// Contains the main structure of the algorithm
// The actual logic used to solve each query is in AlgoSteps.h
#pragma once

#include "Agarwal.h"
#include "AgarwalProg.h"
#include "EqualTimeDistance.h"
#include "DiHash.h"
#include "Query.h"
#include "CDFQueued.h"
#include "CDFQShortcuts.h"
#include "settings.h"


#include <thread>
#include <chrono>
#include <mutex>
#include <sstream>
#include <iomanip>
#include <iostream>


// All data needed by the algorithm to solve a specific query file
// Also contains structures needed for preprocessing
struct AlgoData {
	std::vector<Query> *queries;
	std::vector<Trajectory*> *trajectories;
	std::vector<std::string> *trajectoryNames;
	int numTrajectories;
	DiHash* diHash;
	FileIO fio;
	BoundingBox* boundingBox;
	volatile int startedSolving = 0;
	volatile int startedSimplifying = 0;
	int numWorkers;


};

// Data needed by each worker thread executing the algorithm
// Each thread has seperate algorithm objects so they can use
// their own buffers to store intermediate results.
struct AlgorithmObjects {
	std::ostringstream results;
	std::vector<Trajectory*> candidates;
	BoundingBox bbox;

	FileIO fio;
	AgarwalSimplification agarwal;
	ProgressiveAgarwal agarwalProg;
	CDFQueued cdfq;
	CDFQShortcuts cdfqs;

};

// TODO: Included here to avoid include problem
#include "AlgoSteps.h"

// Does all needed preprocessing for the given the dataset
void preprocessDataSet(AlgoData *a) {
	constructSimplifications(*a);
	addPtsToDiHash(*a);
}

// Solves a single query, calls functions in AlgoSteps.h
// Before solving, also loads the query trajectory (since it may 
// not be present in the dataset) and constructs simplifications for it.
void solveQuery(AlgoData *a, Query &q, AlgorithmObjects *algo) {
	Trajectory *queryTrajectory = algo->fio.parseTrajectoryFile(q.queryTrajectoryFilename, -1);


	double diagonal = queryTrajectory->boundingBox->getDiagonal();
	makeSourceSimplificationsForTrajectory(*queryTrajectory, *queryTrajectory, diagonal, *algo, numSimplifications);
	// for query trajectories, we also simplify the simplifications. Not because we use them directly, but because
	// we use their freespace jumps
	for (int i = 1; i < numSimplifications; i++) {
		makeSourceSimplificationsForTrajectory(*queryTrajectory->simplifications[i], *queryTrajectory, diagonal, *algo, i-1);
	}



#if WRITE_OUTPUT_TO_QUERY 
	std::ostringstream stringStream;
	stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
	std::string filename = stringStream.str();
	std::ofstream outfile(filename);

	if (!outfile.is_open()) {
		std::cout << "Failed to open: " << filename << "\n";
		exit(1);
	}
#endif

	// statistics
	int results = 0;
	int dihash = 0;
	int simp = 0;
	int et = 0;

	const std::function< void(Trajectory*) >& result = [&](Trajectory *t) -> void{
		results++;
#if WRITE_OUTPUT_TO_QUERY 
		outfile << t->name << "\n";
#endif
	};

	// the actual pruning steps of the algorithm
	collectDiHashPoints(a, q, algo, *queryTrajectory, [&](Trajectory *t) -> void {
	//for (Trajectory *t : *a->trajectories) {
		dihash++;
		pruneWithSimplifications(a, q, algo, *queryTrajectory, t, [&](Trajectory *t) -> void {
			simp++;
			pruneWithEqualTime(a, q, algo, *queryTrajectory, t, [&](Trajectory *t) -> void {
				et++;
				pruneWithDecisionFrechet(a, q, algo, *queryTrajectory, t, result);
			}, result);
		}, result);
	//}
	});

#if WRITE_OUTPUT_TO_QUERY 
	outfile.close();
#endif

	// TODO: cleanup queryTrajectory, doesn't work from destructor somehow
	for (int i = 0; i < queryTrajectory->simplifications.size(); i++) {
		Trajectory *s = queryTrajectory->simplifications[i];
		for (int j = 0; j < s->simplifications.size(); j++) {
			delete s->simplifications[j];
		}
		delete queryTrajectory->simplifications[i];
	}
	delete queryTrajectory;

	return;

}

// Mutex guarding access to the queryset from the worker threads
std::mutex queryMtx;
// Number of queries allocated to a worker as one 'job'
int querySteps = 20;

// Returns a query index for a worker to solve, locking the query set
int getConcurrentQuery(AlgoData *a) {
	queryMtx.lock();
	if (a->startedSolving > a->queries->size()) {
		queryMtx.unlock();
		return -1;
	}
	int returnQuery = a->startedSolving;
	a->startedSolving += querySteps;
	queryMtx.unlock();
	if (returnQuery % 100 == 0) {
		std::cout << " --- Solving: " << returnQuery << "\n";
	}
	return returnQuery;
}

// Function executed by the worker threads, obtains a query index
// solves a fixed number of queries from that index, then
// tries to obtain a new index. When no queries are left, it exits
// and merges its statistics with the complete statistics.
void worker(AlgoData *a, AlgorithmObjects *algo) {
	int current = getConcurrentQuery(a);
	std::vector<Query> &queries = *a->queries;
	while (current != -1) {
		int limit = querySteps;
		if (current + querySteps > queries.size()) {
			limit = queries.size() - current;
		}
		for (int step = 0; step < limit; step++) {
			Query &c = queries[current + step];
			solveQuery(a, c, algo);
		}
		current = getConcurrentQuery(a);
	}
	delete algo;
}


void printMS(std::string msg, long ms) {
	double timeSec = ms / 1000.0;
	std::cout << msg << ": " << timeSec << " sec \n";
}

void print(std::string msg, double value) {
	std::cout << msg << ": " << value << "\n";
}


// Datastructure containing worker threads
std::vector<std::thread*> threads;

// Spins up all worker threads, waits for them to complete,
// then prints statistics.
void solveQueries(AlgoData *a) {
	for (int i = 0; i < a->numWorkers; i++) {
		std::thread *t = new std::thread(worker, a, new AlgorithmObjects());
		threads.push_back(t);
	}
	for (int i = 0; i < a->numWorkers; i++) {
		(*threads[i]).join();
		delete threads[i];
	}
}

void cleanup(AlgoData *a) {
	//TODO: improve code by moving deallocations here
	//TODO: not strictly necessary because program exits
}

// Entrypoint for the algorithm
void runAlgorithm(AlgoData *a) {
	long timeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	std::cout << " - Preprocess\n";
	preprocessDataSet(a);
	long ptimeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	std::cout << " - Solve\n";
	solveQueries(a);
	std::cout << " - Cleanup\n";
	cleanup(a);
	long total = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1) - timeMS;
	long totalp = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1) - ptimeMS;
	printMS("TOTAL", total);
	printMS("TOTAL_SOLVE", totalp);
	printMS("PREPROCESSING", total - totalp);

}


