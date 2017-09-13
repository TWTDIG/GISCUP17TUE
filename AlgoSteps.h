// Contains the logic used to solve each query,
// called from the main algorithm at Algorithm.h
// Each step changes the set of candidate and result sets,
// where the candidate set indicates a trajectory that MAY
// become a result, but we aren't sure yet, and the result set 
// contains the actual results.
#pragma once

#include "DiHash.h"
#include "Algorithm.h"
#include "Trajectory.h"
#include "FileIO.h"

#include <unordered_set>
#include <map>

// pre-processing steps --------------------------------------------------------------

int slotsPerDimension = 500;
double tolerance = 0.00001;

// Preprocessing step. Inserts start and endpoints in a regular grid so they can be used
// for range queries later.
void addPtsToDiHash(AlgoData &a) {
	a.diHash = new DiHash(*a.boundingBox, slotsPerDimension, tolerance);
	for (Trajectory *t : *a.trajectories) {
		if (t != nullptr) {
			a.diHash->addPoint(t->vertices[0]);
			a.diHash->addPoint(t->vertices[t->size - 1]);
		}
	}
}

// number of simplification steps constructed for each trajectory
int numSimplifications = 4;
int avgs[4];
int count = 0;


// the averages of all learned simplifications
// TODO: is not properly handled for parallel execution
double avgsBBRatio[4];


// Calculates numSimplification trajectory simplifications for one trajectory
void makeSimplificationsForTrajectory(Trajectory &t, double diagonal, AlgorithmObjects &algo, int size) {
	// target ratio of input vertices for simps
	double targets[4] = {.07, .19, .24, .32};

	// convert to integers
	int targetCounts[4];
	for (int i = 0; i < 4; i++) {
		targetCounts[i] = t.size * targets[i];
		targetCounts[i] = std::max(20, targetCounts[i]);
	}
	targetCounts[0] = std::min(18, targetCounts[0]);//start simple in case dihash is useless

	// construct upper and lowerbounds for bsearching epsilon
	double diag = t.boundingBox->getDiagonal();
	double lowerBound = diagonal / 100000;
	double upperBound = diagonal / 2;
	// number of bsearch steps
	int numIterations = 10;

	// constructs all simplifications for t
	for (int i = 0; i < size; i++) {
		double targetVertexPercentage = targets[i];
		TrajectorySimplification* simp = nullptr;
		int tries = 0;
		double newUpperbound = 0;
		// does bsearch on epsilon space
		binaryDoubleSearch(
			[&](double value) -> int {
				newUpperbound = value;
				if (simp != nullptr) delete simp;
				simp = algo.agarwal.simplify(t, value);
				tries++;
				if (tries == 10) {
					return -1;
				}
				else {
					return simp->size > targetCounts[i];
				}
			},
			upperBound,
			lowerBound
		);
		// adjust bsearch parameters to be faster
		upperBound = newUpperbound;
		numIterations -= 2;
		double ratio = simp->size/(double)t.size;
		// apply epsilon learning for query trajectories
		avgsBBRatio[i] += newUpperbound / diagonal;
		t.simplifications.push_back(simp);
	}
	count++;

	/*
	// debug code used to check how close the bsearch is to the target vertex ratio
	for (int i = 0; i < size; i++) {
		TrajectorySimplification* ts = algo.agarwal.simplify(t, simplificationEpsilon / (6 * (i + .25)));
		int diff = ts->size - t.simplifications[i]->size;
		avgs[i] += diff;
	}
	double avg = avgs[2] / (double)count;
	std::cout << "avg " << avg << "\n";
	*/

	// compile portals from simplifications to source trajectory
	for (int i = 0; i < size; i++) {
		for (Portal &p : t.simplifications[i]->portals) {
			// check if it is a useful portal
			if (p.destination - p.source != 1) {
				// check if it is not a duplicate
				bool found = false;
				for (Portal &q : t.simpPortals[p.source]) {
					if (q.destination == p.destination) {
						found = true;
					}
				}
				if (!found) {
					t.simpPortals[p.source].push_back(p);
				}
			}
		}
	}
	//sort portals from small to large (these are many sorts but on small sets)
	for (std::map<int, std::vector<Portal>>::iterator iter = t.simpPortals.begin(); iter != t.simpPortals.end(); ++iter)
	{	
		std::vector<Portal> &k = iter->second;
		std::sort(k.begin(), k.end(), portalCompare);
	}
}

// Calculates numSimplification trajectory simplifications for one trajectory, using guesswork instead of binary search
void makeSourceSimplificationsForTrajectory(Trajectory &t, Trajectory &source, double diagonal, AlgorithmObjects &algo, int size) {
	// apply learned ratio from avgsBBRatio
	for (int i = 0; i < size; i++) {
		double eps = diagonal * (avgsBBRatio[i]/count);
		t.simplifications.push_back(algo.agarwalProg.simplify(t, source, eps));
	}
	// compile portals
	for (int i = 0; i < size; i++) {
		for (Portal &p : t.simplifications[i]->portals) {
			// check if it is a useful portal
			if (p.destination - p.source != 1) {
				// check if it is not a duplicate
				bool found = false;
				for (Portal &q : t.simpPortals[p.source]) {
					if (q.destination == p.destination) {
						found = true;
					}
				}
				if (!found) {
					t.simpPortals[p.source].push_back(p);
				}
			}
		}
	}
	//sort portals from small to large
	for (std::map<int, std::vector<Portal>>::iterator iter = t.simpPortals.begin(); iter != t.simpPortals.end(); ++iter)
	{
		std::vector<Portal> &k = iter->second;
		std::sort(k.begin(), k.end(), portalCompare);
	}
}

void makeSimplificationsForTrajectory(Trajectory &t, AlgorithmObjects &algo) {
	makeSimplificationsForTrajectory(t, t.boundingBox->getDiagonal(), algo, numSimplifications);
}

void loadAndSimplifyTrajectory(std::string &tname, int tIndex, AlgorithmObjects &algo, AlgoData &a) {
	Trajectory *t = algo.fio.parseTrajectoryFile(tname, tIndex);
	if (t->size == 1) {
		delete t;
		// ugly but necessary
		a.trajectories->at(tIndex) = nullptr;
		return;
	}
	// add to bbox used by dihash
	algo.bbox.addPoint(t->boundingBox->minx, t->boundingBox->miny);
	algo.bbox.addPoint(t->boundingBox->maxx, t->boundingBox->maxy);
	makeSimplificationsForTrajectory(*t, algo);
	a.trajectories->at(tIndex) = t;

}


// Mutex guarding access to the dataset from the worker threads
std::mutex simplificationMTX;
// Number of loads/simplifications allocated to a worker as one 'job'
int simplificationBatchSize = 20;

// Returns a trajectory index for the worker to load/simplify, locking the trajectory set
int getConcurrentTrajectory(AlgoData *a) {
	simplificationMTX.lock();
	if (a->startedSimplifying > a->trajectories->size()) {
		simplificationMTX.unlock();
		return -1;
	}
	int returnTrajectory = a->startedSimplifying;
	a->startedSimplifying += simplificationBatchSize;
	simplificationMTX.unlock();
	return returnTrajectory;
}

// obtains a piece of the dataset which it loads and simplifies
void simplificationWorker(AlgoData *a, AlgorithmObjects *algo) {
	int current = getConcurrentTrajectory(a);
	std::vector<std::string> &trajectories = *a->trajectoryNames;
	while (current != -1) {
		int limit = simplificationBatchSize;
		if (current + simplificationBatchSize > trajectories.size()) {
			limit = trajectories.size() - current;
		}
		for (int step = 0; step < limit; step++) {
			std::string &t = trajectories[current + step];
			loadAndSimplifyTrajectory(t, current + step, *algo, *a);
		}
		current = getConcurrentTrajectory(a);
	}
}


std::vector<std::thread*> simplificationThreads;
std::vector<AlgorithmObjects*> algos;


// Preprocessing step. Calculates simplifications for all trajectories in the dataset
void constructSimplifications(AlgoData &a) {
	a.trajectories = new std::vector<Trajectory*>();
	a.trajectories->resize(a.numTrajectories);
	// spawn workers which load/simplify
	for (int i = 0; i < a.numWorkers; i++) {
		AlgorithmObjects *algo = new AlgorithmObjects();
		std::thread *t = new std::thread(simplificationWorker, &a, algo);
		simplificationThreads.push_back(t);
		algos.push_back(algo);
	}
	// kill all workers
	for (int i = 0; i < a.numWorkers; i++) {
		(*simplificationThreads[i]).join();
		AlgorithmObjects *algo = algos[i];
		a.boundingBox->addPoint(algo->bbox.minx, algo->bbox.miny);
		a.boundingBox->addPoint(algo->bbox.maxx, algo->bbox.maxy);
		delete simplificationThreads[i];
		delete algo;
	}

	//cleanup
	algos.clear();
	simplificationThreads.clear();
	std::cout << "done";
}






// runtime steps ---------------------------------------------------------------------

// to reduce bookkeeping time, and to increase memory locality, we use 
// all pruning step on a single trajectory pair before considering another,
// we implement this by giving each pruning step emittor functions indicating
// success, failure and unknown results

// Query step. Does rangequeries for start/endpoints of dataset. Adds all found trajectories
// to candidates.
void collectDiHashPoints(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, const std::function< void(Trajectory*) >& emit) {
	Vertex start = queryTrajectory.vertices[0];
	Vertex end = queryTrajectory.vertices[queryTrajectory.size - 1];

	a->diHash->neighborsWithCallback(start, end, q.queryDelta, *a->trajectories, [&](Trajectory* t) -> void{
		emit(t);
	});
}


// Query step. For each trajectory T in the dataset and query trajectory Q, this step
// compares successive simplifications of T and Q with continuous decision frechet.
// Each comparison can result in YES, NO, or MAYBE.
// YES   -> remove from candidates, add to results
// NO    -> remove from candidates
// MAYBE -> try next simplification, if none are left, continue to next algorithm step
void pruneWithSimplifications(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t, 
	const std::function< void(Trajectory*) >& maybe, const std::function< void(Trajectory*) >& result) {
	bool broke = false;
	for (int i = 0; i < numSimplifications; i++) {

		// construct epsilons for tri ineq.
		double decisionEpsilonLower = q.queryDelta
			- queryTrajectory.simplifications[i]->simplificationEpsilon
			- t->simplifications[i]->simplificationEpsilon;

		double decisionEpsilonUpper = q.queryDelta
			+ queryTrajectory.simplifications[i]->simplificationEpsilon
			+ t->simplifications[i]->simplificationEpsilon;

		double dist = equalTimeDistance(*t->simplifications[i], *queryTrajectory.simplifications[i]);

		// do ETD greedy check
		if (dist < decisionEpsilonLower) {
			result(t);
			broke = true;
			break;
		}

		// do lower frechet check
		if (decisionEpsilonLower > 0) {
			bool r = algo->cdfqs.calculate(*queryTrajectory.simplifications[i], *t->simplifications[i], decisionEpsilonLower, q.queryDelta);
			if (r) {
				result(t);
				broke = true;
				break;
			}
		}

		// do upper frechet check
		if (decisionEpsilonUpper > 0) {
			bool r = algo->cdfqs.calculate(*queryTrajectory.simplifications[i], *t->simplifications[i], decisionEpsilonUpper, q.queryDelta);
			if (!r) {
				broke = true;
				break;
			}
		}
	}
	if (!broke) {
		maybe(t);
	}
}

// Query step. Uses equal time distance as an upperbound for the actual frechet distance
// If ETD(P, Q) <= queryDelta then CDF(P,Q) <= queryDelta. With P in dataset and Q query trajectory.
void pruneWithEqualTime(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t,
	const std::function< void(Trajectory*) >& maybe, const std::function< void(Trajectory*) >& result) {
	double dist = equalTimeDistance(*t, queryTrajectory);
	if (dist < q.queryDelta) {
		result(t);
	}
	else {
		maybe(t);
	}
}


// Query step. The final step for each query is to do a full decision frechet computation.
// This step contains no additional smart optimization, and so is very slow.
void pruneWithDecisionFrechet(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t,
	const std::function< void(Trajectory*) >& result) {
	bool r = algo->cdfqs.calculate(queryTrajectory, *t, q.queryDelta);
	if (r) {
		result(t);
	}
	// last step, conclusive, no maybe
}