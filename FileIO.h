// Contains wrappers for all file IO needed for the algorithm
// Such as loading trajectories, query files, dataset files, etc
#pragma once

#include "Trajectory.h"
#include "Query.h"
#include "Settings.h"

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <iomanip>

class FileIO {
	std::vector<Vertex> vertexBuffer;
	std::vector<double> distanceBuffer;
	std::vector<double> totalBuffer;
	std::vector<int> sourceIndex;


public:

	// Parses trajectory file, also computes trajectory metrics
	// TODO: move temporary buffer allocation out of function
	Trajectory* parseTrajectoryFileStreams(std::string filename, int trajectoryNumber) {
		std::ifstream infile(filename);



		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}
		std::string line;
		std::getline(infile, line);
		double x, y, z, w;
		vertexBuffer.clear();
		distanceBuffer.clear();
		totalBuffer.clear();
		sourceIndex.clear();

		distanceBuffer.push_back(0);
		totalBuffer.push_back(0);

		BoundingBox *b = new BoundingBox();

		bool start = true;

		Trajectory *t = new Trajectory();
		t->name = filename;
		Vertex v;

		while (infile >> x >> y >> z >> w)
		{
			v.x = x;
			v.y = y;
			v.trajectoryNumber = trajectoryNumber;
			v.isStart = start;
			start = false;
			//update boundingbox
			b->addPoint(v.x, v.y);
			if (vertexBuffer.empty()) {
				vertexBuffer.push_back(v);
				sourceIndex.push_back(sourceIndex.size());
			}
			else {
				// ignore duplicate verts, they are annoying
				int prevIndex = vertexBuffer.size() - 1;
				Vertex &prev = vertexBuffer[prevIndex];
				if (prev.x != v.x || prev.y != v.y) {
					double dx = v.x - prev.x;
					double dy = v.y - prev.y;
					double dist = sqrt(dx*dx + dy*dy);
					distanceBuffer.push_back(dist);
					totalBuffer.push_back(totalBuffer[prevIndex] + dist);
					vertexBuffer.push_back(v);
					sourceIndex.push_back(sourceIndex.size());
				}
			}
		}
		infile.close();

		t->size = vertexBuffer.size();
		t->totalLength = totalBuffer[t->size - 1];
		t->boundingBox = b;

		t->vertices = vertexBuffer;
		t->distances = distanceBuffer;
		t->totals = totalBuffer;
		t->sourceIndex = sourceIndex;



		return t;
	}

	#define BUFFER_SIZE (1024 * 1024 * 1024)//1Mb
	#define READ_SIZE (1024 * 1024 * 5)//5kb

	//std::vector<char> buffer = std::vector<char>(BUFFER_SIZE);
	char* buffer = nullptr;
	Trajectory* parseTrajectoryFileFast(std::string filename, int trajectoryNumber) {
		if (buffer == nullptr) {
			buffer = new char[BUFFER_SIZE];
		}
#if USE_FOPEN_S
		FILE* file;
		int err = fopen_s(&file, filename.c_str(), "rb");
		if (err != 0) {
			std::cout << "Failed to open (fopens): " << filename << "\n";
			exit(1);
		}
#else
		FILE* file = fopen(filename.c_str(), "rb");
		if (file == NULL) {
			std::cout << "Failed to open (fopen): " << filename << "\n";
			exit(1);
		}
#endif


		double x, y, z, w;
		vertexBuffer.clear();
		distanceBuffer.clear();
		totalBuffer.clear();
		sourceIndex.clear();

		distanceBuffer.push_back(0);
		totalBuffer.push_back(0);

		BoundingBox *b = new BoundingBox();

		bool start = true;

		Trajectory *t = new Trajectory();
		t->name = filename;
		Vertex v;

		int line = 0;

		int offset = 0;
		int reads = 0;

		while (true)
		{
			size_t readLength = fread(&buffer[READ_SIZE*reads], 1, READ_SIZE, file);
			reads++;
			while (true) {
				if (offset == readLength) break;
				line++;
				int eol = -1;
				for (int i = offset; i < offset + readLength; i++) {
					if (buffer[i] == '\n') {
						eol = i;
						break;
					}
				}
				int prevOffset = offset;
				int nextOffset = eol + 1;
				if (eol == -1) break;
				if (line == 1) {
					offset = eol + 1;
					continue;
				}
				// vertex from 0 to nl
				int space = -1;
				for (int i = offset; i < eol; i++) {
					if (buffer[i] == ' ') {
						space = i;
						break;
					}
				}
				offset = nextOffset;
				if (space == -1 || space == prevOffset) break;// bad line
				char* pEnd;
				double x = strtod(&buffer[prevOffset], &pEnd);
				double y = strtod(pEnd, NULL);
				v.x = x;
				v.y = y;
				v.trajectoryNumber = trajectoryNumber;
				v.isStart = start; start = false;
				//update boundingbox
				b->addPoint(v.x, v.y);
				if (vertexBuffer.empty()) {
					vertexBuffer.push_back(v);
					sourceIndex.push_back(sourceIndex.size());
				}
				else {
					// ignore duplicate verts, they are annoying
					int prevIndex = vertexBuffer.size() - 1;
					Vertex &prev = vertexBuffer[prevIndex];
					if (prev.x != v.x || prev.y != v.y) {
						double dx = v.x - prev.x;
						double dy = v.y - prev.y;
						double dist = sqrt(dx*dx + dy*dy);
						distanceBuffer.push_back(dist);
						totalBuffer.push_back(totalBuffer[prevIndex] + dist);
						vertexBuffer.push_back(v);
						sourceIndex.push_back(sourceIndex.size());
					}
				}
			}

			//EOF
			if (readLength != BUFFER_SIZE) {
				break;
			}
		}
		fclose(file);


		t->size = vertexBuffer.size();
		t->totalLength = totalBuffer[t->size - 1];
		t->boundingBox = b;

		t->vertices = vertexBuffer;
		t->distances = distanceBuffer;
		t->totals = totalBuffer;
		t->sourceIndex = sourceIndex;



		return t;
	}

	// delegating function for file loading
	Trajectory* parseTrajectoryFile(std::string filename, int trajectoryNumber) {
		#if USE_FAST_IO
			return parseTrajectoryFileFast(TRAJECTORY_FILES_OFFSET + filename, trajectoryNumber);
		#else
			return parseTrajectoryFileStreams(TRAJECTORY_FILES_OFFSET + filename, trajectoryNumber);
		#endif
	}


	// Parses query file, does not load query trajectories
	std::vector<Query>* parseQueryFile(char* filename) {
		std::ifstream infile(filename);

		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		double eps;
		std::string queryTrajectoryFileName;

		std::vector<Query> *queries = new std::vector<Query>;

		int queryNumber = 0;
		while (infile >> queryTrajectoryFileName >> eps) {
			Query q;
			q.queryNumber = queryNumber;
			q.queryDelta = eps;
			q.queryTrajectoryFilename = queryTrajectoryFileName;

			queries->push_back(q);
			queryNumber++;
		}
		infile.close();

		return queries;

	}

	// Parses dataset file, does not translate filenames based on dataset location
	// This behavior is consistent with the conversations on the mailing list
	std::vector<std::string>* parseDatasetFile(std::string filename) {
		std::ifstream infile(filename);

		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		std::string trajectoryFileName;
		std::vector<std::string> * trajectoryNames = new std::vector<std::string>();

		int trajectoryNumber = 0;
		while (infile >> trajectoryFileName) {
			trajectoryNames->push_back(trajectoryFileName);
		}
		infile.close();

		return trajectoryNames;
	}

	// Writes result trajectories for one query to a file called result-XXXXX.txt 
	void writeQueryOutputFile(Query &q, std::vector<std::string> results) {
		std::ostringstream stringStream;
		stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
		std::string filename = stringStream.str();
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		for (std::string trajectoryName : results) {
			outfile << trajectoryName << "\n";
		}
		outfile.close();
	}

	void writeQueryOutputFile(Query &q, std::string results) {
		std::ostringstream stringStream;
		stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
		std::string filename = stringStream.str();
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}
		outfile << results;
		outfile.close();
	}
};