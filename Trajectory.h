
#pragma once

#include "Vertex.h"
#include "BoundingBox.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>




struct Portal {
	int source;
	int destination;
	double distance;
};

bool portalCompare(Portal &lhs, Portal &rhs) { return lhs.destination < rhs.destination; }

class TrajectorySimplification;

class Trajectory {
public: 
	std::string name;
	
	std::vector<Vertex> vertices;// all vertices of traj
	std::vector<double> distances;// distance between vertices
	std::vector<double> totals;// total length at vertex x
	std::vector<int> sourceIndex;// mapping back to the source trajectory if this is a simplification
	std::map<int, std::vector<Portal>> simpPortals;// freespace jumps

	int size;
	int uniqueIDInDataset;
	double totalLength;

	BoundingBox *boundingBox;
	std::vector<TrajectorySimplification*> simplifications;

	~Trajectory() {
		delete boundingBox;
	}

	void print() {

		std::cout << "Trajectory: " << name << "\n";
		for (int i = 0; i < size; i++) {
			Vertex v = vertices[i];
			std::cout << i << " " << v.x << " " << v.y << " " << distances[i] << " " << totals[i] << "\n";
		}
	}
};

class TrajectorySimplification : public Trajectory {

public:
	Trajectory* source;
	std::vector<Portal> portals;
	double simplificationEpsilon;// epsilon used in agarwal to make traj
};