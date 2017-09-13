

#pragma once

#include "Trajectory.h"
#include <string>

struct Vertex {
	double x;
	double y;

	// unique ID of trajectory in dataset
	int trajectoryNumber;

	// true -> this vertex is the start of trajectory
	// probably not the most elegant
	bool isStart;
};