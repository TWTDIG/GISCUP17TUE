#pragma once

#include "Trajectory.h"

#include <string>

// Represents a query in a queryset file
// The queryNumber is the index in the query file,
// and should be used when outputting results
struct Query {
	std::string queryTrajectoryFilename;
	double queryDelta;
	int queryNumber;
};