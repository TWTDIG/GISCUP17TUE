#pragma once

#include "Vertex.h"
#include <cfloat>
#include <algorithm>
#include <cmath>

// Represents a boundingbox around a trajectory or group of trajectories
class BoundingBox {
public:
	double minx = DBL_MAX;
	double miny = DBL_MAX;

	double maxx = -DBL_MAX;
	double maxy = -DBL_MAX;

	void addPoint(double x, double y) {
		if (x < minx) minx = x;
		if (y < miny) miny = y;

		if (x > maxx) maxx = x;
		if (y > maxy) maxy = y;
	}

	double getDiagonal() {
		double width = maxx - minx;
		double height = maxy - miny;

		return sqrt(width*width + height*height);
	}
};