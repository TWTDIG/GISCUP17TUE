// Contains various geometric utility functions used by frechet algorithms.
#pragma once

#include "Vertex.h"

#include <algorithm>

struct Range {
	double start;
	double end;
};

Range emptyRange = { 0,0 };
Range dontCare = { 0, 0 };

bool isEmpty(Range &r) {
	return r.start == r.end;
}

bool isComplete(Range &r) {
	return r.start == 0.0 && r.end == 1.0;
}

void setRange(Range &r, Range &s) {
	r.start = s.start;
	r.end = s.end;
}

void setRange(Range &r, double start, double end) {
	r.start = start;
	r.end = end;
}

double distSQ(Vertex &p, Vertex &q) {
	double dx = p.x - q.x;
	double dy = p.y - q.y;
	return dx*dx + dy*dy;
}

double clamp01(double input) {
	return input > 1 ? 1 : (input < 0 ? 0 : input);
}

// by WM
inline bool computeInterval(Vertex &a, Vertex &b1, Vertex &b2, double eps, Range &r) {
	// compute the interval along b1-b2 that is within eps-range of a

	// L(t) = b1 + t * (b2 - b1)
	// D(t) = |L(t) - a|
	// for which t, does D(t) = eps hold?
	// square it to make it easier
	// (b1_x + t * (b2_x - b1_x) - a_x)^2
	//    + (b1_y + t * (b2_y - b1_y) - a_y)^2 = eps^2
	// so, we get:
	// A * t^2 + B * t + C = 0 with
	// A = (b2_x - b1_x)^2 + (b2_y - b1_y)^2
	// B = 2 * ((b2_x - b1_x)*(b1_x - a_x) + (b2_y - b1_y)*(b1_y - a_y));
	// C = (b1_x - a_x)^2 + (b1_y - a_y)^2 - eps^2

	// pull out some of the identical computations in here
	double b2m1x = b2.x - b1.x;
	double b2m1y = b2.y - b1.y;
	double b1max = b1.x - a.x;
	double b1may = b1.y - a.y;

	double A = b2m1x * b2m1x + b2m1y * b2m1y;
	double B = 2 * ((b2m1x) * (b1max) + (b2m1y) * (b1may));
	double C = b1max*b1max + b1may * b1may - eps*eps;

	double D = B * B - 4 * A * C;
	if (D < 0) {
		// no solution
		return false;
	}
	else {
		// pull out the Sqrt(D)
		double sqrtD = sqrt(D);
		double t1 = (-B + sqrtD) / (2 * A);
		double t2 = (-B - sqrtD) / (2 * A);
		// rather than doing a swap, if may be faster to check the sign of A before doing the assignment, OR just use min/max, etc 
		double tempt1 = t1;
		t1 = std::min(t1, t2);
		t2 = std::max(tempt1, t2);

		if (t2 < 0 || t1 > 1) {
			return false;
		}
		else {
			r.start = std::max(0.0, t1);
			r.end = std::min(1.0, t2);
			return true;
		}
	}
};



