#pragma once


#include "Vertex.h"
#include "Trajectory.h"
#include "DoubleNSearch.h"
#include "EqualTimeDistance.h"

#include <vector>

// This class contains all logic need to compute a simplification
// of any input trajectory using agarwal simplification with
// double & search. The algorithm uses EqualTimeDistance.h to
// satisfy the agarwal constraints.
class AgarwalSimplification {
	int doubleNSearchBase = 2;
	double doubleNSearchExponentStep = 1;

	// temp data so we dont have to reallocate
	std::vector<Vertex> simpBuffer;
	std::vector<double> simpDistances;
	std::vector<double> simpTotals;

public:
	// wrapper for the simplify function
	TrajectorySimplification* simplify(Trajectory &t, double simplificationEpsilon) {
		TrajectorySimplification* simplified = new TrajectorySimplification();
		simplified->name = t.name + "[simplified]";
		simplified->simplificationEpsilon = simplificationEpsilon;
		simplified->source = &t;

		simplify(t, *simplified, simplificationEpsilon);

		return simplified;
	}

	// uses agarwal to simplify (t) into (simplification) with agarwal epsilon (simplificationEpsilon)
	void simplify(Trajectory &t, TrajectorySimplification &simplification, double simplificationEpsilon) {
		// reset temp data
		simpBuffer.clear();
		simpDistances.clear();
		simpTotals.clear();

		// initialize first vertex
		std::vector<Vertex> &P = t.vertices;
		simpBuffer.push_back(P[0]);
		simpDistances.push_back(0);
		simpTotals.push_back(0);

		int simpSize = 1;

		int rangeStart = 1;
		int prevk = 0;
		while (true) {
			// find last index of (t) satisfying (simplificationEpsilon), given current index (prevk)
			int k = findLastFrechetMatch(P, simpBuffer, t.totals, simpTotals, t.distances, simpDistances, simpSize, rangeStart, t.size, prevk, simplificationEpsilon, simplification.portals);
			// put vertex (k) of (t) into simplification
			simpSize++;
			simpBuffer[simpSize - 1] = P[k];
			// check if we reached the end
			if (k == t.size - 1) {
				break;
			}
			prevk = k;
			rangeStart = k + 1;
		}

		// copy temp data into simplification
		simplification.size = simpBuffer.size();
		simplification.vertices = simpBuffer;
		simplification.distances = simpDistances;
		simplification.totals = simpTotals;
	}

private:
	// Finds index k of last vertex v that still satisfies 
	// the simplification epsilon 
	int findLastFrechetMatch(
		std::vector<Vertex> &P, std::vector<Vertex> &simp,
		std::vector<double> &ptotals, std::vector<double> &simptotals,
		std::vector<double> &pdists, std::vector<double> &simpdists,
		int simpSize,
		int start, int end, int prevk, double epsilon,
		std::vector<Portal> &portals) {

		// add temporary data which we modify with double & search
		simp.push_back(P[0]);
		simpdists.push_back(0);
		simpTotals.push_back(0);

		// Use lambda's to easily double & search the function from (start) to (end)
		return doubleNsearch(
			[&](int index) -> bool {
				// update temporary data to reflect double&search choice
				simp[simpSize] = P[index];
				double dx = simp[simpSize].x - simp[simpSize - 1].x;
				double dy = simp[simpSize].y - simp[simpSize - 1].y;
				double d = sqrt(dx*dx + dy*dy);
				simpdists[simpSize] = d;
				simptotals[simpSize] = simptotals[simpSize - 1] + d;
				// calculate upper bound to subtrajectory frechet with ETD
				double dist = equalTimeDistance(
					P, simp,
					ptotals, simptotals,
					pdists, simpdists,
					index + 1, simpSize + 1,
					prevk, simpSize - 1
				);
				// tell double & search whether its a match
				return dist <= epsilon;
			},
			start,
			end,
			doubleNSearchBase,
			doubleNSearchExponentStep
		);
	}


};