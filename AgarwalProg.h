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
class ProgressiveAgarwal {
	int doubleNSearchBase = 2;
	double doubleNSearchExponentStep = 1;

	// temp data so we dont have to reallocate
	std::vector<Vertex> simpBuffer;
	std::vector<double> simpDistances;
	std::vector<double> simpTotals;
	std::vector<int> sourceIndex;

public:
	// wrapper for the simplify function
	TrajectorySimplification* simplify(Trajectory &parent, Trajectory &sourceTrajectory, double simplificationEpsilon) {
		TrajectorySimplification* simplified = new TrajectorySimplification();
		simplified->name = parent.name + "[simplified]";
		simplified->simplificationEpsilon = simplificationEpsilon;
		simplified->source = &parent;

		simplify(parent, *simplified, sourceTrajectory, simplificationEpsilon);
		
		return simplified;
	}

	// uses agarwal to simplify (t) into (simplification) with agarwal epsilon (simplificationEpsilon)
	void simplify(Trajectory &parent, TrajectorySimplification &simplification, Trajectory &sourceTrajectory, double simplificationEpsilon) {
		// reset temp data
		simpBuffer.clear();
		simpDistances.clear();
		simpTotals.clear();
		sourceIndex.clear();

		// initialize first vertex
		std::vector<Vertex> &P = parent.vertices;
		simpBuffer.push_back(P[0]);
		simpDistances.push_back(0);
		simpTotals.push_back(0);
		sourceIndex.push_back(0);

		int simpSize = 1;

		int rangeStart = 1;
		int prevk = 0;
		while (true) {
			// find last index of (t) satisfying (simplificationEpsilon), given current index (prevk)
			int k = findLastFrechetMatch(P, simpBuffer, parent.totals, simpTotals, parent.distances, simpDistances, simpSize, rangeStart, parent.size, prevk, simplificationEpsilon, parent.sourceIndex, sourceTrajectory, simplification.portals);
			// put vertex (k) of (t) into simplification
			simpSize++;
			simpBuffer[simpSize - 1] = P[k];
			sourceIndex.push_back(parent.sourceIndex[k]);
			// check if we reached the end
			if (k == parent.size - 1) {
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
		simplification.sourceIndex = sourceIndex;
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
		std::vector<int> &parentSourceIndices,
		Trajectory &sourceTrajectory,
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
				int end = 0;
				int start = parentSourceIndices[prevk];
				if (index + 1 >= parentSourceIndices.size()) {
					end = sourceTrajectory.size;
				}
				else {
					end = parentSourceIndices[index + 1];
				}
				// calculate upper bound to subtrajectory frechet with ETD
				double dist = equalTimeDistance(
					sourceTrajectory.vertices, simp,
					sourceTrajectory.totals, simptotals,
					sourceTrajectory.distances, simpdists,
					end, simpSize + 1,
					start, simpSize - 1
				);
				// construct freespace jump (portal) from data
				Portal p;
				p.source = prevk;
				p.destination = index;
				p.distance = dist;
				portals.push_back(p);
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