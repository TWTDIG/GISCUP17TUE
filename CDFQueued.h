#pragma once


#include "Vertex.h"
#include "FrechetUtil.h"

#include <algorithm>
#include <vector>
#include <stdio.h>



// Only computes reachable part of freespace diagram
// unused in algorithm, but same as CDFQueued without jumping administration
// used for performance comparisons and correctness testing
class CDFQueued {

private:
	struct QEntry {
		int row_index;
		double lowest_right;
	};

	std::vector<QEntry> queue[2];
	int queueSize[2];

	double dist(Vertex p, Vertex q) {
		double dx = p.x - q.x;
		double dy = p.y - q.y;
		return sqrt(dx*dx + dy*dy);
	}


public:

	int numRows = 0;
	bool calculate(
		std::vector<Vertex> &P, std::vector<Vertex> &Q,
		int offset_p, int offset_q,
		int size_p, int size_q,
		double queryDelta
	) {

		double startDist = dist(P[offset_p], Q[offset_q]);
		double endDist = dist(P[size_p - 1], Q[size_q - 1]);
		if (startDist > queryDelta || endDist > queryDelta) return false;
		if (size_p <= offset_p + 1 || size_q <= offset_q + 1) return false;//TODO: do we need this? // WM: added offset

		int first = 0;
		int second = 1;

		Range Rf;
		Range Tf;

		// ensure queue capacity
		// WM: added offsets
		int max = size_p - offset_p;
		if (size_q - offset_q > max) max = size_q - offset_q;
		if (queue[0].size() < max) {
			for (int i = queue[0].size(); i < max; i++) {
				queue[0].push_back({ 0,0 });
				queue[1].push_back({ 0,0 });
			}
		}

		// setup
		//WM: this is always free space! by check on startDist
		//bool LFree = computeInterval(Q[0], P[0], P[1], queryDelta, Rf); 
		queue[first][0].row_index = 0;
		queue[first][0].lowest_right = 0;

		queueSize[first] = 1;
		queueSize[second] = 0;


		// For each column
		for (int column = offset_q; column < size_q - 1; column++) {
			if (queueSize[first] == 0) {
				// nothing reachable anymore
				return false;
			}
			queueSize[second] = 0;
			int row = queue[first][0].row_index;
			int qIndex = 0;
			// while there's reachable cells left in the queue
			while (qIndex < queueSize[first]) {
				double left_most_top = 2;
				// start at reachable cell at the head of the queue, and continue until
				// reachability cannot propagate, consuming the queue as we progress
				do {
					// tracks whether we overshoot the queue
					bool outsideQueue = qIndex >= queueSize[first];
					// Right edge stored in Rf, RFree = false means not free
					bool RFree = computeInterval(Q[column + 1], P[row], P[row + 1], queryDelta, Rf);
					if (RFree) {
						if (left_most_top <= 1) {
							// push to queue
							queue[second][queueSize[second]].row_index = row;
							queue[second][queueSize[second]].lowest_right = Rf.start;
							queueSize[second]++;
						}
						else {
							// WM: think you should be checking row here as well
							if (!outsideQueue && row == queue[first][qIndex].row_index && queue[first][qIndex].lowest_right <= Rf.end) {
								// push to queue
								queue[second][queueSize[second]].row_index = row;
								queue[second][queueSize[second]].lowest_right = std::max(queue[first][qIndex].lowest_right, Rf.start);
								queueSize[second]++;
							}
						}
					}
					// Top edge stored in Tf, TFree = false means not free
					bool TFree = computeInterval(P[row + 1], Q[column], Q[column + 1], queryDelta, Tf);
					if (!outsideQueue && row == queue[first][qIndex].row_index) {
						// consume the first queue
						qIndex++;
						if (TFree) {
							left_most_top = Tf.start;
						}
						else {
							left_most_top = 2;
						}
					}
					else if (TFree && left_most_top <= Tf.end) {
						left_most_top = std::max(left_most_top, Tf.start);
					}
					else {
						left_most_top = 2;
					}
					// propagated reachability by one cell, so look at next row
					row++;
					numRows++;
				} while (left_most_top <= 1 && row < size_p - 1);
			}

			// swap first and second column
			int temp = first;
			first = second;
			second = temp;
		}

		int endIndex = queueSize[first] - 1;
		if(endIndex<0) return false;
		return queue[first][endIndex].row_index == size_p - 2 && queue[first][endIndex].lowest_right <= 1;
	}

	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta);
	}
};
