#pragma once


#include "Vertex.h"
#include "FrechetUtil.h"

#include <algorithm>
#include <vector>
#include <stdio.h>



// Only computes reachable part of freespace diagram, uses shortcuts to skip columns in reachable part
class CDFQShortcuts {

private:
	// entry in the queue
	struct QEntry {
		int start_row_index;
		int end_row_index;
		double lowest_right;
	};

	// queues for column-based frechet algo
	std::vector<QEntry> queue[2];
	int queueSize[2];

	// wrapper distance function
	double dist(Vertex p, Vertex q) {
		double dx = p.x - q.x;
		double dy = p.y - q.y;
		return sqrt(dx*dx + dy*dy);
	}

	// compute frechet distance of point vs line
	inline double computeSegmentFrechet(
		Portal &p,
		int q,
		std::vector<Vertex> &p_array, std::vector<Vertex> &q_array
	) {
		Vertex &pstart = p_array[p.source];
		Vertex &pend = p_array[p.destination];

		Vertex &qstart = q_array[q];
		Vertex &qend = q_array[q];

		return computeSegmentFrechet(pstart, pend, qstart, qend);
	}

	inline double computeSegmentFrechet(
		Vertex &pstart, Vertex &pend,
		Vertex &qstart, Vertex &qend) {
		double startdx = pstart.x - qstart.x;
		double startdy = pstart.y - qstart.y;
		double enddx = pend.x - qend.x;
		double enddy = pend.y - qend.y;
		double startdist = startdx*startdx + startdy*startdy;
		double enddist = enddx*enddx + enddy*enddy;
		return sqrt(std::max(startdist, enddist));
	}



public:

	// debug variable
	int numRows = 0;

	// calculate frechet decision between P, and Q given queryDelta
	bool calculate(
		std::vector<Vertex> &P, std::vector<Vertex> &Q,
		int offset_p, int offset_q,
		int size_p, int size_q,
		double queryDelta,
		double baseQueryDelta,
		std::map<int, std::vector<Portal>> &portals
	) {
		double startDist = dist(P[offset_p], Q[offset_q]);
		double endDist = dist(P[size_p - 1], Q[size_q - 1]);
		if (startDist > queryDelta || endDist > queryDelta) return false;
		if (size_p <= offset_p + 1 || size_q <= offset_q + 1) return false;//TODO: do we need this? // WM: added offset

		int first = 0;
		int second = 1;

		Range Rf;
		Range Tf;
		Portal choice;

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
		queue[first][0].start_row_index = 0;
		queue[first][0].end_row_index = 0;
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
			int row = queue[first][0].start_row_index;
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
							double newLR = Rf.start;
							if (isComplete(Rf) && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row - 1) {
								// complete reachable right means increase previous queue entry to span this cell
								queue[second][queueSize[second] - 1].end_row_index = row;
							}
							else {
								// push to queue
								queue[second][queueSize[second]].start_row_index = row;
								queue[second][queueSize[second]].end_row_index = row;
								queue[second][queueSize[second]].lowest_right = newLR;
								queueSize[second]++;
							}
						}
						else {
							// WM: think you should be checking row here as well
							if (!outsideQueue && row >= queue[first][qIndex].start_row_index && row <= queue[first][qIndex].end_row_index) {
								if (!(row == queue[first][qIndex].start_row_index && queue[first][qIndex].lowest_right > Rf.end)) {
									double prevR = row == queue[first][qIndex].start_row_index ? queue[first][qIndex].lowest_right : 0.0;
									double newLR = std::max(prevR, Rf.start);
									if (isComplete(Rf) && newLR == 0.0 && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row - 1) {
										// complete reachable right means increase previous queue entry to span this cell
										queue[second][queueSize[second] - 1].end_row_index = row;
									}
									else {
										// push to queue
										queue[second][queueSize[second]].start_row_index = row;
										queue[second][queueSize[second]].end_row_index = row;
										queue[second][queueSize[second]].lowest_right = newLR;
										queueSize[second]++;
									}
								}
							}
						}
					}
					// Top edge stored in Tf, TFree = false means not free
					bool TFree = computeInterval(P[row + 1], Q[column], Q[column + 1], queryDelta, Tf);
					if (!outsideQueue && row <= queue[first][qIndex].end_row_index && row >= queue[first][qIndex].start_row_index) {
						if (row == queue[first][qIndex].end_row_index) {
							// consume the first queue
							qIndex++;
						}
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
					//try and jump
					if (!outsideQueue && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row && Rf.end == 1) {
						// jump-off point possible
						// check if minimum jump distance is big enough
						int gapSize = queue[first][qIndex].end_row_index - queue[first][qIndex].start_row_index;
						if (gapSize > 1) {
							std::vector<Portal> &ports = portals[row];
							choice.source = -1;
							for (Portal &p : ports) {
								int jumpSize = p.destination - p.source;
								// check if jump within range
								if (p.destination <= queue[first][qIndex].end_row_index) {
									// check if jump distance fits
									double segmentFrechet = computeSegmentFrechet(p, column, P, Q);
									// we use adjusted querydelta to avoid false negatives due to epsilons from simplification steps
									if (segmentFrechet + p.distance <= baseQueryDelta) {
										choice = p;
									}
								}
								else {
									// can't reach this jump, jumps are sorted, so break
									break;
								}
							}
							// JUMP!
							if (choice.source != -1) {
								row = choice.destination - 1;// - 1 to counter ++ later
								queue[second][queueSize[second] - 1].end_row_index = row;
							}
						}
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

		// figure out what constitutes success decision and return it
		int endIndex = queueSize[first] - 1;
		if (endIndex == -1) return false;
		bool exit = queue[first][endIndex].start_row_index == size_p - 2 && queue[first][endIndex].lowest_right <= 1;
		return exit || (queue[first][endIndex].end_row_index == size_p - 2 && queue[first][endIndex].start_row_index != size_p - 2);
	}

	// wrapper
	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta, double baseQueryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta, baseQueryDelta, P.simpPortals);
	}

	// wrapper
	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta, queryDelta, P.simpPortals);
	}
};
