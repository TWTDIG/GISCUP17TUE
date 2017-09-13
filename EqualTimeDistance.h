
#include "Vertex.h"
#include <vector>
#include <math.h>
#include <algorithm>

#pragma once



// Implements Equal Time Distance algorithm between two trajectories.
// The ETD algorithm computes an approximation of frechet distance by
// taking the 'dog leash' length when traversing two trajectories at 
// the same speed. Used by agarwal and simplification step.
// If found relevant, this function can be optimized using SIMD instructions
static double equalTimeDistance(
	std::vector<Vertex> &pverts, std::vector<Vertex> &qverts,
	std::vector<double> &ptotals, std::vector<double> &qtotals,
	std::vector<double> &pdistances, std::vector<double> &qdistances,
	int psize, int qsize,
	int pstart, int qstart) {

	double pdistOffset = ptotals[pstart];
	double qdistOffset = qtotals[qstart];
	double pdist = ptotals[psize - 1] - pdistOffset;
	double qdist = qtotals[qsize - 1] - qdistOffset;
	double pscale = qdist / pdist;
	int p_ptr = pstart + 1;
	int q_ptr = qstart + 1;

	//startpoints
	double dx = pverts[pstart].x - qverts[qstart].x;
	double dy = pverts[pstart].y - qverts[qstart].y;
	double smax = dx*dx + dy*dy;
	dx = pverts[psize - 1].x - qverts[qsize - 1].x;
	dy = pverts[psize - 1].y - qverts[qsize - 1].y;
	double emax = dx*dx + dy*dy;
	if (qdist == 0 || pdist == 0) return sqrt(std::max(emax, smax));
	double position = 0; // from 0 to 1

	Vertex p_pt;
	Vertex q_pt;

	// start traversing diagonal
	while (
			!(
				p_ptr == psize - 1
				&&
				q_ptr == qsize - 1
			)
		){
		// figure out which cell edge we hit next on the diagonal, vertical or horizontal edge
		double posP = position * pdist;
		double posQ = position * qdist;
		double nextDistP = ptotals[p_ptr] - pdistOffset - posP;
		double nextDistQ = qtotals[q_ptr] - qdistOffset - posQ;

		// safety net to ensure this function terminates if there's still a bug in there
		if (p_ptr == psize - 1) nextDistP = DBL_MAX;
		if (q_ptr == qsize - 1) nextDistQ = DBL_MAX;



		if (nextDistP * pscale < nextDistQ) { // treat P first
			p_pt.x = pverts[p_ptr].x;
			p_pt.y = pverts[p_ptr].y;
			position = (ptotals[p_ptr] - pdistOffset) / pdist;
			double scale = (position * qdist - (qtotals[q_ptr - 1] - qdistOffset)) / qdistances[q_ptr];
			dx = qverts[q_ptr].x - qverts[q_ptr - 1].x;
			dy = qverts[q_ptr].y - qverts[q_ptr - 1].y;
			q_pt.x = qverts[q_ptr - 1].x + dx * scale;
			q_pt.y = qverts[q_ptr - 1].y + dy * scale;
			p_ptr++;
		}
		else { // treat Q first
			q_pt.x = qverts[q_ptr].x;
			q_pt.y = qverts[q_ptr].y;
			position = (qtotals[q_ptr] - qdistOffset) / qdist;
			double scale = (position * pdist - (ptotals[p_ptr - 1] - pdistOffset)) / pdistances[p_ptr];
			dx = pverts[p_ptr].x - pverts[p_ptr - 1].x;
			dy = pverts[p_ptr].y - pverts[p_ptr - 1].y;
			p_pt.x = pverts[p_ptr - 1].x + dx * scale;
			p_pt.y = pverts[p_ptr - 1].y + dy * scale;
			q_ptr++;
		}

		// figure out distance we need for this point on the diagonal
		dx = p_pt.x - q_pt.x;
		dy = p_pt.y - q_pt.y;
		double nm = dx*dx + dy*dy;
		// overwrite current frechet if dist is larger
		// skipping sqrts for speed
		if (nm > smax) {
			smax = nm;
		}
	}

	//endpoints
	dx = pverts[p_ptr].x - qverts[q_ptr].x;
	dy = pverts[p_ptr].y - qverts[q_ptr].y;
	double nm = dx*dx + dy*dy;
	if (nm > smax) {
		smax = nm;
	}
	// finally sqrt whatever max is found
	return sqrt(smax);
}

double equalTimeDistance(Trajectory &p, Trajectory &q) {
	return equalTimeDistance(p.vertices, q.vertices, p.totals, q.totals, p.distances, q.distances, p.size, q.size, 0, 0);
}
