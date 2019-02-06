#pragma once

#include "Vertex.h"
#include "BoundingBox.h"

#include <stdio.h>
#include <vector>
#include <unordered_set>

// Adapted from implementation of Yago Diez, could be sped up but is not a bottleneck
class DiHash
{
public:

	int slotsPerDimension; //number of equally spaced subdivisions in each dimension
	double limits[2][2]; // two rows (x,y) and two columns(min,max), keeps the information on limits in each dimension, for 3D, just add one row.

	double tol; // tolerance to prevent numerical representation errors

	std::vector<Vertex>** elements;

	DiHash(BoundingBox &boundingBox, int numC, double iTol) {

		slotsPerDimension = numC;
		tol = iTol;

		// First, locate numerichal limits

		// initialize limits at a three-component two doubles (max and min) vector

		// set initial values for extreme values
		limits[0][0] = boundingBox.minx; //min in X
		limits[0][1] = boundingBox.maxx; //max in X
		limits[1][0] = boundingBox.miny; //min in Y
		limits[1][1] = boundingBox.maxy; //max in Y

		// Now that we have extreme dimensions, set points in its corresponding part of the grid

		elements = new std::vector<Vertex>*[slotsPerDimension];
		for (int i = 0; i < slotsPerDimension; i++) {
			elements[i] = new std::vector<Vertex>[slotsPerDimension];
		}
	}


	void addPoint(Vertex pActual) {
		int x, y;

		x = findSlot(pActual.x, 'x', false);
		y = findSlot(pActual.y, 'y', false);

		elements[x][y].push_back(pActual);

	}

	int findSlot(double val, char type, bool allowOverflow) {
		double min, max;
		int retorn;

		//cout<<"DiHash::findSlot limits "<<endl<<"x: ("<<limits[0][0]<<" , "<<limits[0][1]<<")"<<endl;
		//cout<<"y: ("<<limits[1][0]<<" , "<<limits[1][1]<<")"<<endl;

		switch (type) {
		case 'x':
			//cout<<"DiHash::findSlot x "<<endl;
			min = limits[0][0];
			max = limits[0][1];
			break;
		case 'y':
			//cout<<"DiHash::findSlot y"<<endl;
			min = limits[1][0];
			max = limits[1][1];
			break;
		default:
			std::cout << "DiHash::findSlot(double val, char type) wrong slot ";
			break;
		}
		// check that we are not exactly at the maximum or minimum value
        if (fabs(min - val) < tol) retorn = 0;
        else if (fabs(max - val) < tol) retorn = slotsPerDimension - 1;
		else {
			double pas = (fabs(max - min) / slotsPerDimension);

			//cout<<"DiHash::findSlot computing "normal" return value Math.fabs(max-min)"<<Math.fabs(max-min)<<" pas "<<pas<<endl;

			retorn = (int)((val - min) / pas);
		}

		//cout<<"DiHash::findSlot return valure before checking strange cases (if asked to): "<<retorn<<endl;

		if ((retorn >= slotsPerDimension) || (retorn < 0)) {
			if (!allowOverflow) {
				std::cout << "DiHash::findSlot(double val, char type) wrong return value ";
			}
			else    // in this case, used for searches, we set values out of the values outside the extreme to the last slot
			{
				if (retorn >= slotsPerDimension) retorn = slotsPerDimension - 1;
				else retorn = 0;
			}
		}


		return retorn;

	}
	// Given a type of search (x or y) and a slot number c, return slot c's lower bound
	double slotLowerBound(int c, char type) {
		double min, max;

		switch (type) {
		case 'x':
			min = limits[0][0];
			max = limits[0][1];
			break;
		case 'y':
			min = limits[1][0];
			max = limits[1][1];
			break;
		default:
			std::cout << "DiHash::slotLowerBound wrong slot";
		}

		double pas = (fabs(max - min) / slotsPerDimension);

		return (min + pas * c);
	}

	// Given a type of search (x or y) and a search range (min, max(, return the two indexes of the first and last slots affected by the search
	int* slotsTouched(double min, double max, char type) {
		// first position of the following vector is for the minimum slot affected and second for the maximum
		int* retorn = new int[2];

		retorn[0] = findSlot(min, type, true);
		retorn[1] = findSlot(max, type, true);

		return retorn;
	}

	// Return neighbors at range strictly less than distance for a given point, does not return the query point if present
	void neighbors(Vertex &p, double eps, std::unordered_set<int> &retorn) {
		int* limitsX = slotsTouched(p.x - eps, p.x + eps, 'x');
		int* limitsY = slotsTouched(p.y - eps, p.y + eps, 'y');

		for (int i = limitsX[0]; i <= limitsX[1]; i++) {
			for (int j = limitsY[0]; j <= limitsY[1]; j++) {
				for (Vertex &pActual : elements[i][j]) {
					if (pActual.isStart == p.isStart) {
						double dx = p.x - pActual.x;
						double dy = p.y - pActual.y;
						double distSQ = dx * dx + dy * dy;
						if (distSQ < eps * eps) {
							retorn.insert(pActual.trajectoryNumber);
						}
					}
				}

			}
		}
		delete[] limitsX;
		delete[] limitsY;
	}

	// Emits neighbors at range strictly less than distance for a given point, does not return the query point if present
	// also checks endpoints directly
	void neighborsWithCallback(Vertex &p, Vertex& end, double eps, std::vector<Trajectory*> &trajectories, const std::function< void(Trajectory*) >& emit) {
		int* limitsX = slotsTouched(p.x - eps, p.x + eps, 'x');
		int* limitsY = slotsTouched(p.y - eps, p.y + eps, 'y');

		for (int i = limitsX[0]; i <= limitsX[1]; i++) {
			for (int j = limitsY[0]; j <= limitsY[1]; j++) {
				for (Vertex &pActual : elements[i][j]) {
					if (pActual.isStart == p.isStart) {
						double dx = p.x - pActual.x;
						double dy = p.y - pActual.y;
						double distSQ = dx * dx + dy * dy;
						if (distSQ < eps * eps) {
							Trajectory *t = trajectories[pActual.trajectoryNumber];
							Vertex &pend = t->vertices[t->size - 1];
							double dex = end.x - pend.x;
							double dey = end.y - pend.y;
							double disteSQ = dex * dex + dey * dey;
							if (disteSQ < eps * eps) {
								emit(t);
							}
						}
					}
				}

			}
		}
		delete[] limitsX;
		delete[] limitsY;
	}
	


	~DiHash();
};

