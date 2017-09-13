#pragma once

#include <functional>
#include <algorithm>

// Does binary search on integer range (lowerbound, upperbound),
// accepts lambda function returning whether the given search index
// satisfies the search criterion.
int binaryIntSearch(const std::function< bool(int) >& f, int upperbound, int lowerbound) {
	int rangeLength = upperbound - lowerbound;
	if (rangeLength <= 1) {
		return lowerbound;
	}
	int middle = lowerbound + (rangeLength) / 2;
	bool result = f(middle);
	if (result) {
		return binaryIntSearch(f, upperbound, middle);
	}
	else {
		return binaryIntSearch(f, middle, lowerbound);
	}
}

// does binary search on double range (lowerbound, upperbound),
// accepts lambda function returning whether the given search index
// satisfies the search criterion.
void binaryDoubleSearch(const std::function< int(double) >& f, double upperbound, double lowerbound) {
	double rangeLength = upperbound - lowerbound;
	double avg = lowerbound + (rangeLength) / 2;
	int result = f(avg);
	if (result == 1) {
		binaryDoubleSearch(f, upperbound, avg);
	}
	else if (result == 0){
		binaryDoubleSearch(f, avg, lowerbound);
	}
	else {
		return;
	}
}


// Does double & search on integer range (lowerbound, upperbound),
// accepts lambda function returning whether the given search index
// satisfies the search criterion.
int doubleNsearch(const std::function< bool(int) >& f, int start, int end, int doubleNSearchBase, double doubleNSearchExponentStep) {
	int k = start;
	int prevk = start;
	int iteration = 0;
	while (true) {
		//double
		if (k > end - 1) {
			k = end - 1;
		}
		bool epsValid = f(k);
		if (!epsValid) {
			//binary search
			k = binaryIntSearch(
				f,
				k,
				prevk
				);
			return k;
		}
		else {
			if (k == end - 1) {
				return k;
			}
			prevk = k;
			k += (int)floor(pow(doubleNSearchBase, doubleNSearchExponentStep*iteration));
			iteration++;
		}
	}
}

