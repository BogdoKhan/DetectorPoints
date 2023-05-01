#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {
	//std::vector<double> isochrones = {4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0};
	//std::vector<double> isochrones = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	std::vector<double> isochrones = { 4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0 };
	std::vector<Detector> detArray;

	GetDetectorHitData(isochrones, detArray, 5.0);
	size_t counter = 1;
	vector <pt3d> hitPointsByLines = GetHitPointsByLines(detArray);
	////pt3d hitPointsErr = {2*10.0/sqrt(12), 2 * 10.0 / sqrt(12) , 2 * 10.0 / sqrt(12) };
	counter = 1;
	cout << "Hit points, new counter: \n";
	for (pt3d& i : hitPointsByLines) {
		cout << "No. " << counter << " " << i << endl;
		counter++;
	}

	return 0;
}