#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {
	std::vector<double> isochrones = {
		4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0
	};
	std::vector<Detector> detArray;
	//uint8_t word = 0b0000'0000;
	//GetWorkingDetectors(word);

	GetDetectorHitData(0b01010101, isochrones, detArray, 5.0);
	vector <pt3d> hitPoints = GetHitPoints(detArray);
	////pt3d hitPointsErr = {2*10.0/sqrt(12), 2 * 10.0 / sqrt(12) , 2 * 10.0 / sqrt(12) };
	size_t counter = 1;
	for (pt3d& i : hitPoints) {
		cout << "Point: " << counter << " {" 
			<< i.x << 
			", " << i.y  <<
			", " << i.z << "} \n";
		counter++;
	}

	return 0;
}