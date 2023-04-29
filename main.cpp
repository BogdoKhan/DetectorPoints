#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {
	//
	std::vector<double> isochrones = {
		4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0
	};
	std::vector<Detector> detArray;
	//uint8_t word = 0b0000'0000;
	//GetWorkingDetectors(word);

	GetDetectorHitData(isochrones, detArray, 5.0);
	vector <pt3d> hitPoints = GetHitPoints(detArray);
	////pt3d hitPointsErr = {2*10.0/sqrt(12), 2 * 10.0 / sqrt(12) , 2 * 10.0 / sqrt(12) };
	size_t counter = 1;
	cout << "Hit points, old counter: \n";
	for (pt3d& i : hitPoints) {
		cout << "Point: " << counter << " {" 
			<< i.x << 
			", " << i.y  <<
			", " << i.z << "} \n";
		counter++;
	}

	vector <pt3d> hitPointsByLines = GetHitPointsByLines(detArray);
	////pt3d hitPointsErr = {2*10.0/sqrt(12), 2 * 10.0 / sqrt(12) , 2 * 10.0 / sqrt(12) };
	counter = 1;
	cout << "Hit points, new counter: \n";
	for (pt3d& i : hitPointsByLines) {
		cout << "Point: " << counter << " {"
			<< i.x <<
			", " << i.y <<
			", " << i.z << "} \n";
		counter++;
	}
	//std::vector<double> iso = {
	//	0.0, 0.0, 0.0, 0.0, 3.2, 0.0, 1.6, 0.0
	//};
	//uint8_t word = GetWordFromIsochrones(iso);
	//std::cout << (std::bitset<8>)word << std::endl;

	return 0;
}