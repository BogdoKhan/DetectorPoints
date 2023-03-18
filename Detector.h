#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include "Tangents.h"

const double DetectorRadius = 1.0;

class Detector {
public:
	Detector();
	Detector(std::string emplacement, bool isHit,
		double isochroneRadius);
	cylinder _isoCylinder = cylinder();
private:
	std::string _emplacement = "NULL";
	bool _isHit = false;
	double _isochroneRadius = 0.0;
	double _centerPosition = 0.0;

};

Detector::Detector() {
	
};

Detector::Detector(std::string emplacement, bool isHit,
	double isochroneRadius) :
	_emplacement(emplacement),
	_isHit(isHit),
	_isochroneRadius(isochroneRadius)
{
	std::cout << "Detector at " << _emplacement;
	if (_isHit) {
		std::cout << " with a hit and isochrone with a radius of " << std::setprecision(3) << _isochroneRadius;
	}
	else {
		std::cout << " with no hit";
	}
	std::cout << " is created\n";
};