#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include "Tangents.h"



class Detector {
public:
	Detector();
	Detector(std::string emplacement, bool isHit);
	cylinder _isoCylinder = cylinder();
	double DetectorRadius = 1.0;
	bool isHit() const;
private:
	std::string _emplacement = "NULL";
	bool _isHit = false;
	double _centerPosition = 0.0;

};
