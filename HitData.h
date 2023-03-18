#pragma once
#include <iostream>
#include <bitset>
#include <map>
#include <vector>
#include <string>
#include <math.h>

#include "Detector.h"
#include "Tangents.h"


//(x-x1)^2 + (y-y1)^2 == R1^2;
//dx + fy + e = 0
//returns coordinates of tangence point between line & circle1
//(x_tang, y_tang)
pt tangentCrd(double _x1, double _y1, double _r1,
	double _d, double _e, double _f) {
	double coef = 1.0e6; //round to 6th digit, else goes to nan
	double x1 = std::round(_x1 * coef) / coef;
	double y1 = std::round(_y1 * coef) / coef;
	double r1 = std::round(_r1 * coef) / coef;
	double d = std::round(_d * coef) / coef;
	double e = std::round(_e * coef) / coef;
	double f = std::round(_f * coef) / coef;

	double resultX = 0.0;
	double resultY = 0.0;

	//magic equations for x and y 
	//calcualted with Wolfram Mathematica
	//Solve[{(x - x1)^2 + (y - y1)^2 == r1^2, 
	//y == -((d*x + e)/f)}, {x, y}]
	//this is the same as
	//Solve[{(x - a)^2 + (y - b)^2 == r^2, 
	//d*x + f*y + e == 0}, {x, y}]
	//X coordinate
	resultX = (1 / (sqr(d) + sqr(f)) * 0.5 *
	(-2 * d * e - 2 * y1 * d * f + 2 * x1 * sqr(f) -
		sqrt(sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f)) -
			4 * (sqr(d) + sqr(f)) *
			(sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
				sqr(y1) * sqr(f) - sqr(f) * sqr(r1)
				)
			)
		)
	);
	//Y coordinate
	resultY = (1 / f) * 
		(-e + (sqr(d) * e) / (sqr(d) + sqr(f)) +
		(y1 * sqr(d) * f) / (sqr(d) + sqr(f)) -
		(x1 * d * sqr(f)) / (sqr(d) + sqr(f)) +
		1 / (2 * (sqr(d) + sqr(f))) * d *
		sqrt(sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f)) -
			4 * (sqr(d) + sqr(f)) *
			(sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
				sqr(y1) * sqr(f) - sqr(f) * sqr(r1))
			)
		);

	pt result{ resultX, resultY };
	return result;
}


plane TangPlaneToCylinder(pt3d tangPoint, cylinder c1) {
	double coefX = 2 * (tangPoint.x - c1.x);
	double coefY = 2 * (tangPoint.y - c1.y);
	double coefZ = 2 * (tangPoint.z - c1.z);
	double coefD = -coefX * tangPoint.x - coefY * tangPoint.y
		- coefZ * tangPoint.z;
	plane tangPlane = plane(coefX, coefY, coefZ, coefD);
	return tangPlane;
}

pt3d HitPointOnDefinedPlane(plane& pl1, plane& pl2, double& z_plane) {
	pt3d pt{ 0.0,0.0,0.0 };
	return pt;
}


//0b|y22|y21|y12|y11|'|x22|x21|x12|x11|
void GetDetectorHitData(const uint8_t& word ) {
	std::vector<std::string> DetectorNames = {
		"X11", "X12", "X21", "X22", "Y11", "Y12", "Y21", "Y22"
	};
	std::vector<double> isochrones = {
		5.2, 0.0, 2.4, 0.0, 3.2, 0.0, 7.6, 0.0
	};
	std::vector<Detector> detArray;
	for (size_t i = 0b0000'0001, j = 0; i <= 0b1000'0000; i <<= 0b0000'0001, j++) {
		std::cout << (std::bitset<8>)word << " "
			<< (std::bitset<8>)i << " " << ((std::bitset<8>)(word & i)).any() << " " << DetectorNames.at(j) << std::endl;
		Detector det(DetectorNames.at(j), ((std::bitset<8>)(word & i)).any(), 5.32);
		det._isoCylinder = cylinder();
		if (DetectorNames.at(j).at(0) == 'X') det._isoCylinder.orientation = posPlane::XZ;
		else if (DetectorNames.at(j).at(0) == 'Y') det._isoCylinder.orientation = posPlane::YZ;
		else det._isoCylinder.orientation = posPlane::NA;
		detArray.push_back(det);
	}

}