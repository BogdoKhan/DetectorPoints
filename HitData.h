#pragma once
#include <iostream>
#include <bitset>
#include <map>
#include <vector>
#include <string>
#include <math.h>

#include "Detector.h"
#include "Tangents.h"

void ShowDetArrayData(std::vector<Detector>& detArray);

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
	pt result = {};
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
	result = { resultX, resultY };

	if (isnan(resultX) || isnan(resultY)) {
		 x1 = std::ceil(_x1 * coef) / coef;
		 y1 = std::ceil(_y1 * coef) / coef;
		 r1 = std::ceil(_r1 * coef) / coef;
		 d = std::ceil(_d * coef) / coef;
		 e = std::ceil(_e * coef) / coef;
		 f = std::ceil(_f * coef) / coef;
		 result = tangentCrd(x1, y1, r1, d, e, f);
		 if (isnan(result.crd1) || isnan(result.crd2)) {
			 x1 = std::floor(_x1 * coef) / coef;
			 y1 = std::floor(_y1 * coef) / coef;
			 r1 = std::floor(_r1 * coef) / coef;
			 d = std::floor(_d * coef) / coef;
			 e = std::floor(_e * coef) / coef;
			 f = std::floor(_f * coef) / coef;
			 result = tangentCrd(x1, y1, r1, d, e, f);
		 }
	}
	return result;
}

pt3d tangentCrd(cylinder cyl,
	line lin) {
	double coef = 1.0e6; //round to 6th digit, else goes to nan
	double x1 = std::round(cyl.crd1 * coef) / coef;
	double y1 = std::round(cyl.crd2 * coef) / coef;
	double r1 = std::round(cyl.r * coef) / coef;
	double d = std::round(lin.a * coef) / coef;
	double e = std::round(lin.b * coef) / coef;
	double f = std::round(lin.c * coef) / coef;

	pt tangPt2d = tangentCrd(x1, y1, r1, d, e, f);

	if (cyl.orientation == posPlane::XZ) return { tangPt2d.crd1, 0.0, tangPt2d.crd2 };
	else if (cyl.orientation == posPlane::YZ) return { 0.0, tangPt2d.crd1, tangPt2d.crd2 };
	pt3d result{ tangPt2d.crd1, 0.0, tangPt2d.crd2 };
	return result;
}

//create equation of plane tangent to 1 cylinder in
//a pre-defined tangence point based on 
//lines tangent to cross-sectional circles
//this plane is tangent to 2 cylinders
plane TangPlaneToCylinder(pt3d tangPoint, cylinder c1) {
	double coefX = 2 * (tangPoint.x - c1.x);
	double coefY = 2 * (tangPoint.y - c1.y);
	double coefZ = 2 * (tangPoint.z - c1.z);
	double coefD = -coefX * tangPoint.x - coefY * tangPoint.y
		- coefZ * tangPoint.z;
	plane tangPlane = plane(coefX, coefY, coefZ, coefD);
	return tangPlane;
}
//coordinates of intersection point between the plane parallel to XY
//and tangent line built on 2 intersecting planes
pt3d HitPointOnDefinedPlane(plane& pl1, plane& pl2, double& z_plane) {
	pt3d pt{ 0.0,0.0,0.0 };
	pt.z = z_plane;
	//for XZ plane: ax + by + cz + d = 0 => b==0; x = -(cz+d)/a
	pt.x = - (pl1.c * pt.z + pl1.d) / pl1.a;
	//for YZ plane: ax + by + cz + d = 0 => b==0; y = -(cz+d)/b
	pt.y = -(pl2.c * pt.z + pl2.d) / pl2.b;
	return pt;
}

double SetIsochroneRadius(double& isoRadius, const double& tubeRadius) {
	if (isnan(isoRadius) || isoRadius >= tubeRadius || isoRadius <= 0) return tubeRadius;
	return isoRadius;
}

//0b|y22|y21|y12|y11|'|x22|x21|x12|x11|
//here we create the grid for our detectors field
void GetDetectorHitData(const uint8_t& word, 
	std::vector<double>& isochrones, 
	std::vector<Detector>& detArray,
	const double& tubeRadius ) {
	std::vector<std::string> DetectorNames = {
		"X11", "X12", "X21", "X22", "Y11", "Y12", "Y21", "Y22"
	};

	//std::vector<Detector> detArray;
	for (size_t i = 0b0000'0001, j = 0; i <= 0b1000'0000; i <<= 0b0000'0001, j++) {
		std::cout << (std::bitset<8>)word << " "
			<< (std::bitset<8>)i << " " << ((std::bitset<8>)(word & i)).any() << " " << DetectorNames.at(j) << std::endl;

		Detector det(DetectorNames.at(j), ((std::bitset<8>)(word & i)).any()); //create detector: its emplacement(X,Y); was it hit? (counts 1s in binary corresp to its mask);
		det._isoCylinder = cylinder();
		det.DetectorRadius = tubeRadius;

		if (DetectorNames.at(j).at(0) == 'X') {
			det._isoCylinder.orientation = posPlane::XZ;

			det._isoCylinder.x = j / 2 * tubeRadius + j % 2 * 2 * tubeRadius;
			det._isoCylinder.crd1 = det._isoCylinder.x;
			//j = 0, 1, 2, 3, 4, 5, 6, 7
			//x = 0 +0 = 0 || 0 + 1*2 = 2 || 1 + 0 = 1 || 1 + 2 = 3...
			det._isoCylinder.z = j / 2 * tubeRadius;
			if (j / 2 < 1) det._isoCylinder.z;
			else det._isoCylinder.z;
			det._isoCylinder.crd2 = det._isoCylinder.z;

			det._isoCylinder.y = 0.0;
			det._isoCylinder.r = SetIsochroneRadius(isochrones.at(j), tubeRadius);
		}

		else if (DetectorNames.at(j).at(0) == 'Y') {
			det._isoCylinder.orientation = posPlane::YZ;

			det._isoCylinder.x = 0.0;
			det._isoCylinder.y = (j - 4) / 2 * tubeRadius + (j - 4) % 2 * 2 * tubeRadius;

			det._isoCylinder.crd1 = det._isoCylinder.y;
			det._isoCylinder.crd2 = det._isoCylinder.z;

			det._isoCylinder.z = j / 2 * tubeRadius;
			det._isoCylinder.r = SetIsochroneRadius(isochrones.at(j), tubeRadius);
		}
		else det._isoCylinder.orientation = posPlane::NA;

		detArray.push_back(det);
	}
	ShowDetArrayData(detArray);
}

void ShowDetArrayData(std::vector<Detector>& detArray) {
	int counter = 0;
	for (const auto& item : detArray) {
		std::cout << "Detector placed at ";
		if (item._isoCylinder.orientation == posPlane::XZ) { std::cout << "XZ plane "; }
		else if (item._isoCylinder.orientation == posPlane::YZ) { std::cout << "YZ plane "; }
		else { std::cout << "nowhere "; }
		std::cout << "at X= " << item._isoCylinder.x << " Y= " << item._isoCylinder.y << " Z= " << item._isoCylinder.z;
		std::cout << " with r = " << item._isoCylinder.r;
		if (item.isHit()) { std::cout << ", hit"; }
		else std::cout << ", no hit";
		std::cout << "\n";
	}
}

vector<pt3d> GetHitPoints(std::vector<Detector>& detArray) {
	vector<pt3d> hitPoints = {};
	vector<Detector> workingDets = {};
	map<posPlane, vector<plane>> tangPlanesAllDir = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };
	double z_planePosition = 3.0 * detArray.at(0).DetectorRadius;
	for (Detector& det : detArray) {
		if (det.isHit()) workingDets.push_back(det);
	}
	if (workingDets.size() >= 4) {
		if (workingDets.at(0)._isoCylinder.orientation == posPlane::XZ && workingDets.at(1)._isoCylinder.orientation == posPlane::XZ) {
			vector<line> tangLines1 = tangents(workingDets.at(0)._isoCylinder, workingDets.at(1)._isoCylinder);
			for (line& line_ : tangLines1) {
				pt3d tangPt = tangentCrd(workingDets.at(0)._isoCylinder, line_);
				plane tangPlane = TangPlaneToCylinder(tangPt, workingDets.at(0)._isoCylinder);
				tangPlanesAllDir.at(posPlane::XZ).push_back(tangPlane);
			}
		}
		if (workingDets.at(2)._isoCylinder.orientation == posPlane::YZ && workingDets.at(3)._isoCylinder.orientation == posPlane::YZ) {
			vector<line> tangLines1 = tangents(workingDets.at(2)._isoCylinder, workingDets.at(3)._isoCylinder);
			for (line& line_ : tangLines1) {
				pt3d tangPt = tangentCrd(workingDets.at(2)._isoCylinder, line_);
				plane tangPlane = TangPlaneToCylinder(tangPt, workingDets.at(2)._isoCylinder);
				tangPlanesAllDir.at(posPlane::YZ).push_back(tangPlane);
			}
		}
	}
	for (plane& plXZ : tangPlanesAllDir.at(posPlane::XZ)) {
		for (plane& plYZ : tangPlanesAllDir.at(posPlane::YZ)) {
			hitPoints.push_back(HitPointOnDefinedPlane(plXZ, plYZ, z_planePosition));
		}
	}
	return hitPoints;
}