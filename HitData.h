#pragma once
#include <iostream>
#include <bitset>
#include <map>
#include <vector>
#include <string>
#include <math.h>

#include "Detector.h"
#include "Tangents.h"
#include "RandGen.h"
#include "GSL_Solvers.h"

#define PI 3.14159265

void ShowDetArrayData(std::vector<Detector>& detArray);

pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
	const double& _d, const double& _e, const double& _f, bool flag_round);

pt3d tangentCrd(const cylinder& cyl,
	const line& lin);

plane TangPlaneToCylinder(const pt3d& tangPoint, const cylinder& c1);

pt3d HitPointOnDefinedPlane(const plane& pl1, const plane& pl2, const double& z_plane);

double SetIsochroneRadius(double& isoRadius, const double& tubeRadius);

void GetDetectorHitData(
	std::vector<double>& isochrones,
	std::vector<Detector>& detArray,
	const double& tubeRadius);

void ShowDetArrayData(std::vector<Detector>& detArray);

vector<pt3d> GetHitPoints(std::vector<Detector>& detArray);

pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
	const double& _d, const double& _e, const double& _f, bool flag_round) {
	const gsl_multiroot_fsolver_type* T; //solver type
	gsl_multiroot_fsolver* s;			 //solver as object

	int status;							//GSL_CONTINUE or GSL_SUCCESS
	size_t iter = 0;
	const size_t n = 2;

	struct params_tangentCrd p =
	{ _x1, _y1, _r1, _d, _e, _f }; // structure containing parameters

	gsl_multiroot_function f;
	f.f = &circleAndLine;	//our function
	f.n = n;				//max dimensions
	f.params = &p;			//structure containing parameters

	double x_init[2] = { _x1 + _r1, _y1 + _r1 }; //non-zero initial conditions to prevent algorithm from being stuck

	gsl_vector* x = gsl_vector_alloc(n);		//initial guess for variables
	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);

	T = gsl_multiroot_fsolver_hybrids;				//solver type
	s = gsl_multiroot_fsolver_alloc(T, 2);
	gsl_multiroot_fsolver_set(s, &f, x);			//set solver
	print_state(iter, s);
	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(iter, s);
		if (status) /* check if solver is stuck */
			break;
		status =
			gsl_multiroot_test_residual(s->f, 1e-7);
	} while (status == GSL_CONTINUE && iter < 1000);

	printf("status = %s\n", gsl_strerror(status));
	std::cout << "crd1: " << gsl_vector_get(s->x, 0) << " crd2: "
		<< gsl_vector_get(s->x, 1) << std::endl;

	return { gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1) };
}

pt3d tangentCrd(const cylinder& cyl,
	const line& lin) {
	double coef = 1.0e6; //round to 6th digit, else goes to nan
	double x1 = std::round(cyl.crd1 * coef) / coef;
	double y1 = std::round(cyl.crd2 * coef) / coef;
	double r1 = std::round(cyl.r * coef) / coef;
	double d = std::round(lin.a * coef) / coef;
	double e = std::round(lin.b * coef) / coef; //fixed swapped values
	double f = std::round(lin.c * coef) / coef; 

	pt tangPt2d = tangentCrd(x1, y1, r1, d, e, f, false);

	if (cyl.orientation == posPlane::XZ) return { tangPt2d.crd1, 0.0, tangPt2d.crd2 };
	else if (cyl.orientation == posPlane::YZ) return { 0.0, tangPt2d.crd1, tangPt2d.crd2 };
	pt3d result{ tangPt2d.crd1, 0.0, tangPt2d.crd2 };
	return result;
}

//create equation of plane tangent to 1 cylinder in
//a pre-defined tangence point based on 
//lines tangent to cross-sectional circles
//this plane is tangent to 2 cylinders
plane TangPlaneToCylinder(const pt3d& tangPoint, const cylinder& c1) {
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
pt3d HitPointOnDefinedPlane(const plane& pl1, const plane& pl2, const double& z_plane) {
	pt3d pt{ 0.0,0.0,0.0 };
	pt.z = z_plane;
	//for XZ plane: ax + by + cz + d = 0 => b==0; x = -(cz+d)/a
	pt.x = -(pl1.c * pt.z + pl1.d) / pl1.a;
	//for YZ plane: ax + by + cz + d = 0 => b==0; y = -(cz+d)/b
	pt.y = -(pl2.c * pt.z + pl2.d) / pl2.b;
	return pt;
}
//set constraints on isochrone radius, it must be above 0 and below tube radius
double SetIsochroneRadius(double& isoRadius, const double& tubeRadius) {
	if (isnan(isoRadius) || isoRadius >= tubeRadius || isoRadius <= 0) return tubeRadius;
	return isoRadius;
}
//form binary word from isochrones vector
//0b|y22|y21|y12|y11|'|x22|x21|x12|x11|
//starts from 0b0000'0001, corresponding to x11
//please note that x11 corresponds to the first element of isochrones vector
uint8_t GetWordFromIsochrones(const std::vector<double>& isochrones) {
	uint8_t word = 0b0000'0000;
	for (size_t i = 0b0000'0001, j = 0; i <= 0b1000'0000; i <<= 0b0000'0001, j++) {
		if (isochrones.at(j) != 0.0 && !isnan(isochrones.at(j))) {
			word += i;
		}
	}
	return word;
}

//here we create the grid for our detectors field
void GetDetectorHitData(
	std::vector<double>& isochrones,
	std::vector<Detector>& detArray,
	const double& tubeRadius) {
	std::vector<std::string> DetectorNames = {
		"X11", "X12", "X21", "X22", "Y11", "Y12", "Y21", "Y22"
	};
	const uint8_t word = GetWordFromIsochrones(isochrones);
	//std::vector<Detector> detArray;
	for (size_t i = 0b0000'0001, j = 0; i <= 0b1000'0000; i <<= 0b0000'0001, j++) {
		std::cout << (std::bitset<8>)word << " "
			<< (std::bitset<8>)i << " " << ((std::bitset<8>)(word & i)).any() << " " << DetectorNames.at(j) << std::endl; //status info

		Detector det(DetectorNames.at(j), ((std::bitset<8>)(word & i)).any()); //create detector: its emplacement(X,Y); was it hit? (counts 1s in binary corresp to its mask);
		det._isoCylinder = cylinder();											//create cylinder as base of detector
		det.DetectorRadius = tubeRadius;										//detector tube radius is defined now
		det._isoCylinder.circle::orientation = det._isoCylinder.orientation;	//set orientation of detector and place it on the grid

		if (DetectorNames.at(j).at(0) == 'X') {									//if name starts with "X"
			det._isoCylinder.orientation = posPlane::XZ;

			//j = 0, 1, 2, 3, 4, 5, 6, 7
			//x = 0 +0 = 0 || 0 + 1*2 = 2 || 1 + 0 = 1 || 1 + 2 = 3...
			det._isoCylinder.x = floor(j / 2) * tubeRadius + j % 2 * 2 * tubeRadius;	//0; 10; 5; 15
			det._isoCylinder.y = 0.0;
			det._isoCylinder.z = floor(j / 2) * 2 * tubeRadius * cos(30 * PI / 180);

			det._isoCylinder.circle::pt::crd1 = det._isoCylinder.x;
			det._isoCylinder.circle::pt::crd2 = det._isoCylinder.z;
			det._isoCylinder.r = SetIsochroneRadius(isochrones.at(j), tubeRadius);
		}

		else if (DetectorNames.at(j).at(0) == 'Y') {
			det._isoCylinder.orientation = posPlane::YZ;

			det._isoCylinder.x = 0.0;
			det._isoCylinder.y = floor((j - 4) / 2) * tubeRadius + (j - 4) % 2 * 2 * tubeRadius; //0;10;5;15
			det._isoCylinder.z = 2 * tubeRadius + floor((j - 2) / 2) * 2 * tubeRadius * cos(30 * PI / 180);

			det._isoCylinder.circle::pt::crd1 = det._isoCylinder.y;
			det._isoCylinder.circle::pt::crd2 = det._isoCylinder.z;
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
	double z_planePosition = detArray.at(0).DetectorRadius * (1 + 2*cos(30*PI/180)); //R + 2Rcos(30*) is the position of plane dividing XZ and YZ detectors
	//define 4 working detectors
	for (const Detector& det : detArray) {
		if (det.isHit()) workingDets.push_back(det);
	}
	//check if all 4 working detectors for each layer are chosen propely
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
	//get coordinates of hit points at the plane between XZ and YZ detectors
	for (plane& plXZ : tangPlanesAllDir.at(posPlane::XZ)) {
		for (plane& plYZ : tangPlanesAllDir.at(posPlane::YZ)) {
			hitPoints.push_back(HitPointOnDefinedPlane(plXZ, plYZ, z_planePosition));
		}
	}
	return hitPoints;
}

void GetWorkingDetectors(uint8_t word) {
	if (((std::bitset<8>)(~word & 0b1111'0000)).count() == 4) {

	}
	if (((std::bitset<8>)(~word & 0b0000'1111)).count() == 4) {

	}
	//if (((std::bitset<8>)(word & 0b0000)).count() > 0) {
		std::cout << ((std::bitset<8>)(~word)) << " " <<
			((std::bitset<8>)(0b1111'0000)) << " "
			<< ((std::bitset<8>)(~word & 0b1111'0000)).count()
			<< " "
			<< ((std::bitset<8>)(~word & 0b0000'1111)).count() << std::endl;
	//}

	for (size_t i = 0b1100'0000; i >= 0b0000'0010; i >>= 0b0000'0010) {
		//if 4 in row are without hit, skip

		if ((((std::bitset<8>)(word & i)).count()) != 1) {
			if (randNum()) {
				std::cout << ((std::bitset<8>)(i & (i << 0b0000'0001))) << std::endl;
			}
			else {
				std::cout << ((std::bitset<8>)(i & (i >> 0b0000'0001))) << std::endl;
			}
		}
		std::cout << (std::bitset<8>)word << " "
			<< (std::bitset<8>)i << " " << ((std::bitset<8>)(word & i)).count() << " " << randNum() << std::endl;
	}
}