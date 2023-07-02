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

pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
	const double& _d, const double& _e, const double& _f, bool showOutput);

pt3d tangentCrd(const cylinder& cyl,
	const line& lin, bool showOutput);

plane TangPlaneToCylinder(const pt3d& tangPoint, const cylinder& c1);

pt3d HitPointOnDefinedPlane(const plane& pl1, const plane& pl2, const double& z_plane);

double SetIsochroneRadius(double& isoRadius, const double& tubeRadius);

uint8_t GetWordFromIsochrones(const std::vector<double>& isochrones);

void GetDetectorHitData(
	std::vector<double>& isochrones,
	std::vector<Detector>& detArray,
	const double& tubeRadius);

void ShowDetArrayData(std::vector<Detector>& detArray);

double geomVectorLength(const pt3d& pt1, const pt3d& pt2);

double getIncidenceAngle(const pt3d& intersectPt, const pt3d& tangPt, const pt3d& tangPtPerp);

void GetIntersectionPoint_case1(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	const double& z_planePosition);

void GetIntersectionPoint_case2(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	const double& z_planePosition);

vector<pt3d> GetHitPointsByLines(std::vector<Detector>& detArray);


pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
	const double& _d, const double& _e, const double& _f, bool showOutput) {
	const gsl_multiroot_fsolver_type* T; //solver type
	gsl_multiroot_fsolver* s;			 //solver as object

	int status;							//GSL_CONTINUE or GSL_SUCCESS
	size_t iter = 0;
	const size_t n = 2;

	struct params_tangentCrd p =
	{ _x1, _y1, _r1, _d, _e, _f }; // structure containing parameters

	gsl_multiroot_function function;
	function.f = &circleAndLine;	//our function
	function.n = n;				//max dimensions
	function.params = &p;			//structure containing parameters

	double x_init[2] = { _x1 + _r1, _y1 + _r1 }; //non-zero initial conditions to prevent algorithm from being stuck

	gsl_vector* x = gsl_vector_alloc(n);		//initial guess for variables
	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);

	T = gsl_multiroot_fsolver_hybrids;				//solver type
	s = gsl_multiroot_fsolver_alloc(T, 2);
	gsl_multiroot_fsolver_set(s, &function, x);			//set solver
	if (showOutput) print_state(iter, s);
	do
	{
		iter++;										//iterate while required tolerance is not reached
		status = gsl_multiroot_fsolver_iterate(s);
		if (showOutput) print_state(iter, s);
		if (status)									// check if solver is stuck
			break;
		status =
			gsl_multiroot_test_residual(s->f, 1e-6);
	} while (status == GSL_CONTINUE && iter < 1000);

	if (showOutput) {
		printf("status = %s\n", gsl_strerror(status));
		std::cout << "crd1: " << gsl_vector_get(s->x, 0) << " crd2: "
			<< gsl_vector_get(s->x, 1) << std::endl;
	}
	if (status == GSL_SUCCESS || (gsl_vector_get(s->f, 0) < 5e-4 &&
		(gsl_vector_get(s->f, 1) < 5e-4))) {
		if (showOutput) std::cout << gsl_vector_get(s->x, 0) << " " << gsl_vector_get(s->x, 1) << std::endl;
	}
	else {
		if (showOutput) std::cout << "Roots not found\n";
		return { NAN, NAN };
	}

	return { gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1) };
}


pt3d tangentCrd(const cylinder& cyl,
	const line& lin, bool showOutput) {
	double coef = 1.0e6; //round to 6th digit, else goes to nan
	double x1 = std::round(cyl.crd1 * coef) / coef;
	double y1 = std::round(cyl.crd2 * coef) / coef;
	double r1 = std::round(cyl.r * coef) / coef;
	double d = std::round(lin.a * coef) / coef;
	double e = std::round(lin.b * coef) / coef; //fixed swapped values
	double f = std::round(lin.c * coef) / coef;

	pt tangPt2d = tangentCrd(x1, y1, r1, d, e, f, showOutput);

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


double geomVectorLength(const pt3d& pt1, const pt3d& pt2) {
	return sqrt(gsl_pow_2(pt2.x - pt1.x) + gsl_pow_2(pt2.y - pt1.y) + gsl_pow_2(pt2.z - pt1.z));
}

double getIncidenceAngle(const pt3d& intersectPt, const pt3d& tangPt, const pt3d& tangPtPerp) {
	double tangentLength = geomVectorLength(tangPt, intersectPt);
	double tangPerpLength = geomVectorLength(tangPt, tangPtPerp);
	double incidenceAngle = (asin(tangPerpLength / tangentLength) / PI * 180);
	std::cout << "Angle is: " << incidenceAngle << " degrees\n";
	return incidenceAngle;
}

void GetIntersectionPoint_case1(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	const double& z_planePosition) {
	size_t i_start = 0;
	if (planePosition == posPlane::XZ) i_start = 0;
	else if (planePosition == posPlane::YZ) i_start = 4;
	size_t i_end = i_start + 4;

	double thresholdAngle = 30.0;

	map<size_t, vector<line>> tangentLines;
	map<posPlane, vector<pt3d>> tangencyPt;
	map<posPlane, vector<pt3d>> intersectionPtCandidates;
	for (size_t i = i_start; i < i_end; i++) {
		if (detArray.at(i).isHit()) {
			tangentLines[i] = {};						//skips the calculation for the only detector hit, we use it as anchor to find tangents to non-active detectors
			continue;
		}
		else {
			tangentLines[i] = tangents(workingDets.at(planePosition).at(0)._isoCylinder, detArray.at(i)._isoCylinder);	//iterate over detectors array to build tangents between working one and non-active dets
			for (line& line_ : tangentLines.at(i)) {
				double intersect_crd = -(line_.b * z_planePosition + line_.c) / line_.a;
				pt3d intersectPt{ 0.0, 0.0, z_planePosition };
				if (planePosition == posPlane::XZ) intersectPt.x = intersect_crd;
				else if (planePosition == posPlane::YZ) intersectPt.y = intersect_crd;			//X or Y intersection coordinates

				pt3d tangPt = tangentCrd(workingDets.at(planePosition).at(0)._isoCylinder, line_, true); //tangence point
				pt3d tangPtPerp = tangPt;
				tangPtPerp.z = z_planePosition;															//point coming from tangence point perpendicular to Z plane
				double incidenceAngle = getIncidenceAngle(intersectPt, tangPt, tangPtPerp);				//calculate angle of incidence wrt Z axis
				
				if (incidenceAngle < thresholdAngle) continue;
				//place here the implementation for the line tangent @ specified angle

				bool noIntersectWithNonActiveDets = false;												//flag: if the tangent intersects any non-hit detector, which it is not tangent to, reject the line
				for (size_t item = i_start; item < i_end; item++) {
					if (detArray.at(item).isHit() || item == i) continue;								//skip if crosses the pair of detectors to which it is tangent
					noIntersectWithNonActiveDets = false;
					if (isnan(tangentCrd(detArray.at(item)._isoCylinder, line_, false).z))				//if no intersections with non-hit detector were found, the line is OK
					{
						noIntersectWithNonActiveDets = true;
					}
					else break;																			//if crosses any non-hit detector, breaks the loop, line is rejected
				}
				if (noIntersectWithNonActiveDets) {
					//intersectionPoints[planePosition].push_back(intersectPt);							//if line is OK, add to the intersection points

					intersectionPtCandidates[planePosition].push_back(intersectPt);
					tangencyPt[planePosition].push_back(tangPt);

				}
			}

		}
	}


	double workingDetX = workingDets.at(posPlane::XZ).at(0)._isoCylinder.x;
	double workingDetY = workingDets.at(posPlane::YZ).at(0)._isoCylinder.y;
	double intersectPtCenter = 0.0;
	for (const pair<posPlane, vector<pt3d>>& det : tangencyPt) {
		if (intersectionPtCandidates[planePosition].size() == 1) {
			intersectionPoints[planePosition] = intersectionPtCandidates[planePosition];
		}
		else {
			for (size_t i = 0; i < det.second.size(); i++) {
				for (size_t j = i; j < det.second.size(); j++) {

					if (i == j) continue;

					if (det.first == posPlane::XZ) {
						bool lhsTangents = ((det.second.at(i).x < workingDetX) &&
							(det.second.at(j).x < workingDetX));
						bool rhsTangents = ((det.second.at(i).x > workingDetX) &&
							(det.second.at(j).x > workingDetX));
						if (lhsTangents || rhsTangents) {
							intersectPtCenter = (intersectionPtCandidates.at(det.first).at(i).x + intersectionPtCandidates.at(det.first).at(j).x) / 2;
							intersectionPoints[posPlane::XZ].push_back({ intersectPtCenter, 0.0, z_planePosition });
						}
						else continue;
					}

					else if (det.first == posPlane::YZ) {
						bool lhsTangents = ((det.second.at(i).y < workingDetY) &&
							(det.second.at(j).y < workingDetY));
						bool rhsTangents = ((det.second.at(i).y > workingDetY) &&
							(det.second.at(j).y > workingDetY));
						if (lhsTangents || rhsTangents) {
							intersectPtCenter = (intersectionPtCandidates.at(det.first).at(i).y + intersectionPtCandidates.at(det.first).at(j).y) / 2;
							intersectionPoints[posPlane::YZ].push_back({ 0.0, intersectPtCenter, z_planePosition });
						}
						else continue;
					}

				}
			}

		}
	}

}

void GetIntersectionPoint_case2(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	const double& z_planePosition) {
	size_t i_start = 0;
	if (planePosition == posPlane::XZ) i_start = 0;
	else if (planePosition == posPlane::YZ) i_start = 4;
	size_t i_end = i_start + 4;

	vector<line> tangLines1 = tangents(workingDets.at(planePosition).at(0)._isoCylinder,
		workingDets.at(planePosition).at(1)._isoCylinder);
	for (line& line_ : tangLines1) {
																						//f(crd): a*crd + b*z +c == 0
		double intersect_crd = -(line_.b * z_planePosition + line_.c) / line_.a;		//required coordinate of the intersection point between tangent line and specific line at z_planePosition (@XZ plane)
		pt3d intersectPt{ 0.0, 0.0, z_planePosition };
		if (planePosition == posPlane::XZ) intersectPt.x = intersect_crd;
		else if (planePosition == posPlane::YZ) intersectPt.y = intersect_crd;			//intersection point at XZ plane, y coordinate is set to 0.0
		pt3d tangPt = tangentCrd(workingDets.at(planePosition).at(0)._isoCylinder, line_, true);//tangence point coordinate
		pt3d tangPtPerp = tangPt;													//projection of tangence point to z_planePosition 
		tangPtPerp.z = z_planePosition;												//so we can build the normal from tangPT to z_planePosition and build a triangle (tangPt, tangPtPerp, intersectPt)
		double incidenceAngle = getIncidenceAngle(intersectPt, tangPt, tangPtPerp);	//trajectory projection incidence angle is found as ASIN[(tangence line)/(perpendicular)] 

		bool noIntersectWithNonActiveDets = false;									//here we check if the found trajectory does not cross the cylinders of non-active detectors
		for (size_t item = i_start; item < i_end; item++) {
			if (detArray.at(item).isHit()) continue;								//we skip intesection with active detector, as they are obviously hit
			noIntersectWithNonActiveDets = false;
			if (isnan(tangentCrd(detArray.at(item)._isoCylinder, line_, false).z))	//if common point coordinates for the line and non-active detector are {NAN, NAN}
			{
				noIntersectWithNonActiveDets = true;								//there is no intersection between our line and detector cylinder, we can check next detector
			}
			else break;																//else it crosses the cylinder, and we have to break the cycle
		}
		if (noIntersectWithNonActiveDets) {
			intersectionPoints[planePosition].push_back(intersectPt);				//all legit points are stored in vector corresponding to the detector working plane
		}
	}
}

vector<pt3d> GetHitPointsByLines(std::vector<Detector>& detArray) {
	vector<pt3d> hitPoints = {};
	map<posPlane, vector<Detector>> workingDets = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };
	map<posPlane, vector<pt3d>> intersectionPoints;
	double z_planePosition = detArray.at(0).DetectorRadius * (1 + 2 * cos(30 * PI / 180)); //R + 2Rcos(30*) is the position of plane dividing XZ and YZ detectors
	//define working detectors
	for (const Detector& det : detArray) {
		if (det.isHit()) workingDets[det._isoCylinder.orientation].push_back(det);
	}
	//check if all 4 working detectors for each layer are chosen propely
	if (workingDets[posPlane::XZ].size() > 0) {
		switch (workingDets.at(posPlane::XZ).size()) {
			case 1:
			{
				GetIntersectionPoint_case1(detArray, posPlane::XZ, workingDets,
					intersectionPoints, z_planePosition);
				break;
			}
			case 2:
			{
				GetIntersectionPoint_case2(detArray, posPlane::XZ, workingDets,
					intersectionPoints, z_planePosition);
				break;
			}
		}
	}
	else if (workingDets[posPlane::XZ].size() == 0) {
		intersectionPoints[posPlane::XZ].push_back({ NAN, NAN, z_planePosition });
	}

	if (workingDets[posPlane::YZ].size() > 0) {
		switch (workingDets.at(posPlane::YZ).size()) {
			case 1:
			{
				GetIntersectionPoint_case1(detArray, posPlane::YZ, workingDets,
					intersectionPoints, z_planePosition);
				break;
			}
			case 2:
			{
				GetIntersectionPoint_case2(detArray, posPlane::YZ, workingDets,
					intersectionPoints, z_planePosition);
				break;
			}
		}
	}
	else if (workingDets[posPlane::YZ].size() == 0) {
		intersectionPoints[posPlane::YZ].push_back({ NAN, NAN, z_planePosition });
	}
	//get coordinates of hit points at the plane between XZ and YZ detectors
	for (auto& itemXZ : intersectionPoints.at(posPlane::XZ)) {
		for (auto& itemYZ : intersectionPoints.at(posPlane::YZ)) {
			hitPoints.push_back({ itemXZ.x, itemYZ.y, itemXZ.z });
		}
	}

	return hitPoints;
}


//----------------------------------TODO:
//this part is responsible for the calcualtion of line tangent to circle
//at the specified threshold angle

//if (incidenceAngle < thresholdAngle) {			//reject the line with unrealistic incidence angle
//	double central_z = workingDets.at(planePosition).at(0)._isoCylinder.z;
//	double R = workingDets.at(planePosition).at(0)._isoCylinder.r;
//	double perp_a = R * sin((90 - thresholdAngle) / 180.0 * PI);
//	double parall_b = R * cos((90 - thresholdAngle) / 180.0 * PI);
//	double perp_d = perp_a + (z_planePosition - central_z);
//	double parall_f = perp_d / tan((thresholdAngle) / 180.0 * PI);
//	double parall_g = parall_b;
// 
//	//add case for +thrAngle and -thrAngle!
//	//lower: (z_planePosition - central_z + perp_a)

//	if (planePosition == posPlane::XZ) {
//		double central_crd = workingDets.at(planePosition).at(0)._isoCylinder.x;
//		if (tangPt.x < central_crd) {
//			tangPtPerp = { (central_crd - parall_g), 0.0, z_planePosition };
//			tangPt = { (central_crd - parall_g), 0.0, (z_planePosition - perp_d) };
//			intersectPt = { (central_crd - (parall_f + parall_g)), 0.0, z_planePosition };
//		}
//		else {
//			tangPtPerp = { (central_crd + parall_g), 0.0, z_planePosition };
//			tangPt = { (central_crd + parall_g), 0.0, (z_planePosition - perp_d) };
//			intersectPt = { (central_crd + (parall_f + parall_g)), 0.0, z_planePosition };
//		}
//		line_ = {
//			tangPt.z - intersectPt.z,
//			intersectPt.x - tangPt.x,
//			-(intersectPt.x * (tangPt.z - intersectPt.z) +
//			intersectPt.z * (intersectPt.x - tangPt.x))
//		};
//	}
//	else if (planePosition == posPlane::YZ) {
//		double central_crd = workingDets.at(planePosition).at(0)._isoCylinder.y;
//		if (tangPt.y < central_crd) {
//			tangPtPerp = { 0.0, (central_crd - parall_g), z_planePosition };
//			tangPt = { 0.0, (central_crd - parall_g), (z_planePosition + (perp_a + R)) };
//			intersectPt = { 0.0, (central_crd - (parall_f + parall_g)), z_planePosition };
//		}
//		else {
//			tangPtPerp = { 0.0, (central_crd + parall_g), z_planePosition };
//			tangPt = { 0.0, (central_crd + parall_g), (z_planePosition + (perp_a + R)) };
//			intersectPt = { 0.0, (central_crd + (parall_f + parall_g)), z_planePosition };
//		}
//		line_ = {
//			tangPt.z - intersectPt.z,
//			intersectPt.y - tangPt.y,
//			-(intersectPt.y * (tangPt.z - intersectPt.z) +
//			intersectPt.z * (intersectPt.y - tangPt.y))
//		};
//	}
//}
//