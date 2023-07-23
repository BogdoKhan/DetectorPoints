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

double getAngleInRadians(const double& angle);

double getAngleInDegrees(const double& angle);

double getIncidenceAngle(const pt3d& intersectPt, const pt3d& tangPt, const pt3d& tangPtPerp);

void RecalcTangLineWithAngle(const double& incidenceAngle, const double& thresholdAngle, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, const double& R, const double& z_planePosition, pt3d& tangPt,
	pt3d& tangPtPerp, pt3d& intersectPt, line& line_);

void GetIntersectionPoint_case1(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, const double& z_planePosition, const double& thresholdAngle);

void GetIntersectionPoint_case2(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, const double& z_planePosition);

void GetIntersectionPoint_Multi(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, const double& z_planePosition);

void GetHitPointsByLines(std::vector<Detector>& detArray, const double& thresholdAngle, vector<pt3d>& hitPoints,
	vector<pt3d>& hitPointsUncert);


