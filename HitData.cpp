#include "HitData.h"
#include "Logger.h"

pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
	const double& _d, const double& _e, const double& _f, FileLogger& logger) {
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

	if (logger.getVerbosity() >= LogVerbosity::FULL) {
		print_state(iter, s, logger._info);
	}

	do
	{
		iter++;										//iterate while required tolerance is not reached
		status = gsl_multiroot_fsolver_iterate(s);

		if (logger.getVerbosity() >= LogVerbosity::FULL) {
			print_state(iter, s, logger._info);
		}

		if (status)									// check if solver is stuck
			break;
		status =
			gsl_multiroot_test_residual(s->f, 1e-6);
	} while (status == GSL_CONTINUE && iter < 1000);

	if (logger.getVerbosity() >= LogVerbosity::FULL) {
		logger._info << "\nStatus = " << gsl_strerror(status) << "\n";
		logger._info << "crd1: " << gsl_vector_get(s->x, 0) << " crd2: "
			<< gsl_vector_get(s->x, 1) << std::endl;
	}
	if (status == GSL_SUCCESS || (gsl_vector_get(s->f, 0) < 5e-4 &&
		(gsl_vector_get(s->f, 1) < 5e-4))) {

		if (logger.getVerbosity() >= LogVerbosity::FULL) {
			logger._info << gsl_vector_get(s->x, 0) << " " << gsl_vector_get(s->x, 1) << std::endl;
		}
	}
	else {
		if (logger.getVerbosity() >= LogVerbosity::FULL) {
			logger._info << "Roots not found\n";
		}
		return { NAN, NAN };
	}

	return { gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1) };
}


pt3d tangentCrd(const cylinder& cyl,
	const line& lin, FileLogger& logger) {
	double coef = 1.0e6; //round to 6th digit, else goes to nan
	double x1 = std::round(cyl.crd1 * coef) / coef;
	double y1 = std::round(cyl.crd2 * coef) / coef;
	double r1 = std::round(cyl.r * coef) / coef;
	double d = std::round(lin.a * coef) / coef;
	double e = std::round(lin.b * coef) / coef; //fixed swapped values
	double f = std::round(lin.c * coef) / coef;

	pt tangPt2d = tangentCrd(x1, y1, r1, d, e, f, logger);

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
	for (uint8_t i = 0b0000'0001, j = 0; (i <= 0b1000'0000 && j <= 7); i <<= 0b0000'0001, j++) {
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

double getAngleInRadians(const double& angle) {
	return angle * PI / 180;
}

double getAngleInDegrees(const double& angle) {
	return angle * 180 / PI;
}

double getIncidenceAngle(const pt3d& intersectPt, const pt3d& tangPt, const pt3d& tangPtPerp, FileLogger& logger) {
	double tangentLength = geomVectorLength(tangPt, intersectPt);
	double tangPerpLength = geomVectorLength(tangPt, tangPtPerp);
	double incidenceAngle = (asin(tangPerpLength / tangentLength) / PI * 180);
	std::cout << "\nAngle is: " << incidenceAngle << " degrees\n";
	if (logger.getVerbosity() >= LogVerbosity::FULL) logger._info << "\nAngle is: " << incidenceAngle << " degrees\n";
	return incidenceAngle;
}

//this recalculates the tangent line and corresponding points for the lines that incident at unrealistic angles
//a new tangent line is built according to the chosen threshold angle value
void RecalcTangLineWithAngle(const double& incidenceAngle, const double& thresholdAngle, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, const double& R, const double& z_planePosition, pt3d& tangPt,
	pt3d& tangPtPerp, pt3d& intersectPt, line& line_) {

	double zPt = std::fabs((tangPt.z - z_planePosition));
	double z0 =
		std::fabs(workingDets.at(planePosition).at(0)._isoCylinder.z - z_planePosition);


	if (planePosition == posPlane::XZ) {
		double x0 = workingDets.at(planePosition).at(0)._isoCylinder.x;

		double perpAB = R * std::cos(getAngleInRadians(thresholdAngle));
		double parOB = R * std::sin(getAngleInRadians(thresholdAngle));

		double XperpLen = 0;
		double ZparLen = 0;
		double XinterLen = 0;

		if (tangPt.x < x0) {						//if x < x0
			if (zPt > z0) {							//if z > z0
				ZparLen = z0 + perpAB;
				XperpLen = x0 - parOB;
				XinterLen = XperpLen - ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
			else {
				ZparLen = z0 - perpAB;
				XperpLen = x0 - parOB;
				XinterLen = XperpLen + ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
		}
		else {
			if (zPt > z0) {
				ZparLen = z0 + perpAB;
				XperpLen = x0 + parOB;
				XinterLen = XperpLen + ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
			else {
				ZparLen = z0 - perpAB;
				XperpLen = x0 + parOB;
				XinterLen = XperpLen - ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
		}
		tangPtPerp = { XperpLen, 0.0, z_planePosition };
		tangPt = { XperpLen, 0.0, (z_planePosition - ZparLen) };
		intersectPt = { XinterLen, 0.0, z_planePosition };

		line_ = {
			tangPt.z - intersectPt.z,
			intersectPt.x - tangPt.x,
			-(intersectPt.x * (tangPt.z - intersectPt.z) +
			intersectPt.z * (intersectPt.x - tangPt.x))
		};
	}

	else if (planePosition == posPlane::YZ) {
		double y0 = workingDets.at(planePosition).at(0)._isoCylinder.y;

		double perpAB = R * std::cos(getAngleInRadians(thresholdAngle));
		double parOB = R * std::sin(getAngleInRadians(thresholdAngle));

		double YperpLen = 0;
		double ZparLen = 0;
		double YinterLen = 0;

		if (tangPt.y < y0) {						//if y < y0
			if (zPt > z0) {							//if z > z0
				ZparLen = z0 + perpAB;
				YperpLen = y0 - parOB;
				YinterLen = YperpLen - ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
			else {
				ZparLen = z0 - perpAB;
				YperpLen = y0 - parOB;
				YinterLen = YperpLen + ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
		}
		else {
			if (zPt > z0) {
				ZparLen = z0 + perpAB;
				YperpLen = y0 + parOB;
				YinterLen = YperpLen + ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
			else {
				ZparLen = z0 - perpAB;
				YperpLen = y0 + parOB;
				YinterLen = YperpLen - ZparLen / std::tan(getAngleInRadians(thresholdAngle));
			}
		}
		tangPtPerp = { 0.0, YperpLen, z_planePosition };
		tangPt = { 0.0, YperpLen, (z_planePosition + ZparLen) };
		intersectPt = { 0.0, YinterLen, z_planePosition };

		line_ = {
			tangPt.z - intersectPt.z,
			intersectPt.y - tangPt.y,
			-(intersectPt.y * (tangPt.z - intersectPt.z) +
			intersectPt.z * (intersectPt.y - tangPt.y))
		};
	}

}

void GetIntersectionPoint_case1(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, map<posPlane, vector<pt3d>>& tangencyPoints,
	const double& z_planePosition, const double& thresholdAngle, FileLogger& logger) {
	size_t i_start = 0;
	if (planePosition == posPlane::XZ) i_start = 0;
	else if (planePosition == posPlane::YZ) i_start = 4;
	size_t i_end = i_start + 4;

	map<size_t, vector<line>> tangentLines;
	map<posPlane, vector<pt3d>> tangencyPt;
	map<posPlane, vector<pt3d>> intersectionPtCandidates;
	for (size_t i = i_start; i < i_end; i++) {
		if (detArray.at(i).isHit() &&
			detArray.at(i)._isoCylinder.x == workingDets.at(planePosition).at(0)._isoCylinder.x) {
			tangentLines[i] = {};						//skips the calculation for the only detector hit, we use it as anchor to find tangents to non-active detectors
			continue;
		}
		else {
			tangentLines[i] = tangents(workingDets.at(planePosition).at(0)._isoCylinder, detArray.at(i)._isoCylinder);	//iterate over detectors array to build tangents between working one and non-active dets
			for (line& line_ : tangentLines.at(i)) {

				//Logging tangent line equations
				logger._info.str("");
				logger._info << "\n";
				if (logger.getVerbosity() >= LogVerbosity::GEOM) PrintTangentLine(logger._info, line_, workingDets.at(planePosition).at(0)._isoCylinder);

				double intersect_crd = -(line_.b * z_planePosition + line_.c) / line_.a;
				pt3d intersectPt{ 0.0, 0.0, z_planePosition };
				if (planePosition == posPlane::XZ) intersectPt.x = intersect_crd;
				else if (planePosition == posPlane::YZ) intersectPt.y = intersect_crd;			//X or Y intersection coordinates

				pt3d tangPt = tangentCrd(workingDets.at(planePosition).at(0)._isoCylinder, line_, logger); //tangence point
				pt3d tangPtPerp = tangPt;
				tangPtPerp.z = z_planePosition;															//point coming from tangence point perpendicular to Z plane
				double incidenceAngle = getIncidenceAngle(intersectPt, tangPt, tangPtPerp, logger);				//calculate angle of incidence wrt Z axis
				double R = workingDets.at(planePosition).at(0)._isoCylinder.r;

				if (incidenceAngle < thresholdAngle)
				{
					logger._info << "\nIncidence angle is " << incidenceAngle << " deg., but threshold angle is "
						<< thresholdAngle << " deg. Recalculating tangent line.\n";
					RecalcTangLineWithAngle(incidenceAngle, thresholdAngle, planePosition, workingDets, R,
						z_planePosition, tangPt, tangPtPerp, intersectPt, line_);
					if (logger.getVerbosity() >= LogVerbosity::GEOM) PrintTangentLine(logger._info, line_, workingDets.at(planePosition).at(0)._isoCylinder);
					//prints recalculated line equation
				}

				bool noIntersectWithNonActiveDets = false;												//flag: if the tangent intersects any non-hit detector, which it is not tangent to, reject the line
				for (size_t item = i_start; item < i_end; item++) {
					if (detArray.at(item).isHit() ||													//for detectors where tangence line was recalculated, do not skip check of intersection with this detector
						((incidenceAngle >= thresholdAngle) && (item == i))) continue;					//skip if crosses the pair of detectors to which it is tangent
					noIntersectWithNonActiveDets = false;
					if (isnan(tangentCrd(detArray.at(item)._isoCylinder, line_, logger).z))				//if no intersections with non-hit detector were found, the line is OK
					{
						noIntersectWithNonActiveDets = true;
					}
					else break;																			//if crosses any non-hit detector, breaks the loop, line is rejected
				}
				if (noIntersectWithNonActiveDets) {
					logger._info << " ...OK!";
					//intersectionPoints[planePosition].push_back(intersectPt);							//if line is OK, add to the intersection points
					intersectionPtCandidates[planePosition].push_back(intersectPt);
					tangencyPt[planePosition].push_back(tangPt);

				}
				else logger._info << " ...NOT OK!";

				logger._info << "\n";
				logger.out("", logger._info.str(), LogVerbosity::GEOM);
			}

		}
	}
	//central point position calculation
	double workingDetX = 0.0;
	double workingDetY = 0.0;
	double intersectPtCenter = 0.0;
	double intersectPtUnc = 0.0;

	if (workingDets.at(posPlane::XZ).size() > 0) {
		workingDetX = workingDets.at(posPlane::XZ).at(0)._isoCylinder.x;
	}
	if (workingDets.at(posPlane::YZ).size() > 0) {
		workingDetY = workingDets.at(posPlane::YZ).at(0)._isoCylinder.y;
	}

	for (const pair<posPlane, vector<pt3d>>& det : tangencyPt) {
		if (intersectionPtCandidates[planePosition].size() == 1) {
			intersectionPoints[planePosition] = intersectionPtCandidates[planePosition];
		}
		else {
			for (size_t i = 0; i < det.second.size(); i++) {
				for (size_t j = i; j < det.second.size(); j++) {
					double R = workingDets.at(planePosition).at(0)._isoCylinder.r;
					if (i == j) continue;
					

					if (det.first == posPlane::XZ) {
						if (std::fabs(det.second.at(i).x - workingDetX) > R ||			//for case of 2 working detectors in the same layer, skip if points 
							std::fabs(det.second.at(j).x - workingDetX) > R) continue;	//do not belong to the current working detector

						bool lhsTangents = ((det.second.at(i).x < workingDetX) &&
							(det.second.at(j).x < workingDetX));
						bool rhsTangents = ((det.second.at(i).x > workingDetX) &&
							(det.second.at(j).x > workingDetX));
						if (lhsTangents || rhsTangents) {
							intersectPtCenter = (intersectionPtCandidates.at(det.first).at(i).x + intersectionPtCandidates.at(det.first).at(j).x) / 2;
							intersectPtUnc = std::fabs(intersectPtCenter - intersectionPtCandidates.at(det.first).at(i).x);
							intersectionPoints[posPlane::XZ].push_back({ intersectPtCenter, 0.0, z_planePosition });
							intersectionPointsUnc[posPlane::XZ].push_back({ intersectPtUnc, 0.0, (R / sqrt(12)) });


							pt3d tangPt = { 
								(det.second.at(i).x + det.second.at(j).x) / 2,
								0.0,
								(det.second.at(i).z + det.second.at(j).z) / 2 };
							logger._info.str("");
							tangencyPoints[planePosition].push_back(tangPt);
							logger._info << "\nTangency point: " << tangPt;
							logger._info << " at plane XZ\n";
							logger.out("", logger._info.str(), LogVerbosity::GEOM);
						}
						else continue;
					}

					else if (det.first == posPlane::YZ) {
						if (std::fabs(det.second.at(i).y - workingDetY) > R ||			//for case of 2 working detectors in the same layer, skip if points 
							std::fabs(det.second.at(j).y - workingDetY) > R) continue;	//do not belong to the current working detector

						bool lhsTangents = ((det.second.at(i).y < workingDetY) &&
							(det.second.at(j).y < workingDetY));
						bool rhsTangents = ((det.second.at(i).y > workingDetY) &&
							(det.second.at(j).y > workingDetY));
						if (lhsTangents || rhsTangents) {
							intersectPtCenter = (intersectionPtCandidates.at(det.first).at(i).y + intersectionPtCandidates.at(det.first).at(j).y) / 2;
							intersectPtUnc = std::fabs(intersectPtCenter - intersectionPtCandidates.at(det.first).at(i).y);
							intersectionPoints[posPlane::YZ].push_back({ 0.0, intersectPtCenter, z_planePosition });
							intersectionPointsUnc[posPlane::YZ].push_back({ 0.0, intersectPtUnc, (R / sqrt(12)) });

							pt3d tangPt = { 0.0, 
								(det.second.at(i).y + det.second.at(j).y) / 2, 
								(det.second.at(i).z + det.second.at(j).z) / 2 };
							logger._info.str("");
							tangencyPoints[planePosition].push_back(tangPt);
							logger._info << "\nTangency point: " << tangPt;
							logger._info << " at plane YZ\n";
							logger.out("", logger._info.str(), LogVerbosity::GEOM);
						}
						else continue;
					}

				}
			}

		}
	}
	//end of central point position calculation
}

void GetIntersectionPoint_case2(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, map<posPlane, vector<pt3d>>& tangencyPoints, 
	const double& z_planePosition, FileLogger& logger) {
	size_t i_start = 0;
	if (planePosition == posPlane::XZ) i_start = 0;
	else if (planePosition == posPlane::YZ) i_start = 4;
	size_t i_end = i_start + 4;
	double R = workingDets.at(planePosition).at(0)._isoCylinder.r;

	vector<line> tangLines1 = tangents(workingDets.at(planePosition).at(0)._isoCylinder,
		workingDets.at(planePosition).at(1)._isoCylinder);
	for (line& line_ : tangLines1) {

		//Logging tangent line equations
		logger._info.str("");
		logger._info << "\n";
		if (logger.getVerbosity() >= LogVerbosity::GEOM) PrintTangentLine(logger._info, line_, workingDets.at(planePosition).at(0)._isoCylinder);

		//f(crd): a*crd + b*z +c == 0

		double intersect_crd = -(line_.b * z_planePosition + line_.c) / line_.a;		//required coordinate of the intersection point between tangent line and specific line at z_planePosition (@XZ plane)
		pt3d intersectPt{ 0.0, 0.0, z_planePosition };
		if (planePosition == posPlane::XZ) intersectPt.x = intersect_crd;
		else if (planePosition == posPlane::YZ) intersectPt.y = intersect_crd;			//intersection point at XZ plane, y coordinate is set to 0.0
		pt3d tangPt = tangentCrd(workingDets.at(planePosition).at(0)._isoCylinder, line_, logger);//tangence point coordinate
		pt3d tangPtPerp = tangPt;													//projection of tangence point to z_planePosition 
		tangPtPerp.z = z_planePosition;												//so we can build the normal from tangPT to z_planePosition and build a triangle (tangPt, tangPtPerp, intersectPt)
		double incidenceAngle = getIncidenceAngle(intersectPt, tangPt, tangPtPerp, logger);	//trajectory projection incidence angle is found as ASIN[(tangence line)/(perpendicular)] 

		bool noIntersectWithNonActiveDets = false;									//here we check if the found trajectory does not cross the cylinders of non-active detectors
		for (size_t item = i_start; item < i_end; item++) {
			if (detArray.at(item).isHit()) continue;								//we skip intesection with active detector, as they are obviously hit
			noIntersectWithNonActiveDets = false;
			if (isnan(tangentCrd(detArray.at(item)._isoCylinder, line_, logger).z))	//if common point coordinates for the line and non-active detector are {NAN, NAN}
			{
				noIntersectWithNonActiveDets = true;								//there is no intersection between our line and detector cylinder, we can check next detector
			}
			else break;																//else it crosses the cylinder, and we have to break the cycle
		}

		if (noIntersectWithNonActiveDets) {
			logger._info << " ...OK!";
			intersectionPoints[planePosition].push_back(intersectPt);				//all legit points are stored in vector corresponding to the detector working plane
			
			tangencyPoints[planePosition].push_back(tangPt);
			logger._info << "\nTangency point: " << tangPt;

			if (planePosition == posPlane::XZ) {
				logger._info << " at plane XZ\n";
				intersectionPointsUnc[planePosition].push_back({ (R / sqrt(12)), 0.0, (R / sqrt(12)) });
			}
			else if (planePosition == posPlane::YZ) {
				logger._info << " at plane YZ\n";
				intersectionPointsUnc[planePosition].push_back({ 0.0, (R / sqrt(12)), (R / sqrt(12)) });
			}
		}
		else logger._info << " ...NOT OK!";

		logger._info << "\n";
		logger.out("", logger._info.str(), LogVerbosity::GEOM);
	}
}

void GetIntersectionPoint_Multi(const std::vector<Detector>& detArray, const posPlane& planePosition,
	const map<posPlane, vector<Detector>>& workingDets, map<posPlane, vector<pt3d>>& intersectionPoints,
	map<posPlane, vector<pt3d>>& intersectionPointsUnc, map<posPlane, vector<pt3d>>& tangencyPoints,
	const double& z_planePosition, FileLogger& logger) {
	size_t i_start = 0;
	if (planePosition == posPlane::XZ) i_start = 0;
	else if (planePosition == posPlane::YZ) i_start = 4;
	size_t i_end = i_start + 4;
	double R = workingDets.at(planePosition).at(0)._isoCylinder.r;


	vector<line> tangLines1 = {};
	for (size_t i = 0; i < workingDets.at(planePosition).size(); i++) {
		for (size_t item = i; item < workingDets.at(planePosition).size(); item++) {
			if (i == item) continue;
			if (workingDets.at(planePosition).at(i)._isoCylinder.z ==
				workingDets.at(planePosition).at(item)._isoCylinder.z) continue;

			else {
				tangLines1 = tangents(workingDets.at(planePosition).at(i)._isoCylinder,
					workingDets.at(planePosition).at(item)._isoCylinder);


				for (line& line_ : tangLines1) {
					//Logging tangent line equations
					logger._info.str("");
					logger._info << "\n";
					if (logger.getVerbosity() >= LogVerbosity::GEOM) PrintTangentLine(logger._info, line_, workingDets.at(planePosition).at(0)._isoCylinder);

					//f(crd): a*crd + b*z +c == 0
					double intersect_crd = -(line_.b * z_planePosition + line_.c) / line_.a;		//required coordinate of the intersection point between tangent line and specific line at z_planePosition (@XZ plane)
					pt3d intersectPt{ 0.0, 0.0, z_planePosition };
					if (planePosition == posPlane::XZ) intersectPt.x = intersect_crd;
					else if (planePosition == posPlane::YZ) intersectPt.y = intersect_crd;			//intersection point at XZ plane, y coordinate is set to 0.0
					pt3d tangPt = tangentCrd(workingDets.at(planePosition).at(0)._isoCylinder, line_, logger);//tangence point coordinate
					pt3d tangPtPerp = tangPt;													//projection of tangence point to z_planePosition 
					tangPtPerp.z = z_planePosition;												//so we can build the normal from tangPT to z_planePosition and build a triangle (tangPt, tangPtPerp, intersectPt)
					double incidenceAngle = getIncidenceAngle(intersectPt, tangPt, tangPtPerp, logger);	//trajectory projection incidence angle is found as ASIN[(tangence line)/(perpendicular)] 

					bool noIntersectWithNonActiveDets = false;									//here we check if the found trajectory does not cross the cylinders of non-active detectors
					for (size_t item = i_start; item < i_end; item++) {
						if (detArray.at(item).isHit()) continue;								//we skip intersection with active detectors, as they are obviously hit
						noIntersectWithNonActiveDets = false;
						if (isnan(tangentCrd(detArray.at(item)._isoCylinder, line_, logger).z))	//if common point coordinates for the line and non-active detector are {NAN, NAN}
						{
							noIntersectWithNonActiveDets = true;								//there is no intersection between our line and detector cylinder, we can check next detector
						}
						else break;																//else it crosses the cylinder, and we have to break the cycle
					}

					if (workingDets.at(planePosition).size() == 4) noIntersectWithNonActiveDets = true; //all detectors in the layers are active

					if (noIntersectWithNonActiveDets) {
						logger._info << " ...OK!";
						intersectionPoints[planePosition].push_back(intersectPt);				//all legit points are stored in vector corresponding to the detector working plane

						tangencyPoints[planePosition].push_back(tangPt);
						logger._info << "\nTangency point: " << tangPt;

						if (planePosition == posPlane::XZ) {
							logger._info << " at plane XZ\n";
							intersectionPointsUnc[planePosition].push_back({ (R / sqrt(12)), 0.0, (R / sqrt(12)) });
						}
						else if (planePosition == posPlane::YZ) {
							logger._info << " at plane YZ\n";
							intersectionPointsUnc[planePosition].push_back({ 0.0, (R / sqrt(12)), (R / sqrt(12)) });
						}

					}
					else logger._info << " ...NOT OK!";

					logger._info << "\n";
					logger.out("", logger._info.str(), LogVerbosity::GEOM);
				}

			}

		}
	}

}

void GetHitPointsByLines(std::vector<Detector>& detArray, const double& thresholdAngle, vector<pt3d>& hitPoints,
	vector<pt3d>& hitPointsUncert, vector<pt3d>& directionVectors, FileLogger& logger) {
	hitPoints = {};
	map<posPlane, vector<Detector>> workingDets = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };
	map<posPlane, vector<pt3d>> intersectionPoints = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };
	map<posPlane, vector<pt3d>> intersectionPointsUnc = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };
	map<posPlane, vector<pt3d>> tangencyPoints = { {posPlane::XZ, {}}, {posPlane::YZ, {}} };


	double z_planePosition = detArray.at(0).DetectorRadius * (1 + 2 * cos(30 * PI / 180)); //R + 2Rcos(30*) is the position of plane dividing XZ and YZ detectors
	//define working detectors
	size_t det_id = 1;
	for (const Detector& det : detArray) {

		logger._info.str("");
		logger._info << "\nDetector No. " << det_id << " with R = " << det.DetectorRadius;
		if (det.isHit()) {
			logger._info << " is HIT,";
			logger._info << " isochrone radius is " << det._isoCylinder.r << "\n";

			workingDets[det._isoCylinder.orientation].push_back(det);
		}
		else logger._info << " is NOT HIT\n";

		if (logger.getVerbosity() >= LogVerbosity::GEOM) logger._info << det;
		logger.out("", logger._info.str(), LogVerbosity::GEOM);
		det_id++;
	}

	logger._info.str("");
	logger._info << "\nControl XY plane is placed at z = " << z_planePosition << "\n";
	logger.out("", logger._info.str(), LogVerbosity::GEOM);

	//check if all 4 working detectors for each layer are chosen propely
	if (workingDets[posPlane::XZ].size() > 0) {
		switch (workingDets.at(posPlane::XZ).size()) {
		case 1:
		{
			GetIntersectionPoint_case1(detArray, posPlane::XZ, workingDets,
				intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);
			break;
		}
		case 2:
		{
			//check if 2 workind detectors are in the same layer
			if (workingDets.at(posPlane::XZ).at(0)._isoCylinder.z ==
				workingDets.at(posPlane::XZ).at(1)._isoCylinder.z) {
				map<posPlane, vector<Detector>> workingDetsSameLayer = workingDets;
				workingDetsSameLayer.at(posPlane::XZ).clear();
				workingDetsSameLayer.at(posPlane::XZ).push_back(workingDets.at(posPlane::XZ).at(0));
				GetIntersectionPoint_case1(detArray, posPlane::XZ, workingDetsSameLayer,
					intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);

				workingDetsSameLayer.at(posPlane::XZ).clear();
				workingDetsSameLayer.at(posPlane::XZ).push_back(workingDets.at(posPlane::XZ).at(1));
				GetIntersectionPoint_case1(detArray, posPlane::XZ, workingDetsSameLayer,
					intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);
			}
			else {
				GetIntersectionPoint_case2(detArray, posPlane::XZ, workingDets,
					intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, logger);
			}
			break;
		}
		default: {
			GetIntersectionPoint_Multi(detArray, posPlane::XZ, workingDets,
				intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, logger);
		}
		}
	}
	if (workingDets[posPlane::XZ].size() == 0 || intersectionPoints[posPlane::XZ].size() == 0) {
		intersectionPoints[posPlane::XZ].push_back({ NAN, NAN, z_planePosition });
	}

	if (workingDets[posPlane::YZ].size() > 0) {
		switch (workingDets.at(posPlane::YZ).size()) {
		case 1:
		{
			GetIntersectionPoint_case1(detArray, posPlane::YZ, workingDets,
				intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);
			break;
		}
		case 2:
		{
			//check if 2 workind detectors are in the same layer
			if (workingDets.at(posPlane::YZ).at(0)._isoCylinder.z ==
				workingDets.at(posPlane::YZ).at(1)._isoCylinder.z) {
				map<posPlane, vector<Detector>> workingDetsSameLayer = workingDets;
				workingDetsSameLayer.at(posPlane::YZ).clear();
				workingDetsSameLayer.at(posPlane::YZ).push_back(workingDets.at(posPlane::YZ).at(0));
				GetIntersectionPoint_case1(detArray, posPlane::YZ, workingDetsSameLayer,
					intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);

				workingDetsSameLayer.at(posPlane::YZ).clear();
				workingDetsSameLayer.at(posPlane::YZ).push_back(workingDets.at(posPlane::YZ).at(1));
				GetIntersectionPoint_case1(detArray, posPlane::YZ, workingDetsSameLayer,
					intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, thresholdAngle, logger);
			}
			else GetIntersectionPoint_case2(detArray, posPlane::YZ, workingDets,
				intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, logger);
			break;
		}
		default: {
			GetIntersectionPoint_Multi(detArray, posPlane::YZ, workingDets,
				intersectionPoints, intersectionPointsUnc, tangencyPoints, z_planePosition, logger);
		}
		}
	}
	if (workingDets[posPlane::YZ].size() == 0 || intersectionPoints[posPlane::YZ].size() == 0) {
		intersectionPoints[posPlane::YZ].push_back({ NAN, NAN, z_planePosition });
	}

	if (workingDets[posPlane::XZ].size() == 0 &&
		workingDets[posPlane::YZ].size() == 0) {
		hitPoints = {};
		hitPointsUncert = {};
		return;
	}

	//get coordinates of hit points at the plane between XZ and YZ detectors
	for (auto& itemXZ : intersectionPoints.at(posPlane::XZ)) {
		for (auto& itemYZ : intersectionPoints.at(posPlane::YZ)) {
			hitPoints.push_back({ itemXZ.x, itemYZ.y, itemXZ.z });
		}
	}
	if (intersectionPointsUnc.at(posPlane::XZ).size() > 0 &&						//TODO: change for case of defined uncertainties in 1 layer
		intersectionPointsUnc.at(posPlane::YZ).size() > 0) {
		for (auto& itemXZ : intersectionPointsUnc.at(posPlane::XZ)) {
			for (auto& itemYZ : intersectionPointsUnc.at(posPlane::YZ)) {
				hitPointsUncert.push_back({ itemXZ.x, itemYZ.y, itemXZ.z });
			}
		}
	}

	//this code executes the calculation of 
	if (hitPoints.size() > 0) {
		size_t i = 0;
		for (auto& itemXZ : tangencyPoints.at(posPlane::XZ)) {
			for (auto& itemYZ : tangencyPoints.at(posPlane::YZ)) {
				pt3d intersectPt = hitPoints.at(i);
				double x_1 = itemXZ.x;
				double x_2 = intersectPt.x;
				double y_1 = itemYZ.y;
				double y_2 = intersectPt.y;
				double zx_1 = itemXZ.z;
				double z_2 = intersectPt.z;
				double zy_1 = itemYZ.z;

				double x_on_y = ((x_2 - x_1) * zy_1 + (x_1 * (z_2 - zx_1) - zx_1 * (x_2 - x_1))) / (z_2 - zx_1);
				double y_on_x = ((y_2 - y_1) * zx_1 + (y_1 * (z_2 - zy_1) - zy_1 * (y_2 - y_1))) / (z_2 - zy_1);

				pt3d tangPtX = { x_1, y_on_x, zx_1 };
				pt3d tangPtY = { x_on_y, y_1, zy_1 };

				pt3d directVector = tangPtY - tangPtX;

				logger._info.str("");
				logger._info << "\nHit point No. " << i << " related tangent line reconstruction points: "
					<< "\nX layer: " << tangPtX
					<< "\nZ-plane: " << intersectPt
					<< "\nY layer: " << tangPtY << "\n";
				if (logger.getVerbosity() >= LogVerbosity::FULL) {
					logger._info << "Direction vector: { "
						<< directVector.x << ", " << directVector.y << ", " << directVector.z << " }\n";
				}
				logger.out("", logger._info.str(), LogVerbosity::GEOM);
				
				directionVectors.push_back(directVector);
				i++;
			}
		}
	}

	return;
}
