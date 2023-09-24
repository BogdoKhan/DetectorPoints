#include <cstdint>
#include <iostream>
#include <sstream>

#include "HitData.h"
#include "Logger.h"

//TODO: move the solver output and add info about the created lines, circles, planes, points into the log files
//only required info (hit points +- uncertainty) must be pushed to the cout

//add visualization? for example, add function that fills in the gnuplot script, so you can execute it
//all geometrical primitives may be presented as the set of points

int main() {
	//std::vector<double> isochrones = {4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0};
	//std::vector<double> isochrones = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	LogVerbosity verbosity = LogVerbosity::GEOM;		//set verbosity level, should be user-defined
	std::ofstream logfile;								//create output log file
	stringstream info;									//buffer for complex objects before logging them
	FileLogger logger(verbosity, logfile, info);				//create logger that works with the log file
			
	std::vector<double> isochrones = { 4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0 };
	double detectorRadius = 5.0;
	std::vector<Detector> detArray;
	std::vector<pt3d> hitPointsByLines = {};
	std::vector<pt3d> hitPointsUncert = {};
	double thresholdAngle = 30.0;

	GetDetectorHitData(isochrones, detArray, detectorRadius);
	size_t counter = 1;
	GetHitPointsByLines(detArray, thresholdAngle, hitPointsByLines, hitPointsUncert, logger);
	////pt3d hitPointsErr = {2*10.0/sqrt(12), 2 * 10.0 / sqrt(12) , 2 * 10.0 / sqrt(12) };
	//counter = 1;
	//cout << "Hit points, new counter: \n";
	//for (pt3d& i : hitPointsByLines) {
	//	cout << "No. " << counter << " " << i << endl;
	//	counter++;
	//}

	BufferDetectorHitsPoints(hitPointsByLines,
		hitPointsUncert, counter, std::cout);

	logger._info.str("");
	BufferDetectorHitsPoints(hitPointsByLines,
		hitPointsUncert, counter, logger._info);
	logger.out("Hits data:\n ", logger._info.str(), LogVerbosity::BASIC);

	return 0;
}