#include "Logger.h"


Logger::Logger():
	_verbosity(LogVerbosity::NONE)
{}

LogVerbosity& Logger::getVerbosity() {
	return _verbosity;
}

void Logger::setVerbosity(const LogVerbosity& verb) {
	_verbosity = verb;
}

FileLogger::FileLogger(const LogVerbosity& verb, std::ofstream& outFile, std::stringstream& buffer):
	_verbosity(verb), _outFile(outFile), _info(buffer)
{
	if (_verbosity > LogVerbosity::NONE) {
		_outFile.open("log_run.txt", std::ofstream::out | std::ofstream::trunc);		//open log file
		out("Starting the log file...\n", "", LogVerbosity::BASIC);
	}
}

FileLogger::~FileLogger() {
	if (_outFile.is_open()) {
		out("\nClosing the log file...", "", LogVerbosity::BASIC);
		_outFile.close();
	}
}

LogVerbosity& FileLogger::getVerbosity() {
	return _verbosity;
}


void BufferDetectorHitsPoints(
	const std::vector<pt3d>& hitPointsByLines,
	const std::vector<pt3d>& hitPointsUncert,
	const std::vector<pt3d>& directionVectors,
	size_t counter,
	ostream& ss){
	for (size_t i = 0; i < hitPointsByLines.size(); i++) {
		ss << "No. " << counter << " " << hitPointsByLines.at(i) << endl;
		if (hitPointsUncert.size() > 0) {
			ss << "Uncert +-" << hitPointsUncert.at(i) << endl;
		}
		if (directionVectors.size() > 0 && (directionVectors.size() == hitPointsByLines.size())) {
			ss << "Direction vector: { "
				<< directionVectors.at(i).x << ", " << directionVectors.at(i).y
				<< ", " << directionVectors.at(i).z << " }\n";
		}
		ss << "--------------------------------------------------------------\n";
		counter++;
	}
}