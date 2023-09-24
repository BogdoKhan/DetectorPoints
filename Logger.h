#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

#include "Tangents.h"

enum LogVerbosity {
	NONE = 0,
	BASIC = 1,
	GEOM = 2,
	FULL = 3
};

class Logger {
public:
	Logger();

	virtual LogVerbosity& getVerbosity();
	virtual void setVerbosity(const LogVerbosity& verb);

private:
	LogVerbosity _verbosity = LogVerbosity::NONE;
};


class FileLogger : public Logger {
public:
	FileLogger(const LogVerbosity& verb, std::ofstream& outFile, std::stringstream& buffer);
	~FileLogger();

	LogVerbosity& getVerbosity();

	template<typename T, typename U>
	std::ofstream& out(const U& caption, const T& data, const LogVerbosity& requiredVerb) {
		if (_outFile.is_open() && (_verbosity >= requiredVerb)) {
			_outFile << caption;
			_outFile << data;
		}
		else if (!_outFile.is_open()) std::cout << "File is not accessible or not opened\n";
		return _outFile;
	}

	std::ofstream& _outFile;
	std::stringstream& _info;
private:
	LogVerbosity _verbosity;
};

void BufferDetectorHitsPoints(
	const std::vector<pt3d>& hitPointsByLines,
	const std::vector<pt3d>& hitPointsUncert,
	size_t counter,
	ostream& ss);