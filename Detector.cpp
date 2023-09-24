#include "Detector.h"

Detector::Detector() {

};

Detector::Detector(std::string emplacement, bool isHit) :
	_emplacement(emplacement),
	_isHit(isHit)
{
	std::cout << "Detector at " << _emplacement;
	if (_isHit) {
		std::cout << " with a hit";
	}
	else {
		std::cout << " with no hit";
	}
	std::cout << " is created\n";
};

bool Detector::isHit() const { return _isHit; }

bool operator< (const Detector& lhs, const Detector& rhs) {
	if (lhs._isoCylinder.crd1 < rhs._isoCylinder.crd1) return true;
	if (lhs._isoCylinder.orientation < rhs._isoCylinder.orientation) return true;
	else return false;
}

std::ostream& operator<< (std::ostream& out, const Detector& det) {
	pt3d detCenter = { det._isoCylinder.x, det._isoCylinder.y, det._isoCylinder.z };
	out << "Detector center is placed at: " << detCenter << "\n";
	out << "Isochrone equation is: "; 
	if (det._isoCylinder.orientation == posPlane::XZ) {
		(det._isoCylinder.x >= 0) ? out << "(x-" << det._isoCylinder.x : out << "(x+" << det._isoCylinder.x; out << ")^2 + ";
		(det._isoCylinder.z >= 0) ? out << "(z-" << det._isoCylinder.z : out << "(z+" << det._isoCylinder.z; out << ")^2 - ";
		out << pow(det._isoCylinder.r, 2) << " == 0\n";
	}
	else {
		(det._isoCylinder.y >= 0) ? out << "(y-" << det._isoCylinder.y : out << "(y+" << det._isoCylinder.y; out << ")^2 + ";
		(det._isoCylinder.z >= 0) ? out << "(z-" << det._isoCylinder.z : out << "(z+" << det._isoCylinder.z; out << ")^2 - ";
		out << pow(det._isoCylinder.r, 2) << " == 0\n";
	}
	
	return out;
}