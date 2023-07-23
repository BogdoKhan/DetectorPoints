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