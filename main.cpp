#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {
	std::vector<double> isochrones = {
		5.2, 0.0, 2.4, 0.0, 3.2, 0.0, 7.6, 0.0
	};
	std::vector<Detector> detArray;
	//GetDetectorHitData(0b0110'1010);
	//circle a1;
	//a1.r = 5.2;
	//a1.crd1 = 0.0;
	//a1.crd2 = 0.0;
	//a1.orientation = posPlane::XZ;
	//circle a2;
	//a2.r = 2.4;
	//a2.crd1 = 10.0;
	//a2.crd2 = 10.0;
	//a2.orientation = posPlane::YZ;

	//vector<line> ll1 = tangents(a1, a2);
	//for (auto& i : ll1) {
	//	cout << i.a << " " << i.c << " " << i.b << endl;
	//}
	//pt res = tangentCrd(a1.crd1, a1.crd2, a1.r,
	//	ll1[0].a,
	//	ll1[0].c,
	//	ll1[0].b);
	//cylinder cyl1(a1, 0.0);
	//for (auto& i : ll1) {
	//	pt res1 = tangentCrd(a1.crd1, a1.crd2, a1.r,
	//		i.a,
	//		i.c,
	//		i.b);
	//	pt3d tangPt = { res1.crd1, 0.0, res1.crd2 };
	//	cout << "Tangence point: {" << tangPt.x << ", " << tangPt.y << ", " << tangPt.z << "}" << endl;
	//	plane pp = TangPlaneToCylinder(tangPt, cyl1);
	//	cout << pp << endl;
	//}

	GetDetectorHitData(0b01010101, isochrones, detArray, 10.0);
	vector <pt3d> hitPoints = GetHitPoints(detArray);
	size_t counter = 1;
	for (pt3d& i : hitPoints) {
		cout << "Point: " << counter << " {" << i.x << ", " << i.y << ", " << i.z << "} \n";
		counter++;
	}

	//std::cout << res.x << " " << res.y << endl;
	return 0;
}