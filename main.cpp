#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {

	//GetDetectorHitData(0b0110'1010);
	circle a1;
	a1.r = 5.2;
	a1.crd1 = 0.3;
	a1.crd2 = 2.0;
	a1.orientation = posPlane::XZ;
	circle a2;
	a2.r = 2.4;
	a2.crd1 = 11.6;
	a2.crd2 = 4.1;
	a2.orientation = posPlane::YZ;

	vector<line> ll1 = tangents(a1, a2);
	for (auto& i : ll1) {
		cout << i.a << " " << i.c << " " << i.b << endl;
	}
	pt res = tangentCrd(a1.crd1, a1.crd2, a1.r,
		ll1[0].a,
		ll1[0].c,
		ll1[0].b);
	cylinder cyl1(a1, 0.0);
	for (auto& i : ll1) {
		pt res1 = tangentCrd(a1.crd1, a1.crd2, a1.r,
			i.a,
			i.c,
			i.b);
		pt3d tangPt = { res1.crd1, 0.0, res1.crd2 };
		cout << res1.crd1 << " " << res1.crd2 << endl;
		plane pp = TangPlaneToCylinder(tangPt, cyl1);
		cout << "Plane equation: " << pp.a;
		(pp.b >= 0) ? cout << "x + " : cout << "x "; cout << pp.b; 
		(pp.c >= 0) ? cout << "y + " : cout << "y "; cout << pp.c; 
		(pp.d >= 0) ? cout << "z + " : cout << "z "; cout << pp.d << endl;
	}

	GetDetectorHitData(0b01010101, 10.0);


	//std::cout << res.x << " " << res.y << endl;
	return 0;
}