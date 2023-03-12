#include <cstdint>
#include <iostream>

#include "HitData.h"


int main() {

	//GetDetectorHitData(0b0110'1010);
	circle a1;
	a1.r = 5.2;
	a1.x = 0.3;
	a1.y = 2.0;
	circle a2;
	a2.r = 2.4;
	a2.x = 11.6;
	a2.y = 4.1;
	vector<line> ll1 = tangents(a1, a2);
	for (auto i : ll1) {
		cout << i.a << " " << i.c << " " << i.b << endl;
	}
	pt res = tangentCrd(a1.x, a1.y, a1.r,
		ll1[1].a,
		ll1[1].c,
		ll1[1].b);
	for (auto i : ll1) {
		pt res1 = tangentCrd(a1.x, a1.y, a1.r,
			i.a,
			i.c,
			i.b);
		cout << res1.x << " " << res1.y << endl;
	}
	//std::cout << res.x << " " << res.y << endl;
	return 0;
}