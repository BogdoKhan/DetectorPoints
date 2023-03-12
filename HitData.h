#pragma once
#include <iostream>
#include <bitset>
#include <map>
#include <vector>
#include <string>
#include <math.h>

#include "Detector.h"

//0b|y22|y21|y12|y11|'|x22|x21|x12|x11|

using namespace std;

struct pt {
	double x, y;

	pt operator- (pt p) {
		pt res = { x - p.x, y - p.y };
		return res;
	}
};

struct circle : pt {
	double r;
};

struct pt3d {
	double x,y,z;
	pt3d operator- (pt3d p) {
		pt3d res = { x - p.x, y - p.y, z - p.z };
		return res;
	}
};

struct line {
	double a, b, c;
};

const double EPS = 1E-9;

double sqr(double a) {
	return a * a;
}

void tangents(pt c, double r1, double r2, vector<line>& ans) {
	double r = r2 - r1;
	double z = sqr(c.x) + sqr(c.y);
	double d = z - sqr(r);
	if (d < -EPS)  return;
	d = sqrt(abs(d));
	line l;
	l.a = (c.x * r + c.y * d) / z;
	l.b = (c.y * r - c.x * d) / z;
	l.c = r1;
	ans.push_back(l);
}

vector<line> tangents(circle a, circle b) {
	vector<line> ans;
	for (int i = -1; i <= 1; i += 2)
		for (int j = -1; j <= 1; j += 2)
			tangents(b - a, a.r * i, b.r * j, ans);
	for (size_t i = 0; i < ans.size(); ++i)
		ans[i].c -= ans[i].a * a.x + ans[i].b * a.y;
	return ans;
}

//(x-x1)^2 + (y-y1)^2 == R1^2;
//dx + fy + e = 0
//returns coordinates of tangence point between line & circle1
//(x_tang, y_tang)
pt tangentCrd(double _x1, double _y1, double _r1,
	double _d, double _e, double _f) {
	int coef = 1e6; //round to 6th digit, else goes to nan
	double x1 = std::round(_x1 * coef) / coef;
	double y1 = std::round(_y1 * coef) / coef;
	double r1 = std::round(_r1 * coef) / coef;
	double d = std::round(_d * coef) / coef;
	double e = std::round(_e * coef) / coef;
	double f = std::round(_f * coef) / coef;

	double resultX = 0.0;
	double resultY = 0.0;

	//magic equations for x and y 
	//calcualted with Wolfram Mathematica
	//Solve[{(x - x1)^2 + (y - y1)^2 == r1^2, 
	//y == -((d*x + e)/f)}, {x, y}]
	//this is the same as
	//Solve[{(x - a)^2 + (y - b)^2 == r^2, 
	//d*x + f*y + e == 0}, {x, y}]
	//X coordinate
	resultX = (1 / (sqr(d) + sqr(f)) * 0.5 *
	(-2 * d * e - 2 * y1 * d * f + 2 * x1 * sqr(f) -
		sqrt(sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f)) -
			4 * (sqr(d) + sqr(f)) *
			(sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
				sqr(y1) * sqr(f) - sqr(f) * sqr(r1)
				)
			)
		)
	);
	//Y coordinate
	resultY = (1 / f) * 
		(-e + (sqr(d) * e) / (sqr(d) + sqr(f)) +
		(y1 * sqr(d) * f) / (sqr(d) + sqr(f)) -
		(x1 * d * sqr(f)) / (sqr(d) + sqr(f)) +
		1 / (2 * (sqr(d) + sqr(f))) * d *
		sqrt(sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f)) -
			4 * (sqr(d) + sqr(f)) *
			(sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
				sqr(y1) * sqr(f) - sqr(f) * sqr(r1))
			)
		);

	pt result{ resultX, resultY };
	return result;
}

void GetDetectorHitData(const uint8_t& word ) {
	std::vector<std::string> DetectorNames = {
		"X11", "X12", "X21", "X22", "Y11", "Y12", "Y21", "Y22"
	};
	std::vector<Detector> detArray;
	for (size_t i = 0b0000'0001, j = 0; i <= 0b1000'0000; i <<= 0b0000'0001, j++) {
		std::cout << (std::bitset<8>)word << " "
			<< (std::bitset<8>)i << " " << ((std::bitset<8>)(word & i)).any() << " " << DetectorNames.at(j) << std::endl;
		detArray.push_back(Detector(DetectorNames.at(j), ((std::bitset<8>)(word & i)).any(), 5.32));
	}

}