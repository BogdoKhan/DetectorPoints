#pragma once
#include <vector>
#include <iostream>

using namespace std;

enum class posPlane {
	XZ, YZ, XY, NA
};

//tangent to two circles and corresponding structures
//equation is taken from https://e-maxx.ru/algo/circle_tangents
struct pt {
	double crd1, crd2; //coordinate 1 & coordinate 2, not X, Y, etc
						//common names because in XZ and YZ cases coords are different
	pt operator- (pt p);
};
//a circle, has radius r and central point at crd1, crd2
struct circle : pt {
public:
	double r = 0.0;
	posPlane orientation = posPlane::NA;	//orientation of detector cross section: XZ, YZ, XY
	circle(double crd1, double crd2, double r);
	circle();
};

struct pt3d {
	double x, y, z;	//coordinates in space based on cross sectional orientation
	pt3d();
	pt3d(double x_, double y_, double z_);
	pt3d operator- (pt3d p);

};
//cylinder is based on circle
struct cylinder : circle {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
public:
	cylinder(circle c, double crd3);
	cylinder();
};
//common equation of line
struct line {
	double a, b, c;
};

const double EPS = 1E-9;//epsilon, defines how accurate we want to be

double sqr(double a);
//finds inner and outer tangents to 2 circles
void tangents(pt c, double r1, double r2, vector<line>& ans);
//creates vector of 4 tangents
vector<line> tangents(circle a, circle b);

vector<line> tangents(cylinder cyl1, cylinder cyl2);
//common equation of plane
struct plane {
	double a;
	double b;
	double c;
	double d;

	plane() : a(0), b(0), c(0), d(0) {};
	plane(double _a, double _b, double _c, double _d) :
		a(_a), b(_b), c(_c), d(_d) {};

};

std::ostream& operator<< (std::ostream& out, const plane& pp);

std::ostream& operator<< (std::ostream& out, const pt3d& pt);