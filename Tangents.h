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
	pt operator- (pt p) {
		pt res = { crd1 - p.crd1, crd2 - p.crd2 };
		return res;
	}
};
//a circle, has radius r and central point at crd1, crd2
struct circle : pt {
public:
	double r = 0.0;
	posPlane orientation = posPlane::NA;	//orientation of detector cross section: XZ, YZ, XY
	circle(double crd1, double crd2, double r) {
		pt::crd1 = crd1;
		pt::crd2 = crd2;
		circle::r = r;
	};
	circle() {
		pt::crd1 = 0.0;
		pt::crd2 = 0.0;
		circle::r = 0.0;
	}
};

struct pt3d {
	double x, y, z;	//coordinates in space based on cross sectional orientation
	pt3d() {
		x = 0; y = 0; z = 0;
	}
	pt3d(double x_, double y_, double z_) {
		x = x_; y = y_; z = z_;
	}
	pt3d operator- (pt3d p) {
		pt3d res( (x - p.x), (y - p.y), (z - p.z) );
		return res;
	}

};
//cylinder is based on circle
struct cylinder : circle {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
public:
	cylinder(circle c, double crd3) {
		circle::r = c.r;
		if (c.orientation == posPlane::XZ) {	//for XZ orientation
			x = c.crd1;
			y = crd3;
			z = c.crd2;
			circle::orientation = c.orientation;
		}
		else if (c.orientation == posPlane::YZ) { //for YZ orientation
			x = crd3;
			y = c.crd1;
			z = c.crd2;
			circle::orientation = c.orientation;
		}
	};
	cylinder() {
		circle::r = 0.0;
		x = 0.0;
		y = 0.0;
		z = 0.0;
		circle::orientation = posPlane::NA;
	}
};
//common equation of line
struct line {
	double a, b, c;
};

const double EPS = 1E-9;//epsilon, defines how accurate we want to be

double sqr(double a) {
	return a * a;
}
//finds inner and outer tangents to 2 circles
void tangents(pt c, double r1, double r2, vector<line>& ans) {
	double r = r2 - r1;
	double z = sqr(c.crd1) + sqr(c.crd2);
	double d = z - sqr(r);
	if (abs(d) < -EPS)  return;
	d = sqrt(abs(d));
	line l;
	l.a = (c.crd1 * r + c.crd2 * d) / z;
	l.b = (c.crd2 * r - c.crd1 * d) / z;
	l.c = r1;
	ans.push_back(l);
}
//creates vector of 4 tangents
vector<line> tangents(circle a, circle b) {
	vector<line> ans;
	for (int i = -1; i <= 1; i += 2)
		for (int j = -1; j <= 1; j += 2)
			tangents(b - a, a.r * i, b.r * j, ans);
	for (size_t i = 0; i < ans.size(); ++i)
		ans[i].c -= (ans[i].a * a.crd1 + ans[i].b * a.crd2);
	return ans;
}

vector<line> tangents(cylinder cyl1, cylinder cyl2) {
	circle a = circle(cyl1.crd1, cyl1.crd2, cyl1.r);
	circle b = circle(cyl2.crd1, cyl2.crd2, cyl2.r);
	vector<line> ans = tangents(a, b);
	return ans;
}
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

std::ostream& operator<< (std::ostream& out, const plane& pp) {
	out << "Plane equation: " << pp.a;
	(pp.b >= 0) ? out << "x + " : out << "x "; out << pp.b;
	(pp.c >= 0) ? out << "y + " : out << "y "; out << pp.c;
	(pp.d >= 0) ? out << "z + " : out << "z "; out << pp.d;
	return out;
}

std::ostream& operator<< (std::ostream& out, const pt3d& pt) {
	out << "point at cooridnates: {" << pt.x;
	out << ", " << pt.y << ", " << pt.z << "}";
	return out;
}