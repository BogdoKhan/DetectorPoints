#pragma once
#include <vector>

using namespace std;

enum class posPlane {
	XZ, YZ, XY, NA
};

//tangent to two circles and corresponding structures
//equation is taken from https://e-maxx.ru/algo/circle_tangents
struct pt {
	double crd1, crd2;

	pt operator- (pt p) {
		pt res = { crd1 - p.crd1, crd2 - p.crd2 };
		return res;
	}
};

struct circle : pt {
public:
	double r = 0.0;
	posPlane orientation = posPlane::NA;
};

struct pt3d {
	double x, y, z;
	pt3d operator- (pt3d p) {
		pt3d res = { x - p.x, y - p.y, z - p.z };
		return res;
	}
};

struct cylinder : circle {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
public:
	cylinder(circle c, double crd3) {
		circle::r = c.r;
		if (c.orientation == posPlane::XZ) {
			x = c.crd1;
			y = crd3;
			z = c.crd2;
			circle::orientation = c.orientation;
		}
		else if (c.orientation == posPlane::YZ) {
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

struct line {
	double a, b, c;
};

const double EPS = 1E-9;

double sqr(double a) {
	return a * a;
}

void tangents(pt c, double r1, double r2, vector<line>& ans) {
	double r = r2 - r1;
	double z = sqr(c.crd1) + sqr(c.crd2);
	double d = z - sqr(r);
	if (d < -EPS)  return;
	d = sqrt(abs(d));
	line l;
	l.a = (c.crd1 * r + c.crd2 * d) / z;
	l.b = (c.crd2 * r - c.crd1 * d) / z;
	l.c = r1;
	ans.push_back(l);
}

vector<line> tangents(circle a, circle b) {
	vector<line> ans;
	for (int i = -1; i <= 1; i += 2)
		for (int j = -1; j <= 1; j += 2)
			tangents(b - a, a.r * i, b.r * j, ans);
	for (size_t i = 0; i < ans.size(); ++i)
		ans[i].c -= ans[i].a * a.crd1 + ans[i].b * a.crd2;
	return ans;
}

struct plane {
	double a;
	double b;
	double c;
	double d;

	plane() : a(0), b(0), c(0), d(0) {};
	plane(double _a, double _b, double _c, double _d) :
		a(_a), b(_b), c(_c), d(_d) {};
};
