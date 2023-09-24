#include "Tangents.h"
#include "Logger.h"

pt pt::operator- (pt p) {
	pt res = { crd1 - p.crd1, crd2 - p.crd2 };
	return res;
}

circle::circle(double crd1, double crd2, double r) {
	pt::crd1 = crd1;
	pt::crd2 = crd2;
	circle::r = r;
};

circle::circle() {
	pt::crd1 = 0.0;
	pt::crd2 = 0.0;
	circle::r = 0.0;
}

pt3d::pt3d() {
	x = 0; y = 0; z = 0;
}

pt3d::pt3d(double x_, double y_, double z_) {
	x = x_; y = y_; z = z_;
}

pt3d pt3d::operator- (pt3d p) {
	pt3d res((x - p.x), (y - p.y), (z - p.z));
	return res;
}

cylinder::cylinder(circle c, double crd3) {
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

cylinder::cylinder() {
	circle::r = 0.0;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	circle::orientation = posPlane::NA;
}

double sqr(double a) {
	return a * a;
}

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

std::ostream& operator<< (std::ostream& out, const plane& pp) {
	out << "Plane equation: " << pp.a;
	(pp.b >= 0) ? out << "x + " : out << "x "; out << pp.b;
	(pp.c >= 0) ? out << "y + " : out << "y "; out << pp.c;
	(pp.d >= 0) ? out << "z + " : out << "z "; out << pp.d;
	return out;
}

std::ostream& operator<< (std::ostream& out, const pt3d& pt) {
	out << "Point at cooridnates: {" << pt.x;
	out << ", " << pt.y << ", " << pt.z << "}";
	return out;
}

std::ostream& PrintTangentLine (std::ostream& out, const line& line, const cylinder& cyl) {

	if (cyl.orientation == posPlane::XZ) {
		out << "Line equation: " << line.a;
		(line.b >= 0) ? out << "x+" : out << "x"; out << line.b;
		(line.c >= 0) ? out << "z+" : out << "z"; out << line.c;
		out << " == 0";
	}
	else {
		out << "Line equation: " << line.a;
		(line.b >= 0) ? out << "y+" : out << "y"; out << line.b;
		(line.c >= 0) ? out << "z+" : out << "z"; out << line.c;
		out << " == 0";
	}

	return out;
}