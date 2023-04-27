////There are parts of 
////magic equations for coord1 (let it be x) and coord2 (let it be y)
////calcualted with Wolfram Mathematica
////Solve[{(x - x1)^2 + (y - y1)^2 == r1^2, 
////y == -((d*x + e)/f)}, {x, y}]
////this is the same as
////Solve[{(x - a)^2 + (y - b)^2 == r^2, 
////d*x + f*y + e == 0}, {x, y}]
//
////First coordinate (x or Y) for tangence point equation
//double getCoord1_TangPt(const double& x1, const double& y1, const double& r1,
//	const double& d, const double& e, const double& f) {
//	double x_part1 = 1 / (sqr(d) + sqr(f)) * 0.5;
//	double x_part2 = -2 * d * e - 2 * y1 * d * f + 2 * x1 * sqr(f);
//	double x_part3 = sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f));
//	double x_part4 = 4 * (sqr(d) + sqr(f));
//	double x_part5 = (sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
//		sqr(y1) * sqr(f) - sqr(f) * sqr(r1));
//	double value = (x_part1 * (x_part2 - sqrt(abs(x_part3 - x_part4 * x_part5))));
//	return value;
//}
////2nd coordiante (typically Z) for tangence point equation
//double getCoord2_TangPt(const double& x1, const double& y1, const double& r1,
//	const double& d, const double& e, const double& f) {
//	double y_part1 = (1 / f);
//	double y_part2 = -e;
//	double y_part3 = (sqr(d) * e) / (sqr(d) + sqr(f));
//	double y_part4 = (y1 * sqr(d) * f) / (sqr(d) + sqr(f));
//	double y_part5 = (x1 * d * sqr(f)) / (sqr(d) + sqr(f));
//	double y_part6 = 1 / (2 * (sqr(d) + sqr(f))) * d;
//	double y_part7 = sqr(2 * d * e + 2 * y1 * d * f - 2 * x1 * sqr(f));
//	double y_part8 = 4 * (sqr(d) + sqr(f));
//	double y_part9 = (sqr(e) + 2 * y1 * e * f + sqr(x1) * sqr(f) +
//		sqr(y1) * sqr(f) - sqr(f) * sqr(r1));
//	double value = y_part1 * (y_part2 + y_part3 + y_part4 - y_part5 + y_part6 *
//		sqrt(abs(y_part7 - y_part8 * y_part9)));
//	return value;
//}
//
////(x-x1)^2 + (y-y1)^2 == R1^2;
////dx + fy + e = 0
////returns coordinates of tangence point between line & circle1
////(x_tang, y_tang)
////!!DELETE flag_round, check
////pt tangentCrd(const double& _x1, const double& _y1, const double& _r1,
////	const double& _d, const double& _e, const double& _f, bool flag_round) {
////	double coef = 1.0e6; //round to 6th digit, else goes to nan
////	double x1 = std::round(_x1 * coef) / coef;
////	double y1 = std::round(_y1 * coef) / coef;
////	double r1 = std::round(_r1 * coef) / coef;
////	double d = std::round(_d * coef) / coef;
////	double e = std::round(_e * coef) / coef;
////	double f = std::round(_f * coef) / coef;
////
////	double resultX = 0.0;
////	double resultY = 0.0;
////	pt result = {};
////
////	// 
////	//tangence point coordinates
////	result = { getCoord1_TangPt(x1, y1, r1, d, e, f), getCoord2_TangPt(x1, y1, r1, d, e, f) };
////
////	//if we get NaN - try to get solution with ceil rounded values
////	//MAYBE NOT NEEDED
////	if ((isnan(resultX) || isnan(resultY)) && !flag_round) {
////		x1 = std::ceil(_x1 * coef) / coef;
////		y1 = std::ceil(_y1 * coef) / coef;
////		r1 = (std::round(_r1 * coef) / coef);
////		d = std::ceil(_d * coef) / coef;
////		e = std::ceil(_e * coef) / coef;
////		f = std::ceil(_f * coef) / coef;
////		flag_round = true;
////		result = tangentCrd(x1, y1, r1, d, e, f, true);
////
////	}
////	//DEPRECATED, CHECK AND DELETE IF NOT NECESSARY
////	//if (isnan(result.crd1) || isnan(result.crd2)) {
////	//	x1 = std::floor(_x1 * coef) / coef;
////	//	y1 = std::floor(_y1 * coef) / coef;
////	//	r1 = std::floor(_r1 * coef) / coef;
////	//	d = std::floor(_d * coef) / coef;
////	//	e = std::floor(_e * coef) / coef;
////	//	f = std::floor(_f * coef) / coef;
////	//	result = tangentCrd(x1, y1, r1, d, e, f, true);
////	//}
////	return result;
////}