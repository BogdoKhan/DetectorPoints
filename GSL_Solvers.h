#pragma once

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>

#include <iostream>

#include "Logger.h"

int circleAndLine(const gsl_vector* variables, void* params,
	gsl_vector* f);
int print_state(size_t iter, gsl_multiroot_fsolver* s, std::ostream& out);
int lineAndPlane(const gsl_vector* variables, void* params,
	gsl_vector* f);
int print_state_linePl(size_t iter, gsl_multiroot_fsolver* s, std::ostream& out);

struct params_tangentCrd
{
	double crd1;
	double crd2;
	double r;
	double a;
	double b;
	double c;
};


struct params_linPlane {
	double a;
	double b;
	double c;
	double d;
	double zz;
};

