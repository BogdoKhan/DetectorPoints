#pragma once

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>

int circleAndLine(const gsl_vector* variables, void* params,
	gsl_vector* f);
int print_state(size_t iter, gsl_multiroot_fsolver* s);
int lineAndPlane(const gsl_vector* variables, void* params,
	gsl_vector* f);
int print_state_linePl(size_t iter, gsl_multiroot_fsolver* s);

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

