#pragma once

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>

int circleAndLine(const gsl_vector* variables, void* params,
	gsl_vector* f);
int print_state(size_t iter, gsl_multiroot_fsolver* s);

struct params_tangentCrd
{
	double crd1;
	double crd2;
	double r;
	double a;
	double b;
	double c;
};

int circleAndLine(const gsl_vector* variables, void* params,
	gsl_vector* f) {

	double crd1 = ((struct params_tangentCrd*)params)->crd1;
	double crd2 = ((struct params_tangentCrd*)params)->crd2;
	double r = ((struct params_tangentCrd*)params)->r;
	double a = ((struct params_tangentCrd*)params)->a;
	double b = ((struct params_tangentCrd*)params)->b;
	double c = ((struct params_tangentCrd*)params)->c;
	const double x = gsl_vector_get(variables, 0);
	const double y = gsl_vector_get(variables, 1);
	const double f0 = gsl_pow_2(x - crd1) + gsl_pow_2(y - crd2) - gsl_pow_2(r);
	const double f1 = a * x + b * y + c;
	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	return GSL_SUCCESS;
}

int print_state(size_t iter, gsl_multiroot_fsolver* s)
{
	printf("iter = %3u x = % .6f % .6f "
		"f(x) = % .6e % .6e\n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, 1));
	return 0;
}

struct params_linPlane {
	double a;
	double b;
	double c;
	double d;
	double zz;
};

int lineAndPlane(const gsl_vector* variables, void* params,
	gsl_vector* f) {
	double a = ((struct params_linPlane*)params)->a;
	double b = ((struct params_linPlane*)params)->b;
	double c = ((struct params_linPlane*)params)->c;
	double d = ((struct params_linPlane*)params)->d;
	double zz = ((struct params_linPlane*)params)->zz;
	const double x = gsl_vector_get(variables, 0);
	const double y = gsl_vector_get(variables, 1);
	const double z = gsl_vector_get(variables, 2);
	const double f0 = a*x + b*y + c*z + d;
	const double f1 = z - zz;
	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	return GSL_SUCCESS;
}

int print_state_linePl(size_t iter, gsl_multiroot_fsolver* s)
{
	printf("iter = %3u x = % .6f % .6f "
		"f(x) = % .6e % .6e\n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->x, 2),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, 1));
	return 0;
}

