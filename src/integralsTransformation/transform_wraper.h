#ifndef TRANSFORM_WRAPER
#define TRANSFORM_WRAPER

#include <iostream>
#include <cstdio>
#include <transform.h>

void print_matrix( char* desc, int m, int n, double* a, int lda );

#ifdef __cplusplus
extern "C" {
#endif

void c_test(double *coeff, double *ints, int nao);

#ifdef __cplusplus
}
#endif

#endif
