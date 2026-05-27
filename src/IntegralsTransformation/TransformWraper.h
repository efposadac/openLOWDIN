/********************************************************-
  This code is part of LOWDIN Quantum chemistry package

  this program has been developed under direction of:

  Prof. A REYES' Lab. Universidad Nacional de Colombia
    http://www.qcc.unal.edu.co
  Prof. R. FLORES' Lab. Universidad de Guadalajara
    http://www.cucei.udg.mx/~robertof

    Todos los derechos reservados, 2013

**********************************************************/

#ifndef TRANSFORMWRAPER
#define TRANSFORMWRAPER

#include <iostream>
#include <cstdio>
#include <Transform.h>

#ifdef __cplusplus
extern "C" {
#endif

void c_integrals_transform(double *coeff, double *ints, int nao, int lp, int up,
                           int lq, int uq, int lr, int ur, int ls, int us);

#ifdef __cplusplus
}
#endif

#endif
