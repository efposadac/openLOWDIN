/**
 * @file transform
 * @author:  Ruben Dario Guerrero <rdguerrerom@unal.edu.co>
 * @version 7.0
 *
 *     This file is part of the HELIOS ab-initio dynamics package
 *
 * @section LICENSE
 *     this program is being developed under the advice of:
 *     Prof. A. Reyes(a) and Prof. J. Vanicek(b)
 *      (a) Universidad Nacional de Colombia
 *      (b) Swiss Federal Institute of Technology in Lausanne - EPFL
 *
 *     Copyright R.D Guerrero 2016
 *     No portion of this file may be distributed, modified, or used without
 express
 *     written permission from  the author.
 *
 * @section DESCRIPTION

 */

#ifndef TRANSFORM
#define TRANSFORM

#include <Eigen/Eigenvalues>
#include <assert.h>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <limits>
#include <stdlib.h>

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    mat;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vec;

inline int Multi_Index(int i, int j, int k, int l);

/// This routine works perfectly!
// May  be easier try cholesky decomposition here.
void four_index_trans(int nao, mat &C, vec &EXC_ERIS);
void four_index_trans2(int nao, double C[], double ERIS[], int i_lower,
                       int i_upper, int j_lower, int j_upper, int k_lower,
                       int k_upper, int l_lower, int l_upper);

void four_index_trans(int nao, double C[], double ERIS[]);

void four_index_trans_inter(int nao, int onao, double C[], double OC[], double ERIS[]);

#ifdef __cplusplus
extern "C" {
#endif

// TODO: It does not work for now.
void c_integrals_transform_partial(double *coeff, double *ints, int nao, int lp,
                                   int up, int lq, int uq, int lr, int ur,
                                   int ls, int us);

void c_integrals_transform_all(double *coeff, double *ints, int nao);

void c_integrals_transform_inter_all(double *coeff, double *ocoeff, double *ints, int nao, int onao);

void dgemm_(char *transA, char *transB, int *m, int *n, int *k, double *alpha,
            double *A, int *lda, double *B, int *ldb, double *beta, double *C,
            int *ldc);

void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda,
            double *X, int *incx, double *beta, double *Y, int *incy);

void Similarity_Transform(unsigned n, double Op[], double C[],
                                 double Work[]);

#ifdef __cplusplus
}
#endif

#endif
