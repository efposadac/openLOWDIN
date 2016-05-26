#ifndef SCF_TRANSFORM
#define SCF_TRANSFORM
#include <iostream>
#include <Lapack.h>
#include <armadillo>
#include <fstream>
#include <limits>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <cstring>
using namespace std;
using namespace arma;

inline int Multi_Index(int i, int j, int k, int l);

///This routine works perfectly!
//May  be easier try cholesky decomposition here.
void four_index_trans(int nao, mat &C, vec &EXC_ERIS);

void four_index_trans2(int nao, double C[], double ERIS[], int i_lower, int i_upper, int j_lower, int j_upper, int k_lower, int k_upper, int l_lower, int l_upper);

void four_index_trans(int nao, double C[], double ERIS[]);


#endif
