/*
file: Libint2Iface.h
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#ifndef LIBINT2IFACE_H
#define LIBINT2IFACE_H

// #include <libint2.h>
#include <libint2.hpp>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdlib.h>
#include <thread>
#include <unordered_map>

#ifdef _OMP
#include <omp.h>
#endif

#include "Iterators.h"

/*
Raw libint2 interface
*/

#define VECLEN LIBINT2_MAX_VECLEN

typedef struct {
  double LIBINT_T_SS_EREP_SS0;
  double LIBINT_T_SS_EREP_SS1;
  double LIBINT_T_SS_EREP_SS2;
  double LIBINT_T_SS_EREP_SS3;
  double LIBINT_T_SS_EREP_SS4;
  double LIBINT_T_SS_EREP_SS5;
  double LIBINT_T_SS_EREP_SS6;
  double LIBINT_T_SS_EREP_SS7;
  double LIBINT_T_SS_EREP_SS8;
  double LIBINT_T_SS_EREP_SS9;
  double LIBINT_T_SS_EREP_SS10;
  double LIBINT_T_SS_EREP_SS11;
  double LIBINT_T_SS_EREP_SS12;
  double LIBINT_T_SS_EREP_SS13;
  double LIBINT_T_SS_EREP_SS14;
  double LIBINT_T_SS_EREP_SS15;
  double LIBINT_T_SS_EREP_SS16;
  double LIBINT_T_SS_EREP_SS17;
  double LIBINT_T_SS_EREP_SS18;
  double LIBINT_T_SS_EREP_SS19;
  double LIBINT_T_SS_EREP_SS20;

  // Prefactors for recurrence relations from Weber and Daul, Comp. Phys. Comm.
  // 158, 1 (2004).

  double LIBINT_T_SS_K0G12_SS_0;
  double LIBINT_T_SS_K2G12_SS_0;
  double LIBINT_T_SS_Km1G12_SS0;
  double LIBINT_T_SS_Km1G12_SS1;
  double LIBINT_T_SS_Km1G12_SS2;
  double LIBINT_T_SS_Km1G12_SS3;
  double LIBINT_T_SS_Km1G12_SS4;
  double LIBINT_T_SS_Km1G12_SS5;
  double LIBINT_T_SS_Km1G12_SS6;
  double LIBINT_T_SS_Km1G12_SS7;
  double LIBINT_T_SS_Km1G12_SS8;
  double LIBINT_T_SS_Km1G12_SS9;
  double LIBINT_T_SS_Km1G12_SS10;
  double LIBINT_T_SS_Km1G12_SS11;
  double LIBINT_T_SS_Km1G12_SS12;
  double LIBINT_T_SS_Km1G12_SS13;
  double LIBINT_T_SS_Km1G12_SS14;
  double LIBINT_T_SS_Km1G12_SS15;
  double LIBINT_T_SS_Km1G12_SS16;

  // LRL1991, Eq. 30, prefactor in front of (a0|c0)
  double TwoPRepITR_pfac0_0_x;
  double TwoPRepITR_pfac0_0_y;
  double TwoPRepITR_pfac0_0_z;
  double TwoPRepITR_pfac0_1_x;
  double TwoPRepITR_pfac0_1_y;
  double TwoPRepITR_pfac0_1_z;

  // LRL1991, Eq. 30, prefactor in front of (a0|c+10)
  double TwoPRepITR_pfac1_0;
  double TwoPRepITR_pfac1_1;

  // WD2004, Eq. 30, prefactor in front of (a0|k|c0)
  double R12kG12_pfac0_0_x;
  double R12kG12_pfac0_0_y;
  double R12kG12_pfac0_0_z;
  double R12kG12_pfac0_1_x;
  double R12kG12_pfac0_1_y;
  double R12kG12_pfac0_1_z;

  // WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0)
  double R12kG12_pfac1_0;
  double R12kG12_pfac1_1;

  // WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0)
  double R12kG12_pfac2;

  // WD2004, Eq. 30, prefactor in front of curly brakets (excludes k)
  double R12kG12_pfac3_0;
  double R12kG12_pfac3_1;

  // WD2004, Eq. 30, prefactor in front of (a0|k-2|c0)
  double R12kG12_pfac4_0_x;
  double R12kG12_pfac4_0_y;
  double R12kG12_pfac4_0_z;
  double R12kG12_pfac4_1_x;
  double R12kG12_pfac4_1_y;
  double R12kG12_pfac4_1_z;

  // Exponents
  double zeta_A;
  double zeta_B;
  double zeta_C;
  double zeta_D;

  // Squared exponents
  double zeta_A_2;
  double zeta_B_2;
  double zeta_C_2;
  double zeta_D_2;

  // Appear in OS RR for ERIs

  // One over 2.0*zeta
  double oo2z;

  // One over 2.0*eta
  double oo2e;

  // One over 2.0*(zeta+eta)
  double oo2ze;

  // rho over zeta
  double roz;

  // rho over eta
  double roe;

  // Appear in standard OS RR for ERI and almost all other recurrence relations
  double WP_x, WP_y, WP_z;
  double WQ_x, WQ_y, WQ_z;
  double PA_x, PA_y, PA_z;
  double QC_x, QC_y, QC_z;
  double AB_x, AB_y, AB_z;
  double CD_x, CD_y, CD_z;

} lowdin_t;

/*
Namespace extension for threads
*/

namespace libint2 {
int nthreads;

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda> void parallel_do(Lambda &lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else // use C++11 threads
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  } // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}

struct QuartetBuffer {
  std::vector<int> p;
  std::vector<int> q;
  std::vector<int> r;
  std::vector<int> s;
  std::vector<double> val;
};

} // end libint2 namespace

/*
Type definitions
*/
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix; // import dense, dynamically sized Matrix type from Eigen;
            // this is a matrix with row-major storage
            // (http://en.wikipedia.org/wiki/Row-major_order)

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

/*
Auxiliary functions
*/

template <libint2::Operator Kernel = libint2::Operator::coulomb>
Matrix compute_schwartz_ints(
    std::vector<libint2::Shell> bs1,
    const std::vector<libint2::Shell> &bs2 = std::vector<libint2::Shell>(),
    bool use_2norm = false, // use infty norm by default
    typename libint2::operator_traits<Kernel>::oper_params_type params =
        libint2::operator_traits<Kernel>::default_params());

shellpair_list_t compute_shellpair_list(
    const std::vector<libint2::Shell> &bs1,
    const std::vector<libint2::Shell> &_bs2 = std::vector<libint2::Shell>(),
    double threshold = 1e-12);

__inline__ void write_buffer(const libint2::QuartetBuffer &buffer,
                             const int &s_size, std::ofstream &outfile);
/*
Libint interface class for a single species
*/
class LibintInterface {
private:
  size_t max_nprim;
  size_t nbasis;
  int s_size; // stack size
  int max_l;
  int speciesID;
  bool is_electron;
  shellpair_list_t obs_shellpair_list;
  std::vector<libint2::Atom> atoms;
  std::vector<libint2::Shell> shells;
  std::vector<double> norma;
  Matrix compute_shellblock_norm(const Matrix &A);

public:

  LibintInterface(const int stack_size, const int id, const bool el);

  ~LibintInterface() { libint2::finalize(); };

  void add_particle(const int z, const double *center);

  void add_shell(double *alpha, double *coeff, double *origin, int l,
                 int nprim);

  Matrix compute_1body_ints(libint2::Operator obtype);

  void init_2body_ints();

  void
  compute_2body_disk(const char *filename, const Matrix &D,
                     const Matrix &Schwartz,
                     double precision = std::numeric_limits<double>::epsilon());

  Matrix compute_2body_direct(
      const Matrix &D, const Matrix &Schwartz,
      double precision = std::numeric_limits<double>::epsilon());

  void compute_coupling_disk(
      LibintInterface &other, const char *filename,
      double precision = std::numeric_limits<double>::epsilon());

  Matrix compute_coupling_direct(
      LibintInterface &other, const Matrix &D, const bool permuted,
      double precision = std::numeric_limits<double>::epsilon());

  void compute_g12_disk(const char *filename,
                                         const double *coefficients,
                                         const double *exponents,
                                         const int pot_size);

  std::vector<size_t> map_shell_to_basis_function();

  std::vector<libint2::Shell> get_shells() { return shells; };

  double get_norma(int index) { return norma[index]; };
  size_t get_max_nprim() { return max_nprim; };
  size_t get_nbasis() { return nbasis; };
  int get_max_l() { return max_l; };
  int get_speciesID() { return speciesID; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
Fortran interface routines.
*/
LibintInterface *LibintInterface_new(const int stack_size, const int id,
                                     const bool el);

void LibintInterface_del(LibintInterface *lint);

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center);

void LibintInterface_add_shell(LibintInterface *lint, double *alpha,
                               double *coeff, double *origin, int l, int nprim);

int LibintInterface_get_nbasis(LibintInterface *lint);

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result);

void LibintInterface_init_2body_ints(LibintInterface *lint);

void LibintInterface_compute_2body_direct(LibintInterface *lint, double *dens,
                                          double *result);

void LibintInterface_compute_2body_disk(LibintInterface *lint,
                                        const char *filename, double *dens);

void LibintInterface_compute_coupling_direct(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *result);

void LibintInterface_compute_coupling_disk(LibintInterface *lint,
                                           LibintInterface *olint,
                                           const char *filename);

void libintinterface_compute_g12_disk(LibintInterface *lint,
                                      const char *filename,
                                      const double *coefficients,
                                      const double *exponents,
                                      const int pot_size);

void libintinterface_buildg12_(int *, int *, int *, int *, int *, int *,
                               lowdin_t *, double *);

void LibintInterface_setLibint(Libint_t *, lowdin_t *);

#ifdef __cplusplus
}
#endif

#endif
