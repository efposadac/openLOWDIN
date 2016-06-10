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
    const std::vector<libint2::Shell> bs1,
    const std::vector<libint2::Shell> &_bs2 = std::vector<libint2::Shell>(),
    double threshold = 1e-12);

/*
Libint interface class
*/
class LibintInterface {
private:
  size_t max_nprim;
  size_t nbasis;
  int s_size; // stack size
  int max_l;
  shellpair_list_t obs_shellpair_list;
  std::vector<libint2::Atom> atoms;
  std::vector<libint2::Shell> shells;
  std::vector<double> norma;
  std::vector<size_t> map_shell_to_basis_function();
  Matrix compute_shellblock_norm(const Matrix &A);

public:
  LibintInterface(const int stack_size);
  ~LibintInterface() { libint2::finalize(); };
  void add_particle(const int z, const double *center);
  void add_shell(double *alpha, double *coeff, double *origin, int l,
                 int nprim);
  void init_2body_ints();
  void
  compute_2body_ints(const char *filename, const Matrix &D,
                     const Matrix &Schwartz,
                     double precision = std::numeric_limits<double>::epsilon());
  Matrix compute_1body_ints(libint2::Operator obtype);
  Matrix
  compute_2body_fock(const Matrix &D, const Matrix &Schwartz,
                     double precision = std::numeric_limits<double>::epsilon());
  std::vector<libint2::Shell> get_shells() { return shells; };
  size_t get_max_nprim() { return max_nprim; };
  size_t get_nbasis() { return nbasis; };
  int get_max_l() { return max_l; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
Fortran interface routines.
*/
LibintInterface *LibintInterface_new(int stack_size);

void LibintInterface_del(LibintInterface *lint);

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center);

void LibintInterface_add_shell(LibintInterface *lint, double *alpha,
                               double *coeff, double *origin, int l, int nprim);

int LibintInterface_get_nbasis(LibintInterface *lint);

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result);

void LibintInterface_init_2body_ints(LibintInterface *lint);

void LibintInterface_compute_2body_fock(LibintInterface *lint, double *dens,
                                        double *result);

void LibintInterface_compute_2body_ints(LibintInterface *lint,
                                        const char *filename, double *dens);

#ifdef __cplusplus
}
#endif

#endif