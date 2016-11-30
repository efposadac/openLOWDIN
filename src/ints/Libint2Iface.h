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
Namespace extension for threads
*/

namespace libint2 {
int nthreads;

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda &lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else  // use C++11 threads
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  }  // threads_id
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

}  // end libint2 namespace

/*
Type definitions
*/
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;  // import dense, dynamically sized Matrix type from Eigen;
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
    bool use_2norm = false,  // use infty norm by default
    typename libint2::operator_traits<Kernel>::oper_params_type params =
        libint2::operator_traits<Kernel>::default_params());

shellpair_list_t compute_shellpair_list(
    const std::vector<libint2::Shell> & bs1,
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
  int s_size;  // stack size
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

  void compute_2body_disk(
      const char *filename, const Matrix &D, const Matrix &Schwartz,
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
LibintInterface *LibintInterface_new(const int stack_size, const int id, const bool el);

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

#ifdef __cplusplus
}
#endif

#endif