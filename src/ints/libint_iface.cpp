/*
file: libint_iface.cpp
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "libint_iface.h"

/*
Fortran interface
*/

LibintInterface *LibintInterface_new() { return new LibintInterface(); }

void LibintInterface_del(LibintInterface *lint) { lint->~LibintInterface(); }

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center) {
  lint->add_particle(z, center);
}

void LibintInterface_add_shell(LibintInterface *lint, double *alpha,
                               double *coeff, double *origin, int l,
                               int nprim) {
  lint->add_shell(alpha, coeff, origin, l, nprim);
}

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result) {
  // Default case
  libint2::Operator obtype = libint2::Operator::overlap;

  switch (integral_kind) {
  case 1: // Overlap
    obtype = libint2::Operator::overlap;
    break;
  case 2: // Kinetic
    obtype = libint2::Operator::kinetic;
    break;
  case 3: // Nuclear
    obtype = libint2::Operator::nuclear;
    break;
  }

  auto tmp = lint->compute_1body_ints(obtype);

  for (int i = 0; i < tmp.size(); ++i) {
    result[i] = tmp.array()(i);
  }
}

void LibintInterface_init_2body_ints(LibintInterface *lint) {
  lint->init_2body_ints();
}

void LibintInterface_compute_2body_fock(LibintInterface *lint, double *dens,
                                        double *result) {
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto K = compute_schwartz_ints<>(lint->get_shells());
  auto G = lint->compute_2body_fock(D, K);

  for (int i = 0; i < G.size(); ++i) {
    result[i] = G.array()(i);
  }
}

void LibintInterface_compute_2body_ints(LibintInterface *lint,
                                        const char *filename) {
  lint->compute_2body_ints(filename);
}

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
}

/*
LibintInterface class implementation
*/

LibintInterface::LibintInterface() {
  // set up thread pool
  {
    using libint2::nthreads;
    auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
    nthreads = 1;
    if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
      std::istringstream iss(nthreads_cstr);
      iss >> nthreads;
      if (nthreads > 1 << 16 || nthreads <= 0)
        nthreads = 1;
    }

#if defined(_OPENMP)
    omp_set_num_threads(nthreads);
#endif
    std::cout << "Will scale over " << nthreads
#if defined(_OPENMP)
              << " OpenMP"
#else
              << " C++11"
#endif
              << " threads" << std::endl;
  }

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();
};

void LibintInterface::init_2body_ints() {
  // compute OBS non-negligible shell-pair list
  // {
  //   obs_shellpair_list = compute_shellpair_list(shells);
  //   size_t nsp = 0;
  //   for (auto &sp : obs_shellpair_list) {
  //     nsp += sp.second.size();
  //   }
  //   std::cout << "# of {all,non-negligible} shell-pairs = {"
  //             << shells.size() * (shells.size() + 1) / 2 << "," << nsp << "}"
  //             << std::endl;
  // }
}

void LibintInterface::add_particle(const int z, const double *center) {
  // add atom information
  libint2::Atom atom;
  atom.atomic_number = z;
  atom.x = center[0];
  atom.y = center[1];
  atom.z = center[2];
  atoms.push_back(atom);
}

void LibintInterface::add_shell(double *alpha, double *coeff, double *origin,
                                int l, int nprim) {
  std::vector<double> exponents(nprim);
  std::vector<double> coefficients(nprim);

  exponents.assign(alpha, alpha + nprim);
  coefficients.assign(coeff, coeff + nprim);

  shells.push_back({{exponents},
                    {
                        {l, false, {coefficients}},
                    },
                    {{origin[0], origin[1], origin[2]}}});

  nbasis = 0;
  max_nprim = 0;
  max_l = 0;
  for (const auto &shell : shells) {
    nbasis += shell.size();
    max_nprim = std::max(shell.nprim(), max_nprim);
    for (auto c : shell.contr)
      max_l = std::max(c.l, max_l);
  }
}

std::vector<size_t> LibintInterface::map_shell_to_basis_function() {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell : shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

Matrix LibintInterface::compute_shellblock_norm(const Matrix &A) {
  const auto nsh = shells.size();
  Matrix Ash(nsh, nsh);

  auto shell2bf = map_shell_to_basis_function();
  for (size_t s1 = 0; s1 != nsh; ++s1) {
    const auto &s1_first = shell2bf[s1];
    const auto &s1_size = shells[s1].size();
    for (size_t s2 = 0; s2 != nsh; ++s2) {
      const auto &s2_first = shell2bf[s2];
      const auto &s2_size = shells[s2].size();

      Ash(s1, s2) = A.block(s1_first, s2_first, s1_size, s2_size)
                        .lpNorm<Eigen::Infinity>();
    }
  }

  return Ash;
}

Matrix LibintInterface::compute_1body_ints(libint2::Operator obtype) {
  const auto n = nbasis;
  Matrix result(n, n);

  // construct the integrals engine
  libint2::Engine engine(obtype, max_nprim, max_l, 0);

  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<double, std::array<double, 3>>> q;
    for (const auto &atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for (unsigned int s1 = 0; s1 != shells.size(); ++s1) {
    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for (unsigned int s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      const auto *buf = engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the
      // corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                    // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  return result;
}

void LibintInterface::compute_2body_ints(const char *filename) {

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  std::chrono::duration<double> time_elapsed =
      std::chrono::duration<double>::zero();

  std::cout << filename << std::endl;
  std::ofstream outfile;
  outfile.open(filename, std::ios::binary);

  const auto n = nbasis;
  Matrix G = Matrix::Zero(n, n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim, max_l, 0);

  auto shell2bf = map_shell_to_basis_function();

  // loop over permutationally-unique set of shells
  for (auto s1 = 0; s1 < shells.size(); ++s1) {
    auto n = s1;
    auto n1 = shells[s1].size(); // number of basis functions in this shell

    for (auto s2 = s1; s2 < shells.size(); ++s2) {
      auto u = s2;
      auto n2 = shells[s2].size();

      for (auto s3 = n; s3 < shells.size(); ++s3) {
        auto n3 = shells[s3].size();

        for (auto s4 = u; s4 < shells.size(); ++s4) {
          auto n4 = shells[s4].size();

          const auto tstart = std::chrono::high_resolution_clock::now();

          const auto *buf =
              engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

          const auto tstop = std::chrono::high_resolution_clock::now();

          time_elapsed += tstop - tstart;
          auto buf_size = n1 * n2 * n3 * n4;

          outfile.write((char *)&buf[0], buf_size * sizeof(buf_size));
        }
        u = s3 + 1;
      }
    }
  }
  outfile.close();
  std::cout<<"Time elapsed: "<<time_elapsed.count()<<std::endl;
}

Matrix LibintInterface::compute_2body_fock(const Matrix &D,
                                           const Matrix &Schwartz,
                                           double precision) {
  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  std::chrono::duration<double> time_elapsed =
      std::chrono::duration<double>::zero();

  const auto n = nbasis;
  Matrix G = Matrix::Zero(n, n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim, max_l, 0);

  auto shell2bf = map_shell_to_basis_function();

  // std::ofstream outfile;
  // outfile.open("lowdin.txt");
  // std::cout << "Dens: "<< D.sum()<<std::endl;
  // outfile << std::setprecision(12) << D << " ";

  // loop over permutationally-unique set of shells
  for (auto s1 = 0; s1 != shells.size(); ++s1) {
    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for (auto s2 = 0; s2 <= s1; ++s2) {
      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for (auto s3 = 0; s3 <= s1; ++s3) {
        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for (auto s4 = 0; s4 <= s4_max; ++s4) {
          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the
          // given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          const auto tstart = std::chrono::high_resolution_clock::now();

          const auto *buf =
              engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

          // std::cout << n1 << ' ' << n2 << ' ' << n3 << ' ' << n4 << ' '
          //           << std::endl;
          // std::cout << bf1_first << ' ' << bf2_first << ' '<< bf3_first << '
          // '<< bf4_first << std::endl;

          // for (int i = 0; i < n1 * n2 * n3 * n4; ++i) {
          //   std::cout << buf[i] << std::endl;
          // }

          const auto tstop = std::chrono::high_resolution_clock::now();
          time_elapsed += tstop - tstart;

          for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for (auto f2 = 0; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for (auto f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                  G(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                  G(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                  G(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                  G(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                  G(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                }
              }
            }
          }
        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  std::cout<<"Time elapsed: "<<time_elapsed.count()<<std::endl;
  return 0.25 * (G + Gt);
}

// Matrix LibintInterface::compute_2body_fock(const Matrix &D,
//                                            const Matrix &Schwartz,
//                                            double precision) {
//   const auto n = nbasis;
//   const auto nshells = shells.size();

//   using libint2::nthreads;

//   std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

//   const auto do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() !=
//   0;
//   Matrix D_shblk_norm;  // matrix of infty-norms of shell blocks
//   if (do_schwartz_screen) {
//     D_shblk_norm = compute_shellblock_norm(D);
//   }

//   auto fock_precision = precision;
//   // engine precision controls primitive truncation, assume
//   worst-casescenario
//   // (all primitive combinations add up constructively)
//   auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
//   auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
//                                    std::numeric_limits<double>::epsilon()) /
//                           max_nprim4;

//   // construct the 2-electron repulsion integrals engine pool
//   using libint2::Engine;
//   std::vector<Engine> engines(nthreads);
//   engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
//   engines[0].set_precision(engine_precision);  // shellset-dependent
//   precision
//                                                // control will likely break
//                                                // positive definiteness
//                                                // stick with this simple
//                                                recipe
//   std::cout << "compute_2body_fock:precision = " << precision << std::endl;
//   std::cout << "Engine::precision = " << engines[0].precision() << std::endl;
//   for (size_t i = 1; i != nthreads; ++i) {
//     engines[i] = engines[0];
//   }
//   std::atomic<size_t> num_ints_computed{0};

// #if defined(REPORT_INTEGRAL_TIMINGS)
//   std::vector<libint2::Timers<1>> timers(nthreads);
// #endif

//   auto shell2bf = map_shell_to_basis_function();

//   auto lambda = [&](int thread_id) {

//     auto &engine = engines[thread_id];
//     auto &g = G[thread_id];

// #if defined(REPORT_INTEGRAL_TIMINGS)
//     auto &timer = timers[thread_id];
//     timer.clear();
//     timer.set_now_overhead(25);
// #endif

//     // loop over permutationally-unique set of shells
//     for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
//       auto bf1_first = shell2bf[s1];  // first basis function in this shell
//       auto n1 = shells[s1].size();    // number of basis functions in this
//       shell

//       for (const auto &s2 : obs_shellpair_list[s1]) {
//         auto bf2_first = shell2bf[s2];
//         auto n2 = shells[s2].size();

//         const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

//         for (auto s3 = 0; s3 <= s1; ++s3) {
//           auto bf3_first = shell2bf[s3];
//           auto n3 = shells[s3].size();

//           const auto Dnorm123 =
//               do_schwartz_screen
//                   ? std::max(D_shblk_norm(s1, s3),
//                              std::max(D_shblk_norm(s2, s3), Dnorm12))
//                   : 0.;

//           const auto s4_max = (s1 == s3) ? s2 : s3;
//           for (const auto &s4 : obs_shellpair_list[s3]) {
//             if (s4 > s4_max)
//               break;  // for each s3, s4 are stored in monotonically
//               increasing
//                       // order

//             if ((s1234++) % nthreads != thread_id) continue;

//             const auto Dnorm1234 =
//                 do_schwartz_screen
//                     ? std::max(
//                           D_shblk_norm(s1, s4),
//                           std::max(D_shblk_norm(s2, s4),
//                                    std::max(D_shblk_norm(s3, s4), Dnorm123)))
//                     : 0.;

//             if (do_schwartz_screen &&
//                 Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
//                     fock_precision)
//               continue;

//             auto bf4_first = shell2bf[s4];
//             auto n4 = shells[s4].size();

//             num_ints_computed += n1 * n2 * n3 * n4;

//             // compute the permutational degeneracy (i.e. # of equivalents)
//             of
//             // the given shell set
//             auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
//             auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
//             auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
//             auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

// #if defined(REPORT_INTEGRAL_TIMINGS)
//             timer.start(0);
// #endif

//             const auto *buf = engine.compute2<libint2::Operator::coulomb,
//                                               libint2::BraKet::xx_xx, 0>(
//                 shells[s1], shells[s2], shells[s3], shells[s4]);

// #if defined(REPORT_INTEGRAL_TIMINGS)
//             timer.stop(0);
// #endif

//             // 1) each shell set of integrals contributes up to 6 shell sets
//             of
//             // the Fock matrix:
//             //    F(a,b) += (ab|cd) * D(c,d)
//             //    F(c,d) += (ab|cd) * D(a,b)
//             //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
//             //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
//             //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
//             //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
//             // 2) each permutationally-unique integral (shell set) must be
//             // scaled by its degeneracy,
//             //    i.e. the number of the integrals/sets equivalent to it
//             // 3) the end result must be symmetrized
//             for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
//               const auto bf1 = f1 + bf1_first;
//               for (auto f2 = 0; f2 != n2; ++f2) {
//                 const auto bf2 = f2 + bf2_first;
//                 for (auto f3 = 0; f3 != n3; ++f3) {
//                   const auto bf3 = f3 + bf3_first;
//                   for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
//                     const auto bf4 = f4 + bf4_first;

//                     const auto value = buf[f1234];

//                     const auto value_scal_by_deg = value * s1234_deg;

//                     g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
//                     g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
//                     g(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
//                     g(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
//                     g(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
//                     g(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }

//   };  // end of lambda

//   libint2::parallel_do(lambda);

//   // accumulate contributions from all threads
//   for (size_t i = 1; i != nthreads; ++i) {
//     G[0] += G[i];
//   }

// #if defined(REPORT_INTEGRAL_TIMINGS)
//   double time_for_ints = 0.0;
//   for (auto &t : timers) {
//     time_for_ints += t.read(0);
//   }
//   std::cout << "time for integrals = " << time_for_ints << std::endl;
//   for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
// #endif

//   Matrix GG = 0.25 * (G[0] + G[0].transpose());

//   std::cout << "# of integrals = " << num_ints_computed << std::endl;

//   // symmetrize the result and return
//   return GG;
// }

/*
Auxiliary functions
*/
shellpair_list_t compute_shellpair_list(const std::vector<libint2::Shell> bs1,
                                        const std::vector<libint2::Shell> &_bs2,
                                        const double threshold) {
  const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  libint2::BasisSet obs_aux;
  std::vector<Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(libint2::Operator::overlap,
                       std::max(obs_aux.max_nprim(bs1), obs_aux.max_nprim(bs2)),
                       std::max(obs_aux.max_l(bs1), obs_aux.max_l(bs2)), 0);

  for (size_t i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&](int thread_id) {

    // auto &engine = engines[thread_id];

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      mx.lock();
      if (result.find(s1) == result.end())
        result.insert(std::make_pair(s1, std::vector<size_t>()));
      mx.unlock();

      auto n1 = bs1[s1].size(); // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id)
          continue;

        auto on_same_center = (bs1[s1].O == bs2[s2].O);
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = bs2[s2].size();
          const auto *buf = engines[thread_id].compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          result[s1].emplace_back(s2);
          mx.unlock();
        }
      }
    }
  }; // end of compute

  libint2::parallel_do(compute);

  // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
  // if s1 < s2
  auto sort = [&](int thread_id) {
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % nthreads == thread_id) {
        auto &list = result[s1];
        std::sort(list.begin(), list.end());
      }
    }
  }; // end of sort

  libint2::parallel_do(sort);

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return result;
}

template <libint2::Operator Kernel>
Matrix compute_schwartz_ints(
    const std::vector<libint2::Shell> bs1,
    const std::vector<libint2::Shell> &_bs2, bool use_2norm,
    typename libint2::operator_traits<Kernel>::oper_params_type params) {
  const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1, nsh2);

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  using libint2::nthreads;
  std::vector<Engine> engines(nthreads);

  // !!! very important: cannot screen primitives in Schwartz computation !!!
  libint2::BasisSet obs_aux;
  auto epsilon = 0.;
  engines[0] = Engine(
      Kernel, std::max(obs_aux.max_nprim(bs1), obs_aux.max_nprim(bs2)),
      std::max(obs_aux.max_l(bs1), obs_aux.max_l(bs2)), 0, epsilon, params);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  // std::cout << "computing Schwartz bound prerequisites (kernel=" <<
  // (int)Kernel
  // << ") ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  auto lambda = [&](int thread_id) {

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      auto n1 = bs1[s1].size(); // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id)
          continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1 * n2;

        const auto *buf =
            engines[thread_id].compute(bs1[s1], bs2[s2], bs1[s1], bs2[s2]);

        // the diagonal elements are the Schwartz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf, n12, n12);
        auto norm2 = use_2norm ? buf_mat.diagonal().norm()
                               : buf_mat.diagonal().lpNorm<Eigen::Infinity>();
        K(s1, s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2)
          K(s2, s1) = K(s1, s2);
      }
    }
  }; // thread lambda

  libint2::parallel_do(lambda);

  timer.stop(0);
  // std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return K;
}
