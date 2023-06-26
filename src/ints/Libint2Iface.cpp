/*
file: Libint2Iface.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "Libint2Iface.h"
#include <libint2/util/small_vector.h>

/*
LibintInterface class implementation
*/

LibintInterface::LibintInterface(const int stack_size, const int id,
                                 const bool el, const bool parallel)
    : max_nprim(0), nbasis(0), s_size(stack_size), max_l(0), speciesID(id),
      is_electron(el) {
  // set up thread pool
  using libint2::nthreads;
  nthreads = 1;
  if (parallel) {
    auto nthreads_cstr = getenv("OMP_NUM_THREADS");
    if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
      std::istringstream iss(nthreads_cstr);
      iss >> nthreads;
      if (nthreads > 1 << 16 || nthreads <= 0)
        nthreads = 1;
    }
  }
  {
#if defined(_OPENMP)
    omp_set_num_threads(nthreads);
#endif
    // std::cout << "Will scale over " << nthreads
    // #if defined(_OPENMP)
    //               << " OpenMP"
    // #else
    //               << " C++11"
    // #endif
    //               << " threads" << std::endl;
  }

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();
  name = "Thread:"+(std::to_string(omp_get_thread_num()+1))+"Species"+(std::to_string(id));

  // std::cout << "libint2::initialized() in constructor" << libint2::initialized() << "with name" << name << std::endl;
  // std::cout << std::flush;
};

void LibintInterface::init_2body_ints() {

  // compute OBS non-negligible shell-pair list

  obs_shellpair_list = compute_shellpair_list(shells);
  size_t nsp = 0;
  for (auto &sp : obs_shellpair_list) {
    nsp += sp.second.size();
  }

  // std::cout << "# of {all,non-negligible} shell-pairs = {"
  //           << shells.size() * (shells.size() + 1) / 2 << "," << nsp << "}"
  //           << std::endl;
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
  libint2::Shell::do_enforce_unit_normalization(false);

  libint2::svector<double> exponents(nprim);
  libint2::svector<double> coefficients(nprim);

  // exponents.assign(alpha, alpha + nprim);
  // coefficients.assign(coeff, coeff + nprim);

  for(auto a=0; a<nprim; ++a) {
    exponents[a]=alpha[a];
    coefficients[a]=coeff[a];
  }
  
  shells.push_back(
		   libint2::Shell{ //
		     exponents,
		     {
		       {l, false, coefficients}
		     },
		     {{origin[0], origin[1], origin[2]}}
		   }
		   );

  nbasis = 0;
  max_nprim = 0;
  max_l = 0;
  
  for (const auto &shell : shells) {
    nbasis += shell.size();
    max_nprim = std::max(shell.nprim(), max_nprim);
    for (auto c : shell.contr)
      max_l = std::max(c.l, max_l);
  }

  // Renormalize
  const auto &shell = shells.back();

  // std::cout << "libint2::initialized() in add shell" << libint2::initialized() << "with name" << name << std::endl;
  // std::cout << std::flush;
 
  libint2::Engine engine(libint2::Operator::overlap, shell.nprim(), max_l, 0);
  const auto &buf = engine.results();
  engine.compute(shell, shell);
  Eigen::Map<const Matrix> buf_mat(buf[0], shell.size(), shell.size());
  for (int i = 0; i < shell.size(); ++i) {
    norma.push_back(1 / sqrt(buf_mat(i, i)));
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
  const auto &buf = engine.results();

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
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the
      // corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                    // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  // normalize
  for (int i = 0; i < result.cols(); ++i) {
    for (int j = 0; j < result.rows(); ++j) {
      result(i, j) *= norma[i] * norma[j];
    }
  }

  return result;
}

void LibintInterface::compute_2body_disk(const char *filename, const Matrix &D,
                                         const Matrix &Schwartz,
                                         double precision) {
  const auto nshells = shells.size();

  using libint2::nthreads;

  bool do_schwartz_screen = is_electron;
  if (is_electron) {
    do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  }

  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = 0.0;
  if (do_schwartz_screen) {
    engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                std::numeric_limits<double>::epsilon()) /
                       max_nprim4;
  } else {
    engine_precision =
        std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
        max_nprim4;
  }
  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_disk:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    const auto &buf = engines[thread_id].results();
    auto &engine = engines[thread_id];
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(s_size);
    buffer.q.reserve(s_size);
    buffer.r.reserve(s_size);
    buffer.s.reserve(s_size);
    buffer.val.reserve(s_size);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif
            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);


            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                buffer.p[counter] = intIter.i() + 1;
                buffer.q[counter] = intIter.j() + 1;
                buffer.r[counter] = intIter.k() + 1;
                buffer.s[counter] = intIter.l() + 1;
                buffer.val[counter] = buf[0][intIter.index()] *
                                      norma[intIter.i()] * norma[intIter.j()] *
                                      norma[intIter.k()] * norma[intIter.l()];

                // printf("\t(%2d %2d | %2d %2d) = %20.15f\n", p[counter],
                // q[counter], r[counter], s[counter], integral[counter]);

                ++num_ints_computed;
                ++counter;
              }

              if (counter == s_size) {
                write_buffer(buffer, s_size, outfile);
                counter = 0;
              }
            }
          }
        }
      }
    }

    buffer.p[counter] = -1;
    write_buffer(buffer, s_size, outfile);
    outfile.close();

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intra-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  std::cout << " Number of unique integrals for species: " << speciesID << " = "
            << num_ints_computed << std::endl;
}

Matrix LibintInterface::compute_2body_direct(const Matrix &D,
                                             const Matrix &Schwartz,
					     double &factor, 
                                             double precision) {
  const auto n = nbasis;
  const auto nshells = shells.size();

  using libint2::nthreads;

  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  bool do_schwartz_screen = is_electron;
  if (is_electron) {
    do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  }

  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = 0.0;
  if (do_schwartz_screen) {
    engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                std::numeric_limits<double>::epsilon()) /
                       max_nprim4;
  } else {
    engine_precision =
        std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
        max_nprim4;
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    auto &g = G[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf[0][f1234];

                    const auto value_scal_by_deg = value * s1234_deg *
		      norma[bf1] * norma[bf2] *
		      norma[bf3] * norma[bf4];

                    g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                    g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                    g(bf1, bf3) += 0.5 * factor * D(bf2, bf4) * value_scal_by_deg;
                    g(bf2, bf4) += 0.5 * factor * D(bf1, bf3) * value_scal_by_deg;
                    g(bf1, bf4) += 0.5 * factor * D(bf2, bf3) * value_scal_by_deg;
                    g(bf2, bf3) += 0.5 * factor * D(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  }; // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intras-pecies integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // = "
  //           << num_ints_computed << std::endl;
  // symmetrize the result and return
  Matrix GG = 0.25 * (G[0] + G[0].transpose());
  return GG;
}

void LibintInterface::compute_2body_directIT(const Matrix &D, const Matrix &C,
                                             const Matrix &Schwartz,
					     int &p, 
					     double *A,
                                             double precision) {
  const auto n = nbasis;
  const auto nshells = shells.size();

  using libint2::nthreads;

	// allocate a 3D temporary array GG
  double*** GG = new double**[n];
  for ( int i = 0; i < n; i++)
  {
		GG[i] = new double*[n];
  		for ( int j = 0; j < n; j++) 
				GG[i][j] = new double[n];
  }
	// initialize GG
  for ( int i = 0; i < n; i++ )
	{
  	for ( int j = 0; j < n; j++ )
		{
  		for ( int k = 0; k < n; k++ )
			{
				GG[i][j][k] = 0;
			}
		}
	}

  //printf("%e\n", C(0,1));
  //printf("ppp %i\n", p);

  // std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));
  bool do_schwartz_screen = is_electron;
  if (is_electron) {
    do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  }

  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = 0.0;
  if (do_schwartz_screen) {
    engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                std::numeric_limits<double>::epsilon()) /
                       max_nprim4;
  } else {
    engine_precision =
        std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
        max_nprim4;
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  // std::cout << "libint2::initialized() in compute_2body" << libint2::initialized() << "with name" << name << std::endl;
  // std::cout << std::flush;

  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  int m = 0;
  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    // auto &g = G[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            //auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            //auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            //auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            //auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);

            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                auto bf1 = intIter.i() ;
                auto bf2 = intIter.j() ;
                auto bf3 = intIter.k() ;
                auto bf4 = intIter.l() ;
                const auto value_scal_by_deg = buf[0][intIter.index()] *
                                      norma[intIter.i()] * norma[intIter.j()] *
                                      norma[intIter.k()] * norma[intIter.l()];

								#pragma omp critical
								{
          	  	  if (bf3 < bf4)
									{
	          		  	//printf ("A1 1234 %i %i %i %i %f %f %f\n", bf1+1, bf2+1,bf3+1,bf4+1,GG[bf2][bf3][bf4],value_scal_by_deg, C(p, bf1));
	            			GG[bf2][bf3][bf4] = GG[bf2][bf3][bf4] + value_scal_by_deg * C(p, bf1);
	            		} else {
	          			  //printf ("A2 1234 %i %i %i %i %f %f %f\n", bf1+1, bf2+1,bf4+1,bf3+1,GG[bf2][bf4][bf3],value_scal_by_deg, C(p, bf1));
	            	    GG[bf2][bf4][bf3] = GG[bf2][bf4][bf3] + value_scal_by_deg * C(p, bf1);
		    					}

              	  if (bf1 != bf2)
									{
	            	  	if (bf3 < bf4)	
										{
	           					//printf ("B1 2134 %i %i %i %i %f %f %f\n", bf2+1, bf1+1,bf3+1,bf4+1,GG[bf1][bf3][bf4],value_scal_by_deg, C(p, bf2));
	            	      GG[bf1][bf3][bf4] = GG[bf1][bf3][bf4] + value_scal_by_deg * C(p, bf2);
										} else {
	           					//printf ("B2 2134 %i %i %i %i %f %f %f\n", bf2+1, bf1+1,bf4+1,bf3+1,GG[bf1][bf4][bf3],value_scal_by_deg, C(p, bf2));
	            	      GG[bf1][bf4][bf3] = GG[bf1][bf4][bf3] + value_scal_by_deg * C(p, bf2);
		        				}
 
	            		}
              	  if (bf1 != bf3 || bf2 != bf4 )
									{

	            	    if (bf1 < bf2)
										{
	            	      //printf ("C1 3412 %i %i %i %i %f %f %f\n", bf3+1, bf4+1,bf1+1,bf2+1,GG[bf4][bf1][bf2],value_scal_by_deg, C(p, bf3));
	            		  	GG[bf4][bf1][bf2] = GG[bf4][bf1][bf2] + value_scal_by_deg * C(p, bf3);
										} else {
	            	      // printf ("C2 3412 %i %i %i %i %f %f %f\n", bf3+1, bf4+1,bf2+1,bf1+1,GG[bf4][bf2][bf1],value_scal_by_deg, C(p, bf3));
	          		   	  GG[bf4][bf2][bf1] = GG[bf4][bf2][bf1] + value_scal_by_deg * C(p, bf3);
										}

										if (bf3 != bf4)
										{
		          				if (bf1 < bf2)
											{
	            	      	//printf ("D1 4312 %i %i %i %i %f %f %f\n", bf4+1, bf3+1,bf1+1,bf2+1,GG[bf3][bf1][bf2],value_scal_by_deg, C(p, bf4));
		          	  			GG[bf3][bf1][bf2] = GG[bf3][bf1][bf2] + value_scal_by_deg * C(p, bf4);
			  							} else {
	            	      	//printf ("D2 4312 %i %i %i %i %f %f %f\n", bf4+1, bf3+1,bf2+1,bf1+1,GG[bf3][bf2][bf1],value_scal_by_deg, C(p, bf4));
	            		  		GG[bf3][bf2][bf1] = GG[bf3][bf2][bf1] + value_scal_by_deg * C(p, bf4);
			  							}
										}
	            		}

								} //end critical
	         		} //end if >10
	       		} //end for

          }
        }
      }
    }

  }; // end of lambda

  libint2::parallel_do(lambda);
 
	// symmetrize
  for ( int i = 0; i < n; i++ )
  {
  	for ( int j = 0; j < n; j++ )
    {
  		for ( int k = j; k < n; k++ )
      {
      	GG[i][k][j] = GG[i][j][k] ;
      }
    }
  }
  
  m = 0;

  for ( int i = 0; i < n; i++ )
  {
  	for ( int j = 0; j < n; j++ )
    {
  		for ( int k = 0; k < n; k++ )
      {
      	A[m] = GG[i][j][k] ;
        //printf ("%i %i %i %2.8f\n", i+1,j+1,k+1,A[m] );
        m = m + 1 ;
      }
    }
  }

  // // accumulate contributions from all threads
  // for (size_t i = 1; i != nthreads; ++i) {
  //   G[0] += G[i];
  // }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intras-pecies integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // do I really have to free this memory?
  for ( int i = 0; i < n; i++)
    {
      for ( int j = 0; j < n; j++) 
	delete(GG[i][j]);
      delete(GG[i]);
    }
  delete(GG);
  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // = "
  //           << num_ints_computed << std::endl;
  // symmetrize the result and return
  //Matrix GG = 0.25 * (G[0] + G[0].transpose());
  //return GG;
}


void LibintInterface::compute_coupling_disk(LibintInterface &other,
                                            const char *filename,
                                            double precision) {
  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision = 1.0e-40;
      //std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      //max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};
  int screened;

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  screened = 0;
  auto lambda = [&](int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(s_size);
    buffer.q.reserve(s_size);
    buffer.r.reserve(s_size);
    buffer.s.reserve(s_size);
    buffer.val.reserve(s_size);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr){
              screened+=1;
              continue; // all screened out
            }

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      buffer.p[counter] = bf1 + 1;
                      buffer.q[counter] = bf2 + 1;
                      buffer.r[counter] = obf1 + 1;
                      buffer.s[counter] = obf2 + 1;
                      buffer.val[counter] =
                          buf[0][f1212] * get_norma(bf1) * get_norma(bf2) *
                          other.get_norma(obf1) * other.get_norma(obf2);

                      // printf("(%d %d | %d %d) %lf \n", buffer.p[counter],
                      // buffer.q[counter], buffer.r[counter],
                      // buffer.s[counter], buffer.val[counter] *
                      //                  get_norma(bf1) * get_norma(bf2) *
                      //                  other.get_norma(obf1) *
                      //                  other.get_norma(obf2));

                      ++num_ints_computed;
                      ++counter;

                      if (counter == s_size) {
                        write_buffer(buffer, s_size, outfile);
                        counter = 0;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

    buffer.p[counter] = -1;
    write_buffer(buffer, s_size, outfile);
    outfile.close();

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  std::cout << " Number of unique integrals for species: " << speciesID << " / "
            << other.get_speciesID() << " = " << num_ints_computed << " Screened: "<< screened<<std::endl;
}

Matrix LibintInterface::compute_coupling_direct(LibintInterface &other,
                                                const Matrix &D,
                                                const bool permuted,
                                                double precision) {
  const auto n = permuted ? other.get_nbasis() : get_nbasis();

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  std::vector<Matrix> B(nthreads, Matrix::Zero(n, n));

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engine_precision = 1.0e-40;
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

// printf("(%ld %ld | %ld %ld) \n", s1, s2, os1, os2);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      // printf("(%ld %ld | %ld %ld) %lf \n", bf1, bf2, obf1,
                      // obf2, buf[0][f1212] *
                      //                  get_norma(bf1) * get_norma(bf2) *
                      //                  other.get_norma(obf1) *
                      //                  other.get_norma(obf2));

                      if (permuted) {
                        auto coulomb = D(bf1, bf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        b(obf2, obf1) += coulomb;
                        if (bf1 != bf2)
                          b(obf2, obf1) += coulomb;
                      } else {
                        auto coulomb = D(obf1, obf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        // printf("%lf\n", coulomb);
                        b(bf2, bf1) += coulomb;
                        if (obf1 != obf2)
                          b(bf2, bf1) += coulomb;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    B[0] += B[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // / "
  //           << other.get_speciesID() << " = " << num_ints_computed <<
  //           std::endl;

  return B[0];
}

void LibintInterface::compute_coupling_directIT(LibintInterface &other,
                                                const Matrix &D, const Matrix &C, 
						int &p, double *A, 
                                                const bool permuted,
                                                double precision) {
  const auto n = permuted ? other.get_nbasis() : get_nbasis();
  const auto on = permuted ? get_nbasis() : other.get_nbasis()  ;

	//if ( permuted ) {
	//	printf ("permuted\n");
	//}
	//printf ("second n: %i on: %i \n", n, on);
  //printf ("p index: %i \n", p);

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

	// allocate a 3D temporary array GG

  int nn1 = 0;
  int nn2 = 0;
  int nn3 = 0;

  if ( permuted ) {
    nn1 = n;
    nn2 = on;
    nn3 = on;
  } else {
    nn1 = n;
    nn2 = on;
    nn3 = on;
  }

	//printf ("%i %i %i\n",nn1,nn2,nn3);
  double*** GG = new double**[nn1];
  for ( int i = 0; i < nn1; i++)
  {
		GG[i] = new double*[nn2];
  		for ( int j = 0; j < nn2; j++) 
				GG[i][j] = new double[nn3];
  }
	// initialize GG
  for ( int i = 0; i < nn1; i++ )
	{
  	for ( int j = 0; j < nn2; j++ )
		{
  		for ( int k = 0; k < nn3; k++ )
			{
				GG[i][j][k] = 0;
			}
		}
	}

  // std::vector<Matrix> B(nthreads, Matrix::Zero(n, n));

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engine_precision = 1.0e-40;
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    // auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

// printf("(%ld %ld | %ld %ld) \n", s1, s2, os1, os2);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      // printf("(%ld %ld | %ld %ld) %lf \n", bf1, bf2, obf1,
                      // obf2, buf[0][f1212] *
                      //                  get_norma(bf1) * get_norma(bf2) *
                      //                  other.get_norma(obf1) *
                      //                  other.get_norma(obf2));

                      ///if (permuted) {
                      ///  auto coulomb = D(bf1, bf2) * buf[0][f1212] *
                      ///                 get_norma(bf1) * get_norma(bf2) *
                      ///                 other.get_norma(obf1) *
                      ///                 other.get_norma(obf2);

                      ///  b(obf2, obf1) += coulomb;
                      ///  if (bf1 != bf2)
                      ///    b(obf2, obf1) += coulomb;
                      ///} else {
                      ///  auto coulomb = D(obf1, obf2) * buf[0][f1212] *
                      ///                 get_norma(bf1) * get_norma(bf2) *
                      ///                 other.get_norma(obf1) *
                      ///                 other.get_norma(obf2);

                      ///  // printf("%lf\n", coulomb);
                      ///  b(bf2, bf1) += coulomb;
                      ///  if (obf1 != obf2)
                      ///    b(bf2, bf1) += coulomb;
											///}

                      auto coulomb = buf[0][f1212] *
                                     get_norma(bf1) * get_norma(bf2) *
                                     other.get_norma(obf1) *
                                     other.get_norma(obf2);

											//printf ("%i %i %i %i %f\n", bf1, bf2, obf1, obf2, coulomb);

								      #pragma omp critical
								      {
                      	if (permuted) {
          	  	          if (bf1 < bf2)
								      	  {
	            	      	  	GG[obf2][bf1][bf2] = GG[obf2][bf1][bf2] + coulomb * C(p, obf1);
	            	      	  } else {
	          		      	    //printf ("A2 1234 %i %i %i %i %f %f %f\n", obf1+1, obf2+1,bf2+1,bf1+1,GG[obf2][bf2][bf1],coulomb, C(p, obf1));
	            	            GG[obf2][bf2][bf1] = GG[obf2][bf2][bf1] + coulomb * C(p, obf1);
		    				      	  }

              	          if (obf1 != obf2)
								      	  {
	            	          	if (bf1 < bf2)	
								      	  	{
	           		      	  		//printf ("B1 2134 %i %i %i %i %f %f %f\n", obf2+1, obf1+1,bf1+1,bf2+1,GG[obf1][bf1][bf2],coulomb, C(p, obf2));
	            	              GG[obf1][bf1][bf2] = GG[obf1][bf1][bf2] + coulomb * C(p, obf2);
								      	  	} else {
	           		      	  		//printf ("B2 2134 %i %i %i %i %f %f %f\n", obf2+1, obf1+1,bf2+1,bf1+1,GG[obf1][bf2][bf1],coulomb, C(p, obf2));
	            	              GG[obf1][bf2][bf1] = GG[obf1][bf2][bf1] + coulomb * C(p, obf2);
		        		      	  	}
 
	            	      	  }

								      	} else {

          	  	          if (obf1 < obf2)
								      	  {
	          		          	//printf ("A1 1234 %i %i %i %i %f %f %f\n", bf1+1, bf2+1,obf1+1,obf2+1,GG[bf2][obf1][obf2],coulomb, C(p, bf1));
	            	      	  	GG[bf2][obf1][obf2] = GG[bf2][obf1][obf2] + coulomb * C(p, bf1);
	            	      	  } else {
	          		      	    //printf ("A2 1234 %i %i %i %i %f %f %f\n", bf1+1, bf2+1,obf2+1,obf1+1,GG[bf2][obf2][obf1],coulomb, C(p, bf1));
	            	            GG[bf2][obf2][obf1] = GG[bf2][obf2][obf1] + coulomb * C(p, bf1);
		    				      	  }

              	          if (bf1 != bf2)
								      	  {
	            	          	if (obf1 < obf2)	
								      	  	{
	           		      	  		//printf ("B1 2134 %i %i %i %i %f %f %f\n", bf2+1, bf1+1,obf1+1,obf2+1,GG[bf1][obf1][obf2],coulomb, C(p, bf2));
	            	              GG[bf1][obf1][obf2] = GG[bf1][obf1][obf2] + coulomb * C(p, bf2);
								      	  	} else {
	           		      	  		//printf ("B2 2134 %i %i %i %i %f %f %f\n", bf2+1, bf1+1,obf2+1,obf1+1,GG[bf1][obf2][obf1],coulomb, C(p, bf2));
	            	              GG[bf1][obf2][obf1] = GG[bf1][obf2][obf1] + coulomb * C(p, bf2);
		        		      	  	}
 
	            	      	  }
								      
								      	}

								      } //end critical


                    } // end if
                  }
                }
              }
            } // end cartesian loop

          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

	// symmetrize
  for ( int i = 0; i < nn1; i++ )
  {
  	for ( int j = 0; j < nn2; j++ )
    {
  		for ( int k = j; k < nn3; k++ )
      {
      	GG[i][k][j] = GG[i][j][k] ;
      }
    }
  }
  
  int m = 0;
  m = 0;

  for ( int i = 0; i < nn1; i++ )
  {
  	for ( int j = 0; j < nn2; j++ )
    {
  		for ( int k = 0; k < nn3; k++ )
      {
      	A[m] = GG[i][j][k] ;
        //printf ("%i %i %i %2.8f\n", i+1,j+1,k+1,A[m] );
        m = m + 1 ;
      }
    }
  }

  // accumulate contributions from all threads
  // for (size_t i = 1; i != nthreads; ++i) {
  //   B[0] += B[i];
  // }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // do I really have to free this memory?
  for ( int i = 0; i < nn1; i++)
    {
      for ( int j = 0; j < nn2; j++) 
	delete(GG[i][j]);
      delete(GG[i]);
    }
  delete(GG);

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // / "
  //           << other.get_speciesID() << " = " << num_ints_computed <<
  //           std::endl;

  //return B[0];
}



Matrix LibintInterface::compute_alphabeta_direct(LibintInterface &other,
                                                const Matrix &D,
                                                const Matrix &oD,
                                                const bool permuted,
                                                double precision) {
  const auto n = permuted ? other.get_nbasis() : get_nbasis();

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  //std::vector<double> B(nthreads);
  std::vector<Matrix> B(nthreads, Matrix::Zero(1, 1));

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engine_precision = 1.0e-40;
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

// printf("(%ld %ld | %ld %ld) \n", s1, s2, os1, os2);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      // printf("(%ld %ld | %ld %ld) %lf \n", bf1, bf2, obf1,
                      // obf2, buf[0][f1212] *
                      //                  get_norma(bf1) * get_norma(bf2) *
                      //                  other.get_norma(obf1) *
                      //                  other.get_norma(obf2));

                      if (permuted) {
                        auto coulomb = D(bf1, bf2) * buf[0][f1212] * oD(obf1, obf2) *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                         //printf("%lf\n", coulomb);
                        b(0,0) += coulomb;
                        if (bf1 != bf2)
                          b(0,0) += coulomb;
                        if (obf1 != obf2)
                          b(0,0) += coulomb;
                        if (bf1 != bf2 && obf1 != obf2)
                          b(0,0) += coulomb;

                      } else {
                        auto coulomb = buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        //printf("%lf\n", coulomb );

                        coulomb = coulomb * oD(obf1, obf2) * D(bf1,bf2) ;

                        b(0,0) += coulomb;
                        if (obf1 != obf2)
                          b(0,0) += coulomb;
                        if (bf1 != bf2)
                          b(0,0) += coulomb;
                        if (bf1 != bf2 && obf1 != obf2)
                          b(0,0) += coulomb;

                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    B[0] += B[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // / "
  //           << other.get_speciesID() << " = " << num_ints_computed <<
  //           std::endl;

  return B[0];
}



void LibintInterface::compute_g12_disk(const char *filename,
                                       const double *coefficients,
                                       const double *exponents,
                                       const int pot_size) {

  const auto nshells = shells.size();

  using libint2::nthreads;

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  // engines[0].set_precision(engine_precision); // shellset-dependentprecision
  // control will likely break
  // positive definiteness
  // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_disk:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    const auto &buf = engines[thread_id].results();
    auto &engine = engines[thread_id];
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(s_size);
    buffer.q.reserve(s_size);
    buffer.r.reserve(s_size);
    buffer.s.reserve(s_size);
    buffer.val.reserve(s_size);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif
            if (buf[0] == nullptr)
              continue; // all screened out

            // printf("(%2d %2d| %2d %2d)\n", s1, s2, s3, s4);

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);

            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                buffer.p[counter] = intIter.i() + 1;
                buffer.q[counter] = intIter.j() + 1;
                buffer.r[counter] = intIter.k() + 1;
                buffer.s[counter] = intIter.l() + 1;
                buffer.val[counter] = buf[0][intIter.index()] *
                                      norma[intIter.i()] * norma[intIter.j()] *
                                      norma[intIter.k()] * norma[intIter.l()];

                //     printf("counter: %d buffer_size: %d Integral: %lf\n",
                //     counter, s_size, buf[0][intIter.index()]);
                //     printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
                //     buffer.p[counter],
                //            buffer.q[counter], buffer.r[counter],
                //            buffer.s[counter], buffer.val[counter]);

                ++num_ints_computed;
                ++counter;
              }

              if (counter == s_size) {
                write_buffer(buffer, s_size, outfile);
                counter = 0;
              }
            }
          }
        }
      }
    }

    buffer.p[counter] = -1;
    write_buffer(buffer, s_size, outfile);
    outfile.close();

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intra-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  std::cout << " Number of unique integrals for species: " << speciesID << " = "
            << num_ints_computed << std::endl;
}

void LibintInterface::compute_g12inter_disk(LibintInterface &other,
                                            const char *filename,
                                            const double *coefficients,
                                            const double *exponents,
                                            const int pot_size) {

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto aux = std::max(max_nprim, other.max_nprim);
  auto max_nprim = aux;
  auto aux2 = std::max(max_l, other.max_l);
  auto max_l = aux2;

  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_disk:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    const auto &buf = engines[thread_id].results();
    auto &engine = engines[thread_id];
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(s_size);
    buffer.q.reserve(s_size);
    buffer.r.reserve(s_size);
    buffer.s.reserve(s_size);
    buffer.val.reserve(s_size);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], oshells[os1], oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      buffer.p[counter] = bf1 + 1;
                      buffer.q[counter] = bf2 + 1;
                      buffer.r[counter] = obf1 + 1;
                      buffer.s[counter] = obf2 + 1;
                      buffer.val[counter] =
                          buf[0][f1212] * get_norma(bf1) * get_norma(bf2) *
                          other.get_norma(obf1) * other.get_norma(obf2);

                      ++num_ints_computed;
                      ++counter;

                      if (counter == s_size) {
                        write_buffer(buffer, s_size, outfile);
                        counter = 0;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

    buffer.p[counter] = -1;
    write_buffer(buffer, s_size, outfile);
    outfile.close();

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  std::cout << " Number of unique integrals for species: " << speciesID << " / "
            << other.get_speciesID() << " = " << num_ints_computed << std::endl;
}
Matrix LibintInterface::compute_g12_direct(const Matrix &D,
					   double &factor, 
					   const double *coefficients,
					   const double *exponents,
					   const int pot_size) {

  const auto n = nbasis;
  const auto nshells = shells.size();

  using libint2::nthreads;

  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  // engines[0].set_precision(engine_precision); // shellset-dependentprecision
  // control will likely break
  // positive definiteness
  // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    auto &g = G[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif
            if (buf[0] == nullptr)
              continue; // all screened out

            // printf("(%2d %2d| %2d %2d)\n", s1, s2, s3, s4);

            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf[0][f1234];

                    const auto value_scal_by_deg = value * s1234_deg *
		      norma[bf1] * norma[bf2] *
		      norma[bf3] * norma[bf4];

                    g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                    g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                    g(bf1, bf3) += 0.5 * factor * D(bf2, bf4) * value_scal_by_deg;
                    g(bf2, bf4) += 0.5 * factor * D(bf1, bf3) * value_scal_by_deg;
                    g(bf1, bf4) += 0.5 * factor * D(bf2, bf3) * value_scal_by_deg;
                    g(bf2, bf3) += 0.5 * factor * D(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }
    
  }; // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }
  
#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intra-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // = "
  //           << num_ints_computed << std::endl;
  // symmetrize the result and return
  Matrix GG = 0.25 * (G[0] + G[0].transpose());
  return GG;
}

Matrix LibintInterface::compute_g12inter_direct(LibintInterface &other,
						const Matrix &D,
                                                const bool permuted,
						const double *coefficients,
						const double *exponents,
						const int pot_size) {

  const auto n = permuted ? other.get_nbasis() : get_nbasis();

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  std::vector<Matrix> B(nthreads, Matrix::Zero(n, n));

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto aux = std::max(max_nprim, other.max_nprim);
  auto max_nprim = aux;
  auto aux2 = std::max(max_l, other.max_l);
  auto max_l = aux2;

  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], oshells[os1], oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;
                      if (permuted) {
                        auto coulomb = D(bf1, bf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        b(obf2, obf1) += coulomb;
                        if (bf1 != bf2)
                          b(obf2, obf1) += coulomb;
                      } else {
                        auto coulomb = D(obf1, obf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        // printf("%lf\n", coulomb);
                        b(bf2, bf1) += coulomb;
                        if (obf1 != obf2)
                          b(bf2, bf1) += coulomb;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

  for (size_t i = 1; i != nthreads; ++i) {
    B[0] += B[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  return B[0];
}

void LibintInterface::compute_2body_directAll(const Matrix &D, 
					      const Matrix &Schwartz,
					      double *A,
					      double precision) {
  const auto n = nbasis;
  const auto nshells = shells.size();

  using libint2::nthreads;

  bool do_schwartz_screen = is_electron;
  if (is_electron) {
    do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  }

  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = 0.0;
  if (do_schwartz_screen) {
    engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                std::numeric_limits<double>::epsilon()) /
                       max_nprim4;
  } else {
    engine_precision =
        std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
        max_nprim4;
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  // std::cout << "libint2::initialized() in compute_2body" << libint2::initialized() << "with name" << name << std::endl;
  // std::cout << std::flush;

  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  int m = 0;
  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    // auto &g = G[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            //auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            //auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            //auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            //auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);

            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                auto bf1 = intIter.i() ;
                auto bf2 = intIter.j() ;
                auto bf3 = intIter.k() ;
                auto bf4 = intIter.l() ;
                const auto value_scal_by_norm = buf[0][intIter.index()] *
		  norma[intIter.i()] * norma[intIter.j()] *
		  norma[intIter.k()] * norma[intIter.l()];

		A[bf1*n*n*n+bf2*n*n+bf3*n+bf4] = value_scal_by_norm;
		A[bf2*n*n*n+bf1*n*n+bf3*n+bf4] = value_scal_by_norm;
		A[bf1*n*n*n+bf2*n*n+bf4*n+bf3] = value_scal_by_norm;
		A[bf2*n*n*n+bf1*n*n+bf4*n+bf3] = value_scal_by_norm;
		A[bf3*n*n*n+bf4*n*n+bf1*n+bf2] = value_scal_by_norm;
		A[bf3*n*n*n+bf4*n*n+bf2*n+bf1] = value_scal_by_norm;
		A[bf4*n*n*n+bf3*n*n+bf1*n+bf2] = value_scal_by_norm;
		A[bf4*n*n*n+bf3*n*n+bf2*n+bf1] = value_scal_by_norm;
		  
	      } //end if >10
	    } //end for
	    
          }
        }
      }
    }

  }; // end of lambda

  libint2::parallel_do(lambda);
 
#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intras-pecies integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

}

void LibintInterface::compute_coupling_directAll(LibintInterface &other,
                                                const Matrix &D, double *A, 
                                                const bool permuted,
                                                double precision) {
  const auto n = permuted ? other.get_nbasis() : get_nbasis();
  const auto on = permuted ? get_nbasis() : other.get_nbasis()  ;
  
	//if ( permuted ) {
	//	printf ("permuted\n");
	//}
	//printf ("second n: %i on: %i \n", n, on);
  //printf ("p index: %i \n", p);

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engine_precision = 1.0e-40;
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    // auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

// printf("(%ld %ld | %ld %ld) \n", s1, s2, os1, os2);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      auto coulomb = buf[0][f1212] *
			get_norma(bf1) * get_norma(bf2) *
			other.get_norma(obf1) *
			other.get_norma(obf2);

		      // printf ("%i %i %i %i %f\n", obf2, obf1, bf1, bf2 , coulomb);
		      // fflush(stdout);
		      A[bf1*n*on*on+bf2*on*on+obf1*on+obf2]=coulomb;
		      A[bf2*n*on*on+bf1*on*on+obf1*on+obf2]=coulomb;
		      A[bf1*n*on*on+bf2*on*on+obf2*on+obf1]=coulomb;
		      A[bf2*n*on*on+bf1*on*on+obf2*on+obf1]=coulomb;

                    } // end if
                  }
                }
              }
            } // end cartesian loop

          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

}

void LibintInterface::compute_g12_directAll(const Matrix &D,
					    double *A,
					    const double *coefficients,
					    const double *exponents,
					    const int pot_size) {

  const auto n = nbasis;
  const auto nshells = shells.size();

  using libint2::nthreads;

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  // engines[0].set_precision(engine_precision); // shellset-dependentprecision
  // control will likely break
  // positive definiteness
  // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            // auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            // auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            // auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            // auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif
            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);
	    
            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                auto bf1 = intIter.i() ;
                auto bf2 = intIter.j() ;
                auto bf3 = intIter.k() ;
                auto bf4 = intIter.l() ;
                const auto value_scal_by_norm = buf[0][intIter.index()] *
		  norma[intIter.i()] * norma[intIter.j()] *
		  norma[intIter.k()] * norma[intIter.l()];

		A[bf1*n*n*n+bf2*n*n+bf3*n+bf4] = value_scal_by_norm;
		A[bf2*n*n*n+bf1*n*n+bf3*n+bf4] = value_scal_by_norm;
		A[bf1*n*n*n+bf2*n*n+bf4*n+bf3] = value_scal_by_norm;
		A[bf2*n*n*n+bf1*n*n+bf4*n+bf3] = value_scal_by_norm;
		A[bf3*n*n*n+bf4*n*n+bf1*n+bf2] = value_scal_by_norm;
		A[bf3*n*n*n+bf4*n*n+bf2*n+bf1] = value_scal_by_norm;
		A[bf4*n*n*n+bf3*n*n+bf1*n+bf2] = value_scal_by_norm;
		A[bf4*n*n*n+bf3*n*n+bf2*n+bf1] = value_scal_by_norm;
		  
	      } //end if >10
	    } //end for
            // printf("(%2d %2d| %2d %2d)\n", s1, s2, s3, s4);

          }
        }
      }
    }
    
  }; // end of lambda

  libint2::parallel_do(lambda);
  
#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intra-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << speciesID << "
  // = "
  //           << num_ints_computed << std::endl;
  // symmetrize the result and return
}

void LibintInterface::compute_g12inter_directAll(LibintInterface &other,
						 const Matrix &D,
						 double *A,
						 const bool permuted,
						 const double *coefficients,
						 const double *exponents,
						 const int pot_size) {

  const auto n = permuted ? other.get_nbasis() : get_nbasis();
  const auto on = permuted ? get_nbasis() : other.get_nbasis()  ;

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using libint2::nthreads;

  auto fock_precision = std::numeric_limits<double>::epsilon();

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto aux = std::max(max_nprim, other.max_nprim);
  auto max_nprim = aux;
  auto aux2 = std::max(max_l, other.max_l);
  auto max_l = aux2;

  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // Setting CGTG params
  std::vector<std::pair<double, double>> cgtg_params;
  cgtg_params.reserve(pot_size);

  for (int i = 0; i < pot_size; ++i) {
    cgtg_params.push_back({exponents[i], coefficients[i]});
  }

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;

  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::cgtg, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe
  engines[0].set_params(cgtg_params);

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  using libint2::QuartetBuffer;
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (auto s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (auto os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (auto os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute(shells[s1], shells[s2], oshells[os1], oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (auto f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (auto of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (auto of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      auto coulomb = buf[0][f1212] *
			get_norma(bf1) * get_norma(bf2) *
			other.get_norma(obf1) *
			other.get_norma(obf2);

		      // printf ("%i %i %i %i %f\n", obf2, obf1, bf1, bf2 , coulomb);
		      // fflush(stdout);
		      A[bf1*n*on*on+bf2*on*on+obf1*on+obf2]=coulomb;
		      A[bf2*n*on*on+bf1*on*on+obf1*on+obf2]=coulomb;
		      A[bf1*n*on*on+bf2*on*on+obf2*on+obf1]=coulomb;
		      A[bf2*n*on*on+bf1*on*on+obf2*on+obf1]=coulomb;
		      
                    } // end if
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

  }; // end of lambda

  libint2::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species G12 integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

}

/*
Auxiliary functions
*/
shellpair_list_t compute_shellpair_list(const std::vector<libint2::Shell> &bs1,
                                        const std::vector<libint2::Shell> &_bs2,
                                        const double threshold) {
  const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  using libint2::nthreads;

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  std::vector<Engine> engines;

  engines.reserve(nthreads);

  // std::cout << "libint2::initialized() in add compute_shellpair" << libint2::initialized() << "in thread" << omp_get_thread_num()+1 << std::endl;
  // std::cout << std::flush;

  engines.emplace_back(libint2::Operator::overlap,
                       std::max(max_nprim(bs1), max_nprim(bs2)),
                       std::max(max_l(bs1), max_l(bs2)), 0);

  for (size_t i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  // std::cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&](int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();

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
          engine.compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
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
  // std::cout << "done (" << timer.read(0) << " s)" << std::endl;

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
  auto epsilon = 0.;
  engines[0] = Engine(
      Kernel, std::max(max_nprim(bs1), max_nprim(bs2)),
      std::max(max_l(bs1), max_l(bs2)), 0, epsilon, params);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  // std::cout << "computing Schwartz bound prerequisites (kernel=" <<
  // (int)Kernel
  //           << ") ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  auto lambda = [&](int thread_id) {
    const auto &buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      auto n1 = bs1[s1].size(); // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id)
          continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1 * n2;

        engines[thread_id].compute(bs1[s1], bs2[s2], bs1[s1], bs2[s2]);

        // the diagonal elements are the Schwartz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
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

__inline__ void write_buffer(const libint2::QuartetBuffer &buffer,
                             const int &s_size, std::ofstream &outfile) {
  outfile.write(reinterpret_cast<const char *>(&buffer.p[0]),
                s_size * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.q[0]),
                s_size * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.r[0]),
                s_size * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.s[0]),
                s_size * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.val[0]),
                s_size * sizeof(double));
}

size_t nbasis(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

/*
Fortran interface
*/

LibintInterface *LibintInterface_new(const int stack_size, const int id,
                                     const bool el, const bool parallel) {
  return new LibintInterface(stack_size, id, el, parallel);
  // printf("%s\n", "LibintInterface_new");
}

void LibintInterface_del(LibintInterface *lint) {
  // printf("%s\n", "LibintInterface_del");
  // lint->~LibintInterface();
  // std::cout << "deleting object with name" << lint->get_name() << std::endl;
  delete(lint);
}

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center) {
  // printf("%s\n", "LibintInterface_add_particle");
  lint->add_particle(z, center);
}

void LibintInterface_add_shell(LibintInterface *lint, double *alpha,
                               double *coeff, double *origin, int l,
                               int nprim) {
  // printf("%s\n", "LibintInterface_add_shell");
  lint->add_shell(alpha, coeff, origin, l, nprim);
}

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result) {
  // Default case
  // printf("%s\n", "LibintInterface_compute_1body_ints");
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
  // printf("%s\n", "LibintInterface_init_2body_ints");
  lint->init_2body_ints();
}

void LibintInterface_compute_2body_direct(LibintInterface *lint, double *dens,
                                          double *result, double &factor) {
  // printf("%s\n", "LibintInterface_compute_2body_direct");
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto K = compute_schwartz_ints<>(lint->get_shells());
  auto G = lint->compute_2body_direct(D, K, factor);

  for (int i = 0; i < G.size(); ++i) {
    result[i] = G.array()(i);
  }
}

void LibintInterface_compute_2body_directIT(LibintInterface *lint, double *dens, double *coeff,
                                          double *result, int &p) {
  // printf("%s\n", "LibintInterface_compute_2body_direct");
  auto n = lint->get_nbasis();
  Matrix D(n, n);
  Matrix C(n, n);

  for (int i = 0; i < C.size(); ++i) {
    C.array()(i) = coeff[i];
  }
  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto K = compute_schwartz_ints<>(lint->get_shells());
  //auto G = lint->compute_2body_directIT(C, K, p);
  lint->compute_2body_directIT(D, C, K, p, result);

  //for (int i = 0; i < G.size(); ++i) {
  //  result[i] = G.array()(i);
  //}
}

void LibintInterface_compute_2body_disk(LibintInterface *lint,
                                        const char *filename, double *dens) {
  // printf("%s\n", "LibintInterface_compute_2body_disk");
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto K = compute_schwartz_ints<>(lint->get_shells());

  lint->compute_2body_disk(filename, D, K);
}

void LibintInterface_compute_coupling_direct(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *result) {
  // printf("%s\n", "LibintInterface_compute_coupling_direct");
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = olint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto B = permuted ? olint->compute_coupling_direct(*lint, D, permuted)
                    : lint->compute_coupling_direct(*olint, D, permuted);

  for (int i = 0; i < B.size(); ++i) {
    result[i] = B.array()(i);
  }
}

void LibintInterface_compute_coupling_directIT(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *coeff, double *results, int &p) {
  // printf("%s\n", "LibintInterface_compute_coupling_direct");
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = lint->get_nbasis();
  auto on = olint->get_nbasis();

	//printf ("s: %i os: %i \n", sid, osid);
	//printf ("first n: %i on: %i \n", n, on);
  Matrix D(n, n);
  Matrix C(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }
  for (int i = 0; i < C.size(); ++i) {
    C.array()(i) = coeff[i];
  }
 
  permuted ? olint->compute_coupling_directIT(*lint, D, C, p, results, permuted)
                    : lint->compute_coupling_directIT(*olint, D, C, p, results, permuted);

  //for (int i = 0; i < B.size(); ++i) {
  //  result[i] = B.array()(i);
  //}
}



void LibintInterface_compute_alphabeta_direct(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *otherdens, double *result) {
  // printf("%s\n", "LibintInterface_compute_alphabeta_direct");
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = lint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto on = olint->get_nbasis();

  Matrix oD(on, on);

  for (int i = 0; i < oD.size(); ++i) {
    oD.array()(i) = otherdens[i];
  }



  auto B = permuted ? olint->compute_alphabeta_direct(*lint, D, oD, permuted)
                    : lint->compute_alphabeta_direct(*olint, oD, D, permuted);

  for (int i = 0; i < B.size(); ++i) {
    result[i] = B.array()(i);
  }
  //for (int i = 0; i < B.size(); ++i) {
  //  result[i] = B.array()(i);
  //}
}

void LibintInterface_compute_coupling_disk(LibintInterface *lint,
                                           LibintInterface *olint,
                                           const char *filename) {
  // printf("%s\n", "LibintInterface_compute_coupling_disk");
  lint->compute_coupling_disk(*olint, filename);
}

void libintinterface_compute_g12_disk(LibintInterface *lint,
                                      const char *filename,
                                      const double *coefficients,
                                      const double *exponents,
                                      const int pot_size) {

  lint->compute_g12_disk(filename, coefficients, exponents, pot_size);
}

void libintinterface_compute_g12inter_disk(
    LibintInterface *lint, LibintInterface *olint, const char *filename,
    const double *coefficients, const double *exponents, const int pot_size) {

  lint->compute_g12inter_disk(*olint, filename, coefficients, exponents,
                              pot_size);
}

void libintinterface_compute_g12_direct(LibintInterface *lint,
					double *dens,
					double *result,
					double &factor,
					const double *coefficients,
					const double *exponents,
					const int pot_size) {
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto G = lint->compute_g12_direct(D, factor, coefficients, exponents, pot_size); 
  
  for (int i = 0; i < G.size(); ++i) {
    result[i] = G.array()(i);
  }
}

void libintinterface_compute_g12inter_direct(LibintInterface *lint,
					     LibintInterface *olint,
					     double *dens,
					     double *result,
					     const double *coefficients,
					     const double *exponents,
					     const int pot_size) {
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = olint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto B = permuted ? olint->compute_g12inter_direct(*lint, D, permuted, coefficients, exponents, pot_size)
                    : lint->compute_g12inter_direct(*olint, D, permuted, coefficients, exponents, pot_size);
  
  for (int i = 0; i < B.size(); ++i) {
    result[i] = B.array()(i);
  }
}

void LibintInterface_compute_2body_directAll(LibintInterface *lint, double *dens,
                                          double *results) {
  // printf("%s\n", "LibintInterface_compute_2body_direct");
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto K = compute_schwartz_ints<>(lint->get_shells());

  lint->compute_2body_directAll(D, K, results);

}

void LibintInterface_compute_coupling_directAll(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *results) {
  // printf("%s\n", "LibintInterface_compute_coupling_direct");
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = lint->get_nbasis();
  auto on = olint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }
  
  permuted ? olint->compute_coupling_directAll(*lint, D, results, permuted)
    : lint->compute_coupling_directAll(*olint, D, results, permuted);

}

void libintinterface_compute_g12_directAll(LibintInterface *lint,
					double *dens,
					double *result,
					const double *coefficients,
					const double *exponents,
					const int pot_size) {
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  lint->compute_g12_directAll(D, result, coefficients, exponents, pot_size); 
  
}

void libintinterface_compute_g12inter_directAll(LibintInterface *lint,
					     LibintInterface *olint,
					     double *dens,
					     double *results,
					     const double *coefficients,
					     const double *exponents,
					     const int pot_size) {
  auto sid = lint->get_speciesID();
  auto osid = olint->get_speciesID();

  auto permuted = sid > osid;

  auto n = olint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  permuted ? olint->compute_g12inter_directAll(*lint, D, results, permuted, coefficients, exponents, pot_size)
    : lint->compute_g12inter_directAll(*olint, D, results, permuted, coefficients, exponents, pot_size);
  
}

void libintinterface_buildg12_(int *nd, int *nc, int *nb, int *na, int *max_am,
                               int *size, lowdin_t *data, double *output) {
  printf("TEST CALCULATING:\n");
  Libint_t inteval;

  LIBINT2_PREFIXED_NAME(libint2_init_r12kg12)(&inteval, *max_am, 0);

  LibintInterface_setLibint(&inteval, data);

  inteval.contrdepth = 1;
  // Calculate integrals
  LIBINT2_PREFIXED_NAME(libint2_build_r12kg12)[*na][*nb][*nc][*nd](&inteval);

  for (int i = 0; i < *size; ++i) {
    output[i] = inteval.targets[0][i];
    // printf("output %d %d %d %d %f \n", *na, *nb, *nc, *nd,
    //        inteval.targets[0][i]);
  }

  //   this releases all memory that was allocated for this object
  // LIBINT2_PREFIXED_NAME( libint2_cleanup_r12kg12)(&inteval);
}

/*
Set values for Libint_t
*/
void LibintInterface_setLibint(Libint_t *erieval, lowdin_t *data) {
/** Appear in standard OS RR for ERI and almost all other recurrence relations
 */
#if LIBINT2_DEFINED(eri, WP_x)
  erieval->WP_x[0] = data->WP_x;
#endif
#if LIBINT2_DEFINED(eri, WP_y)
  erieval->WP_y[0] = data->WP_y;
#endif
#if LIBINT2_DEFINED(eri, WP_z)
  erieval->WP_z[0] = data->WP_z;
#endif

#if LIBINT2_DEFINED(eri, WQ_x)
  erieval->WQ_x[0] = data->WQ_x;
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
  erieval->WQ_y[0] = data->WQ_y;
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
  erieval->WQ_z[0] = data->WQ_z;
#endif

#if LIBINT2_DEFINED(eri, PA_x)
  erieval->PA_x[0] = data->PA_x;
#endif
#if LIBINT2_DEFINED(eri, PA_y)
  erieval->PA_y[0] = data->PA_y;
#endif
#if LIBINT2_DEFINED(eri, PA_z)
  erieval->PA_z[0] = data->PA_z;
#endif

#if LIBINT2_DEFINED(eri, QC_x)
  erieval->QC_x[0] = data->QC_x;
#endif
#if LIBINT2_DEFINED(eri, QC_y)
  erieval->QC_y[0] = data->QC_y;
#endif
#if LIBINT2_DEFINED(eri, QC_z)
  erieval->QC_z[0] = data->QC_z;
#endif

#if LIBINT2_DEFINED(eri, AB_x)
  erieval->AB_x[0] = data->AB_x;
#endif
#if LIBINT2_DEFINED(eri, AB_y)
  erieval->AB_y[0] = data->AB_y;
#endif
#if LIBINT2_DEFINED(eri, AB_z)
  erieval->AB_z[0] = data->AB_z;
#endif

#if LIBINT2_DEFINED(eri, CD_x)
  erieval->CD_x[0] = data->CD_x;
#endif
#if LIBINT2_DEFINED(eri, CD_y)
  erieval->CD_y[0] = data->CD_y;
#endif
#if LIBINT2_DEFINED(eri, CD_z)
  erieval->CD_z[0] = data->CD_z;
#endif

/** Appear in OS RR for ERIs */
#if LIBINT2_DEFINED(eri, oo2z)
  /** One over 2.0*zeta */
  erieval->oo2z[0] = data->oo2z;
#endif
#if LIBINT2_DEFINED(eri, oo2e)
  /** One over 2.0*eta */
  erieval->oo2e[0] = data->oo2e;
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
  /** One over 2.0*(zeta+eta) */
  erieval->oo2ze[0] = data->oo2ze;
#endif
#if LIBINT2_DEFINED(eri, roz)
  /** rho over zeta */
  erieval->roz[0] = data->roz;
#endif
#if LIBINT2_DEFINED(eri, roe)
  /** rho over eta */
  erieval->roe[0] = data->roe;
#endif

/** Exponents */
#if LIBINT2_DEFINED(eri, zeta_A)
  erieval->zeta_A[0] = data->zeta_A;
#endif
#if LIBINT2_DEFINED(eri, zeta_B)
  erieval->zeta_B[0] = data->zeta_B;
#endif
#if LIBINT2_DEFINED(eri, zeta_C)
  erieval->zeta_C[0] = data->zeta_C;
#endif
#if LIBINT2_DEFINED(eri, zeta_D)
  erieval->zeta_D[0] = data->zeta_D;
#endif
/** Squared exponents */
#if LIBINT2_DEFINED(eri, zeta_A_2)
  erieval->zeta_A_2[0] = data->zeta_A_2;
#endif
#if LIBINT2_DEFINED(eri, zeta_B_2)
  erieval->zeta_B_2[0] = data->zeta_B_2;
#endif
#if LIBINT2_DEFINED(eri, zeta_C_2)
  erieval->zeta_C_2[0] = data->zeta_C_2;
#endif
#if LIBINT2_DEFINED(eri, zeta_D_2)
  erieval->zeta_D_2[0] = data->zeta_D_2;
#endif

// Prefactors for interelecttron transfer relation
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_x)
  erieval->TwoPRepITR_pfac0_0_x[0] = data->TwoPRepITR_pfac0_0_x;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_y)
  erieval->TwoPRepITR_pfac0_0_y[0] = data->TwoPRepITR_pfac0_0_y;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_z)
  erieval->TwoPRepITR_pfac0_0_z[0] = data->TwoPRepITR_pfac0_0_z;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_x)
  erieval->TwoPRepITR_pfac0_1_x[0] = data->TwoPRepITR_pfac0_1_x;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_y)
  erieval->TwoPRepITR_pfac0_1_y[0] = data->TwoPRepITR_pfac0_1_y;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_z)
  erieval->TwoPRepITR_pfac0_1_z[0] = data->TwoPRepITR_pfac0_1_z;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac1_0)
  erieval->TwoPRepITR_pfac1_0[0] = data->TwoPRepITR_pfac1_0;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac1_1)
  erieval->TwoPRepITR_pfac1_1[0] = data->TwoPRepITR_pfac1_1;
#endif

/** WD2004, Eq. 30, prefactor in front of (a0|k|c0) */
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_0_x)
  erieval->R12kG12_pfac0_0_x[0] = data->R12kG12_pfac0_0_x;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_0_y)
  erieval->R12kG12_pfac0_0_y[0] = data->R12kG12_pfac0_0_y;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_0_z)
  erieval->R12kG12_pfac0_0_z[0] = data->R12kG12_pfac0_0_z;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_1_x)
  erieval->R12kG12_pfac0_1_x[0] = data->R12kG12_pfac0_1_x;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_1_y)
  erieval->R12kG12_pfac0_1_y[0] = data->R12kG12_pfac0_1_y;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac0_1_z)
  erieval->R12kG12_pfac0_1_z[0] = data->R12kG12_pfac0_1_z;
#endif

/** WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0) */
#if LIBINT2_DEFINED(eri, R12kG12_pfac1_0)
  erieval->R12kG12_pfac1_0[0] = data->R12kG12_pfac1_0;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac1_1)
  erieval->R12kG12_pfac1_1[0] = data->R12kG12_pfac1_1;
#endif

/** WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0) */
#if LIBINT2_DEFINED(eri, R12kG12_pfac2)
  erieval->R12kG12_pfac2[0] = data->R12kG12_pfac2;
#endif

/** WD2004, Eq. 30, prefactor in front of curly brakets (excludes k) */
#if LIBINT2_DEFINED(eri, R12kG12_pfac3_0)
  erieval->R12kG12_pfac3_0[0] = data->R12kG12_pfac3_0;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac3_1)
  erieval->R12kG12_pfac3_1[0] = data->R12kG12_pfac3_1;
#endif

// * WD2004, Eq. 30, prefactor in front of (a0|k-2|c0)
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_0_x)
  erieval->R12kG12_pfac4_0_x[0] = data->R12kG12_pfac4_0_x;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_0_y)
  erieval->R12kG12_pfac4_0_y[0] = data->R12kG12_pfac4_0_y;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_0_z)
  erieval->R12kG12_pfac4_0_z[0] = data->R12kG12_pfac4_0_z;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_1_x)
  erieval->R12kG12_pfac4_1_x[0] = data->R12kG12_pfac4_1_x;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_1_y)
  erieval->R12kG12_pfac4_1_y[0] = data->R12kG12_pfac4_1_y;
#endif
#if LIBINT2_DEFINED(eri, R12kG12_pfac4_1_z)
  erieval->R12kG12_pfac4_1_z[0] = data->R12kG12_pfac4_1_z;
#endif

// using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = data->LIBINT_T_SS_EREP_SS0;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = data->LIBINT_T_SS_EREP_SS1;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = data->LIBINT_T_SS_EREP_SS2;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = data->LIBINT_T_SS_EREP_SS3;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = data->LIBINT_T_SS_EREP_SS4;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = data->LIBINT_T_SS_EREP_SS5;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = data->LIBINT_T_SS_EREP_SS6;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = data->LIBINT_T_SS_EREP_SS7;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = data->LIBINT_T_SS_EREP_SS8;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = data->LIBINT_T_SS_EREP_SS9;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = data->LIBINT_T_SS_EREP_SS10;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = data->LIBINT_T_SS_EREP_SS11;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = data->LIBINT_T_SS_EREP_SS12;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = data->LIBINT_T_SS_EREP_SS13;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = data->LIBINT_T_SS_EREP_SS14;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = data->LIBINT_T_SS_EREP_SS15;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = data->LIBINT_T_SS_EREP_SS16;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = data->LIBINT_T_SS_EREP_SS17;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = data->LIBINT_T_SS_EREP_SS18;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = data->LIBINT_T_SS_EREP_SS19;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = data->LIBINT_T_SS_EREP_SS20;
#endif

#if LIBINT2_DEFINED(eri, LIBINT_T_SS_K0G12_SS_0)
  erieval->LIBINT_T_SS_K0G12_SS_0[0] = data->LIBINT_T_SS_K0G12_SS_0;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_K2G12_SS_0)
  erieval->LIBINT_T_SS_K2G12_SS_0[0] = data->LIBINT_T_SS_K2G12_SS_0;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(0))
  erieval->LIBINT_T_SS_Km1G12_SS(0)[0] = data->LIBINT_T_SS_Km1G12_SS0;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(1))
  erieval->LIBINT_T_SS_Km1G12_SS(1)[0] = data->LIBINT_T_SS_Km1G12_SS1;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(2))
  erieval->LIBINT_T_SS_Km1G12_SS(2)[0] = data->LIBINT_T_SS_Km1G12_SS2;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(3))
  erieval->LIBINT_T_SS_Km1G12_SS(3)[0] = data->LIBINT_T_SS_Km1G12_SS3;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(4))
  erieval->LIBINT_T_SS_Km1G12_SS(4)[0] = data->LIBINT_T_SS_Km1G12_SS4;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(5))
  erieval->LIBINT_T_SS_Km1G12_SS(5)[0] = data->LIBINT_T_SS_Km1G12_SS5;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(6))
  erieval->LIBINT_T_SS_Km1G12_SS(6)[0] = data->LIBINT_T_SS_Km1G12_SS6;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(7))
  erieval->LIBINT_T_SS_Km1G12_SS(7)[0] = data->LIBINT_T_SS_Km1G12_SS7;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(8))
  erieval->LIBINT_T_SS_Km1G12_SS(8)[0] = data->LIBINT_T_SS_Km1G12_SS8;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(9))
  erieval->LIBINT_T_SS_Km1G12_SS(9)[0] = data->LIBINT_T_SS_Km1G12_SS9;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(10))
  erieval->LIBINT_T_SS_Km1G12_SS(10)[0] = data->LIBINT_T_SS_Km1G12_SS10;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(11))
  erieval->LIBINT_T_SS_Km1G12_SS(11)[0] = data->LIBINT_T_SS_Km1G12_SS11;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(12))
  erieval->LIBINT_T_SS_Km1G12_SS(12)[0] = data->LIBINT_T_SS_Km1G12_SS12;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(13))
  erieval->LIBINT_T_SS_Km1G12_SS(13)[0] = data->LIBINT_T_SS_Km1G12_SS13;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(14))
  erieval->LIBINT_T_SS_Km1G12_SS(14)[0] = data->LIBINT_T_SS_Km1G12_SS14;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(15))
  erieval->LIBINT_T_SS_Km1G12_SS(15)[0] = data->LIBINT_T_SS_Km1G12_SS15;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_Km1G12_SS(16))
  erieval->LIBINT_T_SS_Km1G12_SS(16)[0] = data->LIBINT_T_SS_Km1G12_SS16;
#endif
}
