/********************************************************-
  This code is part of LOWDIN Quantum chemistry package

  this program has been developed under direction of:

  Prof. A REYES' Lab. Universidad Nacional de Colombia
    http://www.qcc.unal.edu.co
  Prof. R. FLORES' Lab. Universidad de Guadalajara
    http://www.cucei.udg.mx/~robertof

    Todos los derechos reservados, 2013

**********************************************************/

#include "TransformWraper.h"

using namespace std;

void c_integrals_transform(double *coeff, double *ints, int nao, int lp, int up,
                           int lq, int uq, int lr, int ur, int ls, int us) {
  auto start = std::chrono::system_clock::now();
  // four_index_trans(nao, coeff, ints);
  // four_index_trans2(nao, coeff, ints, lp, up, lq, uq, lr, ur, ls, us)
  four_index_trans2(nao, coeff, ints, lp, up, lq, uq, lr, ur, ls, us);
  auto end = std::chrono::system_clock::now();
  double elapsed_seconds =
      std::chrono::duration_cast<std::chrono::duration<double> >(end - start)
          .count();
  printf("End of the computation. \n");
  std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << endl;
}
