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

#include <IntTransfD.h>

inline int Multi_Index(int i, int j, int k, int l) {
  int tmp, ij, kl;
  if (i < j) {
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k < l) {
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i * (i + 1) / 2 + j;
  kl = k * (k + 1) / 2 + l;
  if (ij < kl) {
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij * (ij + 1) / 2 + kl;
}

// Fernando Posada
inline int Multi_Index_Inter(int i, int j, int k, int l, int w) {
  int tmp, ij, kl;
  if (i < j) {
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k < l) {
    tmp = k;
    k = l;
    l = tmp;
  }

  ij = i * (i + 1) / 2 + j;
  kl = k * (k + 1) / 2 + l;

  return ij * w + kl;
}

void four_index_trans(int nao, mat &C, vec &ERIS) {
  int i, j, k, l, ij, kl;
  mat X(nao, nao);
  mat Y(nao, nao);
  mat TMP(((nao + 1) * (nao + 2) / 2), ((nao + 1) * (nao + 2) / 2));
  TMP.setZero();
  for (i = 0, ij = 0; i < nao; i++)
    for (j = 0; j <= i; j++, ij++) {
      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          X(k, l) = X(l, k) = ERIS[Multi_Index(i, j, k, l)];
        }
      Y.setZero();
      Y = C.transpose() * X;
      X = Y * C;

      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          TMP(kl, ij) = X(k, l);
        }
    }

  ERIS.setZero();
  X.setZero();
  Y.setZero();
  for (k = 0, kl = 0; k < nao; k++)
    for (l = 0; l <= k; l++, kl++) {
      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++)
          X(i, j) = X(j, i) = TMP(kl, ij);
      Y = C.transpose() * X;
      X = Y * C;

      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++) {
          ERIS[Multi_Index(k, l, i, j)] = X(i, j);
        }
    }
}

static void Similarity_Transform(unsigned n, double Op[], double C[],
                                 double Work[]) {
  char *not_t = "N";
  char *yes_t = "T";
  int n_copy = n;
  double one = 1.0;

  double zero = 0;

  memset(Work, 0.0, n * n * sizeof(double));

  dgemm_(yes_t, not_t, &n_copy, &n_copy, &n_copy, &one, C, &n_copy, Op, &n_copy,
         &zero, Work, &n_copy);

  dgemm_(not_t, not_t, &n_copy, &n_copy, &n_copy, &one, Work, &n_copy, C,
         &n_copy, &zero, Op, &n_copy);
}

void four_index_trans(int nao, double C[], double ERIS[]) {
  /*Size of the ERIs tensor*/
  int m = nao * (nao + 1) / 2;
  int eris_size = m * (m + 1) / 2;
  int i, j, k, l, ij, kl;

  double *X;
  size_t x_size = nao * nao;
  X = (double *)malloc(x_size * sizeof(double));
  assert((X != NULL) && "We have problems allocating X!");
  double *Scratch;
  Scratch = (double *)malloc(x_size * sizeof(double));
  assert((Scratch != NULL) && "We have problems allocating Scratch!");

  int mm = (nao + 1) * (nao + 2) / 2;
  int tmp_size = mm * mm;
  double *TMP;
  TMP = (double *)malloc(tmp_size * sizeof(double));
  assert((TMP != NULL) && "We have problems allocating TMP!");

  memset(TMP, 0.0, tmp_size * sizeof(double));

  /* First half-transformation */
  for (i = 0, ij = 0; i < nao; i++)
    for (j = 0; j <= i; j++, ij++) {
      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          X[k * nao + l] = X[l * nao + k] = ERIS[Multi_Index(i, j, k, l)];
        }

      Similarity_Transform(nao, &X[0], &C[0], &Scratch[0]);

      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          TMP[(kl * mm + ij)] = X[(k * nao + l)];
        }
    }
  /* Second half-transformation */
  memset(ERIS, 0.0, eris_size * sizeof(double));
  memset(X, 0.0, nao * nao * sizeof(double));

  for (k = 0, kl = 0; k < nao; k++)
    for (l = 0; l <= k; l++, kl++) {
      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++)
          X[i * nao + j] = X[j * nao + i] = TMP[kl * mm + ij];

      Similarity_Transform(nao, &X[0], &C[0], &Scratch[0]);

      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++) {
          ERIS[Multi_Index(k, l, i, j)] = X[i * nao + j];
        }
    }
  free(X);
  free(Scratch);
}

void four_index_trans2(int nao, double C[], double ERIS[], int i_lower,
                       int i_upper, int j_lower, int j_upper, int k_lower,
                       int k_upper, int l_lower, int l_upper) {
  /*Size of the ERIs tensor*/
  int m = nao * (nao + 1) / 2;
  int eris_size = m * (m + 1) / 2;
  int i, j, k, l, ij, kl;

  double *X;
  size_t x_size = nao * nao;
  X = (double *)malloc(x_size * sizeof(double));
  assert((X != NULL) && "We have problems allocating X!");
  double *Scratch;
  Scratch = (double *)malloc(x_size * sizeof(double));
  assert((Scratch != NULL) && "We have problems allocating Scratch!");

  int mm = (nao + 1) * (nao + 2) / 2;
  int tmp_size = mm * mm;
  double *TMP;
  TMP = (double *)malloc(tmp_size * sizeof(double));
  assert((TMP != NULL) && "We have problems allocating TMP!");

  memset(TMP, 0.0, tmp_size * sizeof(double));

  /* First half-transformation */
  for (i = i_lower, ij = 0; i <= i_upper; i++)
    for (j = j_lower; j <= i; j++, ij++) {
      for (k = k_lower, kl = 0; k <= k_upper; k++)
        for (l = l_lower; l <= k; l++, kl++) {
          X[k * nao + l] = X[l * nao + k] = ERIS[Multi_Index(i, j, k, l)];
        }

      Similarity_Transform(nao, &X[0], &C[0], &Scratch[0]);

      for (k = k_lower, kl = 0; k <= k_upper; k++)
        for (l = l_lower; l <= k; l++, kl++) {
          TMP[(kl * mm + ij)] = X[(k * nao + l)];
        }
    }
  /* Second half-transformation */
  memset(ERIS, 0.0, eris_size * sizeof(double));
  memset(X, 0.0, nao * nao * sizeof(double));

  for (k = k_lower, kl = 0; k <= k_upper; k++)
    for (l = l_lower; l <= k; l++, kl++) {
      for (i = i_lower, ij = 0; i <= k_upper; i++)
        for (j = j_lower; j <= i; j++, ij++)
          X[i * nao + j] = X[j * nao + i] = TMP[kl * mm + ij];

      Similarity_Transform(nao, &X[0], &C[0], &Scratch[0]);

      for (i = i_lower, ij = 0; i <= k_upper; i++)
        for (j = j_lower; j <= i; j++, ij++) {
          ERIS[Multi_Index(k, l, i, j)] = X[i * nao + j];
        }
    }
  free(X);
  free(Scratch);
  free(TMP);
}

// Fernando Posada
void four_index_trans_inter(int nao, int onao, double C[], double OC[],
                            double ERIS[]) {

  /*Size of the ERIs tensor*/
  int m = nao * (nao + 1) / 2;
  int om = onao * (onao + 1) / 2;
  int eris_size = m * om;
  int i, j, k, l, ij, kl;

  // printf("%d %d %d\n", eris_size, nao, onao);
  
  double *OX;
  size_t ox_size = onao * onao;
  OX = (double *)malloc(ox_size * sizeof(double));
  assert((OX != NULL) && "We have problems allocating X!");

  double *OScratch;
  OScratch = (double *)malloc(ox_size * sizeof(double));
  assert((OScratch != NULL) && "We have problems allocating Scratch!");

  double *TMP;
  TMP = (double *)malloc(eris_size * sizeof(double));
  assert((TMP != NULL) && "We have problems allocating TMP!");

  memset(TMP, 0.0, eris_size * sizeof(double));
  memset(OX, 0.0, nao * nao * sizeof(double));

  /* First half-transformation */
  for (i = 0, ij = 0; i < nao; i++)
    for (j = 0; j <= i; j++, ij++) {
      for (k = 0, kl = 0; k < onao; k++)
        for (l = 0; l <= k; l++, kl++) {
          OX[k * onao + l] = OX[l * onao + k] =
              ERIS[Multi_Index_Inter(i, j, k, l, om)];
        }

      Similarity_Transform(onao, &OX[0], &OC[0], &OScratch[0]);

      for (k = 0, kl = 0; k < onao; k++)
        for (l = 0; l <= k; l++, kl++) {
          TMP[(ij * om + kl)] = OX[(k * onao + l)];
        }
    }

  free(OX);
  free(OScratch);

  /* Second half-transformation */

  double *X;
  size_t x_size = nao * nao;
  X = (double *)malloc(x_size * sizeof(double));
  assert((X != NULL) && "We have problems allocating X!");

  double *Scratch;
  Scratch = (double *)malloc(x_size * sizeof(double));
  assert((Scratch != NULL) && "We have problems allocating Scratch!");

  memset(ERIS, 0.0, eris_size * sizeof(double));
  memset(X, 0.0, nao * nao * sizeof(double));

  for (k = 0, kl = 0; k < onao; k++)
    for (l = 0; l <= k; l++, kl++) {
      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++)
          X[i * nao + j] = X[j * nao + i] = TMP[ij * om + kl];

      Similarity_Transform(nao, &X[0], &C[0], &Scratch[0]);

      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++) {
          ERIS[Multi_Index_Inter(i, j, k, l, om)] = X[i * nao + j];
        }
    }

  free(X);
  free(Scratch);
  free(TMP);
}

// Fortran interface

void c_integrals_transform_partial(double *coeff, double *ints, int nao, int lp,
                                   int up, int lq, int uq, int lr, int ur,
                                   int ls, int us) {

  // auto start = std::chrono::system_clock::now();
  four_index_trans2(nao, coeff, ints, lp, up, lq, uq, lr, ur, ls, us);
  // auto end = std::chrono::system_clock::now();

  // double elapsed_seconds =
  //     std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
  //         .count();

  // printf("End of the computation. \n");
  // std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << endl;
}

void c_integrals_transform_all(double *coeff, double *ints, int nao) {
  // auto start = std::chrono::system_clock::now();
  four_index_trans(nao, coeff, ints);
  // auto end = std::chrono::system_clock::now();

  // double elapsed_seconds =
  //     std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
  //         .count();

  // printf("End of the computation. \n");
  // std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << endl;
}

void c_integrals_transform_inter_all(double *coeff, double *ocoeff,
                                     double *ints, int nao, int onao) {
  // auto start = std::chrono::system_clock::now();
  four_index_trans_inter(nao, onao, coeff, ocoeff, ints);
  // auto end = std::chrono::system_clock::now();

  // double elapsed_seconds =
  //     std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
  //         .count();

  // printf("End of the computation. \n");
  // std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << endl;
}