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

#include <Transform.h>

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

void four_index_trans(int nao, mat &C, vec &ERIS) {
  int i, j, k, l, ij, kl;
  mat X(nao, nao);
  mat Y(nao, nao);
  mat TMP(((nao + 1) * (nao + 2) / 2), ((nao + 1) * (nao + 2) / 2));
  TMP.zeros();
  for (i = 0, ij = 0; i < nao; i++)
    for (j = 0; j <= i; j++, ij++) {
      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          X(k, l) = X(l, k) = ERIS[Multi_Index(i, j, k, l)];
        }
      Y.zeros();
      Y = C.t() * X;
      X = Y * C;

      for (k = 0, kl = 0; k < nao; k++)
        for (l = 0; l <= k; l++, kl++) {
          TMP(kl, ij) = X(k, l);
        }
    }

  ERIS.zeros();
  X.zeros();
  Y.zeros();
  for (k = 0, kl = 0; k < nao; k++)
    for (l = 0; l <= k; l++, kl++) {
      for (i = 0, ij = 0; i < nao; i++)
        for (j = 0; j <= i; j++, ij++) X(i, j) = X(j, i) = TMP(kl, ij);
      Y = C.t() * X;
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
