#ifndef LAPACK_H
#define LAPACK_H


extern "C" {

   void dgemm_(char *transA,char *transB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *X, int *incx, double *beta, double *Y, int *incy);

}
#endif
