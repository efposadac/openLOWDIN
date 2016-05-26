#include "transform_wraper.h"

using namespace std;

void print_matrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}

void c_test(double *coeff, double *ints, int nao) {
  auto start = std::chrono::system_clock::now();
  //four_index_trans(nao, coeff, ints);
  four_index_trans2(nao, coeff, ints, 0, 20, 2, 50, 0, 10, 0,nao-1);
  auto end = std::chrono::system_clock::now();
  double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
  printf("End of the computation. \n");
  std::cout <<"Elapsed time: "<< elapsed_seconds << " seconds"<<endl;
  
}
