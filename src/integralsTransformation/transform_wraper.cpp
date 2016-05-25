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

	print_matrix("Coefficients", nao, nao, coeff, nao);

}