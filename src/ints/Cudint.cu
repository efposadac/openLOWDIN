#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>

const int numberOfThreads = 64;

__global__ void intssss(int N, int *primIndices_d, int *contIndices_d, double *integralValues_d)
{
  int threadIndex = threadIdx.x + threadIdx.y*blockDim.x;
  int blockIndex = blockIdx.x;
  int global = threadIndex + blockIndex*blockDim.x*blockDim.y;
  
  // if(global<N)
    
}

extern "C" void cuda_int_intraspecies_(int *numberOfContractions,
                                       int *maxLength,
                                       int *maxNumCartesianOrbital,
				       int *primNormalizationSize,
                                       int *contractionId,
                                       int *contractionLength,
                                       int *contractionAngularMoment,
                                       int *contractionNumCartesianOrbital,
                                       int *contractionOwner,
                                       double *contractionOrigin,
                                       double *contractionOrbitalExponents,
                                       double *contractionCoefficients,
                                       double *contractionContNormalization,
                                       double *contractionPrimNormalization)
{
  int N;
  double *integralValues, *integralValues_d;
  int a, b, r, s, u, n;
  int *contLength;
  int contractionsMem, totalPrimitives, unicintegrals, unicintegralsMem, exponentSize;
  int *contIndices, *primIndices;
  double *exponents;
  int *numberOfPPUC;
  int i,j,k,l,m,p;

  //Cuda Arrays
  int *contIndices_d, *primIndices_d;

  unicintegrals = ((*numberOfContractions*(*numberOfContractions+1)/2)+1)*(*numberOfContractions*(*numberOfContractions+1)/2)/2;

  contractionsMem = *numberOfContractions*sizeof(int);
  unicintegralsMem = unicintegrals*sizeof(int);

  //////////////////////////////////////////////////////////////////////
  /// Malloc
  //contLength = Contraction size
  contLength = (int *)malloc(contractionsMem);
  //numberOfPPC = Number of Primitives per Unic Integral Contraction
  numberOfPPUC = (int *)malloc(unicintegralsMem);
  //Unic Integral Contraction Indices
  contIndices = (int *)malloc(4*unicintegralsMem); 
  //////////////////////////////////////////////////////////////////////

  exponentSize = 0;
  for(i=0; i<*numberOfContractions;i++)
    {
      contLength[i] = *(contractionLength+i);
      exponentSize += contLength[i];
      printf("Contraction length: %d %d\n", contLength[i], exponentSize);
    }

  int counter=0;
  exponents = (double *)malloc(exponentSize*sizeof(double));
  m=0;
  for(i=0; i<*numberOfContractions;i++){
    for(j=0; j<*maxLength;j++)
      {
	printf(" (%d, %d) %d, %d %f",i,j, *maxLength, *numberOfContractions, *(contractionOrbitalExponents+(j+i*(*numberOfContractions))));
	m++;	     }
	  printf("\n");
      }

	  m=0;

  for(i=0; i<*numberOfContractions;i++)
    {
      for(j=counter;j<(counter+contLength[i]);j++)
	{
	  exponents[m] = *(contractionOrbitalExponents+(j+i*(*numberOfContractions)));
	  printf("Exponent: %f %d, %d\n", exponents[m], m, j);
	  m++;
	}
      counter += *maxLength;
    }

  m=0;
  totalPrimitives = 0;
  for( a = 1;  a<=*numberOfContractions; a++)
    {
      n = a;
      for( b = a; b<=*numberOfContractions;b++)
  	{
          u = b;
          for( r = n ;r <=*numberOfContractions;r++)
  	    {
  	      for( s = u; s<=*numberOfContractions; s++)
  		{
		  contIndices[m*4] = a;
		  contIndices[m*4+1] = b;
		  contIndices[m*4+2] = r;
		  contIndices[m*4+3] = s;
		  numberOfPPUC[m] = contLength[a-1]*contLength[b-1]*contLength[r-1]*contLength[s-1];
		  totalPrimitives += numberOfPPUC[m];
		  printf("Primitives per contraction: %d, %d, %d, %d, %d\n", numberOfPPUC[m], a, b, r, s );
		  m++;
  		}
  	      u = r+1;
  	    }
  	}
    }

  m=0;
  p=0;
  primIndices = (int *)malloc(totalPrimitives*5*sizeof(int));
  for( a = 1;  a<=*numberOfContractions; a++)
    {
      n = a;
      for( b = a; b<=*numberOfContractions;b++)
  	{
          u = b;
          for( r = n ;r <=*numberOfContractions;r++)
  	    {
  	      for( s = u; s<=*numberOfContractions; s++)
  		{
		  for(i=1;i<=contLength[a-1];i++)
		    for(j=1;j<=contLength[b-1];j++)
		      for(k=1;k<=contLength[r-1];k++)
			for(l=1;l<=contLength[s-1];l++)
			  {
			    primIndices[5*p] = m;
			    primIndices[5*p+1] = i;
			    primIndices[5*p+2] = j;
			    primIndices[5*p+3] = k;
			    primIndices[5*p+4] = l;
			    // printf("%d, %d, %d, %d\n",i,j,k,l);
			    printf("Primitives %d, %d, %d, %d, %d\n", primIndices[5*p],
			    	   primIndices[5*p+1],
			    	   primIndices[5*p+2],
			    	   primIndices[5*p+3],
			    	   primIndices[5*p+4]);
			    p++;
			  }
		  m++;
  		}
  	      u = r+1;
  	    }
  	}
    }


  // printf("Total Primitive: %d\n", totalPrimitives);

  N=totalPrimitives;	  
  integralValues = (double *)malloc(totalPrimitives*sizeof(double));

  ////////////////////////////////////////////////////////////////////////////
  /// CUDA Malloc
  cudaMalloc((void **)&integralValues_d, totalPrimitives*sizeof(double));
  cudaMalloc((void **)&primIndices_d, totalPrimitives*5*sizeof(int));
  cudaMalloc((void **)&contIndices_d, 4*unicintegralsMem);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  ///CUDA copy
  cudaMemcpy(primIndices_d, primIndices, totalPrimitives*5*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(contIndices_d, contIndices, 4*unicintegralsMem, cudaMemcpyHostToDevice);
  //////////////////////////////////////////////////////////////////////////

  dim3 blockSize(8,8,1);
  dim3 gridSize(360,1,1);

  intssss<<<gridSize,blockSize>>>(N, primIndices_d, contIndices_d, integralValues);

  cudaMemcpy(integralValues, integralValues_d, totalPrimitives*sizeof(double),cudaMemcpyDeviceToHost);

  cudaFree(integralValues_d);
  free(integralValues);
}
