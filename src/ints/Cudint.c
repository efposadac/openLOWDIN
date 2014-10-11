#include <stdio.h>
#include <stdlib.h>
#include <math.h>


extern void cuda_int_intraspecies_(int *numberOfContractions,
				   int *maxLength,
				   int *maxNumCartesianOrbital,
				   int *primNormalizationLength,
				   int *maxprimNormalizationLength,
				   int *contractionId,
				   int *contractionLength, 
				   int *contractionAngularMoment, 
				   int *contractionNumCartesianOrbital, 
				   int *contractionOwner,
				   double contractionOrigin[3][*numberOfContractions],
				   double contractionOrbitalExponents[*maxLength][*numberOfContractions],
				   double contractionCoefficients[*maxLength][*numberOfContractions],
				   double contractionContNormalization[*maxNumCartesianOrbital][*numberOfContractions],
				   double contractionPrimNormalization[*maxprimNormalizationLength][*numberOfContractions])
{
  int i;
  int j;
  printf("Estoy entrando a C!!!, %d\n", *numberOfContractions);

  for(i=0; i< *numberOfContractions; i++)
    {
      /* printf("IDs!!! %d %d %d %d %d %f %f %f \n", contractionId[i], contractionLength[i], contractionAngularMoment[i], contractionNumCartesianOrbital[i], contractionOwner[i], contractionOrigin[0][i], contractionOrigin[1][i], contractionOrigin[2][i]); */

      /* printf("\n"); */
      /* printf("Orbital exponents and coefficients in c for contraction %d\n", i); */
      /* for(j=0;j<*maxLength;j++) */
      /* 	{ */
      /* 	  printf("%f %e\n", contractionOrbitalExponents[j][i], contractionCoefficients[j][i]); */
      /* 	} */

      printf("\n");
      printf("Cont Normalization in c for contraction %d\n", i);
      for(j=0;j<primNormalizationLength[i];j++)
	{
	  printf("%f \n", contractionPrimNormalization[j][i]);
	}

    }


  return;

}


