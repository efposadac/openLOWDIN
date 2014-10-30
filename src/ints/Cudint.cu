#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>

const int numberOfThreads = 256;
const double pi = 3.14159265358979323846;

__global__ void intssss(int N, 
			int *primIndices_d,
			int *contIndices_d,
			double *exponents_d,
			double *primNormalization_d,
			double *coefficients_d,
			int *contCounter_d,
			int *contLength_d,
			double *origin_d,
			int *angularMoments_d,
			double *integralValues_d,
			int control,
			int kernelIter)
{
  int threadIndex = threadIdx.x + threadIdx.y*blockDim.x;
  int blockIndex = blockIdx.x;
  int global1 = threadIndex + blockIndex*blockDim.x*blockDim.y;
  int global = global1 + kernelIter; 
  
  int aa, bb, rr, ss, ii, jj, kk, ll;
  int contractionID;
  double exponentII, exponentJJ, exponentKK, exponentLL;
  double coefficientsII, coefficientsJJ, coefficientsKK, coefficientsLL;
  double primNormII, primNormJJ, primNormKK, primNormLL;
  int exponentIterII, exponentIterJJ, exponentIterKK, exponentIterLL;
  double IIx, IIy, IIz, JJx, JJy, JJz, KKx, KKy, KKz, LLx, LLy, LLz;
  double preIntegral, normIntegral; 
  double etha;
  int lAA, lBB, lRR, lSS; // Angular moments of contractions
  int integralCase;

  double A, B, C, D, KIJ, KKL, rPx, rPy, rPz, rQx, rQy, rQz, rPQ, rIJ, rKL, tFunc, tFuncsqrt, F, prefact;

  if(global1< control)
    {
      // ID of unic integrals
      contractionID = primIndices_d[global*5];

      // Contraction Indices
      aa = contIndices_d[contractionID*4];
      bb = contIndices_d[contractionID*4+1];
      rr = contIndices_d[contractionID*4+2];
      ss = contIndices_d[contractionID*4+3];
      // Primitive indices
      ii = primIndices_d[global*5+1];
      jj = primIndices_d[global*5+2];
      kk = primIndices_d[global*5+3];
      ll = primIndices_d[global*5+4];
      
      lAA = angularMoments_d[aa-1];
      lBB = angularMoments_d[bb-1];
      lRR = angularMoments_d[rr-1];
      lSS = angularMoments_d[ss-1];
      
      exponentIterII = contCounter_d[aa-1] + ii - 1;
      exponentIterJJ = contCounter_d[bb-1] + jj - 1;
      exponentIterKK = contCounter_d[rr-1] + kk - 1;
      exponentIterLL = contCounter_d[ss-1] + ll - 1;

      exponentII = exponents_d[exponentIterII];
      exponentJJ = exponents_d[exponentIterJJ];
      exponentKK = exponents_d[exponentIterKK];
      exponentLL = exponents_d[exponentIterLL];

      coefficientsII = coefficients_d[exponentIterII];
      coefficientsJJ = coefficients_d[exponentIterJJ];
      coefficientsKK = coefficients_d[exponentIterKK];
      coefficientsLL = coefficients_d[exponentIterLL];

      primNormII = primNormalization_d[exponentIterII];
      primNormJJ = primNormalization_d[exponentIterJJ];
      primNormKK = primNormalization_d[exponentIterKK];
      primNormLL = primNormalization_d[exponentIterLL];

      IIx = origin_d[(aa*3)-3];
      IIy = origin_d[(aa*3)-2];
      IIz = origin_d[(aa*3)-1];
      JJx = origin_d[(bb*3)-3];
      JJy = origin_d[(bb*3)-2];
      JJz = origin_d[(bb*3)-1];
      KKx = origin_d[(rr*3)-3];
      KKy = origin_d[(rr*3)-2];
      KKz = origin_d[(rr*3)-1];
      LLx = origin_d[(ss*3)-3];
      LLy = origin_d[(ss*3)-2];
      LLz = origin_d[(ss*3)-1];
    
      A = exponentII + exponentJJ;
      B = exponentKK + exponentLL;
      C = exponentII*exponentJJ;
      D = exponentKK*exponentLL;

      etha = (A*B)/(A+B);

      rIJ = (IIx-JJx)*(IIx-JJx) + (IIy-JJy)*(IIy-JJy) + (IIz-JJz)*(IIz-JJz);
      rKL = (KKx-LLx)*(KKx-LLx) + (KKy-LLy)*(KKy-LLy) + (KKz-LLz)*(KKz-LLz);

      KIJ = exp(-(C/A)*rIJ);
      KKL = exp(-(D/B)*rKL);

      prefact = (pi*pi*pi)/(A*B*(sqrt(A+B)));

      rPx =(exponentII*IIx+exponentJJ*JJx)/A;
      rPy =(exponentII*IIy+exponentJJ*JJy)/A;
      rPz =(exponentII*IIz+exponentJJ*JJz)/A;
      rQx = (exponentKK*KKx+exponentLL*LLx)/B;
      rQy = (exponentKK*KKy+exponentLL*LLy)/B;
      rQz = (exponentKK*KKz+exponentLL*LLz)/B;
      
      rPQ = (rPx-rQx)*(rPx-rQx) + (rPy-rQy)*(rPy-rQy) + (rPz-rQz)*(rPz-rQz);

      tFunc = etha*rPQ;

      tFuncsqrt = sqrt(tFunc);

      if(tFunc == 0.0)
	F = 2/(sqrt(pi));
      else
	F = erf(tFuncsqrt)/tFuncsqrt;
      
      integralCase = 64*lAA + 16*lBB + 4*lRR + lSS;

      switch(integralCase)
	{
	case 0: // Integral (s,s|s,s)
	  preIntegral = prefact*KIJ*KKL*F;
	  break;
	case 1: // Integral (s,s|s,p)
	  preIntegral = (2*(exponentKK*((KKx-LLx)+(rPx-rQx))*etha)*(exponentKK*((KKy-LLy)+(rPy-rQy))*etha)*(exponentKK*((KKz-LLz)+(rPz-rQz))*etha)*F)/(B*B*B);
	  break;
	case 4: // Integral (s,s|p,s)
	  preIntegral = -(2*(exponentLL*((KKx-LLx)+(rPx-rQx))*etha)*(exponentLL*((KKy-LLy)+(rPy-rQy))*etha)*(exponentLL*((KKz-LLz)+(rPz-rQz))*etha)*F)/(B*B*B);
	  break;
	case 16: // Integral (s,p|s,s)
	  preIntegral = (2*(exponentII*((IIx-JJx)+(rPx-rQx))*etha)*(exponentII*((IIy-JJy)+(rPy-rQy))*etha)*(exponentII*((IIz-JJz)+(rPz-rQz))*etha)*F)/(A*A*A);
	  break;
	case 64: // Integral (p,s|s,s)
	  preIntegral = (2*(exponentJJ*((IIx-JJx)+(rPx-rQx))*etha)*(exponentJJ*((IIy-JJy)+(rPy-rQy))*etha)*(exponentJJ*((IIz-JJz)+(rPz-rQz))*etha)*F)/(A*A*A);
	  break;
	case 5: // Integral (s,s|p,p)
	  preIntegral = (1/(4*pow(B,6))*(B+2*exponentKK*(KKx-LLx)*(-exponentLL*((KKx-LLx)+(rPx-rQx))*etha)+etha*(-1-2*exponentLL*(KKx-LLx)*(rPx-rQx)+2*(rPx-rQx)*(rPx-rQx)*etha))*(4*exponentKK*exponentKK*(KKy-LLy)*(KKz-LLz)*(exponentLL*((KKy-LLy)-(rPy-rQy))*etha)*(exponentLL*((KKz-LLz)-(rPz-rQz))*etha)+(B+etha*(-1-2*exponentLL*(KKy-LLy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)*etha))*(B+etha*(-1-2*exponentLL*(KKz-LLz)*(rPz-rQz)+2*(rPz-rQz)*(rPz-rQz)*etha))+2*exponentKK*(2*exponentLL*exponentLL*(KKy-LLy)*(KKz-LLz)*((KKz-LLz)*(rPy-rQy)+(KKy-LLy)*(rPz-rQz))*etha-exponentLL*(4*(KKy-LLy)*(KKz-LLz)*(rPy-rQy)*(rPz-rQz)*etha*etha+(KKz-LLz)*(KKz-LLz)*(B+etha*(-1+2*(rPy-rQy)*(rPy-rQy)*etha))+(KKy-LLy)*(KKy-LLy)*(B+etha*(-1+2*(rPz-rQz)*(rPz-rQz)*etha)))+etha*((KKz-LLz)*(rPz-rQz)*(B+etha*(-1+2*(rPy-rQy)*(rPy-rQy)*etha))+(KKy-LLy)*(rPy-rQy)*(B+etha*(-1+2*(rPz-rQz)*(rPz-rQz)*etha))))))*F;
	  break;
	case 17: // Integral (s,p|s,p)
	  preIntegral = (1/(4*A*A*A*A*A*B*B*B*B*exponentKK*(rPx-rQx)*exponentII))*(2*(IIx-JJx)*A*A*B*((KKx-LLx)*(rPx-rQx)+exponentKK*etha)+A*A*B*(rPx-rQx)*exponentII*(exponentKK-2*(KKx-LLx)*(rPx-rQx)-2*exponentKK*(rPx-rQx)*(rPx-rQx)*etha))*(-2*exponentII*(IIy-JJy)*(exponentKK*(KKy-LLy)+(rPy-rQy)*etha)+etha*(-1+2*exponentKK*(KKy-LLy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)*etha))*(-2*exponentII*(IIz-JJz)*(exponentKK*(KKz-LLz)+(rPz-rQz)*etha)+etha*(-1+2*exponentKK*(KKz-LLz)*(rPz-rQz)+2*(rPz-rQz)*(rPz-rQz)*etha))*F;
	  break;
	case 20: // Integral (s,p|p,s)
	  preIntegral = (1/(4*A*A*A*A*A*B*B*B*B*exponentLL*(rPx-rQx)*exponentII))*(2*(IIx-JJx)*A*A*B*(-(KKx-LLx)*(rPx-rQx)+exponentLL*etha)+A*A*B*(rPx-rQx)*exponentII*(exponentLL+2*(KKx-LLx)*(rPx-rQx)-2*exponentLL*(rPx-rQx)*(rPx-rQx)*etha))*(2*exponentII*(IIy-JJy)*(exponentLL*(KKy-LLy)-(rPy-rQy)*etha)+etha*(-1-2*exponentLL*(KKy-LLy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)*etha))*(-2*exponentII*(IIz-JJz)*(exponentLL*(KKz-LLz)-(rPz-rQz)*etha)+etha*(-1-2*exponentLL*(KKz-LLz)*(rPz-rQz)+2*(rPz-rQz)*(rPz-rQz)*etha))*F;
	  break;
	case 65: // Integral (p,s|s,p)
	  preIntegral = -(1/(4*A*A*A*A*A*B*B*B*B*exponentKK*(rPx-rQx)*exponentJJ))*(2*exponentJJ*(IIy-JJy)*(exponentKK*(KKy-LLy)+(rPy-rQy)*etha)+etha*(-1+2*exponentKK*(KKy-LLy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)*etha))*(2*exponentJJ*(IIz-JJz)*(exponentKK*(KKz-LLz)+(rPz-rQz)*etha)+etha*(-1+2*exponentKK*(KKz-LLz)*(rPz-rQz)+2*(rPz-rQz)*(rPz-rQz)*etha))*(2*(IIx-JJx)*A*A*B*((KKx-LLx)*(rPx-rQx)+exponentKK*etha)+A*A*B*(rPx-rQx)*exponentJJ*etha*(2*(KKx-LLx)*(rPx-rQx)+exponentKK*(-1+2*(rPx-rQx)*(rPx-rQx)*etha)))*F;
	  break;
	case 68: // Integral (p,s|p,s)
	  preIntegral = -1/(4*A*A*A*A*exponentLL*A*(rPx-rQx)*B*B*B*B*exponentJJ)*(2*(IIx-JJx)*A*B*((rPx-rQx)*A*(KKx-LLx)-A*exponentLL*etha)+A*(KKx-LLx)*B*exponentJJ*etha*(A*exponentLL+2*(KKx-LLx)*A*(rPx-rQx)-2*A*exponentLL*(rPx-rQx)*(rPx-rQx)*etha))*(exponentJJ*(-2*exponentLL*(IIx-JJx)*(KKx-LLx)+2*(IIy-JJy)*(rPy-rQy)*etha)+etha*(-1-2*exponentLL*(KKy-LLy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)*etha))*(2*exponentJJ*(KKz-LLz)*(exponentLL*(KKz-LLz)-(rPz-rQz)*etha)+etha*(1+2*exponentLL*(KKz-LLz)*(rPz-rQz)-2*(rPz-rQz)*(rPz-rQz)*etha))*F;
	  break;
	case 80: // Integral (p,p|s,s)
	  preIntegral = 1/(4*A*A*A*A*A*A)*(A-2*exponentII*(IIx-JJx)*(exponentJJ*(IIx-JJx)+(rPx-rQx)*etha)+etha*(2*exponentJJ*(IIx-JJx)*(rPx-rQx)-1+2*(rPx-rQx)*(rPx-rQx)*etha)*(4*exponentII*exponentII*(IIy-JJy)*(IIz-JJz)*(exponentJJ*(IIy-JJy)+(rPy-rQy)*etha)*(exponentJJ*(IIz-JJz)+(rPx-rQx)*etha)+(A+etha*(-1+2*exponentJJ*(IIy-JJy)*(rPy-rQy)+2*(rPy-rQy)*(rPy-rQy)))*(A+etha*(-1+2*exponentJJ*(IIz-JJz)*(rPz-rQz)+2*(rPz-rQz)*(rPy-rQy)))-2*exponentII*(2*exponentJJ*exponentJJ*(IIy-JJy)*(IIz-JJz)*((IIz-JJz)*(rPy-rQy)+(rPz-rQz)*((IIy-JJy))*etha)+exponentJJ*(4*(IIy-JJy)*(IIz-JJz)*(rPy-rQy)*(rPz-rQz)*etha*etha+(IIz-JJz)*(IIz-JJz)*(A+etha*(-1+2*(rPy-rQy)*(rPy-rQy)*etha))+(IIy-JJy)*(IIy-JJy)*(A+etha*(-1+2*(rPz-rQz)*(rPz-rQz)*etha)))+etha*((IIz-JJz)*(rPz-rQz)*(A+etha*(-1+2*(rPz-rQz)*(rPy-rQy)*etha))+(IIy-JJy)*(rPy-rQy)*(A+etha*(-1+2*(rPz-rQz)*(rPy-rQy)*etha))))))*F;
	  break;
	case 21: // Integral (s,p|p,p)
	  preIntegral = -1/(4*A*A*A*A*A*B*B*B*B*B*B*exponentKK*exponentLL)*(exponentII*(IIy-JJy)*(-B+2*exponentKK*(KKy-LLy)*(-(rPy-rQy)*(etha)+exponentLL*(KKy-LLy))+etha*(-2*(rPy-rQy)*(rPy-rQy)*etha+1+2*exponentLL*(KKy-LLy)*(rPy-rQy))+etha*(exponentLL*((KKy-LLy)-2*(KKy-LLy)*(rPy-rQy)*(rPy-rQy)*etha)-(KKy-LLy)*exponentKK*(1+2*exponentLL*(KKy-LLy)*(rPy-rQy)-2*(rPy-rQy)*(rPy-rQy)*etha)+(rPy-rQy)*(B+etha*(2*(rPy-rQy)*(rPy-rQy)*etha-3)))))*(exponentII*(IIz-JJz)*(-B+2*exponentKK*(KKz-LLz)*(-etha*(rPz-rQz)+(KKz-LLz)*exponentLL)+etha*(-2*(rPz-rQz)*(rPz-rQz)+1+2*exponentLL*(KKz-LLz)*(rPz-rQz)))+etha*(exponentLL*((KKz-LLz)-2*(KKz-LLz)*(rPz-rQz)*(rPz-rQz)*etha)-exponentKK*(KKz-LLz)*(1+2*exponentLL*(KKz-LLz)*(rPz-rQz)-2*(rPz-rQz)*(rPz-rQz)*etha)+(rPz-rQz)*(B+etha*(2*etha*(rPz-rQz)*(rPz-rQz)-3))))*(A*exponentLL*etha*(exponentKK*(KKx-LLx)*A*exponentKK*(-1+2*(rPx-rQx)*(rPx-rQx))+A*exponentKK*(rPx-rQx)*(B+etha*(2*(rPx-rQx)*(rPx-rQx)*etha-3))+exponentLL*(KKx-LLx)*(A*exponentKK-2*(KKx-LLx)*A*(rPx-rQx)-2*exponentKK*A*(rPx-rQx)*(rPx-rQx)))+exponentII*(IIx-JJx)*(2*exponentLL*(KKx-LLx)*(KKx-LLx)*A*A*exponentLL-2*(KKx-LLx)*A*A*exponentLL*(rPx-rQx)*etha+A*exponentKK*(2*(KKx-LLx)*A*(rPx-rQx)*etha+A*exponentLL*(etha-B-2*(rPx-rQx)*(rPx-rQx)*etha*etha))))*F;
	  break;
	case 69: // Integral (p,s|p,p)
	  preIntegral = 1/(4*A*A*A*A*A*B*B*B*B*B*B*exponentKK*exponentLL)*(exponentJJ*(IIy-JJy)*(B+2*exponentKK*(KKy-LLy)*((rPy-rQy)*(etha)-exponentLL*(KKy-LLy))+etha*(2*(rPy-rQy)*(rPy-rQy)*etha-1-2*exponentLL*(KKy-LLy)*(rPy-rQy))+etha*(exponentLL*((KKy-LLy)-2*(KKy-LLy)*(rPy-rQy)*(rPy-rQy)*etha)-(KKy-LLy)*exponentKK*(1+2*exponentLL*(KKy-LLy)*(rPy-rQy)-2*(rPy-rQy)*(rPy-rQy)*etha)+(rPy-rQy)*(B+etha*(2*(rPy-rQy)*(rPy-rQy)*etha-3)))))*(exponentJJ*(IIz-JJz)*(B+2*exponentKK*(KKz-LLz)*(etha*(rPz-rQz)-(KKz-LLz)*exponentLL)+etha*(2*(rPz-rQz)*(rPz-rQz)-1-2*exponentLL*(KKz-LLz)*(rPz-rQz)))+etha*(exponentLL*((KKz-LLz)-2*(KKz-LLz)*(rPz-rQz)*(rPz-rQz)*etha)-exponentKK*(KKz-LLz)*(1+2*exponentLL*(KKz-LLz)*(rPz-rQz)-2*(rPz-rQz)*(rPz-rQz)*etha)+(rPz-rQz)*(B+etha*(2*etha*(rPz-rQz)*(rPz-rQz)-3))))*(A*exponentLL*etha*(exponentKK*(KKx-LLx)*A*exponentKK*(1-2*(rPx-rQx)*(rPx-rQx))-A*exponentKK*(rPx-rQx)*(B+etha*(2*(rPx-rQx)*(rPx-rQx)*etha-3))+exponentLL*(KKx-LLx)*(2*(KKx-LLx)*A*(rPx-rQx)+exponentKK*A*(2*(rPx-rQx)*(rPx-rQx)*etha-1)))+exponentJJ*(IIx-JJx)*(2*exponentLL*(KKx-LLx)*(KKx-LLx)*A*A*exponentLL-2*(KKx-LLx)*A*A*exponentLL*(rPx-rQx)*etha+A*exponentKK*(2*(KKx-LLx)*A*(rPx-rQx)*etha+A*exponentLL*(etha-B-2*(rPx-rQx)*(rPx-rQx)*etha*etha))))*F;
	  break;
	case 81: // Integral (p,p|s,p)
	  preIntegral = 1/(4*A*A*A*A*A*A*B*B*B*B*B*exponentII*exponentII*exponentJJ*exponentJJ)*(-2*exponentJJ*exponentJJ*(IIx-JJx)*(IIx-JJx)*B*B*(exponentKK*(KKx-LLx) + (rPx-rQx)*etha) + exponentKK*(KKx-LLx)*(A*B*B*exponentII*exponentJJ + etha*B*B*(-exponentII*exponentJJ + 2*(IIx-JJx)*exponentII*(rPx-rQx) - 2*(IIx-JJx)*exponentJJ*(rPx-rQx) + 2*exponentII*exponentJJ*(rPx-rQx)*(rPx-rQx)*etha)) + etha*((IIx-JJx)*B*B*(exponentII-exponentJJ)*(-1 + 2*(rPx-rQx)*(rPx-rQx)*etha) + B*B*exponentII*exponentJJ*(rPx-rQx)*(A + etha*(-3 + 2*(rPx-rQx)*(rPx-rQx)*etha))))*(-2*exponentJJ*exponentJJ*(IIy-JJy)*(IIy-JJy)*B*B*(exponentKK*(KKy-LLy) + (rPy-rQy)*etha) + exponentKK*(KKy-LLy)*(A*B*B*exponentII*exponentJJ + etha*B*B*(-exponentII*exponentJJ + 2*(IIy-JJy)*exponentII*(rPy-rQy) - 2*(IIy-JJy)*exponentJJ*(rPy-rQy) + 2*exponentII*exponentJJ*(rPy-rQy)*(rPy-rQy)*etha)) + etha*((IIy-JJy)*B*B*(exponentII-exponentJJ)*(-1 + 2*(rPy-rQy)*(rPy-rQy)*etha) + B*B*exponentII*exponentJJ*(rPy-rQy)*(A + etha*(-3 + 2*(rPy-rQy)*(rPy-rQy)*etha))))*(exponentKK*(KKz-LLz)*(A + etha*(-1 + 2*exponentJJ*(IIz-JJz)*(rPz-rQz) + 2*(rPz-rQz)*(rPz-rQz)*etha)) - exponentII*(IIz-JJz)*(2*exponentJJ*(IIz-JJz)*(exponentKK*(KKz-LLz) + (rPz-rQz)*etha) + etha*(-1 + 2*exponentKK*(KKz-LLz)*(rPz-rQz) + 2*(rPz-rQz)*(rPz-rQz)*etha)) + etha*(exponentJJ*(IIz-JJz)*(-1 + 2*(rPz-rQz)*(rPz-rQz)*etha) + (rPz-rQz)*(A + etha*(-3 + 2*(rPz-rQz)*(rPz-rQz)*etha))))*F;
	  break;
	case 84: // Integral (p,p|p,s)
	  preIntegral = -1/(4*A*A*A*A*A*A*B*B*B*B*B*exponentII*exponentII*exponentJJ*exponentJJ)*(2*exponentJJ*exponentJJ*(IIx-JJx)*(IIx-JJx)*B*B*(exponentLL*(KKx-LLx)-(rPx-rQx)*etha)+etha*((IIx-JJx)*B*B*(exponentII-exponentJJ)*(2*(rPx-rQx)*(rPx-rQx)*etha-1)+B*B*exponentII*exponentJJ*(rPx-rQx)*(A+etha*(2*(rPx-rQx)*(rPx-rQx)*etha-3)))+exponentLL*(KKx-LLx)*(etha*(2*(IIx-JJx)*B*B*exponentJJ*(rPx-rQx)-A*B*B*exponentII*exponentJJ+B*B*exponentII*(exponentJJ-2*(IIx-JJx)*(rPx-rQx)-2*exponentJJ*(rPx-rQx)*(rPx-rQx)*etha))))*(2*exponentJJ*exponentJJ*(IIy-JJy)*(IIy-JJy)*B*B*(exponentLL*(KKy-LLy)-(rPy-rQy)*etha)+etha*((IIy-JJy)*B*B*(exponentII-exponentJJ)*(2*(rPy-rQy)*(rPy-rQy)*etha-1)+B*B*exponentII*exponentJJ*(rPy-rQy)*(A+etha*(2*(rPy-rQy)*(rPy-rQy)*etha-3)))+exponentLL*(KKy-LLy)*(etha*(2*(IIy-JJy)*B*B*exponentJJ*(rPy-rQy)-A*B*B*exponentII*exponentJJ+B*B*exponentII*(exponentJJ-2*(IIy-JJy)*(rPy-rQy)-2*exponentJJ*(rPy-rQy)*(rPy-rQy)*etha))))*(exponentLL*(KKz-LLz)*(A+etha*(2*exponentJJ*(IIz-JJz)*(rPz-rQz)+2*(rPz-rQz)*(rPz-rQz)*etha-1))-exponentII*(IIz-JJz)*(2*exponentJJ*(IIz-JJz)*(exponentLL*(KKz-LLz)-(rPz-rQz)*etha)+etha*(1+2*exponentLL*(KKz-LLz)*(rPz-rQz)-2*(rPz-rQz)*(rPz-rQz)*etha))+etha*(exponentJJ*((IIz-JJz)-2*(IIz-JJz)*(rPz-rQz)*(rPz-rQz)*etha)-(rPz-rQz)*(A+etha*(2*(rPz-rQz)*(rPz-rQz)*etha-3))))*F;
	  break;
	case 85: // Integral (p,p|p,p)
	  preIntegral = 1/(32*A*A*A*A*A*A*B*B*B*B*B*B)*(A*B + 2*exponentJJ*exponentLL*(IIx-JJx)*(KKx-LLx)*etha - A*etha-B*etha - 2*exponentLL*(KKx-LLx)*A*(rPx-rQx)*etha + 2*exponentJJ*(IIx-JJx)*B*(rPx-rQx)*etha + 3*etha*etha - 6*exponentJJ*(IIx-JJx)*(rPx-rQx)*etha*etha + 6*exponentLL*(KKx-LLx)*(rPx-rQx)*etha*etha - 4*exponentJJ*exponentLL*(IIx-JJx)*(KKx-LLx)*(rPx-rQx)*(rPx-rQx)*etha*etha + 2*A*(rPx-rQx)*(rPx-rQx)*etha*etha + 2*B*(rPx-rQx)*(rPx-rQx)*etha*etha - 12*(rPx-rQx)*(rPx-rQx)*etha*etha*etha + 4*exponentJJ*(IIx-JJx)*(rPx-rQx)*(rPx-rQx)*(rPx-rQx)*etha*etha*etha - 4*exponentLL*(KKx-LLx)*(rPx-rQx)*(rPx-rQx)*(rPx-rQx)*etha*etha*etha + 4*(rPx-rQx)*(rPx-rQx)*(rPx-rQx)*(rPx-rQx)*etha*etha*etha*etha - 2*exponentKK*(KKx-LLx)*(exponentLL*(KKx-LLx)*(A+etha*(-1 + 2*exponentJJ*(IIx-JJx)*(rPx-rQx) + 2*(rPx-rQx)*(rPx-rQx)*etha)) + etha*(exponentJJ*((IIx-JJx) - 2*(IIx-JJx)*(rPx-rQx)*(rPx-rQx)*etha) - (rPx-rQx)*(A + etha*(-3 + 2*(rPx-rQx)*(rPx-rQx)*etha)))) + 2*exponentII*(IIx-JJx)*(exponentJJ*(IIx-JJx)*(-B + 2*exponentKK*(KKx-LLx)*(exponentLL*(KKx-LLx) - (rPx-rQx)*etha) + etha*(1 + 2*exponentLL*(KKx-LLx)*(rPx-rQx) - 2*(rPx-rQx)*(rPx-rQx)*etha)) + etha*(exponentKK*(KKx-LLx)*(1 + 2*exponentLL*(KKx-LLx)*(rPx-rQx) - 2*(rPx-rQx)*(rPx-rQx)*etha) + exponentLL*(KKx-LLx)*(-1 + 2*(rPx-rQx)*(rPx-rQx)*etha) - (rPx-rQx)*(B + etha*(-3 + 2*(rPx-rQx)*(rPx-rQx)*etha)))))*(A*B + 2*exponentJJ*exponentLL*(IIy-JJy)*(KKy-LLy)*etha - A*etha-B*etha - 2*exponentLL*(KKy-LLy)*A*(rPy-rQy)*etha + 2*exponentJJ*(IIy-JJy)*B*(rPy-rQy)*etha + 3*etha*etha - 6*exponentJJ*(IIy-JJy)*(rPy-rQy)*etha*etha + 6*exponentLL*(KKy-LLy)*(rPy-rQy)*etha*etha - 4*exponentJJ*exponentLL*(IIy-JJy)*(KKy-LLy)*(rPy-rQy)*(rPy-rQy)*etha*etha + 2*A*(rPy-rQy)*(rPy-rQy)*etha*etha + 2*B*(rPy-rQy)*(rPy-rQy)*etha*etha - 12*(rPy-rQy)*(rPy-rQy)*etha*etha*etha + 4*exponentJJ*(IIy-JJy)*(rPy-rQy)*(rPy-rQy)*(rPy-rQy)*etha*etha*etha - 4*exponentLL*(KKy-LLy)*(rPy-rQy)*(rPy-rQy)*(rPy-rQy)*etha*etha*etha + 4*(rPy-rQy)*(rPy-rQy)*(rPy-rQy)*(rPy-rQy)*etha*etha*etha*etha - 2*exponentKK*(KKy-LLy)*(exponentLL*(KKy-LLy)*(A+etha*(-1 + 2*exponentJJ*(IIy-JJy)*(rPy-rQy) + 2*(rPy-rQy)*(rPy-rQy)*etha)) + etha*(exponentJJ*((IIy-JJy) - 2*(IIy-JJy)*(rPy-rQy)*(rPy-rQy)*etha) - (rPy-rQy)*(A + etha*(-3 + 2*(rPy-rQy)*(rPy-rQy)*etha)))) + 2*exponentII*(IIy-JJy)*(exponentJJ*(IIy-JJy)*(-B + 2*exponentKK*(KKy-LLy)*(exponentLL*(KKy-LLy) - (rPy-rQy)*etha) + etha*(1 + 2*exponentLL*(KKy-LLy)*(rPy-rQy) - 2*(rPy-rQy)*(rPy-rQy)*etha)) + etha*(exponentKK*(KKy-LLy)*(1 + 2*exponentLL*(KKy-LLy)*(rPy-rQy) - 2*(rPy-rQy)*(rPy-rQy)*etha) + exponentLL*(KKy-LLy)*(-1 + 2*(rPy-rQy)*(rPy-rQy)*etha) - (rPy-rQy)*(B + etha*(-3 + 2*(rPy-rQy)*(rPy-rQy)*etha)))))*(A*B + 2*exponentJJ*exponentLL*(IIz-JJz)*(KKz-LLz)*etha - A*etha-B*etha - 2*exponentLL*(KKz-LLz)*A*(rPz-rQz)*etha + 2*exponentJJ*(IIz-JJz)*B*(rPz-rQz)*etha + 3*etha*etha - 6*exponentJJ*(IIz-JJz)*(rPz-rQz)*etha*etha + 6*exponentLL*(KKz-LLz)*(rPz-rQz)*etha*etha - 4*exponentJJ*exponentLL*(IIz-JJz)*(KKz-LLz)*(rPz-rQz)*(rPz-rQz)*etha*etha + 2*A*(rPz-rQz)*(rPz-rQz)*etha*etha + 2*B*(rPz-rQz)*(rPz-rQz)*etha*etha - 12*(rPz-rQz)*(rPz-rQz)*etha*etha*etha + 4*exponentJJ*(IIz-JJz)*(rPz-rQz)*(rPz-rQz)*(rPz-rQz)*etha*etha*etha - 4*exponentLL*(KKz-LLz)*(rPz-rQz)*(rPz-rQz)*(rPz-rQz)*etha*etha*etha + 4*(rPz-rQz)*(rPz-rQz)*(rPz-rQz)*(rPz-rQz)*etha*etha*etha*etha - 2*exponentKK*(KKz-LLz)*(exponentLL*(KKz-LLz)*(A+etha*(-1 + 2*exponentJJ*(IIz-JJz)*(rPz-rQz) + 2*(rPz-rQz)*(rPz-rQz)*etha)) + etha*(exponentJJ*((IIz-JJz) - 2*(IIz-JJz)*(rPz-rQz)*(rPz-rQz)*etha) - (rPz-rQz)*(A + etha*(-3 + 2*(rPz-rQz)*(rPz-rQz)*etha)))) + 2*exponentII*(IIz-JJz)*(exponentJJ*(IIz-JJz)*(-B + 2*exponentKK*(KKz-LLz)*(exponentLL*(KKz-LLz) - (rPz-rQz)*etha) + etha*(1 + 2*exponentLL*(KKz-LLz)*(rPz-rQz) - 2*(rPz-rQz)*(rPz-rQz)*etha)) + etha*(exponentKK*(KKz-LLz)*(1 + 2*exponentLL*(KKz-LLz)*(rPz-rQz) - 2*(rPz-rQz)*(rPz-rQz)*etha) + exponentLL*(KKz-LLz)*(-1 + 2*(rPz-rQz)*(rPz-rQz)*etha) - (rPz-rQz)*(B + etha*(-3 + 2*(rPz-rQz)*(rPz-rQz)*etha)))))*F;
	  break;
	}

      normIntegral = primNormII*primNormJJ*primNormKK*primNormLL*preIntegral;
      integralValues_d[global1] = coefficientsII*coefficientsJJ*coefficientsKK*coefficientsLL*normIntegral;
    }
}

extern "C" void cuda_int_intraspecies_(int *numberOfContractions,
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
                                       double *contractionPrimNormalization,
				       double *contractionIntegrals,
				       int *contractionIndices)
{
  int N;
  double *integralValues, *integralValues_d;
  int a, b, r, s, u, n;
  int *contLength;
  int contractionsMem, totalPrimitives, unicintegrals, unicintegralsMem, exponentSize;
  int *contIndices, *primIndices, *contCounter;
  double *exponents, *primNormalization, *coefficients, *origin, *contNormalization, *contractedIntegrals, *integralValuesTotal;
  int *angularMoments;
  int *numberOfPPUC, contractionsMemDoub, unicintegralsMemDoub;
  int i,j,k,l,m,p;
  int auxCounter, originSize;

  //Cuda Arrays
  int *contIndices_d, *primIndices_d, *contLength_d, *contCounter_d, *angularMoments_d;
  double *exponents_d, *primNormalization_d, *coefficients_d, *origin_d;

  unicintegrals = ((*numberOfContractions*(*numberOfContractions+1)/2)+1)*(*numberOfContractions*(*numberOfContractions+1)/2)/2;

  //////////////////////////////////////////////////////////////////////
  /// Memory size
  contractionsMem = *numberOfContractions*sizeof(int);
  contractionsMemDoub = *numberOfContractions*sizeof(double);
  unicintegralsMem = unicintegrals*sizeof(int);
  unicintegralsMemDoub = unicintegrals*sizeof(double);
  exponentSize = *primNormalizationSize*sizeof(double);
  originSize = *numberOfContractions*3*sizeof(double);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  /// Malloc
  //contLength = Contraction size
  contLength = (int *)malloc(contractionsMem);
  // Counter for contractions
  contCounter = (int *)malloc(contractionsMem);
  //numberOfPPC = Number of Primitives per Unic Integral Contraction
  numberOfPPUC = (int *)malloc(unicintegralsMem);
  //Unic Integral Contraction Indices
  contIndices = (int *)malloc(4*unicintegralsMem); 
  //Exponents of contractions
  exponents = (double *)malloc(exponentSize);
  //Primitive normalization constants
  primNormalization = (double *)malloc(exponentSize);
  //Coefficients of contractions
  coefficients = (double *)malloc(exponentSize);
  // Origins of contractions
  origin = (double *)malloc(originSize);
  // Contracted Integrals
  contractedIntegrals = (double *)malloc(unicintegralsMemDoub);
  // Normalization constants of contractions
  contNormalization = (double *)malloc(contractionsMemDoub);
  // Angular moments of contractions
  angularMoments = (int *)malloc(contractionsMem);
  //////////////////////////////////////////////////////////////////////

  auxCounter = 0;
  for(i=0; i<*numberOfContractions;i++)
    {
      contNormalization[i] = *(contractionContNormalization+i);
      angularMoments[i] = *(contractionAngularMoment+i);
      // printf("Angular moments: %d\n", angularMoments[i]);
      for(j=0; j<3; j++)
	{
	  origin[j+i*3] = *(contractionOrigin+(j+i*3));
	  // printf("Origin %f \n",*(contractionOrigin+(j+i*3)));
	}
      contLength[i] = *(contractionLength+i);
      contCounter[i] = auxCounter; 
      // printf("Contraction length: %d %d\n", contLength[i], contCounter[i]);
      // printf("Origins: (%f, %f, %f)\n", origin[i*3], origin[i*3+1], origin[i*3+2]);
      auxCounter += contLength[i];
    }

  // printf("Exponents, coefficients and Primitive Normalization constants:\n");
  for(i=0; i<*primNormalizationSize;i++)
      {
	exponents[i] = *(contractionOrbitalExponents+i);
	primNormalization[i] = *(contractionPrimNormalization+i);
	coefficients[i] = *(contractionCoefficients+i);
	// printf(" (%d) %f %f %f\n", i, exponents[i], coefficients[i], primNormalization[i]);
      }

  m=0;
  totalPrimitives = 0;
  printf("NUmber of Contractions Cudint: %d\n", *numberOfContractions);
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
		  printf("Contraction C (%d): (%d,%d|%d,%d) %d\n", m, a, b, r, s, numberOfPPUC[m] );
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
			    // printf("Primitives %d, %d, %d, %d, %d\n", primIndices[5*p],
			    	   // primIndices[5*p+1],
			    	   // primIndices[5*p+2],
			    	   // primIndices[5*p+3],
			    	   // primIndices[5*p+4]);
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
  integralValuesTotal = (double *)malloc(N*sizeof(double));
  ////////////////////////////////////////////////////////////////////                                                                                                                                                                        /// Total threads in GPUs
  // printf("     *** GPU Especifications ***\n");
  int gpu, count;
  cudaDeviceProp prop;
  cudaGetDeviceCount(&count);
  int totalThreads=0;
  for (gpu = 0; gpu < count; gpu++) {
    cudaGetDeviceProperties(&prop,gpu);
    totalThreads+=prop.multiProcessorCount*prop.maxThreadsPerMultiProcessor;
  }
  ////////////////////////////////////////////////////////////////////   
  int numberOfBlocks = totalThreads/numberOfThreads;
  dim3 blockSize(16,16,1);
  dim3 gridSize(numberOfBlocks,1,1);

  ////////////////////////////////////////////////////////////////////////////
  /// CUDA Malloc
  cudaMalloc((void **)&primIndices_d, totalPrimitives*5*sizeof(int));
  cudaMalloc((void **)&contIndices_d, 4*unicintegralsMem);
  cudaMalloc((void **)&exponents_d, exponentSize);
  cudaMalloc((void **)&primNormalization_d, exponentSize);
  cudaMalloc((void **)&coefficients_d, exponentSize);
  cudaMalloc((void **)&contCounter_d, contractionsMem);
  cudaMalloc((void **)&angularMoments_d, contractionsMem);
  cudaMalloc((void **)&contLength_d, contractionsMem);
  cudaMalloc((void **)&origin_d, originSize);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  ///CUDA copy
  cudaMemcpy(primIndices_d, primIndices, totalPrimitives*5*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(contIndices_d, contIndices, 4*unicintegralsMem, cudaMemcpyHostToDevice);
  cudaMemcpy(exponents_d, exponents, exponentSize, cudaMemcpyHostToDevice);
  cudaMemcpy(primNormalization_d, primNormalization, exponentSize, cudaMemcpyHostToDevice);
  cudaMemcpy(coefficients_d, coefficients, exponentSize, cudaMemcpyHostToDevice);
  cudaMemcpy(contCounter_d, contCounter, contractionsMem, cudaMemcpyHostToDevice);
  cudaMemcpy(angularMoments_d, angularMoments, contractionsMem, cudaMemcpyHostToDevice);
  cudaMemcpy(contLength_d, contLength, contractionsMem, cudaMemcpyHostToDevice);
  cudaMemcpy(origin_d, origin, originSize, cudaMemcpyHostToDevice);
  //////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  ///Number of Calls to kernel

  int numberCallkernel = 0;
  
  i=0;
  int kernelIter = 0;
  int control2=0;
  while(control2<=totalPrimitives-1)
    {
      int control = 0;
      kernelIter = control2;
      while(control+numberOfPPUC[i]<=totalThreads && i < unicintegrals)
	{
	  control += numberOfPPUC[i];
          control2 += numberOfPPUC[i];
	  i++;
	  // printf("Control: %d %d\n",i, control);
	}
      numberCallkernel++;
      integralValues = (double *)malloc(control*sizeof(double));
      cudaMalloc((void **)&integralValues_d, control*sizeof(double));

      // printf("Control2: %d %d\n", numberCallkernel, control2);

      //      printf("Kernel Call Number: %d\n", numberCallkernel );
      intssss<<<gridSize,blockSize>>>(N, primIndices_d, contIndices_d, exponents_d, primNormalization_d, coefficients_d, contCounter_d, contLength_d, origin_d, angularMoments_d, integralValues_d, control, kernelIter);

      cudaMemcpy(integralValues, integralValues_d, control*sizeof(double),cudaMemcpyDeviceToHost);
      m=0;

      for(j=kernelIter;j<control2;j++)
	{
	  integralValuesTotal[j] = integralValues[j-kernelIter];    
	  // if(numberCallkernel==3)
	  //    printf("Integral post Kernel: %d, %d -> %f\n", j, j-kernelIter, integralValuesTotal[j]);
	}

      cudaFree(integralValues_d);
      free(integralValues);
    }

  printf("Unic Integrals:%d\n", unicintegrals);
      for(i=0; i<unicintegrals;i++)
	{
	  contractedIntegrals[i] = 0.0;
	  a = contIndices[i*4];
	  b = contIndices[i*4+1];
	  r = contIndices[i*4+2];
	  s = contIndices[i*4+3];
	  for(j=0; j<numberOfPPUC[i];j++)
	    {
	      contractedIntegrals[i] += contNormalization[a-1]*contNormalization[b-1]*contNormalization[r-1]*contNormalization[s-1]*integralValuesTotal[m];
	      *(contractionIntegrals+i) = contractedIntegrals[i];
	      *(contractionIndices+(i*4)) = a;
	      *(contractionIndices+(i*4+1)) = b;
	      *(contractionIndices+(i*4+2)) = r;
	      *(contractionIndices+(i*4+3)) = s;
	      m++;
	    }
	  // printf("%d %f %f %f %f\n", i, contNormalization[a],contNormalization[b],contNormalization[r],contNormalization[s]);
	  printf("(%d,%d|%d,%d) = %f \n", a,b,r,s,contractedIntegrals[i]);
	}

  // for(i=0;i<N;i++)
  //   printf("Integral en Host: %d %f\n", i, integralValues[i]);



  cudaFree(primIndices_d);
  cudaFree(contIndices_d);
  cudaFree(exponents_d);
  cudaFree(primNormalization_d);
  cudaFree(coefficients_d);
  cudaFree(contCounter_d);
  cudaFree(contLength_d);
  cudaFree(origin_d);
  free(integralValuesTotal);
  free(contLength);
  free(contCounter);
  free(numberOfPPUC);
  free(contIndices);
  free(exponents);
  free(primNormalization);
  free(coefficients);
  free(origin);
  free(contractedIntegrals);
  free(contNormalization);

  return;
}
