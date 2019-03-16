!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>@brief Modulo para calculo de integrales de overlap
!!
!! Este modulo contiene los algoritmos necesarios para la evaluacian de integrales
!! de overlap entre pares de funciones gaussianas primitivas (PrimitiveGaussian_),
!! sin normalizar.
!!
!! \f[ (\bf{a} \mid \bf{b}) = \int_{TE} {{\varphi(\bf{r};{\zeta}_a,
!! \bf{a,A})}{\varphi(\bf{r};{\zeta}_b,\bf{b,B})}}\,dr \f]
!!
!! Donde:
!!
!! <table>
!! <tr> <td> \f$ \zeta \f$ : <td> <dfn> exponente orbital. </dfn>
!! <tr> <td> <b> r </b> : <td> <dfn> coordenas espaciales de la funcian. </dfn>
!! <tr> <td> <b> n </b> : <td> <dfn> indice de momento angular. </dfn>
!! <tr> <td> <b> R </b> : <td> <dfn> origen de la funcion gaussiana cartesiana. </dfn>
!! </table>
!!
!! La integral de traslapamiento entre funciones gaussinas primitivas se calcula
!! con el matodo recursivo propuesto por Obara-Sayka, el cual transforma las
!! gaussiana de entrada en una sola mediante la identidad del producto gausiano.
!! La expresian general para el calculos de las integrales es:
!!
!! \f[ ({\bf{a + 1_i}} \parallel {\bf{b}}) = (P_i -A_i) ({\bf{a}} \parallel
!! {\bf{b}} ) + \frac{1}{2 \zeta} N_i(\bf{a}) ({\bf{a - 1_i}} \parallel
!! {\bf{b}} ) + \frac{1}{2 \zeta} N_i(\bf{b}) + ({\bf{a}} \parallel {\bf{b-1_i}})\f]
!!
!! Los parametros <b> P </b> y \f$ \zeta \f$ de la expresian provienen del producto de dos
!! funciones gaussianas.
!!
!! @author E. F. Posada
!!
!! <b> Fecha de creacion : </b> 2010-03-11
!!
!! <b> Historial de modificaciones: </b>
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Implementa los metodos para el calculo de integrales de traslape para cualquier valor de n
!!
module ThreeCOverlapIntegrals_
  use Exception_
  use Math_
  use ContractedGaussian_
  implicit none

  public ::  &
       ThreeCOverlapIntegrals_computeShell
  
  private :: &
       ThreeCOverlapIntegrals_obaraSaikaRecursion, &
       ThreeCOverlapIntegrals_computeContractionPair
  
  
contains

  !>
  !! @brief Calculates overlap integral between two contractions (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2011.02.04: E.F.Posada: chage for usage on opints
  !! @return  output: overlap integral of a shell (all combinations)
  !! @version 1.0
  subroutine ThreeCOverlapIntegrals_computeShell(contractedGaussianA, contractedGaussianB, contractedGaussianC, integral)
    implicit none
    
    type(ContractedGaussian), intent(in) :: contractedGaussianA, contractedGaussianB, contractedGaussianC
    real(8), intent(inout) :: integral(contractedGaussianA%numCartesianOrbital * contractedGaussianB%numCartesianOrbital * &
                                       contractedGaussianC%numCartesianOrbital)

    integer ::  am1(0:3)
    integer ::  am2(0:3)
    integer ::  am3(0:3)
    integer ::  nprim1
    integer ::  nprim2
    integer ::  nprim3
    real(8) ::  A(0:3)
    real(8) ::  B(0:3)
    real(8) ::  C(0:3)
    real(8) ::  exp1(0:contractedGaussianA%length)
    real(8) ::  exp2(0:contractedGaussianB%length)
    real(8) ::  exp3(0:contractedGaussianC%length)
    real(8) ::  coef1(0:contractedGaussianA%length)
    real(8) ::  coef2(0:contractedGaussianB%length)
    real(8) ::  coef3(0:contractedGaussianC%length)
    real(8) ::  nor1(0:contractedGaussianA%length)
    real(8) ::  nor2(0:contractedGaussianB%length)
    real(8) ::  nor3(0:contractedGaussianC%length)
    real(8) :: auxIntegral
    integer, allocatable :: angularMomentIndexA(:,:)
    integer, allocatable :: angularMomentIndexB(:,:)
    integer, allocatable :: angularMomentIndexC(:,:)
    integer ::  i, m, p, q, r
    
    integral = 0.0_8

    if(allocated(angularMomentIndexA)) deallocate(angularMomentIndexA)
    if(allocated(angularMomentIndexB)) deallocate(angularMomentIndexB)
    if(allocated(angularMomentIndexC)) deallocate(angularMomentIndexC)

    allocate(angularMomentIndexA(3, contractedGaussianA%numCartesianOrbital))
    allocate(angularMomentIndexB(3, contractedGaussianB%numCartesianOrbital))
    allocate(angularMomentIndexC(3, contractedGaussianC%numCartesianOrbital))
    
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexA, contractedGaussianA)
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexB, contractedGaussianB)
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexC, contractedGaussianC)

    nprim1 = contractedGaussianA%length
    A(0) = contractedGaussianA%origin(1)
    A(1) = contractedGaussianA%origin(2)
    A(2) = contractedGaussianA%origin(3)
    coef1(0:nprim1-1) =  contractedGaussianA%contractionCoefficients(1:nprim1)
    
    
    nprim2 = contractedGaussianB%length
    B(0) = contractedGaussianB%origin(1)
    B(1) = contractedGaussianB%origin(2)
    B(2) = contractedGaussianB%origin(3)
    coef2(0:nprim2-1) =  contractedGaussianB%contractionCoefficients(1:nprim2)

    nprim3 = contractedGaussianC%length
    C(0) = contractedGaussianC%origin(1)
    C(1) = contractedGaussianC%origin(2)
    C(2) = contractedGaussianC%origin(3)
    coef3(0:nprim3-1) =  contractedGaussianC%contractionCoefficients(1:nprim3)
 
    
    m = 0
    
    do p = 1, contractedGaussianA%numcartesianOrbital
      do q = 1, contractedGaussianB%numcartesianOrbital
        do r = 1, contractedGaussianC%numcartesianOrbital

          m = m + 1
           
          exp1(0:nprim1-1) = contractedGaussianA%orbitalExponents(1:nprim1)
          nor1(0:nprim1-1) = contractedGaussianA%primNormalization(1:nprim1,p)
           
          exp2(0:nprim2-1) = contractedGaussianB%orbitalExponents(1:nprim2)
          nor2(0:nprim2-1) = contractedGaussianB%primNormalization(1:nprim2,q)

          exp3(0:nprim3-1) = contractedGaussianC%orbitalExponents(1:nprim3)
          nor3(0:nprim3-1) = 1!contractedGaussianC%primNormalization(1:nprim3,r)

          am1 = 0
          am2 = 0
          am3 = 0

          am1(0:2) = angularMomentIndexA(1:3, p)
          am2(0:2) = angularMomentIndexB(1:3, q)
          am3(0:2) = angularMomentIndexC(1:3, r)
           
          call ThreeCOverlapIntegrals_computeContractionPair(am1, am2, am3, nprim1, nprim2, nprim3, A, B, C, exp1, exp2, exp3, &
                                                         coef1, coef2, coef3, nor1, nor2, nor3, auxIntegral)

          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q) !&
               !* contractedGaussianC%contNormalization(r)
           
          integral(m) = auxIntegral

        end do
      end do
    end do

  end subroutine ThreeCOverlapIntegrals_computeShell

  !>
  !!@brief Evalua integrales overlap para cualquier momento angular
  !!@author E. F. Posada, 2010
  !!@version 2.0
  !!@return devuelve los valores de integrales de overlap para una capa o individual
  !!@param A, B : origin A and B
  subroutine ThreeCOverlapIntegrals_computeContractionPair(angularMomentIndexA, angularMomentIndexB, angularMomentIndexC, &
                                                       lengthA, lengthB, lengthC, A, B, C, &
       orbitalExponentsA, orbitalExponentsB, orbitalExponentsC, &
       contractionCoefficientsA, contractionCoefficientsB, contractionCoefficientsC, &
       normalizationConstantsA, normalizationConstantsB, normalizationConstantsC, integralValue)
    implicit none

    integer, intent(in) :: angularMomentIndexA(0:3), angularMomentIndexB(0:3), angularMomentIndexC(0:3)
    integer, intent(in) :: lengthA, lengthB, lengthC
    real(8), intent(in) :: A(0:3), B(0:3), C(0:3)
    real(8), intent(in) :: orbitalExponentsA(0:lengthA), orbitalExponentsB(0:lengthB), orbitalExponentsC(0:lengthC)
    real(8), intent(in) :: contractionCoefficientsA(0:lengthA), contractionCoefficientsB(0:lengthB), contractionCoefficientsC(0:lengthC)
    real(8), intent(in) :: normalizationConstantsA(0:lengthA), normalizationConstantsB(0:lengthB), normalizationConstantsC(0:lengthC)
    real(8), intent(out):: integralValue

    real(8), allocatable ::  x(:,:,:), y(:,:,:), z(:,:,:)
    real(8) :: AB2, PC2
    real(8) :: auxExponentA, auxCoefficientA, auxConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxConstantB
    real(8) :: auxExponentC, auxCoefficientC, auxConstantC
    real(8) :: gamma, gammaC, gammaInv, gammaInvC
    real(8) :: P(0:3), G(0:3)
    real(8) :: GA(0:3), GB(0:3), GC(0:3)
    real(8) :: commonPreFactor, commonPreFactorABC
    real(8) :: x0, y0, z0

    integer :: angularMomentA, angularMomentB, angularMomentC
    integer :: maxAngularMoment
    integer :: p1, p2, p3 !< iteradores
    
    angularMomentA = sum(angularMomentIndexA)
    angularMomentB = sum(angularMomentIndexB)
    angularMomentC = sum(angularMomentIndexC)

    integralValue = 0.0_8

    maxAngularMoment = max(angularMomentA, angularMomentB, angularMomentC) + 1 !!check

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2, 0:maxAngularMoment+2)) !!check
    allocate(y(0:maxAngularMoment+2, 0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(z(0:maxAngularMoment+2, 0:maxAngularMoment+2, 0:maxAngularMoment+2))

    x = 0.0_8
    y = 0.0_8
    z = 0.0_8

    AB2 = 0.0_8
    AB2 = AB2 + (A(0) - B(0)) * (A(0) - B(0))
    AB2 = AB2 + (A(1) - B(1)) * (A(1) - B(1))
    AB2 = AB2 + (A(2) - B(2)) * (A(2) - B(2))

    do p1=0, lengthA - 1
       auxExponentA = orbitalExponentsA(p1)
       auxCoefficientA = contractionCoefficientsA(p1)
       auxConstantA = normalizationConstantsA(p1)
       do p2=0, lengthB - 1
          auxExponentB = orbitalExponentsB(p2)
          auxCoefficientB = contractionCoefficientsB(p2)
          auxConstantB = normalizationConstantsB(p2)
          gamma = auxExponentA + auxExponentB
          gammaInv = 1.0/gamma

          P(0) = (auxExponentA*A(0) + auxExponentB*B(0))*gammaInv
          P(1) = (auxExponentA*A(1) + auxExponentB*B(1))*gammaInv
          P(2) = (auxExponentA*A(2) + auxExponentB*B(2))*gammaInv

          commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*gammaInv) &
               * sqrt(Math_PI*gammaInv) * Math_PI * gammaInv &
               * auxCoefficientA * auxCoefficientB * auxConstantA * auxConstantB

          do p3 = 0, lengthC -1 

            auxExponentC = orbitalExponentsC(p3)
            auxCoefficientC = contractionCoefficientsC(p3)
            auxConstantC = normalizationConstantsC(p3) !?

            PC2 = 0.0_8
            PC2 = PC2 + (P(0) - C(0)) * (P(0) - C(0))
            PC2 = PC2 + (P(1) - C(1)) * (P(1) - C(1))
            PC2 = PC2 + (P(2) - C(2)) * (P(2) - C(2))

            gammaC = gamma + auxExponentC
            gammaInvC = 1.0 / gammaC

            G(0) = ( gamma*P(0) + auxExponentC*C(0))*gammaInvC
            G(1) = ( gamma*P(1) + auxExponentC*C(1))*gammaInvC
            G(2) = ( gamma*P(2) + auxExponentC*C(2))*gammaInvC

            GA(0) = G(0) - A(0)
            GA(1) = G(1) - A(1)
            GA(2) = G(2) - A(2)
            GB(0) = G(0) - B(0)
            GB(1) = G(1) - B(1)
            GB(2) = G(2) - B(2)
            GC(0) = G(0) - C(0)
            GC(1) = G(1) - C(1)
            GC(2) = G(2) - C(2)

            commonPreFactorABC = exp(-gamma*auxExponentC*gammaInvC * PC2 ) * &
                                 sqrt(gamma*gammaInvC) * (gamma*gammaInvC) * commonPreFactor * &
                                 auxCoefficientC * auxConstantC 
          !! recursion
          call ThreeCOverlapIntegrals_obaraSaikaRecursion(x, y, z, GA, GB, GC, gammaC, angularMomentA+2, angularMomentB+2, angularMomentC+2)
          
          x0 = x(angularMomentIndexA(0),angularMomentIndexC(0),angularMomentIndexB(0))
          y0 = y(angularMomentIndexA(1),angularMomentIndexC(1),angularMomentIndexB(1))
          z0 = z(angularMomentIndexA(2),angularMomentIndexC(2),angularMomentIndexB(2))

          integralValue = integralValue + commonPreFactorABC*x0*y0*z0

         end do
      end do
    end do
    
    deallocate(x)
    deallocate(y)
    deallocate(z)

  end subroutine ThreeCOverlapIntegrals_computeContractionPair

  !>
  !!@brief Implementation of recursion proposed by Obara-Saika for overlap integrals.
  !!@author Edwin Posada, 2010
  !!@return x, y, z : recursion matrixes
  !!@param PA, PB : reduced origin for gaussian A and B
  !!@param gamma : reduced exponent
  !!@see Gaussian product: if you want to know what reduced exponent and origin is.
  subroutine ThreeCOverlapIntegrals_obaraSaikaRecursion(x, y, z, GA, GB, GC, gamma, angularMoment1, angularMoment2, angularMoment3)
    implicit none

    real(8), intent(inout), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(8), intent(in) :: GA(0:3), GB(0:3), GC(0:3)
    real(8), intent(in) :: gamma
    integer, intent(in) :: angularMoment1, angularMoment2, angularMoment3

    real(8) :: pp
    integer :: i, j, k

    pp = 1/(2*gamma)

    x(0,0,0) = 1.0_8
    y(0,0,0) = 1.0_8
    z(0,0,0) = 1.0_8

    !! 1 up from the bottom.

    x(1,0,0) = GA(0)
    y(1,0,0) = GA(1)
    z(1,0,0) = GA(2)

    x(0,1,0) = GB(0)
    y(0,1,0) = GB(1)
    z(0,1,0) = GB(2)

    x(0,0,1) = GC(0)
    y(0,0,1) = GC(1)
    z(0,0,1) = GC(2)

    !! Begin - Upward recursion in b for a=c=0

    do j=1, angularMoment2 -1
       x(0,0,j+1) = GB(0)*x(0,0,j)
       y(0,0,j+1) = GB(1)*y(0,0,j)
       z(0,0,j+1) = GB(2)*z(0,0,j)
       x(0,0,j+1) = x(0,0,j+1) + j*pp*x(0,0,j-1)
       y(0,0,j+1) = y(0,0,j+1) + j*pp*y(0,0,j-1)
       z(0,0,j+1) = z(0,0,j+1) + j*pp*z(0,0,j-1)
    end do
    !! End - Upward recursion in b for a=c=0

    !! Begin - Upward recursion in c for all b and a=0

    do j=1, angularMoment2 
       x(0,1,j) = GC(0)*x(0,0,j)
       y(0,1,j) = GC(1)*y(0,0,j)
       z(0,1,j) = GC(2)*z(0,0,j)
       x(0,1,j) = x(0,1,j+1) + j*pp*x(0,0,j-1)
       y(0,1,j) = y(0,1,j+1) + j*pp*y(0,0,j-1)
       z(0,1,j) = z(0,1,j+1) + j*pp*z(0,0,j-1)
    end do

    do k=1, angularMoment3 -1
       x(0,k+1,0) = GC(0)*x(0,k,0)
       y(0,k+1,0) = GC(1)*y(0,k,0)
       z(0,k+1,0) = GC(2)*z(0,k,0)
       x(0,k+1,0) = x(0,k+1,0) + k*pp*x(0,k-1,0)
       y(0,k+1,0) = y(0,k+1,0) + k*pp*y(0,k-1,0)
       z(0,k+1,0) = z(0,k+1,0) + k*pp*z(0,k-1,0)

      do j = 1, angularMoment2 

         x(0,k+1,j) = GC(0)*x(0,k,j)
         y(0,k+1,j) = GC(1)*y(0,k,j)
         z(0,k+1,j) = GC(2)*z(0,k,j)
         x(0,k+1,j) = x(0,k+1,j) + k*pp*x(0,k-1,j)
         y(0,k+1,j) = y(0,k+1,j) + k*pp*y(0,k-1,j)
         z(0,k+1,j) = z(0,k+1,j) + k*pp*z(0,k-1,j)
         x(0,k+1,j) = x(0,k+1,j) + j*pp*x(0,k,j-1)
         y(0,k+1,j) = y(0,k+1,j) + j*pp*y(0,k,j-1)
         z(0,k+1,j) = z(0,k+1,j) + j*pp*z(0,k,j-1)
      end do
    end do

    !! End - Upward recursion in c for all b and a = 0

    !! Begin - Upward recursion in a for all b and c
    do j = 1, angularMoment2 
       x(1,0,j) = GA(0)*x(0,0,j)
       y(1,0,j) = GA(1)*y(0,0,j)
       z(1,0,j) = GA(2)*z(0,0,j)
       x(1,0,j) = x(1,0,j) + j*pp*x(0,0,j-1)
       y(1,0,j) = y(1,0,j) + j*pp*y(0,0,j-1)
       z(1,0,j) = z(1,0,j) + j*pp*z(0,0,j-1)
    end do
    do k = 1, angularMoment3
       x(1,k,0) = GA(0)*x(0,k,0)
       y(1,k,0) = GA(1)*y(0,k,0)
       z(1,k,0) = GA(2)*z(0,k,0)
       x(1,k,0) = x(1,k,0) + k*pp*x(0,k-1,0)
       y(1,k,0) = y(1,k,0) + k*pp*y(0,k-1,0)
       z(1,k,0) = z(1,k,0) + k*pp*z(0,k-1,0)
    end do
    do i = 1, angularMoment1 - 1
       x(i+1,0,0) = GA(0)*x(i,0,0)
       y(i+1,0,0) = GA(1)*y(i,0,0)
       z(i+1,0,0) = GA(2)*z(i,0,0)
       x(i+1,0,0) = x(i+1,0,0) + i*pp*x(i-1,0,0)
       y(i+1,0,0) = y(i+1,0,0) + i*pp*y(i-1,0,0)
       z(i+1,0,0) = z(i+1,0,0) + i*pp*z(i-1,0,0)
      do j = 1, angularMoment2 - 1
         x(i+1,0,j) = GA(0)*x(i,0,j)
         y(i+1,0,j) = GA(1)*y(i,0,j)
         z(i+1,0,j) = GA(2)*z(i,0,j)
         x(i+1,0,j) = x(i+1,0,j) + i*pp*x(i-1,0,j)
         y(i+1,0,j) = y(i+1,0,j) + i*pp*y(i-1,0,j)
         z(i+1,0,j) = z(i+1,0,j) + i*pp*z(i-1,0,j)
         x(i+1,0,j) = x(i+1,0,j) + j*pp*x(i,0,j-1)
         y(i+1,0,j) = y(i+1,0,j) + j*pp*y(i,0,j-1)
         z(i+1,0,j) = z(i+1,0,j) + j*pp*z(i,0,j-1)
      end do
      do k = 1, angularMoment3
         x(i+1,k,0) = GA(0)*x(i,k,0)
         y(i+1,k,0) = GA(1)*y(i,k,0)
         z(i+1,k,0) = GA(2)*z(i,k,0)
         x(i+1,k,0) = x(i+1,k,0) + i*pp*x(i-1,k,0)
         y(i+1,k,0) = y(i+1,k,0) + i*pp*y(i-1,k,0)
         z(i+1,k,0) = z(i+1,k,0) + i*pp*z(i-1,k,0)
         x(i+1,k,0) = x(i+1,k,0) + k*pp*x(i,k-1,0)
         y(i+1,k,0) = y(i+1,k,0) + k*pp*y(i,k-1,0)
         z(i+1,k,0) = z(i+1,k,0) + k*pp*z(i,k-1,0)
      end do
    end do
   !! End - Upward recursion in a for all b and c = 0


    do j = 1, angularMoment2
      do k = 1, angularMoment3
         x(1,k,j) = GA(0)*x(0,k,j)
         y(1,k,j) = GA(1)*y(0,k,j)
         z(1,k,j) = GA(2)*z(0,k,j)
         x(1,k,j) = x(1,k,j) + j*pp*x(0,k,j-1)
         y(1,k,j) = y(1,k,j) + j*pp*y(0,k,j-1)
         z(1,k,j) = z(1,k,j) + j*pp*z(0,k,j-1)
         x(1,k,j) = x(1,k,j) + k*pp*x(0,k-1,j)
         y(1,k,j) = y(1,k,j) + k*pp*y(0,k-1,j)
         z(1,k,j) = z(1,k,j) + k*pp*z(0,k-1,j)
      end do
    end do
    !!  Begin - Bring everything together
    do i = 1, angularMoment1 - 1
      do j = 1, angularMoment2
        do k = 1, angularMoment3
           x(i+1,k,j) = GA(0)*x(i,k,j)
           y(i+1,k,j) = GA(1)*y(i,k,j)
           z(i+1,k,j) = GA(2)*z(i,k,j)
           x(i+1,k,j) = x(i+1,k,j) + i*pp*x(i-1,k,j)
           y(i+1,k,j) = y(i+1,k,j) + i*pp*y(i-1,k,j)
           z(i+1,k,j) = z(i+1,k,j) + i*pp*z(i-1,k,j)
           x(i+1,k,j) = x(i+1,k,j) + j*pp*x(i,k,j-1)
           y(i+1,k,j) = y(i+1,k,j) + j*pp*y(i,k,j-1)
           z(i+1,k,j) = z(i+1,k,j) + j*pp*z(i,k,j-1)
           x(i+1,k,j) = x(i+1,k,j) + k*pp*x(i,k-1,j)
           y(i+1,k,j) = y(i+1,k,j) + k*pp*y(i,k-1,j)
           z(i+1,k,j) = z(i+1,k,j) + k*pp*z(i,k-1,j)
        end do
      end do
    end do

  end subroutine ThreeCOverlapIntegrals_obaraSaikaRecursion

  !>
  !!@brief  Handle exceptions
  subroutine ThreeCOverlapIntegrals_exception( typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine ThreeCOverlapIntegrals_exception

end module ThreeCOverlapIntegrals_
