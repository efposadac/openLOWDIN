!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!> @brief This module contains all the routines to Confininig potentials for using in DFTB pseudoatom generation
!! @author I. Ortiz-Verano (ieortizv@unal.edu.co)
!! @version 1.0
!! <b> Creation data : </b> 14-02-2018
!!
!! <b> History change: </b>
!!
!!   - <tt> 14-02-18 </tt>:  I. Ortiz-Verano (ieortizv@unal.edu.co)
!!        -# Creacioon del modulo y metodos basicos

module ConfiningIntegrals_
  use Exception_
  use Math_
  use ContractedGaussian_
  !  use AtomicElement_
  use ConfiningPotential_
  implicit none

contains
  subroutine ConfiningIntegrals_computeShell(contractedGaussianA,contractedGaussianB,integral)
    !      subroutine ConfiningIntegrals_computeShell(contractedGaussianA,contractedGaussianB,integral,covalentRadius) 
    implicit none

    type(contractedGaussian), intent(in) :: contractedGaussianA,contractedGaussianB
    real(8), intent(out) :: integral(contractedGaussianA%numCartesianOrbital * contractedGaussianB%numCartesianOrbital ) 
    !    real(8), intent(in) :: covalentRadius
    integer ::  am1(0:3)
    integer ::  am2(0:3)
    integer ::  nprim1
    integer ::  nprim2
    real(8) ::  A(0:3)
    real(8) ::  B(0:3)
    real(8) ::  exp1(0:contractedGaussianA%length)
    real(8) ::  exp2(0:contractedGaussianB%length)
    real(8) ::  coef1(0:contractedGaussianA%length)
    real(8) ::  coef2(0:contractedGaussianB%length)
    real(8) ::  nor1(0:contractedGaussianA%length)
    real(8) ::  nor2(0:contractedGaussianB%length)
    real(8) :: auxIntegral
    real(8) :: radiusCoefficient
    real(8) :: exponentM
    real(8) :: confPot
    integer, allocatable :: angularMomentIndexA(:,:)
    integer, allocatable :: angularMomentIndexB(:,:)
    integer ::  i, m, p, q

    integral = 0.0_8

    if(allocated(angularMomentIndexA)) deallocate(angularMomentIndexA)
    if(allocated(angularMomentIndexB)) deallocate(angularMomentIndexB)

    allocate(angularMomentIndexA(3, contractedGaussianA%numCartesianOrbital))
    allocate(angularMomentIndexB(3, contractedGaussianB%numCartesianOrbital))

    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexA, contractedGaussianA)
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexB, contractedGaussianB)

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

    m = 0

    do p = 1, contractedGaussianA%numcartesianOrbital
       do q = 1, contractedGaussianB%numcartesianOrbital

          m = m + 1

          exp1(0:nprim1-1) = contractedGaussianA%orbitalExponents(1:nprim1)
          nor1(0:nprim1-1) = contractedGaussianA%primNormalization(1:nprim1,p)

          exp2(0:nprim2-1) = contractedGaussianB%orbitalExponents(1:nprim2)
          nor2(0:nprim2-1) = contractedGaussianB%primNormalization(1:nprim2,q)

          am1 = 0
          am2 = 0

          am1(0:2) = angularMomentIndexA(1:3, p)
          am2(0:2) = angularMomentIndexB(1:3, q)

          call ConfiningIntegrals_computeContractionPair(am1, am2, nprim1, nprim2, A, B, exp1, exp2, coef1, coef2, nor1, nor2, auxIntegral)

          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q)

          integral(m) = auxIntegral

       end do
    end do
    !! DEBUG
    ! call ConfiningPotential_show() 

  end subroutine ConfiningIntegrals_computeShell


  subroutine ConfiningIntegrals_computeContractionPair(angularMomentIndexA, angularMomentIndexB, lengthA, lengthB, A, B, &
       orbitalExponentsA, orbitalExponentsB, &
       contractionCoefficientsA, contractionCoefficientsB, &
       normalizationConstantsA, normalizationConstantsB, integralValue)
    implicit none

    integer, intent(in) :: angularMomentIndexA(0:3), angularMomentIndexB(0:3)
    integer, intent(in) :: lengthA, lengthB
    real(8), intent(in) :: A(0:3), B(0:3)
    real(8), intent(in) :: orbitalExponentsA(0:lengthA), orbitalExponentsB(0:lengthB)
    real(8), intent(in) :: contractionCoefficientsA(0:lengthA), contractionCoefficientsB(0:lengthB)
    real(8), intent(in) :: normalizationConstantsA(0:lengthA), normalizationConstantsB(0:lengthB)
    real(8), intent(out):: integralValue

    real(8), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: AB2  ! cuadrado de la distancia
    real(8) :: auxExponentA, auxCoefficientA, auxConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxConstantB
    real(8) :: gamma, gammaInv
    real(8) :: PA(0:3), PB(0:3), P(0:3)
    real(8) :: commonPreFactor
    real(8) :: x0, y0, z0

    integer :: angularMomentA, angularMomentB
    integer :: maxAngularMoment
    integer :: p1, p2 !< iteradores

    angularMomentA = sum(angularMomentIndexA)
    angularMomentB = sum(angularMomentIndexB)

    integralValue = 0.0_8

    maxAngularMoment = max(angularMomentA, angularMomentB) + 1

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(y(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(z(0:maxAngularMoment+2, 0:maxAngularMoment+2))

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
          gamma = auxExponentA + auxExponentB ! \xi, \nu y \rho de (3) en [Ahlrichs2006]
          gammaInv = 1.0/gamma

          P(0) = (auxExponentA*A(0) + auxExponentB*B(0))*gammaInv ! P y Q de (4) en [Ahlrichs2006]
          P(1) = (auxExponentA*A(1) + auxExponentB*B(1))*gammaInv
          P(2) = (auxExponentA*A(2) + auxExponentB*B(2))*gammaInv
          PA(0) = P(0) - A(0) ! P-Q de (5) en [Ahlrichs2006] (??????????????????)
          PA(1) = P(1) - A(1)
          PA(2) = P(2) - A(2)
          PB(0) = P(0) - B(0)
          PB(1) = P(1) - B(1)
          PB(2) = P(2) - B(2)

          commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*gammaInv) & ! Esta primera línea es (6) en [Ahlrichs2006]
               * sqrt(Math_PI*gammaInv) * Math_PI * gammaInv & ! Fracción entre paréntesis en (7) en [Ahlrichs2006]
               * auxCoefficientA * auxCoefficientB * auxConstantA * auxConstantB

          !! recursion
          call ConfiningIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+2, angularMomentB+2)

          x0 = x(angularMomentIndexA(0),angularMomentIndexB(0))
          y0 = y(angularMomentIndexA(1),angularMomentIndexB(1))
          z0 = z(angularMomentIndexA(2),angularMomentIndexB(2))

          integralValue = integralValue + commonPreFactor*x0*y0*z0

       end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)

  end subroutine ConfiningIntegrals_computeContractionPair



  !>
  !!@brief Evaluates Overlap integral for a couple of gaussian primitives
  !!@author E. F. Posada, 2013
  !!@version 1.0
  !!@return Returns the integral value 
  !!@param A, B : origin A and B
  subroutine ConfiningIntegrals_computePrimitivePair(angularMomentA, angularMomentB, &
       numcartesianOrbitalA, numcartesianOrbitalB, & 
       A, B, &
       contractionCoefficientsA, contractionCoefficientsB, &
       orbitalExponentsA, orbitalExponentsB, &
       normalizationConstantsA, normalizationConstantsB, integralValue)
    implicit none

    integer, intent(in) :: angularMomentA, angularMomentB
    integer, intent(in) :: numcartesianOrbitalA, numcartesianOrbitalB
    real(8), intent(in) :: A(3), B(3)
    real(8), intent(in) :: orbitalExponentsA, orbitalExponentsB
    real(8), intent(in) :: contractionCoefficientsA, contractionCoefficientsB
    real(8), intent(in) :: normalizationConstantsA(numcartesianOrbitalA)
    real(8), intent(in) :: normalizationConstantsB(numcartesianOrbitalB)
    real(8), intent(out):: integralValue(numcartesianOrbitalA*numcartesianOrbitalB)

    real(8), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: AB2
    real(8) :: gamma, gammaInv
    real(8) :: PA(0:3), PB(0:3), P(0:3)
    real(8) :: commonPreFactor
    real(8) :: x0, y0, z0

    integer, allocatable :: angularMomentIndexA(:,:)
    integer, allocatable :: angularMomentIndexB(:,:)

    integer :: maxAngularMoment
    integer :: p1, p2 !< iteradores
    integer :: counter

    integralValue = 0.0_8

    maxAngularMoment = max(angularMomentA, angularMomentB) + 1

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(y(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(z(0:maxAngularMoment+2, 0:maxAngularMoment+2))

    x = 0.0_8
    y = 0.0_8
    z = 0.0_8

    AB2 = 0.0_8
    AB2 = AB2 + (A(1) - B(1)) * (A(1) - B(1))
    AB2 = AB2 + (A(2) - B(2)) * (A(2) - B(2))
    AB2 = AB2 + (A(3) - B(3)) * (A(3) - B(3))

    gamma = orbitalExponentsA + orbitalExponentsB
    gammaInv = 1.0/gamma

    P(0) = (orbitalExponentsA*A(1) + orbitalExponentsB*B(1))*gammaInv
    P(1) = (orbitalExponentsA*A(2) + orbitalExponentsB*B(2))*gammaInv
    P(2) = (orbitalExponentsA*A(3) + orbitalExponentsB*B(3))*gammaInv
    PA(0) = P(0) - A(1)
    PA(1) = P(1) - A(2)
    PA(2) = P(2) - A(3)
    PB(0) = P(0) - B(1)
    PB(1) = P(1) - B(2)
    PB(2) = P(2) - B(3)

    !! recursion
    call ConfiningIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+2, angularMomentB+2)

    if(allocated(angularMomentIndexA)) deallocate(angularMomentIndexA)
    if(allocated(angularMomentIndexB)) deallocate(angularMomentIndexB)

    allocate(angularMomentIndexA(3, numCartesianOrbitalA))
    allocate(angularMomentIndexB(3, numCartesianOrbitalB))

    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexA, angularMoment=angularMomentA)
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexB, angularMoment=angularMomentB)

    counter = 1

    do p1 = 1, numcartesianOrbitalA
       do p2 = 1, numcartesianOrbitalB

          commonPreFactor = exp(-orbitalExponentsA*orbitalExponentsB*AB2*gammaInv) &
               * sqrt(Math_PI*gammaInv) * Math_PI * gammaInv &
               * contractionCoefficientsA * contractionCoefficientsB * &
               normalizationConstantsA(p1) * normalizationConstantsB(p2)

          x0 = x(angularMomentIndexA(1,p1),angularMomentIndexB(1,p2))
          y0 = y(angularMomentIndexA(2,p1),angularMomentIndexB(2,p2))
          z0 = z(angularMomentIndexA(3,p1),angularMomentIndexB(3,p2))

          integralValue(counter) = commonPreFactor*x0*y0*z0

          counter = counter + 1

       end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)

  end subroutine ConfiningIntegrals_computePrimitivePair

  !>
  !!@brief Implementation of recursion proposed by Obara-Saika for overlap integrals.
  !!@author Edwin Posada, 2010
  !!@return x, y, z : recursion matrixes
  !!@param PA, PB : reduced origin for gaussian A and B
  !!@param gamma : reduced exponent
  !!@see Gaussian product: if you want to know what reduced exponent and origin is.
  subroutine ConfiningIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMoment1, angularMoment2)
    implicit none

    real(8), intent(inout), allocatable :: x(:,:), y(:,:), z(:,:)
    real(8), intent(in) :: PA(0:3), PB(0:3)
    real(8), intent(in) :: gamma
    integer, intent(in) :: angularMoment1, angularMoment2

    real(8) :: pp, gammaInv
    integer :: i, j, k

    real(8) :: mExponent
    real(8) :: factor
    real(8) :: rZero
    real(8) :: coefRZero



    ! coefRZero = 2
    ! rZero = coefRZero * covalentRadius
    ! pp = 1/(2*gamma)
    ! gammaInv = 2*pp
    ! factor = (4 * Math_PI ) * sqrt(gammaInv*gammaInv*gammaInv) ! eq (47) de Ahlrichs2006
    ! print *, "Factor: ",factor,i,coefRZero,rZero!,j,k,x,y,z
    !* exp(-PA) !

    x(0,0) = 1.0_8
    y(0,0) = 1.0_8
    z(0,0) = 1.0_8

    !! Upward recursion in j for i=0
    x(0,1) = PB(0)
    y(0,1) = PB(1)
    z(0,1) = PB(2)

    do j=1, angularMoment2 -1
       x(0,j+1) = PB(0)*x(0,j)
       y(0,j+1) = PB(1)*y(0,j)
       z(0,j+1) = PB(2)*z(0,j)
       x(0,j+1) = x(0,j+1) + j*pp*x(0,j-1)
       y(0,j+1) = y(0,j+1) + j*pp*y(0,j-1)
       z(0,j+1) = z(0,j+1) + j*pp*z(0,j-1)
    end do

    !! Upward recursion in i for all j
    x(1,0) = PA(0)
    y(1,0) = PA(1)
    z(1,0) = PA(2)

    do j=1, angularMoment2
       x(1,j) = PA(0)*x(0,j)
       y(1,j) = PA(1)*y(0,j)
       z(1,j) = PA(2)*z(0,j)
       x(1,j) = x(1,j) + j*pp*x(0,j-1)
       y(1,j) = y(1,j) + j*pp*y(0,j-1)
       z(1,j) = z(1,j) + j*pp*z(0,j-1)
    end do

    do i=1, angularMoment1 - 1
       x(i+1,0) = PA(0)*x(i,0)
       y(i+1,0) = PA(1)*y(i,0)
       z(i+1,0) = PA(2)*z(i,0)
       x(i+1,0) = x(i+1,0) + i*pp*x(i-1,0)
       y(i+1,0) = y(i+1,0) + i*pp*y(i-1,0)
       z(i+1,0) = z(i+1,0) + i*pp*z(i-1,0)
       do j=1, angularMoment2
          x(i+1,j) = PA(0)*x(i,j)
          y(i+1,j) = PA(1)*y(i,j)
          z(i+1,j) = PA(2)*z(i,j)
          x(i+1,j) = x(i+1,j) + i*pp*x(i-1,j)
          y(i+1,j) = y(i+1,j) + i*pp*y(i-1,j)
          z(i+1,j) = z(i+1,j) + i*pp*z(i-1,j)
          x(i+1,j) = x(i+1,j) + j*pp*x(i,j-1)
          y(i+1,j) = y(i+1,j) + j*pp*y(i,j-1)
          z(i+1,j) = z(i+1,j) + j*pp*z(i,j-1)
       end do
    end do

  end subroutine ConfiningIntegrals_obaraSaikaRecursion




  ! subroutine ConfiningPotentialForElectrons_constructor(covalentRadious)
  !   implicit none

  ! end subroutine ConfiningPotentialForElectrons_constructor

  !>
  !!@brief  Handle exceptions
  subroutine ConfiningIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine ConfiningIntegrals_exception


end module ConfiningIntegrals_
