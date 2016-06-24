!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package
!!
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Overlap Derivative Module.
!!        This module contains all basic functions of the overlap derivatives calculations
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-02-28
!!
!! <b> History: </b>
!!
!!   - <tt> 2015-03-13 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs,
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module OverlapDerivatives_
  use ContractedGaussian_
  use MolecularSystem_
  use Math_
  implicit none

  public :: &
       OverlapDerivatives_getDerive

  private :: &
       OverlapDerivatives_obaraSaikaRecursion

contains

  subroutine OverlapDerivatives_getDerive(this, indexA, indexB, deriveValue, specieID)
    type(ContractedGaussian), intent(in):: this(:)
    integer, intent(in) :: indexA, indexB, specieID
    integer :: i

    integer :: angularMomentA, angularMomentB
    integer :: lengthA, lengthB
    real(8) :: originA(0:3)
    real(8) :: originB(0:3)
    integer :: ssize
    integer :: center_i, center_j
    real(8) :: AB2
    real(8), allocatable :: deriveValue(:)
    real(8), allocatable :: x(:,:), y(:,:), z(:,:)
    integer :: maxAngularMoment
    integer :: p1, p2
    real(8) :: auxExponentA, auxCoefficientA, auxContConstantA, auxPrimConstantA, c1
    real(8) :: auxExponentB, auxCoefficientB, auxContConstantB, auxPrimConstantB, c2
    real(8) :: gamma, gammaInv
    real(8) :: PA(0:3), PB(0:3), P(0:3)
    real(8) :: commonPreFactor
    integer :: ao12
    integer :: ii, jj, kk, ll
    integer :: l1, m1, n1
    integer :: l2, m2, n2
    real(8) :: Ix, Iy, Iz
    real(8) :: lambda

    angularMomentA = this(indexA)%angularMoment
    angularMomentB = this(indexB)%angularMoment
    lengthA = this(indexA)%length
    lengthB = this(indexB)%length

    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    lambda = 2.0_8!MolecularSystem_getLambda(specieID)

    do i = 1, 3
       originA(i-1) = this(indexA)%origin(i)
       originB(i-1) = this(indexB)%origin(i)
    end do

    ssize = this(indexA)%numCartesianOrbital*this(indexB)%numCartesianOrbital
    center_i = 0
    center_j = 3*ssize

    do i = 0, 2
       AB2 = AB2 + (originA(i) - originB(i))*(originA(i) - originB(i))
    end do

    if(allocated(deriveValue)) deallocate(deriveValue)
    allocate(deriveValue(0:6*ssize))

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(y(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(z(0:maxAngularMoment+2, 0:maxAngularMoment+2))

    deriveValue = 0.0_8

    do p1 = 0, lengthA - 1
       auxExponentA = this(indexA)%orbitalExponents(p1 + 1)
       auxCoefficientA = this(indexA)%contractionCoefficients(p1 + 1)
       auxPrimConstantA = this(indexA)%primNormalization(p1 + 1,1)
       auxContConstantA = this(indexA)%contNormalization(1)
       c1 = auxCoefficientA*auxPrimConstantA*auxContConstantA
       do p2 = 0, lengthB - 1
          auxExponentB = this(indexB)%orbitalExponents(p2 + 1)
          auxCoefficientB = this(indexB)%contractionCoefficients(p2 + 1)
          auxPrimConstantB = this(indexB)%primNormalization(p2 + 1,1)
          auxContConstantB = this(indexB)%contNormalization(1)
          c2 = auxCoefficientB*auxPrimConstantB*auxContConstantB

          gamma = auxExponentA + auxExponentB
          gammaInv = 1.0/gamma

          P(0)  = (auxExponentA*originA(0) + auxExponentB*originB(0))*gammaInv
          P(1)  = (auxExponentA*originA(1) + auxExponentB*originB(1))*gammaInv
          P(2)  = (auxExponentA*originA(2) + auxExponentB*originB(2))*gammaInv
          PA(0) = P(0) - originA(0)
          PA(1) = P(1) - originA(1)
          PA(2) = P(2) - originA(2)
          PB(0) = P(0) - originB(0)
          PB(1) = P(1) - originB(1)
          PB(2) = P(2) - originB(2)

          commonPreFactor = exp(-auxExponentA * auxExponentB * AB2 * gammaInv)
          commonPreFactor = commonPreFactor * sqrt(Math_PI*gammaInv)
          commonPreFactor = commonPreFactor * Math_PI * gammaInv * c1 * c2

          ! write(*,"(A,f17.12, A, I, A, I)") "Common Prefactor: ", commonPreFactor, " am1: ", angularMomentA, " am2: ", angularMomentB

          call OverlapDerivatives_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+1, angularMomentB+1)

          ao12 = 0.0_8
          do ii = 0, angularMomentA
             l1 = angularMomentA - ii
             do jj = 0, ii
                m1 = ii - jj
                n1 = jj
                do kk = 0, angularMomentB
                   l2 = angularMomentB - kk
                   do ll = 0, kk
                      m2 = kk - ll
                      n2 = ll

                      Ix = 0.0_8
                      Iy = 0.0_8
                      Iz = 0.0_8

                      ! x on center i
                      Ix = Ix + lambda*auxExponentA*commonPreFactor*x(l1+1,l2)*y(m1,m2)*z(n1,n2)
                      ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1+1,l2), " y: ", y(m1,m2), " z: ", z(n1,n2), " Ix: ", Ix
                      if (l1>0) then
                         Ix = Ix - l1*commonPreFactor*x(l1-1,l2)*y(m1,m2)*z(n1,n2)
                         ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1-1,l2), "y: ", y(m1,m2), "z: ", z(n1,n2), " Ix: ", Ix
                      end if
                      
                      ! y on center i 
                      Iy = Iy + lambda*auxExponentA*commonPreFactor*x(l1,l2)*y(m1+1,m2)*z(n1,n2)
                      ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1,l2), "y: ", y(m1+1,m2), "z: ", z(n1,n2), " Iy: ", Iy
                      if (m1>0) then
                         Iy = Iy - m1*commonPreFactor*x(l1,l2)*y(m1-1,m2)*z(n1,n2)
                         ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1,l2), "y: ", y(m1-1,m2), "z: ", z(n1,n2), " Iy: ", Iy
                      end if
                      
                      ! z on center i 
                      Iz = Iz + lambda*auxExponentA*commonPreFactor*x(l1,l2)*y(m1,m2)*z(n1+1,n2)
                      ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1,l2), "y: ", y(m1,m2), "z: ", z(n1+1,n2), " Iz: ", Iz
                      if (n1>0) then
                         Iz = Iz - n1*commonPreFactor*x(l1,l2)*y(m1,m2)*z(n1-1,n2)
                         ! write(*,"(A,I17.12,4(A,f17.12))") "am1: ", auxExponentA, "x: ", x(l1,l2), "y: ", y(m1,m2), "z: ", z(n1-1,n2), " Iz: ", Iz
                      end if
                      ! write(*,"(A)") "-----------------------------------------------"
                      
                      deriveValue(center_i + (0*ssize) + ao12) = deriveValue(center_i + (0*ssize) + ao12) + Ix
                      deriveValue(center_j + (0*ssize) + ao12) = deriveValue(center_j + (0*ssize) + ao12) - Ix

                      deriveValue(center_i + (1*ssize) + ao12) = deriveValue(center_i + (1*ssize) + ao12) + Iy
                      deriveValue(center_j + (1*ssize) + ao12) = deriveValue(center_j + (1*ssize) + ao12) - Iy

                      deriveValue(center_i + (2*ssize) + ao12) = deriveValue(center_i + (2*ssize) + ao12) + Iz
                      deriveValue(center_j + (2*ssize) + ao12) = deriveValue(center_j + (2*ssize) + ao12) - Iz


                      ao12 = ao12 + 1

                   end do
                end do
             end do
          end do

       end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)    
    ! write(*,"(A)") "Derivadas de overlap"
   ! do i = 0, 6*ssize - 1
   !    write(*,"(f17.12)") deriveValue(i)
   ! end do

  end subroutine OverlapDerivatives_getDerive

  !>
  !! @brief Implementacion de la recursion propuesta por Obara - Saika para el caso de integrales de energia cinetica
  !! @author E. F. Posada, 2010
  !! @version 1.0
  subroutine OverlapDerivatives_obaraSaikaRecursion(x, y, z, PA, PB, zeta, angularMomentIndexA, angularMomentIndexB)
    implicit none

    real(8), intent(inout), allocatable :: x(:,:), y(:,:), z(:,:)
    real(8), intent(in) :: PA(0:3), PB(0:3)
    real(8), intent(in) :: zeta
    integer, intent(in) :: angularMomentIndexA, angularMomentIndexB

    real(8) :: twoZetaInv
    integer :: i, j, k

    twoZetaInv = 1_8/(2_8*zeta)

    x = 0.0_8
    y = 0.0_8
    z = 0.0_8

    x(0,0) = 1.0_8
    y(0,0) = 1.0_8
    z(0,0) = 1.0_8

    !! Upward recursion in j for i=0
    x(0,1) = PB(0)
    y(0,1) = PB(1)
    z(0,1) = PB(2)

    do j=1, angularMomentIndexB -1
       x(0,j+1) = PB(0)*x(0,j)
       y(0,j+1) = PB(1)*y(0,j)
       z(0,j+1) = PB(2)*z(0,j)
       x(0,j+1) = x(0,j+1) + j*twoZetaInv*x(0,j-1)
       y(0,j+1) = y(0,j+1) + j*twoZetaInv*y(0,j-1)
       z(0,j+1) = z(0,j+1) + j*twoZetaInv*z(0,j-1)
    end do

    !! Upward recursion in i for all j's
    x(1,0) = PA(0)
    y(1,0) = PA(1)
    z(1,0) = PA(2)

    do j=1, angularMomentIndexB
       x(1,j) = PA(0)*x(0,j)
       y(1,j) = PA(1)*y(0,j)
       z(1,j) = PA(2)*z(0,j)
       x(1,j) = x(1,j) + j*twoZetaInv*x(0,j-1)
       y(1,j) = y(1,j) + j*twoZetaInv*y(0,j-1)
       z(1,j) = z(1,j) + j*twoZetaInv*z(0,j-1)
    end do

    do i=1, angularMomentIndexA - 1
       x(i+1,0) = PA(0)*x(i,0)
       y(i+1,0) = PA(1)*y(i,0)
       z(i+1,0) = PA(2)*z(i,0)
       x(i+1,0) = x(i+1,0) + i*twoZetaInv*x(i-1,0)
       y(i+1,0) = y(i+1,0) + i*twoZetaInv*y(i-1,0)
       z(i+1,0) = z(i+1,0) + i*twoZetaInv*z(i-1,0)
       do j=1, angularMomentIndexB
          x(i+1,j) = PA(0)*x(i,j)
          y(i+1,j) = PA(1)*y(i,j)
          z(i+1,j) = PA(2)*z(i,j)
          x(i+1,j) = x(i+1,j) + i*twoZetaInv*x(i-1,j)
          y(i+1,j) = y(i+1,j) + i*twoZetaInv*y(i-1,j)
          z(i+1,j) = z(i+1,j) + i*twoZetaInv*z(i-1,j)
          x(i+1,j) = x(i+1,j) + j*twoZetaInv*x(i,j-1)
          y(i+1,j) = y(i+1,j) + j*twoZetaInv*y(i,j-1)
          z(i+1,j) = z(i+1,j) + j*twoZetaInv*z(i,j-1)
       end do
    end do

  end subroutine OverlapDerivatives_obaraSaikaRecursion

end module OverlapDerivatives_
