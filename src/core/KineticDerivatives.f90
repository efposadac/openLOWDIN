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
!! @brief Kinetic Derivative Module.
!!        This module contains all basic functions of the kinetic derivatives calculations
!! @author  R. Gonzalez
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-02-28
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-02-28 </tt>: Ronald Gonzalez( rogonzalez@unal.edu.co )
!!        -# Basics functions has been created
!!   - <tt> 2015-02-27 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs,
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module KineticDerivatives_
  use ContractedGaussian_
  use MolecularSystem_
  use Math_
  implicit none

  public :: &
       KineticDerivatives_getDerive

  private :: &
       KineticDerivatives_obaraSaikaRecursion, &
       KineticDerivatives_int

contains

  subroutine KineticDerivatives_getDerive(this, indexA, indexB, deriveValue, specieID)
    type(ContractedGaussian), intent(in):: this(:)
    integer, intent(in) :: indexA, indexB, specieID
    !        integer, intent(in) :: cartA, cartB
    !        integer, optional :: nuclei
    !        integer, optional :: component
    !        real(8) :: output
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
    real(8) :: Jx, Jy, Jz
    real(8) :: lambda
    ! real(8) :: commonA, commonB, commonC

    angularMomentA = this(indexA)%angularMoment
    angularMomentB = this(indexB)%angularMoment
    lengthA = this(indexA)%length
    lengthB = this(indexB)%length

    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)

    lambda = MolecularSystem_getLambda(specieID)

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

    allocate(x(0:maxAngularMoment+3, 0:maxAngularMoment+3))
    allocate(y(0:maxAngularMoment+3, 0:maxAngularMoment+3))
    allocate(z(0:maxAngularMoment+3, 0:maxAngularMoment+3))

    !        write(*,"(A,I,I)")"Numero de cartesianos: ", this(indexA)%numCartesianOrbital, this(indexB)%numCartesianOrbital
    !        write(*,"(A,I)")"Tamanio de derive: ", 6*ssize
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

          ! commonA = gammaInv
          ! commonB = auxExponentA
          ! commonC = auxExponentB
          commonPreFactor = exp(-auxExponentA * auxExponentB * AB2 * gammaInv)
          commonPreFactor = commonPreFactor * sqrt(Math_PI*gammaInv)
          commonPreFactor = commonPreFactor * Math_PI * gammaInv * c1 * c2
          ! commonPreFactor = commonPreFactor * auxCoefficientA * auxCoefficientB
          ! commonPreFactor = commonPreFactor * auxConstantA * auxConstantB

          call KineticDerivatives_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+2, angularMomentB+2)

          ao12 = 0.0_8
         ! write(*,"(A)") "----------------------------------------"
         ! write(*,"(A)") "   Valores previos de las derivadas     "
         ! write(*,"(A,f)") "Common prefactor: ", commonPreFactor
         ! write(*,"(A,f)") "Common A: ", commonA
         ! write(*,"(A,f)") "Common B: ", commonB
         ! write(*,"(A,f)") "Common C: ", commonC
         ! write(*,"(A)") "----------------------------------------"
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
                      Jx = 0.0_8
                      Jy = 0.0_8
                      Jz = 0.0_8
                      ! ! x on center i
                      Ix = Ix + lambda*auxExponentA*KineticDerivatives_int(x, y, z, auxExponentA, l1+1, m1, n1, auxExponentB, l2, m2, n2)*commonPreFactor
                      if(l1>0) then
                         Ix = Ix - l1 * KineticDerivatives_int(x, y, z, auxExponentA, l1-1, m1, n1, auxExponentB, l2, m2, n2) * commonPreFactor
                      end if
                      ! write(*,"(A,f12.8)") "Ix: ", Ix
                      ! y on center i
                      Iy = Iy +  lambda*auxExponentA*KineticDerivatives_int(x, y, z, auxExponentA, l1, m1+1, n1, auxExponentB, l2, m2, n2)*commonPreFactor
                      if (m1>0) then
                         Iy = Iy - m1 * KineticDerivatives_int(x, y, z, auxExponentA, l1, m1-1, n1, auxExponentB, l2, m2, n2) * commonPreFactor
                      end if
                      ! write(*,"(A,f12.8)") "Iy: ", Iy
                      ! z on center i
                      Iz = Iz +  lambda * auxExponentA * KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1+1, auxExponentB, l2, m2, n2) * commonPreFactor
                      if (n1>0) then
                         Iz = Iz - n1 * KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1-1, auxExponentB, l2, m2, n2) * commonPreFactor
                      end if
                      ! ! write(*,"(A,f12.8)") "Iz: ", Iz
!---------------------------------------------------------------------------------------------
! Esto es nuevo buscando corregir la cinetica
!---------------------------------------------------------------------------------------------
                      ! Ix = Ix + lambda*auxExponentA*KineticDerivatives_int(x, y, z, auxExponentA, l1+1, m1, n1, auxExponentB, l2, m2, n2)*commonPreFactor
                      ! if(l1>0) then
                      !    Ix = Ix - KineticDerivatives_int(x, y, z, auxExponentA, l1-1, m1, n1, auxExponentB, l2, m2, n2) * commonPreFactor
                      ! end if
                      ! ! write(*,"(A,f12.8)") "Ix: ", Ix
                      ! ! y on center i
                      ! Iy = Iy +  lambda*auxExponentA*KineticDerivatives_int(x, y, z, auxExponentA, l1, m1+1, n1, auxExponentB, l2, m2, n2)*commonPreFactor
                      ! if (m1>0) then
                      !    Iy = Iy - KineticDerivatives_int(x, y, z, auxExponentA, l1, m1-1, n1, auxExponentB, l2, m2, n2) * commonPreFactor
                      ! end if
                      ! ! write(*,"(A,f12.8)") "Iy: ", Iy
                      ! ! z on center i
                      ! Iz = Iz +  lambda*auxExponentA * KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1+1, auxExponentB, l2, m2, n2) * commonPreFactor
                      ! if (n1>0) then
                      !    Iz = Iz - KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1-1, auxExponentB, l2, m2, n2) * commonPreFactor
                      ! end if
                      ! ! write(*,"(A,f12.8)") "Iz: ", Iz

                      ! ! x on center j
                      ! Jx = Jx + lambda*auxExponentB*KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2+1, m2, n2)*commonPreFactor
                      ! if(l2>0) then
                      !    Jx = Jx - KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2-1, m2, n2) * commonPreFactor
                      ! end if
                      ! ! write(*,"(A,f12.8)") "Jx: ", Jx
                      ! ! y on center j
                      ! Jy = Jy +  lambda*auxExponentB*KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2, m2+1, n2)*commonPreFactor
                      ! if (m2>0) then
                      !    Jy = Jy - KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2, m2-1, n2) * commonPreFactor
                      ! end if
                      ! ! write(*,"(A,f12.8)") "Jy: ", Jy
                      ! ! z on center j
                      ! Jz = Jz +  lambda*auxExponentB * KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2, m2, n2+1) * commonPreFactor
                      ! if (n2>0) then
                      !    Jz = Jz - KineticDerivatives_int(x, y, z, auxExponentA, l1, m1, n1, auxExponentB, l2, m2, n2-1) * commonPreFactor
                      ! end if
                      ! write(*,"(A,f12.8)") "Jz: ", Jz


                      ! write(*,"(A)") "----------------------------------------"
                      ! deriveValue(center_i + (0*ssize) + ao12) = deriveValue(center_i + (0*ssize) + ao12) + Ix
                      ! deriveValue(center_j + (0*ssize) + ao12) = deriveValue(center_j + (0*ssize) + ao12) + Jx

                      ! deriveValue(center_i + (1*ssize) + ao12) = deriveValue(center_i + (1*ssize) + ao12) + Iy
                      ! deriveValue(center_j + (1*ssize) + ao12) = deriveValue(center_j + (1*ssize) + ao12) + Jy

                      ! deriveValue(center_i + (2*ssize) + ao12) = deriveValue(center_i + (2*ssize) + ao12) + Iz
                      ! deriveValue(center_j + (2*ssize) + ao12) = deriveValue(center_j + (2*ssize) + ao12) + Jz
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

                      deriveValue(center_i + (0*ssize) + ao12) = deriveValue(center_i + (0*ssize) + ao12) + Ix
                      deriveValue(center_j + (0*ssize) + ao12) = deriveValue(center_j + (0*ssize) + ao12) - Ix

                      deriveValue(center_i + (1*ssize) + ao12) = deriveValue(center_i + (1*ssize) + ao12) + Iy
                      deriveValue(center_j + (1*ssize) + ao12) = deriveValue(center_j + (1*ssize) + ao12) - Iy

                      deriveValue(center_i + (2*ssize) + ao12) = deriveValue(center_i + (2*ssize) + ao12) + Iz
                      deriveValue(center_j + (2*ssize) + ao12) = deriveValue(center_j + (2*ssize) + ao12) - Iz


!                      write(*,"(f17.12)") deriveValue(center_i + ao12 )
!                      write(*,"(f17.12)") deriveValue(center_i + ssize + ao12)
!                      write(*,"(f17.12)") deriveValue(center_i+(2*ssize)+ao12)
!                      write(*,"(f17.12)") deriveValue(center_j + ao12 )
!                      write(*,"(f17.12)") deriveValue(center_j + ssize + ao12)
!                      write(*,"(f17.12)") deriveValue(center_j+(2*ssize)+ao12)
                      ao12 = ao12 + 1

                   end do
                end do
             end do
          end do

       end do
    end do
    ! write(*,"(A)") "----------------------------------------"
    ! write(*,"(A)") "Derivadas de cinetica"
    ! do i = 0, 6*ssize - 1
    !    write(*,"(f12.8)") deriveValue(i)
    ! end do

    deallocate(x)
    deallocate(y)
    deallocate(z)
    
  end subroutine KineticDerivatives_getDerive

  !>
  !! @brief Implementacion de la recursion propuesta por Obara - Saika para el caso de integrales de energia cinetica
  !! @author E. F. Posada, 2010
  !! @version 1.0
  subroutine KineticDerivatives_obaraSaikaRecursion(x, y, z, PA, PB, zeta, angularMomentIndexA, angularMomentIndexB)
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

  end subroutine KineticDerivatives_obaraSaikaRecursion

  !>
  !! @author R. Gonzalez, rogonzalez@unal.edu.co   2014
  function KineticDerivatives_int( x, y, z,  auxExponentA, l1, m1, n1, auxExponentB, l2, m2, n2) result (output)
    implicit none

    integer, intent(in) :: l1, m1, n1
    integer, intent(in) :: l2, m2, n2
    real(8), intent(in) :: auxExponentA, auxExponentB
    real(8), intent(in), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: output

    real(8) :: I1, I2, I3, I4
    real(8) :: Ix, Iy, Iz


    if (l1 == 0 .or. l2 == 0) then
       I1 = 0.0_8
    else
       I1 = x(l1-1,l2-1) * y(m1,m2) * z(n1,n2)
    end if

    I2 = x(l1+1,l2+1) * y(m1,m2) * z(n1,n2)

    if (l2 == 0) then
       I3 = 0.0_8
    else
       I3 = x(l1+1,l2-1) * y(m1,m2) * z(n1,n2)
    end if

    if (l1 == 0) then
       I4 = 0.0_8
    else
       I4 = x(l1-1,l2+1) * y(m1,m2) * z(n1,n2)
    end if

    Ix = 0.5 * l1 * l2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * l2 * I3 - l1 * auxExponentB * I4

    if (m1 == 0 .or. m2 == 0) then
       I1 = 0.0_8
    else
       I1 =   x(l1,l2) * y(m1-1,m2-1) * z(n1,n2)
    end if

    I2 = x(l1,l2) * y(m1+1,m2+1) * z(n1,n2)

    if (m2 == 0) then
       I3 = 0.0_8
    else
       I3 = x(l1,l2) * y(m1+1,m2-1) * z(n1,n2)
    end if

    if (m1 == 0) then
       I4 = 0.0_8
    else
       I4 = x(l1,l2) * y(m1-1,m2+1) * z(n1,n2)
    end if

    Iy = 0.5 * m1 * m2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * m2 * I3 - m1 * auxExponentB * I4

    if (n1 == 0 .or. n2 == 0) then
       I1 = 0.0_8
    else
       I1 = x(l1,l2) * y(m1,m2) * z(n1-1,n2-1)
    end if

    I2 = x(l1,l2) * y(m1,m2) * z(n1+1,n2+1)

    if (n2 == 0) then
       I3 = 0.0_8
    else
       I3 = x(l1,l2) * y(m1,m2) * z(n1+1,n2-1)
    end if

    if (n1 == 0) then
       I4 = 0.0_8
    else
       I4 = x(l1,l2) * y(m1,m2) * z(n1-1,n2+1)
    end if

    Iz = 0.5 * n1 * n2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * n2 * I3 - n1 * auxExponentB * I4
    
    output = Ix + Iy + Iz
   ! write(*,"(A,f)") "Integral: ", output

  end function KineticDerivatives_int

end module KineticDerivatives_
