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
!! @brief Attraction Derivative Module.
!!        This module contains all basic functions of the potential derivatives calculations
!! @author  R. Gonzalez
!! @author  M. Rodriguez
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-02-28
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-02-28 </tt>: Ronald Gonzalez( rogonzalez@unal.edu.co ), Matheus Rodriguez
!!        -# Basics functions has been created
!!   - <tt> 2015-02-27 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs,
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module AttractionDerivatives_
  use ParticleManager_
  use MolecularSystem_
  use ContractedGaussian_
  use Math_
  implicit none

  public :: &
       AttractionDerivatives_getDerive

  private :: &
       AttractionDerivatives_obaraSaikaRecursion!, &
  !      AttractionDerivatives_int

contains

  subroutine AttractionDerivatives_getDerive(this, indexA, indexB, deriveValue, centerA, centerB, specieID)
    type(ContractedGaussian), intent(in):: this(:)
    integer, intent(in) :: indexA, indexB, centerA, centerB, specieID
    !        integer, intent(in) :: cartA, cartB
    !        integer, optional :: nuclei
    !        integer, optional :: component
    !        real(8) :: output
    integer :: i, j

    integer :: angularMomentA, angularMomentB
    integer :: lengthA, lengthB
    real(8) :: originA(0:3)
    real(8) :: originB(0:3)
    integer :: ssize
    integer :: recurSize, zSize
    integer :: center_i, center_j
    real(8) :: AB2
    real(8), allocatable :: deriveValue(:)
    real(8), allocatable :: vi(:,:,:), vx(:,:,:), vy(:,:,:), vz(:,:,:)
    real(8), allocatable :: pointCharges(:,:)
    integer :: maxAngularMoment
    integer :: ixm1, iym1, izm1
    integer :: jxm1, jym1, jzm1
    integer :: nCenters, center, Z
    integer :: p1, p2
    real(8) :: c1, c2
    real(8) :: auxExponentA, auxCoefficientA, auxPrimConstantA, auxContConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxPrimConstantB, auxContConstantB
    real(8) :: gamma, gammaInv, pfac, temp
    real(8) :: PA(0:3), PB(0:3), PC(0:3), P(0:3)
    real(8) :: commonPreFactor
    integer :: ao12
    integer :: ii, jj, kk, ll
    integer :: l1, m1, n1
    integer :: l2, m2, n2
    integer :: iind, jind

    angularMomentA = this(indexA)%angularMoment
    angularMomentB = this(indexB)%angularMoment
    lengthA = this(indexA)%length
    lengthB = this(indexB)%length

    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)

    recurSize = maxAngularMoment + 4
    recurSize = (recurSize - 1)*recurSize*(recurSize + 1) + 1
    zSize = (maxAngularMoment + 3)*2 + 1
    allocate(vi(0:recurSize, 0:recurSize, 0:zSize))
    allocate(vx(0:recurSize, 0:recurSize, 0:zSize))
    allocate(vy(0:recurSize, 0:recurSize, 0:zSize))
    allocate(vz(0:recurSize, 0:recurSize, 0:zSize))

    do i = 1, 3
       originA(i-1) = this(indexA)%origin(i)
       originB(i-1) = this(indexB)%origin(i)
    end do

    ssize = this(indexA)%numCartesianOrbital*this(indexB)%numCartesianOrbital
    center_i = (centerA - 1)*3*ssize
    center_j = (centerB - 1)*3*ssize

    izm1 = 1;
    iym1 = angularMomentA + 1 + 1;  
    ixm1 = iym1 * iym1;
    jzm1 = 1;
    jym1 = angularMomentB + 1 + 1;  
    jxm1 = jym1 * jym1;


    do i = 0, 2
       AB2 = AB2 + (originA(i) - originB(i))*(originA(i) - originB(i))
    end do

    nCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if(allocated(deriveValue)) deallocate(deriveValue)
    allocate(deriveValue(0:3*nCenters*ssize))

    if(allocated(pointCharges)) deallocate(pointCharges)
    allocate(pointCharges(0:nCenters,0:4))
    
    deriveValue = 0.0_8

    j = 0
    do i = 1, size(ParticleManager_instance)
       if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
          pointCharges(j,0) = ParticleManager_instance(i)%particlePtr%charge
          pointCharges(j,1) = ParticleManager_instance(i)%particlePtr%origin(1)
          pointCharges(j,2) = ParticleManager_instance(i)%particlePtr%origin(2)
          pointCharges(j,3) = ParticleManager_instance(i)%particlePtr%origin(3)
          j = j + 1
       end if
    end do

    ! write(*,"(A)") "------------------------------------ "
    do p1 = 0, lengthA - 1
       auxExponentA = this(indexA)%orbitalExponents(p1 + 1)
       auxCoefficientA = this(indexA)%contractionCoefficients(p1 + 1)
       auxPrimConstantA = this(indexA)%primNormalization(p1 + 1,1)
       auxContConstantA = this(indexA)%contNormalization(1)
       c1 = auxCoefficientA*auxPrimConstantA*auxContConstantA
       ! write(*,"(A)") "------------------------------------ "
       do p2 = 0, lengthB - 1
          auxExponentB = this(indexB)%orbitalExponents(p2 + 1)
          auxCoefficientB = this(indexB)%contractionCoefficients(p2 + 1)
          auxPrimConstantB = this(indexB)%primNormalization(p2 + 1,1)
          auxContConstantB = this(indexB)%contNormalization(1)
          c2 = auxCoefficientB*auxPrimConstantB*auxContConstantB
          
          !! Debug Mauricio Rodas
          ! write(*,"(A,f,A,f,A,f)") "coefA: ", auxCoefficientA, " PrimA:", auxPrimConstantA, " contA: ", auxContConstantA
          ! write(*,"(A,f,A,f,A,f)") "coefB: ", auxCoefficientB, " PrimB:", auxPrimConstantB, " contB: ", auxContConstantB
          ! write(*,"(A,I,A,I,A,f,A,f)") "P1: ", p1, " P2:", p2, " c1: ", c1, " c2: ", c2

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
          commonPreFactor = commonPreFactor * Math_PI * gammaInv*c1*c2

          do center=0, nCenters - 1
             
             Z = pointCharges(center,0)

             PC(0) = P(0) - pointCharges(center,1)
             PC(1) = P(1) - pointCharges(center,2)
             PC(2) = P(2) - pointCharges(center,3)
             call AttractionDerivatives_obaraSaikaRecursion(vi, vx, vy, vz, PA, PB, PC, gamma, angularMomentA+1, angularMomentB+1)
          
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

                         iind = l1 * ixm1 + m1 * iym1 + n1 * izm1
                         jind = l2 * jxm1 + m2 * jym1 + n2 * jzm1

                         pfac = commonPreFactor * Z
                         temp = 2.0*auxExponentA*vi(iind+ixm1,jind,0)
                         if (l1>0) then
                            temp = temp - l1*vi(iind-ixm1,jind,0)
                         end if
                         deriveValue(center_i+(0*ssize)+ao12) = deriveValue(center_i+(0*ssize)+ao12) - temp * pfac

                         temp = 2.0*auxExponentB*vi(iind,jind+jxm1,0)
                         if (l2>0) then
                            temp = temp - l2*vi(iind,jind-jxm1,0)
                         end if
                         deriveValue(center_j+(0*ssize)+ao12) = deriveValue(center_j+(0*ssize)+ao12) - temp * pfac
                         
                         deriveValue(3*ssize*center+ao12) = deriveValue(3*ssize*center+ao12) - vx(iind,jind,0) * pfac

                         temp = 2.0*auxExponentA*vi(iind+iym1,jind,0)
                         if (m1>0) then
                            temp = temp - m1*vi(iind-iym1,jind,0)
                         end if
                         deriveValue(center_i+(1*ssize)+ao12) = deriveValue(center_i+(1*ssize)+ao12) - temp * pfac
                         
                         temp = 2.0*auxExponentB*vi(iind,jind+jym1,0)
                         if (m2>0) then
                            temp = temp - m2*vi(iind,jind-jym1,0)
                         end if
                         deriveValue(center_j+(1*ssize)+ao12) = deriveValue(center_j+(1*ssize)+ao12) - temp * pfac
                         
                         deriveValue(3*ssize*center+ssize+ao12) = deriveValue(3*ssize*center+ssize+ao12) - vy(iind,jind,0) * pfac

                         temp = 2.0*auxExponentA*vi(iind+izm1,jind,0)
                         if (n1>0) then
                            temp = temp - n1*vi(iind-izm1,jind,0)
                         end if
                         deriveValue(center_i+(2*ssize)+ao12) = deriveValue(center_i+(2*ssize)+ao12) - temp * pfac

                         temp = 2.0*auxExponentB*vi(iind,jind+jzm1,0)
                         if (n2>0) then
                            temp = temp - n2*vi(iind,jind-jzm1,0)
                         end if
                         deriveValue(center_j+(2*ssize)+ao12) = deriveValue(center_j+(2*ssize)+ao12) - temp * pfac

                         deriveValue(3*ssize*center+2*ssize+ao12) = deriveValue(3*ssize*center+2*ssize+ao12) - vz(iind,jind,0) * pfac

                         ao12 = ao12 + 1
                      end do
                   end do
                end do
             end do
          end do
       end do
       ! write(*,"(A)") "------------------------------------ "
    end do
    ! write(*,"(A)") "------------------------------------ "
    ! !    write(*,"(A)") "Derivadas de cinetica"
    ! !    do i = 0, 6*ssize - 1
    ! !       write(*,"(f12.8)") deriveValue(i)
    ! !    end do

  end subroutine AttractionDerivatives_getDerive

  !>
  !! @brief Implementacion de la recursion propuesta por Obara - Saika para el caso de derivadas integrales de energia potencial
  !! @author J.M. Rodas, 2015
  !! @version 1.0
  subroutine AttractionDerivatives_obaraSaikaRecursion(vi, vx, vy, vz, PA, PB, PC, zeta, angularMomentIndexA, angularMomentIndexB)
    implicit none

    real(8), intent(inout), allocatable :: vi(:,:,:), vx(:,:,:), vy(:,:,:), vz(:,:,:)
    real(8), intent(in) :: PA(0:3), PB(0:3), PC(0:3)
    real(8), intent(in) :: zeta
    integer, intent(in) :: angularMomentIndexA, angularMomentIndexB

    integer :: a, b, m
    integer :: azm
    integer :: aym
    integer :: axm
    integer :: bzm
    integer :: bym
    integer :: bxm
    integer :: ax, ay, az, bx, by, bz
    integer :: aind, bind
    real(8) :: ooz 
    integer :: mmax
    real(8) :: tmp, u
    real(8), allocatable :: F(:)

    azm = 1
    aym = angularMomentIndexA + 1
    axm = aym * aym
    bzm = 1
    bym = angularMomentIndexB + 1
    bxm = bym * bym
    ooz = 1.0/(2.0 * zeta)

    mmax = angularMomentIndexA + angularMomentIndexB

    tmp = sqrt(zeta)*(2.0_8/Math_SQRT_PI)
    u = zeta * (PC(0) * PC(0) + PC(1) * PC(1) + PC(2) * PC(2))
    
    allocate(F(0:mmax+1))

    F = 0.0_8

    call Math_fgamma0(mmax,u,F)

    do m=0, mmax
       vi(0,0,m) = tmp*F(m)
    end do

    do m=0, mmax-1
       vx(0,0,m) = 2.0*zeta*PC(0)*vi(0,0,m+1)
       vy(0,0,m) = 2.0*zeta*PC(1)*vi(0,0,m+1)
       vz(0,0,m) = 2.0*zeta*PC(2)*vi(0,0,m+1)
    end do

    do b=1, angularMomentIndexB
       do bx=0, b
          do by=0, b-bx
             bz = b-bx-by
             bind = bx*bxm + by*bym + bz*bzm
             
             if(bz>0) then
                do m=0, mmax-b
                   vi(0,bind,m) = PB(2)*vi(0,bind-bzm,m)-PC(2)*vi(0,bind-bzm,m+1)
                end do
                do m=0,mmax-b-1
                   vx(0,bind,m) = PB(2)*vx(0,bind-bzm,m)-PC(2)*vx(0,bind-bzm,m+1)
                   vy(0,bind,m) = PB(2)*vy(0,bind-bzm,m)-PC(2)*vy(0,bind-bzm,m+1)
                   vz(0,bind,m) = PB(2)*vz(0,bind-bzm,m)-PC(2)*vz(0,bind-bzm,m+1) + vi(0,bind-bzm,m+1)
                end do
                if(bz>1) then
                   do m=0, mmax-b
                      vi(0,bind,m) = vi(0,bind,m) + ooz * (bz-1) * (vi(0,bind-2*bzm,m)-vi(0,bind-2*bzm,m+1))
                   end do
                   do m=0, mmax-b-1
                      vx(0,bind,m) = vx(0,bind,m) + ooz * (bz-1) * (vx(0,bind-2*bzm,m)-vx(0,bind-2*bzm,m+1))
                      vy(0,bind,m) = vy(0,bind,m) + ooz * (bz-1) * (vy(0,bind-2*bzm,m)-vy(0,bind-2*bzm,m+1))
                      vz(0,bind,m) = vz(0,bind,m) + ooz * (bz-1) * (vz(0,bind-2*bzm,m)-vz(0,bind-2*bzm,m+1))
                   end do
                end if
             else if(by>0) then
                do m=0, mmax-b
                   vi(0,bind,m) = PB(1)*vi(0,bind-bym,m)-PC(1)*vi(0,bind-bym,m+1)
                end do
                do m=0,mmax-b-1
                   vx(0,bind,m) = PB(1)*vx(0,bind-bym,m)-PC(1)*vx(0,bind-bym,m+1)
                   vy(0,bind,m) = PB(1)*vy(0,bind-bym,m)-PC(1)*vy(0,bind-bym,m+1) + vi(0,bind-bym,m+1)
                   vz(0,bind,m) = PB(1)*vz(0,bind-bym,m)-PC(1)*vz(0,bind-bym,m+1)
                end do
                if(by>1) then
                   do m=0, mmax-b
                      vi(0,bind,m) = vi(0,bind,m) + ooz * (by-1) * (vi(0,bind-2*bym,m)-vi(0,bind-2*bym,m+1))
                   end do
                   do m=0, mmax-b-1
                      vx(0,bind,m) = vx(0,bind,m) + ooz * (by-1) * (vx(0,bind-2*bym,m)-vx(0,bind-2*bym,m+1))
                      vy(0,bind,m) = vy(0,bind,m) + ooz * (by-1) * (vy(0,bind-2*bym,m)-vy(0,bind-2*bym,m+1))
                      vz(0,bind,m) = vz(0,bind,m) + ooz * (by-1) * (vz(0,bind-2*bym,m)-vz(0,bind-2*bym,m+1))
                   end do
                end if
             else if(bx>0) then
                do m=0, mmax-b
                   vi(0,bind,m) = PB(0)*vi(0,bind-bxm,m)-PC(0)*vi(0,bind-bxm,m+1)
                end do
                do m=0,mmax-b-1
                   vx(0,bind,m) = PB(0)*vx(0,bind-bxm,m)-PC(0)*vx(0,bind-bxm,m+1) + vi(0,bind-bxm,m+1)
                   vy(0,bind,m) = PB(0)*vy(0,bind-bxm,m)-PC(0)*vy(0,bind-bxm,m+1)
                   vz(0,bind,m) = PB(0)*vz(0,bind-bxm,m)-PC(0)*vz(0,bind-bxm,m+1)
                end do
                if(bx>1) then
                   do m=0, mmax-b
                      vi(0,bind,m) = vi(0,bind,m) + ooz * (bx-1) * (vi(0,bind-2*bxm,m)-vi(0,bind-2*bxm,m+1))
                   end do
                   do m=0, mmax-b-1
                      vx(0,bind,m) = vx(0,bind,m) + ooz * (bx-1) * (vx(0,bind-2*bxm,m)-vx(0,bind-2*bxm,m+1))
                      vy(0,bind,m) = vy(0,bind,m) + ooz * (bx-1) * (vy(0,bind-2*bxm,m)-vy(0,bind-2*bxm,m+1))
                      vz(0,bind,m) = vz(0,bind,m) + ooz * (bx-1) * (vz(0,bind-2*bxm,m)-vz(0,bind-2*bxm,m+1))
                   end do
                end if
             end if
          end do
       end do
    end do

    do b=0, angularMomentIndexB
       do bx=0, b
          do by=0, b-bx
             bz = b-bx-by
             bind = bx*bxm + by*bym + bz*bzm

             do a=1, angularMomentIndexA
                do ax=0, a
                   do ay=0, a-ax
                      az = a-ax-ay;
                      aind = ax*axm + ay*aym + az*azm
                      
                      if (az > 0) then
                         do m=0, mmax-a-b
                            vi(aind,bind,m) = PA(2) * vi(aind-azm,bind,m) - PC(2) * vi(aind-azm,bind,m+1)
                         end do
                         do m=0, mmax-a-b-1
                            vx(aind,bind,m) = PA(2) * vx(aind-azm,bind,m) - PC(2) * vx(aind-azm,bind,m+1)
                            vy(aind,bind,m) = PA(2) * vy(aind-azm,bind,m) - PC(2) * vy(aind-azm,bind,m+1)
                            vz(aind,bind,m) = PA(2) * vz(aind-azm,bind,m) - PC(2) * vz(aind-azm,bind,m+1) + vi(aind-azm,bind,m+1)
                         end do
                         if (az > 1) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * (az-1) * (vi(aind-2*azm,bind,m) - vi(aind-2*azm,bind,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * (az-1) * (vx(aind-2*azm,bind,m) - vx(aind-2*azm,bind,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * (az-1) * (vy(aind-2*azm,bind,m) - vy(aind-2*azm,bind,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * (az-1) * (vz(aind-2*azm,bind,m) - vz(aind-2*azm,bind,m+1))
                            end do
                         end if
                         if (bz > 0) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * bz * (vi(aind-azm,bind-bzm,m) - vi(aind-azm,bind-bzm,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * bz * (vx(aind-azm,bind-bzm,m) - vx(aind-azm,bind-bzm,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * bz * (vy(aind-azm,bind-bzm,m) - vy(aind-azm,bind-bzm,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * bz * (vz(aind-azm,bind-bzm,m) - vz(aind-azm,bind-bzm,m+1))
                            end do
                         end if
                      else if (ay > 0) then
                         do m=0, mmax-a-b
                            vi(aind,bind,m) = PA(1) * vi(aind-aym,bind,m) - PC(1) * vi(aind-aym,bind,m+1)
                         end do
                         do m=0, mmax-a-b-1
                            vx(aind,bind,m) = PA(1) * vx(aind-aym,bind,m) - PC(1) * vx(aind-aym,bind,m+1)
                            vy(aind,bind,m) = PA(1) * vy(aind-aym,bind,m) - PC(1) * vy(aind-aym,bind,m+1) + vi(aind-aym,bind,m+1)
                            vz(aind,bind,m) = PA(1) * vz(aind-aym,bind,m) - PC(1) * vz(aind-aym,bind,m+1)
                         end do
                         if (ay > 1) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * (ay-1) * (vi(aind-2*aym,bind,m) - vi(aind-2*aym,bind,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * (ay-1) * (vx(aind-2*aym,bind,m) - vx(aind-2*aym,bind,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * (ay-1) * (vy(aind-2*aym,bind,m) - vy(aind-2*aym,bind,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * (ay-1) * (vz(aind-2*aym,bind,m) - vz(aind-2*aym,bind,m+1))
                            end do
                         end if
                         if (by > 0) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * by * (vi(aind-aym,bind-bym,m) - vi(aind-aym,bind-bym,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * by * (vx(aind-aym,bind-bym,m) - vx(aind-aym,bind-bym,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * by * (vy(aind-aym,bind-bym,m) - vy(aind-aym,bind-bym,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * by * (vz(aind-aym,bind-bym,m) - vz(aind-aym,bind-bym,m+1))
                            end do
                         end if
                      else if (ax > 0) then
                         do m=0, mmax-a-b
                            vi(aind,bind,m) = PA(0) * vi(aind-axm,bind,m) - PC(0) * vi(aind-axm,bind,m+1)
                         end do
                         do m=0, mmax-a-b-1
                            vx(aind,bind,m) = PA(0) * vx(aind-axm,bind,m) - PC(0) * vx(aind-axm,bind,m+1) + vi(aind-axm,bind,m+1)
                            vy(aind,bind,m) = PA(0) * vy(aind-axm,bind,m) - PC(0) * vy(aind-axm,bind,m+1)
                            vz(aind,bind,m) = PA(0) * vz(aind-axm,bind,m) - PC(0) * vz(aind-axm,bind,m+1)
                         end do
                         if (ax > 1) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * (ax-1) * (vi(aind-2*axm,bind,m) - vi(aind-2*axm,bind,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * (ax-1) * (vx(aind-2*axm,bind,m) - vx(aind-2*axm,bind,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * (ax-1) * (vy(aind-2*axm,bind,m) - vy(aind-2*axm,bind,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * (ax-1) * (vz(aind-2*axm,bind,m) - vz(aind-2*axm,bind,m+1))
                            end do
                         end if
                         if (bx > 0) then
                            do m=0, mmax-a-b
                               vi(aind,bind,m) = vi(aind,bind,m) + ooz * bx * (vi(aind-axm,bind-bxm,m) - vi(aind-axm,bind-bxm,m+1))
                            end do
                            do m=0, mmax-a-b-1
                               vx(aind,bind,m) = vx(aind,bind,m) + ooz * bx * (vx(aind-axm,bind-bxm,m) - vx(aind-axm,bind-bxm,m+1))
                               vy(aind,bind,m) = vy(aind,bind,m) + ooz * bx * (vy(aind-axm,bind-bxm,m) - vy(aind-axm,bind-bxm,m+1))
                               vz(aind,bind,m) = vz(aind,bind,m) + ooz * bx * (vz(aind-axm,bind-bxm,m) - vz(aind-axm,bind-bxm,m+1))
                            end do
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    
    deallocate(F)

  end subroutine AttractionDerivatives_obaraSaikaRecursion

  ! !>
  ! !! @author R. Gonzalez, rogonzalez@unal.edu.co   2014
  ! function AttractionDerivatives_int( x, y, z,  auxExponentA, l1, m1, n1, auxExponentB, l2, m2, n2) result (output)
  !   implicit none

  !   integer, intent(in) :: l1, m1, n1
  !   integer, intent(in) :: l2, m2, n2
  !   real(8), intent(in) :: auxExponentA, auxExponentB
  !   real(8), intent(in), allocatable ::  x(:,:), y(:,:), z(:,:)
  !   real(8) :: output

  !   real(8) :: I1, I2, I3, I4
  !   real(8) :: Ix, Iy, Iz


  !   if (l1 == 0 .or. l2 == 0) then
  !      I1 = 0.0_8
  !   else
  !      I1 = x(l1-1,l2-1) * y(m1,m2) * z(n1,n2)
  !   end if

  !   I2 = x(l1+1,l2+1) * y(m1,m2) * z(n1,n2)

  !   if (l2 == 0) then
  !      I3 = 0.0_8
  !   else
  !      I3 = x(l1+1,l2-1) * y(m1,m2) * z(n1,n2)
  !   end if

  !   if (l1 == 0) then
  !      I4 = 0.0_8
  !   else
  !      I4 = x(l1-1,l2+1) * y(m1,m2) * z(n1,n2)
  !   end if

  !   Ix = 0.5 * l1 * l2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * l2 * I3 - l1 * auxExponentB * I4

  !   if (m1 == 0 .or. m2 == 0) then
  !      I1 = 0.0_8
  !   else
  !      I1 =   x(l1,l2) * y(m1-1,m2-1) * z(n1,n2)
  !   end if

  !   I2 = x(l1,l2) * y(m1+1,m2+1) * z(n1,n2)

  !   if (m2 == 0) then
  !      I3 = 0.0_8
  !   else
  !      I3 = x(l1,l2) * y(m1+1,m2-1) * z(n1,n2)
  !   end if

  !   if (m1 == 0) then
  !      I4 = 0.0_8
  !   else
  !      I4 = x(l1,l2) * y(m1-1,m2+1) * z(n1,n2)
  !   end if

  !   Iy = 0.5 * m1 * m2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * m2 * I3 - m1 * auxExponentB * I4

  !   if (n1 == 0 .or. n2 == 0) then
  !      I1 = 0.0_8
  !   else
  !      I1 = x(l1,l2) * y(m1,m2) * z(n1-1,n2-1)
  !   end if

  !   I2 = x(l1,l2) * y(m1,m2) * z(n1+1,n2+1)

  !   if (n2 == 0) then
  !      I3 = 0.0_8
  !   else
  !      I3 = x(l1,l2) * y(m1,m2) * z(n1+1,n2-1)
  !   end if

  !   if (n1 == 0) then
  !      I4 = 0.0_8
  !   else
  !      I4 = x(l1,l2) * y(m1,m2) * z(n1-1,n2+1)
  !      Iz = 0.5 * n1 * n2 * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * n2 * I3 - n1 * auxExponentB * I4
  !   end if

  !   output = Ix + Iy + Iz
  !   !    write(*,"(A,f)") "Integral: ", output

  ! end function AttractionDerivatives_int

end module AttractionDerivatives_
