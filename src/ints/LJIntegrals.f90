!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief   Modulo para calculo de integrales tipo <\phi | (1-EXP[-alpha*r^2]^n)* 1/r^k | \phi>
!!          donde k= 6 y 12 
!!
module LJIntegrals_
  use Exception_
  use Math_
  use ContractedGaussian_
  use LJPotential_
  implicit none

!  !>
!  !! Point  atributes
!  type, public :: pointsorigin
!     real(16) :: x
!     real(16) :: y
!     real(16) :: z
!  end type pointsorigin

  public ::  &
     !  LJIntegrals_computeShell, &
     !  LJIntegrals_computePrimitive
     !  LJIntegral_computesix, &
       LJIntegrals_computeShell
  private :: &
       LJIntegrals_obaraSaikaRecursion, &
       LJIntegrals_primitive6, &
       LJIntegrals_primitive12

contains
  !>
  !! Calculates point charge - quantum specie integral between two contractions (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2013.02.04: E.F.Posada: change for use in opints
  !! @return  output: attraction integral of a shell (all combinations)
  !! @version 1.0
!  subroutine LJIntegrals_computeShell(contractedGaussianA, contractedGaussianB, pointorigin, npoints, integral)
!    implicit none
!
!    type(ContractedGaussian), intent(in) :: contractedGaussianA, contractedGaussianB
!    type(pointorigin), intent(in), allocatable :: pointorigin(:)
!
!    integer, intent(in) :: npoints
!    real(8), intent(inout) :: integral(contractedGaussianA%numCartesianOrbital * contractedGaussianB%numCartesianOrbital)
!
!    integer ::  am1(0:3)
!    integer ::  am2(0:3)
!    integer ::  nprim1
!    integer ::  nprim2
!    real(8) ::  A(0:3)
!    real(8) ::  B(0:3)
!    real(8) ::  exp1(0:contractedGaussianA%length)
!    real(8) ::  exp2(0:contractedGaussianB%length)
!    real(8) ::  coef1(0:contractedGaussianA%length)
!    real(8) ::  coef2(0:contractedGaussianB%length)
!    real(8) ::  nor1(0:contractedGaussianA%length)
!    real(8) ::  nor2(0:contractedGaussianB%length)
!    real(8) :: auxIntegral
!    integer, allocatable :: angularMomentIndexA(:,:)
!    integer, allocatable :: angularMomentIndexB(:,:)
!    integer ::  i, m, p, q
!
!    integral = 0.0_8
!    auxIntegral = 0.0_8
!
!    if(allocated(angularMomentIndexA)) deallocate(angularMomentIndexA)
!    if(allocated(angularMomentIndexB)) deallocate(angularMomentIndexB)
!
!    allocate(angularMomentIndexA(3, contractedGaussianA%numCartesianOrbital))
!    allocate(angularMomentIndexB(3, contractedGaussianB%numCartesianOrbital))
!
!    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexA, contractedGaussianA)
!    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexB, contractedGaussianB)
!
!
!
!    nprim1 = contractedGaussianA%length
!    A(0) = contractedGaussianA%origin(1)
!    A(1) = contractedGaussianA%origin(2)
!    A(2) = contractedGaussianA%origin(3)
!    coef1(0:nprim1-1) =  contractedGaussianA%contractionCoefficients(1:nprim1)
!
!
!    nprim2 = contractedGaussianB%length
!    B(0) = contractedGaussianB%origin(1)
!    B(1) = contractedGaussianB%origin(2)
!    B(2) = contractedGaussianB%origin(3)
!    coef2(0:nprim2-1) =  contractedGaussianB%contractionCoefficients(1:nprim2)
!
!    m = 0
!
!
!    do p = 1, contractedGaussianA%numcartesianOrbital
!       do q = 1, contractedGaussianB%numcartesianOrbital
!
!          m = m + 1
!
!          exp1(0:nprim1-1) = contractedGaussianA%orbitalExponents(1:nprim1)
!          nor1(0:nprim1-1) = contractedGaussianA%primNormalization(1:nprim1,p)
!
!          exp2(0:nprim2-1) = contractedGaussianB%orbitalExponents(1:nprim2)
!          nor2(0:nprim2-1) = contractedGaussianB%primNormalization(1:nprim2,q)
!
!          am1 = 0
!          am2 = 0
!
!          am1(0:2) = angularMomentIndexA(1:3, p)
!          am2(0:2) = angularMomentIndexB(1:3, q)
!
!          call LJIntegrals_computePrimitive(am1, am2, nprim1, nprim2, npoints, A, B, exp1, exp2, coef1, coef2, nor1, nor2, point, auxintegral)
!
!
!          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
!               * contractedGaussianB%contNormalization(q)
!
!          integral(m) = auxIntegral
!
!       end do
!       end do
!
!  end subroutine LJIntegrals_computeShell

  !>
  !! @brief Evaluates attraction - repulsion integrals for any angular momentum
  !! @author Edwin Posada, 2010
  !! @return integral values for all shell (all possibles combinations if angular momentum index) (output)
  subroutine LJIntegrals_computeShell(contractedGaussianA, contractedGaussianB, contractedGaussianC, atomsCenter, numberOfPoints, ljparameters, integralValue)

    implicit none 
     
    type(ContractedGaussian) :: contractedGaussianA, contractedGaussianB, contractedGaussianC
!   type(pointsorigin) :: pointorigin(1:numberOfPoints)
    type(LJPot) :: potential
    integer :: potId
    real(16), allocatable :: atomsCenter(:,:)
    real(16), allocatable :: ljparameters(:,:)
    integer, intent(in) :: numberOfPoints 
    real(16) :: cexp6, cexp12, cp2, cpap, cp2_sqrt, cp2cpap
    real(16) :: w, c   
    real(16) :: exp1, exp2, exp3, classicr6, tmp(32), primIntegralValue6  
    integer :: i, ia, ib
    integer :: atom
  
    integer :: iq, ic, l, f 
    real(16) :: integralValue6, integralValue12
    real(16), allocatable :: auxValue6(:,:), classicValues6(:,:), auxValue12(:,:), classicValues12(:,:)
    real(16), intent(inout) :: integralValue(contractedGaussianA%numCartesianOrbital * contractedGaussianB%numCartesianOrbital)
   
    cexp6 = 1.0_16 !EXPONENTE DEL CUTOFF PARA r6
    cexp12 = 4.0_16 !EXPONENTE DEL CUTOFF PARA r12
     
        call ContractedGaussian_product(contractedGaussianA, contractedGaussianB, contractedGaussianC)
        
!                  print*, "primNormalizationA", contractedGaussianA%primNormalization
!                  print*, "contNormalizationA", contractedGaussianA%contNormalization
!                  print*, "primNormalizationB", contractedGaussianB%primNormalization
!                  print*, "contNormalizationB", contractedGaussianB%contNormalization
!        do f=1, contractedGaussianC%length
!            print*, "origin", contractedGaussianC%primOrigin(f,:)
!        end do
        
!                print*, "ORIGEN BASE A", contractedGaussianA%origin(:)
!                print*, "ORIGEN BASE B", contractedGaussianB%origin(:)

        if(allocated(auxValue6)) deallocate(auxValue6)
        allocate(auxValue6(numberOfPoints,contractedGaussianC%length))

        if(allocated(classicValues6)) deallocate(classicValues6)
        allocate(classicValues6(numberOfPoints,contractedGaussianC%length))
        
        if(allocated(auxValue12)) deallocate(auxValue12)
        allocate(auxValue12(numberOfPoints,contractedGaussianC%length))

        if(allocated(classicValues12)) deallocate(classicValues12)
        allocate(classicValues12(numberOfPoints,contractedGaussianC%length))
        
        auxValue6 = 0.0_8
        auxValue12 = 0.0_8

        do atom = 1, numberOfPoints

            do ic = 1, contractedGaussianC%length

                !CONVENCIONES: cp2 = CP^2
                cp2 = 0.0_16
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,1) - atomsCenter(atom,1)) * (contractedGaussianC%primOrigin(ic,1) - atomsCenter(atom,1))
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,2) - atomsCenter(atom,2)) * (contractedGaussianC%primOrigin(ic,2) - atomsCenter(atom,2))
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,3) - atomsCenter(atom,3)) * (contractedGaussianC%primOrigin(ic,3) - atomsCenter(atom,3))
                
!                print*, "ORIGEN BASE C", contractedGaussianC%primOrigin(:,:)
!                print*, "CENTRO ATOMO", atomsCenter(:,:)
                cp2_sqrt = sqrt(cp2)

                c=real(cp2_sqrt,16)
                w= real(contractedGaussianC%orbitalExponents(ic),16)
                
!                print*, "CP", c
!                print*, "Exponent", w

                call LJIntegrals_primitive6(c, w, cexp6, auxValue6(atom,ic))

                classicValues6(atom,ic) = ((1.0_8-EXP(-cexp6*c**2))**3)*(1.0_8/c**6)

!                print*, "auxValue6", auxValue6(atom,ic)
!                print*, "classicValue6", classicValues6(atom,ic)

                call LJIntegrals_primitive12(c, w, cexp12, auxValue12(atom,ic))

                classicValues12(atom,ic) = ((1.0_8-EXP(-cexp12*c**2))**6)*(1.0_8/c**12)

!                print*, "auxValue12", auxValue12(atom,ic)
!                print*, "classicValue12", classicValues12(atom,ic)
            end do
             
        end do

      do l=1, contractedGaussianC%numCartesianOrbital
          do atom=1, numberOfPoints
              do ic =1, contractedGaussianC%length

                cp2 = 0.0_16
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,1) - atomsCenter(atom,1)) * (contractedGaussianC%primOrigin(ic,1) - atomsCenter(atom,1))
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,2) - atomsCenter(atom,2)) * (contractedGaussianC%primOrigin(ic,2) - atomsCenter(atom,2))
                cp2 = cp2 + (contractedGaussianC%primOrigin(ic,3) - atomsCenter(atom,3)) * (contractedGaussianC%primOrigin(ic,3) - atomsCenter(atom,3))
                
                cp2_sqrt = sqrt(cp2)

                c=real(cp2_sqrt,16)
                w= real(contractedGaussianC%orbitalExponents(ic),16)
                
!                  print*, "A PARAMETER", ljparameters(atom,1)
!                  print*, "B PARAMETER", ljparameters(atom,2)

                  If (auxValue6(atom,ic) < 1.0E-10_8) then 
!                        print *, "Warning, the value integral for r6 is less than 1E-10, so shall be calculated classic form"
                        integralValue(l) = integralValue(l) - classicValues6(atom,ic)*ljparameters(atom,2) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  else if (w>=500.0_16 .and. c>=5.0_16) then
!                          print*, "Exponent > 500 and c > 5"
                        integralValue(l) = integralValue(l) - classicValues6(atom,ic)*ljparameters(atom,2) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  else if (c>=50.0_16) then
                        integralValue(l) = integralValue(l) - classicValues6(atom,ic)*ljparameters(atom,2) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  else       
                        integralValue(l) = integralValue(l) - auxValue6(atom,ic)*ljparameters(atom,2) * contractedGaussianC%primNormalization(ic, l) *  contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  end if

                  !!Condiciones de evaluación para 1/r^12                

                  If (auxValue12(atom,ic) < 1.0E-10_8) then 
!                      print *, "Warning, the value integral for r12 is less than 1E-10, so shall be calculated classic form"
                      integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)!
                 ! else if (w>=1.0_16) then
                 !     print*, "Exponent is greater than 1, so calculated with r6^2"
                 !     integralValue(l) = integralValue(l) +((auxValue6(atom,ic)*contractedGaussianC%primNormalization(ic,l))**2)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 !  integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)!
                 else if (c>=50.0_16) then
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 else if(w >= 50.0_16) then
!                     print*, "Exponent is greater than 50, change for Classic Value 12"
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 else if (w>= 20.0_16 .and. c>=2.4_16) then  
!                     print*, "Exponent is greater than 20 and c is greater than 2"
                     !integralValue(l) = integralValue(l) +((auxValue6(atom,ic)*contractedGaussianC%primNormalization(ic,l))**2)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 else if (w>= 10.0_16 .and. c>=4.0_16) then  
!                     print*, "Exponent is greater than 10 and c is greater than 3"
                     !integralValue(l) = integralValue(l) +((auxValue6(atom,ic)*contractedGaussianC%primNormalization(ic,l))**2)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 else if (w>= 5.0_16 .and. c>=5.0_16) then  
!                     print*, "Exponent is greater than 5 and c is greater than 4"
                     !integralValue(l) = integralValue(l) +((auxValue6(atom,ic)*contractedGaussianC%primNormalization(ic,l))**2)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                 else if (w>= 1.0_16 .and. c>=20.0_16) then
!                     print*, "Exponent is greater than 1 and c is greater than 20"
                     !integralValue(l) = integralValue(l) +((auxValue6(atom,ic)*contractedGaussianC%primNormalization(ic,l))**2)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                     integralValue(l) = integralValue(l) + classicValues12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  else 
                      integralValue(l) = integralValue(l) + auxValue12(atom,ic)*ljparameters(atom,1) * contractedGaussianC%primNormalization(ic, l) *  contractedGaussianC%contNormalization(l) * contractedGaussianC%contractionCoefficients(ic)
                  end if

              end do
          end do
      end do
      
!      print*, "Integral Value", integralValue

  end subroutine LJIntegrals_computeShell

  subroutine LJIntegrals_primitive6(c, w, cexp6, auxValue6)

      implicit none
      
      real(16) :: w, c, exp1, exp2, exp3, cexp6, tmp(32) 
      real(16) :: auxValue6

      exp1=EXP(-(c**2*w) + (c**2*w**2)/(cexp6 + w))
      exp2=EXP(-(c**2*w) + (c**2*w**2)/(2.0_16*cexp6 + w))
      exp3=EXP(-(c**2*w) + (c**2*w**2)/(3.0_16*cexp6 + w))

      tmp = 0.0_16
      tmp(1) = (5.0_16*w**1.5)/3.0_16 
      tmp(2) = -(2.0_16*c**2*w**2.5)/3.0_16 
      tmp(3) = -5.0_16*exp1*cexp6*Sqrt(cexp6 + w) 
      tmp(4) = -5.0_16*exp1*w*Sqrt(cexp6 + w) 
      tmp(5) = 2.0_16*c**2*exp1*w**2*Sqrt(cexp6 + w) 
      tmp(6) = 10.0_16*exp2*cexp6*Sqrt(2.0_16*cexp6 + w) 
      tmp(7) = 5.0_16*exp2*w*Sqrt(2.0_16*cexp6 + w) 
      tmp(8) = -2.0_16*c**2*exp2*w**2*Sqrt(2.0_16*cexp6 + w) 
      tmp(9) = - 5.0_16*exp3*cexp6*Sqrt(3.0_16*cexp6 + w) 
      tmp(10)= -(5.0_16*exp3*w*Sqrt(3.0_16*cexp6 + w))/3.0_16 
      tmp(11) = (2.0_16*c**2*exp3*w**2*Sqrt(3.0_16*cexp6 + w))/3.0_16
      tmp(12) = (w*Math_Dawson(c*Sqrt(w)))/c 
      tmp(13) = -4.0_16*c*w**2*Math_Dawson(c*Sqrt(w)) 
      tmp(14) = (4.0_16*c**3*w**3*Math_Dawson(c*Sqrt(w)))/3.0_16 
      tmp(15) = -(6.0_16*exp1*cexp6*Math_Dawson((c*w)/Sqrt(cexp6 + w)))/c 
      tmp(16) = -(3.0_16*exp1*cexp6**2*Math_Dawson((c*w)/Sqrt(cexp6 + w)))/(c*w) 
      tmp(17) = -(3.0_16*exp1*w*Math_Dawson((c*w)/Sqrt(cexp6 + w)))/c 
      tmp(18) = 12.0_16*c*exp1*cexp6*w*Math_Dawson((c*w)/Sqrt(cexp6 + w)) 
      tmp(19) = 12.0_16*c*exp1*w**2*Math_Dawson((c*w)/Sqrt(cexp6 + w)) 
      tmp(20) = -4.0_16*c**3*exp1*w**3*Math_Dawson((c*w)/Sqrt(cexp6 + w)) 
      tmp(21) = (12.0_16*exp2*cexp6*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)))/c 
      tmp(22) = (12.0_16*exp2*cexp6**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)))/(c*w) 
      tmp(23) = (3.0_16*exp2*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)))/c 
      tmp(24) = -24.0_16*c*exp2*cexp6*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)) 
      tmp(25) = -12.0_16*c*exp2*w**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)) 
      tmp(26) = 4.0_16*c**3*exp2*w**3*Math_Dawson((c*w)/Sqrt(2.0_16*cexp6 + w)) 
      tmp(27) = -(6.0_16*exp3*cexp6*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)))/c 
      tmp(28) = -(9.0_16*exp3*cexp6**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)))/(c*w) 
      tmp(29) = -(exp3*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)))/c 
      tmp(30) = 12.0_16*c*exp3*cexp6*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)) 
      tmp(31) = 4.0_16*c*exp3*w**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)) 
      tmp(32) = -(4.0_16*c**3*exp3*w**3*Math_Dawson((c*w)/Sqrt(3.0_16*cexp6 + w)))/3.0_16
   
      call sort(tmp, int(size(tmp), 8), int(size(tmp), 8))
    
      auxValue6 = sum_kahan(tmp) * (Math_PI_MOD**1.5)

  end subroutine LJIntegrals_primitive6



  subroutine LJIntegrals_primitive12(c, w, cexp12, auxValue12)

      implicit none
      
      real(16) :: w, c, exp1, exp2, exp3, exp4, exp5, exp6,  cexp12, tmp(227) 
      real(16) :: auxValue12


      exp1=EXP(-(c**2*w) + (c**2*w**2)/(cexp12 + w))
      exp2=EXP(-(c**2*w) + (c**2*w**2)/(2*cexp12 + w))
      exp3=EXP(-(c**2*w) + (c**2*w**2)/(3*cexp12 + w))
      exp4=EXP(-(c**2*w) + (c**2*w**2)/(4*cexp12 + w))
      exp5=EXP(-(c**2*w) + (c**2*w**2)/(5*cexp12 + w))
      exp6=EXP(-(c**2*w) + (c**2*w**2)/(6*cexp12 + w))
      
      tmp = 0.0_16
      tmp(1) = (-193.0_16*w**4.5)/3780.0_16
      tmp(2) = (88.0_16*c**2*w**5.5)/945.0_16
      tmp(3) = - (28.0_16*c**4*w**6.5)/675.0_16
      tmp(4) = (88.0_16*c**6*w**7.5)/14175.0_16
      tmp(5) = - (4.0_16*c**8*w**8.5)/14175.0_16 
      tmp(6) = (193.0_16*exp1*cexp12**4*Sqrt(cexp12 +w))/630.0_16
      tmp(7) = (386.0_16*exp1*cexp12**3*w*Sqrt(cexp12 + w))/315.0_16
      tmp(8) = (193.0_16*exp1*cexp12**2*w**2*Sqrt(cexp12 + w))/105.0_16
      tmp(9) = -(176.0_16*c**2*exp1*cexp12**3*w**2*Sqrt(cexp12 + w))/315.0_16
      tmp(10) =  (386.0_16*exp1*cexp12*w**3*Sqrt(cexp12 + w))/315.0_16
      tmp(11) = -(176.0_16*c**2*exp1*cexp12**2*w**3*Sqrt(cexp12 + w))/105.0_16
      tmp(12) =  (193.0_16*exp1*w**4*Sqrt(cexp12 + w))/630.0_16
      tmp(13) = -(176.0_16*c**2*exp1*cexp12*w**4*Sqrt(cexp12 + w))/105.0_16
      tmp(14) =  (56.0_16*c**4*exp1*cexp12**2*w**4*Sqrt(cexp12 + w))/225.0_16
      tmp(15) = -(176.0_16*c**2*exp1*w**5*Sqrt(cexp12 + w))/315.0_16
      tmp(16) =  (112.0_16*c**4*exp1*cexp12*w**5*Sqrt(cexp12 + w))/225.0_16
      tmp(17) = (56.0_16*c**4*exp1*w**6*Sqrt(cexp12 + w))/225.0_16
      tmp(18) = - (176.0_16*c**6*exp1*cexp12*w**6*Sqrt(cexp12 + w))/4725.0_16
      tmp(19) = -(176.0_16*c**6*exp1*w**7*Sqrt(cexp12 + w))/4725.0_16
      tmp(20) =  (8.0_16*c**8*exp1*w**8*Sqrt(cexp12 + w))/4725.0_16
      tmp(21) = -(772.0_16*exp2*cexp12**4*Sqrt(2.0_16*cexp12 + w))/63.0_16
      tmp(22) = - (1544.0_16*exp2*cexp12**3*w*Sqrt(2.0_16*cexp12 + w))/63.0_16
      tmp(23) = -(386.0_16*exp2*cexp12**2*w**2*Sqrt(2.0_16*cexp12 + w))/21.0_16
      tmp(24) =  (704.0_16*c**2*exp2*cexp12**3*w**2*Sqrt(2.0_16*cexp12 + w))/63.0_16
      tmp(25) = -(386.0_16*exp2*cexp12*w**3*Sqrt(2.0_16*cexp12 + w))/63.0_16
      tmp(26) =  (352.0_16*c**2*exp2*cexp12**2*w**3*Sqrt(2.0_16*cexp12 + w))/21.0_16
      tmp(27) = -(193.0_16*exp2*w**4*Sqrt(2.0_16*cexp12 + w))/252.0_16
      tmp(28) =  (176.0_16*c**2*exp2*cexp12*w**4*Sqrt(2.0_16*cexp12 + w))/21.0_16
      tmp(29) = -(112.0_16*c**4*exp2*cexp12**2*w**4*Sqrt(2.0_16*cexp12 + w))/45.0_16
      tmp(30) =  (88.0_16*c**2*exp2*w**5*Sqrt(2.0_16*cexp12 + w))/63.0_16
      tmp(31) = -(112.0_16*c**4*exp2*cexp12*w**5*Sqrt(2.0_16*cexp12 + w))/45.0_16
      tmp(32) = - (28.0_16*c**4*exp2*w**6*Sqrt(2.0_16*cexp12 + w))/45.0_16
      tmp(33) = (176.0_16*c**6*exp2*cexp12*w**6*Sqrt(2.0_16*cexp12 + w))/945.0_16
      tmp(34) =  (88.0_16*c**6*exp2*w**7*Sqrt(2.0_16*cexp12 + w))/945.0_16
      tmp(35) = -(4.0_16*c**8*exp2*w**8*Sqrt(2.0_16*cexp12 + w))/945.0_16
      tmp(36) =  (579.0_16*exp3*cexp12**4*Sqrt(3.0_16*cexp12 + w))/7.0_16
      tmp(37) = (772.0_16*exp3*cexp12**3*w*Sqrt(3.0_16*cexp12 + w))/7.0_16
      tmp(38) =  (386.0_16*exp3*cexp12**2*w**2*Sqrt(3.0_16*cexp12 + w))/7.0_16
      tmp(39) = -(352.0_16*c**2*exp3*cexp12**3*w**2*Sqrt(3.0_16*cexp12 + w))/7.0_16
      tmp(40) =  (772.0_16*exp3*cexp12*w**3*Sqrt(3.0_16*cexp12 + w))/63.0_16
      tmp(41) = -(352.0_16*c**2*exp3*cexp12**2*w**3*Sqrt(3.0_16*cexp12 + w))/7.0_16
      tmp(42) =  (193.0_16*exp3*w**4*Sqrt(3.0_16*cexp12 + w))/189.0_16
      tmp(43) = -(352.0_16*c**2*exp3*cexp12*w**4*Sqrt(3.0_16*cexp12 + w))/21.0_16
      tmp(44) =  (112.0_16*c**4*exp3*cexp12**2*w**4*Sqrt(3.0_16*cexp12 + w))/15.0_16
      tmp(45) = -(352.0_16*c**2*exp3*w**5*Sqrt(3.0_16*cexp12 + w))/189.0_16
      tmp(46) =  (224.0_16*c**4*exp3*cexp12*w**5*Sqrt(3.0_16*cexp12 + w))/45.0_16
      tmp(47) = (112.0_16*c**4*exp3*w**6*Sqrt(3.0_16*cexp12 + w))/135.0_16
      tmp(48) = - (352.0_16*c**6*exp3*cexp12*w**6*Sqrt(3.0_16*cexp12 + w))/945.0_16
      tmp(49) = -(352.0_16*c**6*exp3*w**7*Sqrt(3.0_16*cexp12 + w))/2835.0_16
      tmp(50) =  (16.0_16*c**8*exp3*w**8*Sqrt(3.0_16*cexp12 + w))/2835.0_16
      tmp(51) = -(12352.0_16*exp4*cexp12**4*Sqrt(4.0_16*cexp12 + w))/63.0_16 
      tmp(52) = -(12352.0_16*exp4*cexp12**3*w*Sqrt(4.0_16*cexp12 + w))/63.0_16
      tmp(53) = - (1544.0_16*exp4*cexp12**2*w**2*Sqrt(4.0_16*cexp12 + w))/21.0_16
      tmp(54) = (5632.0_16*c**2*exp4*cexp12**3*w**2*Sqrt(4.0_16*cexp12 + w))/63.0_16
      tmp(55) = - (772.0_16*exp4*cexp12*w**3*Sqrt(4.0_16*cexp12 + w))/63.0_16
      tmp(56) = (1408.0_16*c**2*exp4*cexp12**2*w**3*Sqrt(4.0_16*cexp12 + w))/21.0_16
      tmp(57) = - (193.0_16*exp4*w**4*Sqrt(4.0_16*cexp12 + w))/252.0_16
      tmp(58) = (352.0_16*c**2*exp4*cexp12*w**4*Sqrt(4.0_16*cexp12 + w))/21.0_16
      tmp(59) = - (448.0_16*c**4*exp4*cexp12**2*w**4*Sqrt(4.0_16*cexp12 + w))/45.0_16
      tmp(60) = (88.0_16*c**2*exp4*w**5*Sqrt(4.0_16*cexp12 + w))/63.0_16
      tmp(61) = - (224.0_16*c**4*exp4*cexp12*w**5*Sqrt(4.0_16*cexp12 + w))/45.0_16
      tmp(62) = -(28.0_16*c**4*exp4*w**6*Sqrt(4.0_16*cexp12 + w))/45.0_16
      tmp(63) =  (352.0_16*c**6*exp4*cexp12*w**6*Sqrt(4.0_16*cexp12 + w))/945.0_16
      tmp(64) = (88.0_16*c**6*exp4*w**7*Sqrt(4.0_16*cexp12 + w))/945.0_16
      tmp(65) = - (4.0_16*c**8*exp4*w**8*Sqrt(4.0_16*cexp12 + w))/945.0_16
      tmp(66) = (24125.0_16*exp5*cexp12**4*Sqrt(5.0_16*cexp12 + w))/126.0_16
      tmp(67) =  (9650.0_16*exp5*cexp12**3*w*Sqrt(5.0_16*cexp12 + w))/63.0_16
      tmp(68) = (965.0_16*exp5*cexp12**2*w**2*Sqrt(5.0_16*cexp12 + w))/21.0_16
      tmp(69) = - (4400.0_16*c**2*exp5*cexp12**3*w**2*Sqrt(5.0_16*cexp12 + w))/63.0_16 
      tmp(70) = (386.0_16*exp5*cexp12*w**3*Sqrt(5.0_16*cexp12 + w))/63.0_16
      tmp(71) = - (880.0_16*c**2*exp5*cexp12**2*w**3*Sqrt(5.0_16*cexp12 + w))/21.0_16
      tmp(72) = (193.0_16*exp5*w**4*Sqrt(5.0_16*cexp12 + w))/630.0_16
      tmp(73) = - (176.0_16*c**2*exp5*cexp12*w**4*Sqrt(5.0_16*cexp12 + w))/21.0_16
      tmp(74) = (56.0_16*c**4*exp5*cexp12**2*w**4*Sqrt(5.0_16*cexp12 + w))/9.0_16
      tmp(75) = - (176.0_16*c**2*exp5*w**5*Sqrt(5.0_16*cexp12 + w))/315.0_16
      tmp(76) = (112.0_16*c**4*exp5*cexp12*w**5*Sqrt(5.0_16*cexp12 + w))/45.0_16
      tmp(77) =  (56.0_16*c**4*exp5*w**6*Sqrt(5.0_16*cexp12 + w))/225.0_16
      tmp(78) = -(176.0_16*c**6*exp5*cexp12*w**6*Sqrt(5.0_16*cexp12 + w))/945.0_16
      tmp(79) = - (176.0_16*c**6*exp5*w**7*Sqrt(5.0_16*cexp12 + w))/4725.0_16
      tmp(80) = (8.0_16*c**8*exp5*w**8*Sqrt(5.0_16*cexp12 + w))/4725.0_16
      tmp(81) = - (2316.0_16*exp6*cexp12**4*Sqrt(6.0_16*cexp12 + w))/35.0_16
      tmp(82) = -(1544.0_16*exp6*cexp12**3*w*Sqrt(6.0_16*cexp12 + w))/35.0_16
      tmp(83) = - (386.0_16*exp6*cexp12**2*w**2*Sqrt(6.0_16*cexp12 + w))/35.0_16
      tmp(84) = (704.0_16*c**2*exp6*cexp12**3*w**2*Sqrt(6.0_16*cexp12 + w))/35.0_16
      tmp(85) = - (386.0_16*exp6*cexp12*w**3*Sqrt(6.0_16*cexp12 + w))/315.0_16
      tmp(86) = (352.0_16*c**2*exp6*cexp12**2*w**3*Sqrt(6.0_16*cexp12 + w))/35.0_16
      tmp(87) = - (193.0_16*exp6*w**4*Sqrt(6.0_16*cexp12 + w))/3780.0_16
      tmp(88) = (176.0_16*c**2*exp6*cexp12*w**4*Sqrt(6.0_16*cexp12 + w))/105.0_16
      tmp(89) = - (112.0_16*c**4*exp6*cexp12**2*w**4*Sqrt(6.0_16*cexp12 + w))/75.0_16
      tmp(90) = (88.0_16*c**2*exp6*w**5*Sqrt(6.0_16*cexp12 + w))/945.0_16
      tmp(91) = - (112.0_16*c**4*exp6*cexp12*w**5*Sqrt(6.0_16*cexp12 + w))/225.0_16
      tmp(92) = -(28.0_16*c**4*exp6*w**6*Sqrt(6.0_16*cexp12 + w))/675.0_16
      tmp(93) =  (176.0_16*c**6*exp6*cexp12*w**6*Sqrt(6.0_16*cexp12 + w))/4725.0_16
      tmp(94) = (88.0_16*c**6*exp6*w**7*Sqrt(6.0_16*cexp12 + w))/14175.0_16
      tmp(95) = - (4.0_16*c**8*exp6*w**8*Sqrt(6.0_16*cexp12 + w))/14175.0_16
      tmp(96) = -(w**4*Math_Dawson(c*Sqrt(w)))/(60.0_16*c)
      tmp(97) =  (c*w**5*Math_Dawson(c*Sqrt(w)))/6.0_16
      tmp(98) = -(2.0_16*c**3*w**6*Math_Dawson(c*Sqrt(w)))/9.0_16
      tmp(99) =  (4.0_16*c**5*w**7*Math_Dawson(c*Sqrt(w)))/45.0_16
      tmp(100) = -(4.0_16*c**7*w**8*Math_Dawson(c*Sqrt(w)))/315.0_16
      tmp(101) =  (8.0_16*c**9*w**9*Math_Dawson(c*Sqrt(w)))/14175.0_16
      tmp(102) = (exp1*cexp12**4*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/(2.0_16*c)
      tmp(103) = (exp1*cexp12**5*Math_Dawson((c*w)/Sqrt(cexp12 +w)))/(10.0_16*c*w)
      tmp(104) = (exp1*cexp12**3*w*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/c
      tmp(105) = -c*exp1*cexp12**4*w*Math_Dawson((c*w)/Sqrt(cexp12 + w)) 
      tmp(106) = (exp1*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(cexp12 +w)))/c
      tmp(107) = - 4.0_16*c*exp1*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(cexp12 + w))
      tmp(108) = (exp1*cexp12*w**3*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/(2.0_16*c) 
      tmp(109) = -6.0_16*c*exp1*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(cexp12 +w)) 
      tmp(110) = (4.0_16*c**3*exp1*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/3.0_16
      tmp(111) =  (exp1*w**4*Math_Dawson((c*w)/Sqrt(cexp12+ w)))/(10.0_16*c) 
      tmp(112) = - 4.0_16*c*exp1*cexp12*w**4*Math_Dawson((c*w)/Sqrt(cexp12 + w))
      tmp(113) = 4.0_16*c**3*exp1*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(cexp12 + w))
      tmp(114) = - c*exp1*w**5*Math_Dawson((c*w)/Sqrt(cexp12 +w)) 
      tmp(115) =  4.0_16*c**3*exp1*cexp12*w**5*Math_Dawson((c*w)/Sqrt(cexp12 + w))
      tmp(116) = -(8.0_16*c**5*exp1*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/15.0_16
      tmp(117) = (4.0_16*c**3*exp1*w**6*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/3.0_16
      tmp(118) = -(16.0_16*c**5*exp1*cexp12*w**6*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/15.0_16
      tmp(119) = -(8.0_16*c**5*exp1*w**7*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/15.0_16
      tmp(120) = (8.0_16*c**7*exp1*cexp12*w**7*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/105.0_16
      tmp(121) = (8.0_16*c**7*exp1*w**8*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/105.0_16
      tmp(122) = -(16.0_16*c**9*exp1*w**9*Math_Dawson((c*w)/Sqrt(cexp12 + w)))/4725.0_16
      tmp(123) = -(20.0_16*exp2*cexp12**4*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/c
      tmp(124) = - (8.0_16*exp2*cexp12**5*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12+ w)))/(c*w)
      tmp(125) = - (20.0_16*exp2*cexp12**3*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/c
      tmp(126) = 40.0_16*c*exp2*cexp12**4*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w))
      tmp(127) = -(10.0_16*exp2*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/c 
      tmp(128) = 80.0_16*c*exp2*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w))
      tmp(129) = -(5.0_16*exp2*cexp12*w**3*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/(2.0_16*c)
      tmp(130) = 60.0_16*c*exp2*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w))
      tmp(131) = -(80.0_16*c**3*exp2*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/3.0_16
      tmp(132) = -(exp2*w**4*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/(4.0_16*c)
      tmp(133) = 20.0_16*c*exp2*cexp12*w**4*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)) 
      tmp(134) = -40.0_16*c**3*exp2*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)) 
      tmp(135) = (5.0_16*c*exp2*w**5*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/2.0_16
      tmp(136) = -20.0_16*c**3*exp2*cexp12*w**5*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w))
      tmp(137) = (16.0_16*c**5*exp2*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/3.0_16
      tmp(138) = -(10.0_16*c**3*exp2*w**6*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/3.0_16
      tmp(139) = (16.0_16*c**5*exp2*cexp12*w**6*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/3.0_16
      tmp(140) = (4.0_16*c**5*exp2*w**7*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/3.0_16
      tmp(141) = -(8.0_16*c**7*exp2*cexp12*w**7*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/21.0_16
      tmp(142) = -(4.0_16*c**7*exp2*w**8*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/21.0_16
      tmp(143) = (8.0_16*c**9*exp2*w**9*Math_Dawson((c*w)/Sqrt(2.0_16*cexp12 + w)))/945.0_16
      tmp(144) = (135.0_16*exp3*cexp12**4*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/c
      tmp(145) = (81.0_16*exp3*cexp12**5*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/(c*w)
      tmp(146) = (90.0_16*exp3*cexp12**3*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/c
      tmp(147) = -270.0_16*c*exp3*cexp12**4*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)) 
      tmp(148) = (30.0_16*exp3*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/c
      tmp(149) = -360.0_16*c*exp3*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(150) = (5.0_16*exp3*cexp12*w**3*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/c 
      tmp(151) = -180.0_16*c*exp3*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(152) = 120.0_16*c**3*exp3*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(153) = (exp3*w**4*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/(3.0_16*c) 
      tmp(154) = -40.0_16*c*exp3*cexp12*w**4*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(155) = 120.0_16*c**3*exp3*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)) 
      tmp(156) = -(10.0_16*c*exp3*w**5*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/3.0_16
      tmp(157) = 40.0_16*c**3*exp3*cexp12*w**5*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(158) = -16.0_16*c**5*exp3*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w))
      tmp(159) = (40.0_16*c**3*exp3*w**6*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/9.0_16
      tmp(160) = -(32.0_16*c**5*exp3*cexp12*w**6*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/3.0_16
      tmp(161) = -(16.0_16*c**5*exp3*w**7*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/9.0_16
      tmp(162) = (16.0_16*c**7*exp3*cexp12*w**7*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/21.0_16
      tmp(163) = (16.0_16*c**7*exp3*w**8*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/63.0_16
      tmp(164) = -(32.0_16*c**9*exp3*w**9*Math_Dawson((c*w)/Sqrt(3.0_16*cexp12 + w)))/2835.0_16
      tmp(165) = -(320.0_16*exp4*cexp12**4*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/c 
      tmp(166) = -(256.0_16*exp4*cexp12**5*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/(c*w) 
      tmp(167) = -(160.0_16*exp4*cexp12**3*w*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/c 
      tmp(168) = 640.0_16*c*exp4*cexp12**4*w*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)) 
      tmp(169) = -(40.0_16*exp4*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/c 
      tmp(170) = 640.0_16*c*exp4*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w))
      tmp(171) = -(5.0_16*exp4*cexp12*w**3*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/c
      tmp(172) = 240.0_16*c*exp4*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)) 
      tmp(173) = -(640.0_16*c**3*exp4*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/3.0_16
      tmp(174) = -(exp4*w**4*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/(4.0_16*c)
      tmp(175) = 40.0_16*c*exp4*cexp12*w**4*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w))
      tmp(176) = -160.0_16*c**3*exp4*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w))
      tmp(177) = (5.0_16*c*exp4*w**5*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/2.0_16
      tmp(178) = -40.0_16*c**3*exp4*cexp12*w**5*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)) 
      tmp(179) =  (64.0_16*c**5*exp4*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/3.0_16
      tmp(180) = -(10.0_16*c**3*exp4*w**6*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/3.0_16
      tmp(181) = (32.0_16*c**5*exp4*cexp12*w**6*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/3.0_16
      tmp(182) = (4.0_16*c**5*exp4*w**7*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/3.0_16
      tmp(183) = -(16.0_16*c**7*exp4*cexp12*w**7*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/21.0_16
      tmp(184) = -(4.0_16*c**7*exp4*w**8*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/21.0_16
      tmp(185) = (8.0_16*c**9*exp4*w**9*Math_Dawson((c*w)/Sqrt(4.0_16*cexp12 + w)))/945.0_16
      tmp(186) = (625.0_16*exp5*cexp12**4*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/(2.0_16*c)
      tmp(187) = (625.0_16*exp5*cexp12**5*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/(2.0_16*c*w)
      tmp(188) = (125.0_16*exp5*cexp12**3*w*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/c 
      tmp(189) = -625.0_16*c*exp5*cexp12**4*w*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)) 
      tmp(190) = (25.0_16*exp5*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/c
      tmp(191) = -500.0_16*c*exp5*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w))
      tmp(192) = (5.0_16*exp5*cexp12*w**3*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/(2.0_16*c) 
      tmp(193) = -150.0_16*c*exp5*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)) 
      tmp(194) = (500.0_16*c**3*exp5*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/3.0_16
      tmp(195) = (exp5*w**4*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/(10.0_16*c) 
      tmp(196) = -20.0_16*c*exp5*cexp12*w**4*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w))
      tmp(197) = 100.0_16*c**3*exp5*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w))
      tmp(198) = -c*exp5*w**5*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)) 
      tmp(199) =  20.0_16*c**3*exp5*cexp12*w**5*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 +w))
      tmp(200) = - (40.0_16*c**5*exp5*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/3.0_16
      tmp(201) = (4.0_16*c**3*exp5*w**6*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/3.0_16
      tmp(202) = -(16.0_16*c**5*exp5*cexp12*w**6*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/3.0_16
      tmp(203) = -(8.0_16*c**5*exp5*w**7*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/15.0_16
      tmp(204) = (8.0_16*c**7*exp5*cexp12*w**7*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/21.0_16
      tmp(205) = (8.0_16*c**7*exp5*w**8*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/105.0_16
      tmp(206) = -(16.0_16*c**9*exp5*w**9*Math_Dawson((c*w)/Sqrt(5.0_16*cexp12 + w)))/4725.0_16
      tmp(207) = -(108.0_16*exp6*cexp12**4*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/c
      tmp(208) = -(648.0_16*exp6*cexp12**5*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/(5.0_16*c*w)
      tmp(209) = -(36.0_16*exp6*cexp12**3*w*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/c
      tmp(210) = 216.0_16*c*exp6*cexp12**4*w*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(211) = -(6.0_16*exp6*cexp12**2*w**2*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/c 
      tmp(212) = 144.0_16*c*exp6*cexp12**3*w**2*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(213) = -(exp6*cexp12*w**3*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/(2.0_16*c)
      tmp(214) = 36.0_16*c*exp6*cexp12**2*w**3*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(215) = -48.0_16*c**3*exp6*cexp12**3*w**3*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(216) = - (exp6*w**4*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12+ w)))/(60.0_16*c)
      tmp(217) =  4.0_16*c*exp6*cexp12*w**4*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(218) = - 24.0_16*c**3*exp6*cexp12**2*w**4*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w))
      tmp(219) =  (c*exp6*w**5*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/6.0_16
      tmp(220) = - 4.0_16*c**3*exp6*cexp12*w**5*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)) 
      tmp(221) =  (16.0_16*c**5*exp6*cexp12**2*w**5*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/5.0_16
      tmp(222) = - (2.0_16*c**3*exp6*w**6*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/9.0_16
      tmp(223) =  (16.0_16*c**5*exp6*cexp12*w**6*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/15.0_16
      tmp(224) =  (4.0_16*c**5*exp6*w**7*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/45.0_16
      tmp(225) = - (8.0_16*c**7*exp6*cexp12*w**7*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/105.0_16
      tmp(226) = - (4.0_16*c**7*exp6*w**8*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/315.0_16
      tmp(227) =  (8.0_16*c**9*exp6*w**9*Math_Dawson((c*w)/Sqrt(6.0_16*cexp12 + w)))/14175.0_16
   
      call sort(tmp, int(size(tmp), 8), int(size(tmp), 8))
    
      auxValue12 = sum_kahan(tmp) * (Math_PI_MOD**1.5)

  end subroutine LJIntegrals_primitive12

!  subroutine LJintegral_computesix(A, B, D, primIntegralValue6)
!
!    implicit none
!    type(ContractedGaussian) :: A, B, D
!    type(pointsorigin) :: pointorigin(:)
!
!    real(8) :: cexp, cp2, cp2_sqrt
!    real(16) :: w, c   
!    real(16) :: exp1, exp2, exp3, tmp(32), primIntegralValue6
!
!    call ContractedGaussian_product(A, B, D)
!
!    call ContractedGaussian_showInSimpleForm(D,6)
!    
!    cexp = 1.0_16 !EXPONENTE DEL CUTOFF
!    
!    !CONVENCIONES: cp2 = CP^2
!    cp2 = sum(D%origin-pointorigin)**2
!    cp2_sqrt = sqrt(cp2)
!        
!    c=real(cp2_sqrt,16)
!    w=real(D%orbitalExponents,16)
!                print*, "CENTRO ATOMO", atomsCenter(:,:)
!
!    exp1=EXP(-(c**2*w) + (c**2*w**2)/(cexp + w))
!    exp2=EXP(-(c**2*w) + (c**2*w**2)/(2.0_16*cexp + w))
!    exp3=EXP(-(c**2*w) + (c**2*w**2)/(3.0_16*cexp + w))
!
!    tmp = 0.0_16
!    tmp(1) = (5.0_16*w**1.5)/3.0_16 
!    tmp(2) = -(2.0_16*c**2*w**2.5)/3.0_16 
!    tmp(3) = -5.0_16*exp1*cexp*Sqrt(cexp + w) 
!    tmp(4) = -5.0_16*exp1*w*Sqrt(cexp + w) 
!    tmp(5) = 2.0_16*c**2*exp1*w**2*Sqrt(cexp + w) 
!    tmp(6) = 10.0_16*exp2*cexp*Sqrt(2.0_16*cexp + w) 
!    tmp(7) = 5.0_16*exp2*w*Sqrt(2.0_16*cexp + w) 
!    tmp(8) = -2.0_16*c**2*exp2*w**2*Sqrt(2.0_16*cexp + w) 
!    tmp(9) = - 5.0_16*exp3*cexp*Sqrt(3.0_16*cexp + w) 
!    tmp(10)= -(5.0_16*exp3*w*Sqrt(3.0_16*cexp + w))/3.0_16 
!    tmp(11) = (2.0_16*c**2*exp3*w**2*Sqrt(3.0_16*cexp + w))/3.0_16
!    tmp(12) = (w*Math_Dawson(c*Sqrt(w)))/c 
!    tmp(13) = -4.0_16*c*w**2*Math_Dawson(c*Sqrt(w)) 
!    tmp(14) = (4.0_16*c**3*w**3*Math_Dawson(c*Sqrt(w)))/3.0_16 
!    tmp(15) = -(6.0_16*exp1*cexp*Math_Dawson((c*w)/Sqrt(cexp + w)))/c 
!    tmp(16) = -(3.0_16*exp1*cexp**2*Math_Dawson((c*w)/Sqrt(cexp + w)))/(c*w) 
!    tmp(17) = -(3.0_16*exp1*w*Math_Dawson((c*w)/Sqrt(cexp + w)))/c 
!    tmp(18) = 12.0_16*c*exp1*cexp*w*Math_Dawson((c*w)/Sqrt(cexp + w)) 
!    tmp(19) = 12.0_16*c*exp1*w**2*Math_Dawson((c*w)/Sqrt(cexp + w)) 
!    tmp(20) = -4.0_16*c**3*exp1*w**3*Math_Dawson((c*w)/Sqrt(cexp + w)) 
!    tmp(21) = (12.0_16*exp2*cexp*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)))/c 
!    tmp(22) = (12.0_16*exp2*cexp**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)))/(c*w) 
!    tmp(23) = (3.0_16*exp2*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)))/c 
!    tmp(24) = -24.0_16*c*exp2*cexp*w*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)) 
!    tmp(25) = -12.0_16*c*exp2*w**2*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)) 
!    tmp(26) = 4.0_16*c**3*exp2*w**3*Math_Dawson((c*w)/Sqrt(2.0_16*cexp + w)) 
!    tmp(27) = -(6.0_16*exp3*cexp*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)))/c 
!    tmp(28) = -(9.0_16*exp3*cexp**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)))/(c*w) 
!    tmp(29) = -(exp3*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)))/c 
!    tmp(30) = 12.0_16*c*exp3*cexp*w*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)) 
!    tmp(31) = 4.0_16*c*exp3*w**2*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)) 
!    tmp(32) = -(4.0_16*c**3*exp3*w**3*Math_Dawson((c*w)/Sqrt(3.0_16*cexp + w)))/3.0_16
!   
!    call sort(tmp, int(size(tmp), 8), int(size(tmp), 8))
!    
!    primIntegralValue6 = sum_kahan(tmp) * (Math_PI_MOD**1.5)
!    
!  end subroutine LJintegral_computesix
!





  !> @brief Obara-Saika recursion for nucleo-electron attraction integrals. Supports all angular momentum numbers
  !! @author E. F. Posada, 2010
  !! @version 1.0  
  subroutine LJIntegrals_obaraSaikaRecursion(AI0, PA, PB, PC, zeta, sumAngularMoment, angularMomentA, angularMomentB)
    implicit none
    real(8), intent(inout), allocatable :: AI0(:,:,:)
    real(8), intent(in) :: PA(0:3)
    real(8), intent(in) :: PB(0:3)
    real(8), intent(in) :: PC(0:3)
    real(8), intent(in) :: zeta
    integer, intent(in) :: sumAngularMoment
    integer, intent(in) :: angularMomentA
    integer, intent(in) :: angularMomentB

    real(8), dimension(0: sumAngularMoment) :: F
    real(8) :: twoZetaInv
    real(8) :: tmp
    real(8) :: u
    integer :: a, b, m
    integer :: izm, iym, ixm
    integer :: jzm, jym, jxm
    integer :: ix,iy,iz,jx,jy,jz
    integer :: iind,jind

    izm = 1
    iym = angularMomentA + 1
    ixm = iym * iym
    jzm = 1
    jym = angularMomentB + 1
    jxm = jym * jym
    twoZetaInv = 1 / (2 * zeta)
    tmp = sqrt(zeta)*(2.0_8/Math_SQRT_PI)
    u = zeta*(PC(0) * PC(0) + PC(1) * PC(1) + PC(2) * PC(2))

    call Math_fgamma0(sumAngularMoment,u,F)

    do m = 0, sumAngularMoment
       AI0(0, 0, m) = tmp * F(m)
    end do

    !! Upward recursion in j with i=0
    do b = 1, angularMomentB
       do jx = 0, b
          do jy=0, b - jx
             jz = b - jx - jy
             jind = jx * jxm + jy * jym + jz * jzm
             if (jz > 0) then
                do m=0, sumAngularMoment - b	!! Electrostatic potential integrals
                   AI0(0,jind,m) = PB(2)*AI0(0,jind-jzm,m) - PC(2)*AI0(0, jind-jzm, m+1)
                end do

                if (jz > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jz-1)*(AI0(0,jind-2*jzm,m) - AI0(0,jind-2*jzm,m+1))
                   end do
                end if

             else if (jy > 0) then
                do m=0, sumAngularMoment-b
                   AI0(0,jind,m) = PB(1)*AI0(0,jind-jym,m) -	PC(1)*AI0(0,jind-jym,m+1)
                end do

                if (jy > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jy-1)*(AI0(0,jind-2*jym,m) - AI0(0,jind-2*jym,m+1))
                   end do
                end if

             else if (jx > 0) then
                do m=0, sumAngularMoment-b
                   AI0(0,jind,m) = PB(0)*AI0(0,jind-jxm,m) -	PC(0)*AI0(0,jind-jxm,m+1)
                end do

                if (jx > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jx-1)*(AI0(0,jind-2*jxm,m) - AI0(0,jind-2*jxm,m+1))
                   end do
                end if

             else
                call LJIntegrals_exception( ERROR, "LJIntegrals in obaraSaikaRecursion Function", &
                     "There's some error in the obaraSaika algorithm")
             end if
          end do
       end do
    end do

    !! The following fragment cannot be vectorized easily, I guess :-)
    !! Upward recursion in i with all possible j's
    do b=0, angularMomentB
       do jx=0, b
          do jy=0, b-jx
             jz = b-jx-jy
             jind = jx*jxm + jy*jym + jz*jzm
             do a=1, angularMomentA
                do ix=0, a
                   do iy=0, a-ix
                      iz = a-ix-iy
                      iind = ix*ixm + iy*iym + iz*izm
                      if (iz > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(2)*AI0(iind-izm,jind,m) - PC(2)*AI0(iind-izm,jind,m+1)
                         end do

                         if (iz > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(iz-1)*(AI0(iind-2*izm,jind,m) - AI0(iind-2*izm,jind,m+1))
                            end do
                         end if

                         if (jz > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jz*(AI0(iind-izm,jind-jzm,m) - AI0(iind-izm,jind-jzm,m+1))
                            end do
                         end if

                      else if (iy > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(1)*AI0(iind-iym,jind,m) - PC(1)*AI0(iind-iym,jind,m+1)
                         end do

                         if (iy > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(iy-1)*(AI0(iind-2*iym,jind,m) - AI0(iind-2*iym,jind,m+1))
                            end do
                         end if

                         if (jy > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jy*(AI0(iind-iym,jind-jym,m) - AI0(iind-iym,jind-jym,m+1))
                            end do
                         end if

                      else if (ix > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(0)*AI0(iind-ixm,jind,m) - PC(0)*AI0(iind-ixm,jind,m+1)
                         end do

                         if (ix > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(ix-1)*(AI0(iind-2*ixm,jind,m) - AI0(iind-2*ixm,jind,m+1))
                            end do
                         end if

                         if (jx > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jx*(AI0(iind-ixm,jind-jxm,m) - AI0(iind-ixm,jind-jxm,m+1))
                            end do
                         end if
                      else
                         call LJIntegrals_exception( ERROR, "LJIntegrals in obaraSaikaRecursion Function", &
                              "There's some error in the obaraSaika algorithm")
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine LJIntegrals_obaraSaikaRecursion

  !>
  !! @brief  Handle exceptions of this module
  subroutine LJIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine LJIntegrals_exception

end module LJIntegrals_
