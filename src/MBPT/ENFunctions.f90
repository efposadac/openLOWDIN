!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Epstein-Nesbet and APMO-Epstein-Nesbet program.
!!        This module allows to make calculations in the APMO-Epstein-Nesbet framework
!! @author  
!!
!! <b> Creation date : </b> 2019-06-19
!!
!! <b> History: </b>
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module ENFunctions_
  use CONTROL_
  use MolecularSystem_
  use IndexMap_
  use Exception_
  use Vector_
  use ReadTransformedIntegrals_
  use String_
  implicit none

  !< enum EpsteinNesbet_correctionFlags {
  integer, parameter :: FIRST_ORDER = 1
  integer, parameter :: SECOND_ORDER = 2
  integer, parameter :: THIRD_ORDER = 3
  !< }

  type :: EpsteinNesbetTP
    
        character(50) :: name
        real(8) :: energyHF
        integer :: orderOfCorrection
        integer :: numberOfSpecies
        integer :: frozenCoreBoundary
        real(8) :: totalEnergy(3)
        real(8) :: totalCorrection(3)
        real(8) :: secondOrderCorrection(3)
        real(8) :: thirdOrderCorrection
        !! Vectores para almacenamiento de correcciones a la energia.para cada especie
        type(Matrix) :: energyCorrectionOfSecondOrder
        type(Matrix) :: energyOfCouplingCorrectionOfSecondOrder

        logical :: isInstanced

  end type EpsteinNesbetTP

       type(EpsteinNesbetTP), target :: EpsteinNesbet_instance
!       character(50) :: job
       private :: &
            EpsteinNesbet_secondOrderCorrection, &
            EpsteinNesbet_thirdOrderCorrection

       public :: &
            EpsteinNesbet_constructor, &
            EpsteinNesbet_destructor, &
            EpsteinNesbet_show, &
            EpsteinNesbet_run, &
            EpsteinNesbet_getTotalEnergy, &
            EpsteinNesbet_getEnergyCorrection, &
            EpsteinNesbet_getSpecieCorrection


contains
  !**
  ! Define el constructor para la clase
  !
  !**
  subroutine EpsteinNesbet_constructor( orderOfCorrection )
    implicit none
    integer, intent(in) :: orderOfCorrection

    integer :: i
    type(Exception) :: ex

!   if( .not.EpsteinNesbet_instance%isInstanced ) then

      EpsteinNesbet_instance%orderOfCorrection = orderOfCorrection
      EpsteinNesbet_instance%numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


      if ( EpsteinNesbet_instance%orderOfCorrection >= 2 ) then
        !!call Vector_constructor( EpsteinNesbet_instance%energyCorrectionOfSecondOrder, EpsteinNesbet_instance%numberOfSpecies)
        call Matrix_constructor( EpsteinNesbet_instance%energyCorrectionOfSecondOrder, int(EpsteinNesbet_instance%numberOfSpecies,8), 3_8)

        i = EpsteinNesbet_instance%numberOfSpecies * ( EpsteinNesbet_instance%numberOfSpecies-1 ) / 2



        !!call Vector_constructor( EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder, i)
        call Matrix_constructor( EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder, int(i,8), 3_8)
      end if

      if ( EpsteinNesbet_instance%orderOfCorrection >= 3 ) then

        call Exception_constructor( ex , ERROR )
        call Exception_setDebugDescription( ex, "Class object EpsteinNesbet in the constructor() function" )
        call Exception_setDescription( ex, "This order correction hasn't been implemented" )
        call Exception_show( ex )

      end if

      EpsteinNesbet_instance%isInstanced =.true.

!   end if

  end subroutine EpsteinNesbet_constructor

  !**
  ! Define el destructor para clase
  !
  !**
  subroutine EpsteinNesbet_destructor()
    implicit none

               EpsteinNesbet_instance%isInstanced =.false.

  end subroutine EpsteinNesbet_destructor


  !**
  ! @brief Muestra informacion asociada con la correccion MPn
  !**
  subroutine EpsteinNesbet_show()
    implicit none

    integer :: i
    integer :: j
    integer :: k

    if ( EpsteinNesbet_instance%isInstanced )  then

        print *,""
        print *," POST HARTREE-FOCK CALCULATION"
        print *," MANY-BODY PERTURBATION THEORY:"
        print *,"=============================="
        print *,""

        write (6,"(T10,A25)") "MOLLER-PLESSET FORMALISM "
        write (6,"(T10,A23, I5)") "ORDER OF CORRECTION = ",EpsteinNesbet_instance%orderOfCorrection

        print *,""
        write (6,'(T10,A15,ES20.12)') "E(0) + E(1) = ", EpsteinNesbet_instance%energyHF
        write (6,'(T10,A15,ES20.12)') "E(2) = ", EpsteinNesbet_instance%secondOrderCorrection(1)
        write (6,'(T25,A20)') "________________________"
        write (6,'(T10,A15,ES25.17)') "E(MP2)= ", EpsteinNesbet_instance%totalEnergy(1)
        print *,""
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        write ( 6,'(T10,A15,A20)') " E(n){ Specie } ","   E(n) / Hartree "
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        print *,""

        do i=1, EpsteinNesbet_instance%numberOfSpecies
       write (*,'(T10,A5,A8,A8,ES16.8)') "E(2){", trim(  MolecularSystem_getNameOfSpecie( i ) ),"   } = ", &
               EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(i,1)
        end do

        print *,""
        k=0
        do i=1, EpsteinNesbet_instance%numberOfSpecies
       do j=i+1,EpsteinNesbet_instance%numberOfSpecies
           k=k+1
           write (*,'(T10,A5,A16,A4,ES16.8)') "E(2){", &
               trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
               "} = ", EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder%values(k,1)
       end do
     end do

        print *,""
        write (6,"(T10,A37)") "NON-SINGULAR EPSTEIN-NESBET FORMALISM "
        write (6,"(T10,A23, I5)") "ORDER OF CORRECTION = ",EpsteinNesbet_instance%orderOfCorrection

        print *,""
        write (6,'(T10,A15,ES20.12)') "E(0) + E(1) = ", EpsteinNesbet_instance%energyHF
        write (6,'(T10,A15,ES20.12)') "E(2) = ", EpsteinNesbet_instance%secondOrderCorrection(2)
        write (6,'(T25,A20)') "________________________"
        write (6,'(T10,A15,ES25.17)') "E(NS-EN2)= ", EpsteinNesbet_instance%totalEnergy(2)
        print *,""
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        write ( 6,'(T10,A15,A20)') " E(n){ Specie } ","   E(n) / Hartree "
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        print *,""

        do i=1, EpsteinNesbet_instance%numberOfSpecies
       write (*,'(T10,A5,A8,A8,ES16.8)') "E(2){", trim(  MolecularSystem_getNameOfSpecie( i ) ),"   } = ", &
               EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(i,2)
        end do

        print *,""
        k=0
        do i=1, EpsteinNesbet_instance%numberOfSpecies
       do j=i+1,EpsteinNesbet_instance%numberOfSpecies
           k=k+1
           write (*,'(T10,A5,A16,A4,ES16.8)') "E(2){", &
               trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
               "} = ", EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder%values(k,2)
       end do
     end do

        print *,""
        write (6,"(T10,A25)") "EPSTEIN-NESBET FORMALISM "
        write (6,"(T10,A23, I5)") "ORDER OF CORRECTION = ",EpsteinNesbet_instance%orderOfCorrection

        print *,""
        write (6,'(T10,A15,ES20.12)') "E(0) + E(1) = ", EpsteinNesbet_instance%energyHF
        write (6,'(T10,A15,ES20.12)') "E(2) = ", EpsteinNesbet_instance%secondOrderCorrection(3)
        write (6,'(T25,A20)') "________________________"
        write (6,'(T10,A15,ES25.17)') "E(EN2)= ", EpsteinNesbet_instance%totalEnergy(3)
        print *,""
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        write ( 6,'(T10,A15,A20)') " E(n){ Specie } ","   E(n) / Hartree "
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        print *,""

        do i=1, EpsteinNesbet_instance%numberOfSpecies
       write (*,'(T10,A5,A8,A8,ES16.8)') "E(2){", trim(  MolecularSystem_getNameOfSpecie( i ) ),"   } = ", &
               EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(i,3)
        end do

        print *,""
        k=0
        do i=1, EpsteinNesbet_instance%numberOfSpecies
       do j=i+1,EpsteinNesbet_instance%numberOfSpecies
           k=k+1
           write (*,'(T10,A5,A16,A4,ES16.8)') "E(2){", &
               trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
               "} = ", EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder%values(k,3)
       end do
     end do



    end if

  end subroutine EpsteinNesbet_show

  !**
  ! @brief realiza la correccion de energia MPn orden dado
  !**
  subroutine EpsteinNesbet_run()
    implicit none

    type(Exception) :: ex
    integer :: i

  character(50) :: wfnFile
  integer :: wfnUnit
  wfnFile = "lowdin.wfn"
  wfnUnit = 20

  
  if ( EpsteinNesbet_instance%isInstanced ) then

     do i=2, EpsteinNesbet_instance%orderOfCorrection

        select case( i )

        case(  SECOND_ORDER )

           call EpsteinNesbet_secondOrderCorrection()

           EpsteinNesbet_instance%totalCorrection(1) = EpsteinNesbet_instance%secondOrderCorrection(1)
           EpsteinNesbet_instance%totalCorrection(2) = EpsteinNesbet_instance%secondOrderCorrection(2)
           EpsteinNesbet_instance%totalCorrection(3) = EpsteinNesbet_instance%secondOrderCorrection(3)

        case(  THIRD_ORDER )

           call EpsteinNesbet_thirdOrderCorrection()

           EpsteinNesbet_instance%totalCorrection(1) = EpsteinNesbet_instance%totalCorrection(1) + &
                EpsteinNesbet_instance%thirdOrderCorrection

        case default

           call Exception_constructor( ex , ERROR )
           call Exception_setDebugDescription( ex, "Class object EpsteinNesbet in run() function" )
           call Exception_setDescription( ex, "This order correction hasn't been implemented" )
           call Exception_show( ex )

        end select

     end do
  !! Open file for wavefunction
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  !! Load results...

   call Vector_getFromFile(unit=wfnUnit, binary=.true., value=EpsteinNesbet_instance%energyHF, arguments=["TOTALENERGY"])

      EpsteinNesbet_instance%totalEnergy(1) = EpsteinNesbet_instance%energyHF + EpsteinNesbet_instance%totalCorrection(1)
      EpsteinNesbet_instance%totalEnergy(2) = EpsteinNesbet_instance%energyHF + EpsteinNesbet_instance%totalCorrection(2)
      EpsteinNesbet_instance%totalEnergy(3) = EpsteinNesbet_instance%energyHF + EpsteinNesbet_instance%totalCorrection(3)

close(wfnUnit)

    else

      call Exception_constructor( ex , ERROR )
      call Exception_setDebugDescription( ex, "Class object EpsteinNesbet in run() function" )
      call Exception_setDescription( ex, "You should to instance EpsteinNesbet module before use this function" )
      call Exception_show( ex )

end if


  end subroutine EpsteinNesbet_run


  !**
  ! @ Retorna el la correccion a la energia
  !**
  function EpsteinNesbet_getEnergyCorrection() result( output)
    implicit none
    real(8) :: output

    output = EpsteinNesbet_instance%totalCorrection(1) !!MP2

  end function EpsteinNesbet_getEnergyCorrection

  !**
  ! @ Retorna el la correccion a la energia
  !**
  function EpsteinNesbet_getSpecieCorrection( specieName) result( output)
    implicit none
    character(*) :: specieName
    real(8) :: output

    integer :: i

    do i=1, EpsteinNesbet_instance%numberOfSpecies

      if ( trim(specieName) == trim(  MolecularSystem_getNameOfSpecie( i ) ) ) then

        output = EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(i,1) !!MP2
        return

      end if

    end do

  end function EpsteinNesbet_getSpecieCorrection



  !**
  ! @ Retorna la energia final con correccion  de orden dado
  !**
  function EpsteinNesbet_getTotalEnergy() result(output)
    implicit none
    real(8) :: output

    output = EpsteinNesbet_instance%totalEnergy(1) !!MP2

  end function EpsteinNesbet_getTotalEnergy


 !<
 !! @brief Realiza la correccion de segundo orden a la enegia
 !>
 subroutine EpsteinNesbet_secondOrderCorrection()
   implicit none
   
   
   integer :: a
   integer :: b
   integer :: r
   integer :: s
   integer :: p
   integer :: t
   integer :: i
   integer :: j
   integer :: m
   integer :: k
   integer :: is, js
   integer :: specieID
   integer :: otherSpecieID
   integer :: electronsID
   integer :: ocupationNumber
   integer :: ocupationNumberOfOtherSpecie
   integer :: numberOfContractions
   integer :: numberOfContractionsOfOtherSpecie
   integer(8) :: auxIndex
    character(10) :: order
   character(10) :: nameOfSpecie
   character(10) :: nameOfOtherSpecie
   type(Vector) :: eigenValues
   type(Vector) :: eigenValuesOfOtherSpecie
   type(Matrix) :: auxMatrix, auxMatrixA, auxMatrixB
   integer :: occupation, otherOccupation
!   type(TransformIntegrals) :: repulsionTransformer
   real(8) :: lambda
   real(8) :: lambdaOfOtherSpecie
   real(8) :: independentEnergyCorrection(3), E1, E2
   real(8) :: couplingEnergyCorrection(3)
   real(8) :: deltaE, deltaD, couplingEnergyCorrection2, auxfactor, totalEN
   real(8) :: auxVal_A
   real(8) :: auxVal_B
   real(8) :: auxVal_C(12)
   type(Matrix) :: eigenVec
   type(Matrix) :: eigenVecOtherSpecie 
   character(50) :: wfnFile
   character(50) :: arguments(2)
   integer :: wfnUnit

   integer :: i1, i2, j1, j2, nVirtualOrbitals
   integer :: ii(CONTROL_instance%INTEGRAL_STACK_SIZE)
   integer :: jj(CONTROL_instance%INTEGRAL_STACK_SIZE)
   integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
   integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
   real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
   integer :: integralStackSize
   character(255) :: prefixOfFile
   integer :: unidOfOutputForIntegrals
   real(8), allocatable :: intArray(:,:,:)
   real(8) :: auxIntegral
   integer :: errorValue

   !! TypeOfIntegrals 
   integer, parameter :: ONE_SPECIE     = 0
   integer, parameter :: TWO_SPECIES      = 1

   wfnFile = "lowdin.wfn"
   wfnUnit = 20
   

   electronsID = MolecularSystem_getSpecieID(  nameOfSpecie="E-" )

   !! Define file names
!   call TransformIntegrals_constructor( repulsionTransformer )

   !!*******************************************************************************************
   !! Calculo de correcciones de segundo orden para particulas de la misma especie
   !!

   open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 
   rewind(wfnUnit)
   totalEN = 0

   do is=1, EpsteinNesbet_instance%numberOfSpecies

      
      EpsteinNesbet_instance%frozenCoreBoundary = 1

      if ( is == electronsID )  EpsteinNesbet_instance%frozenCoreBoundary = &
           CONTROL_instance%MP_FROZEN_CORE_BOUNDARY
      
      nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( is ) )
      
      independentEnergyCorrection = 0.0_8
      
      if( trim(nameOfSpecie)=="E-" .or. .not.CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION ) then
         
         numberOfContractions = MolecularSystem_getTotalNumberOfContractions(is)
         arguments(2) = MolecularSystem_getNameOfSpecie(is)

         arguments(1) = "COEFFICIENTS"
         eigenVec= Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
              columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
   
         arguments(1) = "ORBITALS"
         call Vector_getFromFile( elementsNum = numberOfContractions, &
              unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
              output = eigenValues )     
         
         
         specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
         ocupationNumber = MolecularSystem_getOcupationNumber( is )
         numberOfContractions = MolecularSystem_getTotalNumberOfContractions( is )
         lambda = MolecularSystem_instance%species(is)%lambda
         
         !! Read transformed integrals from file
         !if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

           call ReadTransformedIntegrals_readOneSpecies( specieID, auxMatrix)

           couplingEnergyCorrection2 = 0
           do i=EpsteinNesbet_instance%frozenCoreBoundary, ocupationNumber
             do j=EpsteinNesbet_instance%frozenCoreBoundary,ocupationNumber
               do a=ocupationNumber+1, numberOfContractions
                 do b=ocupationNumber+1, numberOfContractions
                 !do b=a, numberOfContractions
                   auxVal_C = 0

                   auxIndex = IndexMap_tensorR4ToVectorB(int(i,8),int(a,8),int(j,8),int(b,8), &
                                                        int(numberOfContractions,8) )

                   auxVal_A= auxMatrix%values(auxIndex, 1)

                   !! iijj
                   auxIndex = IndexMap_tensorR4ToVectorB(int(i,8),int(i,8),int(j,8),int(j,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(1) = - auxMatrix%values(auxIndex, 1)
                   !if ( b == a ) then
                   !if ( MolecularSystem_getNumberOfParticles(is) > lambda ) then
                   !if ( i == j ) then
                   if ( i /= j ) then
                   auxIndex = IndexMap_tensorR4ToVectorB(int(i,8),int(j,8),int(j,8),int(i,8), &
                                                       int(numberOfContractions,8) )
                   auxVal_C(1+6) =  auxMatrix%values(auxIndex, 1) !add if
                   end if
                   !end if
 
                   !! aabb
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(a,8),int(b,8),int(b,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(2) = - auxMatrix%values(auxIndex, 1)
                   !if ( MolecularSystem_getNumberOfParticles(is) > lambda ) then
                   !if ( a /= b ) then
                   if ( a /= b ) then
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(b,8),int(b,8),int(a,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(2+6) = auxMatrix%values(auxIndex, 1) !! add iff
                   end if
                   !end if

                   !! aaii
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(a,8),int(i,8),int(i,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(3) = auxMatrix%values(auxIndex, 1)
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(i,8),int(i,8),int(a,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(3+6) = - auxMatrix%values(auxIndex, 1)

                   !! aajj
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(a,8),int(j,8),int(j,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(4) = auxMatrix%values(auxIndex, 1)
                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(j,8),int(j,8),int(a,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(4+6) = - auxMatrix%values(auxIndex, 1) !! add if

                   !! bbii
                   auxIndex = IndexMap_tensorR4ToVectorB(int(b,8),int(b,8),int(i,8),int(i,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(5) = auxMatrix%values(auxIndex, 1)
                   auxIndex = IndexMap_tensorR4ToVectorB(int(b,8),int(i,8),int(i,8),int(b,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(5+6) = - auxMatrix%values(auxIndex, 1) !! add if

                   !! bbjj
                   auxIndex = IndexMap_tensorR4ToVectorB(int(b,8),int(b,8),int(j,8),int(j,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(6) = auxMatrix%values(auxIndex, 1)
                   auxIndex = IndexMap_tensorR4ToVectorB(int(b,8),int(j,8),int(j,8),int(b,8), &
                                                        int(numberOfContractions,8) )
                   auxVal_C(6+6) = - auxMatrix%values(auxIndex, 1)

                    deltaE = (eigenValues%values(i) + eigenValues%values(j) &
                         -eigenValues%values(a)-eigenValues%values(b))
                    deltaD = sum(auxVal_C(1:6)) + auxVal_C(3+6) + auxVal_C(6+6)


                    if ( ( deltaE + deltaD ) < 0 ) then
                      auxfactor = -1.0
                    else 
                      auxfactor = 1.0
                    end if

                   if (  dabs( auxVal_A)  > 1.0E-20_8 ) then

                   !write(*,"(A2, I3, I3, I3, I3, F10.6, 12F10.6)")  "AA", i, a, j, b, auxVal_A, auxVal_C 

                     if ( b /= a ) then

                       if (i==j) then

                         if( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                           independentEnergyCorrection(1) = independentEnergyCorrection(1) + auxVal_A**2.0  &
                                                        * ( lambda  -  1.0_8 ) / ( deltaE )

                           !deltaD = sum(auxVal_C) - auxVal_C(1+6)

                           independentEnergyCorrection(2) = independentEnergyCorrection(2) + &
                           0.5 * ( auxfactor * sqrt (  4.0*auxVal_A **2.0_8*(lambda - 1.0_8 ) + ( deltaE + deltaD)**2.0_8 ) - &
                           ( deltaE + deltaD ) )

                           independentEnergyCorrection(3) = independentEnergyCorrection(3) + auxVal_A**2.0  &
                                                        * ( lambda  -  1.0_8 ) / ( deltaE + deltaD  )

                         end if

                       else

                         deltaD = deltaD + auxVal_C(1+6)
                         deltaD = deltaD + auxVal_C(2+6)
                         deltaD = deltaD + auxVal_C(4+6)
                         deltaD = deltaD + auxVal_C(5+6)
 
                         if ( ( deltaE + deltaD ) < 0 ) then
                           auxfactor = -1.0
                         else 
                           auxfactor = 1.0
                         end if

                         auxIndex = IndexMap_tensorR4ToVector(a, j, b, i, numberOfContractions )
                         auxVal_B= auxMatrix%values(auxIndex, 1)

                         independentEnergyCorrection(1) = independentEnergyCorrection(1) + auxVal_A  &
                                                     * ( lambda * auxVal_A  - auxVal_B ) / ( deltaE )

                         !deltaD = sum(auxVal_C) 

                         independentEnergyCorrection(2) = independentEnergyCorrection(2) + &
                         0.5 * ( auxfactor * sqrt ( ( 4.0*auxVal_A * (lambda *auxVal_A - auxVal_B )) + ( deltaE + deltaD)**2.0_8 ) - &
                         ( deltaE + deltaD ) )

                         independentEnergyCorrection(3) = independentEnergyCorrection(3) + auxVal_A  &
                                                     * ( lambda * auxVal_A  - auxVal_B ) / ( deltaE + deltaD )


                    !deltaD = sum(auxVal_C(1:6)) + auxVal_C(3+6) + auxVal_C(6+6)
                       end if

                     else if ( i==j .and. a==b ) then

                       if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                         independentEnergyCorrection(1) = independentEnergyCorrection(1) +  auxVal_A**2.0_8  &
                                                      * ( lambda - 1.0_8 ) /  ( deltaE )

                         !deltaD = sum(auxVal_C) - auxVal_C(1+6) - auxVal_C(2+6)

                         independentEnergyCorrection(2) = independentEnergyCorrection(2) + &
                         0.5 * ( auxfactor * sqrt ( 4.0*( auxVal_A )**2.0_8*(lambda - 1.0_8 ) + ( deltaE + deltaD)**2.0_8 ) - &
                         ( deltaE + deltaD ) )

                         independentEnergyCorrection(3) = independentEnergyCorrection(3) +  auxVal_A**2.0_8  &
                                                      * ( lambda - 1.0_8 ) /  ( deltaE + deltaD )

                       end if

                     else

                       if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                         independentEnergyCorrection(1) = independentEnergyCorrection(1) + auxVal_A**2.0  &
                                                      * ( lambda  - 1.0_8 ) / ( deltaE )

                         !if ( b == a ) then
                         !  deltaD = sum(auxVal_C) - auxVal_C(2+6)
                         !else 
                         !  deltaD = sum(auxVal_C)
                         !end if

                         independentEnergyCorrection(2) = independentEnergyCorrection(2) + &
                         0.5 * ( auxfactor * sqrt ( 4.0*( auxVal_A )**2.0_8*(lambda - 1.0_8 ) + ( deltaE + deltaD)**2.0_8 ) - &
                         ( deltaE + deltaD ) )

                         independentEnergyCorrection(3) = independentEnergyCorrection(3) + auxVal_A**2.0  &
                                                      * ( lambda  - 1.0_8 ) / ( deltaE + deltaD )

                       end if

                     end if
                   end if

                 end do
               end do
             end do
           end do

                     
       !end if        
     end if        
                     
    do k = 1, 3
      EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(is,k) = independentEnergyCorrection(k) &
        * ( ( MolecularSystem_getCharge( specieID ) )**4.0_8 )
                     
      if ( nameOfSpecie == "E-ALPHA" .or. nameOfSpecie == "E-BETA" ) then
                       
        EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(is,k) = &
          EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(is,k) / ( 2.0 )
        couplingEnergyCorrection2 = couplingEnergyCorrection2 / (2.0)

      else             
                       
        EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(is,k) = &
          EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(is,k) / &
          ( MolecularSystem_getParticlesFraction ( specieID ) * 2.0 )

        couplingEnergyCorrection2 = couplingEnergyCorrection2 / &
          ( MolecularSystem_getParticlesFraction ( specieID ) * 2.0 )

      end if           
    end do !k

                     
                     
    call Matrix_destructor(auxMatrix)
    !!
    !!**************************************************************************

  end do


  !! Suma las correcciones de energia para especies independientes
  do k = 1, 3
    EpsteinNesbet_instance%secondOrderCorrection(k) = sum( EpsteinNesbet_instance%energyCorrectionOfSecondOrder%values(:,k) )
  end do

  !!
  !!*******************************************************************************************


  !!*******************************************************************************************
  !! Calculo de correcion se segundo orden para interaccion entre particulas de especie diferente
  !!
  if ( EpsteinNesbet_instance%numberOfSpecies > 1 ) then

   couplingEnergyCorrection = 0.0_8
   m=0

   do is = 1 , EpsteinNesbet_instance%numberOfSpecies

     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(is)
     arguments(2) = trim(MolecularSystem_getNameOfSpecie(is))

     arguments(1) = "COEFFICIENTS"
     eigenVec = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "ORBITALS"
     call Vector_getFromFile( elementsNum = numberOfContractions, &
          unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = eigenValues )     

     nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( is ) )
     specieID =MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecie) )
     ocupationNumber = MolecularSystem_getOcupationNumber( is )
     numberOfContractions = MolecularSystem_getTotalNumberOfContractions( is )
     lambda = MolecularSystem_getEta( is )

     do js = is + 1 , EpsteinNesbet_instance%numberOfSpecies
        m = m + 1


        numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( js )

        arguments(2) = trim(MolecularSystem_getNameOfSpecie(js))

        arguments(1) = "COEFFICIENTS"
        eigenVecOtherSpecie = &
             Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractionsOfOtherSpecie,4), &
             columns= int(numberOfContractionsOfOtherSpecie,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = numberOfContractionsofOtherSpecie, &
             unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
             output = eigenValuesOfOtherSpecie )     


        nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( js ) )
        otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
        ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( js )

        lambdaOfOtherSpecie = MolecularSystem_instance%species(js)%lambda
        couplingEnergyCorrection = 0.0_8

!        if ( .not.CONTROL_instance%OPTIMIZE ) then
!           write (6,"(T10,A)") "INTER-SPECIES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)//"/"//trim(nameOfOtherSpecie)
!           print *,""
!        end if

         !if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

         !! Read transformed integrals from file

        call ReadTransformedIntegrals_readOneSpecies( specieID, auxMatrixA)
        call ReadTransformedIntegrals_readOneSpecies( otherSpecieID, auxMatrixB)
        call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, auxMatrix)

!        call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!             eigenVec, eigenVecOtherSpecie, &
!             auxMatrix, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )

        auxMatrix%values = auxMatrix%values * MolecularSystem_getCharge( specieID ) &
             * MolecularSystem_getCharge( otherSpecieID )
        auxMatrixA%values = auxMatrixA%values * MolecularSystem_getCharge( specieID ) 
        auxMatrixB%values = auxMatrixB%values * MolecularSystem_getCharge( otherSpecieID ) 

        auxVal_C = 0

        !print *, "iiiiiiiiiiiiii"

        do i=1, ocupationNumber
           do j=1,ocupationNumberOfOtherSpecie
              do a=ocupationNumber+1, numberOfContractions
                 do b=ocupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie

                    ! ii II 
                    auxIndex = IndexMap_tensorR4ToVector(i,i,j,j, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                    auxVal_C(1) = -auxMatrix%values(auxIndex,1) 

                    ! aa AA
                    auxIndex = IndexMap_tensorR4ToVector(a,a,b,b, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                    auxVal_C(2) = -auxMatrix%values(auxIndex,1) 

                    if ( MolecularSystem_instance%species(is)%isElectron .or. &
                         CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE .or. &
                         MolecularSystem_getNumberOfParticles(is) > 1 ) then

                    ! aaii
                    auxIndex = IndexMap_tensorR4ToVector(i,i,a,a, numberOfContractions )
                    auxVal_C(3) = -auxMatrixA%values(auxIndex,1) 

                    auxIndex = IndexMap_tensorR4ToVector(i,a,a,i, numberOfContractions )
                    auxVal_C(3) = auxVal_C(3) + auxMatrixA%values(auxIndex,1) 
                    end if

                    ! aa II
                    auxIndex = IndexMap_tensorR4ToVector(a,a,j,j, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                    auxVal_C(4) = auxMatrix%values(auxIndex,1) 
 
                    ! ii AA
                    auxIndex = IndexMap_tensorR4ToVector(i,i,b,b, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                    auxVal_C(5) = auxMatrix%values(auxIndex,1) 
 
                    if ( MolecularSystem_instance%species(js)%isElectron .or. &
                         CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE .or. &
                         MolecularSystem_getNumberOfParticles(js) > 1 ) then
                    ! IIAA
                    auxIndex = IndexMap_tensorR4ToVector(j,j,b,b, numberOfContractionsOfOtherSpecie )
                    auxVal_C(6) = -auxMatrixB%values(auxIndex,1) 

                    auxIndex = IndexMap_tensorR4ToVector(j,b,b,j, numberOfContractionsOfOtherSpecie )
                    auxVal_C(6) = auxVal_C(6) + auxMatrixB%values(auxIndex,1) 
                    end if


                    auxIndex = IndexMap_tensorR4ToVector(i,a,j,b, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                   if (  dabs( auxMatrix%values(auxIndex,1))  > 1.0E-20_8 ) then

                   !write(*,"(A2, I3, I3, I3, I3, F10.6, 12F10.6)")  "AB", i, a, j, b, auxMatrix%values(auxIndex,1), auxVal_C 
                   !write(*,"(A2, I3, I3, I3, I3, F14.10, F14.10)")  "AB", i, a, j, b, abs(auxMatrix%values(auxIndex,1)), abs(sum(auxVal_C ))

                    deltaE = (eigenValues%values(i) + eigenValuesOfOtherSpecie%values(j) &
                         -eigenValues%values(a)-eigenValuesOfOtherSpecie%values(b))
                    deltaD = sum(auxVal_C)

                    if ( ( deltaE + deltaD ) < 0 ) then
                      auxfactor = -1.0
                    else 
                      auxfactor = 1.0
                    end if

                    couplingEnergyCorrection(1) = couplingEnergyCorrection(1) +  &
                         ( ( auxMatrix%values(auxIndex,1) )**2.0_8 ) / ( deltaE )

                    couplingEnergyCorrection(2) = couplingEnergyCorrection(2) +  &
                         0.5 * ( auxfactor * sqrt ( ( 2.0*auxMatrix%values(auxIndex,1) )**2.0_8 + ( deltaE + deltaD)**2.0_8 ) - &
                         ( deltaE + deltaD ) )

                    couplingEnergyCorrection(3) = couplingEnergyCorrection(3) +  &
                         ( ( auxMatrix%values(auxIndex,1) )**2.0_8 ) / ( deltaE + deltaD )

                   !write(*,"(A2, I3, I3, I3, I3, F14.10, F14.10, F14.10)")  "AB", i, a, j, b, abs(auxMatrix%values(auxIndex,1)), &
                   !        abs(sum(auxVal_C )), couplingEnergyCorrection(3)

                   end if 

                 end do
              end do
           end do
        end do

      !end if

        do k = 1, 3
          EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder%values(m,k)= &
               ( lambda * lambdaOfOtherSpecie * couplingEnergyCorrection(k) )  

        end do 
     end do
  end do


    call Matrix_destructor(auxMatrix)
    !! Adiciona la correccion del termino de acoplamiento
    do k = 1, 3
      EpsteinNesbet_instance%secondOrderCorrection(k) = EpsteinNesbet_instance%secondOrderCorrection(k) + &
       sum( EpsteinNesbet_instance%energyOfCouplingCorrectionOfSecondOrder%values(:,k) )
    end do
  end if
!!
!!*******************************************************************************************

close(wfnUnit)

end subroutine EpsteinNesbet_secondOrderCorrection

!<
!! @brief Realiza la correccion de tercer orden a la enegia
!>
subroutine EpsteinNesbet_thirdOrderCorrection()
  implicit none

  call EpsteinNesbet_exception( ERROR, "Class object EpsteinNesbet in thirdOrderCorrection() function" , &
       "This order correction hasn't been implemented" )

end subroutine EpsteinNesbet_thirdOrderCorrection

!>
!! @brief  Maneja excepciones de la clase
!<
subroutine EpsteinNesbet_exception( typeMessage, description, debugDescription)
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

end subroutine EpsteinNesbet_exception


end module ENFunctions_
