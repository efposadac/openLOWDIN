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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!   - <tt> 2014-08-26 </tt>: Jorge Charry (jacharrym@unal.edu.co)
!!        -# Adapts the module to read the transformed integrals from the
!!         integralsTransformation program.
!!        -# Enable the UMP2 method.
!!        -# Consider the particle fraction in the intraspecies term.
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module MPFunctions_
  use CONTROL_
  use MolecularSystem_
  use IndexMap_
  use Exception_
  use Vector_
  use ReadTransformedIntegrals_
  use String_
  implicit none

  !< enum MollerPlesset_correctionFlags {
  integer, parameter :: FIRST_ORDER = 1
  integer, parameter :: SECOND_ORDER = 2
  integer, parameter :: THIRD_ORDER = 3
  !< }

  type :: MollerP
    
        character(50) :: name
        real(8) :: energyHF
        integer :: orderOfCorrection
        integer :: numberOfSpecies
        integer :: frozenCoreBoundary
        real(8) :: totalEnergy
        real(8) :: totalCorrection
        real(8) :: secondOrderCorrection
        real(8) :: thirdOrderCorrection
        !! Vectores para almacenamiento de correcciones a la energia.para cada especie
        type(Vector) :: energyCorrectionOfSecondOrder
        type(Vector) :: energyOfCouplingCorrectionOfSecondOrder

        logical :: isInstanced

  end type MollerP

       type(MollerP), target :: MollerPlesset_instance
!       character(50) :: job
       private :: &
            MollerPlesset_secondOrderCorrection, &
            MollerPlesset_thirdOrderCorrection

       public :: &
            MollerPlesset_constructor, &
            MollerPlesset_destructor, &
            MollerPlesset_show, &
            MollerPlesset_run, &
            MollerPlesset_getTotalEnergy, &
            MollerPlesset_getEnergyCorrection, &
            MollerPlesset_getSpecieCorrection


contains
  !**
  ! Define el constructor para la clase
  !
  !**
  subroutine MollerPlesset_constructor( orderOfCorrection )
    implicit none
    integer, intent(in) :: orderOfCorrection

    integer :: i
    type(Exception) :: ex

!   if( .not.MollerPlesset_instance%isInstanced ) then

      MollerPlesset_instance%orderOfCorrection = orderOfCorrection
      MollerPlesset_instance%numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


      if ( MollerPlesset_instance%orderOfCorrection >= 2 ) then
        call Vector_constructor( MollerPlesset_instance%energyCorrectionOfSecondOrder, MollerPlesset_instance%numberOfSpecies)

        i = MollerPlesset_instance%numberOfSpecies * ( MollerPlesset_instance%numberOfSpecies-1 ) / 2



        call Vector_constructor( MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder, i)
      end if

      if ( MollerPlesset_instance%orderOfCorrection >= 3 ) then

        call Exception_constructor( ex , ERROR )
        call Exception_setDebugDescription( ex, "Class object MollerPlesset in the constructor() function" )
        call Exception_setDescription( ex, "This order correction hasn't been implemented" )
        call Exception_show( ex )

      end if

      MollerPlesset_instance%isInstanced =.true.

!   end if

  end subroutine MollerPlesset_constructor

  !**
  ! Define el destructor para clase
  !
  !**
  subroutine MollerPlesset_destructor()
    implicit none

               MollerPlesset_instance%isInstanced =.false.

  end subroutine MollerPlesset_destructor


  !**
  ! @brief Muestra informacion asociada con la correccion MPn
  !**
  subroutine MollerPlesset_show()
    implicit none

    integer :: i
    integer :: j
    integer :: k

               if ( MollerPlesset_instance%isInstanced )  then

                print *,""
          print *," POST HARTREE-FOCK CALCULATION"
          print *," MANY-BODY PERTURBATION THEORY:"
          print *,"=============================="
          print *,""
          write (6,"(T10,A25)") "MOLLER-PLESSET FORMALISM "
          write (6,"(T10,A23, I5)") "ORDER OF CORRECTION = ",MollerPlesset_instance%orderOfCorrection


        print *,""
        write (6,'(T10,A15,ES20.12)') "E(0) + E(1) = ", MollerPlesset_instance%energyHF
        write (6,'(T10,A15,ES20.12)') "E(2) = ", MollerPlesset_instance%secondOrderCorrection
        write (6,'(T25,A20)') "________________________"
        write (6,'(T10,A15,ES25.17)') "E(MP2)= ", MollerPlesset_instance%totalEnergy
        print *,""
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        write ( 6,'(T10,A15,A20)') " E(n){ Specie } ","   E(n) / Hartree "
        write ( 6,'(T10,A35)') "-----------------------------------------------"
        print *,""

        do i=1, MollerPlesset_instance%numberOfSpecies
       write (*,'(T10,A5,A8,A8,ES16.8)') "E(2){", trim(  MolecularSystem_getNameOfSpecie( i ) ),"   } = ", &
               MollerPlesset_instance%energyCorrectionOfSecondOrder%values(i)
        end do

        print *,""
        k=0
        do i=1, MollerPlesset_instance%numberOfSpecies
       do j=i+1,MollerPlesset_instance%numberOfSpecies
           k=k+1
           write (*,'(T10,A5,A16,A4,ES16.8)') "E(2){", &
               trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
               "} = ", MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder%values(k)
       end do
        end do

               end if

  end subroutine MollerPlesset_show

  !**
  ! @brief realiza la correccion de energia MPn orden dado
  !**
  subroutine MollerPlesset_run()
    implicit none

    type(Exception) :: ex
    integer :: i

  character(50) :: wfnFile
  integer :: wfnUnit
  wfnFile = "lowdin.wfn"
  wfnUnit = 20

  
  if ( MollerPlesset_instance%isInstanced ) then

     do i=2, MollerPlesset_instance%orderOfCorrection

        select case( i )

        case(  SECOND_ORDER )

           call MollerPlesset_secondOrderCorrection()

           MollerPlesset_instance%totalCorrection= MollerPlesset_instance%secondOrderCorrection
        case(  THIRD_ORDER )

           call MollerPlesset_thirdOrderCorrection()

           MollerPlesset_instance%totalCorrection= MollerPlesset_instance%totalCorrection +&
                MollerPlesset_instance%thirdOrderCorrection

        case default

           call Exception_constructor( ex , ERROR )
           call Exception_setDebugDescription( ex, "Class object MollerPlesset in run() function" )
           call Exception_setDescription( ex, "This order correction hasn't been implemented" )
           call Exception_show( ex )

        end select

     end do
  !! Open file for wavefunction
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  !! Load results...

   call Vector_getFromFile(unit=wfnUnit, binary=.true., value=MollerPlesset_instance%energyHF, arguments=["TOTALENERGY"])

      MollerPlesset_instance%totalEnergy = MollerPlesset_instance%energyHF + MollerPlesset_instance%totalCorrection

close(wfnUnit)

    else

      call Exception_constructor( ex , ERROR )
      call Exception_setDebugDescription( ex, "Class object MollerPlesset in run() function" )
      call Exception_setDescription( ex, "You should to instance MollerPlesset module before use this function" )
      call Exception_show( ex )

end if


  end subroutine MollerPlesset_run


  !**
  ! @ Retorna el la correccion a la energia
  !**
  function MollerPlesset_getEnergyCorrection() result( output)
    implicit none
    real(8) :: output

    output = MollerPlesset_instance%totalCorrection

  end function MollerPlesset_getEnergyCorrection

  !**
  ! @ Retorna el la correccion a la energia
  !**
  function MollerPlesset_getSpecieCorrection( specieName) result( output)
    implicit none
    character(*) :: specieName
    real(8) :: output

    integer :: i

    do i=1, MollerPlesset_instance%numberOfSpecies

      if ( trim(specieName) == trim(  MolecularSystem_getNameOfSpecie( i ) ) ) then

        output = MollerPlesset_instance%energyCorrectionOfSecondOrder%values(i)
        return

      end if

    end do

  end function MollerPlesset_getSpecieCorrection



  !**
  ! @ Retorna la energia final con correccion Moller-Plesset de orrden dado
  !**
  function MollerPlesset_getTotalEnergy() result(output)
    implicit none
    real(8) :: output

    output = MollerPlesset_instance%totalEnergy

  end function MollerPlesset_getTotalEnergy


 !<
 !! @brief Realiza la correccion de segundo orden a la enegia
 !>
 subroutine MollerPlesset_secondOrderCorrection()
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
   integer :: is
   integer :: specieID
   integer :: otherSpecieID
   integer :: electronsID
   integer :: ocupationNumber
   integer :: ocupationNumberOfOtherSpecie
   integer :: numberOfContractions
   integer :: numberOfContractionsOfOtherSpecie
   integer(8) :: auxIndex
   character(10) :: nameOfSpecie
   character(10) :: nameOfOtherSpecie
   type(Vector) :: eigenValues
   type(Vector) :: eigenValuesOfOtherSpecie
   type(Matrix) :: auxMatrix
!   type(TransformIntegrals) :: repulsionTransformer
   real(8) :: lambda
   real(8) :: lambdaOfOtherSpecie
   real(8) :: independentEnergyCorrection, E1, E2
   real(8) :: couplingEnergyCorrection
   real(8) :: auxVal_A
   real(8) :: auxVal_B
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

   do is=1, MollerPlesset_instance%numberOfSpecies

      
      MollerPlesset_instance%frozenCoreBoundary = 1

      if ( is == electronsID )  MollerPlesset_instance%frozenCoreBoundary = &
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

         if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

           call ReadTransformedIntegrals_readOneSpecies( specieID, auxMatrix)

!           call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!                eigenVec, auxMatrix, specieID, trim(nameOfSpecie) )

                  !!**************************************************************************
         !         !!  Calcula la correccion de segundo orden a la energia
                  !!****
           do a=MollerPlesset_instance%frozenCoreBoundary, ocupationNumber
             do b=MollerPlesset_instance%frozenCoreBoundary,ocupationNumber
               do r=ocupationNumber+1, numberOfContractions
                 do s=r, numberOfContractions

                   auxIndex = IndexMap_tensorR4ToVectorB(int(a,8),int(r,8),int(b,8),int(s,8), &
                                                        int(numberOfContractions,8) )

                   auxVal_A= auxMatrix%values(auxIndex, 1)
                   if (  dabs( auxVal_A)  > 1.0E-10_8 ) then

                     !if ( b >= a) print *, "A", a, b, r, s, auxVal_A
                     !if ( b < a) print *, "B", a, b, r, s, auxVal_A

                     if ( s>r ) then

                       if (a==b) then

                         if( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                           independentEnergyCorrection = independentEnergyCorrection + 2.0_8 *  auxVal_A**2.0  &
                                                        * ( lambda  -  1.0_8 ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                        - eigenValues%values(r) - eigenValues%values(s) )

                         end if

                       else

                         auxIndex = IndexMap_tensorR4ToVector(r, b, s, a, numberOfContractions )
                         auxVal_B= auxMatrix%values(auxIndex, 1)

                         independentEnergyCorrection = independentEnergyCorrection + 2.0_8 *  auxVal_A  &
                                                     * ( lambda * auxVal_A  - auxVal_B ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                     - eigenValues%values(r) - eigenValues%values(s) )

                       end if

                     else if ( a==b .and. r==s ) then

                       if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                         independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0_8  &
                                                      * ( lambda - 1.0_8 ) / ( 2.0_8*( eigenValues%values(a)-eigenValues%values(r)))

                       end if


                     else

                       if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                         independentEnergyCorrection = independentEnergyCorrection + auxVal_A**2.0  &
                                                      * ( lambda  - 1.0_8 ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                      - eigenValues%values(r) - eigenValues%values(s) )

                       end if

                     end if
                   end if

                 end do
               end do
             end do
           end do

         else !! DIRECT

           integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE
           prefixOfFile =""//trim(nameOfSpecie)
           unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
           nVirtualOrbitals = numberOfContractions - ocupationNumber 

           allocate (intArray(numberOfContractions, numberOfContractions, numberOfContractions))
           intArray = 0
           i1 = 1 !! replace to frozen core

           !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
           open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
             status='old',access='sequential', form='unformatted' )

           readIntegralsC : do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) ii, aa, jj, bb, shellIntegrals

             do p = 1, integralStackSize

               if ( ii(p) /= i1 ) then
                 i = i1
                 i1 = ii(p)

                 !print *, i,j
                 do j = 1,ocupationNumber
                   do a = ocupationNumber+1, numberOfContractions
                     do b = a, numberOfContractions
                     
                       !!auxIndex = IndexMap_tensorR4ToVectorB( int(pp(i),8), int(qq(i),8), int(rr(i),8), int(ss(i),8), int(numberOfContractions,8 ))
                       !!matrixContainer%values( auxIndex, 1 ) = auxIntegrals(i)
                       auxVal_A = intArray(j,a,b)
                       !if ( abs(auxVal_A ) > 1E-10) print *, "A", i, j, a, b, auxVal_A
                       if ( abs(auxVal_A ) > 1E-10) then
                     
                         if ( b > a ) then
                      
                           if ( i == j ) then
                      
                             if( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                      
                               independentEnergyCorrection = independentEnergyCorrection + 2.0_8 *  auxVal_A**2.0  &
                                                         * ( lambda  -  1.0_8 ) / ( eigenValues%values(i) + eigenValues%values(j) &
                                                         - eigenValues%values(a) - eigenValues%values(b) )
                      
                             end if
                      
                           else
                      
                             !auxIndex = IndexMap_tensorR4ToVector(r, b, s, a, numberOfContractions )
                             !auxVal_B= auxMatrix%values(auxIndex, 1)
                             auxVal_B = intArray(j,b,a)
                      
                             independentEnergyCorrection = independentEnergyCorrection + 2.0_8 * auxVal_A  &
                                                      * ( lambda * auxVal_A  - auxVal_B ) / ( eigenValues%values(i) + eigenValues%values(j) &
                                                      - eigenValues%values(a) - eigenValues%values(b) )
                           end if
                      
                         else if ( i == j .and. a == b ) then
                      
                           if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                      
                           independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0_8  &
                                                       * ( lambda - 1.0_8 ) / ( 2.0_8*( eigenValues%values(i)-eigenValues%values(a)))

                           end if
                      
                         else
                      
                           if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                             independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0  &
                                                       * ( lambda  - 1.0_8 ) / ( eigenValues%values(i) + eigenValues%values(j) &
                                                       - eigenValues%values(a) - eigenValues%values(b) )

                           end if
                      
                         end if
                       end if !! >1E-10
                     end do !! b

                     do b = 1, a

                       if ( j > i ) then
                         auxVal_A = intArray(j,a,b)
                         !if ( abs(auxVal_A ) > 1E-10) print *, "B", i, j, a, b, auxVal_A
                     
                         if ( a > b ) then
                     
                           !auxIndex = IndexMap_tensorR4ToVector(r, b, s, a, numberOfContractions )
                           !auxVal_B= auxMatrix%values(auxIndex, 1)
                           auxVal_B = intArray(j,b,a)
                     
                           independentEnergyCorrection = independentEnergyCorrection + 2.0_8 * auxVal_A  &
                                                       * ( lambda * auxVal_A  - auxVal_B ) / ( eigenValues%values(i) + eigenValues%values(j) &
                                                       - eigenValues%values(a) - eigenValues%values(b) )
                      
                         else
                           if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                             independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0  &
                                                      * ( lambda  - 1.0_8 ) / ( eigenValues%values(i) + eigenValues%values(j) &
                                                      - eigenValues%values(a) - eigenValues%values(b) )

                            end if
                     
                         end if
                       end if
                     end do
                   end do
                 end do

                 intArray = 0
               end if  
                     
               if ( ii(p) == -1_8 ) exit readIntegralsC
                 !print *, "===", ii(p), jj(p),aa(p),bb(p), shellIntegrals(p)

               if ( ii(p) == i1 ) then
                 intArray(jj(p),aa(p),bb(p))  = shellIntegrals(p)
                  !intArray(jj(p),bb(p),aa(p))  = shellIntegrals(p)
                     
               end if
                     
             end do  
                     
           end do readIntegralsC

           close(unidOfOutputForIntegrals)
           deallocate(intArray)
                     
                     
         end if !!direct
                     
       end if        
                     
    MollerPlesset_instance%energyCorrectionOfSecondOrder%values(is) = independentEnergyCorrection &
      * ( ( MolecularSystem_getCharge( specieID ) )**4.0_8 )
                     
    if ( nameOfSpecie == "E-ALPHA" .or. nameOfSpecie == "E-BETA" ) then
                     
      MollerPlesset_instance%energyCorrectionOfSecondOrder%values(is) = &
        MollerPlesset_instance%energyCorrectionOfSecondOrder%values(is) / ( 2.0 )
                     
    else             
                     
      MollerPlesset_instance%energyCorrectionOfSecondOrder%values(is) = &
        MollerPlesset_instance%energyCorrectionOfSecondOrder%values(is) / &
        ( MolecularSystem_getParticlesFraction ( specieID ) * 2.0 )
    end if           
                     
                     
    call Matrix_destructor(auxMatrix)
    !!
    !!**************************************************************************

  end do


  !! Suma las correcciones de energia para especies independientes
  MollerPlesset_instance%secondOrderCorrection = sum( MollerPlesset_instance%energyCorrectionOfSecondOrder%values )

  !!
  !!*******************************************************************************************


  !!*******************************************************************************************
  !! Calculo de correcion se segundo orden para interaccion entre particulas de especie diferente
  !!
  if ( MollerPlesset_instance%numberOfSpecies > 1 ) then

   couplingEnergyCorrection = 0.0_8
   m=0

   do i = 1 , MollerPlesset_instance%numberOfSpecies

     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(i)
     arguments(2) = trim(MolecularSystem_getNameOfSpecie(i))

     arguments(1) = "COEFFICIENTS"
     eigenVec = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "ORBITALS"
     call Vector_getFromFile( elementsNum = numberOfContractions, &
          unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = eigenValues )     

     nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
     specieID =MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecie) )
     ocupationNumber = MolecularSystem_getOcupationNumber( i )
     numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
     lambda = MolecularSystem_getEta( i )

     do j = i + 1 , MollerPlesset_instance%numberOfSpecies
        m = m + 1


        numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )

        arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))

        arguments(1) = "COEFFICIENTS"
        eigenVecOtherSpecie = &
             Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractionsOfOtherSpecie,4), &
             columns= int(numberOfContractionsOfOtherSpecie,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = numberOfContractionsofOtherSpecie, &
             unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
             output = eigenValuesOfOtherSpecie )     


        nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
        otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
        ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )

        lambdaOfOtherSpecie = MolecularSystem_instance%species(j)%lambda
        couplingEnergyCorrection = 0.0_8

!        if ( .not.CONTROL_instance%OPTIMIZE ) then
!           write (6,"(T10,A)") "INTER-SPECIES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)//"/"//trim(nameOfOtherSpecie)
!           print *,""
!        end if

         !! Read transformed integrals from file
        call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, auxMatrix)

!        call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!             eigenVec, eigenVecOtherSpecie, &
!             auxMatrix, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )

        auxMatrix%values = auxMatrix%values * MolecularSystem_getCharge( specieID ) &
             * MolecularSystem_getCharge( otherSpecieID )

        do a=1, ocupationNumber
           do p=1,ocupationNumberOfOtherSpecie
              do r=ocupationNumber+1, numberOfContractions
                 do t=ocupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie

                    auxIndex = IndexMap_tensorR4ToVector(a,r,p,t, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )
                    couplingEnergyCorrection = couplingEnergyCorrection +  &
                         ( ( auxMatrix%values(auxIndex,1) )**2.0_8 ) &
                         / (eigenValues%values(a) + eigenValuesOfOtherSpecie%values(p) &
                         -eigenValues%values(r)-eigenValuesOfOtherSpecie%values(t) )

                 end do
              end do
           end do
        end do

        MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder%values(m)= &
             ( lambda * lambdaOfOtherSpecie * couplingEnergyCorrection )  

     end do
  end do

  call Matrix_destructor(auxMatrix)
  !! Adiciona la correccion del termino de acoplamiento
  MollerPlesset_instance%secondOrderCorrection = MollerPlesset_instance%secondOrderCorrection + &
       sum( MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder%values )

end if
!!
!!*******************************************************************************************

close(wfnUnit)

end subroutine MollerPlesset_secondOrderCorrection

!<
!! @brief Realiza la correccion de tercer orden a la enegia
!>
subroutine MollerPlesset_thirdOrderCorrection()
  implicit none

  call MollerPlesset_exception( ERROR, "Class object MollerPlesset in thirdOrderCorrection() function" , &
       "This order correction hasn't been implemented" )

end subroutine MollerPlesset_thirdOrderCorrection

!>
!! @brief  Maneja excepciones de la clase
!<
subroutine MollerPlesset_exception( typeMessage, description, debugDescription)
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

end subroutine MollerPlesset_exception


end module MPFunctions_
