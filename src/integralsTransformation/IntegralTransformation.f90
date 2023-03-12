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
!! @brief Integrals transformation program.
!!        This module allows to make calculations in the APMO framework
!! @author  J. A. Charry, J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2014-08-26
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
! this is a change
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!   - <tt> 2014-08-26 </tt>: Jorge Charry (jacharrym@unal.edu.co)
!!        -# Write this program to calculate the transformed integrals in Lowdin2
!!           independently from MP2 program.
!! @warning This programs only works linked to lowdincore library,
!!          provided by LOWDIN quantum chemistry package
!!
program IntegralsTransformation
  use CONTROL_
  use MolecularSystem_
  use InputCI_
  use IndexMap_
  use Matrix_
  use Exception_
  use Vector_
  use TransformIntegralsA_
  use TransformIntegralsB_
  use TransformIntegralsC_
  use TransformIntegralsD_
  use TransformIntegralsE_
  use String_
  implicit none

  character(50) :: job
  integer :: i, j, z
  integer :: speciesID, otherSpeciesID
  integer :: numberOfContractions
  integer :: numberOfContractionsOfOtherSpecie
  integer :: occupation, otherOccupation
  character(10) :: nameOfSpecies
  character(10) :: nameOfOtherSpecie
  type(Vector) :: eigenValues
  type(Vector) :: eigenValuesOfOtherSpecie
  type(Matrix) :: auxMatrix
  type(TransformIntegralsA) :: repulsionTransformer
  type(TransformIntegralsB) :: transformInstanceB
  type(TransformIntegralsC) :: transformInstanceC
  type(TransformIntegralsD) :: transformInstanceD
  type(TransformIntegralsE) :: transformInstanceE
  type(Matrix) :: eigenVec
  type(Matrix) :: densityMatrix
  type(Matrix) :: otherdensityMatrix
  type(Matrix) :: eigenVecOtherSpecie 
  character(50) :: wfnFile
  character(50) :: arguments(2)
  character(50) :: partialTransform
  integer :: wfnUnit
  integer :: numberOfQuantumSpecies
  logical :: transformThisSpecies
  logical :: transformTheseSpecies
  real(8) :: timeA, timeB
  
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )
  
  !Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  ! call MolecularSystem_showInformation()  
  ! call MolecularSystem_showParticlesInformation()
  ! call MolecularSystem_showCartesianMatrix()

  wfnFile = "lowdin.wfn"
  wfnUnit = 20

  call InputCI_constructor( )
  call InputCI_load( MolecularSystem_getNumberOfQuantumSpecies() )
  
  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  !!Start time
  timeA = omp_get_wtime()

  !! Checks if all integrals are needed
  if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION == 2 .and. &
       CONTROL_instance%PT_ORDER == 0 .and. &
       CONTROL_instance%EPSTEIN_NESBET_CORRECTION == 0 .and. &
       CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL == "NONE") then
     partialTransform="MP2"
  else if ( CONTROL_instance%PT_ORDER == 2 .and. &
       CONTROL_instance%MOLLER_PLESSET_CORRECTION == 0 .and. &
       CONTROL_instance%EPSTEIN_NESBET_CORRECTION == 0 .and. &
       CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL == "NONE") then
     partialTransform="PT2"
  else if ( CONTROL_instance%PT_ORDER == 2 .and. &
       CONTROL_instance%MOLLER_PLESSET_CORRECTION == 2 .and. &
       CONTROL_instance%EPSTEIN_NESBET_CORRECTION == 0 .and. &
       CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL == "NONE") then
     partialTransform="MP2-PT2"
  else if (CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL .ne. "NONE") then
     ! partialTransform="ALLACTIVE"
     partialTransform="ALL"
  else
     partialTransform="BOUNDS"
  end if
  
  write(*,*) ""
  write(*,*) "BEGIN FOUR-INDEX INTEGRALS TRANFORMATION:"
  write(*,*) "========================================"
  write(*,*) "Selected method: ",CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD 
  write(*,*) ""


  select case( CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD )

    case ( "A" ) 

      call TransformIntegralsA_show
      call TransformIntegralsA_constructor( repulsionTransformer )

    case ( "B" ) 

      call TransformIntegralsB_show
      call TransformIntegralsB_constructor( transformInstanceB )

    case ( "C" ) 

      call TransformIntegralsC_show
      call TransformIntegralsC_constructor( transformInstanceC,partialTransform )

    case ( "D" )

      call TransformIntegralsD_show
      call TransformIntegralsD_constructor( transformInstanceD,partialTransform )

    case ( "E" )

      call TransformIntegralsE_show
      call TransformIntegralsE_constructor( transformInstanceE,partialTransform )

  end select

  !!*******************************************************************************************
  !! BEGIN

  !! get the number of nonzero integrals

  numberOfQuantumSpecies = MolecularSystem_getNumberOfQuantumSpecies()            
  
    do i=1, numberOfQuantumSpecies
  
        nameOfSpecies = trim( MolecularSystem_getNameOfSpecie( i ) )

        !! For PT = 2 there is no need to transform integrals for all species"
        if ( partialTransform == "PT2" .and. CONTROL_instance%IONIZE_SPECIE(1) /= "NONE" ) then
          transformThisSpecies = .False.
          do z = 1, size(CONTROL_instance%IONIZE_SPECIE )
            if ( nameOfSpecies == CONTROL_instance%IONIZE_SPECIE(z) )  then
              transformThisSpecies = .True.
            end if
          end do
        else 
          transformThisSpecies = .True.
        end if

          if ( .not.CONTROL_instance%OPTIMIZE .and. transformThisSpecies) then
              write (*,*) ""
              write (6,"(T2,A)")"Integrals transformation for: "//trim(nameOfSpecies)
              write (*,*) ""
           end if

           !! Reading the coefficients
           open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(i)
           occupation = MolecularSystem_getOcupationNumber( i )
           arguments(2) = MolecularSystem_getNameOfSpecie(i)
  
           arguments(1) = "COEFFICIENTS"
           eigenVec= Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                columns= int(max(numberOfContractions,occupation),4), binary=.true., arguments=arguments(1:2))

           arguments(1) = "DENSITY"
           densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
 
          arguments(1) = "ORBITALS"
           call Vector_getFromFile( elementsNum = numberOfContractions, &
                unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                output = eigenValues )     
           close(wfnUnit)

           speciesID = i

          !! Transforms integrals for one species

          if ( transformThisSpecies .eqv. .True. ) then

            select case( CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD )
          
              case ( "A" ) 
          
                call TransformIntegralsA_atomicToMolecularOfOneSpecie( repulsionTransformer, &
                 eigenVec, auxMatrix, speciesID, trim(nameOfSpecies) )
          
              case ( "B" ) 
          
                call TransformIntegralsB_atomicToMolecularOfOneSpecie(  transformInstanceB, &
                       eigenVec, auxMatrix, speciesID, trim(nameOfSpecies) )

              case ( "C" ) 

                call TransformIntegralsC_atomicToMolecularOfOneSpecie(  transformInstanceC, &
                       densityMatrix, eigenVec, speciesID, trim(nameOfSpecies) ) 

              case ( "D" )

                call TransformIntegralsD_atomicToMolecularOfOneSpecie( transformInstanceD, &
                       eigenVec, speciesID,  trim(nameOfSpecies))

              case ( "E" )

                call TransformIntegralsE_atomicToMolecularOfOneSpecie(  transformInstanceE, &
                       eigenVec, auxMatrix, speciesID, trim(nameOfSpecies) ) 

            end select

          end if

          !!*******************************************************************************************
          !! Two species
          !!
          if ( numberOfQuantumSpecies > 1 ) then
                  do j = i + 1 , numberOfQuantumSpecies
                          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )

                          !! For PT = 2 there is no need to transform integrals for all species"
                          if ( partialTransform == "PT2" .and. CONTROL_instance%IONIZE_SPECIE(1) /= "NONE" ) then
                            transformTheseSpecies = .False.
                            do z = 1, size(CONTROL_instance%IONIZE_SPECIE )
                              if ( nameOfSpecies == CONTROL_instance%IONIZE_SPECIE(z) .or. &
                                   nameOfOtherSpecie == CONTROL_instance%IONIZE_SPECIE(z) )  then
                                transformTheseSpecies = .True.
                              end if
                            end do
                          else 
                            transformTheseSpecies = .True.
                          end if

                  
                          if ( .not.CONTROL_instance%OPTIMIZE .and. transformTheseSpecies ) then
                             write (*,*) ""
                             write (6,"(T2,A)") "Inter-species integrals transformation for: "//trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecie)
                             write (*,*) ""
                          end if

                          !! Reading the coefficients                          
                          open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 
                          numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
                          otherOccupation = MolecularSystem_getOcupationNumber( j )

                          arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))

                          arguments(1) = "COEFFICIENTS"
                          eigenVecOtherSpecie = &
                                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractionsOfOtherSpecie,4), &
                                  columns= int(max(numberOfContractionsOfOtherSpecie,otherOccupation),4), binary=.true., arguments=arguments(1:2))

                          arguments(1) = "ORBITALS"
                          call Vector_getFromFile( elementsNum = numberOfContractionsofOtherSpecie, &
                               unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                               output = eigenValuesOfOtherSpecie )     

                          arguments(1) = "DENSITY"
                          otherdensityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractionsOfOtherSpecie,4), &
                          columns= int(numberOfContractionsOfOtherSpecie,4), binary=.true., arguments=arguments(1:2))
                          close(wfnUnit)
           
                          otherSpeciesID = j
  
                          !! Transforms integrals for two species

                          if ( transformTheseSpecies .eqv. .True. ) then

                            select case( CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD )
                          
                              case ( "A" ) 

                                call TransformIntegralsA_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
                                 eigenVec, eigenVecOtherSpecie, &
                                 auxMatrix, speciesID, nameOfSpecies, otherSpeciesID, nameOfOtherSpecie )
                          
                               case ("B")

                              call TransformIntegralsB_atomicToMolecularOfTwoSpecies(transformInstanceB, &
                                 eigenVec, eigenVecOtherSpecie, &
                                 auxMatrix, speciesID, nameOfSpecies, otherSpeciesID, nameOfOtherSpecie )

                              case ( "C" ) 

                              if ( occupation <= otherOccupation ) then

                                call TransformIntegralsC_atomicToMolecularOfTwoSpecies(transformInstanceC, &
                                 densityMatrix, eigenVec, eigenVecOtherSpecie, &
                                 auxMatrix, speciesID, nameOfSpecies, otherSpeciesID, nameOfOtherSpecie)

                              else 

                                 call TransformIntegralsC_atomicToMolecularOfTwoSpecies(transformInstanceC, &
                                   otherdensityMatrix, eigenVecOtherSpecie, eigenVec, &
                                   auxMatrix, otherSpeciesID, nameOfOtherSpecie,  speciesID, nameOfSpecies)

                              end if

                            case ( "D" )
                              
                              call TransformIntegralsD_atomicToMolecularOfTwoSpecies( transformInstanceD, &
                                 eigenVec, eigenVecOtherSpecie, &
                                 speciesID, nameOfSpecies, &
                                 otherSpeciesID, nameOfOtherSpecie)

                            case ( "E" )
                              
                              call TransformIntegralsE_atomicToMolecularOfTwoSpecies(transformInstanceE, &
                                 eigenVec, eigenVecOtherSpecie, &
                                 auxMatrix, speciesID, nameOfSpecies, otherSpeciesID, nameOfOtherSpecie )

                            end select

                          end if 
                  end do
           end if
     end do

  timeB = omp_get_wtime()
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for integrals transformation : ", timeB - timeA ," (s)"
  write(*, *) ""
  close(30)

close(wfnUnit)

end program IntegralsTransformation

