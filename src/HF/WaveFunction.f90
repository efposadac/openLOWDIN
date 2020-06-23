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
!! @brief Modulo para definicion de funciones de onda RHF
!!  Este modulo define funciones de onda para diferentes especies cuanticas dentro del formalismo OMNE
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-30
!! <b> Historial de modificaciones: </b>
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y metodos
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
module WaveFunction_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use MolecularSystem_
  use CosmoCore_
  implicit none


     !!**************************************************************
     !! Matrices asociadas al calculo de funciones de onda RHF
     !!
     type(vector) :: energyofmolecularorbital
     !! Cosmo Things
     type(Matrix) :: cosmo1
     type(Matrix) :: cosmo2
     type(Matrix) :: cosmo4
     type(Matrix) :: cosmoCoupling
     !!**************************************************************



  end type WaveFunction

  type(WaveFunction), public, allocatable, target :: WaveFunction_instance(:)

  public


contains



  !   function WaveFunction_getValueForOrbitalAt( nameOfSpecie, orbitalNum, coordinate ) result(output)
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     integer :: orbitalNum
  !     real(8) :: coordinate(3)
  !     real(8) :: output

  !     integer :: specieID
  !     character(30) :: nameOfSpecieSelected
  !     integer :: numberOfContractions
  !     integer :: totalNumberOfContractions
  !     integer :: particleID
  !     integer :: contractionID
  !     integer :: i, j
  !     real(8), allocatable :: auxVal(:)


  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )

  !     output=0.0_8
  !     do i=1,numberOfContractions

  !        particleID = MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%particleID
  !        contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%contractionIDInParticle

  !        totalNumberOfContractions = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital

  !        if( allocated(auxVal)) deallocate(auxVal)
  !        allocate(auxVal(totalNumberOfContractions))

  !        auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( specieID,  numberOfContraction=i ), coordinate )

  !        do j = 1, totalNumberOfContractions

  !           output = output + auxVal(j) * WaveFunction_instance( specieID )%waveFunctionCoefficients%values(j,orbitalNum)

  !        end do

  !     end do


  !   end function WaveFunction_getValueForOrbitalAt
  !   !
  !   !
  !   subroutine WaveFunction_draw2DOrbital( nameOfSpecie, orbitalNum, flags )
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     integer :: orbitalNum
  !     integer :: flags


  !     character(30) :: nameOfSpecieSelected
  !     character(50) :: fileName
  !     character(50) :: xRange
  !     integer :: specieID
  !     integer :: numberOfContractions
  !     integer :: j
  !     integer :: i
  !     integer :: numOfGraphs
  !     integer :: auxInitOrbitalNum
  !     integer :: auxLastOrbitalNum

  !     ! 	nameOfSpecieSelected = "e-"
  !     ! 	if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     ! 	specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     ! 	numberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
  !     ! 	fileName=trim(CONTROL_instance%INPUT_FILE)//'orbital.'//trim(String_convertIntegerToString(orbitalNum))//'.'//trim(nameOfSpecieSelected)

  !     ! 	select case( flags )

  !     ! 		case(ORBITAL_ALONE)

  !     ! 			open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
  !     ! 			do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
  !     ! 				CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
  !     ! 			write (5,"(F20.10,F20.10)") j*CONTROL_instance%STEP_OF_GRAPHS, &
  !     ! 			WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected,&
  !     ! 			orbitalNum, [0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) !&
  !     ! !			!!! + (WaveFunction_instance( specieID )%molecularOrbitalsEnergy%values(orbitalNum) * CM_NEG1)
  !     ! 			end do
  !     ! 			close(5)
  !     ! 			call InputOutput_make2DGraph(trim(fileName),&
  !     ! 				"Nuclear Wave Function",&
  !     ! 				"r / Bohr",&
  !     ! 				"U / a.u.", &
  !     ! 				y_format="%.2e")

  !     ! 		case(ORBITAL_WITH_POTENTIAL)

  !     ! 			auxInitOrbitalNum=orbitalNum
  !     ! 			auxLastOrbitalNum=orbitalNum
  !     ! 			numOfGraphs=2
  !     ! 			if(orbitalNum==0) then
  !     ! 				auxInitOrbitalNum=1
  !     ! 				auxLastOrbitalNum=numberOfContractions
  !     ! 				numOfGraphs=numberOfContractions+1
  !     ! 			end if
  !     ! 		open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
  !     ! 		do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
  !     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
  !     ! 			write (5,"(2ES20.10$)") &
  !     ! 			j*CONTROL_instance%STEP_OF_GRAPHS, &
  !     ! 			ExternalPotential_getPotential(ExternalPotential_Manager_instance%externalsPots(1),&
  !     ! 			[j*CONTROL_instance%STEP_OF_GRAPHS,0.0_8,0.0_8])*CM_NEG1
  !     ! 			do i=auxInitOrbitalNum,auxLastOrbitalNum
  !     ! 				write (5,"(ES20.10$)") &
  !     ! 				CONTROL_instance%WAVE_FUNCTION_SCALE&
  !     ! 				*WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected, i,&
  !     ! 					[0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) &
  !     ! 				+ (WaveFunction_instance( specieID )%molecularOrbitalsEnergy%values(i) * CM_NEG1)
  !     ! 			end do
  !     ! 			write (5,"(A)") ""
  !     ! 		end do
  !     ! 		close(5)

  !     ! 		xRange=trim(adjustl(String_convertRealToString(real(&
  !     ! 			-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
  !     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))//':'//trim(adjustl(String_convertRealToString(real(&
  !     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
  !     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))

  !     ! 			call InputOutput_make2DGraph(trim(fileName),&
  !     ! 			"Nuclear Wave Function in potential ",&
  !     ! 			"r / Bohr",&
  !     ! 			"U / cm-1",&
  !     ! 			y_format="%.2e",numOfGraphs=numOfGraphs,x_range=trim(xRange))

  !     ! 	end select

  !   end subroutine WaveFunction_draw2DOrbital

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine WaveFunction_exception( typeMessage, description, debugDescription)
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
  end subroutine WaveFunction_exception

end module WaveFunction_
