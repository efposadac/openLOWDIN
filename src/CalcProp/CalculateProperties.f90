
!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module CalculateProperties_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use Units_
  use Exception_
  use WaveFunction_
  use ContractedGaussian_
	implicit none

	!>
	!!
	!!  Este modulo define una seudoclase para calculo de propiedades derivadas de
	!! la funcion de onda como cargas, dipolos, polarizabilidades, etc.
	!!
	!! @author Sergio A. Gonzalez Monico
	!!
	!! <b> Fecha de creacion : </b> 2007-09-18
	!!
	!! <b> Historial de modificaciones: </b>
	!!
	!!   - <tt> 2007-09-18 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
	!!        -# Creacion de modulo y metodos basicos.
	!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
	!!   - <tt> 2011-11-23 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
	!!        -# Adds numerical integration properties, ADPT calculations and brings population analyses 
	!!   - <tt> 2014-01-23 </tt>: Matheus Rodriguez ( matrodriguezalv@unal.edu.co )
        !!        -# Reescribe y adapta el modulo de Calculate properties en Lowdin2
	!<

	type, public :: CalculateProperties
        character(30) :: name
       ! type(Matrix) :: contributionsOfdipoleMoment
       ! type(Matrix) :: expectedPositions
       ! type(Vector) :: expectedR2
       ! type(Matrix) :: polarizabilityTensor
       ! type(Matrix) :: hyperPolarizabilityTensor(3)
       ! type(Matrix) :: interparticleDistances
       ! type(Matrix) :: interparticleDistancesErrors
        type(Matrix) :: interparticleOverlap
       ! type(Vector) :: volume
       ! type(Vector) :: cumulativeDensity
       ! type(Cube), allocatable :: densityCube(:)
       ! type(Cube), allocatable :: orbitalCube(:)
        type(Matrix) :: negativeFukui
        type(Matrix) :: positiveFukui
        type(Matrix) :: overlapMatrix   !!! JORGE
        type(Matrix) :: densityMatrix
	end type

	type(CalculateProperties) :: CalculateProperties_instance

       integer, parameter, public :: MULLIKEN  =  1
       integer, parameter, public :: LOWDIN    =  2



	!private :: &
		!CalculateProperties_getDipoleOfPuntualCharges
		! CalculateProperties_getDipoleOfQuantumSpecie

	public :: &
!		CalculateProperties_constructor, &
!		CalculateProperties_destructor, &
	!	CalculateProperties_dipole, &
         !       CalculateProperties_expectedPosition, &
         !       CalculateProperties_expectedR2, &
	!	CalculateProperties_polarizability, &
	!	CalculateProperties_showContributionsToElectrostaticMoment, &
	!	CalculateProperties_showExpectedPositions, &
!		CalculateProperties_showExpectedR2, &
	!	CalculateProperties_showPolarizabilityTensor, &
        !        CalculateProperties_interparticleDistance,  &
!                CalculateProperties_interparticleOverlap, &
         !       CalculateProperties_distanceToPoint, &
!                CalculateProperties_buildDensityCubesLimits, &
 !               CalculateProperties_buildDensityCubes, &
         !       CalculateProperties_volumes, &
		CalculateProperties_showPopulationAnalyses, &
		CalculateProperties_getPopulation
!		CalculateProperties_getPartialCharges, &
!		CalculateProperties_showIonizationPotentials, &
!		CalculateProperties_showCharges
 !               CalculateProperties_showVolumes, &
!                CalculateProperties_getFukuiAt
		
contains



	subroutine CalculateProperties_showPopulationAnalyses()
		implicit none
               ! type (CalculateProperties) :: this ! por medio de este this accedo a todo lo que este en la estructura o type
                                                   ! calculate properties 

		real(8) :: total
		character(10) :: specieName
		integer :: i,specieID
		logical :: showPopulations 


		specieID=1
		!! Recorre las especies buscando electrones
                                        
		search_specie: do i = 1, MolecularSystem_getNumberOfQuantumSpecies()
			specieName=""
			specieName = trim(MolecularSystem_getNameOfSpecie(i))

			if( scan(trim(specieName),"E")==1 ) then
				if( scan(trim(specieName),"-")>1 ) then
					showPopulations=.true.
					specieID=i
					exit search_specie
				end if
			else
				showPopulations=.false.
			end if

		end do search_specie

		if( showPopulations ) then
			!!Obtiene Poblaciones de Mulliken
			print *,""
			print *," POPULATION ANALYSES: "
			print *,"===================="
			print *,""
			print *, " Mulliken Population: "
			print *,"---------------------"
			print *,""
			call Vector_show( CalculateProperties_getPopulation(MULLIKEN, total,trim(specieName)),&
				flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )

			write (6,"(T25,A10)") "__________"
			write (6,"(T10,A15,F10.6)") "Total = ", total
			print *,""
			print *,"...end of Mulliken Population"
			print *,""
			print *, " Lowdin Population:"
			print *,"---------------------"
			print *,""
			call Vector_show( CalculateProperties_getPopulation( LOWDIN, total,trim(specieName)),&
				flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )
			write (6,"(T25,A10)") "__________"
			write (6,"(T10,A15,F10.6)") "Total = ", total
			print *,""
			print *,"...end of Lowdin Population"
			print *,""
			print *,"END POPULATION ANALYSES "
			print *,""
		end if

	end subroutine CalculateProperties_showPopulationAnalyses


 !<
 !! @brief Retorna la poblacion de Mulliken o Lowdin del sistema molecular
 !>
 function CalculateProperties_getPopulation( typeOfPopulation, totalSum, nameOfSpecie, fukuiType )  result( output )
   implicit none
  ! type (CalculateProperties) :: this 
   integer :: typeOfPopulation
   real(8), optional, intent(out) :: totalSum
   character(*),optional  :: nameOfSpecie
   character(*), optional :: fukuiType
   type(Vector) :: output
   
   type(Matrix) :: densityMatrix
   type(Matrix) :: overlapMatrix
   type(Matrix) :: auxMatrix
   type(Matrix) :: auxMatrixB
   character(10) :: auxNameOfSpecie
   integer :: numberOfcontractions
   integer :: speciesID
   integer :: i
   character(50) :: wfnFile
   integer :: wfnUnit
   character(50) :: arguments(20)
   character(50) ::  integralsFile
   integer ::  integralsUnit


  integralsFile = "lowdin.opints"
  integralsUnit = 30


  wfnFile = "lowdin.wfn"
  wfnUnit = 20


  !! Open file for wavefunction                                                                                                                                             
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted") 


  auxNameOfSpecie="E-" 
  if (present( nameOfSpecie ) )	then
     auxNameOfSpecie = trim(nameOfSpecie)
  end if

  !		if ( MolecularSystem_isSet() ) then
  if ( .not. present( fukuiType) .or. (present(fukuiType) .and. trim(auxNameofSpecie) .eq. "E-") ) then
     speciesID =MolecularSystem_getSpecieID (  nameOfSpecie = trim(auxNameOfSpecie) )
     numberOfcontractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )
     call Matrix_constructor( auxMatrix, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )
     call Vector_constructor( output, numberOfcontractions   )


     arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
     arguments(1) = "DENSITY"
     !  WaveFunction_instance(speciesID)%densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

!     print *, "CalculateProperties_getPopulation 0"

     densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
!!! Solo para esta subrutina

     arguments(1) = "OVERLAP"
 

     !!! Abrir el archivo lowdin.opints

     !  WaveFunction_instance(speciesID)%overlapMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))


!     print *, "CalculateProperties_getPopulation 1"

     overlapMatrix = Matrix_getFromFile(unit=integralsUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
  !!! General


!     print *, "CalculateProperties_getPopulation 2"

     select case( typeOfPopulation )

     case( MULLIKEN )
        ! Reto por que cambia bastante en lowdin 2
        ! Leer el modulo wave function de HF en lowdin 2
        ! Leer todo el directorio HF de lowdin 2    
    !!    overlapMatrix = WaveFunction_instance( speciesID )%overlapMatrix
        !  if (trim(fukuiType)=="positive" .and. trim(auxNameOfSpecie) == "E-") then
        !    densityMatrix = this%positiveFukui
        !  else if  (trim(fukuiType)=="negative" .and. trim(auxNameOfSpecie) == "E-") then
        !     densityMatrix = this%negativeFukui
        ! else              
    !!    densityMatrix =  WaveFunction_instance( speciesID )%densityMatrix      
        ! end if
        auxMatrix%values = matmul(densityMatrix%values, overlapMatrix%values )
        
     case (LOWDIN)
        
     !!   overlapMatrix = WaveFunction_instance( speciesID )%overlapMatrix
     !!   densityMatrix =  WaveFunction_instance( speciesID )%densityMatrix          
        auxMatrix%values = matmul(densityMatrix%values, overlapMatrix%values )
        
        auxMatrix = Matrix_pow( overlapMatrix, 0.5_8 )
        auxMatrixB = auxMatrix
        auxMatrix%values = matmul( matmul( auxMatrixB%values , densityMatrix%values), auxMatrixB%values )
        
        call Matrix_destructor(auxMatrixB)
        
     case default
        
     end select

     
     
     do i=1, numberOfcontractions
        output%values(i) = auxMatrix%values(i,i)
     end do
                     
     !print*,"auxMatrix%values", auxMatrix%values      
     if ( present( totalSum ) ) totalSum = sum(output%values)

     call Matrix_destructor(overlapMatrix)
     call Matrix_destructor(auxMatrix)
     call Matrix_destructor(auxMatrixB)
     call Matrix_destructor(densityMatrix)
     
  end if
                       
!else

  ! call CalculateProperties_exception(ERROR, "You should set the molecular system before use this function", &
  !      "Class object CalculateProperties in the getPopulation function" )
   
   close(wfnUnit)
!end if

	end function CalculateProperties_getPopulation



	subroutine CalculateProperties_exception( typeMessage, description, debugDescription)
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
	
	end subroutine CalculateProperties_exception

end module CalculateProperties_
