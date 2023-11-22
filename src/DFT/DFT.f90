!This code is part of LOWDIN Quantum chemistry package                 
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
!! @brief DFT program.
!!        This module allows calculations in the APMO-DFT framework
!! @author F. Moncada
!!
!! <b> Creation date : </b> 2017
!!
!! <b> History: </b>
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program DFT
  use CONTROL_
  use MolecularSystem_
  use DensityFunctionalTheory_
  use GridManager_
  use String_
  use Matrix_
  use Exception_
  use omp_lib
  implicit none

  character(50) :: job
  character(100) :: densFile
  type(Matrix), allocatable :: densityMatrix(:)
  type(Matrix), allocatable :: exchangeCorrelationMatrix(:)
  type(Matrix) :: exchangeCorrelationEnergy
  real(8), allocatable :: numberOfParticles(:)
  character(100) :: excFile
  character(50) ::  labels(2)
  integer :: densUnit, excUnit
  integer :: numberOfContractions
  integer :: numberOfSpecies
  integer :: speciesID, otherSpeciesID
  
  job = ""
  call get_command_argument(1,value=job)
  job = trim(String_getUppercase(job))


  densFile = ""
  call get_command_argument(2,value=densFile)

  ! write(*,"(A,A)") trim(job), trim(densFile)
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  call Functional_createFunctionals( )

  !!!Building grids jobs
  select case ( job ) 
  case ("BUILD_SCF_GRID")
     call DensityFunctionalTheory_buildSCFGrid()
     return
  case ("BUILD_FINAL_GRID" )
     call DensityFunctionalTheory_buildFinalGrid()
     return
  end select

  !!!Computing energy and potential jobs  
  numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
  
  allocate( densityMatrix(numberOfSpecies) , numberOfParticles(numberOfSpecies), &
       exchangeCorrelationMatrix(numberOfSpecies))

  call Matrix_constructor(exchangeCorrelationEnergy,int(numberOfSpecies,8),int(numberOfSpecies,8),0.0_8)

  do speciesID = 1 , numberOfSpecies
     ! Read density matrices    
     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
     
     densUnit = 78
     open(unit = densUnit, file=trim(densFile), status="old", form="unformatted")
     
     labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
     labels(1) = "DENSITY-MATRIX"
     densityMatrix(speciesID) =Matrix_getFromFile(unit=densUnit, rows= int(numberOfContractions,4), &
          columns=int(numberOfContractions,4), binary=.true., arguments=labels)
     
     close(unit=densUnit)
  end do
  
  select case ( job ) 
  case ("SCF_DFT")
     call DensityFunctionalTheory_SCFDFT(densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
  case ("FINAL_DFT")
     !read scf information for comparison
     do speciesID = 1 , numberOfSpecies
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
        excUnit = 79
        excFile = trim(densfile)//".exc"
        open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")

        labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
        labels(1) = "NUMBER-OF-PARTICLES"
        call Vector_getFromFile(unit=excUnit, binary=.true., value=numberOfParticles(speciesID), arguments= labels )
        labels(1) = "EXCHANGE-CORRELATION-MATRIX"
        exchangeCorrelationMatrix(speciesID)=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
             binary=.true., arguments=labels(1:2))

        do otherSpeciesID = speciesID, numberOfSpecies
           labels(2) = trim(MolecularSystem_getNameOfSpecies(speciesID))//trim(MolecularSystem_getNameOfSpecies(otherSpeciesID))
           labels(1) = "EXCHANGE-CORRELATION-ENERGY"
           call Vector_getFromFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy%values(speciesID,otherSpeciesID), arguments= labels )
        end do

        close(unit=excUnit)
     end do

     call DensityFunctionalTheory_finalDFT(densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
  ! case default 
  !    write(*,*) "USAGE: lowdin-DFT.x job "
  !    write(*,*) "Where job can be: "
  !    write(*,*) "  initialize"
  !    write(*,*) "  build_matrices"
  !    STOP "ERROR At DFT program, requested an unknown job type"
  end select

  excUnit = 79
  excFile = trim(densfile)//".exc"
  open(unit = excUnit, file=trim(excFile), status="replace", form="unformatted")
  ! Write results to file
  do speciesID = 1 , numberOfSpecies

     ! print *, trim(MolecularSystem_getNameOfSpecies(speciesID)), numberOfParticles(speciesID)
     ! call Matrix_show(exchangeCorrelationMatrix(speciesID))

     labels(2) = trim(MolecularSystem_getNameOfSpecies(speciesID))
     labels(1) = "NUMBER-OF-PARTICLES"
     call Vector_writeToFile(unit=excUnit, binary=.true., value=numberOfParticles(speciesID), arguments= labels(1:2) )

     labels(1) = "EXCHANGE-CORRELATION-MATRIX"
     call Matrix_writeToFile( exchangeCorrelationMatrix(speciesID), unit=excUnit, binary=.true., arguments = labels(1:2) )

     do otherSpeciesID = speciesID, numberOfSpecies
        labels(2) = trim(MolecularSystem_getNameOfSpecies(speciesID))//trim(MolecularSystem_getNameOfSpecies(otherSpeciesID))
        labels(1) = "EXCHANGE-CORRELATION-ENERGY"
        call Vector_writeToFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy%values(speciesID,otherSpeciesID), arguments= labels(1:2) )
     end do

  end do
  close(unit=excUnit)

  
end program DFT


