!!This code is part of LOWDIN Quantum chemistry package                 
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
!!        This module allows to make calculations in the APMO-DFT framework
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
  use GridManager_
  use String_
  use Exception_
  use omp_lib
  implicit none

  character(50) :: job
  type(Matrix), allocatable :: densityMatrix(:), overlapMatrix(:)
  type(Matrix), allocatable :: exchangeCorrelationMatrix(:)
  real(8), allocatable :: exchangeCorrelationEnergy(:), numberOfParticles(:)
  character(50) :: densFile, dftFile, excFile, wfnFile,labels(2)
  integer :: densUnit, dftUnit, excUnit, wfnUnit

  integer :: numberOfSpecies
  integer :: numberOfContractions
  integer :: speciesID, otherSpeciesID, otherElectronID
  character(50) :: nameOfSpecies, nameOfOtherSpecies

  integer :: i, u, dir
  real(8) :: sumCheck, auxEnergy, otherAuxEnergy, otherElectronAuxEnergy
  real(8) :: time1, time2, time3
  
  job = ""
  call get_command_argument(1,value=job)
  job = trim(String_getUppercase(job))


  densFile = ""
  call get_command_argument(2,value=densFile)
  
  ! write(*,"(A)") trim(job)
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

  if (trim(job).eq."INITIALIZE")  then

     !!Start time
     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)
     
     print *, ""
     print *, "--------------------------------------------------------------------------------------"
     print *, "|---------------------Building DFT Integration Grids---------------------------------|"
     print *, "--------------------------------------------------------------------------------------"
     print *, "Euler-Maclaurin radial grids - Lebedev angular grids"
     print *, ""


     call GridManager_buildGrids( "INITIAL" )
     call Functional_createFunctionals( )
     call Functional_show( )
     call GridManager_writeGrids( )
     call GridManager_atomicOrbitals( "WRITE" )
     ! call Stopwatch_stop(lowdin_stopwatch)     
     ! write(*,"(A,F10.3,A4)") "** Building and writing grids and atomic orbitals:", lowdin_stopwatch%enlapsetTime ," (s)"

  else
     call Functional_createFunctionals( )

     allocate( densityMatrix(numberOfSpecies) , numberOfParticles(numberOfSpecies), exchangeCorrelationMatrix(numberOfSpecies), exchangeCorrelationEnergy(numberOfSpecies))
     
     if(trim(job).eq."BUILD_MATRICES") then
        call GridManager_readGrids( )
        call GridManager_atomicOrbitals( "READ" )
        
     else if(trim(job).eq."FINAL_GRID") then
        print *, ""
        print *, "--------------------------------------------------------------------------------------"
        print *, "|---------------------------- DFT Final Integration ---------------------------------|"
        print *, "--------------------------------------------------------------------------------------"
        print *, ""
        call GridManager_buildGrids( "FINAL" )
        call GridManager_atomicOrbitals( "GET" )

        do speciesID = 1 , numberOfSpecies
           excUnit = 79
           excFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".excmatrix"
           open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")

           labels(2) = Grid_instance(speciesID)%nameOfSpecies
           labels(1) = "NUMBER-OF-PARTICLES"
           call Vector_getFromFile(unit=excUnit, binary=.true., value=numberOfParticles(speciesID), arguments= labels )
           labels(1) = "EXCHANGE-CORRELATION-ENERGY"
           call Vector_getFromFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(speciesID), arguments= labels )

           close(unit=excUnit)
           write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the SCF grid: ", numberOfParticles(speciesID)
        end do

        print *, ""
        write (*,"(A50, F15.8)") "Exchange correlation energy with the SCF grid: ", sum(exchangeCorrelationEnergy)
        print *, ""

     else
        write(*,*) "USAGE: lowdin-DFT.x job "
        write(*,*) "Where job can be: "
        write(*,*) "  initialize"
        write(*,*) "  build_matrices"
        STOP "ERROR At DFT program, requested an unknown job type"
     end if

    !!Start time
     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)

     do speciesID = 1 , numberOfSpecies

        ! Read density matrices    
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

        densUnit = 78
        ! densFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".densmatrix"
        open(unit = densUnit, file=trim(densFile), status="old", form="unformatted")

        labels(2) = Grid_instance(speciesID)%nameOfSpecies
        labels(1) = "DENSITY-MATRIX"
        densityMatrix(speciesID) =Matrix_getFromFile(unit=densUnit, rows= int(numberOfContractions,4), &
             columns=int(numberOfContractions,4), binary=.true., arguments=labels)

        close(unit=densUnit)

        ! Calculate density and gradients

        !Initialize
        call Vector_Constructor( Grid_instance(speciesID)%potential, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%density, Grid_instance(speciesID)%totalSize, 0.0_8)
        do dir=1,3
           call Vector_Constructor( Grid_instance(speciesID)%gradientPotential(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
           call Vector_Constructor( Grid_instance(speciesID)%densityGradient(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
        end do
        exchangeCorrelationEnergy(speciesID)=0.0_8

        call GridManager_getDensityGradientAtGrid( speciesID, densityMatrix(speciesID), Grid_instance(speciesID)%density, Grid_instance(speciesID)%densityGradient)

        ! Check density and gradient in z
        numberOfParticles(speciesID)=0.0_8
        do i=1,Grid_instance(speciesID)%totalSize
           numberOfParticles(speciesID)=numberOfParticles(speciesID)+Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4)
           ! if (Grid_instance(speciesID)%points%values(i,1) .eq. 0.0 .and. Grid_instance(speciesID)%points%values(i,2) .eq. 0.0)  then
           !    print*, Grid_instance(speciesID)%points%values(i,3), Grid_instance(speciesID)%density%values(i), Grid_instance(speciesID)%densityGradient(3)%values(i)
           ! end if
        end do

        if(trim(job).eq."FINAL_GRID") then
           write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the final grid: ", numberOfParticles(speciesID)
        else
           write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the SCF grid: ", numberOfParticles(speciesID)
        end if
     end do

     ! call Stopwatch_stop(lowdin_stopwatch)     
     ! write(*,"(A,F10.3,A4)") "** Calculating density and gradient:", lowdin_stopwatch%enlapsetTime ," (s)"

     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)
     
     ! Calculate energy density and potential for one species
     do speciesID = 1 , numberOfSpecies
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

        if( nameOfSpecies .eq. "E-"  ) then 
           call GridManager_getElectronicEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID))

        elseif( nameOfSpecies .eq. "E-ALPHA"  ) then !El potencial de BETA se calcula simultaneamente con ALPHA
           otherSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie="E-BETA" )
           call GridManager_getElectronicEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                otherSpeciesID, exchangeCorrelationEnergy(otherSpeciesID) )

        elseif (nameOfSpecies .eq. "E-BETA") then
           !Todo se hizo en el paso anterior

        else
           !There aren't more same species functionals implemented so far
        end if

        if(trim(job).eq."FINAL_GRID") then
           write (*,"(A50, F15.8)") trim(MolecularSystem_getNameOfSpecie(speciesID))//" Exchange correlation contribution: ", exchangeCorrelationEnergy(speciesID)
        end if

     end do

     ! Calculate energy density and potential for two species
     do speciesID = 1 , numberOfSpecies-1
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

        do otherSpeciesID = speciesID+1 , numberOfSpecies
           nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)

           auxEnergy=0.0_8
           otherAuxEnergy=0.0_8
           otherElectronAuxEnergy=0.0_8
           if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then
              !Nada, todo se hace como si fuera una sola especie

              cycle

           elseif ( nameOfSpecies .eq. "E-" .and. &
                (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
              !Closed shell electron and other species terms

              call GridManager_getInterspeciesEnergyAndPotentialAtGrid( speciesID, auxEnergy, &
                   otherSpeciesID, otherAuxEnergy )

              exchangeCorrelationEnergy(speciesID)=exchangeCorrelationEnergy(speciesID)+auxEnergy/2
              exchangeCorrelationEnergy(otherSpeciesID)=exchangeCorrelationEnergy(otherSpeciesID)+otherAuxEnergy/2

              if(trim(job).eq."FINAL_GRID") write (*,"(A50, F15.8)") trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ", auxEnergy

           elseif ( nameOfSpecies .eq. "E-ALPHA" .and. &
                (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then

              otherElectronID=MolecularSystem_getSpecieID("E-BETA")

              call GridManager_getInterspeciesEnergyAndPotentialAtGrid( speciesID, auxEnergy, &
                   otherSpeciesID, otherAuxEnergy, &
                   otherElectronID, otherElectronAuxEnergy )
              
              exchangeCorrelationEnergy(speciesID)=exchangeCorrelationEnergy(speciesID)+auxEnergy/2
              exchangeCorrelationEnergy(otherSpeciesID)=exchangeCorrelationEnergy(otherSpeciesID)+otherAuxEnergy/2
              exchangeCorrelationEnergy(otherElectronID)=exchangeCorrelationEnergy(otherElectronID)+otherElectronAuxEnergy/2
              
              if(trim(job).eq."FINAL_GRID") then
                 write (*,"(A50, F15.8)") trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ",  auxEnergy
                 write (*,"(A50, F15.8)") trim("E-BETA")//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ",  otherElectronAuxEnergy
              end if
              
           else

              cycle
              !There aren't more different species functionals implemented so far
           end if
        end do
     end do

     if(trim(job).eq."FINAL_GRID") then
        print *, ""
        write (*,"(A50, F15.8)") "Exchange correlation energy with the final grid: ", sum(exchangeCorrelationEnergy)

        print *, ""
        print *, "Contact density in the final grid"
        print *, ""

        do speciesID = 1 , numberOfSpecies-1
           nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
           do otherSpeciesID = speciesID+1 , numberOfSpecies
              nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)

              if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then

                 !Nada

              elseif ( nameOfSpecies .eq. "E-" .and. &
                   (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
                 !Closed shell electron and other species terms

                 call GridManager_getContactDensity( speciesID, otherSpeciesID )

              elseif ( nameOfSpecies .eq. "E-ALPHA" .and. &
                   (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
                 !Open shell Electron and other species terms

                 otherElectronID=MolecularSystem_getSpecieID("E-BETA")

                 call GridManager_getContactDensity( speciesID, otherSpeciesID, otherElectronID )

              end if

           end do
        end do

     end if
        
     ! call Stopwatch_stop(lowdin_stopwatch)    
     ! write(*,"(A,F10.3,A4)") "** Calculating energy and potential:", lowdin_stopwatch%enlapsetTime ," (s)"

     do speciesID = 1 , numberOfSpecies
        numberOfContractions=MolecularSystem_getTotalNumberOfContractions( speciesID )
        call Matrix_constructor(exchangeCorrelationMatrix(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
        call GridManager_buildExchangeCorrelationMatrix(speciesID, exchangeCorrelationMatrix(speciesID))
     end do
     
     ! Write results to file
     do speciesID = 1 , numberOfSpecies

        excUnit = 79
        excFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".excmatrix"
        open(unit = excUnit, file=trim(excFile), status="replace", form="unformatted")

        ! print *, Grid_instance(speciesID)%nameOfSpecies
        ! print *, speciesID, numberOfParticles(speciesID), exchangeCorrelationEnergy(speciesID)
        ! call Matrix_show(exchangeCorrelationMatrix(speciesID))

        labels(2) = Grid_instance(speciesID)%nameOfSpecies
        labels(1) = "NUMBER-OF-PARTICLES"
        call Vector_writeToFile(unit=excUnit, binary=.true., value=numberOfParticles(speciesID), arguments= labels )

        labels(1) = "EXCHANGE-CORRELATION-ENERGY"
        call Vector_writeToFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(speciesID), arguments= labels )

        labels(1) = "EXCHANGE-CORRELATION-MATRIX"
        call Matrix_writeToFile( exchangeCorrelationMatrix(speciesID), unit=excUnit, binary=.true., arguments = labels(1:2) )

        close(unit=excUnit)

     end do
     
  end if
     
end program DFT

  
