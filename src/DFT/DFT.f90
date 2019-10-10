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
!! @brief HF and APMO-HF program.
!!        This module allows to make calculations in the APMO-HF framework
!! @author  E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-05-09
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-09-17 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos (RHF.f90)
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adpata el modulo para incluirlo en LOWDIN (RHF.90)
!!   - <tt> 2013-05-09 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrite the module as a program, handles closed and open shell systems
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
  integer :: speciesID, otherSpeciesID
  character(50) :: nameOfSpecies, nameOfOtherSpecies

  integer :: i
  real(8) :: sumCheck, auxEnergy, otherAuxEnergy
  real(8) :: time1, time2, time3
  
  job = ""
  call get_command_argument(1,value=job)
  job = trim(String_getUppercase(job))

  ! write(*,"(A)") trim(job)
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  select case(trim(job))         

  case("INITIALIZE")

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
     call GridManager_writeAtomicOrbitals( )

     ! call Stopwatch_stop(lowdin_stopwatch)     
     ! write(*,"(A,F10.3,A4)") "** Building and writing grids and atomic orbitals:", lowdin_stopwatch%enlapsetTime ," (s)"


!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
          
  case("BUILD_MATRICES")

     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)

     call Functional_createFunctionals( )
     call GridManager_readGrids( )

     numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
     allocate( densityMatrix(numberOfSpecies) , numberOfParticles(numberOfSpecies), exchangeCorrelationMatrix(numberOfSpecies), exchangeCorrelationEnergy(numberOfSpecies)) !overlapMatrix(numberOfSpecies) 

     ! Read density matrices    
     do speciesID = 1 , MolecularSystem_getNumberOfQuantumSpecies()

        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

        densUnit = 78
        densFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".densmatrix"
        open(unit = densUnit, file=trim(densFile), status="old", form="unformatted")

        labels(2) = Grid_instance(speciesID)%nameOfSpecies
        labels(1) = "DENSITY-MATRIX"
        densityMatrix(speciesID) =Matrix_getFromFile(unit=densUnit, rows= int(numberOfContractions,4), &
             columns=int(numberOfContractions,4), binary=.true., arguments=labels)

        close(unit=densUnit)

        ! print *, "density matrix for", speciesID
        ! call Matrix_show(densityMatrix(speciesID))
               
        ! wfnUnit = 300
        ! wfnFile = "lowdin.wfn"

        ! open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

        ! labels(1) = "OVERLAP"
        ! overlapMatrix(speciesID) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
        !      columns= int(numberOfContractions,4), binary=.true., arguments=labels)

        ! close(unit=wfnUnit)

     end do

     !stop time
     ! call Stopwatch_stop(lowdin_stopwatch)     
     ! write(*,"(A,F10.3,A4)") "** Reading matrices:", lowdin_stopwatch%enlapsetTime ," (s)"

     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)

    ! Calculate density and gradients
     do speciesID = 1 , numberOfSpecies

        !Initialize
        call Vector_Constructor( Grid_instance(speciesID)%potential, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%sigmaPotential, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%density, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(1), Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(2), Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(3), Grid_instance(speciesID)%totalSize, 0.0_8)
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

     end do

     ! call Stopwatch_stop(lowdin_stopwatch)     
     ! write(*,"(A,F10.3,A4)") "** Calculating density and gradient:", lowdin_stopwatch%enlapsetTime ," (s)"

     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)

     
     ! Calculate energy density and potential for one species
     do speciesID = 1 , numberOfSpecies
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
        if( nameOfSpecies .eq. "E-ALPHA"  ) then !El potencial de BETA se calcula simultaneamente con ALPHA
           otherSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie="E-BETA" )

           call GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential,&
                otherSpeciesID, exchangeCorrelationEnergy(otherSpeciesID),&
                Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%sigmaPotential )
           
        elseif (nameOfSpecies .eq. "E-BETA") then

           !Todo se hizo en el paso anterior

        else
           call GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential )

        end if

     end do
    
     ! Calculate energy density and potential for two species
     do speciesID = 1 , numberOfSpecies
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
        do otherSpeciesID = speciesID+1 , numberOfSpecies
           nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)
           if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then

              !Nada, todo se hace como si fuera una sola especie

           elseif ( (nameOfSpecies .eq. "E-" .or. nameOfSpecies .eq. "E-ALPHA" .or. nameOfSpecies .eq. "E-BETA") .and. &
                (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
              !Electron and other species terms

              call GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                   Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential,&
                   otherSpeciesID, exchangeCorrelationEnergy(otherSpeciesID), &
                   Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%sigmaPotential )

           end if

        end do
     end do


     ! call Stopwatch_stop(lowdin_stopwatch)    
     ! write(*,"(A,F10.3,A4)") "** Calculating energy and potential:", lowdin_stopwatch%enlapsetTime ," (s)"

     do speciesID = 1 , numberOfSpecies
        numberOfContractions=MolecularSystem_getTotalNumberOfContractions( speciesID )
        call Matrix_constructor(exchangeCorrelationMatrix(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
     end do

     do speciesID = 1 , numberOfSpecies
        numberOfContractions=MolecularSystem_getTotalNumberOfContractions( speciesID )
        call GridManager_buildExchangeCorrelationMatrix(speciesID, exchangeCorrelationMatrix(speciesID))
     end do
     
     ! Write results to file
     do speciesID = 1 , numberOfSpecies

        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

        excUnit = 79
        excFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".excmatrix"
        open(unit = excUnit, file=trim(excFile), status="replace", form="unformatted")


        ! print *, Grid_instance(speciesID)%nameOfSpecies
        ! print *, speciesID, exchangeCorrelationEnergy(speciesID)
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

!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
     
  case("FINAL_GRID")

     !!! IN THIS ROUTINE WE DO NOT UPDATE THE EXCHANGE CORRELATION MATRIX BECAUSE IT IS INTENDED TO BE A POST-SCF CALL
    print *, ""
    print *, "--------------------------------------------------------------------------------------"
    print *, "|---------------------------- DFT Final Integration ---------------------------------|"
    print *, "--------------------------------------------------------------------------------------"
    print *, ""

    !!Start time
     ! call Stopwatch_constructor(lowdin_stopwatch)
     ! call Stopwatch_start(lowdin_stopwatch)

     call GridManager_buildGrids( "FINAL" )
     call Functional_createFunctionals( )
     call GridManager_writeAtomicOrbitals( )

     numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
     allocate( densityMatrix(numberOfSpecies) , numberOfParticles(numberOfSpecies), exchangeCorrelationMatrix(numberOfSpecies), exchangeCorrelationEnergy(numberOfSpecies)) !overlapMatrix(numberOfSpecies) 

     ! Read density and exchange correlation matrices    
     do speciesID = 1 , MolecularSystem_getNumberOfQuantumSpecies()

        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

        densUnit = 78
        densFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".densmatrix"
        open(unit = densUnit, file=trim(densFile), status="old", form="unformatted")

        labels(2) = Grid_instance(speciesID)%nameOfSpecies
        labels(1) = "DENSITY-MATRIX"
        densityMatrix(speciesID) =Matrix_getFromFile(unit=densUnit, rows= int(numberOfContractions,4), &
             columns=int(numberOfContractions,4), binary=.true., arguments=labels)

        close(unit=densUnit)

        excUnit = 79
        excFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".excmatrix"
        open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")

        labels(1) = "NUMBER-OF-PARTICLES"
        call Vector_getFromFile(unit=excUnit, binary=.true., value=numberOfParticles(speciesID), arguments= labels )

        labels(1) = "EXCHANGE-CORRELATION-ENERGY"
        call Vector_getFromFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(speciesID), arguments= labels )

        labels(1) = "EXCHANGE-CORRELATION-MATRIX"
        exchangeCorrelationMatrix(speciesID)=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
             binary=.true., arguments=labels)

        close(unit=excUnit)

        write (*,"(A50 F12.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the SCF grid: ", numberOfParticles(speciesID)

     end do

     print *, ""
     write (*,"(A50, F12.8)") "Exchange correlation energy with the SCF grid: ", sum(exchangeCorrelationEnergy)
     print *, ""
     
     ! Calculate density and gradients
     do speciesID = 1 , numberOfSpecies

        !Initialize
        call Vector_Constructor( Grid_instance(speciesID)%potential, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%sigmaPotential, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%density, Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(1), Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(2), Grid_instance(speciesID)%totalSize, 0.0_8)
        call Vector_Constructor( Grid_instance(speciesID)%densityGradient(3), Grid_instance(speciesID)%totalSize, 0.0_8)
        exchangeCorrelationEnergy(speciesID)=0.0_8

        call GridManager_getDensityGradientAtGrid( speciesID, densityMatrix(speciesID), Grid_instance(speciesID)%density, Grid_instance(speciesID)%densityGradient)

        numberOfParticles(speciesID)=0.0_8
        do i=1,Grid_instance(speciesID)%totalSize
           numberOfParticles(speciesID)=numberOfParticles(speciesID)+Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4)
        end do
        
        write (*,"(A50 F12.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the final grid: ",  numberOfParticles(speciesID)

  end do
  

  !    print *, ""
     
     ! Calculate energy density and potential for one species
     do speciesID = 1 , numberOfSpecies
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
        if( nameOfSpecies .eq. "E-ALPHA"  ) then !El potencial de BETA se calcula simultaneamente con ALPHA
           otherSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie="E-BETA" )

           call GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential,&
                otherSpeciesID, exchangeCorrelationEnergy(otherSpeciesID),&
                Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%sigmaPotential )
           
        elseif (nameOfSpecies .eq. "E-BETA") then

           !Todo se hizo en el paso anterior

        else
           call GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy(speciesID), &
                Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential )

        end if

        write (*,"(A50, F12.8)") trim(MolecularSystem_getNameOfSpecie(speciesID))//" Exchange correlation contribution: ", exchangeCorrelationEnergy(speciesID)
        
     end do

    
     ! Calculate energy density and potential for two species
     do speciesID = 1 , numberOfSpecies-1
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
        do otherSpeciesID = speciesID+1 , numberOfSpecies
           nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)

           auxEnergy=0.0_8
           otherAuxEnergy=0.0_8

           if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then

              !Nada, todo se hace como si fuera una sola especie

           elseif ( (nameOfSpecies .eq. "E-" .or. nameOfSpecies .eq. "E-ALPHA" .or. nameOfSpecies .eq. "E-BETA") .and. &
                (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
              !Electron and other species terms

              call GridManager_getEnergyAndPotentialAtGrid( speciesID, auxEnergy, &
                   Grid_instance(speciesID)%potential, Grid_instance(speciesID)%sigmaPotential,&
                   otherSpeciesID, otherAuxEnergy, &
                   Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%sigmaPotential )

           end if
           
           write (*,"(A50, F12.8)") &
                trim(MolecularSystem_getNameOfSpecie(speciesID))//"/"//trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))//" Correlation contribution: ",  auxEnergy+otherAuxEnergy

           exchangeCorrelationEnergy(speciesID)=exchangeCorrelationEnergy(speciesID)+auxEnergy
           exchangeCorrelationEnergy(otherSpeciesID)=exchangeCorrelationEnergy(otherSpeciesID)+otherAuxEnergy

        end do
     end do

  !    print *, ""
     
  !    write (*,"(A50, F12.8)") "Exchange correlation energy with the final grid: ", sum(exchangeCorrelationEnergy)

     
     ! Write results to file
     ! do speciesID = 1 , numberOfSpecies

     !    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

     !    excUnit = 79
     !    excFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".excmatrix"
     !    open(unit = excUnit, file=trim(excFile), status="replace", form="unformatted")

     !    labels(2) = Grid_instance(speciesID)%nameOfSpecies
     !    labels(1) = "EXCHANGE-CORRELATION-ENERGY"
     !    call Vector_writeToFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(speciesID), arguments= labels )

        ! labels(1) = "EXCHANGE-CORRELATION-MATRIX"
        ! call Matrix_writeToFile( exchangeCorrelationMatrix(speciesID), unit=excUnit, binary=.true., arguments = labels(1:2) )

        ! write (*,"(A50, A10)") "Final Exchange correlation matrix for:", Grid_instance(speciesID)%nameOfSpecies
        ! call Matrix_show(exchangeCorrelationMatrix(speciesID))
        
     !    close(unit=excUnit)

     ! end do

!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************
!!*********************************************************************************************************************************

  case default

     write(*,*) "USAGE: lowdin-DFT.x job "
     write(*,*) "Where job can be: "
     write(*,*) "  initialize"
     write(*,*) "  build_matrices"
     stop "ERROR"

  end select

end program DFT

  
