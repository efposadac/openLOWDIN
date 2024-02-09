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
!! @brief DFT callable routines for memory and direct calculations.
!!        This module allows calculations in the APMO-DFT framework
!! @author F. Moncada
!!
!! <b> Creation date : </b> 2023
!!
!! <b> History: </b>
!!
!!
module DensityFunctionalTheory_
  use CONTROL_
  use MolecularSystem_
  use GridManager_
  use String_
  use Exception_
  use omp_lib
  implicit none

  public :: &
       DensityFunctionalTheory_buildSCFGrid

contains

  !>
  !! @brief Builds a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine DensityFunctionalTheory_buildSCFGrid(exactExchangeFractions)
    implicit none
    real(8), optional :: exactExchangeFractions(*)
    integer :: speciesID

    !!Start time
    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       print *, ""
       print *, "--------------------------------------------------------------------------------------"
       print *, "|---------------------Building SCF Integration Grids---------------------------------|"
       print *, "--------------------------------------------------------------------------------------"
       print *, "Euler-Maclaurin radial grids - Lebedev angular grids"
       print *, ""
    end if

    call GridManager_buildGrids( "INITIAL" )
    if(CONTROL_instance%GRID_STORAGE .ne. "DISK") call Functional_createFunctionals( )
    if(CONTROL_instance%PRINT_LEVEL .gt. 0) call Functional_show( )
    if(CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_writeGrids( "INITIAL" )
       call GridManager_atomicOrbitals( "WRITE","INITIAL" )
    else
       call GridManager_atomicOrbitals( "COMPUTE","INITIAL" )
    end if

    if(present(exactExchangeFractions)) then
       do speciesID=1, MolecularSystem_getNumberOfQuantumSpecies()
          exactExchangeFractions(speciesID)=Functional_getExchangeFraction(speciesID)
       end do
    end if
    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** Building and writing grids and atomic orbitals:", lowdin_stopwatch%enlapsetTime ," (s)"

  end subroutine DensityFunctionalTheory_buildSCFGrid

  subroutine DensityFunctionalTheory_buildFinalGrid()
    implicit none
    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       print *, ""
       print *, "--------------------------------------------------------------------------------------"
       print *, "|-------------------Building Final Integration Grids---------------------------------|"
       print *, "--------------------------------------------------------------------------------------"
       print *, "Euler-Maclaurin radial grids - Lebedev angular grids"
       print *, ""
    end if
    call GridManager_buildGrids( "FINAL" )
    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_writeGrids( "FINAL" )
       call GridManager_atomicOrbitals( "WRITE","FINAL" )
    else
       call GridManager_atomicOrbitals( "COMPUTE","FINAL" )
    end if

  end subroutine DensityFunctionalTheory_buildFinalGrid

  subroutine DensityFunctionalTheory_SCFDFT(densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
    implicit none
    type(Matrix), intent(in) :: densityMatrix(*) !IN
    type(Matrix) :: exchangeCorrelationMatrix(*) !OUT
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    real(8) :: numberOfParticles(*) !OUT

    type(Matrix), allocatable :: overlapMatrix(:)
    character(50) ::  labels(2)
    integer :: densUnit, excUnit

    integer :: numberOfSpecies
    integer :: numberOfContractions
    integer :: speciesID, otherSpeciesID, otherElectronID
    character(50) :: nameOfSpecies, nameOfOtherSpecies

    integer :: i, u, dir
    real(8) :: sumCheck, auxEnergy, otherAuxEnergy, otherElectronAuxEnergy
    real(8) :: time1, time2, time3

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_readGrids( "INITIAL")
       call GridManager_atomicOrbitals( "READ", "INITIAL" )
    end if
    !!Start time
    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    call DensityFunctionalTheory_calculateDensityAndGradients(densityMatrix,numberOfParticles)

    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** Calculating density and gradient:", lowdin_stopwatch%enlapsetTime ," (s)"

    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    call DensityFunctionalTheory_calculateEnergyDensity(exchangeCorrelationEnergy)

    !!In the final iteration we don't update the exchange correlation matrix to save time
    do speciesID = 1 , numberOfSpecies
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions( speciesID )
       call Matrix_constructor(exchangeCorrelationMatrix(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
       call GridManager_buildExchangeCorrelationMatrix(speciesID, exchangeCorrelationMatrix(speciesID))
    end do

    ! call Stopwatch_stop(lowdin_stopwatch)    
    ! write(*,"(A,F10.3,A4)") "** Calculating energy and potential:", lowdin_stopwatch%enlapsetTime ," (s)"
  end subroutine DensityFunctionalTheory_SCFDFT

  subroutine DensityFunctionalTheory_finalDFT(densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
    implicit none
    type(Matrix) :: densityMatrix(*) !IN
    type(Matrix) :: exchangeCorrelationMatrix(*) !OUT
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    real(8) :: numberOfParticles(*) !OUT

    type(Matrix), allocatable :: overlapMatrix(:)
    character(100) :: excFile
    character(50) ::  labels(2)
    integer :: densUnit, excUnit

    integer :: numberOfSpecies
    integer :: numberOfContractions
    integer :: speciesID, otherSpeciesID, otherElectronID
    character(50) :: nameOfSpecies, nameOfOtherSpecies

    integer :: i, u, dir
    real(8) :: sumCheck, auxEnergy, otherAuxEnergy, otherElectronAuxEnergy
    real(8) :: time1, time2, time3

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

    !print scf grid information for comparison
    if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) then
       do speciesID = 1 , numberOfSpecies
          write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the SCF grid: ", numberOfParticles(speciesID)
       end do
       print *, ""
       write (*,"(A50, F15.8)") "Exchange-correlation energy with the SCF grid: ", sum(exchangeCorrelationEnergy%values)
       print *, ""
    end if

    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_readGrids( "FINAL" )
       call GridManager_atomicOrbitals( "READ", "FINAL" )
    end if

    call DensityFunctionalTheory_calculateDensityAndGradients(densityMatrix,numberOfParticles)

    !!Start time
    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    if(CONTROL_instance%BETA_FUNCTION .eq. "PsBeta" ) then
       CONTROL_instance%BETA_FUNCTION = "PsBetaMax"
       if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) &
            print *, "We are changing the PsBeta with cutoff to the PsBeta with max(pe,pp)"
    end if

    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** Calculating density and gradient:", lowdin_stopwatch%enlapsetTime ," (s)"

    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    call DensityFunctionalTheory_calculateEnergyDensity(exchangeCorrelationEnergy)

    !print scf grid information for comparison
    if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) then
       do speciesID = 1 , numberOfSpecies
          write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecie(speciesID))//" particles in the final grid: ", numberOfParticles(speciesID)
       end do
       print *, ""
       write (*,"(A50, F15.8)") "Exchange-correlation energy with the final grid: ", sum(exchangeCorrelationEnergy%values)
       print *, ""
    end if
    do speciesID = 1 , numberOfSpecies-1
       nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
       do otherSpeciesID = speciesID+1 , numberOfSpecies
          nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)

          if ( nameOfSpecies .eq. "E-" ) then
             ! if ( nameOfSpecies .eq. "E-" .and. nameOfOtherSpecies .eq. "POSITRON" ) then
             !Closed shell electron and other species terms

             call GridManager_getContactDensity( speciesID, otherSpeciesID )

          elseif ( nameOfSpecies .eq. "E-ALPHA" ) then
             ! elseif ( nameOfSpecies .eq. "E-ALPHA" .and. nameOfOtherSpecies .eq. "POSITRON" ) then
             !Open shell Electron and other species terms

             otherElectronID=MolecularSystem_getSpecieID("E-BETA")

             call GridManager_getContactDensity( speciesID, otherSpeciesID, otherElectronID )

          end if

       end do
    end do

    ! print *, ""
    ! print *, "Density-point charges expected distances"
    ! print *, ""
    ! do speciesID = 1 , numberOfSpecies
    !    call GridManager_getExpectedDistances( speciesID )
    ! end do
    ! print *, ""


    ! call Stopwatch_stop(lowdin_stopwatch)    
    ! write(*,"(A,F10.3,A4)") "** Calculating energy and potential:", lowdin_stopwatch%enlapsetTime ," (s)"
    if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) print *, "END DFT FINAL GRID INTEGRATION"

  end subroutine DensityFunctionalTheory_finalDFT

  subroutine DensityFunctionalTheory_calculateDensityAndGradients(densityMatrix,numberOfParticles)
    implicit none
    type(Matrix) :: densityMatrix(*) !IN
    real(8) :: numberOfParticles(*) !OUT
    integer :: numberOfSpecies
    integer :: speciesID
    integer :: i,dir

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    do speciesID = 1 , numberOfSpecies

       ! Calculate density and gradients

       !Initialize
       call Vector_Constructor( Grid_instance(speciesID)%potential, Grid_instance(speciesID)%totalSize, 0.0_8)
       call Vector_Constructor( Grid_instance(speciesID)%density, Grid_instance(speciesID)%totalSize, 0.0_8)
       do dir=1,3
          call Vector_Constructor( Grid_instance(speciesID)%gradientPotential(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
          call Vector_Constructor( Grid_instance(speciesID)%densityGradient(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
       end do

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

  end subroutine DensityFunctionalTheory_calculateDensityAndGradients

  subroutine DensityFunctionalTheory_calculateEnergyDensity(exchangeCorrelationEnergy)
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    integer :: numberOfSpecies
    integer :: speciesID, otherSpeciesID, otherElectronID
    character(50) :: nameOfSpecies, nameOfOtherSpecies

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    exchangeCorrelationEnergy%values(:,:)=0.0_8
    ! Calculate energy density and potential for one species
    do speciesID = 1 , numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

       if( nameOfSpecies .eq. "E-"  ) then 
          call GridManager_getElectronicEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy%values(speciesID,speciesID))

       elseif( nameOfSpecies .eq. "E-ALPHA"  ) then !El potencial de BETA se calcula simultaneamente con ALPHA
          otherSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie="E-BETA" )
          call GridManager_getElectronicEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy%values(speciesID,speciesID), &
               otherSpeciesID, exchangeCorrelationEnergy%values(otherSpeciesID,otherSpeciesID) )

       elseif (nameOfSpecies .eq. "E-BETA") then
          !Todo se hizo en el paso anterior

       else
          !There aren't more same species functionals implemented so far
       end if

       ! write (*,"(A50, F15.8)") trim(MolecularSystem_getNameOfSpecie(speciesID))//" Exchange-correlation contribution: ", exchangeCorrelationEnergy(speciesID,speciesID)

    end do

    ! Calculate energy density and potential for two species
    do speciesID = 1 , numberOfSpecies-1
       nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

       do otherSpeciesID = speciesID+1 , numberOfSpecies
          nameOfOtherSpecies=MolecularSystem_getNameOfSpecie(otherSpeciesID)

          if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then
             !Nada, todo se hace como si fuera una sola especie

             cycle

          elseif ( nameOfSpecies .eq. "E-" .and. &
               (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
             !Closed shell electron and other species terms

             call GridManager_getInterspeciesEnergyAndPotentialAtGrid( speciesID, otherSpeciesID, exchangeCorrelationEnergy%values(speciesID,otherSpeciesID) )

             ! write (*,"(A50, F15.8)") trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ", exchangeCorrelationEnergy(speciesID,otherSpeciesID)


          elseif ( nameOfSpecies .eq. "E-ALPHA" .and. &
               (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then

             otherElectronID=MolecularSystem_getSpecieID("E-BETA")

             call GridManager_getInterspeciesEnergyAndPotentialAtGrid( speciesID, otherSpeciesID, exchangeCorrelationEnergy%values(speciesID,otherSpeciesID), &
                  otherElectronID, exchangeCorrelationEnergy%values(otherElectronID,otherSpeciesID) )

             ! write (*,"(A50, F15.8)") trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ",  exchangeCorrelationEnergy(speciesID,otherSpeciesID)
             ! write (*,"(A50, F15.8)") trim("E-BETA")//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ",  exchangeCorrelationEnergy(otherElectronID,otherSpeciesID)
          else

             cycle
             !There aren't more different species functionals implemented so far
          end if
       end do
    end do
  end subroutine DensityFunctionalTheory_calculateEnergyDensity

end module DensityFunctionalTheory_
