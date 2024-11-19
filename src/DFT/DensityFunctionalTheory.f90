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
  use Functional_
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
  subroutine DensityFunctionalTheory_buildSCFGrid(scfGrids,scfGridsCommonPoints,exactExchangeFractions,system)
    implicit none
    type(Grid) :: scfGrids(:), scfGridsCommonPoints(:,:)
    real(8), optional :: exactExchangeFractions(*)
    type(MolecularSystem), optional, target :: system
    
    type(Functional), allocatable :: Functionals(:,:)
    type(MolecularSystem), pointer :: molSys
    integer :: speciesID,numberOfSpecies

    !!Start time
    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(molSys)

    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       print *, ""
       print *, "--------------------------------------------------------------------------------------"
       print *, "|---------------------Building SCF Integration Grids---------------------------------|"
       print *, "--------------------------------------------------------------------------------------"
       print *, "Euler-Maclaurin radial grids - Lebedev angular grids"
       print *, ""
    end if

    call GridManager_buildGrids(scfGrids,scfGridsCommonPoints,"INITIAL",molSys )

    allocate(Functionals(numberOfSpecies,numberOfSpecies))
    call Functional_createFunctionals(Functionals,numberOfSpecies,molSys)

    if(CONTROL_instance%PRINT_LEVEL .gt. 0) call Functional_show(Functionals)
    if(CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_writeGrids(scfGrids,scfGridsCommonPoints,Functionals,"INITIAL")
       call GridManager_atomicOrbitals(scfGrids,scfGridsCommonPoints,"WRITE","INITIAL" )
    else
       call GridManager_atomicOrbitals(scfGrids,scfGridsCommonPoints,"COMPUTE","INITIAL" )
    end if

    if(present(exactExchangeFractions)) then
       do speciesID=1, MolecularSystem_getNumberOfQuantumSpecies(molSys)
          exactExchangeFractions(speciesID)=Functional_getExchangeFraction(Functionals,speciesID)
       end do
    end if
    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** Building and writing grids and atomic orbitals:", lowdin_stopwatch%enlapsetTime ," (s)"

  end subroutine DensityFunctionalTheory_buildSCFGrid

  subroutine DensityFunctionalTheory_buildFinalGrid(finalGrids,finalGridsCommonPoints,system)
    implicit none
    type(Grid) :: finalGrids(:), finalGridsCommonPoints(:,:)
    type(MolecularSystem), optional, target :: system

    type(Functional), allocatable :: Functionals(:,:)
    type(MolecularSystem), pointer :: molSys
    integer :: numberOfSpecies

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(molSys)
    
    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       print *, ""
       print *, "--------------------------------------------------------------------------------------"
       print *, "|-------------------Building Final Integration Grids---------------------------------|"
       print *, "--------------------------------------------------------------------------------------"
       print *, "Euler-Maclaurin radial grids - Lebedev angular grids"
       print *, ""
    end if
    call GridManager_buildGrids(finalGrids,finalGridsCommonPoints,"FINAL",molSys)

    allocate(Functionals(numberOfSpecies,numberOfSpecies))
    call Functional_createFunctionals(Functionals,numberOfSpecies,molSys)

    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_writeGrids(finalGrids,finalGridsCommonPoints,Functionals,"FINAL" )
       call GridManager_atomicOrbitals(finalGrids,finalGridsCommonPoints,"WRITE","FINAL" )
    else
       call GridManager_atomicOrbitals(finalGrids,finalGridsCommonPoints,"COMPUTE","FINAL" )
    end if

  end subroutine DensityFunctionalTheory_buildFinalGrid

  subroutine DensityFunctionalTheory_SCFDFT(scfGrids,scfGridsCommonPoints,densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
    implicit none
    type(Grid) :: scfGrids(:), scfGridsCommonPoints(:,:)
    type(Matrix), intent(in) :: densityMatrix(*) !IN
    type(Matrix) :: exchangeCorrelationMatrix(*) !OUT
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    real(8) :: numberOfParticles(*) !OUT

    type(Functional), allocatable :: Functionals(:,:)
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

    numberOfSpecies=size(scfGrids(:))

    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_readGrids(scfGrids,scfGridsCommonPoints,"INITIAL")
       call GridManager_atomicOrbitals(scfGrids,scfGridsCommonPoints,"READ", "INITIAL" )
    end if
    !!Start time
    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)
    allocate(Functionals(numberOfSpecies,numberOfSpecies))
    call Functional_createFunctionals(Functionals,numberOfSpecies,scfGrids(1)%molSys)

    call DensityFunctionalTheory_calculateDensityAndGradients(scfGrids,scfGridsCommonPoints,densityMatrix,numberOfParticles)

    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** Calculating density and gradient:", lowdin_stopwatch%enlapsetTime ," (s)"

    ! call Stopwatch_constructor(lowdin_stopwatch)
    ! call Stopwatch_start(lowdin_stopwatch)

    call DensityFunctionalTheory_calculateEnergyDensity(scfGrids,scfGridsCommonPoints,Functionals,exchangeCorrelationEnergy)

    !!In the final iteration we don't update the exchange correlation matrix to save time
    do speciesID = 1 , numberOfSpecies
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions( speciesID, scfGrids(speciesID)%molSys )
       call Matrix_constructor(exchangeCorrelationMatrix(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
       call GridManager_buildExchangeCorrelationMatrix(scfGrids,scfGridsCommonPoints,speciesID, exchangeCorrelationMatrix(speciesID))
    end do

    ! call Stopwatch_stop(lowdin_stopwatch)    
    ! write(*,"(A,F10.3,A4)") "** Calculating energy and potential:", lowdin_stopwatch%enlapsetTime ," (s)"
  end subroutine DensityFunctionalTheory_SCFDFT

  subroutine DensityFunctionalTheory_finalDFT(finalGrids,finalGridsCommonPoints,densityMatrix, exchangeCorrelationMatrix, exchangeCorrelationEnergy, numberOfParticles)
    implicit none
    type(Grid) :: finalGrids(:), finalGridsCommonPoints(:,:)
    type(Matrix) :: densityMatrix(*) !IN
    type(Matrix) :: exchangeCorrelationMatrix(*) !OUT
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    real(8) :: numberOfParticles(*) !OUT

    type(Functional), allocatable :: Functionals(:,:)
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

    numberOfSpecies=size(finalGrids(:))

    !print scf grid information for comparison
    if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) then
       do speciesID = 1 , numberOfSpecies
          write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecies(speciesID, finalGrids(speciesID)%molSys))//" particles in the SCF grid: ", numberOfParticles(speciesID)
       end do
       print *, ""
       write (*,"(A50, F15.8)") "Exchange-correlation energy with the SCF grid: ", sum(exchangeCorrelationEnergy%values)
       print *, ""
    end if

    if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
       call GridManager_readGrids(finalGrids,finalGridsCommonPoints,"FINAL" )
       call GridManager_atomicOrbitals(finalGrids,finalGridsCommonPoints,"READ", "FINAL" )
    end if

    allocate(Functionals(numberOfSpecies,numberOfSpecies))
    call Functional_createFunctionals(Functionals,numberOfSpecies,finalGrids(1)%molSys)
    
    call DensityFunctionalTheory_calculateDensityAndGradients(finalGrids,finalGridsCommonPoints,densityMatrix,numberOfParticles)

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

    call DensityFunctionalTheory_calculateEnergyDensity(finalGrids,finalGridsCommonPoints,Functionals,exchangeCorrelationEnergy)

    !print scf grid information for comparison
    if(CONTROL_instance%PRINT_LEVEL .gt. 0 ) then
       do speciesID = 1 , numberOfSpecies
          write (*,"(A50 F15.8)") "Number of "//trim(MolecularSystem_getNameOfSpecies(speciesID, finalGrids(speciesID)%molSys))//" particles in the final grid: ", numberOfParticles(speciesID)
       end do
       print *, ""
       write (*,"(A50, F15.8)") "Exchange-correlation energy with the final grid: ", sum(exchangeCorrelationEnergy%values)
       print *, ""
    end if
    do speciesID = 1 , numberOfSpecies-1
       nameOfSpecies=finalGrids(speciesID)%nameOfSpecies
       do otherSpeciesID = speciesID+1 , numberOfSpecies
          nameOfOtherSpecies=finalGrids(otherSpeciesID)%nameOfSpecies

          if ( nameOfSpecies .eq. "E-" ) then
             ! if ( nameOfSpecies .eq. "E-" .and. nameOfOtherSpecies .eq. "POSITRON" ) then
             !Closed shell electron and other species terms

             call GridManager_getContactDensity(finalGrids,finalGridsCommonPoints,speciesID, otherSpeciesID )

          elseif ( nameOfSpecies .eq. "E-ALPHA" ) then
             ! elseif ( nameOfSpecies .eq. "E-ALPHA" .and. nameOfOtherSpecies .eq. "POSITRON" ) then
             !Open shell Electron and other species terms

             otherElectronID=MolecularSystem_getSpecieID("E-BETA",finalGrids(speciesID)%molSys)

             call GridManager_getContactDensity(finalGrids,finalGridsCommonPoints,speciesID, otherSpeciesID, otherElectronID )

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

  subroutine DensityFunctionalTheory_calculateDensityAndGradients(Grid_instance,GridCommonPoints,densityMatrix,numberOfParticles)
    implicit none
    type(Grid) :: Grid_instance(:)
    type(Grid) :: GridCommonPoints(:,:)
    type(Matrix) :: densityMatrix(*) !IN
    real(8) :: numberOfParticles(*) !OUT
    integer :: numberOfSpecies
    integer :: speciesID
    integer :: i,dir

    numberOfSpecies=size(Grid_instance(:))
    do speciesID = 1 , numberOfSpecies

       ! Calculate density and gradients

       !Initialize
       call Vector_Constructor( Grid_instance(speciesID)%potential, Grid_instance(speciesID)%totalSize, 0.0_8)
       call Vector_Constructor( Grid_instance(speciesID)%density, Grid_instance(speciesID)%totalSize, 0.0_8)
       do dir=1,3
          call Vector_Constructor( Grid_instance(speciesID)%gradientPotential(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
          call Vector_Constructor( Grid_instance(speciesID)%densityGradient(dir), Grid_instance(speciesID)%totalSize, 0.0_8)
       end do

       call GridManager_getDensityGradientAtGrid(Grid_instance,GridCommonPoints, speciesID, densityMatrix(speciesID), Grid_instance(speciesID)%density, Grid_instance(speciesID)%densityGradient)

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

  subroutine DensityFunctionalTheory_calculateEnergyDensity(Grid_instance,GridCommonPoints,Functionals,exchangeCorrelationEnergy)
    type(Grid) :: Grid_instance(:)
    type(Grid) :: GridCommonPoints(:,:)
    type(Functional) :: Functionals(:,:)
    type(Matrix) :: exchangeCorrelationEnergy !OUT
    integer :: numberOfSpecies
    integer :: speciesID, otherSpeciesID, otherElectronID
    character(50) :: nameOfSpecies, nameOfOtherSpecies

    numberOfSpecies=size(Grid_instance(:))
    exchangeCorrelationEnergy%values(:,:)=0.0_8
    ! Calculate energy density and potential for one species
    do speciesID = 1 , numberOfSpecies
       nameOfSpecies=Grid_instance(speciesID)%nameOfSpecies

       if( nameOfSpecies .eq. "E-"  ) then 
          call GridManager_getElectronicEnergyAndPotentialAtGrid(Grid_instance,GridCommonPoints,Functionals, speciesID, exchangeCorrelationEnergy%values(speciesID,speciesID))

       elseif( nameOfSpecies .eq. "E-ALPHA"  ) then !El potencial de BETA se calcula simultaneamente con ALPHA
          otherSpeciesID = MolecularSystem_getSpecieID( "E-BETA", Grid_instance(speciesID)%molSys )
          call GridManager_getElectronicEnergyAndPotentialAtGrid(Grid_instance,GridCommonPoints,Functionals, speciesID, exchangeCorrelationEnergy%values(speciesID,speciesID), &
               otherSpeciesID, exchangeCorrelationEnergy%values(otherSpeciesID,otherSpeciesID) )

       elseif (nameOfSpecies .eq. "E-BETA") then
          !Todo se hizo en el paso anterior

       else
          !There aren't more same species functionals implemented so far
       end if

       ! write (*,"(A50, F15.8)") trim(MolecularSystem_getNameOfSpecies(speciesID))//" Exchange-correlation contribution: ", exchangeCorrelationEnergy(speciesID,speciesID)

    end do

    ! Calculate energy density and potential for two species
    do speciesID = 1 , numberOfSpecies-1
       nameOfSpecies=Grid_instance(speciesID)%nameOfSpecies

       do otherSpeciesID = speciesID+1 , numberOfSpecies
          nameOfOtherSpecies=Grid_instance(otherSpeciesID)%nameOfSpecies

          if (nameOfSpecies .eq. "E-ALPHA" .and. nameOfSpecies .eq. "E-BETA") then
             !Nada, todo se hace como si fuera una sola especie

             cycle

          elseif ( nameOfSpecies .eq. "E-" .and. &
               (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then
             !Closed shell electron and other species terms

             call GridManager_getInterspeciesEnergyAndPotentialAtGrid(Grid_instance, GridCommonPoints,Functionals, speciesID, otherSpeciesID, exchangeCorrelationEnergy%values(speciesID,otherSpeciesID) )

             ! write (*,"(A50, F15.8)") trim(nameOfSpecies)//"/"//trim(nameOfOtherSpecies)//" Correlation contribution: ", exchangeCorrelationEnergy(speciesID,otherSpeciesID)


          elseif ( nameOfSpecies .eq. "E-ALPHA" .and. &
               (nameOfOtherSpecies .ne. "E-" .and. nameOfOtherSpecies .ne. "E-ALPHA" .and. nameOfOtherSpecies .ne. "E-BETA") ) then

             otherElectronID=MolecularSystem_getSpecieID("E-BETA",Grid_instance(speciesID)%molSys)

             call GridManager_getInterspeciesEnergyAndPotentialAtGrid(Grid_instance, GridCommonPoints,Functionals, speciesID, otherSpeciesID, exchangeCorrelationEnergy%values(speciesID,otherSpeciesID), &
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
