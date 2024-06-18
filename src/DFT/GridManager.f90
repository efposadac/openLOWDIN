!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	PROF. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module manages the orbital and density represented in the DFT grids.
!! @author F. Moncada, 2017
module GridManager_

  use Matrix_
  use Grid_
  use Exception_
  use String_
  use MolecularSystem_
  use Functional_

  implicit none

  public :: &
       GridManager_buildGrids, &
       GridManager_writeGrids, &
       GridManager_readGrids, &
       GridManager_getDensityGradientAtGrid, &
       GridManager_getElectronicEnergyAndPotentialAtGrid, & !, &       GridManager_getEnergyFromGrid
       GridManager_getInterspeciesEnergyAndPotentialAtGrid, & !, &       GridManager_getEnergyFromGrid
       GridManager_atomicOrbitals, &
       GridManager_buildExchangeCorrelationMatrix, &
       GridManager_findCommonPoints, &
       GridManager_getContactDensity
  ! GridManager_getOrbitalGradientAtPoint, &
  ! GridManager_getOrbitalGradientMatrix!, &
  ! GridManager_getOrbitalMatrix
  ! GridManager_getOrbitalAtGrid, &
  ! GridManager_getDensityAtGrid, &


contains

  !>
  !! @brief Builds a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_buildGrids( type )
    implicit none
    character(*) :: type
    integer :: numberOfSpecies
    integer :: speciesID,otherSpeciesID
    character(50) :: labels(2) 
    character(100) ::   dftFile
    integer :: dftUnit


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate memory.
    if (.not. allocated(Grid_instance)) allocate(Grid_instance(numberOfSpecies))
    if (.not. allocated(GridsCommonPoints)) allocate(GridsCommonPoints(numberOfSpecies,numberOfSpecies))

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       call Grid_constructor(Grid_instance(speciesID), speciesID , type )
       
    end do

    do speciesID = 1, numberOfSpecies-1
       do otherSpeciesID = speciesID+1, numberOfSpecies

          call GridManager_findCommonPoints(Grid_instance(speciesID)%points,Grid_instance(speciesID)%totalSize,&
               Grid_instance(otherSpeciesID)%points,Grid_instance(otherSpeciesID)%totalSize,&
               GridsCommonPoints(speciesID,otherSpeciesID)%points,GridsCommonPoints(speciesID,otherSpeciesID)%totalSize)
       end do
    end do

  end subroutine GridManager_buildGrids
  !>
  !! @brief Writes a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_writeGrids( type )
    implicit none
    character(*) :: type
    integer :: numberOfSpecies
    integer :: speciesID,otherSpeciesID
    character(50) :: labels(2)
    character(100) ::   dftFile
    integer :: dftUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       !! Open file for dft
       dftUnit = 77
       if( trim(type) .eq. "INITIAL" ) then
          dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".grid"
       else if( trim(type) .eq. "FINAL" ) then
          dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".finalGrid"
       else
          STOP "ERROR At DFT program, requested an unknown grid type to writeGrids at GridManager"
       end if
       
       open(unit = dftUnit, file=trim(dftFile), status="replace", form="unformatted")

       labels(2) = Grid_instance(speciesID)%nameOfSpecies
       labels(1) = "GRID-SIZE"

       call Vector_writeToFile(unit=dftUnit, binary=.true., value=real(Grid_instance(speciesID)%totalSize,8), arguments= labels )

       !! This goes here for convenience only
       labels(1) = "EXACT-EXCHANGE-FRACTION" 

       call Vector_writeToFile(unit=dftUnit, binary=.true., value=Functional_getExchangeFraction(speciesID), arguments= labels )

       labels(1) = "INTEGRATION-GRID"
       call Matrix_writeToFile(Grid_instance(speciesID)%points, unit=dftUnit, binary=.true., arguments = labels(1:2) )

       ! print *, "Grid_instance(speciesID)%points", speciesID
       ! call Matrix_show (Grid_instance(speciesID)%points)

       close(unit=dftUnit)
    end do

    !! Writes common points
    do speciesID = 1, numberOfSpecies-1

       do otherSpeciesID = speciesID+1, numberOfSpecies

          dftUnit = 77
          if( trim(type) .eq. "INITIAL" ) then
             dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)//".commonGrid"
          else if( trim(type) .eq. "FINAL" ) then
             dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)//".commonFinalGrid"
          else
             STOP "ERROR At DFT program, requested an unknown grid type to writeGrids at GridManager"
          end if
          open(unit = dftUnit, file=trim(dftFile), status="replace", form="unformatted")

          labels(2) = trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)
          labels(1) = "GRID-SIZE"

          call Vector_writeToFile(unit=dftUnit, binary=.true., value=real(GridsCommonPoints(speciesID,otherSpeciesID)%totalSize,8), arguments= labels )

          labels(1) = "COMMON-POINTS"
          call Matrix_writeToFile(GridsCommonPoints(speciesID,otherSpeciesID)%points, unit=dftUnit, binary=.true., arguments = labels(1:2) )

          close(unit=dftUnit)
       end do

    end do
  end subroutine GridManager_writeGrids


  !>
  !! @brief Reads a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_readGrids( type )
    implicit none
    character(*) :: type
    integer :: numberOfSpecies
    integer :: speciesID,otherSpeciesID
    character(50) :: labels(2)
    character(100) ::   dftFile
    integer :: dftUnit
    real(8) :: auxVal

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate memory.
    if (.not. allocated(Grid_instance)) allocate(Grid_instance(numberOfSpecies))
    if (.not. allocated(GridsCommonPoints)) allocate(GridsCommonPoints(numberOfSpecies,numberOfSpecies))

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       Grid_instance(speciesID)%nameOfSpecies=trim(MolecularSystem_getNameOfSpecie(speciesID))
       !! Open file for dft
       dftUnit = 77
       if( trim(type) .eq. "INITIAL" ) then
          dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".grid"
       else if( trim(type) .eq. "FINAL" ) then
          dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".finalGrid"
       else
          STOP "ERROR At DFT program, requested an unknown grid type to readGrids at GridManager"
       end if

       open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

       labels(2) = Grid_instance(speciesID)%nameOfSpecies
       labels(1) = "GRID-SIZE"
       call Vector_getFromFile(unit=dftUnit, binary=.true., value=auxVal, arguments=labels)
       Grid_instance(speciesID)%totalSize=int(auxVal)

       labels(1) = "INTEGRATION-GRID"
       Grid_instance(speciesID)%points=Matrix_getFromFile(unit=dftUnit, rows= int(Grid_instance(speciesID)%totalSize,4), &
            columns=int(4,4), binary=.true., arguments=labels)

       ! print *, "grid recien leida"
       ! print *, size(Grid_instance(speciesID)%points%values)
       ! call Matrix_show (Grid_instance(speciesID)%points)

       close(unit=dftUnit)
    end do
    
    do speciesID = 1, numberOfSpecies-1
       
       do otherSpeciesID = speciesID+1, numberOfSpecies

          dftUnit = 77
          if( trim(type) .eq. "INITIAL" ) then
             dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)//".commonGrid"
          else if( trim(type) .eq. "FINAL" ) then
             dftFile = "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)//".commonFinalGrid"
          else
             STOP "ERROR At DFT program, requested an unknown grid type to readGrids at GridManager"
          end if
          open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

          labels(2) = trim(Grid_instance(speciesID)%nameOfSpecies)//trim(Grid_instance(otherSpeciesID)%nameOfSpecies)
          labels(1) = "GRID-SIZE"

          call Vector_getFromFile(unit=dftUnit, binary=.true., value=auxVal, arguments=labels)
          GridsCommonPoints(speciesID,otherSpeciesID)%totalSize=int(auxVal)

          labels(1) = "COMMON-POINTS"
          GridsCommonPoints(speciesID,otherSpeciesID)%points=Matrix_getFromFile(unit=dftUnit, rows= int(GridsCommonPoints(speciesID,otherSpeciesID)%totalSize,4), &
            columns=int(2,4), binary=.true., arguments=labels)
          
          close(unit=dftUnit)
       end do

    end do

  end subroutine GridManager_readGrids


  !>
  !! @brief Writes the values of all the atomic orbitals and their gradients in a set of coordinates to a file
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_atomicOrbitals( action, type )
    implicit none
    character(*) action !read, compute or write
    character(*) type

    integer :: numberOfSpecies
    integer :: totalNumberOfContractions
    integer :: speciesID
    integer :: gridSize
    integer :: mu,nu, point, index

    character(50) :: labels(2)
    character(100) ::   orbsFile
    integer :: orbsUnit

    type(Matrix) :: auxMatrix(4)
    integer :: i, j, k, g, u
    integer :: numberOfCartesiansOrbitals

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do speciesID = 1 , numberOfSpecies
       gridSize = Grid_instance(speciesID)%totalSize
       totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
       orbsUnit = 78
       if( trim(type) .eq. "INITIAL" ) then
          write( orbsFile, "(A,I0.4)") "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbitals"
       else if( trim(type) .eq. "FINAL" ) then
          write( orbsFile, "(A,I0.4)") "lowdin."//trim(Grid_instance(speciesID)%nameOfSpecies)//".finalOrbitals"
       else
          STOP "ERROR At DFT program, requested an unknown grid type to orbitals at GridManager"
       end if

       if (allocated(Grid_instance(speciesID)%orbitalsWithGradient)) then
          do u=1, size(Grid_instance(speciesID)%orbitalsWithGradient(:))
             call Matrix_destructor(Grid_instance(speciesID)%orbitalsWithGradient(u))
          end do
          deallocate(Grid_instance(speciesID)%orbitalsWithGradient)
       end if
       allocate(Grid_instance(speciesID)%orbitalsWithGradient(totalNumberOfContractions))

       do u=1, totalNumberOfContractions
          call Matrix_Constructor( Grid_instance(speciesID)%orbitalsWithGradient(u), int(gridSize,8), int(4,8), 0.0_8)
       end do

       if( trim(action) .eq. "READ") then
          open(unit = orbsUnit, file=trim(orbsFile), status="old", form="unformatted")
          do u = 1, totalNumberOfContractions
             write( labels(1), "(A,I0.4)") "ORBITAL_", u
             labels(2) = Grid_instance(speciesID)%nameOfSpecies

             Grid_instance(speciesID)%orbitalsWithGradient(u)=Matrix_getFromFile(unit=orbsUnit, rows= int(gridSize,4), &
                  columns= int(4,4), binary=.true., arguments=labels)

          end do
          close(unit=orbsUnit)

       else if (trim(action) .eq. "COMPUTE" .or. trim(action) .eq. "WRITE") then

          if(trim(action) .eq. "WRITE") open(unit = orbsUnit, file=trim(orbsFile), status="replace", form="unformatted")
          k=0
          do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
             do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)
                numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital

                call Matrix_constructor( auxMatrix(1), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !orbital
                call Matrix_constructor( auxMatrix(2), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dx
                call Matrix_constructor( auxMatrix(3), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dy
                call Matrix_constructor( auxMatrix(4), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dz
                
                call ContractedGaussian_getGradientAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), &
                     Grid_instance(speciesID)%points, gridSize, auxMatrix(1), auxMatrix(2), auxMatrix(3), auxMatrix(4))

                ! call Matrix_show(auxMatrix(1))

                ! print *, "wololooooo"

                do j = 1, numberOfCartesiansOrbitals
                   k=k+1

                   do point = 1 , gridSize
                      Grid_instance(speciesID)%orbitalsWithGradient(k)%values(point,1) = auxMatrix(1)%values(point,j)
                      Grid_instance(speciesID)%orbitalsWithGradient(k)%values(point,2) = auxMatrix(2)%values(point,j)
                      Grid_instance(speciesID)%orbitalsWithGradient(k)%values(point,3) = auxMatrix(3)%values(point,j)
                      Grid_instance(speciesID)%orbitalsWithGradient(k)%values(point,4) = auxMatrix(4)%values(point,j)
                   end do

                   ! print *, "viveee"
                   ! call Matrix_show(orbitalAndGradientInGrid)

                   if(trim(action) .eq. "WRITE") then
                      write( labels(1), "(A,I0.4)") "ORBITAL_", k
                      labels(2) = Grid_instance(speciesID)%nameOfSpecies                   
                      call Matrix_writeToFile(Grid_instance(speciesID)%orbitalsWithGradient(k), unit=orbsUnit, binary=.true., arguments = labels(1:2) )
                   end if
                end do
             end do
          end do

          call Matrix_destructor(auxMatrix(1))
          call Matrix_destructor(auxMatrix(2))
          call Matrix_destructor(auxMatrix(3))
          call Matrix_destructor(auxMatrix(4))

          if(trim(action) .eq. "WRITE") close(unit=orbsUnit)

       end if

    end do

  end subroutine GridManager_atomicOrbitals

  !>
  !! @brief Returns the values of the density in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getDensityGradientAtGrid( speciesID, densityMatrix, densityInGrid, gradientInGrid)
    implicit none
    integer :: speciesID
    type(Matrix) :: densityMatrix
    type(Vector) :: densityInGrid
    type(Vector) :: gradientInGrid(*)

    integer :: gridSize
    integer :: numberOfCartesiansOrbitalsU
    integer :: numberOfCartesiansOrbitalsV
    integer :: point
    integer :: i, j, u, g
    integer :: ii, jj, v, gg
    integer :: s, ss
    integer :: numberOfContractions
    real(8) :: auxreal

    integer :: n, nproc
    real :: time1,time2
    type(Vector),allocatable :: nodeDensityInGrid(:)
    type(Vector),allocatable :: nodeGradientInGrid(:,:)

    
    gridSize = Grid_instance(speciesID)%totalSize
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
    nproc=omp_get_max_threads()

    if(.not. allocated(nodeDensityInGrid) ) allocate( nodeDensityInGrid(nproc),nodeGradientInGrid(nproc,3) )
    if(.not. allocated(nodeGradientInGrid) ) allocate( nodeGradientInGrid(nproc,3) )

    do n=1, nproc
       call Vector_Constructor( nodeDensityInGrid(n), gridSize, 0.0_8)
       call Vector_Constructor( nodeGradientInGrid(n,1), gridSize, 0.0_8)
       call Vector_Constructor( nodeGradientInGrid(n,2), gridSize, 0.0_8)
       call Vector_Constructor( nodeGradientInGrid(n,3), gridSize, 0.0_8)
    end do


    time1=omp_get_wtime()
    !$omp parallel private(n,u,v,point)
    n = omp_get_thread_num() +1
    !$omp do schedule (dynamic)
    do u = 1, numberOfContractions
       do point = 1 , gridSize

          nodeDensityInGrid(n)%values(point)=nodeDensityInGrid(n)%values(point)+densityMatrix%values(u,u)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)**2
          nodeGradientInGrid(n,1)%values(point)=nodeGradientInGrid(n,1)%values(point)+2*densityMatrix%values(u,u)&
               *Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,2)
          nodeGradientInGrid(n,2)%values(point)=nodeGradientInGrid(n,2)%values(point)+2*densityMatrix%values(u,u)&
               *Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,3)
          nodeGradientInGrid(n,3)%values(point)=nodeGradientInGrid(n,3)%values(point)+2*densityMatrix%values(u,u)&
               *Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,4)

       end do

       do v = u+1, numberOfContractions
          do point = 1 , gridSize

             nodeDensityInGrid(n)%values(point)=nodeDensityInGrid(n)%values(point)+2*densityMatrix%values(u,v)*&
                  Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1)

             nodeGradientInGrid(n,1)%values(point)=nodeGradientInGrid(n,1)%values(point)+2*densityMatrix%values(u,v)*&
                  (Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,2)+&
                  Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,2)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))

             nodeGradientInGrid(n,2)%values(point)=nodeGradientInGrid(n,2)%values(point)+2*densityMatrix%values(u,v)*&
                  (Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,3)+&
                  Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,3)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))

             nodeGradientInGrid(n,3)%values(point)=nodeGradientInGrid(n,3)%values(point)+2*densityMatrix%values(u,v)*&
                  (Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,4)+&
                  Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,4)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))

          end do
       end do
    end do
    !$omp end do 
    !$omp end parallel

    do n=1, nproc
       densityInGrid%values=densityInGrid%values+nodeDensityInGrid(n)%values
       gradientInGrid(1)%values=gradientInGrid(1)%values+nodeGradientInGrid(n,1)%values
       gradientInGrid(2)%values=gradientInGrid(2)%values+nodeGradientInGrid(n,2)%values
       gradientInGrid(3)%values=gradientInGrid(3)%values+nodeGradientInGrid(n,3)%values
       call Vector_Destructor( nodeDensityInGrid(n))
       call Vector_Destructor( nodeGradientInGrid(n,1))
       call Vector_Destructor( nodeGradientInGrid(n,2))
       call Vector_Destructor( nodeGradientInGrid(n,3))
    end do

    !!Check for negative values, which should be numerical mistakes
    auxreal=1.0E8
    do point = 1 , gridSize
       if(densityInGrid%values(point).lt.auxreal) auxreal=densityInGrid%values(point)
       if(densityInGrid%values(point).lt.0.0) densityInGrid%values(point)=0.0
    end do
    if(-auxreal.gt.CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD) print *, "found significative negative density, up to", auxreal, ", for species", speciesID
    
    
    time2=omp_get_wtime()
    ! write(*,"(A,F10.3,A4)") "**getDensityGradientAtGrid:", time2-time1 ," (s)"
    
    ! call Vector_show(gradientInGrid(1))
    ! call Vector_show(gradientInGrid(2))
    ! call Vector_show(gradientInGrid(3))

  end subroutine GridManager_getDensityGradientAtGrid

  !>
  !! @brief Returns the values of the exchange correlation potential for a specie in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getElectronicEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy, &
    otherSpeciesID, otherExchangeCorrelationEnergy)
    implicit none
    integer :: speciesID
    real(8) :: exchangeCorrelationEnergy
    integer, optional :: otherSpeciesID
    real(8), optional :: otherExchangeCorrelationEnergy

    character(50) :: nameOfSpecies, otherNameOfSpecies
    integer :: gridSize
    type(Vector) :: energyDensity
    type(Vector) :: sigma, sigmaPotential
    type(Vector) :: densityAB, potentialAB, sigmaAB, sigmaPotentialAB
    integer :: i, index, dir

    
    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )
    if( present(otherSpeciesID) )  otherNameOfSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )

    gridSize = Grid_instance(speciesID)%totalSize

    call Vector_Constructor(energyDensity, gridSize, 0.0_8)

    if (CONTROL_instance%CALL_LIBXC) then
       !Closed Shell
       if (nameOfSpecies=="E-" .and. .not. present(otherSpeciesID) ) then

          index=Functional_getIndex(speciesID)

          call Vector_Constructor(sigma, gridSize, 0.0_8)         
          call Vector_Constructor(sigmaPotential, gridSize, 0.0_8)
          !$omp parallel private(i)
          !$omp do schedule (dynamic)
          do i=1, gridSize

             if(Grid_instance(speciesID)%density%values(i) .gt. CONTROL_instance%ELECTRON_DENSITY_THRESHOLD) then !
                !libxc works with the gradient squared - sigma
                sigma%values(i)=(Grid_instance(speciesID)%densityGradient(1)%values(i)**2&
                     +Grid_instance(speciesID)%densityGradient(2)%values(i)**2&
                     +Grid_instance(speciesID)%densityGradient(3)%values(i)**2)

                !evaluates energy density, potential and sigma potential
                call Functional_libxcEvaluate(Functionals(index), 1, Grid_instance(speciesID)%density%values(i), &
                     sigma%values(i), energyDensity%values(i) , &
                     Grid_instance(speciesID)%potential%values(i), sigmaPotential%values(i) )
             end if
            
          end do
          !$omp end do 
          !$omp end parallel
          
          do i=1, gridSize
             !energy integral
             exchangeCorrelationEnergy=exchangeCorrelationEnergy&
                  +energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4)

             !convert to gradient potential
             do dir=1,3
                Grid_instance(speciesID)%gradientPotential(dir)%values(i)=Grid_instance(speciesID)%gradientPotential(dir)%values(i)&
                     +2.0*sigmaPotential%values(i)*Grid_instance(speciesID)%densityGradient(dir)%values(i)
             end do

          end do

          call Vector_Destructor(sigma)

          call Vector_Destructor(sigmaPotential)

          ! print *, "electronicEXC RKS", exchangeCorrelationEnergy
       elseif (nameOfSpecies=="E-ALPHA" .and. otherNameOfSpecies=="E-BETA") then

          call Vector_Constructor(densityAB, 2*gridSize, 0.0_8)         
          call Vector_Constructor(sigmaAB, 3*gridSize, 0.0_8)         
          call Vector_Constructor(potentialAB, 2*gridSize, 0.0_8)         
          call Vector_Constructor(sigmaPotentialAB, 3*gridSize, 0.0_8)         

          index=Functional_getIndex(speciesID)

          !$omp parallel private(i)
          !$omp do schedule (dynamic)
          do i=1, gridSize

             if(Grid_instance(speciesID)%density%values(i)+Grid_instance(otherSpeciesID)%density%values(i) .gt. CONTROL_instance%ELECTRON_DENSITY_THRESHOLD) then
                !libxc expects a single array of alpha and beta densities 
                densityAB%values(2*i-1)=Grid_instance(speciesID)%density%values(i)

                densityAB%values(2*i)=Grid_instance(otherSpeciesID)%density%values(i)

                !libxc works with the gradient squared - sigma
                sigmaAB%values(3*i-2)=(Grid_instance(speciesID)%densityGradient(1)%values(i)**2 + Grid_instance(speciesID)%densityGradient(2)%values(i)**2 + &
                     Grid_instance(speciesID)%densityGradient(3)%values(i)**2)

                sigmaAB%values(3*i-1)=Grid_instance(speciesID)%densityGradient(1)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(1)%values(i)+&
                     Grid_instance(speciesID)%densityGradient(2)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(2)%values(i)+&
                     Grid_instance(speciesID)%densityGradient(3)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(3)%values(i)

                sigmaAB%values(3*i)=(Grid_instance(otherSpeciesID)%densityGradient(1)%values(i)**2 + Grid_instance(otherSpeciesID)%densityGradient(2)%values(i)**2 + &
                     Grid_instance(otherSpeciesID)%densityGradient(3)%values(i)**2)


                !evaluates energy density, potential and sigma potential
                call Functional_libxcEvaluate(Functionals(index), 1, densityAB%values(2*i-1:2*i), sigmaAB%values(3*i-2:3*i), &
                     energyDensity%values(i) , potentialAB%values(2*i-1:2*i), sigmaPotentialAB%values(3*i-2:3*i) )

                !potential assignment
                Grid_instance(speciesID)%potential%values(i)=Grid_instance(speciesID)%potential%values(i)+potentialAB%values(2*i-1)

                Grid_instance(otherSpeciesID)%potential%values(i)=Grid_instance(otherSpeciesID)%potential%values(i)+potentialAB%values(2*i)

             end if
          end do
          !$omp end do 
          !$omp end parallel

          do i=1, gridSize
             !energy integrals
             exchangeCorrelationEnergy=exchangeCorrelationEnergy+&
                  energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4)

             otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+&
                  energyDensity%values(i)*Grid_instance(otherSpeciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 


             !convert to gradient potential
             do dir=1,3
                Grid_instance(speciesID)%gradientPotential(dir)%values(i)=Grid_instance(speciesID)%gradientPotential(dir)%values(i)&
                     +2.0*sigmaPotentialAB%values(3*i-2)*Grid_instance(speciesID)%densityGradient(dir)%values(i)&
                     +Grid_instance(otherSpeciesID)%densityGradient(dir)%values(i)*sigmaPotentialAB%values(3*i-1)
                
                Grid_instance(otherSpeciesID)%gradientPotential(dir)%values(i)=Grid_instance(otherSpeciesID)%gradientPotential(dir)%values(i)&
                     +2.0*sigmaPotentialAB%values(3*i)*Grid_instance(otherSpeciesID)%densityGradient(dir)%values(i)&
                     +Grid_instance(speciesID)%densityGradient(dir)%values(i)*sigmaPotentialAB%values(3*i-1)

             end do

             ! print *, i, densityAB%values(2*i-1), densityAB%values(2*i), energyDensity%values(i), potentialAB%values(2*i-1), potentialAB%values(2*i)
             ! print *, sigmaAB%values(3*i-2), sigmaAB%values(3*i-1), sigmaAB%values(3*i), sigmaPotentialAB%values(3*i-2), sigmaPotentialAB%values(3*i-1), sigmaPotentialAB%values(3*i)
             
          end do
          
          ! print *, "density", densityAB%values(1), densityAB%values(2*gridSize)
          ! call Vector_show(densityAB)
          call Vector_Destructor(densityAB)

          ! print *, "sigma"
          ! call Vector_show(sigmaAB)
          call Vector_Destructor(sigmaAB)

          ! print *, "electronicEXC UKS", exchangeCorrelationEnergy, otherExchangeCorrelationEnergy
          call Vector_Destructor(potentialAB)         
          call Vector_Destructor(sigmaPotentialAB)         

       end if
    else

       index=Functional_getIndex(speciesID)
       if ( Functionals(index)%name .eq. "exchange:Slater-correlation:VWN5") then

          if (nameOfSpecies=="E-" .and. .not. present(otherSpeciesID) ) then

             call Functional_LDAEvaluate(gridSize, Grid_instance(speciesID)%density%values/2, Grid_instance(speciesID)%density%values/2, &
                  energyDensity%values, Grid_instance(speciesID)%potential%values )
             do i=1, gridSize
                exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%points%values(i,4) 
             end do

             ! print *, "electronicEXC RKS", exchangeCorrelationEnergy
          elseif (nameOfSpecies=="E-ALPHA" .and. otherNameOfSpecies=="E-BETA") then

             call Functional_LDAEvaluate(gridSize, Grid_instance(speciesID)%density%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, Grid_instance(speciesID)%potential%values, Grid_instance(otherSpeciesID)%potential%values )

             do i=1, gridSize
                exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
                otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(otherSpeciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
             end do

             ! print *, "electronicEXC UKS", exchangeCorrelationEnergy, otherExchangeCorrelationEnergy
          end if

       end if
    end if
 
    call Vector_Destructor(energyDensity)

    ! do i=1, gridSize
    !    print *, densityInGrid%values(i), exchange%values(i), correlationA%values(i)
    ! end do

  end subroutine GridManager_getElectronicEnergyAndPotentialAtGrid

  !>
  !! @brief Returns the values of the exchange correlation potential for a specie in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getInterspeciesEnergyAndPotentialAtGrid( speciesID, otherSpeciesID, exchangeCorrelationEnergy, &
       otherElectronID, otherElectronExchangeCorrelationEnergy) 

    implicit none
    integer :: speciesID
    integer :: otherSpeciesID
    real(8) :: exchangeCorrelationEnergy
    integer, optional :: otherElectronID
    real(8), optional :: otherElectronExchangeCorrelationEnergy

    character(50) :: nameOfSpecies, otherNameOfSpecies, auxstring
    integer :: gridSize, otherGridSize
    type(Vector) :: energyDensity
    type(Vector) :: sigma
    type(Vector) :: densityAB, potentialAB, sigmaAB, sigmaPotentialAB
    type(Vector) :: electronicDensityAtOtherGrid, electronicGradientAtOtherGrid(3), electronicPotentialAtOtherGrid, electronicGradientPotentialAtOtherGrid(3)
    integer :: i, j, k, index, dir

    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )
    otherNameOfSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )

    if ( (nameOfSpecies.ne."E-" .and. nameOfSpecies.ne."E-ALPHA") ) return

    gridSize = Grid_instance(speciesID)%totalSize
    otherGridSize = Grid_instance(otherSpeciesID)%totalSize

    index=Functional_getIndex(speciesID, otherSpeciesID)

    !Nuclear electron correlation
    !"E-BETA" is treated at the same time as E-ALPHA

    call Vector_constructor(energyDensity, otherGridSize, 0.0_8)
    call Vector_constructor(electronicDensityAtOtherGrid, otherGridSize, 0.0_8)
    call Vector_constructor(electronicPotentialAtOtherGrid, otherGridSize, 0.0_8)
    do dir=1,3
       call Vector_constructor(electronicGradientAtOtherGrid(dir), otherGridSize, 0.0_8)
       call Vector_constructor(electronicGradientPotentialAtOtherGrid(dir), otherGridSize, 0.0_8)
    end do

    !!This adds E-BETA density and gradient
    call GridManager_getElectronicDensityInOtherGrid(speciesID, otherSpeciesID, &
         GridsCommonPoints(speciesID,otherSpeciesID)%totalSize, int(GridsCommonPoints(speciesID,otherSpeciesID)%points%values), &
         electronicDensityAtOtherGrid, electronicGradientAtOtherGrid)

    if (trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL) .ne. "NONE") then
       auxstring=trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
    else if(trim(CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL) .ne. "NONE") then
       auxstring=trim(CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL)
    else
       auxstring="NONE"
    end if
!!!These routines return the electronic energy density
    select case (trim(auxstring) )

    case ("EXPCS-GGA")
       call Functional_expCSGGAEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid, electronicGradientAtOtherGrid, &
            Grid_instance(otherSpeciesID)%density, Grid_instance(otherSpeciesID)%densityGradient, &
            energyDensity, electronicPotentialAtOtherGrid, electronicGradientPotentialAtOtherGrid, &
            Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%gradientPotential )

    case ("EXPCS-GGA-NOA")
       call Functional_expCSGGAEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid, electronicGradientAtOtherGrid, &
            Grid_instance(otherSpeciesID)%density, Grid_instance(otherSpeciesID)%densityGradient, &
            energyDensity, electronicPotentialAtOtherGrid, electronicGradientPotentialAtOtherGrid, &
            Grid_instance(otherSpeciesID)%potential, Grid_instance(otherSpeciesID)%gradientPotential )

    case ("EPC17-1")
       call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("EPC17-2")
       call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("IKN-NSF")
       call Functional_IKNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("MLCS-FIT")
       call Functional_MLCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("MLCS-A")
       call Functional_MLCSAEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("MLCS-AN")
       call Functional_MLCSANEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("CS-MYFIT")
       call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("IMAMURA-MYFIT")
       call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("MEJIA-MYFIT")
       call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("MEJIAA-MYFIT")
       call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("EXPCS-A")
       call Functional_expCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("EXPCS-NOA")
       call Functional_expCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("PSN")
       call Functional_PSNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("PSNAP")
       call Functional_PSNAPEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
            electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
            energyDensity%values, electronicPotentialAtOtherGrid%values, Grid_instance(otherSpeciesID)%potential%values  )

    case ("NONE")

    case default
       print *, trim(auxstring)
       STOP "The "//otherNameOfSpecies//"electron functional chosen is not implemented"

    end select

    !!Adds the nuclear electron potential to the relevant points in the electronic grid
!!!We partition the energy considering the electronic density
    do k=1, GridsCommonPoints(speciesID,otherSpeciesID)%totalSize
       !electron index
       i=int(GridsCommonPoints(speciesID,otherSpeciesID)%points%values(k,1))
       !nuclear index
       j=int(GridsCommonPoints(speciesID,otherSpeciesID)%points%values(k,2)) 

       exchangeCorrelationEnergy=exchangeCorrelationEnergy+&
            energyDensity%values(j)*Grid_instance(speciesID)%density%values(i)*Grid_instance(otherSpeciesID)%points%values(j,4)

       Grid_instance(speciesID)%potential%values(i) = Grid_instance(speciesID)%potential%values(i) + electronicPotentialAtOtherGrid%values(j)

       do dir=1,3
          Grid_instance(speciesID)%gradientPotential(dir)%values(i) = Grid_instance(speciesID)%gradientPotential(dir)%values(i) + electronicGradientPotentialAtOtherGrid(dir)%values(j)
       end do

       if(nameOfSpecies .eq. "E-ALPHA") then

          otherElectronExchangeCorrelationEnergy=otherElectronExchangeCorrelationEnergy+&
               energyDensity%values(j)*Grid_instance(otherElectronID)%density%values(i)*Grid_instance(otherSpeciesID)%points%values(j,4)

          Grid_instance(otherElectronID)%potential%values(i) = Grid_instance(otherElectronID)%potential%values(i) + electronicPotentialAtOtherGrid%values(j)
          do dir=1,3
             Grid_instance(otherElectronID)%gradientPotential(dir)%values(i) = Grid_instance(otherElectronID)%gradientPotential(dir)%values(i)&
                  + electronicGradientPotentialAtOtherGrid(dir)%values(j)
          end do

       end if

    end do

    call Vector_Destructor(energyDensity)
    call Vector_destructor(electronicDensityAtOtherGrid)
    call Vector_destructor(electronicPotentialAtOtherGrid)
    do dir=1,3
       call Vector_destructor(electronicGradientAtOtherGrid(dir))
       call Vector_destructor(electronicGradientPotentialAtOtherGrid(dir))
    end do
    ! do i=1, gridSize
    !    print *, densityInGrid%values(i), exchange%values(i), correlationA%values(i)
    ! end do

    ! call Vector_Destructor(correlationA)
    ! call Vector_Destructor(correlationB)


  end subroutine GridManager_getInterspeciesEnergyAndPotentialAtGrid

  !>
  !! @brief Builds the exchange correlation for a species
  ! Felix Moncada, 2017
  subroutine GridManager_buildExchangeCorrelationMatrix( speciesID, exchangeCorrelationMatrix)
    implicit none
    integer :: speciesID
    type(Matrix) :: exchangeCorrelationMatrix

    integer :: gridSize
    integer :: numberOfContractions

    integer :: numberOfCartesiansOrbitalsU
    integer :: numberOfCartesiansOrbitalsV
    integer :: u, v, point

    real(8) :: time1, time2
    integer :: n, nproc
    type(Matrix),allocatable :: nodeExchangeCorrelationMatrix(:)

    gridSize = Grid_instance(speciesID)%totalSize
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    nproc=omp_get_max_threads()

    if(.not. allocated(nodeExchangeCorrelationMatrix))allocate( nodeExchangeCorrelationMatrix(nproc) )

    do n=1, nproc
       call Matrix_Constructor( nodeExchangeCorrelationMatrix(n), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
    end do

    ! time1=omp_get_wtime()
    !$omp parallel private(n,u,v,point)
    n = omp_get_thread_num() +1
    !$omp do schedule (dynamic)
    do u = 1, numberOfContractions

       do point = 1 , gridSize
          nodeExchangeCorrelationMatrix(n)%values(u,u)=&
               nodeExchangeCorrelationMatrix(n)%values(u,u)&
               +(Grid_instance(speciesID)%potential%values(point)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)**2&
               +2.0*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)&
               *(Grid_instance(speciesID)%gradientPotential(1)%values(point)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,2)&
               +Grid_instance(speciesID)%gradientPotential(2)%values(point)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,3)&
               +Grid_instance(speciesID)%gradientPotential(3)%values(point)*Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,4))&
               )*Grid_instance(speciesID)%points%values(point,4)
       end do

       do v = u+1, numberOfContractions

          do point = 1 , gridSize
             nodeExchangeCorrelationMatrix(n)%values(u,v)=&
                  nodeExchangeCorrelationMatrix(n)%values(u,v)&
                  +(Grid_instance(speciesID)%potential%values(point)&
                  *Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1)&
                  +Grid_instance(speciesID)%gradientPotential(1)%values(point)&
                  *(Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,2)&
                  +Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,2)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))&
                  +Grid_instance(speciesID)%gradientPotential(2)%values(point)&
                  *(Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,3)&
                  +Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,3)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))&
                  +Grid_instance(speciesID)%gradientPotential(3)%values(point)&
                  *(Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,1)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,4)&
                  +Grid_instance(speciesID)%orbitalsWithGradient(u)%values(point,4)*Grid_instance(speciesID)%orbitalsWithGradient(v)%values(point,1))&
                  )*Grid_instance(speciesID)%points%values(point,4)

          end do


          ! call Stopwatch_stop(lowdin_stopwatch)     
          ! time2=time2+lowdin_stopwatch%enlapsetTime

       end do
       
    end do
    !$omp end do 
    !$omp end parallel

    do n=1, nproc
       exchangeCorrelationMatrix%values=exchangeCorrelationMatrix%values+nodeExchangeCorrelationMatrix(n)%values
       call Matrix_Destructor( nodeExchangeCorrelationMatrix (n))
    end do
    
    do u=1, numberOfContractions
       do v=u+1, numberOfContractions
          exchangeCorrelationMatrix%values(v,u)=exchangeCorrelationMatrix%values(u,v)
       end do
    end do

    ! time2=omp_get_wtime()
    ! write(*,"(A,F10.3,A4)") "**buildExchangeCorrelationMatrix:", time2-time1 ," (s)"

    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** reading orbital files:", time1 ," (s)"
    ! write(*,"(A,F10.3,A4)") "** integrating over the grid:", time2 ," (s)"

  end subroutine GridManager_buildExchangeCorrelationMatrix



  subroutine GridManager_getElectronicDensityInOtherGrid(electronicID,otherSpeciesID, commonGridSize, commonPoints, electronicDensityAtOtherGrid, electronicGradientAtOtherGrid )
    integer :: electronicID, otherSpeciesID
    integer :: commonGridSize
    integer :: commonPoints(commonGridSize,2)
    type(Vector) :: electronicDensityAtOtherGrid
    type(Vector) :: electronicGradientAtOtherGrid(3)
    
    character(50) :: nameOfElectron
    integer :: otherElectronicID
    integer :: electronicGridSize, otherGridSize
    integer :: i,j,k
    real :: time1,time2

    electronicGridSize=Grid_instance(electronicID)%totalSize
    otherGridSize=Grid_instance(otherSpeciesID)%totalSize

    nameOfElectron=MolecularSystem_getNameOfSpecie(electronicID)
    if (nameOfElectron .eq. "E-ALPHA") otherElectronicID=MolecularSystem_getSpecieID( "E-BETA" )
    if (nameOfElectron .eq. "E-BETA") otherElectronicID=MolecularSystem_getSpecieID( "E-ALPHA" )
       
    
    call Vector_constructor(electronicDensityAtOtherGrid, otherGridSize, 1.0E-12_8)
    
    time1=omp_get_wtime()
    do k=1, commonGridSize
       !here we are assuming that the electron came in the first position
       i=commonPoints(k,1)
       j=commonPoints(k,2)
       if(nameOfElectron .eq. "E-ALPHA" .or. nameOfElectron .eq. "E-BETA") then
          electronicDensityAtOtherGrid%values(j)=Grid_instance(electronicID)%density%values(i)+Grid_instance(otherElectronicID)%density%values(i)
          electronicGradientAtOtherGrid(1)%values(j)=Grid_instance(electronicID)%densityGradient(1)%values(i)+Grid_instance(otherElectronicID)%densityGradient(1)%values(i)
          electronicGradientAtOtherGrid(2)%values(j)=Grid_instance(electronicID)%densityGradient(2)%values(i)+Grid_instance(otherElectronicID)%densityGradient(2)%values(i)
          electronicGradientAtOtherGrid(3)%values(j)=Grid_instance(electronicID)%densityGradient(3)%values(i)+Grid_instance(otherElectronicID)%densityGradient(3)%values(i)
       else
          electronicDensityAtOtherGrid%values(j)=Grid_instance(electronicID)%density%values(i)
          electronicGradientAtOtherGrid(1)%values(j)=Grid_instance(electronicID)%densityGradient(1)%values(i)
          electronicGradientAtOtherGrid(2)%values(j)=Grid_instance(electronicID)%densityGradient(2)%values(i)
          electronicGradientAtOtherGrid(3)%values(j)=Grid_instance(electronicID)%densityGradient(3)%values(i)
       end if
    end do
    time2=omp_get_wtime()
    ! write(*,"(A,F10.3,A4)") "**getElectronicDensityInOtherGrid:", time2-time1 ," (s)"

  end subroutine GridManager_getElectronicDensityInOtherGrid

  subroutine GridManager_findCommonPoints(grid,gridSize,otherGrid,otherGridSize,commonPoints,commonSize)
    type(Matrix) :: grid,otherGrid
    integer :: gridSize,otherGridSize,commonSize
    type(Matrix) :: commonPoints

    type(Matrix) :: auxFinder
    integer :: i,j, k,last, point
    real :: time1,time2

    !The two grids must have been generated from the same amount of atomic angular and radial points
    call Matrix_constructor(auxFinder, int(max(gridSize,otherGridSize),8), int(2,8), 0.0_8)

    ! print *, "matrix 1 size is", size(grid%values,1), gridSize
    ! print *, "matrix 2 size is", size(otherGrid%values,1), otherGridSize
    ! print *, "auxFinder size is", size(auxFinder%values,1), max(gridSize,otherGridSize)
    
    time1=omp_get_wtime()
    point=0
    k=0
    do i=1, gridSize
       do j=1, otherGridSize
          k=k+1
          if(  (grid%values(i,1) .eq. otherGrid%values(k,1)) .and. (grid%values(i,2) .eq. otherGrid%values(k,2)) .and. (grid%values(i,3) .eq. otherGrid%values(k,3)) ) then
             point=point+1
             auxFinder%values(point,1)=i
             auxFinder%values(point,2)=k
             ! print *, i, j, k
             ! print *, grid%values(i,1:3), otherGrid%values(k,1:3)
             exit
          end if
          if(k.ge.otherGridSize) k=0
       end do
       if(point .ge. max(gridSize,otherGridSize)) exit
    end do

    call Matrix_constructor(commonPoints, int(point,8), int(2,8), 0.0_8)

    do i=1, point
       commonPoints%values(i,1)=auxFinder%values(i,1)
       commonPoints%values(i,2)=auxFinder%values(i,2)
    end do
    commonSize=point
    
    ! print *, "number of common points", point

    time2=omp_get_wtime()
    ! write(*,"(A,F10.3,A4)") "**FindCommonPoints:", time2-time1 ," (s)"
    ! Call Matrix_show(commonPoints)

  end subroutine GridManager_FindCommonPoints

  !>
  !! @brief Builds the exchange correlation for a species
  ! Felix Moncada, 2017
  subroutine GridManager_getContactDensity( speciesID, otherSpeciesID, otherElectronID)
    implicit none
    integer :: speciesID, otherSpeciesID
    integer, optional :: otherElectronID
    
    integer :: gridSize
    integer :: point
    real(8) :: kf,a0n,a1n,a2n,a3n,a4n,a0d,a1d,a2d,p,q0,q2,q4,Eab,Eab2
    real(8) :: rhoE,rhoP, rhoTot, rhoDif, npos, densityThreshold
    real(8) :: beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP

    integer :: n, nproc
    real(8) :: contactDensity, overlapDensity
    character(50) :: auxstring
    type(Vector) :: electronicDensityAtOtherGrid, electronicGradientAtOtherGrid(3), gfactor

    gridSize =Grid_instance(otherSpeciesID)%totalSize

    call Vector_constructor(electronicDensityAtOtherGrid, gridSize, 0.0_8)
    call Vector_constructor(electronicGradientAtOtherGrid(1), gridSize, 0.0_8)
    call Vector_constructor(electronicGradientAtOtherGrid(2), gridSize, 0.0_8)
    call Vector_constructor(electronicGradientAtOtherGrid(3), gridSize, 0.0_8)

    !electrons go on the first position
    call GridManager_getElectronicDensityInOtherGrid(speciesID, otherSpeciesID, &
         GridsCommonPoints(speciesID,otherSpeciesID)%totalSize, int(GridsCommonPoints(speciesID,otherSpeciesID)%points%values), electronicDensityAtOtherGrid, electronicGradientAtOtherGrid )

    call Vector_constructor(gfactor, gridSize, 0.0_8)
    gfactor%values=0.0_8
    
    nproc=omp_get_max_threads()

    contactDensity=0.0
    overlapDensity=0.0
    ! do n=1, nproc
    !    call Matrix_Constructor( nodeExchangeCorrelationMatrix(n), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
    ! end do

    if (trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL) .ne. "NONE") then
       auxstring=trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
    else if(trim(CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL) .ne. "NONE") then
       auxstring=trim(CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL)
    else
       auxstring="NONE"
    end if

    
    if(MolecularSystem_getMass(otherSpeciesID) .lt. 2.0 .and. (auxstring.eq."expCS-A" .or. auxstring.eq."expCS-GGA")) then !positron

       kf=2.2919886876120283056
       a0n=0.3647813291441602
       a1n=0.04801434878972582
       a2n=1.6987053215381047
       a3n=0.428189835287642
       a4n=1.0
       a0d=1.0
       a1d=0.011817504796076739
       a2d=0.8862269254527579

       densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD

       !$omp parallel private(n,beta, rhoE,rhoP, rhoTot, rhoDif, point), shared(gfactor)
       n = omp_get_thread_num() +1
       !$omp do schedule (dynamic)
       do point = 1 , gridSize
          rhoE=electronicDensityAtOtherGrid%values(point)
          rhoP=Grid_instance(otherSpeciesID)%density%values(point)

          call Functional_getBeta( rhoE, rhoP, MolecularSystem_getMass(otherSpeciesID), kf, beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP)

          if( rhoE .gt. densityThreshold .and. rhoP .gt. densityThreshold ) then !
             gfactor%values(point)=(a0n+a1n*beta+a2n*beta**2+a3n*beta**3+a4n*beta**4)/(a0d*beta**3+a1d*beta**4+a2d*beta**5)
          else if( (rhoE .gt. densityThreshold) .or. (rhoP .gt. densityThreshold)  ) then 
             gfactor%values(point)=a0n/(a0d*beta**3)
          else
             gfactor%values(point)=0.0
          end if
          ! print *, point, rhoE, rhoP, beta, gfactor%values(point)
       end do
       !$omp end do 
       !$omp end parallel
    end if

    do point = 1 , gridSize

       overlapDensity=overlapDensity+&
            electronicDensityAtOtherGrid%values(point)*&
            Grid_instance(otherSpeciesID)%density%values(point)*&
            Grid_instance(otherSpeciesID)%points%values(point,4)

       contactDensity=contactDensity+&
            electronicDensityAtOtherGrid%values(point)*&
            Grid_instance(otherSpeciesID)%density%values(point)*&
            (1.0+gfactor%values(point))*&
            Grid_instance(otherSpeciesID)%points%values(point,4) 
    end do

    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       print *, ""
       print *, "Contact density between ", trim(MolecularSystem_getNameOfSpecie(speciesID)),"-", trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))
       if(present(otherElectronID) ) print *, "Including contact density between ", trim(MolecularSystem_getNameOfSpecie(otherElectronID)),"-", trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))
       
       if(auxstring.eq."expCS-A" .or. auxstring.eq."expCS-GGA" ) then
          print *, "As the integral of rhoA*rhoB(1+g[beta])"
          print *, "With g[beta] from the expCS-A functional"
       else
          print *, "As the integral of rhoA*rhoB"
       end if
       write (*,"(A10,F20.10)") "overlap=", overlapDensity
       write (*,"(A10,F20.10)") "pep=", contactDensity
    end if
    ! npos=0.0
    ! do point = 1 , gridSize

    !    npos=npos+&
    !         Grid_instance(otherSpeciesID)%density%values(point)*&
    !         Grid_instance(otherSpeciesID)%points%values(point,4) 
    ! end do

    ! print *, "npos=", npos
    
  end subroutine GridManager_getContactDensity

    !>
  !! @brief Builds the exchange correlation for a species
  ! Felix Moncada, 2017
  subroutine GridManager_getExpectedDistances( speciesID)
    implicit none
    integer :: speciesID
    
    integer :: point,gridSize
    integer :: center,numberOfCenters
    real(8), allocatable :: distances(:)


    gridSize =Grid_instance(speciesID)%totalSize
    numberOfCenters=MolecularSystem_instance%numberOfPointCharges

    if(.not. allocated(distances))allocate(distances(numberOfCenters))

    distances(:)=0.0
    do center = 1, numberOfCenters
       
       do point = 1 , gridSize

          distances(center)=distances(center) + &
               Grid_instance(speciesID)%density%values(point)* &
               sqrt((MolecularSystem_instance%pointCharges(center)%origin(1)-Grid_instance(speciesID)%points%values(point,1) )**2+ &
               (MolecularSystem_instance%pointCharges(center)%origin(2)-Grid_instance(speciesID)%points%values(point,2) )**2+ &
               (MolecularSystem_instance%pointCharges(center)%origin(3)-Grid_instance(speciesID)%points%values(point,3) )**2 ) * &
               Grid_instance(speciesID)%points%values(point,4) 

       end do
       
       write (*,"(A10,A10,F20.10)") trim(Grid_instance(speciesID)%nameOfSpecies), trim(MolecularSystem_instance%pointCharges(center)%nickname), distances(center)
       
    end do


    
  end subroutine GridManager_getExpectedDistances

  
end module GridManager_

  !>
  !! @brief Returns the values of the gradient of a contracted atomic shell in a set of coordinates
!!! Felix Moncada, 2017
  !<
  ! subroutine GridManager_getOrbitalAtGrid( this, gridA, gridSize, output)
  !   implicit none
  !   type(ContractedGaussian) , intent(in) :: this
  !   type(Matrix) :: gridA
  !   type(Matrix) :: output
  !   integer :: gridSize
    
  !   integer :: h
  !   integer :: nx, ny, nz !< indices de momento angular
  !   integer :: i, j, m, xx
  !   integer :: point
  !   real(8) :: coordinate(3)
  !   real(8) :: exponential
  !   real(8) :: auxOutput(this%numCartesianOrbital)

  !   do point=1, gridSize
  !      coordinate(1)=gridA%values(point,1)
  !      coordinate(2)=gridA%values(point,2)
  !      coordinate(3)=gridA%values(point,3)
  !      do h=1, this%length
  !         exponential=dexp(-this%orbitalExponents(h) &
  !              *(  (this%origin(1)-coordinate(1))**2 &
  !              +(this%origin(2)-coordinate(2))**2 &
  !              +(this%origin(3)-coordinate(3))**2) )
  !         m = 0
  !         do i = 0 , this%angularMoment
  !            nx = this%angularMoment - i
  !            do j = 0 , i
  !               ny = i - j
  !               nz = j
  !               m = m + 1
  !               auxOutput(m) = this%contNormalization(m) &
  !                    * this%primNormalization(h,m) &
  !                    * (coordinate(1)-this%origin(1))** nx &
  !                    * (coordinate(2)-this%origin(2))** ny &
  !                    * (coordinate(3)-this%origin(3))** nz &
  !                    * exponential 
  !            end do
  !         end do

  !         do xx=1, m
  !            output%values(point,xx) = output%values(point,xx) + auxOutput(xx) * this%contractionCoefficients(h)
  !         end do
  !      end do
  !   end do
  ! end subroutine GridManager_getOrbitalAtGrid

!   >
!   @brief Returns the values of all atomic orbital isn a set of coordinates
! Felix Moncada, 2017
!   <
  ! subroutine GridManager_getOrbitalMatrix( speciesID, grid, gridSize, orbitalsInGrid)
  !   implicit none
  !   integer :: speciesID
  !   type(Matrix) :: grid, orbitalsInGrid
  !   integer :: gridSize

  !   type(Matrix) :: auxMatrix
  !   integer :: numberOfCartesiansOrbitals
  !   integer :: totalNumberOfContractions
  !   integer :: point
  !   integer :: i, j, k, g

  !   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

  !   k=1
  !   do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
  !      do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)
  !         numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital

  !         call Matrix_constructor( auxMatrix, int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8)

  !         call GridManager_getOrbitalAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), grid, gridSize, auxMatrix)
          
  !         do j = 1, numberOfCartesiansOrbitals
  !            do point = 1 , gridSize
  !               orbitalsInGrid%values(point,k) = auxMatrix%values(point,j)
  !            end do
  !            k=k+1
  !         end do
  !      end do
  !   end do
    
  ! end subroutine GridManager_getOrbitalMatrix

  !>
  !! @brief Returns the values of a contracted atomic shell in a set of coordinates
!!! Felix Moncada, 2017
  !<
  ! subroutine GridManager_getOrbitalGradientAtPoint( this, gridPoint, output)
  !   implicit none
  !   type(ContractedGaussian) , intent(in) :: this
  !   real(8) :: gridPoint(3)
  !   integer :: gridSize
  !   real(8) :: output(this%numCartesianOrbital,4)

  !   real(8) :: coordinate(3)
  !   integer :: h
  !   integer :: nx, ny, nz !< indices de momento angular
  !   integer :: i, j, m, w
  !   integer :: point
  !   real(8) :: exponential, dx, dy, dz

  !   coordinate(1)=gridPoint(1)-this%origin(1)
  !   coordinate(2)=gridPoint(2)-this%origin(2)
  !   coordinate(3)=gridPoint(3)-this%origin(3)
  !   do h=1, this%length
  !      exponential= this%contractionCoefficients(h) * dexp(-this%orbitalExponents(h)*(coordinate(1)**2 + coordinate(2)**2 +coordinate(3)**2) )
  !      m = 0
  !      do i = 0 , this%angularMoment
  !         nx = this%angularMoment - i
  !         do j = 0 , i
  !            ny = i - j
  !            nz = j
  !            m = m + 1

  !            !!Orbital
  !            output(m,4) = this%contNormalization(m) &
  !                 * this%primNormalization(h,m) &
  !                 * coordinate(1)** nx &
  !                 * coordinate(2)** ny &
  !                 * coordinate(3)** nz &
  !                 * exponential 

  !            dx=-2*this%orbitalExponents(h) &
  !                 * coordinate(1)** (nx+1) &
  !                 * coordinate(2)** ny &
  !                 * coordinate(3)** nz 

  !            ! Orbital derivative
  !            if( nx .ge. 1 ) then
  !               dx= dx + &
  !                    nx*coordinate(1)** (nx-1) &
  !                    * coordinate(2)** ny &
  !                    * coordinate(3)** nz 
  !            end if

  !            dy=-2*this%orbitalExponents(h) &
  !                 * coordinate(1)** nx &
  !                 * coordinate(2)** (ny+1) &
  !                 * coordinate(3)** nz 

  !            if( ny .ge. 1 ) then
  !               dy= dy + &
  !                    coordinate(1)** nx &
  !                    *ny*coordinate(2)** (ny-1) &
  !                    * coordinate(3)** nz 
  !            end if

  !            dz=-2*this%orbitalExponents(h) &
  !                 * coordinate(1)** nx &
  !                 * coordinate(2)** ny &
  !                 * coordinate(3)** (nz+1) 

  !            if( nz .ge. 1 ) then
  !               dz= dz+ &
  !                    coordinate(1)** nx &
  !                    *coordinate(2)** ny &
  !                    *nz*coordinate(3)** (nz-1) 
  !            end if

  !            output(m,1) = this%contNormalization(m) &
  !                 *this%primNormalization(h,m) &
  !                 *exponential*dx

  !            output(m,2) = this%contNormalization(m) &
  !                 *this%primNormalization(h,m) &
  !                 *exponential*dy

  !            output(m,3) = this%contNormalization(m) &
  !                 *this%primNormalization(h,m) &
  !                 *exponential*dz

  !         end do
  !      end do
  !   end do

  ! end subroutine GridManager_getOrbitalGradientAtPoint

  
