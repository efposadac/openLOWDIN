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
       GridManager_getOrbitalGradientAtGrid, &
       GridManager_getDensityGradientAtGrid, &
       GridManager_getEnergyAndPotentialAtGrid, & !, &       GridManager_getEnergyFromGrid
       GridManager_writeAtomicOrbitals, &
       GridManager_buildExchangeCorrelationMatrix!, &
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
    integer :: speciesID
    character(50) :: labels(2), dftFile
    integer :: dftUnit


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate memory.
    allocate(Grid_instance(numberOfSpecies))

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       call Grid_constructor(Grid_instance(speciesID), speciesID , type )

    end do


  end subroutine GridManager_buildGrids
  !>
  !! @brief Writes a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_writeGrids( )
    implicit none
    integer :: numberOfSpecies
    integer :: speciesID
    character(50) :: labels(2), dftFile
    integer :: dftUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       !! Open file for dft
       dftUnit = 77
       dftFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".grid"
       open(unit = dftUnit, file=trim(dftFile), status="replace", form="unformatted")

       labels(2) = Grid_instance(speciesID)%nameOfSpecies
       labels(1) = "GRID-SIZE"

       call Vector_writeToFile(unit=dftUnit, binary=.true., value=real(Grid_instance(speciesID)%totalSize,8), arguments= labels )

       !! This goes here for convenience only
       labels(1) = "EXACT-EXCHANGE-FRACTION" 

       call Vector_writeToFile(unit=dftUnit, binary=.true., value=Functional_getExchangeFraction(speciesID), arguments= labels )

       labels(1) = "INTEGRATION-GRID"
       call Matrix_writeToFile(Grid_instance(speciesID)%points, unit=dftUnit, binary=.true., arguments = labels(1:2) )

       ! call Matrix_show (Grid_instance(speciesID)%points)

       close(unit=dftUnit)
    end do

  end subroutine GridManager_writeGrids


  !>
  !! @brief Reads a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_readGrids( )
    implicit none
    integer :: numberOfSpecies
    integer :: speciesID
    character(50) :: labels(2), dftFile
    integer :: dftUnit
    real(8) :: auxVal

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate memory.
    allocate(Grid_instance(numberOfSpecies))

    !! Build and write species grids
    do speciesID = 1, numberOfSpecies

       Grid_instance(speciesID)%nameOfSpecies=trim(MolecularSystem_getNameOfSpecie(speciesID))
       !! Open file for dft
       dftUnit = 77
       dftFile = trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".grid"
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

  end subroutine GridManager_readGrids


  !>
  !! @brief Writes the values of all the atomic orbitals and their gradients in a set of coordinates to a file
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_writeAtomicOrbitals( )
    implicit none

    integer :: numberOfSpecies
    integer :: totalNumberOfContractions
    integer :: speciesID
    integer :: gridSize
    integer :: mu,nu, point, index
    type(Matrix) :: grid
    type(Matrix) :: orbitalAndGradientInGrid

    character(50) :: labels(2), dftFile
    integer :: dftUnit

    type(Matrix) :: auxMatrix(4)
    integer :: i, j, k, g
    integer :: numberOfCartesiansOrbitals

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do speciesID = 1, numberOfSpecies

       totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
       gridSize = Grid_instance(speciesID)%totalSize

       k=1
       do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
          do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)
             numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital

             ! print *, "holaaaaa"

             call Matrix_constructor( auxMatrix(1), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !orbital
             call Matrix_constructor( auxMatrix(2), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dx
             call Matrix_constructor( auxMatrix(3), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dy
             call Matrix_constructor( auxMatrix(4), int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8) !d orbital/dz

             call GridManager_getOrbitalGradientAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), &
                  Grid_instance(speciesID)%points, gridSize, auxMatrix(1), auxMatrix(2), auxMatrix(3), auxMatrix(4))

             ! call Matrix_show(auxMatrix(1))

             ! print *, "wololooooo"

             do j = 1, numberOfCartesiansOrbitals

                call Matrix_constructor(orbitalAndGradientInGrid, int(gridSize,8), int(4,8), 0.0_8)

                do point = 1 , gridSize
                   orbitalAndGradientInGrid%values(point,1) = auxMatrix(1)%values(point,j)
                   orbitalAndGradientInGrid%values(point,2) = auxMatrix(2)%values(point,j)
                   orbitalAndGradientInGrid%values(point,3) = auxMatrix(3)%values(point,j)
                   orbitalAndGradientInGrid%values(point,4) = auxMatrix(4)%values(point,j)
                end do

                ! print *, "viveee"
                ! call Matrix_show(orbitalAndGradientInGrid)

                dftUnit = 77

                write( dftFile, "(A,I0.4)") trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbital_", k

                open(unit = dftUnit, file=trim(dftFile), status="replace", form="unformatted")

                write( labels(1), "(A,I0.4)") "ORBITAL_", k
                labels(2) = Grid_instance(speciesID)%nameOfSpecies

                call Matrix_writeToFile(orbitalAndGradientInGrid, unit=dftUnit, binary=.true., arguments = labels(1:2) )

                close(unit=dftUnit)
                k=k+1
             end do
          end do
       end do

       call Matrix_destructor(auxMatrix(1))
       call Matrix_destructor(auxMatrix(2))
       call Matrix_destructor(auxMatrix(3))
       call Matrix_destructor(auxMatrix(4))
       call Matrix_destructor(orbitalAndGradientInGrid)

    end do

    ! close(unit=dftUnit)

  end subroutine GridManager_writeAtomicOrbitals

  !>
  !! @brief Returns the values of a contracted atomic shell in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getOrbitalGradientAtGrid( this, grid, gridSize, orbital, orbitaldX, orbitaldY, orbitaldZ)
    implicit none
    type(ContractedGaussian) , intent(in) :: this
    type(Matrix) :: grid
    type(Matrix) :: orbital    
    type(Matrix) :: orbitaldX, orbitaldY, orbitaldZ
    integer :: gridSize

    integer :: h
    integer :: nx, ny, nz !< indices de momento angular
    integer :: i, j, m, w
    integer :: point
    real(8) :: coordinate(3)
    real(8) :: exponential, dx, dy, dz
    real(8) :: auxOutput(this%numCartesianOrbital,4)

    do point=1, gridSize
       coordinate(1)=grid%values(point,1)-this%origin(1)
       coordinate(2)=grid%values(point,2)-this%origin(2)
       coordinate(3)=grid%values(point,3)-this%origin(3)
       do h=1, this%length
          exponential=dexp(-this%orbitalExponents(h)*(coordinate(1)**2 + coordinate(2)**2 +coordinate(3)**2) )
          m = 0
          do i = 0 , this%angularMoment
             nx = this%angularMoment - i
             do j = 0 , i
                ny = i - j
                nz = j
                m = m + 1

                !!Orbital
                auxOutput(m,4) = this%contNormalization(m) &
                     * this%primNormalization(h,m) &
                     * coordinate(1)** nx &
                     * coordinate(2)** ny &
                     * coordinate(3)** nz &
                     * exponential 

                dx=-2*this%orbitalExponents(h) &
                     * coordinate(1)** (nx+1) &
                     * coordinate(2)** ny &
                     * coordinate(3)** nz 

                ! Orbital derivative
                if( nx .ge. 1 ) then
                   dx= dx + &
                        nx*coordinate(1)** (nx-1) &
                        * coordinate(2)** ny &
                        * coordinate(3)** nz 
                end if

                dy=-2*this%orbitalExponents(h) &
                     * coordinate(1)** nx &
                     * coordinate(2)** (ny+1) &
                     * coordinate(3)** nz 

                if( ny .ge. 1 ) then
                   dy= dy + &
                        coordinate(1)** nx &
                        *ny*coordinate(2)** (ny-1) &
                        * coordinate(3)** nz 
                end if

                dz=-2*this%orbitalExponents(h) &
                     * coordinate(1)** nx &
                     * coordinate(2)** ny &
                     * coordinate(3)** (nz+1) 

                if( nz .ge. 1 ) then
                   dz= dz+ &
                        coordinate(1)** nx &
                        *coordinate(2)** ny &
                        *nz*coordinate(3)** (nz-1) 
                end if

                auxOutput(m,1) = this%contNormalization(m) &
                     *this%primNormalization(h,m) &
                     *exponential*dx

                auxOutput(m,2) = this%contNormalization(m) &
                     *this%primNormalization(h,m) &
                     *exponential*dy

                auxOutput(m,3) = this%contNormalization(m) &
                     *this%primNormalization(h,m) &
                     *exponential*dz

             end do
          end do

          auxOutput = auxOutput * this%contractionCoefficients(h)

          do w=1, m
             orbital%values(point,w) =   orbital%values(point,w)   + auxOutput(w,4) 
             orbitaldX%values(point,w) = orbitaldX%values(point,w) + auxOutput(w,1) 
             orbitaldY%values(point,w) = orbitaldY%values(point,w) + auxOutput(w,2) 
             orbitaldZ%values(point,w) = orbitaldZ%values(point,w) + auxOutput(w,3) 
          end do
       end do
    end do
  end subroutine GridManager_getOrbitalGradientAtGrid


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
    type(Matrix) :: orbitalUAndGradientInGrid, orbitalVAndGradientInGrid
    integer :: point
    integer :: i, j, u, g
    integer :: ii, jj, v, gg
    integer :: s, ss
    real(8) :: sum
    character(50) ::  dftFile,labels(2)
    integer ::  dftUnit
    integer :: numberOfContractions

    gridSize = Grid_instance(speciesID)%totalSize
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    do u = 1, numberOfContractions

       call Stopwatch_constructor(lowdin_stopwatch)
       call Stopwatch_start(lowdin_stopwatch)

       dftUnit = 77
       write( dftFile, "(A,I0.4)") trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbital_", u

       open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

       write( labels(1), "(A,I0.4)") "ORBITAL_", u
       labels(2) = Grid_instance(speciesID)%nameOfSpecies

       orbitalUAndGradientInGrid= Matrix_getFromFile(unit=dftUnit, rows= int(gridSize,4), &
            columns= int(4,4), binary=.true., arguments=labels)

       ! print *, "orbital", u
       ! call Matrix_show( orbitalUAndGradientInGrid)

       close(unit=dftUnit)

       do point = 1 , gridSize

          densityInGrid%values(point)=densityInGrid%values(point)+densityMatrix%values(u,u)*orbitalUAndGradientInGrid%values(point,1)**2
          gradientInGrid(1)%values(point)=gradientInGrid(1)%values(point)+2*densityMatrix%values(u,u)&
               *orbitalUAndGradientInGrid%values(point,1)*orbitalUAndGradientInGrid%values(point,2)
          gradientInGrid(2)%values(point)=gradientInGrid(2)%values(point)+2*densityMatrix%values(u,u)&
               *orbitalUAndGradientInGrid%values(point,1)*orbitalUAndGradientInGrid%values(point,3)
          gradientInGrid(3)%values(point)=gradientInGrid(3)%values(point)+2*densityMatrix%values(u,u)&
               *orbitalUAndGradientInGrid%values(point,1)*orbitalUAndGradientInGrid%values(point,4)

       end do


       do v = u+1, numberOfContractions

          dftUnit = 77
          write( dftFile, "(A,I0.4)") trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbital_", v

          open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

          write( labels(1), "(A,I0.4)") "ORBITAL_", v
          labels(2) = Grid_instance(speciesID)%nameOfSpecies

          orbitalVAndGradientInGrid= Matrix_getFromFile(unit=dftUnit, rows= int(gridSize,4), &
               columns= int(4,4), binary=.true., arguments=labels(1:2))

          close(unit=dftUnit)

          do point = 1 , gridSize

             densityInGrid%values(point)=densityInGrid%values(point)+2*densityMatrix%values(u,v)*&
                  orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,1)

             gradientInGrid(1)%values(point)=gradientInGrid(1)%values(point)+2*densityMatrix%values(u,v)*&
                  (orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,2)+&
                  orbitalUAndGradientInGrid%values(point,2)*orbitalVAndGradientInGrid%values(point,1))

             gradientInGrid(2)%values(point)=gradientInGrid(2)%values(point)+2*densityMatrix%values(u,v)*&
                  (orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,3)+&
                  orbitalUAndGradientInGrid%values(point,3)*orbitalVAndGradientInGrid%values(point,1))

             gradientInGrid(3)%values(point)=gradientInGrid(3)%values(point)+2*densityMatrix%values(u,v)*&
                  (orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,4)+&
                  orbitalUAndGradientInGrid%values(point,4)*orbitalVAndGradientInGrid%values(point,1))

          end do
       end do
    end do

    ! call Vector_show(gradientInGrid(1))
    ! call Vector_show(gradientInGrid(2))
    ! call Vector_show(gradientInGrid(3))

  end subroutine GridManager_getDensityGradientAtGrid

  !>
  !! @brief Returns the values of the exchange correlation potential for a specie in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getEnergyAndPotentialAtGrid( speciesID, exchangeCorrelationEnergy, potentialInGrid, sigmaPotentialInGrid,&
       otherSpeciesID, otherExchangeCorrelationEnergy, otherPotentialInGrid, otherSigmaPotentialInGrid)
    implicit none
    integer :: speciesID
    real(8) :: exchangeCorrelationEnergy
    type(Vector) :: potentialInGrid
    type(Vector) :: sigmaPotentialInGrid
    integer, optional :: otherSpeciesID
    real(8), optional :: otherExchangeCorrelationEnergy
    type(Vector), optional :: otherPotentialInGrid
    type(Vector), optional :: otherSigmaPotentialInGrid

    character(50) :: nameOfSpecies, otherNameOfSpecies
    integer :: gridSize, otherGridSize
    type(Vector) :: energyDensity
    type(Vector) :: sigma
    type(Vector) :: densityAB, potentialAB, sigmaAB, sigmaPotentialAB
    type(Vector) :: commonPoints, electronicDensityAtOtherGrid, electronicPotentialAtOtherGrid, holdNuclearPotential
    integer :: i, index
    real(8) :: nuclearElectronCorrelationEnergy

    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )
    if( present(otherSpeciesID) )     otherNameOfSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )

    gridSize = Grid_instance(speciesID)%totalSize

    call Vector_Constructor(energyDensity, gridSize, 0.0_8)

    !Closed Shell
    if (nameOfSpecies=="E-" .and. .not. present(otherSpeciesID) ) then

       if (CONTROL_instance%CALL_LIBXC) then

          index=Functional_getIndex(speciesID)

          !libxc works with the gradient squared - sigma
          call Vector_Constructor(sigma, gridSize, 0.0_8)         
          do i=1, gridSize
             sigma%values(i)=(Grid_instance(speciesID)%densityGradient(1)%values(i)**2 + Grid_instance(speciesID)%densityGradient(2)%values(i)**2 + Grid_instance(speciesID)%densityGradient(3)%values(i)**2)
          end do

          call Functional_libxcEvaluate(Functionals(index), gridSize, Grid_instance(speciesID)%density%values, sigma%values, energyDensity%values , potentialInGrid%values, sigmaPotentialInGrid%values )

          ! print *, "sigma"
          ! call Vector_Show(sigma)
          call Vector_Destructor(sigma)

          do i=1, gridSize
             exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
          end do

          ! print *, "electronicEXC RKS", exchangeCorrelationEnerg

       else

          index=Functional_getIndex(speciesID)

          if ( Functionals(index)%name .eq. "exchange:Slater-correlation:VWN5") then
             call Functional_LDAEvaluate(gridSize, Grid_instance(speciesID)%density%values/2, Grid_instance(speciesID)%density%values/2, &
                  energyDensity%values, potentialInGrid%values )

             do i=1, gridSize
                exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%points%values(i,4) 
             end do
          end if

          ! print *, "electronicEXC RKS", exchangeCorrelationEnergy
       end if

    elseif (nameOfSpecies=="E-ALPHA" .and. otherNameOfSpecies=="E-BETA") then

       if (CONTROL_instance%CALL_LIBXC) then

          call Vector_Constructor(densityAB, 2*gridSize, 0.0_8)         
          call Vector_Constructor(sigmaAB, 3*gridSize, 0.0_8)         
          call Vector_Constructor(potentialAB, 2*gridSize, 0.0_8)         
          call Vector_Constructor(sigmaPotentialAB, 3*gridSize, 0.0_8)         


          do i=1, gridSize
             densityAB%values(2*i-1)=Grid_instance(speciesID)%density%values(i)
             
             densityAB%values(2*i)=Grid_instance(otherSpeciesID)%density%values(i)

             sigmaAB%values(3*i-2)=(Grid_instance(speciesID)%densityGradient(1)%values(i)**2 + Grid_instance(speciesID)%densityGradient(2)%values(i)**2 + &
                  Grid_instance(speciesID)%densityGradient(3)%values(i)**2)

             sigmaAB%values(3*i-1)=Grid_instance(speciesID)%densityGradient(1)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(1)%values(i)+&
                  Grid_instance(speciesID)%densityGradient(2)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(2)%values(i)+&
                  Grid_instance(speciesID)%densityGradient(3)%values(i)*Grid_instance(otherSpeciesID)%densityGradient(3)%values(i)

             sigmaAB%values(3*i)=(Grid_instance(otherSpeciesID)%densityGradient(1)%values(i)**2 + Grid_instance(otherSpeciesID)%densityGradient(2)%values(i)**2 + &
                  Grid_instance(otherSpeciesID)%densityGradient(3)%values(i)**2)

          end do

          index=Functional_getIndex(speciesID)

          call Functional_libxcEvaluate(Functionals(index), gridSize, densityAB%values, sigmaAB%values, energyDensity%values , potentialAB%values, sigmaPotentialAB%values )

          print *, "density", densityAB%values(1), densityAB%values(2*gridSize)
          call Vector_show(densityAB)
          call Vector_Destructor(densityAB)

          ! print *, "sigma"
          ! call Vector_show(sigmaAB)
          call Vector_Destructor(sigmaAB)

          
          do i=1, gridSize
             potentialInGrid%values(i)=potentialInGrid%values(i)+potentialAB%values(2*i-1)
             
             otherPotentialInGrid%values(i)=otherPotentialInGrid%values(i)+potentialAB%values(2*i)
             
             sigmaPotentialInGrid%values(i)=sigmaPotentialInGrid%values(i)+sigmaPotentialAB%values(3*i-2)!+sigmaPotentialAB%values(3*i-1)
             
             otherSigmaPotentialInGrid%values(i)=otherSigmaPotentialInGrid%values(i)+sigmaPotentialAB%values(3*i)!+sigmaPotentialAB%values(3*i-1)

             exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4)
             
             otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(otherSpeciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
          end do

          ! print *, "electronicEXC UKS", exchangeCorrelationEnergy, otherExchangeCorrelationEnergy
          call Vector_Destructor(potentialAB)         
          call Vector_Destructor(sigmaPotentialAB)         
          
       else
          index=Functional_getIndex(speciesID)

          if ( Functionals(index)%name .eq. "exchange:Slater-correlation:VWN5") then

             call Functional_LDAEvaluate(gridSize, Grid_instance(speciesID)%density%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, potentialInGrid%values, otherPotentialInGrid%values )

             do i=1, gridSize
                exchangeCorrelationEnergy=exchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(speciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
                otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+energyDensity%values(i)*Grid_instance(otherSpeciesID)%density%values(i)*Grid_instance(speciesID)%points%values(i,4) 
             end do

             ! exchangeCorrelationEnergy=exchangeCorrelationEnergy*(MolecularSystem_getNumberOfParticles( speciesID ))/(MolecularSystem_getNumberOfParticles( speciesID )+MolecularSystem_getNumberOfParticles( otherSpeciesID ))
             ! otherExchangeCorrelationEnergy=exchangeCorrelationEnergy*(MolecularSystem_getNumberOfParticles( otherSpeciesID ))/(MolecularSystem_getNumberOfParticles( speciesID )+MolecularSystem_getNumberOfParticles( otherSpeciesID ))
             ! exchangeCorrelationEnergy=exchangeCorrelationEnergy/2
             ! otherExchangeCorrelationEnergy=exchangeCorrelationEnergy/2

             ! print *, "electronicEXC UKS", exchangeCorrelationEnergy, otherExchangeCorrelationEnergy
          end if

       end if

       !Closed shell nuclear electron correlation
    elseif ( (nameOfSpecies=="E-") .and. present(otherSpeciesID)  ) then

       index=Functional_getIndex(speciesID, otherSpeciesID)
       
       if(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL .ne. "NONE" &
            .and. .not. (otherNameOfSpecies=="E-" .or. otherNameOfSpecies=="E-ALPHA" .or. otherNameOfSpecies=="E-BETA")) then

          otherGridSize = Grid_instance(otherSpeciesID)%totalSize
          otherNameOfSpecies =trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))

          call Vector_Constructor(energyDensity, otherGridSize, 0.0_8)
          call Vector_constructor(electronicDensityAtOtherGrid, otherGridSize, 0.0_8)
          call Vector_constructor(electronicPotentialAtOtherGrid, otherGridSize, 0.0_8)
          call Vector_constructor(commonPoints, otherGridSize, 0.0_8) !!To build the electronic potential with the correct size
          
          call GridManager_getElectronicDensityInOtherGrid(speciesID, otherSpeciesID, commonPoints, electronicDensityAtOtherGrid )

          select case (trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL) )
             
          case ("epc17-1")
             call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("epc17-2")
             call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("ikn-nsf")
             call Functional_IKNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("mlcs-fit")
             call Functional_MLCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("mlcs-a")
             call Functional_MLCSAEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("mlcs-an")
             call Functional_MLCSANEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("CS-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("Imamura-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("Mejia-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("MejiaA-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("expCS-A")
             call Functional_expCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("psn")
             call Functional_PSNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case default
             ! print *, trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
             STOP "The nuclear electron functional chosen is not implemented"

          end select
          
          ! call Functional_lowLimitEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
          !      electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
          !      energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          !!Adds the nuclear electron potential to the relevant points in the electronic grid
          ! print *, "i, electronicPotentialAtOtherGrid%values(i), Grid_instance(otherSpeciesID)%density%values(i), energyDensity%values(i)"
          do i=1, otherGridSize
             index=int(commonPoints%values(i),4)
             potentialInGrid%values(index) = potentialInGrid%values(index) + electronicPotentialAtOtherGrid%values(i)
             ! print *, i, electronicDensityAtOtherGrid%values(i), Grid_instance(otherSpeciesID)%density%values(i), energyDensity%values(i)
          end do
          
          nuclearElectronCorrelationEnergy=0.0_8

          do i=1, otherGridSize
             nuclearElectronCorrelationEnergy=nuclearElectronCorrelationEnergy+energyDensity%values(i)*Grid_instance(otherSpeciesID)%points%values(i,4) 
          end do

          ! print *, "nuclear electron correlation energy", nuclearElectronCorrelationEnergy

          exchangeCorrelationEnergy=exchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2
          otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2
          
          ! STOP "trolololoooooo"
          
       end if

    elseif ( (nameOfSpecies=="E-ALPHA" .or. nameOfSpecies=="E-BETA") .and. present(otherSpeciesID)  ) then

       index=Functional_getIndex(speciesID, otherSpeciesID)
       
       if(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL .ne. "NONE" &
            .and. .not. (otherNameOfSpecies=="E-" .or. otherNameOfSpecies=="E-ALPHA" .or. otherNameOfSpecies=="E-BETA")) then

          otherGridSize = Grid_instance(otherSpeciesID)%totalSize
          otherNameOfSpecies =trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))

          call Vector_Constructor(energyDensity, otherGridSize, 0.0_8)
          call Vector_constructor(electronicDensityAtOtherGrid, otherGridSize, 0.0_8)
          call Vector_constructor(electronicPotentialAtOtherGrid, otherGridSize, 0.0_8)

          call Vector_constructor(commonPoints, otherGridSize, 0.0_8) !!To build the electronic potential with the correct size

          call GridManager_getElectronicDensityInOtherGrid(speciesID, otherSpeciesID, commonPoints, electronicDensityAtOtherGrid )

          !!!This is a very dirty way of preventing double counting of the nuclear potential
          if(nameOfSpecies .eq. "E-BETA") call Vector_copyConstructor (holdNuclearPotential,otherPotentialInGrid)

          select case (trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL) )
             
          case ("epc17-1")
             call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("epc17-2")
             call Functional_EPCEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("ikn-nsf")
             call Functional_IKNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("mlcs-fit")
             call Functional_MLCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("mlcs-a")
             call Functional_MLCSAEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("mlcs-an")
             call Functional_MLCSANEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
                  electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
                  energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("CS-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("Imamura-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("Mejia-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("MejiaA-myfit")
             call Functional_myCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case ("expCS-A")
             call Functional_expCSEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
             
          case ("psn")
             call Functional_PSNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
               electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
               energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          case default
             print *, trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
             STOP "The nuclear electron functional chosen is not implemented"

          end select

          
          ! call Functional_lowLimitEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
          !      electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
          !      energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )
          
          ! call Functional_PSNEvaluate(Functionals(index), MolecularSystem_getMass( otherSpeciesID ), otherGridSize, &
          !      electronicDensityAtOtherGrid%values, Grid_instance(otherSpeciesID)%density%values, &
          !      energyDensity%values, electronicPotentialAtOtherGrid%values, otherPotentialInGrid%values  )

          if(nameOfSpecies .eq. "E-BETA") call Vector_copyConstructor (otherPotentialInGrid, holdNuclearPotential)

          !!Adds the nuclear electron potential to the relevant points in the electronic grid
          ! print *, "i, electronicPotentialAtOtherGrid%values(i), Grid_instance(otherSpeciesID)%density%values(i), energyDensity%values(i)"
          do i=1, otherGridSize
             index=int(commonPoints%values(i),4)
             if(i .le. gridSize) potentialInGrid%values(index) = potentialInGrid%values(index) + electronicPotentialAtOtherGrid%values(i)
             ! print *, i, electronicDensityAtOtherGrid%values(i), Grid_instance(otherSpeciesID)%density%values(i), energyDensity%values(i)
          end do
          
          nuclearElectronCorrelationEnergy=0.0_8

          do i=1, otherGridSize
             nuclearElectronCorrelationEnergy=nuclearElectronCorrelationEnergy+energyDensity%values(i)*Grid_instance(otherSpeciesID)%points%values(i,4) 
          end do

          ! print *, "nuclear electron correlation energy", nuclearElectronCorrelationEnergy

          !!!This is also a dirty way of preventing double counting of the nuclear electron correlation energy
          if(nameOfSpecies .eq. "E-ALPHA") then
             exchangeCorrelationEnergy=exchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2/2          
             otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2/2

          else if(nameOfSpecies .eq. "E-BETA") then
             exchangeCorrelationEnergy=exchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2/2          
             otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2/2
          else
             exchangeCorrelationEnergy=exchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2        
             otherExchangeCorrelationEnergy=otherExchangeCorrelationEnergy+nuclearElectronCorrelationEnergy/2
          end if
          ! STOP "trolololoooooo"

       end if
       
    end if


    call Vector_Destructor(energyDensity)

    ! do i=1, gridSize
    !    print *, densityInGrid%values(i), exchange%values(i), correlationA%values(i)
    ! end do

    ! call Vector_Destructor(correlationA)
    ! call Vector_Destructor(correlationB)


  end subroutine GridManager_getEnergyAndPotentialAtGrid

  !>
  !! @brief Builds the exchange correlation for a species
  ! Felix Moncada, 2017
  subroutine GridManager_buildExchangeCorrelationMatrix( speciesID, exchangeCorrelationMatrix  )
    implicit none
    integer :: speciesID
    type(Matrix) :: exchangeCorrelationMatrix

    integer :: gridSize
    integer :: numberOfContractions

    type(Matrix) :: orbitalUAndGradientInGrid, orbitalVAndGradientInGrid
    integer :: numberOfCartesiansOrbitalsU
    integer :: numberOfCartesiansOrbitalsV
    integer :: u, v, point

    character(50) :: dftFile, labels(2)
    integer :: dftUnit

    real(8) :: time1, time2

    time1=0.0_8
    time2=0.0_8

    gridSize = Grid_instance(speciesID)%totalSize
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    do u = 1, numberOfContractions

       ! call Stopwatch_constructor(lowdin_stopwatch)
       ! call Stopwatch_start(lowdin_stopwatch)

       dftUnit = 77
       write( dftFile, "(A,I0.4)") trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbital_", u

       open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

       write( labels(1), "(A,I0.4)") "ORBITAL_", u
       labels(2) = Grid_instance(speciesID)%nameOfSpecies

       orbitalUAndGradientInGrid= Matrix_getFromFile(unit=dftUnit, rows= int(gridSize,4), &
            columns= int(4,4), binary=.true., arguments=labels)

       ! print *, "orbital", u
       ! call Matrix_show( orbitalUAndGradientInGrid)

       close(unit=dftUnit)

       ! call Stopwatch_stop(lowdin_stopwatch)     
       ! time1=time1+lowdin_stopwatch%enlapsetTime

       ! call Stopwatch_constructor(lowdin_stopwatch)
       ! call Stopwatch_start(lowdin_stopwatch)

       do point = 1 , gridSize
          exchangeCorrelationMatrix%values(u,u)=&
               exchangeCorrelationMatrix%values(u,u)&
               +(Grid_instance(speciesID)%potential%values(point)*orbitalUAndGradientInGrid%values(point,1)**2&
               +4*Grid_instance(speciesID)%sigmaPotential%values(point)*orbitalUAndGradientInGrid%values(point,1)&
               *(Grid_instance(speciesID)%densityGradient(1)%values(point)*orbitalUAndGradientInGrid%values(point,2)&
               +Grid_instance(speciesID)%densityGradient(2)%values(point)*orbitalUAndGradientInGrid%values(point,3)&
               +Grid_instance(speciesID)%densityGradient(3)%values(point)*orbitalUAndGradientInGrid%values(point,4))&
               )*Grid_instance(speciesID)%points%values(point,4)
       end do

       call Stopwatch_stop(lowdin_stopwatch)     
       ! time2=time2+lowdin_stopwatch%enlapsetTime

       do v = u+1, numberOfContractions

          call Stopwatch_constructor(lowdin_stopwatch)
          call Stopwatch_start(lowdin_stopwatch)

          dftUnit = 77
          write( dftFile, "(A,I0.4)") trim(CONTROL_instance%INPUT_FILE)//trim(Grid_instance(speciesID)%nameOfSpecies)//".orbital_", v

          open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

          write( labels(1), "(A,I0.4)") "ORBITAL_", v
          labels(2) = Grid_instance(speciesID)%nameOfSpecies

          orbitalVAndGradientInGrid= Matrix_getFromFile(unit=dftUnit, rows= int(gridSize,4), &
               columns= int(4,4), binary=.true., arguments=labels(1:2))

          close(unit=dftUnit)

          ! call Stopwatch_stop(lowdin_stopwatch)     
          ! time1=time1+lowdin_stopwatch%enlapsetTime

          ! call Stopwatch_constructor(lowdin_stopwatch)
          ! call Stopwatch_start(lowdin_stopwatch)

          do point = 1 , gridSize
             exchangeCorrelationMatrix%values(u,v)=&
                  exchangeCorrelationMatrix%values(u,v)&
                  +(Grid_instance(speciesID)%potential%values(point)&
                  *orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,1)&
                  +2*Grid_instance(speciesID)%sigmaPotential%values(point)*&
                  (Grid_instance(speciesID)%densityGradient(1)%values(point)&
                  *(orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,2)&
                  +orbitalUAndGradientInGrid%values(point,2)*orbitalVAndGradientInGrid%values(point,1))&
                  +Grid_instance(speciesID)%densityGradient(2)%values(point)&
                  *(orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,3)&
                  +orbitalUAndGradientInGrid%values(point,3)*orbitalVAndGradientInGrid%values(point,1))&
                  +Grid_instance(speciesID)%densityGradient(3)%values(point)&
                  *(orbitalUAndGradientInGrid%values(point,1)*orbitalVAndGradientInGrid%values(point,4)&
                  +orbitalUAndGradientInGrid%values(point,4)*orbitalVAndGradientInGrid%values(point,1))&
                  ))*Grid_instance(speciesID)%points%values(point,4)

          end do

          ! call Stopwatch_stop(lowdin_stopwatch)     
          ! time2=time2+lowdin_stopwatch%enlapsetTime

       end do
    end do

    do u=1, numberOfContractions
       do v=u+1, numberOfContractions
          exchangeCorrelationMatrix%values(v,u)=exchangeCorrelationMatrix%values(u,v)
       end do
    end do

    ! call Stopwatch_stop(lowdin_stopwatch)     
    ! write(*,"(A,F10.3,A4)") "** reading orbital files:", time1 ," (s)"
    ! write(*,"(A,F10.3,A4)") "** integrating over the grid:", time2 ," (s)"

  end subroutine GridManager_buildExchangeCorrelationMatrix



  subroutine GridManager_getElectronicDensityInOtherGrid(electronicID,otherSpeciesID, commonPoints, electronicDensityAtOtherGrid )
    integer :: electronicID, otherSpeciesID
    type(Vector) :: commonPoints, electronicDensityAtOtherGrid

    character(50) :: nameOfElectron
    integer :: otherElectronicID
    integer :: electronicGridSize, otherGridSize
    integer :: i,j

    electronicGridSize=Grid_instance(electronicID)%totalSize
    otherGridSize= Grid_instance(otherSpeciesID)%totalSize

    nameOfElectron=MolecularSystem_getNameOfSpecie(electronicID)
    if (nameOfElectron .eq. "E-ALPHA") otherElectronicID=MolecularSystem_getSpecieID( "E-BETA" )
    if (nameOfElectron .eq. "E-BETA") otherElectronicID=MolecularSystem_getSpecieID( "E-ALPHA" )
       
    
    call Vector_constructor(commonPoints, otherGridSize, 0.0_8)
    call Vector_constructor(electronicDensityAtOtherGrid, otherGridSize, 1.0E-12_8)

    !The other grid must be a subset of the electronic grid
    !FELIX: This is a problem for positron calculations
    do i=1, electronicGridSize
       do j=1, otherGridSize
          if(Grid_instance(electronicID)%points%values(i,1) .eq. Grid_instance(otherSpeciesID)%points%values(j,1) .and. &
               Grid_instance(electronicID)%points%values(i,2) .eq. Grid_instance(otherSpeciesID)%points%values(j,2) .and. &
               Grid_instance(electronicID)%points%values(i,3) .eq. Grid_instance(otherSpeciesID)%points%values(j,3) ) then 
             commonPoints%values(j)=i
             if(nameOfElectron .eq. "E-ALPHA" .or. nameOfElectron .eq. "E-BETA") then
                electronicDensityAtOtherGrid%values(j)=Grid_instance(electronicID)%density%values(i)+Grid_instance(otherElectronicID)%density%values(i)
             else
                electronicDensityAtOtherGrid%values(j)=Grid_instance(electronicID)%density%values(i)
             end if
          end if
       end do
    end do
  end subroutine GridManager_getElectronicDensityInOtherGrid
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
!   @brief Returns the values of all atomic orbitals in a set of coordinates
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

  
