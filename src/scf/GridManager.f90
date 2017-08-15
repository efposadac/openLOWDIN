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
!! @brief This module manages the orbital and density represented in the DFT grids. Partially based on R. Flores-Moreno Parakata's modules
!! @author F. Moncada, 2017
module GridManager_
  use Matrix_
  use Exception_
  use String_
  use MolecularSystem_
  use Grid_
  use Functional_
  implicit none

  public :: &
       GridManager_buildGrids, &
       GridManager_atomicOrbitals, &
       GridManager_getOrbitalMatrix, &
       GridManager_getOrbitalAtGrid, &
       GridManager_getDensityAtGrid, &
       GridManager_createFunctionals, &
       GridManager_getPotentialAtGrid, &
       GridManager_getEnergyFromGrid

contains

  !>
  !! @brief Builds a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine GridManager_buildGrids( )
    implicit none
    type(Matrix) :: atomicGrid
    integer :: numberOfSpecies, numberOfCenters, atomGridSize
    integer :: speciesID, particleID, particleID2, particleID3, point, i
    character(50) :: labels(2), dftFile
    integer :: dftUnit
      
    real(8) :: cutoff, sum, r, w, mu
    real(8), allocatable :: origins(:,:), distance(:),factor(:)

    !! Open file for dft
    ! dftUnit = 77
    ! dftFile = trim(CONTROL_instance%INPUT_FILE)//"dft"
    ! open(unit = dftUnit, file=trim(dftFile), status="new", form="formatted")

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate memory.
    allocate(Grid_instance(numberOfSpecies))
    
    !! Allocate memory for specie in system and load some matrices.
    do speciesID = 1, numberOfSpecies

       call Grid_constructor(Grid_instance(speciesID), speciesID )

    end do
    
  end subroutine GridManager_buildGrids


  subroutine GridManager_atomicOrbitals( )
    implicit none
    
    integer :: numberOfSpecies
    integer :: totalNumberOfContractions
    integer :: speciesID
    integer :: gridSize
    integer :: mu,nu, point, index
    type(Matrix) :: grid
    type(Matrix) :: orbitalsInGrid

    character(50) :: labels(2), dftFile
    integer :: dftUnit

    !! Open file for dft
    ! dftUnit = 77
    ! dftFile = trim(CONTROL_instance%INPUT_FILE)//"dft"
    ! open(unit = dftUnit, file=trim(dftFile), status="new", access='sequential', form="formatted")

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do speciesID = 1, numberOfSpecies

       totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
       grid=Grid_instance(speciesID)%points
       gridSize = Grid_instance(speciesID)%totalSize

       call Matrix_Constructor( Grid_instance(speciesID)%orbitals, int(gridSize,8), int(totalNumberOfContractions,8), 0.0_8)

       call GridManager_getOrbitalMatrix( speciesID, grid, gridSize, Grid_instance(speciesID)%orbitals)

      !!!!Esto es un abuso de memoria!!! escala como por N^5
       ! call Matrix_Constructor( Grid_instance(speciesID)%orbitalsProduct, int(gridSize,8), int(totalNumberOfContractions*(totalNumberOfContractions+1)/2,8), 0.0_8)

       ! index=1
       ! do mu=1, totalNumberOfContractions
       !    do nu=mu, totalNumberOfContractions
       !       do point=1, gridSize
       !          Grid_instance(speciesID)%orbitalsProduct%values(point,index)=orbitalsInGrid%values(point,mu)*orbitalsInGrid%values(point,nu)
       !       end do
       !       index=index+1
       !    end do
       ! end do

       !Write to disk
       ! labels(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
       ! labels(1) = "ORBITAL-GRID"
       ! call Matrix_show(orbitalsInGrid)
       ! call Matrix_writeToFile( orbitalsInGrid, unit=dftUnit, binary=.false., arguments = labels(1:2) )
       ! call Matrix_destructor( orbitalsInGrid)
       ! labels(1) = "ORBITALPRODUCT-GRID"
       ! call Matrix_show(orbitalsProductInGrid)
       ! call Matrix_writeToFile(orbitalsProductInGrid, unit=dftUnit, binary=.false., arguments = labels(1:2) )
       ! call Matrix_destructor(orbitalsProductInGrid)

    end do

    ! close(unit=dftUnit)
    
  end subroutine GridManager_atomicOrbitals
  
    subroutine GridManager_createFunctionals()
    implicit none

    integer :: numberOfSpecies
    integer :: speciesID, otherSpeciesID, i

    !! Open file for dft
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate(Functionals(numberOfSpecies*2+numberOfSpecies*(numberOfSpecies-1)/2))

    i=1
    do speciesID=1, numberOfSpecies
       call Functional_constructor(Functionals(i),"exchange",speciesID,speciesID)
       i=i+1
       call Functional_constructor(Functionals(i),"correlation",speciesID,speciesID)
       i=i+1
    end do

    do speciesID=1, numberOfSpecies-1
       do otherSpeciesID=speciesID+1, numberOfSpecies  
          call Functional_constructor(Functionals(i),"correlation",speciesID,otherSpeciesID)
          i=i+1
       end do
    end do

    print *, ""
    print *, "-------------------------------------------------"
    print *, "|-----------Functionals summary -----------------|"
    print *, "-------------------------------------------------"
    print *, ""

    i=1
    do i=1, size(Functionals)
       call Functional_show(Functionals(i))
    end do

    print *, "-------------------------------------------------"
    
  end subroutine GridManager_createFunctionals
  

  !>
  !! @brief Returns the values of all the atomic orbitals in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getOrbitalMatrix( speciesID, grid, gridSize, orbitalsInGrid)
    implicit none
    integer :: speciesID
    type(Matrix) :: grid, orbitalsInGrid
    integer :: gridSize

    type(Matrix) :: auxMatrix
    integer :: numberOfCartesiansOrbitals
    integer :: totalNumberOfContractions
    integer :: point
    integer :: i, j, k, g

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    k=1
    do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
       do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)
          numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital

          call Matrix_constructor( auxMatrix, int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8)

          call GridManager_getOrbitalAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), grid, gridSize, auxMatrix)
          
          do j = 1, numberOfCartesiansOrbitals
             do point = 1 , gridSize
                orbitalsInGrid%values(point,k) = auxMatrix%values(point,j)
             end do
             k=k+1
          end do
       end do
    end do
  end subroutine GridManager_getOrbitalMatrix

  !>
  !! @brief Returns the values of a contracted atomic shell in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getOrbitalAtGrid( this, grid, gridSize, output)
    implicit none
    type(ContractedGaussian) , intent(in) :: this
    type(Matrix) :: grid
    type(Matrix) :: output
    integer :: gridSize
    
    integer :: h
    integer :: nx, ny, nz !< indices de momento angular
    integer :: i, j, m, xx
    integer :: point
    real(8) :: coordinate(3)
    real(8) :: exponential
    real(8) :: auxOutput(this%numCartesianOrbital)

    do point=1, gridSize
       coordinate(1)=grid%values(point,1)
       coordinate(2)=grid%values(point,2)
       coordinate(3)=grid%values(point,3)
       do h=1, this%length
          exponential=dexp(-this%orbitalExponents(h) &
               *(  (this%origin(1)-coordinate(1))**2 &
               +(this%origin(2)-coordinate(2))**2 &
               +(this%origin(3)-coordinate(3))**2) )
          m = 0
          do i = 0 , this%angularMoment
             nx = this%angularMoment - i
             do j = 0 , i
                ny = i - j
                nz = j
                m = m + 1
                auxOutput(m) = this%contNormalization(m) &
                     * this%primNormalization(h,m) &
                     * (coordinate(1)-this%origin(1))** nx &
                     * (coordinate(2)-this%origin(2))** ny &
                     * (coordinate(3)-this%origin(3))** nz &
                     * exponential 
             end do
          end do

          do xx=1, m
             output%values(point,xx) = output%values(point,xx) + auxOutput(xx) * this%contractionCoefficients(h)
          end do
       end do
    end do
  end subroutine GridManager_getOrbitalAtGrid

  
      !>
  !! @brief Returns the values of the density in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getDensityAtGrid( speciesID, gridSize, densityMatrix, densityInGrid)
    implicit none
    integer :: speciesID
    integer :: gridSize
    type(Matrix) :: densityMatrix
    type(Vector) :: densityInGrid

    character(50) :: nameOfSpecies
    character(50) :: labels(2), dftFile
    ! type(Matrix) :: orbitalsProductInGrid
    integer :: numberOfContractions
    integer :: u, v, i, index
    integer :: dftUnit, wfnUnit
    real(8) :: sum
    
    !! Open file for dft

    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

    !Read grid and orbital products from disk
    ! dftUnit = 77
    ! dftFile = trim(CONTROL_instance%INPUT_FILE)//"dft"

    ! labels(2) = trim(nameOfSpecies)
    ! labels(1) = "ORBITALPRODUCT-GRID"
    ! open( UNIT=dftUnit,FILE=dftFile, status='old', access='sequential', form='formatted')
    ! orbitalsProductInGrid=Matrix_getFromFile(unit=dftUnit, rows=gridSize, columns=numberOfContractions*(numberOfContractions+1)/2, binary=.false., arguments=labels)
    ! close( dftUnit)

    ! call Matrix_show(orbitalsProductInGrid)
    ! call Matrix_show(densityMatrix)
    
    do u = 1 , numberOfContractions
       do v = u , numberOfContractions
          if ( u .eq. v) then
             do i=1,gridSize
                densityInGrid%values(i)=densityInGrid%values(i)+densityMatrix%values(u,v)*&
                     Grid_instance(speciesID)%orbitals%values(i,u)*Grid_instance(speciesID)%orbitals%values(i,v)
             end do
          else
             do i=1,gridSize
                densityInGrid%values(i)=densityInGrid%values(i)+2*densityMatrix%values(u,v)*&
                     Grid_instance(speciesID)%orbitals%values(i,u)*Grid_instance(speciesID)%orbitals%values(i,v)
             end do
          end if
       end do
    end do

    !Check normalization
    ! sum=0
    ! do i=1,gridSize
    !    print *, i, densityInGrid%values(i) 
    !    sum=sum+densityInGrid%values(i)*Grid_instance(speciesID)%points%values(i,4)
    ! end do
    
    ! print *, "particles", sum
    
  end subroutine GridManager_getDensityAtGrid


  !>
  !! @brief Returns the values of the exchange correlation potential for a specie in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getPotentialAtGrid( speciesID, gridSize, densityInGrid, potentialInGrid)
    implicit none
    integer :: speciesID
    type(Vector) :: densityInGrid
    type(Vector) :: potentialInGrid
    integer :: gridSize
    character(50) :: nameOfSpecies

    type(Matrix) :: grid
    type(Vector) :: exchange,correlationA,correlationB
    integer :: dftUnit
    integer :: numberOfContractions
    integer :: i
    integer :: u, v

    !! Open file for dft

    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )

    call Vector_Constructor(exchange, gridSize, 0.0_8)
    call Vector_Constructor(correlationA, gridSize, 0.0_8)
    call Vector_Constructor(correlationB, gridSize, 0.0_8)
    
    if (nameOfSpecies=="E-") then

       call vxdirac(densityInGrid%values/2, exchange%values, gridSize)

       ! call Vector_show(exchange)

       call vcvwn(densityInGrid%values/2, densityInGrid%values/2 , correlationA%values, correlationB%values, gridSize)

       ! call Vector_show(correlationA)

       potentialInGrid%values=correlationA%values+exchange%values

       ! call Vector_show(potentialInGrid)
       
    end if

    ! do i=1, gridSize
    !    print *, densityInGrid%values(i), exchange%values(i), correlationA%values(i)
    ! end do
    
  end subroutine GridManager_getPotentialAtGrid

  
  !>
  !! @brief Returns the values of the exchange correlation potential for a specie in a set of coordinates
!!! Felix Moncada, 2017
  !<
  subroutine GridManager_getEnergyFromGrid( speciesID, gridSize, densityInGrid, exchangeCorrelationEnergy)
    implicit none
    integer :: speciesID
    type(Vector) :: densityInGrid
    real(8) :: exchangeCorrelationEnergy
    integer :: gridSize
    character(50) :: nameOfSpecies

    real(8) :: exchangeEnergy,correlationEnergy,particles
    type(Matrix) :: grid
    type(Vector) :: exchange,correlation
    integer :: dftUnit
    integer :: numberOfContractions
    integer :: i
    integer :: u, v

    !! Open file for dft
 
    grid=Grid_instance(speciesID)%points
    
    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )

    call Vector_Constructor(exchange, gridSize, 0.0_8)
    call Vector_Constructor(correlation, gridSize, 0.0_8)

    if (nameOfSpecies=="E-") then

       call exdirac(densityInGrid%values/2, exchange%values, gridSize)

       ! call Vector_show(exchange)

       call ecvwn(densityInGrid%values/2, densityInGrid%values/2 , correlation%values, gridSize)

       ! call Vector_show(correlationA)
       
    end if

    particles=0.0
    exchangeEnergy=0.0
    correlationEnergy=0.0

    do i=1, gridSize
       particles=particles+densityInGrid%values(i)*grid%values(i,4)
       exchangeEnergy=exchangeEnergy+exchange%values(i)*grid%values(i,4)
       correlationEnergy=correlationEnergy+correlation%values(i)*grid%values(i,4)
    end do
    
    exchangeCorrelationEnergy=exchangeEnergy+correlationEnergy

    print *, "particles", particles 
    print *, "exchange", exchangeEnergy
    print *, "correlation", correlationEnergy
    print *, "total", exchangeCorrelationEnergy
    
  end subroutine GridManager_getEnergyFromGrid

  
end module GridManager_


  
!   !<
!   !! @brief  Calculates density at one point
!   !>
!   function CalculateWaveFunction_getDensityAt ( nameOfSpecie, coordinate, densityMatrix ) result( output )
!   implicit none
!   character(*), optional, intent(in):: nameOfSpecie
!   real(8) :: coordinate(3)
!   real(8) :: output
!   type(Matrix) :: densityMatrix
  
!   integer :: specieID
!   character(30) :: nameOfSpecieSelected
!   integer :: numberOfContractions
!   integer :: numberOfCartesiansOrbitals
!   integer :: totalNumberOfContractions
!   integer :: particleID
!   integer :: contractionID
!   integer :: i, j, k, u, v, g
!   real(8), allocatable :: auxVal(:)
!   real(8), allocatable :: basisSetValue(:)
  
!   character(50) :: wfnFile
!   character(50) :: arguments(20)
!   integer :: wfnUnit
  
!      nameOfSpecieSelected = "e-"
!      if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
!      specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
!      numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )
!      totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
         
!      if( allocated(basisSetValue)) deallocate(basisSetValue)
!      allocate(basisSetValue(totalNumberOfContractions))
!      output=0.0_8
!      k=1
!      do g = 1, size(MolecularSystem_instance%species(specieID)%particles)
!        do i = 1, size(MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction)
!          numberOfCartesiansOrbitals = MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction(i)%numCartesianOrbital
!          if( allocated(auxVal)) deallocate(auxVal)
!          allocate(auxVal(numberOfCartesiansOrbitals))
!          auxVal = CalculateWaveFunction_getDensityValueAt( &
!                  MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction(i), coordinate )
!          do j = 1, numberOfCartesiansOrbitals
!            basisSetValue(k) = auxVal(j) 
!            k=k+1
!          end do
!        end do
!      end do
!      do u=1,totalNumberOfContractions
!         do v=1,totalNumberOfContractions
!           output=output + densityMatrix%values(u,v)*basisSetValue(u)*basisSetValue(v)
!         end do
!      end do
  
!   end function CalculateWaveFunction_getDensityAt

!     !!**************************************************************

!     do h = 1, this%numCartesianOrbital
!        output(h) = output(h) * this%contNormalization(h)
!     end do

!   end function  CalculateWaveFunction_getDensityValueAt


! ! 	!<
! !	!! @brief  Calculates gradient density at one point
! !	!>
! !	 function CalculateWaveFunction_getGradientDensityAt ( nameOfSpecie, coordinate ) result( output )
! !		implicit none
! !		character(*), optional, intent(in):: nameOfSpecie
! !		real(8) :: coordinate(3)
! !		real(8) :: output
! !
! !		integer :: specieID
! !		character(30) :: nameOfSpecieSelected
! !		integer :: numberOfContractions
! !                integer :: numberOfCartesiansOrbitals
! !		integer :: totalNumberOfContractions
! !		integer :: particleID
! !		integer :: contractionID
! !		integer :: i, j, k, u, v
! !		real(8), allocatable :: auxVal(:)
! !		real(8), allocatable :: basisSetValue(:)
! !
! !                if ( CalculateWaveFunction_isSet() ) then
! !                   nameOfSpecieSelected = "e-"
! !                   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
! !                   specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
! !                   numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )
! !                   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
! !                   if( allocated(basisSetValue)) deallocate(basisSetValue)
! !                   allocate(basisSetValue(totalNumberOfContractions))
! !                   output=0.0_8
! !                   k=1
! !                   do i=1,numberOfContractions
! !                      particleID = MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%particleID
! !                      contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%contractionIDInParticle
! !                      numberOfCartesiansOrbitals = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital
! !                      if( allocated(auxVal)) deallocate(auxVal)
! !                      allocate(auxVal(numberOfCartesiansOrbitals))
! !                      auxVal = ContractedGaussian_getGradientAt(MolecularSystem_getContractionPtr( specieID,  numberOfContraction=i ), coordinate )
! !                      do j = 1, numberOfCartesiansOrbitals
! !
! !                         basisSetValue(k) = auxVal(j) 
! !                         k=k+1
! !                      end do
! !                   end do
! !                   do u=1,totalNumberOfContractions
! !                      do v=1,totalNumberOfContractions
! !                         output=output + CalculateWaveFunction_instance%densityMatrix(specieID)%values(u,v)*basisSetValue(u)*basisSetValue(v)
! !                      end do
! !                   end do
! !                else
! !                   call CalculateWaveFunction_exception(ERROR, "You should set the molecular system before use this function", &
! !                        "Class object Calculate Properties in the getDensityAt function" )
! !                end if
! !
! !	end function CalculateWaveFunction_getGradientDensityAt
! !
!   function CalculateWaveFunction_getOrbitalValueAt ( nameOfSpecie, orbitalNum, coordinate ) result( output )
!      implicit none
!      character(*), optional, intent(in) :: nameOfSpecie
!      integer :: orbitalNum
!      real(8) :: coordinate(3)
!      real(8) :: output

!      type(Matrix) ::  coefficientsofcombination
!      integer :: specieID
!      character(30) :: nameOfSpecieSelected
!      integer :: numberOfContractions
!      integer :: numberOfCartesiansOrbitals
!      integer :: totalNumberOfContractions
!      integer :: i, j, k, u, v, g
!      real(8), allocatable :: auxVal(:)
!      real(8), allocatable :: basisSetValue(:)
!      character(50) :: wfnFile
!      character(50) :: arguments(20)
!      integer :: wfnUnit

!      nameOfSpecieSelected = "e-"
!      if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
!      specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
!      numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )
!      totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
  
!      wfnFile = "lowdin.wfn"
!      wfnUnit = 20
  
!      !! Open file for wavefunction
!      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
  
!      arguments(2) = MolecularSystem_getNameOfSpecie(specieID)
!      arguments(1) = "COEFFICIENTS"
  
!      coefficientsofcombination = &
!      Matrix_getFromFile(unit=wfnUnit, rows= int(totalNumberOfContractions,4), &
!      columns= int(totalNumberOfContractions,4), binary=.true., arguments=arguments(1:2))
  
!      if( allocated(basisSetValue)) deallocate(basisSetValue)
!      allocate(basisSetValue(totalNumberOfContractions))
!      output=0.0_8
!      k=1
!      do g = 1, size(MolecularSystem_instance%species(specieID)%particles)
!        do i = 1, size(MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction)
!          numberOfCartesiansOrbitals = MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction(i)%numCartesianOrbital
!          if( allocated(auxVal)) deallocate(auxVal)
!          allocate(auxVal(numberOfCartesiansOrbitals))
!          auxVal = CalculateWaveFunction_getDensityValueAt( &
!               MolecularSystem_instance%species(specieID)%particles(g)%basis%contraction(i), coordinate )
!          do j = 1, numberOfCartesiansOrbitals
!            basisSetValue(k) = auxVal(j) 
!            k=k+1
!          end do
!        end do
!      end do
!      do u=1,totalNumberOfContractions
!           output=output + coefficientsofcombination%values(u,orbitalNum)*basisSetValue(u)
!      end do
!      close (wfnUnit)

!   end function CalculateWaveFunction_getOrbitalValueAt

! !
! !	 function CalculateWaveFunction_getFukuiFunctionAt ( nameOfSpecie, fukuiType ,coordinate ) result( output )
! !		implicit none
! !		character(*), optional, intent(in):: nameOfSpecie
! !		real(8) :: coordinate(3)
! !		character(*) :: fukuiType
! !		real(8) :: output
! !
! !		integer :: specieID
! !		character(30) :: nameOfSpecieSelected
! !		integer :: numberOfContractions
! !                integer :: numberOfCartesiansOrbitals
! !		integer :: totalNumberOfContractions
! !		integer :: particleID
! !		integer :: contractionID
! !		integer :: i, j, k, u, v
! !		real(8), allocatable :: auxVal(:)
! !		real(8), allocatable :: basisSetValue(:)
! !
! !                if ( CalculateWaveFunction_isSet() ) then
! !                   nameOfSpecieSelected = "e-"
! !                   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
! !                   specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
! !                   numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )
! !                   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
! !                   if( allocated(basisSetValue)) deallocate(basisSetValue)
! !                   allocate(basisSetValue(totalNumberOfContractions))
! !                   output=0.0_8
! !                   k=1
! !                   do i=1,numberOfContractions
! !                      particleID = MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%particleID
! !                      contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%contractionIDInParticle
! !                      numberOfCartesiansOrbitals = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital
! !                      if( allocated(auxVal)) deallocate(auxVal)
! !                      allocate(auxVal(numberOfCartesiansOrbitals))
! !                      auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( specieID,  numberOfContraction=i ), coordinate )
! !                      do j = 1, numberOfCartesiansOrbitals
! !
! !                         basisSetValue(k) = auxVal(j) 
! !                         k=k+1
! !                      end do
! !                   end do
! !                   if ( trim(fukuiType) .eq. "positive" ) then
! !                      do u=1,totalNumberOfContractions
! !                         do v=1,totalNumberOfContractions
! !                            ! output=output + CalculateProperties_instance%positiveFukui%values(u,v)*basisSetValue(u)*basisSetValue(v)
! !                         end do
! !                      end do
! !                   end if
! !                   if ( trim(fukuiType) .eq. "negative" ) then
! !                      do u=1,totalNumberOfContractions
! !                         do v=1,totalNumberOfContractions
! !                            ! output=output + CalculateProperties_instance%negativeFukui%values(u,v)*basisSetValue(u)*basisSetValue(v)
! !                         end do
! !                      end do
! !                   end if
! !                else
! !                   call CalculateWaveFunction_exception(ERROR, "You should set the molecular system before use this function", &
! !                        "Class object Calculate Properties in the getDensityAt function" )
! !                end if
! !
! !	end function CalculateWaveFunction_getFukuiFunctionAt


