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

module CalculateWaveFunction_
  use Matrix_
  use Exception_
  use String_
  use MolecularSystem_
  implicit none

  !>
  !! @brief Modulo para definicion de datos de elemetos atomicos.
  !!
  !!  Este modulo define una seudoclase para manejo de datos atomicos, correspondientes
  !!  a elementos atomicos.
  !!
  !! @author Sergio A. Gonzalez Monico
  !!
  !! <b> Fecha de creacion : </b> 2008-08-05
  !!
  !! <b> Historial de modificaciones: </b>
  !!
  !!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
  !!        -# Propuso estandar de codificacion.
  !!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
  !!        -# Se adapto al estandar de codificacion propuesto.
  !!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
  !!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
  !!
  !! @see XMLParser_
  !<

  public :: &
       CalculateWaveFunction_getDensityAt, &
       CalculateWaveFunction_getOrbitalValueAt
  !		CalculateWaveFunction_getFukuiFunctionAt

contains

  !<
  !! @brief  Calculates density at one point
  !>
  subroutine CalculateWaveFunction_getDensityAt ( speciesID, coordinates, densityMatrix, output )
    implicit none
    integer :: speciesID
    type(Matrix) :: coordinates !npoints,3
    type(Matrix) :: densityMatrix
    type(Vector) :: output !npoints

    integer :: totalNumberOfContractions, gridSize
    integer :: u, v
    type(Matrix) :: basisSetValues

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    gridSize=size(coordinates%values(:,1))
    call Vector_constructor( output, gridSize, 0.0_8) 
    call Matrix_constructor( basisSetValues, int(gridSize,8), int(totalNumberOfContractions,8), 0.0_8 )

    call CalculateWaveFunction_getBasisValueAt(speciesID, coordinates, gridSize, basisSetValues)

    do u=1,totalNumberOfContractions
       output%values=output%values + densityMatrix%values(u,u)*basisSetValues%values(:,u)**2
       do v=u+1,totalNumberOfContractions
          output%values=output%values + 2.0*densityMatrix%values(u,v)*basisSetValues%values(:,u)*basisSetValues%values(:,v)
       end do
    end do

  end subroutine CalculateWaveFunction_getDensityAt

  subroutine CalculateWaveFunction_getOrbitalValueAt ( speciesID, orbitalNum, coordinates, output )
    implicit none
    integer :: speciesID
    integer :: orbitalNum
    type(Matrix) :: coordinates
    type(Vector) :: output

    type(Matrix) :: coefficientsofcombination
    type(Matrix) :: basisSetValues
    integer :: totalNumberOfContractions, gridSize
    integer :: u
    integer :: wfnUnit
    character(50) :: wfnFile
    character(50) :: arguments(20)

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    gridSize=size(coordinates%values(:,1))
    call Vector_constructor( output, gridSize, 0.0_8) 
    call Matrix_constructor( basisSetValues, int(gridSize,8), int(totalNumberOfContractions,8), 0.0_8 )

    call CalculateWaveFunction_getBasisValueAt(speciesID, coordinates, gridSize, basisSetValues)

    wfnFile = "lowdin.wfn"
    wfnUnit = 20  
    !! Open file for wavefunction
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    arguments(2) = MolecularSystem_getNameOfSpecies(speciesID)
    arguments(1) = "COEFFICIENTS"

    coefficientsofcombination = &
         Matrix_getFromFile(unit=wfnUnit, rows= int(totalNumberOfContractions,4), &
         columns= int(totalNumberOfContractions,4), binary=.true., arguments=arguments(1:2))
    close (wfnUnit)

    do u=1,totalNumberOfContractions
       output%values=output%values + coefficientsofcombination%values(u,orbitalNum)*basisSetValues%values(:,u)
    end do

  end subroutine CalculateWaveFunction_getOrbitalValueAt


  subroutine CalculateWaveFunction_getBasisValueAt ( speciesID, coordinates, gridSize, basisSetValues )
    integer :: speciesID
    type(Matrix) :: coordinates
    integer :: gridSize
    type(Matrix) :: basisSetValues

    integer :: numberOfCartesiansOrbitals
    integer :: i, j, k, g
    type(Matrix) :: auxVal

    k=1
    do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
       do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)
          numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital
          call Matrix_constructor(auxVal, int(gridSize,8), int(numberOfCartesiansOrbitals,8), 0.0_8)
          call ContractedGaussian_getValuesAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), &
               coordinates, gridSize, auxVal)
          do j = 1, numberOfCartesiansOrbitals
             basisSetValues%values(:,k) = auxVal%values(:,j) 
             k=k+1
          end do
       end do
    end do

  end subroutine CalculateWaveFunction_getBasisValueAt

end module CalculateWaveFunction_
!>
!! @brief Retorna el valor de la funcion en la coordenada especificada
!<
! function  CalculateWaveFunction_getDensityValueAt( this, coordinate ) result(output)
!   implicit none
!   type(ContractedGaussian) , intent(in) :: this
!   real(8) :: coordinate(3)
!   real(8) :: output(this%numCartesianOrbital)
!   integer :: h
!   integer :: nx, ny, nz !< indices de momento angular
!   integer :: i, j, m
!   real(8) :: auxOutput(this%numCartesianOrbital)

!   output = 0.0_8

!   do h=1, this%length
!       m = 0
!       do i = 0 , this%angularMoment
!               nx = this%angularMoment - i
!               do j = 0 , i
!                       ny = i - j
!                       nz = j
!                       m = m + 1

!                       auxOutput(m) = this%primNormalization(h,m) &
!                               * dexp(-this%orbitalExponents(h) &
!                               *((this%origin(1)-coordinate(1))**2 &
!                               + (this%origin(2)-coordinate(2))**2 &
!                               + (this%origin(3)-coordinate(3))**2) ) &
!                               * ( (coordinate(1)-this%origin(1))** nx ) &
!                               * ( (coordinate(2)-this%origin(2))** ny ) &
!                              * ( (coordinate(3)-this%origin(3))** nz )
!                end do
!         end do

!       output = output + auxOutput * this%contractionCoefficients(h)

!   end do

!   !!**************************************************************

!   do h = 1, this%numCartesianOrbital
!      output(h) = output(h) * this%contNormalization(h)
!   end do

! end function  CalculateWaveFunction_getDensityValueAt


! 	!<
!	!! @brief  Calculates gradient density at one point
!	!>
!	 function CalculateWaveFunction_getGradientDensityAt ( nameOfSpecie, coordinate ) result( output )
!		implicit none
!		character(*), optional, intent(in):: nameOfSpecie
!		real(8) :: coordinate(3)
!		real(8) :: output
!
!		integer :: speciesID
!		character(30) :: nameOfSpecieSelected
!		integer :: numberOfContractions
!                integer :: numberOfCartesiansOrbitals
!		integer :: totalNumberOfContractions
!		integer :: particleID
!		integer :: contractionID
!		integer :: i, j, k, u, v
!		real(8), allocatable :: auxVal(:)
!		real(8), allocatable :: basisSetValues(:)
!
!                if ( CalculateWaveFunction_isSet() ) then
!                   nameOfSpecieSelected = "e-"
!                   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
!                   speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
!                   numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
!                   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
!                   if( allocated(basisSetValues)) deallocate(basisSetValues)
!                   allocate(basisSetValues(totalNumberOfContractions))
!                   output=0.0_8
!                   k=1
!                   do i=1,numberOfContractions
!                      particleID = MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%particleID
!                      contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%contractionIDInParticle
!                      numberOfCartesiansOrbitals = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital
!                      if( allocated(auxVal)) deallocate(auxVal)
!                      allocate(auxVal(numberOfCartesiansOrbitals))
!                      auxVal = ContractedGaussian_getGradientAt(MolecularSystem_getContractionPtr( speciesID,  numberOfContraction=i ), coordinate )
!                      do j = 1, numberOfCartesiansOrbitals
!
!                         basisSetValues(k) = auxVal(j) 
!                         k=k+1
!                      end do
!                   end do
!                   do u=1,totalNumberOfContractions
!                      do v=1,totalNumberOfContractions
!                         output=output + CalculateWaveFunction_instance%densityMatrix(speciesID)%values(u,v)*basisSetValues(u)*basisSetValues(v)
!                      end do
!                   end do
!                else
!                   call CalculateWaveFunction_exception(ERROR, "You should set the molecular system before use this function", &
!                        "Class object Calculate Properties in the getDensityAt function" )
!                end if
!
!	end function CalculateWaveFunction_getGradientDensityAt
!

!
!	 function CalculateWaveFunction_getFukuiFunctionAt ( nameOfSpecie, fukuiType ,coordinate ) result( output )
!		implicit none
!		character(*), optional, intent(in):: nameOfSpecie
!		real(8) :: coordinate(3)
!		character(*) :: fukuiType
!		real(8) :: output
!
!		integer :: speciesID
!		character(30) :: nameOfSpecieSelected
!		integer :: numberOfContractions
!                integer :: numberOfCartesiansOrbitals
!		integer :: totalNumberOfContractions
!		integer :: particleID
!		integer :: contractionID
!		integer :: i, j, k, u, v
!		real(8), allocatable :: auxVal(:)
!		real(8), allocatable :: basisSetValues(:)
!
!                if ( CalculateWaveFunction_isSet() ) then
!                   nameOfSpecieSelected = "e-"
!                   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
!                   speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )
!                   numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
!                   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
!                   if( allocated(basisSetValues)) deallocate(basisSetValues)
!                   allocate(basisSetValues(totalNumberOfContractions))
!                   output=0.0_8
!                   k=1
!                   do i=1,numberOfContractions
!                      particleID = MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%particleID
!                      contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%contractionIDInParticle
!                      numberOfCartesiansOrbitals = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital
!                      if( allocated(auxVal)) deallocate(auxVal)
!                      allocate(auxVal(numberOfCartesiansOrbitals))
!                      auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( speciesID,  numberOfContraction=i ), coordinate )
!                      do j = 1, numberOfCartesiansOrbitals
!
!                         basisSetValues(k) = auxVal(j) 
!                         k=k+1
!                      end do
!                   end do
!                   if ( trim(fukuiType) .eq. "positive" ) then
!                      do u=1,totalNumberOfContractions
!                         do v=1,totalNumberOfContractions
!                            ! output=output + CalculateProperties_instance%positiveFukui%values(u,v)*basisSetValues(u)*basisSetValues(v)
!                         end do
!                      end do
!                   end if
!                   if ( trim(fukuiType) .eq. "negative" ) then
!                      do u=1,totalNumberOfContractions
!                         do v=1,totalNumberOfContractions
!                            ! output=output + CalculateProperties_instance%negativeFukui%values(u,v)*basisSetValues(u)*basisSetValues(v)
!                         end do
!                      end do
!                   end if
!                else
!                   call CalculateWaveFunction_exception(ERROR, "You should set the molecular system before use this function", &
!                        "Class object Calculate Properties in the getDensityAt function" )
!                end if
!
!	end function CalculateWaveFunction_getFukuiFunctionAt

