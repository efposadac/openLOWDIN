!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
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
!! @brief Module to handle direct integral calculations
!! @author Jorge Charry
!! @version 1.0
!! <b> Fecha de creacion : </b> 2015-12-28
!!
!! <b> Historial de modificaciones: </b>
!!
module DirectIntegralManager_
  use CONTROL_
  use MolecularSystem_
  use Libint2Interface_
  use RysQuadrature_
  use Matrix_
  use Stopwatch_
  !# use RysQInts_  !! Please do not remove this line

  implicit none

  public :: &
       DirectIntegralManager_getDirectIntraRepulsionIntegrals

contains

  !> 
  !! @brief Calculate Intra-species repulsion integrals directly
  !! @author J. A. Charry, 2015
  !! @version 1.0
  !! @par History
  !!    
  recursive subroutine DirectIntegralManager_getDirectIntraRepulsionIntegrals(speciesID, scheme, &
       densityMatrix, coeffMatrix, matrixA, p )
    implicit none

    integer :: speciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: coeffMatrix
    real(8), allocatable, target :: matrixA(:,:,:)
    integer :: p

    integer :: numberOfContractions

    real(8), allocatable, target :: coefficients(:,:)
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    ssize = size(coeffMatrix%values, DIM=1)
    allocate(coefficients(ssize, ssize))
    coefficients = coeffMatrix%values

    allocate(density(ssize, ssize))
    density = densityMatrix%values

    select case (trim(String_getUppercase(trim(scheme))))

       !     case("RYS")
       !        call RysQuadrature_directIntraSpecies( speciesID, "ERIS", starting, ending, int( process ) , &
       !               densityMatrix, & 
       !               twoParticlesMatrix, factor)
    case("LIBINT")
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p )

       !     ! case("CUDINT")
       !     !    call CudintInterface_computeIntraSpecies(speciesID)
    case default
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p )
    end select

    deallocate(density)
  end subroutine DirectIntegralManager_getDirectIntraRepulsionIntegrals


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly
  !! @author E. F. Posada 2016
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionIntegrals(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, couplingMatrix )
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: couplingMatrix(:,:)

    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    select case (trim(String_getUppercase(trim(scheme))))

       !case("RYS")
       ! Not implemented

    case("LIBINT")
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix)
    case default
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix)
    end select


    deallocate(density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionIntegrals

end module DirectIntegralManager_
