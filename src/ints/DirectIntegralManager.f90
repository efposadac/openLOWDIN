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
!! <b> Historial de modificaciones: </b> Merge DirectIntegralManagers from different programs Felix 2022-4
!!
module DirectIntegralManager_
  use CONTROL_  
  use MolecularSystem_
  use OverlapIntegrals_
  use ThreeCOverlapIntegrals_
  use AttractionIntegrals_
  use KineticIntegrals_
  use MomentIntegrals_
  use HarmonicIntegrals_
  use Libint2Interface_
  use RysQuadrature_
  use Matrix_
  use Stopwatch_
  use GTFPotential_
  use String_
  !# use RysQInts_  !! Please do not remove this line

  implicit none

  public :: &
       DirectIntegralManager_getDirectIntraRepulsionMatrix, &
       DirectIntegralManager_getDirectInterRepulsionMatrix, &
       DirectIntegralManager_getDirectIntraRepulsionG12Matrix, &
       DirectIntegralManager_getDirectInterRepulsionG12Matrix, &
       DirectIntegralManager_getDirectIntraRepulsionFirstQuarter, &
       DirectIntegralManager_getDirectInterRepulsionFirstQuarter, &
       DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix, &
       DirectIntegralManager_constructor, &
       DirectIntegralManager_destructor, &
       DirectIntegralManager_getOverlapIntegrals, &
       DirectIntegralManager_getKineticIntegrals, &
       DirectIntegralManager_getAttractionIntegrals, &
       DirectIntegralManager_getMomentIntegrals, &
       DirectIntegralManager_getElectricFieldIntegrals, &
       DirectIntegralManager_getHarmonicIntegrals, &
       DirectIntegralManager_getExternalPotentialIntegrals, &
       DirectIntegralManager_getDirectIntraRepulsionIntegralsAll, &
       DirectIntegralManager_getDirectInterRepulsionIntegralsAll
contains

  !> 
  !! @brief Calculate Intra-species repulsion integrals directly and use it to obtain two particles matrix
  !! @author J. A. Charry, 2015
  !! @version 1.0
  !! @par History
  !!    
  recursive subroutine DirectIntegralManager_getDirectIntraRepulsionMatrix(speciesID, scheme, densityMatrix, twoParticlesMatrix, factor, system, Libint2Local )
    implicit none

    integer :: speciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: twoParticlesMatrix(:,:)
    real(8) :: factor
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2Local(:)

    type(MolecularSystem), pointer :: molSys
    ! integer :: numberOfContractions
    ! integer(8) :: integralsByProcess
    ! integer(8) :: nprocess
    ! integer(8) :: process
    ! integer(8) :: starting
    ! integer(8) :: ending
    real(8), allocatable, target :: density(:,:)
    integer(8) :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    ssize = size(densityMatrix%values, DIM=1)
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    !     numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)

    !     ssize = (numberOfContractions * (numberOfContractions + 1))/2
    !     ssize = (ssize * (ssize + 1))/2

    !     integralsByProcess = ceiling( real(ssize,8)/real(nprocess,8) )

    !     ending = process * integralsByProcess
    !     starting = ending - integralsByProcess + 1

    !     if( starting > ssize ) return

    !     if( ending > ssize ) ending = ssize

    !! Calculate integrals
    if (trim(String_getUppercase(trim(scheme))) .ne. "LIBINT") STOP "The integral method selected has not been implemented"
       !     case("RYS")
       !        call RysQuadrature_directIntraSpecies( speciesID, "ERIS", starting, ending, int( process ) , &
       !               densityMatrix, & 
       !               twoParticlesMatrix, factor)
       !     case("CUDINT")
       !        call CudintInterface_computeIntraSpecies(speciesID)
    if( present(Libint2Local) ) then
       if (.not. Libint2Local(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(speciesID), molSys, speciesID)
       call Libint2Interface_compute2BodyIntraspecies_direct(speciesID, density, twoParticlesMatrix, factor, molSys, Libint2Local(speciesID) )
    else
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molSys%species)))
       if (.not. Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), molSys, speciesID)
       call Libint2Interface_compute2BodyIntraspecies_direct(speciesID, density, twoParticlesMatrix, factor, molSys, Libint2Instance(speciesID) )
    end if

    
    deallocate(density)
  end subroutine DirectIntegralManager_getDirectIntraRepulsionMatrix


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly and use it to obtain coupling matrix
  !! @author E. F. Posada 2016
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionMatrix(speciesID, OtherSpeciesID, scheme, densityMatrix, couplingMatrix, system, Libint2Local )
    implicit none
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: couplingMatrix(:,:)
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2Local(:)

    type(MolecularSystem), pointer :: molSys
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if
    
    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    if (trim(String_getUppercase(trim(scheme))) .ne. "LIBINT") STOP "The integral method selected has not been implemented"

    if( present(Libint2Local) ) then
       if (.not. Libint2Local(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(speciesID), molSys, speciesID)
       if (.not. Libint2Local(otherSpeciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(otherSpeciesID), molSys, otherSpeciesID)
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix, molSys, Libint2Local(speciesID), Libint2Local(otherSpeciesID))
    else
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molSys%species)))
       if (.not. Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), molSys, speciesID)
       if (.not. Libint2Instance(otherSpeciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), molSys, otherSpeciesID)
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix, molSys, Libint2Instance(speciesID), Libint2Instance(otherSpeciesID))
    end if

    deallocate(density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionMatrix

  !> 
  !! @brief Calculate Intra-species G12 repulsion integrals directly and use it to obtain two particles matrix
  !! @author F. Moncada 2023
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectIntraRepulsionG12Matrix(speciesID, densityMatrix, twoParticlesMatrix, factor, system, Libint2Local )
    implicit none

    integer :: speciesID
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: twoParticlesMatrix(:,:)
    real(8) :: factor
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2Local(:)

    type(MolecularSystem), pointer :: molSys
    real(8), allocatable, target :: density(:,:)
    integer(8) :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if
    
    ssize = size(densityMatrix%values, DIM=1)
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    if( present(Libint2Local) ) then
       if (.not. Libint2Local(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(speciesID), molSys, speciesID)
       call Libint2Interface_computeG12Intraspecies_direct(speciesID, density, twoParticlesMatrix, factor, molSys, Libint2Local(speciesID) )
    else
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molSys%species)))
       if (.not. Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), molSys, speciesID)
       call Libint2Interface_computeG12Intraspecies_direct(speciesID, density, twoParticlesMatrix, factor, molSys, Libint2Instance(speciesID) )
    end if
    
    deallocate(density)
  end subroutine DirectIntegralManager_getDirectIntraRepulsionG12Matrix


  !> 
  !! @brief Calculate Inter-species G12 repulsion integrals directly and use it to obtain coupling matrix
  !! @author F. Moncada 2023
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionG12Matrix(speciesID, OtherSpeciesID, densityMatrix, couplingMatrix, system, Libint2Local)
    implicit none
    integer :: speciesID
    integer :: otherSpeciesID
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: couplingMatrix(:,:)
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2Local(:)

    type(MolecularSystem), pointer :: molSys
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if
    
    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    ! Initialize libint objects
    if( present(Libint2Local)) then
       if (.not. Libint2Local(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(speciesID), molSys, speciesID)
       if (.not. Libint2Local(otherSpeciesID)%isInstanced) call Libint2Interface_constructor(Libint2Local(otherSpeciesID), molSys, otherSpeciesID)
       call Libint2Interface_computeG12Interspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix, molSys, Libint2Local(speciesID), Libint2Local(otherSpeciesID))
    else
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molSys%species)))
       if (.not. Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), molSys, speciesID)
       if (.not. Libint2Instance(otherSpeciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), molSys, otherSpeciesID)
       call Libint2Interface_computeG12Interspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix, molSys, Libint2Instance(speciesID), Libint2Instance(otherSpeciesID))
    end if

    deallocate(density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionG12Matrix

  !> 
  !! @brief Calculate Intra-species repulsion integrals directly and transform the first index from atomic to molecular orbitals
  !! @author J. A. Charry, 2015
  !! @version 1.0
  !! @par History
  !!    
  recursive subroutine DirectIntegralManager_getDirectIntraRepulsionFirstQuarter(speciesID, scheme, &
       densityMatrix, coeffMatrix, matrixA, p, system, Libint2LocalForSpecies )
    implicit none

    integer :: speciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: coeffMatrix
    real(8), allocatable, target :: matrixA(:,:,:)
    integer :: p
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2LocalForSpecies

    type(MolecularSystem), pointer :: molSys

    real(8), allocatable, target :: coefficients(:,:)
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    ssize = size(coeffMatrix%values, DIM=1)
    allocate(coefficients(ssize, ssize))
    coefficients = coeffMatrix%values

    allocate(density(ssize, ssize))
    density = densityMatrix%values

    if( present(Libint2LocalForSpecies) ) then
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p, molSys, Libint2LocalForSpecies )
    else
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p, molSys, Libint2Instance(speciesID) )
    end if
    ! select case (trim(String_getUppercase(trim(scheme))))

    !     case("RYS")
    !        call RysQuadrature_directIntraSpecies( speciesID, "ERIS", starting, ending, int( process ) , &
    !               densityMatrix, & 
    !               twoParticlesMatrix, factor)
    ! case("LIBINT")
    !    call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p, molSys )

    !     ! case("CUDINT")
    !     !    call CudintInterface_computeIntraSpecies(speciesID)
    ! case default
    ! end select

    deallocate(coefficients,density)

  end subroutine DirectIntegralManager_getDirectIntraRepulsionFirstQuarter


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly and transform the first index from atomic to molecular orbitals
  !! @author E. F. Posada 2016
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionFirstQuarter(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, coeffMatrix, couplingMatrix, p, system, Libint2LocalForSpecies, Libint2LocalForOtherSpecies )
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: coeffMatrix
    real(8), allocatable, target :: couplingMatrix(:,:,:)
    integer :: p
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2LocalForSpecies, Libint2LocalForOtherSpecies

    type(MolecularSystem), pointer :: molSys
    real(8), allocatable, target :: coefficients(:,:)
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    ssize = size(coeffMatrix%values, DIM=1)

    allocate(coefficients(ssize, ssize))
    coefficients = coeffMatrix%values

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    ! select case (trim(String_getUppercase(trim(scheme))))

    !case("RYS")
    ! Not implemented

    ! case("LIBINT")
    !    call Libint2Interface_compute2BodyInterspecies_direct_IT(speciesID, otherSpeciesID, density, coefficients, couplingMatrix, p, molSys)
    ! case default
    if( present(Libint2LocalForSpecies) .and. present(Libint2LocalForOtherSpecies) ) then
       call Libint2Interface_compute2BodyInterspecies_direct_IT(speciesID, otherSpeciesID, density, coefficients, couplingMatrix, p, &
            molSys, Libint2LocalForSpecies, Libint2LocalForOtherSpecies)
    else
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molSys%species)))
       if (.not. Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), molSys, speciesID)
       if (.not. Libint2Instance(otherSpeciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), molSys, otherSpeciesID)
       call Libint2Interface_compute2BodyInterspecies_direct_IT(speciesID, otherSpeciesID, density, coefficients, couplingMatrix, p, &
            molSys, Libint2Instance(speciesID), Libint2Instance(otherSpeciesID))
    end if
    ! end select

    deallocate(coefficients,density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionFirstQuarter

  subroutine DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, otherdensityMatrix, coupling )
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: otherdensityMatrix
    real(8), allocatable, target :: coupling(:,:)
    real(8), allocatable, target :: density(:,:)
    real(8), allocatable, target :: otherdensity(:,:)
    integer :: ssize

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    if(allocated(density))deallocate(density)
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    ssize = size(otherdensityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    if(allocated(otherdensity))deallocate(otherdensity)
    allocate(otherdensity(ssize, ssize))
    otherdensity = otherdensityMatrix%values


    select case (trim(String_getUppercase(trim(scheme))))

       !case("RYS")
       ! Not implemented

    case("LIBINT")
       call Libint2Interface_compute2BodyAlphaBeta_direct(speciesID, otherSpeciesID, density, otherdensity, coupling)
    case default
       call Libint2Interface_compute2BodyAlphaBeta_direct(speciesID, otherSpeciesID, density, otherdensity, coupling)
    end select


    deallocate(density)

  end subroutine DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix

  !> 
  !! @brief Create libint interface objects, to be computed simultaneously
  !! @author F. Moncada 2022
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_constructor(Libint2LocalInstance, molSys)
    implicit none
    type(Libint2Interface) :: Libint2LocalInstance(:)
    type(molecularSystem) :: molSys
    integer speciesID
    logical(1) :: parallel

    ! if (allocated(Libint2LocalInstance)) then
    !    call DirectIntegralManager_destructor(Libint2LocalInstance)
    !    deallocate(Libint2LocalInstance,stat=allocation)
    !    if(allocation .gt. 0) STOP "libint failed deallocation in DirectIntegralManager_constructor"
    ! end if

    ! allocate(Libint2LocalInstance(molSys%numberOfQuantumSpecies),stat=allocation)
    ! if(allocation .gt. 0) STOP "libint failed allocation in DirectIntegralManager_constructor"

    parallel=.false.
    do speciesID=1, size(Libint2LocalInstance(:))
       call Libint2Interface_constructor(Libint2LocalInstance(speciesID), molSys, speciesID, parallel)
       call c_LibintInterface_init2BodyInts(Libint2LocalInstance(speciesID)%this)
    end do

  end subroutine DirectIntegralManager_constructor

  !> 
  !! @brief Destroy libint interface objects, to update the molecular system  
  !! @author F. Moncada 2022
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_destructor(Libint2LocalInstance)
    implicit none
    type(Libint2Interface), optional :: Libint2LocalInstance(:)
    integer speciesID

    if( present(Libint2LocalInstance) ) then
       do speciesID=1, size(Libint2LocalInstance(:))
          if(Libint2LocalInstance(speciesID)%isInstanced) call Libint2Interface_destructor(Libint2LocalInstance(speciesID)) 
       end do
    else
       do speciesID=1, size(Libint2Instance(:))
          if(Libint2Instance(speciesID)%isInstanced) call Libint2Interface_destructor(Libint2Instance(speciesID))
       end do
    end if

  end subroutine DirectIntegralManager_destructor

  !> 
  !! @brief Calculate overlap integrals and return them as matrices
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  subroutine DirectIntegralManager_getOverlapIntegrals(molSystem,speciesID,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    Type(Matrix), intent(out) :: integralsMatrix

    integer :: g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    ! real(8) :: maxOverlap

    !!Overlap Integrals for one species    

    if(allocated(labels)) deallocate(labels)
    allocate(labels(molSystem%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(molSystem%species(speciesID))

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)

          hh = h

          ii = ii + 1
          jj = ii - 1

          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !! allocating memory Integrals for shell
                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                integralValue = 0.0_8

                !! Calculating integrals for shell
                call OverlapIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j), integralValue)

                !! saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1

                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)

                      ! if( (integralValue(m) .gt. maxOverlap) .and. (l .ne. k) ) maxOverlap=integralValue(m)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do

    !!Depuration block
    ! print*, "Overlap Matrix for species: ", speciesID
    ! call Matrix_show(integralsMatrix)

  end subroutine DirectIntegralManager_getOverlapIntegrals

  !> 
  !! @brief Calculate kinetic energy integrals and return them as a matrix
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  subroutine DirectIntegralManager_getKineticIntegrals(molSystem,speciesID,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    Type(Matrix), intent(out) :: integralsMatrix

    integer :: g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)

    !!Kinetic Integrals for one species

    if(allocated(labels)) deallocate(labels)
    allocate(labels(molSystem%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(molSystem%species(speciesID))

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)

          hh = h

          ii = ii + 1
          jj = ii - 1

          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !! allocating memory Integrals for shell
                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                !!Calculating integrals for shell
                call KineticIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j), integralValue)

                !!saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1

                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do

    !!Depuration block
    ! print*, "Kinetic Matrix for specie: ", speciesID
    ! call Matrix_show(integralsMatrix)

  end subroutine DirectIntegralManager_getKineticIntegrals


  !> 
  !! @brief Calculate point charge - quantum particle attraction integrals, return them as a matrix
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: reads point charge information from lowdin.bas file
  subroutine DirectIntegralManager_getAttractionIntegrals(molSystem,speciesID,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    Type(Matrix), intent(out) :: integralsMatrix
    !
    integer :: p, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer :: numberOfPointCharges
    ! integer :: numberOfSurfaceSegments
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    type(pointCharge), allocatable :: point(:)
    ! character(20) :: colNum
    character(50) :: symbolOfSpecies

    symbolOfSpecies = molSystem%species(speciesID)%symbol

    numberOfPointCharges = molSystem%numberOfPointCharges

    !! Allocating memory for point charges objects
    if(allocated(point)) deallocate(point)
    allocate(point(0:numberOfPointCharges - 1))

    do p = 0, numberOfPointCharges - 1
       point(p)%charge = molSystem%pointCharges(p+1)%charge
       point(p)%x  = molSystem%pointCharges(p+1)%origin(1)
       point(p)%y  = molSystem%pointCharges(p+1)%origin(2)
       point(p)%z  = molSystem%pointCharges(p+1)%origin(3)
       point(p)%qdoCenterOf  = molSystem%pointCharges(p+1)%qdoCenterOf
    end do

    if(allocated(labels)) deallocate(labels)
    allocate(labels(molSystem%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(molSystem%species(speciesID))

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)

          hh = h
          ii = ii + 1
          jj = ii - 1

          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !! allocating memory Integrals for shell

                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))


                !!Calculating integrals for shell
                call AttractionIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j), point, numberOfPointCharges, integralValue, speciesID, symbolOfSpecies)

                !!saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1


                      ! write(*,*)"lowdin integrals:f, m,k,l, integral value",f,m,k,l,integralValue(m)
                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do

    !!Depuration block
    ! print*, "Attraction  Matrix for specie: ", speciesID
    ! call Matrix_show(integralsMatrix)

  end subroutine DirectIntegralManager_getAttractionIntegrals

  !> 
  !! @brief Calculate moment integrals
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  subroutine DirectIntegralManager_getMomentIntegrals(molSystem,speciesID,component,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    integer :: component  !! components 1=x, 2=y, 3=z    
    Type(Matrix), intent(out) :: integralsMatrix

    integer :: g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)

    !!Moment Integrals for one species, one component

    if(allocated(labels)) deallocate(labels)
    allocate(labels(molSystem%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(molSystem%species(speciesID))

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    !if(component.gt.3) return !????
    
    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)

          hh = h

          ii = ii + 1
          jj = ii - 1

          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !! allocating memory Integrals for shell
                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                !!Calculating integrals for shell
                call MomentIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j), [0.0_8, 0.0_8, 0.0_8], component, integralValue)

                !!saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1

                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do

    !!Depuration block
    ! print*, "Moment Matrix for specie: ", speciesID, " and component ", component
    ! call Matrix_show(integralsMatrix)

  end subroutine DirectIntegralManager_getMomentIntegrals

  ! !>
  ! !! @brief Return real labels for integrals
  ! !! @autor E. F. Posada, 2011
  ! !! @version 1.0
  function DirectIntegralManager_getLabels(specieSelected) result(labelsOfContractions)
    implicit none

    type(species) :: specieSelected
    integer:: labelsOfContractions(specieSelected%basisSetSize)

    integer:: auxLabelsOfContractions
    integer:: i, j, k

    auxLabelsOfContractions = 1

    ! write(*,*)"labels data from integral manager"

    k = 0
    ! write(*,*)size(specieSelected%particles)
    do i = 1, size(specieSelected%particles)

       do j = 1, size(specieSelected%particles(i)%basis%contraction)

          k = k + 1

          !!position for cartesian contractions
          labelsOfContractions(k) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital
          ! write(*,*)specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital

       end do
    end do


  end function DirectIntegralManager_getLabels

  subroutine DirectIntegralManager_getExternalPotentialIntegrals(molSystem,speciesID,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    Type(Matrix) :: integralsMatrix

    integer :: g, h, i
    integer :: j, k, l, m, r
    integer :: p, potID
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:), auxintegralValue(:)

    !!Overlap Integrals for one species    
    potID = 0
    
    do i= 1, ExternalPotential_instance%ssize
       !if( trim(potential(i)%specie)==trim(interactNameSelected) ) then ! This does not work for UHF
       ! if ( String_findSubstring(trim( molSystem%species(speciesID)%name  ), &
       !      trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%species)))) == 1 ) then
       if ( trim( molSystem%species(speciesID)%symbol) == trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%species))) ) then
          potID=i
          exit
       end if
    end do

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    ! Safety check, If there is no potential, do not try to fill the matrix
    if(potID .eq. 0) return

    if(allocated(labels)) deallocate(labels)
    allocate(labels(molSystem%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(molSystem%species(speciesID))
    
    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)

          hh = h

          ii = ii + 1
          jj = ii - 1

          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !! allocating memory Integrals for shell
                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                if(allocated(auxintegralValue)) deallocate(auxintegralValue)
                allocate(auxintegralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                integralValue = 0.0_8
                do p = 1, ExternalPotential_instance%potentials(potID)%numOfComponents
                   auxintegralValue = 0.0_8

                   do r = 1, ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%numcartesianOrbital
                      ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%primNormalization( &
                           1:ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%length,r) = 1

                      ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%contNormalization(r) = 1
                   end do

                   !! Calculating integrals for shell
                   call ThreeCOverlapIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                        molSystem%species(speciesID)%particles(i)%basis%contraction(j), &
                        ExternalPotential_instance%potentials(potID)%gaussianComponents(p), auxintegralValue)
                   integralValue = integralValue + auxintegralValue

                end do !! potential
                !! saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1

                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do

  end subroutine DirectIntegralManager_getExternalPotentialIntegrals
 
  !> 
  !! @brief Calculate Intra-species repulsion integrals directly 
  !! @author F.M. 3/2023
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectIntraRepulsionIntegralsAll(speciesID, &
       densityMatrix, ints, system, Libint2LocalForSpecies )
    implicit none

    integer :: speciesID
    type(matrix) :: densityMatrix
    real(8), target :: ints(:,:,:,:)
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2LocalForSpecies

    type(MolecularSystem), pointer :: molSys
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    ssize = size(densityMatrix%values, DIM=1)

    allocate(density(ssize, ssize))
    density = densityMatrix%values

    if( present(Libint2LocalForSpecies) ) then
       call Libint2Interface_compute2BodyIntraspecies_direct_all(speciesID, density, ints, molSys, Libint2LocalForSpecies )
    else
       call Libint2Interface_compute2BodyIntraspecies_direct_all(speciesID, density, ints, molSys, Libint2Instance(speciesID) )
    end if

    deallocate(density)

  end subroutine DirectIntegralManager_getDirectIntraRepulsionIntegralsAll


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly
  !! @author F.M. 3/2023
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionIntegralsAll(speciesID, OtherSpeciesID, &
       densityMatrix, ints, system, Libint2LocalForSpecies, Libint2LocalForOtherSpecies )
    integer :: speciesID
    integer :: otherSpeciesID
    type(matrix) :: densityMatrix
    real(8), target :: ints(:,:,:,:)
    type(MolecularSystem), optional, target :: system
    type(Libint2Interface), optional :: Libint2LocalForSpecies, Libint2LocalForOtherSpecies

    type(MolecularSystem), pointer :: molSys

    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    if( present(system) ) then
       molSys=>system
    else
       molSys=>MolecularSystem_instance
    end if

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    if( present(Libint2LocalForSpecies) .and. present(Libint2LocalForOtherSpecies) ) then
       call Libint2Interface_compute2BodyInterspecies_direct_all(speciesID, otherSpeciesID, density, ints, &
            molSys, Libint2LocalForSpecies, Libint2LocalForOtherSpecies)
    else
       call Libint2Interface_compute2BodyInterspecies_direct_all(speciesID, otherSpeciesID, density, ints, &
            molSys, Libint2Instance(speciesID), Libint2Instance(otherSpeciesID))
    end if

    deallocate(density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionIntegralsAll

  !! @author F.M. 3/2025
  subroutine DirectIntegralManager_getElectricFieldIntegrals(molSystem,speciesID,integralsMatrices)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    Type(Matrix) :: integralsMatrices(1:3)

    integer :: component

    do component=1, 3
       call DirectIntegralManager_getMomentIntegrals(molSystem,speciesID,component,integralsMatrices(component))
    end do
    
  end subroutine DirectIntegralManager_getElectricFieldIntegrals

  !! @author F.M. 3/2025
  subroutine DirectIntegralManager_getHarmonicIntegrals(molSystem,speciesID,origin,integralsMatrix)
    implicit none
    type(MolecularSystem) :: molSystem
    integer :: speciesID
    real(8) :: origin(3)
    Type(Matrix) :: integralsMatrix

    integer :: g, h, i
    integer :: j, k, l, m, r
    integer :: p, potID
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:), auxintegralValue(:)

    if(allocated(labels)) deallocate(labels)
    allocate(labels(MolecularSystem_instance%species(speciesID)%basisSetSize))
    labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(speciesID))

    call Matrix_constructor(integralsMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), &
         int(MolecularSystem_getTotalNumberOfContractions(speciesID,molSystem),8), 0.0_8)

    ii = 0
    do g = 1, size(molSystem%species(speciesID)%particles)
       do h = 1, size(molSystem%species(speciesID)%particles(g)%basis%contraction)
          
          hh = h
          
          ii = ii + 1
          jj = ii - 1
          
          do i = g, size(molSystem%species(speciesID)%particles)
             do j = hh, size(molSystem%species(speciesID)%particles(i)%basis%contraction)
                
                jj = jj + 1
                
                !! allocating memory Integrals for shell
                if(allocated(integralValue)) deallocate(integralValue)
                allocate(integralValue(molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                !!Calculating integrals for shell
                call HarmonicIntegrals_computeShell( molSystem%species(speciesID)%particles(g)%basis%contraction(h), &
                     molSystem%species(speciesID)%particles(i)%basis%contraction(j), integralValue, origin)
                
                !!saving integrals on Matrix
                m = 0
                do k = labels(ii), labels(ii) + (molSystem%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (molSystem%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1
                      integralsMatrix%values(k, l) = integralValue(m)
                      integralsMatrix%values(l, k) = integralsMatrix%values(k, l)
                      
                   end do
                end do
                
             end do
             hh = 1
          end do
          
       end do
    end do
    
  end subroutine DirectIntegralManager_getHarmonicIntegrals

end module DirectIntegralManager_
