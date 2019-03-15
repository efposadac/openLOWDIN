!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Libint2 interface
!!
!! @author E. F. Posada
!!
!! <b> Creation data : </b> 28-04-2016
!!
!! <b> History change: </b>
!!
!!   - <tt> 28-04-16 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# Creation of the module. This module is intended to be the interface with libint2 library

module Libint2Interface_
  use, intrinsic :: iso_c_binding
  use MolecularSystem_
  use InterPotential_
  use ContractedGaussian_
  ! use Matrix_

  implicit none

  !> @brief Type for libint2/c++ library
  type, public :: Libint2Interface
     type(c_ptr) :: this = C_NULL_ptr
     logical :: isInstanced = .false.
  end type Libint2Interface

  !<
  !! libint.h structure
  type, bind(c) :: lib_int
     type(c_ptr)     :: int_stack
     type(c_ptr)     :: PrimQuartet
     real(c_double)  :: AB(3)
     real(c_double)  :: CD(3)
     type(c_ptr)     :: vrr_classes(9,9)
     type(c_ptr)     :: vrr_stack
  end type lib_int

  !<
  !! Type structure for libint2
  !>
  type, bind(c) :: libint2
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS0
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS1
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS2
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS3
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS4
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS5
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS6
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS7
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS8
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS9
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS10
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS11
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS12
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS13
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS14
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS15
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS16
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS17
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS18
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS19
     real(kind=c_double)  :: LIBINT_T_SS_EREP_SS20
     !! Prefactors for recurrence relations from Weber and Daul, Comp. Phys. Comm. 158, 1 (2004).
     real(kind=c_double)  :: LIBINT_T_SS_K0G12_SS_0
     real(kind=c_double)  :: LIBINT_T_SS_K2G12_SS_0
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS0
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS1
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS2
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS3
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS4
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS5
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS6
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS7
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS8
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS9
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS10
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS11
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS12
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS13
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS14
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS15
     real(kind=c_double)  :: LIBINT_T_SS_Km1G12_SS16
     !! LRL1991, Eq. 30, prefactor in front of (a0|c0)
     real(kind=c_double)  :: TwoPRepITR_pfac0_0_x
     real(kind=c_double)  :: TwoPRepITR_pfac0_0_y
     real(kind=c_double)  :: TwoPRepITR_pfac0_0_z
     real(kind=c_double)  :: TwoPRepITR_pfac0_1_x
     real(kind=c_double)  :: TwoPRepITR_pfac0_1_y
     real(kind=c_double)  :: TwoPRepITR_pfac0_1_z
     !! LRL1991, Eq. 30, prefactor in front of (a0|c+10)
     real(kind=c_double)  :: TwoPRepITR_pfac1_0
     real(kind=c_double)  :: TwoPRepITR_pfac1_1
     !! WD2004, Eq. 30, prefactor in front of (a0|k|c0)
     real(kind=c_double)  :: R12kG12_pfac0_0_x
     real(kind=c_double)  :: R12kG12_pfac0_0_y
     real(kind=c_double)  :: R12kG12_pfac0_0_z
     real(kind=c_double)  :: R12kG12_pfac0_1_x
     real(kind=c_double)  :: R12kG12_pfac0_1_y
     real(kind=c_double)  :: R12kG12_pfac0_1_z
     !! WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0)
     real(kind=c_double)  :: R12kG12_pfac1_0
     real(kind=c_double)  :: R12kG12_pfac1_1
     !! WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0)
     real(kind=c_double)  :: R12kG12_pfac2
     !! WD2004, Eq. 30, prefactor in front of curly brakets (excludes k)
     real(kind=c_double)  :: R12kG12_pfac3_0
     real(kind=c_double)  :: R12kG12_pfac3_1
     !! WD2004, Eq. 30, prefactor in front of (a0|k-2|c0)
     real(kind=c_double)  :: R12kG12_pfac4_0_x
     real(kind=c_double)  :: R12kG12_pfac4_0_y
     real(kind=c_double)  :: R12kG12_pfac4_0_z
     real(kind=c_double)  :: R12kG12_pfac4_1_x
     real(kind=c_double)  :: R12kG12_pfac4_1_y
     real(kind=c_double)  :: R12kG12_pfac4_1_z
     !! Exponents
     real(kind=c_double)  :: zeta_A
     real(kind=c_double)  :: zeta_B
     real(kind=c_double)  :: zeta_C
     real(kind=c_double)  :: zeta_D
     !! Squared exponents
     real(kind=c_double)  :: zeta_A_2
     real(kind=c_double)  :: zeta_B_2
     real(kind=c_double)  :: zeta_C_2
     real(kind=c_double)  :: zeta_D_2

     !! Appear in OS RR for ERIs

     !! One over 2.0*zeta
     real(kind=c_double)  :: oo2z
     !! One over 2.0*eta
     real(kind=c_double)  :: oo2e
     !! One over 2.0*(zeta+eta)
     real(kind=c_double)  :: oo2ze
     !! rho over zeta
     real(kind=c_double)  :: roz
     !! rho over eta
     real(kind=c_double)  :: roe

     !! Appear in standard OS RR for ERI and almost all other recurrence relations
     real(kind=c_double)  :: WP_x, WP_y, WP_z
     real(kind=c_double)  :: WQ_x, WQ_y, WQ_z
     real(kind=c_double)  :: PA_x, PA_y, PA_z
     real(kind=c_double)  :: QC_x, QC_y, QC_z
     real(kind=c_double)  :: AB_x, AB_y, AB_z
     real(kind=c_double)  :: CD_x, CD_y, CD_z

  end type libint2

  !>
  !!Interface to libint_iface.cpp
  interface

    function c_LibintInterface_new (stack_size, id, el) result(this) bind(C,name="LibintInterface_new")
     use, intrinsic :: iso_c_binding
     implicit none
     integer(c_int), value :: stack_size
     integer(c_int), value :: id
     logical(c_bool), value :: el
     type(c_ptr) :: this
    end function c_LibintInterface_new

     subroutine c_LibintInterface_del(this) bind(C, name="LibintInterface_del")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this
     end subroutine c_LibintInterface_del

     subroutine c_LibintInterface_addParticle(this, z, origin) bind(C, name="LibintInterface_add_particle")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this
       integer(c_int), value :: z
       type(c_ptr), value :: origin
     end subroutine c_LibintInterface_addParticle

     subroutine c_LibintInterface_addShell(this, alpha, coeff, origin, l, nprim) bind(C, name="LibintInterface_add_shell")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this
       type(c_ptr), value :: alpha
       type(c_ptr), value :: coeff
       type(c_ptr), value :: origin
       integer(c_int), value :: l
       integer(c_int), value :: nprim
     end subroutine c_LibintInterface_addShell

     subroutine c_LibintInterface_compute1BodyInts(this, integral_kind, result) bind(C, name="LibintInterface_compute_1body_ints")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this
       integer(c_int), value :: integral_kind
       type(c_ptr), value :: result

     end subroutine c_LibintInterface_compute1BodyInts

     subroutine c_LibintInterface_init2BodyInts(this) bind(C, name="LibintInterface_init_2body_ints")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this
     end subroutine c_LibintInterface_init2BodyInts

     subroutine c_LibintInterface_compute2BodyDirect(this, density, result) bind(C, name="LibintInterface_compute_2body_direct")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       type(c_ptr), value :: density
       type(c_ptr), value :: result
     end subroutine c_LibintInterface_compute2BodyDirect

     subroutine c_LibintInterface_compute2BodyDisk(this, filename, density) bind(C, name="LibintInterface_compute_2body_disk")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       character(c_char) :: filename(*)
       type(c_ptr), value :: density
     end subroutine c_LibintInterface_compute2BodyDisk

     subroutine c_LibintInterface_computeCouplingDirect(this, othis, density, result) bind(C, name="LibintInterface_compute_coupling_direct")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       type(c_ptr), value :: othis        
       type(c_ptr), value :: density
       type(c_ptr), value :: result
     end subroutine c_LibintInterface_computeCouplingDirect

     subroutine c_LibintInterface_computeCouplingDisk(this, othis, filename) bind(C, name="LibintInterface_compute_coupling_disk")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       type(c_ptr), value :: othis        
       character(c_char) :: filename(*)
     end subroutine c_LibintInterface_computeCouplingDisk

     subroutine c_LibintInterface_computeG12Disk(this, filename, coefficients, exponents, pot_size) bind(C, name="libintinterface_compute_g12_disk")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       character(c_char) :: filename(*)
       type(c_ptr), value :: coefficients
       type(c_ptr), value :: exponents
       integer(c_int), value :: pot_size
       
     end subroutine c_LibintInterface_computeG12Disk
     
     subroutine c_LibintInterface_computeG12InterDisk(this, othis, filename, coefficients, exponents, pot_size) bind(C, name="libintinterface_compute_g12inter_disk")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: this        
       type(c_ptr), value :: othis        
       character(c_char) :: filename(*)
       type(c_ptr), value :: coefficients
       type(c_ptr), value :: exponents
       integer(c_int), value :: pot_size

     end subroutine c_LibintInterface_computeG12InterDisk

  end interface

  type(Libint2Interface), allocatable, dimension(:) :: Libint2Instance

contains

  !>
  !! Initialize the objects to calculate integrals using libint2 library
  !! Uses the C++ API.
  subroutine Libint2Interface_constructor(this, speciesID)
    implicit none
    type(Libint2Interface) :: this
    integer :: speciesID

    type(Particle) :: particle_tmp
    type(ContractedGaussian) :: contraction_tmp
    type(c_ptr) :: origin_ptr, alpha_ptr, coeff_ptr

    real(8), target :: origin(3)
    real(8), target, allocatable :: coefficients(:), exponents(:)

    integer :: p, c

    ! Create Libint object
    this%this = c_LibintInterface_new(CONTROL_instance%INTEGRAL_STACK_SIZE, speciesID, &
                  MolecularSystem_instance%species(speciesID)%isElectron)

    ! Iterate over particles
    do p = 1, size(MolecularSystem_instance%species(speciesID)%particles)
       particle_tmp = MolecularSystem_instance%species(speciesID)%particles(p)

       ! Add particle to the object
       origin = particle_tmp%origin
       origin_ptr = c_loc(origin(1))

       call c_LibintInterface_addParticle(this%this, int(-particle_tmp%totalCharge), origin_ptr)

       ! Add basis-set to the object
       do c = 1, size(particle_tmp%basis%contraction)
          contraction_tmp = particle_tmp%basis%contraction(c)

          allocate(exponents(contraction_tmp%length))
          exponents = contraction_tmp%orbitalExponents 
          alpha_ptr = c_loc(exponents(1))

          allocate(coefficients(contraction_tmp%length))
          coefficients = contraction_tmp%contractionCoefficients
          coeff_ptr = c_loc(coefficients(1))

          origin = contraction_tmp%origin
          origin_ptr = c_loc(origin(1))

          call c_LibintInterface_addShell(&
               this%this, alpha_ptr, coeff_ptr, origin_ptr, contraction_tmp%angularMoment, contraction_tmp%length &
               ) 

          deallocate(exponents)
          deallocate(coefficients)

       end do
    end do

    this%isInstanced = .true.

  end subroutine Libint2Interface_constructor

  !>
  !! clear the object
  subroutine Libint2Interface_destructor(this)
    implicit none
    type(Libint2Interface) :: this

    call c_LibintInterface_del(this%this)
    this%isInstanced = .false.

  end subroutine Libint2Interface_destructor

  !>
  !! Compute all 1-body integrals and store them as a matrix
  subroutine Libint2Interface_compute1BodyInts(integral_kind)
    implicit none

    !! 1: Overlap, 2: Kinetic, 3: Attraction
    integer :: integral_kind

    type(c_ptr) :: matrix_ptr
    character(20), dimension(3) :: job = ['OVERLAP   ', 'KINETIC   ', 'ATTRACTION']
    real(8), allocatable, target :: integralsMatrix(:,:)
    integer :: nspecies
    integer :: s

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    open(unit=1000, file='libint.txt')

    do s = 1, nspecies
       ! Prepare matrix
       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = s), &
            MolecularSystem_getTotalNumberOfContractions(specieID = s)))
       matrix_ptr = c_loc(integralsMatrix(1,1))

       ! Initialize libint objects
       if (.not. Libint2Instance(s)%isInstanced) then
          call Libint2Interface_constructor(Libint2Instance(s), s)
       endif

       write(30) job(integral_kind)
       write(30) MolecularSystem_instance%species(s)%name

       call c_LibintInterface_compute1BodyInts(Libint2Instance(s)%this, integral_kind, matrix_ptr)

       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix
       write(1000, *) integralsMatrix

       ! ! Delete Libint object
       ! call Libint2Interface_destructor(Libint2Instance(s))

    end do
    close(1000)

  end subroutine Libint2Interface_compute1BodyInts

  !>
  !! Compute  2-body integrals and computes the G matrix
  subroutine Libint2Interface_compute2BodyIntraspecies_direct(speciesID, density, twoBody)
    implicit none

    integer :: speciesID
    real(8), allocatable, target :: density(:,:)
    real(8), allocatable, target :: twoBody(:,:)

    type(c_ptr) :: density_ptr
    type(c_ptr) :: twoBody_ptr

    integer :: nspecies

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    ! Prepare matrix
    if(allocated(twoBody)) deallocate(twoBody)
    allocate(twoBody(MolecularSystem_getTotalNumberOfContractions(specieID = speciesID), &
         MolecularSystem_getTotalNumberOfContractions(specieID = speciesID)))

    twoBody_ptr = c_loc(twoBody(1,1))
    density_ptr = c_loc(density(1,1))

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    endif

    call c_LibintInterface_init2BodyInts(Libint2Instance(speciesID)%this)
    call c_LibintInterface_compute2BodyDirect(Libint2Instance(speciesID)%this, density_ptr, twoBody_ptr)

  end subroutine Libint2Interface_compute2BodyIntraspecies_direct

  !>
  !! Compute 2-body integrals and store them on disk
  subroutine Libint2Interface_compute2BodyIntraspecies_disk(speciesID)
    implicit none

    integer :: speciesID
    real(8), allocatable, target :: density(:,:)

    character(50) :: filename, wfnFile
    character(30) :: labels(2)
    integer :: nspecies
    integer :: wfnUnit
    integer :: ssize
    integer :: numberOfContractions

    type(Matrix) :: aux_dens
    type(c_ptr) :: density_ptr

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    !! Get dens
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

    wfnUnit = 30
    wfnFile = "lowdin.wfn"

    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    labels(1) = "DENSITY"
    labels(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
    aux_dens = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
         columns= int(numberOfContractions,4), binary=.true., arguments=labels)

    ssize = size(aux_dens%values, DIM=1)
    allocate(density(ssize, ssize))
    density = aux_dens%values

    !! filename for integrals
    filename = C_CHAR_""//trim(MolecularSystem_instance%species(speciesID)%name)//".ints"//C_NULL_CHAR
    density_ptr = c_loc(density(1,1))

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    end if

    call c_LibintInterface_init2BodyInts(Libint2Instance(speciesID)%this)
    call c_LibintInterface_compute2BodyDisk(Libint2Instance(speciesID)%this, filename, density_ptr)

  end subroutine Libint2Interface_compute2BodyIntraspecies_disk

  !>
  !! Compute  2-body integrals and computes the G matrix
  subroutine Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, coupling)
    implicit none

    integer :: speciesID
    integer :: otherSpeciesID
    real(8), allocatable, target :: density(:,:)
    real(8), allocatable, target :: coupling(:,:)

    type(c_ptr) :: density_ptr
    type(c_ptr) :: coupling_ptr

    integer :: nspecies

    nspecies = size(MolecularSystem_instance%species)

    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    ! Prepare matrix
    if(allocated(coupling)) deallocate(coupling)
    allocate(coupling(MolecularSystem_getTotalNumberOfContractions(specieID = speciesID), &
         MolecularSystem_getTotalNumberOfContractions(specieID = speciesID)))

    coupling_ptr = c_loc(coupling(1,1))
    density_ptr = c_loc(density(1,1))

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    endif

    if (.not. Libint2Instance(otherSpeciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), otherSpeciesID)
    endif


    call c_LibintInterface_computeCouplingDirect(&
         Libint2Instance(speciesID)%this, Libint2Instance(otherSpeciesID)%this, density_ptr, coupling_ptr)

  end subroutine Libint2Interface_compute2BodyInterSpecies_direct


  !>
  !! Compute  2-body integrals and and write them to disk
  subroutine Libint2Interface_compute2BodyInterspecies_disk(speciesID, otherSpeciesID)
    implicit none

    integer :: speciesID
    integer :: otherSpeciesID

    character(50) :: filename

    integer :: nspecies

    nspecies = size(MolecularSystem_instance%species)

    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    !! filename for integrals
    filename = C_CHAR_""//trim(MolecularSystem_instance%species(speciesID)%name)//"."&
         //trim(MolecularSystem_instance%species(otherSpeciesID)%name)//".ints"//C_NULL_CHAR

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    endif

    if (.not. Libint2Instance(otherSpeciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), otherSpeciesID)
    endif


    call c_LibintInterface_computeCouplingDisk(&
         Libint2Instance(speciesID)%this, Libint2Instance(otherSpeciesID)%this, filename)

  end subroutine Libint2Interface_compute2BodyInterSpecies_disk


 !>
  !! Compute 2-body integrals and store them on disk
  subroutine Libint2Interface_computeG12Intraspecies_disk(speciesID)
    implicit none

    integer :: speciesID

    character(50) :: filename
    integer :: nspecies
    integer :: i, potID, pot_size

    real(8), allocatable, target :: coefficients(:)
    real(8), allocatable, target :: exponents(:)

    type(ContractedGaussian), pointer :: contractionG12
    type(c_ptr) :: coefficients_ptr
    type(c_ptr) :: exponents_ptr

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    !Get potential ID
    do i=1, InterPotential_instance%ssize
       if ( trim( MolecularSystem_instance%species(speciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. &
            trim( MolecularSystem_instance%species(speciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie))) ) then
          potID=i
          exit
       end if
    end do

    pot_size = size(InterPotential_instance%Potentials(potID)%gaussianComponents)
    allocate(coefficients(pot_size), exponents(pot_size))

    do i=1, pot_size
       contractionG12 => InterPotential_instance%Potentials(potID)%gaussianComponents(i)
       exponents(i) = contractionG12%orbitalExponents(1)
       coefficients(i) = contractionG12%contractionCoefficients(1)
    end do

    coefficients_ptr = c_loc(coefficients(1))
    exponents_ptr = c_loc(exponents(1))

    !! filename for integrals
    filename = C_CHAR_""//trim(MolecularSystem_instance%species(speciesID)%name)//".ints"//C_NULL_CHAR

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    end if

    call c_LibintInterface_init2BodyInts(Libint2Instance(speciesID)%this)

    call c_LibintInterface_computeG12Disk(Libint2Instance(speciesID)%this, filename, coefficients_ptr, exponents_ptr, pot_size)

  end subroutine Libint2Interface_computeG12Intraspecies_disk

  !! Compute 2-body integrals and store them on disk
  subroutine Libint2Interface_computeG12Interspecies_disk(speciesID,otherSpeciesID)
    implicit none

    integer :: speciesID, otherSpeciesID

    character(50) :: filename
    integer :: nspecies
    integer :: i, potID, pot_size

    real(8), allocatable, target :: coefficients(:)
    real(8), allocatable, target :: exponents(:)

    type(ContractedGaussian), pointer :: contractionG12
    type(c_ptr) :: coefficients_ptr
    type(c_ptr) :: exponents_ptr

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    !Get potential ID
    do i=1, InterPotential_instance%ssize
       if ( (trim(MolecularSystem_instance%species(speciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. &
            trim(MolecularSystem_instance%species(otherSpeciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie)) ) &
            ) .or. &
            (trim( MolecularSystem_instance%species(otherSpeciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. &
            trim( MolecularSystem_instance%species(speciesID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie)) ) &
            ) &
            ) then
          potID=i
          exit
       end if
    end do

    pot_size = size(InterPotential_instance%Potentials(potID)%gaussianComponents)
    allocate(coefficients(pot_size), exponents(pot_size))

    do i=1, pot_size
       contractionG12 => InterPotential_instance%Potentials(potID)%gaussianComponents(i)
       exponents(i) = contractionG12%orbitalExponents(1)
       coefficients(i) = contractionG12%contractionCoefficients(1)
    end do

    coefficients_ptr = c_loc(coefficients(1))
    exponents_ptr = c_loc(exponents(1))

    !! filename for integrals
    filename = C_CHAR_""//trim(MolecularSystem_instance%species(speciesID)%name)//"."&
         //trim(MolecularSystem_instance%species(otherSpeciesID)%name)//".ints"//C_NULL_CHAR

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    endif

    if (.not. Libint2Instance(otherSpeciesID)%isInstanced) then
       call Libint2Interface_constructor(Libint2Instance(otherSpeciesID), otherSpeciesID)
    endif

    
    call c_LibintInterface_computeG12InterDisk(&
         Libint2Instance(speciesID)%this, Libint2Instance(otherSpeciesID)%this, filename, coefficients_ptr, exponents_ptr, pot_size)

    
  end subroutine Libint2Interface_computeG12Interspecies_disk



end module Libint2Interface_
