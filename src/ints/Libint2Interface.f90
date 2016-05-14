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
  ! use Matrix_

  implicit none

  !> @brief Type for libint2/c++ library
  type, public :: LibintInterface
     type(c_ptr) :: this = C_NULL_ptr
     logical :: isInstanced = .false.
  end type LibintInterface

  !>
  !!Interface to libint_iface.cpp
  interface

    function c_LibintInterface_new () result(this) bind(C,name="LibintInterface_new")
     use, intrinsic :: iso_c_binding
     implicit none
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

    subroutine c_LibintInterface_compute2BodyMatrix(this, density, result) bind(C, name="LibintInterface_compute_2body_fock")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: this        
      type(c_ptr), value :: density
      type(c_ptr), value :: result
    end subroutine c_LibintInterface_compute2BodyMatrix

    subroutine c_LibintInterface_compute2BodyInts(this, filename) bind(C, name="LibintInterface_compute_2body_ints")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: this        
      character(c_char) :: filename(*)
    end subroutine c_LibintInterface_compute2BodyInts
  end interface

  type(LibintInterface), allocatable, dimension(:) :: Libint2Instance

contains

  !>
  !! Initialize the objects to calculate integrals using libint2 library
  !! Uses the C++ API.
  subroutine Libint2Interface_constructor(this, speciesID)
    implicit none
    type(LibintInterface) :: this
    integer :: speciesID

    type(Particle) :: particle_tmp
    type(ContractedGaussian) :: contraction_tmp
    type(c_ptr) :: origin_ptr, alpha_ptr, coeff_ptr

    real(8), target :: origin(3)
    real(8), target, allocatable :: coefficients(:), exponents(:)

    integer :: p, c

    ! Create Libint object
    this%this = c_LibintInterface_new()

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
    type(LibintInterface) :: this

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
    call c_LibintInterface_compute2BodyMatrix(Libint2Instance(speciesID)%this, density_ptr, twoBody_ptr)

  end subroutine Libint2Interface_compute2BodyIntraspecies_direct

  !>
  !! Compute 2-body integrals and store them on disk
  subroutine Libint2Interface_compute2BodyIntraspecies_disk(speciesID)
    implicit none

    integer :: speciesID

    character(50) :: filename
    integer :: nspecies

    nspecies = size(MolecularSystem_instance%species)
    if (.not. allocated(Libint2Instance)) then
       allocate(Libint2Instance(nspecies))  
    endif

    !! open file for integrals
    filename = C_CHAR_""//trim(MolecularSystem_instance%species(speciesID)%name)//".ints"//C_NULL_CHAR

    ! Initialize libint objects
    if (.not. Libint2Instance(speciesID)%isInstanced) then
      call Libint2Interface_constructor(Libint2Instance(speciesID), speciesID)
    end if

    call c_LibintInterface_init2BodyInts(Libint2Instance(speciesID)%this)
    call c_LibintInterface_compute2BodyInts(Libint2Instance(speciesID)%this, filename)

  end subroutine Libint2Interface_compute2BodyIntraspecies_disk

end module Libint2Interface_
