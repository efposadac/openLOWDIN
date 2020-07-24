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
!! @brief This class performs integrals tranformation using Ruben's implementation.
!!
!<
module TransformIntegralsD_
  use, intrinsic :: iso_c_binding
  use MolecularSystem_
  use InputManager_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
  use ReadIntegrals_
  implicit none

  type, public :: TransformIntegralsD

     character(30) :: name
     character(255) :: fileForCoefficients
     character(255) :: fileForIntegrals
     character(255) :: prefixOfFile
     integer :: numberOfContractions
     integer :: otherNumberOfContractions
     integer :: bias
     integer :: specieID
     integer :: otherSpecieID
     integer :: unidOfOutputForCoefficients
     integer :: unidOfOutputForIntegrals
     integer :: nproc
     integer :: integralStackSize

     integer :: p_lowerOrbital, p_upperOrbital
     integer :: q_lowerOrbital, q_upperOrbital
     integer :: r_lowerOrbital, r_upperOrbital
     integer :: s_lowerOrbital, s_upperOrbital

     character(50) :: partialTransform

  end type TransformIntegralsD

  interface
     subroutine TransformIntegralsD_integralsTransform_partial(coeff, ints, nao, lp, up, lq, uq, lr, ur, ls, us) &
          bind(C, name="c_integrals_transform_partial")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: coeff
       type(c_ptr), value :: ints
       integer(c_int), value :: nao
       integer(c_int), value :: lp, up, lq, uq, lr, ur, ls, us
     end subroutine TransformIntegralsD_integralsTransform_partial

     subroutine TransformIntegralsD_integralsTransform_all(coeff, ints, nao) &
          bind(C, name="c_integrals_transform_all")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: coeff
       type(c_ptr), value :: ints
       integer(c_int), value :: nao
     end subroutine TransformIntegralsD_integralsTransform_all

     subroutine TransformIntegralsD_integralsTransformInter_all(coeff, ocoeff, ints, nao, onao) &
          bind(C, name="c_integrals_transform_inter_all")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: coeff
       type(c_ptr), value :: ocoeff
       type(c_ptr), value :: ints
       integer(c_int), value :: nao
       integer(c_int), value :: onao
     end subroutine TransformIntegralsD_integralsTransformInter_all

     subroutine TransformIntegralsD_similarityTransform(nao, Op, C, Work) &
          bind(C, name='Similarity_Transform')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: nao
       type(c_ptr), value :: Op
       type(c_ptr), value :: C
       type(c_ptr), value :: Work
     end subroutine TransformIntegralsD_similarityTransform

  end interface


  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE = 0
  integer, parameter :: TWO_SPECIES = 1
  !! }

  public :: &
       TransformIntegralsD_constructor, &
       TransformIntegralsD_destructor, &
       TransformIntegralsD_show, &
       TransformIntegralsD_atomicToMolecularOfOneSpecie
  !TransformIntegralsD_readIntegralsTransformed

contains

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsD_constructor(this,partial)
    implicit none
    type(TransformIntegralsD) :: this
    character(*) :: partial
    
    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"

    this%partialTransform=trim(partial)

  end subroutine TransformIntegralsD_constructor

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsD_destructor(this)
    implicit none
    type(TransformIntegralsD) :: this

  end subroutine TransformIntegralsD_destructor

  !>
  !! @brief show
  !<
  subroutine TransformIntegralsD_show()
    implicit none
    write(*,  "(A)")  " IN-CORE TRANSFORMATION OF INTEGRALS                 " 
    write(*, "(A)")   " Implementation V. 1.0   Guerrero R. D.  2016         "
    write(*, "(A)")   " Literature:         "
    write(*, "(A)")   " Sherrill, C. David, and Henry F. Schaefer.        "
    write(*, "(A)")   " The configuration interaction method: Advances in highly correlated approaches.         "
    write(*, "(A)")   " Advances in quantum chemistry 34 (1999): 143-269.         "
    write(*, "(A)")   " ----------------------------------------------------------------------"

  end subroutine TransformIntegralsD_show

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !! 		a integrales moleculares.
  !<
  subroutine TransformIntegralsD_atomicToMolecularOfOneSpecie(this, coefficients, speciesID, nameOfSpecies)
    implicit none
    type(TransformIntegralsD) :: this
    type(Matrix) :: coefficients
    integer :: speciesID
    character(*) :: nameOfSpecies

    integer :: nao, sze
    integer :: j, k
    integer :: p, q, r, s, s_max
    integer :: reclen
    real(8) :: integral
    real(8), allocatable, target :: ints(:)
    real(8), allocatable, target :: coeff(:, :)

    type(c_ptr) :: coeff_ptr, ints_ptr
    logical :: direct = .true.

    this%prefixOfFile =""//trim(nameOfSpecies)
    this%fileForCoefficients =""//trim(nameOfSpecies)//"mo.values"
    this%specieID = speciesID

    ! Setting up intervals for transformation
    call TransformIntegralsD_checkMOIntegralType(speciesID, this)

    !! Read Integrals
    nao = MolecularSystem_getTotalNumberOfContractions(speciesID)
    sze = nao * (nao + 1) / 2
    sze = sze * (sze + 1) / 2

    this%numberOfContractions = nao

    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(nao, nao))

    do j = 1, nao
       do k = 1, nao
          coeff(j, k) = coefficients%values(j, k)
       end do
    end do

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    if (.not. direct) then

       if(allocated(ints)) deallocate(ints)
       allocate(ints(sze))
       ints = 0.0_8

       call ReadIntegrals_intraSpecies(trim(nameOfSpecies), ints)

       ! Calling C function
       coeff_ptr = c_loc(coeff(1, 1))
       ints_ptr = c_loc(ints(1))

       call TransformIntegralsD_integralsTransform_all(coeff_ptr, ints_ptr, nao)

       ! call TransformIntegralsD_integralsTransform_partial(coeff_ptr, ints_ptr, nao, &
       !                          this%p_lowerOrbital, this%p_upperOrbital, &
       !                          this%q_lowerOrbital, this%q_upperOrbital, &
       !                          this%r_lowerOrbital, this%r_upperOrbital, &
       !                          this%s_lowerOrbital, this%s_upperOrbital)

       do p = 1, nao
          do q = 1, p
             do r = 1 , p
                s_max = r
                if(p == r) s_max = q 
                do s = 1,  s_max

                   !! TODO: Use chunks instead.
                   write(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p, q, r, s, ints(ReadIntegrals_index4Intra(p, q, r, s))

                end do
             end do
          end do
       end do


    else

       call ReadIntegrals_intraSpecies(trim(nameOfSpecies), ints)

       inquire(iolength=reclen) integral
       open(unit=50,FILE=trim(nameOfSpecies)//".dints",ACCESS="direct",FORM="Unformatted",RECL=reclen, STATUS="old")

       call TransformIntegralsD_transformIntegralsIntra(50, nao, nameOfSpecies, coeff)


       do p = 1, nao
          do q = 1, p
             do r = 1 , p
                s_max = r
                if(p == r) s_max = q 
                do s = 1,  s_max

                   !! TODO: Use chunks instead.
                   read(50, rec=ReadIntegrals_index4Intra(p, q, r, s)) integral
                   write(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p, q, r, s, integral

                end do
             end do
          end do
       end do

       close(50)

    end if

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0.0_8 
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

    print *, "Non zero transformed repulsion integrals: ", sze

  end subroutine TransformIntegralsD_atomicToMolecularOfOneSpecie


  subroutine TransformIntegralsD_transformIntegralsIntra(unit, nao, nameOfSpecies, coeff)
    implicit none

    integer :: unit
    integer :: nao
    character(*) :: nameOfSpecies
    real(8), allocatable, target :: coeff(:, :)


    real(8) :: integral
    integer :: reclen
    integer :: i, j, k, l, ij, kl, index
    real(8), allocatable, target :: X(:,:)
    real(8), allocatable, target :: W(:,:)

    type(c_ptr) :: coeff_ptr, X_ptr, W_ptr

    ! First half-transformation
    inquire(iolength=reclen) integral
    open(60,FILE=trim(nameOfSpecies)//".tmp",ACCESS="direct",FORM="Unformatted",RECL=reclen, STATUS="replace")

    if(allocated(W))deallocate(W)
    allocate(W(nao, nao))

    if(allocated(X))deallocate(X)
    allocate(X(nao, nao))

    coeff_ptr = c_loc(coeff(1,1))
    W_ptr = c_loc(W(1,1))
    X_ptr = c_loc(X(1,1))

    ij = 0
    do i = 1, nao
       do j = 1, i
          ij = ij + 1
          do k = 1, nao
             do l = 1, k
                index = ReadIntegrals_index4Intra(i, j, k, l)
                read(unit, rec=index) X(k, l)
                X(l, k) = X(k, l)
             end do
          end do

          call TransformIntegralsD_similarityTransform(nao, X_ptr, coeff_ptr, W_ptr)

          kl = 0
          do k = 1, nao
             do l = 1, k
                kl = kl + 1
                index = ReadIntegrals_index2(kl, ij)
                write(60, rec=index) X(k, l);
             end do
          end do
       end do
    end do

    ! Second half-transformation
    X = 0
    kl = 0
    do k = 1, nao
       do l = 1, k
          kl = kl + 1
          ij = 0
          do i = 1, nao
             do j = 1, i
                ij = ij + 1
                index = ReadIntegrals_index2(kl, ij)
                read(60, rec=index) X(i, j)
                X(j, i) = X(i, j)
             end do
          end do

          call TransformIntegralsD_similarityTransform(nao, X_ptr, coeff_ptr, W_ptr)

          do i = 1, nao
             do j = 1, i
                index = ReadIntegrals_index4Intra(i, j, k, l)
                write(unit, rec=index) X(i, j)
             end do
          end do
       end do
    end do

    close(60)


  end subroutine TransformIntegralsD_transformIntegralsIntra

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !!    a integrales moleculares.
  subroutine TransformIntegralsD_atomicToMolecularOfTwoSpecies( this, &
       coefficients, otherCoefficients, &
       speciesID, nameOfSpecies, &
       otherSpeciesID, nameOfOtherSpecies )

    implicit none

    type(TransformIntegralsD) :: this
    type(Matrix) :: coefficients
    type(Matrix) :: otherCoefficients
    integer :: speciesID, otherSpeciesID
    character(*) :: nameOfSpecies, nameOfOtherSpecies

    real(8), allocatable, target :: ints(:)
    real(8), allocatable, target :: coeff(:, :), ocoeff(:, :)
    real(8) :: integral
    integer :: reclen
    integer :: nao, onao, sze, osze, tsze
    integer :: j, k
    integer :: p, q, r, s

    type(c_ptr) :: coeff_ptr, ocoeff_ptr, ints_ptr
    logical :: direct = .true.

    this%prefixOfFile =""//trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)
    this%fileForCoefficients =""//trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)//"mo.values"

    nao = size(coefficients%values,dim=1)
    onao = size(otherCoefficients%values,dim=1)

    this%numberOfContractions = nao 
    this%otherNumberOfContractions = onao
    this%specieID = speciesID
    this%numberOfContractions = nao

    call TransformIntegralsD_checkInterMOIntegralType(speciesID, otherSpeciesID, this)

    !! Read Integrals
    sze = nao * (nao + 1) / 2
    osze = onao * (onao + 1) / 2
    tsze = sze * osze

    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(nao, nao))

    do j = 1, nao
       do k = 1, nao
          coeff(j, k) = coefficients%values(j, k)
       end do
    end do

    if (allocated(ocoeff)) deallocate(ocoeff)
    allocate(ocoeff(onao, onao))

    do j = 1, onao
       do k = 1, onao
          ocoeff(j, k) = otherCoefficients%values(j, k)
       end do
    end do

  
    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    if (.not. direct) then

      if(allocated(ints)) deallocate(ints)
      allocate(ints(tsze))

      ints = 0.0_8
      call ReadIntegrals_interSpecies(trim(nameOfSpecies), trim(nameOfOtherSpecies), osze, ints)

      ! Calling C function
      coeff_ptr = c_loc(coeff(1, 1))
      ocoeff_ptr = c_loc(ocoeff(1, 1))
      ints_ptr = c_loc(ints(1))

      call TransformIntegralsD_integralsTransformInter_all(coeff_ptr, ocoeff_ptr, ints_ptr, nao, onao)

      do p = 1, nao
         do q = p, nao
            do r = 1 , onao
               do s = r,  onao
                  !! TODO: Use chunks instead.
                  write(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p, q, r, s, ints(ReadIntegrals_index4Inter(p, q, r, s, osze))

               end do
            end do
         end do
      end do

    else

       call ReadIntegrals_interSpecies(trim(nameOfSpecies), trim(nameOfOtherSpecies), osze, ints)

       inquire(iolength=reclen) integral
       open(unit=50,FILE=trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)//".dints",ACCESS="direct",FORM="Unformatted",RECL=reclen, STATUS="old")

       call TransformIntegralsD_transformIntegralsInter(50, nao, onao,  nameOfSpecies, nameOfOtherSpecies, coeff, ocoeff, osze)


       do p = 1, nao
          do q = p, nao
             do r = 1 , nao
                do s = r,  nao

                   !! TODO: Use chunks instead.
                   read(50, rec=ReadIntegrals_index4Inter(p, q, r, s, osze)) integral
                   write(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p, q, r, s, integral

                end do
             end do
          end do
       end do

       close(50)

    end if
    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0.0_8 

    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

  end subroutine TransformIntegralsD_atomicToMolecularOfTwoSpecies


  subroutine TransformIntegralsD_transformIntegralsInter(unit, nao, onao, nameOfSpecies, nameOfOtherSpecies, coeff, ocoeff, om)
    implicit none

    integer :: unit
    integer :: nao, onao, om
    character(*) :: nameOfSpecies, nameOfOtherSpecies
    real(8), allocatable, target :: coeff(:, :), ocoeff(:, :)


    real(8) :: integral
    integer :: reclen
    integer :: i, j, k, l, ij, kl, index
    real(8), allocatable, target :: X(:,:), OX(:,:)
    real(8), allocatable, target :: W(:,:), OW(:,:)

    type(c_ptr) :: coeff_ptr, ocoeff_ptr, X_ptr, W_ptr, OX_ptr, OW_ptr

    ! First half-transformation
    inquire(iolength=reclen) integral
    open(60,FILE="scratch_transform.tmp",ACCESS="direct",FORM="Unformatted",RECL=reclen, STATUS="replace")

    if(allocated(OW))deallocate(OW)
    allocate(OW(onao, onao))

    if(allocated(OX))deallocate(OX)
    allocate(OX(onao, onao))

    ocoeff_ptr = c_loc(ocoeff(1,1))
    OW_ptr = c_loc(OW(1,1))
    OX_ptr = c_loc(OX(1,1))

    ij = 0
    do i = 1, nao
       do j = 1, i
          ij = ij + 1
          do k = 1, onao
             do l = 1, k
                index = ReadIntegrals_index4Inter(i, j, k, l, om)
                read(unit, rec=index) OX(k, l)
                OX(l, k) = OX(k, l)
             end do
          end do

          call TransformIntegralsD_similarityTransform(onao, OX_ptr, ocoeff_ptr, OW_ptr)

          kl = 0
          do k = 1, onao
             do l = 1, k
                kl = kl + 1
                index = ReadIntegrals_index2(kl, ij)
                write(60, rec=index) OX(k, l);
             end do
          end do
       end do
    end do


    if(allocated(W))deallocate(W)
    allocate(W(nao, nao))

    if(allocated(X))deallocate(X)
    allocate(X(nao, nao))

    coeff_ptr = c_loc(coeff(1,1))
    W_ptr = c_loc(W(1,1))
    X_ptr = c_loc(X(1,1))


    ! Second half-transformation
    X = 0
    kl = 0
    do k = 1, onao
       do l = 1, k
          kl = kl + 1
          ij = 0
          do i = 1, nao
             do j = 1, i
                ij = ij + 1
                index = ReadIntegrals_index2(kl, ij)
                read(60, rec=index) X(i, j)
                X(j, i) = X(i, j)
             end do
          end do

          call TransformIntegralsD_similarityTransform(nao, X_ptr, coeff_ptr, W_ptr)

          do i = 1, nao
             do j = 1, i
                index = ReadIntegrals_index4Inter(i, j, k, l, om)
                write(unit, rec=index) X(i, j)
             end do
          end do
       end do
    end do

    close(60)

  end subroutine TransformIntegralsD_transformIntegralsInter

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsD_checkMOIntegralType(speciesID, this)
    implicit none
    integer :: speciesID
    type(TransformIntegralsD) :: this
    integer :: totalOccupation 
    integer :: totalNumberOfContractions

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )
    totalNumberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID)

    !! All orbitals. Default
    this%p_lowerOrbital = 0
    this%p_upperOrbital = totalNumberOfContractions - 1 
    this%q_lowerOrbital = 0
    this%q_upperOrbital = totalNumberOfContractions - 1 
    this%r_lowerOrbital = 0
    this%r_upperOrbital = totalNumberOfContractions - 1
    this%s_lowerOrbital = 0
    this%s_upperOrbital = totalNumberOfContractions - 1


    !! only the (ia|jb) integrals will be transformed
    if ( trim(this%partialTransform)=="MP2"  ) then

       this%p_lowerOrbital = 0
       this%p_upperOrbital = totalOccupation - 1
       this%q_lowerOrbital = totalOccupation
       this%q_upperOrbital = totalNumberOfContractions - 1
       this%r_lowerOrbital = 0
       this%r_upperOrbital = totalOccupation - 1
       this%s_lowerOrbital = totalOccupation 
       this%s_upperOrbital = totalNumberOfContractions - 1

    end if

    !!    !! only the (ia|bc) integrals will be transformed
    !!    if ( CONTROL_instance%PT_ORDER == 2 .and.  CONTROL_instance%IONIZE_MO <= totalOCcupation ) then
    !!
    !!      this%p_lowerOrbital = 1
    !!      this%p_upperOrbital = totalOccupation
    !!      this%q_lowerOrbital = 1
    !!      this%q_upperOrbital = totalNumberOfContractions
    !!      this%r_lowerOrbital = 1
    !!      this%r_upperOrbital = totalNumberOfContractions
    !!      this%s_lowerOrbital = 1
    !!      this%s_upperOrbital = totalNumberOfContractions
    !!
    !!    end if
    write(*,"(T15,A)") "Transformation boundaries "
    write(*,"(T15,A10,A6,A6)") "orbital","lower", "upper"
    write(*,"(T20,A5,I6,I6)") "p", this%p_lowerOrbital, this%p_upperOrbital
    write(*,"(T20,A5,I6,I6)") "q", this%q_lowerOrbital, this%q_upperOrbital
    write(*,"(T20,A5,I6,I6)") "r", this%r_lowerOrbital, this%r_upperOrbital
    write(*,"(T20,A5,I6,I6)") "s", this%s_lowerOrbital, this%s_upperOrbital
    print *, ""

  end subroutine TransformIntegralsD_checkMOIntegralType


  subroutine TransformIntegralsD_checkInterMOIntegralType(speciesID, otherSpeciesID, this)
    implicit none
    integer :: speciesID, otherSpeciesID
    type(TransformIntegralsD) :: this
    integer :: totalOccupation, otherTotalOccupation
    integer :: totalNumberOfContractions, otherTotalNumberOfContractions

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )
    totalNumberOfContractions =  MolecularSystem_getTotalNumberOfContractions ( speciesID )
    otherTotalOccupation = MolecularSystem_getOcupationNumber( otherSpeciesID )
    otherTotalNumberOfContractions =  MolecularSystem_getTotalNumberOfContractions ( otherSpeciesID )


    !! All orbitals. Default
    this%p_lowerOrbital = 1
    this%p_upperOrbital = totalNumberOfContractions
    this%q_lowerOrbital = 1
    this%q_upperOrbital = totalNumberOfContractions
    this%r_lowerOrbital = 1
    this%r_upperOrbital = otherTotalNumberOfContractions
    this%s_lowerOrbital = 1
    this%s_upperOrbital = otherTotalNumberOfContractions


    !! only the (ia|jb) integrals will be transformed
    if ( trim(this%partialTransform) .eq. "MP2"  ) then

       this%p_lowerOrbital = 1
       this%p_upperOrbital = totalOccupation
       this%q_lowerOrbital = totalOccupation + 1
       this%q_upperOrbital = totalNumberOfContractions
       this%r_lowerOrbital = 1
       this%r_upperOrbital = otherTotalOccupation
       this%s_lowerOrbital = otherTotalOccupation + 1
       this%s_upperOrbital = otherTotalNumberOfContractions

    end if

    write(*,"(T15,A)") "Transformation boundaries "
    write(*,"(T15,A10,A6,A6)") "orbital","lower", "upper"
    write(*,"(T20,A5,I6,I6)") "p", this%p_lowerOrbital, this%p_upperOrbital
    write(*,"(T20,A5,I6,I6)") "q", this%q_lowerOrbital, this%q_upperOrbital
    write(*,"(T20,A5,I6,I6)") "r", this%r_lowerOrbital, this%r_upperOrbital
    write(*,"(T20,A5,I6,I6)") "s", this%s_lowerOrbital, this%s_upperOrbital
    print *, ""
    
  end subroutine TransformIntegralsD_checkInterMOIntegralType


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine TransformIntegralsD_exception( typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine TransformIntegralsD_exception

end module TransformIntegralsD_
