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
  use Interface_
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

  end type TransformIntegralsD

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
  subroutine TransformIntegralsD_constructor(this)
    implicit none
    type(TransformIntegralsD) :: this

    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"

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
    real(8), allocatable, target :: ints(:)
    real(8), allocatable, target :: coeff(:, :)

    type(c_ptr) :: coeff_ptr, ints_ptr

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


    if(allocated(ints)) deallocate(ints)
    allocate(ints(sze))

    ints = 0.0_8
    call ReadIntegrals_intraSpecies(trim(nameOfSpecies), ints, CONTROL_instance%NUMBER_OF_CORES)


    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(nao, nao))

    do j = 1, nao
      do k = 1, nao
        coeff(j, k) = coefficients%values(j, k)
      end do
    end do


    ! Calling C function
    coeff_ptr = c_loc(coeff(1, 1))
    ints_ptr = c_loc(ints(1))

    call Interface_integralsTransform(coeff_ptr, ints_ptr, nao, &
                             this%p_lowerOrbital, this%p_upperOrbital, &
                             this%q_lowerOrbital, this%q_upperOrbital, &
                             this%r_lowerOrbital, this%r_upperOrbital, &
                             this%s_lowerOrbital, this%s_upperOrbital)

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    do p = 1, nao
       do q = 1, p
          do r = 1 , p
            s_max = r
            if(p == r) s_max = q 
             do s = 1,  s_max

              !! TODO: Use chunks instead.
              write(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p, q, r, s, ints(ReadIntegrals_index4(p, q, r, s))

           end do
         end do
       end do
    end do

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0.0_8 

    print *, "Non zero transformed repulsion integrals: ", sze

    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

  end subroutine TransformIntegralsD_atomicToMolecularOfOneSpecie

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
    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION == 2  ) then

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

  end subroutine TransformIntegralsD_checkMOIntegralType

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
