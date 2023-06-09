!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!> @brief This module contains all the routines to handle external potentials
!! @author E. F. Posada (efposadac@unal.edu.co)
!! @version 2.0
!! <b> Creation data : </b> 06-08-10
!!
!! <b> History change: </b>
!!
!!   - <tt> 06-08-10 </tt>:  Sergio A. Gonzalez ( sergmonic@gmail.com )
!!        -# Creacioon del modulo y metodos basicos
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
module ExternalPotential_
  use ContractedGaussian_
  use String_
  use Matrix_
  use Units_
  use Exception_
  implicit none

  type, public :: ExternalPot
    character(20) :: name
    character(50) :: specie
    character(50) :: ttype
    character(50) :: units
    integer :: numOfComponents
    integer :: iter
    type(ContractedGaussian), allocatable :: gaussianComponents(:)
  end type

  type, public :: ExternalPotential
    integer :: ssize
    type(ExternalPot), allocatable :: potentials(:)
    logical :: isInstanced
  end type

  type(ExternalPotential), public, target :: ExternalPotential_instance

contains


  !>
  !! @brief Constructor by default
  !! @param this
  subroutine ExternalPotential_constructor(numberOfPotentials)
    implicit none

    integer :: numberOfPotentials

    ExternalPotential_instance%ssize = numberOfPotentials
    allocate(ExternalPotential_instance%potentials(numberOfPotentials))
    ExternalPotential_instance%isInstanced = .true.

  end subroutine ExternalPotential_constructor

  !>
  !! @brief Destroys the object
  !! @param this
  subroutine ExternalPotential_destructor()
      implicit none
      
      integer :: i

      do i = 1, ExternalPotential_instance%ssize
        if (allocated(ExternalPotential_instance%potentials(i)%gaussianComponents)) deallocate(ExternalPotential_instance%potentials(i)%gaussianComponents)
      end do

      if (allocated(ExternalPotential_instance%potentials) ) deallocate(ExternalPotential_instance%potentials)
      ExternalPotential_instance%isInstanced=.false.

  end subroutine ExternalPotential_destructor

  !>
  !! @brief Shows information of the object
  !! @param this
  subroutine ExternalPotential_show()
      implicit none
      type(ExternalPot), pointer :: this
      integer ::  potId, i

      do potId = 1, ExternalPotential_instance%ssize
        this => ExternalPotential_instance%potentials(potId)

        print *,""
        print *,"======="
        print *, "External Potential for ", trim(this%specie), " : ", trim(this%name)
        print *, "Type : ", trim(this%ttype)
        write(6,"(T10,A20,A10,A10,A10,A10,A20)") "Exponent", "l", "R_x", "R_y", "R_z", "Factor"

        do i=1,this%numOfComponents
          write(6,"(T10,F20.10,I10,F10.5,F10.5,F10.5,F20.10)") &
          this%gaussianComponents(i)%orbitalExponents(1), &
          this%gaussianComponents(i)%angularMoment, this%gaussianComponents(i)%origin(:), &
          this%gaussianComponents(i)%contractionCoefficients(1)
        end do
      end do

  end subroutine ExternalPotential_show

  !>
  !! @brief loads information from the input file
  !! @param this
  !! @author E. F. Posada, 2013
  subroutine ExternalPotential_load(potId, name, species)
    implicit none
    integer :: potId
    character(*) :: name
    character(*) :: species
    
    type(ExternalPot), pointer :: this
    integer :: status, i, j
    character(150) :: fileName
    character(20) :: token
    character(10) :: symbol
    logical :: existFile, found

    this => ExternalPotential_instance%potentials(potId)

    this%name= trim(name)
    this%specie= trim(species)
    this%ttype=""
    this%units="bohr"
    this%numOfComponents=0
    this%iter=1

    fileName = trim( trim( CONTROL_instance%DATA_DIRECTORY ) // &
                trim(CONTROL_instance%POTENTIALS_DATABASE)// String_getUppercase(trim(name)))

    inquire(file=trim(fileName), exist = existFile)   
    if(existFile) then
          
    !! Open File
    open(unit=30, file=trim(fileName), status="old",form="formatted")
    rewind(30)
          
    found = .false.
          
    !! Open element and Find the proper potential
    do while(found .eqv. .false.)
      read(30,*, iostat=status) token
      symbol = token(3:)

      !! Some debug information in case of error!
      if (status > 0 ) then

        call ExternalPotential_exception(ERROR, &
          "ERROR reading ExternalPotential file: "//trim(this%name)//&
          " Please check that file!","ExternalPotential module at Load function.")

      end if

      if (status == -1 ) then

        call ExternalPotential_exception(ERROR, &
          "The ExternalPotential: "//trim(this%name)//&
          " for: "//trim(species)//&
          " was not found!","ExternalPotential module at Load function.")

      end if

      if(trim(token(1:2)) == "O-") then
        if(trim(symbol) == trim(species)) then
          found = .true.

        end if

      end if

    end do

    !! Neglect any comment
    token = "#"
    do while(trim(token(1:1)) == "#")

      read(30,*) token

    end do

    !! Start reading Potential
    backspace(30)

    read(30,*, iostat=status) this%numOfComponents

    !! Some debug information in case of error!
    if (status > 0 ) then

      call ExternalPotential_exception(ERROR, &
        "ERROR reading ExternalPotential file: "//trim(this%name)//&
        " Please check that file!","ExternalPotential module at Load function.")

    end if

    allocate(this%gaussianComponents(this%numOfComponents))

    do i = 1, this%numOfComponents

      read(30,*,iostat=status) this%gaussianComponents(i)%id, &
      this%gaussianComponents(i)%angularMoment
      this%gaussianComponents(i)%length = 1
      
      !! Some debug information in case of error!
      if (status > 0 ) then

        call ExternalPotential_exception(ERROR, &
          "ERROR reading ExternalPotential file: "//trim(this%name)//&
          " Please check that file!","ExternalPotential module at Load function.")

      end if

      allocate(this%gaussianComponents(i)%orbitalExponents(this%gaussianComponents(i)%length))
      allocate(this%gaussianComponents(i)%contractionCoefficients(this%gaussianComponents(i)%length))

      do j = 1, this%gaussianComponents(i)%length

        read(30,*,iostat=status) this%gaussianComponents(i)%orbitalExponents(j), &
        this%gaussianComponents(i)%contractionCoefficients(j)
        read(30,*,iostat=status) this%gaussianComponents(i)%origin

        !! Some debug information in case of error!
        if (status > 0 ) then

          call ExternalPotential_exception(ERROR, &
            "ERROR reading ExternalPotential file: "//trim(this%name)//&
            " Please check that file!","ExternalPotential module at Load function.")

        end if

      end do

 
      !! Calculates the number of Cartesian orbitals, by dimensionality
      select case(CONTROL_instance%DIMENSIONALITY)
        case(3)
          this%gaussianComponents(i)%numCartesianOrbital = ( ( this%gaussianComponents(i)%angularMoment + 1_8 )*( this%gaussianComponents(i)%angularMoment + 2_8 ) ) / 2_8
        case(2)
          this%gaussianComponents(i)%numCartesianOrbital = ( ( this%gaussianComponents(i)%angularMoment + 1_8 ) )
        case(1)
          this%gaussianComponents(i)%numCartesianOrbital = 1 
        case default
          call ExternalPotential_exception( ERROR, &
            "Class object ExternalPotential in load function",&
            "This Dimensionality is not available") 
      end select

      !! Normalize
      allocate(this%gaussianComponents(i)%contNormalization(this%gaussianComponents(i)%numCartesianOrbital))
      allocate(this%gaussianComponents(i)%primNormalization(this%gaussianComponents(i)%length, &
        this%gaussianComponents(i)%length*this%gaussianComponents(i)%numCartesianOrbital))

      this%gaussianComponents(i)%contNormalization = 1.0_8
      this%gaussianComponents(i)%primNormalization = 1.0_8

      call ContractedGaussian_normalizePrimitive(this%gaussianComponents(i))
      call ContractedGaussian_normalizeContraction(this%gaussianComponents(i))

      !! DEBUG
      ! call ContractedGaussian_showInCompactForm(ExternalPotential_instance%potentials(potId)%gaussianComponents(i))

    end do

    close(30)

    !!DONE

    else

      call ExternalPotential_exception(ERROR, &
        "The ExternalPotential file: "//trim(name)//&
        " was not found!","ExternalPotential module at Load function.")

    end if
      
  end subroutine ExternalPotential_load

!     !>
!     !! @brief 
!     !! @param this 
!     function ExternalPotential_getPotential( this, coords )  result(output)
!         implicit none
!         type(ExternalPotential) :: this
!         real(8) :: coords(3)
!         real(8) :: output

!     !   integer :: i

!     !   output=0.0

!     !   do i=1, this%gaussianComponents%length
!     !       output = output+( this%gaussianComponents%contractionCoefficients(i)* &
!     !       exp(-this%gaussianComponents%primitives(i)%orbitalExponent*( dot_product(coords,coords) ) ) )
!     !   end do

!     end function ExternalPotential_getPotential


! !   function ExternalPotential_getInteractionMtx(this, contractions) result(output)
!     subroutine ExternalPotential_getInteractionMtx( this, contractions )
!         implicit none
!         type(ExternalPotential) :: this
!         type(ContractedGaussian) :: contractions(:)
!         type(Matrix) :: output

! !       integer :: i, j, k, l, m, a, b
! !       integer :: numContractions
! !       real(8), allocatable :: auxVal(:)
! !       type(ContractedGaussian) :: auxContract

! !       do i = 1, size(contractions)
! !           numContractions = numContractions + contractions(i)%numCartesianOrbital
! !       end do

! !       call Matrix_constructor(output,int(numContractions,8),int(numContractions,8))

! !       a = 1
! !       b = 1

! !       do i=1, size(contractions)

! !           call ContractedGaussian_product(contractions(i), &
! !                   this%gaussianComponents, auxContract)

! !           do j=1, size(contractions)

! !               call ContractedGaussian_overlapIntegral( auxContract, contractions(j), auxVal)

! !               m = 0
! !               do k = a, auxContract%numCartesianOrbital - 1
! !                   do l = b, contractions(j)%numCartesianOrbital - 1
! !                       m = m + 1
! !                       output%values(k,l)= auxVal(m)
! !                   end do
! !               end do
! !               b = b + contractions(j)%numCartesianOrbital
! !           end do
! !           a = a + auxContract%numCartesianOrbital
! !           call ContractedGaussian_destructor(auxContract)
! !       end do

! !         call Matrix_show(output)

! !         stop "ExternalPotential_getInteractionMtx"

! ! ! end function ExternalPotential_getInteractionMtx
!     end subroutine ExternalPotential_getInteractionMtx


    !>
    !! @brief  Maneja excepciones de la clase
    subroutine ExternalPotential_exception( typeMessage, description, debugDescription)
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

    end subroutine ExternalPotential_exception

end module ExternalPotential_
