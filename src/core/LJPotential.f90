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

!> @brief This module contains all the routines to handle LJ potentials
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
module LJPotential_
  use String_
  use Matrix_
  use Units_
  use Exception_
  implicit none
 
  type, public :: LJPot
    character(20) :: name
    character(50) :: specie
    character(50) :: ttype
    character(50) :: units
    character(50) :: id
    integer :: numOfatoms
    integer :: iter
    real(16), allocatable :: ljparameters(:,:)
    real(16), allocatable :: atomsCenter(:,:)
  end type

  type, public :: LJPotential
    integer :: ssize
    type(LJPot), allocatable :: potentials(:)
    logical :: isInstanced
  end type

  type(LJPotential), public, target :: LJPotential_instance

contains

  !>
  !! @brief Constructor by default
  !! @param this
  subroutine LJPotential_constructor(numberOfPotentials)
    implicit none

    integer :: numberOfPotentials

    LJPotential_instance%ssize = numberOfPotentials
    allocate(LJPotential_instance%potentials(numberOfPotentials))
    LJPotential_instance%isInstanced = .true.

  end subroutine LJPotential_constructor

  !>
  !! @brief Destroys the object
  !! @param this
  subroutine LJPotential_destructor()
      implicit none
      
      integer :: i

      do i = 1, LJPotential_instance%ssize
        if (allocated(LJPotential_instance%potentials(i)%atomsCenter)) deallocate(LJPotential_instance%potentials(i)%atomsCenter)
        if (allocated(LJPotential_instance%potentials(i)%ljparameters)) deallocate(LJPotential_instance%potentials(i)%ljparameters)
      end do
     
        if (allocated(LJPotential_instance%potentials) ) deallocate(LJPotential_instance%potentials)
      LJPotential_instance%isInstanced=.false.

  end subroutine LJPotential_destructor

  !>
  !! @brief Shows information of the object
  !! @param this
  subroutine LJPotential_show()
      implicit none
      type(LJPot), pointer :: this
      integer ::  potId, i

      do potId = 1, LJPotential_instance%ssize
        this => LJPotential_instance%potentials(potId)
 
        print *,""
        print *,"======="
        print *, "LJ Potential for ", trim(this%specie), " : ", trim(this%name)
        print *, "Type : ", trim(this%ttype)
        write(6,"(T15,A12,1X,A12,A10,A15,A15)") "A_parameter", "B_parameter", "R_x", "R_y", "R_z"

        do i=1,this%numOfatoms
          write(6,"(T10,F15.5,3X,F10.5,F15.5,F15.5,F15.5)") &
          this%ljparameters(i,:), &
          this%atomsCenter(i,:)
        end do
      end do

  end subroutine LJPotential_show

  !>
  !! @brief loads information from the input file
  !! @param this
  !! @author E. F. Posada, 2013
  subroutine LJPotential_load(potId, name, species)
    implicit none
    integer :: potId
    character(*) :: name
    character(*) :: species
    
    type(LJPot), pointer :: this
    integer :: status, i, j
    character(150) :: fileName
    character(20) :: token
    character(10) :: symbol
    logical :: existFile, found
 
    this => LJPotential_instance%potentials(potId)

    this%name= trim(name)
    this%specie= trim(species)
    this%ttype=""
    this%units="bohr"
    this%numOfatoms=0
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

        call LJPotential_exception(ERROR, &
          "ERROR reading LJPotential file: "//trim(this%name)//&
          " Please check that file!","LJPotential module at Load function.")

      end if

      if (status == -1 ) then

        call LJPotential_exception(ERROR, &
          "The LJPotential: "//trim(this%name)//&
          " for: "//trim(species)//&
          " was not found!","LJPotential module at Load function.")

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

    read(30,*, iostat=status) this%numOfatoms

    !! Some debug information in case of error!
    if (status > 0 ) then

      call LJPotential_exception(ERROR, &
        "ERROR reading LJPotential file: "//trim(this%name)//&
        " Please check that file!","LJPotential module at Load function.")

    end if

    allocate(this%atomsCenter(this%numOfatoms,3))
    allocate(this%ljparameters(this%numOfatoms,2))

    
    do i = 1, this%numOfatoms

      read(30,*,iostat=status) this%id, this%ljparameters(i,:), this%atomsCenter(i,:)
      
      !! Some debug information in case of error!
      if (status > 0 ) then

        call LJPotential_exception(ERROR, &
          "ERROR reading LJPotential file: "//trim(this%name)//&
          " Please check that file!","LJPotential module at Load function.")

      end if

    end do

    close(30)

    !!DONE

    else

      call LJPotential_exception(ERROR, &
        "The LJPotential file: "//trim(name)//&
        " was not found!","LJPotential module at Load function.")

    end if
      
  end subroutine LJPotential_load

    !>
    !! @brief  Maneja excepciones de la clase
    subroutine LJPotential_exception( typeMessage, description, debugDescription)
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

    end subroutine LJPotential_exception

end module LJPotential_
