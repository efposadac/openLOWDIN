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

!> @brief This module contains all the routines to handle inter particle potentials
!! @author E. F. Posada (efposadac@unal.edu.co)
!! @version 2.0

module InterPotential_
  use ContractedGaussian_
  use String_
  use Exception_
  implicit none

  type, public :: InterPot
      character(20) :: name
      character(50) :: specie
      character(50) :: otherSpecie
      character(50) :: ttype
      character(50) :: units
      integer :: numOfComponents
      integer :: iter
      type(ContractedGaussian), allocatable :: gaussianComponents(:)
  end type

  type, public :: InterPotential
    integer :: ssize
    type(InterPot), allocatable :: Potentials(:)
    logical :: isInstanced
  end type

  type(InterPotential), public, target :: InterPotential_instance

contains
  
  !>
  !! @brief Initializes the class
  !! @param this
  !! @author E. F. Posada, 2013
  subroutine InterPotential_constructor(numberOfPotentials)
      implicit none
      integer :: numberOfPotentials

      InterPotential_instance%ssize = numberOfPotentials
      allocate(InterPotential_instance%potentials(numberOfPotentials))
      InterPotential_instance%isInstanced = .true.
      
  end subroutine InterPotential_constructor

  !>
  !! @brief destroy the class
  !! @param this
  !! @author E. F. Posada, 2013
  subroutine InterPotential_destructor()
      implicit none

      integer :: i

      do i = 1, InterPotential_instance%ssize
        if (allocated(InterPotential_instance%potentials(i)%gaussianComponents) ) deallocate(InterPotential_instance%potentials(i)%gaussianComponents)
      end do

      if (allocated(InterPotential_instance%potentials) )deallocate(InterPotential_instance%potentials)

  end subroutine InterPotential_destructor

  !>
  !! @brief Shows information of the object
  !! @param this 
  subroutine InterPotential_show()
    implicit none
    integer :: i, j
    type(InterPot), pointer :: this

    do i=1,InterPotential_instance%ssize
      this => InterPotential_instance%potentials(i)
      print *,""
      print *,"======="
      print *, "InterParticle Potential for ", trim(this%specie) ," and ",  trim(this%otherSpecie), " : ", trim(this%name)
      print *, "Type : ", trim(this%ttype)
      write(6,"(T10,A20,A10,A10,A10,A20)") "Exponent", "l", "Factor"
      do j=1,this%numOfComponents
        write(6,"(T10,F20.15,I10,F20.15)") this%gaussianComponents(j)%orbitalExponents, &
        this%gaussianComponents(j)%angularMoment,  this%gaussianComponents(j)%contractionCoefficients(1)
      end do
    end do

  end subroutine InterPotential_show


  !>
  !! @brief loads information from the input file
  !! @param this
  !! @author E. F. Posada, 2015
  subroutine InterPotential_load(potId, name, species, otherSpecies)
    implicit none
    integer :: potId
    character(*) :: name
    character(*) :: species
    character(*) :: otherSpecies

    type(InterPot), pointer :: this
    integer :: status, i, j
    character(150) :: fileName
    character(20) :: token
    character(20) :: symbol
    logical :: existFile, found
    
    this => InterPotential_instance%potentials(potId)

    this%name= trim(name)
    this%specie= trim(species)
    this%otherSpecie= trim(otherSpecies)
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

        call InterPotential_exception(ERROR, &
          "ERROR reading InterPotential file: "//trim(this%name)//&
          " Please check that file!","InternalPotential module at Load function.")

      end if

      if (status == -1 ) then

        call InterPotential_exception(ERROR, &
          "The InterPotential: "//trim(this%name)//&
          " for: "//trim(species)//trim(otherSpecies)//&
          " was not found!","InternalPotential module at Load function.")

      end if

      if(trim(token(1:2)) == "O-") then
        if(trim(symbol) == trim(species)//trim(otherSpecies)) then
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

      call InterPotential_exception(ERROR, &
        "ERROR reading InternalPotential file: "//trim(this%name)//&
        " Please check that file!","InternalPotential module at Load function.")

    end if

    allocate(this%gaussianComponents(this%numOfComponents))

    do i = 1, this%numOfComponents

      read(30,*,iostat=status) this%gaussianComponents(i)%id, &
      this%gaussianComponents(i)%angularMoment
      this%gaussianComponents(i)%length = 1
      
      !! Some debug information in case of error!
      if (status > 0 ) then

        call InterPotential_exception(ERROR, &
          "ERROR reading InternalPotential file: "//trim(this%name)//&
          " Please check that file!","InternalPotential module at Load function.")

      end if

      allocate(this%gaussianComponents(i)%orbitalExponents(this%gaussianComponents(i)%length))
      allocate(this%gaussianComponents(i)%contractionCoefficients(this%gaussianComponents(i)%length))

      do j = 1, this%gaussianComponents(i)%length

        read(30,*,iostat=status) this%gaussianComponents(i)%orbitalExponents(j), &
        this%gaussianComponents(i)%contractionCoefficients(j)
        read(30,*,iostat=status) this%gaussianComponents(i)%origin

        !! Some debug information in case of error!
        if (status > 0 ) then

          call InterPotential_exception(ERROR, &
            "ERROR reading InternalPotential file: "//trim(this%name)//&
            " Please check that file!","InternalPotential module at Load function.")

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
          call InterPotential_exception( ERROR, &
            "Class object InternalPotential in load function",&
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
      ! call ContractedGaussian_showInCompactForm(InterPotential_instance%potentials(potId)%gaussianComponents(i))

    end do

    close(30)

    !!DONE

    else

      call InterPotential_exception(ERROR, &
        "The InternalPotential file: "//trim(name)//&
        " was not found!","InternalPotential module at Load function.")

    end if

  end subroutine InterPotential_load


  !>
  !! @brief  Handles class exceptions
  !<
  subroutine InterPotential_exception( typeMessage, description, debugDescription)
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

  end subroutine InterPotential_exception
end module InterPotential_
