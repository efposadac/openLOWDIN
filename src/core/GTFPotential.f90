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

!> @brief This module contains all the routines to handle external and interal GTF potentials
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
!!   - <tt> 2024-11-26 </tt>: Felix
!!        -# Merges ExternalPotential and InternalPotentials modules into a single file (GTFPotential)
module GTFPotential_
  use ContractedGaussian_
  use String_
  use Matrix_
  use Units_
  use Exception_
  implicit none

  type, public :: GaussPot
     character(50) :: name
     character(50) :: species
     character(50) :: otherSpecies
     character(50) :: ttype
     character(50) :: units
     integer :: numOfComponents
     integer :: iter
     type(ContractedGaussian), allocatable :: gaussianComponents(:)
  end type GaussPot

  type, public :: GTFPotential
     integer :: ssize
     type(GaussPot), allocatable :: potentials(:)
     character(50) :: type
     logical :: isInstanced
  end type GTFPotential

  type(GTFPotential), public, target :: ExternalPotential_instance, InterPotential_instance

contains

  !>
  !! @brief Initializes the class
  !! @param this, n
  !! @author E. F. Posada, 2013
  subroutine GTFPotential_constructor(this,numberOfPotentials,type)
    implicit none
    type(GTFPotential) :: this
    integer :: numberOfPotentials
    character(*) :: type


    this%ssize = numberOfPotentials
    allocate(this%potentials(numberOfPotentials))
    this%isInstanced = .true.
    this%type = type

  end subroutine GTFPotential_constructor

  !>
  !! @brief Destroys the object
  !! @param this
  subroutine GTFPotential_destructor(this)
    implicit none
    type(GTFPotential) :: this

    integer :: i

    do i = 1, this%ssize
       if (allocated(this%potentials(i)%gaussianComponents)) deallocate(this%potentials(i)%gaussianComponents)
    end do

    if (allocated(this%potentials) ) deallocate(this%potentials)
    this%isInstanced=.false.

  end subroutine GTFPotential_destructor

  !>
  !! @brief loads information from the input file
  !! @param this
  !! @author E. F. Posada, 2013
  subroutine GTFPotential_load(this, type, potId, name, species, otherSpecies)
    implicit none
    type(GTFPotential) :: this
    character(*) :: type
    integer :: potId
    character(*) :: name
    character(*) :: species
    character(*) :: otherSpecies

    if(trim(type) .eq. "INTERNAL") then
       call GaussPot_load(this%potentials(potId), name, trim(type), trim(species), trim(otherSpecies))
    else if(trim(type) .eq. "EXTERNAL") then
       call GaussPot_load(this%potentials(potId), name, trim(type), trim(species), "")
    else
       call Exception_stopError("Unknown GTFPotential potential type requested: "//trim(type),"GTFPotential module at Load function.")
    end if
  end subroutine GTFPotential_load

  !>
  !! @brief Shows information of the object
  !! @param this 
  subroutine GTFPotential_show(this)
    implicit none
    type(GTFPotential) :: this
    integer :: i, j

    do i=1,this%ssize       
       if( this%potentials(i)%ttype .eq. "INTERNAL") then
          write(*,*) ""
          write(*,*)"======="
          write(*,"(A30,A)") "GTF Interparticle potential: ", trim(this%potentials(i)%name)
          write(*,"(A4,A10,A5,A10)") "for ", trim(this%potentials(i)%species) ," and ",  trim(this%potentials(i)%otherSpecies)
          write(*,"(T10,A10,A10)") "Units:", trim(this%potentials(i)%units)
          write(*,"(T10,A16,A16)") "Exponent", "Factor"
          do j=1,this%potentials(i)%numOfComponents
             write(*,"(T10,E16.8,E16.8)") this%potentials(i)%gaussianComponents(j)%orbitalExponents, &
                  this%potentials(i)%gaussianComponents(j)%contractionCoefficients(1)
          end do
       else if( this%potentials(i)%ttype .eq. "EXTERNAL") then
          write(*,*) ""
          write(*,*) "======="
          write(*,"(A25,A20,A5,A10)") "GTF External potential: ", trim(this%potentials(i)%name), " for ", trim(this%potentials(i)%species)
          write(*,"(T10,A10,A10)") "Units:", trim(this%potentials(i)%units)
          write(*,"(T10,A16,A10,A10,A10,A16)") "Exponent", "R_x", "R_y", "R_z", "Factor"

          do j=1,this%potentials(i)%numOfComponents
             write(*,"(T10,E16.8,F10.5,F10.5,F10.5,E16.8)") &
                  this%potentials(i)%gaussianComponents(j)%orbitalExponents(1), &
                  this%potentials(i)%gaussianComponents(j)%origin(:), &
                  this%potentials(i)%gaussianComponents(j)%contractionCoefficients(1)
          end do
       end if
    end do

  end subroutine GTFPotential_show

  !>
  !! @brief loads information from the input file
  !! @param this
  !! @author Felix, 2024
  subroutine GaussPot_load(this, name, type, species, otherSpecies)
    type(GaussPot) :: this
    character(*) :: name
    character(*) :: type
    character(*) :: species
    character(*) :: otherSpecies

    integer :: status, i, j
    character(150) :: fileName
    character(20) :: token
    character(50) :: symbol
    logical :: existFile, found

    this%name= trim(name)
    this%species= trim(species)
    this%ttype=type
    this%units="BOHR"
    this%numOfComponents=0
    this%iter=1

    if(trim(this%ttype) .eq. "EXTERNAL" ) then
       this%otherSpecies=""
    else
       this%otherSpecies=trim(otherSpecies)
    end if

    fileName = trim( trim( CONTROL_instance%DATA_DIRECTORY ) // &
         trim(CONTROL_instance%POTENTIALS_DATABASE)// String_getUppercase(trim(name)))

    !! Open Potential file from library
    inquire(file=trim(fileName), exist = existFile)
    if(existFile) then
       open(unit=30, file=trim(fileName), status="old",form="formatted")
    else
       !! Open Potential file from directory
       inquire(file=trim(this%name), exist = existFile)
       if(existFile) then
          open(unit=30, file=trim(this%name), status="old",form="formatted")
       else
          !! File not found
          call Exception_stopError("The GTFPotential file: "//trim(this%name)//&
               " was not found!","GTFPotential module at Load function.")
       end if
    end if

    !! Open File
    rewind(30)

    found = .false.

    !! Open element and Find the proper potential
    do while(found .eqv. .false.)
       read(30,*, iostat=status) token
       symbol = token(3:)

       !! Some debug information in case of error!
       if (status > 0 ) call Exception_stopError("ERROR reading InterPotential file: "//trim(this%name)//&
            " Please check that file!","GTFPotential module at Load function.")

       if (status == -1 ) &
            call Exception_stopError("The "//trim(this%ttype)//" Potential: "//trim(this%name)//&
            " for: "//trim(species)//trim(otherSpecies)//&
            " was not found!","GTFPotential module at Load function.")

       if(trim(token(1:2)) == "O-") then
          if(this%ttype=="EXTERNAL" .and. trim(symbol) == trim(species)) found = .true.             
          if(this%ttype=="INTERNAL" .and. trim(symbol) == trim(species)//trim(otherSpecies)) found = .true.             
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
    if (status > 0 ) call Exception_stopError("ERROR reading ExterPotential file: "//trim(this%name)//&
         " Please check that file!","GTFPotential module at Load function.")

    allocate(this%gaussianComponents(this%numOfComponents))

    do i = 1, this%numOfComponents
       read(30,*,iostat=status) this%gaussianComponents(i)%id, &
            this%gaussianComponents(i)%angularMoment

       if(this%gaussianComponents(i)%angularMoment .gt. 0) then
          call Exception_sendWarning("you provided a non-zero angular momentum in GTFpotential "//trim(this%name)//&
               ". This feature is not yet implemented, will be ignored and set to zero", "GTFPotential module at Load function.")
          this%gaussianComponents(i)%angularMoment=0
       end if

       this%gaussianComponents(i)%length = 1

       !! Some debug information in case of error!
       if (status > 0 ) call Exception_stopError("ERROR reading InternalPotential file: "//trim(this%name)//&
            " Please check that file!","GTFPotential module at Load function.")

       allocate(this%gaussianComponents(i)%orbitalExponents(this%gaussianComponents(i)%length))
       allocate(this%gaussianComponents(i)%contractionCoefficients(this%gaussianComponents(i)%length))

       do j = 1, this%gaussianComponents(i)%length

          read(30,*,iostat=status) this%gaussianComponents(i)%orbitalExponents(j), &
               this%gaussianComponents(i)%contractionCoefficients(j)
          read(30,*,iostat=status) this%gaussianComponents(i)%origin

          !! Some debug information in case of error!
          if (status > 0 ) call Exception_stopError("ERROR reading InternalPotential file: "//trim(this%name)//&
               " Please check that file!","GTFPotential module at Load function.")

       end do

       if(this%ttype=="INTERNAL" .and. sum(this%gaussianComponents(i)%origin(:)**2) .gt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
          call Exception_sendWarning("you provided a non-zero origin for interpotential "//trim(this%name)//&
               ". This feature is not yet implemented, will be ignored and set to zero", "GTFPotential module at Load function.")
          this%gaussianComponents(i)%origin=0.0_8
       end if

       !! Calculates the number of Cartesian orbitals, by dimensionality
       select case(CONTROL_instance%DIMENSIONALITY)
       case(3)
          this%gaussianComponents(i)%numCartesianOrbital = ( ( this%gaussianComponents(i)%angularMoment + 1_8 )*( this%gaussianComponents(i)%angularMoment + 2_8 ) ) / 2_8
       case(2)
          this%gaussianComponents(i)%numCartesianOrbital = ( ( this%gaussianComponents(i)%angularMoment + 1_8 ) )
       case(1)
          this%gaussianComponents(i)%numCartesianOrbital = 1 
       case default
          call Exception_stopError("Class object InternalPotential in load function",&
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

    end do

    close(30)

    !!DONE
  end subroutine GaussPot_load

end module GTFPotential_
