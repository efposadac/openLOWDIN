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

!> @brief This module contains all the routines to confining Potentials
!! @author I. Ortiz-Verano (ieortizv@unal.edu.co)
!! @version 2.0
!! <b> Creation data : </b> 20-02-2018
!!
!! <b> History change: </b>
!!
!!   - <tt> 20-02-1018 </tt>:  Ismael Ortiz-Verano ( ieortizv@unal.edu.co )
!!        -# Creacioon del modulo y metodos basicos
module ConfiningPotential_
  use String_
  use Exception_
  use AtomicElement_
  use MolecularSystem_

  implicit none

  type, public :: ConfPot
     real(8) :: exponent
     real(8) :: radius
     real(8) :: rZero
     real(8) :: confiningCoefficient
     character(30) :: particle
     character(5) :: elementSymbol
  end type ConfPot

contains

  !>
  !! @brief Shows information relative to confining potential
  !! @param this
  subroutine ConfiningPotential_show()
    implicit none
    type(ConfPot),allocatable :: this(:)
    integer :: quantumSpecies, i
    character(5) :: elementSymbol

    quantumSpecies = MolecularSystem_instance%numberOfQuantumSpecies
    allocate(this(quantumSpecies))

    call ConfiningPotential_loadParticles(elementSymbol)
    this(1)%elementSymbol=elementSymbol

    print *,""
    print *,"======="
    print *, "Confining Potential for atom ",this(1)%elementSymbol
    write(6,"(T10,A10,A10,A20,A20,A20,A20)") "Element", "Particle", "Covalent radius", "Conf. exponent", "Conf. coefficient","Conf. radius"


    do i = 1, quantumSpecies
       this%particle = MolecularSystem_getNameOfSpecie(i)
       call ConfiningPotential_constructPotential(this(i))
       write(6,"(T10,A10,A10,F20.10,F20.10,F20.10,F20.10)") &
            this(i)%elementSymbol,this(i)%particle,this(i)%radius,this(i)%exponent, this(i)%confiningCoefficient, this(i)%rZero
    end do

  end subroutine ConfiningPotential_show


  !>
  !! @brief construct the confining potential
  !! @param this
  subroutine ConfiningPotential_constructPotential(this)
    implicit none
    type(ConfPot):: this
    type(AtomicElement) :: element
    real(8) :: confiningCoefficient
    character(5) :: elementSymbol
    integer :: aux,mass

    call ConfiningPotential_loadParticles(elementSymbol)
    this%elementSymbol=elementSymbol

    aux = index(this%particle,"_")
    if (aux == 0) then
       this%elementSymbol=elementSymbol
       confiningCoefficient = 1.85_8
       this%confiningCoefficient = confiningCoefficient
       this%exponent = 2.0_8
       call AtomicElement_load(element, trim(elementSymbol), mass)
       this%radius = element%covalentRadius
       this%rZero = this%radius * confiningCoefficient
    else
       this%elementSymbol=elementSymbol
       confiningCoefficient = 5.85_8
       this%confiningCoefficient = confiningCoefficient
       this%exponent = 1.5_8
       call AtomicElement_load(element, trim(elementSymbol), mass)
       this%radius = 0.1_8*element%covalentRadius
       this%rZero = this%radius * confiningCoefficient
    end if
  end subroutine ConfiningPotential_constructPotential

  !>
  !! @brief loads some necesary information from the input file
  !! @param this
  !! @author I. Ortiz-Verano, 2018
  subroutine ConfiningPotential_loadParticles(elementSymbol)
    implicit none
    character(100) :: auxFile
    logical :: existFile
    logical :: foundElement
    character(20) :: token, ttoken,tttoken
    integer :: status,aux
    character(5), intent(out) :: elementSymbol


    auxFile=trim(CONTROL_instance%INPUT_FILE)//"aux"

    inquire(FILE = trim(auxFile), EXIST = existFile )


    if ( existFile ) then

       open(unit=99, file=trim(auxFile), status="old",form="formatted")
       rewind(99)

       foundElement = .false.

       do while(foundElement .eqv. .false.)

          read(99,*, iostat=status) token

          if(trim(token) == "InputParticle_name") then

             foundElement = .true.

             backspace(99)
             read(99,*, iostat=status) token, ttoken, tttoken
             !! Debug
             !! print*,"cuántos imprime? ",tttoken
             elementSymbol = tttoken (scan( tttoken, "[" ) : scan( tttoken, "]" ))

             aux = index(elementSymbol,"]")
             elementSymbol = trim(elementSymbol(2:aux-1))
             elementSymbol = trim(String_getUppercase(trim(elementSymbol)))

             !! Debug
             !!print*,"Símbolo químico!!!!!!!!!!!!!!!!!!!!!!!!!",elementSymbol


          end if
       end do

    else
       print*,"*_*_*_*_*_*_*_*_*_*_   ERROR   _*_*_*_*_*_*_*_*_*_*"
       call ConfiningPotential_exception(ERROR, "ERROR reading chemical symbol (Confining Potential) from file: "//trim(auxFile)//" This file doesn't exists!","ConfiningPotential module at loadParticles function.")
    end if

  end subroutine ConfiningPotential_loadParticles


  subroutine ConfiningPotential_exception( typeMessage, description, debugDescription)
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

  end subroutine ConfiningPotential_exception

end module ConfiningPotential_
