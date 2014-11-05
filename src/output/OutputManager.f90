!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!		R. Flores (roberto.floresmoreno.qt@gmail.com)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module OutputManager_
  use Exception_
  use OutputBuilder_
  use Vector_
	implicit none

	!>
	!! @brief Description
	!!
	!! @author felix
	!!
	!! <b> Creation data : </b> 08-04-11
	!!
	!! <b> History change: </b>
	!!
	!!   - <tt> 08-04-11 </tt>:  felix ( email@server )
	!!        -# description.
        !!   - <tt> 10-31-2014 </tt>:  Jorge Charry ( jacharry@unal.edu.co )
	!!        -# Adapts this module to Lowdin2
	!!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
	!!        -# description
	!!
	!<
	type, public :: OutputManager
		character(20) :: name
		logical :: isInstanced
		integer :: numberOfOutputs
		type(OutputBuilder), allocatable :: outputs(:)
	end type

	type(OutputManager), public :: OutputManager_instance


	public :: &
		OutputManager_constructor, &
		OutputManager_destructor, &
		OutputManager_buildOutputs, &
		OutputManager_show
		
private		
contains


	!>
	!! @brief Constructor por omision
	!!
	!! @param this
	!<
	subroutine OutputManager_constructor( this, &
                                               outputType, specie, orbital, dimensions, cubeSize, point1, point2, point3)
                implicit none
		type(OutputManager) :: this
		character(*) :: outputType(:)
		character(*) :: specie(:)
		integer :: orbital(:)
		integer :: dimensions(:)
                real(8) :: cubeSize(:)
                type(Vector) :: point1(:)
                type(Vector) :: point2(:)
                type(Vector) :: point3(:)

                integer :: i

		this%numberOfOutputs= size(outputType)
		this%isInstanced = .true.

		allocate(this%outputs(this%numberOfOutputs))

		do i=1,this%numberOfOutputs

                    call OutputBuilder_constructor( this%outputs(i), i, trim(outputType(i)), trim(specie(i)), orbital(i), dimensions(i), cubeSize(i), point1(i), point2(i), point3(i)  )
		end do

	end subroutine OutputManager_constructor


	!>
	!! @brief Destructor por omision
	!!
	!! @param this
	!<
	subroutine OutputManager_destructor(this)
		implicit none
		type(OutputManager) :: this

		this%isInstanced = .false.

	end subroutine OutputManager_destructor

	!>
	!! @brief Muestra informacion del objeto
	!!
	!! @param this 
	!<
	subroutine OutputManager_show(this)
		implicit none 
		type(OutputManager) :: this
                integer :: i
                   print*, "OUTPUT FILES INFORMATION"
                   print*, "========================"
                   print*, ""

                do i=1,this%numberOfOutputs
                    call OutputBuilder_show( this%outputs(i))
		end do


	end subroutine OutputManager_show

	!!>
	!! @brief Build Different Types of Output Files
	!!
	!<
	subroutine OutputManager_buildOutputs( this )
                implicit none
		type(OutputManager) :: this
                integer :: i

		do i=1,this%numberOfOutputs
                    call OutputBuilder_buildOutput( this%outputs(i) )
		end do

	end subroutine OutputManager_buildOutputs


	!!>
	!! @brief Indica si el objeto ha sido instanciado o no
	!!
	!<
	function OutputManager_isInstanced( this ) result( output )
		implicit  none
		type(OutputManager), intent(in) :: this
		logical :: output
		
		output = this%isInstanced
	
	end function OutputManager_isInstanced

	!>
	!! @brief  Maneja excepciones de la clase
	!<
	subroutine OutputManager_exception( typeMessage, description, debugDescription)
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
	
	end subroutine OutputManager_exception

end module OutputManager_
