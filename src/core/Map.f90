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
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module Map_
  use CONTROL_
  use Exception_
  implicit none

	!>
	!! @brief Defincion de clase para mapas
	!!
	!! Define una seudoclase que alamcena elementos promados por una combinacion clave-valor
	!!
	!! @author Sergio A. Gonzalez Monico
	!!
	!! <b> Fecha de creacion : </b> 2008-09-03
	!!
	!! <b> Historial de modificaciones: </b>
	!!
	!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
	!!        -# Creacion de modulo y metodos
	!!   - <tt> 2014-05-14 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin2
	!!
	!! @todo Implementar mapa mediante una arbol binario, emplenado punteros
	!<

	type, public :: Map

		character(30) :: name
		character(30), allocatable :: key(:)
		real(8) , allocatable :: value(:)
		integer :: first
		integer :: last

	end type Map

	type(Map), target, public :: currentMap

	public :: &
		Map_constructor , &
 		Map_destructor, &
		! Map_show, &
		Map_insert, &
		! Map_erase, &
		Map_find, &
		! Map_swap, &
		Map_getValue
		! Map_getKey, &
		! Map_size, &
		! Map_begin, &
		! Map_end, &
		! Map_clear, &
		! Map_isEmpty, &
		! Map_doeskeyExist

contains

	!<
	!! Define el constructor para la clase
	!!
	!>
	subroutine Map_constructor( this, ssize, name )
		implicit none
		type(Map), intent(inout) :: this
		integer, optional :: ssize
		character(*), optional :: name

		integer :: auxSize

		auxSize = 0
		if ( present(ssize) ) auxSize = ssize
		if ( present( name) )  this%name = trim(name)

		if ( allocated( this%key ) ) deallocate( this%key )
		allocate( this%key( auxSize ) )

		if ( allocated( this%value ) ) deallocate( this%value )
		allocate( this%value( auxSize ) )


		this%first = 1
		this%last = auxSize

	end subroutine Map_constructor


	!<
	!! Define el destructor para clase
	!!
	!>
	subroutine Map_destructor( this)
		implicit none
		type(Map), intent(inout) :: this

		if ( allocated( this%key ) ) deallocate( this%key )
		if ( allocated( this%value ) ) deallocate( this%value )
		this%name = "undefined"
		this%first = 0
		this%last = 0

	end subroutine Map_destructor

	!<
	!! Muestra el contenido del mapa especificado
	!!
	!>
	! subroutine Map_show( this)
	! 	implicit none
	! 	type(Map), intent(in) :: this

	! 	integer :: i

	! 	if ( allocated( this%key ) ) then
	! 		print *,""
	! 		print *,"===="
	! 		print *,"Map: "
	! 		print *,"===="
	! 		print *,""
	! 		do i=1, this%last
	! 			write (6,"(T10,A15,F10.2)") this%key(i), this%value(i)
	! 		end do
	! 		print *,""

	! 	end if

	! end subroutine Map_show


	! !<
	! !! Adiciona un elemento al mapa especificado
	! !!
	! !>
	subroutine Map_insert( this, key, value )
		implicit none
		type(Map), intent(inout), target :: this
		character(*), intent(in) :: key
		real(8), intent(in) :: value

		type(Map) :: auxMap
		integer :: i
		logical :: auxLogical
		type(Exception) :: ex

		if ( Map_find( this, trim(key) ) == 0 ) then

			call Map_constructor( auxMap,  this%last )

			auxMap%key = this%key
			auxMap%value = this%value

			call Map_constructor( this,  this%last + 1 )

			this%key(1 : auxMap%last ) = auxMap%key
			this%value(1 : auxMap%last) = auxMap%value

			this%key( this%last ) = trim( key )
			this%value( this%last ) = value

			call Map_destructor( auxMap )

		else

			call Map_exception( ERROR, "Class object Map in the insert() function", &
								"The key: "//trim(key)//" already was defined " )

		end if

	end subroutine Map_insert

	! !<
	! !! Borra un elemento del mapa
	! !!
	! !>
	! subroutine Map_erase( this, key, iterator )
	! 	implicit none
	! 	type(Map), intent(inout), target :: this
	! 	character(*), intent(in), optional :: key
	! 	integer, intent(in), optional :: iterator

	! 	type(Map) :: auxMap
	! 	integer :: i,j



	! 	if ( Map_find( this, trim(key) ) /= 0 ) then

	! 		call Map_constructor( auxMap, this%last - 1)

	! 		j=0
	! 		do  i = 1, this%last
	! 			if ( trim( this%key(i) ) /= trim( key ) ) then
	! 				j=j+1
	! 				auxMap%key(j) = this%key(i)
	! 				auxMap%value(j) = this%value(i)
	! 			end if
	! 		end do

	! 		call Map_constructor( this, this%last - 1)

	! 		this%key = auxMap%key
	! 		this%value = auxMap%value

	! 		call Map_destructor( auxMap )

	! 	end if

	! end subroutine Map_erase

	! !<
	! !! Invierte la psocion de dos elementos
	! !!
	! !>
	! subroutine Map_swap( this, keyA, keyB )
	! 	implicit none
	! 	type(Map), intent(inout), target :: this
	! 	character(*), intent(in) :: keyA
	! 	character(*), intent(in) :: keyB

	! 	integer :: i
	! 	integer :: j
	! 	integer :: auxValue
	! 	character(30) :: auxKey

	! 	i = Map_find( this, trim(keyA) )
	! 	J = Map_find( this, trim(keyB) )

	! 	if ( ( i /= 0 ) .and. ( j /=0 )) then

	! 		auxKey = trim( this%key(i))
	! 		this%key(i) = trim( this%key(j) )
	! 		this%key(j) = trim( auxKey )

	! 		auxValue= this%value(i)
	! 		this%value(i) = this%value(j)
	! 		this%value(j) = auxValue

	! 	end if

	! end subroutine Map_swap

	!>
	!! @brief Busca un elemento del mapa y retorna su iterador
	!!
	!<
	function Map_find( this, key ) result( output )
		implicit none
		type(Map), intent(in), target :: this
		character(*), intent(in), optional :: key
		integer :: output

		integer :: i

		do  i = 1, this%last
			if ( trim( this%key(i) ) == trim( key ) ) exit
		end do


		if ( i > this%last ) then
			output = 0
		else
			output = i
		end if



	end function Map_find

	! function Map_doeskeyExist( this, key ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	character(*), intent(in), optional :: key
	! 	integer :: output

	! 	integer :: i
	! 	output = 0
	! 	do  i = 1, this%last
	! 		if ( trim( this%key(i) ) == trim( key ) )then
	! 			output = 1
	! 			exit
	! 		end if

	! 	end do
	! end function Map_doeskeyExist

	! !<
	! !! Busca un elemento del mapa y retorna su valor
	! !!
	! !>
	function Map_getValue( this, key, iterator ) result( output )
		implicit none
		type(Map), intent(in), target :: this
		character(*), intent(in), optional :: key
		integer, intent(in), optional :: iterator
		real(8) :: output

		integer :: i

		if ( present(key) ) then

			do  i = 1, this%last
				if ( trim( this%key(i) ) == trim( key ) ) exit
			end do

		else if ( present(iterator) ) then

			i= iterator

		end if

		if ( i >  this%last ) then
			output = this%value( this%last )
		else
			output = this%value( i )
		end if


	end function Map_getValue

	! !<
	! !! @brief Retorna la clave aociada a un valor dado
	! !!
	! !>
	! function Map_getKey( this, value, iterator ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	real(8), intent(in), optional :: value
	! 	integer, intent(in), optional :: iterator
	! 	character(30) :: output

	! 	integer :: i

	! 	if ( present(value) ) then

	! 		do  i = 1, this%last
	! 			if (  abs(this%value(i) - value)  < Parameters%DOUBLE_ZERO_THRESHOLD ) exit
	! 		end do

	! 	else if ( present(iterator) ) then

	! 		i= iterator

	! 	end if

	! 	if ( i >  this%last ) then
	! 		output =  this%key( this%last)
	! 	else
	! 		output = this%key( i )
	! 	end if


	! end function Map_getKey

	! !<
	! !! Devuelve el tamano del mapa
	! !!
	! !>
	! function Map_size( this ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	integer :: output

	! 	output = this%last

	! end function Map_size

	! !<
	! !! Devuelve un iterador al primer elemento del mapa
	! !!
	! !>
	! function Map_begin( this ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	integer :: output

	! 	output = this%first

	! end function Map_begin

	! !<
	! !! Devuelve un iterador al �ltimo elemento del mapa
	! !!
	! !>
	! function Map_end( this ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	integer :: output

	! 	output = this%last

	! end function Map_end

	! !<
	! !! Limpiar el mapa
	! !!
	! !>
	! subroutine Map_clear( this )
	! 	implicit none
	! 	type(Map), intent(inout), target :: this
	! 	integer :: output

	! 	call Map_destructor( this )

	! end subroutine Map_clear

	! !<
	! !! Indica si el arreglo esta vacio o no
	! !!
	! !>
	! function Map_isEmpty( this ) result( output )
	! 	implicit none
	! 	type(Map), intent(in), target :: this
	! 	logical :: output
		
	! 	output = .false.
	! 	if ( ( this%last - this%first ) <= 0 )	output = .true.

	! end function Map_isEmpty

	!>
	!! @brief  Maneja excepciones de la clase
	!<
	subroutine Map_exception( typeMessage, description, debugDescription)
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
	
	end subroutine Map_exception

end module Map_
