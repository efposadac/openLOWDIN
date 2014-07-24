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

module ListInteger_
  use Exception_
  implicit none

	!>
	!! @brief  Modulo para implementacion de listas de enteros
	!!		Esta es un implementacion temporal, ver List para detalles
	!!
	!! @author Sergio A. Gonzalez
	!!
	!! <b> Fecha de creacion : </b> 2008-09-19
	!!   - <tt> 2007-08-19 </tt>: Sergio A. Gonz�lez ( sagonzalezm@unal.edu.co )
	!!        -# Creacion del archivo y las funciones basicas
	!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
	!!   - <tt> 2014-05-14 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin2
	!!
	!! @todo implementar empleando punteros como una lista enlazada de estos
	!! @todo Extender el m\'etodo para aceptar lista de vetores y matrices
	!<

type, public :: ListInteger

		character(30) :: name
		integer, allocatable :: data(:)
		integer, pointer :: current
		integer :: iterator
		integer :: maxSize
		logical :: isSizeUndefined

	end type ListInteger


	public :: &
		ListInteger_constructor, &
		ListInteger_assign, &
		ListInteger_back, &
		ListInteger_begin, &
		ListInteger_clear, &
		ListInteger_empty, &
		ListInteger_end, &
		ListInteger_erase, &
		ListInteger_front, &
		ListInteger_insert, &
		ListInteger_max_size, &
		ListInteger_merge, &
		ListInteger_pop_back, &
		ListInteger_pop_front, &
		ListInteger_remove, &
		ListInteger_push_back, &
		ListInteger_push_front, &
		ListInteger_resize	, &
		ListInteger_reverse, &
		ListInteger_size, &
		ListInteger_sort, &
		ListInteger_swap , &
		ListInteger_unique, &
		ListInteger_iterate, &
		ListInteger_current

		private :: &
			ListInteger_exception

contains

	!>
	!! @brief Constructor
	!!	Crea una lista y la inicializa
	!<
	subroutine ListInteger_constructor( this, name, ssize )
		implicit none
		type(ListInteger), target, intent(inout) :: this
		character(*), optional :: name
		integer, optional, intent(in) :: ssize

		this%maxSize =10
		this%name = "undefined"
		this%isSizeUndefined = .false.

		if ( present(name) ) this%name = trim(name)

		if (present(ssize)) then

			this%maxSize =ssize
			if( ssize < 0 ) then
				this%isSizeUndefined = .true.
				this%maxSize =1
			end if

		end if

		allocate( this%data(1) )
		this%data = 0
		this%iterator = 0
		this%current => null()

	end  subroutine ListInteger_constructor

	!>
	!! @brief Destructor
	!!            Elimina una lista previamente creada y libera la memoria asociada
	!!
	!! @ Verificar que la memoria priviamente separada sea liberada
	!<
	subroutine ListInteger_destructor( this )
		implicit none
		type(ListInteger), intent(inout) :: this

		integer :: i
		integer :: currentSize

		if ( allocated(this%data) ) then

			this%current => null()
			this%maxSize = 0
			this%iterator = 0
			deallocate(this%data)

		end if

	end  subroutine ListInteger_destructor

	!>
	!! @brief Asigna un grupo elementos a la lista
	!!
	!! @todo Falta implementar
	!<
	subroutine ListInteger_assign(this, first, last, value)
		implicit none
		type(ListInteger), intent(inout) :: this
		integer :: first
		integer :: last
		integer :: value

	end subroutine ListInteger_assign

	!>
	!! @brief Retorna un puntero al ultimo elemento de la lista
	!!
	!<
	function ListInteger_back( this ) result(output)
		implicit none
		type(ListInteger), target, intent(inout) :: this
		integer, pointer :: output

		if ( allocated( this%data) ) then

			output => null()
			this%iterator = size( this%data )
			output => this%data(this%iterator)


		else

			call ListInteger_exception( ERROR, "The list "// trim(this%name) //" is empty", "Class object ListInteger in the back() function")

		end if

	end function ListInteger_back

	!>
	!! @brief Coloca el puntero del iterador al inicio de la lista
	!
	!<
	subroutine ListInteger_begin( this)
		type(ListInteger), target, intent(inout) :: this

		if ( allocated(this%data) ) then

			this%current => null()
			this%current => this%data(1)
			this%iterator = 1

		else

			call ListInteger_exception(ERROR, "The list "// trim(this%name) //" is empty", "Class object ListInteger in the begin() function" )

		end if

	end subroutine ListInteger_begin

	!>
	!! @brief Remueve todo los elementos del la lista
	!!
	!<
	subroutine ListInteger_clear( this)
		type(ListInteger), intent(inout) :: this

		if ( allocated(this%data) ) then

			this%current=> null()
			deallocate( this%data )
			allocate( this%data(1) )
			this%data = 0.0_8
			this%iterator = 0


		end if

	end subroutine ListInteger_clear

	!>
	!! @brief Indica si la lista tiene o no elementos
	!<
	function ListInteger_empty( this ) result(output)
		implicit none
		type(ListInteger), intent(in) :: this
		logical :: output

		output=.false.

		if( .not.allocated(this%data) .or. this%iterator == 0 ) output =.true.

	end function ListInteger_empty

	!>
	!! @brief Coloca el puntero del iterador al final de la lista
	!!
	!<
	subroutine ListInteger_end( this)
		type(ListInteger), target, intent(inout) :: this

		if ( allocated( this%data ) ) then

			this%current=>null()
			this%iterator = size(this%data)
			this%current => this%data( this%iterator )

		else

			call ListInteger_exception( ERROR, "The list "// trim(this%name) //" is empty", "Class object ListInteger in the end() function" )

		end if

	end subroutine ListInteger_end

	!>
	!! @brief Remueve todo los elementos del la lista
	!!
	!! @todo Falta por implementar
	!<
	subroutine ListInteger_erase( this)
		type(ListInteger), intent(inout) :: this


	end subroutine ListInteger_erase


	!>
	!! @brief Retorna una puntero al primer elemento de la lista
	!!
	!<
	function ListInteger_front( this ) result(output)
		implicit none
		type(ListInteger), target, intent(inout) :: this
		integer, pointer :: output

		if ( allocated( this%data ) ) then

			output => null()
			output => this%data(1)

		else

			call ListInteger_exception( ERROR, "The list "// trim(this%name) //" is empty", "Class object ListInteger in the back() function" )

		end if

	end function ListInteger_front

	!>
	!! @brief Inserta un elemento en la posicion especificada
	!!
	!! @todo Falta implementacion completa
	!<
	subroutine ListInteger_insert( this, position)
		type(ListInteger), intent(inout) :: this
		integer :: position

	end subroutine ListInteger_insert

	!>
	!! @brief Retorna el maximo numero de elementos que se puden almacenar
	!! 		en la lista.
	!<
	function ListInteger_max_size( this) result(output)
		implicit none
		type(ListInteger),  intent(in) :: this
		integer :: output

		output = this%maxSize

	end function ListInteger_max_size

	!>
	!! @brief Fusiona dos lista en una borrando un de ellas
	!!
	!! @todo Falta implementacion completa
	!<
	subroutine ListInteger_merge( this, otherThisToMerge )
		type(ListInteger), intent(inout) :: this
		type(ListInteger), intent(inout) :: otherThisToMerge


	end subroutine ListInteger_merge


	!>
	!! @brief Elimina el �ltimo elemento de la lista, el primer elemento no es eliminado
	!<
	subroutine ListInteger_pop_back( this )
		implicit none
		type(ListInteger), target, intent(inout) :: this

		integer, allocatable :: auxData(:)
		integer :: ssize

		if ( allocated(this%data)  ) then

			if ( size(this%data) > 1 ) then

				this%current  => null()
				ssize = size(this%data) - 1
				allocate( auxData(ssize) )
				auxData = this%data(1:ssize)
				deallocate( this%data )
				allocate( this%data(ssize) )
				this%data = auxData
				this%iterator = ssize
				this%current => this%data(ssize)
				deallocate( auxData )

			end if

		else

			call ListInteger_exception( WARNING, "The list is empty", "Class object ListInteger in the pop_back() function" )

		end if

	end subroutine ListInteger_pop_back

	!>
	!! @brief Elimina el primer elemento de la lista
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_pop_front( this )
		implicit none
		type(ListInteger), intent(inout) :: this

	end subroutine ListInteger_pop_front

	!>
	!! @brief Adiciona un elemento al final a la lista
	!!
	!! @ Verificar que la memoria de cada nodo sea liberada cuando �ste se elimine
	!<
	subroutine ListInteger_push_back( this, data )
		implicit none
		type(ListInteger), target, intent(inout) :: this
		integer, intent(in) :: data

		integer, allocatable :: auxData(:)
		integer :: ssize

		if ( allocated(this%data)  ) then
		
			ssize = size( this%data )

			if ( this%isSizeUndefined ) this%maxSize = this%maxSize + 1

			if ( this%iterator /= 0 .and. ssize < this%maxSize ) then

				this%current => null()
				allocate( auxData( ssize ) )
				auxData = this%data
				deallocate(this%data)
				allocate( this%data( ssize + 1 ) )
				this%data(1:ssize) = auxData
				this%data(ssize+1) = data
				this%iterator = ssize+1
				this%current => this%data(ssize+1)
				deallocate(auxData)

			else if ( ssize== this%maxSize) then

				this%current => null()
				allocate( auxData( ssize ) )
				auxData = this%data
				this%data(1:ssize-1) = auxData(2:ssize)
				this%data(ssize) = data
				this%current => this%data(ssize)
				deallocate(auxData)

			else

				this%data(1) = data
				this%iterator = 1
				this%current => this%data(1)

			end if

		else

			call ListInteger_exception(ERROR, "The list "// trim(this%name) //" is empty", "Class object ListInteger in the begin() function" )

		end if

	end subroutine ListInteger_push_back

	!>
	!! @brief Adiciona un elemento al inicio de la lista
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_push_front( this )
		implicit none
		type(ListInteger), intent(inout) :: this

	end subroutine ListInteger_push_front

	!>
	!! @brief Remueve el nodo especificado de la lista
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_remove( this, position )
		implicit none
		type(ListInteger), intent(inout) :: this
		integer :: position

	end subroutine ListInteger_remove

	!>
	!! @brief Redefine el tama�o maximo de la lista
	!!
	!<
	subroutine ListInteger_resize( this, resize )
		implicit none
		type(ListInteger), intent(inout) :: this
		integer :: resize

		this%maxSize = resize

	end subroutine ListInteger_resize

	!>
	!! @brief Invierte la lista
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_reverse( this)
		implicit none
		type(ListInteger), intent(inout) :: this

	end subroutine ListInteger_reverse

	!>
	!! @brief Retorna el tama�o actual o ocupacion de la pila
	!<
	function ListInteger_size( this) result(output)
		implicit none
		type(ListInteger),  intent(in) :: this
		integer :: output

		output=0

		if( allocated(this%data) ) then

			output =size( this%data )

		end if

	end function ListInteger_size

	!>
	!! @brief Ordena los elementos de la lista en orden ascendente
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_sort( this)
		implicit none
		type(ListInteger), intent(inout) :: this


	end subroutine ListInteger_sort

	!>
	!! @brief Intercambia los elementos entre dos listas
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_swap( this, otherThis)
		implicit none
		type(ListInteger), intent(inout) :: this
		type(ListInteger), intent(inout) :: otherThis


	end subroutine ListInteger_swap

	!>
	!! @brief remueve los elemetos repetidos de la lista siempre que
	!!		sean consecutivos
	!!
	!! @todo Falta toda la implemtacion
	!<
	subroutine ListInteger_unique( this)
		implicit none
		type(ListInteger), intent(inout) :: this


	end subroutine ListInteger_unique

	!>
	!! @brief Mueve el iterador de la lista el numero de nodo que se especifique
	!!
	!<
	subroutine ListInteger_iterate( this, iterator)
		implicit none
		type(ListInteger), target, intent(inout) :: this
		integer, optional :: iterator

		integer :: auxIterator
		integer :: i

		auxIterator=1
		if ( present(iterator) ) auxIterator = iterator

		if ( allocated(this%data) ) then

			if( ( ( this%iterator + auxIterator ) >= 1 ) .and.  ( ( this%iterator + auxIterator) <= this%maxSize  )  ) then

				do i=1, abs(auxIterator)

					if ( auxIterator > 0) then

						this%iterator = this%iterator + 1
						this%current => this%data(this%iterator)

					else

						this%iterator = this%iterator - 1
						this%current => this%data(this%iterator)

					end if

				end do

			else

				call ListInteger_exception( ERROR, "The current iterator is greater than  the size of the ListInteger", &
						 "Class object ListInteger in the iterate(i) function" )

			end if

		end if

	end subroutine ListInteger_iterate


	!>
	!! @brief Retorna un puntero (referencia) al valor del nodo donde se encuentra el iterador
	!!
	!! @todo Hacer que el procedimineto salga mediante una excepcion.
	!<
	function ListInteger_current( this ) result( output )
		implicit none
		type(ListInteger), target, intent(inout) :: this
		integer, pointer :: output
		
	
		if ( associated(this%current) ) then

			output => this%data( this%iterator )

		else

			call ListInteger_exception(ERROR, " in the ListInteger_current function ", " is empty, current iterator is pointing to null()")

		end if

	end function ListInteger_current


	!>
	!! @brief  Maneja excepciones de la clase
	!<
	subroutine ListInteger_exception(typeMessage, description, debugDescription)
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
	
	end subroutine ListInteger_exception

end module ListInteger_
