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
!! @brief  Modulo para implementacion de listas de reales, esta implemtacion es temporal
!!		en remplazo a la depuracion de la lista enlzadas la cual presenta problemas
!!		con el uso de memoria. Ver archivo list.f90.bkp. Al  completaarse
!! 		tamano de la lista, esta se comporta como una pila -adicion realizada por conveniencia-!!!!!!
!!
!! Modulo para manejo de pilas de datos, esta implementacion esta basada en
!! una estructura de datos tipo lista enlazada.
!!
!! @author Sergio A. Gonzalez
!!
!! <b> Fecha de creacion : </b> 2008-09-19
!!   - <tt> 2007-08-19 </tt>: Sergio A. Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
!!
!! @todo verificar la liberacion de memoria al eliminar elementos de la pila
module List_
  use Exception_
  implicit none

  
  type, public :: List
     
     character(30) :: name
     real(8), allocatable :: data(:)
     real(8), pointer :: current
     integer :: iterator
     integer :: maxSize
     logical :: isSizeUndefined
     
  end type List
  
  
  public :: &
       List_constructor, &
       List_assign, &
       List_back, &
       List_begin, &
       List_clear, &
       List_empty, &
       List_end, &
       List_erase, &
       List_front, &
       List_insert, &
       List_max_size, &
       List_merge, &
       List_pop_back, &
       List_pop_front, &
       List_remove, &
       List_push_back, &
       List_push_front, &
       List_resize, &
       List_reverse, &
       List_size, &
       List_sort, &
       List_swap , &
       List_unique, &
       List_iterate, &
       List_current
  
  
contains
  
  !>
  !! @brief Crea una lista y la inicializa
  subroutine List_constructor( this, name, ssize )
    implicit none
    type(List), target, intent(inout) :: this
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
    this%data = 0.0_8
    this%iterator = 0
    this%current => null()
    
  end  subroutine List_constructor
  
  !>
  !! @brief Elimina una lista previamente creada y libera la memoria asociada
  !! @todo Verificar que la memoria previamente separada sea liberada
  subroutine List_destructor( this )
    implicit none
    type(List), intent(inout) :: this
    
    integer :: i
    integer :: currentSize
    
    if ( allocated(this%data) ) then
       
       this%current => null()
       this%maxSize = 0
       this%iterator = 0
       deallocate(this%data)
       
    end if
    
  end  subroutine List_destructor
  
  !>
  !! @brief Asigna un grupo elementos a la lista
  !! @todo Falta implementar
  subroutine List_assign(this, first, last, value)
    implicit none
    type(List), intent(inout) :: this
    integer :: first
    integer :: last
    real(8) :: value
    
  end subroutine List_assign
  
  !>
  !! @brief Retorna un puntero al ultimo elemento de la lista
  function List_back( this ) result(output)
    implicit none
    type(List), target, intent(inout) :: this
    real(8) ,pointer :: output
    
    type(Exception) :: ex
    
    if ( allocated( this%data) ) then
       
       output => null()
       this%iterator = size( this%data )
       output => this%data(this%iterator)
       
       
    else
       call List_exception( ERROR, "The list "// trim(this%name) //" is empty" , &
            "Class object List in the back() function" )
    end if
    
  end function List_back
  
  !>
  !! @brief Coloca el puntero del iterador al inicio de la lista
  subroutine List_begin( this)
    type(List), target, intent(inout) :: this
    
    type(Exception) :: ex
    
    if ( allocated(this%data) ) then
       
       this%current => null()
       this%current => this%data(1)
       this%iterator = 1
       
    else
       call List_exception( ERROR, "The list "// trim(this%name) //" is empty" , &
            "Class object List in the begin() function")
       
    end if
    
  end subroutine List_begin
  
  !>
  !! @brief Remueve todo los elementos del la lista
  subroutine List_clear( this)
    type(List), intent(inout) :: this
    
    if ( allocated(this%data) ) then
       
       this%current=> null()
       deallocate( this%data )
       allocate( this%data(1) )
       this%data = 0.0_8
       this%iterator = 0
       
       
    end if
    
  end subroutine List_clear
  
  !>
  !! @brief Indica si la lista tiene o no elementos
  function List_empty( this ) result(output)
    implicit none
    type(List), intent(in) :: this
    logical :: output
    
    output=.false.
    
    if( .not.allocated(this%data) .or. this%iterator == 0 ) output =.true.
    
  end function List_empty
  
  !>
  !! @brief Coloca el puntero del iterador al final de la lista
  subroutine List_end( this)
    type(List), target, intent(inout) :: this
    
    type(Exception) :: ex
    
    if ( allocated( this%data ) ) then
       
       this%current=>null()
       this%iterator = size(this%data)
       this%current => this%data( this%iterator )
       
       
    else
       call List_exception( ERROR, "The list "// trim(this%name) //" is empty" , &
            "Class object List in the end() function")
       
    end if
    
  end subroutine List_end

  !>
  !! @brief Remueve todo los elementos del la lista
  subroutine List_erase( this)
    type(List), intent(inout) :: this
    
    type(Exception) :: ex
    
    
  end subroutine List_erase
  
  
  !>
  !! @brief Retorna una puntero al primer elemento de la lista
  function List_front( this ) result(output)
    implicit none
    type(List), target, intent(inout) :: this
    real(8),pointer :: output
    
    type(Exception) :: ex
    
    if ( allocated( this%data ) ) then
       
       output => null()
       output => this%data(1)
       
    else
       call List_exception( ERROR, "The list "// trim(this%name) //" is empty" , &
            "Class object List in the front() function")
       
    end if
    
  end function List_front
  
  !>
  !! @brief Inserta un elemento en la posicion especificada
  !! @todo Falta implementacion completa
  subroutine List_insert( this, position)
    type(List), intent(inout) :: this
    integer :: position
    
  end subroutine List_insert
  
  !>
  !! @brief Retorna el maximo numero de elementos que se puden almacenar
  !! 		en la lista.
  function List_max_size( this) result(output)
    implicit none
    type(List),  intent(in) :: this
    integer :: output
    
    output = this%maxSize
    
  end function List_max_size
  
  !>
  !! @brief Fusiona dos lista en una borrando un de ellas
  !! @todo Falta implementacion completa
  subroutine List_merge( this, otherThisToMerge )
    type(List), intent(inout) :: this
    type(List), intent(inout) :: otherThisToMerge
    
    
  end subroutine List_merge
  
  !>
  !! @brief Elimina el ultimo elemento de la lista, el primer elemento no es eliminado
  subroutine List_pop_back( this )
    implicit none
    type(List), target, intent(inout) :: this
    
    type(Exception) :: ex
    real(8), allocatable :: auxData(:)
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
       call List_exception(WARNING, "The list is empty" , &
            "Class object List in the pop_back() function")
    end if
    
    
  end subroutine List_pop_back
  
  !>
  !! @brief Elimina el primer elemento de la lista
  !! @todo Falta toda la implemtacion
  subroutine List_pop_front( this )
    implicit none
    type(List), intent(inout) :: this
    
  end subroutine List_pop_front
  
  !>
  !! @brief Adiciona un elemento al final a la lista
  !! @ Verificar que la memoria de cada nodo sea liberada cuando este se elimine
  subroutine List_push_back( this, data )
    implicit none
    type(List), target, intent(inout) :: this
    real(8),intent(in) :: data
    
    real(8), allocatable :: auxData(:)
    integer :: ssize
    type(Exception) :: ex
    
    if ( allocated(this%data)  ) then
       
       ssize = size( this%data )
       
       if ( this%isSizeUndefined  ) this%maxSize = this%maxSize + 1
       
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
       call List_exception(ERROR, "The list "// trim(this%name) //" is empty" , &
            "Class object List in the push back() function")
       
    end if
    
  end subroutine List_push_back
  
  !>
  !! @brief Adiciona un elemento al inicio de la lista
  !! @todo Falta toda la implemtacion
  subroutine List_push_front( this )
    implicit none
    type(List), intent(inout) :: this
    
  end subroutine List_push_front
  
  !>
  !! @brief Remueve el nodo especificado de la lista
  !! @todo Falta toda la implemtacion
  subroutine List_remove( this, position )
    implicit none
    type(List), intent(inout) :: this
    integer :: position
    
  end subroutine List_remove
  
  !>
  !! @brief Redefine el tamano maximo de la lista
  subroutine List_resize( this, resize )
    implicit none
    type(List), intent(inout) :: this
    integer :: resize
    
    this%maxSize = resize
    
  end subroutine List_resize
  
  !>
  !! @brief Invierte la lista
  !! @todo Falta toda la implemtacion
  subroutine List_reverse( this)
    implicit none
    type(List), intent(inout) :: this
    
  end subroutine List_reverse
  
  !>
  !! @brief Retorna el tamano  actual o ocupacion de la pila
  function List_size( this) result(output)
    implicit none
    type(List),  intent(in) :: this
    integer :: output
    
    output=0
    
    if( allocated(this%data) ) then
       
       output =size( this%data )
       
    end if
    
  end function List_size
  
  !>
  !! @brief Ordena los elementos de la lista en orden ascendente
  !! @todo Falta toda la implemtacion
  subroutine List_sort( this)
    implicit none
    type(List), intent(inout) :: this
    
    
  end subroutine List_sort
  
  !>
  !! @brief Intercambia los elementos entre dos listas
  !! @todo Falta toda la implemtacion
  subroutine List_swap( this, otherThis)
    implicit none
    type(List), intent(inout) :: this
    type(List), intent(inout) :: otherThis
    
    
  end subroutine List_swap
  
  !>
  !! @brief remueve los elemetos repetidos de la lista siempre que
  !!		sean consecutivos
  !! @todo Falta toda la implemtacion
  subroutine List_unique( this)
    implicit none
    type(List), intent(inout) :: this
    
    
  end subroutine List_unique
  
  !>
  !! @brief Mueve el iterador de la lista el numero de nodo que se especifique
  subroutine List_iterate( this, iterator)
    implicit none
    type(List), target, intent(inout) :: this
    integer, optional :: iterator
    
    type(Exception) :: ex
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
          call List_exception(ERROR, "The current iterator is greater than  the size of the List", &
               "Class object List in the iterate(i) function")
          
       end if
       
    end if
    
  end subroutine List_iterate
  
  !>
  !! @brief Retorna un puntero (referencia) al valor del nodo donde se encuentra el iterador
  !! @todo Hacer que el procedimineto salga mediante una excepcion.
  function List_current( this ) result( output )
    implicit none
    type(List), target, intent(inout) :: this
    real(8), pointer :: output
    
    
    if ( associated(this%current) ) then
       
       output => this%data( this%iterator )
       
    else
       
       call List_exception(ERROR, "in the List_current function", "is empty, current iterator is pointing to null()")
       
    end if
    
  end function List_current
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine List_exception(typeMessage, description, debugDescription)
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
    
  end subroutine List_exception
  
end module List_
