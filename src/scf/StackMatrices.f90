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
!! @brief Modulo para implementacion de pilas
!! Modulo para manejo de pilas de datos, esta implementacion esta basada en
!! una estructura de datos tipo lista enlazada.
!!
!! @author Sergio A. Gonzalez
!!
!! <b> Fecha de creacion : </b> 2008-09-19
!!   - <tt> 2007-08-19 </tt>: Sergio A. Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# adapta el modulo para su inclusion en Lowdin
!!
!! @todo verificar la liberacion de memoria al eliminar nodos de la pila
module StackMatrices_
  use Matrix_
  use Exception_
  implicit none

  type, public :: StackMatrices
     character(30) :: name
     type(Matrix), allocatable :: node(:)
     real(8), pointer :: current(:,:)
     integer :: iterator
     integer :: maxSize
     logical :: isIstanced
  end type StackMatrices

  
  public :: &
       StackMatrices_constructor, &
       StackMatrices_destructor, &
       StackMatrices_isInstanced, &
       StackMatrices_size, &
       StackMatrices_empty, &
       StackMatrices_push, &
       StackMatrices_pop, &
       StackMatrices_popFront  
  
contains
  
  !>
  !! @brief Define el constructor para la clase
  subroutine StackMatrices_constructor( this, name, ssize )
    implicit none
    type(StackMatrices), intent(inout) :: this
    character(*),optional :: name
    integer, optional, intent(in) :: ssize
    
    this%maxSize =10
    this%name = "undefined"
    
    if ( present(name) ) this%name = trim(name)
    if (present(ssize)) this%maxSize =ssize
    this%isIstanced =.true.
    
    
  end  subroutine StackMatrices_constructor
  
  !>
  !! @brief Define el destructor para la clase
  !! @ Verificar que la memoria priviamente separada sea liberada
  subroutine StackMatrices_destructor( this )
    implicit none
    type(StackMatrices), intent(inout) :: this
    
    integer i
    
    if (allocated(this%node)) then
       
       this%current => null()
       
       do i=1,size(this%node)
          
          call Matrix_destructor( this%node(i) )
          
       end do
       
       deallocate( this%node )
       this%name=""
       this%iterator = 0
       this%maxSize =0
       this%isIstanced =.false.
    end if
    
  end  subroutine StackMatrices_destructor
  
  !>
  !! @brief Indica si la pila esta instanciadad
  !>
  function StackMatrices_isInstanced( this) result(output)
    implicit none
    type(StackMatrices), intent(in) :: this
    logical :: output
    
    output = this%isIstanced
    
  end function StackMatrices_isInstanced
  
  !>
  !! @brief Indica si la pila tiene o no nodos
  function StackMatrices_empty( this) result(output)
    implicit none
    type(StackMatrices), intent(in) :: this
    logical :: output
    
    output=.false.
    
    if( .not.allocated(this%node) .or. this%iterator == 0 ) output =.true.
    
  end function StackMatrices_empty
  
  !>
  !! @brief Retorna el tamano actual de la pila
  function StackMatrices_size( this) result(output)
    implicit none
    type(StackMatrices),  intent(in) :: this
    integer :: output
    
    output=0
    
    if( allocated(this%node) ) output =size(this%node)
    
  end function StackMatrices_size
  
  
  !>
  !! @brief Adiciona un nodeo a la pila
  subroutine StackMatrices_push( this, node )
    implicit none
    type(StackMatrices), intent(inout), target :: this
    real(8),intent(in) :: node(:,:)
    
    type(Exception) :: ex
    type(Matrix), allocatable ::  auxStack(:)
    integer :: currentSize
    integer :: columns
    integer :: rows
    integer :: i
    
    rows =size( node, 1)
    columns =size( node, 2)
    
    if ( allocated(this%node) ) then
       
       currentSize = size(this%node)
       if(  currentSize < this%maxSize ) then
          
          allocate( auxStack( currentSize ) )
          this%current => null()
          
          do i=1, currentSize
             auxStack(i)=  this%node(i)
             call Matrix_destructor( this%node(i) )
          end do
          
          deallocate( this%node )
          allocate( this%node(currentSize+1) )
          this%iterator = currentSize+1
          
          do i=1, currentSize
             this%node(i) = auxStack(i)
             call Matrix_destructor( auxStack(i) )
          end do
          
          deallocate( auxStack )
          call Matrix_constructor( this%node(currentSize + 1), int(rows,8), int(columns,8) )
          this%node(currentSize + 1)%values=node
          this%current => this%node(currentSize + 1)%values
          
       else
          
          allocate( auxStack( this%maxSize ) )
          this%current => null()
          
          do i=1, this%maxSize
             auxStack(i)=  this%node(i)
             call Matrix_destructor( this%node(i) )
          end do
          
          do i=2, this%maxSize
             this%node(i-1) = auxStack(i)
             call Matrix_destructor( auxStack(i) )
          end do
          
          call Matrix_destructor( auxStack(1) )
          
          deallocate( auxStack )
          call Matrix_constructor( this%node(this%maxSize), int(rows,8), int(columns,8) )
          this%node(this%maxSize)%values=node
          this%current => this%node(this%maxSize )%values
          
       end if
       
    else
       
       allocate( this%node(1) )
       call Matrix_constructor( this%node(1), int(rows,8), int(columns,8) )
       this%node(1)%values=node
       this%current => this%node(1)%values
       this%iterator = 1
       
    end if
    
  end subroutine StackMatrices_push
  
  !>
  !! @brief Elimina el utimo nodeo de la pila
  subroutine StackMatrices_pop(this)
    implicit none
    type(StackMatrices), intent(inout) :: this
    
    type(Exception) :: ex
    
  end subroutine StackMatrices_pop
  
  !>
  !! @brief Elimina el primer elemento de la pila
  subroutine StackMatrices_popFront( this )
    implicit none
    type(StackMatrices), intent(inout), target :: this
    
    type(Matrix), allocatable ::  auxStack(:)
    integer :: currentSize
    integer :: i
    
    currentSize = size(this%node)
    
    if(  currentSize >= 2 ) then
       
       allocate( auxStack( currentSize-1 ) )
       this%current => null()
       
       do i=2, currentSize
          auxStack(i-1)=  this%node(i)
          call Matrix_destructor( this%node(i) )
       end do
       
       call Matrix_destructor( this%node(1) )
       
       deallocate( this%node )
       allocate( this%node(currentSize-1) )
       this%iterator = currentSize-1
       currentSize=this%iterator
       
       do i=1, currentSize
          this%node(i) = auxStack(i)
          call Matrix_destructor( auxStack(i) )
       end do
       
       deallocate( auxStack )
       
       this%current => this%node(currentSize)%values
       
    end if
    
  end subroutine StackMatrices_popFront
  
  !>
  !! @brief Elimina el ultimo nodeo de la pila
  function StackMatrices_getValuesPtr( this, iterator ) result(output)
    implicit none
    type(StackMatrices), intent(inout), target :: this
    integer, intent(in) :: iterator
    real(8), pointer :: output(:,:)
    
    output => null()
    output => this%node(iterator)%values
    
  end function StackMatrices_getValuesPtr
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine StackMatrices_exception( typeMessage, description, debugDescription)
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
    
  end subroutine StackMatrices_exception
  
end module StackMatrices_
