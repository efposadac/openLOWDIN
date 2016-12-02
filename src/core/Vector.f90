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
!! @brief Clase encargada de manipular todo lo relacionado con vectores
!!
!! Esta clase manipula todo lo relacionado con matrices de tipo numerico, ademas
!! de servir como una interface transparente para el uso de LAPACK en cuanto a metodos
!! que involuran algebra lineal
!!
!! @author Nestor Aguirre
!!
!! <b> Fecha de creacion : </b> 2008-08-19
!!   - <tt> 2007-08-19 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2007-08-26 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Acoplamiento con el primer nivel de BLAS (daxpy,ddot,dnrm2,idamax)
!!   - <tt> 2008-10-20 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Modifico el formato del la funcion show
!!   - <tt> 2009-04-27 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Corrigio las funciones que retornan el maximo o  minimo de un Vector%
!!   - <tt> 2009-05-27 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Creo operadores y funciones para operaciones basicas entre vectores
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Apadta el Modulo a LOWDIN, no usa Blas
module Vector_
  use CONTROL_
  use Exception_
  implicit none

  interface assignment(=)
     module procedure Vector_copyConstructor
  end interface

  interface operator ( / )
     module procedure Vector_scalarDiv
  end interface


  !< enum Vector_printFormatFlags {
  integer, parameter :: HORIZONTAL = 1
  integer, parameter :: VERTICAL = 2
  integer, parameter :: WITH_KEYS = 4
  !< }
  
  type, public :: Vector
     character(50) :: name
     real(8) , allocatable :: values(:)
  end type Vector

  type, public :: Vector8
     character(50) :: name
     real(8) , allocatable :: values(:)
  end type Vector8

  type, public :: IVector
     integer , allocatable :: values(:)
  end type IVector
  
  type, public :: IVector8
     integer(8) , allocatable :: values(:)
  end type IVector8
 
  
  public :: &
       Vector_constructor, &
       Vector_constructor8, &
       Vector_copyConstructor, &
       Vector_copyConstructor8, &
       Vector_destructor, &
       Vector_destructor8, &
       Vector_show, &
       Vector_writeToFile, &
       Vector_getPtr, &
       Vector_sortElements, &
       Vector_reverseSortElements, &
       Vector_reverseSortElements8, &
       Vector_swapElements, &
       Vector_getSize, &
       Vector_getElement, &
       Vector_setElement, &
       Vector_setIdentity, &
       Vector_setNull, &
       Vector_getMax, &
       Vector_getMin, &
       Vector_isNull, &
       Vector_plus, &
       Vector_scalarDiv, &
       Vector_dot, &
       Vector_cross, &
       Vector_norm, &
       Vector_removeElement, &
       Vector_constructorInteger, &
       Vector_constructorInteger8, &
       Vector_destructorInteger, &
       Vector_swapIntegerElements
  
contains
  
  !>
  !! @brief Constructor por omision
  subroutine Vector_constructor( this, ssize, value, values, name )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(in) :: ssize
    real(8), optional, intent(in) :: value
    real(8), optional, intent(in) :: values(:)
    character(50), optional :: name
    
    real(8) :: valueTmp
    character(50) :: auxName
    
    valueTmp = 0.0_8
    
    if ( allocated( this%values ) ) then
       deallocate( this%values )
       
    end if
    
    allocate( this%values( ssize ) )
    
    auxName = "none"
    
    if( present( name )) then
       
       auxName = trim(name)
       
    end if
    
    if( present(value) ) then
       
       valueTmp = value
       this%values = valueTmp
       
    end if
    
    if( present(values) ) then
       
       this%values = values
       
    end if
    
  end subroutine Vector_constructor

  !>
  !! @brief Constructor por omision
  subroutine Vector_constructor8( this, ssize, value, values, name )
    implicit none
    type(Vector8), intent(inout) :: this
    integer(8), intent(in) :: ssize
    real(8), optional, intent(in) :: value
    real(8), optional, intent(in) :: values(:)
    character(50), optional :: name
    
    real(8) :: valueTmp
    character(50) :: auxName
    
    valueTmp = 0.0_8
    
    if ( allocated( this%values ) ) then
       deallocate( this%values )
       
    end if
    
    allocate( this%values( ssize ) )
    
    auxName = "none"
    
    if( present( name )) then
       
       auxName = trim(name)
       
    end if
    
    if( present(value) ) then
       
       valueTmp = value
       this%values = valueTmp
       
    end if
    
    if( present(values) ) then
       
       this%values = values
       
    end if
    
  end subroutine Vector_constructor8

 

  !>
  !! @brief Constructor por omision
  subroutine Vector_constructorInteger( this, ssize, value, values )
    implicit none
    type(IVector), intent(inout) :: this
    integer, intent(in) :: ssize
    integer, optional, intent(in) :: value
    integer, optional, intent(in) :: values(:)
    
    integer :: valueTmp
    
    valueTmp = 0
    
    if ( allocated( this%values ) ) then
       deallocate( this%values )
       
    end if
    
    allocate( this%values( ssize ) )
    
    if( present(value) ) then
       
       valueTmp = value
       this%values = valueTmp
       
    end if
    
    if( present(values) ) then
       
       this%values = values
       
    end if
    
  end subroutine Vector_constructorInteger

  !>
  !! @brief Constructor por omision
  subroutine Vector_constructorInteger8( this, ssize, value, values )
    implicit none
    type(IVector8), intent(inout) :: this
    integer(8), intent(in) :: ssize
    integer(8), optional, intent(in) :: value
    integer(8), optional, intent(in) :: values(:)
    
    integer :: valueTmp
    
    valueTmp = 0
    
    if ( allocated( this%values ) ) then
       deallocate( this%values )
       
    end if
    
    allocate( this%values( ssize ) )
    
    if( present(value) ) then
       
       valueTmp = value
       this%values = valueTmp
       
    end if
    
    if( present(values) ) then
       
       this%values = values
       
    end if
    
  end subroutine Vector_constructorInteger8
 

  !>
  !! @brief Constructor de copia
  !! Reserva la memoria necesaria para otherMatrix y le asigna los valores de this
  subroutine Vector_copyConstructorInteger( this, otherVector )
    implicit none
    type(IVector), intent(inout) :: this
    type(IVector), intent(in) :: otherVector
    
    if ( allocated( this%values ) ) deallocate( this%values )
    allocate( this%values( size(otherVector%values, DIM=1) ) )
    
    this%values = otherVector%values
    
  end subroutine Vector_copyConstructorInteger


  
  !>
  !! @brief Constructor de copia
  !! Reserva la memoria necesaria para otherMatrix y le asigna los valores de this
  subroutine Vector_copyConstructor( this, otherVector )
    implicit none
    type(Vector), intent(inout) :: this
    type(Vector), intent(in) :: otherVector
    
    if ( allocated( this%values ) ) deallocate( this%values )
    allocate( this%values( size(otherVector%values, DIM=1) ) )
    
    this%values = otherVector%values
    
  end subroutine Vector_copyConstructor

  !>
  !! @brief Constructor de copia
  !! Reserva la memoria necesaria para otherMatrix y le asigna los valores de this
  subroutine Vector_copyConstructor8( this, otherVector )
    implicit none
    type(Vector8), intent(inout) :: this
    type(Vector8), intent(in) :: otherVector
    
    if ( allocated( this%values ) ) deallocate( this%values )
    allocate( this%values( size(otherVector%values, DIM=1) ) )
    
    this%values = otherVector%values
    
  end subroutine Vector_copyConstructor8

  !>
  !! @brief Destructor
  subroutine Vector_destructor( this )
    implicit none
    type(Vector), intent(inout) :: this
    
    if( allocated(this%values) ) deallocate( this%values )
    
  end subroutine Vector_destructor

  !>
  !! @brief Destructor
  subroutine Vector_destructor8( this )
    implicit none
    type(Vector8), intent(inout) :: this
    
    if( allocated(this%values) ) deallocate( this%values )
    
  end subroutine Vector_destructor8



  !>
  !! @brief Destructor
  subroutine Vector_destructorInteger( this )
    implicit none
    type(IVector), intent(inout) :: this
    
    if( allocated(this%values) ) deallocate( this%values )
    
  end subroutine Vector_destructorInteger

  
  !>
  !! @brief Imprime a salida estandar la matriz realizando cambio de linea
  !! con un maximo de "CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS" columnas
  subroutine Vector_show( this, flags, keys )
    implicit none
    type(Vector), intent(in) :: this
    integer, intent(in), optional :: flags
    character(*), intent(in), optional :: keys(:)
    
    integer :: columns
    integer :: ssize
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: u
    integer :: tmpFlags
    character(50) :: formatSize
    
    type(Exception) :: ex
    
    tmpFlags = HORIZONTAL
    if( present(flags) ) then
       if( flags == WITH_KEYS )  then
          tmpFlags = HORIZONTAL + WITH_KEYS
       else
          tmpFlags = flags
       end if
    end if
    
    ssize = size( this%values , DIM=1 )
    
    if( tmpFlags == HORIZONTAL .or. tmpFlags == HORIZONTAL + WITH_KEYS ) then
       
       columns = ssize
       write(formatSize,*) ssize
       
       do k=1, ceiling( (ssize * 1.0)/(CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * 1.0 ) )
          
          l = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k - 1 ) + 1
          u = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k )
          
          if( u > ssize ) then
             columns = l + CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS*( 1 - k ) +  ssize - 1
             u = columns
          end if
          
          if( ( tmpFlags - HORIZONTAL ) == WITH_KEYS ) then
             
             if( present( keys ) ) then
                write (6,"("//trim(formatSize)//"A18)") ( trim(keys(i)), i = l, u )
             else
                write (6,"("//trim(formatSize)//"I15)") ( i, i = l, u )
             end if
             
          end if
          
          print *,""
          write (6,"("//trim(formatSize)//"F15.6)") ( this%values(i), i = l, u )
          print *,""
          
       end do
       
    else if( tmpFlags == VERTICAL .or. tmpFlags == VERTICAL + WITH_KEYS ) then
       
       if( ( tmpFlags - VERTICAL ) == WITH_KEYS ) then
          
          if( present( keys ) ) then
             do i=1, ssize
                
                write (6,"(A18,F15.6)") trim(keys(i)), this%values(i)
             end do
          else
             do i=1, ssize
                write (6,"(I5,F15.6)") i, this%values(i)
             end do
          end if
          
       else
          
          do i=1, ssize
             write (6,"(F15.6)") this%values(i)
          end do
          
       end if
       
    else
       
       call Exception_constructor( ex , WARNING )
       call Exception_setDebugDescription( ex, "Class object Vector in the show() function" )
       call Exception_setDescription( ex, "Bad flags selected" )
       call Exception_show( ex )
       
    end if
    
  end subroutine Vector_show

  !>
  !! @brief Imprime a salida estandar la matriz realizando cambio de linea
  !! con un maximo de "CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS" columnas
  subroutine Vector_showInteger( this, flags, keys )
    implicit none
    type(IVector), intent(in) :: this
    integer, intent(in), optional :: flags
    character(*), intent(in), optional :: keys(:)
    
    integer :: columns
    integer :: ssize
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: u
    integer :: tmpFlags
    character(50) :: formatSize
    
    type(Exception) :: ex
    
    tmpFlags = HORIZONTAL
    if( present(flags) ) then
       if( flags == WITH_KEYS )  then
          tmpFlags = HORIZONTAL + WITH_KEYS
       else
          tmpFlags = flags
       end if
    end if
    
    ssize = size( this%values , DIM=1 )
    
    if( tmpFlags == HORIZONTAL .or. tmpFlags == HORIZONTAL + WITH_KEYS ) then
       
       columns = ssize
       write(formatSize,*) ssize
       
       do k=1, ceiling( (ssize * 1.0)/(CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * 1.0 ) )
          
          l = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k - 1 ) + 1
          u = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k )
          
          if( u > ssize ) then
             columns = l + CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS*( 1 - k ) +  ssize - 1
             u = columns
          end if
          
          if( ( tmpFlags - HORIZONTAL ) == WITH_KEYS ) then
             
             if( present( keys ) ) then
                write (6,"("//trim(formatSize)//"A18)") ( trim(keys(i)), i = l, u )
             else
                write (6,"("//trim(formatSize)//"I15)") ( i, i = l, u )
             end if
             
          end if
          
          print *,""
          write (6,"("//trim(formatSize)//"F15.6)") ( this%values(i), i = l, u )
          print *,""
          
       end do
       
    else if( tmpFlags == VERTICAL .or. tmpFlags == VERTICAL + WITH_KEYS ) then
       
       if( ( tmpFlags - VERTICAL ) == WITH_KEYS ) then
          
          if( present( keys ) ) then
             do i=1, ssize
                
                write (6,"(A18,F15.6)") trim(keys(i)), this%values(i)
             end do
          else
             do i=1, ssize
                write (6,"(I5,F15.6)") i, this%values(i)
             end do
          end if
          
       else
          
          do i=1, ssize
             write (6,"(F15.6)") this%values(i)
          end do
          
       end if
       
    else
       
       call Exception_constructor( ex , WARNING )
       call Exception_setDebugDescription( ex, "Class object Vector in the show() function" )
       call Exception_setDescription( ex, "Bad flags selected" )
       call Exception_show( ex )
       
    end if
    
  end subroutine Vector_showInteger
  
  
  !>
  !! @brief Escribe un vector en el lugar especificado (en binario o con formato)
  subroutine Vector_writeToFile(vvector, unit, file, binary, value, arguments)
    implicit none
    
    type(Vector),optional :: vvector
    integer, optional :: unit
    character(*), optional :: file
    logical, optional :: binary
    real(8), optional :: value
    character(*), optional :: arguments(:)
    
    integer :: elementsNum
    integer :: status
    integer :: n
    
    character(20) :: auxSize
    
    logical :: bbinary
    logical :: existFile
    
    bbinary = .false.
    if(present(binary)) bbinary = binary
    
    
    if ( present( unit ) ) then
       !! It is assumed that the unit y conected to any file (anyways will check)
       inquire(unit=unit, exist=existFile)
       
       if(existFile) then
          
          if(present(arguments)) then
             
             do n = 1, size(arguments)
                
                write(unit) arguments(n)
                
             end do
          end if
             
          if(present(value)) then
             write(unit) 1_8
             write(unit) value
             
          else
             
             write(unit) int(size(vvector%values), 8)
             write(unit) vvector%values
             
          end if

       else
          
          call Vector_exception( ERROR, "Unit file no connected!",&
               "Class object Matrix  in the writeToFile() function" )
             
       end if
       
       
    else if ( present(file) ) then
       if(bbinary) then
          open ( 4,FILE=trim(file),STATUS='REPLACE',ACTION='WRITE', FORM ='UNFORMATTED')
          write(4) int(size(vvector%values), 8)
          write(4) vvector%values
          close(4)
          
       else
          
          open ( 4,FILE=trim(file),STATUS='REPLACE',ACTION='WRITE')
          elementsNum = size( vvector%values )
          write(auxSize,*) elementsNum
          write (4,"("//trim(auxSize)//"ES15.8)") (vvector%values(n), n=1 , elementsNum)
          close(4)
       end if
       
    end if
    
  end subroutine Vector_writeToFile
  
  !>
  !! @brief Obtiene un vector  del lugar especificado
  subroutine Vector_getFromFile(elementsNum, unit, file, binary, value, arguments, output )
    implicit none
    integer, optional, intent(in) :: elementsNum
    integer, optional :: unit
    character(*), optional :: file
    logical, optional :: binary
    character(*), optional :: arguments(:)
    real(8), optional :: value
    type(Vector), optional, intent(out) :: output
    
    type(Exception) :: ex
    character(5000) :: line
    character(20) :: auxSize
    integer :: status
    integer :: n
    integer(8) :: totalSize
    logical :: bbinary
    logical :: existFile
    logical :: found


    if (present(elementsNum)) write(auxSize,*) elementsNum
    
    bbinary = .false.
    existFile = .false.
    
    if(present(binary)) bbinary = binary
    
    if ( present( unit ) ) then
       
       !! check file
       inquire(unit=unit, exist=existFile)
       
       if(existFile) then
          
          rewind(unit)
          
          found = .false.
          line = ""
             
          if(present(arguments)) then

             do                   
                read(unit, iostat = status) line (1:len_trim(arguments(1)))

                if(status == -1) then
                      
                   call vector_exception( ERROR, "End of file!",&
                        "Class object Vector in the getfromFile() function" )
                end if

                if(trim(line) == trim(arguments(1))) then
                   
                   found = .true.                   
                   
                end if
                
                if(found) then
                      
                   backspace(unit)
                   
                   do n = 1, size(arguments)
                      
                      found = .false.
                      read(unit, iostat = status) line
                         
                      if(trim(line) == trim(arguments(n))) then
                            
                         found = .true.
                                                     
                      end if
                         
                   end do
                      
                end if
                   
                if(found) exit
                      
             end do


          end if
             
          !! check size
          read(unit) totalSize
          
          if(present(value)) then
             
             read(unit) value
             
          else
             
             if(totalSize == int(elementsNum,8)) then
                
                if(.not. allocated(output%values)) then                   
                   
                   call Vector_constructor( output, elementsNum )
                   
                end if
                
                read(unit) output%values
                
                ! call Vector_show(output)
                
             else
                
                call Vector_exception( ERROR, "The dimensions of the matrix "//trim(file)//" are wrong ",&
                     "Class object Matrix  in the getFromFile() function"  )
                
             end if

          end if
             
       else

          call Vector_exception( ERROR, "Unit file no connected!",&
               "Class object Matrix  in the getFromFile() function" )
             
       end if

    else if ( present(file) ) then
       
       inquire( FILE = trim(file), EXIST = existFile )
       
       if ( existFile ) then
          
          if(bbinary) then
             open( 4, FILE=trim(file), ACTION='read', FORM='unformatted', STATUS='old' )
             
             read(4) totalSize
             if(totalSize == int(elementsNum,8)) then
                
                call Vector_constructor( output, elementsNum )
                
                read(4) output%values
                close(4)
                
                !! print*, output%values
                
             else
                close(4)
                
                call Vector_exception(ERROR, "The dimensions of the vector in "//trim(file)//" are wrong ", &
                     "Class object Vector_  in the getFromFile() function")
             end if
          else
             
             call Vector_constructor( output, elementsNum )
             open ( 4,FILE=trim(file),STATUS='unknown',ACCESS='sequential' )
             !! verifica el tamahno de las matrice
             read (4,"(A)") line
             if ( len_trim(line) /=elementsNum*15) then
                close(4)
                
                !! call tracebackqq(string="Probando trace ",USER_EXIT_CODE=-1,EPTR=%LOC(ExceptionInfo))
                
                call Vector_exception(ERROR, "The dimensions of the vector in "//trim(file)//" are wrong ", &
                     "Class object Vector_  in the getFromFile() function" )
                
             else
                rewind(4)
             end if
             read ( 4,"("//trim(auxSize)//"ES15.8)") (output%values(n), n=1 ,elementsNum)
             close(4)
             
          end if
       else
          
          call Vector_exception(ERROR, "The file "//trim(file)//" don't exist " , &
               "Class object Vector_  in the getFromFile() function" )
          
       end if
       
    end if
    
  end subroutine Vector_getFromFile
  
  !>
  !! @brief Devuelve un apuntador a la matrix solicitada
  !! @param this matrix de m x n
  !! @return Apuntador a la matriz solicitada.
  !! @todo No ha sido probada
  function Vector_getPtr( this ) result( output )
    implicit none
    type(Vector) , target , intent(in) :: this
    real(8) , pointer :: output(:)
    
    output => null()
    output => this%values
    
  end function Vector_getPtr
  
  subroutine Vector_sortElements(this, factor)
    type(Vector) :: this
    
    integer i,j,n
    integer, optional :: factor
    
    n = Vector_getSize(this)
    if ( .not. present (factor) ) then
      do i=1,n
         do j=i+1,n
            if (this%values(j).gt.this%values(i)) then
               call Vector_swapElements( this, i, j )
            end if
         end do
      end do
    else 
      factor = 0
      do i=1,n
         do j=i+1,n
            if (this%values(j).gt.this%values(i)) then
               factor = factor + 1
               call Vector_swapElements( this, i, j )
            end if
         end do
      end do
    end if

  end subroutine Vector_sortElements

  subroutine Vector_reverseSortElements(this,indexVector,m)
    type(Vector) :: this
    type(IVector), optional :: indexVector
    integer, optional :: m
    integer i,j,n
    
    n = Vector_getSize(this)
    if ( .not. present (indexVector) ) then
      do i=1,n
         do j=i+1,n
            if (this%values(j).lt.this%values(i)) then
               call Vector_swapElements( this, i, j )
            end if
         end do
      end do
    else
    
      if ( .not. present (m) ) then

        do i=1,n
          indexVector%values(i) = i
        end do 

        do i=1,n
           do j=i+1,n
              if (this%values(j).lt.this%values(i)) then
                 call Vector_swapElements( this, i, j )
                 call Vector_swapIntegerElements( indexVector, i, j )
              end if
           end do
        end do
      else

        do i=1,n
          indexVector%values(i) = i
        end do 

        do i=1,m
           do j=i+1,n
              if (this%values(j).lt.this%values(i)) then
                 call Vector_swapElements( this, i, j )
                 call Vector_swapIntegerElements( indexVector, i, j )
              end if
           end do
        end do
      end if
    end if

  end subroutine Vector_reverseSortElements

  subroutine Vector_reverseSortElements8(this,indexVector,m)
    type(Vector8) :: this
    type(IVector8), optional :: indexVector
    integer(8), optional :: m
    integer(8) i,j,n
    
    n = Vector_getSize8(this)
    if ( .not. present (indexVector) ) then
      do i=1,n
         do j=i+1,n
            if (this%values(j).lt.this%values(i)) then
               call Vector_swapElements8( this, i, j )
            end if
         end do
      end do
    else
    
      if ( .not. present (m) ) then

        do i=1,n
          indexVector%values(i) = i
        end do 

        do i=1,n
           do j=i+1,n
              if (this%values(j).lt.this%values(i)) then
                 call Vector_swapElements8( this, i, j )
                 call Vector_swapIntegerElements8( indexVector, i, j )
              end if
           end do
        end do
      else

        do i=1,n
          indexVector%values(i) = i
        end do 

        do i=1,m
           do j=i+1,n
              if (this%values(j).lt.this%values(i)) then
                 call Vector_swapElements8( this, i, j )
                 call Vector_swapIntegerElements8( indexVector, i, j )
              end if
           end do
        end do
      end if
    end if

  end subroutine Vector_reverseSortElements8
  
  !>
  !! @brief Intercambia los elementos i y j el vector
  subroutine Vector_swapElements( this, i, j )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    
    real(8) :: value1
    real(8) :: value2
    
    value1 = this%values( i )
    value2 = this%values( j )
    
    this%values( i ) = value2
    this%values( j ) = value1
    
  end subroutine Vector_swapElements

  !>
  !! @brief Intercambia los elementos i y j el vector
  subroutine Vector_swapElements8( this, i, j )
    implicit none
    type(Vector8), intent(inout) :: this
    integer(8), intent(in) :: i
    integer(8), intent(in) :: j
    
    real(8) :: value1
    real(8) :: value2
    
    value1 = this%values( i )
    value2 = this%values( j )
    
    this%values( i ) = value2
    this%values( j ) = value1
    
  end subroutine Vector_swapElements8

  !>
  !! @brief Intercambia los elementos i y j el vector
  subroutine Vector_swapIntegerElements( this, i, j )
    implicit none
    type(IVector), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    
    integer :: value1
    integer :: value2
    
    value1 = this%values( i )
    value2 = this%values( j )
    
    this%values( i ) = value2
    this%values( j ) = value1
    
  end subroutine Vector_swapIntegerElements


  !>
  !! @brief Intercambia los elementos i y j el vector
  subroutine Vector_swapIntegerElements8( this, i, j )
    implicit none
    type(IVector8), intent(inout) :: this
    integer(8), intent(in) :: i
    integer(8), intent(in) :: j
    
    integer(8) :: value1
    integer(8) :: value2
    
    value1 = this%values( i )
    value2 = this%values( j )
    
    this%values( i ) = value2
    this%values( j ) = value1
    
  end subroutine Vector_swapIntegerElements8


  
  !>
  !! @brief Retorna el tamano del vector
  function Vector_getSize( this ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    integer :: output
    
    output = size( this%values , DIM=1 )
    
  end function Vector_getSize
  
  !>
  !! @brief Retorna el tamano del vector
  function Vector_getSize8( this ) result ( output )
    implicit none
    type(Vector8), intent(inout) :: this
    integer(8) :: output
    
    output = size( this%values , DIM=1 )
    
  end function Vector_getSize8
  
  !>
  !! @brief Retorna el elemento i-esimo del vector
  function Vector_getElement( this, i ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(in) :: i
    real(8) :: output
    
    output = this%values( i )
    
  end function Vector_getElement
  
  !>
  !! @brief Selecciona el valor del elemento i-esimo del vector
  subroutine Vector_setElement( this, i, value )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(in) :: i
    real(8), intent(in) :: value
    
    this%values( i ) = value
    
  end subroutine Vector_setElement

  !>
  !! @brief Selecciona todos los elementos del vector a 1.0
  subroutine Vector_setIdentity( this )
    implicit none
    type(Vector), intent(inout) :: this
    
    this%values = 1.0_8
    
  end subroutine Vector_setIdentity

  !>
  !! @brief Selecciona a cero todos los elementos del vector
  subroutine Vector_setNull( this )
    implicit none
    type(Vector), intent(inout) :: this
    
    this%values = 0.0_8
    
  end subroutine Vector_setNull
  
  !>
  !! @brief Retorna el maximo valor absoluto y la posicion de este
  function Vector_getMax( this, pos ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(inout), optional :: pos
    real(8) :: output
    
    integer :: posTmp
    
    if( present(pos) ) then
       pos = maxloc( this%values, DIM=1 )
    else
       posTmp = maxloc( this%values, DIM=1 )
    end if
    
    output = this%values(pos)
    
  end function Vector_getMax
  
  !>
  !! @brief Retorna el valor minimo valor absoluto y la posicion de este
  function Vector_getMin( this, pos ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    integer, intent(inout), optional :: pos
    real(8) :: output
    
    integer :: posTmp
    
    if( present(pos) ) then
       pos = minloc( this%values, DIM=1 )
    else
       posTmp = minloc( this%values, DIM=1 )
    end if
    
    output = this%values(pos)
    
  end function Vector_getMin

  !>
  !! @brief Retorna true si todos lo elementos del vector sun iguales a
  !! cero, false de otra manera
  !! @todo Falta implementar, depende de la declaracion de un umbral de comparacion de enteros en la clase APMO
  function Vector_isNull( this ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    logical :: output
    
    output = .false.
    
  end function Vector_isNull
  
  !>
  !! @brief Suma dos vectores
  !! @warning No verifica que los dos vectores sean del mismo tamano
  function Vector_plus( this, otherVector ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    type(Vector), intent(in) :: otherVector
    type(Vector) :: output
    
    call Vector_copyConstructor( output, otherVector )
    output%values = this%values + otherVector%values
    
  end function Vector_plus
  
  function Vector_scalarDiv( this, scalar ) result( output )
    implicit none
    type(vector), intent(in) :: this
    real(8), intent(in) :: scalar
    type(Vector) :: output
    
    output=this
    output%values=output%values/scalar
    
  end function Vector_scalarDiv

  !>
  !! @brief Calcula el producto punto de dos vectores
  !! @warning No verifica que los dos vectores sean del mismo tamano
  function Vector_dot(this, otherVector) result(output)
    implicit none
    type(Vector), intent(inout) :: this
    type(Vector), intent(in) :: otherVector
    real(8) :: output
    
    output = dot_product( this%values, otherVector%values )
    
  end function Vector_dot
  
  !>
  !! @brief Calcula el producto vectorial tridimedional
  !! @warning No verifica que los dos vectores sean del mismo tamano
  !! y asume que se trata de vectores tridimesionales
  function Vector_cross( this, otherVector ) result ( output )
    implicit none
    type(Vector), intent(in) :: this
    type(Vector), intent(in) :: otherVector
    type(Vector) :: output
    
    allocate( output%values(3) )
    
    output%values = [this%values(2)*otherVector%values(3) -  this%values(3)*otherVector%values(2), &
         this%values(3)*otherVector%values(1) -  this%values(1)*otherVector%values(3), &
         this%values(1)*otherVector%values(2) -  this%values(2)*otherVector%values(1)  ]
    
  end function Vector_cross
  
  !>
  !! @brief Calcula el producto vectorial tridimedional entre vectores nativos de fortran
  function Vector_Fortran_Cross( fotranVector, otherFotranVector ) result ( output )
    implicit none
    real(8), intent(in) :: fotranVector(3)
    real(8), intent(in) :: otherFotranVector(3)
    real(8) :: output(3)
    
    
    
    output = [fotranVector(2)*otherFotranVector(3) -  fotranVector(3)*otherFotranVector(2), &
         fotranVector(3)*otherFotranVector(1) -  fotranVector(1)*otherFotranVector(3), &
         fotranVector(1)*otherFotranVector(2) -  fotranVector(2)*otherFotranVector(1)  ]
    
  end function Vector_Fortran_Cross
  
  
  !>
  !! @brief Calcula la norma euclideana de un vector
  function Vector_norm( this ) result ( output )
    implicit none
    type(Vector), intent(inout) :: this
    real(8) :: output
    
    output = sqrt( dot_product( this%values, this%values ) )
    
  end function Vector_norm
  
  !>
  !! @brief Remueve la fila especificada de una matriz
  subroutine Vector_removeElement( this, numberOfElement )
    implicit none
    type(Vector) :: this
    integer, intent(in) :: numberOfElement
    
    real(8), allocatable :: auxArray(:)
    integer :: numberOfElements
    
    numberOfElements = size( this%values )
    
    if (numberOfElement <= numberOfElements ) then
       
       allocate( auxArray(numberOfElements-1) )
       auxArray(1:numberOfElement-1) = this%values(1:numberOfElement-1)
       auxArray(numberOfElement:numberOfElements-1) = this%values(numberOfElement+1:numberOfElements)
       deallocate( this%values )
       allocate( this%values(numberOfElements-1) )
       this%values = auxArray
       deallocate( auxArray )
       
    end if
    
  end subroutine Vector_removeElement

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine Vector_exception( typeMessage, description, debugDescription)
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
    
  end subroutine Vector_exception
  
end module Vector_
