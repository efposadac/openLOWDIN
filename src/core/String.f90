!!*****************************************************************************
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
!!*****************************************************************************
!>
!! @brief  Modulo para manejo de cadenas de caracateres
!!
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creaciacion : </b> 2008-18-10
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2008-18-10 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
module String_
  use Exception_
  use Vector_
  implicit none
  
  public :: &
       String_getLowercase, &
       String_getUppercase, &
       String_findSubstring
  
contains

  !>
  !! @brief Devuelve una cadena de entrada en minuscula
  function String_getLowercase( inputString, from, to ) result( output )
    implicit none
    character(*):: inputString
    integer, optional :: from
    integer, optional :: to
    character(255):: output
    
    integer:: ffrom
    integer:: tto
    integer:: i
    
    ffrom = 1
    tto =  len_trim(inputString)
    
    if(present(from)) ffrom = from
    if(present(to)) tto = to
    
    output=trim(inputString)
    
    do i=ffrom, tto
       output(i:i) = char( IOR( ichar(inputString (i:i) ),32) )
    end do
    
  end function String_getLowercase
  
  !>
  !! @brief Devuelve una cadena de entrada en mayuscula
  function String_getUppercase( inputString, from, to ) result( output )
    implicit none
    character(*):: inputString
    integer, optional :: from
    integer, optional :: to
    character(255):: output
    
    integer:: ffrom
    integer:: tto
    integer:: i
    integer::j
    
    ffrom = 1
    tto =  len_trim(inputString)
    
    if(present(from)) ffrom = from
    if(present(to)) tto = to
    
    output=trim(inputString)
    do i=ffrom, tto
       j=ichar(inputString (i:i))
       if ( j>96 ) output(i:i) = char( XOR( j,32) )
       
    end do
    
  end function String_getUppercase
  
  
  !>
  !! @brief encuentra la posicion de primera coincidencia de una cadena especificada
  !!        e indica el posible numero de coincidencias en la cadena.
  function String_findSubstring( inputString, substring, probabilityOfCoincidence ) result(output)
    implicit none
    character(*):: inputString
    character(*):: substring
    integer, optional, intent(inout) :: probabilityOfCoincidence
    integer:: output
    
    integer:: i
    integer::j
    integer::k
    integer, allocatable :: mark(:)
    logical:: areCoincidence
    integer:: probOfCoincidence
    
    output=0
    k=len(inputString)
    
    !! Separa memoria para el maximo numero de posible coincidencias entre cadenas
    allocate( mark( k ) )
    mark = 0
    
    
    !!*************************************************************
    !! Busca los indices de las posibles posiciones de coincidencia
    !!****
    i=0
    do j=1,len(inputString)
       
       if( inputString(j:j) == substring(1:1) ) then
          i=i+1
          mark(i) = j
       end if
    end do
    
    !!*************************************************************
    probOfCoincidence=i
    if (present(probabilityOfCoincidence)) probabilityOfCoincidence=probOfCoincidence
    
    !!*************************************************************
    !!  Rechaza los posible indice que no podrian servir
    !!*******
    !????
    
    !!*************************************************************
    !! inicia la busqueda de la primera coincidencias
    !!*******
    do  i=1, probOfCoincidence
       
       areCoincidence=.true.
       do j=2,len(substring)
          k=mark(i)+j-1
          if( inputString(k:k)/=substring(j:j) ) then
             areCoincidence=.false.
             exit
          end if
       end do
       
       if(areCoincidence) exit
       
    end do
    
    !!*************************************************************
    
    output=mark(i)
    
    deallocate(mark)
    
  end function String_findSubstring
  
  !>
  !! @brief Convierte un vector de reales de doble precision en una cadena de caracteres
  function String_convertVectorOfRealsToString( vectorOfReals ) result( output )
    implicit none
    type(Vector), intent(in) :: vectorOfReals
    character(20), allocatable :: output(:)
    
    integer:: orderOfMatrix
    integer:: i
    
    orderOfMatrix=size(vectorOfReals%values )
    
    if ( allocated( output ) ) deallocate( output )
    allocate( output( orderOfMatrix ) )
    
    do i=1, orderOfMatrix

       write (output(i)," (F12.4)") vectorOfReals%values(i)
       
    end do
    
  end function String_convertVectorOfRealsToString

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine String_exception( typeMessage, description, debugDescription)
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
    
  end subroutine String_exception
  
end module String_
