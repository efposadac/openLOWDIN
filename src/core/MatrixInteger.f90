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

module MatrixInteger_
  use CONTROL_
  use Matrix_
  use Exception_
  use Vector_
  use LapackInterface_
  implicit none

	!>
	!! @brief Clase encargada de manipular matrices de enteros
	!!
	!! Esta clase manipula todo lo relacionado con matrices de tipo numerico
	!!
	!! @author Sergio Gonzalez
	!!
	!! <b> Fecha de creacion : </b> 2009-04-24
	!!   - <tt> 2009-04-24 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
	!!        -# Creacion del archivo y las funciones basicas
	!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
	!<

	type, public :: MatrixInteger
		integer , allocatable :: values(:,:)
	end type MatrixInteger

	interface assignment(=)
		module procedure MatrixInteger_copyConstructor
	end interface


	public :: &
		MatrixInteger_constructor, &
		MatrixInteger_copyConstructor, &
		! MatrixInteger_diagonalConstructor, &
		! MatrixInteger_randomElementsConstructor, &
		MatrixInteger_destructor, &
		! MatrixInteger_show, &
		! MatrixInteger_getPtr, &
		! MatrixInteger_swapRows, &
		! MatrixInteger_swapColumns, &
		! MatrixInteger_getNumberOfRows, &
		! MatrixInteger_getNumberOfColumns, &
		! MatrixInteger_getElement, &
		! MatrixInteger_setElement, &
		! MatrixInteger_getColumn, &
		! MatrixInteger_getRow, &
		! MatrixInteger_setIdentity, &
		! MatrixInteger_setNull, &
		! MatrixInteger_getMax, &
		! MatrixInteger_getMin, &
		! MatrixInteger_getDeterminant, &
		! MatrixInteger_isNull, &
		! MatrixInteger_getTransPose, &
		! MatrixInteger_trace, &
		! MatrixInteger_plus, &
		MatrixInteger_removeRow
		! MatrixInteger_removeColumn

	private :: &
		MatrixInteger_exception
		
contains

	!>
	!! @brief Constructor
	!! Constructor por omisi�n
	!<
	subroutine MatrixInteger_constructor( this, rows, cols, value )
		implicit none
		type(MatrixInteger), intent(inout) :: this
		integer, intent(in) :: rows
		integer, intent(in) :: cols
		integer, optional, intent(in) :: value

		integer :: valueTmp

		valueTmp = 0
		if( present(value) ) valueTmp = value

		if (allocated(this%values)) deallocate(this%values)
		allocate( this%values( rows, cols ) )

		this%values = valueTmp

	end subroutine MatrixInteger_constructor

	!>
	!! @brief Constructor de copia
	!! Reserva la memoria necesaria para this y le asigna los valores de otherMatrixInteger
	!<
	subroutine MatrixInteger_copyConstructor( this, otherMatrixInteger )
		implicit none
		type(MatrixInteger), intent(inout) :: this
		type(MatrixInteger), intent(in) :: otherMatrixInteger


		if ( allocated(  otherMatrixInteger%values ) ) then

			if ( allocated(  this%values ) ) then

				if ( ( size(this%values, DIM=1) /= size(otherMatrixInteger%values, DIM=1) ) .or. &
					(size(this%values, DIM=2) /= size(otherMatrixInteger%values, DIM=2) )  ) then

					deallocate( this%values )
					allocate( this%values( size(otherMatrixInteger%values, DIM=1), size(otherMatrixInteger%values, DIM=2) ) )

				end if

			else

				allocate( this%values( size(otherMatrixInteger%values, DIM=1), size(otherMatrixInteger%values, DIM=2) ) )

			end if

			this%values = otherMatrixInteger%values

		else

			call MatrixInteger_exception(WARNING, "The original matrix wasn't allocated ", "Class object MatrixInteger in the copyConstructor() function")

		end if



	end subroutine MatrixInteger_copyConstructor

	! !>
	! !! @brief Constructor
	! !!    Reserva la memoria necesaria para this y le asigna los valores del vector
	! !!    diagonalVector a sus elementos de la diagonal y al resto cero
	! !<
	! subroutine MatrixInteger_diagonalConstructor( this, diagonalVector )
	! 	implicit none
	! 	type(MatrixInteger), intent(inout) :: this
	! 	type(Vector), intent(in) :: diagonalVector

	! 	integer :: i

	! 	if( allocated( this%values ) ) deallocate( this%values )

	! 	if( allocated( diagonalVector%values ) ) then

	! 		allocate( this%values( size(diagonalVector%values), size(diagonalVector%values) ) )

	! 		call MatrixInteger_setIdentity( this )

	! 		do i=1, size( this%values, DIM=1 )
	! 			this%values(i, i) = diagonalVector%values(i)
	! 		end do

	! 	else

	! 		call MatrixInteger_exception( WARNING, "The original matrix wasn't allocated ", "Class object MatrixInteger in the copyConstructor() function" )

	! 	end if

	! end subroutine MatrixInteger_diagonalConstructor

	! !>
	! !! @brief Constructor
	! !! Reserva la memoria necesaria para this y le asigna los valores aleatorios
	! !! a sus elementos
	! !!
	! !! @param symmetric Si su valor es .true. genera una matriz simetrica, de lo
	! !!                 contrario genera una matriz no simetrica
	! !<
	! subroutine MatrixInteger_randomElementsConstructor( this, rows, cols, symmetric )
	! 	implicit none
	! 	type(MatrixInteger), intent(inout) :: this
	! 	integer, intent(in) :: rows
	! 	integer, intent(in) :: cols
	! 	logical, intent(in), optional :: symmetric

	! 	integer :: value
	! 	integer :: j
	! 	integer :: i                  ! Counts random numbers
	! 	integer(4) :: timeArray(3)    ! Holds the hour, minute, and second
	! 	logical :: symmetricTmp

	! 	if( allocated( this%values ) ) deallocate( this%values )

	! 	symmetricTmp = .false.
	! 	if( present(symmetric) ) symmetricTmp = symmetric

	! 	call itime(timeArray)     ! Get the current time
	! 	i = rand ( timeArray(1)+timeArray(2)+timeArray(3) )

	! 	allocate( this%values( rows, cols ) )

	! 	call MatrixInteger_setIdentity( this )

	! 	if( symmetricTmp ) then

	! 		do i=1, rows
	! 			do j=i, cols
	! 				this%values(i, j) = irand(0)
	! 				this%values(j, i) = this%values(i, j)
	! 			end do
	! 		end do

	! 	else

	! 		do i=1, rows
	! 			do j=1, cols
	! 				this%values(i, j) = irand(0)
	! 			end do
	! 		end do

	! 	end if

	! end subroutine MatrixInteger_randomElementsConstructor


	!>
	!! @brief  Destructor
	!<
	subroutine MatrixInteger_destructor( this )
		implicit none
		type(MatrixInteger), intent(inout) :: this

		if( allocated(  this%values ) ) deallocate( this%values )

	end subroutine MatrixInteger_destructor

	!>
	!! @brief  Imprime a salida estandar la matriz realizando cambio de linea
	!!              con un maximo de "Parameters%FORMAT_NUMBER_OF_COLUMNS" columnas
	!<
! 	subroutine MatrixInteger_show( this, rowKeys, columnKeys, flags )
! 		implicit none
! 		type(MatrixInteger), intent(in) :: this
! 		character(*), intent(in), optional :: rowKeys(:)
! 		character(*), intent(in), optional :: columnKeys(:)
! 		integer, intent(in), optional :: flags

! 		character(20) :: colNum
! 		integer :: auxColNum
! 		integer :: columns
! 		integer :: rows
! 		integer :: i
! 		integer :: j
! 		integer :: k
! 		integer :: lowerLimit
! 		integer :: upperLimit
! 		integer :: tmpFlags

! 		tmpFlags = WITHOUT_KEYS
! 		if( present(flags) ) then
! 			tmpFlags = flags
! 		end if

! 		rows = size( this%values, DIM=1 )
! 		columns = size( this%values, DIM=2 )

! 		if( present( rowKeys ) ) then
! 			if( size( rowKeys ) < rows ) then

! 				call MatrixInteger_exception(WARNING, "The size of row keys is low than number of matrix rows", &
! 					 "Class object MatrixInteger in the show() function" )

! 			end if
! 		end if

! 		if( present( columnKeys ) ) then
! 			if( size( columnKeys ) < columns ) then

! 				call MatrixInteger_exception(WARNING, "The size of column keys is low than number of matrix columns", &
! 					"Class object MatrixInteger in the show() function" )

! 			end if
! 		end if

! 		do k=1, ceiling( (columns * 1.0)/(Parameters%FORMAT_NUMBER_OF_COLUMNS * 1.0 ) )

! 			lowerLimit = Parameters%FORMAT_NUMBER_OF_COLUMNS * ( k - 1 ) + 1
! 			upperLimit = Parameters%FORMAT_NUMBER_OF_COLUMNS * ( k )
! 			auxColNum = Parameters%FORMAT_NUMBER_OF_COLUMNS

! 			if ( upperLimit > columns ) then
! 				auxColNum =  Parameters%FORMAT_NUMBER_OF_COLUMNS -  upperLimit + columns
! 				upperLimit = columns
! 			end if

! 			write(colNum,*) auxColNum

! 			if( present( columnKeys ) ) then

! 				if( tmpFlags == WITH_COLUMN_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
! 					write (6,"(21X,"//trim(colNum)//"A15)") ( columnKeys(i), i = lowerLimit, upperLimit )
! 				end if

! 			else

! 				if( tmpFlags /= WITHOUT_KEYS ) then
! 					if( tmpFlags == WITH_COLUMN_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
! 						write (6,"(5X,"//trim(colNum)//"I15)") ( i,i=lowerLimit,upperLimit )
! 					end if
! 				end if

! 			end if

! 			print *,""

! 			if( present( rowKeys ) ) then

! 				if( tmpFlags == WITH_ROW_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
! 					write (6,"(A18,"//trim(colNum)//"I15)") ( rowKeys(i), ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
! 				else
! 					write (6,"(5X,"//trim(colNum)//"I15)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
! 				end if

! 			else
! 				if( tmpFlags /= WITHOUT_KEYS ) then

! 					if( ( tmpFlags == WITH_ROW_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) .and. tmpFlags /= WITHOUT_KEYS ) then
! 						write (6,"(I5,"//trim(colNum)//"I15)") ( i, ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
! 					else
! 						write (6,"(5X,"//trim(colNum)//"I15)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
! 					end if

! 				else

! 					write (6,"(5X,"//trim(colNum)//"I15)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )

! 				end if
! 			end if

! 			print *,""

! 		end do

! 	end subroutine MatrixInteger_show

! 	!>
! 	!! @brief Devuelve un apuntador a la matrix solicitada
! 	!!
! 	!! @param this matrix de m x n
! 	!! @return Apuntador a la matriz solicitada.
! 	!! @todo No ha sido probada
! 	!<
! 	function MatrixInteger_getPtr( this ) result( output )
! 		implicit none
! 		type(MatrixInteger) , target , intent(in) :: this
! 		integer , pointer :: output(:,:)

! 		output => null()
! 		output => this%values

! 	end function MatrixInteger_getPtr

! 	!>
! 	!! @brief  Intercambia las filas i y j
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	subroutine MatrixInteger_swapRows( this, i, j )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		integer, intent(in) :: j

! 	end subroutine MatrixInteger_swapRows

! 	!>
! 	!! @brief  Intercambia las columnsa i y j
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	subroutine MatrixInteger_swapColumns( this, i, j )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		integer, intent(in) :: j

! 	end subroutine MatrixInteger_swapColumns

! 	!>
! 	!! @brief  Retorna el n�mero de filas de la matriz
! 	!<
! 	function MatrixInteger_getNumberOfRows( this ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer :: output

! 		output = size( this%values , DIM=1 )

! 	end function MatrixInteger_getNumberOfRows

! 	!>
! 	!! @brief  Retorna el n�mero de columnas de la matriz
! 	!<
! 	function MatrixInteger_getNumberOfColumns( this ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer :: output

! 		output = size( this%values , DIM=2 )

! 	end function MatrixInteger_getNumberOfColumns

! 	!>
! 	!! @brief Retorna el elemento de la columna i-esima y la fila j-esima
! 	!<
! 	function MatrixInteger_getElement( this, i, j ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		integer, intent(in) :: j
! 		integer :: output

! 		output = this%values( i, j )

! 	end function MatrixInteger_getElement

! 	!>
! 	!! @brief Selecciona el valor del elemento de la columna i-esima y la fila j-esima
! 	!<
! 	subroutine MatrixInteger_setElement( this, i, j, value )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		integer, intent(in) :: j
! 		integer, intent(in) :: value

! 		this%values( i, j ) = value

! 	end subroutine MatrixInteger_setElement

! 	!>
! 	!! @brief Retorna un vector con los elementos de la columna i-esima
! 	!!
! 	! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getColumn( this, i ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		type(Vector) , pointer :: output(:)

! 		output => null()
! !! 		output => this%values(:,n)

! 		!! A=M(:,n)//columnas

! 	end function MatrixInteger_getColumn

! 	!>
! 	!! @brief  Retorna un vector con los elementos de la fila i-esima
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getRow( this, i ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in) :: i
! 		type(Vector) , pointer :: output(:)

! 		output => null()
! !! 		output => this%values(x,:)

! 		!! A=M(x,:)//filas

! 	end function MatrixInteger_getRow

! 	!>
! 	!! @brief  Convierte la matriz en una matriz identidad
! 	!<
! 	subroutine MatrixInteger_setIdentity( this )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this

! 		integer :: i
! 		integer :: rows

! 		this%values = 0.0_8
! 		rows = size( this%values , DIM=1 )

! 		do i=1, rows

! 			this%values(i,i) = 1

! 		end do

! 	end subroutine MatrixInteger_setIdentity

! 	!>
! 	!! @brief Selecciona todos los valores de la matriz a cero
! 	!<
! 	subroutine MatrixInteger_setNull( this )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this

! 		this%values = 0

! 	end subroutine MatrixInteger_setNull

! 	!>
! 	!! Retorna el maximo valor encontrado dentro de la matriz
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getMax( this, colPos, rowPos ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(inout), optional :: colPos
! 		integer, intent(inout), optional :: rowPos
! 		integer :: output

! 		colPos = 0
! 		rowPos = 0
! 		output = 0

! 	end function MatrixInteger_getMax

! 	!>
! 	!! @brief  Retorna el minimo valor encontrado dentro de la matriz
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getMin( this, colPos, rowPos ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(inout), optional :: colPos
! 		integer, intent(inout), optional :: rowPos
! 		integer :: output

! 		colPos = 0
! 		rowPos = 0
! 		output = 0

! 	end function MatrixInteger_getMin

! 	!>
! 	!! @brief Retorna el determinante de la matriz
! 	!!
! 	!! @param flags Indica las propiedades adicionales de la matriz que
! 	!!              permite optimizar el calculo
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getDeterminant( this, method, flags ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer, intent(in), optional :: method
! 		integer, intent(in), optional :: flags
! 		integer :: output

! 		output = 0

! 	end function MatrixInteger_getDeterminant
! 	!>
! 	!! @brief  Retorna true si la matriz tiene todos sus elementos iguales
! 	!!               a cero, de lo contrario false
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_isNull( this ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		logical :: output

! 		output = .false.

! 	end function MatrixInteger_isNull

! 	!>
! 	!! @brief Retorna la transpuesta de la matriz
! 	!!
! 	!! @todo Falta implementar
! 	!<
! 	function MatrixInteger_getTransPose( this ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		integer :: output

! 		output = 0.0_8

! 	end function MatrixInteger_getTransPose

! 	!>
! 	!! @brief  Retorna la traza de la matriz
! 	!<
! 	function MatrixInteger_trace( this ) result ( output )
! 		implicit none
! 		type(MatrixInteger) , intent(in) :: this
! 		integer :: output

! 		integer :: i

! 		output = 0

! 		do i = 1, size( this%values, DIM=1 )
		
! 			output = output + this%values( i, i )

! 		end do

! 	end function MatrixInteger_trace
	
! 	!>
! 	!! @brief Suma dos matrices
! 	!!
! 	!! @todo Falta meter el soporte blas
! 	!<
! 	function MatrixInteger_plus( this, otherMatrixInteger ) result ( output )
! 		implicit none
! 		type(MatrixInteger), intent(inout) :: this
! 		type(MatrixInteger), intent(in) :: otherMatrixInteger
! 		type(MatrixInteger) :: output

! #ifdef BLAS_MATRIX_SUPPORT
! !! 		type(Vector) :: x
! !! 		integer :: alpha
! !! 		integer :: beta
! !! 		integer :: matrixOrder
! !!
! !! 		call Vector_constructor( x, size(this%values, DIM=1), value=1.0_8 )
! !! 		alpha = 1.0_8
! !! 		beta = 0.0_8
! !! 		matrixOrder = size( this%values, DIM=1 )
! !!
! !! 		call MatrixInteger_constructor( output, size(this%values, DIM=1), size(otherMatrixInteger%values, DIM=2) )
! !!
! !! 		call dsymv( &
! !! 			UPPER_TRIANGLE_IS_STORED, &
! !! 			matrixOrder, &
! !! 			alpha, &
! !! 			output%values, &
! !! 			matrixOrder, &
! !! 			x, &
! !! 			1, &
! !! 			beta, &
! !! 			x, &
! !! 			1 )
! #else
! 		call MatrixInteger_constructor( output, size(this%values, DIM=1), size(otherMatrixInteger%values, DIM=2) )
! 		output%values = this%values + otherMatrixInteger%values
! #endif

! 	end function MatrixInteger_plus


! 	!>
! 	!! @brief Remueve la fila especificada de una matriz
! 	!<
	subroutine MatrixInteger_removeRow( this, numberOfRow )
		implicit none
		type(MatrixInteger) :: this
		integer, intent(in) :: numberOfRow

		integer, allocatable :: auxArray(:,:)
		integer :: rows
		integer :: columns

		rows = size( this%values, dim=1 )

		if (numberOfRow <= rows ) then

			columns = size( this%values, dim=2 )

			allocate( auxArray(rows-1,columns) )
			auxArray(1:numberOfRow-1,:) = this%values(1:numberOfRow-1,:)
			auxArray(numberOfRow:rows-1,:) = this%values(numberOfRow+1:rows,:)
			deallocate( this%values )
			allocate( this%values(rows-1,columns) )
			this%values = auxArray
			deallocate( auxArray )

		end if

	end subroutine MatrixInteger_removeRow

! 	!>
! 	!! @brief Remueve la columna especificada de una matriz
! 	!<
! 	subroutine MatrixInteger_removeColumn( this, numberOfColumn )
! 		implicit none
! 		type(MatrixInteger) :: this
! 		integer, intent(in) :: numberOfColumn

! 		integer, allocatable :: auxArray(:,:)
! 		integer :: rows
! 		integer :: columns

! 		columns = size( this%values, dim=2 )

! 		if (numberOfColumn <= columns ) then

! 			rows = size( this%values, dim=1 )


! 			allocate( auxArray(rows,columns-1) )
! 			auxArray(:,1:numberOfColumn-1) = this%values(:,1:numberOfColumn-1)
! 			auxArray(:,numberOfColumn:columns-1) = this%values(:,numberOfColumn+1:columns)
! 			deallocate( this%values )
! 			allocate( this%values(rows,columns-1) )
! 			this%values = auxArray
! 			deallocate( auxArray )

! 		end if

! 	end subroutine MatrixInteger_removeColumn


	!>
	!! @brief  Maneja excepciones de la clase
	!<
	subroutine MatrixInteger_exception( typeMessage, description, debugDescription)
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

	end subroutine MatrixInteger_exception
	
end module MatrixInteger_
