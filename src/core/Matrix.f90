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
!! @brief Clase encargada de manipular todo lo relacionado con matrices
!!
!! Esta clase manipula todo lo relacionado con matrices de tipo numerico, ademas
!! de servir como una interface transparente para el uso de LAPACK en cuanto a metodos
!! que involuran algebra lineal
!!
!! @author Nestor Aguirre
!!
!! <b> Fecha de creacion : </b> 2008-08-19
!!   - <tt> 2007-08-19 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Define las estructurura del modulo
!!   - <tt> 2007-08-19 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2007-08-26 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Acoplamiento algunas funciones de LAPACK (dsyev)
!!   - <tt> 2007-09-17 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Adiciono metodos complementarios al modulo
!!   - <tt> 2009-04-22 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Establecio relacion de herencia con Exception class para llamado exciones
!!   - <tt> 2009-04-22 </tt>: Sergio Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Define y ajusta el formato de la documentacion para doxygen
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
!!   - <tt> 2014-02-12 </tt> Mauricio Rodas ( jmrodasr@unal.edu.co )
!!        -# Adapta el modulo para usar la libreria MAGMA-CUDA y acelerar en GPUs
!!   - <tt> 2015-02-22 </tt> Jorge Charry ( jacharrym@unal.edu.co )
!!        -# Adapt the module to the Lapack routine DSYEVX. 
!!   - <tt> 2015-03-17 </tt> Mauricio Rodas ( jmrodasr@unal.edu.co )
!!        -# Create the new routine for compute matrix-matrix operation using DGEMM of lapack
!!   - <tt> 2015-04-20 </tt> Jorge Charry ( jacharrym@unal.edu.co )
!!        -# Implement the swap rows routine. 
module Matrix_
  use CONTROL_
  use Exception_
  use Vector_
  use LapackInterface_
  !use CudaInterface_
  use Math_
  use String_
  use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none

  
  type, public :: Matrix
     real(8), allocatable :: values(:,:)
     real(8), allocatable :: eigenVectors(:,:)
     real(8), allocatable :: eigenValues(:)
     integer :: numberOfNegatives
     logical :: isPositiveDefinited
     logical :: isSymmetric
     logical :: isUnitary
     logical :: isInstanced
  end type Matrix

  type, public :: IMatrix8
     integer(8), allocatable :: values(:,:)
     logical :: isInstanced
  end type IMatrix8

  type, public :: IMatrix
     integer, allocatable :: values(:,:)
     logical :: isInstanced
  end type IMatrix

  type, public :: IMatrix1
     integer(1), allocatable :: values(:,:)
  end type IMatrix1
  
  type, public :: FourIndexMatrix
     real(8), allocatable :: values(:,:,:,:)
     logical :: isInstanced
  end type FourIndexMatrix

  interface assignment(=)
     module procedure Matrix_copyConstructor
  end interface
  
  
  !> enum Vector_printFormatFlags {
  integer, parameter, public :: WITH_COLUMN_KEYS = 1
  integer, parameter, public :: WITH_ROW_KEYS = 2
  integer, parameter, public :: WITH_BOTH_KEYS = 3
  integer, parameter, public :: WITHOUT_KEYS = 0
  integer, parameter, public :: WITHOUT_MESSAGES=4
  !> }
  
  !> enum Matrix_jobs {
  character, parameter, public :: COMPUTE_EIGENVALUES = 'N'
  character, parameter, public :: COMPUTE_EIGENVALUES_AND_EIGENVECTORS = 'V'
  !> }

  !> enum Matrix_storedFlags {
  character, parameter, public :: UPPER_TRIANGLE_IS_STORED = 'U'
  character, parameter, public :: LOWER_TRIANGLE_IS_STORED = 'L'
  !> }

  !> enum Matrix_type {
  !! @todo Hay que definir bien cules son los tipos de matrices a utilizar
  !! @todo Hay que intentar que sea posible mezclarlar, por ejemplo utilizando codigos binarios para su definicion
  integer, parameter, public :: SYMMETRIC		= 1
  integer, parameter, public :: DIAGONAL		= 2
  integer, parameter, public :: BIDIAGONAL		= 3
  integer, parameter, public :: TRIDIAGONAL		= 4
  integer, parameter, public :: TRIANGULAR		= 5
  integer, parameter, public :: UNKNOWN			= 6
  !> }

  public :: &
       Matrix_constructor, &
       Matrix_copyConstructor, &
       Matrix_diagonalConstructor, &
       Matrix_randomElementsConstructor, &
       Matrix_destructor, &
       Matrix_show, &
       Matrix_writeToFile, &
       Matrix_getFromFile, &
       Matrix_getPtr, &
       Matrix_swapRows, &
       Matrix_swapBlockOfColumns, &
       Matrix_swapColumns, &
       Matrix_getLinearlyIndependentVectors, &
       Matrix_selectLinearlyIndependentVectors, &
       Matrix_getNumberOfRows, &
       Matrix_getNumberOfColumns, &
       Matrix_getElement, &
       Matrix_setElement, &
       Matrix_getColumn, &
       Matrix_getRow, &
       Matrix_setIdentity, &
       Matrix_setNull, &
       Matrix_getMax, &
       Matrix_getMin, &
       Matrix_getDeterminant, &
       Matrix_inverse, &
       Matrix_orthogonalizeLastVector, &
       Matrix_isPositiveDefinited, &
       Matrix_orthogonalize, &
       Matrix_eigen, &
       Matrix_eigen_select, &
       Matrix_eigen_dsyevr, &
       Matrix_svd, &
       Matrix_symmetrize, &
       Matrix_solveLinearEquation, &
       Matrix_isNull, &
       Matrix_getTranspose, &
       Matrix_factorizeLU, &
       Matrix_trace, &
       Matrix_multiplication, &
       Matrix_product_dgemm, &
       Matrix_product, &
       Matrix_plus, &
       Matrix_pow, &
       Matrix_sqrt, &
       Matrix_log, &
       Matrix_linear, &
       Matrix_log10, &
       Matrix_sin, &
       Matrix_cos, &
       Matrix_tan, &
       Matrix_asin, &
       Matrix_acos, &
       Matrix_atan, &
       Matrix_sinh, &
       Matrix_cosh, &
       Matrix_tanh, &
       Matrix_boundEigenValues, &
       Matrix_standardDeviation, &
       Matrix_addRow, &
       Matrix_removeRow, &
       Matrix_addColumn, &
       Matrix_removeColumn, &
       Matrix_Fortran_orthogonalizeLastVector, &
       Matrix_eigenProperties, &
       diagonalize_matrix, & ! Copiada de Parakata
       Matrix_constructorInteger8, &
       Matrix_constructorInteger, &
       Matrix_constructorInteger1, &
       Matrix_fourIndexConstructor

  private

contains

  !>
  !! @brief Constructor
  !! Constructor por omision
  subroutine Matrix_constructor( this, dim1, dim2, value)
    implicit none
    type(Matrix), intent(inout) :: this
    integer(8), intent(in) :: dim1
    integer(8), intent(in) :: dim2
    real(8), optional, intent(in) :: value

    real(8) :: valueTmp
    this%isInstanced = .true.
    valueTmp = 0.0_8
    if( present(value) ) valueTmp = value
    if (allocated(this%values)) deallocate(this%values)
    allocate( this%values( dim1, dim2 ) )

    this%values = valueTmp
    this%isInstanced = .true.

  end subroutine Matrix_constructor
  
  !>
  !! @brief Constructor
  !! Constructor por omision
  subroutine Matrix_constructorInteger8( this, dim1, dim2, value)
    implicit none
    type(IMatrix8), intent(inout) :: this
    integer(8), intent(in) :: dim1
    integer(8), intent(in) :: dim2
    integer(8), optional, intent(in) :: value

    integer(8) :: valueTmp
    this%isInstanced = .true.
    valueTmp = 0.0_8
    if( present(value) ) valueTmp = value

    if (allocated(this%values)) deallocate(this%values)
    allocate( this%values( dim1, dim2 ) )

    this%values = valueTmp
    this%isInstanced = .true.

  end subroutine Matrix_constructorInteger8

  !>
  !! @brief Constructor
  !! Constructor por omision
  subroutine Matrix_constructorInteger( this, dim1, dim2, value)
    implicit none
    type(IMatrix), intent(inout) :: this
    integer(8), intent(in) :: dim1
    integer(8), intent(in) :: dim2
    integer, optional, intent(in) :: value

    integer :: valueTmp

    valueTmp = 0.0_8
    if( present(value) ) valueTmp = value

    if (allocated(this%values)) deallocate(this%values)
    allocate( this%values( dim1, dim2 ) )

    this%values = valueTmp
    this%isInstanced = .true.

  end subroutine Matrix_constructorInteger

  !>
  !! @brief Constructor
  !! Constructor por omision
  subroutine Matrix_constructorInteger1( this, dim1, dim2, value)
    implicit none
    type(IMatrix1), intent(inout) :: this
    integer(8), intent(in) :: dim1
    integer(8), intent(in) :: dim2
    integer(1), optional, intent(in) :: value
    integer(1) :: valueTmp

    valueTmp = 0
    if( present(value) ) valueTmp = value

    if (allocated(this%values)) deallocate(this%values)
    allocate( this%values( dim1, dim2 ) )

    this%values = valueTmp

  end subroutine Matrix_constructorInteger1

  !>
  !! @brief Constructor
  !! Constructor por omision
  subroutine Matrix_fourIndexConstructor( this, dim1, dim2, dim3, dim4, value)
    implicit none
    type(FourIndexMatrix), intent(inout) :: this
    integer(8), intent(in) :: dim1
    integer(8), intent(in) :: dim2
    integer(8), intent(in) :: dim3
    integer(8), intent(in) :: dim4
    real(8), optional, intent(in) :: value

    real(8) :: valueTmp
    this%isInstanced = .true.
    valueTmp = 0.0_8
    if( present(value) ) valueTmp = value
    if (allocated(this%values)) deallocate(this%values)
    allocate( this%values( dim1, dim2, dim3, dim4 ) )

    this%values = valueTmp
    this%isInstanced = .true.

  end subroutine Matrix_fourIndexConstructor

  !>
  !! @brief Constructor de copia
  !! Reserva la memoria necesaria para this y le asigna los valores de otherMatrix
  !! @warning Actualmenet solo copia los valores las matrices
  subroutine Matrix_copyConstructor( this, otherMatrix )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix), intent(in) :: otherMatrix


    if ( allocated(  otherMatrix%values ) ) then

       if ( allocated(  this%values ) ) then

          if ( ( size(this%values, DIM=1) /= size(otherMatrix%values, DIM=1) ) .or. &
               (size(this%values, DIM=2) /= size(otherMatrix%values, DIM=2) )  ) then

             deallocate( this%values )
             allocate( this%values( size(otherMatrix%values, DIM=1), size(otherMatrix%values, DIM=2) ) )

          end if

       else

          allocate( this%values( size(otherMatrix%values, DIM=1), size(otherMatrix%values, DIM=2) ) )

       end if

       this%values = otherMatrix%values

    else

       call Matrix_exception(WARNING, "The original matrix wasn't allocated ", "Class object Matrix in the copyConstructor() function")

    end if

  end subroutine Matrix_copyConstructor

  !>
  !! @brief Constructor
  !!    Reserva la memoria necesaria para this y le asigna los valores del vector
  !!    diagonalVector a sus elementos de la diagonal y al resto cero
  subroutine Matrix_diagonalConstructor( this, diagonalVector )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Vector), intent(in) :: diagonalVector

    integer :: i

    if( allocated( this%values ) ) deallocate( this%values )

    if( allocated( diagonalVector%values ) ) then

       allocate( this%values( size(diagonalVector%values), size(diagonalVector%values) ) )

       call Matrix_setIdentity( this )

       do i=1, size( this%values, DIM=1 )
          this%values(i, i) = diagonalVector%values(i)
       end do

    else

       call Matrix_exception( WARNING, "The original matrix wasn't allocated ", &
            "Class object Matrix in the copyConstructor() function" )

    end if

  end subroutine Matrix_diagonalConstructor

  !>
  !! @brief Constructor
  !! Reserva la memoria necesaria para this y le asigna los valores aleatorios
  !! a sus elementos
  !!
  !! @param symmetric Si su valor es .true. genera una matriz simetrica, de lo
  !!                 contrario genera una matriz no simetrica
  subroutine Matrix_randomElementsConstructor( this, rows, cols, symmetric )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    logical, intent(in), optional :: symmetric

    real(8) :: value
    integer :: j
    integer :: i                  ! Counts random numbers
    integer(4) :: timeArray(3)    ! Holds the hour, minute, and second
    logical :: symmetricTmp

    if( allocated( this%values ) ) deallocate( this%values )

    symmetricTmp = .false.
    if( present(symmetric) ) symmetricTmp = symmetric

    call itime(timeArray)     ! Get the current time
    !! i = rand ( timeArray(1)+timeArray(2)+timeArray(3) )

    allocate( this%values( rows, cols ) )

    call Matrix_setIdentity( this )

    if( symmetricTmp ) then

       do i=1, rows
          do j=i, cols
             !! this%values(i, j) = rand(0)
             this%values(j, i) = this%values(i, j)
          end do
       end do

    else

       do i=1, rows
          do j=1, cols
             !! this%values(i, j) = rand(0)
          end do
       end do

    end if

  end subroutine Matrix_randomElementsConstructor


  !>
  !! @brief  Destructor
  subroutine Matrix_destructor( this )
    implicit none
    type(Matrix), intent(inout) :: this
    this%isInstanced = .false.
    if ( allocated(  this%values ) ) deallocate( this%values )
    if ( allocated( this%eigenVectors) ) deallocate( this%eigenVectors )
    if ( allocated( this%eigenValues) ) deallocate( this%eigenValues )

  end subroutine Matrix_destructor

  !>
  !! @brief  Imprime a salida estandar la matriz realizando cambio de linea
  !!              con un maximo de "CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS" columnas
  subroutine Matrix_show( this, rowKeys, columnKeys, flags )
    implicit none
    type(Matrix), intent(in) :: this
    character(*), intent(in), optional :: rowKeys(:)
    character(*), intent(in), optional :: columnKeys(:)
    integer, intent(in), optional :: flags

    integer :: auxColNum
    integer :: columns
    integer :: rows
    integer :: i
    integer :: j
    integer :: k
    integer :: lowerLimit
    integer :: upperLimit
    integer :: tmpFlags
    character(20) :: colNum

    colNum = ""

    tmpFlags = WITHOUT_KEYS
    if( present(flags) ) then
       tmpFlags = flags
    end if

    rows = size( this%values, DIM=1 )
    columns = size( this%values, DIM=2 )

    if( present( rowKeys ) ) then
       if( size( rowKeys ) < rows ) then

          call Matrix_exception(WARNING, "The size of row keys is low than number of matrix rows", &
               "Class object Matrix in the show() function" )

       end if
    end if

    if( present( columnKeys ) ) then
       if( size( columnKeys ) < columns ) then

          call Matrix_exception(WARNING, "The size of column keys is low than number of matrix columns", &
               "Class object Matrix in the show() function" )

       end if
    end if

    do k=1, ceiling( (columns * 1.0)/(CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * 1.0 ) )

       lowerLimit = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k - 1 ) + 1
       upperLimit = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS * ( k )
       auxColNum = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS

       if ( upperLimit > columns ) then
          auxColNum =  CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS -  upperLimit + columns
          upperLimit = columns
       end if

       write(colNum,*) auxColNum

       if( present( columnKeys ) ) then

          if( tmpFlags == WITH_COLUMN_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
             write (6,"(21X,"//trim(colNum)//"A15)") ( columnKeys(i), i = lowerLimit, upperLimit )
          end if

       else

          if( tmpFlags /= WITHOUT_KEYS ) then
             if( tmpFlags == WITH_COLUMN_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
                write (6,"(5X,"//trim(colNum)//"I15)") ( i,i=lowerLimit,upperLimit )
             end if
          end if

       end if

       print *,""

       if( present( rowKeys ) ) then

          if( tmpFlags == WITH_ROW_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) then
             write (6,"(A18,"//trim(colNum)//"F15.6)") ( rowKeys(i), ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
          else
             write (6,"(5X,"//trim(colNum)//"F15.6)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
          end if

       else
          if( tmpFlags /= WITHOUT_KEYS ) then

             if( ( tmpFlags == WITH_ROW_KEYS .or. tmpFlags == WITH_BOTH_KEYS ) .and. tmpFlags /= WITHOUT_KEYS ) then
                write (6,"(I5,"//trim(colNum)//"F15.6)") ( i, ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
             else
                write (6,"(5X,"//trim(colNum)//"F15.6)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )
             end if

          else

             write (6,"(5X,"//trim(colNum)//"F15.6)") ( ( this%values(i,j), j=lowerLimit,upperLimit ), i = 1, rows )

          end if
       end if

       print *,""

    end do

  end subroutine Matrix_show

  !>
  !! @brief Escribe una matriz en el lugar especificado (en binario o con formato)
  !! @warning Arguments option are available only for unit option and for binary format.
  subroutine Matrix_writeToFile( mmatrix, unit, file, binary, arguments)    
    implicit none
    type(Matrix) :: mmatrix
    integer, optional :: unit
    character(*), optional :: file
    logical, optional :: binary
    character(*), optional :: arguments(:)

    integer :: rows
    integer :: columns
    integer :: status
    integer :: m
    integer :: n
    character(50) :: auxSize
    logical :: bbinary
    logical :: existFile

    bbinary = .false.
    if (present(binary)) bbinary = binary

    if (bbinary)  then
       
       if ( present( unit ) ) then
          !! It is assumed that the unit y conected to any file (anyways will check)
          inquire(unit=unit, exist=existFile)
          
          if(existFile) then
             
             if(present(arguments)) then
                
                do m = 1, size(arguments)
                   
                   write(unit) arguments(m)
                   
                end do
             end if
             
             write(unit) int(size(mmatrix%values), 8)
             write(unit) mmatrix%values

          else
             
             call Matrix_exception( ERROR, "Unit file no connected!",&
                  "Class object Matrix  in the writeToFile() function" )
             
          end if

             
       else if ( present(file) ) then

          open ( 4,FILE=trim(file),STATUS='REPLACE',ACTION='WRITE', FORM ='UNFORMATTED')
          write(4) int(size(mmatrix%values), 8)
          write(4) mmatrix%values
          close(4)
       end if

    else

       if ( present( unit ) ) then

          !! Falta implementar

          inquire(unit=unit, exist=existFile)


          if(existFile) then

             if(present(arguments)) then
                
                do m = 1, size(arguments)
                   
                   write(unit,*) arguments(m)
                   
                end do
             end if
            rows = size( mmatrix%values , DIM=1 )
            columns = size( mmatrix%values , DIM=2 )
            write(auxSize,*)  columns
            write(unit,*) int(size(mmatrix%values), 8)
            write (unit,"("//trim(auxSize)//"(1X,ES15.8))") ((mmatrix%values(m,n), n=1 ,columns) , m=1 , rows)
          end if

       else if ( present(file) ) then

          open ( 4,FILE=trim(file),STATUS='REPLACE',ACTION='WRITE')
          rows = size( mmatrix%values , DIM=1 )
          columns = size( mmatrix%values , DIM=2 )
          write(auxSize,*)  columns
          write (4,"("//trim(auxSize)//"ES15.8)") ((mmatrix%values(m,n), n=1 ,columns) , m=1 , rows)
          close(4)


       end if
    end if


  end subroutine Matrix_writeToFile

  !>
  !! @brief Obtiene una matriz del lugar especificado
  !! @warning The arguments options are only available with unit option
  function Matrix_getFromFile(rows, columns, unit, file, binary, arguments, failContinue) result( output )
    implicit none
    integer, intent(in) :: rows
    integer, intent(in) :: columns
    integer, optional :: unit
    character(*), optional :: file
    logical, optional :: binary
    character(*), optional :: arguments(:)
    logical, optional :: failContinue

    type(Matrix) :: output

    integer :: failAction
    character(5000) :: line
    character(20) :: auxSize
    real(8), allocatable :: values(:)
    integer(8) :: totalSize
    integer :: status
    integer :: m, i
    integer :: n, j
    logical :: bbinary
    logical :: existFile
    logical :: found
    
    write(auxSize,*) columns

    bbinary = .false.
    if(present(binary)) bbinary = binary

    if(present(failContinue)) then
       failAction=WARNING
    else
       failAction=ERROR
    end if

    if(bbinary) then
       
       !! it is assumed that if you want to load a file for unit parameter, that file must be connected to unit "unit".
       if ( present( unit ) ) then
          
          !! check file
          inquire(unit=unit, exist=existFile)
          
          if(existFile) then
             
             rewind(unit)
             
             found = .false.
             
             if(present(arguments)) then
                
                do          
                   line = ""
                   read(unit, iostat = status) line (1:len_trim(arguments(1)))
                
                   if(status == -1) then
                      call Matrix_exception( failAction, "End of file!",&
                           "Class object Matrix in the getfromFile() function "//trim(arguments(1))//" "//trim(arguments(2)) )
                      return
                   end if
                   
                   if(trim(line) == trim(arguments(1))) then
                      
                      found = .true.                   
                      
                   end if
                   
                   if(found) then
                      
                      backspace(unit)
                      
                      do i = 1, size(arguments)
                         
                         found = .false.
                         read(unit, iostat = status) line
                         
                         if(trim(line) == trim(arguments(i))) then
                            
                            found = .true.
                                                     
                         end if
                         
                      end do
                      
                   end if
                   
                   if(found) exit
                      
                end do
                
             end if
             
             !! check size
             read(unit) totalSize
             
             if(totalSize == int(columns*rows,8)) then

                if(.not. allocated(output%values)) then                   
                
                   call Matrix_constructor( output, int(rows,8), int(columns,8) )
                   
                end if
                
                if (allocated(values)) deallocate(values)
                allocate( values(totalSize) )
                
                read(unit) values
                
                m = 1

                ! do i=1,rows
                !    do j=1,columns
                      ! output%values(j, i) = values(m)
                      ! FELIX: We may need to check this again to read asymmetrical matrices
                do j=1,columns
                   do i=1,rows
                      output%values(i, j) = values(m)
                      m = m + 1
                   end do
                end do
                
                !! call Matrix_show(output)
                
                deallocate(values)

             else
                
                call Matrix_exception( failAction, "The dimensions of the matrix "//trim(arguments(1))//" "//trim(arguments(2))//" are wrong ",&
                     "Class object Matrix  in the getFromFile() function"  )
                return
             end if
             
          else
             call Matrix_exception( failAction, "Unit file no connected!",&
                  "Class object Matrix  in the getFromFile() function" )
             return

          end if
          
       else if ( present(file) ) then
          
          inquire( FILE = trim(file), EXIST = existFile )
          
          if ( existFile ) then
             
             open( 4, FILE=trim(file), ACTION='read', FORM='unformatted', STATUS='old' )
             
             !!Comprueba el tamanho de las matrices
             read(4) totalSize
             
             if(totalSize == int(columns*rows,8)) then
                
                call Matrix_constructor( output, int(rows,8), int(columns,8) )
                
                if (allocated(values)) deallocate(values)
                allocate( values(totalSize) )
                
                read(4) values
                close(4)
                
                m = 1
                do i=1,rows
                   do j=1,columns
                      output%values(j, i) = values(m)
                      m = m + 1
                   end do
                end do
                
                deallocate(values)
                
             else
                
                close(4)
                
                call Matrix_exception( failAction, "The dimensions of the matrix "//trim(arguments(1))//" "//trim(arguments(2))//" are wrong ",&
                     "Class object Matrix  in the getFromFile() function"  )
                return

             end if
             
          else
             call Matrix_exception( failAction, "The file "//trim(file)//" don't exist ",&
                  "Class object Matrix  in the getFromFile() function" )
             return
                                   
          end if
       end if
       
    else
       
       if ( present( unit ) ) then
          
          !! check file
          inquire(unit=unit, exist=existFile)
          
          if(existFile) then
             
             rewind(unit)
             
             found = .false.
             
             if(present(arguments)) then
                
                do          
                   line = ""
                   read(unit,*, iostat = status) line !(1:len_trim(arguments(1)))
                
                   if(status == -1) then
                      
                      call Matrix_exception( failAction, "End of file!",&
                           "Class object Matrix in the getfromFile() function" )
                      return

                   end if
                   
                   if(trim(line) == trim(arguments(1))) then
                      
                      found = .true.                   
                      
                   end if
                   
                   if(found) then
                      
                      backspace(unit)
                      
                      do i = 1, size(arguments)
                         
                         found = .false.
                         read(unit,*, iostat = status) line
                         
                         if(trim(line) == trim(arguments(i))) then
                            
                            found = .true.
                                                     
                         end if
                         
                      end do
                      
                   end if
                   
                   if(found) exit
                      
                end do
                
             end if
             
             !! check size
             read(unit,*) totalSize
             
             if(totalSize == int(columns*rows,8)) then

                if(.not. allocated(output%values)) then                   
                
                   call Matrix_constructor( output, int(rows,8), int(columns,8) )
                   
                end if
                
                if (allocated(values)) deallocate(values)
                allocate( values(totalSize) )
                
                read(unit,*) values
                
                m = 1
                
                do i=1,rows
                   do j=1,columns
!                      output%values(j, i) = values(m) !! ???
                      output%values(i, j) = values(m)
                      m = m + 1
                   end do
                end do
                
                !! call Matrix_show(output)
                
                deallocate(values)

            else
                
                call Matrix_exception( failAction, "The dimensions of the matrix "//trim(arguments(1))//" "//trim(arguments(2))//" are wrong ",&
                     "Class object Matrix  in the getFromFile() function"  )
                return

                
             end if
             
          else
             call Matrix_exception( failAction, "Unit file no connected!",&
                  "Class object Matrix  in the getFromFile() function" )
             return

             
          end if

          
       else if ( present(file) ) then
          
          inquire( FILE = trim(file), EXIST = existFile )
          
          if ( existFile ) then
             
             call Matrix_constructor( output, int(rows,8), int(columns,8) )
             
             open ( 4,FILE=trim(file),STATUS='unknown',ACCESS='sequential' )
             
             !! verifica el tamahno de las matrices
             read  (4,"(A)") line
             
             if ( len_trim(line) /=columns*15) then
                
                close(4)
                
                call Matrix_exception( failAction, "The dimensions of the matrix "//trim(arguments(1))//" "//trim(arguments(2))//" are wrong ",&
                     "Class object Matrix  in the getFromFile() function"  )
                return

                
             else
                
                rewind(4)
                
             end if
             
             read ( 4,"("//trim(auxSize)//"ES15.8)") ( (output%values(m,n), n=1 ,columns) , m=1, rows )
             close(4)
             
          else
             
             call Matrix_exception( failAction, "The file "//trim(file)//" don't exist ",&
                  "Class object Matrix  in the getFromFile() function" )
             return

             
          end if
          
       end if
       
    end if

  end function Matrix_getFromFile

  !>
  !! @brief Devuelve un apuntador a la matrix solicitada
  !! @param this matrix de m x n
  !! @return Apuntador a la matriz solicitada.
  !! @todo No ha sido probada
  function Matrix_getPtr( this ) result( output )
    implicit none
    type(Matrix) , target , intent(in) :: this
    real(8) , pointer :: output(:,:)

    output => null()
    output => this%values

  end function Matrix_getPtr

  !>
  !! @brief  Intercambia las filas i y j
  !! @todo Falta implementar
  subroutine Matrix_swapRows( this, i, j )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(8), allocatable :: value1(:)
    real(8), allocatable :: value2(:)

    if (allocated(value1)) deallocate(value1)
    allocate (value1(size(this%values,dim=2)) )
    if (allocated(value2)) deallocate(value2)
    allocate (value2(size(this%values,dim=2)) )


    value1 = this%values(i,:) 
    value2 = this%values(j,:) 

    this%values(i,:) = value2
    this%values(j,:) = value1

  end subroutine Matrix_swapRows

  !>
  !! @brief  Intercambia dos bloques de columnas especificados por los rangos A y B
  !! @warning Coloca el bloque de columnas especificado al inicio de la matriz
  !!		el resto de columnas al final de la misma
  !! @warning Actualmente no se soportan rangos intermedios, solamente rangos contains
  !!			abiertos que incluyan el elemento terminal
  !! @todo Dar soporte a rangos no consecutivos
  subroutine Matrix_swapBlockOfColumns( this, rangeSpecification )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: rangeSpecification(2)

    real(8), allocatable :: auxMatrix(:,:)

    allocate( auxMatrix(size(this%values,dim=1), size(this%values,dim=2) ) )
    auxMatrix=this%values

    this%values(:, 1: rangeSpecification(2)-rangeSpecification(1)+1) = auxMatrix(:, rangeSpecification(1):rangeSpecification(2) )
    this%values(:, rangeSpecification(2)-rangeSpecification(1)+2:size(this%values,dim=2) ) = auxMatrix(:,1:rangeSpecification(1)-1)

    deallocate(auxMatrix)

  end subroutine Matrix_swapBlockOfColumns


  !>
  !! @brief  Intercambia las columnsa i y j
  !! @todo Falta implementar
  subroutine Matrix_swapColumns( this, i, j )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j


  end subroutine Matrix_swapColumns

  !>
  !! @brief devuelve una matrix ortogonal con columnas correspondientes
  !! 	a la columnas linealmentes independientes de la matriz pasada como
  !! 	parametro.
  function Matrix_getLinearlyIndependentVectors(this) result(output)
    implicit none
    type(Matrix) :: this
    type(Matrix) :: output

    type(Matrix) :: auxMatrix
    type(Matrix) :: linearlyIndependentVectors
    type(Vector) :: auxVector
    real(8) :: auxVal
    integer :: ssize
    integer ::numOfLinearlyDependent
    integer :: numOfLinearlyIndependent
    integer :: i
    integer :: j

    auxMatrix=this
    ssize = size(this%values,dim=1)
    call Vector_constructor(auxVector, ssize, 0.0_8)
    linearlyIndependentVectors=this
    linearlyIndependentVectors%values=0.0

    ssize = size(this%values,dim=2)
    i=1
    numOfLinearlyDependent = 0
    numOfLinearlyIndependent=0

    do while( i <= ( ssize- numOfLinearlyDependent ) )
       auxVector%values=auxMatrix%values(:,i)
       do j=1,numOfLinearlyIndependent
          auxVal=-dot_product(auxMatrix%values(:,i),linearlyIndependentVectors%values(:,j))
          auxVector%values=auxVal*linearlyIndependentVectors%values(:,j)+auxVector%values
       end do
       auxVal=sqrt( dot_product(auxVector%values, auxVector%values) )
       if( abs(auxVal) < 1.0D-4) then
          auxVector%values=auxMatrix%values(:,i)
          auxMatrix%values(:,i)=auxMatrix%values(:,ssize-numOfLinearlyDependent)
          auxMatrix%values(:,ssize-numOfLinearlyDependent) = auxVector%values
          numOfLinearlyDependent=numOfLinearlyDependent+1
       else
          auxVector%values = (1.0/auxVal) * auxVector%values
          numOfLinearlyIndependent= numOfLinearlyIndependent+1
          linearlyIndependentVectors%values(:,numOfLinearlyIndependent)=auxVector%values
          i=i+1
       end if
    end do

    call Matrix_constructor(output,int(numOfLinearlyIndependent,8),int(numOfLinearlyIndependent,8))
    output%values = linearlyIndependentVectors%values(:,1:numOfLinearlyIndependent)
    call Matrix_destructor(linearlyIndependentVectors)
    call Vector_destructor(auxVector)
    call Matrix_destructor(auxMatrix)


  end function Matrix_getLinearlyIndependentVectors


  !>
  !! @brief devuelve una matrix en la cual se reordenan  las columnas
  !! 	linealmente independientes (ortogonalizadas) al principio de la misma.
  subroutine Matrix_selectLinearlyIndependentVectors(this, numOfLinearlyIndependent)
    implicit none
    type(Matrix) :: this
    integer, intent(inout) :: numOfLinearlyIndependent

    type(Matrix) :: linearlyIndependentVectors
    type(Vector) :: auxVector
    real(8) :: auxVal
    integer :: ssize
    integer ::numOfLinearlyDependent
    integer :: i
    integer :: j

    ssize = size(this%values,dim=1)
    call Vector_constructor(auxVector, ssize, 0.0_8)
    linearlyIndependentVectors=this
    linearlyIndependentVectors%values=0.0

    ssize = size(this%values,dim=2)
    i=1

    numOfLinearlyDependent = 0
    numOfLinearlyIndependent=0

    do while( i <= ( ssize- numOfLinearlyDependent ) )
       auxVector%values=this%values(:,i)
       do j=1,numOfLinearlyIndependent
          auxVal=-dot_product(this%values(:,i),linearlyIndependentVectors%values(:,j))
          auxVector%values=auxVal*linearlyIndependentVectors%values(:,j)+auxVector%values
       end do
       auxVal=sqrt( dot_product(auxVector%values, auxVector%values) )
       if( abs(auxVal) < 1.0D-4) then
          auxVector%values=this%values(:,i)
          this%values(:,i)=this%values(:,ssize-numOfLinearlyDependent)
          this%values(:,ssize-numOfLinearlyDependent) = auxVector%values
          numOfLinearlyDependent=numOfLinearlyDependent+1
       else
          auxVector%values = (1.0/auxVal) * auxVector%values
          numOfLinearlyIndependent= numOfLinearlyIndependent+1
          linearlyIndependentVectors%values(:,numOfLinearlyIndependent)=auxVector%values
          i=i+1
       end if
    end do

    !! Ortogonaliza los vectores linealmente independientes

    do i=1, numOfLinearlyIndependent
       do j=1,i-1
          auxVal=-dot_product( this%values(:,i), this%values(:,j) )
          this%values(:,i)=auxVal*this%values(:,j)+this%values(:,i)
       end do
       auxVal=1.0/sqrt( dot_product( this%values(:,i), this%values(:,i) ))
       this%values(:,i) = this%values(:,i)*auxVal
    end do

  end subroutine Matrix_selectLinearlyIndependentVectors


  !>
  !! @brief  Retorna el numero de filas de la matriz
  function Matrix_getNumberOfRows( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer :: output

    output = size( this%values , DIM=1 )

  end function Matrix_getNumberOfRows

  !>
  !! @brief  Retorna el numero de columnas de la matriz
  function Matrix_getNumberOfColumns( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer :: output

    output = size( this%values , DIM=2 )

  end function Matrix_getNumberOfColumns

  !>
  !! @brief Retorna el elemento de la columna i-esima y la fila j-esima
  function Matrix_getElement( this, i, j ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(8) :: output

    output = this%values( i, j )

  end function Matrix_getElement

  !>
  !! @brief Selecciona el valor del elemento de la columna i-esima y la fila j-esima
  subroutine Matrix_setElement( this, i, j, value )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(8), intent(in) :: value

    this%values( i, j ) = value

  end subroutine Matrix_setElement

  !>
  !! @brief Retorna un vector con los elementos de la columna i-esima
  !! @todo Falta implementar
  function Matrix_getColumn( this, i ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    type(Vector) , pointer :: output(:)

    output => null()
    !! 		output => this%values(:,n)

    !! A=M(:,n)//columnas

  end function Matrix_getColumn

  !>
  !! @brief  Retorna un vector con los elementos de la fila i-esima
  !! @todo Falta implementar
  function Matrix_getRow( this, i ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in) :: i
    type(Vector) , pointer :: output(:)

    output => null()
    !! 		output => this%values(x,:)

    !! A=M(x,:)//filas

  end function Matrix_getRow

  !>
  !! @brief  Convierte la matriz en una matriz identidad
  subroutine Matrix_setIdentity( this )
    implicit none
    type(Matrix), intent(inout) :: this

    integer :: i
    integer :: rows

    this%values = 0.0_8
    rows = size( this%values , DIM=1 )

    do i=1, rows

       this%values(i,i) = 1.0_8

    end do

  end subroutine Matrix_setIdentity

  !>
  !! @brief Selecciona todos los valores de la matriz a cero
  subroutine Matrix_setNull( this )
    implicit none
    type(Matrix), intent(inout) :: this

    this%values = 0.0_8

  end subroutine Matrix_setNull

  !>
  !! Retorna el maximo valor encontrado dentro de la matriz
  !! @todo Falta implementar
  function Matrix_getMax( this, colPos, rowPos ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(inout), optional :: colPos
    integer, intent(inout), optional :: rowPos
    real(8) :: output

    colPos = 0
    rowPos = 0
    output = 0.0_8

  end function Matrix_getMax

  !>
  !! @brief  Retorna el minimo valor encontrado dentro de la matriz
  !! @todo Falta implementar
  function Matrix_getMin( this, colPos, rowPos ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(inout), optional :: colPos
    integer, intent(inout), optional :: rowPos
    real(8) :: output

    colPos = 0
    rowPos = 0
    output = 0.0_8

  end function Matrix_getMin

  !>
  !! @brief Retorna el determinante de la matriz
  !! @param method se calcula a partir de una descomposicion SVD (valor absoluto) o LU
  subroutine Matrix_getDeterminant( this, determinant, method)
    implicit none
    type(Matrix), intent(inout) :: this
    real(8) :: determinant
    character(*), optional :: method

    type(Matrix) :: range, nullSpace, singular
    type(Matrix) :: U, auxMatrix
    integer :: dim, i
    integer, allocatable :: pivotIndices(:)
    character(10) :: selectedMethod

    if( present ( method ) ) then
       selectedMethod=trim(method)
    else
       selectedMethod="LU"
    end if
    dim=size(this%values, dim=1)
    
    select case( trim(selectedMethod) )
    case("SVD")

       call Matrix_constructor(range, int(dim,8), int(dim,8), 0.0_8)
       call Matrix_constructor(nullSpace, int(dim,8), int(dim,8), 0.0_8)
       call Matrix_constructor(singular, int(dim,8), int(dim,8), 0.0_8)

       call Matrix_svd( this, range, nullSpace, singular )

       determinant=1.0
       do i=1,dim
          determinant=determinant*singular%values(i,i)
       end do

    case("LU")

       call Matrix_constructor(U, int(dim,8), int(dim,8), 0.0_8)
       allocate( pivotIndices( dim ))
       
       auxMatrix=Matrix_factorizeLU( this, pivotIndices=pivotIndices, U=U )
       
       determinant=1.0
       do i=1,dim
          determinant=determinant*U%values(i,i)
          if(pivotIndices(i) .ne. i) determinant=-determinant
       end do

    case default
       call Matrix_exception(ERROR, "The selected method to compute the determinant is not implemented", "Class object Matrix in the getDeterminant() function")
    end select

  end subroutine Matrix_getDeterminant

  !>
  !! @brief Retorna la matriz inversa de la matriz
  !! @param flags Indica las propiedades adicionales de la matriz que
  !!              permite optimizar el calculo
  !! @todo Falta implementar
  function Matrix_inverse( this, method, flags,printFormatFlags, info ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    integer, intent(in), optional :: method
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: printFormatFlags
    integer, optional,intent(out) :: info
    type(Matrix) :: output

    integer :: matrixSize
    integer :: infoProcess
    real(8), allocatable :: workSpace(:)
    integer, allocatable :: pivotIndices(:)
    logical :: printError

    printError=.true.
    if(present(printFormatFlags)) then
       if(printFormatFlags/=WITHOUT_MESSAGES) printError=.false.
    end if
    
    !! Detemina variables y parametros requeridos para el calculo
    matrixSize = size( this%values, DIM=1 )

    allocate( pivotIndices( matrixSize ))
    allocate( workSpace( matrixSize ) )

    call Matrix_copyConstructor( output, this )

    if (present(printFormatFlags)) then
       output = Matrix_factorizeLU( output, pivotIndices=pivotIndices, printFormatFlags=printFormatFlags )
    else
       output = Matrix_factorizeLU( output, pivotIndices=pivotIndices)
    end if

    !! Invierte la matriz de entrada
    call dgetri( &
         matrixSize, &
         output%values, &
         matrixSize, &
         pivotIndices, &
         workSpace, &
         matrixSize, &
         infoProcess )

    !! Determina la ocurrencia de errores
    if ( infoProcess /= 0 .and. printError )  then
       call Matrix_exception(WARNING, "Factorization failed", "Class object Matrix in the factorizeLU() function" )
    end if

    if (present(info)) info=infoProcess

    !! libera memoria
    deallocate(workSpace)
    deallocate(pivotIndices)


  end function Matrix_inverse

  !>
  !! @brief Ortogonaliza el ultimo vector asumiendo que los anteriores ya estan ortogonalizados
  !! @todo Implementar el metodo modificado de Gram-Schmidt
  function Matrix_orthogonalizeLastVector( this ) result(output)
    implicit none
    type(Matrix) :: this
    type(Matrix) :: output

    integer :: i
    integer :: last
    real(8) :: squareNorm
    real(8) :: projectionOverOrthogonalizedBasis
    real(8) :: auxValue

    last = size(this%values,dim=2)
    allocate( output%values( size(this%values,dim=1), last ) )
    output%values = this%values

    !!***********************************************************************************
    !! Realiza de ortogonalizacion sobre los last-1 vectores, previamente ortogonalizados.
    !!
    do i=1,last-1
       squareNorm = dot_product( output%values(:,i), output%values(:,i) )
       projectionOverOrthogonalizedBasis=dot_product( output%values(:,i),output%values(:,last) )
       if ( squareNorm>CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
          output%values( :, last ) = output%values( :, last ) - projectionOverOrthogonalizedBasis/squareNorm*output%values(:,i)
       end if
    end do
    !!
    !!***********************************************************************************

  end function Matrix_orthogonalizeLastVector

  !>
  !! @brief Ortogonaliza el ultimo vector asumiendo que los anteriores ya estan ortogonalizados
  !! @todo Implementar el metodo modificado de Gram-Schmidt
  function Matrix_Fortran_orthogonalizeLastVector( thisFortran ) result(output)
    implicit none
    real(8) :: thisFortran(:,:)
    real(8), allocatable :: output(:,:)

    integer :: i
    integer :: last
    real(8) :: squareNorm
    real(8) :: projectionOverOrthogonalizedBasis
    real(8) :: auxValue

    last = size(thisFortran,dim=2)
    allocate( output( size(thisFortran,dim=1), last ) )
    output = thisFortran

    !!***********************************************************************************
    !! Realiza de ortogonalizacion sobre los last-1 vectores, previamente ortogonalizados.
    !!
    do i=1,last-1
       squareNorm = dot_product( output(:,i), output(:,i) )
       projectionOverOrthogonalizedBasis=dot_product( output(:,i),output(:,last) )
       if ( squareNorm>CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
          output( :, last ) = output( :, last ) - projectionOverOrthogonalizedBasis/squareNorm*output(:,i)
       end if
    end do
    !!
    !!***********************************************************************************

  end function Matrix_Fortran_orthogonalizeLastVector


  !>
  !! @brief Ortogonaliza mediante proceso Gram-Schmidt, las columnas de la matriz pasada como parametro
  function Matrix_orthogonalize( this ) result(output)
    implicit none
    type(Matrix) :: this
    type(Matrix) :: output

    integer :: i
    integer :: last
    real(8) :: squareNorm
    real(8) :: projectionOverOrthogonalizedBasis
    real(8) :: auxValue


    last = size(this%values,dim=2)
    allocate( output%values( size(this%values,dim=1), last ) )
    output%values = 0.0_8
    output%values = this%values

    !!
    !! Realiza de ortogonalizacion consecutiva de cada uno de los vectores
    !! presentes en la matriz
    !!
    do i=2,last

       output%values(:,1:i)=Matrix_Fortran_orthogonalizeLastVector( output%values(:,1:i) )

    end do

  end function Matrix_orthogonalize


  !>
  !! @brief  Calcula los values propios y opcionalmente los vectores propios de una matriz
  !! @param De manera opcional se pueden almacenar los vectores propios en la
  !!        variable eigenVectors
  !! @return Retorna los valores propios
  !! @todo Hay que agregar las rutinas de Lapack para los diferentes tipos de matrices
  !!       como NONSYMMETRIC, entre otras

! #ifdef CUDA

!   subroutine Matrix_eigen( this, eigenValues, eigenVectors, flags, m, dm )
!     implicit none
!     type(Matrix), intent(in) :: this
!     type(Vector), intent(inout) :: eigenValues
!     type(Matrix), intent(inout), optional :: eigenVectors
!     integer, intent(in), optional :: flags
!     integer, intent(in), optional :: dm
!     real(8), intent(in), optional :: m(:,:)
!     real(8), allocatable :: vectorAux(:)

!     integer :: lengthWorkSpace
!     integer :: matrixSize
!     integer :: infoProcess
!     real(8), allocatable :: workSpace(:)
!     type(Matrix) :: eigenVectorsTmp
!     integer :: i
!     integer :: matrixOrder


!     COMMON/NCB/ matrixSize, matrixOrder, lengthWorkSpace, infoProcess

!     matrixSize = size( this%values, DIM=1 )
!     matrixOrder = matrixSize

!     allocate(vectorAux (matrixSize))

!     if( flags == SYMMETRIC ) then

!        !! Determina la longitud adecuada del vector de trabajo
!        lengthWorkSpace=3*matrixSize-1

!        !! Crea el vector de trabajo
!        allocate( workSpace( lengthWorkSpace ) )

!        if( present( eigenVectors ) ) then

!           if (present ( dm ) ) then
!              eigenvectors%values = m
!           else	

!              eigenVectors%values=this%values

!           end if

!           !! Calcula valores propios de la matriz de entrada

!           call eigen_vec( &
!                COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
!                UPPER_TRIANGLE_IS_STORED, &
!                matrixSize, &
!                eigenVectors%values, &
!                matrixSize, &
!                vectorAux, &
!                workSpace, &
!                lengthWorkSpace, &
!                infoProcess )

!           eigenValues%values=vectorAux
!           deallocate(vectorAux)

!        else
!           !! Crea la matriz que almacenara los vectores propios
!           call Matrix_copyConstructor( eigenVectorsTmp, this )

!           !! Calcula valores propios de la matriz de entrada
!           call dsyev( &
!                COMPUTE_EIGENVALUES, &
!                UPPER_TRIANGLE_IS_STORED, &
!                matrixSize, &
!                eigenVectorsTmp%values, &
!                matrixSize, &
!                eigenValues%values, &
!                workSpace, &
!                lengthWorkSpace, &
!                infoProcess )

!           call Matrix_destructor( eigenVectorsTmp )

!        end if

!        !! Determina la ocurrencia de errores
!        if ( infoProcess /= 0 )  then

!           call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")

!        end if

!        do i=1,size(eigenValues%values)
!           if( eigenValues%values(i) == Math_NaN ) then

!              call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")
!           end if
!        end do


!        !! libera memoria separada para vector de trabajo
!        deallocate(workSpace)

!     end if

!   end subroutine Matrix_eigen

! #else

 subroutine Matrix_eigen( this, eigenValues, eigenVectors, flags, m, dm )
    implicit none
    type(Matrix), intent(in) :: this
    type(Vector), intent(inout) :: eigenValues
    type(Matrix), intent(inout), optional :: eigenVectors
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: dm
    real(8), intent(in), optional :: m(:,:)

    integer :: lengthWorkSpace
    integer :: matrixSize
    integer :: infoProcess
    real(8), allocatable :: workSpace(:)
    type(Matrix) :: eigenVectorsTmp
    integer :: i

    
    matrixSize = size( this%values, DIM=1 )

    if( flags == SYMMETRIC ) then

       !! Determina la longitud adecuada del vector de trabajo
       lengthWorkSpace=3*matrixSize-1

       !! Crea el vector de trabajo
       allocate( workSpace( lengthWorkSpace ) )

       if( present( eigenVectors ) ) then

          if (present ( dm ) ) then
             eigenvectors%values = m
          else	

             eigenVectors%values=this%values
             
          end if

          !! Calcula valores propios de la matriz de entrada
          call dsyev( &
               COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
               UPPER_TRIANGLE_IS_STORED, &
               matrixSize, &
               eigenVectors%values, &
               matrixSize, &
               eigenValues%values, &
               workSpace, &
               lengthWorkSpace, &
               infoProcess )

       else
          !! Crea la matriz que almacenara los vectores propios
          call Matrix_copyConstructor( eigenVectorsTmp, this )

          !! Calcula valores propios de la matriz de entrada
          call dsyev( &
               COMPUTE_EIGENVALUES, &
               UPPER_TRIANGLE_IS_STORED, &
               matrixSize, &
               eigenVectorsTmp%values, &
               matrixSize, &
               eigenValues%values, &
               workSpace, &
               lengthWorkSpace, &
               infoProcess )

          call Matrix_destructor( eigenVectorsTmp )

       end if

       !! Determina la ocurrencia de errores
       if ( infoProcess /= 0 )  then

          call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")
          print *, "Info Process: ", infoProcess

       end if

       do i=1,size(eigenValues%values)
          if( eigenValues%values(i) == Math_NaN ) then

             call Matrix_exception(WARNING, "Diagonalization failed, Math_NaN", "Class object Matrix in the getEigen() function")
          end if
       end do


       !! libera memoria separada para vector de trabajo
       deallocate(workSpace)

    end if

  end subroutine Matrix_eigen

 subroutine Matrix_eigen2stage( this, eigenValues, eigenVectors, flags, m, dm )
    implicit none
    type(Matrix), intent(in) :: this
    type(Vector), intent(inout) :: eigenValues
    type(Matrix), intent(inout), optional :: eigenVectors
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: dm
    real(8), intent(in), optional :: m(:,:)

    integer :: lengthWorkSpace
    integer :: matrixSize
    integer :: infoProcess
    real(8), allocatable :: workSpace(:)
    type(Matrix) :: eigenVectorsTmp
    integer :: i

    
    matrixSize = size( this%values, DIM=1 )

    if( flags == SYMMETRIC ) then

       !! Determina la longitud adecuada del vector de trabajo
       lengthWorkSpace=3*matrixSize-1

       !! Crea el vector de trabajo
       allocate( workSpace( lengthWorkSpace ) )

       if( present( eigenVectors ) ) then

          if (present ( dm ) ) then
             eigenvectors%values = m
          else	

             eigenVectors%values=this%values
             
          end if

          !! Calcula valores propios de la matriz de entrada
          call dsyev_2stage( &
               COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
               UPPER_TRIANGLE_IS_STORED, &
               matrixSize, &
               eigenVectors%values, &
               matrixSize, &
               eigenValues%values, &
               workSpace, &
               lengthWorkSpace, &
               infoProcess )

       else
          !! Crea la matriz que almacenara los vectores propios
          call Matrix_copyConstructor( eigenVectorsTmp, this )

          !! Calcula valores propios de la matriz de entrada
          call dsyev_2stage( &
               COMPUTE_EIGENVALUES, &
               UPPER_TRIANGLE_IS_STORED, &
               matrixSize, &
               eigenVectorsTmp%values, &
               matrixSize, &
               eigenValues%values, &
               workSpace, &
               lengthWorkSpace, &
               infoProcess )

          call Matrix_destructor( eigenVectorsTmp )

       end if

       !! Determina la ocurrencia de errores
       if ( infoProcess /= 0 )  then

          call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")
          print *, "Info Process: ", infoProcess

       end if

       do i=1,size(eigenValues%values)
          if( eigenValues%values(i) == Math_NaN ) then

             call Matrix_exception(WARNING, "Diagonalization failed, Math_NaN", "Class object Matrix in the getEigen() function")
          end if
       end do


       !! libera memoria separada para vector de trabajo
       deallocate(workSpace)

    end if

  end subroutine Matrix_eigen2stage
  
!>
!! @brief Matrix_eigen_select
!!  -- LAPACK driver routine (version 3.2) --
!!  DSYEVX computes selected eigenvalues and, optionally, eigenvectors
!!  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
!!  selected by specifying either a range of values or a range of indices
!!  for the desired eigenvalues. 
!! For further information, visit http://www.netlib.org/lapack/double/dsyevx.f
!! @warning Implemented only to compute eigenvalues. 

  subroutine Matrix_eigen_select( this, eigenValues, smallestEigenValue, largestEigenValue, eigenVectors, flags, m, dm, method )
    implicit none
    type(Matrix), intent(in) :: this
    type(Vector8), intent(inout) :: eigenValues
    type(Matrix), intent(inout), optional :: eigenVectors
    integer(4), intent(in) :: smallestEigenValue, largestEigenValue  !! The indices (in ascending order)
    !! of the eigenvalues to be computed. Start at 1
    integer(4), intent(in), optional :: flags
    integer(4), intent(in), optional :: dm
    real(8), intent(in), optional :: m(:,:)
    integer(4), intent(in), optional :: method
    integer(4) :: lengthWorkSpace
    integer(4) :: matrixSize
    integer(4) :: infoProcess
    real(8), allocatable :: workSpace(:)
    type(Matrix) :: eigenVectorsTmp
    integer(4) :: i
    real(8) :: vl, vu  !! If RANGE='V', the lower and upper bounds of the interval to be searched for eigenvalues.
    real(8) :: abstol
    integer(4) :: m_dsyevx   !! The total number of eigenvalues found.
    integer(4), allocatable :: iFail(:), iwork(:)

    !!Negative ABSTOL means using the default value
    abstol = -1.0

    matrixSize = size( this%values, DIM=1 )
    if (allocated (iFail) ) deallocate (iFail)
    allocate (iFail(matrixSize))

    if (allocated (iwork) ) deallocate (iwork)
    allocate( iwork( 5*matrixSize ) )

    call omp_set_num_threads(omp_get_max_threads())
!    call omp_set_num_threads (OMP_GET_NUM_THREADS())

    if( flags == SYMMETRIC ) then

      !! Determina la longitud adecuada del vector de trabajo
      lengthWorkSpace=3*matrixSize-1

      !! Crea el vector de trabajo
      allocate( workSpace( lengthWorkSpace ) )

      if( present( eigenVectors ) ) then

         if (present ( dm ) ) then
            eigenvectors%values = m
         else

            eigenVectors%values=this%values

         end if

         !! Crea la matriz que almacenara los vectores propios
         !! call Matrix_copyConstructor( eigenVectorsTmp, this )

         lengthWorkSpace = -1
         !! calculates the optimal size of the WORK array
         call dsyevx( &
              COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectors%values, &
              matrixSize, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              ifail, infoProcess )

         lengthWorkSpace = int(workSpace(1))

         !! Crea el vector de trabajo
         if (allocated(workSpace)) deallocate(workSpace)
         allocate( workSpace( lengthWorkSpace ) )

         !! Calcula valores propios de la matriz de entrada
         call dsyevx( &
              COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectors%values, &
              matrixSize, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              ifail, infoProcess )

       else !! only eigenvalues
       
         !! Crea la matriz que almacenara los vectores propios
         call Matrix_copyConstructor( eigenVectorsTmp, this )

         lengthWorkSpace = -1
         !! calculates the optimal size of the WORK array
         call dsyevx( &
              COMPUTE_EIGENVALUES, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectorsTmp%values, &
              matrixSize, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              ifail, infoProcess )

         lengthWorkSpace = int(workSpace(1))

         !! Crea el vector de trabajo
         if (allocated(workSpace)) deallocate(workSpace)
         allocate( workSpace( lengthWorkSpace ) )

         !! Calcula valores propios de la matriz de entrada
         call dsyevx( &
              COMPUTE_EIGENVALUES, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectorsTmp%values, &
              matrixSize, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              ifail, infoProcess )

         call Matrix_destructor( eigenVectorsTmp )
       
       end if !!  present( eigenVectors ) 

      !! Determina la ocurrencia de errores
      if ( infoProcess /= 0 )  then

         call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")

      end if

      do i=1,size(eigenValues%values)
        if( eigenValues%values(i) == Math_NaN ) then

         call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")
         end if
      end do


      !! libera memoria separada para vector de trabajo
      deallocate(workSpace)

    end if
    
  end subroutine Matrix_eigen_select

  subroutine Matrix_eigen_dsyevr( this, eigenValues, smallestEigenValue, largestEigenValue, eigenVectors, flags, m, dm, method )
    implicit none
    type(Matrix), intent(in) :: this
    type(Vector8), intent(inout) :: eigenValues
    type(Matrix), intent(inout), optional :: eigenVectors
    integer, intent(in) :: smallestEigenValue, largestEigenValue  !! The indices (in ascending order)
    !! of the eigenvalues to be computed. Start at 1
    integer, intent(in), optional :: flags
    integer, intent(in), optional :: dm
    real(8), intent(in), optional :: m(:,:)
    integer, intent(in), optional :: method
    integer :: lengthWorkSpace, lengthiwork
    integer :: matrixSize
    integer :: infoProcess
    real(8), allocatable :: workSpace(:)
    type(Matrix) :: eigenVectorsTmp
    integer :: i
    real(8) :: vl, vu  !! If RANGE='V', the lower and upper bounds of the interval to be searched for eigenvalues.
    real(8) :: abstol
    integer :: m_dsyevx   !! The total number of eigenvalues found.
    integer, allocatable :: iFail(:), iwork(:)

    !!Negative ABSTOL means using the default value
    abstol = -1.0

    matrixSize = size( this%values, DIM=1 )
    if (allocated (iFail) ) deallocate (iFail)
    allocate (iFail(matrixSize))

    if (allocated (iwork) ) deallocate (iwork)
    allocate( iwork( matrixSize*5 ) )

    call omp_set_num_threads(omp_get_max_threads())
!    call omp_set_num_threads (OMP_GET_NUM_THREADS())

    if( flags == SYMMETRIC ) then

      !! Determina la longitud adecuada del vector de trabajo
      lengthWorkSpace=3*matrixSize-1

      !! Crea el vector de trabajo
      allocate( workSpace( lengthWorkSpace ) )

      if( present( eigenVectors ) ) then

         if (present ( dm ) ) then
            eigenvectors%values = m
         else

            eigenVectors%values=this%values

         end if
         
         lengthWorkSpace = -1
         lengthiwork = -1

         !! calculates the optimal size of the WORK array
         call dsyevr( &
              COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectors%values, matrixSize, &
              largestEigenValue - smallestEigenValue + 1, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              lengthiwork, &
              infoProcess )

         lengthWorkSpace = int(workSpace(1))
         lengthiwork = int(iwork(1))
         
         !! Crea el vector de trabajo
         if (allocated(workSpace)) deallocate(workSpace)
         allocate( workSpace( lengthWorkSpace ) )

         !! Crea el vector de trabajo
         if (allocated(iwork)) deallocate(iwork)
         allocate( iwork( lengthiwork ) )


         !! Calcula valores propios de la matriz de entrada
         call dsyevr( &
              COMPUTE_EIGENVALUES_AND_EIGENVECTORS, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectors%values, matrixSize, &
              largestEigenValue - smallestEigenValue + 1, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              lengthiwork, &
              infoProcess )

      else

         !! Crea la matriz que almacenara los vectores propios
         call Matrix_copyConstructor( eigenVectorsTmp, this )

         lengthWorkSpace = -1
         lengthiwork = -1

         !! calculates the optimal size of the WORK array
         call dsyevr( &
              COMPUTE_EIGENVALUES, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectorsTmp%values, &
              matrixSize, &
              largestEigenValue - smallestEigenValue + 1, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              lengthiwork, &
              infoProcess )

         lengthWorkSpace = int(workSpace(1))
         lengthiwork = int(iwork(1))

         !! Crea el vector de trabajo
         if (allocated(workSpace)) deallocate(workSpace)
         allocate( workSpace( lengthWorkSpace ) )

         !! Crea el vector de trabajo
         if (allocated(iwork)) deallocate(iwork)
         allocate( iwork( lengthiwork ) )


         !! Calcula valores propios de la matriz de entrada
         call dsyevr( &
              COMPUTE_EIGENVALUES, &
              "I", &
              UPPER_TRIANGLE_IS_STORED, &
              matrixSize, &
              this%values, &
              matrixSize, &
              vl, vu, &
              smallestEigenValue, largestEigenValue, &
              abstol, &
              m_dsyevx, &
              eigenValues%values, &
              eigenVectorsTmp%values, &
              matrixSize, &
              largestEigenValue - smallestEigenValue + 1, &
              workSpace, &
              lengthWorkSpace, &
              iwork, &
              lengthiwork, &
              infoProcess )
            
        call Matrix_destructor( eigenVectorsTmp )
             
      end if

      !! Determina la ocurrencia de errores
      if ( infoProcess /= 0 )  then

         call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")

      end if

      do i=1,size(eigenValues%values)
        if( eigenValues%values(i) == Math_NaN ) then

            call Matrix_exception(WARNING, "Diagonalization failed", "Class object Matrix in the getEigen() function")
        end if
      end do

      !! libera memoria separada para vector de trabajo
      deallocate(workSpace)

    end if

  end subroutine Matrix_eigen_dsyevr

! #endif

  !>
  !! @brief  Calcula la descomposicion en valores simples de la matriz especificada
  subroutine Matrix_svd( this, basisOfRange, basisOfNullSpace, singularValues )
    implicit none
    type(Matrix), intent(in) :: this
    type(Matrix), intent(inout) :: basisOfRange
    type(Matrix), intent(inout) :: basisOfNullSpace
    type(Matrix), intent(inout) :: singularValues

    real(8), allocatable :: singularValuesVector(:)
    real(8), allocatable :: workSpace(:)
    real(8) :: dummy(1,1)
    integer :: lengthWorkSpace
    integer :: numberOfRows
    integer :: numberOfColumns
    integer :: infoProcess
    integer :: i

    basisOfRange = this
    numberOfRows = size(this%values, dim=1)
    numberOfColumns=size(this%values, dim=2)
    lengthWorkSpace = numberOfRows+4*numberOfColumns+64*(numberOfColumns+numberOfRows)
    if ( allocated(singularValues%values) ) deallocate(singularValues%values)
    allocate( singularValues%values(min(numberOfRows,numberOfColumns),min(numberOfRows,numberOfColumns)) )
    if ( allocated(basisOfNullSpace%values) ) deallocate(basisOfNullSpace%values)
    allocate( basisOfNullSpace%values( min(numberOfColumns,numberOfRows),numberOfColumns) )

    !! Crea el vector de trabajo
    allocate( workSpace( lengthWorkSpace ) )
    allocate( singularValuesVector(min( numberOfRows, numberOfColumns)) )

    call  dgesvd( &
         'O',&
         'S',&
         numberOfRows, &
         numberOfColumns,&
         basisOfRange%values, &
         numberOfRows, &
         singularValuesVector, &
         dummy, &
         1, &
         basisOfNullSpace%values, &
         min(numberOfColumns,numberOfRows), &
         workSpace, &
         lengthWorkSpace, &
         infoProcess )

    singularValues%values = 0.0_8
    do i=1,size(singularValuesVector)
       singularValues%values(i,i) = singularValuesVector(i)
    end do

    !! libera memoria separada para vector de trabajo
    deallocate( workSpace )
    deallocate( singularValuesVector)

  end subroutine Matrix_svd

  !>
  !! @brief Simetriza una matriz triangular inferior
  !! @param flags Indica si la matriz de entrada es triangular inferior (L) o triangular superior(U)
  subroutine Matrix_symmetrize( this, flags )
    implicit none
    type(Matrix) :: this
    character(*) :: flags

    integer :: i
    integer :: numberOfRows
    integer :: numberOfColumns

    numberOfRows = size(this%values,dim=1)
    numberOfColumns = size(this%values,dim=2)

    select case( trim(flags) )

    case("U") !!Copia la mitad superior en la mitad inferior

       do i=1, numberOfColumns
          this%values(i+1:numberOfRows,i)=this%values(i,i+1:numberOfColumns)
       end do

    case("L") !! Copia la mitad inferior en la mitad superior

       do i=1, numberOfRows
          this%values(i,i+1:numberOfColumns)=this%values(i+1:numberOfRows,i)
       end do

    end select

  end subroutine Matrix_symmetrize

  !>
  !! @brief  resulve un sistema de euaciones lineales
  function Matrix_solveLinearEquation( this, rightSideCoefficients ) result( output )
    implicit none
    type(Matrix), intent(in)  :: this
    real(8), intent(in) :: rightSideCoefficients(:)
    real(8), allocatable :: output(:)

    integer :: matrixSize
    integer :: infoProcess
    real(8), allocatable :: workSpace(:)
    integer, allocatable :: pivotIndices(:)
    integer :: i
    integer :: j
    type(Matrix) :: auxMatrix

    !! Detemina variables y parametros requeridos para el calculo
    matrixSize = size( this%values, DIM=1 )

    allocate( pivotIndices( matrixSize ))
    allocate( workSpace( matrixSize ) )


    if ( allocated(output) ) deallocate( output)
    allocate( output(matrixSize) )

    call Matrix_copyConstructor( auxMatrix, this )

    auxMatrix = Matrix_factorizeLU( auxMatrix, pivotIndices=pivotIndices )

    !! Invierte la matriz de entrada
    call dgetri( &
         matrixSize, &
         auxMatrix%values, &
         matrixSize, &
         pivotIndices, &
         workSpace, &
         matrixSize, &
         infoProcess )

    output = 0.0_8
    do i=1, matrixSize
       do j=1, matrixSize
          output(i)= output(i) + rightSideCoefficients(j) * auxMatrix%values(i,j)
       end do
    end do

    !! Determina la ocurrencia de errores
    if ( infoProcess /= 0 )  then

       call Matrix_exception(WARNING, "Procedure failed", "Class object Matrix in the solveLinearEquation() function")

    end if

    !! libera memoria
    deallocate(workSpace)
    deallocate(pivotIndices)

  end function Matrix_solveLinearEquation

  !>
  !! @brief  Retorna true si la matriz tiene todos sus elementos iguales
  !!               a cero, de lo contrario false
  !! @todo Falta implementar
  function Matrix_isNull( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    logical :: output

    output = .false.

  end function Matrix_isNull

  !>
  !! @brief Retorna la transpuesta de la matriz
  !! @author F.M. 2025
  function Matrix_getTranspose(this) result(output)
    implicit none
    type(Matrix), intent(in) :: this
    type(Matrix) :: output
    integer :: i,j

    
    call Matrix_constructor(output, int(size(this%values,DIM=2),8), int(size(this%values,DIM=1),8))

    do i=1, size(this%values,DIM=2)
       do j=1, size(this%values,DIM=1)
          output%values(i,j) = this%values(j,i)
       end do
    end do

  end function Matrix_getTranspose

  !>
  !! @brief Factoriza la matriz
  !! @todo Hay que tener cuidado con el valor de numberOfRows y numberOfColumns
  !!       ya que en el caso de matrices cuadradas no hay problema ( DIM=1 o DIM=2 ?)
  function Matrix_factorizeLU( this, L, U, pivotIndices, printFormatFlags ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output
    type(Matrix), intent(inout), optional :: L
    type(Matrix), intent(inout), optional :: U
    integer, allocatable, intent(inout), optional :: pivotIndices(:)
    integer,intent(in),optional :: printFormatFlags

    integer :: numberOfRows
    integer :: numberOfColumns
    integer :: infoProcess
    integer, allocatable :: pivotIndicesTmp(:)
    integer :: methodTmp
    integer :: i
    integer :: j
    logical :: printError
    
    printError=.true.
    if(present(printFormatFlags)) then
       if(printFormatFlags/=WITHOUT_MESSAGES) printError=.false.
    end if
    
    !! Determina variables y parametros requeridos para el calculo
    numberOfRows = size( this%values, DIM=1 )
    numberOfColumns = size( this%values, DIM=1 )
    
    call Matrix_copyConstructor( output, this )

    if( .not. present( pivotIndices ) ) then
       allocate( pivotIndicesTmp( min( numberOfRows, numberOfColumns )  )  )

       call dgetrf( &
            numberOfRows, &
            numberOfColumns, &
            output%values, &
            numberOfRows, &
            pivotIndicesTmp, &
            infoProcess )

       deallocate( pivotIndicesTmp )
    else
       call dgetrf( &
            numberOfRows, &
            numberOfColumns, &
            output%values, &
            numberOfRows, &
            pivotIndices, &
            infoProcess )
    end if

    if( present(L) ) then
       call Matrix_setNull( L )

       do i=1, numberOfRows
          do j=1, i-1
             L%values(i, j) = output%values(i, j)
          end do

          L%values(i, i) = 1.0_8
       end do
    end if

    if( present(U) ) then
       call Matrix_setNull( U )

       do i=1, numberOfRows
          do j=i, numberOfColumns
             U%values(i, j) = output%values(i, j)
          end do
       end do
    end if

    !! Determina la ocurrencia de errores
    if ( infoProcess /= 0 .and. printError )  then
       call Matrix_exception(WARNING, "Factorization failed", "Class object Matrix in the factorizeLU() function" )
    end if

  end function Matrix_factorizeLU

  !>
  !! @brief  Retorna la traza de la matriz
  function Matrix_trace( this ) result ( output )
    implicit none
    type(Matrix) , intent(in) :: this
    real(8) :: output

    integer :: i

    output = 0.0_8

    do i = 1, size( this%values, DIM=1 )

       output = output + this%values( i, i )

    end do

  end function Matrix_trace

  !>
  !! @brief perform one of the matrix-matrix operations
  !! @author J.M. Rodas, 2015
  !! C = alpha*op( A )*op( B ) + beta*C, using DGEMM of Lapack 
  function Matrix_multiplication( TRANSA, TRANSB, M, N, K, ALPHA, matrixA, LDA, matrixB, LDB, BETA, LDC ) result ( matrixC )
    implicit none
    type(Matrix), intent(in) :: matrixA 
    type(Matrix), intent(in) :: matrixB
    character, intent(in) :: TRANSA, TRANSB
    integer, intent(in) :: M, N, K, LDA, LDB, LDC
    real(8), intent(in) :: ALPHA, BETA
    type(Matrix) :: matrixC
    integer :: i

    call Matrix_constructor(matrixC, int(LDC,8), int(N,8))

    call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, matrixA%values, LDA, matrixB%values, LDB, BETA, matrixC%values, LDC)
   
  end function Matrix_multiplication

  !>
  !! @brief Multiplica dos matrices, funciona para matrices grandes
  !! @author F.M. 2025
  function Matrix_product_dgemm( this, otherThis ) result ( output )
    implicit none
    type(Matrix), intent(in) :: this
    type(Matrix), intent(in) :: otherThis
    type(Matrix) :: output

    integer :: M, N, K
    real(8) :: ALPHA, BETA

    M=size( this%values , DIM=1 )
    N=size( otherThis%values , DIM=2 )
    K=size( this%values , DIM=2 )
  
    call Matrix_constructor(output, int(M,8), int(N,8))

    call dgemm("N", "N", M, N, K, 1.0_8, this%values, M, otherThis%values, K, 0.0_8, output%values, M)
		
  end function Matrix_product_dgemm

  !>
  !! @brief Multiplica dos matrices
  function Matrix_product( this, otherMatrix ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix), intent(in) :: otherMatrix
    type(Matrix) :: output

    call Matrix_constructor( output, int(size(this%values, DIM=1),8), int(size(otherMatrix%values, DIM=2),8) )
    output%values = matmul( this%values, otherMatrix%values )
		
  end function Matrix_product

  !>
  !! @brief Suma dos matrices
  function Matrix_plus( this, otherMatrix ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix), intent(in) :: otherMatrix
    type(Matrix) :: output

    call Matrix_constructor( output, int( size(this%values, DIM=1), 8 ), int( size(otherMatrix%values, DIM=2), 8 ) )
    output%values = this%values + otherMatrix%values

  end function Matrix_plus

  !>
  !! @brief  Calcula una funcion generica de la matriz this
  function Matrix_function( this, func, param ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    character(*), intent(in) :: func
    type(Matrix) :: output
    real(8), intent(in), optional :: param


    type(Vector) :: eigenValues
    type(Matrix) :: eigenVectors
    type(Matrix) :: eigenVectorsInverted
    type(Matrix) :: diagonal
    integer :: i

    call Vector_constructor( eigenValues, Matrix_getNumberOfColumns( this ) )
    call Matrix_copyConstructor( eigenVectors, this )
    call Matrix_copyConstructor( eigenVectorsInverted, this )

    call Matrix_eigen( this, eigenValues, eigenVectors=eigenVectors, flags=SYMMETRIC )
    call Matrix_diagonalConstructor( diagonal, eigenValues )

    eigenVectorsInverted = Matrix_inverse( eigenVectors )

    do i=1, size( diagonal%values, DIM=1 )

       select case ( trim(func) )

       case ( "pow" )
          diagonal%values(i, i) = diagonal%values(i, i)**param

       case ( "sqrt" )
          diagonal%values(i, i) = sqrt( diagonal%values(i, i) )

       case ( "log" )
          diagonal%values(i, i) = log( diagonal%values(i, i) )

       case ( "log10" )
          diagonal%values(i, i) = log10( diagonal%values(i, i) )

       case ( "sin" )
          diagonal%values(i, i) = sin( diagonal%values(i, i) )

       case ( "cos" )
          diagonal%values(i, i) = cos( diagonal%values(i, i) )

       case ( "tan" )
          diagonal%values(i, i) = tan( diagonal%values(i, i) )

       case ( "asin" )
          diagonal%values(i, i) = asin( diagonal%values(i, i) )

       case ( "acos" )
          diagonal%values(i, i) = acos( diagonal%values(i, i) )

       case ( "atan" )
          diagonal%values(i, i) = atan( diagonal%values(i, i) )

       case ( "sinh" )
          diagonal%values(i, i) = sinh( diagonal%values(i, i) )

       case ( "cosh" )
          diagonal%values(i, i) = cosh( diagonal%values(i, i) )

       case ( "tanh" )
          diagonal%values(i, i) = tanh( diagonal%values(i, i) )

       end select

    end do

    output = Matrix_product( diagonal, eigenVectorsInverted )
    output = Matrix_product( eigenVectors, output )

  end function Matrix_function

  !>
  !! @brief  Calcula una funcion generica de la matriz this
  function Matrix_function_svd( this, func, param ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    character(*), intent(in) :: func
    type(Matrix) :: output
    real(8), intent(in), optional :: param


    type(Matrix) :: U, VT, singular
    integer :: i
    integer :: dim

    dim=Matrix_getNumberOfColumns(this)
    
    call Matrix_constructor(U, int(dim,8), int(dim,8), 0.0_8)
    call Matrix_constructor(VT, int(dim,8), int(dim,8), 0.0_8)
    call Matrix_constructor(singular, int(dim,8), int(dim,8), 0.0_8)
    call Matrix_svd( this, U, VT, singular )
    
    ! vectorsInverted = Matrix_inverse( range )

    do i=1, size( singular%values, DIM=1 )

       select case ( trim(func) )

       case ( "pow" )
          singular%values(i, i) = singular%values(i, i)**param

       case ( "sqrt" )
          singular%values(i, i) = sqrt( singular%values(i, i) )

       case ( "log" )
          singular%values(i, i) = log( singular%values(i, i) )

       case ( "log10" )
          singular%values(i, i) = log10( singular%values(i, i) )

       case ( "sin" )
          singular%values(i, i) = sin( singular%values(i, i) )

       case ( "cos" )
          singular%values(i, i) = cos( singular%values(i, i) )

       case ( "tan" )
          singular%values(i, i) = tan( singular%values(i, i) )

       case ( "asin" )
          singular%values(i, i) = asin( singular%values(i, i) )

       case ( "acos" )
          singular%values(i, i) = acos( singular%values(i, i) )

       case ( "atan" )
          singular%values(i, i) = atan( singular%values(i, i) )

       case ( "sinh" )
          singular%values(i, i) = sinh( singular%values(i, i) )

       case ( "cosh" )
          singular%values(i, i) = cosh( singular%values(i, i) )

       case ( "tanh" )
          singular%values(i, i) = tanh( singular%values(i, i) )

       end select

    end do

    output = Matrix_product_dgemm( singular, VT )
    output = Matrix_product_dgemm( U, output )

  end function Matrix_function_svd
  
  !>
  !! @brief  Calcula la potencia de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_pow( this, eexponent, method ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    real(8), intent(in) :: eexponent
    character(*), optional :: method    
    type(Matrix) :: output
    character(10) :: selectedMethod
    
    if( present ( method ) ) then
       selectedMethod=trim(method)
    else
       selectedMethod="eigen"
    end if

    select case( trim(selectedMethod) )
    case("SVD")
       output = Matrix_function_svd( this, "pow", eexponent )

    case("eigen")
       output = Matrix_function( this, "pow", eexponent )

    case default
       call Matrix_exception(ERROR, "The selected method to compute the matrix power is not implemented", "Class object Matrix in the Matrix_pow() function")
    end select

  end function Matrix_pow

  !>
  !! @brief  Calcula la raiz cuadrada de la matriz this
  !! @warning La matriz tiene que ser diagonalizable y tener valores propios mayores que cero
  function Matrix_sqrt( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "sqrt" )

  end function Matrix_sqrt

  !>
  !! @brief  Calcula el logaritmo natural de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_log( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "log" )

  end function Matrix_log

  !>
  !! @brief  Calcula el logaritmo en base 10 de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_log10( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "log10" )

  end function Matrix_log10

  !>
  !! @brief Calcula la funcion seno de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_sin( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "sin" )

  end function Matrix_sin

  !>
  !! @brief  Calcula la funcion coseno de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_cos( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "cos" )

  end function Matrix_cos

  !>
  !! @brief  Calcula la funcion tangente de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_tan( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "tan" )

  end function Matrix_tan

  !>
  !! @brief  Calcula la funcion arcoseno de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_asin( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "asin" )

  end function Matrix_asin

  !>
  !! @brief  Calcula la funcion arcocoseno de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_acos( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "acos" )

  end function Matrix_acos

  !>
  !! @brief  Calcula la funcion arcotangente de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_atan( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "atan" )

  end function Matrix_atan

  !>
  !! @brief Calcula la funcion seno hiperbolico de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_sinh( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "sinh" )

  end function Matrix_sinh

  !>
  !! @brief  Calcula la funcion coseno hiperbolico de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_cosh( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "cosh" )

  end function Matrix_cosh

  !>
  !! @brief  Calcula la fincion tangente hiperbolico de la matriz this
  !! @warning La matriz tiene que ser diagonalizable
  function Matrix_tanh( this ) result ( output )
    implicit none
    type(Matrix), intent(inout) :: this
    type(Matrix) :: output

    output = Matrix_function( this, "tanh" )

  end function Matrix_tanh

  !>
  !! @todo Debe adicionarse la posiblidad de tratar matrices no simetricas y rectangulares
  subroutine Matrix_eigenProperties(this)
    implicit none
    type(Matrix) :: this

    integer :: rows
    integer :: columns
    type(Vector) :: eigenValues
    type(Matrix) :: eigenVectors

    rows = size( this%values, dim = 1 )
    columns = size( this%values, dim = 2 )

    if( columns == rows ) then

       if ( .not.allocated(this%eigenValues) ) then

          allocate( this%eigenVectors(rows,rows))
          this%eigenVectors = 0.0_8

          allocate( this%eigenValues(rows))
          this%eigenValues = 0.0_8

       else &
            if ( allocated(this%eigenValues) .and. ( size(this%eigenValues) /= rows ) ) then

          deallocate(this%eigenVectors)
          allocate( this%eigenVectors(rows,rows))
          this%eigenVectors = 0.0_8

          deallocate(this%eigenValues)
          allocate( this%eigenValues(rows))
          this%eigenValues = 0.0_8

       end if


       call Vector_constructor( eigenValues, rows)
       call Matrix_constructor( eigenVectors, int(rows,8), int(rows,8) )


       !! Obtiene vectores y valores propios de la matriz de entrada
       call Matrix_eigen( this, eigenValues, eigenvectors, SYMMETRIC ) !! Verificar que esta matriz sea simetrica

       this%eigenVectors=eigenVectors%values
       this%eigenValues=eigenValues%values

       call Vector_destructor( eigenValues)
       call Matrix_destructor( eigenVectors )

    else
       call Matrix_exception(ERROR, "The matrix is not square", "Class object Matrix in the eigenProperties() function")
    end if

  end subroutine Matrix_eigenProperties

  !>
  !! @brief Indica si la matriz es definida positiva
  function Matrix_isPositiveDefinited(this, eigenValueBound ) result(output)
    implicit none
    type(Matrix) :: this
    logical :: output
    real(8), optional :: eigenValueBound

    integer :: i
    real(8) :: auxEigenValueBound

    if( .not.allocated(this%eigenValues) ) call Matrix_eigenProperties(this)

    !! Acota los valores propios de la hessiana
    auxEigenValueBound = CONTROL_instance%DOUBLE_ZERO_THRESHOLD
    if(present(eigenValueBound)) auxEigenValueBound= eigenValueBound

    do i=1,size(this%eigenValues)
       if( abs( this%eigenValues(i) ) < auxEigenValueBound ) this%eigenValues(i)=0.0
    end do

    this%numberOfNegatives = 0
    do i=1,size(this%eigenValues)
       if ( ( this%eigenValues(i) < 0.0_8 )  ) this%numberOfNegatives = this%numberOfNegatives + 1
    end do

    this%isPositiveDefinited=.true.
    if ( this%numberOfNegatives > 0 ) 	this%isPositiveDefinited=.false.

    output = this%isPositiveDefinited

  end function Matrix_isPositiveDefinited

  !>
  !! @brief  reviza el espectro de valores propios de la matriz especificada
  !!	(asumida simetrica) acontando sus valores entre las cotas especidicadas
  !! @param lowerBond Cota inferior
  !! @param lowerBond Cota superior
  !! @warning La matriz tiene que ser diagonalizable
  subroutine Matrix_boundEigenValues( this, lowerBond, upperBond )
    implicit none
    type(Matrix) :: this
    real(8), intent(in) :: lowerBond
    real(8), intent(in) :: upperBond

    integer :: i

    if( .not.allocated(this%eigenValues) ) call Matrix_eigenProperties(this)

    !!***********************************************************
    !! Acota los valores propios entre el rango establecido
    !!**
    do i=1,size( this%eigenValues )
       if ( abs( this%eigenValues(i) ) < lowerBond ) then
          this%eigenValues(i) = sign( lowerBond,this%eigenValues(i) )
       else &
            if ( abs( this%eigenValues(i) ) > upperBond ) then
          this%eigenValues(i) = sign( upperBond,this%eigenValues(i) )
       end if
    end do

  end subroutine Matrix_boundEigenValues

  !>
  !! @brief Retorna los valores propios de la matriz especificada
  function Matrix_getEigenValues(this) result(output)
    implicit none
    type(Matrix) :: this
    type(Vector) :: output

    if( .not.allocated(this%eigenValues) ) call Matrix_eigenProperties(this)
    call Vector_constructor(output,size(this%eigenValues))
    output%values=this%eigenValues

  end function Matrix_getEigenValues

  !>
  !! @brief Retorna los vectores propios de la matriz especificada
  function Matrix_getEigenVectors(this) result(output)
    implicit none
    type(Matrix) :: this
    type(Matrix) :: output

    if( .not.allocated(this%eigenValues) ) call Matrix_eigenProperties(this)
    call Matrix_constructor(output,int(size(this%eigenValues),8),int(size(this%eigenValues),8))
    output%values=this%eigenVectors

  end function Matrix_getEigenVectors

  !>
  !! @brief  Devuelve la desviacion estandar de los elementos de un par de matrices
  function Matrix_standardDeviation( this , otherThis ) result ( output)
    implicit none
    type(Matrix), intent(in) :: this
    type(Matrix), intent(in) :: otherThis
    real(8) :: output

    integer :: i
    integer :: j
    integer :: orderMatrix

    orderMatrix = size( this%values , dim = 1 )
    output = 0.0_8

    !!****************************************************************
    !! Obtiene la desviacion estandar de los elementos de la matriz
    !!
    do i = 1 , orderMatrix
       do j = 1 , orderMatrix
          output = output + ( this%values(i,j) - otherthis%values( i,j ) ) **  2.0_8
       end do
    end do

    output = dabs (dsqrt ( ( orderMatrix **(-2.0_8) ) * output ) )

    !!****************************************************************

  end function Matrix_standardDeviation

  !>
  !! @brief Adiciona una fila al final de la matriz y la inicializa con cero
  subroutine Matrix_addRow( this )
    implicit none
    type(Matrix) :: this

    real(8), allocatable :: auxArray(:,:)
    integer :: rows
    integer :: columns

    columns = size( this%values, dim=2 )
    rows = size( this%values, dim=1 )

    allocate( auxArray(rows+1,columns) )
    auxArray(1:rows,:) = this%values
    auxArray(rows+1,:) = 0.0_8
    deallocate( this%values )
    allocate( this%values(rows+1,columns) )
    this%values = auxArray
    deallocate( auxArray )

  end subroutine Matrix_addRow


  !>
  !! @brief Remueve la fila especificada de una matriz
  subroutine Matrix_removeRow( this, numberOfRow )
    implicit none
    type(Matrix) :: this
    integer, intent(in) :: numberOfRow

    real(8), allocatable :: auxArray(:,:)
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

  end subroutine Matrix_removeRow

  !>
  !! @brief Adiciona una columna al final de la matriz y la inicializa con cero
  subroutine Matrix_addColumn( this )
    implicit none
    type(Matrix) :: this

    real(8), allocatable :: auxArray(:,:)
    integer :: rows
    integer :: columns

    columns = size( this%values, dim=2 )
    rows = size( this%values, dim=1 )

    allocate( auxArray(rows,columns+1) )
    auxArray(:,1:columns) = this%values
    auxArray(:,columns + 1) = 0.0_8
    deallocate( this%values )
    allocate( this%values(rows,columns+1) )
    this%values = auxArray
    deallocate( auxArray )

  end subroutine Matrix_addColumn

  !>
  !! @brief Remueve la columna especificada de una matriz
  subroutine Matrix_removeColumn( this, numberOfColumn )
    implicit none
    type(Matrix) :: this
    integer, intent(in) :: numberOfColumn

    real(8), allocatable :: auxArray(:,:)
    integer :: rows
    integer :: columns

    columns = size( this%values, dim=2 )

    if (numberOfColumn <= columns ) then

       rows = size( this%values, dim=1 )


       allocate( auxArray(rows,columns-1) )
       auxArray(:,1:numberOfColumn-1) = this%values(:,1:numberOfColumn-1)
       auxArray(:,numberOfColumn:columns-1) = this%values(:,numberOfColumn+1:columns)
       deallocate( this%values )
       allocate( this%values(rows,columns-1) )
       this%values = auxArray
       deallocate( auxArray )

    end if

  end subroutine Matrix_removeColumn

  subroutine diagonalize_matrix(m,eig,dm)
  ! Diagonalize
  ! Roberto Flores-Moreno, Aug 2008
  implicit none
    integer dm
    real(8) m(dm,dm),eig(dm)

    type(Matrix) :: tmp1
    type(Matrix) :: eigenVectors
    type(Vector) :: eigenValues

    call Matrix_constructor( tmp1, int(dm,8), int(dm,8 ) )
    call Matrix_constructor( eigenVectors, int(dm,8), int(dm,8 ) )
    call Vector_constructor( eigenValues, dm )

    tmp1%values = m

    call Matrix_eigen(tmp1, eigenValues, eigenVectors, SYMMETRIC, m,dm )
   	m = eigenVectors%values	

    eig(1:dm) = eigenValues%values(1:dm)

    call Matrix_destructor( tmp1 )
    call Matrix_destructor( eigenVectors )
    call Vector_destructor( eigenValues )

  end subroutine



  !>
  !! @brief  resuelve ecuaciones lineales con lapack
  !<
  subroutine Matrix_linear( N, columnsB, A, Asize, B, rowsB, Xsolv, info )
    implicit none
    type(Matrix), intent(in) :: A
    integer, intent(in) :: columnsB
    integer, intent(in) :: N
    integer, intent(in) :: Asize
    real(8), allocatable, intent(in) :: B(:)
    integer, intent(in) :: rowsB
    integer, intent(out) :: info
    real(8), allocatable, intent(out) :: Xsolv(:) 
    integer :: pivot(N)
    type(Matrix) :: Baux
    integer :: i
    integer(8) :: Bsize

    Bsize = rowsB

    call Matrix_constructor( Baux, int(1,8), Bsize)
    
    do i=1, rowsB
       Baux%values(i,1) = B(i)
    end do

    call dgesv( &
         N, &
         columnsB, &
         A%values, &
         ASize, &
         pivot, &
         Baux%values, &
         rowsB, &
         info )

    allocate(Xsolv(rowsB))
    do i=1, rowsB
       Xsolv(i) = Baux%values(i,1)
    end do

  end subroutine Matrix_linear

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine Matrix_exception( typeMessage, description, debugDescription)
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

  end subroutine Matrix_exception



end module Matrix_
