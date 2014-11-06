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
!! @brief Clase encargada de manipular elementos de espacios Vectoriales lineales
!!
!! @author Sergio Gonzalez
!!
!! <b> Fecha de creacion : </b> 2009-05-28
!!   - <tt> 2009-05-28 </tt>: Sergio Gonzalez ( sagonzalez@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
module LinearVectorialSpace_
  use Matrix_
  use Vector_
  use Exception_
  implicit none
  
  type , public :: LinearVectorialSpace
     
     character(30) :: name
     type(Matrix) :: elements
     integer :: numberOfElements
     
  end type LinearVectorialSpace
  
  public :: &
       LinearVectorialSpace_constructor, &
       LinearVectorialSpace_destructor, &
       LinearVectorialSpace_getElements, &
       LinearVectorialSpace_orthogonalize
  
  private
  
contains
  
  !>
  !! @brief Define el constructor para la clase
  subroutine LinearVectorialSpace_constructor(  this, mmatrix, dimensionOfSpace, vectorDimension )
    implicit none
    type(LinearVectorialSpace) :: this
    real(8), optional :: mmatrix(:,:)
    integer, optional :: dimensionOfSpace
    integer, optional :: vectorDimension
    
    this%name = "undefined"
    
    if ( present(mmatrix)  ) then
       call Matrix_constructor(this%elements, int(size(mmatrix, dim=1), 8), int(size(mmatrix, dim=2), 8) )
       this%elements%values= mmatrix
       this%numberOfElements = size(mmatrix, dim=2)
       return
    end if
    
  end subroutine LinearVectorialSpace_constructor
  
  !>
  !! @brief Define el destructor para clase
  !! @param thisPointer Funcion base
  subroutine LinearVectorialSpace_destructor( this )
    implicit none
    type(LinearVectorialSpace) :: this
    
    call Matrix_destructor( this%elements)
    this%numberOfElements = 0
    
  end subroutine LinearVectorialSpace_destructor
  
  
  subroutine LinearVectorialSpace_show( this )
    implicit none
    type(LinearVectorialSpace) :: this
    
    print *,""
    print *, "VECTORIAL SPACE: "
    print *, "================="
    call Matrix_show( this%elements)
    
  end subroutine LinearVectorialSpace_show
  
  function LinearVectorialSpace_getNumberOfElements( this ) result( output )
    implicit none
    type(LinearVectorialSpace) :: this
    integer :: output
    
    output = this%numberOfElements
    
  end function LinearVectorialSpace_getNumberOfElements
  
  function LinearVectorialSpace_getElements(this) result( output )
    implicit none
    type(LinearVectorialSpace) :: this
    type(Matrix) :: output
    
    output = this%elements
    
  end function LinearVectorialSpace_getElements
  
  
  !>
  !! @brief Adiciana un nuevo elemento al espacio vectorial existente
  subroutine LinearVectorialSpace_addElement( this, fortranVector )
    implicit none
    type(LinearVectorialSpace) :: this
    real(8) :: fortranVector(:)
    
    if ( this%numberOfElements == 0 ) then
       call Matrix_constructor( this%elements, int(size(fortranVector),8),1_8 )
       this%elements%values(:,1) = fortranVector
    else
       call Matrix_addColumn( this%elements )
       this%elements%values(:, this%numberOfElements+1) = fortranVector
    end if
    
    this%numberOfElements = this%numberOfElements + 1
    
  end subroutine LinearVectorialSpace_addElement
  
  !>
  !! @brief Expande el espacio vectorial en el numero de elementos especificados
  !!		El espacio resultante se ortogonaliza mediante un proceso Gram-Schmidt
  !! @param this Espacio vectorial
  !! @param expansionSize tama\~no de la expansi\'on
  !! @todo Falta implementacion completa
  subroutine LinearVectorialSpace_expadSpace( this, expansionSize )
    implicit none
    type(LinearVectorialSpace) :: this
    integer :: expansionSize
    
    
  end subroutine LinearVectorialSpace_expadSpace
  
  !>
  !! @brief Proyecta un numero dado de elementos sobre el resto de los elementos del espacio vectorial
  !!		y ortogonaliza el espacio resultante mediante un proceso Gram-Schmidt
  !! @param this Espacio vectorial
  function LinearVectorialSpace_projectLastElements( vectorialSpace ) result(output)
    implicit none
    real(8) :: vectorialSpace(:,:)
    real(8), allocatable :: output(:,:)
    
    integer :: i
    integer :: last
    real(8) :: squareNorm
    real(8) :: projectionOverOrthogonalizedBasis
    
    last = size(vectorialSpace,dim=2)
    allocate( output( size(vectorialSpace,dim=1), last ) )
    output = vectorialSpace
    
    !!***********************************************************************************
    !! Realiza de ortogonalizacion sobre los last-1 vectores, previamente ortogonalizados.
    !!
    do i=1,last-1
       squareNorm = dot_product( output(:,i), output(:,i) )
       
       projectionOverOrthogonalizedBasis=dot_product( output(:,i),vectorialSpace(:,last) )
       
       if ( squareNorm>1.0D-12 ) then
          
          output( :, last ) = output( :, last ) - projectionOverOrthogonalizedBasis/sqrt(squareNorm)*output(:,i)
          
       end if
    end do
    squareNorm = dot_product( output(:,last), output(:,last) )
    output( :, last )=output( :, last )/sqrt(squareNorm)
    
    !!
    !!******************************************************************
    
  end function LinearVectorialSpace_projectLastElements
  
  !>
  !! @brief Ajusta las componentes de un elemento vetorial al valor especificado
  !! @todo falta por implementar
  subroutine LinearVectorialSpace_setElement(this)
    implicit none
    type(LinearVectorialSpace) :: this
    
    
  end subroutine LinearVectorialSpace_setElement
  
  
  !>
  !! @brief Ortogonaliza la componentes de espacio vectorial
  !! @todo falta por implementar
  subroutine LinearVectorialSpace_orthogonalize(this)
    implicit none
    type(LinearVectorialSpace) :: this
    
    integer :: i
    integer :: last
    real(8) :: norm
    
    last = size(this%elements%values,dim=2)
    norm=sqrt(dot_product(this%elements%values(:,1),this%elements%values(:,1)))
    
    this%elements%values(:,1)=this%elements%values(:,1)/norm
    
    !!
    !! Realiza de ortogonalizacion consecutiva de cada uno de los vectores
    !! presentes en la matriz
    !!
    do i=2,last
       
       this%elements%values(:,1:i)=LinearVectorialSpace_projectLastElements( this%elements%values(:,1:i) )
       
    end do
    
    !! Reortonormaliza para asegurar la ortonormalizacion
    do i=2,last
       
       this%elements%values(:,1:i)=LinearVectorialSpace_projectLastElements( this%elements%values(:,1:i) )
       
    end do
    
  end subroutine LinearVectorialSpace_orthogonalize
  
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine LinearVectorialSpace_exception( typeMessage, description, debugDescription)
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
    
  end subroutine LinearVectorialSpace_exception
  
end module LinearVectorialSpace_
