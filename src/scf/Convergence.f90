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
!! @brief Modulo para implementacion de metodos de acelaracion convergencia SCF
!!
!! Este modulo contiene las implementaciones requeridas para realizar
!! procedimientos de aseguramineto o aceleracion de convergencia de metodo SCF.
!! Normalmente estos metodos alteran la matriz de Fock y la matriz de densidad
!! para forzar la convergencia de un procedimineto SCF no convergente.
!!
!! @author Sergio A. Gonzalez
!!
!! <b> Fecha de creacion : </b> 2008-09-20
!!   - <tt> 2007-08-20 </tt>: Sergio A. Gonzalez ( sagonzalezm@unal.edu.co )
!!        -# Creacion del modulo, y metodo "optimal damping"
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# adapta el modulo para su inclusion en Lowdin
!! @todo implementar otros metdos como DIIS y Level shifting  
module Convergence_
  use Exception_
  use Matrix_
  use Vector_
  use CONTROL_
  use MolecularSystem_
  use StackMatrices_
  implicit none
  
  type, public :: Convergence
     
     character(30) :: name
     
     type(Matrix) :: initialDensityMatrix
     type(Matrix) :: initialFockMatrix
     
     type(Matrix), pointer :: newFockMatrixPtr
     type(Matrix), pointer :: newDensityMatrixPtr
     type(Matrix) :: overlapMatrix
     type(Matrix) :: coefficientMatrix
     
     !! Matrices requeidas en DIIS
     type(StackMatrices) :: componentsOfError
     type(StackMatrices) :: fockMatrices
     type(Matrix) :: diisEquationsSystem
     type(Matrix) :: orthonormalizationMatrix
     
     real(8) :: diisError
     integer :: methodType
     integer :: speciesID
     integer :: diisDimensionality
     integer :: currentsize
     logical :: isInstanced
     logical :: isLastComponent
     logical :: hadLowDiisError
     
  end type Convergence
  
  !< enum Matrix_type {
  integer, parameter, public :: SCF_CONVERGENCE_DEFAULT			=  -1
  integer, parameter, public :: SCF_CONVERGENCE_NONE			=  0
  integer, parameter, public :: SCF_CONVERGENCE_DAMPING		=  1
  integer, parameter, public :: SCF_CONVERGENCE_DIIS				=  2
  integer, parameter, public :: SCF_CONVERGENCE_LEVEL_SHIFTING	=  3
  integer, parameter, public :: SCF_CONVERGENCE_MIXED			=  4
  !< }
  
  public ::&
       Convergence_constructor, &
       Convergence_destructor, &
       Convergence_setName, &
       Convergence_setInitialDensityMatrix, &
       Convergence_setInitialFockMatrix, &
       Convergence_getName, &
       Convergence_isInstanced, &
       Convergence_getInitialDensityMatrix, &
       Convergence_getInitialFockMatrix, &
       Convergence_getDiisError, &
       Convergence_setMethodType, &
       Convergence_setMethod, &
       Convergence_run, &
       Convergence_reset
  
  private
  
contains
  
  !>
  !! @brief Define el constructor para la clase
  subroutine Convergence_constructor( this, name ,methodType )
    implicit none
    type(Convergence), intent(inout) :: this
    character(*),optional :: name
    integer, optional :: methodType
    
    this%name = "undefined"
    if ( present(name) ) this%name = trim(name)
    this%methodType = SCF_CONVERGENCE_DEFAULT
    if ( present(methodType) ) this%methodType = methodType
    this%isInstanced = .true.
    this%diisDimensionality = CONTROL_instance%DIIS_DIMENSIONALITY
    this%diisError = 0.0_8
    this%isLastComponent = .false.
    this%hadLowDiisError = .false.
    
  end subroutine Convergence_constructor
  
  !>
  !! @brief Define el destructor para la clase
  subroutine Convergence_destructor( this )
    implicit none
    type(Convergence), intent(inout) :: this
    
    this%newFockMatrixPtr => null()
    this%newDensityMatrixPtr => null()
    if ( allocated(this%initialFockMatrix%values) ) call Matrix_destructor( this%initialFockMatrix)
    if ( allocated(this%initialDensityMatrix%values) ) call Matrix_destructor( this%initialDensityMatrix)
			if ( allocated(this%overlapMatrix%values) ) call Matrix_destructor( this%overlapMatrix)
   if ( StackMatrices_isInstanced( this%fockMatrices ) ) then
      call StackMatrices_destructor( this%fockMatrices)
      call StackMatrices_destructor( this%componentsOfError)
      call Matrix_destructor( this%diisEquationsSystem )
      call Matrix_destructor( this%orthonormalizationMatrix )
   end if
   
   this%methodType =SCF_CONVERGENCE_DEFAULT
   this%isInstanced = .false.
   this%isLastComponent = .false.
   this%hadLowDiisError = .false.
   
 end subroutine Convergence_destructor
 
 !>
 !! @brief Permite asiganar un nombre al metodo
 subroutine Convergence_setName( this, name )
   implicit none
   type(Convergence), intent(inout) :: this
   character(*), intent(in) :: name

   this%name = trim(name)
   
 end subroutine Convergence_setName
 
 !>
 !! @brief Ajusta la dimension del subespacio empleado por metodo DIIS
 subroutine Convergence_setDiisDimensionality( this, diisDimensionality )
   implicit none
   type(Convergence), intent(inout) :: this
   integer, intent(in) :: diisDimensionality
   
   this%diisDimensionality = diisDimensionality
   
 end subroutine Convergence_setDiisDimensionality
 
 !>
 !! @brief Ajusta  la matriz de densidad inicial
 subroutine Convergence_setInitialDensityMatrix( this, otherMatrix )
   implicit none
   type(Convergence), intent(inout) :: this
   type(Matrix), intent(in) :: otherMatrix
   
   if (  this%isInstanced ) then
      
      call Matrix_copyConstructor( this%initialDensityMatrix, otherMatrix )
      
   else
      
      call Convergence_exception( ERROR, "The object has not beeb instanced",&
           "Class object Convergence in setInitialDensityMatrix() function" )
      
   end if
   
 end subroutine Convergence_setInitialDensityMatrix
 
 
 !>
 !! @brief Ajusta la matrix del Fock inicial
 subroutine Convergence_setInitialFockMatrix( this, otherMatrix )
   implicit none
   type(Convergence), intent(inout) :: this
   type(Matrix), intent(in) :: otherMatrix
   
   if ( allocated( this%initialDensityMatrix%values ) ) then
      
      call Matrix_copyConstructor( this%initialFockMatrix, otherMatrix )
      
   else
      
      call Convergence_exception( ERROR, "The object has not beeb instanced",&
           "Class object Convergence in setInitialFockMatrix() function"  )
      
   end if
   
 end subroutine Convergence_setInitialFockMatrix
 
 
 !>
 !! @brief Define el metodo de convergencia SCF que debe emplear
 subroutine Convergence_setMethodType( this,  methodType )
   implicit none
   type(Convergence), intent(inout), target :: this
   integer,optional :: methodType
   
   if ( present(methodType) ) this%methodType = methodType
   
 end subroutine Convergence_setMethodType

 !>
 !! @brief Ajusta los parametros requeridos por el metodo convergencia
 subroutine Convergence_setMethod( this, newFockMatrix, newDensityMatrix, overlapMatrix, &
      methodType, coefficientMatrix, speciesID )
   implicit none
   type(Convergence), intent(inout), target :: this
   type(Matrix),target :: newFockMatrix
   type(Matrix),target :: newDensityMatrix
   type(Matrix),target, optional :: overlapMatrix
   type(Matrix),target, optional :: coefficientMatrix
   integer, optional :: methodType
   integer, optional :: speciesID

   type(Matrix) :: eigenVectors
   type(Vector) :: eigenValues
   integer :: orderOfMatrix
   integer :: i
   integer :: j
   
   this%newFockMatrixPtr => newFockMatrix
   this%newDensityMatrixPtr => newDensityMatrix
   
   if ( present(methodType) ) this%methodType = methodType
   
   if( present(overlapMatrix) .and. ( this%methodType == 2 .or. this%methodType == 4 .or. &
        CONTROL_instance%DIIS_ERROR_IN_DAMPING ) ) then
      
      call Matrix_copyConstructor(this%overlapMatrix, overlapMatrix)
      
      if( .not.allocated(this%orthonormalizationMatrix%values) ) then
         
         !!*****************************************************
         !! Separa espacio de memoria para sistema de ecuaciones DIIS
         !!*****
         call Matrix_constructor( this%diisEquationsSystem, int(this%diisDimensionality+1,8), int(this%diisDimensionality+1,8) )
         this%diisEquationsSystem%values= 0.0_8
         this%diisEquationsSystem%values(1,:) = -1.0_8
         this%diisEquationsSystem%values(:,1) = -1.0_8
         this%diisEquationsSystem%values(1,1)= 0.0_8
         
         orderOfMatrix = size(this%newFockMatrixPtr%values, 1)
         this%orthonormalizationMatrix= this%overlapMatrix
         
         call Vector_constructor( eigenValues, orderOfMatrix )
         call Matrix_constructor( eigenVectors, int(orderOfMatrix,8), int(orderOfMatrix,8) )
         
         call Matrix_eigen( this%orthonormalizationMatrix, eigenValues, eigenVectors, SYMMETRIC  )
         
         do i = 1 , orderOfMatrix
            do j = 1 , orderOfMatrix
               
               if ( eigenValues%values(j) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                  this%orthonormalizationMatrix%values(i,j) = &
                       eigenVectors%values(i,j)/sqrt( eigenValues%values(j) )
               else
                  
                  call Convergence_exception( ERROR, "The matrix of orthonormalization is singular",&
                       "Class object Convergence diis() function" )
                  
               end if
            end do
         end do
         
         !!
         !! Esta matriz oermite obtener un vector de error balanceado e'=A~*e*A al usar una base ortonormal
         !!	donde A es la matriz de ortogonalizacion.
         !!
         this%orthonormalizationMatrix%values =	matmul( eigenVectors%values,&
              matmul( this%orthonormalizationMatrix%values, &
              transpose(eigenVectors%values) ) )
         
         call Matrix_destructor(eigenVectors)
         call Vector_destructor(eigenValues)
         
      end if
      
   end if
   
   if( present(overlapMatrix) .and. present(coefficientMatrix) .and. ( this%methodType == 3 )) then
      
      call Matrix_copyConstructor(this%overlapMatrix, overlapMatrix)
      call Matrix_copyConstructor(this%coefficientMatrix, coefficientMatrix)
      this%speciesID=speciesID  
      
   end if
   
 end subroutine Convergence_setMethod
 
 !>
 !! @brief Indica si el objeto ha sidi instanciado
 function Convergence_isInstanced( this ) result( output )
   implicit none
   type(Convergence), intent(in) :: this
   logical :: output
   
   output =this%isInstanced
   
 end function Convergence_isInstanced
 
 !>
 !! @brief devuelve el nombre del metodo
 function Convergence_getName( this ) result( output )
   implicit none
   type(Convergence), intent(in) :: this
   character(30) :: output
   
   output =trim(this%name)
   
 end function Convergence_getName
 
 !>
 !! @brief devuelve el nombre del metodo
 function Convergence_getDiisError( this ) result( output )
   implicit none
   type(Convergence), intent(in) :: this
   real(8) :: output
   
   output = this%diisError
   
 end function Convergence_getDiisError
 
 !>
 !! @brief devuelve la matrix de densidad inicial
 function Convergence_getInitialDensityMatrix( this ) result( output )
   implicit none
   type(Convergence), intent(in) :: this
   type(Matrix) :: output
   
   if ( allocated( this%initialDensityMatrix%values ) ) then
      
      call Matrix_copyConstructor( output, this%initialDensityMatrix )
      
   else
      
      call Convergence_exception( ERROR, "The object has not beeb instanced",&
           "Class object Convergence in getInitialDensityMatrix() function" )
      
   end if
   
 end function Convergence_getInitialDensityMatrix
 
 !>
 !! @brief devuelve la matrix de Fock inicial
 function Convergence_getInitialFockMatrix( this ) result( output )
   implicit none
   type(Convergence), intent(in) :: this
   type(Matrix) :: output
   
   if ( allocated( this%initialFockMatrix%values ) ) then
      
      call Matrix_copyConstructor( output, this%initialFockMatrix )
      
   else
      
      call Convergence_exception( ERROR, "The object has not beeb instanced",&
           "Class object Convergence in getInitialFockMatrix() function" )
      
   end if
   
 end function Convergence_getInitialFockMatrix
 
 !>
 !! @brief Calcula un mejor estimado de la matriz de Fock, mediante el
 !!		metodo de convergencia seleccionado
 subroutine Convergence_run(this)
   implicit none
   type(Convergence) :: this
   
   if ( this%isInstanced ) then
      
      select case(this%methodType)
         
      case(SCF_CONVERGENCE_NONE)
         
         call Convergence_noneMethod( this)
         
      case(SCF_CONVERGENCE_DAMPING)
         
         call Convergence_dampingMethod( this )
         
      case(SCF_CONVERGENCE_DIIS)
         
         if ( StackMatrices_isInstanced( this%fockMatrices ) .eqv. .false. ) then
            
            call StackMatrices_constructor( this%fockMatrices, &
                 ssize = this%diisDimensionality )
            call StackMatrices_constructor( this%componentsOfError, &
                 ssize = this%diisDimensionality )
         end if
         
         call Convergence_diisMethod( this )
         
      case(SCF_CONVERGENCE_LEVEL_SHIFTING)
         
         call Convergence_dampingMethod( this )
         
         call Convergence_levelShifting( this )
         
      case( SCF_CONVERGENCE_MIXED )
         
         if ( StackMatrices_isInstanced( this%fockMatrices ) .eqv. .false. ) then
            
            call StackMatrices_constructor( this%fockMatrices, &
                 ssize = this%diisDimensionality )
            call StackMatrices_constructor( this%componentsOfError, &
                 ssize = this%diisDimensionality )
         end if
         
         call Convergence_diisMethod( this )
         if ( abs(this%diisError) > CONTROL_instance%DIIS_SWITCH_THRESHOLD .and. .not.this%hadLowDiisError ) then
            call Convergence_dampingMethod( this )
         else
            CONTROL_instance%DIIS_SWITCH_THRESHOLD=0.5
         end if
         
      case default
         
         call Convergence_dampingMethod( this )
         
      end select
      
   else
      
      call Convergence_exception( ERROR, "The object has not beeb instanced",&
           "Class object Convergence in run() function" )
      
   end if
   
 end subroutine Convergence_run
 
 !>
 !! @brief No realiza ninguna modificacion a la matriz de Fock o a la
 !! 		matriz de densidad asociada
 subroutine Convergence_noneMethod(  this )
   implicit none
   type(Convergence) :: this
   
 end subroutine Convergence_noneMethod
 
 
 !>
 !! @brief Calcula una nueva matriz de Fock con el metodo de amortiguamiento
 subroutine Convergence_dampingMethod( this )
   implicit none
   type(Convergence) :: this
   
   real(8) :: densityEffect
   real(8) :: fockAndDensityEffect
   real(8) :: dampingFactor
   type(Matrix) :: auxMatrix
   
   if( CONTROL_instance%DIIS_ERROR_IN_DAMPING ) then
      
      call Convergence_obtainDiisError(this)
   end if
   !!********************************************************************************************
   !! Evalua el efecto del cambio de la matriz densidad en la energia
   !!
   
   call Matrix_copyConstructor( auxMatrix, this%initialFockMatrix)
   auxMatrix%values = matmul( auxMatrix%values, (this%newDensityMatrixPtr%values &
        - this%initialDensityMatrix%values) )
   
   densityEffect = -0.5_8 * Matrix_trace( auxMatrix )
   
   !!
   !!********************************************************************************************
   
   !!********************************************************************************************
   !! Evalua el efecto de los cambios de la matriz densidad y  la matiz de fock en energia
   !!
   auxMatrix%values = matmul( (this%newFockMatrixPtr%values &
        - this%initialFockMatrix%values), &
        ( this%newDensityMatrixPtr%values &
        - this%initialDensityMatrix%values) )
   
   fockAndDensityEffect = Matrix_trace( auxMatrix )
   
   !!
   !!********************************************************************************************
   
   if ( fockAndDensityEffect <=  densityEffect ) then
      
      this%initialFockMatrix%values = this%newFockMatrixPtr%values
      this%initialDensityMatrix%values = this%newDensityMatrixPtr%values
      
   else
      dampingFactor = densityEffect / fockAndDensityEffect
      this%initialFockMatrix%values = this%initialFockMatrix%values &
           + dampingFactor * ( this%newFockMatrixPtr%values &
           - this%initialFockMatrix%values )
      
      !! Modifica la matrix de Fock con una matriz amortiguada
      !! Modificacion del damping factor para prueba
      !                                         dampingFactor=0.1_8
      this%newFockMatrixPtr%values = this%initialFockMatrix%values

      this%initialDensityMatrix%values = this%initialDensityMatrix%values &
           + dampingFactor * ( this%newDensityMatrixPtr%values &
           - this%initialDensityMatrix%values )
      
   end if
   
   call Matrix_destructor( auxMatrix )
   
 end subroutine Convergence_dampingMethod
 
 !>
 !! @brief Calcula una nueva matriz de densidad con la componentes de un
 !!		subespacio definido con la matrices de Fock anteriores
 !! 		y a la matriz de densidad asociada
 subroutine Convergence_diisMethod( this)
   implicit none
   type(Convergence) :: this
   
   type(Matrix) :: currentError
   real(8), pointer :: currentErrorPtr(:,:)
   real(8), pointer :: beforeErrorPtr(:,:)
   real(8), pointer :: fockMatrix(:,:)
   real(8), allocatable :: rightSideCoefficients(:)
   real(8), allocatable :: coefficientsOfCombination(:)
   type(Matrix) :: diisEquationsSystem !! Debe almacenarse una para cada especie
   type(Matrix) :: rightMatrix
   type(Matrix) :: leftMatrix
   type(Matrix) :: singularValues
   type(Matrix) :: auxMatrix
   integer :: orderOfMatrix
   integer :: i
   integer :: j
   
   orderOfMatrix = size(this%newFockMatrixPtr%values, 1)
   
   call Matrix_constructor( currentError, int(orderOfMatrix,8), int(orderOfMatrix,8 ) )
   
   !! Calcula FDS
   currentError%values = matmul( matmul( this%newFockMatrixPtr%values, this%newDensityMatrixPtr%values), this%overlapMatrix%values )
   
   !! Calcula la matriz de error FDS-SDF
   currentError%values = currentError%values - transpose(currentError%values)
   
   !! Ortonormaliza la componente de error
   currentError%values = matmul(transpose(this%orthonormalizationMatrix%values), &
        matmul( currentError%values, this%orthonormalizationMatrix%values) )
   
   !! Determina el error para la nueva componente
   this%diisError = abs( maxval(currentError%values) )
   
   if( this%diisError < CONTROL_instance%DIIS_SWITCH_THRESHOLD )  then
      
      !! Almacena matriz de error y matriz de Fock
      call StackMatrices_push( this%componentsOfError, currentError%values )
      call StackMatrices_push( this%fockMatrices,  this%newFockMatrixPtr%values)
      
      this%currentSize = StackMatrices_size( this%fockMatrices )
      
      if ( this%currentsize >= 2 ) then
         
         !! Vector con valores de la constantes  del sistema lineal
         allocate( rightSideCoefficients( this%currentSize + 1 ) )
         rightSideCoefficients=0.0_8
         rightSideCoefficients(1)= -1.0_8
         
         allocate( coefficientsOfCombination( this%currentSize + 1 ) )
         call Matrix_constructor(auxMatrix, int(this%currentsize+1,8), int(this%currentsize+1,8) )
         
         !! Elimina la primera fila de la matriz de error y mueve las siguiente en un nivel
         if ( this%isLastComponent ) then
            
            this%diisEquationsSystem%values(2:this%currentsize,2:this%currentSize) = &
                 this%diisEquationsSystem%values(3:this%currentSize+1,3:this%currentSize+1)
            
         end if
         
         if ( this%currentSize <= this%diisDimensionality ) then
            
            !!
            !! Calcula la ultima fila de la matriz de error
            !!
            currentErrorPtr => StackMatrices_getValuesPtr( this%componentsOfError, this%currentSize )
            do i=this%currentSize,1,-1
               beforeErrorPtr  => StackMatrices_getValuesPtr( this%componentsOfError, i )
               this%diisEquationsSystem%values(i+1,this%currentSize+1) = sum(currentErrorPtr*beforeErrorPtr)
               this%diisEquationsSystem%values(this%currentSize+1,i+1)  =  &
                    this%diisEquationsSystem%values(i+1,this%currentSize+1)
               beforeErrorPtr => null()
            end do
            currentErrorPtr => null()
            if( this%currentSize==this%diisDimensionality ) this%isLastComponent=.true.
         end if
         
         auxMatrix%values = this%diisEquationsSystem%values(1:this%currentSize+1,1:this%currentSize+1)
         
         !! Diagonaliza la matiz por SVD de la forma A=UWV^T
         call Matrix_svd( auxMatrix, leftMatrix, rightMatrix, singularValues)
         
         !! Invierte la matriz de valores singulares
         do i=1,size(singularValues%values,dim=1)
            if( abs(singularValues%values(i,i)) > 1.0D-4 ) then
               singularValues%values(i,i) = 1/singularValues%values(i,i)
            else
               singularValues%values(i,i) = singularValues%values(i,i)
            end if
         end do
         
         coefficientsOfCombination=matmul( &
              matmul(matmul(transpose(rightMatrix%values), &
              singularValues%values),transpose(leftMatrix%values)),rightSideCoefficients )
         
         call Matrix_destructor( auxMatrix )
         call Matrix_destructor( singularValues )
         call Matrix_destructor( rightMatrix )
         call Matrix_destructor( leftMatrix )
         
         !! Construye una nueva matriz de Fock con una combinacion lineal de elementos del sub-espacio de F
         this%newFockMatrixPtr%values = 0.0_8
         
         do i=1, this%currentSize
            fockMatrix => StackMatrices_getValuesPtr( this%fockMatrices, i )
            this%newFockMatrixPtr%values =  &
                 this%newFockMatrixPtr%values + coefficientsOfCombination(i+1) * fockMatrix
            fockMatrix => null()
         end do
         
         deallocate( rightSideCoefficients )
         deallocate( coefficientsOfCombination )
         
      else &
           if ( this%currentsize == 1 ) then
         
         currentErrorPtr => StackMatrices_getValuesPtr( this%componentsOfError, 1 )
         beforeErrorPtr  => StackMatrices_getValuesPtr( this%componentsOfError, 1 )
         this%diisEquationsSystem%values(2,2) = sum(currentErrorPtr*beforeErrorPtr)
         currentErrorPtr => null()
         beforeErrorPtr => null()
         
      end if
      this%hadLowDiisError = .true.
   end if
   
   call Matrix_destructor(currentError)
   
 end subroutine Convergence_diisMethod
 
 !>
 !! @brief Retorna el error diis de una matriz de fock
 subroutine Convergence_obtainDiisError(this)
   implicit none
   type(Convergence) :: this
   
   type(Matrix) :: currentError
   integer :: orderOfMatrix
   
   orderOfMatrix = size(this%newFockMatrixPtr%values, 1)
   call Matrix_constructor( currentError, int(orderOfMatrix,8), int(orderOfMatrix,8 ) )
   
   !! Calcula FDS
   currentError%values = matmul( matmul( this%newFockMatrixPtr%values,&
        this%newDensityMatrixPtr%values), this%overlapMatrix%values )
   
   !! Calcula la matriz de error FDS-SDF
   currentError%values = currentError%values - transpose(currentError%values)
   
   !! Ortonormaliza la componente de error
   currentError%values = matmul(transpose(this%orthonormalizationMatrix%values), &
        matmul( currentError%values, this%orthonormalizationMatrix%values) )
   
   this%diisError = abs( maxval(currentError%values) )
   
   call Matrix_destructor( currentError)
   
 end subroutine Convergence_obtainDiisError
 
 !>
 !! @brief Modifica la Matriz de Fock empleando el metodo level shifting
 !!
 !! @todo Falta la implementacion completa
 subroutine Convergence_levelShifting(this)
   implicit none
   type(Convergence) :: this
   
   CONTROL_instance%ACTIVATE_LEVEL_SHIFTING=.true.
   
   ! real(8) :: factor
   ! type(Matrix) :: auxMatrix
   ! type(Matrix) :: auxOverlapMatrix
   ! integer :: occupiedOrbitals
   ! integer :: totalOrbitals
   ! integer :: i
   
   ! if ( (this%name).ne."e-" ) then
   !     factor=-0.1_8
   !     else
   !     factor=1.0_8
   ! end if
   
   ! call Matrix_copyConstructor( auxMatrix, this%initialFockMatrix)
   
   ! auxMatrix%values = &
   ! 	matmul( matmul( transpose( this%coefficientMatrix%values ) , &
   ! 		auxMatrix%values), this%coefficientMatrix%values )
   
   ! totalOrbitals =Matrix_getNumberOfRows(auxMatrix)
   ! occupiedOrbitals=MolecularSystem_getOcupationNumber( this%speciesID )
   
   ! do i= occupiedOrbitals + 1, totalOrbitals
   !   auxMatrix%values(i,i) = factor*CONTROL_instance%LEVEL_SHIFTING + auxMatrix%values(i,i)
   ! end do
   
   ! auxMatrix%values = &
   ! 	matmul( matmul(this%coefficientMatrix%values, &
   ! 		auxMatrix%values), transpose(this%coefficientMatrix%values ))
   
   ! call Matrix_copyConstructor(auxOverlapMatrix,this%overlapMatrix)
   
   ! auxMatrix%values = &
   ! 	matmul( matmul(auxOverlapMatrix%values, &
   ! 		auxMatrix%values), transpose(auxOverlapMatrix%values ))
   
   ! 	!! Changing the fock matrix
   
   ! 		this%newFockMatrixPtr%values = auxMatrix%values
   
 end subroutine Convergence_levelShifting
 
 
 subroutine Convergence_reset()
   implicit none
   
   CONTROL_instance%DIIS_SWITCH_THRESHOLD = CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP
   
 end subroutine Convergence_reset
 
 !>
 !! @brief  Maneja excepciones de la clase
 subroutine Convergence_exception( typeMessage, description, debugDescription)
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
   
 end subroutine Convergence_exception
 
end module Convergence_
