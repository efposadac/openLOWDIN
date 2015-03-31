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
!! @brief Module to minimize using Trust Region method
!! @author  J.M. Rodas
!! @author  S.A. Gonzalez
!!
!! <b> Creation date : </b> 2009-06-11
!!
!! <b> History: </b>
!!
!!   - <tt> 2009-06-11 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Basics functions and functions has been created using trust region method
!!   - <tt> 2015-02-25 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0 and prepare the module for new optimizers
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module TrustRegionOptimizer_
  use CONTROL_
  use Matrix_
  use Vector_
  use Math_
  use Exception_
  implicit none

  type, public :: TrustRegionOptimizer
     character(20) :: name
     ! integer :: outputUnid
     integer :: numberOfVariables
     integer :: realNumberOfIteration
     integer :: numberOfIteration
     integer :: oldNumberOfIteration
     integer :: principalMode
     integer :: maximumNumberOfIterations
     integer :: maximizationMode
     type(Vector)  :: variables
     type(Vector)  :: gradient
     type(Vector)  :: oldVariables
     type(Vector)  :: oldGradient
     type(Vector)  :: step
     type(Vector) :: vibrationalMode
     type(Vector)  :: gradientProjectedOnExtDegrees
     type(Vector)  :: gradientProjectedOnHessiane
     type(Matrix) :: hessiane
     type(Matrix) :: hessianeProjected
     real(8) :: functionValue
     real(8) :: oldFunctionValue
     real(8) :: trustRadius
     real(8) :: lagrangeMultiplier
     real(8) :: oldLagrangeMultiplier
     real(8) :: xlambda
     real(8) :: oldXlambda
     real(8) :: squareStepSize
     real(8) :: minimumStepSize
     real(8) :: maximumStepSize
     real(8) :: minimumTrustRadius
     real(8) :: maximumTrustRadius
     real(8) :: optimizationTolerance
     real(8) :: predictedChangeOfFunction
     real(8) :: realChangeOfFunction
     real(8) :: adaptativeMaximumTrustRadius
     real(8) :: ratioOfChange
     real(8) :: overlapInterModes
     logical :: hasConverged
     logical :: isHardCase
     logical :: isSuitableStep
     logical :: hasRestarted
  end type TrustRegionOptimizer

  !< enum Vector_printFormatFlags {
  real(8), parameter :: MAXIMUM_TRUST_RADIO = 0.3_8
  real(8), private, parameter :: MINIMUM_TRUST_RADIO = 0.001_8
  real(8), private, parameter :: ZERO_BOUND = 1.0D-12
  real(8), private, parameter :: epsilonTol=1.0D-16
  !< }


  public :: &
       TrustRegionOptimizer_constructor, &
       TrustRegionOptimizer_destructor, &
       ! TrustRegionOptimizer_getHessiane,&
       TrustRegionOptimizer_getGradient, &
       ! TrustRegionOptimizer_getNumberOfIterations,&
       ! TrustRegionOptimizer_getRealNumberOfIterations,&
       TrustRegionOptimizer_setInitialHessian, &
       ! TrustRegionOptimizer_setMaximumIterations, &
       TrustRegionOptimizer_iterate, &
       TrustRegionOptimizer_isMinimum!, &
       ! TrustRegionOptimizer_restart


  private
contains

  subroutine TrustRegionOptimizer_constructor( this, initialPoint )
    implicit none
    type(TrustRegionOptimizer) :: this
    real(8) :: initialPoint(:)

    this%numberOfVariables = size( initialPoint )
    call Vector_constructor( this%variables, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%oldVariables, this%numberOfVariables, 0.0_8 )
    this%variables%values = initialPoint
    call Vector_constructor( this%gradient, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%oldGradient, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%step, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%gradientProjectedOnExtDegrees, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%vibrationalMode, this%numberOfVariables, 0.0_8 )
    call Vector_constructor( this%gradientProjectedOnHessiane, this%numberOfVariables, 0.0_8 )
    call Matrix_constructor( this%hessiane, int(this%numberOfVariables,8), int(this%numberOfVariables,8), 0.0_8 )
    call Matrix_constructor( this%hessianeProjected, int(this%numberOfVariables,8), int(this%numberOfVariables,8), 0.0_8 )
    this%trustRadius = 0.0_8
    this%optimizationTolerance = CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT
    this%realNumberOfIteration = 0
    this%numberOfIteration = 0
    this%principalMode = 1
    this%maximumNumberOfIterations = CONTROL_instance%MINIMIZATION_MAX_ITERATION
    this%hasConverged = .false.
    this%isSuitableStep = .false.
    this%hasRestarted = .false.
    ! * = CONTROL_instance%UNID_FOR_OUTPUT_FILE
    this%adaptativeMaximumTrustRadius= MAXIMUM_TRUST_RADIO
    this%minimumTrustRadius = 0.05_8
    this%maximumTrustRadius = 0.5_8
  end subroutine TrustRegionOptimizer_constructor

  subroutine TrustRegionOptimizer_destructor( this)
    implicit none
    type(TrustRegionOptimizer) :: this

    call Vector_destructor( this%variables )
    call Vector_destructor( this%oldVariables )
    call Vector_destructor( this%gradient )
    call Vector_destructor( this%oldGradient )
    call Vector_destructor( this%step )
    call Vector_destructor( this%gradientProjectedOnExtDegrees )
    call Vector_destructor( this%gradientProjectedOnHessiane )
    call Matrix_destructor( this%hessiane)
    call Matrix_destructor( this%hessianeProjected)
    this%trustRadius = 0.1_8
    this%realNumberOfIteration = 0
    this% numberOfIteration = 0
    this%numberOfVariables = 0
    this%hasConverged = .false.
    this%isSuitableStep = .false.
    this%hasRestarted = .false.

  end subroutine TrustRegionOptimizer_destructor

!   function TrustRegionOptimizer_getHessiane( this ) result(output)
!     implicit none
!     type(TrustRegionOptimizer) :: this
!     type(Matrix) :: output

!     output = this%hessiane

!   end function TrustRegionOptimizer_getHessiane

  function TrustRegionOptimizer_getGradient( this ) result(output)
    implicit none
    type(TrustRegionOptimizer) :: this
    real(8) :: output

    output = sqrt(dot_product(this%gradient%values,this%gradient%values)/this%numberOfVariables)

  end function TrustRegionOptimizer_getGradient



!   !>
!   !! @brief Retorna el numero de iteraciones realizadas hasta el momento
!   !<
!   function TrustRegionOptimizer_getNumberOfIterations( this ) result(output)
!     implicit none
!     type(TrustRegionOptimizer) :: this
!     integer :: output

!     output = this%numberOfIteration

!   end function TrustRegionOptimizer_getNumberOfIterations

!   !>
!   !! @brief Retorna el numero de iteraciones realizadas hasta el momento
!   !<
!   function TrustRegionOptimizer_getRealNumberOfIterations( this ) result(output)
!     implicit none
!     type(TrustRegionOptimizer) :: this
!     integer :: output

!     output = this%realNumberOfIteration

!   end function TrustRegionOptimizer_getRealNumberOfIterations

  function TrustRegionOptimizer_isMinimum( this ) result(output)
    implicit none
    type(TrustRegionOptimizer) :: this
    logical :: output

    call TrustRegionOptimizer_checkOptimizationCriteria( this )
    output = this%hasConverged

  end function TrustRegionOptimizer_isMinimum



  subroutine TrustRegionOptimizer_setInitialHessian( this, initialHessianMatrix)
    implicit none
    type(TrustRegionOptimizer) :: this
    real(8) :: initialHessianMatrix(:,:)

    this%hessiane%values = initialHessianMatrix

  end subroutine TrustRegionOptimizer_setInitialHessian

!   subroutine TrustRegionOptimizer_setMaximumIterations( this, numberOfIterations)
!     implicit none
!     type(TrustRegionOptimizer) :: this
!     integer :: numberOfIterations

!     this%maximumNumberOfIterations = numberOfIterations

!   end subroutine TrustRegionOptimizer_setMaximumIterations


  subroutine TrustRegionOptimizer_iterate( this, functionValue, gradientOfFunction, &
       projectGradientFunction, &
       projectHessianeFunction, &
       showFunction )
    implicit none
    type(TrustRegionOptimizer) :: this

    integer :: i
    interface

       function functionValue( pointOfEvaluation ) result( output )
         real(8) :: pointOfEvaluation(:)
         real(8) :: output
       end function functionValue

       subroutine gradientOfFunction( pointOfEvaluation, gradient )
         real(8) :: pointOfEvaluation(:)
         real(8) :: gradient(:)
       end subroutine gradientOfFunction

       subroutine projectGradientFunction( gradient, gradientProjectedOnExtDegrees )
         real(8) :: gradient(:)
         real(8) :: gradientProjectedOnExtDegrees(:)
       end subroutine projectGradientFunction

       subroutine projectHessianeFunction( hessiane,hessianeProjected )
         real(8) :: hessiane(:,:)
         real(8) :: hessianeProjected(:,:)
       end subroutine projectHessianeFunction

       subroutine showFunction( pointOfEvaluation, gradient, functionValue )
         real(8) :: pointOfEvaluation(:)
         real(8) :: gradient(:)
         real(8) :: functionValue
       end subroutine showFunction

    end interface

    if ( this%realNumberOfIteration == 0  ) then

       !! Acota el radio de confianza inicial
       this%trustRadius = 0.3
       !! Calcula el valor de la funcion y su gradiente en el punto inicial
       this%oldFunctionValue=0.0
       this%functionValue= functionValue( this%variables%values )

       call gradientOfFunction(this%variables%values, this%gradient%values )
       this%hasRestarted=.false.

       !! Muestra el estimado actual de la minimizacion
       write(*,"(T10,A21,I5)") "NUMBER OF ITERATION: ", 0
       write(*,"(T10,A20)") "-----------------------------------------------------"
       write(*,"(A1)") " "
       call showFunction( this%variables%values, this%gradient%values, this%functionValue )
       write(*,"(T10,A19,F12.7)") "RMS GRADIENT     = ", TrustRegionOptimizer_getGradient(this)
       write(*,*) ""

    end if

    if ( this%realNumberOfIteration == 0 .or. this%hasRestarted ) this%trustRadius = min(MAXIMUM_TRUST_RADIO,this%trustRadius)

    this%realNumberOfIteration = this%realNumberOfIteration + 1
    this%numberOfIteration = this%numberOfIteration +1

    write (*,"(T10,A21,I5)") "NUMBER OF ITERATION: ", this%realNumberOfIteration
    write(*,"(T10,A20)") "-----------------------------------------------------------"
    write(*,"(A1)") " "

    if ( this%realNumberOfIteration <= this%maximumNumberOfIterations ) then

       if( this%realNumberOfIteration > 1) &
            call TrustRegionOptimizer_updateHessiane( this )

       !! Proyecta los grados de libetad externos del gradiente
       call projectGradientFunction( this%gradient%values, this%gradientProjectedOnExtDegrees%values )

       ! !! Proyecta los grados de libetad externos de la hesiana
       call projectHessianeFunction( this%hessiane%values, this%hessianeProjected%values )

       ! !! Calcula la nueva direccion de busqueda
       call TrustRegionOptimizer_calculateStep(this)
       call TrustRegionOptimizer_predictedChangeOfFunction( this )
       call TrustRegionOptimizer_updateTrustRadio( this )

       !           if (this%isSuitableStep) then
       !! Almacena valores actuales del punto sobre la superficie
       this%oldVariables%values = this%variables%values
       this%oldGradient%values = this%gradient%values
       this%oldFunctionValue   = this%functionValue

       !! Ejecuta el paso
       this%variables%values=this%variables%values + this%step%values

       !! Calcula el valor de la funcion y el gradiente
       this%functionValue = functionValue( this%variables%values )
       call gradientOfFunction( this%variables%values, this%gradient%values )

       !! Muestra el estimado actual de la minimizacion
       call showFunction( this%variables%values, this%gradient%values, this%functionValue )
       !           else
       !               print *,"The current no is valid"
       !           end if

    else

       call TrustRegionOptimizer_exception( ERROR, "The maximum number of iterations was exceded", &
            "Class object TrustRegionOptimizer in iterate() function")

    end if

  end subroutine TrustRegionOptimizer_iterate

  subroutine TrustRegionOptimizer_checkOptimizationCriteria(this)
    implicit none
    type(TrustRegionOptimizer) :: this

    real(8) :: rmsOfGradient
    real(8) :: rmsOfDisplacement
    real(8) :: maximumGradient
    real(8) :: maximumDisplacement
    real(8) :: thresholdForRmsGradient
    real(8) :: thresholdForRmsDisplacement
    real(8) :: thresholdForSingleGradient
    real(8) :: thresholdForSingleDisplacement
    real(8) :: thresholdForFunctionChange
    real(8) :: changeOfFunction
    logical :: isRmsOfGradientConverged
    logical :: isRmsOfDisplacementConverged
    logical :: isMaximumGradientConverged
    logical :: isMaximumDisplacementConverged
    logical :: isFunctionChangeConverged
    character(4) :: stateOfConvergenceOfRmsGradient
    character(4) :: stateOfConvergenceOfRmsDisplacement
    character(4) :: stateOfConvergenceOfSingleGradient
    character(4) :: stateOfConvergenceOfSingleDisplacement
    character(4) :: stateOfConvergenceOfFunctionChange

    !!******************************************************
    !! define los umbrales de convergencia para gradiente y desplazamiento
    thresholdForSingleGradient = 5.0*this%optimizationTolerance
    thresholdForRmsGradient = 1.0*this%optimizationTolerance
    thresholdForSingleDisplacement = 6.0*this%optimizationTolerance
    thresholdForRmsDisplacement = 4.0*this%optimizationTolerance
    thresholdForFunctionChange = 1.0D-5
    !!
    !!******************************************************


    !!*******************************************************
    !! Determina convergencia en gradiente RMS
    !!****
    isRmsOfGradientConverged=.false.
    stateOfConvergenceOfRmsGradient = "--"
    rmsOfGradient = sqrt(dot_product(this%gradient%values,this%gradient%values)/this%numberOfVariables)
    if (rmsOfGradient <= thresholdForRmsGradient) then
       stateOfConvergenceOfRmsGradient = "OK  "
       isRmsOfGradientConverged=.true.
    end if
    !!
    !!*******************************************************

    !!*******************************************************
    !! Determina convergencia en desplazamiento RMS
    !!****
    isRmsOfDisplacementConverged =.false.
    stateOfConvergenceOfRmsDisplacement = "--  "
    rmsOfDisplacement = sqrt(dot_product(this%step%values,this%step%values)/this%numberOfVariables)
    if (rmsOfDisplacement <= thresholdForRmsDisplacement) then
       stateOfConvergenceOfRmsDisplacement = "OK  "
       isRmsOfDisplacementConverged =.true.
    end if
    !!
    !!*******************************************************

    !!*******************************************************
    !! Determina convergencia para la componente mas grande del gradiente
    !!****
    isMaximumGradientConverged=.false.
    stateOfConvergenceOfSingleGradient = "--  "
    maximumGradient = max(maxval(this%gradient%values),abs(minval(this%gradient%values)))
    if ( maximumGradient <= thresholdForSingleGradient ) then
       stateOfConvergenceOfSingleGradient = "OK  "
       isMaximumGradientConverged=.true.
    end if
    !!
    !!*******************************************************

    !!*******************************************************
    !! Determina convergencia para la componente mas grande del desplazamiento
    !!****
    isMaximumDisplacementConverged =.false.
    stateOfConvergenceOfSingleDisplacement="--  "
    maximumDisplacement = max(maxval(this%step%values),abs(minval(this%step%values)))

    if (maximumDisplacement <= thresholdForSingleDisplacement) then
       stateOfConvergenceOfSingleDisplacement="OK  "
       isMaximumDisplacementConverged =.true.
    end if
    !!
    !!*******************************************************


    !!*******************************************************
    !!  Determina el cambio de la funcion entre las dos ultimas iteraciones
    !!****
    changeOfFunction = abs(this%functionValue - this%oldFunctionValue)
    changeOfFunction = min(abs(this%oldFunctionValue),changeOfFunction)
    isFunctionChangeConverged = .false.
    stateOfConvergenceOfFunctionChange="--  "
    if ( changeOfFunction  < thresholdForFunctionChange ) then
       stateOfConvergenceOfFunctionChange="OK  "
       isFunctionChangeConverged =.true.
    end if
    !!
    !!*******************************************************

    write(*,"(T10,A19,F12.7,A13,F12.7,A1,A4)") "RMS GRADIENT     = ",rmsOfGradient, &
         " THRESHOLD = ",thresholdForRmsGradient," ",trim(stateOfConvergenceOfRmsGradient)
    write(*,"(T10,A19,F12.7,A13,F12.7,A1,A4)") "MAX GRADIENT     = ",maximumGradient, &
         " THRESHOLD = ",thresholdForSingleGradient," ",trim(stateOfConvergenceOfSingleGradient)
    write(*,"(T10,A19,F12.7,A13,F12.7,A1,A4)") "RMS DISPLACEMENT = ",rmsOfDisplacement,&
         " THRESHOLD = ",thresholdForRmsDisplacement," ",trim(stateOfConvergenceOfRmsDisplacement)
    write(*,"(T10,A19,F12.7,A13,F12.7,A1,A4)") "MAX DISPLACEMENT = ",maximumDisplacement,&
         " THRESHOLD = ",thresholdForSingleDisplacement," ",trim(stateOfConvergenceOfSingleDisplacement)
    write(*,"(T10,A19,F12.7,A13,F12.7,A1,A4)") "FUNCTION  CHANGE = ",changeOfFunction,&
         " THRESHOLD = ",thresholdForFunctionChange," ",trim(stateOfConvergenceOfFunctionChange)
    write(*,*) ""

    if(maximumGradient> 1000) then

       call TrustRegionOptimizer_exception( ERROR, "The maximum gradient is very high", &
            "Class object TrustRegionOptimizer in checkOptimizationCriteria() function")

    end if

    !!*********************************************************
    !! Analiza los criterios de convergencia
    !!****
    if ( isMaximumGradientConverged .and. isRmsOfGradientConverged ) then


       if ( isRmsOfDisplacementConverged .and. isMaximumDisplacementConverged ) then
          this%hasConverged = .true.
          write (*,"(T10,A57)") "THE THRESHOLD FOR GRADIENT AND DISPLACEMENT WERE ACHIEVED"

       else if ( isFunctionChangeConverged ) then
          this%hasConverged = .true.
          write (*,"(T10,A63)") "THE THRESHOLD FOR GRADIENT AND CHANGE OF FUNCTION WERE ACHIEVED"

       else if (this%realNumberOfIteration == 1) then
          this%hasConverged = .true.
          write (*,"(T15,A32)") "THE CURRET POSITION IS A MINIMUM"
       end if

    else
       this%hasConverged = .false.
       write (*,"(T2,A25)") "---BEGINNING NEW SEARCH---"
    end if
    !!
    !!*********************************************************


    if ( this%hasConverged ) then
       write(*,*) ""
       write(*,"(T15,A35)") "//////////////////////////////////////////////////////////////////////////////////////////////////////////////"
       write(*,"(T15,A23,I5,A7)") "MINIMUM LOCATED AFTER ",this%realNumberOfIteration," STEPS"
       write(*,"(T15,A35)") "//////////////////////////////////////////////////////////////////////////////////////////////////////////////"
       write(*,*) ""
    end if

  end subroutine TrustRegionOptimizer_checkOptimizationCriteria


  subroutine TrustRegionOptimizer_calculateStep( this)
    implicit none
    type(TrustRegionOptimizer) :: this

    real(8) :: tau
    integer :: i
    integer :: j
    integer :: numberOfTransAndRot
    integer :: lowest
    integer :: negativeLowest
    integer :: numberOfSearches
    real(8), allocatable :: auxMatrix(:,:)
    real(8), allocatable :: auxVector(:)
    real(8) :: gradientThreshold
    real(8) :: auxVal
    real(8) :: auxEigenVal
    real(8) :: step
    real(8) :: FL,FU,FM,BU,BL,BULI
    logical :: flagA,flagB,flagC


    allocate( auxMatrix(this%numberOfVariables,this%numberOfVariables) )
    allocate( auxVector(this%numberOfVariables) )
    call Matrix_symmetrize( this%hessianeProjected, 'L' )

    !! Calcula vectores y valores propios, acotando los ultimos entre el rango (1E-8,1E3)
    call Matrix_eigenProperties(this%hessianeProjected)

    !!********************************************************************************
    !! Ordena los vectores propios de manera que los grados translacionales y rotaciones
    !! queden al final
    !!***********************

    j=0
    numberOfTransAndRot=0
    auxMatrix=0.0
    auxVector=0.0
    do i=1, this%numberOfVariables
       if( abs(this%hessianeProjected%eigenValues(i) ) <= CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
          numberOfTransAndRot=numberOfTransAndRot+1

          if(numberOfTransAndRot>6) then
             call TrustRegionOptimizer_exception(WARNING,"There are more of six null eigenValues  ", &
                  "Class object TrustRegionOptimizer in calculateStep() function")
          end if

          auxMatrix(:,i)=this%hessianeProjected%eigenVectors(:,i)
          auxVector(i)=this%hessianeProjected%eigenValues(i)
       else
          j=j+1
          this%hessianeProjected%eigenVectors(:,j)=this%hessianeProjected%eigenVectors(:,i)
          this%hessianeProjected%eigenValues(j)=this%hessianeProjected%eigenValues(i)
       end if
    end do
    this%hessianeProjected%eigenVectors(:,j+1:this%numberOfVariables)=auxMatrix(:,1:numberOfTransAndRot)
    this%hessianeProjected%eigenValues(j+1:this%numberOfVariables)=auxVector(1:numberOfTransAndRot)

    !! Determina si la hesiana es definida positiva, acotando a los valores propios a 1E-6
    this%hessianeProjected%isPositiveDefinited=Matrix_isPositiveDefinited( this%hessianeProjected, 1.0D-6 )

    deallocate(auxMatrix)
    deallocate(auxVector)
    !
    !!********************************************************************************

    !! Proyecta el gradiente sobre los grados de libertad externos
    this%gradientProjectedOnHessiane%values= &
         matmul( transpose(this%hessianeProjected%eigenVectors), &
         this%gradientProjectedOnExtDegrees%values )

    do i=1,this%numberOfVariables
       if( abs(this%hessianeProjected%eigenValues(i))< 1.0D-8 )&
            this%gradientProjectedOnHessiane%values(i)=0.0_8
    end do

    if (this%numberOfIteration > 1) then
       this%realChangeOfFunction= this%functionValue-this%oldFunctionValue
       this%ratioOfChange = this%realChangeOfFunction/this%predictedChangeOfFunction

       this%adaptativeMaximumTrustRadius=this%trustRadius
       if( this%ratioOfChange <= 0.1 .or. this%ratioOfChange >= 3.0) &
            this%adaptativeMaximumTrustRadius=this%adaptativeMaximumTrustRadius/2.0
       if( this%ratioOfChange >= 0.75 .and. this%ratioOfChange <= (4.0/3.0) ) &
            this%adaptativeMaximumTrustRadius=this%adaptativeMaximumTrustRadius *sqrt(2.0)

       if ( abs(this%ratioOfChange-1.0) < 0.1 ) &
            this%adaptativeMaximumTrustRadius=this%adaptativeMaximumTrustRadius * sqrt(2.0)
       this%adaptativeMaximumTrustRadius = max(this%adaptativeMaximumTrustRadius, this%minimumTrustRadius)
       this%adaptativeMaximumTrustRadius = min(this%adaptativeMaximumTrustRadius,this%maximumTrustRadius)

    end if

    write(*, "(T10,A25,F12.8)") "ACTUAL FUNCTION CHANGE = ",this%realChangeOfFunction
    write(*, "(T10,A28,F12.8)") "PREDICTED FUNCTION CHANGE = ",this%predictedChangeOfFunction
    write(*, "(T10,A18,F12.8)") "RATIO OF CHANGE = ", this%ratioOfChange


    !!***************************************************************
    !! Inicia proceso de calculo de paso
    !!***************
    negativeLowest = 0
    auxEigenVal=0.0
    flagC = .false.
    if ( .not.this%hessianeProjected%isPositiveDefinited ) then

       if( this%numberOfIteration == 1 ) then
          negativeLowest = this%principalMode
       else

          !! Calcula el acoplamiento entre modos
          negativeLowest = 1
          this%overlapInterModes=dot_product(this%vibrationalMode%values, &
               this%hessianeProjected%eigenVectors(:,1) )
          do i=2,this%numberOfVariables
             auxVal=abs( dot_product(this%vibrationalMode%values, &
                  this%hessianeProjected%eigenVectors(:,i) ) )
             if( auxVal > this%overlapInterModes) then
                this%overlapInterModes = auxVal
                negativeLowest = i
             end if
          end do

       end if

       do i=1, this%numberOfVariables
          this%vibrationalMode%values(i) = &
               this%hessianeProjected%eigenVectors(i,negativeLowest)
       end do
       this%maximizationMode = negativeLowest

       !! Reajusta el radio de confianza cuando el acoplamiento entre modos
       !! es muy pequenho
       if( this%numberOfIteration > 1 ) then
          if( this%overlapInterModes < 0.9) this%adaptativeMaximumTrustRadius= &
               this%trustRadius/sqrt(2.0)
          if( this%overlapInterModes < 0.8) this%adaptativeMaximumTrustRadius= &
               this%adaptativeMaximumTrustRadius/sqrt(2.0)
          if( this%overlapInterModes < 0.7) this%adaptativeMaximumTrustRadius= &
               this%adaptativeMaximumTrustRadius/sqrt(2.0)
          this%adaptativeMaximumTrustRadius=max(this%adaptativeMaximumTrustRadius,this%minimumTrustRadius)
          this%adaptativeMaximumTrustRadius=min(this%adaptativeMaximumTrustRadius,this%maximumTrustRadius)
       end if

       this%principalMode= this%maximizationMode
       negativeLowest = this%principalMode
       auxEigenVal=this%hessianeProjected%eigenValues(negativeLowest)
       if( abs(auxEigenVal) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
          call TrustRegionOptimizer_exception(ERROR,"The current mode will be a translational or a rotational mode, try a new geometry","")
       end if
    end if

    !!***************************************************************************
    !! Busca y verifica que el modo mas bajo el cual no puede ser ni una rotacion ni
    !! una translacion, tenga un gradiente, con el cual se pueda obtener un
    !! multiplicador de lagrange adecuado
    !!*****************
    lowest=1

    lowest_search : do

       do i=lowest,this%numberOfVariables

          if ( this%hessianeProjected%eigenValues(i) ==Math_NaN ) then

             call TrustRegionOptimizer_exception(ERROR,"The hessiane is wrong", &
                  "Class object TrustRegionOptimizer in calculateStep() function")

          else if ( i /= negativeLowest .and. abs(this%hessianeProjected%eigenValues(i)) >CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
             lowest = i
             exit
          end if

       end do

       !! Acota el tamahno del paso
       this%minimumStepSize = max( abs(this%hessianeProjected%eigenValues(lowest))*epsilonTol,10*epsilonTol )
       this%maximumStepSize=1000.0*max(1000.0,abs(this%hessianeProjected%eigenValues(lowest)))

       gradientThreshold=sqrt( this%minimumStepSize * abs(this%hessianeProjected%eigenValues(lowest)))

       if( this%hessianeProjected%eigenValues(lowest) < 0.0 .and. &
            abs(this%gradientProjectedOnHessiane%values(lowest)) < gradientThreshold ) then

          call TrustRegionOptimizer_exception(WARNING," The gradient for the more lowest mode is more low, using the next", &
               "Class object TrustRegionOptimizer in calculateStep() function")

          lowest=lowest+1
          if ( lowest > this%numberOfVariables ) then
             call TrustRegionOptimizer_exception(ERROR,"All gradients are too low", &
                  "Class object TrustRegionOptimizer in calculateStep() function")
          end if
          cycle lowest_search
       else
          exit lowest_search
       end if
    end do lowest_search
    !!
    !!***************************************************************************
    
    this%oldLagrangeMultiplier = 0.0
    this%lagrangeMultiplier=0.0
    flagA = .false.
    flagB = .false.

    auxVal = 0.0
    do i=1,this%numberOfVariables
       if( i /= negativeLowest .and. this%hessianeProjected%eigenValues(i) < 0.0 ) &
            auxVal= this%hessianeProjected%eigenValues(i)
    end do

    !! FALTA ADICINAR RESTRICCIONES  (A.2) AQUI

    if( .not.this%hessianeProjected%isPositiveDefinited ) then

       if( auxEigenVal < 0.0 .and. auxVal >= 0.0) goto 100

    else
       if(this%hessianeProjected%eigenValues(lowest) >= 0.0 ) goto 100
    end if

    do j=1,this%numberOfVariables
       print *,"valor propio: ",j,this%hessianeProjected%eigenValues(j)
    end do

    goto 200

    !!***************************************************************************
    !!   Calculo el tamanho del paso
    !!*************

100 step_search : do

       this%step%values=0.0
       do i=1,this%numberOfVariables
          if( abs(this%lagrangeMultiplier) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
               abs( this%hessianeProjected%eigenValues(i) ) < 1.0D-6 ) then
             auxVal = 0.0
          else
             if( abs( this%LagrangeMultiplier - this%hessianeProjected%eigenValues(i) ) < 1.0D-6) then
                auxVal=0.0
             else
                auxVal=this%gradientProjectedOnHessiane%values(i) &
                     /( this%LagrangeMultiplier - this%hessianeProjected%eigenValues(i) )
             end if
          end if

          if (  i == negativeLowest) then
             if ( abs(this%oldLagrangeMultiplier- auxEigenVal) < 1.0D-6) then
                auxVal = 0.0
             else
                auxVal=this%gradientProjectedOnHessiane%values(negativeLowest) &
                     /( this%oldLagrangeMultiplier - this%hessianeProjected%eigenValues(negativeLowest) )
             end if
          end if

          do j=1,this%numberOfVariables
             this%step%values(j)=this%step%values(j)+auxVal*this%hessianeProjected%eigenVectors(j,i)
          end do

       end do

       auxVal=sqrt(dot_product(this%step%values,this%step%values))

       if (auxVal < ( this%adaptativeMaximumTrustRadius + 1.0D-6 ) ) then
          this%xlambda=this%lagrangeMultiplier
          this%oldXlambda=this%oldLagrangeMultiplier
          this%trustRadius=auxVal
          return
       end if

       if( abs(this%lagrangeMultiplier) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD &
            .and. abs(this%oldLagrangeMultiplier) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then

          if( flagC ) then
             this%xlambda=this%lagrangeMultiplier
             this%oldXlambda=this%oldLagrangeMultiplier
             this%trustRadius=auxVal
             return
          end if

       end if
       
200    this%lagrangeMultiplier = 0.0
       flagA = .false.
       flagB = .false.
       BULI = this%hessianeProjected%eigenValues(lowest)
       step = 0.001
       if (  -auxEigenVal <  BULI ) BULI=-auxEigenVal
       if ( BULI > 0.0 ) BULI= 0.0
       if ( BULI <= 0.0 ) this%lagrangeMultiplier= BULI- step
       this%lagrangeMultiplier = BULI - step
       BL= this%lagrangeMultiplier - step
       BU= this%lagrangeMultiplier + 0.5 * step

       numberOfSearches = 0

       !!*******************************************************************************
       !! Realiza el calculo del multiplicador de Lagrange
       !!************
       lagrangeMultiplier_search : do
          FL= 0.0
          FU= 0.0

          do i=1, this%numberOfVariables

             if ( i /= negativeLowest) then
                FL = FL + ( this%gradientProjectedOnHessiane%values(i) / &
                     (BL - this%hessianeProjected%eigenValues(i) ) )**2.0
                FU = FU + ( this%gradientProjectedOnHessiane%values(i) / &
                     (BU - this%hessianeProjected%eigenValues(i) ) )**2.0
             end if

          end do

          if ( .not.this%hessianeProjected%isPositiveDefinited ) then
             FL = FL + ( this%gradientProjectedOnHessiane%values(negativeLowest) / &
                  (BL + this%hessianeProjected%eigenValues(negativeLowest) ) )**2.0
             FU = FU + ( this%gradientProjectedOnHessiane%values(negativeLowest) / &
                  (BU + this%hessianeProjected%eigenValues(negativeLowest) ) )**2.0
          end if

          FL = FL - this%adaptativeMaximumTrustRadius**2
          FU = FU - this%adaptativeMaximumTrustRadius**2

          if ( FL*FU < 0.0 ) then

             write (*, "(T10,A40,I5,A11)") "THE LAGRANGE MULTIPLIER CONVERGED AFTER OF ", numberOfSearches," ITERATIONS"
             exit lagrangeMultiplier_search
          else

             BL = BL - ( this%hessianeProjected%eigenValues(lowest) - BL )
             BU = BU + 0.5 * ( BULI - BU)

             if (  BL <= -this%maximumStepSize ) then
                BL= -this%maximumStepSize
                flagA=.true.
             end if

             if ( abs(BULI-BU) <= this%minimumStepSize ) then
                BU = BULI - this%minimumStepSize
                flagB = .true.
             end if

             if( flagA .and. flagB ) then
                this%lagrangeMultiplier = this%hessianeProjected%eigenValues(lowest) - 0.01
                this%oldLagrangeMultiplier = auxEigenVal + 0.01
                cycle step_search
             end if
             numberOfSearches = numberOfSearches +1

             if (numberOfSearches < 1000 )    then
                cycle lagrangeMultiplier_search
             else
                call TrustRegionOptimizer_exception(ERROR,"Lagrange Multiplier search has failed","Class object TrustRegionOptimizer in calculateStep() function")
             end if
          end if
       end do lagrangeMultiplier_search
       !!
       !!*******************************************************************************
       
       numberOfSearches =0
       this%xlambda = 0.0

       lambda_search : &
            do
       FL = 0.0
       FU = 0.0
       FM = 0.0
       this%lagrangeMultiplier = 0.5 * (BL + BU)

       do i=1,this%numberOfVariables
          if ( i /= negativeLowest) then
             FL = FL + ( this%gradientProjectedOnHessiane%values(i) / &
                  ( BL - this%hessianeProjected%eigenValues(i)) )** 2.0
             FU = FU + ( this%gradientProjectedOnHessiane%values(i) / &
                  ( BU - this%hessianeProjected%eigenValues(i)) )** 2.0
             FM = FM + ( this%gradientProjectedOnHessiane%values(i) / &
                  ( this%lagrangeMultiplier - this%hessianeProjected%eigenValues(i)) )** 2.0
          end if
       end do

       if ( .not.this%hessianeProjected%isPositiveDefinited) then
          FL = FL + ( this%gradientProjectedOnHessiane%values(negativeLowest) / &
               ( BL + this%hessianeProjected%eigenValues(negativeLowest)) )** 2.0
          FU = FU + ( this%gradientProjectedOnHessiane%values(negativeLowest) / &
               ( BU + this%hessianeProjected%eigenValues(negativeLowest)) )** 2.0
          FM = FM + ( this%gradientProjectedOnHessiane%values(negativeLowest) / &
               ( this%lagrangeMultiplier + this%hessianeProjected%eigenValues(negativeLowest)) )** 2.0
       end if

       FL = FL - (this%adaptativeMaximumTrustRadius**2)
       FU = FU - (this%adaptativeMaximumTrustRadius**2)
       FM = FM - (this%adaptativeMaximumTrustRadius**2)

       if( abs(this%xlambda-this%lagrangeMultiplier) < 1.0D-10) exit lambda_search

       numberOfSearches = numberOfSearches +1
       if( numberOfSearches < 1000 )  then
          this%xlambda = this%lagrangeMultiplier
          if ( FM*FU < 0.0) BL =this%lagrangeMultiplier
          if ( FM*FL < 0.0) BU =this%lagrangeMultiplier
          cycle lambda_search
       else

          call TrustRegionOptimizer_exception(ERROR,"Lagrange Multiplier search has failed", &
               "Class object TrustRegionOptimizer in calculateStep() function")

       end if

    end do lambda_search

    this%oldLagrangeMultiplier=-this%lagrangeMultiplier
    flagC = .true.
    cycle step_search
 end do step_search

end subroutine TrustRegionOptimizer_calculateStep

!>
!! @brief Calcula una prediccion al cambio de la funcion
!<
subroutine TrustRegionOptimizer_predictedChangeOfFunction( this )
  implicit none
  type(TrustRegionOptimizer) :: this

  real(8), allocatable :: auxVector(:)
  real(8) :: scaleFactor
  real(8) :: step
  real(8) :: auxVal
  integer :: i

  scaleFactor = 1.0
  if( this%trustRadius > this%adaptativeMaximumTrustRadius) &
       scaleFactor= this%adaptativeMaximumTrustRadius/this%trustRadius

  this%predictedChangeOfFunction = 0.0
  do i=1,this%numberOfVariables
     auxVal= this%xlambda
     if ( .not.this%hessianeProjected%isPositiveDefinited .and. i == this%principalMode ) auxVal= this%oldXlambda

     if ( abs( auxVal-this%hessianeProjected%eigenValues(i) ) < 1.0D-6 ) then
        step = 0.0
     else
        step=this%gradientProjectedOnHessiane%values(i)/(auxVal-this%hessianeProjected%eigenValues(i))
     end if
     step = step*scaleFactor
     auxVal=step*this%gradientProjectedOnHessiane%values(i) + 0.5*(step**2)*this%hessianeProjected%eigenValues(i)
     this%predictedChangeOfFunction=this%predictedChangeOfFunction+auxVal
  end do

  this%oldFunctionValue = this%functionValue
  this%oldNumberOfIteration = this%numberOfIteration

end subroutine TrustRegionOptimizer_predictedChangeOfFunction

!>
!! @brief Actualiza en valor de radio de confianza
!<
subroutine TrustRegionOptimizer_updateTrustRadio( this )
  implicit none
  type(TrustRegionOptimizer) :: this

  real(8), parameter :: upperBondForRatio = 0.75_8
  real(8), parameter :: lowerBondForRatio = 0.25_8
  real(8), parameter :: radiusFactor = 2.0_8
  real(8) :: scaleFactor

  this%isSuitableStep = .true.

  !! No altera el tamano del radio para cambio pequenos de energia
  !       if ( ( this%realChangeOfFunction < 0.0_8) .and. (abs(this%realChangeOfFunction)< 1.0D-5) ) return

  !!**********************************************
  !! Modifica el valor del radio de confianza
  !!****
  if ( this%ratioOfChange < 0.0_8 ) this%isSuitableStep = .false.

  this%trustRadius = sqrt( dot_product(this%step%values,this%step%values) )
  if( this%trustRadius > this%adaptativeMaximumTrustRadius+ 1.0D-06 ) then
     scaleFactor = this%adaptativeMaximumTrustRadius/this%trustRadius
     this%trustRadius = this%adaptativeMaximumTrustRadius
     this%step%values= scaleFactor*this%step%values
  end if


  !       if ( this%ratioOfChange > upperBondForRatio ) then
  !           if ( sqrt(this%squareStepSize) > ( 0.8_8 * this%trustRadius ) ) this%trustRadius = radiusFactor * this%trustRadius
  !       else &
  !       if ( ( this%ratioOfChange >= lowerBondForRatio) .and. ( this%ratioOfChange <= upperBondForRatio ) ) then
  !           this%trustRadius = this%trustRadius
  !       else
  !           this%trustRadius = this%trustRadius / radiusFactor
  !       end if
  !!
  !!**********************************************

  !       !! Establece limite inferior para cambio de radio de confianza
  !       if ( (this%trustRadius < MINIMUM_TRUST_RADIO ) .and. ( .not.this%isSuitableStep ) ) then
  !           this%isSuitableStep = .true.
  !           write(*, *) "Can not reduce the trust radius anymore, forcing to accept it"
  !       end if
  !
  !       !! Acota el radio de confianza
  !       this%trustRadius = max( this%trustRadius, MINIMUM_TRUST_RADIO )
  !       this%trustRadius = min( this%trustRadius, MAXIMUM_TRUST_RADIO )
  write(*, "(T10,A15,F12.8)") "TRUST RADIUS = ",this%trustRadius

end subroutine TrustRegionOptimizer_updateTrustRadio

  !>
  !! @brief Actualiza la matriz hesiana empleando la
  !! formula de Broyden-Fletcher-Goldfab-Shano (BFGS)
  !<
  subroutine TrustRegionOptimizer_updateHessiane( this )
    implicit none
    type( TrustRegionOptimizer ) :: this

    type(Vector) :: changeOfGradient
    type(Vector) :: auxVector
    real(8) :: factor(4)
    integer :: i
    integer :: j

    call Vector_constructor( changeOfGradient, this%numberOfVariables )
    call Vector_constructor( auxVector, this%numberOfVariables )

    changeOfGradient%values= this%gradient%values - this%oldGradient%values

    auxVector%values=matmul(this%hessiane%values,this%step%values)
    changeOfGradient%values=this%gradient%values-this%oldGradient%values


    factor(1)=dot_product(changeOfGradient%values,this%step%values)
    factor(2)=dot_product(auxVector%values,this%step%values)
    factor(3)=dot_product(changeOfGradient%values,changeOfGradient%values)
    factor(4)=dot_product(this%step%values,this%step%values)

    if ( factor(1) < sqrt(3.0D-8*factor(3)*factor(4) ) ) then
       write(*, "(T10,A33)") "THE HESSIANE HAS NOT BEEN UPDATED"
    else
       if( factor(1) < 0.0 ) then
          call TrustRegionOptimizer_exception(WARNING,"POSITIVE DEFINITENESS NOT IS WARRANTED","")
       end if

       do i=1, this%numberOfVariables
          do j=1,i
             this%hessiane%values(i,j)=this%hessiane%values(i,j) + &
                  ( changeOfGradient%values(i)*changeOfGradient%values(j)/factor(1) ) - &
                  ( auxVector%values(i) * auxVector%values(j)/ factor(2) )
             if( abs(this%hessiane%values(i,j)) < 3.0D-6 ) this%hessiane%values(i,j)=0.0
             this%hessiane%values(j,i)=this%hessiane%values(i,j)
          end do
       end do
    end if
    call Vector_destructor( changeOfGradient )
    call Vector_destructor( auxVector )
  end subroutine TrustRegionOptimizer_updateHessiane

! !>
! !! @brief Reinicia el iterador tomando como punto de partida el estimado actual del minimo
! !!      y la hessiana actual
! !>
! subroutine TrustRegionOptimizer_restart( this )
! implicit none
! type(TrustRegionOptimizer) :: this



! type(Vector) :: initialpoint
! type(Vector) :: auxGradient
! integer :: currentIterations


! initialPoint= this%variables
! auxGradient = this%gradient
! currentIterations=this%realNumberOfIteration

! write (*,*) ""
! write (*,"(T10,A52)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
! write (*,"(T10,A52)") "!!!!  RESTARTING THE ALGORITHM  OF MINIMIZATION !!!!"
! write (*,"(T10,A52)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
! write (*,*) ""

! call TrustRegionOptimizer_destructor( this )
! call TrustRegionOptimizer_constructor(this, initialPoint%values)
! !!      call Matrix_setIdentity( this%hessiane )
! this%hasRestarted=.true.
! this%realNumberOfIteration=currentIterations
! this%gradient=auxGradient
! call Vector_destructor( initialPoint )
! call Vector_destructor( auxGradient )


! end subroutine TrustRegionOptimizer_restart

!>
!! @brief  Maneja excepciones de la clase
!<
subroutine TrustRegionOptimizer_exception( typeMessage, description, debugDescription)
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

end subroutine TrustRegionOptimizer_exception


end module TrustRegionOptimizer_
