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
!! @brief Energy Gradients Module.
!!        This module contains all basic functions of the energy gradients calculations
!! @author  J.M. Rodas
!! 
!! <b> Creation date : </b> 2015-02-27
!!
!! <b> History: </b>
!!
!!   - <tt> 2015-02-27 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0 and prepare the module for new methods
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!

module EnergyGradients_
#ifdef intel
  use IFPORT
#endif
  use CONTROL_
  use Matrix_
  use MatrixInteger_
  use Vector_
  use MolecularSystem_
  use ParticleManager_
  use ContractedGaussian_
  use CosmoCore_
  use MecanicProperties_
  use DerivativeManager_
  use Math_
  use Exception_
  implicit none
  ! use omp_lib

  type, public :: Gradients

     real(8), allocatable :: kinetic(:)
     real(8), allocatable :: attraction(:)
     real(8), allocatable :: overlap(:)
     real(8), allocatable :: coulomb(:)
     real(8), allocatable :: exchange(:)
     real(8), allocatable :: nuclear(:)
     real(8), allocatable :: coupling(:)
     real(8), allocatable :: total(:)

  end type Gradients


  type, public :: EnergyGradients

     character(30) :: name
     logical :: isAnalytic
     logical :: isIntanced
     real(8) :: deltaX
     real(8) :: output
     integer :: order
     type(Vector) pointsOfEvaluation
     integer, allocatable :: components(:)
     logical :: derivativeWithMP
     type(Gradients) :: gradients 

  end type EnergyGradients

  type(EnergyGradients), public :: EnergyGradients_instance

  private :: &
       EnergyGradients_calculateAnalyticUncoupledFirstDerivative, &
       EnergyGradients_calculateAnalyticCouplingFirstDerivative, &
       EnergyGradients_calculateFistDerivativeOfPuntualEnergy

  public :: &
       EnergyGradients_constructor, &
       EnergyGradients_destructor, &
                                ! EnergyGradients_setName, &
                                ! EnergyGradients_setAnalytic, &
                                ! EnergyGradients_setNumeric, &
                                ! EnergyGradients_setNumericStepSize, &
       EnergyGradients_getDerivative, &
       EnergyGradients_getAnalyticDerivative, &
       EnergyGradients_getNumericGradient, &
       EnergyGradients_show
  ! EnergyGradients_getNumericDerivative, &
  ! EnergyGradients_getTypeGradient, &
  ! EnergyGradients_isInstanced

contains
  !**
  ! @brief Define el constructor
  !**
  subroutine EnergyGradients_constructor()
    implicit none


    EnergyGradients_instance%name = ""
    EnergyGradients_instance%isIntanced = .true.
    EnergyGradients_instance%isAnalytic = CONTROL_instance%ANALYTIC_GRADIENT
    EnergyGradients_instance%deltaX = CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA
    EnergyGradients_instance%derivativeWithMP = CONTROL_instance%OPTIMIZE_WITH_MP
    EnergyGradients_instance%order = 1

  end subroutine EnergyGradients_constructor

  !**
  ! @brief Define el destructor
  !**
  subroutine EnergyGradients_destructor()
    implicit none

    integer :: i
    integer :: status
    EnergyGradients_instance%isIntanced = .false.

  end subroutine EnergyGradients_destructor

  subroutine EnergyGradients_show()
    implicit none
    integer :: numberOfOptimizationCenters
    integer :: i,j,k

    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()


    write(*,"(A)") "----------------------------------------------------"
    write(*,"(A)") "            ENERGY GRADIENTS(HARTREE/BOHR)          "
    write(*,"(A)") "        dE/dx            dE/dy            dE/dz     "
    write(*,"(A)") "----------------------------------------------------"
    do i=0, numberOfOptimizationCenters - 1
       j = i*3 + 1
       k = i*3 + 3
       write(*,"(3F17.12)") EnergyGradients_instance%gradients%total(j:k)
    end do

    if(CONTROL_instance%AMBER_FILE) then
       !! Results for Amber package
       !!

       !! open file
       open(unit=45, file="low2amber.dat", status="replace", form="formatted")

       !!save all options
       write(45,*) "Este es un archivo de prueba"
       close(45)
    end if


  end subroutine EnergyGradients_show

  !         !**
  !         ! @brief Ajusta el nombre del minimizador empleado
  !         !**
  !         subroutine EnergyGradients_setName()
  !             implicit none


  !         end subroutine EnergyGradients_setName

  !         !**
  !         ! @brief Ajusta calculo del gradiente de forma numerica
  !         !**
  !         subroutine EnergyGradients_setNumeric()
  !             implicit none

  !             EnergyGradients_instance%isAnalytic = .false.

  !         end subroutine EnergyGradients_setNumeric

  !         !**
  !         ! @brief Ajusta calculo del gradiente de forma analitic
  !         !**
  !         subroutine EnergyGradients_setAnalytic()
  !             implicit none

  !             EnergyGradients_instance%isAnalytic = .true.

  !         end subroutine EnergyGradients_setAnalytic

  !         !**
  !         ! @brief Ajusta calculo del gradiente de forma analitic
  !         !**
  !         subroutine EnergyGradients_setNumericStepSize(stepSize)
  !             implicit none
  !             real(8) :: stepSize

  !             EnergyGradients_instance%deltaX = stepSize

  !         end subroutine EnergyGradients_setNumericStepSize

  !         !**
  !         ! @brief Indica si el gradiente es numerico o analitico
  !         !**
  !         function EnergyGradients_getTypeGradient() result(output)
  !             implicit none
  !             character(30) :: output

  !             if ( EnergyGradients_instance%isAnalytic .eqv. .true. ) then
  !                 output = "ANALYTIC"
  !             else
  !                 output = "NUMERICAL"
  !             end if

  !         end function EnergyGradients_getTypeGradient


  !**
  ! @brief    Retorna el valor de la derivada de orden n en un punto dado
  !       de las variables independientes
  !**
  subroutine EnergyGradients_getDerivative( evaluationPoint, functionValue, gradientVector, order )
    implicit none
    type(Vector), target, intent(in) :: evaluationPoint
    real(8), intent(inout) :: gradientVector(:)
    integer, intent(in),optional :: order
    integer :: i, j

    interface
       function functionValue( pointOfEvaluation ) result( output )
         real(8) :: pointOfEvaluation(:)
         real(8) :: output
       end function functionValue
    end interface

    EnergyGradients_instance%order = 1
    if ( present(order) ) EnergyGradients_instance%order = order
    EnergyGradients_instance%pointsOfEvaluation = evaluationPoint
    !    allocate ( EnergyGradients_instance%components( size(components) ) )
    !    EnergyGradients_instance%components = components

    !! Debug Mauricio Rodas
    ! write(*,"(A)") "Componenetes de entrada"
    ! do i=1, size(components)
    !    write(*,"(I)") EnergyGradients_instance%components(i)
    ! write(*,"(A)") "-------------------------"
    ! do j=1, size(evaluationPoint%values)
    !    write(*,"(f12.8)") evaluationPoint%values(j)
    ! end do
    ! end do
    ! write(*,"(A)") "-------------------------"

    gradientVector=0.0_8

    ! if ( EnergyGradients_instance%isAnalytic ) then
    call EnergyGradients_getAnalyticDerivative(gradientVector)

    ! else
    ! if( EnergyGradients_instance%derivativeWithMP .and. CONTROL_instance%MOLLER_PLESSET_CORRECTION==2) then
    !    output = EnergyGradients_gradientWithMP()
    ! else
    ! output = EnergyGradients_getNumericDerivative( functionValue )
    ! end if
    ! end if

    !    if ( abs(output) < 1.0E-8 ) gradientVector=0.0D+0

    !! Deja el sistema en su posicion original
    call ParticleManager_setParticlesPositions( evaluationPoint )
    !    deallocate(EnergyGradients_instance%components)

  end subroutine EnergyGradients_getDerivative

  function EnergyGradients_getNumericGradient( evaluationPoint, functionValue ) result(output)
    implicit none
    type(Vector), target, intent(in) :: evaluationPoint
    real(8), allocatable :: output(:)

    interface
       function functionValue( pointOfEvaluation ) result( output )
         real(8) :: pointOfEvaluation(:)
         real(8) :: output
       end function functionValue
    end interface

    real(8) :: auxVal
    real(8) :: auxValue(2)
    type(Matrix) :: projector
    type(Matrix) :: gradientProjected
    integer :: infoProcess
    integer :: ssize
    integer :: numberOfSymmetricModes
    integer :: nucleiAndComponent
    integer :: i

    EnergyGradients_instance%pointsOfEvaluation = evaluationPoint

    allocate( output( size(evaluationPoint%values) ) )

    projector = MecanicProperties_getHandyProjector( MolecularSystem_instance%mechanicalProp, infoProcess )

    ssize=size(projector%values, dim=1)

    do i=1,ssize
       auxVal=sqrt( dot_product(projector%values(i,:),projector%values(i,:)) )
       if( auxVal > 1.0D-4 )  projector%values(i,:)=projector%values(i,:)/auxVal
    end do

    !! Verifica la independencia lineal de los vetores
    projector%values=transpose(projector%values)
    call Matrix_selectLinearlyIndependentVectors(projector, numberOfSymmetricModes)
    projector%values=transpose(projector%values)

    output=0.0
    do i=1, numberOfSymmetricModes

       ! if( CONTROL_instance%OPTIMIZE_WITH_MP ) then

       !    !! Ajusta punto x+h
       !    EnergyGradients_instance%pointsOfEvaluation%values = EnergyGradients_instance%deltaX*projector%values(i,:) &
       !         + evaluationPoint%values

       !    call ParticleManager_setParticlesPositions( EnergyGradients_instance%pointsOfEvaluation )
       !    call EnergyGradients_writeInput(trim(CONTROL_instance%INPUT_FILE)//"0.apmo")

       !    !! Ajusta punto x-h
       !    EnergyGradients_instance%pointsOfEvaluation%values = -EnergyGradients_instance%deltaX*projector%values(i,:) &
       !         + evaluationPoint%values

       !    call ParticleManager_setParticlesPositions( EnergyGradients_instance%pointsOfEvaluation )
       !    call EnergyGradients_writeInput(trim(CONTROL_instance%INPUT_FILE)//"1.apmo")

       !    call EnergyGradients_calculateFile(trim(CONTROL_instance%INPUT_FILE),auxValue)

       ! else

       EnergyGradients_instance%pointsOfEvaluation%values = EnergyGradients_instance%deltaX*projector%values(i,:) &
            + evaluationPoint%values

       auxValue(1) = functionValue( EnergyGradients_instance%pointsOfEvaluation%values )

       EnergyGradients_instance%pointsOfEvaluation%values = -EnergyGradients_instance%deltaX*projector%values(i,:) &
            + evaluationPoint%values

       auxValue(2) = functionValue( EnergyGradients_instance%pointsOfEvaluation%values )

       ! end if

       ! print "(A1$)",":"
       output(i) = ( auxValue(1) - auxValue(2) ) / ( 2.0_8 * EnergyGradients_instance%deltaX )

    end do

    output=matmul( transpose(projector%values), output)

    do i=1,size(output)
       if( abs(output(i)) < 1.0D-8 ) output(i)=0.0
    end do

    call ParticleManager_setParticlesPositions( evaluationPoint )
    call Matrix_destructor(projector)

  end function EnergyGradients_getNumericGradient

  !**
  ! @brief    Retorna el valor de la derivada de orden n en un punto dado
  !       de las variables independientes
  !**
  subroutine EnergyGradients_getAnalyticDerivative(gradientVector)
    implicit none
    real(8), intent(inout), optional :: gradientVector(:)
    ! real(8), allocatable :: auxGradient(:)
    type(Exception) :: ex
    integer :: numberOfOptimizationCenters
    type(Vector) :: geometry
    integer :: i,j,k 
    type(surfaceSegment) :: surface_aux
    ! real(8), allocatable :: auxNuclear(:,:), auxCoupling(:,:)

    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if(CONTROL_instance%COSMO)then
       call CosmoCore_lines(surface_aux)
       call CosmoCore_filler(surface_aux)
    end if


    ! geometry = ParticleManager_getPositionOfCenterOfOptimizacion()

    ! write(*,"(A,I)") "Number of opt: ", numberOfOptimizationCenters
    ! write(*,"(3F17.12)") geometry%values(1:3)
    ! write(*,"(3F17.12)") geometry%values(4:6)
    ! write(*,"(3F17.12)") geometry%values(7:9)
    ! write(*,"(3F17.12)") geometry%values(10:12)
    ! write(*,"(3F17.12)") geometry%values(13:15)

    ! write(*,"(A,I)") "Entrando a calcular Gradientes ..."
    if(allocated(EnergyGradients_instance%gradients%coupling)) deallocate(EnergyGradients_instance%gradients%coupling)
    allocate(EnergyGradients_instance%gradients%coupling(numberOfOptimizationCenters*3))

    if(allocated(EnergyGradients_instance%gradients%total)) deallocate(EnergyGradients_instance%gradients%total)
    allocate(EnergyGradients_instance%gradients%total(numberOfOptimizationCenters*3))

    EnergyGradients_instance%gradients%total = 0.0_8
    EnergyGradients_instance%gradients%coupling = 0.0_8

    if ( EnergyGradients_instance%order == 1 ) then
      
       if(CONTROL_instance%COSMO)then
          call EnergyGradients_calculateAnalyticUncoupledFirstDerivative(surface_aux)
       else
          call EnergyGradients_calculateAnalyticUncoupledFirstDerivative()
       end if

       if ( MolecularSystem_instance%numberOfQuantumSpecies > 1) then
          call EnergyGradients_calculateAnalyticCouplingFirstDerivative()
       end if
       
       call EnergyGradients_calculateFistDerivativeOfPuntualEnergy()

       k=1
       do i=1, numberOfOptimizationCenters
          do j=1,3
             EnergyGradients_instance%gradients%total(k) = EnergyGradients_instance%gradients%total(k) + &
                  EnergyGradients_instance%gradients%nuclear(k) + EnergyGradients_instance%gradients%coupling(k)
             if(present(gradientVector)) then
                gradientVector(k) = EnergyGradients_instance%gradients%total(k)
             end if
             k = k + 1
          end do
       end do

       if(present(gradientVector)) then
          ! write(*,"(A)") "----------------------------------------------------------------"
          ! write(*,"(A)") " Gradientes Totales"
          ! write(*,"(A)") "----------------------------------------------------------------"
          ! write(*,"(3(f17.12))") gradientVector(1), gradientVector(2), gradientVector(3)
          ! write(*,"(3(f17.12))") gradientVector(4), gradientVector(5), gradientVector(6)
          ! write(*,"(3(f17.12))") gradientVector(7), gradientVector(8), gradientVector(9)
          ! write(*,"(3(f17.12))") gradientVector(10), gradientVector(11), gradientVector(12)
          ! write(*,"(3(f17.12))") gradientVector(13), gradientVector(14), gradientVector(15)
          ! write(*,"(A)") "----------------------------------------------------------------"
       else
          call EnergyGradients_show()
       end if

       deallocate(EnergyGradients_instance%gradients%total)
       deallocate(EnergyGradients_instance%gradients%nuclear)
       deallocate(EnergyGradients_instance%gradients%kinetic)
       deallocate(EnergyGradients_instance%gradients%attraction)
       deallocate(EnergyGradients_instance%gradients%overlap)
       deallocate(EnergyGradients_instance%gradients%coupling)
       deallocate(EnergyGradients_instance%gradients%coulomb)
       deallocate(EnergyGradients_instance%gradients%exchange)
    else

       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object EnergyGradients in the getAnalyticDerivative function" )
       call Exception_setDescription( ex, "This order isn't implemented" )
       call Exception_show( ex )

    end if

  end subroutine EnergyGradients_getAnalyticDerivative

  ! !**
  ! ! @brief    Retorna el valor de la derivada de orden n en un punto dado
  ! !       de las variables independientes
  ! !**
  ! function EnergyGradients_getNumericDerivative( functionValue ) result( output )
  !   implicit none
  !   real(8) :: output

  !   interface
  !      function functionValue( pointOfEvaluation ) result( output )
  !        real(8) :: pointOfEvaluation(:)
  !        real(8) :: output
  !      end function functionValue
  !   end interface

  !   type(Exception) :: ex
  !   type(Vector) :: auxEvaluationPoint
  !   real(8), allocatable :: auxValue(:)
  !   real(8) :: energyDeltaLess

  !   integer :: nucleiAndComponent(2)

  !   output=0.0_8

  !   if ( EnergyGradients_instance%order <= 2 ) then

  !      select case( EnergyGradients_instance%order )

  !      case(1)

  !         nucleiAndComponent= ParticleManager_getCenterOfOptimization(EnergyGradients_instance%components(1) )
  !         if ( ParticleManager_isComponentFixed( nucleiAndComponent(1),nucleiAndComponent(2) ) ) return

  !         !                       if( abs(EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) ) < &
  !         !                           1.0D-10 ) return

  !         allocate( auxValue(2) )

  !         !! Ajusta punto x+h
  !         EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !              EnergyGradients_instance%pointsOfEvaluation%values( &
  !              EnergyGradients_instance%components(1) ) &
  !              + EnergyGradients_instance%deltaX

  !         auxValue(1) = functionValue( EnergyGradients_instance%pointsOfEvaluation%values )

  !         !! Ajusta punto x-h
  !         EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !              EnergyGradients_instance%pointsOfEvaluation%values( &
  !              EnergyGradients_instance%components(1) ) &
  !              - 2.0_8 * EnergyGradients_instance%deltaX

  !         auxValue(2) = functionValue( &
  !              EnergyGradients_instance%pointsOfEvaluation%values )

  !         !! calcula primera derivada de la energia
  !         output = ( auxValue(1) - auxValue(2) ) / ( 2.0_8 * EnergyGradients_instance%deltaX )

  !         deallocate( auxValue )

  !      case(2)


  !         if ( EnergyGradients_instance%components(1) == EnergyGradients_instance%components(2) ) then

  !            allocate( auxValue(5) )

  !            !! Ajusta en punto x-2h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 -2.0_8 * EnergyGradients_instance%deltaX

  !            auxValue(1) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x-h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 + EnergyGradients_instance%deltaX

  !            auxValue(2) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto inicial
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 + EnergyGradients_instance%deltaX

  !            auxValue(3) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x+h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 + EnergyGradients_instance%deltaX

  !            auxValue(4) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x+2*h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 + EnergyGradients_instance%deltaX

  !            auxValue(5) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Calcula segunda derivada d^2/dx^2
  !            output = ( 16.0_8 * ( auxValue(4) - auxValue(2) ) - 30.0_8 * auxValue(3) - auxValue(1) - auxValue(5) ) &
  !                 / ( 12.0_8 * ( EnergyGradients_instance%deltaX**2.0_8 ) )

  !         else

  !            allocate( auxValue(4) )

  !            !! Ajusta en punto x-h, y-h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 - EnergyGradients_instance%deltaX

  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(2) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(2) ) &
  !                 - EnergyGradients_instance%deltaX

  !            auxValue(1) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x-h, y+h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(2) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(2) ) &
  !                 + 2.0_8 * EnergyGradients_instance%deltaX

  !            auxValue(2) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x+h, y+h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(1) ) &
  !                 + 2.0_8 * EnergyGradients_instance%deltaX

  !            auxValue(3) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Ajusta en punto x+h, y-h
  !            EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(2) ) = &
  !                 EnergyGradients_instance%pointsOfEvaluation%values( &
  !                 EnergyGradients_instance%components(2) ) &
  !                 - 2.0_8 * EnergyGradients_instance%deltaX

  !            auxValue(4) = functionValue( &
  !                 EnergyGradients_instance%pointsOfEvaluation%values )

  !            !! Calcula segunda derivada d^2/dxy
  !            output =( auxValue(3) - auxValue(4) - auxValue(2) + auxValue(1) ) &
  !                 / ( 4.0_8 * (EnergyGradients_instance%deltaX**2.0_8) )

  !         end if

  !      end select

  !   else

  !      call Exception_constructor( ex , ERROR )
  !      call Exception_setDebugDescription( ex, "Class object EnergyGradients in the getNumericDerivative function" )
  !      call Exception_setDescription( ex, "This order isn't implemented" )
  !      call Exception_show( ex )

  !   end if

  ! end function EnergyGradients_getNumericDerivative

  !         function EnergyGradients_gradientWithMP() result(output)
  !             implicit none
  !             real(8) :: output

  !             integer :: status
  !             real(8) :: auxValue(2)
  !             integer :: nucleiAndComponent(2)

  !             output = 0.0_8

  !             nucleiAndComponent= ParticleManager_getCenterOfOptimization(EnergyGradients_instance%components(1) )
  !             if ( ParticleManager_isComponentFixed( nucleiAndComponent(1),nucleiAndComponent(2) ) ) return

  !             !! Ajusta punto x+h
  !             EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                                         EnergyGradients_instance%pointsOfEvaluation%values( &
  !                                         EnergyGradients_instance%components(1) ) &
  !                                         + EnergyGradients_instance%deltaX

  !             call ParticleManager_setParticlesPositions( EnergyGradients_instance%pointsOfEvaluation )

  !             call EnergyGradients_writeInput(trim(CONTROL_instance%INPUT_FILE)//"0.apmo")

  !             !! Ajusta punto x-h
  !             EnergyGradients_instance%pointsOfEvaluation%values( EnergyGradients_instance%components(1) ) = &
  !                                         EnergyGradients_instance%pointsOfEvaluation%values( &
  !                                         EnergyGradients_instance%components(1) ) &
  !                                         - 2.0_8 * EnergyGradients_instance%deltaX

  !             call ParticleManager_setParticlesPositions( EnergyGradients_instance%pointsOfEvaluation )

  !             call EnergyGradients_writeInput(trim(CONTROL_instance%INPUT_FILE)//"1.apmo")

  !             call EnergyGradients_calculateFile(trim(CONTROL_instance%INPUT_FILE),auxValue)

  !             !! calcula primera derivada de la energia
  !             output = ( auxValue(1) - auxValue(2) ) / ( 2.0_8 * EnergyGradients_instance%deltaX )

  !         end function EnergyGradients_gradientWithMP

  !         subroutine EnergyGradients_writeInput(nameOfFile)
  !             implicit none
  !             character(*) :: nameOfFile
  !             integer :: i

  !             open( UNIT=34,FILE=trim(nameOfFile),STATUS='REPLACE', &
  !                 ACCESS='SEQUENTIAL', FORM='FORMATTED' )

  !             write (34,*) "SYSTEM_DESCRIPTION='' "
  !             write (34,*) ""
  !             write (34,*) "GEOMETRY"
  !                 do i=1,size(  ParticleManager_instance%particles )
  ! #ifdef intel
  !                     write(34,"(A10,A10,<3>F15.10)")     trim(ParticleManager_instance%particles(i)%nickname), &
  !                                         trim(ParticleManager_instance%particles(i)%basisSetName), &
  !                                         ParticleManager_instance%particles(i)%origin*0.52917724924_8
  ! #else
  !                     write(34,"(A10,A10,3F15.10)")     trim(ParticleManager_instance%particles(i)%nickname), &
  !                                         trim(ParticleManager_instance%particles(i)%basisSetName), &
  !                                         ParticleManager_instance%particles(i)%origin*0.52917724924_8
  ! #endif
  !                 end do
  !             write (34,*) "END GEOMETRY"
  !             write (34,*) ""
  !             write (34,*) "TASKS"
  !                 write (34,*) "     method=RHF"
  !                 write (34,*) "     mollerPlessetCorrection=2"
  !             write (34,*) "END TASKS"
  !             write (34,*) ""
  !             write (34,*) "CONTROL"
  !                 write (34,*) "     scfnonelectronicenergytolerance =",CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE
  !                 write (34,*) "     scfelectronicenergytolerance =",CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE
  !                 write (34,*) "     nonelectronicdensitymatrixtolerance =",CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
  !                 write (34,*) "     electronicdensitymatrixtolerance =",CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
  !                 write (34,*) "     totalenergytolerance =",CONTROL_instance%TOTAL_ENERGY_TOLERANCE
  !                 write (34,*) "     strongenergytolerance =",CONTROL_instance%STRONG_ENERGY_TOLERANCE
  !                 write (34,*) "     densityfactorthreshold =",CONTROL_instance%DENSITY_FACTOR_THRESHOLD
  !                 write (34,*) "     scfnonelectronicmaxiterations =",CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
  !                 write (34,*) "     scfelectronicmaxiterations =",CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
  !                 write (34,*) "     scfmaxiterations =",CONTROL_instance%SCF_MAX_ITERATIONS
  !                 write (34,*) "     scfinterspeciesmaximumiterations =",CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS
  !                 write (34,*) "     diisswitchthreshold =",CONTROL_instance%DIIS_SWITCH_THRESHOLD
  !                 write (34,*) "     convergencemethod =",CONTROL_instance%CONVERGENCE_METHOD
  !             write (34,*) "END CONTROL"
  !             close(34)

  !         end subroutine EnergyGradients_writeInput

  !         subroutine EnergyGradients_calculateFile(fileName, values)
  !             implicit none
  !             character(*) :: fileName
  !             real(8) :: values(2)

  !             real(8) :: aux
  !             integer :: status
  !             integer :: i
  !             integer :: auxValue
  !             character :: auxChar
  !             character(255) :: line


  !             do i=1, ParticleManager_getNumberOfQuantumSpecies()
  !                 if ( ParticleManager_instance%particles(i)%isQuantum ) then
  !                     status=system("cp "//trim(fileName)//trim(ParticleManager_getNameOfSpecie( i ))//".vec " &
  !                             //trim(fileName)//"0."//trim(ParticleManager_getNameOfSpecie( i ))//".vec ")
  !                     status=system("cp "//trim(fileName)//trim(ParticleManager_getNameOfSpecie( i ))//".vec " &
  !                             //trim(fileName)//"1."//trim(ParticleManager_getNameOfSpecie( i ))//".vec ")
  !                 end if
  !             end do

  !             call omp_set_num_threads(2)
  !             !$OMP PARALLEL private(status)

  !                 if ( omp_get_thread_num()==0) then
  !                     status=system("apmo "//trim(fileName)//"0.apmo")
  !                 else if ( omp_get_thread_num()==1) then
  !                     status=system("apmo "//trim(fileName)//"1.apmo")
  !                 end if

  !             !$OMP END PARALLEL

  !             !!     lee el valor de la energia

  !             open(UNIT=70,FILE=trim(fileName)//"0.out",STATUS='unknown',ACCESS='sequential')
  !             read(70,'(A)') line
  !             line=trim(adjustl(trim(line)))

  !             do while((line(1:13)) /= "Enlapsed Time")
  !                 read(70,'(A)',end=10) line
  !                 line=trim(adjustl(trim(line)))

  ! #ifdef intel
  !                 if ( (line(1:7)) == "E(MP2)=" )  &
  !                     values(1) = dnum(line(scan(trim(line),"=" )+1:len_trim(line)))
  ! #else
  !                 if ( (line(1:7)) == "E(MP2)=" )  then
  !                     write(line(scan(trim(line),"=" )+1:len_trim(line)), *) aux
  !                     values(1) = aux
  !                 end if
  ! #endif

  !             end do
  !             close(70)
  !             status=system("rm "//trim(fileName)//"0.*")
  !             open(UNIT=70,FILE=trim(fileName)//"1.out",STATUS='unknown',ACCESS='sequential')
  !             read(70,'(A)') line
  !             line=trim(adjustl(trim(line)))

  !             do while((line(1:13)) /= "Enlapsed Time")
  !                 read(70,'(A)',end=10) line
  !                 line=trim(adjustl(trim(line)))

  ! #ifdef intel
  !                 if ( (line(1:7)) == "E(MP2)=" )  &
  !                     values(2) = dnum(line(scan(trim(line),"=" )+1:len_trim(line)))
  ! #else
  !                 if ( (line(1:7)) == "E(MP2)=" )  then
  !                     write(line(scan(trim(line),"=" )+1:len_trim(line)), *) aux
  !                     values(2) = aux
  !                 end if
  ! #endif

  !             end do
  !             close(70)
  !             status=system("rm "//trim(fileName)//"1.*")

  !             return

  ! 10          call EnergyGradients_exception(ERROR, "The single point calculation has failed", &
  !                 "Class object EnergyGradients in calculateFile() function")

  !         end subroutine EnergyGradients_calculateFile

  !         !**
  !         ! @brief    Indica si el derivadoor a sido o instanciado
  !         !
  !         !**
  !         function EnergyGradients_isInstanced() result(output)
  !             implicit none
  !             logical :: output

  !             output=EnergyGradients_instance%isIntanced


  !         end function EnergyGradients_isInstanced

  !**
  ! @brief  Calcula la primera derivada analitica para la matriz de Fock sin incluir el termino de
  !               acoplamiento interspecie
  ! 16/07/2012: revisada para implementacion de derivadas analiticas para cualquier momento angular (Edwin Posada)
  !**
  subroutine EnergyGradients_calculateAnalyticUncoupledFirstDerivative(surface)
    implicit none
    ! real(8), allocatable :: gradientVector(:)
    character(30) :: nameOfSpecie
    type(ContractedGaussian), allocatable :: contractions(:)
    type(surfaceSegment), intent(in), optional :: surface
    integer :: numberOfContractions
    real(8), allocatable :: auxVector(:), auxVector2(:), auxVector3(:), auxVector4(:)
    integer :: specieIterator
    integer :: ocupationNumber
    integer :: orderOfMatrix
    type(MatrixInteger) :: shellPairs
    type(Matrix) :: matrixOfEigenvectors
    type(Matrix) :: densityMatrix
    type(Vector) :: vectorOfEigenvalues
    type(Matrix) :: auxWeightDensity
    type(Matrix) :: weightDensityMatrix
    type(Matrix) :: twoParticles
    integer(8) :: numberOfIntegrals
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: u
    integer :: v
    integer :: x
    integer :: y
    integer :: npairs, npairs2, pqrsIter
    integer :: P, Q, R, S, PQ, RS, pIter, qIter, rIter, sIter, auxIter, A, B, stride
    integer :: numCartesianP, numCartesianQ, numCartesianR, numCartesianS
    integer :: centerP, centerQ, centerR, centerS, center
    integer :: numberOfOptimizationCenters, delta
    integer :: Zi
    real(8) :: perm, Vval, JKval, mass
    real(8) :: Duv, Dxy
    real(8) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz
    integer, allocatable :: labelsOfContractions(:)
    real(8), allocatable :: auxKinetic(:,:), auxPotential(:,:), auxOverlap(:,:), auxCoulomb(:,:), auxExchange(:,:)
    real(8), allocatable :: auxCOSMO(:,:), auxCOSMO2(:,:), auxCOSMO3(:,:,:)
    integer, allocatable :: auxOwnerId(:)
    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20)
    real(8) :: charge
    real(8) :: eta, lambda, kappa
    integer :: deltasum
    real(8) :: matrixSum, derivativesum
    real(8), allocatable :: qTotal(:)
    real(8) :: distance, cubeDistance
    real(8) :: deltaOrigin(3)
    real(8) :: cosmoEpsilon
    integer :: AA
    real(8) :: dax, day, daz
    real(8) :: constantDer
    integer :: stat
    integer :: icharges


    wfnFile = "lowdin.wfn"
    wfnUnit = 20


    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if(allocated(auxOwnerId)) deallocate(auxOwnerId)
    allocate(auxOwnerId(size(ParticleManager_instance)))

    auxOwnerId = 0

    j = 1
    do i = 1, size(ParticleManager_instance)
       if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
          auxOwnerId(i) = j
          j = j + 1
       else
          auxOwnerId(i) = 0
       end if
    end do

    if(allocated(EnergyGradients_instance%gradients%kinetic)) deallocate(EnergyGradients_instance%gradients%kinetic)
    allocate(EnergyGradients_instance%gradients%kinetic(numberOfOptimizationCenters*3))

    if(allocated(EnergyGradients_instance%gradients%attraction)) deallocate(EnergyGradients_instance%gradients%attraction)
    allocate(EnergyGradients_instance%gradients%attraction(numberOfOptimizationCenters*3))

    if(allocated(EnergyGradients_instance%gradients%overlap)) deallocate(EnergyGradients_instance%gradients%overlap)
    allocate(EnergyGradients_instance%gradients%overlap(numberOfOptimizationCenters*3))

    if(allocated(EnergyGradients_instance%gradients%coulomb)) deallocate(EnergyGradients_instance%gradients%coulomb)
    allocate(EnergyGradients_instance%gradients%coulomb(numberOfOptimizationCenters*3))

    if(allocated(EnergyGradients_instance%gradients%exchange)) deallocate(EnergyGradients_instance%gradients%exchange)
    allocate(EnergyGradients_instance%gradients%exchange(numberOfOptimizationCenters*3))

    ! if(allocated(EnergyGradients_instance%gradients%total)) deallocate(EnergyGradients_instance%gradients%total)
    ! allocate(EnergyGradients_instance%gradients%total(numberOfOptimizationCenters*3))

    ! EnergyGradients_instance%gradients%total = 0.0_8
    EnergyGradients_instance%gradients%kinetic = 0.0_8
    EnergyGradients_instance%gradients%attraction = 0.0_8
    EnergyGradients_instance%gradients%overlap = 0.0_8
    EnergyGradients_instance%gradients%coulomb = 0.0_8
    EnergyGradients_instance%gradients%exchange = 0.0_8


    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    do specieIterator=1, MolecularSystem_instance%numberOfQuantumSpecies

       if(allocated(auxKinetic)) deallocate(auxKinetic)
       allocate(auxKinetic(numberOfOptimizationCenters,3))

       if(allocated(auxPotential)) deallocate(auxPotential)
       allocate(auxPotential(numberOfOptimizationCenters,3))

       if(CONTROL_instance%COSMO) then
          if(allocated(auxCOSMO)) deallocate(auxCOSMO)
          allocate(auxCOSMO(numberOfOptimizationCenters,3))
          if(allocated(auxCOSMO2)) deallocate(auxCOSMO2)
          allocate(auxCOSMO2(numberOfOptimizationCenters,3))
          if(allocated(auxCOSMO3)) deallocate(auxCOSMO3)
          allocate(auxCOSMO3(surface%sizeSurface,numberOfOptimizationCenters,3))
          if(allocated(qTotal)) deallocate(qTotal)
          allocate(qTotal(surface%sizeSurface))
       end if

       if(allocated(auxOverlap)) deallocate(auxOverlap)
       allocate(auxOverlap(numberOfOptimizationCenters,3))

       if(allocated(auxCoulomb)) deallocate(auxCoulomb)
       allocate(auxCoulomb(numberOfOptimizationCenters,3))

       if(allocated(auxExchange)) deallocate(auxExchange)
       allocate(auxExchange(numberOfOptimizationCenters,3))

       auxKinetic = 0.0_8
       auxPotential = 0.0_8
       auxOverlap = 0.0_8
       auxCoulomb = 0.0_8
       auxExchange = 0.0_8       

       call MolecularSystem_getBasisSet(specieIterator, contractions)

       numberOfContractions = MolecularSystem_getNumberOfContractions(specieIterator)

       if(allocated(labelsOfContractions)) deallocate(labelsOfContractions)
       allocate(labelsOfContractions(numberOfContractions))

       j = 1
       do i = 1, numberOfContractions
          labelsOfContractions(i) = j
          j = j + contractions(i)%numCartesianOrbital
       end do

       nameOfSpecie = trim(MolecularSystem_instance%species(specieIterator)%symbol)
       orderOfMatrix = MolecularSystem_getTotalNumberOfContractions(specieIterator)

       ! write(*,"(A,A,A,I)")"Especie: ", trim(nameOfSpecie), " size: ", numberOfContractions
       ! write(*,"(A,I)") "Order: ", orderOfMatrix
       ! do i = 1, numberOfContractions
       !    write(*,"(A,I)") "ID: ", contractions(i)%id
       !    write(*,"(A,I)") "Length: ", contractions(i)%length
       !    write(*,"(A,I)") "am: ", contractions(i)%angularMoment
       !    write(*,"(A,I)") "nc: ", contractions(i)%numCartesianOrbital
       !    write(*,"(A,I)") "owner: ", contractions(i)%owner
       ! end do
       ocupationNumber = MolecularSystem_getOcupationNumber( specieIterator )

       arguments(2) = MolecularSystem_getNameOfSpecie(specieIterator)

       arguments(1) = "DENSITY"
       densityMatrix = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       ! write(*,*) "Matriz densidad"
       ! call Matrix_show(densityMatrix)
       ! densityMatrix%values = transpose(densityMatrix%values) 

       arguments(1) = "COEFFICIENTS"
       matrixOfEigenvectors = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       ! arguments(1) = "COUPLING"
       ! twoParticles = &
       !      Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
       !      columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = int(orderOfMatrix,4), &
            unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
            output = vectorOfEigenvalues )

       call Matrix_constructor(auxWeightDensity, int(orderOfMatrix,8), int(ocupationNumber,8))

       matrixSum = 0.0_8
       do i=1, orderOfMatrix
          do j=1, orderOfMatrix
             matrixSum = matrixSum + densityMatrix%values(i,j)
          end do
       end do
       ! write(*,"(A,F17.12)") "Suma de la matriz densidad: ", matrixSum



       do i=1, orderOfMatrix
          do j=1, ocupationNumber
             auxWeightDensity%values(i,j) = matrixOfEigenvectors%values(i,j)*vectorOfEigenvalues%values(j)
          end do
       end do

       call Matrix_constructor(weightDensityMatrix, int(orderOfMatrix,8), int(orderOfMatrix,8))

       !! OJO esto solo es para RHF, para UHF hay que separar los orbitales alpha y beta
       weightDensityMatrix = Matrix_multiplication("N","T",&
            orderOfMatrix,&
            orderOfMatrix,&
            ocupationNumber,&
            2.0_8,&
            matrixOfEigenvectors,&
            orderOfMatrix,&
            auxWeightDensity,&
            orderOfMatrix,&
            1.0_8,&
            orderOfMatrix)

       mass = 0.0_8
       mass = MolecularSystem_getMass(specieIterator)
       charge = MolecularSystem_getCharge(specieIterator)
       eta = MolecularSystem_getEta(specieIterator)
       lambda = MolecularSystem_getLambda(specieIterator)
       kappa = MolecularSystem_getKappa(specieIterator)


       ! write(*,"(A,F12.8)") "lambda: ", lambda
       ! write(*,"(A,F12.8)") "kappa: ", kappa
       ! write(*,"(A,F12.8)") "eta: ", eta

       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes Cineticos"
       ! write(*,"(A)") "----------------------------------------------------------------"
       ! do i=1, numberOfContractions
       !    write(*,"(A,I)") "owner: ", owner
       ! end do
       ! Kinetic Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerP = auxOwnerId(contractions(P)%owner) !! aP
             centerQ = auxOwnerId(contractions(Q)%owner) !! aQ

             if (P == Q) then
                perm = 1.0
             else
                perm = 2.0
             end if

             perm = (perm/mass)



             ! write(*,"(A1,I1,A1,I1,A1)") "(", P, "|", Q, ")"
             ! write(*,"(A,I1,A,I1)") "CenterP: ", centerP, " centerQ: ", centerQ
             ! write(*,"(A,f12.8)") "Scal: ", 1.0/mass

             call DerivativeManager_getElement( KINETIC_DERIVATIVES, auxVector, i=P, j=Q, nameOfSpecie=nameOfSpecie )

             !write(*,*) 'Derivadas Cineticas'
             ! !write(*,*) auxVector(:)

             i = 0
             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerP,1) = auxKinetic(centerP,1) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Px: ", auxKinetic(centerP,1), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerP,2) = auxKinetic(centerP,2) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Py: ", auxKinetic(centerP,2), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerP,3) = auxKinetic(centerP,3) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Pz: ", auxKinetic(centerP,3), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerQ,1) = auxKinetic(centerQ,1) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Qx: ", auxKinetic(centerQ,1), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerQ,2) = auxKinetic(centerQ,2) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Qy: ", auxKinetic(centerQ,2), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerQ,3) = auxKinetic(centerQ,3) + perm*densityMatrix%values(u,v)*auxVector(i)
                   !write(*,"(A,3f17.12)") "Kinetic Vector Qz: ", auxKinetic(centerQ,3), auxVector(i), densityMatrix%values(u,v)
                   i = i + 1
                end do
             end do
             ! write(*,"(A)") "----------------------------------------------------------------"
             ! write(*,"(A,3f17.12)") "Kinetic Vector final: ", auxKinetic(centerQ,:)
             ! write(*,"(A)") "----------------------------------------------------------------"
          end do
       end do

       ! write(*,"(3(f17.12))") auxKinetic(1,1), auxKinetic(1,2), auxKinetic(1,3)
       ! write(*,"(3(f17.12))") auxKinetic(2,1), auxKinetic(2,2), auxKinetic(2,3)
       ! write(*,"(3(f17.12))") auxKinetic(3,1), auxKinetic(3,2), auxKinetic(3,3)
       ! write(*,"(A)") "----------------------------------------------------------------"


       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes de Potencial"
       ! write(*,"(A)") "----------------------------------------------------------------"
       !! Potential Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerP = auxOwnerId(contractions(P)%owner) !! aP
             centerQ = auxOwnerId(contractions(Q)%owner) !! aQ

             if (P == Q) then
                perm = 1.0
             else
                perm = 2.0
             end if

             ! write(*,"(A1,I1,A1,I1,A1)") "(", P, "|", Q, ")"
             ! write(*,"(A,I,A,I)") "Centro P: ", centerP, " owner: ", contractions(P)%owner
             ! write(*,"(A,I,A,I)") "Centro Q: ", centerQ, " owner: ", contractions(Q)%owner
             call DerivativeManager_getElement( ATTRACTION_DERIVATIVES, &
                  auxVector2, i=P, j=Q, nameOfSpecie=nameOfSpecie, A=centerP, B=centerQ )


             ! write(*,"(A)") "---------------------------------------------------------------------------------------"
             ! write(*,*) "Derivadas Potencial: "
             ! write(*,*) P, Q, " | ", auxVector2(:)

             i = 0
             j = 0
             k = 0
             ! center = 1
             ! A = 1
             ! do l = 1, size(ParticleManager_instance)
             !    if(ParticleManager_instance(l)%particlePtr%isCenterOfOptimization) then
             !       if(ParticleManager_instance(l)%particlePtr%isQuantum) then
             !          auxPotential(center,1) = 0.0_8
             !          auxPotential(center,2) = 0.0_8
             !          auxPotential(center,3) = 0.0_8
             !          center = center + 1
             !          A = A + 1
             !       else
             do A=1, numberOfOptimizationCenters
                i = 3*(A-1)*numCartesianP*numCartesianQ + 0*numCartesianP*numCartesianQ
                j = 3*(A-1)*numCartesianP*numCartesianQ + 1*numCartesianP*numCartesianQ
                k = 3*(A-1)*numCartesianP*numCartesianQ + 2*numCartesianP*numCartesianQ
                ! write(*,"(A,I2,I2,3I3)") "center: ", center, A, i, j, k
                do pIter = 0, numCartesianP - 1
                   do qIter = 0, numCartesianQ - 1
                      u = pIter + labelsOfContractions(P)
                      v = qIter + labelsOfContractions(Q)

                      Vval = perm*densityMatrix%values(u,v)
                      ! write(*,*) "Label P or ", labelsOfContractions(P)
                      ! write(*,*) "Label Q or ", labelsOfContractions(Q)
                      ! write(*,*) "(u|v) ", u, " | ", v, " : ", densityMatrix%values(u,v)
                      ! write(*,*) "(v,u) ", v, " | ", u, " : ", densityMatrix%values(v,u)

                      auxPotential(A,1) = auxPotential(A,1) + Vval*auxVector2(i)
                      !write(*,"(A,3f17.12)") "Potencial Vector x: ", auxPotential(center,1), Vval, auxVector2(i)

                      i = i + 1
                      auxPotential(A,2) = auxPotential(A,2) + Vval*auxVector2(j)
                      !write(*,"(A,3f17.12)") "Potencial Vector y: ", auxPotential(center,2), Vval, auxVector2(j)

                      j = j + 1
                      auxPotential(A,3) = auxPotential(A,3) + Vval*auxVector2(k)
                      !write(*,"(A,3f17.12)") "Potencial Vector z: ", auxPotential(center,3), Vval, auxVector2(k)

                      k = k + 1

                   end do
                end do
                ! A = A + 1
                ! center = center + 1
                !    end if
                ! end if
             end do
             ! write(*,"(A)") "================================================================"
             ! write(*,"(3(f17.12))") auxPotential(1,1), auxPotential(1,2), auxPotential(1,3)
             ! write(*,"(A)") "================================================================"
          end do
       end do

       ! write(*,"(3(f17.12))") auxPotential(1,1), auxPotential(1,2), auxPotential(1,3)
       ! write(*,"(3(f17.12))") auxPotential(2,1), auxPotential(2,2), auxPotential(2,3)
       ! write(*,"(3(f17.12))") auxPotential(3,1), auxPotential(3,2), auxPotential(3,3)
       ! write(*,"(A)") "----------------------------------------------------------------"

       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes de Overlap"
       ! write(*,"(A)") "----------------------------------------------------------------"
       !! Overlap Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             call DerivativeManager_getElement( OVERLAP_DERIVATIVES, auxVector3, i=P, j=Q, nameOfSpecie=nameOfSpecie )

             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerP = auxOwnerId(contractions(P)%owner) !! aP
             centerQ = auxOwnerId(contractions(Q)%owner) !! aQ

             if (P == Q) then
                perm = 1.0
             else
                perm = 2.0
             end if

             i = 0
             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerP,1) = auxOverlap(centerP,1) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerP,2) = auxOverlap(centerP,2) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerP,3) = auxOverlap(centerP,3) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerQ,1) = auxOverlap(centerQ,1) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerQ,2) = auxOverlap(centerQ,2) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerQ,3) = auxOverlap(centerQ,3) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

          end do
       end do

       ! write(*,"(3(f17.12))") auxOverlap(1,1), auxOverlap(1,2), auxOverlap(1,3)
       ! write(*,"(3(f17.12))") auxOverlap(2,1), auxOverlap(2,2), auxOverlap(2,3)
       ! write(*,"(3(f17.12))") auxOverlap(3,1), auxOverlap(3,2), auxOverlap(3,3)
       ! write(*,"(A)") "----------------------------------------------------------------"

       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes de Coulomb"
       ! write(*,"(A)") "----------------------------------------------------------------"
       !! Coulomb and Exchange Gradients
       call EnergyGradients_getShellPairs(numberOfContractions, shellPairs)
       npairs = size(shellPairs%values,DIM=1)
       npairs2 = npairs*npairs

       deltasum = 0
       do pqrsIter = 0, npairs2 - 1
          PQ = pqrsIter/npairs
          RS = mod(pqrsIter,npairs)

          if (RS.gt.PQ) cycle

          P = shellPairs%values(PQ+1,1)
          Q = shellPairs%values(PQ+1,2)
          R = shellPairs%values(RS+1,1)
          S = shellPairs%values(RS+1,2)
          ! write(*,"(A)") "Contraida:"
          ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",P,",",Q,"|",R,",",S,")"
          ! write(*,"(A)") "Primitivas:"
          call DerivativeManager_getElement( REPULSION_DERIVATIVES, auxVector4, i=P, j=Q, k=R, l=S, nameOfSpecie=nameOfSpecie )

          numCartesianP = contractions(P)%numCartesianOrbital
          numCartesianQ = contractions(Q)%numCartesianOrbital
          numCartesianR = contractions(R)%numCartesianOrbital
          numCartesianS = contractions(S)%numCartesianOrbital

          centerP = auxOwnerId(contractions(P)%owner)
          centerQ = auxOwnerId(contractions(Q)%owner)
          centerR = auxOwnerId(contractions(R)%owner)
          centerS = auxOwnerId(contractions(S)%owner)

          perm = 1.0_8
          if (P /= Q) then
             perm = perm*2.0_8
          end if
          if (R /= S) then
             perm = perm*2.0_8
          end if
          if (PQ /= RS) then
             perm = perm*2.0_8
          end if

          deltasum = deltasum + perm
          perm = perm*0.5_8

          stride = numCartesianP*numCartesianQ*numCartesianR*numCartesianS

          Ax = 0.0_8
          Ay = 0.0_8
          Az = 0.0_8
          Bx = 0.0_8
          By = 0.0_8
          Bz = 0.0_8
          Cx = 0.0_8
          Cy = 0.0_8
          Cz = 0.0_8
          Dx = 0.0_8
          Dy = 0.0_8
          Dz = 0.0_8
          delta = 0


          ! derivativesum = 0.0_8
          ! do i=1, stride*9
          !    derivativesum = derivativesum + auxVector4(i)
          ! end do

          ! write(*,"(A,F17.12)") "Suma derivadas: ", derivativesum

          do pIter=0, numCartesianP-1
             do qIter=0, numCartesianQ-1
                do rIter=0, numCartesianR-1
                   do sIter=0, numCartesianS-1
                      u = pIter + labelsOfContractions(P)
                      v = qIter + labelsOfContractions(Q)
                      x = rIter + labelsOfContractions(R)
                      y = sIter + labelsOfContractions(S)
                      Duv = densityMatrix%values(u,v)
                      Dxy = densityMatrix%values(x,y)
                      JKval = perm*Duv*Dxy*2.0_8!eta
                      !write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",pIter,",",qIter,"|",rIter,",",sIter,")"
                      !write(*,"(A3,I2,A4,I2,A4,I2,A4,I2)") "u: ", u, " v: ", v," x: ", x," y: ", y
                      !write(*,"(2(A,F17.12))") "Duv: ", Duv, " Dxy: ", Dxy
                      Ax = Ax + JKval * auxVector4(0 * stride + delta)
                      !write(*,"(A,2F17.12)") "Ax: ", Ax, auxVector4(0 * stride + delta)
                      Ay = Ay + JKval * auxVector4(1 * stride + delta)
                      !write(*,"(A,2F17.12)") "Ay: ", Ay, auxVector4(1 * stride + delta)
                      Az = Az + JKval * auxVector4(2 * stride + delta)
                      !write(*,"(A,2F17.12)") "Az: ", Az, auxVector4(2 * stride + delta)
                      Cx = Cx + JKval * auxVector4(3 * stride + delta)
                      !write(*,"(A,2F17.12)") "Cx: ", Cx, auxVector4(3 * stride + delta)
                      Cy = Cy + JKval * auxVector4(4 * stride + delta)
                      !write(*,"(A,2F17.12)") "Cy: ", Cy, auxVector4(4 * stride + delta)
                      Cz = Cz + JKval * auxVector4(5 * stride + delta)
                      !write(*,"(A,2F17.12)") "Cz: ", Cz, auxVector4(5 * stride + delta)
                      Dx = Dx + JKval * auxVector4(6 * stride + delta)
                      !write(*,"(A,2F17.12)") "Dx: ", Dx, auxVector4(6 * stride + delta)
                      Dy = Dy + JKval * auxVector4(7 * stride + delta)
                      !write(*,"(A,2F17.12)") "Dy: ", Dy, auxVector4(7 * stride + delta)
                      Dz = Dz + JKval * auxVector4(8 * stride + delta)
                      !write(*,"(A,2F17.12)") "Dz: ", Dz, auxVector4(8 * stride + delta)
                      !write(*,"(A)") "---------------------------"
                      delta = delta + 1
                   end do
                end do
             end do
          end do
          Bx = -(Ax + Cx + Dx)
          By = -(Ay + Cy + Dy)
          Bz = -(Az + Cz + Dz)

          auxCoulomb(centerP,1) = auxCoulomb(centerP,1) + Ax
          auxCoulomb(centerP,2) = auxCoulomb(centerP,2) + Ay
          auxCoulomb(centerP,3) = auxCoulomb(centerP,3) + Az
          auxCoulomb(centerQ,1) = auxCoulomb(centerQ,1) + Bx
          auxCoulomb(centerQ,2) = auxCoulomb(centerQ,2) + By
          auxCoulomb(centerQ,3) = auxCoulomb(centerQ,3) + Bz
          auxCoulomb(centerR,1) = auxCoulomb(centerR,1) + Cx
          auxCoulomb(centerR,2) = auxCoulomb(centerR,2) + Cy
          auxCoulomb(centerR,3) = auxCoulomb(centerR,3) + Cz
          auxCoulomb(centerS,1) = auxCoulomb(centerS,1) + Dx
          auxCoulomb(centerS,2) = auxCoulomb(centerS,2) + Dy
          auxCoulomb(centerS,3) = auxCoulomb(centerS,3) + Dz

          ! write(*,"(A,F17.12)") "valor: ", auxCoulomb(1,3)

          ! write(*,"(A,A4,3F17.12)") "A: ", " -> ", Ax, Ay, Az
          ! write(*,"(A,A4,3F17.12)") "B: ", " -> ", Bx, By, Bz
          ! write(*,"(A,A4,3F17.12)") "C: ", " -> ", Cx, Cy, Cz
          ! write(*,"(A,A4,3F17.12)") "D: ", " -> ", Dx, Dy, Dz
          ! write(*,"(A)") "---------------------------"



          ! deltasum = deltasum + delta
          ! write(*,"(A)") "---------------------------"
          ! write(*,"(A,I1,A4,3F17.12)") "Centro: ", centerP, " -> ", auxCoulomb(centerP,:)
          ! write(*,"(A,I1,A4,3F17.12)") "Centro: ", centerQ, " -> ", auxCoulomb(centerQ,:)
          ! write(*,"(A,I1,A4,3F17.12)") "Centro: ", centerR, " -> ", auxCoulomb(centerR,:)
          ! write(*,"(A,I1,A4,3F17.12)") "Centro: ", centerS, " -> ", auxCoulomb(centerS,:)
          ! write(*,"(A,I)") "delta: ", delta
          ! write(*,"(A)") "---------------------------"

          Ax = 0.0_8
          Ay = 0.0_8
          Az = 0.0_8
          Bx = 0.0_8
          By = 0.0_8
          Bz = 0.0_8
          Cx = 0.0_8
          Cy = 0.0_8
          Cz = 0.0_8
          Dx = 0.0_8
          Dy = 0.0_8
          Dz = 0.0_8
          delta = 0
          do pIter=0, numCartesianP-1
             do qIter=0, numCartesianQ-1
                do rIter=0, numCartesianR-1
                   do sIter=0, numCartesianS-1
                      JKval = 0.0_8
                      u = pIter + labelsOfContractions(P)
                      v = qIter + labelsOfContractions(Q)
                      x = rIter + labelsOfContractions(R)
                      y = sIter + labelsOfContractions(S)
                      Duv = densityMatrix%values(u,x)*0.5_8
                      Dxy = densityMatrix%values(v,y)*0.5_8
                      JKval = JKval + perm*Duv*Dxy
                      Duv = densityMatrix%values(u,y)*0.5_8
                      Dxy = densityMatrix%values(v,x)*0.5_8
                      JKval = JKval + perm*Duv*Dxy
                      Duv = densityMatrix%values(u,x)*0.5_8
                      Dxy = densityMatrix%values(v,y)*0.5_8
                      JKval = JKval + perm*Duv*Dxy
                      Duv = densityMatrix%values(u,y)*0.5_8
                      Dxy = densityMatrix%values(v,x)*0.5_8
                      JKval = JKval + perm*Duv*Dxy
                      Ax = Ax + JKval * auxVector4(0 * stride + delta)
                      Ay = Ay + JKval * auxVector4(1 * stride + delta)
                      Az = Az + JKval * auxVector4(2 * stride + delta)
                      Cx = Cx + JKval * auxVector4(3 * stride + delta)
                      Cy = Cy + JKval * auxVector4(4 * stride + delta)
                      Cz = Cz + JKval * auxVector4(5 * stride + delta)
                      Dx = Dx + JKval * auxVector4(6 * stride + delta)
                      Dy = Dy + JKval * auxVector4(7 * stride + delta)
                      Dz = Dz + JKval * auxVector4(8 * stride + delta)
                      delta = delta + 1
                   end do
                end do
             end do
          end do
          Bx = -(Ax + Cx + Dx)
          By = -(Ay + Cy + Dy)
          Bz = -(Az + Cz + Dz)

          auxExchange(centerP,1) = auxExchange(centerP,1) + Ax
          auxExchange(centerP,2) = auxExchange(centerP,2) + Ay
          auxExchange(centerP,3) = auxExchange(centerP,3) + Az
          auxExchange(centerQ,1) = auxExchange(centerQ,1) + Bx
          auxExchange(centerQ,2) = auxExchange(centerQ,2) + By
          auxExchange(centerQ,3) = auxExchange(centerQ,3) + Bz
          auxExchange(centerR,1) = auxExchange(centerR,1) + Cx
          auxExchange(centerR,2) = auxExchange(centerR,2) + Cy
          auxExchange(centerR,3) = auxExchange(centerR,3) + Cz
          auxExchange(centerS,1) = auxExchange(centerS,1) + Dx
          auxExchange(centerS,2) = auxExchange(centerS,2) + Dy
          auxExchange(centerS,3) = auxExchange(centerS,3) + Dz

       end do

       auxCoulomb = auxCoulomb*charge*charge*0.5
       auxExchange = auxExchange*charge*charge*kappa*0.5

       ! write(*,"(3(f17.12))") auxCoulomb(1,1), auxCoulomb(1,2), auxCoulomb(1,3)
       ! write(*,"(3(f17.12))") auxCoulomb(2,1), auxCoulomb(2,2), auxCoulomb(2,3)
       ! write(*,"(3(f17.12))") auxCoulomb(3,1), auxCoulomb(3,2), auxCoulomb(3,3)
       ! write(*,"(A)") "----------------------------------------------------------------"

       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes de Intercambio"
       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(3(f17.12))") auxExchange(1,1), auxExchange(1,2), auxExchange(1,3)
       ! write(*,"(3(f17.12))") auxExchange(2,1), auxExchange(2,2), auxExchange(2,3)
       ! write(*,"(3(f17.12))") auxExchange(3,1), auxExchange(3,2), auxExchange(3,3)
       ! write(*,"(A)") "----------------------------------------------------------------"

       ! write(*,"(A,I)") "delta: ", deltasum
       ! write(*,"(A)") "----------------------------------------------------------------"
       ! write(*,"(A)") " Gradientes de Potencial COSMO"
       ! write(*,"(A)") "----------------------------------------------------------------"
       !! Potential Gradients
       if(CONTROL_instance%COSMO) then
          cosmoEpsilon=(CONTROL_instance%COSMO_SOLVENT_DIELECTRIC+CONTROL_instance%COSMO_SCALING)/(CONTROL_instance%COSMO_SOLVENT_DIELECTRIC-1)

          open(unit=77, file="qTotalCosmo.charges", status="unknown",form="unformatted")
          read(77)(qTotal(i),i=1,surface%sizeSurface)

          close(unit=77)

          write(*,"(A)") "---------------------------------------------------------------------------------------"
          write(*,*) "Derivadas: "

          do P = 1, numberOfContractions
             do Q = 1, P
                numCartesianP = contractions(P)%numCartesianOrbital  !! nP
                numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
                centerP = auxOwnerId(contractions(P)%owner) !! aP
                centerQ = auxOwnerId(contractions(Q)%owner) !! aQ

                if (P == Q) then
                   perm = 1.0
                else
                   perm = 2.0
                end if

                ! write(*,"(A1,I1,A1,I1,A1)") "(", P, "|", Q, ")"
                ! write(*,"(A,I,A,I)") "Centro P: ", centerP, " owner: ", contractions(P)%owner
                ! write(*,"(A,I,A,I)") "Centro Q: ", centerQ, " owner: ", contractions(Q)%owner
                call DerivativeManager_getElement( ATTRACTION_DERIVATIVES, &
                     auxVector2, surface, i=P, j=Q, nameOfSpecie=nameOfSpecie, A=centerP, B=centerQ )

                ! write(*,*) P, Q, " | ", auxVector2(:)

                i = 0
                j = 0
                k = 0
                ! center = 1


                do A=1, surface%sizeSurface
                   i = 3*(A-1)*numCartesianP*numCartesianQ + 0*numCartesianP*numCartesianQ
                   j = 3*(A-1)*numCartesianP*numCartesianQ + 1*numCartesianP*numCartesianQ
                   k = 3*(A-1)*numCartesianP*numCartesianQ + 2*numCartesianP*numCartesianQ
                   ! write(*,"(A,I2,I2,3I3)") "center: ", center, A, i, j, k
                   do pIter = 0, numCartesianP - 1
                      do qIter = 0, numCartesianQ - 1
                         u = pIter + labelsOfContractions(P)
                         v = qIter + labelsOfContractions(Q)
                         Vval = perm*densityMatrix%values(u,v)
                         ! write(*,*) "Label P", labelsOfContractions(P), " Center P: ", centerP
                         ! write(*,*) "Label Q", labelsOfContractions(Q), " Center Q: ", centerQ
                         ! write(*,*) "(u|v) ", u, " | ", v, " : ", densityMatrix%values(u,v)
                         ! write(*,*) "(v,u) ", v, " | ", u, " : ", densityMatrix%values(v,u)

                         do center=1, numberOfOptimizationCenters
                            ! if((centerP.EQ.center).and.(surface%atoms(A).NE.center))then
                            !    auxCOSMO(center,1) = auxCOSMO(center,1) + Vval*auxVector2(i)!*qTotal(A)
                            !    i=i+1
                            !    auxCOSMO(center,2) = auxCOSMO(center,2) + Vval*auxVector2(j)!*qTotal(A)
                            !    j=j+1
                            !    auxCOSMO(center,3) = auxCOSMO(center,3) + Vval*auxVector2(k)!*qTotal(A)
                            !    k=k+1
                            ! else if((centerP.NE.center).and.(surface%atoms(A).EQ.center))then
                            !    auxCOSMO(center,1) = auxCOSMO(center,1) - Vval*auxVector2(i)!*qTotal(A)
                            !    i=i+1
                            !    auxCOSMO(center,2) = auxCOSMO(center,2) - Vval*auxVector2(j)!*qTotal(A)
                            !    j=j+1
                            !    auxCOSMO(center,3) = auxCOSMO(center,3) - Vval*auxVector2(k)!*qTotal(A)
                            !    k=k+1
                            ! end if
                            auxCOSMO(center,1) = auxCOSMO(center,1) + Vval*auxVector2(i)!*qTotal(A)
                            i=i+1
                            auxCOSMO(center,2) = auxCOSMO(center,2) + Vval*auxVector2(j)!*qTotal(A)
                            j=j+1
                            auxCOSMO(center,3) = auxCOSMO(center,3) + Vval*auxVector2(k)!*qTotal(A)
                            k=k+1

                         end do

                      end do

                   end do
                   do center=1, numberOfOptimizationCenters
                      do m=1, 3
                         auxCOSMO3(A,center,m)= auxCOSMO3(A,center,m) + auxCOSMO(center,m) - auxCOSMO2(center,m)
                      end do
                   end do
                   auxCOSMO2=auxCOSMO

                   ! A = A + 1
                   !    end if
                   ! end if
                end do
                ! write(*,"(A)") "================================================================"
                ! write(*,"(3(f17.12))") auxPotential(1,1), auxPotential(1,2), auxPotential(1,3)
                ! write(*,"(A)") "================================================================"

             end do
          end do
          do center=1, numberOfOptimizationCenters
             do A=1, surface%sizeSurface

                write(*,*) auxCOSMO3(A,center,:), center, A

             end do
          end do
          write(*,*)"1X",auxCOSMO(1,1) 
          write(*,*)"2X",auxCOSMO(2,1) 
          write(*,*)"3X",auxCOSMO(3,1) 

          deltaOrigin=0.0_8
          do A=1, surface%sizeSurface
             do center=1, numberOfOptimizationCenters
                Zi = ParticleManager_instance(center)%particlePtr%charge
                deltaOrigin(1) = (ParticleManager_instance(center)%particlePtr%origin(1) - surface%xs(A)) 
                deltaOrigin(2) = (ParticleManager_instance(center)%particlePtr%origin(2) - surface%ys(A)) 
                deltaOrigin(3) = (ParticleManager_instance(center)%particlePtr%origin(3) - surface%zs(A)) 
                distance = sqrt( sum( deltaOrigin**2.0_8 ) )
                cubeDistance = distance*distance*distance
                do centerP=1, numberOfOptimizationCenters

                   if((centerP.EQ.center).and.(surface%atoms(A).NE.center))then
                      auxCOSMO(center,1) = auxCOSMO(center,1)  - Zi*qTotal(A)*deltaOrigin(1)/cubeDistance
                      auxCOSMO(center,2) = auxCOSMO(center,2)  - Zi*qTotal(A)*deltaOrigin(2)/cubeDistance
                      auxCOSMO(center,3) = auxCOSMO(center,3)  - Zi*qTotal(A)*deltaOrigin(3)/cubeDistance
                   else if((centerP.NE.center).and.(surface%atoms(A).EQ.center))then
                      auxCOSMO(center,1) = auxCOSMO(center,1)  + Zi*qTotal(A)*deltaOrigin(1)/cubeDistance
                      auxCOSMO(center,2) = auxCOSMO(center,2)  + Zi*qTotal(A)*deltaOrigin(2)/cubeDistance
                      auxCOSMO(center,3) = auxCOSMO(center,3)  + Zi*qTotal(A)*deltaOrigin(3)/cubeDistance
                   end if
                end do
             end do
          end do
          deltaOrigin=0.0_8
	  constantDer=-0.5_8*cosmoEpsilon*1.07_8*Math_SQRT_PI
          open(unit=78, file=trim(CONTROL_instance%INPUT_FILE)//"der", status="old")
	  read(78,*) AA, centerP, dax, day, daz
   ! write(*,*) "lectura de .der ", AA, centerP, dax, day, daz  	
          do A=1, surface%sizeSurface
             do B=1, surface%sizeSurface
                deltaOrigin(1) = surface%xs(A) - surface%xs(B) 
                deltaOrigin(2) = surface%ys(A) - surface%ys(B) 
                deltaOrigin(3) = surface%zs(A) - surface%zs(B) 
                distance = sqrt( sum( deltaOrigin**2.0_8 ) )
                cubeDistance = distance*distance*distance
                do center=1, numberOfOptimizationCenters
                   if(A.EQ.B)then
                      if (AA.EQ.A) then
                         if (centerP.EQ.center) then


                            auxCOSMO(center,1) = auxCOSMO(center,1)  +constantDer*qTotal(A)*qTotal(A)*dax/sqrt(surface%area(A)*surface%area(A)*surface%area(A))
                            auxCOSMO(center,2) = auxCOSMO(center,2)  +constantDer*qTotal(A)*qTotal(A)*day/sqrt(surface%area(A)*surface%area(A)*surface%area(A))
                            auxCOSMO(center,3) = auxCOSMO(center,3)  +constantDer*qTotal(A)*qTotal(A)*daz/sqrt(surface%area(A)*surface%area(A)*surface%area(A))
                            read(78,*,iostat=stat) AA, centerP, dax, day, daz
                            ! write(*,*) "lectura de .der ", AA, centerP, dax, day, daz  	
                         end if
                      end if
                   else if ((surface%atoms(A).EQ.center).and.(surface%atoms(B).NE.center)) then 
                      auxCOSMO(center,1) = auxCOSMO(center,1)  -0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(1)/cubeDistance
                      auxCOSMO(center,2) = auxCOSMO(center,2)  -0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(2)/cubeDistance
                      auxCOSMO(center,3) = auxCOSMO(center,3)  -0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(3)/cubeDistance
                   else if ((surface%atoms(B).EQ.center).and.(surface%atoms(A).NE.center)) then 
                      auxCOSMO(center,1) = auxCOSMO(center,1)  +0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(1)/cubeDistance
                      auxCOSMO(center,2) = auxCOSMO(center,2)  +0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(2)/cubeDistance
                      auxCOSMO(center,3) = auxCOSMO(center,3)  +0.5_8*(cosmoEpsilon)*qTotal(A)*qTotal(B)*deltaOrigin(3)/cubeDistance
                   end if
                end do
             end do
          end do
          close(78)
          write(*,*)"gradiente cosmo"
          write(*,*)auxCOSMO(:,:)	
       end if


       k=1
       do i=1, numberOfOptimizationCenters
          do j=1,3
             EnergyGradients_instance%gradients%total(k) = EnergyGradients_instance%gradients%total(k) + auxKinetic(i,j) + auxPotential(i,j) + auxOverlap(i,j) + auxCoulomb(i,j) + auxExchange(i,j)
             EnergyGradients_instance%gradients%kinetic(k) = EnergyGradients_instance%gradients%kinetic(k) + auxKinetic(i,j)
             EnergyGradients_instance%gradients%attraction(k) = EnergyGradients_instance%gradients%attraction(k) + auxPotential(i,j)
             EnergyGradients_instance%gradients%overlap(k) = EnergyGradients_instance%gradients%overlap(k) + auxOverlap(i,j)
             EnergyGradients_instance%gradients%coulomb(k) = EnergyGradients_instance%gradients%coulomb(k) + auxCoulomb(i,j)
             EnergyGradients_instance%gradients%exchange(k) = EnergyGradients_instance%gradients%Exchange(k) + auxExchange(i,j)
             k = k + 1
          end do
       end do

       deallocate(auxKinetic) 
       deallocate(auxPotential)
       deallocate(auxOverlap)
       deallocate(auxCoulomb)
       deallocate(auxExchange)
       deallocate(auxVector) 
       deallocate(auxVector2)
       deallocate(auxVector3)
       deallocate(auxVector4)
       deallocate(labelsOfContractions)

       call Matrix_destructor(matrixOfEigenvectors)
       call Matrix_destructor(densityMatrix)
       call Vector_destructor(vectorOfEigenvalues)
       call Matrix_destructor(auxWeightDensity)
       call Matrix_destructor(weightDensityMatrix)

    end do

    deallocate(auxOwnerId)

    close(wfnUnit)

  end subroutine EnergyGradients_calculateAnalyticUncoupledFirstDerivative

  !**
  ! @brief  Calcula la primera derivada analitica para la matriz de Fock sin incluir el termino de
  !               acoplamiento interspecie
  !**
  subroutine EnergyGradients_calculateAnalyticCouplingFirstDerivative()
    implicit none
    real(8), allocatable :: output(:,:)
    real(8), allocatable :: auxoutput(:,:)
    integer :: specieIterator, otherSpecieIterator
    type(ContractedGaussian), allocatable :: contractions(:)
    type(ContractedGaussian), allocatable :: otherContractions(:)
    integer :: numberOfContractions, otherNumberOfContractions
    integer :: P, Q, R, S
    character(30) :: nameOfSpecie, otherNameOfSpecie
    real(8), allocatable :: auxVector(:)
    integer :: numberOfOptimizationCenters
    integer :: i, j, k
    integer, allocatable :: auxOwnerId(:)
    integer :: pIter, qIter, rIter, sIter, stride, delta
    integer :: u, v, x, y
    integer :: numCartesianP, numCartesianQ, numCartesianR, numCartesianS
    integer :: centerP, centerQ, centerR, centerS
    real(8) :: perm, JKval
    real(8) :: Duv, Dxy
    real(8) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz
    integer, allocatable :: labelsOfContractions(:)
    integer, allocatable :: otherLabelsOfContractions(:)
    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20), otherArguments(20)
    integer :: orderOfMatrix, otherOrderOfMatrix
    type(Matrix) :: densityMatrix, otherDensityMatrix
    real(8) :: charge, otherCharge
    real(8) :: lambda, otherLambda

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if(allocated(output)) deallocate(output)
    allocate(output(numberOfOptimizationCenters,3))

    if(allocated(auxoutput)) deallocate(auxoutput)
    allocate(auxoutput(numberOfOptimizationCenters,3))

    if(allocated(auxOwnerId)) deallocate(auxOwnerId)
    allocate(auxOwnerId(size(ParticleManager_instance)))

    output = 0.0_8
    auxOwnerId = 0

    j = 1
    do i = 1, size(ParticleManager_instance)
       if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
          auxOwnerId(i) = j
          j = j + 1
       else
          auxOwnerId(i) = 0
       end if
    end do


    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    do specieIterator = 1 , MolecularSystem_instance%numberOfQuantumSpecies

       call MolecularSystem_getBasisSet(specieIterator, contractions)
       numberOfContractions = MolecularSystem_getNumberOfContractions(specieIterator)
       nameOfSpecie = trim(MolecularSystem_instance%species(specieIterator)%symbol)

       if(allocated(labelsOfContractions)) deallocate(labelsOfContractions)
       allocate(labelsOfContractions(numberOfContractions))

       j = 1
       do i = 1, numberOfContractions
          labelsOfContractions(i) = j
          j = j + contractions(i)%numCartesianOrbital
       end do

       orderOfMatrix = MolecularSystem_getTotalNumberOfContractions(specieIterator)

       arguments(2) = MolecularSystem_getNameOfSpecie(specieIterator)

       arguments(1) = "DENSITY"
       densityMatrix = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       charge = MolecularSystem_getCharge(specieIterator)
       lambda = MolecularSystem_getLambda(specieIterator)

       do otherSpecieIterator = specieIterator + 1, MolecularSystem_instance%numberOfQuantumSpecies

          auxoutput = 0.0_8
          ! write(*,*) 'Iteradores acoplamiento: ', otherSpecieIterator, specieIterator

          call MolecularSystem_getBasisSet(otherSpecieIterator, otherContractions)
          otherNumberOfContractions = MolecularSystem_getNumberOfContractions(otherSpecieIterator)
          otherNameOfSpecie = trim(MolecularSystem_instance%species(otherSpecieIterator)%symbol)

          if(allocated(otherLabelsOfContractions)) deallocate(otherLabelsOfContractions)
          allocate(otherLabelsOfContractions(otherNumberOfContractions))

          j = 1
          do i = 1, otherNumberOfContractions
             otherLabelsOfContractions(i) = j
             j = j + otherContractions(i)%numCartesianOrbital
          end do

          otherOrderOfMatrix = MolecularSystem_getTotalNumberOfContractions(otherSpecieIterator)

          otherArguments(2) = MolecularSystem_getNameOfSpecie(otherSpecieIterator)

          otherArguments(1) = "DENSITY"
          otherDensityMatrix = &
               Matrix_getFromFile(unit=wfnUnit, rows= int(otherOrderOfMatrix,4), &
               columns= int(otherOrderOfMatrix,4), binary=.true., arguments=otherArguments(1:2))

          otherCharge = MolecularSystem_getCharge(otherSpecieIterator)
          otherLambda = MolecularSystem_getLambda(otherSpecieIterator)

          do P = 1, numberOfContractions
             do Q = P, numberOfContractions
                do R = 1 , otherNumberOfContractions
                   do S = R,  otherNumberOfContractions
                      call DerivativeManager_getElement(&
                           TWOPARTICLE_REPULSION_DERIVATIVES, &
                           auxVector, &
                           i=P, j=Q, k=R, l=S, &
                           nameOfSpecie=nameOfSpecie, &
                           otherNameOfSpecie=otherNameOfSpecie )

                      ! write(*,*) 'shell: ', P,Q,R,S

                      numCartesianP = contractions(P)%numCartesianOrbital
                      numCartesianQ = contractions(Q)%numCartesianOrbital
                      numCartesianR = otherContractions(R)%numCartesianOrbital
                      numCartesianS = otherContractions(S)%numCartesianOrbital

                      centerP = auxOwnerId(contractions(P)%owner)
                      centerQ = auxOwnerId(contractions(Q)%owner)
                      centerR = auxOwnerId(otherContractions(R)%owner)
                      centerS = auxOwnerId(otherContractions(S)%owner)

                      ! write(*,*) 'CENTROS'
                      ! write(*,*) centerP, centerQ, centerR, centerS

                      perm = 1.0
                      if (P /= Q) then
                         perm = perm*2.0
                      end if
                      if (R /= S) then
                         perm = perm*2.0
                      end if
                      ! if (PQ /= RS) then
                      !    perm = perm*2.0
                      ! end if

                      stride = numCartesianP*numCartesianQ*numCartesianR*numCartesianS

                      Ax = 0.0_8
                      Ay = 0.0_8
                      Az = 0.0_8
                      Bx = 0.0_8
                      By = 0.0_8
                      Bz = 0.0_8
                      Cx = 0.0_8
                      Cy = 0.0_8
                      Cz = 0.0_8
                      Dx = 0.0_8
                      Dy = 0.0_8
                      Dz = 0.0_8
                      delta = 0
                      do pIter=0, numCartesianP-1
                         do qIter=0, numCartesianQ-1
                            do rIter=0, numCartesianR-1
                               do sIter=0, numCartesianS-1
                                  u = pIter + labelsOfContractions(P)
                                  v = qIter + labelsOfContractions(Q)
                                  x = rIter + otherLabelsOfContractions(R)
                                  y = sIter + otherLabelsOfContractions(S)
                                  Duv = densityMatrix%values(u,v)
                                  Dxy = otherDensityMatrix%values(x,y)
                                  JKval = perm*Duv*Dxy
                                  Ax = Ax + JKval * auxVector(0 * stride + delta)
                                  Ay = Ay + JKval * auxVector(1 * stride + delta)
                                  Az = Az + JKval * auxVector(2 * stride + delta)
                                  Cx = Cx + JKval * auxVector(3 * stride + delta)
                                  Cy = Cy + JKval * auxVector(4 * stride + delta)
                                  Cz = Cz + JKval * auxVector(5 * stride + delta)
                                  Dx = Dx + JKval * auxVector(6 * stride + delta)
                                  Dy = Dy + JKval * auxVector(7 * stride + delta)
                                  Dz = Dz + JKval * auxVector(8 * stride + delta)
                                  delta = delta + 1
                               end do
                            end do
                         end do
                      end do

                      Bx = -(Ax + Cx + Dx)
                      By = -(Ay + Cy + Dy)
                      Bz = -(Az + Cz + Dz)

                      ! write(*,*) 'x', Ax, Bx, Cx, Dx
                      ! write(*,*) 'y', Ay, By, Cy, Dy
                      ! write(*,*) 'z', Az, Bz, Cz, Dz

                      auxoutput(centerP,1) = auxoutput(centerP,1) + Ax
                      auxoutput(centerP,2) = auxoutput(centerP,2) + Ay
                      auxoutput(centerP,3) = auxoutput(centerP,3) + Az
                      auxoutput(centerQ,1) = auxoutput(centerQ,1) + Bx
                      auxoutput(centerQ,2) = auxoutput(centerQ,2) + By
                      auxoutput(centerQ,3) = auxoutput(centerQ,3) + Bz
                      auxoutput(centerR,1) = auxoutput(centerR,1) + Cx
                      auxoutput(centerR,2) = auxoutput(centerR,2) + Cy
                      auxoutput(centerR,3) = auxoutput(centerR,3) + Cz
                      auxoutput(centerS,1) = auxoutput(centerS,1) + Dx
                      auxoutput(centerS,2) = auxoutput(centerS,2) + Dy
                      auxoutput(centerS,3) = auxoutput(centerS,3) + Dz

                   end do
                end do
             end do
          end do

          auxoutput = auxoutput*charge*otherCharge!*lambda*otherLambda

          ! write(*,"(A)") "----------------------------------------------------------------"
          ! write(*,"(A)") " Gradientes de Acoplamiento"
          ! write(*,*) specieIterator, otherSpecieIterator
          ! write(*,"(A)") "----------------------------------------------------------------"
          ! write(*,"(3(f17.12))") auxoutput(1,1), auxoutput(1,2), auxoutput(1,3)
          ! write(*,"(3(f17.12))") auxoutput(2,1), auxoutput(2,2), auxoutput(2,3)
          ! write(*,"(3(f17.12))") auxoutput(3,1), auxoutput(3,2), auxoutput(3,3)
          ! write(*,"(A)") "----------------------------------------------------------------"


          output = output + auxoutput

       end do
    end do

    k=1
    do i=1, numberOfOptimizationCenters
       do j=1,3
          EnergyGradients_instance%gradients%coupling(k) = output(i,j)
          k = k + 1
       end do
    end do

    ! write(*,"(A)") "----------------------------------------------------------------"
    ! write(*,"(A)") " Gradientes de Acoplamiento"
    ! write(*,"(A)") "----------------------------------------------------------------"
    ! write(*,"(3(f17.12))") output(1,1), output(1,2), output(1,3)
    ! write(*,"(3(f17.12))") output(2,1), output(2,2), output(2,3)
    ! write(*,"(3(f17.12))") output(3,1), output(3,2), output(3,3)
    ! write(*,"(A)") "----------------------------------------------------------------"

    close(wfnUnit)


  end subroutine EnergyGradients_calculateAnalyticCouplingFirstDerivative

  !**
  ! @brief Retorna la componente de la derivada anlitica asociada a las particulas puntuales del sistema
  !
  !**
  subroutine EnergyGradients_calculateFistDerivativeOfPuntualEnergy()
    implicit none
    real(8), allocatable :: output(:,:)
    integer :: i,j,k
    integer, allocatable :: ownerId(:)
    real(8) :: distance, cubeDistance
    real(8) :: deltaOrigin(3)
    real(8) :: Zi, Zj
    integer :: numberOfOptimizationCenters
    integer :: auxIter
    integer :: numberOfPointCharges

    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()
    numberOfPointCharges = size(MolecularSystem_instance%pointCharges)

    if(allocated(ownerId)) deallocate(ownerId)
    allocate(ownerId(size(ParticleManager_instance)))

    if(allocated(output)) deallocate(output)
    allocate(output(numberOfOptimizationCenters,3))

    if(allocated(EnergyGradients_instance%gradients%nuclear)) deallocate(EnergyGradients_instance%gradients%nuclear)
    allocate(EnergyGradients_instance%gradients%nuclear(numberOfOptimizationCenters*3))


    EnergyGradients_instance%gradients%nuclear = 0.0_8

    j = 1

    do i = 1, size(ParticleManager_instance)
       if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
          ownerId(i) = j
          j = j + 1
       else
          ownerId(i) = 0
       end if
    end do

    !print*, size(ParticleManager_instance), numberOfPointCharges

    ! write(*,"(A)") "----------------------------------------------------------------"
    ! write(*,"(A)") " Gradientes Nucleares"
    ! write(*,"(A)") "----------------------------------------------------------------"
    !! Nuclear Gradients
    output = 0.0_8
    do i = 1, size(ParticleManager_instance)
       if(ParticleManager_instance(i)%particlePtr%isQuantum) cycle
       if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
          auxIter = ownerId(i)
          do j = 1, size(ParticleManager_instance)
             if(ParticleManager_instance(j)%particlePtr%isQuantum) cycle
             if(ParticleManager_instance(j)%particlePtr%isCenterOfOptimization) then
                if(i .NE. j) then
                   deltaOrigin = ParticleManager_instance(i)%particlePtr%origin - ParticleManager_instance(j)%particlePtr%origin
                   distance = sqrt( sum( deltaOrigin**2.0_8 ) )
                   cubeDistance = distance*distance*distance
                   Zi = ParticleManager_instance(i)%particlePtr%charge
                   Zj = ParticleManager_instance(j)%particlePtr%charge

                   output(auxIter,1) = output(auxIter,1) - (deltaOrigin(1)*Zi*Zj)/(cubeDistance)
                   output(auxIter,2) = output(auxIter,2) - (deltaOrigin(2)*Zi*Zj)/(cubeDistance)
                   output(auxIter,3) = output(auxIter,3) - (deltaOrigin(3)*Zi*Zj)/(cubeDistance)
                end if
             end if
          end do
       end if
    end do

    k=1
    do i=1, numberOfOptimizationCenters
       do j=1,3
          EnergyGradients_instance%gradients%nuclear(k) = output(i,j)
          k = k + 1
       end do
    end do

    ! write(*,"(3(f17.12))") output(1,1), output(1,2), output(1,3)
    ! write(*,"(3(f17.12))") output(2,1), output(2,2), output(2,3)
    ! write(*,"(3(f17.12))") output(3,1), output(3,2), output(3,3)
    ! write(*,"(A)") "----------------------------------------------------------------"

  end subroutine EnergyGradients_calculateFistDerivativeOfPuntualEnergy

  !>
  !! @brief This routine returns the derivative indices of ket (a,b|
  !! @author Mauricio Rodas 2015
  !<
  subroutine EnergyGradients_getShellPairs(numberOfContractions, output)
    implicit none
    integer, intent(in) :: numberOfContractions
    type(MatrixInteger) :: output
    integer :: P, Q, pairsCount, iter

    call MatrixInteger_constructor(output, 1, 2)

    pairsCount = 0
    do P = 1, numberOfContractions
       do Q = 1, P
          if(pairsCount == 0) then
             output%values(1,1) = P
             output%values(1,2) = Q
          else
             iter = 0
             iter = pairsCount + 1
             call MatrixInteger_addRow(output)
             output%values(iter,1) = P
             output%values(iter,2) = Q
          end if
          pairsCount = pairsCount + 1
       end do
    end do

  end subroutine EnergyGradients_getShellPairs

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine EnergyGradients_exception( typeMessage, description, debugDescription)
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

  end subroutine EnergyGradients_exception

end module EnergyGradients_
