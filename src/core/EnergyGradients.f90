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
  use Vector_
  use MolecularSystem_
  use ParticleManager_
  use ContractedGaussian_
  use MecanicProperties_
  use DerivativeManager_
  use Math_
  use Exception_
  implicit none
  ! use omp_lib

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

  end type EnergyGradients

  type(EnergyGradients), public :: EnergyGradients_instance

  private :: &
       EnergyGradients_calculateAnalyticUncoupledFirstDerivative!, &
  !     EnergyGradients_calculateAnalyticCouplingFirstDerivative, &
  !     EnergyGradients_calculateFistDerivativeOfPuntualEnergy

  public :: &
       EnergyGradients_constructor, &
       EnergyGradients_destructor, &
                                ! EnergyGradients_setName, &
                                ! EnergyGradients_setAnalytic, &
                                ! EnergyGradients_setNumeric, &
                                ! EnergyGradients_setNumericStepSize, &
       EnergyGradients_getDerivative, &
       EnergyGradients_getAnalyticDerivative, &
       EnergyGradients_getNumericGradient!, &
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
    !    write(*,"(A)") "-------------------------"
    !    do j=1, size(evaluationPoint%values)
    !       write(*,"(f12.8)") evaluationPoint%values(j)
    !    end do
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
    real(8), intent(inout) :: gradientVector(:)

    type(Exception) :: ex

    if ( EnergyGradients_instance%order <= 1 ) then

       call EnergyGradients_calculateAnalyticUncoupledFirstDerivative(gradientVector)

       ! if ( ParticleManager_getNumberOfQuantumSpecies() > 1) &
       !      output = output + EnergyGradients_calculateAnalyticCouplingFirstDerivative()

       ! output = output + EnergyGradients_calculateFistDerivativeOfPuntualEnergy()

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
  subroutine EnergyGradients_calculateAnalyticUncoupledFirstDerivative(gradientVector)
    implicit none
    real(8), intent(inout) :: gradientVector(:)
    character(30) :: nameOfSpecie
    type(ContractedGaussian), allocatable :: contractions(:)
    integer :: numberOfContractions
    real(8), allocatable :: auxVector(:), auxVector2(:), auxVector3(:)
    integer :: specieIterator
    integer :: ocupationNumber
    integer :: orderOfMatrix
    type(Matrix) :: matrixOfEigenvectors
    type(Matrix) :: densityMatrix
    type(Vector) :: vectorOfEigenvalues
    type(Matrix) :: auxWeightDensity
    type(Matrix) :: weightDensityMatrix
    integer(8) :: numberOfIntegrals
    integer :: i
    integer :: j
    integer :: k
    integer :: u
    integer :: v
    integer :: r
    integer :: s
    integer :: P, Q, pIter, qIter, auxIter, A
    integer :: numCartesianP, numCartesianQ
    integer :: centerA, centerB
    integer :: numberOfOptimizationCenters
    real(8) :: perm, distance, cubeDistance, Zi, Zj, Vval
    real(8) :: deltaOrigin(3)
    integer, allocatable ::  labelsOfContractions(:)
    real(8), allocatable :: auxNuclear(:,:), auxKinetic(:,:), auxPotential(:,:), auxOverlap(:,:), auxCoulomb(:,:), auxExchange(:,:)
    integer, allocatable :: auxOwnerId(:)
    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20)

    wfnFile = "lowdin.wfn"
    wfnUnit = 20


    numberOfOptimizationCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if(allocated(auxNuclear)) deallocate(auxNuclear)
    allocate(auxNuclear(numberOfOptimizationCenters,3))

    if(allocated(auxKinetic)) deallocate(auxKinetic)
    allocate(auxKinetic(numberOfOptimizationCenters,3))

    if(allocated(auxPotential)) deallocate(auxPotential)
    allocate(auxPotential(numberOfOptimizationCenters,3))

    if(allocated(auxOverlap)) deallocate(auxOverlap)
    allocate(auxOverlap(numberOfOptimizationCenters,3))

    if(allocated(auxCoulomb)) deallocate(auxCoulomb)
    allocate(auxCoulomb(numberOfOptimizationCenters,3))

    if(allocated(auxExchange)) deallocate(auxExchange)
    allocate(auxExchange(numberOfOptimizationCenters,3))

    if(allocated(auxOwnerId)) deallocate(auxOwnerId)
    allocate(auxOwnerId(size(ParticleManager_instance)))

    auxNuclear = 0.0_8
    auxKinetic = 0.0_8
    auxPotential = 0.0_8
    auxOverlap = 0.0_8
    auxCoulomb = 0.0_8
    auxExchange = 0.0_8
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

    do specieIterator=1, MolecularSystem_instance%numberOfQuantumSpecies

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

       ocupationNumber = MolecularSystem_getOcupationNumber( specieIterator )

       arguments(2) = MolecularSystem_getNameOfSpecie(specieIterator)

       arguments(1) = "DENSITY"
       densityMatrix = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       arguments(1) = "COEFFICIENTS"
       matrixOfEigenvectors = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = int(orderOfMatrix,4), &
            unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
            output = vectorOfEigenvalues )

       call Matrix_constructor(auxWeightDensity, int(orderOfMatrix,8), int(ocupationNumber,8))

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

       gradientVector = 0.0_8

       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes Nucleares"
       write(*,"(A)") "----------------------------------------------------------------"
       !! Nuclear Gradients
       auxNuclear = 0.0_8
       do i = 1, size(ParticleManager_instance)
          if(ParticleManager_instance(i)%particlePtr%isCenterOfOptimization) then
             auxIter = auxOwnerId(i)

             do j = 1, size(ParticleManager_instance)
                if(ParticleManager_instance(j)%particlePtr%isCenterOfOptimization) then
                   if(i .NE. j) then
                      deltaOrigin = ParticleManager_instance(i)%particlePtr%origin - ParticleManager_instance(j)%particlePtr%origin
                      distance = sqrt( sum( deltaOrigin**2.0_8 ) )
                      cubeDistance = distance*distance*distance
                      Zi = ParticleManager_instance(i)%particlePtr%charge
                      Zj = ParticleManager_instance(j)%particlePtr%charge

                      auxNuclear(auxIter,1) = auxNuclear(auxIter,1) - (deltaOrigin(1)*Zi*Zj)/(cubeDistance)
                      auxNuclear(auxIter,2) = auxNuclear(auxIter,2) - (deltaOrigin(2)*Zi*Zj)/(cubeDistance)
                      auxNuclear(auxIter,3) = auxNuclear(auxIter,3) - (deltaOrigin(3)*Zi*Zj)/(cubeDistance)
                   end if
                end if
             end do
          end if
       end do

       write(*,"(3(f17.12))") auxNuclear(1,1), auxNuclear(1,2), auxNuclear(1,3)
       write(*,"(3(f17.12))") auxNuclear(2,1), auxNuclear(2,2), auxNuclear(2,3)
       write(*,"(3(f17.12))") auxNuclear(3,1), auxNuclear(3,2), auxNuclear(3,3)
       write(*,"(A)") "----------------------------------------------------------------"

       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes Cineticos"
       write(*,"(A)") "----------------------------------------------------------------"
       !! Kinetic Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             call DerivativeManager_getElement( KINETIC_DERIVATIVES, auxVector, i=P, j=Q, nameOfSpecie=nameOfSpecie )

             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerA = auxOwnerId(contractions(P)%owner) !! aP
             centerB = auxOwnerId(contractions(Q)%owner) !! aQ

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
                   auxKinetic(centerA,1) = auxKinetic(centerA,1) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerA,2) = auxKinetic(centerA,2) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerA,3) = auxKinetic(centerA,3) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerB,1) = auxKinetic(centerB,1) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerB,2) = auxKinetic(centerB,2) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxKinetic(centerB,3) = auxKinetic(centerB,3) + perm*densityMatrix%values(u,v)*auxVector(i)
                   i = i + 1
                end do
             end do

          end do
       end do

       write(*,"(3(f17.12))") auxKinetic(1,1), auxKinetic(1,2), auxKinetic(1,3)
       write(*,"(3(f17.12))") auxKinetic(2,1), auxKinetic(2,2), auxKinetic(2,3)
       write(*,"(3(f17.12))") auxKinetic(3,1), auxKinetic(3,2), auxKinetic(3,3)
       write(*,"(A)") "----------------------------------------------------------------"


       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes de Potencial"
       write(*,"(A)") "----------------------------------------------------------------"
       !! Potential Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerA = auxOwnerId(contractions(P)%owner) !! aP
             centerB = auxOwnerId(contractions(Q)%owner) !! aQ

             if (P == Q) then
                perm = 1.0
             else
                perm = 2.0
             end if

             call DerivativeManager_getElement( ATTRACTION_DERIVATIVES, &
                  auxVector2, i=P, j=Q, nameOfSpecie=nameOfSpecie, A=centerA, B=centerB )

             i = 0
             j = 0
             k = 0
             do A=1, numberOfOptimizationCenters
                i = 3*(A-1)*numCartesianP*numCartesianQ + 0*numCartesianP*numCartesianQ
                j = 3*(A-1)*numCartesianP*numCartesianQ + 1*numCartesianP*numCartesianQ
                k = 3*(A-1)*numCartesianP*numCartesianQ + 2*numCartesianP*numCartesianQ
                do pIter = 0, numCartesianP - 1
                   do qIter = 0, numCartesianQ - 1
                      u = pIter + labelsOfContractions(P)
                      v = qIter + labelsOfContractions(Q)
                      Vval = perm*densityMatrix%values(u,v)

                      auxPotential(A,1) = auxPotential(A,1) + Vval*auxVector2(i)

                      i = i + 1
                      auxPotential(A,2) = auxPotential(A,2) + Vval*auxVector2(j)

                      j = j + 1
                      auxPotential(A,3) = auxPotential(A,3) + Vval*auxVector2(k)

                      k = k + 1

                   end do
                end do
             end do
          end do
       end do

       write(*,"(3(f17.12))") auxPotential(1,1), auxPotential(1,2), auxPotential(1,3)
       write(*,"(3(f17.12))") auxPotential(2,1), auxPotential(2,2), auxPotential(2,3)
       write(*,"(3(f17.12))") auxPotential(3,1), auxPotential(3,2), auxPotential(3,3)
       write(*,"(A)") "----------------------------------------------------------------"

       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes de Overlap"
       write(*,"(A)") "----------------------------------------------------------------"
       !! Overlap Gradients
       do P = 1, numberOfContractions
          do Q = 1, P
             call DerivativeManager_getElement( OVERLAP_DERIVATIVES, auxVector3, i=P, j=Q, nameOfSpecie=nameOfSpecie )

             numCartesianP = contractions(P)%numCartesianOrbital  !! nP
             numCartesianQ = contractions(Q)%numCartesianOrbital  !! nQ
             centerA = auxOwnerId(contractions(P)%owner) !! aP
             centerB = auxOwnerId(contractions(Q)%owner) !! aQ

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
                   auxOverlap(centerA,1) = auxOverlap(centerA,1) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerA,2) = auxOverlap(centerA,2) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerA,3) = auxOverlap(centerA,3) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerB,1) = auxOverlap(centerB,1) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerB,2) = auxOverlap(centerB,2) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

             do pIter = 0, numCartesianP - 1
                do qIter = 0, numCartesianQ - 1

                   u = pIter + labelsOfContractions(P)
                   v = qIter + labelsOfContractions(Q)
                   auxOverlap(centerB,3) = auxOverlap(centerB,3) - perm*weightDensityMatrix%values(u,v)*auxVector3(i)
                   i = i + 1
                end do
             end do

          end do
       end do

       write(*,"(3(f17.12))") auxOverlap(1,1), auxOverlap(1,2), auxOverlap(1,3)
       write(*,"(3(f17.12))") auxOverlap(2,1), auxOverlap(2,2), auxOverlap(2,3)
       write(*,"(3(f17.12))") auxOverlap(3,1), auxOverlap(3,2), auxOverlap(3,3)
       write(*,"(A)") "----------------------------------------------------------------"

       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes de Coulomb"
       write(*,"(A)") "----------------------------------------------------------------"


       write(*,"(3(f17.12))") auxCoulomb(1,1), auxCoulomb(1,2), auxCoulomb(1,3)
       write(*,"(3(f17.12))") auxCoulomb(2,1), auxCoulomb(2,2), auxCoulomb(2,3)
       write(*,"(3(f17.12))") auxCoulomb(3,1), auxCoulomb(3,2), auxCoulomb(3,3)
       write(*,"(A)") "----------------------------------------------------------------"

       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(A)") " Gradientes de Intercambio"
       write(*,"(A)") "----------------------------------------------------------------"
       write(*,"(3(f17.12))") auxExchange(1,1), auxExchange(1,2), auxExchange(1,3)
       write(*,"(3(f17.12))") auxExchange(2,1), auxExchange(2,2), auxExchange(2,3)
       write(*,"(3(f17.12))") auxExchange(3,1), auxExchange(3,2), auxExchange(3,3)
       write(*,"(A)") "----------------------------------------------------------------"


    end do

    close(wfnUnit)

  end subroutine EnergyGradients_calculateAnalyticUncoupledFirstDerivative

  !         !**
  !         ! @brief  Calcula la primera derivada analitica para la matriz de Fock sin incluir el termino de
  !         !               acoplamiento interspecie
  !         !**
  !         function EnergyGradients_calculateAnalyticCouplingFirstDerivative() result(output)
  !             implicit none
  !             real(8) :: output

  !             character(30) :: nameOfSpecie
  !             character(30) :: nameOfOtherSpecie
  !             integer :: nucleiAndComponent(2)
  !             integer :: specieIterator
  !             integer :: otherSpecieIterator
  !             integer :: ocupationNumber
  !             integer :: otherOcupationNumber
  !             integer :: orderOfMatrix
  !             integer :: orderOfOtherMatrix
  !             real(8) :: couplingEnergyDerivative
  !             real(8) :: auxCouplingEnergyDerivative
  !             type(Matrix), allocatable :: derivativeOfRepulsionVector(:)
  !             real(8) :: lambda, otherLambda
  !             type(Matrix) :: matrixOfEigenvectors
  !             type(Matrix) :: otherMatrixOfEigenvectors
  !             integer(8) :: numberOfIntegrals
  !             integer :: iteratorOfMatrix
  !             integer :: a
  !             integer :: b
  !             integer :: u
  !             integer :: v
  !             integer :: r
  !             integer :: s

  !             nucleiAndComponent = ParticleManager_getCenterOfOptimization( EnergyGradients_instance%components(1) )

  !             if ( ParticleManager_isComponentFixed( nucleiAndComponent(1),nucleiAndComponent(2) ) ) then
  !                 output=0.0_8
  !                 return
  !             end if

  !             iteratorOfMatrix =  ParticleManager_getNumberOfQuantumSpecies() &
  !                                 * ( ParticleManager_getNumberOfQuantumSpecies() - 1 ) / 2

  !             allocate( derivativeOfRepulsionVector(iteratorOfMatrix) )

  !             couplingEnergyDerivative = 0.0_8
  !             iteratorOfMatrix = 0

  !             do specieIterator = 1 , ParticleManager_getNumberOfQuantumSpecies()

  !                 !! Obtiene vectores propios para una de las especies
  !                 nameOfSpecie =  ParticleManager_getNameOfSpecie( specieIterator )

  !                 matrixOfEigenvectors = MolecularSystem_instance%getWaveFunctionCoefficients( trim(nameOfSpecie) )

  !                 orderOfMatrix = ParticleManager_getNumberOfContractions( specieIterator )
  !                 ocupationNumber = ParticleManager_getOcupationNumber( specieIterator )
  !                 lambda =ParticleManager_getLambda( specieIterator )


  !                 do otherSpecieIterator = specieIterator + 1, ParticleManager_getNumberOfQuantumSpecies()

  !                     iteratorOfMatrix = iteratorOfMatrix + 1

  !                     !! Obtiene vectores propios para otra especie
  !                     nameOfOtherSpecie =  ParticleManager_getNameOfSpecie( otherSpecieIterator )

  !                     otherMatrixOfEigenvectors = MolecularSystem_instance%getWaveFunctionCoefficients( trim(nameOfOtherSpecie) )

  !                     orderOfOtherMatrix = ParticleManager_getNumberOfContractions( otherSpecieIterator )
  !                     otherOcupationNumber = ParticleManager_getOcupationNumber( otherSpecieIterator )
  !                     otherLambda =ParticleManager_getLambda( otherSpecieIterator )
  !                     auxCouplingEnergyDerivative=0.0_8

  !                     do a=1, ocupationNumber
  !                         do b=1,otherOcupationNumber

  !                 !!CONVERTIR EN UN UNICO CICLO
  !                         !!**************************************************************
  !                         !! Suma derivadas de  la energia iterativamente sobre pares especies.
  !                         !!
  !                         do u = 1 , orderOfMatrix
  !                             do v = 1 , orderOfMatrix
  !                                 do r = 1 , orderOfOtherMatrix
  !                                     do s = 1 , orderOfOtherMatrix

  !                                         auxCouplingEnergyDerivative = auxCouplingEnergyDerivative &
  !                                             + ( matrixOfEigenvectors%values(u,a) * matrixOfEigenvectors%values(v,a) &
  !                                             * otherMatrixOfEigenvectors%values(r,b) * otherMatrixOfEigenvectors%values(s,b) &
  !                                             * dRepulsion( u,v,r,s )  )

  !                                     end do
  !                                 end do
  !                             end do
  !                         end do
  !                         !!**************************************************************

  !                         end do
  !                     end do

  !                     couplingEnergyDerivative = couplingEnergyDerivative + auxCouplingEnergyDerivative &
  !                         * ParticleManager_getCharge( specieID=specieIterator ) &
  !                         * ParticleManager_getCharge( specieID=otherSpecieIterator )  * lambda * otherLambda

  !                 end do

  !             end do
  !             !!***************************************************************************

  !             output = couplingEnergyDerivative

  !             do a=1, size(derivativeOfRepulsionVector)
  !                 call Matrix_destructor( derivativeOfRepulsionVector(a) )
  !             end do
  !             deallocate( derivativeOfRepulsionVector )

  !             contains

  !                 function dRepulsion(index_u, index_v, index_r, index_s) result( output )
  !                     implicit none
  !                     integer :: index_u
  !                     integer :: index_v
  !                     integer :: index_r
  !                     integer :: index_s
  !                     real(8) :: output
  !                     integer(8) :: numberOfIntegrals

  !                     integer(8) :: indexOfIntegral

  !                     if ( .not.allocated( derivativeOfRepulsionVector(iteratorOfMatrix)%values ) ) then

  !                         numberOfIntegrals =     ( orderOfMatrix    *   ( orderOfMatrix + 1.0_8) / 2.0_8 ) * &
  !                                             ( orderOfOtherMatrix * ( orderOfOtherMatrix + 1.0_8 ) / 2.0_8 )
  !                         call Matrix_constructor( derivativeOfRepulsionVector(iteratorOfMatrix), numberOfIntegrals, 1_8, Math_NaN )

  !                     end if

  !                     indexOfIntegral =   IndexMap_tensorR4ToVector(index_u,index_v,index_r,index_s, orderOfMatrix, orderOfOtherMatrix )

  !                     if ( isNaN( derivativeOfRepulsionVector( iteratorOfMatrix )%values(indexOfIntegral, 1 ) ) ) then

  !                         output = ContractedGaussian_repulsionDerivative( &
  ! ParticleManager_instance%particlesPtr(ParticleManager_instance%idsOfContractionsForSpecie(specieIterator)%contractionID(index_u)%particleID)%basis%contractions(ParticleManager_instance%idsOfContractionsForSpecie(specieIterator)%contractionID(index_u)%contractionIDInParticle), &
  !                 ParticleManager_instance%particlesPtr(ParticleManager_instance%idsOfContractionsForSpecie(specieIterator)%contractionID(index_v)%particleID)%basis%contractions(ParticleManager_instance%idsOfContractionsForSpecie(specieIterator)%contractionID(index_v)%contractionIDInParticle), &
  !                 ParticleManager_instance%particlesPtr(ParticleManager_instance%idsOfContractionsForSpecie(otherSpecieIterator)%contractionID(index_r)%particleID)%basis%contractions(ParticleManager_instance%idsOfContractionsForSpecie(otherSpecieIterator)%contractionID(index_r)%contractionIDInParticle), &
  !                 ParticleManager_instance%particlesPtr(ParticleManager_instance%idsOfContractionsForSpecie(otherSpecieIterator)%contractionID(index_s)%particleID)%basis%contractions(ParticleManager_instance%idsOfContractionsForSpecie(otherSpecieIterator)%contractionID(index_s)%contractionIDInParticle), &
  !                             nucleiAndComponent(1), nucleiAndComponent(2) )

  !                         derivativeOfRepulsionVector( iteratorOfMatrix )%values(indexOfIntegral, 1 ) = output

  !                     else

  !                         output = derivativeOfRepulsionVector( iteratorOfMatrix )%values(indexOfIntegral, 1 )

  !                     end if

  !                 end function dRepulsion

  !         end function EnergyGradients_calculateAnalyticCouplingFirstDerivative


  !         !**
  !         ! @brief Retorna la componente de la derivada anlitica asociada a las particulas puntuales del sistema
  !         !
  !         !**
  !         function EnergyGradients_calculateFistDerivativeOfPuntualEnergy() result( output )
  !             implicit none
  !             real(8) :: output

  !             integer :: i
  !             integer :: j
  !             integer :: center
  !             integer :: otherCenter
  !             integer :: auxKronecker
  !             integer :: nucleiAndComponent(2)
  !             real(8) :: originComponent(3)


  !             nucleiAndComponent= ParticleManager_getCenterOfOptimization(EnergyGradients_instance%components(1) )

  !             if ( ParticleManager_isComponentFixed( nucleiAndComponent(1),nucleiAndComponent(2) ) ) then
  !                 output=0.0_8
  !                 return
  !             end if

  !             output = 0.0_8

  !             do i = 1 , ParticleManager_getNumberOfPuntualParticles()
  !                 do j = i + 1 , ParticleManager_getNumberOfPuntualParticles()

  !                     center= ParticleManager_getOwnerOfPuntualParticle( i )
  !                     otherCenter= ParticleManager_getOwnerOfPuntualParticle( j )

  !                     auxKronecker =  Math_kroneckerDelta( nucleiAndComponent(1), center ) &
  !                                     - Math_kroneckerDelta( nucleiAndComponent(1), otherCenter )

  !                     if (abs( auxKronecker ) > 0 ) then

  !                         originComponent=    ParticleManager_getOriginOfPuntualParticle(i) &
  !                                         -ParticleManager_getOriginOfPuntualParticle( j)

  !                         output = output &
  !                             +  ( (-1.0_8 * ParticleManager_getChargeOfPuntualParticle( i )  &
  !                             *  ParticleManager_getChargeOfPuntualParticle( j )&
  !                             / ( sum( ( ParticleManager_getOriginOfPuntualParticle(i) &
  !                             - ParticleManager_getOriginOfPuntualParticle(j) )**2.0_8 )  )**(1.5_8)  ) &
  !                             * originComponent( nucleiAndComponent(2) ) * auxKronecker )
  !                     end if

  !                 end do
  !             end do

  !         end function EnergyGradients_calculateFistDerivativeOfPuntualEnergy

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
