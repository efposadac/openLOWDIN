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
!! @brief Geometry Optimizer Module.
!!        This module contains all basic functions of the geometry optimizer
!! @author  J.M. Rodas
!! @author  S.A. Gonzalez
!!
!! <b> Creation date : </b> 2008-10-10
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-10-10 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Basics functions and functions has been created using trust region method
!!   - <tt> 2015-02-23 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0 and prepare the module for new optimizers
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module GeometryOptimizer_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use EnergyGradients_
  use Hessians_
  use InternalCoordinates_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  ! use Input_Parameters_
  ! use Units_
  ! use CalculateProperties_
  ! use InternalCoordinates_
  use MecanicProperties_
  use TrustRegionOptimizer_
  ! use TrustRegionIterative_
  ! use RHF_
  ! use MollerPlesset_
  ! use Input_Parsing_
  ! use ExternalSoftware_
  use Solver_
  use Exception_
  implicit none

  type, public :: GeometryOptimizer
     
     character(30) :: name
     integer :: numberOfIterations
     integer :: numberOfIndependVariables
     type(Vector) :: valuesOfIndependentVariables
     type(Vector) :: gradient
     real(8) :: functionValue
     real(8) :: stepSize
     ! type(InternalCoordinates) :: primitivesCoordinates
     type(TrustRegionOptimizer) :: trustRegion
     ! type(TrustRegionIterative) :: trustRegionIterative
     logical :: isInitialExecution
     logical :: isInstanced
     ! type(Input_Parameters) :: LOWDIN_bkp
     ! type(Solver), pointer :: solverPtr

  end type GeometryOptimizer

  type(GeometryOptimizer), pointer , private :: this_pointer
  type(GeometryOptimizer), target :: GeometryOptimizer_instance

  public :: &
       GeometryOptimizer_constructor, &
       GeometryOptimizer_destructor, &
       GeometryOptimizer_run, &
  ! GeometryOptimizer_setName, &
  ! GeometryOptimizer_getNumberOfIterations, &
       GeometryOptimizer_getFunctionValue, &
       GeometryOptimizer_getGradient

contains
  !**
  ! @brief Define el constructor
  !**
  subroutine GeometryOptimizer_constructor( this )
    implicit none
    type(GeometryOptimizer) :: this
    type(Matrix) :: initialHessian
    type(Vector) :: initialGeometry
    type(Hessians) :: hessiansInfo
    integer :: i
    integer :: numberOfCenterofOptimization

    this%isInstanced =.true.
    ! CONTROL_instance%LAST_STEP=.false.

    numberOfCenterofOptimization = ParticleManager_getNumberOfCentersOfOptimization()

    this%name = ""
    this%isInitialExecution=.true.
    initialGeometry = ParticleManager_getPositionOfCenterOfOptimizacion()

    this%numberOfIndependVariables = size( initialGeometry%values )

    ! call Input_Parameters_copyConstructor(this%LOWDIN_bkp, Parameters)
    call Hessians_constructor( hessiansInfo )
    initialHessian = Hessians_getEmpirical( hessiansInfo , MolecularSystem_instance )

    !! Debug Mauricio Rodas
    ! do i=1, size(initialHessian%values, dim=1)
    !    write(*,"(A,2x,f12.8)") "Diagonal hesiana: ", initialHessian%values(i,i)
    ! end do

    ! this%solverPtr => ssolver

    select case ( trim(CONTROL_instance%MINIMIZATION_METHOD))

        case("TR")
            !!by default
            call TrustRegionOptimizer_constructor( this%trustRegion, initialGeometry%values )
            call TrustRegionOptimizer_setInitialHessian( this%trustRegion, initialHessian%values )

    !     case("TRI")

    !         call TrustRegionIterative_constructor( this%trustRegionIterative, initialGeometry%values )
    !         call TrustRegionIterative_setInitialHessian( this%trustRegionIterative, initialHessian%values )

        case default

            call TrustRegionOptimizer_constructor( this%trustRegion, initialGeometry%values )
            call TrustRegionOptimizer_setInitialHessian( this%trustRegion, initialHessian%values )

    end select

    call EnergyGradients_constructor()
    call Matrix_destructor(initialHessian)
    call Vector_destructor(initialGeometry)
    call Hessians_destructor( hessiansInfo)

  end subroutine GeometryOptimizer_constructor

        !**
        ! @brief Define el destructor
        !**
  subroutine GeometryOptimizer_destructor(this)
    implicit none

    type(GeometryOptimizer) :: this

    ! this%solverPtr => null()

    select case ( trim(CONTROL_instance%MINIMIZATION_METHOD))
    case("TR")
       call TrustRegionOptimizer_destructor( this%trustRegion )
    ! case("TRI")
    !    call TrustRegionIterative_destructor( this%trustRegionIterative )
    end select

    call EnergyGradients_destructor()

  end subroutine GeometryOptimizer_destructor

  !>
  !! @brief Muestra en la la  unidad de salida el estimado del minimo
  !!
  !!  @warning Anque no se define como privada de la clase por propositos de conveniencia
  !!          no debe ser empleada fuera de la clase
  !<
  subroutine GeometryOptimizer_showMinimum( evaluationPoint, gradient, functionValue )
    implicit none
    real(8) :: evaluationPoint(:)
    real(8) :: gradient(:)
    real(8) :: functionValue

    integer :: i
    integer :: k
    real(8) :: origin(3)
    integer :: totalNumberOfParticles

    totalNumberOfParticles = size( ParticleManager_instance )

    write (6,*) ""
    write (6,"(T20,A25)") "COORDINATES: "//trim(CONTROL_instance%UNITS)
    write (6,"(A10,A16,A20,A20)") "Particle","<x>","<y>","<z>"
    do i = 1, totalNumberOfParticles!ParticleManager_getTotalNumberOfParticles()

       if ( ParticleManager_isCenterOfOptimization( i ) ) then
          origin = ParticleManager_getOrigin( iterator = i )
          if ( CONTROL_instance%UNITS=="ANGSTROMS") then
             origin = origin * AMSTRONG
          end if
#ifdef intel
          write (6,"(A10,<3>F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), origin(1), origin(2), origin(3)
#else
          write (6,"(A10,3F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), origin(1), origin(2), origin(3)
#endif
       end if

    end do


    ! print *,""
    write (6,"(A)") ""
    write (6,"(T20,A25)") "GRADIENT: HARTREE/BOHR"
    write (6,"(A10,A16,A20,A20)") "Particle","dE/dx","dE/dy","dE/dz"
    k=1
    do i = 1, totalNumberOfParticles!ParticleManager_getTotalNumberOfParticles()

       if ( ParticleManager_isCenterOfOptimization( i ) .and. k < size(gradient)) then

          write (6,"(A10,3F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), &
               gradient(k), gradient(k+1),gradient(k+2)
          k=k+3

       end if

    end do
    write (6,*) ""
    write (6,"(T10,A24,F20.10)") "TOTAL ENERGY (Hartree) =", functionValue

  end subroutine GeometryOptimizer_showMinimum

  !>
  !! @brief Ejecuta el minimizador hasta encontrar un minimo de la estructura
  !<
!!!!!!!16/07/2012
  subroutine GeometryOptimizer_run(this )
    implicit none
    type(GeometryOptimizer), target :: this

    real(8) :: auxValue
    logical :: isMinumum
    type(Matrix) :: initialHessian
    type(Hessians) :: hessiansInfo
    logical :: lastStep

    print   *,""
    print *, " BEGIN GEOMETRY OPTIMIZATION: "
    print *,"============================"
    print   *,""
    isMinumum = .false.
    lastStep = .false.

    open(unit=40, file="lowdin.dat", status="replace", form="formatted")
    
    !!save all options
    call CONTROL_save(40, lastStep)

    close(40)

    CONTROL_instance%SCF_CONVERGENCE_CRITERIUM="energy"
    if( CONTROL_instance%MINIMIZATION_WITH_SINGLE_POINT ) CONTROL_instance%SCF_CONVERGENCE_CRITERIUM = "energy"

    this_pointer => this

    !!************************************************************
    !! Inicia proceso de minimizacion
    !!***

    select case ( trim(CONTROL_instance%MINIMIZATION_METHOD))

    case("TR")

       do while ( .not.isMinumum )

          write (6,*) ""

          call TrustRegionOptimizer_iterate( this%trustRegion, &
               GeometryOptimizer_getFunctionValue, &
               GeometryOptimizer_getGradient, &
               GeometryOptimizer_projectGradient, &
               GeometryOptimizer_projectHessiane, &
               GeometryOptimizer_showMinimum )

          isMinumum = TrustRegionOptimizer_isMinimum( this%trustRegion )

          ! if  ( ( TrustRegionOptimizer_getNumberOfIterations(this%trustRegion) > &
          !      this%numberOfIndependVariables ) .and. (CONTROL_instance%RESTART_OPTIMIZATION) ) then

          !    call TrustRegionOptimizer_restart(this%trustRegion)
          !    call Hessians_constructor(hessiansInfo)
          !    initialHessian =Hessians_getCartesianForceConstants( hessiansInfo , MolecularSystem_instance )
          !    call TrustRegionOptimizer_setInitialHessian( this%trustRegion, initialHessian%values )
          !    call Matrix_destructor(initialHessian)
          !    call Hessians_destructor(hessiansInfo)

          ! end if

       end do

       write(6,*) ""
       write (6,"(T20,A30)") " FINAL GEOMETRY: AMSTRONG"
       write (6,"(T18,A35)") "------------------------------------------"
       write(6,*) ""
       write (6,"(T20,A25,ES20.8)") "FINAL GRADIENT (RMS):",TrustRegionOptimizer_getGradient(this%trustRegion)

    ! case("TRI")

    !    do while ( .not.isMinumum )

    !       write (6,*) ""

    !       call TrustRegionIterative_iterate( this%trustRegionIterative, &
    !            GeometryOptimizer_getFunctionValue, &
    !            GeometryOptimizer_getGradient, &
    !            GeometryOptimizer_projectGradient, &
    !            GeometryOptimizer_projectHessiane, &
    !            GeometryOptimizer_showMinimum )

    !       isMinumum = TrustRegionIterative_isMinimum( this%trustRegionIterative )

    !       if  ( ( TrustRegionIterative_getNumberOfIterations(this%trustRegionIterative) > &
    !            this%numberOfIndependVariables ) .and. (CONTROL_instance%RESTART_OPTIMIZATION) ) then

    !          call TrustRegionIterative_restart(this%trustRegionIterative)
    !          call Hessians_constructor(hessiansInfo)
    !          initialHessian =Hessians_getCartesianForceConstants( hessiansInfo, &
    !               MolecularSystem_instance )
    !          call TrustRegionIterative_setInitialHessian( this%trustRegionIterative, initialHessian%values )
    !          call Matrix_destructor(initialHessian)
    !          call Hessians_destructor(hessiansInfo)

    !       end if

    !    end do

    !    write(6,*) ""
    !    write (6,"(T20,A30)") " FINAL GEOMETRY: AMSTRONG"
    !    write (6,"(T18,A35)") "------------------------------------------"
    !    write(6,*) ""
    !    write (6,"(T20,A25,ES20.8)") "FINAL GRADIENT (RMS):",TrustRegionIterative_getGradient(this%trustRegionIterative)

    end select

    write(6,*) ""
    call MolecularSystem_showCartesianMatrix()
    call MolecularSystem_showDistanceMatrix()
    call MolecularSystem_showZMatrix( MolecularSystem_instance )

    ! !! Realiza un calculo final en la geometria de equilibrio
    CONTROL_instance%OPTIMIZE=.false.
    CONTROL_instance%SCF_CONVERGENCE_CRITERIUM="density"
    lastStep = .true.

    open(unit=40, file="lowdin.dat", status="replace", form="formatted")
    
    !!save all options
    call CONTROL_save(40, lastStep)

    close(40)

    ! CONTROL_instance%LAST_STEP = .true.
    ! call SCF_Global_setStopingThreshold( CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE, &
    !      CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE )
    ! call Solver_run()

    print *,""
    print *,"END GEOMETRY OPTIMIZATION "
    print *,""

    ! this%solverPtr%withProperties=.true.
    call Solver_run()

  end subroutine GeometryOptimizer_run


!         !>
!         !! @brief Ajusta el nombre del minimizador empleado
!         !!
!         !! @todo Falta implementacion completa
!         !<
!         subroutine GeometryOptimizer_setName(this)
!             implicit none
!             type(GeometryOptimizer) :: this

!         end subroutine GeometryOptimizer_setName


!         !>
!         !! @brief retorna el numero de iteraciones actual del minimizador
!         !!
!         !! @todo Falta por implementar
!         !<
!         function GeometryOptimizer_getNumberOfIterations(this) result(output)
!             implicit none
!             type(GeometryOptimizer) :: this
!             integer :: output

!             output = 0

!         end function GeometryOptimizer_getNumberOfIterations

  !>
  !! @brief retorna el valor la funcion a minimiza en el punto especificado
  !!      esta funcion debe ser llamada privia asigacion de apuntador this_pointer
  !!
  !! @warning Anque no se define como privada de la clase por propositos de conveniencia
  !!          no debe ser empleada fuera de la clase
  !<
  function GeometryOptimizer_getFunctionValue( evaluationPoint ) result(output)
    implicit none
    real(8):: evaluationPoint(:)
    real(8) :: output
    real(8) :: totalEnergy
    character(50) :: wfnFile
    integer :: wfnUnit

    type(Vector) :: valuesOfIndependentVariables

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    call Vector_constructor( valuesOfIndependentVariables, size(evaluationPoint) )
    valuesOfIndependentVariables%values = evaluationPoint

    !! Ajusta el origen de las particulas presentes en el sistema
    call ParticleManager_setParticlesPositions( valuesOfIndependentVariables)
    call Vector_destructor( valuesOfIndependentVariables )

    ! select case ( trim(CONTROL_instance%ENERGY_CALCULATOR) )

    ! case("internal")

       ! if ( this_pointer%isInitialExecution ) then
       !    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION >= 2) call EnergyGradients_setNumeric()
       ! else
       !    call Solver_reset( this_pointer%solverPtr )
       ! end if

       !! Restituye las variables originales para el control de programa
       ! call Input_Parameters_copyConstructor(Parameters, this_pointer%LOWDIN_bkp)

       ! this_pointer%solverPtr%withProperties=.false.
    call Solver_run()
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])
    close(wfnUnit)
   
    output = totalEnergy
    ! output = this_pointer%solverPtr%energy

    this_pointer%isInitialExecution=.false.
    ! print "(A1$)","."

    ! case("external")

    !    call EnergyGradients_setNumeric()
    !    if ( .not.external_instance%isInstanced ) &
    !         call ExternalSoftware_constructor( external_instance )

    !    output = ExternalSoftware_getEnergy( external_instance, &
    !         withCounterPoiseCorrection = CONTROL_instance%OPTIMIZE_WITH_CP_CORRECTION )
    ! end select

  end function GeometryOptimizer_getFunctionValue

  ! >
  ! @brief Retorna el gradiente obtenido en la ultima iteracion
  !      Esta funcion solo debe emplease previa asignacion del apuntador this_pointer
  
  ! @warning Anque no se define como privada de la clase por propositos de conveniencia
  !          no debe ser empleada fuera de la clase
  ! <
  subroutine GeometryOptimizer_getGradient( evaluationPoint, gradientVector )
    implicit none
    real(8) :: evaluationPoint(:)
    real(8) :: gradientVector(:)
    type(Vector) :: valuesOfIndependentVariables
    integer :: i
    ! integer :: sizeGradients !Only for debug

    call Vector_constructor( valuesOfIndependentVariables, size(evaluationPoint) )
    valuesOfIndependentVariables%values = evaluationPoint

    if( .not. CONTROL_instance%ANALYTIC_GRADIENT ) then
       gradientVector = EnergyGradients_getNumericGradient( valuesOfIndependentVariables, GeometryOptimizer_getFunctionValue )
    else
!       do i=1, this_pointer%numberOfIndependVariables
!          gradientVector(i) = EnergyGradients_getDerivative(valuesOfIndependentVariables , [i], &
!               GeometryOptimizer_getFunctionValue )
!       end do
        call EnergyGradients_getDerivative(valuesOfIndependentVariables,  GeometryOptimizer_getFunctionValue, gradientVector )
    end if

    ! Debug Mauricio Rodas
    ! sizeGradients = size(gradientVector) 
    ! write(*,"(A)") "Gradientes numericos"
    ! do i=1, sizeGradients
    !    write(*,"(f20.12)") gradientVector(i)
    ! end do
    ! print *,""
    call Vector_destructor( valuesOfIndependentVariables )

  end subroutine GeometryOptimizer_getGradient

!         subroutine GeometryOptimizer_builtInputFile( evaluationPoint )
!             real(8) :: evaluationPoint(:)

!             open( UNIT=37,FILE="tmp.aux",STATUS='REPLACE', &
!                 ACCESS='SEQUENTIAL', FORM='FORMATTED' )

!                 call Input_Parsing_writeNameList(37,"InputLowdinParameters")

!             close(37)

!         end subroutine GeometryOptimizer_builtInputFile

  !>
  !! @brief Define un proyector para el gradiente cartesiano
  !!
  !!  @warning Anque no se define como privada de la clase por propositos de conveniencia
  !!          no debe ser empleada fuera de la clase
  !<
  subroutine GeometryOptimizer_projectGradient( gradientVector, projectedGradient )
    implicit none
    real(8) :: gradientVector(:)
    real(8) :: projectedGradient(:)

    type(vector) :: pprojectedGradient


    !           if ( CONTROL_instance%PROJECT_HESSIANE ) then
    !
    !               MolecularSystem_instance(fragment)%mecProperties.transformationMatrix = &
    !                   MecanicProperties_getTransformationMatrix( MolecularSystem_instance(fragment)%mecProperties )
    !               call Vector_constructor(pprojectedGradient,size(gradientVector))
    !               pprojectedGradient%values = gradientVector
    !               pprojectedGradient =MecanicProperties_getGradientProjectedOnExternalDegrees( MolecularSystem_instance(fragment)%mecProperties, &
    !                   pprojectedGradient )
    !               projectedGradient=pprojectedGradient%values
    !               call Vector_destructor(pprojectedGradient)
    !
    !           else

    projectedGradient = gradientVector

    !           end if

  end subroutine GeometryOptimizer_projectGradient

  !>
  !! @brief Define un proyector para la hesiana cartesiano
  !!
  !!  @warning Anque no se define como privada de la clase por propositos de conveniencia
  !!          no debe ser empleada fuera de la clase
  !<
  subroutine GeometryOptimizer_projectHessiane( hessianeMatrix, projectedHessiane)
    implicit none
    real(8) :: hessianeMatrix(:,:)
    real(8) :: projectedHessiane(:,:)

    type(Matrix) :: projector
    real(8) :: auxValue
    integer :: infoProcess
    integer :: i
    integer :: j

    if ( CONTROL_instance%PROJECT_HESSIANE .or. CONTROL_instance%MOLLER_PLESSET_CORRECTION >= 2 ) then

       projector = MecanicProperties_getHandyProjector( MolecularSystem_instance%mechanicalProp, infoProcess)

       if(infoProcess==0) then
          projectedHessiane=matmul(projector%values,matmul(hessianeMatrix,projector%values))

          do i=1,size(hessianeMatrix,dim=1)
             do  j=1,i
                auxValue=0.5*( projectedHessiane(i,j) + projectedHessiane(j,i) )
                projectedHessiane(j,i) = auxValue
                projectedHessiane(i,j) = auxValue
             end do
          end do

       else
          projectedHessiane = hessianeMatrix
          print *, "HESSIANE WILL NOT PROJECTED, THE INERTIA MATRIX IS BADLY CONDITIONED"
          print *,""

       end if
       call Matrix_destructor(projector)

    else
       projectedHessiane = hessianeMatrix
    end if

  end subroutine GeometryOptimizer_projectHessiane

!         !>
!         !! @brief  Maneja excepciones de la clase
!         !<
!         subroutine GeometryOptimizer_exception( this, typeMessage, description, debugDescription)
!             implicit none
!             type(GeometryOptimizer) :: this
!             integer :: typeMessage
!             character(*) :: description
!             character(*) :: debugDescription

!             type(Exception) :: ex

!             call Exception_constructor( ex , typeMessage )
!             call Exception_setDebugDescription( ex, debugDescription )
!             call Exception_setDescription( ex, description )
!             call Exception_show( ex )
!             call Exception_destructor( ex )

!         end subroutine GeometryOptimizer_exception

end module GeometryOptimizer_
