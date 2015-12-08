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
  use GSLInterface_
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
     integer :: method ! using by MR
     ! integer :: minimizationType  ! using by MR
     integer :: numberOfIterations ! using by MR
     real(8) :: toleranceGradient ! using by MR
     real(8) :: toleranceDx ! using by MR
     integer :: numberOfIndependVariables ! using by MR
     type(Vector) :: valuesOfIndependentVariables ! using by MR
     type(Vector) :: gradient
     real(8) :: functionValue
     real(8) :: stepSize ! using by MR
     logical :: isInitialExecution ! using by MR
     logical :: isInstanced ! using by MR

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
    ! type(Vector) :: initialGeometry
    integer :: i
    integer :: numberOfCenterofOptimization

    this%isInstanced =.true.

    numberOfCenterofOptimization = ParticleManager_getNumberOfCentersOfOptimization()
    this%valuesOfIndependentVariables = ParticleManager_getPositionOfCenterOfOptimizacion()
    ! initialGeometry = ParticleManager_getPositionOfCenterOfOptimizacion()

    this%name = ""
    this%isInitialExecution=.true.
    this%numberOfIndependVariables = size( this%valuesOfIndependentVariables%values )
    this%method = CONTROL_instance%MINIMIZATION_METHOD
    this%toleranceGradient =  CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT
    this%toleranceDx = CONTROL_instance%MINIMIZATION_LINE_TOLERANCE
    this%stepSize = CONTROL_instance%MINIMIZATION_INITIAL_STEP_SIZE
    this%numberOfIterations = CONTROL_instance%MINIMIZATION_MAX_ITERATION
    ! this%minimizationType = 1

    ! write(*,*) "MAURICIO RODAS: "
    ! write(*,*) "Size: ", this%numberOfIndependVariables
    ! write(*,*) "Method: ", this%method
    ! write(*,*) "TOL GRAD: ", this%toleranceGradient
    ! write(*,*) "TOL Dx: ", this%toleranceDx
    ! write(*,*) "Step size: ", this%stepSize
    ! write(*,*) "Max Iterations: ", this%numberOfIterations
    ! write(*,*) "Initial geometry: "
    ! do i=1, this%numberOfIndependVariables
    !    write(*,*) this%valuesOfIndependentVariables%values(i)
    ! end do

    call EnergyGradients_constructor()
    ! call Vector_destructor(initialGeometry)

  end subroutine GeometryOptimizer_constructor

        !**
        ! @brief Define el destructor
        !**
  subroutine GeometryOptimizer_destructor(this)
    implicit none
    type(GeometryOptimizer) :: this

    call EnergyGradients_destructor()

  end subroutine GeometryOptimizer_destructor

  !>
  !! @brief Ejecuta el minimizador hasta encontrar un minimo de la estructura
  !<
  subroutine GeometryOptimizer_run(this )
    implicit none
    type(GeometryOptimizer), target :: this
    real(8), allocatable :: coordinates(:)
    integer :: infoError
    real(8) :: energy
    type(Vector) :: geometry
    real(8) :: auxValue
    logical :: isMinumum
    type(Matrix) :: initialHessian
    type(Hessians) :: hessiansInfo
    logical :: lastStep
    logical :: firstStep

    allocate(coordinates(this%numberOfIndependVariables))

    coordinates = this%valuesOfIndependentVariables%values

    print   *,""
    print *, " BEGIN GEOMETRY OPTIMIZATION: "
    print *,"------------------------------------------------------------"
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

    infoError = tolow_minimize(this%method, this%numberOfIndependVariables, coordinates(1),&
         this%stepSize, this%toleranceGradient, this%toleranceDx, this%numberOfIterations, &
         GeometryOptimizer_calculatePoint, GeometryOptimizer_showIterationInfo, energy)

    if (infoError /= 0) then
       write(6,*) "Error occurred during the GSL minimization procedure"
       stop
    else
       write(6,*) ""
       write (6,"(T20,A30)") " FINAL GEOMETRY: AMSTRONG"
       write (6,"(T18,A35)") "------------------------------------------"
       write(6,*) ""

       call Vector_constructor( geometry, size(coordinates) )
       geometry%values = coordinates

       !! Ajusta el origen de las particulas presentes en el sistema
       call ParticleManager_setParticlesPositions(geometry)
       call Vector_destructor( geometry )
       call MolecularSystem_showCartesianMatrix()
       call MolecularSystem_showDistanceMatrix()
       call MolecularSystem_saveToFile()
       ! call MolecularSystem_showZMatrix( MolecularSystem_instance )
       CONTROL_instance%OPTIMIZE=.false.
       CONTROL_instance%SCF_CONVERGENCE_CRITERIUM="density"
       lastStep = .true.

       open(unit=40, file="lowdin.dat", status="replace", form="formatted")

       !!save all options
       call CONTROL_save(40, lastStep)

       close(40)

       print *,""
       print *,"END GEOMETRY OPTIMIZATION "
       print *,""

       call Solver_run()

    end if

    ! !! Realiza un calculo final en la geometria de equilibrio

  end subroutine GeometryOptimizer_run

  subroutine GeometryOptimizer_calculatePoint(size, coordinates, functionValue, getgrad, gradients)
    integer, intent(in)    :: size
    real(8), intent(in)    :: coordinates(size)
    real(8), intent(inout) :: functionValue
    integer, intent(in)    :: getgrad
    real(8), intent(inout) :: gradients(size)


    ! write(6,*) "Entre a la funcion"
    functionValue = GeometryOptimizer_getFunctionValue(coordinates)
    ! write(6,*) "Energia: ", functionValue

    if(getgrad .eq. 1) then
       ! write(6,*) "Voy a calcular gradientes"
       call GeometryOptimizer_getGradient( coordinates, gradients )
       ! write(6,*) "Saliendo de calcular gradientes"
    end if

  end subroutine GeometryOptimizer_calculatePoint

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
    logical :: lastStep

    type(Vector) :: valuesOfIndependentVariables

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    call Vector_constructor( valuesOfIndependentVariables, size(evaluationPoint) )
    valuesOfIndependentVariables%values = evaluationPoint

    !! Ajusta el origen de las particulas presentes en el sistema
    call ParticleManager_setParticlesPositions( valuesOfIndependentVariables)
    call Vector_destructor( valuesOfIndependentVariables )
    call MolecularSystem_saveToFile()

    lastStep = .false.

    open(unit=40, file="lowdin.dat", status="replace", form="formatted")
    
    !!save all options
    call CONTROL_save(40, lastStep)

    close(40)
    
    CONTROL_instance%SCF_CONVERGENCE_CRITERIUM="energy"

    call Solver_run()
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])
    close(wfnUnit)
   
    output = totalEnergy

    ! this_pointer%isInitialExecution=.false.

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

    ! write(*,"(A)") "Dentro de GeometryOptimizer_getGradient"
    ! do i=1, size(evaluationPoint)
    !    write(*,"(F17.12)") valuesOfIndependentVariables%values(i)
    ! end do

    if( .not. CONTROL_instance%ANALYTIC_GRADIENT ) then
       gradientVector = EnergyGradients_getNumericGradient( valuesOfIndependentVariables, GeometryOptimizer_getFunctionValue )
    else
       call EnergyGradients_getDerivative(valuesOfIndependentVariables,  GeometryOptimizer_getFunctionValue, gradientVector )
    end if

    ! Debug Mauricio Rodas
    ! sizeGradients = size(gradientVector) 
    ! write(*,"(A)") "Gradientes"
    ! do i=1, sizeGradients
    !    write(*,"(f20.12)") gradientVector(i)
    ! end do
    ! print *,""

    call Vector_destructor( valuesOfIndependentVariables )

  end subroutine GeometryOptimizer_getGradient

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
    do i = 1, totalNumberOfParticles

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
  !! @brief Muestra en la informacion en un punto de iteracion
  !!
  !<
  subroutine GeometryOptimizer_showIterationInfo(iterationPoint, size, energy, maxdx, maxdf, coordinates)
    implicit none
    integer, intent(in) :: iterationPoint
    integer, intent(in) :: size
    real(8), intent(in) :: energy, maxdx, maxdf
    real(8), intent(in) :: coordinates(size)

    ! real(8) :: evaluationPoint(:)
    ! real(8) :: gradient(:)
    ! real(8) :: functionValue

    integer :: i
    ! integer :: k
    ! real(8) :: origin(3)
    ! integer :: totalNumberOfParticles

    ! totalNumberOfParticles = size( ParticleManager_instance )

    write(6,"(T20,A65)") "-----------------------------------------------------------------"
    write(6,"(T20,A29,I4)") "GEOMETRY OPTIMIZATION POINT: ", iterationPoint
    write(6,"(T20,A65)") "-----------------------------------------------------------------"
    write(6,"(T20,A9,F20.10)") "ENERGY = ", energy
    write(6,"(T20,A22,F20.10)") "MAX POSITION CHANGE = ", maxdx
    write(6,"(T20,A22,F20.10)") "MAX GRADIENT CHANGE = ", maxdf
    write(6,"(T20,A65)") "-----------------------------------------------------------------"
    write(6,*) ""
    ! write(6,"(T20,A25)") "COORDINATES: "//trim(CONTROL_instance%UNITS)
    ! do i = 1, size
    !    write(6,"(T20,F20.10)") coordinates(i)
    ! end do
!     write (6,"(A10,A16,A20,A20)") "Particle","<x>","<y>","<z>"
!     do i = 1, totalNumberOfParticles!ParticleManager_getTotalNumberOfParticles()

!        if ( ParticleManager_isCenterOfOptimization( i ) ) then
!           origin = ParticleManager_getOrigin( iterator = i )
!           if ( CONTROL_instance%UNITS=="ANGSTROMS") then
!              origin = origin * AMSTRONG
!           end if
! #ifdef intel
!           write (6,"(A10,<3>F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), origin(1), origin(2), origin(3)
! #else
!           write (6,"(A10,3F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), origin(1), origin(2), origin(3)
! #endif
!        end if

!     end do


!     ! print *,""
!     write (6,"(A)") ""
!     write (6,"(T20,A25)") "GRADIENT: HARTREE/BOHR"
!     write (6,"(A10,A16,A20,A20)") "Particle","dE/dx","dE/dy","dE/dz"
!     k=1
!     do i = 1, totalNumberOfParticles

!        if ( ParticleManager_isCenterOfOptimization( i ) .and. k < size(gradient)) then

!           write (6,"(A10,3F20.10)") trim( ParticleManager_getSymbol( iterator = i ) ), &
!                gradient(k), gradient(k+1),gradient(k+2)
!           k=k+3

!        end if

!     end do
!     write (6,*) ""
!     write (6,"(T10,A24,F20.10)") "TOTAL ENERGY (Hartree) =", functionValue

  end subroutine GeometryOptimizer_showIterationInfo

  !>
  !! @brief Define un proyector para el gradiente cartesiano
  !!
  !!  @warning Anque no se define como privada de la clase por propositos de conveniencia
  !!          no debe ser empleada fuera de la clase
  !<
  ! subroutine GeometryOptimizer_projectGradient( gradientVector, projectedGradient )
  !   implicit none
  !   real(8) :: gradientVector(:)
  !   real(8) :: projectedGradient(:)

  !   type(vector) :: pprojectedGradient


  !   !           if ( CONTROL_instance%PROJECT_HESSIANE ) then
  !   !
  !   !               MolecularSystem_instance(fragment)%mecProperties.transformationMatrix = &
  !   !                   MecanicProperties_getTransformationMatrix( MolecularSystem_instance(fragment)%mecProperties )
  !   !               call Vector_constructor(pprojectedGradient,size(gradientVector))
  !   !               pprojectedGradient%values = gradientVector
  !   !               pprojectedGradient =MecanicProperties_getGradientProjectedOnExternalDegrees( MolecularSystem_instance(fragment)%mecProperties, &
  !   !                   pprojectedGradient )
  !   !               projectedGradient=pprojectedGradient%values
  !   !               call Vector_destructor(pprojectedGradient)
  !   !
  !   !           else

  !   projectedGradient = gradientVector

  !   !           end if

  ! end subroutine GeometryOptimizer_projectGradient

  ! !>
  ! !! @brief Define un proyector para la hesiana cartesiano
  ! !!
  ! !!  @warning Anque no se define como privada de la clase por propositos de conveniencia
  ! !!          no debe ser empleada fuera de la clase
  ! !<
  ! subroutine GeometryOptimizer_projectHessiane( hessianeMatrix, projectedHessiane)
  !   implicit none
  !   real(8) :: hessianeMatrix(:,:)
  !   real(8) :: projectedHessiane(:,:)

  !   type(Matrix) :: projector
  !   real(8) :: auxValue
  !   integer :: infoProcess
  !   integer :: i
  !   integer :: j

  !   if ( CONTROL_instance%PROJECT_HESSIANE .or. CONTROL_instance%MOLLER_PLESSET_CORRECTION >= 2 ) then

  !      projector = MecanicProperties_getHandyProjector( MolecularSystem_instance%mechanicalProp, infoProcess)

  !      if(infoProcess==0) then
  !         projectedHessiane=matmul(projector%values,matmul(hessianeMatrix,projector%values))

  !         do i=1,size(hessianeMatrix,dim=1)
  !            do  j=1,i
  !               auxValue=0.5*( projectedHessiane(i,j) + projectedHessiane(j,i) )
  !               projectedHessiane(j,i) = auxValue
  !               projectedHessiane(i,j) = auxValue
  !            end do
  !         end do

  !      else
  !         projectedHessiane = hessianeMatrix
  !         print *, "HESSIANE WILL NOT PROJECTED, THE INERTIA MATRIX IS BADLY CONDITIONED"
  !         print *,""

  !      end if
  !      call Matrix_destructor(projector)

  !   else
  !      projectedHessiane = hessianeMatrix
  !   end if

  ! end subroutine GeometryOptimizer_projectHessiane

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
