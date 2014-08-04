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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Vertex_
  use CONTROL_
  use MolecularSystem_
  use Matrix_
  use MMCommons_
  use MatrixInteger_
  use AtomTypeUFF_
  use ParticleManager_
  use UFFParameters_
  use Exception_
  implicit none

  type , public :: Vertex

     integer :: numberOfVertices
     character(10), allocatable :: symbol(:)
     character(10), allocatable :: type(:)
     real(8), allocatable :: charges(:)
     type(Matrix) :: cartesianMatrix
     real(8), allocatable :: bondValence(:)
     real(8), allocatable :: angleValence(:)
     real(8), allocatable :: distanceVdW(:)
     real(8), allocatable :: energyVdW(:)
     real(8), allocatable :: scaleVdW(:)
     real(8), allocatable :: effectiveCharge(:)
     real(8), allocatable :: torsionalBarrier(:)
     real(8), allocatable :: torsionalConstant(:)
     real(8), allocatable :: electronegativityGMP(:)
     real(8), allocatable :: hard(:)
     real(8), allocatable :: radius(:)
     integer, allocatable :: connectivity(:)

  end type Vertex


       public :: &
            Vertex_constructor

contains

  subroutine Vertex_constructor( this, forcefield )
    implicit none
    type(Vertex) :: this
    character(50), intent(in) :: forcefield
    character(10), allocatable :: ffAtomType(:)
    type(UFFParameters) :: atomType
    integer :: i
    type(Exception) :: ex
    type(MatrixInteger) :: connectivityMatrix

    call MMCommons_constructor( MolecularSystem_instance )

    this%numberOfVertices = ParticleManager_getNumberOfCentersOfOptimization()
    allocate( this%symbol( this%numberOfVertices ) )
    this%symbol = ParticleManager_getLabelsOfCentersOfOptimization()
    allocate( this%type( this%numberOfVertices ) ) 
    allocate( this%charges( this%numberOfVertices ) )
    this%charges = ParticleManager_getChargesOfCentersOfOptimization()

    allocate( this%bondValence( this%numberOfVertices ) ) 
    allocate( this%angleValence( this%numberOfVertices ) )
    allocate( this%distanceVdW( this%numberOfVertices ) )
    allocate( this%energyVdW( this%numberOfVertices ) )
    allocate( this%scaleVdW( this%numberOfVertices ) )
    allocate( this%effectiveCharge( this%numberOfVertices ) )
    allocate( this%torsionalBarrier( this%numberOfVertices ) )
    allocate( this%torsionalConstant( this%numberOfVertices ) )
    allocate( this%electronegativityGMP( this%numberOfVertices ) )
    allocate( this%hard( this%numberOfVertices ) )
    allocate( this%radius( this%numberOfVertices ) )
    allocate( this%connectivity( this%numberOfVertices ) )

    this%cartesianMatrix = ParticleManager_getCartesianMatrixOfCentersOfOptimization()

    call MatrixInteger_constructor( connectivityMatrix, this%numberOfVertices, 2 )
    call MMCommons_getConnectivityMatrix( MolecularSystem_instance, this%numberOfVertices, connectivityMatrix )

    if ( forcefield == "UFF" ) then
       call AtomTypeUFF_run(ffAtomType)
       do i=1,this%numberOfVertices
          this%type(i) = ffAtomType(i)
          call UFFParameters_load( atomType, trim(this%type(i)) )
          this%bondValence(i) = atomType%bond
          this%angleValence(i) = atomType%angle
          this%distanceVdW(i) = atomType%distanceVdW
          this%energyVdW(i) = atomType%energyVdW
          this%scaleVdW(i) = atomType%scaleVdW
          this%effectiveCharge(i) = atomType%effectiveCharge
          this%torsionalBarrier(i) = atomType%torsionalBarrier
          this%torsionalConstant(i) = atomType%torsionalConstant
          this%electronegativityGMP(i) = atomType%electronegativityGMP
          this%hard(i) = atomType%hard
          this%radius(i) = atomType%radius
          this%connectivity(i) = connectivityMatrix%values(i,2)
       end do
    else
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object Vertex in constructor() function" )
       call Exception_setDescription( ex, "This Force Field hasn't been implemented" )
       call Exception_show( ex )
    end if

  end subroutine Vertex_constructor

  subroutine Vertex_exception( typeMessage, description, debugDescription)
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

  end subroutine Vertex_exception

end module Vertex_
