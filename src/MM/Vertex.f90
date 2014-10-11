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
!! @brief Molecular Mechanics program.
!!        This module creates a class with the information of the vertices in the system
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Vertex_
  use CONTROL_
  use MolecularSystem_
  use Matrix_
  use Rings_
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
     integer, allocatable :: hybridization(:)
     type(Matrix), allocatable :: ionizationPotential(:) 
     integer, allocatable :: connectivity(:)

  end type Vertex


       public :: &
            Vertex_constructor

contains

  !>
  !! @brief Defines the class constructor, call the atom typing module a charge the parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the vertices
  !! @param [in] forcefield CHARACTER Force Field selected by the user
  !! @param [in] ring Class with the information of the rings
  !! @param numberOfVertices INTEGER number of vertices in the system
  !! @param symbol CHARACTER ARRAY symbols of the vertices
  !! @param type CHARACTER ARRAY types of atoms, those are specific for the force field
  !! @param charges REAL(8) ARRAY Nuclear charges of all atoms
  !! @param cartesianMatrix REAL(8) ARRAY cartesian coordinates for atoms
  !! @param bondValence REAL(8) ARRAY valence radius for all atoms (UFF parameter)
  !! @param angleValence REAL(8) ARRAY ideal angle for all atoms (UFF parameter)
  !! @param distanceVdW REAL(8) ARRAY Van der Waals distance for all atoms (UFF parameter)
  !! @param energyVdW REAL(8) ARRAY Van der Waals energy for all atoms (UFF parameter)
  !! @param scaleVdW REAL(8) ARRAY Van der Waals scaling factor for all atoms (UFF parameter)
  !! @param effectiveCharge REAL(8) ARRAY Nuclear effective charges for all atoms (UFF parameter)
  !! @param torsionalBarrier REAL(8) ARRAY torsional energy barrier for all atoms (UFF parameter)
  !! @param torsionalConstant REAL(8) ARRAY torsional constant for all atoms (UFF parameter)
  !! @param electronegativityGMP REAL(8) ARRAY electronegativity for all atoms (UFF parameter)
  !! @param hard REAL(8) ARRAY hardness for all atoms (UFF parameter)
  !! @param radius REAL(8) ARRAY radius for all atoms (UFF parameter)
  !! @param hybridization REAL(8) ARRAY hybridization for all atoms (UFF parameter)
  !! @param ionizationPotential REAL(8) ARRAY experimental Ionization Potential for atoms
  !! @param connectivity REAL(8) ARRAY connectivity for all atoms (UFF parameter)
  !! @see atomtypeuff_::atomtypeuff_run
  !! @see uffparameters_::uffparameters_load
  subroutine Vertex_constructor( this, forcefield, ring )
    implicit none
    type(Vertex), intent(in out) :: this
    character(50), intent(in) :: forcefield
    type(Rings), intent(in) :: ring
    character(10), allocatable :: ffAtomType(:)
    type(UFFParameters) :: atomType
    integer :: i, j
    type(Exception) :: ex
    type(MatrixInteger) :: connectivityMatrix
    type(Matrix) :: auxIonizationPotentials
    integer, allocatable :: ionizationSize(:)
    integer(8) :: size1, size2


    call MMCommons_constructor( MolecularSystem_instance )
    
    this%numberOfVertices = ParticleManager_getNumberOfCentersOfOptimization()
    size1 = this%numberOfVertices
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
    allocate( this%hybridization( this%numberOfVertices ) )
    allocate( this%connectivity( this%numberOfVertices ) )

    
    call Matrix_constructor( auxIonizationPotentials, size1, 9 )

    this%cartesianMatrix = ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    
    write(*,"(A)") "Cartesianas dentro del MM"
    do i=1, this%numberOfVertices
       write(*,"(3(F12.5))") this%cartesianMatrix%values(i,:)
    end do
    call MatrixInteger_constructor( connectivityMatrix, this%numberOfVertices, 2 )
    call MMCommons_getConnectivityMatrix( MolecularSystem_instance, this%numberOfVertices, connectivityMatrix )

    if ( forcefield == "UFF" ) then
       call AtomTypeUFF_run(ffAtomType, ring)
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
          this%hybridization(i) = atomType%hybridization
          auxIonizationPotentials%values(i,:) = atomType%ionizationPotential(:)
          this%connectivity(i) = connectivityMatrix%values(i,2)
       end do
       allocate( ionizationSize ( this%numberOfVertices )) 
       do i=1,this%numberOfVertices
          ionizationSize(i) = 9
          do j=9,2,-1
             if(auxIonizationPotentials%values(i,j) == 0.0) then
                ionizationSize(i) = ionizationSize(i) - 1
             end if
          end do
       end do

       allocate( this%ionizationPotential( this%numberOfVertices ))

       do i=1,this%numberOfVertices
          size2 = ionizationSize(i)
          call Matrix_constructor( this%ionizationPotential(i), 1, size2 )
          do j=1,ionizationSize(i)
             this%ionizationPotential(i)%values(1,j) = auxIonizationPotentials%values(i,j)
          end do
       end do
    else
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object Vertex in constructor() function" )
       call Exception_setDescription( ex, "This Force Field hasn't been implemented" )
       call Exception_show( ex )
    end if

  end subroutine Vertex_constructor

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
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
