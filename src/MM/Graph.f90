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
!!        This module creates a graph with the whole information about the system
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions using Universal Force Field has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Graph_
  use CONTROL_
  use Matrix_
  use MatrixInteger_
  use Rings_
  use Edges_
  use Vertex_
  use Angles_
  use Torsions_
  use VDWaals_
  use Electrostatic_
  use Inversions_
  use Exception_
  implicit none
  
  type , public :: Graph
     
     type(Rings) :: rings
     type(Edges) :: edges
     type(Vertex) :: vertex
     type(Angles) :: angles
     type(Torsions) :: torsions
     type(VDWaals) :: vdwaals
     type(Electrostatic) :: electrostatic
     type(Inversions) :: inversions

  end type Graph

  public :: &
       Graph_initialize,&
       Graph_destructor

  !>Singleton
  type(Graph), public, target :: Graph_instance

contains

  !>
  !! @brief Initializes the graph calling all constructors
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] forcefield CHARACTER Force Field selected by the user, for now only is implemented UFF
  !! @param [in] electrostatic LOGICAL evaluates if the user requires Electrostatic Energy
  subroutine Graph_initialize( forcefield, electrostatic )
    implicit none
    character(50), intent(in) :: forcefield
    logical, intent(in) :: electrostatic
 
    call Rings_constructor( Graph_instance%rings )
    call Vertex_constructor( Graph_instance%vertex, forcefield, Graph_instance%rings )
    call Edges_constructor( Graph_instance%edges, Graph_instance%vertex, Graph_instance%rings )
    call Angles_constructor( Graph_instance%angles, Graph_instance%vertex, Graph_instance%edges )
    call Torsions_constructor( Graph_instance%torsions, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )
    call VDWaals_constructor( Graph_instance%vdwaals, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )
    if(electrostatic) then
       call Electrostatic_constructor(Graph_instance%electrostatic, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )
    end if
    call Inversions_constructor(Graph_instance%inversions, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles)

  end subroutine Graph_initialize

  subroutine Graph_destructor(this, electrostatic)
    implicit none
    type(Graph), intent(inout) :: this
    logical, intent(in):: electrostatic
    integer :: i
    
    ! Deallocate Rings
    if(this%rings%hasRings) then
       deallocate(this%rings%ringSize)
       deallocate(this%rings%aromaticity)
       do i=1, this%rings%numberOfRings
          call MatrixInteger_destructor(this%rings%connectionMatrix(i))
       end do
       deallocate(this%rings%connectionMatrix)
    end if
    ! Deallocate Vertex
    deallocate(this%vertex%symbol)
    deallocate(this%vertex%type)
    deallocate(this%vertex%charges)
    call Matrix_destructor(this%vertex%cartesianMatrix)
    deallocate(this%vertex%bondValence)
    deallocate(this%vertex%angleValence)
    deallocate(this%vertex%distanceVdW)
    deallocate(this%vertex%energyVdW)
    deallocate(this%vertex%scaleVdW)
    deallocate(this%vertex%effectiveCharge)
    deallocate(this%vertex%torsionalBarrier)
    deallocate(this%vertex%torsionalConstant)
    deallocate(this%vertex%electronegativityGMP)
    deallocate(this%vertex%hard)
    deallocate(this%vertex%radius)
    deallocate(this%vertex%hybridization)
    do i=1, this%vertex%numberOfVertices
       call Matrix_destructor(this%vertex%ionizationPotential(i))
    end do
    deallocate(this%vertex%ionizationPotential)
    deallocate(this%vertex%connectivity)
    ! Deallocate Edges
    call MatrixInteger_destructor(this%edges%connectionMatrix)
    deallocate(this%edges%distance)
    deallocate(this%edges%bondOrder)
    deallocate(this%edges%idealDistance)
    deallocate(this%edges%forceConstant)
    deallocate(this%edges%stretchingEnergy) 
    deallocate(this%edges%stretchingEnergyKJ)
    ! Deallocate Angles
     call MatrixInteger_destructor(this%angles%connectionMatrix)
     deallocate(this%angles%theta)
     deallocate(this%angles%idealTheta)
     deallocate(this%angles%forceConstant)
     deallocate(this%angles%cosTheta)
     deallocate(this%angles%cosIdealTheta)
     deallocate(this%angles%sinTheta)
     deallocate(this%angles%sinIdealTheta)
     deallocate(this%angles%bendingEnergy)
     deallocate(this%angles%bendingEnergyKJ)
    ! Deallocate Torsions
     if(this%torsions%hasTorsion) then
        call MatrixInteger_destructor(this%torsions%connectionMatrix)
        deallocate(this%torsions%phi)
        deallocate(this%torsions%rotationalBarrier)
        deallocate(this%torsions%idealPhi)
        deallocate(this%torsions%order)
        deallocate(this%torsions%torsionEnergy)
        deallocate(this%torsions%torsionEnergyKJ)
     end if
    ! Deallocate Inversions
     if(this%inversions%hasInversions) then
        call MatrixInteger_destructor(this%inversions%connectionMatrix)
        deallocate(this%inversions%omega)
        deallocate(this%inversions%C0)
        deallocate(this%inversions%C1)
        deallocate(this%inversions%C2)
        deallocate(this%inversions%forceConstant)
        deallocate(this%inversions%inversionEnergy)
        deallocate(this%inversions%inversionEnergyKJ)
     end if
    ! Deallocate Van der Waals
     if(this%vdwaals%VDW) then
        call MatrixInteger_destructor(this%vdwaals%connectionMatrix)
        deallocate(this%vdwaals%distance)
        deallocate(this%vdwaals%idealDistance)
        deallocate(this%vdwaals%wellDepth)
        deallocate(this%vdwaals%VDWEnergy)
        deallocate(this%vdwaals%VDWEnergyKJ)
     end if
    ! Deallocate Electrostatics
     if(electrostatic) then
        call MatrixInteger_destructor(this%electrostatic%connectionMatrix)
        deallocate(this%electrostatic%partialCharge)
        deallocate(this%electrostatic%distance)
        deallocate(this%electrostatic%electrostaticEnergy)
        deallocate(this%electrostatic%electrostaticEnergyKJ)
     end if

  end subroutine Graph_destructor

end module Graph_
