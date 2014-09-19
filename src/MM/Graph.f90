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
       Graph_initialize

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

end module Graph_
