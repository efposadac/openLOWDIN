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
!! @brief This module handles all molecular system
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-14
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapts module to inclusion on LOWDIN package
!!   - <tt> 2013-04-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrites the module to avoid XML dependence and for implement new LOWDIN standard.
!!
module Graph_
  use CONTROL_
  use Edges_
  use Vertex_
  use Angles_
  use Torsions_
  use VDWaals_
  use Electrostatic_
  use Exception_
  implicit none
  
  type , public :: Graph
     
     type(Edges) :: edges
     type(Vertex) :: vertex
     type(Angles) :: angles
     type(Torsions) :: torsions
     type(VDWaals) :: vdwaals
     type(Electrostatic) :: electrostatic

  end type Graph

  public :: &
       Graph_initialize

  !>Singleton
  type(Graph), public, target :: Graph_instance

contains

  subroutine Graph_initialize( forcefield )
    implicit none
    character(50), intent(in) :: forcefield

    call Vertex_constructor( Graph_instance%vertex, forcefield )
    call Edges_constructor( Graph_instance%edges, Graph_instance%vertex )
    call Angles_constructor( Graph_instance%angles, Graph_instance%vertex, Graph_instance%edges )
    call Torsions_constructor( Graph_instance%torsions, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )
    call VDWaals_constructor( Graph_instance%vdwaals, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )
    call Electrostatic_constructor(Graph_instance%electrostatic, Graph_instance%vertex, Graph_instance%edges, Graph_instance%angles )

  end subroutine Graph_initialize

end module Graph_
