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
module Electrostatic_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Vertex_
  use Edges_
  use Angles_
  use MMCommons_
  use Matrix_
  use Vector_
  use ChargesEQeq_
  use Exception_
  implicit none

  type , public :: Electrostatic
     
     real(8), allocatable :: partialCharge(:)
     real(8), allocatable :: electrostaticEnergy(:) !! Kcal/mol
     real(8), allocatable :: electrostaticEnergyKJ(:) !! KJ/mol

  end type Electrostatic


       public :: &
            Electrostatic_constructor

contains

  subroutine Electrostatic_constructor( this, vertices )
    implicit none
    type(Electrostatic), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    real(8), allocatable :: partials(:)

    call ChargesEQeq_getCharges(partials, vertices)
    allocate( this%partialCharge( vertices%numberOfVertices ) )

    do i=1,vertices%numberOfVertices
       this%partialCharge(i) = partials(i)
    end do

  end subroutine Electrostatic_constructor

  subroutine Electrostatic_exception( typeMessage, description, debugDescription)
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

  end subroutine Electrostatic_exception

end module Electrostatic_
