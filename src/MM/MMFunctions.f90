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
module MMFunctions_
  use CONTROL_
  use Graph_
  ! use ParticleManager_
  ! use Edges_
  ! use AtomTypeUFF_
  ! use EnergyUFF_
  use Exception_
  implicit none


	type :: MolecularMechanics
    
       	character(50) :: ffmethod
        logical :: isInstanced

	end type MolecularMechanics

       type(MolecularMechanics), target :: MolecularMechanics_instance

       public :: &
            MolecularMechanics_constructor, &
            MolecularMechanics_destructor, &
            MolecularMechanics_run

contains
	!**
	! Define el constructor para la clase
	!
	!**
  subroutine MolecularMechanics_constructor()
    implicit none
  
    MolecularMechanics_instance%isInstanced =.true.

  end subroutine MolecularMechanics_constructor

	!**
	! Define el destructor para clase
	!
	!**
  subroutine MolecularMechanics_destructor()
    implicit none

    MolecularMechanics_instance%isInstanced =.false.

  end subroutine MolecularMechanics_destructor

  subroutine MolecularMechanics_run( ffmethod )
    implicit none
    character(50), intent(in) :: ffmethod
    type(Exception) :: ex
!! Parametros para impresion borrar luego
    integer :: atomAIdx, AtomBIdx
    character(10) :: atomA, AtomB
    integer :: i

    
    MolecularMechanics_instance%ffmethod = ffmethod

    if ( MolecularMechanics_instance%isInstanced ) then
    
       if ( MolecularMechanics_instance%ffmethod == "UFF" ) then
          call Graph_initialize(MolecularMechanics_instance%ffmethod)

          write(*,"(T20,A)") ""
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write(*,"(T60,A)") "INITIAL GEOMETRY: AMSTRONG"
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write (*,"(T20,A,T30,A,T40,A,T50,A,T68,A,T83,A,T98,A)") "Idx", &
               "Atom", "Type", "Charge(Z)", &
               "<x>","<y>","<z>"
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          do i=1,Graph_instance%vertex%numberOfVertices
                write(*,"(T10,I,T30,A,T40,A,T50,F8.5,T60,F12.5,T75,F12.5,T90,F12.5)") i, &
                     trim(Graph_instance%vertex%symbol(i)), &
                     trim( Graph_instance%vertex%type(i) ), &
                     Graph_instance%vertex%charges(i), &
                     Graph_instance%vertex%cartesianMatrix%values(i,1), &
                     Graph_instance%vertex%cartesianMatrix%values(i,2), &
                     Graph_instance%vertex%cartesianMatrix%values(i,3)
          end do
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write(*,"(T20,A)") ""

          write(*,"(T20,A)") ""
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write(*,"(T60,A)") "UFF Parameters"
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write (*,"(T20,A,T30,A,T40,A,T50,A,T61,A,T78,A,T91,A)") "Idx", &
               "Atom", "Type", "Charge(Z)", &
               "Eff. Charge(Z*)","Bond Valence","Angle Valence"
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          do i=1,Graph_instance%vertex%numberOfVertices
                write(*,"(T10,I,T30,A,T40,A,T50,F8.5,T60,F12.5,T75,F12.5,T90,F12.5)") i, &
                     trim(Graph_instance%vertex%symbol(i)), &
                     trim( Graph_instance%vertex%type(i) ), &
                     Graph_instance%vertex%charges(i), &
                     Graph_instance%vertex%effectiveCharge(i), &
                     Graph_instance%vertex%bondValence(i), &
                     Graph_instance%vertex%angleValence(i)
          end do
          write(*,"(T20,A)") "----------------------------------------------------------------------------------------------------"
          write(*,"(T20,A)") ""

          
          write(*,"(T20,A)") ""
          write(*,"(T20,A)") "--------------------------------------------------------------------------------------------"
          write(*,"(T20,A)") "                Informacion completa de enlaces "
          write(*,"(T20,A)") "--------------------------------------------------------------------------------------------"
          do i=1,Graph_instance%edges%numberOfEdges
             atomAIdx=Graph_instance%edges%connectionMatrix%values(i,1)
             Write( atomA, '(i10)' ) atomAIdx
             atomA = adjustl(trim(atomA))
             atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
             atomBIdx=Graph_instance%edges%connectionMatrix%values(i,2)
             Write( atomB, '(i10)' ) atomBIdx
             atomB = adjustl(trim(atomB))
             atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
             write(*,"(T20,I5,2x,2A,2x,F8.5,2x,F8.5)") i, atomA, atomB, &
                  Graph_instance%edges%bondOrder%values(i), Graph_instance%edges%distance%values(i)
          end do
          write(*,"(T20,A)") "--------------------------------------------------------------------------------------------"
          write(*,"(T20,A)") ""

          ! write(*,"(T20,A,I,A)") "Voy a construir un grafo con ", Graph_instance%numberOfVertex, " vertices" 
          ! call Edges_getBondOrders(bondOrders)
          ! call AtomTypeUFF_run(ffAtomType)
          ! call EnergyUFF_run(ffAtomType, bondOrders)
       else
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object MolecularMechanics in run() function" )
          call Exception_setDescription( ex, "This Force Field hasn't been implemented" )
          call Exception_show( ex )
       end if
    else
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object MolecularMechanics in run() function" )
       call Exception_setDescription( ex, "You should to instance MolecularMechanics module before use this function" )
       call Exception_show( ex )
    end if

  end subroutine MolecularMechanics_run

subroutine MolecularMechanics_exception( typeMessage, description, debugDescription)
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

end subroutine MolecularMechanics_exception

  
end module MMFunctions_
