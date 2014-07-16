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
!! @brief Module for atomic elements definitions
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-05
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapto al estandar de codificacion propuesto.
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.cn
!!        -# Reescribe y adapta el modulo  para su inclusion en Lowdin
!!        -# Elminates XML dependence. The module is rewritten.
module EnergyUFF_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use UFFParameters_
  use Exception_
  implicit none
  
  
  public :: &
       EnergyUFF_run
       ! EnergyUFF_show

contains
    
  !>
  !! @brief Loads an atomic element from library.
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine EnergyUFF_run( this )
    implicit none
    character(10), allocatable, intent(in) :: this(:)
    integer :: ffAtomTypeSize, i
    type(UFFParameters) :: atomType
    character(10) :: type

          ffAtomTypeSize = size(this)
          do i=1,ffAtomTypeSize
             type = trim(this(i))
             call UFFParameters_load( atomType, trim(type) )
          end do

  end subroutine EnergyUFF_run
  
  !<
  !! @brief Define el destructor para clase
  ! subroutine EnergyUFF_show( this )
  !   implicit none
  !   type(EnergyUFF) , intent(in) :: this

    
  !   print *,""
  !   print *,"---------------------------------------------------"
  !   print *,"  Atom Type Parameters   "
  !   print *,"---------------------------------------------------"
  !   print *,""
  !   write (6,"(T10,A22,A12,A10)")	"Type                   = ",this%type,""
  !   write (6,"(T10,A22,F12.5,A10)")	"Bond                   = ",this%bond," Angstroms"
  !   write (6,"(T10,A22,F12.5,A10)")	"Angle                  = ",this%angle," Degrees"
  !   write (6,"(T10,A22,F12.5,A10)")	"VdW distance           = ",this%distanceVdW," Angstroms"
  !   write (6,"(T10,A22,F12.5,A10)")	"VdW energy             = ",this%energyVdW," Kcal/mol"
  !   write (6,"(T10,A22,F12.5,A10)")	"VdW scale              = ",this%scaleVdW,""
  !   write (6,"(T10,A22,F12.5,A10)")	"Effective Charge       = ",this%effectiveCharge,""
  !   write (6,"(T10,A22,F12.5,A10)")	"Torsional Barrier      = ",this%torsionalBarrier," kcal/mol"
  !   write (6,"(T10,A22,F12.5,A10)")	"Torsional Constant     = ",this%torsionalConstant," kcal/mol"
  !   write (6,"(T10,A22,F12.5,A10)")	"Electronegativity GMP  = ",this%electronegativityGMP," (pauling)"
  !   write (6,"(T10,A22,F12.5,A10)")	"Hard                   = ",this%hard,""
  !   write (6,"(T10,A22,F12.5,A10)")	"Radius                 = ",this%radius," Angstroms"
  !   print *,""
    
  ! end subroutine EnergyUFF_show

  ! function EnergyUFF_getCovalentRadius( symbolOfElement ) result( output )
  !       	implicit none
  !       	character(*),  intent( in ) :: symbolOfElement
  !       	real(8)  :: output

  !       	type(EnergyUFF) :: element
  !       	character(10) :: auxSymbol

  !       	! call EnergyUFF_constructor( element )
  !       	auxSymbol=trim( symbolOfElement )

  !       	if (  scan( auxSymbol, "_" ) /= 0 ) then
  !       		auxSymbol = trim(auxSymbol(1: scan( auxSymbol, "_" ) - 1 ) )
  !       	end if

  !       	call EnergyUFF_load( element, auxSymbol, 0 )
  !       	output = element%covalentRadius

  !       	! call EnergyUFF_destructor(element)

  ! end function EnergyUFF_getCovalentRadius

  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine EnergyUFF_exception( typeMessage, description, debugDescription)
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
    
  end subroutine EnergyUFF_exception
  
end module EnergyUFF_
