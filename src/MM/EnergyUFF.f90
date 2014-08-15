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
  use Graph_
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
    type(Graph) :: this
    real(8) :: totalStretchingEnergy
    real(8) :: totalStretchingEnergyKJ
    real(8) :: totalBendingEnergy
    real(8) :: totalBendingEnergyKJ
    real(8) :: totalTorsionEnergy
    real(8) :: totalTorsionEnergyKJ
    real(8) :: totalVDWEnergy
    real(8) :: totalVDWEnergyKJ
    real(8) :: totalInversionEnergy
    real(8) :: totalInversionEnergyKJ
    real(8) :: totalElectrostaticEnergy
    real(8) :: totalElectrostaticEnergyKJ
    real(8) :: totalEnergy
    real(8) :: totalEnergyKJ
    integer :: i


    totalStretchingEnergy = 0.0
    totalBendingEnergy = 0.0
    totalTorsionEnergy = 0.0
    totalVDWEnergy = 0.0
    totalInversionEnergy = 0.0
    totalElectrostaticEnergy = 0.0

    do i=1, this%edges%numberOfEdges
       totalStretchingEnergy = totalStretchingEnergy + this%edges%stretchingEnergy(i)
    end do

    do i=1, this%angles%numberOfAngles
       totalBendingEnergy = totalBendingEnergy + this%angles%bendingEnergy(i)
    end do

    do i=1, this%torsions%numberOfTorsions
       totalTorsionEnergy = totalTorsionEnergy + this%torsions%torsionEnergy(i)
    end do

    do i=1, this%vdwaals%numberOfVDWaals
       totalVDWEnergy = totalVDWEnergy + this%vdwaals%VDWEnergy(i)
    end do

    totalStretchingEnergyKJ = totalStretchingEnergy*4.1868
    totalBendingEnergyKJ = totalBendingEnergy*4.1868
    totalTorsionEnergyKJ = totalTorsionEnergy*4.1868
    totalVDWEnergyKJ = totalVDWEnergy*4.1868
    totalInversionEnergyKJ = totalInversionEnergy*4.1868
    totalElectrostaticEnergyKJ = totalElectrostaticEnergy*4.1868

    totalEnergy = totalStretchingEnergy + totalBendingEnergy + totalTorsionEnergy + totalVDWEnergy + totalInversionEnergy + totalElectrostaticEnergy
    totalEnergyKJ = totalEnergy*4.1868

    write(*,"(T5,A)") ""
    write(*,"(T5,A)") ""
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T12,A)") "Summary of total energies"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T7,A,T25,A,T39,A)") "Type", "kcal/mol", "kJ/mol"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Stretching", totalStretchingEnergy, totalStretchingEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Bending", totalBendingEnergy, totalBendingEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Torsional", totalTorsionEnergy, totalTorsionEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Van der Waals", totalVDWEnergy, totalVDWEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Out of Plane", totalInversionEnergy, totalInversionEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Electrostatic", totalElectrostaticEnergy, totalElectrostaticEnergyKJ
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "TOTAL", totalEnergy, totalEnergyKJ
    write(*,"(T5,A)") "------------------------------------------"

    
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
