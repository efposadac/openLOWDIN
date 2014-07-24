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
module UFFParameters_
  use CONTROL_
  use Exception_
  implicit none
  
  type , public :: UFFParameters
     character(30) :: type
     real(8) :: bond
     real(8) :: angle
     real(8) :: distanceVdW
     real(8) :: energyVdW
     real(8) :: scaleVdW
     real(8) :: effectiveCharge
     real(8) :: torsionalBarrier
     real(8) :: torsionalConstant
     real(8) :: electronegativityGMP
     real(8) :: hard
     real(8) :: radius
     logical :: isInstanced
  end type UFFParameters
  
  public :: &
       UFFParameters_load, &
       UFFParameters_show
       ! UFFParameters_getCovalentRadius

contains
    
  !>
  !! @brief Loads an atomic element from library.
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine UFFParameters_load( this, typeSelected )
    implicit none
    
    type(UFFParameters), intent(inout) :: this
    character(*) :: typeSelected
    logical :: existFile
    integer :: stat
    integer :: i
    
    !! Namelist definition
    character(30) :: type
    real(8) :: bond
    real(8) :: angle
    real(8) :: distanceVdW
    real(8) :: energyVdW
    real(8) :: scaleVdW
    real(8) :: effectiveCharge
    real(8) :: torsionalBarrier
    real(8) :: torsionalConstant
    real(8) :: electronegativityGMP
    real(8) :: hard
    real(8) :: radius


    NAMELIST /atomtype/ &
         type, &
         bond, &
         angle, &
         distanceVdW, &
         energyVdW, &
         scaleVdW, &
         effectiveCharge, &
         torsionalBarrier, &
         torsionalConstant, &
         electronegativityGMP, &
         hard, &
         radius

    !! Looking for library    
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%UFF_PARAMETERS_DATABASE), exist=existFile)
    
    if ( existFile ) then
       
       !! Open library
       open(unit=10, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%UFF_PARAMETERS_DATABASE), status="old", form="formatted" )
              
       !! Read information
       type = "NONE"
       stat = 0

       rewind(10)
       
       do while(trim(type) /= trim(typeSelected))
       
          !! Setting defaults
          bond = 0
          angle = 0
          distanceVdW = 0
          energyVdW = 0
          scaleVdW = 0
          effectiveCharge = 0
          torsionalBarrier = 0
          torsionalConstant = 0
          electronegativityGMP = 0
          hard = 0
          radius = 0
         
          if (stat == -1 ) then
             
             this%isInstanced = .false.
             return
             
          end if
          
          read(10,NML=atomtype, iostat=stat)
 
          if (stat > 0 ) then
             
             print*, "ERROR!!! ", stat
             call UFFParameters_exception( ERROR, "Failed reading uffParameters.lib file!! please check this file.", "In UFFParameters at load function.")
             
          end if

       end do
       
       !! Set object variables
       this%type = type
       this%bond = bond
       this%angle = angle
       this%distanceVdW = distanceVdW
       this%energyVdW = energyVdW
       this%scaleVdW = scaleVdW
       this%effectiveCharge = effectiveCharge
       this%torsionalBarrier = torsionalBarrier
       this%torsionalConstant = torsionalConstant
       this%electronegativityGMP = electronegativityGMP
       this%hard = hard
       this%radius = radius
       this%isInstanced = .true.
       
       
       !! Debug information.
       ! call UFFParameters_show(this)
       
       close(10)
       
    else 

       call UFFParameters_exception( ERROR, "LOWDIN library not found!! please export lowdinvars.sh file.", "In UFFParameters at load function.")

    end if 

    !! Done
    
  end subroutine UFFParameters_load
  
  !<
  !! @brief Define el destructor para clase
  subroutine UFFParameters_show( this )
    implicit none
    type(UFFParameters) , intent(in) :: this

    
    print *,""
    print *,"---------------------------------------------------"
    print *,"  Atom Type Parameters   "
    print *,"---------------------------------------------------"
    print *,""
    write (6,"(T10,A22,A12,A10)")	"Type                   = ",this%type,""
    write (6,"(T10,A22,F12.5,A10)")	"Bond                   = ",this%bond," Angstroms"
    write (6,"(T10,A22,F12.5,A10)")	"Angle                  = ",this%angle," Degrees"
    write (6,"(T10,A22,F12.5,A10)")	"VdW distance           = ",this%distanceVdW," Angstroms"
    write (6,"(T10,A22,F12.5,A10)")	"VdW energy             = ",this%energyVdW," Kcal/mol"
    write (6,"(T10,A22,F12.5,A10)")	"VdW scale              = ",this%scaleVdW,""
    write (6,"(T10,A22,F12.5,A10)")	"Effective Charge       = ",this%effectiveCharge,""
    write (6,"(T10,A22,F12.5,A10)")	"Torsional Barrier      = ",this%torsionalBarrier," kcal/mol"
    write (6,"(T10,A22,F12.5,A10)")	"Torsional Constant     = ",this%torsionalConstant," kcal/mol"
    write (6,"(T10,A22,F12.5,A10)")	"Electronegativity GMP  = ",this%electronegativityGMP," (pauling)"
    write (6,"(T10,A22,F12.5,A10)")	"Hard                   = ",this%hard,""
    write (6,"(T10,A22,F12.5,A10)")	"Radius                 = ",this%radius," Angstroms"
    print *,""
    
  end subroutine UFFParameters_show

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine UFFParameters_exception( typeMessage, description, debugDescription)
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
    
  end subroutine UFFParameters_exception
  
end module UFFParameters_
