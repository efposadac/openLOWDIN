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
!!        This module loads the UFF parameters
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
     integer :: hybridization
     real(8) :: ionizationPotential(9)
     logical :: isInstanced
  end type UFFParameters
  
  public :: &
       UFFParameters_load, &
       UFFParameters_show
       ! UFFParameters_getCovalentRadius

contains
    
  !>
  !! @brief Loads UFF parameters from library.
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the atom type information
  !! @param [in] typeSelected CHARACTER atom type to evaluate
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
    integer :: hybridization
    real(8) :: ionizationPotential(9)


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
         radius, &
         hybridization, &
         ionizationPotential


    !! Looking for library    
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%UFF_PARAMETERS_DATABASE), exist=existFile)
    
    if ( existFile ) then
       !! Open library
       open(unit=20, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%UFF_PARAMETERS_DATABASE), status="old", form="formatted" )
              
       !! Read information
       type = "NONE"
       stat = 0

       rewind(20)
       
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
          hybridization = 0
          ionizationPotential = 0

        
          if (stat == -1 ) then
             
             this%isInstanced = .false.
             return
             
          end if
          
          read(20,NML=atomtype, iostat=stat)

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
       this%hybridization = hybridization
       this%isInstanced = .true.

       do i = 1, size(ionizationPotential)
          this%ionizationPotential(i) = ionizationPotential(i)
       end do
       
       !! Debug information.
       ! call UFFParameters_show(this)
       
       close(20)
       
    else 

       call UFFParameters_exception( ERROR, "LOWDIN library not found!! please export lowdinvars.sh file.", "In UFFParameters at load function.")

    end if 

    !! Done
    
  end subroutine UFFParameters_load
  
  !<
  !! @brief This routine show the information for an atom type
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param this Class with informations of the atom type
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
    write (6,"(T10,A22,I5,A10)")	"Hybridization          = ",this%hybridization,""
    print *,""
    
  end subroutine UFFParameters_show

  !>
  !! @brief  Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
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
