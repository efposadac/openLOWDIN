!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Selecciona el metodo apropiado para el calculo solicitado
!! @author Sergio Gonzalez
!! <b> Creation data : </b> 02-15-11
!! <b> History change: </b>
!!   - <tt> 02-15-11 </tt>:  fernando ( sagonzalez@unal.edu.co )
!!        -# Creacion del Modulo
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
module Solver_
  use CONTROL_
  use MolecularSystem_
  use String_
  use Exception_
  use InputManager_
  implicit none

  public :: &
       Solver_run
  
contains

  !>
  !! @brief Run the properly programs depending of the requested tasks
  subroutine Solver_run( )
    implicit none
    character(100) :: auxString

    !Check cosmo
    if (CONTROL_instance%COSMO) call system("lowdin-cosmo.x")

    !Do SCF
    select case ( trim(CONTROL_instance%METHOD) )

    case('MM')
       call system("lowdin-MolecularMechanics.x CONTROL_instance%FORCE_FIELD")
    case('RHF')
       call system("lowdin-SCF.x RHF")
    case('UHF')
       call system("lowdin-SCF.x UHF")
    case('RKS')
       call system("lowdin-SCF.x RKS")
    case('UKS')
       call system("lowdin-SCF.x UKS")
    case default
       call Solver_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented", &
            "At Solver module in run function")
    end select

    !!calculate HF/KS HF/KS properties
    call system ("lowdin-CalcProp.x")

    !Check for inconsistent methods
    if ( (CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 .or. &
         CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 .or. &
         CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. &
         CONTROL_instance%PT_ORDER /= 0) .and. &
         (trim(CONTROL_instance%METHOD) .eq. "RKS" .or. trim(CONTROL_instance%METHOD) .eq. "UKS")) then
       print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       call Solver_exception(WARNING, "You have selected a post-HF calculation that probably doesn't make sense with a KS reference."// &
            " The calculation will proceed but be mindful of the results", "At Solver module in run function")
       print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    end if

    if ( (CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. CONTROL_instance%PT_ORDER .ge. 3)  .and. &
         trim(CONTROL_instance%METHOD) .ne. "UHF" ) then
       print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       call Solver_exception(WARNING, "CI calculations have been tested only for UHF. You have selected "//trim(CONTROL_instance%METHOD)//&
            " The calculation will proceed but be mindful of the results", "At Solver module in run function")
       print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    end if
    
    !Post SCF corrections
    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 .or. &
         CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 .or. &
         CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. &
         CONTROL_instance%PT_ORDER /= 0) call system("lowdin-integralsTransformation.x")

    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 ) call system("lowdin-MBPT.x CONTROL_instance%MOLLER_PLESSET_CORRECTION")

    if ( CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 ) call system("lowdin-MBPT.x CONTROL_instance%EPSTEIN_NESBET_CORRECTION")

    if ( CONTROL_instance%PT_ORDER /= 0 ) call system("lowdin-PT.x CONTROL_instance%PT_ORDER")
    
    if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" ) then
       write(auxString,"(I10)") Input_instance%numberOfSpeciesInCI
       call system("lowdin-CI.x" //trim(auxString))
       !!calculate CI density properties
       if (CONTROL_instance%CI_STATES_TO_PRINT .ge. 1) call system ("lowdin-CalcProp.x")
    end if

    if ( CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION ) then
       call system("lowdin-NOCI.x POSTSCF")
       !!calculate CI density properties
       if ( CONTROL_instance%CI_STATES_TO_PRINT .ge. 1 .and. &
            .not. (CONTROL_instance%COMPUTE_ROCI_FORMULA .or. CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS)) &
            call system ("lowdin-CalcProp.x")
    end if

    ! if(optimization) then
    !    call system("lowdin-Optimizer.x")
    ! else
    ! end if
     
   end subroutine Solver_run
  
  !>
  !! @brief Manejo de excepciones
  subroutine Solver_exception(typeMessage, description, debugDescription)
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
    
  end subroutine Solver_exception
  
end module Solver_
