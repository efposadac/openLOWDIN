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

  type, public :: Solver
     logical :: withProperties
  end type Solver
  
  public :: &
       Solver_run
  
  !> Singleton lock
  type(Solver), public :: lowdin_solver
    
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

    !Post SCF corrections
    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 .or. &
         CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 .or. &
         CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. &
         CONTROL_instance%PT_ORDER /= 0) then
       call system("lowdin-integralsTransformation.x")
    end if

    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 ) then
       call system("lowdin-MBPT.x CONTROL_instance%MOLLER_PLESSET_CORRECTION")
    end if

    if ( CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 ) then
       call system("lowdin-MBPT.x CONTROL_instance%EPSTEIN_NESBET_CORRECTION")
    end if
        
    if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" ) then
       write(auxString,"(I10)") Input_instance%numberOfSpeciesInCI
       call system("lowdin-CI.x" //trim(auxString))
    end if

    if ( CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION ) then
       call system("lowdin-CI.x NOCI")
    end if

    if ( CONTROL_instance%PT_ORDER /= 0 ) then
       call system("lowdin-PT.x CONTROL_instance%PT_ORDER")
    end if
   
    ! lowdin_solver%withProperties = .false.
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
