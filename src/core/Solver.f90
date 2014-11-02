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
  implicit none

  type, public :: Solver
     character(20) :: methodName
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
    
    if( String_findSubstring( trim(CONTROL_instance%METHOD), "-" ) > 1 ) then

       lowdin_solver%methodName = CONTROL_instance%METHOD(1: scan(CONTROL_instance%METHOD, "-") - 1)

    else
       
       lowdin_solver%methodName = CONTROL_instance%METHOD
       
    end if
    
    lowdin_solver%withProperties = .false.
    
    select case ( trim(lowdin_solver%methodName) )
       
      case('RHF')
        call Solver_RHFRun( )
      case('UHF')
        call Solver_UHFRun( )
      case('RKS')
        call Solver_RKSRun( )
      case('UKS')
        call Solver_UKSRun( )
      case default
       
        call Solver_exception(ERROR, "The method: "//trim(lowdin_solver%methodName)//" is not implemented", &
             "At Solver module in run function")
        
     end select
     
   end subroutine Solver_run
  
  !> @brief run RHF-based calculation
  subroutine Solver_RHFRun( )
    implicit none

    !! Run HF program in RHF mode
    call system("lowdin-HF.x RHF")
    
    select case(CONTROL_instance%METHOD)
              
    case("RHF")

       call system("lowdin-HF.x RHF")
       
    case ("RHF-COSMO")
       
       call system("lowdin-cosmo.x")
       call system("lowdin-HF.x RHF")
       write(*,*) CONTROL_instance%METHOD
       
    case('RHF-MP2')

       call system("lowdin-integralsTransformation.x")

       call system("lowdin-MollerPlesset.x CONTROL_instance%MOLLER_PLESSET_CORRECTION")

    case('RHF-CI')

    case('RHF-PT')

       call system("lowdin-integralsTransformation.x")
       call system("lowdin-PT.x CONTROL_instance%PT_ORDER")
       
    case default

       call Solver_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented", &
            "At Solver module in RHFrun function")

    end select
    
!     call RHF_run()
!     if ( this%withProperties ) then
!        call CalculateProperties_dipole( CalculateProperties_instance )
!        call CalculateProperties_expectedPosition( CalculateProperties_instance )
!        if (Parameters%CALCULATE_INTERPARTICLE_DISTANCES .or. Parameters%CALCULATE_DENSITY_VOLUME ) &
!             call CalculateProperties_buildDensityCubesLimits( CalculateProperties_instance )
!        if (Parameters%CALCULATE_DENSITY_VOLUME) then
!           call CalculateProperties_buildDensityCubes( CalculateProperties_instance )
!           call CalculateProperties_volumes( CalculateProperties_instance )
!        end if
!        if (Parameters%CALCULATE_INTERPARTICLE_DISTANCES) call CalculateProperties_interparticleDistance( CalculateProperties_instance )
!        ! call CalculateProperties_interparticleOverlap( CalculateProperties_instance )
!        ! call CalculateProperties_expectedR2( CalculateProperties_instance )
!     end if
    
!     this%energy = MolecularSystem_getTotalEnergy()

  end subroutine Solver_RHFRun

  !> @brief run ROHF-based calculation
  subroutine Solver_ROHFRun( )
    implicit none
    
    select case(CONTROL_instance%METHOD)

    case("ROHF")

    case('ROHF-MP2')

    case('ROHF-CI')

    case('ROHF-PT')

    case default

       call Solver_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented", &
            "At Solver module in ROHFrun function")

    end select
    
!     call RHF_run()
!     if ( this%withProperties ) then
!        call CalculateProperties_dipole( CalculateProperties_instance )
!        call CalculateProperties_expectedPosition( CalculateProperties_instance )
!        if (Parameters%CALCULATE_INTERPARTICLE_DISTANCES .or. Parameters%CALCULATE_DENSITY_VOLUME ) &
!             call CalculateProperties_buildDensityCubesLimits( CalculateProperties_instance )
!        if (Parameters%CALCULATE_DENSITY_VOLUME) then
!           call CalculateProperties_buildDensityCubes( CalculateProperties_instance )
!           call CalculateProperties_volumes( CalculateProperties_instance )
!        end if
!        if (Parameters%CALCULATE_INTERPARTICLE_DISTANCES) call CalculateProperties_interparticleDistance( CalculateProperties_instance )
!        ! call CalculateProperties_interparticleOverlap( CalculateProperties_instance )
!        ! call CalculateProperties_expectedR2( CalculateProperties_instance )
!     end if
    
!     this%energy = MolecularSystem_getTotalEnergy()
    
  end subroutine Solver_ROHFRun

  !> @brief run UHF-based calculation
  subroutine Solver_UHFRun( )
    implicit none

    call system("lowdin-HF.x RHF")

    select case(CONTROL_instance%METHOD)
       
    case("UHF")
       
       !! Run HF program in RHF mode
       call system("lowdin-HF.x RHF")
       
    case('UHF-MP2')
       call system("lowdin-integralsTransformation.x")
       call system("lowdin-MollerPlesset.x CONTROL_instance%MOLLER_PLESSET_CORRECTION")
       
    case('UHF-PT')
       
       ! <<<<<<< head
       !        call system("lowdin-HF.x UHF") ???
       !        call system("lowdin-MOERI.x UHF") ???
       ! ===========
       
       call system("lowdin-integralsTransformation.x")
       call system("lowdin-PT.x CONTROL_instance%PT_ORDER")
       ! >>>>>>> master
       !rfm call system("lowdin-EPT.x UHF")
       
    case default
       
       call Solver_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented", &
            "At Solver module in UHFrun function")
       
    end select
    
    
!     type(Solver) :: this
    
!     call UHF_run()
!     if ( this%withProperties ) then
!        call CalculateProperties_dipole( CalculateProperties_instance )
!        call CalculateProperties_expectedPosition( CalculateProperties_instance )
!     end if
    
!     this%energy = MolecularSystem_getTotalEnergy()
    
  end subroutine Solver_UHFRun

  !> @brief run RKS-based calculation
  subroutine Solver_RKSRun( )
    implicit none
!     type(Solver) :: this
    
!     call RKS_run()
!     if ( this%withProperties ) then
!        call CalculateProperties_dipole( CalculateProperties_instance )
!        call CalculateProperties_expectedPosition( CalculateProperties_instance )
!        if (Parameters%POLARIZATION_ORDER>1) then
!           call CalculateProperties_polarizability( CalculateProperties_instance )
!        end if
!        if (Parameters%FUKUI_FUNCTIONS) then
!           call CalculateProperties_fukuiFunctions (CalculateProperties_instance)
!        end if
!     end if
    
!     this%energy = MolecularSystem_getTotalEnergy()
    
  end subroutine Solver_RKSRun
  
  !> @brief run UKS-based calculation
  subroutine Solver_UKSRun( )
    implicit none
!     type(Solver) :: this
    
!     call UKS_run()
!     if ( this%withProperties ) then
!        call CalculateProperties_expectedPosition( CalculateProperties_instance )
!        call CalculateProperties_dipole( CalculateProperties_instance )
!     end if
    
!     this%energy = MolecularSystem_getTotalEnergy()
    
  end subroutine Solver_UKSRun
  
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
