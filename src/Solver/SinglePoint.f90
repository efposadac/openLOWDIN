module SinglePoint_
  use CONTROL_
  use MolecularSystem_
  use String_
  use Exception_
  use InputManager_
  use COSMO_ ,                  only : COSMO_main
  use SCF_ ,                    only : SCF_main
  use IntegralsTransformation_, only : IntegralsTransformation_main
  use MolecularMechanics_ ,     only : MolecularMechanics_main
  use CalcProp_ ,               only : CalcProp_main
  use MBPT_ ,                   only : MBPT_main
  use PT_ ,                     only : PT_main
  use CI_ ,                     only : CI_main
  use NOCI_ ,                   only : NOCI_main

  implicit none

  public SinglePoint_run

contains

!>
  !! @brief Run the properly programs depending of the requested tasks
  subroutine SinglePoint_run()
    implicit none
    character(100) :: auxString

    !Check cosmo
    if (CONTROL_instance%COSMO) call COSMO_main()

    !Do SCF
    select case (trim(CONTROL_instance%METHOD))

    case ('MM')
      call MolecularMechanics_main( CONTROL_instance%FORCE_FIELD )
    case ('RHF')
      call SCF_main( "RHF" )
    case ('UHF')
      call SCF_main( "UHF" )
    case ('RKS')
      call SCF_main( "RKS" )
    case ('UKS')
      call SCF_main( "UKS" )
    case default
      call SinglePoint_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented", &
                            "At Solver module in run function")
    end select

    !!calculate HF/KS HF/KS properties
    call CalcProp_main( "lowdin" )

    !Check for inconsistent methods
    if ((CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 .or. &
         CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 .or. &
         CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. &
         CONTROL_instance%PT_ORDER /= 0) .and. &
        (trim(CONTROL_instance%METHOD) .eq. "RKS" .or. trim(CONTROL_instance%METHOD) .eq. "UKS")) then
      print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      call SinglePoint_exception(WARNING, "You have selected a post-HF calculation that probably doesn't make sense with a KS reference."// &
                            " The calculation will proceed but be mindful of the results", "At Solver module in run function")
      print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    end if

    if ((CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. CONTROL_instance%PT_ORDER .ge. 3) .and. &
        trim(CONTROL_instance%METHOD) .ne. "UHF") then
      print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      call SinglePoint_exception(WARNING, "CI calculations have been tested only for UHF. You have selected "//trim(CONTROL_instance%METHOD)// &
                            " The calculation will proceed but be mindful of the results", "At Solver module in run function")
      print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    end if

    !Post SCF corrections
    if (CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0 .or. &
        CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0 .or. &
        CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" .or. &
        CONTROL_instance%PT_ORDER /= 0) call IntegralsTransformation_main()

    if (CONTROL_instance%MOLLER_PLESSET_CORRECTION /= 0) call MBPT_main()

    if (CONTROL_instance%EPSTEIN_NESBET_CORRECTION /= 0) call MBPT_main()

    if (CONTROL_instance%PT_ORDER /= 0) call PT_main()

    if (CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE") then
      write (auxString, "(I10)") Input_instance%numberOfSpeciesInCI
      call CI_main(trim(auxString))
       !!calculate CI density properties
      if (CONTROL_instance%CI_STATES_TO_PRINT .ge. 1) call CalcProp_main( "lowdin" )
    end if

    if (CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION) then
      call NOCI_main("POSTSCF")
       !!calculate CI density properties
      if (CONTROL_instance%CI_STATES_TO_PRINT .ge. 1 .and. &
          .not. (CONTROL_instance%COMPUTE_ROCI_FORMULA .or. CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS)) &
        call CalcProp_main( "lowdin" )
    end if

    ! if(optimization) then
    !    call system("lowdin-Optimizer.x")
    ! else
    ! end if

  end subroutine SinglePoint_run

  !>
  !! @brief Manejo de excepciones
  subroutine SinglePoint_exception(typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor(ex, typeMessage)
    call Exception_setDebugDescription(ex, debugDescription)
    call Exception_setDescription(ex, description)
    call Exception_show(ex)
    call Exception_destructor(ex)

  end subroutine SinglePoint_exception

end module SinglePoint_


