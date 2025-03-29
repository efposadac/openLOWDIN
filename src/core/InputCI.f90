!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!  Prof. G. MERINO's Lab. Universidad de Guanajuato
!!    http://quimera.ugto.mx/qtc/gmerino.html
!!
!!  Authors:
!!    E. F. Posada (efposadac@unal.edu.co)
!!    R. Flores (roberto.floresmoreno.qt@gmail.com)
!!
!!  Contributors:
!!
!!    Todos los derechos reservados, 2011
!!
!!******************************************************************************

module InputCI_
  use Exception_
  use CONTROL_
  use MolecularSystem_
  implicit none

  !>
  !! @brief Description
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 08-04-11
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 08-04-11 </tt>:  felix ( email@server )
  !!        -# description.
  !!   - <tt> 10-31-2014 </tt>:  Jorge Charry ( jacharry@unal.edu.co )
  !!        -# Adapts this module to Lowdin2
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<
  type, public :: InputCI
    character(50) :: species
    integer :: coreOrbitals
    integer :: activeOrbitals
    integer :: excitationType
		logical :: isInstanced
  end type

  character(50) :: InputCI_species
  integer :: InputCI_core
  integer :: InputCI_active
  integer :: InputCI_excitation

  NAMELIST /InputCINamelist/ &
            InputCI_species, &
            InputCI_core, &
            InputCI_active, &
            InputCI_excitation

  public :: &
    InputCI_constructor, &
    InputCI_destructor, &
    InputCI_show, &
    InputCI_load
    
private  

  !<Singleton
  type(InputCI), allocatable, public :: InputCI_Instance(:)

contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine InputCI_constructor()
    implicit none
    integer :: ssize
    integer :: i
    character(50) :: fileName
    integer :: fileUnit

    ssize = MolecularSystem_getNumberOfQuantumSpecies()
    fileUnit = 4
    fileName = trim(CONTROL_instance%INPUT_FILE)//"aux"

    if(.not.allocated(InputCI_Instance)) then
       allocate(InputCI_Instance(ssize))
       InputCI_Instance%species = ""
       InputCI_Instance%coreOrbitals = 0
       InputCI_Instance%activeOrbitals = 0
       InputCI_Instance%excitationType = 0
    end if
    open (unit=fileUnit, file= fileName, status="old")

  end subroutine InputCI_constructor

  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine InputCI_destructor()
    implicit none
    integer :: i 

    if(allocated(InputCI_Instance)) then
      deallocate(InputCI_Instance)
    end if

  end subroutine InputCI_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine InputCI_show(this)
    implicit none
    type(InputCI) :: this
  end subroutine InputCI_show

  !>
  !! @brief Carga la informacion de potenciales externos desde el input
  !!
  !! @param this
  !<
  subroutine  InputCI_load( numberOfSpeciesInCIinput )
    implicit none
    integer :: i,newID
    integer :: stat
    integer :: numberOfSpeciesInCIinput
    character(50) :: wfnFile,labels(2)
    integer :: wfnUnit
    real(8) :: removedOrbitals

    if ( allocated(InputCI_Instance) ) then
       rewind(4)

       do i=1, numberOfSpeciesInCIinput
          InputCI_species=""
          InputCI_core=0
          InputCI_active=0
          InputCI_excitation=0

          read(4,NML=InputCINamelist, iostat=stat)

          if(stat > 0 ) then

             call InputCI_exception( ERROR, "Class object InputCI in the load function", &
                  "check the INPUT_CI block in your input file")
          else
             if(trim(InputCI_species) .ne. "" ) then
             
                newID = MolecularSystem_getSpeciesIDFromSymbol (InputCI_species)

                if ( newID .ne. 0 ) then
                   InputCI_Instance(newID)%species = trim(InputCI_species)
                   InputCI_Instance(newID)%coreOrbitals = InputCI_core
                   InputCI_Instance(newID)%activeOrbitals = InputCI_active
                   InputCI_Instance(newID)%excitationType = InputCI_excitation
                   
                else
                   call InputCI_exception( ERROR, "Class object InputCI in the load function", &
                     "check the name of the species in the INPUT_CI block of your input file")
                end if
             else
                InputCI_Instance(i)%species = MolecularSystem_getNameOfSpecies(i)          
                InputCI_excitation=0
                if( CONTROL_instance%MP_FROZEN_CORE_BOUNDARY .ne. 0 &
                     .and. (trim(InputCI_Instance(i)%species) .eq. "E-" .or. trim(InputCI_Instance(i)%species) .eq. "E-ALPHA" .or. trim(InputCI_Instance(i)%species) .eq. "E-BETA")) &
                     InputCI_Instance(i)%coreOrbitals=CONTROL_instance%MP_FROZEN_CORE_BOUNDARY

                !! Check for orbitals removed in the SCF
                ! if ( .not. CONTROL_instance%SUBSYSTEM_EMBEDDING) then
                wfnFile = "lowdin.wfn"
                ! else
                !    wfnFile = "lowdin-subsystemA.wfn"
                ! end if
                wfnUnit = 20
                open(unit = wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

                labels(2) = InputCI_Instance(i)%species
                labels(1) = "REMOVED-ORBITALS"
                call Vector_getFromFile(unit=wfnUnit, binary=.true., value=removedOrbitals, arguments=labels)
                close(wfnUnit)

                if(removedOrbitals .ne. 0.0) InputCI_Instance(i)%activeOrbitals = MolecularSystem_getTotalNumberOfContractions(i)- int(removedOrbitals) 

             end if
          end if
       end do
       close(4)

       !      call OutputManager_constructor( OutputManager_instance, &
       !                                                         InputCI_Instance%type, &
       !                                                         InputCI_Instance%specie, & 
       !                                                         InputCI_Instance%orbital, &
       !                                                         InputCI_Instance%dimensions, &
       !                                                         InputCI_Instance%cubeSize, &
       !                                                         InputCI_Instance%point1, & 
       !                                                         InputCI_Instance%point2, &
       !                                                         InputCI_Instance%point3  )
       !                        
    else

       call InputCI_exception( ERROR, "Class object InputCI in the load function", &
            "The Input_Parsing module wasn't instanced")
    end if

  end subroutine InputCI_load


  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  !<
  function InputCI_isInstanced( this ) result( output )
    implicit  none
    type(InputCI), intent(in) :: this
    logical :: output
    
    output = this%isInstanced
  
  end function InputCI_isInstanced

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine InputCI_exception( typeMessage, description, debugDescription)
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
  
  end subroutine InputCI_exception

end module InputCI_
