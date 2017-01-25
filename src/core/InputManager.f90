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
!! @brief Carga el archivo de entrada *.aux, All the characters variables for the system must be uppercase
!! @author E. F. Posada, 2013
!! @version 1.0
module InputManager_
  use Units_
  use String_
  use CONTROL_
  use Exception_
  use Particle_
  use MolecularSystem_
  implicit none

  
  type, public :: InputManager
     !!System
     character(100) :: fileName
     character(100) :: systemName
     character(255) :: systemDescription
     integer :: numberOfParticles
     integer :: numberOfExternalPots
     integer :: numberOfInterPots
     integer :: numberOfOutputs
     integer :: numberOfSpeciesInCI
     integer :: numberOfFragments
     !!Task
     character(50) :: method
     integer :: mollerPlessetCorrection
     character(20) :: configurationInteractionLevel
     integer :: propagatorTheoryCorrection
     logical :: optimizeGeometry
     logical :: TDHF
     logical :: cosmo
     logical :: ecp

  end type InputManager

  !> Singleton
  type(InputManager), public :: input_instance

  public :: &
       InputManager_loadSystem, &
       InputManager_loadControl, &
       InputManager_loadTask, &
       InputManager_loadGeometry, &
       InputManager_loadEcp

contains


  !>
  !! @brief load system definitions
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine InputManager_loadSystem()
    implicit none

    logical :: existFile
    integer :: status

    character(255):: InputSystem_description
    integer:: InputSystem_numberOfParticles
    integer:: InputSystem_numberOfExternalPots
    integer:: InputSystem_numberOfInterPots
    integer:: InputSystem_numberOfOutputs
    integer:: InputSystem_numberOfSpeciesInCI
    
    NAMELIST /InputSystem/ &
         InputSystem_numberOfParticles, &
         InputSystem_numberOfExternalPots, &
         InputSystem_numberOfInterPots, &
         InputSystem_numberOfOutputs, &
         InputSystem_numberOfSpeciesInCI, & 
         InputSystem_description
    
    !! Get input name
    call get_command_argument (1,value=Input_instance%fileName)

    !! Clear extension .aux --- first check the rigth extension name.
    if(Input_instance%fileName(len_trim(Input_instance%fileName)-3:len_trim(Input_instance%fileName)) /= ".aux") then       
       call InputManager_exception( ERROR, "Invalid input file!. It must be a *.aux file", "InputManager loadSystem function")
    end if
    
    Input_instance%fileName = Input_instance%fileName(1:len_trim(Input_instance%fileName)-3)
    Input_instance%systemName = trim(Input_instance%fileName(1:len_trim(Input_instance%fileName)-1))

    !! Setting name in global parameters
    
    
    !! Exist that file?
    inquire(FILE = trim(Input_instance%fileName)//"aux", EXIST = existFile )
       
    if ( existFile ) then
       
       !! Open file
       open (unit=4, file=trim(Input_instance%fileName)//"aux")
       
       !! Setting defaults
       InputSystem_numberOfExternalPots=0
       InputSystem_numberOfInterPots=0
       InputSystem_numberOfOutputs=0
       InputSystem_numberOfSpeciesInCI=0
       InputSystem_numberOfParticles=0
       InputSystem_description=trim(Input_instance%systemName)

       !! Load inputSystem namelist
       rewind(4)
       read(4,NML=InputSystem, iostat=status)
       
       if(status /= 0) then
          
          call InputManager_exception( ERROR, "Error reading "// trim(Input_instance%fileName)//"aux" //" file!", "Class object InputManager in the constructor function")
          
       end if
              
       Input_instance%systemDescription = trim(InputSystem_description)       
       Input_instance%numberOfParticles = InputSystem_numberOfParticles
       Input_instance%numberOfExternalPots = InputSystem_numberOfExternalPots
       Input_instance%numberOfInterPots = InputSystem_numberOfInterPots
       Input_instance%numberOfOutputs = InputSystem_numberOfOutputs
       Input_instance%numberOfSpeciesInCI = InputSystem_numberOfSpeciesInCI

       !!done
       
    else
       
       call InputManager_exception( ERROR, "The input file "// trim(Input_instance%fileName)//"aux" //" not found!", "Class object InputManager in the constructor function")
       
    end if
    
  end subroutine InputManager_loadSystem

  !>
  !! @brief Load Control block
  !! @author E. F. Posada
  !! @version 1.0
  subroutine InputManager_loadControl()
    implicit none
    
    !! Setting all defaults
    call CONTROL_start()

    !! Loads Control Block
    call CONTROL_load()

    CONTROL_instance%INPUT_FILE = trim(input_instance%fileName)
    
  end subroutine InputManager_loadControl

  !>
  !! @brief Load tasks
  !! @author E. F. Posada
  !! @version 1.0
  subroutine InputManager_loadTask()
    implicit none
    
    integer :: stat

    !! Namelist definition
    character(50):: InputTasks_method
    character(20):: InputTasks_configurationInteractionLevel
    integer:: InputTasks_mollerPlessetCorrection
    integer:: InputTasks_propagatorTheoryCorrection
    logical:: InputTasks_optimizeGeometry
    logical:: InputTasks_TDHF
    logical:: InputTasks_cosmo
    logical:: InputTasks_ecp

    
    NAMELIST /InputTasks/ &
         InputTasks_method, &
         InputTasks_configurationInteractionLevel, &
         InputTasks_mollerPlessetCorrection, &
         InputTasks_propagatorTheoryCorrection, &
         InputTasks_optimizeGeometry, &
         InputTasks_TDHF, &
         InputTasks_cosmo, &
         InputTasks_ecp

    
    !! Setting defaults    
    InputTasks_method = "NONE"
    InputTasks_mollerPlessetCorrection = 0
    InputTasks_configurationInteractionLevel = "NONE"
    InputTasks_propagatorTheoryCorrection = 0
    InputTasks_optimizeGeometry = .false.
    InputTasks_TDHF = .false.
    InputTasks_cosmo= .false.
    InputTasks_ecp = .false.
    
    !! reload input file
    rewind(4)
    
    !! Read InputTask namelist from input file
    read(4,NML=InputTasks, iostat=stat)
    
    if( stat > 0 ) then       
       call InputManager_exception( ERROR, "check the TASKS block in your input file", "InputManager loadTask function" )       
    end if
    
    !! all uppercase! Mandatory for ALL character variables
    Input_instance%method = trim(String_getUppercase(trim(InputTasks_method)))
    Input_instance%mollerPlessetCorrection = InputTasks_mollerPlessetCorrection
    Input_instance%configurationInteractionLevel = trim(String_getUppercase(trim(InputTasks_configurationInteractionLevel)))
    Input_instance%propagatorTheoryCorrection = InputTasks_propagatorTheoryCorrection
    Input_instance%optimizeGeometry = InputTasks_optimizeGeometry
    Input_instance%TDHF = InputTasks_TDHF
    Input_instance%cosmo = InputTasks_cosmo
    Input_instance%ecp = InputTasks_ecp
    
    !! If the method is for open shell systems
    if ( trim(Input_instance%method) == "UHF" .or. trim(Input_instance%method) == "ROHF" .or. & 
         trim(Input_instance%method) == "UKS"  .or. trim(Input_instance%method) == "ROKS") then
       
       CONTROL_instance%IS_OPEN_SHELL = .true.
       
    end if
    
    !! If there is no method in the input file
    if( Input_instance%method == "NONE" .or. Input_instance%method == "" ) then       
       call InputManager_exception( ERROR, "check the TASKS block in  input file, what method you want to use? I dont have super cow powers!", "InputManager loadTask function" )       
    end if
    
    !! Setting some parameters-object variables
    CONTROL_instance%METHOD = trim(input_instance%method)
    CONTROL_instance%MOLLER_PLESSET_CORRECTION = input_instance%mollerPlessetCorrection
    CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL = input_instance%configurationInteractionLevel
    CONTROL_instance%PT_ORDER = input_instance%propagatorTheoryCorrection

    if ( input_instance%mollerPlessetCorrection /= 0 ) then
       CONTROL_instance%METHOD=trim(CONTROL_instance%METHOD)//"-MP2"
    end if
        
    if ( input_instance%configurationInteractionLevel /= "NONE" ) then
       CONTROL_instance%METHOD=trim(CONTROL_instance%METHOD)//"-CI"
    end if

    if ( input_instance%propagatorTheoryCorrection /= 0 ) then
       CONTROL_instance%METHOD=trim(CONTROL_instance%METHOD)//"-PT"
    end if
    
    if ( input_instance%optimizeGeometry ) then 
       CONTROL_instance%OPTIMIZE = .true.
    end if
    
    if ( input_instance%TDHF ) then 
       CONTROL_instance%TDHF = .true.
    end if    
    
    if ( input_instance%ecp ) then 
       CONTROL_instance%ecp = .true.
    end if    
    
    if (input_instance%cosmo) then 	
       CONTROL_instance%cosmo = .true.	
       CONTROL_instance%METHOD=trim(CONTROL_instance%METHOD)//"-COSMO"
    end if
    
    if( Input_instance%numberOfExternalPots > 0) then    
       CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL=.true.
    end if
        
    if(Input_instance%numberOfInterPots > 0) then
       CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL=.true.
    end if
           
    if(Input_instance%numberOfOutputs > 0) then
       CONTROL_instance%IS_THERE_OUTPUT=.true.
    end if
        
  end subroutine InputManager_loadTask  

  !>
  !! @brief Load all particles
  !! @author E. F. Posada
  !! @version 1.0
  subroutine InputManager_loadGeometry()
    implicit none
    
    character(15), allocatable :: quantumSpeciesName(:)
    integer, allocatable :: numberOfParticlesForSpecies(:)
    integer, allocatable :: particlesID(:)

    character(15) :: atomName
    integer :: stat
    integer :: numberOfQuantumSpecies
    integer :: numberOfPointCharges
    integer :: speciesID
    integer :: i, j, counter
    logical :: isNewSpecies
    logical :: isElectron
    !! Namelist definition
    character(15):: InputParticle_name
    character(30):: InputParticle_basisSetName
    real(8):: InputParticle_origin(3)
    real(8) :: InputParticle_charge
    character(3):: InputParticle_fixedCoordinates
    integer:: InputParticle_addParticles
    real(8):: InputParticle_multiplicity
    character(3):: InputParticle_ecp
    
    NAMELIST /InputParticle/ &
         InputParticle_name, &
         InputParticle_basisSetName, &
         InputParticle_charge, &
         InputParticle_origin, &
         InputParticle_fixedCoordinates, &
         InputParticle_multiplicity, &
         InputParticle_addParticles, &
         InputParticle_ecp
    
    !! Allocate memory for buffer.
    allocate(quantumSpeciesName(Input_instance%numberOfParticles))
    allocate(numberOfParticlesForSpecies(Input_instance%numberOfParticles))

    !! Initializes some variables.
    quantumSpeciesName = ""
    numberOfParticlesForSpecies = 0
    numberOfQuantumSpecies = 0
    numberOfPointCharges = 0
    counter = 0

    !! Reload input file
    rewind(4)
    
    !!*****************************************************************************
    !! Explore GEOMETRY block 
    !! (determines the number of quantum species, the number of point charges 
    !! and the number of particles by specie)
    do i=1, Input_instance%numberOfParticles
       
       !! Setting defaults
       InputParticle_name = "NONE"
       InputParticle_basisSetName = "NONE"
       
       !! Reads namelist from input file
       read(4,NML = InputParticle, iostat = stat)
       
       if( stat > 0 ) then          
          call InputManager_exception( ERROR, "check the GEOMETRY block in your input file", "InputManager loadParticles function")          
       end if
       
       !! All uppercase
       InputParticle_name = trim(String_getUppercase(trim(InputParticle_name)))
       InputParticle_basisSetName = trim(String_getUppercase(trim(InputParticle_basisSetName)))
       
       !! Electrons case
       isElectron = .false.
       if( String_findSubstring( trim(InputParticle_name), "E-") == 1 ) then
          isElectron = .true.
          if(CONTROL_instance%IS_OPEN_SHELL) then
             InputParticle_name = "E-ALPHA"
          else
             InputParticle_name = "E-"
          end if
       end if
       
       !! identifies if this species is a new species
       isNewSpecies = .true.
       do j = 1, Input_instance%numberOfParticles
          !! Search for this particle in the buffer
          if(trim(InputParticle_name) == quantumSpeciesName(j)) then
             !! Counters for old species
             isNewSpecies = .false.
             !! Only are species those who have basis-set
             if((trim(InputParticle_basisSetName) /= "MM") .and. &
                  (trim(InputParticle_basisSetName) /= "DIRAC") .and. (trim(InputParticle_basisSetName) /= "") &
                  .and. (trim(InputParticle_basisSetName) /= "ECP")) then
                
                !! For open-shell case electrons are alpha and beta
                if(isElectron .and. CONTROL_instance%IS_OPEN_SHELL) then
                   !! Alpha
                   numberOfParticlesForSpecies(j) = numberOfParticlesForSpecies(j) + 1
                   !! Beta
                   numberOfParticlesForSpecies(j+1) = numberOfParticlesForSpecies(j+1) + 1
                else
                   numberOfParticlesForSpecies(j) = numberOfParticlesForSpecies(j) + 1
                end if

             else
                
                numberOfPointCharges = numberOfPointCharges + 1
                
             end if
             
             exit             
          end if          
       end do
       
       !! Counters for new species
       if(isNewSpecies) then          
          !! Look for Quantum species
          if((trim(InputParticle_basisSetName) /= "MM") .and. &
               (trim(InputParticle_basisSetName) /= "DIRAC") .and. (trim(InputParticle_basisSetName) /= "") &
                  .and. (trim(InputParticle_basisSetName) /= "ECP")) then 
             
             !! For open-shell case electrons are alpha and beta
             if(isElectron .and. CONTROL_instance%IS_OPEN_SHELL) then
                !! Aplha
                counter = counter + 1
                !! Store new species in the buffer (to avoid count again)
                quantumSpeciesName(counter) = "E-ALPHA"
                !! counter to calculate number Of Particles For this species
                numberOfParticlesForSpecies(counter) = numberOfParticlesForSpecies(counter) + 1
                numberOfQuantumSpecies = numberOfQuantumSpecies + 1
                !! Beta
                counter = counter + 1
                !! Store new species in the buffer (to avoid count again)
                quantumSpeciesName(counter) = "E-BETA"
                !! counter to calculate number Of Particles For this species
                numberOfParticlesForSpecies(counter) = numberOfParticlesForSpecies(counter) + 1
                numberOfQuantumSpecies = numberOfQuantumSpecies + 1
             else
                !! Other quantum species
                counter = counter + 1
                !! Store new species in the buffer (to avoid count again)
                quantumSpeciesName(counter) = trim(InputParticle_name)
                !! counter to calculate number Of Particles For this species
                numberOfParticlesForSpecies(counter) = numberOfParticlesForSpecies(counter) + 1
                numberOfQuantumSpecies = numberOfQuantumSpecies + 1
             end if
                
          else !! Point Charges
             
             numberOfPointCharges = numberOfPointCharges + 1
             
          end if
          
       end if
    end do
    
    !! Reshape arrays
    quantumSpeciesName = reshape(quantumSpeciesName, (/numberOfQuantumSpecies/))
    numberOfParticlesForSpecies = reshape(numberOfParticlesForSpecies, (/numberOfQuantumSpecies/))
    

    !!*****************************************************************************
    !! LOAD GEOMETRY block 
    !!
    !! Allocate particlesID
    allocate(particlesID(numberOfQuantumSpecies))
    particlesID = 0
    counter = 0
    
    !! Initializes the molecular system object
    call MolecularSystem_initialize(numberOfQuantumSpecies, numberOfPointCharges, numberOfParticlesForSpecies, quantumSpeciesName, &
         trim(input_instance%systemName), trim(input_instance%systemDescription))
    
    !! Reload input file
    rewind(4)
    
    do i=1, Input_instance%numberOfParticles
       
       !! Setting defaults
       InputParticle_name = "NONE"
       InputParticle_basisSetName = "NONE"
       InputParticle_charge=0.0_8
       InputParticle_origin=0.0_8
       InputParticle_fixedCoordinates = "NONE"
       InputParticle_multiplicity = 1.0_8
       InputParticle_addParticles = 0
       InputParticle_ecp = "NONE"
       
       !! Reads namelist from input file
       read(4,NML = InputParticle, iostat = stat)
       if( stat > 0 ) then          
          call InputManager_exception( ERROR, "check the GEOMETRY block in your input file", "InputManager loadParticles function")          
       end if
       
       !! Fix units ( converts ANGS to BOHR )
       if ( trim(CONTROL_instance%UNITS) == "ANGS") then
          InputParticle_origin= InputParticle_origin / AMSTRONG
       end if
       
       !! All uppercase
       InputParticle_name = trim(String_getUppercase(trim(InputParticle_name)))
       InputParticle_basisSetName = trim(String_getUppercase(trim(InputParticle_basisSetName)))
       InputParticle_fixedCoordinates = trim(String_getUppercase(trim(InputParticle_fixedCoordinates)))
       
       !!***************************************************************
       !! Start loading particles 
       
       !! Electrons case
       isElectron = .false.
       if( String_findSubstring( trim(InputParticle_name), "E-") == 1 ) then          
          isElectron = .true.
          !! Extract atom's symbol
          atomName = InputParticle_name (scan( InputParticle_name, "[" ) : scan( InputParticle_name, "]" )) 
          atomName = trim(adjustl(atomName))
          
          if(CONTROL_instance%IS_OPEN_SHELL) then
             InputParticle_name = "E-ALPHA"
          else
             InputParticle_name = "E-"
          end if
          
       end if
       
       !! Load quantum species
       if((trim(InputParticle_basisSetName) /= "MM") .and. &
            (trim(InputParticle_basisSetName) /= "DIRAC") .and. (trim(InputParticle_basisSetName) /= "") &
                  .and. (trim(InputParticle_basisSetName) /= "ECP")) then

          !! Locate specie
          do j = 1, numberOfQuantumSpecies          
             
             if(trim(InputParticle_name) == trim(quantumSpeciesName(j))) then             
                speciesID = j
                exit
             end if
             
          end do
          
          !! Load open-shell case electrons
          if(isElectron .and. CONTROL_instance%IS_OPEN_SHELL) then

             !!ALPHA SET
             particlesID(speciesID) = particlesID(speciesID) + 1
             InputParticle_name = "E-ALPHA-"//trim(atomName)
             
             !! Loads Particle
             call Particle_load( MolecularSystem_instance%species(speciesID)%particles(particlesID(speciesID)), &
                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
                  multiplicity=inputParticle_multiplicity, spin="ALPHA", id = particlesID(speciesID), ecp = inputParticle_ecp )
             
             !!BETA SET
             speciesID = speciesID + 1
             
             particlesID(speciesID) = particlesID(speciesID) + 1
             InputParticle_name = "E-BETA-"//trim(atomName)             
             
             !! Loads Particle
             call Particle_load( MolecularSystem_instance%species(speciesID)%particles(particlesID(speciesID)), &
                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
                  multiplicity=inputParticle_multiplicity, spin="BETA", id = particlesID(speciesID), ecp = inputParticle_ecp )
             
          else 

             !! Set name for particle, in electrons case.
             if(isElectron) then
                InputParticle_name = "E-"//trim(atomName)
             end if
             
             !! loads electrons for Closed-shell case and other quantum particles.
             particlesID(speciesID) = particlesID(speciesID) + 1
             
             !! Loads Particle
             call Particle_load( MolecularSystem_instance%species(speciesID)%particles(particlesID(speciesID)),&
                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
                  multiplicity=inputParticle_multiplicity, id = particlesID(speciesID), ecp = inputParticle_ecp)


             
          end if

       else
          
          !! Loads Point charges   
          counter = counter + 1
          !! Loads Molecular Mechanics Particle 
          if(trim(InputParticle_basisSetName) == "MM") then
             call Particle_load( MolecularSystem_instance%pointCharges(counter),&
                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
                  multiplicity=inputParticle_multiplicity, id = counter, charge = inputParticle_charge)
          !! Loads ECP particle (cores)
          else if(trim(InputParticle_basisSetName) == "ECP") then
!                  call Particle_load( MolecularSystem_instance%pointCharges(counter),&
!                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
!                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
!                  multiplicity=inputParticle_multiplicity, id = counter, charge = inputParticle_charge)
                call InputManager_exception( ERROR, "Aquí va la información de ECP")       
          else
          !! Loads Particle
             call Particle_load( MolecularSystem_instance%pointCharges(counter),&
                  name = trim(InputParticle_name), baseName = trim(InputParticle_basisSetName), &
                  origin = inputParticle_origin, fix=trim(inputParticle_fixedCoordinates), addParticles=inputParticle_addParticles, &
                  multiplicity=inputParticle_multiplicity, id = counter)
          end if
       end if
    end do

  end subroutine InputManager_loadGeometry

  !>
  !! @brief Retorna la descripcion del sistema
  function InputManager_getSystemDescription() result(output)
    implicit none
    character(255) :: output
    
    output = trim(Input_instance%systemDescription)
    
  end function InputManager_getSystemDescription
  
  !>
  !! @brief Retorna el nombre del sistema
  function InputManager_getSystemName() result(output)
    implicit none
    character(100) :: output
    
    output = trim(Input_instance%systemName)
    
  end function InputManager_getSystemName
  
  !>
  !! @brief  Handle class exceptions
  subroutine InputManager_exception( typeMessage, description, debugDescription)
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

  end subroutine InputManager_exception

end module InputManager_
