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
!! @brief Module to load constans of coupling for elemental particles
!! @author E. F. Posada, 2013
!! @version 1.0
module ConstantsOfCoupling_
  use CONTROL_
  use Exception_
  implicit none
  
  type , public :: ConstantsOfCoupling
     character(30) :: name
     character(30) :: symbol
     real(8) :: kappa
     real(8) :: eta
     real(8) :: lambda
     real(8) :: particlesFraction
     logical :: isInstanced
  end type ConstantsOfCoupling
  
  public :: &
       ConstantsOfCoupling_load, &
       ConstantsOfCoupling_show
  
contains
    
  !>
  !! @brief Loads constants of coupling from library.
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine ConstantsOfCoupling_load( this, symbolSelected )
    implicit none
    
    type(ConstantsOfCoupling), intent(inout) :: this
    character(*) :: symbolSelected

    logical :: existFile
    integer :: stat
    
    !! Namelist definition
    character(30) :: name
    character(30) :: symbol
    real(8) :: kappa
    real(8) :: eta
    real(8) :: lambda
    real(8) :: particlesFraction
    
    NAMELIST /specie/ &
         name, &
         symbol, &
         kappa, &
         eta, &
         lambda, &
         particlesFraction

    this%isInstanced = .true.
    
    !! Looking for library    
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//"/dataBases/constantsOfCoupling.lib", exist=existFile)
    
    if ( existFile ) then
       
       !! Open library
       open(unit=10, file=trim(CONTROL_instance%DATA_DIRECTORY)//"/dataBases/constantsOfCoupling.lib", status="old", form="formatted" )
       
       !! Read information
       symbol = "NONE"
       stat = 0
       
       do while(trim(symbol) /= trim(symbolSelected))
       
          !! Setting defaults
          name = "NONE"
          kappa = 0
          eta = 0
          particlesFraction = 1
          
          if (stat == -1 ) then
             
             call ConstantsOfCoupling_exception( ERROR, "Elemental particle: "//trim(symbolSelected)//" NOT found!!", "In ConstantsOfCoupling at load function.")
             this%isInstanced = .false.

          end if
          
          read(10,NML=specie, iostat=stat)
 
          if (stat > 0 ) then
             
             call ConstantsOfCoupling_exception( ERROR, "Failed reading ConstantsOfCouplings.lib file!! please check this file.", "In ConstantsOfCoupling at load function.")
             
          end if

       end do

       !! Set object variables
       this%name = name
       this%symbol = symbol
       this%kappa = kappa
       this%eta = eta
       this%lambda = lambda
       this%particlesFraction = particlesFraction
       
       !! Debug information.
       !! call ConstantsOfCoupling_show(this)
       
       close(10)
       
    else 

       call ConstantsOfCoupling_exception( ERROR, "LOWDIN library not found!! please export lowdinvars.sh file.", "In ConstantsOfCoupling at load function.")

    end if 

    !! Done
    
  end subroutine ConstantsOfCoupling_load
  
  !<
  !! @brief print out the object
  subroutine ConstantsOfCoupling_show( this )
    implicit none
    type(ConstantsOfCoupling) , intent(in) :: this
    
    
    print *,""
    print *,"====================="
    print *,"  Particle Properties"
    print *,"====================="
    print *,""
    write (6,"(T10,A10,A12)") "Name    = ",this%name
    write (6,"(T10,A10,A12)") "Symbol  = ",this%symbol
    write (6,"(T10,A10,F12.5)") "Mass    = ",this%kappa
    write (6,"(T10,A10,F12.5)") "Charge  = ",this%eta
    write (6,"(T10,A10,F12.5)") "Spin    = ",this%lambda
    write (6,"(T10,A10,F12.5)") "Spin    = ",this%particlesFraction
    print *,""
    
  end subroutine ConstantsOfCoupling_show
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine ConstantsOfCoupling_exception( typeMessage, description, debugDescription)
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
    
  end subroutine ConstantsOfCoupling_exception
  
end module ConstantsOfCoupling_
