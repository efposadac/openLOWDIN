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
!! @brief Module for elemental particles definitions
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
module ElementalParticle_
  use CONTROL_
  use Exception_
  implicit none
  
  type , public :: ElementalParticle
     character(30) :: name
     character(30) :: symbol
     character(30) :: category
     real(8) :: mass
     real(8) :: charge
     real(8) :: spin
  end type ElementalParticle
  
  public :: &
       ElementalParticle_load, &
       ElementalParticle_show
  
contains
    
  !>
  !! @brief Loads an elemental particles from library.
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine ElementalParticle_load( this, symbolSelected )
    implicit none
    
    type(ElementalParticle), intent(inout) :: this
    character(*) :: symbolSelected

    logical :: existFile
    integer :: stat
    integer :: i
    
    !! Namelist definition
    character(30) :: name
    character(30) :: symbol
    character(30) :: category
    real(8) :: charge
    real(8) :: mass
    real(8) :: spin

    NAMELIST /particle/ &
         name, &
         symbol, &
         category, &
         charge, &
         mass, &
         spin
    
    !! Looking for library    
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%ELEMENTAL_PARTICLES_DATABASE), exist=existFile)
    
    if ( existFile ) then
       
       !! Open library
       open(unit=10, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%ELEMENTAL_PARTICLES_DATABASE), status="old", form="formatted" )
              
       !! Read information
       symbol = "NONE"
       stat = 0
       
       do while(trim(symbol) /= trim(symbolSelected))
       
          !! Setting defaults
          name = "NONE"
          category = "NONE"
          mass = -1
          charge = 0
          spin = 0
          
          if (stat == -1 ) then
             
             call ElementalParticle_exception( ERROR, "Elemental particle: "//trim(symbolSelected)//" NOT found!!", "In ElementalParticle at load function.")

          end if
          
          read(10,NML=particle, iostat=stat)
 
          if (stat > 0 ) then
             
             call ElementalParticle_exception( ERROR, "Failed reading ElementalParticles.lib file!! please check this file.", "In ElementalParticle at load function.")
             
          end if

       end do
       
       !! Set object variables
       this%name = name
       this%symbol = symbol
       this%category = category 
       this%mass = mass
       this%charge = charge
       this%spin = spin
       
       !! Debug information.
       !! call ElementalParticle_show(this)
       
       close(10)
       
    else 

       call ElementalParticle_exception( ERROR, "LOWDIN library not found!! please export lowdinvars.sh file.", "In ElementalParticle at load function.")

    end if 

    !! Done
    
  end subroutine ElementalParticle_load
  
  !<
  !! @brief print out the object
  subroutine ElementalParticle_show( this )
    implicit none
    type(ElementalParticle) , intent(in) :: this
    
    
    print *,""
    print *,"====================="
    print *,"  Particle Properties"
    print *,"====================="
    print *,""
    write (6,"(T10,A10,A12)") "Name    = ",this%name
    write (6,"(T10,A10,A12)") "Symbol  = ",this%symbol
    write (6,"(T10,A10,F12.5)") "Mass    = ",this%mass
    write (6,"(T10,A10,F12.5)") "Charge  = ",this%charge
    write (6,"(T10,A10,F12.5)") "Spin    = ",this%spin
    print *,""
    
  end subroutine ElementalParticle_show
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine ElementalParticle_exception( typeMessage, description, debugDescription)
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
    
  end subroutine ElementalParticle_exception
  
end module ElementalParticle_
