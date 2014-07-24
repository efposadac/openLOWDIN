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
module AtomicElement_
  use CONTROL_
  use Exception_
  implicit none
  
  type , public :: AtomicElement
     character(30) :: name
     character(30) :: symbol
     real(8) :: atomicNumber
     real(8) :: massicNumber
     real(8) :: abundance
     real(8) :: meltingPoint
     real(8) :: boilingPoint
     real(8) :: density
     real(8) :: electronAffinity
     real(8) :: electronegativity
     real(8) :: ionizationEnergy1
     real(8) :: covalentRadius
     real(8) :: atomicRadio
     real(8) :: vanderWaalsRadio
     real(8) :: atomicWeight
     real(8) :: nuclearSpin     
     !!< Variables con propositos de conveniencia
     logical :: isMassNumberPresent
     logical :: isInstanced
  end type AtomicElement
  
  public :: &
       AtomicElement_load, &
       AtomicElement_show
  
contains
    
  !>
  !! @brief Loads an atomic element from library.
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine AtomicElement_load( this, symbolSelected, massicNumber)
    implicit none
    
    type(AtomicElement), intent(inout) :: this
    character(*) :: symbolSelected
    integer :: massicNumber
    logical :: existFile
    integer :: stat
    integer :: auxMassicNumber
    integer :: i
    
    !! Namelist definition
    character(30) :: name
    character(30) :: symbol
    real(8) :: atomicNumber
    real(8) :: mass
    real(8) :: meltingPoint
    real(8) :: boilingPoint
    real(8) :: density
    real(8) :: electronAffinity
    real(8) :: ionizationEnergy1
    real(8) :: electronegativity    
    real(8) :: covalent
    real(8) :: atomic
    real(8) :: vanderWaals
    real(8) :: isotopes(4,30) !! asumming that no one atom have more than 30 isotopes...

    NAMELIST /element/ &
         name, &
         symbol, &
         atomicNumber, &
         mass, &
         meltingPoint, &
         boilingPoint, &
         density, &
         electronAffinity, & 
         ionizationEnergy1, & 
         electronegativity, &
         covalent, & 
         atomic, & 
         vanderWaals, &
         isotopes

    !! Looking for library    
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%ATOMIC_ELEMENTS_DATABASE), exist=existFile)
    
    if ( existFile ) then
       
       !! Open library
       open(unit=10, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%ATOMIC_ELEMENTS_DATABASE), status="old", form="formatted" )
              
       !! Read information
       symbol = "NONE"
       stat = 0

       rewind(10)
       
       do while(trim(symbol) /= trim(symbolSelected))
       
          !! Setting defaults
          name = "NONE"
          atomicNumber = 0
          mass = -1
          meltingPoint = 0
          boilingPoint = 0
          density = 0
          electronAffinity = 0
          ionizationEnergy1 = 0
          electronegativity = 0
          covalent = 0
          atomic = 0
          vanderWaals = 0
          isotopes = -1
          
          if (stat == -1 ) then
             
             this%isInstanced = .false.
             return
             
          end if
          
          read(10,NML=element, iostat=stat)
 
          if (stat > 0 ) then
             
             print*, "ERROR!!! ", stat
             call AtomicElement_exception( ERROR, "Failed reading AtomicElements.lib file!! please check this file.", "In AtomicElement at load function.")
             
          end if

       end do
       
       !! Set object variables
       this%name = name
       this%symbol = symbol
       this%atomicNumber = atomicNumber       
       this%meltingPoint = meltingPoint       
       this%boilingPoint = boilingPoint
       this%density = density
       this%electronAffinity =	electronAffinity
       this%electronegativity = electronegativity
       this%ionizationEnergy1 = ionizationEnergy1
       this%covalentRadius = covalent
       this%atomicRadio =atomic
       this%vanderWaalsRadio =	vanderWaals
       this%isMassNumberPresent = .false.
       this%abundance = 0
       this%isInstanced = .true.
       
       !! Set mass weight: The most abundant isotope if mass specification is not given.
       if ( massicNumber /= 0 ) then
          do i = 1, size(isotopes,DIM=2)
             if(isotopes(1,i) == massicNumber) then
                this%massicNumber = isotopes(1,i)
                this%atomicWeight = isotopes(2,i)
                this%abundance = isotopes(3,i)
                this%nuclearSpin = isotopes(4,i)
                this%isMassNumberPresent = .true.          
                exit
             end if
          end do

       else

          do i = 1, size(isotopes,DIM=2)
             if( isotopes(3,i) > this%abundance ) then
                this%massicNumber = isotopes(1,i)
                this%atomicWeight = isotopes(2,i)
                this%abundance = isotopes(3,i)
                this%nuclearSpin = isotopes(4,i)
                this%isMassNumberPresent = .false.
             end if
          end do

       end if
       
       !! Debug information.
       !! call AtomicElement_show(this)
       
       close(10)
       
    else 

       call AtomicElement_exception( ERROR, "LOWDIN library not found!! please export lowdinvars.sh file.", "In AtomicElement at load function.")

    end if 

    !! Done
    
  end subroutine AtomicElement_load
  
  !<
  !! @brief Define el destructor para clase
  subroutine AtomicElement_show( this )
    implicit none
    type(AtomicElement) , intent(in) :: this
    
    
    print *,""
    print *,"======================="
    print *,"  Element Properties   "
    print *,"======================="
    print *,""
    write (6,"(T10,A22,A12,A10)")	"Name                = ",this%name,""
    write (6,"(T10,A22,A12,A10)")	"Symbol              = ",this%symbol,""
    write (6,"(T10,A22,F12.5,A10)")	"Atomic number       = ",this%atomicNumber,""
    write (6,"(T10,A22,F12.5,A10)")	"masic number        = ",this%massicNumber,""
    write (6,"(T10,A22,F12.5,A10)")	"Isotopic abundance  = ",this%abundance,""
    write (6,"(T10,A22,F12.5,A10)")	"Melting point       = ",this%meltingpoint," K"
    write (6,"(T10,A22,F12.5,A10)")	"Boiling point       = ",this%boilingpoint," K"
    write (6,"(T10,A22,F12.5,A10)")	"Density             = ",this%density," g/cm^3"
    write (6,"(T10,A22,F12.5,A10)")	"Electron affinity   = ",this%electronaffinity," kJ/mol"
    write (6,"(T10,A22,F12.5,A10)")	"1.Ionization energy = ",this%ionizationenergy1," kJ/mol"
    write (6,"(T10,A22,F12.5,A10)")	"Electronegativity   = ",this%electronegativity," (pauling)"
    write (6,"(T10,A22,F12.5,A10)")	"Covalent radius     = ",this%covalentRadius," pm"
    write (6,"(T10,A22,F12.5,A10)")	"Atomic radius       = ",this%atomicRadio," pm"
    write (6,"(T10,A22,F12.5,A10)")	"Van der Waals radius= ",this%vanderwaalsRadio," pm"
    write (6,"(T10,A22,F12.5,A10)")	"Atomic weight       = ",this%atomicWeight,"u.m.a"
    write (6,"(T10,A22,F12.5,A10)")	"Nuclear Spin        = ",this%nuclearSpin,""
    print *,""
    
  end subroutine AtomicElement_show
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine AtomicElement_exception( typeMessage, description, debugDescription)
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
    
  end subroutine AtomicElement_exception
  
end module AtomicElement_
