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
!! @brief Module as "quantum species administrator"
!!
!!  This module handles all kind of both classic and quantum species for any system.
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-14
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapts module to inclusion on LOWDIN package
!!   - <tt> 2013-04-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrites the module to avoid XML dependence and for implement new LOWDIN standard.
!!
module Species_
  use CONTROL_
  use Exception_
  use Particle_
  use ConstantsOfCoupling_
  implicit none


  type , public :: Species
     
     character(10) :: name                     !< Name of species: electron, muons, etc.
     character(10) :: symbol                   !< Symbol of species: e-, u-, etc.
     character(10) :: statistics               !< Boson / fermion
     real(8) :: charge        !< Carga asociada a la especie.
     real(8) :: mass       !< Masa asociada a la particula.
     real(8) :: spin       !< Especifica el espin de la especie
     real(8) :: totalCharge       !< Carga total asociada a la especie.
     real(8) :: kappa
     real(8) :: eta
     real(8) :: lambda
     real(8) :: particlesFraction
     real(8) :: ocupationNumber
     real(8) :: multiplicity
     integer :: internalSize       !< Numero de particulas si se trata de una particula estructurada
     integer :: basisSetSize       
     integer :: speciesID
     logical(1) :: isElectron
     
     type(Particle), allocatable :: particles(:)
     
  end type Species

contains
  
  !>
  !! @brief Updates the information of the species from charged particles.
  !! @author E. F. Posada, 2013
  subroutine Species_setSpecie( this, speciesID)
    implicit none
    
    type(species) :: this
    integer :: speciesID
    
    type(ConstantsOfCoupling) :: couplingConstants    
    integer :: i

    !! Adjust coupling constants
    call ConstantsOfCoupling_load( couplingConstants, trim(this%symbol) )
    
    !! If the species was found on library
    if(couplingConstants%isInstanced) then
       
       this%kappa = couplingConstants%kappa
       this%eta = couplingConstants%eta
       this%lambda = couplingConstants%lambda
       this%particlesFraction = couplingConstants%particlesFraction
       
    else
       !! if not... Setting defaults by statistics
       if(trim(this%statistics) == "FERMION") then
          
          this%kappa = -1.0_8
          this%eta = 2.0_8
          this%lambda = 2.0_8
          this%particlesFraction = 0.5_8
          
       else 
          
          this%kappa = -1.0_8
          this%eta = 1.0_8
          this%lambda = 1.0_8
          this%particlesFraction = 1.0_8
          
       end if
       
    end if
    
    this%multiplicity = 0
    this%totalCharge = 0
    this%basisSetSize = 0
    this%internalSize = 0
    this%ocupationNumber = 0    
    
    do i = 1, size(this%particles)
       
       !! Set multiplicity
       this%multiplicity = this%multiplicity + ( this%particles(i)%multiplicity - 1 )
       !! Set totalCharge       
       this%totalCharge = this%totalCharge + this%particles(i)%totalCharge
       !! Calculate basis set size       
       this%basisSetSize = this%basisSetSize + size(this%particles(i)%basis%contraction)
       
       if( .not. this%particles(i)%isDummy ) then
          !! Set internal size
          this%internalSize = this%internalSize + this%particles(i)%internalSize
          !! Set ocupation number
          this%ocupationNumber = this%ocupationNumber + this%particles(i)%internalSize
       end if
       
    end do
    
    !! Set speciesID
    this%speciesID = speciesID
    !! Adjust name of species
    this%name = trim(this%particles(1)%name)
    !! Adjust symbol of species
    this%symbol = trim(this%particles(1)%symbol)    
    !! Adjust statistics
    this%statistics = trim(this%particles(1)%statistics)
    !! Adjust charge
    this%charge = this%particles(1)%charge
    !! Adjust mass
    this%mass = this%particles(1)%mass
    !! Adjust spin
    this%spin = this%particles(1)%spin
    !! Adjust multiplicity
    this%multiplicity = this%multiplicity + 1
    !! Adjust Occupation number
    this%ocupationNumber = this%ocupationNumber * this%particlesFraction
    
    !! are electrons?
    this%isElectron = .false.
    if(trim(this%particles(1)%symbol) == "E-") this%isElectron = .true.    
    if(trim(this%particles(1)%symbol) == "E-ALPHA") this%isElectron = .true.    
    if(trim(this%particles(1)%symbol) == "E-BETA") this%isElectron = .true.    
    
  end subroutine Species_setSpecie

  !>
  !! @brief Saves the info of the this species to file.
  !! @author E. F. Posada, 2013
  subroutine Species_saveToFile(this, unit)
    implicit none
    
    type(species), intent(in) :: this
    integer :: unit
    
    integer :: i
    
    write(unit,*) this%name
    write(unit,*) this%symbol
    write(unit,*) this%statistics
    write(unit,*) this%charge
    write(unit,*) this%mass
    write(unit,*) this%spin
    write(unit,*) this%totalCharge
    write(unit,*) this%kappa
    write(unit,*) this%eta
    write(unit,*) this%lambda
    write(unit,*) this%particlesFraction
    write(unit,*) this%ocupationNumber
    write(unit,*) this%multiplicity
    write(unit,*) this%internalSize
    write(unit,*) this%basisSetSize
    write(unit,*) this%speciesID
    write(unit,*) this%isElectron
    write(unit,*) size(this%particles)
    
    do i = 1, size(this%particles)
       call Particle_saveToFile(this%particles(i), unit)
    end do
    
  end subroutine Species_saveToFile

  !>
  !! @brief Load the info of the this species from file.
  !! @author E. F. Posada, 2013
  subroutine Species_loadFromFile(this, unit)
    implicit none
    
    type(species) :: this
    integer :: unit
    
    integer :: i
    integer :: numberOfParticles
    
    read(unit,*) this%name
    read(unit,*) this%symbol
    read(unit,*) this%statistics
    read(unit,*) this%charge
    read(unit,*) this%mass
    read(unit,*) this%spin
    read(unit,*) this%totalCharge
    read(unit,*) this%kappa
    read(unit,*) this%eta
    read(unit,*) this%lambda
    read(unit,*) this%particlesFraction
    read(unit,*) this%ocupationNumber
    read(unit,*) this%multiplicity
    read(unit,*) this%internalSize
    read(unit,*) this%basisSetSize
    read(unit,*) this%speciesID
    read(unit,*) this%isElectron
    read(unit,*) numberOfParticles

    allocate(this%particles(numberOfParticles))
    
    do i = 1, size(this%particles)
       call Particle_loadFromFile(this%particles(i), unit)
    end do
    
  end subroutine Species_loadFromFile  

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine Species_exception( typeMessage, description, debugDescription)
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

  end subroutine Species_exception

end module Species_
