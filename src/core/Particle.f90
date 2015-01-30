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
!! @brief Module to particles definition
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2007-02-06
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapto al estandar de codificacion propuesto.
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrittes the code to be included on LOWDIN package.
!!   - <tt> 2013-04-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrittes the code to eliminate the XML dependence
module Particle_
  use Exception_
  use CONTROL_
  use PhysicalConstants_
  use AtomicElement_
  use ElementalParticle_
  use BasisSet_  
  implicit none
  
  type, public :: Particle
     type(BasisSet) :: basis            !< Basis set for particle (if any)
     character(10) :: name              !< Name of particle: electron, muons, etc.
     character(10) :: symbol            !< Symbol of particle: e-, u-, etc.
     character(50) :: nickname          !< Name in input file: e-(H), U-, He_4, etc.
     character(10) :: statistics        !< Boson / fermion
     character(20) :: basisSetName      !< basis set name
     real(8) :: origin(3)		!< Posicion espacial
     real(8) :: charge 			!< Carga asociada a la particula.
     real(8) :: mass			!< Masa asociada a la particula.
     real(8) :: spin			!< Especifica el espin de la particula
     real(8) :: totalCharge		!< Carga total asociada a la particula.
     logical :: isQuantum		!< Indica comportamiento cuantico/puntual
     logical :: isDummy			!< Permite verificar si se trata de una particula virtual
     logical :: fixComponent(3)		!< Indica cual coordenada sera un parametro. (is fixed)
     logical :: isCenterOfOptimization	!< Especifica si la particula sera centro de optimizacion -Atributo requerido por conveniencia-
     integer :: multiplicity
     integer :: id			!< Indice de particula dentro del sistema
     integer :: internalSize		!< Numero de particulas si se trata de una particula estructurada
     integer :: owner			!< asocia un indice a la particula que la indentifica como un centro de referencia.
     integer :: basisSetSize		!< Este atributo es adicionado por conveniencia
     integer, allocatable :: childs(:)	!< Cuando la particula es un centro de optimizacion (es padre), almacena los Ids de sus hijas
  end type Particle

  public :: &
       Particle_load
  
contains

  !>
  !! @brief loads all particle information, particles can be elements, nuclei or elemental particles.
  !! @author S. A. Gonz\'alez
  !! @par History
  !!      -Adapted for open shell systems, 2011. E. F. Posada
  !!      -Re-written and  Verified, 2013. E. F. Posada
  !! @version 2.0
  subroutine Particle_load( this, name, baseName, origin, fix, multiplicity, addParticles, spin, id, charge )
    implicit none
    type(particle) :: this
    character(*), intent(in) :: name
    character(*), intent(in), optional :: baseName
    character(*), intent(in), optional :: fix
    character(*), intent(in), optional :: spin
    real(8), intent(in), optional :: origin(3)
    real(8), intent(in), optional :: multiplicity
    integer, intent(in), optional :: addParticles
    integer, intent(in) :: id
    real(8), intent(in), optional :: charge
    
    type(AtomicElement) :: element
    type(ElementalParticle) :: eparticle
    character(3) :: varsToFix
    character(5) :: massNumberString
    character(5) :: elementSymbol
    integer :: i
    integer :: j
    integer :: massNumber
    integer :: auxAdditionOfParticles
    real(8) :: auxOrigin(3)
    real(8) :: auxMultiplicity
    logical :: isDummy
    logical :: isElectron
    
    varsToFix=""
    if ( present(fix) ) varsToFix =trim(fix)

    auxOrigin=[0.0_8,0.0_8,0.0_8]
    if   ( present(origin) ) auxOrigin= origin
    
    auxMultiplicity=1.0_8
    if   ( present(multiplicity) ) auxMultiplicity=multiplicity
    
    auxAdditionOfParticles=0
    if   ( present(addParticles) ) auxAdditionOfParticles= addParticles

    !! Initialize some variables
    isDummy = .false.
    isElectron = .false.
    massNumber = 0
    elementSymbol = trim(name)
    massNumberString=""
    
    !!*******************************************************************************************
    !! Identify what kind of particle is
    !!*******************************************************************************************
    
    !! Scan if the particle is dummy. This must have * at name. ie. (He*)
    if (  scan( name, "*" ) /= 0 ) then
       elementSymbol = name(1: scan( name, "*" ) - 1 )
       isDummy =.true.
       CONTROL_instance%ARE_THERE_DUMMY_ATOMS = .true.
    end if
    
    !! Scan if electrons of an atom. At input file they must be apear as e-(H), H is any atom symbol.
    if ( scan( name, "E-" ) == 1 ) then
       !! Just to be sure that it is an electron.
       if(  scan( name, "[" ) /= 0 ) then          
          elementSymbol = name(scan( name, "[" ) + 1: scan( name, "]" ) - 1)
          isElectron = .true.          
       else !! Maybe it is another kind of particle, you never know.          
          elementSymbol=trim(name)          
       end if       
    end if
    
    !! Scan if nucleous with mass especification. ie, isotope specification. He_4 for instance.
    if (  scan( name, "_" ) /= 0 ) then       
       elementSymbol = name(1: scan( name, "_" ) - 1 )       
       !! Obtains the mass number
       massNumberString = name( scan( name, "_" ) + 1 : len( trim(name) ) )
       read( massNumberString , * ) massNumber
    end if

    !!*******************************************************************************************
    !! Load particles 
    !!*******************************************************************************************
    
    this%multiplicity = auxMultiplicity
    
    !! Obtains information of the atom., if any...
    call AtomicElement_load ( element, trim(elementSymbol), massicNumber=massNumber)
    
    !! Load information related with atoms
    if(element%isInstanced) then  !! it is true when elementSymbol variable is found as an atomic element.
       
       !! Load electrons
       if (isElectron) then

          !! Loads Particle information
          call Particle_build( this=this, &
               isQuantum= .true., &
               origin=auxOrigin, &
               basisSetName=trim(baseName), &
               elementSymbol=trim(elementSymbol), &
               isDummy= isDummy, &
               owner=id, &
               nickname=trim(element%symbol) )
          
          !! Setting remaining variables...
          if (present(spin)) then
             
             !! Alpha electrons
             if ( trim(spin) == "ALPHA") then
                
                this%name = "E-ALPHA"
                this%symbol = "E-ALPHA"
                this%charge = PhysicalConstants_ELECTRON_CHARGE
                this%mass = PhysicalConstants_ELECTRON_MASS
                this%totalCharge = -element%atomicNumber
                this%internalSize = element%atomicNumber + (auxMultiplicity-1)  + auxAdditionOfParticles             
                this%statistics="FERMION"
                this%spin = PhysicalConstants_SPIN_ELECTRON
                this%id = id    
                
             !! Beta electrons   
             else if ( trim(spin) == "BETA") then
                
                this%name = "E-BETA"
                this%symbol = "E-BETA"
                this%charge = PhysicalConstants_ELECTRON_CHARGE
                this%mass = PhysicalConstants_ELECTRON_MASS
                this%totalCharge = -element%atomicNumber
                this%internalSize = element%atomicNumber - (auxMultiplicity-1)  + auxAdditionOfParticles             
                this%statistics="FERMION"
                this%spin = PhysicalConstants_SPIN_ELECTRON
                this%id = id    
                
             end if
             
          else 
             
             !! Restricted case
             this%name = "E-"
             this%symbol = "E-"
             this%charge = PhysicalConstants_ELECTRON_CHARGE
             this%mass = PhysicalConstants_ELECTRON_MASS
             this%totalCharge = -element%atomicNumber
             this%internalSize = element%atomicNumber  + auxAdditionOfParticles
             this%statistics="FERMION"
             this%spin = PhysicalConstants_SPIN_ELECTRON
             this%id = id    
             
          end if
          
          call Particle_setComponentFixed( this, varsToFix )
          
       !! Load quantum nuclei (of an atomic element)
       else if( present(baseName) .and. trim(baseName) /= "DIRAC" .and. trim(baseName) /= "MM") then

          call Particle_build( this = this, &
               isQuantum = .true., &
               origin = auxOrigin,&
               basisSetName = trim(baseName), &
               elementSymbol = trim(elementSymbol), &
               isDummy = isDummy, &
               owner = id, &
               massNumber = int(element%massicNumber),&
               nickname = trim(name))

          call Particle_setComponentFixed( this, varsToFix )
          
          !! Setting remaining variables...
          this%name = name
          this%symbol = trim(name)
          this%charge = element%atomicNumber
          this%mass = element%massicNumber * PhysicalConstants_NEUTRON_MASS &
               + element%atomicNumber * (PhysicalConstants_PROTON_MASS - PhysicalConstants_NEUTRON_MASS)          
          this%totalCharge = element%atomicNumber
          this%internalSize = 1  + auxAdditionOfParticles
          this%id = id    
          
          !! Identify statistics for this particle
          if ( Math_isEven( INT( element%massicNumber ) ) ) then
             
             this%statistics = "BOSON"

             if ( Math_isEven( INT( element%atomicNumber ) ) ) then
                
                this%spin = 0.0_8
                
             else
                
                this%spin = element%nuclearSpin
                
             end if

          else

             this%statistics = "FERMION"
             this%spin = element%nuclearSpin
             
          end if
       
       !! Load Point charges (of an atomic element)
       else

          if(present(baseName) .and. trim(baseName) == "MM" ) then
             
             call Particle_build( this, &
                  isQuantum=.false., &
                  name=trim(elementSymbol), &
                  origin=auxOrigin, &
                  owner=id, &
                  nickname=trim(element%symbol))

             call Particle_setComponentFixed(this, varsToFix )

             this%basisSetName="MM"
             this%name = trim(name)
             this%symbol = trim(elementSymbol)

             if(CONTROL_instance%CHARGES_MM) then
                this%charge = charge
                this%totalCharge = charge
             else
                this%charge = element%atomicNumber
                this%totalCharge = element%atomicNumber
             end if
             this%mass = element%massicNumber * PhysicalConstants_NEUTRON_MASS &
                  + element%atomicNumber * (PhysicalConstants_PROTON_MASS - PhysicalConstants_NEUTRON_MASS)

             this%internalSize = 1 + auxAdditionOfParticles
             this%id = id

             if ( Math_isEven( INT( element%massicNumber ) ) ) then

                this%statistics = "BOSON"
                if ( Math_isEven( INT( element%atomicNumber ) ) ) then
                   this%spin = 0.0_8
                else
                   this%spin = element%nuclearSpin
                end if

             else

                this%statistics = "FERMION"
                this%spin = element%nuclearSpin

             end if


          else
          
             call Particle_build( this, &
                  isQuantum=.false., &
                  name=trim(elementSymbol), &
                  origin=auxOrigin, &
                  owner=id, &
                  nickname=trim(element%symbol))

             call Particle_setComponentFixed(this, varsToFix )

             this%basisSetName="DIRAC"
             this%name = trim(name)
             this%symbol = trim(elementSymbol)
             this%charge = element%atomicNumber
             this%mass = element%massicNumber * PhysicalConstants_NEUTRON_MASS &
                  + element%atomicNumber * (PhysicalConstants_PROTON_MASS - PhysicalConstants_NEUTRON_MASS)
             this%totalCharge = element%atomicNumber
             this%internalSize = 1 + auxAdditionOfParticles
             this%id = id

             if ( Math_isEven( INT( element%massicNumber ) ) ) then

                this%statistics = "BOSON"
                if ( Math_isEven( INT( element%atomicNumber ) ) ) then
                   this%spin = 0.0_8
                else
                   this%spin = element%nuclearSpin
                end if

             else

                this%statistics = "FERMION"
                this%spin = element%nuclearSpin

             end if
          
          end if
       end if
       
    !! Load quantum elemental particles (is not an atomic element)
    else if ( present(baseName) .and. trim(baseName) /= "DIRAC" .and. trim(baseName) /= "MM") then
       
       call ElementalParticle_load( eparticle, trim(name) )
       
       call Particle_build( this = this, &
            isQuantum = .true., &
            origin = auxOrigin, &
            basisSetName = trim(baseName), &
            elementSymbol = trim(elementSymbol), &
            isDummy = isDummy, &
            owner = id, &
            nickname = trim(eparticle%symbol))
       
       !! Setting remaining variables       
       this%name = eParticle%name
       this%symbol = eParticle%symbol
       this%charge = eParticle%charge
       this%mass = eParticle%mass
       this%totalCharge = eParticle%charge
       this%internalSize = 1  + auxAdditionOfParticles
       this%spin = eParticle%spin
       this%id = id
       
       !! Identifica la estadistica de la particula
       if ( abs( ceiling( abs( eParticle%spin ) ) - abs( eParticle%spin ) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
          this%statistics = "BOSON"
       else
          this%statistics = "FERMION"
       end if
       
   !! Load Point charges (elemental Particles)       
    else

       if(present(baseName) .and. trim(baseName) == "MM" ) then
          call ElementalParticle_load( eparticle, trim( name ) )

          call Particle_build( this, &
               isQuantum=.false., &
               name=trim(name), &
               origin=auxOrigin, &
               owner=id, &
               nickname=trim(eparticle%symbol))

          call Particle_setComponentFixed(this, varsToFix )

          this%basisSetName="MM"
          this%name = eParticle%name
          this%symbol = eParticle%symbol
          if(CONTROL_instance%CHARGES_MM) then
             this%charge = charge
             this%totalCharge = charge
          end if
          this%mass = eParticle%mass
          this%internalSize = 1 + auxAdditionOfParticles
          this%spin = eParticle%spin
          this%id = id
          
       else
          
          call ElementalParticle_load( eparticle, trim( name ) )

          call Particle_build( this, &
               isQuantum=.false., &
               name=trim(name), &
               origin=auxOrigin, &
               owner=id, &
               nickname=trim(eparticle%symbol))

          call Particle_setComponentFixed(this, varsToFix )

          this%basisSetName="DIRAC"
          this%name = eParticle%name
          this%symbol = eParticle%symbol
          this%charge = eParticle%charge
          this%mass = eParticle%mass
          this%totalCharge = eParticle%charge
          this%internalSize = 1 + auxAdditionOfParticles
          this%spin = eParticle%spin
          this%id = id
          
       end if
       
    end if

    !!*******************************************************************
    !! Debug information
    ! call Particle_show(this)
    
  end subroutine Particle_load
  
  
  !>
  !! @brief Build and load particle information.
  !! @param this Quantum particle or point charge
  !! @author S. A. Gonzalez (before known as Particle_constructor)
  subroutine Particle_build( this, name, symbol, basisSetName, elementSymbol, nickname, &
       origin, mass, charge, totalCharge, spin, &
       owner, massNumber, isQuantum, isDummy)

    implicit none
    
    type(Particle) , intent(inout) :: this
    character(*), optional, intent(in) :: name
    character(*), optional, intent(in) :: symbol
    character(*), optional, intent(in) :: basisSetName
    character(*), optional, intent(in) :: elementSymbol
    character(*), optional, intent(in) :: nickname
    real(8), optional, intent(in) :: origin(3)
    real(8), optional, intent(in) :: mass
    real(8), optional, intent(in) :: charge
    real(8), optional, intent(in) :: totalCharge
    real(8), optional, intent(in) :: spin
    integer, optional, intent(in) :: owner
    integer, optional, intent(in) :: massNumber
    logical, optional, intent(in) :: isQuantum
    logical, optional, intent(in) :: isDummy
    
    integer :: i
    type(Exception) :: ex              
              
    !!***************************************************************
    !! Setting defaults    
    this%isQuantum = .false.
    this%origin = 0.0_8
    this%mass = PhysicalConstants_ELECTRON_MASS
    this%charge = PhysicalConstants_ELECTRON_CHARGE
    this%totalCharge = PhysicalConstants_ELECTRON_CHARGE
    this%name = "ELECTRON"
    this%symbol = "E-"

    this%statistics = "FERMION"
    this%basisSetName = "DIRAC"
    this%spin = PhysicalConstants_SPIN_ELECTRON
    this%fixComponent =.false.
    this%owner = 0
    this%isDummy = .false.
    this%internalSize = 0
    this%isCenterOfOptimization = .true.
    this%nickname="NONE"

    !! Loads optional information    
    if ( present(origin) ) this%origin = origin
    if ( present( mass ) ) this%mass = mass
    if ( present(charge) ) this%charge = charge    
    if ( present(totalCharge)) then
       this%totalCharge = totalCharge
    else
       this%totalCharge = this%charge
    end if    
    if ( present(spin) )  this%spin = spin
    if ( present(isQuantum) ) this%isQuantum = isQuantum
    if ( present(isDummy) ) this%isDummy = isDummy
    if ( present(name) ) this%name=trim(name)
    if ( present(symbol) ) this%symbol=trim(symbol)
    if ( present(owner) ) this%owner=owner
    if ( present(basisSetName) ) this%basisSetName = trim(basisSetName)
    if ( present(nickname) ) this%nickname=trim(nickname)
    this%id = this%owner

    !!*******************************************************************
    !! Load Basis set
    if ( .not.present(basisSetName) .and. (this%isQuantum .eqv. .true.) ) then
       
       call Particle_exception( ERROR, "The basis set for this particle wasn't specified!!", &
            "Particle at load function")
       
    else if ( (this%isQuantum .eqv. .true.) .and. present(basisSetName) ) then
       
       if ( present(massNumber) ) then

          call BasisSet_load ( this%basis, "DEMON2K", trim(basisSetName), trim(nickName), origin=this%origin)
          
       else
          
          call BasisSet_load ( this%basis, "DEMON2K", trim(basisSetName), trim(elementSymbol), origin=this%origin)
          
       end if
       
       this%basisSetSize = this%basis%length
       
       if(this%isDummy) this%isCenterOfOptimization = .false.
       
    end if
        
  end subroutine Particle_build

  

  
  !>
  !! @brief Destroy the object
  !! @param this Quantum particle or point charge
  !! @author S. A. Gonzalez (before known as Particle_destructor)
  subroutine Particle_destroy( this )
    implicit none
    
    type(particle) :: this
    
    if(allocated(this%basis%contraction)) deallocate(this%basis%contraction)
    if(allocated(this%childs)) deallocate(this%childs)
    this%name = "NONE"
    this%symbol = "NONE"      
    this%nickname = "NONE"    
    this%statistics = "NONE"  
    this%basisSetName = "NONE"
    this%origin = 0.0_8
    this%charge = 0	
    this%mass = 0.0_8	
    this%spin = 0.0_8	
    this%totalCharge = 0	
    this%isQuantum = .false.
    this%isDummy = .false.	
    this%fixComponent = .false.	
    this%isCenterOfOptimization = .false.
    this%id = 0
    this%internalSize = 0
    this%owner = 0
    this%basisSetSize = 0
    
  end subroutine Particle_destroy

  !>
  !! @brief Muestra los atributos de la particula dada (this).
  !! @param this Particula cuantica o puntual.
  subroutine Particle_show( this )
    implicit none
    type(Particle) , intent(in) :: this
    integer :: i
    character(20) :: childsNumber
    
    if ( this%isQuantum ) then
       print *," Quantum particle"
    else
       print *," Puntual particle"
    end if
    
    print *,"================="
    print *,""
    write (6,"(T10,A16,A8)") "Name          : ",trim( this%name )
    write (6,"(T10,A16,A8)") "Symbol        : ",trim( this%symbol )
    write (6,"(T10,A16,A8)") "Nickname      : ",trim( this%nickname )
    write (6,"(T10,A16,I8)") "Particle ID   : ",this%id
    write (6,"(T10,A16,I8)") "Owner         : ",this%owner
    write (6,"(T10,A16,F8.2)") "Charge        : ",this%charge
    write (6,"(T10,A16,F8.2)") "Mass          : ",this%mass
    write (6,"(T10,A16,F8.2)") "Spin          : ",this%spin
    write (6,"(T10,A16,F8.2)") "Total charge  : ",this%totalCharge
    write (6,"(T10,A16,I8)") "Internal size : ",this%internalSize
    write (6,"(T10,A16,L8)") "Is dummy      : ",this%isDummy
    write (6,"(T10,A16,A8)") "Statistics    : ",trim( this%statistics )
    
    if ( this%isQuantum ) then
       write (6,"(T10,A16,A8)") "Basis set name: ",trim( this%basisSetName )
       write (6,"(T10,A16,I8)") "Basis set size: ", this%basisSetSize
    end if
    
    write (6,"(T10,A16,F8.2,F8.2,F8.2)") "Origin        : ",this%origin(1),this%origin(2),this%origin(3)
    
!     if ( allocated(this%childs) ) then
!        write(childsNumber, *)  size(this%childs)
!        write (6,"(T10,A16,"//trim(childsNumber)//"I8)") "Childs        : ",( this%childs(i),i=1,size(this%childs) )
!     end if
    
    print *,""
    
  end subroutine Particle_show
  
  !>
  !! @brief Muestra los atributos de una lista de particulas
  !! @param this lista de particulas
  subroutine Particle_showStructure( this )
    implicit none
    type(Particle) , intent(in) :: this
    
    call Particle_show( this )
    
    if ( this%isQuantum	) then
       
       print *,""
       print *,"INFORMATION OF BASIS SET: "
       
       call BasisSet_showInCompactForm( this%basis, trim(this%nickname))
       
    end if
    
    
  end subroutine Particle_showStructure
  
  !>
  !! @brief Saves the particle structure to file.
  !! @param this Particula cuantica o puntual.
  !! @author E. F. Posada, 2013
  subroutine Particle_saveToFile( this, unit )
    implicit none
    
    type(Particle) , intent(in) :: this
    integer :: unit

    integer :: i
    logical :: childs

    
    write(unit,*) this%name
    write(unit,*) this%symbol
    write(unit,*) this%nickname
    write(unit,*) this%statistics
    write(unit,*) this%basisSetName
    write(unit,*) this%origin
    write(unit,*) this%charge
    write(unit,*) this%mass
    write(unit,*) this%spin
    write(unit,*) this%totalCharge
    write(unit,*) this%isQuantum
    write(unit,*) this%isDummy
    write(unit,*) this%fixComponent
    write(unit,*) this%isCenterOfOptimization
    write(unit,*) this%multiplicity
    write(unit,*) this%id
    write(unit,*) this%internalSize
    write(unit,*) this%owner
    write(unit,*) this%basisSetSize
    
    if ( allocated(this%childs) ) then
       childs = .true.
       write(unit,*) childs
       write(unit,*) size(this%childs)
       write(unit,*) this%childs
    else
       childs = .false.
       write(unit,*) childs
    end if
    
    !!Save basis set
    if(this%isQuantum) then
       call BasisSet_saveToFile(this%basis, unit)
    end if
    
  end subroutine Particle_saveToFile


  !>
  !! @brief Loads the particle structure from file.
  !! @param this Particula cuantica o puntual.
  !! @author E. F. Posada, 2013
  subroutine Particle_loadFromFile( this, unit )
    implicit none
    
    type(Particle) :: this
    integer :: unit

    integer :: i
    integer :: childsSize
    logical :: childs
    
    
    
    read(unit,*) this%name
    read(unit,*) this%symbol
    read(unit,*) this%nickname
    read(unit,*) this%statistics
    read(unit,*) this%basisSetName
    read(unit,*) this%origin
    read(unit,*) this%charge
    read(unit,*) this%mass
    read(unit,*) this%spin
    read(unit,*) this%totalCharge
    read(unit,*) this%isQuantum
    read(unit,*) this%isDummy
    read(unit,*) this%fixComponent
    read(unit,*) this%isCenterOfOptimization
    read(unit,*) this%multiplicity
    read(unit,*) this%id
    read(unit,*) this%internalSize
    read(unit,*) this%owner
    read(unit,*) this%basisSetSize
    read(unit,*) childs
    
    if (childs ) then
       
       read(unit,*) childsSize
       allocate(this%childs(childsSize))
       read(unit,*) this%childs
       
    end if
    
    !!loads basis set
    if(this%isQuantum) then
       call BasisSet_loadFromFile(this%basis, unit)
    end if
    
  end subroutine Particle_loadFromFile

  !>
  !! Devuelve un apuntador a una particula solicitada.
  !! @param this Lista de particulas
  function Particle_getParticlePtr(this) result(output)
    implicit none
    type(Particle) , target :: this
    type(Particle) , pointer :: output
    
    output => null()
    
    !! Devuelve el apuntador a la particula solicitada.
    output => this
    
  end function Particle_getParticlePtr

  !>
  !! Devuelve el name asignado a una particula solicitada.
  !! @param this Lista de particulas
  function Particle_getName(this) result(output)
    implicit none
    type(Particle) , intent(in) :: this
    character(30) :: output
    
    output = trim( this%name )
    
  end function Particle_getName

  !>
  !! Devuelve el indice asignado a la particula
  function Particle_getOwner(this) result(output)
    implicit none
    type(Particle) , intent(in) :: this
    integer :: output

    output= this%owner
    
  end function Particle_getOwner


  !>
  !! Devuelve la componente cartesiana del origen de la particula considerada
  !! @param this Particula cuantica o carga puntual.
  !! @param component componente x, y, o z del origen.
  function Particle_getOrigin(this) result(output)
    implicit none
    type(Particle) , target :: this
    real(8) :: output(3)
    
    output = this%origin
    
  end function Particle_getOrigin


  !>
  !! Devuelve la componente cartesiana del origen de la particula considerada
  !! @param this Particula cuantica o carga puntual.
  !! @param component componente x, y, o z del origen.
  function Particle_getOriginComponent(this,component) result(output)
    implicit none
    type(Particle) , target :: this
    integer :: component
    real(8) :: output
    
    output = this%origin(component)
    
  end function Particle_getOriginComponent
  
  
  !>
  !! Devuelve el comportamiento asignado a la particula, (TRUE/FALSE).
  !! @param this Particula cuantica o carga puntual.
  function Particle_isQuantum(this) result(output)
    implicit none
    type(Particle) , intent(in) :: this
    logical :: output
    
    output = this%isQuantum
    
  end function Particle_isQuantum

  !>
  !! Devuelve la masa de la particula solicitada.
  !! @param this Particula cuantica o carga puntual.
  function Particle_getMass(this) result(output)
    implicit none
    type(Particle) , intent(in) :: this
    real(8) :: output
    
    output = this%mass
    
  end function Particle_getMass

  !>
  !! Devuelve la carga de la particula solicitada
  !! @param this Particula cuantica o carga puntual.
  function Particle_getCharge( this ) result(output)
    implicit none
    type(Particle) , intent(in) :: this
    real(8) :: output
    
    output = this%charge
    
  end function Particle_getCharge
  
  
  !>
  !! Devuelve el espin de la particula solicitada.
  !! @param this Particula cuantica o carga puntual.
  function Particle_getSpin( this ) result( output )
    implicit none
    type(Particle) , intent(in) :: this
    real(8) :: output
    
    if (this%isQuantum .eqv. .true.) then
       
       output=this%spin
    end if
    
  end function Particle_getSpin

  !>
  !! Devuelve el identificador unico de la particula solicitada.
  !! @param this Particula cuantica o carga puntual.
  function Particle_getId( this ) result( output )
    implicit none
    type(Particle) , intent(in) :: this
    integer output
    
    output=this%id
    
  end function Particle_getId


  !>
  !! Devuelve el numero de contracciones asociadas a una particula.
  !! @param this Particula cuantica o carga puntual.
  function Particle_getNumberOfContractions( this ) result( output )
    implicit none
    type(Particle) , intent(in) :: this
    integer :: output
    
    output = 0
    
    if ( this%isQuantum .eqv. .true.)  then
       
       output = this%basisSetSize
       
    end if
    
  end function Particle_getNumberOfContractions
  
  !>
  !! Devuelve un apuntador asociado a una contraccion, de una particula cuantica.
  !! @param this Particula cuantica.
  function Particle_getContractionPtr( this, id ) result( output )
    implicit none
    type(Particle), target, intent(inout) :: this
    integer :: id
    type(ContractedGaussian) , pointer :: output
    type(Exception) :: ex
    
    
    if ( this%isQuantum ) then
       
       output => null()
       
       if ( allocated(this%basis%contraction) ) then
          
          output => this%basis%contraction(id)
          
       end if
       
    end if
    
    
  end function Particle_getContractionPtr
  
    
  !>
  !! Indica si la coordenada especificada se ha definido como variable
  !! o como parametro.
  !! @param this Particula cuantica o carga puntual.
  !! @param component Componente cartesiana x=1, y=2, z =3.
  function Particle_hasFixedComponent(this , component) result( output )
    implicit none
    type(Particle) , intent(inout) :: this
    integer , intent(in) :: component
    logical :: output
    
    output = this%fixComponent(component)
    
  end function Particle_hasFixedComponent

  
  !>
  !! Ajusta el id de las particulas hijas.
  !! @param this Particula cuantica o carga puntual.
  !>
  subroutine Particle_setChild( this, child )
    implicit none
    type(Particle) , intent(inout) :: this
    integer  :: child
    
    integer, allocatable :: auxChilds(:)
    integer :: ssize
    
    
    if ( allocated(this%childs) ) then
       ssize =size( this%childs )
       allocate( auxChilds( ssize  ) )
       auxChilds=this%childs
       deallocate(this%childs)
       allocate( this%childs( ssize +1  ) )
       this%childs(1:ssize)=auxChilds
       this%childs(ssize+1)=child
       deallocate(auxChilds)
    else
       allocate( this%childs(1))
       this%childs(1) = child
    end if
    
  end subroutine Particle_setChild
  

  !>
  !! @brief Remueve los Id de las particulas hijas
  !! @param this Particula cuantica o carga puntual.
  subroutine Particle_removeChilds( this )
    implicit none
    type(Particle) , intent(inout) :: this
    
    
    if ( allocated(this%childs) ) then
       deallocate(this%childs)
    end if
    
  end subroutine Particle_removeChilds

  !>
  !! @brief Ajusta el propieario de la particula especificada.
  !! @param this Particula cuantica o carga puntual.
  !>
  subroutine Particle_setOwner( this, owner )
    implicit none
    type(Particle) , intent(inout) :: this
    integer  :: owner

    integer :: i
    
    if (this%isQuantum .eqv. .true.) then
       
       do i = 1, size(this%basis%contraction)
          this%basis%contraction(i)%owner = owner
       end do

       this%owner = owner
       
    else
       
       this%owner = owner
       
    end if
    
  end subroutine Particle_setOwner

  
  !>
  !! @brief Permiter ajustar la componente de una particula como parametro o como
  !! variable durante la optimizacion de geometria. Por defecto lo
  !! ajusta como parametro.
  !! @param this Particula cuyo atributo se desea modificar
  !! @param component Componente cartesiana que se fija como parametro
  subroutine Particle_setComponentFixed(this , component)
    implicit none
    type(Particle) , intent(inout) :: this
    character(*) , intent(in) :: component
    
    select case ( trim( component ) )
       
       !!
       !! Ajusta componente x como parametro
       !!
    case ("x")
       this%fixComponent(1)= .true.
       
       !!
       !! Ajusta componente y como parametro
       !!
    case ("y")
       this%fixComponent(2)= .true.
       
       !!
       !! Ajusta componente z como parametro
       !!
    case ("z")
       this%fixComponent(3)= .true.
       
       !!
       !! Ajusta componente x como parametro
       !!
    case ("xy")
       this%fixComponent(1)= .true.
       this%fixComponent(2)= .true.
       
       !!
       !! Ajusta componente x como parametro
       !!
    case ("xz")
       this%fixComponent(1)= .true.
       this%fixComponent(3)= .true.
       
       !!
       !! Ajusta componente x como parametro
       !!
    case ("yz")
       this%fixComponent(2)= .true.
       this%fixComponent(3)= .true.
       
       !!
       !! Ajusta todas las componente cartesianas como parametro
       !!
    case ("xyz")
       this%fixComponent= .true.
       
       !!
       !! Deja que las coordenadas de la particula cambien durante
       !! la optimizacion de geometria.
       !!
    case default
       this%fixComponent= .false.
       
    end select
    
  end subroutine Particle_setComponentFixed  
  
  !>
  !! @brief  Handles exceptions for this class
  subroutine Particle_exception( typeMessage, description, debugDescription)
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
    
  end subroutine Particle_exception
  
end module Particle_

