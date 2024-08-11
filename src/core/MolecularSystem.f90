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
!! @brief This module handles all molecular system
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
module MolecularSystem_
  use CONTROL_
  use Exception_
  use Units_
  use Particle_
  use Species_
  use MecanicProperties_
  use Matrix_
  use Vector_
  use InternalCoordinates_
  use ExternalPotential_
  use InterPotential_
  implicit none
  
  type , public :: MolecularSystem
     
     character(100) :: name
     character(100) :: description
     integer :: numberOfParticles
     integer :: numberOfPointCharges
     integer :: numberOfQuantumParticles
     integer :: numberOfQuantumSpecies
     integer :: charge
          
     type(Species), allocatable :: species(:)
     type(particle), allocatable :: pointCharges(:)
     type(particleManager), allocatable :: allParticles(:)
     type(InternalCoordinates) :: intCoordinates

     type(MecanicProperties) :: mechanicalProp
     
  end type MolecularSystem

  public :: &
       MolecularSystem_initialize, &
       MolecularSystem_build, &
       MolecularSystem_destroy, &
       MolecularSystem_moveToCenterOfMass, &
       MolecularSystem_changeOriginOfSystem, &
       MolecularSystem_rotateOnPrincipalAxes, &
       MolecularSystem_showInformation, &
       MolecularSystem_showParticlesInformation, &
       MolecularSystem_showCartesianMatrix, &
       MolecularSystem_showZMatrix, &
       MolecularSystem_showDistanceMatrix, &
       MolecularSystem_saveToFile, &
       MolecularSystem_loadFromFile, &
       MolecularSystem_getNumberOfQuantumSpecies, &
       MolecularSystem_getNumberOfParticles, &
       MolecularSystem_getNumberOfContractions, &
       MolecularSystem_getTotalNumberOfContractions, &
       MolecularSystem_getBasisSet, &
       MolecularSystem_getMaxAngularMoment, &
       MolecularSystem_getMaxNumberofPrimitives, &
       MolecularSystem_getMaxNumberofCartesians, &
       MolecularSystem_getOcupationNumber, &
       MolecularSystem_getCharge, &
       MolecularSystem_getMass, &
       MolecularSystem_getEta, &
       MolecularSystem_getLambda, &
       MolecularSystem_getKappa, &
       MolecularSystem_getMultiplicity, &
       MolecularSystem_getParticlesFraction, &
       MolecularSystem_getFactorOfExchangeIntegrals, &
       MolecularSystem_getNameOfSpecie, &
       MolecularSystem_getNameOfSpecies, &
       MolecularSystem_getSpecieID, &
       MolecularSystem_getSpecieIDFromSymbol, &
       MolecularSystem_getPointChargesEnergy, &
       MolecularSystem_getMMPointChargesEnergy, &
       MolecularSystem_getlabelsofcontractions, &
       MolecularSystem_changeOrbitalOrder, &
       MolecularSystem_readFchk, &
       MolecularSystem_copyConstructor, &
       MolecularSystem_mergeTwoSystems, &
       MolecularSystem_GetTwoSystemsDisplacement, &
       MolecularSystem_checkParticleEquivalence,&
       MolecularSystem_getTotalMass
  
  !>Singleton
  type(MolecularSystem), public, target :: MolecularSystem_instance
  
contains
  
  !>
  !! @brief initializes the molecular system.
  !! @author E. F. Posada, 2013
  subroutine MolecularSystem_initialize(numberOfQuantumSpecies, numberOfPointCharges, numberOfParticlesForSpecies, quantumSpeciesName, systemName, systemDescription)
    implicit none
    
    integer :: numberOfQuantumSpecies
    integer :: numberOfPointCharges
    integer, allocatable :: numberOfParticlesForSpecies(:)
    character(15), allocatable :: quantumSpeciesName(:)
    character(*) :: systemName
    character(*) :: systemDescription

    integer :: i, j
    integer :: counter

    !! Setting defaults
    MolecularSystem_instance%name = trim(systemName)
    MolecularSystem_instance%description = trim(systemDescription)
    MolecularSystem_instance%numberOfPointCharges = numberOfPointCharges
    MolecularSystem_instance%numberOfQuantumSpecies = numberOfQuantumSpecies
    MolecularSystem_instance%numberOfParticles = sum(numberOfParticlesForSpecies) + numberOfPointCharges

    !! Allocate memory for particles
    allocate(MolecularSystem_instance%species(MolecularSystem_instance%numberOfQuantumSpecies))

    do i = 1, numberOfQuantumSpecies
       allocate(MolecularSystem_instance%species(i)%particles(numberOfParticlesForSpecies(i)))
       MolecularSystem_instance%species(i)%symbol = trim(quantumSpeciesName(i))
    end do

    allocate(MolecularSystem_instance%pointCharges(MolecularSystem_instance%numberOfPointCharges))
    allocate(molecularSystem_instance%allParticles(MolecularSystem_instance%numberOfParticles))
    
    
    !! Set the particles manager (all pointers)        
    counter = 1
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       do j = 1, size(MolecularSystem_instance%species(i)%particles)
          
          molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%species(i)%particles(j)
          counter = counter + 1

       end do
    end do

    do i = 1, MolecularSystem_instance%numberOfPointCharges
       
       molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%pointCharges(i)
       counter = counter + 1

    end do

    particleManager_instance => molecularSystem_instance%allParticles
    
  end subroutine MolecularSystem_initialize
  
  !>
  !! @brief Builds the molecular system.
  !! @author E. F. Posada, 2013
  subroutine MolecularSystem_build()
    implicit none
    
    integer :: i
    
    !! Setting quantum species
    MolecularSystem_instance%numberOfQuantumParticles = 0
    MolecularSystem_instance%charge = 0
   
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       
       call Species_setSpecie(MolecularSystem_instance%species(i), speciesID = i)
       
       if( MolecularSystem_instance%species(i)%isElectron .and. CONTROL_instance%IS_OPEN_SHELL ) then
          
          MolecularSystem_instance%charge = MolecularSystem_instance%charge + int(MolecularSystem_instance%species(i)%totalCharge / 2)
          MolecularSystem_instance%numberOfQuantumParticles = MolecularSystem_instance%numberOfQuantumParticles + int(MolecularSystem_instance%species(i)%internalSize / 2)
          
       else
          
          MolecularSystem_instance%charge = MolecularSystem_instance%charge + MolecularSystem_instance%species(i)%totalCharge
          MolecularSystem_instance%numberOfQuantumParticles = MolecularSystem_instance%numberOfQuantumParticles + MolecularSystem_instance%species(i)%internalSize
          
       end if

       !!Check for input errors in the number of particles
       if( (abs(int(MolecularSystem_instance%species(i)%ocupationNumber)-MolecularSystem_instance%species(i)%ocupationNumber) .gt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD)) then
          print *, "species ", trim(MolecularSystem_getNameOfSpecie(i)) , "has fractional ocupation number ", &
               MolecularSystem_instance%species(i)%ocupationNumber, "please check your input addParticles and multiplicity"
          call MolecularSystem_exception(ERROR, "Fractional ocupation number, imposible combination of charge and multiplicity","MolecularSystem module at build function.")
       end if
    end do
    
    do i = 1, MolecularSystem_instance%numberOfPointCharges
       
       MolecularSystem_instance%charge = MolecularSystem_instance%charge + MolecularSystem_instance%pointCharges(i)%totalCharge

    end do
    
    call ParticleManager_setOwner()
    
    !! Debug
    !! call ParticleManager_show()

  end subroutine MolecularSystem_build

  !>
  !! @brief Destry the molecular system object.
  !! @author E. F. Posada, 2013
  subroutine MolecularSystem_destroy()
    implicit none
    
    MolecularSystem_instance%name = "NONE"
    MolecularSystem_instance%description = "NONE"
    MolecularSystem_instance%numberOfParticles = 0
    MolecularSystem_instance%numberOfPointCharges = 0
    MolecularSystem_instance%numberOfQuantumParticles = 0
    MolecularSystem_instance%numberOfQuantumSpecies = 0
    MolecularSystem_instance%charge = 0
    
    if(allocated(MolecularSystem_instance%species)) deallocate(MolecularSystem_instance%species)
    if(allocated(MolecularSystem_instance%pointCharges)) deallocate(MolecularSystem_instance%pointCharges)
    if(allocated(MolecularSystem_instance%allParticles)) deallocate(MolecularSystem_instance%allParticles)
    
    call MecanicProperties_destructor(MolecularSystem_instance%mechanicalProp)

    call ExternalPotential_destructor()
    call InterPotential_destructor()


  end subroutine MolecularSystem_destroy

  !>
  !! @brief Traslada el origen del sistema molecular al centro de masa
  subroutine MolecularSystem_moveToCenterOfMass()
    implicit none

    call MecanicProperties_getCenterOfMass( MolecularSystem_instance%mechanicalProp )
    call MolecularSystem_changeOriginOfSystem( MolecularSystem_instance%mechanicalProp%centerOfMass )
    
  end subroutine MolecularSystem_moveToCenterOfMass

  !>
  !! @brief Traslada el origen del sistema molecular al centro de masa  
  subroutine MolecularSystem_changeOriginOfSystem( origin )
    implicit none
    real(8) :: origin(3)
    real(8) :: auxOrigin(3)
    
    integer :: i
    
    do i=1, size(MolecularSystem_instance%allParticles)
       
       auxOrigin = MolecularSystem_instance%allParticles(i)%particlePtr%origin
       auxOrigin=auxOrigin-origin
       call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(i)%particlePtr, auxOrigin )
       
    end do

  end subroutine MolecularSystem_changeOriginOfSystem
  
  !>
  !! @brief Traslada el origen del sistema molecular al centro de masa
  subroutine MolecularSystem_rotateOnPrincipalAxes()
    implicit none
    
    type(Matrix) :: matrixOfCoordinates
    type(Vector) :: vectorOfCoordinates
    real(8), allocatable :: auxVector(:)
    real(8) :: coordinates(3)
    integer :: i
    integer :: numOfzeros(2)
    
    
    MolecularSystem_instance%mechanicalProp%molecularInertiaTensor = MecanicProperties_getMolecularInertiaTensor( MolecularSystem_instance%mechanicalProp )
    matrixOfCoordinates = ParticleManager_getCartesianMatrixOfCentersOfOptimization()

    do i=1, size(matrixOfCoordinates%values,dim=1)
       
       !!
       !! Proyecta las coordenadas cartesianas sobre el tensor de inercia molecular
       !!
       coordinates(1)=dot_product(matrixOfCoordinates%values(i,:),MolecularSystem_instance%mechanicalProp%molecularInertiaTensor%values(:,3) )
       coordinates(2)=dot_product(matrixOfCoordinates%values(i,:),MolecularSystem_instance%mechanicalProp%molecularInertiaTensor%values(:,2) )
       coordinates(3)=dot_product(matrixOfCoordinates%values(i,:),MolecularSystem_instance%mechanicalProp%molecularInertiaTensor%values(:,1) )
       matrixOfCoordinates%values(i,:) = coordinates
    end do

    numOfzeros=0
    do i=1, size(matrixOfCoordinates%values,dim=1)
       if( abs(matrixOfCoordinates%values(i,3)) < 1.0D-6  ) numOfzeros(1) = numOfzeros(1)+1
       if( abs(matrixOfCoordinates%values(i,2)) < 1.0D-6  ) numOfzeros(2) = numOfzeros(2)+1
    end do

    if ( numOfzeros(1) > numOfzeros(2) ) then
       allocate( auxVector( size(matrixOfCoordinates%values,dim=1) ) )
       auxVector=matrixOfCoordinates%values(:,2)
       matrixOfCoordinates%values(:,2) = matrixOfCoordinates%values(:,3)
       matrixOfCoordinates%values(:,3) = auxVector
       deallocate( auxVector )
    end if
    
    !! Invierte la orientacion sobre el eje Z
    matrixOfCoordinates%values(:,3)=-1.0*matrixOfCoordinates%values(:,3)
    
    call Vector_constructor( vectorOfCoordinates, size( matrixOfCoordinates%values,dim=1)*3 )
    do i=1, size(matrixOfCoordinates%values,dim=1)
       vectorOfCoordinates%values(3*i-2:3*i)=matrixOfCoordinates%values(i,:)
    end do
    
    !! Realiza la rotacion de las pariculas
    call ParticleManager_setParticlesPositions( vectorOfCoordinates )
    
    call Vector_destructor( vectorOfCoordinates )
    call Matrix_destructor( matrixOfCoordinates )
    
  end subroutine MolecularSystem_rotateOnPrincipalAxes
  
  !>
  !! @brief shows general information of the molecular system.
  !! @author S. A. Gonzalez
  subroutine MolecularSystem_showInformation(this)
    type(MolecularSystem), optional, target :: this
    
    type(MolecularSystem), pointer :: system

    if( present(this) ) then
       system=>this
    else
       system=>MolecularSystem_instance
    end if
        
    print *,""
    print *," MOLECULAR SYSTEM: ",trim(system%name)
    print *,"-----------------"
    print *,""
    write (6,"(T5,A16,A)") "DESCRIPTION   : ", trim( system%description )
    write (6,"(T5,A16,I3)") "CHARGE        : ",system%charge
    write (6,"(T5,A16,A4)") "PUNTUAL GROUP : ", "NONE"
    print *,""
    
 
  end subroutine MolecularSystem_showInformation

  !>
  !! @brief Muestra los atributos de todas las particulas en el Administrador de particulas
  !! @author S. A. Gonzalez
  !! @par changes : 
  !!     - rewritten, E. F. Posada. 2013
  subroutine MolecularSystem_showParticlesInformation(this)
    implicit none
    type(MolecularSystem), optional, target :: this
    
    type(MolecularSystem), pointer :: system
    integer :: i, j

    if( present(this) ) then
       system=>this
    else
       system=>MolecularSystem_instance
    end if

    print *,""
    print *," INFORMATION OF PARTICLES :"
    print *,"========================="
    write (6,"(T10,A29,I8.0)") "Total number of particles   =", system%numberOfQuantumParticles + system%numberOfPointCharges
    write (6,"(T10,A29,I8.0)") "Number of quantum particles =", system%numberOfQuantumParticles
    write (6,"(T10,A29,I8.0)") "Number of puntual charges   =", system%numberOfPointCharges
    write (6,"(T10,A29,I8.0)") "Number of quantum species   =", system%numberOfQuantumSpecies

    !!***********************************************************************
    !! Imprime iformacion sobre masa, carga y numero de particulas encontradas
    !!
    print *,""
    print *,"                INFORMATION OF QUANTUM SPECIES "
    write (6,"(T5,A70)") "---------------------------------------------------------------------------------------------"
    write (6,"(T10,A2,A4,A8,A12,A4,A5,A6,A5,A6,A5,A4,A5,A12)") "ID", " ","Symbol", " ","mass", " ","charge", " ","omega","","spin","","multiplicity"
    write (6,"(T5,A70)") "---------------------------------------------------------------------------------------------"

    do i = 1, system%numberOfQuantumSpecies
       write (6,'(T8,I3.0,A5,A10,A5,F10.4,A5,F5.2,A5,F5.2,A5,F5.2,A5,F5.2)') &
            i, " ", &
            trim(system%species(i)%symbol)," ",&
            system%species(i)%mass," ",&
            system%species(i)%charge, " ",&
            system%species(i)%omega," ",&
            system%species(i)%spin, "",&
            system%species(i)%multiplicity
    end do

    print *,""
    print *,"                  CONSTANTS OF COUPLING "
    write (6,"(T7,A60)") "------------------------------------------------------------"
    write (6,"(T10,A11,A11,A11,A11,A11)") "Symbol", "kappa","eta","lambda","occupation"
    write (6,"(T7,A60)") "------------------------------------------------------------"
    
    do i = 1, system%numberOfQuantumSpecies
       write (6,'(T10,A11,F11.2,F11.2,F11.2,F11.2)') &
       trim(system%species(i)%symbol),&
            system%species(i)%kappa,&
            system%species(i)%eta,&
            system%species(i)%lambda,&
            system%species(i)%ocupationNumber
    end do

    print *,""

    print *,"                  BASIS SET FOR SPECIES "
    write (6,"(T7,A60)") "------------------------------------------------------------"
    write (6,"(T10,A8,A10,A8,A5,A12,A5,A9)") "Symbol", " ","N. Basis", " ","N. Particles"," ","Basis Set"
    write (6,"(T7,A60)") "------------------------------------------------------------"
    
    !! Only shows the basis-set name of the first particle by specie.
    do i = 1, system%numberOfQuantumSpecies
       
       if( system%species(i)%isElectron .and. CONTROL_instance%IS_OPEN_SHELL ) then

          write (6,'(T10,A10,A5,I8,A5,I12,A5,A10)') &
               trim(system%species(i)%symbol)," ",&
               !MolecularSystem_getTotalNumberOfContractions(i)," ",&
               system%species(i)%basisSetSize," ",&
               int(system%species(i)%internalSize / 2), " ",&
               trim(system%species(i)%particles(1)%basis%name)
          ! write(*,*)MolecularSystem_getTotalNumberOfContractions(i)
          
       else

          write (6,'(T10,A10,A5,I8,A5,I12,A5,A10)') &
               trim(system%species(i)%symbol)," ",&
               !MolecularSystem_getTotalNumberOfContractions(i)," ",&
               system%species(i)%basisSetSize," ",&
               system%species(i)%internalSize, " ",&
               trim(system%species(i)%particles(1)%basis%name)
          ! write(*,*)MolecularSystem_getTotalNumberOfContractions(i)
          
       end if
       
    end do

    write (6,*) ""
    write (6,"(T10,A35)")"                     ATOMIC BASIS"
    write (6,"(T10,A60)") "------------------------------------------------------------"
    write (6,"(T10,A11,A9,A20,A20)") " PRIMITIVE ", "  SHELL  "," EXPONENT "," COEFFICIENT "
    write (6,"(T10,A60)") "------------------------------------------------------------"
    
    do i = 1, system%numberOfQuantumSpecies
       
       !! Avoid print twice basis in open-shell case
       if(trim(system%species(i)%name) == "E-BETA" ) cycle
       
       write (6,*) ""
       write( 6, "(T5,A32,A5)") "BEGIN DESCRIPTION OF BASIS FOR: ", trim(system%species(i)%symbol)
       write (6,"(T5,A30)") "================================"
       write (6,*) ""
       
       do j = 1, size( system%species(i)%particles )
          
          call BasisSet_showInCompactForm( system%species(i)%particles(j)%basis,&
               trim(system%species(i)%particles(j)%nickname ))
          
       end do
       
       write (6,"(T5,A28)") "... END DESCRIPTION OF BASIS"
       write (6,*) ""
    end do
    
    print *,""
    print *," END INFORMATION OF PARTICLES"
    print *,""

    !!***********************************************************************
    !! Prints information about number of occupied orbitales and basis set size
    !!
    print *,""
    print *," INFORMATION OF THE QUANTUM SYSTEM "
    write (6,"(T5,A70)") "---------------------------------------------------------------------"
    write (6,"(T10,A10,A5,A17,A5,A10)") "ID", " ", "Occupied Orbitals", " " ,"Basis size"
    write (6,"(T5,A70)") "---------------------------------------------------------------------"

    do i = 1, system%numberOfQuantumSpecies
          write (6,'(T10,A10,A5,I8,A5,I12)') &
               trim(system%species(i)%symbol)," ",&
                MolecularSystem_getOcupationNumber( i )," ",&
               MolecularSystem_getTotalNumberOfContractions(i)
    end do

    print *,""
    print *," END INFORMATION OF QUANTUM SYSTEM"
    print *,""


    if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
      print *,""
      print *," INFORMATION OF EXTERNAL POTENTIALS "
      call ExternalPotential_show()
      print *,""
      print *," END INFORMATION OF EXTERNAL POTENTIALS"
      print *,""
    end if

    if(CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL) then
      print *,""
      print *," INFORMATION OF INTER-PARTICLE POTENTIALS "
      call InterPotential_show()
      print *,""
      print *," END INFORMATION OF INTER-PARTICLE POTENTIALS"
      print *,""
    end if
 
  end subroutine MolecularSystem_showParticlesInformation

  !>
  !! @brief Muestra una matriz cartesianas de las particulas del sistema
  subroutine MolecularSystem_showCartesianMatrix(this,fragmentNumber,unit)
    implicit none
    type(MolecularSystem), optional, target :: this
    integer,optional :: fragmentNumber
    integer,optional :: unit
    
    type(MolecularSystem), pointer :: system

    integer :: i, j, outUnit
    real(8) :: origin(3)

    if( present(this) ) then
       system=>this
    else
       system=>MolecularSystem_instance
    end if
    
    
    outUnit=6
    if(present(unit)) outUnit=unit

    write (outUnit,"(A10,A16,A20,A20)") " ","<x>","<y>","<z>"
    
    !! Print quatum species information
    do i = 1, system%numberOfQuantumSpecies
       
       !! Avoid print twice basis in open-shell case
       if(trim(system%species(i)%name) == "E-BETA" ) cycle

       do j = 1, size(system%species(i)%particles)

          origin = system%species(i)%particles(j)%origin * AMSTRONG

          if(present(fragmentNumber) .and. (system%species(i)%particles(j)%subsystem .ne. fragmentNumber )) cycle
          
          if(system%species(i)%isElectron) then
             write (outUnit,"(A10,3F20.10)") trim( system%species(i)%particles(j)%symbol )//trim(system%species(i)%particles(j)%nickname),&
                  origin(1), origin(2), origin(3)
          else
             write (outUnit,"(A10,3F20.10)") trim(system%species(i)%particles(j)%nickname), origin(1), origin(2), origin(3)
          end if
          
       end do
    end do
    
    !! Print Point charges information
    do i = 1, system%numberOfPointCharges
       
       origin = system%pointCharges(i)%origin * AMSTRONG
       write (outUnit,"(A10,3F20.10)") trim(system%pointCharges(i)%nickname), origin(1), origin(2), origin(3)
       
    end do
    
  end subroutine MolecularSystem_showCartesianMatrix

  
  !>                      
  !! @Construye y mustra una matriz Z de coordenadas internas
  !> 
  subroutine MolecularSystem_showZMatrix( this )
    implicit none
    type(MolecularSystem) :: this

    call InternalCoordinates_constructor(  this%intCoordinates )
    call InternalCoordinates_obtainCoordinates( this%intCoordinates )
    call InternalCoordinates_destructor(  this%intCoordinates )

  end subroutine MolecularSystem_showZMatrix

  !>
  !! @brief Imprime una matriz de distancias entre particulas presentes en el sistema
  subroutine MolecularSystem_showDistanceMatrix()
    implicit none
    
    type(Matrix) :: auxMatrix
    write (6,*) ""
    write (6,"(T20,A30)") " MATRIX OF DISTANCE: ANGSTROM"
    write (6,"(T18,A35)") "------------------------------------------"
    
    auxMatrix=ParticleManager_getDistanceMatrix()
    auxMatrix%values = auxMatrix%values * AMSTRONG
    call Matrix_show(auxMatrix, rowKeys=ParticleManager_getLabelsOfCentersOfOptimization( LABELS_NUMERATED ),&
         columnKeys = ParticleManager_getLabelsOfCentersOfOptimization( LABELS_NUMERATED), flags=WITH_BOTH_KEYS  )
    
    call Matrix_destructor(auxMatrix)
    
  end subroutine MolecularSystem_showDistanceMatrix
  
  !>
  !! @brief Saves all system information.
  !! @author E. F. Posada, 2013
  subroutine MolecularSystem_saveToFile(targetFilePrefix)
    implicit none
    
    character(*), optional :: targetFilePrefix
    character(50) :: fileName, prefix    
    character(100) :: title
    integer i, j, k

    prefix="lowdin"
    if ( present( targetFilePrefix ) ) prefix=trim(targetFilePrefix)
    
    !!****************************************************************************
    !! CONTROL parameters on file.
    !!
    
    !! open file
    filename=trim(prefix)//".dat"
    open(unit=40, file=filename, status="replace", form="formatted")
    
    !!save all options
    call CONTROL_save(40)
    close(40)
    
    !!****************************************************************************
    !!Save the molecular system on file.
    !!

    !!Open file
    filename=trim(prefix)//".sys"
    open(unit=40, file=filename, status="replace", form="formatted")
    
    !! Saving general information.
    write(40,*) MolecularSystem_instance%name
    write(40,'(A100)') MolecularSystem_instance%description    
    write(40,*) MolecularSystem_instance%charge
    
    !! Saving quantum species.
    write(40,*) MolecularSystem_instance%numberOfQuantumSpecies
    !! Saving particles for each species.
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies       
       call Species_saveToFile(MolecularSystem_instance%species(i), unit=40)       
    end do    
    !! Saving Point charges
    write(40,*) MolecularSystem_instance%numberOfPointCharges

    !! Saving info of each point charge
    do i = 1, MolecularSystem_instance%numberOfPointCharges
       call Particle_saveToFile(MolecularSystem_instance%pointCharges(i), unit=40)
    end do

    !! Saving the total of particles on the system
    write(40,*) MolecularSystem_instance%numberOfParticles
    write(40,*) MolecularSystem_instance%numberOfQuantumParticles

    ! Saving External/Inter-particle potentials information
    if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
      write(40,*) ExternalPotential_instance%ssize 
      do i = 1, ExternalPotential_instance%ssize 
        write(40,*) i 
        write(40,*) ExternalPotential_instance%potentials(i)%name
        write(40,*) ExternalPotential_instance%potentials(i)%specie 
      end do

    end if

    if(CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL) then
      write(40,*) InterPotential_instance%ssize 
      do i = 1, InterPotential_instance%ssize 
        write(40,*) i 
        write(40,*) InterPotential_instance%potentials(i)%name
        write(40,*) InterPotential_instance%potentials(i)%specie 
        write(40,*) InterPotential_instance%potentials(i)%otherSpecie 
      end do

    end if

    close(40)
    
    !!****************************************************************************
    !! Saving info on the lowdin.bas format
    !!
    
    !!open file
    filename=trim(prefix)//".bas"
    open(unit=40, file=filename, status="replace", form="formatted")
    
    write(40,*) MolecularSystem_instance%numberOfQuantumSpecies    
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       write(40,*) MolecularSystem_instance%species(i)%name
       write(40,*) size(MolecularSystem_instance%species(i)%particles)
       do j = 1, size(MolecularSystem_instance%species(i)%particles)
          write(40,*) MolecularSystem_instance%species(i)%particles(j)%nickname
          write(40,*) size(MolecularSystem_instance%species(i)%particles(j)%basis%contraction)
          do k = 1, size(MolecularSystem_instance%species(i)%particles(j)%basis%contraction)
             write(40,*) &
                  MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%angularMoment, &
                  MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%length
             write(40,*) &
                  MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%origin
             write(40,*) &
                  MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%orbitalExponents
             write(40,*) &
                  MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%contractionCoefficients
          end do
       end do
    end do

    write(40,*) MolecularSystem_instance%numberOfPointCharges
    do i = 1, MolecularSystem_instance%numberOfPointCharges
       write(40,*) MolecularSystem_instance%pointCharges(i)%charge
       write(40,*) MolecularSystem_instance%pointCharges(i)%origin
    end do
    
    close(40)
    
    !!****************************************************************************
    !! Saving info for gepol program
    !!
  
    call get_command_argument (1,value=title)
    150 format (4(F10.5))
    open(unit=41, file="gepol.xyzr",status="replace",form="formatted")

      do i = 1,MolecularSystem_instance%numberOfQuantumSpecies  
        if (MolecularSystem_instance%species(i)%isElectron .eqv. .true.) then
        write(41,"(I8)") size(MolecularSystem_instance%species(i)%particles)
        do j = 1, size(MolecularSystem_instance%species(i)%particles)
            write(41,150)&
              MolecularSystem_instance%species(i)%particles(j)%origin(1)*AMSTRONG, &
              MolecularSystem_instance%species(i)%particles(j)%origin(2)*AMSTRONG, &
              MolecularSystem_instance%species(i)%particles(j)%origin(3)*AMSTRONG, &
              MolecularSystem_instance%species(i)%particles(j)%vanderwaalsRadio
        end do
        end if
      end do
    close(41)
    
    160 format (A,A)
    open(unit=42, file="gepol.inp",status="replace", form="formatted")
      write(42,160)"TITL=",trim(title)
      write(42,160)"COOF=gepol.xyzr"
      write(42,160)"VECF=vectors.vec"
      write(42,160)"LPRIN"
      write(42,160)"NDIV=5"
      ! write(42,160)"ESURF"
    
    close(42)
    
  end subroutine MolecularSystem_saveToFile

  !>
  !! @brief loads all system information from file
  !! @author E. F. Posada, 2013
  subroutine MolecularSystem_loadFromFile( form, targetPrefix )
    implicit none
    
    character(*) :: form
    character(*), optional :: targetPrefix
    character(50) :: filePrefix,fileName

    integer :: auxValue
    integer :: counter
    integer :: i, j
    logical :: existFile
    character(20) :: name
    character(50) :: species
    character(50) :: otherSpecies

    filePrefix="lowdin"
    if ( present( targetPrefix ) ) filePrefix=trim(targetPrefix)

    existFile = .false.

    select case (trim(form))
       
    case("LOWDIN.BAS")

       !!****************************************************************************
       !! Loading info from the lowdin.bas format
       !!

       !!open file
       fileName=trim(filePrefix)//".bas"
       inquire(file=fileName, exist=existFile)
       if( .not. existFile) call MolecularSystem_exception(ERROR, "The file: "//trim(fileName)//" was not found!","MolecularSystem module at LoadFromFile function.")

       !! Destroy the molecular system if any
       ! call MolecularSystem_destroy()
       if(allocated(MolecularSystem_instance%pointCharges)) deallocate(MolecularSystem_instance%pointCharges)
       if(allocated(MolecularSystem_instance%allParticles)) deallocate(MolecularSystem_instance%allParticles)

       open(unit=40, file=fileName, status="old", form="formatted")

       read(40,*) MolecularSystem_instance%numberOfQuantumSpecies
       if(.not. allocated(MolecularSystem_instance%species)) allocate(MolecularSystem_instance%species(MolecularSystem_instance%numberOfQuantumSpecies))

       MolecularSystem_instance%numberOfQuantumParticles = 0

       do i = 1, MolecularSystem_instance%numberOfQuantumSpecies

          if(allocated(MolecularSystem_instance%species(i)%particles)) deallocate(MolecularSystem_instance%species(i)%particles)

          read(40,*) MolecularSystem_instance%species(i)%name             
          read(40,*) auxValue

          allocate(MolecularSystem_instance%species(i)%particles(auxValue))

          do j = 1, size(MolecularSystem_instance%species(i)%particles)

             read(40,*) MolecularSystem_instance%species(i)%particles(j)%nickname

             MolecularSystem_instance%numberOfQuantumParticles = MolecularSystem_instance%numberOfQuantumParticles + 1
             call BasisSet_load(MolecularSystem_instance%species(i)%particles(j)%basis, filename, unit = 40)

          end do

          MolecularSystem_instance%species(i)%basisSetSize = MolecularSystem_getNumberOfContractions(i)

       end do

       read(40,*) MolecularSystem_instance%numberOfPointCharges
       allocate(MolecularSystem_instance%pointCharges(MolecularSystem_instance%numberOfPointCharges))

       do i = 1, MolecularSystem_instance%numberOfPointCharges
          read(40,*) MolecularSystem_instance%pointCharges(i)%charge
          read(40,*) MolecularSystem_instance%pointCharges(i)%origin
       end do

       close(40)

       !! Set the particles manager (all pointers)
       MolecularSystem_instance%numberOfParticles = MolecularSystem_instance%numberOfQuantumParticles + MolecularSystem_instance%numberOfPointCharges
       allocate(molecularSystem_instance%allParticles(MolecularSystem_instance%numberOfParticles ))

       counter = 1          
       do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
          do j = 1, size(MolecularSystem_instance%species(i)%particles)
             molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%species(i)%particles(j)
             counter = counter + 1
          end do
       end do

       do i = 1, MolecularSystem_instance%numberOfPointCharges
          molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%pointCharges(i)
          counter = counter + 1
       end do

       particleManager_instance => molecularSystem_instance%allParticles
       
    case("LOWDIN.DAT")
       
       !!****************************************************************************
       !! Load CONTROL parameters from file.
       !!
       
       !! open file
       fileName=trim(filePrefix)//".dat"
       inquire(file=fileName, exist=existFile)
       if( .not. existFile) call MolecularSystem_exception(ERROR, "The file: "//trim(fileName)//" was not found!","MolecularSystem module at LoadFromFile function.")
       
       open(unit=40, file=fileName, status="old", form="formatted")
       
       call CONTROL_start()
       call CONTROL_load(unit = 40)
       
       close(40)
                 
    case("LOWDIN.SYS")
       
       !!****************************************************************************
       !! Load  the molecular system from file.
       !!
       
       !! Destroy the molecular system if any
       call MolecularSystem_destroy()
    
       !! open file
       fileName=trim(filePrefix)//".sys"
       inquire(file=fileName, exist=existFile)
       if( .not. existFile) call MolecularSystem_exception(ERROR, "The file: "//trim(fileName)//" was not found!","MolecularSystem module at LoadFromFile function.")

       !!Open file
       open(unit=40, file=fileName, status="old", form="formatted")

       !! read general information.
       read(40,*) MolecularSystem_instance%name
       read(40,'(A100)') MolecularSystem_instance%description    
       read(40,*) MolecularSystem_instance%charge

       !! load quantum species.
       read(40,*) MolecularSystem_instance%numberOfQuantumSpecies
       allocate(MolecularSystem_instance%species(MolecularSystem_instance%numberOfQuantumSpecies))

       !! load particles for each species.
       do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
          call Species_loadFromFile(MolecularSystem_instance%species(i), unit=40)       
       end do

       !! load Point charges
       read(40,*) MolecularSystem_instance%numberOfPointCharges
       allocate(MolecularSystem_instance%pointCharges(MolecularSystem_instance%numberOfPointCharges))

       !! load info of each point charge
       do i = 1, MolecularSystem_instance%numberOfPointCharges
          call Particle_loadFromFile(MolecularSystem_instance%pointCharges(i), unit=40)
       end do

       !! read the total of particles on the system
       read(40,*) MolecularSystem_instance%numberOfParticles
       read(40,*) MolecularSystem_instance%numberOfQuantumParticles


       !! Set the particles manager (all pointers)              
       allocate(molecularSystem_instance%allParticles(MolecularSystem_instance%numberOfParticles ))

       counter = 1
       do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
          do j = 1, size(MolecularSystem_instance%species(i)%particles)

             molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%species(i)%particles(j)
             counter = counter + 1

          end do
       end do

       do i = 1, MolecularSystem_instance%numberOfPointCharges

          molecularSystem_instance%allParticles(counter)%particlePtr => MolecularSystem_instance%pointCharges(i)
          counter = counter + 1

       end do

       particleManager_instance => molecularSystem_instance%allParticles


       !! Loading External/Inter-particle potentials information
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

          read(40,*) auxValue
          call ExternalPotential_constructor(auxValue)

          !! FELIX TODO: create function to get potential ID

          do j = 1, ExternalPotential_instance%ssize 
             read(40,*) i 
             read(40,*) name
             read(40,*) species

             call ExternalPotential_load(i, name, species)

          end do

       end if

       if(CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL) then

          read(40,*) auxValue 
          call InterPotential_constructor(auxValue)

          do j = 1, InterPotential_instance%ssize 
             read(40,*) i 
             read(40,*) name
             read(40,*) species
             read(40,*) otherSpecies

             call InterPotential_load(i, name, species, otherSpecies)

          end do

       end if

       close(40)
          
    end select
    
  end subroutine MolecularSystem_loadFromFile

  !>
  !! @brief Returns the number of quantum species in the system.
  !! @author E. F. Posada, 2013
  function MolecularSystem_getNumberOfQuantumSpecies() result( output )
    implicit none

    integer :: output
    
    output = MolecularSystem_instance%numberOfQuantumSpecies
    
  end function MolecularSystem_getNumberOfQuantumSpecies

  !>
  !! @brief Returns the number of particles of speciesID.
  !! @author E. F. Posada, 2013
  function MolecularSystem_getNumberOfParticles(speciesID) result(output)
    implicit none

    integer :: speciesID
    integer :: output

    output = MolecularSystem_instance%species(speciesID)%internalSize
    
  end function MolecularSystem_getNumberOfParticles
  
  !>
  !! @brief Returns the number of shells for specie.
  !! @author E. F. Posada, 2013
  function MolecularSystem_getNumberOfContractions( specieID ) result( output )
    implicit none
    integer :: specieID
    integer :: output
    
    integer :: i, j

    output = 0

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies

       if ( specieID ==  i ) then

          do j = 1, size(MolecularSystem_instance%species(i)%particles)
             
             output = output + size(MolecularSystem_instance%species(i)%particles(j)%basis%contraction)
             
          end do
          
       end if
       
    end do
    
  end function MolecularSystem_getNumberOfContractions
  
  !>
  !! @brief Returns the number of cartesian shells for specie.
  !! @author E. F. Posada, 2013
  function MolecularSystem_getTotalNumberOfContractions( specieID, this ) result( output )
    implicit none
    integer :: specieID
    type(MolecularSystem), optional, target :: this

    type(MolecularSystem), pointer :: system
    integer :: output
    integer :: j, k

    output = 0

    if( present(this) ) then
       system=>this
    else
       system=>MolecularSystem_instance
    end if
    
    do j = 1, size(system%species(specieID)%particles)

       do k = 1, size(system%species(specieID)%particles(j)%basis%contraction)
          
          output = output + system%species(specieID)%particles(j)%basis%contraction(k)%numCartesianOrbital
          
       end do
       
    end do
          
  end function MolecularSystem_getTotalNumberOfContractions

  !> @brief gets all basis set for one specie as an array of contractions
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine MolecularSystem_getBasisSet(specieID, output)
    implicit none
    
    integer :: specieID
    type(ContractedGaussian), allocatable, intent(inout) :: output(:)
    integer :: i, j
    integer :: counter
    
    if (allocated(output)) deallocate(output)
    allocate(output(MolecularSystem_instance%species(specieID)%basisSetSize))
    
    counter = 0
    
    do i = 1, size(MolecularSystem_instance%species(specieID)%particles)
       do j = 1, size(MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction)
          
          counter = counter + 1
          output(counter)= MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction(j)
          
       end do
    end do

   end subroutine  MolecularSystem_getBasisSet
   
   !> @brief find de maximun angular moment for specie specieID
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getMaxAngularMoment(specieID) result(output)
     implicit none
     
     integer :: specieID
     integer :: output
     
     integer :: i, j
     
     output = -1
     
     do i = 1, size(MolecularSystem_instance%species(specieID)%particles)
        do j = 1, size(MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction)
           
           output = max(output, MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction(j)%angularMoment)
           
        end do
     end do
     
   end function MolecularSystem_getMaxAngularMoment

   !> @brief find the maximun number of primitives for specieID, necessary for derive with libint
   !! @author J.M. Rodas 2015
   !! @version 1.0
   function MolecularSystem_getMaxNumberofPrimitives(specieID) result(output)
     implicit none
     
     integer :: specieID
     integer :: output
     
     integer :: i, j
     
     output = -1
     
     do i = 1, size(MolecularSystem_instance%species(specieID)%particles)
        do j = 1, size(MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction)
           
           output = max(output, MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction(j)%length)
           
        end do
     end do

   end function MolecularSystem_getMaxNumberofPrimitives

   !> @brief find de maximun number of primitives for specieID, necessary for derive with libint
   !! @author J.M. Rodas 2015
   !! @version 1.0
   function MolecularSystem_getMaxNumberofCartesians(specieID) result(output)
     implicit none
     
     integer :: specieID
     integer :: output
     
     integer :: i, j
     
     output = -1
     
     do i = 1, size(MolecularSystem_instance%species(specieID)%particles)
        do j = 1, size(MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction)
           
           output = max(output, MolecularSystem_instance%species(specieID)%particles(i)%basis%contraction(j)%numCartesianOrbital)
           
        end do
     end do

   end function MolecularSystem_getMaxNumberofCartesians
     
   !> @brief Returns the occupation number of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getOcupationNumber(speciesID,this) result(output)
     implicit none
     integer :: speciesID
     type(MolecularSystem), optional, target :: this
     integer :: output

     type(MolecularSystem), pointer :: system

     if( present(this) ) then
        system=>this
     else
        system=>MolecularSystem_instance
     end if

     output = -1
     output = system%species(speciesID)%ocupationNumber
          
   end function MolecularSystem_getOcupationNumber

   !> @brief Returns the eta parameter of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getEta(speciesID,this) result(output)
     implicit none
     
     integer :: speciesID
     type(MolecularSystem), optional, target :: this
     integer :: output

     type(MolecularSystem), pointer :: system

     if( present(this) ) then
        system=>this
     else
        system=>MolecularSystem_instance
     end if
     
     output = -1
     output = system%species(speciesID)%eta
          
   end function MolecularSystem_getEta

   function MolecularSystem_getLambda(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     integer :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%lambda
          
   end function MolecularSystem_getLambda


   function MolecularSystem_getKappa(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     integer :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%kappa
          
   end function MolecularSystem_getKappa

   function MolecularSystem_getMultiplicity(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     integer :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%spin
          
   end function MolecularSystem_getMultiplicity




   function MolecularSystem_getParticlesFraction(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     real(8) :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%particlesFraction
          
   end function MolecularSystem_getParticlesFraction


   !> @brief Returns the charge of speciesID
   !! @author E. F. Posada, 2013
   !! @version 1.0   
   function MolecularSystem_getCharge( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%charge
     
   end function MolecularSystem_getCharge

   !> @brief Returns the omega frequency of speciesID. Why we have these functions??
   function MolecularSystem_getOmega( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%omega
     
   end function MolecularSystem_getOmega

   !> @brief Returns the mass of speciesID
   !! @author E. F. Posada, 2013
   !! @version 1.0   
   function MolecularSystem_getMass( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%mass
     
   end function MolecularSystem_getMass

   !> @brief Returns QDO center of quantum species
   function MolecularSystem_getQDOcenter( speciesID ) result( origin )
     implicit none
     integer :: speciesID
     integer :: i
     logical :: centerFound
     real(8) :: origin(3)

     centerFound = .False.
     do i = 1 , size( MolecularSystem_instance%pointCharges )
        if ( trim(MolecularSystem_instance%pointCharges(i)%qdoCenterOf) == trim(MolecularSystem_instance%species(speciesID)%symbol) ) then 
          origin = MolecularSystem_instance%pointCharges(i)%origin 
          centerFound = .True.
          exit
        end if
    end do
    if ( .not. centerFound ) then
        call MolecularSystem_exception(ERROR, "No QDO center for species: "//MolecularSystem_instance%species(speciesID)%symbol, "MolecularSystem_getQDOcenter"   )
    end if
 
     
   end function MolecularSystem_getQDOCenter

   
   !> @brief Returns the Factor Of Exchange Integrals
   !! @author E. F. Posada, 2013
   !! @version 1.0   
   function MolecularSystem_getFactorOfExchangeIntegrals( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%kappa / MolecularSystem_instance%species(speciesID)%eta
     
   end function MolecularSystem_getFactorOfExchangeIntegrals

   !> @brief Returns the name of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getNameOfSpecie(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     character(30) :: output
     
     output = MolecularSystem_instance%species(speciesID)%name
          
   end function MolecularSystem_getNameOfSpecie

   !> @brief Returns the name of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getNameOfSpecies(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     character(30) :: output
     
     output = MolecularSystem_instance%species(speciesID)%name
          
   end function MolecularSystem_getNameOfSpecies
   
   !> @brief Returns the name of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getSpecieID( nameOfSpecie ) result(output)
     implicit none
     
     character(*) :: nameOfSpecie
     integer :: output
     integer i 
     
     output = 0

     do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
        if( trim(MolecularSystem_instance%species(i)%name) == trim(nameOfSpecie)) output = i
     end do

   end function MolecularSystem_getSpecieID

      !> @brief Returns the name of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getSpecieIDFromSymbol( symbolOfSpecie ) result(output)
     implicit none
     
     character(*) :: symbolOfSpecie
     integer :: output
     integer i 
     
     output = 0

     do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
        if( trim(MolecularSystem_instance%species(i)%symbol) == trim(symbolOfSpecie)) output = i
     end do

   end function MolecularSystem_getSpecieIDFromSymbol

   !>
   !! @brief calcula la energia total para una especie especificada
   function MolecularSystem_getPointChargesEnergy() result( output )
     implicit none
     real(8) :: output
     
     integer :: i
     integer :: j
     real(8) :: deltaOrigin(3)
     
     output =0.0_8
     
     do i=1, size( MolecularSystem_instance%pointCharges )      
        do j = i + 1 , size( MolecularSystem_instance%pointCharges )
            
           deltaOrigin = MolecularSystem_instance%pointCharges(i)%origin &
                - MolecularSystem_instance%pointCharges(j)%origin
           
           output=output + ( ( MolecularSystem_instance%pointCharges(i)%charge &
                * MolecularSystem_instance%pointCharges(j)%charge )&
                / sqrt( sum( deltaOrigin**2.0_8 ) ) )
           
        end do
     end do
    
    !! Point charge potential with the external electric field
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      do i=1, size( MolecularSystem_instance%pointCharges )      
        output = output + sum(CONTROL_instance%ELECTRIC_FIELD(:) * MolecularSystem_instance%pointCharges(i)%origin(:) )*  MolecularSystem_instance%pointCharges(i)%charge 
      end do
    end if

     
   end function MolecularSystem_getPointChargesEnergy

   function MolecularSystem_getMMPointChargesEnergy() result( output )
     implicit none
     real(8) :: output

     integer :: i
     integer :: j
     real(8) :: deltaOrigin(3)

     output =0.0_8
     
     do i=1, size( MolecularSystem_instance%pointCharges )
        if(trim(MolecularSystem_instance%pointCharges(i)%nickname) == "PC") then
           do j = i + 1 , size( MolecularSystem_instance%pointCharges )

              deltaOrigin = MolecularSystem_instance%pointCharges(i)%origin &
                   - MolecularSystem_instance%pointCharges(j)%origin

              output=output + ( ( MolecularSystem_instance%pointCharges(i)%charge &
                   * MolecularSystem_instance%pointCharges(j)%charge )&
                   / sqrt( sum( deltaOrigin**2.0_8 ) ) )

           end do
        end if
     end do

   end function MolecularSystem_getMMPointChargesEnergy
   
   !>
   !! @brief returns an array of labels of all basis set of speciesID
   function MolecularSystem_getlabelsofcontractions(speciesID) result(output)
     implicit none
     
     integer :: speciesID     
     character(19),allocatable :: output(:)

     integer :: i, j, k
     integer :: counter
     character(9), allocatable :: shellCode(:)

     if(allocated(output)) deallocate(output)
     allocate(output(MolecularSystem_getTotalNumberOfContractions(speciesID)))
     
     output = ""
     counter = 1
     
     do i = 1, size(MolecularSystem_instance%species(speciesID)%particles)
        do j = 1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
           
           if(allocated(shellCode)) deallocate(shellCode)
           allocate(shellCode(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital))
           shellCode = ""

           shellCode = ContractedGaussian_getShellCode(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j))
           
           do k = 1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
              
              write (output(counter),"(I5,A1,A6,A1,A6)") counter, " ", &
                   trim(MolecularSystem_instance%species(speciesID)%particles(i)%nickname), " ", &
                   trim(shellCode(k))//" "

              counter = counter + 1 
              
           end do

        end do
     end do
     
   end function MolecularSystem_getlabelsofcontractions

   !>
   !! @brief  Change from Lowdin order to Molden/Gaussian or Gamess order
   !! 
   subroutine MolecularSystem_changeOrbitalOrder( coefficientsOfCombination, speciesID, actualFormat, desiredFormat )
     implicit none
     type(Matrix), intent(inout) :: coefficientsOfCombination
     integer, intent(in) :: speciesID
     character(*), intent(in) :: actualFormat
     character(*), intent(in) :: desiredFormat
     character(19) , allocatable :: labelsOfContractions(:)
     integer :: numberOfContractions
     character(6) :: nickname
     character(6) :: shellCode
     character(1) :: space
     integer :: k, counter, auxcounter

     numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID)
     !! Build a vector of labels of contractions
     if(allocated(labelsOfContractions)) deallocate(labelsOfContractions)
     allocate(labelsOfContractions(numberOfContractions))

     labelsOfContractions =  MolecularSystem_getlabelsofcontractions(speciesID)

     if( (actualFormat.eq."LOWDIN" .and. desiredFormat.eq."MOLDEN") ) then
        !! Swap some columns according to the molden format
        do k=1,numberOfContractions
           !! Take the shellcode
           read (labelsOfContractions(k), "(I5,A1,A6,A1,A6)") counter, space, nickname, space, shellcode 

           !! Reorder the D functions
           !! counter:  0,  1,  2,  3,  4,  5
           !! Lowdin:  XX, XY, XZ, YY, YZ, ZZ
           !! Molden:  XX, YY, ZZ, XY, XZ, YZ 
           !!  1-1, 2-4, 3-5, 4-2, 5-6, 6-3
           !!  2-4, 3-5, 5-6

           if ( adjustl(shellcode) == "Dxx" ) then 
              auxcounter = counter
              !! Swap XY and YY
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+3)
              !! Swap XZ and ZZ
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+5)
              !! Swap YZ and XZ'
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
           end if

           !! Reorder the F functions
           !! counter:   0,   1,   2,   3,   4,   5,   6,   7,   8    9
           !! Lowdin:  XXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ
           !! Molden:  XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ

           if ( adjustl(shellcode) == "Fxxx" ) then 
              auxcounter = counter
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)

           end if

           !! Reorder the G functions
           !! counter:   0,  1,   2,   3,   4,   5,   6,   7,   8    9,   10,  11,  12,  13,  14
           !! Lowdin:  XXXX,XXXY,XXXZ,XXYY,XXYZ,XXZZ,XYYY,XYYZ,XYZZ,XZZZ,YYYY,YYYZ,YYZZ,YZZZ,ZZZZ
           !! erkale-FCHK:  ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

           ! Molden 15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy
           if ( adjustl(shellcode) == "Gxxxx" ) then
              auxcounter = counter
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter   , auxcounter+14)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+13)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+12)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+3 , auxcounter+11)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+10)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+8)

           end if

           if ( adjustl(shellcode) == "Hxxxxx" ) then 
              call MolecularSystem_exception(WARNING, "The order of the coefficients only works until G orbitals", "MolecularSystem_changeOrbitalOrder" )
           end if

        end do

     else if( ( actualFormat.eq."LOWDIN" .and. desiredFormat.eq."GAMESS") ) then
        !! Swap some columns according to the Gamess format
        do k=1,numberOfContractions
           !! Take the shellcode
           read (labelsOfContractions(k), "(I5,A1,A6,A1,A6)") counter, space, nickname, space, shellcode 

           !! Reorder the D functions
           !! counter:  1,  2,  3,  4,  5,  6
           !! Lowdin:  XX, XY, XZ, YY, YZ, ZZ
           !! Molden:  XX, YY, ZZ, XY, XZ, YZ 
           !!  1-1, 2-4, 3-5, 4-2, 5-6, 6-3
           !!  2-4, 3-5, 5-6

           if ( adjustl(shellcode) == "Dxx" ) then 
              auxcounter = counter
              !! Swap XY and YY
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+3)
              !! Swap XZ and ZZ
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+5)
              !! Swap YZ and XZ'
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
           end if

           !! Reorder the F functions
           !! counter:   1,   2,   3,   4,   5,   6,   7,   8    9,  10
           !! Lowdin:  XXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ
           !! Molden:  XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
           !! Gamess:  XXX, YYY, ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ

           if ( adjustl(shellcode) == "Fxxx" ) then 
              auxcounter = counter
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+3 , auxcounter+4)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+8 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)
           end if

           if ( adjustl(shellcode) == "Gxxxx" ) then 
              call MolecularSystem_exception(WARNING, "The order of the coefficients only works until F orbitals", "MolecularSystem_changeOrbitalOrder" )
           end if

           
        end do

     else if( actualFormat.eq."MOLDEN" .and. desiredFormat.eq."LOWDIN") then
        !! Swap some columns according to the molden format
        do k=1,numberOfContractions
           !! Take the shellcode
           read (labelsOfContractions(k), "(I5,A1,A6,A1,A6)") counter, space, nickname, space, shellcode 

           !! Reorder the D functions
           !! counter:  1,  2,  3,  4,  5,  6
           !! Molden:  XX, YY, ZZ, XY, XZ, YZ 
           !! Lowdin:  XX, XY, XZ, YY, ZZ, YZ
           !!  1-1, 2-4, 3-5, 4-2, 5-6, 6-3
           !!  2-4, 3-5, 5-6

           if ( adjustl(shellcode) == "Dxx" ) then 
              auxcounter = counter
              !! Swap YY and XY
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+3)
              !! Swap ZZ and XZ
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+4)
              !! Swap ZZ and YZ'
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
           end if

           !! Reorder the F functions
           !! counter:   1,   2,   3,   4,   5,   6,   7,   8    9,  10
           !! Molden:  XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
           !! Lowdin:  XXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ

           if ( adjustl(shellcode) == "Fxxx" ) then 
              auxcounter = counter
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+4)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+5)             
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)
           end if
           !! Reorder the G functions
           !! counter:   0,  1,   2,   3,   4,   5,   6,   7,   8    9,   10,  11,  12,  13,  14
           !! erkale-FCHK:  ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
           !! Lowdin:  XXXX,XXXY,XXXZ,XXYY,XXYZ,XXZZ,XYYY,XYYZ,XYZZ,XZZZ,YYYY,YYYZ,YYZZ,YZZZ,ZZZZ
           if ( adjustl(shellcode) == "Gxxxx" ) then 
              auxcounter = counter
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter   , auxcounter+14)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+13)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+12)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+3 , auxcounter+11)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+10)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+8)

           end if

           if ( adjustl(shellcode) == "Hxxxxx" ) then 
              call MolecularSystem_exception(WARNING, "The order of the coefficients only works until G orbitals", "MolecularSystem_changeOrbitalOrder" )
           end if

           
        end do

     else

        call MolecularSystem_exception(ERROR, "The desired format change from "//actualFormat//" to "//desiredFormat//"has not been implemented","MolecularSystem module at changeOrbitalOrder function.")

     end if

   end subroutine MolecularSystem_changeOrbitalOrder

   

  !>
  !! @brief Lee la matriz de densidad y los orbitales de un archivo fchk tipo Gaussian
  subroutine MolecularSystem_readFchk( fileName, coefficients, densityMatrix, nameOfSpecies )
    implicit none

    character(*), intent(in) :: fileName
    type(Matrix), intent(inout) :: coefficients
    type(Matrix), intent(inout) :: densityMatrix
    character(*) :: nameOfSpecies

    integer :: speciesID
    integer :: numberOfContractions
    integer :: i, j, k
    character(100) :: info
    character(40) :: object
    character(4) :: type
    character(6) :: array
    integer :: size
    integer :: fchkUnit, io
    integer :: auxInteger(6)
    integer, allocatable :: integerArray(:)
    real(8) :: auxReal(5)
    real(8), allocatable :: realArray(:)
    logical :: existFchk
    

    speciesID=MolecularSystem_getSpecieID(nameOfSpecies)
    numberOfContractions=MolecularSystem_getTotalnumberOfContractions( speciesID )
    inquire(FILE = trim(fileName), EXIST = existFchk )
    if ( .not. existFchk ) call MolecularSystem_exception( ERROR, "I did not find any .fchk coefficients file", "At MolecularSystem_readFchk")
    

    fchkUnit = 50

    open(unit=fchkUnit, file=filename, status="old", form="formatted", access='sequential', action='read')

    print *, ""
    print *, "reading FCHK orbitals from ", fileName
    !The first two lines don't matter
    read(fchkUnit,"(A100)",iostat=io) info
    print *, info
    read(fchkUnit,"(A100)",iostat=io) info
    print *, info

    !Read line by line, the first 40 characters determine what's being read
    do 
       read(fchkUnit,"(A40,A4,A6,I14)",iostat=io) object,type,array,size
       if (io.ne.0) exit
       ! print *, object,type,array,size
       if( trim(adjustl(array)).eq."N=" .and. trim(adjustl(type)).eq."I") then
          ! print *, "voy a leer el arreglo ", object , " de enteros de ", size, "en lineas ", ceiling( size/6.0 )
          allocate(integerArray(size))
          k=0
          do i=1, ceiling( size/6.0 )
             read(fchkUnit,"(6I12)") auxInteger(1:6)
             do j=1,6
                k=k+1
                if(k.gt.size)exit
                integerArray(k)=auxInteger(j)
             end do
          end do
          ! print *, integerArray   
          deallocate(integerArray)

       else if( trim(adjustl(array)).eq."N=" .and. trim(adjustl(type)).eq."R") then
          ! print *, "voy a leer el arreglo ", object , " de reales", size, "en lineas ", ceiling( size/5.0 )
          allocate(realArray(size))
          k=0
          do i=1, ceiling( size/5.0 )
             read(fchkUnit,"(5E16.8)") auxReal(1:5)
             do j=1,5
                k=k+1
                if(k.gt.size)exit
                realArray(k)=auxReal(j)
             end do
          end do
          ! if(trim(adjustl(object)).eq."Total SCF Density") then
          !    k=0
          !    do i=1, numberOfContractions
          !       do j=1,i
          !          k=k+1
          !          densityMatrix%values(i,j)=realArray(k)
          !          densityMatrix%values(j,i)=realArray(k)
          !       end do
          !    end do
          !    print *, "density matrix read"
          !    call Matrix_show(densityMatrix)
          ! end if
          if(trim(adjustl(object)).eq."Alpha MO coefficients") then
             k=0
             do i=1, numberOfContractions
                do j=1, numberOfContractions
                   k=k+1
                   coefficients%values(j,i)=realArray(k)
                end do
             end do
          end if
          ! print *, realArray   
          deallocate(realArray)
       end if
    end do

    call MolecularSystem_changeOrbitalOrder( coefficients, speciesID, "MOLDEN", "LOWDIN" )
    ! print *, "coefficients read"
    ! call Matrix_show(coefficients)

    !Build density matrix with the new order
    densityMatrix%values=0.0
    do i=1, numberOfContractions
       do j=1, numberOfContractions
          do k=1, MolecularSystem_getOcupationNumber(speciesID)
             densityMatrix%values(i,j) =  &
                  densityMatrix%values( i,j ) + &
                  coefficients%values(i,k)*coefficients%values(j,k)*MolecularSystem_getEta(speciesID)
          end do
       end do
    end do
    ! print *, "density matrix from orbitals read"
    ! call Matrix_show(densityMatrix)

    
    close(fchkUnit)
    
  end subroutine MolecularSystem_readFchk

   
  !>
  !! @brief Copies a molecular system into another molecular system object
  !! @author F. M. Moncada 2022
  subroutine MolecularSystem_copyConstructor(this,originalThis)
    implicit none

    type(MolecularSystem), intent(out), target :: this
    type(MolecularSystem), intent(in) :: originalThis

    integer :: i, j
    integer :: counter

    !! Destroy the molecular system if any
    ! call MolecularSystem_destroy()
    if(allocated(this%pointCharges)) deallocate(this%pointCharges)
    if(allocated(this%allParticles)) deallocate(this%allParticles)

    !! Start copy atribute by atribute
    !! Non structured information 
    this%name=originalThis%name 
    this%description=originalThis%description
    this%numberOfParticles=originalThis%numberOfParticles
    this%numberOfPointCharges=originalThis%numberOfPointCharges 
    this%numberOfQuantumParticles=originalThis%numberOfQuantumParticles
    this%numberOfQuantumSpecies=originalThis%numberOfQuantumSpecies
    this%charge=originalThis%charge

    !! Allocate memory for species and particles
    allocate(this%species(originalThis%numberOfQuantumSpecies))

    !!Copies species information
    do i = 1, originalThis%numberOfQuantumSpecies
       allocate(this%species(i)%particles( size(originalThis%species(i)%particles)))
       this%species(i)%name = originalThis%species(i)%name
       this%species(i)%symbol = originalThis%species(i)%symbol
       this%species(i)%statistics = originalThis%species(i)%statistics
       this%species(i)%charge = originalThis%species(i)%charge
       this%species(i)%mass = originalThis%species(i)%mass
       this%species(i)%omega = originalThis%species(i)%omega
       this%species(i)%spin = originalThis%species(i)%spin
       this%species(i)%totalCharge = originalThis%species(i)%totalCharge
       this%species(i)%kappa = originalThis%species(i)%kappa
       this%species(i)%eta = originalThis%species(i)%eta
       this%species(i)%lambda = originalThis%species(i)%lambda
       this%species(i)%particlesFraction = originalThis%species(i)%particlesFraction
       this%species(i)%ocupationNumber = originalThis%species(i)%ocupationNumber
       this%species(i)%multiplicity = originalThis%species(i)%multiplicity
       this%species(i)%internalSize = originalThis%species(i)%internalSize
       this%species(i)%basisSetSize = originalThis%species(i)%basisSetSize
       this%species(i)%speciesID = originalThis%species(i)%speciesID
       this%species(i)%isElectron = originalThis%species(i)%isElectron
    end do

    allocate(this%pointCharges(originalThis%numberOfPointCharges))
    allocate(this%allParticles(originalThis%numberOfParticles))
       
    !! Set the particles manager (all pointers)        
    counter = 1
    do i = 1, originalThis%numberOfQuantumSpecies
       do j = 1, size(originalThis%species(i)%particles)        
          this%allParticles(counter)%particlePtr => this%species(i)%particles(j)
          counter = counter + 1
       end do
    end do

    do i = 1, originalThis%numberOfPointCharges       
       this%allParticles(counter)%particlePtr => this%pointCharges(i)
       counter = counter + 1
    end do
    
    
    !!Copies particles information
    do i = 1, originalThis%numberOfParticles
       this%allParticles(i)%particlePtr%name=originalThis%allParticles(i)%particlePtr%name
       this%allParticles(i)%particlePtr%symbol=originalThis%allParticles(i)%particlePtr%symbol
       this%allParticles(i)%particlePtr%nickname=originalThis%allParticles(i)%particlePtr%nickname
       this%allParticles(i)%particlePtr%statistics=originalThis%allParticles(i)%particlePtr%statistics
       this%allParticles(i)%particlePtr%basisSetName=originalThis%allParticles(i)%particlePtr%basisSetName
       this%allParticles(i)%particlePtr%origin=originalThis%allParticles(i)%particlePtr%origin
       this%allParticles(i)%particlePtr%charge=originalThis%allParticles(i)%particlePtr%charge
       this%allParticles(i)%particlePtr%mass=originalThis%allParticles(i)%particlePtr%mass
       this%allParticles(i)%particlePtr%omega=originalThis%allParticles(i)%particlePtr%omega
       this%allParticles(i)%particlePtr%qdoCenterOf=originalThis%allParticles(i)%particlePtr%qdoCenterOf
       this%allParticles(i)%particlePtr%spin=originalThis%allParticles(i)%particlePtr%spin
       this%allParticles(i)%particlePtr%totalCharge=originalThis%allParticles(i)%particlePtr%totalCharge
       this%allParticles(i)%particlePtr%klamt=originalThis%allParticles(i)%particlePtr%klamt
       this%allParticles(i)%particlePtr%vanderWaalsRadio=originalThis%allParticles(i)%particlePtr%vanderWaalsRadio
       this%allParticles(i)%particlePtr%isQuantum=originalThis%allParticles(i)%particlePtr%isQuantum
       this%allParticles(i)%particlePtr%isDummy=originalThis%allParticles(i)%particlePtr%isDummy
       this%allParticles(i)%particlePtr%fixComponent=originalThis%allParticles(i)%particlePtr%fixComponent
       this%allParticles(i)%particlePtr%isCenterOfOptimization=originalThis%allParticles(i)%particlePtr%isCenterOfOptimization
       this%allParticles(i)%particlePtr%multiplicity=originalThis%allParticles(i)%particlePtr%multiplicity
       this%allParticles(i)%particlePtr%subsystem=originalThis%allParticles(i)%particlePtr%subsystem
       this%allParticles(i)%particlePtr%translationCenter=originalThis%allParticles(i)%particlePtr%translationCenter
       this%allParticles(i)%particlePtr%rotationPoint=originalThis%allParticles(i)%particlePtr%rotationPoint
       this%allParticles(i)%particlePtr%rotateAround=originalThis%allParticles(i)%particlePtr%rotateAround
       this%allParticles(i)%particlePtr%id=originalThis%allParticles(i)%particlePtr%id
       this%allParticles(i)%particlePtr%internalSize=originalThis%allParticles(i)%particlePtr%internalSize
       this%allParticles(i)%particlePtr%owner=originalThis%allParticles(i)%particlePtr%owner
       this%allParticles(i)%particlePtr%basisSetSize=originalThis%allParticles(i)%particlePtr%basisSetSize

       if ( allocated(originalThis%allParticles(i)%particlePtr%childs) ) then
          allocate(this%allParticles(i)%particlePtr%childs( size(originalThis%allParticles(i)%particlePtr%childs)))
          this%allParticles(i)%particlePtr%childs=originalThis%allParticles(i)%particlePtr%childs
       end if

       !! Copies basis set information
       if ( this%allParticles(i)%particlePtr%isQuantum ) then

          this%allParticles(i)%particlePtr%basis%name=originalThis%allParticles(i)%particlePtr%basis%name
          this%allParticles(i)%particlePtr%basis%origin=originalThis%allParticles(i)%particlePtr%basis%origin
          this%allParticles(i)%particlePtr%basis%length=originalThis%allParticles(i)%particlePtr%basis%length
          this%allParticles(i)%particlePtr%basis%ttype=originalThis%allParticles(i)%particlePtr%basis%ttype
          this%allParticles(i)%particlePtr%basis%contractionLength=originalThis%allParticles(i)%particlePtr%basis%contractionLength
          this%allParticles(i)%particlePtr%basis%numberOfPrimitives=originalThis%allParticles(i)%particlePtr%basis%numberOfPrimitives

          allocate(this%allParticles(i)%particlePtr%basis%contraction(originalThis%allParticles(i)%particlePtr%basis%length))
          
          do j=1, originalThis%allParticles(i)%particlePtr%basis%length
             this%allParticles(i)%particlePtr%basis%contraction(j)%id=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%id
             this%allParticles(i)%particlePtr%basis%contraction(j)%length=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%length
             this%allParticles(i)%particlePtr%basis%contraction(j)%angularMoment=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%angularMoment
             this%allParticles(i)%particlePtr%basis%contraction(j)%numCartesianOrbital=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%numCartesianOrbital
             this%allParticles(i)%particlePtr%basis%contraction(j)%owner=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%owner
             this%allParticles(i)%particlePtr%basis%contraction(j)%subsystem=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%subsystem
             this%allParticles(i)%particlePtr%basis%contraction(j)%origin=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%origin

             allocate(this%allParticles(i)%particlePtr%basis%contraction(j)%orbitalExponents(&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%length))
             allocate(this%allParticles(i)%particlePtr%basis%contraction(j)%contractionCoefficients(&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%length))

             this%allParticles(i)%particlePtr%basis%contraction(j)%orbitalExponents=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%orbitalExponents
             this%allParticles(i)%particlePtr%basis%contraction(j)%contractionCoefficients=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%contractionCoefficients

             allocate(this%allParticles(i)%particlePtr%basis%contraction(j)%contNormalization(&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%numCartesianOrbital))
             allocate(this%allParticles(i)%particlePtr%basis%contraction(j)%primNormalization(&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%length,&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%numCartesianOrbital))

             this%allParticles(i)%particlePtr%basis%contraction(j)%contNormalization=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%contNormalization
             this%allParticles(i)%particlePtr%basis%contraction(j)%primNormalization=&
                  originalThis%allParticles(i)%particlePtr%basis%contraction(j)%primNormalization                  

          end do
       end if

    end do
    
    ! particleManager_instance => this%allParticles
    
  end subroutine MolecularSystem_copyConstructor

  !>
  !! @brief Stack together the basis sets of two molecular systems for Non Orthogonal CI calculations
  !! Adds system B particles to system A, sysBasisList indicate the position shifts of the basis functions in the merged molecular system
  !! @author F. M. Moncada 2022
  subroutine MolecularSystem_mergeTwoSystems(mergedThis,thisA,thisB,sysAbasisList, sysBbasisList, reorder)
    type(MolecularSystem), intent(out), target :: mergedThis
    type(MolecularSystem), intent(in) :: thisA, thisB
    type(IVector) :: sysAbasisList(*), sysBbasisList(*) !length = numberOfSpecies
    logical, optional :: reorder !reorder system A to put common basis functions first
    
    integer :: i, j, k, l, jj, speciesID, mu, nu
    integer :: counter, common
    integer, allocatable :: notCommonParticles(:), notCommonBasisSize(:), auxMuPositions(:), auxNuPositions(:)
    type(IVector1), allocatable :: notCommonListA(:), notCommonListB(:) !1=not common, 0=common
    logical :: reorderA
    
    !! Destroy the molecular system if any
    ! call MolecularSystem_destroy()
    if(allocated(mergedThis%pointCharges)) deallocate(mergedThis%pointCharges)
    if(allocated(mergedThis%allParticles)) deallocate(mergedThis%allParticles)

    reorderA=.true.
    if(present(reorder)) reorderA=reorder

    !! Start copy atribute by atribute
    !! Non structured information 
    mergedThis%name=thisA%name 
    mergedThis%description=thisA%description
    mergedThis%charge=thisA%charge
    mergedThis%numberOfPointCharges=thisA%numberOfPointCharges 
    mergedThis%numberOfQuantumParticles=thisA%numberOfQuantumParticles
    mergedThis%numberOfQuantumSpecies=thisA%numberOfQuantumSpecies
    mergedThis%numberOfParticles=thisA%numberOfPointCharges 

    !! Allocate memory for particles
    allocate(mergedThis%species(mergedThis%numberOfQuantumSpecies))

    !! Count the number of basis functions of A that are also in B
    !! a comparison of origin, basis set names etc    
    allocate(notCommonParticles(mergedThis%numberOfQuantumSpecies),&
         notCommonBasisSize(mergedThis%numberOfQuantumSpecies),&
         notCommonListA(mergedThis%numberOfQuantumSpecies),&
         notCommonListB(mergedThis%numberOfQuantumSpecies))

    do i = 1, mergedThis%numberOfQuantumSpecies
       notCommonParticles(i)=0
       notCommonBasisSize(i)=0
       call Vector_constructorInteger1(notCommonListA(i), int(size(thisA%species(i)%particles),8) , 1_1) 
       call Vector_constructorInteger1(notCommonListB(i), int(size(thisB%species(i)%particles),8) , 0_1) 
       Bloop: do j = 1, size(thisB%species(i)%particles)        
          Aloop: do k= 1, size(thisA%species(i)%particles)
          ! if( thisA%species(i)%particles(j)%translationCenter .ne. 0 .or. thisA%species(i)%particles(j)%rotateAround .ne. 0) then
          !    notCommonParticles(i)=notCommonParticles(i)+1
          !    notCommonBasisSize(i)=notCommonBasisSize(i)+thisA%species(i)%particles(j)%basisSetSize
          ! end if
             if ( MolecularSystem_checkParticleEquivalence(thisA%species(i)%particles(k),thisB%species(i)%particles(j)) .eqv. .true.) then
                notCommonListA(i)%values(k)=0
                cycle Bloop
             end if
          end do Aloop
          notCommonParticles(i)=notCommonParticles(i)+1
          notCommonBasisSize(i)=notCommonBasisSize(i)+thisB%species(i)%particles(j)%basisSetSize
          notCommonListB(i)%values(j)=1 
       end do Bloop
       ! print *, "species", i, "has", notCommonParticles(i), "notCommonParticles and", notCommonBasisSize(i), "notCommonBasisSize"
       ! print *, notCommonListB(i)%values(:)
       if(.not. reorderA) notCommonListA(i)%values=0
    end do

    
    !!Copies species information
    !!A basis set plus B basis functions not included in A
    do i = 1, mergedThis%numberOfQuantumSpecies
       allocate(mergedThis%species(i)%particles( size(thisA%species(i)%particles) + notCommonParticles(i) ))
       mergedThis%species(i)%name = thisA%species(i)%name
       mergedThis%species(i)%symbol = thisA%species(i)%symbol
       mergedThis%species(i)%statistics = thisA%species(i)%statistics
       mergedThis%species(i)%charge = thisA%species(i)%charge
       mergedThis%species(i)%mass = thisA%species(i)%mass
       mergedThis%species(i)%omega = thisA%species(i)%omega
       mergedThis%species(i)%spin = thisA%species(i)%spin
       mergedThis%species(i)%totalCharge = thisA%species(i)%totalCharge
       mergedThis%species(i)%kappa = thisA%species(i)%kappa
       mergedThis%species(i)%eta = thisA%species(i)%eta
       mergedThis%species(i)%lambda = thisA%species(i)%lambda
       mergedThis%species(i)%particlesFraction = thisA%species(i)%particlesFraction
       mergedThis%species(i)%ocupationNumber = thisA%species(i)%ocupationNumber+thisB%species(i)%ocupationNumber
       mergedThis%species(i)%multiplicity = thisA%species(i)%multiplicity
       mergedThis%species(i)%internalSize = thisA%species(i)%internalSize
       mergedThis%species(i)%basisSetSize = thisA%species(i)%basisSetSize+notCommonBasisSize(i)
       mergedThis%species(i)%speciesID = thisA%species(i)%speciesID
       mergedThis%species(i)%isElectron = thisA%species(i)%isElectron
       mergedThis%numberOfParticles=mergedThis%numberOfParticles+size(mergedThis%species(i)%particles)
    end do

    allocate(mergedThis%pointCharges(mergedThis%numberOfPointCharges))
    allocate(mergedThis%allParticles(mergedThis%numberOfParticles))
       
    !! Set the particles manager (all pointers)        
    counter = 1
    do i = 1, mergedThis%numberOfQuantumSpecies
       do j = 1, size(mergedThis%species(i)%particles)        
          mergedThis%allParticles(counter)%particlePtr => mergedThis%species(i)%particles(j)
          counter = counter + 1
       end do
    end do

    do i = 1, mergedThis%numberOfPointCharges       
       mergedThis%allParticles(counter)%particlePtr => mergedThis%pointCharges(i)
       counter = counter + 1
    end do

    !!Copies point charges particles information
    do i = 1, mergedThis%numberOfPointCharges
       mergedThis%pointCharges(i)%name=thisA%pointCharges(i)%name
       mergedThis%pointCharges(i)%symbol=thisA%pointCharges(i)%symbol
       mergedThis%pointCharges(i)%nickname=thisA%pointCharges(i)%nickname
       mergedThis%pointCharges(i)%statistics=thisA%pointCharges(i)%statistics
       mergedThis%pointCharges(i)%basisSetName=thisA%pointCharges(i)%basisSetName
       mergedThis%pointCharges(i)%origin=thisA%pointCharges(i)%origin
       mergedThis%pointCharges(i)%charge=thisA%pointCharges(i)%charge
       mergedThis%pointCharges(i)%mass=thisA%pointCharges(i)%mass
       mergedThis%pointCharges(i)%omega=thisA%pointCharges(i)%omega
       mergedThis%pointCharges(i)%qdoCenterOf=thisA%pointCharges(i)%qdoCenterOf
       mergedThis%pointCharges(i)%spin=thisA%pointCharges(i)%spin
       mergedThis%pointCharges(i)%totalCharge=thisA%pointCharges(i)%totalCharge
       mergedThis%pointCharges(i)%klamt=thisA%pointCharges(i)%klamt
       mergedThis%pointCharges(i)%vanderWaalsRadio=thisA%pointCharges(i)%vanderWaalsRadio
       mergedThis%pointCharges(i)%isQuantum=thisA%pointCharges(i)%isQuantum
       mergedThis%pointCharges(i)%isDummy=thisA%pointCharges(i)%isDummy
       mergedThis%pointCharges(i)%fixComponent=thisA%pointCharges(i)%fixComponent
       mergedThis%pointCharges(i)%isCenterOfOptimization=thisA%pointCharges(i)%isCenterOfOptimization
       mergedThis%pointCharges(i)%multiplicity=thisA%pointCharges(i)%multiplicity
       mergedThis%pointCharges(i)%subsystem=thisA%pointCharges(i)%subsystem
       mergedThis%pointCharges(i)%translationCenter=thisA%pointCharges(i)%translationCenter
       mergedThis%pointCharges(i)%rotationPoint=thisA%pointCharges(i)%rotationPoint
       mergedThis%pointCharges(i)%rotateAround=thisA%pointCharges(i)%rotateAround
       mergedThis%pointCharges(i)%id=thisA%pointCharges(i)%id
       mergedThis%pointCharges(i)%internalSize=thisA%pointCharges(i)%internalSize
       mergedThis%pointCharges(i)%owner=thisA%pointCharges(i)%owner
       mergedThis%pointCharges(i)%basisSetSize=thisA%pointCharges(i)%basisSetSize

       if ( allocated(thisA%pointCharges(i)%childs) ) then
          allocate(mergedThis%pointCharges(i)%childs( size(thisA%pointCharges(i)%childs)))
          mergedThis%pointCharges(i)%childs=thisA%pointCharges(i)%childs
       end if

    end do

    do i = 1, mergedThis%numberOfQuantumSpecies
       j=0
       !!Copies common quantum particles information first
       do common=0,1
          !!Copies quantum particles information from system A
          do jj = 1, size(thisA%species(i)%particles)
             if( notCommonListA(i)%values(jj) .eq. common) then
                j=j+1
                mergedThis%species(i)%particles(j)%name=thisA%species(i)%particles(jj)%name
                mergedThis%species(i)%particles(j)%symbol=thisA%species(i)%particles(jj)%symbol
                mergedThis%species(i)%particles(j)%nickname=thisA%species(i)%particles(jj)%nickname
                mergedThis%species(i)%particles(j)%statistics=thisA%species(i)%particles(jj)%statistics
                mergedThis%species(i)%particles(j)%basisSetName=thisA%species(i)%particles(jj)%basisSetName
                mergedThis%species(i)%particles(j)%origin=thisA%species(i)%particles(jj)%origin
                mergedThis%species(i)%particles(j)%charge=thisA%species(i)%particles(jj)%charge
                mergedThis%species(i)%particles(j)%mass=thisA%species(i)%particles(jj)%mass
                mergedThis%species(i)%particles(j)%omega=thisA%species(i)%particles(jj)%omega
                mergedThis%species(i)%particles(j)%qdoCenterOf=thisA%species(i)%particles(jj)%qdoCenterOf
                mergedThis%species(i)%particles(j)%spin=thisA%species(i)%particles(jj)%spin
                mergedThis%species(i)%particles(j)%totalCharge=thisA%species(i)%particles(jj)%totalCharge
                mergedThis%species(i)%particles(j)%klamt=thisA%species(i)%particles(jj)%klamt
                mergedThis%species(i)%particles(j)%vanderWaalsRadio=thisA%species(i)%particles(jj)%vanderWaalsRadio
                mergedThis%species(i)%particles(j)%isQuantum=thisA%species(i)%particles(jj)%isQuantum
                mergedThis%species(i)%particles(j)%isDummy=thisA%species(i)%particles(jj)%isDummy
                mergedThis%species(i)%particles(j)%fixComponent=thisA%species(i)%particles(jj)%fixComponent
                mergedThis%species(i)%particles(j)%isCenterOfOptimization=thisA%species(i)%particles(jj)%isCenterOfOptimization
                mergedThis%species(i)%particles(j)%multiplicity=thisA%species(i)%particles(jj)%multiplicity
                mergedThis%species(i)%particles(j)%subsystem=thisA%species(i)%particles(jj)%subsystem
                mergedThis%species(i)%particles(j)%translationCenter=thisA%species(i)%particles(jj)%translationCenter
                mergedThis%species(i)%particles(j)%rotationPoint=thisA%species(i)%particles(jj)%rotationPoint
                mergedThis%species(i)%particles(j)%rotateAround=thisA%species(i)%particles(jj)%rotateAround
                mergedThis%species(i)%particles(j)%id=thisA%species(i)%particles(jj)%id
                mergedThis%species(i)%particles(j)%internalSize=thisA%species(i)%particles(jj)%internalSize
                mergedThis%species(i)%particles(j)%owner=thisA%species(i)%particles(jj)%owner
                mergedThis%species(i)%particles(j)%basisSetSize=thisA%species(i)%particles(jj)%basisSetSize

                if ( allocated(thisA%species(i)%particles(jj)%childs) ) then
                   allocate(mergedThis%species(i)%particles(j)%childs( size(thisA%species(i)%particles(jj)%childs)))
                   mergedThis%species(i)%particles(j)%childs=thisA%species(i)%particles(jj)%childs
                end if

                !!Basis information 
                mergedThis%species(i)%particles(j)%basis%name=thisA%species(i)%particles(jj)%basis%name
                mergedThis%species(i)%particles(j)%basis%origin=thisA%species(i)%particles(jj)%basis%origin
                mergedThis%species(i)%particles(j)%basis%length=thisA%species(i)%particles(jj)%basis%length
                mergedThis%species(i)%particles(j)%basis%ttype=thisA%species(i)%particles(jj)%basis%ttype
                mergedThis%species(i)%particles(j)%basis%contractionLength=thisA%species(i)%particles(jj)%basis%contractionLength
                mergedThis%species(i)%particles(j)%basis%numberOfPrimitives=thisA%species(i)%particles(jj)%basis%numberOfPrimitives

                allocate(mergedThis%species(i)%particles(j)%basis%contraction(thisA%species(i)%particles(jj)%basis%length))

                do k=1, thisA%species(i)%particles(jj)%basis%length
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%id=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%id
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%length=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%length
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%angularMoment=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%angularMoment
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%numCartesianOrbital=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%owner=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%owner
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%subsystem=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%subsystem
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%origin=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%origin

                   allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%orbitalExponents(&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%length))
                   allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%contractionCoefficients(&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%length))

                   mergedThis%species(i)%particles(j)%basis%contraction(k)%orbitalExponents=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%orbitalExponents
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%contractionCoefficients=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%contractionCoefficients

                   allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%contNormalization(&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital))
                   allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%primNormalization(&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%length,&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital))

                   mergedThis%species(i)%particles(j)%basis%contraction(k)%contNormalization=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%contNormalization
                   mergedThis%species(i)%particles(j)%basis%contraction(k)%primNormalization=&
                        thisA%species(i)%particles(jj)%basis%contraction(k)%primNormalization                  
                end do
             end if
          end do
       end do

       !!Copies quantum particles information from system B if the particle coordinates were displaced
       do jj=1, size(thisB%species(i)%particles)        
          ! if( thisB%species(i)%particles(jj)%translationCenter .ne. 0 .or. thisB%species(i)%particles(jj)%rotateAround .ne. 0) then
          if( notCommonListB(i)%values(jj) .eq. 1) then
             j=j+1
             mergedThis%species(i)%particles(j)%name=thisB%species(i)%particles(jj)%name
             mergedThis%species(i)%particles(j)%symbol=thisB%species(i)%particles(jj)%symbol
             mergedThis%species(i)%particles(j)%nickname=thisB%species(i)%particles(jj)%nickname
             mergedThis%species(i)%particles(j)%statistics=thisB%species(i)%particles(jj)%statistics
             mergedThis%species(i)%particles(j)%basisSetName=thisB%species(i)%particles(jj)%basisSetName
             mergedThis%species(i)%particles(j)%origin=thisB%species(i)%particles(jj)%origin
             mergedThis%species(i)%particles(j)%charge=thisB%species(i)%particles(jj)%charge
             mergedThis%species(i)%particles(j)%mass=thisB%species(i)%particles(jj)%mass
             mergedThis%species(i)%particles(j)%omega=thisB%species(i)%particles(jj)%omega
             mergedThis%species(i)%particles(j)%qdoCenterOf=thisB%species(i)%particles(jj)%qdoCenterOf
             mergedThis%species(i)%particles(j)%spin=thisB%species(i)%particles(jj)%spin
             mergedThis%species(i)%particles(j)%totalCharge=thisB%species(i)%particles(jj)%totalCharge
             mergedThis%species(i)%particles(j)%klamt=thisB%species(i)%particles(jj)%klamt
             mergedThis%species(i)%particles(j)%vanderWaalsRadio=thisB%species(i)%particles(jj)%vanderWaalsRadio
             mergedThis%species(i)%particles(j)%isQuantum=thisB%species(i)%particles(jj)%isQuantum
             mergedThis%species(i)%particles(j)%isDummy=thisB%species(i)%particles(jj)%isDummy
             mergedThis%species(i)%particles(j)%fixComponent=thisB%species(i)%particles(jj)%fixComponent
             mergedThis%species(i)%particles(j)%isCenterOfOptimization=thisB%species(i)%particles(jj)%isCenterOfOptimization
             mergedThis%species(i)%particles(j)%multiplicity=thisB%species(i)%particles(jj)%multiplicity
             mergedThis%species(i)%particles(j)%subsystem=thisB%species(i)%particles(jj)%subsystem
             mergedThis%species(i)%particles(j)%translationCenter=thisB%species(i)%particles(jj)%translationCenter
             mergedThis%species(i)%particles(j)%rotationPoint=thisB%species(i)%particles(jj)%rotationPoint
             mergedThis%species(i)%particles(j)%rotateAround=thisB%species(i)%particles(jj)%rotateAround
             mergedThis%species(i)%particles(j)%id=thisB%species(i)%particles(jj)%id
             mergedThis%species(i)%particles(j)%internalSize=thisB%species(i)%particles(jj)%internalSize
             mergedThis%species(i)%particles(j)%owner=thisB%species(i)%particles(jj)%owner
             mergedThis%species(i)%particles(j)%basisSetSize=thisB%species(i)%particles(jj)%basisSetSize

             if ( allocated(thisB%species(i)%particles(jj)%childs) ) then
                allocate(mergedThis%species(i)%particles(j)%childs( size(thisB%species(i)%particles(jj)%childs)))
                mergedThis%species(i)%particles(j)%childs=thisB%species(i)%particles(jj)%childs
             end if

             !!Basis information 
             mergedThis%species(i)%particles(j)%basis%name=thisB%species(i)%particles(jj)%basis%name
             mergedThis%species(i)%particles(j)%basis%origin=thisB%species(i)%particles(jj)%basis%origin
             mergedThis%species(i)%particles(j)%basis%length=thisB%species(i)%particles(jj)%basis%length
             mergedThis%species(i)%particles(j)%basis%ttype=thisB%species(i)%particles(jj)%basis%ttype
             mergedThis%species(i)%particles(j)%basis%contractionLength=thisB%species(i)%particles(jj)%basis%contractionLength
             mergedThis%species(i)%particles(j)%basis%numberOfPrimitives=thisB%species(i)%particles(jj)%basis%numberOfPrimitives

             allocate(mergedThis%species(i)%particles(j)%basis%contraction(thisB%species(i)%particles(jj)%basis%length))

             do k=1, thisB%species(i)%particles(jj)%basis%length
                mergedThis%species(i)%particles(j)%basis%contraction(k)%id=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%id
                mergedThis%species(i)%particles(j)%basis%contraction(k)%length=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%length
                mergedThis%species(i)%particles(j)%basis%contraction(k)%angularMoment=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%angularMoment
                mergedThis%species(i)%particles(j)%basis%contraction(k)%numCartesianOrbital=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital
                mergedThis%species(i)%particles(j)%basis%contraction(k)%owner=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%owner
                mergedThis%species(i)%particles(j)%basis%contraction(k)%subsystem=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%subsystem
                mergedThis%species(i)%particles(j)%basis%contraction(k)%origin=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%origin

                allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%orbitalExponents(&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%length))
                allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%contractionCoefficients(&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%length))

                mergedThis%species(i)%particles(j)%basis%contraction(k)%orbitalExponents=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%orbitalExponents
                mergedThis%species(i)%particles(j)%basis%contraction(k)%contractionCoefficients=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%contractionCoefficients

                allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%contNormalization(&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital))
                allocate(mergedThis%species(i)%particles(j)%basis%contraction(k)%primNormalization(&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%length,&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%numCartesianOrbital))

                mergedThis%species(i)%particles(j)%basis%contraction(k)%contNormalization=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%contNormalization
                mergedThis%species(i)%particles(j)%basis%contraction(k)%primNormalization=&
                     thisB%species(i)%particles(jj)%basis%contraction(k)%primNormalization                  
             end do
          end if
       end do
    end do

    ! particleManager_instance => mergedThis%allParticles

    ! if( (.not. present(sysAbasisList)) .and. (.not. present(sysBbasisList)) ) return
    
    !!Fill the basis set lists
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
       call Vector_constructorInteger(sysAbasisList(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID,mergedThis), 0 )
       call Vector_constructorInteger(sysBbasisList(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID,mergedThis), 0 )
       
       if(allocated(auxMuPositions)) deallocate(auxMuPositions)
       allocate(auxMuPositions(size(mergedThis%species(speciesID)%particles)+1))
       auxMuPositions(:)=0
       !saving merged system basis lengths
       do i = 1, size(mergedThis%species(speciesID)%particles)
          auxMuPositions(i+1)=auxMuPositions(i)
          do l = 1, size(mergedThis%species(speciesID)%particles(i)%basis%contraction)
             auxMuPositions(i+1)=auxMuPositions(i+1)+mergedThis%species(speciesID)%particles(i)%basis%contraction(l)%numCartesianOrbital
          end do
          ! print *, "merged particle", i, "start", auxMuPositions(i), "end", auxMuPositions(i+1)
       end do

       !SysI loop
       if(allocated(auxNuPositions)) deallocate(auxNuPositions)
       allocate(auxNuPositions(size(thisA%species(speciesID)%particles)+1))
       auxNuPositions(:)=0
       do i = 1, size(thisA%species(speciesID)%particles)
          auxNuPositions(i+1)=auxNuPositions(i)
          do l = 1, size(thisA%species(speciesID)%particles(i)%basis%contraction)
             auxNuPositions(i+1)=auxNuPositions(i+1)+thisA%species(speciesID)%particles(i)%basis%contraction(l)%numCartesianOrbital
          end do
          ! print *, "sysA particle", i, "start", auxNuPositions(i), "end", auxNuPositions(i+1)
       end do

       !Assign equivalence positions
       do i = 1, size(mergedThis%species(speciesID)%particles)
          do j = 1, size(thisA%species(speciesID)%particles)
             if ( MolecularSystem_checkParticleEquivalence( &
                  mergedThis%species(speciesID)%particles(i), thisA%species(speciesID)%particles(j)) .eqv. .true. ) then
                mu=auxMuPositions(i)
                do nu=auxNuPositions(j)+1,auxNuPositions(j+1)
                   mu=mu+1
                   sysAbasisList(speciesID)%values(mu)=nu
                   ! print *, "sysA", nu, "is equivalent to merged", mu
                end do
             end if
          end do
       end do

       !SysII loop
       if(allocated(auxNuPositions)) deallocate(auxNuPositions)
       allocate(auxNuPositions(size(thisB%species(speciesID)%particles)+1))
       auxNuPositions(:)=0
       do i = 1, size(thisB%species(speciesID)%particles)
          auxNuPositions(i+1)=auxNuPositions(i)
          do l = 1, size(thisB%species(speciesID)%particles(i)%basis%contraction)
             auxNuPositions(i+1)=auxNuPositions(i+1)+thisB%species(speciesID)%particles(i)%basis%contraction(l)%numCartesianOrbital
          end do
          ! print *, "sysB particle", i, "start", auxNuPositions(i), "end", auxNuPositions(i+1)
       end do

       !Assign equivalence positions
       do i = 1, size(mergedThis%species(speciesID)%particles)
          do j = 1, size(thisB%species(speciesID)%particles)
             if ( MolecularSystem_checkParticleEquivalence( &
                  mergedThis%species(speciesID)%particles(i), thisB%species(speciesID)%particles(j)) .eqv. .true. ) then
                mu=auxMuPositions(i)
                do nu=auxNuPositions(j)+1,auxNuPositions(j+1)
                   mu=mu+1
                   sysBbasisList(speciesID)%values(mu)=nu
                   ! print *, "sysB", nu, "is equivalent to merged", mu
                end do
             end if
          end do
       end do

       ! print *, "species", speciesID, "systemA list"
       ! print *, sysAbasisList(speciesID)%values
       ! print *, "species", speciesID, "systemB list"
       ! print *, sysBbasisList(speciesID)%values
    end do
    
  end subroutine MolecularSystem_mergeTwoSystems

  !>
  !! @brief Computes the displacement between equivalent particles (basis set centers) of two molecular systems
  !! @param thisA,thisB: molecular system, distanceVector: displacement of each species particles
  !! @author F. M. Moncada 2022
  subroutine MolecularSystem_getTwoSystemsDisplacement(thisA,thisB,displacementVector)
    type(MolecularSystem), intent(in) :: thisA, thisB
    type(Vector) :: displacementVector(*)
    real(8) :: distance

    type(IVector1) :: skip !avoid computing distance to the same particle of B
    real(8) :: minDistance
    integer :: minIndex
    
    integer :: nparticles
    integer :: speciesID, i, j

    
    ! print *, "max distance between equivalent particles of systems", thisA%description, thisB%description
    do speciesID=1, thisA%numberOfQuantumSpecies
       nparticles=size(thisA%species(speciesID)%particles)
       call Vector_constructor(displacementVector(speciesID),nparticles,0.0_8)
       call Vector_constructorInteger1(skip,int(nparticles,8),0_1)
       do i=1, size(thisA%species(speciesID)%particles) !ParticlesSystemA          
          minDistance=1.0E8
          minIndex=0
          do j=1, size(thisB%species(speciesID)%particles) !ParticlesSystemB 
             if(skip%values(j) .eq. 1) cycle
             
             if (thisA%species(speciesID)%particles(i)%nickname .ne. thisB%species(speciesID)%particles(j)%nickname .and. &
                  thisA%species(speciesID)%particles(i)%basisSetName .ne. thisB%species(speciesID)%particles(j)%basisSetName ) cycle
             
             distance= sqrt((thisA%species(speciesID)%particles(i)%origin(1) - thisB%species(speciesID)%particles(j)%origin(1))**2+&
                  (thisA%species(speciesID)%particles(i)%origin(2) - thisB%species(speciesID)%particles(j)%origin(2))**2+&
                  (thisA%species(speciesID)%particles(i)%origin(3) - thisB%species(speciesID)%particles(j)%origin(3))**2)

             if(distance .lt. minDistance) then
                minDistance=distance
                minIndex=j
             end if
          end do
          ! print *, "speciesID, i, closestParticle, distance", speciesID, i, minIndex, minDistance
          displacementVector(speciesID)%values(i)=minDistance
          skip%values(minIndex)=1
       end do
       ! print *, "speciesID", speciesID , "displacementVector", displacementVector(speciesID)%values
    end do
  end subroutine MolecularSystem_GetTwoSystemsDisplacement

  !>
  !! @brief Check if two particles are in the same position, are of the same speciers and have the same basis set
  !! @author F. M. Moncada, 2022
  function MolecularSystem_checkParticleEquivalence(ParticleI,ParticleII) result( output )
    implicit none
    logical :: output
    Type(Particle) :: ParticleI,ParticleII

    output =.false.
    if((ParticleI%origin(1) .eq. ParticleII%origin(1)) .and. &
         (ParticleI%origin(2) .eq. ParticleII%origin(2)) .and. &
         (ParticleI%origin(3) .eq. ParticleII%origin(3)) .and. &
         (ParticleI%basisSetName .eq. ParticleII%basisSetName) .and. &
         (ParticleI%basisSetSize .eq. ParticleII%basisSetSize) .and. &
         (ParticleI%symbol .eq. ParticleII%symbol) .and. &
         (ParticleI%nickname .eq. ParticleII%nickname) ) output =.true.

  end function MolecularSystem_checkParticleEquivalence

  !>
  !! @brief Retorna la masa total del sistema en las unidades solicitadas
  function MolecularSystem_getTotalMass( this, unid ) result( output )
    implicit none

    type(MolecularSystem), optional, target :: this
    character(*), optional :: unid
    real(8) :: output

    type(MolecularSystem), pointer :: system
    integer :: i
    
    if( present(this) ) then
       system=>this
    else
       system=>MolecularSystem_instance
    end if
    
    output = 0.0_8

    do i=1, size( molecularSystem_instance%allParticles)
       output = output + molecularSystem_instance%allParticles(i)%particlePtr%mass *  &
            molecularSystem_instance%allParticles(i)%particlePtr%internalSize
    end do

    if ( present(unid) ) then

       select case( trim(unid) )
       case ("AU")

       case("SI")
          output = output * kg

       case("AMU")
          output = output * AMU

       case default

       end select
    end if

  end function MolecularSystem_getTotalMass

  
  !>
   !! @brief  Maneja excepciones de la clase
   subroutine MolecularSystem_exception( typeMessage, description, debugDescription)
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
     
   end subroutine MolecularSystem_exception
   
 end module MolecularSystem_
