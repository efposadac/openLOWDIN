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
       MolecularSystem_getFactorOfInterchangeIntegrals, &
       MolecularSystem_getNameOfSpecie, &
       MolecularSystem_getSpecieID, &
       MolecularSystem_getPointChargesEnergy, &
       MolecularSystem_getMMPointChargesEnergy, &
       MolecularSystem_getlabelsofcontractions
  
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
  subroutine MolecularSystem_showInformation()
    implicit none
    
    type(Exception) :: ex
    
    print *,""
    print *," MOLECULAR SYSTEM: ",trim(MolecularSystem_instance%name)
    print *,"-----------------"
    print *,""
    write (6,"(T5,A16,A)") "DESCRIPTION   : ", trim( MolecularSystem_instance%description )
    write (6,"(T5,A16,I3)") "CHARGE        : ",MolecularSystem_instance%charge
    write (6,"(T5,A16,A4)") "PUNTUAL GROUP : ", "NONE"
    print *,""
    
 
  end subroutine MolecularSystem_showInformation

  !>
  !! @brief Muestra los atributos de todas las particulas en el Administrador de particulas
  !! @author S. A. Gonzalez
  !! @par changes : 
  !!     - rewritten, E. F. Posada. 2013
  subroutine MolecularSystem_showParticlesInformation()
    implicit none
    integer :: i
    integer :: j
    
    print *,""
    print *," INFORMATION OF PARTICLES :"
    print *,"========================="
    write (6,"(T10,A29,I8.0)") "Total number of particles   =", MolecularSystem_instance%numberOfQuantumParticles + MolecularSystem_instance%numberOfPointCharges
    write (6,"(T10,A29,I8.0)") "Number of quantum particles =", MolecularSystem_instance%numberOfQuantumParticles
    write (6,"(T10,A29,I8.0)") "Number of puntual charges   =", MolecularSystem_instance%numberOfPointCharges
    write (6,"(T10,A29,I8.0)") "Number of quantum species   =", MolecularSystem_instance%numberOfQuantumSpecies

    !!***********************************************************************
    !! Imprime iformacion sobre masa, carga y numero de particulas encontradas
    !!
    print *,""
    print *,"                INFORMATION OF QUANTUM SPECIES "
    write (6,"(T5,A70)") "---------------------------------------------------------------------"
    write (6,"(T10,A2,A4,A8,A10,A4,A5,A6,A5,A4,A5,A12)") "ID", " ","Symbol", " ","mass", " ","charge", " ","spin","","multiplicity"
    write (6,"(T5,A70)") "---------------------------------------------------------------------"

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       write (6,'(T8,I3.0,A5,A10,A5,F7.1,A5,F5.1,A5,F5.2,A5,F5.2)') &
            i, " ", &
            trim(MolecularSystem_instance%species(i)%symbol)," ",&
            MolecularSystem_instance%species(i)%mass," ",&
            MolecularSystem_instance%species(i)%charge, " ",&
            MolecularSystem_instance%species(i)%spin, "",&
            MolecularSystem_instance%species(i)%multiplicity
    end do

    print *,""
    print *,"                  CONSTANTS OF COUPLING "
    write (6,"(T7,A60)") "------------------------------------------------------------"
    write (6,"(T10,A8,A10,A5,A5,A4,A5,A6,A5,A9)") "Symbol", " ","kappa", " ","eta", " ","lambda","","ocupation"
    write (6,"(T7,A60)") "------------------------------------------------------------"
    
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       write (6,'(T10,A10,A5,F7.1,A5,F5.2,A5,F5.2,A5,F5.2,A5,F5.2)'), &
       trim(MolecularSystem_instance%species(i)%symbol)," ",&
            MolecularSystem_instance%species(i)%kappa, " ", &
            MolecularSystem_instance%species(i)%eta, " ",&
            MolecularSystem_instance%species(i)%lambda," ",&
            MolecularSystem_instance%species(i)%ocupationNumber
    end do

    print *,""

    print *,"                  BASIS SET FOR SPECIES "
    write (6,"(T7,A60)") "------------------------------------------------------------"
    write (6,"(T10,A8,A10,A8,A5,A12,A5,A9)") "Symbol", " ","N. Basis", " ","N. Particles"," ","Basis Set"
    write (6,"(T7,A60)") "------------------------------------------------------------"
    
    !! Only shows the basis-set name of the first particle by specie.
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       
       if( MolecularSystem_instance%species(i)%isElectron .and. CONTROL_instance%IS_OPEN_SHELL ) then

          write (6,'(T10,A10,A5,I8,A5,I12,A5,A10)') &
               trim(MolecularSystem_instance%species(i)%symbol)," ",&
               !MolecularSystem_getTotalNumberOfContractions(i)," ",&
               MolecularSystem_instance%species(i)%basisSetSize," ",&
               int(MolecularSystem_instance%species(i)%internalSize / 2), " ",&
               trim(MolecularSystem_instance%species(i)%particles(1)%basis%name)
          ! write(*,*)MolecularSystem_getTotalNumberOfContractions(i)
          
       else

          write (6,'(T10,A10,A5,I8,A5,I12,A5,A10)') &
               trim(MolecularSystem_instance%species(i)%symbol)," ",&
               !MolecularSystem_getTotalNumberOfContractions(i)," ",&
               MolecularSystem_instance%species(i)%basisSetSize," ",&
               MolecularSystem_instance%species(i)%internalSize, " ",&
               trim(MolecularSystem_instance%species(i)%particles(1)%basis%name)
          ! write(*,*)MolecularSystem_getTotalNumberOfContractions(i)
          
       end if
       
    end do

    write (6,*) ""
    write (6,"(T10,A35)")"                     ATOMIC BASIS"
    write (6,"(T10,A60)") "------------------------------------------------------------"
    write (6,"(T10,A11,A9,A20,A20)") " PRIMITIVE ", "  SHELL  "," EXPONENT "," COEFFICIENT "
    write (6,"(T10,A60)") "------------------------------------------------------------"
    
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       
       !! Avoid print twice basis in open-shell case
       if(trim(MolecularSystem_instance%species(i)%name) == "E-BETA" ) cycle
       
       write (6,*) ""
       write( 6, "(T5,A32,A5)") "BEGIN DESCRIPTION OF BASIS FOR: ", trim(MolecularSystem_instance%species(i)%symbol)
       write (6,"(T5,A30)") "================================"
       write (6,*) ""
       
       do j = 1, size( MolecularSystem_instance%species(i)%particles )
          
          call BasisSet_showInCompactForm( MolecularSystem_instance%species(i)%particles(j)%basis,&
               trim(MolecularSystem_instance%species(i)%particles(j)%nickname ))
          
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

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
          write (6,'(T10,A10,A5,I8,A5,I12)') &
               trim(MolecularSystem_instance%species(i)%symbol)," ",&
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
  subroutine MolecularSystem_showCartesianMatrix()
    implicit none
    
    integer :: i, j
    real(8) :: origin(3)
    
    write (6,"(A10,A16,A20,A20)") " ","<x>","<y>","<z>"
    
    !! Print quatum species information
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       
       !! Avoid print twice basis in open-shell case
       if(trim(MolecularSystem_instance%species(i)%name) == "E-BETA" ) cycle

       do j = 1, size(MolecularSystem_instance%species(i)%particles)
       
          origin = MolecularSystem_instance%species(i)%particles(j)%origin * AMSTRONG
          
          if(MolecularSystem_instance%species(i)%isElectron) then
             write (6,"(A10,3F20.10)") trim( MolecularSystem_instance%species(i)%particles(j)%symbol )//trim(MolecularSystem_instance%species(i)%particles(j)%nickname),&
                  origin(1), origin(2), origin(3)
          else
             write (6,"(A10,3F20.10)") trim(MolecularSystem_instance%species(i)%particles(j)%nickname), origin(1), origin(2), origin(3)
          end if
          
       end do
    end do
    
    !! Print Point charges information
    do i = 1, MolecularSystem_instance%numberOfPointCharges
       
       origin = MolecularSystem_instance%pointCharges(i)%origin * AMSTRONG
       write (6,"(A10,3F20.10)") trim(MolecularSystem_instance%pointCharges(i)%nickname), origin(1), origin(2), origin(3)
       
    end do
    
    print *," "
    
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
  subroutine MolecularSystem_saveToFile()
    implicit none
    
    integer i, j, k
    character(100) :: title
    
    !!****************************************************************************
    !! CONTROL parameters on file.
    !!
    
    !! open file
    open(unit=40, file="lowdin.dat", status="replace", form="formatted")
    
    !!save all options
    call CONTROL_save(40)
    close(40)
    
    !!****************************************************************************
    !!Save the molecular system on file.
    !!

    !!Open file
    open(unit=40, file="lowdin.sys", status="replace", form="formatted")
    
    !! Saving general information.
    write(40,*) MolecularSystem_instance%name
    write(40,*) MolecularSystem_instance%description    
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
    open(unit=40, file="lowdin.bas", status="replace", form="formatted")
    
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
  subroutine MolecularSystem_loadFromFile( form )
    implicit none
    
    character(*) :: form

    integer :: auxValue
    integer :: counter
    integer :: i, j
    logical :: existFile
    character(20) :: name
    character(50) :: species
    character(50) :: otherSpecies

    select case (trim(form))
       
    case("LOWDIN.BAS")
       
       !!****************************************************************************
       !! Loading info from the lowdin.bas format
       !!
       
       !!open file
       existFile = .false.
       inquire(file="lowdin.bas", exist=existFile)
       
       if(existFile) then

          !! Destroy the molecular system if any
          call MolecularSystem_destroy()

          open(unit=40, file="lowdin.bas", status="old", form="formatted")
          
          read(40,*) MolecularSystem_instance%numberOfQuantumSpecies
          allocate(MolecularSystem_instance%species(MolecularSystem_instance%numberOfQuantumSpecies))

          MolecularSystem_instance%numberOfQuantumParticles = 0

          do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
             
             read(40,*) MolecularSystem_instance%species(i)%name             
             read(40,*) auxValue
             
             allocate(MolecularSystem_instance%species(i)%particles(auxValue))
             
             do j = 1, size(MolecularSystem_instance%species(i)%particles)
                
                read(40,*) MolecularSystem_instance%species(i)%particles(j)%nickname
                
                MolecularSystem_instance%numberOfQuantumParticles = MolecularSystem_instance%numberOfQuantumParticles + 1
                call BasisSet_load(MolecularSystem_instance%species(i)%particles(j)%basis, "LOWDIN.BAS", unit = 40)
                
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

       else
          
          call MolecularSystem_exception(ERROR, "The file: lowdin.bas  was not found!","MolecularSystem module at LoadFromFile function.")
          
       end if
       
    case("LOWDIN.DAT")
       
       !!****************************************************************************
       !! Load CONTROL parameters from file.
       !!
       
       !! open file
       existFile = .false.
       inquire(file="lowdin.dat", exist=existFile)
       
       if(existFile) then
          
          open(unit=40, file="lowdin.dat", status="old", form="formatted")
          
          call CONTROL_start()
          call CONTROL_load(unit = 40)
        
          close(40)
          
       else

          call MolecularSystem_exception(ERROR, "The file: lowdin.dat  was not found!","MolecularSystem module at LoadFromFile function.")
          
       end if
       
    case("LOWDIN.SYS")
       
       !!****************************************************************************
       !! Load  the molecular system from file.
       !!
       
       !! Destroy the molecular system if any
       call MolecularSystem_destroy()
    
       !! open file
       existFile = .false.
       inquire(file="lowdin.sys", exist=existFile)
       
       if(existFile) then
          !!Open file
          open(unit=40, file="lowdin.sys", status="old", form="formatted")
          
          !! read general information.
          read(40,*) MolecularSystem_instance%name
          read(40,*) MolecularSystem_instance%description    
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

       else
          
          call MolecularSystem_exception(ERROR, "The file: lowdin.sys  was not found!","MolecularSystem module at LoadFromFile function.")
          
       end if
          
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
    
    integer :: i, j, k

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
  function MolecularSystem_getTotalNumberOfContractions( specieID ) result( output )
    implicit none
    integer :: specieID
    integer :: output
    
    integer :: i, j, k

    output = 0

    
    do j = 1, size(MolecularSystem_instance%species(specieID)%particles)

       do k = 1, size(MolecularSystem_instance%species(specieID)%particles(j)%basis%contraction)
          
          output = output + MolecularSystem_instance%species(specieID)%particles(j)%basis%contraction(k)%numCartesianOrbital
          
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
   function MolecularSystem_getOcupationNumber(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     integer :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%ocupationNumber
          
   end function MolecularSystem_getOcupationNumber

   !> @brief Returns the eta parameter of a species
   !! @author E. F. Posada, 2013
   !! @version 1.0
   function MolecularSystem_getEta(speciesID) result(output)
     implicit none
     
     integer :: speciesID
     integer :: output
     
     output = -1
     output = MolecularSystem_instance%species(speciesID)%eta
          
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

   !> @brief Returns the mass of speciesID
   !! @author E. F. Posada, 2013
   !! @version 1.0   
   function MolecularSystem_getMass( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%mass
     
   end function MolecularSystem_getMass
   
   !> @brief Returns the Factor Of Interchange Integrals
   !! @author E. F. Posada, 2013
   !! @version 1.0   
   function MolecularSystem_getFactorOfInterchangeIntegrals( speciesID ) result( output )
     implicit none
     integer :: speciesID
     
     real(8) :: output
     
     output = MolecularSystem_instance%species(speciesID)%kappa / MolecularSystem_instance%species(speciesID)%eta
     
   end function MolecularSystem_getFactorOfInterchangeIntegrals

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
   function MolecularSystem_getSpecieID( nameOfSpecie ) result(output)
     implicit none
     
     character(*) nameOfSpecie
     integer :: output
     integer i 
     
     output = 0
     
     do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
        if( trim(MolecularSystem_instance%species(i)%name) == trim(nameOfSpecie)) then
           output = i
        end if
     end do
          
   end function MolecularSystem_getSpecieID
   
   !>
   !! @brief calcula la energia total para una especie especificada
   function MolecularSystem_getPointChargesEnergy() result( output )
     implicit none
     real(8) :: output
     
     integer :: i
     integer :: j
     real(8) :: deltaOrigin(3)
     real(8) :: tmp
     
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
     
   end function MolecularSystem_getPointChargesEnergy

   function MolecularSystem_getMMPointChargesEnergy() result( output )
     implicit none
     real(8) :: output

     integer :: i
     integer :: j
     real(8) :: deltaOrigin(3)
     real(8) :: tmp

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
              
              write (output(counter),"(I5,A2,A6,A2,A4)") counter, "  ", &
                   trim(MolecularSystem_instance%species(speciesID)%particles(i)%nickname), "  ", &
                   trim(shellCode(k))//" "

              counter = counter + 1 
              
           end do

        end do
     end do
     
   end function MolecularSystem_getlabelsofcontractions

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
