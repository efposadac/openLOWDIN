!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	PROF. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module handles all matrices for SCF program
!! @author E. F. Posada, 2013
!! @author F. Moncada, 2020
module WaveFunction_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use List_
  use Convergence_
  use MolecularSystem_
  use CosmoCore_
  use DirectIntegralManager_
  use DensityFunctionalTheory_
  use Libint2Interface_

  implicit none

  !< enum Matrix_type {
  integer, parameter :: CANONICAL_ORTHOGONALIZATION = 1
  integer, parameter :: SYMMETRIC_ORTHOGONALIZATION = 2
  !< }

  !< enum type of orbital graph {
  integer, parameter, public :: ORBITAL_ALONE = 1
  integer, parameter, public :: ORBITAL_WITH_POTENTIAL = 2
  !< }


  type, public :: WaveFunction

     !!Identity
     character(30) :: name
     integer :: species
     type(MolecularSystem), pointer :: molSys

     !!**************************************************************
     !! Matrices requeridas y alteradas en la realizacion del ciclo SCF
     !!
     type(Matrix) :: overlapMatrix     
     type(Matrix) :: transformationMatrix
     type(Matrix) :: kineticMatrix
     type(Matrix) :: puntualInteractionMatrix
     type(Matrix) :: fockMatrix     
     type(Matrix) :: densityMatrix
     type(Matrix) :: hcoreMatrix
     type(Matrix) :: twoParticlesMatrix !! Coulomb+ExchangeHF interactions between one species
     type(Matrix) :: couplingMatrix !!Sum of coulomb interaction between one species and all the different species 
     type(Matrix), allocatable :: hartreeMatrix(:) !!Coulomb Interaction between species
     type(Matrix) :: exchangeHFMatrix
     type(Matrix) :: exchangeCorrelationMatrix !!Kohn-Sham contributions. Felix: Separate this into exchange and correlation contributions
     type(Matrix) :: beforeDensityMatrix
     type(Matrix) :: waveFunctionCoefficients
     type(Vector) :: molecularOrbitalsEnergy     
     type(FourIndexMatrix), allocatable :: fourCenterIntegrals(:) !!Coulomb Interaction between species

     !! External potential contributions
     type(Matrix) :: externalPotentialMatrix
     type(Matrix) :: electricField(3)
     type(Matrix) :: harmonic
     
     !! Cosmo Things

     type(Matrix) :: cosmo1
     type(Matrix) :: cosmo2
     type(Matrix) :: cosmo4
     type(Matrix) :: cosmoCoupling
     real(8) :: cosmoCharge
     real(8) :: cosmoChargeValue

     !!**************************************************************
     !!  Variables y objetos asociados al metodo SCF
     !!
     integer :: numberOfIterations
     type(List) :: energySCF
     type(List) :: standardDesviationOfDensityMatrixElements
     type(List) :: diisError
     type(Convergence) :: convergenceMethod

     !!**************************************************************
     !! Variable por conveniencia
     real(8) :: exactExchangeFraction
     real(8) :: particlesInGrid
     integer :: removedOrbitals

     !!**************************************************************
     !! Atributos asociados a los valores de energia al final de
     !! la funcion de onda RHF
     real(8) :: totalEnergyForSpecies
     real(8) :: independentSpeciesEnergy
     real(8) :: kineticEnergy
     real(8) :: puntualInteractionEnergy
     real(8) :: independentParticleEnergy
     real(8) :: twoParticlesEnergy
     real(8) :: exchangeHFEnergy
     real(8) :: couplingEnergy
     real(8), allocatable :: hartreeEnergy(:)
     real(8) :: externalPotentialEnergy
     real(8), allocatable :: exchangeCorrelationEnergy(:)
     !! Cosmo Things
     real(8) :: cosmoEnergy
     !!**************************************************************
     
  end type WaveFunction

  type(WaveFunction), public, allocatable :: WaveFunction_instance(:)

contains

  !>
  !! @brief Define el constructor para la clase
  subroutine WaveFunction_constructor(these,nspecies,molsystem)
    implicit none
    type(WaveFunction) :: these(nspecies)
    integer :: nspecies
    type(MolecularSystem), optional, target :: molsystem

    integer :: speciesID, otherSpeciesID
    integer(8) :: numberOfContractions, otherNumberOfContractions
        
    !! Allocate memory for specie in system and load some matrices.
    do speciesID = 1, nspecies
       if( present(molsystem) ) then
          these(speciesID)%molSys=>molsystem
       else
          these(speciesID)%molSys=>MolecularSystem_instance
       end if

       these(speciesID)%species=speciesID
       these(speciesID)%name=trim(MolecularSystem_getNameOfSpecies(speciesID,these(speciesID)%molSys))

       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID,these(speciesID)%molSys)


       if(allocated(these(speciesID)%hartreeMatrix)) deallocate(these(speciesID)%hartreeMatrix)
       if(allocated(these(speciesID)%hartreeEnergy)) deallocate(these(speciesID)%hartreeEnergy)
       if(allocated(these(speciesID)%exchangeCorrelationEnergy)) deallocate(these(speciesID)%exchangeCorrelationEnergy)

       allocate(these(speciesID)%hartreeMatrix( nspecies))
       allocate(these(speciesID)%hartreeEnergy( nspecies))
       allocate(these(speciesID)%exchangeCorrelationEnergy( nspecies))


       !! Parametros Asociados con el SCF
       call List_constructor( these( speciesID )%energySCF,"energy",CONTROL_instance%LISTS_SIZE )
       call List_constructor( these( speciesID )%diisError,"diisError",CONTROL_instance%LISTS_SIZE )
       call List_constructor( these( speciesID )%standardDesviationOfDensityMatrixElements, "densitySD",CONTROL_instance%LISTS_SIZE )

       !! Instancia un objeto para manejo de aceleracion y convergencia del metodo SCF
       call Convergence_constructor(these( speciesID )%convergenceMethod, &
            these(speciesID)%name,CONTROL_instance%CONVERGENCE_METHOD,these(speciesID)%molSys)

       !! Set defaults
       these(speciesID)%totalEnergyForSpecies = 0.0_8
       these(speciesID)%independentSpeciesEnergy =0.0_8
       these(speciesID)%numberOfIterations = 0 
       these(speciesID)%kineticEnergy = 0.0_8
       these(speciesID)%puntualInteractionEnergy = 0.0_8
       these(speciesID)%independentParticleEnergy = 0.0_8
       these(speciesID)%twoParticlesEnergy = 0.0_8
       these(speciesID)%exchangeHFEnergy = 0.0_8
       these(speciesID)%externalPotentialEnergy = 0.0_8
       these(speciesID)%couplingEnergy = 0.0_8
       these(speciesID)%hartreeEnergy(:) = 0.0_8
       these(speciesID)%exchangeCorrelationEnergy(:) = 0.0_8

       !! Build some matrices
       call Matrix_constructor( these(speciesID)%overlapMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%transformationMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%kineticMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%puntualInteractionMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%fockMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%densityMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%beforeDensityMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%hcoreMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%twoParticlesMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%exchangeHFMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%couplingMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%externalPotentialMatrix, numberOfContractions, numberOfContractions, 0.0_8 )

       do otherSpeciesID = 1, nspecies
         call Matrix_constructor( these(speciesID)%hartreeMatrix(otherSpeciesID), numberOfContractions, numberOfContractions, 0.0_8 )
      end do
      
       call Matrix_constructor( these(speciesID)%exchangeCorrelationMatrix, numberOfContractions, numberOfContractions, 0.0_8 )

       call Matrix_constructor( these(speciesID)%waveFunctionCoefficients,numberOfContractions, numberOfContractions, 0.0_8 )
       call Vector_constructor( these(speciesID)%molecularOrbitalsEnergy, int(numberOfContractions) )

       !! Cosmo things
       call Matrix_constructor( these(speciesID)%cosmo1, numberOfContractions, numberOfContractions, 0.0_8 )     
       call Matrix_constructor( these(speciesID)%cosmo2, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%cosmo4, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( these(speciesID)%cosmoCoupling, numberOfContractions, numberOfContractions, 0.0_8 )

       these(speciesID)%exactExchangeFraction = 1.0_8
       these(speciesID)%particlesInGrid = 0.0
       these(speciesID)%removedOrbitals = 0

       !!Allocate arrays for integrals memory
       if (CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then
          if(allocated(these(speciesID)%fourCenterIntegrals)) deallocate(these(speciesID)%fourCenterIntegrals)

          allocate(these(speciesID)%fourCenterIntegrals(nspecies))
          !its not necessary to allocate all the species
          do otherSpeciesID=speciesID, nspecies
             otherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,these(speciesID)%molSys)
             call Matrix_fourIndexConstructor(these(speciesID)%fourCenterIntegrals(otherSpeciesID),&
                  otherNumberOfContractions,otherNumberOfContractions,numberOfContractions,numberOfContractions,0.0_8)
          end do
       end if
    end do

    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      write (*,"(T2,A15,3F12.8)") "ELECTRIC FIELD:", CONTROL_instance%ELECTRIC_FIELD
    end if

  end subroutine WaveFunction_constructor

  !>
  !! @brief Lee la matrix de overlap.
  subroutine WaveFunction_readOverlapMatrix(this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(1) = "OVERLAP"
    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    this%overlapMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
     print *,"Matriz de overlap: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
     call Matrix_show(this%overlapMatrix)
    end if
  end subroutine WaveFunction_readOverlapMatrix

  !>
  !! @brief Lee la matrix de energia cinetica.
  subroutine WaveFunction_readKineticMatrix(this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(1) = "KINETIC"
    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    this%kineticMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
     print *,"Matriz de kinetic: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
     call Matrix_show(this%kineticMatrix)
    end if
  end subroutine WaveFunction_readKineticMatrix

  !>
  !! @brief Lee la matrix de interaccion con cargas puntuales.
  subroutine WaveFunction_readPuntualInteractionMatrix(this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(1) = "ATTRACTION"
    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    this%puntualInteractionMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Matriz de puntual interaction: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
       call Matrix_show(this%puntualInteractionMatrix)
    end if
  end subroutine WaveFunction_readPuntualInteractionMatrix
  
  !>
  !! @brief Lee la matrix de interaccion con cargas puntuales.
  subroutine WaveFunction_readElectricFieldMatrices(this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    arguments(1) = "MOMENTX"
    this%electricField(1) = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    
    arguments(1) = "MOMENTY"
    this%electricField(2) = Matrix_getFromFile(rows=totalNumberOfContractions, & 
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    
    arguments(1) = "MOMENTZ"
    this%electricField(3) = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"External electric field Matrix: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
       call Matrix_show(this%electricField(1))
       call Matrix_show(this%electricField(2))
       call Matrix_show(this%electricField(3))
    end if
  end subroutine WaveFunction_readElectricFieldMatrices

  !>
  !! @brief Lee la matrix de interaccion con cargas puntuales.
  subroutine WaveFunction_readHarmonicOscillatorMatrix( this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    arguments(1) = "HARMONIC"
    this%harmonic = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Harmonic oscillator Matrix: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
       call Matrix_show(this%harmonic )
    end if

  end subroutine WaveFunction_readHarmonicOscillatorMatrix

  !! @brief Contruye la matrix de de transformacion.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_buildTransformationMatrix(this, typeOfOrthogonalization )
    implicit none
    type(WaveFunction) :: this
    integer, optional, intent(in) :: typeOfOrthogonalization

    type(Matrix) :: eigenVectors
    type(Vector) :: eigenValues
    integer(8) :: numberOfContractions
    integer :: i, j

    !! Numero de contracciones "totales"
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)

    if ( numberOfContractions > 1) then
       call Vector_constructor( eigenValues, int(numberOfContractions) )
       call Matrix_constructor( eigenVectors, numberOfContractions, numberOfContractions)
       !!****************************************************************
       !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
       !!          
       call Matrix_eigen( this%overlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )      
       ! do i = 1 , numberOfContractions
       !   print *, eigenvalues%values(i) 
       ! end do
       do i = 1 , numberOfContractions
          do j = 1 , numberOfContractions
             if ( abs(eigenValues%values(j)) >= CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
                this%transformationMatrix%values(i,j) = &
                     eigenVectors%values(i,j)/sqrt( eigenvalues%values(j) )
             else
                this%transformationMatrix%values(i,j) = 0
             end if
          end do
       end do
       do i = 1 , numberOfContractions
          if ( abs(eigenValues%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) &
               this%removedOrbitals=this%removedOrbitals+1
       end do
       if (this%removedOrbitals .gt. 0 .and. CONTROL_instance%PRINT_LEVEL .gt. 0) &
            write(*,"(A,I5,A,A,A,ES10.3)") "Removed ", this%removedOrbitals , " orbitals for species ", &
            trim(MolecularSystem_getSymbolOfSpecies(this%species,this%molSys)), " with overlap eigen threshold of ", CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
       !!
       !!****************************************************************

       !!****************************************************************
       !! Calcula matriz de transformacion
       !!
       select case (typeOfOrthogonalization)
          !! Ortogonalizacion canonica
       case (CANONICAL_ORTHOGONALIZATION)
          this%transformationMatrix%values = this%transformationMatrix%values
          !!Ortogonalizacion simetrica
       case (SYMMETRIC_ORTHOGONALIZATION)
          this%transformationMatrix%values  = &
               matmul(this%transformationMatrix%values, transpose(eigenVectors%values))
       case default
          this%transformationMatrix%values  = &
               matmul(this%transformationMatrix%values, transpose(eigenVectors%values))
       end select

       call Vector_destructor( eigenValues )
       call Matrix_destructor( eigenVectors )
       !!
       !!****************************************************************
    else
       this%transformationMatrix%values= 1.0_8
    end if

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Matriz de transformation: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
       call Matrix_show(this%transformationMatrix)
    end if

  end subroutine WaveFunction_buildTransformationMatrix

  !>
  !! @brief Contruye la matrix de particula independiente.
  subroutine WaveFunction_buildHCoreMatrix(this)
    implicit none
    type(WaveFunction) :: this
    integer :: k, l, r, s
    integer :: ParticleID, ParticleID_2
    integer :: contractionID, contractionID_2
    integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    integer :: owner, owner_2
    real(8) :: auxCharge
    real(8) :: auxOmega

    !! Incluiding mass effect       
    if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
       this%kineticMatrix%values =  &
            this%kineticMatrix%values * &
            ( 1.0_8/MolecularSystem_getMass(this%species,this%molSys) -1.0_8 / MolecularSystem_getTotalMass(this%molSys) )
    else
       this%kineticMatrix%values = &
            this%kineticMatrix%values / &
            MolecularSystem_getMass(this%species,this%molSys)
    end if

    !! Finite Nuclear Mass Correction
    if ( CONTROL_instance%FINITE_MASS_CORRECTION ) then
       k=1
       do particleID = 1, size(this%molSys%species(this%species)%particles)
          do contractionID = 1, size(this%molSys%species(this%species)%particles(particleID)%basis%contraction)

             numberOfCartesiansOrbitals = this%molSys%species(this%species)%particles(particleID)%basis%contraction(contractionID)%numCartesianOrbital
             owner = this%molSys%species(this%species)%particles(particleID)%basis%contraction(contractionID)%owner

             do s = 1, numberOfCartesiansOrbitals
                l=k

                do particleID_2 = 1, size(this%molSys%species(this%species)%particles)
                   do contractionID_2 = 1, size(this%molSys%species(this%species)%particles(particleID_2)%basis%contraction)

                      numberOfCartesiansOrbitals_2 = this%molSys%species(this%species)%particles(particleID_2)%basis%contraction(contractionID_2)%numCartesianOrbital
                      owner_2 = this%molSys%species(this%species)%particles(particleID_2)%basis%contraction(contractionID_2)%owner

                      do r = 1, numberOfCartesiansOrbitals_2

                         if ( owner .eq. owner_2) then
                            this%kineticMatrix%values(k,l)=&
                                 this%kineticMatrix%values(k,l)*&
                                 ( 1 + MolecularSystem_getMass(this%species,this%molSys) / this%molSys%species(this%species)%particles(particleID)%mass  )

                            this%kineticMatrix%values(l,k)=&
                                 this%kineticMatrix%values(k,l)
                         end if
                      end do

                      l=l+1
                   end do
                end do
                k=k+1
             end do
          end do
       end do
    end if

    !! Incluiding charge effect
    auxcharge = MolecularSystem_getCharge(this%species,this%molSys)

    this%puntualInteractionMatrix%values = &
         this%puntualInteractionMatrix%values * (-auxCharge)

    this%HCoreMatrix%values = &
         this%kineticMatrix%values + &
         this%puntualInteractionMatrix%values

    !! Add GTF external potential
    if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
       this%HCoreMatrix%values = this%HCoreMatrix%values + this%externalPotentialMatrix%values
    end if

    !! Add electric field F_i < \mu | e_i | \nu >
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
       this%HCoreMatrix%values = this%HCoreMatrix%values + auxcharge * &
            (CONTROL_instance%ELECTRIC_FIELD(1)*this%electricField(1)%values + &
            CONTROL_instance%ELECTRIC_FIELD(2)*this%electricField(2)%values + &
            CONTROL_instance%ELECTRIC_FIELD(3)*this%electricField(3)%values )
       this%externalPotentialMatrix%values = this%externalPotentialMatrix%values + auxcharge * &
            (CONTROL_instance%ELECTRIC_FIELD(1)*this%electricField(1)%values + &
            CONTROL_instance%ELECTRIC_FIELD(2)*this%electricField(2)%values + &
            CONTROL_instance%ELECTRIC_FIELD(3)*this%electricField(3)%values )
    end if


    !! Add harmonic oscillator potential 1/2 m omega**2 < \mu | r**2 | \nu >
    auxOmega = MolecularSystem_getOmega(this%species,this%molSys)

    if ( auxOmega .ne. 0.0_8 ) then
       CONTROL_instance%ARE_THERE_QDO_POTENTIALS=.true.
       this%HCoreMatrix%values = this%HCoreMatrix%values + &                                                  
            (1.0/2.0) * MolecularSystem_getMass(this%species,this%molSys) * auxOmega**2 * this%harmonic%values
       this%externalPotentialMatrix%values = this%externalPotentialMatrix%values + &
            (1.0/2.0) * MolecularSystem_getMass(this%species,this%molSys) * auxOmega**2 * this%harmonic%values
       !! Setting the kinetic and potential energy relative to the free QDO results (omega/4)
       if(CONTROL_instance%SET_QDO_ENERGY_ZERO) then
          this%HCoreMatrix%values=this%HCoreMatrix%values-auxOmega*3.0/2.0*this%overlapMatrix%values
          this%externalPotentialMatrix%values=this%externalPotentialMatrix%values-auxOmega*3.0/4.0*this%overlapMatrix%values
          this%kineticMatrix%values=this%kineticMatrix%values-auxOmega*3.0/4.0*this%overlapMatrix%values
       end if
    end if


    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Matriz de hcore: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
       call Matrix_show(this%HCoreMatrix)
    end if

  end subroutine WaveFunction_BuildHCoreMatrix

  !>
  !! @brief Calcula las componentes de energia para la especie especificada
  !! @warning Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  subroutine WaveFunction_obtainEnergyComponentsForSpecies(this)
    implicit none
    type(WaveFunction) :: this

    integer :: otherSpeciesID
    real(8) :: auxCharge

    auxcharge = MolecularSystem_getCharge(this%species,this%molSys)

    !! Remove the electric field matrix to calculate the energy components
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      this%HCoreMatrix%values = &
        this%HCoreMatrix%values - &
                auxcharge * &
        (CONTROL_instance%ELECTRIC_FIELD(1)*this%electricField(1)%values + &
         CONTROL_instance%ELECTRIC_FIELD(2)*this%electricField(2)%values + &
         CONTROL_instance%ELECTRIC_FIELD(3)*this%electricField(3)%values )
    end if
    
    !! Calcula la energia de dos particulas
    this%twoParticlesEnergy = 0.5_8 * &
         sum( transpose( this%densityMatrix%values ) * &
         this%twoParticlesMatrix%values )

    !! Calcula la energia de repulsion intraespecie
    this%hartreeEnergy(this%species) = 0.5_8 * &
         sum( transpose( this%densityMatrix%values ) * &
         this%hartreeMatrix(this%species)%values )

    !! Calcula la energia de intercambio HF
    this%exchangeHFEnergy = 0.5_8 * &
         sum( transpose( this%densityMatrix%values ) * &
         this%exchangeHFMatrix%values )

    !! Calcula energia de particula independiente
    this%independentParticleEnergy = &
         sum( transpose( this%densityMatrix%values ) * &
         this%hcoreMatrix%values )

    !! Calcula energia cinetica para la especie dada
    this%kineticEnergy = &
         sum( transpose(this%densityMatrix%values) * &
         this%kineticMatrix%values )

    !! Calcula energia de potencial externo para la especie dada
    this%externalPotentialEnergy = &
         sum( transpose(this%densityMatrix%values) * &
         this%externalPotentialMatrix%values )

    !! Calcula energia de interaccion entre particulas puntuales y cuanticas
    this%puntualInteractionEnergy =  &
         sum( transpose(this%densityMatrix%values) * &
         this%puntualInteractionMatrix%values )

    !! Calula enegia de especie independiente (  sin considerar el termino de acoplamiento )
    this%independentSpeciesEnergy = &
         sum( transpose(this%densityMatrix%values) * &
         (  ( this%hcoreMatrix%values ) + &
         0.5_8 * this%twoParticlesMatrix%values + &
         this%externalPotentialMatrix%values))

    !! Calcula energia de acoplamiento por especies
    do otherSpeciesID=1, this%molSys%numberOfQuantumSpecies
       if (this%species .ne. otherSpeciesID) then
          this%hartreeEnergy( otherSpeciesID ) = &
               sum( transpose( this%densityMatrix%values ) * &
               this%hartreeMatrix( otherSpeciesID )%values )
       end if
    end do
    
    !! Calcula energia de acoplamiento en caso de mas de una especie presente
    this%couplingEnergy = &
         sum( transpose( this%densityMatrix%values ) * &
         this%couplingMatrix%values )    

    !! Total energy for species
    this%totalEnergyForSpecies = &
         this%independentSpeciesEnergy +  &
         this%couplingEnergy + &
         sum(this%exchangeCorrelationEnergy(:))


    !! Calcula la energia COSMO	

    if(CONTROL_instance%COSMO)then

       write(*,*)"COSMO energy contributions"

       write(*,*)"Species = ",trim(this%molSys%species(this%species)%symbol)

       this%cosmoEnergy =  &
            0.5_8* (sum( transpose( this%densitymatrix%values ) * &
            this%cosmo1%values )) +0.5_8 *  &
            (sum( transpose( this%densityMatrix%values ) * &
            this%cosmo4%values )) + &
            0.5_8*(sum( transpose( this%densityMatrix%values ) * &
            this%cosmo2%values ) + &
            sum( transpose( this%densityMatrix%values ) * &
            this%cosmoCoupling%values )) 

       write(*,*)"COSMO energy 1 = ",0.5_8*(sum( transpose( this%densitymatrix%values ) * this%cosmo1%values )) 
       write(*,*)"COSMO energy 2 = ",0.5_8*(sum( transpose( this%densitymatrix%values ) * this%cosmo2%values )) 
       write(*,*)"COSMO energy 4 = ",0.5_8*(sum( transpose( this%densitymatrix%values ) * this%cosmo4%values )) 
       write(*,*)"COSMO coupling energy = ",0.5_8*(sum( transpose( this%densitymatrix%values ) * this%cosmoCoupling%values )) 
       write(*,*)"COSMO total energy = ",this%cosmoEnergy


    end if

    !! Put back the electric field matrix to the Hcore matrix
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      this%HCoreMatrix%values = &
        this%HCoreMatrix%values + &
                auxcharge * &
        (CONTROL_instance%ELECTRIC_FIELD(1)*this%electricField(1)%values + &
         CONTROL_instance%ELECTRIC_FIELD(2)*this%electricField(2)%values + &
         CONTROL_instance%ELECTRIC_FIELD(3)*this%electricField(3)%values )
    end if



    if (  CONTROL_instance%DEBUG_SCFS) then   
       print *, "__________________ ENERGY COMPONENTS _______________________"
       print *, "	Specie                       ", MolecularSystem_getNameOfSpecies(this%species,this%molSys)
       print *, "	Total Energy                =", this%totalEnergyForSpecies
       print *, "	Indepent Specie Energy      =", this%independentSpeciesEnergy
       print *, "	Kinetic Energy              =",this%kineticEnergy
       print *, "	Puntual Interaction Energy  =",this%puntualInteractionEnergy
       print *, "	Independent Particle Energy =",this%independentParticleEnergy
       print *, "	Repultion Energy            =",this%twoParticlesEnergy
       print *, "	Coupling Energy             =", this%couplingEnergy
       print *, "____________________________________________________________"
    end if
  end subroutine WaveFunction_obtainEnergyComponentsForSpecies

  !!cosmo matrices construction

  subroutine WaveFunction_cosmoHCoreMatrix(this,file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    ! integer :: k, l, r, s
    ! integer :: ParticleID, ParticleID_2
    ! integer :: contractionID, contractionID_2
    ! integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    ! integer :: owner, owner_2
    ! integer :: auxCharge
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    !! Open file
    unit = 44
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))  


    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)          


    !! Load electron potential vs clasical charges cosmo matrix
    arguments(1) = "COSMO1"    
    this%cosmo1 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)



    ! DEBUG
    ! print *,"Matriz cosmo1: ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    ! call Matrix_show( this%cosmo1 )

    !! Load clasical potential vs quantum charges cosmo matrix
    arguments(1) = "COSMO4"


    this%cosmo4 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    


    !! DEBUG
    ! print *,"Matriz cosmo 4 ", trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    ! call Matrix_show( this%cosmo4 )

    close(44)    


  end subroutine WaveFunction_cosmoHcoreMatrix

  !<
  !! @brief Lee una matriz de interaccion con un potencial externo
  !!
  !>
  subroutine WaveFunction_readExternalPotentialMatrix(this, file)
    implicit none
    type(WaveFunction) :: this
    character(*), intent(in) :: file

    integer :: unit
    integer :: totalNumberOfContractions
    character(50) :: arguments(2)

    arguments(1) = "EXTERNAL_POTENTIAL"
    arguments(2) = trim(MolecularSystem_getNameOfSpecies(this%species,this%molSys))
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    !! Get number of shells and number of cartesian contractions
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)          
    this%externalPotentialMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments(1:2))
    close(34)

    !!if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )  then  !! Numerical integration
    !!if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
    !!  this%externalPotentialMatrix = &
    !!    IntegralManager_getNumericalInteractionWithPotentialMatrix( &
    !!    ExternalPotential_Manager_instance%externalsPots, this%species, integralName="external" )
    
    !!else     !! From xml file
    !!  this%externalPotentialMatrix = &
    !!    IntegralManager_getNumericalPotentialMatrixFromXml( &
    !!    ExternalPotential_Manager_instance%externalsPots, this%species, integralName="external" )
    !!end if
    !!else    !! Analytical Integration  
    
    !!end if
    
    !! Debug
    if (  CONTROL_instance%DEBUG_SCFS) then
      print *,"EXTERNAL POTENTIAL MATRIX FOR: ", arguments(2)
      call Matrix_show(this%externalPotentialMatrix)
    end if

  end subroutine WaveFunction_readExternalPotentialMatrix
  
  
  !>
  !! @brief Builds two-particles matrix.
  subroutine WaveFunction_buildTwoParticlesMatrix( this, densityMatrixIN, factorIN, twoParticlesMatrixOUT, Libint2Objects )
    implicit none
    type(WaveFunction) :: this
    type(Matrix), optional :: densityMatrixIN
    real(8), optional :: factorIN
    type(Matrix), optional :: twoParticlesMatrixOUT
    type(Libint2Interface), optional :: Libint2Objects(:)

    real(8) :: coulomb
    real(8) :: exchange
    real(8) :: factor
    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8), allocatable, target :: tmpTwoParticlesMatrix(:,:)
    real(8), allocatable :: tmpArray(:,:)

    integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: totalNumberOfContractions
    integer :: status
    integer :: u, v, r, s, i

    type(Matrix) :: densityMatrix
    type(Matrix) :: twoParticlesMatrix

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)

    call Matrix_constructor(densityMatrix, int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), 0.0_8 )
    if ( present(densityMatrixIN)) then
       densityMatrix%values = densityMatrixIN%values
    else       
       densityMatrix%values = this%densityMatrix%values
    end if

    if ( present(factorIN)) then
       factor = MolecularSystem_getFactorOfExchangeIntegrals(this%species,this%molSys)*factorIN
    else
       factor = MolecularSystem_getFactorOfExchangeIntegrals(this%species,this%molSys)*this%exactExchangeFraction
    end if
    
    call Matrix_constructor(twoParticlesMatrix, int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), 0.0_8 )
    
    !! This matrix is only calculated if there are more than one particle for this%species or if the user want to calculate it.
    if ( MolecularSystem_getNumberOfParticles(this%species,this%molSys) > 1 .or.  CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE ) then

       if ( CONTROL_instance%INTEGRAL_STORAGE == "DISK" ) then

          !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, aa, bb, rr, ss, shellIntegrals, i, coulomb, exchange, tmpArray)
          nthreads = OMP_GET_NUM_THREADS()
          threadid =  OMP_GET_THREAD_NUM()
          unitid = 40 + threadid

          write(fileid,*) threadid
          fileid = trim(adjustl(fileid))

          if(CONTROL_instance%IS_OPEN_SHELL .and. this%molSys%species(this%species)%isElectron) then
             open( UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.ints", status='old', access='stream', form='Unformatted')
          else
             open( UNIT=unitid,FILE=trim(fileid)//trim(this%name)//".ints", status='old', access='stream', form='Unformatted')
          end if

          allocate(tmpArray(totalNumberOfContractions,totalNumberOfContractions))
          tmpArray = 0.0_8

          loadintegrals : do

             read(UNIT=unitid, iostat=status) &
                  ss(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  rr(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  bb(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  aa(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  shellIntegrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

             if(status == -1 ) then
                print*, "end of file! file: ",trim(fileid)//"E-ALPHA.ints"
                exit loadintegrals
             end if

             buildmatrix : do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                if( ss(i) == -1 ) exit loadintegrals
                ! print*, ss(i), rr(i),  bb(i), aa(i), shellIntegrals(i)
                coulomb = densityMatrix%values(rr(i),ss(i)) * shellIntegrals(i)

                !!*****************************************************************************
                !! Adds coulomb operator contributions

                if( aa(i) == rr(i) .and. bb(i) == ss(i) ) then

                   tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   if( rr(i) /= ss(i) ) then

                      tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   end if

                else

                   tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   if( rr(i) /= ss(i) ) then

                      tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   end if

                   coulomb = densityMatrix%values(aa(i),bb(i))*shellIntegrals(i)

                   tmpArray(rr(i),ss(i)) = tmpArray(rr(i),ss(i)) + coulomb

                   if ( aa(i) /= bb(i) ) then

                      tmpArray( rr(i), ss(i) ) = tmpArray( rr(i), ss(i) ) + coulomb

                   end if

                end if

                !!
                !!*****************************************************************************

                !!*****************************************************************************
                !! Adds exchange operator contributions
                !! 
                if ( abs(factor) .gt. 0.0_8 ) then

                   if( rr(i) /= ss(i) ) then

                      exchange =densityMatrix%values(bb(i),ss(i)) * shellIntegrals(i) * factor

                      tmpArray( aa(i), rr(i) ) = tmpArray( aa(i), rr(i) ) + exchange

                      if( aa(i) == rr(i) .and. bb(i) /= ss(i) ) then

                         tmpArray( aa(i), rr(i) ) = tmpArray( aa(i), rr(i) ) + exchange

                      end if

                   end if

                   if ( aa(i) /= bb(i) ) then

                      exchange = densityMatrix%values(aa(i),rr(i)) * shellIntegrals(i) * factor

                      if( bb(i) > ss(i) ) then

                         tmpArray( ss(i), bb(i) ) = tmpArray( ss(i), bb(i)) + exchange

                      else

                         tmpArray( bb(i), ss(i) ) = tmpArray( bb(i), ss(i) ) + exchange

                         if( bb(i)==ss(i) .and. aa(i) /= rr(i) ) then

                            tmpArray( bb(i), ss(i) ) = tmpArray( bb(i), ss(i) ) + exchange

                         end if

                      end if

                      if ( rr(i) /= ss(i) ) then

                         exchange = densityMatrix%values(aa(i),ss(i)) * shellIntegrals(i) * factor

                         if( bb(i) <= rr(i) ) then

                            tmpArray( bb(i), rr(i) ) = tmpArray( bb(i), rr(i) ) + exchange

                            if( bb(i) == rr(i) ) then

                               tmpArray( bb(i), rr(i) ) = tmpArray( bb(i), rr(i) ) + exchange

                            end if

                         else

                            tmpArray( rr(i), bb(i) ) = tmpArray( rr(i), bb(i)) + exchange

                            if( aa(i) == rr(i) .and. ss(i) == bb(i) ) cycle buildmatrix

                         end if

                      end if

                   end if

                   exchange = densityMatrix%values(bb(i),rr(i))*shellIntegrals(i) * factor

                   tmpArray( aa(i), ss(i) ) = tmpArray( aa(i), ss(i) ) + exchange

                end if
                !!
                !!*****************************************************************************

             end do buildmatrix

          end do loadintegrals

          close(unitid)
         
          do u = 1, totalNumberOfContractions
             do v = 1, totalNumberOfContractions
                !$OMP ATOMIC
                twoParticlesMatrix%values(u,v) = &
                     twoParticlesMatrix%values(u,v) + tmpArray(u,v) 
             end do
          end do
          
          deallocate(tmpArray)

          !$OMP END PARALLEL

          do u = 1 , totalNumberOfContractions
             do v = u + 1 , totalNumberOfContractions

                twoParticlesMatrix%values(v,u) = twoParticlesMatrix%values(v,u) + &
                     twoParticlesMatrix%values(u,v)

                twoParticlesMatrix%values(u,v) = twoParticlesMatrix%values(v,u)

             end do
          end do
          
          if ( .not. InterPotential_instance%isInstanced) &
               twoParticlesMatrix%values=twoParticlesMatrix%values*(MolecularSystem_getCharge(this%species,this%molSys))**2.0_8
          
       else if ( CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then

          !coulomb loop
          do s = 1 , totalNumberOfContractions
             !diagonal p2
             do v = 1 , totalNumberOfContractions
                !triangular unique p1
                do u = v , totalNumberOfContractions
                   twoParticlesMatrix%values(u,v)=twoParticlesMatrix%values(u,v)+&
                        densityMatrix%values(s,s)*this%fourCenterIntegrals(this%species)%values(u,v,s,s)
                end do
             end do
             !off diagonal p2
             do r = s+1 , totalNumberOfContractions
                do v = 1 , totalNumberOfContractions
                   !triangular unique p1
                   do u = v , totalNumberOfContractions
                      twoParticlesMatrix%values(u,v)=twoParticlesMatrix%values(u,v)+&
                           2.0*densityMatrix%values(r,s)*this%fourCenterIntegrals(this%species)%values(u,v,r,s)
                   end do
                end do
             end do
          end do
          !exchange loop
          do v = 1 , totalNumberOfContractions
             do r = 1 , totalNumberOfContractions
                do s = 1 , totalNumberOfContractions
                   !triangular unique p1
                   do u = v , totalNumberOfContractions
                      twoParticlesMatrix%values(u,v)=twoParticlesMatrix%values(u,v)+&
                           factor*densityMatrix%values(r,s)*this%fourCenterIntegrals(this%species)%values(u,s,r,v)
                   end do
                end do
             end do
          end do
          !symmetrize
          do v = 1 , totalNumberOfContractions
             do u = v+1 , totalNumberOfContractions
                twoParticlesMatrix%values(v,u)=twoParticlesMatrix%values(u,v)
             end do
          end do
          if ( .not. InterPotential_instance%isInstanced) &
               twoParticlesMatrix%values=twoParticlesMatrix%values*(MolecularSystem_getCharge(this%species,this%molSys))**2.0_8
          
       else !! Direct

          if ( .not. InterPotential_instance%isInstanced) then !!regular integrals
             if( present(Libint2Objects) ) then
                call DirectIntegralManager_getDirectIntraRepulsionMatrix(&
                     this%species, &
                     trim(CONTROL_instance%INTEGRAL_SCHEME), &
                     densityMatrix, & 
                     tmpTwoParticlesMatrix, &
                     factor, &
                     this%molSys, &
                     Libint2Objects(:))
             else
                call DirectIntegralManager_getDirectIntraRepulsionMatrix(&
                     this%species, &
                     trim(CONTROL_instance%INTEGRAL_SCHEME), &
                     densityMatrix, & 
                     tmpTwoParticlesMatrix, &
                     factor, &
                     this%molSys, &
                     Libint2Instance)
             end if

             tmpTwoParticlesMatrix = &
                  tmpTwoParticlesMatrix * ( MolecularSystem_getCharge(this%species,this%molSys) )**2.0_8

          else !! G12 integrals
             if( present(Libint2Objects) ) then
                call DirectIntegralManager_getDirectIntraRepulsionG12Matrix(&
                     this%species, &
                     densityMatrix, & 
                     tmpTwoParticlesMatrix, &
                     factor,&
                     this%molSys, &
                     Libint2Objects(:))
             else
                call DirectIntegralManager_getDirectIntraRepulsionG12Matrix(&
                     this%species, &
                     densityMatrix, & 
                     tmpTwoParticlesMatrix, &
                     factor,&
                     this%molSys, &
                     Libint2Instance)
             end if

          end if
          twoParticlesMatrix%values = tmpTwoParticlesMatrix
          deallocate(tmpTwoParticlesMatrix)

       end if

    end if

    if ( .not. ( present(twoParticlesMatrixOUT) )) then
       this%twoParticlesMatrix%values=twoParticlesMatrix%values
    else
       twoParticlesMatrixOUT%values=twoParticlesMatrix%values
    end if
    
    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "two particle matrix for: ", trim(this%name)
       call Matrix_show(this%twoParticlesMatrix)
    end if

  end subroutine WaveFunction_buildTwoParticlesMatrix

  !>
  !! @brief Builds the coupling matrix for the selected speciesID.
  subroutine WaveFunction_buildCouplingMatrix( these, speciesID, densityMatricesIN, couplingMatrixOUT, hartreeMatricesOUT, Libint2Objects)
    implicit none
    type(WaveFunction) :: these(*)
    integer :: speciesID
    type(Matrix), optional :: densityMatricesIN(*)
    type(Matrix), optional :: couplingMatrixOUT
    type(Matrix), optional :: hartreeMatricesOUT(*)
    type(Libint2Interface), optional :: Libint2Objects(:)

    character(30) :: nameOfSpecies
    character(30) :: nameOfOtherSpecies
    integer :: numberOfSpecies
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: otherSpeciesID
    integer :: speciesIterator
    integer :: i, j, u, v, rr, ss
    real(8), allocatable, target :: auxMatrix(:,:)
    real(8) :: coulomb

    integer :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)    

    type(Matrix), allocatable :: densityMatrices(:)
    type(Matrix), allocatable :: hartreeMatrices(:)
    type(Matrix) :: couplingMatrix
    
    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: status

    nameOfSpecies=these(speciesID)%name

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies(these(1)%molSys)
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID,these(1)%molSys)

    allocate(densityMatrices(numberOfSpecies))
    allocate(hartreeMatrices(numberOfSpecies))

    !!Initialize matrices
    do speciesIterator = 1, numberOfSpecies
       call Matrix_constructor(densityMatrices(speciesIterator), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
       call Matrix_constructor(hartreeMatrices(speciesIterator), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

       if (present(densityMatricesIN) ) then
          densityMatrices(speciesIterator)%values=densityMatricesIN(speciesIterator)%values
       else
          densityMatrices(speciesIterator)%values=these(speciesIterator)%densityMatrix%values          
       end if
    end do

    call Matrix_constructor(couplingMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

    if( numberOfSpecies > 1 ) then
       
       if ( CONTROL_instance%INTEGRAL_STORAGE == "DISK" ) then

          !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, a, b, r, s, integral, u, i, j, coulomb, auxMatrix, speciesIterator, &
          !$OMP& otherSpeciesID, nameofOtherSpecies, otherNumberOfContractions)

          nthreads = OMP_GET_NUM_THREADS()
          threadid = OMP_GET_THREAD_NUM()
          unitid = 40 + threadid

          write(fileid,*) threadid
          fileid = trim(adjustl(fileid))

          allocate(auxMatrix(numberOfContractions, numberOfContractions))
          auxMatrix=0.0_8                

          do otherSpeciesID = 1, numberOfSpecies

             nameOfOtherSpecies = MolecularSystem_getNameOfSpecies(otherSpeciesID,these(otherSpeciesID)%molSys)          
             OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,these(otherSpeciesID)%molSys)

             !! Restringe suma de terminos repulsivos de la misma especie.
             if ( otherSpeciesID .eq. speciesID ) cycle
                
             ! hartreeMatrices(otherSpeciesID)%values = 0.0_8             
             if( speciesID > otherSpeciesID) then  

                auxMatrix = 0.0_8


                !! open file for integrals
                if(CONTROL_instance%IS_OPEN_SHELL .and. &
                     these(speciesID)%molSys%species(speciesID)%isElectron .and. &
                     these(speciesID)%molSys%species(otherSpeciesID)%isElectron ) then
                   open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.E-BETA.ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else if(CONTROL_instance%IS_OPEN_SHELL .and. these(speciesID)%molSys%species(otherSpeciesID)%isElectron) then
                   open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA."//trim(nameOfSpecies)//".ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else if(CONTROL_instance%IS_OPEN_SHELL .and. these(speciesID)%molSys%species(speciesID)%isElectron) then
                   open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfOtherSpecies)//".E-ALPHA.ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else
                   open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfOtherSpecies)//"."//trim(nameOfSpecies)//".ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                end if

                readIntegrals1 : do

                   read(unitid, iostat=status)  a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   if(status == -1 ) then
                      print*, "end of file! file: ",trim(fileid)//trim(nameOfOtherSpecies)//"."//trim(nameOfSpecies)//".ints"
                      exit readIntegrals1
                   end if

                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE                      
                      if (a(u) == -1) exit readIntegrals1

                      coulomb = densityMatrices(otherSpeciesID)%values(a(u),b(u))*integral(u)

                      auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

                      if( a(u) /= b(u) ) auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

                   end do

                end do readIntegrals1

                close(unitid)

                if ( .not. InterPotential_instance%isInstanced) &
                     auxMatrix = auxMatrix * MolecularSystem_getCharge(speciesID,these(speciesID)%molSys ) * MolecularSystem_getCharge( otherSpeciesID,these(otherSpeciesID)%molSys )

                do i = 1 , numberOfContractions
                   do j = i , numberOfContractions
                      !$OMP ATOMIC
                      hartreeMatrices(otherSpeciesID)%values(i,j) = &
                           hartreeMatrices(otherSpeciesID)%values(i,j) + auxMatrix(i,j)
                   end do
                end do

             else

                auxMatrix=0.0_8

                !! open file for integrals
                if(CONTROL_instance%IS_OPEN_SHELL .and. &
                     these(speciesID)%molSys%species(speciesID)%isElectron .and. &
                     these(speciesID)%molSys%species(otherSpeciesID)%isElectron ) then
                   open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.E-BETA.ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else if(CONTROL_instance%IS_OPEN_SHELL .and. these(speciesID)%molSys%species(otherSpeciesID)%isElectron) then
                   open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//".E-ALPHA.ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else if(CONTROL_instance%IS_OPEN_SHELL .and. these(speciesID)%molSys%species(speciesID)%isElectron) then
                   open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA."//trim(nameOfOtherSpecies)//".ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                else
                   open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)//".ints", &
                        STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                end if

                readIntegrals2 : do

                   read(unitid, iostat=status) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   if(status == -1 ) then
                      print*, "end of file! file: ", trim(fileid)//trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)//".ints"
                      exit readIntegrals2
                   end if

                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
                      if (a(u) == -1) exit readIntegrals2

                      coulomb = densityMatrices(otherSpeciesID)%values(r(u),s(u))*integral(u)

                      auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

                      if( r(u) /= s(u) ) auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

                   end do

                end do readIntegrals2

                close(unitid)

                if ( .not. InterPotential_instance%isInstanced) &
                     auxMatrix = auxMatrix * MolecularSystem_getCharge(speciesID,these(speciesID)%molSys ) * MolecularSystem_getCharge( otherSpeciesID,these(otherSpeciesID)%molSys )

                do i = 1 , numberOfContractions
                   do j = i , numberOfContractions
                      !$OMP ATOMIC
                      hartreeMatrices(otherSpeciesID)%values(i,j) = &
                           hartreeMatrices(otherSpeciesID)%values(i,j) + auxMatrix(i,j)
                   end do
                end do

             end if

          end do

          deallocate(auxMatrix)

          !$OMP END PARALLEL

       else if ( CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then
          do otherSpeciesID = 1, numberOfSpecies
             if ( otherSpeciesID .eq. speciesID ) cycle
             OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,these(otherSpeciesID)%molSys)
             !integral storage order
             if( speciesID < otherSpeciesID) then  
                do v = 1 , numberOfContractions
                   !triangular unique species
                   do u = v , numberOfContractions
                      do ss = 1 , OtherNumberOfContractions
                         !diagonal other species
                         hartreeMatrices(otherSpeciesID)%values(u,v)=hartreeMatrices(otherSpeciesID)%values(u,v)+&
                              densityMatrices(otherSpeciesID)%values(ss,ss)*these(speciesID)%fourCenterIntegrals(otherSpeciesID)%values(ss,ss,u,v)
                         !off diagonal other species
                         do rr = ss+1 , OtherNumberOfContractions
                            hartreeMatrices(otherSpeciesID)%values(u,v)=hartreeMatrices(otherSpeciesID)%values(u,v)+&
                                 2.0*densityMatrices(otherSpeciesID)%values(rr,ss)*these(speciesID)%fourCenterIntegrals(otherSpeciesID)%values(rr,ss,u,v)
                         end do
                      end do
                   end do
                end do
             else
                do ss = 1 , OtherNumberOfContractions
                   !diagonal other species
                   do v = 1 , numberOfContractions
                      !triangular unique species
                      do u = v , numberOfContractions
                         hartreeMatrices(otherSpeciesID)%values(u,v)=hartreeMatrices(otherSpeciesID)%values(u,v)+&
                              densityMatrices(otherSpeciesID)%values(ss,ss)*these(otherSpeciesID)%fourCenterIntegrals(speciesID)%values(u,v,ss,ss)
                      end do
                   end do
                   !off diagonal other species
                   do rr = ss+1 , OtherNumberOfContractions
                      do v = 1 , numberOfContractions
                         !triangular unique species
                         do u = v , numberOfContractions
                            hartreeMatrices(otherSpeciesID)%values(u,v)=hartreeMatrices(otherSpeciesID)%values(u,v)+&
                                 2.0*densityMatrices(otherSpeciesID)%values(rr,ss)*these(otherSpeciesID)%fourCenterIntegrals(speciesID)%values(u,v,rr,ss)
                         end do
                      end do
                   end do
                end do
             end if
             !symmetrize
             do v = 1 , numberOfContractions
                do u = v+1 , numberOfContractions
                   hartreeMatrices(otherSpeciesID)%values(v,u)=hartreeMatrices(otherSpeciesID)%values(u,v)
                end do
             end do
             if ( .not. InterPotential_instance%isInstanced) &
                  hartreeMatrices(otherSpeciesID)%values=hartreeMatrices(otherSpeciesID)%values*&
                  MolecularSystem_getCharge(speciesID,these(speciesID)%molSys)*MolecularSystem_getCharge(otherSpeciesID,these(otherSpeciesID)%molSys)
          end do

          
       else !! Direct

          do otherSpeciesID = 1, numberOfSpecies

             !! Restringe suma de terminos repulsivos de la misma especie.
             if ( otherSpeciesID .eq. speciesID ) cycle

             if ( .not. InterPotential_instance%isInstanced) then !!regular integrals
                if(present(Libint2Objects)) then
                   call DirectIntegralManager_getDirectInterRepulsionMatrix(&
                        speciesID, OtherSpeciesID, &
                        trim(CONTROL_instance%INTEGRAL_SCHEME), &
                        densityMatrices(otherSpeciesID), &
                        auxMatrix, &
                        these(speciesID)%molSys, &
                        Libint2Objects(:))
                else
                   call DirectIntegralManager_getDirectInterRepulsionMatrix(&
                        speciesID, OtherSpeciesID, &
                        trim(CONTROL_instance%INTEGRAL_SCHEME), &
                        densityMatrices(otherSpeciesID), &
                        auxMatrix, &
                        these(speciesID)%molSys, &
                        Libint2Instance)
                end if

                auxMatrix = auxMatrix *MolecularSystem_getCharge(speciesID,these(speciesID)%molSys)*MolecularSystem_getCharge(otherSpeciesID,these(otherSpeciesID)%molSys)

             else !! G12 integrals
                if(present(Libint2Objects)) then
                   call DirectIntegralManager_getDirectInterRepulsionG12Matrix(&
                        speciesID, OtherSpeciesID, &
                        densityMatrices(otherSpeciesID), &
                        auxMatrix, &
                        these(speciesID)%molSys, &
                        Libint2Objects(:))
                else
                   call DirectIntegralManager_getDirectInterRepulsionG12Matrix(&
                        speciesID, OtherSpeciesID, &
                        densityMatrices(otherSpeciesID), &
                        auxMatrix, &
                        these(speciesID)%molSys, &
                        Libint2Instance)
                end if
             end if
                
             hartreeMatrices(otherSpeciesID)%values = &
                  hartreeMatrices(otherSpeciesID)%values + auxMatrix
             
             deallocate(auxMatrix)
             
          end do
        
       end if

       ! Symmetrize
       do otherSpeciesID = 1, numberOfSpecies

          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpeciesID .eq. speciesID ) cycle
          
          do i = 1 , numberOfContractions
             do j = i , numberOfContractions
                hartreeMatrices(otherSpeciesID)%values(j,i) = &
                     hartreeMatrices(otherSpeciesID)%values(i,j)
             end do
          end do
                                       
          if ( MolecularSystem_getSymbolOfSpecies( otherSpeciesID,these(otherSpeciesID)%molSys ) .ne. CONTROL_instance%SCF_GHOST_SPECIES ) &
               couplingMatrix%values = couplingMatrix%values + hartreeMatrices(otherSpeciesID)%values 

       end do
       
    end if

    if ( .not. ( present(hartreeMatricesOUT) .or. present(couplingMatrixOUT) ) ) then
       these(speciesID)%couplingMatrix%values=couplingMatrix%values
       do otherSpeciesID = 1, numberOfSpecies
          if ( otherSpeciesID .eq. speciesID ) cycle
          these(speciesID)%hartreeMatrix(otherSpeciesID)%values=hartreeMatrices(otherSpeciesID)%values
       end do
    else
       couplingMatrixOUT%values=couplingMatrix%values
       do otherSpeciesID = 1, numberOfSpecies
          if ( otherSpeciesID .eq. speciesID ) cycle
          hartreeMatricesOUT(otherSpeciesID)%values=hartreeMatrices(otherSpeciesID)%values
       end do
    end if
    
    
    if (  CONTROL_instance%DEBUG_SCFS) then
       do otherSpeciesID = 1, numberOfSpecies
          if ( otherSpeciesID .eq. speciesID ) cycle
          nameOfOtherSpecies = MolecularSystem_getNameOfSpecies( otherSpeciesID,these(otherSpeciesID)%molSys )          
          write(*,*) "Hartree Matrix for: ", trim(nameOfSpecies), trim(nameOfOtherSpecies)
          call Matrix_show( these(speciesID)%hartreeMatrix(otherSpeciesID) )
       end do
       write(*,*) "Coupling Matrix: ", trim(nameOfSpecies)
       call Matrix_show( these(speciesID)%couplingMatrix )
    end if

 end subroutine WaveFunction_buildCouplingMatrix

  !>
  !! @brief Builds exchange correlation contributions Matrix for DFT calculations (FELIX)
 subroutine WaveFunction_readExchangeCorrelationMatrix( this, excFileIN, &
      exchangeCorrelationMatrixOUT, exchangeCorrelationEnergyOUT, particlesInGridOUT  )       
    implicit none
    type(WaveFunction) :: this
    character(*), optional :: excFileIN
    type(Matrix), optional :: exchangeCorrelationMatrixOUT
    real(8), optional :: exchangeCorrelationEnergyOUT(*)
    real(8), optional :: particlesInGridOUT

    character(50) :: excFile
    type(Matrix) :: exchangeCorrelationMatrix
    real(8),allocatable :: exchangeCorrelationEnergy(:)
    real(8) :: particlesInGrid
    
    integer :: numberOfContractions, numberOfSpecies
    integer :: otherSpeciesID
    character(50) :: otherNameOfSpecies
    
    character(50) :: labels(2)
    integer :: excUnit
    
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(this%molSys)
    allocate(exchangeCorrelationEnergy(numberOfSpecies))

    !! Open file from dft and read matrices
    excUnit = 79
    if(present(excFileIN) ) then
       excFile = excFileIN
    else
       excFile = "lowdin.densmatrix.exc"
    end if
    
    open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")

    labels(2) = trim(this%name)
    labels(1) = "NUMBER-OF-PARTICLES"
    call Vector_getFromFile(unit=excUnit, binary=.true., value=particlesInGrid, arguments= labels(1:2) )

    labels(1) = "EXCHANGE-CORRELATION-MATRIX"
    exchangeCorrelationMatrix=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
         binary=.true., arguments=labels(1:2))

    do otherSpeciesID = this%species, numberOfSpecies
       otherNameOfSpecies=trim(MolecularSystem_getNameOfSpecies(otherSpeciesID,this%molSys))
       labels(1) = "EXCHANGE-CORRELATION-ENERGY"
       labels(2) = trim(this%name)//trim(otherNameOfSpecies)
       call Vector_getFromFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(otherSpeciesID), arguments= labels )
    end do
    
    close(unit=excUnit)

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *, "Exc. Corr. Energy ", this%name
       print *, exchangeCorrelationEnergy
       print *, "Exc. Corr. Matrix for species ", this%name
       call Matrix_show(exchangeCorrelationMatrix)
    end if

    if(.not. present(excFileIN)) then
       this%particlesInGrid=particlesInGrid
       this%exchangeCorrelationMatrix%values=exchangeCorrelationMatrix%values
       do otherSpeciesID = this%species, numberOfSpecies
          this%exchangeCorrelationEnergy(otherSpeciesID)=exchangeCorrelationEnergy(otherSpeciesID)
          Wavefunction_instance(otherSpeciesID)%exchangeCorrelationEnergy(this%species)=exchangeCorrelationEnergy(otherSpeciesID)
       end do
    else
       if(present(particlesInGridOUT)) particlesInGridOUT=particlesInGrid
       if(present(exchangeCorrelationEnergyOUT)) exchangeCorrelationEnergyOUT(1:numberOfSpecies)=exchangeCorrelationEnergy(1:numberOfSpecies)
       if(present(exchangeCorrelationMatrixOUT)) exchangeCorrelationMatrixOUT=exchangeCorrelationMatrix       
    end if
  end subroutine WaveFunction_readExchangeCorrelationMatrix

  !>
  !! @brief Builds exchange correlation contributions Matrix for DFT calculations (FELIX)
 subroutine WaveFunction_getDFTContributions( these, DFTGrids, DFTGridCommonPoints, status, densityMatricesIN, &
      exchangeCorrelationMatricesOUT, exchangeCorrelationEnergyOUT, particlesInGridOUT  )       
    implicit none
    type(WaveFunction) :: these(*)
    type(Grid) :: DFTGrids(:), DFTGridCommonPoints(:,:)
    character(*) :: status
    type(Matrix), optional :: densityMatricesIN(*)
    type(Matrix), optional :: exchangeCorrelationMatricesOUT(*)
    type(Matrix), optional :: exchangeCorrelationEnergyOUT
    real(8), optional :: particlesInGridOUT(*)

    type(Matrix), allocatable :: densityMatrices(:), exchangeCorrelationMatrices(:)
    real(8), allocatable :: particlesInGrid(:)
    type(Matrix) :: energyMatrix
    integer :: numberOfSpecies, i,j
    
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(these(1)%molSys)
    allocate(densityMatrices(numberOfSpecies), exchangeCorrelationMatrices(numberOfSpecies), particlesInGrid(numberOfSpecies))
    
    if ( present(densityMatricesIN)) then
       do i=1, numberOfSpecies
          call Matrix_copyConstructor(densityMatrices(i),densityMatricesIN(i))
       end do
    else       
       do i=1, numberOfSpecies
          call Matrix_copyConstructor(densityMatrices(i),these(i)%densityMatrix)
       end do
    end if

    call Matrix_constructor(energyMatrix, int(numberOfSpecies,8), int(numberOfSpecies,8), 0.0_8 )
    
    if (status .eq. "SCF" ) then
       call DensityFunctionalTheory_SCFDFT(DFTGrids, DFTGridCommonPoints, densityMatrices, &
            exchangeCorrelationMatrices, &
            energyMatrix, &
            particlesInGrid)
    else
       do i=1, numberOfSpecies
          call Matrix_copyConstructor(exchangeCorrelationMatrices(i),these(i)%exchangeCorrelationMatrix)
          particlesInGrid(i)=these(i)%particlesInGrid
          energyMatrix%values(i,i)=these(i)%exchangeCorrelationEnergy(i)
          do j=i+1, numberOfSpecies
             energyMatrix%values(i,j)=these(i)%exchangeCorrelationEnergy(j)
          end do
       end do
       call DensityFunctionalTheory_buildFinalGrid(DFTGrids, DFTGridCommonPoints, these(1)%molSys)

       call DensityFunctionalTheory_finalDFT(DFTGrids, DFTGridCommonPoints,densityMatrices, &
            exchangeCorrelationMatrices, &
            energyMatrix, &
            particlesInGrid)
    end if

    if ( present(exchangeCorrelationMatricesOUT)) then
       do i=1, numberOfSpecies
          call Matrix_copyConstructor(exchangeCorrelationMatricesOUT(i),exchangeCorrelationMatrices(i))
       end do
    else       
       do i=1, numberOfSpecies
          call Matrix_copyConstructor(these(i)%exchangeCorrelationMatrix,exchangeCorrelationMatrices(i))
       end do
    end if

    if ( present(exchangeCorrelationEnergyOUT)) then
       do i=1, numberOfSpecies
          exchangeCorrelationEnergyOUT%values(i,i)=energyMatrix%values(i,i)
          do j=i+1, numberOfSpecies
             exchangeCorrelationEnergyOUT%values(i,j)=energyMatrix%values(i,j)
             exchangeCorrelationEnergyOUT%values(j,i)=energyMatrix%values(i,j)
          end do
       end do
    else       
       do i=1, numberOfSpecies
          these(i)%exchangeCorrelationEnergy(i)=energyMatrix%values(i,i)
          do j=i, numberOfSpecies
             these(i)%exchangeCorrelationEnergy(j)=energyMatrix%values(i,j)
             these(j)%exchangeCorrelationEnergy(i)=energyMatrix%values(i,j)
          end do
       end do
    end if

    if ( present(particlesInGridOUT)) then
       do i=1, numberOfSpecies
          particlesInGridOUT(i)=particlesInGrid(i)
       end do
    else       
       do i=1, numberOfSpecies
          these(i)%particlesInGrid=particlesInGrid(i)
       end do
    end if
    
    do i=1, numberOfSpecies
       if (  CONTROL_instance%DEBUG_SCFS) then
          print *, "Exc. Corr. Energy ", these(i)%name
          print *, these(i)%exchangeCorrelationEnergy
          print *, "Exc. Corr. Matrix for species ", these(i)%name
          call Matrix_show(these(i)%exchangeCorrelationMatrix)
       end if
    end do

    call Matrix_destructor(energyMatrix)
    do i=1, numberOfSpecies
       call Matrix_destructor(densityMatrices(i))
       call Matrix_destructor(exchangeCorrelationMatrices(i))
    end do
    
    deallocate(densityMatrices, exchangeCorrelationMatrices, particlesInGrid)

  end subroutine WaveFunction_getDFTContributions
  
  !>
  !! @brief Writes density matrices prior to a call to DFT (FELIX)
  subroutine WaveFunction_writeDensityMatricesToFile( these, densityFileOUT, densityMatricesIN )       
    implicit none
    type(WaveFunction) :: these(*)

    character(*) :: densityFileOUT
    type(Matrix), optional :: densityMatricesIN(*)

    type(Matrix), allocatable :: densityMatrices(:)
    
    integer :: speciesID, numberOfSpecies
    
    integer :: densUnit
    character(50) :: labels(2)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(these(1)%molSys)

    allocate(densityMatrices(numberOfSpecies))
    densUnit = 78

    do speciesID=1,numberOfSpecies
       if(present(densityMatricesIN)) then
          call Matrix_copyConstructor(densityMatrices(speciesID),densityMatricesIN(speciesID))
       else
          call Matrix_copyConstructor(densityMatrices(speciesID),these(speciesID)%densityMatrix)
       end if
    end do

    open(unit = densUnit, file=trim(densityFileOUT), status="replace", form="unformatted")
    labels(1) = "DENSITY-MATRIX"
    do speciesID = 1, numberOfSpecies
       labels(2) = MolecularSystem_getNameOfSpecies(speciesID,these(speciesID)%molSys)
       call Matrix_writeToFile(densityMatrices(speciesID), unit=densUnit, binary=.true., arguments = labels )
       call Matrix_destructor(densityMatrices(speciesID))
    end do
    close (densUnit)
    
    deallocate(densityMatrices)
        
  end subroutine WaveFunction_writeDensityMatricesToFile

  !>
  !! @brief Builds fock Matrix
  subroutine WaveFunction_buildFockMatrix(this)       
    implicit none
    type(WaveFunction) :: this

    ! type(Matrix)::cosmoContribution

    this%fockMatrix%values = this%hcoreMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 1 (hcore): "//trim(this%name)
       call Matrix_show(this%fockMatrix)
    end if

    !! cosmo fock matrix
    !!full coupling
    this%fockMatrix%values = this%fockMatrix%values + &
         0.5_8*(this%cosmo1%values + &
         this%cosmo4%values)+ &
         this%cosmo2%values + &
         this%cosmoCoupling%values 

    if (  CONTROL_instance%DEBUG_SCFS .and. CONTROL_instance%COSMO) then
       print *,"MATRIZ DE FOCK 1.5 (cosmo): "//trim(this%name)
       call Matrix_show(this%fockMatrix)
    end if
    !!half coupling
    ! this%fockMatrix%values = this%fockMatrix%values + &
    !      0.5_8*(this%cosmo1%values + &
    !      this%cosmo4%values)+ &
    !      this%cosmo2%values +0.5_8*( &
    !      this%cosmoCoupling%values) 

    !!without coupling
    ! this%fockMatrix%values = this%fockMatrix%values + &
    !      0.5_8*(this%cosmo1%values + &
    !      this%cosmo4%values)+ &
    !      this%cosmo2%values

    this%fockMatrix%values = this%fockMatrix%values + this%twoParticlesMatrix%values
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 2 (+ two particles): "//trim(this%name)
       call Matrix_show(this%fockMatrix)
    end if

    this%fockMatrix%values = this%fockMatrix%values + this%couplingMatrix%values
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 3 (+ coupling): "//trim(this%name)
       call Matrix_show(this%fockMatrix)
    end if

    !!!FELIX, agrega la matriz para hacer calculo DFT
    this%fockMatrix%values = this%fockMatrix%values + this%exchangeCorrelationMatrix%values
      
    if (  CONTROL_instance%DEBUG_SCFS .and. (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS")) then
       print *,"MATRIZ DE FOCK 4 (+ exchangeCorrelation): "//trim(this%name)
       call Matrix_show(this%fockMatrix)
    end if
    
  end subroutine WaveFunction_buildFockMatrix

  !>
  !! @brief Calcula la matriz de densidad para una especie especificada
  subroutine WaveFunction_buildDensityMatrix(this)
    implicit none
    type(WaveFunction) :: this

    type(Matrix) :: auxMatrix
    integer :: orderMatrix
    integer :: ocupationNumber
    integer :: i
    integer :: j
    integer :: k

    orderMatrix = size( this%densityMatrix%values, DIM = 1 )

    ocupationNumber = MolecularSystem_getOcupationNumber(this%species,this%molSys)

    this%densityMatrix%values = 0.0_8

    call Matrix_copyConstructor(auxMatrix, this%waveFunctionCoefficients)
    
    !! Segment for fractional occupations: introduce fractional occupation
    if (trim(MolecularSystem_getSymbolOfSpecies(this%species,this%molSys)) == trim(CONTROL_instance%IONIZE_SPECIES(1)) ) then
       do i=1,size(CONTROL_instance%IONIZE_MO)
          if(CONTROL_instance%IONIZE_MO(i) .gt. 0 .and. CONTROL_instance%MO_FRACTION_OCCUPATION(i) .lt. 1.0_8) &
               auxMatrix%values(:,CONTROL_instance%IONIZE_MO(i)) = &
               auxMatrix%values(:,CONTROL_instance%IONIZE_MO(i))*sqrt(CONTROL_instance%MO_FRACTION_OCCUPATION(i))
       end do
    end if
    
    do i = 1 , orderMatrix
       do j = 1 , orderMatrix
          do k = 1 , ocupationNumber

             this%densityMatrix%values(i,j) =  &
                  this%densityMatrix%values( i,j ) + &
                  ( auxMatrix%values(i,k) &
                  * auxMatrix%values(j,k) )
          end do
       end do
    end do

    this%densityMatrix%values =  MolecularSystem_getEta(this%species,this%molSys)  * this%densityMatrix%values

    !!DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Density Matrix ", trim(this%name)
       call Matrix_show(this%densityMatrix)
    end if

  end subroutine WaveFunction_buildDensityMatrix

  !>
  !! @brief Calculates total energy for one species
  subroutine WaveFunction_obtainTotalEnergyForSpecies(this)
    implicit none
    type(WaveFunction) :: this
    
    this%totalEnergyForSpecies = &
         sum(  transpose(this%densityMatrix%values) &
         *  (( this%hcoreMatrix%values ) &
         + 0.5_8 *this%twoParticlesMatrix%values &
         + this%couplingMatrix%values)) &
         + sum(this%exchangeCorrelationEnergy(:))

    if(CONTROL_instance%COSMO)then

       this%totalEnergyForSpecies =this%totalEnergyForSpecies + 0.5_8 * &
            (sum( transpose( this%densityMatrix%values ) * &
            this%cosmo1%values )+ &
            sum( transpose( this%densityMatrix%values ) * &
            this%cosmo2%values ) + &
            sum( transpose( this%densityMatrix%values ) * &
            this%cosmo4%values ) + &
            sum( transpose( this%densityMatrix%values ) * &
            this%cosmoCoupling%values ))
    end if

    if( allocated(this%externalPotentialMatrix%values) ) then

       this%totalEnergyForSpecies =this%totalEnergyForSpecies + &
       sum( transpose(this%densityMatrix%values) &
            *  this%externalPotentialMatrix%values)
    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Total energy for "// trim(this%name) //"= ", this%totalEnergyForSpecies
       print *,"Core:            ", sum(transpose(this%densityMatrix%values)*this%hcoreMatrix%values )
       print *,"Two particles:   ", sum(transpose(this%densityMatrix%values)*0.5_8 *this%twoParticlesMatrix%values )
       print *,"Coupling:        ", sum(transpose(this%densityMatrix%values)*this%couplingMatrix%values )
       print *,"DFT exch-corr:   ", sum(this%exchangeCorrelationEnergy(:))
    end if

  end subroutine WaveFunction_obtainTotalEnergyForSpecies

  !>
  !! @brief Calcula la energia total para el sistema estudiado  
  subroutine WaveFunction_obtainTotalEnergy( these, totalEnergy, totalCouplingEnergy,  cosmo3Energy)
    implicit none
    type(WaveFunction) :: these(*)
    real(8) :: totalEnergy
    real(8) :: totalCouplingEnergy
    real(8) :: cosmo3Energy

    integer :: speciesID
    integer :: otherSpeciesID

    !! cosmo
    type(surfaceSegment) :: surface_aux2

    totalCouplingEnergy = 0.0_8
    cosmo3Energy = 0.0_8

    !! Adicionado energia de interaccion entre particulas puntuales
    totalEnergy = MolecularSystem_getPointChargesEnergy(these(1)%molSys)
    !! cosmo potential nuclei-charges nuclei

    if(CONTROL_instance%COSMO)then
       call CosmoCore_lines(surface_aux2)
       call CosmoCore_filler(surface_aux2)
       call CosmoCore_nucleiPotentialNucleiCharges(surface_aux2,cosmo3Energy)
       totalEnergy=totalEnergy+cosmo3Energy
    end if

    do speciesID = 1, these(1)%molSys%numberOfQuantumSpecies

       !! Calula enegia de especie independiente ( sin considerar el termino de acoplamiento )
       these(speciesID)%independentSpeciesEnergy = &
            sum(  transpose(these(speciesID)%densityMatrix%values) &
            *  (  ( these(speciesID)%hcoreMatrix%values ) &
            + 0.5_8 * these(speciesID)%twoParticlesMatrix%values))

       these(speciesID)%independentSpeciesEnergy = &
            these(speciesID)%independentSpeciesEnergy + &
            these(speciesID)%exchangeCorrelationEnergy(speciesID)
       
       if(CONTROL_instance%COSMO)then
          
          these(speciesID)%independentSpeciesEnergy =these(speciesID)%independentSpeciesEnergy + 0.5_8 * &
               (sum( transpose( these(speciesID)%densityMatrix%values ) * &
               these(speciesID)%cosmo1%values )+ &
               sum( transpose( these(speciesID)%densityMatrix%values ) * &
               these(speciesID)%cosmo2%values ) +  &
               sum( transpose( these(speciesID)%densityMatrix%values ) * &
               these(speciesID)%cosmo4%values ) + &
               sum( transpose( these(speciesID)%densityMatrix%values ) * &
               these(speciesID)%cosmoCoupling%values))

          ! these(speciesID)%independentSpeciesEnergy =these(speciesID)%independentSpeciesEnergy + 0.5_8 * &
          !      (sum( transpose( these(speciesID)%densityMatrix%values ) * &
          !      these(speciesID)%cosmo1%values )+ &
          !      sum( transpose( these(speciesID)%densityMatrix%values ) * &
          !      these(speciesID)%cosmo2%values ) +  &
          !      sum( transpose( these(speciesID)%densityMatrix%values ) * &
          !      these(speciesID)%cosmo4%values)) 
       end if

       totalEnergy = totalEnergy + these(speciesID)%independentSpeciesEnergy

    end do

    !! Adicionar  energia de acoplamiento y recalcula matrices de acoplamiento, including E-ALPHA/E-BETA
    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies(these(1)%molSys)
      do otherSpeciesID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies(these(1)%molSys)
          totalCouplingEnergy = totalCouplingEnergy + (sum(  transpose(these(speciesID)%densityMatrix%values) &
              * (these(speciesID)%hartreeMatrix(otherSpeciesID)%values))) 
          totalEnergy = totalEnergy+these(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)
      end do
    end do

    !! Adds inter-electron species coupling energy
    ! electronicRepulsionEnergy = WaveFunction_getAlphaBetaRepulsion()

    !! Remove the UHF "repulsion" energy
    ! totalCouplingEnergy = totalCouplingEnergy - electronicRepulsionEnergy

    !! Total Energy
    totalEnergy = totalEnergy +  totalCouplingEnergy

  end subroutine WaveFunction_obtainTotalEnergy


  !>
  !! @brief calcula la energia total de acoplamiento para una especie especificada
!   function WaveFunction_getAlphaBetaRepulsion() result( output )
!     implicit none

!     real(8) :: output

!     character(30) :: nameOfSpecie
!     character(30) :: nameOfOtherSpecie
!     real(8) :: auxValue
!     ! real(8) :: auxRepulsion
!     real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
!     integer :: numberOfContractions
!     integer :: numberOfContractionsOfOtherSpecie
!     integer :: numberOfTotalContractions
!     integer :: numberOfTotalContractionsOfOtherSpecie
!     integer :: otherSpeciesID
!     ! integer :: outFile
!     integer :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
!     integer :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
!     integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
!     integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
!     integer :: u
!     integer :: m
!     real(8), allocatable, target :: auxMatrix(:,:)
! !    real(8), allocatable :: auxMatrix2(:,:)
!     ! integer :: arrayNumber

!     !! OpenMP related variables
!     character(50) :: fileid
!     integer :: nthreads
!     integer :: threadid
!     integer :: unitid
!     integer :: status

!     output =0.0_8

!     do this%species = 1, MolecularSystem_getNumberOfQuantumSpecies()
!        do otherSpeciesID = this%species+1, MolecularSystem_getNumberOfQuantumSpecies()

!           !! Restringe suma de terminos repulsivos de la misma especie.
!           if ( otherSpeciesID /= this%species ) then

!              nameOfSpecie = MolecularSystem_getNameOfSpecies(this%species)
 
!              numberOfContractions = MolecularSystem_getNumberOfContractions(this%species)
!              numberOfTotalContractions = MolecularSystem_getTotalNumberOfContractions(this%species)

!              nameOfOtherSpecie = MolecularSystem_getNameOfSpecies( otherSpeciesID )
!              numberOfContractionsOfOtherSpecie = MolecularSystem_getNumberOfContractions( otherSpeciesID )
!              numberOfTotalContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( otherSpeciesID )

!              !Restringe la suma a solo electrones
!              if(trim(nameOfSpecie)=="E-ALPHA" .and. trim(nameOfOtherSpecie)=="E-BETA") then

!               if ( .not. CONTROL_instance%INTEGRAL_STORAGE == "DIRECT" ) then
               
!                !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, auxValue, m, a, b, r, s, integral, u)

!                nthreads = OMP_GET_NUM_THREADS()
!                threadid =  OMP_GET_THREAD_NUM()
!                unitid = 40 + threadid

!                write(fileid,*) threadid
!                fileid = trim(adjustl(fileid))

!                !! open file for integrals
!                open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
!                     STATUS='OLD', ACCESS='stream', FORM='Unformatted')

!                auxValue = 0.0_8
!                m = 0

!                readIntegrals : do

!                   read(unitid, iostat=status) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

!                   if(status == -1 ) then
!                      print*, "end of file! file: ",trim(fileid)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints"
!                      exit readIntegrals
!                   end if

!                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

!                      if (a(u) == -1) exit readIntegrals

!                      m = m + 1
!                      !print *, integral(u)


!                      auxValue = auxValue +&
!                           (  this%densityMatrix%values(b(u),a(u)) &
!                           * this( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
!                           *  integral(u))

!                      if(b(u) /= a(u)) then

!                         m = m + 1

!                         auxValue = auxValue +&
!                              (  this%densityMatrix%values(b(u),a(u)) &
!                              * this( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
!                              *  integral(u))
!                      end if

!                      if(s(u) /= r(u)) then

!                         m = m + 1

!                         auxValue = auxValue +&
!                              (  this%densityMatrix%values(b(u),a(u)) &
!                              * this( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
!                              *  integral(u))
!                      end if

!                      if(b(u) /= a(u) .and. s(u) /= r(u)) then

!                         m = m + 1

!                         auxValue = auxValue +&
!                              (  this%densityMatrix%values(b(u),a(u)) &
!                              * this( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
!                              *  integral(u))
!                      end if

!                   end do

!                end do readIntegrals

!                auxValue = auxValue *  MolecularSystem_getCharge( this%species=this%species ) &
!                     * MolecularSystem_getCharge( this%species=otherSpeciesID )

!                !$OMP ATOMIC
!                output = output + auxValue

!                close(unitid)                

!                !$OMP END PARALLEL
!               else !! DIRECT

!                ! if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
!                !  call WaveFunction_exception(ERROR, "Direct integrals are not implemented in DFT yet", "trololo")
!                ! end if
          
!                call DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix(&
!                       this%species, OtherSpeciesID, &
!                       trim(CONTROL_instance%INTEGRAL_SCHEME), &
!                       this%densityMatrix, &
!                       wavefunction_instance(otherSpeciesID)%densityMatrix, &
!                       auxMatrix )

              
!                auxMatrix = auxMatrix * MolecularSystem_getCharge(this%species) * MolecularSystem_getCharge( otherSpeciesID )
!                output = output + (sum( (auxMatrix))) 

!                deallocate(auxMatrix)

!               end if !! 

!                 output = output + (sum( transpose ( this%densityMatrix%values ) * &
!                                 this%hartreeMatrix(otherSpeciesID)%values ))

!              end if
!           end if


!        end do
!     end do

!   end function WaveFunction_getAlphaBetaRepulsion

  !**
  ! @brief Retorna la matriz  de coeficientes de combinacion
  !
  !**
  function WaveFunction_getWaveFunctionCoefficients(this) result(output)
    implicit none
    type(WaveFunction) :: this
    type(Matrix) ::  output

    if ( allocated( this%waveFunctionCoefficients%values) ) then
       call Matrix_copyConstructor( output, this%waveFunctionCoefficients )
    end if

  end function WaveFunction_getWaveFunctionCoefficients

  !**
  ! @brief Retorna valores propios del sistema molecular
  !
  !**
  function WaveFunction_getMolecularOrbitalsEnergy(this) result( output )
    implicit none
    type(WaveFunction) :: this
    type(Vector) ::  output

    if ( allocated( this%molecularOrbitalsEnergy%values) ) then
       call Vector_copyConstructor( output, this%molecularOrbitalsEnergy )
    end if

  end function WaveFunction_getMolecularOrbitalsEnergy

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine WaveFunction_exception( typeMessage, description, debugDescription)
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
  end subroutine WaveFunction_exception

  !! build cosmo 2 matrix 

  subroutine WaveFunction_buildCosmo2Matrix(this)
    type(WaveFunction) :: this

    type(species) :: specieSelected

    integer, allocatable :: labels(:)
    ! real(8), allocatable :: cosmo_int(:)
    real(8), allocatable :: ints_mat_aux(:,:)
    real(8), allocatable :: cosmo2_aux(:,:)

    ! real(8), allocatable :: qe(:)
    real(8) :: cosmo_int


    integer :: g,i,ii,h,hh,j,jj,k,l,m,o,p
    integer :: iii,jjj,hhh,gg,ll,pp,oo

    integer:: auxLabelsOfContractions
    integer:: a, b, c

    specieSelected=this%molSys%species(this%species)

    open(unit=110, file=trim(this%name)//"_qq.inn", status='old', form="unformatted")
    read(110)m


    ! if(allocated(cosmo_int)) deallocate(cosmo_int)
    ! allocate(cosmo_int(m))

    ! close(unit=110)


    if(allocated(labels)) deallocate(labels)
    allocate(labels(this%molSys%species(this%species)%basisSetSize))

    if(allocated(ints_mat_aux)) deallocate(ints_mat_aux)
    allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys), MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)))


    if(allocated(cosmo2_aux)) deallocate(cosmo2_aux)
    allocate(cosmo2_aux(MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys), MolecularSystem_getTotalNumberOfContractions(this%species,this%molSys)))


    auxLabelsOfContractions = 1

    c = 0
    do a = 1, size(specieSelected%particles)
       do b = 1, size(specieSelected%particles(a)%basis%contraction)

          c = c + 1

          !!position for cartesian contractions

          labels(c) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital


       end do
    end do


    ! call Matrix_show(this%densityMatrix)

    m = 0

    ii = 0
    do g = 1, size(this%molSys%species(this%species)%particles)
       do h = 1, size(this%molSys%species(this%species)%particles(g)%basis%contraction)

          hh = h
          ii = ii + 1
          jj = ii - 1

          do i = g, size(this%molSys%species(this%species)%particles)
             do j = hh, size(this%molSys%species(this%species)%particles(i)%basis%contraction)

                jj = jj + 1

                !!saving integrals on Matrix
                do k = labels(ii), labels(ii) + (this%molSys%species(this%species)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (this%molSys%species(this%species)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      iii=0
                      do gg = 1, size(this%molSys%species(this%species)%particles)
                         do ll = 1, size(this%molSys%species(this%species)%particles(gg)%basis%contraction)

                            hhh = ll
                            iii = iii + 1
                            jjj = iii - 1
                            do p = gg, size(this%molSys%species(this%species)%particles)
                               do o = hhh, size(this%molSys%species(this%species)%particles(p)%basis%contraction)
                                  jjj = jjj + 1

                                  !!saving integrals on Matrix
                                  do pp = labels(iii), labels(iii) + (this%molSys%species(this%species)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                     do oo = labels(jjj), labels(jjj) + (this%molSys%species(this%species)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
                                        m = m + 1

                                        read(110)cosmo_int

                                        ints_mat_aux(pp, oo) =(this%densityMatrix%values(pp,oo))* cosmo_int
                                        ints_mat_aux(oo, pp) = ints_mat_aux(pp, oo)

                                     end do
                                  end do
                               end do
                               hhh = 1
                            end do

                         end do
                      end do
                      cosmo2_aux(k,l)=0.0_8
                      do pp=1,size(ints_mat_aux,DIM=1)
                         do oo=1,size(ints_mat_aux,DIM=1)
                            cosmo2_aux(k,l)=cosmo2_aux(k,l)+ints_mat_aux(pp,oo)
                            this%cosmo2%values(k,l)=cosmo2_aux(k,l)
                            this%cosmo2%values(l,k)=this%cosmo2%values(k,l)
                         end do
                      end do
                   end do
                end do
             end do
             hh = 1
          end do
       end do
    end do

    close(unit=110)

    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "COSMO 2 matrix for: ", trim(this%name)
       call Matrix_show(this%cosmo2)
    end if

  end subroutine WaveFunction_buildCosmo2Matrix


  subroutine WaveFunction_buildCosmoCoupling(this)
    type(WaveFunction) :: this

    type(species) :: specieSelected
    type(species) :: otherSpecieSelected
    character(30) :: nameOfOtherSpecie

    integer, allocatable :: labels(:)
    integer, allocatable :: otherLabels(:)
    ! real(8), allocatable :: cosmo_int(:)
    real(8), allocatable :: ints_mat_aux(:,:)
    real(8), allocatable :: cosmoCoup_aux(:,:)

    ! real(8), allocatable :: auxMatrix(:,:)

    ! real(8), allocatable :: qe(:)
    real(8) :: cosmo_int


    integer :: currentSpeciesID
    integer :: otherSpeciesID
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: speciesIterator

    integer :: g,i,ii,h,hh,j,jj,k,l,m,o,p
    integer :: iii,jjj,hhh,gg,ll,pp,oo

    integer:: auxLabelsOfContractions
    integer:: otherAuxLabelsOfContractions
    integer:: a, b, c


    currentSpeciesID = MolecularSystem_getSpeciesID(this%name,this%molSys)
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpeciesID,this%molSys)
    specieSelected=this%molSys%species(currentSpeciesID)

    if(allocated(labels)) deallocate(labels)
    allocate(labels(this%molSys%species(currentSpeciesID)%basisSetSize))

    this%cosmoCoupling%values(:,:)=0.0_8

    auxLabelsOfContractions = 1



    c = 0
    do a = 1, size(specieSelected%particles)
       do b = 1, size(specieSelected%particles(a)%basis%contraction)

          c = c + 1

          !!position for cartesian contractions

          labels(c) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital


       end do
    end do


    if( MolecularSystem_getNumberOfQuantumSpecies(this%molSys) > 1 ) then

       this%cosmoCoupling%values = 0.0_8

       do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies(this%molSys)

          otherSpeciesID = speciesIterator
          nameOfOtherSpecie = MolecularSystem_getNameOfSpecies(otherSpeciesID,this%molSys)          
          OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,this%molSys)
          otherSpecieSelected=this%molSys%species(otherSpeciesID)

          if ( otherSpeciesID /= currentSpeciesID ) then

             ! write(*,*)"hola other and current", otherSpeciesID,currentSpeciesID 

             !      wavefunction_instance(currentSpeciesID)%cosmoCoupling%values = 0.0_8

             open(unit=110, file=trim(nameOfOtherSpecie)//trim(this%name)//"_qq.cup", status='old', form="unformatted")
             ! open(unit=110, file=trim(this%name)//trim(nameOfOtherSpecie)//"_qq.cup", status='old', form="unformatted")
             read(110)m

             ! if(allocated(cosmo_int)) deallocate(cosmo_int)
             ! allocate(cosmo_int(m))


             if(allocated(otherLabels)) deallocate(otherLabels)
             allocate(otherLabels(this%molSys%species(otherSpeciesID)%basisSetSize))

             otherAuxLabelsOfContractions=1

             c = 0

             do a = 1, size(otherSpecieSelected%particles)
                do b = 1, size(otherSpecieSelected%particles(a)%basis%contraction)

                   c = c + 1

                   !!position for cartesian contractions

                   otherlabels(c) = otherAuxLabelsOfContractions
                   otherAuxLabelsOfContractions = otherAuxLabelsOfContractions + otherSpecieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital

                end do
             end do


             if(allocated(ints_mat_aux)) deallocate(ints_mat_aux)
             allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,this%molSys), MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,this%molSys)))

             ints_mat_aux=0.0_8                


             if(allocated(cosmoCoup_aux)) deallocate(cosmoCoup_aux)
             allocate(cosmoCoup_aux(MolecularSystem_getTotalNumberOfContractions(currentSpeciesID,this%molSys), MolecularSystem_getTotalNumberOfContractions(currentSpeciesID,this%molSys)))


             m = 0

             ii = 0
             do g = 1, size(this%molSys%species(currentSpeciesID)%particles)
                do h = 1, size(this%molSys%species(currentSpeciesID)%particles(g)%basis%contraction)

                   hh = h
                   ii = ii + 1
                   jj = ii - 1

                   do i = g, size(this%molSys%species(currentSpeciesID)%particles)
                      do j = hh, size(this%molSys%species(currentSpeciesID)%particles(i)%basis%contraction)

                         jj = jj + 1

                         !!saving integrals on Matrix
                         do k = labels(ii), labels(ii) + (this%molSys%species(currentSpeciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                            do l = labels(jj), labels(jj) + (this%molSys%species(currentSpeciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                               iii=0
                               do gg = 1, size(this%molSys%species(otherSpeciesID)%particles)
                                  do ll = 1, size(this%molSys%species(otherSpeciesID)%particles(gg)%basis%contraction)

                                     hhh = ll
                                     iii = iii + 1
                                     jjj = iii - 1

                                     do p = gg, size(this%molSys%species(otherSpeciesID)%particles)
                                        do o = hhh, size(this%molSys%species(otherSpeciesID)%particles(p)%basis%contraction)
                                           jjj = jjj + 1

                                           !!saving integrals on Matrix
                                           do pp = otherlabels(iii), otherlabels(iii) + (this%molSys%species(otherSpeciesID)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                              do oo = otherlabels(jjj), otherlabels(jjj) + (this%molSys%species(otherSpeciesID)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
                                                 m = m + 1

                                                 ! write(*,*)"m,cosmo_int(m),P_element,pp,oo",m,cosmo_int(m),wavefunction_instance(otherSpeciesID)%densityMatrix%values(pp,oo),pp,oo
                                                 read(110)cosmo_int
                                                 ints_mat_aux(pp, oo) =(wavefunction_instance(otherSpeciesID)%densityMatrix%values(pp,oo))* cosmo_int
                                                 ints_mat_aux(oo, pp) = ints_mat_aux(pp, oo)

                                              end do
                                           end do

                                        end do
                                        hhh = 1
                                     end do

                                  end do
                               end do
                               ! write(*,*)"m ", m
                               cosmoCoup_aux(k,l)=0.0_8
                               do pp=1,size(ints_mat_aux,DIM=1)
                                  do oo=1,size(ints_mat_aux,DIM=1)
                                     cosmoCoup_aux(k,l)=cosmoCoup_aux(k,l)+ints_mat_aux(pp,oo)
                                  end do
                               end do
                            end do
                         end do
                      end do
                      hh = 1
                   end do
                end do
             end do
             do k=1,size(cosmoCoup_aux,DIM=1)
                do l=k,size(cosmoCoup_aux,DIM=1)
                   wavefunction_instance(currentSpeciesID)%cosmoCoupling%values(k,l)=cosmoCoup_aux(k,l)+wavefunction_instance(currentSpeciesID)%cosmoCoupling%values(k,l)
                   wavefunction_instance(currentSpeciesID)%cosmoCoupling%values(l,k)=wavefunction_instance(currentSpeciesID)%cosmoCoupling%values(k,l)
                end do
             end do



             !! debug

             if (  CONTROL_instance%DEBUG_SCFS) then

                write(*,*)"cosmo Coupling = "//trim(this%name)

                call Matrix_show(wavefunction_instance(currentSpeciesID)%cosmoCoupling)

                write(*,*)"cosmo density matrix used = "//trim(nameOfOtherSpecie)

                call Matrix_show(wavefunction_instance(otherSpeciesID)%densityMatrix)

             end if
             close(unit=110)
          end if
       end do


    end if



  end subroutine WaveFunction_buildCosmoCoupling

  subroutine WaveFunction_cosmoQuantumCharge(molSys)
    type(MolecularSystem) :: molSys
    integer :: f,g,a,c,b
    integer :: m,k,l
    integer :: h,hh,i,ii,jj,j
    integer :: auxLabelsOfContractions
    integer :: numberOfPointCharges
    integer :: orderOfMatrix
    integer :: numberOfSpecies

    integer, allocatable :: labels(:)
    real(8), allocatable :: qTotalCosmo(:)
    real(8), allocatable :: qiCosmo(:)
    real(8), allocatable :: qiDensityCosmo(:,:,:)

    character(100) :: charges_file
    character(50) :: arguments(20)

    type(species) :: specieSelected
    type(Matrix) :: densityMatrix

    character(50) :: wfnFile
    integer :: wfnUnit
    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    open(unit=100,file="cosmo.clasical",status="old",form="unformatted")
    read(100)numberOfPointCharges

    if(allocated(qTotalCosmo)) deallocate(qTotalCosmo)
    allocate(qTotalCosmo(numberOfPointCharges))

    if(allocated(qiCosmo)) deallocate(qiCosmo)
    allocate(qiCosmo(numberOfPointCharges))
    qTotalCosmo(:)=0.0_8
    read(100)(qTotalCosmo(i),i=1,numberOfPointCharges)
    close(100)
    ! write(*,*)"Cosmo Clasical Charges : ", qTotalCosmo(:)
    ! write(*,*)"sum Cosmo Clasical Charges : ", sum(qTotalCosmo(:))

    numberOfSpecies = molSys%numberOfQuantumSpecies

    do f = 1, numberOfSpecies

       specieSelected=molSys%species(f)

       if(allocated(labels)) deallocate(labels)
       allocate(labels(molSys%species(f)%basisSetSize))

       orderOfMatrix = MolecularSystem_getTotalNumberOfContractions(f,molSys)

       arguments(2) = MolecularSystem_getNameOfSpecies(f,molSys)

       arguments(1) = "DENSITY"
       densityMatrix = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       auxLabelsOfContractions = 1

       c = 0
       do a = 1, size(specieSelected%particles)
          do b = 1, size(specieSelected%particles(a)%basis%contraction)

             c = c + 1

             !!position for cartesian contractions

             labels(c) = auxLabelsOfContractions
             auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital


          end do
       end do

       charges_file="cosmo"//trim( MolecularSystem_getNameOfSpecies( f,molSys ) )//".charges"
       open(unit=100, file=trim(charges_file), status='old', form="unformatted")
       read(100)m

       if(allocated(qiDensityCosmo)) deallocate(qiDensityCosmo)
       allocate(qiDensityCosmo(orderOfMatrix, orderOfMatrix,numberOfPointCharges))
       ii = 0
       do g = 1, size(molSys%species(f)%particles)
          do h = 1, size(molSys%species(f)%particles(g)%basis%contraction)
             hh = h
             ii = ii + 1
             jj = ii - 1
             do i = g, size(molSys%species(f)%particles)
                do j = hh, size(molSys%species(f)%particles(i)%basis%contraction)
                   jj = jj + 1
                   do k = labels(ii), labels(ii) + (molSys%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (molSys%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                         read(100)(qiCosmo(m),m=1,numberOfPointCharges)
                         do m=1, numberOfPointCharges
                            qiDensityCosmo(k, l, m) = densityMatrix%values(k,l)*qiCosmo(m)
                            qiDensityCosmo(l, k, m) = qiDensityCosmo(k, l, m) 
                         end do
                      end do
                   end do
                end do
                hh = 1
             end do
          end do
       end do

       close(100)

       qiCosmo(:)=0.0_8

       do m=1, numberOfPointCharges
          do k=1, orderOfMatrix
             do l=1, orderOfMatrix
                qiCosmo(m)=qiCosmo(m)+qiDensityCosmo(k, l, m)
             end do
          end do
          qTotalCosmo(m)=qTotalCosmo(m)+qiCosmo(m)
       end do

       ! write(*,*)"Cosmo Quantum Charges : ", qiCosmo(:)
       write(*,*) "COSMO Charges for ",MolecularSystem_getSymbolOfSpecies( f,molSys )," = ", sum(qiCosmo(:))

    end do
    close(wfnUnit)


    charges_file="qTotalCosmo.charges"
    open(unit=100, file=trim(charges_file), status='replace', form="unformatted")
    write(100) qTotalCosmo(:)
    close(100)
    ! Debug
    ! write(*,*)"qTotalCosmo ", qTotalCosmo(:)
    ! write(*,*)"sumqTotalCosmo ", sum(qTotalCosmo(:))

  end subroutine WaveFunction_cosmoQuantumCharge
  
  !>
  !! @brief Remove orbitals from basis functions deleted by the overlap criteria selected
  subroutine Wavefunction_removeOrbitalsBelowEigenThreshold(this)
    implicit none
    type(WaveFunction) :: this

    integer(8) :: numberOfContractions
    type(Vector) :: normCheck, normCheckSorted
    real(8) :: threshold
    integer :: i, j, mu, nu, index

    if(this%removedOrbitals .eq. 0) return
    
    numberOfContractions = MolecularSystem_getTotalnumberOfContractions(this%species,this%molSys)
    call Vector_constructor(normCheck, int(numberOfContractions,4), 0.0_8)

    do i=1, numberOfContractions
       do mu = 1 , numberOfContractions
          do nu = 1 , numberOfContractions
             normCheck%values(i)=normCheck%values(i)+this%waveFunctionCoefficients%values(mu,i)*&
                  this%waveFunctionCoefficients%values(nu,i)*&
                  this%overlapMatrix%values(mu,nu)
          end do
       end do
    end do

    normCheckSorted=normCheck
    call Vector_reverseSortElements(normCheckSorted)
    threshold=normCheckSorted%values(this%removedOrbitals)

    i=0
    do index = 1 , numberOfContractions
       i=i+1
       if ( normCheck%values(i) .le. threshold) then
          if (  CONTROL_instance%DEBUG_SCFS) &
               print *, "shifting eigenvector no.", i, "with normCheck", normCheck%values(i), "to the end of the coefficients matrix"
          do j = i , numberOfContractions-1
             this%molecularOrbitalsEnergy%values(j)=this%molecularOrbitalsEnergy%values(j+1)
             this%waveFunctionCoefficients%values(1:numberOfContractions,j) = this%waveFunctionCoefficients%values(1:numberOfContractions,j+1)
             normCheck%values(j)=normCheck%values(j+1)
          end do
          i=i-1
          ! Make eigenenergy a very large number
          this%molecularOrbitalsEnergy%values(numberOfContractions)=1.0E+308_8
          this%waveFunctionCoefficients%values(1:numberOfContractions,numberOfContractions)=0.0_8
          normCheck%values(numberOfContractions)=1.0E+308_8
       end if
    end do
    
  end subroutine Wavefunction_removeOrbitalsBelowEigenThreshold


end module WaveFunction_

