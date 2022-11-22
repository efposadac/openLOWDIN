!!******************************************************************************
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

     character(30) :: name

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
     type(Matrix) :: externalPotentialMatrix
     type(Matrix) :: beforeDensityMatrix
     type(Matrix) :: waveFunctionCoefficients
     type(Vector) :: molecularOrbitalsEnergy     

     !! Cosmo Things

     type(Matrix) :: cosmo1
     type(Matrix) :: cosmo2
     type(Matrix) :: cosmo4
     type(Matrix) :: cosmoCoupling
     type(Matrix) :: electricField(3)
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
     real(8) :: totalEnergyForSpecie
     real(8) :: independentSpecieEnergy
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
  subroutine WaveFunction_constructor( )
    implicit none

    integer :: speciesID, i, otherSpeciesID
    integer :: statusSystem    
    integer(8) :: numberOfContractions
    character(50) :: labels(2)
    character(50) :: dftFile, vecFile
    integer :: dftUnit, vecUnit
    logical :: existFile

    !! Allocate memory.
    allocate(WaveFunction_instance(MolecularSystem_instance%numberOfQuantumSpecies))

    !! Allocate memory for specie in system and load some matrices.
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       labels = ""
       labels(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)


       allocate(WaveFunction_instance(speciesID)%hartreeMatrix( MolecularSystem_instance%numberOfQuantumSpecies))
       allocate(WaveFunction_instance(speciesID)%hartreeEnergy( MolecularSystem_instance%numberOfQuantumSpecies))
       allocate(WaveFunction_instance(speciesID)%exchangeCorrelationEnergy( MolecularSystem_instance%numberOfQuantumSpecies))


       !! Parametros Asociados con el SCF
       call List_constructor( WaveFunction_instance( speciesID )%energySCF,"energy",CONTROL_instance%LISTS_SIZE )
       call List_constructor( WaveFunction_instance( speciesID )%diisError,"diisError",CONTROL_instance%LISTS_SIZE )
       call List_constructor( WaveFunction_instance( speciesID )%standardDesviationOfDensityMatrixElements, "densitySD",CONTROL_instance%LISTS_SIZE )

       !! Instancia un objeto para manejo de aceleracion y convergencia del metodo SCF
       call Convergence_constructor(WaveFunction_instance( speciesID )%convergenceMethod, &
            WaveFunction_instance( speciesID )%name,CONTROL_instance%CONVERGENCE_METHOD)

       !! Set defaults
       WaveFunction_instance( speciesID )%totalEnergyForSpecie = 0.0_8
       WaveFunction_instance( speciesID )%independentSpecieEnergy =0.0_8
       WaveFunction_instance( speciesID )%numberOfIterations = 0 
       WaveFunction_instance( speciesID )%kineticEnergy = 0.0_8
       WaveFunction_instance( speciesID )%puntualInteractionEnergy = 0.0_8
       WaveFunction_instance( speciesID )%independentParticleEnergy = 0.0_8
       WaveFunction_instance( speciesID )%twoParticlesEnergy = 0.0_8
       WaveFunction_instance( speciesID )%exchangeHFEnergy = 0.0_8
       WaveFunction_instance( speciesID )%externalPotentialEnergy = 0.0_8
       WaveFunction_instance( speciesID )%couplingEnergy = 0.0_8
       WaveFunction_instance( speciesID )%hartreeEnergy(:) = 0.0_8
       WaveFunction_instance( speciesID )%exchangeCorrelationEnergy(:) = 0.0_8

       !! Cosmo things
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo1, numberOfContractions, numberOfContractions, 0.0_8 )     
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo4,numberOfContractions, numberOfContractions, 0.0_8 )

       call Matrix_constructor( WaveFunction_instance(speciesID)%externalPotentialMatrix, numberOfContractions, numberOfContractions, 0.0_8 )

       !! Build some matrices
       call Matrix_constructor( WaveFunction_instance(speciesID)%overlapMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%transformationMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%kineticMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%puntualInteractionMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%fockMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%densityMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%beforeDensityMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%hcoreMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%twoParticlesMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%exchangeHFMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%couplingMatrix, numberOfContractions, numberOfContractions, 0.0_8 )


       do otherSpeciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
         call Matrix_constructor( WaveFunction_instance(speciesID)%hartreeMatrix(otherSpeciesID), numberOfContractions, numberOfContractions, 0.0_8 )
      end do
      
       call Matrix_constructor( WaveFunction_instance(speciesID)%exchangeCorrelationMatrix, numberOfContractions, numberOfContractions, 0.0_8 )

       call Matrix_constructor( WaveFunction_instance(speciesID)%waveFunctionCoefficients,numberOfContractions, numberOfContractions, 0.0_8 )
       call Vector_constructor( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, int(numberOfContractions) )

       !!cosmo things
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo2, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmoCoupling, numberOfContractions, numberOfContractions, 0.0_8 )

       WaveFunction_instance(speciesID)%exactExchangeFraction = 1.0_8
       WaveFunction_instance(speciesID)%particlesInGrid = 0.0
       WaveFunction_instance(speciesID)%removedOrbitals = 0
       
       !! read the coefficients again from the ".vec" file again
       ! if ( CONTROL_instance%READ_COEFFICIENTS) then

       !    labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
       !    labels(1) = "COEFFICIENTS"

       !    vecUnit=77
          
       !    vecFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
       !    inquire(FILE = vecFile, EXIST = existFile )

       !    if ( existFile) then
       !       open(unit=vecUnit, file=trim(vecFile), status="old", form="formatted")

       !       WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(unit=vecUnit, &
       !            rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.false.,  & 
       !            arguments=labels(1:2))

       !       close(vecUnit)

       !    else 
       !       vecFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
       !       inquire(FILE = vecFile, EXIST = existFile )

       !       if ( existFile) then
       !          open(unit=vecUnit, file=trim(vecFile), status="old", form="unformatted")

       !          WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(unit=vecUnit, &
       !               rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., & 
       !               arguments=labels(1:2))

       !          close(vecUnit)

       !       else
       !          call  Wavefunction_exception( ERROR, "I did not find any .vec coefficients file", "At SCF program, at Wavefunction_constructor")
       !       end if

       !    end if
       ! end if
       
       
    end do
    !!Initialize DFT: Calculate Grids and build functionals
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then

       call system ("lowdin-DFT.x BUILD_SCF_GRID")

       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
          dftUnit = 77
          dftFile = "lowdin."//trim(MolecularSystem_getNameOfSpecie(speciesID))//".grid"
          open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

          labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
          labels(1) = "EXACT-EXCHANGE-FRACTION"

          call Vector_getFromFile(unit=dftUnit, binary=.true., value=WaveFunction_instance(speciesID)%exactExchangeFraction, arguments=labels)
          close(unit=dftUnit)

          ! print *, "el tormento tuyo", speciesID, WaveFunction_instance(speciesID)%exactExchangeFraction
          
       end do

    end if

  end subroutine WaveFunction_constructor

  !>
  !! @brief Contruye la matrix de overlap.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_buildOverlapMatrix(file, speciesID)
    implicit none

    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    arguments(1) = "OVERLAP"
    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))

    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    WaveFunction_instance( speciesID )%overlapMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)
    close(34)

    !! DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
     print *,"Matriz de overlap: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
     call Matrix_show(WaveFunction_instance( speciesID )%overlapMatrix)
    end if

  end subroutine WaveFunction_buildOverlapMatrix

  !>
  !! @brief Contruye la matrix de de transformacion.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_buildTransformationMatrix(file, speciesID, typeOfOrthogonalization )
    implicit none

    character(*), intent(in) :: file
    integer, intent(in) :: speciesID    
    integer, optional, intent(in) :: typeOfOrthogonalization

    type(Matrix) :: eigenVectors
    type(Vector) :: eigenValues
    integer(8) :: numberOfContractions
    integer :: i, j, removed

    !! Numero de contracciones "totales"
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    if ( .not. allocated( WaveFunction_instance( speciesID )%overlapMatrix%values ) ) then

       call WaveFunction_buildOverlapMatrix( trim(file), speciesID )

    end if

    if ( .not. allocated(WaveFunction_instance( speciesID )%transformationMatrix%values)) then

       call Matrix_constructor( WaveFunction_instance( speciesID )%transformationMatrix, &
            int(numberOfContractions,8), int(numberOfContractions,8), Math_NaN )

    end if

    if ( numberOfContractions > 1) then

       call Vector_constructor( eigenValues, int(numberOfContractions) )

       call Matrix_constructor( eigenVectors, numberOfContractions, numberOfContractions)

       !!****************************************************************
       !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
       !!          
       call Matrix_eigen( WaveFunction_instance( speciesID )%overlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )
       
       ! do i = 1 , numberOfContractions
       !   print *, eigenvalues%values(i) 
       ! end do

       do i = 1 , numberOfContractions
          do j = 1 , numberOfContractions
             if ( abs(eigenValues%values(j)) >= CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
                WaveFunction_instance( speciesID )%transformationMatrix%values(i,j) = &
                     eigenVectors%values(i,j)/sqrt( eigenvalues%values(j) )
             else
                WaveFunction_instance( speciesID )%transformationMatrix%values(i,j) = 0
             end if
          end do
       end do

       do i = 1 , numberOfContractions
          if ( abs(eigenValues%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) &
               WaveFunction_instance( speciesID )%removedOrbitals=WaveFunction_instance( speciesID )%removedOrbitals+1
       end do

       if (WaveFunction_instance( speciesID )%removedOrbitals .gt. 0) &
            write(*,"(A,I5,A,A,A,ES9.3)") "Removed ", WaveFunction_instance( speciesID )%removedOrbitals , " orbitals for species ", &
            trim(MolecularSystem_getNameOfSpecie(speciesID)), " with overlap eigen threshold of ", CONTROL_instance%OVERLAP_EIGEN_THRESHOLD

      !!
       !!****************************************************************

       !!****************************************************************
       !! Calcula matriz de transformacion
       !!
       select case (typeOfOrthogonalization)

          !! Ortogonalizacion canonica
       case (CANONICAL_ORTHOGONALIZATION)
          WaveFunction_instance( speciesID )%transformationMatrix%values = &
               WaveFunction_instance( speciesID )%transformationMatrix%values

          !!Ortogonalizacion simetrica
       case (SYMMETRIC_ORTHOGONALIZATION)
          WaveFunction_instance( speciesID )%transformationMatrix%values  = &
               matmul(WaveFunction_instance( speciesID )%transformationMatrix%values, transpose(eigenVectors%values))

       case default
          
          WaveFunction_instance( speciesID )%transformationMatrix%values  = &
               matmul(WaveFunction_instance( speciesID )%transformationMatrix%values, transpose(eigenVectors%values))

       end select

       call Vector_destructor( eigenValues )
       call Matrix_destructor( eigenVectors )
       !!
       !!****************************************************************

    else

       WaveFunction_instance( speciesID )%transformationMatrix%values= 1.0_8

    end if

    !! Debug
    ! print *,"Matriz de transformacion: ", trim(MolecularSystem_instance%species(speciesID)%symbol)
    ! call Matrix_show( WaveFunction_instance( speciesID )%transformationMatrix )

  end subroutine WaveFunction_buildTransformationMatrix

  !>
  !! @brief Contruye la matrix de particula independiente.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_HCoreMatrix(file, speciesID)
    implicit none

    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    integer :: k, l, r, s
    integer :: ParticleID, ParticleID_2
    integer :: contractionID, contractionID_2
    integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    integer :: owner, owner_2
    real(8) :: auxCharge
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))    

    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          

    !! Load Kinetic Matrix
    arguments(1) = "KINETIC"    

    WaveFunction_instance( speciesID )%kineticMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)

    !! Incluiding mass effect       
    if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
       WaveFunction_instance( speciesID )%kineticMatrix%values =  &
            WaveFunction_instance( speciesID )%kineticMatrix%values * &
            ( 1.0_8/MolecularSystem_getMass( speciesID ) -1.0_8 / ParticleManager_getTotalMass() )

    else
       WaveFunction_instance( speciesID )%kineticMatrix%values = &
            WaveFunction_instance( speciesID )%kineticMatrix%values / &
            MolecularSystem_getMass( speciesID )
    end if

    !! Finite Nuclear Mass Correction
    if ( CONTROL_instance%FINITE_MASS_CORRECTION ) then

       k=1
       do particleID = 1, size(MolecularSystem_instance%species(speciesID)%particles)
          do contractionID = 1, size(MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction)

             numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction(contractionID)%numCartesianOrbital
             owner = MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction(contractionID)%owner

             do s = 1, numberOfCartesiansOrbitals
                l=k

                do particleID_2 = 1, size(MolecularSystem_instance%species(speciesID)%particles)
                   do contractionID_2 = 1, size(MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction)

                      numberOfCartesiansOrbitals_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%numCartesianOrbital
                      owner_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%owner

                      do r = 1, numberOfCartesiansOrbitals_2

                         if ( owner .eq. owner_2) then
                            WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)=&
                                 WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)*&
                                 ( 1 + MolecularSystem_getMass( speciesID ) / MolecularSystem_instance%species(speciesID)%particles(particleID)%mass  )

                            WaveFunction_instance( speciesID )%kineticMatrix%values(l,k)=&
                                 WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)
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

    !! DEBUG
    !!   print *,"Matriz de energia cinetica: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !!    call Matrix_show( WaveFunction_instance(speciesID)%kineticMatrix )

    !! Load N-Q- Attraction  Matrix
    arguments(1) = "ATTRACTION"

    WaveFunction_instance( speciesID )%puntualInteractionMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    

    !! Incluiding charge effect
    auxcharge = MolecularSystem_getCharge( speciesID )

    WaveFunction_instance( speciesID )%puntualInteractionMatrix%values = &
         WaveFunction_instance( speciesID )%puntualInteractionMatrix%values * (-auxCharge)

    !! DEBUG
    ! print *,"Matriz de interaccion n-quantum: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%puntualInteractionMatrix )

    !! Build Hcore Matrix
    if ( .not.allocated(WaveFunction_instance( speciesID )%HcoreMatrix%values ) ) then

       call Matrix_constructor( WaveFunction_instance( speciesID )%HcoreMatrix, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )

    end if

    WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%kineticMatrix%values + &
        WaveFunction_instance(speciesID)%puntualInteractionMatrix%values

    !! Add electric field F_i < \mu | e_i | \nu >
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      write (*,"(T2,A15,3F12.8)") "ELECTRIC FIELD:", CONTROL_instance%ELECTRIC_FIELD

      arguments(1) = "MOMENTX"
      WaveFunction_instance(speciesID)%electricField(1) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                           columns=totalNumberOfContractions, &
                                                            unit=unit, binary=.true., arguments=arguments)    
      arguments(1) = "MOMENTY"
      WaveFunction_instance(speciesID)%electricField(2) = Matrix_getFromFile(rows=totalNumberOfContractions, & 
                                                            columns=totalNumberOfContractions, &
                                                            unit=unit, binary=.true., arguments=arguments)    
      arguments(1) = "MOMENTZ"
      WaveFunction_instance(speciesID)%electricField(3) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                            columns=totalNumberOfContractions, &
                                                            unit=unit, binary=.true., arguments=arguments)    

      WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%HCoreMatrix%values + &
                auxcharge * &
        (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(speciesID)%electricField(1)%values + &
         CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(speciesID)%electricField(2)%values + &
         CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(speciesID)%electricField(3)%values )
    end if

    close(34)    
         !WaveFunction_instance(speciesID)%externalPotentialMatrix%values 

    !! DEBUG
    !!   print *,"Matriz de hcore: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !!   call Matrix_show( WaveFunction_instance( speciesID )%HcoreMatrix )

  end subroutine WaveFunction_HCoreMatrix

  !>
  !! @brief Ajusta la matriz de densidad para una especie espcificada
  subroutine WaveFunction_setCoefficientsMatrix(coefficientsMatrix, speciesID )
    implicit none

    type(Matrix), intent(in) :: coefficientsMatrix
    integer :: speciesID

    integer :: totalNumberOfContractions

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

    if( .not. allocated(WaveFunction_instance( speciesID )%waveFunctionCoefficients%values )) then
       call Matrix_constructor( WaveFunction_instance( speciesID )%waveFunctionCoefficients, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
    end if

    call Matrix_copyConstructor( WaveFunction_instance( speciesID )%waveFunctionCoefficients, coefficientsMatrix )

    !! Debug
    ! print*, "Matriz de coefficients inicial ", MolecularSystem_getNameOfSpecie(speciesID)
    ! call Matrix_show(WaveFunction_instance( speciesID )%waveFunctionCoefficients)

  end subroutine WaveFunction_setCoefficientsMatrix



  !>
  !! @brief Calcula las componentes de energia para la especie especificada
  !! @warning Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  subroutine WaveFunction_obtainEnergyComponents(speciesID )
    implicit none

    integer :: speciesID
    integer :: otherSpeciesID
    real(8) :: auxCharge

    auxcharge = MolecularSystem_getCharge( speciesID )

    !! Remove the electric field matrix to calculate the energy components
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%HCoreMatrix%values - &
                auxcharge * &
        (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(speciesID)%electricField(1)%values + &
         CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(speciesID)%electricField(2)%values + &
         CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(speciesID)%electricField(3)%values )
    end if
    
    !! Calcula la energia de dos particulas
    WaveFunction_instance( speciesID )%twoParticlesEnergy = 0.5_8 * &
         sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
         WaveFunction_instance( speciesID )%twoParticlesMatrix%values )

    !! Calcula la energia de repulsion intraespecie
    WaveFunction_instance( speciesID )%hartreeEnergy(speciesID) = 0.5_8 * &
         sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
         WaveFunction_instance( speciesID )%hartreeMatrix(speciesID)%values )

    !! Calcula la energia de intercambio HF
    WaveFunction_instance( speciesID )%exchangeHFEnergy = 0.5_8 * &
         sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
         WaveFunction_instance( speciesID )%exchangeHFMatrix%values )

    !! Calcula energia de particula independiente
    WaveFunction_instance( speciesID )%independentParticleEnergy = &
         sum( transpose( WaveFunction_instance(speciesID)%densityMatrix%values ) * &
         WaveFunction_instance( speciesID )%hcoreMatrix%values )

    !! Calcula energia cinetica para la especie dada
    WaveFunction_instance( speciesID )%kineticEnergy = &
         sum( transpose(WaveFunction_instance( speciesID )%densityMatrix%values) * &
         WaveFunction_instance( speciesID )%kineticMatrix%values )

    !! Calcula energia de potencial externo para la especie dada
    WaveFunction_instance( speciesID )%externalPotentialEnergy = &
         sum( transpose(WaveFunction_instance( speciesID )%densityMatrix%values) * &
         WaveFunction_instance( speciesID )%externalPotentialMatrix%values )

    !! Calcula energia de interaccion entre particulas puntuales y cuanticas
    WaveFunction_instance( speciesID )%puntualInteractionEnergy =  &
         WaveFunction_instance( speciesID )%independentParticleEnergy - &
         WaveFunction_instance( speciesID )%kineticEnergy

    !! Calula enegia de especie independiente (  sin considerar el termino de acoplamiento )
    WaveFunction_instance( speciesID )%independentSpecieEnergy = &
         sum( transpose(WaveFunction_instance( speciesID )%densityMatrix%values) * &
         (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) + &
         0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values + &
         WaveFunction_instance( speciesID )%externalPotentialMatrix%values))

    WaveFunction_instance( speciesID )%independentSpecieEnergy = &
         WaveFunction_instance( speciesID )%independentSpecieEnergy 

    !! Calcula energia de acoplamiento por especies
    do otherSpeciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
       if (speciesID .ne. otherSpeciesID) then
          WaveFunction_instance( speciesID )%hartreeEnergy( otherSpeciesID ) = &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               WaveFunction_instance( speciesID )%hartreeMatrix( otherSpeciesID )%values )
       end if
    end do
    
    !! Calcula energia de acoplamiento en caso de mas de una especie presente
    WaveFunction_instance( speciesID )%couplingEnergy = &
         sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
         WaveFunction_instance( speciesID )%couplingMatrix%values )    

    !! Total energy for species
    WaveFunction_instance( speciesID )%totalEnergyForSpecie = &
         WaveFunction_instance( speciesID )%independentSpecieEnergy +  &
         WaveFunction_instance( speciesID )%couplingEnergy + &
         sum(WaveFunction_instance( speciesID )%exchangeCorrelationEnergy(:))


    !! Calcula la energia COSMO	

    if(CONTROL_instance%COSMO)then

       write(*,*)"COSMO energy contributions"

       write(*,*)"Especie = ",trim(MolecularSystem_instance%species(speciesID)%name)

       WaveFunction_instance( speciesID )%cosmoEnergy =  &
            0.5_8* (sum( transpose( wavefunction_instance( speciesID )%densitymatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo1%values )) +0.5_8 *  &
            (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            WaveFunction_instance( speciesID )%cosmo4%values )) + &
            0.5_8*(sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            WaveFunction_instance( speciesID )%cosmo2%values ) + &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            WaveFunction_instance( speciesID )%cosmoCoupling%values )) 

       write(*,*)"COSMO energy 1 = ",0.5_8*(sum( transpose( wavefunction_instance( speciesID )%densitymatrix%values ) * wavefunction_instance( speciesID )%cosmo1%values )) 
       write(*,*)"COSMO energy 2 = ",0.5_8*(sum( transpose( wavefunction_instance( speciesID )%densitymatrix%values ) * wavefunction_instance( speciesID )%cosmo2%values )) 
       write(*,*)"COSMO energy 4 = ",0.5_8*(sum( transpose( wavefunction_instance( speciesID )%densitymatrix%values ) * wavefunction_instance( speciesID )%cosmo4%values )) 
       write(*,*)"COSMO coupling energy = ",0.5_8*(sum( transpose( wavefunction_instance( speciesID )%densitymatrix%values ) * wavefunction_instance( speciesID )%cosmoCoupling%values )) 
       write(*,*)"COSMO total energy = ",WaveFunction_instance( speciesID )%cosmoEnergy


    end if

    !! Put back the electric field matrix to the Hcore matrix
    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
      WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%HCoreMatrix%values + &
                auxcharge * &
        (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(speciesID)%electricField(1)%values + &
         CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(speciesID)%electricField(2)%values + &
         CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(speciesID)%electricField(3)%values )
    end if




    ! print *, "__________________ ENERGY COMPONENTS _______________________"
    ! print *, "	Specie                       ", MolecularSystem_getNameOfSpecie( speciesID )
    ! print *, "	Total Energy                =", WaveFunction_instance( speciesID )%totalEnergyForSpecie
    ! print *, "	Indepent Specie Energy      =", WaveFunction_instance( speciesID )%independentSpecieEnergy
    ! print *, "	Kinetic Energy              =",WaveFunction_instance( speciesID )%kineticEnergy
    ! print *, "	Puntual Interaction Energy  =",WaveFunction_instance( speciesID )%puntualInteractionEnergy
    ! print *, "	Independent Particle Energy =",WaveFunction_instance( speciesID )%independentParticleEnergy
    ! print *, "	Repultion Energy            =",WaveFunction_instance( speciesID )%twoParticlesEnergy
    ! print *, "	Coupling Energy             =", WaveFunction_instance( speciesID )%couplingEnergy
    ! print *, "____________________________________________________________"
    !
  end subroutine WaveFunction_obtainEnergyComponents

  !!cosmo matrices construction

  subroutine WaveFunction_cosmoHCoreMatrix(file, speciesID)
    implicit none

    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    ! integer :: k, l, r, s
    ! integer :: ParticleID, ParticleID_2
    ! integer :: contractionID, contractionID_2
    ! integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    ! integer :: owner, owner_2
    ! integer :: auxCharge
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    !! Open file
    unit = 44
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))  


    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          


    !! Load electron potential vs clasical charges cosmo matrix
    arguments(1) = "COSMO1"    
    WaveFunction_instance( speciesID )%cosmo1 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)



    ! DEBUG
    ! print *,"Matriz cosmo1: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%cosmo1 )

    !! Load clasical potential vs quantum charges cosmo matrix
    arguments(1) = "COSMO4"


    WaveFunction_instance( speciesID )%cosmo4 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    


    !! DEBUG
    ! print *,"Matriz cosmo 4 ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%cosmo4 )

    close(44)    


  end subroutine WaveFunction_cosmoHcoreMatrix

  !<
  !! @brief Contruye una matriz de interaccion con un potencial externo
  !!
  !! @param nameOfSpecie nombre de la especie seleccionada.
  !>
  subroutine WaveFunction_buildExternalPotentialMatrix( file, speciesID )
    implicit none
    character(30) :: nameOfSpecieSelected
    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(50) :: arguments(2)

    arguments(1) = "EXTERNAL_POTENTIAL"
    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))

    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    WaveFunction_instance(speciesID)%externalPotentialMatrix%values = 0.0_8

    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          
    WaveFunction_instance( speciesID )%externalPotentialMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments(1:2))
    close(34)

      !!if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )  then  !! Numerical integration
      !!if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
      !!  WaveFunction_instance(speciesID)%externalPotentialMatrix = &
      !!    IntegralManager_getNumericalInteractionWithPotentialMatrix( &
      !!    ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )

      !!else     !! From xml file
      !!  WaveFunction_instance(speciesID)%externalPotentialMatrix = &
      !!    IntegralManager_getNumericalPotentialMatrixFromXml( &
      !!    ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )
      !!end if
      !!else    !! Analytical Integration  

      !!end if

    !! Debug
    if (  CONTROL_instance%DEBUG_SCFS) then
      print *,"EXTERNAL POTENTIAL MATRIX FOR: ", arguments(2)
      call Matrix_show(WaveFunction_instance(speciesID)%externalPotentialMatrix)
    end if

  end subroutine WaveFunction_buildExternalPotentialMatrix
  
  
  !>
  !! @brief Builds two-particles matrix.
  subroutine WaveFunction_buildTwoParticlesMatrix( nameOfSpecies, densityMatrixIN, factorIN, twoParticlesMatrixOUT )
    implicit none

    character(*)  :: nameOfSpecies
    type(Matrix), optional :: densityMatrixIN
    real(8), optional :: factorIN
    type(Matrix), optional :: twoParticlesMatrixOUT
    
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
    integer :: speciesID
    integer :: status
    integer :: u, v, i

    type(Matrix) :: densityMatrix
    type(Matrix) :: twoParticlesMatrix

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecies ) )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

    call Matrix_constructor(densityMatrix, int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), 0.0_8 )
    if ( present(densityMatrixIN)) then
       densityMatrix%values = densityMatrixIN%values
    else       
       densityMatrix%values = wavefunction_instance(speciesID)%densityMatrix%values
    end if

    if ( present(factorIN)) then
       factor = MolecularSystem_getFactorOfExchangeIntegrals( speciesID)*factorIN
    else
       factor = MolecularSystem_getFactorOfExchangeIntegrals( speciesID)*WaveFunction_instance(speciesID)%exactExchangeFraction
    end if

    call Matrix_constructor(twoParticlesMatrix, int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), 0.0_8 )
    
    !! This matrix is only calculated if there are more than one particle for speciesID or if the user want to calculate it.
    if ( MolecularSystem_getNumberOfParticles( speciesID ) > 1 .or.  CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE ) then

       if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

          !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, aa, bb, rr, ss, shellIntegrals, i, coulomb, exchange, tmpArray)
          nthreads = OMP_GET_NUM_THREADS()
          threadid =  OMP_GET_THREAD_NUM()
          unitid = 40 + threadid

          write(fileid,*) threadid
          fileid = trim(adjustl(fileid))

          if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(speciesID)%isElectron) then

             open( UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.ints", status='old', access='stream', form='Unformatted')

          else

             open( UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//".ints", status='old', access='stream', form='Unformatted')

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

       else !! Direct

          call DirectIntegralManager_getDirectIntraRepulsionMatrix(&
               speciesID, &
               trim(CONTROL_instance%INTEGRAL_SCHEME), &
               densityMatrix, & 
               tmpTwoParticlesMatrix, &
               factor)

          twoParticlesMatrix%values = tmpTwoParticlesMatrix

          deallocate(tmpTwoParticlesMatrix)

       end if

       if ( .not. InterPotential_instance%isInstanced) then
          twoParticlesMatrix%values = &
              twoParticlesMatrix%values * ( MolecularSystem_getCharge(speciesID=speciesID ) )**2.0_8
       end if

    end if

    if ( .not. ( present(twoParticlesMatrixOUT) )) then
       wavefunction_instance(speciesID)%twoParticlesMatrix%values=twoParticlesMatrix%values
    else
       twoParticlesMatrixOUT%values=twoParticlesMatrix%values
    end if
    
    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "two particle matrix for: ", trim(nameOfSpecies)
       call Matrix_show(wavefunction_instance(speciesID)%twoParticlesMatrix)
    end if

  end subroutine WaveFunction_buildTwoParticlesMatrix

  !>
  !! @brief Builds the coupling matrix.
  subroutine WaveFunction_buildCouplingMatrix( nameOfSpecies, initialSpeciesIterator, densityMatricesIN, couplingMatrixOUT, hartreeMatricesOUT )
    implicit none
    
    character(*), optional :: nameOfSpecies
    integer, optional :: initialSpeciesIterator
    integer :: initialSpeciesIteratorSelected
    type(Matrix), optional :: densityMatricesIN(*)
    type(Matrix), optional :: couplingMatrixOUT
    type(Matrix), optional :: hartreeMatricesOUT(*)


    character(30) :: nameOfSpeciesSelected
    character(30) :: nameOfOtherSpecies
    integer :: numberOfSpecies
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: currentSpeciesID
    integer :: otherSpeciesID
    integer :: speciesIterator
    integer :: ssize
    integer :: i, j, u
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

    nameOfSpeciesSelected = "E-"    

    initialSpeciesIteratorSelected = 1

    if ( present( initialSpeciesIterator ) ) initialSpeciesIteratorSelected= initialSpeciesIterator
    if ( present( nameOfSpecies ) )  nameOfSpeciesSelected= trim( nameOfSpecies )

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    currentSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpeciesSelected )
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpeciesID)

    allocate(densityMatrices(numberOfSpecies))
    allocate(hartreeMatrices(numberOfSpecies))

    !!Initialize matrices
    do speciesIterator = initialSpeciesIteratorSelected, numberOfSpecies
       call Matrix_constructor(densityMatrices(speciesIterator), int(MolecularSystem_getTotalNumberOfContractions(speciesIterator),8), &
            int(MolecularSystem_getTotalNumberOfContractions(speciesIterator),8), 0.0_8 )
       call Matrix_constructor(hartreeMatrices(speciesIterator), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

       if (present(densityMatricesIN) ) then
          densityMatrices(speciesIterator)%values=densityMatricesIN(speciesIterator)%values
       else
          densityMatrices(speciesIterator)%values=wavefunction_instance(speciesIterator)%densityMatrix%values          
       end if

    end do

    call Matrix_constructor(couplingMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

    
    if( numberOfSpecies > 1 ) then
       
       ssize = size(couplingMatrix%values,dim=1)

       if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

          !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, a, b, r, s, integral, u, i, j, coulomb, auxMatrix, speciesIterator, &
          !$OMP& otherSpeciesID, nameofOtherSpecies, otherNumberOfContractions)

          nthreads = OMP_GET_NUM_THREADS()
          threadid = OMP_GET_THREAD_NUM()
          unitid = 40 + threadid

          write(fileid,*) threadid
          fileid = trim(adjustl(fileid))

          allocate(auxMatrix(ssize, ssize))
          auxMatrix=0.0_8                

          do speciesIterator = initialSpeciesIteratorSelected, numberOfSpecies

             otherSpeciesID = speciesIterator
             nameOfOtherSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )          
             OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID)

             !! Restringe suma de terminos repulsivos de la misma especie.
             if ( otherSpeciesID /= currentSpeciesID ) then
                
                ! hartreeMatrices(otherSpeciesID)%values = 0.0_8

                if( currentSpeciesID > otherSpeciesID) then  

                   auxMatrix = 0.0_8

                   
                   !! open file for integrals
                   if(CONTROL_instance%IS_OPEN_SHELL .and. &
                        MolecularSystem_instance%species(currentSpeciesID)%isElectron .and. &
                        MolecularSystem_instance%species(otherSpeciesID)%isElectron ) then
                      open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.E-BETA.ints", &
                           STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                   else if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(otherSpeciesID)%isElectron) then
                      open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA."//trim(nameOfSpecies)//".ints", &
                           STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                   else if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(currentSpeciesID)%isElectron) then
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

                   auxMatrix = auxMatrix * MolecularSystem_getCharge(currentSpeciesID ) * MolecularSystem_getCharge( otherSpeciesID )

                   do i = 1 , ssize
                     do j = i , ssize
                        !$OMP ATOMIC
                       hartreeMatrices(otherSpeciesID)%values(i,j) = &
                       hartreeMatrices(otherSpeciesID)%values(i,j) + auxMatrix(i,j)
                     end do
                   end do

                else

                   auxMatrix=0.0_8

                   !! open file for integrals
                   if(CONTROL_instance%IS_OPEN_SHELL .and. &
                        MolecularSystem_instance%species(currentSpeciesID)%isElectron .and. &
                        MolecularSystem_instance%species(otherSpeciesID)%isElectron ) then
                      open(UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.E-BETA.ints", &
                           STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                   else if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(otherSpeciesID)%isElectron) then
                      open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//".E-ALPHA.ints", &
                           STATUS='OLD', ACCESS='stream', FORM='Unformatted')
                   else if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(currentSpeciesID)%isElectron) then
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

                   auxMatrix = auxMatrix * MolecularSystem_getCharge(currentSpeciesID ) * MolecularSystem_getCharge( otherSpeciesID )

                   do i = 1 , ssize
                     do j = i , ssize
                        !$OMP ATOMIC
                       hartreeMatrices(otherSpeciesID)%values(i,j) = &
                       hartreeMatrices(otherSpeciesID)%values(i,j) + auxMatrix(i,j)
                     end do
                   end do
                   
                end if

             end if

          end do

          deallocate(auxMatrix)

          !$OMP END PARALLEL

       else !! Direct

          do speciesIterator = initialSpeciesIteratorSelected, numberOfSpecies

             otherSpeciesID = speciesIterator
             
             !! Restringe suma de terminos repulsivos de la misma especie.
             if ( otherSpeciesID /= currentSpeciesID ) then

                   call DirectIntegralManager_getDirectInterRepulsionMatrix(&
                        currentSpeciesID, OtherSpeciesID, &
                        trim(CONTROL_instance%INTEGRAL_SCHEME), &
                        densityMatrices(otherSpeciesID), &
                        auxMatrix )

                auxMatrix = auxMatrix * MolecularSystem_getCharge(currentSpeciesID ) * MolecularSystem_getCharge( otherSpeciesID )

                hartreeMatrices(otherSpeciesID)%values = &
                hartreeMatrices(otherSpeciesID)%values + auxMatrix

                !!wavefunction_instance(currentSpeciesID)%couplingMatrix%values = wavefunction_instance(currentSpeciesID)%couplingMatrix%values + auxMatrix
                deallocate(auxMatrix)

             end if
          end do
        
       end if

       ! Symmetric
       do speciesIterator = initialSpeciesIteratorSelected, numberOfSpecies

          otherSpeciesID = speciesIterator
          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpeciesID /= currentSpeciesID ) then

             do i = 1 , ssize
                do j = i , ssize

                   hartreeMatrices(otherSpeciesID)%values(j,i) = &
                        hartreeMatrices(otherSpeciesID)%values(i,j)

                end do
             end do

             nameOfOtherSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )          

             if ( nameOfOtherSpecies .ne. CONTROL_instance%SCF_GHOST_SPECIES ) &
                  couplingMatrix%values = couplingMatrix%values + hartreeMatrices(otherSpeciesID)%values 

          end if
       end do

    end if

    if ( .not. ( present(hartreeMatricesOUT) .or. present(couplingMatrixOUT) ) ) then
       wavefunction_instance(currentSpeciesID)%couplingMatrix%values=couplingMatrix%values
       do otherSpeciesID = initialSpeciesIteratorSelected, numberOfSpecies
          if ( otherSpeciesID /= currentSpeciesID ) &
          wavefunction_instance(currentSpeciesID)%hartreeMatrix(otherSpeciesID)%values=hartreeMatrices(otherSpeciesID)%values
       end do
    else
       couplingMatrixOUT%values=couplingMatrix%values
       do otherSpeciesID = initialSpeciesIteratorSelected, numberOfSpecies
          if ( otherSpeciesID /= currentSpeciesID ) &
          hartreeMatricesOUT(otherSpeciesID)%values=hartreeMatrices(otherSpeciesID)%values
       end do
    end if
    
    
    if (  CONTROL_instance%DEBUG_SCFS) then
       do otherSpeciesID = initialSpeciesIteratorSelected, numberOfSpecies
          if ( otherSpeciesID /= currentSpeciesID ) then
             write(*,*) "Hartree Matrix for: ", trim(nameOfSpeciesSelected), MolecularSystem_getNameOfSpecie( otherSpeciesID )          
             call Matrix_show( wavefunction_instance(currentSpeciesID)%hartreeMatrix(otherSpeciesID) )
          end if
       end do
       write(*,*) "Coupling Matrix: ", trim(nameOfSpeciesSelected)
       call Matrix_show( wavefunction_instance(currentSpeciesID)%couplingMatrix )
    end if

 end subroutine WaveFunction_buildCouplingMatrix

  !>
  !! @brief Builds exchange correlation contributions Matrix for DFT calculations (FELIX)
  subroutine WaveFunction_readExchangeCorrelationMatrix( nameOfSpecies, excFileIN, exchangeCorrelationMatrixOUT, exchangeCorrelationEnergyOUT, particlesInGridOUT  )       
    implicit none
    character(*) :: nameOfSpecies
    character(*), optional :: excFileIN
    type(Matrix), optional :: exchangeCorrelationMatrixOUT
    real(8), optional :: exchangeCorrelationEnergyOUT(*)
    real(8), optional :: particlesInGridOUT

    character(50) :: excFile
    type(Matrix) :: exchangeCorrelationMatrix
    real(8),allocatable :: exchangeCorrelationEnergy(:)
    real(8) :: particlesInGrid
    
    integer :: numberOfContractions, numberOfSpecies
    integer :: speciesID
    integer :: otherSpeciesID
    character(50) :: otherNameOfSpecies
    
    character(50) :: labels(2)
    integer :: excUnit
    
    speciesID = MolecularSystem_getSpecieID( nameOfSpecies )

    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(exchangeCorrelationEnergy(numberOfSpecies))

    !! Open file from dft and read matrices
    excUnit = 79
    if(present(excFileIN) ) then
       excFile = excFileIN
    else
       excFile = "lowdin.densmatrix.exc"
    end if
    
    open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")

    labels(2) = trim(nameOfSpecies)
    labels(1) = "NUMBER-OF-PARTICLES"
    call Vector_getFromFile(unit=excUnit, binary=.true., value=particlesInGrid, arguments= labels(1:2) )

    labels(1) = "EXCHANGE-CORRELATION-MATRIX"
    exchangeCorrelationMatrix=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
         binary=.true., arguments=labels(1:2))

    do otherSpeciesID = speciesID, numberOfSpecies
       otherNameOfSpecies=trim(MolecularSystem_getNameOfSpecie(otherSpeciesID))
       labels(1) = "EXCHANGE-CORRELATION-ENERGY"
       labels(2) = trim(nameOfSpecies)//trim(otherNameOfSpecies)
       call Vector_getFromFile(unit=excUnit, binary=.true., value=exchangeCorrelationEnergy(otherSpeciesID), arguments= labels )
    end do
    
    close(unit=excUnit)

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *, "Exc. Corr. Matrix for species ", nameOfSpecies 
       print *, exchangeCorrelationEnergy
       call Matrix_show(exchangeCorrelationMatrix)
    end if

    if(.not. present(excFileIN)) then
       Wavefunction_instance(speciesID)%particlesInGrid=particlesInGrid
       Wavefunction_instance(speciesID)%exchangeCorrelationMatrix%values=exchangeCorrelationMatrix%values
       do otherSpeciesID = speciesID, numberOfSpecies
          Wavefunction_instance(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)=exchangeCorrelationEnergy(otherSpeciesID)
          Wavefunction_instance(otherSpeciesID)%exchangeCorrelationEnergy(speciesID)=exchangeCorrelationEnergy(otherSpeciesID)
       end do
    else
       if(present(particlesInGridOUT)) particlesInGridOUT=particlesInGrid
       if(present(exchangeCorrelationEnergyOUT)) exchangeCorrelationEnergyOUT(1:numberOfSpecies)=exchangeCorrelationEnergy(1:numberOfSpecies)
       if(present(exchangeCorrelationMatrixOUT)) exchangeCorrelationMatrixOUT=exchangeCorrelationMatrix       
    end if
  end subroutine WaveFunction_readExchangeCorrelationMatrix

  !>
  !! @brief Writes density matrices prior to a call to DFT (FELIX)
  subroutine WaveFunction_writeDensityMatricesToFile( densityFileOUT, densityMatricesIN )       
    implicit none
    character(*) :: densityFileOUT
    type(Matrix), optional :: densityMatricesIN(*)

    type(Matrix), allocatable :: densityMatrices(:)
    
    integer :: numberOfSpecies
    integer :: speciesID
    
    integer :: densUnit
    character(50) :: labels(2)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate(densityMatrices(numberOfSpecies))
    densUnit = 78

    do speciesID=1,numberOfSpecies
       if(present(densityMatricesIN)) then
          call Matrix_copyConstructor(densityMatrices(speciesID),densityMatricesIN(speciesID))
       else
          call Matrix_copyConstructor(densityMatrices(speciesID),WaveFunction_instance(speciesID)%densityMatrix)
       end if
    end do
    
    

    open(unit = densUnit, file=trim(densityFileOUT), status="replace", form="unformatted")
    labels(1) = "DENSITY-MATRIX"
    do speciesID = 1, numberOfSpecies
       labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
       call Matrix_writeToFile(densityMatrices(speciesID), unit=densUnit, binary=.true., arguments = labels )
       call Matrix_destructor(densityMatrices(speciesID))
    end do
    close (densUnit)

    
    deallocate(densityMatrices)
    
    
  end subroutine WaveFunction_writeDensityMatricesToFile

  !>
  !! @brief Builds fock Matrix
  subroutine WaveFunction_buildFockMatrix( nameOfSpecie )       
    implicit none

    character(*), optional :: nameOfSpecie

    character(30) :: nameOfSpecieSelected
    integer :: speciesID
    ! type(Matrix)::cosmoContribution

    nameOfSpecieSelected = "E-"    
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%hcoreMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 1 (hcore): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    !! cosmo fock matrix

    !!full coupling
    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
         0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
         wavefunction_instance(speciesID)%cosmo4%values)+ &
         wavefunction_instance(speciesID)%cosmo2%values + &
         wavefunction_instance(speciesID)%cosmoCoupling%values 

    !!half coupling
    ! wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
    !      0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
    !      wavefunction_instance(speciesID)%cosmo4%values)+ &
    !      wavefunction_instance(speciesID)%cosmo2%values +0.5_8*( &
    !      wavefunction_instance(speciesID)%cosmoCoupling%values) 

    !!without coupling
    ! wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
    !      0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
    !      wavefunction_instance(speciesID)%cosmo4%values)+ &
    !      wavefunction_instance(speciesID)%cosmo2%values


    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + wavefunction_instance(speciesID)%twoParticlesMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 2 (+ two particles): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + wavefunction_instance(speciesID)%couplingMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 3 (+ coupling): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    !!!FELIX, agrega la matriz para hacer calculo DFT
    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + wavefunction_instance(speciesID)%exchangeCorrelationMatrix%values
   
    
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 3.1 (+ exchangeCorrelation): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if
    

    if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
      wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
         wavefunction_instance(speciesID)%externalPotentialMatrix%values
      if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 4 (+ external potential): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)

      end if
    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK: "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

  end subroutine WaveFunction_buildFockMatrix

  !>
  !! @brief Calcula la matriz de densidad para una especie especificada
  subroutine WaveFunction_buildDensityMatrix( nameOfSpecies )
    implicit none

    character(*), optional :: nameOfSpecies

    character(30) :: nameOfSpeciesSelected
    integer :: orderMatrix
    integer :: speciesID
    integer :: ocupationNumber
    integer :: i
    integer :: j
    integer :: k
    character(50) :: densFile, labels(2)
    integer :: densUnit

    nameOfSpeciesSelected = "E-"
    if ( present( nameOfSpecies ) )  nameOfSpeciesSelected= trim( nameOfSpecies )

    speciesID = MolecularSystem_getSpecieID( trim(nameOfSpeciesSelected) )

    orderMatrix = size( wavefunction_instance(speciesID)%densityMatrix%values, DIM = 1 )

    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )

    wavefunction_instance(speciesID)%densityMatrix%values = 0.0_8

    !! Segment for fractional occupations: 1
    if (CONTROL_instance%IONIZE_MO /= 0 .and. trim(nameOfSpeciesSelected) == trim(CONTROL_instance%IONIZE_SPECIE(1)) ) then
       wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO) = &
            wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO)*sqrt(CONTROL_instance%MO_FRACTION_OCCUPATION)

    end if

    do i = 1 , orderMatrix
       do j = 1 , orderMatrix
          do k = 1 , ocupationNumber

             wavefunction_instance(speciesID)%densityMatrix%values(i,j) =  &
                  wavefunction_instance(speciesID)%densityMatrix%values( i,j ) + &
                  ( wavefunction_instance(speciesID)%waveFunctionCoefficients%values(i,k) &
                  * wavefunction_instance(speciesID)%waveFunctionCoefficients%values(j,k) )
          end do
       end do
    end do

    wavefunction_instance(speciesID)%densityMatrix%values =  MolecularSystem_getEta( speciesID )  * wavefunction_instance(speciesID)%densityMatrix%values

    !! Segment for fractional occupations: 1
    if (CONTROL_instance%IONIZE_MO /= 0 .and. trim(nameOfSpeciesSelected) == trim(CONTROL_instance%IONIZE_SPECIE(1))) then
       wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO) = &
            wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO)/sqrt(CONTROL_instance%MO_FRACTION_OCCUPATION)
    end if

    !!DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Density Matrix ", trim(nameOfSpeciesSelected)
       call Matrix_show(wavefunction_instance(speciesID)%densityMatrix)
    end if

  end subroutine WaveFunction_buildDensityMatrix

  !>
  !! @brief Calculates total energy for one species
  subroutine WaveFunction_obtainTotalEnergyForSpecie( nameOfSpecie )

    implicit none

    character(*), optional :: nameOfSpecie

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    wavefunction_instance(speciesID)%totalEnergyForSpecie = &
         sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
         *  (( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
         + 0.5_8 *wavefunction_instance(speciesID)%twoParticlesMatrix%values &
         + wavefunction_instance(speciesID)%couplingMatrix%values)) &
         + sum(wavefunction_instance(speciesID)%exchangeCorrelationEnergy(:))

    if(CONTROL_instance%COSMO)then

       wavefunction_instance( speciesID )%totalEnergyForSpecie =wavefunction_instance( speciesID )%totalEnergyForSpecie + 0.5_8 * &
            (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo1%values )+ &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo2%values ) + &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo4%values ) + &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmoCoupling%values ))
    end if

    if( allocated(wavefunction_instance(speciesID)%externalPotentialMatrix%values) ) then

       wavefunction_instance( speciesID )%totalEnergyForSpecie =wavefunction_instance( speciesID )%totalEnergyForSpecie + &
       sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
            *  wavefunction_instance(speciesID)%externalPotentialMatrix%values)
    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Total energy for "// trim(nameOfSpecieSelected) //"= ", wavefunction_instance(speciesID)%totalEnergyForSpecie
    end if

  end subroutine WaveFunction_obtainTotalEnergyForSpecie

  !>
  !! @brief Calcula la energia total para el sistema estudiado  
  subroutine WaveFunction_obtainTotalEnergy( totalEnergy, totalCouplingEnergy,  electronicRepulsionEnergy, cosmo3Energy)
    implicit none

    real(8) :: totalEnergy
    real(8) :: totalCouplingEnergy
    real(8) :: electronicRepulsionEnergy

    ! character(30) :: nameOfSpecieSelected
    integer :: speciesID
    integer :: otherSpeciesID

    !! cosmo

    type(surfaceSegment) :: surface_aux2
    real(8) :: cosmo3Energy

    totalEnergy = 0.0_8
    totalCouplingEnergy = 0.0_8
    cosmo3Energy = 0.0_8



    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       !! Calula enegia de especie independiente ( sin considerar el termino de acoplamiento )
       if( .not. allocated( WaveFunction_instance( speciesID )%externalPotentialMatrix%values ) ) then

          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values))

       else if(CONTROL_instance%COSMO)then


          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values))


          wavefunction_instance( speciesID )%independentSpecieEnergy =wavefunction_instance( speciesID )%independentSpecieEnergy + 0.5_8 * &
               (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo1%values )+ &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo2%values ) +  &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo4%values ) + &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmoCoupling%values))

          ! wavefunction_instance( speciesID )%independentSpecieEnergy =wavefunction_instance( speciesID )%independentSpecieEnergy + 0.5_8 * &
          !      (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo1%values )+ &
          !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo2%values ) +  &
          !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo4%values)) 


       else

          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values &
               + WaveFunction_instance( speciesID )%externalPotentialMatrix%values))

       end if

       WaveFunction_instance( speciesID )%independentSpecieEnergy = &
            WaveFunction_instance( speciesID )%independentSpecieEnergy + &
            WaveFunction_instance( speciesID )%exchangeCorrelationEnergy(speciesID)

       totalEnergy = totalEnergy + WaveFunction_instance( speciesID )%independentSpecieEnergy

    end do

    !! Adicionado energia de interaccion entre particulas puntuales
    totalEnergy = totalEnergy + MolecularSystem_getPointChargesEnergy()


    !! cosmo potential nuclei-charges nuclei

    if(CONTROL_instance%COSMO)then
       call CosmoCore_lines(surface_aux2)
       call CosmoCore_filler(surface_aux2)

       call CosmoCore_nucleiPotentialNucleiCharges(surface_aux2,cosmo3Energy)
       totalEnergy=totalEnergy+cosmo3Energy

    end if

    !! Adicionar  energia de acoplamiento y recalcula matrices de acoplamiento, including E-ALPHA/E-BETA

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
      do otherSpeciesID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies()

          totalCouplingEnergy = totalCouplingEnergy + (sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
              * (wavefunction_instance(speciesID)%hartreeMatrix(otherSpeciesID)%values))) 

          totalEnergy = totalEnergy+wavefunction_instance(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)

      end do

    end do

    !! Adds inter-electron species coupling energy
    ! electronicRepulsionEnergy = WaveFunction_getAlphaBetaRepulsion()

    !! Remove the UHF "repulsion" energy
    ! totalCouplingEnergy = totalCouplingEnergy - electronicRepulsionEnergy

    !! Total Energy
    totalEnergy = totalEnergy +  totalCouplingEnergy + electronicRepulsionEnergy 

  end subroutine WaveFunction_obtainTotalEnergy


  !>
  !! @brief calcula la energia total de acoplamiento para una especie especificada
  function WaveFunction_getAlphaBetaRepulsion() result( output )
    implicit none

    real(8) :: output

    character(30) :: nameOfSpecie
    character(30) :: nameOfOtherSpecie
    real(8) :: auxValue
    ! real(8) :: auxRepulsion
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    integer :: numberOfTotalContractions
    integer :: numberOfTotalContractionsOfOtherSpecie
    integer :: speciesID
    integer :: otherSpeciesID
    ! integer :: outFile
    integer :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: u
    integer :: m
    real(8), allocatable, target :: auxMatrix(:,:)
!    real(8), allocatable :: auxMatrix2(:,:)
    ! integer :: arrayNumber

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: status

    output =0.0_8

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
       do otherSpeciesID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies()

          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpeciesID /= speciesID ) then

             nameOfSpecie = MolecularSystem_getNameOfSpecie( speciesID )
             numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
             numberOfTotalContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

             nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID )
             numberOfContractionsOfOtherSpecie = MolecularSystem_getNumberOfContractions( otherSpeciesID )
             numberOfTotalContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( otherSpeciesID )

             !Restringe la suma a solo electrones
             if(trim(nameOfSpecie)=="E-ALPHA" .and. trim(nameOfOtherSpecie)=="E-BETA") then

              if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then
               
               !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, auxValue, m, a, b, r, s, integral, u)

               nthreads = OMP_GET_NUM_THREADS()
               threadid =  OMP_GET_THREAD_NUM()
               unitid = 40 + threadid

               write(fileid,*) threadid
               fileid = trim(adjustl(fileid))

               !! open file for integrals
               open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
                    STATUS='OLD', ACCESS='stream', FORM='Unformatted')

               auxValue = 0.0_8
               m = 0

               readIntegrals : do

                  read(unitid, iostat=status) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                       r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                       integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                  if(status == -1 ) then
                     print*, "end of file! file: ",trim(fileid)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints"
                     exit readIntegrals
                  end if

                  do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                     if (a(u) == -1) exit readIntegrals

                     m = m + 1
                     !print *, integral(u)


                     auxValue = auxValue +&
                          (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                          * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
                          *  integral(u))

                     if(b(u) /= a(u)) then

                        m = m + 1

                        auxValue = auxValue +&
                             (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                             * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
                             *  integral(u))
                     end if

                     if(s(u) /= r(u)) then

                        m = m + 1

                        auxValue = auxValue +&
                             (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                             * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
                             *  integral(u))
                     end if

                     if(b(u) /= a(u) .and. s(u) /= r(u)) then

                        m = m + 1

                        auxValue = auxValue +&
                             (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                             * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
                             *  integral(u))
                     end if

                  end do

               end do readIntegrals

               auxValue = auxValue *  MolecularSystem_getCharge( speciesID=speciesID ) &
                    * MolecularSystem_getCharge( speciesID=otherSpeciesID )

               !$OMP ATOMIC
               output = output + auxValue

               close(unitid)                

               !$OMP END PARALLEL
              else !! DIRECT

               ! if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
               !  call WaveFunction_exception(ERROR, "Direct integrals are not implemented in DFT yet", "trololo")
               ! end if
          
               call DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix(&
                      speciesID, OtherSpeciesID, &
                      trim(CONTROL_instance%INTEGRAL_SCHEME), &
                      wavefunction_instance(speciesID)%densityMatrix, &
                      wavefunction_instance(otherSpeciesID)%densityMatrix, &
                      auxMatrix )

              
               auxMatrix = auxMatrix * MolecularSystem_getCharge(speciesID ) * MolecularSystem_getCharge( otherSpeciesID )
               output = output + (sum( (auxMatrix))) 

               deallocate(auxMatrix)

              end if !! 

                output = output + (sum( transpose ( wavefunction_instance(speciesID)%densityMatrix%values ) * &
                                WaveFunction_instance(speciesID)%hartreeMatrix(otherSpeciesID)%values ))

             end if
          end if


       end do
    end do

  end function WaveFunction_getAlphaBetaRepulsion

  !**
  ! @brief Retorna la matriz  de coeficientes de combinacion
  !
  !**
  function WaveFunction_getWaveFunctionCoefficients( nameOfSpecie ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    type(Matrix) ::  output

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    if ( allocated( WaveFunction_instance(speciesID)%waveFunctionCoefficients%values) ) then

       call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%waveFunctionCoefficients )

    end if

  end function WaveFunction_getWaveFunctionCoefficients

  !**
  ! @brief Retorna valores propios del sistema molecular
  !
  !**
  function WaveFunction_getMolecularOrbitalsEnergy( nameOfSpecie ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    type(Vector) ::  output

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    if ( allocated( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values) ) then

       call Vector_copyConstructor( output, WaveFunction_instance(speciesID)%molecularOrbitalsEnergy )

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

  subroutine WaveFunction_buildCosmo2Matrix(nameOfSpecie)

    character(*), optional :: nameOfSpecie
    type(species) :: specieSelected
    character(30) :: nameOfSpecieSelected
    integer :: speciesID


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


    nameOfSpecieSelected = "E-"
    if (present(nameOfSpecie))  nameOfSpecieSelected= trim(nameOfSpecie)
    speciesID = MolecularSystem_getSpecieID(nameOfSpecie=trim(nameOfSpecieSelected))
    specieSelected=MolecularSystem_instance%species(speciesID)


    open(unit=110, file=trim(nameOfSpecieSelected)//"_qq.inn", status='old', form="unformatted")
    read(110)m


    ! if(allocated(cosmo_int)) deallocate(cosmo_int)
    ! allocate(cosmo_int(m))

    ! close(unit=110)


    if(allocated(labels)) deallocate(labels)
    allocate(labels(MolecularSystem_instance%species(speciesID)%basisSetSize))

    if(allocated(ints_mat_aux)) deallocate(ints_mat_aux)
    allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID)))


    if(allocated(cosmo2_aux)) deallocate(cosmo2_aux)
    allocate(cosmo2_aux(MolecularSystem_getTotalNumberOfContractions(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID)))


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


    ! call Matrix_show(wavefunction_instance(speciesID)%densityMatrix)

    m = 0

    ii = 0
    do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
       do h = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)

          hh = h
          ii = ii + 1
          jj = ii - 1

          do i = g, size(MolecularSystem_instance%species(speciesID)%particles)
             do j = hh, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !!saving integrals on Matrix
                do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      iii=0
                      do gg = 1, size(MolecularSystem_instance%species(speciesID)%particles)
                         do ll = 1, size(MolecularSystem_instance%species(speciesID)%particles(gg)%basis%contraction)

                            hhh = ll
                            iii = iii + 1
                            jjj = iii - 1
                            do p = gg, size(MolecularSystem_instance%species(speciesID)%particles)
                               do o = hhh, size(MolecularSystem_instance%species(speciesID)%particles(p)%basis%contraction)
                                  jjj = jjj + 1

                                  !!saving integrals on Matrix
                                  do pp = labels(iii), labels(iii) + (MolecularSystem_instance%species(speciesID)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                     do oo = labels(jjj), labels(jjj) + (MolecularSystem_instance%species(speciesID)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
                                        m = m + 1

                                        read(110)cosmo_int

                                        ints_mat_aux(pp, oo) =(wavefunction_instance(speciesID)%densityMatrix%values(pp,oo))* cosmo_int
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
                            wavefunction_instance(speciesID)%cosmo2%values(k,l)=cosmo2_aux(k,l)
                            wavefunction_instance(speciesID)%cosmo2%values(l,k)=wavefunction_instance(speciesID)%cosmo2%values(k,l)
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
       write(*,*) "COSMO 2 matrix for: ", trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%cosmo2)
    end if



  end subroutine WaveFunction_buildCosmo2Matrix


  subroutine WaveFunction_buildCosmoCoupling(nameOfSpecie)

    character(*), optional :: nameOfSpecie
    type(species) :: specieSelected
    type(species) :: otherSpecieSelected
    character(30) :: nameOfSpecieSelected
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


    nameOfSpecieSelected = "E-"
    if (present(nameOfSpecie))  nameOfSpecieSelected= trim(nameOfSpecie)

    currentSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpeciesID)
    specieSelected=MolecularSystem_instance%species(currentSpeciesID)

    if(allocated(labels)) deallocate(labels)
    allocate(labels(MolecularSystem_instance%species(currentSpeciesID)%basisSetSize))

    wavefunction_instance(currentSpeciesID)%cosmoCoupling%values(:,:)=0.0_8

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


    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

       wavefunction_instance(currentSpeciesID)%cosmoCoupling%values = 0.0_8


       do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies()

          otherSpeciesID = speciesIterator
          nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID )          
          OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID)
          otherSpecieSelected=MolecularSystem_instance%species(otherSpeciesID)

          if ( otherSpeciesID /= currentSpeciesID ) then

             ! write(*,*)"hola other and current", otherSpeciesID,currentSpeciesID 

             !      wavefunction_instance(currentSpeciesID)%cosmoCoupling%values = 0.0_8

             open(unit=110, file=trim(nameOfOtherSpecie)//trim(nameOfSpecieSelected)//"_qq.cup", status='old', form="unformatted")
             ! open(unit=110, file=trim(nameOfSpecieSelected)//trim(nameOfOtherSpecie)//"_qq.cup", status='old', form="unformatted")
             read(110)m

             ! if(allocated(cosmo_int)) deallocate(cosmo_int)
             ! allocate(cosmo_int(m))


             if(allocated(otherLabels)) deallocate(otherLabels)
             allocate(otherLabels(MolecularSystem_instance%species(otherSpeciesID)%basisSetSize))

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
             allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(otherSpeciesID), MolecularSystem_getTotalNumberOfContractions(otherSpeciesID)))

             ints_mat_aux=0.0_8                


             if(allocated(cosmoCoup_aux)) deallocate(cosmoCoup_aux)
             allocate(cosmoCoup_aux(MolecularSystem_getTotalNumberOfContractions(currentSpeciesID), MolecularSystem_getTotalNumberOfContractions(currentSpeciesID)))


             m = 0

             ii = 0
             do g = 1, size(MolecularSystem_instance%species(currentSpeciesID)%particles)
                do h = 1, size(MolecularSystem_instance%species(currentSpeciesID)%particles(g)%basis%contraction)

                   hh = h
                   ii = ii + 1
                   jj = ii - 1

                   do i = g, size(MolecularSystem_instance%species(currentSpeciesID)%particles)
                      do j = hh, size(MolecularSystem_instance%species(currentSpeciesID)%particles(i)%basis%contraction)

                         jj = jj + 1

                         !!saving integrals on Matrix
                         do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(currentSpeciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                            do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(currentSpeciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                               iii=0
                               do gg = 1, size(MolecularSystem_instance%species(otherSpeciesID)%particles)
                                  do ll = 1, size(MolecularSystem_instance%species(otherSpeciesID)%particles(gg)%basis%contraction)

                                     hhh = ll
                                     iii = iii + 1
                                     jjj = iii - 1

                                     do p = gg, size(MolecularSystem_instance%species(otherSpeciesID)%particles)
                                        do o = hhh, size(MolecularSystem_instance%species(otherSpeciesID)%particles(p)%basis%contraction)
                                           jjj = jjj + 1

                                           !!saving integrals on Matrix
                                           do pp = otherlabels(iii), otherlabels(iii) + (MolecularSystem_instance%species(otherSpeciesID)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                              do oo = otherlabels(jjj), otherlabels(jjj) + (MolecularSystem_instance%species(otherSpeciesID)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
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

                write(*,*)"cosmo Coupling = "//trim(nameofSpecieSelected)

                call Matrix_show(wavefunction_instance(currentSpeciesID)%cosmoCoupling)

                write(*,*)"cosmo density matrix used = "//trim(nameOfOtherSpecie)

                call Matrix_show(wavefunction_instance(otherSpeciesID)%densityMatrix)

             end if
             close(unit=110)
          end if
       end do


    end if



  end subroutine WaveFunction_buildCosmoCoupling

  subroutine WaveFunction_cosmoQuantumCharge()

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

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

    do f = 1, numberOfSpecies

       specieSelected=MolecularSystem_instance%species(f)

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))

       orderOfMatrix = MolecularSystem_getTotalNumberOfContractions(f)

       arguments(2) = MolecularSystem_getNameOfSpecie(f)

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

       charges_file="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".charges"
       open(unit=100, file=trim(charges_file), status='old', form="unformatted")
       read(100)m

       if(allocated(qiDensityCosmo)) deallocate(qiDensityCosmo)
       allocate(qiDensityCosmo(orderOfMatrix, orderOfMatrix,numberOfPointCharges))
       ii = 0
       do g = 1, size(MolecularSystem_instance%species(f)%particles)
          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)
             hh = h
             ii = ii + 1
             jj = ii - 1
             do i = g, size(MolecularSystem_instance%species(f)%particles)
                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)
                   jj = jj + 1
                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
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
       write(*,*) "COSMO Charges for ",MolecularSystem_getNameOfSpecie( f )," = ", sum(qiCosmo(:))

    end do


    charges_file="qTotalCosmo.charges"
    open(unit=100, file=trim(charges_file), status='replace', form="unformatted")
    write(100) qTotalCosmo(:)
    close(100)
    ! Debug
    ! write(*,*)"qTotalCosmo ", qTotalCosmo(:)
    ! write(*,*)"sumqTotalCosmo ", sum(qTotalCosmo(:))

  end subroutine WaveFunction_cosmoQuantumCharge



end module WaveFunction_



  !   !>
  !   !! @brief   Contruye la matrix de acoplamineto para la especie especificada (C)
  !   !!      la matrix resultante sera tenida en cuenta en la construccion de la matriz de Fock
  !   !!
  !   !>
  !   subroutine WaveFunction_buildCouplingMatrixElectronFree( nameOfSpecie, output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(matrix) :: output

  !     character(30) :: nameOfSpecieSelected
  !     character(30) :: nameOfOtherSpecie
  !     integer :: numberOfContractions
  !     integer :: otherNumberOfContractions
  !     integer :: currentSpeciesID
  !     integer :: otherSpeciesID
  !     integer :: speciesIterator
  !     integer(8) :: ssize
  !     integer :: i, j
  !     integer :: k, l
  !     integer :: arrayNumber
  !     real(8), allocatable :: auxMatrix(:,:)
  !     real(8) :: coulomb

  !     integer :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer :: u

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     currentSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
  !     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpeciesID)

  !     if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !        ssize = size(WaveFunction_instance( currentSpeciesID )%wavefunction_instance(speciesID)%couplingMatrix%values,dim=1)
  !        allocate(auxMatrix(ssize, ssize))
  !        call Matrix_constructor(output, ssize, ssize)
  !        output%values=0.0_8

  !        do speciesIterator = MolecularSystem_beginSpecie(), MolecularSystem_endSpecie()

  !           otherSpeciesID = MolecularSystem_getSpecieID( iteratorOfSpecie = speciesIterator )
  !           nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID )

  !           if(trim(nameOfSpecieSelected) == "e-ALPHA" .and. trim(nameOfOtherSpecie) == "e-BETA" ) cycle
  !           if(trim(nameOfSpecieSelected) == "e-BETA" .and. trim(nameOfOtherSpecie) == "e-ALPHA" ) cycle

  !           OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID)

  !           !! Restringe suma de terminos repulsivos de la misma especie.
  !           if ( otherSpeciesID /= currentSpeciesID ) then

  !              !! ALL IN MEMORY
  !              if( .not. IntegralManager_instance%toDisk ) then

  !                 if( currentSpeciesID > otherSpeciesID) then

  !                    call IntegralManager_interspecieRepulsionIntegral (otherSpeciesID, currentSpeciesID , isInterSpecies=.true., arrayNumber=arrayNumber)

  !                    auxMatrix=0.0_8
  !                    u = 0

  !                    do i=1,otherNumberOfContractions
  !                       do j=i, otherNumberOfContractions
  !                          do k=1,numberOfContractions
  !                             do l=k,numberOfContractions

  !                                u = u + 1

  !                                coulomb = WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values(i,j) * &
  !                                     IntegralManager_instance%interSpecieRepulsionIntegrals(arrayNumber)%values(u)

  !                                auxMatrix(k,l) = auxMatrix(k,l) + coulomb

  !                                if( i /= j ) auxMatrix(k,l) = auxMatrix(k,l) + coulomb

  !                             end do
  !                          end do
  !                       end do
  !                    end do

  !                 else

  !                    auxMatrix=0.0_8
  !                    u = 0

  !                    call IntegralManager_interspecieRepulsionIntegral (currentSpeciesID, otherSpeciesID, isInterSpecies=.true., arrayNumber=arrayNumber)

  !                    do i=1,numberOfContractions
  !                       do j=i, numberOfContractions
  !                          do k=1,otherNumberOfContractions
  !                             do l=k,otherNumberOfContractions

  !                                u = u + 1

  !                                coulomb = WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values(k,l)* &
  !                                     IntegralManager_instance%interSpecieRepulsionIntegrals(arrayNumber)%values(u)

  !                                auxMatrix(i,j) = auxMatrix(i,j) + coulomb

  !                                if( k /= l ) auxMatrix(i,j) = auxMatrix(i,j) + coulomb

  !                             end do
  !                          end do
  !                       end do
  !                    end do

  !                 end if

  !                 !!ALL IN DISK
  !              else

  !                 if( currentSpeciesID > otherSpeciesID) then

  !                    call IntegralManager_interspecieRepulsionIntegral (otherSpeciesID, currentSpeciesID , isInterSpecies=.true.)

  !                    auxMatrix=0.0_8

  !                    !! open file for integrals
  !                    open(UNIT=34,FILE=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)//".ints", &
  !                         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

  !                    do
  !                       read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !                       do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

  !                          if (a(u) == -1) goto 30

  !                          coulomb = WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values(a(u),b(u))*integral(u)

  !                          auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

  !                          if( a(u) /= b(u) ) auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

  !                       end do

  !                    end do

  ! 30                 continue

  !                    close(34)

  !                 else

  !                    call IntegralManager_interspecieRepulsionIntegral (currentSpeciesID, otherSpeciesID, isInterSpecies=.true.)

  !                    auxMatrix=0.0_8

  !                    !! open file for integrals
  !                    open(UNIT=34,FILE=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
  !                         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

  !                    do
  !                       read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !                       do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

  !                          if (a(u) == -1) goto 20

  !                          coulomb = WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values(r(u),s(u))*integral(u)

  !                          auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

  !                          if( r(u) /= s(u) ) auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

  !                       end do

  !                    end do

  ! 20                 continue

  !                    close(34)

  !                 end if

  !              end if

  !              auxMatrix= auxMatrix * MolecularSystem_getCharge( speciesID=currentSpeciesID ) &
  !                   * MolecularSystem_getCharge( speciesID=otherSpeciesID )

  !              output%values= &
  !                   output%values + auxMatrix



  !           end if

  !        end do

  !        deallocate(auxMatrix)

  !        !! Simetriza la matriz de Acoplamineto
  !        do i = 1 , ssize
  !           do j = i , ssize
  !              output%values(j,i) = &
  !                   output%values(i,j)
  !           end do
  !        end do

  !        WaveFunction_instance( currentSpeciesID )%addCouplingMatrix =.true.
  !        WaveFunction_instance( currentSpeciesID )%wasBuiltFockMatrix = .false.

  !     end if

  !     !       print *,"Matriz de acoplamiento: ", trim(nameOfSpecieSelected)
  !     !       call Matrix_show( output )

  !   end subroutine WaveFunction_buildCouplingMatrixElectronFree

  !   !! Add nuclear-electron correlation with ADFT (this could become useful)
  !   subroutine WaveFunction_buildInterParticleCorrMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     character(30) :: nameOfOtherSpecie
  !     integer :: numberOfContractions
  !     integer :: otherNumberOfContractions
  !     integer :: currentSpeciesID
  !     integer :: otherSpeciesID
  !     integer :: speciesIterator
  !     real(8) :: coulomb
  !     real(8) :: startTime, endTime
  !     real(8) :: correlationEnergy

  !     call cpu_time(startTime)

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     currentSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpeciesID)

  !     WaveFunction_instance( currentSpeciesID )%addInterParticleCorrMatrix =.false.

  !     if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !        otherSpeciesID = MolecularSystem_getSpecieID( nameOfSpecie = "e-" ) ! Only for electron and nuclei
  !        nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID )
  !        OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID)

  !        !! Restringe suma de terminos repulsivos de la misma especie.
  !        if ( otherSpeciesID /= currentSpeciesID ) then
  !           if ( nameofOtherSpecie .ne. nameOfSpecieSelected ) then

  !              if ( .not.allocated(WaveFunction_instance( otherSpeciesID)%wavefunction_instance(speciesID)%densityMatrix%values) ) then


  !                 call WaveFunction_exception(ERROR, "Class object WaveFunction_RHF in the builtCouplingMatrix(" &
  !                      // trim(nameOfOtherSpecie) //") function", &
  !                      "Density matrix for "// trim(nameOfOtherSpecie) //" specie, hasn't been defined." )

  !              end if

  !              if ( CONTROL_instance%CALL_DFT ) then

  !                 call bld_aux_ks_c_mat(6, &
  !                      WaveFunction_instance( otherSpeciesID )%bridge%system, &
  !                      WaveFunction_instance( otherSpeciesID )%bridge%system, &
  !                      WaveFunction_instance( currentSpeciesID )%bridge%system, &
  !                      WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( currentSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      WaveFunction_instance( currentSpeciesID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      otherNumberOfContractions,NumberOfContractions, &
  !                      correlationEnergy,.false.)

  ! 		WaveFunction_instance(currentSpeciesID)%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy=correlationEnergy
  ! 		      WaveFunction_instance( otherSpeciesID )%addInterParticleCorrMatrix =.true.
  ! 		      WaveFunction_instance( currentSpeciesID )%addInterParticleCorrMatrix =.true.
  ! 		      WaveFunction_instance( currentSpeciesID )%wasBuiltFockMatrix = .false.

  !              else 
  !                 WaveFunction_instance(currentSpeciesID)%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy= 0.0
  !              end if
  !           end if
  !        end if


  !     end if

  !   end subroutine WaveFunction_buildInterParticleCorrMatrix

  !   !<
  !   !! @brief Contruye una matriz de interaccion con un potencial externo
  !   !!
  !   !! @param nameOfSpecie nombre de la especie seleccionada.
  !   !>
  !   subroutine WaveFunction_buildExternalPotentialMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( nameOfspecie /= "e-BETA" ) then

  ! 	if( WaveFunction_instance(speciesID)%isThereExternalPotential ) then
  ! 		WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix%values = 0.0_8


  ! 	if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )	then	!! Numerical integration
  ! 		if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
  ! 			WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 				IntegralManager_getNumericalInteractionWithPotentialMatrix( &
  ! 				ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )

  ! 		else 		!! From xml file
  ! 			WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 				IntegralManager_getNumericalPotentialMatrixFromXml( &
  ! 				ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )
  ! 		end if
  ! 	else		!! Analytical Integration	

  ! 		WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 		IntegralManager_getInteractionWithPotentialMatrix( &
  ! 		ExternalPotential_Manager_instance%externalsPots, speciesID, "external" )

  !           end if

  !        end if

  !    else !! Use the same matrix for e-beta and e-alpha

  !        WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  !             WaveFunction_instance( MolecularSystem_getSpecieID( nameOfSpecie="e-ALPHA" ))%wavefunction_instance(speciesID)%externalPotentialMatrix

  !     end if

  !     			print *,"EXTERNAL POTENTIAL MATRIX FOR: ", nameOfSpecie
  !     			call Matrix_show(WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix)

  !   end subroutine WaveFunction_buildExternalPotentialMatrix





  !   subroutine WaveFunction_buildPuntualParticleMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if(.not. associated(WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr)) then
  !        call IntegralManager_buildMatrix( ATTRACTION_INTEGRALS, trim(nameOfSpecieSelected  ) )

  !        WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr => 	&
  !             IntegralManager_getMatrixPtr(ATTRACTION_INTEGRALS, nameOfSpecieSelected)

  !     end if

  !     WaveFunction_instance(speciesID)%puntualParticleMatrix%values = &
  !          WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr%values

  !   end subroutine WaveFunction_buildPuntualParticleMatrix





  !   !**
  !   ! @brief Calcula las componentes de energia para la especie especificada
  !   !
  !   ! @warning 	Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  !   !			obtainTotalEnergyForSpecie
  !   !**
  !   subroutine WaveFunction_obtainEnergyComponents( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID
  !     integer :: otherSpeciesID
  !     type(Matrix) :: auxMatrix


  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     if( trim(nameOfSpecieSelected) == "e-ALPHA") then

  !        otherSpeciesID = MolecularSystem_getSpecieID(nameOfSpecie="e-BETA")

  !     else if (trim(nameOfSpecieSelected) == "e-BETA") then

  !        otherSpeciesID = MolecularSystem_getSpecieID(nameOfSpecie="e-ALPHA")

  !     end if

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     ! 		if( nameOfSpecieSelected == "e-") then
  !     ! 		    auxMatrix=IntegralManager_getInteractionBetweenQuantumClassicMtx(speciesID,3)
  !     ! 		    print *,"ATRACION H-e:",sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values )&
  !     ! 			 * auxMatrix%values )
  !     ! 		end if

  !     !! Calcula energia de repulsion para la especie dada                
  !     !                if ( trim(nameOfSpecieSelected) == "e-ALPHA" .or. trim(nameOfSpecieSelected) == "e-BETA") then
  !     !
  !     !                   WaveFunction_instance( speciesID )%twoParticlesEnergy= (0.5_8 * sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !     !                        * wavefunction_instance(speciesID)%twoParticlesMatrix%values )) + &
  !     !                        (0.5_8 * sum( transpose( WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%densityMatrix%values ) &
  !     !                        * WaveFunction_instance( otherSpeciesID )%wavefunction_instance(speciesID)%twoParticlesMatrix%values ))

  !     !                else 

  !     WaveFunction_instance( speciesID )%twoParticlesEnergy= 0.5_8 * sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !          * wavefunction_instance(speciesID)%twoParticlesMatrix%values )
  !     !                end if

  !     !! Calcula energia de particula independiente
  !     WaveFunction_instance( speciesID )%independentParticleEnergy = sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values )&
  !          * WaveFunction_instance( speciesID )%hcoreMatrix%values )

  !     !! Calcula energia cinetica para la especie dada
  !     WaveFunction_instance( speciesID )%kineticEnergy = sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values)  &
  !          * WaveFunction_instance( speciesID )%kineticMatrixValuesPtr%values )

  !     !! Calcula energia de potencial externo para la especie dada
  !     if( WaveFunction_instance( speciesID )%isThereExternalPotential ) &
  !          WaveFunction_instance( speciesID )%externalPotentialEnergy = sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values)  &
  !          * WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%externalPotentialMatrix%values )

  !     !! Calcula energia de interaccion entre particulas puntuales y cuanticas
  !     WaveFunction_instance( speciesID )%puntualInteractionEnergy =  WaveFunction_instance( speciesID )%independentParticleEnergy &
  !          - WaveFunction_instance( speciesID )%kineticEnergy

  !     !! Calula enegia de especie independiente (  sin considerar el termino de acoplamiento)
  !     if( .not.WaveFunction_instance( speciesID )%isThereExternalPotential ) then

  !        WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !             sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !             *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
  !             + 0.5_8 * wavefunction_instance(speciesID)%twoParticlesMatrix%values))
  !     else

  !        WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !             sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !             *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
  !             + 0.5_8 * wavefunction_instance(speciesID)%twoParticlesMatrix%values &
  !             + WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%externalPotentialMatrix%values))

  !     end if

  !     WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !          WaveFunction_instance( speciesID )%independentSpecieEnergy + &
  !          WaveFunction_instance( SpeciesID )%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy




  !     !! Calcula energia de acoplamiento en caso de mas de una especie presente
  !     if ( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) &
  !          WaveFunction_instance( speciesID )%couplingEnergy = sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !          * WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%couplingMatrix%values )

  !     !print *, "__________________ ENERGY COMPONENTS _______________________"
  !     !print *, "	Specie                       ", MolecularSystem_getNameOfSpecie( speciesID )
  !     !print *, "	Total Energy                =", WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%totalEnergyForSpecie
  !     !print *, "	Indepent Specie Energy      =", WaveFunction_instance( speciesID )%independentSpecieEnergy
  !     !print *, "	Kinetic Energy              =",WaveFunction_instance( speciesID )%kineticEnergy
  !     !print *, "	Puntual Interaction Energy  =",WaveFunction_instance( speciesID )%puntualInteractionEnergy
  !     !print *, "	Independent Particle Energy =",WaveFunction_instance( speciesID )%independentParticleEnergy
  !     !print *, "	N-E Correlation Energy      =",WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy
  !     !print *, "	Repultion Energy            =",WaveFunction_instance( speciesID )%twoParticlesEnergy
  !     !print *, "	Coupling Energy             =", WaveFunction_instance( speciesID )%couplingEnergy
  !     !print *, "____________________________________________________________"

  !   end subroutine WaveFunction_obtainEnergyComponents

  !   !**
  !   ! @brief indica el objeto ha sido instanciado
  !   !
  !   !**
  !   function WaveFunction_isInstanced() result(output)
  !     implicit none
  !     logical :: output

  !     if ( allocated( WaveFunction_instance ) ) then
  !        output = .true.
  !     else
  !        output = .false.
  !     end if

  !   end function WaveFunction_isInstanced


  !   !**
  !   ! @brief Retorna la matrix de Overlap
  !   !
  !   ! @param nameOfSpecie nombre de la especie seleccionada.
  !   !**
  !   function WaveFunction_getOverlapMatrix( nameOfSpecie, sspeciesID ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     integer, optional :: sspeciesID
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     if ( present(sspeciesID) ) then

  !        speciesID = sspeciesID

  !     else

  !        nameOfSpecieSelected = "e-"

  !        if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !        speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     end if

  !     if ( .not.associated(WaveFunction_instance(speciesID)%overlapMatrixValuesPtr) ) then

  !        call WaveFunction_buildOverlapMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%overlapMatrixValuesPtr )

  !   end function WaveFunction_getOverlapMatrix





  !   !**
  !   ! @brief Retorna la matrix de dos particulas para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getTwoParticlesMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( wavefunction_instance(speciesID)%twoParticlesMatrix%values) ) then

  !        call WaveFunction_buildTwoParticlesMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, wavefunction_instance(speciesID)%twoParticlesMatrix )

  !   end function WaveFunction_getTwoParticlesMatrix

  !   !**
  !   ! @brief Retorna la matrix de acoplamiento para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getCouplingMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     if(trim(nameOfSpecieSelected) == "e-ALPHA" .or. trim(nameOfSpecieSelected) == "e-BETA" ) then
  !        call WaveFunction_buildCouplingMatrixElectronFree( nameOfSpecieSelected, output )
  !        return
  !     end if

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%couplingMatrix%values) ) then

  !        call WaveFunction_buildCouplingMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%couplingMatrix )

  !   end function WaveFunction_getCouplingMatrix

  !   !**
  !   ! @brief Retorna la matrix de interaccion con un potencial externo
  !   !
  !   !**
  !   function WaveFunction_getExternalPotentialMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix%values) ) then

  !        call WaveFunction_buildExternalPotentialMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix )

  !   end function WaveFunction_getExternalPotentialMatrix


  !   !**
  !   ! @brief Contruye la matrix de Fock para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getFockMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%fockMatrix%values) ) then

  !        call WaveFunction_buildFockMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%fockMatrix )

  !   end function WaveFunction_getFockMatrix


  !   !**
  !   ! @brief Retorna la matriz  la ultima matriz de densidad calculada para una especie especificada
  !   !
  !   !**
  !   function WaveFunction_getDensityMatrix( nameOfSpecie, sspeciesID ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     integer, optional :: sspeciesID
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     if ( present(sspeciesID) ) then

  !        speciesID = sspeciesID

  !     else

  !        nameOfSpecieSelected = "e-"

  !        if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !        speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     end if

  !     if ( allocated( wavefunction_instance(speciesID)%densityMatrix%values) ) then

  !        call Matrix_copyConstructor( output, wavefunction_instance(speciesID)%densityMatrix )

  !     end if

  !   end function WaveFunction_getDensityMatrix


  !   function WaveFunction_getKineticMatrix( nameOfSpecie ) result(output)
  !     implicit none

  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( allocated( WaveFunction_instance(speciesID)%kineticMatrix%values) ) then

  !        call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%kineticMatrix)

  !     end if

  !   end function WaveFunction_getKineticMatrix


  !   !**
  !   ! @brief Retorna la matrix de acoplamiento para la especie especificada
  !   !        respecto a las particula puntuales
  !   !
  !   !**
  !   function WaveFunction_getPuntualParticleMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     call WaveFunction_buildPuntualParticleMatrix( nameOfSpecieSelected )

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%puntualParticleMatrix )

  !   end function WaveFunction_getPuntualParticleMatrix





  !   function WaveFunction_getValueForOrbitalAt( nameOfSpecie, orbitalNum, coordinate ) result(output)
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     integer :: orbitalNum
  !     real(8) :: coordinate(3)
  !     real(8) :: output

  !     integer :: speciesID
  !     character(30) :: nameOfSpecieSelected
  !     integer :: numberOfContractions
  !     integer :: totalNumberOfContractions
  !     integer :: particleID
  !     integer :: contractionID
  !     integer :: i, j
  !     real(8), allocatable :: auxVal(:)


  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )

  !     output=0.0_8
  !     do i=1,numberOfContractions

  !        particleID = MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%particleID
  !        contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%contractionIDInParticle

  !        totalNumberOfContractions = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital

  !        if( allocated(auxVal)) deallocate(auxVal)
  !        allocate(auxVal(totalNumberOfContractions))

  !        auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( speciesID,  numberOfContraction=i ), coordinate )

  !        do j = 1, totalNumberOfContractions

  !           output = output + auxVal(j) * wavefunction_instance(speciesID)%waveFunctionCoefficients%values(j,orbitalNum)

  !        end do

  !     end do


  !   end function WaveFunction_getValueForOrbitalAt
  !   !
  !   !
  !   subroutine WaveFunction_draw2DOrbital( nameOfSpecie, orbitalNum, flags )
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     integer :: orbitalNum
  !     integer :: flags


  !     character(30) :: nameOfSpecieSelected
  !     character(50) :: fileName
  !     character(50) :: xRange
  !     integer :: speciesID
  !     integer :: numberOfContractions
  !     integer :: j
  !     integer :: i
  !     integer :: numOfGraphs
  !     integer :: auxInitOrbitalNum
  !     integer :: auxLastOrbitalNum

  !     ! 	nameOfSpecieSelected = "e-"
  !     ! 	if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     ! 	speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     ! 	numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
  !     ! 	fileName=trim(CONTROL_instance%INPUT_FILE)//'orbital.'//trim(String_convertIntegerToString(orbitalNum))//'.'//trim(nameOfSpecieSelected)

  !     ! 	select case( flags )

  !     ! 		case(ORBITAL_ALONE)

  !     ! 			open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
  !     ! 			do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
  !     ! 				CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
  !     ! 			write (5,"(F20.10,F20.10)") j*CONTROL_instance%STEP_OF_GRAPHS, &
  !     ! 			WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected,&
  !     ! 			orbitalNum, [0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) !&
  !     ! !			!!! + (WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values(orbitalNum) * CM_NEG1)
  !     ! 			end do
  !     ! 			close(5)
  !     ! 			call InputOutput_make2DGraph(trim(fileName),&
  !     ! 				"Nuclear Wave Function",&
  !     ! 				"r / Bohr",&
  !     ! 				"U / a.u.", &
  !     ! 				y_format="%.2e")

  !     ! 		case(ORBITAL_WITH_POTENTIAL)

  !     ! 			auxInitOrbitalNum=orbitalNum
  !     ! 			auxLastOrbitalNum=orbitalNum
  !     ! 			numOfGraphs=2
  !     ! 			if(orbitalNum==0) then
  !     ! 				auxInitOrbitalNum=1
  !     ! 				auxLastOrbitalNum=numberOfContractions
  !     ! 				numOfGraphs=numberOfContractions+1
  !     ! 			end if
  !     ! 		open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
  !     ! 		do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
  !     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
  !     ! 			write (5,"(2ES20.10$)") &
  !     ! 			j*CONTROL_instance%STEP_OF_GRAPHS, &
  !     ! 			ExternalPotential_getPotential(ExternalPotential_Manager_instance%externalsPots(1),&
  !     ! 			[j*CONTROL_instance%STEP_OF_GRAPHS,0.0_8,0.0_8])*CM_NEG1
  !     ! 			do i=auxInitOrbitalNum,auxLastOrbitalNum
  !     ! 				write (5,"(ES20.10$)") &
  !     ! 				CONTROL_instance%WAVE_FUNCTION_SCALE&
  !     ! 				*WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected, i,&
  !     ! 					[0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) &
  !     ! 				+ (WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values(i) * CM_NEG1)
  !     ! 			end do
  !     ! 			write (5,"(A)") ""
  !     ! 		end do
  !     ! 		close(5)

  !     ! 		xRange=trim(adjustl(String_convertRealToString(real(&
  !     ! 			-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
  !     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))//':'//trim(adjustl(String_convertRealToString(real(&
  !     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
  !     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))

  !     ! 			call InputOutput_make2DGraph(trim(fileName),&
  !     ! 			"Nuclear Wave Function in potential ",&
  !     ! 			"r / Bohr",&
  !     ! 			"U / cm-1",&
  !     ! 			y_format="%.2e",numOfGraphs=numOfGraphs,x_range=trim(xRange))

  !     ! 	end select

  !   end subroutine WaveFunction_draw2DOrbital




  !   !**
  !   ! @brief Resetea los atributos de clase
  !   !**
  !   subroutine WaveFunction_reset()
  !     implicit none


  !     integer :: speciesIterator
  !     integer :: speciesID

  !     do speciesIterator = MolecularSystem_beginSpecie(), MolecularSystem_endSpecie()

  !        speciesID = MolecularSystem_getSpecieID( iteratorOfSpecie=speciesIterator )
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%totalEnergyForSpecie = 0.0_8
  !        WaveFunction_instance( speciesID )%independentSpecieEnergy =0.0_8
  !        WaveFunction_instance( speciesID )%kineticEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%puntualInteractionEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%independentParticleEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%twoParticlesEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%couplingEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%externalPotentialEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%addTwoParticlesMatrix = .false.
  !        WaveFunction_instance( speciesID )%addCouplingMatrix = .false.
  !        WaveFunction_instance( speciesID )%addInterParticleCorrMatrix = .false.
  !        WaveFunction_instance( speciesID )%wasBuiltFockMatrix = .false.

  !        wavefunction_instance(speciesID)%waveFunctionCoefficients%values = 0.0_8
  !        WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values = 0.0_8
  !        WaveFunction_instance( speciesID )%hcoreMatrix%values = 0.0_8
  !        wavefunction_instance(speciesID)%densityMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%transformationMatrix%values = 0.0_8
  !        wavefunction_instance(speciesID)%twoParticlesMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%couplingMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%fockMatrix%values = 0.0_8

  !        if ( associated( WaveFunction_instance(speciesID)%kineticMatrixValuesPtr ) )  &
  !             WaveFunction_instance(speciesID)%kineticMatrixValuesPtr => null()

  !        if ( associated( WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr )) &
  !             WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr => null()

  !        if ( associated( WaveFunction_instance(speciesID)%overlapMatrixValuesPtr )) &
  !             WaveFunction_instance(speciesID)%overlapMatrixValuesPtr => null()


  !     end do

  !   end subroutine WaveFunction_reset


  !   !**
  !   ! @brief Asigna los coeficientes de la función de onda
  !   !**
  !   subroutine WaveFunction_setWaveFunctionCoefficients(wavefunction_instance(speciesID)%waveFunctionCoefficients, nameOfSpecie)
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     character(30) :: nameOfSpecieSelected
  !     type(Matrix), intent(in) :: wavefunction_instance(speciesID)%waveFunctionCoefficients
  !     integer :: speciesID
  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )


  !     call Matrix_copyConstructor( wavefunction_instance(speciesID)%waveFunctionCoefficients, wavefunction_instance(speciesID)%waveFunctionCoefficients)

  !   end subroutine WaveFunction_setWaveFunctionCoefficients


  !   !<
  !   !! Write out matrices in MOLPRO'S format
  !   !>
  !   subroutine WaveFunction_writeMatrices(nameOfSpecie)
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) :: A
  !     type(Matrix) :: P, Pnew
  !     type(Matrix) :: C, Cnew
  !     type(Matrix) :: S, Snew
  !     type(Matrix) :: K, Knew
  !     type(Matrix) :: Pot, PotNew
  !     type(Matrix) :: Jcoup, JcoupNew
  !     type(Matrix) :: Icoup, IcoupNew
  !     type(Exception) :: ex

  !     character(20) :: fileName
  !     character(50) :: nameOfSpecieSelected

  !     integer :: numberOfContractions
  !     integer :: totalNumberOfContractions
  !     integer :: speciesID
  !     integer :: st(1), pt(3), dt(6), ft(10), gt(15)
  !     integer :: angularMoment,numCartesianOrbital
  !     integer, allocatable :: trans(:)
  !     integer :: i, j
  !     integer :: u, v
  !     real(8) :: Epn, Ecn
  !     real(8) :: aVal, bVal

  !     !!Esta rutina no se usa mas... resucitoo....!!!
  !     !return

  !     nameOfSpecieSelected="e-"
  !     if(present(nameOfSpecie)) nameOfSpecieSelected=trim(nameOfSpecie)

  !     fileName = CONTROL_instance%INPUT_FILE

  !     !!Id de la especie seleccionada (por defecto e-)
  !     speciesID = MolecularSystem_getSpecieID(nameOfSpecie = trim(nameOfSpecieSelected))

  !     !!Numero total de contracciones
  !     numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)
  !     totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

  !     !!Tamano del arreglo de nuevas etiquetas
  !     if(allocated(trans)) deallocate(trans)
  !     allocate(trans(totalNumberOfContractions))

  !     !! Reglas de transformacion de indices
  !     st(1)    = 1
  !     pt(1:3)  = [1, 2, 3]
  !     dt(1:6)  = [1, 4, 5, 2, 6, 3]
  !     ft(1:10) = [1, 4, 5, 6, 10, 8, 2, 7, 9, 3]
  !     gt(1:15) = [1, 4, 5, 10, 13, 11, 6, 14, 15, 8, 2, 7, 12, 9, 3]

  !     write(*,*) "RULES FOR INDEX TRANSFORMATION"
  !     write(*,*) "================================"
  !     write(*,*) ""
  !     write(*,"(' s =',I3)") st(1)
  !     write(*,"(' p =',3I3)") ( pt(i), i=1,3 )
  !     write(*,"(' d =',6I3)") ( dt(i), i=1,6 )
  !     write(*,"(' f =',10I3)") ( ft(i), i=1,10 )
  !     write(*,"(' g =',15I3)") ( gt(i), i=1,15 )
  !     write(*,*) ""

  !     trans = 0
  !     u = 1
  !     v = 0

  !     do i = 1, numberOfContractions
  !        angularMoment = MolecularSystem_instance%particlesPtr(MolecularSystem_instance%idsOfContractionsForSpecie(&
  !             speciesID)%contractionID(i)%particleID)%basis%contractions( &
  !             MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(&
  !             i)%contractionIDInParticle)%angularMoment

  !        select case(CONTROL_instance%DIMENSIONALITY)
  !        case(3)
  !           numCartesianOrbital = ((angularMoment + 1_8)*(angularMoment + 2_8))/2_8
  !        case(2)
  !           numCartesianOrbital = ((angularMoment + 1_8))
  !        case(1)
  !           numCartesianOrbital = 1 
  !        case default
  !           call Exception_constructor( ex , WARNING )
  !           call Exception_setDebugDescription( ex, "Class object WaveFunction in the writeMatrices function" )
  !           call Exception_setDescription( ex, "this Dimensionality is not implemented (D>3)" )
  !           call Exception_show( ex )

  !           numCartesianOrbital = 0
  !        end select

  !        select case(angularMoment)
  !        case(0)
  !           trans(u) = st(1) + v
  !           u = u + 1

  !        case(1)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = pt(j) + v
  !              u = u + 1
  !           end do

  !        case(2)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = dt(j) + v
  !              u = u + 1
  !           end do
  !        case(3)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = ft(j) + v
  !              u = u + 1
  !           end do
  !        case(4)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = gt(j) + v
  !              u = u + 1
  !           end do

  !        case default
  !           call Exception_constructor( ex , WARNING )
  !           call Exception_setDebugDescription( ex, "Class object WaveFunction in the writeMatrices function" )
  !           call Exception_setDescription( ex, "this angular moment is not implemented (l>4)" )
  !           call Exception_show( ex )

  !           return

  !        end select
  !        v = u - 1
  !     end do
  !     i = 1
  !     !write(*,"(' trans =',<totalNumberOfContractions>I3)") ( trans(i), i=1,totalNumberOfContractions )

  !     !Epn = MolecularSystem_instance%twoParticlesEnergy(2)
  !     !Ecn = MolecularSystem_instance%kineticEnergy(2)

  !     !write(*,*) "Epn = ", Epn
  !     !write(*,*) "Ecn = ", Ecn

  !     P = WaveFunction_getDensityMatrix( trim(nameOfSpecieSelected) )
  !     C = WaveFunction_getWaveFunctionCoefficients( trim(nameOfSpecieSelected) )
  !     S = WaveFunction_getOverlapMatrix( trim(nameOfSpecieSelected) )
  !     K = WaveFunction_getKineticMatrix( trim(nameOfSpecieSelected) )
  !     Pot= WaveFunction_getExternalPotentialMatrix(trim(nameOfSpecieSelected)) 
  !     Jcoup = WaveFunction_getCouplingMatrix( trim(nameOfSpecieSelected) )
  !     Icoup = WaveFunction_getPuntualParticleMatrix( trim(nameOfSpecieSelected) )

  !     call Matrix_constructor( Pnew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Cnew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Snew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Knew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( JcoupNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( IcoupNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )    
  !     call Matrix_constructor( PotNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )

  !     do i=1,totalNumberOfContractions
  !        do j=1,totalNumberOfContractions
  !           Pnew%values( trans(i), trans(j) ) = P%values( i, j )
  !           Cnew%values( trans(i), trans(j) ) = C%values( i, j )
  !           Snew%values( trans(i), trans(j) ) = S%values( i, j )
  !           Knew%values( trans(i), trans(j) ) = K%values( i, j )
  !           JcoupNew%values( trans(i), trans(j) ) = Jcoup%values( i, j )
  !           IcoupNew%values( trans(i), trans(j) ) = Icoup%values( i, j )
  !           PotNew%values( trans(i), trans(j) ) = Pot%values( i, j )
  !        end do
  !     end do

  !     ! 		do i=1,totalNumberOfContractions
  !     ! 			aVal = Pnew.values( 3, i )
  !     ! 			bVal = Pnew.values( 6, i )
  !     ! 			Pnew.values( 6, i ) = aVal
  !     ! 			Pnew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = Cnew.values( 3, i )
  !     ! 			bVal = Cnew.values( 6, i )
  !     ! 			Cnew.values( 6, i ) = aVal
  !     ! 			Cnew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = Snew.values( 3, i )
  !     ! 			bVal = Snew.values( 6, i )
  !     ! 			Snew.values( 6, i ) = aVal
  !     ! 			Snew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = JcoupNew.values( 3, i )
  !     ! 			bVal = JcoupNew.values( 6, i )
  !     ! 			JcoupNew.values( 6, i ) = aVal
  !     ! 			JcoupNew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = IcoupNew.values( 3, i )
  !     ! 			bVal = IcoupNew.values( 6, i )
  !     ! 			IcoupNew.values( 6, i ) = aVal
  !     ! 			IcoupNew.values( 3, i ) = bVal
  !     ! 		end do


  !     write(*,"(3A)", advance="no") " Saving coefficients matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"coeff"), " ) ... "
  !     open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"coeff")), action='write', form='unformatted' )
  !     write(20) int(size(Cnew%values), 8)
  !     write(20) Cnew%values
  !     close(20)
  !     write(*,*) "OK"

  !     if(trim(nameOfSpecieSelected) /= "e-BETA") then 

  !        write(*,"(3A)", advance="no") " Saving potential matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"pot"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"pot")), action='write', form='unformatted')
  !        write(20) int(size(PotNew%values), 8)
  !        write(20) PotNew%values
  !        close(20)
  !        write(*,*) "OK"


  !        write(*,"(3A)", advance="no") " Saving density matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"dens"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"dens")), action='write', form='unformatted' )
  !        write(20) int(size(Pnew%values), 8)
  !        write(20) Pnew%values
  !        close(20)
  !        write(*,*) "OK"

  !        !		write(*,*) "Jcoup ="
  !        !		call Matrix_show( JcoupNew )

  !        !		write(*,*) "Icoup ="
  !        !		call Matrix_show( IcoupNew )

  !        write(*,"(3A)", advance="no") " Saving overlap matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"over"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"over")), action='write', form='unformatted' )
  !        write(20) int(size(Snew%values), 8)
  !        write(20) Snew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving kinetic matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"kin"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"kin")), action='write', form='unformatted' )
  !        write(20) int(size(Snew%values), 8)
  !        write(20) Knew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving coupling matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"jcoup"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"jcoup")), action='write', form='unformatted' )
  !        write(20) int(size(JcoupNew%values), 8)
  !        write(20) JcoupNew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving fixed interaction matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"icoup"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"icoup")), action='write', form='unformatted' )
  !        write(20) int(size(IcoupNew%values), 8)
  !        write(20) IcoupNew%values
  !        close(20)
  !        write(*,*) "OK"

  !     end if

  !   end subroutine WaveFunction_writeMatrices

  !>
  !! @brief calcula la energia total de acoplamiento para una especie especificada
  ! function WaveFunction_getTotalCouplingEnergy() result( output )
  !   implicit none
  !   real(8) :: output

  !   character(30) :: nameOfSpecie
  !   character(30) :: nameOfOtherSpecie
  !   integer :: speciesID
  !   integer :: otherSpeciesID
  !   real(8) :: factor

  !   output = 0.0_8
  !   factor = 0.5

  !   do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

  !     nameOfSpecie = MolecularSystem_getNameOfSpecie( speciesID ) 
  !     call WaveFunction_buildCouplingMatrix(nameOfSpecie)

  !     if ( nameOfSpecie == CONTROL_instance%SCF_GHOST_SPECIES ) factor = 1.0

  !     do otherSpeciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

  !       if ( otherSpeciesID /= speciesID ) then 

  !         nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID ) 
  !         !if ( nameOfOtherSpecie == CONTROL_instance%SCF_GHOST_SPECIES  ) factor = 0

  !         output = output + factor*(sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !             * (wavefunction_instance(speciesID)%hartreeMatrix(otherSpeciesID)%values))) 

  !       end if
  !     end do

  !   end do

  !   ! real(8) :: auxValue
  !   ! real(8) :: auxRepulsion
  !   ! real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !   ! integer :: numberOfContractions
  !   ! integer :: numberOfContractionsOfOtherSpecie
  !   ! integer :: numberOfTotalContractions
  !   ! integer :: numberOfTotalContractionsOfOtherSpecie
  !   ! integer :: outFile
  !   ! integer :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !   ! integer :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !   ! integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !   ! integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !   ! integer :: u, v
  !   ! integer :: k, l, m
  !   ! integer :: arrayNumber
  !      !!output=output+ 0.5*(sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !      !!     * (wavefunction_instance(speciesID)%couplingMatrix%values))) 


  !      !       do otherSpeciesID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies()


  !      ! !! Restringe suma de terminos repulsivos de la misma especie.
  !      ! if ( otherSpeciesID /= speciesID ) then

  !      !    nameOfSpecie = MolecularSystem_getNameOfSpecie( speciesID )
  !      !    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
  !      !    numberOfTotalContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

  !      !    nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpeciesID )
  !      !    numberOfContractionsOfOtherSpecie = MolecularSystem_getNumberOfContractions( otherSpeciesID )
  !      !    numberOfTotalContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( otherSpeciesID )

  !      !    !Restringe la suma de terminos repulsivos electronicos
  !      !    if(trim(nameOfSpecie)=="E-ALPHA" .and. trim(nameOfOtherSpecie)=="E-BETA") cycle

  !      !    auxValue = 0.0_8
  !      !    m = 0

  !      !    !! open file for integrals
  !      !    open(UNIT=34,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
  !      !         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

  !      !    auxValue=0.0_8

  !      !    readIntegrals : do

  !      !       read(34) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !      !            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !      !            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !      !       do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

  !      !          if (a(u) == -1) exit readIntegrals

  !      !          m = m + 1

  !      !          auxValue = auxValue +&
  !      !               (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
  !      !               * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
  !      !               *  integral(u))

  !      !          if(b(u) /= a(u)) then

  !      !             m = m + 1

  !      !             auxValue = auxValue +&
  !      !                  (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
  !      !                  * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
  !      !                  *  integral(u))
  !      !          end if

  !      !          if(s(u) /= r(u)) then

  !      !             m = m + 1

  !      !             auxValue = auxValue +&
  !      !                  (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
  !      !                  * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
  !      !                  *  integral(u))
  !      !          end if

  !      !          if(b(u) /= a(u) .and. s(u) /= r(u)) then

  !      !             m = m + 1

  !      !             auxValue = auxValue +&
  !      !                  (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
  !      !                  * WaveFunction_instance( otherSpeciesID)%densityMatrix%values(r(u),s(u)) &
  !      !                  *  integral(u))
  !      !          end if


  !      !       end do

  !      !    end do readIntegrals

  !      !    auxValue = auxValue *  MolecularSystem_getCharge( speciesID=speciesID ) &
  !      !         * MolecularSystem_getCharge( speciesID=otherSpeciesID )

  !      !    output = output + auxValue

  !      !    close(34)

  !      ! end if
  ! end function WaveFunction_getTotalCouplingEnergy

  !>
  !! @brief retorna la matrix de particula independiente.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  ! function WaveFunction_getHcoreMatrix( nameOfSpecie ) result(output)
  !   implicit none
  !   character(*), optional :: nameOfSpecie
  !   type(Matrix) :: output

  !   integer :: speciesID
  !   character(30) :: nameOfSpecieSelected

  !   nameOfSpecieSelected = "E-"
  !   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !   speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !   if ( .not. allocated(WaveFunction_instance( speciesID )%HcoreMatrix%values ) ) then

  !      call WaveFunction_exception(ERROR, "You need to build the Hcore matrix before to use it.", "At HF program, at WaveFunction_getHcoreMatrix function.")

  !   end if

  !   call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%HcoreMatrix )

  ! end function WaveFunction_getHcoreMatrix

  !>
  !! @brief Retorna la matrix de transformationMatrix
  ! function WaveFunction_getTransformationMatrix( nameOfSpecie ) result( output )
  !   implicit none
  !   character(*), optional :: nameOfSpecie
  !   type(Matrix) ::  output

  !   character(30) :: nameOfSpecieSelected
  !   integer :: speciesID

  !   nameOfSpecieSelected = "e-"
  !   if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !   speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !   if ( .not.allocated(WaveFunction_instance(speciesID)%transformationMatrix%values) ) then

  !      call WaveFunction_exception(ERROR, "You need build the transformation Matrix before to use it.", "At HF program, at WaveFunction_getTransformationMatrix function.")

  !   end if

  !   call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%transformationMatrix )

  ! end function WaveFunction_getTransformationMatrix

  !>
  !! @brief Ajusta la matriz de densidad para una especie espcificada
  ! subroutine WaveFunction_setDensityMatrix( densityMatrix, speciesID )
  !   implicit none

  !   type(Matrix), intent(in) :: densityMatrix
  !   integer :: speciesID

  !   integer :: totalNumberOfContractions
  !   character(50) :: densFile, labels(2)
  !   integer :: densUnit

  !   totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

  !   if( .not. allocated(WaveFunction_instance( speciesID )%densityMatrix%values )) then
  !      call Matrix_constructor( WaveFunction_instance( speciesID )%densityMatrix, &
  !           int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
  !   end if

  !   call Matrix_copyConstructor( WaveFunction_instance( speciesID )%densityMatrix, densityMatrix )

  !   !! Debug
  !   ! print*, "Matriz de densidad inicial ", MolecularSystem_getNameOfSpecie(speciesID)
  !   ! call Matrix_show(WaveFunction_instance( speciesID )%densityMatrix)

  !   !!Save this matrix for DFT calculations, because reasons
  !   if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       
  !      densUnit = 78
  !      densFile = trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//".densmatrix"
  !      open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")

  !      labels(1) = "DENSITY-MATRIX"
  !      labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
     
  !      call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )

  !      close (78)
  !   end if

  ! end subroutine WaveFunction_setDensityMatrix

