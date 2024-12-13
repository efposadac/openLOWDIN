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
!! @brief Modulo para realizacion de iteraciones SCF sobre diferentes especies cuanticas
!!
!!  Este modulo define una seudoclase para realizacion de ciclos SCF sobre todas las
!! especies cuanticas.
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-30
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y metodos
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe el modulo para su inclusion en Lowdin
!!   - <tt> 2023-03-10 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
!!        -# Rereescribe el modulo para simplificarlo
module MultiSCF_
  use Exception_
  use Stopwatch_
  use CONTROL_
  use Matrix_
  use List_  
  use MolecularSystem_
  use WaveFunction_
  use DensityFunctionalTheory_
  use SingleSCF_
  use omp_lib
  use DensityMatrixSCFGuess_
  use OrbitalLocalizer_
  use Convergence_
  use Libint2Interface_
  
  implicit none    

  !< For effective Fock matrix for restricted open-shell SCF
  ! logical :: isROHF
  !>

  type, public :: MultiSCF

     type(MolecularSystem), pointer :: molSys
     type(List) :: energyOMNE
     character(100) :: name
     integer :: numberOfIterations
     integer :: status
     integer, allocatable :: singleMaxIterations(:)
     real(8), allocatable :: singleEnergyTolerance(:)
     real(8), allocatable :: singleDensityTolerance(:)

     !! Global energies
     real(8) :: totalEnergy
     real(8) :: totalCouplingEnergy
     real(8) :: totalPotentialEnergy
     real(8) :: totalKineticEnergy

     !! Global Density Standard Deviation
     real(8) :: totalDensityMatrixStandardDeviation

     !! Cosmo
     real(8) :: cosmo3Energy

     !!
     logical :: printSCFiterations

     !!
     type(Grid), allocatable :: DFTGrids(:), DFTGridsCommonPoints(:,:)

  end type MultiSCF

  type(MultiSCF), public, target :: MultiSCF_instance

  public :: &
       MultiSCF_constructor, &
       MultiSCF_destructor, &
       MultiSCF_getNumberOfIterations, &
       MultiSCF_getLastEnergy, &
       MultiSCF_iterate, &
       MultiSCF_restart, &
       MultiSCF_solveHartreeFockRoothan, &
       MultiSCF_obtainFinalEnergy, &
       MultiSCF_reset

contains


  !>
  !! @brief Define el constructor para la clase
  subroutine MultiSCF_constructor(this,wfObjects,iterationScheme,molsystem)
    implicit none
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)
    integer :: iterationScheme
    type(MolecularSystem), target :: molsystem
    
    integer :: i, nspecies, speciesID
    integer :: dftUnit
    character(50) :: labels(2)
    character(50) :: dftFile
    
    this%molSys=>molsystem

    nspecies=MolecularSystem_getNumberOfQuantumSpecies(this%molSys)
    
    ! isROHF = .false.
    select case(iterationScheme)
    case(0)
       this%name="NONELECTRONIC_FULLY_CONVERGED_BY_ELECTRONIC_ITERATION"
    case(1)
       this%name="ELECTRONIC_FULLY_CONVERGED_BY_NONELECTRONIC_ITERATION"
    case(2)
       this%name="SPECIES_FULLY_CONVERGED_INDIVIDUALLY"
    case(3)
       this%name="SIMULTANEOUS"
    case default
       call MultiSCF_exception( ERROR, "selected an iteration scheme not implemented (>3)", "at SCF program, this_constructor")
    end select

    call List_constructor( this%energyOMNE,"ENERGY", CONTROL_instance%LISTS_SIZE)
    this%numberOfIterations = 0
    this%status = 0

    allocate(this%singleEnergyTolerance(nspecies),&
         this%singleDensityTolerance(nspecies),&
         this%singleMaxIterations(nspecies))

    do i = 1, nspecies
       if(this%molSys%species(i)%isElectron ) then
          this%singleEnergyTolerance(i)=CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE
          this%singleDensityTolerance(i)=CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
       else
          this%singleEnergyTolerance(i)=CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE
          this%singleDensityTolerance(i)=CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
       end if

       !The maximum number of iterations for each species in selected according to the SCF scheme chosen
       select case(iterationScheme)

       case( 0 )          
          ! we perform single species iterations for nonelectrons
          if(this%molSys%species(i)%isElectron ) then
             this%singleMaxIterations(i)=1
          else
             this%singleMaxIterations(i)=CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
          end if

       case( 1 )
          ! we perform single species iterations for nelectrons
          if(this%molSys%species(i)%isElectron ) then
             this%singleMaxIterations(i)=CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
          else
             this%singleMaxIterations(i)=1
          end if


       case( 2 )
          ! we perform single species  for all species
          if(this%molSys%species(i)%isElectron ) then
             this%singleMaxIterations(i)=CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
          else
             this%singleMaxIterations(i)=CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
          end if

       case ( 3 )
          ! we do not perform single species SCF
          if(this%molSys%species(i)%isElectron ) then
             this%singleMaxIterations(i)=1
          else
             this%singleMaxIterations(i)=1
          end if
       end select

    end do

    this%totalEnergy = 0.0_8
    this%cosmo3Energy = 0.0_8
    this%totalCouplingEnergy = 0.0_8

    !!Print information - only if the method requires few (one) SCF
    this%printSCFiterations=.true.
    if(CONTROL_instance%OPTIMIZE) this%printSCFiterations=.false.
    if(CONTROL_instance%PRINT_LEVEL .eq. 0) this%printSCFiterations=.false. 
    if(CONTROL_instance%DEBUG_SCFS) this%printSCFiterations=.true.

    !! Start the wavefunction object   
    call WaveFunction_constructor(wfObjects,nspecies,this%molSys)

    !!Initialize DFT: Calculate Grids and build functionals
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          call system("lowdin-DFT.x BUILD_SCF_GRID")
          do speciesID = 1, nspecies
             dftUnit = 77
             dftFile = "lowdin."//trim(MolecularSystem_getNameOfSpecies(speciesID,this%molSys))//".grid"
             open(unit = dftUnit, file=trim(dftFile), status="old", form="unformatted")

             labels(2) = MolecularSystem_getNameOfSpecies(speciesID,this%molSys)
             labels(1) = "EXACT-EXCHANGE-FRACTION"

             call Vector_getFromFile(unit=dftUnit, binary=.true., value=wfObjects(speciesID)%exactExchangeFraction, arguments=labels)
             close(unit=dftUnit)
             ! print *, "el tormento tuyo", speciesID, these(speciesID)%exactExchangeFraction
          end do
       else     !! Allocate DFT grids memory.
          if(allocated(this%DFTGrids)) deallocate(this%DFTGrids)
          allocate(this%DFTGrids(nspecies))

          if (allocated(this%DFTGridsCommonPoints)) deallocate(this%DFTGridsCommonPoints)
          allocate(this%DFTGridsCommonPoints(nspecies,nspecies))
          
          call DensityFunctionalTheory_buildSCFGrid(this%DFTGrids,this%DFTGridsCommonPoints,wfObjects(1:nspecies)%exactExchangeFraction,this%molSys)
       end if
    end if

    
    !! Start the orbital localizer object
    if (CONTROL_instance%LOCALIZE_ORBITALS) call OrbitalLocalizer_constructor( )

  end subroutine MultiSCF_constructor

  !>
  !! @brief Define el destructor para clase
  subroutine MultiSCF_destructor(this)
    implicit none
    type(MultiSCF) :: this

    deallocate(this%singleEnergyTolerance,this%singleDensityTolerance,this%singleMaxIterations)
    call List_destructor( this%energyOMNE )

  end subroutine MultiSCF_destructor

  !>
  !! @brief Retorna la energia obtenida en la ultima iteracion iter-especies
  function MultiSCF_getNumberOfIterations(this) result( output )
    implicit none

    type(MultiSCF) :: this
    integer :: output

    output= this%numberOfIterations

  end function MultiSCF_getNumberOfIterations

  !>
  !! @brief Retorna la energia obtenida en la ultima iteracion iter-especies
  function MultiSCF_getLastEnergy(this) result( output )
    implicit none

    type(MultiSCF) :: this
    real(8) :: output

    call List_end( this%energyOMNE )
    output= List_current( this%energyOMNE )

  end function MultiSCF_getLastEnergy

  !>
  !! @brief Realiza esquema de iteracion SCF para todas las especies cuanticas presentes
  !! @param
  !!
  subroutine MultiSCF_iterate(this, wfObjects, libint2Objects, iterationScheme)
    implicit none
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)
    type(Libint2Interface) :: libint2Objects(:)    
    integer, intent(in) :: iterationScheme

    integer :: i,j
    integer :: numberOfSpecies
    character(50) :: nameOfSpecies
    character(50) :: densFile
    integer :: singleIterator
    real(8) :: oldEnergy
    real(8) :: deltaEnergy

!!!We start with an update of the global energy and matrices

    this%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies(this%molSys)

    densFile="lowdin.densmatrix"

    !!DFT calculations are performed only in global iterations
    !!Save density matrices to file for DFT calculations
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          call WaveFunction_writeDensityMatricesToFile(wfObjects, densFile)
          call system("lowdin-DFT.x SCF_DFT "//trim(densFile))
       else
          call WaveFunction_getDFTContributions(wfObjects,this%DFTGrids,this%DFTGridsCommonPoints,"SCF")
       end if
    end if

    !Build Fock matrices and calculate energy with the initial densities
    !Coupling Matrix is only updated in global SCF cycles
    do i = 1, numberOfSpecies

       call WaveFunction_buildTwoParticlesMatrix(wfObjects(i),libint2Objects=libint2Objects(:))

       call WaveFunction_buildCouplingMatrix(wfObjects,i,libint2Objects=libint2Objects(:))

       if ( (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") .and. CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          call WaveFunction_readExchangeCorrelationMatrix(wfObjects(i))
       end if

       if (CONTROL_instance%COSMO) then
          call WaveFunction_buildCosmo2Matrix(wfObjects(i))
          call WaveFunction_buildCosmoCoupling(wfObjects(i))
       end if

       call WaveFunction_buildFockMatrix(wfObjects(i))

       call WaveFunction_obtainTotalEnergyForSpecies(wfObjects(i))

       !Save initial matrices for convergence comparisons
       if(this%numberOfIterations .eq. 0 ) then
          call Convergence_setInitialDensityMatrix(wfObjects(i)%convergenceMethod, &
               wfObjects(i)%densityMatrix )
          call Convergence_setInitialFockMatrix(wfObjects(i)%convergenceMethod, &
               wfObjects(i)%fockMatrix )
       end if
       
    end do

    !Get total energy
    call WaveFunction_obtainTotalEnergy(wfObjects,&
         this%totalEnergy, &
         this%totalCouplingEnergy, &
         this%cosmo3Energy)

    call List_push_back( this%energyOMNE, this%totalEnergy)             
    this%numberOfIterations = this%numberOfIterations + 1

!!!Now we procede to update each species density matrices according to the iteration scheme selected
    this%totalDensityMatrixStandardDeviation=0.0
    do i = 1, numberOfSpecies
       nameOfSpecies = MolecularSystem_getNameOfSpecies(i,this%molSys)
       oldEnergy=wfObjects(i)%totalEnergyForSpecies
       deltaEnergy=1.0E16_8
       singleIterator=0

       !The maximum number of iterations for each species in selected according to the SCF scheme chosen

       do while( (abs(deltaEnergy) .gt. this%singleEnergyTolerance(i) .or. &
            Matrix_standardDeviation( wfObjects(i)%beforeDensityMatrix, wfObjects(i)%densityMatrix ) .gt. &
            this%singleDensityTolerance(i)) .and. &
            singleIterator .lt. this%singleMaxIterations(i) )

          singleIterator=singleIterator+1
          !! Obtains the new eigenvectors
          call SingleSCF_iterate(wfObjects(i))

          !! Determina la desviacion estandar de los elementos de la matriz de densidad
          call Matrix_copyConstructor( wfObjects(i)%beforeDensityMatrix, wfObjects(i)%densityMatrix )
          call WaveFunction_buildDensityMatrix(wfObjects(i))
          call List_push_back( wfObjects(i)%standardDesviationOfDensityMatrixElements, &
               Matrix_standardDeviation( wfObjects(i)%beforeDensityMatrix, wfObjects(i)%densityMatrix ) )

          !! Calcula la energy de la especie con la nueva densidad
          call WaveFunction_obtainTotalEnergyForSpecies(wfObjects(i))
          call List_push_back( wfObjects(i)%energySCF, wfObjects(i)%totalEnergyForSpecies )
          call List_push_back( wfObjects(i)%diisError, Convergence_getDiisError( wfObjects(i)%convergenceMethod) )

          if(this%singleMaxIterations(i).gt.1) then
             !! Updates two particle matrix
             call WaveFunction_buildTwoParticlesMatrix(wfObjects(i),libint2Objects=libint2Objects(:))

             if (CONTROL_instance%COSMO) then
                call WaveFunction_buildCosmo2Matrix(wfObjects(i))
             end if
             !!Builds new fock Matrix
             call WaveFunction_buildFockMatrix(wfObjects(i))

             !!Checks energy convergence
             deltaEnergy = oldEnergy -wfObjects(i)%totalEnergyForSpecies
             oldEnergy = wfObjects(i)%totalEnergyForSpecies

             !!Prints iteration results
             if (  CONTROL_instance%DEBUG_SCFS) then
                write(*,"(A10,I5,F20.12,F20.12,F20.12)") trim(nameOfSpecies), singleIterator , &
                     wfObjects(i)%totalEnergyForSpecies, deltaEnergy, &
                     List_current(wfObjects(i)%standardDesviationOfDensityMatrixElements)
             end if

             if ( this%printSCFiterations ) then
                !!Prints convergence messages
                if(singleIterator .ge. this%singleMaxIterations(i) ) &
                     write(*,"(T35,A10,A30,I4,A)") trim(nameOfSpecies), " Max. subcycles reached(", this%singleMaxIterations(i),")"

                if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "DENSITY" .and. &
                     Matrix_standardDeviation( wfObjects(i)%beforeDensityMatrix, wfObjects(i)%densityMatrix ) .lt. &
                     this%singleDensityTolerance(i)) & 
                     write(*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Density converged in", singleIterator ," subcycles"

                if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "ENERGY" .and. &
                     abs(deltaEnergy) .lt. this%singleEnergyTolerance(i) )&
                     write(*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Energy converged in", singleIterator ," subcycles"

                if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "BOTH" .and. & 
                     Matrix_standardDeviation( wfObjects(i)%beforeDensityMatrix, wfObjects(i)%densityMatrix ) .lt. &
                     this%singleDensityTolerance(i) .and. & 
                     abs(deltaEnergy) .lt. this%singleEnergyTolerance(i) )&
                     write(*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Energy-density converged in", singleIterator ," subcycles"
             end if
          end if

       end do
       this%totalDensityMatrixStandardDeviation=this%totalDensityMatrixStandardDeviation+&
            sqrt(List_current(wfObjects(i)%standardDesviationOfDensityMatrixElements)**2)
    end do

    if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
         (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
       i=MolecularSystem_getSpecieIDFromSymbol(trim("E-ALPHA"),this%molSys)
       j=MolecularSystem_getSpecieIDFromSymbol(trim("E-BETA"),this%molSys)

       if(MolecularSystem_getNumberOfParticles(i,this%molSys) .eq. MolecularSystem_getNumberOfParticles(j,this%molSys) ) then
          wfObjects(j)%waveFunctionCoefficients%values= wfObjects(i)%waveFunctionCoefficients%values
          wfObjects(j)%densityMatrix%values= wfObjects(i)%densityMatrix%values
       end if

    end if


  end subroutine MultiSCF_iterate

  ! >
  ! @brief Prueba si la energia del la ultima iteracion a sufrido un cambio por debajo de cierta tolerancia
  ! function MultiSCF_checkConvergence(  tolerace ) result( output )
  !   implicit none
  !   real(8), intent(in) :: tolerace
  !   integer :: output

  ! character(30) :: nameOfSpecie
  ! real(8) :: deltaEnergy
  ! real(8) :: finalEnergy
  ! real(8) :: toleraceOfSpecie
  ! type(Exception) :: ex
  ! integer :: speciesID

  ! if ( this%numberOfIterations > CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS ) then

  !    output = SCF_GLOBAL_CONVERGENCE_SUCCESS
  !    call List_end( this%energyOMNE )
  !    finalEnergy= List_current( this%energyOMNE )

  !    !! Obtiene el valor  de energia anterior a la ultima iteracion
  !    call List_iterate( this%energyOMNE, -1 )

  !    !! Obtiene el cambio de energia en las ultimas dos iteraciones
  !    deltaEnergy = finalEnergy - List_current( this%energyOMNE )
  !    print *,""
  !    write(*,"(A20,A20)") "   Current Energy   ", "   Current Change   "
  !    print *,"-----------------------------------------------"
  !    write(*,"(F20.10,F20.10)") finalEnergy, deltaEnergy

  !    print *,  "The number of Iterations was exceded, the convergence had failed"

  ! else

  !    if ( this%numberOfIterations > 1 ) then

  !       auxVar=.true.
  !       do speciesID = 1, nspecies

  !          nameOfSpecie = MolecularSystem_getNameOfSpecies(speciesID)

  !          toleraceOfSpecie = this%electronicTolerance

  !          if (.not. this%molSys%species(speciesID)%isElectron ) then
  !             toleraceOfSpecie = this%nonelectronicTolerance
  !          end if

  !          if ( SingleSCF_testDensityMatrixChange( nameOfSpecie, &
  !               toleraceOfSpecie ) == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then
  !             auxVar = auxVar .and. .true.
  !          else
  !             auxVar = auxVar .and. .false.
  !          end if

  !       end do

  !       call List_end( this%energyOMNE )
  !       deltaEnergy= List_current( this%energyOMNE )

  !       !! Obtiene el valor  de energia anterior a la ultima iteracion
  !       call List_iterate( this%energyOMNE, -1 )

  !       !! Obtiene el cambio de energia en las ultimas dos iteraciones
  !       deltaEnergy = deltaEnergy - List_current( this%energyOMNE )

  !       if( ( ( abs( deltaEnergy ) < tolerace ) .and. auxVar ) .or. abs(deltaEnergy) < CONTROL_instance%TOTAL_ENERGY_TOLERANCE ) then
  !          output = SCF_GLOBAL_CONVERGENCE_SUCCESS
  !       else
  !          output = SCF_GLOBAL_CONVERGENCE_CONTINUE
  !       end if

  !    else

  !       output = SCF_GLOBAL_CONVERGENCE_CONTINUE

  !    end if

  ! end if

  ! end function MultiSCF_checkConvergence



  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine MultiSCF_restart(this)
    implicit none
    type(MultiSCF) :: this  

    this%numberOfIterations = 0

    call List_clear( this%energyOMNE )

  end subroutine MultiSCF_restart

  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine MultiSCF_reset(this,wfObjects)
    implicit none
    type(MultiSCF) :: this  
    type(WaveFunction) :: wfObjects(*)

    integer :: speciesIterator

    this%numberOfIterations = 0

    call List_clear( this%energyOMNE )
    this%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

    do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies(this%molSys)
       call SingleSCF_reset(wfObjects(speciesIterator))
    end do

  end subroutine MultiSCF_reset

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine MultiSCF_exception( typeMessage, description, debugDescription)
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

  end subroutine MultiSCF_exception

  !>
  !! @brief create hcore operators and use them to obtain a guess HC=eSC for the SCF equations, or read previous coefficients
  subroutine MultiSCF_buildHcore(this,wfObjects)
    type(MultiSCF) :: this  
    type(WaveFunction) :: wfObjects(*)

    integer :: numberOfSpecies
    integer :: wfnUnit, densUnit
    integer :: speciesID
    integer :: integralsUnit
    character(50) :: wfnFile, densFile
    character(50) :: integralsFile
    logical :: existFile

    !!cosmo things
    character(50) :: cosmoIntegralsFile

    if ( trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DISK" ) then
       !! Open file for wfn
       wfnUnit = 300
       wfnFile = "lowdin.wfn"

       integralsUnit = 30
       integralsFile = "lowdin.opints"

       densUnit=78
       densFile="lowdin.densmatrix"

       !**************************************************************************************************
       !! Builds the fock operator
       !!
       !! Check the one-particle integrals file  
       existFile = .false.     
       inquire(file=trim(integralsFile), exist=existFile)
       if( .not. existFile ) &
            call MolecularSystem_exception(ERROR,"lowdin.opints file not found!", "In SCF.f90 at main program")
       open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")
       read(integralsUnit) numberOfSpecies
       if(this%molSys%numberOfQuantumSpecies /= numberOfSpecies ) &
            call MolecularSystem_exception( ERROR, "Bad "//trim(integralsFile)//" file!", "In SCF.f90 at main program")
       close(integralsUnit)
       do speciesID = 1, this%molSys%numberOfQuantumSpecies
          call WaveFunction_readOverlapMatrix(wfObjects(speciesID), trim(integralsFile))
          call WaveFunction_readKineticMatrix(wfObjects(speciesID), trim(integralsFile))
          call WaveFunction_readPuntualInteractionMatrix(wfObjects(speciesID), trim(integralsFile))
          if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
             call WaveFunction_readExternalPotentialMatrix(wfObjects(speciesID), trim(integralsFile))
          end if
          if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
             call WaveFunction_readElectricFieldMatrices(wfObjects(speciesID), trim(integralsFile))
          end if

          if ( MolecularSystem_getOmega(speciesID,this%molSys) .ne. 0.0_8) then
            call WaveFunction_readHarmonicOscillatorMatrix(wfObjects(speciesID), trim(integralsFile))
          end if
          !! Builds Cosmo hcore integrals
          if(CONTROL_instance%COSMO)then
             cosmoIntegralsFile="cosmo.opints"
             call WaveFunction_cosmoHcoreMatrix(wfObjects(speciesID), trim(cosmoIntegralsFile))
          end if
       end do
    else !!DIRECT or MEMORY
       numberOfSpecies = this%molSys%numberOfQuantumSpecies

       do speciesID = 1, numberOfSpecies
          call DirectIntegralManager_getOverlapIntegrals(this%molSys,speciesID,&
               wfObjects(speciesID)%overlapMatrix)
          call DirectIntegralManager_getKineticIntegrals(this%molSys,speciesID,&
               wfObjects(speciesID)%kineticMatrix)
          call DirectIntegralManager_getAttractionIntegrals(this%molSys,speciesID,&
               wfObjects(speciesID)%puntualInteractionMatrix)
          if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
             call DirectIntegralManager_getExternalPotentialIntegrals(this%molSys,speciesID,&
                  wfObjects(speciesID)%externalPotentialMatrix)
          end if
       end do
    end if

    !!**********************************************************
    do speciesID = 1, this%molSys%numberOfQuantumSpecies
       !! Transformation Matrix
       call WaveFunction_buildTransformationMatrix(wfObjects(speciesID), 2)
       !! Hcore Matrix
       call WaveFunction_buildHCoreMatrix(wfObjects(speciesID))
    end do

  end subroutine MultiSCF_buildHcore

  !>
  !! @brief create hcore operators and use them to obtain a guess HC=eSC for the SCF equations, or read previous coefficients
  subroutine MultiSCF_getInitialGuess(this,wfObjects)
    type(MultiSCF) :: this  
    type(WaveFunction) :: wfObjects(*)

    integer :: speciesID, otherSpeciesID, i
    real(8) :: expectedOccupation,normCheck

    !!**********************************************************
    !! Build Guess and first density matrix
    !!
    do speciesID = 1, this%molSys%numberOfQuantumSpecies
       call DensityMatrixSCFGuess_getGuess( speciesID, wfObjects(speciesID)%HcoreMatrix, &
            wfObjects(speciesID)%transformationMatrix, &
            wfObjects(speciesID)%densityMatrix,&
            wfObjects(speciesID)%waveFunctionCoefficients, &
            this%printSCFiterations, &
            this%molSys)
       normCheck=sum( transpose(wfObjects(speciesID)%densityMatrix%values)*wfObjects(speciesID)%overlapMatrix%values)

       if ( this%printSCFiterations ) &
            write(*,"(A15,A10,A40,F12.6)") "number of ", trim(MolecularSystem_getNameOfSpecies(speciesID,this%molSys)) , &
            " particles in guess density matrix: ", normCheck

       expectedOccupation=MolecularSystem_getEta(speciesID,this%molSys)*this%molSys%species(speciesID)%ocupationNumber
       if (trim(MolecularSystem_getNameOfSpecies(speciesID,this%molSys)) .eq. trim(CONTROL_instance%IONIZE_SPECIES(1))) then
          do i=1,size(CONTROL_instance%IONIZE_MO)
             if(CONTROL_instance%IONIZE_MO(i) .gt. 0 .and. CONTROL_instance%MO_FRACTION_OCCUPATION(i) .lt. 1.0_8) &
                  expectedOccupation=expectedOccupation-MolecularSystem_getEta(speciesID,this%molSys)*(1.0-CONTROL_instance%MO_FRACTION_OCCUPATION(i))
          end do
       end if
       
       if (abs(expectedOccupation/normCheck-1.0) .gt. 1.0E-3 ) then
          ! wfObjects(speciesID)%densityMatrix%values=wfObjects(speciesID)%densityMatrix%values*expectedOccupation/normCheck
          if ( this%printSCFiterations ) &
               write(*,"(A65,F12.6)") "Warning!, density matrix with deviation from number of particles", expectedOccupation
               ! write(*,"(A65,F12.6)") "renormalized density matrix with deviation from number of particles", &
               ! sum( transpose(wfObjects(speciesID)%densityMatrix%values)*wfObjects(speciesID)%overlapMatrix%values)
       end if
       
       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "Initial Density Matrix ", trim(MolecularSystem_getNameOfSpecies(speciesID,this%molSys))
          call Matrix_show(wfObjects(speciesID)%densityMatrix)
       end if

    end do

    !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
    if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
         (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
       speciesID=MolecularSystem_getSpecieIDFromSymbol(trim("E-ALPHA"),this%molSys)
       otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol(trim("E-BETA"),this%molSys)

       if(MolecularSystem_getNumberOfParticles(speciesID,this%molSys) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID,this%molSys)) then
          wfObjects(otherSpeciesID)%waveFunctionCoefficients%values= wfObjects(speciesID)%waveFunctionCoefficients%values
          wfObjects(otherSpeciesID)%densityMatrix%values= wfObjects(speciesID)%densityMatrix%values
       end if

       if ( this%printSCFiterations ) &
            print *, "E-ALPHA AND E-BETA COEFFICIENTS ARE FORCED TO BE EQUAL IN THIS RUN"
    end if

  end subroutine MultiSCF_getInitialGuess

  !>
  !! @brief solve multcomponent FC=eSC SCF equations, store the coefficients in wfObjects, use the libint2Objects to compute the integrals in direct calculations
  subroutine MultiSCF_solveHartreeFockRoothan(this,wfObjects,libint2Objects)
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)
    type(Libint2Interface) :: libint2Objects(:)

    real(8) :: oldEnergy
    real(8) :: deltaEnergy
    integer :: numberOfSpecies
    integer :: wfnUnit, densUnit
    character(50) :: wfnFile, densFile
    character(50) :: integralsFile
    character(100) :: convergenceMessage
    integer :: integralsUnit

    logical :: GLOBAL_SCF_CONTINUE

    !! Open file for wfn
    wfnUnit = 300
    wfnFile = "lowdin.wfn"

    integralsUnit = 30
    integralsFile = "lowdin.opints"

    densUnit=78
    densFile="lowdin.densmatrix"

    numberOfSpecies = this%molSys%numberOfQuantumSpecies

    !!
    !!***************************************************************************************************************
    !! Begin Multi-species SCF

    if ( this%printSCFiterations ) then
       write(*,"(A)") "INFO: RUNNING SCHEME "//trim(this%name)
       write(*,"(A)")" "
       write(*,*) "Begin Multi-Species SCF calculation:"
       write(*,*) ""
       write(*,*) "--------------------------------------------------------------------------------------------------"
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          write(*,"(A15,A25,A25,A25,A25)") "Iteration", "Energy","Energy Change","Density Change","ParticlesInGrid" 
       else
          write(*,"(A15,A25,A25,A25)") "Iteration", "Energy","Energy Change","Density Change"
       end if
       write(*,*) "--------------------------------------------------------------------------------------------------"
    end if

    oldEnergy=0.0_8
    deltaEnergy=1.0E16_8
    this%totalDensityMatrixStandardDeviation=1.0E16_8
    GLOBAL_SCF_CONTINUE=.true.

    do while(GLOBAL_SCF_CONTINUE)

       call MultiSCF_iterate(this, wfObjects, libint2Objects(:), CONTROL_instance%ITERATION_SCHEME )
       deltaEnergy = oldEnergy -MultiSCF_getLastEnergy(this)
       oldEnergy = MultiSCF_getLastEnergy(this)

!!!!Print iteration results
       if ( this%printSCFiterations ) then
          if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
             write(*,"(I15,F25.12,F25.12,F25.12,F25.12)") MultiSCF_getNumberOfIterations(this), &
                  MultiSCF_getLastEnergy(this), deltaEnergy, &
                  this%totalDensityMatrixStandardDeviation ,&
                  sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%particlesInGrid)
          else
             write(*,"(I15,F25.12,F25.12,F25.12)") MultiSCF_getNumberOfIterations(this), &
                  MultiSCF_getLastEnergy(this), deltaEnergy, &
                  this%totalDensityMatrixStandardDeviation
          end if
       end if

!!!!Check convergence and print messages
       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "energyContinue", abs(deltaEnergy) .gt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE
          print *, "densityContinue",  this%totalDensityMatrixStandardDeviation .gt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE
          print *, "iterationContinue", MultiSCF_getNumberOfIterations(this) .lt. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS
       end if

       if(MultiSCF_getNumberOfIterations(this) .ge. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS) then
          write(convergenceMessage,"(A,I4,A)")  "The number of Iterations was exceded, the convergence had failed after", MultiSCF_getNumberOfIterations(this), "global iterations"
          GLOBAL_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "DENSITY" .and. &
            this%totalDensityMatrixStandardDeviation .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
          write(convergenceMessage,"(A,I4,A)") "Total density converged after", MultiSCF_getNumberOfIterations(this) ," global iterations"
          GLOBAL_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "ENERGY" .and. &
            abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE) then
          write(convergenceMessage,"(A,I4,A)") "Total energy converged after", MultiSCF_getNumberOfIterations(this) ," global iterations"
          GLOBAL_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "BOTH" .and. & 
            abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE .and. &
            this%totalDensityMatrixStandardDeviation .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
          write(convergenceMessage,"(A,I4,A)") "Total energy and density converged after", MultiSCF_getNumberOfIterations(this) ," global iterations"
          GLOBAL_SCF_CONTINUE=.false.
       end if

    end do

    if ( this%printSCFiterations ) then
       print *, convergenceMessage       
       print *,""
       print *,"...end Multi-Species SCF calculation"
       print *,""
    end if
    !!**************************************************************************************************************

    !! Multi-species SCF if HPG was instanced
    ! if (CONTROL_instance%HARTREE_PRODUCT_GUESS) then

    !    CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = .true.

    !    if ( this%printSCFiterations ) then
    !       print *,""
    !       print *,"Begin Second Multi-Species SCF calculation:"
    !       print *,""
    !       print *,"---------------------------------------------------------"
    !       write(*,"(A10,A12,A25)") "Iteration", "Energy","Energy Change"
    !       print *,"---------------------------------------------------------"
    !    end if

    !    print *,""
    !    print *,"...end Second Multi-Species SCF calculation"
    !    print *,""
    ! end if

    call MultiSCF_obtainFinalEnergy(this,wfObjects,libint2Objects(:))
    
  end subroutine MultiSCF_solveHartreeFockRoothan

    !>
  !! @brief solve multcomponent FC=eSC SCF equations, store the coefficients in wfObjects, use the libint2Objects to compute the integrals in direct calculations
  subroutine MultiSCF_obtainFinalEnergy(this,wfObjects,libint2Objects,method)
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)
    type(Libint2Interface) :: libint2Objects(:)
    character(*), optional :: method
    
    integer :: numberOfSpecies
    integer :: wfnUnit, densUnit
    integer :: speciesID, otherSpeciesID, i
    character(50) :: wfnFile, densFile
    character(30) :: labels(2)
    character(50) :: integralsFile
    integer :: integralsUnit

    if( .not. present(method) ) method=CONTROL_instance%METHOD
    
    !! Open file for wfn
    wfnUnit = 300
    wfnFile = "lowdin.wfn"

    integralsUnit = 30
    integralsFile = "lowdin.opints"

    densUnit=78
    densFile="lowdin.densmatrix"

    numberOfSpecies = this%molSys%numberOfQuantumSpecies

    if (CONTROL_instance%LOCALIZE_ORBITALS) then

       !! write coefficients and orbitals required in fchk files
       open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
       rewind(wfnUnit)
       do speciesID=1, numberOfSpecies
          labels(2) = MolecularSystem_getNameOfSpecies(speciesID,this%molSys)
          labels(1) = "COEFFICIENTS"
          call Matrix_writeToFile(wfObjects(speciesID)%waveFunctionCoefficients, unit=wfnUnit, binary=.true., arguments = labels )
          labels(1) = "ORBITALS"
          call Vector_writeToFile(wfObjects(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )
       end do
       close(wfnUnit)

       !! Build fchk files with Lowdin results
       call system("lowdin-output.x FCHK")

       !! Erkale Orbital Localization calls
       !! Orbital localization should not change density matrices
       write(*,*) ""
       print *, " ERKALE ORBITAL LOCALIZATION"
       write(*,*) "=============================="
       write(*,*) ""
       do speciesID=1, numberOfSpecies
          if(MolecularSystem_getMass(speciesID,this%molSys) .lt. 10.0 .and. MolecularSystem_getOcupationNumber(speciesID,this%molSys) .gt. 1) then !We assume that heavy particle orbitals are naturally localized
             call OrbitalLocalizer_erkaleLocal(speciesID,&
                  wfObjects( speciesID )%densityMatrix,&
                  wfObjects( speciesID )%fockMatrix, &
                  wfObjects( speciesID )%waveFunctionCoefficients, &
                  wfObjects( speciesID )%molecularOrbitalsEnergy)
          end if
       end do
    end if

    
    ! Final energy evaluation - larger integration grid for DFT
    do speciesID=1, numberOfSpecies
       call WaveFunction_buildDensityMatrix(wfObjects(speciesID))
    end do

    ! Reorder coefficients after creating the last density matrix
    do i=1,size(CONTROL_instance%IONIZE_MO)
       if(CONTROL_instance%IONIZE_MO(i) .gt. 0 .and. CONTROL_instance%MO_FRACTION_OCCUPATION(i) .eq. 0.0_8) then
          call MultiSCF_reorderIonizedCoefficients(this,wfObjects)
          exit
       end if
    end do

    !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
    if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
         (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
       speciesID=MolecularSystem_getSpecieIDFromSymbol(trim("E-ALPHA"),this%molSys)
       otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol(trim("E-BETA"),this%molSys)

       if(MolecularSystem_getNumberOfParticles(speciesID,this%molSys) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID,this%molSys)) then
          wfObjects(otherSpeciesID)%waveFunctionCoefficients%values= wfObjects(speciesID)%waveFunctionCoefficients%values
          wfObjects(otherSpeciesID)%densityMatrix%values= wfObjects(speciesID)%densityMatrix%values
       end if

       print *, "E-ALPHA AND E-BETA COEFFICIENTS ARE FORCED TO BE EQUAL IN THIS RUN"
    end if

    if (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") then
       call WaveFunction_writeDensityMatricesToFile(wfObjects, densFile)
       if (CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          call system("lowdin-DFT.x BUILD_FINAL_GRID "//trim(densFile))
          call system("lowdin-DFT.x FINAL_DFT "//trim(densFile))
       else
          call WaveFunction_getDFTContributions(wfObjects,this%DFTGrids,this%DFTGridsCommonPoints,"FINAL")
       end if
    end if

    do speciesID = 1, numberOfSpecies

       if (CONTROL_instance%COSMO) then
          call WaveFunction_buildCosmo2Matrix(wfObjects(speciesID))
          call WaveFunction_buildCosmoCoupling(wfObjects(speciesID))
       end if

       if ( (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") .and. CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          call WaveFunction_readExchangeCorrelationMatrix(wfObjects(speciesID))
       end if

       call WaveFunction_buildTwoParticlesMatrix(wfObjects(speciesID),libint2Objects=libint2Objects(:))
       !Separate coulomb and exchange contributions to two particles matrix

       call WaveFunction_buildTwoParticlesMatrix(wfObjects(speciesID), &
            twoParticlesMatrixOUT=wfObjects(speciesID)%hartreeMatrix(speciesID), factorIN=0.0_8, libint2Objects=libint2Objects(:) )     

       wfObjects(speciesID)%exchangeHFMatrix%values= wfObjects(speciesID)%twoParticlesMatrix%values &
            -wfObjects(speciesID)%hartreeMatrix(speciesID)%values

       call WaveFunction_buildCouplingMatrix(wfObjects,speciesID, libint2Objects=libint2Objects(:))       

       call WaveFunction_buildFockMatrix(wfObjects(speciesID))

       !!Obtain energy components for species
       call WaveFunction_obtainEnergyComponentsForSpecies(wfObjects(speciesID))

    end do

    !! Obtain energy compotents for whole system
    call WaveFunction_obtainTotalEnergy(wfObjects,&
         this%totalEnergy, &
         this%totalCouplingEnergy, &
         this%cosmo3Energy)

  end subroutine MultiSCF_obtainFinalEnergy

  subroutine MultiSCF_showResults(this,wfObjects)
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)

    integer :: i,j
    integer :: speciesID, otherSpeciesID
    integer :: numberOfContractions
    integer :: numberOfIterations
    real(8) :: diisError
    real(8) :: totalHartreeEnergy
    real(8) :: totalExchangeHFEnergy
    real(8) :: totalExchangeCorrelationEnergy
    real(8) :: totalQuantumPuntualInteractionEnergy
    real(8) :: totalExternalPotentialEnergy
    real(8) :: puntualInteractionEnergy
    real(8) :: puntualMMInteractionEnergy
    real(8) :: totalCosmoEnergy
    character :: convergenceType
    character(30) :: nameOfSpecies
    type(Matrix) :: coefficientsShow

    !! Show results
    !! Shows iterations by species

    if ( this%printSCFiterations ) then

       if(.not. CONTROL_instance%ELECTRONIC_WaveFunction_ANALYSIS ) then

          do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies(this%molSys)

             nameOfSpecies =  MolecularSystem_getNameOfSpecies(speciesID,this%molSys)                 
             numberOfIterations = List_size( wfObjects(speciesID)%energySCF )

             call List_begin( wfObjects(speciesID)%energySCF )
             call List_begin( wfObjects(speciesID)%diisError )
             call List_begin( wfObjects(speciesID)%standardDesviationOfDensityMatrixElements )

             print *,""
             print *,"Begin SCF calculation by: ",trim(nameOfSpecies)
             print *,"-------------------------"
             print *,""
             print *,"-----------------------------------------------------------------------------------"
             write(*,"(A10,A25,A25,A20)") "Iteration", "Energy", " Density Change","DIIS Error "
             print *,"-----------------------------------------------------------------------------------"

             do i=1, numberOfIterations-1

                call List_iterate( wfObjects(speciesID)%energySCF )
                call List_iterate( wfObjects(speciesID)%standardDesviationOfDensityMatrixElements )
                call List_iterate( wfObjects(speciesID)%diisError )
                diisError = List_current( wfObjects(speciesID)%diisError )

                convergenceType = ""

                if ( diisError > CONTROL_instance%DIIS_SWITCH_THRESHOLD ) convergenceType = "*"

                if (abs(diisError) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                   write(*,"(I5,F25.12,F25.12,A19,A1)") i,  List_current( wfObjects(speciesID)%energySCF ),&
                        List_current( wfObjects(speciesID)%standardDesviationOfDensityMatrixElements ), &
                        "--",convergenceType
                else
                   write(*,"(I5,F25.12,F25.12,F25.12,A1)") i,  List_current( wfObjects(speciesID)%energySCF ),&
                        List_current( wfObjects(speciesID)%standardDesviationOfDensityMatrixElements ), &
                        diisError,convergenceType
                end if

             end do
             print *,""
             print *,"... end SCF calculation"

          end do

       end if

    end if

    write(*,*) ""
    write(*,*) " EIGENVALUES AND EIGENVECTORS: "
    write(*,*) "=============================="
    write(*,*) ""

    if ( CONTROL_instance%HF_PRINT_EIGENVALUES ) then
       do speciesID = 1, this%molSys%numberOfQuantumSpecies                
          write(*,*) ""
          write(*,*) " Eigenvalues for: ", trim( this%molSys%species(speciesID)%name )
          write(*,*) "-----------------"
          write(*,*) ""
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID,this%molSys)
          do i = 1 , numberOfContractions 
             write(6,"(T2,I4,F25.12)") i,wfObjects(speciesID)%molecularOrbitalsEnergy%values(i)
          end do
          write(*,*) ""
       end do
       write(*,*) " end of eigenvalues "
    end if

    if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL" .or. trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then

       do speciesID = 1, this%molSys%numberOfQuantumSpecies      

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID,this%molSys)

          if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL") then

             write(*,*) ""
             write(*,*) " Eigenvectors for: ", trim( this%molSys%species(speciesID)%name )
             write(*,*) "-----------------"
             write(*,*) ""

             call Matrix_constructor(coefficientsShow,int(numberOfContractions,8),&
                  int(numberOfContractions-wfObjects(speciesID)%removedOrbitals,8),0.0_8)
             do i=1, numberOfContractions
                do j=1, numberOfContractions-wfObjects(speciesID)%removedOrbitals
                   coefficientsShow%values(i,j)=wfObjects(speciesID)%waveFunctionCoefficients%values(i,j)
                end do
             end do

          else if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then

             write(*,*) ""
             write(*,*) " Occupied Eigenvectors for: ", trim( this%molSys%species(speciesID)%name )
             write(*,*) "--------------------------- "
             write(*,*) ""

             call Matrix_constructor(coefficientsShow,int(numberOfContractions,8),int(MolecularSystem_getOcupationNumber(speciesID,this%molSys),8),0.0_8)
             do i=1, numberOfContractions
                do j=1, MolecularSystem_getOcupationNumber(speciesID,this%molSys)
                   coefficientsShow%values(i,j)=wfObjects(speciesID)%waveFunctionCoefficients%values(i,j)
                end do
             end do

          end if

          call Matrix_show(coefficientsShow , &
               rowkeys = MolecularSystem_getlabelsofcontractions(speciesID,this%molSys), &
               columnkeys = string_convertvectorofrealstostring(wfObjects(speciesID)%molecularOrbitalsEnergy ),&
               flags=WITH_BOTH_KEYS)

          call Matrix_destructor(coefficientsShow)
       end do

       write(*,*) ""
       write(*,*) " end of eigenvectors "

    end if

    write(*,*) ""
    write(*,*) " END OF EIGENVALUES AND EIGENVECTORS"
    write(*,*) ""

    !!Shows Energy components
    write(*,*) ""             
    write(*,*) " COMPONENTS OF KINETIC ENERGY: "
    write(*,*) "-----------------------------"
    write(*,*) ""             

    do speciesID = 1, this%molSys%numberOfQuantumSpecies                
       write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
            " Kinetic energy = ", wfObjects(speciesID)%kineticEnergy
    end do
    this%totalKineticEnergy = sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%kineticEnergy)             

    write(*,"(T38,A25)") "___________________________"
    write(*,"(A38,F25.12)") "Total kinetic energy = ", this%totalKineticEnergy

    write(*,*) ""
    write(*,*) " COMPONENTS OF POTENTIAL ENERGY: "
    write(*,*) "-------------------------------"
    write(*,*) ""

    puntualInteractionEnergy = MolecularSystem_getPointChargesEnergy(this%molSys)
    write(*,"(A38,F25.12)") "Fixed potential energy    = ", puntualInteractionEnergy

    puntualMMInteractionEnergy = MolecularSystem_getMMPointChargesEnergy(this%molSys)
    if(CONTROL_instance%CHARGES_MM) then
       write(*,"(A38,F25.12)") "Self MM potential energy   = ", puntualMMInteractionEnergy
    end if

    write(*,*) ""
    write(*,*) " Quantum/Fixed interaction energy: "
    write(*,*) "----------------------------------"
    write(*,*) ""

    do speciesID = 1, this%molSys%numberOfQuantumSpecies                
       write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
            "/Fixed interact. energy = ", wfObjects(speciesID)%puntualInteractionEnergy
    end do
    totalQuantumPuntualInteractionEnergy = sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%puntualInteractionEnergy )
    write(*,"(T38,A25)") "___________________________"
    write(*,"(A38,F25.12)") "Total Q/Fixed energy = ", totalQuantumPuntualInteractionEnergy

    write(*,*) ""
    write(*,*) " Coulomb energy: "
    write(*,*) "------------------"
    write(*,*) ""
    totalHartreeEnergy=0.0
    do speciesID = 1, this%molSys%numberOfQuantumSpecies                
       write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
            "/"//trim( this%molSys%species(speciesID)%name ) // &
            " Hartree energy = ", wfObjects(speciesID)%hartreeEnergy(speciesID)
       totalHartreeEnergy=totalHartreeEnergy+wfObjects(speciesID)%hartreeEnergy(speciesID)
    end do
    do speciesID = 1, this%molSys%numberOfQuantumSpecies                
       do otherSpeciesID = speciesID + 1, this%molSys%numberOfQuantumSpecies                
          write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
               "/"//trim( this%molSys%species(otherSpeciesID)%name ) // &
               " Hartree energy = ", wfObjects(speciesID)%hartreeEnergy(otherSpeciesID)
          totalHartreeEnergy=totalHartreeEnergy+wfObjects(speciesID)%hartreeEnergy(otherSpeciesID)
       end do
    end do
    write(*,"(T38,A25)") "___________________________"
    write(*,"(A38,F25.12)") "Total Hartree energy = ", totalHartreeEnergy

    write(*,*) ""
    write(*,*) " Exchange(HF) energy: "
    write(*,*) "----------------------"
    write(*,*) ""
    do speciesID = 1, this%molSys%numberOfQuantumSpecies                
       write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
            " Exchange energy = ", wfObjects(speciesID)%exchangeHFEnergy
    end do
    totalExchangeHFEnergy=sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%exchangeHFEnergy)
    write(*,"(T38,A25)") "___________________________"
    write(*,"(A38,F25.12)") "Total Exchange energy = ", totalExchangeHFEnergy


    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       write(*,*) ""
       write(*,*) " Exchange-Correlation(DFT) energy: "
       write(*,*) "-----------------------------------"
       write(*,*) "" 
       totalExchangeCorrelationEnergy=0.0
       do speciesID = 1, this%molSys%numberOfQuantumSpecies                
          write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
               " Exc.Corr. energy = ", wfObjects(speciesID)%exchangeCorrelationEnergy(speciesID)
          totalExchangeCorrelationEnergy=totalExchangeCorrelationEnergy+wfObjects(speciesID)%exchangeCorrelationEnergy(speciesID)
       end do
       do speciesID = 1, this%molSys%numberOfQuantumSpecies                
          do otherSpeciesID = speciesID + 1, this%molSys%numberOfQuantumSpecies                
             write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name ) // &
                  "/"//trim( this%molSys%species(otherSpeciesID)%name ) // &
                  " Corr. energy = ", wfObjects(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)
             totalExchangeCorrelationEnergy=totalExchangeCorrelationEnergy+wfObjects(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)
          end do
       end do
       write(*,"(T38,A25)") "___________________________"
       write(*,"(A38,F25.12)") "Total Exchange Correlation energy = ", totalExchangeCorrelationEnergy
    end if


    if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

       write(*,*) ""
       write(*,*) " External Potential energy: "
       write(*,*) "----------------------------"
       write(*,*) ""

       do speciesID = 1, this%molSys%numberOfQuantumSpecies                
          write(*,"(A38,F25.12)") trim( this%molSys%species(speciesID)%name) // &
               " Ext Pot energy = ", wfObjects(speciesID)%externalPotentialEnergy
       end do
       totalExternalPotentialEnergy=sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%externalPotentialEnergy)
       write(*,"(T38,A25)") "___________________________"
       write(*,"(A38,F25.12)") "Total External Potential energy = ", totalExternalPotentialEnergy             

    end if

    if(CONTROL_instance%COSMO) then
       write(*,*) ""
       write(*,*) " COSMO ENERGY: "
       write(*,*) "--------------"
       write(*,*) ""
       totalCosmoEnergy = sum(wfObjects(1:this%molSys%numberOfQuantumSpecies)%cosmoEnergy)
       write(*,"(A38,F25.12)") "Total Cosmo Energy    = ", totalCosmoEnergy
       write(*,"(A38,F25.12)") "Cosmo 3 Energy    = ", this%cosmo3Energy
    end if

    this%totalPotentialEnergy = puntualInteractionEnergy &
         + totalQuantumPuntualInteractionEnergy &
         + totalHartreeEnergy &
         + totalExchangeHFEnergy &
         + totalExchangeCorrelationEnergy &
         + totalExternalPotentialEnergy &
         + totalCosmoEnergy &
         + this%cosmo3Energy

    write(*,*) ""
    write(*,*) " TOTAL ENERGY COMPONENTS: "
    write(*,*) "=========================="
    write(*,*) ""
    write(*,"(A38,F25.12)") "TOTAL KINETIC ENERGY      = ", this%totalKineticEnergy
    write(*,"(A38,F25.12)") "TOTAL POTENTIAL ENERGY    = ", this%totalPotentialEnergy
    write(*,"(T38,A25)") "___________________________"
    write(*,"(A38,F25.12)") "TOTAL ENERGY = ", this%totalEnergy             
    write(*,*) ""
    write(*,"(A38,F25.12)") "VIRIAL RATIO (V/T) = ", - ( this%totalPotentialEnergy / this%totalKineticEnergy)
    write(*,*) ""
    write(*,*) ""
    write(*,*) " END ENERGY COMPONENTS"
    write(*,*) ""  


    if(CONTROL_instance%COSMO) then
       write(*,*) ""
       write(*,*) " COSMO CHARGE: "
       write(*,*) "--------------"
       write(*,*) ""
       call WaveFunction_cosmoQuantumCharge(this%molSys)
    end if

  end subroutine MultiSCF_showResults

  !!**********************************************************
  !! Save matrices to lowdin.wfn file
  !!
  subroutine MultiSCF_saveWfn(this,wfObjects)
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)

    integer :: numberOfSpecies
    integer :: wfnUnit, vecUnit, densUnit
    integer :: speciesID
    character(50) :: wfnFile, vecFile, densFile
    character(30) :: labels(2)

    character(50) :: integralsFile
    integer :: integralsUnit

    !! Open file for wfn
    wfnUnit = 300
    wfnFile = "lowdin.wfn"

    integralsUnit = 30
    integralsFile = "lowdin.opints"

    densUnit=78
    densFile="lowdin.densmatrix"


    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    rewind(wfnUnit)

    labels = ""

    numberOfSpecies = this%molSys%numberOfQuantumSpecies

    do speciesID = 1, numberOfSpecies

       labels(2) = MolecularSystem_getNameOfSpecies(speciesID,this%molSys)

       labels(1) = "REMOVED-ORBITALS"
       call Vector_writeToFile(unit=wfnUnit, binary=.true., value=real(wfObjects(speciesID)%removedOrbitals,8), arguments= labels )

       labels(1) = "TWOPARTICLES"
       call Matrix_writeToFile(wfObjects(speciesID)%twoParticlesMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       labels(1) = "COUPLING"
       call Matrix_writeToFile(wfObjects(speciesID)%couplingMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       labels(1) = "EXCHANGE-CORRELATION"
       call Matrix_writeToFile(wfObjects(speciesID)%exchangeCorrelationMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       labels(1) = "EXCHANGE-CORRELATION-ENERGY"
       call Vector_writeToFile(unit=wfnUnit, binary=.true., value=sum(wfObjects(speciesID)%exchangeCorrelationEnergy(:)), arguments= labels )

       labels(1) = "COEFFICIENTS"
       call Matrix_writeToFile(wfObjects(speciesID)%waveFunctionCoefficients, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "DENSITY"
       call Matrix_writeToFile(wfObjects(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "KINETIC"
       call Matrix_writeToFile(wfObjects(speciesID)%kineticMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "ATTRACTION"
       call Matrix_writeToFile(wfObjects(speciesID)%puntualInteractionMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "EXTERNAL"
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
            call Matrix_writeToFile(wfObjects(speciesID)%externalPotentialMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "HCORE"
       call Matrix_writeToFile(wfObjects(speciesID)%hcoreMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "ORBITALS"
       call Vector_writeToFile(wfObjects(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "FOCK"
       call Matrix_writeToFile(wfObjects(speciesID)%fockMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "OVERLAP"
       call Matrix_writeToFile(wfObjects(speciesID)%overlapMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "TRANSFORMATION"
       call Matrix_writeToFile(wfObjects(speciesID)%transformationMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       if (CONTROL_instance%COSMO) then
          labels(1) = "COSMO1"
          call Matrix_writeToFile(wfObjects(speciesID)%cosmo1, unit=wfnUnit, binary=.true., arguments = labels )
          labels(1) = "COSMO2"
          call Matrix_writeToFile(wfObjects(speciesID)%cosmo2, unit=wfnUnit, binary=.true., arguments = labels )  
          labels(1) = "COSMOCOUPLING"
          call Matrix_writeToFile(wfObjects(speciesID)%cosmoCoupling, unit=wfnUnit, binary=.true., arguments = labels ) 
          labels(1) = "COSMO4"
          call Matrix_writeToFile(wfObjects(speciesID)%cosmo4, unit=wfnUnit, binary=.true., arguments = labels )
       end if

    end do

    labels = ""
    !! Open file for vec
    vecUnit = 36
    if ( CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then
       vecFile = trim(CONTROL_instance%INPUT_FILE)//"vec"
       open(unit=vecUnit, file=trim(vecFile), form="unformatted", status='replace')
       do speciesID = 1, numberOfSpecies
          labels(2) = MolecularSystem_getNameOfSpecies(speciesID,this%molSys)
          labels(1) = "COEFFICIENTS"
          call Matrix_writeToFile(wfObjects(speciesID)%waveFunctionCoefficients, &
               unit=vecUnit, binary=.true., arguments = labels)

          labels(1) = "ORBITALS"
          call Vector_writeToFile(wfObjects(speciesID)%molecularOrbitalsEnergy, & 
               unit=vecUnit, binary=.true., arguments = labels )
       end do

    else
       vecFile = trim(CONTROL_instance%INPUT_FILE)//"plainvec"
       open(unit=vecUnit, file=trim(vecFile), form="formatted", status='replace')

       do speciesID = 1, numberOfSpecies
          labels(2) = MolecularSystem_getNameOfSpecies(speciesID,this%molSys)
          labels(1) = "COEFFICIENTS"
          call Matrix_writeToFile(wfObjects(speciesID)%waveFunctionCoefficients, &
               unit=vecUnit, binary=.false., arguments = labels)

          labels(1) = "ORBITALS"
          call Vector_writeToFile(wfObjects(speciesID)%molecularOrbitalsEnergy, & 
               unit=vecUnit, binary=.false., arguments = labels )
       end do

    end if
    close (vecUnit)


    !!**********************************************************
    !! Save Some energies
    !!
    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=this%totalEnergy, arguments=["TOTALENERGY"])

    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=this%cosmo3Energy, arguments=["COSMO3ENERGY"])

    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=this%totalCouplingEnergy, arguments=["COUPLINGENERGY"])

    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MolecularSystem_getPointChargesEnergy(this%molSys), arguments=["PUNTUALINTERACTIONENERGY"])

    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=- ( this%totalPotentialEnergy / this%totalKineticEnergy) , arguments=["VIRIAL"])

    close(wfnUnit)

  end subroutine MultiSCF_saveWfn

!!Reorders full ionized orbitals
  subroutine MultiSCF_reorderIonizedCoefficients(this,wfObjects)
    implicit none
    type(MultiSCF) :: this
    type(WaveFunction) :: wfObjects(*)
    type(Matrix) :: auxMatrix
    type(Vector) :: auxVector
    integer :: occupationNumber, newOccupationNumber, i, j, speciesID

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies(this%molSys)
       if (trim(wfObjects(speciesID)%name) .eq. trim(CONTROL_instance%IONIZE_SPECIES(1)) ) then
          occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molSys)
          newOccupationNumber=occupationNumber
          call Matrix_copyConstructor(auxMatrix,wfObjects(speciesID)%waveFunctionCoefficients)
          call Vector_copyConstructor(auxVector,wfObjects(speciesID)%molecularOrbitalsEnergy)
          do i= 1, size(CONTROL_instance%IONIZE_MO)
             if(CONTROL_instance%IONIZE_MO(i) .gt. 0 .and. CONTROL_instance%MO_FRACTION_OCCUPATION(i) .eq. 0.0_8 ) then
                wfObjects(speciesID)%waveFunctionCoefficients%values(:,occupationNumber)=auxMatrix%values(:,CONTROL_instance%IONIZE_MO(i))
                wfObjects(speciesID)%molecularOrbitalsEnergy%values(occupationNumber)=auxVector%values(CONTROL_instance%IONIZE_MO(i))
                do j=CONTROL_instance%IONIZE_MO(i),occupationNumber-1
                   wfObjects(speciesID)%waveFunctionCoefficients%values(:,j)=auxMatrix%values(:,j+1)
                   wfObjects(speciesID)%molecularOrbitalsEnergy%values(j)=auxVector%values(j+1)
                end do
                newOccupationNumber=newOccupationNumber-1
             end if             
          end do
          this%molSys%species(speciesID)%ocupationNumber=newOccupationNumber
          if(CONTROL_instance%DEBUG_SCFS) then
             print *, "newOccupationNumber for",  trim(wfObjects(speciesID)%name), this%molSys%species(speciesID)%ocupationNumber
             call Matrix_show(wfObjects(speciesID)%waveFunctionCoefficients)
          end if
       end if
    end do
       
  end subroutine MultiSCF_reorderIonizedCoefficients
  
  
end module MultiSCF_

