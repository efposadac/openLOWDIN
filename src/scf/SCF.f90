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

!> @brief This program handles booth single and multi species SCF procedure in the APMO approach 
!! @author E. F. Posada, 2013.
!! @info This code is based on the code of S. A. Gonzalez (APMO, 2010), but mostly code has been completly rewritten to fit 
!!       new LOWDIN coding standard. 2013
!! @info All iterations schemes have been tested. But if you want to add new matrix, have to fix all shchemes 
!!       to support this matrix.
!! This program needs lowdincore library to compile, all functions of molecular system are extensively used.
program SCF
  use Stopwatch_
  use CONTROL_ 
  use WaveFunction_
  use MolecularSystem_
  use MultiSCF_
  use DensityMatrixSCFGuess_
  use String_
  use Exception_
  use omp_lib
  use OrbitalLocalizer_
  implicit none

  real(8) :: oldEnergy
  real(8) :: deltaEnergy
  real(8) :: diisError
  integer :: i,j
  integer :: numberOfSpecies
  integer :: status
  integer :: wfnUnit, vecUnit, densUnit
  integer :: speciesID, otherSpeciesID
  integer :: numberOfIterations
  character(50) :: wfnFile, vecFile, densFile
  character(30) :: nameOfSpecies
  character(30) :: labels(2)
  character(100) :: iterationScheme(0:4)
  character :: convergenceType

  type(Matrix) :: auxDensity
  type(Matrix) :: coefficientsShow
  character(50) :: job
  character(50) :: integralsFile
  character(50) :: arguments(20)
  integer :: integralsUnit
  logical :: existFile
  integer :: numberOfContractions
  real(8) :: totalEnergy
  real(8) :: totalHartreeEnergy
  real(8) :: totalExchangeHFEnergy
  real(8) :: totalExchangeCorrelationEnergy
  real(8) :: totalKineticEnergy
  real(8) :: totalRepulsionEnergy
  real(8) :: totalQuantumPuntualInteractionEnergy
  real(8) :: totalExternalPotentialEnergy
  real(8) :: electronicRepulsionEnergy
  real(8) :: puntualInteractionEnergy
  real(8) :: puntualMMInteractionEnergy
  real(8) :: potentialEnergy

  !!cosmo things
  character(50) :: cosmoIntegralsFile
  real(8) :: totalCosmoEnergy
  real(8) :: cosmo3Energy

  logical :: GLOBAL_SCF_CONTINUE

  
  job = ""
  call get_command_argument(1,value=job)
  job = trim(String_getUppercase(job))

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  write(*,"(A)")"--------------------------------------------------------------------------------"
  write(*,"(A)")"** PROGRAM: SCF (Hartree Fock).     Author: S.A. Gonzalez, E. Posada, F. Moncada"
  write(*,"(A)")"--------------------------------------------------------------------------------"
  write(*,"(A)") "INFO: RUNNING IN "//trim(job)//" MODE."
  write(*,"(A)")" "

  if(.not. CONTROL_instance%FIRST_STEP) then
     write (6,"(T20,A30)") " TEST GEOMETRY: AMSTRONG"
     write (6,"(T18,A35)") "------------------------------------------"
     call MolecularSystem_showCartesianMatrix()
     call MolecularSystem_showDistanceMatrix()
  end if


  !! Open file for wfn
  wfnUnit = 300
  wfnFile = "lowdin.wfn"

  integralsUnit = 30
  integralsFile = "lowdin.opints"

  densUnit=78
  densFile="lowdin.densmatrix"

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !! Start the MultiSCF object
  call MultiSCF_constructor(CONTROL_instance%ITERATION_SCHEME)

  !! Start the wavefunction object
  call WaveFunction_constructor( )

  !! Start the orbital localizer object
  if (CONTROL_instance%LOCALIZE_ORBITALS) call OrbitalLocalizer_constructor( )

  numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

  !****************************************************************************************************
  !! Builds the fock operator
  !!

  !! Calculate one-particle integrals  
  call system("lowdin-ints.x ONE_PARTICLE")

  !! Check the one-particle integrals file  
  existFile = .false.     
  inquire(file=trim(integralsFile), exist=existFile)

  if( existFile ) then
     open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

     read(integralsUnit) numberOfSpecies

     if(MolecularSystem_instance%numberOfQuantumSpecies /= numberOfSpecies ) then

        call MolecularSystem_exception( ERROR, "Bad "//trim(integralsFile)//" file!", "In SCF.f90 at main program")

     end if

     close(integralsUnit)

  else

     call MolecularSystem_exception(ERROR,"lowdin.opints file not found!", "In SCF.f90 at main program")

  end if

  open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")

  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

     !!**********************************************************
     !! Builds Hcore
     !!
     !! Overlap Matrix
     call WaveFunction_buildOverlapMatrix(trim(integralsFile), speciesID)

     !! Transformation Matrix
     call WaveFunction_buildTransformationMatrix( trim(integralsFile), speciesID, 2 )

     if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
        call WaveFunction_buildExternalPotentialMatrix(trim(integralsFile), speciesID)
     end if

     !! Hcore Matrix
     call WaveFunction_HCoreMatrix(trim(integralsFile), speciesID)

     !!**********************************************************
     !! Builds Cosmo hcore integrals
     !!
     !! 
     if(CONTROL_instance%COSMO)then
        cosmoIntegralsFile="cosmo.opints"
        call WaveFunction_cosmoHcoreMatrix(trim(cosmoIntegralsFile), speciesID)
     end if

     !!**********************************************************
     !! Build Guess and first density matrix
     !!

     call DensityMatrixSCFGuess_getGuess( speciesID, WaveFunction_instance(speciesID)%densityMatrix, WaveFunction_instance(speciesID)%waveFunctionCoefficients )

     write(*,"(A,A,A,F7.3)") "number of ", trim(MolecularSystem_getNameOfSpecie( speciesID )) ," particles in guess density matrix: ",  &
         sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values)*WaveFunction_instance(speciesID)%overlapMatrix%values)
  end do
    
  !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
  if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
     (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
     speciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-ALPHA")  )
     otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-BETA")  )

     if(MolecularSystem_getNumberOfParticles(speciesID) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID) ) then
        WaveFunction_instance(otherSpeciesID)%waveFunctionCoefficients%values= WaveFunction_instance(speciesID)%waveFunctionCoefficients%values
        WaveFunction_instance(otherSpeciesID)%densityMatrix%values= WaveFunction_instance(speciesID)%densityMatrix%values
     end if

     print *, "E-ALPHA AND E-BETA COEFFICIENTS ARE FORCED TO BE EQUAL IN THIS RUN"
  end if
   
  !!**********************************************************
  !! Save matrices to lowdin.wfn file required by ints program
  !!
  
  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
     labels(2) = MolecularSystem_getNameOfSpecie(speciesID)

     labels(1) = "DENSITY"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )

     if ( CONTROL_instance%DEBUG_SCFS ) then
        print *, "Initial Density Matrix ", trim(MolecularSystem_getNameOfSpecie( speciesID ))
        call Matrix_show(WaveFunction_instance(speciesID)%densityMatrix)
     end if
  end do
  close(wfnUnit)
     

  !!**************************************************************************************************************************
  !! Calculate two-particle integrals (not building 2 particles and coupling matrix... those matrices updated at each SCF cycle)
  !!

  if( CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL ) then

     call system(" lowdin-ints.x TWO_PARTICLE_G12")

  else        

     call system(" lowdin-ints.x TWO_PARTICLE_R12")

  end if

  !!
  !!***************************************************************************************************************
 
  !! Begin SCF calculation...  

  write(*,"(A)") "INFO: RUNNING SCHEME "//trim(MultiSCF_instance%name)
  write(*,"(A)")" "

  !!**************************************************************************************************************
  !! Multi-species SCF

  if ( .not. CONTROL_instance%OPTIMIZE ) then
     write(*,*) "Begin Multi-Species SCF calculation:"
     write(*,*) ""
     write(*,*) "-------------------------------------------------------------------------"
     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        write(*,"(A20,A12,A20,A20,A20)") "Iteration", "Energy","Energy Change","Density Change","ParticlesInGrid" 
     else
        write(*,"(A20,A12,A20,A20)") "Iteration", "Energy","Energy Change","Density Change"
     end if
     write(*,*) "-------------------------------------------------------------------------"
  end if

  oldEnergy=0.0_8
  deltaEnergy=1.0E16_8
  MultiSCF_instance%totalDensityMatrixStandardDeviation=1.0E16_8
  GLOBAL_SCF_CONTINUE=.true.
  
  do while(GLOBAL_SCF_CONTINUE)
          
     call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )
     deltaEnergy = oldEnergy -MultiSCF_getLastEnergy()
     oldEnergy = MultiSCF_getLastEnergy()

     !!!!Print iteration results
     if ( .not.CONTROL_instance%OPTIMIZE ) then
        if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
           write (*,"(I15,F20.12,F20.12,F20.12,F20.12)") MultiSCF_getNumberOfIterations(), &
                MultiSCF_getLastEnergy(), deltaEnergy, &
                MultiSCF_instance%totalDensityMatrixStandardDeviation ,&
                sum(WaveFunction_instance(:)%particlesInGrid)
        else
           write (*,"(I15,F20.12,F20.12,F20.12)") MultiSCF_getNumberOfIterations(), &
                MultiSCF_getLastEnergy(), deltaEnergy, &
                MultiSCF_instance%totalDensityMatrixStandardDeviation
        end if
     end if
     
     !!!!Check convergence and print messages
     if ( CONTROL_instance%DEBUG_SCFS ) then
        print *, "energyContinue", abs(deltaEnergy) .gt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE, "densityContinue", &
             MultiSCF_instance%totalDensityMatrixStandardDeviation .gt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE, &
             "iterationContinue", MultiSCF_getNumberOfIterations() .lt. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS
     end if

     if(MultiSCF_getNumberOfIterations() .ge. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS) then
        write (*,"(A,I4,A)")  "The number of Iterations was exceded, the convergence had failed after", MultiSCF_getNumberOfIterations() ," global iterations"
        GLOBAL_SCF_CONTINUE=.false.
     end if
          
     if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "DENSITY" .and. &
          MultiSCF_instance%totalDensityMatrixStandardDeviation .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
        write (*,"(A,I4,A)") "Total density converged after", MultiSCF_getNumberOfIterations() ," global iterations"
        GLOBAL_SCF_CONTINUE=.false.
     end if
     
     if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "ENERGY" .and. &
          abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE) then
        write (*,"(A,I4,A)") "Total energy converged after", MultiSCF_getNumberOfIterations() ," global iterations"
        GLOBAL_SCF_CONTINUE=.false.
     end if

     if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "BOTH" .and. & 
          abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE .and. &
          MultiSCF_instance%totalDensityMatrixStandardDeviation .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
        write (*,"(A,I4,A)") "Total energy and density converged after", MultiSCF_getNumberOfIterations() ," global iterations"
        GLOBAL_SCF_CONTINUE=.false.
     end if
     
  end do
  
  print *,""
  print *,"...end Multi-Species SCF calculation"
  print *,""
  !!**************************************************************************************************************

  !! Multi-species SCF if HPG was instanced
  if (CONTROL_instance%HARTREE_PRODUCT_GUESS) then

     CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = .true.

     if ( .not. CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
        print *,""
        print *,"Begin Second Multi-Species SCF calculation:"
        print *,""
        print *,"---------------------------------------------------------"
        write (*,"(A10,A12,A25)") "Iteration", "Energy","Energy Change"
        print *,"---------------------------------------------------------"
     end if

     !     do i = 1, numberOfSpecies
     !        nameOfSpecies = MolecularSystem_getNameOfSpecie(i)
     !        call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies), MultiSCF_instance%nproc )
     !        call WaveFunction_buildFockMatrix( trim(nameOfSpecies) )
     !
     !        if (CONTROL_instance%COSMO) then
     !           call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecies))
     !        end if
     !     end do

     ! do while( status == SCF_GLOBAL_CONVERGENCE_CONTINUE .and. &
     !      MultiSCF_getNumberOfIterations() <= CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS )

     !    call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )
     !    deltaEnergy = oldEnergy-MultiSCF_getLastEnergy()

     !    if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
     !       write (6,"(I5,F20.12,F20.12)") MultiSCF_getNumberOfIterations(), &
     !            MultiSCF_getLastEnergy(),deltaEnergy
     !    end if

     !    ! status = MultiSCF_testEnergyChange(CONTROL_instance%TOTAL_ENERGY_TOLERANCE  )
     !    oldEnergy=MultiSCF_getLastEnergy()

     ! end do

     print *,""
     print *,"...end Second Multi-Species SCF calculation"
     print *,""

  end if

  !! Shows iterations by species
  if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then

     if(.not. CONTROL_instance%ELECTRONIC_WaveFunction_ANALYSIS ) then

        do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

           nameOfSpecies =  MolecularSystem_getNameOfSpecie(speciesID)                 
           numberOfIterations = List_size( WaveFunction_instance(speciesID)%energySCF )

           call List_begin( WaveFunction_instance(speciesID)%energySCF )
           call List_begin( WaveFunction_instance(speciesID)%diisError )
           call List_begin( WaveFunction_instance(speciesID)%standardDesviationOfDensityMatrixElements )

           print *,""
           print *,"Begin SCF calculation by: ",trim(nameOfSpecies)
           print *,"-------------------------"
           print *,""
           print *,"-----------------------------------------------------------------"
           write (*,"(A10,A12,A25,A20)") "Iteration", "Energy", " Density Change","         DIIS Error "
           print *,"-----------------------------------------------------------------"

           do i=1, numberOfIterations-1

              call List_iterate( WaveFunction_instance(speciesID)%energySCF )
              call List_iterate( WaveFunction_instance(speciesID)%standardDesviationOfDensityMatrixElements )
              call List_iterate( WaveFunction_instance(speciesID)%diisError )
              diisError = List_current( WaveFunction_instance(speciesID)%diisError )

              convergenceType = ""

              if ( diisError > CONTROL_instance%DIIS_SWITCH_THRESHOLD ) convergenceType = "*"

              if (abs(diisError) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                 write (6,"(I5,F20.12,F20.12,A20,A1)") i,  List_current( WaveFunction_instance(speciesID)%energySCF ),&
                      List_current( WaveFunction_instance(speciesID)%standardDesviationOfDensityMatrixElements ), &
                      "         --         ",convergenceType
              else
                 write (6,"(I5,F20.12,F20.12,F20.12,A1)") i,  List_current( WaveFunction_instance(speciesID)%energySCF ),&
                      List_current( WaveFunction_instance(speciesID)%standardDesviationOfDensityMatrixElements ), &
                      diisError,convergenceType
              end if

           end do
           print *,""
           print *,"... end SCF calculation"

        end do

     end if

  end if

  ! else

  !    call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )


  ! end if


  if (CONTROL_instance%LOCALIZE_ORBITALS) then

     !! write coefficients and orbitals required in fchk files
     open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
     rewind(wfnUnit)
     do speciesID=1, numberOfSpecies
        labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
        labels(1) = "COEFFICIENTS"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, unit=wfnUnit, binary=.true., arguments = labels )
        labels(1) = "ORBITALS"
        call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )
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
        if(MolecularSystem_getMass( speciesID ) .lt. 10.0 .and. MolecularSystem_getOcupationNumber( speciesID ) .gt. 1) then !We assume that heavy particle orbitals are naturally localized
           call OrbitalLocalizer_erkaleLocal(speciesID,&
                WaveFunction_instance( speciesID )%densityMatrix,&
                WaveFunction_instance( speciesID )%fockMatrix, &
                WaveFunction_instance( speciesID )%waveFunctionCoefficients, &
                WaveFunction_instance( speciesID )%molecularOrbitalsEnergy)
        end if
     end do
  end if

  !! Show results
  if(CONTROL_instance%LAST_STEP) then
     write(*,*) ""
     write(*,*) " EIGENVALUES AND EIGENVECTORS: "
     write(*,*) "=============================="
     write(*,*) ""

     if ( CONTROL_instance%HF_PRINT_EIGENVALUES ) then
        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
           write(*,*) ""
           write(*,*) " Eigenvalues for: ", trim( MolecularSystem_instance%species(speciesID)%name )
           write(*,*) "-----------------"
           write(*,*) ""
           numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
           do i = 1 , numberOfContractions 
              write(6,"(T2,I4,F20.12)") i,WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)
           end do
           write(*,*) ""
        end do
        write(*,*) " end of eigenvalues "
     end if

     if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL" .or. trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then

        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies      

           numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
          
           if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL") then

              write(*,*) ""
              write(*,*) " Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
              write(*,*) "-----------------"
              write(*,*) ""
              
              call Matrix_constructor(coefficientsShow,int(numberOfContractions,8),int(numberOfContractions-WaveFunction_instance(speciesID)%removedOrbitals,8),0.0_8)
              do i=1, numberOfContractions
                 do j=1, numberOfContractions-WaveFunction_instance(speciesID)%removedOrbitals
                    coefficientsShow%values(i,j)=WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(i,j)
                 end do
              end do
              
           else if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then
              
              write(*,*) ""
              write(*,*) " Occupied Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
              write(*,*) "--------------------------- "
              write(*,*) ""
              
              call Matrix_constructor(coefficientsShow,int(numberOfContractions,8),int(MolecularSystem_getOcupationNumber(speciesID),8),0.0_8)
              do i=1, numberOfContractions
                 do j=1, MolecularSystem_getOcupationNumber(speciesID)
                    coefficientsShow%values(i,j)=WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(i,j)
                 end do
              end do
              
           end if

           call Matrix_show(coefficientsShow , &
                rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
                columnkeys = string_convertvectorofrealstostring( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy ),&
                flags=WITH_BOTH_KEYS)

           call Matrix_destructor(coefficientsShow)
        end do
        
        write(*,*) ""
        write(*,*) " end of eigenvectors "

     end if

     write(*,*) ""
     write(*,*) " END OF EIGENVALUES AND EIGENVECTORS"
     write(*,*) ""
  end if


  ! Final energy evaluation - larger integration grid for DFT
  do i=1, numberOfSpecies
     nameOfSpecies = MolecularSystem_getNameOfSpecie(i)
     call WaveFunction_buildDensityMatrix( trim(nameOfSpecies) )
  end do

  !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
  if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
     (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
     speciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-ALPHA")  )
     otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-BETA")  )

     if(MolecularSystem_getNumberOfParticles(speciesID) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID) ) then
        WaveFunction_instance(otherSpeciesID)%waveFunctionCoefficients%values= WaveFunction_instance(speciesID)%waveFunctionCoefficients%values
        WaveFunction_instance(otherSpeciesID)%densityMatrix%values= WaveFunction_instance(speciesID)%densityMatrix%values
     end if

     print *, "E-ALPHA AND E-BETA COEFFICIENTS ARE FORCED TO BE EQUAL IN THIS RUN"
  end if
  
  if (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") then
     call system("lowdin-DFT.x BUILD_FINAL_GRID "//trim(densFile))
     call WaveFunction_writeDensityMatricesToFile(trim(densFile))
     call system("lowdin-DFT.x FINAL_DFT "//trim(densFile))
  end if

  do speciesID = 1, numberOfSpecies
     nameOfSpecies = MolecularSystem_getNameOfSpecie(speciesID)

     if (CONTROL_instance%COSMO) then
        call WaveFunction_buildCosmo2Matrix( trim(nameOfSpecies) )
        call WaveFunction_buildCosmoCoupling( trim(nameOfSpecies) )
     end if

     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        call WaveFunction_readExchangeCorrelationMatrix( trim(nameOfSpecies))
     end if

     call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies) )
     !Separate coulomb and exchange contributions to two particles matrix

     call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies), &
          twoParticlesMatrixOUT=WaveFunction_instance(speciesID)%hartreeMatrix(speciesID), factorIN=0.0_8 )     

     WaveFunction_instance(speciesID)%exchangeHFMatrix%values= WaveFunction_instance(speciesID)%twoParticlesMatrix%values &
          -WaveFunction_instance(speciesID)%hartreeMatrix(speciesID)%values

     call WaveFunction_buildCouplingMatrix( trim(nameOfSpecies) )       

     call WaveFunction_buildFockMatrix( trim(nameOfSpecies) )

     !!Obtain energy components for species
     call WaveFunction_obtainEnergyComponents(speciesID)

  end do

  !! Obtain energy compotents for whole system
  call WaveFunction_obtainTotalEnergy(&
       MultiSCF_instance%totalEnergy, &
       MultiSCF_instance%totalCouplingEnergy, &
       MultiSCF_instance%electronicRepulsionEnergy, &
       MultiSCF_instance%cosmo3Energy)


  !!Shows Energy components
  write(*,*) ""             
  write(*,*) " COMPONENTS OF KINETIC ENERGY: "
  write(*,*) "-----------------------------"
  write(*,*) ""             

  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
          " Kinetic energy = ", WaveFunction_instance(speciesID)%kineticEnergy
  end do
  totalKineticEnergy = sum( WaveFunction_instance(:)%kineticEnergy)             

  write (6,"(T10,A50)") "_____________________"
  write (6,"(A38,F20.12)") "Total kinetic energy = ", totalKineticEnergy

  write(*,*) ""
  write(*,*) " COMPONENTS OF POTENTIAL ENERGY: "
  write(*,*) "-------------------------------"
  write(*,*) ""

  write(*,*) ""
  write(*,*) " Quantum/Fixed interaction energy: "
  write(*,*) "----------------------------------"
  write(*,*) ""

  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
          "/Fixed interact. energy = ", WaveFunction_instance(speciesID)%puntualInteractionEnergy
  end do
  totalQuantumPuntualInteractionEnergy = sum ( WaveFunction_instance(:)%puntualInteractionEnergy )
  write (6,"(T10,A50)") "_____________________"
  write (6,"(A38,F20.12)") "Total Q/Fixed energy = ", totalQuantumPuntualInteractionEnergy

  write(*,*) ""
  write(*,*) " Coulomb energy: "
  write(*,*) "------------------"
  write(*,*) ""
  totalHartreeEnergy=0.0
  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
          "/"//trim( MolecularSystem_instance%species(speciesID)%name ) // &
          " Hartree energy = ", WaveFunction_instance(speciesID)%hartreeEnergy(speciesID)
     totalHartreeEnergy=totalHartreeEnergy+WaveFunction_instance(speciesID)%hartreeEnergy(speciesID)
  end do
  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     do otherSpeciesID = speciesID + 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
             "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
             " Hartree energy = ", WaveFunction_instance(speciesID)%hartreeEnergy(otherSpeciesID)
        totalHartreeEnergy=totalHartreeEnergy+WaveFunction_instance(speciesID)%hartreeEnergy(otherSpeciesID)
     end do
  end do
  write (6,"(T10,A50)") "_____________________"
  write (6,"(A38,F20.12)") "Total Hartree energy = ", totalHartreeEnergy

  write(*,*) ""
  write(*,*) " Exchange(HF) energy: "
  write(*,*) "----------------------"
  write(*,*) ""
  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
          " Exchange energy = ", WaveFunction_instance(speciesID)%exchangeHFEnergy
  end do
  totalExchangeHFEnergy=sum(WaveFunction_instance(:)%exchangeHFEnergy)
  write (6,"(T10,A50)") "_____________________"
  write (6,"(A38,F20.12)") "Total Exchange energy = ", totalExchangeHFEnergy


  if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
     write(*,*) ""
     write(*,*) " Exchange-Correlation(DFT) energy: "
     write(*,*) "-----------------------------------"
     write(*,*) "" 
     totalExchangeCorrelationEnergy=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
             " Exc.Corr. energy = ", WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(speciesID)
        totalExchangeCorrelationEnergy=totalExchangeCorrelationEnergy+WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(speciesID)
     end do
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        do otherSpeciesID = speciesID + 1, MolecularSystem_instance%numberOfQuantumSpecies                
           write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                " Corr. energy = ", WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)
           totalExchangeCorrelationEnergy=totalExchangeCorrelationEnergy+WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)
        end do
     end do
     write (6,"(T10,A50)") "_____________________"
     write (6,"(A38,F20.12)") "Total Exchange Correlation energy = ", totalExchangeCorrelationEnergy
  end if


  if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

     write(*,*) ""
     write(*,*) " External Potential energy: "
     write(*,*) "----------------"
     write(*,*) ""

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A38,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name) // &
             " Ext Pot energy = ", WaveFunction_instance(speciesID)%externalPotentialEnergy
     end do
     totalExternalPotentialEnergy=sum(WaveFunction_instance(:)%externalPotentialEnergy)
     write (6,"(A38,F20.12)") "Total External Potential energy = ", totalExternalPotentialEnergy             

  end if

  write(*,*) ""
  puntualInteractionEnergy = MolecularSystem_getPointChargesEnergy()
  write (6,"(A38,F20.12)") "Fixed potential energy    = ", puntualInteractionEnergy

  puntualMMInteractionEnergy = MolecularSystem_getMMPointChargesEnergy()
  if(CONTROL_instance%CHARGES_MM) then
     write (6,"(A38,F20.12)") "Self MM potential energy   = ", puntualMMInteractionEnergy
  end if

  potentialEnergy = puntualInteractionEnergy &
       + totalQuantumPuntualInteractionEnergy &
       + totalHartreeEnergy &
       + totalExchangeHFEnergy &
       + totalExchangeCorrelationEnergy &
       + totalExternalPotentialEnergy

  totalCosmoEnergy = sum( WaveFunction_instance(:)%cosmoEnergy)

  if(CONTROL_instance%COSMO) then
     write(*,*)"totalCosmoEnergy",totalCosmoEnergy
     write(*,*)"cosmo3energy",cosmo3Energy

     potentialEnergy=potentialEnergy+totalCosmoEnergy+cosmo3Energy

  end if

  write(*,*) ""
  write(*,*) " TOTAL ENERGY COMPONENTS: "
  write(*,*) "=========================="
  write(*,*) ""
  write (6,"(A38,F20.12)") "TOTAL KINETIC ENERGY      = ", totalKineticEnergy
  write (6,"(A38,F20.12)") "TOTAL POTENTIAL ENERGY    = ", potentialEnergy
  write (6,"(T10,A50)") "_____________________"
  write (6,"(A38,F20.12)") "TOTAL ENERGY = ", MultiSCF_instance%totalEnergy             
  write(*,*) ""
  write (6,"(A38,F20.12)") "VIRIAL RATIO (V/T) = ", - ( potentialEnergy / totalKineticEnergy)
  write(*,*) ""
  write(*,*) ""
  write(*,*) " END ENERGY COMPONENTS"
  write(*,*) ""


  !stop time
  ! call Stopwatch_stop(lowdin_stopwatch)



  !!**********************************************************
  !! Save matrices to lowdin.wfn file
  !!
  open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
  rewind(wfnUnit)

  labels = ""



  do speciesID = 1, numberOfSpecies

     labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
     
     labels(1) = "REMOVED-ORBITALS"
     call Vector_writeToFile(unit=wfnUnit, binary=.true., value=real(WaveFunction_instance(speciesID)%removedOrbitals,8), arguments= labels )

     labels(1) = "TWOPARTICLES"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%twoParticlesMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "COUPLING"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%couplingMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "EXCHANGE-CORRELATION"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%exchangeCorrelationMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "EXCHANGE-CORRELATION-ENERGY"
     call Vector_writeToFile(unit=wfnUnit, binary=.true., value=sum(WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(:)), arguments= labels )

     labels(1) = "COEFFICIENTS"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, unit=wfnUnit, binary=.true., arguments = labels )

     labels(1) = "DENSITY"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )

     labels(1) = "HCORE"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%hcoreMatrix, unit=wfnUnit, binary=.true., arguments = labels )

     labels(1) = "ORBITALS"
     call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )

     labels(1) = "FOCK"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%fockMatrix, unit=wfnUnit, binary=.true., arguments = labels )

     labels(1) = "OVERLAP"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%overlapMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     labels(1) = "TRANSFORMATION"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%transformationMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )
     
     if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
        labels(1) = "EXTERNAL_POTENTIAL"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%externalPotentialMatrix, unit=wfnUnit, binary=.true., arguments = labels )
     end if

     if (CONTROL_instance%COSMO) then
        labels(1) = "COSMO1"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo1, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )
        labels(1) = "COSMO2"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo2, unit=wfnUnit, binary=.true., arguments = labels )  
        labels(1) = "COSMOCOUPLING"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmoCoupling, unit=wfnUnit, binary=.true., arguments = labels ) 
        labels(1) = "COSMO4"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo4, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )
     end if

  end do

  labels = ""
  !! Open file for vec
  vecUnit = 36
  if ( CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then
     vecFile = trim(CONTROL_instance%INPUT_FILE)//"vec"
     open(unit=vecUnit, file=trim(vecFile), form="unformatted", status='replace')
     do speciesID = 1, numberOfSpecies
        labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
        labels(1) = "COEFFICIENTS"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
             unit=vecUnit, binary=.true., arguments = labels)

        labels(1) = "ORBITALS"
        call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, & 
             unit=vecUnit, binary=.true., arguments = labels )
     end do

  else
     vecFile = trim(CONTROL_instance%INPUT_FILE)//"plainvec"
     open(unit=vecUnit, file=trim(vecFile), form="formatted", status='replace')

     do speciesID = 1, numberOfSpecies
        labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
        labels(1) = "COEFFICIENTS"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
             unit=vecUnit, binary=.false., arguments = labels)

        labels(1) = "ORBITALS"
        call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, & 
             unit=vecUnit, binary=.false., arguments = labels )
     end do

  end if
  close (vecUnit)


  !!**********************************************************
  !! Save Some energies
  !!
  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalEnergy, arguments=["TOTALENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%cosmo3Energy, arguments=["COSMO3ENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalCouplingEnergy, arguments=["COUPLINGENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%electronicRepulsionEnergy, arguments=["COUPLING-E-"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MolecularSystem_getPointChargesEnergy(), arguments=["PUNTUALINTERACTIONENERGY"])
  
  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=- ( potentialEnergy / totalKineticEnergy) , arguments=["VIRIAL"])

  close(wfnUnit)
  
  !stop time
  call Stopwatch_stop(lowdin_stopwatch)

  if(CONTROL_instance%LAST_STEP) then
     write(*, *) ""
     write(*,"(A,F10.3,A4)") "** TOTAL CPU Time SCF : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time SCF : ", lowdin_stopwatch%elapsetWTime ," (s)"
     write(*, *) ""
  end if


  if( .not.CONTROL_instance%OPTIMIZE .and. CONTROL_instance%GET_GRADIENTS ) then        

     call system("lowdin-ints.x GET_GRADIENTS")

  end if

  
  if (CONTROL_instance%COSMO) then
     call WaveFunction_cosmoQuantumCharge()
  end if


  !!calculate HF/KS properties
  call system ("lowdin-CalcProp.x")

  if (CONTROL_instance%SUBSYSTEM_EMBEDDING) then

     print *, ""
     print *, "-------------------------------------------"
     print *, "STARTING SUBSYSTEM EMBEDDED SCF CALCULATION"
     print *, "-------------------------------------------"
     print *, ""

     !! Asign orbitals to fragments and creates the subsystem fock matrix
     call OrbitalLocalizer_levelShiftSubSystemOrbitals()

  end if

end program SCF



  !   vecUnit = 36
  !   if ( CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then

  !      open(unit=vecUnit, file=trim(vecFile), status="replace", form="unformatted")

  !      do speciesID = 1, numberOfSpecies

  !         labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
  !         labels(1) = "COEFFICIENTS"
  !         call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
  !              unit=vecUnit, binary=.true., arguments = labels)

  !      end do

  !      close (vecUnit)

  !   else
  !      labels = ""
  !      !! Open file for wfn
  !      vecFile = trim(CONTROL_instance%INPUT_FILE)//"lowdin-plain.vec"
  !      open(unit=vecUnit, file=trim(vecFile), status="replace", form="formatted")

  !      if ( .not. CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY .or. .not. CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then

  !         labels = ""
  !     !! Open file for vec
  !     vecUnit = 36
  !     vecFile = "lowdin-plain.vec"
  !     open(unit=vecUnit, file=trim(vecFile), form="formatted", status='unknown')

  !     if ( .not. CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then
  !       do speciesID = 1, numberOfSpecies

  !         labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
  !         labels(1) = "COEFFICIENTS"
  !         call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
  !              unit=vecUnit, binary=.false., arguments = labels)
  !       end do
  !     end if

