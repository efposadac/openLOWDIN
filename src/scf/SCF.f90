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

  real(8) :: auxValue
  real(8) :: deltaEnergy
  real(8) :: diisError
  integer :: i
  integer :: numberOfSpecies
  integer :: status
  integer :: wfnUnit, vecUnit, densUnit
  integer :: speciesID
  integer :: numberOfIterations
  character(50) :: wfnFile, vecFile, densFile
  character(30) :: nameOfSpecies
  character(30) :: labels(2)
  character(100) :: iterationScheme(0:4)
  character :: convergenceType

  type(Matrix) :: auxDensity
  character(50) :: job
  character(50) :: integralsFile
  character(50) :: arguments(20)
  integer :: integralsUnit
  logical :: existFile
  integer :: numberOfContractions
  real(8) :: totalEnergy
  real(8) :: totalCouplingEnergy
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

  job = ""
  call get_command_argument(1,value=job)
  job = trim(String_getUppercase(job))

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  if(CONTROL_instance%LAST_STEP) then
     write(*,"(A)")"----------------------------------------------------------------------"
     write(*,"(A)")"** PROGRAM: SCF (Hartree Fock).      Author: S.A. Gonzalez, E. Posada  "
     write(*,"(A)")"----------------------------------------------------------------------"

     write(*,"(A)") "INFO: RUNNING IN "//trim(job)//" MODE."
     write(*,"(A)")" "
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
  
  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  print *, "integralsFile trololoooooooooooooooooooooooooo", integralsFile

  !! Start the MultiSCF object
  call MultiSCF_constructor()

  print *, "integralsFile trololooooooooooo", integralsFile
  !! Start the wavefunction object
  call WaveFunction_constructor( )

  print *, "integralsFile trololoooooo", integralsFile
  !! Start the orbital localizer object
  call OrbitalLocalizer_constructor( )

  numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

  !****************************************************************************************************
  !! Builds the fock operator
  !!

  print *, "integralsFile trololo", integralsFile

  !! Calculate one-particle integrals  
  call system("lowdin-ints.x ONE_PARTICLE")

  !! Check the one-particle integrals file  
  existFile = .false.     
  inquire(file=trim(integralsFile), exist=existFile)

  print *, "integralsFile", integralsFile
  
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


  !! Open file for wavefunction
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
     if ( MolecularSystem_instance%species(speciesID)%isElectron ) then

        auxDensity=DensityMatrixSCFGuess_getGuess( CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS, speciesID )

        call WaveFunction_setDensityMatrix(  auxDensity, speciesID )                 
        call Matrix_destructor(auxDensity)

     else


        auxDensity=DensityMatrixSCFGuess_getGuess( CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS, speciesID )

        call WaveFunction_setDensityMatrix(  auxDensity, speciesID )
        call Matrix_destructor(auxDensity)


     end if


     !!**********************************************************
     !! Save matrices to lowdin.wfn file
     !!

     arguments = ""
     arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

     arguments(1) = "OVERLAP"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%overlapMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     arguments(1) = "HCORE"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%HcoreMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     arguments(1) = "DENSITY"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     arguments(1) = "TRANSFORMATION"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%transformationMatrix, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
       arguments(1) = "EXTERNAL_POTENTIAL"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%externalPotentialMatrix, unit=wfnUnit, binary=.true., &
               arguments = arguments(1:2) )
     end if

     if(CONTROL_instance%COSMO)then

        arguments(1) = "COSMO1"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo1, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

        arguments(1) = "COSMO4"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo4, unit=wfnUnit, binary=.true., arguments = arguments(1:2) )

     end if

  end do


  close(wfnUnit)


  !!**************************************************************************************************************************
  !! Calculate two-particle integrals (not building 2 particles and coupling matrix... those matrices are done by SCF program)
  !!

  if( CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL ) then

     call system(" lowdin-ints.x TWO_PARTICLE_G12")

  else        

     call system(" lowdin-ints.x TWO_PARTICLE_R12")

  end if

  if(CONTROL_instance%LAST_STEP) then
     write(*,*) "DONE!"
  end if

  !!
  !!***************************************************************************************************************

  !! Begin SCF calculation...  
  
  ! if( numberOfSpecies > 1 ) then

  iterationScheme(0) = "NONELECRONIC FULLY CONVERGED BY ELECTRONIC ITERATION"
  iterationScheme(1) = "ELECTRONIC FULLY CONVERGED BY NONELECRONIC ITERATION"
  iterationScheme(2) = "SPECIES FULLY CONVERGED INDIVIDIALLY"
  iterationScheme(3) = "SIMULTANEOUS"
  iterationScheme(4) = "SIMULTANEOUS_TRIAL"

  write(*,"(A)") "INFO: RUNNING "//trim(iterationScheme(CONTROL_instance%ITERATION_SCHEME))//" SCHEME."
  write(*,"(A)")" "

  !! Multi-species SCF
  status = SCF_GLOBAL_CONVERGENCE_CONTINUE      

  if ( .not. CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
     write(*,*) "Begin Multi-Species SCF calculation:"
     write(*,*) ""
     write(*,*) "---------------------------------------------------------"
     write (6,"(A10,A12,A25)") "Iteration", "Energy","Energy Change"
     write(*,*) "---------------------------------------------------------"
  end if

  if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
     !!Save density matrices to file for DFT calculations
     densUnit = 78
     densFile = trim(CONTROL_instance%INPUT_FILE)//"densmatrix"
     open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
     labels(1) = "DENSITY-MATRIX"
     do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
        labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
     end do
     close (densUnit)
     call system("lowdin-DFT.x BUILD_MATRICES "//trim(densFile))
  end if

  do i = 1, numberOfSpecies
     nameOfSpecies = MolecularSystem_getNameOfSpecie(i)

     call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies))

     call WaveFunction_buildCouplingMatrix( trim(nameOfSpecies))

     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        call WaveFunction_buildExchangeCorrelationMatrix( trim(nameOfSpecies))
     end if

     if (CONTROL_instance%COSMO) then
        call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecies))
     end if

     call WaveFunction_buildFockMatrix( trim(nameOfSpecies) )

  end do

  auxValue = 0.0_8 
  do while( status == SCF_GLOBAL_CONVERGENCE_CONTINUE .and. &
       MultiSCF_getNumberOfIterations() <= CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS )

     call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )
     deltaEnergy = auxValue-MultiSCF_getLastEnergy()

     if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
        write (6,"(I5,F20.10,F20.10)") MultiSCF_getNumberOfIterations(), &
             MultiSCF_getLastEnergy(),deltaEnergy
     end if

     status = MultiSCF_testEnergyChange(CONTROL_instance%TOTAL_ENERGY_TOLERANCE  )
     auxValue = MultiSCF_getLastEnergy()

  end do

  print *,""
  print *,"...end Multi-Species SCF calculation"
  print *,""

  !! Multi-species SCF if HPG was instanced
  if (CONTROL_instance%HARTREE_PRODUCT_GUESS) then

     status = SCF_GLOBAL_CONVERGENCE_CONTINUE

     CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = .true.

     if ( .not. CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
        print *,""
        print *,"Begin Second Multi-Species SCF calculation:"
        print *,""
        print *,"---------------------------------------------------------"
        write (6,"(A10,A12,A25)") "Iteration", "Energy","Energy Change"
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

     do while( status == SCF_GLOBAL_CONVERGENCE_CONTINUE .and. &
          MultiSCF_getNumberOfIterations() <= CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS )

        call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )
        deltaEnergy = auxValue-MultiSCF_getLastEnergy()

        if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
           write (6,"(I5,F20.10,F20.10)") MultiSCF_getNumberOfIterations(), &
                MultiSCF_getLastEnergy(),deltaEnergy
        end if

        status = MultiSCF_testEnergyChange(CONTROL_instance%TOTAL_ENERGY_TOLERANCE  )
        auxValue=MultiSCF_getLastEnergy()

     end do

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
           call List_begin( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )

           print *,""
           print *,"Begin SCF calculation by: ",trim(nameOfSpecies)
           print *,"-------------------------"
           print *,""
           print *,"-----------------------------------------------------------------"
           write (*,"(A10,A12,A25,A20)") "Iteration", "Energy", " Density Change","         DIIS Error "
           print *,"-----------------------------------------------------------------"

           do i=1, numberOfIterations-1

              call List_iterate( WaveFunction_instance(speciesID)%energySCF )
              call List_iterate( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )
              call List_iterate( WaveFunction_instance(speciesID)%diisError )
              diisError = List_current( WaveFunction_instance(speciesID)%diisError )

              convergenceType = ""

              if ( diisError > CONTROL_instance%DIIS_SWITCH_THRESHOLD ) convergenceType = "*"

              if (abs(diisError) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                 write (6,"(I5,F20.10,F20.10,A20,A1)") i,  List_current( WaveFunction_instance(speciesID)%energySCF ),&
                      List_current( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements ), &
                      "         --         ",convergenceType
              else
                 write (6,"(I5,F20.10,F20.10,F20.10,A1)") i,  List_current( WaveFunction_instance(speciesID)%energySCF ),&
                      List_current( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements ), &
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


  ! Final energy evaluation - larger integration grid 
  do i=1, numberOfSpecies
     nameOfSpecies = MolecularSystem_getNameOfSpecie(i)
     call SingleSCF_actualizeDensityMatrix( trim(nameOfSpecies) )
  end do
  
  if ( (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") .and. &
       ( CONTROL_instance%FINAL_GRID_ANGULAR_POINTS*CONTROL_instance%FINAL_GRID_RADIAL_POINTS  .gt. &
       CONTROL_instance%GRID_ANGULAR_POINTS*CONTROL_instance%GRID_RADIAL_POINTS ) ) then
     !!Save density matrices to file for DFT calculations
     densUnit = 78
     densFile = trim(CONTROL_instance%INPUT_FILE)//"densmatrix"
     open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
     labels(1) = "DENSITY-MATRIX"
     do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
        labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
     end do
     close (densUnit)
     call system("lowdin-DFT.x FINAL_GRID "//trim(densFile))
  end if

  do i = 1, numberOfSpecies
     nameOfSpecies = MolecularSystem_getNameOfSpecie(i)

     if (CONTROL_instance%COSMO) then
        call WaveFunction_buildCosmo2Matrix( trim(nameOfSpecies) )
        call WaveFunction_buildCosmoCoupling( trim(nameOfSpecies) )
     end if
     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        call WaveFunction_buildExchangeCorrelationMatrix( trim(nameOfSpecies))
     end if

     call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies) )
     call WaveFunction_buildCouplingMatrix( trim(nameOfSpecies) )       

     call WaveFunction_buildFockMatrix( trim(nameOfSpecies) )
     WaveFunction_instance(i)%wasBuiltFockMatrix = .true.

     !!Obtain energy components for species
     call WaveFunction_obtainEnergyComponents(i)

  end do

  !! Obtain energy compotents for whole system
  call WaveFunction_obtainTotalEnergy(&
       MultiSCF_instance%totalEnergy, &
       MultiSCF_instance%totalCouplingEnergy, &
       MultiSCF_instance%electronicRepulsionEnergy, &
       MultiSCF_instance%cosmo3Energy)

  totalKineticEnergy = sum( WaveFunction_instance(:)%kineticEnergy)             
  totalRepulsionEnergy = sum( WaveFunction_instance(:)%repulsionEnergy ) + electronicRepulsionEnergy                          
  totalQuantumPuntualInteractionEnergy = sum ( WaveFunction_instance(:)%puntualInteractionEnergy )
  totalCouplingEnergy=MultiSCF_instance%totalCouplingEnergy-electronicRepulsionEnergy !! sum(WaveFunction_instance(:)%couplingEnergy)/2.0
  totalExchangeCorrelationEnergy = sum ( WaveFunction_instance(:)%exchangeCorrelationEnergy )             
  totalExternalPotentialEnergy = sum ( WaveFunction_instance(:)%externalPotentialEnergy )             
  puntualInteractionEnergy = MolecularSystem_getPointChargesEnergy()
  puntualMMInteractionEnergy = MolecularSystem_getMMPointChargesEnergy()
  potentialEnergy = totalRepulsionEnergy &
       + puntualInteractionEnergy &
       + totalQuantumPuntualInteractionEnergy &
       + totalCouplingEnergy &
       + totalExchangeCorrelationEnergy &
       + totalExternalPotentialEnergy
  totalCosmoEnergy = sum( WaveFunction_instance(:)%cosmoEnergy)

  if(CONTROL_instance%COSMO) then
     write(*,*)"totalCosmoEnergy",totalCosmoEnergy
     write(*,*)"cosmo3energy",cosmo3Energy

     potentialEnergy=potentialEnergy+totalCosmoEnergy+cosmo3Energy

  end if

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
     print *, "ERKALE ORBITAL LOCALIZATION"
     do speciesID=1, numberOfSpecies
        print *, "density before"
        call OrbitalLocalizer_erkaleLocal(speciesID,&
             WaveFunction_instance( speciesID )%densityMatrix,&
             WaveFunction_instance( speciesID )%fockMatrix, &
             WaveFunction_instance( speciesID )%waveFunctionCoefficients, &
             WaveFunction_instance( speciesID )%molecularOrbitalsEnergy)
     end do
  end if
  
  !! Show results
  if(CONTROL_instance%LAST_STEP) then
     write(*,*) ""
     write(*,*) " EIGENVALUES AND EIGENVECTORS: "
     write(*,*) "=============================="
     write(*,*) ""

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies      

        write(*,*) ""
        write(*,*) " Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
        write(*,*) "-----------------"
        write(*,*) ""

        call Matrix_show( WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
             rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
             columnkeys = string_convertvectorofrealstostring( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy ),&
             flags=WITH_BOTH_KEYS)

        write(*,*) ""
        write(*,*) " end of eigenvectors "

     end do

     write(*,*) ""
     write(*,*) " END OF EIGENVALUES AND EIGENVECTORS"
     write(*,*) ""

     !!Shows Energy components
     write(*,*) ""
     write(*,*) " ENERGY COMPONENTS: "
     write(*,*) "=================="
     write(*,*) ""
     write (6,"(T10,A28,F20.12)") "TOTAL KINETIC ENERGY      = ", sum(WaveFunction_instance(:)%kineticEnergy)
     write (6,"(T10,A28,F20.12)") "TOTAL POTENTIAL ENERGY    = ", potentialEnergy
     write (6,"(T10,A50)") "________________"
     write (6,"(T10,A28,F20.12)") "TOTAL ENERGY = ", MultiSCF_instance%totalEnergy             
     write(*,*) ""
     write (6,"(T10,A28,F20.12)") "VIRIAL RATIO (V/T) = ", - ( potentialEnergy / totalKineticEnergy)
     write(*,*) ""             
     write(*,*) " COMPONENTS OF KINETIC ENERGY: "
     write(*,*) "-----------------------------"
     write(*,*) ""             

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

        write (6,"(T10,A8,A20,F20.10)") trim( MolecularSystem_instance%species(speciesID)%name ), &
             " Kinetic energy   = ", WaveFunction_instance(speciesID)%kineticEnergy
     end do

     write (6,"(T10,A50)") "________________"
     write (6,"(T10,A28,F20.10)") "Total kinetic energy = ", totalKineticEnergy

     write(*,*) ""
     write(*,*) " COMPONENTS OF POTENTIAL ENERGY: "
     write(*,*) "-------------------------------"
     write(*,*) ""
     write (6,"(T10,A30,F20.10)") "Fixed potential energy     = ", puntualInteractionEnergy
     if(CONTROL_instance%CHARGES_MM) then
     write (6,"(T10,A30,F20.10)") "Self MM potential energy   = ", puntualMMInteractionEnergy
     end if
     write (6,"(T10,A30,F20.10)") "Q/Fixed potential energy   = ", totalQuantumPuntualInteractionEnergy
     write (6,"(T10,A30,F20.10)") "Coupling energy            = ", totalCouplingEnergy
     write (6,"(T10,A30,F20.10)") "Repulsion energy           = ", totalRepulsionEnergy
     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
     write (6,"(T10,A30,F20.10)") "Exchange Correlation energy= ", totalExchangeCorrelationEnergy
     end if
     write (6,"(T10,A30,F20.10)") "External Potential energy  = ", totalExternalPotentialEnergy             
     write (6,"(T10,A50)") "________________"
     write (6,"(T10,A30,F20.10)") "Total potential energy     = ", potentialEnergy

     write(*,*) ""
     write(*,*) " Repulsion energy: "
     write(*,*) "------------------"
     write(*,*) ""

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

        write (6,"(T10,A26,A2,F20.10)") trim( MolecularSystem_instance%species(speciesID)%name ) // "/" // &
             trim(MolecularSystem_instance%species(speciesID)%name ) // &
             " Repulsion energy  ","= ", WaveFunction_instance(speciesID)%repulsionEnergy
     end do

     if(CONTROL_instance%IS_OPEN_SHELL) then

        write (6,"(T10,A26,A2,F20.10)") "E-ALPHA" // "/" // &
             "E-BETA" // " Repulsion energy  ","= ", electronicRepulsionEnergy
     end if

     write(*,*) ""
     write(*,*) " Coupling energy: "
     write(*,*) "----------------"
     write(*,*) ""

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

        write (6,"(T10,A26,A2,F20.10)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
             " Coupling energy  ","= ", WaveFunction_instance(speciesID)%couplingEnergy
     end do

     if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

        write(*,*) ""
        write(*,*) " External Potential energy: "
        write(*,*) "----------------"
        write(*,*) ""

        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

           write (6,"(T10,A26,A2,F20.10)") trim( MolecularSystem_instance%species(speciesID)%name) // &
                " Ext Pot energy  ","= ", WaveFunction_instance(speciesID)%externalPotentialEnergy
        end do

     end if

     write(*,*) ""
     write(*,*) " Quantum/Fixed interaction energy: "
     write(*,*) "-----------------------"
     write(*,*) ""

     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

        write (6,"(T10,A26,A2,F20.10)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
             "/Fixed interact. energy ","= ", WaveFunction_instance(speciesID)%puntualInteractionEnergy
     end do

     write(*,*) ""
     write(*,*) " END ENERGY COMPONENTS"
     write(*,*) ""


     if ( CONTROL_instance%HF_PRINT_EIGENVALUES ) then
       write(*,*) "BEGIN EIGENVALUES"
       write(*,*) ""
       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                

         write (6,"(T2,A12)") trim( MolecularSystem_instance%species(speciesID)%name) 
         write(*,*) ""

         numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

         do i = 1 , numberOfContractions 
           write(6,"(T2,I4,F20.10)") i,WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)
         end do
         write(*,*) ""
       end do
       write(*,*) "END OF EIGENVALUES"
    end if
  end if

  !!save virial
  open(unit=wfnUnit, file=trim(wfnFile), status="old", position="append", form="unformatted")

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=- ( potentialEnergy / totalKineticEnergy) , arguments=["VIRIAL"])

  close(wfnUnit)

  !stop time
  call Stopwatch_stop(lowdin_stopwatch)


  if( .not.CONTROL_instance%OPTIMIZE .and. CONTROL_instance%GET_GRADIENTS ) then        

     call system("lowdin-ints.x GET_GRADIENTS")

  end if
  
  
  !!**********************************************************
  !! Save matrices to lowdin.wfn file
  !!
  open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
  rewind(wfnUnit)

  labels = ""
  
  
  
  do speciesID = 1, numberOfSpecies

     labels(2) = MolecularSystem_getNameOfSpecie(speciesID)

     labels(1) = "TWOPARTICLES"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%twoParticlesMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "COUPLING"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%couplingMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "EXCHANGE-CORRELATION"
     call Matrix_writeToFile(WaveFunction_instance(speciesID)%exchangeCorrelationMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

     labels(1) = "EXCHANGE-CORRELATION-ENERGY"
     call Vector_writeToFile(unit=wfnUnit, binary=.true., value=WaveFunction_instance(speciesID)%exchangeCorrelationEnergy, arguments= labels )
     
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

     if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
       labels(1) = "EXTERNAL_POTENTIAL"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%externalPotentialMatrix, unit=wfnUnit, binary=.true., arguments = labels )
     end if

     if (CONTROL_instance%COSMO) then
        labels(1) = "COSMO2"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo2, unit=wfnUnit, binary=.true., arguments = labels )  
        labels(1) = "COSMOCOUPLING"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmoCoupling, unit=wfnUnit, binary=.true., arguments = labels ) 
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


  !!**********************************************************
  !! Save Some energies
  !!
  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalEnergy, arguments=["TOTALENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%cosmo3Energy, arguments=["COSMO3ENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalCouplingEnergy, arguments=["COUPLINGENERGY"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%electronicRepulsionEnergy, arguments=["COUPLING-E-"])

  call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MolecularSystem_getPointChargesEnergy(), arguments=["PUNTUALINTERACTIONENERGY"])

  !stop time
  call Stopwatch_stop(lowdin_stopwatch)

  if(CONTROL_instance%LAST_STEP) then
     write(*, *) ""
     write(*,"(A,F10.3,A4)") "** TOTAL CPU Time SCF : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time SCF : ", lowdin_stopwatch%elapsetWTime ," (s)"
     write(*, *) ""
  end if

  close(wfnUnit)

  if (CONTROL_instance%COSMO) then
     call WaveFunction_cosmoQuantumCharge()
  end if


  !!calculate HF/KS properties
  call system ("lowdin-CalcProp.x")
  
  
  if (CONTROL_instance%LOCALIZE_ORBITALS) then

     !! Asign orbitals to fragments and creates the subsystem fock matrix
     call OrbitalLocalizer_levelShiftSubSystemOrbitals()

  end if


  
end program SCF


