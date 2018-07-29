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
  implicit none

  real(8) :: auxValue
  real(8) :: deltaEnergy
  real(8) :: diisError
  integer :: i
  integer :: numberOfSpecies
  integer :: status
  integer :: wfnUnit, vecUnit
  integer :: speciesID
  integer :: numberOfIterations
  character(50) :: wfnFile, vecFile
  character(30) :: nameOfSpecie
  character(30) :: labels(2)
  character(100) :: iterationScheme(0:3)
  character :: convergenceType
  integer :: statusSystem

  !! Open file for wfn
  wfnUnit = 300
  wfnFile = "lowdin.wfn"

  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  !! Start the MultiSCF object
  call MultiSCF_constructor()

  !! Start the wavefunction object
  call WaveFunction_constructor( wfnUnit )

  numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

  if( numberOfSpecies > 1 ) then

     iterationScheme(0) = "NONELECRONIC FULLY CONVERGED BY ELECTRONIC ITERATION"
     iterationScheme(1) = "ELECTRONIC FULLY CONVERGED BY NONELECRONIC ITERATION"
     iterationScheme(2) = "SPECIES FULLY CONVERGED INDIVIDIALLY"
     iterationScheme(3) = "SIMULTANEOUS"

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
        statusSystem = system ("lowdin-DFT.x BUILD_MATRICES")
     end if
     
     do i = 1, numberOfSpecies
        nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
        
        call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie))

        !call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie))

        if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
           call WaveFunction_buildExchangeCorrelationMatrix( trim(nameOfSpecie))
        end if
        
        call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )

        if (CONTROL_instance%COSMO) then
           call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecie))
        end if
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
!        nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
!        call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), MultiSCF_instance%nproc )
!        call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )
!
!        if (CONTROL_instance%COSMO) then
!           call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecie))
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

              nameOfSpecie =  MolecularSystem_getNameOfSpecie(speciesID)                 
              numberOfIterations = List_size( WaveFunction_instance(speciesID)%energySCF )

              call List_begin( WaveFunction_instance(speciesID)%energySCF )
              call List_begin( WaveFunction_instance(speciesID)%diisError )
              call List_begin( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )

              print *,""
              print *,"Begin SCF calculation by: ",trim(nameOfSpecie)
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

  else
     
     call MultiSCF_iterate( CONTROL_instance%ITERATION_SCHEME )


  end if

  close(wfnUnit)

  ! Final integration grid goes here
  if ( (CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS") .and. &
       ( CONTROL_instance%FINAL_GRID_ANGULAR_POINTS*CONTROL_instance%FINAL_GRID_RADIAL_POINTS  .gt. &
       CONTROL_instance%GRID_ANGULAR_POINTS*CONTROL_instance%GRID_RADIAL_POINTS ) ) then
     statusSystem = system ("lowdin-DFT.x FINAL_GRID")
     do speciesID = 1, numberOfSpecies
        call WaveFunction_buildExchangeCorrelationMatrix( trim(MolecularSystem_getNameOfSpecie(speciesID)) )
     end do
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

  if ( .not. CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY .or. .not. CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then

     labels = ""
     !! Open file for vec
     vecUnit = 36
     vecFile = "lowdin-plain.vec"
     open(unit=vecUnit, file=trim(vecFile), form="formatted", status='unknown')

     if ( .not. CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then
        do speciesID = 1, numberOfSpecies

           labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
           labels(1) = "COEFFICIENTS"
           call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
                unit=vecUnit, binary=.false., arguments = labels)
        end do
     end if

     if ( .not. CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY ) then

        do speciesID = 1, numberOfSpecies

           labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
           labels(1) = "ORBITALS"
           call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, & 
                unit=vecUnit, binary=.false., arguments = labels )
        end do
     end if
     close (vecUnit)
  end if

!   vecUnit = 36
!   if ( CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY ) then
     
!      vecFile = trim(CONTROL_instance%INPUT_FILE)//"lowdin.vec"
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

!     if ( .not. CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY ) then

!       do speciesID = 1, numberOfSpecies

! <<<<<<< HEAD
!         labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
!         labels(1) = "ORBITALS"
!         call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, & 
!              unit=vecUnit, binary=.false., arguments = labels )
!       end do
!     end if
!     close (vecUnit)
! =======

  
! >>>>>>> 65f17bb51f2cb9cb348777354352854d08e936e2
! end if


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
     write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time SCF : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*, *) ""
  end if

  close(wfnUnit)

  if (CONTROL_instance%COSMO) then
     call WaveFunction_cosmoQuantumCharge()
  end if

end program SCF

