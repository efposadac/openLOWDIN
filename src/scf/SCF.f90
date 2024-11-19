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
  use String_
  use Exception_
  use omp_lib
  use OrbitalLocalizer_
  implicit none

  integer :: speciesID, otherSpeciesID
  integer :: wfnUnit
  character(50) :: job
  character(50) :: wfnFile
  character(30) :: labels(2)

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
     call MolecularSystem_showCartesianMatrix(molecularSystem_instance)
     call MolecularSystem_showDistanceMatrix()
  end if


  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !! Start the MultiSCF object
  allocate(WaveFunction_instance(MolecularSystem_instance%numberOfQuantumSpecies))
  call MultiSCF_constructor(MultiSCF_instance,WaveFunction_instance,CONTROL_instance%ITERATION_SCHEME,molecularSystem_instance)

  !! Calculate one-particle integrals  
  if ( CONTROL_instance%INTEGRAL_STORAGE == "DISK" ) &
       call system("lowdin-ints.x ONE_PARTICLE")

  !! Build hcore operators and use them to get guess (or read previous coefficients)
  call MultiSCF_buildHcore(MultiSCF_instance,WaveFunction_instance)

  call MultiSCF_getInitialGuess(MultiSCF_instance,WaveFunction_instance)

  !!**************************************************************************************************************************
  !! Calculate two-particle integrals (not building 2 particles and coupling matrix... those matrices updated at each SCF cycle)
  !!
  if ( CONTROL_instance%INTEGRAL_STORAGE == "DISK" ) then
     !! Save matrices to lowdin.wfn file required by ints program
     wfnUnit = 300
     wfnFile = "lowdin.wfn"
     open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
        labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
        labels(1) = "DENSITY"
        call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )
     end do
     close(wfnUnit)

     if( CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL ) then       
        call system(" lowdin-ints.x TWO_PARTICLE_G12")        
     else        
        call system(" lowdin-ints.x TWO_PARTICLE_R12")
     end if
  else if (CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then
     allocate(Libint2Instance(MolecularSystem_instance%numberOfQuantumSpecies))
     call DirectIntegralManager_constructor(Libint2Instance,MolecularSystem_instance)
     do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
        call DirectIntegralManager_getDirectIntraRepulsionIntegralsAll(&
             speciesID, &
             WaveFunction_instance(speciesID)%densityMatrix, & 
             WaveFunction_instance(speciesID)%fourCenterIntegrals(speciesID)%values, &
             MolecularSystem_instance,Libint2Instance(speciesID))
     end do

     do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies-1
        do otherSpeciesID=speciesID+1, MolecularSystem_instance%numberOfQuantumSpecies
           call DirectIntegralManager_getDirectInterRepulsionIntegralsAll(&
                speciesID, otherSpeciesID, &
                WaveFunction_instance(speciesID)%densityMatrix, & 
                WaveFunction_instance(speciesID)%fourCenterIntegrals(otherSpeciesID)%values, &
                MolecularSystem_instance,Libint2Instance(speciesID),Libint2Instance(otherSpeciesID))
        end do
     end do

  end if

  !!**********************************************************
  !! Obtain the energy and wavefunction coefficients iteratively
  !!
  call MultiSCF_solveHartreeFockRoothan(MultiSCF_instance,WaveFunction_instance,Libint2Instance)

  call MultiSCF_showResults(MultiSCF_instance,WaveFunction_instance)

  !!**********************************************************
  !! Save matrices to lowdin.wfn file
  !!
  call MultiSCF_saveWfn(MultiSCF_instance,WaveFunction_instance)

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

  if (CONTROL_instance%SUBSYSTEM_EMBEDDING) then
     !!calculate HF/KS properties for the full system
     call system ("lowdin-CalcProp.x")

     print *, ""
     print *, "-------------------------------------------"
     print *, "STARTING SUBSYSTEM EMBEDDED SCF CALCULATION"
     print *, "-------------------------------------------"
     print *, ""

     !! Asign orbitals to fragments and creates the subsystem fock matrix
     call OrbitalLocalizer_levelShiftSubSystemOrbitals()

  end if

end program SCF


