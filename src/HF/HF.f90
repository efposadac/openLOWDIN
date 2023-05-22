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
!! @brief HF and APMO-HF program.
!!        This module allows to make calculations in the APMO-HF framework
!! @author  E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-05-09
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-09-17 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos (RHF.f90)
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adpata el modulo para incluirlo en LOWDIN (RHF.90)
!!   - <tt> 2013-05-09 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrite the module as a program, handles closed and open shell systems
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program HF
  use CONTROL_
  use WaveFunction_
  use MolecularSystem_
  use DensityMatrixSCFGuess_
  use String_
  use Exception_
  use omp_lib
  use OrbitalLocalizer_
  implicit none




  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)


  call WaveFunction_constructor()

  integralsFile = "lowdin.opints"
  integralsUnit = 30
  wfnFile = "lowdin.wfn"
  wfnUnit = 20


  call system("lowdin-SCF.x ")

  !! Check the one-particle integrals file  
  existFile = .false.     
  inquire(file=trim(integralsFile), exist=existFile)
  
  if( existFile ) then
     
     open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")
     
     read(integralsUnit) numberOfSpecies
     
     if(MolecularSystem_instance%numberOfQuantumSpecies /= numberOfSpecies ) then
        
        call MolecularSystem_exception( ERROR, "Bad "//trim(integralsFile)//" file!", "In HF.f90 at main program")

     end if
     
     close(integralsUnit)
     
  else
     
     call MolecularSystem_exception(ERROR,"lowdin.opints file not found!", "In HF.f90 at main program")
     
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
     
     !! Hcore Matrix
     call WaveFunction_HCoreMatrix(trim(integralsFile), speciesID)
     
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
     
  end do
  
  close(wfnUnit)
  
  !!**************************************************************************************************************************
  !! Calculate two-particle integrals (not building 2 particles and coupling matrix... those matrices are done by SCF program)
  !!
  select case (trim(String_getUppercase(trim(CONTROL_instance%INTEGRAL_SCHEME))))
     
  case("RYS")
     write(*,  "(A)")  " RYS QUADRATURE SCHEME                 " 
     write(*, "(A)")   " LOWDIN-RYS Implementation V. 1.0   Guerrero R. D. ; Posada E. F. 2013 "
     write(*, "(A)")   " ----------------------------------------------------------------------"
     
  case("LIBINT")
     write(*,  "(A)")  " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
     write(*, "(A)")   " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2011   "
     write(*, "(A)")   " ----------------------------------------------------------------------"
     
  case default
     write(*,  "(A)")  " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
     write(*, "(A)")   " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2011   "
     write(*, "(A)")   " ----------------------------------------------------------------------"
     
  end select

  if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL ) then        
     
     call system(" lowdin-ints.x TWO_PARTICLE_F12")
     
  else        
     
     nproc = CONTROL_instance%NUMBER_OF_CORES
     
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
        
        !$OMP PARALLEL firstprivate(arguments)
        write(arguments(1),*) nproc
        !$OMP DO 
        do i = nproc, 1, -1
           
           write(arguments(i+1),*) i
           call system(" lowdin-ints.x TWO_PARTICLE_R12_INTRA "//trim(arguments(1))//" "//trim(arguments(i+1))//" "//trim(MolecularSystem_getNameOfSpecie(speciesID)))          
           
        end do
        !$OMP END DO
        !$OMP END PARALLEL            

     end do
     
     call system(" lowdin-ints.x TWO_PARTICLE_R12_INTER")
        
  end if
  
  write(*,*) "DONE!"

  !!
  !!***************************************************************************************************************

  !! Begin SCF calculation...  

  write(arguments(1),*) nproc
  call system("lowdin-SCF.x "//trim(arguments(1)))

  !! Open file for wavefunction
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  !! Load results...
  call Vector_getFromFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])
  call Vector_getFromFile(unit=wfnUnit, binary=.true., value=cosmo3Energy, arguments=["COSMO3ENERGY"])
  call Vector_getFromFile(unit=wfnUnit, binary=.true., value=totalCouplingEnergy, arguments=["COUPLINGENERGY"])
  call Vector_getFromFile(unit=wfnUnit, binary=.true., value=electronicRepulsionEnergy, arguments=["COUPLING-E-"])

  !!Load matrices
  do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
     arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

     arguments(1) = "DENSITY"
     WaveFunction_instance(speciesID)%densityMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "COEFFICIENTS"
     WaveFunction_instance(speciesID)%waveFunctionCoefficients = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "COUPLING"
     WaveFunction_instance(speciesID)%couplingMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "EXCHANGE-CORRELATION"
     WaveFunction_instance(speciesID)%exchangeCorrelationMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "EXCHANGE-CORRELATION-ENERGY"
     call Vector_getFromFile( value=WaveFunction_instance(speciesID)%exchangeCorrelationEnergy, unit=wfnUnit, &
          binary=.true., arguments=arguments(1:2))
          
     arguments(1) = "TWOPARTICLES"
     WaveFunction_instance(speciesID)%twoParticlesMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "ORBITALS"
     call Vector_getFromFile( elementsNum = numberOfContractions, &
          unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = WaveFunction_instance(speciesID)%molecularOrbitalsEnergy )     

     if(CONTROL_instance%COSMO)then

        arguments(1) = "COSMO2"
        WaveFunction_instance(speciesID)%cosmo2 = &
             Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "COSMOCOUPLING"
        WaveFunction_instance(speciesID)%cosmoCoupling = &
             Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     end if


  end do



  close(wfnUnit)




  
  !! not necessary for now...
  select case(trim(job))         

  case("RHF")

  case("UHF")

  case("RKS")

  case("UKS")

  case default

     write(*,*) "USAGE: lowdin-HF.x job "
     write(*,*) "Where job can be: "
     write(*,*) "  RHF"
     write(*,*) "  UHF"
     write(*,*) "  RKS"
     write(*,*) "  UKS"
     stop "ERROR"

  end select

end program HF
