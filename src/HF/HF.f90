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
  implicit none

  type(Matrix) :: auxDensity
  character(50) :: job
  character(50) :: integralsFile
  character(50) :: wfnFile
  character(50) :: arguments(20)
  integer :: wfnUnit
  integer :: integralsUnit
  logical :: existFile
  integer :: numberOfSpecies
  integer :: numberOfContractions
  integer :: speciesID, i
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

  ! write(*,"(A)") trim(job)
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  if(CONTROL_instance%LAST_STEP) then
     write(*,"(A)")"----------------------------------------------------------------------"
     write(*,"(A)")"** PROGRAM: HF (Hartree Fock).      Author: S.A. Gonzalez, E. Posada  "
     write(*,"(A)")"----------------------------------------------------------------------"

     write(*,"(A)") "INFO: RUNNING IN "//trim(job)//" MODE."
     write(*,"(A)")" "
     write (6,"(T20,A30)") " TEST GEOMETRY: AMSTRONG"
     write (6,"(T18,A35)") "------------------------------------------"
     call MolecularSystem_showCartesianMatrix()
     call MolecularSystem_showDistanceMatrix()
  end if

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)


  call WaveFunction_constructor()

  integralsFile = "lowdin.opints"
  integralsUnit = 30
  wfnFile = "lowdin.wfn"
  wfnUnit = 20

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

  call system("lowdin-SCF.x ")

  if( .not.CONTROL_instance%OPTIMIZE .and. CONTROL_instance%GET_GRADIENTS ) then        

     call system("lowdin-ints.x GET_GRADIENTS")

  end if
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
     WaveFunction_instance(speciesID)%coefficientsofcombination = &
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
          output = WaveFunction_instance(speciesID)%energyofmolecularorbital )     

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

     !!Obtain energy components for species
     call WaveFunction_obtainEnergyComponents(speciesID)

  end do


  !! Obtain energy compotents for whole system
  totalKineticEnergy = sum( WaveFunction_instance(:)%kineticEnergy)             
  totalRepulsionEnergy = sum( WaveFunction_instance(:)%repulsionEnergy ) + electronicRepulsionEnergy                          
  totalQuantumPuntualInteractionEnergy = sum ( WaveFunction_instance(:)%puntualInteractionEnergy )
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

  close(wfnUnit)

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

        !!call Matrix_show( WaveFunction_instance(speciesID)%coefficientsofcombination, &
        !!     rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
        !!     columnkeys = string_convertvectorofrealstostring( WaveFunction_instance(speciesID)%energyofmolecularorbital ),&
        !!     flags=WITH_BOTH_KEYS)

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
     write (6,"(T10,A28,F20.10)") "TOTAL KINETIC ENERGY      = ", sum(WaveFunction_instance(:)%kineticEnergy)
     write (6,"(T10,A28,F20.10)") "TOTAL POTENTIAL ENERGY    = ", potentialEnergy
     write (6,"(T10,A50)") "________________"
     write (6,"(T10,A28,F20.12)") "TOTAL ENERGY = ", totalEnergy             
     write(*,*) ""
     write (6,"(T10,A28,F20.10)") "VIRIAL RATIO (V/T) = ", - ( potentialEnergy / totalKineticEnergy)
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

        write (6,"(T10,A26,A2,F20.10)") "e-ALPHA" // "/" // &
             "e-BETA" // " Repulsion energy  ","= ", electronicRepulsionEnergy
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
           write(6,"(T2,I4,F20.10)") i,WaveFunction_instance(speciesID)%energyofmolecularorbital%values(i)
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

  if(CONTROL_instance%LAST_STEP) then
     write(*, *) ""
     write(*,"(A,F10.3,A4)") "** TOTAL CPU Time HF : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time HF : ", lowdin_stopwatch%elapsetWTime ," (s)"
     write(*, *) ""
  end if

  close(30)

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
