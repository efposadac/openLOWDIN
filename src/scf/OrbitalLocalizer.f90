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
!! @brief This module handles orbital localization routines to be used during or after SCF calculations
!! @author F. Moncada 2020

module OrbitalLocalizer_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use MolecularSystem_
  use ParticleManager_
  use BasisSet_
  use WaveFunction_
  use Convergence_
  use omp_lib

  implicit none


  type, public :: OrbitalLocalizer
     real(8) :: levelShiftingValue !control parameter
     real(8) :: orbitalThreshold !control parameter
     real(8) :: basisThreshold !control parameter
     type(Matrix) :: projectionMatrix
     type(Matrix) :: hcoreMatrixA
     type(Matrix) :: fockMatrixA
     type(Matrix) :: waveFunctionCoefficientsA
     type(Vector) :: molecularOrbitalsEnergyA
     integer,allocatable :: subsystemList(:)
     integer,allocatable :: orbitalSubsystem(:)
     integer,allocatable :: basisSubsystem(:)
     integer :: occupiedOrbitalsA
     integer :: virtualOrbitalsA
     integer :: occupiedOrbitalsB
     integer :: virtualOrbitalsB
     integer :: removedOrbitalsA

  end type OrbitalLocalizer

  type(OrbitalLocalizer), public, allocatable :: OrbitalLocalizer_instance(:)

contains
 
  subroutine OrbitalLocalizer_constructor()
    integer :: speciesID, otherSpeciesID, numberOfSpecies
    integer :: numberOfContractions, otherNumberOfContractions, occupationNumber
    integer :: mu, i, j, k

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies
    allocate(OrbitalLocalizer_instance(numberOfSpecies))

    do speciesID=1, numberOfSpecies

       OrbitalLocalizer_instance(speciesID)%levelShiftingValue = CONTROL_INSTANCE%SUBSYSTEM_LEVEL_SHIFTING
       OrbitalLocalizer_instance(speciesID)%orbitalThreshold = CONTROL_INSTANCE%SUBSYSTEM_ORBITAL_THRESHOLD
       OrbitalLocalizer_instance(speciesID)%basisThreshold = CONTROL_INSTANCE%SUBSYSTEM_BASIS_THRESHOLD

       !Gets atomic subsystem classification
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       occupationNumber = MolecularSystem_getOcupationNumber( speciesID )
       
       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%hcoreMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%fockMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Vector_constructor(OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA, numberOfContractions, 0.0_8)
       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%projectionMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

       
       allocate(OrbitalLocalizer_instance(speciesID)%subsystemList(numberOfContractions))
       allocate(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(numberOfContractions))
       allocate(OrbitalLocalizer_instance(speciesID)%basisSubsystem(numberOfContractions))

       mu=0
       do i = 1, size(MolecularSystem_instance%species(speciesID)%particles)
          do j = 1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             do k = 1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
                mu=mu+1
                if(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%subsystem .eq. 1 ) then
                   OrbitalLocalizer_instance(speciesID)%subsystemList(mu)=1
                else
                   OrbitalLocalizer_instance(speciesID)%subsystemList(mu)=2
                end if
             end do
          end do
       end do       
    end do

  end subroutine OrbitalLocalizer_constructor

  subroutine OrbitalLocalizer_erkaleLocal(speciesID,densityMatrix,fockMatrix,orbitalCoefficients,orbitalEnergies)
    integer :: speciesID
    type(Matrix) :: densityMatrix
    type(Matrix) :: fockMatrix
    type(Matrix) :: orbitalCoefficients
    type(Vector) :: orbitalEnergies

    character(30) :: nameOfSpecies, symbolOfSpecies
    integer :: statusSystem

    integer :: numberOfContractions
    integer :: mu,nu, i

    
    !! Convert lowdin fchk files to erkale chk files
    nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
    symbolOfSpecies=MolecularSystem_getSymbolOfSpecies(speciesID)
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    open(unit=30, file="erkale.read", status="replace", form="formatted")

    write(30,*) "LoadFChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".fchk"
    write(30,*) "SaveChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".chk"
    write(30,*) "Reorthonormalize true"
    
    close(30)

    call system("erkale_fchkpt erkale.read")

    !! Localize orbitals
    open(unit=30, file="erkale.local", status="replace", form="formatted")
    
    write(30,*) "LoadChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".chk"
    write(30,*) "SaveChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".local.chk"
    write(30,*) "Method ", trim(CONTROL_instance%ERKALE_LOCALIZATION_METHOD)
    write(30,*) "Virtual false"
    write(30,*) "Maxiter 5000"
    write(30,*) "StartingPoint CAN"
    write(30,*) "FThreshold 1E-8"
    write(30,*) "GThreshold 1E-8"
    

    close(30)

    call system("erkale_loc_omp erkale.local")
    
    !!Convert erkale chk files to lowdin fchk files
    open(unit=30, file="erkale.write", status="replace", form="formatted")
    
    write(30,*) "LoadChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".local.chk"
    write(30,*) "SaveFChk ", trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".local.fchk"
    
    close(30)

    call system("erkale_fchkpt erkale.write")
    
    !! Read orbital coefficients from fchk files
    call MolecularSystem_readFchk(trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".local.fchk",  orbitalCoefficients, densityMatrix, nameOfSpecies )

    orbitalEnergies%values=0.0
    !! Molecular orbital fock operator expected value
    do i=1, numberOfContractions
       do mu=1, numberOfContractions
          do nu=1, numberOfContractions
             orbitalEnergies%values(i)=&
                  orbitalEnergies%values(i)+&
                  fockMatrix%values(mu,nu)*&
                  orbitalCoefficients%values(mu,i)*&
                  orbitalCoefficients%values(nu,i)
          end do
       end do
    end do
    
  end subroutine OrbitalLocalizer_erkaleLocal

  
  subroutine OrbitalLocalizer_levelShiftSubsystemOrbitals()
    integer :: speciesID, otherSpeciesID
    type(Matrix) :: coefficients
    type(Matrix) :: fockMatrix
    type(Matrix) :: overlapMatrix
    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: newDensityMatrix
    real(8) :: sumSubsystem(2), sumAB, maxSumAB, normCheck
    real(8) :: totalEnergy, newTotalEnergy, deltaEnergy
    integer :: shellSize
    
    type(BasisSet) :: reducedBasisA
    
    type(Matrix),allocatable :: auxMatrix(:)
    type(Matrix),allocatable :: densityMatrixA(:)
    type(Matrix),allocatable :: densityMatrixB(:)
    type(Matrix),allocatable :: densityMatrixAB(:)
    type(Matrix),allocatable :: hartreeMatrixA(:,:)
    type(Matrix),allocatable :: hartreeMatrixB(:,:)
    type(Matrix),allocatable :: hartreeMatrixAB(:,:)
    type(Matrix),allocatable :: couplingMatrixA(:)
    type(Matrix),allocatable :: couplingMatrixB(:)
    type(Matrix),allocatable :: couplingMatrixAB(:)
    type(Matrix),allocatable :: exchangeCorrelationMatrixA(:)
    type(Matrix),allocatable :: exchangeCorrelationMatrixB(:)
    type(Matrix),allocatable :: exchangeCorrelationMatrixAB(:)
    type(Matrix),allocatable :: exchangeHFMatrixA(:)
    type(Matrix),allocatable :: exchangeHFMatrixB(:)
    type(Matrix),allocatable :: exchangeHFMatrixAB(:)
    type(Matrix),allocatable :: populationMatrixA(:)
    type(Matrix),allocatable :: populationMatrixB(:)
    type(Vector),allocatable :: atomPopulationA(:)
    type(Vector),allocatable :: atomPopulationB(:)
    type(Vector),allocatable :: shellPopulationA(:)
    type(Vector),allocatable :: shellPopulationB(:)
    real(8),allocatable :: excCorrEnergyA(:,:)
    real(8),allocatable :: excCorrEnergyB(:,:)
    real(8),allocatable :: excCorrEnergyAB(:,:)
    real(8),allocatable :: particlesInGridA(:)
    real(8),allocatable :: particlesInGridB(:)
    real(8),allocatable :: particlesInGridAB(:)
    real(8),allocatable :: densStd(:)
    real(8),allocatable :: auxEnergy(:)
    character(9), allocatable :: shellCode(:)
    character(9), allocatable :: holdPointChargeSymbol(:)
    real(8), allocatable :: holdPointChargeMass(:)
    
    real(8) :: totalDensStd

    logical :: includeAtom
    logical :: SUBSYSTEM_SCF_CONTINUE
    
    character(30) :: nameOfSpecies, nickname

    integer :: numberOfSpecies,numberOfContractions, numberOfCenters, numberOfShells, occupationNumber
    integer :: subsystemMu, subsystemNu, orbitalsInA, orbitalsInB
    integer :: mu,nu,alpha,beta, i, ii, j, jj, k, l, iter
    integer :: oldMu, oldNu

    character(50) :: labels(2)
    character(100) :: densFileA,densFileB,densFileAB, wfnFile, integralsFile, vecFile
    integer :: densUnitA,densUnitB,densUnitAB, wfnUnit, vecUnit, integralsUnit

    real(8) :: levelShiftingFactor
    real(8) :: totalKineticEnergyA,totalKineticEnergyB, puntualInteractionEnergy
    real(8) :: totalQuantumPuntualInteractionEnergyA, totalQuantumPuntualInteractionEnergyB
    real(8) :: totalExternalPotentialEnergyA, totalExternalPotentialEnergyB
    real(8) :: totalHartreeEnergyA, totalHartreeEnergyB, totalHartreeEnergyAB
    real(8) :: totalExchangeHFEnergyA, totalExchangeHFEnergyB, totalExchangeHFEnergyAB 
    real(8) :: totalExchangeCorrelationEnergyA, totalExchangeCorrelationEnergyB, totalExchangeCorrelationEnergyAB
    real(8) :: totalEmbeddingPotentialEnergyA, totalProjectionCorrectionA=0.0
    real(8) :: totalPotentialEnergy, totalKineticEnergy
    real(8) :: kAB, pAB
    
    densUnitA = 77
    densUnitB = 78
    densUnitAB = 79
    
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

    allocate(hartreeMatrixA(numberOfSpecies,numberOfSpecies))
    allocate(hartreeMatrixB(numberOfSpecies,numberOfSpecies))
    allocate(hartreeMatrixAB(numberOfSpecies,numberOfSpecies))

    allocate(auxMatrix(numberOfSpecies))
    allocate(densityMatrixA(numberOfSpecies))
    allocate(densityMatrixB(numberOfSpecies))
    allocate(densityMatrixAB(numberOfSpecies))
    allocate(couplingMatrixA(numberOfSpecies))
    allocate(couplingMatrixB(numberOfSpecies))
    allocate(couplingMatrixAB(numberOfSpecies))
    allocate(exchangeCorrelationMatrixA(numberOfSpecies))
    allocate(exchangeCorrelationMatrixB(numberOfSpecies))
    allocate(exchangeCorrelationMatrixAB(numberOfSpecies))
    allocate(exchangeHFMatrixA(numberOfSpecies))
    allocate(exchangeHFMatrixB(numberOfSpecies))
    allocate(exchangeHFMatrixAB(numberOfSpecies))
    allocate(populationMatrixA(numberOfSpecies))
    allocate(populationMatrixB(numberOfSpecies))
    allocate(atomPopulationA(numberOfSpecies))
    allocate(atomPopulationB(numberOfSpecies))
    allocate(shellPopulationA(numberOfSpecies))
    allocate(shellPopulationB(numberOfSpecies))
    allocate(excCorrEnergyA(numberOfSpecies,numberOfSpecies))
    allocate(excCorrEnergyB(numberOfSpecies,numberOfSpecies))
    allocate(excCorrEnergyAB(numberOfSpecies,numberOfSpecies))
    allocate(particlesInGridA(numberOfSpecies))
    allocate(particlesInGridB(numberOfSpecies))
    allocate(particlesInGridAB(numberOfSpecies))
    allocate(densStd(numberOfSpecies))
    allocate(auxEnergy(numberOfSpecies))

    densFileA ="lowdin.densmatrixA"
    densFileB ="lowdin.densmatrixB"
    densFileAB ="lowdin.densmatrixAB"
    
    do speciesID=1, numberOfSpecies

       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       numberOfCenters = size(MolecularSystem_instance%species(speciesID)%particles)
       occupationNumber = MolecularSystem_getOcupationNumber( speciesID )
       overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
       coefficients=WaveFunction_instance( speciesID )%waveFunctionCoefficients
       fockMatrix=WaveFunction_instance( speciesID )%fockMatrix
       call Matrix_constructor (densityMatrixA(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (densityMatrixB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (densityMatrixAB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (couplingMatrixA(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (couplingMatrixB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (couplingMatrixAB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeCorrelationMatrixA(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeCorrelationMatrixB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeCorrelationMatrixAB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeHFMatrixA(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeHFMatrixB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeHFMatrixAB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (populationMatrixA(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (populationMatrixB(speciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Vector_constructor (atomPopulationA(speciesID), numberOfCenters, 0.0_8)
       call Vector_constructor (atomPopulationB(speciesID), numberOfCenters, 0.0_8)
       call Vector_constructor (shellPopulationA(speciesID), MolecularSystem_getNumberOfContractions(speciesID), 0.0_8)
       call Vector_constructor (shellPopulationB(speciesID), MolecularSystem_getNumberOfContractions(speciesID), 0.0_8)

       do otherSpeciesID=1, numberOfSpecies
          call Matrix_constructor(hartreeMatrixA(speciesID,otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
          call Matrix_constructor(hartreeMatrixB(speciesID,otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
          call Matrix_constructor(hartreeMatrixAB(speciesID,otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       end do

       excCorrEnergyA(speciesID,:)=0.0
       excCorrEnergyB(speciesID,:)=0.0
       excCorrEnergyAB(speciesID,:)=0.0
       
       !Assigns occupied molecular orbital to each subsystem
       OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA=0
       OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB=0
       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA=0
       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=0

       print *, ""
       print *, "Occupied orbitals subsystem separation for ", nameOfSpecies
       print *, ""
       write(*,"(A15,A15,A15)") "Orbital", "Subsystem", "A Population"
       
       ! print *, "orbital, assigned to, sumSubsystem(1), sumSubsystem(1), sumAB, normCheck"
       !!Occupied
       do k = 1 , occupationNumber
          normCheck=0.0
          sumSubsystem(:)=0.0
          sumAB=0.0
          do mu = 1 , numberOfContractions
             do nu = 1 , numberOfContractions
                if( OrbitalLocalizer_instance(speciesID)%subsystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subsystemList(nu) ) then
                   sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))=&
                        sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))+&
                        overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
                ! else
                !    sumAB=sumAB+&
                !         overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
                end if
                ! normCheck=normCheck+coefficients%values(mu,k)*&
                !      coefficients%values(nu,k)*&
                !      WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
             end do
          end do
          if(sumSubsystem(1) .gt. OrbitalLocalizer_instance(speciesID)%orbitalThreshold) then
             OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+1
             OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 1
             write(*,"(I15,A15,F12.6)") k, "A    ", sumSubsystem(1)
          else
             OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB+1
             OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 2
             write(*,"(I15,A15,F12.6)") k, "B    ", sumSubsystem(1)
          end if
          ! print *, k, OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k), sumSubsystem(1), sumSubsystem(2), sumAB, normCheck
       end do
       print *, ""

       !!Virtual
       ! do k = occupationNumber+1 , numberOfContractions
       !    normCheck=0.0
       !    sumSubsystem(:)=0.0
       !    sumAB=0.0
       !    do mu = 1 , numberOfContractions
       !       do nu = 1 , numberOfContractions
       !          if( OrbitalLocalizer_instance(speciesID)%subsystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subsystemList(nu) ) then
       !             sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))=&
       !                  sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))+&
       !                  overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
       !          else
       !             sumAB=sumAB+&
       !                  overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
       !          end if
       !          normCheck=normCheck+coefficients%values(mu,k)*&
       !               coefficients%values(nu,k)*&
       !               WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
       !       end do
       !    end do
       !    if(sumAB**2 .gt. (sumSubsystem(1)**2+sumSubsystem(2)**2) ) then
       !       OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals=OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 0
       !    else if(sumSubsystem(1)**2 .gt. 0.1) then
       !       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 1
       !    else
       !       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 2
       !    end if
       !    print *, k, OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k), sumSubsystem(1), sumSubsystem(2), sumAB, normCheck
       ! end do
       ! print *, "virtualOrbitalsA, virtualOrbitalsB, virtualDelocalizedOrbitals", &
       !      OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA, OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB, OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals


       !Decomposes atomic density matrix
       do k = 1 , occupationNumber
          if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k).eq.1) then
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   densityMatrixA(speciesID)%values(mu,nu)=&
                        densityMatrixA(speciesID)%values(mu,nu)&
                        +MolecularSystem_getEta( speciesID )*coefficients%values(mu,k)*coefficients%values(nu,k)
                end do
             end do
          else
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   densityMatrixB(speciesID)%values(mu,nu)=&
                        densityMatrixB(speciesID)%values(mu,nu)&
                        +MolecularSystem_getEta( speciesID )*coefficients%values(mu,k)*coefficients%values(nu,k)
                end do
             end do
          end if
       end do

       !Reorders coefficients Matrix, first occupied orbitals in A
       i=0
       j=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA
       do k = 1 , numberOfContractions
          if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k).eq.1) then
             i=i+1
             OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,i)=WaveFunction_instance( speciesID )%waveFunctionCoefficients%values(:,k)
          else
             j=j+1
             OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,j)=WaveFunction_instance( speciesID )%waveFunctionCoefficients%values(:,k)
          end if
       end do
          
       
       !Forces equal density matrices for E-ALPHA and E-BETA in open shell calculations
       if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
            (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") .and. &
            trim(nameOfSpecies) .eq. "E-BETA"  ) then
          otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-ALPHA")  )
          if(MolecularSystem_getNumberOfParticles(speciesID) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID) ) then
             densityMatrixB(speciesID)%values=densityMatrixB(otherSpeciesID)%values
          end if
       end if

       
       densityMatrixAB(speciesID)%values=densityMatrixA(speciesID)%values+densityMatrixB(speciesID)%values
       
       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "densityMatrixA", nameOfSpecies
          Call Matrix_show(densityMatrixA(speciesID))
          print *, "densityMatrixB", nameOfSpecies
          Call Matrix_show(densityMatrixB(speciesID))
       end if

       !Build Proyection Matrix to force orbital orthogonality between A and B occupied orbitals
       do alpha = 1 , numberOfContractions
          do beta = 1 , numberOfContractions
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(alpha,beta)=&
                        OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(alpha,beta)&
                        +overlapMatrix%values(alpha,mu)*overlapMatrix%values(nu,beta)&
                        *densityMatrixB(speciesID)%values(mu,nu)                
                end do
             end do
          end do
       end do
       
       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "projectionMatrix", nameOfSpecies
          Call Matrix_show(OrbitalLocalizer_instance(speciesID)%projectionMatrix)
       end if
       
       
    end do
    
    !!Reduce basis set loops
    do speciesID=1, numberOfSpecies

       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       numberOfCenters = size(MolecularSystem_instance%species(speciesID)%particles)
       overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
       
       !Calculates atomic Mulliken population on each fragment
       if( trim(nameOfSpecies) .eq. "E-ALPHA"  ) then
          otherSpeciesID=speciesID+1
          
          populationMatrixA(speciesID)%values= matmul(densityMatrixA(speciesID)%values, overlapMatrix%values )+&
               matmul(densityMatrixA(otherSpeciesID)%values, overlapMatrix%values )
          populationMatrixB(speciesID)%values= matmul(densityMatrixB(speciesID)%values, overlapMatrix%values )+&
               matmul(densityMatrixA(otherSpeciesID)%values, overlapMatrix%values )
          nickname="E-"
          
       else if( trim(nameOfSpecies) .eq. "E-BETA" ) then
          otherSpeciesID=speciesID-1
         
          shellPopulationA(speciesID)%values=shellPopulationA(otherSpeciesID)%values
          shellPopulationB(speciesID)%values=shellPopulationA(otherSpeciesID)%values
          
          populationMatrixA(speciesID)%values=populationMatrixA(otherSpeciesID)%values
          populationMatrixB(speciesID)%values=populationMatrixA(otherSpeciesID)%values

          OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=OrbitalLocalizer_instance(otherSpeciesID)%virtualOrbitalsB
          OrbitalLocalizer_instance(speciesID)%basisSubsystem=OrbitalLocalizer_instance(otherSpeciesID)%basisSubsystem
          
          cycle
       else
          populationMatrixA(speciesID)%values= matmul(densityMatrixA(speciesID)%values, overlapMatrix%values )
          populationMatrixB(speciesID)%values= matmul(densityMatrixB(speciesID)%values, overlapMatrix%values )
          nickname=nameOfSpecies
       end if
       
       print *, ""
       print *, "Subsystem Max. Shell Mulliken Populations: ", nickname
       print *, ""
       write(*,"(A10,A10,A10,A10)") "Atom ","Shell A.M." , "A   ", "B   "
       mu=0
       nu=0
       OrbitalLocalizer_instance(speciesID)%basisSubsystem(:)=1
       do i=1, numberOfCenters
          atomPopulationA(speciesID)%values(i)=0
          atomPopulationB(speciesID)%values(i)=0
          do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             nu=nu+1
             shellPopulationA(speciesID)%values(nu)=0
             shellPopulationB(speciesID)%values(nu)=0
             shellSize=MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
             do k=1, shellSize
                mu=mu+1
                if(abs(populationMatrixA(speciesID)%values(mu,mu)) .gt. shellPopulationA(speciesID)%values(nu) ) &
                     shellPopulationA(speciesID)%values(nu)=abs(populationMatrixA(speciesID)%values(mu,mu))     
                   
                if(abs(populationMatrixB(speciesID)%values(mu,mu)) .gt. shellPopulationB(speciesID)%values(nu) ) &
                     shellPopulationB(speciesID)%values(nu)=abs(populationMatrixB(speciesID)%values(mu,mu))     

                atomPopulationA(speciesID)%values(i)=atomPopulationA(speciesID)%values(i)+populationMatrixA(speciesID)%values(mu,mu)
                atomPopulationB(speciesID)%values(i)=atomPopulationB(speciesID)%values(i)+populationMatrixB(speciesID)%values(mu,mu)
             end do
             if (abs(shellPopulationA(speciesID)%values(nu)) .lt. OrbitalLocalizer_instance(speciesID)%basisThreshold .and. &
                  OrbitalLocalizer_instance(speciesID)%subsystemList(mu) .ne. 1 ) then
                write(*,"(A1,I3,A6,I10,F10.6,F10.6)") "*",i,trim(ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner)),&
                     MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%angularMoment,&
                     shellPopulationA(speciesID)%values(nu), shellPopulationB(speciesID)%values(nu)
                OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu-shellSize+1:mu)=2
                OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=&
                     OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB+shellSize
             else
                write(*,"(I4,A6,I10,F10.6,F10.6)") i,trim(ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner)),&
                     MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%angularMoment,&
                     shellPopulationA(speciesID)%values(nu), shellPopulationB(speciesID)%values(nu)
                OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu-shellSize+1:mu)=1
             end if
          end do
          ! write(*,"(T20,A20)") "____________________"
          ! write(*,"(A10,A10,F10.6,F10.6)") trim(ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner)),"total:",&
          !      atomPopulationA(speciesID)%values(i), atomPopulationB(speciesID)%values(i)
          ! print *, ""
       end do
       ! write(*,"(T20,A20)") "--------------------"
       ! write(*,"(A20,F10.6,F10.6)") "Subsystem Total", sum(atomPopulationA(speciesID)%values(:)), sum(atomPopulationB(speciesID)%values(:))

       print *, "* will be projected out from A basis set"
    end do
       
    do speciesID=1, numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

       !Adds diagonal proyection elements to orbitals with small contributions to A orbitals - small mulliken population
       !These will correspond to virtual orbitals with high energy that can be removed from post-SCF calculations
       ! print *, "basis", "subsystem"
       do mu=1,numberOfContractions
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
          ! print *, mu, OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu)
          if (OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu) .ne. 1 ) then
             OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)=&
                  OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)+1.0
          end if
       end do

       !Adds diagonal proyection elements to atoms with small contributions to A orbitals - small mulliken population
       !These will correspond to virtual orbitals with high energy that can be removed from post-SCF calculations
       ! print *, "Basis Subsystem"
       ! OrbitalLocalizer_instance(speciesID)%basisSubsystem(:)=1
       ! if(OrbitalLocalizer_instance(speciesID)%reduceSubsystemBasis) then
       !    mu=0
       !    OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=0
       !    do i=1, numberOfCenters
       !       do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
       !          do k=1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
       !             mu=mu+1
       !             if (abs(atomPopulationA(speciesID)%values(i)) .lt. OrbitalLocalizer_instance(speciesID)%orbitalThreshold ) then
       !                OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=&
       !                     OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB+1
       !                OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)=&
       !                     OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)+1.0
       !                OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu)=2
       !             else
       !                OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu)=1
       !             end if
       !             ! print *, mu, OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu)
       !          end do
       !       end do
       !    end do
       ! end if
       
       OrbitalLocalizer_instance(speciesID)%removedOrbitalsA=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB&
            +OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB

       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA=numberOfContractions&
            -OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB

       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "projectionMatrix-2", nameOfSpecies
          Call Matrix_show(OrbitalLocalizer_instance(speciesID)%projectionMatrix)
       end if
       

       !!Manual fixes for some problematic cases
       if(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA .eq.0 ) OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA = 0
       if(OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA .lt. OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA ) then
          OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB = OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB - &
               (OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA - OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA)
          OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA = OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA
       end if
       if(OrbitalLocalizer_instance(speciesID)%removedOrbitalsA .gt. numberOfContractions ) OrbitalLocalizer_instance(speciesID)%removedOrbitalsA=numberOfContractions

       if(OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA .eq. 0 ) then
          OrbitalLocalizer_instance(speciesID)%removedOrbitalsA = WaveFunction_instance(speciesID)%removedOrbitals
          OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA = OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA - OrbitalLocalizer_instance(speciesID)%removedOrbitalsA
       end if

    end do
    
  !!Modify and write molecular system for subsystem A
    print *, ""
    print *, "**************************************************************************"
    print *, "*************Molecular subsystem A****************************************"
    print *, "**************************************************************************"

    !For now, we only remove the particles from system B
    !We will rebuild the whole molecular system later
    do speciesID=1, numberOfSpecies
       MolecularSystem_instance%species(speciesID)%ocupationNumber=&
            MolecularSystem_instance%species(speciesID)%ocupationNumber&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB
       MolecularSystem_instance%species(speciesID)%internalSize=&
            MolecularSystem_instance%species(speciesID)%internalSize&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB&
            *MolecularSystem_getEta(speciesID)
    end do

    print *, ""
    print *, "Orbital Subsystem Distribution"
    print *, ""
    write(*,"(A15,A15,A15,A15,A15)") "Species", "Occupied A","Virtual A","Occupied B", "Virtual B"
    write(*,"(A75)") "---------------------------------------------------------------------------"
    do speciesID=1, numberOfSpecies
       write(*,"(A15,I7,F7.1,A,I7,F7.1,A,I7,F7.1,A,I7,F7.1,A)") trim(MolecularSystem_getNameOfSpecies(speciesID)), &
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA,&
            100.0_8*OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA/(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB),&
            "%",&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA,&
            100.0_8*OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA/(OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA+OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB),&
            "%",&
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB,&
            100.0_8*OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB/(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB),&
            "%",&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB,&
            100.0_8*OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB/(OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA+OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB),&
            "%"

    end do
   
    !!Save density matrix B and run DFT for the subsystem B   
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileB),densityMatrixB(:))
       call system("lowdin-DFT.x SCF_DFT "//trim(densFileB))
    end if

    !!!Build subsystem B matrices - these do not change in the second SCF cycle
    !Two particles and coupling matrices
    do speciesID=1, numberOfSpecies     
       
       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       !Only Coulomb (factor=0.0)
       call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
         densityMatrixIN=densityMatrixB(speciesID),&
         factorIN=0.0_8,&
         twoParticlesMatrixOUT=hartreeMatrixB(speciesID,speciesID) )

       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "CoulombMatrixB", nameOfSpecies
          Call Matrix_show(hartreeMatrixB(speciesID,speciesID))
       end if

       !Exchange-HF matrix with the HF fraction used in the global SCF calculation       
       call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
         densityMatrixIN=densityMatrixB(speciesID),&
         twoParticlesMatrixOUT=exchangeHFMatrixB(speciesID))

       exchangeHFMatrixB(speciesID)%values=&
            exchangeHFMatrixB(speciesID)%values-hartreeMatrixB(speciesID,speciesID)%values

       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "exchangeHFMatrixB", nameOfSpecies
          Call Matrix_show(exchangeHFMatrixB(speciesID))
       end if
       
       do otherSpeciesID=1, numberOfSpecies
          call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       end do

       call WaveFunction_buildCouplingMatrix(WaveFunction_instance,speciesID,&
         densityMatricesIN=densityMatrixB(1:numberOfSpecies),&
         couplingMatrixOUT=couplingMatrixB(speciesID),&
         hartreeMatricesOUT=auxMatrix(1:numberOfSpecies)) 

       do otherSpeciesID=1, numberOfSpecies
          if(otherSpeciesID.ne.speciesID) hartreeMatrixB(speciesID,otherSpeciesID)%values=auxMatrix(otherSpeciesID)%values
          call Matrix_destructor (auxMatrix(otherSpeciesID))
       end do

       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "couplingMatrixB", nameOfSpecies
          Call Matrix_show(couplingMatrixB(speciesID))
       end if

       
       auxEnergy(:)=0.0
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
               excFileIN=trim(densFileB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridB(speciesID) )

          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "exchangeCorrelationMatrixB", nameOfSpecies
             Call Matrix_show(exchangeCorrelationMatrixB(speciesID))
          end if
          
       end if
       do otherSpeciesID=1, numberOfSpecies
          excCorrEnergyB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
       end do
       
    end do

    SUBSYSTEM_SCF_CONTINUE=.true.
    iter=0
    totalEnergy=0.0
    deltaEnergy=1.0E16
    ! do speciesID=1,numberOfSpecies
    !    totalEnergy=totalEnergy+WaveFunction_instance( speciesID )%exchangeCorrelationEnergy&
    !         +sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
    !         *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
    !         + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values))&
    !         +0.5_8*(sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
    !         * (wavefunction_instance(speciesID)%couplingMatrix%values))) 
    ! end do

    ! print *, "Complete system energy", totalEnergy

    if ( .not. CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
       write(*,*) ""
       write(*,*) "Begin Multi-Species Subsystem A SCF calculation:"
       write(*,*) ""
       write(*,*) "---------------------------------------------------------"
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          write(*,"(A20,A12,A20,A20,A20)") "Iteration", "EnergyAB","Energy Change","DensityChangeA","ParticlesInGridAB" 
       else
          write(*,"(A20,A12,A20,A20)") "Iteration", "EnergyAB","Energy Change","DensityChangeA"
       end if
       write(*,*) "---------------------------------------------------------"
    end if

    newTotalEnergy=0.0

    do while( SUBSYSTEM_SCF_CONTINUE )
       iter=iter+1
      
       !Run DFT with the total density and the subsystem density 
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then

          !!Save density matrix A, A+B
          call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileA),densityMatrixA(:))
          call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileAB),densityMatrixAB(:))
          
          !!Run DFT for subsystem A
          call system("lowdin-DFT.x SCF_DFT "//trim(densFileA))
          !!Run DFT for complete system A+B
          call system("lowdin-DFT.x SCF_DFT "//trim(densFileAB))
          
       end if
       
       !Calculates the fock matrix with the new subsystem density 
       do speciesID=1, numberOfSpecies
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
          nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
          
          !Updates two particles matrix - only Coulomb (factor=0.0)
          call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
               densityMatrixIN=densityMatrixA(speciesID),&
               factorIN=0.0_8,&
               twoParticlesMatrixOUT=hartreeMatrixA(speciesID,speciesID))

          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "CoulombMatrixA", nameOfSpecies
             Call Matrix_show(hartreeMatrixA(speciesID,speciesID))
          end if

          
          !Obtains exchange-correlation matrix from a full HF calculation (factor=1.0)
          call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
               densityMatrixIN=densityMatrixA(speciesID),&
               factorIN=1.0_8,&
               twoParticlesMatrixOUT=exchangeHFMatrixA(speciesID))

          exchangeHFMatrixA(speciesID)%values=&
               exchangeHFMatrixA(speciesID)%values-hartreeMatrixA(speciesID,speciesID)%values

          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "ExchangeMatrixA", nameOfSpecies
             Call Matrix_show(exchangeHFMatrixA(speciesID))
          end if
          
          !Updates coupling matrices
          do otherSpeciesID=1, numberOfSpecies
             call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
          end do

          call WaveFunction_buildCouplingMatrix(WaveFunction_instance,speciesID,&
               densityMatricesIN=densityMatrixA(:),&
               couplingMatrixOUT=couplingMatrixA(speciesID),&
               hartreeMatricesOUT=auxMatrix(:)) 

          do otherSpeciesID=1, numberOfSpecies
             if(otherSpeciesID.ne.speciesID) hartreeMatrixA(speciesID,otherSpeciesID)%values=auxMatrix(otherSpeciesID)%values
             call Matrix_destructor (auxMatrix(otherSpeciesID))
          end do

          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "CouplingMatrixA", nameOfSpecies
             Call Matrix_show(couplingMatrixA(speciesID))
          end if

          
          !Updates exchange correlation matrices
          if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
             auxEnergy=0.0
             call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
                  excFileIN=trim(densFileA)//".exc",&
                  exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixA(speciesID),&
                  exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
                  particlesInGridOUT=particlesInGridA(speciesID) )
             do otherSpeciesID=1, numberOfSpecies
                excCorrEnergyA(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
             end do
             
             if ( CONTROL_instance%DEBUG_SCFS ) then
                print *, "exchangeCorrelationMatrixA", nameOfSpecies
                Call Matrix_show(exchangeCorrelationMatrixA(speciesID))
             end if

             auxEnergy=0.0
             call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
                  excFileIN=trim(densFileAB)//".exc",&
                  exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixAB(speciesID),&
                  exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
                  particlesInGridOUT=particlesInGridAB(speciesID) ) 
             do otherSpeciesID=1, numberOfSpecies
                excCorrEnergyAB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
             end do
                  
             if ( CONTROL_instance%DEBUG_SCFS ) then
                print *, "exchangeCorrelationMatrixAB", nameOfSpecies
                Call Matrix_show(exchangeCorrelationMatrixAB(speciesID))
             end if
             
          end if
          
       !!Updates Fock Matrix
          OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values=WaveFunction_instance(speciesID)%hcoreMatrix%values&
               +WaveFunction_instance(speciesID)%externalPotentialMatrix%values&
               +hartreeMatrixB(speciesID,speciesID)%values&
               +exchangeHFMatrixB(speciesID)%values&
               +couplingMatrixB(speciesID)%values&
               +(exchangeCorrelationMatrixAB(speciesID)%values-exchangeCorrelationMatrixA(speciesID)%values)&
               +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
               
          OrbitalLocalizer_instance(speciesID)%fockMatrixA%values = &
               OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values &
               +hartreeMatrixA(speciesID,speciesID)%values&
               +exchangeHFMatrixA(speciesID)%values&
               +couplingMatrixA(speciesID)%values
                  
          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "fockMatrixA", nameOfSpecies
             Call Matrix_show(OrbitalLocalizer_instance(speciesID)%fockMatrixA)
          end if
          
          if ( iter > 1 )then
             call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
                  OrbitalLocalizer_instance(speciesID)%fockMatrixA, densityMatrixA(speciesID), &
                  WaveFunction_instance(speciesID)%OverlapMatrix, &
                  methodType=SCF_CONVERGENCE_DAMPING, &
                  coefficientMatrix=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, speciesID=speciesID )

             call Convergence_run( WaveFunction_instance(speciesID)%convergenceMethod )

             if ( CONTROL_instance%DEBUG_SCFS ) then
                print *, "fockMatrixA-damped", nameOfSpecies
                Call Matrix_show(OrbitalLocalizer_instance(speciesID)%fockMatrixA)
             end if

          end if
          
          !! cosmo fock matrix
          ! OrbitalLocalizer_instance(speciesID)%fockMatrixA%values= OrbitalLocalizer_instance(speciesID)%fockMatrixA%values+&
          !      0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
          !      wavefunction_instance(speciesID)%cosmo4%values)+ &
          !      wavefunction_instance(speciesID)%cosmo2%values + &
          !      wavefunction_instance(speciesID)%cosmoCoupling%values
          ! if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          !    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
          !         wavefunction_instance(speciesID)%externalPotentialMatrix%values

       end do

       !Calculates the total energy with the new subsystem density 
       newTotalEnergy=MolecularSystem_getPointChargesEnergy()

       do speciesID=1,numberOfSpecies
          !constants, depend on B, are not in Fock Matrix
          newTotalEnergy= newTotalEnergy&
               +sum(transpose(densityMatrixB(speciesID)%values)*WaveFunction_instance(speciesID)%hCoreMatrix%values) &
               +0.5*sum(transpose(densityMatrixB(speciesID)%values)*hartreeMatrixB(speciesID,speciesID)%values)&
               +0.5*sum(transpose(densityMatrixB(speciesID)%values)*exchangeHFMatrixB(speciesID)%values)&
               +0.5*sum(transpose(densityMatrixB(speciesID)%values)*couplingMatrixB(speciesID)%values)
          
          !variables, depend on A, are in Fock Matrix
          newTotalEnergy= newTotalEnergy&
               +sum(transpose(densityMatrixA(speciesID)%values)*WaveFunction_instance(speciesID)%hCoreMatrix%values) &
               +sum(transpose(densityMatrixA(speciesID)%values)*hartreeMatrixB(speciesID,speciesID)%values)&
               +sum(transpose(densityMatrixA(speciesID)%values)*exchangeHFMatrixB(speciesID)%values)&
               +sum(transpose(densityMatrixA(speciesID)%values)*couplingMatrixB(speciesID)%values)&
               +sum(transpose(densityMatrixA(speciesID)%values)*OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)&
               +0.5*sum(transpose(densityMatrixA(speciesID)%values)*hartreeMatrixA(speciesID,speciesID)%values)&               
               +0.5*sum(transpose(densityMatrixA(speciesID)%values)*exchangeHFMatrixA(speciesID)%values)&
               +0.5*sum(transpose(densityMatrixA(speciesID)%values)*couplingMatrixA(speciesID)%values)
          
          ! +(exchangeCorrelationMatrixAB(speciesID)%values-exchangeCorrelationMatrixA(speciesID)%values) is replaced by
          do otherSpeciesID=speciesID, numberOfSpecies
             newTotalEnergy= newTotalEnergy+(excCorrEnergyAB(speciesID,otherSpeciesID)-excCorrEnergyA(speciesID,otherSpeciesID))
          end do

          ! newTotalEnergy= newTotalEnergy+&
          !      sum(transpose(densityMatrixAB(speciesID)%values)* WaveFunction_instance(speciesID)%hCoreMatrix%values) &
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values)&
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)&
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* couplingMatrixA(speciesID)%values)&
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values)&
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)&
          !      +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* couplingMatrixB(speciesID)%values)&
          !      +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*sum(transpose(densityMatrixA(speciesID)%values)*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)
               ! +(sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixAB(speciesID)%values)&
               ! -sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixA(speciesID)%values))&
          ! +excCorrEnergyAB(speciesID,speciesID)-excCorrEnergyA(speciesID,speciesID)&

       end do
       ! print *, "total energy", newTotalEnergy

       !Calculates the subsystem orbitals with the new fock matrix
       do speciesID=1, numberOfSpecies
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
          nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)

          call Matrix_copyConstructor( fockMatrixTransformed, OrbitalLocalizer_instance(speciesID)%fockMatrixA )
          
          !!**********************************************************************************************
          !! Level Shifting Convergence Method       
          !!
          if ( CONTROL_instance%ACTIVATE_LEVEL_SHIFTING .eqv. .true. ) then
             
             if ( MolecularSystem_instance%species(speciesID)%isElectron ) then
                levelShiftingFactor=CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING
             else
                levelShiftingFactor=CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING
             end if
             
             fockMatrixTransformed%values = &
                  matmul( matmul( transpose( OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values ) , &
                  fockMatrixTransformed%values), OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values )

             do i=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+1, numberOfContractions
                   fockMatrixTransformed%values(i,i) = levelShiftingFactor + fockMatrixTransformed%values(i,i)
             end do

             fockMatrixTransformed%values = &
                  matmul( matmul( OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values, &
                  fockMatrixTransformed%values), transpose(OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values ))

             fockMatrixTransformed%values = &
                  matmul( matmul(WaveFunction_instance( speciesID )%OverlapMatrix%values, &
                  fockMatrixTransformed%values), transpose(WaveFunction_instance( speciesID )%OverlapMatrix%values ))                    

          end if
          !!**********************************************************************************************

          
          fockMatrixTransformed%values = &
               matmul( matmul( transpose( WaveFunction_instance(speciesID)%transformationMatrix%values ) , &
               fockMatrixTransformed%values), WaveFunction_instance(speciesID)%transformationMatrix%values )

          !! Calcula valores y vectores propios de matriz de Fock transformada.
          call Matrix_eigen( fockMatrixTransformed, OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA, &
               OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, SYMMETRIC )

          !! Calcula los  vectores propios para matriz de Fock       
          OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values = &
               matmul( WaveFunction_instance(speciesID)%transformationMatrix%values, &
               OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values )

          !! Removes contributions from basis sets of subsystem B and renormalizes orbital
          i=0
          do k = 1 , OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA + OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA
             i=i+1
             ! do mu = 1 , numberOfContractions
             !    if(OrbitalLocalizer_instance(speciesID)%basisSubsystem(mu).ne.1) &
             !         OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,i)=0.0
             ! end do
             
             normCheck=0.0
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   normCheck=normCheck+OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,i)*&
                        OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(nu,i)*&
                        WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
                end do
             end do
             ! print *, i, "normCheck", normCheck

             !! Filtra los orbitales eliminados por el umbral de solapamiento
             if ( normCheck .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD) then
                ! Shift orbital coefficients to the end of the matrix and Make energy a very large number
                do j = i , numberOfContractions-1
                   OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA%values(j)=OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA%values(j+1)
                   OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,j) = OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,j+1)
                end do
                OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA%values(numberOfContractions)=1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
                OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,numberOfContractions)=0.0
                i=i-1
             else
                !! Renormalize
                ! OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,i)=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,i)/normCheck
             end if
          end do
          
          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "waveFunctionCoefficientsA", nameOfSpecies
             Call Matrix_show(OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA)
          end if

          !Updates atomic density matrix
          call Matrix_constructor (newDensityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

          do k = 1 , OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   newDensityMatrix%values(mu,nu)=&
                        newDensityMatrix%values(mu,nu)&
                        +MolecularSystem_getEta( speciesID )*&
                        OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,k)*&
                        OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(nu,k)
                end do
             end do
          end do
          densStd(speciesID)= Matrix_standardDeviation(newDensityMatrix,densityMatrixA(speciesID))

          densityMatrixA(speciesID)=newDensityMatrix
          densityMatrixAB(speciesID)%values=densityMatrixA(speciesID)%values+densityMatrixB(speciesID)%values

          if ( CONTROL_instance%DEBUG_SCFS ) then
             print *, "densityMatrixA", nameOfSpecies
             Call Matrix_show(densityMatrixA(speciesID))

             print *, "densityMatrixAB", nameOfSpecies
             Call Matrix_show(densityMatrixAB(speciesID))
             
          end if

          
       end do
       !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
       if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
            (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
          speciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-ALPHA")  )
          otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-BETA")  )

          if(MolecularSystem_getNumberOfParticles(speciesID) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID) ) then
             OrbitalLocalizer_instance(otherSpeciesID)%waveFunctionCoefficientsA%values = OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values
             densityMatrixA(otherSpeciesID)%values=densityMatrixA(speciesID)%values
             densityMatrixAB(otherSpeciesID)%values=densityMatrixAB(speciesID)%values
          end if
       end if
       
       !Check convergence and print messages
       totalDensStd=sqrt(sum(densStd(:)**2))
       deltaEnergy=totalEnergy-newTotalEnergy
       !Writes iteration results
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
           write (*,"(I15,F20.12,F20.12,F20.12,F20.12)") iter, newTotalEnergy, deltaEnergy , totalDensStd, sum(particlesInGridAB(:)) 
       else
           write (*,"(I15,F20.12,F20.12,F20.12)") iter, newTotalEnergy, deltaEnergy, totalDensStd 
       end if

       totalEnergy=newTotalEnergy

       if ( CONTROL_instance%DEBUG_SCFS ) then
          print *, "energyContinue", deltaEnergy .gt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE, "densityContinue", &
               totalDensStd .gt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE, &
               "iterationContinue", iter .lt. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS
       end if

       if(iter .ge. CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS) then
          write (*,"(A,I4,A)")  "The number of Iterations was exceded, the convergence had failed after", iter ," global iterations"
          SUBSYSTEM_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "DENSITY" .and. &
            totalDensStd .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
          write (*,"(A,I4,A)") "Total density converged after", iter ," global iterations"
          SUBSYSTEM_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "ENERGY" .and. &
            abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE) then
          write (*,"(A,I4,A)") "Total energy converged after", iter ," global iterations"
          SUBSYSTEM_SCF_CONTINUE=.false.
       end if

       if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "BOTH" .and. & 
            abs(deltaEnergy) .lt. CONTROL_instance%TOTAL_ENERGY_TOLERANCE .and. &
            totalDensStd .lt. CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE) then
          write (*,"(A,I4,A)") "Total energy and density converged after", iter ," global iterations"
          SUBSYSTEM_SCF_CONTINUE=.false.
       end if

       
    end do

    !Forces equal coefficients for E-ALPHA and E-BETA in open shell calculations
    if ( CONTROL_instance%FORCE_CLOSED_SHELL .and. &
         (CONTROL_instance%METHOD .eq. "UKS" .or. CONTROL_instance%METHOD .eq. "UHF") ) then
       speciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-ALPHA")  )
       otherSpeciesID=MolecularSystem_getSpecieIDFromSymbol( trim("E-BETA")  )

       if(MolecularSystem_getNumberOfParticles(speciesID) .eq. MolecularSystem_getNumberOfParticles(otherSpeciesID) ) then
          OrbitalLocalizer_instance(otherSpeciesID)%waveFunctionCoefficientsA%values = OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values
          densityMatrixA(otherSpeciesID)%values=densityMatrixA(speciesID)%values
          densityMatrixAB(otherSpeciesID)%values=densityMatrixAB(speciesID)%values
       end if
    end if

    ! do speciesID=1, numberOfSpecies     

    !    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    !    call Matrix_constructor(auxMatrix(speciesID),int(numberOfContractions,8),&
    !         int(numberOfContractions-OrbitalLocalizer_instance(speciesID)%removedOrbitalsA,8),0.0_8)

    !    do i=1, numberOfContractions
    !       do j=1, numberOfContractions-OrbitalLocalizer_instance(speciesID)%removedOrbitalsA
    !          auxMatrix(speciesID)%values(i,j)=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(i,j)
    !       end do
    !    end do

    !    write(*,*) ""
    !    write(*,*) " Subsystem Eigenvectors for: ", trim(MolecularSystem_instance%species(speciesID)%name )
    !    write(*,*) "-----------------------------"
    !    write(*,*) ""

    !    call Matrix_show( auxMatrix(speciesID), &
    !         rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
    !         columnkeys = string_convertvectorofrealstostring( OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA ),&
    !         flags=WITH_BOTH_KEYS)

    ! end do

  ! Final energy evaluation - larger integration grid for DFT 

    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       !!Save density matrix A, B, A+B
       call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileA),densityMatrixA(:))
       call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileB),densityMatrixB(:))
       call WaveFunction_writeDensityMatricesToFile(WaveFunction_instance,trim(densFileAB),densityMatrixAB(:))

       write(*,*) " FINAL GRID DFT EVALUATION FOR SUBSYSTEM A: "
       write(*,*) "-----------------------------"
       call system("lowdin-DFT.x FINAL_DFT "//trim(densFileA))

       write(*,*) " FINAL GRID DFT EVALUATION FOR SUBSYSTEM B: "
       write(*,*) "-----------------------------"
       call system("lowdin-DFT.x FINAL_DFT "//trim(densFileB))

       write(*,*) " FINAL GRID DFT EVALUATION FOR SYSTEM AB: "
       write(*,*) "-----------------------------"
       call system("lowdin-DFT.x FINAL_DFT "//trim(densFileAB))
    end if

    !Calculates the final matrices with the final subsystem density 
    do speciesID=1, numberOfSpecies
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

       !Updates two particles matrix - only Coulomb (factor=0.0)
       call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
            densityMatrixIN=densityMatrixA(speciesID),&
            factorIN=0.0_8,&
            twoParticlesMatrixOUT=hartreeMatrixA(speciesID,speciesID))

       !Obtains exchange-correlation matrix from a full HF calculation (factor=1.0)
       call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),&
            densityMatrixIN=densityMatrixA(speciesID),&
            factorIN=1.0_8,&
            twoParticlesMatrixOUT=exchangeHFMatrixA(speciesID))

       exchangeHFMatrixA(speciesID)%values=&
            exchangeHFMatrixA(speciesID)%values-hartreeMatrixA(speciesID,speciesID)%values

       !Updates coupling matrices
       do otherSpeciesID=1, numberOfSpecies
          call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       end do

       call WaveFunction_buildCouplingMatrix(WaveFunction_instance,speciesID,&
            densityMatricesIN=densityMatrixA(:),&
            couplingMatrixOUT=couplingMatrixA(speciesID),&
            hartreeMatricesOUT=auxMatrix(:)) 

       do otherSpeciesID=1, numberOfSpecies
          if(otherSpeciesID.ne.speciesID) hartreeMatrixA(speciesID,otherSpeciesID)%values=auxMatrix(otherSpeciesID)%values
          call Matrix_destructor (auxMatrix(otherSpeciesID))
       end do

       !Updates exchange correlation matrices
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          auxEnergy=0.0
          call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
               excFileIN=trim(densFileA)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixA(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridA(speciesID) )
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyA(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

          auxEnergy=0.0
          call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
               excFileIN=trim(densFileB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridB(speciesID) )
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

          auxEnergy=0.0
          call WaveFunction_readExchangeCorrelationMatrix(WaveFunction_instance(speciesID),&
               excFileIN=trim(densFileAB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixAB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridAB(speciesID) ) 
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyAB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

       end if

       !!Updates Fock Matrix
       OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values=WaveFunction_instance(speciesID)%hcoreMatrix%values&
            +WaveFunction_instance(speciesID)%externalPotentialMatrix%values&
            +hartreeMatrixB(speciesID,speciesID)%values&
            +exchangeHFMatrixB(speciesID)%values&
            +couplingMatrixB(speciesID)%values&
            +(exchangeCorrelationMatrixAB(speciesID)%values-exchangeCorrelationMatrixA(speciesID)%values)&
            +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
       
       OrbitalLocalizer_instance(speciesID)%fockMatrixA%values = &
            OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values &
            +hartreeMatrixA(speciesID,speciesID)%values&
            +couplingMatrixA(speciesID)%values&
            +exchangeHFMatrixA(speciesID)%values

    end do

  !! Obtain energy compotents for whole system
     write(*,*) ""             
     write(*,*) "---------------------------------------------------------------------------"
     write(*,"(A)") " ENERGY CONTRIBUTIONS FOR THE EMBEDDED SUBSYSTEM:"
     write(*,*) "---------------------------------------------------------------------------"
     write(*,*) ""             


     write(*,*) ""             
     write(*,"(A35,A20,A20)") " COMPONENTS OF KINETIC ENERGY: ", "Subsystem A", "Subsystem B"
     write(*,*) "---------------------------------------------------------------------------"
     write(*,*) ""             

     totalKineticEnergyA=0.0
     totalKineticEnergyB=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A35,F20.12,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // " Kinetic energy = ", &
             sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%kineticMatrix%values), &
             sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%kineticMatrix%values)
        totalKineticEnergyA=totalKineticEnergyA+&
             sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%kineticMatrix%values)
        totalKineticEnergyB=totalKineticEnergyB+&
             sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%kineticMatrix%values)
     end do
     write (6,"(T10,A70)") "_________________________________________"
     write (6,"(A35,F20.12,F20.12)") "Total kinetic energy = ", totalKineticEnergyA, totalKineticEnergyB

     write(*,*) ""
     write(*,*) " COMPONENTS OF POTENTIAL ENERGY: "
     write(*,*) "---------------------------------"
     write(*,*) ""
     ! puntualMMInteractionEnergy = MolecularSystem_getMMPointChargesEnergy()
     ! if(CONTROL_instance%CHARGES_MM) then
     !    write (6,"(T10,A28,F20.12)") "Self MM potential energy   = ", puntualMMInteractionEnergy
     ! end if

     write(*,*) ""
     write(*,"(A35,A20,A20)") " Quantum/Fixed interaction energy: ", "Subsystem A", "Subsystem B"
     write(*,*) "---------------------------------------------------------------------------"
     write(*,*) ""

     totalQuantumPuntualInteractionEnergyA=0.0
     totalQuantumPuntualInteractionEnergyB=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A35,F20.12,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // "/Fixed interact. energy = ",  &
             sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%puntualInteractionMatrix%values),&
             sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%puntualInteractionMatrix%values)
        totalQuantumPuntualInteractionEnergyA=totalQuantumPuntualInteractionEnergyA+&
             sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%puntualInteractionMatrix%values)
        totalQuantumPuntualInteractionEnergyB=totalQuantumPuntualInteractionEnergyB+&
             sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%puntualInteractionMatrix%values)
     end do
     write (6,"(T10,A70)") "_________________________________________"
     write (6,"(A35,F20.12,F20.12)") "Total Q/Fixed energy = ", totalQuantumPuntualInteractionEnergyA, totalQuantumPuntualInteractionEnergyB
     
     write(*,*) ""
     write(*,"(A35,A20,A20,A20)") " Coulomb energy: ", "Subsystem A", "Subsystem B", "A-B Interaction"
     write(*,*) "---------------------------------------------------------------------------------------------"
     write(*,*) ""
     totalHartreeEnergyA=0.0
     totalHartreeEnergyB=0.0
     totalHartreeEnergyAB=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        write (6,"(A35,F20.12,F20.12,F20.12)") &
             trim( MolecularSystem_instance%species(speciesID)%name ) // &
             "/"//trim( MolecularSystem_instance%species(speciesID)%name ) // &
             " Hartree energy = ", &
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values), &
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values), &
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values)+ &
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values)
        
        totalHartreeEnergyA=totalHartreeEnergyA+0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values)
        totalHartreeEnergyB=totalHartreeEnergyB+0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values)
        totalHartreeEnergyAB=totalHartreeEnergyAB+0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values)+&
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values)        
     end do
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
        do otherSpeciesID = speciesID + 1, MolecularSystem_instance%numberOfQuantumSpecies                
           write (6,"(A35,F20.12,F20.12,F20.12)") &
                trim( MolecularSystem_instance%species(speciesID)%name ) // &
                "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                " Hartree energy = ", &
                sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values), &
                sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values), &
                sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values)+ &
                sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values)

           totalHartreeEnergyA=totalHartreeEnergyA+sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values)
           totalHartreeEnergyB=totalHartreeEnergyB+sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values)
           totalHartreeEnergyAB=totalHartreeEnergyAB+sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values)+&
                sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values)                   
        end do
     end do
     write (6,"(T10,A90)") "_____________________________________________________________"
     write (6,"(A35,F20.12,F20.12,F20.12)") "Total Hartree energy = ", totalHartreeEnergyA, totalHartreeEnergyB, totalHartreeEnergyAB


     write(*,*) ""
     write(*,"(A35,A20,A20,A20)") " Exchange(HF) energy: ", "Subsystem A", "Subsystem B", "A-B Interaction"
     write(*,*) "---------------------------------------------------------------------------------------------"
     write(*,*) ""
     totalExchangeHFEnergyA=0.0
     totalExchangeHFEnergyB=0.0
     totalExchangeHFEnergyAB=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
           write (6,"(A35,F20.12,F20.12,F20.12)") &
             trim( MolecularSystem_instance%species(speciesID)%name ) // &
             " Exchange energy = ", &
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixA(speciesID)%values), &
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixB(speciesID)%values), &
             sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)
           !0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
           
           totalExchangeHFEnergyA=totalExchangeHFEnergyA+0.5*sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
           totalExchangeHFEnergyB=totalExchangeHFEnergyB+0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)
           totalExchangeHFEnergyAB=totalExchangeHFEnergyAB+sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)! + &
                ! 0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
        end do

     write (6,"(T10,A90)") "_____________________________________________________________"
     write (6,"(A35,F20.12,F20.12,F20.12)") "Total Exchange energy = ", totalExchangeHFEnergyA, totalExchangeHFEnergyB, totalExchangeHFEnergyAB

     write(*,*) ""

     totalExchangeCorrelationEnergyA=0.0
     totalExchangeCorrelationEnergyB=0.0
     totalExchangeCorrelationEnergyAB=0.0
     totalEmbeddingPotentialEnergyA=0.0
     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        write(*,*) ""
        write(*,"(A35,A20,A20,A20)") " Exchange-Correlation(DFT) energy: ", "(Subsystem A)*", "Subsystem B", "A-B Interaction"
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
           write (6,"(A35,F20.12,F20.12,F20.12)") &
                trim( MolecularSystem_instance%species(speciesID)%name ) // &
                " Exc.Corr. energy = ", &
                excCorrEnergyA(speciesID,speciesID), excCorrEnergyB(speciesID,speciesID), &
                excCorrEnergyAB(speciesID,speciesID)-excCorrEnergyA(speciesID,speciesID)-excCorrEnergyB(speciesID,speciesID)

           totalExchangeCorrelationEnergyA=totalExchangeCorrelationEnergyA+excCorrEnergyA(speciesID,speciesID)
           totalExchangeCorrelationEnergyB=totalExchangeCorrelationEnergyB+excCorrEnergyB(speciesID,speciesID)
           totalExchangeCorrelationEnergyAB=totalExchangeCorrelationEnergyAB+&
                excCorrEnergyAB(speciesID,speciesID)-excCorrEnergyA(speciesID,speciesID)-excCorrEnergyB(speciesID,speciesID)
        end do
        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
           do otherSpeciesID = speciesID + 1, MolecularSystem_instance%numberOfQuantumSpecies                
              write (6,"(A35,F20.12,F20.12,F20.12)") &
                   trim( MolecularSystem_instance%species(speciesID)%name ) // &
                   "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                   " Corr. energy = ", &
                   excCorrEnergyA(speciesID,otherSpeciesID), excCorrEnergyB(speciesID,otherSpeciesID), &
                   excCorrEnergyAB(speciesID,otherSpeciesID)-excCorrEnergyA(speciesID,otherSpeciesID)-excCorrEnergyB(speciesID,otherSpeciesID)

              totalExchangeCorrelationEnergyA=totalExchangeCorrelationEnergyA+excCorrEnergyA(speciesID,otherSpeciesID)
              totalExchangeCorrelationEnergyB=totalExchangeCorrelationEnergyB+excCorrEnergyB(speciesID,otherSpeciesID)
              totalExchangeCorrelationEnergyAB=totalExchangeCorrelationEnergyAB+&
                   excCorrEnergyAB(speciesID,otherSpeciesID)-excCorrEnergyA(speciesID,otherSpeciesID)-excCorrEnergyB(speciesID,otherSpeciesID)
           end do
        end do
        
        write (6,"(T10,A90)") "_____________________________________________________________"
        write (6,"(A35,F20.12,F20.12,F20.12)") "Total Exchange Correlation energy = ", totalExchangeCorrelationEnergyA, totalExchangeCorrelationEnergyB, totalExchangeCorrelationEnergyAB
        write(*,*) ""
        write(*,*) "* Total energy does not include subsystem A DFT exchange-correlation"
        write(*,*) ""
        write(*,"(A35,A20,A20,A20)") " A Embedding in B Potential Energy: ", "(A-B Interaction)*"
        write(*,*) "-----------------------------------------------------------"
        write(*,*) ""
        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
           write (6,"(A35,F20.12)") &
                trim( MolecularSystem_instance%species(speciesID)%name ) // &
                " Embedding energy = ", &
                sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixAB(speciesID)%values)&
                -sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixA(speciesID)%values)

           totalEmbeddingPotentialEnergyA=totalEmbeddingPotentialEnergyA&
                +sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixAB(speciesID)%values)&
                -sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixA(speciesID)%values)

        end do
        write (6,"(T10,A50)") "_____________________"
        write (6,"(A35,F20.12,F20.12,F20.12)") "Total Embedding Potential energy = ", totalEmbeddingPotentialEnergyA
        write(*,*) ""
        write(*,*) "* This energy term is added to the core operator and removed from the total energy"
        write(*,*) ""

     end if

     write(*,*) ""
     write(*,"(A35,A20,A20,A20)") " A/B Projection Operator Correction: ", "A-B Interaction"
     write(*,*) "-----------------------------------------------------------"
     write(*,*) ""
     totalProjectionCorrectionA=0.0
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
        write (6,"(A35,F20.12)") &
             trim( MolecularSystem_instance%species(speciesID)%name ) // &
             " Projection energy = ", &
             OrbitalLocalizer_instance(speciesID)%levelShiftingValue*&
             sum(transpose(densityMatrixA(speciesID)%values)*&
             OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)

        totalProjectionCorrectionA=totalProjectionCorrectionA+&
             OrbitalLocalizer_instance(speciesID)%levelShiftingValue*&
             sum(transpose(densityMatrixA(speciesID)%values)*&
             OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)

     end do
     write (6,"(T10,A50)") "_____________________"
     write (6,"(A35,F20.12,F20.12,F20.12)") "Total Projection energy correction  = ", totalProjectionCorrectionA


     totalExternalPotentialEnergyA=0.0
     totalExternalPotentialEnergyB=0.0
     if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

        write(*,*) ""
        write(*,"(A35,A20,A20)") " External Potential energy: ", "Subsystem A", "Subsystem B"
        write(*,*) "---------------------------------------------------------------------------"
        write(*,*) ""

        do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
           write (6,"(A35,F20.12,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // "/Ext. Pot. energy = ",  &
                sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%externalPotentialMatrix%values),&
                sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%externalPotentialMatrix%values)
           totalExternalPotentialEnergyA=totalExternalPotentialEnergyA+&
                sum(transpose(densityMatrixA(speciesID)%values)* WaveFunction_instance( speciesID )%externalPotentialMatrix%values)
           totalExternalPotentialEnergyB=totalExternalPotentialEnergyB+&
                sum(transpose(densityMatrixB(speciesID)%values)* WaveFunction_instance( speciesID )%externalPotentialMatrix%values)
        end do
        write (6,"(T10,A70)") "_________________________________________"
        write (6,"(A35,F20.12,F20.12)") "Total External Potential energy = ", totalExternalPotentialEnergyA, totalExternalPotentialEnergyB        
     end if

     ! potentialEnergy = totalRepulsionEnergy &
     !      + puntualInteractionEnergy &
     !      + totalQuantumPuntualInteractionEnergy &
     !      + totalHartreeEnergy &
     !      + totalExchangeHFEnergy &
     !      + totalExchangeCorrelationEnergy &
     !      + totalExternalPotentialEnergy

     ! totalCosmoEnergy = sum( WaveFunction_instance(:)%cosmoEnergy)

     ! if(CONTROL_instance%COSMO) then
     !    write(*,*)"totalCosmoEnergy",totalCosmoEnergy
     !    write(*,*)"cosmo3energy",cosmo3Energy

     !    potentialEnergy=potentialEnergy+totalCosmoEnergy+cosmo3Energy

     ! end if

     
     
     write(*,*) ""
     write(*,*) " TOTAL ENERGY COMPONENTS: "
     write(*,*) "=================="
     write(*,*) ""
     puntualInteractionEnergy = MolecularSystem_getPointChargesEnergy()
     write(*,"(A35,F20.12)") "Fixed potential energy = ", puntualInteractionEnergy
     write(*,"(A35,F20.12)") "TOTAL SUBSYSTEM A ENERGY = ", totalKineticEnergyA+totalQuantumPuntualInteractionEnergyA+&
          totalExternalPotentialEnergyA+&
          totalHartreeEnergyA+totalExchangeHFEnergyA
     write(*,"(A35,F20.12)") "TOTAL SUBSYSTEM B ENERGY = ", totalKineticEnergyB+totalQuantumPuntualInteractionEnergyB+&
          totalExternalPotentialEnergyB+&
          totalHartreeEnergyB+totalExchangeHFEnergyB+totalExchangeCorrelationEnergyB
     write(*,"(A35,F20.12)") "TOTAL A-B INTERACTION ENERGY = ", totalHartreeEnergyAB+totalExchangeHFEnergyAB+totalExchangeCorrelationEnergyAB+totalProjectionCorrectionA
     !+totalEmbeddingPotentialEnergyA

     totalKineticEnergy=totalKineticEnergyA+totalKineticEnergyB
     totalPotentialEnergy=puntualInteractionEnergy+&
          totalQuantumPuntualInteractionEnergyA+&
          totalQuantumPuntualInteractionEnergyB+&
          totalExternalPotentialEnergyA+&
          totalExternalPotentialEnergyB+&
          totalHartreeEnergyA+&
          totalHartreeEnergyB+&
          totalHartreeEnergyAB+&
          totalExchangeHFEnergyA+&
          totalExchangeHFEnergyB+&
          totalExchangeHFEnergyAB+&
          totalExchangeCorrelationEnergyB+&
          totalExchangeCorrelationEnergyAB+&          
          totalProjectionCorrectionA

     totalEnergy=totalKineticEnergy+totalPotentialEnergy

     ! totalEmbeddingPotentialEnergyA+&

     write(*,"(T10,A50)") "_____________________"
     write(*,"(A35,F20.12)") "TOTAL ENERGY = ", totalEnergy
     write(*,*) ""
     write(*,"(A35,F20.12)") "TOTAL KINETIC ENERGY = ", totalKineticEnergy
     write(*,"(A35,F20.12)") "TOTAL POTENTIAL ENERGY = ", totalPotentialEnergy
     write(*,*) ""
     write(*,"(A35,F20.12)") "TOTAL VIRIAL RATIO (V/T) = ", - (totalPotentialEnergy / totalKineticEnergy)
     write(*,*) ""
     write(*,*) ""
     write(*,*) " END ENERGY COMPONENTS"
     write(*,*) ""

    write(*,*) ""
    write(*,*) "--------------------------------------"
    write(*,*) "Max Deviation From A-B Orthogonality"
    write(*,*) ""
    write(*,"(A10,A20,A20,A10)") "Species", "Original B orbital", "Subsystem A orbital", "<A|B>"
    do speciesID=1, numberOfSpecies
       maxSumAB=0.0
       ii=0
       jj=0
       occupationNumber = OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       if(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA .gt. 0) then
          do j = 1 , occupationNumber
             if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(j).ne.1) then
                do i = 1, OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA
                   sumAB=0.0
                   do mu = 1 , numberOfContractions
                      do nu = 1 , numberOfContractions
                         sumAB=sumAB+&
                              OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,i)*&
                              WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(nu,j)*&
                              WaveFunction_instance(speciesID)%OverlapMatrix%values(mu,nu)
                      end do
                   end do
                   if( (abs(sumAB) .gt. abs(maxSumAB))) then
                      ii=i
                      jj=j
                      maxSumAB=sumAB
                   end if
                end do
             end if
          end do
          write(*,"(A10,I20,I20,ES15.3)") trim(MolecularSystem_instance%species(speciesID)%name ), jj, ii, maxSumAB
       end if
       ! kAB=0.0
       ! pAB=0.0
       ! if(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA .gt. 0) then
       !    do j = 1 , occupationNumber
       !       if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(j).ne.1) then
       !          do i = 1, OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA
       !             sumAB=0.0
       !             do mu = 1 , numberOfContractions
       !                do nu = 1 , numberOfContractions
       !                   kAB=kAB+&
       !                        OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,i)*&
       !                        WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(nu,j)*&
       !                        WaveFunction_instance(speciesID)%kineticMatrix%values(mu,nu)
       !                   pAB=pAB+&
       !                        OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(mu,i)*&
       !                        WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(nu,j)*&
       !                        WaveFunction_instance(speciesID)%puntualInteractionMatrix%values(mu,nu)
       !                end do
       !             end do
       !          end do
       !       end if
       !    end do
       !    write(*,"(A10,ES15.3,A10,ES15.3)") "kAB", kAB, "pAB", pAB
       ! end if
    end do
    write(*,*) "--------------------------------------"

    write(*,*) ""             
    write(*,*) "---------------------------------------------------------------------------"
    write(*,"(A)") " EMBEDDED SUBSYSTEM INFORMATION "
    write(*,*) "---------------------------------------------------------------------------"
    write(*,*) ""             


    !Destroy unrequired basis centers and shells from the molecular system
    !Write a LOWDIN.BAS with the reduced basis set file
    
    open(unit=40, file="subsystem.bas", status="replace", form="formatted")

    write(40,*) numberOfSpecies

    do speciesID=1, numberOfSpecies
       write(40,*) MolecularSystem_instance%species(speciesID)%name

       numberOfCenters=0
       nu=0
       do i=1, size(MolecularSystem_instance%species(speciesID)%particles)
          includeAtom=.false.
          do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             nu=nu+1
             if(shellPopulationA(speciesID)%values(nu) .gt. OrbitalLocalizer_instance(speciesID)%basisThreshold .or. &
                  MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%subsystem .eq. 1 ) includeAtom=.true.
          end do
          if( includeAtom ) numberOfCenters=numberOfCenters+1
       end do

       write(40,*) numberOfCenters
       nu=0
       do i=1, size(MolecularSystem_instance%species(speciesID)%particles)

          oldNu=nu
          numberOfShells=0
          includeAtom=.false.
          do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             nu=nu+1
             if(shellPopulationA(speciesID)%values(nu) .gt. OrbitalLocalizer_instance(speciesID)%basisThreshold .or. &
                  MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%subsystem .eq. 1 ) then
                includeAtom=.true.
                numberOfShells=numberOfShells+1
             end if
          end do
          nu=oldNu
          
          if( includeAtom ) then
             write(40,*) MolecularSystem_instance%species(speciesID)%particles(i)%nickname
             write(40,*) numberOfShells
             do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
                nu=nu+1
                if(shellPopulationA(speciesID)%values(nu) .gt. OrbitalLocalizer_instance(speciesID)%basisThreshold .or. &
                     MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%subsystem .eq. 1 ) then
                   write(40,*) &
                        MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%angularMoment, &
                        MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%length
                   write(40,*) &
                        MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%origin
                   write(40,*) &
                        MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%orbitalExponents
                   write(40,*) &
                        MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%contractionCoefficients
                end if
             end do

          else
             nu=nu+size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
          end if

       end do
    end do
    
    allocate(holdPointChargeSymbol(MolecularSystem_instance%numberOfPointCharges),holdPointChargeMass(MolecularSystem_instance%numberOfPointCharges))
    write(40,*) MolecularSystem_instance%numberOfPointCharges
    do i = 1, MolecularSystem_instance%numberOfPointCharges
       holdPointChargeSymbol(i)=trim(MolecularSystem_instance%pointCharges(i)%nickname)
       holdPointChargeMass(i)=MolecularSystem_instance%pointCharges(i)%mass
       write(40,*) MolecularSystem_instance%pointCharges(i)%charge
       write(40,*) MolecularSystem_instance%pointCharges(i)%origin
    end do
    
    close(40)

    call MolecularSystem_loadFromFile( "LOWDIN.BAS", "subsystem" )

    !! Fix missing information from file loading

    do speciesID=1, numberOfSpecies
       
       MolecularSystem_instance%species(speciesID)%particles(:)%basis%name="REDUCED"      

       do i=1, size(MolecularSystem_instance%species(speciesID)%particles)
          MolecularSystem_instance%species(speciesID)%particles(i)%name=MolecularSystem_instance%species(speciesID)%name
          MolecularSystem_instance%species(speciesID)%particles(i)%symbol=MolecularSystem_instance%species(speciesID)%symbol
          MolecularSystem_instance%species(speciesID)%particles(i)%statistics=MolecularSystem_instance%species(speciesID)%statistics
          MolecularSystem_instance%species(speciesID)%particles(i)%basisSetName="REDUCED"
          MolecularSystem_instance%species(speciesID)%particles(i)%origin=MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(1)%origin
          MolecularSystem_instance%species(speciesID)%particles(i)%charge=MolecularSystem_instance%species(speciesID)%charge
          MolecularSystem_instance%species(speciesID)%particles(i)%mass=MolecularSystem_instance%species(speciesID)%mass
          MolecularSystem_instance%species(speciesID)%particles(i)%spin=MolecularSystem_instance%species(speciesID)%spin
          MolecularSystem_instance%species(speciesID)%particles(i)%totalCharge=0
          MolecularSystem_instance%species(speciesID)%particles(i)%internalSize=0
          MolecularSystem_instance%species(speciesID)%particles(i)%isQuantum=.true.
          MolecularSystem_instance%species(speciesID)%particles(i)%isDummy=.false.
          MolecularSystem_instance%species(speciesID)%particles(i)%fixComponent=.false.
          MolecularSystem_instance%species(speciesID)%particles(i)%isCenterOfOptimization=.true.
          MolecularSystem_instance%species(speciesID)%particles(i)%multiplicity=MolecularSystem_instance%species(speciesID)%multiplicity
          MolecularSystem_instance%species(speciesID)%particles(i)%id=i
          MolecularSystem_instance%species(speciesID)%particles(i)%owner=i
          MolecularSystem_instance%species(speciesID)%particles(i)%basisSetSize=size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
          MolecularSystem_instance%species(speciesID)%particles(i)%subsystem=1
          ! MolecularSystem_instance%species(speciesID)%particles(i)%childs=.false.
          
          ! do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
          !    print *, "particle", i, "contraction", j
          !    print *, "angularMoment ", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%angularMoment
          !    print *, "numCartesianOrbitals ", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
          !    print *, "length ", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%length
          !    print *, "origin ", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%origin
          !    print *, "exponents", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%orbitalExponents
          !    print *, "coefficients", MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%contractionCoefficients
          ! end do
       end do
    end do

    do i = 1, MolecularSystem_instance%numberOfPointCharges
       MolecularSystem_instance%pointCharges(i)%name=trim(holdPointChargeSymbol(i))
       MolecularSystem_instance%pointCharges(i)%symbol=trim(holdPointChargeSymbol(i))
       MolecularSystem_instance%pointCharges(i)%statistics="BOSON"
       MolecularSystem_instance%pointCharges(i)%basisSetName="DIRAC"
       MolecularSystem_instance%pointCharges(i)%nickname=trim(holdPointChargeSymbol(i))
       MolecularSystem_instance%pointCharges(i)%mass=holdPointChargeMass(i)
       MolecularSystem_instance%pointCharges(i)%totalCharge=MolecularSystem_instance%pointCharges(i)%charge
       MolecularSystem_instance%pointCharges(i)%isQuantum=.false.
    end do

    call ParticleManager_setOwner()

    do speciesID=1, numberOfSpecies
       do i=1, size(MolecularSystem_instance%species(speciesID)%particles)
          MolecularSystem_instance%species(speciesID)%particles(i)%totalCharge=&
               molecularSystem_instance%allParticles( MolecularSystem_instance%species(speciesID)%particles(i)%owner )%particlePtr%charge
          MolecularSystem_instance%species(speciesID)%particles(i)%internalSize=&
               abs(molecularSystem_instance%allParticles( MolecularSystem_instance%species(speciesID)%particles(i)%owner )%particlePtr%charge)
       end do
    end do
    ! call MolecularSystem_showInformation()
    call MolecularSystem_showParticlesInformation()

    call MolecularSystem_showCartesianMatrix(molecularSystem_instance,1)
    
    !Save molecular system to file
    call MolecularSystem_saveToFile()
    
    !! Recalculate one particle integrals  
    call system("rm lowdin.opints")
    call system("lowdin-ints.x ONE_PARTICLE")

    !! Start the wavefunction object
    deallocate(WaveFunction_instance)
    call WaveFunction_constructor(WaveFunction_instance,numberOfSpecies)
       
    do speciesID=1, numberOfSpecies
       call WaveFunction_readOverlapMatrix(WaveFunction_instance(speciesID), "lowdin.opints")
    end do
    
    !!Now, we are resizing the matrices

    do speciesID=1, numberOfSpecies
            
       numberOfContractions=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+&
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB+&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA+&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB

       k=0
       do i = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          k=k+1
          ! OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA+&
          !   OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA

          wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(k)=OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA%values(k)

          mu=0
          do oldMu = 1 , numberOfContractions
             if( OrbitalLocalizer_instance(speciesID)%basisSubsystem(oldMu) .eq. 1) then
                mu=mu+1
                wavefunction_instance(speciesID)%waveFunctionCoefficients%values(mu,k)=&
                     OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(oldMu,k)
             end if
          end do
          
          normCheck=0.0
          do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
             do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                normCheck=normCheck+wavefunction_instance(speciesID)%waveFunctionCoefficients%values(mu,k)*&
                     wavefunction_instance(speciesID)%waveFunctionCoefficients%values(nu,k)*&
                     WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
             end do
          end do
          ! print *, k, "normCheck", normCheck

          if ( normCheck .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD) then
             ! Shift orbital coefficients to the end of the matrix and Make energy a very large number
             do j = k , numberOfContractions-1
                wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(j)=OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA%values(j+1)
                wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,j) = OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(:,j+1)
             end do
             wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(numberOfContractions)=1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
             wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,numberOfContractions)=0.0
             k=k-1
          else
          !! Renormalize
             wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,k)=wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,k)/normCheck
          end if
       end do
      
       do k = 1 , MolecularSystem_getOcupationNumber(speciesID)
          do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
             do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                wavefunction_instance(speciesID)%densityMatrix%values(mu,nu)=&
                     wavefunction_instance(speciesID)%densityMatrix%values(mu,nu)&
                     +MolecularSystem_getEta(speciesID)*&
                     wavefunction_instance(speciesID)%waveFunctionCoefficients%values(mu,k)*&
                     wavefunction_instance(speciesID)%waveFunctionCoefficients%values(nu,k)
             end do
          end do
       end do

       mu=0
       do oldMu = 1 , numberOfContractions
          if(OrbitalLocalizer_instance(speciesID)%basisSubsystem(oldMu) .eq. 1) then
             mu=mu+1
             nu=0
             do oldNu = 1 , numberOfContractions
                if( OrbitalLocalizer_instance(speciesID)%basisSubsystem(oldNu) .eq. 1) then
                   nu=nu+1
                   WaveFunction_instance(speciesID)%hcoreMatrix%values(mu,nu)=&
                        OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values(oldMu,oldNu)
                   ! WaveFunction_instance(speciesID)%twoParticlesMatrix%values(mu,nu)=&
                   !      hartreeMatrixA(speciesID,speciesID)%values(oldMu,oldNu)+exchangeHFMatrixA(speciesID)%values(oldMu,oldNu)
                   ! WaveFunction_instance(speciesID)%couplingMatrix%values(mu,nu)=&
                   !      couplingMatrixA(speciesID)%values(oldMu,oldNu)
                   ! WaveFunction_instance(speciesID)%fockMatrix%values(mu,nu)=&
                   !      OrbitalLocalizer_instance(speciesID)%fockMatrixA%values(oldMu,oldNu)
                end if
             end do
          end if
       end do      
    end do
        

    wfnUnit = 300
    wfnFile = "lowdin.wfn"
    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    rewind(wfnUnit)
    do speciesID = 1, numberOfSpecies
       labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
       labels(1) = "DENSITY"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )
    end do
    close(wfnUnit)
    
    if ( trim(CONTROL_instance%INTEGRAL_STORAGE) == "DIRECT" ) then
       do speciesID=1, numberOfSpecies
          if (Libint2Instance(speciesID)%isInstanced) call Libint2Interface_constructor(Libint2Instance(speciesID), MolecularSystem_instance, speciesID)
       end do
    else
       call system("rm  *.ints")   
       if( CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL ) then
          call system(" lowdin-ints.x TWO_PARTICLE_G12")
       else        
          call system(" lowdin-ints.x TWO_PARTICLE_R12")
       end if
    end if
    
    !! Recalculate matrices - not really needed, but do it anyway to be sure
    do speciesID=1, numberOfSpecies
       call WaveFunction_readOverlapMatrix(WaveFunction_instance(speciesID), "lowdin.opints")
       call WaveFunction_buildTransformationMatrix(WaveFunction_instance(speciesID), 2 )
       call WaveFunction_buildTwoParticlesMatrix(WaveFunction_instance(speciesID),factorIN=1.0_8)
       call WaveFunction_buildCouplingMatrix(WaveFunction_instance,speciesID)
       call WaveFunction_buildFockMatrix(WaveFunction_instance(speciesID))
    end do

    !! Molecular orbital fock operator expected value
    do speciesID=1, numberOfSpecies
       ! print *, "Species", speciesID, "eigenvalues" 
       do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
          if(wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(i) .eq. 1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) cycle
          wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)=0.0
          do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
             do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
                wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)=&
                     wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)+&
                     wavefunction_instance(speciesID)%fockMatrix%values(mu,nu)*&
                     wavefunction_instance(speciesID)%waveFunctionCoefficients%values(mu,i)*&
                     wavefunction_instance(speciesID)%waveFunctionCoefficients%values(nu,i)
             end do
          end do
          ! print *, i, wavefunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)
       end do
    end do

    
    write(*,*) ""
    write(*,*) " SUBSYSTEM EIGENVALUES AND EIGENVECTORS: "
    write(*,*) "======================================== "
    write(*,*) ""

    ! if ( CONTROL_instance%HF_PRINT_EIGENVALUES ) then
    do speciesID=1, numberOfSpecies
       write(*,*) ""
       write(*,*) " Subsystem Eigenvalues for: ", trim( MolecularSystem_instance%species(speciesID)%name )
       write(*,*) "--------------------------- "
       write(*,*) ""
       WaveFunction_instance(speciesID)%removedOrbitals=0
       do i = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          if(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i) .gt. sqrt(OrbitalLocalizer_instance(speciesID)%levelShiftingValue) ) then
             WaveFunction_instance(speciesID)%removedOrbitals=WaveFunction_instance(speciesID)%removedOrbitals+1
             write(*,"(T1,A1,I4,F24.12)") "*",i,WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)
          else
             write(*,"(T2,I4,F24.12)") i,WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)
          end if
       end do
       write(*,"(I10,A30,A10)") WaveFunction_instance(speciesID)%removedOrbitals, " eigenvectors were removed for ", trim(MolecularSystem_instance%species(speciesID)%name )
    end do
    write(*,*) " end of eigenvalues "
    ! end if
    
    if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL" .or. trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then
       
       do speciesID=1, numberOfSpecies

          if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "ALL") then

             write(*,*) ""
             write(*,*) " Subsystem Eigenvectors for: ", trim(MolecularSystem_instance%species(speciesID)%name )
             write(*,*) "-----------------------------"
             write(*,*) ""

             call Matrix_constructor(auxMatrix(speciesID),int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),&
                  int(MolecularSystem_getTotalNumberOfContractions(speciesID)-WaveFunction_instance(speciesID)%removedOrbitals,8),0.0_8)

             do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
                do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)-WaveFunction_instance(speciesID)%removedOrbitals
                   auxMatrix(speciesID)%values(mu,nu)=wavefunction_instance(speciesID)%waveFunctionCoefficients%values(mu,nu)
                end do
             end do

             call Matrix_show( auxMatrix(speciesID), &
                  rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
                  columnkeys = string_convertvectorofrealstostring( wavefunction_instance(speciesID)%molecularOrbitalsEnergy ),&
                  flags=WITH_BOTH_KEYS)

           else if ( trim(CONTROL_instance%HF_PRINT_EIGENVECTORS) .eq. "OCCUPIED" ) then
             
              write(*,*) ""
              write(*,*) " Subsystem Occupied Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
              write(*,*) "------------------------------------- "
              write(*,*) ""
              
             call Matrix_constructor(auxMatrix(speciesID),int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),&
                  int(MolecularSystem_getOcupationNumber(speciesID),8),0.0_8)
              do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
                 do j=1, MolecularSystem_getOcupationNumber(speciesID)
                    auxMatrix(speciesID)%values(i,j)=WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(i,j)
                 end do
              end do
              
              call Matrix_show( auxMatrix(speciesID), &
                   rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
                   columnkeys = string_convertvectorofrealstostring( wavefunction_instance(speciesID)%molecularOrbitalsEnergy ),&
                   flags=WITH_BOTH_KEYS)

           end if
        end do
        write(*,*) ""
        write(*,*) " end of eigenvectors "
     end if
    write(*,*) ""
    write(*,*) " Subsystem energy summary in the reduced basis: "
    write(*,*) "------------------------------------------------"
    write(*,*) ""
    write(*,"(A20,A20,A20,A20)") "Species:", "embedded hcore", "two particles", "coupling"
    write(*,"(A20,A20,A20,A20)") "--------", "--------------", "-------------", "--------"
    write(*,*) ""

    totalEnergy=0.0
    do speciesID=1, numberOfSpecies
       
       write (*,"(A20,F20.12,F20.12,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name )//":", &
            sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)* wavefunction_instance(speciesID)%hcoreMatrix%values), &
            0.5*sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)* WaveFunction_instance(speciesID)%twoParticlesMatrix%values), &
            sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)* WaveFunction_instance(speciesID)%couplingMatrix%values)
       
       totalEnergy=totalEnergy+sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)* wavefunction_instance(speciesID)%hcoreMatrix%values)+&
            0.5*sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)* WaveFunction_instance(speciesID)%twoParticlesMatrix%values)

       do otherSpeciesID=speciesID+1, numberOfSpecies
          totalEnergy=totalEnergy+ sum(transpose(wavefunction_instance(speciesID)%densityMatrix%values)*wavefunction_instance(speciesID)%hartreeMatrix(otherSpeciesID)%values)
       end do
    end do
    
    write(*,*) ""
    write(*,"(A35,F20.12)") "Fixed potential energy = ", puntualInteractionEnergy
    write(*,"(A35,F20.12)") "Subsystem A energy = ", totalEnergy
    write(*,"(A35,F20.12)") "Subsystem B energy = ", totalKineticEnergyB+totalQuantumPuntualInteractionEnergyB+&
         totalHartreeEnergyB+totalExchangeHFEnergyB+totalExchangeCorrelationEnergyB
    write(*,"(A35,F20.12)") "Embedding Correction Energy A = ", totalExchangeCorrelationEnergyAB-totalEmbeddingPotentialEnergyA
    totalEnergy=totalEnergy+puntualInteractionEnergy+totalKineticEnergyB+totalQuantumPuntualInteractionEnergyB+&
         totalHartreeEnergyB+totalExchangeHFEnergyB+totalExchangeCorrelationEnergyB+totalExchangeCorrelationEnergyAB-totalEmbeddingPotentialEnergyA
    write(*,"(T10,A50)") "_____________________"
    write(*,"(A35,F20.12)") "TOTAL ENERGY = ", totalEnergy
    write(*,*) ""
    
   
    !!**********************************************************
    !! Save matrices to subsystem lowdin.wfn file
    !!
    wfnUnit = 300
    wfnFile = "lowdin.wfn"

    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    rewind(wfnUnit)

    labels = ""

    do speciesID = 1, numberOfSpecies

       labels(2) = MolecularSystem_getNameOfSpecies(speciesID)

       labels(1) = "REMOVED-ORBITALS"
       call Vector_writeToFile(unit=wfnUnit, binary=.true., value=real(WaveFunction_instance(speciesID)%removedOrbitals,8), arguments= labels )

       labels(1) = "TWOPARTICLES"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%twoParticlesMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       labels(1) = "COUPLING"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%couplingMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       ! labels(1) = "EXCHANGE-CORRELATION"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%exchangeCorrelationMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       ! labels(1) = "EXCHANGE-CORRELATION-ENERGY"
       ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=WaveFunction_instance(speciesID)%exchangeCorrelationEnergy, arguments= labels )

       labels(1) = "COEFFICIENTS"
       call Matrix_writeToFile(wavefunction_instance(speciesID)%waveFunctionCoefficients, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "DENSITY"
       call Matrix_writeToFile(wavefunction_instance(speciesID)%densityMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "HCORE"      
       call Matrix_writeToFile(wavefunction_instance(speciesID)%hcoreMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "ORBITALS"
       call Vector_writeToFile(wavefunction_instance(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "FOCK"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%fockMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "OVERLAP"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%overlapMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "TRANSFORMATION"
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%transformationMatrix, unit=wfnUnit, binary=.true., arguments = labels )

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
       vecFile = trim(CONTROL_instance%INPUT_FILE)//"subvec"
       open(unit=vecUnit, file=trim(vecFile), form="unformatted", status='replace')
       do speciesID = 1, numberOfSpecies
          labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
          labels(1) = "COEFFICIENTS"
          call Matrix_writeToFile(WaveFunction_instance(speciesID)%waveFunctionCoefficients, &
               unit=vecUnit, binary=.true., arguments = labels)

          labels(1) = "ORBITALS"
          call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, & 
               unit=vecUnit, binary=.true., arguments = labels )
       end do

    else
       vecFile = trim(CONTROL_instance%INPUT_FILE)//"plainsubvec"
       open(unit=vecUnit, file=trim(vecFile), form="formatted", status='replace')

       do speciesID = 1, numberOfSpecies
          labels(2) = MolecularSystem_getNameOfSpecies(speciesID)
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
    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%cosmo3Energy, arguments=["COSMO3ENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalCouplingEnergy, arguments=["COUPLINGENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%electronicRepulsionEnergy, arguments=["COUPLING-E-"])

    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MolecularSystem_getPointChargesEnergy()&
         +totalKineticEnergyB&
         +totalQuantumPuntualInteractionEnergyB&
         +totalHartreeEnergyB&
         +totalExchangeHFEnergyB&
         +totalExchangeCorrelationEnergyB&
         , arguments=["PUNTUALINTERACTIONENERGY"])
    close(wfnUnit)       
    
    
    !!calculate HF/KS properties
    call system ("lowdin-CalcProp.x")
    
  end subroutine OrbitalLocalizer_levelShiftSubsystemOrbitals

  ! subroutine OrbitalLocalizer_reorderSubsystemOrbitals(speciesID,coefficients,fockMatrix,eigenvalues)
  !   integer :: speciesID
  !   type(Matrix) :: coefficients
  !   type(Matrix) :: fockMatrix
  !   type(Vector) :: eigenvalues
    
  !   type(Matrix) :: overlapMatrix
  !   type(Matrix) :: densityMatrix
  !   type(Matrix) :: miniDensityMatrix
  !   type(Vector) :: energyContribution
  !   real(8) :: sumSubsystem(2), sumAB, normCheck, holdEnergy
  !   real(8),allocatable ::holdValues(:)

  !   integer :: numberOfContractions, occupationNumber
  !   integer :: index1,index2,holdIndex,holdSystem
  !   integer :: mu,nu,alpha,beta, i, ii, j, jj, k

  !   !Removes the projection matrix from the fock matrix
  !   fockMatrix%values=fockMatrix%values-OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
  !   print *, "transformed fockMatrix"
  !   Call Matrix_show(fockMatrix)

    
  !   numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
  !   occupationNumber = MolecularSystem_getOcupationNumber( speciesID )
  !   overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
  !   !Asigns molecular orbital to each subsystem and calculates its energy contribution
  !   call Matrix_constructor (densityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
  !   call Vector_constructor (energyContribution, numberOfContractions, 0.0_8)
  !   allocate(holdValues(numberOfContractions))

  !   print *, "k, eigenvalue, energyContribution, subsystem"
  !   do k = 1 , numberOfContractions
  !      normCheck=0.0
  !      sumSubsystem(:)=0.0
  !      sumAB=0.0
  !      call Matrix_constructor (miniDensityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)       
  !      do mu = 1 , numberOfContractions
  !         do nu = 1 , numberOfContractions
  !            miniDensityMatrix%values(mu,nu)=miniDensityMatrix%values(mu,nu)+&
  !                 coefficients%values(mu,k)*coefficients%values(nu,k)
  !            if( OrbitalLocalizer_instance(speciesID)%subsystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subsystemList(nu) ) then
  !               sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))=&
  !                    sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))+&
  !                    overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
  !            else
  !               sumAB=sumAB+&
  !                    overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
  !            end if
  !         end do
  !      end do

  !      ! if(k .le. occupationNumber) densityMatrix%values=densityMatrix%values+miniDensityMatrix%values
       
  !      energyContribution%values(k)= sum(  transpose(miniDensityMatrix%values) &
  !           *  (  ( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
  !           + wavefunction_instance(speciesID)%twoParticlesMatrix%values &
  !           + wavefunction_instance(speciesID)%couplingMatrix%values &
  !           + WaveFunction_instance(speciesID)%exchangeCorrelationMatrix%values ))
       
  !      if(sumSubsystem(1)**2 .gt. sumSubsystem(2)**2) then
  !         OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 1
  !      else
  !         OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 2
  !      end if
  !      print *, k, eigenvalues%values(k), energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k)
  !   end do

  !   !Reorders coefficients matrix and eigenvalues according to the energy contribution
        
  !   do index1 = 1 , numberOfContractions-1
  !      holdIndex=index1
  !      holdEnergy=energyContribution%values(index1)
  !      holdValues=coefficients%values(:,index1)
  !      holdSystem=OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(index1) 

  !      do index2 = index1+1 , numberOfContractions
  !         if (energyContribution%values(index2).lt.energyContribution%values(holdIndex)) then
  !            print *, index1, index2, "switching", energyContribution%values(index2), energyContribution%values(holdIndex)
  !            holdIndex=index2
  !         end if
  !      end do
  !      coefficients%values(:,index1)=coefficients%values(:,holdIndex)
  !      coefficients%values(:,holdIndex)=holdValues
  !      energyContribution%values(index1)=energyContribution%values(holdIndex)
  !      energyContribution%values(holdIndex)=holdEnergy
  !      OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(index1)=OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(holdIndex)
  !      OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(holdIndex)=holdSystem

  !   end do
    
  !   print *, "sorted orbitals, subsystem"
  !   do k=1, numberOfContractions
  !      print *, k, energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k)
  !   end do
  !   eigenvalues%values=energyContribution%values

  !   deallocate(holdValues)
  !   call Matrix_destructor(densityMatrix)
  !   call Matrix_destructor(miniDensityMatrix)
  !   call Vector_destructor(energyContribution)
  !   !    end if
       
  !   !    do k = 1 , numberOfContractions
  !   !       print *, k, eigenvalues%values(k)
  !   ! end do
  !   ! print *, "subsystem 2 orbitals"
  !   ! do k = 1 , numberOfContractions
  !   !    if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k).eq.2) &
  !   !       print *, k, eigenvalues%values(k)
  !   ! end do

  !   !    i=i+1
  !   !    normCheck=0.0
  !   !    if ( abs(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
  !   !       do mu = 1 , numberOfContractions
  !   !          do nu = 1 , numberOfContractions
  !   !             normCheck=normCheck+WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(mu,i)*&
  !   !                  WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(nu,i)*&
  !   !                  WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
  !   !          end do
  !   !       end do
  !   !       ! print *, "eigenvalue", i, WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i), "normCheck", normCheck

  !   !       if ( normCheck .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD) then
  !   !          ! Shift orbital coefficients to the end of the matrix and Make energy a very large number
  !   !          do j = i , numberOfContractions-1
  !   !             WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(j)=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(j+1)
  !   !             WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,j) = WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,j+1)
  !   !          end do
  !   !          WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(numberOfContractions)=1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
  !   !          WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,numberOfContractions)=0.0
  !   !          i=i-1
  !   !       end if

  !   !    end if
  !   ! end do
    
  ! end subroutine OrbitalLocalizer_reorderSubsystemOrbitals
  


end module OrbitalLocalizer_
