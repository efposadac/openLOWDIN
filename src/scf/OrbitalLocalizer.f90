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
  use WaveFunction_
  use Convergence_

  implicit none


  type, public :: OrbitalLocalizer
     real(8) :: levelShiftingValue !control parameter
     real(8) :: populationThreshold !control parameter
     logical :: reduceSubsystemBasis !control parameter
     type(Matrix) :: projectionMatrix
     type(Matrix) :: hcoreMatrixA
     type(Matrix) :: fockMatrixA
     type(Matrix) :: waveFunctionCoefficientsA
     type(Vector) :: molecularOrbitalsEnergyA
     integer,allocatable :: subsystemList(:)
     integer,allocatable :: orbitalSubsystem(:)
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
       OrbitalLocalizer_instance(speciesID)%populationThreshold = CONTROL_INSTANCE%SUBSYSTEM_POPULATION_THRESHOLD
       OrbitalLocalizer_instance(speciesID)%reduceSubsystemBasis = CONTROL_INSTANCE%REDUCE_SUBSYSTEM_BASIS
       

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

    character(30) :: nameOfSpecies
    integer :: statusSystem

    integer :: numberOfContractions
    integer :: mu,nu, i

    
    !! Convert lowdin fchk files to erkale chk files
    nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    open(unit=30, file="erkale.read", status="replace", form="formatted")

    write(30,*) "LoadFChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".fchk"
    write(30,*) "SaveChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".chk"
    write(30,*) "Reorthonormalize true"
    
    close(30)

    statusSystem = system ("erkale_fchkpt erkale.read")

    !! Localize orbitals
    open(unit=30, file="erkale.local", status="replace", form="formatted")
    
    write(30,*) "LoadChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".chk"
    write(30,*) "SaveChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".local.chk"
    write(30,*) "Method MU"
    write(30,*) "Virtual false"
    write(30,*) "Maxiter 500"

    close(30)

    statusSystem = system ("erkale_loc_omp erkale.local")
    
    !!Convert erkale chk files to lowdin fchk files
    open(unit=30, file="erkale.write", status="replace", form="formatted")
    
    write(30,*) "LoadChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".local.chk"
    write(30,*) "SaveFChk ", trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".local.fchk"
    
    close(30)

    statusSystem = system ("erkale_fchkpt erkale.write")
    
    !! Read orbital coefficients from fchk files
    call MolecularSystem_readFchk(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".local.fchk",  orbitalCoefficients, densityMatrix, nameOfSpecies )

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
    real(8),allocatable :: excCorrEnergyA(:,:)
    real(8),allocatable :: excCorrEnergyB(:,:)
    real(8),allocatable :: excCorrEnergyAB(:,:)
    real(8),allocatable :: particlesInGridA(:)
    real(8),allocatable :: particlesInGridB(:)
    real(8),allocatable :: particlesInGridAB(:)
    real(8),allocatable :: densStd(:)
    real(8),allocatable :: auxEnergy(:)
    
    real(8) :: totalDensStd
    logical :: SUBSYSTEM_SCF_CONTINUE
    
    character(30) :: nameOfSpecies
    
    integer :: numberOfSpecies,numberOfContractions, numberOfCenters, occupationNumber
    integer :: subsystemMu, subsystemNu, orbitalsInA, orbitalsInB
    integer :: mu,nu,alpha,beta, i, ii, j, jj, k, l, iter

    character(50) :: labels(2)
    character(100) :: densFileA,densFileB,densFileAB, wfnFile
    integer :: densUnitA,densUnitB,densUnitAB, wfnUnit

    real(8) :: totalKineticEnergyA,totalKineticEnergyB, puntualInteractionEnergy
    real(8) :: totalQuantumPuntualInteractionEnergyA, totalQuantumPuntualInteractionEnergyB
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

       nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
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
          if(sumSubsystem(1) .gt. OrbitalLocalizer_instance(speciesID)%populationThreshold) then
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
       densityMatrixAB(speciesID)%values=densityMatrixA(speciesID)%values+densityMatrixB(speciesID)%values
       ! print *, "densityMatrixA"
       ! Call Matrix_show(densityMatrixA(speciesID))
       ! print *, "densityMatrixB"
       ! Call Matrix_show(densityMatrixB(speciesID))

       !Calculates atomic Mulliken population on each fragment
       populationMatrixA(speciesID)%values= matmul(densityMatrixA(speciesID)%values, overlapMatrix%values )
       populationMatrixB(speciesID)%values= matmul(densityMatrixB(speciesID)%values, overlapMatrix%values )
       print *, ""
       print *, "Subsystem Atomic Mulliken Population for: ", nameOfSpecies
       print *, ""
       write(*,"(A20,A10,A10)") "Atom:      ", "A   ", "B   "
       mu=0
       do i=1, numberOfCenters
          atomPopulationA(speciesID)%values(i)=0
          atomPopulationB(speciesID)%values(i)=0
          do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             do k=1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
                mu=mu+1
                atomPopulationA(speciesID)%values(i)=atomPopulationA(speciesID)%values(i)+populationMatrixA(speciesID)%values(mu,mu)
                atomPopulationB(speciesID)%values(i)=atomPopulationB(speciesID)%values(i)+populationMatrixB(speciesID)%values(mu,mu)
             end do
          end do
          if (abs(atomPopulationA(speciesID)%values(i)) .lt. OrbitalLocalizer_instance(speciesID)%populationThreshold ) then
             write(*,"(A20,F10.6,F10.6)") trim(ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner))//"*",&
                  atomPopulationA(speciesID)%values(i), atomPopulationB(speciesID)%values(i)
          else
             write(*,"(A20,F10.6,F10.6)") trim(ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner)),&
                  atomPopulationA(speciesID)%values(i), atomPopulationB(speciesID)%values(i)
          end if
       end do
       write(*,"(T20,A20)") "____________________"
       write(*,"(A20,F10.6,F10.6)") "Total", sum(atomPopulationA(speciesID)%values(:)), sum(atomPopulationB(speciesID)%values(:))
       print *, "* will be projected out from A basis set"
       
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
       
       ! print *, "projectionMatrix"
       ! Call Matrix_show(OrbitalLocalizer_instance(speciesID)%projectionMatrix)

       !Adds diagonal proyection elements to atoms with small contributions to A orbitals - small mulliken population
       !These will correspond to virtual orbitals with high energy that can be removed from post-SCF calculations
       if(OrbitalLocalizer_instance(speciesID)%reduceSubsystemBasis) then
          mu=0
          OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=0
          do i=1, numberOfCenters
             do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
                do k=1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
                   mu=mu+1
                   if (abs(atomPopulationA(speciesID)%values(i)) .lt. OrbitalLocalizer_instance(speciesID)%populationThreshold ) then
                      OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB=&
                           OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB+1
                      OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)=&
                           OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(mu,mu)+1.0
                   end if
                end do
             end do
          end do
       end if
       OrbitalLocalizer_instance(speciesID)%removedOrbitalsA=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB&
            +OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB

       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA=numberOfContractions&
            -OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB

       if(OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA .eq.0 ) OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA = 0
       if(OrbitalLocalizer_instance(speciesID)%removedOrbitalsA .gt. numberOfContractions ) OrbitalLocalizer_instance(speciesID)%removedOrbitalsA=numberOfContractions
       
       ! print *, "projectionMatrix-2"
       ! Call Matrix_show(OrbitalLocalizer_instance(speciesID)%projectionMatrix)

    end do
    
  !!Modify and write molecular system for subsystem A
    print *, ""
    print *, "**************************************************************************"
    print *, "*************Molecular subsystem A****************************************"
    print *, "**************************************************************************"

    !For now, we only remove the particles from system B
    do speciesID=1, numberOfSpecies
       MolecularSystem_instance%species(speciesID)%ocupationNumber=&
            MolecularSystem_instance%species(speciesID)%ocupationNumber&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB
       MolecularSystem_instance%species(speciesID)%internalSize=&
            MolecularSystem_instance%species(speciesID)%internalSize&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB&
            *MolecularSystem_getEta(speciesID)
    end do
    !To avoid trouble caused by removing and adding particles
    ! CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE=.true.
    call MolecularSystem_saveToFile( "lowdin-subsystemA" )

    call MolecularSystem_destroy()

    !!Load CONTROL Parameters
    call MolecularSystem_loadFromFile( "LOWDIN.DAT", "lowdin-subsystemA" )

    !!Load the system in lowdin.sys format
    call MolecularSystem_loadFromFile( "LOWDIN.SYS", "lowdin-subsystemA" )

    ! call MolecularSystem_showInformation()  
    ! call MolecularSystem_showParticlesInformation()

    print *, ""
    print *, "Orbital Subsystem Distribution"
    print *, ""
    write(*,"(A15,A15,A15,A15,A15)") "Species", "Occupied A","Virtual A","Occupied B", "Virtual B"
    write(*,"(A75)") "---------------------------------------------------------------------------"
    do speciesID=1, numberOfSpecies
       write(*,"(A15,I15,I15,I15,I15)") trim(MolecularSystem_getNameOfSpecie(speciesID)), &
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsA,&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsA,&
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsB,&
            OrbitalLocalizer_instance(speciesID)%virtualOrbitalsB
    end do

    
    print *, ""
    print *, "Subsystem A atoms"
    print *, ""
    call MolecularSystem_showCartesianMatrix(1)

    !!Save density matrix B and run DFT for the subsystem B
    
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       call WaveFunction_writeDensityMatricesToFile(trim(densFileB),densityMatrixB(:))
       call system("lowdin-DFT.x SCF_DFT "//trim(densFileB))
    end if
    !!!Build subsystem B matrices - these do not change in the second SCF cycle
    !Two particles and coupling matrices
    do speciesID=1, numberOfSpecies     
       
       !Only Coulomb (factor=0.0)
       call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
         densityMatrixIN=densityMatrixB(speciesID),&
         factorIN=0.0_8,&
         twoParticlesMatrixOUT=hartreeMatrixB(speciesID,speciesID) )
      
       !Exchange-HF matrix with the HF fraction used in the global SCF calculation       
       call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
         densityMatrixIN=densityMatrixB(speciesID),&
         twoParticlesMatrixOUT=exchangeHFMatrixB(speciesID))

       exchangeHFMatrixB(speciesID)%values=&
            exchangeHFMatrixB(speciesID)%values-hartreeMatrixB(speciesID,speciesID)%values

       do otherSpeciesID=1, numberOfSpecies
          call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       end do
       
       call WaveFunction_buildCouplingMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
         densityMatricesIN=densityMatrixB(1:numberOfSpecies),&
         couplingMatrixOUT=couplingMatrixB(speciesID),&
         hartreeMatricesOUT=auxMatrix(1:numberOfSpecies)) 

       do otherSpeciesID=1, numberOfSpecies
          if(otherSpeciesID.ne.speciesID) hartreeMatrixB(speciesID,otherSpeciesID)%values=auxMatrix(otherSpeciesID)%values
          call Matrix_destructor (auxMatrix(otherSpeciesID))
       end do

       auxEnergy(:)=0.0
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               excFileIN=trim(densFileB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridB(speciesID) )
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
          write(*,"(A10,A20,A20,A20,A20)") "Iteration", "EnergyAB","EnergyChange","DensityChangeA", "ParticlesInGridAB"
       else
          write(*,"(A10,A20,A20,A20)") "Iteration", "EnergyAB","EnergyChange","DensityChangeA"
       end if
       write(*,*) "---------------------------------------------------------"
    end if


    newTotalEnergy=0.0
    do while( SUBSYSTEM_SCF_CONTINUE )
       iter=iter+1
    
       !Run DFT with the total density and the subsystem density 
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then

          !!Save density matrix A, A+B
          call WaveFunction_writeDensityMatricesToFile(trim(densFileA),densityMatrixA(:))
          call WaveFunction_writeDensityMatricesToFile(trim(densFileAB),densityMatrixAB(:))
          
          !!Run DFT for subsystem A
          call system("lowdin-DFT.x SCF_DFT "//trim(densFileA))
          !!Run DFT for complete system A+B
          call system("lowdin-DFT.x SCF_DFT "//trim(densFileAB))
          
       end if
       
       !Calculates the fock matrix with the new subsystem density 
       do speciesID=1, numberOfSpecies
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
          
          !Updates two particles matrix - only Coulomb (factor=0.0)
          call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               densityMatrixIN=densityMatrixA(speciesID),&
               factorIN=0.0_8,&
               twoParticlesMatrixOUT=hartreeMatrixA(speciesID,speciesID))

          !Obtains exchange-correlation matrix from a full HF calculation (factor=1.0)
          call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               densityMatrixIN=densityMatrixA(speciesID),&
               factorIN=1.0_8,&
               twoParticlesMatrixOUT=exchangeHFMatrixA(speciesID))

          exchangeHFMatrixA(speciesID)%values=&
               exchangeHFMatrixA(speciesID)%values-hartreeMatrixA(speciesID,speciesID)%values

          !Updates coupling matrices
          do otherSpeciesID=1, numberOfSpecies
             call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
          end do

          call WaveFunction_buildCouplingMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
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
             call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
                  excFileIN=trim(densFileA)//".exc",&
                  exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixA(speciesID),&
                  exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
                  particlesInGridOUT=particlesInGridA(speciesID) )
             do otherSpeciesID=1, numberOfSpecies
                excCorrEnergyA(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
             end do

             auxEnergy=0.0
             call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
                  excFileIN=trim(densFileAB)//".exc",&
                  exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixAB(speciesID),&
                  exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
                  particlesInGridOUT=particlesInGridAB(speciesID) ) 
             do otherSpeciesID=1, numberOfSpecies
                excCorrEnergyAB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
             end do
                  
          end if
          
                  
          ! Updates Fock Matrix
          OrbitalLocalizer_instance(speciesID)%fockMatrixA%values = &
               wavefunction_instance(speciesID)%hcoreMatrix%values &
               +hartreeMatrixA(speciesID,speciesID)%values&
               +couplingMatrixA(speciesID)%values&
               +exchangeHFMatrixA(speciesID)%values&
               +hartreeMatrixB(speciesID,speciesID)%values&
               +exchangeHFMatrixB(speciesID)%values&
               +couplingMatrixB(speciesID)%values&
               +(exchangeCorrelationMatrixAB(speciesID)%values&
               -exchangeCorrelationMatrixA(speciesID)%values)&
               +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values

          call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
               OrbitalLocalizer_instance(speciesID)%fockMatrixA, densityMatrixA(speciesID), &
               WaveFunction_instance(speciesID)%OverlapMatrix, &
               methodType=SCF_CONVERGENCE_DAMPING, &
               coefficientMatrix=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, speciesID=speciesID )

          call Convergence_run( WaveFunction_instance(speciesID)%convergenceMethod )


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
          newTotalEnergy= newTotalEnergy+&
               sum(transpose(densityMatrixAB(speciesID)%values)* WaveFunction_instance( speciesID )%hCoreMatrix%values) &
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* hartreeMatrixA(speciesID,speciesID)%values)&
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)&
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* couplingMatrixA(speciesID)%values)&
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* hartreeMatrixB(speciesID,speciesID)%values)&
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)&
               +0.5*sum(transpose(densityMatrixAB(speciesID)%values)* couplingMatrixB(speciesID)%values)&
               -(sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixAB(speciesID)%values)&
               -sum(transpose(densityMatrixA(speciesID)%values)*exchangeCorrelationMatrixA(speciesID)%values))&
               +excCorrEnergyAB(speciesID,speciesID)-excCorrEnergyA(speciesID,speciesID) +&
               OrbitalLocalizer_instance(speciesID)%levelShiftingValue*sum(transpose(densityMatrixA(speciesID)%values)*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)
          do otherSpeciesID=speciesID+1, numberOfSpecies
             newTotalEnergy= newTotalEnergy+excCorrEnergyAB(speciesID,otherSpeciesID)
          end do
       end do
       ! print *, "total energy", newTotalEnergy

       !Calculates the subsystem orbitals with the new fock matrix
       do speciesID=1, numberOfSpecies
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

          call Matrix_copyConstructor( fockMatrixTransformed, OrbitalLocalizer_instance(speciesID)%fockMatrixA )

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
          densityMatrixAB(speciesID)%values=densityMatrixA(speciesID)%values+&
               densityMatrixB(speciesID)%values

       end do

       !Check convergence and print messages
       totalDensStd=sqrt(sum(densStd(:)**2))
       deltaEnergy=newTotalEnergy-totalEnergy
       !Writes iteration results
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          write (*,"(I10,F20.10,F20.10,F20.10,F20.10)") iter, newTotalEnergy, deltaEnergy , totalDensStd, sum(particlesInGridAB(:)) 
       else
          write (*,"(I10,F20.10,F20.10,F20.10)") iter, newTotalEnergy, deltaEnergy, totalDensStd 
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

    do speciesID=1, numberOfSpecies
       write(*,*) ""
       write(*,*) " Subsystem Eigenvectors for: ", trim(MolecularSystem_instance%species(speciesID)%name )
       write(*,*) "-----------------------------"
       write(*,*) ""

       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       call Matrix_constructor(auxMatrix(speciesID),int(numberOfContractions,8),&
            int(numberOfContractions-OrbitalLocalizer_instance(speciesID)%removedOrbitalsA,8),0.0_8)

       do i=1, numberOfContractions
          do j=1, numberOfContractions-OrbitalLocalizer_instance(speciesID)%removedOrbitalsA
             auxMatrix(speciesID)%values(i,j)=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA%values(i,j)
          end do
       end do

       call Matrix_show( auxMatrix(speciesID), &
            rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
            columnkeys = string_convertvectorofrealstostring( OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA ),&
            flags=WITH_BOTH_KEYS)

    end do

  ! Final energy evaluation - larger integration grid for DFT 

    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       !!Save density matrix A, B, A+B
       call WaveFunction_writeDensityMatricesToFile(trim(densFileA),densityMatrixA(:))
       call WaveFunction_writeDensityMatricesToFile(trim(densFileB),densityMatrixB(:))
       call WaveFunction_writeDensityMatricesToFile(trim(densFileAB),densityMatrixAB(:))

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
       call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
            densityMatrixIN=densityMatrixA(speciesID),&
            factorIN=0.0_8,&
            twoParticlesMatrixOUT=hartreeMatrixA(speciesID,speciesID))

       !Obtains exchange-correlation matrix from a full HF calculation (factor=1.0)
       call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
            densityMatrixIN=densityMatrixA(speciesID),&
            factorIN=1.0_8,&
            twoParticlesMatrixOUT=exchangeHFMatrixA(speciesID))

       exchangeHFMatrixA(speciesID)%values=&
            exchangeHFMatrixA(speciesID)%values-hartreeMatrixA(speciesID,speciesID)%values

       !Updates coupling matrices
       do otherSpeciesID=1, numberOfSpecies
          call Matrix_constructor (auxMatrix(otherSpeciesID), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       end do

       call WaveFunction_buildCouplingMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
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
          call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               excFileIN=trim(densFileA)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixA(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridA(speciesID) )
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyA(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

          auxEnergy=0.0
          call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               excFileIN=trim(densFileB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridB(speciesID) )
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

          auxEnergy=0.0
          call WaveFunction_readExchangeCorrelationMatrix(MolecularSystem_getNameOfSpecie(speciesID),&
               excFileIN=trim(densFileAB)//".exc",&
               exchangeCorrelationMatrixOUT=exchangeCorrelationMatrixAB(speciesID),&
               exchangeCorrelationEnergyOUT=auxEnergy(1:numberOfSpecies),&
               particlesInGridOUT=particlesInGridAB(speciesID) ) 
          do otherSpeciesID=1, numberOfSpecies
             excCorrEnergyAB(speciesID,otherSpeciesID)=auxEnergy(otherSpeciesID)
          end do

       end if

       !!Updates Fock Matrix
       OrbitalLocalizer_instance(speciesID)%fockMatrixA%values = &
            wavefunction_instance(speciesID)%hcoreMatrix%values &
            +hartreeMatrixA(speciesID,speciesID)%values&
            +couplingMatrixA(speciesID)%values&
            +exchangeHFMatrixA(speciesID)%values&
            +hartreeMatrixB(speciesID,speciesID)%values&
            +exchangeHFMatrixB(speciesID)%values&
            +couplingMatrixB(speciesID)%values&
            +(exchangeCorrelationMatrixAB(speciesID)%values&
            -exchangeCorrelationMatrixA(speciesID)%values)&
            +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values

    end do
       
  !! Obtain energy compotents for whole system
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
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values), &
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values), &
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* hartreeMatrixB(speciesID,otherSpeciesID)%values)+ &
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* hartreeMatrixA(speciesID,otherSpeciesID)%values)

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
             0.5*sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)+ &             
             0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
           
           totalExchangeHFEnergyA=totalExchangeHFEnergyA+0.5*sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
           totalExchangeHFEnergyB=totalExchangeHFEnergyB+0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)
           totalExchangeHFEnergyAB=totalExchangeHFEnergyAB+0.5*sum(transpose(densityMatrixA(speciesID)%values)* exchangeHFMatrixB(speciesID)%values)+ &
                0.5*sum(transpose(densityMatrixB(speciesID)%values)* exchangeHFMatrixA(speciesID)%values)
        end do

     write (6,"(T10,A90)") "_____________________________________________________________"
     write (6,"(A35,F20.12,F20.12,F20.12)") "Total Exchange energy = ", totalExchangeHFEnergyA, totalExchangeHFEnergyB, totalExchangeHFEnergyAB

     write(*,*) ""

     if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
        write(*,*) ""
        write(*,"(A35,A20,A20,A20)") " Exchange-Correlation(DFT) energy: ", "(Subsystem A)*", "Subsystem B", "A-B Interaction"
        write(*,*) "---------------------------------------------------------------------------------------------"
        write(*,*) ""
        totalExchangeCorrelationEnergyA=0.0
        totalExchangeCorrelationEnergyB=0.0
        totalExchangeCorrelationEnergyAB=0.0
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
        write(*,"(A35,A20,A20,A20)") " A Embedding in B Potential Energy: ", "A-B Interaction"
        write(*,*) "-----------------------------------------------------------"
        write(*,*) ""
        totalEmbeddingPotentialEnergyA=0.0
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


     ! if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then

     !    write(*,*) ""
     !    write(*,*) " External Potential energy: "
     !    write(*,*) "----------------"
     !    write(*,*) ""

     !    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
     !       write (6,"(T10,A26,A2,F20.12)") trim( MolecularSystem_instance%species(speciesID)%name) // &
     !            " Ext Pot energy  ","= ", WaveFunction_instance(speciesID)%externalPotentialEnergy
     !    end do
     !    totalExternalPotentialEnergy=sum(WaveFunction_instance(:)%externalPotentialEnergy)
     !    write (6,"(T10,A28,F20.12)") "External Potential energy  = ", totalExternalPotentialEnergy             
        
     ! end if

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
          totalHartreeEnergyA+totalExchangeHFEnergyA
     write(*,"(A35,F20.12)") "TOTAL SUBSYSTEM B ENERGY = ", totalKineticEnergyB+totalQuantumPuntualInteractionEnergyB+&
          totalHartreeEnergyB+totalExchangeHFEnergyB+totalExchangeCorrelationEnergyB
     write(*,"(A35,F20.12)") "TOTAL A-B INTERACTION ENERGY = ", totalHartreeEnergyAB+totalExchangeHFEnergyAB+totalExchangeCorrelationEnergyAB+&
          totalEmbeddingPotentialEnergyA+totalProjectionCorrectionA

     totalKineticEnergy=totalKineticEnergyA+totalKineticEnergyB
     totalPotentialEnergy=puntualInteractionEnergy+&
          totalQuantumPuntualInteractionEnergyA+&
          totalQuantumPuntualInteractionEnergyB+&
          totalHartreeEnergyA+&
          totalHartreeEnergyB+&
          totalHartreeEnergyAB+&
          totalExchangeHFEnergyA+&
          totalExchangeHFEnergyB+&
          totalExchangeHFEnergyAB+&
          totalExchangeCorrelationEnergyB+&
          totalExchangeCorrelationEnergyAB+&          
          totalEmbeddingPotentialEnergyA+&
          totalProjectionCorrectionA
     write(*,"(T10,A50)") "_____________________"
     write(*,"(A35,F20.12)") "TOTAL ENERGY = ", totalKineticEnergy+totalPotentialEnergy
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


      
    !!**********************************************************
    !! Save matrices to subsystem lowdin.wfn file
    !!
    wfnUnit = 300
    wfnFile = "lowdin-subsystemA.wfn"

    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    rewind(wfnUnit)

    labels = ""

    do speciesID = 1, numberOfSpecies

       labels(2) = MolecularSystem_getNameOfSpecie(speciesID)

       labels(1) = "REMOVED-ORBITALS"

       call Vector_writeToFile(unit=wfnUnit, binary=.true., value=real(OrbitalLocalizer_instance(speciesID)%removedOrbitalsA,8), arguments= labels )
       ! labels(1) = "TWOPARTICLES"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%twoParticlesMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       ! labels(1) = "COUPLING"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%couplingMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       ! labels(1) = "EXCHANGE-CORRELATION"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%exchangeCorrelationMatrix, unit=wfnUnit, binary=.true., arguments = labels )  

       ! labels(1) = "EXCHANGE-CORRELATION-ENERGY"
       ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=WaveFunction_instance(speciesID)%exchangeCorrelationEnergy, arguments= labels )

       labels(1) = "COEFFICIENTS"
       call Matrix_writeToFile(OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "DENSITY"
       call Matrix_writeToFile(densityMatrixA(speciesID), unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "HCORE"
       OrbitalLocalizer_instance(speciesID)%hcoreMatrixA%values=WaveFunction_instance(speciesID)%hcoreMatrix%values&
            +hartreeMatrixB(speciesID,speciesID)%values&
            +exchangeHFMatrixB(speciesID)%values&
            +couplingMatrixB(speciesID)%values&
            +(exchangeCorrelationMatrixAB(speciesID)%values-exchangeCorrelationMatrixA(speciesID)%values)&
            +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values

       ! print *, "mod hcore matrix ", speciesID
       ! call Matrix_show(OrbitalLocalizer_instance(speciesID)%hcoreMatrixA)
       
       call Matrix_writeToFile(OrbitalLocalizer_instance(speciesID)%hcoreMatrixA, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "ORBITALS"
       call Vector_writeToFile(OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA, unit=wfnUnit, binary=.true., arguments = labels )

       ! labels(1) = "FOCK"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%fockMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       ! if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
       !    labels(1) = "EXTERNAL_POTENTIAL"
       !    call Matrix_writeToFile(WaveFunction_instance(speciesID)%externalPotentialMatrix, unit=wfnUnit, binary=.true., arguments = labels )
       ! end if

       ! if (CONTROL_instance%COSMO) then
       !    labels(1) = "COSMO2"
       !    call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmo2, unit=wfnUnit, binary=.true., arguments = labels )  
       !    labels(1) = "COSMOCOUPLING"
       !    call Matrix_writeToFile(WaveFunction_instance(speciesID)%cosmoCoupling, unit=wfnUnit, binary=.true., arguments = labels ) 
       ! end if

    end do

    !!**********************************************************
    !! Save Some energies
    !!
    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=newTotalEnergy, arguments=["TOTALENERGY"])

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
    call system ("lowdin-CalcProp.x lowdin-subsystemA")

    
  end subroutine OrbitalLocalizer_levelShiftSubsystemOrbitals

  subroutine OrbitalLocalizer_reorderSubsystemOrbitals(speciesID,coefficients,fockMatrix,eigenvalues)
    integer :: speciesID
    type(Matrix) :: coefficients
    type(Matrix) :: fockMatrix
    type(Vector) :: eigenvalues
    
    type(Matrix) :: overlapMatrix
    type(Matrix) :: densityMatrix
    type(Matrix) :: miniDensityMatrix
    type(Vector) :: energyContribution
    real(8) :: sumSubsystem(2), sumAB, normCheck, holdEnergy
    real(8),allocatable ::holdValues(:)

    integer :: numberOfContractions, occupationNumber
    integer :: index1,index2,holdIndex,holdSystem
    integer :: mu,nu,alpha,beta, i, ii, j, jj, k

    !Removes the projection matrix from the fock matrix
    fockMatrix%values=fockMatrix%values-OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
    print *, "transformed fockMatrix"
    Call Matrix_show(fockMatrix)

    
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    occupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
    !Asigns molecular orbital to each subsystem and calculates its energy contribution
    call Matrix_constructor (densityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
    call Vector_constructor (energyContribution, numberOfContractions, 0.0_8)
    allocate(holdValues(numberOfContractions))

    print *, "k, eigenvalue, energyContribution, subsystem"
    do k = 1 , numberOfContractions
       normCheck=0.0
       sumSubsystem(:)=0.0
       sumAB=0.0
       call Matrix_constructor (miniDensityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)       
       do mu = 1 , numberOfContractions
          do nu = 1 , numberOfContractions
             miniDensityMatrix%values(mu,nu)=miniDensityMatrix%values(mu,nu)+&
                  coefficients%values(mu,k)*coefficients%values(nu,k)
             if( OrbitalLocalizer_instance(speciesID)%subsystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subsystemList(nu) ) then
                sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))=&
                     sumSubsystem(OrbitalLocalizer_instance(speciesID)%subsystemList(mu))+&
                     overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
             else
                sumAB=sumAB+&
                     overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
             end if
          end do
       end do

       ! if(k .le. occupationNumber) densityMatrix%values=densityMatrix%values+miniDensityMatrix%values
       
       energyContribution%values(k)= sum(  transpose(miniDensityMatrix%values) &
            *  (  ( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
            + wavefunction_instance(speciesID)%twoParticlesMatrix%values &
            + wavefunction_instance(speciesID)%couplingMatrix%values &
            + WaveFunction_instance(speciesID)%exchangeCorrelationMatrix%values ))
       
       if(sumSubsystem(1)**2 .gt. sumSubsystem(2)**2) then
          OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 1
       else
          OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k) = 2
       end if
       print *, k, eigenvalues%values(k), energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k)
    end do

    !Reorders coefficients matrix and eigenvalues according to the energy contribution
        
    do index1 = 1 , numberOfContractions-1
       holdIndex=index1
       holdEnergy=energyContribution%values(index1)
       holdValues=coefficients%values(:,index1)
       holdSystem=OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(index1) 

       do index2 = index1+1 , numberOfContractions
          if (energyContribution%values(index2).lt.energyContribution%values(holdIndex)) then
             print *, index1, index2, "switching", energyContribution%values(index2), energyContribution%values(holdIndex)
             holdIndex=index2
          end if
       end do
       coefficients%values(:,index1)=coefficients%values(:,holdIndex)
       coefficients%values(:,holdIndex)=holdValues
       energyContribution%values(index1)=energyContribution%values(holdIndex)
       energyContribution%values(holdIndex)=holdEnergy
       OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(index1)=OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(holdIndex)
       OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(holdIndex)=holdSystem

    end do
    
    print *, "sorted orbitals, subsystem"
    do k=1, numberOfContractions
       print *, k, energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k)
    end do
    eigenvalues%values=energyContribution%values

    deallocate(holdValues)
    call Matrix_destructor(densityMatrix)
    call Matrix_destructor(miniDensityMatrix)
    call Vector_destructor(energyContribution)
    !    end if
       
    !    do k = 1 , numberOfContractions
    !       print *, k, eigenvalues%values(k)
    ! end do
    ! print *, "subsystem 2 orbitals"
    ! do k = 1 , numberOfContractions
    !    if(OrbitalLocalizer_instance(speciesID)%orbitalSubsystem(k).eq.2) &
    !       print *, k, eigenvalues%values(k)
    ! end do

    !    i=i+1
    !    normCheck=0.0
    !    if ( abs(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
    !       do mu = 1 , numberOfContractions
    !          do nu = 1 , numberOfContractions
    !             normCheck=normCheck+WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(mu,i)*&
    !                  WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(nu,i)*&
    !                  WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
    !          end do
    !       end do
    !       ! print *, "eigenvalue", i, WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(i), "normCheck", normCheck

    !       if ( normCheck .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD) then
    !          ! Shift orbital coefficients to the end of the matrix and Make energy a very large number
    !          do j = i , numberOfContractions-1
    !             WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(j)=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(j+1)
    !             WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,j) = WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,j+1)
    !          end do
    !          WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(numberOfContractions)=1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
    !          WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,numberOfContractions)=0.0
    !          i=i-1
    !       end if

    !    end if
    ! end do
    
  end subroutine OrbitalLocalizer_reorderSubsystemOrbitals
  


end module OrbitalLocalizer_
