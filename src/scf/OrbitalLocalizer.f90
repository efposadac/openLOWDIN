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
     real(8) :: levelShiftingValue !should be a value in Control
     type(Matrix) :: projectionMatrix
     type(Matrix) :: fockMatrixA
     type(Matrix) :: waveFunctionCoefficientsA
     type(Vector) :: molecularOrbitalsEnergyA
     type(Matrix) :: densityMatrixA
     type(Matrix) :: densityMatrixB
     integer,allocatable :: subSystemList(:)
     integer,allocatable :: orbitalSubSystem(:)
     integer :: occupiedOrbitalsInA
     integer :: virtualOrbitalsInA
     integer :: occupiedOrbitalsInB
     integer :: virtualOrbitalsInB
     integer :: occupiedDelocalizedOrbitals
     integer :: virtualDelocalizedOrbitals

  end type OrbitalLocalizer

  type(OrbitalLocalizer), public, allocatable :: OrbitalLocalizer_instance(:)

contains
 
  subroutine OrbitalLocalizer_constructor()
    integer :: speciesID, numberOfSpecies
    integer :: numberOfContractions, ocupationNumber
    integer :: mu, i, j, k

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies
    allocate(OrbitalLocalizer_instance(numberOfSpecies))

    do speciesID=1, numberOfSpecies
       OrbitalLocalizer_instance(speciesID)%levelShiftingValue = 1.0E6
       ! OrbitalLocalizer_instance(speciesID)%levelShiftingValue = 0.0
       !Gets atomic subsystem classification
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )

       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor(OrbitalLocalizer_instance(speciesID)%fockMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Vector_constructor(OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA, numberOfContractions, 0.0_8)

       
       allocate(OrbitalLocalizer_instance(speciesID)%subSystemList(numberOfContractions))
       allocate(OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(numberOfContractions))

       print *, "subSystemList"
       mu=0
       do i = 1, size(MolecularSystem_instance%species(speciesID)%particles)
          do j = 1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
             do k = 1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
                mu=mu+1
                OrbitalLocalizer_instance(speciesID)%subSystemList(mu)=MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%subSystem
                print *, mu, OrbitalLocalizer_instance(speciesID)%subSystemList(mu)
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
    call MolecularSystem_readFchk( orbitalCoefficients, densityMatrix, nameOfSpecies )

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

  
  subroutine OrbitalLocalizer_levelShiftSubSystemOrbitals()
    integer :: speciesID
    type(Matrix) :: coefficients
    type(Matrix) :: fockMatrix
    type(Matrix) :: overlapMatrix
    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: newDensityMatrix
    type(Matrix) :: twoParticlesMatrixA, twoParticlesMatrixB
    type(Matrix) :: exchangeCorrelationMatrixA, exchangeCorrelationMatrixB
    real(8) :: sumSubSystem(2), sumAB, normCheck
    real(8) :: excCorrEnergyA,excCorrEnergyB
    real(8) :: totalEnergy
    
    real(8) :: densStd
    logical :: SUBSYSTEM_SCF_CONTINUE
    
    character(30) :: nameOfSpecies
    
    integer :: numberOfSpecies,numberOfContractions, ocupationNumber
    integer :: subSystemMu, subSystemNu, orbitalsInA, orbitalsInB
    integer :: mu,nu,alpha,beta, i, ii, j, jj, k, iter

    character(50) :: labels(2), excFile, densFile, wfnFile
    integer :: excUnit, densUnit, wfnUnit

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    
    do speciesID=1, numberOfSpecies

       nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
       ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
       overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
       coefficients=WaveFunction_instance( speciesID )%waveFunctionCoefficients
       fockMatrix=WaveFunction_instance( speciesID )%fockMatrix
       call Matrix_constructor (newDensityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (twoParticlesMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (twoParticlesMatrixB, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeCorrelationMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (exchangeCorrelationMatrixB, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (OrbitalLocalizer_instance(speciesID)%densityMatrixA, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (OrbitalLocalizer_instance(speciesID)%densityMatrixB, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
       call Matrix_constructor (OrbitalLocalizer_instance(speciesID)%projectionMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

       !Asigns molecular orbital to each subsystem
       OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInA=0
       OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB=0
       ! OrbitalLocalizer_instance(speciesID)%occupiedDelocalizedOrbitals=0
       ! OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInA=0
       ! OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInB=0
       ! OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals=0

       print *, "orbital, assigned to, sumSubSystem(1), sumSubSystem(1), sumAB, normCheck"
       print *, "occupied"
       !!Occupied
       do k = 1 , ocupationNumber
          normCheck=0.0
          sumSubSystem(:)=0.0
          sumAB=0.0
          do mu = 1 , numberOfContractions
             do nu = 1 , numberOfContractions
                if( OrbitalLocalizer_instance(speciesID)%subSystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subSystemList(nu) ) then
                   sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))=&
                        sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))+&
                        overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
                else
                   sumAB=sumAB+&
                        overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
                end if
                normCheck=normCheck+coefficients%values(mu,k)*&
                     coefficients%values(nu,k)*&
                     WaveFunction_instance(speciesID)%overlapMatrix%values(mu,nu)
             end do
          end do
          ! if(sumAB**2 .gt. (sumSubSystem(1)**2+sumSubSystem(2)**2) ) then
          !    OrbitalLocalizer_instance(speciesID)%occupiedDelocalizedOrbitals=OrbitalLocalizer_instance(speciesID)%occupiedDelocalizedOrbitals+1
          !    OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 0
          if(sumSubSystem(1) .gt. 0.1) then
             OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInA=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInA+1
             OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 1
          else
             OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB=OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB+1
             OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 2
          end if
          print *, k, OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k), sumSubSystem(1), sumSubSystem(2), sumAB, normCheck
       end do
       print *, "occupiedOrbitalsInA, occupiedOrbitalsInB", &
            OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInA, OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB

       !!Virtual
       ! do k = ocupationNumber+1 , numberOfContractions
       !    normCheck=0.0
       !    sumSubSystem(:)=0.0
       !    sumAB=0.0
       !    do mu = 1 , numberOfContractions
       !       do nu = 1 , numberOfContractions
       !          if( OrbitalLocalizer_instance(speciesID)%subSystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subSystemList(nu) ) then
       !             sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))=&
       !                  sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))+&
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
       !    if(sumAB**2 .gt. (sumSubSystem(1)**2+sumSubSystem(2)**2) ) then
       !       OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals=OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 0
       !    else if(sumSubSystem(1)**2 .gt. 0.1) then
       !       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInA=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInA+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 1
       !    else
       !       OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInB=OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInB+1
       !       OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 2
       !    end if
       !    print *, k, OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k), sumSubSystem(1), sumSubSystem(2), sumAB, normCheck
       ! end do
       ! print *, "virtualOrbitalsInA, virtualOrbitalsInB, virtualDelocalizedOrbitals", &
       !      OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInA, OrbitalLocalizer_instance(speciesID)%virtualOrbitalsInB, OrbitalLocalizer_instance(speciesID)%virtualDelocalizedOrbitals


       !Decomposes atomic density matrix
       do k = 1 , ocupationNumber
          if(OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k).eq.1) then
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   OrbitalLocalizer_instance(speciesID)%densityMatrixA%values(mu,nu)=&
                        OrbitalLocalizer_instance(speciesID)%densityMatrixA%values(mu,nu)&
                        +MolecularSystem_getEta( speciesID )*coefficients%values(mu,k)*coefficients%values(nu,k)
                end do
             end do
          else
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   OrbitalLocalizer_instance(speciesID)%densityMatrixB%values(mu,nu)=&
                        OrbitalLocalizer_instance(speciesID)%densityMatrixB%values(mu,nu)&
                        +MolecularSystem_getEta( speciesID )*coefficients%values(mu,k)*coefficients%values(nu,k)
                end do
             end do
          end if
       end do


       print *, "densityMatrixA"
       Call Matrix_show(OrbitalLocalizer_instance(speciesID)%densityMatrixA)
       print *, "densityMatrixB"
       Call Matrix_show(OrbitalLocalizer_instance(speciesID)%densityMatrixB)
       call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID), densityMatrixIN=OrbitalLocalizer_instance(speciesID)%densityMatrixB, twoParticlesMatrixOUT=twoParticlesMatrixB)

       !Build Proyection Matrix to force orbital orthogonality between A and B
       do alpha = 1 , numberOfContractions
          do beta = 1 , numberOfContractions
             do mu = 1 , numberOfContractions
                do nu = 1 , numberOfContractions
                   OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(alpha,beta)=&
                        OrbitalLocalizer_instance(speciesID)%projectionMatrix%values(alpha,beta)&
                        +overlapMatrix%values(alpha,mu)*overlapMatrix%values(nu,beta)&
                        *OrbitalLocalizer_instance(speciesID)%densityMatrixB%values(mu,nu)                
                end do
             end do
          end do
       end do

       print *, "projectionMatrix"
       Call Matrix_show(OrbitalLocalizer_instance(speciesID)%projectionMatrix)

       !Adds the projection matrix to the fock matrix
       OrbitalLocalizer_instance(speciesID)%fockMatrixA%values=fockMatrix%values+OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
       print *, "shifted fockMatrix"
       Call Matrix_show(OrbitalLocalizer_instance(speciesID)%fockMatrixA)

    end do

  !!Modify and write molecular system for subsystem A

    print *, "**************************************************************************"
    print *, "*************Molecular subsystem A****************************************"
    print *, "**************************************************************************"

    !For now, we only remove the particles from system B
    do speciesID=1, numberOfSpecies
       MolecularSystem_instance%species(speciesID)%ocupationNumber=&
            MolecularSystem_instance%species(speciesID)%ocupationNumber&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB
       MolecularSystem_instance%species(speciesID)%internalSize=&
            MolecularSystem_instance%species(speciesID)%internalSize&
            -OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInB*2
    end do
    
    call MolecularSystem_saveToFile( "lowdin-subsystemA" )

    call MolecularSystem_destroy()

    !!Load CONTROL Parameters
    call MolecularSystem_loadFromFile( "LOWDIN.DAT", "lowdin-subsystemA" )

    !!Load the system in lowdin.sys format
    call MolecularSystem_loadFromFile( "LOWDIN.SYS", "lowdin-subsystemA" )

    call MolecularSystem_showInformation()  
    call MolecularSystem_showParticlesInformation()

    
    SUBSYSTEM_SCF_CONTINUE=1
    iter=0
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) call system ("lowdin-DFT.x INITIALIZE")
    print *, "Starting subsystem A SCF"
    do while( SUBSYSTEM_SCF_CONTINUE )
       iter=iter+1
    
       do speciesID=1, numberOfSpecies
       
       !Calculates the subsystem orbitals with the new fock matrix

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
          newDensityMatrix%values=0.0
          do k = 1 , OrbitalLocalizer_instance(speciesID)%occupiedOrbitalsInA
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
          densStd= Matrix_standardDeviation(newDensityMatrix,OrbitalLocalizer_instance(speciesID)%densityMatrixA)
          if (densStd .lt. CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE) SUBSYSTEM_SCF_CONTINUE=0
          if (iter .gt. CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) SUBSYSTEM_SCF_CONTINUE=0

          OrbitalLocalizer_instance(speciesID)%densityMatrixA=newDensityMatrix
          print *, "subsystem A iteration", iter, "densityMatrix difference", densStd

       end do
       
       do speciesID=1, numberOfSpecies

          ! print *, "densityMatrixA after"
          ! Call Matrix_show(newDensityMatrix)


          newDensityMatrix%values=OrbitalLocalizer_instance(speciesID)%densityMatrixA%values+OrbitalLocalizer_instance(speciesID)%densityMatrixB%values
          !Updates two particles matrix
          call WaveFunction_buildTwoParticlesMatrix(MolecularSystem_getNameOfSpecie(speciesID), densityMatrixIN=newDensityMatrix, twoParticlesMatrixOUT=twoParticlesMatrixA)
          !Updates coupling matrices
          
          
          !Updates exchange correlation matrix
          if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
             densUnit = 78
             densFile = trim(CONTROL_instance%INPUT_FILE)//"localdensmatrix"
             excUnit = 79
             excFile = trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".excmatrix"
             !!Save density matrix A+B
             open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
             do i = 1, numberOfSpecies
                labels(1) = "DENSITY-MATRIX"
                labels(2) = MolecularSystem_getNameOfSpecie(i)
                if (i .eq. speciesID) then
                   call Matrix_writeToFile(newDensityMatrix, unit=densUnit, binary=.true., arguments = labels )
                else
                   call Matrix_writeToFile(WaveFunction_instance(i)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
                end if
             end do
             close (densUnit)
             !!Run DFT for complete system A+B
             call system("lowdin-DFT.x BUILD_MATRICES "//trim(densFile))

             !!Get results
             open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")
             labels(2) = nameOfSpecies
             labels(1) = "EXCHANGE-CORRELATION-ENERGY"
             call Vector_getFromFile(unit=excUnit, binary=.true., value=excCorrEnergyA, arguments= labels )

             labels(1) = "EXCHANGE-CORRELATION-MATRIX"
             exchangeCorrelationMatrixA=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
                  binary=.true., arguments=labels)

             close(unit=excUnit)


          end if

          !! Updates Fock Matrix
          OrbitalLocalizer_instance(speciesID)%fockMatrixA%values = wavefunction_instance(speciesID)%hcoreMatrix%values &
               +twoParticlesMatrixA%values&
               +wavefunction_instance(speciesID)%couplingMatrix%values&
               +exchangeCorrelationMatrixA%values&
               +OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values

          call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
               OrbitalLocalizer_instance(speciesID)%fockMatrixA, OrbitalLocalizer_instance(speciesID)%densityMatrixA, &
               WaveFunction_instance(speciesID)%OverlapMatrix, &
               methodType=SCF_CONVERGENCE_DAMPING, &
               coefficientMatrix=OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, speciesID=speciesID )

          call Convergence_run( WaveFunction_instance(speciesID)%convergenceMethod )

          print *, "core energy", sum(transpose(newDensityMatrix%values)* WaveFunction_instance( speciesID )%hCoreMatrix%values)
          print *, "repulsion energy", 0.5*sum(transpose(newDensityMatrix%values)* twoParticlesMatrixA%values)
          print *, "excCorrEnergy", excCorrEnergyA
          print *, "energy correction", OrbitalLocalizer_instance(speciesID)%levelShiftingValue*&
               sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixA%values)*&
               OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)

          totalEnergy= sum(transpose(newDensityMatrix%values)* WaveFunction_instance( speciesID )%hCoreMatrix%values) &
               +0.5*sum(transpose(newDensityMatrix%values)* twoParticlesMatrixA%values)&
               +excCorrEnergyA &
               +MolecularSystem_getPointChargesEnergy()

          print *, "uncorrected energy", totalEnergy

          ! print *, "shifted fockMatrix"
          ! Call Matrix_show(OrbitalLocalizer_instance(speciesID)%fockMatrixA)


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
    end do

    do speciesID=1, numberOfSpecies
       write(*,*) ""
       write(*,*) " Subsystem Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
       write(*,*) "-----------------"
       write(*,*) ""

       call Matrix_show( OrbitalLocalizer_instance(speciesID)%waveFunctionCoefficientsA, &
            rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
            columnkeys = string_convertvectorofrealstostring( OrbitalLocalizer_instance(speciesID)%molecularOrbitalsEnergyA ),&
            flags=WITH_BOTH_KEYS)



       print *, "core energy A", sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixA%values)* WaveFunction_instance( speciesID )%hCoreMatrix%values)
       print *, "core energy B", sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixB%values)* WaveFunction_instance( speciesID )%hCoreMatrix%values)

       print *, "repulsion energy A", 0.5*sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixA%values)* twoParticlesMatrixA%values)
       print *, "repulsion energy B", 0.5*sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixB%values)* twoParticlesMatrixB%values)
       print *, "repulsion energy AB", sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixA%values)* twoParticlesMatrixB%values)


       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          densUnit = 78
          densFile = trim(CONTROL_instance%INPUT_FILE)//"localdensmatrix"
          excUnit = 79
          excFile = trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".excmatrix"
          !!Save density matrix A
          open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
          do i = 1, numberOfSpecies
             labels(1) = "DENSITY-MATRIX"
             labels(2) = MolecularSystem_getNameOfSpecie(i)
             if (i .eq. speciesID) then
                call Matrix_writeToFile(OrbitalLocalizer_instance(speciesID)%densityMatrixA, unit=densUnit, binary=.true., arguments = labels )
             else
                call Matrix_writeToFile(WaveFunction_instance(i)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
             end if
          end do
          close (densUnit)
          !!Run DFT for subsystem A
          call system("lowdin-DFT.x BUILD_MATRICES "//trim(densFile))

          !!Get results
          open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")
          labels(2) = nameOfSpecies
          labels(1) = "EXCHANGE-CORRELATION-ENERGY"
          call Vector_getFromFile(unit=excUnit, binary=.true., value=excCorrEnergyA, arguments= labels )

          labels(1) = "EXCHANGE-CORRELATION-MATRIX"
          exchangeCorrelationMatrixA=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
               binary=.true., arguments=labels)

          close(unit=excUnit)

          print *, "excCorrEnergyA", excCorrEnergyA

          !!Save density matrix B
          open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
          do i = 1, numberOfSpecies
             labels(1) = "DENSITY-MATRIX"
             labels(2) = MolecularSystem_getNameOfSpecie(i)
             if (i .eq. speciesID) then
                call Matrix_writeToFile(OrbitalLocalizer_instance(speciesID)%densityMatrixB, unit=densUnit, binary=.true., arguments = labels )
             else
                call Matrix_writeToFile(WaveFunction_instance(i)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
             end if
          end do
          close (densUnit)

          !!Run DFT for subsystem B
          call system("lowdin-DFT.x BUILD_MATRICES "//trim(densFile))

          !!Get results
          open(unit = excUnit, file=trim(excFile), status="old", form="unformatted")
          labels(2) = nameOfSpecies
          labels(1) = "EXCHANGE-CORRELATION-ENERGY"
          call Vector_getFromFile(unit=excUnit, binary=.true., value=excCorrEnergyB, arguments= labels )

          labels(1) = "EXCHANGE-CORRELATION-MATRIX"
          exchangeCorrelationMatrixB=Matrix_getFromFile(unit=excUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4),&
               binary=.true., arguments=labels)

          close(unit=excUnit)

          print *, "excCorrEnergyB", excCorrEnergyB

          print *, "energy correction mu*densityMatrix*projectionMatrix"

          print *, OrbitalLocalizer_instance(speciesID)%levelShiftingValue*&
               sum(transpose(OrbitalLocalizer_instance(speciesID)%densityMatrixA%values)*&
               OrbitalLocalizer_instance(speciesID)%projectionMatrix%values)

       end if

    end do
  
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
       call Matrix_writeToFile(OrbitalLocalizer_instance(speciesID)%densityMatrixA, unit=wfnUnit, binary=.true., arguments = labels )

       ! labels(1) = "HCORE"
       ! call Matrix_writeToFile(WaveFunction_instance(speciesID)%hcoreMatrix, unit=wfnUnit, binary=.true., arguments = labels )

       labels(1) = "ORBITALS"
       call Vector_writeToFile(WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, unit=wfnUnit, binary=.true., arguments = labels )

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
    call Vector_writeToFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%cosmo3Energy, arguments=["COSMO3ENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%totalCouplingEnergy, arguments=["COUPLINGENERGY"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MultiSCF_instance%electronicRepulsionEnergy, arguments=["COUPLING-E-"])

    ! call Vector_writeToFile(unit=wfnUnit, binary=.true., value=MolecularSystem_getPointChargesEnergy(), arguments=["PUNTUALINTERACTIONENERGY"])

    close(wfnUnit)       

    !!calculate HF/KS properties
    call system ("lowdin-CalcProp.x lowdin-subsystemA")

    
  end subroutine OrbitalLocalizer_levelShiftSubSystemOrbitals

  subroutine OrbitalLocalizer_reorderSubSystemOrbitals(speciesID,coefficients,fockMatrix,eigenvalues)
    integer :: speciesID
    type(Matrix) :: coefficients
    type(Matrix) :: fockMatrix
    type(Vector) :: eigenvalues
    
    type(Matrix) :: overlapMatrix
    type(Matrix) :: densityMatrix
    type(Matrix) :: miniDensityMatrix
    type(Vector) :: energyContribution
    real(8) :: sumSubSystem(2), sumAB, normCheck, holdEnergy
    real(8),allocatable ::holdValues(:)

    integer :: numberOfContractions, ocupationNumber
    integer :: index1,index2,holdIndex,holdSystem
    integer :: mu,nu,alpha,beta, i, ii, j, jj, k

    !Removes the projection matrix from the fock matrix
    fockMatrix%values=fockMatrix%values-OrbitalLocalizer_instance(speciesID)%levelShiftingValue*OrbitalLocalizer_instance(speciesID)%projectionMatrix%values
    print *, "transformed fockMatrix"
    Call Matrix_show(fockMatrix)

    
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    overlapMatrix=WaveFunction_instance( speciesID )%OverlapMatrix
    !Asigns molecular orbital to each subsystem and calculates its energy contribution
    call Matrix_constructor (densityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
    call Vector_constructor (energyContribution, numberOfContractions, 0.0_8)
    allocate(holdValues(numberOfContractions))

    print *, "k, eigenvalue, energyContribution, subsystem"
    do k = 1 , numberOfContractions
       normCheck=0.0
       sumSubSystem(:)=0.0
       sumAB=0.0
       call Matrix_constructor (miniDensityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)       
       do mu = 1 , numberOfContractions
          do nu = 1 , numberOfContractions
             miniDensityMatrix%values(mu,nu)=miniDensityMatrix%values(mu,nu)+&
                  coefficients%values(mu,k)*coefficients%values(nu,k)
             if( OrbitalLocalizer_instance(speciesID)%subSystemList(mu) .eq. OrbitalLocalizer_instance(speciesID)%subSystemList(nu) ) then
                sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))=&
                     sumSubSystem(OrbitalLocalizer_instance(speciesID)%subSystemList(mu))+&
                     overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
             else
                sumAB=sumAB+&
                     overlapMatrix%values(mu,nu)*coefficients%values(mu,k)*coefficients%values(nu,k)
             end if
          end do
       end do

       ! if(k .le. ocupationNumber) densityMatrix%values=densityMatrix%values+miniDensityMatrix%values
       
       energyContribution%values(k)= sum(  transpose(miniDensityMatrix%values) &
            *  (  ( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
            + wavefunction_instance(speciesID)%twoParticlesMatrix%values &
            + wavefunction_instance(speciesID)%couplingMatrix%values &
            + WaveFunction_instance(speciesID)%exchangeCorrelationMatrix%values ))
       
       if(sumSubSystem(1)**2 .gt. sumSubSystem(2)**2) then
          OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 1
       else
          OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k) = 2
       end if
       print *, k, eigenvalues%values(k), energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k)
    end do

    !Reorders coefficients matrix and eigenvalues according to the energy contribution
        
    do index1 = 1 , numberOfContractions-1
       holdIndex=index1
       holdEnergy=energyContribution%values(index1)
       holdValues=coefficients%values(:,index1)
       holdSystem=OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(index1) 

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
       OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(index1)=OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(holdIndex)
       OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(holdIndex)=holdSystem

    end do
    
    print *, "sorted orbitals, subSystem"
    do k=1, numberOfContractions
       print *, k, energyContribution%values(k), OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k)
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
    !    if(OrbitalLocalizer_instance(speciesID)%orbitalSubSystem(k).eq.2) &
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
    
  end subroutine OrbitalLocalizer_reorderSubSystemOrbitals
  


end module OrbitalLocalizer_
