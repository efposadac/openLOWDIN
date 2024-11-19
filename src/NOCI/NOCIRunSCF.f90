!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                
!!	
!!	this program has been developed under direction of:
!!
!!	  UNIVERSIDAD NACIONAL DE COLOMBIA"
!!	  PROF. ANDRES REYES GROUP"
!!	  http://www.qcc.unal.edu.co"
!!	
!!	  UNIVERSIDAD DE GUADALAJARA"
!!	  PROF. ROBERTO FLORES GROUP"
!!	  http://www.cucei.udg.mx/~robertof"
!!	
!!	AUTHORS
!!		E.F. POSADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		S.A. GONZALEZ. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		F.S. MONCADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		J. ROMERO. UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!	CONTRIBUTORS
!!		N.F.AGUIRRE. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		GABRIEL MERINO. UNIVERSIDAD DE GUANAJUATO
!!   		J.A. CHARRY UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module NOCIRunSCF_
  use NOCIBuild_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use DirectIntegralManager_
  use MultiSCF_
  use WaveFunction_
  use omp_lib
  implicit none

  !>
  !! @brief non Orthogonal Configuration Interaction Module. APMO implementation of Skone et al 2005 10.1063/1.2039727
  !!
  !! @author Felix
  !!
  !! <b> Creation data : </b> 02-22
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 02-22 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# Creation of the module.
  !!
  !<

  public :: &
       NOCIRunSCF_runHFs

  private

contains
  
  !>
  !! @brief Run a Hartree-Fock calculation at displaced geometries and fill CI matrix diagonals 
  !!
  !! @param this -> NOCI instance
  !<
  subroutine NOCIRunSCF_runHFs(this)
    implicit none
    type(NonOrthogonalCI) :: this

    integer, allocatable :: sysIbatch(:)
    type(MultiSCF), allocatable :: MultiSCFParallelInstance(:)
    type(WaveFunction), allocatable :: WaveFunctionParallelInstance(:,:)
    type(Libint2Interface), allocatable :: Libint2ParallelInstance(:,:)
    integer :: speciesID, otherSpeciesID, nspecies
    integer :: sysI,me,mySysI
    integer :: ncores, batchSize
    integer :: coordsUnit
    real(8) :: timeA
    character(100) :: coordsFile
    character(50) :: auxmethod

    !$  timeA = omp_get_wtime()
    !!Read HF energy of the non displaced SCF calculation 
    ! print *, "HF reference energy is ", hfEnergy
    nspecies=molecularSystem_instance%numberOfQuantumSpecies 

    allocate(this%HFCoefficients(this%numberOfDisplacedSystems,nspecies))
    allocate(this%systemLabels(this%numberOfDisplacedSystems))

    call Matrix_constructor(this%configurationHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Matrix_constructor(this%configurationOverlapMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Vector_constructor(this%configurationCorrelationEnergies, this%numberOfDisplacedSystems, 0.0_8)
    do speciesID=1, nspecies
         call Matrix_constructor(this%configurationKineticMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationPuntualMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationExternalMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationExchangeMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         do otherSpeciesID=speciesID, nspecies
            call Matrix_constructor(this%configurationHartreeMatrix(speciesID,otherSpeciesID), &
                 int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
            call Matrix_constructor(this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID), &
                 int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         end do
    end do

    coordsUnit=333
    coordsFile=trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords"
    open(unit=coordsUnit, file=trim(coordsFile), status="replace", form="formatted")

    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       print *, "running KS calculations at the displaced geometries ... saving results on file ", coordsFile
    else
       print *, "running HF calculations at the displaced geometries ... saving results on file ", coordsFile
    end if


    write (coordsUnit,'(A25,I20)') "numberOfDisplacedSystems ", this%numberOfDisplacedSystems

    !Allocate objets to distribute in parallel
    ncores=CONTROL_instance%NUMBER_OF_CORES
    batchSize=min(ncores,this%numberOfDisplacedSystems)
    print *, "ncores", ncores, "batchsize", batchSize

    call MolecularSystem_destroy()
    
    !Skip priting scf iterations
    CONTROL_instance%PRINT_LEVEL=0
    CONTROL_instance%NUMBER_OF_CORES=1
    
    allocate(sysIbatch(batchSize),&
         MultiSCFParallelInstance(batchSize),&
         WaveFunctionParallelInstance(nspecies,batchSize),&
         Libint2ParallelInstance(nspecies,batchSize))
    
    sysI=0
    systemLoop: do while(sysI.le.this%numberOfDisplacedSystems)
       !In serial, prepare systems
       sysIbatch(:)=0
       me=0
       mySysI=sysI

       do while(me.lt.batchSize)
          mySysI=mySysI+1
          if(mySysI .gt. this%numberOfDisplacedSystems) exit
          me=me+1
          sysIbatch(me)=mySysI

          write(this%systemLabels(mySysI), '(A)') trim(this%molecularSystems(mySysI)%description)
          call MultiSCF_constructor(MultiSCFParallelInstance(me),WaveFunctionParallelInstance(1:nspecies,me),CONTROL_instance%ITERATION_SCHEME,this%molecularSystems(mySysI))
          call MultiSCF_buildHcore(MultiSCFParallelInstance(me),WaveFunctionParallelInstance(1:nspecies,me))       
          call MultiSCF_getInitialGuess(MultiSCFParallelInstance(me),WaveFunctionParallelInstance(1:nspecies,me))
          call DirectIntegralManager_constructor(Libint2ParallelInstance(1:nspecies,me),this%molecularSystems(mySysI))
          
       end do
       ! STOP "NOCI runs only work with CONTROL_instance%INTEGRAL_STORAGE == MEMORY"

       !In parallel, run SCF calculations without calling lowdin-scf.x      
       call OMP_set_num_threads(ncores)
       !$omp parallel&
       !$omp& private(mySysI,auxmethod,speciesID,otherSpeciesID)
       !$omp do schedule(dynamic)
       procs: do me=1, batchSize
          mySysI=sysIbatch(me)
          if(mySysI .eq. 0) cycle procs
          
          if (CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then
             do speciesID=1, nspecies
                call DirectIntegralManager_getDirectIntraRepulsionIntegralsAll(&
                     speciesID, &
                     WaveFunctionParallelInstance(speciesID,me)%densityMatrix, & 
                     WaveFunctionParallelInstance(speciesID,me)%fourCenterIntegrals(speciesID)%values, &
                     this%molecularSystems(mySysI),Libint2ParallelInstance(speciesID,me))
             end do

             do speciesID=1, nspecies-1
                do otherSpeciesID=speciesID+1,nspecies
                   call DirectIntegralManager_getDirectInterRepulsionIntegralsAll(&
                        speciesID, otherSpeciesID, &
                        WaveFunctionParallelInstance(speciesID,me)%densityMatrix, & 
                        WaveFunctionParallelInstance(speciesID,me)%fourCenterIntegrals(otherSpeciesID)%values, &
                        this%molecularSystems(mySysI),Libint2ParallelInstance(speciesID,me),Libint2ParallelInstance(otherSpeciesID,me))
                end do
             end do
          end if

          call MultiSCF_solveHartreeFockRoothan(MultiSCFParallelInstance(me),WaveFunctionParallelInstance(1:nspecies,me),Libint2ParallelInstance(1:nspecies,me))

          !Save HF results
          ! call MultiSCF_saveWfn(MultiSCF_instance,WaveFunction_instance)
          ! call MolecularSystem_copyConstructor(this%molecularSystems(sysI),molecularSystem_instance)
          this%configurationHamiltonianMatrix%values(mySysI,mySysI)=MultiSCFParallelInstance(me)%totalEnergy
          
          do speciesID = 1, nspecies
             this%HFCoefficients(mySysI,speciesID) = WaveFunctionParallelInstance(speciesID,me)%waveFunctionCoefficients
             this%configurationKineticMatrix(speciesID)%values(mySysI,mySysI)=WaveFunctionParallelInstance(speciesID,me)%kineticEnergy
             this%configurationPuntualMatrix(speciesID)%values(mySysI,mySysI)=WaveFunctionParallelInstance(speciesID,me)%puntualInteractionEnergy
             this%configurationExternalMatrix(speciesID)%values(mySysI,mySysI)=WaveFunctionParallelInstance(speciesID,me)%externalPotentialEnergy
             this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysI)=WaveFunctionParallelInstance(speciesID,me)%exchangeHFEnergy
             do otherSpeciesID = speciesID, this%molecularSystems(mySysI)%numberOfQuantumSpecies
                this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysI)=&
                     WaveFunctionParallelInstance(speciesID,me)%hartreeEnergy(otherSpeciesID)
                this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysI)=&
                     WaveFunctionParallelInstance(speciesID,me)%exchangeCorrelationEnergy(otherSpeciesID)                  
             end do
          end do

          ! Compute HF energy with KS determinants
          if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
             if(CONTROL_instance%METHOD.eq."RKS") then
                auxmethod="RHF"
             else
                auxmethod="UHF"
             end if

             do speciesID = 1, nspecies
                WaveFunctionParallelInstance(speciesID,me)%exchangeCorrelationEnergy=0.0_8
                WaveFunctionParallelInstance(speciesID,me)%exchangeCorrelationMatrix%values=0.0_8
                this%exactExchangeFraction(speciesID)=WaveFunctionParallelInstance(speciesID,me)%exactExchangeFraction
                WaveFunctionParallelInstance(speciesID,me)%exactExchangeFraction=1.0_8
             end do
             call MultiSCF_obtainFinalEnergy(MultiSCFParallelInstance(me),WaveFunctionParallelInstance(1:nspecies,me),Libint2ParallelInstance(1:nspecies,me),auxmethod)
             !Difference between HF and KS energies
             this%configurationCorrelationEnergies%values(mySysI)=this%configurationHamiltonianMatrix%values(mySysI,mySysI)-MultiSCFParallelInstance(me)%totalEnergy
          end if
       end do procs
       !$omp end do nowait
       !$omp end parallel

       !In serial, free memory and print
       do me=1, batchSize
          mySysI=sysIbatch(me)
          if(mySysI .eq. 0) exit systemLoop

          write (coordsUnit,'(A10,I10,A10,ES20.12,A20,ES20.12)') "Geometry ", mySysI, "Energy", this%configurationHamiltonianMatrix%values(mySysI,mySysI), &
               "Correlation energy", this%configurationCorrelationEnergies%values(mySysI)               
          call MolecularSystem_showCartesianMatrix(this%molecularSystems(mySysI),unit=coordsUnit)

          if (this%numberOfDisplacedSystems .le. this%printMatrixThreshold) then 
             write (*,'(A10,I10,A10,ES20.12,A20,ES20.12)') "Geometry ", mySysI, "Energy", this%configurationHamiltonianMatrix%values(mySysI,mySysI), &
                  "Correlation energy", this%configurationCorrelationEnergies%values(mySysI)               
             call MolecularSystem_showCartesianMatrix(this%molecularSystems(mySysI))
             ! do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
             !    print *, "sysI", sysI, "speciesID", speciesID, "occupation number", MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(mySysI))
             ! end do
          end if

          call DirectIntegralManager_destructor(Libint2ParallelInstance(1:nspecies,me))
          call MultiSCF_destructor(MultiSCFParallelInstance(me))

          sysI=mySysI

       end do

       CONTROL_instance%NUMBER_OF_CORES=ncores

       !!Screen geometries with high energies
       ! if( CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD .ne. 0.0 .and. &
       !      testEnergy .gt. this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD) then
       !    write (coordsUnit,"(A,F20.12)") "Skipping system with high energy", testEnergy
       !    this%numberOfEnergyRejectedSystems=this%numberOfEnergyRejectedSystems+1                      
       ! else
       ! if(this%numberOfEnergyRejectedSystems .gt. 0) &
       !      write (*,'(A10,I10,A,F18.12)') "Rejected ", this%numberOfEnergyRejectedSystems, &
       !      " geometries with energy higher than", this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD       

    end do systemLoop

    close(coordsUnit)
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for HF calculations at displaced geometries : ", omp_get_wtime() - timeA ," (s)"
    
  end subroutine NOCIRunSCF_runHFs
  
end module NOCIRunSCF_

