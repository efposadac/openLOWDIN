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

module NOCIMatrices_
  use NOCIBuild_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use DirectIntegralManager_
  use Libint2Interface_
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
       NOCIMatrices_buildOverlapAndHamiltonian,&
       NOCIMatrices_diagonalize,&
       NOCIMatrices_mergeCoefficients

  private

contains
  
  !>
  !! @brief Computes overlap and hamiltonian non orthogonal CI matrices for previously calculated molecular systems at different geometries
  !!
  !! @param this
  !<
  subroutine NOCIMatrices_buildOverlapAndHamiltonian(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem), allocatable :: mergedMolecularSystem(:)
    type(Libint2Interface), allocatable :: Libint2ParallelInstance(:,:)
    integer, allocatable :: sysIbatch(:), sysIIbatch(:)
    integer :: sysI,sysII,me,mySysI,mySysII
    type(Matrix), allocatable :: mergedCoefficients(:), inverseOverlapMatrices(:)
    type(IVector), allocatable :: sysIbasisList(:,:),sysIIbasisList(:,:)
    real(8) :: overlapUpperBound
    integer :: prescreenedElements, overlapScreenedElements

    integer :: speciesID, otherSpeciesID
    integer :: nspecies
    integer :: ncores, batchSize, upperBound
    
    integer :: matrixUnit
    character(100) :: matrixFile
    real(8) :: empiricalScaleFactor
    
    real(8) :: timeMerging, timePrescreen, timeOverlap, timeTwoIntegrals
    real(8) :: timeA
    real(8) :: timeB

    timePrescreen=0.0
    timeOverlap=0.0
    timeTwoIntegrals=0.0
        
    print *, ""
    print *, "A prescreening of the overlap matrix elements is performed for the heavy species"
    write (*,'(A,ES8.1)') "Overlap and Hamiltonian matrix elements are saved for pairs with overlap higher than",&
         CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD
    print *, "For pairs with lower overlap, setting H(I,II)=0, S(I,II)=0"
    print *, ""

    prescreenedElements=0
    overlapScreenedElements=0

    matrixUnit=290
    matrixFile= trim(CONTROL_instance%INPUT_FILE)//"NOCI-Matrix.ci"

    print *, "computing NOCI overlap and hamiltonian matrices... saving them to ", trim(matrixFile)

    open(unit=matrixUnit, file=trim(matrixFile), status="replace", form="formatted")

    write (matrixUnit,'(A20,I20)') "MatrixSize", this%numberOfDisplacedSystems
    write (matrixUnit,'(A10,A10,A20,A20)') "Conf. ", "Conf. ", "Overlap ","Hamiltonian "
    !Save diagonal elements
    do sysI=1,this%numberOfDisplacedSystems
       this%configurationOverlapMatrix%values(sysI,sysI)=1.0
       write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, sysI, &
            this%configurationOverlapMatrix%values(sysI,sysI), this%configurationHamiltonianMatrix%values(sysI,sysI)               
    end do
    
    !Allocate objets to distribute in parallel
    nspecies=this%molecularSystems(1)%numberOfQuantumSpecies 
    ncores=CONTROL_instance%NUMBER_OF_CORES
    batchSize=this%numberOfDisplacedSystems
    print *, "ncores", ncores, "batchsize", batchSize

    allocate(mergedMolecularSystem(batchSize),&
         mergedCoefficients(nspecies),&
         inverseOverlapMatrices(nspecies),&
         Libint2ParallelInstance(nspecies,batchSize),&
         sysIbatch(batchSize),&
         sysIIbatch(batchSize),&
         sysIbasisList(nspecies,batchSize),&
         sysIIbasisList(nspecies,batchSize))

    if(CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS) then
       upperBound=1
       this%printMatrixThreshold=this%numberOfDisplacedSystems       
    else
       upperBound=this%numberOfDisplacedSystems       
    end if
    ! print *, "upperBound", upperBound

    sysI=1
    sysII=1
    
    systemLoop: do while((sysI.le.upperBound .and. sysII.le.this%numberOfDisplacedSystems))
       ! print *, "distributing sysI ", sysI,  " sysII ", sysII, " into", batchSize, " batches"          
       !In serial, prepare systems
       sysIbatch(:)=0
       sysIIbatch(:)=0
       me=0
       mySysI=sysI
       mySysII=sysII

       do while(me.lt.batchSize)
          mySysII=mySysII+1
          if(mySysII .gt. this%numberOfDisplacedSystems) then
             mySysI=mySysI+1
             mySysII=mySysI+1
             if(mySysI .gt. upperBound .or. mySysII .gt. this%numberOfDisplacedSystems) exit
          end if

          ! print *, "checking prescreening of elements", mySysI, mySysII
          !$  timeA = omp_get_wtime()
          !Estimates overlap with a 1s-1s integral approximation
          call NOCIMatrices_prescreenOverlap(this,mySysI,mySysII,overlapUpperBound)

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               overlapUpperBound .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             ! print *, "preskipping elements", mySysI, mySysII, "with overlap estimated as", overlapUpperBound
             prescreenedElements=prescreenedElements+1
          else
             !$  timeB = omp_get_wtime()
             !$  timePrescreen=timePrescreen+(timeB - timeA)
             me=me+1
             sysIbatch(me)=mySysI
             sysIIbatch(me)=mySysII
             !$  timeA = omp_get_wtime()
             !This generates a new molecular system
             ! print *, "Merging systems from geometries ", mySysI, mySysII
             call MolecularSystem_mergeTwoSystems(mergedMolecularSystem(me), &
                  this%molecularSystems(mySysI), this%molecularSystems(mySysII),sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me))
             ! call MolecularSystem_showInformation()  
             ! call MolecularSystem_showParticlesInformation()
             ! call MolecularSystem_showCartesianMatrix(mergedMolecularSystem)
             call DirectIntegralManager_constructor(Libint2ParallelInstance(1:nspecies,me),mergedMolecularSystem(me))
             !$  timeB = omp_get_wtime()
             !$  timeMerging=timeMerging+(timeB - timeA)
             
          end if
       end do

       !In parallel, fill matrices

       call OMP_set_num_threads(ncores)
       !$omp parallel & 
       !$omp& private(mySysI,mySysII,mergedCoefficients,inverseOverlapMatrices),&
       !$omp& shared(this,sysI,sysII,matrixUnit,prescreenedElements,overlapScreenedElements,sysIbasisList,sysIIbasisList,mergedMolecularSystem,Libint2ParallelInstance,nspecies,batchSize)
       !$omp do schedule(dynamic,10)
       procs: do me=1, batchSize
          mySysI=sysIbatch(me)
          mySysII=sysIIbatch(me)          
          if(mySysI .eq. 0 .or. mySysII .eq. 0) cycle procs
          
          ! print *, "evaluating S and H elements for", mySysI, mySysII

          !! Merge occupied coefficients into a single matrix 
          call NOCIMatrices_mergeCoefficients(this%HFCoefficients(mySysI,:),this%HFCoefficients(mySysII,:),&
               this%molecularSystems(mySysI),this%molecularSystems(mySysII),mergedMolecularSystem(me),&
               sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me),mergedCoefficients)
          !$  timeA = omp_get_wtime()

          call NOCIMatrices_computeOverlapAndHCoreElements(this,mySysI,mySysII,mergedMolecularSystem(me),mergedCoefficients,&
               sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me),inverseOverlapMatrices)
          !$  timeB = omp_get_wtime()
          !$  timeOverlap=timeOverlap+(timeB - timeA)

          !! SKIP ENERGY EVALUATION IF OVERLAP IS TOO LOW

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               abs(this%configurationOverlapMatrix%values(mySysI,mySysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             ! print *, "screening elements", mySysI, mySysII, "with overlap", this%configurationOverlapMatrix%values(mySysI,mySysII)
             this%configurationOverlapMatrix%values(mySysI,mySysII)=0.0
             this%configurationHamiltonianMatrix%values(mySysI,mySysII)=0.0
             !$OMP ATOMIC
             overlapScreenedElements=overlapScreenedElements+1
          else

             !$  timeA = omp_get_wtime()
             ! print *, "evaluating twoParticlesContributions for", mySysI, mySysII
             call NOCIMatrices_twoParticlesContributions(this,mySysI,mySysII,mergedMolecularSystem(me),&
                  inverseOverlapMatrices,mergedCoefficients,Libint2ParallelInstance(1:nspecies,me))

             if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
                !    !DFT energy correction for off diagonal elements following Gao2016 - scaled average of the diagonal elements

                do speciesID = 1, nspecies
                   this%configurationHamiltonianMatrix%values(mySysI,mySysII)=this%configurationHamiltonianMatrix%values(mySysI,mySysII)-&
                        this%configurationOverlapMatrix%values(mySysI,mySysII)/2.0*&
                        (1/this%exactExchangeFraction(speciesID)-1)*&
                        (this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysI)+&
                        this%configurationExchangeMatrix(speciesID)%values(mySysII,mySysII))

                   this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)=this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)-&
                        this%configurationOverlapMatrix%values(mySysI,mySysII)/2.0*&
                        (1/this%exactExchangeFraction(speciesID)-1)*&
                        (this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysI)+&
                        this%configurationExchangeMatrix(speciesID)%values(mySysII,mySysII))

                   do otherSpeciesID = speciesID, nspecies
                      this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)=this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)+&
                           this%configurationOverlapMatrix%values(mySysI,mySysII)/2.0*&
                           (this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysI)+&
                           this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysII,mySysII))
                      this%configurationHamiltonianMatrix%values(mySysI,mySysII)=this%configurationHamiltonianMatrix%values(mySysI,mySysII)+&
                           this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)
                   end do
                end do

                ! this%configurationHamiltonianMatrix%values(mySysI,mySysII)=this%configurationHamiltonianMatrix%values(mySysI,mySysII)+&
                !      this%configurationOverlapMatrix%values(mySysI,mySysII)/2.0*&
                !      (this%configurationCorrelationEnergies%values(mySysI)+&
                !      this%configurationCorrelationEnergies%values(mySysII))
                !    !DFT energy correction for off diagonal elements
                !    call NOCIMatrices_getOffDiagonalDensityMatrix(this,mySysI,mySysII,mergedCoefficients,mergedMolecularSystem(me),this%configurationOverlapMatrix%values(mySysI,mySysII),&
                !         inverseOverlapMatrices,sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me))
             end if

             !$  timeB = omp_get_wtime()
             !$  timeTwoIntegrals=timeTwoIntegrals+(timeB - timeA)
          end if

          ! print *, "thread", omp_get_thread_num()+1,"me", me, "mySysI", " mySysII", mySysI, mySysII, "S", this%configurationOverlapMatrix%values(mySysI,mySysII), "H", this%configurationHamiltonianMatrix%values(mySysI,mySysII)
       end do procs
       !$omp end do nowait
       !$omp end parallel

       !In serial, symmetrize, free memory and print
       do me=1, batchSize
          mySysI=sysIbatch(me)
          mySysII=sysIIbatch(me)

          if(mySysI .eq. 0 .or. mySysII .eq. 0) exit systemLoop

          !Yu2020 magical empirical correction
          if(CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION .and. &
               abs(this%configurationOverlapMatrix%values(mySysI,mySysII)) .gt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             empiricalScaleFactor=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A*&
                  abs(this%configurationOverlapMatrix%values(mySysI,mySysII))**CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B/&
                  abs(this%configurationOverlapMatrix%values(mySysI,mySysII))
             this%configurationOverlapMatrix%values(mySysI,mySysII)=&
                  this%configurationOverlapMatrix%values(mySysI,mySysII)*empiricalScaleFactor
             this%configurationHamiltonianMatrix%values(mySysI,mySysII)=&
                  this%configurationHamiltonianMatrix%values(mySysI,mySysII)*empiricalScaleFactor                     
             do speciesID=1, nspecies
                this%configurationKineticMatrix(speciesID)%values(mySysI,mySysII)=&
                     this%configurationKineticMatrix(speciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                this%configurationPuntualMatrix(speciesID)%values(mySysI,mySysII)=&
                     this%configurationPuntualMatrix(speciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                this%configurationExternalMatrix(speciesID)%values(mySysI,mySysII)=&
                     this%configurationExternalMatrix(speciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)=&
                     this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                do otherSpeciesID=speciesID, nspecies
                   this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)=&
                        this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                   this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)=&
                        this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)*empiricalScaleFactor
                end do
             end do
          end if

          !Symmetrize
          this%configurationOverlapMatrix%values(mySysII,mySysI)=this%configurationOverlapMatrix%values(mySysI,mySysII)
          this%configurationHamiltonianMatrix%values(mySysII,mySysI)=this%configurationHamiltonianMatrix%values(mySysI,mySysII)

          do speciesID=1, nspecies
             this%configurationKineticMatrix(speciesID)%values(mySysII,mySysI)=this%configurationKineticMatrix(speciesID)%values(mySysI,mySysII)
             this%configurationPuntualMatrix(speciesID)%values(mySysII,mySysI)=this%configurationPuntualMatrix(speciesID)%values(mySysI,mySysII)
             this%configurationExternalMatrix(speciesID)%values(mySysII,mySysI)=this%configurationExternalMatrix(speciesID)%values(mySysI,mySysII)
             this%configurationExchangeMatrix(speciesID)%values(mySysII,mySysI)=this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)
             do otherSpeciesID=speciesID, nspecies
                this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysII,mySysI)=&
                     this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)
                this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysII,mySysI)=&
                     this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)
             end do
          end do

          write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') mySysI, mySysII, &
               this%configurationOverlapMatrix%values(mySysI,mySysII), this%configurationHamiltonianMatrix%values(mySysI,mySysII)               

          if (this%numberOfDisplacedSystems .le. this%printMatrixThreshold) then 
             write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, "Overlap element = ", this%configurationOverlapMatrix%values(mySysI,mySysII)
             do speciesID = 1, nspecies                
                write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     " Kinetic element = ", this%configurationKineticMatrix(speciesID)%values(mySysI,mySysII)
                write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     " Puntual element = ", this%configurationPuntualMatrix(speciesID)%values(mySysI,mySysII)
                if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                     write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     " External element = ", this%configurationExternalMatrix(speciesID)%values(mySysI,mySysII)
                write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     "/"//trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     " Hartree element = ", this%configurationHartreeMatrix(speciesID,speciesID)%values(mySysI,mySysII)
                write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                     " Exchange element = ", this%configurationExchangeMatrix(speciesID)%values(mySysI,mySysII)
             end do
             do speciesID=1, nspecies-1
                do otherSpeciesID=speciesID+1, nspecies
                   write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                        "/"//trim( this%molecularSystems(mySysI)%species(otherSpeciesID)%name ) // &
                        " Hartree element = ", this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)
                end do
             end do
             if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
                do speciesID=1, nspecies
                   write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                        "/"//trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                        " DFTcorrelation element = ", this%configurationDFTcorrelationMatrix(speciesID,speciesID)%values(mySysI,mySysII)
                end do
                do speciesID=1, nspecies
                   do otherSpeciesID=speciesID+1, nspecies-1
                      write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, trim( this%molecularSystems(mySysI)%species(speciesID)%name ) // &
                           "/"//trim( this%molecularSystems(mySysI)%species(otherSpeciesID)%name ) // &
                           " DFTcorrelation element = ", this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysI,mySysII)
                   end do
                end do
                
                ! write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, "Total DFT Correlation element = ", this%configurationOverlapMatrix%values(mySysI,mySysII)/2.0*&
                !      (this%configurationCorrelationEnergies%values(mySysI)+&
                !      this%configurationCorrelationEnergies%values(mySysII))
             end if
             write (*,'(I10,I10,A38,ES20.12)') mySysI, mySysII, "Hamiltonian element = ", this%configurationHamiltonianMatrix%values(mySysI,mySysII)
             print *, ""

          end if

          call DirectIntegralManager_destructor(Libint2ParallelInstance(1:nspecies,me))

          sysI=mySysI
          sysII=mySysII

       end do

    end do systemLoop

    close(matrixUnit)

    print *, ""
    print *, "Configuration pairs skipped by overlap prescreening: ", prescreenedElements
    print *, "Configuration pairs skipped by overlap    screening: ", overlapScreenedElements
    if( .not. CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS) then
       print *, "Overlap integrals computed for    ", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
            -prescreenedElements, "configuration pairs"
       print *, "Four center integrals computed for", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
            -prescreenedElements-overlapScreenedElements, "configuration pairs"
    else
       print *, "Overlap integrals computed for    ", this%numberOfDisplacedSystems&
            -prescreenedElements, "configuration pairs"
       print *, "Four center integrals computed for", this%numberOfDisplacedSystems&
            -prescreenedElements-overlapScreenedElements, "configuration pairs"
    end if
    print *, ""

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for overlap prescreening : ", timePrescreen ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for merging systems      : ", timeMerging ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for two index integrals  : ", timeOverlap ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for four index integrals : ", timeTwoIntegrals ," (s)"

    print *, ""

    deallocate(mergedMolecularSystem,&
         mergedCoefficients,&
         inverseOverlapMatrices,&
         Libint2ParallelInstance,&
         sysIbatch,&
         sysIIbatch,&
         sysIbasisList,&
         sysIIbasisList)
    
    ! integer :: symmetryEquivalentElements
    ! timeSymmetry=0.0
    ! symmetryEquivalentElements=0
    ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
    !    write (matrixUnit,'(A10,A10,A10,A20,A20)') "Conf. ", "Conf. ", "Type ", "Overlap ","Hamiltonian "
    ! else
    ! end if
    ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
       !    write (matrixUnit,'(I10,I10,I10,ES20.12,ES20.12)') sysI, sysI, this%configurationPairTypes%values(sysI,sysI), &
       !         this%configurationOverlapMatrix%values(sysI,sysI), this%configurationHamiltonianMatrix%values(sysI,sysI)               
       ! else
       ! end if
    ! write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for element symmetry     : ", timeSymmetry ," (s)"
          ! !$  timeA = omp_get_wtime()
          ! !!Check symmetry of the element
          ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
          !    call NOCIMatrices_classifyConfigurationPair(this,sysI,sysII,newPairFlag)
          !    !$  timeB = omp_get_wtime()
          !    !$  timeSymmetry=timeSymmetry+(timeB - timeA)

          !    !!Copy results from previously computed equivalent elements
          !    if (newPairFlag .eqv. .false.) then
          !       do preSysI=1, sysI
          !          do preSysII=preSysI+1, sysII                   
          !             if(this%configurationPairTypes%values(preSysI,preSysII) .eq. this%configurationPairTypes%values(sysI,sysII)) then
          !                this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(preSysI,preSysII)
          !                this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)
          !                this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(preSysI,preSysII)
          !                this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)
          !                symmetryEquivalentElements=symmetryEquivalentElements+1

          !                if( this%configurationOverlapMatrix%values(sysI,sysII) .ne. 0.0) &
          !                     write (*,'(A,I10,I10,A,I10,A,ES20.12,ES20.12)') "Pair ",sysI, sysII," is type ", &
          !                     this%configurationPairTypes%values(sysI,sysII), " Overlap and Hamiltonian elements", &
          !                     this%configurationOverlapMatrix%values(sysI,sysII), this%configurationHamiltonianMatrix%values(sysI,sysII)

          !                cycle systemII
          !             end if
          !          end do
          !       end do
          !    end if
          ! end if
          !!This is a symmetry test, assume positive phase
          ! if( this%configurationOverlapMatrix%values(sysI,sysII) .lt. 0.0) then
          !    this%configurationOverlapMatrix%values(sysI,sysII)=-this%configurationOverlapMatrix%values(sysI,sysII)
          !    this%configurationOverlapMatrix%values(sysII,sysI)=-this%configurationOverlapMatrix%values(sysII,sysI)
          !    this%configurationHamiltonianMatrix%values(sysI,sysII)=-this%configurationHamiltonianMatrix%values(sysI,sysII)
          !    this%configurationHamiltonianMatrix%values(sysII,sysI)=-this%configurationHamiltonianMatrix%values(sysII,sysI)
          ! end if
          
          ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
          !    write (matrixUnit,'(I10,I10,I10,ES20.12,ES20.12)') sysI, sysII, this%configurationPairTypes%values(sysI,sysII), &
          !         this%configurationOverlapMatrix%values(sysI,sysII), this%configurationHamiltonianMatrix%values(sysI,sysII)               
          ! else
    ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
    !      print *, "Configuration pairs skipped by symmetry equivalence: ", symmetryEquivalentElements

  end subroutine NOCIMatrices_buildOverlapAndHamiltonian

  
  !>
  !! @brief Merges the occupied orbitals coefficients from two systems
  !! @param occupationI and occupationII: Number of orbitals to merge from each matrix. 
  !! sysBasisList: array indicating which basis functions of the merged molecular system belong to sysI and sysII Merged Coefficients: Matrices for output.
  !<
  subroutine NOCIMatrices_mergeCoefficients(coefficientsI,coefficientsII,molecularSystemI,molecularSystemII,mergedMolecularSystem,&
       sysIbasisList,sysIIbasisList,mergedCoefficients)
    type(Matrix), intent(in) :: coefficientsI(*), coefficientsII(*)
    type(MolecularSystem), intent(in) :: molecularSystemI, molecularSystemII, mergedMolecularSystem
    type(IVector), intent(in) :: sysIbasisList(*), sysIIbasisList(*)
    type(Matrix), intent(out) :: mergedCoefficients(*)
    
    ! character(100) :: wfnFile
    ! character(50) :: arguments(2)
    ! integer :: wfnUnit
    integer :: speciesID, i, j, mu
    
    !! Mix coefficients of occupied orbitals of both systems    
    !!Create a dummy density matrix to lowdin.wfn file
    ! wfnUnit = 500
    ! wfnFile = "lowdin.wfn"
    ! open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    do speciesID = 1, mergedMolecularSystem%numberOfQuantumSpecies
       
       ! arguments(2) = mergedMolecularSystem%species(speciesID)%name

       ! arguments(1) = "COEFFICIENTS"

       !    !Max: to make the matrix square for the integral calculations for configuration pairs, and rectangular for the merged coefficients of all systems 
       call Matrix_constructor(mergedCoefficients(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8), &
            int(max(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),MolecularSystem_getOcupationNumber(speciesID,mergedMolecularSystem)),8), 0.0_8 )

       ! print *, "sysI coefficients for ", speciesID
       ! call Matrix_show(coefficientsI(speciesID))
       ! print *, "sysII coefficients for ", speciesID
       ! call Matrix_show(coefficientsII(speciesID))

       !sysI orbitals on the left columns, sysII on the right columns
       !sysI coefficients
       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
          if((sysIbasisList(speciesID)%values(mu) .ne. 0) ) then
             do i=1, MolecularSystem_getOcupationNumber(speciesID,molecularSystemI)!sysI
                mergedCoefficients(speciesID)%values(mu,i)=coefficientsI(speciesID)%values(sysIbasisList(speciesID)%values(mu),i)
                ! print *, "sys I", mu, i, mergedCoefficients(speciesID)%values(mu,i)
             end do
          end if
       end do

       ! !sysII coefficients
       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
          if((sysIIbasisList(speciesID)%values(mu) .ne. 0) ) then
             do i=1, MolecularSystem_getOcupationNumber(speciesID,molecularSystemII)!sysII
                j=MolecularSystem_getOcupationNumber(speciesID,molecularSystemI)+i !column
                mergedCoefficients(speciesID)%values(mu,j)=coefficientsII(speciesID)%values(sysIIbasisList(speciesID)%values(mu),i)
                ! print *, "sys II", mu, j, mergedCoefficients(speciesID)%values(mu,j)
             end do
          end if
       end do       
              
       ! print *, "Merged coefficients matrix for ", speciesID
       ! call Matrix_show(mergedCoefficients(speciesID))

       ! call Matrix_writeToFile(mergedCoefficients(speciesID), unit=wfnUnit, binary=.true., arguments = arguments      
       ! call Matrix_writeToFile(auxMatrix, unit=wfnUnit, binary=.true., arguments = arguments )

       ! arguments(1) = "ORBITALS"
       ! call Vector_constructor(auxVector, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem), 0.0_8 )

       ! call Vector_writeToFile(auxVector, unit=wfnUnit, binary=.true., arguments = arguments )

       ! Only occupied orbitals are going to be transformed - handled in integral transformation program
       ! print *, "removed", MolecularSystem_getTotalNumberOfContractions(speciesID)-MolecularSystem_getOcupationNumber(speciesID)
       ! arguments(1) = "REMOVED-ORBITALS"
       ! call Vector_writeToFile(unit=wfnUnit, binary=.true., &
       !      value=real(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)-MolecularSystem_getOcupationNumber(speciesID,mergedMolecularSystem),8),&
       !      arguments= arguments )

    end do
    ! close(wfnUnit)
    
  end subroutine NOCIMatrices_mergeCoefficients

  !>
  !! @brief Merges the occupied orbitals coefficients from two systems
  !! @param occupationI and occupationII: Number of orbitals to merge from each matrix. 
  !! sysBasisList: array indicating which basis functions of the merged molecular system belong to sysI and sysII Merged Coefficients: Matrices for output.
  !<
  subroutine NOCIMatrices_getOffDiagonalDensityMatrix(this,sysI,sysII,mergedCoefficients,mergedMolecularSystem,overlapElement,inverseOverlapMatrices,&
       sysIbasisList,sysIIbasisList)
    type(NonOrthogonalCI), intent(inout) :: this
    integer, intent(in) :: sysI, sysII
    type(Matrix), intent(in) :: mergedCoefficients(*), inverseOverlapMatrices(*)
    type(MolecularSystem), intent(in) :: mergedMolecularSystem
    real(8), intent(in) :: overlapElement
    type(IVector), intent(in) :: sysIbasisList(*), sysIIbasisList(*)

    type(Matrix), allocatable :: mergedDensityMatrix(:)
    type(Matrix), allocatable :: exchangeCorrelationMatrices(:)
    type(Matrix) :: dftEnergyMatrix
    real(8), allocatable :: particlesInGrid(:)
   
    integer :: speciesID, otherSpeciesID, i, j, ii, jj, mu, nu
    integer :: numberOfSpecies, particlesPerOrbital, occupationNumber, numberOfContractions

    !!"Non Diagonal" terms - system pairs
    if( abs(overlapElement) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) return

    numberOfSpecies=mergedMolecularSystem%numberOfQuantumSpecies
    allocate(mergedDensityMatrix(numberOfSpecies))

    call MolecularSystem_copyConstructor(MolecularSystem_instance,mergedMolecularSystem)
    
    ! Compute density contributions
    do speciesID=1, numberOfSpecies
       particlesPerOrbital=MolecularSystem_getEta(speciesID,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(speciesID,mergedMolecularSystem)/2
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
       call Matrix_constructor(mergedDensityMatrix(speciesID),int(numberOfContractions,8),int(numberOfContractions,8),0.0_8)
       do mu = 1 , numberOfContractions
          if(sysIbasisList(speciesID)%values(mu) .eq. 0) cycle
          do nu = 1 , numberOfContractions
             if(sysIIbasisList(speciesID)%values(nu) .eq. 0) cycle
             do i = 1 , occupationNumber 
                ii= i
                do j = 1 , occupationNumber
                   jj=occupationNumber + j
                   mergedDensityMatrix(speciesID)%values(mu,nu) =  mergedDensityMatrix(speciesID)%values(mu,nu) + &
                        inverseOverlapMatrices(speciesID)%values(j,i)*&
                        mergedCoefficients(speciesID)%values(mu,ii)*&
                        mergedCoefficients(speciesID)%values(nu,jj)
                   mergedDensityMatrix(speciesID)%values(nu,mu) = mergedDensityMatrix(speciesID)%values(nu,mu) + &
                        inverseOverlapMatrices(speciesID)%values(j,i)*&
                        mergedCoefficients(speciesID)%values(mu,ii)*&
                        mergedCoefficients(speciesID)%values(nu,jj)
                end do
             end do
          end do
       end do
       mergedDensityMatrix(speciesID)%values=0.5*particlesPerOrbital*mergedDensityMatrix(speciesID)%values
       ! print *, "off diagonal matrix for", speciesID
       ! call Matrix_show(mergedDensityMatrix(speciesID))
    end do


    ! if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
    !    print *, "Superposed DFT energies:"

    !    allocate(exchangeCorrelationMatrices(numberOfSpecies), &
    !         particlesInGrid(numberOfSpecies))
    !    call DensityFunctionalTheory_buildFinalGrid()
    !    call Matrix_constructor(dftEnergyMatrix, int(numberOfSpecies,8), &
    !         int(numberOfSpecies,8), 0.0_8 )
    !    do speciesID=1, numberOfSpecies
    !       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
    !       call Matrix_constructor(exchangeCorrelationMatrices(speciesID), int(numberOfContractions,8), &
    !            int(numberOfContractions,8), 0.0_8)
    !    end do
    !    call DensityFunctionalTheory_finalDFT(mergedDensityMatrix(1:numberOfSpecies), &
    !         exchangeCorrelationMatrices, &
    !         dftEnergyMatrix, &
    !         particlesInGrid)

    !    do speciesID = 1, numberOfSpecies
    !       write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
    !            " Particles in grid = ", particlesInGrid(speciesID)
    !    end do

    !    do speciesID = 1, numberOfSpecies
    !       do otherSpeciesID = speciesID, numberOfSpecies
    !          write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
    !               "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
    !               " DFT Corr. energy = ", dftEnergyMatrix%values(speciesID,otherSpeciesID)
    !          this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,sysII)=dftEnergyMatrix%values(speciesID,otherSpeciesID)*overlapElement
    !       end do
    !    end do
    ! end if

    
    do speciesID=1, numberOfSpecies
       call Matrix_destructor(mergedDensityMatrix(speciesID))
    end do

    deallocate(mergedDensityMatrix)
    
  end subroutine NOCIMatrices_getOffDiagonalDensityMatrix


  !>
  !! @brief Computes an upper bound of the overlap between two configurations, based on the max distance between particles of the same species and the lowest exponent of the basis set functions. Assumes a localized hartree product for the heaviest species
  !!
  !! @param sysI and sysII: molecular system indices. estimatedOverlap: output value
  !<
  subroutine NOCIMatrices_prescreenOverlap(this,sysI,sysII,estimatedOverlap)
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII !Indices of the systems to screen
    real(8) :: estimatedOverlap

    type(Vector), allocatable :: displacementVector(:)    
    integer :: speciesID, k, l, m
    real(8) :: massThreshold, minExponent, speciesOverlap
    
    !displacement vectors contains the max distance between equivalent basis function centers 
    allocate(displacementVector(this%molecularSystems(sysI)%numberOfQuantumSpecies))
  
    call MolecularSystem_GetTwoSystemsDisplacement(this%molecularSystems(sysI), this%molecularSystems(sysII),displacementVector(:))

    estimatedOverlap=1.0

    !only compute for heavy particles, maybe should be a control parameter
    massThreshold=10.0   
    
    do speciesID = 1, this%molecularSystems(sysI)%numberOfQuantumSpecies
       if(this%molecularSystems(sysI)%species(speciesID)%mass .lt. massThreshold) cycle
       speciesOverlap=1.0
       !!get smallest exponent of the basis set
       do k = 1, size(this%molecularSystems(sysI)%species(speciesID)%particles)
          minExponent=1.0E8
          do l = 1, size(this%molecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction)             
             do m = 1, size(this%molecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents)
                if(this%molecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents(m).lt.minExponent) &
                     minExponent=this%molecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents(m)
                   !Assume a 1S GTF
                   ! normCoefficients(speciesID)=(2.0*minExponents(speciesID)/Math_PI)**(3.0/4.0)
             end do
          end do
          !!Compute an hipothetical overlap between two 1S functions with the lowest orbital exponent separated at the distance between systems 
          speciesOverlap=speciesOverlap*exp(-minExponent*displacementVector(speciesID)%values(k)**2/2.0)
       end do

       ! print *, "sysI", sysI, "sysII", sysII, "species", speciesID,"overlap approx", speciesOverlap
       estimatedOverlap=estimatedOverlap*speciesOverlap
    end do

    deallocate(displacementVector)
    
  end subroutine NOCIMatrices_prescreenOverlap

  !>
  !! @brief Classify the sysI and sysII pair according to their distance matrix
  !!
  !! @param sysI and sysII: molecular system indices.
  !<
  ! subroutine NOCIMatrices_classifyConfigurationPair(this,currentSysI,currentSysII,newPairFlag)
  !   implicit none
  !   type(NonOrthogonalCI) :: this
  !   integer :: currentSysI, currentSysII !Indices of the systems to classify
  !   logical :: newPairFlag
    
  !   type(MolecularSystem) :: currentMolecularSystem
  !   type(Matrix) :: currentDistanceMatrix,previousDistanceMatrix
    
  !   integer :: sysI, sysII, i, checkingType
  !   logical :: match
    
  !   call MolecularSystem_copyConstructor(currentMolecularSystem, molecularSystem_instance)
  !   newPairFlag=.true.
  !   currentDistanceMatrix=ParticleManager_getDistanceMatrix()

  !   ! print *, "Current distance matrix"
  !   ! call Matrix_show(currentDistanceMatrix)

  !   types: do checkingType=1, this%numberOfUniquePairs
  !      ! print *, "checkingType", checkingType
  !      systemI: do sysI=1, currentSysI
  !         systemII: do sysII=sysI+1, currentSysII

  !            if(sysI .eq. currentSysI .and. sysII .eq. currentSysII ) cycle types

  !            if((this%configurationPairTypes%values(sysI,sysII) .eq. checkingType) .and. &
  !               (this%systemTypes%values(sysI) .eq. this%systemTypes%values(currentSysI)) .and. & 
  !               (this%systemTypes%values(sysII) .eq. this%systemTypes%values(currentSysII))) then

  !               ! call MolecularSystem_mergeTwoSystems(molecularSystem_instance, this%MolecularSystems(sysI), this%MolecularSystems(sysII))
                
  !               previousDistanceMatrix=ParticleManager_getDistanceMatrix()

  !               ! print *, "Comparing with previous distance matrix", checkingType
  !               ! call Matrix_show(previousDistanceMatrix)          
          
  !               match=.true.
  !               do i=1, size(currentDistanceMatrix%values(:,1))
  !                  if(sum(abs(currentDistanceMatrix%values(i,:) - previousDistanceMatrix%values(i,:))) .gt. &
  !                       CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) then
  !                     match=.false.
  !                     exit
  !                  end if
  !               end do
             
  !               if(match) then
  !                  newPairFlag=.false.
  !                  this%configurationPairTypes%values(currentSysI,currentSysII)=this%configurationPairTypes%values(sysI,sysII)
  !                  exit types
  !               else
  !                  cycle types
  !               end if
  !            end if
  !         end do systemII
  !      end do systemI
  !   end do types

  !   if(newPairFlag) then
  !      this%numberOfUniquePairs=this%numberOfUniquePairs+1
  !      this%configurationPairTypes%values(currentSysI,currentSysII)=this%numberOfUniquePairs
  !   end if

  !   if(this%configurationPairTypes%values(currentSysI,currentSysII).eq.0) then
  !      print *, "newPairFlag", newPairFlag
  !      print *, currentSysI, currentSysII, this%configurationPairTypes%values(currentSysI,currentSysII)
  !      STOP "I found a type zero"
  !   end if
  !   call MolecularSystem_copyConstructor(molecularSystem_instance, currentMolecularSystem)

  ! end subroutine NOCIMatrices_classifyConfigurationPair

  
  !>
  !! @brief Computes overlap matrix element between two configurations along with one particle energy contributions
  !!
  !! @param sysI and sysII: molecular system indices. Merged Molecular System: Union of objects from sysI and sysII. Merged Coefficients: Mixed molecular system coefficients. Sys basis list indicate the basis functions of each sysI and sysII in the merged molecular system. inverseOverlapMatrices: output required for two particle contributions
  !<
  subroutine NOCIMatrices_computeOverlapAndHCoreElements(this,sysI,sysII,mergedMolecularSystem,mergedCoefficients, &
       sysIbasisList, sysIIbasisList,inverseOverlapMatrices)

    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: mergedMolecularSystem
    integer :: sysI, sysII 
    type(Matrix) :: mergedCoefficients(*), inverseOverlapMatrices(*) 
    type(IVector) :: sysIbasisList(*), sysIIbasisList(*)

    integer :: speciesID
    integer :: a,b,bb,mu,nu    
    integer :: numberOfContractions,occupationNumber,particlesPerOrbital
    type(Matrix) :: molecularOverlapMatrix
    type(Matrix), allocatable :: auxOverlapMatrix(:), auxKineticMatrix(:), auxAttractionMatrix(:), auxExternalPotMatrix(:)
    type(Matrix), allocatable :: molecularKineticMatrix(:), molecularAttractionMatrix(:), molecularExternalMatrix(:)
    type(Vector) :: overlapDeterminant
    real(8) :: oneParticleKineticEnergy,oneParticleAttractionEnergy,oneParticleExternalEnergy
    
    allocate(auxOverlapMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxKineticMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxAttractionMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxExternalPotMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         molecularKineticMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         molecularAttractionMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         molecularExternalMatrix(mergedMolecularSystem%numberOfQuantumSpecies))

    !!Initialize overlap
    this%configurationOverlapMatrix%values(sysI,sysII)=1.0
       
    call Vector_constructor(overlapDeterminant, mergedMolecularSystem%numberOfQuantumSpecies, 0.0_8)        
    
!!!!Overlap first
    do speciesID = 1, this%MolecularSystems(sysI)%numberOfQuantumSpecies

       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
       particlesPerOrbital=MolecularSystem_getEta(speciesID,mergedMolecularSystem)
       !! Calculate one- particle integrals  
       call DirectIntegralManager_getOverlapIntegrals(mergedMolecularSystem,speciesID,&
            auxOverlapMatrix(speciesID))

       !!Test 

       ! print *, "auxOverlapMatrix", speciesID
       ! call Matrix_show(auxOverlapMatrix(speciesID))

       call Matrix_constructor(molecularOverlapMatrix, int(occupationNumber,8), &
            int(occupationNumber,8), 0.0_8 )

       do mu=1, numberOfContractions!sysI
          if(sysIbasisList(speciesID)%values(mu) .eq. 0 ) cycle
          do nu=1, numberOfContractions  !sysII
             if(sysIIbasisList(speciesID)%values(nu) .eq. 0) cycle
             do a=1, occupationNumber !sysI
                do b=occupationNumber+1, 2*occupationNumber
                     bb=b-occupationNumber
                     ! print *, "a, b, mu, nu, coefI, coefII", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b),auxOverlapMatrix(speciesID)%values(mu,nu)
                   
                     molecularOverlapMatrix%values(a,bb)=molecularOverlapMatrix%values(a,bb)+&
                          mergedCoefficients(speciesID)%values(mu,a)*&
                          mergedCoefficients(speciesID)%values(nu,b)*&
                          auxOverlapMatrix(speciesID)%values(mu,nu)
                end do
             end do
          end do
       end do

       ! print *, "molecularOverlapMatrix sysI, sysII, speciesID", sysI, sysII, speciesID
       ! call Matrix_show(molecularOverlapMatrix)

       !Sometimes we run calculations for systems with ghost species
       if(occupationNumber .ne. 0) then
          inverseOverlapMatrices(speciesID)=Matrix_inverse(molecularOverlapMatrix)
          ! print *, "inverseOverlapMatrices sysI, sysII", speciesID, sysI, sysII
          ! call Matrix_show(inverseOverlapMatrices(speciesID))
          call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant%values(speciesID),method="LU")             
          ! print *, "OverlapDeterminantLU speciesID, sysI, sysII", speciesID, sysI, sysII, overlapDeterminant%values(speciesID)
       else
          overlapDeterminant%values(speciesID)=1.0
       end if
       
       this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*overlapDeterminant%values(speciesID)**particlesPerOrbital


    end do

    ! print *, "total overlap", this%configurationOverlapMatrix%values(sysI,sysII)
    
    !!Skip the rest of the evaluation if the overlap is smaller than the threshold
    if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
         abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) return

    !!Point charge-Point charge repulsion
    this%configurationHamiltonianMatrix%values(sysI,sysII)=MolecularSystem_getPointChargesEnergy(this%molecularSystems(sysI))*&
         this%configurationOverlapMatrix%values(sysI,sysII)
    ! print *, "Point charge-Point charge repulsion", MolecularSystem_getPointChargesEnergy()

    !!Compute hcore if overlap is significant
    do speciesID = 1, this%molecularSystems(sysI)%numberOfQuantumSpecies

       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
       particlesPerOrbital=MolecularSystem_getEta(speciesID,mergedMolecularSystem)
       
       call Matrix_constructor(auxKineticMatrix(speciesID),&
            int(numberOfContractions,8),int(numberOfContractions,8),0.0_8)
       call Matrix_constructor(auxAttractionMatrix(speciesID),&
            int(numberOfContractions,8),int(numberOfContractions,8),0.0_8)
       call Matrix_constructor(auxExternalPotMatrix(speciesID),&
            int(numberOfContractions,8),int(numberOfContractions,8),0.0_8)

       call DirectIntegralManager_getKineticIntegrals(mergedMolecularSystem,speciesID,auxKineticMatrix(speciesID))
       call DirectIntegralManager_getAttractionIntegrals(mergedMolecularSystem,speciesID,auxAttractionMatrix(speciesID))
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
            call DirectIntegralManager_getExternalPotentialIntegrals(mergedMolecularSystem,speciesID,auxExternalPotMatrix(speciesID))
       
       !! Incluiding mass effect       
       if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
          auxKineticMatrix(speciesID)%values =  &
            auxKineticMatrix(speciesID)%values * &
            ( 1.0_8/MolecularSystem_getMass( speciesID,this%molecularSystems(sysI) ) -1.0_8 / MolecularSystem_getTotalMass(this%molecularSystems(sysI)) )
       else
          auxKineticMatrix(speciesID)%values =  &
            auxKineticMatrix(speciesID)%values / &
            MolecularSystem_getMass( speciesID,this%molecularSystems(sysI) )
       end if

       !! Including charge
       auxAttractionMatrix(speciesID)%values=auxAttractionMatrix(speciesID)%values*(-MolecularSystem_getCharge(speciesID,this%molecularSystems(sysI)))                         
       
       call Matrix_constructor(molecularKineticMatrix(speciesID), int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
       call Matrix_constructor(molecularAttractionMatrix(speciesID), int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
       call Matrix_constructor(molecularExternalMatrix(speciesID), int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )

       !!Test 
       ! print *, "auxKineticMatrix", speciesID
       ! call Matrix_show(auxKineticMatrix(speciesID))
       ! print *, "auxAttractionMatrix", speciesID
       ! call Matrix_show(auxAttractionMatrix(speciesID))

       do mu=1, numberOfContractions !sysI
          if(sysIbasisList(speciesID)%values(mu) .eq. 0) cycle
          do nu=1, numberOfContractions !sysII
             if(sysIIbasisList(speciesID)%values(nu) .eq. 0) cycle
             do a=1, occupationNumber !sysI
                do b=occupationNumber+1, 2*occupationNumber
                   bb=b-occupationNumber

                   ! print *, "hcore", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
                   !      auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
                   !      auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID))

                   molecularKineticMatrix(speciesID)%values(a,bb)=molecularKineticMatrix(speciesID)%values(a,bb)+&
                        mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                        auxKineticMatrix(speciesID)%values(mu,nu)                        

                   molecularAttractionMatrix(speciesID)%values(a,bb)=molecularAttractionMatrix(speciesID)%values(a,bb)+&
                        mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                        auxAttractionMatrix(speciesID)%values(mu,nu)                         
                   
                   if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                        molecularExternalMatrix(speciesID)%values(a,bb)=molecularExternalMatrix(speciesID)%values(a,bb)+&
                        mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                        auxExternalPotMatrix(speciesID)%values(mu,nu)                         
                end do
             end do
          end do
       end do
       molecularKineticMatrix(speciesID)%values=particlesPerOrbital*molecularKineticMatrix(speciesID)%values
       molecularAttractionMatrix(speciesID)%values=particlesPerOrbital*molecularAttractionMatrix(speciesID)%values
       molecularExternalMatrix(speciesID)%values=particlesPerOrbital*molecularExternalMatrix(speciesID)%values
       !!End test                          
    end do

    !!One Particle Terms
    do speciesID=1, this%molecularSystems(sysI)%numberOfQuantumSpecies
       oneParticleKineticEnergy=0.0
       oneParticleAttractionEnergy=0.0
       oneParticleExternalEnergy=0.0
       occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
       do a=1, occupationNumber !sysI
          do b=1, occupationNumber !sysII
             oneParticleKineticEnergy=oneParticleKineticEnergy+ molecularKineticMatrix(speciesID)%values(a,b)*&
                  inverseOverlapMatrices(speciesID)%values(b,a)
             oneParticleAttractionEnergy=oneParticleAttractionEnergy+ molecularAttractionMatrix(speciesID)%values(a,b)*&
                  inverseOverlapMatrices(speciesID)%values(b,a)
             oneParticleExternalEnergy=oneParticleExternalEnergy+ molecularExternalMatrix(speciesID)%values(a,b)*&
                  inverseOverlapMatrices(speciesID)%values(b,a)
          end do          
       end do
       this%configurationKineticMatrix(speciesID)%values(sysI,sysII)=oneParticleKineticEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
       this%configurationPuntualMatrix(speciesID)%values(sysI,sysII)=oneParticleAttractionEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
       this%configurationExternalMatrix(speciesID)%values(sysI,sysII)=oneParticleExternalEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
       
       this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+&
            (oneParticleKineticEnergy+oneParticleAttractionEnergy+oneParticleExternalEnergy)*this%configurationOverlapMatrix%values(sysI,sysII)
       ! print *, "sysI, sysII", sysI, sysII, "oneParticleEnergy for species", speciesID, oneParticleEnergy
    end do

    deallocate(auxOverlapMatrix, auxKineticMatrix, auxAttractionMatrix, auxExternalPotMatrix, &
         molecularKineticMatrix, molecularAttractionMatrix, molecularExternalMatrix)
    
  end subroutine NOCIMatrices_computeOverlapAndHCoreElements
  !>
  !! @brief Computes the two particles contributions to the non diagonal elements of the hamiltonian matrix
  !!
  !! @param this, sysI,sysII: system indexes, inverseOverlapMatrices, mergedCoefficients are required to evaluate the elements
  !<
  subroutine NOCIMatrices_twoParticlesContributions(this,sysI,sysII,mergedMolecularSystem,inverseOverlapMatrices,mergedCoefficients,Libint2LocalInstance)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII 
    type(MolecularSystem) :: mergedMolecularSystem
    type(Matrix) :: inverseOverlapMatrices(*) 
    type(Matrix) :: mergedCoefficients(*) 
    type(Libint2Interface) :: Libint2LocalInstance(*)
   
    type(matrix), allocatable :: fourCenterIntegrals(:,:)
    type(imatrix), allocatable :: twoIndexArray(:),fourIndexArray(:)
    integer :: numberOfContractions,occupationNumber,particlesPerOrbital
    integer :: otherNumberOfContractions,otherOccupationNumber,otherParticlesPerOrbital
    integer :: ssize1, auxIndex, auxIndex1
    integer :: a,b,bb,c,d,dd,i,j
    real(8) :: hartreeEnergy, exchangeEnergy 

    allocate(fourCenterIntegrals(mergedMolecularSystem%numberOfQuantumSpecies,mergedMolecularSystem%numberOfQuantumSpecies), &
         twoIndexArray(mergedMolecularSystem%numberOfQuantumSpecies), &
         fourIndexArray(mergedMolecularSystem%numberOfQuantumSpecies))

    !!Fill indexes arrays
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies       
       ! print *, "reading integrals species", i
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(i,mergedMolecularSystem)
       !!Two particle integrals indexes
       call Matrix_constructorInteger(twoIndexArray(i),  &
            int(max(numberOfContractions,occupationNumber),8), &
            int(max(numberOfContractions,occupationNumber),8), 0 )

       c = 0
       do a=1,max(numberOfContractions,occupationNumber)
          do b=a, max(numberOfContractions,occupationNumber)
             c = c + 1
             twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
             twoIndexArray(i)%values(b,a) = twoIndexArray(i)%values(a,b)
          end do
       end do

       ssize1 = max(numberOfContractions,occupationNumber)
       ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2

       call Matrix_constructorInteger(fourIndexArray(i), int( ssize1,8), int( ssize1,8) , 0 )

       c = 0
       do a = 1, ssize1
          do b = a, ssize1
             c = c + 1
             fourIndexArray(i)%values(a,b) = c! IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
             fourIndexArray(i)%values(b,a) = fourIndexArray(i)%values(a,b)
          end do
       end do
    end do
    
    !! Calculate two- particle integrals  
    call NOCIMatrices_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, &
         twoIndexArray, fourIndexArray, fourCenterIntegrals, Libint2LocalInstance)

!!!Add charges
    if ( .not. InterPotential_instance%isInstanced) then
       do i=1, mergedMolecularSystem%numberOfQuantumSpecies
          fourCenterIntegrals(i,i)%values = &
               fourCenterIntegrals(i,i)%values * mergedMolecularSystem%species(i)%charge**2.0

          do j = i+1 , mergedMolecularSystem%numberOfQuantumSpecies
             fourCenterIntegrals(i,j)%values = &
                  fourCenterIntegrals(i,j)%values * mergedMolecularSystem%species(i)%charge * mergedMolecularSystem%species(j)%charge
          end do
       end do
    end if
!!!Compute Hamiltonian Matrix element between displaced geometries

    ! !!Point charge-Point charge repulsion
    ! !!One Particle Terms
    ! !!Have already been computed

    !!Same species repulsion
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(i,this%molecularSystems(sysI))
       particlesPerOrbital=MolecularSystem_getEta(i,mergedMolecularSystem)
       hartreeEnergy=0.0
       exchangeEnergy=0.0
       do a=1,occupationNumber !sysI
          do b=occupationNumber+1, 2*occupationNumber !sysII
             bb=b-occupationNumber
             do c=1, occupationNumber !sysI
                do d=occupationNumber+1, 2*occupationNumber !sysII
                   dd=d-occupationNumber
                   auxIndex = fourIndexArray(i)%values(twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d) )
                   hartreeEnergy=hartreeEnergy+0.5*fourCenterIntegrals(i,i)%values(auxIndex, 1)*&
                        inverseOverlapMatrices(i)%values(bb,a)*inverseOverlapMatrices(i)%values(dd,c)*particlesPerOrbital**2 !coulomb
                   exchangeEnergy=exchangeEnergy-0.5*fourCenterIntegrals(i,i)%values(auxIndex, 1)*&
                        inverseOverlapMatrices(i)%values(dd,a)*inverseOverlapMatrices(i)%values(bb,c)*particlesPerOrbital !exchange
                   ! print *, a, b, c, d, twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d), fourIndexArray(i)%values( &
                   !      twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d)), 
                end do
             end do
          end do
       end do
       this%configurationHartreeMatrix(i,i)%values(sysI,sysII)=hartreeEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
       this%configurationExchangeMatrix(i)%values(sysI,sysII)=exchangeEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
       this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+&
            (hartreeEnergy+exchangeEnergy)*this%configurationOverlapMatrix%values(sysI,sysII)
       ! print *, "same species interactionEnergy for species", i, hartreeEnergy, exchangeEnergy
    end do

    !!Interspecies repulsion
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies-1
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem)
       occupationNumber=MolecularSystem_getOcupationNumber(i,this%molecularSystems(sysI))
       particlesPerOrbital=MolecularSystem_getEta(i,mergedMolecularSystem)
       do j=i+1, mergedMolecularSystem%numberOfQuantumSpecies
          otherNumberOfContractions=MolecularSystem_getTotalNumberOfContractions(j,mergedMolecularSystem)
          otherOccupationNumber=MolecularSystem_getOcupationNumber(j,mergedMolecularSystem)
          otherParticlesPerOrbital=MolecularSystem_getEta(j,mergedMolecularSystem)
          hartreeEnergy=0.0
          ssize1 = max(otherNumberOfContractions,otherOccupationNumber)
          ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2
          otherOccupationNumber=MolecularSystem_getOcupationNumber(j,this%molecularSystems(sysI))
          do a=1, occupationNumber !sysI
             do b=occupationNumber+1, 2*occupationNumber !sysII
                bb=b-MolecularSystem_getOcupationNumber(i,this%molecularSystems(sysI))
                auxIndex1 = ssize1 * (twoIndexArray(i)%values(a,b) - 1 ) 
                do c=1, otherOccupationNumber  !sysI
                   do d=otherOccupationNumber+1,2*otherOccupationNumber !sysII
                      dd=d-otherOccupationNumber
                      auxIndex = auxIndex1  + twoIndexArray(j)%values(c,d) 
                      hartreeEnergy=hartreeEnergy+fourCenterIntegrals(i,j)%values(auxIndex, 1)*&
                           inverseOverlapMatrices(i)%values(bb,a)*inverseOverlapMatrices(j)%values(dd,c)*&
                           particlesPerOrbital*otherParticlesPerOrbital
                           ! print *, a, b, c, d,  fourCenterIntegrals(i,j)%values(auxIndex, 1), inverseOverlapMatrices(i)%values(bb,a), inverseOverlapMatrices(j)%values(dd,c)
                   end do
                end do
             end do
          end do
          this%configurationHartreeMatrix(i,j)%values(sysI,sysII)=hartreeEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
          this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+&
               hartreeEnergy*this%configurationOverlapMatrix%values(sysI,sysII)
          ! print *, "interspecies hartreeEnergy for species", i, j, hartreeEnergy
       end do
    end do

    deallocate(fourCenterIntegrals,twoIndexArray,fourIndexArray)

  end subroutine NOCIMatrices_twoParticlesContributions

  !>
  !! @brief Solves the NOCI matrix equation
  !!
  !! @param this
  !<
  subroutine NOCIMatrices_diagonalize(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(Matrix) :: transformationMatrix,transformedHamiltonianMatrix,eigenVectors,auxMatrix
    type(Vector) :: eigenValues
    integer :: removedStates
    integer :: speciesID,otherSpeciesID,sysI,sysII,state,i,j
    real(8) :: auxEnergy
    real(8) :: timeA
    
    !$  timeA = omp_get_wtime()

    call Matrix_constructor(this%configurationCoefficients, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Vector_constructor(this%statesEigenvalues, this%numberOfDisplacedSystems, 0.0_8)

    ! print *, "non orthogonal CI overlap Matrix "
    ! call Matrix_show(this%configurationOverlapMatrix)

    ! print *, "non orthogonal CI Hamiltionian Matrix "
    ! call Matrix_show(this%configurationHamiltonianMatrix)
    ! 
    print *, ""
    print *, "Transforming non orthogonal CI Hamiltonian Matrix..."
    
    call Matrix_constructor(transformationMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8) , 0.0_8)
    call Matrix_constructor(transformedHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8) , 0.0_8)

    call Vector_constructor( eigenValues, this%numberOfDisplacedSystems )
    call Matrix_constructor( eigenVectors,int(this%numberOfDisplacedSystems,8),int(this%numberOfDisplacedSystems,8))

    !!****************************************************************
    !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
    !!          
    call Matrix_eigen( this%configurationOverlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )

    ! print *,"Overlap eigenvectors "
    ! call Matrix_show( eigenVectors )

    ! print *,"Overlap eigenvalues "
    ! call Vector_show( eigenValues )
    
    !! Remove states from configurations with linear dependencies
    do i = 1 , this%numberOfDisplacedSystems
       do j = 1 , this%numberOfDisplacedSystems
          if ( abs(eigenValues%values(j)) >= CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
             transformationMatrix%values(i,j) = &
                  eigenVectors%values(i,j)/sqrt( eigenvalues%values(j) )
          else
             transformationMatrix%values(i,j) = 0
          end if
       end do
    end do

    removedStates=0
    do i = 1 , this%numberOfDisplacedSystems
       if ( abs(eigenValues%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) &
            removedStates=removedStates+1
    end do

    if (removedStates .gt. 0) &
         write(*,"(A,I5,A,ES9.3)") "Removed ", removedStates , &
         " states from the CI transformation Matrix with overlap eigen threshold of ", CONTROL_instance%OVERLAP_EIGEN_THRESHOLD


    !!Ortogonalizacion simetrica
    transformationMatrix%values  = &
         matmul(transformationMatrix%values, transpose(eigenVectors%values))

    ! print *,"Matriz de transformacion "
    ! call Matrix_show( transformationMatrix )

    !!**********************************************************************************************
    !! Transform configuration hamiltonian matrix
    !!
    transformedHamiltonianMatrix%values = &
         matmul( matmul( transpose( transformationMatrix%values ) , &
         this%configurationHamiltonianMatrix%values), transformationMatrix%values )

    ! print *,"transformed Hamiltonian Matrix "
    ! call Matrix_show( this%configurationHamiltonianMatrix )

    print *, "Diagonalizing non orthogonal CI Hamiltonian Matrix..."
    !! Calcula valores y vectores propios de matriz de CI transformada.
    call Matrix_eigen( transformedHamiltonianMatrix, this%statesEigenvalues, this%configurationCoefficients, SYMMETRIC )

    !! Calcula los  vectores propios para matriz de CI       
    this%configurationCoefficients%values = matmul( transformationMatrix%values, this%configurationCoefficients%values )

    ! print *,"non orthogonal CI eigenvalues "
    ! call Vector_show( this%statesEigenvalues )

    ! print *,"configuration Coefficients"
    ! call Matrix_show( this%configurationCoefficients )

    write(*,"(A)") ""
    write(*,"(A)") " MIXED HARTREE-FOCK CALCULATION"
    write(*,"(A)") " NON ORTHOGONAL CONFIGURATION INTERACTION"
    write(*,"(A)") " EIGENVALUES AND EIGENVECTORS: "
    write(*,"(A)") "========================================="
    write(*,"(A)") ""
    do state = 1, min(CONTROL_instance%NUMBER_OF_CI_STATES,this%numberOfDisplacedSystems)
       write (*,"(A)")  ""
       write (*,"(T9,A17,I3,A10, F25.12)") "STATE: ", state, " ENERGY = ", this%statesEigenvalues%values(state)
       write (*,"(A38)")  "Components: "
       write(*,"(A38,F25.12)") " Point charges energy = ", MolecularSystem_getPointChargesEnergy(this%molecularSystems(1))
       do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
          auxEnergy=0
          do sysI=1, this%numberOfDisplacedSystems
             auxEnergy= auxEnergy+ &
                  this%configurationCoefficients%values(sysI,state)**2*&
                  this%configurationKineticMatrix(speciesID)%values(sysI,sysI)
             do sysII=sysI+1, this%numberOfDisplacedSystems
                auxEnergy= auxEnergy+ &
                     2.0_8*this%configurationCoefficients%values(sysI,state)*&
                     this%configurationCoefficients%values(sysII,state)*&
                     this%configurationKineticMatrix(speciesID)%values(sysI,sysII)
             end do
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Kinetic energy = ", auxEnergy
       end do
       do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
          auxEnergy=0
          do sysI=1, this%numberOfDisplacedSystems
             auxEnergy= auxEnergy+ &
                  this%configurationCoefficients%values(sysI,state)**2*&
                  this%configurationPuntualMatrix(speciesID)%values(sysI,sysI)
             do sysII=sysI+1, this%numberOfDisplacedSystems
                auxEnergy= auxEnergy+ &
                     2.0_8*this%configurationCoefficients%values(sysI,state)*&
                     this%configurationCoefficients%values(sysII,state)*&
                     this%configurationPuntualMatrix(speciesID)%values(sysI,sysII)
             end do
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Puntual energy = ", auxEnergy
       end do
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
             auxEnergy=0
             do sysI=1, this%numberOfDisplacedSystems
                auxEnergy= auxEnergy+ &
                     this%configurationCoefficients%values(sysI,state)**2*&
                     this%configurationExternalMatrix(speciesID)%values(sysI,sysI)
                do sysII=sysI+1, this%numberOfDisplacedSystems
                   auxEnergy= auxEnergy+ &
                        2.0_8*this%configurationCoefficients%values(sysI,state)*&
                        this%configurationCoefficients%values(sysII,state)*&
                        this%configurationExternalMatrix(speciesID)%values(sysI,sysII)
                end do
             end do
             write(*,"(A38,F25.12)") trim(this%molecularSystems(1)%species(speciesID)%name ) // &
                  " External energy = ", auxEnergy
          end do
       end if
       do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
          auxEnergy=0
          do sysI=1, this%numberOfDisplacedSystems
             auxEnergy= auxEnergy+ &
                  this%configurationCoefficients%values(sysI,state)**2*&
                  this%configurationExchangeMatrix(speciesID)%values(sysI,sysI)
             do sysII=sysI+1, this%numberOfDisplacedSystems
                auxEnergy= auxEnergy+ &
                     2.0_8*this%configurationCoefficients%values(sysI,state)*&
                     this%configurationCoefficients%values(sysII,state)*&
                     this%configurationExchangeMatrix(speciesID)%values(sysI,sysII)
             end do
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                  " Exchange energy = ", auxEnergy

          do otherSpeciesID=speciesID, this%molecularSystems(1)%numberOfQuantumSpecies
             auxEnergy=0
             do sysI=1, this%numberOfDisplacedSystems
                auxEnergy= auxEnergy+ &
                     this%configurationCoefficients%values(sysI,state)**2*&
                     this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,sysI)
                do sysII=sysI+1, this%numberOfDisplacedSystems
                   auxEnergy= auxEnergy+ &
                        2.0_8*this%configurationCoefficients%values(sysI,state)*&
                        this%configurationCoefficients%values(sysII,state)*&
                        this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,sysII)
                end do
             end do
             write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                  "/"//trim( this%molecularSystems(1)%species(otherSpeciesID)%name ) // &
                  " Hartree energy = ", auxEnergy
          end do
       end do
    end do
    write(*,"(A)") ""

    call Matrix_constructor(auxMatrix,int(this%numberOfDisplacedSystems,8),&
         int(CONTROL_instance%CI_STATES_TO_PRINT,8),0.0_8)
    do i=1, this%numberOfDisplacedSystems
       do j=1, CONTROL_instance%CI_STATES_TO_PRINT
          auxMatrix%values(i,j)=this%configurationCoefficients%values(i,j)
       end do
    end do

    
    write(*,"(I5,A)") CONTROL_instance%CI_STATES_TO_PRINT, " LOWEST LYING STATES CONFIGURATION COEFFICIENTS"
    write(*,"(A)") ""
    call Matrix_show(auxMatrix , &
         rowkeys = this%systemLabels, &
         columnkeys = string_convertvectorofrealstostring( this%statesEigenvalues ),&
         flags=WITH_BOTH_KEYS)
    write(*,"(A)") ""

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI matrix diagonalization : ", omp_get_wtime() - timeA ," (s)"
    
  end subroutine NOCIMatrices_diagonalize

  !>
  !! @brief Calculate and Transform the four center integrals in one sweep without writing anything to disk
  !!
  !! @param molecularSystem, HFCoefficients: species array with the atomic coefficients, fourCenterIntegrals: species*species array to save integrals
  !<
  subroutine NOCIMatrices_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, twoIndexArray, fourIndexArray, fourCenterIntegrals, Libint2LocalInstance)
    implicit none
    type(MolecularSystem), intent(in) :: mergedMolecularSystem
    type(Matrix), intent(in) :: mergedCoefficients(mergedMolecularSystem%numberOfQuantumSpecies)
    type(iMatrix), intent(in) :: twoIndexArray(mergedMolecularSystem%numberOfQuantumSpecies)
    type(iMatrix), intent(in) :: fourIndexArray(mergedMolecularSystem%numberOfQuantumSpecies)
    type(Matrix), intent(out) :: fourCenterIntegrals(mergedMolecularSystem%numberOfQuantumSpecies,mergedMolecularSystem%numberOfQuantumSpecies)
    type(Libint2Interface) :: Libint2LocalInstance(mergedMolecularSystem%numberOfQuantumSpecies)

    real(8), allocatable, target :: ints(:,:,:,:)
    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer :: p, p_l, p_u
    integer :: q, q_l, q_u
    integer :: r, r_l, r_u
    integer :: s, s_l, s_u
    integer :: ssize, ssizeb, auxIndex, auxIndexA
    integer :: n,u, mu,nu, lambda,sigma
    real(8) :: auxTransformedTwoParticlesIntegral
    
    type(Matrix) :: densityMatrix
    integer :: speciesID, otherSpeciesID
    integer :: numberOfOrbitals, otherNumberOfOrbitals
    integer(8) :: numberOfIntegrals
    
    do speciesID=1, mergedMolecularSystem%numberOfQuantumSpecies
       numberOfOrbitals = max( MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem), &
            MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem ))
       numberOfIntegrals= int( ( (  numberOfOrbitals * ( numberOfOrbitals + 1.0_8 ) / 4.0_8 ) * &
            ( (  numberOfOrbitals * (  numberOfOrbitals + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )

       call Matrix_constructor( fourCenterIntegrals(speciesID,speciesID), numberOfIntegrals, 1_8, 0.0_8 )

       p_l = 1
       p_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2 
       q_l = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2+1
       q_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )

       r_l = 1
       r_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2 
       s_l = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2+1
       s_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )
       
       ssize = MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
       call Matrix_constructor( densityMatrix, int(ssize,8), int(ssize,8), 1.0_8 ) !Test filling with values later
       
       ! Prepare matrix
       if(allocated(ints)) deallocate(ints)
       if(allocated(tempA)) deallocate (tempA)
       if(allocated(tempB)) deallocate (tempB)
       if(allocated(tempC)) deallocate (tempC)      
       allocate (ints ( ssize, ssize, ssize, ssize ), &
            tempA ( ssize, ssize, ssize ), &
            tempB ( ssize, ssize ), &
            tempC ( ssize ))
       ints = 0

       call DirectIntegralManager_getDirectIntraRepulsionIntegralsAll(&
            speciesID, &
            densityMatrix, & 
            ints, mergedMolecularSystem, Libint2LocalInstance(speciesID) )

       
       do p = p_l, p_u
          tempA = 0
          n = p

       !    !First quarter transformation happens here
          do mu = 1, ssize
             !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
             tempA(:,:,:) = tempA(:,:,:) + mergedCoefficients(speciesID)%values( mu, p )* &
                  ints(:,:,:,mu)
          end do
          
          do q = p, q_u
             u = q
             tempB = 0

             if ( q < q_l ) cycle
             !! second quarter
             do nu = 1, ssize
                !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
                tempB(:,:) = tempB(:,:) + mergedCoefficients(speciesID)%values( nu, q )* &
                     tempA(:,:,nu)
             end do

             do r = n, r_u

                tempC = 0

                !Why??
                !if ( r <  this%r_l  ) cycle

                !! third quarter
                do lambda = 1, ssize
                   !if ( abs(coefficientsOfAtomicOrbitals%values( lambda, r )) < 1E-10 ) cycle
                   tempC(:) = tempC(:) + mergedCoefficients(speciesID)%values( lambda, r )* &
                        tempB(:,lambda)
                end do

                do s = u, s_u
                   auxTransformedTwoParticlesIntegral = 0

                   if ( s < s_l ) cycle
                   !! fourth quarter
                   do sigma = 1, ssize
                      auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                           mergedCoefficients(speciesID)%values( sigma, s )* &
                           tempC(sigma)

                   end do
                   auxIndex = fourIndexArray(speciesID)%values(twoIndexArray(speciesID)%values(p,q), twoIndexArray(speciesID)%values(r,s) )
                   fourCenterIntegrals(speciesID,speciesID)%values(auxIndex, 1) = auxTransformedTwoParticlesIntegral
                   ! print *, speciesID, p, q, r, s, auxIndex, auxTransformedTwoParticlesIntegral
                end do
                u = r + 1
             end do
          end do
       end do       
    end do
    
    do speciesID=1, mergedMolecularSystem%numberOfQuantumSpecies-1
       do otherSpeciesID=speciesID+1, mergedMolecularSystem%numberOfQuantumSpecies

          numberOfOrbitals = max( MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem), &
               MolecularSystem_getOcupationNumber(speciesID,mergedMolecularSystem))
          otherNumberOfOrbitals = max( MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,mergedMolecularSystem), &
               MolecularSystem_getOcupationNumber(otherSpeciesID,mergedMolecularSystem))

          numberOfIntegrals = int((numberOfOrbitals*((numberOfOrbitals+1.0_8)/2.0_8)) * &
                              (otherNumberOfOrbitals*(otherNumberOfOrbitals+1.0_8)/2.0_8),8)

          call Matrix_constructor( fourCenterIntegrals(speciesID,otherSpeciesID), numberOfIntegrals, 1_8, 0.0_8 )

          ssize = MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
          ssizeb = MolecularSystem_getTotalNumberOfContractions(otherSpeciesID,mergedMolecularSystem)

          call Matrix_constructor( densityMatrix, int(ssize,8), int(ssize,8), 1.0_8 ) !Test filling with values later

          p_l = 1
          p_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2 
          q_l = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )/2+1
          q_u = MolecularSystem_getOcupationNumber( speciesID, mergedMolecularSystem )

          r_l = 1
          r_u = MolecularSystem_getOcupationNumber( otherSpeciesID, mergedMolecularSystem )/2 
          s_l = MolecularSystem_getOcupationNumber( otherSpeciesID, mergedMolecularSystem )/2+1
          s_u = MolecularSystem_getOcupationNumber( otherSpeciesID, mergedMolecularSystem )
          
          ! Prepare matrix
       ! Prepare matrix
          if(allocated(ints)) deallocate(ints)
          if(allocated(tempA)) deallocate (tempA)
          if(allocated(tempB)) deallocate (tempB)
          if(allocated(tempC)) deallocate (tempC)      
          allocate (ints ( ssizeb, ssizeb, ssize, ssize ), &
               tempA ( ssizeb, ssizeb, ssize ), &
               tempB ( ssizeb, ssizeb ), &
               tempC ( ssizeb ))
          ints = 0

          call DirectIntegralManager_getDirectInterRepulsionIntegralsAll(&
               speciesID, otherSpeciesID, &
               densityMatrix, & 
               ints, mergedMolecularSystem, Libint2LocalInstance(speciesID), Libint2LocalInstance(otherSpeciesID) )

          ! do mu = 1, ssize
          !    do nu = 1, ssize
          !       do lambda = 1, ssizeb
          !          do sigma = 1, ssizeb
          !             print *, mu, nu, lambda, sigma, ints(lambda,sigma,nu,mu)
          !          end do
          !       end do
          !    end do
          ! end do
          do p = p_l, p_u
             tempA = 0
             !First quarter transformation happens here
             do mu = 1, ssize
                !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
                tempA(:,:,:) = tempA(:,:,:) + mergedCoefficients(speciesID)%values( mu, p )* &
                     ints(:,:,:,mu)
             end do

             do q = q_l, q_u
                tempB = 0

                ! if ( q < p ) cycle
                !! second quarter
                do nu = 1, ssize
                   !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle

                   tempB(:,:) = tempB(:,:) + mergedCoefficients(speciesID)%values( nu, q )* &
                        tempA(:,:,nu)
                end do

                auxIndexA = (otherNumberOfOrbitals*(otherNumberOfOrbitals+1))/2 * (twoIndexArray(speciesID)%values(p,q) - 1 ) 

                do r = r_l , r_u

                   tempC = 0

                   !! third quarter
                   do lambda = 1, ssizeb

                      tempC(:) = tempC(:) + mergedCoefficients(otherSpeciesID)%values( lambda, r )* &
                           tempB(:,lambda)

                   end do
                   do s = s_l, s_u
                      auxTransformedTwoParticlesIntegral = 0

                      ! if ( s < r ) cycle
                      !! fourth quarter
                      do sigma = 1, ssizeb
                         auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                              mergedCoefficients(otherSpeciesID)%values( sigma, s )* &
                              tempC(sigma)

                      end do

                      auxIndex = auxIndexA  + twoIndexArray(otherSpeciesID)%values(r,s) 

                      fourCenterIntegrals(speciesID,otherSpeciesID)%values(auxIndex, 1) = auxTransformedTwoParticlesIntegral

                      ! print *, speciesID,otherSpeciesID, p, q, r, s, auxIndex, auxTransformedTwoParticlesIntegral

                   end do
                end do
             end do
          end do
          
       end do
    end do

    ! call DirectIntegralManager_destructor(Libint2LocalInstance)
    
  end subroutine NOCIMatrices_transformIntegralsMemory

end module NOCIMatrices_

