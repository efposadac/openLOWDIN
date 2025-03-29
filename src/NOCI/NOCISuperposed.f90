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

module NOCISuperposed_
  use NOCIBuild_
  use NOCIMatrices_
  use MolecularSystem_
  use ParticleManager_
  use Matrix_
  use Vector_
  use DirectIntegralManager_
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
       NOCISuperposed_generateSuperposedSystem,&
       NOCISuperposed_buildDensityMatrix,&
       NOCISuperposed_getNaturalOrbitals,&
       NOCISuperposed_saveToFile

  private

contains
  
  !>
  !! @brief Generates one molecular system combining all the displaced geometries and coefficients
  !!
  !! @param this
  !<
  subroutine NOCISuperposed_generateSuperposedSystem(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: auxMolecularSystem
    type(Matrix), allocatable :: auxCoefficients(:)
    type(IVector), allocatable :: auxBasisList(:)
    
    integer :: i, sysI, speciesID
    integer :: numberOfSpecies
    
    real(8) :: timeA

    !$  timeA = omp_get_wtime()
    
    if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    numberOfSpecies=this%molecularSystems(1)%numberOfQuantumSpecies
    
    allocate(this%sysBasisList(this%numberOfDisplacedSystems,numberOfSpecies),&
         auxCoefficients(numberOfSpecies),&         
         auxBasisList(numberOfSpecies))
    
    !Create a super molecular system
    !!!Merge coefficients from system 1 and system 2
    call MolecularSystem_mergeTwoSystems(this%mergedMolecularSystem, this%molecularSystems(1), this%molecularSystems(2), &
         this%sysBasisList(1,:),this%sysBasisList(2,:))
    
    call NOCIMatrices_mergeCoefficients(numberOfSpecies,this%HFCoefficients(1,:),this%HFCoefficients(2,:),&
         this%molecularSystems(1),this%molecularSystems(2),this%mergedMolecularSystem,&
         this%sysBasisList(1,:),this%sysBasisList(2,:),this%mergedCoefficients(:))

    ! do speciesID=1, numberOfSpecies
       ! print *, "2", speciesID, "ocupationNumber", MolecularSystem_getOcupationNumber(speciesID,this%mergedMolecularSystem)
       ! print *, "2", speciesID, "mergedCoefficients"
       ! call Matrix_show(this%mergedCoefficients(speciesID))
    ! end do
    ! 
    !! Loop other systems expanding the merged coefficients matrix 
    do sysI=3, this%numberOfDisplacedSystems
       call MolecularSystem_copyConstructor(auxMolecularSystem,this%mergedMolecularSystem)
       do speciesID=1, numberOfSpecies
          call Matrix_copyConstructor(auxCoefficients(speciesID), this%mergedCoefficients(speciesID))       
       end do
       call MolecularSystem_mergeTwoSystems(this%mergedMolecularSystem, auxMolecularSystem, this%molecularSystems(sysI), &
            auxBasisList,this%sysBasisList(sysI,:),reorder=.false.)
       call NOCIMatrices_mergeCoefficients(numberOfSpecies,auxCoefficients,this%HFCoefficients(sysI,:),&
            auxMolecularSystem,this%molecularSystems(sysI),this%mergedMolecularSystem,&
            auxBasisList,this%sysBasisList(sysI,:),this%mergedCoefficients(:))
       ! do speciesID=1, numberOfSpecies
             ! print *, sysI, speciesID, "ocupationNumber", MolecularSystem_getOcupationNumber(speciesID,this%mergedMolecularSystem)
             ! print *, sysI, speciesID, "mergedCoefficients"
             ! call Matrix_show(this%mergedCoefficients(speciesID))
       ! end do
    end do
    
    !!!Fix basis list size
    do sysI=1, this%numberOfDisplacedSystems
       do speciesID=1, numberOfSpecies
          call Vector_copyConstructorInteger(auxBasisList(speciesID),this%sysBasisList(sysI,speciesID))
          call Vector_constructorInteger(this%sysBasisList(sysI,speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID,this%mergedMolecularSystem), 0)           
          do i=1, size(auxBasisList(speciesID)%values)
             this%sysBasisList(sysI,speciesID)%values(i)=auxBasisList(speciesID)%values(i)
          end do
          ! print *, "sysI", sysI, "speciesID", speciesID, "after list"
          ! call Vector_showInteger(this%sysBasisList(sysI,speciesID))
       end do
    end do
    
    write(*,*) ""
    print *, "Superposed molecular system geometry"
    write(*,*) "---------------------------------- "
    ! call MolecularSystem_showInformation()  
    ! call MolecularSystem_showParticlesInformation()
    call MolecularSystem_copyConstructor(molecularSystem_instance,this%mergedMolecularSystem)
    call MolecularSystem_showCartesianMatrix(molecularSystem_instance)
    particleManager_instance => molecularSystem_instance%allParticles
    call ParticleManager_setOwner()
    call MolecularSystem_saveToFile()

    ! do speciesID=1, numberOfSpecies
    !    write(*,*) ""
    !    write(*,*) " Merged Occupied Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%symbol )
    !    write(*,*) "---------------------------------- "
    !    write(*,*) ""
    !    print *, "contractions", speciesID, int(MolecularSystem_getTotalNumberOfContractions(speciesID),8)
    !    print *, "ocupation", speciesID, int(MolecularSystem_getOcupationNumber(speciesID),8)
    !    call Matrix_constructor(auxCoefficients(speciesID),int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),&
    !         int(MolecularSystem_getOcupationNumber(speciesID),8),0.0_8)
    !    do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
    !       do j=1, MolecularSystem_getOcupationNumber(speciesID)
    !          auxCoefficients(speciesID)%values(i,j)=mergedCoefficients(speciesID)%values(i,j)
    !       end do
    !    end do
    !    call Matrix_show(auxCoefficients(speciesID))                    
    ! end do

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time creating supermolecular system : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
        
    deallocate(auxCoefficients,&
         auxBasisList)
      
    return
    
  end subroutine NOCISuperposed_generateSuperposedSystem

  !>
  !! @brief Generates the NOCI density matrix in the superposed molecular system
  !!
  !! @param this
  !<
  subroutine NOCISuperposed_buildDensityMatrix(this)
    implicit none
    type(NonOrthogonalCI) :: this
   
    type(Matrix) :: molecularOverlapMatrix 
    type(Matrix), allocatable :: inverseOverlapMatrix(:) !,kineticMatrix(:), attractionMatrix(:), externalPotMatrix(:)
    integer :: state
    integer :: i,ii,j,jj,mu,nu, sysI, sysII, speciesID, otherSpeciesID
    integer :: particlesPerOrbital
    integer :: numberOfSpecies
    
    integer :: densUnit
    character(100) :: densFile
    character(50) :: arguments(2), auxString
    type(Matrix), allocatable :: exchangeCorrelationMatrices(:)
    type(Matrix) :: dftEnergyMatrix
    real(8), allocatable :: particlesInGrid(:)
    real(8) :: timeA

    !$  timeA = omp_get_wtime()
    
    if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    numberOfSpecies=molecularSystem_instance%numberOfQuantumSpecies
                
    allocate(InverseOverlapMatrix(numberOfSpecies))

    print *, ""
    print *, "Computing overlap integrals for the superposed systems..."
    print *, ""
    do speciesID = 1, numberOfSpecies
       call DirectIntegralManager_getOverlapIntegrals(molecularSystem_instance,speciesID,this%mergedOverlapMatrix(speciesID))
    end do
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for supermolecular 1-body integrals : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
    
    print *, ""
    print *, "Building merged density matrices for the superposed systems..."
    print *, ""
    !!Build the merged density matrix
    do state=1, CONTROL_instance%CI_STATES_TO_PRINT
       do speciesID=1, numberOfSpecies
          call Matrix_constructor(this%mergedDensityMatrix(state,speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8)
       end do
    end do
    !!Fill the merged density matrix
    ! "Diagonal" terms - same system
    do sysI=1, this%numberOfDisplacedSystems
       do speciesID=1, numberOfSpecies
          particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(sysI))
          do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
             if(this%sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
             do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                if(this%sysBasisList(sysI,speciesID)%values(nu) .eq. 0) cycle
                do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
                   ii=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))*(sysI-1)+i
                   do state=1, CONTROL_instance%CI_STATES_TO_PRINT
                      this%mergedDensityMatrix(state,speciesID)%values(mu,nu) =  this%mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                           this%configurationCoefficients%values(sysI,state)**2*&
                           this%mergedCoefficients(speciesID)%values(mu,ii)*&
                           this%mergedCoefficients(speciesID)%values(nu,ii)*&
                           particlesPerOrbital
                   end do
                end do
             end do
          end do
       end do
    end do
    !!"Non Diagonal" terms - system pairs
    do sysI=1, this%numberOfDisplacedSystems
       do sysII=sysI+1, this%numberOfDisplacedSystems
          if( abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) cycle
          !!Compute molecular overlap matrix and its inverse
          do speciesID=1, numberOfSpecies
             call Matrix_constructor(molecularOverlapMatrix, &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII)),8), 0.0_8 )
             call Matrix_constructor(inverseOverlapMatrix(speciesID), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII)),8), 0.0_8 )
                
             do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysI
                if(this%sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                   if(this%sysBasisList(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
                      ii=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))*(sysI-1)+i
                      do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII))
                         jj=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII))*(sysII-1)+j
                         ! print *, "i, j, mu, nu, coefI, coefII", i,j,mu,nu,mergedCoefficients(speciesID)%values(mu,ii),mergedCoefficients(speciesID)%values(nu,jj)
                         molecularOverlapMatrix%values(i,j)=molecularOverlapMatrix%values(i,j)+&
                              this%mergedCoefficients(speciesID)%values(mu,ii)*&
                              this%mergedCoefficients(speciesID)%values(nu,jj)*&
                              this%mergedOverlapMatrix(speciesID)%values(mu,nu)
                      end do
                   end do
                end do
             end do
             ! print *, "molecularOverlapMatrix sysI, sysII, speciesID", sysI, sysII, speciesID
             ! call Matrix_show(molecularOverlapMatrix)
             if(MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)) .ne. 0) &
                  inverseOverlapMatrix(speciesID)=Matrix_inverse(molecularOverlapMatrix)
          end do
          
          ! Compute density contributions
          do speciesID=1, numberOfSpecies
             particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(sysI))
             do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                if(this%sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                   if(this%sysBasisList(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do state=1, CONTROL_instance%CI_STATES_TO_PRINT
                      do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
                         ii=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))*(sysI-1)+i
                         do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII))
                            jj=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysII))*(sysII-1)+j
                            this%mergedDensityMatrix(state,speciesID)%values(mu,nu) =  this%mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                                 this%configurationCoefficients%values(sysI,state)*&
                                 this%configurationCoefficients%values(sysII,state)*&
                                 this%configurationOverlapMatrix%values(sysI,sysII)*&
                                 inverseOverlapMatrix(speciesID)%values(j,i)*&
                                 this%mergedCoefficients(speciesID)%values(mu,ii)*&
                                 this%mergedCoefficients(speciesID)%values(nu,jj)*&
                                 particlesPerOrbital
                            this%mergedDensityMatrix(state,speciesID)%values(nu,mu) = this%mergedDensityMatrix(state,speciesID)%values(nu,mu) + &
                                 this%configurationCoefficients%values(sysI,state)*&
                                 this%configurationCoefficients%values(sysII,state)*&
                                 this%configurationOverlapMatrix%values(sysI,sysII)*&
                                 inverseOverlapMatrix(speciesID)%values(j,i)*&
                                 this%mergedCoefficients(speciesID)%values(mu,ii)*&
                                 this%mergedCoefficients(speciesID)%values(nu,jj)*&
                                 particlesPerOrbital
                         end do
                      end do
                   end do
                end do
             end do
          end do
          !!symmetrize
          ! do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !    do nu = mu+1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !       this%mergedDensityMatrix(state,speciesID)%values(nu,mu) = this%mergedDensityMatrix(state,speciesID)%values(mu,nu)  
          !    end do
          ! end do
       end do
    end do
    
    !! Open file - to write density matrices
    densUnit = 29
       
    densFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    open(unit = densUnit, file=trim(densFile), status="replace", form="formatted")
    do state=1, CONTROL_instance%CI_STATES_TO_PRINT
       do speciesID=1, numberOfSpecies
          ! print *, "this%mergedDensityMatrix", state, trim( MolecularSystem_instance%species(speciesID)%symbol )
          ! call Matrix_show(this%mergedDensityMatrix(state,speciesID))
          write(auxString,*) state
          arguments(2) = trim(MolecularSystem_instance%species(speciesID)%symbol)
          arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxString)) 
          call Matrix_writeToFile ( this%mergedDensityMatrix(state,speciesID), densUnit , arguments=arguments(1:2) )
       end do
    end do

    ! if(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL.ne."NONE" .or. &
    !      CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL.ne."NONE") then
    !    print *, "Superposed DFT energies:"

    !    allocate(exchangeCorrelationMatrices(numberOfSpecies), &
    !         particlesInGrid(numberOfSpecies))
    !    call DensityFunctionalTheory_buildFinalGrid()
    !    do state=1, CONTROL_instance%CI_STATES_TO_PRINT
    !       call Matrix_constructor(dftEnergyMatrix, int(numberOfSpecies,8), &
    !            int(numberOfSpecies,8), 0.0_8 )
    !       do speciesID=1, numberOfSpecies
    !          call Matrix_constructor(exchangeCorrelationMatrices(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
    !               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8)
    !       end do
    !       call DensityFunctionalTheory_finalDFT(this%mergedDensityMatrix(state,1:numberOfSpecies), &
    !            exchangeCorrelationMatrices, &
    !            dftEnergyMatrix, &
    !            particlesInGrid)

    !       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
    !          do otherSpeciesID = speciesID, MolecularSystem_instance%numberOfQuantumSpecies                
    !             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%symbol ) // &
    !                  "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%symbol ) // &
    !                  " DFT Corr. energy = ", dftEnergyMatrix%values(speciesID,otherSpeciesID)
    !          end do
    !       end do
    !    end do
    ! end if
        
    close(densUnit)

    deallocate(inverseOverlapMatrix)

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for merging density matrices : ", omp_get_wtime() - timeA ," (s)"
    
    return

    ! allocate(kineticMatrix(numberOfSpecies),&
    !      attractionMatrix(numberOfSpecies),&
    !      externalPotMatrix(numberOfSpecies))
    ! do speciesID = 1, numberOfSpecies
       ! call DirectIntegralManager_getKineticIntegrals(molecularSystem_instance,speciesID,kineticMatrix(speciesID))
       ! if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
       !    kineticMatrix(speciesID)%values =  &
       !      kineticMatrix(speciesID)%values * &
       !      ( 1.0_8/MolecularSystem_getMass( speciesID ) -1.0_8 / MolecularSystem_getTotalMass() )
       ! else
       !    kineticMatrix(speciesID)%values =  &
       !      kineticMatrix(speciesID)%values / &
       !      MolecularSystem_getMass( speciesID )
       ! end if

       ! call DirectIntegralManager_getAttractionIntegrals(molecularSystem_instance,speciesID,attractionMatrix(speciesID))
       ! attractionMatrix(speciesID)%values=attractionMatrix(speciesID)%values*(-MolecularSystem_getCharge(speciesID))                         

       ! if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
    !      call DirectIntegralManager_getExternalPotentialIntegrals(molecularSystem_instance,speciesID,externalPotMatrix(speciesID))
    ! end do
    ! write(*,*) ""
    ! write(*,*) "=========================================================="
    ! write(*,*) " ONE BODY ENERGY CONTRIBUTIONS OF THE SUPERPOSED SYSTEMS: "
    ! write(*,*) ""
    ! do state=1, CONTROL_instance%CI_STATES_TO_PRINT
    !    write(*,*) " STATE: ", state
    !    do speciesID=1, numberOfSpecies
    !       write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%symbol ) // &
    !            " Kinetic energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*kineticMatrix(speciesID)%values)
    !       write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%symbol ) // &
    !            "/Fixed interact. energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*attractionMatrix(speciesID)%values)
    !       if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
    !            write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%symbol) // &
    !            " Ext Pot energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*externalPotMatrix(speciesID)%values)
    !       print *, ""
    !    end do
    !    print *, ""
    ! end do
    ! deallocate(kineticMatrix,&
    !      attractionMatrix,&
    !      externalPotMatrix)  
    
  end subroutine NOCISuperposed_buildDensityMatrix

  !>
  !! @brief Generates the NOCI natural orbitals from the NOCI density matrix in the superposed molecular system
  !!
  !! @param this
  !<
  subroutine NOCISuperposed_getNaturalOrbitals(this)
    implicit none
    type(NonOrthogonalCI) :: this

    type(Matrix) :: auxMatrix, densityEigenVectors, auxdensityEigenVectors
    type(Vector) :: auxVector, densityEigenValues, auxdensityEigenValues

    integer :: state
    integer :: i,j,k,speciesID
    integer :: numberOfSpecies

    integer :: densUnit
    character(100) :: densFile
    character(50) :: arguments(2), auxString
    real(8) :: timeA

    !$  timeA = omp_get_wtime()
    if(.not. CONTROL_instance%CI_NATURAL_ORBITALS) return
    if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    numberOfSpecies=molecularSystem_instance%numberOfQuantumSpecies

    write(*,*) ""
    write(*,*) "============================================="
    write(*,*) " NATURAL ORBITALS OF THE SUPERPOSED SYSTEMS: "
    write(*,*) ""
    !! Open file - to write density matrices
    densUnit = 29

    densFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    open(unit = densUnit, file=trim(densFile), status="old", form="formatted", position="append")

    do state=1, CONTROL_instance%CI_STATES_TO_PRINT

       write(*,*) " STATE: ", state

       do speciesID=1, numberOfSpecies

          write(*,*) ""
          write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%symbol )
          write(*,*) "--------------------------------------------------------------"

          call Vector_constructor ( densityEigenValues, &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),4),  0.0_8 )
          call Matrix_constructor ( densityEigenVectors, &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),  0.0_8 )

          call Vector_constructor ( auxdensityEigenValues, &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),4),  0.0_8 )
          call Matrix_constructor ( auxdensityEigenVectors, &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),  0.0_8 )

          !! Lowdin orthogonalization of the density matrix
          auxMatrix = Matrix_pow(this%mergedOverlapMatrix(speciesID), 0.5_8, method="SVD" )

          auxMatrix%values=matmul(matmul(auxMatrix%values,this%mergedDensityMatrix(state,speciesID)%values),auxMatrix%values)

          ! print *, "Diagonalizing non orthogonal CI density Matrix..."

          !! Calcula valores y vectores propios de matriz de densidad CI ortogonal.
          call Matrix_eigen(auxMatrix , auxdensityEigenValues, auxdensityEigenVectors, SYMMETRIC )

          !! Transform back to the atomic basis
          auxMatrix = Matrix_pow(this%mergedOverlapMatrix(speciesID), -0.5_8, method="SVD" )

          auxdensityEigenVectors%values=matmul(auxMatrix%values,auxdensityEigenVectors%values)

          ! reorder and count significant occupations
          k=0
          do i = 1, MolecularSystem_getTotalNumberOfContractions(speciesID)
             densityEigenValues%values(i) =  auxdensityEigenValues%values(MolecularSystem_getTotalNumberOfContractions(speciesID) - i + 1)
             densityEigenVectors%values(:,i) = auxdensityEigenVectors%values(:,MolecularSystem_getTotalNumberOfContractions(speciesID) - i + 1)
             if(abs(densityEigenValues%values(i)) .ge. 1.0E-4_8 ) k=k+1
          end do
          if(k .eq. 0) k=1
          ! Print eigenvectors with occupation larger than 0.01
          call Vector_constructor(auxVector,k,0.0_8)
          call Matrix_constructor(auxMatrix,int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),int(k,8),0.0_8)
          k=0
          do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
             if(abs(densityEigenValues%values(i)) .ge. 1.0E-4_8 ) then
                k=k+1
                auxVector%values(k)=densityEigenValues%values(i)
                do j=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
                   auxMatrix%values(j,k)=densityEigenVectors%values(j,i)
                end do
             end if
          end do
          !densityEigenVectors
          call Matrix_show( auxMatrix , &
               rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
               columnkeys = string_convertvectorofrealstostring( auxVector ),&
               flags=WITH_BOTH_KEYS)

          write(*,"(A10,A10,A20,I5,A15,F17.12)") "number of ", trim(MolecularSystem_getSymbolOfSpecies( speciesID )) ," particles in state", state , &
               " density matrix: ", sum( transpose(this%mergedDensityMatrix(state,speciesID)%values)*this%mergedOverlapMatrix(speciesID)%values)
          write(*,"(A10,A10,A40,F17.12)") "sum of ", trim(MolecularSystem_getSymbolOfSpecies( speciesID )) , "natural orbital occupations", sum(densityEigenValues%values)

          ! density matrix check
          ! auxMatrix%values=0.0
          ! do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
          !    do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
          !       do k=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
          !          auxMatrix%values(mu,nu)=auxMatrix%values(mu,nu)+densityEigenValues%values(k)*&
          !               densityEigenVectors%values(mu,k)*densityEigenVectors%values(nu,k)
          !       end do
          !    end do
          ! end do          
          ! print *, "atomicDensityMatrix again"
          ! call Matrix_show(auxMatrix)

          write(auxString,*) state

          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%symbol )
          arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 

          call Matrix_writeToFile ( densityEigenVectors, densUnit , arguments=arguments(1:2) )

          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%symbol )
          arguments(1) = "OCCUPATIONS"//trim(adjustl(auxstring))

          call Vector_writeToFile( densityEigenValues, densUnit, arguments=arguments(1:2) )

          write(*,*) " End of natural orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%symbol )
       end do
    end do

    write(*,*) ""
    write(*,*) " END OF NATURAL ORBITALS"
    write(*,*) "=============================="
    write(*,*) ""

    close(densUnit)

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI natural orbitals : ", omp_get_wtime() - timeA ," (s)"

    return

  end subroutine NOCISuperposed_getNaturalOrbitals
  
  !>
  !! @brief Save NOCI results to file
  !!
  !! @param 
  !<
  subroutine NOCISuperposed_saveToFile(this)
    type(NonOrthogonalCI) :: this
    integer :: nociUnit, speciesID, numberOfSpecies, sysI
    character(100) :: prefix, nociFile
    character(50) :: arguments(2), auxString

    !Save merged molecular system
    call MolecularSystem_copyConstructor(molecularSystem_instance,this%mergedMolecularSystem)

    prefix=trim(CONTROL_instance%INPUT_FILE)//"NOCI"
    call MolecularSystem_saveToFile(prefix)

    numberOfSpecies=molecularSystem_instance%numberOfQuantumSpecies

    nociUnit=123
    nociFile=trim(prefix)//".states"
    open(unit = nociUnit, file=trim(nociFile), status="replace", form="unformatted")

    arguments(1:1) = "NOCI-NUMBEROFDISPLACEDSYSTEMS"
    call Vector_writeToFileInteger(unit=nociUnit, binary=.true., value=this%numberOfDisplacedSystems, arguments=arguments(1:1) )

    arguments(1:1) = "NOCI-NUMBEROFSPECIES"
    call Vector_writeToFileInteger(unit=nociUnit, binary=.true., value=numberOfSpecies, arguments=arguments(1:1) )

    arguments(1:1) = "NOCI-CONFIGURATIONCOEFFICIENTS"
    call Matrix_writeToFile ( this%configurationCoefficients, nociUnit , binary=.true., arguments=arguments(1:1) )       

    arguments(1:1) = "NOCI-CONFIGURATIONENERGIES"
    call Vector_writeToFile ( this%statesEigenvalues, nociUnit , binary=.true., arguments=arguments(1:1) )       
    
    arguments(1) = "MERGEDCOEFFICIENTS"
    do speciesID=1, numberOfSpecies
       arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
       call Matrix_writeToFile ( this%mergedCoefficients(speciesID), nociUnit, binary=.true. , arguments=arguments(1:2) )       
    end do

    do sysI=1, this%numberOfDisplacedSystems
       do speciesID=1, numberOfSpecies
          write(auxString,*) sysI
          arguments(1) = "SYSBASISLIST"//trim(auxString) 
          arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
          call Vector_writeToFileInteger(this%sysBasisList(sysI,speciesID), nociUnit, binary=.true., arguments=arguments(1:2) )
       end do
    end do

    ! do state=1, min(CONTROL_instance%NUMBER_OF_CI_STATES,this%numberOfDisplacedSystems)
    ! end do
    
    ! do state=1, CONTROL_instance%CI_STATES_TO_PRINT
    !       write(auxString,*) state
    !       call Matrix_writeToFile ( this%mergedDensityMatrix(state,speciesID), densUnit , arguments=arguments(1:2) )
    !    end do
    ! end do

    close(nociUnit)
    
  end subroutine NOCISuperposed_saveToFile

  
end module NOCISuperposed_

