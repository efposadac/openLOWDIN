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

module NOCIFranckCondon_
  use NOCIBuild_
  use NOCIMatrices_
  use MolecularSystem_
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
       NOCIFranckCondon_computeFranckCondon

  private

contains
  
  !>
  !! @brief Compute Franck-Condon factors from the current NOCI calculations and previous results read from file
  !!
  !! @param 
  !<
  subroutine NOCIFranckCondon_computeFranckCondon(this)
    type(NonOrthogonalCI) :: this
    integer :: nociUnit, numberOfSpecies, occupationNumber,numberOfDisplacedSystems, numberOfContractions, dim2
    character(100) :: nociFile
    type(Matrix) :: ciCoefficients
    type(Vector) :: ciEnergies
    type(Matrix), allocatable :: auxCoefficients(:), superMergedCoefficients(:)
    type(IVector), allocatable :: sysListCur(:,:), sysListRef(:,:), orbListI(:), orbListII(:)
    type(IVector) :: auxIVector
    type(MolecularSystem) :: superMergedMolecularSystem
    logical :: existFile
    type(Matrix) :: molecularOverlapMatrix
    type(Matrix), allocatable :: superOverlapMatrix(:), superMomentMatrix(:,:), inverseOverlapMatrix(:), molecularMomentMatrix(:,:) !,attractionMatrix(:), externalPotMatrix(:)
    integer :: stateI, stateII
    integer :: i,ii,j,jj,k,mu,nu,mumu,nunu,sysI, sysII, speciesID, otherSpeciesID
    integer :: particlesPerOrbital
    real(8) :: overlapDeterminant, trololo, trolololo(3), pointchargesdipole(3)
    
    integer :: densUnit
    character(100) :: densFile
    character(50) :: arguments(2), auxString
    type(Matrix), allocatable :: franckCondonMatrix(:), transitionDipoleMatrix(:,:), refCurOverlapMatrix(:), refCurMomentMatrix(:,:)
    type(Matrix) :: refCurTotalOverlap
    real(8) :: timeA

    !$  timeA = omp_get_wtime()
    
    existFile=.false.
    
    nociFile = trim(CONTROL_instance%INPUT_FILE)//"refNOCI"
    inquire( FILE = trim(nociFile)//".sys", EXIST = existFile )

    if(.not. existFile) return
    print *, "Found a reference molecular system for NOCI calculations ",  trim(nociFile)//".sys"

    pointchargesdipole=0.0
    do i=1, size( MolecularSystem_instance%pointCharges )      
       pointchargesdipole = pointchargesdipole + MolecularSystem_instance%pointCharges(i)%origin(:) * MolecularSystem_instance%pointCharges(i)%charge
    end do

    
    call MolecularSystem_loadFromFile("LOWDIN.SYS",nociFile)
    call MolecularSystem_showInformation()  
    call MolecularSystem_showParticlesInformation()
    call MolecularSystem_showCartesianMatrix()

    nociFile = trim(CONTROL_instance%INPUT_FILE)//"refNOCI.states"
    inquire( FILE = trim(nociFile), EXIST = existFile )

    if(.not. existFile) then
       print *, "Did not find reference states for NOCI calculations ",  nociFile
       return
    end if
    print *, "Found reference states for NOCI calculations ",  nociFile
    print *, "Computing the Franck-Condon factors with respect to that system" 
    
    nociUnit=123
    open(unit = nociUnit, file=trim(nociFile), status="old", form="unformatted")

    arguments(1) = "NOCI-NUMBEROFSPECIES"
    call Vector_getFromFileInteger(1,unit=nociUnit, binary=.true., value=numberOfSpecies, arguments=arguments(1:1) )

    arguments(1) = "NOCI-NUMBEROFDISPLACEDSYSTEMS"
    call Vector_getFromFileInteger(1,unit=nociUnit, binary=.true., value=numberOfDisplacedSystems, arguments=arguments(1:1) )

    allocate(auxCoefficients(numberOfSpecies))           
    allocate(sysListCur(numberOfDisplacedSystems,numberOfSpecies),sysListRef(numberOfDisplacedSystems,numberOfSpecies))
    allocate(orbListI(numberOfDisplacedSystems),orbListII(numberOfDisplacedSystems))
    allocate(superMergedCoefficients(numberOfSpecies))           
    allocate(superOverlapMatrix(numberOfSpecies), superMomentMatrix(numberOfSpecies,3),inverseOverlapMatrix(numberOfSpecies),molecularMomentMatrix(numberOfSpecies,3))
    allocate(franckCondonMatrix(numberOfSpecies),transitionDipoleMatrix(numberOfSpecies+1,3),refCurOverlapMatrix(numberOfSpecies),refCurMomentMatrix(numberOfSpecies,3))
    
    arguments(1) = "NOCI-CONFIGURATIONCOEFFICIENTS"
    ciCoefficients = Matrix_getFromFile(numberOfDisplacedSystems,numberOfDisplacedSystems,nociUnit,binary=.true.,arguments=arguments(1:1) )       
    
    arguments(1:1) = "NOCI-CONFIGURATIONENERGIES"
    call Vector_getFromFile(numberOfDisplacedSystems, nociUnit, output=ciEnergies, binary=.true., arguments=arguments(1:1) )
    
    arguments(1) = "MERGEDCOEFFICIENTS"
    do speciesID=1, numberOfSpecies
       numberOfContractions=molecularSystem_getTotalNumberOfContractions(speciesID)
       dim2=max(MolecularSystem_getTotalNumberOfContractions(speciesID),MolecularSystem_getOcupationNumber(speciesID))
       arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
       auxCoefficients(speciesID) = Matrix_getFromFile(numberOfContractions,dim2,nociUnit,binary=.true.,arguments=arguments(1:2) )       
    end do

    do sysI=1, numberOfDisplacedSystems
       do speciesID=1, numberOfSpecies
          numberOfContractions=molecularSystem_getTotalNumberOfContractions(speciesID)
          write(auxString,*) sysI
          arguments(1) = "SYSBASISLIST"//trim(auxString) 
          arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
          call Vector_getFromFileInteger(numberOfContractions, nociUnit, output=sysListRef(sysI,speciesID), binary=.true., arguments=arguments(1:2) )
       end do
    end do

    close(nociUnit)
    
    !Create a super-mega molecular system
    !Merge coefficients from NOCI calculation and reference system

    print *, "super-mega molecular system"
    call MolecularSystem_mergeTwoSystems(superMergedMolecularSystem, this%mergedMolecularSystem, MolecularSystem_instance, &
         orbListI(:),orbListII(:), reorder=.false.)
    call MolecularSystem_showInformation(superMergedMolecularSystem)  
    call MolecularSystem_showParticlesInformation(superMergedMolecularSystem)
    call MolecularSystem_showCartesianMatrix(superMergedMolecularSystem)

    call NOCIMatrices_mergeCoefficients(this%mergedMolecularSystem,MolecularSystem_instance,superMergedMolecularSystem,&
         this%mergedCoefficients(:),auxCoefficients(:),&
         orbListI(:),orbListII(:),superMergedCoefficients(:))

    ! do speciesID=1, numberOfSpecies
    !    print *, "superMergedCoefficients", speciesID
    !    call Matrix_show(superMergedCoefficients(speciesID))
    ! end do
    
    !Fix basis list size
    do speciesID=1, numberOfSpecies
       ! print *, "orbListI", "speciesID", speciesID
       ! call Vector_showInteger(orbListI(speciesID))
       do sysI=1, this%numberOfDisplacedSystems
          call Vector_copyConstructorInteger(auxIVector,this%sysBasisList(sysI,speciesID))
          call Vector_constructorInteger(sysListCur(sysI,speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem), 0)           
          do i=1, size(auxIVector%values)
             if(orbListI(speciesID)%values(i) .eq. 0) cycle
             sysListCur(sysI,speciesID)%values(i)=auxIVector%values(orbListI(speciesID)%values(i))
          end do
          ! print *, "sysListCur", "sysI", sysI, "speciesID", speciesID
          ! call Vector_showInteger(sysListCur(sysI,speciesID))
       end do
    end do

    do speciesID=1, numberOfSpecies
       occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(1)) !not using the merged molecular systems
       ! print *, "orbListII", "speciesID", speciesID
       ! call Vector_showInteger(orbListII(speciesID))
       do sysII=1, numberOfDisplacedSystems
          call Vector_copyConstructorInteger(auxIVector,sysListRef(sysII,speciesID))
          call Vector_constructorInteger(sysListRef(sysII,speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem), 0)           
          do i=1, size(orbListII(speciesID)%values)
             if(orbListII(speciesID)%values(i) .eq. 0) cycle
             sysListRef(sysII,speciesID)%values(i)=auxIVector%values(orbListII(speciesID)%values(i))
          end do
          ! print *, "sysListRef", "sysII", sysII, "speciesID", speciesID
          ! call Vector_showInteger(sysListRef(sysII,speciesID))
       end do
    end do
    
    ! if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    ! numberOfSpecies=molecularSystem_instance%numberOfQuantumSpecies
                

    print *, ""
    print *, "Computing overlap and moment integrals for the super-mega system..."
    print *, ""
    do speciesID = 1, numberOfSpecies
       call DirectIntegralManager_getOverlapIntegrals(superMergedMolecularSystem,speciesID,superOverlapMatrix(speciesID))
       call DirectIntegralManager_getMomentIntegrals(superMergedMolecularSystem,speciesID,1,superMomentMatrix(speciesID,1))
       call DirectIntegralManager_getMomentIntegrals(superMergedMolecularSystem,speciesID,2,superMomentMatrix(speciesID,2))
       call DirectIntegralManager_getMomentIntegrals(superMergedMolecularSystem,speciesID,3,superMomentMatrix(speciesID,3))
    end do
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for supermolecular 1-body integrals : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
    
    print *, ""
    print *, "Self overlap matrices for the supermegaposed systems..."
    print *, ""

    do speciesID=1, numberOfSpecies
       call Matrix_constructor(refCurOverlapMatrix(speciesID), int(this%numberOfDisplacedSystems,8), &
            int(numberOfDisplacedSystems,8), 1.0_8)
    end do
    !!Fill the merged density matrix
    !!"Non Diagonal" terms - system pairs
    do sysI=1, numberOfDisplacedSystems !computed
       do sysII=1, numberOfDisplacedSystems !reference
          ! if( abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) cycle
          !!Compute molecular overlap matrix and its inverse
          do speciesID=1, numberOfSpecies
             occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)) !not using the merged molecular systems
             particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(sysI))
             call Matrix_constructor(molecularOverlapMatrix, int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
             ! call Matrix_constructor(inverseOverlapMatrix, int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
                
             do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem) !sysI
                if(sysListRef(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem)  !sysII
                   if(sysListRef(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do i = 1 , occupationNumber
                      ii=occupationNumber*(sysI-1)+i+MolecularSystem_getOcupationNumber(speciesID,superMergedMolecularSystem)/2
                      do j = 1 , occupationNumber
                         jj=occupationNumber*(sysII-1)+j+MolecularSystem_getOcupationNumber(speciesID,superMergedMolecularSystem)/2
                         ! print *, "i, j, mu, nu, coefI, coefII, overlap", i,j,mu,nu,superMergedCoefficients(speciesID)%values(mu,ii),&
                         !      superMergedCoefficients(speciesID)%values(nu,jj),&
                         !      superOverlapMatrix(speciesID)%values(mu,nu)
                         molecularOverlapMatrix%values(i,j)=molecularOverlapMatrix%values(i,j)+&
                              superMergedCoefficients(speciesID)%values(mu,ii)*&
                              superMergedCoefficients(speciesID)%values(nu,jj)*&
                              superOverlapMatrix(speciesID)%values(mu,nu)
                      end do
                   end do
                end do
             end do
             if(occupationNumber .ne. 0) then
                ! inverseOverlapMatrix=Matrix_inverse(molecularOverlapMatrix)
                ! print *, "inverseOverlapMatrices sysI, sysII", speciesID, sysI, sysII
                ! call Matrix_show(inverseOverlapMatrices(speciesID))
                call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant,method="LU")             
                ! print *, "OverlapDeterminantLU speciesID, sysI, sysII", speciesID, sysI, sysII, overlapDeterminant
             else
                overlapDeterminant=1.0
             end if
             refCurOverlapMatrix(speciesID)%values(sysI,sysII)=refCurOverlapMatrix(speciesID)%values(sysI,sysII)*overlapDeterminant**particlesPerOrbital
          end do
          
       end do
    end do
    
    do speciesID=1, numberOfSpecies
       print *, "Reference Overlap Matrix for", speciesID
       call Matrix_show(refCurOverlapMatrix(speciesID))
    end do

    print *, ""
    print *, "Building Franck-Condon matrices for the superposed systems..."
    print *, ""

    do speciesID=1, numberOfSpecies
       call Matrix_constructor(refCurOverlapMatrix(speciesID), int(this%numberOfDisplacedSystems,8), &
            int(numberOfDisplacedSystems,8), 1.0_8)
       do k=1,3
          call Matrix_constructor(refCurMomentMatrix(speciesID,k), int(this%numberOfDisplacedSystems,8), &
               int(numberOfDisplacedSystems,8), 0.0_8)
       end do
    end do
    call Matrix_constructor(refCurTotalOverlap, int(this%numberOfDisplacedSystems,8), &
         int(numberOfDisplacedSystems,8), 1.0_8)
    
    !!Fill the merged density matrix
    !!"Non Diagonal" terms - system pairs
    do sysI=1, this%numberOfDisplacedSystems !computed
       do sysII=1, numberOfDisplacedSystems !reference
          ! if( abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) cycle
          !!Compute molecular overlap matrix and its inverse
          do speciesID=1, numberOfSpecies
             occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)) !not using the merged molecular systems
             particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(sysI))
             call Matrix_constructor(molecularOverlapMatrix, int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
             do k=1,3
                call Matrix_constructor(molecularMomentMatrix(speciesID,k), int(occupationNumber,8), int(occupationNumber,8), 0.0_8 )
             end do
                
             do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem) !sysI
                if(sysListCur(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,superMergedMolecularSystem)  !sysII
                   if(sysListRef(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do i = 1 , occupationNumber
                      ii=occupationNumber*(sysI-1)+i
                      do j = 1 , occupationNumber
                         jj=occupationNumber*(sysII-1)+j+MolecularSystem_getOcupationNumber(speciesID,superMergedMolecularSystem)/2
                         ! print *, "i, j, mu, nu, coefI, coefII, overlap", i,j,mu,nu,superMergedCoefficients(speciesID)%values(mu,ii),&
                         !      superMergedCoefficients(speciesID)%values(nu,jj),&
                         !      superOverlapMatrix(speciesID)%values(mu,nu)
                         molecularOverlapMatrix%values(i,j)=molecularOverlapMatrix%values(i,j)+&
                              superMergedCoefficients(speciesID)%values(mu,ii)*&
                              superMergedCoefficients(speciesID)%values(nu,jj)*&
                              superOverlapMatrix(speciesID)%values(mu,nu)
                         do k=1,3
                            molecularMomentMatrix(speciesID,k)%values(i,j)=molecularMomentMatrix(speciesID,k)%values(i,j)+&
                                 superMergedCoefficients(speciesID)%values(mu,ii)*&
                                 superMergedCoefficients(speciesID)%values(nu,jj)*&
                                 superMomentMatrix(speciesID,k)%values(mu,nu)
                         end do
                      end do
                   end do
                end do
             end do
             if(occupationNumber .ne. 0) then
                inverseOverlapMatrix(speciesID)=Matrix_inverse(molecularOverlapMatrix)
                ! print *, "inverseOverlapMatrices sysI, sysII", speciesID, sysI, sysII
                ! call Matrix_show(inverseOverlapMatrices(speciesID))
                call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant,method="LU")             
                ! print *, "OverlapDeterminantLU speciesID, sysI, sysII", speciesID, sysI, sysII, overlapDeterminant
                refCurOverlapMatrix(speciesID)%values(sysI,sysII)=overlapDeterminant**particlesPerOrbital
             else
                overlapDeterminant=1.0
             end if
             refCurTotalOverlap%values(sysI,sysII)=refCurTotalOverlap%values(sysI,sysII)*refCurOverlapMatrix(speciesID)%values(sysI,sysII)
          end do

          do speciesID=1, numberOfSpecies
             occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI)) !not using the merged molecular systems
             particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(sysI))
             do i = 1 , occupationNumber
                do j = 1 , occupationNumber
                   do k=1,3
                      refCurMomentMatrix(speciesID,k)%values(sysI,sysII)=refCurMomentMatrix(speciesID,k)%values(sysI,sysII)+&
                           molecularMomentMatrix(speciesID,k)%values(i,j)*&
                           inverseOverlapMatrix(speciesID)%values(j,i)
                      end do
                   end do
                end do
                do k=1,3
                   refCurMomentMatrix(speciesID,k)%values(sysI,sysII)=refCurMomentMatrix(speciesID,k)%values(sysI,sysII)*refCurTotalOverlap%values(sysI,sysII)*particlesPerOrbital
                end do
          end do
       end do
    end do

    do speciesID=1, numberOfSpecies
       print *, "refCurOverlapMatrix(speciesID)", speciesID
       call Matrix_show(refCurOverlapMatrix(speciesID))
       call Matrix_constructor(franckCondonMatrix(speciesID), int(CONTROL_instance%CI_STATES_TO_PRINT,8), int(CONTROL_instance%CI_STATES_TO_PRINT,8), 0.0_8)
    end do

    !+1 For point charges
    do speciesID=1, numberOfSpecies+1
       do k=1,3
          call Matrix_constructor(transitionDipoleMatrix(speciesID,k), int(CONTROL_instance%CI_STATES_TO_PRINT,8), int(CONTROL_instance%CI_STATES_TO_PRINT,8), 0.0_8)
       end do
    end do
    
    do stateII=1, CONTROL_instance%CI_STATES_TO_PRINT
       print *, "Reference state:", stateII
       do stateI=1, CONTROL_instance%CI_STATES_TO_PRINT
          print *, "            current state:", stateI
          do speciesID=1, numberOfSpecies
             occupationNumber=MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(1)) !not using the merged molecular systems
             print *, "occupationNumber", occupationNumber
             particlesPerOrbital=MolecularSystem_getEta(speciesID,this%molecularSystems(1))
             trololo=0
             do sysI=1, this%numberOfDisplacedSystems !computed
                do sysII=1, numberOfDisplacedSystems !reference
                   do i = 1 , occupationNumber
                      do j = 1 , occupationNumber
                         trololo = trololo + &
                              inverseOverlapMatrix(speciesID)%values(j,i)*&
                              this%configurationCoefficients%values(sysI,stateI)*&
                              ciCoefficients%values(sysII,stateII)*& !!reference
                              refCurOverlapMatrix(speciesID)%values(sysI,sysII)*&
                              particlesPerOrbital
                      end do
                   end do
                              ! refCurTotalOverlap%values(sysI,sysII)*&
                   ! franckCondonMatrix(speciesID)%values(stateI,stateII)+&

                   do k=1,3
                      transitionDipoleMatrix(speciesID,k)%values(stateI,stateII) = transitionDipoleMatrix(speciesID,k)%values(stateI,stateII) + &
                           molecularsystem_getcharge( speciesID )*&
                           this%configurationCoefficients%values(sysI,stateI)*&
                           ciCoefficients%values(sysII,stateII)*& !!reference
                           refCurMomentMatrix(speciesID,k)%values(sysI,sysII)
                   end do

                end do
             end do
             print *, "speciesID", speciesID, "trololo", trololo
             franckCondonMatrix(speciesID)%values(stateI,stateII)=trololo
             franckCondonMatrix(speciesID)%values(stateI,stateII)=franckCondonMatrix(speciesID)%values(stateI,stateII)/(occupationNumber*particlesPerOrbital)
             print *, "                      F.C. factor for ", molecularSystem_getNameOfSpecies(speciesID),&
                  franckCondonMatrix(speciesID)%values(stateI,stateII)
          end do
          do sysI=1, this%numberOfDisplacedSystems !computed
             do sysII=1, numberOfDisplacedSystems !reference
                do k=1,3
                   transitionDipoleMatrix(numberOfSpecies+1,k)%values(stateI,stateII) = transitionDipoleMatrix(numberOfSpecies+1,k)%values(stateI,stateII) + &
                        pointchargesdipole(k)*&
                        this%configurationCoefficients%values(sysI,stateI)*&
                        ciCoefficients%values(sysII,stateII)*& !!reference
                        refCurTotalOverlap%values(sysI,sysII)
                end do
             end do
          end do
          ! trololo=1
          ! do speciesID=1, numberOfSpecies
          !         trololo=trololo*franckCondonMatrix(speciesID)%values(stateI,stateII)
          ! end do
          ! print *, "                      F.C. factor product ", trololo
          ! trololo=0
          ! do speciesID=1, numberOfSpecies
          !         trololo=trololo+franckCondonMatrix(speciesID)%values(stateI,stateII)
          ! end do
          ! print *, "                      F.C. factor sum ", trololo
          ! trololo=0
          ! do sysI=1, this%numberOfDisplacedSystems !computed
          !    do sysII=1, numberOfDisplacedSystems !reference
          !          trololo = trololo + &
          !               this%configurationCoefficients%values(sysI,stateI)*&
          !               ciCoefficients%values(sysII,stateII)*& !!reference
          !               refCurTotalOverlap%values(sysI,sysII)
          !    end do
          ! end do
          ! print *, "                      total overlap ", trololo
       end do
    end do

    print *, "Dipole approximation spectrum"    
    do stateII=1, CONTROL_instance%CI_STATES_TO_PRINT
       print *, "Reference state:", stateII
       do stateI=1, CONTROL_instance%CI_STATES_TO_PRINT
          trolololo=0
          print *, "current state:", stateI
          do speciesID=1, numberOfSpecies
             do k=1,3
                trolololo(k)=trolololo(k)+transitionDipoleMatrix(speciesID,k)%values(stateI,stateII)
             end do
             print *, "                      T.D. integrals for ", molecularSystem_getNameOfSpecies(speciesID),&
                  transitionDipoleMatrix(speciesID,1)%values(stateI,stateII),&
                  transitionDipoleMatrix(speciesID,2)%values(stateI,stateII),&
                  transitionDipoleMatrix(speciesID,3)%values(stateI,stateII)
          end do
          do k=1,3
             trolololo(k)=trolololo(k)+transitionDipoleMatrix(numberOfSpecies+1,k)%values(stateI,stateII)
          end do
          print *, "                      T.D. integrals point charges ", &
               transitionDipoleMatrix(numberOfSpecies+1,1)%values(stateI,stateII),&
               transitionDipoleMatrix(numberOfSpecies+1,2)%values(stateI,stateII),&
               transitionDipoleMatrix(numberOfSpecies+1,3)%values(stateI,stateII)
          print *, "energy dif", ciEnergies%values(stateII)-this%statesEigenvalues%values(stateI), "total components", trolololo(1:3) ,"intensity", sqrt(sum(trolololo(1:3)**2))
       end do
    end do
    
    close(densUnit)

    deallocate(auxCoefficients,&           
         sysListCur,sysListRef,&
         orbListI,orbListII,&
         superMergedCoefficients,&          
         superOverlapMatrix,&
         franckCondonMatrix)
    
  end subroutine NOCIFranckCondon_computeFranckCondon

  
end module NOCIFranckCondon_

