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

module NonOrthogonalCI_
  use Exception_
  use MolecularSystem_
  use Matrix_
  use ReadTransformedIntegrals_
  use Lebedev_
  use Matrix_
  use Vector_
  use Solver_
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

  type, public :: NonOrthogonalCI
     logical :: isInstanced
     integer :: numberOfDisplacedSystems
     integer :: numberOfRejectedSystems
     integer :: numberOfTransformedCenters
     integer :: numberOfIndividualTransformations
     integer, allocatable :: rotationCenterList(:,:)
     type(Matrix) :: configurationOverlapMatrix, configurationHamiltonianMatrix, configurationCoefficients
     type(Vector) :: statesEigenvalues
     type(Matrix), allocatable :: HFCoefficients(:,:)
     type(Matrix), allocatable :: HCoreMatrices(:,:,:)
     type(Matrix), allocatable :: inverseOverlapMatrices(:,:,:)
     character(50) :: transformationType
     character(15),allocatable :: systemLabels(:)
     real(8) :: refEnergy
  end type NonOrthogonalCI

  type(NonOrthogonalCI), public :: NonOrthogonalCI_instance

  public :: &
       NonOrthogonalCI_constructor,&
       NonOrthogonalCI_displaceGeometries,&
       NonOrthogonalCI_firstImplementation,& 
       NonOrthogonalCI_configurationsOverlapAndOneParticle,&
       NonOrthogonalCI_diagonalizeCImatrix,&
       NonOrthogonalCI_plotDensities

  private

contains

  !>
  !! @brief Allocates memory and run HF calculations to be used in the construction of the NOCI matrix
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_constructor(this)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: numberOfRotationCenters, numberOfTranslationCenters
    integer :: p,q,r
    
    print *, "-------------------------------------------------------------"
    print *, "STARTING NON ORTHOGONAL CONFIGURATION INTERACTION CALCULATION"
    print *, "-------------------------------------------------------------"
    print *, ""
    this%isInstanced=.true.
    this%numberOfDisplacedSystems=0
    this%numberOfRejectedSystems=0

    numberOfTranslationCenters=0
    numberOfRotationCenters=0

    allocate(this%rotationCenterList(size(MolecularSystem_instance%allParticles),2))
    !For rotations, 0,0: leave alone, N,M: rotation center number to be rotated around point M
    this%rotationCenterList=0
    
    !!Translation count
    do p = 1, size(MolecularSystem_instance%allParticles)

       if(MolecularSystem_instance%allParticles(p)%particlePtr%translationCenter.gt.numberOfTranslationCenters) &
            numberOfTranslationCenters=MolecularSystem_instance%allParticles(p)%particlePtr%translationCenter

       if(MolecularSystem_instance%allParticles(p)%particlePtr%translationCenter.ne.0) &
            write (*,"(A,A10,A,3F9.5,A)") "Particle ", trim(ParticleManager_getSymbol(p)), &
            " basis functions at ", MolecularSystem_instance%allParticles(p)%particlePtr%origin(1:3), &
            " are going to be displaced"
    end do

    !!Rotation count
    do p = 1, size(MolecularSystem_instance%allParticles)
       if(MolecularSystem_instance%allParticles(p)%particlePtr%rotationPoint.ne.0) then
          write (*,"(A,A10,A,3F9.5,A,I5)") "Particle ", trim(ParticleManager_getSymbol(p)), &
               " located at ", MolecularSystem_instance%allParticles(p)%particlePtr%origin(1:3), &
               " is center of rotation number", MolecularSystem_instance%allParticles(p)%particlePtr%rotationPoint

          do q = 1, size(MolecularSystem_instance%allParticles)
             if(MolecularSystem_instance%allParticles(q)%particlePtr%rotateAround .eq. &
                  MolecularSystem_instance%allParticles(p)%particlePtr%rotationPoint) then
                write (*,"(A,A10,A,3F9.5,A,I5)") "Particle ", trim(ParticleManager_getSymbol(q)), &
                     " basis functions at ", MolecularSystem_instance%allParticles(q)%particlePtr%origin(1:3), &
                     " are going to be rotated around center ", MolecularSystem_instance%allParticles(q)%particlePtr%rotateAround

                if(q .eq. MolecularSystem_instance%allParticles(q)%particlePtr%owner) then
                   !in the case of several species with the same center, rotate them as one
                   numberOfRotationCenters=numberOfRotationCenters+1               
                   this%rotationCenterList(q,1)=numberOfRotationCenters
                   !find childs
                   do r=1,size(MolecularSystem_instance%allParticles(q)%particlePtr%childs)
                      this%rotationCenterList( MolecularSystem_instance%allParticles(q)%particlePtr%childs(r),1)=numberOfRotationCenters
                   end do
                end if
                this%rotationCenterList(q,2)=p
             end if
          end do
       end if
    end do
    print *, ""

    print *, "this%rotationCenterList" 
    do p=1, size(MolecularSystem_instance%allParticles)
       print *, "Particle ", trim(ParticleManager_getSymbol(p)),this%rotationCenterList(p,1), this%rotationCenterList(p,2)
    end do
    
    if(numberOfTranslationCenters.ne.0) then

       print *, ""
       write (*,"(A,I5,A,I3,A,I3,A,I3,A)") "Displacing coordinates of ", numberOfTranslationCenters, " centers", &
            CONTROL_instance%TRANSLATION_SCAN_GRID(1),"*",&
            CONTROL_instance%TRANSLATION_SCAN_GRID(2),"*",&
            CONTROL_instance%TRANSLATION_SCAN_GRID(3)," times"
       print *, ""

       this%transformationType="TRANSLATION"
       this%numberOfTransformedCenters=numberOfTranslationCenters
       this%numberOfIndividualTransformations=&
            CONTROL_instance%TRANSLATION_SCAN_GRID(1)*CONTROL_instance%TRANSLATION_SCAN_GRID(2)*CONTROL_instance%TRANSLATION_SCAN_GRID(3)

    else if(numberOfRotationCenters.ne.0) then
       print *, ""
       write (*,"(A,I5,A,I5,A,I5,A)") "Rotating coordinates of ", numberOfRotationCenters, " centers", CONTROL_instance%ROTATIONAL_SCAN_GRID, &
            " times in ", CONTROL_instance%NESTED_ROTATIONAL_GRIDS, " nested grids"
            
       print *, ""

       this%transformationType="ROTATION"
       this%numberOfTransformedCenters=numberOfRotationCenters
       this%numberOfIndividualTransformations=CONTROL_instance%ROTATIONAL_SCAN_GRID*CONTROL_instance%NESTED_ROTATIONAL_GRIDS

    end if
        
  end subroutine NonOrthogonalCI_constructor
  !>
  !! @brief Generates different geometries and runs HF calculations at each 
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_displaceGeometries(this)
    implicit none
    type(NonOrthogonalCI) :: this

    type(MolecularSystem) :: originalMolecularSystem
    real(8) :: testEnergy
    character(50) :: wfnFile, molFile
    character(50) :: arguments(2)
    integer, allocatable :: transformationCounter(:)
    integer :: wfnUnit
    integer :: i,j
    integer :: sysI, speciesID
    real(8) :: timeA
    real(8) :: timeB

    !$  timeA = omp_get_wtime()
    
    !!Read HF energy of the non displaced SCF calculation 
    wfnUnit=300
    wfnFile="lowdin.wfn"
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=this%refEnergy, arguments=["TOTALENERGY"])
    close(300)
    ! print *, "HF reference energy is ", hfEnergy

    call MolecularSystem_copyConstructor(originalMolecularSystem, molecularSystem_instance)

    allocate(transformationCounter(this%numberOfTransformedCenters))

    transformationCounter(1:this%numberOfTransformedCenters)=1

    this%numberOfDisplacedSystems=0

    print *, "this%numberOfTransformedCenters", this%numberOfTransformedCenters
    print *, "this%numberOfIndividualTransformations", this%numberOfIndividualTransformations
!!!!! clock type iterations to form all the possible combination of modified geometries
    do while (.true.)

       !Apply the transformation given by transformationCounter to each center, the result is saved in molecularSystemInstance
       call NonOrthogonalCI_transformCoordinates(this,transformationCounter(1:this%numberOfTransformedCenters),originalMolecularSystem)
       
       write (*,"(A,A)") " Running system at transformation", MolecularSystem_instance%description

       call MolecularSystem_showCartesianMatrix()

!!!Run a HF calculation at each displaced geometry
       call MolecularSystem_saveToFile()
       call NonOrthogonalCI_runHF(testEnergy)

       !!Screen geometries with high energies
       if( CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD .ne. 0.0 .and. &
            testEnergy .gt. this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD) then
          print *, " Skipping system with energy", testEnergy
          this%numberOfRejectedSystems=this%numberOfRejectedSystems+1                      
       else
          this%numberOfDisplacedSystems=this%numberOfDisplacedSystems+1
          write (*,"(A,I10,A,F20.12)") "Saving system with ID ", this%numberOfDisplacedSystems, " and energy", testEnergy
          call NonOrthogonalCI_saveSystem(this%numberOfDisplacedSystems)
       end if
       print *, ""

       !Determine the next movement like a clock iteration
       transformationCounter(1)=transformationCounter(1)+1
       do i=1,this%numberOfTransformedCenters-1
          if(transformationCounter(i) .gt. this%numberOfIndividualTransformations) then
             j=i+1 
             transformationCounter(j)=transformationCounter(j)+1
             transformationCounter(1:i)=1
          end if
       end do
       if(transformationCounter(this%numberOfTransformedCenters) .gt. this%numberOfIndividualTransformations) exit
    end do

    print *, ""
    write (*,'(A10,I10,A)') "Mixing ", this%numberOfDisplacedSystems, " HF calculations at different geometries"
    if(CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD.ne.0.0) &
         write (*,'(A10,I10,A,F18.12)') "Rejected ", this%numberOfRejectedSystems, &
         " geometries with energy higher than", this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD     
    print *, ""
    
    allocate(this%HFCoefficients(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    allocate(this%systemLabels(this%numberOfDisplacedSystems))
    call Matrix_constructor(this%configurationHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    !!Read HF energies and coefficients and fill CI matrix diagonals 
    ! minEnergy=0.0
    
    do sysI=1, this%numberOfDisplacedSystems
       write(wfnFile, '(I16)') sysI
       wfnUnit=300
       wfnFile="lowdin-"//trim(adjustl(wfnFile))//".wfn"
       open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
       call Vector_getFromFile(unit=wfnUnit, binary=.true., value=this%configurationHamiltonianMatrix%values(sysI,sysI), arguments=["TOTALENERGY"])

       do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
          arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
          arguments(1) = "COEFFICIENTS"
          this%HFCoefficients(sysI,speciesID) = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               unit=wfnUnit, binary=.true., arguments=arguments)
       end do
       
       close(wfnUnit)
       
       write(molFile, '(I16)') sysI
       molFile="LOWDIN-"//trim(adjustl(molFile))
       call MolecularSystem_loadFromFile('LOWDIN.SYS', molFile)
       write(this%systemLabels(sysI), '(A)') trim(MolecularSystem_instance%description)

       ! if(this%configurationHamiltonianMatrix%values(sysI,sysI) .lt. minEnergy) then
       ! print *, sysI, "newMinEnergy", this%configurationHamiltonianMatrix%values(sysI,sysI)
       ! minEnergy=this%configurationHamiltonianMatrix%values(sysI,sysI)
       ! end if
    end do
    ! print *, "minimum energy between configurations", minEnergy
    
!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for HF calculations at displaced geometries : ", timeB - timeA ," (s)"
    print *, ""
    
  end subroutine NonOrthogonalCI_displaceGeometries


  !>
  !! @brief Apply the transformation (translation or rotation) given by transformationCounter to each center, based in the originalMolecularSystemPositions the result is saved in molecularSystemInstance 
  !! @param this,transformationCounter,originalMolecularSystem
  !<
  subroutine NonOrthogonalCI_transformCoordinates(this,transformationCounter,originalMolecularSystem)
    type(NonOrthogonalCI) :: this
    integer :: transformationCounter(*)
    type(MolecularSystem) :: originalMolecularSystem

    real(8) :: centerX, centerY, centerZ, displacedOrigin(3), distanceToCenter
    integer :: center, displacementId
    real(8),allocatable :: X(:), Y(:), Z(:), W(:)
    integer :: i,j,k,p

    MolecularSystem_instance%description=""
    do i=1,this%numberOfTransformedCenters
       write(MolecularSystem_instance%description, '(A,I3)') trim(MolecularSystem_instance%description), transformationCounter(i)
    end do
    
    if(trim(this%transformationType).eq."TRANSLATION") then

       do center=1, this%numberOfTransformedCenters
          do p=1, size(originalMolecularSystem%allParticles)
             if(center.eq.originalMolecularSystem%allParticles(p)%particlePtr%translationCenter) then
                centerX=originalMolecularSystem%allParticles(p)%particlePtr%origin(1)
                centerY=originalMolecularSystem%allParticles(p)%particlePtr%origin(2)
                centerZ=originalMolecularSystem%allParticles(p)%particlePtr%origin(3)
             end if
          end do

          !!These loops update the molecular system file for each displaced geometry
          !!ADD DIFFERENT AXIS DISPLACEMENTS!
          displacementId=0
          do i=1,CONTROL_instance%TRANSLATION_SCAN_GRID(1)
             do j=1,CONTROL_instance%TRANSLATION_SCAN_GRID(2)
                do k=1,CONTROL_instance%TRANSLATION_SCAN_GRID(3)
                   displacementId=displacementId+1
                   if(displacementId .eq. transformationCounter(center) ) then
                      displacedOrigin(1)=centerX+CONTROL_instance%TRANSLATION_STEP*(i-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0)
                      displacedOrigin(2)=centerY+CONTROL_instance%TRANSLATION_STEP*(j-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0)
                      displacedOrigin(3)=centerZ+CONTROL_instance%TRANSLATION_STEP*(k-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0)

                      do p=1, size(MolecularSystem_instance%allParticles)
                         if(center.eq.MolecularSystem_instance%allParticles(p)%particlePtr%translationCenter) then
                            call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                         end if
                      end do

                      ! write(*, '(F4.1,A,F4.1,A,F4.1)') &
                      !      (i-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0)," ", &
                      !      (j-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0)," ", &
                      !      (k-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0)
                   end if
                end do
             end do
          end do
       end do
       
    else if(trim(this%transformationType).eq."ROTATION") then

       allocate(X(CONTROL_instance%ROTATIONAL_SCAN_GRID),&
            Y(CONTROL_instance%ROTATIONAL_SCAN_GRID),&
            Z(CONTROL_instance%ROTATIONAL_SCAN_GRID),&
            W(CONTROL_instance%ROTATIONAL_SCAN_GRID))

       call Lebedev_angularGrid(X(:),Y(:),Z(:),W(:),CONTROL_instance%ROTATIONAL_SCAN_GRID)

       do center=1, this%numberOfTransformedCenters
          displacementId=0
          do p=1, size(originalMolecularSystem%allParticles)
             if(this%rotationCenterList(p,1) .eq. center ) then
                centerX=originalMolecularSystem%allParticles(this%rotationCenterList(p,2))%particlePtr%origin(1)
                centerY=originalMolecularSystem%allParticles(this%rotationCenterList(p,2))%particlePtr%origin(2)
                centerZ=originalMolecularSystem%allParticles(this%rotationCenterList(p,2))%particlePtr%origin(3)
             end if
          end do
       
          do i=1,CONTROL_instance%ROTATIONAL_SCAN_GRID
             displacementId=displacementId+1
             if(displacementId .eq. transformationCounter(center) ) then
                do p=1, size(MolecularSystem_instance%allParticles)
                   if(this%rotationCenterList(p,1).eq. center ) then

                      distanceToCenter=sqrt((MolecularSystem_instance%allParticles(p)%particlePtr%origin(1)-centerX)**2 &
                           +(MolecularSystem_instance%allParticles(p)%particlePtr%origin(2)-centerY)**2 &
                           +(MolecularSystem_instance%allParticles(p)%particlePtr%origin(3)-centerZ)**2)

                      displacedOrigin(1)=centerX+X(i)*distanceToCenter
                      displacedOrigin(2)=centerY+Y(i)*distanceToCenter
                      displacedOrigin(3)=centerZ+Z(i)*distanceToCenter

                      call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                   end if
                end do
             end if
          end do
       end do
    end if
                   
  end subroutine NonOrthogonalCI_transformCoordinates
  !>
  !! @brief Run a Hartree-Fock calculation at a displaced geometry
  !!
  !! @param testEnergy -> HF energy obtained
  !<
  subroutine NonOrthogonalCI_runHF(testEnergy)
    implicit none
    real(8) :: testEnergy
    character(50) :: wfnFile
    integer :: wfnUnit
    !Do SCF
    select case ( trim(CONTROL_instance%METHOD) )

       ! case('MM')
       !    call system("lowdin-MolecularMechanics.x CONTROL_instance%FORCE_FIELD")
       ! case('RHF')
       !    call system("lowdin-SCF.x RHF")
    case('UHF')
       call system("lowdin-SCF.x UHF >> extendedOutput.eut")
       ! case('RKS')
       !    call system("lowdin-SCF.x RKS")
       ! case('UKS')
       !    call system("lowdin-SCF.x UKS")
    case default
       call Solver_exception(ERROR, "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented for non orthogonal CI", &
            "At NOCI module ")
    end select
    
    !!Read HF energy
    wfnUnit=300
    wfnFile="lowdin.wfn"
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=testEnergy, arguments=["TOTALENERGY"])
    close(300)

  end subroutine NonOrthogonalCI_runHF

  !>
  !! @brief Saves molecular system and wfn files for a displaced system 
  !!
  !! @param systemID
  !<
  subroutine NonOrthogonalCI_saveSystem(systemID)
    implicit none
    integer :: systemID
    character(50) :: wfnFile, molFile, auxString
    integer :: wfnUnit
    !!Save the molecular system and WF results to separate files for each geometry 
    write(auxString, '(I16)') systemID
    
    molFile="LOWDIN-"//trim(adjustl(auxString))
    wfnFile="lowdin-"//trim(adjustl(auxString))//".wfn"

    ! print *, "molFile", molFile, "wfnFile", wfnFile

    call MolecularSystem_saveToFile(molFile)
    call system("cp lowdin.wfn "//wfnFile)
    
  end subroutine NonOrthogonalCI_saveSystem

  !>
  !! @brief Computes overlap matrix element between two configurations along with one particle energy contributions
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_configurationsOverlapAndOneParticle(this)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: sysI,sysII
    type(Matrix), allocatable :: mergedCoefficients(:), molecularHCoreMatrix(:)
    integer :: integralsUnit
    character(50) :: integralsFile    
    character(50) :: arguments(2)
    integer :: speciesID
    integer :: a,b,bb,mu,nu    
    type(Matrix) :: auxMatrix, auxOverlapMatrix, auxKineticMatrix, auxAttractionMatrix
    type(Matrix) :: molecularOverlapMatrix
    type(Vector) :: overlapDeterminant
    real(8) :: oneParticleEnergy
    integer :: skippedElements
    
    real(8) :: timeA
    real(8) :: timeB

    call Matrix_constructor(this%configurationOverlapMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 1.0_8)
    
    allocate(this%inverseOverlapMatrices(this%numberOfDisplacedSystems,this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    allocate(this%HCoreMatrices(this%numberOfDisplacedSystems,this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    
    print *, ""
    print *, "Overlap matrix elements are computed for all configurations pairs"
    write (*,'(A,ES8.1)') " Hamiltonian matrix elements are computed for pairs with overlap higher than",&
         CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD
    print *, "For pairs with lower overlap, aproximating H(I,II)=S(I,II)*[E(I)+E(II)]/2"    
    print *, ""

    allocate(mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies))
    allocate(molecularHCoreMatrix(molecularSystem_instance%numberOfQuantumSpecies))
    call Vector_constructor(overlapDeterminant, molecularSystem_instance%numberOfQuantumSpecies, 0.0_8)        

    skippedElements=0
    !$  timeA = omp_get_wtime()
    do sysI=1, this%numberOfDisplacedSystems       
       do sysII=sysI+1, this%numberOfDisplacedSystems !sysI+1

          !This generates a new molecular system
          call NonOrthogonalCI_mergeBasisAndCoefficients(this,sysI,sysII,mergedCoefficients)

          !!!!Overlap
          integralsUnit = 30
          integralsFile = "lowdin.opints"
          
          !! Calculate one- particle integrals  
          call system("lowdin-ints.x ONE_PARTICLE >> extendedOutput.eut")

          open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

          do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

             arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
             arguments(1) = "OVERLAP"

             auxOverlapMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  unit=integralsUnit, binary=.true., arguments=arguments)

             !!Test 

             ! print *, "auxOverlapMatrix", speciesID
             ! call Matrix_show(auxOverlapMatrix)

             call Matrix_constructor(molecularOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID),8), 0.0_8 )

             do a=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
                do b=MolecularSystem_getOcupationNumber(speciesID)+1, MolecularSystem_getOcupationNumber(speciesID)*2 !sysII
                   bb=b-MolecularSystem_getOcupationNumber(speciesID)
                   do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
                      do nu= MolecularSystem_getTotalNumberOfContractions(speciesID)/2+1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                         ! print *, "overlap", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
                         !      auxOverlapMatrix%values(mu,nu)

                         molecularOverlapMatrix%values(a,bb)=molecularOverlapMatrix%values(a,bb)+&
                              mergedCoefficients(speciesID)%values(mu,a)*&
                              mergedCoefficients(speciesID)%values(nu,b)*&
                              auxOverlapMatrix%values(mu,nu)

                      end do
                   end do
                end do
             end do

             ! print *, "molecularOverlapMatrix", speciesID
             ! call Matrix_show(molecularOverlapMatrix)

             this%inverseOverlapMatrices(sysI,sysII,speciesID)=Matrix_inverse(molecularOverlapMatrix)

             ! print *, "this%inverseOverlapMatrices", speciesID
             ! call Matrix_show(this%inverseOverlapMatrices(speciesID))

             call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant%values(speciesID),method="LU")             

             ! print *, "OverlapDeterminantLU", speciesID, overlapDeterminant%values(speciesID)

             this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*overlapDeterminant%values(speciesID)

          end do

          close(integralsUnit)

          this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)
          !! SKIP ENERGY EVALUATION IF OVERLAP IS TOO LOW

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*&
                  (this%configurationHamiltonianMatrix%values(sysI,sysI)+this%configurationHamiltonianMatrix%values(sysII,sysII))/2.0
             this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)             
             skippedElements=skippedElements+1
             cycle
          end if

          print *, "Overlap determinant product for", sysI, sysII, this%configurationOverlapMatrix%values(sysI,sysII)           
          
          open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

          !!!!HCore
          do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

             ! print *, "coefficientsI", speciesID
             ! call Matrix_show(coefficientsI(speciesID))

             ! print *, "coefficientsII", speciesID
             ! call Matrix_show(coefficientsII(speciesID))

             arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))

             arguments(1) = "KINETIC"

             auxKineticMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  unit=integralsUnit, binary=.true., arguments=arguments)

             arguments(1) = "ATTRACTION"

             auxAttractionMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
                  unit=integralsUnit, binary=.true., arguments=arguments)

             !!Test 
             call Matrix_constructor(molecularHCoreMatrix(speciesID), int(MolecularSystem_getOcupationNumber(speciesID),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID),8), 0.0_8 )

             ! print *, "auxKineticMatrix", speciesID
             ! call Matrix_show(auxKineticMatrix)
             ! print *, "auxAttractionMatrix", speciesID
             ! call Matrix_show(auxAttractionMatrix)

             do a=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
                do b=MolecularSystem_getOcupationNumber(speciesID)+1, MolecularSystem_getOcupationNumber(speciesID)*2 !sysII
                   bb=b-MolecularSystem_getOcupationNumber(speciesID)
                   do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
                      do nu= MolecularSystem_getTotalNumberOfContractions(speciesID)/2+1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                         ! print *, "hcore", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
                         !      auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
                         !      auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID))

                         molecularHCoreMatrix(speciesID)%values(a,bb)=molecularHCoreMatrix(speciesID)%values(a,bb)+&
                              mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                              (auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
                              auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID)))                         
                      end do
                   end do
                end do
             end do

             ! print *, "molecularHCoreMatrix", speciesID
             ! call Matrix_show(molecularHCoreMatrix(speciesID))
             !!End test                          
          end do

          close(integralsUnit)


          !!Point charge-Point charge repulsion
          this%configurationHamiltonianMatrix%values(sysI,sysII)=MolecularSystem_getPointChargesEnergy()
          ! print *, "Point charge-Point charge repulsion", MolecularSystem_getPointChargesEnergy()

          !!One Particle Terms
          do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
             oneParticleEnergy=0.0
             do a=1, int(MolecularSystem_getOcupationNumber(speciesID),8) !sysI
                do b=1, int(MolecularSystem_getOcupationNumber(speciesID),8) !sysII
                   oneParticleEnergy=oneParticleEnergy+ molecularHCoreMatrix(speciesID)%values(a,b)*&
                        this%inverseOverlapMatrices(sysI,sysII,speciesID)%values(b,a)
                end do
             end do
             this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+oneParticleEnergy
             ! print *, "oneParticleEnergy for species", i, oneParticleEnergy
          end do
          
       end do
    end do

    print *, ""
    print *, "Configuration pairs skipped by overlap threshold: ", skippedElements
    print *, "Four center integrals will be computed for", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2-skippedElements, "configuration pairs"
    print *, ""
    !$  timeB = omp_get_wtime()
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Overlap and one particle integrals : ", timeB-timeA ," (s)"
    print *, ""
    
    deallocate(mergedCoefficients,molecularHCoreMatrix)

  end subroutine NonOrthogonalCI_configurationsOverlapAndOneParticle

  !>
  !! @brief Generates a new molecular system from the combination of sysI and sysII geometries along with their combined coefficients
  !!
  !! @param sysI and sysII: molecular system indices. Merged Coefficients: Matrices for output
  !<
  subroutine NonOrthogonalCI_mergeBasisAndCoefficients(this,sysI,sysII,mergedCoefficients)
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII !Indices of the systems to merge
    type(Matrix) :: mergedCoefficients(*)
    
    type(MolecularSystem), target :: molecularSystemI, molecularSystemII
    character(50) :: molFileI, molFileII, wfnFile
    character(50) :: arguments(2)
    integer :: wfnUnit
    integer :: speciesID, i, j, mu, nu
    type(Matrix) :: auxMatrix
    type(Vector) :: auxVector
    
    !!Read molecular system I
    write(molFileI, '(I16)') sysI
    molFileI="LOWDIN-"//trim(adjustl(molFileI))
    call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileI)
    call MolecularSystem_copyConstructor(molecularSystemI, molecularSystem_instance)

    !!Read molecular system II
    write(molFileII, '(I16)') sysII
    molFileII="LOWDIN-"//trim(adjustl(molFileII))
    call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileII)
    call MolecularSystem_copyConstructor(molecularSystemII, molecularSystem_instance)

    ! print *, "Merging systems from geometries ", sysI, sysII
    call MolecularSystem_mergeTwoSystems(molecularSystem_instance, molecularSystemI, molecularSystemII)

    ! call MolecularSystem_showInformation()  
    ! call MolecularSystem_showParticlesInformation()
    ! call MolecularSystem_showCartesianMatrix()

    call MolecularSystem_saveToFile()

    !!***********************************************************************************************************
    !! Mix coefficients of occupied orbitals of both systems and create a dummy density matrix to lowdin.wfn file
    wfnUnit = 500
    wfnFile = "lowdin.wfn"
    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

       arguments(1) = "COEFFICIENTS"
       call Matrix_constructor(mergedCoefficients(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
            int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8 )

       do i=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
          do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
             mergedCoefficients(speciesID)%values(mu,i)=this%HFCoefficients(sysI,speciesID)%values(mu,i)
          end do
       end do

       do i=1, MolecularSystem_getOcupationNumber(speciesID) !sysII
          do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2  !sysII
             j=MolecularSystem_getOcupationNumber(speciesID)+i
             nu=MolecularSystem_getTotalNumberOfContractions(speciesID)/2+mu
             mergedCoefficients(speciesID)%values(nu,j)=this%HFCoefficients(sysII,speciesID)%values(mu,i)
          end do
       end do

       ! print *, "Merged coefficients matrix for ", speciesID
       ! call Matrix_show(mergedCoefficients(speciesID))

       call Matrix_writeToFile(mergedCoefficients(speciesID), unit=wfnUnit, binary=.true., arguments = arguments )

       arguments(1) = "DENSITY"
       call Matrix_constructor(auxMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
            int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 1.0_8 )

       ! do i = 1 , MolecularSystem_getOcupationNumber(speciesID)*2 !!double size A+B
       !    do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
       !       do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
       !          auxMatrix%values(mu,nu)=auxMatrix%values(mu,nu)&
       !               +MolecularSystem_getEta(speciesID)*mergedCoefficients(speciesID)%values(mu,i)*mergedCoefficients(speciesID)%values(nu,i)
       !       end do
       !    end do
       ! end do

       ! print *, "auxDensity", speciesID
       ! call Matrix_show(auxMatrix)

       call Matrix_writeToFile(auxMatrix, unit=wfnUnit, binary=.true., arguments = arguments )

       arguments(1) = "ORBITALS"
       call Vector_constructor(auxVector, MolecularSystem_getTotalNumberOfContractions(speciesID), 0.0_8 )

       call Vector_writeToFile(auxVector, unit=wfnUnit, binary=.true., arguments = arguments )

       ! Only 2*occupied orbitals are going to be transformed
       arguments(1) = "REMOVED-ORBITALS"
       call Vector_writeToFile(unit=wfnUnit, binary=.true., &
            value=real(MolecularSystem_getTotalNumberOfContractions(speciesID)-MolecularSystem_getOcupationNumber(speciesID)*2,8),&
            arguments= arguments )

    end do
    close(wfnUnit)

  end subroutine NonOrthogonalCI_mergeBasisAndCoefficients

  subroutine NonOrthogonalCI_firstImplementation(this)
    implicit none
    type(NonOrthogonalCI) :: this

    integer :: a,b,bb,c,d,dd,g,i,j,k,l,p,mu,nu,sysI,sysII,speciesID
    character(50) :: molFileI, wfnFileI, molFileII, wfnFileII, auxString
    character(50) :: wfnFile, integralsFile, outFile
    integer :: wfnUnit,wfnUnitI,wfnUnitII,integralsUnit,outUnit
    type(MolecularSystem), target :: molecularSystemI, molecularSystemII
    character(50) :: arguments(2)
    type(Vector) :: auxVector
    type(Vector) :: overlapDeterminant
    type(Matrix) :: auxMatrix, auxOverlapMatrix, auxKineticMatrix, auxAttractionMatrix
    type(Matrix) :: atomicOverlapMatrix, atomicHCoreMatrix, atomicKineticMatrix, atomicAttractionMatrix, molecularOverlapMatrix
    type(Matrix), allocatable :: mergedCoefficients(:), molecularHCoreMatrix(:)
    logical :: existFile
    real(8) :: oneParticleEnergy
    type(matrix), allocatable :: fourCenterIntegrals(:,:)
    type(imatrix), allocatable :: twoIndexArray(:),fourIndexArray(:)
    integer :: ssize1, auxIndex, auxIndex1
    real(8) :: interactionEnergy
    real(8) :: timeA
    real(8) :: timeB
    real(8) :: timeSetup, timeOverlap, timeHCore, timeTwoIntegrals, timeIntegralTransformation, timeMatrixElement

    timeSetup=0.0
    timeOverlap=0.0
    timeHCore=0.0
    timeTwoIntegrals=0.0
    timeMatrixElement=0.0
    
    allocate(mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies))
    ! allocate(molecularHCoreMatrix(molecularSystem_instance%numberOfQuantumSpecies))

    ! call Vector_constructor(overlapDeterminant, molecularSystem_instance%numberOfQuantumSpecies, 0.0_8)        

    do sysI=1, this%numberOfDisplacedSystems       
       ! !!Read molecular system I
       ! write(auxString, '(I16)') sysI
       ! molFileI="LOWDIN-"//trim(adjustl(auxString))
       ! call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileI)
       ! call MolecularSystem_copyConstructor(molecularSystemI, molecularSystem_instance)

       ! write(this%systemLabels(sysI), '(A)') trim(MolecularSystem_instance%description)

       do sysII=sysI+1, this%numberOfDisplacedSystems !sysI+1

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) cycle
          
          !This generates a new molecular system
          call NonOrthogonalCI_mergeBasisAndCoefficients(this,sysI,sysII,mergedCoefficients)

          ! !!Read molecular system II
          ! write(auxString, '(I16)') sysII
          ! molFileII="LOWDIN-"//trim(adjustl(auxString))
          ! call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileII)
          ! call MolecularSystem_copyConstructor(molecularSystemII, molecularSystem_instance)

          ! print *, "Merging systems from geometries ", sysI, sysII
          ! call MolecularSystem_mergeTwoSystems(molecularSystem_instance, molecularSystemI, molecularSystemII)

          ! ! call MolecularSystem_showInformation()  
          ! ! call MolecularSystem_showParticlesInformation()
          ! ! call MolecularSystem_showCartesianMatrix()

          ! call MolecularSystem_saveToFile()

          ! !!***********************************************************************************************************
          ! !! Mix coefficients of occupied orbitals of both systems and create a dummy density matrix to lowdin.wfn file
          ! wfnUnit = 500
          ! wfnFile = "lowdin.wfn"
          ! open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
          ! do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

          !    arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

          !    arguments(1) = "COEFFICIENTS"
          !    call Matrix_constructor(mergedCoefficients(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
          !         int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8 )

          !    do i=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
          !       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
          !          mergedCoefficients(speciesID)%values(mu,i)=this%HFCoefficients(sysI,speciesID)%values(mu,i)
          !       end do
          !    end do

          !    do i=1, MolecularSystem_getOcupationNumber(speciesID) !sysII
          !       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2  !sysII
          !          j=MolecularSystem_getOcupationNumber(speciesID)+i
          !          nu=MolecularSystem_getTotalNumberOfContractions(speciesID)/2+mu
          !          mergedCoefficients(speciesID)%values(nu,j)=this%HFCoefficients(sysII,speciesID)%values(mu,i)
          !       end do
          !    end do

          !    ! print *, "Merged coefficients matrix for ", speciesID
          !    ! call Matrix_show(mergedCoefficients(speciesID))

          !    call Matrix_writeToFile(mergedCoefficients(speciesID), unit=wfnUnit, binary=.true., arguments = arguments )

          !    arguments(1) = "DENSITY"
          !    call Matrix_constructor(auxMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
          !         int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 1.0_8 )

          !    ! do i = 1 , MolecularSystem_getOcupationNumber(speciesID)*2 !!double size A+B
          !    !    do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !    !       do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !    !          auxMatrix%values(mu,nu)=auxMatrix%values(mu,nu)&
          !    !               +MolecularSystem_getEta(speciesID)*mergedCoefficients(speciesID)%values(mu,i)*mergedCoefficients(speciesID)%values(nu,i)
          !    !       end do
          !    !    end do
          !    ! end do

          !    ! print *, "auxDensity", speciesID
          !    ! call Matrix_show(auxMatrix)

          !    call Matrix_writeToFile(auxMatrix, unit=wfnUnit, binary=.true., arguments = arguments )

          !    arguments(1) = "ORBITALS"
          !    call Vector_constructor(auxVector, MolecularSystem_getTotalNumberOfContractions(speciesID), 0.0_8 )

          !    call Vector_writeToFile(auxVector, unit=wfnUnit, binary=.true., arguments = arguments )

          !    ! Only 2*occupied orbitals are going to be transformed
          !    arguments(1) = "REMOVED-ORBITALS"
          !    call Vector_writeToFile(unit=wfnUnit, binary=.true., &
          !         value=real(MolecularSystem_getTotalNumberOfContractions(speciesID)-MolecularSystem_getOcupationNumber(speciesID)*2,8),&
          !         arguments= arguments )

          ! end do
          ! close(wfnUnit)

          ! integralsUnit = 30
          ! integralsFile = "lowdin.opints"

          ! !$  timeB = omp_get_wtime()

          ! !$  timeSetup=timeSetup+(timeB - timeA)

          ! !$  timeA = omp_get_wtime()
          
          ! !! Calculate one- particle integrals  
          ! call system("lowdin-ints.x ONE_PARTICLE >> extendedOutput.eut")

          ! open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

          ! do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

          !    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
          !    arguments(1) = "OVERLAP"

          !    auxOverlapMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         unit=integralsUnit, binary=.true., arguments=arguments)

          !    !!Test 

          !    ! print *, "auxOverlapMatrix", speciesID
          !    ! call Matrix_show(auxOverlapMatrix)

          !    call Matrix_constructor(molecularOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID),8), &
          !         int(MolecularSystem_getOcupationNumber(speciesID),8), 0.0_8 )

          !    do a=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
          !       do b=MolecularSystem_getOcupationNumber(speciesID)+1, MolecularSystem_getOcupationNumber(speciesID)*2 !sysII
          !          bb=b-MolecularSystem_getOcupationNumber(speciesID)
          !          do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
          !             do nu= MolecularSystem_getTotalNumberOfContractions(speciesID)/2+1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
          !                ! print *, "overlap", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
          !                !      auxOverlapMatrix%values(mu,nu)

          !                molecularOverlapMatrix%values(a,bb)=molecularOverlapMatrix%values(a,bb)+&
          !                     mergedCoefficients(speciesID)%values(mu,a)*&
          !                     mergedCoefficients(speciesID)%values(nu,b)*&
          !                     auxOverlapMatrix%values(mu,nu)

          !             end do
          !          end do
          !       end do
          !    end do

          !    ! print *, "molecularOverlapMatrix", speciesID
          !    ! call Matrix_show(molecularOverlapMatrix)

          !    this%inverseOverlapMatrices(sysI,sysII,speciesID)=Matrix_inverse(molecularOverlapMatrix)

          !    ! print *, "this%inverseOverlapMatrices", speciesID
          !    ! call Matrix_show(this%inverseOverlapMatrices(speciesID))

          !    call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant%values(speciesID),method="LU")             

          !    ! print *, "OverlapDeterminantLU", overlapDeterminant%values(speciesID)

          !    this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*overlapDeterminant%values(speciesID)

          ! end do

          ! close(integralsUnit)

          ! !$  timeB = omp_get_wtime()

          ! !$  timeOverlap=timeOverlap+(timeB - timeA)

          ! !$  timeA = omp_get_wtime()

          ! !! SKIP ENERGY EVALUATION IF OVERLAP IS TOO LOW
          ! print *, "Overlap determinant product for", sysI, sysII, this%configurationOverlapMatrix%values(sysI,sysII)           
          ! this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)

          ! if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
          !      abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
          !    this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*&
          !         (this%configurationHamiltonianMatrix%values(sysI,sysI)+this%configurationHamiltonianMatrix%values(sysII,sysII))/2.0
          !    this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)             
          !    cycle
          ! end if

          ! open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

          ! do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

          !    ! print *, "coefficientsI", speciesID
          !    ! call Matrix_show(coefficientsI(speciesID))

          !    ! print *, "coefficientsII", speciesID
          !    ! call Matrix_show(coefficientsII(speciesID))

          !    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))

          !    arguments(1) = "KINETIC"

          !    auxKineticMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         unit=integralsUnit, binary=.true., arguments=arguments)

          !    arguments(1) = "ATTRACTION"

          !    auxAttractionMatrix = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
          !         unit=integralsUnit, binary=.true., arguments=arguments)

          !    !!Test 
          !    call Matrix_constructor(molecularHCoreMatrix(speciesID), int(MolecularSystem_getOcupationNumber(speciesID),8), &
          !         int(MolecularSystem_getOcupationNumber(speciesID),8), 0.0_8 )

          !    ! print *, "auxKineticMatrix", speciesID
          !    ! call Matrix_show(auxKineticMatrix)
          !    ! print *, "auxAttractionMatrix", speciesID
          !    ! call Matrix_show(auxAttractionMatrix)

          !    do a=1, MolecularSystem_getOcupationNumber(speciesID) !sysI
          !       do b=MolecularSystem_getOcupationNumber(speciesID)+1, MolecularSystem_getOcupationNumber(speciesID)*2 !sysII
          !          bb=b-MolecularSystem_getOcupationNumber(speciesID)
          !          do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)/2 !sysI
          !             do nu= MolecularSystem_getTotalNumberOfContractions(speciesID)/2+1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
          !                ! print *, "hcore", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
          !                !      auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
          !                !      auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID))

          !                molecularHCoreMatrix(speciesID)%values(a,bb)=molecularHCoreMatrix(speciesID)%values(a,bb)+&
          !                     mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
          !                     (auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
          !                     auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID)))                         
          !             end do
          !          end do
          !       end do
          !    end do

          !    ! print *, "molecularHCoreMatrix", speciesID
          !    ! call Matrix_show(molecularHCoreMatrix(speciesID))
          !    !!End test                          
          ! end do

          ! close(integralsUnit)

          ! !$  timeB = omp_get_wtime()

          ! !$  timeHCore=timeHCore+(timeB - timeA)

          !$  timeA = omp_get_wtime()
          
          !! Calculate two- particle integrals  
          call system("lowdin-ints.x TWO_PARTICLE_R12 >> extendedOutput.eut")

          !$  timeB = omp_get_wtime()

          !$  timeTwoIntegrals=timeTwoIntegrals+(timeB - timeA)

          !$  timeA = omp_get_wtime()

          !! Transform and load integrals
          call system("lowdin-integralsTransformation.x >> extendedOutput.eut")

          allocate(fourCenterIntegrals(MolecularSystem_instance%numberOfQuantumSpecies,MolecularSystem_instance%numberOfQuantumSpecies), &
               twoIndexArray(MolecularSystem_instance%numberOfQuantumSpecies), &
               fourIndexArray(MolecularSystem_instance%numberOfQuantumSpecies))

          do i=1, MolecularSystem_instance%numberOfQuantumSpecies

             ! print *, "reading integrals species", i
             !!Two particle integrals indexes
             call Matrix_constructorInteger(twoIndexArray(i), int(MolecularSystem_getTotalNumberOfContractions(i),8), &
                  int(MolecularSystem_getTotalNumberOfContractions(i),8) , 0 )

             c = 0
             do a=1,MolecularSystem_getTotalNumberOfContractions(i)
                do b=a, MolecularSystem_getTotalNumberOfContractions(i)
                   c = c + 1
                   twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
                   twoIndexArray(i)%values(b,a) = twoIndexArray(i)%values(a,b)
                end do
             end do

             ssize1 = MolecularSystem_getTotalNumberOfContractions( i )
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

             call ReadTransformedIntegrals_readOneSpecies( i , fourCenterIntegrals(i,i)   )
             fourCenterIntegrals(i,i)%values = &
                  fourCenterIntegrals(i,i)%values * MolecularSystem_getCharge(i)**2.0

             ! print *, "transformed integrals for species", i
             ! do a=1, MolecularSystem_getOcupationNumber(i)*2
             !    do b=1, MolecularSystem_getOcupationNumber(i)*2
             !       do c=1, MolecularSystem_getOcupationNumber(i)*2
             !          do d=1, MolecularSystem_getOcupationNumber(i)*2
             !             auxIndex = fourIndexArray(i)%values( &
             !                  twoIndexArray(i)%values(a,b), &
             !                  twoIndexArray(i)%values(c,d) )
             !             ! print *, a, b, c, d, twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d), fourIndexArray(i)%values( &
             !             !      twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d)), fourCenterIntegrals(i,i)%values(auxIndex, 1)
             !          end do
             !       end do
             !    end do
             ! end do
          end do
          
          do i=1, MolecularSystem_instance%numberOfQuantumSpecies-1
             do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
                ! print *, "reading integrals species", i, j

                call ReadTransformedIntegrals_readTwoSpecies( i, j, fourCenterIntegrals(i,j) )

                fourCenterIntegrals(i,j)%values = &
                     fourCenterIntegrals(i,j)%values * MolecularSystem_getCharge(i) * MolecularSystem_getCharge(j)

                ! ssize1 = MolecularSystem_getTotalNumberOfContractions( j )
                ! ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2

                ! print *, "transformed integrals for species", i, j
                ! do a=1, MolecularSystem_getOcupationNumber(i)*2
                !    do b=1, MolecularSystem_getOcupationNumber(i)*2
                !       auxIndex1 = ssize1 * (twoIndexArray(i)%values(a,b) - 1 ) 
                !       do c=1, MolecularSystem_getOcupationNumber(j)*2
                !          do d=1, MolecularSystem_getOcupationNumber(j)*2
                !             auxIndex = auxIndex1  + twoIndexArray(j)%values(c,d) 
                !             print *, a, b, c, d, auxIndex1, twoIndexArray(j)%values(c,d), auxIndex, fourCenterIntegrals(i,j)%values(auxIndex, 1)
                !          end do
                !       end do
                !    end do
                ! end do

             end do

          end do

          !$  timeB = omp_get_wtime()

          !$  timeIntegralTransformation=timeIntegralTransformation+(timeB - timeA)

          !$  timeA = omp_get_wtime()

          
!!!Compute Hamiltonian Matrix element between displaced geometries

          ! !!Point charge-Point charge repulsion
          ! this%configurationHamiltonianMatrix%values(sysI,sysII)=MolecularSystem_getPointChargesEnergy()
          ! ! print *, "Point charge-Point charge repulsion", MolecularSystem_getPointChargesEnergy()

          ! !!One Particle Terms
          ! do i=1, MolecularSystem_instance%numberOfQuantumSpecies
          !    oneParticleEnergy=0.0
          !    do a=1, int(MolecularSystem_getOcupationNumber(i),8) !sysI
          !       do b=1, int(MolecularSystem_getOcupationNumber(i),8) !sysII
          !          oneParticleEnergy=oneParticleEnergy+ molecularHCoreMatrix(i)%values(a,b)*&
          !               this%inverseOverlapMatrices(sysI,sysII,i)%values(b,a)
          !       end do
          !    end do
          !    this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+oneParticleEnergy
          !    ! print *, "oneParticleEnergy for species", i, oneParticleEnergy
          ! end do


          !!Same species repulsion
          do i=1, MolecularSystem_instance%numberOfQuantumSpecies
             interactionEnergy=0.0
             do a=1, MolecularSystem_getOcupationNumber(i) !sysI
                do b=MolecularSystem_getOcupationNumber(i)+1, MolecularSystem_getOcupationNumber(i)*2 !sysII
                   bb=b-MolecularSystem_getOcupationNumber(i)
                   do c=1, MolecularSystem_getOcupationNumber(i) !sysI
                      do d=MolecularSystem_getOcupationNumber(i)+1, MolecularSystem_getOcupationNumber(i)*2 !sysII
                         dd=d-MolecularSystem_getOcupationNumber(i)
                         auxIndex = fourIndexArray(i)%values(twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d) )
                         interactionEnergy=interactionEnergy+0.5*fourCenterIntegrals(i,i)%values(auxIndex, 1)*&
                              (this%inverseOverlapMatrices(sysI,sysII,i)%values(bb,a)*this%inverseOverlapMatrices(sysI,sysII,i)%values(dd,c)-&
                              this%inverseOverlapMatrices(sysI,sysII,i)%values(dd,a)*this%inverseOverlapMatrices(sysI,sysII,i)%values(bb,c))
                         ! print *, a, b, c, d, twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d), fourIndexArray(i)%values( &
                         !      twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d)), 
                      end do
                   end do
                end do
             end do
             this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+interactionEnergy
             ! print *, "same species interactionEnergy for species", i, interactionEnergy
          end do

          !!Interspecies repulsion
          do i=1, MolecularSystem_instance%numberOfQuantumSpecies-1
             do j=i+1, MolecularSystem_instance%numberOfQuantumSpecies
                interactionEnergy=0.0
                ssize1 = MolecularSystem_getTotalNumberOfContractions( j )
                ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2
                do a=1, MolecularSystem_getOcupationNumber(i) !sysI
                   do b=MolecularSystem_getOcupationNumber(i)+1, MolecularSystem_getOcupationNumber(i)*2  !sysII
                      bb=b-MolecularSystem_getOcupationNumber(i)
                      auxIndex1 = ssize1 * (twoIndexArray(i)%values(a,b) - 1 ) 
                      do c=1, MolecularSystem_getOcupationNumber(j)  !sysI
                         do d=MolecularSystem_getOcupationNumber(j)+1, MolecularSystem_getOcupationNumber(j)*2  !sysII
                            dd=d-MolecularSystem_getOcupationNumber(j)
                            auxIndex = auxIndex1  + twoIndexArray(j)%values(c,d) 
                            interactionEnergy=interactionEnergy+fourCenterIntegrals(i,j)%values(auxIndex, 1)*&
                                 this%inverseOverlapMatrices(sysI,sysII,i)%values(bb,a)*this%inverseOverlapMatrices(sysI,sysII,j)%values(dd,c)
                            ! print *, a, b, c, d,  fourCenterIntegrals(i,j)%values(auxIndex, 1), this%inverseOverlapMatrices(i)%values(bb,a), this%inverseOverlapMatrices(j)%values(dd,c)
                         end do
                      end do
                   end do
                end do
                this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+interactionEnergy
                ! print *, "interspecies interactionEnergy for species", i, j, interactionEnergy
             end do
          end do

          deallocate(fourCenterIntegrals,twoIndexArray,fourIndexArray)

          print *, "Unscaled and overlap scaled hamiltonian element for", sysI,sysII, this%configurationHamiltonianMatrix%values(sysI,sysII), &
               this%configurationHamiltonianMatrix%values(sysI,sysII)*this%configurationOverlapMatrix%values(sysI,sysII)

          this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)*&
                  this%configurationOverlapMatrix%values(sysI,sysII)

          this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)
          
          !$  timeB = omp_get_wtime()

          !$  timeMatrixElement=timeMatrixElement+(timeB - timeA)

          !$  timeA = omp_get_wtime()


       end do !Displacements loops
    end do !Displacements loops

    print *, ""
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for atomic four index integrals : ", timeTwoIntegrals ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for four index integrals transformation : ", timeIntegralTransformation ," (s)"
    print *, ""
    
  end subroutine NonOrthogonalCI_firstImplementation

  subroutine NonOrthogonalCI_diagonalizeCImatrix(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(Matrix) :: transformationMatrix,eigenVectors,auxMatrix
    type(Vector) :: eigenValues
    integer :: removedStates
    integer :: i,j
    real(8) :: timeA
    real(8) :: timeB
    
    !$  timeA = omp_get_wtime()

    call Matrix_constructor(this%configurationCoefficients, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Vector_constructor(this%statesEigenvalues, this%numberOfDisplacedSystems, 0.0_8)

    ! print *, "non orthogonal CI overlap Matrix "
    ! call Matrix_show(this%configurationOverlapMatrix)

    ! print *, "non orthogonal CI Hamiltionian Matrix "
    ! call Matrix_show(this%configurationHamiltonianMatrix)

    print *, ""
    print *, "Transforming non orthogonal CI Hamiltionian Matrix..."
    
    call Matrix_constructor(transformationMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8) , 0.0_8)

    call Vector_constructor( eigenValues, this%numberOfDisplacedSystems )
    call Matrix_constructor( eigenVectors,int(this%numberOfDisplacedSystems,8),int(this%numberOfDisplacedSystems,8))

    !!****************************************************************
    !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
    !!          
    call Matrix_eigen( this%configurationOverlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )

    !! Remove states from configurations with linear dependencies
    if ( CONTROL_instance%OVERLAP_EIGEN_THRESHOLD .gt. 0.0 ) then

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

    end if

    !!Ortogonalizacion simetrica
    transformationMatrix%values  = &
         matmul(transformationMatrix%values, transpose(eigenVectors%values))

    call Vector_destructor( eigenValues )
    call Matrix_destructor( eigenVectors )

    ! print *,"Matriz de transformacion "
    ! call Matrix_show( transformationMatrix )

    !!**********************************************************************************************
    !! Transform configuration hamiltonian matrix
    !!
    this%configurationHamiltonianMatrix%values = &
         matmul( matmul( transpose( transformationMatrix%values ) , &
         this%configurationHamiltonianMatrix%values), transformationMatrix%values )

    ! print *,"transformed Hamiltonian Matrix "
    ! call Matrix_show( this%configurationHamiltonianMatrix )

    print *, "Diagonalizing non orthogonal CI Hamiltionian Matrix..."
    !! Calcula valores y vectores propios de matriz de CI transformada.
    call Matrix_eigen( this%configurationHamiltonianMatrix, this%statesEigenvalues, this%configurationCoefficients, SYMMETRIC )

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
    do i = 1, size(this%statesEigenvalues%values)
       write (*,"(T8,A17,I3,A10, F18.12)") "STATE: ", i, " ENERGY = ", this%statesEigenvalues%values(i)
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

    !$  timeB = omp_get_wtime()
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI matrix diagonalization : ", timeB - timeA ," (s)"
    
  end subroutine NonOrthogonalCI_diagonalizeCImatrix

  subroutine NonOrthogonalCI_plotDensities(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(Matrix) :: plotPoints, auxMatrix
    type(Matrix), allocatable :: orbitalsInGrid(:,:)
    type(Vector), allocatable :: densityInGrid(:)
    integer :: gridSize, state
    real(8) :: densityIntegral, auxValue
    integer :: a,b,g,i,j,k,p, sysI, sysII, speciesID
    character(50) :: molFileI, outFile, auxString
    integer :: wfnUnit,wfnUnitI,wfnUnitII,integralsUnit,outUnit
    real(8) :: timeA
    real(8) :: timeB
    
!$  timeA = omp_get_wtime()

    write(*,"(A)") ""
    write(*,"(A)") " DENSITY PLOTS"
    write(*,"(A)") "========================================="
    write(*,"(A)") ""

    !Grid for density plots

    gridSize=1001
    call Matrix_constructor( plotPoints,int(gridSize**2,8),3_8)
    p=0
    do i=0,gridSize-1
       do j=0,0!gridSize-1
          do k=0,gridSize-1
             p=p+1
             plotPoints%values(p,1)=0+(-500.0+i)/250.0
             plotPoints%values(p,2)=0!0+(-50.0+j)/50.0
             plotPoints%values(p,3)=0+(-500.0+k)/250.0
          end do
       end do
    end do
    gridSize=gridSize**2

    !Loading molecular orbitals to memory
    if(allocated(orbitalsInGrid)) deallocate(orbitalsInGrid)
    allocate(orbitalsInGrid(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))

    do sysI=1, this%numberOfDisplacedSystems       
       write(auxString, '(I16)') sysI
       molFileI="LOWDIN-"//trim(adjustl(auxString))
       call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileI)

       do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
          call Matrix_constructor(orbitalsInGrid(sysI,speciesID),int(MolecularSystem_getOcupationNumber(speciesID),8),&
               int(gridSize,8), 0.0_8)

          !Get atomic orbital values
          k=0
          do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
             do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)

                call Matrix_constructor( auxMatrix, int(gridSize,8), &
                     int(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital,8), 0.0_8) !orbital

                call ContractedGaussian_getValuesAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), &
                     plotPoints, gridSize, auxMatrix)

                ! call Matrix_show(auxMatrix)

                do j = 1, MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital
                   k=k+1

                   do a=1,MolecularSystem_getOcupationNumber(speciesID)
                      do p = 1 , gridSize
                         ! print *, speciesID, k, a, coefficients(sysI,speciesID)%values(k,a), plotPoints%values(p,1:3), auxMatrix%values(p,j)
                         orbitalsInGrid(sysI,speciesID)%values(a,p)=orbitalsInGrid(sysI,speciesID)%values(a,p)+&
                              this%HFCoefficients(sysI,speciesID)%values(k,a)*auxMatrix%values(p,j)

                      end do
                   end do
                end do

             end do
          end do
       end do
    end do


    !Computing NOCI densities
    if(allocated(densityInGrid)) deallocate(densityInGrid)
    allocate(densityInGrid(molecularSystem_instance%numberOfQuantumSpecies))

    do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies

       if(speciesID .ne. 3 ) cycle

       do state = 1, CONTROL_instance%CI_STATES_TO_PRINT

          call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)

          !Diagonal contributions
          do sysI=1, this%numberOfDisplacedSystems       
             ! call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)
             do a=1,MolecularSystem_getOcupationNumber(speciesID)
                do b=1,MolecularSystem_getOcupationNumber(speciesID)
                   do p = 1 , gridSize
                      densityInGrid(speciesID)%values(p)=densityInGrid(speciesID)%values(p)+&
                           this%configurationCoefficients%values(sysI,state)**2.0*&
                           orbitalsInGrid(sysI,speciesID)%values(a,p)*&
                           orbitalsInGrid(sysI,speciesID)%values(b,p)
                   end do
                end do
             end do
          end do

          !Off Diagonal contributions
          do sysI=1, this%numberOfDisplacedSystems       
             do sysII=sysI+1, this%numberOfDisplacedSystems !sysI+1
                ! call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)
                do a=1,MolecularSystem_getOcupationNumber(speciesID)
                   do b=1,MolecularSystem_getOcupationNumber(speciesID)
                      do p = 1 , gridSize
                         densityInGrid(speciesID)%values(p)=densityInGrid(speciesID)%values(p)+&
                              2.0*this%configurationCoefficients%values(sysI,state)*&
                              this%configurationCoefficients%values(sysII,state)*&
                              this%inverseOverlapMatrices(sysI,sysII,speciesID)%values(b,a)*&
                              this%configurationOverlapMatrix%values(sysI,sysII)*&
                              orbitalsInGrid(sysI,speciesID)%values(a,p)*&
                              orbitalsInGrid(sysII,speciesID)%values(b,p)
                      end do
                   end do
                end do
             end do
          end do

          !density check
          ! densityIntegral=0.0
          ! do p = 1 , gridSize
          ! densityIntegral=densityIntegral+densityInGrid(speciesID)%values(p)*0.02**3.0   
          ! end do
          ! print *, "densityIntegral", densityIntegral

          !1D plot
          outUnit=700
          write(auxString, '(I16)') state
          outFile=trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//"-"//trim(adjustl(auxString))//".2D.dens"

          print *, "printing...", outFile

          open(unit=outUnit, file=trim(outFile), status="replace", form="formatted")          
          write (outUnit,*) "# Z dens for speciesID", speciesID
          do p = 1 , gridSize
             if(plotPoints%values(p,1) .eq. 0.0 .and. plotPoints%values(p,2) .eq. 0.0) &
                  write (outUnit,"(T10,F20.8,E20.8)") plotPoints%values(p,3), densityInGrid(speciesID)%values(p)
          end do
          close(outUnit)

          !2D plot
          outFile=trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//"-"//trim(adjustl(auxString))//".3D.dens"
          print *, "printing...", outFile

          open(unit=outUnit, file=trim(outFile), status="replace", form="formatted")          
          write (outUnit,*) "# X Z dens for speciesID", speciesID

          auxValue=plotPoints%values(1,1)
          do p = 1 , gridSize
             if(plotPoints%values(p,2) .eq. 0.0) then
                if(plotPoints%values(p,1) .gt. auxValue ) then             
                   write (outUnit,"(T10)")
                   auxValue=plotPoints%values(p,1)
                end if
                write (outUnit,"(T10,F20.8,F20.8,E20.8)") plotPoints%values(p,1), plotPoints%values(p,3), densityInGrid(speciesID)%values(p)
             end if
          end do

       end do
    end do

    !$  timeB = omp_get_wtime()
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI density plots : ", timeB - timeA ," (s)"

  end subroutine NonOrthogonalCI_plotDensities

  
end module NonOrthogonalCI_

