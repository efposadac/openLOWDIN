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
  use Math_
  use Exception_
  use MolecularSystem_
  use ParticleManager_
  use Matrix_
  use ReadTransformedIntegrals_
  use Lebedev_
  use Matrix_
  use Vector_
  use Solver_
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

  type, public :: NonOrthogonalCI
     logical :: isInstanced
     integer :: numberOfDisplacedSystems
     integer :: numberOfRejectedSystems
     integer :: numberOfTransformedCenters
     integer :: numberOfIndividualTransformations
     integer :: numberOfUniqueSystems !sort of symmetry
     integer :: numberOfUniquePairs !sort of symmetry
     type(IVector) :: systemTypes  !sort of symmetry
     integer, allocatable :: rotationCenterList(:,:)
     type(Matrix) :: configurationOverlapMatrix, configurationHamiltonianMatrix, configurationCoefficients
     type(IMatrix) :: configurationPairTypes !, uniqueOverlapElements, uniqueHamiltonianElements
     type(Vector) :: statesEigenvalues
     type(Matrix), allocatable :: HFCoefficients(:,:)
     type(MolecularSystem), allocatable :: MolecularSystems(:)
     type(MolecularSystem), allocatable :: uniqueMolecularSystems(:)
     ! type(Matrix), allocatable :: HCoreMatrices(:,:,:)
     ! type(Matrix), allocatable :: inverseOverlapMatrices(:,:,:)
     character(50) :: transformationType
     character(15),allocatable :: systemLabels(:)
     real(8) :: refEnergy
  end type NonOrthogonalCI

  type(NonOrthogonalCI), public :: NonOrthogonalCI_instance

  public :: &
       NonOrthogonalCI_constructor,&
       NonOrthogonalCI_displaceGeometries,&
       NonOrthogonalCI_buildOverlapAndHamiltonianMatrix,&
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
    this%numberOfUniqueSystems=0
    this%numberOfUniquePairs=0

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

    ! print *, "this%rotationCenterList" 
    ! do p=1, size(MolecularSystem_instance%allParticles)
    !    print *, "Particle ", trim(ParticleManager_getSymbol(p)),this%rotationCenterList(p,1), this%rotationCenterList(p,2)
    ! end do
    
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

    !The size here is an upper bound
    allocate(this%MolecularSystems(this%numberOfIndividualTransformations**this%numberOfTransformedCenters), &
         this%systemLabels(this%numberOfIndividualTransformations**this%numberOfTransformedCenters))

    call Vector_constructorInteger(this%systemTypes,this%numberOfIndividualTransformations**this%numberOfTransformedCenters,0)
    
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
    real(8) :: displacement
    real(8) :: testEnergy
    character(50) :: wfnFile, molFile
    character(50) :: arguments(2)
    integer, allocatable :: transformationCounter(:)
    integer :: wfnUnit
    integer :: i,j
    integer :: sysI, speciesID
    integer :: closestSystem
    integer :: systemType
    logical :: newSystemFlag
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
       
       write (*,"(A,A)") "Transformation counter: ", MolecularSystem_instance%description

       call MolecularSystem_showCartesianMatrix()

       !Check if the new system is not to close to previous calculated systems - duplicate protection
       call NonOrthogonalCI_checkNewSystemDisplacement(this,closestSystem,displacement) 

       !Classify the system according to its distance matrix (symmetry) 
       if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
            call NonOrthogonalCI_classifyNewSystem(this,systemType, newSystemFlag) 
       
       if(displacement .lt. CONTROL_instance%CONFIGURATION_DISPLACEMENT_THRESHOLD ) then
          print *, " Skipping system with distance ", displacement , "a.u. from system ", closestSystem
          this%numberOfRejectedSystems=this%numberOfRejectedSystems+1                      
       else
!!!Run a HF calculation at each displaced geometry 
          call MolecularSystem_saveToFile()          
          call NonOrthogonalCI_runHF(testEnergy)

          !!Screen geometries with high energies
          if( CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD .ne. 0.0 .and. &
               testEnergy .gt. this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD) then
             print *, " Skipping system with high energy", testEnergy
             this%numberOfRejectedSystems=this%numberOfRejectedSystems+1                      
          else
             this%numberOfDisplacedSystems=this%numberOfDisplacedSystems+1
             if(newSystemFlag) then
                this%numberOfUniqueSystems=this%numberOfUniqueSystems+1
                this%systemTypes%values(this%numberOfDisplacedSystems)=this%numberOfUniqueSystems
             else
                this%systemTypes%values(this%numberOfDisplacedSystems)=systemType
             end if

             if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
                write (*,"(A,I5,A,I10,A,F20.12)") "Saving system of type ",  this%systemTypes%values(this%numberOfDisplacedSystems) , &
                     " with ID ", this%numberOfDisplacedSystems, " and energy", testEnergy
             else
                write (*,"(A,I10,A,F20.12)") "Saving system with ID ", this%numberOfDisplacedSystems, " and energy", testEnergy
             end if
             !!Copy the molecular system to the NonOrthogonalCI object
             call NonOrthogonalCI_saveSystem(this,this%numberOfDisplacedSystems)
          end if
       end if
       
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
    if(this%numberOfRejectedSystems .gt. 0) &
         write (*,'(A10,I10,A,F18.12,A,F18.12)') "Rejected ", this%numberOfRejectedSystems, &
         " geometries with energy higher than", this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD, &     
         " or with distance to other systems lower than", CONTROL_instance%CONFIGURATION_DISPLACEMENT_THRESHOLD
    print *, ""
    
    allocate(this%HFCoefficients(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    call Matrix_constructor(this%configurationHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Matrix_constructor(this%configurationOverlapMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 1.0_8)
    call Matrix_constructorInteger(this%configurationPairTypes,  int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0)
    
    !!Read HF energies and coefficients and fill CI matrix diagonals 
    ! minEnergy=0.0    
    do sysI=1, this%numberOfDisplacedSystems
       write(this%systemLabels(sysI), '(A)') &
            trim(this%MolecularSystems(sysI)%description)

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
       
    end do
    
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
    integer :: i,j,k,p,q,mu

    MolecularSystem_instance%description=""
    do i=1,this%numberOfTransformedCenters
       write(MolecularSystem_instance%description, '(A,I3)') trim(MolecularSystem_instance%description), transformationCounter(i)
    end do

    particleManager_instance => MolecularSystem_instance%allParticles
    
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
                            ! call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                            MolecularSystem_instance%allParticles(p)%particlePtr%origin=displacedOrigin
                            do mu = 1, MolecularSystem_instance%allParticles(p)%particlePtr%basis%length
                               MolecularSystem_instance%allParticles(p)%particlePtr%basis%contraction(mu)%origin = displacedOrigin
                            end do
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
       
          do i=1,CONTROL_instance%ROTATIONAL_SCAN_GRID
             do j=1,CONTROL_instance%NESTED_ROTATIONAL_GRIDS
                displacementId=displacementId+1
                if(displacementId .eq. transformationCounter(center) ) then
                   do p=1, size(MolecularSystem_instance%allParticles)
                      if(this%rotationCenterList(p,1).eq. center ) then

                         do q=1, size(originalMolecularSystem%allParticles)
                            if(this%rotationCenterList(q,1) .eq. center ) then
                               centerX=originalMolecularSystem%allParticles(this%rotationCenterList(q,2))%particlePtr%origin(1)
                               centerY=originalMolecularSystem%allParticles(this%rotationCenterList(q,2))%particlePtr%origin(2)
                               centerZ=originalMolecularSystem%allParticles(this%rotationCenterList(q,2))%particlePtr%origin(3)
                            end if
                         end do

                         distanceToCenter=sqrt((originalMolecularSystem%allParticles(p)%particlePtr%origin(1)-centerX)**2 &
                              +(originalMolecularSystem%allParticles(p)%particlePtr%origin(2)-centerY)**2 &
                              +(originalMolecularSystem%allParticles(p)%particlePtr%origin(3)-centerZ)**2)

                         distanceToCenter=distanceToCenter+&
                              CONTROL_instance%NESTED_GRIDS_DISPLACEMENT*(j-(CONTROL_instance%NESTED_ROTATIONAL_GRIDS+1)/2.0)
                                                  
                         displacedOrigin(1)=centerX+X(i)*distanceToCenter
                         displacedOrigin(2)=centerY+Y(i)*distanceToCenter
                         displacedOrigin(3)=centerZ+Z(i)*distanceToCenter

                         ! call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                         MolecularSystem_instance%allParticles(p)%particlePtr%origin=displacedOrigin
                         do mu = 1, MolecularSystem_instance%allParticles(p)%particlePtr%basis%length
                            MolecularSystem_instance%allParticles(p)%particlePtr%basis%contraction(mu)%origin = displacedOrigin
                         end do
                      end if
                   end do
                end if
             end do
          end do
       end do
    end if
                   
  end subroutine NonOrthogonalCI_transformCoordinates

  !>
  !! @brief Computes the distance between the particles of latest generated molecular system with all the previous saved ones 
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NonOrthogonalCI_checkNewSystemDisplacement(this,closestSystem,displacement) 
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: closestSystem
    real(8) :: displacement

    integer :: sysI, i
    type(Vector), allocatable :: displacementVector(:)
    real(8) :: dispSum
    
    displacement=1.0E8
    
    allocate(displacementVector(molecularSystem_instance%numberOfQuantumSpecies))
        
    do sysI=1, this%numberOfDisplacedSystems

       call MolecularSystem_GetTwoSystemsDisplacement(this%MolecularSystems(sysI), molecularSystem_instance, displacementVector(:))

       dispSum=0.0
       do i=1, molecularSystem_instance%numberOfQuantumSpecies
          dispSum=dispSum+sum(displacementVector(i)%values(:))
       end do
       if(dispSum/molecularSystem_instance%numberOfQuantumSpecies &
            .lt. CONTROL_instance%CONFIGURATION_DISPLACEMENT_THRESHOLD ) then
          displacement=dispSum/molecularSystem_instance%numberOfQuantumSpecies
          closestSystem=sysI
          exit
       end if
    end do

  end subroutine NonOrthogonalCI_checkNewSystemDisplacement

  !>
  !! @brief Classify the new system by comparing its distance matrix to previosly saved systems
  !!
  !! @param this, systemType: integer defining system equivalence type,  newSystemFlag: returns if the system is new or not
  !<
  subroutine NonOrthogonalCI_classifyNewSystem(this, systemType, newSystemFlag) 
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: systemType
    logical :: newSystemFlag
    
    type(MolecularSystem) :: currentMolecularSystem
    type(Matrix) :: currentDistanceMatrix,previousDistanceMatrix

    integer :: sysI, i, checkingType
    logical :: match
    
    call MolecularSystem_copyConstructor(currentMolecularSystem, molecularSystem_instance)
    systemType=0
    newSystemFlag=.true.
    currentDistanceMatrix=ParticleManager_getDistanceMatrix()

    ! print *, "Current distance matrix"
    ! call Matrix_show(currentDistanceMatrix)

    types: do checkingType=1, this%numberOfUniqueSystems
       ! print *, "checkingType", checkingType
       systems: do sysI=1, this%numberOfDisplacedSystems

          if(this%systemTypes%values(sysI) .eq. checkingType) then
             call MolecularSystem_copyConstructor(molecularSystem_instance, this%MolecularSystems(sysI))

             previousDistanceMatrix=ParticleManager_getDistanceMatrix()

             ! print *, "Comparing with previous distance matrix", checkingType
             ! call Matrix_show(previousDistanceMatrix)          
          
             match=.true.
             do i=1, size(currentDistanceMatrix%values(:,1))
                if(sum(abs(currentDistanceMatrix%values(i,:) - previousDistanceMatrix%values(i,:))) .gt. &
                     CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) then
                   match=.false.
                   exit
                end if
             end do
             
             ! print *, "match?", match

             if(match) then
                systemType=this%systemTypes%values(sysI)
                newSystemFlag=.false.
                exit types
             else
                cycle types
             end if
          end if
       end do systems
    end do types
    
    ! print *, "newSystemFlag", newSystemFlag
    
    call MolecularSystem_copyConstructor(molecularSystem_instance, currentMolecularSystem)

  end subroutine NonOrthogonalCI_classifyNewSystem

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

  ! >
  ! @brief Saves molecular system and wfn files for a displaced system 
  
  ! @param systemID
  ! <
  subroutine NonOrthogonalCI_saveSystem(this,systemID)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: systemID
    character(50) :: wfnFile, molFile, auxString
    integer :: wfnUnit
    
    !!Copy the molecular system to the NonOrthogonalCI object
    call MolecularSystem_copyConstructor(this%MolecularSystems(systemID), molecularSystem_instance)    
    
    !!Save the molecular system and WF results to separate files for each geometry     

    write(auxString, '(I16)') systemID
    wfnFile="lowdin-"//trim(adjustl(auxString))//".wfn"
    call system("cp lowdin.wfn "//wfnFile)
    ! molFile="LOWDIN-"//trim(adjustl(auxString))
    ! call MolecularSystem_saveToFile(molFile)
    ! print *, "molFile", molFile, "wfnFile", wfnFile
    
  end subroutine NonOrthogonalCI_saveSystem

  !>
  !! @brief Computes overlap and hamiltonian non orthogonal CI matrices for previously calculated molecular systems at different geometries
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_buildOverlapAndHamiltonianMatrix(this)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: sysI,sysII, preSysI, preSysII
    type(Matrix), allocatable :: mergedCoefficients(:), inverseOverlapMatrices(:)
    type(IVector), allocatable :: sysIbasisList(:),sysIIbasisList(:)
    real(8) :: overlapUpperBound
    integer :: symmetryEquivalentElements, prescreenedElements, overlapScreenedElements
    logical :: newPairFlag

    integer :: matrixUnit
    character(50) :: matrixFile
    
    real(8) :: timePrescreen, timeSymmetry, timeOverlap, timeTwoIntegrals
    real(8) :: timeA
    real(8) :: timeB

    timePrescreen=0.0
    timeSymmetry=0.0
    timeOverlap=0.0
    timeTwoIntegrals=0.0
        
    print *, ""
    print *, "A prescreening of the overlap matrix elements is performed for the heavy species"
    write (*,'(A,ES8.1)') "Overlap and Hamiltonian matrix elements are saved for pairs with overlap higher than",&
         CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD
    print *, "For pairs with lower overlap, setting H(I,II)=0, S(I,II)=0"
    print *, ""

    allocate(mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         inverseOverlapMatrices(molecularSystem_instance%numberOfQuantumSpecies))
    allocate(sysIbasisList(molecularSystem_instance%numberOfQuantumSpecies),&
         sysIIbasisList(molecularSystem_instance%numberOfQuantumSpecies))

    prescreenedElements=0
    symmetryEquivalentElements=0
    overlapScreenedElements=0

    matrixUnit=290
    matrixFile= trim(CONTROL_instance%INPUT_FILE)//"NOCI-Matrix.ci"

    print *, "computing NOCI overlap and hamiltonian matrices... saving them to ", trim(matrixFile)

    open(unit=matrixUnit, file=trim(matrixFile), status="replace", form="formatted")

    write (matrixUnit,'(A20,I20)') "MatrixSize", this%numberOfDisplacedSystems

    if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
       write (matrixUnit,'(A10,A10,A10,A20,A20)') "Conf. ", "Conf. ", "Type ", "Overlap ","Hamiltonian "
    else
       write (matrixUnit,'(A10,A10,A20,A20)') "Conf. ", "Conf. ", "Overlap ","Hamiltonian "
    end if
    
    systemI: do sysI=1, this%numberOfDisplacedSystems       
       !Save diagonal elements
       if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
          write (matrixUnit,'(I10,I10,I10,ES20.12,ES20.12)') sysI, sysI, this%configurationPairTypes%values(sysI,sysI), &
               this%configurationOverlapMatrix%values(sysI,sysI), this%configurationHamiltonianMatrix%values(sysI,sysI)               
       else
          write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, sysI, &
               this%configurationOverlapMatrix%values(sysI,sysI), this%configurationHamiltonianMatrix%values(sysI,sysI)               
       end if
       
       systemII: do sysII=sysI+1, this%numberOfDisplacedSystems !sysI+1
          !initialize overlap
          this%configurationOverlapMatrix%values(sysI,sysII)=1.0
          
          !This generates a new molecular system
          ! print *, "Merging systems from geometries ", sysI, sysII
          call MolecularSystem_mergeTwoSystems(molecularSystem_instance, this%MolecularSystems(sysI), this%MolecularSystems(sysII), &
               sysIbasisList,sysIIbasisList)

          ! call MolecularSystem_showInformation()  
          ! call MolecularSystem_showParticlesInformation()
          ! call MolecularSystem_showCartesianMatrix()
                   
          !$  timeA = omp_get_wtime()
          !Estimates overlap with a 1s-1s integral approximation
          call NonOrthogonalCI_prescreenOverlap(this,sysI,sysII,overlapUpperBound)
          !$  timeB = omp_get_wtime()
          !$  timePrescreen=timePrescreen+(timeB - timeA)

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               overlapUpperBound .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             ! print *, "preskipping elements", sysI, sysII, "with overlap estimated as", overlapUpperBound
             this%configurationOverlapMatrix%values(sysI,sysII)=0.0
             this%configurationOverlapMatrix%values(sysII,sysI)=0.0
             this%configurationHamiltonianMatrix%values(sysI,sysII)=0.0
             this%configurationHamiltonianMatrix%values(sysII,sysI)=0.0
             prescreenedElements=prescreenedElements+1
             cycle systemII
          end if

          !$  timeA = omp_get_wtime()
          !!Check symmetry of the element
          if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
             call NonOrthogonalCI_classifyConfigurationPair(this,sysI,sysII,newPairFlag)
             !$  timeB = omp_get_wtime()
             !$  timeSymmetry=timeSymmetry+(timeB - timeA)

             !!Copy results from previously computed equivalent elements
             if (newPairFlag .eqv. .false.) then
                do preSysI=1, sysI
                   do preSysII=preSysI+1, sysII                   
                      if(this%configurationPairTypes%values(preSysI,preSysII) .eq. this%configurationPairTypes%values(sysI,sysII)) then
                         this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(preSysI,preSysII)
                         this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)
                         this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(preSysI,preSysII)
                         this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)
                         symmetryEquivalentElements=symmetryEquivalentElements+1

                         if( this%configurationOverlapMatrix%values(sysI,sysII) .ne. 0.0) &
                              write (*,'(A,I10,I10,A,I10,A,ES20.12,ES20.12)') "Pair ",sysI, sysII," is type ", &
                              this%configurationPairTypes%values(sysI,sysII), " Overlap and Hamiltonian elements", &
                              this%configurationOverlapMatrix%values(sysI,sysII), this%configurationHamiltonianMatrix%values(sysI,sysII)

                         cycle systemII
                      end if
                   end do
                end do
             end if
          end if
         
          !! Save to file to compute integrals
          if ( trim(CONTROL_instance%INTEGRAL_STORAGE) .ne. "DIRECT" ) call MolecularSystem_saveToFile()

          !! Merge occupied coefficients into a single matrix 
          call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(sysI,:),this%HFCoefficients(sysII,:),&
               this%MolecularSystems(sysI),this%MolecularSystems(sysII),&
               sysIbasisList,sysIIbasisList,mergedCoefficients)
          !$  timeA = omp_get_wtime()

          call NonOrthogonalCI_computeOverlapAndHCoreElements(this,sysI,sysII,mergedCoefficients,sysIbasisList,sysIIbasisList,inverseOverlapMatrices)
          !$  timeB = omp_get_wtime()
          !$  timeOverlap=timeOverlap+(timeB - timeA)
          
          !! SKIP ENERGY EVALUATION IF OVERLAP IS TOO LOW

          if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
               abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
             ! print *, "screening elements", sysI, sysII, "with overlap", this%configurationOverlapMatrix%values(sysI,sysII)
             this%configurationOverlapMatrix%values(sysI,sysII)=0.0
             this%configurationOverlapMatrix%values(sysII,sysI)=0.0
             this%configurationHamiltonianMatrix%values(sysI,sysII)=0.0
             this%configurationHamiltonianMatrix%values(sysII,sysI)=0.0
             ! this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*&
             !      (this%configurationHamiltonianMatrix%values(sysI,sysI)+this%configurationHamiltonianMatrix%values(sysII,sysII))/2.0
             ! this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)             
             overlapScreenedElements=overlapScreenedElements+1
             cycle systemII
          end if

          !$  timeA = omp_get_wtime()
          call NonOrthogonalCI_twoParticlesContributions(this,sysI,sysII,inverseOverlapMatrices,mergedCoefficients)
          !$  timeB = omp_get_wtime()
          !$  timeTwoIntegrals=timeTwoIntegrals+(timeB - timeA)

          !!This is a symmetry test, assume positive phase
          ! if( this%configurationOverlapMatrix%values(sysI,sysII) .lt. 0.0) then
          !    this%configurationOverlapMatrix%values(sysI,sysII)=-this%configurationOverlapMatrix%values(sysI,sysII)
          !    this%configurationOverlapMatrix%values(sysII,sysI)=-this%configurationOverlapMatrix%values(sysII,sysI)
          !    this%configurationHamiltonianMatrix%values(sysI,sysII)=-this%configurationHamiltonianMatrix%values(sysI,sysII)
          !    this%configurationHamiltonianMatrix%values(sysII,sysI)=-this%configurationHamiltonianMatrix%values(sysII,sysI)
          ! end if
          
          if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
             write (matrixUnit,'(I10,I10,I10,ES20.12,ES20.12)') sysI, sysII, this%configurationPairTypes%values(sysI,sysII), &
                  this%configurationOverlapMatrix%values(sysI,sysII), this%configurationHamiltonianMatrix%values(sysI,sysII)               
          else
             write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, sysII, &
                  this%configurationOverlapMatrix%values(sysI,sysII), this%configurationHamiltonianMatrix%values(sysI,sysII)               
          end if
       end do systemII
    end do systemI

    close(matrixUnit)

    print *, ""
    print *, "Configuration pairs skipped by overlap prescreening: ", prescreenedElements
    if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
         print *, "Configuration pairs skipped by symmetry equivalence: ", symmetryEquivalentElements
    print *, "Configuration pairs skipped by overlap    screening: ", overlapScreenedElements
    print *, "Overlap integrals computed for    ", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
         -prescreenedElements-symmetryEquivalentElements, "configuration pairs"
    print *, "Four center integrals computed for", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
         -prescreenedElements-symmetryEquivalentElements-overlapScreenedElements, "configuration pairs"
    print *, ""

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for overlap prescreening : ", timePrescreen ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for element symmetry     : ", timeSymmetry ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for two index integrals  : ", timeOverlap ," (s)"
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for four index integrals : ", timeTwoIntegrals ," (s)"
    print *, ""
    
  end subroutine NonOrthogonalCI_buildOverlapAndHamiltonianMatrix

  !>
  !! @brief Merges the occupied orbitals coefficients from two systems
  !! molecularsystem_instance brings the merged molecular system
  !! @param occupationI and occupationII: Number of orbitals to merge from each matrix. 
  !! sysBasisList: array indicating which basis functions of the merged molecular system belong to sysI and sysII Merged Coefficients: Matrices for output.
  !<
  subroutine NonOrthogonalCI_mergeCoefficients(coefficientsI,coefficientsII,molecularSystemI,molecularSystemII,&
       sysIbasisList,sysIIbasisList,mergedCoefficients)
    type(Matrix), intent(in) :: coefficientsI(*), coefficientsII(*)
    type(MolecularSystem), intent(in) :: molecularSystemI, molecularSystemII
    type(IVector), intent(in) :: sysIbasisList(*), sysIIbasisList(*)
    type(Matrix), intent(out) :: mergedCoefficients(*)
    
    character(50) :: wfnFile
    character(50) :: arguments(2)
    integer :: wfnUnit
    integer :: speciesID, i, j, k, l, m, mu, nu, notCommonBasis
    type(Matrix) :: auxMatrix
    type(Vector) :: auxVector
    
    !! Mix coefficients of occupied orbitals of both systems    
    !!Create a dummy density matrix to lowdin.wfn file
    wfnUnit = 500
    wfnFile = "lowdin.wfn"
    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
       
       arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

       arguments(1) = "COEFFICIENTS"

       !    !Max: to make the matrix square for the integral calculations for configuration pairs, and rectangular for the merged coefficients of all systems 
       call Matrix_constructor(mergedCoefficients(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
            int(max(MolecularSystem_getTotalNumberOfContractions(speciesID),MolecularSystem_getOcupationNumber(speciesID)),8), 0.0_8 )

       ! print *, "sysI coefficients for ", speciesID
       ! call Matrix_show(coefficientsI(speciesID))
       ! print *, "sysII coefficients for ", speciesID
       ! call Matrix_show(coefficientsII(speciesID))

       !sysI orbitals on the left columns, sysII on the right columns

       !sysI coefficients
       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
          if((sysIbasisList(speciesID)%values(mu) .ne. 0) ) then
             do i=1, MolecularSystem_getOcupationNumber(speciesID,molecularSystemI)!sysI
                mergedCoefficients(speciesID)%values(mu,i)=coefficientsI(speciesID)%values(sysIbasisList(speciesID)%values(mu),i)
                ! print *, "sys I", mu, i, mergedCoefficients(speciesID)%values(mu,i)
             end do
          end if
       end do

       ! !sysII coefficients
       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
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

       call Matrix_writeToFile(mergedCoefficients(speciesID), unit=wfnUnit, binary=.true., arguments = arguments )
       
       arguments(1) = "DENSITY"
       call Matrix_constructor(auxMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
            int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8 )

       auxMatrix%values=1.0
       
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

       ! Only occupied orbitals are going to be transformed - handled in integral transformation program
       ! print *, "removed", MolecularSystem_getTotalNumberOfContractions(speciesID)-MolecularSystem_getOcupationNumber(speciesID)
       arguments(1) = "REMOVED-ORBITALS"
       call Vector_writeToFile(unit=wfnUnit, binary=.true., &
            value=real(MolecularSystem_getTotalNumberOfContractions(speciesID)-MolecularSystem_getOcupationNumber(speciesID),8),&
            arguments= arguments )

    end do
    close(wfnUnit)
    
  end subroutine NonOrthogonalCI_mergeCoefficients


  !>
  !! @brief Computes an upper bound of the overlap between two configurations, based on the max distance between particles of the same species and the lowest exponent of the basis set functions. Assumes a localized hartree product for the heaviest species
  !!
  !! @param sysI and sysII: molecular system indices. estimatedOverlap: output value
  !<
  subroutine NonOrthogonalCI_prescreenOverlap(this,sysI,sysII,estimatedOverlap)
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII !Indices of the systems to screen
    real(8) :: estimatedOverlap

    type(Vector), allocatable :: displacementVector(:)    
    integer :: speciesID, k, l, m
    real(8) :: massThreshold, minExponent, speciesOverlap
    
    !displacement vectors contains the max distance between equivalent basis function centers 
    allocate(displacementVector(this%MolecularSystems(sysI)%numberOfQuantumSpecies))
  
    call MolecularSystem_GetTwoSystemsDisplacement(this%MolecularSystems(sysI), this%MolecularSystems(sysII),displacementVector(:))

    estimatedOverlap=1.0

    !only compute for heavy particles, maybe should be a control parameter
    massThreshold=10.0   
    
    do speciesID = 1, this%MolecularSystems(sysI)%numberOfQuantumSpecies
       if(this%MolecularSystems(sysI)%species(speciesID)%mass .lt. massThreshold) cycle
       speciesOverlap=1.0
       !!get smallest exponent of the basis set
       do k = 1, size(this%MolecularSystems(sysI)%species(speciesID)%particles)
          minExponent=1.0E8
          do l = 1, size(this%MolecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction)             
             do m = 1, size(this%MolecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents)
                if(this%MolecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents(m).lt.minExponent) &
                     minExponent=this%MolecularSystems(sysI)%species(speciesID)%particles(k)%basis%contraction(l)%orbitalExponents(m)
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

  end subroutine NonOrthogonalCI_prescreenOverlap

  !>
  !! @brief Classify the sysI and sysII pair according to their distance matrix
  !!
  !! @param sysI and sysII: molecular system indices.
  !<
  subroutine NonOrthogonalCI_classifyConfigurationPair(this,currentSysI,currentSysII,newPairFlag)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: currentSysI, currentSysII !Indices of the systems to classify
    logical :: newPairFlag
    
    type(MolecularSystem) :: currentMolecularSystem
    type(Matrix) :: currentDistanceMatrix,previousDistanceMatrix
    
    integer :: sysI, sysII, i, checkingType
    logical :: match
    
    call MolecularSystem_copyConstructor(currentMolecularSystem, molecularSystem_instance)
    newPairFlag=.true.
    currentDistanceMatrix=ParticleManager_getDistanceMatrix()

    ! print *, "Current distance matrix"
    ! call Matrix_show(currentDistanceMatrix)

    types: do checkingType=1, this%numberOfUniquePairs
       ! print *, "checkingType", checkingType
       systemI: do sysI=1, currentSysI
          systemII: do sysII=sysI+1, currentSysII

             if(sysI .eq. currentSysI .and. sysII .eq. currentSysII ) cycle types

             if((this%configurationPairTypes%values(sysI,sysII) .eq. checkingType) .and. &
                (this%systemTypes%values(sysI) .eq. this%systemTypes%values(currentSysI)) .and. & 
                (this%systemTypes%values(sysII) .eq. this%systemTypes%values(currentSysII))) then

                ! call MolecularSystem_mergeTwoSystems(molecularSystem_instance, this%MolecularSystems(sysI), this%MolecularSystems(sysII))
                
                previousDistanceMatrix=ParticleManager_getDistanceMatrix()

                ! print *, "Comparing with previous distance matrix", checkingType
                ! call Matrix_show(previousDistanceMatrix)          
          
                match=.true.
                do i=1, size(currentDistanceMatrix%values(:,1))
                   if(sum(abs(currentDistanceMatrix%values(i,:) - previousDistanceMatrix%values(i,:))) .gt. &
                        CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) then
                      match=.false.
                      exit
                   end if
                end do
             
                if(match) then
                   newPairFlag=.false.
                   this%configurationPairTypes%values(currentSysI,currentSysII)=this%configurationPairTypes%values(sysI,sysII)
                   exit types
                else
                   cycle types
                end if
             end if
          end do systemII
       end do systemI
    end do types

    if(newPairFlag) then
       this%numberOfUniquePairs=this%numberOfUniquePairs+1
       this%configurationPairTypes%values(currentSysI,currentSysII)=this%numberOfUniquePairs
    end if

    if(this%configurationPairTypes%values(currentSysI,currentSysII).eq.0) then
       print *, "newPairFlag", newPairFlag
       print *, currentSysI, currentSysII, this%configurationPairTypes%values(currentSysI,currentSysII)
       STOP "I found a type zero"
    end if
    call MolecularSystem_copyConstructor(molecularSystem_instance, currentMolecularSystem)

  end subroutine NonOrthogonalCI_classifyConfigurationPair


  
  !>
  !! @brief Computes overlap matrix element between two configurations along with one particle energy contributions
  !!
  !! @param sysI and sysII: molecular system indices. Merged Coefficients: Mixed molecular system coefficients. Sys basis list indicate the basis functions of each sysI and sysII in the merged molecular system. inverseOverlapMatrices: output required for two particle contributions
  !<
  subroutine NonOrthogonalCI_computeOverlapAndHCoreElements(this,sysI,sysII,mergedCoefficients, &
       sysIbasisList, sysIIbasisList,inverseOverlapMatrices)

    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII 
    type(Matrix) :: mergedCoefficients(*), inverseOverlapMatrices(*) 
    type(IVector) :: sysIbasisList(*), sysIIbasisList(*)

    integer :: integralsUnit
    character(50) :: integralsFile    
    character(50) :: arguments(2)
    integer :: speciesID
    integer :: a,b,bb,mu,nu    
    type(Matrix) :: auxMatrix
    type(Matrix) :: molecularOverlapMatrix
    type(Matrix), allocatable :: auxOverlapMatrix(:), auxKineticMatrix(:), auxAttractionMatrix(:), molecularHCoreMatrix(:)
    type(Vector) :: overlapDeterminant
    real(8) :: oneParticleEnergy
    
    allocate(auxOverlapMatrix(molecularSystem_instance%numberOfQuantumSpecies), &
         auxKineticMatrix(molecularSystem_instance%numberOfQuantumSpecies), &
         auxAttractionMatrix(molecularSystem_instance%numberOfQuantumSpecies), &
         molecularHCoreMatrix(molecularSystem_instance%numberOfQuantumSpecies))
    
    !! Calculate one- particle integrals  
!!!!Overlap first
    if ( trim(CONTROL_instance%INTEGRAL_STORAGE) .eq. "DIRECT" ) then
       call DirectIntegralManager_getOverlapIntegrals(auxOverlapMatrix)
    else
       integralsUnit = 30
       integralsFile = "lowdin.opints"

       call system("lowdin-ints.x ONE_PARTICLE >> extendedOutput.eut")

       open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")
       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

          arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
          arguments(1) = "OVERLAP"

          auxOverlapMatrix(speciesID) = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               unit=integralsUnit, binary=.true., arguments=arguments)
       end do
       close(integralsUnit)
    end if
    
    call Vector_constructor(overlapDeterminant, molecularSystem_instance%numberOfQuantumSpecies, 0.0_8)        
    
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       !!Test 

       ! print *, "auxOverlapMatrix", speciesID
       ! call Matrix_show(auxOverlapMatrix)

       call Matrix_constructor(molecularOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
            int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )

       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysI
          if(sysIbasisList(speciesID)%values(mu) .ne. 0 ) then
             do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                if(sysIIbasisList(speciesID)%values(nu) .ne. 0) then
                   do a=1, MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)) !sysI
                      do b=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+1, &
                           MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+ &
                           MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)) !sysII
                         bb=b-MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                         ! print *, "a, b, mu, nu, coefI, coefII", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b)
                         !      auxOverlapMatrix%values(mu,nu)

                         molecularOverlapMatrix%values(a,bb)=molecularOverlapMatrix%values(a,bb)+&
                              mergedCoefficients(speciesID)%values(mu,a)*&
                              mergedCoefficients(speciesID)%values(nu,b)*&
                              auxOverlapMatrix(speciesID)%values(mu,nu)
                      end do
                   end do
                end if
             end do
          end if
       end do

       ! print *, "molecularOverlapMatrix sysI, sysII, speciesID", sysI, sysII, speciesID
       ! call Matrix_show(molecularOverlapMatrix)

       inverseOverlapMatrices(speciesID)=Matrix_inverse(molecularOverlapMatrix)

       ! print *, "inverseOverlapMatrices sysI, sysII", speciesID, sysI, sysII
       ! call Matrix_show(inverseOverlapMatrices(speciesID))

       call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant%values(speciesID),method="LU")             

       ! print *, "OverlapDeterminantLU speciesID, sysI, sysII", speciesID, sysI, sysII, overlapDeterminant%values(speciesID)
       
       this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*overlapDeterminant%values(speciesID)

    end do

    this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)
    
    !!Skip the rest of the evaluation if the overlap is smaller than the threshold
    if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
         abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) return

    !!Compute hcore if overlap is significant
    if ( trim(CONTROL_instance%INTEGRAL_STORAGE) .eq. "DIRECT" ) then
       call DirectIntegralManager_getKineticIntegrals(auxKineticMatrix)
       call DirectIntegralManager_getAttractionIntegrals(auxAttractionMatrix)
    else
       open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")
       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
          arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
          arguments(1) = "KINETIC"
          auxKineticMatrix(speciesID) = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               unit=integralsUnit, binary=.true., arguments=arguments)
          arguments(1) = "ATTRACTION"
          auxAttractionMatrix(speciesID) = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
               unit=integralsUnit, binary=.true., arguments=arguments)
       end do
       close(integralsUnit)
    end if

!!!!HCore
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       !! Incluiding mass effect       
       if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
          auxKineticMatrix(speciesID)%values =  &
            auxKineticMatrix(speciesID)%values * &
            ( 1.0_8/MolecularSystem_getMass( speciesID ) -1.0_8 / ParticleManager_getTotalMass() )
       else
          auxKineticMatrix(speciesID)%values =  &
            auxKineticMatrix(speciesID)%values / &
            MolecularSystem_getMass( speciesID )
       end if

       !! Including charge
       auxAttractionMatrix(speciesID)%values=auxAttractionMatrix(speciesID)%values*(-MolecularSystem_getCharge(speciesID))                         
       
       !!Test 
       call Matrix_constructor(molecularHCoreMatrix(speciesID), int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
            int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )

       ! print *, "auxKineticMatrix", speciesID
       ! call Matrix_show(auxKineticMatrix)
       ! print *, "auxAttractionMatrix", speciesID
       ! call Matrix_show(auxAttractionMatrix)

       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysI
          if(sysIbasisList(speciesID)%values(mu) .ne. 0) then
             do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysII
                if(sysIIbasisList(speciesID)%values(nu) .ne. 0) then
                   do a=1, MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)) !sysI
                      do b=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+1, &
                           MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+&
                           MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)) !sysII
                         bb=b-MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))

                         ! print *, "hcore", a, b, mu, nu, mergedCoefficients(speciesID)%values(mu,a), mergedCoefficients(speciesID)%values(nu,b), &
                         !      auxKineticMatrix%values(mu,nu)/MolecularSystem_getMass(speciesID)+&
                         !      auxAttractionMatrix%values(mu,nu)*(-MolecularSystem_getCharge(speciesID))

                         molecularHCoreMatrix(speciesID)%values(a,bb)=molecularHCoreMatrix(speciesID)%values(a,bb)+&
                              mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                              (auxKineticMatrix(speciesID)%values(mu,nu)+auxAttractionMatrix(speciesID)%values(mu,nu))                         
                      end do
                   end do
                end if
             end do
          end if
       end do
       
       ! print *, "molecularHCoreMatrix", speciesID
       ! call Matrix_show(molecularHCoreMatrix(speciesID))
       !!End test                          
    end do


    !!Point charge-Point charge repulsion
    this%configurationHamiltonianMatrix%values(sysI,sysII)=MolecularSystem_getPointChargesEnergy()
    ! print *, "Point charge-Point charge repulsion", MolecularSystem_getPointChargesEnergy()

    !!One Particle Terms
    do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
       oneParticleEnergy=0.0
       do a=1, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8) !sysI
          do b=1, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8) !sysII
             oneParticleEnergy=oneParticleEnergy+ molecularHCoreMatrix(speciesID)%values(a,b)*&
                  inverseOverlapMatrices(speciesID)%values(b,a)
          end do
       end do
       this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+oneParticleEnergy
       ! print *, "oneParticleEnergy for species", speciesID, oneParticleEnergy
    end do
        
  end subroutine NonOrthogonalCI_computeOverlapAndHCoreElements
  !>
  !! @brief Computes the two particles contributions to the non diagonal elements of the hamiltonian matrix
  !!
  !! @param this, sysI,sysII: system indexes, inverseOverlapMatrices, mergedCoefficients are required to evaluate the elements
  !<
  subroutine NonOrthogonalCI_twoParticlesContributions(this,sysI,sysII,inverseOverlapMatrices,mergedCoefficients)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII 
    type(Matrix) :: inverseOverlapMatrices(*) 
    type(Matrix) :: mergedCoefficients(*) 

    type(matrix), allocatable :: fourCenterIntegrals(:,:)
    type(imatrix), allocatable :: twoIndexArray(:),fourIndexArray(:)
    integer :: ssize1, auxIndex, auxIndex1
    integer :: a,b,bb,c,d,dd,i,j
    real(8) :: interactionEnergy

    allocate(fourCenterIntegrals(MolecularSystem_instance%numberOfQuantumSpecies,MolecularSystem_instance%numberOfQuantumSpecies), &
         twoIndexArray(MolecularSystem_instance%numberOfQuantumSpecies), &
         fourIndexArray(MolecularSystem_instance%numberOfQuantumSpecies))

    !!Fill indexes arrays
    do i=1, MolecularSystem_instance%numberOfQuantumSpecies       
       ! print *, "reading integrals species", i
       !!Two particle integrals indexes
       call Matrix_constructorInteger(twoIndexArray(i),  &
            int(max(MolecularSystem_getTotalNumberOfContractions(i),MolecularSystem_getOcupationNumber(i)),8), &
            int(max(MolecularSystem_getTotalNumberOfContractions(i),MolecularSystem_getOcupationNumber(i)),8) , 0 )

       c = 0
       do a=1,max(MolecularSystem_getTotalNumberOfContractions(i),MolecularSystem_getOcupationNumber(i))
          do b=a, max(MolecularSystem_getTotalNumberOfContractions(i),MolecularSystem_getOcupationNumber(i))
             c = c + 1
             twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
             twoIndexArray(i)%values(b,a) = twoIndexArray(i)%values(a,b)
          end do
       end do

       ssize1 = max(MolecularSystem_getTotalNumberOfContractions( i ),MolecularSystem_getOcupationNumber(i))
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
    if ( trim(CONTROL_instance%INTEGRAL_STORAGE) .eq. "DIRECT" ) then
       call NonOrthogonalCI_transformIntegralsMemory(MolecularSystem_instance, mergedCoefficients, &
            twoIndexArray, fourIndexArray, fourCenterIntegrals)

    else
       call system("lowdin-ints.x TWO_PARTICLE_R12 >> extendedOutput.eut")
       ! call system("lowdin-ints.x TWO_PARTICLE_R12")

       !! Transform integrals
       call system("lowdin-integralsTransformation.x >> extendedOutput.eut")
       ! call system("lowdin-integralsTransformation.x")

       ! Load integrals      
       do i=1, MolecularSystem_instance%numberOfQuantumSpecies
          call ReadTransformedIntegrals_readOneSpecies( i , fourCenterIntegrals(i,i)   )
       end do

       do i=1, MolecularSystem_instance%numberOfQuantumSpecies-1
          do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
             call ReadTransformedIntegrals_readTwoSpecies( i, j, fourCenterIntegrals(i,j) )
          end do
       end do

    end if

!!!Add charges
    do i=1, MolecularSystem_instance%numberOfQuantumSpecies
       fourCenterIntegrals(i,i)%values = &
            fourCenterIntegrals(i,i)%values * MolecularSystem_getCharge(i)**2.0
       do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
          fourCenterIntegrals(i,j)%values = &
               fourCenterIntegrals(i,j)%values * MolecularSystem_getCharge(i) * MolecularSystem_getCharge(j)
       end do
    end do

!!!Print integrals - debug loops
    ! do i=1, MolecularSystem_instance%numberOfQuantumSpecies
    !    print *, "transformed integrals for species", i
    !    do a=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
    !       do b=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
    !            MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
    !            MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII)) !sysII
    !          bb=b-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
    !          do c=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
    !             do d=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
    !                  MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
    !                  MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII)) !sysII
    !                dd=d-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
    !                auxIndex = fourIndexArray(i)%values(twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d) )
    !                print *, a, b, c, d, twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d), fourIndexArray(i)%values( &
    !                     twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d)), fourCenterIntegrals(i,i)%values(auxIndex, 1)
    !             end do
    !          end do
    !       end do
    !    end do
    ! end do

    ! do i=1, MolecularSystem_instance%numberOfQuantumSpecies-1
    !    do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
    !       print *, "transformed integrals for species", i, j

    !       ssize1 = max(MolecularSystem_getTotalNumberOfContractions( j ),MolecularSystem_getOcupationNumber(j))
    !       ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2
    !       do a=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
    !          do b=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
    !               MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
    !               MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII))  !sysII
    !             bb=b-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
    !             auxIndex1 = ssize1 * (twoIndexArray(i)%values(a,b) - 1 ) 
    !             do c=1, MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))  !sysI
    !                do d=MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))+1, &
    !                     MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))+&
    !                     MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))  !sysII
    !                   dd=d-MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))
    !                   auxIndex = auxIndex1  + twoIndexArray(j)%values(c,d) 
    !                   print *, a, b, c, d, auxIndex1, twoIndexArray(j)%values(c,d), auxIndex, fourCenterIntegrals(i,j)%values(auxIndex, 1)
    !                end do
    !             end do
    !          end do
    !       end do
    !    end do
    ! end do
   
    
!!!Compute Hamiltonian Matrix element between displaced geometries

    ! !!Point charge-Point charge repulsion
    ! !!One Particle Terms
    ! !!Have already been computed     
    !!Same species repulsion
    do i=1, MolecularSystem_instance%numberOfQuantumSpecies
       interactionEnergy=0.0
       do a=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
          do b=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
               MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
               MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII)) !sysII
             bb=b-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
             do c=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
                do d=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
                     MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
                     MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII)) !sysII
                   dd=d-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
                   auxIndex = fourIndexArray(i)%values(twoIndexArray(i)%values(a,b), twoIndexArray(i)%values(c,d) )
                   interactionEnergy=interactionEnergy+0.5*fourCenterIntegrals(i,i)%values(auxIndex, 1)*&
                        (inverseOverlapMatrices(i)%values(bb,a)*inverseOverlapMatrices(i)%values(dd,c)-&
                        inverseOverlapMatrices(i)%values(dd,a)*inverseOverlapMatrices(i)%values(bb,c))
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
          ssize1 = max(MolecularSystem_getTotalNumberOfContractions( j ),MolecularSystem_getOcupationNumber(j))
          ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2
          do a=1, MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI)) !sysI
             do b=MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+1, &
                  MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))+&
                  MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysII))  !sysII
                bb=b-MolecularSystem_getOcupationNumber(i,this%MolecularSystems(sysI))
                auxIndex1 = ssize1 * (twoIndexArray(i)%values(a,b) - 1 ) 
                do c=1, MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))  !sysI
                   do d=MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))+1, &
                        MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))+&
                        MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))  !sysII
                      dd=d-MolecularSystem_getOcupationNumber(j,this%MolecularSystems(sysI))
                      auxIndex = auxIndex1  + twoIndexArray(j)%values(c,d) 
                      interactionEnergy=interactionEnergy+fourCenterIntegrals(i,j)%values(auxIndex, 1)*&
                           inverseOverlapMatrices(i)%values(bb,a)*inverseOverlapMatrices(j)%values(dd,c)
                      ! print *, a, b, c, d,  fourCenterIntegrals(i,j)%values(auxIndex, 1), inverseOverlapMatrices(i)%values(bb,a), inverseOverlapMatrices(j)%values(dd,c)
                   end do
                end do
             end do
          end do
          this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)+interactionEnergy
          ! print *, "interspecies interactionEnergy for species", i, j, interactionEnergy
       end do
    end do

    deallocate(fourCenterIntegrals,twoIndexArray,fourIndexArray)

    this%configurationHamiltonianMatrix%values(sysI,sysII)=this%configurationHamiltonianMatrix%values(sysI,sysII)*&
         this%configurationOverlapMatrix%values(sysI,sysII)

    this%configurationHamiltonianMatrix%values(sysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,sysII)

  end subroutine NonOrthogonalCI_twoParticlesContributions

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

    print *, "non orthogonal CI overlap Matrix "
    call Matrix_show(this%configurationOverlapMatrix)

    print *, "non orthogonal CI Hamiltionian Matrix "
    call Matrix_show(this%configurationHamiltonianMatrix)

    print *, ""
    print *, "Transforming non orthogonal CI Hamiltionian Matrix..."
    
    call Matrix_constructor(transformationMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8) , 0.0_8)

    call Vector_constructor( eigenValues, this%numberOfDisplacedSystems )
    call Matrix_constructor( eigenVectors,int(this%numberOfDisplacedSystems,8),int(this%numberOfDisplacedSystems,8))

    !!****************************************************************
    !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
    !!          
    call Matrix_eigen( this%configurationOverlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )

    print *,"Overlap eigenvectors "
    call Matrix_show( eigenVectors )

    print *,"Overlap eigenvalues "
    call Vector_show( eigenValues )
    
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

    print *,"Matriz de transformacion "
    call Matrix_show( transformationMatrix )

    !!**********************************************************************************************
    !! Transform configuration hamiltonian matrix
    !!
    this%configurationHamiltonianMatrix%values = &
         matmul( matmul( transpose( transformationMatrix%values ) , &
         this%configurationHamiltonianMatrix%values), transformationMatrix%values )

    print *,"transformed Hamiltonian Matrix "
    call Matrix_show( this%configurationHamiltonianMatrix )

    print *, "Diagonalizing non orthogonal CI Hamiltionian Matrix..."
    !! Calcula valores y vectores propios de matriz de CI transformada.
    call Matrix_eigen( this%configurationHamiltonianMatrix, this%statesEigenvalues, this%configurationCoefficients, SYMMETRIC )

    !! Calcula los  vectores propios para matriz de CI       
    this%configurationCoefficients%values = matmul( transformationMatrix%values, this%configurationCoefficients%values )

    print *,"non orthogonal CI eigenvalues "
    call Vector_show( this%statesEigenvalues )

    print *,"configuration Coefficients"
    call Matrix_show( this%configurationCoefficients )

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
    type(MolecularSystem) :: auxMolecularSystem
    type(Matrix), allocatable :: auxCoefficients(:),mergedCoefficients(:), overlapMatrix(:)
    type(IVector), allocatable :: sysBasisList(:,:),auxBasisList(:)
    type(Matrix), allocatable :: mergedDensityMatrix(:,:)

    type(Matrix) :: auxMatrix, densityEigenVectors, auxdensityEigenVectors
    type(Vector) :: densityEigenValues, auxdensityEigenValues
    
    type(Matrix) :: plotPoints, molecularOverlapMatrix, inverseOverlapMatrix
    real(8) :: overlapDeterminant
    type(Matrix), allocatable :: orbitalsInGrid(:,:)
    type(Vector), allocatable :: densityInGrid(:)
    integer :: gridSize, state
    real(8) :: densityIntegral, auxValue
    integer :: a,b,g,i,ii,j,jj,k,p,mu,nu, sysI, sysII, speciesID

    integer :: integralsUnit, densUnit
    character(50) :: integralsFile, densFile
    character(50) :: arguments(2), auxString
    character(50) :: molFileI, outFile
    integer :: wfnUnit,wfnUnitI,wfnUnitII,outUnit
    real(8) :: timeA
    real(8) :: timeB
    logical :: existFile
    
    !$  timeA = omp_get_wtime()

    if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    allocate(mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         auxCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         sysBasisList(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies),&
         auxBasisList(molecularSystem_instance%numberOfQuantumSpecies))

    allocate(overlapMatrix(molecularSystem_instance%numberOfQuantumSpecies))
    allocate(mergedDensityMatrix(CONTROL_instance%CI_STATES_TO_PRINT,molecularSystem_instance%numberOfQuantumSpecies))
    
    !Create a super molecular system
    !!!Merge coefficients from system 1 and system 2
    call MolecularSystem_mergeTwoSystems(molecularSystem_instance, this%MolecularSystems(1), this%MolecularSystems(2), &
         sysBasisList(1,:),sysBasisList(2,:))

    call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(1,:),this%HFCoefficients(2,:),&
         this%MolecularSystems(1),this%MolecularSystems(2),&
         sysBasisList(1,:),sysBasisList(2,:),mergedCoefficients)

    ! do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
    !    print *, "2", speciesID, "ocupationNumber", MolecularSystem_getOcupationNumber(speciesID)
    !    print *, "2", speciesID, "mergedCoefficients"
    !    call Matrix_show(mergedCoefficients(speciesID))
    ! end do
    ! 
    !! Loop other systems expanding the merged coefficients matrix 
    do sysI=3, this%numberOfDisplacedSystems
       call MolecularSystem_copyConstructor(auxMolecularSystem,molecularSystem_instance)
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          call Matrix_copyConstructor(auxCoefficients(speciesID), mergedCoefficients(speciesID))       
       end do
       call MolecularSystem_mergeTwoSystems(molecularSystem_instance, auxMolecularSystem, this%MolecularSystems(sysI), &
            auxBasisList,sysBasisList(sysI,:),reorder=.false.)
       call NonOrthogonalCI_mergeCoefficients(auxCoefficients,this%HFCoefficients(sysI,:),&
            auxMolecularSystem,this%MolecularSystems(sysI),&
            auxBasisList,sysBasisList(sysI,:),mergedCoefficients)
       ! do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          !    print *, sysI, speciesID, "ocupationNumber", MolecularSystem_getOcupationNumber(speciesID)
          !    print *, sysI, speciesID, "mergedCoefficients"
          !    call Matrix_show(mergedCoefficients(speciesID))
       ! end do
    end do
    
    !!!Fix basis list size
    do sysI=1, this%numberOfDisplacedSystems
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          call Vector_copyConstructorInteger(auxBasisList(speciesID),sysBasisList(sysI,speciesID))
          call Vector_constructorInteger(sysBasisList(sysI,speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID), 0)           
          do i=1, size(auxBasisList(speciesID)%values)
             sysBasisList(sysI,speciesID)%values(i)=auxBasisList(speciesID)%values(i)
          end do
          ! print *, "sysI", sysI, "speciesID", speciesID, "after list"
          ! call Vector_showInteger(sysBasisList(sysI,speciesID))
       end do
    end do
    
    write(*,*) ""
    print *, "Superposed molecular system geometry"
    write(*,*) "---------------------------------- "
    ! call MolecularSystem_showInformation()  
    ! call MolecularSystem_showParticlesInformation()
    call MolecularSystem_showCartesianMatrix()

    ! do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
       ! write(*,*) ""
       ! write(*,*) " Merged Occupied Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
       ! write(*,*) "---------------------------------- "
       ! write(*,*) ""
       ! print *, "contractions", speciesID, int(MolecularSystem_getTotalNumberOfContractions(speciesID),8)
       ! print *, "ocupation", speciesID, int(MolecularSystem_getOcupationNumber(speciesID),8)
       ! call Matrix_constructor(auxCoefficients(speciesID),int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),&
       !      int(MolecularSystem_getOcupationNumber(speciesID),8),0.0_8)
       ! do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
       !    do j=1, MolecularSystem_getOcupationNumber(speciesID)
       !       auxCoefficients(speciesID)%values(i,j)=mergedCoefficients(speciesID)%values(i,j)
       !    end do
       ! end do
       ! call Matrix_show(auxCoefficients(speciesID))                    
    ! end do

    print *, ""
    print *, "Computing overlap integrals in for the superposed systems..."
    print *, ""
    !!Compute overlap integrals 
    call MolecularSystem_saveToFile()          

    integralsUnit = 30
    integralsFile = "lowdin.opints"
    call system("lowdin-ints.x ONE_PARTICLE >> extendedOutput.eut")

    existFile = .false.     
    inquire(file=trim(integralsFile), exist=existFile)
    if( .not. existFile ) STOP "opints does not exist"

    open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted")

    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
       arguments(1) = "OVERLAP"
       
       overlapMatrix(speciesID) = Matrix_getFromFile(rows=MolecularSystem_getTotalNumberOfContractions(speciesID), &
            columns=MolecularSystem_getTotalNumberOfContractions(speciesID), &
            unit=integralsUnit, binary=.true., arguments=arguments)

    end do
    close(integralsUnit)

    print *, ""
    print *, "Building merged density matrices for the superposed systems..."
    print *, ""
    !!Fill the merged density matrix
    do state=1, CONTROL_instance%CI_STATES_TO_PRINT
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          call Matrix_constructor(mergedDensityMatrix(state,speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8)
          ! "Diagonal" terms - same system
          do sysI=1, this%numberOfDisplacedSystems
             do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                   if(sysBasisList(sysI,speciesID)%values(mu) .ne. 0) then
                      do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                         if(sysBasisList(sysI,speciesID)%values(nu) .ne. 0) then
                            mergedDensityMatrix(state,speciesID)%values(mu,nu) =  mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                                 this%configurationCoefficients%values(sysI,state)**2*&
                                 mergedCoefficients(speciesID)%values(mu,ii)*&
                                 mergedCoefficients(speciesID)%values(nu,ii)
                         end if
                      end do
                   end if
                end do
             end do
          end do

          !!"Non Diagonal" terms - system pairs
          do sysI=1, this%numberOfDisplacedSystems
             do sysII=sysI+1, this%numberOfDisplacedSystems
             ! do sysII=1, this%numberOfDisplacedSystems
             !    if(sysI .eq. sysII) cycle
                
                if( abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) cycle

                call Matrix_constructor(molecularOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
                     int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )
                call Matrix_constructor(inverseOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
                     int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )
                
                !!Compute molecular overlap matrix and its inverse
                do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                   ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                   do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))
                      jj=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))*(sysII-1)+j
                      do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysI
                         do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                            if((sysBasisList(sysI,speciesID)%values(mu) .ne. 0) .and. (sysBasisList(sysII,speciesID)%values(nu) .ne. 0)) then
                               ! print *, "i, j, mu, nu, coefI, coefII", i,j,mu,nu,mergedCoefficients(speciesID)%values(mu,ii),mergedCoefficients(speciesID)%values(nu,jj)
                               molecularOverlapMatrix%values(i,j)=molecularOverlapMatrix%values(i,j)+&
                                    mergedCoefficients(speciesID)%values(mu,ii)*&
                                    mergedCoefficients(speciesID)%values(nu,jj)*&
                                    overlapMatrix(speciesID)%values(mu,nu)
                            end if
                         end do
                      end do
                   end do
                end do
                ! print *, "molecularOverlapMatrix sysI, sysII, speciesID", sysI, sysII, speciesID
                ! call Matrix_show(molecularOverlapMatrix)
                ! call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant,method="LU")             
                ! print *, "overlapDeterminant sysI, sysII, speciesID", sysI, sysII, speciesID, overlapDeterminant
                ! if( overlapDeterminant .gt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD ) then
                inverseOverlapMatrix=Matrix_inverse(molecularOverlapMatrix)
                ! else
                !    print *, "skipping inverse sysI, sysII, speciesID",  sysI, sysII, speciesID
                !    cycle
                ! end if
                ! Compute density contributions
                do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                   ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                   do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))
                      jj=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))*(sysII-1)+j
                      do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                         do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                            if((sysBasisList(sysI,speciesID)%values(mu) .ne. 0) .and. (sysBasisList(sysII,speciesID)%values(nu) .ne. 0)) then
                               mergedDensityMatrix(state,speciesID)%values(mu,nu) =  mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                                    this%configurationCoefficients%values(sysI,state)*&
                                    this%configurationCoefficients%values(sysII,state)*&
                                    this%configurationOverlapMatrix%values(sysI,sysII)*&
                                    inverseOverlapMatrix%values(j,i)*&
                                    mergedCoefficients(speciesID)%values(mu,ii)*&
                                    mergedCoefficients(speciesID)%values(nu,jj)
                               mergedDensityMatrix(state,speciesID)%values(nu,mu) = mergedDensityMatrix(state,speciesID)%values(nu,mu) + &
                                    this%configurationCoefficients%values(sysI,state)*&
                                    this%configurationCoefficients%values(sysII,state)*&
                                    this%configurationOverlapMatrix%values(sysI,sysII)*&
                                    inverseOverlapMatrix%values(j,i)*&
                                    mergedCoefficients(speciesID)%values(mu,ii)*&
                                    mergedCoefficients(speciesID)%values(nu,jj)
                            end if
                         end do
                      end do
                   end do
                end do
             end do
          end do
          !!symmetrize
          ! do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !    do nu = mu+1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
          !       mergedDensityMatrix(state,speciesID)%values(nu,mu) = mergedDensityMatrix(state,speciesID)%values(mu,nu)  
          !    end do
          ! end do
          ! call Matrix_show(mergedDensityMatrix(state,speciesID))

       end do
    end do
    
    !! Natural orbitals
    !! Open file - to write density matrices
    densUnit = 29
       
    densFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    open(unit = densUnit, file=trim(densFile), status="replace", form="formatted")

    write(*,*) ""
    write(*,*) "============================================="
    write(*,*) " NATURAL ORBITALS OF THE SUPERPOSED SYSTEMS: "
    write(*,*) ""

    do state=1, CONTROL_instance%CI_STATES_TO_PRINT

       write(*,*) " STATE: ", state

       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
             
          write(*,*) ""
          write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%name )
          write(*,*) "-----------------"

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

          ! print *,"Matriz de overlap "
          ! call Matrix_show( overlapMatrix(speciesID) )

          !! Lowdin orthogonalization of the density matrix
          auxMatrix = Matrix_pow( overlapMatrix(speciesID), 0.5_8 )

          auxMatrix%values=matmul(matmul(auxMatrix%values,mergedDensityMatrix(state,speciesID)%values),auxMatrix%values)
          
          print *, "Diagonalizing non orthogonal CI density Matrix..."

          !! Calcula valores y vectores propios de matriz de densidad CI ortogonal.
          call Matrix_eigen(auxMatrix , auxdensityEigenValues, auxdensityEigenVectors, SYMMETRIC )

          !! Transform back to the atomic basis
          auxMatrix = Matrix_pow( overlapMatrix(speciesID), -0.5_8 )
         
          auxdensityEigenVectors%values=matmul(auxMatrix%values,auxdensityEigenVectors%values)
          
          ! reorder and count significant occupations
          k=0
          do i = 1, MolecularSystem_getTotalNumberOfContractions(speciesID)
             densityEigenValues%values(i) =  auxdensityEigenValues%values(MolecularSystem_getTotalNumberOfContractions(speciesID) - i + 1)
             densityEigenVectors%values(:,i) = auxdensityEigenVectors%values(:,MolecularSystem_getTotalNumberOfContractions(speciesID) - i + 1)
             if(densityEigenValues%values(i) .ge. 5.0E-5 ) k=k+1
          end do
          ! Print eigenvectors with occupation larger than 5.0E-5
          call Matrix_constructor(auxMatrix,int(MolecularSystem_getTotalNumberOfContractions(speciesID),8),int(k,8),0.0_8)
          do i=1, MolecularSystem_getTotalNumberOfContractions(speciesID)
             do j=1, k
                auxMatrix%values(i,j)=densityEigenVectors%values(i,j)
             end do
          end do
          !densityEigenVectors
          call Matrix_show( auxMatrix , &
             rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
             columnkeys = string_convertvectorofrealstostring( densityEigenValues ),&
             flags=WITH_BOTH_KEYS)

          write(*,"(A10,A10,A20,I5,A15,F17.12)") "number of ", trim(MolecularSystem_getNameOfSpecie( speciesID )) ," particles in state", state , &
               " density matrix: ", sum( transpose(mergedDensityMatrix(state,speciesID)%values)*overlapMatrix(speciesID)%values)
          write(*,"(A10,A10,A40,F17.12)") "sum of ", trim(MolecularSystem_getNameOfSpecie( speciesID )) , "natural orbital occupations", sum(densityEigenValues%values)

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
          arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
          arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxString)) 

          call Matrix_writeToFile ( mergedDensityMatrix(state,speciesID), densUnit , arguments=arguments(1:2) )
          
          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%name )
          arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 
             
          call Matrix_writeToFile ( densityEigenVectors, densUnit , arguments=arguments(1:2) )

          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%name )
          arguments(1) = "OCCUPATIONS"//trim(adjustl(auxstring))

          call Vector_writeToFile( densityEigenValues, densUnit, arguments=arguments(1:2) )

        write(*,*) " End of natural orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%name )
        end do
      end do

      close(densUnit)

      write(*,*) ""
      write(*,*) " END OF NATURAL ORBITALS"
      write(*,*) "=============================="
      write(*,*) ""

    !$  timeB = omp_get_wtime()
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI density plots : ", timeB - timeA ," (s)"
    
    return
    
    ! write(*,"(A)") ""
    ! write(*,"(A)") " DENSITY PLOTS"
    ! write(*,"(A)") "========================================="
    ! write(*,"(A)") ""

    ! !Grid for density plots

    ! gridSize=1001
    ! call Matrix_constructor( plotPoints,int(gridSize**2,8),3_8)
    ! p=0
    ! do i=0,gridSize-1
    !    do j=0,0!gridSize-1
    !       do k=0,gridSize-1
    !          p=p+1
    !          plotPoints%values(p,1)=0+(-500.0+i)/250.0
    !          plotPoints%values(p,2)=0!0+(-50.0+j)/50.0
    !          plotPoints%values(p,3)=0+(-500.0+k)/250.0
    !       end do
    !    end do
    ! end do
    ! gridSize=gridSize**2

    ! !Loading molecular orbitals to memory
    ! if(allocated(orbitalsInGrid)) deallocate(orbitalsInGrid)
    ! allocate(orbitalsInGrid(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))

    !    write(auxString, '(I16)') sysI
    !    molFileI="LOWDIN-"//trim(adjustl(auxString))
    !    call MolecularSystem_loadFromFile('LOWDIN.SYS', molFileI)

    !    do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
    !       call Matrix_constructor(orbitalsInGrid(sysI,speciesID),int(MolecularSystem_getOcupationNumber(speciesID),8),&
    !            int(gridSize,8), 0.0_8)

    !       !Get atomic orbital values
    !       k=0
    !       do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
    !          do i = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)

    !             call Matrix_constructor( auxMatrix, int(gridSize,8), &
    !                  int(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital,8), 0.0_8) !orbital

    !             call ContractedGaussian_getValuesAtGrid( MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i), &
    !                  plotPoints, gridSize, auxMatrix)

    !             ! call Matrix_show(auxMatrix)

    !             do j = 1, MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(i)%numCartesianOrbital
    !                k=k+1

    !                do a=1,MolecularSystem_getOcupationNumber(speciesID)
    !                   do p = 1 , gridSize
    !                      ! print *, speciesID, k, a, coefficients(sysI,speciesID)%values(k,a), plotPoints%values(p,1:3), auxMatrix%values(p,j)
    !                      orbitalsInGrid(sysI,speciesID)%values(a,p)=orbitalsInGrid(sysI,speciesID)%values(a,p)+&
    !                           this%HFCoefficients(sysI,speciesID)%values(k,a)*auxMatrix%values(p,j)

    !                   end do
    !                end do
    !             end do

    !          end do
    !       end do
    !    end do
    ! end do


    ! !Computing NOCI densities
    ! if(allocated(densityInGrid)) deallocate(densityInGrid)
    ! allocate(densityInGrid(molecularSystem_instance%numberOfQuantumSpecies))

    ! do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies

    !    if(speciesID .ne. 3 ) cycle

    !    do state = 1, CONTROL_instance%CI_STATES_TO_PRINT

    !       call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)

    !       !Diagonal contributions
    !       do sysI=1, this%numberOfDisplacedSystems       
    !          ! call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)
    !          do a=1,MolecularSystem_getOcupationNumber(speciesID)
    !             do b=1,MolecularSystem_getOcupationNumber(speciesID)
    !                do p = 1 , gridSize
    !                   densityInGrid(speciesID)%values(p)=densityInGrid(speciesID)%values(p)+&
    !                        this%configurationCoefficients%values(sysI,state)**2.0*&
    !                        orbitalsInGrid(sysI,speciesID)%values(a,p)*&
    !                        orbitalsInGrid(sysI,speciesID)%values(b,p)
    !                end do
    !             end do
    !          end do
    !       end do

    !       !Off Diagonal contributions
    !       do sysI=1, this%numberOfDisplacedSystems       
    !          do sysII=sysI+1, this%numberOfDisplacedSystems !sysI+1
    !             ! call Vector_constructor(densityInGrid(speciesID), gridSize, 0.0_8)
    !             do a=1,MolecularSystem_getOcupationNumber(speciesID)
    !                do b=1,MolecularSystem_getOcupationNumber(speciesID)
    !                   do p = 1 , gridSize
    !                      ! densityInGrid(speciesID)%values(p)=densityInGrid(speciesID)%values(p)+&
    !                      !      2.0*this%configurationCoefficients%values(sysI,state)*&
    !                      !      this%configurationCoefficients%values(sysII,state)*&
    !                      !      inverseOverlapMatrices(speciesID)%values(b,a)*&
    !                      !      this%configurationOverlapMatrix%values(sysI,sysII)*&
    !                      !      orbitalsInGrid(sysI,speciesID)%values(a,p)*&
    !                      !      orbitalsInGrid(sysII,speciesID)%values(b,p)
    !                   end do
    !                end do
    !             end do
    !          end do
    !       end do

    !       !density check
    !       ! densityIntegral=0.0
    !       ! do p = 1 , gridSize
    !       ! densityIntegral=densityIntegral+densityInGrid(speciesID)%values(p)*0.02**3.0   
    !       ! end do
    !       ! print *, "densityIntegral", densityIntegral

    !       !1D plot
    !       outUnit=700
    !       write(auxString, '(I16)') state
    !       outFile=trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//"-"//trim(adjustl(auxString))//".2D.dens"

    !       print *, "printing...", outFile

    !       open(unit=outUnit, file=trim(outFile), status="replace", form="formatted")          
    !       write (outUnit,*) "# Z dens for speciesID", speciesID
    !       do p = 1 , gridSize
    !          if(plotPoints%values(p,1) .eq. 0.0 .and. plotPoints%values(p,2) .eq. 0.0) &
    !               write (outUnit,"(T10,F20.8,E20.8)") plotPoints%values(p,3), densityInGrid(speciesID)%values(p)
    !       end do
    !       close(outUnit)

    !       !2D plot
    !       outFile=trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//"-"//trim(adjustl(auxString))//".3D.dens"
    !       print *, "printing...", outFile

    !       open(unit=outUnit, file=trim(outFile), status="replace", form="formatted")          
    !       write (outUnit,*) "# X Z dens for speciesID", speciesID

    !       auxValue=plotPoints%values(1,1)
    !       do p = 1 , gridSize
    !          if(plotPoints%values(p,2) .eq. 0.0) then
    !             if(plotPoints%values(p,1) .gt. auxValue ) then             
    !                write (outUnit,"(T10)")
    !                auxValue=plotPoints%values(p,1)
    !             end if
    !             write (outUnit,"(T10,F20.8,F20.8,E20.8)") plotPoints%values(p,1), plotPoints%values(p,3), densityInGrid(speciesID)%values(p)
    !          end if
    !       end do

    !    end do
    ! end do

  end subroutine NonOrthogonalCI_plotDensities

  !>
  !! @brief Calculate and Transform the four center integrals in one sweep without writing anything to disk
  !!
  !! @param molecularSystem, HFCoefficients: species array with the atomic coefficients, fourCenterIntegrals: species*species array to save integrals
  !<
  subroutine NonOrthogonalCI_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, twoIndexArray, fourIndexArray, fourCenterIntegrals)
    implicit none
    type(MolecularSystem), intent(in) :: mergedMolecularSystem
    type(Matrix), intent(in) :: mergedCoefficients(mergedMolecularSystem%numberOfQuantumSpecies)
    type(iMatrix), intent(in) :: twoIndexArray(mergedMolecularSystem%numberOfQuantumSpecies)
    type(iMatrix), intent(in) :: fourIndexArray(mergedMolecularSystem%numberOfQuantumSpecies)
    type(Matrix), intent(out) :: fourCenterIntegrals(mergedMolecularSystem%numberOfQuantumSpecies,mergedMolecularSystem%numberOfQuantumSpecies)

    real(8), allocatable, target :: auxtempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer :: p, p_l, p_u
    integer :: q, q_l, q_u
    integer :: r, r_l, r_u
    integer :: s, s_l, s_u
    integer :: ssize, ssizeb, auxIndex, auxIndexA
    integer :: n,u, nu, lambda,sigma
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
       
       if ( allocated (tempB)) deallocate (tempB )
       allocate (tempB ( ssize, ssize ) )
       tempB = 0

       if ( allocated (tempC)) deallocate (tempC )
       allocate (tempC ( ssize ) )
       tempC = 0

       do p = p_l, p_u
          n = p

          !First quarter transformation happens here
          call DirectIntegralManager_getDirectIntraRepulsionFirstQuarter(&
               speciesID, &
               trim(CONTROL_instance%INTEGRAL_SCHEME), &
               densityMatrix, & 
               mergedCoefficients(speciesID), &
               auxtempA, p-1 )

          do q = p, q_u
             u = q
             tempB = 0

             if ( q < q_l ) cycle
             !! second quarter
             do nu = 1, ssize
                !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
                tempB(:,:) = tempB(:,:) + mergedCoefficients(speciesID)%values( nu, q )* &
                     auxtempA(:,:,nu)
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
          
          if ( allocated (tempB)) deallocate (tempB )
          allocate (tempB ( ssizeb, ssizeb ) )
          tempB = 0

          if ( allocated (tempC)) deallocate (tempC )
          allocate (tempC ( ssizeb ) )
          tempC = 0

          do p = p_l, p_u

             !First quarter transformation happens here
             call DirectIntegralManager_getDirectInterRepulsionFirstQuarter(&
                  speciesID, otherSpeciesID, &
                  trim(CONTROL_instance%INTEGRAL_SCHEME), &
                  densityMatrix, & 
                  mergedCoefficients(speciesID), &
                  auxtempA, p-1 )

             do q = q_l, q_u
                tempB = 0

                if ( q < p ) cycle
                !! second quarter
                do nu = 1, ssize
                   !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle

                   tempB(:,:) = tempB(:,:) + mergedCoefficients(speciesID)%values( nu, q )* &
                        auxtempA(:,:,nu)
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

                      if ( s < r ) cycle
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

    call DirectIntegralManager_destructor()
    
  end subroutine NonOrthogonalCI_transformIntegralsMemory

  
end module NonOrthogonalCI_

