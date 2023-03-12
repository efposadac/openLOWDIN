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
  use Libint2Interface_
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

  type, public :: NonOrthogonalCI
     logical :: isInstanced
     integer :: numberOfDisplacedSystems
     integer :: numberOfEnergyRejectedSystems
     integer :: numberOfEllipsoidRejectedSystems
     integer :: numberOfPPdistanceRejectedSystems
     integer :: numberOfNPdistanceRejectedSystems
     integer :: numberOfEquivalentSystems
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
       NonOrthogonalCI_runHFs,&
       NonOrthogonalCI_buildOverlapAndHamiltonianMatrix,&
       NonOrthogonalCI_diagonalizeCImatrix,&
       NonOrthogonalCI_generateDensities

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
    this%numberOfEnergyRejectedSystems=0
    this%numberOfEllipsoidRejectedSystems=0
    this%numberOfPPdistanceRejectedSystems=0
    this%numberOfNPdistanceRejectedSystems=0
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
       if(MolecularSystem_instance%allParticles(p)%particlePtr%rotationPoint.eq.0) cycle
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
                if ( allocated(MolecularSystem_instance%allParticles(q)%particlePtr%childs) ) then
                   do r=1,size(MolecularSystem_instance%allParticles(q)%particlePtr%childs)
                      this%rotationCenterList( MolecularSystem_instance%allParticles(q)%particlePtr%childs(r),1)=numberOfRotationCenters
                   end do
                end if
             end if
             this%rotationCenterList(q,2)=p
          end if
       end do
    end do
    print *, ""

    ! print *, "this%rotationCenterList" 
    ! do p=1, size(MolecularSystem_instance%allParticles)
    !    print *, "Particle ", trim(ParticleManager_getSymbol(p)),this%rotationCenterList(p,1), this%rotationCenterList(p,2)
    ! end do
    
    if(numberOfTranslationCenters.ne.0) then

       this%transformationType="TRANSLATION"
       this%numberOfTransformedCenters=numberOfTranslationCenters
       this%numberOfIndividualTransformations=&
            CONTROL_instance%TRANSLATION_SCAN_GRID(1)*CONTROL_instance%TRANSLATION_SCAN_GRID(2)*CONTROL_instance%TRANSLATION_SCAN_GRID(3)&
            +(CONTROL_instance%TRANSLATION_SCAN_GRID(1)-1)*(CONTROL_instance%TRANSLATION_SCAN_GRID(2)-1)*(CONTROL_instance%TRANSLATION_SCAN_GRID(3)-1)

              print *, ""
       write (*,"(A,I5,A,I10,A)") "Displacing coordinates of ", numberOfTranslationCenters, " centers", &
            this%numberOfIndividualTransformations," times"
       print *, ""

    else if(numberOfRotationCenters.ne.0) then
       print *, ""
       write (*,"(A,I5,A,I5,A,I5,A)") "Rotating coordinates of ", numberOfRotationCenters, " centers", CONTROL_instance%ROTATIONAL_SCAN_GRID, &
            " times in ", CONTROL_instance%NESTED_ROTATIONAL_GRIDS, " nested grids"
            
       print *, ""

       this%transformationType="ROTATION"
       this%numberOfTransformedCenters=numberOfRotationCenters
       this%numberOfIndividualTransformations=CONTROL_instance%ROTATIONAL_SCAN_GRID*CONTROL_instance%NESTED_ROTATIONAL_GRIDS

    end if

    !!Dynamically allocated through the displacement routine
    allocate(this%MolecularSystems(0))

    ! call Vector_constructorInteger(this%systemTypes,this%numberOfIndividualTransformations**this%numberOfTransformedCenters,0)
    
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
    type(MolecularSystem) :: displacedMolecularSystem
    real(8) :: displacement
    character(100) :: coordsFile
    integer, allocatable :: transformationCounter(:)
    integer :: coordsUnit
    integer :: i,j
    integer :: sysI, speciesID
    integer :: closestSystem
    integer :: systemType
    logical :: newSystemFlag
    logical :: skip
    real(8) :: timeA
    
    !$  timeA = omp_get_wtime()
    
    call MolecularSystem_copyConstructor(originalMolecularSystem, molecularSystem_instance)

    allocate(transformationCounter(this%numberOfTransformedCenters))

    transformationCounter(1:this%numberOfTransformedCenters)=1
    transformationCounter(1)=0

    this%numberOfDisplacedSystems=0

    coordsUnit=333
    coordsFile=trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords"

    print *, "generating NOCI displaced geometries and HF wavefunctions... saving coords to ", trim(coordsFile)

    open(unit=coordsUnit, file=trim(coordsFile), status="replace", form="formatted")

!!!!! clock type iterations to form all the possible combination of modified geometries
    do while (.true.)

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
       
       write (coordsUnit,"(A)",advance="no") "Transformation counter: "
       do i=1,this%numberOfTransformedCenters
          write (coordsUnit,"(I10)",advance="no") transformationCounter(i)
       end do
       write (coordsUnit,*) ""

       skip=.false.
       !Apply the transformation given by transformationCounter to each center, the result is saved in molecularSystemInstance
       call NonOrthogonalCI_transformCoordinates(this,transformationCounter(1:this%numberOfTransformedCenters),originalMolecularSystem,displacedMolecularSystem,skip)       
       
       write(displacedMolecularSystem%description, '(10I6)') transformationCounter(:)
       
       call MolecularSystem_showCartesianMatrix(displacedMolecularSystem,unit=coordsUnit)

       !Classify the system according to its distance matrix (symmetry) 
       if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
            call NonOrthogonalCI_classifyNewSystem(this,systemType, newSystemFlag) 

       !Check if the new system is not beyond the max displacement
       if(skip) then
          write (coordsUnit,"(A)") "Skipping system beyond the ellipsoids boundaries"
          this%numberOfEllipsoidRejectedSystems=this%numberOfEllipsoidRejectedSystems+1                      
          cycle
       end if

       !Check if the separation between particles of the same charge is not too small
       call NonOrthogonalCI_checkSameChargesDistance(displacedMolecularSystem,displacement,skip)

       if(skip) then
          write (coordsUnit,"(A,F20.12)") "Skipping system with same charge particle separation", displacement
          this%numberOfPPdistanceRejectedSystems=this%numberOfPPdistanceRejectedSystems+1                      
          cycle
       end if

       !Check if the separation between positive and negative particles is not too big
       call NonOrthogonalCI_checkOppositeChargesDistance(displacedMolecularSystem,displacement,skip)

       if(skip) then
          write (coordsUnit,"(A,F20.12)") "Skipping system with positive and negative particle separation", displacement
          this%numberOfNPdistanceRejectedSystems=this%numberOfNPdistanceRejectedSystems+1                      
          cycle
       end if
       
       !Check if the new system is not to close to previous calculated systems - duplicate protection
       call NonOrthogonalCI_checkNewSystemDisplacement(this,displacedMolecularSystem,closestSystem,displacement) 

       if(displacement .lt. CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) then
          write (coordsUnit,"(A,F20.12,A,I10)") "Skipping system with distance ", displacement , "a.u. from system ", closestSystem
          skip=.true.
          this%numberOfEquivalentSystems=this%numberOfEquivalentSystems+1                      
          cycle
       end if

       !!Copy the molecular system to the NonOrthogonalCI object
       ! if(newSystemFlag) then
       !    this%numberOfUniqueSystems=this%numberOfUniqueSystems+1
       !    this%systemTypes%values(this%numberOfDisplacedSystems)=this%numberOfUniqueSystems
       ! else
       !    this%systemTypes%values(this%numberOfDisplacedSystems)=systemType
       ! end if

       ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) then
       !    write (coordsUnit,"(A,I5,A,I10,A,F20.12)") "Saving system of type ",  this%systemTypes%values(this%numberOfDisplacedSystems) , &
       !         " with ID ", this%numberOfDisplacedSystems, " and energy", testEnergy
       ! else<
       if(skip .eqv. .false.) then
          call NonOrthogonalCI_saveSystem(this,displacedMolecularSystem)
          write (coordsUnit,"(A,I10)") "Saving system with ID ", this%numberOfDisplacedSystems
       end if
    end do

    close(coordsUnit)

    print *, ""
    write (*,'(A10,I10,A)') "Mixing ", this%numberOfDisplacedSystems, " HF calculations at different geometries"

    if(this%numberOfEllipsoidRejectedSystems .gt. 0) &
         write (*,'(A10,I10,A)') "Rejected ", this%numberOfEllipsoidRejectedSystems, &
         " geometries outside the ellipsoids area"

    if(this%numberOfPPdistanceRejectedSystems .gt. 0) &
         write (*,'(A10,I10,A,F18.12,A,F18.12)') "Rejected ", this%numberOfPPdistanceRejectedSystems, &
         " geometries with separation between same charge basis sets smaller than", CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE, &
         " or larger than", CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE

    if(this%numberOfNPdistanceRejectedSystems .gt. 0) &
         write (*,'(A10,I10,A,F18.12)') "Rejected ", this%numberOfNPdistanceRejectedSystems, &
         " geometries with separation between positive and negative basis sets larger than", CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE

    if(this%numberOfEquivalentSystems .gt. 0) &
         write (*,'(A10,I10,A)') "Rejected ", this%numberOfEquivalentSystems, &
         " duplicated geometries after permutations"
    
    print *, ""
    
    call Matrix_constructor(this%configurationHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Matrix_constructor(this%configurationOverlapMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)

    if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
         call Matrix_constructorInteger(this%configurationPairTypes,int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0)
    

    ! minEnergy=0.0    
    
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time displacing coordinates : ", omp_get_wtime() - timeA ," (s)"
    print *, ""
    
  end subroutine NonOrthogonalCI_displaceGeometries


  !>
  !! @brief Apply the transformation (translation or rotation) given by transformationCounter to each center, based in the originalMolecularSystemPositions the result is saved in molecularSystemInstance 
  !! @param this,transformationCounter,originalMolecularSystem
  !<
  subroutine NonOrthogonalCI_transformCoordinates(this,transformationCounter,originalMolecularSystem,displacedMolecularSystem,skip)
    type(NonOrthogonalCI) :: this
    integer :: transformationCounter(*)
    type(MolecularSystem) :: originalMolecularSystem
    type(MolecularSystem), target :: displacedMolecularSystem
    logical, intent(out) :: skip

    real(8) :: centerX, centerY, centerZ, displacedOrigin(3), distanceCheck, distanceToCenter
    integer :: center, displacementId
    real(8),allocatable :: X(:), Y(:), Z(:), W(:)
    integer :: i,j,k,p,q,mu
    
    skip=.false.

    call MolecularSystem_copyConstructor(displacedMolecularSystem, originalMolecularSystem)
    
    displacedMolecularSystem%description=""
    do i=1,this%numberOfTransformedCenters
       write(MolecularSystem_instance%description, '(A,I6)') trim(originalMolecularSystem%description), transformationCounter(i)
    end do

    particleManager_instance => displacedMolecularSystem%allParticles
    
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
          !Body centered cube
          do i=1,CONTROL_instance%TRANSLATION_SCAN_GRID(1)*2-1
             do j=1,CONTROL_instance%TRANSLATION_SCAN_GRID(2)*2-1
                do k=1,CONTROL_instance%TRANSLATION_SCAN_GRID(3)*2-1

                   if( (mod(i,2) .eq. mod(j,2)) .and. (mod(i,2) .eq. mod(k,2)) ) then
                      displacementId=displacementId+1

                      if(displacementId .eq. transformationCounter(center) ) then

                         distanceCheck= &
                              (CONTROL_instance%TRANSLATION_STEP*((i+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT(1)**2+&
                              (CONTROL_instance%TRANSLATION_STEP*((j+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT(2)**2+&
                              (CONTROL_instance%TRANSLATION_STEP*((k+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT(3)**2

                         if(distanceCheck .gt. 1.0) then
                            skip=.true.
                            ! return
                         end if
                         
                         distanceCheck= &
                              (CONTROL_instance%TRANSLATION_STEP*((i+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT(1)**2+&
                              (CONTROL_instance%TRANSLATION_STEP*((j+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT(2)**2+&
                              (CONTROL_instance%TRANSLATION_STEP*((k+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0))**2/&
                              CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT(3)**2

                         if(distanceCheck .lt. 1.0) then
                            skip=.true.
                            ! return
                         end if

                         displacedOrigin(1)=centerX+CONTROL_instance%TRANSLATION_STEP*((i+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0)
                         displacedOrigin(2)=centerY+CONTROL_instance%TRANSLATION_STEP*((j+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0)
                         displacedOrigin(3)=centerZ+CONTROL_instance%TRANSLATION_STEP*((k+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0)
                         
                         do p=1, size(displacedMolecularSystem%allParticles)
                            if(center.eq.displacedMolecularSystem%allParticles(p)%particlePtr%translationCenter) then
                               ! call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                               displacedMolecularSystem%allParticles(p)%particlePtr%origin=displacedOrigin
                               do mu = 1, displacedMolecularSystem%allParticles(p)%particlePtr%basis%length
                                  displacedMolecularSystem%allParticles(p)%particlePtr%basis%contraction(mu)%origin = displacedOrigin
                               end do
                            end if
                         end do
                         
                         ! write(*, '(3I5,F4.1,A,F4.1,A,F4.1)') i,j,k, &
                         !      (i+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(1)+1)/2.0," ", &
                         !      (j+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(2)+1)/2.0," ", &
                         !      (k+1)/2.0-(CONTROL_instance%TRANSLATION_SCAN_GRID(3)+1)/2.0
                      end if
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
                   do p=1, size(displacedMolecularSystem%allParticles)
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
                         displacedMolecularSystem%allParticles(p)%particlePtr%origin=displacedOrigin
                         do mu = 1, displacedMolecularSystem%allParticles(p)%particlePtr%basis%length
                            displacedMolecularSystem%allParticles(p)%particlePtr%basis%contraction(mu)%origin = displacedOrigin
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
  subroutine NonOrthogonalCI_checkNewSystemDisplacement(this,newMolecularSystem,closestSystem,displacement) 
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: newMolecularSystem
    integer :: closestSystem
    real(8) :: displacement

    integer :: sysI, i
    type(Vector), allocatable :: displacementVector(:)
    real(8) :: dispSum
    
    displacement=1.0E8
    
    allocate(displacementVector(newMolecularSystem%numberOfQuantumSpecies))
        
    do sysI=1, this%numberOfDisplacedSystems

       call MolecularSystem_GetTwoSystemsDisplacement(this%MolecularSystems(sysI), newMolecularSystem, displacementVector)

       dispSum=0.0
       do i=1, newMolecularSystem%numberOfQuantumSpecies
          dispSum=dispSum+sum(displacementVector(i)%values(:))
       end do
       if(dispSum .lt. displacement ) then
          displacement=dispSum
          closestSystem=sysI
          if(displacement .lt. CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) exit
       end if
    end do

    deallocate(displacementVector)
    
  end subroutine NonOrthogonalCI_checkNewSystemDisplacement


  !>
  !! @brief Finds the maximum of the distances between the basis set center of a particle to its closest neighbour with opposite charge
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NonOrthogonalCI_checkOppositeChargesDistance(molSys,minNPDistance,skip) 
    implicit none
    type(MolecularSystem) :: molSys
    real(8) :: minNPDistance
    logical :: skip

    integer :: p,q
    real(8) :: npDistance
    

    minNPDistance=1E8
    do p=1, size(molSys%allParticles)-1
       if(.not.(molSys%allParticles(p)%particlePtr%translationCenter .ne. 0 .or. &
            molSys%allParticles(p)%particlePtr%rotateAround .ne. 0)) cycle
       do q=p+1, size(molSys%allParticles)
          if(.not.(molSys%allParticles(q)%particlePtr%translationCenter .ne. 0 .or. &
               molSys%allParticles(q)%particlePtr%rotateAround .ne. 0)) cycle
          if( molSys%allParticles(p)%particlePtr%charge*molSys%allParticles(q)%particlePtr%charge .gt. 0.0 ) cycle
          npDistance=sqrt(&
               (molSys%allParticles(p)%particlePtr%origin(1)-molSys%allParticles(q)%particlePtr%origin(1))**2+&
               (molSys%allParticles(p)%particlePtr%origin(2)-molSys%allParticles(q)%particlePtr%origin(2))**2+&
               (molSys%allParticles(p)%particlePtr%origin(3)-molSys%allParticles(q)%particlePtr%origin(3))**2)
          if(npDistance .lt. minNPDistance) minNPDistance=npDistance
       end do
    end do

    if(minNPDistance .gt. CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE) skip=.true.
    
  end subroutine NonOrthogonalCI_checkOppositeChargesDistance

  !>
  !! @brief Finds the maximum of the distances between the basis set center of a particle to its closest neighbour with the same charge
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NonOrthogonalCI_checkSameChargesDistance(molSys,distance,skip) 
    implicit none
    type(MolecularSystem) :: molSys
    real(8) :: distance
    logical :: skip
    
    real(8) :: minPPDistance

    integer :: p,q
    real(8) :: ppDistance
    

    minPPDistance=1.0E8
    do p=1, size(molSys%allParticles)-1
       if(.not.(molSys%allParticles(p)%particlePtr%translationCenter .ne. 0 .or. &
            molSys%allParticles(p)%particlePtr%rotateAround .ne. 0)) cycle
       do q=p+1, size(molSys%allParticles)
          if(.not.(molSys%allParticles(q)%particlePtr%translationCenter .ne. 0 .or. &
               molSys%allParticles(q)%particlePtr%rotateAround .ne. 0)) cycle
          if( molSys%allParticles(p)%particlePtr%charge*molSys%allParticles(q)%particlePtr%charge .lt. 0.0 ) cycle

          ppDistance=sqrt(&
               (molSys%allParticles(p)%particlePtr%origin(1)-molSys%allParticles(q)%particlePtr%origin(1))**2+&
               (molSys%allParticles(p)%particlePtr%origin(2)-molSys%allParticles(q)%particlePtr%origin(2))**2+&
               (molSys%allParticles(p)%particlePtr%origin(3)-molSys%allParticles(q)%particlePtr%origin(3))**2)
          if(ppDistance .lt. minPPDistance) minPPDistance=ppDistance

       end do
    end do

    if(minPPDistance .gt. CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE) skip=.true.
    if(minPPDistance .lt. CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE) skip=.true.

  end subroutine NonOrthogonalCI_checkSameChargesDistance
  
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
  !! @brief Run a Hartree-Fock calculation at displaced geometries and fill CI matrix diagonals 
  !!
  !! @param this -> NOCI instance
  !<
  subroutine NonOrthogonalCI_runHFs(this)
    implicit none
    type(NonOrthogonalCI) :: this

    character(50) :: wfnFile
    character(50) :: arguments(2)
    integer :: wfnUnit
    integer :: sysI, speciesID
    real(8) :: timeA

    !$  timeA = omp_get_wtime()
    !!Read HF energy of the non displaced SCF calculation 
    ! print *, "HF reference energy is ", hfEnergy

    allocate(this%HFCoefficients(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    allocate(this%systemLabels(this%numberOfDisplacedSystems))
        
    print *, "running HF calculations at the displaced geometries ..."
    
    do sysI=1, this%numberOfDisplacedSystems
       write(this%systemLabels(sysI), '(A)') trim(this%MolecularSystems(sysI)%description)

       !!Do SCF without calling lowdin-scf.x
       call MolecularSystem_copyConstructor(molecularSystem_instance, this%MolecularSystems(sysI))

       allocate(WaveFunction_instance(molecularSystem_instance%numberOfQuantumSpecies))
       
       call MultiSCF_constructor(MultiSCF_instance,WaveFunction_instance,CONTROL_instance%ITERATION_SCHEME)
       MultiSCF_instance%printSCFiterations=.false.

       call MultiSCF_buildHcore(MultiSCF_instance,WaveFunction_instance)

       call MultiSCF_getInitialGuess(MultiSCF_instance,WaveFunction_instance)

       call MultiSCF_solveHartreeFockRoothan(MultiSCF_instance,WaveFunction_instance,Libint2Instance)
       
       this%configurationHamiltonianMatrix%values(sysI,sysI)=MultiSCF_instance%totalEnergy
       do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
          this%HFCoefficients(sysI,speciesID) = WaveFunction_instance(speciesID)%waveFunctionCoefficients
       end do

       call DirectIntegralManager_destructor(Libint2Instance)
       call MultiSCF_destructor(MultiSCF_instance)
       deallocate(WaveFunction_instance)
       
       !!Screen geometries with high energies
       ! if( CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD .ne. 0.0 .and. &
       !      testEnergy .gt. this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD) then
       !    write (coordsUnit,"(A,F20.12)") "Skipping system with high energy", testEnergy
       !    this%numberOfEnergyRejectedSystems=this%numberOfEnergyRejectedSystems+1                      
       ! else
       ! if(this%numberOfEnergyRejectedSystems .gt. 0) &
       !      write (*,'(A10,I10,A,F18.12)') "Rejected ", this%numberOfEnergyRejectedSystems, &
       !      " geometries with energy higher than", this%refEnergy+CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD       
       
    end do
    
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for HF calculations at displaced geometries : ", omp_get_wtime() - timeA ," (s)"
    
  end subroutine NonOrthogonalCI_runHFs

  ! >
  ! @brief Saves molecular system and wfn files for a displaced system 
  
  ! @param systemID
  ! <
  subroutine NonOrthogonalCI_saveSystem(this, newSystem)
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: newSystem
    
    type(MolecularSystem), allocatable :: tempMolecularSystems(:)
    integer :: i

    !!Increase the size of the molecular systems array by 1
    this%numberOfDisplacedSystems=this%numberOfDisplacedSystems+1

    allocate(tempMolecularSystems(size(this%MolecularSystems)))

    do i=1, size(this%MolecularSystems)
       call MolecularSystem_copyConstructor(tempMolecularSystems(i),this%MolecularSystems(i))
    end do

    deallocate(this%MolecularSystems)
    allocate(this%MolecularSystems(this%numberOfDisplacedSystems))

    do i=1, size(tempMolecularSystems)
       call MolecularSystem_copyConstructor(this%MolecularSystems(i),tempMolecularSystems(i))    
    end do

    deallocate(tempMolecularSystems)
    !!Copy the molecular system to the NonOrthogonalCI object

    call MolecularSystem_copyConstructor(this%MolecularSystems(this%numberOfDisplacedSystems), newSystem)    
        
  end subroutine NonOrthogonalCI_saveSystem

  !>
  !! @brief Computes overlap and hamiltonian non orthogonal CI matrices for previously calculated molecular systems at different geometries
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_buildOverlapAndHamiltonianMatrix(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem), allocatable :: mergedMolecularSystem(:)
    type(Libint2Interface), allocatable :: Libint2ParallelInstance(:,:)
    integer, allocatable :: sysIIbatch(:)
    integer :: sysI,sysII,me,mySysII, preSysI, preSysII
    type(Matrix), allocatable :: mergedCoefficients(:), inverseOverlapMatrices(:)
    type(IVector), allocatable :: sysIbasisList(:,:),sysIIbasisList(:,:)
    real(8) :: overlapUpperBound
    integer :: n
    integer :: prescreenedElements, overlapScreenedElements
    logical :: newPairFlag

    integer :: nspecies
    integer :: ncores, batchSize

    integer :: matrixUnit
    character(100) :: matrixFile
    
    real(8) :: timeMerging, timePrescreen, timeSymmetry, timeOverlap, timeTwoIntegrals
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

    !Allocate objets to distribute in parallel
    nspecies=molecularSystem_instance%numberOfQuantumSpecies 
    ncores=OMP_get_max_threads()
    batchSize=ncores*10
    print *, "ncores", ncores, "batchsize", batchSize
    
    allocate(mergedMolecularSystem(batchSize),&
         mergedCoefficients(nspecies),&
         inverseOverlapMatrices(nspecies),&
         Libint2ParallelInstance(nspecies,batchSize),&
         sysIIbatch(batchSize),&
         sysIbasisList(nspecies,batchSize),&
         sysIIbasisList(nspecies,batchSize))
    
    systemI: do sysI=1, this%numberOfDisplacedSystems       

       this%configurationOverlapMatrix%values(sysI,sysI)=1.0
       !Save diagonal elements
       write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, sysI, &
            this%configurationOverlapMatrix%values(sysI,sysI), this%configurationHamiltonianMatrix%values(sysI,sysI)               

       sysII=sysI
       systemII: do while(sysII.lt.this%numberOfDisplacedSystems)

          ! print *, "distributing sysII", sysII, "into", batchSize, "batches"
          
          !In serial, prepare systems
          sysIIbatch(:)=0
          me=0
          mySysII=sysII
          do while(me.lt.batchSize)
             mySysII=mySysII+1
             if(mySysII .gt. this%numberOfDisplacedSystems) exit

             !$  timeA = omp_get_wtime()
             !Estimates overlap with a 1s-1s integral approximation
             call NonOrthogonalCI_prescreenOverlap(this,sysI,mySysII,overlapUpperBound)

             if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
                  overlapUpperBound .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
                ! this%configurationOverlapMatrix%values(sysI,mySysII)=0.0
                ! this%configurationOverlapMatrix%values(mySysII,sysI)=0.0
                ! this%configurationHamiltonianMatrix%values(sysI,mySysII)=0.0
                ! this%configurationHamiltonianMatrix%values(mySysII,sysI)=0.0
                ! print *, "preskipping elements", sysI, mySysII, "with overlap estimated as", overlapUpperBound
                prescreenedElements=prescreenedElements+1
             else
                !$  timeB = omp_get_wtime()
                !$  timePrescreen=timePrescreen+(timeB - timeA)
                me=me+1
                sysIIbatch(me)=mySysII
                !$  timeA = omp_get_wtime()
                !This generates a new molecular system
                ! print *, "Merging systems from geometries ", sysI, mySysII
                call MolecularSystem_mergeTwoSystems(mergedMolecularSystem(me), &
                     this%MolecularSystems(sysI), this%MolecularSystems(mySysII),sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me))
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
          !$omp& private(mySysII,mergedCoefficients,inverseOverlapMatrices),&
          !$omp& shared(this,sysI,sysII,matrixUnit,prescreenedElements,overlapScreenedElements,sysIbasisList,sysIIbasisList,mergedMolecularSystem,Libint2ParallelInstance,nspecies,batchSize)
          !$omp do schedule(dynamic,10)
          procs: do me=1, batchSize
             mySysII=sysIIbatch(me)

             ! if(mySysII .gt. this%numberOfDisplacedSystems) cycle procs
             if(mySysII .eq. 0) cycle procs

             ! print *, "evaluating S and H elements for", sysI, mySysII

             ! cycle systemII
             !! Merge occupied coefficients into a single matrix 
             call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(sysI,:),this%HFCoefficients(mySysII,:),&
                  this%MolecularSystems(sysI),this%MolecularSystems(mySysII),mergedMolecularSystem(me),&
                  sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me),mergedCoefficients)
             !$  timeA = omp_get_wtime()

             call NonOrthogonalCI_computeOverlapAndHCoreElements(this,sysI,mySysII,mergedMolecularSystem(me),mergedCoefficients,&
                  sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me),inverseOverlapMatrices)
             !$  timeB = omp_get_wtime()
             !$  timeOverlap=timeOverlap+(timeB - timeA)

             !! SKIP ENERGY EVALUATION IF OVERLAP IS TOO LOW

             if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
                  abs(this%configurationOverlapMatrix%values(sysI,mySysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
                ! print *, "screening elements", sysI, mySysII, "with overlap", this%configurationOverlapMatrix%values(sysI,mySysII)
                this%configurationOverlapMatrix%values(sysI,mySysII)=0.0
                this%configurationOverlapMatrix%values(mySysII,sysI)=0.0
                this%configurationHamiltonianMatrix%values(sysI,mySysII)=0.0
                this%configurationHamiltonianMatrix%values(mySysII,sysI)=0.0
                !$OMP ATOMIC
                overlapScreenedElements=overlapScreenedElements+1
                ! cycle systemII
             else
                !$  timeA = omp_get_wtime()
                call NonOrthogonalCI_twoParticlesContributions(this,sysI,mySysII,mergedMolecularSystem(me),&
                     inverseOverlapMatrices,mergedCoefficients,Libint2ParallelInstance(1:nspecies,me))
                !$  timeB = omp_get_wtime()
                !$  timeTwoIntegrals=timeTwoIntegrals+(timeB - timeA)
             end if

             ! print *, "thread", omp_get_thread_num()+1,"me", me, "sysI", " mySysII", sysI, mySysII, "S", this%configurationOverlapMatrix%values(sysI,mySysII), "H", this%configurationHamiltonianMatrix%values(sysI,mySysII)
          end do procs
          !$omp end do nowait
          !$omp end parallel

          !In serial, free memory and print
          do me=1, batchSize
             mySysII=sysIIbatch(me)
             ! if(mySysII .gt. this%numberOfDisplacedSystems) exit
             if(mySysII .eq. 0) exit systemII
             write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, mySysII, &
                  this%configurationOverlapMatrix%values(sysI,mySysII), this%configurationHamiltonianMatrix%values(sysI,mySysII)               
             call DirectIntegralManager_destructor(Libint2ParallelInstance(1:nspecies,me))
          end do
          sysII=mySysII
          
       end do systemII
    end do systemI

    close(matrixUnit)

    print *, ""
    print *, "Configuration pairs skipped by overlap prescreening: ", prescreenedElements
    print *, "Configuration pairs skipped by overlap    screening: ", overlapScreenedElements
    print *, "Overlap integrals computed for    ", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
         -prescreenedElements, "configuration pairs"
    print *, "Four center integrals computed for", this%numberOfDisplacedSystems*(this%numberOfDisplacedSystems-1)/2&
         -prescreenedElements-overlapScreenedElements, "configuration pairs"
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
          !    call NonOrthogonalCI_classifyConfigurationPair(this,sysI,sysII,newPairFlag)
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

  end subroutine NonOrthogonalCI_buildOverlapAndHamiltonianMatrix

  
  !>
  !! @brief Merges the occupied orbitals coefficients from two systems
  !! @param occupationI and occupationII: Number of orbitals to merge from each matrix. 
  !! sysBasisList: array indicating which basis functions of the merged molecular system belong to sysI and sysII Merged Coefficients: Matrices for output.
  !<
  subroutine NonOrthogonalCI_mergeCoefficients(coefficientsI,coefficientsII,molecularSystemI,molecularSystemII,mergedMolecularSystem,&
       sysIbasisList,sysIIbasisList,mergedCoefficients)
    type(Matrix), intent(in) :: coefficientsI(*), coefficientsII(*)
    type(MolecularSystem), intent(in) :: molecularSystemI, molecularSystemII, mergedMolecularSystem
    type(IVector), intent(in) :: sysIbasisList(*), sysIIbasisList(*)
    type(Matrix), intent(out) :: mergedCoefficients(*)
    
    ! character(100) :: wfnFile
    ! character(50) :: arguments(2)
    ! integer :: wfnUnit
    integer :: speciesID, i, j, k, l, m, mu, nu, notCommonBasis
    type(Matrix) :: auxMatrix
    type(Vector) :: auxVector
    
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

       ! call Matrix_writeToFile(mergedCoefficients(speciesID), unit=wfnUnit, binary=.true., arguments = arguments )
       
       ! arguments(1) = "DENSITY"
       ! call Matrix_constructor(auxMatrix, int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8), &
       !      int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8), 0.0_8 )

       ! auxMatrix%values=1.0
       
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

    deallocate(displacementVector)
    
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
  !! @param sysI and sysII: molecular system indices. Merged Molecular System: Union of objects from sysI and sysII. Merged Coefficients: Mixed molecular system coefficients. Sys basis list indicate the basis functions of each sysI and sysII in the merged molecular system. inverseOverlapMatrices: output required for two particle contributions
  !<
  subroutine NonOrthogonalCI_computeOverlapAndHCoreElements(this,sysI,sysII,mergedMolecularSystem,mergedCoefficients, &
       sysIbasisList, sysIIbasisList,inverseOverlapMatrices)

    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: mergedMolecularSystem
    integer :: sysI, sysII 
    type(Matrix) :: mergedCoefficients(*), inverseOverlapMatrices(*) 
    type(IVector) :: sysIbasisList(*), sysIIbasisList(*)

    integer :: speciesID
    integer :: a,b,bb,mu,nu    
    type(Matrix) :: auxMatrix
    type(Matrix) :: molecularOverlapMatrix
    type(Matrix), allocatable :: auxOverlapMatrix(:), auxKineticMatrix(:), auxAttractionMatrix(:), auxExternalPotMatrix(:), molecularHCoreMatrix(:)
    type(Vector) :: overlapDeterminant
    real(8) :: oneParticleEnergy
    
    allocate(auxOverlapMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxKineticMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxAttractionMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         auxExternalPotMatrix(mergedMolecularSystem%numberOfQuantumSpecies), &
         molecularHCoreMatrix(mergedMolecularSystem%numberOfQuantumSpecies))

    !!Initialize overlap
    this%configurationOverlapMatrix%values(sysI,sysII)=1.0
       
    call Vector_constructor(overlapDeterminant, mergedMolecularSystem%numberOfQuantumSpecies, 0.0_8)        
    
!!!!Overlap first
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
       !! Calculate one- particle integrals  
       call DirectIntegralManager_getOverlapIntegrals(mergedMolecularSystem,speciesID,&
            auxOverlapMatrix(speciesID))

       !!Test 

       ! print *, "auxOverlapMatrix", speciesID
       ! call Matrix_show(auxOverlapMatrix(speciesID))

       call Matrix_constructor(molecularOverlapMatrix, int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
            int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )

       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem) !sysI
          if(sysIbasisList(speciesID)%values(mu) .eq. 0 ) cycle
          do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)  !sysII
             if(sysIIbasisList(speciesID)%values(nu) .eq. 0) cycle
             do a=1, MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)) !sysI
                do b=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+1, &
                     MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))+ &
                     MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)) !sysII
                   bb=b-MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
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

       inverseOverlapMatrices(speciesID)=Matrix_inverse(molecularOverlapMatrix)

       ! print *, "inverseOverlapMatrices sysI, sysII", speciesID, sysI, sysII
       ! call Matrix_show(inverseOverlapMatrices(speciesID))

       call Matrix_getDeterminant(molecularOverlapMatrix,overlapDeterminant%values(speciesID),method="LU")             

       ! print *, "OverlapDeterminantLU speciesID, sysI, sysII", speciesID, sysI, sysII, overlapDeterminant%values(speciesID)
       
       this%configurationOverlapMatrix%values(sysI,sysII)=this%configurationOverlapMatrix%values(sysI,sysII)*overlapDeterminant%values(speciesID)

    end do

    ! print *, "total overlap", this%configurationOverlapMatrix%values(sysI,sysII)
    this%configurationOverlapMatrix%values(sysII,sysI)=this%configurationOverlapMatrix%values(sysI,sysII)
    
    !!Skip the rest of the evaluation if the overlap is smaller than the threshold
    if( CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0 .and. &
         abs(this%configurationOverlapMatrix%values(sysI,sysII)) .lt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) return

    !!Compute hcore if overlap is significant
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       call Matrix_constructor(auxKineticMatrix(speciesID),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),0.0_8)
       call Matrix_constructor(auxAttractionMatrix(speciesID),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),0.0_8)
       call Matrix_constructor(auxExternalPotMatrix(speciesID),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),&
            int(MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem),8),0.0_8)

       call DirectIntegralManager_getKineticIntegrals(mergedMolecularSystem,speciesID,auxKineticMatrix(speciesID))
       call DirectIntegralManager_getAttractionIntegrals(mergedMolecularSystem,speciesID,auxAttractionMatrix(speciesID))
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
            call DirectIntegralManager_getExternalPotentialIntegrals(mergedMolecularSystem,speciesID,auxExternalPotMatrix(speciesID))
       
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

       do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem) !sysI
          if(sysIbasisList(speciesID)%values(mu) .eq. 0) cycle
          do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem) !sysII
             if(sysIIbasisList(speciesID)%values(nu) .eq. 0) cycle
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

                   if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                        molecularHCoreMatrix(speciesID)%values(a,bb)=molecularHCoreMatrix(speciesID)%values(a,bb)+&
                        mergedCoefficients(speciesID)%values(mu,a)*mergedCoefficients(speciesID)%values(nu,b)*&
                        auxExternalPotMatrix(speciesID)%values(mu,nu)                         
                end do
             end do
          end do
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

    deallocate(auxOverlapMatrix, auxKineticMatrix, auxAttractionMatrix, auxExternalPotMatrix, molecularHCoreMatrix)
    
  end subroutine NonOrthogonalCI_computeOverlapAndHCoreElements
  !>
  !! @brief Computes the two particles contributions to the non diagonal elements of the hamiltonian matrix
  !!
  !! @param this, sysI,sysII: system indexes, inverseOverlapMatrices, mergedCoefficients are required to evaluate the elements
  !<
  subroutine NonOrthogonalCI_twoParticlesContributions(this,sysI,sysII,mergedMolecularSystem,inverseOverlapMatrices,mergedCoefficients,Libint2LocalInstance)
    implicit none
    type(NonOrthogonalCI) :: this
    integer :: sysI, sysII 
    type(MolecularSystem) :: mergedMolecularSystem
    type(Matrix) :: inverseOverlapMatrices(*) 
    type(Matrix) :: mergedCoefficients(*) 
    type(Libint2Interface) :: Libint2LocalInstance(*)
   
    type(matrix), allocatable :: fourCenterIntegrals(:,:)
    type(imatrix), allocatable :: twoIndexArray(:),fourIndexArray(:)
    integer :: ssize1, auxIndex, auxIndex1
    integer :: a,b,bb,c,d,dd,i,j
    real(8) :: interactionEnergy

    allocate(fourCenterIntegrals(mergedMolecularSystem%numberOfQuantumSpecies,mergedMolecularSystem%numberOfQuantumSpecies), &
         twoIndexArray(mergedMolecularSystem%numberOfQuantumSpecies), &
         fourIndexArray(mergedMolecularSystem%numberOfQuantumSpecies))

    !!Fill indexes arrays
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies       
       ! print *, "reading integrals species", i
       !!Two particle integrals indexes
       call Matrix_constructorInteger(twoIndexArray(i),  &
            int(max(MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem),MolecularSystem_getOcupationNumber(i,mergedMolecularSystem)),8), &
            int(max(MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem),MolecularSystem_getOcupationNumber(i,mergedMolecularSystem)),8) , 0 )

       c = 0
       do a=1,max(MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem),MolecularSystem_getOcupationNumber(i,mergedMolecularSystem))
          do b=a, max(MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem),MolecularSystem_getOcupationNumber(i,mergedMolecularSystem))
             c = c + 1
             twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
             twoIndexArray(i)%values(b,a) = twoIndexArray(i)%values(a,b)
          end do
       end do

       ssize1 = max(MolecularSystem_getTotalNumberOfContractions(i,mergedMolecularSystem),MolecularSystem_getOcupationNumber(i,mergedMolecularSystem))
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
    call NonOrthogonalCI_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, &
         twoIndexArray, fourIndexArray, fourCenterIntegrals, Libint2LocalInstance)

!!!Add charges
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies
       fourCenterIntegrals(i,i)%values = &
            fourCenterIntegrals(i,i)%values * mergedMolecularSystem%species(i)%charge**2.0

       do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
          fourCenterIntegrals(i,j)%values = &
               fourCenterIntegrals(i,j)%values * mergedMolecularSystem%species(i)%charge * mergedMolecularSystem%species(j)%charge
       end do
    end do
    
!!!Compute Hamiltonian Matrix element between displaced geometries

    ! !!Point charge-Point charge repulsion
    ! !!One Particle Terms
    ! !!Have already been computed     
    !!Same species repulsion
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies
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
    do i=1, mergedMolecularSystem%numberOfQuantumSpecies-1
       do j=i+1, mergedMolecularSystem%numberOfQuantumSpecies
          interactionEnergy=0.0
          ssize1 = max(MolecularSystem_getTotalNumberOfContractions(j,mergedMolecularSystem),MolecularSystem_getOcupationNumber(j,mergedMolecularSystem))
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
    this%configurationHamiltonianMatrix%values = &
         matmul( matmul( transpose( transformationMatrix%values ) , &
         this%configurationHamiltonianMatrix%values), transformationMatrix%values )

    ! print *,"transformed Hamiltonian Matrix "
    ! call Matrix_show( this%configurationHamiltonianMatrix )

    print *, "Diagonalizing non orthogonal CI Hamiltonian Matrix..."
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
    do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
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

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI matrix diagonalization : ", omp_get_wtime() - timeA ," (s)"
    
  end subroutine NonOrthogonalCI_diagonalizeCImatrix

  subroutine NonOrthogonalCI_generateDensities(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(MolecularSystem) :: auxMolecularSystem
    type(Matrix), allocatable :: auxCoefficients(:),mergedCoefficients(:), overlapMatrix(:)
    type(IVector), allocatable :: sysBasisList(:,:),auxBasisList(:)
    type(Matrix), allocatable :: mergedDensityMatrix(:,:)

    type(Matrix) :: auxMatrix, densityEigenVectors, auxdensityEigenVectors
    type(Vector) :: densityEigenValues, auxdensityEigenValues
    
    type(Matrix) :: plotPoints, molecularOverlapMatrix 
    type(Matrix), allocatable :: inverseOverlapMatrix(:)
    real(8) :: overlapDeterminant
    type(Matrix), allocatable :: orbitalsInGrid(:,:)
    type(Vector), allocatable :: densityInGrid(:)
    integer :: gridSize, state
    real(8) :: densityIntegral, auxValue
    integer :: a,b,g,i,ii,j,jj,k,p,mu,nu, sysI, sysII, speciesID

    integer :: densUnit
    character(100) :: densFile
    character(50) :: arguments(2), auxString
    character(100) :: molFileI, outFile
    integer :: outUnit
    real(8) :: timeA
    logical :: existFile
    
    !$  timeA = omp_get_wtime()

    if(CONTROL_instance%CI_STATES_TO_PRINT .eq. 0) return

    allocate(mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         auxCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         sysBasisList(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies),&
         auxBasisList(molecularSystem_instance%numberOfQuantumSpecies),&
         overlapMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         mergedDensityMatrix(CONTROL_instance%CI_STATES_TO_PRINT,molecularSystem_instance%numberOfQuantumSpecies))
    
    !Create a super molecular system
    !!!Merge coefficients from system 1 and system 2
    call MolecularSystem_mergeTwoSystems(molecularSystem_instance, this%MolecularSystems(1), this%MolecularSystems(2), &
         sysBasisList(1,:),sysBasisList(2,:))

    call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(1,:),this%HFCoefficients(2,:),&
         this%MolecularSystems(1),this%MolecularSystems(2),molecularSystem_instance,&
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
            auxMolecularSystem,this%MolecularSystems(sysI),molecularSystem_instance,&
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
    call MolecularSystem_showCartesianMatrix(molecularSystem_instance)
    
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

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time creating supermolecular system : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
    
    print *, ""
    print *, "Computing overlap integrals in for the superposed systems..."
    print *, ""
    !!Compute overlap integrals 
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
       call DirectIntegralManager_getOverlapIntegrals(molecularSystem_instance,speciesID,overlapMatrix(speciesID))
    end do
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for supermolecular overlap : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
    
    print *, ""
    print *, "Building merged density matrices for the superposed systems..."
    print *, ""
    allocate(InverseOverlapMatrix(molecularSystem_instance%numberOfQuantumSpecies))
    !!Build the merged density matrix
    do state=1, CONTROL_instance%CI_STATES_TO_PRINT
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          call Matrix_constructor(mergedDensityMatrix(state,speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
               int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8)
       end do
    end do
    !!Fill the merged density matrix
    ! "Diagonal" terms - same system
    do sysI=1, this%numberOfDisplacedSystems
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
          do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
             if(sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
             do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                if(sysBasisList(sysI,speciesID)%values(nu) .eq. 0) cycle
                do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                   ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                   do state=1, CONTROL_instance%CI_STATES_TO_PRINT
                      mergedDensityMatrix(state,speciesID)%values(mu,nu) =  mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                           this%configurationCoefficients%values(sysI,state)**2*&
                           mergedCoefficients(speciesID)%values(mu,ii)*&
                           mergedCoefficients(speciesID)%values(nu,ii)
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
          do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
             call Matrix_constructor(molecularOverlapMatrix, &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )
             call Matrix_constructor(inverseOverlapMatrix(speciesID), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI)),8), &
                  int(MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII)),8), 0.0_8 )
                
             do mu=1, MolecularSystem_getTotalNumberOfContractions(speciesID) !sysI
                if(sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu=1, MolecularSystem_getTotalNumberOfContractions(speciesID)  !sysII
                   if(sysBasisList(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                      ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                      do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))
                         jj=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))*(sysII-1)+j
                         ! print *, "i, j, mu, nu, coefI, coefII", i,j,mu,nu,mergedCoefficients(speciesID)%values(mu,ii),mergedCoefficients(speciesID)%values(nu,jj)
                         molecularOverlapMatrix%values(i,j)=molecularOverlapMatrix%values(i,j)+&
                              mergedCoefficients(speciesID)%values(mu,ii)*&
                              mergedCoefficients(speciesID)%values(nu,jj)*&
                              overlapMatrix(speciesID)%values(mu,nu)
                      end do
                   end do
                end do
             end do
             ! print *, "molecularOverlapMatrix sysI, sysII, speciesID", sysI, sysII, speciesID
             ! call Matrix_show(molecularOverlapMatrix)
             inverseOverlapMatrix(speciesID)=Matrix_inverse(molecularOverlapMatrix)
          end do
          
          ! Compute density contributions
          do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
             do mu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                if(sysBasisList(sysI,speciesID)%values(mu) .eq. 0) cycle
                do nu = 1 , MolecularSystem_getTotalNumberOfContractions(speciesID)
                   if(sysBasisList(sysII,speciesID)%values(nu) .eq. 0) cycle
                   do state=1, CONTROL_instance%CI_STATES_TO_PRINT
                      do i = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))
                         ii=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysI))*(sysI-1)+i
                         do j = 1 , MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))
                            jj=MolecularSystem_getOcupationNumber(speciesID,this%MolecularSystems(sysII))*(sysII-1)+j
                            mergedDensityMatrix(state,speciesID)%values(mu,nu) =  mergedDensityMatrix(state,speciesID)%values(mu,nu) + &
                                 this%configurationCoefficients%values(sysI,state)*&
                                 this%configurationCoefficients%values(sysII,state)*&
                                 this%configurationOverlapMatrix%values(sysI,sysII)*&
                                 inverseOverlapMatrix(speciesID)%values(j,i)*&
                                 mergedCoefficients(speciesID)%values(mu,ii)*&
                                 mergedCoefficients(speciesID)%values(nu,jj)
                            mergedDensityMatrix(state,speciesID)%values(nu,mu) = mergedDensityMatrix(state,speciesID)%values(nu,mu) + &
                                 this%configurationCoefficients%values(sysI,state)*&
                                 this%configurationCoefficients%values(sysII,state)*&
                                 this%configurationOverlapMatrix%values(sysI,sysII)*&
                                 inverseOverlapMatrix(speciesID)%values(j,i)*&
                                 mergedCoefficients(speciesID)%values(mu,ii)*&
                                 mergedCoefficients(speciesID)%values(nu,jj)
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
          ! print *, "mergedDensityMatrix", state, trim( MolecularSystem_instance%species(speciesID)%name )
          ! call Matrix_show(mergedDensityMatrix(state,speciesID))
       end do
    end do

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for merging density matrices : ", omp_get_wtime() - timeA ," (s)"
    !$  timeA = omp_get_wtime()
    
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
     
    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI natural orbitals : ", omp_get_wtime() - timeA ," (s)"

    deallocate(mergedCoefficients,&
         auxCoefficients,&
         sysBasisList,&
         auxBasisList,&
         overlapMatrix,&
         inverseOverlapMatrix,&
         mergedDensityMatrix)
      
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

  end subroutine NonOrthogonalCI_generateDensities

  !>
  !! @brief Calculate and Transform the four center integrals in one sweep without writing anything to disk
  !!
  !! @param molecularSystem, HFCoefficients: species array with the atomic coefficients, fourCenterIntegrals: species*species array to save integrals
  !<
  subroutine NonOrthogonalCI_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, twoIndexArray, fourIndexArray, fourCenterIntegrals, Libint2LocalInstance)
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
    
  end subroutine NonOrthogonalCI_transformIntegralsMemory

  
end module NonOrthogonalCI_

