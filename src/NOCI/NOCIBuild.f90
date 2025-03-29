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

module NOCIBuild_
  use Math_
  use MolecularSystem_
  use ParticleManager_
  use Lebedev_
  use Matrix_
  use Vector_
  use String_
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
     integer :: printMatrixThreshold
     integer, allocatable :: rotationCenterList(:,:)
     type(Matrix) :: configurationOverlapMatrix, configurationHamiltonianMatrix, configurationCoefficients
     type(Matrix), allocatable :: configurationKineticMatrix(:), configurationPuntualMatrix(:), configurationExternalMatrix(:), configurationExchangeMatrix(:)  
     type(Matrix), allocatable :: configurationHartreeMatrix(:,:), configurationDFTcorrelationMatrix(:,:)
     type(Vector) :: configurationCorrelationEnergies, statesEigenvalues
     type(IVector), allocatable :: sysBasisList(:,:)
     type(Matrix), allocatable :: HFCoefficients(:,:)
     type(Matrix), allocatable :: mergedCoefficients(:)
     type(Matrix), allocatable :: mergedOverlapMatrix(:)
     type(Matrix), allocatable :: mergedDensityMatrix(:,:)
     type(MolecularSystem), allocatable :: molecularSystems(:)
     type(MolecularSystem) :: mergedMolecularSystem
     character(50) :: transformationType
     character(15),allocatable :: systemLabels(:)
     real(8) :: refEnergy
     real(8), allocatable :: exactExchangeFraction(:)
     ! integer :: numberOfUniqueSystems !sort of symmetry
     ! integer :: numberOfUniquePairs !sort of symmetry
     ! type(IVector) :: systemTypes  !sort of symmetry
     ! type(IMatrix) :: configurationPairTypes !, uniqueOverlapElements, uniqueHamiltonianElements
     ! type(MolecularSystem), allocatable :: uniqueMolecularSystems(:)
  end type NonOrthogonalCI

  type(NonOrthogonalCI), public :: NOCI_instance

  public :: &
       NOCIBuild_constructor,&
       NOCIBuild_displaceGeometries,&
       NOCIBuild_readGeometries

  private

contains

  !>
  !! @brief Allocates memory and run HF calculations to be used in the construction of the NOCI matrix
  !!
  !! @param this
  !<
  subroutine NOCIBuild_constructor(this)
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
    ! this%numberOfUniqueSystems=0
    ! this%numberOfUniquePairs=0
    this%printMatrixThreshold=30
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
    else if(CONTROL_instance%ROTATION_AROUND_Z_STEP.ne.0) then
       print *, ""
       write (*,"(A)") "Rotating around the z axis the basis centers of the quantum particles "

       if(CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE .eq. 360 ) then
          this%numberOfIndividualTransformations=int(CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE/CONTROL_instance%ROTATION_AROUND_Z_STEP)
       else
          this%numberOfIndividualTransformations=int(CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE/CONTROL_instance%ROTATION_AROUND_Z_STEP)+1
       end if
       
       write (*,"(A,I5,A,I5,A)") "From 0 to ", CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE ," degrees in ", this%numberOfIndividualTransformations, " steps"
       print *, ""

       this%transformationType="ROTATION_AROUND_Z"
       this%numberOfTransformedCenters=1
    else if(CONTROL_instance%READ_NOCI_GEOMETRIES) then
       this%transformationType="READ_GEOMETRIES"
       write (*,"(A)") "Reading input geometries from "//trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords file"
    else
       STOP "To perform a NOCI calculation, please provide either instructions for a geometry transformation or a NOCI.coords file"
    end if

    ! call Vector_constructorInteger(this%systemTypes,this%numberOfIndividualTransformations**this%numberOfTransformedCenters,0)

    
    allocate(this%mergedDensityMatrix(CONTROL_instance%CI_STATES_TO_PRINT,molecularSystem_instance%numberOfQuantumSpecies),&
         this%mergedOverlapMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         this%mergedCoefficients(molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationKineticMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationPuntualMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationExternalMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationExchangeMatrix(molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationHartreeMatrix(molecularSystem_instance%numberOfQuantumSpecies,molecularSystem_instance%numberOfQuantumSpecies),&
         this%configurationDFTcorrelationMatrix(molecularSystem_instance%numberOfQuantumSpecies,molecularSystem_instance%numberOfQuantumSpecies),&
         this%exactExchangeFraction(molecularSystem_instance%numberOfQuantumSpecies))
    
    this%exactExchangeFraction(molecularSystem_instance%numberOfQuantumSpecies)=1.0_8
    
  end subroutine NOCIBuild_constructor
  !>
  !! @brief Generates different geometries and runs HF calculations at each 
  !!
  !! @param this
  !<
  subroutine NOCIBuild_displaceGeometries(this)
    implicit none
    type(NonOrthogonalCI) :: this

    type(MolecularSystem) :: originalMolecularSystem
    type(MolecularSystem) :: displacedMolecularSystem
    real(8) :: displacement
    character(100) :: coordsFile
    integer, allocatable :: transformationCounter(:)
    integer :: coordsUnit
    integer :: i,j
    integer :: closestSystem
    logical :: skip
    real(8) :: timeA
    
    !$  timeA = omp_get_wtime()
    
    call MolecularSystem_copyConstructor(originalMolecularSystem, molecularSystem_instance)

    !!Dynamically allocated through the displacement routine
    allocate(this%molecularSystems(0))
    
    allocate(transformationCounter(this%numberOfTransformedCenters))

    transformationCounter(1:this%numberOfTransformedCenters)=1
    transformationCounter(1)=0

    this%numberOfDisplacedSystems=0

    coordsUnit=333
    coordsFile=trim(CONTROL_instance%INPUT_FILE)//"trial.coords"

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
       call NOCIBuild_transformCoordinates(this,transformationCounter(1:this%numberOfTransformedCenters),originalMolecularSystem,displacedMolecularSystem,skip)       
       
       call MolecularSystem_showCartesianMatrix(displacedMolecularSystem,unit=coordsUnit)

       !Classify the system according to its distance matrix (symmetry) 
       ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
       !      call NOCIBuild_classifyNewSystem(this,systemType, newSystemFlag) 

       !Check if the new system is not beyond the max displacement
       if(skip) then
          write (coordsUnit,"(A)") "Skipping system beyond the ellipsoids boundaries"
          this%numberOfEllipsoidRejectedSystems=this%numberOfEllipsoidRejectedSystems+1                      
          cycle
       end if

       !Check if the separation between particles of the same charge is not too small
       call NOCIBuild_checkSameChargesDistance(displacedMolecularSystem,displacement,skip)

       if(skip) then
          write (coordsUnit,"(A,F20.12)") "Skipping system with same charge particle separation", displacement
          this%numberOfPPdistanceRejectedSystems=this%numberOfPPdistanceRejectedSystems+1                      
          cycle
       end if

       !Check if the separation between positive and negative particles is not too big
       call NOCIBuild_checkOppositeChargesDistance(displacedMolecularSystem,displacement,skip)

       if(skip) then
          write (coordsUnit,"(A,F20.12)") "Skipping system with positive and negative particle separation", displacement
          this%numberOfNPdistanceRejectedSystems=this%numberOfNPdistanceRejectedSystems+1                      
          cycle
       end if
       
       !Check if the new system is not to close to previous calculated systems - duplicate protection
       call NOCIBuild_checkNewSystemDisplacement(this,displacedMolecularSystem,closestSystem,displacement) 

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
          call NOCIBuild_saveSystem(this,displacedMolecularSystem)
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
         write (*,'(A10,I10,A,ES18.8,A,ES18.8)') "Rejected ", this%numberOfPPdistanceRejectedSystems, &
         " geometries with separation between same charge basis sets smaller than", CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE, &
         " or larger than", CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE

    if(this%numberOfNPdistanceRejectedSystems .gt. 0) &
         write (*,'(A10,I10,A,ES18.8)') "Rejected ", this%numberOfNPdistanceRejectedSystems, &
         " geometries with separation between positive and negative basis sets larger than", CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE

    if(this%numberOfEquivalentSystems .gt. 0) &
         write (*,'(A10,I10,A)') "Rejected ", this%numberOfEquivalentSystems, &
         " duplicated geometries after permutations"
    
    print *, ""
    
    ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
    !      call Matrix_constructorInteger(this%configurationPairTypes,int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0)   
    ! minEnergy=0.0    
    
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time displacing coordinates : ", omp_get_wtime() - timeA ," (s)"
    print *, ""
    
  end subroutine NOCIBuild_displaceGeometries

  !>
  !! @brief Read different geometries 
  !!
  !! @param this
  !<
  subroutine NOCIBuild_readGeometries(this)
    implicit none
    type(NonOrthogonalCI) :: this

    type(MolecularSystem) :: originalMolecularSystem
    real(8) :: origin(3)
    character(100) :: string,coordsFile
    integer :: coordsUnit
    integer :: sysI,i,ii,j,mu
    real(8) :: timeA
    logical :: readSuccess

    !$  timeA = omp_get_wtime()

    call MolecularSystem_copyConstructor(originalMolecularSystem, molecularSystem_instance)

    coordsUnit=333
    coordsFile=trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords"
    readSuccess=.false.

    inquire(FILE = coordsFile, EXIST = readSuccess )
    if(.not. readSuccess) then
       print *, "Didn't find the file ", trim(coordsFile)
       STOP "Please provide one or turn the readNOCIGeometries flag off!" 
    end if

    open(unit=coordsUnit, file=trim(coordsFile), status="old", form="formatted")

    read(coordsUnit,*) string, this%numberOfDisplacedSystems
    print *, "reading ", this%numberOfDisplacedSystems, " systems"

    allocate(this%molecularSystems(this%numberOfDisplacedSystems))
    
    do sysI = 1, this%numberOfDisplacedSystems
       call MolecularSystem_copyConstructor(molecularSystem_instance, originalMolecularSystem)
       write(molecularSystem_instance%description,"(I10)") sysI       
       read(coordsUnit,*) string !skip line
       read(coordsUnit,*) string !skip line 

       !! Print quatum species information
       do i = 1, molecularSystem_instance%numberOfQuantumSpecies

          !! Copy origins in open-shell case
          if(trim(molecularSystem_instance%species(i)%name) .eq. "E-BETA" ) then
             do ii = 1, i-1
                if(trim(molecularSystem_instance%species(ii)%name) .ne. "E-ALPHA" ) cycle
                do j = 1, size(molecularSystem_instance%species(i)%particles)
                   molecularSystem_instance%species(i)%particles(j)%origin = &
                        molecularSystem_instance%species(ii)%particles(j)%origin
                   do mu = 1, molecularSystem_instance%species(i)%particles(j)%basis%length
                      molecularSystem_instance%species(i)%particles(j)%basis%contraction(mu)%origin = &
                           molecularSystem_instance%species(i)%particles(j)%origin
                   end do
                end do
             end do
             cycle !skip the rest of the read
          end if

          do j = 1, size(molecularSystem_instance%species(i)%particles)
             read(coordsUnit,*) string, origin(1), origin(2), origin(3)

             molecularSystem_instance%species(i)%particles(j)%origin = origin/ANGSTROM
             do mu = 1, molecularSystem_instance%species(i)%particles(j)%basis%length
                molecularSystem_instance%species(i)%particles(j)%basis%contraction(mu)%origin = &
                     molecularSystem_instance%species(i)%particles(j)%origin
             end do
          end do
       end do

       !! Point charges information       
       do i = 1, molecularSystem_instance%numberOfPointCharges
          read(coordsUnit,*) string, origin(1), origin(2), origin(3)
          
          molecularSystem_instance%pointCharges(i)%origin = origin/ANGSTROM
       end do
       call MolecularSystem_copyConstructor(this%molecularSystems(sysI), molecularSystem_instance)

    end do

    close(unit=coordsUnit)

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time reading coordinates : ", omp_get_wtime() - timeA ," (s)"
  end subroutine NOCIBuild_readGeometries


  !>
  !! @brief Apply the transformation (translation or rotation) given by transformationCounter to each center, based in the originalMolecularSystemPositions the result is saved in molecularSystemInstance 
  !! @param this,transformationCounter,originalMolecularSystem
  !<
  subroutine NOCIBuild_transformCoordinates(this,transformationCounter,originalMolecularSystem,displacedMolecularSystem,skip)
    type(NonOrthogonalCI) :: this
    integer :: transformationCounter(*)
    type(MolecularSystem) :: originalMolecularSystem
    type(MolecularSystem), target :: displacedMolecularSystem
    logical, intent(out) :: skip

    real(8) :: centerX, centerY, centerZ, displacedOrigin(3), distanceCheck, distanceToCenter, angle, maxAngle
    integer :: center, displacementId
    real(8),allocatable :: X(:), Y(:), Z(:), W(:)
    integer :: i,j,k,p,q,s,mu, nsteps
    character(200) :: description
    
    skip=.false.

    call MolecularSystem_copyConstructor(displacedMolecularSystem, originalMolecularSystem)

    write(displacedMolecularSystem%description, '(I10)') transformationCounter(1)
    do i=2,this%numberOfTransformedCenters
       write(description, '(A)') adjustl(adjustr(displacedMolecularSystem%description)//"-"//adjustl(String_convertIntegerToString(transformationCounter(i))))
       displacedMolecularSystem%description=trim(description)
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
    else if(trim(this%transformationType).eq."ROTATION_AROUND_Z") then

       do center=1, this%numberOfTransformedCenters
          displacementId=0

          maxAngle=CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE
          nsteps=this%numberOfIndividualTransformations

          do i=1, nsteps
             ! do j=1,CONTROL_instance%NESTED_ROTATIONAL_GRIDS
             displacementId=displacementId+1
             if(displacementId .eq. transformationCounter(center) ) then
                angle=(i-1)*CONTROL_instance%ROTATION_AROUND_Z_STEP*Math_PI/180
                do s = 1, originalMolecularSystem%numberOfQuantumSpecies
                   do p = 1, size(originalMolecularSystem%species(s)%particles)

                      centerX=originalMolecularSystem%species(s)%particles(p)%origin(1)
                      centerY=originalMolecularSystem%species(s)%particles(p)%origin(2)
                      centerZ=originalMolecularSystem%species(s)%particles(p)%origin(3)

                      ! distanceToCenter=sqrt((originalMolecularSystem%allParticles(p)%particlePtr%origin(1)-centerX)**2 &
                      !      +(originalMolecularSystem%allParticles(p)%particlePtr%origin(2)-centerY)**2)

                      displacedOrigin(1)=centerX*cos(angle) - centerY*sin(angle) 
                      displacedOrigin(2)=centerX*sin(angle) + centerY*cos(angle) 
                      displacedOrigin(3)=centerZ

                      ! distanceToCenter=distanceToCenter+&
                      !      CONTROL_instance%NESTED_GRIDS_DISPLACEMENT*(j-(CONTROL_instance%NESTED_ROTATIONAL_GRIDS+1)/2.0)

                      ! call ParticleManager_setOrigin( MolecularSystem_instance%allParticles(p)%particlePtr, displacedOrigin )
                      displacedMolecularSystem%species(s)%particles(p)%origin=displacedOrigin
                      do mu = 1, displacedMolecularSystem%species(s)%particles(p)%basis%length
                         displacedMolecularSystem%species(s)%particles(p)%basis%contraction(mu)%origin = displacedOrigin
                      end do
                   end do
                end do
             end if
             ! end do
          end do
       end do
    end if
                   
  end subroutine NOCIBuild_transformCoordinates

  !>
  !! @brief Computes the distance between the particles of latest generated molecular system with all the previous saved ones 
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NOCIBuild_checkNewSystemDisplacement(this,newMolecularSystem,closestSystem,displacement) 
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

       call MolecularSystem_GetTwoSystemsDisplacement(this%molecularSystems(sysI), newMolecularSystem, displacementVector)

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
    
  end subroutine NOCIBuild_checkNewSystemDisplacement


  !>
  !! @brief Finds the maximum of the distances between the basis set center of a particle to its closest neighbour with opposite charge
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NOCIBuild_checkOppositeChargesDistance(molSys,minNPDistance,skip) 
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
    
  end subroutine NOCIBuild_checkOppositeChargesDistance

  !>
  !! @brief Finds the maximum of the distances between the basis set center of a particle to its closest neighbour with the same charge
  !!
  !! @param this, output: closestSystem: ID of previous system closest to the new one, displacement: sum of the distances between particles 
  !<
  subroutine NOCIBuild_checkSameChargesDistance(molSys,distance,skip) 
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

  end subroutine NOCIBuild_checkSameChargesDistance
  
  !>
  !! @brief Classify the new system by comparing its distance matrix to previosly saved systems
  !!
  !! @param this, systemType: integer defining system equivalence type,  newSystemFlag: returns if the system is new or not
  !<
  ! subroutine NOCIBuild_classifyNewSystem(this, systemType, newSystemFlag) 
  !   implicit none
  !   type(NonOrthogonalCI) :: this
  !   integer :: systemType
  !   logical :: newSystemFlag
    
  !   type(MolecularSystem) :: currentMolecularSystem
  !   type(Matrix) :: currentDistanceMatrix,previousDistanceMatrix

  !   integer :: sysI, i, checkingType
  !   logical :: match
    
  !   call MolecularSystem_copyConstructor(currentMolecularSystem, molecularSystem_instance)
  !   systemType=0
  !   newSystemFlag=.true.
  !   currentDistanceMatrix=ParticleManager_getDistanceMatrix()

  !   ! print *, "Current distance matrix"
  !   ! call Matrix_show(currentDistanceMatrix)

  !   types: do checkingType=1, this%numberOfUniqueSystems
  !      ! print *, "checkingType", checkingType
  !      systems: do sysI=1, this%numberOfDisplacedSystems

  !         if(this%systemTypes%values(sysI) .eq. checkingType) then
  !            call MolecularSystem_copyConstructor(molecularSystem_instance, this%molecularSystems(sysI))

  !            previousDistanceMatrix=ParticleManager_getDistanceMatrix()

  !            ! print *, "Comparing with previous distance matrix", checkingType
  !            ! call Matrix_show(previousDistanceMatrix)          
          
  !            match=.true.
  !            do i=1, size(currentDistanceMatrix%values(:,1))
  !               if(sum(abs(currentDistanceMatrix%values(i,:) - previousDistanceMatrix%values(i,:))) .gt. &
  !                    CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE ) then
  !                  match=.false.
  !                  exit
  !               end if
  !            end do
             
  !            ! print *, "match?", match

  !            if(match) then
  !               systemType=this%systemTypes%values(sysI)
  !               newSystemFlag=.false.
  !               exit types
  !            else
  !               cycle types
  !            end if
  !         end if
  !      end do systems
  !   end do types
    
  !   ! print *, "newSystemFlag", newSystemFlag
    
  !   call MolecularSystem_copyConstructor(molecularSystem_instance, currentMolecularSystem)

  ! end subroutine NOCIBuild_classifyNewSystem


  ! >
  ! @brief Saves molecular system and wfn files for a displaced system 
  
  ! @param systemID
  ! <
  subroutine NOCIBuild_saveSystem(this, newSystem)
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
        
  end subroutine NOCIBuild_saveSystem

  !>
  !! @brief Classify the sysI and sysII pair according to their distance matrix
  !!
  !! @param sysI and sysII: molecular system indices.
  !<
  ! subroutine NOCIBuild_classifyConfigurationPair(this,currentSysI,currentSysII,newPairFlag)
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

  ! end subroutine NOCIBuild_classifyConfigurationPair
  
end module NOCIBuild_

