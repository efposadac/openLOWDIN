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
     ! integer :: numberOfUniqueSystems !sort of symmetry
     ! integer :: numberOfUniquePairs !sort of symmetry
     ! type(IVector) :: systemTypes  !sort of symmetry
     ! type(IMatrix) :: configurationPairTypes !, uniqueOverlapElements, uniqueHamiltonianElements
     ! type(MolecularSystem), allocatable :: uniqueMolecularSystems(:)
  end type NonOrthogonalCI

  type(NonOrthogonalCI), public :: NonOrthogonalCI_instance

  public :: &
       NonOrthogonalCI_constructor,&
       NonOrthogonalCI_displaceGeometries,&
       NonOrthogonalCI_readGeometries,&
       NonOrthogonalCI_runHFs,&
       NonOrthogonalCI_buildOverlapAndHamiltonianMatrix,&
       NonOrthogonalCI_diagonalizeCImatrix,&
       NonOrthogonalCI_generateSuperposedSystem,&
       NonOrthogonalCI_buildDensityMatrix,&
       NonOrthogonalCI_getNaturalOrbitals,&
       NonOrthogonalCI_saveToFile,&
       NonOrthogonalCI_computeFranckCondon

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
    else if(CONTROL_instance%READ_NOCI_GEOMETRIES) then
       this%transformationType="READ_GEOMETRIES"
       write (*,"(A)") "Reading input geometries from "//trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords file"
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
         this%configurationDFTcorrelationMatrix(molecularSystem_instance%numberOfQuantumSpecies,molecularSystem_instance%numberOfQuantumSpecies))  
    
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
       call NonOrthogonalCI_transformCoordinates(this,transformationCounter(1:this%numberOfTransformedCenters),originalMolecularSystem,displacedMolecularSystem,skip)       
       
       call MolecularSystem_showCartesianMatrix(displacedMolecularSystem,unit=coordsUnit)

       !Classify the system according to its distance matrix (symmetry) 
       ! if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY .eqv. .true.) &
       !      call NonOrthogonalCI_classifyNewSystem(this,systemType, newSystemFlag) 

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
    
  end subroutine NonOrthogonalCI_displaceGeometries

  !>
  !! @brief Read different geometries 
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_readGeometries(this)
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

             molecularSystem_instance%species(i)%particles(j)%origin = origin/AMSTRONG
             do mu = 1, molecularSystem_instance%species(i)%particles(j)%basis%length
                molecularSystem_instance%species(i)%particles(j)%basis%contraction(mu)%origin = &
                     molecularSystem_instance%species(i)%particles(j)%origin
             end do
          end do
       end do

       !! Point charges information       
       do i = 1, molecularSystem_instance%numberOfPointCharges
          read(coordsUnit,*) string, origin(1), origin(2), origin(3)
          
          molecularSystem_instance%pointCharges(i)%origin = origin/AMSTRONG
       end do
       call MolecularSystem_copyConstructor(this%molecularSystems(sysI), molecularSystem_instance)

    end do

    close(unit=coordsUnit)

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time reading coordinates : ", omp_get_wtime() - timeA ," (s)"
  end subroutine NonOrthogonalCI_readGeometries


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
  ! subroutine NonOrthogonalCI_classifyNewSystem(this, systemType, newSystemFlag) 
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

  ! end subroutine NonOrthogonalCI_classifyNewSystem

  !>
  !! @brief Run a Hartree-Fock calculation at displaced geometries and fill CI matrix diagonals 
  !!
  !! @param this -> NOCI instance
  !<
  subroutine NonOrthogonalCI_runHFs(this)
    implicit none
    type(NonOrthogonalCI) :: this

    integer :: sysI, speciesID, otherSpeciesID
    integer :: coordsUnit
    real(8) :: timeA
    character(100) :: coordsFile

    !$  timeA = omp_get_wtime()
    !!Read HF energy of the non displaced SCF calculation 
    ! print *, "HF reference energy is ", hfEnergy

    allocate(this%HFCoefficients(this%numberOfDisplacedSystems,molecularSystem_instance%numberOfQuantumSpecies))
    allocate(this%systemLabels(this%numberOfDisplacedSystems))

    call Matrix_constructor(this%configurationHamiltonianMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Matrix_constructor(this%configurationOverlapMatrix, int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
    call Vector_constructor(this%configurationCorrelationEnergies, this%numberOfDisplacedSystems, 0.0_8)
    do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
         call Matrix_constructor(this%configurationKineticMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationPuntualMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationExternalMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         call Matrix_constructor(this%configurationExchangeMatrix(speciesID), &
              int(this%numberOfDisplacedSystems,8), int(this%numberOfDisplacedSystems,8), 0.0_8)
         do otherSpeciesID=speciesID, MolecularSystem_instance%numberOfQuantumSpecies
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
    
    do sysI=1, this%numberOfDisplacedSystems
       write(this%systemLabels(sysI), '(A)') trim(this%molecularSystems(sysI)%description)

       !!Do SCF without calling lowdin-scf.x
       call MolecularSystem_copyConstructor(molecularSystem_instance, this%molecularSystems(sysI))
       CONTROL_instance%PRINT_LEVEL=0
       
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) &
            call MolecularSystem_saveToFile()
       
       if(allocated(WaveFunction_instance)) deallocate(WaveFunction_instance)
       allocate(WaveFunction_instance(molecularSystem_instance%numberOfQuantumSpecies))
       
       call MultiSCF_constructor(MultiSCF_instance,WaveFunction_instance,CONTROL_instance%ITERATION_SCHEME)
       
       call MultiSCF_buildHcore(MultiSCF_instance,WaveFunction_instance)

       call MultiSCF_getInitialGuess(MultiSCF_instance,WaveFunction_instance)

       if (CONTROL_instance%INTEGRAL_STORAGE == "MEMORY" ) then
          if(allocated(Libint2Instance)) deallocate(Libint2Instance)
          allocate(Libint2Instance(MolecularSystem_instance%numberOfQuantumSpecies))
          call DirectIntegralManager_constructor(Libint2Instance,MolecularSystem_instance)
          do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
             call DirectIntegralManager_getDirectIntraRepulsionIntegralsAll(&
                  speciesID, &
                  WaveFunction_instance(speciesID)%densityMatrix, & 
                  WaveFunction_instance(speciesID)%fourCenterIntegrals(speciesID)%values, &
                  MolecularSystem_instance,Libint2Instance(speciesID))
          end do

          do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies-1
             do otherSpeciesID=speciesID+1, MolecularSystem_instance%numberOfQuantumSpecies
                call DirectIntegralManager_getDirectInterRepulsionIntegralsAll(&
                     speciesID, otherSpeciesID, &
                     WaveFunction_instance(speciesID)%densityMatrix, & 
                     WaveFunction_instance(speciesID)%fourCenterIntegrals(otherSpeciesID)%values, &
                     MolecularSystem_instance,Libint2Instance(speciesID),Libint2Instance(otherSpeciesID))
             end do
          end do
       end if

       call MultiSCF_solveHartreeFockRoothan(MultiSCF_instance,WaveFunction_instance,Libint2Instance)

       !Save HF results
       ! call MultiSCF_saveWfn(MultiSCF_instance,WaveFunction_instance)
       call MolecularSystem_copyConstructor(this%molecularSystems(sysI),molecularSystem_instance)
       this%configurationHamiltonianMatrix%values(sysI,sysI)=MultiSCF_instance%totalEnergy

       do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
          this%HFCoefficients(sysI,speciesID) = WaveFunction_instance(speciesID)%waveFunctionCoefficients

          this%configurationKineticMatrix(speciesID)%values(sysI,sysI)=WaveFunction_instance(speciesID)%kineticEnergy
          this%configurationPuntualMatrix(speciesID)%values(sysI,sysI)=WaveFunction_instance(speciesID)%puntualInteractionEnergy
          this%configurationExternalMatrix(speciesID)%values(sysI,sysI)=WaveFunction_instance(speciesID)%externalPotentialEnergy
          this%configurationExchangeMatrix(speciesID)%values(sysI,sysI)=WaveFunction_instance(speciesID)%exchangeHFEnergy
          do otherSpeciesID = speciesID, molecularSystem_instance%numberOfQuantumSpecies
             this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,sysI)=&
                  WaveFunction_instance(speciesID)%hartreeEnergy(otherSpeciesID)
             this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,sysI)=&
                  WaveFunction_instance(speciesID)%exchangeCorrelationEnergy(otherSpeciesID)                  
          end do
       end do

       ! Compute HF energy with KS determinants
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          if(CONTROL_instance%METHOD.eq."RKS") then
             CONTROL_instance%METHOD="RHF"
          else
             CONTROL_instance%METHOD="UHF"
          end if
          
          do speciesID = 1, molecularSystem_instance%numberOfQuantumSpecies
             WaveFunction_instance(speciesID)%exchangeCorrelationEnergy=0.0_8
             WaveFunction_instance(speciesID)%exchangeCorrelationMatrix%values=0.0_8
             WaveFunction_instance(speciesID)%exactExchangeFraction=1.0_8
          end do
          call MultiSCF_obtainFinalEnergy(MultiSCF_instance,WaveFunction_instance,Libint2Instance)
          !Difference between HF and KS energies
          this%configurationCorrelationEnergies%values(sysI)=this%configurationHamiltonianMatrix%values(sysI,sysI)-MultiSCF_instance%totalEnergy

          if(CONTROL_instance%METHOD.eq."RHF") then
             CONTROL_instance%METHOD="RKS"
          else
             CONTROL_instance%METHOD="UKS"
          end if
       end if

       write (coordsUnit,'(A10,I10,A10,ES20.12,A20,ES20.12)') "Geometry ", sysI, "Energy", this%configurationHamiltonianMatrix%values(sysI,sysI), &
            "Correlation energy", this%configurationCorrelationEnergies%values(sysI)               
       call MolecularSystem_showCartesianMatrix(MolecularSystem_instance,unit=coordsUnit)

       if (this%numberOfDisplacedSystems .le. this%printMatrixThreshold) then 
          write (*,'(A10,I10,A10,ES20.12,A20,ES20.12)') "Geometry ", sysI, "Energy", this%configurationHamiltonianMatrix%values(sysI,sysI), &
               "Correlation energy", this%configurationCorrelationEnergies%values(sysI)               
          call MolecularSystem_showCartesianMatrix(MolecularSystem_instance)
          ! do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
          !    print *, "sysI", sysI, "speciesID", speciesID, "occupation number", MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(sysI))
          ! end do
       end if

       call DirectIntegralManager_destructor(Libint2Instance)
       call MultiSCF_destructor(MultiSCF_instance)
       
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

    close(coordsUnit)
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
    integer, allocatable :: sysIbatch(:), sysIIbatch(:)
    integer :: sysI,sysII,me,mySysII
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

    !Allocate objets to distribute in parallel
    nspecies=molecularSystem_instance%numberOfQuantumSpecies 
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
    
    systemI: do sysI=1, upperBound

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
                     this%molecularSystems(sysI), this%molecularSystems(mySysII),sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me))
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

             if(mySysII .eq. 0) cycle procs

             ! print *, "evaluating S and H elements for", sysI, mySysII

             ! cycle systemII
             !! Merge occupied coefficients into a single matrix 
             call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(sysI,:),this%HFCoefficients(mySysII,:),&
                  this%molecularSystems(sysI),this%molecularSystems(mySysII),mergedMolecularSystem(me),&
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
                this%configurationHamiltonianMatrix%values(sysI,mySysII)=0.0
                !$OMP ATOMIC
                overlapScreenedElements=overlapScreenedElements+1
                ! cycle systemII
             else

                if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
                   !DFT energy correction for off diagonal elements
                   call NonOrthogonalCI_getOffDiagonalDensityMatrix(this,sysI,mySysII,mergedCoefficients,mergedMolecularSystem(me),this%configurationOverlapMatrix%values(sysI,mySysII),&
                        inverseOverlapMatrices,sysIbasisList(1:nspecies,me),sysIIbasisList(1:nspecies,me))
                   this%configurationHamiltonianMatrix%values(sysI,mySysII)=this%configurationHamiltonianMatrix%values(sysI,mySysII)+&
                        this%configurationOverlapMatrix%values(sysI,mySysII)/2.0*&
                        (this%configurationCorrelationEnergies%values(sysI)+&
                        this%configurationCorrelationEnergies%values(mySysII))
                end if

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

          !In serial, symmetrize, free memory and print
          do me=1, batchSize
             mySysII=sysIIbatch(me)

             if(mySysII .eq. 0) exit systemII
             
             !Yu2020 magical empirical correction
             if(CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION .and. &
                  abs(this%configurationOverlapMatrix%values(sysI,mySysII)) .gt. CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD) then
                empiricalScaleFactor=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A*&
                     abs(this%configurationOverlapMatrix%values(sysI,mySysII))**CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B/&
                     abs(this%configurationOverlapMatrix%values(sysI,mySysII))
                this%configurationOverlapMatrix%values(sysI,mySysII)=&
                     this%configurationOverlapMatrix%values(sysI,mySysII)*empiricalScaleFactor
                this%configurationHamiltonianMatrix%values(sysI,mySysII)=&
                     this%configurationHamiltonianMatrix%values(sysI,mySysII)*empiricalScaleFactor                     
                do speciesID=1, nspecies
                   this%configurationKineticMatrix(speciesID)%values(sysI,mySysII)=&
                        this%configurationKineticMatrix(speciesID)%values(sysI,mySysII)*empiricalScaleFactor
                   this%configurationPuntualMatrix(speciesID)%values(sysI,mySysII)=&
                        this%configurationPuntualMatrix(speciesID)%values(sysI,mySysII)*empiricalScaleFactor
                   this%configurationExternalMatrix(speciesID)%values(sysI,mySysII)=&
                        this%configurationExternalMatrix(speciesID)%values(sysI,mySysII)*empiricalScaleFactor
                   this%configurationExchangeMatrix(speciesID)%values(sysI,mySysII)=&
                        this%configurationExchangeMatrix(speciesID)%values(sysI,mySysII)*empiricalScaleFactor
                   do otherSpeciesID=speciesID, nspecies
                      this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)=&
                           this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)*empiricalScaleFactor
                      this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)=&
                           this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)*empiricalScaleFactor
                   end do
                end do
             end if
             
             !Symmetrize
             this%configurationOverlapMatrix%values(mySysII,sysI)=this%configurationOverlapMatrix%values(sysI,mySysII)
             this%configurationHamiltonianMatrix%values(mySysII,sysI)=this%configurationHamiltonianMatrix%values(sysI,mySysII)

             do speciesID=1, nspecies
                this%configurationKineticMatrix(speciesID)%values(mySysII,sysI)=this%configurationKineticMatrix(speciesID)%values(sysI,mySysII)
                this%configurationPuntualMatrix(speciesID)%values(mySysII,sysI)=this%configurationPuntualMatrix(speciesID)%values(sysI,mySysII)
                this%configurationExternalMatrix(speciesID)%values(mySysII,sysI)=this%configurationExternalMatrix(speciesID)%values(sysI,mySysII)
                this%configurationExchangeMatrix(speciesID)%values(mySysII,sysI)=this%configurationExchangeMatrix(speciesID)%values(sysI,mySysII)
                do otherSpeciesID=speciesID, nspecies
                   this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(mySysII,sysI)=&
                        this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)
                   this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(mySysII,sysI)=&
                        this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)
                end do
             end do
             
             write (matrixUnit,'(I10,I10,ES20.12,ES20.12)') sysI, mySysII, &
                  this%configurationOverlapMatrix%values(sysI,mySysII), this%configurationHamiltonianMatrix%values(sysI,mySysII)               

             if (this%numberOfDisplacedSystems .le. this%printMatrixThreshold) then 
                write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, "Overlap element = ", this%configurationOverlapMatrix%values(sysI,mySysII)
                do speciesID = 1, nspecies                
                   write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( this%molecularSystems(sysI)%species(speciesID)%name ) // &
                        " Kinetic element = ", this%configurationKineticMatrix(speciesID)%values(sysI,mySysII)
                   write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( this%molecularSystems(sysI)%species(speciesID)%name ) // &
                        " Puntual element = ", this%configurationPuntualMatrix(speciesID)%values(sysI,mySysII)
                   if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                        write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( this%molecularSystems(sysI)%species(speciesID)%name ) // &
                        " External element = ", this%configurationExternalMatrix(speciesID)%values(sysI,mySysII)
                end do
                do speciesID=1, nspecies
                   write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( this%molecularSystems(sysI)%species(speciesID)%name ) // &
                        " Exchange element = ", this%configurationExchangeMatrix(speciesID)%values(sysI,mySysII)
                   do otherSpeciesID=speciesID, nspecies
                      write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( MolecularSystem_instance%species(speciesID)%name ) // &
                           "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                           " Hartree element = ", this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)
                   end do
                end do
                if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
                   do speciesID=1, nspecies
                      do otherSpeciesID=speciesID, nspecies
                         write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, trim( MolecularSystem_instance%species(speciesID)%name ) // &
                              "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                              " DFTcorrelation element = ", this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,mySysII)
                      end do
                   end do
                   write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, "Total DFT Correlation element = ", this%configurationOverlapMatrix%values(sysI,mySysII)/2.0*&
                        (this%configurationCorrelationEnergies%values(sysI)+&
                        this%configurationCorrelationEnergies%values(mySysII))
                end if
                write (*,'(I10,I10,A38,ES20.12)') sysI, mySysII, "Hamiltonian element = ", this%configurationHamiltonianMatrix%values(sysI,mySysII)
                print *, ""

             end if
             
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
    
  end subroutine NonOrthogonalCI_mergeCoefficients

  !>
  !! @brief Merges the occupied orbitals coefficients from two systems
  !! @param occupationI and occupationII: Number of orbitals to merge from each matrix. 
  !! sysBasisList: array indicating which basis functions of the merged molecular system belong to sysI and sysII Merged Coefficients: Matrices for output.
  !<
  subroutine NonOrthogonalCI_getOffDiagonalDensityMatrix(this,sysI,sysII,mergedCoefficients,mergedMolecularSystem,overlapElement,inverseOverlapMatrices,&
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


    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       print *, "Superposed DFT energies:"

       allocate(exchangeCorrelationMatrices(numberOfSpecies), &
            particlesInGrid(numberOfSpecies))
       call DensityFunctionalTheory_buildFinalGrid()
       call Matrix_constructor(dftEnergyMatrix, int(numberOfSpecies,8), &
            int(numberOfSpecies,8), 0.0_8 )
       do speciesID=1, numberOfSpecies
          numberOfContractions=MolecularSystem_getTotalNumberOfContractions(speciesID,mergedMolecularSystem)
          call Matrix_constructor(exchangeCorrelationMatrices(speciesID), int(numberOfContractions,8), &
               int(numberOfContractions,8), 0.0_8)
       end do
       call DensityFunctionalTheory_finalDFT(mergedDensityMatrix(1:numberOfSpecies), &
            exchangeCorrelationMatrices, &
            dftEnergyMatrix, &
            particlesInGrid)

       do speciesID = 1, numberOfSpecies
          write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
               " Particles in grid = ", particlesInGrid(speciesID)
       end do

       do speciesID = 1, numberOfSpecies
          do otherSpeciesID = speciesID, numberOfSpecies
             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                  "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                  " DFT Corr. energy = ", dftEnergyMatrix%values(speciesID,otherSpeciesID)
             this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(sysI,sysII)=dftEnergyMatrix%values(speciesID,otherSpeciesID)*overlapElement
          end do
       end do
    end if

    
    do speciesID=1, numberOfSpecies
       call Matrix_destructor(mergedDensityMatrix(speciesID))
    end do

    deallocate(mergedDensityMatrix)
    
  end subroutine NonOrthogonalCI_getOffDiagonalDensityMatrix


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
    
  end subroutine NonOrthogonalCI_prescreenOverlap

  !>
  !! @brief Classify the sysI and sysII pair according to their distance matrix
  !!
  !! @param sysI and sysII: molecular system indices.
  !<
  ! subroutine NonOrthogonalCI_classifyConfigurationPair(this,currentSysI,currentSysII,newPairFlag)
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

  ! end subroutine NonOrthogonalCI_classifyConfigurationPair

  
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
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

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
    this%configurationHamiltonianMatrix%values(sysI,sysII)=MolecularSystem_getPointChargesEnergy()*&
         this%configurationOverlapMatrix%values(sysI,sysII)
    ! print *, "Point charge-Point charge repulsion", MolecularSystem_getPointChargesEnergy()

    !!Compute hcore if overlap is significant
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

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
            ( 1.0_8/MolecularSystem_getMass( speciesID ) -1.0_8 / MolecularSystem_getTotalMass() )
       else
          auxKineticMatrix(speciesID)%values =  &
            auxKineticMatrix(speciesID)%values / &
            MolecularSystem_getMass( speciesID )
       end if

       !! Including charge
       auxAttractionMatrix(speciesID)%values=auxAttractionMatrix(speciesID)%values*(-MolecularSystem_getCharge(speciesID))                         
       
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
    do speciesID=1, MolecularSystem_instance%numberOfQuantumSpecies
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
    call NonOrthogonalCI_transformIntegralsMemory(mergedMolecularSystem, mergedCoefficients, &
         twoIndexArray, fourIndexArray, fourCenterIntegrals, Libint2LocalInstance)

!!!Add charges
    if ( .not. InterPotential_instance%isInstanced) then
       do i=1, mergedMolecularSystem%numberOfQuantumSpecies
          fourCenterIntegrals(i,i)%values = &
               fourCenterIntegrals(i,i)%values * mergedMolecularSystem%species(i)%charge**2.0

          do j = i+1 , MolecularSystem_instance%numberOfQuantumSpecies
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

  end subroutine NonOrthogonalCI_twoParticlesContributions

  !>
  !! @brief Solves the NOCI matrix equation
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_diagonalizeCImatrix(this)
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
       write(*,"(A38,F25.12)") " Point charges energy = ", MolecularSystem_getPointChargesEnergy()
       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
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
          write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
               " Kinetic energy = ", auxEnergy
       end do
       do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
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
          write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
               " Puntual energy = ", auxEnergy
       end do
       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
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
             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                  " External energy = ", auxEnergy
          end do
       end if
       do speciesID=1, molecularSystem_instance%numberOfQuantumSpecies
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
          write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                  " Exchange energy = ", auxEnergy

          do otherSpeciesID=speciesID, molecularSystem_instance%numberOfQuantumSpecies
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
             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                  "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
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
    
  end subroutine NonOrthogonalCI_diagonalizeCImatrix

  !>
  !! @brief Generates one molecular system combining all the displaced geometries and coefficients
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_generateSuperposedSystem(this)
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

    numberOfSpecies=molecularSystem_instance%numberOfQuantumSpecies
    
    allocate(this%sysBasisList(this%numberOfDisplacedSystems,numberOfSpecies),&
         auxCoefficients(numberOfSpecies),&         
         auxBasisList(numberOfSpecies))
    
    !Create a super molecular system
    !!!Merge coefficients from system 1 and system 2
    call MolecularSystem_mergeTwoSystems(this%mergedMolecularSystem, this%molecularSystems(1), this%molecularSystems(2), &
         this%sysBasisList(1,:),this%sysBasisList(2,:))
    
    call NonOrthogonalCI_mergeCoefficients(this%HFCoefficients(1,:),this%HFCoefficients(2,:),&
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
       call NonOrthogonalCI_mergeCoefficients(auxCoefficients,this%HFCoefficients(sysI,:),&
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
    !    write(*,*) " Merged Occupied Eigenvectors for: ", trim( MolecularSystem_instance%species(speciesID)%name )
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
    
  end subroutine NonOrthogonalCI_generateSuperposedSystem

  !>
  !! @brief Generates the NOCI density matrix in the superposed molecular system
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_buildDensityMatrix(this)
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
          ! print *, "this%mergedDensityMatrix", state, trim( MolecularSystem_instance%species(speciesID)%name )
          ! call Matrix_show(this%mergedDensityMatrix(state,speciesID))
          write(auxString,*) state
          arguments(2) = trim(MolecularSystem_instance%species(speciesID)%name)
          arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxString)) 
          call Matrix_writeToFile ( this%mergedDensityMatrix(state,speciesID), densUnit , arguments=arguments(1:2) )
       end do
    end do

    if(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL.ne."NONE" .or. &
         CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL.ne."NONE") then
       print *, "Superposed DFT energies:"

       allocate(exchangeCorrelationMatrices(numberOfSpecies), &
            particlesInGrid(numberOfSpecies))
       call DensityFunctionalTheory_buildFinalGrid()
       do state=1, CONTROL_instance%CI_STATES_TO_PRINT
          call Matrix_constructor(dftEnergyMatrix, int(numberOfSpecies,8), &
               int(numberOfSpecies,8), 0.0_8 )
          do speciesID=1, numberOfSpecies
             call Matrix_constructor(exchangeCorrelationMatrices(speciesID), int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), &
                  int(MolecularSystem_getTotalNumberOfContractions(speciesID),8), 0.0_8)
          end do
          call DensityFunctionalTheory_finalDFT(this%mergedDensityMatrix(state,1:numberOfSpecies), &
               exchangeCorrelationMatrices, &
               dftEnergyMatrix, &
               particlesInGrid)

          do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies                
             do otherSpeciesID = speciesID, MolecularSystem_instance%numberOfQuantumSpecies                
                write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
                     "/"//trim( MolecularSystem_instance%species(otherSpeciesID)%name ) // &
                     " DFT Corr. energy = ", dftEnergyMatrix%values(speciesID,otherSpeciesID)
             end do
          end do
       end do
    end if
        
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
    !       write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
    !            " Kinetic energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*kineticMatrix(speciesID)%values)
    !       write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name ) // &
    !            "/Fixed interact. energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*attractionMatrix(speciesID)%values)
    !       if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
    !            write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(speciesID)%name) // &
    !            " Ext Pot energy = ", sum(transpose(this%mergedDensityMatrix(state,speciesID)%values)*externalPotMatrix(speciesID)%values)
    !       print *, ""
    !    end do
    !    print *, ""
    ! end do
    ! deallocate(kineticMatrix,&
    !      attractionMatrix,&
    !      externalPotMatrix)  
    
  end subroutine NonOrthogonalCI_buildDensityMatrix

  !>
  !! @brief Generates the NOCI natural orbitals from the NOCI density matrix in the superposed molecular system
  !!
  !! @param this
  !<
  subroutine NonOrthogonalCI_getNaturalOrbitals(this)
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
          write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%name )
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

          write(*,"(A10,A10,A20,I5,A15,F17.12)") "number of ", trim(MolecularSystem_getNameOfSpecie( speciesID )) ," particles in state", state , &
               " density matrix: ", sum( transpose(this%mergedDensityMatrix(state,speciesID)%values)*this%mergedOverlapMatrix(speciesID)%values)
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

          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%name )
          arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 

          call Matrix_writeToFile ( densityEigenVectors, densUnit , arguments=arguments(1:2) )

          arguments(2) = trim( MolecularSystem_instance%species(speciesID)%name )
          arguments(1) = "OCCUPATIONS"//trim(adjustl(auxstring))

          call Vector_writeToFile( densityEigenValues, densUnit, arguments=arguments(1:2) )

          write(*,*) " End of natural orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(speciesID)%name )
       end do
    end do

    write(*,*) ""
    write(*,*) " END OF NATURAL ORBITALS"
    write(*,*) "=============================="
    write(*,*) ""

    close(densUnit)

    !$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for NOCI natural orbitals : ", omp_get_wtime() - timeA ," (s)"

    return

  end subroutine NonOrthogonalCI_getNaturalOrbitals
  
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


  !>
  !! @brief Save NOCI results to file
  !!
  !! @param 
  !<
  subroutine NonOrthogonalCI_saveToFile(this)
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
    
  end subroutine NonOrthogonalCI_saveToFile

  !>
  !! @brief Compute Franck-Condon factors from the current NOCI calculations and previous results read from file
  !!
  !! @param 
  !<
  subroutine NonOrthogonalCI_computeFranckCondon(this)
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

    call NonOrthogonalCI_mergeCoefficients(this%mergedCoefficients(:),auxCoefficients(:),&
         this%mergedMolecularSystem,MolecularSystem_instance,superMergedMolecularSystem,&
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
    
  end subroutine NonOrthogonalCI_computeFranckCondon

  
end module NonOrthogonalCI_

