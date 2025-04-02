!!******************************************************************************
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

module Configuration_
  use Vector_
  use Exception_
  use MolecularSystem_
  use InputCI_
  implicit none

  !>
  !! @brief Description
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 07-24-12
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 07-24-12 </tt>:  felix ( email@server )
  !!        -# description.
  !!   - <tt> 07-09-16 </tt>: Jorge Charry ( jacharrym@unal.edu.co )
  !!        -# Add CIS, and Fix CISD.
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<
  type, public :: Configuration
     !character(20) :: name
     !logical(2) :: isInstanced
     !type(Vector) :: coefficients
     !real(8) :: auxEnergy
     !real(8) :: coefficient
     !integer(2) :: nDeterminants
     !integer :: id
     !type(IVector), allocatable :: occupations(:,:)
     !integer, allocatable :: order(:) !! 1=single, 2=double, 3=triple, etc
     !type(IVector), allocatable :: excitations(:,:,:) !! nexcitations (order), ndeterminants, occ -> vir
     !integer, allocatable :: excitations(:,:,:) !! nexcitations (order), ndeterminants, occ -> vir
     integer(2), allocatable :: occupations(:,:) !! species, occ
     !integer(1), allocatable :: orbitals(:,:) !! species, occ
  end type Configuration

  type, public :: GlobalConfiguration

     logical :: isInstanced
     type(ivector) :: numberOfOccupiedOrbitals
     type(ivector) :: numberOfCoreOrbitals
     type(ivector) :: numberOfOrbitals
     type(ivector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     type(ivector) :: excitationType
  end type GlobalConfiguration

  type(GlobalConfiguration) :: GlobalConfiguration_instance

  public :: &
        Configuration_globalConstructor, &
        Configuration_constructor, &
        Configuration_constructorB, &
        Configuration_copyConstructor, &
!        Configuration_checkMaximumCoincidence, &
        Configuration_destructor, &
        Configuration_show, &
!        Configuration_checkTwoConfigurations, &
!        Configuration_checkCoincidence, &
        Configuration_checkCoincidenceB, &
!        Configuration_setAtMaximumCoincidence, &
        Configuration_setAtMaximumCoincidence, &
        Configuration_setAtMaximumCoincidenceB, &
        Configuration_setAtMaximumCoincidenceC
  private
contains


  subroutine Configuration_globalConstructor()
    implicit none
    integer :: i,j,c
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    call Vector_constructorInteger ( GlobalConfiguration_instance%numberOfOccupiedOrbitals, numberOfSpecies, 0 )
    call Vector_constructorInteger ( GlobalConfiguration_instance%numberOfCoreOrbitals, numberOfSpecies, 0 )
    call Vector_constructorInteger ( GlobalConfiguration_instance%numberOfOrbitals, numberOfSpecies, 0 )
    call Vector_constructorInteger ( GlobalConfiguration_instance%lambda, numberOfSpecies, 0 )
    call Vector_constructorInteger ( GlobalConfiguration_instance%excitationType, numberOfSpecies, 0 )

    
    do i=1, numberOfSpecies

       GlobalConfiguration_instance%lambda%values(i) = MolecularSystem_getLambda(i)
       GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(i) = MolecularSystem_getOcupationNumber(i) * &
           GlobalConfiguration_instance%lambda%values(i)
       GlobalConfiguration_instance%numberOfCoreOrbitals%values(i) = 0
       GlobalConfiguration_instance%numberOfOrbitals%values(i) = MolecularSystem_getTotalNumberOfContractions(i) * &
           GlobalConfiguration_instance%lambda%values(i)


      if ( InputCI_Instance(i)%coreOrbitals /= 0 ) then
        GlobalConfiguration_instance%numberOfCoreOrbitals%values(i) = InputCI_Instance(i)%coreOrbitals 
      end if

      if ( InputCI_Instance(i)%activeOrbitals /= 0 ) then
        GlobalConfiguration_instance%numberOfOrbitals%values(i) = InputCI_Instance(i)%activeOrbitals * &
                                    GlobalConfiguration_instance%lambda%values(i) + &
                                    GlobalConfiguration_instance%numberOfCoreOrbitals%values(i) 

     end if

     !if ( InputCI_Instance(i)%excitationType /= 0 ) then
     !  GlobalConfiguration_instance%excitationType%values(i) = InputCI_Instance(i)%excitationType
     !end if

    end do

  end subroutine
  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_constructor(this,occupiedCode,unoccupiedCode,order,numberOfConfigurations, c)
    implicit none
    type(Configuration) :: this
    type(IVector) :: order
    type(Vector), allocatable :: occupiedCode(:)
    type(Vector), allocatable :: unoccupiedCode(:)
    integer(8) :: numberOfConfigurations, c

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

!    call Vector_constructor( this%coefficients, numberOfConfigurations, 0.0_8 ) ???????????
!    call Vector_copyConstructorInteger( this%order, order )

!    this%auxEnergy = 0.0_8

!    this%id = c

!    this%coefficient = 1.0_8


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

!    allocate ( this%order ( numberOfSpecies ) )
!    this%order = order%values

!    this%nDeterminants = 1
!    allocate ( this%occupations(numberOfSpecies,this%nDeterminants))

!    if (allocated ( this%excitations )) deallocate ( this%excitations ) 
!    allocate ( this%excitations(2,maxval(this%order),numberOfSpecies ) )

    if (allocated ( this%occupations )) deallocate ( this%occupations ) 
    allocate ( this%occupations( maxval(GlobalConfiguration_instance%numberOfOccupiedOrbitals%values), numberOfSpecies ) )
    this%occupations = 0

!    if (allocated ( this%orbitals )) deallocate ( this%orbitals ) 
!    allocate ( this%orbitals( maxval(GlobalConfiguration_instance%numberOfOrbitals%values), numberOfSpecies ) )
!    this%orbitals = 0


    !allocate ( this%excitations(numberOfSpecies,this%nDeterminants,2,maxval(this%order%values) ) )
    !this%excitations = 0
    do i=1, numberOfSpecies

!      if ( this%order%values(i) == 2 )  this%coefficient = this%coefficient * ( 1.0_8 / 4.0_8 )
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda
       numberOfOrbitals = GlobalConfiguration_instance%numberOfOrbitals%values(i) 

!       call Vector_constructorInteger ( this%occupations(i,this%nDeterminants), numberOfOrbitals , 0 )
!       if ( this%order%values(i) > 0 ) then
!         call Vector_constructorInteger ( this%excitations(i,this%nDeterminants,1), this%order%values(i), 0)  !! nexcitations (order), occ -> vir
!         call Vector_constructorInteger ( this%excitations(i,this%nDeterminants,2), this%order%values(i), 0)  !! nexcitations (order), occ -> vir

!       end if

       do j=1, numberOfOccupiedOrbitals
          this%occupations(j,i)=j
!          this%orbitals(j,i) = 1
       end do

       !div1= int(occupiedCode%values(i))
       !div2= int(unoccupiedCode%values(i))

       do j= int(order%values(i)), 1, -1 

          div1= int(occupiedCode(i)%values(j))
          div2= int(unoccupiedCode(i)%values(j))

          !this%excitations(i,this%nDeterminants,1)%values(j) = mod( div1, 1024)
          !this%excitations(i,this%nDeterminants,2)%values(j) = mod( div2, 1024)

          !this%excitations(1,j,i) = div1
          !this%excitations(2,j,i) = div2

          this%occupations(div1,i)=div2
!          this%orbitals(div1,i)=0
!          this%orbitals(div2,i)=1

!          this%excitations(i,this%nDeterminants,1,j) = div1
!          this%excitations(i,this%nDeterminants,2,j) = div2



!          this%occupations(i,this%nDeterminants)%values( MOD ( div1, 1024 ) ) = &
!                                                       this%occupations(i,this%nDeterminants)%values( MOD ( div1, 1024 ) ) - 1
!          div1= div1/1024
!          this%occupations(i,this%nDeterminants)%values( MOD ( div2, 1024 ) ) = & 
!                                                       this%occupations(i,this%nDeterminants)%values( MOD ( div2, 1024 ) ) + 1
!          div2= div2/1024

!          this%occupations(i,this%nDeterminants)%values( div1 ) = &
!                                                       this%occupations(i,this%nDeterminants)%values( div1 ) - 1
!          this%occupations(i,this%nDeterminants)%values( div2  ) = & 
!                                                       this%occupations(i,this%nDeterminants)%values( div2 ) + 1

       end do

       !call Configuration_reverseSortElements(this%occupations(:,i),this%factor)
       !!if ( this%order%values(i) > 0 ) then
       !!  print *, "excitation"
       !!  call vector_show ( this%excitations(i,1,1))
       !!  call vector_show ( this%excitations(i,1,2))
       !!  call vector_show ( this%excitations(i,2,1))
       !!  call vector_show ( this%excitations(i,2,2))

       !!end if
       
    end do

    !this%isInstanced = .true.

  end subroutine Configuration_constructor

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_constructorB(this,orbitals,occupiedCode,unoccupiedCode,i,k,order)
    implicit none
    type(imatrix) :: this
    type(imatrix1) :: orbitals
    type(IVector) :: order
    type(Vector), allocatable :: occupiedCode(:)
    type(Vector), allocatable :: unoccupiedCode(:)

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j,jj
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital
    integer(8) :: k

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !spin orbitals not spatial orbitals
    lambda=MolecularSystem_getLambda(i)
    numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda
    !numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda - 7
    numberOfOrbitals = GlobalConfiguration_instance%numberOfOrbitals%values(i) 

    do j=1, numberOfOccupiedOrbitals
      !this%values(j,k) = 1_1
      this%values(j,k)=j
      orbitals%values(j,k) = 1
    end do

    do j= int(order%values(i)), 1, -1 

       div1= int(occupiedCode(i)%values(j))
       div2= int(unoccupiedCode(i)%values(j))

       this%values(div1,k) = div2

       orbitals%values(div1,k) = 0_1
       orbitals%values(div2,k) = 1_1

       !this%values(div1,k) = 0_1
       !this%values(div2,k) = 1_1

    end do

    jj = 0 
    do j=1, numberOfOrbitals
       if ( orbitals%values(j,k) == 1 ) then
         jj = jj + 1
         this%values(jj,k) = j  
       end if
    end do

  end subroutine Configuration_constructorB


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_copyConstructor(this,otherThis)
    implicit none
    type(Configuration) :: this, otherThis

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !call Vector_copyConstructorInteger( otherThis%order, this%order )
!    if (allocated( otherThis%order )) deallocate ( otherThis%order )
!    allocate ( otherThis%order(numberOfSpecies) )
!    otherThis%order = this%order

!    otherThis%nDeterminants = this%nDeterminants

!    otherThis%Coefficient = this%coefficient

!    if (allocated( otherThis%excitations )) deallocate ( otherThis%excitations )
!    allocate ( otherThis%excitations(2, maxval(this%order), numberOfSpecies ) )
!    allocate ( otherThis%excitations(numberOfSpecies,otherThis%nDeterminants, 2, maxval(this%order%values)) )

    if (allocated ( otherThis%occupations )) deallocate ( otherThis%occupations ) 
    allocate ( otherThis%occupations( maxval(GlobalConfiguration_instance%numberOfOccupiedOrbitals%values), numberOfSpecies ) )

    do i=1, numberOfSpecies
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
      ! numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda
       numberOfOrbitals = GlobalConfiguration_instance%numberOfOrbitals%values(i) 

!       do j = 1, this%nDeterminants
!         call Vector_constructorInteger ( otherThis%occupations(i,j), numberOfOrbitals , 0 )

!         otherThis%occupations(i,j)%values = this%occupations(i,j)%values

!         call Vector_copyConstructorInteger ( otherThis%excitations(i,j,1), this%excitations(i,j,1) )
!         call Vector_copyConstructorInteger ( otherThis%excitations(i,j,2), this%excitations(i,j,2) )
         !otherThis%excitations(i,j,1)%values = this%excitations(i,j,1)%values !! copy the occ
         !otherThis%excitations(i,j,2)%values = this%excitations(i,j,2)%values !! copy the occ

!       end do 

    end do
    !otherThis%excitations = this%excitations !! copy the occ
    otherThis%occupations = this%occupations !! copy the occ

    !otherThis%isInstanced = .true.

  end subroutine Configuration_copyConstructor


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
!!  subroutine Configuration_checkTwoConfigurations(thisA,thisB,sameConfiguration, numberOfSpecies)
!!    implicit none
!!    type(Configuration) :: thisA, thisB, auxthisB
!!    logical :: sameConfiguration
!!    integer :: numberOfSpecies
!!    
!!    integer :: lambda
!!    integer :: numberOfOccupiedOrbitals 
!!    integer :: numberOfOrbitals, numberOfSpatialOrbitals
!!    integer :: i,j,s,ii
!!    type(IVector), allocatable :: auxOccupations (:)
!!    type(IVector), allocatable :: auxExcitations (:,:)
!!    real(8), allocatable :: spatialOrbitalA(:)
!!    real(8), allocatable :: spatialOrbitalB(:)
!!    integer :: spin, spatial 
!!
!!!!    allocate ( auxthisB%occupations(numberOfSpecies,1) )
!!
!!    sameConfiguration = .false.  
!!
!!    do s = 1, numberOfSpecies
!!
!!      lambda=MolecularSystem_getLambda(s)
!!      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
!!      numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(s)*lambda
!!      numberOfSpatialOrbitals=MolecularSystem_getTotalNumberOfContractions(s)
!!
!!      if ( allocated (  spatialOrbitalA ) ) deallocate ( spatialOrbitalA )
!!      allocate ( spatialOrbitalA ( numberOfSpatialOrbitals ) )
!!      spatialOrbitalA = 0.0_8
!!
!!      if ( allocated (  spatialOrbitalB ) ) deallocate ( spatialOrbitalB )
!!      allocate ( spatialOrbitalB ( numberOfSpatialOrbitals ) )
!!      spatialOrbitalB = 0.0_8
!!
!!
!!      !!print *, "conf A"
!!      !!call vector_show (thisA%occupations(s,1))
!!      !!print *, "conf B"
!!      !!call vector_show (thisB%occupations(s,1))
!!
!!      if ( lambda > 1 ) then
!!
!!        do i = 1, numberOfOrbitals
!!          if ( thisA%occupations(s,thisA%nDeterminants)%values(i) > 0 ) then
!!            spin = mod (i,lambda)
!!            spatial = int((i+spin)/lambda)
!!            spatialOrbitalA(spatial) = spatialOrbitalA(spatial) + 1
!!          end if
!!        end do    
!!
!!        do i = 1, numberOfOrbitals
!!          if ( thisB%occupations(s,thisB%nDeterminants)%values(i) > 0 ) then
!!            spin = mod (i,lambda)
!!            spatial = int((i+spin)/lambda)
!!            spatialOrbitalB(spatial) = spatialOrbitalB(spatial) + 1
!!          end if
!!        end do    
!!
!!      else ! 
!!
!!        sameConfiguration = .false.
!!        return
!!
!!      end if
!!
!!
!!!!      print *, "spatial A", spatialOrbitalA
!!!!      print *, "spatial B", spatialOrbitalB
!!
!!     !! if ( lambda > 1 ) then
!!
!!     !!     call Vector_constructor ( auxthisB%occupations(s,1), numberOfOrbitals , 0.0_8 )
!!    
!!     !!     do i = 1, numberOfOrbitals
!!    
!!     !!       ii = mod(i,lambda) ! 1 alpha, 0 beta
!!     !!       if ( ii == 1 ) then 
!!     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 0 
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 1 
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 0
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 1 
!!     !!               end if
!!     !!       else if ( ii == 0 ) then
!!     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 0 
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 1 
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 0
!!     !!               end if
!!     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
!!     !!                  auxthisB%occupations(s,1)%values(i) = 1 
!!     !!               end if
!!     !!        end if
!!    
!!     !!     end do
!!     !!     !print *, "conf C"
!!     !!     !call vector_show (auxthisB%occupations(s,1))
!!
!!     !!   else ! 
!!
!!     !!     sameConfiguration = .false.
!!     !!     return
!!
!!     !!   end if
!!
!!    
!!!!      end do !! 
!!        
!!
!!      !print *, "ndetAB", thisA%nDeterminants, thisB%nDeterminants 
!!
!!      if ( allocated ( auxOccupations )) deallocate ( auxOccupations )
!!      allocate ( auxOccupations(thisA%nDeterminants))
!!      if ( allocated ( auxExcitations )) deallocate ( auxExcitations )
!!      allocate ( auxExcitations(thisA%nDeterminants,2))
!!
!!
!!
!!!!      do s = 1, numberOfSpecies
!!
!!     if ( sum(abs ( thisA%occupations(s,1)%values - thisB%occupations(s,1)%values ) ) > 0 ) then
!!
!!!!        if ( sum(abs ( thisA%occupations(s,1)%values - auxthisB%occupations(s,1)%values ) ) == 0 ) then
!!        if ( sum(abs ( spatialOrbitalA - spatialOrbitalB )  ) == 0 .and. &
!!             sum(abs ( thisA%occupations(s,thisA%nDeterminants)%values - thisB%occupations(s,thisB%nDeterminants)%values ) ) > 0  ) then !! avoid equal conf...
!!!!          print *, "true eqruiv"
!!          sameConfiguration = .true.  
!!
!!          !! copy from determinant 1 of B to 2 of A
!!          do i = 1, thisA%nDeterminants
!!            call Vector_copyConstructorInteger ( auxOccupations(i), thisA%occupations(s,i) ) 
!!            call Vector_copyConstructorInteger ( auxExcitations(i,1), thisA%excitations(s,i,1) ) 
!!            call Vector_copyConstructorInteger ( auxExcitations(i,2), thisA%excitations(s,i,2) ) 
!!          end do 
!!
!!          thisA%nDeterminants = thisA%nDeterminants + 1
!!
!!          deallocate ( thisA%occupations )
!!          allocate ( thisA%occupations(numberOfSpecies,thisA%nDeterminants)) !! species!
!!          deallocate ( thisA%excitations )
!!          allocate ( thisA%excitations(numberOfSpecies,thisA%nDeterminants,2)) !! species!
!!
!!
!!          do i = 1, thisA%nDeterminants-1
!!            call Vector_copyConstructorInteger ( thisA%occupations(s,i),  auxOccupations(i) ) 
!!            call Vector_copyConstructorInteger ( thisA%excitations(s,i,1), auxExcitations(i,1) ) 
!!            call Vector_copyConstructorInteger ( thisA%excitations(s,i,2), auxExcitations(i,2) )  
!!          end do 
!!
!!          call Vector_copyConstructorInteger ( thisA%occupations(s,thisA%nDeterminants), thisB%occupations(s,1) )
!!
!!          if ( thisA%order%values(s) > 0 .and. thisB%order%values(s) > 0 ) then
!!            call Vector_copyConstructorInteger (thisA%excitations(s,thisA%nDeterminants,1), thisB%excitations(s,1,1) ) !! copy the occ
!!            call Vector_copyConstructorInteger (thisA%excitations(s,thisA%nDeterminants,2), thisB%excitations(s,1,2) ) !! copy the occ
!!          end if
!!        else 
!!          sameConfiguration = .false.  
!!        end if
!!
!!      ! if ( thisA%order%values(s) > 0 ) then
!!      !   print *, "excitation ndet 1"
!!      !   call vector_show ( thisA%excitations(s,1,1))
!!      !   call vector_show ( thisA%excitations(s,1,2))
!!      !   print *, "excitation ndet 2"
!!      !   call vector_show ( thisA%excitations(s,2,1))
!!      !   call vector_show ( thisA%excitations(s,2,2))
!!      ! end if
!!
!!       end if
!!
!!      end do
!!
!!
!!  end subroutine Configuration_checkTwoConfigurations
!!
!!  !>
!!  !! @brief Constructor por omision
!!  !!
!!  !! @param this
!!  !<
!!  subroutine Configuration_checkMaximumCoincidence(thisA,thisB,factor, ia, ib, numberOfSpecies, occupiedOrbitals)
!!    implicit none
!!    type(Configuration) :: thisA, thisB
!!    real(8) :: factor
!!    integer :: numberOfSpecies
!!    
!!    real(8) :: lambda
!!    integer :: numberOfOccupiedOrbitals 
!!    integer :: numberOfOrbitals 
!!    integer :: i,j,ia,ib,s,ii,jj
!!    real(8), allocatable :: diffVector(:), occupiedOrbitals(:,:)
!!    integer, allocatable :: scoreMatrix(:,:)
!!    integer :: diagonal
!!
!!    factor = 1.0_8
!!
!!    do s = 1, numberOfSpecies
!!      if (allocated(diffVector) ) deallocate (diffVector)
!!      allocate (diffVector(size(thisA%occupations(s,ia)%values)))
!!
!!      lambda=MolecularSystem_getLambda(s)
!!      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
!!
!!      if (allocated(occupiedOrbitals) ) deallocate (occupiedOrbitals)
!!      allocate (occupiedOrbitals(numberOfOccupiedOrbitals,2))
!!      occupiedOrbitals = 0.0_8
!!
!!      if (allocated(scoreMatrix) ) deallocate (scoreMatrix)
!!      allocate (scoreMatrix(numberOfOccupiedOrbitals,numberOfOccupiedOrbitals))
!!      scoreMatrix = 0
!!
!!      ii = 0
!!      do i = 1, numberOfOccupiedOrbitals
!!        ii = ii + 1 
!!        if (thisA%occupations(s,ia)%values(i) > 0 ) then
!!          occupiedOrbitals(ii,1) = i 
!!        end if
!!      end do 
!!
!!      do i = numberOfOccupiedOrbitals+1, size(thisA%occupations(s,ia)%values )
!!        if (thisA%occupations(s,ia)%values(i) > 0 ) then
!!          do j = 1, numberOfOccupiedOrbitals
!!            if ( occupiedOrbitals(j,1) == 0 .and. mod(i,2) == mod(j,2) ) then
!!              occupiedOrbitals(j,1) = i
!!              exit
!!            end if
!!          end do
!!        end if
!!      end do 
!!
!!      ii = 0
!!      do i = 1, numberOfOccupiedOrbitals
!!        ii = ii + 1 
!!        if (thisB%occupations(s,ib)%values(i) > 0 ) then
!!          occupiedOrbitals(ii,2) = i 
!!        end if
!!      end do 
!!
!!      do i = numberOfOccupiedOrbitals+1, size(thisB%occupations(s,ib)%values )
!!        if (thisB%occupations(s,ib)%values(i) > 0 ) then
!!          do j = 1, numberOfOccupiedOrbitals
!!            if ( occupiedOrbitals(j,2) == 0 .and. mod(i,2) == mod(j,2)) then
!!              occupiedOrbitals(j,2) = i
!!              exit
!!            end if
!!          end do
!!        end if
!!      end do 
!!
!!      diagonal = 0
!!      do i = 1, numberOfOccupiedOrbitals 
!!        do j = 1, numberOfOccupiedOrbitals
!!           if ( occupiedOrbitals(i,1) == occupiedOrbitals(j,2) ) then
!!             scoreMatrix(i,j) = 1
!!           else  
!!             scoreMatrix(i,j) = 0
!!           end if 
!!           if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
!!        end do 
!!      end do 
!!
!!      !print *, "scoreMatrix total", sum(scoreMatrix)
!!      !print *, "diagonal", diagonal
!!
!!      if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 1 ) factor = -1.0_8*factor
!!      if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 0 ) factor = 1.0_8*factor
!!
!!      !print *, "factor", factor
!!
!!    end do
!!
!!  end subroutine Configuration_checkMaximumCoincidence

!  function Configuration_checkCoincidence(thisA,thisB, numberOfSpecies) result ( numberOfDiffOrbitals )
!    implicit none
!    type(Configuration), intent(in) :: thisA, thisB
!    integer, intent(in) :: numberOfSpecies
!    integer :: numberOfDiffOrbitals
!    
!    integer :: numberOfOccupiedOrbitals 
!    integer :: i,j,s
!    !type(IVector), allocatable :: occupiedOrbitals(:,:) !! spescies, confA confB
!    integer, allocatable :: occupiedOrbitalsA(:),occupiedOrbitalsB(:) !! spescies, confA confB
!    !integer, allocatable :: scoreMatrix(:,:)
!    integer :: score
!
!!    allocate (occupiedOrbitals(numberOfSpecies,2))
!    numberOfDiffOrbitals = 0 
!    do s = 1, numberOfSpecies
!  
!      !numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
!      numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
!  
!      !!call Vector_constructor (occupiedOrbitals, numberOfOccupiedOrbitals,2))
!      !call Vector_constructorInteger (occupiedOrbitals(s,1), numberOfOccupiedOrbitals, 0 )
!      !call Vector_constructorInteger (occupiedOrbitals(s,2), numberOfOccupiedOrbitals, 0 )
!
!      allocate (occupiedOrbitalsA ( numberOfOccupiedOrbitals ) )
!      occupiedOrbitalsA = 0
!      allocate (occupiedOrbitalsB ( numberOfOccupiedOrbitals ) )
!      occupiedOrbitalsB = 0
!
!  
!!      allocate (scoreMatrix(numberOfOccupiedOrbitals,numberOfOccupiedOrbitals))
!!      scoreMatrix = 0
!      score = 0
!  
!      do i = 1, numberOfOccupiedOrbitals
!        occupiedOrbitalsA(i) = i
!        occupiedOrbitalsB(i) = i
!      end do 
!  
!      do j= thisA%order(s), 1, -1 
!        occupiedOrbitalsA( thisA%excitations(1,j,s) ) = thisA%excitations(2,j,s)
!      end do
!
!      do j= thisB%order(s), 1, -1 
!        occupiedOrbitalsB( thisB%excitations(1,j,s) ) = thisB%excitations(2,j,s)
!      end do
!  
!      do i = 1, numberOfOccupiedOrbitals 
!          do j = 1, numberOfOccupiedOrbitals
!             if ( occupiedOrbitalsA(i) == occupiedOrbitalsB(j) ) then
!               !scoreMatrix(i,j) = 1
!               score = score + 1
!             end if 
!          end do 
!        end do 
!  
!       numberOfDiffOrbitals = numberOfDiffOrbitals + (numberOfOccupiedOrbitals- score )
!
! !     deallocate (scoreMatrix)
!      deallocate (occupiedOrbitalsB )
!      deallocate (occupiedOrbitalsA )
!
! !     call Vector_destructorInteger (occupiedOrbitals(s,1) )
! !     call Vector_destructorInteger (occupiedOrbitals(s,2) )
!
!    end do 
! !   deallocate (occupiedOrbitals)
!
!  end function Configuration_checkCoincidence

  function Configuration_checkCoincidenceB(thisA,thisB, numberOfSpecies) result ( numberOfDiffOrbitals )
    implicit none
    !type(Configuration) :: thisA, thisB
    integer(2), intent(in) :: thisA(:,:), thisB(:,:)
    integer, intent(in) :: numberOfSpecies
    integer(2) :: numberOfDiffOrbitals
    
    integer :: numberOfOccupiedOrbitals 
    integer :: i,j,s
    !integer, allocatable :: occupiedOrbitalsA(:),occupiedOrbitalsB(:) !! spescies, confA confB
    integer :: score

    numberOfDiffOrbitals = 0 
    do s = 1, numberOfSpecies
  
      numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
  
      score = 0
  
      do i = 1, numberOfOccupiedOrbitals 
          do j = 1, numberOfOccupiedOrbitals
             if ( thisA(i,s) == thisB(j,s) ) then
               score = score + 1
             end if 
          end do 
        end do 
  
       numberOfDiffOrbitals = numberOfDiffOrbitals + (numberOfOccupiedOrbitals- score )

    end do 

  end function Configuration_checkCoincidenceB


  subroutine Configuration_setAtMaximumCoincidence(thisA,thisB, numberOfSpecies, factor)
  ! subroutine Configuration_setAtMaximumCoincidenceB(numberOfSpecies, factor)
  !function Configuration_setAtMaximumCoincidenceB(thisA, thisB, numberOfSpecies) result ( factor)
    implicit none
    !type(Configuration) :: thisA, thisB
    integer(2) :: thisA(:,:), thisB(:,:)
    integer :: factor
    integer, intent(in) :: numberOfSpecies
    
    integer :: numberOfOccupiedOrbitals 
    integer :: i,j,k,s
    !type(IVector), allocatable, intent(out) :: occupiedOrbitals(:,:) !! spescies, confA confB
    !integer, allocatable, intent(out) :: occupiedOrbitalsA(:,:),occupiedOrbitalsB(:,:) !! spescies, confA confB
    integer :: auxOcc
    integer(2) :: score, auxscore
    integer(2) :: diagonal
    logical(1) :: swap

    factor = 1
  
    do s = 1, numberOfSpecies
        numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
  
        score = 0
        auxscore = 0
        diagonal = 0

        do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals
               if ( thisA(i,s) == thisB(j,s) ) then
                 auxscore = 1
                if ( i /= j ) then
                   factor = -1*factor
!                  if ( mod(abs(i - j ),2) == 0 ) factor = factor
                end if

               else 
                 auxscore = 0
               end if 
               if ( i == j ) diagonal = diagonal + auxscore
               score = score + auxscore
            end do 
         end do 
         ! print *, thisA(:,s), thisB(:,s)
        
     end do

  end subroutine Configuration_setAtMaximumCoincidence

  subroutine Configuration_setAtMaximumCoincidenceB(thisA,thisB, numberOfSpecies, factor)
  ! subroutine Configuration_setAtMaximumCoincidenceB(numberOfSpecies, factor)
  !function Configuration_setAtMaximumCoincidenceB(thisA, thisB, numberOfSpecies) result ( factor)
    implicit none
    !type(Configuration) :: thisA, thisB
    integer(2) :: thisA(:,:), thisB(:,:)
    integer :: factor
    integer, intent(in) :: numberOfSpecies
    
    integer :: numberOfOccupiedOrbitals 
    integer :: i,j,k,s, m ,n
    !type(IVector), allocatable, intent(out) :: occupiedOrbitals(:,:) !! spescies, confA confB
    !integer, allocatable, intent(out) :: occupiedOrbitalsA(:,:),occupiedOrbitalsB(:,:) !! spescies, confA confB
    integer :: auxOcc
    integer(2) :: score, auxscore
    integer(2) :: diagonal
    logical(1) :: swap

    factor = 1
  
    do s = 1, numberOfSpecies
        numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
  
        score = 0
        auxscore = 0
        diagonal = 0

        do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals
               if ( thisA(i,s) == thisB(j,s) ) then
                 auxscore = 1
               else 
                 auxscore = 0
               end if 
               if ( i == j ) diagonal = diagonal + auxscore
               score = score + auxscore
            end do 
         end do 
  
        do while ( (score) > diagonal )
            if ((score) > diagonal ) then
            swap = .false. 
            do i = 1, numberOfOccupiedOrbitals 
              do j = 1, numberOfOccupiedOrbitals
                if ( i /= j ) then
                  if ( thisA(i,s) == thisB(j,s) ) then

                    auxOcc = thisA(i,s)
                    thisA(i,s) = thisA(j,s)
                    thisA(j,s) = auxOcc
                    swap = .true.  
                  end if
                end if

                if ( swap .eqv. .true. ) exit
              end do 
              if ( swap .eqv. .true. ) exit
            end do 

          diagonal = 0
          score = 0
          do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals

               if ( thisA(i,s) == thisB(j,s) ) then
                  auxscore = 1
               else  
                 auxscore = 0
               end if 
               if ( i == j ) diagonal = diagonal + auxscore
               score = score + auxscore
            end do 
          end do 
            factor = -1 * factor 
          end if
        end do! while

     end do

  end subroutine Configuration_setAtMaximumCoincidenceB

  subroutine Configuration_setAtMaximumCoincidenceC(thisA,thisB, m, n, numberOfSpecies, factor)
  ! subroutine Configuration_setAtMaximumCoincidenceB(numberOfSpecies, factor)
  !function Configuration_setAtMaximumCoincidenceB(thisA, thisB, numberOfSpecies) result ( factor)
    implicit none
    !type(Configuration) :: thisA, thisB
    integer(2) :: thisA(m,n)
    integer(2), intent(in) :: thisB(:,:)
    integer :: factor
    integer, intent(in) :: numberOfSpecies
    
    integer :: numberOfOccupiedOrbitals 
    integer :: i,j,k,s, m ,n
    !type(IVector), allocatable, intent(out) :: occupiedOrbitals(:,:) !! spescies, confA confB
    !integer, allocatable, intent(out) :: occupiedOrbitalsA(:,:),occupiedOrbitalsB(:,:) !! spescies, confA confB
    integer :: auxOcc
    integer(2) :: score, auxscore
    integer(2) :: diagonal
    logical(1) :: swap

    factor = 1
  
    do s = 1, numberOfSpecies
        numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
  
        score = 0
        auxscore = 0
        diagonal = 0

        do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals
               if ( thisA(i,s) == thisB(j,s) ) then
                 auxscore = 1
               else 
                 auxscore = 0
               end if 
               if ( i == j ) diagonal = diagonal + auxscore
               score = score + auxscore
            end do 
         end do 
  
        do while ( (score) > diagonal )
            if ((score) > diagonal ) then
            swap = .false. 
            do i = 1, numberOfOccupiedOrbitals 
              do j = 1, numberOfOccupiedOrbitals
                if ( i /= j ) then
                  if ( thisA(i,s) == thisB(j,s) ) then

                    auxOcc = thisA(i,s)
                    thisA(i,s) = thisA(j,s)
                    thisA(j,s) = auxOcc
                    swap = .true.  
                  end if
                end if

                if ( swap .eqv. .true. ) exit
              end do 
              if ( swap .eqv. .true. ) exit
            end do 

          diagonal = 0
          score = 0
          do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals

               if ( thisA(i,s) == thisB(j,s) ) then
                  auxscore = 1
               else  
                 auxscore = 0
               end if 
               if ( i == j ) diagonal = diagonal + auxscore
               score = score + auxscore
            end do 
          end do 
            factor = -1 * factor 
          end if
        end do! while

     end do

  end subroutine Configuration_setAtMaximumCoincidenceC
!  
!  end function Configuration_setAtMaximumCoincidenceB


!  subroutine Configuration_setAtMaximumCoincidence(thisA,thisB, numberOfSpecies, occupiedOrbitals, factor)
!    implicit none
!    type(Configuration), intent(in) :: thisA, thisB
!    integer :: factor
!    integer, intent(in) :: numberOfSpecies
!    
!    integer :: numberOfOccupiedOrbitals 
!    integer :: i,j,k,s
!    type(IVector), allocatable, intent(out) :: occupiedOrbitals(:,:) !! spescies, confA confB
!    integer, allocatable :: occupiedOrbitalsA(:,:),occupiedOrbitalsB(:,:) !! spescies, confA confB
!    integer :: auxOcc
!    integer, allocatable :: scoreMatrix(:,:)
!    integer :: diagonal
!    logical(2) :: swap
!
!      factor = 1
!      if (allocated(occupiedOrbitals) ) deallocate (occupiedOrbitals)
!      allocate (occupiedOrbitals(numberOfSpecies,2))
!  
!    do s = 1, numberOfSpecies
!        numberOfOccupiedOrbitals = GlobalConfiguration_instance%numberOfOccupiedOrbitals%values(s) 
!  
!        !!call Vector_constructor (occupiedOrbitals, numberOfOccupiedOrbitals,2))
!        call Vector_constructorInteger (occupiedOrbitals(s,1), numberOfOccupiedOrbitals, 0 )
!        call Vector_constructorInteger (occupiedOrbitals(s,2), numberOfOccupiedOrbitals, 0 )
!  
!        do i = 1, numberOfOccupiedOrbitals
!          occupiedOrbitals(s,1)%values(i) = i
!          occupiedOrbitals(s,2)%values(i) = i
!        end do 
!  
!        do j= thisA%order(s), 1, -1 
!          occupiedOrbitals(s,1)%values( thisA%excitations(1,j,s) ) = thisA%excitations(2,j,s)
!        end do
!  
!        do j= thisB%order(s), 1, -1 
!          occupiedOrbitals(s,2)%values( thisB%excitations(1,j,s) ) = thisB%excitations(2,j,s)
!        end do
!
!        allocate (scoreMatrix(numberOfOccupiedOrbitals,numberOfOccupiedOrbitals))
!        scoreMatrix = 0
!
!        diagonal = 0
!        do i = 1, numberOfOccupiedOrbitals 
!            do j = 1, numberOfOccupiedOrbitals
!               if ( occupiedOrbitals(s,1)%values(i) == occupiedOrbitals(s,2)%values(j) ) then
!                 scoreMatrix(i,j) = 1
!               end if 
!               if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
!            end do 
!         end do 
!  
!        do while ( sum(scoreMatrix) > diagonal )
!
!            if (sum(scoreMatrix) > diagonal ) then
!              
!            swap = .false. 
!            do i = 1, numberOfOccupiedOrbitals 
!              do j = 1, numberOfOccupiedOrbitals
!                if ( i /= j ) then
!                  if ( scoreMatrix(i,j) == 1 ) then
!
!                    auxOcc = occupiedOrbitals(s,1)%values(i)
!                    occupiedOrbitals(s,1)%values(i) = occupiedOrbitals(s,1)%values(j)
!                    occupiedOrbitals(s,1)%values(j) = auxOcc
!                    swap = .true.  
!                  end if
!                end if
!
!                if ( swap .eqv. .true. ) exit
!              end do 
!              if ( swap .eqv. .true. ) exit
!            end do 
!
!          diagonal = 0
!
!          do i = 1, numberOfOccupiedOrbitals 
!            do j = 1, numberOfOccupiedOrbitals
!
!               if ( occupiedOrbitals(s,1)%values(i) == occupiedOrbitals(s,2)%values(j) ) then
!                 scoreMatrix(i,j) = 1
!               else  
!                 scoreMatrix(i,j) = 0
!               end if 
!               if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
!            end do 
!          end do 
! 
!            factor = -1 * factor 
!          end if
!        end do! while
!
!        deallocate (scoreMatrix)
!     end do
!
!  end subroutine Configuration_setAtMaximumCoincidence


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_destructor(this)
    implicit none
    type(Configuration) :: this
    integer :: i, j, numberOfSpecies

    !call Vector_destructor( this%coefficients )
    !call Vector_destructorInteger( this%order )
    !deallocate ( this%order )


    !numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

!    do i=1, numberOfSpecies
!        do j = 1, this%nDeterminants
!!!          call Vector_destructorInteger ( this%occupations(i,j) )
!          call Vector_destructorInteger ( this%excitations(i,j,1) )
!          call Vector_destructorInteger ( this%excitations(i,j,2) )
!        end do
!    end do

    !deallocate(this%excitations)
    deallocate(this%occupations)

!!    deallocate ( this%occupations )
!    this%isInstanced = .false.

  end subroutine Configuration_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine Configuration_show(this)
    implicit none
    type(Configuration) :: this

    integer :: i, j, numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !print *, "Configuration", this%id
    print *, "-------------"

    
    do i=1, numberOfSpecies
      print *, "For specie ", MolecularSystem_getSymbolOfSpecies( i )
!      print *, "Excitations: ", this%order(i)
!      print *, "Ndeterminants: ",this%nDeterminants
!      print *, "Occupations"
      
!!      do j = 1, this%nDeterminants
!!        call Vector_showInteger ( this%occupations(i,j) )
!!      end do 

!      if ( this%order(i) > 0 ) then

        !do j = 1, this%nDeterminants
        !  print *, "excitation ndet ", j
        !  print *, this%excitations(i,j,1,:)
        !  print *, this%excitations(i,j,2,:)
        !end do   
!          print *, this%excitations(i,1,:)
!          print *, this%excitations(i,2,:)

        do j = 1, int( MolecularSystem_instance%species(i)%ocupationNumber)
          print *, this%occupations(j,i)
        end do
      !end if

    end do

  end subroutine Configuration_show

!!  subroutine Configuration_reverseSortElements(this,m)
!!    integer :: this(:)
!!    integer(1) :: m
!!    integer i,j,n
!!    
!!    n = size(this)
!!    m = 1 
!!    do i=1,n
!!       do j=i+1,n
!!          if (this(j).lt.this(i)) then
!!             call Configuration_swapElements( this, i, j )
!!             m = -1 * m 
!!          end if
!!       end do
!!    end do
!!
!!  end subroutine Configuration_reverseSortElements
!!
!!  
!!  !>
!!  !! @brief Intercambia los elementos i y j el vector
!!  subroutine Configuration_swapElements( this, i, j )
!!    implicit none
!!    integer, intent(inout) :: this(:)
!!    integer, intent(in) :: i
!!    integer, intent(in) :: j
!!    integer :: value1
!!    integer :: value2
!!    
!!    value1 = this( i )
!!    value2 = this( j )
!!    
!!    this( i ) = value2
!!    this( j ) = value1
!!    
!!  end subroutine Configuration_swapElements



  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  !<
  function Configuration_isInstanced( this ) result( output )
    implicit  none
    type(Configuration), intent(in) :: this
    logical :: output

    !output = this%isInstanced

  end function Configuration_isInstanced

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine Configuration_exception( typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine Configuration_exception

end module Configuration_
