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
     character(20) :: name
     logical :: isInstanced
     type(Vector) :: coefficients
     real(8) :: auxEnergy
     real(8) :: coefficient
     integer :: nDeterminants
     integer :: id
     type(Vector), allocatable :: occupations(:,:)
     type(Vector) :: order !! 1=single, 2=double, 3=triple, etc
     type(Vector), allocatable :: excitations(:,:,:) !! nexcitations (order), ndeterminants, occ -> vir
  end type Configuration

  public :: &
        Configuration_constructor, &
        Configuration_copyConstructor, &
        Configuration_checkMaximumCoincidence, &
        Configuration_destructor, &
        Configuration_show, &
        Configuration_checkTwoConfigurations, &
        Configuration_setAtMaximumCoincidence
  private		
contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_constructor(this,occupiedCode,unoccupiedCode,order,numberOfConfigurations, c)
    implicit none
    type(Configuration) :: this
    type(Vector) :: order
    type(Vector) :: occupiedCode
    type(Vector) :: unoccupiedCode
    integer :: numberOfConfigurations

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j,c
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

    call Vector_constructor( this%coefficients, numberOfConfigurations, 0.0_8 )
    call Vector_copyConstructor( this%order, order )

    this%auxEnergy = 0.0_8

    this%id = c

    this%coefficient = 1.0_8


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    this%nDeterminants = 1
    allocate ( this%occupations(numberOfSpecies,this%nDeterminants))

    if (allocated ( this%excitations )) deallocate ( this%excitations ) 
    allocate ( this%excitations(numberOfSpecies,this%nDeterminants,2) )

    do i=1, numberOfSpecies

      if ( this%order%values(i) == 2 )  this%coefficient = this%coefficient * ( 1.0_8 / 4.0_8 )
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda
       numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda

       call Vector_constructor ( this%occupations(i,this%nDeterminants), numberOfOrbitals , 0.0_8 )
       !! print *, "order ", this%order%values(i), int( this%order%values(i))
       if ( this%order%values(i) > 0 ) then
         call Vector_constructor ( this%excitations(i,this%nDeterminants,1), int( this%order%values(i)), 0.0_8)  !! nexcitations (order), occ -> vir
         call Vector_constructor ( this%excitations(i,this%nDeterminants,2), int( this%order%values(i)), 0.0_8)  !! nexcitations (order), occ -> vir
!!         call Vector_constructor ( this%excitations(i,2,1), int( this%order%values(i)), 0.0_8)  !! nexcitations (order), occ -> vir
!!         call Vector_constructor ( this%excitations(i,2,2), int( this%order%values(i)), 0.0_8)  !! nexcitations (order), occ -> vir

       end if
       !print *, "occ"
       !call Vector_show(occupiedCode)
       !print *, "ucc"
       !call Vector_show(unoccupiedCode)

       do j=1, numberOfOccupiedOrbitals
          this%occupations(i,this%nDeterminants)%values(j)=1
       end do

       div1= int(occupiedCode%values(i))
       div2= int(unoccupiedCode%values(i))

       do j= int(this%order%values(i)), 1, -1 

          this%excitations(i,this%nDeterminants,1)%values(j) = mod( div1, 1024)
          this%excitations(i,this%nDeterminants,2)%values(j) = mod( div2, 1024)

          this%occupations(i,this%nDeterminants)%values( MOD ( div1, 1024 ) ) = &
                                                       this%occupations(i,this%nDeterminants)%values( MOD ( div1, 1024 ) ) - 1
          div1= div1/1024
          this%occupations(i,this%nDeterminants)%values( MOD ( div2, 1024 ) ) = & 
                                                       this%occupations(i,this%nDeterminants)%values( MOD ( div2, 1024 ) ) + 1
          div2= div2/1024
       end do

       !print *, "conf"
       !call vector_show (this%occupations(i))

       !!if ( this%order%values(i) > 0 ) then
       !!  print *, "excitation"
       !!  call vector_show ( this%excitations(i,1,1))
       !!  call vector_show ( this%excitations(i,1,2))
       !!  call vector_show ( this%excitations(i,2,1))
       !!  call vector_show ( this%excitations(i,2,2))

       !!end if
       
    end do

    this%isInstanced = .true.

  end subroutine Configuration_constructor

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_copyConstructor(this,otherThis)
    implicit none
    type(Configuration) :: this, otherThis
    type(Vector) :: order

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

    call Vector_copyConstructor( otherThis%order, this%order )

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    otherThis%nDeterminants = this%nDeterminants

    otherThis%Coefficient = this%coefficient

    if (allocated( otherThis%occupations )) deallocate ( otherThis%occupations ) !...
    if (allocated( otherThis%excitations )) deallocate ( otherThis%excitations )

    allocate ( otherThis%occupations(numberOfSpecies,otherThis%nDeterminants ) )
    allocate ( otherThis%excitations(numberOfSpecies,otherThis%nDeterminants, 2) )

    do i=1, numberOfSpecies
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda

       do j = 1, this%nDeterminants
         call Vector_constructor ( otherThis%occupations(i,j), numberOfOrbitals , 0.0_8 )

         otherThis%occupations(i,j)%values = this%occupations(i,j)%values

         call Vector_copyConstructor ( otherThis%excitations(i,j,1), this%excitations(i,j,1) )
         call Vector_copyConstructor ( otherThis%excitations(i,j,2), this%excitations(i,j,2) )
         !otherThis%excitations(i,j,1)%values = this%excitations(i,j,1)%values !! copy the occ
         !otherThis%excitations(i,j,2)%values = this%excitations(i,j,2)%values !! copy the occ
       end do 

    end do

    otherThis%isInstanced = .true.

  end subroutine Configuration_copyConstructor


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_checkTwoConfigurations(thisA,thisB,sameConfiguration, numberOfSpecies)
    implicit none
    type(Configuration) :: thisA, thisB, auxthisB
    logical :: sameConfiguration
    integer :: numberOfSpecies
    
    integer :: lambda
    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals, numberOfSpatialOrbitals
    integer :: i,j,s,ii
    type(Vector), allocatable :: auxOccupations (:)
    type(Vector), allocatable :: auxExcitations (:,:)
    real(8), allocatable :: spatialOrbitalA(:)
    real(8), allocatable :: spatialOrbitalB(:)
    integer :: spin, spatial 

    allocate ( auxthisB%occupations(numberOfSpecies,1) )

    sameConfiguration = .false.  

    do s = 1, numberOfSpecies

      lambda=MolecularSystem_getLambda(s)
      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
      numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(s)*lambda
      numberOfSpatialOrbitals=MolecularSystem_getTotalNumberOfContractions(s)

      if ( allocated (  spatialOrbitalA ) ) deallocate ( spatialOrbitalA )
      allocate ( spatialOrbitalA ( numberOfSpatialOrbitals ) )
      spatialOrbitalA = 0.0_8

      if ( allocated (  spatialOrbitalB ) ) deallocate ( spatialOrbitalB )
      allocate ( spatialOrbitalB ( numberOfSpatialOrbitals ) )
      spatialOrbitalB = 0.0_8


      !!print *, "conf A"
      !!call vector_show (thisA%occupations(s,1))
      !!print *, "conf B"
      !!call vector_show (thisB%occupations(s,1))

      if ( lambda > 1 ) then

        do i = 1, numberOfOrbitals
          if ( thisA%occupations(s,thisA%nDeterminants)%values(i) > 0.0_8 ) then
            spin = mod (i,lambda)
            spatial = int((i+spin)/lambda)
            spatialOrbitalA(spatial) = spatialOrbitalA(spatial) + 1
          end if
        end do    

        do i = 1, numberOfOrbitals
          if ( thisB%occupations(s,thisB%nDeterminants)%values(i) > 0.0_8 ) then
            spin = mod (i,lambda)
            spatial = int((i+spin)/lambda)
            spatialOrbitalB(spatial) = spatialOrbitalB(spatial) + 1
          end if
        end do    

      else ! 

        sameConfiguration = .false.
        return

      end if


!!      print *, "spatial A", spatialOrbitalA
!!      print *, "spatial B", spatialOrbitalB

     !! if ( lambda > 1 ) then

     !!     call Vector_constructor ( auxthisB%occupations(s,1), numberOfOrbitals , 0.0_8 )
    
     !!     do i = 1, numberOfOrbitals
    
     !!       ii = mod(i,lambda) ! 1 alpha, 0 beta
     !!       if ( ii == 1 ) then 
     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 0 
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 1 
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 0
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 1 
     !!               end if
     !!       else if ( ii == 0 ) then
     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 0 
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 1 
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 0
     !!               end if
     !!               if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
     !!                  auxthisB%occupations(s,1)%values(i) = 1 
     !!               end if
     !!        end if
    
     !!     end do
     !!     !print *, "conf C"
     !!     !call vector_show (auxthisB%occupations(s,1))

     !!   else ! 

     !!     sameConfiguration = .false.
     !!     return

     !!   end if

    
!!      end do !! 
        

      !print *, "ndetAB", thisA%nDeterminants, thisB%nDeterminants 

      if ( allocated ( auxOccupations )) deallocate ( auxOccupations )
      allocate ( auxOccupations(thisA%nDeterminants))
      if ( allocated ( auxExcitations )) deallocate ( auxExcitations )
      allocate ( auxExcitations(thisA%nDeterminants,2))



!!      do s = 1, numberOfSpecies

     if ( sum(abs ( thisA%occupations(s,1)%values - thisB%occupations(s,1)%values ) ) > 0 ) then

!!        if ( sum(abs ( thisA%occupations(s,1)%values - auxthisB%occupations(s,1)%values ) ) == 0 ) then
        if ( sum(abs ( spatialOrbitalA - spatialOrbitalB )  ) == 0 .and. &
             sum(abs ( thisA%occupations(s,thisA%nDeterminants)%values - thisB%occupations(s,thisB%nDeterminants)%values ) ) > 0  ) then !! avoid equal conf...
!!          print *, "true eqruiv"
          sameConfiguration = .true.  

          !! copy from determinant 1 of B to 2 of A
          do i = 1, thisA%nDeterminants
            call Vector_copyConstructor ( auxOccupations(i), thisA%occupations(s,i) ) 
            call Vector_copyConstructor ( auxExcitations(i,1), thisA%excitations(s,i,1) ) 
            call Vector_copyConstructor ( auxExcitations(i,2), thisA%excitations(s,i,2) ) 
          end do 

          thisA%nDeterminants = thisA%nDeterminants + 1

          deallocate ( thisA%occupations )
          allocate ( thisA%occupations(numberOfSpecies,thisA%nDeterminants)) !! species!
          deallocate ( thisA%excitations )
          allocate ( thisA%excitations(numberOfSpecies,thisA%nDeterminants,2)) !! species!


          do i = 1, thisA%nDeterminants-1
            call Vector_copyConstructor ( thisA%occupations(s,i),  auxOccupations(i) ) 
            call Vector_copyConstructor ( thisA%excitations(s,i,1), auxExcitations(i,1) ) 
            call Vector_copyConstructor ( thisA%excitations(s,i,2), auxExcitations(i,2) )  
          end do 

          call Vector_copyConstructor ( thisA%occupations(s,thisA%nDeterminants), thisB%occupations(s,1) )

          if ( thisA%order%values(s) > 0 .and. thisB%order%values(s) > 0 ) then
            call Vector_copyConstructor (thisA%excitations(s,thisA%nDeterminants,1), thisB%excitations(s,1,1) ) !! copy the occ
            call Vector_copyConstructor (thisA%excitations(s,thisA%nDeterminants,2), thisB%excitations(s,1,2) ) !! copy the occ
          end if
        else 
          sameConfiguration = .false.  
        end if

      ! if ( thisA%order%values(s) > 0 ) then
      !   print *, "excitation ndet 1"
      !   call vector_show ( thisA%excitations(s,1,1))
      !   call vector_show ( thisA%excitations(s,1,2))
      !   print *, "excitation ndet 2"
      !   call vector_show ( thisA%excitations(s,2,1))
      !   call vector_show ( thisA%excitations(s,2,2))
      ! end if

       end if

      end do


  end subroutine Configuration_checkTwoConfigurations

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_checkMaximumCoincidence(thisA,thisB,factor, ia, ib, numberOfSpecies, occupiedOrbitals)
    implicit none
    type(Configuration) :: thisA, thisB
    real(8) :: factor
    integer :: numberOfSpecies
    
    real(8) :: lambda
    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j,ia,ib,s,ii,jj
    real(8), allocatable :: diffVector(:), occupiedOrbitals(:,:)
    integer, allocatable :: scoreMatrix(:,:)
    integer :: diagonal

    factor = 1.0_8

    do s = 1, numberOfSpecies
      if (allocated(diffVector) ) deallocate (diffVector)
      allocate (diffVector(size(thisA%occupations(s,ia)%values)))

      lambda=MolecularSystem_getLambda(s)
      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda

      if (allocated(occupiedOrbitals) ) deallocate (occupiedOrbitals)
      allocate (occupiedOrbitals(numberOfOccupiedOrbitals,2))
      occupiedOrbitals = 0.0_8

      if (allocated(scoreMatrix) ) deallocate (scoreMatrix)
      allocate (scoreMatrix(numberOfOccupiedOrbitals,numberOfOccupiedOrbitals))
      scoreMatrix = 0

      ii = 0
      do i = 1, numberOfOccupiedOrbitals
        ii = ii + 1 
        if (thisA%occupations(s,ia)%values(i) > 0 ) then
          occupiedOrbitals(ii,1) = i 
        end if
      end do 

      do i = numberOfOccupiedOrbitals+1, size(thisA%occupations(s,ia)%values )
        if (thisA%occupations(s,ia)%values(i) > 0 ) then
          do j = 1, numberOfOccupiedOrbitals
            if ( occupiedOrbitals(j,1) == 0 .and. mod(i,2) == mod(j,2) ) then
              occupiedOrbitals(j,1) = i
              exit
            end if
          end do
        end if
      end do 

      ii = 0
      do i = 1, numberOfOccupiedOrbitals
        ii = ii + 1 
        if (thisB%occupations(s,ib)%values(i) > 0 ) then
          occupiedOrbitals(ii,2) = i 
        end if
      end do 

      do i = numberOfOccupiedOrbitals+1, size(thisB%occupations(s,ib)%values )
        if (thisB%occupations(s,ib)%values(i) > 0 ) then
          do j = 1, numberOfOccupiedOrbitals
            if ( occupiedOrbitals(j,2) == 0 .and. mod(i,2) == mod(j,2)) then
              occupiedOrbitals(j,2) = i
              exit
            end if
          end do
        end if
      end do 

      diagonal = 0
      do i = 1, numberOfOccupiedOrbitals 
        do j = 1, numberOfOccupiedOrbitals
           if ( occupiedOrbitals(i,1) == occupiedOrbitals(j,2) ) then
             scoreMatrix(i,j) = 1
           else  
             scoreMatrix(i,j) = 0
           end if 
           if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
        end do 
      end do 

      !print *, "scoreMatrix total", sum(scoreMatrix)
      !print *, "diagonal", diagonal

      if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 1 ) factor = -1.0_8*factor
      if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 0 ) factor = 1.0_8*factor

      !print *, "factor", factor

    end do

  end subroutine Configuration_checkMaximumCoincidence

  subroutine Configuration_setAtMaximumCoincidence(thisA,thisB, ia, ib, numberOfSpecies, occupiedOrbitals, factor)
    implicit none
    type(Configuration) :: thisA, thisB
    real(8) :: factor
    integer :: numberOfSpecies
    
    real(8) :: lambda
    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j,k,ia,ib,s,ii,jj
    real(8), allocatable :: diffVector(:)
    type(Vector), allocatable :: occupiedOrbitals(:,:) !! spescies, confA confB
    real(8) :: auxOcc
    integer, allocatable :: scoreMatrix(:,:)
    integer :: diagonal
    logical :: swap

    factor = 1.0_8

    if (allocated(occupiedOrbitals) ) deallocate (occupiedOrbitals)
    allocate (occupiedOrbitals(numberOfSpecies,2))


    do s = 1, numberOfSpecies
      if (allocated(diffVector) ) deallocate (diffVector)
      allocate (diffVector(size(thisA%occupations(s,ia)%values)))

      lambda=MolecularSystem_getLambda(s)
      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
      !print *, "id ab", thisA%id, thisB%id

      !!
      !!call Vector_constructor (occupiedOrbitals, numberOfOccupiedOrbitals,2))
      call Vector_constructor (occupiedOrbitals(s,1), numberOfOccupiedOrbitals, 0.0_8 )
      call Vector_constructor (occupiedOrbitals(s,2), numberOfOccupiedOrbitals, 0.0_8 )

      if (allocated(scoreMatrix) ) deallocate (scoreMatrix)
      allocate (scoreMatrix(numberOfOccupiedOrbitals,numberOfOccupiedOrbitals))
      scoreMatrix = 0

      ii = 0
      do i = 1, numberOfOccupiedOrbitals
        ii = ii + 1 
        if (thisA%occupations(s,ia)%values(i) > 0 ) then
          !!occupiedOrbitals(ii,1) = i 
          occupiedOrbitals(s,1)%values(ii) = i 
        end if
      end do 

      do i = numberOfOccupiedOrbitals+1, size(thisA%occupations(s,ia)%values )
        if (thisA%occupations(s,ia)%values(i) > 0 ) then
          do j = 1, numberOfOccupiedOrbitals
            !!if ( occupiedOrbitals(j,1) == 0 .and. mod(i,2) == mod(j,2) ) then
            if ( occupiedOrbitals(s,1)%values(j) == 0 .and. mod(i,int(lambda)) == mod(j,int(lambda)) ) then
              !!occupiedOrbitals(j,1) = i
              occupiedOrbitals(s,1)%values(j) = i
              exit
            end if
          end do
        end if
      end do 

      ii = 0
      do i = 1, numberOfOccupiedOrbitals
        ii = ii + 1 
        if (thisB%occupations(s,ib)%values(i) > 0 ) then
          !!occupiedOrbitals(ii,2) = i 
          occupiedOrbitals(s,2)%values(ii) = i 
        end if
      end do 

      do i = numberOfOccupiedOrbitals+1, size(thisB%occupations(s,ib)%values )
        if (thisB%occupations(s,ib)%values(i) > 0 ) then
          do j = 1, numberOfOccupiedOrbitals
            if ( occupiedOrbitals(s,2)%values(j) == 0 .and. mod(i,int(lambda)) == mod(j,int(lambda))) then
            !!if ( occupiedOrbitals(j,2) == 0 .and. mod(i,2) == mod(j,2)) then
              !!occupiedOrbitals(j,2) = i
              occupiedOrbitals(s,2)%values(j) = i
              exit
            end if
          end do
        end if
      end do 
!      print *, "occ i1", occupiedOrbitals(s,1)%values(:)
!      call vector_show (occupiedOrbitals(s,1))

!      print *, "occ i2", occupiedOrbitals(s,2)%values(:)
 !     call vector_show (occupiedOrbitals(s,2))


      diagonal = 0
      do i = 1, numberOfOccupiedOrbitals 
          do j = 1, numberOfOccupiedOrbitals
             !!if ( occupiedOrbitals(i,1) == occupiedOrbitals(j,2) ) then
             if ( occupiedOrbitals(s,1)%values(i) == occupiedOrbitals(s,2)%values(j) ) then
               scoreMatrix(i,j) = 1
             else  
               scoreMatrix(i,j) = 0
             end if 
             if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
          end do 
        end do 
        !print *, "scoreMatrix total0", sum(scoreMatrix)
        !print *, "diagonal0", diagonal


      do while ( sum(scoreMatrix) > diagonal )
        !print *, "while"
            if (sum(scoreMatrix) > diagonal ) then
            
          swap = .false. 
          do i = 1, numberOfOccupiedOrbitals 
            do j = 1, numberOfOccupiedOrbitals
              if ( i /= j ) then
                if ( scoreMatrix(i,j) == 1 ) then

                  !auxOcc = occupiedOrbitals(i,1)
                  !occupiedOrbitals(i,1) = occupiedOrbitals(j,1)
                  !occupiedOrbitals(j,1) = auxOcc
  
                  auxOcc = occupiedOrbitals(s,1)%values(i)
                  occupiedOrbitals(s,1)%values(i) = occupiedOrbitals(s,1)%values(j)
                  occupiedOrbitals(s,1)%values(j) = auxOcc
                  swap = .true.  
                end if
              end if

              if ( swap .eqv. .true. ) exit
            end do 
            if ( swap .eqv. .true. ) exit
          end do 

          !print *, "occ 1", occupiedOrbitals(:,1) 
          !print *, "occ 2", occupiedOrbitals(:,2) 
        diagonal = 0
        do i = 1, numberOfOccupiedOrbitals 
          do j = 1, numberOfOccupiedOrbitals
             !!if ( occupiedOrbitals(i,1) == occupiedOrbitals(j,2) ) then
             if ( occupiedOrbitals(s,1)%values(i) == occupiedOrbitals(s,2)%values(j) ) then

               scoreMatrix(i,j) = 1
             else  
               scoreMatrix(i,j) = 0
             end if 
             if ( i == j ) diagonal = diagonal + scoreMatrix(i,j)
          end do 
        end do 
 
          factor = -1.0_8 * factor 
          !print *, "scoreMatrix total1", sum(scoreMatrix)
          !print *, "diagonal", diagonal
    
          !if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 1 ) factor = -1.0_8*factor
          !if (  mod( sum(scoreMatrix) - diagonal, 2 ) == 0 ) factor = 1.0_8*factor
        end if
      end do! while
      !print *, "factorF", factor
    end do

  end subroutine Configuration_setAtMaximumCoincidence

  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_destructor(this)
    implicit none
    type(Configuration) :: this
    integer :: i, numberOfSpecies

    call Vector_destructor( this%coefficients )
    call Vector_destructor( this%order )

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do i=1, numberOfSpecies
       call Vector_destructor ( this%occupations(i,1) )
      
    end do

    deallocate ( this%occupations )

    this%isInstanced = .false.

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

    print *, "Configuration", this%id
    print *, "-------------"

    
    do i=1, numberOfSpecies
      print *, "For specie ", MolecularSystem_getNameOfSpecie ( i )
      print *, "Excitations: ", this%order%values(i)
      print *, "Ndeterminants: ",this%nDeterminants
      print *, "Occupations"
      
      do j = 1, this%nDeterminants
        call Vector_show ( this%occupations(i,j) )
      end do 

      if ( this%order%values(i) > 0 ) then

        do j = 1, this%nDeterminants
          print *, "excitation ndet ", j
          call vector_show ( this%excitations(i,j,1))
          call vector_show ( this%excitations(i,j,2))
        end do   
      end if

    end do

  end subroutine Configuration_show

  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  !<
  function Configuration_isInstanced( this ) result( output )
    implicit  none
    type(Configuration), intent(in) :: this
    logical :: output

    output = this%isInstanced

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
