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
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<
  type, public :: Configuration
     character(20) :: name
     logical :: isInstanced
     type(Vector) :: coefficients
     real(8) :: auxEnergy
     integer :: nDeterminants
     type(Vector), allocatable :: occupations(:,:)
     type(Vector) :: order !! 1=single, 2=double, 3=triple, etc
  end type Configuration

  public :: &
       Configuration_constructor, &
       Configuration_copyConstructor, &
       Configuration_checkMaximumCoincidence, &
       Configuration_destructor, &
       Configuration_show, &
       Configuration_checkTwoConfigurations

  private		
contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_constructor(this,occupiedCode,unoccupiedCode,order,numberOfConfigurations)
    implicit none
    type(Configuration) :: this
    type(Vector) :: order
    type(Vector) :: occupiedCode
    type(Vector) :: unoccupiedCode
    integer :: numberOfConfigurations

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

    call Vector_constructor( this%coefficients, numberOfConfigurations, 0.0_8 )
    call Vector_copyConstructor( this%order, order )
    this%auxEnergy = 0.0_8

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    this%nDeterminants = 1
    allocate ( this%occupations(numberOfSpecies,2) )

    do i=1, numberOfSpecies
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda
       numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda

       call Vector_constructor ( this%occupations(i,1), numberOfOrbitals , 0.0_8 )
       call Vector_constructor ( this%occupations(i,2), numberOfOrbitals , 0.0_8 )
      !  print *, "occ"
       !call Vector_show(occupiedCode)
       ! print *, "ucc"
       !call Vector_show(unoccupiedCode)

       do j=1, numberOfOccupiedOrbitals
          this%occupations(i,1)%values(j)=1
       end do

       div1= int(occupiedCode%values(i))
       div2= int(unoccupiedCode%values(i))

       do j=this%order%values(i), 1, -1 
          this%occupations(i,1)%values( MOD ( div1, 1024 ) ) = this%occupations(i,1)%values( MOD ( div1, 1024 ) ) - 1
          div1= div1/1024
          this%occupations(i,1)%values( MOD ( div2, 1024 ) ) = this%occupations(i,1)%values( MOD ( div2, 1024 ) ) + 1
          div2= div2/1024
       end do

       !print *, "conf"
       !call vector_show (this%occupations(i))

       
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
    allocate ( otherThis%occupations(numberOfSpecies,2) )

    do i=1, numberOfSpecies
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda

       call Vector_constructor ( otherThis%occupations(i,1), numberOfOrbitals , 0.0_8 )
       call Vector_constructor ( otherThis%occupations(i,2), numberOfOrbitals , 0.0_8 )

       otherThis%occupations(i,1)%values = this%occupations(i,1)%values
       otherThis%occupations(i,2)%values = this%occupations(i,2)%values

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
    integer :: numberOfOrbitals 
    integer :: i,j,s,ii


    allocate ( auxthisB%occupations(numberOfSpecies,1) )

    do s = 1, numberOfSpecies

      lambda=MolecularSystem_getLambda(s)
      numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(s)*lambda
      numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(s)*lambda

      !print *, "conf A"
      !call vector_show (thisA%occupations(s,1))
      !print *, "conf B"
      !call vector_show (thisB%occupations(s,1))

      if ( lambda > 1 ) then

          call Vector_constructor ( auxthisB%occupations(s,1), numberOfOrbitals , 0.0_8 )
    
          do i = 1, numberOfOrbitals
    
            ii = mod(i,lambda) ! 1 alpha, 0 beta
            if ( ii == 1 ) then 
                    if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
                       auxthisB%occupations(s,1)%values(i) = 0 
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
                       auxthisB%occupations(s,1)%values(i) = 1 
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 0 ) then
                       auxthisB%occupations(s,1)%values(i) = 0
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i+1) == 1 ) then
                       auxthisB%occupations(s,1)%values(i) = 1 
                    end if
            else if ( ii == 0 ) then
                    if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
                       auxthisB%occupations(s,1)%values(i) = 0 
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 0 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
                       auxthisB%occupations(s,1)%values(i) = 1 
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 0 ) then
                       auxthisB%occupations(s,1)%values(i) = 0
                    end if
                    if ( thisB%occupations(s,1)%values(i) == 1 .and. thisB%occupations(s,1)%values(i-1) == 1 ) then
                       auxthisB%occupations(s,1)%values(i) = 1 
                    end if
             end if
    
          end do
          !print *, "conf C"
          !call vector_show (auxthisB%occupations(s,1))

        else ! 

          sameConfiguration = .false.
          return

        end if

    
      end do
        
      sameConfiguration = .false.  
      do s = 1, numberOfSpecies
        if ( sum(abs ( thisA%occupations(s,1)%values - auxthisB%occupations(s,1)%values ) ) == 0 ) then
          !print *, "true eqruiv"
          sameConfiguration = .true.  
          thisA%nDeterminants = 2
          thisA%occupations(s,2)%values =  thisB%occupations(s,1)%values
        else 
          sameConfiguration = .false.  
        end if
      end do


  end subroutine Configuration_checkTwoConfigurations

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Configuration_checkMaximumCoincidence(thisA,thisB,factor, ia, ib, numberOfSpecies)
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

      !call vector_show(thisA%occupations(s,ia))
      !call vector_show(thisB%occupations(s,ib))

!      ii = 0
!      do i = 1, size(thisA%occupations(s,ia)%values )
!        if (thisA%occupations(s,ia)%values(i) > 0 ) then
!          ii = ii + 1 
!          occupiedOrbitals(ii,1) = i 
!        end if
!      end do 
!
!      ii = 0
!      do i = 1, size(thisB%occupations(s,ib)%values )
!        if (thisB%occupations(s,ib)%values(i) > 0 ) then
!          ii = ii + 1 
!          occupiedOrbitals(ii,2) = i 
!        end if
!      end do 

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

      print *, "occ 1", occupiedOrbitals(:,1) 
      print *, "occ 2", occupiedOrbitals(:,2) 

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

      print *, "factor", factor

    end do

  end subroutine Configuration_checkMaximumCoincidence

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

    integer :: i, numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    print *, "Configuration"
    print *, "-------------"

    do i=1, numberOfSpecies
       print *, "For specie ", MolecularSystem_getNameOfSpecie ( i )
       print *, "Excitations: ", this%order%values(i)
       print *, "Occupations"
       call Vector_show ( this%occupations(i,1) )
       call Vector_show ( this%occupations(i,2) )
       print *, "Ndeterminants: ",this%nDeterminants
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
