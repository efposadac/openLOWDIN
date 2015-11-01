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

module Coupled_
  use Vector_
  use Exception_
  use MolecularSystem_
  implicit none

  !>
  !! @brief Description
  !!
  !! @author Alejandro
  !!
  !! <b> Creation data : </b> 2015
  !!
  !! <b> History change: </b>
  !!
  !<
  type, public :: Coupled
     character(20) :: name
     logical :: isInstanced
     type(Vector) :: coefficients
     type(Vector), allocatable :: occupations(:)
     type(Vector) :: order !! 1=single, 2=double, 3=triple, etc
  end type Coupled

  public :: &
       Coupled_constructor, &
       Coupled_destructor, &
       Coupled_show

  private		
contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine Coupled_constructor(this,occupiedCode,unoccupiedCode,order,numberOfCouplings)
    implicit none
    type(Coupled) :: this
    type(Vector) :: order
    type(Vector) :: occupiedCode
    type(Vector) :: unoccupiedCode
    integer :: numberOfCouplings

    integer :: numberOfOccupiedOrbitals 
    integer :: numberOfOrbitals 
    integer :: i,j
    integer :: numberOfSpecies
    integer :: div1
    integer :: div2
    integer :: lambda !Ocupation per orbital

    call Vector_constructor( this%coefficients, numberOfCouplings, 0.0_8 )
    call Vector_copyConstructor( this%order, order )

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate ( this%occupations(numberOfSpecies) )

    do i=1, numberOfSpecies
       !spin orbitals not spatial orbitals
       lambda=MolecularSystem_getLambda(i)
       numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda
       numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(i)*lambda

       call Vector_constructor ( this%occupations(i), numberOfOrbitals , 0.0_8 )

       do j=1, numberOfOccupiedOrbitals
          this%occupations(i)%values(j)=1
       end do

       div1= int(occupiedCode%values(i))
       div2= int(unoccupiedCode%values(i))

       do j=this%order%values(i), 1, -1 
          this%occupations(i)%values( MOD ( div1, 1024 ) ) = this%occupations(i)%values( MOD ( div1, 1024 ) ) - 1
          div1= div1/1024
          this%occupations(i)%values( MOD ( div2, 1024 ) ) = this%occupations(i)%values( MOD ( div2, 1024 ) ) + 1
          div2= div2/1024
       end do
       
    end do

    this%isInstanced = .true.

  end subroutine Coupled_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine Coupled_destructor(this)
    implicit none
    type(Coupled) :: this
    integer :: i, numberOfSpecies

    call Vector_destructor( this%coefficients )
    call Vector_destructor( this%order )

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do i=1, numberOfSpecies
       call Vector_destructor ( this%occupations(i) )
      
    end do

    deallocate ( this%occupations )

    this%isInstanced = .false.

  end subroutine Coupled_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine Coupled_show(this)
    implicit none
    type(Coupled) :: this

    integer :: i, numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    print *, "Coupled"
    print *, "-------------"

    do i=1, numberOfSpecies
       print *, "For specie ", MolecularSystem_getNameOfSpecie ( i )
       print *, "Excitations: ", this%order%values(i)
       print *, "Occupations"
       call Vector_show ( this%occupations(i) )
    end do

  end subroutine Coupled_show

  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  !<
  function Coupled_isInstanced( this ) result( output )
    implicit  none
    type(Coupled), intent(in) :: this
    logical :: output

    output = this%isInstanced

  end function Coupled_isInstanced

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine Coupled_exception( typeMessage, description, debugDescription)
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

  end subroutine Coupled_exception

end module Coupled_
