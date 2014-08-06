!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Torsions_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Vertex_
  use Edges_
  use Angles_
  use MMCommons_
  use MatrixInteger_
  use Vector_
  use Exception_
  implicit none

  type , public :: Torsions

     integer :: numberOfTorsions
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: phi(:)
     integer, allocatable :: type(:)
     ! real(8), allocatable :: angleType(:)
     ! real(8), allocatable :: idealTheta(:)
     ! real(8), allocatable :: forceConstant(:)
     ! real(8), allocatable :: cosTheta(:)
     ! real(8), allocatable :: cosIdealTheta(:)
     ! real(8), allocatable :: sinTheta(:)
     ! real(8), allocatable :: sinIdealTheta(:)
     ! real(8), allocatable :: bendingEnergy(:) !! Kcal/mol
     ! real(8), allocatable :: bendingEnergyKJ(:) !! KJ/mol

  end type Torsions


       public :: &
            Torsions_constructor, &
            Torsions_getType
            ! Torsions_getBendingEnergies

contains

  subroutine Torsions_constructor( this, vertices, bonds, angle )
    implicit none
    type(Torsions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    integer :: i, j


    
    call MMCommons_constructor( MolecularSystem_instance )
    this%numberOfTorsions = size(MolecularSystem_instance%intCoordinates%dihedralsAngleValue%values)

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfTorsions, 4 )
    allocate( this%phi( this%numberOfTorsions ) )
    allocate( this%type( this%numberOfTorsions ) )
    ! allocate( this%idealTheta( this%numberOfTorsions ) )
    ! allocate( this%cosTheta( this%numberOfTorsions ) )
    ! allocate( this%cosIdealTheta( this%numberOfTorsions ) )
    ! allocate( this%sinTheta( this%numberOfTorsions ) )
    ! allocate( this%sinIdealTheta( this%numberOfTorsions ) )

    do i=1,this%numberOfTorsions
       do j=1,4
          this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForDihedrals%values(i,j)
       end do
       this%type(i) = Torsions_getType(this%connectionMatrix%values(i,2), this%connectionMatrix%values(i,3), vertices)
       this%phi(i) = MolecularSystem_instance%intCoordinates%dihedralsAngleValue%values(i)
    !    this%idealTheta(i) = vertices%angleValence(this%connectionMatrix%values(i,2))
    !    this%cosTheta(i) = cos(this%theta(i)*0.01745329251)
    !    this%cosIdealTheta(i) = cos(this%idealTheta(i)*0.01745329251)
    !    this%sinTheta(i) = sin(this%theta(i)*0.01745329251)
    !    this%sinIdealTheta(i) = sin(this%idealTheta(i)*0.01745329251)
    end do

    ! call Torsions_getForceConstants(this, vertices, bonds)

    ! call Torsions_getBendingEnergies(this, vertices)

  end subroutine Torsions_constructor

  function Torsions_getType(atomB, atomC, vertices) result(output)
    implicit none
    integer, intent(in) :: atomB, atomC
    type(Vertex), intent(in) :: vertices
    integer :: output

    output = 1

  end function Torsions_getType

!   subroutine Torsions_getForceConstants(this, vertices, bonds)
!     implicit none
!     type(Torsions), intent(in out) :: this
!     type(Vertex), intent(in) :: vertices
!     type(Edges), intent(in) :: bonds
!     integer :: i
!     real(8) :: Zi, Zk !! Effective charges
!     real(8) :: rij, rjk, rik
!     real(8) :: lambda

! !! se calcula usando la correccion de openbabel
!     lambda = 664.12

!     allocate( this%forceConstant( this%numberOfTorsions ) )

!     do i=1, this%numberOfTorsions
!        rij = Edges_getDistance(bonds,this%connectionMatrix%values(i,1),this%connectionMatrix%values(i,2))
!        rjk = Edges_getDistance(bonds,this%connectionMatrix%values(i,2),this%connectionMatrix%values(i,3))
!        rik = sqrt((rij**2)+(rjk**2)-2*rij*rjk*this%cosIdealTheta(i))
!        Zi=vertices%effectiveCharge(this%connectionMatrix%values(i,1))
!        Zk=vertices%effectiveCharge(this%connectionMatrix%values(i,3))
!        !! equation 13
!        this%forceConstant(i) = lambda*((Zi*Zk)/(rik**5))*(3*rij*rjk*(1-(this%cosIdealTheta(i)**2))-(rik**2)*this%cosIdealTheta(i))
!     end do

!   end subroutine Torsions_getForceConstants

!   subroutine Torsions_getBendingEnergies(this, vertices)
!     implicit none
!     type(Torsions), intent(in out) :: this
!     type(Vertex), intent(in) :: vertices
!     integer :: i
!     integer :: centralAtom
!     real :: coeff0, coeff1, coeff2


!     allocate( this%bendingEnergy( this%numberOfTorsions ) )
!     allocate( this%bendingEnergyKJ( this%numberOfTorsions ) )

!     do i=1, this%numberOfTorsions
!        centralAtom = this%connectionMatrix%values(i,2)
!        !! Caso lineal
!        if( vertices%connectivity(centralAtom) == 2 .AND. vertices%angleValence(centralAtom) == 180.0 ) then
!           this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))
!        !! Caso Trigonal plana
!        else if( vertices%connectivity(centralAtom) == 3 .AND. vertices%angleValence(centralAtom) == 120.0 ) then
!           this%bendingEnergy(i) = (this%forceConstant(i)/4.5)*(1.0 + (1.0 + this%cosTheta(i))*(4.0*this%cosTheta(i)))
!        !! Caso cuadrado planar y octaedrico
!        else if( (vertices%connectivity(centralAtom) == 4 .AND. vertices%angleValence(centralAtom) == 90.0) .OR. &
!             (vertices%connectivity(centralAtom) == 6 .AND. vertices%angleValence(centralAtom) == 90.0) .OR. &
!             (vertices%connectivity(centralAtom) == 6 .AND. vertices%angleValence(centralAtom) == 180.0)) then
!           this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))*this%cosTheta(i)*this%cosTheta(i)
!        !! Caso bipiramidal pentagonal (IF7)
!        else if( vertices%connectivity(centralAtom) == 7 ) then
!           this%bendingEnergy(i) = this%forceConstant(i)*1.0*&
!                (this%cosTheta(i)-0.30901699)*(this%cosTheta(i)-0.30901699)*&
!                (this%cosTheta(i)+0.80901699)*(this%cosTheta(i)+0.80901699)
!       !! Caso general sp3
!       else 
!          coeff2 = 1.0 / (4.0 * this%sinIdealTheta(i) * this%sinIdealTheta(i))
!          coeff1 = -4.0 * coeff2 * this%cosIdealTheta(i)
!          coeff0 = coeff2*(2.0*this%cosIdealTheta(i)*this%cosIdealTheta(i) + 1.0)
!          this%bendingEnergy(i) = this%forceConstant(i)*(coeff0 + coeff1*this%cosTheta(i) + coeff2*(2.0*this%cosTheta(i)*this%cosTheta(i)-1.0)) !! identidad cos2*Theta  
!        end if
!        this%bendingEnergyKJ(i) = this%bendingEnergy(i)*4.1868 
!      end do

!   end subroutine Torsions_getBendingEnergies

  subroutine Torsions_exception( typeMessage, description, debugDescription)
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

  end subroutine Torsions_exception

end module Torsions_
