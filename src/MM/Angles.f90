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
module Angles_
  use CONTROL_
  use MolecularSystem_
  use Vertex_
  use Edges_
  use MMCommons_
  use MatrixInteger_
  use Exception_
  implicit none

  type , public :: Angles

     integer :: numberOfAngles
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: theta(:)
     real(8), allocatable :: idealTheta(:)
     real(8), allocatable :: forceConstant(:)
     real(8), allocatable :: cosTheta(:)
     real(8), allocatable :: cosIdealTheta(:)
     real(8), allocatable :: sinTheta(:)
     real(8), allocatable :: sinIdealTheta(:)
     real(8), allocatable :: bendingEnergy(:) !! Kcal/mol
     real(8), allocatable :: bendingEnergyKJ(:) !! KJ/mol

  end type Angles


       public :: &
            Angles_constructor, &
            Angles_getBendingEnergies

contains

  subroutine Angles_constructor( this, vertices, bonds )
    implicit none
    type(Angles), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    integer :: i, j
    
    call MMCommons_constructor( MolecularSystem_instance )
    this%numberOfAngles = size(MolecularSystem_instance%intCoordinates%angleOfBondValue%values)

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfAngles, 3 )
    allocate( this%theta( this%numberOfAngles ) )
    allocate( this%idealTheta( this%numberOfAngles ) )
    allocate( this%cosTheta( this%numberOfAngles ) )
    allocate( this%cosIdealTheta( this%numberOfAngles ) )
    allocate( this%sinTheta( this%numberOfAngles ) )
    allocate( this%sinIdealTheta( this%numberOfAngles ) )

    do i=1,this%numberOfAngles
       do j=1,3
          this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForAngles%values(i,j)
       end do
       this%theta(i) = MolecularSystem_instance%intCoordinates%angleOfBondValue%values(i)
       this%idealTheta(i) = vertices%angleValence(this%connectionMatrix%values(i,2))
       this%cosTheta(i) = cos(this%theta(i)*0.01745329251)
       this%cosIdealTheta(i) = cos(this%idealTheta(i)*0.01745329251)
       this%sinTheta(i) = sin(this%theta(i)*0.01745329251)
       this%sinIdealTheta(i) = sin(this%idealTheta(i)*0.01745329251)
    end do

    call Angles_getForceConstants(this, vertices, bonds)

    call Angles_getBendingEnergies(this, vertices)

  end subroutine Angles_constructor

  subroutine Angles_getForceConstants(this, vertices, bonds)
    implicit none
    type(Angles), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    integer :: i
    real(8) :: Zi, Zk !! Effective charges
    real(8) :: rij, rjk, rik
    real(8) :: lambda

!! se calcula usando la correccion de openbabel
    lambda = 664.12

    allocate( this%forceConstant( this%numberOfAngles ) )

    do i=1, this%numberOfAngles
       rij = Edges_getDistance(bonds,this%connectionMatrix%values(i,1),this%connectionMatrix%values(i,2))
       rjk = Edges_getDistance(bonds,this%connectionMatrix%values(i,2),this%connectionMatrix%values(i,3))
       rik = sqrt((rij**2)+(rjk**2)-2*rij*rjk*this%cosIdealTheta(i))
       Zi=vertices%effectiveCharge(this%connectionMatrix%values(i,1))
       Zk=vertices%effectiveCharge(this%connectionMatrix%values(i,3))
       !! equation 13
       this%forceConstant(i) = lambda*((Zi*Zk)/(rik**5))*(3*rij*rjk*(1-(this%cosIdealTheta(i)**2))-(rik**2)*this%cosIdealTheta(i))
    end do

  end subroutine Angles_getForceConstants

  subroutine Angles_getBendingEnergies(this, vertices)
    implicit none
    type(Angles), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    integer :: centralAtom
    real :: coeff0, coeff1, coeff2


    allocate( this%bendingEnergy( this%numberOfAngles ) )
    allocate( this%bendingEnergyKJ( this%numberOfAngles ) )

    do i=1, this%numberOfAngles
       centralAtom = this%connectionMatrix%values(i,2)
       !! Caso lineal
       if( vertices%connectivity(centralAtom) == 2 .AND. vertices%angleValence(centralAtom) == 180.0 ) then
          this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))
       !! Caso Trigonal plana
       else if( vertices%connectivity(centralAtom) == 3 .AND. vertices%angleValence(centralAtom) == 120.0 ) then
          this%bendingEnergy(i) = (this%forceConstant(i)/4.5)*(1.0 + (1.0 + this%cosTheta(i))*(4.0*this%cosTheta(i)))
       !! Caso cuadrado planar y octaedrico
       else if( (vertices%connectivity(centralAtom) == 4 .AND. vertices%angleValence(centralAtom) == 90.0) .OR. &
            (vertices%connectivity(centralAtom) == 6 .AND. vertices%angleValence(centralAtom) == 90.0) .OR. &
            (vertices%connectivity(centralAtom) == 6 .AND. vertices%angleValence(centralAtom) == 180.0)) then
          this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))*this%cosTheta(i)*this%cosTheta(i)
       !! Caso bipiramidal pentagonal (IF7)
       else if( vertices%connectivity(centralAtom) == 7 ) then
          this%bendingEnergy(i) = this%forceConstant(i)*1.0*&
               (this%cosTheta(i)-0.30901699)*(this%cosTheta(i)-0.30901699)*&
               (this%cosTheta(i)+0.80901699)*(this%cosTheta(i)+0.80901699)
      !! Caso general sp3
      else 
         coeff2 = 1.0 / (4.0 * this%sinIdealTheta(i) * this%sinIdealTheta(i))
         coeff1 = -4.0 * coeff2 * this%cosIdealTheta(i)
         coeff0 = coeff2*(2.0*this%cosIdealTheta(i)*this%cosIdealTheta(i) + 1.0)
         this%bendingEnergy(i) = this%forceConstant(i)*(coeff0 + coeff1*this%cosTheta(i) + coeff2*(2.0*this%cosTheta(i)*this%cosTheta(i)-1.0)) !! identidad cos2*Theta  
       end if
       this%bendingEnergyKJ(i) = this%bendingEnergy(i)*4.1868 
     end do

  end subroutine Angles_getBendingEnergies

  subroutine Angles_exception( typeMessage, description, debugDescription)
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

  end subroutine Angles_exception

end module Angles_
