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
!! @brief Molecular Mechanics program.
!!        This module creates a class with the information of the angles in the system
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Angles_
  use CONTROL_
  use MolecularSystem_
  use Vertex_
  use Edges_
  use Matrix_
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

  !>
  !! @brief Defines the class constructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the angles
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param numberOfAngles INTEGER number of angles in the system
  !! @param connectionMatrix INTEGER ARRAY with the information about the vertices in a angle
  !! @param theta REAL ARRAY with the angles(Degrees) of the system
  !! @param idealTheta REAL ARRAY with the ideal angles(Degrees) of the system
  !! @param forceConstant REAL ARRAY with the force constants of the system
  !! @param cosTheta REAL ARRAY with cos(theta) of the system
  !! @param cosIdealTheta REAL ARRAY with cos(idealTheta) of the system
  !! @param sinTheta REAL ARRAY with sin(theta) of the system
  !! @param sinIdealTheta REAL ARRAY with sin(idealTheta) of the system
  !! @param bendingEnergy REAL ARRAY with the bending energies (kcal/mol) of the system
  !! @param stretchingEnergyKJ REAL ARRAY with the bending energies (kJ/mol) of the system
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

  !>
  !! @brief This routine calculates the force constants using UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the angles
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @note The Force constants are calculated using the equation 13 in Rappe et. al. paper (1992)
  !! with the correction proposed by Marcus G. Martin and implemented on TOWHEE (http://towhee.sourceforge.net/forcefields/uff.html) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[ 
  !! k_{ijk}=664.12\frac{Z_{i}^{*}Z_{k}^{*}}{r_{ik}^{5}}(3r_{ij}r_{jk}(1-\cos^{2}\theta_{0})-r_{ik}^{2}\cos\theta_{0})
  !! \f]
  !! where: \n
  !! - \f$Z_{i}^{*}\f$ and \f$Z_{j}^{*}\f$ are the effective atomic charges, parameters in the UFF
  !! - \f$r_{ij}\f$, \f$r_{ik}\f$, \f$r_{jk}\f$ are the ideal distances, parameters in the UFF
  !! - \f$\theta_{0}\f$ is the ideal angle, parameter in the UFF
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

  !>
  !! @brief This routine calculates the bending energies using UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the angles
  !! @param [in] vertices Class with the information of the vertices
  !! @note The Bending energies are calculated using the equation 10-12 in Rappe et. al. paper (1992)
  !! with the correction proposed by Marcus G. Martin and implemented on TOWHEE (http://towhee.sourceforge.net/forcefields/uff.html) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \n
  !! <b>Linear case:</b>
  !! \f[ 
  !! E_{\theta} = K_{ijk}(1+\cos\theta)
  !! \f]
  !! <b>Trigonal-planar case:</b>
  !! \f[ 
  !! E_{\theta} = \frac{K_{ijk}}{4.5}(1+(1+\cos\theta)(4\cos\theta))
  !! \f]
  !! <b>Square-planar and octahedral case:</b>
  !! \f[ 
  !! E_{\theta} = K_{ijk}(1+\cos\theta)\cos^{2}\theta
  !! \f]
  !! <b> General case:</b>
  !! \f[ 
  !! E_{\theta} = K_{ijk}(C_{0}+C_{1}\cos\theta+C_{2}\cos2\theta)
  !! \f]
  !! \f[ 
  !! C_{2}= 1/(4\sin^{2}\theta_{0})
  !! \f]
  !! \f[ 
  !! C_{1}= -4C_{2}\cos\theta_{0}$
  !! \f]
  !! \f[ 
  !! C_{0}= C_{2}(2\cos^{2}\theta_{0}+1)
  !! \f]
  !! where: \n
  !! - \f$K_{ijk}\f$ is the force constant
  !! - \f$\theta\f$ is the angle
  !! - \f$\theta_{0}\f$ is the ideal angle, parameter in the UFF
  subroutine Angles_getBendingEnergies(this, vertices)
    implicit none
    type(Angles), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    integer :: centralAtom
    real :: coeff0, coeff1, coeff2
    real :: angle
    integer :: conectivity


    allocate( this%bendingEnergy( this%numberOfAngles ) )
    allocate( this%bendingEnergyKJ( this%numberOfAngles ) )

    do i=1, this%numberOfAngles
       centralAtom = this%connectionMatrix%values(i,2)
       conectivity = vertices%connectivity(centralAtom)
       angle = vertices%angleValence(centralAtom)
       !! Caso lineal
       if( conectivity == 2 .AND. angle == 180.0 ) then
          this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))
       !! Caso Trigonal plana
       else if( conectivity == 3 .AND. angle == 120.0 ) then
          this%bendingEnergy(i) = (this%forceConstant(i)/4.5)*(1.0 + (1.0 + this%cosTheta(i))*(4.0*this%cosTheta(i)))
       !! Caso Trigonal plana del Nitrogeno
       else if( conectivity == 3 .AND. angle == 111.2 ) then
          this%bendingEnergy(i) = (this%forceConstant(i)/4.5)*(1.0 + (1.0 + this%cosTheta(i))*(4.0*this%cosTheta(i)))
       !! Caso cuadrado planar y octaedrico
       else if( (conectivity == 4 .AND. angle == 90.0) .OR. &
            (conectivity == 6 .AND. angle == 90.0) .OR. &
            (conectivity == 6 .AND. angle == 180.0)) then
          this%bendingEnergy(i) = this%forceConstant(i)*(1.0 + this%cosTheta(i))*this%cosTheta(i)*this%cosTheta(i)
       !! Caso bipiramidal pentagonal (IF7)
       else if( conectivity == 7 ) then
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

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
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
