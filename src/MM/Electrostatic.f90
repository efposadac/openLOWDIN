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
!!        This module creates a class with the information of the Electrostatic interactions in the system, 
!! only the EQeq(Extended Charge Equilibration) approach for to calculate partial charges is implemented,
!! this energy calculation is optional and must be activated in the input like this:
!! <BLOCKQUOTE>
!! CONTROL \n
!! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;electrostaticMM = T \n
!! END CONTROL \n
!! </BLOCKQUOTE>
!! @note For reference see: \n
!! <b>Original Charge Equilibration approach (QEq):</b> \n
!! \n
!! Rappe, A.K., Goddard III, W.A., <b>Charge Equilibration for Molecular Dynamics Simulations</b>,
!! J. Phys. Chem., 95, 3358--3363, 1991 \n
!! \n
!! <b>Extended Charge Equilibration approach (EQeq)</b> \n
!! \n
!! Wilmer, C.E., Kim, K.C., Snurr, R.Q., <b>An Extended Charge Equilibration Method</b>,
!! J. Phys. Chem. Lett, 3, 2506--2511, 2012 
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!! @todo Is necessary to implement an iterative method for to calculate the partial charges
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Electrostatic_
  use CONTROL_
  use Vertex_
  use Edges_
  use Angles_
  use MatrixInteger_
  use List_
  use ListInteger_
  use ChargesEQeq_
  use Exception_
  implicit none

  type , public :: Electrostatic
     
     integer :: numberOfElectrostatics
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: partialCharge(:)
     real(8), allocatable :: distance(:)
     real(8), allocatable :: electrostaticEnergy(:) !! Kcal/mol
     real(8), allocatable :: electrostaticEnergyKJ(:) !! KJ/mol
     logical :: isElectrostatic

  end type Electrostatic


       public :: &
            Electrostatic_constructor, &
            Electrostatic_getDistance, &
            Electrostatic_isElectrostatic

contains

  !>
  !! @brief Defines the class constructor
  !! The Electrostatic energy excludes the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Electrostatic interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  !! @param numberOfElectrostatics INTEGER number of Electrostatic interactions
  !! @param connectionMatrix INTEGER ARRAY with the information about the vertices with Electrostatic interaction
  !! @param partialCharge REAL(8) ARRAY with partial charges of the system
  !! @param distance REAL(8) ARRAY with the Electrostatic distance interaction
  !! @param electrostaticEnergy REAL(8) ARRAY with the Electrostatic energies (kcal/mol) of the system
  !! @param electrostaticEnergyKJ REAL(8) ARRAY with the Electrostatic energies (kJ/mol) of the system
  !! @param isElectrostatic LOGICAL returns .true. if the system has Electrostatic interactions
  subroutine Electrostatic_constructor( this, vertices, bonds, angle )
    implicit none
    type(Electrostatic), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    real(8), allocatable :: partials(:)
    
    call Electrostatic_getDistance(this, vertices, bonds, angle)
    call ChargesEQeq_getCharges(partials, vertices)
    call Electrostatic_getEnergies(this, partials, vertices)

  end subroutine Electrostatic_constructor
  
  !>
  !! @brief This routine calculates the Electrostatic energies excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Electrostatic interactions
  !! @param [in] partials REAL(8) ARRAY with the partial charges of the system
  !! @param [in] vertices Class with the information of the vertices
  !! @note The energies are calculated using the equation 43 in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[
  !! E_{el} = \frac{1}{4\pi\epsilon_{0}}\left(\frac{Q_{i}Q_{j}}{\epsilon R_{ij}}\right)
  !! \f]
  !! where: \n
  !! - \f$\epsilon_{0}\f$ is the vacuum permittivity
  !! - \f$\epsilon\f$ is the relative permittivity
  !! - \f$R_{ij}\f$ is the distance between the charges
  !! - \f$Q_{i}\f$ and \f$Q_{j}\f$ are the partial charges
  subroutine Electrostatic_getEnergies(this, partials, vertices)
    implicit none
    type(Electrostatic), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    real(8), allocatable, intent(in) :: partials(:)
    integer :: i, atomA, atomB
    real(8) :: K

    K = 332.0637

    allocate( this%partialCharge( vertices%numberOfVertices ) )

    do i=1,vertices%numberOfVertices
       this%partialCharge(i) = partials(i)
    end do


    allocate( this%electrostaticEnergy( this%numberOfElectrostatics ) )
    allocate( this%electrostaticEnergyKJ( this%numberOfElectrostatics ) )

    do i=1,this%numberOfElectrostatics
       atomA = this%connectionMatrix%values(i,1)
       atomB = this%connectionMatrix%values(i,2)
       this%electrostaticEnergy(i) = K*((this%partialCharge(atomA)*this%partialCharge(atomB))/this%distance(i))
    end do

    do i=1,this%numberOfElectrostatics
       this%electrostaticEnergyKJ(i) = this%electrostaticEnergy(i)*4.1868
    end do

  end subroutine Electrostatic_getEnergies

  !>
  !! @brief This routine calculates the Electrostatic distance excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Electrostatic interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  subroutine Electrostatic_getDistance(this, vertices, bonds, angle)
    implicit none
    type(Electrostatic), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    type(List) :: electrostaticDistance
    type(ListInteger) :: atomA, atomB
    integer :: i, j
    logical :: isElect
    real(8) :: separationOfCenters 
    integer :: numberofElectrostaticDistances


    call List_constructor( electrostaticDistance, ssize=-1 )
    call ListInteger_constructor( atomA, ssize=-1 )
    call ListInteger_constructor( atomB, ssize=-1 )

    do i=1,vertices%numberOfVertices
       do j=i+1,vertices%numberOfVertices
          isElect = Electrostatic_isElectrostatic(bonds, angle, i, j)
          if( isElect ) then
             separationOfCenters = dsqrt( sum( ( vertices%cartesianMatrix%values(i,:) - vertices%cartesianMatrix%values(j,:))**2.0 ) )
             call ListInteger_push_back(atomA, i)
             call ListInteger_push_back(atomB, j)
             call List_push_back( electrostaticDistance, separationOfCenters )
          end if
       end do
    end do

    numberofElectrostaticDistances = ListInteger_size(atomA)
    
    this%isElectrostatic = .false.
    if(numberofElectrostaticDistances > 1) then
       this%isElectrostatic = .true.
    end if

    allocate(this%distance(numberofElectrostaticDistances))
    this%distance = electrostaticDistance%data * AMSTRONG

    call MatrixInteger_constructor( this%connectionMatrix, numberofElectrostaticDistances, 2 )
    this%connectionMatrix%values(:,1) = atomA%data(:)
    this%connectionMatrix%values(:,2) = atomB%data(:)

    this%numberOfElectrostatics = numberofElectrostaticDistances
    
    call List_destructor(electrostaticDistance)
    call ListInteger_destructor(atomA)
    call ListInteger_destructor(atomB)

  end subroutine Electrostatic_getDistance

  !>
  !! @brief This function evaluates if two atoms has Electrostatic interactions
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  !! @param [in] atomA INTEGER first atom to evaluate
  !! @param [in] atomB INTEGER second atom to evaluate
  !! @return [out] output LOGICAL return .true. if the atoms has Electrostatic interactions
  function Electrostatic_isElectrostatic(bonds, angle, atomA, atomB) result(output)
    implicit none
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    integer, intent(in) :: atomA, atomB
    logical :: isBond, isAngle, output
    integer :: i
    
    output = .false.
    isBond = .false.
    isAngle = .false.

    !! evaluando si es un enlace
    do i=1,bonds%numberOfEdges
       if(bonds%connectionMatrix%values(i,1) == atomA .and. bonds%connectionMatrix%values(i,2) == atomB) then
          isBond = .true.
       else if(bonds%connectionMatrix%values(i,1) == atomB .and. bonds%connectionMatrix%values(i,2) == atomA) then
          isBond = .true.
       end if
    end do

    !! evaluando si hacen parte de un angulo
    do i=1,angle%numberOfAngles
       if( angle%connectionMatrix%values(i,1) == atomA .and. angle%connectionMatrix%values(i,2) == atomB ) then
          isAngle = .true.
       else if( angle%connectionMatrix%values(i,1) == atomA .and. angle%connectionMatrix%values(i,3) == atomB ) then
          isAngle = .true.
       else if( angle%connectionMatrix%values(i,2) == atomA .and. angle%connectionMatrix%values(i,1) == atomB ) then
          isAngle = .true.
       else if( angle%connectionMatrix%values(i,2) == atomA .and. angle%connectionMatrix%values(i,3) == atomB ) then
          isAngle = .true.          
       else if( angle%connectionMatrix%values(i,3) == atomA .and. angle%connectionMatrix%values(i,1) == atomB ) then
          isAngle = .true.
       else if( angle%connectionMatrix%values(i,3) == atomA .and. angle%connectionMatrix%values(i,2) == atomB ) then
          isAngle = .true.
       end if
    end do

    if( isBond ) then
       output = .false.
    else if( isAngle ) then
       output = .false.
    else
       output = .true.
    end if

  end function Electrostatic_isElectrostatic

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine Electrostatic_exception( typeMessage, description, debugDescription)
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

  end subroutine Electrostatic_exception

end module Electrostatic_
