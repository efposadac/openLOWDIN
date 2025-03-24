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
!!        This module creates a class with the information of the Van der Waals interactions in the system
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
module VDWaals_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Vertex_
  use Edges_
  use Angles_
  use MMCommons_
  use MatrixInteger_
  use List_
  use ListInteger_
  use Vector_
  use Exception_
  implicit none

  type , public :: VDWaals

     integer :: numberOfVDWaals
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: distance(:)
     real(8), allocatable :: idealDistance(:)
     real(8), allocatable :: wellDepth(:)
     real(8), allocatable :: VDWEnergy(:) !! Kcal/mol
     real(8), allocatable :: VDWEnergyKJ(:) !! KJ/mol
     logical :: VDW

  end type VDWaals


       public :: &
            VDWaals_constructor, &
            VDWaals_getDistance, &
            VDWaals_getIdealDistance, &
            VDWaals_getWellDepth, &
            VDWaals_isVDW, &
            VDWaals_getVDWEnergies

contains

  !>
  !! @brief Defines the class constructor
  !! The Van der Waals energy excludes the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Van der Waals interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  !! @param numberOfVDWaals INTEGER number of Van der Waals interactions
  !! @param connectionMatrix INTEGER ARRAY with the information about the vertices with Van der Waals interaction
  !! @param distance REAL(8) ARRAY with the Van der Waals distance interaction
  !! @param idealDistance REAL(8) ARRAY with the ideal Van der Waals distance interaction
  !! @param wellDepth REAL(8) ARRAY with the Van der Waals well depth
  !! @param VDWEnergy REAL(8) ARRAY with the Van der Waals energies (kcal/mol) of the system
  !! @param VDWEnergyKJ REAL(8) ARRAY with the Van der Waals energies (kJ/mol) of the system
  !! @param VDW LOGICAL returns .true. if the system has Van der Waals interactions
  subroutine VDWaals_constructor( this, vertices, bonds, angle )
    implicit none
    type(VDWaals), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    integer :: i, j
 
    call VDWaals_getDistance(this, vertices, bonds, angle)
    call VDWaals_getIdealDistance(this, vertices)
    call VDWaals_getWellDepth(this, vertices)
    call VDWaals_getVDWEnergies(this)

  end subroutine VDWaals_constructor

  !>
  !! @brief This routine calculates the Van der Waals distance excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Van der Waals interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  subroutine VDWaals_getDistance(this, vertices, bonds, angle)
    implicit none
    type(VDWaals), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    type(List) :: vdwDistance
    type(ListInteger) :: atomA, atomB
    integer :: i, j
    logical :: isVDW
    real(8) :: separationOfCenters 
    integer :: numberofVDWdistances


    call List_constructor( vdwDistance, ssize=-1 )
    call ListInteger_constructor( atomA, ssize=-1 )
    call ListInteger_constructor( atomB, ssize=-1 )

    do i=1,vertices%numberOfVertices
       do j=i+1,vertices%numberOfVertices
          isVDW = VDWaals_isVDW(bonds, angle, i, j)
          if( isVDW ) then
             separationOfCenters = dsqrt( sum( ( vertices%cartesianMatrix%values(i,:) - vertices%cartesianMatrix%values(j,:))**2.0 ) )
             call ListInteger_push_back(atomA, i)
             call ListInteger_push_back(atomB, j)
             call List_push_back( vdwDistance, separationOfCenters )
          end if
       end do
    end do


    numberofVDWdistances = ListInteger_size(atomA)
    
    this%VDW = .false.
    if(numberofVDWdistances > 1) then
       this%VDW = .true.
    end if

    allocate(this%distance(numberofVDWdistances))
    this%distance = vdwDistance%data * ANGSTROM

    call MatrixInteger_constructor( this%connectionMatrix, numberofVDWdistances, 2 )
    this%connectionMatrix%values(:,1) = atomA%data(:)
    this%connectionMatrix%values(:,2) = atomB%data(:)

    this%numberOfVDWaals = numberofVDWdistances

    call List_destructor(vdwDistance)
    call ListInteger_destructor(atomA)
    call ListInteger_destructor(atomB)

  end subroutine VDWaals_getDistance

  !>
  !! @brief This function evaluates if two atoms has Van der Waals interactions
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  !! @param [in] atomA INTEGER first atom to evaluate
  !! @param [in] atomB INTEGER second atom to evaluate
  !! @return [out] output LOGICAL return .true. if the atoms has Van der Waals interactions
  function VDWaals_isVDW(bonds, angle, atomA, atomB) result(output)
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

  end function VDWaals_isVDW

  !>
  !! @brief This routine calculates the ideal Van der Waals distance excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Van der Waals interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @note The ideal distances are calculated using the equation 21b in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[
  !! x_{ij}=\sqrt{x_{i}{x_{j}}}
  !! \f]
  !! where: \n
  !! - \f$x_{i}\f$ and \f$x_{j}\f$ are the atomic Van der Waals distances, are parameters in UFF
  subroutine VDWaals_getIdealDistance(this, vertices)
    implicit none
    type(VDWaals), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    integer :: atomA, atomB
    real(8) :: Xi, Xj !! atomic Van der Waals distance

    allocate(this%idealDistance(this%numberOfVDWaals))

    do i=1,this%numberOfVDWaals
       atomA = this%connectionMatrix%values(i,1)
       atomB = this%connectionMatrix%values(i,2)
       Xi = vertices%distanceVdW(atomA)
       Xj = vertices%distanceVdW(atomB)
       this%idealDistance(i) = sqrt(Xi*Xj)
    end do

  end subroutine VDWaals_getIdealDistance

  !>
  !! @brief This routine calculates the well depth for Van der Waals interactions excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Van der Waals interactions
  !! @param [in] vertices Class with the information of the vertices
  !! @note The well depths are calculated using the equation 22 in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[
  !! D_{ij}=\sqrt{D_{i}{D_{j}}}
  !! \f]
  !! where: \n
  !! - \f$D_{i}\f$ and \f$D_{j}\f$ are the atomic Van der Waals energy, are parameters in UFF
  subroutine VDWaals_getWellDepth(this, vertices)
    implicit none
    type(VDWaals), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    integer :: atomA, atomB
    real(8) :: Di, Dj !! atomic Van der Waals energy

    allocate(this%wellDepth(this%numberOfVDWaals))

    do i=1,this%numberOfVDWaals
       atomA = this%connectionMatrix%values(i,1)
       atomB = this%connectionMatrix%values(i,2)
       Di = vertices%energyVdW(atomA)
       Dj = vertices%energyVdW(atomB)
       this%wellDepth(i) = sqrt(Di*Dj)
    end do

  end subroutine VDWaals_getWellDepth

  !>
  !! @brief This routine calculates the Van der Waals energies excluding the interactions for atoms
  !! that are bonded to each other (1,2 interactions) and bonded to a common atom (1,3 interactions)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the Van der Waals interactions
  !! @note The energies are calculated using the equation 20 in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[
  !! E_{vdw} = D_{ij}\left[\left(\frac{x_{ij}}{x}\right)^{12}-2\left(\frac{x_{ij}}{x}\right)^{6}\right]
  !! \f]
  !! where: \n
  !! - \f$D_{ij}\f$ is the well depth energy
  !! - \f$x\f$ is the Van der Waals distance
  !! - \f$x_{ij}\f$ is the ideal Van der Waals distance
  subroutine VDWaals_getVDWEnergies(this)
    implicit none
    type(VDWaals), intent(in out) :: this
    real(8) :: twelvepow, sixpow
    integer :: i

    allocate(this%VDWEnergy(this%numberOfVDWaals))
    allocate(this%VDWEnergyKJ(this%numberOfVDWaals))

    do i=1,this%numberOfVDWaals
       sixpow = (this%idealDistance(i)/this%distance(i))**6
       twelvepow = sixpow**2
       this%VDWEnergy(i) = this%wellDepth(i)*(twelvepow - 2*sixpow)
       this%VDWEnergyKJ(i) = this%VDWEnergy(i)*4.1868
    end do

  end subroutine VDWaals_getVDWEnergies

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine VDWaals_exception( typeMessage, description, debugDescription)
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

  end subroutine VDWaals_exception

end module VDWaals_
