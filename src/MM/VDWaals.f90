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
            ! VDWaals_getConstants, &
            ! VDWaals_getTorsionEnergies

contains

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
    this%distance = vdwDistance%data * AMSTRONG

    call MatrixInteger_constructor( this%connectionMatrix, numberofVDWdistances, 2 )
    this%connectionMatrix%values(:,1) = atomA%data(:)
    this%connectionMatrix%values(:,2) = atomB%data(:)

    this%numberOfVDWaals = numberofVDWdistances

    call List_destructor(vdwDistance)
    call ListInteger_destructor(atomA)
    call ListInteger_destructor(atomB)

  end subroutine VDWaals_getDistance

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
