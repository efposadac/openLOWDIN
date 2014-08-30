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
