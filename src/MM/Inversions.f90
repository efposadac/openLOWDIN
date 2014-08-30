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
module Inversions_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Vertex_
  use Edges_
  use Angles_
  use ListInteger_
  use MatrixInteger_
  use Vector_
  use Exception_
  implicit none

  type , public :: Inversions

     integer :: numberOfInversions
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: omega(:)
     real(8), allocatable :: C0(:)
     real(8), allocatable :: C1(:)
     real(8), allocatable :: C2(:)
     real(8), allocatable :: forceConstant(:)
     real(8), allocatable :: inversionEnergy(:) !! Kcal/mol
     real(8), allocatable :: inversionEnergyKJ(:) !! KJ/mol
     logical :: hasInversions

  end type Inversions


       public :: &
            Inversions_constructor, &
            Inversions_getConnectionMatrix, &
            Inversions_getNeighbors, &
            Inversions_getOmega, &
            Inversions_getPhi, &
            Inversions_getConstants

contains

  subroutine Inversions_constructor( this, vertices, bonds, angle )
    implicit none
    type(Inversions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    type(Angles), intent(in) :: angle
    integer :: i, j


    call Inversions_getConnectionMatrix(this, vertices, bonds)

    this%hasInversions= .false.
    if(this%numberOfInversions > 0) then
       this%hasInversions= .true.
    end if

    if(this%hasInversions) then
       call Inversions_getOmega(this, vertices, angle)
       call Inversions_getConstants(this, vertices)
       call Inversions_getInversionEnergies(this)
    end if

  end subroutine Inversions_constructor

  subroutine Inversions_getConnectionMatrix(this, vertices, bonds)
    implicit none
    type(Inversions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    integer, allocatable :: neighbors(:)
    integer :: i, inversionsCounter
    type(ListInteger) :: centralAtom

    call ListInteger_constructor( centralAtom, ssize=-1 )


    inversionsCounter = 0
    do i=1,vertices%numberOfVertices
       if(trim(vertices%type(i)) == "C_R" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "C_2" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "N_R" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "N_3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "N_2" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "P_3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "As3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "Sb3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       else if(trim(vertices%type(i)) == "Bi3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
          inversionsCounter = inversionsCounter + 1
       end if
    end do

    this%numberOfInversions = inversionsCounter
    
    if(this%numberOfInversions>0) then
       call MatrixInteger_constructor( this%connectionMatrix, this%numberOfInversions, 4 )

       this%connectionMatrix%values(:,1) = centralAtom%data(:)

       do i=1,this%numberOfInversions
          call Inversions_getNeighbors(neighbors, centralAtom%data(i), bonds)
          this%connectionMatrix%values(i,2) = neighbors(1)
          this%connectionMatrix%values(i,3) = neighbors(2)
          this%connectionMatrix%values(i,4) = neighbors(3)
       end do
    end if

  end subroutine Inversions_getConnectionMatrix

  subroutine Inversions_getNeighbors(neighbors, AtomI, bonds)
    implicit none
    integer, allocatable, intent(out) :: neighbors(:)
    integer, intent(in) :: AtomI
    type(Edges), intent(in) :: bonds
    integer :: i, j

    allocate( neighbors( 3 ) )

    j = 1
    do i=1,bonds%numberOfEdges
       if(bonds%connectionMatrix%values(i,1) == AtomI) then
          neighbors(j) = bonds%connectionMatrix%values(i,2)
          j = j + 1
       else if(bonds%connectionMatrix%values(i,2) == AtomI) then
          neighbors(j) = bonds%connectionMatrix%values(i,1)
          j = j + 1
       end if
    end do

  end subroutine Inversions_getNeighbors

  subroutine Inversions_getOmega(this, vertices, angle)
    implicit none
    type(Inversions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Angles), intent(in) :: angle
    integer :: i, atom1, atom2, atom3, atom4
    type(Vector) :: R1, R2, R3, R4
    type(Vector) :: R12, R13, R14
    type(Vector) :: unitR12, unitR13, unitR14
    type(Vector) :: crossVector1, crossVector2, crossVector3, perpenU1, perpenU2, perpenU3
    real(8) :: phi1, phi2, phi3, norm12, norm13, norm14
    real(8) :: sinPhi1, sinPhi2, sinPhi3, sinOmega1, sinOmega2, sinOmega3
    real(8) :: omega1, omega2, omega3


    ! 3
    !  \
    !   1----2   plane = 3-1-2
    !  /
    ! 4
    !
    ! R1 = [x1,y1,z1]
    ! R2 = [x2,y2,z2]
    ! R3 = [x3,y3,z3]
    ! R4 = [x4,y4,z4]    
    ! R12 = R2 - R1
    ! R13 = R3 - R1
    ! R14 = R4 - R1
    ! norm12 = sqrt(R12(x)**R12(x)+R12(y)**R12(y)+R12(z)**R12(z))
    ! norm13 = sqrt(R13(x)**R13(x)+R13(y)**R13(y)+R13(z)**R13(z))
    ! norm14 = sqrt(R14(x)**R14(x)+R14(y)**R14(y)+R14(z)**R14(z))
    ! unitR12 = R12/norm12
    ! unitR13 = R13/norm13
    ! unitR14 = R14/nomr14
    ! phi1 = angle 3-1-4
    ! phi2 = angle 2-1-4
    ! phi3 = angle 2-1-3
    ! crossVector1 = unitR13 X unitR14
    ! crossVector2 = unitR14 X unitR12
    ! crossVector3 = unitR12 X unitR13
    ! perpenU1 = crossVector1/sin(phi1)
    ! perpenU2 = crossVector2/sin(phi2)
    ! perpenU3 = crossVector3/sin(phi3)
    ! sinOmega1 = (unitR12)@perpenU1
    ! sinOmega2 = (unitR13)@perpenU2
    ! sinOmega3 = (unitR14)@perpenU3
    ! @ = dot product
    ! X = cross product

    call Vector_constructor(R1, 3)
    call Vector_constructor(R2, 3)
    call Vector_constructor(R3, 3)
    call Vector_constructor(R4, 3)
    call Vector_constructor(R12, 3)
    call Vector_constructor(R13, 3)
    call Vector_constructor(R14, 3)
    call Vector_constructor(unitR12, 3)
    call Vector_constructor(unitR13, 3)
    call Vector_constructor(unitR14, 3)
    call Vector_constructor(crossVector1, 3)
    call Vector_constructor(crossVector2, 3)
    call Vector_constructor(crossVector3, 3)
    call Vector_constructor(perpenU1, 3)
    call Vector_constructor(perpenU2, 3)
    call Vector_constructor(perpenU3, 3)
    
    
    if(allocated(this%omega)) deallocate(this%omega)
    allocate( this%omega( this%numberOfInversions ) )
    
    do i=1,this%numberOfInversions
       atom1 = this%connectionMatrix%values(i,1)
       atom2 = this%connectionMatrix%values(i,2)
       atom3 = this%connectionMatrix%values(i,3)
       atom4 = this%connectionMatrix%values(i,4)

       R1%values = vertices%cartesianMatrix%values(atom1,:) * AMSTRONG
       R2%values = vertices%cartesianMatrix%values(atom2,:) * AMSTRONG
       R3%values = vertices%cartesianMatrix%values(atom3,:) * AMSTRONG
       R4%values = vertices%cartesianMatrix%values(atom4,:) * AMSTRONG

       R12%values = R2%values - R1%values
       R13%values = R3%values - R1%values
       R14%values = R4%values - R1%values

       norm12 = sqrt(R12%values(1)*R12%values(1) + R12%values(2)*R12%values(2) + R12%values(3)*R12%values(3))
       norm13 = sqrt(R13%values(1)*R13%values(1) + R13%values(2)*R13%values(2) + R13%values(3)*R13%values(3))
       norm14 = sqrt(R14%values(1)*R14%values(1) + R14%values(2)*R14%values(2) + R14%values(3)*R14%values(3))
       
       unitR12%values = R12%values/norm12
       unitR13%values = R13%values/norm13
       unitR14%values = R14%values/norm14
       
       phi1 = Inversions_getPhi(angle, atom3, atom1, atom4)
       phi2 = Inversions_getPhi(angle, atom2, atom1, atom4)
       phi3 = Inversions_getPhi(angle, atom2, atom1, atom3)

       crossVector1 = Vector_cross(unitR13, unitR14)
       crossVector2 = Vector_cross(unitR14, unitR12)
       crossVector3 = Vector_cross(unitR12, unitR13)

       sinPhi1 = sin(phi1*0.01745329251)
       sinPhi2 = sin(phi2*0.01745329251)
       sinPhi3 = sin(phi3*0.01745329251)

       perpenU1 = Vector_scalarDiv(crossVector1, sinPhi1)
       perpenU2 = Vector_scalarDiv(crossVector2, sinPhi2)
       perpenU3 = Vector_scalarDiv(crossVector3, sinPhi3)

       sinOmega1 = Vector_dot(unitR12, perpenU1)
       sinOmega2 = Vector_dot(unitR13, perpenU2)
       sinOmega3 = Vector_dot(unitR14, perpenU3)

       omega1 = asin(sinOmega1)*57.2957795
       omega2 = asin(sinOmega2)*57.2957795
       omega3 = asin(sinOmega3)*57.2957795
       
       this%omega(i) = (omega1 + omega2 + omega3)/3
    end do

  end subroutine Inversions_getOmega

  function Inversions_getPhi(angle, atomA, atomB, atomC) result(output)
    implicit none
    type(Angles), intent(in) :: angle
    integer, intent(in) :: AtomB ! Central atom
    integer, intent(in) :: AtomA, AtomC
    real(8) :: output
    integer :: i

    do i=1,angle%numberOfAngles
       if(angle%connectionMatrix%values(i,1) == atomA .and. &
            angle%connectionMatrix%values(i,2) == atomB .and. &
            angle%connectionMatrix%values(i,3) == atomC) then
          output = angle%theta(i)
       else if(angle%connectionMatrix%values(i,1) == atomC .and. &
            angle%connectionMatrix%values(i,2) == atomB .and. &
            angle%connectionMatrix%values(i,3) == atomA) then
          output = angle%theta(i)
       end if
    end do

  end function Inversions_getPhi

  subroutine Inversions_getConstants(this, vertices)
    implicit none
    type(Inversions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i, centralAtom, neighborA, neighborB, neighborC
    real(8) :: idealOmega


    allocate( this%C0( this%numberOfInversions ))
    allocate( this%C1( this%numberOfInversions ))
    allocate( this%C2( this%numberOfInversions ))
    allocate( this%forceConstant( this%numberOfInversions ))

    do i=1,this%numberOfInversions
       centralAtom = this%connectionMatrix%values(i,1)
       if(trim(vertices%type(centralAtom)) == "N_3" .or. &
            trim(vertices%type(centralAtom)) == "N_2" .or. &
            trim(vertices%type(centralAtom)) == "N_R") then
          this%C0(i) = 1.0
          this%C1(i) = -1.0
          this%C2(i) = 0.0          
          this%forceConstant(i) = 6.0
       else if(trim(vertices%type(centralAtom)) == "C_R" .or. &
            trim(vertices%type(centralAtom)) == "C_2") then
          neighborA = this%connectionMatrix%values(i,2)
          neighborB = this%connectionMatrix%values(i,3)
          neighborC = this%connectionMatrix%values(i,4)
          this%C0(i) = 1.0
          this%C1(i) = -1.0
          this%C2(i) = 0.0
          if(trim(vertices%type(neighborA)) == "O_2" .or. &
               trim(vertices%type(neighborB)) == "O_2" .or. &
               trim(vertices%type(neighborC)) == "O_2") then
             this%forceConstant(i) = 50.0
          else
             this%forceConstant(i) = 6.0
          end if
       else if(trim(vertices%type(centralAtom)) == "P_3+3") then
          idealOmega = 84.4339*0.01745329251
          this%C1(i) = -4.0*cos(idealOmega*0.01745329251)
          this%C2(i) = 1.0
          this%C0(i) = -1.0*this%C1(i)*cos(idealOmega) + this%C2(i)*cos(2.0*idealOmega)
          this%forceConstant(i) = 22.0
       else if(trim(vertices%type(centralAtom)) == "As3+3") then
          idealOmega = 86.9735*0.01745329251
          this%C1(i) = -4.0*cos(idealOmega*0.01745329251)
          this%C2(i) = 1.0
          this%C0(i) = -1.0*this%C1(i)*cos(idealOmega) + this%C2(i)*cos(2.0*idealOmega)
          this%forceConstant(i) = 22.0
       else if(trim(vertices%type(centralAtom)) == "Sb3+3") then
          idealOmega = 87.7047*0.01745329251
          this%C1(i) = -4.0*cos(idealOmega*0.01745329251)
          this%C2(i) = 1.0
          this%C0(i) = -1.0*this%C1(i)*cos(idealOmega) + this%C2(i)*cos(2.0*idealOmega)
          this%forceConstant(i) = 22.0
       else if(trim(vertices%type(centralAtom)) == "Bi3+3") then
          idealOmega = 90.000*0.01745329251
          this%C1(i) = -4.0*cos(idealOmega*0.01745329251)
          this%C2(i) = 1.0
          this%C0(i) = -1.0*this%C1(i)*cos(idealOmega) + this%C2(i)*cos(2.0*idealOmega)
          this%forceConstant(i) = 22.0
       end if
    end do
    
  end subroutine Inversions_getConstants

  subroutine Inversions_getInversionEnergies(this)
    implicit none
    type(Inversions), intent(in out) :: this
    integer :: i

    allocate( this%inversionEnergy( this%numberOfInversions ))
    allocate( this%inversionEnergyKJ( this%numberOfInversions ))

    do i=1,this%numberOfInversions
       this%inversionEnergy(i) = this%forceConstant(i)*(this%C0(i) + &
            this%C1(i)*cos(this%omega(i)*0.01745329251) + &
            this%C2(i)*cos(2*this%omega(i)*0.01745329251))
       this%inversionEnergyKJ(i) = this%inversionEnergy(i)**4.1868
    end do
  end subroutine Inversions_getInversionEnergies

  subroutine Inversions_exception( typeMessage, description, debugDescription)
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

  end subroutine Inversions_exception

end module Inversions_
