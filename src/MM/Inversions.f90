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
!!        This module creates a class with the information of the inversion angles(out of plane) in the system
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

  !>
  !! @brief Defines the class constructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the inversion angles
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
  !! @param [in] angle Class with the information of the angles
  !! @param numberOfInversions INTEGER number of inversion angles in the system
  !! @param connectionMatrix INTEGER ARRAY with the information about the vertices in a inversion angle
  !! @param omega REAL ARRAY with the inversion angles(Degrees) of the system
  !! @param C0 REAL ARRAY with parameter of the UFF
  !! @param C1 REAL ARRAY with parameter of the UFF
  !! @param C2 REAL ARRAY with parameter of the UFF
  !! @param forceConstant REAL ARRAY with the force constants in the system
  !! @param inversionEnergy REAL ARRAY with the inversion energies (kcal/mol) of the system
  !! @param inversionEnergyKJ REAL ARRAY with the inversion energies (kJ/mol) of the system
  !! @param hasInversions LOGICAL returns .true. if the system has inversion angles
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

  !>
  !! @brief This routine evaluates if the system has inversion angles and return the members of the angle
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the inversion angles
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] bonds Class with the information of the edges
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

  !>
  !! @brief This routine searches the neighbors of an atom in the inversion angle
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [out] neighbors INTEGER ARRAY with the neighbors of the atom
  !! @param [in] AtomI INTEGER atom to evaluate
  !! @param [in] bonds Class with the information of the edges
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

  !>
  !! @brief This routine calculates the inversion angles (omega) of the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the inversion angles
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] angle Class with the information of the angle
  !! @note The inversion angle (\f$\omega\f$) is calculated using the method proposed by Lee et al. (1999) \n
  !! Lee, S.H., Palmo, K., Krimm, S., <b>New out-of-plane angle and bond angle internal coordinates 
  !! and related potential energy functions for molecular mechanics and dynamics simulations</b>,
  !! J. Comput. Chem., 20, 10, 1067--1084, 1999 \n
  !! \f[
  !! \sin\omega_{1} = e_{41}\cdot\left(\frac{e_{42}\times e_{43}}{\sin\phi_{1}}\right)
  !! \f]
  !! 3 (\f$e_{43}\f$) \n
  !! &nbsp;&nbsp;\ \n
  !! &nbsp;&nbsp;&nbsp;4 \f$\cdots\f$ 1 (\f$e_{41}\f$)&nbsp;&nbsp; plane = 3-4-2 and \f$\phi\f$ angle of 3-4-2\n
  !! &nbsp;&nbsp;/ \n
  !! 2 (\f$e_{42}\f$)\n
  !! \n
  !! where: \n
  !! - \f$e_{43}\f$, \f$e_{41}\f$ and \f$e_{42}\f$ are the unit vectors
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


  !>
  !! @brief This function searches the angle (phi) for the plane A-B-C
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] angle Class with the information of the angles
  !! @param [in] atomA INTEGER first atom to evaluate
  !! @param [in] atomB INTEGER second atom(central atom) to evaluate
  !! @param [in] atomC INTEGER third atom to evaluate
  !! @return [out] output REAL(8) angle between A-B-C
  function Inversions_getPhi(angle, atomA, atomB, atomC) result(output)
    implicit none
    type(Angles), intent(in) :: angle
    integer, intent(in) :: atomB ! Central atom
    integer, intent(in) :: atomA, atomC
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

  !>
  !! @brief This routine calculates the constants with the UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the inversion angles
  !! @param [in] vertices Class with the information of the vertices
  !! @note The constants are calculated using the parameters in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \n
  !! For P, As, Sb and Bi the constants are calculated using the correction implemented on OpenBabel (http://openbabel.org/wiki/Main_Page): \n
  !! \f[
  !! C_{1} = -4.0\cos(\omega_{0})
  !! \f]
  !! \f[
  !! C_{2} = 1.0
  !! \f]
  !! \f[
  !! C_{0} = -C_{1}\cos(\omega_{0}) + C_{2}\cos(2.0\omega_{0})
  !! \f]
  !! where: \n
  !! - \f$\omega_{0}\f$ = 84.4339 for P
  !! - \f$\omega_{0}\f$ = 86.9735 for As
  !! - \f$\omega_{0}\f$ = 87.7047 for Sb
  !! - \f$\omega_{0}\f$ = 90.0000 for Bi
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

  !>
  !! @brief This routine calculates the inversion energies with the UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the inversion angles
  !! @note The inversion energies are calculated using the equation 18 in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[
  !! E_{\omega} = K_{ijkl}(C_{0}+C_{1}\cos\omega_{ijkl}+C_{2}\cos2\omega_{ijkl})
  !! \f]
  !! where: \n
  !! - \f$\omega_{ijkl}\f$ is the inversion angle
  !! - \f$C_{0}\f$, \f$C_{1}\f$ and \f$C_{2}\f$ are parameters in UFF
  !! - \f$K_{ijkl}\f$ is the force constant and it is a parameter in the UFF
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

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
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
