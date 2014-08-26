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
  ! use MMCommons_
  use MatrixInteger_
  ! use Vector_
  use Exception_
  implicit none

  type , public :: Inversions

     integer :: numberOfInversions
     type(MatrixInteger) :: connectionMatrix
     ! real(8), allocatable :: phi(:)
     ! real(8), allocatable :: rotationalBarrier(:)
     ! real(8), allocatable :: idealPhi(:)
     ! real(8), allocatable :: order(:)
     ! real(8), allocatable :: torsionEnergy(:) !! Kcal/mol
     ! real(8), allocatable :: torsionEnergyKJ(:) !! KJ/mol
     logical :: hasInversions

  end type Inversions


       public :: &
            Inversions_constructor, &
            Inversions_getConnectionMatrix, &
            Inversions_getNeighbors

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

    ! call MatrixInteger_constructor( this%connectionMatrix, this%numberOfInversions, 4 )
    ! allocate( this%phi( this%numberOfInversions ) )

    ! do i=1,this%numberOfInversions
    !    do j=1,4
    !       this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForDihedrals%values(i,j)
    !    end do
    !    this%phi(i) = MolecularSystem_instance%intCoordinates%dihedralsAngleValue%values(i)
    ! end do

    ! call Inversions_getConstants(this, vertices, bonds)

    ! call Inversions_getTorsionEnergies(this)

  end subroutine Inversions_constructor

  subroutine Inversions_getConnectionMatrix(this, vertices, bonds)
    implicit none
    type(Inversions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    integer, allocatable :: neighbors(:)
    integer :: i
    type(ListInteger) :: centralAtom

    call ListInteger_constructor( centralAtom, ssize=-1 )

    do i=1,vertices%numberOfVertices
       if(trim(vertices%type(i)) == "C_R" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "C_2" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "N_R" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "N_3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "N_2" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "P_3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "As3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "Sb3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       else if(trim(vertices%type(i)) == "Bi3+3" .and. vertices%connectivity(i) == 3)then
          call ListInteger_push_back(centralAtom, i)
       end if
    end do

    this%numberOfInversions = ListInteger_size(centralAtom)

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfInversions, 4 )

    this%connectionMatrix%values(:,1) = centralAtom%data(:)

    do i=1,this%numberOfInversions
       call Inversions_getNeighbors(neighbors, centralAtom%data(i), bonds)
       this%connectionMatrix%values(i,2) = neighbors(1)
       this%connectionMatrix%values(i,3) = neighbors(2)
       this%connectionMatrix%values(i,4) = neighbors(3)
       write(*,"(T20,5I)") i, this%connectionMatrix%values(i,:)
    end do

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

  ! subroutine Inversions_getConstants(this, vertices, bonds)
  !   implicit none
  !   type(Inversions), intent(in out) :: this
  !   type(Vertex), intent(in) :: vertices
  !   type(Edges), intent(in) :: bonds
  !   integer :: i, atomA, atomB, atomC, atomD
  !   logical :: isBgroupSixMember, isCgroupSixMember
  !   real(8) :: Vj, Vk !! Torsional barriers
  !   real(8) :: Uj,Uk !! TorsionalConstant
  !   real(8) :: bondOrder

  !   isBgroupSixMember = .false.
  !   isCgroupSixMember = .false.

  !   allocate(this%order(this%numberOfInversions))
  !   allocate(this%idealPhi(this%numberOfInversions))
  !   allocate(this%rotationalBarrier(this%numberOfInversions))

  !   do i=1,this%numberOfInversions
  !      atomA = this%connectionMatrix%values(i,1)
  !      atomB = this%connectionMatrix%values(i,2)
  !      atomC = this%connectionMatrix%values(i,3)
  !      atomD = this%connectionMatrix%values(i,4)
  !      if(vertices%hybridization(atomB)==3 .AND. vertices%hybridization(atomC)==3) then
  !         isBgroupSixMember = Inversions_isGroupSixMember(atomB, vertices)
  !         isCgroupSixMember = Inversions_isGroupSixMember(atomC, vertices)
  !         if( isBgroupSixMember .AND. isCgroupSixMember ) then
  !            if(vertices%charges(atomB)==8.0 .AND. vertices%charges(atomC)==8.0) then
  !               Vj = 2.0
  !               Vk = 2.0
  !            else if (vertices%charges(atomB)==8.0 .AND. vertices%charges(atomC)/=8.0) then
  !               Vj = 2.0
  !               Vk = 6.8
  !            else if (vertices%charges(atomB)/=8.0 .AND. vertices%charges(atomC)==8.0) then
  !               Vj = 6.8
  !               Vk = 2.0
  !            else 
  !               Vj = 6.8
  !               Vk = 6.8
  !            end if
  !            this%order(i) = 2.0
  !            this%idealPhi(i) = 90.0
  !         else
  !            Vj = vertices%torsionalBarrier(atomB)
  !            Vk = vertices%torsionalBarrier(atomC)
  !            this%order(i) = 3.0
  !            this%idealPhi(i) = 60.0
  !         end if
  !         this%rotationalBarrier(i) = sqrt(Vj*Vk)
  !      else if(vertices%hybridization(atomB)==2 .AND. vertices%hybridization(atomC)==2) then
  !         this%order(i) = 2.0
  !         this%idealPhi(i) = 180.0
  !         bondOrder = Edges_getOrder(bonds, atomB, atomC)
  !         Uj = vertices%torsionalConstant(atomB) 
  !         Uk = vertices%torsionalConstant(atomC)
  !         this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
  !      else if(vertices%hybridization(atomB)==3 .AND. vertices%hybridization(atomC)==2) then
  !         isBgroupSixMember = Inversions_isGroupSixMember(atomB, vertices)
  !         isCgroupSixMember = Inversions_isGroupSixMember(atomB, vertices)
  !         if(vertices%hybridization(atomD)==2) then
  !            this%order(i) = 3.0
  !            this%idealPhi(i) = 180.0
  !            this%rotationalBarrier(i) = 2.0
  !         else if(isBgroupSixMember .and. (.not.isCgroupSixMember)) then
  !            this%order(i) = 2.0
  !            this%idealPhi(i) = 90.0
  !            bondOrder = Edges_getOrder(bonds, atomB, atomC)
  !            Uj = vertices%torsionalConstant(atomB) 
  !            Uk = vertices%torsionalConstant(atomC)
  !            this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
  !         else
  !            this%order(i) = 6.0
  !            this%idealPhi(i) = 0.0
  !            this%rotationalBarrier(i) = 1.0
  !         end if
  !      else if(vertices%hybridization(atomB)==2 .AND. vertices%hybridization(atomC)==3) then
  !         isBgroupSixMember = Inversions_isGroupSixMember(atomB, vertices)
  !         isCgroupSixMember = Inversions_isGroupSixMember(atomB, vertices)
  !         if(vertices%hybridization(atomA)==2) then
  !            this%order(i) = 3.0
  !            this%idealPhi(i) = 180.0
  !            this%rotationalBarrier(i) = 2.0
  !         else if(isCgroupSixMember .and. (.not.isBgroupSixMember)) then
  !            this%order(i) = 2.0
  !            this%idealPhi(i) = 90.0
  !            bondOrder = Edges_getOrder(bonds, atomB, atomC)
  !            Uj = vertices%torsionalConstant(atomB) 
  !            Uk = vertices%torsionalConstant(atomC)
  !            this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
  !         else
  !            this%order(i) = 6.0
  !            this%idealPhi(i) = 0.0
  !            this%rotationalBarrier(i) = 1.0
  !         end if
  !      else
  !         this%idealPhi(i) = 0.0
  !         this%rotationalBarrier(i) = 0.0
  !      end if
  !   end do
  ! end subroutine Inversions_getConstants

  ! function Inversions_isGroupSixMember(atom, vertices) result(output)
  !   implicit none
  !   integer, intent(in) :: atom
  !   type(Vertex), intent(in) :: vertices
  !   logical :: output

  !   output = .false.
  !   if( vertices%charges(atom)==8.0 .OR. &
  !        vertices%charges(atom)==16.0 .OR. &
  !        vertices%charges(atom)==34.0 .OR. &
  !        vertices%charges(atom)==52.0 .OR. &
  !        vertices%charges(atom)==84.0 ) then
  !      output = .true.
  !   end if

  ! end function Inversions_isGroupSixMember

  ! subroutine Inversions_getTorsionEnergies(this)
  !   implicit none
  !   type(Inversions), intent(in out) :: this
  !   integer :: i
  !   real(8) :: cosIdealTorsion, cosTorsion

  !   allocate(this%torsionEnergy(this%numberOfInversions))
  !   allocate(this%torsionEnergyKJ(this%numberOfInversions))

  !   do i=1,this%numberOfInversions
  !      cosIdealTorsion = cos(this%order(i)*this%idealPhi(i)*0.01745329251)
  !      cosTorsion = cos(this%order(i)*this%phi(i)*0.01745329251)
  !      this%torsionEnergy(i) = 0.5*this%rotationalBarrier(i)*(1-cosIdealTorsion*cosTorsion)
  !      this%torsionEnergyKJ(i) = this%torsionEnergy(i)*4.1868
  !   end do
  ! end subroutine Inversions_getTorsionEnergies

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
