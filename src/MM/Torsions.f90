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
     real(8), allocatable :: rotationalBarrier(:)
     real(8), allocatable :: idealPhi(:)
     real(8), allocatable :: order(:)
     real(8), allocatable :: torsionEnergy(:) !! Kcal/mol
     real(8), allocatable :: torsionEnergyKJ(:) !! KJ/mol
     logical :: hasTorsion

  end type Torsions


       public :: &
            Torsions_constructor, &
            Torsions_getConstants, &
            Torsions_getTorsionEnergies

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

    this%hasTorsion= .false.
    if(this%numberOfTorsions > 0) then
       this%hasTorsion= .true.
    end if

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfTorsions, 4 )
    allocate( this%phi( this%numberOfTorsions ) )

    do i=1,this%numberOfTorsions
       do j=1,4
          this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForDihedrals%values(i,j)
       end do
       this%phi(i) = MolecularSystem_instance%intCoordinates%dihedralsAngleValue%values(i)
    end do

    call Torsions_getConstants(this, vertices, bonds)

    call Torsions_getTorsionEnergies(this)

  end subroutine Torsions_constructor

  subroutine Torsions_getConstants(this, vertices, bonds)
    implicit none
    type(Torsions), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Edges), intent(in) :: bonds
    integer :: i, atomA, atomB, atomC, atomD
    logical :: isBgroupSixMember, isCgroupSixMember
    real(8) :: Vj, Vk !! Torsional barriers
    real(8) :: Uj,Uk !! TorsionalConstant
    real(8) :: bondOrder

    isBgroupSixMember = .false.
    isCgroupSixMember = .false.

    allocate(this%order(this%numberOfTorsions))
    allocate(this%idealPhi(this%numberOfTorsions))
    allocate(this%rotationalBarrier(this%numberOfTorsions))

    do i=1,this%numberOfTorsions
       atomA = this%connectionMatrix%values(i,1)
       atomB = this%connectionMatrix%values(i,2)
       atomC = this%connectionMatrix%values(i,3)
       atomD = this%connectionMatrix%values(i,4)
       if(vertices%hybridization(atomB)==3 .AND. vertices%hybridization(atomC)==3) then
          isBgroupSixMember = Torsions_isGroupSixMember(atomB, vertices)
          isCgroupSixMember = Torsions_isGroupSixMember(atomC, vertices)
          if( isBgroupSixMember .AND. isCgroupSixMember ) then
             if(vertices%charges(atomB)==8.0 .AND. vertices%charges(atomC)==8.0) then
                Vj = 2.0
                Vk = 2.0
             else if (vertices%charges(atomB)==8.0 .AND. vertices%charges(atomC)/=8.0) then
                Vj = 2.0
                Vk = 6.8
             else if (vertices%charges(atomB)/=8.0 .AND. vertices%charges(atomC)==8.0) then
                Vj = 6.8
                Vk = 2.0
             else 
                Vj = 6.8
                Vk = 6.8
             end if
             this%order(i) = 2.0
             this%idealPhi(i) = 90.0
          else
             Vj = vertices%torsionalBarrier(atomB)
             Vk = vertices%torsionalBarrier(atomC)
             this%order(i) = 3.0
             this%idealPhi(i) = 60.0 
          end if
          this%rotationalBarrier(i) = sqrt(Vj*Vk)
       else if(vertices%hybridization(atomB)==2 .AND. vertices%hybridization(atomC)==2) then
          this%order(i) = 2.0
          this%idealPhi(i) = 180.0
          bondOrder = Edges_getOrder(bonds, atomB, atomC)
          Uj = vertices%torsionalConstant(atomB) 
          Uk = vertices%torsionalConstant(atomC)
          this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
       else if(vertices%hybridization(atomB)==3 .AND. vertices%hybridization(atomC)==2) then
          isBgroupSixMember = Torsions_isGroupSixMember(atomB, vertices)
          isCgroupSixMember = Torsions_isGroupSixMember(atomB, vertices)
          if(vertices%hybridization(atomD)==2) then
             this%order(i) = 3.0
             this%idealPhi(i) = 180.0
             this%rotationalBarrier(i) = 2.0
          else if(isBgroupSixMember .and. (.not.isCgroupSixMember)) then
             this%order(i) = 2.0
             this%idealPhi(i) = 90.0
             bondOrder = Edges_getOrder(bonds, atomB, atomC)
             Uj = vertices%torsionalConstant(atomB) 
             Uk = vertices%torsionalConstant(atomC)
             this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
          else
             this%order(i) = 6.0
             this%idealPhi(i) = 0.0
             this%rotationalBarrier(i) = 1.0
          end if
       else if(vertices%hybridization(atomB)==2 .AND. vertices%hybridization(atomC)==3) then
          isBgroupSixMember = Torsions_isGroupSixMember(atomB, vertices)
          isCgroupSixMember = Torsions_isGroupSixMember(atomB, vertices)
          if(vertices%hybridization(atomA)==2) then
             this%order(i) = 3.0
             this%idealPhi(i) = 180.0
             this%rotationalBarrier(i) = 2.0
          else if(isCgroupSixMember .and. (.not.isBgroupSixMember)) then
             this%order(i) = 2.0
             this%idealPhi(i) = 90.0
             bondOrder = Edges_getOrder(bonds, atomB, atomC)
             Uj = vertices%torsionalConstant(atomB) 
             Uk = vertices%torsionalConstant(atomC)
             this%rotationalBarrier(i) = 5.0*sqrt(Uj*Uk)*(1.0+4.18*log(bondOrder))
          else
             this%order(i) = 6.0
             this%idealPhi(i) = 0.0
             this%rotationalBarrier(i) = 1.0
          end if
       else
          this%idealPhi(i) = 0.0
          this%rotationalBarrier(i) = 0.0
       end if
    end do
  end subroutine Torsions_getConstants

  function Torsions_isGroupSixMember(atom, vertices) result(output)
    implicit none
    integer, intent(in) :: atom
    type(Vertex), intent(in) :: vertices
    logical :: output

    output = .false.
    if( vertices%charges(atom)==8.0 .OR. &
         vertices%charges(atom)==16.0 .OR. &
         vertices%charges(atom)==34.0 .OR. &
         vertices%charges(atom)==52.0 .OR. &
         vertices%charges(atom)==84.0 ) then
       output = .true.
    end if

  end function Torsions_isGroupSixMember

  subroutine Torsions_getTorsionEnergies(this)
    implicit none
    type(Torsions), intent(in out) :: this
    integer :: i
    real(8) :: cosIdealTorsion, cosTorsion

    allocate(this%torsionEnergy(this%numberOfTorsions))
    allocate(this%torsionEnergyKJ(this%numberOfTorsions))

    do i=1,this%numberOfTorsions
       cosIdealTorsion = cos(this%order(i)*this%idealPhi(i)*0.01745329251)
       cosTorsion = cos(this%order(i)*this%phi(i)*0.01745329251)
       this%torsionEnergy(i) = 0.5*this%rotationalBarrier(i)*(1-cosIdealTorsion*cosTorsion)
       this%torsionEnergyKJ(i) = this%torsionEnergy(i)*4.1868
    end do
  end subroutine Torsions_getTorsionEnergies

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
