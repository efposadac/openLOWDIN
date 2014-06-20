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
module MMCommons_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use MatrixInteger_
  use Matrix_
  use Vector_
  use InternalCoordinates_
  use Exception_
  implicit none


  public :: &
       MMCommons_constructor, &
       MMCommons_getConnectivity, &
       MMCommons_getConnectivityMatrix, &
       MMCommons_getAngleAverage

contains

  subroutine MMCommons_constructor(this)
    implicit none
    type(MolecularSystem) :: this

    call InternalCoordinates_constructor(  this%intCoordinates )
    call InternalCoordinates_obtainCoordinates(  this%intCoordinates )

  end subroutine MMCommons_constructor

  function MMCommons_getConnectivity( this, atomIdx ) result( output )
    implicit none
    type(MolecularSystem) :: this
    integer, intent(in) :: atomIdx
    integer :: i
    integer :: j
    integer :: output

    output = 0

    do i=1,size(this%intCoordinates%distanceBondValue%values)
       do j=1,2
          if ( this%intCoordinates%connectionMatrixForBonds%values(i,j) == atomIdx ) then
             output=output+1
          end if
       end do
    end do
    
  end function MMCommons_getConnectivity

  function MMCommons_getAngleAverage ( this, atomIdx ) result( output )
    implicit none
    type(MolecularSystem) :: this
    integer, intent(in) :: atomIdx
    integer :: i
    integer :: numberOfAnglesForAtom
    real(8) :: output
    real(8) :: angle

    numberOfAnglesForAtom = 0
    angle = 0.00000000

    do i=1,size(this%intCoordinates%angleOfBondValue%values)
       if( this%intCoordinates%connectionMatrixForAngles%values(i,2) == atomIdx ) then
          angle = angle + this%intCoordinates%angleOfBondValue%values(i)
          numberOfAnglesForAtom = numberOfAnglesForAtom + 1
       end if
    end do

    output = angle/numberOfAnglesForAtom

  end function MMCommons_getAngleAverage


  recursive subroutine MMCommons_sort(a,na)

    ! DUMMY ARGUMENTS
    integer, intent(in) :: nA
    integer, dimension(nA,2), intent (in out) :: A

    ! LOCAL VARIABLES
    integer :: left, right
    real :: random
    real :: pivot
    integer, dimension(1,2) :: temp
    integer :: marker

    if (nA > 1) then

       call random_number(random)
       pivot = A(int(random*real(nA-1))+1,2)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = nA + 1
       do while (left < right)
          right = right - 1
          do while (A(right,2) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (A(left,2) < pivot)
             left = left + 1
          end do
          if (left < right) then
             temp(1,:) = A(left,:)
             A(left,:) = A(right,:)
             A(right,:) = temp(1,:)
          end if
       end do
       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if


       call MMCommons_sort(A(:marker-1,:),marker-1)
       call MMCommons_sort(A(marker:,:),nA-marker+1 )

    end if

  end subroutine MMCommons_sort

  subroutine MMCommons_getConnectivityMatrix ( this, outputSize, output )
    implicit none
    type(MolecularSystem) :: this
    integer, intent(in) :: outputSize
    type(MatrixInteger), intent(in out) :: output
    integer :: i
    integer :: connectivity
    
    call MatrixInteger_constructor( output, outputSize, 2 )

    do i=1, outputSize
       connectivity = MMCommons_getConnectivity( this , i )
       output%values(i,1) = i
       output%values(i,2) = connectivity
    end do

  end subroutine MMCommons_getConnectivityMatrix

  function MMCommons_searchNeighbor( edges, edgesSize, Idx ) result(output)
    implicit none
    integer, intent(in) :: edgesSize
    type(MatrixInteger), allocatable :: edges(:)
    integer :: Idx
    integer :: output
    integer :: i
    integer :: j

    do i=1, edgesSize
       do j=1, 2
          if ( edges(i)%values(1,j) == Idx ) then
             if ( j == 1 ) then
                output = edges(i)%values(1,2)
             else 
                output = edges(i)%values(1,1)
             end if
          end if
       end do
    end do

  end function MMCommons_searchNeighbor

  subroutine MMCommons_searchEdgesRow( edges, edgesSize, Idx, output )
    implicit none
    integer, intent(in) :: edgesSize
    type(MatrixInteger), allocatable :: edges(:)
    integer :: Idx
    integer, allocatable, intent(out) :: output(:)
    integer :: i
    integer :: j
    integer :: outputSize
    integer :: position

    position = 1
    outputSize = 0
    do i=1, edgesSize
       do j=1, 2
          if ( edges(i)%values(1,j) == Idx ) then
                outputSize = outputSize+1
          end if
       end do
    end do
    
    allocate( output( outputSize ) )

    do i=1, edgesSize
       do j=1, 2
          if ( edges(i)%values(1,j) == Idx ) then
             output(position) = i
             position = position + 1
          end if
       end do
    end do

  end subroutine MMCommons_searchEdgesRow

  subroutine MMCommons_pruningGraph(this, numberOfCenters, edges, connectivityMatrix, bonds)
    implicit none
    type(MolecularSystem) :: this
    integer, intent(in) :: numberOfCenters
    type(MatrixInteger), intent(out)  :: connectivityMatrix
    type(MatrixInteger), allocatable, intent(out) :: edges(:)
    type(Vector), intent(out) :: bonds
    integer :: i
    integer :: j
    integer :: atomToRemove
    integer :: neighbor
    integer :: edgesSize
    integer :: minimum
    integer :: connectivitySize
    integer, allocatable :: edgesRow(:)
    integer :: row
    
    call MatrixInteger_constructor( connectivityMatrix, numberOfCenters, 2 )
    call MMCommons_getConnectivityMatrix( this, numberOfCenters, connectivityMatrix )

    edgesSize = size(this%intCoordinates%distanceBondValue%values)

    allocate( edges( edgesSize )  )
    do i=1,edgesSize
       call MatrixInteger_constructor( edges(i), 1, 2 )
    end do

    do i=1,edgesSize
       do j=1,2
          edges(i)%values(1,j) = this%intCoordinates%connectionMatrixForBonds%values(i,j)
       end do
    end do

    call Vector_constructor( bonds, edgesSize )
    do i=1,edgesSize
       do j=1,2
          bonds%values(i) = this%intCoordinates%distanceBondValue%values(i) * AMSTRONG
       end do
    end do

    call MMCommons_sort( connectivityMatrix%values, numberOfCenters )   

    minimum = minval( connectivityMatrix%values(:,2), numberOfCenters )
    do while ( minimum == 1 )
       atomToRemove = connectivityMatrix%values(1,1)
       neighbor = MMCommons_searchNeighbor( edges, edgesSize, atomToRemove )
       call MMCommons_searchEdgesRow( edges, edgesSize, atomToRemove, edgesRow )
       row = edgesRow(1)
       call MatrixInteger_removeRow( connectivityMatrix, 1 )
       call MMCommons_removeEdge( edges, row, edgesSize )       
       call Vector_removeElement( bonds, row )
       connectivitySize = size(connectivityMatrix%values)/2
       edgesSize = size(edges)
       do i=1,connectivitySize
          if ( connectivityMatrix%values(i,1) == neighbor ) then
             connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
          end if
       end do
       call MMCommons_sort( connectivityMatrix%values, connectivitySize )
       minimum = minval( connectivityMatrix%values(:,2), connectivitySize )
    end do

  end subroutine MMCommons_pruningGraph

  subroutine MMCommons_removeEdge( this, numberOfElement, edgesSize )
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    integer, intent(in) :: numberOfElement
    integer, intent(in) :: edgesSize
    type(Matrixinteger), allocatable :: auxArray(:)
    integer :: i
    integer :: j

    if (numberOfElement <= edgesSize ) then

       allocate( auxArray(edgesSize-1) )
       do i=1,edgesSize-1
          call MatrixInteger_constructor( auxArray(i), 1, 2 )
       end do

       do i=1,numberOfElement-1
          do j=1,2
             auxArray(i)%values(1,j) = this(i)%values(1,j)
          end do
       end do

       do i=numberOfElement,edgesSize-1
          do j=1,2
             auxArray(i)%values(1,j) = this(i+1)%values(1,j)
          end do
       end do

       deallocate( this )
       allocate( this(edgesSize-1) )
       do i=1,edgesSize-1
          call MatrixInteger_constructor( this(i), 1, 2 )
       end do
       
       do i=1,edgesSize-1
          do j=1,2
             this(i)%values(1,j) = auxArray(i)%values(1,j)
          end do
       end do

       deallocate( auxArray )

    end if
    
  end subroutine MMCommons_removeEdge

end module MMCommons_
