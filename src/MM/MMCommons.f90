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
!!        This module contains common functions for Molecular Mechanics calculations 
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions using Universal Force Field has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module MMCommons_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use AtomicElement_
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
       MMCommons_isOrganometallic, &
       MMCommons_getAngleAverage, &
       MMCommons_getValences, &
       MMCommons_sort, &
       MMCommons_removeEdge, &
       MMCommons_pruningGraph, &
       MMCommons_searchEdgesRow, &
       MMCommons_searchNeighbor

contains

  !>
  !! @brief Defines the class constructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine MMCommons_constructor(this)
    implicit none
    type(MolecularSystem) :: this

    call InternalCoordinates_constructor(  this%intCoordinates )
    call InternalCoordinates_obtainCoordinates(  this%intCoordinates )

  end subroutine MMCommons_constructor

  !>
  !! @brief This function search for the connectivity for one atom
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class MolecularSystem with all information about the system
  !! @param [in] atomIdx INTEGER Id that indentify the atom to check his connectivity
  !! @return [out] output INTEGER connectivity of the atom
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

  !>
  !! @brief This function calculate the angle average for one atom
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class MolecularSystem with all information about the system
  !! @param [in] atomIdx INTEGER Id that indentify the atom to check his connectivity
  !! @return [out] output REAL(8) Angle average of the atom
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

  !>
  !! @brief Recursive subroutine realizes an ascending sort of an integer array
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] a INTEGER Array to sort, change on the output
  !! @param [in] na INTEGER with the dimension of the array A
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

  !>
  !! @brief Get connectivity Matrix for whole system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class MolecularSystem with all information about the system
  !! @param [in] outputSize INTEGER size of the output array
  !! @param [in,out] output INTEGER ARRAY with connectivity information of the system
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

  !>
  !! @brief Search the neighbors of one atom in the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] edges INTEGER ARRAY with the edges(bonds) information
  !! @param [in] edgesSize INTEGER size of the edges array
  !! @param [in] outputSize INTEGER size of the output array, generally is the connectivity of the atom
  !! @param [in] Idx ID of the atom 
  !! @param [in,out] output INTEGER ARRAY with neighbors information of the atom
  subroutine MMCommons_searchNeighbor( edges, edgesSize, outputSize, Idx, output ) 
    implicit none
    integer, intent(in) :: edgesSize
    integer, intent(in) :: outputSize
    type(MatrixInteger), allocatable :: edges(:)
    integer :: Idx
    integer, allocatable, intent(out) :: output(:)
    integer :: i
    integer :: j
    integer :: row

    allocate( output( outputSize ) )

    row=1
    do i=1, edgesSize
       do j=1, 2
          if ( edges(i)%values(1,j) == Idx ) then
             if ( j == 1 ) then
                output(row) = edges(i)%values(1,2)
                row = row + 1
             else 
                output(row) = edges(i)%values(1,1)
                row = row + 1
             end if
          end if
       end do
    end do

  end subroutine MMCommons_searchNeighbor

  !>
  !! @brief Search the rows in edges array where is the bond information of an specific atom
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] edges INTEGER ARRAY with the edges(bonds) information
  !! @param [in] edgesSize INTEGER size of the edges array
  !! @param [in] Idx ID of the atom 
  !! @return [out] output INTEGER ARRAY with number of each row
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
    integer :: numberOfColumns

    position = 1
    outputSize = 0
    do i=1, edgesSize
       numberOfColumns = size(edges(i)%values)
       do j=1, numberOfColumns
          if ( edges(i)%values(1,j) == Idx ) then
                outputSize = outputSize+1
          end if
       end do
    end do
    
    allocate( output( outputSize ) )

    do i=1, edgesSize
       numberOfColumns = size(edges(i)%values)
       do j=1, numberOfColumns
          if ( edges(i)%values(1,j) == Idx ) then
             output(position) = i
             position = position + 1
          end if
       end do
    end do

  end subroutine MMCommons_searchEdgesRow

  !>
  !! @brief This subroutine cretes a pruned graph, all atoms with connectivity == 1 are cut of the original graph, the connectivity is updated in each iteration step
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this MolecularSystem class with all information of the system
  !! @param [in] numberOfCenters INTEGER number of atoms in the system
  !! @return [out] edges INTEGER ARRAY with edges information of the pruned graph
  !! @return [out] connectivityMatrix INTEGER ARRAY with connectivity information of the pruned graph
  !! @return [out] bonds REAL(8) ARRAY with distances information of the pruned graph
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
    integer, allocatable :: neighbor(:)
    integer :: edgesSize
    integer :: numberOfColumns
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
          bonds%values(i) = this%intCoordinates%distanceBondValue%values(i) * ANGSTROM
       end do
    end do

    call MMCommons_sort( connectivityMatrix%values, numberOfCenters )   

    minimum = minval( connectivityMatrix%values(:,2), numberOfCenters )
    do while ( minimum == 1 )
       atomToRemove = connectivityMatrix%values(1,1)
       call MMCommons_searchNeighbor( edges, edgesSize, 1, atomToRemove, neighbor )
       call MMCommons_searchEdgesRow( edges, edgesSize, atomToRemove, edgesRow )
       row = edgesRow(1)
       call MatrixInteger_removeRow( connectivityMatrix, 1 )
       call MMCommons_removeEdge( edges, row, edgesSize )       
       call Vector_removeElement( bonds, row )
       connectivitySize = size(connectivityMatrix%values)/2
       edgesSize = size(edges)
       do i=1,connectivitySize
          if ( connectivityMatrix%values(i,1) == neighbor(1) ) then
             connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
          end if
       end do
       call MMCommons_sort( connectivityMatrix%values, connectivitySize )
       minimum = minval( connectivityMatrix%values(:,2), connectivitySize )
    end do
   
  end subroutine MMCommons_pruningGraph

  !>
  !! @brief This subroutine removes an edge information of one element of the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this INTEGER ARRAY with edges information, on input has size = edgesSize, on output has size = edgesSize - 1
  !! @param [in] numberOfElement INTEGER number of the element to remove
  !! @param [in] edgesSize INTEGER size of the edges array
  subroutine MMCommons_removeEdge( this, numberOfElement, edgesSize )
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    integer, intent(in) :: numberOfElement
    integer, intent(in) :: edgesSize
    type(Matrixinteger), allocatable :: auxArray(:)
    integer :: i
    integer :: j
    integer :: numberOfColumns


    if (numberOfElement <= edgesSize ) then

       allocate( auxArray(edgesSize-1) )
       ! do i=1,edgesSize-1
       !    call MatrixInteger_constructor( auxArray(i), 1, 2 )
       ! end do

       do i=1,numberOfElement-1
          numberOfColumns = size(this(i)%values)
          call MatrixInteger_constructor( auxArray(i), 1, numberOfColumns )
          do j=1,numberOfColumns
             auxArray(i)%values(1,j) = this(i)%values(1,j)
          end do
       end do

       do i=numberOfElement,edgesSize-1
          numberOfColumns = size(this(i+1)%values)
          call MatrixInteger_constructor( auxArray(i), 1, numberOfColumns )
          do j=1,numberOfColumns
             auxArray(i)%values(1,j) = this(i+1)%values(1,j)
          end do
       end do

       deallocate( this ) 
      allocate( this(edgesSize-1) )
       do i=1,edgesSize-1
          numberOfColumns = size(auxArray(i)%values)
          call MatrixInteger_constructor( this(i), 1, numberOfColumns )
       end do
       
       do i=1,edgesSize-1
          numberOfColumns = size(auxArray(i)%values)
          do j=1,numberOfColumns
             this(i)%values(1,j) = auxArray(i)%values(1,j)
          end do
       end do

       deallocate( auxArray )

    end if
    
  end subroutine MMCommons_removeEdge

  !>
  !! @brief This function evaluates if the neigbors of one atom are metals
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this MolecularSystem class with all information of the system
  !! @param [in] atomIdx INTEGER atom to evaluate if is organometallic
  !! @param [in] connectivity INTEGER connectivity of the atom
  !! @param [in] labelOfCenters CHARACTER ARRAY with the symbols of the atoms in the system
  !! @return [out] output LOGICAL .true. = the atom is organometallic, .false. = atom is not organometallic
  function MMCommons_isOrganometallic( this, atomIdx, connectivity, labelOfCenters ) result( output )
    implicit none
    type(MolecularSystem) :: this
    integer, intent(in) :: atomIdx
    integer, intent(in) :: connectivity
    integer :: i
    integer :: j
    logical :: output
    type(Vector) :: neighbor
    integer :: neighborIdx
    integer :: edgesSize, row
    type(AtomicElement) :: element
    character(10), allocatable :: labelOfCenters(:)

    output = .false.

    edgesSize = size(this%intCoordinates%distanceBondValue%values)

    call Vector_constructor( neighbor, connectivity )

    row = 1
    do i=1,edgesSize
       if ( this%intCoordinates%connectionMatrixForBonds%values(i,1) == atomIdx ) then
          neighbor%values(row) = this%intCoordinates%connectionMatrixForBonds%values(i,2)
          row = row +1
       else if ( this%intCoordinates%connectionMatrixForBonds%values(i,2) == atomIdx ) then
          neighbor%values(row) = this%intCoordinates%connectionMatrixForBonds%values(i,1)
          row = row +1
       end if
    end do

    do j=1,connectivity
       neighborIdx = neighbor%values(j)
       call AtomicElement_load ( element, trim( labelOfCenters(neighborIdx) ), 0 )
       if (element%atomicNumber >= 21 .and. element%atomicNumber <= 31) then
          output = .true.
       else if (element%atomicNumber >= 39 .and. element%atomicNumber <= 50) then
          output = .true.
       else if (element%atomicNumber >= 57 .and. element%atomicNumber <= 83) then
          output = .true.
       else if (element%atomicNumber >= 89) then
          output = .true.
       end if
    end do

  end function MMCommons_isOrganometallic

  !>
  !! @brief This subroutine obtain the valences of all atom in the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this CHARACTER ARRAY with the symbols of the atoms in the system
  !! @param [in] numberOfCenters INTEGER number of atom in the system
  !! @param [in,out] valences REAL(8) ARRAY with the valences of all atoms
  subroutine MMCommons_getValences(this, numberOfCenters, valences)
    implicit none
    character(10), allocatable, intent(in) :: this(:)
    integer, intent(in) :: numberOfCenters
    real(8), allocatable, intent(in out) :: valences(:)
    integer :: i
    type(Exception) :: ex

    allocate( valences( numberOfCenters ) )

    do i=1, numberOfCenters
!!******************************************************************************
!! Se evalua el grupo 1
!!******************************************************************************
       if( trim( this(i) ) == "H" &
            .OR. trim( this(i) ) == "LI" &
            .OR. trim( this(i) ) == "NA" &
            .OR. trim( this(i) ) == "K" &
            .OR. trim( this(i) ) == "RB" &
            .OR. trim( this(i) ) == "CS" &
            .OR. trim( this(i) ) == "FR" ) then
          valences(i) = 1
!!******************************************************************************
!! Se evalua el grupo 2
!!******************************************************************************
       else if( trim( this(i) ) == "BE" &
            .OR. trim( this(i) ) == "MG" &
            .OR. trim( this(i) ) == "CA" &
            .OR. trim( this(i) ) == "SR" &
            .OR. trim( this(i) ) == "BA" &
            .OR. trim( this(i) ) == "RA" ) then
          valences(i) = 2
!!******************************************************************************
!! Se evalua el grupo 13
!!******************************************************************************
       else if( trim( this(i) ) == "B" &
            .OR. trim( this(i) ) == "AL" &
            .OR. trim( this(i) ) == "GA" &
            .OR. trim( this(i) ) == "IN" &
            .OR. trim( this(i) ) == "TL" ) then
          valences(i) = 3
!!******************************************************************************
!! Se evalua el grupo 14
!!******************************************************************************
       else if( trim( this(i) ) == "C" &
            .OR. trim( this(i) ) == "SI" &
            .OR. trim( this(i) ) == "GE" &
            .OR. trim( this(i) ) == "SN" &
            .OR. trim( this(i) ) == "PB") then
          valences(i) = 4
!!******************************************************************************
!! Se evalua el grupo 15
!!******************************************************************************
       else if( trim( this(i) ) == "N" &
            .OR. trim( this(i) ) == "P" &
            .OR. trim( this(i) ) == "AS" &
            .OR. trim( this(i) ) == "SB" &
            .OR. trim( this(i) ) == "BI") then
          valences(i) = 5
!!******************************************************************************
!! Se evalua el grupo 16
!!******************************************************************************
       else if( trim( this(i) ) == "O" ) then
          valences(i) = 2
       else if( trim( this(i) ) == "S" &
            .OR. trim( this(i) ) == "SE" &
            .OR. trim( this(i) ) == "TE" &
            .OR. trim( this(i) ) == "PO") then
          valences(i) = 6
!!******************************************************************************
!! Se evalua el grupo 17
!!******************************************************************************
       else if( trim( this(i) ) == "F" ) then
          valences(i) = 1
       else if( trim( this(i) ) == "CL" &
            .OR. trim( this(i) ) == "BR" &
            .OR. trim( this(i) ) == "I" &
            .OR. trim( this(i) ) == "AT") then
          valences(i) = 7
!!******************************************************************************
!! Se evalua el grupo 18
!!******************************************************************************
       else if( trim( this(i) ) == "HE" &
            .OR. trim( this(i) ) == "NE") then
          valences(i) = 0
       else if( trim( this(i) ) == "AR" &
            .OR. trim( this(i) ) == "KR") then
          valences(i) = 2
       else if( trim( this(i) ) == "XE" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "RN" ) then
          valences(i) = 6
!!******************************************************************************
!! Se evalua la primera serie de transicion 
!!******************************************************************************
       else if( trim( this(i) ) == "SC" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "TI" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "V" ) then
          valences(i) = 5
       else if( trim( this(i) ) == "CR" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "MN" ) then
          valences(i) = 7
       else if( trim( this(i) ) == "FE" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "CO" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "NI" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "CU" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "ZN" ) then
          valences(i) = 2
!!******************************************************************************
!! Se evalua la segunda serie de transicion 
!!******************************************************************************
       else if( trim( this(i) ) == "Y" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "ZR" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "NB" ) then
          valences(i) = 5
       else if( trim( this(i) ) == "MO" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "TC" ) then
          valences(i) = 7
       else if( trim( this(i) ) == "RU" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "RH" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "PD" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "AG" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "CD" ) then
          valences(i) = 2
!!******************************************************************************
!! Se evalua la tercera serie de transicion y lantanidos
!!******************************************************************************
       else if( trim( this(i) ) == "LA" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "CE" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "PR" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "ND" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "PM" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "SM" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "EU" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "GD" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "TB" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "DY" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "HO" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "ER" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "TM" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "YB" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "LU" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "HF" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "TA" ) then
          valences(i) = 5
       else if( trim( this(i) ) == "W" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "RE" ) then
          valences(i) = 7
       else if( trim( this(i) ) == "OS" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "IR" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "PT" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "AU" ) then
          valences(i) = 5
       else if( trim( this(i) ) == "HG" ) then
          valences(i) = 4
!!******************************************************************************
!! Se evalua los Actinidos
!!******************************************************************************
       else if( trim( this(i) ) == "AC" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "TH" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "PA" ) then
          valences(i) = 5
       else if( trim( this(i) ) == "U" ) then
          valences(i) = 6
       else if( trim( this(i) ) == "NP" ) then
          valences(i) = 7
       else if( trim( this(i) ) == "PU" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "AM" ) then
          valences(i) = 7
       else if( trim( this(i) ) == "CM" ) then
          valences(i) = 8
       else if( trim( this(i) ) == "BK" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "CF" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "ES" ) then
          valences(i) = 4
       else if( trim( this(i) ) == "FM" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "MD" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "NO" ) then
          valences(i) = 3
       else if( trim( this(i) ) == "LR" ) then
          valences(i) = 3
!!******************************************************************************
       else
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object MMCommons in getValences() function" )
          call Exception_setDescription( ex, "This Atom hasn't been implemented" )
          call Exception_show( ex )
       end if
    end do

  end subroutine MMCommons_getValences

  !>
  !! @brief Defines the exception of the class
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine MMCommons_exception( typeMessage, description, debugDescription)
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

  end subroutine MMCommons_exception

end module MMCommons_
