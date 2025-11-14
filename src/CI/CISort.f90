module CISort_
!! Author: Jorge Charry
!! Collection of sorting functions adapted for SCI program.
!! Inspired from CENCALC code https://github.com/dimassuarez/cencalc_quicksort/ 
  use MolecularSystem_
  use CIcore_
  use Vector_
  use Matrix_
  implicit none

  type, public :: CISort
    integer(1), allocatable :: tmp_Vector (:) ! species
    integer(1), allocatable :: pivot (:) ! species, orbs
    integer :: numberOfSpecies
    integer :: combinedNumberOfOrbitals
  end type CISort

  type(CISort) :: CISort_instance

contains

  !! allocating global array, to avoid construting them in OMP
  subroutine CISort_constructor()
    implicit none
    integer :: numberOfSpecies, spi, m
   
    numberOfSpecies = CIcore_instance%numberOfSpecies 

    CISort_instance%numberOfSpecies = numberOfSpecies

    !! auxiliary variables to map orbitals from vector to array location
    CISort_instance%combinedNumberOfOrbitals = sum(CIcore_instance%numberOfOrbitals%values(:))
  
    allocate ( CISort_instance%tmp_Vector(  CISort_instance%combinedNumberOfOrbitals ) )
    allocate ( CISort_instance%pivot( CISort_instance%combinedNumberOfOrbitals ) )

    CISort_instance%pivot = -1_1

  end subroutine CISort_constructor

  !! allocating global array, to avoid construting them in OMP
  subroutine CISort_destructor()
    implicit none
   
    deallocate ( CISort_instance%pivot )
    deallocate ( CISort_instance%tmp_Vector )

  end subroutine CISort_destructor

  !! Modified Recursive Quicksort routine
  recursive subroutine CISort_quicksort(array, indices, left_end, right_end)
    implicit none
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8) :: pivot
    integer, parameter :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) THEN
      !! Use interchange sort for small lists
      call CISort_interchange_sort(array, indices, left_end, right_end)

    else

      if ( left_end < right_end ) then
          !! Partition using the indexed version
          call CISort_partition_hoare(array, indices, left_end, right_end, pivot)
  
          !! Recursively sort the sub-arrays
          call CISort_quicksort(array, indices, left_end, pivot)
          call CISort_quicksort(array, indices, pivot + 1, right_end )
      end if
    endif

  end subroutine CISort_quicksort

  !! Modified Partition scheme to use swap_both
  subroutine CISort_partition_hoare(array, indices, left_end, right_end, pivot_index )
    implicit none
    integer(1), intent(inout) :: array(:,:) ! original array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8), intent(out) :: pivot_index
    integer(8) :: i, j
    integer :: spi, orb
    logical :: is_more    
    integer(8) :: mid_index

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( left_end + right_end ) /2
    CISort_instance%pivot(:) = array(:, mid_index )

    i = left_end - 1 ! Index of the smaller element
    j = right_end + 1 ! Index of the higher element

    do 
      do 
        i = i + 1
        is_more = .false.
        do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( CISort_instance%pivot(orb) > array(orb,i) ) then
            is_more = .true.
            exit
          else if ( CISort_instance%pivot(orb) < array(orb,i) ) then
            is_more = .false. 
            exit
          endif
        enddo
        if ( .not. is_more ) exit
      enddo

      do 
        j = j - 1
        is_more = .false.
        do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( CISort_instance%pivot(orb) < array(orb,j) ) then
            is_more = .true.
            exit
          else if ( CISort_instance%pivot(orb) > array(orb,j) ) then
            is_more = .false. 
            exit
          endif
        enddo
        if ( .not. is_more ) exit
      enddo

      if ( i >= j ) exit 

      !! Swap both the vector and the index
      call CISort_swapVector( array(:,i), array(:,j) )
      call CISort_swapIndex( indices%values(i), indices%values(j) )
    enddo
    pivot_index = j 

  end subroutine CISort_partition_hoare
 
  !! Lomuto partition scheme
  subroutine CISort_partition_lomuto(array, indices, left_end, right_end, pivot_index )
    implicit none
    integer(1), intent(inout) :: array(:,:) ! original array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8), intent(out) :: pivot_index
    integer(8) :: i, j
    integer :: spi, orb
    logical :: is_more    
    integer(8) :: mid_index

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( left_end + right_end ) /2
    ! swap mid point to the end, really? why?
    call CISort_swapVector( array(:,mid_index), array(:,right_end) )
    call CISort_swapIndex( indices%values(mid_index), indices%values(right_end ) )

    !! Select the last element (and its index) as the pivot
    CISort_instance%pivot(:) = array(:,right_end )

    i = left_end - 1 ! Index of the smaller element

    do j = left_end, right_end - 1
      is_more = .false.

      if ( array(1,j) == -1_1 ) cycle !! avoid zero, it means no configuration

       do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( array(orb,j) > CISort_instance%pivot(orb) ) then
            is_more = .true.
            exit
          else if ( array(orb,j) < CISort_instance%pivot(orb) ) then
            is_more = .false. 
            exit
          endif
       enddo
      
      if ( is_more ) then
        i = i + 1
        !! Swap both the vector and the index
        call CISort_swapVector( array(:,i), array(:,j) )
        call CISort_swapIndex( indices%values(i), indices%values(j) )
      end if
    end do
    
    !! Place the pivot in the correct sorted position (swap with arr(i+1) and index(i+1))
    call CISort_swapVector( array(:,i + 1), array(:,right_end) )
    call CISort_swapIndex( indices%values(i + 1), indices%values(right_end) )

    pivot_index = i + 1

  end subroutine CISort_partition_lomuto

  !! Sorting for small arrays, no more pivot
  subroutine CISort_interchange_sort(array, indices, left_end, right_end) 
    implicit none
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer :: i, j
    logical :: is_more    
    integer :: orb

    do i = left_end, right_end - 1
      if ( array(1,i) == -1_1 ) cycle !! avoid zero, it means no configuration
      do j = i+1, right_end
        is_more = .false.

        do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( array(orb,j) < array(orb,i) ) then
            is_more = .true.
            exit
          else if ( array(orb,j) > array(orb,i) ) then
            is_more = .false. 
            exit
          endif
       enddo
        
        if ( is_more ) then
          !! Swap both the vector and the index
          call CISort_swapVector( array(:,i), array(:,j) )
          call CISort_swapIndex( indices%values(i), indices%values(j) )
        end if

      enddo 
    enddo 

  end subroutine CISort_interchange_sort

  !! Swap both data and index elements
  subroutine CISort_swapVector(A, B)
    implicit none
    integer(1), intent(inout) :: A(:), B(:)

    !! Swap data
    CISort_instance%tmp_Vector(:) = A(:)
    A(:) = B(:)
    B(:) = CISort_instance%tmp_Vector(:)

  end subroutine CISort_swapVector

  !! Swap both data and index elements
  subroutine CISort_swapIndex(IndexA, IndexB)
    implicit none
    integer(8), intent(inout) :: IndexA, IndexB
    integer(8) :: tmpIndex

    !! Swap index
    tmpIndex = IndexA
    IndexA = IndexB
    IndexB = tmpIndex

  end subroutine CISort_swapIndex

  !! (The 'is_more' function for vector_t remains unchanged)
  function CISort_is_more(A, B ) result (is_more)
    implicit none
    logical :: is_more
    real(8), intent(in) :: A(:), B(:)
    integer :: i

    !! logic for lexicographical comparison
    do i = 1, CISort_instance%numberOfSpecies 
      if ( A(i) > B(i) ) then
        is_more = .true.
        return
      else if (A(i) < B(i) ) then
        is_more = .false.
        return
      end if
    end do
    is_more = .false.
    return

  end function CISort_is_more

  !! (The 'is_more' function for vector_t remains unchanged)
  function CISort_is_equal(A, B ) result (is_equal)
    implicit none
    logical :: is_equal
    real(8), intent(in) :: A(:), B(:)
    integer :: i

    is_equal = .true.
    !! logic for lexicographical comparison
    do i = 1, CISort_instance%numberOfSpecies 
      if ( A(i) == B(i) ) then
        is_equal = .true.
      else 
        is_equal = .false.
        return
      end if
    end do
    return

  end function CISort_is_equal
 
  !! Modified Recursive Quicksort routine
  subroutine CISort_mergeDuplicates( amplitudes, array, indices, arraySize )
    implicit none
    type(vector8), intent(inout) :: amplitudes
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: arraySize
    integer :: i
    integer :: spi, orb
    logical :: is_equal

    do i = 1, arraySize
      if ( amplitudes%values(i) == 0.0_8 ) cycle

      is_equal = .false.
      do orb = 1, CISort_instance%combinedNumberOfOrbitals 
        if ( array(orb,i) == array(orb,i+1) ) then
          is_equal = .true.
        else
          is_equal = .false.
          exit
        endif
      enddo

      if ( is_equal ) then

        !! merge in the next element 
        amplitudes%values(i+1) = amplitudes%values(i+1) + amplitudes%values(i)  
        !! clean the current element
        amplitudes%values(i) = 0.0_8
        indices%values(i) = 0_8
        array(:,i) = -1_1

      endif
    enddo

  end subroutine CISort_mergeDuplicates

  !! Reordering an array based on sorted indices in another array
  !! Algorithm taken from Zecong Hu
  !! https://stackoverflow.com/questions/60917343/

  subroutine CISort_sortArrayByIndex( array, auxindex_array, vectorSize  )
    implicit none
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: auxindex_array
    integer(8) :: vectorSize
    !real(8), allocatable :: temp(:)
    integer :: i, x, y

    CISort_instance%tmp_Vector(:) = -1_1

    do i = 1, vectorSize
      if ( auxindex_array%values(i) == -1 ) cycle 

        CISort_instance%tmp_Vector(:) = array(:,i)
        x = i
        y = auxindex_array%values(i) 

      do while ( y /= i )
        auxindex_array%values(x) = -1
        array(:,x) = array(:,y) 
        x = y
        y = auxindex_array%values(x)
      enddo 

      array(:,x) = CISort_instance%tmp_Vector(:) 
      auxindex_array%values(x) = -1

    enddo

  endsubroutine CISort_sortArrayByIndex

  subroutine CISort_sortVectorByIndex( targetVector, auxindex_array, vectorSize  )
    implicit none
    type(vector8), intent(inout) :: targetVector
    type(ivector8), intent(inout) :: auxindex_array
    integer(8) :: vectorSize
    real(8) :: temp
    integer :: i, x, y

    temp = 0.0_8

    do i = 1, vectorSize
      if ( auxindex_array%values(i) == -1 ) cycle 
      temp = targetVector%values(i) 
      x = i
      y = auxindex_array%values(i) 
      do while ( y /= i )
        auxindex_array%values(x) = -1
        targetVector%values(x) = targetVector%values(y)
        x = y
        y = auxindex_array%values(x)
      enddo 
      targetVector%values(x) = temp
      auxindex_array%values(x) = -1
    enddo

  endsubroutine CISort_sortVectorByIndex



end module CISort_
