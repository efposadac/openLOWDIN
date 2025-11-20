module CISort_
!! Author: Jorge Charry
!! Collection of sorting functions adapted for SCI program.
!! Inspired from CENCALC code https://github.com/dimassuarez/cencalc_quicksort/ 
  use CIcore_
  use Vector_
  use Matrix_
  implicit none

  type, public :: CISort
    integer(1), allocatable :: tmp_Vector (:,:) ! species, n threads
    integer(1), allocatable :: pivot (:,:) ! species, orbs
    integer :: numberOfSpecies
    integer :: combinedNumberOfOrbitals
  end type CISort

  type(CISort) :: CISort_instance

contains

  !! allocating global array, to avoid construting them in OMP
  subroutine CISort_constructor()
    implicit none
    integer :: numberOfSpecies, spi, m
    integer :: nproc, n
   
    numberOfSpecies = CIcore_instance%numberOfSpecies 

    CISort_instance%numberOfSpecies = numberOfSpecies
    nproc = CIcore_instance%nproc 

    !! auxiliary variables to map orbitals from vector to array location
    CISort_instance%combinedNumberOfOrbitals = sum(CIcore_instance%numberOfOrbitals%values(:))
  
    allocate ( CISort_instance%tmp_Vector(  CISort_instance%combinedNumberOfOrbitals, nproc + 1 ) )
    allocate ( CISort_instance%pivot( CISort_instance%combinedNumberOfOrbitals, nproc + 1 ) )

    CISort_instance%tmp_Vector = 0_1
    CISort_instance%pivot = -1_1

  end subroutine CISort_constructor

  !! allocating global array, to avoid construting them in OMP
  subroutine CISort_destructor()
    implicit none
   
    deallocate ( CISort_instance%pivot )
    deallocate ( CISort_instance%tmp_Vector )

  end subroutine CISort_destructor

  !! Modified Recursive Quicksort routine
  recursive subroutine CISort_quicksort_matrix(matrix, indices, left_end, right_end, n)
    implicit none
    integer(1), intent(inout) :: matrix(:,:)
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer, intent(in) :: n
    integer(8) :: pivot_index
    integer, parameter :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) THEN
      !! Use interchange sort for small lists
      call CISort_interchange_sort_matrix(matrix, indices, left_end, right_end, n)

    else

      if ( left_end < right_end ) then
          !! Partition using the indexed version
          call CISort_partition_hoare_matrix(matrix, indices, left_end, right_end, pivot_index, n)
  
          !! Recursively sort the sub-arrays
          call CISort_quicksort_matrix(matrix, indices, left_end, pivot_index, n)
          call CISort_quicksort_matrix(matrix, indices, pivot_index + 1, right_end, n )
      end if
    endif

  end subroutine CISort_quicksort_matrix

  !! Modified Partition scheme to use swap_both
  subroutine CISort_partition_hoare_matrix(matrix, indices, left_end, right_end, pivot_index, n )
    implicit none
    integer(1), intent(inout) :: matrix(:,:) ! original array
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8), intent(out) :: pivot_index
    integer, intent(in) :: n
    integer(8) :: i, j
    integer :: spi, orb
    logical :: is_more    
    integer(8) :: mid_index

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( left_end + right_end ) /2
    CISort_instance%pivot(:,n) = matrix(:, mid_index )

    i = left_end - 1 ! Index of the smaller element
    j = right_end + 1 ! Index of the higher element

    do 
      do 
        i = i + 1
        is_more = .false.
        do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( CISort_instance%pivot(orb,n) > matrix(orb,i) ) then
            is_more = .true.
            exit
          else if ( CISort_instance%pivot(orb,n) < matrix(orb,i) ) then
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
          if ( CISort_instance%pivot(orb,n) < matrix(orb,j) ) then
            is_more = .true.
            exit
          else if ( CISort_instance%pivot(orb,n) > matrix(orb,j) ) then
            is_more = .false. 
            exit
          endif
        enddo
        if ( .not. is_more ) exit
      enddo

      if ( i >= j ) exit 

      !! Swap both the vector and the index
      call CISort_swapVector( matrix(:,i), matrix(:,j), n )
      call CISort_swapIndex( indices(i), indices(j) )
    enddo
    pivot_index = j 

  end subroutine CISort_partition_hoare_matrix
 
  !! Lomuto partition scheme
  subroutine CISort_partition_lomuto(matrix, indices, left_end, right_end, pivot_index, n )
    implicit none
    integer(1), intent(inout) :: matrix(:,:) ! original array
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8), intent(out) :: pivot_index
    integer, intent(in) :: n
    integer(8) :: i, j
    integer :: spi, orb
    logical :: is_more    
    integer(8) :: mid_index

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( left_end + right_end ) /2
    ! swap mid point to the end, really? why?
    call CISort_swapVector( matrix(:,mid_index), matrix(:,right_end), n )
    call CISort_swapIndex( indices(mid_index), indices(right_end ) )

    !! Select the last element (and its index) as the pivot
    CISort_instance%pivot(:,n) = matrix(:,right_end )

    i = left_end - 1 ! Index of the smaller element

    do j = left_end, right_end - 1
      is_more = .false.

      if ( matrix(1,j) == -1_1 ) cycle !! avoid zero, it means no configuration

       do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( matrix(orb,j) > CISort_instance%pivot(orb,n) ) then
            is_more = .true.
            exit
          else if ( matrix(orb,j) < CISort_instance%pivot(orb,n) ) then
            is_more = .false. 
            exit
          endif
       enddo
      
      if ( is_more ) then
        i = i + 1
        !! Swap both the vector and the index
        call CISort_swapVector( matrix(:,i), matrix(:,j), n )
        call CISort_swapIndex( indices(i), indices(j) )
      end if
    end do
    
    !! Place the pivot in the correct sorted position (swap with arr(i+1) and index(i+1))
    call CISort_swapVector( matrix(:,i + 1), matrix(:,right_end), n )
    call CISort_swapIndex( indices(i + 1), indices(right_end) )

    pivot_index = i + 1

  end subroutine CISort_partition_lomuto

  !! Sorting for small arrays, no more pivot
  subroutine CISort_interchange_sort_matrix(matrix, indices, left_end, right_end, n) 
    implicit none
    integer(1), intent(inout) :: matrix(:,:)
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer, intent(in) :: n
    integer :: i, j
    logical :: is_more    
    integer :: orb

    do i = left_end, right_end - 1
      !if ( matrix(1,i) == -1_1 ) cycle !! avoid zero, it means no configuration
      do j = i+1, right_end
        is_more = .false.

        do orb = 1, CISort_instance%combinedNumberOfOrbitals 
          if ( matrix(orb,j) < matrix(orb,i) ) then
            is_more = .true.
            exit
          else if ( matrix(orb,j) > matrix(orb,i) ) then
            is_more = .false. 
            exit
          endif
       enddo
        
        if ( is_more ) then
          !! Swap both the vector and the index
          call CISort_swapVector( matrix(:,i), matrix(:,j), n )
          call CISort_swapIndex( indices(i), indices(j) )
        end if

      enddo 
    enddo 

  end subroutine CISort_interchange_sort_matrix

  !! Swap both data and index elements
  subroutine CISort_swapVector(A, B, n)
    implicit none
    integer(1), intent(inout) :: A(:), B(:)
    integer, intent(in) :: n

    !! Swap data
    CISort_instance%tmp_Vector(:,n) = A(:)
    A(:) = B(:)
    B(:) = CISort_instance%tmp_Vector(:,n )

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
  recursive subroutine CISort_quicksort_vector(vector, indices, left_end, right_end)
    implicit none
    real(8), intent(inout) :: vector(:)
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8) :: pivot_index
    integer, parameter :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) THEN
      !! Use interchange sort for small lists
      call CISort_interchange_sort_vector(vector, indices, left_end, right_end)

    else

      if ( left_end < right_end ) then
          !! Partition using the indexed version
          call CISort_partition_hoare_vector(vector, indices, left_end, right_end, pivot_index )
  
          !! Recursively sort the sub-arrays
          call CISort_quicksort_vector(vector, indices, left_end, pivot_index )
          call CISort_quicksort_vector(vector, indices, pivot_index + 1, right_end )
      end if
    endif

  end subroutine CISort_quicksort_vector

  !! Modified Partition scheme to use swap_both
  subroutine CISort_partition_hoare_vector(vector, indices, left_end, right_end, pivot_index )
    implicit none
    real(8), intent(inout) :: vector(:) ! original array
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer(8), intent(out) :: pivot_index
    integer(8) :: i, j
    logical :: is_less
    integer(8) :: mid_index
    real(8) :: pivot

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( left_end + right_end ) /2
    pivot = vector( mid_index )

    i = left_end - 1 ! Index of the smaller element
    j = right_end + 1 ! Index of the higher element

    do 
      do 
        i = i + 1
        is_less = .false.
        if ( abs(pivot) < abs(vector(i)) ) then
          is_less = .true.
        else if ( abs(pivot) > abs(vector(i)) ) then
          is_less = .false. 
        endif
        if ( .not. is_less ) exit
      enddo

      do 
        j = j - 1
        is_less = .false.
        if ( abs(pivot) > abs(vector(j)) ) then
          is_less = .true.
        else if ( abs(pivot) < abs(vector(j)) ) then
          is_less = .false. 
        endif
        if ( .not. is_less ) exit
      enddo

      if ( i >= j ) exit 

      !! Swap both the element and the index
      call CISort_swapElement( vector(i), vector(j) )
      call CISort_swapIndex( indices(i), indices(j) )
    enddo
    pivot_index = j 

  end subroutine CISort_partition_hoare_vector

  !! Sorting for small arrays, no more pivot
  subroutine CISort_interchange_sort_vector(vector, indices, left_end, right_end) 
    implicit none
    real(8), intent(inout) :: vector(:)
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: left_end, right_end
    integer :: i, j
    logical :: is_less

    do i = left_end, right_end - 1
      !if ( vector(i) == 0.0_8 ) cycle !! avoid zero, it means no configuration
      do j = i+1, right_end
        is_less = .false.

        if ( abs(vector(j)) > abs(vector(i)) ) then
          is_less = .true.
        else if ( abs(vector(j)) < abs(vector(i)) ) then
          is_less = .false. 
        endif
        
        if ( is_less ) then
          !! Swap both the element and the index
          call CISort_swapElement( vector(i), vector(j) )
          call CISort_swapIndex( indices(i), indices(j) )
        end if

      enddo 
    enddo 

  end subroutine CISort_interchange_sort_vector

  !! Swap both data and index elements. I know it's an overkill, but just to keep the structure
  subroutine CISort_swapElement(A, B)
    implicit none
    real(8), intent(inout) :: A, B
    real(8) :: tmp_element

    !! Swap data
    tmp_element = A
    A = B
    B = tmp_element

  end subroutine CISort_swapElement
 
  !! Modified Recursive Quicksort routine
  subroutine CISort_mergeDuplicates( amplitudes, matrix, indices, arraySize )
    implicit none
    real(8), intent(inout) :: amplitudes(:)
    integer(1), intent(inout) :: matrix(:,:)
    integer(8), intent(inout) :: indices(:) ! The index array
    integer(8), intent(in) :: arraySize
    integer :: i
    integer :: spi, orb
    logical :: is_equal

    do i = 1, arraySize -1 
      if ( amplitudes(i) == 0.0_8 ) cycle

      is_equal = .false.
      do orb = 1, CISort_instance%combinedNumberOfOrbitals 
        if ( matrix(orb,i) == matrix(orb,i+1) ) then
          is_equal = .true.
        else
          is_equal = .false.
          exit
        endif
      enddo

      if ( is_equal ) then

        !! merge in the next element 
        amplitudes(i+1) = amplitudes(i+1) + amplitudes(i)  
        !! clean the current element
        amplitudes(i) = 0.0_8
        indices(i) = 0_8
        matrix(:,i) = -1_1

      endif
    enddo

  end subroutine CISort_mergeDuplicates

  !! Reordering an array based on sorted indices in another array
  !! Algorithm taken from Zecong Hu
  !! https://stackoverflow.com/questions/60917343/

  subroutine CISort_sortArrayByIndex( matrix, auxindex_array, vectorSize, n  )
    implicit none
    integer(1), intent(inout) :: matrix(:,:)
    integer(8), intent(inout) :: auxindex_array(:)
    integer(8) :: vectorSize
    integer, intent(in) :: n
    !real(8), allocatable :: temp(:)
    integer :: i, x, y

    CISort_instance%tmp_Vector(:,n) = -1_1

    do i = 1, vectorSize
      if ( auxindex_array(i) == -1 ) cycle 

        CISort_instance%tmp_Vector(:,n) = matrix(:,i)
        x = i
        y = auxindex_array(i) 

      do while ( y /= i )
        auxindex_array(x) = -1
        matrix(:,x) = matrix(:,y) 
        x = y
        y = auxindex_array(x)
      enddo 

      matrix(:,x) = CISort_instance%tmp_Vector(:,n) 
      auxindex_array(x) = -1

    enddo

  endsubroutine CISort_sortArrayByIndex

  subroutine CISort_sortVectorByIndex( targetVector, auxindex_array, vectorSize  )
    implicit none
    real(8), intent(inout) :: targetVector(:)
    integer(8), intent(inout) :: auxindex_array(:)
    integer(8) :: vectorSize
    real(8) :: temp
    integer :: i, x, y

    temp = 0.0_8

    do i = 1, vectorSize
      if ( auxindex_array(i) == -1 ) cycle 
      temp = targetVector(i) 
      x = i
      y = auxindex_array(i) 
      do while ( y /= i )
        auxindex_array(x) = -1
        targetVector(x) = targetVector(y)
        x = y
        y = auxindex_array(x)
      enddo 
      targetVector(x) = temp
      auxindex_array(x) = -1
    enddo

  endsubroutine CISort_sortVectorByIndex




end module CISort_
