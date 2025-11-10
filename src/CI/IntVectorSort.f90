module IntVectorSort_
  use MolecularSystem_
  implicit none

  type, public :: IntVectorSort
    real(8), allocatable :: tmp_Vector (:) ! species
    real(8), allocatable :: pivot (:) ! species
    integer :: numberOfSpecies
  end type IntVectorSort

  type(IntVectorSort) :: IntVectorSort_instance

  !real(8), PARAMETER :: EPSILON = 1.0E-6

contains

  !! allocating global array, to avoid construting them in OMP
  subroutine IntVectorSort_constructor()
    implicit none
    integer :: numberOfSpecies
   
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    IntVectorSort_instance%numberOfSpecies = numberOfSpecies

    allocate ( IntVectorSort_instance%tmp_Vector( numberOfSpecies) )
    allocate ( IntVectorSort_instance%pivot( numberOfSpecies) )

    IntVectorSort_instance%tmp_Vector = 0.0_8
    IntVectorSort_instance%pivot = 0.0_8

  end subroutine IntVectorSort_constructor

  !! allocating global array, to avoid construting them in OMP
  subroutine IntVectorSort_destructor()
    implicit none
   
    deallocate ( IntVectorSort_instance%pivot )
    deallocate ( IntVectorSort_instance%tmp_Vector )

  end subroutine IntVectorSort_destructor

  !! (The 'is_more' function for vector_t remains unchanged)
  function IntVectorSort_is_more(A, B ) result (is_more)
    implicit none
    logical :: is_more
    real(8), intent(in) :: A(:), B(:)
    integer :: i

    !! logic for lexicographical comparison
    do i = 1, IntVectorSort_instance%numberOfSpecies 
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

  end function IntVectorSort_is_more

  !! (The 'is_more' function for vector_t remains unchanged)
  function IntVectorSort_is_equal(A, B ) result (is_equal)
    implicit none
    logical :: is_equal
    real(8), intent(in) :: A(:), B(:)
    integer :: i

    is_equal = .true.
    !! logic for lexicographical comparison
    do i = 1, IntVectorSort_instance%numberOfSpecies 
      if ( A(i) == B(i) ) then
        is_equal = .true.
      else 
        is_equal = .false.
        return
      end if
    end do
    return

  end function IntVectorSort_is_equal
  
 
  
  !! Swap both data and index elements
  subroutine IntVectorSort_swap(A, B, IndexA, IndexB)
    implicit none
    real(8), intent(inout) :: A(:), B(:)
    integer(8), intent(inout) :: IndexA, IndexB
    integer(8) :: tmpIndex

    !! Swap data
    IntVectorSort_instance%tmp_Vector(:) = A(:)
    A(:) = B(:)
    B(:) = IntVectorSort_instance%tmp_Vector(:)

    !! Swap index
    tmpIndex = IndexA
    IndexA = IndexB
    IndexB = tmpIndex

  end subroutine IntVectorSort_swap


  !! Modified Partition scheme to use swap_both
  subroutine IntVectorSort_partition(array, indices, low, high, pivot_index )
    implicit none
    type(Matrix), intent(inout) :: array ! original array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: low, high
    integer(8), intent(out) :: pivot_index
    integer(8) :: i, j
    
    !! Select the last element (and its index) as the pivot
    IntVectorSort_instance%pivot(:) = array%values(:,high)

    i = low - 1 ! Index of the smaller element

    do j = low, high - 1
      !! Use the custom comparison function (arr(j) < pivot)
      if ( array%values(1,j) == 0.0_8 ) cycle !! avoid zero, it means no configuration
      if ( IntVectorSort_is_more(array%values(:,j), IntVectorSort_instance%pivot(:) ) ) then
        i = i + 1
        !! Swap both the vector and the index
        call  IntVectorSort_swap(array%values(:,i), array%values(:,j), indices%values(i), indices%values(j) )
      end if
    end do
    
    !! Place the pivot in the correct sorted position (swap with arr(i+1) and index(i+1))
    call IntVectorSort_swap(array%values(:,i + 1), array%values(:,high), indices%values(i+1), indices%values(high) )
    pivot_index = i + 1

  end subroutine IntVectorSort_partition

  !! Modified Recursive Quicksort routine
  recursive subroutine IntVectorSort_quicksort(array, indices, low, high)
    implicit none
    type(matrix), intent(inout) :: array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: low, high
    integer(8) :: pivot

    !print *, low, high
    if (low < high) then
        !! Partition using the indexed version
        call IntVectorSort_partition(array, indices, low, high, pivot)

        !! Recursively sort the sub-arrays
        call IntVectorSort_quicksort(array, indices, low, pivot - 1)
        call IntVectorSort_quicksort(array, indices, pivot + 1, high)
    end if
  end subroutine IntVectorSort_quicksort

  !! Modified Recursive Quicksort routine
  subroutine IntVectorSort_mergeDuplicates( amplitudes, array, indices, arraySize )
    implicit none
    type(vector8), intent(inout) :: amplitudes
    type(matrix), intent(inout) :: array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: arraySize
    integer :: i

    do i = 1, arraySize
      if ( amplitudes%values(i) == 0.0_8 ) cycle

      !! check if the configuration is the same
      if ( IntVectorSort_is_equal( array%values(:,i),  array%values(:,i+1) ) )  then
        !! merge in the next element 
        amplitudes%values(i+1) = amplitudes%values(i+1) + amplitudes%values(i)  
        !! clean the current element
        amplitudes%values(i) = 0.0_8
        array%values(:,i) = 0.0_8
        indices%values(i) = 0.0_8
      endif
    enddo

  end subroutine IntVectorSort_mergeDuplicates



end module IntVectorSort_
