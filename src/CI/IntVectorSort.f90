module IntVectorSort_
  use MolecularSystem_
  use CIcore_
  use Vector_
  use Matrix_
  implicit none

  type, public :: IntVectorSort
    integer(1), allocatable :: tmp_Vector (:) ! species
    integer(1), allocatable :: pivot (:) ! species, orbs
    integer :: numberOfSpecies
    integer :: combinedNumberOfOrbitals
  end type IntVectorSort

  type(IntVectorSort) :: IntVectorSort_instance

  !real(8), PARAMETER :: EPSILON = 1.0E-6

contains

  !! allocating global array, to avoid construting them in OMP
  subroutine IntVectorSort_constructor()
    implicit none
    integer :: numberOfSpecies, spi, m
   
    numberOfSpecies = CIcore_instance%numberOfSpecies 

    IntVectorSort_instance%numberOfSpecies = numberOfSpecies

    !! auxiliary variables to map orbitals from vector to array location
    IntVectorSort_instance%combinedNumberOfOrbitals = sum(CIcore_instance%numberOfOrbitals%values(:))
  
    allocate ( IntVectorSort_instance%tmp_Vector(  IntVectorSort_instance%combinedNumberOfOrbitals ) )
    allocate ( IntVectorSort_instance%pivot( IntVectorSort_instance%combinedNumberOfOrbitals ) )

    !do spi = 1, numberOfSpecies
      !call Vector_constructorInteger1(  IntVectorSort_instance%tmp_Vector(spi), int( CIcore_instance%numberOfOrbitals%values(spi), 8), -1_1  )
      !call Vector_constructorInteger1(  IntVectorSort_instance%pivot(spi), int( CIcore_instance%numberOfOrbitals%values(spi), 8), -1_1  )
    !enddo

    IntVectorSort_instance%pivot = -1_1

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
  subroutine IntVectorSort_swapVector(A, B)
    implicit none
    integer(1), intent(inout) :: A(:), B(:)

    !! Swap data
    IntVectorSort_instance%tmp_Vector(:) = A(:)
    A(:) = B(:)
    B(:) = IntVectorSort_instance%tmp_Vector(:)

  end subroutine IntVectorSort_swapVector

  !! Swap both data and index elements
  subroutine IntVectorSort_swapIndex(IndexA, IndexB)
    implicit none
    integer(8), intent(inout) :: IndexA, IndexB
    integer(8) :: tmpIndex

    !! Swap index
    tmpIndex = IndexA
    IndexA = IndexB
    IndexB = tmpIndex

  end subroutine IntVectorSort_swapIndex


  !! Modified Partition scheme to use swap_both
  subroutine IntVectorSort_partition(array, indices, low, high, pivot_index )
    implicit none
    integer(1), intent(inout) :: array(:,:) ! original array
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: low, high
    integer(8), intent(out) :: pivot_index
    integer(8) :: i, j
    integer :: spi, orb
    logical :: is_more    
    integer(8) :: mid_index

    ! select the midpoint as pivot to avoid worst case sceneario
    mid_index = ( low + high ) /2
    ! swap mid point to the end, really? why?
    call IntVectorSort_swapVector( array(:,mid_index), array(:,high) )
    call IntVectorSort_swapIndex( indices%values(mid_index), indices%values(high) )

    !! Select the last element (and its index) as the pivot
      IntVectorSort_instance%pivot(:) = array(:,high )

    i = low - 1 ! Index of the smaller element

    do j = low, high - 1
      is_more = .false.

      if ( array(1,j) == -1_1 ) cycle !! avoid zero, it means no configuration

       do orb = 1, IntVectorSort_instance%combinedNumberOfOrbitals 
          if ( array(orb,j) > IntVectorSort_instance%pivot(orb) ) then
            is_more = .true.
            exit
          else if ( array(orb,j) < IntVectorSort_instance%pivot(orb) ) then
            is_more = .false. 
            exit
          endif
       enddo
      
      if ( is_more ) then
        i = i + 1
        !! Swap both the vector and the index
        call IntVectorSort_swapVector( array(:,i), array(:,j) )
        call IntVectorSort_swapIndex( indices%values(i), indices%values(j) )
      end if
    end do
    
    !! Place the pivot in the correct sorted position (swap with arr(i+1) and index(i+1))
    call IntVectorSort_swapVector( array(:,i + 1), array(:,high) )
    call IntVectorSort_swapIndex( indices%values(i + 1), indices%values(high) )

    pivot_index = i + 1

  end subroutine IntVectorSort_partition

  !! Modified Recursive Quicksort routine
  recursive subroutine IntVectorSort_quicksort(array, indices, low, high)
    implicit none
    integer(1), intent(inout) :: array(:,:)
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
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: indices ! The index array
    integer(8), intent(in) :: arraySize
    integer :: i
    integer :: spi, orb
    logical :: is_equal

    do i = 1, arraySize
      if ( amplitudes%values(i) == 0.0_8 ) cycle

      is_equal = .false.
      !! check if the configuration is the same
      !species: do spi = 1, CIcore_instance%numberOfSpecies 
      !  do orb = 1, CIcore_instance%numberOfOrbitals%values(spi)
      !    if ( array(spi)%values(orb,i) == array(spi)%values(orb,i+1) ) then
      !      is_equal = .true.
      !    else
      !      is_equal = .false.
      !      exit species
      !    endif
      !  enddo
      !enddo species

       do orb = 1, IntVectorSort_instance%combinedNumberOfOrbitals 
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
        !do spi = 1, CIcore_instance%numberOfSpecies 
          array(:,i) = -1_1
        !enddo

      endif
    enddo

  end subroutine IntVectorSort_mergeDuplicates

  !! Reordering an array based on sorted indices in another array
  !! Algorithm taken from Zecong Hu
  !! https://stackoverflow.com/questions/60917343/

  subroutine IntVectorSort_sortArrayByIndex( array, auxindex_array, vectorSize  )
    implicit none
    integer(1), intent(inout) :: array(:,:)
    type(ivector8), intent(inout) :: auxindex_array
    integer(8) :: vectorSize
    !real(8), allocatable :: temp(:)
    integer :: numberOfSpecies, spi
    integer :: i, x, y

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !allocate ( temp ( numberOfSpecies ) )
    !do spi = 1, CIcore_instance%numberOfSpecies 
      IntVectorSort_instance%tmp_Vector(:) = -1_1
    !enddo
    !temp(:) = 0.0_8

    do i = 1, vectorSize
      if ( auxindex_array%values(i) == -1 ) cycle 
      !do spi = 1, CIcore_instance%numberOfSpecies 
        IntVectorSort_instance%tmp_Vector(:) = array(:,i)
        !IntVectorSort_instancetmp_Vector(:) = array%values(:,i) 
      !enddo

      x = i
      y = auxindex_array%values(i) 
      do while ( y /= i )
        auxindex_array%values(x) = -1
        !do spi = 1, CIcore_instance%numberOfSpecies 
          !array%values(:,x) = array%values(:,y)
          array(:,x) = array(:,y) 
        !enddo
        x = y
        y = auxindex_array%values(x)
      enddo 
      !array%values(:,x) = temp(:)
      !do spi = 1, CIcore_instance%numberOfSpecies 
        array(:,x) = IntVectorSort_instance%tmp_Vector(:) 
        !IntVectorSort_instancetmp_Vector(:) = array%values(:,i) 
      !enddo

      auxindex_array%values(x) = -1
    enddo

    !deallocate ( temp ) 

  endsubroutine IntVectorSort_sortArrayByIndex

  subroutine IntVectorSort_sortVectorByIndex( targetVector, auxindex_array, vectorSize  )
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

  endsubroutine IntVectorSort_sortVectorByIndex



end module IntVectorSort_
