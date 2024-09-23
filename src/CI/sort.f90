
module sort_
! Jorge Charry: Modified to evalute absolute value of the input Array as real
!
! This is a slightly-modified FORTRAN95 module containing the OpenMP recursive implementation of
! the quicksort algorithm by David Bal as available in: 
! https://bitbucket.org/daviddbal/multi-threaded-quicksort-mergesort-fortran/src/master/ 
! (Last access date: April 26th 2022)
!
! The only changes made are:
!
!      1) Usage of integer(kind=16) data type for the array A to be sorted
!      2) Inclusion of a companion array (indxA) 
!

use, intrinsic :: iso_fortran_env, only: int64, real64, int32

implicit none

public :: MTSort
private :: QuickSortA, QuickSortD, InsertionSortA, InsertionSortD, Merge8A, Merge8D

contains

! Main sorting subroutine - quicksort and insertion sort - multi-threaded
subroutine MTSort(A, indxA, nA, direction, nt)

! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32
use omp_lib

! DUMMY ARGUMENTS
! nA:  size of input array a
! a: input data to be sorted
! index array accompanying A: resorted as A
! direction: A = ascending, D = descending
! nt: number of threads, if omitted then uses all available
integer (kind=int64), intent(in) :: nA
integer (kind=int64), dimension(nA), intent(in out) :: indxA
real (kind=real64), dimension(nA), intent(in out) :: A
character, intent(in) :: direction
integer, optional, intent(in) :: nt

! LOCAL VARIABLES
! nt1: number of threads available.  Either a copy of input variable nt or omp_get_max_threads() if nt is not present.
! nt2: largest power of 2 number of threads (i.e. 2,4,8,16...)
! t: thread number
! s, f: Chunk indexes (start and finish) for each thread used by the initial quicksort/insertion sort combo.
! i, j: loop counters to merge sorted chunks together in order For example, with 4 threads first merge chunks 1 & 2 and 3 & 4 at the
!       same time.  Next merge the previously made 1+2 chunk and 3+4 chunk to generate the completely sorted result.
! levels: Log(base 2) of nt2 (i.e 8 threads has 3 levels 2**3).  Extent of outer merge loop (level_loop).
! step: number of threads available for each merge.
! gap: number of pieces between left and right chunks to be merged.
! span: number of pieces in each left and right chunk.  As level_loop progresses each left and right chunk spans multiple chunks.
! l1, l2, r1, r2: index locations of left and right chunks to be merged.  1 is start, 2 is finish
! left_part: temp array half the size of chunks being merged
! i_limit: array size limit when quicksort changes to insertion sort.  50 is a good value for small data types.  As derived
!               derived type increases in bytes the i_limit should become lower.
! verbose: T means to output messages, F means to output no messages
integer :: nt1
integer :: nt2
integer :: t=1
integer (kind=int64), dimension(:), allocatable :: s, f
integer :: i, j
integer :: levels
integer :: step
integer :: gap
integer :: span
integer (kind=int64) :: l1, l2, r1, r2
real (kind=real64), dimension(:), allocatable :: left_part    ! temp array for left half of array to be sorted
integer (kind=int64), dimension(:), allocatable :: indx_left_part    ! temp array for left half of array to be sorted
integer, parameter :: i_limit = 50
logical, parameter :: verbose = .false.
real(8) :: timeA, timeB

! ABSTRACT INTERFACE (used by procedure pointer)
abstract interface ! explicit interface block to make sure procedure pointers are given correct arguments

    subroutine QuickSort_Interface(A,indxA,nA,i_limit)
        use, intrinsic :: iso_fortran_env, only: int64, real64, int32
                integer (kind=int64), intent(in) :: nA
        real (kind=real64), dimension(nA), intent(in out) :: A
        integer (kind=int64), dimension(nA), intent(in out) :: indxA
        integer, intent(in) :: i_limit
    end subroutine QuickSort_Interface

    subroutine Merge_Interface(A,indxA,nA,B,indxB,nB,C,indxC,nC)
        use, intrinsic :: iso_fortran_env, only: int64, real64, int32
                integer (kind=int64), intent(in) :: nA, nB, nC
        real (kind=real64), dimension(nA), intent(in) :: A
        integer (kind=int64), dimension(nA), intent(in) :: indxA
        real (kind=real64), dimension(nB), intent(in) :: B
        integer (kind=int64), dimension(nB), intent(in) :: indxB
        real (kind=real64), dimension(nC), intent(out) :: C
        integer (kind=int64), dimension(nC), intent(out) :: indxC
    end subroutine Merge_Interface

end interface

! PROCEDURE POINTER
procedure (QuickSort_Interface), pointer :: QuickSort  ! which sort to use (ascending or descending)
procedure (Merge_Interface), pointer :: Merge8         ! which merge to use (ascending or descending)

!$  timeA = omp_get_wtime()

    ! POINT TO CORRECT SORT & MERGE PROCEDURES
    if (direction == "A" .or. direction == "a") then
        QuickSort => QuickSortA
        Merge8 => Merge8A
    else if (direction == "D" .or. direction == "d") then
        QuickSort => QuickSortD
        Merge8 => Merge8D
    else
        write (*,*) "ERROR: Invalid sort direction: ", direction
        write (*,*) "Valid options are A and D for ascending and descending sort"
        write (*,*) "Can not sort."
        return
    endif

    ! FIND AVAILABLE THREADS
    nt1 = 1  ! default to 1 thread in not using openmp
    !$ nt1 = omp_get_max_threads()  ! get max threads available if using OpenMP
    if (nt1 == 1) then
        if (verbose) then
            write (*,*) "WARNING: Multi-threaded sort requested, but either system has only one CPU core or OpenMP is not enabled."
        end if
    end if
    if (present(nt)) then
         nt1 = nt
    end if

    multithread: if (nA < 100000 .or. nt1 == 1) then

        ! Single-threaded
        if (verbose) write (*,*) "Single threaded"
        call QuickSort(A, indxA, nA, i_limit)

    else multithread

        ! PARALLEL MERGE SORT
        nt2 = 2 ** int(log(real(nt1))/log(2.0)) ! get largest power of 2 number of threads (i.e. 2,4,8,16...)
        if (verbose) then
            write (*,"(A,I3)") "Threads used =", nt1
            if (nt2 /= nt1) write (*,"(A,I3,a)") "Only efficiently using",nt2," threads."
        end if
        allocate (s(nt2),f(nt2))

        ! SORT PIECES
        !$omp parallel &
        !$omp num_threads(nt2) &
        !$omp private(t) &
        !$omp shared(a, s, f, nt2, nA)
            !$ t = omp_get_thread_num() + 1   ! add 1 to make first thread 1 instead of 0
            s(t) = nA * (t-1) / nt2 + 1       ! start
            f(t) = nA * t / nt2               ! finish
            call QuickSort(a(s(t):f(t)), &    ! section to be sorted
                           indxa(s(t):f(t)), &    
            &               f(t)-s(t)+1, &    ! size of section
            &               i_limit)          ! Insertion sort limit (50 is a good for small data elements)
        !$omp end parallel

        ! MERGE SORTED PIECES
        levels = int ( log(real(nt2))/log(2.0) ) 
        level_loop: do i = levels, 1, -1
            step = 2 ** (levels - i + 1)
            gap = 2 ** (levels - i)
            span = 2 ** (levels - i) - 1

            !$omp parallel &
            !$omp num_threads(nt2) &
            !$omp private (l1, l2, r1, r2, left_part, indx_left_part) &
            !$omp shared (a, indxa, s, f, nt2, gap, span, step)
            allocate (left_part(f(ceiling(real(step)/2.0))+1)) ! allocate left_part to max size of first half of pieces
            allocate (indx_left_part(f(ceiling(real(step)/2.0))+1)) 
            !$omp do
            merge_loop: do j = 1, nt2, step

                l1 = s(j)
                l2 = f(j+span)
                r1 = s(j+gap)
                r2 = f(j+gap+span)
                left_part(1:l2-l1+1) = a(l1:l2)
                indx_left_part(1:l2-l1+1) = indxa(l1:l2)
                call Merge8(left_part(1:l2-l1+1), &     ! left part
                            indx_left_part(1:l2-l1+1), &     ! left part
                &           l2-l1+1, &                  ! size of left part
                &           a(r1:r2), &                 ! right part
                &           indxa(r1:r2), &                 ! right part
                &           r2-r1+1, &                  ! size of right part
                &           a(l1:r2), &                 ! output array
                &           indxa(l1:r2), &                 ! output array
                &           r2-l1+1)                    ! size of output array
            enddo merge_loop
            !$omp end do
            deallocate (left_part)
            deallocate (indx_left_part)
            !$omp end parallel
        enddo level_loop

    endif multithread

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for sorting the vector : ", timeB - timeA ," (s)"

end subroutine MTSort

! Ascending Quicksort/Insertion Sort Combo for derived type
recursive subroutine QuickSortA(A,indxA,nA,limit)
! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64 , int32

! DUMMY ARGUMENTS
integer (kind=int64), intent(in) :: nA
real (kind=real64), dimension(nA), intent(in out) :: A
integer (kind=int64), dimension(nA), intent(in out) :: indxA
integer, intent(in) :: limit

! LOCAL VARIABLES
integer (kind=int64) :: left, right
real (kind=real64) :: pivot
real (kind=real64) :: temp
integer (kind=int64) :: marker ,i,j,k
integer (kind=int64) :: eighth
integer (kind=int64) :: temp_indx
real (kind=real64), dimension(9) :: sample

    if (nA > 1) then
        if (nA > limit) then ! Do quicksort for large groups
            ! ************************
            ! 9-SAMPLE PIVOT METHOD
            eighth = nA/8 
            j=0
            do i=1,nA,eighth-1
               j=j+1
               if ( j .gt. 9 ) exit
               sample(j)=abs(a(i))
            enddo
            ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
            if (sample(1) > sample(9)) then; temp = sample(1); sample(1) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(5)) then; temp = sample(1); sample(1) = sample(5); sample(5) = temp; end if
            if (sample(2) > sample(6)) then; temp = sample(2); sample(2) = sample(6); sample(6) = temp; end if
            if (sample(3) > sample(7)) then; temp = sample(3); sample(3) = sample(7); sample(7) = temp; end if
            if (sample(4) > sample(8)) then; temp = sample(4); sample(4) = sample(8); sample(8) = temp; end if
            if (sample(5) > sample(9)) then; temp = sample(5); sample(5) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(3)) then; temp = sample(1); sample(1) = sample(3); sample(3) = temp; end if
            if (sample(2) > sample(4)) then; temp = sample(2); sample(2) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(7)) then; temp = sample(5); sample(5) = sample(7); sample(7) = temp; end if
            if (sample(6) > sample(8)) then; temp = sample(6); sample(6) = sample(8); sample(8) = temp; end if
            if (sample(3) > sample(9)) then; temp = sample(3); sample(3) = sample(9); sample(9) = temp; end if
            if (sample(3) > sample(5)) then; temp = sample(3); sample(3) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(6)) then; temp = sample(4); sample(4) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(9)) then; temp = sample(7); sample(7) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(2)) then; temp = sample(1); sample(1) = sample(2); sample(2) = temp; end if
            if (sample(3) > sample(4)) then; temp = sample(3); sample(3) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(6)) then; temp = sample(5); sample(5) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(8)) then; temp = sample(7); sample(7) = sample(8); sample(8) = temp; end if
            if (sample(2) > sample(9)) then; temp = sample(2); sample(2) = sample(9); sample(9) = temp; end if
            if (sample(2) > sample(5)) then; temp = sample(2); sample(2) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(7)) then; temp = sample(4); sample(4) = sample(7); sample(7) = temp; end if
           !if (sample(6) > sample(9)) then; temp = sample(6); sample(6) = sample(9); sample(9) = temp; end if ! skipped
           !if (sample(2) > sample(3)) then; temp = sample(2); sample(2) = sample(3); sample(3) = temp; end if ! skipped
            if (sample(4) > sample(5)) then; temp = sample(4); sample(4) = sample(5); sample(5) = temp; end if
           !if (sample(6) > sample(7)) then; temp = sample(6); sample(6) = sample(7); sample(7) = temp; end if ! skipped
           !if (sample(8) > sample(9)) then; temp = sample(8); sample(8) = sample(9); sample(9) = temp; end if ! skipped
            pivot = sample(5)
            ! ************************
            left = 0
            right = nA + 1
            do while (left < right)
                right = right - 1
                do while (abs(A(right)) > pivot)
                    right = right - 1
                end do
                left = left + 1
                do while (abs(A(left)) < pivot)
                    left = left + 1
                end do
                if (left < right) then
                    temp = A(left)
                    A(left) = A(right)
                    A(right) = temp
                    temp_indx = indxA(left)
                    indxA(left) =indxA(right)
                    indxA(right) = temp_indx
                end if
            end do

            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if

            call QuickSortA(A(:marker-1),indxA(:marker-1),marker-1,limit)
            call QuickSortA(A(marker:),indxA(marker:),nA-marker+1,limit)

        else
            call InsertionSortA(A,indxA,nA)    ! Insertion sort for small groups is faster than Quicksort
        end if
    end if

end subroutine QuickSortA

subroutine InsertionSortA(A,indxA,nA)
! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32

! DUMMY ARGUMENTS
integer (kind=int64), intent(in) :: nA
real (kind=real64), dimension(nA), intent(in out) :: A
integer (kind=int64), dimension(nA), intent(in out) :: indxA

! LOCAL VARIABLES
real (kind=real64) :: temp
integer (kind=int64) :: temp_indx
integer (kind=int64) :: i, j

    outter: do i = 2, nA
        j = i - 1
        temp = A(i)
        temp_indx = indxA(i)
        inner: do
            if (j == 0) exit inner
            if (abs(a(j)) <= abs(temp)) exit inner
            A(j+1) = A(j)
            indxA(j+1) = indxA(j)
            j = j - 1
        end do inner
        a(j+1) = temp
        indxA(j+1) = temp_indx
    end do outter

end subroutine InsertionSortA

! Descending Quicksort/Insertion Sort Combo for derived type
recursive subroutine QuickSortD(A,indxA,nA,limit)
! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32

! DUMMY ARGUMENTS
integer (kind=int64), intent(in) :: nA
real (kind=real64), dimension(nA), intent(in out) :: A
integer (kind=int64), dimension(nA), intent(in out) :: indxA
integer, intent(in) :: limit

! LOCAL VARIABLES
integer (kind=int64) :: left, right
real (kind=real64) :: pivot
real (kind=real64) :: temp
integer (kind=int64) :: temp_indx
integer (kind=int64) :: marker
integer (kind=int64) :: eighth,i,j
real (kind=real64), dimension(9) :: sample

    if (nA > 1) then
        if (nA > limit) then ! Do quicksort for large groups
            ! ************************
            ! 9-SAMPLE PIVOT METHOD
            eighth = nA/8
!           sample = a(1:nA:nine)
            j=0
            do i=1,nA,eighth-1
               j=j+1
               if ( j .gt. 9 ) exit
               sample(j)=abs(a(i))
            enddo
            ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
            if (sample(1) > sample(9)) then; temp = sample(1); sample(1) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(5)) then; temp = sample(1); sample(1) = sample(5); sample(5) = temp; end if
            if (sample(2) > sample(6)) then; temp = sample(2); sample(2) = sample(6); sample(6) = temp; end if
            if (sample(3) > sample(7)) then; temp = sample(3); sample(3) = sample(7); sample(7) = temp; end if
            if (sample(4) > sample(8)) then; temp = sample(4); sample(4) = sample(8); sample(8) = temp; end if
            if (sample(5) > sample(9)) then; temp = sample(5); sample(5) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(3)) then; temp = sample(1); sample(1) = sample(3); sample(3) = temp; end if
            if (sample(2) > sample(4)) then; temp = sample(2); sample(2) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(7)) then; temp = sample(5); sample(5) = sample(7); sample(7) = temp; end if
            if (sample(6) > sample(8)) then; temp = sample(6); sample(6) = sample(8); sample(8) = temp; end if
            if (sample(3) > sample(9)) then; temp = sample(3); sample(3) = sample(9); sample(9) = temp; end if
            if (sample(3) > sample(5)) then; temp = sample(3); sample(3) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(6)) then; temp = sample(4); sample(4) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(9)) then; temp = sample(7); sample(7) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(2)) then; temp = sample(1); sample(1) = sample(2); sample(2) = temp; end if
            if (sample(3) > sample(4)) then; temp = sample(3); sample(3) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(6)) then; temp = sample(5); sample(5) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(8)) then; temp = sample(7); sample(7) = sample(8); sample(8) = temp; end if
            if (sample(2) > sample(9)) then; temp = sample(2); sample(2) = sample(9); sample(9) = temp; end if
            if (sample(2) > sample(5)) then; temp = sample(2); sample(2) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(7)) then; temp = sample(4); sample(4) = sample(7); sample(7) = temp; end if
           !if (sample(6) > sample(9)) then; temp = sample(6); sample(6) = sample(9); sample(9) = temp; end if ! skipped
           !if (sample(2) > sample(3)) then; temp = sample(2); sample(2) = sample(3); sample(3) = temp; end if ! skipped
            if (sample(4) > sample(5)) then; temp = sample(4); sample(4) = sample(5); sample(5) = temp; end if
           !if (sample(6) > sample(7)) then; temp = sample(6); sample(6) = sample(7); sample(7) = temp; end if ! skipped
           !if (sample(8) > sample(9)) then; temp = sample(8); sample(8) = sample(9); sample(9) = temp; end if ! skipped
            pivot = sample(5)
            ! ************************

            left = 0
            right = nA + 1
            do while (left < right)
                right = right - 1
                do while (abs(A(right)) < pivot)
                    right = right - 1
                end do
                left = left + 1
                do while (abs(A(left)) > pivot)
                    left = left + 1
                end do
                if (left < right) then
                    temp = A(left)
                    A(left) = A(right)
                    A(right) = temp
                    temp_indx = indxA(left)
                    indxA(left) = indxA(right)
                    indxA(right) = temp_indx
                end if
            end do

            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if

            call QuickSortD(A(:marker-1),indxA(:marker-1),marker-1,limit)
            call QuickSortD(A(marker:),indxA(marker:),nA-marker+1,limit)

        else
            call InsertionSortD(A,indxA,nA)    ! Insertion sort for small groups is faster than Quicksort
        end if
    end if

end subroutine QuickSortD

subroutine InsertionSortD(A,indxA,nA)
! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32

! DUMMY ARGUMENTS
integer (kind=int64), intent(in) :: nA
real (kind=real64), dimension(nA), intent(in out) :: A
integer (kind=int64), dimension(nA), intent(in out) :: indxA

! LOCAL VARIABLES
real (kind=real64) :: temp
integer (kind=int64) :: temp_indx
integer (kind=int64) :: i, j

    outter: do i = 2, nA
        j = i - 1
        temp = A(i)
        temp_indx = indxA(i)
        inner: do
            if (j == 0) exit inner
            if (abs(a(j)) >= abs(temp)) exit inner
            A(j+1) = A(j)
            indxA(j+1) = indxA(j)
            j = j - 1
        end do inner
        A(j+1) = temp
        indxA(j+1) = temp_indx
    end do outter

end subroutine InsertionSortD

! Ascending merge (merges 2 ascending sorted lists into 1 ascending sorted list)
subroutine Merge8a(A, indxA, nA, B, indxB, nB, C, indxC, nC)

! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32

! DUMMY ARGUMENTS
integer(kind=int64), intent(in) :: nA, nB, nC         ! Size of arrays.  Normal usage: nC = nA+nB
! ***NOTE: USING EXPLICIT SHAPE DUMMY ARRAYS AS SHOW BELOW SOMETIMES CAUSES SEGMENTATION FAULTS WITH LARGE ARRAY SIZES
! Explicit shape arrays are faster than allocatable arrays, but if problems occur first try to increase the stack memory size.
! For Linux in terminal enter "ulimit -s unlimited" to set stack size to unlimited.  If problems still occur then switch to above
! allocatable variable declarations and uncomment allocate statement below.
real (kind=real64), dimension(nA), intent(in) :: A  ! left part
real (kind=real64), dimension(nB), intent(in) :: B  ! right part
real (kind=real64), dimension(nC), intent(out) :: C ! output array
integer (kind=int64), dimension(nA), intent(in) :: indxA  ! left part
integer (kind=int64), dimension(nB), intent(in) :: indxB  ! right part
integer (kind=int64), dimension(nC), intent(out) :: indxC ! output array
!   integer (kind=16), dimension(:), intent(in) :: A
!   integer (kind=16), dimension(:), intent(in) :: B
!   integer (kind=16), dimension(:), intent(out) :: C
! Note: Under normal usage array A is a copy of the left part of array A.  B is the right part of array C (nA+1:nC) meaning B and C
! share the same memory.  Array A must be a copy to insure the merge doesn't write over, but B doesn't need to be a copy.

! LOCAL VARIABLES
integer(kind=int64) :: i, j, k

    i = 1; j = 1; k = 1
    do
        if (i > nA .or. j > nB) exit
        if (abs(a(i)) <= abs(b(j))) then
            c(k) = a(i)
            indxc(k) = indxa(i)
            i = i + 1
        else
            c(k) = b(j)
            indxc(k) = indxb(j)
            j = j + 1
        endif
        k = k + 1
    enddo
    if (i <= nA) then
        c(k:) = a(i:)
        indxc(k:) = indxa(i:)
        return
    endif
    if (j <= nB) then 
       c(k:) = b(j:)    ! This statement is only necessary for multi-threaded merge
       indxc(k:) = indxb(j:)  
    endif

end subroutine Merge8a

! Descending merge (merges 2 descending sorted lists into 1 descending sorted list)
subroutine Merge8d(A, indxA, nA, B, indxB, nB, C, indxC, nC)

! USED MODULES
use, intrinsic :: iso_fortran_env, only: int64, real64, int32

! DUMMY ARGUMENTS
integer(kind=int64), intent(in) :: nA, nB, nC         ! Size of arrays.  Normal usage: nC = nA+nB
! ***NOTE: USING EXPLICIT SHAPE DUMMY ARRAYS AS SHOW BELOW SOMETIMES CAUSES SEGMENTATION FAULTS WITH LARGE ARRAY SIZES
! Explicit shape arrays are faster than allocatable arrays, but if problems occur first try to increase the stack memory size.
! For Linux in terminal enter "ulimit -s unlimited" to set stack size to unlimited.  If problems still occur then switch to above
! allocatable variable declarations and uncomment allocate statement below.
real (kind=real64), dimension(nA), intent(in) :: A  ! left part
real (kind=real64), dimension(nB), intent(in) :: B  ! right part
real (kind=real64), dimension(nC), intent(out) :: C ! output array
integer (kind=int64), dimension(nA), intent(in) :: indxA  ! left part
integer (kind=int64), dimension(nB), intent(in) :: indxB  ! right part
integer (kind=int64), dimension(nC), intent(out) :: indxC ! output array
!   integer (kind=16), dimension(:), intent(in) :: A
!   integer (kind=16), dimension(:), intent(in) :: B
!   integer (kind=16), dimension(:), intent(out) :: C
! Note: Under normal usage array A is a copy of the left part of array A.  B is the right part of array C (nA+1:nC) meaning B and C
! share the same memory.  Array A must be a copy to insure the merge doesn't write over, but B doesn't need to be a copy.

! LOCAL VARIABLES
integer(kind=int64) :: i, j, k

    i = 1; j = 1; k = 1
    do
        if (i > nA .or. j > nB) exit
        if (abs(a(i)) >= abs(b(j))) then
            c(k) = a(i)
            indxc(k) = indxa(i)
            i = i + 1
        else
            c(k) = b(j)
            indxc(k) = indxb(j)
            j = j + 1
        endif
        k = k + 1

    enddo
    if (i <= nA) then
        c(k:) = a(i:)
        indxc(k:) = indxa(i:)
        return
    endif
    if (j <= nB)  then
         c(k:) = b(j:)    ! THIS STATEMENT IS NOT NECESSARY IN SINGLE THREADED MERGE
         indxc(k:) = indxb(j:)    
    endif

end subroutine Merge8d

end module sort_

! END OF MODULE SORT
