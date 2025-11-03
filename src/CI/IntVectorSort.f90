MODULE Vector_Sort_Mod
    IMPLICIT NONE
    ! ... (VECTOR_LEN and EPSILON definitions remain here)
    INTEGER, PARAMETER :: VECTOR_LEN = 3 
    REAL, PARAMETER :: EPSILON = 1.0E-6

    TYPE :: vector_t
        REAL, DIMENSION(VECTOR_LEN) :: components
    END TYPE vector_t

CONTAINS
    ! (The 'is_less' function for vector_t remains unchanged)
    PURE LOGICAL FUNCTION is_less(A, B)
        TYPE(vector_t), INTENT(IN) :: A, B
        INTEGER :: i
        ! ... (Logic for lexicographical comparison remains here)
        DO i = 1, VECTOR_LEN
            IF (A%components(i) < B%components(i) - EPSILON) THEN
                is_less = .TRUE.; RETURN
            ELSE IF (A%components(i) > B%components(i) + EPSON) THEN
                is_less = .FALSE.; RETURN
            END IF
        END DO
        is_less = .FALSE.; RETURN
    END FUNCTION is_less
    
---

    ! Modified SUBROUTINE to swap both data and index elements
    SUBROUTINE swap_both(A, B, I_A, I_B)
        TYPE(vector_t), INTENT(INOUT) :: A, B
        INTEGER, INTENT(INOUT) :: I_A, I_B
        TYPE(vector_t) :: temp_data
        INTEGER :: temp_index

        ! Swap data
        temp_data = A
        A = B
        B = temp_data

        ! Swap index
        temp_index = I_A
        I_A = I_B
        I_B = temp_index
    END SUBROUTINE swap_both


    ! Modified Partition scheme to use swap_both
    SUBROUTINE partition_indexed(arr, indices, low, high, pivot_index)
        TYPE(vector_t), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER, DIMENSION(:), INTENT(INOUT) :: indices ! The index array
        INTEGER, INTENT(IN) :: low, high
        INTEGER, INTENT(OUT) :: pivot_index
        INTEGER :: i, j
        
        ! Select the last element (and its index) as the pivot
        TYPE(vector_t) :: pivot
        pivot = arr(high)

        i = low - 1 ! Index of the smaller element

        DO j = low, high - 1
            ! Use the custom comparison function (arr(j) < pivot)
            IF (is_less(arr(j), pivot)) THEN
                i = i + 1
                ! Swap both the vector and the index
                CALL swap_both(arr(i), arr(j), indices(i), indices(j))
            END IF
        END DO
        
        ! Place the pivot in the correct sorted position (swap with arr(i+1) and index(i+1))
        CALL swap_both(arr(i + 1), arr(high), indices(i + 1), indices(high))
        pivot_index = i + 1
    END SUBROUTINE partition_indexed


    ! Modified Recursive Quicksort routine
    RECURSIVE SUBROUTINE quicksort_indexed(arr, indices, low, high)
        TYPE(vector_t), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER, DIMENSION(:), INTENT(INOUT) :: indices ! The index array
        INTEGER, INTENT(IN) :: low, high
        INTEGER :: p

        IF (low < high) THEN
            ! Partition using the indexed version
            CALL partition_indexed(arr, indices, low, high, p)

            ! Recursively sort the sub-arrays
            CALL quicksort_indexed(arr, indices, low, p - 1)
            CALL quicksort_indexed(arr, indices, p + 1, high)
        END IF
    END SUBROUTINE quicksort_indexed

END MODULE Vector_Sort_Mod
