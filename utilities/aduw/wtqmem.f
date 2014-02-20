      SUBROUTINE WTQMEM(
     *                  IGO1
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                checks memory size.
!                If OK, IGO1 -> 1.  If NG, IGO1 -> -1,-2,-3.
!  RELEASE :     v.00  gen = sya 1987-06-24 at ims
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHOPT
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(INOUT) :: IGO1
!
!             ...declare variable(s) in code
!
      INTEGER :: K
!
!---------------------------------------------------------------------
!
      IGO1 = 0
      K = NTIMIJ + NDIFRS * NMOMO
      IF( K > M8CORE ) IGO1 = IGO1 - 1
      K = NTIMIJ + NDIFIJ * NRSEXT
      IF( K > M8CORE ) IGO1 = IGO1 - 2
!
      IF( IGO1 == 0 ) THEN
        IGO1 = 1
      END IF
!
      RETURN
      END SUBROUTINE
