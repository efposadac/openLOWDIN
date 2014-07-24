      SUBROUTINE WTMMPT
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets memory partitioning.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTOMP1
!                WTOMP2
!                WTHDRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'sym4tr.h'
!
!             ...define variable(s) in code
!
      INTEGER :: I
!
!----------------------------------------------------------------------
!
      KK(1) = 1
      DO I = 1, KKLAST
        KK(I+1) = KK(I) + KK8(I)
      END DO
!
      RETURN
      END SUBROUTINE
