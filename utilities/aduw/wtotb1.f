      SUBROUTINE WTOTB1(
     *                  ID1BOX
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets control tables for ordering SO basis integrals.
!                The i-th index belongs to ID1BOX(i)-th box.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSOSO),INTENT(OUT) :: ID1BOX
!
!             ...declare variable(s) in code
!
      INTEGER :: I1BOX, ISEQ, ISEQ1, ISEQ2
!
!---------------------------------------------------------------------
!      WRITE(*,*) 'WTOTB1 start'
!
      DO I1BOX = 1, NBOX
        ISEQ1 = NBXDIF*(I1BOX - 1) + 1
        ISEQ2 = ISEQ1 + NBXDIF - 1
        ISEQ2 = MIN(ISEQ2,NSOSO)
        IF(ISEQ1 > ISEQ2) THEN
 !         WRITE(*,*) 'ISEQ1,ISEQ2:',ISEQ1,ISEQ2
          STOP 'WTOTB1 001'
        END IF
        DO ISEQ = ISEQ1, ISEQ2
          ID1BOX(ISEQ) = I1BOX
        END DO
      END DO
!
!     WRITE(*,*) 'WTOTB1 ended'
      RETURN
      END SUBROUTINE
