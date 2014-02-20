      SUBROUTINE WTOTB2(
     *                  ID2BOX
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets tables for 2nd stage of ordering SO integrals.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
!  N.B.
!     ISEQ      : sequential number of rs-index
!     ID2BOX(i) : i-th rs-index should be distributed to ID2BOX(i)-th
!                 box in bin-2.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSOSO),INTENT(OUT) :: ID2BOX
!
!             ...declare variable(s) in code
!
      INTEGER :: ISEQ, I1BOX, MAX2BX, I2BOX
!
!---------------------------------------------------------------------
!      WRITE(*,*) 'WTOTB2 start'
!
      ISEQ = 0
ccc   MAX1BX = NSOSO/NBXDIF
ccc   IF( (MAX1BX*NBXDIF) < NSOSO ) MAX1BX = MAX1BX + 1
ccc   DO I1BOX = 1, MAX1BX
      DO I1BOX = 1, NBOX
        MAX2BX = MIN( (NSOSO-ISEQ),NBXDIF )
        DO I2BOX = 1, MAX2BX
          ISEQ = ISEQ + 1
          ID2BOX(ISEQ) = I2BOX
        END DO
      END DO
!
!      WRITE(*,*) 'WTOTB2 ended'
      RETURN
      END SUBROUTINE
