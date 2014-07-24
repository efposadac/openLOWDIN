      SUBROUTINE WMCHKD(
     *                  DD8BYT,
     *                  BIGARY
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                checks whether MO coefficient matrix is 
!                block-diagonal.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   PREPAR
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DD8BYT
      REAL(KIND=LDREAL),DIMENSION(DD8BYT) :: BIGARY
!
!             ...define variable(s) in code
!
      INTEGER :: KEND, I
      INTEGER,DIMENSION(20) :: KKK, KK
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMCHKD start'
!
!     ...save MO coefficient matrix in binary format
!
      KKK(1) = LDOB72 ! DBUFF(LDOB71)
      KKK(2) = LDOB72 * 2 / LR2LI + N1PAD ! IBUFF(2,LDOB71)
      KKK(3) = NSO * NMO ! CSOMO(NSO,NMO)
      KEND = 3
      KK(1) = 1
      DO I = 1, KEND
        KK(I+1) = KK(I) + KKK(I)
      END DO
      CALL WMCDIA(LDOB72,
     *            BIGARY(KK(1)), BIGARY(KK(2)), BIGARY(KK(3))  ) 
!
!      WRITE(*,*) 'WMCHKD ended'
!
      RETURN
      END SUBROUTINE
