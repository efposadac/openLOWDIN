      SUBROUTINE WMFILE(
     *                  DD8BYT,
     *                  BIGARY, nproc, integralStackSize,
     *     otherNumberOfContractions)
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads SO basis integrals and MO coefficinet matrix
!                from text files and write them to binary files.
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
      integer :: nproc
      integer :: otherNumberOfContractions
      integer :: integralStackSize
!
!             ...define variable(s) in code
!
      INTEGER :: KEND, I
      INTEGER,DIMENSION(20) :: KKK, KK
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMFILE start'
!
!     ...save MO coefficient matrix in binary format
!
      KKK(1) = LDOB72 * (2 / LR2LI) ! IBUFF(2,LDOB71)
      KKK(2) = LDOB72 ! DBUFF(LDOB72)
      KKK(3) = NSO * NMO ! CSOMO(NSO,NMO)
      KEND = 3
      KK(1) = 1
      DO I = 1, KEND
        KK(I+1) = KK(I) + KKK(I)
      END DO
      CALL WMFCOE(LDOB72, 
     *            BIGARY(KK(1)), BIGARY(KK(2)), BIGARY(KK(3))  ) 
!
!     ...save SO basis integrals in binary format
!
      KKK(1) = LDOB73 * (4 / LR2LI) ! IBUFF(4,LDOB73)
      KKK(2) = LDOB73 ! DBUFF(LDOB73)
      KEND = 2
      KK(1) = 1
      DO I = 1, KEND
        KK(I+1) = KK(I) + KKK(I)
      END DO
      CALL WMFINT(LDOB73, 
     *            BIGARY(KK(1)), BIGARY(KK(2)), nproc, 
     *     integralStackSize, otherNumberOfContractions ) 
!
!      WRITE(*,*) 'WMFILE ended'
!
      RETURN
      END SUBROUTINE
