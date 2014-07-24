      SUBROUTINE VTHREV(NM,NMM,IJREV)
!----------------------------------------------------------------------
!  $$$$TRANS5M
!  Set reversed canonical index
!  Get i and j from their canonical index ij
!  ij = (i * (i-1)) / 2 + j
!  VTMAIN-VTHDRV-VTHCL2-VTHREV
!  CREATION = S.YAMAMOTO 1999-01-19
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: NM, NMM
      INTEGER,INTENT(OUT) :: IJREV(2,NMM)
!
!     working variables within this routine
      INTEGER             :: IJ, I, J
!----------------------------------------------------------------------
      IJREV = 0
      IJ = 0
      DO I = 1, NM
        DO J = 1, I
          IJ = IJ + 1
          IJREV(1,IJ) = I
          IJREV(2,IJ) = J
        END DO
      END DO
      RETURN
      END
