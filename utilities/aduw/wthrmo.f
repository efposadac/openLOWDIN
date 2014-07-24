      SUBROUTINE WTHRMO(
     *                  CSOMO,DBUFFM,IBUFFM
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads MO coefficient matrix.
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO) :: CSOMO
      REAL(KIND=LDREAL),DIMENSION(LDOB72) :: DBUFFM
      INTEGER,DIMENSION(2,LDOB72) :: IBUFFM
!
!             ...declare variable(s) in code
!
      INTEGER :: KW, NW, MW
      INTEGER :: IW, ISO, IMO
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTHRMO start'
!
      CSOMO = 0.0D0
      OPEN(UNIT=NFT72,FILE=FN72,STATUS='OLD',
     *     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!
      DO
        READ(NFT72) KW,MW,DBUFFM,IBUFFM
        NW = ABS(KW)
        DO IW = 1, NW
          ISO = IBUFFM(1,IW)
          IMO = IBUFFM(2,IW)
          CSOMO(ISO,IMO) = DBUFFM(IW)
        END DO
        IF( KW <= 0 ) EXIT
      END DO
!
      CLOSE(UNIT=NFT72)
!      IF( NPFLAG > 0 ) CALL WMUPMO(NSO,NMO,CSOMO)
!
!      WRITE(*,*) 'WTHRMO ended'
      RETURN
      END SUBROUTINE
