      SUBROUTINE WMCDIA(
     *                  DDOB72, 
     *                  DBUFF,  IBUFF,  CSOMO
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads MO coefficinet matrix from binary files and 
!                checks it.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WMCHKD
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!             ...declare variable(s) in argument list
!
      INTEGER,INTENT(IN) :: DDOB72
      REAL(KIND=LDREAL),DIMENSION(DDOB72),INTENT(INOUT) :: DBUFF
      INTEGER,DIMENSION(2,DDOB72),INTENT(INOUT) :: IBUFF
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO) :: CSOMO
!
!             ...define constant(s) in code
!
      REAL(KIND=LDREAL) :: SMALL = 1.0D-10
!
!             ...define variable(s) in code
!
      INTEGER :: ISO, IMO, ISYM
      INTEGER :: IS1, IS2, IM1, IM2
      INTEGER :: NW, KW, IW, MW
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMCDIA start'
!
!     ...read MO coefficient
!
      OPEN(UNIT=NFT72,FILE=FN72,STATUS='OLD',ACCESS='SEQUENTIAL',
     *     FORM='UNFORMATTED')
      CSOMO = 0.D0
!
      DO 
        READ(NFT72) KW,MW,DBUFF,IBUFF
        NW = ABS(KW)
        DO IW = 1, NW
          ISO = IBUFF(1,IW)
          IMO = IBUFF(2,IW)
          CSOMO(ISO,IMO) = DBUFF(IW)
        END DO
        IF( KW <= 0 ) EXIT
      END DO
      CLOSE(UNIT=NFT72)
!
!     ...set zero for block diagonal parts
!
      IS1 = 1
      IM1 = 1
      DO ISYM = 1, NSYM
        IF( NSOG(ISYM) == 0 ) CYCLE
        IS2 = IS1 + NSOG(ISYM) - 1
        IF( NMOG(ISYM) == 0 ) CYCLE
        IM2 = IM1 + NMOG(ISYM) - 1
        DO IMO = IM1, IM2
          DO ISO = IS1, IS2
            CSOMO(ISO,IMO) = 0.0D0
          END DO
        END DO
        IS1 = IS2 + 1
        IM1 = IM2 + 1
      END DO
!
!     ...check any value for off block diagonal parts
!
      DO IMO = 1, NMO
        DO ISO = 1, NSO
          IF( ABS(CSOMO(ISO,IMO)) > SMALL ) THEN
            WRITE(*,*) 'off block diagonal not zero'
            WRITE(*,*) ISO,IMO,CSOMO(ISO,IMO)
            STOP 'WMCDIA 001'
          END IF
        END DO
      END DO
!
!      WRITE(*,*) 'WMCDIA ended'
!
      RETURN
      END SUBROUTINE
