      SUBROUTINE WTHPRI(
     *                  DBUFFW, IBUFFW
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                prints transformed integrals and their indices.
!  RELEASE :     v.01  gen = sya 2003-04-18 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      REAL(KIND=LDREAL),DIMENSION(LDOB74) :: DBUFFW
      INTEGER,DIMENSION(4,LDOB74) :: IBUFFW
!
!             ...declare variable(s) in code
!
      INTEGER :: I, IOS
      INTEGER :: KBUF, NBUF, MBUF
      CHARACTER(LEN=60) :: FN61 
      CHARACTER(LEN=60) :: FN62 
      CHARACTER(LEN=60) :: FN63 
      CHARACTER(LEN=60) :: FN71 
      CHARACTER(LEN=60) :: FN72 
      CHARACTER(LEN=60) :: FN73 
      CHARACTER(LEN=60) :: FN74

      COMMON/FILES/FN61,FN62,FN63,FN71,FN72,FN73,FN74

!
!----------------------------------------------------------------------
 1000 FORMAT(1X,'TRANSFORMATED INTEGRALS')
 1010 FORMAT(1X,'KBUF,NBUF=',2I5)
 1020 FORMAT((1X,3(4I4,1X,E18.10,1X)))
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTHPRI start'
!
      OPEN(UNIT=NFT74,FILE=FN74,STATUS='OLD',
     *     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!
      WRITE(*,1000) 
!
      DO 
        READ(UNIT=NFT74,ERR=901,IOSTAT=IOS) KBUF,NBUF,DBUFFW,IBUFFW
!        WRITE(*,1010) KBUF,NBUF
        MBUF = IABS(KBUF)
        IF( KBUF /= 0 ) THEN
          WRITE(*,1020) (IBUFFW(:,I),DBUFFW(I),I=1,MBUF)
        END IF
        IF( KBUF <= 0 ) EXIT
      END DO
!
      CLOSE(NFT74)
!      WRITE(*,*) 'WTHPRI ended'
      RETURN
!
  901 CONTINUE
!      WRITE(*,*) IOS
      STOP 'WTHPRI 001'
      END SUBROUTINE
