      SUBROUTINE WMFCOE(
     *                  DDOB72, 
     *                  IBUFF,  DBUFF,  CSOMO
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads MO coefficinet matrix from text files and 
!                saves it to a binary file.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WMFILE
!----------------------------------------------------------------------
!  structure of MO matrix file of text format (mocoef.txt).
!     MO index, SO index, matrix element (coefficient) per line.
!
!  record structure of binary MO matrix file
!     KW, LDOB72, DBUFF(1:LDOB72), IBUFF(1:2,1:LDOB72)
! 
!  meaning of variable
!     NW : number of elements of MO coefficients in a record
!     KW : = NW, 
!          = 0 or -NW for the last record
!     DBUFF(i) : coefficient value
!     IBUFF(1,i) : SO number
!     IBUFF(2,i) : MO number
!
!     LDOB72 == 1187
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
!             ...define variable(s) in code
!
      REAL(KIND=LDREAL) :: VAL
      INTEGER :: ISO, IMO
      INTEGER :: NW
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMFCOE start'
!
!     ...read MO coefficient
!
      OPEN(UNIT=NFT62,FILE=FN62,STATUS='OLD',ACCESS='SEQUENTIAL',
     *     FORM='FORMATTED')
      OPEN(UNIT=NFT72,FILE=FN72,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     *     FORM='UNFORMATTED')
      CSOMO = 0.D0
!
      NW = 0
      DO 
        READ(NFT62,FMT=*,END=109) IMO,ISO,VAL
        IF( IMO > NMO ) STOP 'WMFCOE 001'
        IF( ISO > NSO ) STOP 'WMFCOE 002' 
        CSOMO(ISO,IMO) = VAL
        NW = NW + 1
        IBUFF(1,NW) = ISO
        IBUFF(2,NW) = IMO
        DBUFF(NW) = VAL
        IF( NW == LDOB72 ) THEN
          WRITE(NFT72) NW,LDOB72,DBUFF,IBUFF
          NW = 0
        END IF
      END DO
  109 CONTINUE
      WRITE(NFT72) -NW,LDOB72,DBUFF,IBUFF
!
      CLOSE(UNIT=NFT62)
      CLOSE(UNIT=NFT72)
!
!      IF( NPFLAG > 0 ) CALL WMUPMO(NSO,NMO,CSOMO)
!
!      WRITE(*,*) 'WMFILE ended'
!
      RETURN
      END SUBROUTINE
