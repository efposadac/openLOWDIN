      SUBROUTINE WTHRDA(
     *                  NREC,   KHTTR, 
     *                  KDIFIJ, KDIFRS, 
     *                  MRSTOP, MRSEND, MIJTOP, MIJEND,
     *                  DIJRS,  DBUFFH, IBUFFH
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads half-transformed integrals from work direct
!                access file.
!  RELEASE :     v.00  gen = sya 1999-01-19 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHCL2
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(INOUT) :: NREC
      INTEGER,INTENT(OUT) :: KHTTR
      INTEGER,INTENT(OUT) :: KDIFIJ, KDIFRS 
      INTEGER,INTENT(OUT) :: MRSTOP, MRSEND, MIJTOP, MIJEND
      REAL(KIND=LDREAL),DIMENSION(NRSEXT,NDIFIJ),INTENT(OUT) :: DIJRS
!
      REAL(KIND=LDREAL),DIMENSION(LDOB95) :: DBUFFH
      INTEGER,DIMENSION(2,LDOB95) :: IBUFFH
!
!             ...declare variable(s) in code
!
      INTEGER :: MT1, KW, NW, IW, IRSALL, IJALL, IJ, IOS
      INTEGER :: IPAD
!----------------------------------------------------------------------
 9010 FORMAT(/1X,'**** ERROR STOP IN SUB.WTHRDA ****',
     *       /1X,'     END OF FILE DURING READING.',
     *       /1X,'     NREC, IOS = ',2I12)
!----------------------------------------------------------------------
      READ(NFT95,REC=NREC,ERR=901,IOSTAT=IOS) 
     *  KHTTR,KDIFIJ,KDIFRS,MRSTOP,MRSEND,MIJTOP,MIJEND,IPAD
!
      MT1 = MIJTOP - 1
      DO 
        NREC = NREC + 1
        READ(NFT95,REC=NREC,ERR=902,IOSTAT=IOS)
     *    KW,IPAD,DBUFFH,IBUFFH
        NW = IABS(KW)
        DO IW = 1, NW
          IJALL = IBUFFH(1,IW)
          IRSALL = IBUFFH(2,IW)
          IJ = IJALL - MT1
          DIJRS(IRSALL,IJ) = DBUFFH(IW)
        END DO
        IF( KW <= 0 ) EXIT
      END DO
!
      RETURN
!
  901 CONTINUE
      WRITE(*,9010) NREC,IOS 
      STOP 'WTHRDA 901'
  902 CONTINUE
      WRITE(*,9010) NREC,IOS
      STOP 'WTHRDA 902'
      END
