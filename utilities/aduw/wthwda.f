      SUBROUTINE WTHWDA(
     *                  NREC,   KHTTR,  KDIFIJ, KDIFRS, 
     *                  MRSTOP, MRSEND, MIJTOP, MIJEND, 
     *                  DIJRS,  DBUFFH, IBUFFH
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                writes half-transformed integrals in a work direct
!                access file.
!  RELEASE :     v.00  gen = sya 1999-01-19 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHCL1
!----------------------------------------------------------------------
!  N.B.
!     (1) Copy only nonzero half-transformed into DBUFFH buffer.
!     (2) Write IBUFFH and DBUFFH in a work direct access file.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(INOUT) :: NREC
      INTEGER,INTENT(INOUT) :: KHTTR(1)
      INTEGER,INTENT(IN) :: KDIFIJ, KDIFRS 
      INTEGER,INTENT(IN) :: MRSTOP, MRSEND, MIJTOP, MIJEND
      REAL(KIND=LDREAL),INTENT(IN) :: DIJRS(NDIFRS,NMOMO)
!
      REAL(KIND=LDREAL),DIMENSION(LDOB95) :: DBUFFH
      INTEGER,DIMENSION(2,LDOB95) :: IBUFFH
!
!             ...declare constant(s) in code
!
      INTEGER :: IPAD = 0
!
!             ...declare variable(s) in code
!
      INTEGER :: IRS, IJALL, IRSALL, NW, IOS
!
!----------------------------------------------------------------------
      NREC = NREC + 1
      WRITE(NFT95,REC=NREC,ERR=901,IOSTAT=IOS) 
     *  KHTTR(1),
     *  KDIFIJ,KDIFRS,MRSTOP,MRSEND,MIJTOP,MIJEND,IPAD
      KHTTR(1) = NREC
!
      NW = 0
      DO IRS = 1, KDIFRS
        IRSALL = IRS + MRSTOP - 1
        DO IJALL = MIJTOP, MIJEND
          IF( ABS(DIJRS(IRS,IJALL)) < THRHLF ) CYCLE
          NW = NW + 1
          IBUFFH(1,NW) = IJALL
          IBUFFH(2,NW) = IRSALL
          DBUFFH(NW) = DIJRS(IRS,IJALL)
          IF( NW  == LDOB95 ) THEN
            NREC = NREC + 1
            WRITE(NFT95,REC=NREC,ERR=902,IOSTAT=IOS)
     *        NW,IPAD,DBUFFH,IBUFFH
            NW = 0
          END IF
        END DO
      END DO
!
!     ...write the last record in this block.
!
      NREC = NREC + 1
      WRITE(NFT95,REC=NREC,ERR=903,IOSTAT=IOS)
     *  -NW,IPAD,DBUFFH,IBUFFH
      RETURN
!
  901 CONTINUE
      WRITE(*,*) 'IOSTAT=',IOS
      WRITE(*,*) 'KHTTR(1),KDIFIJ,KDIFRS'
      WRITE(*,*)  KHTTR(1),KDIFIJ,KDIFRS
      WRITE(*,*) 'MRSTOP,MRSEND,MIJTOP,MIJEND'
      WRITE(*,*)  MRSTOP,MRSEND,MIJTOP,MIJEND
      STOP 'WTHWDA 901'
!
  902 CONTINUE
      WRITE(*,*) 'IOSTAT=',IOS
      WRITE(*,*) 'NREC=',NREC
      STOP 'WTHWDA 902'
!
  903 CONTINUE
      WRITE(*,*) 'IOSTAT=',IOS
      WRITE(*,*) 'NREC=',NREC
      STOP 'WTHWDA 903'
!
      END SUBROUTINE
