      SUBROUTINE WTHCL1(
     *                  IREPSO, IREPMO, S2MRNG, M2SRNG, CSOMO,
     *                  NHTTR,  IRSPAC, DBUFFH, IBUFFH, WIQRS,
     *                  DIJRS,  DBUFFR, IBUFFR
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                transforms SO integral (pq|rs) to (ij|rs).
!  RELEASE :     v.00  gen = sya 1999-01-18 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
!  N.B.
!     IRSPAC contains r and s sequential indexes for compressed rs-pair
!     index.  This compressed rs-pair index differs from normal 
!     canonical index ( (r*(r-1))/2 + s ) of r and s.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
      INCLUDE 'nmltbl.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSO),INTENT(IN) :: IREPSO
      INTEGER,DIMENSION(NMO),INTENT(IN) :: IREPMO
      INTEGER,DIMENSION(2,NSO),INTENT(IN) :: S2MRNG
      INTEGER,DIMENSION(2,NMO),INTENT(IN) :: M2SRNG
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO),INTENT(IN) :: CSOMO
      INTEGER,DIMENSION(NTIMIJ),INTENT(OUT) :: NHTTR
      INTEGER,DIMENSION(2,NRSEXT),INTENT(OUT) :: IRSPAC
      REAL(KIND=LDREAL),DIMENSION(LDOB95) :: DBUFFH
      INTEGER,DIMENSION(2,LDOB95) :: IBUFFH
      REAL(KIND=LDREAL),DIMENSION(NMO,NSO) :: WIQRS
      REAL(KIND=LDREAL),DIMENSION(NDIFRS,NMOMO) :: DIJRS
      REAL(KIND=LDREAL),DIMENSION(LDOB94) :: DBUFFR
      INTEGER,DIMENSION(2,LDOB94) :: IBUFFR
!
!             ...declare variable(s) in code
!
      INTEGER :: IRSALL
      INTEGER :: NREC, IERROR, IOSTAT, KREC, NXTREC 
      INTEGER :: MINT, IR, IS, IP, IQ
      INTEGER :: MIJEND, MIJTOP, KDIFIJ, KDIFRS, MRSEND, MRSTOP
      INTEGER :: ITIMIJ, ITIMRS, IJ, IINT, IDIFRS, I, J
      INTEGER :: SYMRS, SYMIJ
      INTEGER :: IPAD
      REAL(KIND=LDREAL) :: WIJRS
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) " "
!      WRITE(*,*) 'FIRST STAGE OF INTEGRALS TRANSFORMATION'
!      CALL PRCPUT('WTHCL1 start',0)
!

      NHTTR = 0
      IRSPAC = 0
      NREC = 0
!
      REWIND NFT94
      CALL SETDRW(NFT95,LREC95,IERROR,IOSTAT)
!
!     ................................
      DO_ITIMRS: DO ITIMRS = 1, NTIMRS
!     ................................
!
        MRSTOP = (ITIMRS - 1) * NDIFRS + 1
        MRSEND = ITIMRS * NDIFRS
        KDIFRS = NDIFRS
        IF( ITIMRS == NTIMRS ) THEN
          MRSEND = NRSEXT
          KDIFRS = MRSEND - MRSTOP + 1
        ELSE
          IF( MRSEND > NRSEXT ) THEN
            WRITE(*,*) 'ITIMRS,NTIMRS,MRSEND,NRSEXT'
            WRITE(*,*) ITIMRS,NTIMRS,MRSEND,NRSEXT
            STOP 'WTHCL1 001'
          END IF
        END IF
!
        DIJRS(1:KDIFRS,:) = 0.D0
!
!       ................................
        DO_IDIFRS: DO IDIFRS = 1, KDIFRS
!       ................................
!
          WIQRS = 0.D0
!
!     ---- 1st step (p,q,r,s) -> (i,q,r,s) ----
!
          READ(NFT94,END=799) MINT,IR,IS,IPAD,
     *      DBUFFR(1:MINT),IBUFFR(:,1:MINT)
          IRSALL = IDIFRS + MRSTOP - 1
          IRSPAC(1,IRSALL) = IR
          IRSPAC(2,IRSALL) = IS
          SYMRS = NMLTBL(IREPSO(IR),IREPSO(IS))
!
          DO_IINT: DO IINT = 1, MINT
            IP = IBUFFR(1,IINT)
            IQ = IBUFFR(2,IINT)
            IF( IP == IQ ) THEN
              DO I = S2MRNG(1,IP), S2MRNG(2,IP)
                WIQRS(I,IQ) = WIQRS(I,IQ) + DBUFFR(IINT) * CSOMO(IP,I)
              END DO
            ELSE
              DO I = S2MRNG(1,IP), S2MRNG(2,IP)
                WIQRS(I,IQ) = WIQRS(I,IQ) + DBUFFR(IINT) * CSOMO(IP,I)
              END DO
              DO I = S2MRNG(1,IQ), S2MRNG(2,IQ)
                WIQRS(I,IP) = WIQRS(I,IP) + DBUFFR(IINT) * CSOMO(IQ,I)
              END DO
            END IF
          END DO DO_IINT
!
!     ---- 2nd step  (i,q,r,s) -> (i,j,r,s) ----
!
          IJ = 0
          DO I = 1, NMO
            DO J = 1, I
              IJ = IJ + 1
              SYMIJ = NMLTBL(IREPMO(I),IREPMO(J))
              IF( SYMRS /= SYMIJ ) CYCLE
              IF( M2SRNG(1,J) == 0 ) CYCLE
              WIJRS = 0.D0
              DO IQ = M2SRNG(1,J), M2SRNG(2,J)
                WIJRS = WIJRS + CSOMO(IQ,J) * WIQRS(I,IQ)
              END DO
              DIJRS(IDIFRS,IJ) = WIJRS
            END DO
          END DO
!
!       ................
        END DO DO_IDIFRS
!       ................
!
!     ---- SAVE (I,J,R,S) TO DAFILE (NFT95) ----
!
        DO_ITIMIJ: DO ITIMIJ = 1, NTIMIJ
          MIJTOP = NDIFIJ * (ITIMIJ - 1) + 1
          MIJEND = NDIFIJ * ITIMIJ
          KDIFIJ = NDIFIJ
          IF( NMOMO < MIJEND ) THEN
            MIJEND = NMOMO
            KDIFIJ = MIJEND - MIJTOP + 1
          END IF
!
          CALL WTHWDA(
     *      NREC,   NHTTR(ITIMIJ), 
     *      KDIFIJ, KDIFRS, 
     *      MRSTOP, MRSEND, MIJTOP, MIJEND, 
     *      DIJRS,  DBUFFH, IBUFFH)
!
        END DO DO_ITIMIJ
!
!     ................
      END DO DO_ITIMRS
!     ................
!
!
!     ...confirm EOF for safe. (can be skipped)
!
      INQUIRE(UNIT=NFT95,NEXTREC=NXTREC)
      READ(NFT94,END=799) MINT,IR,IS,IPAD,
     *  DBUFFR(1:MINT),IBUFFR(:,1:MINT)
      STOP 'WTHCL1 002'
!
  799 CONTINUE
!
      KREC = NREC + 1
      IF( NXTREC /= KREC ) THEN
        WRITE(*,*) 'NXTREC,NREC',NXTREC,NREC
        STOP 'WTHCL1 003'
      END IF
      IF( MIJEND /= NMOMO ) THEN
        WRITE(*,*) 'MIJEND,NMOMO',MIJEND,NMOMO
        STOP 'WTHCL1 004'
      END IF
!
      IF( NRFLAG == 0 ) THEN
        CLOSE(NFT94,STATUS='DELETE')
      ELSE
        CLOSE(NFT94)
      END IF
!
!      CALL PRCPUT('',0)
!      WRITE (*,*) 'First stage ended ...'
      
      RETURN
!
      WRITE(*,*) 'ITIMIJ,NTIMIJ',ITIMIJ,NTIMIJ
      STOP 'WTHCL1 005'
      END SUBROUTINE
