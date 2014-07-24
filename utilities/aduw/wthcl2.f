      SUBROUTINE WTHCL2(
     *                  IREPMO, S2MRNG, M2SRNG, CSOMO,
     *                  NHTTR,  IRSPAC, DBUFFH, IBUFFH, WIJKS,
     *                  DIJRS,  DBUFFW, IBUFFW, IJREV
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                transforms (ij|rs) to (ij|kl).
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
!  N.B.
!     The following data are passed from subroutine WTHCL1 through 
!     memory.
!       CSOMO, S2MRNG, M2SRNG, NHTTR, IRSPAC
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
      INCLUDE 'nmltbl.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSO),INTENT(IN) :: IREPMO
      INTEGER,DIMENSION(2,NSO),INTENT(IN) :: S2MRNG
      INTEGER,DIMENSION(2,NMO),INTENT(IN) :: M2SRNG
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO),INTENT(IN) :: CSOMO
!
      INTEGER,DIMENSION(NTIMIJ),INTENT(INOUT) :: NHTTR
      INTEGER,DIMENSION(2,NRSEXT),INTENT(INOUT) :: IRSPAC
      REAL(KIND=LDREAL),DIMENSION(LDOB95) :: DBUFFH
      INTEGER,DIMENSION(2,LDOB95) :: IBUFFH
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO) :: WIJKS
      REAL(KIND=LDREAL),DIMENSION(NRSEXT,NDIFIJ) :: DIJRS
      REAL(KIND=LDREAL),DIMENSION(LDOB74) :: DBUFFW
      INTEGER,DIMENSION(4,LDOB74) :: IBUFFW
      INTEGER,DIMENSION(2,NMOMO) :: IJREV
!
!             ...declare variable(s) in code
!
      INTEGER :: I, J, IJ
      INTEGER :: K, L, KL 
      INTEGER :: IR, IS, IRS
      INTEGER :: J01
      INTEGER :: KBUF, KNG
      INTEGER :: ITIMRS, IRDTTR, KHTTR, ITIMIJ
      INTEGER :: KDIFIJ, MIJEND, IDIFIJ
      INTEGER :: MIJTOP, KDIFRS, MRSTOP, MRSEND
      INTEGER :: SYMIJ, SYMKL
      REAL(KIND=LDREAL) :: WIJKL
!
!----------------------------------------------------------------------
 9000 FORMAT(1X,'ij-integral block is remained in SUB.WTHCL2,  ',
     *          'ITIMIJ,NHTTR=',2I8)
!----------------------------------------------------------------------
!
!      WRITE(*,*) ' '
!      WRITE(*,*) 'SECOND STAGE OF INTEGRALS TRANSFORMATION'
!      CALL PRCPUT('WTHCL2 start',0)
!
      CALL VTHREV(NMO,NMOMO,IJREV)
      OPEN(UNIT=NFT74,FILE=FN74,STATUS='UNKNOWN',
     *     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      KBUF = 0
!
!     ...............................
      DO_ITIMIJ: DO ITIMIJ = 1, NTIMIJ
!     ...............................
!
!     ---- 3rd step (i,j,r,s) -> (i,j,k,s) ----
!
        DIJRS = 0.D0
        ITIMRS = 0
!
!       ..........
        RS_ALL: DO
!       ..........
          ITIMRS = ITIMRS + 1
          IF( ITIMRS > NTIMRS ) GOTO 901
          IRDTTR = NHTTR(ITIMIJ)
!
          CALL WTHRDA(
     *                IRDTTR, KHTTR, 
     *                KDIFIJ, KDIFRS, 
     *                MRSTOP, MRSEND, MIJTOP, MIJEND,
     *                DIJRS,  DBUFFH, IBUFFH)
!
          NHTTR(ITIMIJ) = KHTTR
          IF( NHTTR(ITIMIJ) == 0 ) EXIT
!       .............
        END DO RS_ALL
!       .............
!
!       ................................
        DO_IDIFIJ: DO IDIFIJ = 1, KDIFIJ
!       ................................
!
          WIJKS = 0.D0
!
          DO IRS = 1, NRSEXT
            IR = IRSPAC(1,IRS)
            IS = IRSPAC(2,IRS)
!           The next if-statement is effective for speed-up
            IF( ABS(DIJRS(IRS,IDIFIJ)) < THRHLF ) CYCLE
            IF( IR == IS ) THEN
              DO K = S2MRNG(1,IR), S2MRNG(2,IR)
                WIJKS(IS,K) = WIJKS(IS,K) 
     *                      + DIJRS(IRS,IDIFIJ) * CSOMO(IR,K)
              END DO
            ELSE
              DO K = S2MRNG(1,IR), S2MRNG(2,IR)
                WIJKS(IS,K) = WIJKS(IS,K)
     *                      + DIJRS(IRS,IDIFIJ) * CSOMO(IR,K)
              END DO
              DO K = S2MRNG(1,IS), S2MRNG(2,IS)
                WIJKS(IR,K) = WIJKS(IR,K)
     *                      + DIJRS(IRS,IDIFIJ) * CSOMO(IS,K)
              END DO
            END IF
          END DO 
!
!     ---- 4th step (i,j,k,s) -> (i,j,k,l) ----
!
          IJ = MIJTOP + IDIFIJ - 1
          I = IJREV(1,IJ)
          J = IJREV(2,IJ)
          SYMIJ = NMLTBL(IREPMO(I),IREPMO(J))
!
          KL = 0
          DO_K: DO K = 1, NMO
            DO_L: DO L = 1, K
              KL = KL + 1
!             SAVE ONLY WHEN (IJ>=KL)
              IF( IJ < KL ) CYCLE
              SYMKL = NMLTBL(IREPMO(K),IREPMO(L))
              IF( SYMIJ /= SYMKL ) CYCLE
cccc
              IF( M2SRNG(1,L) == 0 ) CYCLE
cccc
              WIJKL = 0.D0
              DO IS = M2SRNG(1,L), M2SRNG(2,L)
                WIJKL = WIJKL + CSOMO(IS,L) * WIJKS(IS,K)
              END DO
!
!     ---- save MO integrals in FT74 (sequential file) ----
!
              IF( ABS(WIJKL) < THRMOI ) CYCLE
              KBUF = KBUF + 1
              IBUFFW(1,KBUF) = I
              IBUFFW(2,KBUF) = J
              IBUFFW(3,KBUF) = K
              IBUFFW(4,KBUF) = L
              DBUFFW(KBUF) = WIJKL
              IF( KBUF == LDOB74 ) THEN
                WRITE(NFT74) KBUF,KBUF,DBUFFW,IBUFFW
                KBUF = 0
              END IF
!
            END DO DO_L
          END DO DO_K
!       ................
        END DO DO_IDIFIJ
!       ................
!
!     ................
      END DO DO_ITIMIJ
!     ................
!
!     ---- save transformed integrals in FT74 (sequential file) ----
!     ---- last record                                          ----
!
      WRITE(NFT74) -KBUF,LDOB74,DBUFFW,IBUFFW
      CLOSE(NFT74)
      IF( NRFLAG == 0 ) THEN
        CLOSE(NFT95,STATUS='DELETE')
      ELSE
        CLOSE(NFT95)
      END IF
!
      KNG = 0
      DO ITIMIJ = 1, NTIMIJ
        IF( NHTTR(ITIMIJ) /= 0 ) THEN
          WRITE(*,9000) ITIMIJ,NHTTR(ITIMIJ)
          KNG = 1
        END IF
      END DO
      IF( KNG == 1 ) GOTO 902
!
!      CALL PRCPUT('Second stage ended',0)
!      WRITE(*,*) 'Integral transformation completed ...'
!      WRITE(*,*) ' '
      RETURN
!
  901 CONTINUE
      WRITE(*,*) ITIMRS,NTIMRS
      STOP 'WTHCL2 001'
  902 CONTINUE
      STOP 'WTHCL2 002'
      END SUBROUTINE
