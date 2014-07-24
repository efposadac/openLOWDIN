      SUBROUTINE WTHDRV(
     *                  DD8BYT,
     *                  AA
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                is a main driver for WTH* module.
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   SYM4TR
!----------------------------------------------------------------------
!  WTHDRV-+-WTHSTB : set tables 
!         +-WTHRMO : read MO coefficients
!         +-WTHCL1 : 1st stage of transformation (PQ--->IJ)
!         +-WTHCL2 : 2nd stage of transformation (RS--->KL)
!----------------------------------------------------------------------
!  Arrays and Variables
!     IREPSO : irrep number of SO                           
!     IREPMO : irrep number of MO 
!     S2MRNG : MO range of valid symmetry for given SO
!              2*NSO*I4
!     M2SRNG : SO range of valid symmetry for given MO
!              2*NMO*I4
!     CSOMO  : MO coefficients (SOMO) / NSO*NMO*R8
!     NHTTR  : record position for backward chaining 
!              MTIMIJ*I4
!     IRSPAC : array to decompose IRSALL into IR and IS.
!            : This can be an integer*2 array.
!              2*NRSEXT
!     IJREV  : reversed canonical index ij to i and j
!            : This can be an integer*2 array.
!     WIQRS  : working area for transformed integral (i,q,r,s) 
!              NMO*NSO*R8
!     DIJRS  : half-transformed integral (i,j,r,s) 
!              NDIFRS*NMOMO*R8
!     WIJKS  : working area for transformed int. (i,j,k,s)
!     DBUFFM : buffer for reading MO matrix
!              NSO*NMO*R8
!     IBUFFM : buffer for reading MO matrix
!              NSO*NMO*2*I4
!     DBUFFR : buffer for reading ordered SO integrals 
!              NPQINT*R8
!     IBUFFR : buffer for reading indices of ordered SO integrals 
!              2*NPQINT*I4
!     DBUFFH : buffer for half-transformed integrals 
!              LDOB95*R8
!     IBUFFH : buffer for indices of half-transformed integrals
!              LDOB95*2*I4
!     DBUFFW : buffer for writing (i,j,r,s)
!     IBUFFW : buffer for reading index of (i,j,r,s)
!----------------------------------------------------------------------
!  Record structure of work DA (FT95)
!     In WTH*** version, only nonzero half-transformed integrals are
!     saved in DA.  Integrals are saved with their ij- and rs-indices.
!
!     NW,NPAD,DBUFFH(1:LDOB95),IBUFFH(1:2,1:LDOB95)
!
!     (4*2 + 8*LDOB95 + 4*2*LDOB95) .LE. (record length (BYTE) of DA).
!     In the last record within an ij-block, NW is set zero or minus.
!----------------------------------------------------------------------
!  Symmetry consideration
!     For Benzene D2h, WTHDRV is faster than VTQDRV by 4 or 5 times.
!     (1)
!       Only nonzero SO integrals (p,q,r,s) are used.
!     (2)
!       S2MRNG and M2SRNG are arrays to reduce the number of 
!       multiplication operation in WTHCL1 and WTHCL2 by symmetry 
!       consideration.  This reduction is based on blocked form of MO 
!       coefficient matrix, i.e., symmetry based molecular orbital.
!----------------------------------------------------------------------
!  I/O
!     input:        NFT94(seq)
!     intermediate: NFT95(da)
!     output:       NFT74(seq)
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DD8BYT
      REAL(KIND=LDREAL),DIMENSION(DD8BYT) :: AA
!
!             ...declare variable(s) in code
!
      INTEGER :: LDYM1, KX, NSYM1
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTHDRV start'
!      WRITE(*,*) 'THRHLF:',THRHLF
!      WRITE(*,*) 'THRMOI:',THRMOI
!
      NMOMO = (NMO * (NMO + 1)) / 2
!
!     ...calculates buffer size for half-transformed integrals (FT95)
!
      LDOB95 = (LDTRAK - LDINTG*2) / (LDINTG*2 + LDREAL)
      LREC95 = LDINTG*2 + (LDINTG*2 + LDREAL) * LDOB95
!
!     --- memory partition for WTHSTB (unit=8-Byte)
!
      KK8( 1) = NSO / LR2LI + N1PAD      ! IREPSO
      KK8( 2) = NMO / LR2LI + N1PAD      ! IREPMO
      KK8( 3) = NSO * 2 / LR2LI + N1PAD  ! S2MRNG
      KK8( 4) = NMO * 2 / LR2LI + N1PAD  ! M2SRNG
      KK8( 5) = NSYM * 2 /LR2LI + N1PAD  ! SOPAIR
      KK8( 6) = NSYM * 2 /LR2LI + N1PAD  ! MOPAIR
      KKLAST  = 6
      CALL WTMMPT
!
      CALL WTHSTB(
     *  AA(KK( 1)),AA(KK( 2)),AA(KK( 3)),AA(KK( 4)),AA(KK( 5)), 
     *  AA(KK( 6))  )
!       IREPSO      IREPMO      S2MRNG      M2SRNG      SOPAIR
!       MOPAIR
!
!     ...memory partition for WTHRMO
!
!"    (1)IREPSO
!"    (2)IREPMO
!"    (3)S2MRNG 
!"    (4)M2SRNG 
      KK8( 5) = NSO * NMO                    ! CSOMO
      KK8( 6) = LDOB72                       ! DBUFFM
      KK8( 7) = LDOB72 * 2 / LR2LI + N1PAD   ! IBUFFM
      KKLAST  = 7
      CALL WTMMPT
      CALL WTHRMO(
     *  AA(KK( 5)),AA(KK( 6)),AA(KK( 7)) )
!       CSOMO      DBUFFM     IBUFFM
!
!     ...optimize the value of NDIFIJ and NDIFRS.
!        determine other control variables such as NTIMIJ, 
!        NTIMRS,etc.
!
      CALL WTHOPT(DD8BYT)
!
!     ...memory partition for WTHCL1 (pq->ij) 
!
!"    KK8( 1) = NSO / LR2LI + N1PAD         ! IREPSO
!"    KK8( 2) = NMO / LR2LI + N1PAD         ! IREPMO
!"    KK8( 3) = NSO * 2 / LR2LI + N1PAD     ! S2MRNG
!"    KK8( 4) = NMO * 2 / LR2LI + N1PAD     ! M2SRNG
!"    KK8( 5) = NSO * NMO                   ! CSOMO
      KK8( 6) = NTIMIJ / LR2LI + N1PAD      ! NHTTR
      KK8( 7) = NRSEXT * 2 / LR2LI + N1PAD  ! IRSPAC
      KK8( 8) = LDOB95                      ! DBUFFH
      KK8( 9) = LDOB95 * 2 / LR2LI + N1PAD  ! IBUFFH
      KK8(10) = NMO * NSO                   ! WIQRS
      KK8(11) = NDIFRS * NMOMO              ! DIJRS
      KK8(12) = LDOB94                      ! DBUFFR
      KK8(13) = LDOB94 * 2 / LR2LI + N1PAD  ! IBUFFR
      KKLAST = 13
      CALL WTMMPT   
      IF( KK(KKLAST+1) > DD8BYT ) THEN
        STOP 'WTHDRV 001'
      END IF
      CALL WTHCL1(
     *  AA(KK( 1)),AA(KK( 2)),AA(KK( 3)),AA(KK( 4)),AA(KK( 5)),
     *  AA(KK( 6)),AA(KK( 7)),AA(KK( 8)),AA(KK( 9)),AA(KK(10)),
     *  AA(KK(11)),AA(KK(12)),AA(KK(13)),AA(KK(14))  )
!       IREPSO     IREPMO     S2MRNG     M2SRNG     CSOMO
!       NHTTR      IRSPAC     DBUFFH     IBUFFH     WIQRS 
!       DIJRS      DBUFFR     IBUFFR
!
!     ...memory partition for WTHCL2 (rs->kl)
!
      KK8(10) = NSO*NMO                     ! WIJKS
      KK8(11) = NRSEXT * NDIFIJ             ! DIJRS
      KK8(12) = LDOB74                      ! DBUFFW
      KK8(13) = LDOB74 * 4 / LR2LI + N1PAD  ! IBUFFW
      KK8(14) = NMOMO  * 2 / LR2LI + N1PAD  ! IJREV
      KKLAST = 14
      CALL WTMMPT   
      IF( KK(KKLAST+1) > DD8BYT ) THEN
        STOP 'WTHDRV 002'
      END IF
!
      CALL WTHCL2(
     *             AA(KK( 2)),AA(KK( 3)),AA(KK( 4)),AA(KK( 5)),
     *  AA(KK( 6)),AA(KK( 7)),AA(KK( 8)),AA(KK( 9)),AA(KK(10)),
     *  AA(KK(11)),AA(KK(12)),AA(KK(13)),AA(KK(14))  )
!                  IREPMO     S2MRNG     M2SRNG     CSOMO
!       NHTTR      IRSPAC     DBUFFH     IBUFFH     WIJKS 
!       DIJRS      DBUFFW     IBUFFW     IJREV
!
!      IF( NPFLAG > 0 ) THEN
!        CALL WTHPRI(
!     *             AA(KK(12)),AA(KK(13))  )
!                  DBUFFW     IBUFFW     
!      END IF
!
!      WRITE(*,*) 'WTHDRV ended'
      RETURN
      END SUBROUTINE
