      SUBROUTINE WTHOPT(
     *                  DD8BYT
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                optimizes NDIFIJ and NDIFRS and determines other 
!                control variables, i.e., NTIMIJ, etc.
!                This version saves only nonzero half-transformed 
!                integrals.
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
!
!  Optimize the value of 'NDIFIJ'.
!
!     *----------------------------------------------> rs-index
!     |
!     |  *----------------------------*     /
!     |  |4          |4          |2   |     | NDIFIJ
!     |  |           |           |    |     |
!     |  |          A|          D|   G|     |
!     |  |-----------+-----------+----|     /
!     |  |4          |4          |2   |
!     |  |           |           |    |
!     |  |          B|          E|   H|
!     |  |-----------+-----------+----|     /
!     |  |3          |3          |1   |     | MDIFIJ
!     |  |          C|          F|   I|     |
!     |  *----------------------------*     /
!     |
!     |  <===========>           <====>
!     |   NDIFRS                  MDIFRS
!     |
!     /
!    ij-index
!
!     write (WTHCL1):  A->B->C---->D->E->F---->G->H->I
!                      <=====>     <=====>     <=====>
!                       on-memory
!
!     read  (WTHCL2):  G->D->A---->H->E->B---->I->F->C
!                      <=====>     <=====>     <=====>
!
!     Second order equation to solve the value of 'NDIFIJ'.
!         X = value of 'NDIFIJ' to be solved
!         C = memory size (counted by 8-Byte unit)
!         A = NMOMO = NMO*(NMO+1)/2
!         R = NRSEXT = existing rs-pair in SO integrals
!         Y = NTIMIJ = A/X <---- NHTTR
!         W = LR2LI (= 2) <---- because real*8 = 2*(4byte)
!
!         for NHTTR: Y/LR2LI
!         for DIJRS in WTHCL2: NRSEXT*NDIFIJ = R*X
!
!         C = Y/W + R*X
!         solve X from above equation
!         X = (C - Y/W) / R 
!         substitute Y by A/X,
!         X = (C - A/(X*W)) / R
!         accordingly,
!         W*R*X**2 - C*W*X + A = 0
!      -->
!             C*W +- SQRT( (C*W)**2 - 4*W*R*A )
!         X = ---------------------------------
!                      2*W*R
!
!----------------------------------------------------------------------
!  N.B.
!     LDOB74 == 791
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DD8BYT
!
!             ...declare constant(s) in code
!
      INTEGER :: NTRY   = 5
      INTEGER :: NFLCT  = 4
!
!             ...declare variable(s) in code
!
      INTEGER :: MA, MB 
      INTEGER :: K8CORE, MEM01, MEM02, MEM12
      INTEGER :: IGO1, IGO2
      INTEGER :: IJDF1, IJDF2
      REAL(KIND=LDREAL) :: A,B,C,D,SQ
!
!----------------------------------------------------------------------
 1000 FORMAT(
     &  1X,'OPTIMIZING CONTROL DATA FOR TRANSFORMATION',
     & /1X,'  NSO           =',I12,
     & /1X,'  NMO,NMOMO     =',2I12,
     & /1X,'  NPQINT,NRSEXT =',2I12,
     & /1X,'  NTIMIJ,NTIMRS =',2I12,
     & /1X,'  NDIFIJ,NDIFRS =',2I12)
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTHOPT start'
!
!     ...buffer for writing transformed MO integrals
!
      LDOB74 = (LDTRAK - LDINTG*2) / (LDINTG*4 + LDREAL)
!
!     ...calculate MEM01 
!        memory size needed for WTHCL1 except for NHTTR, DIJRS
!        KK8(1:5) has been already set in WTHDRV to call WTHRMO.
!
!"    KK8( 1) = NSO/LR2LI + N1PAD       ! IREPSO
!"    KK8( 2) = NMO/LR2LI + N1PAD       ! IREPMO
!"    KK8( 3) = NSO*2/LR2LI + N1PAD     ! S2MRNG
!"    KK8( 4) = NMO*2/LR2LI + N1PAD     ! M2SRNG
!"    KK8( 5) = NSO*NMO                 ! CSOMO (MO coef)
      KK8( 6) = NRSEXT*2/LR2LI + N1PAD  ! IRSPAC 
      KK8( 7) = LDOB95                  ! DBUFFH
      KK8( 8) = LDOB95*2/LR2LI + N1PAD  ! IBUFFH
      KK8( 9) = NSO*NMO                 ! WIQRS 
      KK8(10) = LDOB94                  ! DBUFFR
      KK8(11) = LDOB94*2/LR2LI + N1PAD  ! IBUFFR
      KK8(12) = MAX(NMOMO,NRSEXT)*NFLCT ! flucuation
      KKLAST = 12
      CALL WTMMPT
      MEM01 = KK(KKLAST+1)
!
!     ...calculate MEM02 
!        memory size needed for WTHCL2 except for NHTTR, DIJRS
!
!        IREPSO is not referred in WTHCL2, so KK8(1) is unnecessary,
!        but not deleted in order to keep KK8(2:8).
!
      KK8( 9) = NSO*NMO                 ! WIJKS
      KK8(10) = LDOB74                  ! DBUFFW
      KK8(11) = LDOB74*4/LR2LI + N1PAD  ! IBUFFW
      KK8(12) = NMOMO*2/LR2LI + N1PAD   ! IJREV
      KK8(13) = MAX(NMOMO,NRSEXT)*NFLCT ! flucuation
      KKLAST = 13
      CALL WTMMPT
      MEM02 = KK(KKLAST+1)
!
      MEM12 = MAX(MEM01,MEM02)
      K8CORE = DD8BYT - MEM12 
!
      A = NRSEXT*LR2LI
      B = K8CORE*LR2LI
      C = NMOMO
      D = B*B - 4.D0*A*C
      IF( D < 0.D0 ) THEN
        WRITE(*,*) 'LOGICAL ERROR',A,B,C,D
        STOP 'WTHOPT 001'
      END IF
      SQ = SQRT(D)
!
      MA = (B+SQ)/(A*2.D0)
      IF( B < SQ ) THEN
        MB = 0
      ELSE
        MB = (B-SQ) / (A*2.D0)
      END IF
      NDIFIJ = MAX(MA,MB)
      IF( NDIFIJ > NMOMO ) NDIFIJ = NMOMO
!
      NTIMIJ = NMOMO/NDIFIJ
      IF( MOD(NMOMO,NDIFIJ) /= 0 ) NTIMIJ = NTIMIJ + 1
!
!     ...for DIJRS in WTHCL1: NDIFRS*NMOMO
      NDIFRS = (K8CORE - NTIMIJ/LR2LI) / NMOMO
      IF( NDIFRS > NRSEXT ) NDIFRS = NRSEXT
!
      NTIMRS = NRSEXT/NDIFRS
      IF( MOD(NRSEXT,NDIFRS) /= 0 ) NTIMRS = NTIMRS + 1
      IF( NTIMIJ == NTIMRS ) THEN
        CALL WTQMEM(IGO1)
        IF( IGO1 == 1 ) GOTO 19
      END IF
!
      IJDF1 = MAX(NDIFIJ-NTRY,1)
      IJDF2 = MIN(NDIFIJ+NTRY,NMOMO)
!
      DO NDIFIJ = IJDF1, IJDF2
        CALL WTQOP2(IGO2)
        IF( IGO2 == 1 ) THEN
          CALL WTQMEM(IGO1)
          IF( IGO1 == 1 ) GOTO 19
        END IF
      END DO
!
!     ...failed in determining NDIFIJ and NDIFRS
!
!!      WRITE(*,*) NDIFIJ,NTIMIJ,NDIFRS,NTIMRS
      STOP 'WTHOPT 002'
!
   19 CONTINUE
!
!     ...get MDIFIJ, MDIFRS.
!
      CALL WTQDIF
!!      WRITE(*,1000) NSO,NMO,NMOMO,NPQINT,NRSEXT,
!!     *              NTIMIJ,NTIMRS,NDIFIJ,NDIFRS
!
!      WRITE(*,*) 'WTHOPT ended'
      RETURN
      END SUBROUTINE
