      SUBROUTINE WTHSTB(
     *                  IREPSO, IREPMO, S2MRNG, M2SRNG, SOPAIR, 
     *                  MOPAIR
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets tables (S2MRNG, M2SRNG, etc.) for WTH*.
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHDRV
!----------------------------------------------------------------------
!  Arrays
!     MOPAIR(j,i) : active MO range of i-th symmetry
!     SOPAIR(j,i) : SO range of i-th symmetry
!
!  Variable
!     SENTNL is a sentinel to skip multiplication operation in WTHCL1
!     and WTHCL2 when there is no active orbitals for symmetry.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSO),INTENT(OUT) :: IREPSO
      INTEGER,DIMENSION(NMO),INTENT(OUT) :: IREPMO
      INTEGER,DIMENSION(2,NSO),INTENT(OUT) :: S2MRNG
      INTEGER,DIMENSION(2,NMO),INTENT(OUT) :: M2SRNG
      INTEGER,DIMENSION(2,NSYM) :: SOPAIR, MOPAIR
!
!             ...declare variable(s) in code
!
      INTEGER :: ISYM, ISO, IMO
      INTEGER :: SPRE, MPRE
!
!             ...declare constant(s) in code
!
      INTEGER :: SENTNL = -1
!
!----------------------------------------------------------------------
!      WRITE(*,*) 'WTHSTB start'
!
!     ...SO range for symmetry
!
      SOPAIR = 0
      SPRE = 0
      DO ISYM = 1, NSYM
        IF( NSOG(ISYM) == 0 ) CYCLE
        SOPAIR(1,ISYM) = SPRE + 1
        SOPAIR(2,ISYM) = SOPAIR(1,ISYM) + NSOG(ISYM) - 1
        SPRE = SOPAIR(2,ISYM)
      END DO
!
!     ...MO range for symmetry
!
      MOPAIR = 0
      MPRE = 0
      DO ISYM = 1, NSYM
        IF( NMOG(ISYM) == 0 ) CYCLE
        MOPAIR(1,ISYM) = MPRE + 1
        MOPAIR(2,ISYM) = MOPAIR(1,ISYM) + NMOG(ISYM) - 1
        MPRE = MOPAIR(2,ISYM)
      END DO
!
!     ...S2M range for symmetry
!
      S2MRNG = 0
      IREPSO = 0
      DO ISYM = 1, NSYM
        IF( NSOG(ISYM) == 0 ) CYCLE
        DO ISO = SOPAIR(1,ISYM), SOPAIR(2,ISYM)
          S2MRNG(1,ISO) = MOPAIR(1,ISYM) 
          IF( MOPAIR(1,ISYM) == 0 ) THEN
            S2MRNG(2,ISO) = SENTNL
          ELSE 
            S2MRNG(2,ISO) = MOPAIR(2,ISYM) 
          END IF
          IREPSO(ISO) = ISYM
        END DO
      END DO
!
!     ...M2S range for symmetry
!
      M2SRNG = 0
      IREPMO = 0
      DO ISYM = 1, NSYM
        IF( NMOG(ISYM) == 0 ) CYCLE
        DO IMO = MOPAIR(1,ISYM), MOPAIR(2,ISYM)
          M2SRNG(1,IMO) = SOPAIR(1,ISYM)
          IF( SOPAIR(1,ISYM) == 0 ) THEN
            M2SRNG(2,IMO) = SENTNL
          ELSE 
            M2SRNG(2,IMO) = SOPAIR(2,ISYM)
          END IF
          IREPMO(IMO) = ISYM
        END DO
      END DO
!
!      WRITE(*,*) 'WTHSTB ended'
      RETURN
      END SUBROUTINE
