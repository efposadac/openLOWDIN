      SUBROUTINE WTQOP2(
     *                  IGO2
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                optimizes NDIFRS.  
!                If it succeeds, get NTIMIJ,NTIMRS.
!                If OK -> IGO2=1, NO-GOOD -> IGO2=-1.
!  RELEASE :     v.00  gen = sya 1987-06-24 at ims
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTQOPT
!----------------------------------------------------------------------
!  inequality
!     N: total number of object
!     x: number of object in a group
!     t: number of groups
!
!     1 <= N - (t-1)*x <= N
!     1-N <= - (t-1)*x <= 0
!     N-1 >= (t-1)*x >= 0     (because t>1, if t=1 then NDIFRS=NRSEXT)
!     (N-1)/(t-1) >= x >= 0 
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(INOUT) :: IGO2
!
!             ...declare constant(s) in code
!
      INTEGER :: NTRY = 10
!
!             ...declare variable(s) in code
!
      INTEGER :: K,K1,K2
!
!----------------------------------------------------------------------
!
      IGO2 = -1
      NTIMIJ = NMOMO / NDIFIJ
      IF( MOD(NMOMO,NDIFIJ) /= 0 ) NTIMIJ = NTIMIJ + 1
      NTIMRS = NTIMIJ
      IF( MOD(NRSEXT,NTIMRS) == 0 ) THEN
        NDIFRS = NRSEXT/NTIMRS
        IGO2 = 1
      ELSE
        NDIFRS = (NRSEXT - 1) / (NTIMRS - 1)
        K = NRSEXT - NDIFRS*(NTIMRS-1)
        IF( K > 0 .AND. K < NDIFRS ) THEN
          IGO2 = 1
        ELSE
          K1 = MAX(NDIFRS-NTRY,1)
          K2 = NDIFRS
          DO NDIFRS = K2,K1,-1
            K = NRSEXT - NDIFRS*(NTIMRS-1)
            IF( K > 0 .AND. K < NDIFRS ) THEN
              IGO2 = 1
              EXIT
            END IF
          END DO
        END IF
      END IF
!
!     ...recalculate NTIMRS
!
      NTIMRS = NRSEXT / NDIFRS
      IF( MOD(NRSEXT,NDIFRS) /= 0 ) NTIMRS = NTIMRS + 1
!
      RETURN
      END SUBROUTINE
