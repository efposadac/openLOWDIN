      SUBROUTINE WTSCTL(
     *                  NFT,NBXDIF,LDOB92,LREC92
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads control data for integral ordering and check
!                them with COMMON data.
!                This is to make WTS* module independent from WTO* 
!                module.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTSDRV
!----------------------------------------------------------------------
      IMPLICIT NONE
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: NFT, NBXDIF, LDOB92, LREC92
!
!             ...declare variable(s) in code
!
      INTEGER :: IBXDIF, IOBJ92, IREC92
      INTEGER :: NG
!
!---------------------------------------------------------------------
!
      NG = 0
      READ(NFT) IBXDIF,IOBJ92,IREC92
!
      IF( IBXDIF /= NBXDIF ) THEN
        WRITE(*,*) 'IBXDIF .ne. NBXDIF',IBXDIF,NBXDIF
        NG = NG + 1
      END IF
!
      IF( IOBJ92 /= LDOB92 ) THEN
        WRITE(*,*) 'IOBJ92 .ne. LDOB92',IOBJ92,LDOB92
        NG = NG + 1
      END IF
!
      IF( IREC92 /= LREC92 ) THEN
        WRITE(*,*) 'IREC92 .ne. LREC92',IREC92,LREC92
        NG = NG + 1
      END IF
!
      IF( NG > 0 ) THEN
        STOP 'WTSCTL 001'
      END IF
!
      RETURN
      END SUBROUTINE
