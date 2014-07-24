      SUBROUTINE WTOWB1(
     *                  DDOB91,NFT,IWREC,
     *                  ITTR,ICOUNT,PQRS,IPQRS
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                writes integrals and pq-indices in bin-1 (work direct 
!                access file).
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTOCL1
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DDOB91
      INTEGER,INTENT(IN) :: NFT,IWREC,ITTR,ICOUNT
      REAL(KIND=LDREAL),DIMENSION(DDOB91),INTENT(IN) :: PQRS
      INTEGER,DIMENSION(4,DDOB91),INTENT(IN) :: IPQRS
!----------------------------------------------------------------------
      WRITE(NFT,REC=IWREC) ITTR,ICOUNT,PQRS,IPQRS
      RETURN
      END SUBROUTINE
