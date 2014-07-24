      SUBROUTINE WTOWB2(
     *                  DDOB92,NFT,IWREC,
     *                  ITTR,ICOUNT,
     *                  IR,IS,PQRS,IPQ
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                writes integrals and pq-indices in bin-2 (work direct 
!                access file).
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTOCL2
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DDOB92
      INTEGER,DIMENSION(2,DDOB92),INTENT(IN) :: IPQ
      REAL(KIND=LDREAL),DIMENSION(DDOB92),INTENT(IN) :: PQRS
      INTEGER,INTENT(IN) :: NFT, IWREC, ITTR, ICOUNT
      INTEGER,INTENT(IN) :: IR, IS
!
!             ...declare constant(s) in code
!
      INTEGER :: IPAD = 0
!----------------------------------------------------------------------
      WRITE(NFT,REC=IWREC) ITTR,ICOUNT,DDOB92,IPAD,IR,IS,PQRS,IPQ
      RETURN
      END SUBROUTINE
