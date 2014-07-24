      SUBROUTINE WTORB1(DDOB91,NFT,ITTR,MX,PQRS,IPQRS)
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads data from bin-2 in 2nd stage of SO integral
!                ordering
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DDOB91,NFT
      INTEGER,INTENT(INOUT) :: ITTR
      INTEGER,INTENT(OUT) :: MX
      REAL(KIND=LDREAL),DIMENSION(DDOB91),INTENT(OUT) :: PQRS
      INTEGER,DIMENSION(4,DDOB91),INTENT(OUT) :: IPQRS
!
!             ...declare variable(s) in code
!
      INTEGER :: K1, K2
!
!---------------------------------------------------------------------
      READ(NFT,REC=ITTR) K1,K2,PQRS,IPQRS
      ITTR = K1
      MX = K2
      RETURN
      END SUBROUTINE
