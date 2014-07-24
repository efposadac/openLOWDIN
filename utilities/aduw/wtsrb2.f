      SUBROUTINE WTSRB2(
     *                  ITTR,  KPQINT,
     *                  IR,    IS,     PQRS,  IPQ
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads data from bin-2 and mergers them.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTSDRV
!----------------------------------------------------------------------
!     FT92 direct-access file (bin-2):                 
!       a record:                      <---A--->
!       another record:                <---B--->
!                          
!     FT94 sequential file:
!       a record:                      <---A---><---B---> 
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(INOUT) :: ITTR, KPQINT
      INTEGER,INTENT(OUT) :: IR, IS
      REAL(KIND=LDREAL),DIMENSION(LDOB94),INTENT(INOUT) :: PQRS
      INTEGER,DIMENSION(2,LDOB94),INTENT(INOUT) :: IPQ
!
!             ...declare variable(s) in code
!
      INTEGER :: KTTR, MPQ, KOBJ, IPAD
!
!---------------------------------------------------------------------
!
      READ(NFT92,REC=ITTR) KTTR,MPQ,KOBJ,IPAD,IR,IS,
     *  PQRS(1:LDOB92),IPQ(:,1:LDOB92)
      ITTR = KTTR
      KPQINT = KPQINT + MPQ
!
      RETURN
      END SUBROUTINE
