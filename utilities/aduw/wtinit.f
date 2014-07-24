      SUBROUTINE WTINIT
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads inf.dat file and initializes control data.
!  RELEASE :     v.01  gen = sya 2003-03-30 at chukyo-u
!  CALLED BY :   SYM4TR
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTINIT start'
!
      OPEN(UNIT=NFT71, FILE=FN71, STATUS='OLD',
     *     ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      READ(NFT71) 
     *  LDOB72,LDOB73,NSYM,NSOG,NMOG,NSO,NMO,NPFLAG,NRFLAG
      READ(NFT71) 
     *  THRSOI,THRHLF,THRMOI
      READ(NFT71) 
     *  TITLEC,NAMEPT
      CLOSE(NFT71)
!
!      WRITE(*,*) 'WTINIT ended'
!
      RETURN
      END SUBROUTINE
