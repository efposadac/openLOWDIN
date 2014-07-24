      SUBROUTINE WMINFO
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                check consistency of the input data and save them.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   PREPAR
!----------------------------------------------------------------------
!  Variables
!     LDOB72 == 1187
!     LDOB73 == 791
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!             ...define constant(s) in code
!
      REAL(KIND=LDREAL) :: SMALL = 1.0D-10
!
!             ...define variable(s) in code
!
      INTEGER :: ISYM
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMINFO start'
!
      IF( NMO > NSO ) STOP 'WMINFO 001'
!
      DO ISYM = 1, NSYM
        IF( NMOG(ISYM) > NSOG(ISYM) ) STOP 'WMINFO 011'
      END DO
!
      IF( THRSOI < SMALL ) THRSOI = SMALL
      IF( THRHLF < SMALL ) THRHLF = SMALL
      IF( THRMOI < SMALL ) THRMOI = SMALL
!
!     ...save information in a binary file
!
      OPEN(UNIT=NFT71, FILE=FN71, STATUS='UNKNOWN',
     *     ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      WRITE(NFT71) 
     *  LDOB72,LDOB73,NSYM,NSOG,NMOG,NSO,NMO,NPFLAG,NRFLAG
      WRITE(NFT71) 
     *  THRSOI,THRHLF,THRMOI
      WRITE(NFT71) 
     *  TITLEC,NAMEPT
      CLOSE(NFT71)
!
!      WRITE(*,*) 'WMINFO ended'
!
      RETURN
      END SUBROUTINE
