      SUBROUTINE WMABEL
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                gets NABEL and NSYM.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!----------------------------------------------------------------------
!     NSYM : order of the input Abel point group
!     NABEL : ID number of the input Abel point group
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!     ...constant array(s)
!
!                    ...name of Abel group
!
      CHARACTER(LEN=4),DIMENSION(LDABEL) :: ABEL 
     * = (/'C1  ','CI  ','CS  ','C2  ','C2V ','C2H ','D2  ','D2H '/)
!
!                    ...dimension of Abel group
!
      INTEGER,DIMENSION(LDABEL) :: DMABEL
     * = (/1,2,2,2,4,4,4,8/)
!
!                    ...irrep name
!
      CHARACTER(LEN=4),DIMENSION(LDSYM,LDABEL) :: IREPNM 
     * =RESHAPE( (/
     * 'A   ','    ','    ','    ','    ','    ','    ','    ',
     * 'AG  ','AU  ','    ','    ','    ','    ','    ','    ',
     * 'A*  ','A** ','    ','    ','    ','    ','    ','    ',
     * 'A   ','B   ','    ','    ','    ','    ','    ','    ',
     * 'A1  ','A2  ','B1  ','B2  ','    ','    ','    ','    ',
     * 'AG  ','BG  ','AU  ','BU  ','    ','    ','    ','    ',
     * 'A   ','B1  ','B2  ','B3  ','    ','    ','    ','    ',
     * 'AG  ','B1G ','B2G ','B3G ',
     * 'AU  ','B1U ','B2U ','B3U '/),
     * (/LDSYM,LDABEL/) )
!
!             ...define variable(s) in code
!
      INTEGER :: IABEL, ISYM
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMABEL start'
!
      DO IABEL = 1, LDABEL
        IF( NAMEPT == ABEL(IABEL) ) THEN
          NABEL = IABEL
          GOTO 109
        END IF
      END DO
!      WRITE(*,*) 'Point group not Abel ', NAMEPT
      STOP 'WMABEL 001'
!
  109 CONTINUE
      NSYM = DMABEL(NABEL)
!
!      WRITE(*,"(' Dimension of group: ',I1,' for ',A4)") NSYM,NAMEPT
!      WRITE(*,"(' Name of each irrep:')")
!      WRITE(*,FMT='(8(1X,A4))') (IREPNM(ISYM,NABEL),ISYM=1,NSYM)
!
!      WRITE(*,*) 'WMABEL ended'
      RETURN
      END SUBROUTINE
