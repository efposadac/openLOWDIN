      SUBROUTINE PRCPUT(MESSAG,INITLZ)
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                is a general utility to display message, time-zone, 
!                clock, cpu-time.
!  RELEASE :     v.01  gen = S.Yamamoto 2003-04-2 at chukyo-u
!----------------------------------------------------------------------
!  Arguments and variables
!     MESSAG : message to be printed
!     INITLZ : if 1, initialize cpu-time.
!----------------------------------------------------------------------
!  N.B.
!     CPU_TIME is an intrinsic subroutine of Fortran 95 standard.
!     DATE_AND_TIME is an intrinsic subroutine of Fortran 90 standard.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
!
!             ...declare argument(s) in code
!
      CHARACTER(LEN=*) :: MESSAG
      INTEGER :: INITLZ
!
!             ...declare constant(s) in code
!
      REAL(KIND=LDREAL) :: SIXTY = 60.D0
!
!             ...declare variable(s) in code
!
      CHARACTER(LEN=32) :: OMESSG
      CHARACTER(LEN=8)  :: CDATE
      CHARACTER(LEN=10) :: ODATE
      CHARACTER(LEN=10) :: CTIME
      CHARACTER(LEN=8)  :: OTIME
      CHARACTER(LEN=5)  :: CZONE
      CHARACTER(LEN=3)  :: OZONE
      REAL :: RCPUTM
      REAL(KIND=LDREAL) :: CPUTM0, CPUTM1, CPUDIF, CPUMIN
!
!             ...declare common 
!
      COMMON/CPUTIM/CPUTM0
!
!             ...declare function(s) in code
!
      INTRINSIC CPU_TIME
      INTRINSIC DATE_AND_TIME
!
!----------------------------------------------------------------------
!
      OMESSG = '          '//'          '//'          '//'  '
      OMESSG = MESSAG
!
      CALL DATE_AND_TIME(DATE=CDATE,TIME=CTIME,ZONE=CZONE)
      OZONE(1:3) = CZONE(1:3)
      ODATE(1:4) = CDATE(1:4)
      ODATE(5:5) = '-'
      ODATE(6:7) = CDATE(5:6)
      ODATE(8:8) = '-'
      ODATE(9:10)= CDATE(7:8)
      OTIME(1:2) = CTIME(1:2)
      OTIME(3:3) = ':'
      OTIME(4:5) = CTIME(3:4)
      OTIME(6:6) = ':'
      OTIME(7:8) = CTIME(5:6)
!
      IF( INITLZ == 1 ) THEN
        CPUTM0 = 0.0D0
      END IF
      CALL CPU_TIME(RCPUTM)
      CPUTM1 = RCPUTM
      CPUDIF = CPUTM1 - CPUTM0
      CPUMIN = CPUDIF / SIXTY
!
!      WRITE(*,FMT='(1X,A32," time-zone=",A3," date=",A10," time=",A8)') 
!     *  OMESSG,OZONE,ODATE,OTIME
      WRITE(*,FMT='(1X,"CPU-time=",F20.2,"sec=",F15.2,"min")') 
     *  CPUDIF,CPUMIN
!
      RETURN
      END SUBROUTINE
