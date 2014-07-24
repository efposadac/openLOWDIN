      SUBROUTINE WMINIT
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                initialized control data
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   PREPAR
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMINIT start'
!
      LDOB72 = (LDTRAK - LDINTG*2) / (LDINTG*2 + LDREAL)
      LDOB73 = (LDTRAK - LDINTG*2) / (LDINTG*4 + LDREAL)
!
!      WRITE(*,*) 'WMINIT ended'
!
      RETURN
      END SUBROUTINE
