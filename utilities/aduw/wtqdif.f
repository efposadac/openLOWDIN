      SUBROUTINE WTQDIF
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                calculates MDIFIJ, MDIFRS.
!  RELEASE :     v.00  gen = sya 1999-01-11 at chukyo-u
!                v.01  mod = sya 2003-04-05 at chukyo-u
!  CALLED BY :   WTHOPT
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'sym4tr.h'
!----------------------------------------------------------------------
      KTIMIJ = NMOMO/NDIFIJ
      MDIFIJ = NMOMO - NDIFIJ*KTIMIJ
!
      KTIMRS = NRSEXT/NDIFRS
      MDIFRS = NRSEXT - NDIFRS*KTIMRS
!
      RETURN
      END SUBROUTINE
