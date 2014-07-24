      SUBROUTINE WTOMP2(
     *                  DD8BYT
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets memory partitioning data for 2nd stage of 
!                ordering SO integrals.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare variable(s) in code
!
      INTEGER :: DD8BYT
!
!             ...declare function(s)
!
      INTRINSIC :: SUM
!
!---------------------------------------------------------------------
!
!      WRITE(*,*) 'WTOMP2 ended'
!
!"    KK8( 1) = NSO / LR2LI + N1PAD        ! ICANON (same as WTOMP1)
!"    KK8( 2) = NBOX / LR2LI + N1PAD       ! N1TTR  (same as WTOMP1)
      KK8( 3) = NSOSO / LR2LI + N1PAD      ! ID2BOX
      KK8( 4) = NBXDIF / LR2LI + N1PAD     ! N2TTR
      KK8( 5) = NBXDIF / LR2LI + N1PAD     ! NCOUNT
      KK8( 6) = LDOB91                     ! DBUFF
      KK8( 7) = LDOB91 * 4 / LR2LI + N1PAD ! IBUFF
      KK8( 8) = NBXDIF * 2 / LR2LI + N1PAD ! IRS
      M8CORE = DD8BYT - SUM(KK8(1:8))
!
!     ...make buffer length less than 1 track
!
      LDOB92 = M8CORE / (NBXDIF + NBXDIF*2/LR2LI)
!
      DO
        LREC92 = 6 * LDINTG + LDOB92 * LDREAL + 2 * LDOB92 * LDINTG
        IF( LREC92 < =  LDTRAK ) EXIT
        LDOB92 = LDOB92 - 1
      END DO
!
!     ...memory partitioning
!
      KK8( 9) = LDOB92 * NBXDIF                     ! PQRS
      KK8(10) = LDOB92 * NBXDIF * 2 / LR2LI + N1PAD ! IPQ
      KKLAST  = 10
      CALL WTMMPT
      IF( KK(KKLAST+1) <= DD8BYT ) THEN
!        WRITE(*,*) 'WTOMP2 ended'
        RETURN
      END IF
!
      WRITE(*,*) "Four-Index transformation library:"
      WRITE(*,*) "Memory insufficient"
      STOP 'WTOMP2 001'
      END SUBROUTINE
