      SUBROUTINE WTSDRV(
     *                  DD8BYT,
     *                  AA
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                is main driver for merging bin-2 to a sequential 
!                file.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   SYM4TR
!---------------------------------------------------------------------
!  WTSDRV-+-WTSCTL   read control data from FT93 to make WTS* module
!         |          independent from WTO* module.
!         +-SETDRW
!         +-WTSCAL   read integrals from FT92 and write them to FT94
!----------------------------------------------------------------------
!  I/O
!     input:  NFT92(da), NFT93(seq)
!     output: NFT94(seq)
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DD8BYT
      REAL(KIND=LDREAL),DIMENSION(DD8BYT) :: AA
!
!             ...declare variable(s) in code
!
      INTEGER :: IERROR, IOSTAT
!
!----------------------------------------------------------------------
!      WRITE(*,*) 'WTSDRV start'
!
      REWIND NFT93
      OPEN(UNIT=NFT94, STATUS='UNKNOWN',
     *     ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
!
      CALL SETDRW(NFT92, LREC92, IERROR, IOSTAT)
!
      CALL WTSCTL(NFT93, NBXDIF, LDOB92, LREC92)
      LDOB94 = LDOB92 + NSOSO
!
      KK8(1) = NBXDIF / LR2LI + N1PAD     ! N2TTR
      KK8(2) = LDOB94                     ! PQRS
      KK8(3) = LDOB94 * 2 / LR2LI + N1PAD ! IPQ
      KKLAST = 3
      CALL WTMMPT
      IF( KK(KKLAST+1) > DD8BYT ) GOTO 900
!
      CALL WTSCAL(
     &  AA(KK(1)),AA(KK(2)),AA(KK(3)) )
!       N2TTR     PQRS      IPQ
!
      IF( NRFLAG == 0 ) THEN
        CLOSE(NFT92,STATUS='DELETE')
        CLOSE(NFT93,STATUS='DELETE')
      ELSE
        CLOSE(NFT92)
        CLOSE(NFT93)
      END IF
!
!      WRITE(*,*) 'WTSDRV ended'
      RETURN
!
  900 CONTINUE
      WRITE(*,*) "Four-Index transformation library:"
      WRITE(*,*) "Memory insufficient"
      WRITE(*,*) DD8BYT,KK(1:KKLAST+1)
      STOP 'WTSDRV 001'
      END SUBROUTINE
