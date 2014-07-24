      SUBROUTINE WTOMP1(
     *                  DD8BYT,
     *                  AA
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                sets control data for ordering SO basis integrals
!  RELEASE :     v.00  gen  =  una
!                v.01  mod  =  sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,INTENT(IN) :: DD8BYT
      REAL(KIND=LDREAL),DIMENSION(DD8BYT),INTENT(INOUT) :: AA
!
!             ...declare variable(s) in code
!
      INTEGER :: IMODE
!
!---------------------------------------------------------------------
!
!      WRITE(*,*) 'WTOMP1 start'
!
!     ...set control data
!
      NSOSO = (NSO * (NSO + 1))/2
!
!     ...initialize by maximum value
!
      LDOB91 = (LDTRAk - 2*LDINTG)/(LDREAL + 4*LDINTG)
      LDOB92 = (LDTRAk - 4*LDINTG)/(LDREAL + 2*LDINTG)
      LREC91 = 2*LDINTG + (LDREAL + 4*LDINTG)*LDOB91
      LREC92 = 4*LDINTG + (LDREAL + 2*LDINTG)*LDOB92
!
      KK8(1) = NSO / LR2LI + N1PAD           ! ICANON, 
!                                            ! referred by WTOMP2 also
      KK8(7) = LDOB73                        ! DBUFF(LDOB73)
      KK8(8) = LDOB73 * (4 / LR2LI) + N1PAD  ! IBUFF(4,LDOB73)
      M8CORE = DD8BYT - (KK8(1) + KK8(7) + KK8(8))
!
      IMODE = 2
!
!     ...determine LDOB91,NBOX,NBXDIF for 2-pass sort
!
      DO
        NBOX = M8CORE * LDREAL / (LDOB91*(LDREAL + LDINTG*4))
        NBXDIF = NSOSO / NBOX
        IF( NBXDIF > NBOX ) THEN
          LDOB91 = LDOB91 - 1
          CYCLE
        ENDIF
!
!     ...check possibility of 1-pass sort
!
        IF( NSOSO <= NBOX ) THEN
          IMODE = 1
          NBOX = NSOSO
          NBXDIF = 1
        END IF
!
!     ...memory partitioning
!
        KK8(2) = NBOX / LR2LI + N1PAD        ! N1TTR, 
!                                            ! referred by WTOMP2 also
        KK8(3) = NSOSO / LR2LI + N1PAD       ! ID1BOX
        KK8(4) = NBOX / LR2LI + N1PAD        ! NCOUNT
        KK8(5) = LDOB91 * NBOX               ! PQRS(LDOB91)
        KK8(6) = LDOB91 * NBOX * (4 / LR2LI) ! IPQRS(4,LDOB91)
        KKLAST = 8
        CALL WTMMPT
        IF( KK(KKLAST+1) < DD8BYT ) EXIT
!
ccc     WRITE(*,*) 'try again by decresing LDOB91'
        LDOB91 = LDOB91 - 1
      END DO
!
!      IF( IMODE == 1 ) THEN
!        WRITE(*,*) '1-pass sort is possible, but 2-pass'
!      END IF
!      IF( NPFLAG > 0 ) THEN
!        WRITE(*,*) 'LDOB91:',LDOB91
!        WRITE(*,*) 'NBOX:  ',NBOX
!        WRITE(*,*) 'NBXDIF:',NBXDIF
!      ENDIF
!
      CALL VUTCN2(NSO,AA(KK(1)))
!
!      WRITE(*,*) 'WTOMP1 ended'
      RETURN
      END SUBROUTINE
