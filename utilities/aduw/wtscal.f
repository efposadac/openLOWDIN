      SUBROUTINE WTSCAL(
     *                  N2TTR,PQRS,IPQ
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads ordered SO integrals from direct access file
!                (bin-2) and write them to sequential file.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTSDRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NBXDIF) :: N2TTR
      INTEGER,DIMENSION(2,LDOB94) :: IPQ
      REAL(KIND=LDREAL),DIMENSION(LDOB94) :: PQRS
!
!             ...declare constant(s) in code
!
      INTEGER :: IPAD = 0
!
!             ...declare variable(s) in code
!
      INTEGER :: MBOX, IPAIR, KPQ1, KPQINT
      INTEGER :: ITTR
      INTEGER :: IR, IS
!
!---------------------------------------------------------------------
!
!      WRITE(*,*) 'WTSCAL start'
!
      NPQINT = 0
      NRSEXT = 0
!
      DO
        READ(NFT93,END=199) MBOX,N2TTR(1:MBOX)
!
        DO IPAIR = 1, MBOX 
          ITTR = N2TTR(IPAIR)
          IF( ITTR == 0 ) CYCLE
          KPQINT = 0
          DO 
            KPQ1 = KPQINT + 1
            CALL WTSRB2(ITTR,KPQINT,IR,IS,
     *                  PQRS(KPQ1),IPQ(1,KPQ1)  )
            IF(ITTR == 0) EXIT
          END DO
          IF( KPQINT > NSOSO ) THEN
!            WRITE(*,*) 'KPQINT .GT. NSOSO',KPQINT,NSOSO
            STOP 'WTSCAL 001'
          END IF
          NPQINT = MAX(NPQINT,KPQINT)
          NRSEXT = NRSEXT + 1
          WRITE(NFT94) KPQINT,IR,IS,IPAD,
     *                 PQRS(1:KPQINT),IPQ(:,1:KPQINT)
        END DO 
      END DO 
!
  199 CONTINUE
!      WRITE(*,*) 'NPQINT:',NPQINT
!      WRITE(*,*) 'NRSEXT:',NRSEXT
!      WRITE(*,*) 'WTSCAL ended'
      RETURN
      END SUBROUTINE
