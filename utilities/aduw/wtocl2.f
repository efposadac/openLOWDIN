      SUBROUTINE WTOCL2(
     *                  ICANON, N1TTR,  ID2BOX, N2TTR,  NCOUNT,
     *                  DBUFF,  IBUFF,  IRS,    PQRS,   IPQ
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                is 2nd stage of integral ordering.
!  RELEASE :     v.00  gen = una
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WTODRV
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'sym4tr.h'
!
!             ...declare argument(s) in code
!
      INTEGER,DIMENSION(NSO) :: ICANON
      INTEGER,DIMENSION(NBOX) :: N1TTR
      INTEGER,DIMENSION(NSOSO) :: ID2BOX
      INTEGER,DIMENSION(NBXDIF) :: N2TTR, NCOUNT
      REAL(KIND=LDREAL),DIMENSION(LDOB91) :: DBUFF
      INTEGER,DIMENSION(4,LDOB91) :: IBUFF
      INTEGER,DIMENSION(2,NBXDIF) :: IRS
      REAL(KIND=LDREAL),DIMENSION(LDOB92,NBXDIF) :: PQRS
      INTEGER,DIMENSION(2,LDOB92,NBXDIF) :: IPQ
!
!             ...declare variable(s) in code
!
      INTEGER :: IW, MW, IR, IS
      INTEGER :: I1TTR, I2TTR, I1BOX, I2BOX
      INTEGER :: ICOUNT, IADRS, IWREC
!
!---------------------------------------------------------------------
!
!      WRITE(*,*) 'WTOCL2 start'
!
      IWREC = 0
!
      WRITE(NFT93) NBXDIF,LDOB92,LREC92
!
      DO I1BOX = 1, NBOX
        N2TTR = 0
        NCOUNT = 0
        I1TTR = N1TTR(I1BOX)
        IF( I1TTR == 0 ) CYCLE
!
        DO
          CALL WTORB1(LDOB91,NFT91,I1TTR,MW,DBUFF,IBUFF)
          DO IW = 1, MW
            IR = IBUFF(3,IW)
            IS = IBUFF(4,IW)
            IADRS = ICANON(IR) + IS
            I2BOX = ID2BOX(IADRS)
            NCOUNT(I2BOX) = NCOUNT(I2BOX) + 1
            ICOUNT = NCOUNT(I2BOX)
            PQRS(ICOUNT,I2BOX) = DBUFF(IW)
            IPQ(1,ICOUNT,I2BOX) = IBUFF(1,IW)
            IPQ(2,ICOUNT,I2BOX) = IBUFF(2,IW)
            IRS(1,I2BOX) = IR
            IRS(2,I2BOX) = IS
            IF( ICOUNT == LDOB92 ) THEN
              IWREC = IWREC + 1
              I2TTR = N2TTR(I2BOX)
              N2TTR(I2BOX) = IWREC
              CALL WTOWB2(LDOB92,NFT92,IWREC,I2TTR,ICOUNT,
     *                    IR,IS,PQRS(1,I2BOX),IPQ(1,1,I2BOX) )
              NCOUNT(I2BOX) = 0
            END IF
          END DO
          IF( I1TTR == 0 ) EXIT
        END DO
!
!     ---- concluding 2nd ordering ----
!
        DO I2BOX = 1, NBXDIF
          IF( NCOUNT(I2BOX) > 0 ) THEN
            ICOUNT = NCOUNT(I2BOX)
            IWREC = IWREC + 1
            I2TTR = N2TTR(I2BOX)
            N2TTR(I2BOX) = IWREC
            IR = IRS(1,I2BOX)
            IS = IRS(2,I2BOX)
            CALL WTOWB2(LDOB92,NFT92,IWREC,I2TTR,ICOUNT,
     *                  IR,IS,PQRS(1,I2BOX),IPQ(1,1,I2BOX) )
          END IF
        END DO
!
        WRITE(NFT93) NBXDIF,N2TTR
!
      END DO
!
      IF( NRFLAG == 0 ) THEN
        CLOSE(NFT91,STATUS='DELETE') 
      ELSE
        CLOSE(NFT91) 
      END IF
      CLOSE(NFT92) ! EOF attached
      CLOSE(NFT93) ! EOF attached
      NWTR92 = IWREC
!
!      WRITE(*,*) 'WTOCL2 ended'
      RETURN
      END SUBROUTINE
