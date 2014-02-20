      SUBROUTINE WTOCL1(
     *                  ICANON, N1TTR,  ID1BOX, NCOUNT, PQRS, 
     *                  IPQRS,  DBUFF,  IBUFF
     *                  )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                is the 1-st stage of ordering SO basis integrals.
!                FT73 (sequential) -> FT91 (direct access)
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
      INTEGER,DIMENSION(4,LDOB73) :: IBUFF
      INTEGER,DIMENSION(4,LDOB91,NBOX) :: IPQRS
      INTEGER,DIMENSION(NBOX) :: N1TTR, NCOUNT
      INTEGER,DIMENSION(NSOSO) :: ID1BOX
      REAL(KIND=LDREAL),DIMENSION(LDOB73) :: DBUFF
      REAL(KIND=LDREAL),DIMENSION(LDOB91,NBOX) :: PQRS
!
!             ...declare variable(s) in code
!
      INTEGER :: IADRS1, I1BOX, ICOUNT, I1TTR, IADRS2
      INTEGER :: IW, KW, NW
      INTEGER :: IP, IQ, IR, IS
      INTEGER :: JP, JQ, JR, JS
      INTEGER :: IWREC
!
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WTOCL1 start'
!
!     ...zero clear
!
      IWREC = 0
      N1TTR = 0
      NCOUNT = 0
!
      OPEN(UNIT=NFT73,FILE=FN73,STATUS='OLD',
     *     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!
! .................
      LOOP100: DO
! .................
!
        READ(NFT73) KW,NW,DBUFF,IBUFF
        IF( KW == 0 ) EXIT
        DO IW = 1,NW
          IP = IBUFF(1,IW)
          IQ = IBUFF(2,IW)
          IR = IBUFF(3,IW)
          IS = IBUFF(4,IW)
          JP = MAX(IP,IQ)
          JQ = MIN(IP,IQ)
          JR = MAX(IR,IS)
          JS = MIN(IR,Is)
          IADRS1 = ICANON(JP) + JQ
          I1BOX = ID1BOX(IADRS1)
          NCOUNT(I1BOX) = NCOUNT(I1BOX) + 1
          ICOUNT = NCOUNT(I1BOX)
          PQRS(ICOUNT,I1BOX) = DBUFF(IW)
          IPQRS(1,ICOUNT,I1BOX) = JR
          IPQRS(2,ICOUNT,I1BOX) = JS
          IPQRS(3,ICOUNT,I1BOX) = JP
          IPQRS(4,ICOUNT,I1BOX) = JQ
!
          IF( ICOUNT == LDOB91 ) THEN
            IWREC = IWREC + 1
            I1TTR = N1TTR(I1BOX)
            N1TTR(I1BOX) = IWREC
            CALL WTOWB1(LDOB91,NFT91,IWREC,I1TTR,ICOUNT,
     *                  PQRS(1,I1BOX),IPQRS(1,1,I1BOX))
            NCOUNT(I1BOX) = 0
          END IF
!
          IADRS2 = ICANON(JR) + JS
          IF( IADRS2 /= IADRS1 ) THEN
            I1BOX = ID1BOX(IADRS2)
            NCOUNT(I1BOX) = NCOUNT(I1BOX) + 1
            ICOUNT = NCOUNT(I1BOX)
            PQRS(ICOUNT,I1BOX) = DBUFF(IW)
            IPQRS(1,ICOUNT,I1BOX) = JP
            IPQRS(2,ICOUNT,I1BOX) = JQ
            IPQRS(3,ICOUNT,I1BOX) = JR
            IPQRS(4,ICOUNT,I1BOX) = JS
            IF(ICOUNT == LDOB91) THEN
              IWREC = IWREC+1
              I1TTR = N1TTR(I1BOX)
              N1TTR(I1BOX) = IWREC
              CALL WTOWB1(LDOB91,NFT91,IWREC,I1TTR,ICOUNT,
     *                    PQRS(1,I1BOX),IPQRS(1,1,I1BOX))
              NCOUNT(I1BOX) = 0
            END IF
          END IF
        END DO
        IF( KW < 0 ) EXIT
!
! ....................
      END DO LOOP100
! ....................
!
!     ...NFT73 not referred anymore
!
      CLOSE(NFT73) 
!
!     ...ending 1st ordering
!
      DO I1BOX = 1, NBOX
        ICOUNT = NCOUNT(I1BOX)
        IF( ICOUNT > 0 ) THEN
          IWREC = IWREC + 1
          I1TTR = N1TTR(I1BOX)
          N1TTR(I1BOX) = IWREC
          CALL WTOWB1(LDOB91,NFT91,IWREC,I1TTR,ICOUNT,
     *                PQRS(1,I1BOX),IPQRS(1,1,I1BOX))
        END IF
      END DO
!
      NWTR91 = IWREC
!
!      WRITE(*,*) 'WTOCL1 ended'
      RETURN
      END SUBROUTINE
