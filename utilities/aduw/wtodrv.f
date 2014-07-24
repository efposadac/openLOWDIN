      SUBROUTINE WTODRV(
     *                  DD8BYT,
     *                  AA
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                orders SO basis integrals by 2-pass 2-bin sort.
!  RELEASE :     v.00  gen = Umpei Nagashima (una)
!                v.01  mod = sya 2003-03-29 at chukyo-u
!  CALLED BY :   SYM4TR
!----------------------------------------------------------------------
!  WTODRV-+-WTOMP1           SET NBOX, MEMORY PARTITIONING
!         +-WTOTB1           SET TABLES
!         +-WTOCL1---WTOWT1  1ST STAGE MAIN CALCULATION
!         !
!         +-WTOMP2           SET LBUFW, MEMORY PARTITIONING
!         +-WTOTB2           SET TABLES
!         +-WTOCL2-+-WTORD2  2ND STAGE MAIN CALCULATION
!                  +-WTOWT2
!----------------------------------------------------------------------
!  Algorithm is based on Yoshimine's backward-chaining technique
!     LCORE  = size of available main memory (byte), e.g. 20MB
!     LBOX   = size of a box = I/O unit size = 1 track, e.g. 40KB
!     NBOX   = number of boxes on main memory
!     NSOSO  = number of items to be sorted
!            = NSO*(NSO+1)/2
!     NSTEP  = number of step needed for sorting
!     MDIF   = number of different items in a box
!----------------------------------------------------------------------
!  Example
!     NSO = 400
!
!     NSOSO = 400*401/2 = 80200 --(semiorthogonalization)--> 30000
!
!     NBOX = LCORE / LBOX
!          = 20000 / 40 = 500
!
!     NSTEP = LOG(NSOSO) / LOG(NBOX)
!           = LOG(30000) / LOG(500) = 1.66
!
!     MDIF = (NSOSO + NBOX - 1) / NBOX
!          = (30000 + 500 - 1) / 500 = 60
!----------------------------------------------------------------------
!  Record structure of FT91 direct access work file
!     ITTR, NW, PQRS(LDOB91), I4PQRS(4,LDOB91)
!  Record structure of FT92 direct access work file
!     ITTR, NW, IR, IS, PQRS(LDOB92), I4PQ(2,LDOB92)
!----------------------------------------------------------------------
!  I/O
!     input:        NFT73(seq)
!     intermediate: NFT91(da)
!     output:       NFT92(da), NFT93(seq)
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
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
!
!      WRITE(*,*) 'Ordering SO integrals'
!      WRITE(*,*) 'WTODRV start'
!
!     ----- 1ST STAGE -----
!
      CALL WTOMP1(DD8BYT,AA)
      CALL WTOTB1(AA(KK(3))) ! ID1BOX
      CALL SETDRW(NFT91, LREC91, IERROR, IOSTAT)
      CALL WTOCL1(
     *  AA(KK(1)),AA(KK(2)),AA(KK(3)),AA(KK(4)),AA(KK(5)),
     *  AA(KK(6)),AA(KK(7)),AA(KK(8))
     *           )
!
!     ----- 2ND STAGE -----
!
      CALL WTOMP2(DD8BYT)
      CALL WTOTB2(AA(KK(3))) ! ID2BOX
      CALL SETDRW(NFT92, LREC92, IERROR, IOSTAT)
      CALL WTOCL2(
     *  AA(KK(1)),AA(KK(2)),AA(KK(3)),AA(KK(4)),AA(KK(5)),
     *  AA(KK(6)),AA(KK(7)),AA(KK(8)),AA(KK(9)),AA(KK(10))
     *           )
!
!      WRITE(*,*) 'WTODRV ended'
!      WRITE(*,*) 'SO integral ordering completed'
      RETURN
      END SUBROUTINE
