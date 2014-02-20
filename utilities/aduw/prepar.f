      PROGRAM PREPAR
!----------------------------------------------------------------------
!  DESCRIPTION : This routine 
!                is main program for preparing files for integral
!                transformation.  It makes binary files of MO 
!                coefficient matrix and SO basis integrals from their
!                text format files.
!  AUTHORS :     Shigeyoshi YAMAMOTO
!  RELEASE :     v.01  gen = S.Yamammoto(sya) 2003-03-29 at chukyo-u
!----------------------------------------------------------------------
!  Files
!     FT61 : inf.txt / old / sequential / formatted 
!     FT62 : mocoef.txt / old / sequential / formatted / text
!     FT63 : soint.txt / old / sequential / formatted / text
!     FT71 : inf.dat / unknown / sequential / unformatted / binary
!     FT72 : mocoef.dat / unknown / sequential / unformatted / binary
!     FT73 : soint.dat / unknown / sequential / unformatted / binary
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!     ...memory size
!
      INTEGER,PARAMETER :: LDMB = 800 ! MB
      INTEGER,PARAMETER :: LDBYTE = LDMB * 10**6 ! Byte
      INTEGER,PARAMETER :: LD8BYT = LDBYTE/LDREAL ! 8-Byte
!
!     ...allocate big array by COMMON
!
      REAL(KIND=LDREAL),DIMENSION(LD8BYT) :: BIGARY
      COMMON/BIGARY/BIGARY
!
!----------------------------------------------------------------------
!
      WRITE(*,*) 'PREPAR start'
!
!     ...initializes control data
!
      CALL WMINIT
!
!     ...read input data from text file
!
      CALL WMINPT
!
!     ...check consistency of input data and save them
!
      CALL WMINFO

!     ...set integral file and mo coefficient file
!
      CALL WMFILE(LD8BYT,BIGARY)
!
!     ...check block-diagonality of MO matrix
!
      CALL WMCHKD(LD8BYT,BIGARY)
!
      WRITE(*,*) 'PREPAR ended'
!
      END PROGRAM
