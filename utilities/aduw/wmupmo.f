      SUBROUTINE WMUPMO(
     *                  NSO, NMO, CSOMO
     *                 )
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                prints MO coefficient matrix.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WMFCOE
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
!
!             ...declare variable(s) in argument list
!
      INTEGER,INTENT(IN) :: NSO, NMO
      REAL(KIND=LDREAL),DIMENSION(NSO,NMO),INTENT(IN) :: CSOMO
!
!             ...define constant(s) in code
!
      INTEGER,PARAMETER :: NCOL = 10
      REAL(KIND=LDREAL) :: THR = 1.0D-10
!
!             ...define variable(s) in code
!
      INTEGER :: NTIME, ITIME, ISO, IMO, IM1, IM2, ICOL, JCOL
      CHARACTER(LEN=13),DIMENSION(NCOL) :: VAL
!
!----------------------------------------------------------------------
 1000 FORMAT(4X,I7,10I11)
 2000 FORMAT(I4,10(1X,A10))
!----------------------------------------------------------------------
!
      WRITE(*,*) 'WMUPMO start'
!
!      WRITE(*,*) 'MO coefficient matrix (SOMO)'
      NTIME = (NMO + NCOL - 1) / NCOL
      DO ITIME = 1, NTIME
        IM1 = (ITIME - 1) * NCOL + 1
        IM2 = IM1 + NCOL - 1
        IM2 = MIN(NMO,IM2)
        WRITE(*,1000) (IMO,IMO=IM1,IM2)
        DO ISO = 1, NSO
          VAL = ''
          ICOL = 0
          DO IMO = IM1, IM2
            ICOL = ICOL + 1
            IF( ABS(CSOMO(ISO,IMO)) > THR ) 
     *        WRITE(VAL(ICOL),"(D10.3)") CSOMO(ISO,IMO)
          END DO
          WRITE(*,2000) ISO,(VAL(JCOL),JCOL=1,ICOL)
        END DO
      END DO
!
      WRITE(*,*) 'WMUPMO ended'
!
      RETURN
      END SUBROUTINE
