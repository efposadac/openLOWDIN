      SUBROUTINE WMINPT(numberOfAO, atomicIntegralThreshold,
     * halfTransformedThreshold, molecularIntegralsThreshold)
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads line input data from file.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   PREPAR
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!             ...define variable(s) in code
!
      integer :: numberOfAO
      real(8) :: atomicIntegralThreshold
      real(8) :: halfTransformedThreshold
      real(8) :: molecularIntegralsThreshold
      INTEGER :: J1, ISYM
!
!----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!! INICIO DE RUTINA ORIGINAL  !!!!!!!!!!!!!!!!!!!!!!!!!!
C$$$!      WRITE(*,*) 'WMINPT start'
C$$$!
C$$$      OPEN(UNIT=NFT61, FILE=FN61, STATUS='OLD',
C$$$     *     ACCESS='SEQUENTIAL', FORM='FORMATTED')
C$$$!
C$$$!     ...title comment
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) TITLEC
C$$$!      WRITE(*,"(' Title of this job: ',/,1X,A72)") TITLEC
C$$$!
C$$$!     ...name of point group
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) NAMEPT
C$$$!      WRITE(*,"(' Name of point group: ',A4)") NAMEPT
C$$$      CALL WMABEL
C$$$!
C$$$!     ...NSOG
C$$$!
C$$$      NSOG = 0
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) (NSOG(ISYM),ISYM=1,NSYM)
C$$$!      WRITE(*,"(' NSOG: ',8I4)") (NSOG(ISYM),ISYM=1,NSYM)
C$$$      NSO = 0
C$$$      DO ISYM = 1, NSYM
C$$$        NSO = NSO + NSOG(ISYM)
C$$$      END DO
C$$$!      WRITE(*,"(' Total number of SOs: ',I4)") NSO
C$$$!
C$$$!     ...NMOG
C$$$!
C$$$      NMOG = 0
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) (NMOG(ISYM),ISYM=1,NSYM)
C$$$!      WRITE(*,"(' NMOG: ',8I4)") (NMOG(ISYM),ISYM=1,NSYM)
C$$$      NMO = 0
C$$$      DO ISYM = 1, NSYM
C$$$        NMO = NMO + NMOG(ISYM)
C$$$      END DO
C$$$!      WRITE(*,"(' Total number of active MOs: ',I4)") NMO
C$$$!
C$$$!     ...NPFLAG
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) NPFLAG
C$$$      IF( NPFLAG /= 0 ) NPFLAG = 1
C$$$!      WRITE(*,"(' Print flag: ',I1)") NPFLAG
C$$$!
C$$$!     ...NRFLAG
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) NRFLAG
C$$$      IF( NRFLAG /= 0 ) NRFLAG = 1
C$$$!      WRITE(*,"(' Run flag: ',I1)") NRFLAG
C$$$!
C$$$!     ...THRSOI
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) THRSOI
C$$$!      WRITE(*,"(' THRSOI: ',D20.10)") THRSOI
C$$$!
C$$$!     ...THRHLF
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) THRHLF
C$$$!      WRITE(*,"(' THRHLF: ',D20.10)") THRHLF
C$$$!
C$$$!     ...THRMOI
C$$$!
C$$$      READ(NFT61,"(1X)")
C$$$      READ(NFT61,*) THRMOI
C$$$!      WRITE(*,"(' THRMOI: ',D20.10)") THRMOI
C$$$!
C$$$      CLOSE(UNIT=NFT61)
C$$$!      WRITE(*,*) 'WMINPT ended'
C$$$!
C$$$      RETURN
!!!!!!!!!!!!!!!!!! FIN DE RUTINA ORIGINAL  !!!!!!!!!!!!!!!!!!!!!!!!!!

      TITLEC=" "
      NAMEPT="C1" !! Pueden definise el grupo puntual asociado a la simetria
      CALL WMABEL
      NSOG = 0
      NSOG(1) = numberOfAO 
      NSO = 0
      DO ISYM = 1, NSYM
        NSO = NSO + NSOG(ISYM)
      END DO

      NMOG = 0
      NMOG(1) = numberOfAO 
      NMO = 0
      DO ISYM = 1, NSYM
        NMO = NMO + NMOG(ISYM)
      END DO
      NPFLAG=1
      NRFLAG=0
      THRSOI=atomicIntegralThreshold
      THRHLF=halfTransformedThreshold
      THRMOI=molecularIntegralsThreshold

      RETURN

      END SUBROUTINE
