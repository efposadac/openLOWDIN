!
!----------------------------------------------------------------------
!  DESCRIPTION : This include file
!                sets constants for PREPAR and SYM4TR programs.
!                This file contains control data.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!----------------------------------------------------------------------
!
!     ...constant(s)
!
      INTEGER,PARAMETER :: LDSYM  = 8 ! number of irreps
      INTEGER,PARAMETER :: LDABEL = 8 ! number of Abel point groups
!
!     ...constant(s) but environment dependent
!
      INTEGER,PARAMETER :: LDTRAK = 19000 ! disk track size (Byte)
      INTEGER :: LDOB72 ! number of objects in a record of FT72
                        ! MO coefficient matrix
      INTEGER :: LDOB73 ! number of objects in a record of FT73
                        ! SO integrals
!
      INTEGER,PARAMETER :: NFT61 = 61
      INTEGER,PARAMETER :: NFT62 = 62
      INTEGER,PARAMETER :: NFT63 = 63
      INTEGER,PARAMETER :: NFT71 = 71
      INTEGER,PARAMETER :: NFT72 = 72
      INTEGER,PARAMETER :: NFT73 = 73
!      CHARACTER(LEN=12),PARAMETER :: FN61 = 'inf.txt'
!      CHARACTER(LEN=12),PARAMETER :: FN62 = 'mo.values'
!      CHARACTER(LEN=12),PARAMETER :: FN63 = '.ints'
!      CHARACTER(LEN=12),PARAMETER :: FN71 = 'inf.dat'
!      CHARACTER(LEN=12),PARAMETER :: FN72 = 'mocoef.dat'
!      CHARACTER(LEN=12),PARAMETER :: FN73 = 'soint.dat'    

      CHARACTER(LEN=60) :: FN61 
      CHARACTER(LEN=60) :: FN62 
      CHARACTER(LEN=60) :: FN63 
      CHARACTER(LEN=60) :: FN71 
      CHARACTER(LEN=60) :: FN72 
      CHARACTER(LEN=60) :: FN73 
      CHARACTER(LEN=60) :: FN74

!
!----------------------------------------------------------------------
!
!     ...working variable(s)
!
      CHARACTER(LEN=72) :: TITLEC ! title comment
      CHARACTER(LEN=4) :: NAMEPT ! belonging Abel point group name
      INTEGER :: NABEL ! belonging Abel point group ID
      INTEGER :: NSYM ! number of irreps
      INTEGER :: NSO ! number of SOs
      INTEGER :: NMO ! number of active MOs
      INTEGER :: NPFLAG ! print flag (0=normal, 1=dump)
      INTEGER :: NRFLAG  ! run flag (0=delete, 1=keep work files)
      REAL(KIND=LDREAL) :: THRSOI ! threshold for SO integrals.
      REAL(KIND=LDREAL) :: THRHLF ! threshold for half-tr int.
      REAL(KIND=LDREAL) :: THRMOI ! threshold for transformed int.
!
!     ...working array(s)
!
      INTEGER,DIMENSION(LDSYM) :: NSOG ! number of SO for each irrep
      INTEGER,DIMENSION(LDSYM) :: NMOG ! number of MO for each irrep
!
!----------------------------------------------------------------------
!
!     ...allocate local COMMON
!
      COMMON/PREPAI/LDOB72,LDOB73,NSYM,NSOG,NMOG,NSO,NMO,NPFLAG,NRFLAG
      COMMON/PREPAD/THRSOI,THRHLF,THRMOI
      COMMON/PREPAC/TITLEC,NAMEPT
      COMMON/FILES/FN61,FN62,FN63,FN71,FN72,FN73,FN74
!
!----------------------------------------------------------------------
!
