!
!----------------------------------------------------------------------
!  DESCRIPTION : This include file
!                is an include file for SYM4TR.
!                This file contains control data.
!  RELEASE :     v.01  gen = sya 2003-03-30 at chukyo-u
!----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: NFT74 = 74 ! moint.dat
      INTEGER,PARAMETER :: NFT91 = 91 ! bin-1
      INTEGER,PARAMETER :: NFT92 = 92 ! bin-2
      INTEGER,PARAMETER :: NFT93 = 93 ! control data for ordering
      INTEGER,PARAMETER :: NFT94 = 94 ! file merged from bin-2
      INTEGER,PARAMETER :: NFT95 = 95 ! half-tranformed integral
!      CHARACTER(LEN=12),PARAMETER :: FN74 = 'moint.dat'
! 
!----------------------------------------------------------------------
!
      INTEGER :: LDOB74 ! number of objects in a record of FT74
      INTEGER :: LDOB91 ! number of objects in a record of FT91
      INTEGER :: LDOB92 ! number of objects in a record of FT92
      INTEGER :: LDOB94 ! number of objects in a record of FT94
      INTEGER :: LDOB95 ! number of objects in a record of FT95
      INTEGER :: LREC91 ! record length of FT91 in Byte
      INTEGER :: LREC92 ! record length of FT92 in Byte
      INTEGER :: LREC95 ! record length of FT95 in Byte
      INTEGER :: NWTR91 ! number of records written on FT91
      INTEGER :: NWTR92 ! number of records written on FT92
      INTEGER :: NSOSO  ! NSO*(NSO+1)/2
      INTEGER :: M8CORE ! remaining memory in 8-ByteE
      INTEGER :: NBOX   ! total number of boxes on memory
      INTEGER :: NBXDIF ! number of different rs/pq indices in a box
!
!----------------------------------------------------------------------
!
      INTEGER :: NPQINT ! maximum nuber of SO-integrals (pq!rs)
!                         with a fixed rs-index, set by SUB.WTSWRS
      INTEGER :: NRSEXT ! number of actually existing rs-indices,
!                         set by SUB.WTSWRS
!
!----------------------------------------------------------------------
!
      INTEGER :: NMOMO  ! NMO*(NMO+1)/2
      INTEGER :: NTIMIJ
      INTEGER :: NTIMRS
      INTEGER :: KTIMIJ
      INTEGER :: KTIMRS
      INTEGER :: NDIFIJ
      INTEGER :: NDIFRS
      INTEGER :: MDIFIJ
      INTEGER :: MDIFRS
!
!----------------------------------------------------------------------
!
!     ...local common 
!
      COMMON/WTCOM/LDOB74,LDOB91,LDOB92,LDOB94,LDOB95,
     *             LREC91,LREC92,LREC95,NWTR91,NWTR92,
     *             M8CORE
      COMMON/WTOCOM/NSOSO, NBOX,  NBXDIF
      COMMON/WTSCOM/NPQINT,NRSEXT
      COMMON/WTHCOM/NMOMO, NTIMIJ,NTIMRS,NDIFIJ,NDIFRS,
     *              KTIMIJ,KTIMRS,MDIFIJ,MDIFRS
!
!----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: LDKMAX = 20
!
!                          ...padding for 8-byte boundarary
      INTEGER :: KKLAST
      INTEGER,DIMENSION(LDKMAX) :: KK,KK8
!
!     ...local common for SYM4TR
!
      COMMON/SYMCOM/KKLAST,KK,KK8
!
!----------------------------------------------------------------------
!
