      SUBROUTINE WMFINT(
     *     DDOB73, 
     *     IBUFF,  DBUFF, nproc, integralStackSize,
     *     otherNumberOfContractions)
!----------------------------------------------------------------------
!     DESCRIPTION : This routine
!     reads SO basis integrals from text files and write 
!     them to a binary file.
!     RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!     CALLED BY :   WMFILE
!----------------------------------------------------------------------
!     record structure of binary SO integral file
!     KW, NW, DBUFF(1:LDOB73), IBUFF(1:4,1:LDOB73)
!     
!     meaning of variable
!     NW : number of SO integrals in a record
!     KW : = NW, 
!     = 0 or -NW for the last record
!     IBUFF(1,i) : p
!     IBUFF(2,i) : q
!     IBUFF(2,i) : r
!     IBUFF(2,i) : s
!     DBUFF(i) : integral value
!----------------------------------------------------------------------
!     N.B.
!     DDOB73 == 791
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
      INCLUDE 'omp_lib.h'
!     
!     ...declare variable(s) in argument list
!     
      INTEGER,INTENT(IN) :: DDOB73
      REAL(KIND=LDREAL),DIMENSION(DDOB73),INTENT(INOUT) :: DBUFF
      INTEGER,DIMENSION(4,DDOB73),INTENT(INOUT) :: IBUFF
      integer :: nproc
      integer :: otherNumberOfContractions
      integer :: integralStackSize
      integer :: j
      integer :: nthreads
      integer :: threadid
      integer :: unitid
      character(50) :: str
      logical :: multispecies
!     

!     ...define variable(s) in code
!     
      
      REAL(KIND=LDREAL),DIMENSION(integralStackSize):: VAL
      INTEGER,DIMENSION(integralStackSize) :: P,Q,R,S
      INTEGER :: NW, status
!----------------------------------------------------------------------
!     
!     WRITE(*,*) 'WMFINT start'
!     
!     ...read SO basis integrals
!     
      OPEN(UNIT=NFT73,FILE=FN73,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     *     FORM='UNFORMATTED')

      str = ""

c$OMP PARALLEL 
c$omp& private(j, nthreads, threadid, unitid, str, P, Q, R, S, VAL)
c$omp& private(NW, IBUFF, DBUFF) 
      nthreads = OMP_GET_NUM_THREADS()
      threadid =  OMP_GET_THREAD_NUM()
      unitid = 40 + threadid        
      

      write(str,*) threadid
      str = trim(adjustl(str))

C       print*, nthreads, threadid, unitid, trim(str)//FN63, 
C      >     integralStackSize
      
      OPEN(UNIT=unitid,FILE=trim(str)//FN63,STATUS='OLD',
     *     ACCESS='STREAM',
     *     FORM='UNFORMATTED')

      NW = 0
      IBUFF = 0
      DBUFF = 0

      rewind(unitid)

      loadintegrals_multi : DO
      READ(UNIT=unitid, iostat=status) P(1:integralStackSize),
     *     Q(1:integralStackSize),
     *     R(1:integralStackSize),
     *     S(1:integralStackSize),
     *     VAL(1:integralStackSize)

      if(status == -1 ) then
         print*, "end of file! file: ",trim(str)//FN63
         exit loadintegrals_multi
      end if

      
      DO j=1, integralStackSize
C     PRINT*, P(j),Q(j),R(j)+ otherNumberOfContractions,
C     >           S(j)+ otherNumberOfContractions, 
C     >           VAL(j) 
         if(P(j) == -1) exit loadintegrals_multi
         NW = NW + 1
         IBUFF(1,NW) = P(j)
         IBUFF(2,NW) = Q(j)
         IBUFF(3,NW) = R(j)+otherNumberOfContractions
         IBUFF(4,NW) = S(j)+otherNumberOfContractions
         DBUFF(NW) = VAL(j)
         

         IF( NW == DDOB73 ) THEN
c$OMP CRITICAL
            WRITE(NFT73) NW,NW,DBUFF,IBUFF
            NW = 0
            DBUFF=0
            IBUFF=0
c$OMP END CRITICAL
         END IF
      END DO             
      END DO loadintegrals_multi

c$OMP CRITICAL
      WRITE(NFT73) NW,NW,DBUFF,IBUFF
c$OMP END CRITICAL

c$OMP BARRIER
      if(threadid == 0) then
         NW = 0
         WRITE(NFT73) NW,NW,DBUFF,IBUFF
      end if
      
      CLOSE(UNIT=unitid)

c$OMP END PARALLEL
      

      WRITE(NFT73) -NW,NW,DBUFF,IBUFF      
      CLOSE(UNIT=NFT73)
      
c$$$  WRITE(*,*) 'WMFINT ended'
      
      RETURN
      END SUBROUTINE
      
      
      
