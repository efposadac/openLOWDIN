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
      integer :: i
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
      multispecies = .false.
      if(scan(FN63, '.') /= scan(FN63, '.', .True.)) then
        multispecies = .True.
      end if

      if(.not.multispecies) then
!     $OMP private(i, j, nthreads, threadid, unitid, str, P, Q, R, S, VAL, NW, IBUFF, DBUFF) 
         nthreads = OMP_GET_NUM_THREADS()
         threadid =  OMP_GET_THREAD_NUM()
         unitid = 40 + threadid        
         

         write(str,*) threadid
         str = trim(adjustl(str))
C          print*, nthreads, threadid, unitid, trim(str)//FN63, 
C      >   integralStackSize
         
         OPEN(UNIT=unitid,FILE=trim(str)//FN63,STATUS='OLD',
     *        ACCESS='STREAM',
     *        FORM='UNFORMATTED')

         NW = 0
         IBUFF = 0
         DBUFF = 0

         rewind(unitid)

         loadintegrals : DO
            READ(UNIT=unitid, iostat=status) S(1:integralStackSize)
            if(status == -1 ) then
                print*, "end of file! file: ",trim(str)//FN63
                exit loadintegrals
            end if
            READ(UNIT=unitid) R(1:integralStackSize)
            READ(UNIT=unitid) Q(1:integralStackSize)
            READ(UNIT=unitid) P(1:integralStackSize)
            READ(UNIT=unitid) VAL(1:integralStackSize)
            
            readstack : DO j=1, integralStackSize
C                PRINT*, P(j),Q(j),R(j),S(j), VAL(j), j
               if(S(j) == -1) exit loadintegrals
               NW = NW + 1
               IBUFF(1,NW) = P(j)
               IBUFF(2,NW) = Q(j)
               IBUFF(3,NW) = R(j)
               IBUFF(4,NW) = S(j)
               DBUFF(NW) = VAL(j)
               

               IF( NW == DDOB73 ) THEN
!     $OMP CRITICAL
                  WRITE(NFT73) NW,NW,DBUFF,IBUFF
                  NW = 0
                  DBUFF=0
                  IBUFF=0
!     $OMP END CRITICAL
               END IF
            END DO readstack            
         END DO loadintegrals

!     $OMP CRITICAL
            WRITE(NFT73) -NW,NW,DBUFF,IBUFF
!     $OMP END CRITICAL
            CLOSE(UNIT=unitid)

!     $ END OMP
      else

         OPEN(UNIT=NFT63,FILE=FN63,STATUS='OLD',
     *        ACCESS='SEQUENTIAL',
     *        FORM='UNFORMATTED')
         

         NW = 0
         IBUFF = 0
         DBUFF = 0

         rewind(NFT63)

         loadintegrals_multi : DO
            READ(UNIT=NFT63, iostat=status) P(1:integralStackSize),
     *           Q(1:integralStackSize),
     *           R(1:integralStackSize),
     *           S(1:integralStackSize),
     *           VAL(1:integralStackSize)

            if(status == -1 ) then
                print*, "end of file! file: ",NFT63
                exit loadintegrals_multi
            end if

            DO j=1, integralStackSize
               if(p(j) == -1) goto 110
C                PRINT*, P(j),Q(j),R(j)+ otherNumberOfContractions,
C      >           S(j)+ otherNumberOfContractions, 
C      >           VAL(j) 
               NW = NW + 1
               IBUFF(1,NW) = P(j)
               IBUFF(2,NW) = Q(j)
               IBUFF(3,NW) = R(j)+otherNumberOfContractions
               IBUFF(4,NW) = S(j)+otherNumberOfContractions
               DBUFF(NW) = VAL(j)
               
               IF( NW == DDOB73 ) THEN
                  WRITE(NFT73) NW,NW,DBUFF,IBUFF
                  NW = 0
                  DBUFF=0
                  IBUFF=0
               END IF
            END DO
         END DO loadintegrals_multi
 110     CONTINUE
         
!     IF(i .LT. nproc) THEN
!     WRITE(NFT73) NW,NW,DBUFF,IBUFF
!     END IF 
         
         CLOSE(UNIT=NFT63)
         
      end if
      
      WRITE(NFT73) -NW,NW,DBUFF,IBUFF      
      CLOSE(UNIT=NFT73)
      
c$$$  WRITE(*,*) 'WMFINT ended'
      
      RETURN
      END SUBROUTINE
      
      
