      SUBROUTINE WMFINT(
     *                  DDOB73, 
     *                  IBUFF,  DBUFF, nproc, integralStackSize,
     *     otherNumberOfContractions)
!----------------------------------------------------------------------
!  DESCRIPTION : This routine
!                reads SO basis integrals from text files and write 
!                them to a binary file.
!  RELEASE :     v.01  gen = sya 2003-03-29 at chukyo-u
!  CALLED BY :   WMFILE
!----------------------------------------------------------------------
!  record structure of binary SO integral file
!     KW, NW, DBUFF(1:LDOB73), IBUFF(1:4,1:LDOB73)
! 
!  meaning of variable
!     NW : number of SO integrals in a record
!     KW : = NW, 
!          = 0 or -NW for the last record
!     IBUFF(1,i) : p
!     IBUFF(2,i) : q
!     IBUFF(2,i) : r
!     IBUFF(2,i) : s
!     DBUFF(i) : integral value
!----------------------------------------------------------------------
!  N.B.
!     DDOB73 == 791
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'declar.h'
      INCLUDE 'prepar.h'
!
!             ...declare variable(s) in argument list
!
      INTEGER,INTENT(IN) :: DDOB73
      REAL(KIND=LDREAL),DIMENSION(DDOB73),INTENT(INOUT) :: DBUFF
      INTEGER,DIMENSION(4,DDOB73),INTENT(INOUT) :: IBUFF
      integer :: nproc
      integer :: otherNumberOfContractions
      integer :: integralStackSize
      integer :: i
      integer :: j
      character(50) :: str
      logical :: multispecies
!     

!             ...define variable(s) in code
!     
      
      REAL(KIND=LDREAL),DIMENSION(integralStackSize):: VAL
      INTEGER*2,DIMENSION(integralStackSize) :: P,Q,R,S
      INTEGER :: NW
!----------------------------------------------------------------------
!
!      WRITE(*,*) 'WMFINT start'
!
!     ...read SO basis integrals
!
      OPEN(UNIT=NFT73,FILE=FN73,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     *     FORM='UNFORMATTED')

      str = ""
      multispecies = .false.
      if(nproc==0) then 
         multispecies = .true.
      end if

      if(.not.multispecies) then

         DO i=1, nproc         
 
            write(str,*) i
            str = trim(adjustl(str))
            OPEN(UNIT=NFT63,FILE=trim(str)//FN63,STATUS='OLD',
     *           ACCESS='SEQUENTIAL',
     *           FORM='UNFORMATTED')

            NW = 0
            IBUFF = 0
            DBUFF = 0

            rewind(NFT63)

            DO
               READ(UNIT=NFT63) P(1:integralStackSize),
     *              Q(1:integralStackSize),
     *              R(1:integralStackSize),S(1:integralStackSize),
     *              VAL(1:integralStackSize)
               
               DO j=1, integralStackSize
                  if(p(j) == -1) goto 109
                  !PRINT*, P(j),Q(j),R(j),S(j), VAL(j)
                  NW = NW + 1
                  IBUFF(1,NW) = P(j)
                  IBUFF(2,NW) = Q(j)
                  IBUFF(3,NW) = R(j)
                  IBUFF(4,NW) = S(j)
                  DBUFF(NW) = VAL(j)
                  
                  IF( NW == DDOB73 ) THEN
                     WRITE(NFT73) NW,NW,DBUFF,IBUFF
                     NW = 0
                     DBUFF=0
                     IBUFF=0
                  END IF
               END DO
            END DO
 109        CONTINUE
            
            IF(i .LT. nproc) THEN
               WRITE(NFT73) NW,NW,DBUFF,IBUFF
            END IF      
            CLOSE(UNIT=NFT63)
            
         END DO
         
      else

         OPEN(UNIT=NFT63,FILE=FN63,STATUS='OLD',
     *        ACCESS='SEQUENTIAL',
     *        FORM='UNFORMATTED')
         

         NW = 0
         IBUFF = 0
         DBUFF = 0

         rewind(NFT63)

         DO
            READ(UNIT=NFT63) P(1:integralStackSize),
     *           Q(1:integralStackSize),
     *           R(1:integralStackSize),S(1:integralStackSize),
     *           VAL(1:integralStackSize)
            
            DO j=1, integralStackSize
               if(p(j) == -1) goto 110
c$$$               PRINT*, P(j),Q(j),R(j)+ otherNumberOfContractions,
c$$$     *              S(j)+ otherNumberOfContractions, 
c$$$     *              VAL(j) 
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
         END DO
 110     CONTINUE
  
!         IF(i .LT. nproc) THEN
!            WRITE(NFT73) NW,NW,DBUFF,IBUFF
!         END IF 
     
         CLOSE(UNIT=NFT63)
  
      end if
      
      WRITE(NFT73) -NW,NW,DBUFF,IBUFF      
      CLOSE(UNIT=NFT73)
     
c$$$      WRITE(*,*) 'WMFINT ended'
     
      RETURN
      END SUBROUTINE
      
