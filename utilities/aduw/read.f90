program readvals
  implicit none
       INTEGER :: I,IOS, KBUF, NBUF, MBUF
      REAL(8),DIMENSION(791) :: DBUFFW
      INTEGER,DIMENSION(4,791) :: IBUFFW


     OPEN(UNIT=64,FILE="moint.dat",STATUS='OLD',ACCESS='SEQUENTIAL',&
          FORM='UNFORMATTED')

      DO 
        READ(UNIT=64,ERR=901,IOSTAT=IOS) KBUF,NBUF,DBUFFW,IBUFFW

        MBUF = IABS(KBUF)
        IF( KBUF /= 0 ) THEN
          WRITE(*,*) (IBUFFW(:,I),DBUFFW(I),I=1,MBUF)
        END IF
        IF( KBUF <= 0 ) EXIT
      END DO


      CLOSE(UNIT=64)

  901 CONTINUE
      WRITE(*,*) IOS
      STOP 'WTHPRI 001'

end program readvals
