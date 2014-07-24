      SUBROUTINE SETDRW(NUNIT,N1BYTE,NERROR,IOSTAT)
!----------------------------------------------------------------------
!   OPEN TEMPORARY DIRECT ACCESS FILES
!   VHUOPN-SETDRW
!   GENERATION = S.Yamamoto 1994-08-22
!----------------------------------------------------------------------
!   meaning of variables
!     NUNIT  : logical unit number of the file
!     N1BYTE : record length counted by Byte
!     NERROR : error status, 0 for OK
!     IOSTAT : error code
!     MACHIN : kind of machine, 
!              if 1, recl is specified by Byte.
!                 on HITAC,FACOM,NEC-SX,SUN,HP,IBM
!              if 2, recl is specified by word.
!                 on  DEC-VAX,TITAN,IRIS-4D
!----------------------------------------------------------------------
      IMPLICIT INTEGER(A-Z)
      DATA MACHIN/1/
!----------------------------------------------------------------------
 9000 FORMAT(/1x,'**** ERROR STOP SUB.SETDRW ****',
     *       /1x,'     error on opening work direct-access file',
     *       /1x,'     nunit, iostat, lrecl = ',3i12)
!----------------------------------------------------------------------
      NERROR = 0
      IOSTAT = 0
!
      IF( MACHIN .EQ. 1 ) THEN
          LRECL = N1BYTE 
      ELSE IF( MACHIN .EQ. 2 ) THEN
          LRECL = N1BYTE / 4
      ELSE
          STOP 'SETDRW'
      ENDIF
!
      OPEN(UNIT=NUNIT,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=LRECL,
     *     STATUS='UNKNOWN',IOSTAT=IOSTAT,ERR=95)
      RETURN 
!
   95 CONTINUE
      WRITE(*,9000) NUNIT,IOSTAT,LRECL
      STOP 'SETDRW 001'
      END
