      SUBROUTINE VUTCN2(N,MTR)
C---------------------------------------------------------------------
C  SET TRIANGULAR ADDRESS TABLE MTR
C  ## VECTORIZATION ##
C  CREATION = S.YAMAMOTO 1984-07-04
C---------------------------------------------------------------------
      INTEGER N
      INTEGER MTR(N)
      INTEGER J1
C---------------------------------------------------------------------
      DO 100 J1 = 1, N
        MTR(J1) = J1*(J1-1)/2
  100 CONTINUE
      RETURN
      END
