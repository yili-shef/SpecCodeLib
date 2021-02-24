RECURSIVE SUBROUTINE dtmnnts(A,N,det)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: N
  REAL(DP), INTENT(OUT) :: det
  REAL(DP), DIMENSION(N,N), INTENT(IN) :: A

  REAL(DP), DIMENSION(N-1,N-1) :: Asub

  REAL(DP) :: detsub
  INTEGER  :: i,nt

  det=0.
  IF (N .EQ. 1) THEN
          det=A(1,1)
          RETURN
  ELSE
          nt=1
          DO i=1,N
          CALL subm(A,1,i,Asub,N)
          CALL dtmnnts(Asub,N-1,detsub)
          det=det+detsub*A(1,i)*nt
          nt=nt*(-1)
          END DO
  END IF
END SUBROUTINE dtmnnts

SUBROUTINE subm(A,i,j,Asub,N)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: i,j,N
  REAL(DP), DIMENSION(N,N), INTENT(IN) :: A
  REAL(DP), DIMENSION(N-1,N-1),INTENT(OUT) :: Asub

  IF (i .EQ. 1 .AND. j .EQ. 1) THEN
    Asub=A(2:N,2:N)
  ELSE IF (i .EQ. 1) THEN
    Asub(:,1:j-1)=A(2:N,1:j-1)
    Asub(:,j:N-1)=A(2:N,j+1:N)
  ELSE IF (j .EQ. 1) THEN
    Asub(1:i-1,:)=A(1:i-1,2:N)
    Asub(i:N-1,:)=A(i+1:N,2:N)
  ELSE
    Asub(1:i-1,1:j-1)=A(1:i-1,1:j-1)
    Asub(i:N-1,1:j-1)=A(i+1:N,1:j-1)
    Asub(1:i-1,j:N-1)=A(1:i-1,j+1:N)
    Asub(i:N-1,j:N-1)=A(i+1:N,j+1:N)
  END IF
END SUBROUTINE subm
