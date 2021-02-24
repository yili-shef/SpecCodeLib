SUBROUTINE ainvrnt (A,pa,ps,qa,qs,qw,ra,rs,v2)
  USE mconstant
  IMPLICIT NONE

  REAL(DP), DIMENSION(3,3), INTENT(IN)  :: A
  REAL(DP), INTENT(OUT) :: pa,ps,qa,qs,qw,ra,rs,v2

  INTEGER,  DIMENSION(3,3) :: I
  REAL(DP), DIMENSION(3,3) :: S
  REAL(DP), DIMENSION(3)   :: w,sw

  I=0
  I(1,1)=1; I(2,2)=1; I(3,3)=1
  
  S=.5*(A+TRANSPOSE(A))
  w(1)=A(2,3)-A(3,2)
  w(2)=A(3,1)-A(1,3)
  w(3)=A(1,2)-A(2,1)

  pa=-SUM(A,MASK=(I .EQ. 1))
  ps=-SUM(S,MASK=(I .EQ. 1))
  qa=-.5*SUM(MATMUL(A,A),MASK=(I .EQ. 1))
  qs=-.5*SUM(MATMUL(S,S),MASK=(I .EQ. 1))
  qw=SUM(w*w)/4.
  ra=-(1./3.)*SUM(MATMUL(A,MATMUL(A,A)),MASK=(I .EQ. 1))
  rs=-(1./3.)*SUM(MATMUL(S,MATMUL(S,S)),MASK=(I .EQ. 1))
  sw=MATMUL(S,w)
  v2=SUM(sw*sw)
  
END SUBROUTINE ainvrnt      
