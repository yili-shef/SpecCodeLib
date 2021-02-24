subroutine taylorch(A,ctau,rank)
  use mconstant
  implicit none

  integer, intent(in) :: rank
  real(sp), intent(in) :: ctau  ! equals approximately the norm of A.
  real(sp), dimension(rank,rank), intent(inout) :: A
  
  real(sp), dimension(rank,rank) :: A2
  real(sp) :: I2A, I3A, c0, c1, c2
  integer :: j, ii

  real(sp), parameter :: ca1=1., ca2=ca1/2., ca3=ca2/3., ca4=ca3/4.
  real(sp), parameter :: ca5=ca4/5., ca6=ca5/6., ca7=ca6/7.
  real(sp), parameter :: c00=1., c01=ca3/3., c02=ca5/6., c03=ca6/9., c04=ca7/12.
  real(sp), parameter :: c10=1., c11=ca3/2., c12=ca4/3., c13=ca5/4., c14=ca6/3., c15=ca7
  real(sp), parameter :: c20=ca2, c21=ca4/2., c22=ca5/3., c23=ca6/4., c24=ca7/3.

  j=max(0, 1+floor(log(ctau)/log(2.)))

  A=A/2**j

  A2=matmul(A,A)
  
  I2A=0.
  do ii=1,rank
    I2A=A2(ii,ii)+I2A
  end do
  I3A = sum(A2*transpose(A))
  
  ! Start here: Valid only for 3X3 matrice with zero trace.
  c0=c00+c01*I3A+c02*I2A*I3A+c03*I3A*I3A+c04*I2A*I2A*I3A
  c1=c10+c11*I2A+c12*I3A+c13*I2A*I2A+c14*I2A*I3A+c15*(I3A*I3A-I2A*I2A*I2A)
  c2=c20+c21*I2A+c22*I3A+c23*I2A*I2A+c24*I2A*I3A
  ! End here.
  
  A=c1*A+c2*A2
  do ii=1,rank
    A(ii,ii)=A(ii,ii)+c0
  end do

  do ii=1,j
    A=matmul(A,A)
  end do
end subroutine taylorch      
