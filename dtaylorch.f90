subroutine taylorch(A,ctau,rank)
  use mconstant
  implicit none

  integer, intent(in) :: rank
  real(dp), intent(in) :: ctau  ! equals approximately the norm of A.
  real(dp), dimension(rank,rank), intent(inout) :: A
  
  real(dp), dimension(rank,rank) :: A2
  real(dp) :: IIA, IIIA, c0, c1, c2
  integer :: j, ii, jj

  real(dp), parameter :: ca1=1., ca2=ca1/2., ca3=ca2/3., ca4=ca3/4.
  real(dp), parameter :: ca5=ca4/5., ca6=ca5/6., ca7=ca6/7.
  real(dp), parameter :: c00=1., c01=ca3/3., c02=ca5/6., c03=ca7/12.
  real(dp), parameter :: c10=1., c11=ca3/2., c12=ca4/3., c13=ca5/4., c14=ca6/3., c15=ca7
  real(dp), parameter :: c20=ca2, c21=ca4/2., c22=ca5/3., c23=ca6/4., c24=ca7/3.

  j=max(0, 1+floor(log(ctau)/log(2.)))

  A=A/2**j

  A2=matmul(A,A)
  
  IIA=0.
  do ii=1,rank
    IIA=A2(ii,ii)+IIA
  end do

  IIIA=0.
  do ii=1,rank
  do jj=1,rank
    IIIA=IIIA+A2(ii,jj)*A(jj,ii)
  end do
  end do

  
  ! Start here: Valid only for 3X3 matrice with zero trace.
  c0=c00+c01*IIIA+c02*IIA*IIIA+c03*IIA*IIA*IIIA
  c1=c10+c11*IIA+c12*IIIA+c13*IIA*IIA+c14*IIA*IIIA+c15*(IIIA*IIIA-IIA*IIA*IIA)
  c2=c20+c21*IIA+c22*IIIA+c23*IIA*IIA+c24*IIA*IIIA
  ! End here.
  
  A=c1*A+c2*A2
  do ii=1,rank
    A(ii,ii)=A(ii,ii)+c0
  end do

  do ii=1,j
    A=matmul(A,A)
  end do
end subroutine taylorch      
