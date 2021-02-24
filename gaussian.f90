SUBROUTINE gaussian(u,lx1,ly,lz)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz
  COMPLEX(SP), DIMENSION (lx1,ly,lz), INTENT(INOUT) :: u

  INTEGER  :: i,j,k,idum
  REAL(SP) :: t1,t2
  REAL     :: ran1

  u = 0.0_SP
  do i=1,lx1
   do j=1,ly
    do k=1,lz
      t1 = MAX(ran1(idum),smallest)
      t2 = MAX(ran1(idum),smallest)
      t2 = 2._SP*pi*t2
      u(i,j,k)=SQRT(-2.0_SP*LOG(t1))*CMPLX(COS(t2),SIN(t2))/SQRT(2.0_SP)
      ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
    END do
   END do
  END do

  RETURN
END SUBROUTINE gaussian

real function ran1(idum)
!  Random generator from Numerical Recipes, p.271
implicit none
integer idum,ia,im,iq,ir,ntab,ndiv
real am,eps,rnmx
parameter(ia = 16807, im = 2147483647, am = 1./im, iq = 127773, ir = 2836,  &
          ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7, rnmx = 1. - eps)
integer j,k
integer, dimension(ntab), save :: iv=(/(0,j=1,ntab)/)
integer, save :: iy = 0


if (idum <= 0 .or. iy == 0) then
  idum = max(-idum,1)
  do j = ntab + 8, 1, -1 
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k
    if (idum < 0) idum = idum + im
    if (j <= ntab) iv(j) = idum
  end do
  iy = iv(1)
end if

k = idum/iq
idum = ia*(idum - k*iq) - ir*k
if (idum < 0) idum = idum + im
j = 1 + iy/ndiv
iy = iv(j)
iv(j) = idum
ran1 = min(am*iy,rnmx)

End function Ran1

