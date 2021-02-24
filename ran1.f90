function ran1(idum)
USE mconstant
!  Random generator from Numerical Recipes, p.271
implicit none
integer idum,ia,im,iq,ir,ntab,ndiv
real(sp) :: ran1, am,eps,rnmx
parameter(ia = 16807, im = 2147483647, am = 1._SP/im, iq = 127773, ir = 2836,  &
          ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7_SP, rnmx = 1._SP - eps)
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

