complex(sp) function gasdev(t1,t2)
  use mconstant
  real(sp) :: t1, t2

  gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
  ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
end function gasdev
  

! This subroutine is highly ad hoc and often modified for the problem at hand. 
subroutine initphi(phi,kx,ky,kz,k2,iseed,lx1,ly,lz,kcut)
  use mconstant
  use msymmetrize
  implicit none

  integer,  intent(in) :: lx1,ly,lz,iseed
  real(sp), intent(in) :: kcut
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: phi
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2, kx, ky, kz

  integer  :: idum,ii,jj,kk,lx2
  real(sp) :: ran1,tmp1,tmp2
  complex(sp) :: gasdev

  lx2 = (lx1-1)*(lx1-1)
  idum=iseed

  phi = 0._sp
  do kk=1,3
  do jj=1,3
  do ii=1,3
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    phi(ii,jj,kk)=gasdev(tmp1,tmp2)
  end do
  end do
  end do
 
  call symmetrize(phi,lx1,ly,lz)

  where(k2 .ge. kcut*kcut)
      phi = 0._sp
  endwhere

  tmp1 = real(sum(phi*conjg(phi)),sp)
  phi = sqrt(1.6_sp / tmp1) * phi

end subroutine initphi
