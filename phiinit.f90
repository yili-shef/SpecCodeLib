subroutine phiinit (phi, k2, iseed, lx1, ly, lz, akp, phirms)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1, ly, lz, iseed
  real(sp), intent(in) :: akp, phirms

  complex(sp), dimension(lx1,ly,lz), intent(out) :: phi
  real(sp), dimension(lx1,ly,lz), intent(in) :: k2

  integer  :: idum, ii, jj, kk
  real(sp) :: ran1, tmp1, tmp2, tmp3, k

  idum=iseed
  tmp3 = 0._sp
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    k = sqrt( k2(ii,jj,kk) )

    k = phirms * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * k * exp (-k*k/(akp*akp)) 

    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    phi(ii,jj,kk)=gasdev(tmp1,tmp2) * k

    tmp1 = phi(ii,jj,kk) * conjg(phi(ii,jj,kk))
    if (ii .eq. 1) tmp1 = tmp1 * .5_sp
    tmp3 = tmp3 + tmp1

  end do
  end do
  end do

  call hermitianize(phi,lx1,ly,lz)
  call setzero(phi,k2,lx1,ly,lz)

  tmp2=(1._sp/2._sp)*phirms*phirms
  phi=sqrt(tmp2/tmp3)*phi
  ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

contains

  complex(sp) function gasdev(t1,t2)
    real(sp) :: t1, t2

    gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
    ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
  end function gasdev

end subroutine phiinit 
