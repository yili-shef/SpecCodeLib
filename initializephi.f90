subroutine initializephi(phi,k2,iseed,lx1,ly,lz,akp,phirms)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in) :: lx1,ly,lz,iseed
  real(sp), intent(in) :: akp, phirms
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: phi
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2

  real(sp), dimension(lx1,ly,lz) :: tmp

  integer  :: idum,ii,jj,kk
  real(sp) :: ran1, tmp1, tmp2, ignore_me, const

  
  const = 1./(ly*lz*(lx1-1)*2)

  ! this is used to initialize ran1 with a given seed iseed.
  idum=iseed
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    phi(ii,jj,kk)=gasdev(tmp1,tmp2)
  end do
  end do
  end do

  call symmetrize(phi,k2,lx1,ly,lz)

  ! Following Eswaran and Pope 1988
  where (k2 .ge. akp*akp)
      phi = phirms * phi / sqrt(4.*pi*k2)
  elsewhere
      phi = 0.
  endwhere

  call rfftwnd_f77_one_complex_to_real(c2r3d,phi,ignore_me)
  phi = phi * const

  where ( real(phi) .ge. 0 ) 
      phi = cmplx( 1. , aimag(phi) )
  elsewhere
      phi = cmplx( -1. , aimag(phi) )
  endwhere

  where ( aimag(phi) .ge. 0 ) 
      phi = cmplx( real(phi) , 1. )
  elsewhere
      phi = cmplx( real(phi) , -1. )
  endwhere

  call rfftwnd_f77_one_real_to_complex(r2c3d,phi,ignore_me)

  where ( k2 .ge. (2*akp)**2 ) 
      phi = 0.
  endwhere

  tmp = phi*conjg(phi)
  tmp(1,:,:)=.5_sp*tmp(1,:,:)

  tmp1=sum(tmp)
  tmp2= phirms * phirms  / 2._sp
  phi=sqrt(tmp2/tmp1)*phi
  ! TODO: Rescaling to match rms phi. Mean phi is assumed to be ZERO.

contains

  complex(sp) function gasdev(t1,t2)
    real(sp) :: t1, t2

    gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
    ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
    ! each component is a Gaussian random number

  end function gasdev

end subroutine initializephi

