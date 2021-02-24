complex(sp) function gasdev(t1,t2)
  use mconstant
  real(sp) :: t1, t2

  gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
  ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
end function gasdev
  
subroutine initvel(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
  use mconstant
  use msymmetrize
  implicit none

  integer,  intent(in) :: lx1,ly,lz,iseed
  real(sp), intent(in) :: akp, u0, kcut
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2, kx, ky, kz

  real(sp), dimension(lx1,ly,lz) :: tmp

  integer  :: idum,ii,jj,kk
  real(sp) :: ran1,tmp1,tmp2
  complex(sp) :: gasdev


  ! this is used to initialize ran1 with a given seed iseed.
  idum=iseed
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vx(ii,jj,kk)=gasdev(tmp1,tmp2)
    
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vy(ii,jj,kk)=gasdev(tmp1,tmp2)

    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vz(ii,jj,kk)=gasdev(tmp1,tmp2)

  end do
  end do
  end do

  tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
  vx = cmplx(real(vx) - kx*tmp, aimag(vx))
  vy = cmplx(real(vy) - ky*tmp, aimag(vy))
  vz = cmplx(real(vz) - kz*tmp, aimag(vz))

  tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
  vx = cmplx(real(vx), aimag(vx) - kx*tmp)
  vy = cmplx(real(vy), aimag(vy) - ky*tmp)
  vz = cmplx(real(vz), aimag(vz) - kz*tmp)

  call symmetrize(vx,lx1,ly,lz)
  call symmetrize(vy,lx1,ly,lz)
  call symmetrize(vz,lx1,ly,lz) 

  where(k2 .ge. kcut*kcut)
      vx = 0._sp
      vy = 0._sp
      vz = 0._sp
  endwhere

  ! different initial conditions can be imposed here

  tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 

  ! Another initial condition
  !where ( k2 .le. akp*akp)
  !    tmp = u0
  !else
  !    tmp = 0.
  !endwhere

  vx = vx * tmp
  vy = vy * tmp
  vz = vz * tmp

  tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
  tmp(1,:,:)=.5_sp*tmp(1,:,:)
  tmp1=sum(tmp)
  tmp2=(3._sp/2._sp)*u0*u0
  vx=sqrt(tmp2/tmp1)*vx
  vy=sqrt(tmp2/tmp1)*vy
  vz=sqrt(tmp2/tmp1)*vz
  ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

end subroutine initvel
