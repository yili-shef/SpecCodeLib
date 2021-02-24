subroutine initdealiased(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
  use mconstant
  use msymmetrize
  implicit none

  integer,  intent(in) :: lx1,ly,lz,iseed
  real(sp), intent(in) :: akp, u0, kcut
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2,kx,ky,kz

  integer  :: idum,ii,jj,kk
  real(sp) :: ran1,tmp1,tmp2, kcut2, kxii, kyjj, kzkk
  complex(sp) :: tmpz, vxt, vyt, vzt, gasdev


  kcut2 = kcut * kcut

  ! this is used to initialize ran1 with a given seed iseed.
  idum=iseed
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1

    kzkk=kz(ii,jj,kk)
    kyjj=ky(ii,jj,kk)
    kxii=kx(ii,jj,kk)

    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vxt = gasdev(tmp1,tmp2)
    
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vyt = gasdev(tmp1,tmp2)

    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    vzt = gasdev(tmp1,tmp2)

    tmpz=kxii*vxt+kyjj*vyt+kzkk*vzt
    tmpz=tmpz/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
    vxt=vxt-kxii*tmpz
    vyt=vyt-kyjj*tmpz
    vzt=vzt-kzkk*tmpz

    tmp2 = k2(ii,jj,kk)
    tmp1 = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(tmp2) * exp (-tmp2/(akp*akp)) 

    vx(ii,jj,kk) = vxt * tmp1
    vy(ii,jj,kk) = vyt * tmp1
    vz(ii,jj,kk) = vzt * tmp1
  end do
  end do
  end do

  call symmetrize(vx,lx1,ly,lz)
  call symmetrize(vy,lx1,ly,lz)
  call symmetrize(vz,lx1,ly,lz) 

  where(k2 .ge. kcut2)
      vx = 0.
      vy = 0.
      vz = 0.
  endwhere

  ! Another initial condition
  !where ( k2 .le. akp*akp)
  !    tmp = u0
  !else
  !    tmp = 0.
  !endwhere

  tmp1 = sum( real( vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz) ) )
  tmp2 =.5_sp * sum( real( vx(1,:,:)*conjg(vx(1,:,:)) + vy(1,:,:)*conjg(vy(1,:,:)) &
                    + vz(1,:,:) * conjg(vz(1,:,:)) ) )
  tmp1 = tmp1 - tmp2
  tmp2 = (3._sp/2._sp)*u0*u0
  vx=sqrt(tmp2/tmp1)*vx
  vy=sqrt(tmp2/tmp1)*vy
  vz=sqrt(tmp2/tmp1)*vz
  ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

end subroutine initdealiased

complex(sp) function gasdev(t1,t2)
  use mconstant
  real(sp) :: t1, t2

  gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
  ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
end function gasdev
  
