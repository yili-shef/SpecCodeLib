subroutine initialize_heli(ux,uy,uz,kx,ky,kz,k2,lx1,ly,lz,iseed,hrel,akp,u0)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz,iseed
  real(sp), intent(in) :: akp,u0,hrel

  complex(sp), dimension(lx1,ly,lz), intent(out) :: ux,uy,uz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2
  real(sp),    dimension(lx1),       intent(in)  :: kx
  real(sp),    dimension(ly),        intent(in)  :: ky
  real(sp),    dimension(lz),        intent(in)  :: kz

  complex(sp), dimension(lx1,ly,lz) :: ap,am
  complex(sp) :: hpx,hpy,hpz,hmx,hmy,hmz

  real(sp),    dimension(lx1,ly,lz) :: atmp
  real(sp) :: tmp1,tmp2,k,k2d
  integer  :: ii,jj,kk,idum
  real(sp) :: ran1

  idum=iseed
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    am(ii,jj,kk)=gasdev(tmp1,tmp2)
    tmp1 = ran1(idum)
    tmp2 = ran1(idum)
    tmp2 = 2._sp*pi*tmp2
    ap(ii,jj,kk)=gasdev(tmp1,tmp2)
  end do
  end do
  end do

  do ii=1,lx1-1
    tmp1=sum(ap*conjg(ap),mask=(abs(sqrt(k2)-ii-0.5_sp*oneless).lt..5_sp))
    tmp2=sum(am*conjg(am),mask=(abs(sqrt(k2)-ii-0.5_sp*oneless).lt..5_sp))
    where(abs(sqrt(k2)-ii).lt..5_sp)
      am=am*sqrt((1._sp-hrel)*tmp1/tmp2)
    endwhere
  end do

  ! scaled so that (1/2)<u_i^2>=(3/2)*atmp^2
  do kk=1,lz
  do jj=1,ly
    if (jj .eq. 1 .and. kk .eq. 1) then  ! on the kx axis
      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        tmp1 = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * k * exp (-k*k/(akp*akp)) 
        k2d=max(sqrt(kx(ii)*kx(ii)+ky(jj)*ky(jj)),smallest)
        
        hpx=(sqrt(2._sp)/2._sp)*(-ky(jj)/k2d-eye*kx(ii)*kz(kk)/k/k2d)
        hpy=(sqrt(2._sp)/2._sp)*(kx(ii)/k2d-eye*ky(jj)*kz(kk)/k/k2d)
        hpz=(sqrt(2._sp)/2._sp)*eye*(ky(jj)*ky(jj)+kx(ii)*kx(ii))/k/k2d
        hmx=conjg(hpx)
        hmy=conjg(hpy)
        hmz=conjg(hpz)

        ux(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpx+am(ii,jj,kk)*hmx)*tmp1
        uy(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpy+am(ii,jj,kk)*hmy)*tmp1
        uz(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpz+am(ii,jj,kk)*hmz)*tmp1
      end do
    else
      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        tmp1 = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * k * exp (-k*k/(akp*akp)) 
        k2d=sqrt(kz(kk)*kz(kk)+ky(jj)*ky(jj))
        
        hpx=(sqrt(2._sp)/2._sp)*eye*(ky(jj)*ky(jj)+kz(kk)*kz(kk))/k/k2d
        hpy=(sqrt(2._sp)/2._sp)*(-kz(kk)/k2d-eye*kx(ii)*ky(jj)/k/k2d)
        hpz=(sqrt(2._sp)/2._sp)*(ky(jj)/k2d-eye*kz(kk)*kx(ii)/k/k2d)
 
        hmx=conjg(hpx)
        hmy=conjg(hpy)
        hmz=conjg(hpz)

        ux(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpx+am(ii,jj,kk)*hmx)*tmp1
        uy(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpy+am(ii,jj,kk)*hmy)*tmp1
        uz(ii,jj,kk)=sqrt(3._sp/2._sp)*(ap(ii,jj,kk)*hpz+am(ii,jj,kk)*hmz)*tmp1
      end do
    end if
  end do
  end do

  call symmetrize(ux,k2,lx1,ly,lz)
  call symmetrize(uy,k2,lx1,ly,lz)
  call symmetrize(uz,k2,lx1,ly,lz) 

  atmp=ux*conjg(ux)+uy*conjg(uy)+uz*conjg(uz)
  atmp(1,:,:)=.5_sp*atmp(1,:,:)
  tmp1=sum(atmp)
  tmp2=(3._sp/2._sp)*u0*u0
  ux=sqrt(tmp2/tmp1)*ux
  uy=sqrt(tmp2/tmp1)*uy
  uz=sqrt(tmp2/tmp1)*uz
  ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.
contains

  complex(sp) function gasdev(t1,t2)
    real(sp) :: t1, t2

    gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
    ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
  end function gasdev

end subroutine initialize_heli
