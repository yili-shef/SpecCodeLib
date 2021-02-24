module mheliwave

  interface heliwave
    module procedure heliwavelowmem, heliwavehighmem
  end interface heliwave

contains 

subroutine heliwavehighmem(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2, kx, ky, kz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: hpx,hpy,hpz

  real(sp) :: k, k2d, kxii, kyjj, kzkk
  integer  :: ii,jj,kk

  do kk=1,lz
  kzkk = kz(1,1,kk)
  do jj=1,ly
    kyjj = ky(1,jj,1)

    if (jj .eq. 1 .and. kk .eq. 1) then  ! on the kx axis

      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        kxii = kx(ii,1,1)
        k2d=max(sqrt(kxii*kxii+kyjj*kyjj), mytiny)
        
        hpx(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(-kyjj/k2d-eye*kxii*kzkk/k/k2d)
        hpy(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(kxii/k2d-eye*kyjj*kzkk/k/k2d)
        hpz(ii,jj,kk)=(sqrt(2._sp)/2._sp)*eye*(kyjj*kyjj+kxii*kxii)/k/k2d
      end do

    else

      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        k2d=sqrt(kzkk*kzkk+kyjj*kyjj)
        kxii = kx(ii,1,1)
        
        hpx(ii,jj,kk)=(sqrt(2._sp)/2._sp)*eye*(kyjj*kyjj+kzkk*kzkk)/k/k2d
        hpy(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(-kzkk/k2d-eye*kxii*kyjj/k/k2d)
        hpz(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(kyjj/k2d-eye*kzkk*kxii/k/k2d)
 
      end do

    end if

  end do
  end do

end subroutine heliwavehighmem

subroutine heliwavelowmem(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2
  real(sp),    dimension(lx1),       intent(in)  :: kx
  real(sp),    dimension(ly) ,       intent(in)  :: ky
  real(sp),    dimension(lz) ,       intent(in)  :: kz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: hpx,hpy,hpz

  real(sp) :: k, k2d
  integer  :: ii,jj,kk

  do kk=1,lz
  do jj=1,ly
    if (jj .eq. 1 .and. kk .eq. 1) then  ! on the kx axis
      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        k2d=max(sqrt(kx(ii)*kx(ii)+ky(jj)*ky(jj)),smallest)
        
        hpx(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(-ky(jj)/k2d-eye*kx(ii)*kz(kk)/k/k2d)
        hpy(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(kx(ii)/k2d-eye*ky(jj)*kz(kk)/k/k2d)
        hpz(ii,jj,kk)=(sqrt(2._sp)/2._sp)*eye*(ky(jj)*ky(jj)+kx(ii)*kx(ii))/k/k2d
      end do
    else
      do ii=1,lx1
        k=sqrt(k2(ii,jj,kk))
        k2d=sqrt(kz(kk)*kz(kk)+ky(jj)*ky(jj))
        
        hpx(ii,jj,kk)=(sqrt(2._sp)/2._sp)*eye*(ky(jj)*ky(jj)+kz(kk)*kz(kk))/k/k2d
        hpy(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(-kz(kk)/k2d-eye*kx(ii)*ky(jj)/k/k2d)
        hpz(ii,jj,kk)=(sqrt(2._sp)/2._sp)*(ky(jj)/k2d-eye*kz(kk)*kx(ii)/k/k2d)
 
      end do
    end if
  end do
  end do

end subroutine heliwavelowmem

end module mheliwave
