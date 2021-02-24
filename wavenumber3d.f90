subroutine wavenumber3d(kx,ky,kz,k2,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp), dimension(lx1,ly,lz), intent(out) :: kx, ky, kz, k2
  integer :: ii, jj, kk

  do ii=1,lx1
    kx(ii,:,:)=real(ii-1)
  end do

  do jj=1,ly
    ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
    ! Ref. fftw_manual.pdf page 22
    ky(:,jj,:) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
  end do

  do kk=1,lz
    ! kz=0, 1, ..., lz/2-1, -lz/2, -lz/2+1, ..., -1 
    ! Ref. fftw_manual.pdf page 22
    kz(:,:,kk) = real(mod(kk-1+lz/2,lz)-lz/2,sp)
  end do

  k2 = kx * kx + ky * ky + kz * kz
  k2(1,1,1)=mytiny_sp


  return
end subroutine wavenumber3d
