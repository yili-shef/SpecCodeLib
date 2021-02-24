!===================================================================
subroutine wavenumber(kx,ky,kz,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz

  real(sp), dimension(lx1), intent(out) :: kx
  real(sp), dimension(ly),  intent(out) :: ky
  real(sp), dimension(lz),  intent(out) :: kz

  integer :: ii, jj, kk

  do ii=1,lx1
    kx(ii)=real(ii-1)
  end do

  do jj=1,ly
    ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
    ! Ref. fftw_manual.pdf page 22
    ky(jj) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
  end do

  do kk=1,lz
    ! kz=0, 1, ..., lz/2-1, -lz/2, -lz/2+1, ..., -1 
    ! Ref. fftw_manual.pdf page 22
    kz(kk) = real(mod(kk-1+lz/2,lz)-lz/2,sp)
  end do

  return
end subroutine wavenumber
!==================================================================
