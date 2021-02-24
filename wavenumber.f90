module mwavenumber

  interface wavenumber
      module procedure wavenumberlowmem2d, wavenumberhighmem2d
      module procedure wavenumberlowmem3d, wavenumberhighmem3d
      module procedure wavenumberlowmem3dnok2
      module procedure wavenumberlowmem2dnok2
      module procedure wavenumberk2only2d
      module procedure wavenumberk2only3d
      module procedure wavenumberhighmem3dnok2
  end interface wavenumber

contains
    subroutine wavenumberk2only2d(k2,lx1,ly)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly
      real(sp), dimension(lx1,ly), intent(out) :: k2
      real(sp) :: kxii, kyjj
    
      integer :: ii, jj

      do jj=1,ly
        kyjj = real(mod(jj-1+ly/2,ly)-ly/2,sp)
        do ii=1,lx1
          kxii=real(ii-1)
          k2(ii,jj)=kxii*kxii+kyjj*kyjj
        end do
      end do
      k2(1,1)=mytiny

    end subroutine wavenumberk2only2d
    
    subroutine wavenumberk2only3d(k2,lx1,ly,lz)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly,lz
      real(sp), dimension(lx1,ly,lz), intent(out) :: k2
      real(sp) :: kzkk, kyjj, kxii
    
      integer :: ii, jj, kk

      do kk=1,lz
        kzkk = real(mod(kk-1+lz/2,lz)-lz/2,sp)
        do jj=1,ly
          kyjj = real(mod(jj-1+ly/2,ly)-ly/2,sp)
          do ii=1,lx1
            kxii=real(ii-1)
            k2(ii,jj,kk)=kxii*kxii+kyjj*kyjj+kzkk*kzkk
          end do
        end do
      end do
      k2(1,1,1)=mytiny

    end subroutine wavenumberk2only3d
    
    subroutine wavenumberlowmem3dnok2(kx,ky,kz,lx1,ly,lz)
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
    end subroutine wavenumberlowmem3dnok2
    

    subroutine wavenumberlowmem3d(kx,ky,kz,k2,lx1,ly,lz)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly,lz
    
      real(sp), dimension(lx1), intent(out) :: kx
      real(sp), dimension(ly),  intent(out) :: ky
      real(sp), dimension(lz),  intent(out) :: kz
      real(sp), dimension(lx1,ly,lz), intent(out) :: k2

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
    
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx1
        k2(ii,jj,kk)=kx(ii)*kx(ii)+ky(jj)*ky(jj)+kz(kk)*kz(kk)
      end do
      end do
      end do
      k2(1,1,1)=mytiny_sp
    
    
      return
    end subroutine wavenumberlowmem3d
    
    subroutine wavenumberhighmem3d(kx,ky,kz,k2,lx1,ly,lz)
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
    
    end subroutine wavenumberhighmem3d

    subroutine wavenumberhighmem3dnok2(kx,ky,kz,lx1,ly,lz)
      use mconstant
      implicit none
    
      integer, intent(in) :: lx1,ly,lz
      real(sp), dimension(lx1,ly,lz), intent(out) :: kx, ky, kz
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
    
    end subroutine wavenumberhighmem3dnok2


    subroutine wavenumberlowmem2d(kx,ky,k2,lx1,ly)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly
    
      real(sp), dimension(lx1), intent(out) :: kx
      real(sp), dimension(ly),  intent(out) :: ky
      real(sp), dimension(lx1,ly), intent(out) :: k2

      integer :: ii, jj
    
      do ii=1,lx1
        kx(ii)=real(ii-1)
      end do
    
      do jj=1,ly
        ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
        ! Ref. fftw_manual.pdf page 22
        ky(jj) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
      end do
    
      do jj=1,ly
      do ii=1,lx1
        k2(ii,jj)=kx(ii)*kx(ii)+ky(jj)*ky(jj)
      end do
      end do
      k2(1,1)=mytiny_sp
    
    
      return
    end subroutine wavenumberlowmem2d

    subroutine wavenumberlowmem2dnok2(kx,ky,lx1,ly)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly
    
      real(sp), dimension(lx1), intent(out) :: kx
      real(sp), dimension(ly),  intent(out) :: ky

      integer :: ii, jj
    
      do ii=1,lx1
        kx(ii)=real(ii-1)
      end do
    
      do jj=1,ly
        ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
        ! Ref. fftw_manual.pdf page 22
        ky(jj) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
      end do
    
    end subroutine wavenumberlowmem2dnok2

    subroutine wavenumberhighmem2d(kx,ky,k2,lx1,ly)
      use mconstant
      implicit none

      integer, intent(in) :: lx1,ly
      real(sp), dimension(lx1,ly), intent(out) :: kx, ky, k2

      integer :: ii, jj
    
      do ii=1,lx1
        kx(ii,:)=real(ii-1)
      end do
    
      do jj=1,ly
        ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
        ! Ref. fftw_manual.pdf page 22
        ky(:,jj) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
      end do

      k2 = kx * kx + ky * ky
      k2(1,1)=mytiny_sp
    
      return
    end subroutine wavenumberhighmem2d

end module mwavenumber

