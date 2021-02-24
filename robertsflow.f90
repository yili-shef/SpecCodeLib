program roberts
  use mconstant
  use mfftwplan3d
  implicit none

  integer, parameter :: nx = 256, ny = nx, nz = nx
  integer, parameter :: lx = nx/2, lx1 = lx + 1, ly = ny, lz = nz
  real(sp), parameter :: const = 1./(nx*ny*nz)

  real(sp),    dimension(nx , ny, nz) :: uxr, uyr, uzr
  complex(sp), dimension(lx1, ly, lz) :: ux,  uy,  uz

  real(sp) :: dx, dy, x, y, ignore_me
  integer :: ii, jj, kk

  dx = 2*pi / nx; dy = 2*pi / ny

  do kk = 1, nz
    do jj = 1, ny
      do ii = 1, nx
        x = (ii-1) * dx
        y = (jj-1) * dy
        uxr(ii,jj,kk) =   sin(x) * cos(y)
        uyr(ii,jj,kk) = - cos(x) * sin(y)
        uzr(ii,jj,kk) = sqrt(2.) * sin(x) * sin(y)
      end do
    end do
  end do

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  ux(1:lx,:,:) = cmplx( uxr(1:nx:2,:,:), uxr(2:nx:2,:,:) )
  uy(1:lx,:,:) = cmplx( uyr(1:nx:2,:,:), uyr(2:nx:2,:,:) )
  uz(1:lx,:,:) = cmplx( uzr(1:nx:2,:,:), uzr(2:nx:2,:,:) )

  ux = ux * const; uy = uy * const; uz = uz * const
  call rfftwnd_f77_one_real_to_complex(r2c3d,ux,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,uy,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,uz,ignore_me)

  open(15, file = './out/ux0000.dat', form = 'unformatted')
    write(15) ux
  close(15)
  open(15, file = './out/uy0000.dat', form = 'unformatted')
    write(15) uy
  close(15)
  open(15, file = './out/uz0000.dat', form = 'unformatted')
    write(15) uz
  close(15)

  call destroyplan3d

  write(*,*) 'Finished'

end program roberts      
