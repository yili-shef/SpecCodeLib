program runpadd
  use mconstant
  implicit none

  integer, parameter :: nx = 256, nxb = 384
  integer, parameter :: lx = nx / 2, lx1 = lx + 1, ly = nx, lz = nx
  integer, parameter :: lxb = nxb / 2, lxb1 = lxb + 1, lyb = nxb, lzb = nxb
  complex(sp), dimension(lx1,  ly,  lz ) :: ux
  complex(sp), dimension(lxb1, lyb, lzb) :: uxb

  open(15, file = 'ux212.dat', form = 'unformatted')
    read(15) ux
  close(15)

  call padd(ux, uxb, lx1, ly, lz, lxb1, lyb, lzb)

  open(15, file = './out/ux1.dat', form = 'unformatted')
    write(15) uxb
  close(15)

  open(15, file = 'uy212.dat', form = 'unformatted')
    read(15) ux
  close(15)

  call padd(ux, uxb, lx1, ly, lz, lxb1, lyb, lzb)

  open(15, file = './out/uy1.dat', form = 'unformatted')
    write(15) uxb
  close(15)

  open(15, file = 'uz212.dat', form = 'unformatted')
    read(15) ux
  close(15)

  call padd(ux, uxb, lx1, ly, lz, lxb1, lyb, lzb)

  open(15, file = './out/uz1.dat', form = 'unformatted')
    write(15) uxb
  close(15)

end program runpadd      
