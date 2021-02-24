program rescale
  use mconstant
  implicit none

  integer :: nx, ny, nz, lx1, lx, ly, lz

  complex(sp), allocatable, dimension(:,:,:) :: vx, vy, vz
  real(sp),    allocatable, dimension(:,:,:) :: ek
  real(sp) :: urms, et
  character(80) :: str, infile, outfile

  if ( iargc() .ne. 4) stop "usage: rescale.x nx infile# urms outfile#"

  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(2,infile)
  infile = adjustl(infile)

  call getarg(3,str)
  read(str, '(F15.3)') urms

  call getarg(4,outfile)
  outfile = adjustl(outfile) 

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz
  lx1 = lx + 1

  allocate( vx(lx1,ly,lz), vy(lx1,ly,lz), vz(lx1,ly,lz) )
  allocate( ek(lx1,ly,lz) )

  open(15, file = './out/ux'//trim(infile)//'.dat', form = 'unformatted')
    read(15) vx
  close(15)
  open(15, file = './out/uy'//trim(infile)//'.dat', form = 'unformatted')
    read(15) vy
  close(15)
  open(15, file = './out/uz'//trim(infile)//'.dat', form = 'unformatted')
    read(15) vz
  close(15)

  et = 3*urms*urms/2

  ek = real( vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz), sp )
  ek(1,:,:)=0.5_sp*ek(1,:,:)

  et = sqrt( et/sum(ek) )

  vx = vx * et; vy = vy * et; vz = vz * et

  open(15, file = './out/ux'//trim(outfile)//'.dat', form = 'unformatted')
    write(15) vx
  close(15)
  open(15, file = './out/uy'//trim(outfile)//'.dat', form = 'unformatted')
    write(15) vy
  close(15)
  open(15, file = './out/uz'//trim(outfile)//'.dat', form = 'unformatted')
    write(15) vz
  close(15)

  deallocate(vx, vy, vz, ek)

  write(*,*) 'rescale.x done'

end program rescale
