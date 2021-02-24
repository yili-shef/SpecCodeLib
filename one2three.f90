program onetothree
  implicit none

  character(80) :: str
  integer, parameter :: nx = 256, lx = nx/2, lx1 = lx + 1
  complex(8), dimension(lx1,nx,nx) :: vx, vy, vz
  integer :: ll

  ll = iarg()
  if ( ll .ne. 1 ) stop "need one argument"
  call getarg(1,str)
  str = adjustl(str)

  open(16, file = 'MTLM_vxvyvz0256.bin', form = 'unformatted')
    read(16) vx
    read(16) vy
    read(16) vz
  close(16)

  open(16, file = './out/ux'//str(1:len_trim(str))//'.bin', form = 'unformatted')
    write(16) vx
  close(16)
  open(16, file = './out/uy'//str(1:len_trim(str))//'.bin', form = 'unformatted')
    write(16) vy
  close(16)
  open(16, file = './out/uz'//str(1:len_trim(str))//'.bin', form = 'unformatted')
    write(16) vz
  close(16)

end program onetothree      
