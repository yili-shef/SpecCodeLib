program unfmtted2bin3dr
  use mconstant
  implicit none

  integer :: nx, numarg
  real(sp), allocatable, dimension(:,:,:) :: uu
  character(80) :: str

  numarg = iarg()
  if ( numarg .ne. 2 ) then
    write(*,'(''>>> convert an unformatted data file to binary <<<'')')
    write(*,*) 'Usage: unfmtted2bin3d-real.x nx filename'
    write(*,*) '       nx: resolution of the data file'
    write(*,*) '       filename: data file name'
    stop
  end if

  call getarg(1,str)
  read(str, '(I20)') nx
  call getarg(2,str)
  str = adjustl(str)

  allocate( uu(nx,nx,nx) )

  open(16, file = str(1:len_trim(str)), form = 'unformatted')
    read(16) uu
  close(16)

  open(16, file = 'bin-'//str(1:len_trim(str)), form = 'binary')
    write(16) uu
  close(16)

  deallocate( uu )

  write(*,*) 'finished'

end program unfmtted2bin3dr 
