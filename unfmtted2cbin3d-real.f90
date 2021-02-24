program unfmtted2cbin3dr
  use mconstant
  implicit none

  integer :: nx, numarg
  real(sp), allocatable, dimension(:,:,:) :: uu
  character(80) :: str
  integer :: c3drealdatawriter

  numarg = iarg()
  if ( numarg .ne. 2 ) then
    write(*,*)
    write(*,'(''>>> convert an unformatted data file to c binary output <<<'')')
    write(*,*)
    write(*,*) 'Usage: unfmtted2cbin3d-real.x nx filename'
    write(*,*) '       nx: resolution of the data file'
    write(*,*) '       filename: data file name'
    write(*,*)
    write(*,*) 'Stopped'
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

  numarg = c3drealdatawriter( uu, nx, nx, nx, str(1:len_trim(str)) )

  if ( numarg .ne. 0 ) then 
    write(*,*) 'Error in writing c binary data file. Stopping...'
    stop
  endif

  deallocate( uu )

  write(*,*) 'finished'

end program unfmtted2cbin3dr 
