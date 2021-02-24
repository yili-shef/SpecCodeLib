! Will reading the energy spectrum and calculate the mean energy
! dissipation rate with given viscosity
program meaneps
  implicit none

  integer :: nx, ii, kk
  real :: eps, rnu, dat
  character(80) :: str

  write(*,*) 
  write(*,*) '>>>>>> Mean energy/helicity dissipation from energy/helicity spectrum <<<<<<'
  write(*,*) 

  ii=iargc()
  if (ii .ne. 3) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./meaneps.x nx rnu filename'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        rnu: viscosity'
          write(*,*) '        filename: name of the file for spectrum'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  call getarg(2,str)
  read(str, '(F15.8)') rnu
  ! file number string
  call getarg(3,str)
  str = adjustl(str)

  eps=0.
  open(15,file=str(1:len_trim(str)))
  do ii=1,nx/2
    read(15,*) kk, dat
    eps=eps+rnu*kk*kk*dat
  end do
  close(15)
  eps=2.*eps

  write(*,*) 'mean energy dissipation is : ', eps

end program meaneps      
