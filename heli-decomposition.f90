program helidecomposition
  use mconstant
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz


  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,hpx,hpy,hpz,ut
  real,    allocatable, dimension(:,:,:) :: k2
  real,    allocatable, dimension(:)     :: kx,ky,kz

  character(80) :: str, fnm

  write(*,*) 
  write(*,'(''>>> Helical decomposition of velocity fields <<< '')')
  write(*,*)
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./heli-decomposition.x nx filelist '
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if


  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  deallocate(k2,kx,ky,kz)
  write(*,*) 'after heliwave'

  allocate(ut(lx1,ly,lz))

  open(20,file=fnm(1:len_trim(fnm))//'.list')

    do while ( .not. eof(20)) 
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      open(10,file='./out/ux'//str(1:len_trim(str)), status='unknown', form='unformatted')
        read(10) ux
      close(10)
      open(10,file='./out/uy'//str(1:len_trim(str)), status='unknown', form='unformatted')
        read(10) uy
      close(10)
      open(10,file='./out/uz'//str(1:len_trim(str)), status='unknown', form='unformatted')
        read(10) uz
      close(10)

      ut = ux * conjg(hpx) + uy * conjg(hpy) + uz * conjg(hpz)

      open(10, file = './out/uxph'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * hpx
      close(10)
      open(10, file = './out/uyph'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * hpy
      close(10)
      open(10, file = './out/uzph'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * hpz
      close(10)

      ut = ux * hpx + uy * hpy + uz * hpz

      open(10, file = './out/uxmh'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * conjg(hpx)
      close(10)
      open(10, file = './out/uymh'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * conjg(hpy)
      close(10)
      open(10, file = './out/uzmh'//str(1:len_trim(str)), status = 'unknown', form = 'unformatted')
        write(10) ut * conjg(hpz)
      close(10)


    end do

  close(20)

  deallocate(ux,uy,uz,hpx,hpy,hpz,ut)
  write(*,*) 'finished'

end program helidecomposition      

