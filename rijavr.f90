program rijavr

  use mconstant
  implicit none

  integer :: ii, ll, nx, lx, nfile
  real(sp), allocatable, dimension(:) :: rxxx, rxxy, rxxz
  real(sp), allocatable, dimension(:) :: ryyx, ryyy, ryyz
  real(sp), allocatable, dimension(:) :: rzzx, rzzy, rzzz

  real(sp), allocatable, dimension(:) :: meanrxxx, meanrxxy, meanrxxz
  real(sp), allocatable, dimension(:) :: meanryyx, meanryyy, meanryyz
  real(sp), allocatable, dimension(:) :: meanrzzx, meanrzzy, meanrzzz
  real(sp), allocatable, dimension(:) :: xx
  
  real(sp) :: rxxx0, rxxy0, rxxz0
  real(sp) :: ryyx0, ryyy0, ryyz0
  real(sp) :: rzzx0, rzzy0, rzzz0

  real(sp) :: meanrxxx0, meanrxxy0, meanrxxz0
  real(sp) :: meanryyx0, meanryyy0, meanryyz0
  real(sp) :: meanrzzx0, meanrzzy0, meanrzzz0

  real(sp) :: dx
  character(80) :: str, fnm

  write(*,*) 
  write(*,'(''>>> Averaging velocity correlation functions <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rijavr-s(d)p.x nx filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
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

  lx = nx / 2 
  dx = pi / lx

  allocate( rxxx(lx), ryyx(lx), rzzx(lx) )
  allocate( rxxy(lx), ryyy(lx), rzzy(lx) )
  allocate( rxxz(lx), ryyz(lx), rzzz(lx) )
  allocate( meanrxxx(lx), meanryyx(lx), meanrzzx(lx) )
  allocate( meanrxxy(lx), meanryyy(lx), meanrzzy(lx) )
  allocate( meanrxxz(lx), meanryyz(lx), meanrzzz(lx) )
  allocate( xx(lx) )

  meanrxxx0=0; meanrxxy0=0; meanrxxz0=0
  meanryyx0=0; meanryyy0=0; meanryyz0=0
  meanrzzx0=0; meanrzzy0=0; meanrzzz0=0

  meanrxxx=0; meanrxxy=0; meanrxxz=0
  meanryyx=0; meanryyy=0; meanryyz=0
  meanrzzx=0; meanrzzy=0; meanrzzz=0

  nfile = 0
  open( 30, file = fnm(1:len_trim(fnm))//'.list' )

    do while ( .not. eof(30) )

      read (30,*) str
      write( *,*) 'rij', str(1:len_trim(str))

      open(10, file = './rij'//str(1:len_trim(str)))

        read(10,*)
        read(10,*) str, rxxx0, ryyx0, rzzx0
        read(10,*)
        read(10,*) str, rxxy0, ryyy0, rzzy0
        read(10,*)
        read(10,*) str, rxxz0, ryyz0, rzzz0

        do ii = 1, nx/2
          read(10,*) xx(ii), rxxx(ii), ryyx(ii), rzzx(ii),  &
                             rxxy(ii), ryyy(ii), rzzy(ii),  &
                             rxxz(ii), ryyz(ii), rzzz(ii)
        end do
        rxxx = rxxx * rxxx0
        ryyx = ryyx * ryyx0
        rzzx = rzzx * rzzx0
        rxxy = rxxy * rxxy0
        ryyy = ryyy * ryyy0
        rzzy = rzzy * rzzy0
        rxxz = rxxz * rxxz0
        ryyz = ryyz * ryyz0
        rzzz = rzzz * rzzz0

      close(10)
      meanrxxx = meanrxxx + rxxx 
      meanryyx = meanryyx + ryyx 
      meanrzzx = meanrzzx + rzzx 
      meanrxxy = meanrxxy + rxxy 
      meanryyy = meanryyy + ryyy 
      meanrzzy = meanrzzy + rzzy 
      meanrxxz = meanrxxz + rxxz 
      meanryyz = meanryyz + ryyz 
      meanrzzz = meanrzzz + rzzz 

      meanrxxx0 = meanrxxx0 + rxxx0 
      meanryyx0 = meanryyx0 + ryyx0 
      meanrzzx0 = meanrzzx0 + rzzx0 
      meanrxxy0 = meanrxxy0 + rxxy0 
      meanryyy0 = meanryyy0 + ryyy0 
      meanrzzy0 = meanrzzy0 + rzzy0 
      meanrxxz0 = meanrxxz0 + rxxz0 
      meanryyz0 = meanryyz0 + ryyz0 
      meanrzzz0 = meanrzzz0 + rzzz0 

      nfile = nfile + 1

    end do

  close(30) 
  meanrxxx = meanrxxx / nfile 
  meanryyx = meanryyx / nfile 
  meanrzzx = meanrzzx / nfile 
  meanrxxy = meanrxxy / nfile 
  meanryyy = meanryyy / nfile 
  meanrzzy = meanrzzy / nfile 
  meanrxxz = meanrxxz / nfile 
  meanryyz = meanryyz / nfile 
  meanrzzz = meanrzzz / nfile 

  meanrxxx0 = meanrxxx0 / nfile 
  meanryyx0 = meanryyx0 / nfile 
  meanrzzx0 = meanrzzx0 / nfile 
  meanrxxy0 = meanrxxy0 / nfile 
  meanryyy0 = meanryyy0 / nfile 
  meanrzzy0 = meanrzzy0 / nfile 
  meanrxxz0 = meanrxxz0 / nfile 
  meanryyz0 = meanryyz0 / nfile 
  meanrzzz0 = meanrzzz0 / nfile 

  meanrxxx = meanrxxx / meanrxxx0
  meanryyx = meanryyx / meanryyx0
  meanrzzx = meanrzzx / meanrzzx0
  meanrxxy = meanrxxy / meanrxxy0
  meanryyy = meanryyy / meanryyy0
  meanrzzy = meanrzzy / meanrzzy0
  meanrxxz = meanrxxz / meanrxxz0
  meanryyz = meanryyz / meanryyz0
  meanrzzz = meanrzzz / meanrzzz0

  open(15,file='rijavr-'//fnm(1:len_trim(fnm))//'.dat')
    write(15, *) '# meanrxxx(1) meanryyx(1) meanrzzx(1)'
    write(15, *) '#', meanrxxx0, meanryyx0, meanrzzx0
    write(15, *) '# meanrxxy(1) meanryyy(1) meanrzzy(1)'
    write(15, *) '#', meanrxxy0, meanryyy0, meanrzzy0
    write(15, *) '# meanrxxz(1) meanryyz(1) meanrzzz(1)'
    write(15, *) '#', meanrxxz0, meanryyz0, meanrzzz0
    write(15, *) '# Lxxx    Lyyx     Lzzx'
    write(15, *) '#', dx*sum(meanrxxx), dx*sum(meanryyx), dx*sum(meanrzzx)
    write(15, *) '# Lxxy    Lyyy     Lzzy'
    write(15, *) '#', dx*sum(meanrxxy), dx*sum(meanryyy), dx*sum(meanrzzy)
    write(15, *) '# Lxxz    Lyyz     Lzzz'
    write(15, *) '#', dx*sum(meanrxxz), dx*sum(meanryyz), dx*sum(meanrzzz)
    do ii=1,lx
      write(15,'(20f9.4)') xx(ii), meanrxxx(ii), meanryyx(ii), meanrzzx(ii), & 
                                   meanrxxy(ii), meanryyy(ii), meanrzzy(ii), &
                                   meanrxxz(ii), meanryyz(ii), meanrzzz(ii)
    end do
  close(15)

  deallocate(meanrxxx,meanryyx,meanrzzx,meanrxxy)
  deallocate(meanryyy,meanrzzy,meanrxxz,meanryyz,meanrzzz)

  write(*,*) 'Finished'

end program rijavr      
