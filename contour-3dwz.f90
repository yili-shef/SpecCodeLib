program contour3dwz
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, lx, ly, lz, lx1, ii, jj, kk, ll, ndel, startpnt, endpnt
  complex, allocatable, dimension(:, :, :) :: ux, wz
  real,    allocatable, dimension(:, :, :) :: kx, ky, kz, g

  real :: delta_c, ignore_me
  character(80) :: str, fnm

  write(*,*) 
  write(*,'(''>>> 3D contour plot for filtered wz<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./contour-3dwz.x nx nfile ndel startpnt endpnt'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        nfile: file number'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        startpnt: starting point to plot the cube'
          write(*,*) '        endpnt: end point of the plotted cube'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! number of file
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! starting point for plotting
  call getarg(4,str)
  read(str,'(I20)') startpnt

  ! end point for plotting
  call getarg(5,str)
  read(str, '(I20)') endpnt

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)


  str='-'//str(1:len_trim(str))//'dx-wz-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.dat'


  lx = nx / 2; lx1 = lx + 1
  ly = nx;  lz = nx

  allocate(kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  allocate(ux(lx1,ly,lz), wz(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call fftwplan3de(nx,nx,nx)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  delta_c = ndel*2*pi / nx
  g=exp(-g*delta_c**2/24.)

  open( 10, file = './out/ux'//fnm( 1 : len_trim(fnm) ), form = 'unformatted' )
    read(10) ux   
  close(10)
  open( 10, file = './out/uy'//fnm( 1 : len_trim(fnm) ), form = 'unformatted' )
    read(10) wz
  close(10)

  wz=eye*(kx*wz-ky*ux)*g

  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  open(15, file = 'contour-3d'//str( 1 : len_trim(str) ) )
    ll = endpnt - startpnt + 1
    write(15,'(''zone t="3d enstr contour", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') ll,ll,ll
    do kk = startpnt, endpnt
    do jj = startpnt, endpnt
    do ii = startpnt, endpnt
      if ( mod(ii, 2) .eq. 0 ) then 
              ll = ii/2
              ignore_me = aimag( wz(ll,jj,kk) )
      else
              ll = (ii+1)/2
              ignore_me = real( wz(ll,jj,kk) )
      end if
      write(15,*) ii, jj, kk, ignore_me
    end do
    end do
    end do
  close(15)

  deallocate(kx, ky, kz, g, ux, wz)

  call destroyplan3d

  write(*,*) 'done'

end program contour3dwz
