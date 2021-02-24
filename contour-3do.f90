program contour3do
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, lx, ly, lz, lx1, ii, jj, kk, ll, ndel, startpnt, endpnt
  complex, allocatable, dimension(:, :, :) :: wx, wy, wz, enst
  real,    allocatable, dimension(:)       :: kx, ky, kz
  real,    allocatable, dimension(:, :, :) :: g

  real :: delta_c, ignore_me
  character(80) :: str, fnm

  write(*,*) 
  write(*,'(''>>> 3D contour plot for filtered enstrophy with given threshold value <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./contour-3do.x nx nfile ndel startpnt endpnt'
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


  str='-'//str(1:len_trim(str))//'dx-enst-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.dat'


  lx = nx / 2; lx1 = lx + 1
  ly = nx;  lz = nx

  allocate(kx(lx1), ky(ly), kz(lz))
  allocate(g(lx1,ly,lz), enst(lx1,ly,lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call fftwplan3de(nx,nx,nx)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  delta_c = ndel*2*pi / nx
  g=exp(-g*delta_c**2/24.)

  open( 10, file = './out/ux'//fnm( 1 : len_trim(fnm) ), form = 'unformatted' )
    read(10) enst
  close(10)
  open( 10, file = './out/uy'//fnm( 1 : len_trim(fnm) ), form = 'unformatted' )
    read(10) wz   
  close(10)
  open( 10, file = './out/uz'//fnm( 1 : len_trim(fnm) ), form = 'unformatted' )
    read(10) wy
  close(10)

  ! enst = ux; wz = uy; wy = uz
  ! staggered use to save memory
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wx(ii,jj,kk)=eye*(ky(jj)*  wy(ii,jj,kk)-kz(kk)*  wz(ii,jj,kk))
    wy(ii,jj,kk)=eye*(kz(kk)*enst(ii,jj,kk)-kx(ii)*  wy(ii,jj,kk))
    wz(ii,jj,kk)=eye*(kx(ii)*  wz(ii,jj,kk)-ky(jj)*enst(ii,jj,kk))
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  enst = .5 * cmplx( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz), &
           aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) &
           )


  open(15, file = 'contour-3d'//str( 1 : len_trim(str) ) )
    ll = endpnt - startpnt + 1
    write(15,'(''zone t="3d enstr contour", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') ll,ll,ll
    do kk = startpnt, endpnt
    do jj = startpnt, endpnt
    do ii = startpnt, endpnt
      if ( mod(ii, 2) .eq. 0 ) then 
              ll = ii/2
              ignore_me = aimag( enst(ll,jj,kk) )
      else
              ll = (ii+1)/2
              ignore_me = real( enst(ll,jj,kk) )
      end if
      write(15,*) ii, jj, kk, ignore_me
    end do
    end do
    end do
  close(15)

  deallocate(kx, ky, kz, g, enst, wx, wy, wz)

  call destroyplan3d

  write(*,*) 'done'

end program contour3do      
