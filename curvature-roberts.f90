program croberts
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(4), allocatable, dimension(:,:,:) :: ux, uy, uz

  complex(dp), allocatable, dimension(:,:,:) :: xix, xiy, xiz
  complex(dp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  real(dp),    allocatable, dimension(:,:) :: curv, curv_wxm

  real(sp),    allocatable, dimension(:,:,:) :: g 
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt = 100
  real(sp), parameter :: bndoo = 8
  real(sp), parameter :: binwoo = 2*bndoo/npnt

  real(dp) :: const, meancurv, rmscurv, meanoo, rmsoo
  real(dp) :: meancurv0, rmscurv0, meanoo0, rmsoo0
  real(dp) :: delta_c, ignore_me, tmp, tmp1, x, y, tmpx, tmpy, tmpz

  character(80) :: fnm, str, fpath

  write(*,*) 
  write(*,'(''>>>>>> curvature distribution in Roberts flow <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./curvature-roberts.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: dns data file list'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! filter parameter
  call getarg(3,str)
  read(str, '(I20)') ndel
  ! file list string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  str='-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx
  const = 1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( ux(lx1, ly, lz), uy(lx1, ly, lz), uz(lx1, ly, lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), xix(lx1,ly,lz), xiy(lx1,ly,lz), xiz(lx1,ly,lz) )
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz), g21(lx1,ly,lz), g22(lx1,ly,lz) )
  allocate( g23(lx1,ly,lz), g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  allocate( curv(nx,ny), curv_wxm(nx,ny) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  nfile = 0

  meancurv = 0._dp
  rmscurv = 0._dp

  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)ux
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)uy
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)uz
    close(10)
    write(*,*) 'after reading data files'

    xix = ux * g
    xiy = uy * g
    xiz = uz * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      ! Vorticity 
      g11(ii,jj,kk) = eye * ( ky(jj) * xiz(ii,jj,kk) - kz(kk) * xiy(ii,jj,kk) )
      g22(ii,jj,kk) = eye * ( kz(kk) * xix(ii,jj,kk) - kx(ii) * xiz(ii,jj,kk) )
      g33(ii,jj,kk) = eye * ( kx(ii) * xiy(ii,jj,kk) - ky(jj) * xix(ii,jj,kk) )

      g12(ii,jj,kk) = eye * g11(ii,jj,kk) * ky(jj)
      g13(ii,jj,kk) = eye * g11(ii,jj,kk) * kz(kk)
      g21(ii,jj,kk) = eye * g22(ii,jj,kk) * kx(ii)
      g23(ii,jj,kk) = eye * g22(ii,jj,kk) * kz(kk)
      g31(ii,jj,kk) = eye * g33(ii,jj,kk) * kx(ii)
      g32(ii,jj,kk) = eye * g33(ii,jj,kk) * ky(jj)

      xix(ii,jj,kk) = g11(ii,jj,kk)
      xiy(ii,jj,kk) = g22(ii,jj,kk)
      xiz(ii,jj,kk) = g33(ii,jj,kk)

      g11(ii,jj,kk) = eye * g11(ii,jj,kk) * kx(ii)
      g22(ii,jj,kk) = eye * g22(ii,jj,kk) * ky(jj)
      g33(ii,jj,kk) = eye * g33(ii,jj,kk) * kz(kk)

    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xix,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiz,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

        if ( mod(ii,2) .eq. 1 ) then
            ll = (ii + 1) / 2

            ignore_me =  real(xix(ll,jj,kk)) ** 2 &
                      +  real(xiy(ll,jj,kk)) ** 2 &
                      +  real(xiz(ll,jj,kk)) ** 2

            tmpx = real(g11(ll,jj,kk)) * real(xix(ll,jj,kk)) &
                 + real(g12(ll,jj,kk)) * real(xiy(ll,jj,kk)) &
                 + real(g13(ll,jj,kk)) * real(xiz(ll,jj,kk))
            tmpy = real(g21(ll,jj,kk)) * real(xix(ll,jj,kk)) &
                 + real(g22(ll,jj,kk)) * real(xiy(ll,jj,kk)) &
                 + real(g23(ll,jj,kk)) * real(xiz(ll,jj,kk))
            tmpz = real(g31(ll,jj,kk)) * real(xix(ll,jj,kk)) &
                 + real(g32(ll,jj,kk)) * real(xiy(ll,jj,kk)) &
                 + real(g33(ll,jj,kk)) * real(xiz(ll,jj,kk))

            delta_c = real(xix(ll,jj,kk)) * tmpx &
                    + real(xiy(ll,jj,kk)) * tmpy &
                    + real(xiz(ll,jj,kk)) * tmpz

            tmpx = tmpx - real(xix(ll,jj,kk)) * delta_c / (ignore_me + mytiny)
            tmpy = tmpy - real(xiy(ll,jj,kk)) * delta_c / (ignore_me + mytiny)
            tmpz = tmpz - real(xiz(ll,jj,kk)) * delta_c / (ignore_me + mytiny)

        else
            ll = ii / 2

            ignore_me =  aimag(xix(ll,jj,kk)) ** 2 &
                      +  aimag(xiy(ll,jj,kk)) ** 2 &
                      +  aimag(xiz(ll,jj,kk)) ** 2

            tmpx = aimag(g11(ll,jj,kk)) * aimag(xix(ll,jj,kk)) &
                 + aimag(g12(ll,jj,kk)) * aimag(xiy(ll,jj,kk)) &
                 + aimag(g13(ll,jj,kk)) * aimag(xiz(ll,jj,kk))
            tmpy = aimag(g21(ll,jj,kk)) * aimag(xix(ll,jj,kk)) &
                 + aimag(g22(ll,jj,kk)) * aimag(xiy(ll,jj,kk)) &
                 + aimag(g23(ll,jj,kk)) * aimag(xiz(ll,jj,kk))
            tmpz = aimag(g31(ll,jj,kk)) * aimag(xix(ll,jj,kk)) &
                 + aimag(g32(ll,jj,kk)) * aimag(xiy(ll,jj,kk)) &
                 + aimag(g33(ll,jj,kk)) * aimag(xiz(ll,jj,kk))

            delta_c = aimag(xix(ll,jj,kk)) * tmpx &
                    + aimag(xiy(ll,jj,kk)) * tmpy &
                    + aimag(xiz(ll,jj,kk)) * tmpz

            tmpx = tmpx - aimag(xix(ll,jj,kk)) * delta_c / (ignore_me + mytiny)
            tmpy = tmpy - aimag(xiy(ll,jj,kk)) * delta_c / (ignore_me + mytiny)
            tmpz = tmpz - aimag(xiz(ll,jj,kk)) * delta_c / (ignore_me + mytiny)

        end if


        tmp1 = sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz) / (ignore_me + mytiny)
        curv(ii,jj) = tmp1

        x = (ii-1) * 2 * pi / real(nx, dp)
        y = (jj-1) * 2 * pi / real(ny, dp)
 
        curv_wxm(ii,jj) = (sqrt(-sin(x)**2*sin(y)**8+(sin(x)**4 &
                        - 2*sin(x)**2)*sin(y)**6+(sin(x)**6-4*sin(x)**4 &
                        + 4*sin(x)**2)*sin(y)**4+ (-sin(x)**8-2*sin(x)**6 &
                        + 4*sin(x)**4)*sin(y)**2))/( sqrt(sin(y)**8+ &
                        4*sin(x)**2*sin(y)**6+6*sin(x)**4*sin(y)**4 &
                        +4*sin(x)**6*sin(y)**2+sin(x)**8)+mytiny)
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(20)

  meancurv = meancurv / nfile * const
  rmscurv = rmscurv / nfile * const
  rmscurv = sqrt( rmscurv - meancurv * meancurv )

  write(*,*) 'meancurv is:', meancurv
  write(*,*) 'rmscurv is: ', rmscurv

  ignore_me = sqrt(sum( (curv - curv_wxm)**2 ) /(nx*ny))
  write(*,*) 'Mean relative error:', ignore_me

  ! Ouput format for gnuplot
  open(15,file='curvature-roberts'//str(1:len_trim(str)))
    write(15,*) '# mean curvature =', meancurv
    write(15,*) '# rms  curvature =', rmscurv
    do ii = 1, ny
    do jj = 1, nx
      write(15,'(15E15.5)') curv(ii,jj), curv_wxm(ii,jj)
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(ux, uy, uz, curv, xix, xiy, xiz, curv_wxm)

  call destroyplan3d

  write(*,*) 'Finished'

end program croberts
