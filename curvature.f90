program curvature
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(sp), allocatable, dimension(:,:,:) :: xix, xiy, xiz, g11,g12,g13,g21,g22,g23,g31,g32,g33
  real(sp),    allocatable, dimension(:,:,:) :: g, curv
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt = 320
  real(sp), parameter :: bndcurv = 8
  real(sp), parameter :: binwcurv = 2*bndcurv/npnt
  real(dp), dimension(npnt) :: pdfcurv

  real(dp) :: const, meancurv, rmscurv
  real(sp) :: meancurv0, rmscurv0
  real(sp) :: delta_c, ignore_me, tmp
  character(80) :: fnm, str, fpath

  write(*,*) 
  write(*,'(''>>>>>> Mean rms and PDFs of Curvature of vortex lines <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./curvature.x nx filelist ndel'
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

  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), xix(lx1,ly,lz), xiy(lx1,ly,lz), xiz(lx1,ly,lz) )
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz), g21(lx1,ly,lz), g22(lx1,ly,lz) )
  allocate( g23(lx1,ly,lz), g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  allocate( curv(nx,ny,nz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  nfile = 0
  meancurv = 0._dp
  rmscurv = 0._dp
  pdfcurv = 0._dp
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)g31
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)g32
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)g33
    close(10)
    write(*,*) 'after reading data files'

    g31 = g31 * g
    g32 = g32 * g
    g33 = g33 * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      ! vorticity 
      g11(ii,jj,kk) = eye * ( ky(jj) * g33(ii,jj,kk) - kz(kk) * g32(ii,jj,kk) )
      g22(ii,jj,kk) = eye * ( kz(kk) * g31(ii,jj,kk) - kx(ii) * g33(ii,jj,kk) )
      g33(ii,jj,kk) = eye * ( kx(ii) * g32(ii,jj,kk) - ky(jj) * g31(ii,jj,kk) )
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx
      ignore_me = real(g11(ii,jj,kk)) ** 2 + real(g22(ii,jj,kk)) ** 2 + real(g33(ii,jj,kk)) ** 2
      ignore_me = sqrt(ignore_me)
      delta_c = aimag(g11(ii,jj,kk)) ** 2 + aimag(g22(ii,jj,kk)) ** 2 + aimag(g33(ii,jj,kk)) ** 2
      delta_c = sqrt(delta_c)

      ! direction of vorticity
      g11(ii,jj,kk) = cmplx( real( g11(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g11(ii,jj,kk) ) / (delta_c + mytiny_sp) )
      g22(ii,jj,kk) = cmplx( real( g22(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g22(ii,jj,kk) ) / (delta_c + mytiny_sp) )
      g33(ii,jj,kk) = cmplx( real( g33(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g33(ii,jj,kk) ) / (delta_c + mytiny_sp) )

      ! backup the direction in real space
      xix(ii,jj,kk) = g11(ii,jj,kk)
      xiy(ii,jj,kk) = g22(ii,jj,kk)
      xiz(ii,jj,kk) = g33(ii,jj,kk)
    end do
    end do
    end do

    g11 = g11 * const; g22 = g22 * const; g33 = g33 * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,g11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,g22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,g33,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      g12(ii,jj,kk) = eye * ky(jj) * g11(ii,jj,kk)
      g13(ii,jj,kk) = eye * kz(kk) * g11(ii,jj,kk)

      g21(ii,jj,kk) = eye * kx(ii) * g22(ii,jj,kk)
      g23(ii,jj,kk) = eye * kz(kk) * g22(ii,jj,kk)

      g31(ii,jj,kk) = eye * kx(ii) * g33(ii,jj,kk)
      g32(ii,jj,kk) = eye * ky(jj) * g33(ii,jj,kk)

      g11(ii,jj,kk) = eye * kx(ii) * g11(ii,jj,kk)
      g22(ii,jj,kk) = eye * ky(jj) * g22(ii,jj,kk)
      g33(ii,jj,kk) = eye * kz(kk) * g33(ii,jj,kk)
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

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx
      delta_c = real( g11(ii,jj,kk) ) * real( xix(ii,jj,kk) ) &
              + real( g12(ii,jj,kk) ) * real( xiy(ii,jj,kk) ) &
              + real( g13(ii,jj,kk) ) * real( xiz(ii,jj,kk) )
      ignore_me = real( g21(ii,jj,kk) ) * real( xix(ii,jj,kk) ) &
                + real( g22(ii,jj,kk) ) * real( xiy(ii,jj,kk) ) &
                + real( g23(ii,jj,kk) ) * real( xiz(ii,jj,kk) )
      tmp = real( g31(ii,jj,kk) ) * real( xix(ii,jj,kk) ) &
          + real( g32(ii,jj,kk) ) * real( xiy(ii,jj,kk) ) &
          + real( g33(ii,jj,kk) ) * real( xiz(ii,jj,kk) )

      curv(2*ii-1,jj,kk) = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)

      delta_c = aimag( g11(ii,jj,kk) ) * aimag( xix(ii,jj,kk) ) &
              + aimag( g12(ii,jj,kk) ) * aimag( xiy(ii,jj,kk) ) &
              + aimag( g13(ii,jj,kk) ) * aimag( xiz(ii,jj,kk) )
      ignore_me = aimag( g21(ii,jj,kk) ) * aimag( xix(ii,jj,kk) ) &
                + aimag( g22(ii,jj,kk) ) * aimag( xiy(ii,jj,kk) ) &
                + aimag( g23(ii,jj,kk) ) * aimag( xiz(ii,jj,kk) )
      tmp = aimag( g31(ii,jj,kk) ) * aimag( xix(ii,jj,kk) ) &
          + aimag( g32(ii,jj,kk) ) * aimag( xiy(ii,jj,kk) ) &
          + aimag( g33(ii,jj,kk) ) * aimag( xiz(ii,jj,kk) )

      curv(2*ii,jj,kk) = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)

    end do
    end do
    end do

    if ( nfile .eq. 0 ) then
        meancurv0 = sum( curv ) / (nx*ny*nz)
        rmscurv0  = sum( ( curv - meancurv0 )**2 ) / (nx*ny*nz)
        rmscurv0  = sqrt( rmscurv0 )

        write(*,*) 'Estimated mean of curvature: ', meancurv0
        write(*,*) 'Estimated rms of curvature:  ', rmscurv0
    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      ignore_me = curv(ii,jj,kk)
      meancurv = meancurv + ignore_me
      rmscurv = rmscurv + ignore_me * ignore_me

      delta_c = ( ignore_me - meancurv0 ) / rmscurv0

      mm = floor( (delta_c + bndcurv) / binwcurv ) + 1

      if ( mm .ge. 1 .and. mm .le. npnt) then
          pdfcurv(mm) = pdfcurv(mm) + 1
      end if

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

  pdfcurv = pdfcurv / nfile * const
  write(*,*) 'Check normalizationof pdfcurv:', sum(pdfcurv)

  pdfcurv = pdfcurv / binwcurv

  ! Ouput format for gnuplot
  open(15,file='pdf-curvature'//str(1:len_trim(str)))
    write(15,*) '# mean curvature =', meancurv
    write(15,*) '# rms curvature =', rmscurv
    write(15,*) '# Zone T= "pdf of curvature ", i=', npnt, ', F=point'
    do ii=1,npnt
      delta_c = -bndcurv + (ii - .5) * binwcurv
      write(15,'(15E15.5)') (delta_c * rmscurv0 + meancurv0 - meancurv) / rmscurv, &
                            pdfcurv(ii) * rmscurv / rmscurv0
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(curv, xix, xiy, xiz)

  call destroyplan3d

  write(*,*) 'Finished'
end program curvature
