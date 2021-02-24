program curvcndalpha
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(4), allocatable, dimension(:,:,:) :: xix, xiy, xiz

  complex(dp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  real(dp),    allocatable, dimension(:,:,:) :: curv, alpha

  real(sp),    allocatable, dimension(:,:,:) :: g 
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt = 100
  real(sp), parameter :: bndalf = 8
  real(sp), parameter :: binwalf = 2*bndalf/npnt

  real(dp), dimension(npnt) :: pdfalf, curvcndalf

  real(dp) :: const, meancurv, rmscurv, meanalf, rmsalf, covca
  real(dp) :: meancurv0, rmscurv0, meanalf0, rmsalf0
  real(dp) :: delta_c, ignore_me, tmp, tmp1

  character(80) :: fnm, str, fpath

  write(*,*) 
  write(*,'(''>>>>>> Mean alpha conditioned on curvature <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./curvature-cnd-alpha.x nx filelist ndel'
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
  allocate( curv(nx,ny,nz), alpha(nx,ny,nz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  nfile = 0

  meancurv = 0._dp
  rmscurv = 0._dp

  pdfalf = 0._dp
  meanalf = 0._dp
  rmsalf = 0._dp

  covca = 0._dp

  curvcndalf = 0._dp
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xix
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xiy
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xiz
    close(10)
    write(*,*) 'after reading data files'

    xix = xix * real(g, kind = 4)
    xiy = xiy * real(g, kind = 4)
    xiz = xiz * real(g, kind = 4)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      ! Vorticity 
      g11(ii,jj,kk) = eye * ( ky(jj) * xiz(ii,jj,kk) - kz(kk) * xiy(ii,jj,kk) )
      g22(ii,jj,kk) = eye * ( kz(kk) * xix(ii,jj,kk) - kx(ii) * xiz(ii,jj,kk) )
      g33(ii,jj,kk) = eye * ( kx(ii) * xiy(ii,jj,kk) - ky(jj) * xix(ii,jj,kk) )

      ! Strain rate tensor
      g12(ii,jj,kk) = eye * ( kx(ii) * xiy(ii,jj,kk) + ky(jj) * xix(ii,jj,kk) ) / 2.
      g13(ii,jj,kk) = eye * ( kx(ii) * xiz(ii,jj,kk) + kz(kk) * xix(ii,jj,kk) ) / 2.
      g23(ii,jj,kk) = eye * ( ky(jj) * xiz(ii,jj,kk) + kz(kk) * xiy(ii,jj,kk) ) / 2.

      ! s11
      g21(ii,jj,kk) = eye * kx(ii) * xix(ii,jj,kk)
      ! s22
      g31(ii,jj,kk) = eye * ky(jj) * xiy(ii,jj,kk)
      ! s33
      g32(ii,jj,kk) = eye * kz(kk) * xiz(ii,jj,kk)

    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g23,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,g21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g32,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx
      ignore_me =  real(g11(ii,jj,kk)) ** 2 +  real(g22(ii,jj,kk)) ** 2 +  real(g33(ii,jj,kk)) ** 2
      delta_c   = aimag(g11(ii,jj,kk)) ** 2 + aimag(g22(ii,jj,kk)) ** 2 + aimag(g33(ii,jj,kk)) ** 2
      ignore_me = sqrt(ignore_me)
      delta_c   = sqrt(delta_c)

      ! direction of vorticity
      g11(ii,jj,kk) = cmplx( real( g11(ii,jj,kk) ) / (ignore_me + mytiny), &
                            aimag( g11(ii,jj,kk) ) / (delta_c + mytiny) )
      g22(ii,jj,kk) = cmplx( real( g22(ii,jj,kk) ) / (ignore_me + mytiny), &
                            aimag( g22(ii,jj,kk) ) / (delta_c + mytiny) )
      g33(ii,jj,kk) = cmplx( real( g33(ii,jj,kk) ) / (ignore_me + mytiny), &
                            aimag( g33(ii,jj,kk) ) / (delta_c + mytiny) )

      ! backup the direction in real space
      xix(ii,jj,kk) = cmplx( g11(ii,jj,kk), kind = 4 )
      xiy(ii,jj,kk) = cmplx( g22(ii,jj,kk), kind = 4 )
      xiz(ii,jj,kk) = cmplx( g33(ii,jj,kk), kind = 4 )

      alpha(2*ii-1,jj,kk) = real(g21(ii,jj,kk)) * real(g11(ii,jj,kk))**2 &
                          + real(g31(ii,jj,kk)) * real(g22(ii,jj,kk))**2 &
                          + real(g32(ii,jj,kk)) * real(g33(ii,jj,kk))**2 &
                          + 2. * real(g12(ii,jj,kk)) * real(g11(ii,jj,kk)) * real(g22(ii,jj,kk)) &
                          + 2. * real(g13(ii,jj,kk)) * real(g11(ii,jj,kk)) * real(g33(ii,jj,kk)) &
                          + 2. * real(g23(ii,jj,kk)) * real(g22(ii,jj,kk)) * real(g33(ii,jj,kk)) 

      alpha(2*ii,jj,kk) = aimag(g21(ii,jj,kk)) * aimag(g11(ii,jj,kk))**2 &
                        + aimag(g31(ii,jj,kk)) * aimag(g22(ii,jj,kk))**2 &
                        + aimag(g32(ii,jj,kk)) * aimag(g33(ii,jj,kk))**2 &
                        + 2. * aimag(g12(ii,jj,kk)) * aimag(g11(ii,jj,kk)) * aimag(g22(ii,jj,kk)) &
                        + 2. * aimag(g13(ii,jj,kk)) * aimag(g11(ii,jj,kk)) * aimag(g33(ii,jj,kk)) &
                        + 2. * aimag(g23(ii,jj,kk)) * aimag(g22(ii,jj,kk)) * aimag(g33(ii,jj,kk)) 
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

      ! Gradient of vorticity direction
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

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1 ) then
          ll = (ii+1)/2

          delta_c   =  real( g11(ll,jj,kk) ) *  real( xix(ll,jj,kk) ) &
                    +  real( g12(ll,jj,kk) ) *  real( xiy(ll,jj,kk) ) &
                    +  real( g13(ll,jj,kk) ) *  real( xiz(ll,jj,kk) )
          ignore_me =  real( g21(ll,jj,kk) ) *  real( xix(ll,jj,kk) ) &
                    +  real( g22(ll,jj,kk) ) *  real( xiy(ll,jj,kk) ) &
                    +  real( g23(ll,jj,kk) ) *  real( xiz(ll,jj,kk) )
          tmp       =  real( g31(ll,jj,kk) ) *  real( xix(ll,jj,kk) ) &
                    +  real( g32(ll,jj,kk) ) *  real( xiy(ll,jj,kk) ) &
                    +  real( g33(ll,jj,kk) ) *  real( xiz(ll,jj,kk) )
      else
          ll = ii/2

          delta_c   = aimag( g11(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                    + aimag( g12(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                    + aimag( g13(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          ignore_me = aimag( g21(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                    + aimag( g22(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                    + aimag( g23(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          tmp       = aimag( g31(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                    + aimag( g32(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                    + aimag( g33(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
      end if

      ! Curvature
      tmp1 = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)
      curv(ii,jj,kk) = tmp1

    end do
    end do
    end do

    if ( nfile .eq. 0 ) then
        meancurv0 = sum( curv ) * const
        meanalf0 = sum( alpha ) * const
        rmsalf0 = sum( alpha * alpha ) * const - meanalf0 * meanalf0
        rmsalf0 = sqrt( rmsalf0 )
    end if


    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      ignore_me = curv(ii,jj,kk)
      meancurv = meancurv + ignore_me
      rmscurv = rmscurv + ignore_me * ignore_me

      delta_c = alpha(ii,jj,kk)
      meanalf = meanalf + delta_c
      rmsalf = rmsalf + delta_c * delta_c

      covca = covca + ignore_me * delta_c

      ignore_me = ignore_me / meancurv0
      delta_c = (delta_c - meanalf0) / rmsalf0

      mm = floor( (delta_c + bndalf) / binwalf ) + 1

      if ( mm .ge. 1 .and. mm .le. npnt ) then
          pdfalf(mm) = pdfalf(mm) + 1
          curvcndalf(mm)  = curvcndalf(mm)  + ignore_me
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

  curvcndalf = curvcndalf / (pdfalf + mytiny)

  pdfalf = pdfalf / nfile * const
  write(*,*) 'Check normalizationof pdfalf:', sum(pdfalf)

  pdfalf = pdfalf / binwalf

  meanalf = meanalf / nfile * const
  rmsalf = rmsalf / nfile * const
  rmsalf = sqrt( rmsalf - meanalf * meanalf )

  write(*,*) 'meanalf is:', meanalf
  write(*,*) 'rmsalf is: ', rmsalf

  covca = covca / nfile * const
  covca = covca - meancurv * meanalf
  covca = covca / (rmscurv * rmsalf) 

  ! Ouput format for gnuplot
  open(15,file='curvature-cnd-alpha'//str(1:len_trim(str)))
    write(15,*) '# mean curvature =', meancurv
    write(15,*) '# rms curvature =', rmscurv
    write(15,*) '# mean alpha =', meanalf
    write(15,*) '# rms alpha =', rmsalf
    write(15,*) '# correlation coefficient = ', covca
    write(15,*) '# Zone T= "conditional mean of alpha", i=', npnt, ', F=point'
    do ii=1,npnt
      write(15,'(15E15.5)') ((-bndalf+(ii-.5)*binwalf)*rmsalf0+meanalf0-meanalf) / rmsalf, &
              pdfalf(ii)*rmsalf/rmsalf0, curvcndalf(ii)*meancurv0/meancurv
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(curv, alpha, xix, xiy, xiz)

  call destroyplan3d

  write(*,*) 'Finished'

end program curvcndalpha
