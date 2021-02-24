program tcndc
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(4), allocatable, dimension(:,:,:) :: xix, xiy, xiz

  complex(dp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  complex(dp), allocatable, dimension(:,:,:) :: binormalx, binormaly, binormalz
  real(dp),    allocatable, dimension(:,:,:) :: curv, torsion

  real(sp),    allocatable, dimension(:,:,:) :: g 
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt = 100
  real(sp), parameter :: bndcurv = 12, bndtor = 6
  real(sp), parameter :: binwcurv = bndcurv/npnt, binwtor = 2.*bndtor / npnt

  real(dp), dimension(npnt) :: pdfcurv, torcndc

  real(dp) :: const, meancurv, rmscurv, meantor, rmstor, covct
  real(dp) :: meancurv0, rmscurv0, meantor0, rmstor0
  real(dp) :: delta_c, ignore_me, tmp, tmp1, normalx, normaly, normalz

  character(80) :: fnm, str, fpath

  write(*,*) 
  write(*,'(''>>>>>> Mean Torsion of vortex lines conditioned on curvature <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./torsion-cnd-curvature.x nx filelist ndel'
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
  allocate( binormalx(lx1,ly,lz), binormaly(lx1,ly,lz), binormalz(lx1,ly,lz) )
  allocate( curv(nx,ny,nz), torsion(nx,ny,nz) )
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

  meantor = 0._dp
  rmstor = 0._dp

  covct = 0._dp

  torcndc = 0._dp
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

    xix = xix * g
    xiy = xiy * g
    xiz = xiz * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      ! vorticity 
      g11(ii,jj,kk) = eye * ( ky(jj) * xiz(ii,jj,kk) - kz(kk) * xiy(ii,jj,kk) )
      g22(ii,jj,kk) = eye * ( kz(kk) * xix(ii,jj,kk) - kx(ii) * xiz(ii,jj,kk) )
      g33(ii,jj,kk) = eye * ( kx(ii) * xiy(ii,jj,kk) - ky(jj) * xix(ii,jj,kk) )
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

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1 ) then
          ll = (ii+1)/2
          delta_c = real( g11(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
                  + real( g12(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
                  + real( g13(ll,jj,kk) ) * real( xiz(ll,jj,kk) )
          ignore_me = real( g21(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
                    + real( g22(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
                    + real( g23(ll,jj,kk) ) * real( xiz(ll,jj,kk) )
          tmp = real( g31(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
              + real( g32(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
              + real( g33(ll,jj,kk) ) * real( xiz(ll,jj,kk) )
      else
          ll = ii/2

          delta_c = aimag( g11(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                  + aimag( g12(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                  + aimag( g13(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          ignore_me = aimag( g21(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                    + aimag( g22(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                    + aimag( g23(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          tmp = aimag( g31(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
              + aimag( g32(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
              + aimag( g33(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
      end if

      tmp1 = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)
      curv(ii,jj,kk) = tmp1

      normalx = delta_c / tmp1
      normaly = ignore_me / tmp1
      normalz = tmp / tmp1

      if ( mod(ii,2) .eq. 1 ) then

          ll = (ii+1)/2
          binormalx(ll,jj,kk) = cmplx( real(xiy(ll,jj,kk)) * normalz &
                              - real(xiz(ll,jj,kk)) * normaly, 0.)
          binormaly(ll,jj,kk) = cmplx( real(xiz(ll,jj,kk)) * normalx &
                              - real(xix(ll,jj,kk)) * normalz, 0.)
          binormalz(ll,jj,kk) = cmplx( real(xix(ll,jj,kk)) * normaly &
                              - real(xiy(ll,jj,kk)) * normalx, 0.)
      else

          ll = ii / 2
          binormalx(ll,jj,kk) = cmplx(real(binormalx(ll,jj,kk)), aimag(xiy(ll,jj,kk)) * normalz &
                              - aimag(xiz(ll,jj,kk)) * normaly )
          binormaly(ll,jj,kk) = cmplx(real(binormaly(ll,jj,kk)), aimag(xiz(ll,jj,kk)) * normalx &
                              - aimag(xix(ll,jj,kk)) * normalz )
          binormalz(ll,jj,kk) = cmplx(real(binormalz(ll,jj,kk)), aimag(xix(ll,jj,kk)) * normaly &
                              - aimag(xiy(ll,jj,kk)) * normalx )
      end if

    end do
    end do
    end do

    binormalx = binormalx * const
    binormaly = binormaly * const
    binormalz = binormalz * const

    call rfftwnd_f77_one_real_to_complex(r2c3d,binormalx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,binormaly,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,binormalz,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      g11(ii,jj,kk) = eye * kx(ii) * binormalx(ii,jj,kk)
      g12(ii,jj,kk) = eye * ky(jj) * binormalx(ii,jj,kk)
      g13(ii,jj,kk) = eye * kz(kk) * binormalx(ii,jj,kk)

      g21(ii,jj,kk) = eye * kx(ii) * binormaly(ii,jj,kk)
      g22(ii,jj,kk) = eye * ky(jj) * binormaly(ii,jj,kk)
      g23(ii,jj,kk) = eye * kz(kk) * binormaly(ii,jj,kk)

      g31(ii,jj,kk) = eye * kx(ii) * binormalz(ii,jj,kk)
      g32(ii,jj,kk) = eye * ky(jj) * binormalz(ii,jj,kk)
      g33(ii,jj,kk) = eye * kz(kk) * binormalz(ii,jj,kk)

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

    call rfftwnd_f77_one_complex_to_real(c2r3d,binormalx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,binormaly,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,binormalz,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1 ) then
          ll = (ii+1)/2
          delta_c = real( g11(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
                  + real( g12(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
                  + real( g13(ll,jj,kk) ) * real( xiz(ll,jj,kk) )
          ignore_me = real( g21(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
                    + real( g22(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
                    + real( g23(ll,jj,kk) ) * real( xiz(ll,jj,kk) )
          tmp = real( g31(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
              + real( g32(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
              + real( g33(ll,jj,kk) ) * real( xiz(ll,jj,kk) )

          normalx = real( binormaly(ll,jj,kk) ) * real( xiz(ll,jj,kk) ) &
                  - real( binormalz(ll,jj,kk) ) * real( xiy(ll,jj,kk) )

          normaly = real( binormalz(ll,jj,kk) ) * real( xix(ll,jj,kk) ) &
                  - real( binormalx(ll,jj,kk) ) * real( xiz(ll,jj,kk) ) 

          normalz = real( binormalx(ll,jj,kk) ) * real( xiy(ll,jj,kk) ) &
                  - real( binormaly(ll,jj,kk) ) * real( xix(ll,jj,kk) ) 

      else
          ll = ii/2

          delta_c = aimag( g11(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                  + aimag( g12(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                  + aimag( g13(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          ignore_me = aimag( g21(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                    + aimag( g22(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                    + aimag( g23(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )
          tmp = aimag( g31(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
              + aimag( g32(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
              + aimag( g33(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) )

          normalx = aimag( binormaly(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) ) &
                  - aimag( binormalz(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) )

          normaly = aimag( binormalz(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) &
                  - aimag( binormalx(ll,jj,kk) ) * aimag( xiz(ll,jj,kk) ) 

          normalz = aimag( binormalx(ll,jj,kk) ) * aimag( xiy(ll,jj,kk) ) &
                  - aimag( binormaly(ll,jj,kk) ) * aimag( xix(ll,jj,kk) ) 

      end if

      tmp1 = normalx * delta_c + normaly * ignore_me + normalz * tmp
      torsion(ii,jj,kk) = -tmp1
    end do
    end do
    end do

    if ( nfile .eq. 0 ) then
        meancurv0 = sum( curv ) * const
        meantor0 = sum( torsion ) * const
        rmstor0 = sum( torsion * torsion ) * const - meantor0 * meantor0
        rmstor0 = sqrt( rmstor0 )
    end if


    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      ignore_me = curv(ii,jj,kk)
      meancurv = meancurv + ignore_me
      rmscurv = rmscurv + ignore_me * ignore_me

      delta_c = torsion(ii,jj,kk)
      meantor = meantor + delta_c
      rmstor = rmstor + delta_c * delta_c

      covct = covct + ignore_me * delta_c

      ignore_me = ignore_me / meancurv0
      delta_c = delta_c / rmstor0

      mm = floor( ignore_me / binwcurv ) + 1

      if ( mm .ge. 1 .and. mm .le. npnt ) then
          pdfcurv(mm) = pdfcurv(mm) + 1
          torcndc(mm)  = torcndc(mm)  + delta_c
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

  torcndc = torcndc / (pdfcurv + mytiny)

  pdfcurv = pdfcurv / nfile * const
  write(*,*) 'Check normalizationof pdfcurv:', sum(pdfcurv)

  pdfcurv = pdfcurv / binwcurv

  meantor = meantor / nfile * const
  rmstor = rmstor / nfile * const
  rmstor = sqrt( rmstor - meantor * meantor )

  write(*,*) 'meantor is:', meantor
  write(*,*) 'rmstor is: ', rmstor

  covct = covct / nfile * const
  covct = covct - meancurv * meantor
  covct = covct / (rmscurv * rmstor) 

  ! Ouput format for gnuplot
  open(15,file='torsion-cnd-curvature'//str(1:len_trim(str)))
    write(15,*) '# mean curvature =', meancurv
    write(15,*) '# rms curvature =', rmscurv
    write(15,*) '# mean torsion =', meantor
    write(15,*) '# rms torsion =', rmstor
    write(15,*) '# correlation coefficient = ', covct
    write(15,*) '# Zone T= "conditional mean of torsion", i=', npnt, ', F=point'
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5) * binwcurv * meancurv0 / meancurv, pdfcurv(ii)*meancurv/meancurv0, &
                            torcndc(ii)*rmstor0/rmstor
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(curv, torsion, xix, xiy, xiz)
  deallocate(binormalx, binormaly, binormalz)

  call destroyplan3d

  write(*,*) 'Finished'

end program tcndc
