program fimeans
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, ly, lz, lx1, ii, jj, kk, ll
  real(dp) :: fifi, fieps
  real(dp) :: meanfi1, rmsfi1, skfi1, flfi1
  real(dp) :: meanfi2, rmsfi2, skfi2, flfi2
  real(dp) :: meanfi3, rmsfi3, skfi3, flfi3

  real(sp) :: rnu, pr, ignore_me, tmp1, tmp2

  complex(sp), allocatable, dimension(:,:,:) :: fi, gfi
  real(sp),    allocatable, dimension(:,:,:) :: k2, kx, ky, kz
  character(80) :: str, fnm


  write(*,*) 
  write(*,'(''       >>>>>> Mean parameters <<<<<<'')')
  write(*,'('' Mean dissipation, gradient of scalar statistics'')')
  write(*,*) 

  ii=iargc()
  if (iargc() .ne. 4) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./fimeans.x nx rnu pr filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        rnu: viscosity'
          write(*,*) '        pr: prandtl number'
          write(*,*) '        filelist: list of data files, *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(2,str)
  read(str, '(F15.8)') rnu

  call getarg(3,str)
  read(str, '(F15.8)') pr

  ! file list string
  call getarg(4,fnm)
  fnm = adjustl(fnm)

  ny = nx; nz = nx;
  lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz

  allocate( fi(lx1,ly,lz), gfi(lx1,ly,lz) )
  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz), k2(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open( 14, file = 'fifieps-'//fnm(1:len_trim(fnm))//'.dat' )
  write(14, '(''# nfile fifi fieps'')')
  open( 19, file = 'meanfi-'//fnm(1:len_trim(fnm))//'.dat' )
  write(19,'(''# nfile fi1 fi2 fi3'')')
  open( 15, file = 'rmsfi-'//fnm(1:len_trim(fnm))//'.dat' )
  write(15,'(''# nfile fi1 fi2 fi3'')')
  open( 16, file = 'skewfi-'//fnm(1:len_trim(fnm))//'.dat' )
  write(16,'(''# nfile fi1 fi2 fi3'')')
  open( 17, file = 'flatfi-'//fnm(1:len_trim(fnm))//'.dat' )
  write(17,'(''# nfile fi1 fi2 fi3'')')

  open(20, file = fnm(1:len_trim(fnm))//'.list')

  ll = 0
  do while ( .true. )

    read(20, *, end=999) str
    write(*,*) 'data file: ', str(1:len_trim(str))

    open(10, file = './out/phi'//str(1:len_trim(str)), form = 'unformatted')
      read(10) fi
    close(10)

    gfi = fi*conjg(fi)
    gfi(1,:,:) = 0.5_sp * gfi(1,:,:)

    fifi = sum(real(gfi))
    fieps = (rnu/pr) * sum( k2 * real(gfi) )

    ! gfi1 
    gfi = eye * kx * fi
    call rfftwnd_f77_one_complex_to_real(c2r3d,gfi,ignore_me)

    meanfi1 = 0.d0
    rmsfi1 = 0.d0
    skfi1 = 0.d0
    flfi1 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( gfi( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( gfi( ii/2,jj,kk ) )
      end if
      meanfi1 = meanfi1 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsfi1 = rmsfi1 + tmp1
      tmp1 = tmp1 * ignore_me
      skfi1 = skfi1 + tmp1
      tmp1 = tmp1 * ignore_me
      flfi1 = flfi1 + tmp1
    end do
    end do
    end do
    tmp2 = 1._sp/(nx*ny*nz)
    meanfi1 = meanfi1 * tmp2

    rmsfi1 = rmsfi1 * tmp2
    flfi1 = flfi1 * tmp2
    skfi1 = skfi1 * tmp2
    flfi1 = flfi1 + 6.*meanfi1**2*rmsfi1 - 4*meanfi1*skfi1 &
                   - 3.*meanfi1**4
    skfi1 = skfi1 - 3*rmsfi1*meanfi1 + 2*meanfi1**3

    rmsfi1 = sqrt( rmsfi1 - meanfi1 * meanfi1 )
    skfi1 = skfi1 / rmsfi1**3
    flfi1 = flfi1 / rmsfi1**4
    
    ! gfi2
    gfi = eye * ky * fi
    call rfftwnd_f77_one_complex_to_real(c2r3d,gfi,ignore_me)

    meanfi2 = 0.d0
    rmsfi2 = 0.d0
    skfi2 = 0.d0
    flfi2 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( gfi( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( gfi( ii/2,jj,kk ) )
      end if
      meanfi2 = meanfi2 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsfi2 = rmsfi2 + tmp1
      tmp1 = tmp1 * ignore_me
      skfi2 = skfi2 + tmp1
      tmp1 = tmp1 * ignore_me
      flfi2 = flfi2 + tmp1
    end do
    end do
    end do
    meanfi2 = meanfi2 * tmp2

    rmsfi2 = rmsfi2 * tmp2
    flfi2 = flfi2 * tmp2
    skfi2 = skfi2 * tmp2
    flfi2 = flfi2 + 6.*meanfi2**2*rmsfi2 - 4*meanfi2*skfi2 &
                   - 3.*meanfi2**4
    skfi2 = skfi2 - 3*rmsfi2*meanfi2 + 2*meanfi2**3

    rmsfi2 = sqrt( rmsfi2 - meanfi2 * meanfi2 )
    skfi2 = skfi2 / rmsfi2**3
    flfi2 = flfi2 / rmsfi2**4
    
    
    ! gfi3
    gfi = eye * kz * fi
    call rfftwnd_f77_one_complex_to_real(c2r3d,gfi,ignore_me)

    meanfi3 = 0.d0
    rmsfi3 = 0.d0
    skfi3 = 0.d0
    flfi3 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( gfi( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( gfi( ii/2,jj,kk ) )
      end if
      meanfi3 = meanfi3 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsfi3 = rmsfi3 + tmp1
      tmp1 = tmp1 * ignore_me
      skfi3 = skfi3 + tmp1
      tmp1 = tmp1 * ignore_me
      flfi3 = flfi3 + tmp1
    end do
    end do
    end do
    meanfi3 = meanfi3 * tmp2

    rmsfi3 = rmsfi3 * tmp2
    flfi3 = flfi3 * tmp2
    skfi3 = skfi3 * tmp2
    flfi3 = flfi3 + 6.*meanfi3**2*rmsfi3 - 4*meanfi3*skfi3 &
                   - 3.*meanfi3**4
    skfi3 = skfi3 - 3*rmsfi3*meanfi3 + 2*meanfi3**3

    rmsfi3 = sqrt( rmsfi3 - meanfi3 * meanfi3 )
    skfi3 = skfi3 / rmsfi3**3
    flfi3 = flfi3 / rmsfi3**4

    write(14, '(I5, 20E12.3)') ll, fifi, fieps
    write(19, '(I5, 20E12.3)') ll, meanfi1,meanfi2,meanfi3
    write(15, '(I5, 20E12.3)') ll, rmsfi1,rmsfi2,rmsfi3
    write(16, '(I5, 20E12.3)') ll, skfi1,skfi2,skfi3
    write(17, '(I5, 20E12.3)') ll, flfi1,flfi2,flfi3

    ll = ll + 1

  end do
  999 continue


  close(14)
  close(15)
  close(16)
  close(17)
  close(19)


  deallocate(fi,gfi,kx,ky,kz,k2)

  call destroyplan3d

  write(*,*) 'fimeans.x done'

end program fimeans 
