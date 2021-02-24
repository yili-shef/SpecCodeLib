program means
  use mconstant
  use mwavenumber
  use mfftw3
  implicit none

  integer :: nx, ny, nz, lx, ly, lz, lx1, ii, jj, kk, ll
  real(dp) :: ek, eps
  real(dp) :: meana11, rmsa11, ska11, fla11
  real(dp) :: meana22, rmsa22, ska22, fla22
  real(dp) :: meana33, rmsa33, ska33, fla33
  real(dp) :: meana12, rmsa12, ska12, fla12
  real(dp) :: meana13, rmsa13, ska13, fla13
  real(dp) :: meana21, rmsa21, ska21, fla21
  real(dp) :: meana23, rmsa23, ska23, fla23
  real(dp) :: meana31, rmsa31, ska31, fla31
  real(dp) :: meana32, rmsa32, ska32, fla32
  real(dp) :: meanwx, rmswx, skwx, flwx
  real(dp) :: meanwy, rmswy, skwy, flwy
  real(dp) :: meanwz, rmswz, skwz, flwz

  real(sp) :: rnu, ignore_me, tmp1, tmp2, tmp3, tmp4
  complex(sp) :: a12, a13, a21, a23, a31, a32

  integer(8) :: aijc2r, uxc2r, uyc2r, uzc2r

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, aij
  real(sp),    allocatable, dimension(:,:,:) :: k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, fnm


  write(*,*) 
  write(*,'(''       >>>>>> Mean parameters <<<<<<'')')
  write(*,'('' Mean dissipation, TKE, gradient, vorticity statistics'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 3) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./means.x nx rnu filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        rnu: viscosity'
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
  ! file list string
  call getarg(3,fnm)
  fnm = adjustl(fnm)

  ny = nx; nz = nx;
  lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz), aij(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz), k2(lx1,ly,lz) )

  !call fftwplan3de(nx,ny,nz)
  call dfftwplan3dc2r(aij,nx,aijc2r)
  call dfftwplan3dc2r(ux,nx,uxc2r)
  call dfftwplan3dc2r(uy,nx,uyc2r)
  call dfftwplan3dc2r(uz,nx,uzc2r)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open( 14, file = 'tkeeps-'//fnm(1:len_trim(fnm))//'.dat' )
  write(14, '(''# nfile tke eps'')')
  open( 19, file = 'meanaij-'//fnm(1:len_trim(fnm))//'.dat' )
  write(19,'(''# nfile a11 a22 a33 a12 a13 a21 a23 a31 a32'')')
  open( 15, file = 'rmsaij-'//fnm(1:len_trim(fnm))//'.dat' )
  write(15,'(''# nfile a11 a22 a33 a12 a13 a21 a23 a31 a32'')')
  open( 16, file = 'skewaij-'//fnm(1:len_trim(fnm))//'.dat' )
  write(16,'(''# nfile a11 a22 a33 a12 a13 a21 a23 a31 a32'')')
  open( 17, file = 'flataij-'//fnm(1:len_trim(fnm))//'.dat' )
  write(17,'(''# nfile a11 a22 a33 a12 a13 a21 a23 a31 a32'')')
  open( 18, file = 'omegameans-'//fnm(1:len_trim(fnm))//'.dat' )
  write(18, '(''# nfile meanx meany meanz rmsx rmsy rmsz skewx skewy skewz flatx flaty flatz'')')

  open(20, file = fnm(1:len_trim(fnm))//'.list')

  ll = 0
  do while ( .not. eof(20) )

    read(20,*) str
    write(*,*) 'data file: ', str(1:len_trim(str))

    open(10, file = './out/ux'//str(1:len_trim(str)), form = 'unformatted')
      read(10) ux
    close(10)
    open(10, file = './out/uy'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uy
    close(10)
    open(10, file = './out/uz'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uz
    close(10)

    !write(*,*) 'after reading data'
 
    aij = ux*conjg(ux) + uy*conjg(uy) + uz*conjg(uz)
    aij(1,:,:) = 0.5_sp * aij(1,:,:)

    ek = sum(real(aij))
    eps = 2. * rnu * sum( k2 * real(aij) )

    ! a11 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kx(ii) * ux(ii,jj,kk)
    end do
    end do
    end do
    call dfftw_execute(aijc2r)
    !write(*,*) 'after first fftw'

    meana11 = 0.d0
    rmsa11 = 0.d0
    ska11 = 0.d0
    fla11 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana11 = meana11 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa11 = rmsa11 + tmp1
      tmp1 = tmp1 * ignore_me
      ska11 = ska11 + tmp1
      tmp1 = tmp1 * ignore_me
      fla11 = fla11 + tmp1
    end do
    end do
    end do
    tmp2 = 1./(nx*ny*nz)
    meana11 = meana11 * tmp2

    rmsa11 = rmsa11 * tmp2
    fla11 = fla11 * tmp2
    ska11 = ska11 * tmp2
    fla11 = fla11 + 6.*meana11**2*rmsa11 - 4*meana11*ska11 &
                   - 3.*meana11**4
    ska11 = ska11 - 3*rmsa11*meana11 + 2*meana11**3

    rmsa11 = sqrt( rmsa11 - meana11 * meana11 )
    ska11 = ska11 / rmsa11**3
    fla11 = fla11 / rmsa11**4
    
    ! a22
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * ky(jj) * uy(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana22 = 0.d0
    rmsa22 = 0.d0
    ska22 = 0.d0
    fla22 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana22 = meana22 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa22 = rmsa22 + tmp1
      tmp1 = tmp1 * ignore_me
      ska22 = ska22 + tmp1
      tmp1 = tmp1 * ignore_me
      fla22 = fla22 + tmp1
    end do
    end do
    end do
    meana22 = meana22 * tmp2

    rmsa22 = rmsa22 * tmp2
    fla22 = fla22 * tmp2
    ska22 = ska22 * tmp2
    fla22 = fla22 + 6.*meana22**2*rmsa22 - 4*meana22*ska22 &
                   - 3.*meana22**4
    ska22 = ska22 - 3*rmsa22*meana22 + 2*meana22**3

    rmsa22 = sqrt( rmsa22 - meana22 * meana22 )
    ska22 = ska22 / rmsa22**3
    fla22 = fla22 / rmsa22**4
    
    
    ! a33
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kz(kk) * uz(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana33 = 0.d0
    rmsa33 = 0.d0
    ska33 = 0.d0
    fla33 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana33 = meana33 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa33 = rmsa33 + tmp1
      tmp1 = tmp1 * ignore_me
      ska33 = ska33 + tmp1
      tmp1 = tmp1 * ignore_me
      fla33 = fla33 + tmp1
    end do
    end do
    end do
    meana33 = meana33 * tmp2

    rmsa33 = rmsa33 * tmp2
    fla33 = fla33 * tmp2
    ska33 = ska33 * tmp2
    fla33 = fla33 + 6.*meana33**2*rmsa33 - 4*meana33*ska33 &
                   - 3.*meana33**4
    ska33 = ska33 - 3*rmsa33*meana33 + 2*meana33**3

    rmsa33 = sqrt( rmsa33 - meana33 * meana33 )
    ska33 = ska33 / rmsa33**3
    fla33 = fla33 / rmsa33**4

    ! Other components of gradient, skewness, flatness.

    ! a12
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * ky(jj) * ux(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana12 = 0.d0
    rmsa12 = 0.d0
    ska12 = 0.d0
    fla12 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana12 = meana12 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa12 = rmsa12 + tmp1
      tmp1 = tmp1 * ignore_me
      ska12 = ska12 + tmp1
      tmp1 = tmp1 * ignore_me
      fla12 = fla12 + tmp1
    end do
    end do
    end do
    meana12 = meana12 * tmp2

    rmsa12 = rmsa12 * tmp2
    fla12 = fla12 * tmp2
    ska12 = ska12 * tmp2
    fla12 = fla12 + 6.*meana12**2*rmsa12 - 4*meana12*ska12 &
                   - 3.*meana12**4
    ska12 = ska12 - 3*rmsa12*meana12 + 2*meana12**3

    rmsa12 = sqrt( rmsa12 - meana12 * meana12 )
    ska12 = ska12 / rmsa12**3
    fla12 = fla12 / rmsa12**4


    ! a13
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kz(kk) * ux(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana13 = 0.d0
    rmsa13 = 0.d0
    ska13 = 0.d0
    fla13 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana13 = meana13 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa13 = rmsa13 + tmp1
      tmp1 = tmp1 * ignore_me
      ska13 = ska13 + tmp1
      tmp1 = tmp1 * ignore_me
      fla13 = fla13 + tmp1
    end do
    end do
    end do
    meana13 = meana13 * tmp2

    rmsa13 = rmsa13 * tmp2
    fla13 = fla13 * tmp2
    ska13 = ska13 * tmp2
    fla13 = fla13 + 6.*meana13**2*rmsa13 - 4*meana13*ska13 &
                   - 3.*meana13**4
    ska13 = ska13 - 3*rmsa13*meana13 + 2*meana13**3

    rmsa13 = sqrt( rmsa13 - meana13 * meana13 )
    ska13 = ska13 / rmsa13**3
    fla13 = fla13 / rmsa13**4


    ! a21
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kx(ii) * uy(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)


    meana21 = 0.d0
    rmsa21 = 0.d0
    ska21 = 0.d0
    fla21 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana21 = meana21 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa21 = rmsa21 + tmp1
      tmp1 = tmp1 * ignore_me
      ska21 = ska21 + tmp1
      tmp1 = tmp1 * ignore_me
      fla21 = fla21 + tmp1
    end do
    end do
    end do
    meana21 = meana21 * tmp2

    rmsa21 = rmsa21 * tmp2
    fla21 = fla21 * tmp2
    ska21 = ska21 * tmp2
    fla21 = fla21 + 6.*meana21**2*rmsa21 - 4*meana21*ska21 &
                   - 3.*meana21**4
    ska21 = ska21 - 3*rmsa21*meana21 + 2*meana21**3

    rmsa21 = sqrt( rmsa21 - meana21 * meana21 )
    ska21 = ska21 / rmsa21**3
    fla21 = fla21 / rmsa21**4


    ! a23
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kz(kk) * uy(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana23 = 0.d0
    rmsa23 = 0.d0
    ska23 = 0.d0
    fla23 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana23 = meana23 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa23 = rmsa23 + tmp1
      tmp1 = tmp1 * ignore_me
      ska23 = ska23 + tmp1
      tmp1 = tmp1 * ignore_me
      fla23 = fla23 + tmp1
    end do
    end do
    end do
    meana23 = meana23 * tmp2

    rmsa23 = rmsa23 * tmp2
    fla23 = fla23 * tmp2
    ska23 = ska23 * tmp2
    fla23 = fla23 + 6.*meana23**2*rmsa23 - 4*meana23*ska23 &
                   - 3.*meana23**4
    ska23 = ska23 - 3*rmsa23*meana23 + 2*meana23**3

    rmsa23 = sqrt( rmsa23 - meana23 * meana23 )
    ska23 = ska23 / rmsa23**3
    fla23 = fla23 / rmsa23**4

    ! a31
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * kx(ii) * uz(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)

    meana31 = 0.d0
    rmsa31 = 0.d0
    ska31 = 0.d0
    fla31 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana31 = meana31 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa31 = rmsa31 + tmp1
      tmp1 = tmp1 * ignore_me
      ska31 = ska31 + tmp1
      tmp1 = tmp1 * ignore_me
      fla31 = fla31 + tmp1
    end do
    end do
    end do
    meana31 = meana31 * tmp2

    rmsa31 = rmsa31 * tmp2
    fla31 = fla31 * tmp2
    ska31 = ska31 * tmp2
    fla31 = fla31 + 6.*meana31**2*rmsa31 - 4*meana31*ska31 &
                   - 3.*meana31**4
    ska31 = ska31 - 3*rmsa31*meana31 + 2*meana31**3

    rmsa31 = sqrt( rmsa31 - meana31 * meana31 )
    ska31 = ska31 / rmsa31**3
    fla31 = fla31 / rmsa31**4


    ! a32
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
    aij(ii,jj,kk) = eye * ky(jj) * uz(ii,jj,kk)
    end do
    end do
    end do
    !call rfftwnd_f77_one_complex_to_real(c2r3d,aij,ignore_me)
    call dfftw_execute(aijc2r)


    meana32 = 0.d0
    rmsa32 = 0.d0
    ska32 = 0.d0
    fla32 = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( aij( (ii+1)/2,jj,kk ) )
      else
        ignore_me = aimag( aij( ii/2,jj,kk ) )
      end if
      meana32 = meana32 + ignore_me
      tmp1 = ignore_me * ignore_me
      rmsa32 = rmsa32 + tmp1
      tmp1 = tmp1 * ignore_me
      ska32 = ska32 + tmp1
      tmp1 = tmp1 * ignore_me
      fla32 = fla32 + tmp1
    end do
    end do
    end do
    meana32 = meana32 * tmp2

    rmsa32 = rmsa32 * tmp2
    fla32 = fla32 * tmp2
    ska32 = ska32 * tmp2
    fla32 = fla32 + 6.*meana32**2*rmsa32 - 4*meana32*ska32 &
                   - 3.*meana32**4
    ska32 = ska32 - 3*rmsa32*meana32 + 2*meana32**3

    rmsa32 = sqrt( rmsa32 - meana32 * meana32 )
    ska32 = ska32 / rmsa32**3
    fla32 = fla32 / rmsa32**4



    !write(*,*) 'before vorticity '

    ! Vorticity statistics

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a12 = eye * ky(jj) * ux(ii,jj,kk)
      a13 = eye * kz(kk) * ux(ii,jj,kk)
      a21 = eye * kx(ii) * uy(ii,jj,kk)
      a23 = eye * kz(kk) * uy(ii,jj,kk)
      a31 = eye * kx(ii) * uz(ii,jj,kk)
      a32 = eye * ky(jj) * uz(ii,jj,kk)

      ! ux uy uz are now vorticity components
      ux(ii,jj,kk) = a32 - a23
      uy(ii,jj,kk) = a13 - a31
      uz(ii,jj,kk) = a21 - a12
    end do
    end do
    end do

    !write(*,*) 'after vorticity '

    ! Fourier Transform
    !call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    !call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    !call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
    call dfftw_execute(uxc2r)
    call dfftw_execute(uyc2r)
    call dfftw_execute(uzc2r)
    !write(*,*) 'after fftw of vorticity '

    ! skewness and flatness of vorticity components
    meanwx = 0.d0; meanwy = 0.d0; meanwz = 0.d0
    rmswx = 0.d0; rmswy = 0.d0; rmswz = 0.d0
    skwx = 0.d0; skwy = 0.d0; skwz = 0.d0
    flwx = 0.d0; flwy = 0.d0; flwz = 0.d0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then 
        ignore_me = real( ux( (ii+1)/2,jj,kk ) )
        tmp4 = real( uy( (ii+1)/2,jj,kk ) )
        tmp3 = real( uz( (ii+1)/2, jj, kk ) )
      else
        ignore_me = aimag( ux( ii/2,jj,kk ) )
        tmp4 = aimag( uy( (ii+1)/2,jj,kk ) )
        tmp3 = aimag( uz( (ii+1)/2, jj, kk ) )
      end if

      meanwx = meanwx + ignore_me
      tmp1 = ignore_me * ignore_me
      rmswx = rmswx + tmp1
      tmp1 = tmp1 * ignore_me
      skwx = skwx + tmp1
      tmp1 = tmp1 * ignore_me
      flwx = flwx + tmp1

      meanwy = meanwy + tmp4
      tmp1 = tmp4 * tmp4
      rmswy = rmswy + tmp1
      tmp1 = tmp1 * tmp4
      skwy = skwy + tmp1
      tmp1 = tmp1 * tmp4
      flwy = flwy + tmp1

      meanwz = meanwz + tmp3
      tmp1 = tmp3 * tmp3
      rmswz = rmswz + tmp1
      tmp1 = tmp1 * tmp3
      skwz = skwz + tmp1
      tmp1 = tmp1 * tmp3
      flwz = flwz + tmp1

    end do
    end do
    end do

    !write(*,*) 'after vorticity statistics'
    meanwx = meanwx * tmp2

    rmswx = rmswx * tmp2
    flwx = flwx * tmp2
    skwx = skwx * tmp2
    flwx = flwx + 6.*meanwx**2*rmswx - 4*meanwx*skwx &
                   - 3.*meanwx**4
    skwx = skwx - 3*rmswx*meanwx + 2*meanwx**3

    rmswx = sqrt( rmswx - meanwx * meanwx )
    skwx = skwx / rmswx**3
    flwx = flwx / rmswx**4

    meanwy = meanwy * tmp2

    rmswy = rmswy * tmp2
    flwy = flwy * tmp2
    skwy = skwy * tmp2
    flwy = flwy + 6.*meanwy**2*rmswy - 4*meanwy*skwy &
                   - 3.*meanwy**4
    skwy = skwy - 3*rmswy*meanwy + 2*meanwy**3

    rmswy = sqrt( rmswy - meanwy * meanwy )
    skwy = skwy / rmswy**3
    flwy = flwy / rmswy**4

    meanwz = meanwz * tmp2

    rmswz = rmswz * tmp2
    flwz = flwz * tmp2
    skwz = skwz * tmp2
    flwz = flwz + 6.*meanwz**2*rmswz - 4*meanwz*skwz &
                   - 3.*meanwz**4
    skwz = skwz - 3*rmswz*meanwz + 2*meanwz**3

    rmswz = sqrt( rmswz - meanwz * meanwz )
    skwz = skwz / rmswz**3
    flwz = flwz / rmswz**4

    !str = str(1: len_trim(str) - 4) ! removing .dat part
    !read(str,'(I5)') ll

    write(14, '(I5, 20E12.3)') ll, ek, eps
    write(19, '(I5, 20E12.3)') ll, meana11,meana22,meana33,meana12,meana13,meana21,meana23,meana31,meana32
    write(15, '(I5, 20E12.3)') ll, rmsa11,rmsa22,rmsa33,rmsa12,rmsa13,rmsa21,rmsa23,rmsa31,rmsa32
    write(16, '(I5, 20E12.3)') ll, ska11,ska22,ska33,ska12,ska13,ska21,ska23,ska31,ska32
    write(17, '(I5, 20E12.3)') ll, fla11,fla22,fla33,fla12,fla13,fla21,fla23,fla31,fla32
    write(18, '(I5, 20E12.3)') ll, meanwx,meanwy,meanwz,rmswx,rmswy,rmswz,skwx,skwy,skwz,flwx,flwy,flwz

    ll = ll + 1

  end do


  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)


  deallocate(ux,uy,uz,aij,kx,ky,kz,k2)

  !call destroyplan3d
  write(*,*) 'means.x done'

end program means 
