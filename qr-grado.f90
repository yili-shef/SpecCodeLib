program qrgrado
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(sp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  real(sp), allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  real(sp), dimension(3,3) :: gg

  integer, parameter :: npnt = 160
  real(sp), parameter :: boundq = 8., boundr = 8.
  real(sp), parameter :: binwq = 2*boundq/npnt, binwr = 2*boundr/npnt
  real(dp), dimension(npnt, npnt) :: pqr

  real(dp) :: const
  real :: delta_c, ignore_me, rmsrr, rmsrr0
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> Joint PDF of (r,q) for gradient of vorticity <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./qr-grado.x nx filelist ndel'
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
  allocate( g(lx1,ly,lz) )
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz), g21(lx1,ly,lz), g22(lx1,ly,lz) )
  allocate( g23(lx1,ly,lz), g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pqr = 0._dp
  rmsrr = 0._dp
  nfile = 0
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

      ! Gradient of vorticity 
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

    if ( nfile .eq. 0 ) then
      rmsrr0 = sum( g11(1:lx,:,:) * conjg( g11(1:lx,:,:) ) ) &
            + sum( g12(1:lx,:,:) * conjg( g12(1:lx,:,:) ) ) &
            + sum( g13(1:lx,:,:) * conjg( g13(1:lx,:,:) ) ) &
            + sum( g21(1:lx,:,:) * conjg( g21(1:lx,:,:) ) ) &
            + sum( g22(1:lx,:,:) * conjg( g22(1:lx,:,:) ) ) &
            + sum( g23(1:lx,:,:) * conjg( g23(1:lx,:,:) ) ) &
            + sum( g31(1:lx,:,:) * conjg( g31(1:lx,:,:) ) ) &
            + sum( g32(1:lx,:,:) * conjg( g32(1:lx,:,:) ) ) &
            + sum( g33(1:lx,:,:) * conjg( g33(1:lx,:,:) ) ) 
      rmsrr0 = sqrt( rmsrr0 * const )         

      write(*,*) 'Estimated rms of rr: ', rmsrr0
    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii + 1 )/2
        gg(1,1) = real( g11(ll,jj,kk) )
        gg(1,2) = real( g12(ll,jj,kk) )
        gg(1,3) = real( g13(ll,jj,kk) )
        gg(2,1) = real( g21(ll,jj,kk) )
        gg(2,2) = real( g22(ll,jj,kk) )
        gg(2,3) = real( g23(ll,jj,kk) )
        gg(3,1) = real( g31(ll,jj,kk) )
        gg(3,2) = real( g32(ll,jj,kk) )
        gg(3,3) = real( g33(ll,jj,kk) )
      else
        ll = ii / 2
        gg(1,1) = aimag( g11(ll,jj,kk) )
        gg(1,2) = aimag( g12(ll,jj,kk) )
        gg(1,3) = aimag( g13(ll,jj,kk) )
        gg(2,1) = aimag( g21(ll,jj,kk) )
        gg(2,2) = aimag( g22(ll,jj,kk) )
        gg(2,3) = aimag( g23(ll,jj,kk) )
        gg(3,1) = aimag( g31(ll,jj,kk) )
        gg(3,2) = aimag( g32(ll,jj,kk) )
        gg(3,3) = aimag( g33(ll,jj,kk) )
      end if

      rmsrr = rmsrr + sum( gg * gg )

      ignore_me = -.5 * sum(  gg * transpose(gg) ) / rmsrr0**2
      delta_c = -sum ( matmul(gg,gg) * transpose(gg) ) / 3. /rmsrr0**3

      ll = floor( (ignore_me + boundq) / binwq ) + 1
      mm = floor( (delta_c + boundr) / binwr ) + 1

      if ( ll .ge. 1 .and. ll .le. npnt .and. mm .ge. 1 .and. mm .le. npnt) then
        pqr(mm,ll) = pqr(mm,ll) + 1
      end if

    end do
    end do
    end do
    nfile = nfile + 1
  end do
  close(20)

  rmsrr = rmsrr / nfile * const
  rmsrr = sqrt(rmsrr)
  write(*,*) 'rmsrr is: ', rmsrr

  pqr = pqr / nfile * const
  write(*,*) 'Check normalizationof pqr:', sum(pqr)

  pqr = pqr / binwq / binwr

  open(15,file='jp-grado-rq'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(r,q)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (-boundr + (ii-.5)*binwr)*rmsrr0**3/rmsrr**3, & 
            (-boundq + (jj-.5)*binwq)*rmsrr0**2/rmsrr**2, pqr(ii,jj)
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)

  call destroyplan3d

  write(*,*) 'Finished'
end program qrgrado
