program qrgrado
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex(sp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  complex(sp), allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  real(sp), dimension(3,3) :: gg

  integer, parameter :: npnt = 80
  real(sp), parameter :: boundq = 2., boundr = 2.
  real(sp), parameter :: binwq = 2*boundq/npnt, binwr = 2*boundr/npnt
  real(dp), dimension(npnt, npnt) :: pqrp, pqrm

  integer(8) :: nump, numm
  real(dp) :: const
  real :: delta_c, ignore_me, rmsrr, rmsrr0, pih
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>>>>> Joint PDF of (r,q) for gradient of vorticity <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./qr-grado-cndpih.x nx filelist ndel'
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
  allocate( t11(lx1,ly,lz), t12(lx1,ly,lz), t13(lx1,ly,lz), t22(lx1,ly,lz) )
  allocate( t23(lx1,ly,lz), t33(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pqrp = 0._dp
  pqrm = 0._dp
  rmsrr = 0._dp
  nfile = 0
  nump = 0
  numm = 0
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

    call rsgstauij(g31,g31,t11,g,nx,ny,nz)
    call rsgstauij(g31,g32,t12,g,nx,ny,nz)
    call rsgstauij(g31,g33,t13,g,nx,ny,nz)
    call rsgstauij(g32,g32,t22,g,nx,ny,nz)
    call rsgstauij(g32,g33,t23,g,nx,ny,nz)
    call rsgstauij(g33,g33,t33,g,nx,ny,nz)

    g11 = - (t11+t22+t33)/3.
    t11 = g11 + t11
    t22 = g11 + t22
    t33 = g11 + t33
    
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

        pih = gg(1,1) * real(t11(ll,jj,kk)) + gg(1,2) * real( t12(ll,jj,kk) ) + gg(1,3) * real(t13(ll,jj,kk)) &
            + gg(2,1) * real(t12(ll,jj,kk)) + gg(2,2) * real( t22(ll,jj,kk) ) + gg(2,3) * real(t23(ll,jj,kk)) &
            + gg(3,1) * real(t13(ll,jj,kk)) + gg(3,2) * real( t23(ll,jj,kk) ) + gg(3,3) * real(t33(ll,jj,kk)) 
        pih = -2 * pih
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

        pih = gg(1,1) * aimag(t11(ll,jj,kk)) + gg(1,2) * aimag( t12(ll,jj,kk) ) + gg(1,3) * aimag(t13(ll,jj,kk)) &
            + gg(2,1) * aimag(t12(ll,jj,kk)) + gg(2,2) * aimag( t22(ll,jj,kk) ) + gg(2,3) * aimag(t23(ll,jj,kk)) &
            + gg(3,1) * aimag(t13(ll,jj,kk)) + gg(3,2) * aimag( t23(ll,jj,kk) ) + gg(3,3) * aimag(t33(ll,jj,kk)) 
        pih = -2 * pih
      end if

      rmsrr = rmsrr + sum( gg * gg )

      ignore_me = -.5 * sum(  gg * transpose(gg) ) / rmsrr0**2
      delta_c = -sum ( matmul(gg,gg) * transpose(gg) ) / 3. /rmsrr0**3

      ll = floor( (ignore_me + boundq) / binwq ) + 1
      mm = floor( (delta_c + boundr) / binwr ) + 1

      if ( ll .ge. 1 .and. ll .le. npnt .and. mm .ge. 1 .and. mm .le. npnt) then

        if ( pih .ge. 0 ) then
          pqrp(mm,ll) = pqrp(mm,ll) + 1
          nump = nump + 1
        else
          pqrm(mm,ll) = pqrm(mm,ll) + 1
          numm = numm + 1
        end if

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

  pqrp = pqrp / nump
  write(*,*) 'Check normalizationof pqrp:', sum(pqrp)

  pqrp = pqrp / binwq / binwr

  pqrm = pqrm / numm
  write(*,*) 'Check normalizationof pqrm:', sum(pqrm)

  pqrm = pqrm / binwq / binwr

  open(15,file='qr-grado-cndpih'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(r,q)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (-boundr + (ii-.5)*binwr)*rmsrr0**3/rmsrr**3, & 
            (-boundq + (jj-.5)*binwq)*rmsrr0**2/rmsrr**2, pqrp(ii,jj)*rmsrr**5/rmsrr0**5, &
            pqrm(ii,jj)*rmsrr**5/rmsrr0**5

    end do
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(t11,t12,t13,t22,t23,t33)

  call destroyplan3d

  write(*,*) 'Finished'
end program qrgrado
