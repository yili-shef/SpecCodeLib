program qrrij
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile, nthresh

  complex(sp), allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33,oo
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  real(sp), allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  real(sp), dimension(3,3) :: gg

  integer, parameter :: npnt = 400
  real(sp), parameter :: boundq = 8., boundr = 4.
  real(sp), parameter :: binwq = boundq/npnt, binwr = 2*boundr/npnt
  real(dp), dimension(npnt, npnt) :: pqr

  real(dp) :: const, rmsrr, rmso
  real(sp) :: delta_c, ignore_me, rmsrr0, rmso0, tmpoo

  integer(8) :: numpnt 

  character(80) :: fnm,str,str1,fpath, disslist

  write(*,*) 
  write(*,'(''>>>>>> Joint PDF of (r,q) for rij conditioned on vorticity <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./qr-rij-cndoo.x nx dnsfilelist ndel nthresh'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        dnsfilelist: dns data file list'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*) '        nthresh: threshold value = nthresh * mdissrate '
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

  ! nthresh
  call getarg(4,disslist)
  read(disslist, '(I20)') nthresh
  disslist = adjustl(disslist)

  str='-'//str(1:len_trim(str))//'dx-cndoo'//disslist(1:len_trim(disslist))//'-'//fnm(1:len_trim(fnm))//'.dat'

  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx
  const = 1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( g(lx1,ly,lz) , oo(lx1,ly,lz) )
  allocate( r11(lx1,ly,lz), r12(lx1,ly,lz), r13(lx1,ly,lz), r22(lx1,ly,lz) )
  allocate( r23(lx1,ly,lz), r33(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pqr = 0._dp
  rmso = 0._dp
  rmsrr = 0._dp
  numpnt = 0
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10) r12
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10) r13
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10) r23
    close(10)

    write(*,*) 'after reading data files'
    
    r12 = r12 * g
    r13 = r13 * g
    r23 = r23 * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1

      ! vorticity 
      wx(ii,jj,kk) = eye * ( ky(jj) * r23(ii,jj,kk) - kz(kk) * r13(ii,jj,kk) )
      wy(ii,jj,kk) = eye * ( kz(kk) * r12(ii,jj,kk) - kx(ii) * r23(ii,jj,kk) )
      wz(ii,jj,kk) = eye * ( kx(ii) * r13(ii,jj,kk) - ky(jj) * r12(ii,jj,kk) )

      ! rij
      r12(ii,jj,kk) = .5 * eye * ( ky(jj) * wx(ii,jj,kk) + kx(ii) * wy(ii,jj,kk) )
      r13(ii,jj,kk) = .5 * eye * ( kx(ii) * wz(ii,jj,kk) + kz(kk) * wx(ii,jj,kk) )
      r23(ii,jj,kk) = .5 * eye * ( ky(jj) * wz(ii,jj,kk) + kz(kk) * wy(ii,jj,kk) )
      r11(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
      r22(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
      r33(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)

    end do
    end do
    end do
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    oo = cmplx( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz), &
                aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) )

    if ( nfile .eq. 0 ) then
      rmso0 = sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )
      rmso0 = sqrt( rmso0 / (nx*ny*nz) )
 
      write(*,*) 'Estimated rms of omega: ', rmso0
    end if
    rmso = rmso + sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )

    call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)

    if ( nfile .eq. 0 ) then
      rmsrr0 = sum( r11(1:lx,:,:) * conjg( r11(1:lx,:,:) ) ) &
             + sum( r12(1:lx,:,:) * conjg( r12(1:lx,:,:) ) ) &
             + sum( r13(1:lx,:,:) * conjg( r13(1:lx,:,:) ) ) &
             + sum( r12(1:lx,:,:) * conjg( r12(1:lx,:,:) ) ) &
             + sum( r22(1:lx,:,:) * conjg( r22(1:lx,:,:) ) ) &
             + sum( r23(1:lx,:,:) * conjg( r23(1:lx,:,:) ) ) &
             + sum( r13(1:lx,:,:) * conjg( r13(1:lx,:,:) ) ) &
             + sum( r23(1:lx,:,:) * conjg( r23(1:lx,:,:) ) ) &
             + sum( r33(1:lx,:,:) * conjg( r33(1:lx,:,:) ) ) 
      rmsrr0 = sqrt( rmsrr0 * const )         

      write(*,*) 'Estimated rms of rr: ', rmsrr0

    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii + 1 )/2
        tmpoo = sqrt( real( oo(ll,jj,kk) ) )

        gg(1,1) = real( r11(ll,jj,kk) )
        gg(1,2) = real( r12(ll,jj,kk) )
        gg(1,3) = real( r13(ll,jj,kk) )
        gg(2,2) = real( r22(ll,jj,kk) )
        gg(2,3) = real( r23(ll,jj,kk) )
        gg(3,3) = real( r33(ll,jj,kk) )
      else
        ll = ii / 2
        tmpoo = sqrt( aimag( oo(ll,jj,kk) ) )

        gg(1,1) = aimag( r11(ll,jj,kk) )
        gg(1,2) = aimag( r12(ll,jj,kk) )
        gg(1,3) = aimag( r13(ll,jj,kk) )
        gg(2,2) = aimag( r22(ll,jj,kk) )
        gg(2,3) = aimag( r23(ll,jj,kk) )
        gg(3,3) = aimag( r33(ll,jj,kk) )
      end if
      gg(2,1) = gg(1,2)
      gg(3,1) = gg(1,3)
      gg(3,2) = gg(2,3)

      rmsrr = rmsrr + sum( gg * gg )

      if ( tmpoo .ge. nthresh * rmso0 ) then

        numpnt = numpnt + 1

        ignore_me = -.5 * sum(  gg * transpose(gg) ) / rmsrr0**2
        delta_c = -sum ( matmul(gg,gg) * transpose(gg) ) / 3. /rmsrr0**3
     
        ll = floor( (ignore_me + boundq) / binwq ) + 1
        mm = floor( (delta_c + boundr) / binwr ) + 1
     
        if ( ll .ge. 1 .and. ll .le. npnt .and. mm .ge. 1 .and. mm .le. npnt) then
          pqr(mm,ll) = pqr(mm,ll) + 1
        end if
     
      end if

    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(20)
  close(21)

  rmsrr = rmsrr / nfile * const
  rmsrr = sqrt(rmsrr)
  write(*,*) 'rmsrr is: ', rmsrr

  rmso = rmso / nfile * const
  rmso = sqrt(rmso)
  write(*,*) 'rmso is: ', rmso

  pqr = pqr / numpnt 
  write(*,*) 'Check normalizationof pqr:', sum(pqr)

  write(*,*) 'Number of points: ', numpnt
  write(*,*) 'percentage of points: ', real(numpnt) / nfile * const

  pqr = pqr / binwq / binwr

  open(15,file='jp-rij-rq'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(r,q)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (-boundr + (ii-.5)*binwr)*rmsrr0**3/rmsrr**3, & 
            (-boundq + (jj-.5)*binwq)*rmsrr0**2/rmsrr**2, pqr(ii,jj)
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,g,r11,r12,r13,r22,r23,r33,oo,wx,wy,wz)

  call destroyplan3d

  write(*,*) 'Finished'
end program qrrij
