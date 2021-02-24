program curvaturecndxx
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(sp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  complex(sp), allocatable, dimension(:,:,:) :: oh11,oh12,oh13,oh21,oh22,oh23,oh31,oh32,oh33
  real(sp),    allocatable, dimension(:,:,:) :: g, curv
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer,  parameter :: npnt = 200
  real(sp), parameter :: boundq = 10.
  real(sp), parameter :: binwq = boundq/npnt
  real(dp), dimension(npnt) :: pdfxx, cndoh

  real(dp) :: const, rmsrr
  real(sp) :: delta_c, ignore_me, rmsrr0, tmp
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> Curvature of vortex lines conditioned on q,r of gij <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./curvature-cnd-xx.x nx filelist ndel'
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
  allocate( oh11(lx1,ly,lz), oh12(lx1,ly,lz), oh13(lx1,ly,lz), oh21(lx1,ly,lz), oh22(lx1,ly,lz) )
  allocate( oh23(lx1,ly,lz), oh31(lx1,ly,lz), oh32(lx1,ly,lz), oh33(lx1,ly,lz) )
  allocate( curv(nx,ny,nz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pdfxx = 0._dp
  cndoh = 0._dp
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
      g12(ii,jj,kk) = cmplx( real( g11(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g11(ii,jj,kk) ) / (delta_c + mytiny_sp) )
      g21(ii,jj,kk) = cmplx( real( g22(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g22(ii,jj,kk) ) / (delta_c + mytiny_sp) )
      g31(ii,jj,kk) = cmplx( real( g33(ii,jj,kk) ) / (ignore_me + mytiny_sp), &
                            aimag( g33(ii,jj,kk) ) / (delta_c + mytiny_sp) )

      ! backup the direction in real space
      g13(ii,jj,kk) = g12(ii,jj,kk)
      g23(ii,jj,kk) = g21(ii,jj,kk)
      g32(ii,jj,kk) = g31(ii,jj,kk)
    end do
    end do
    end do

    g12 = g12 * const; g21 = g21 * const; g31 = g31 * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,g12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,g21,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,g31,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      oh11(ii,jj,kk) = eye * kx(ii) * g12(ii,jj,kk)
      oh12(ii,jj,kk) = eye * ky(jj) * g12(ii,jj,kk)
      oh13(ii,jj,kk) = eye * kz(kk) * g12(ii,jj,kk)

      oh21(ii,jj,kk) = eye * kx(ii) * g21(ii,jj,kk)
      oh22(ii,jj,kk) = eye * ky(jj) * g21(ii,jj,kk)
      oh23(ii,jj,kk) = eye * kz(kk) * g21(ii,jj,kk)

      oh31(ii,jj,kk) = eye * kx(ii) * g31(ii,jj,kk)
      oh32(ii,jj,kk) = eye * ky(jj) * g31(ii,jj,kk)
      oh33(ii,jj,kk) = eye * kz(kk) * g31(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,oh11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oh33,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx
      delta_c = real( oh11(ii,jj,kk) ) * real( g13(ii,jj,kk) ) &
              + real( oh12(ii,jj,kk) ) * real( g23(ii,jj,kk) ) &
              + real( oh13(ii,jj,kk) ) * real( g32(ii,jj,kk) )
      ignore_me = real( oh21(ii,jj,kk) ) * real( g13(ii,jj,kk) ) &
                + real( oh22(ii,jj,kk) ) * real( g23(ii,jj,kk) ) &
                + real( oh23(ii,jj,kk) ) * real( g32(ii,jj,kk) )
      tmp = real( oh31(ii,jj,kk) ) * real( g13(ii,jj,kk) ) &
          + real( oh32(ii,jj,kk) ) * real( g23(ii,jj,kk) ) &
          + real( oh33(ii,jj,kk) ) * real( g32(ii,jj,kk) )

      curv(2*ii-1,jj,kk) = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)

      delta_c = aimag( oh11(ii,jj,kk) ) * aimag( g13(ii,jj,kk) ) &
              + aimag( oh12(ii,jj,kk) ) * aimag( g23(ii,jj,kk) ) &
              + aimag( oh13(ii,jj,kk) ) * aimag( g32(ii,jj,kk) )
      ignore_me = aimag( oh21(ii,jj,kk) ) * aimag( g13(ii,jj,kk) ) &
                + aimag( oh22(ii,jj,kk) ) * aimag( g23(ii,jj,kk) ) &
                + aimag( oh23(ii,jj,kk) ) * aimag( g32(ii,jj,kk) )
      tmp = aimag( oh31(ii,jj,kk) ) * aimag( g13(ii,jj,kk) ) &
          + aimag( oh32(ii,jj,kk) ) * aimag( g23(ii,jj,kk) ) &
          + aimag( oh33(ii,jj,kk) ) * aimag( g32(ii,jj,kk) )

      curv(2*ii,jj,kk) = sqrt( delta_c ** 2 + ignore_me ** 2 + tmp ** 2)

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
      ! xi
      g12(ii,jj,kk) = eye * (ky(jj) * g33(ii,jj,kk) - kz(kk) * g22(ii,jj,kk) )
      g13(ii,jj,kk) = eye * (kz(kk) * g11(ii,jj,kk) - kx(ii) * g33(ii,jj,kk) )
      g21(ii,jj,kk) = eye * (kx(ii) * g22(ii,jj,kk) - ky(jj) * g11(ii,jj,kk) )

    end do
    end do
    end do
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g21,ignore_me)


    if ( nfile .eq. 0 ) then
      rmsrr0 = sum( g12(1:lx,:,:) * conjg( g12(1:lx,:,:) ) ) &
            + sum( g13(1:lx,:,:) * conjg( g13(1:lx,:,:) ) ) &
            + sum( g21(1:lx,:,:) * conjg( g21(1:lx,:,:) ) ) 
      rmsrr0 = sqrt( rmsrr0 * const )         

      write(*,*) 'Estimated rms of rr: ', rmsrr0
    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii + 1 )/2

        delta_c = real( g12(ll,jj,kk) )
        ignore_me = real( g13(ll,jj,kk) )
        tmp = real( g21(ll,jj,kk) )

      else
        ll = ii / 2

        delta_c = aimag( g12(ll,jj,kk) )
        ignore_me = aimag( g13(ll,jj,kk) )
        tmp = aimag( g21(ll,jj,kk) )

      end if

      tmp = sqrt( tmp * tmp + delta_c * delta_c + ignore_me * ignore_me )
      rmsrr = rmsrr + tmp * tmp
      tmp = tmp / rmsrr0

      ll = floor( (tmp) / binwq ) + 1

      if ( ll .ge. 1 .and. ll .le. npnt ) then
          pdfxx(ll) = pdfxx(ll) + 1
          cndoh(ll) = cndoh(ll) + curv(ii,jj,kk)
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

  cndoh = cndoh / (pdfxx + mytiny_dp)

  pdfxx = pdfxx / nfile * const
  write(*,*) 'Check normalizationof pdfxx:', sum(pdfxx)

  pdfxx = pdfxx / binwq 


  ! Ouput format for gnuplot
  open(15,file='curvature-cnd-xx'//str(1:len_trim(str)))
    write(15,*) '# Zone T= "curvature | xx", i=', npnt, ', F=point'
    do jj=1,npnt
      write(15,'(15E15.5)') ( (jj-.5)*binwq)*rmsrr0**2/rmsrr**2, &
                pdfxx(jj)*rmsrr**2/rmsrr0**2, cndoh(jj)
    end do
  close(15)

  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(oh11,oh12,oh13,oh21,oh22,oh23,oh31,oh32,oh33)
  deallocate(curv)

  call destroyplan3d

  write(*,*) 'Finished'
end program curvaturecndxx
