program sgshdmoments
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  integer :: lx1,lx,ly,lz,nx,ny,nz,ii,jj,kk,ll,nfile,ndel,ifilter

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, wx, wy, wz, fux, fuy, fuz
  complex(sp), allocatable, dimension(:,:,:) :: rij, tij, helidiss
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, g, k2

  real(dp) :: m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
  real(sp) :: ignore_me, const, delta_c, delta2, kc2

  character(80) :: fnm,str,str1,fpath,strfltr

  write(*,*) 
  write(*,'(''>>> Moments of SGS helicity dissipation <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sgshdmoments.x nx filelist ndel ifilter'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        ifilter: = 1 for Gaussian; = 0 for cutoff'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! filter type
  call getarg(4,strfltr)
  read(strfltr, '(I20)') ifilter

  ! some parameters 
  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  ! filter scale
  delta_c = ndel*2*pi/nx
  delta2 = delta_c * delta_c
  kc2 = (nx/2/ndel)**2

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3de'

  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
  allocate( g(lx1,ly,lz), k2(lx1,ly,lz), helidiss(lx1,ly,lz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( fux(lx1,ly,lz), fuy(lx1,ly,lz), fuz(lx1,ly,lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( rij(lx1,ly,lz), tij(lx1,ly,lz) )
  write(*,*) 'allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if ( ifilter .eq. 1 ) then

      ! Gaussian filter
      g=exp(-k2*delta_c**2/24.)
      strfltr = 'gau'

  else if ( ifilter .eq. 0 ) then

      ! Cutoff 
      where ( k2 .ge. kc2 )
          g = 0
      elsewhere
          g = 1
      endwhere
      strfltr = 'cut'

  else
      stop 'Wrong filter type. Stopping '
  end if


  nfile = 0
  m1 = 0.0_dp; m2 = 0.0_dp; m3 = 0.0_dp; m4 = 0.0_dp; m5 = 0.0_dp
  m6 = 0.0_dp; m7 = 0.0_dp; m8 = 0.0_dp; m9 = 0.0_dp; m10 = 0.0_dp

  open( 20, file = fnm( 1:len_trim(fnm) )//'.list' )
    do while ( .not. eof(20) )
      read(20,*) str1
      write(*,*) str1( 1:len_trim(str1) )

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)ux
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'

      wx = eye * ( ky * uz - kz * uy ) * g
      wy = eye * ( kz * ux - kx * uz ) * g
      wz = eye * ( kx * uy - ky * ux ) * g
      fux = ux * g; fuy = uy * g; fuz = uz * g

      call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)

      call rfftwnd_f77_one_complex_to_real(c2r3d,fux,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,fuy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,fuz,ignore_me)

      rij = eye * kx * wx
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(ux) * real(ux), aimag(ux) * aimag(ux) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fux) * real(fux), aimag(fux) * aimag(fux) )

      helidiss = cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      rij = eye * ky * wy
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(uy) * real(uy), aimag(uy) * aimag(uy) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fuy) * real(fuy), aimag(fuy) * aimag(fuy) )

      helidiss = helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      rij = eye * kz * wz
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(uz) * real(uz), aimag(uz) * aimag(uz) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fuz) * real(fuz), aimag(fuz) * aimag(fuz) )

      helidiss = helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      rij = .5 * eye * ( kx * wy + ky * wx )
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(ux) * real(uy), aimag(ux) * aimag(uy) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fux) * real(fuy), aimag(fux) * aimag(fuy) )

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      rij = .5 * eye * ( kx * wz + kz * wx )
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(ux) * real(uz), aimag(ux) * aimag(uz) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fux) * real(fuz), aimag(fux) * aimag(fuz) )

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      rij = .5 * eye * ( ky * wz + kz * wy )
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      tij = cmplx( real(uy) * real(uz), aimag(uy) * aimag(uz) ) * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - cmplx( real(fuy) * real(fuz), aimag(fuy) * aimag(fuz) )

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )


      helidiss = -2. * helidiss
      
      m1  = m1  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(1./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(1./3))
      m2  = m2  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(2./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(2./3))
      m3  = m3  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) )         ) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) )         )
      m4  = m4  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(4./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(4./3))
      m5  = m5  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(5./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(5./3))
      m6  = m6  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **2     ) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **2     )
      m7  = m7  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(7./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(7./3))
      m8  = m8  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(8./3)) + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(8./3))
      m9  = m9  + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **3)      + &
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **3)
      m10 = m10 + sum( ( abs( real(delta2 * helidiss(1:lx,:,:))) ) **(10./3)) + & 
                  sum( ( abs(aimag(delta2 * helidiss(1:lx,:,:))) ) **(10./3))

      nfile = nfile + 1

    end do
  close(20)

  ignore_me = const / nfile
  m1 = m1 * ignore_me
  m2 = m2 * ignore_me
  m3 = m3 * ignore_me
  m4 = m4 * ignore_me
  m5 = m5 * ignore_me
  m6 = m6 * ignore_me
  m7 = m7 * ignore_me
  m8 = m8 * ignore_me
  m9 = m9 * ignore_me
  m10 = m10 * ignore_me

  strfltr = 'sgshdmoments-'//strfltr(1:len_trim(strfltr))//'-' 
  strfltr = strfltr(1:len_trim(strfltr))//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file = strfltr(1:len_trim(strfltr)) )
    write(15,'(20E15.4)') delta_c, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
  close(15)

  deallocate( kx, ky, kz, k2, g, ux, uy, uz, wx, wy, wz, fux, fuy, fuz, rij, tij, helidiss )

  call destroyplan3d

  write(*,*) 'sgshdmoments.x finished'

end program sgshdmoments
