program vel2p
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  complex(sp), allocatable, dimension(:,:,:) :: p
  real(sp),    allocatable, dimension(:,:,:) :: k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(sp) :: ignore_me, const
  integer  :: ioinfo
  character(90) :: str, strndel, fnm

  write(*,*) 
  write(*,'(''>>> Calculate pressure fields (in Fourier space) from velocity fields <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./vel2p.x nx dnsfilelist'
          write(*,*) '        nx: resolution'
          write(*,*) '        dnsfilelist: dns data file list'
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

  fnm = fnm(1:len_trim(fnm))//'.list'

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1. / (nx * ny * nz)

  allocate( a11(lx1,ly,lz), a21(lx1,ly,lz), a31(lx1,ly,lz) )
  allocate( a12(lx1,ly,lz), a22(lx1,ly,lz), a32(lx1,ly,lz) )
  allocate( a13(lx1,ly,lz), a23(lx1,ly,lz), a33(lx1,ly,lz) )
  allocate( p(lx1,ly,lz), kx(lx1), ky(ly), kz(lz), k2(lx1,ly,lz) )
  write(*,*) 'after allocation'

  call fftwplan3de(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(25, file = fnm(1:len_trim(fnm)) )

    do while ( .not. eof(25) )
      read(25,*) str
      write(*,*) 'files ', str( 1:len_trim(str) )

      open(15, file = './out/p'//str( 1:len_trim(str) ), status = 'new', form = 'unformatted', &
           iostat = ioinfo, err = 99)

      open(10,file='./out/ux'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a13
      close(10)
      open(10,file='./out/uy'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a23
      close(10)
      open(10,file='./out/uz'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a33
      close(10)

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        a11(ii,jj,kk) = eye * a13(ii,jj,kk) * kx(ii)
        a12(ii,jj,kk) = eye * a13(ii,jj,kk) * ky(jj)
        a13(ii,jj,kk) = eye * a13(ii,jj,kk) * kz(kk)
        a21(ii,jj,kk) = eye * a23(ii,jj,kk) * kx(ii)
        a22(ii,jj,kk) = eye * a23(ii,jj,kk) * ky(jj)
        a23(ii,jj,kk) = eye * a23(ii,jj,kk) * kz(kk)
        a31(ii,jj,kk) = eye * a33(ii,jj,kk) * kx(ii)
        a32(ii,jj,kk) = eye * a33(ii,jj,kk) * ky(jj)
        a33(ii,jj,kk) = eye * a33(ii,jj,kk) * kz(kk)
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

      ! - Tr AA
      p = - cmplx( &
            real(a11)*real(a11) + real(a22)*real(a22) + real(a33)*real(a33) &
            + 2 * ( real(a12)*real(a21) + real(a13)*real(a31) + &
            real(a23)*real(a32) ), &
            aimag(a11)*aimag(a11) + aimag(a22)*aimag(a22) + aimag(a33)*aimag(a33) &
            + 2 * ( aimag(a12)*aimag(a21) + aimag(a13)*aimag(a31) + &
            aimag(a23)*aimag(a32) ) )

      call rfftwnd_f77_one_real_to_complex(r2c3d,p,ignore_me)

      p = - p * const / k2       ! Pressure in Fourier space

      write(15) p
      close(15)

      99 if ( ioinfo .ne. 0 ) then
        write(*,*) 'Pressure data file already exists. Skip to next one.'
      end if

    end do
  close(25)

  deallocate( a11, a12, a13, a21, a22, a23, a31, a32, a33, p, kx, ky, kz, k2)

  call destroyplan3d

  write(*,*) 'Done'

end program vel2p
