program ggspec
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nfile

  complex(sp), allocatable, dimension(:,:,:) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz
  real(dp),    allocatable, dimension(:)     :: ekcurv

  real(dp) :: const
  real(sp) :: kxii, kyjj, kzkk
  real(sp) :: delta_c, ignore_me, tmp
  character(80) :: fnm, str, fpath

  write(*,*) 
  write(*,'(''>>>>>> Mean rms and PDFs of Curvature of vortex lines <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./gijgij-spec.x nx filelist ndel'
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

  allocate( kx(lx1), ky(ly), kz(lz), ekcurv(lx) )
  allocate( g(lx1,ly,lz) )
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz), g21(lx1,ly,lz), g22(lx1,ly,lz) )
  allocate( g23(lx1,ly,lz), g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  ekcurv=0._dp
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

      ! gij 
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

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx

      delta_c = real( g11(ii,jj,kk) ) ** 2 &
              + real( g12(ii,jj,kk) ) ** 2 &
              + real( g13(ii,jj,kk) ) ** 2 &
              + real( g21(ii,jj,kk) ) ** 2 &
              + real( g22(ii,jj,kk) ) ** 2 &
              + real( g23(ii,jj,kk) ) ** 2 &
              + real( g31(ii,jj,kk) ) ** 2 &
              + real( g32(ii,jj,kk) ) ** 2 &
              + real( g33(ii,jj,kk) ) ** 2 

      ignore_me = aimag( g11(ii,jj,kk) ) ** 2 &
                + aimag( g12(ii,jj,kk) ) ** 2 &
                + aimag( g13(ii,jj,kk) ) ** 2 &
                + aimag( g21(ii,jj,kk) ) ** 2 &
                + aimag( g22(ii,jj,kk) ) ** 2 &
                + aimag( g23(ii,jj,kk) ) ** 2 &
                + aimag( g31(ii,jj,kk) ) ** 2 &
                + aimag( g32(ii,jj,kk) ) ** 2 &
                + aimag( g33(ii,jj,kk) ) ** 2 

      g11(ii,jj,kk) = cmplx( delta_c, ignore_me )

    end do
    end do
    end do

    g11 = g11 * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,g11,ignore_me)

    do kk=1,lz
    
      kzkk=kz(kk)
      do jj=1,ly
    
        kyjj=ky(jj)
        do ii=1,lx1
    
          kxii=kx(ii)
    
          ignore_me = real( g11(ii,jj,kk) * conjg( g11(ii,jj,kk) ) )

          if ( ii .eq. 1 ) ignore_me = .5 * ignore_me

          ll=floor(sqrt(kxii*kxii+kyjj*kyjj+kzkk*kzkk)+.5)
          if (ll .ge. 1 .and. ll .le. lx) then
                  ekcurv(ll)=ekcurv(ll)+ignore_me
          end if
    
        end do
      end do
    end do

    nfile = nfile + 1
  end do
  close(20)


  ! Ouput format for gnuplot
  open(24,file='spec-gijgij'//str(1:len_trim(str)))
  do ii=1,lx
    write(24,*) ii,ekcurv(ii)/nfile
  end do
  close(24)


  deallocate(kx,ky,kz,g,g11,g12,g13,g21,g22,g23,g31,g32,g33)
  deallocate(ekcurv)

  call destroyplan3d

  write(*,*) 'Finished'
end program ggspec
