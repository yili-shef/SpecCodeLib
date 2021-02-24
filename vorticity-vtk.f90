program otovtk
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: ii, jj, kk, nx, ny, nz, lx, lx1, ly, lz, ll, ndel
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz, 
  real(sp),    allocatable, dimension(:,:,:) :: g, wnorm, rux, ruy, ruz
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  character(80) :: fnm, str, fpath
  real(sp) :: ignore_me, tmp, delta_c
  integer :: vtkopenfile, vtkheader, vtk3dvtr, vtkscalar, vtkclosefile


  write(*,*) 
  write(*,'(''>>>>>> Generate 3d vorticity fields for paraview visualization <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./vorticity-vtk.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: .list file of data files'
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
  ! file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( xix(lx1,ly,lz), xiy(lx1,ly,lz), xiz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), xinorm(nx,ny,nz) )
  allocate( rux(nx,ny,nz), ruy(nx,ny,nz), ruz(nx,ny,nz) )

  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=fnm(1:len_trim(fnm))//'.list')
    do while ( .not. eof(20)) 
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      fpath='./out/ux'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wx
      close(10)
      fpath='./out/uy'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wy
      close(10)
      fpath='./out/uz'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wz
      close(10)
      write(*,*) 'after reading data files'

      wx = wx * g; wy = wy * g; wz = wz * g

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ignore_me = eye*(ky(jj)*wz(ii,jj,kk)-kz(kk)*wy(ii,jj,kk))
        tmp       = eye*(kz(kk)*wx(ii,jj,kk)-kx(ii)*wz(ii,jj,kk))
        wz(ii,jj,kk) = eye*(kx(ii)*wy(ii,jj,kk)-ky(jj)*wx(ii,jj,kk))
        wx(ii,jj,kk) = ignore_me
        wy(ii,jj,kk) = tmp
        ignore_me = eye*(ky(jj)*wz(ii,jj,kk)-kz(kk)*wy(ii,jj,kk))
        tmp       = eye*(kz(kk)*wx(ii,jj,kk)-kx(ii)*wz(ii,jj,kk))
        xiz(ii,jj,kk) = eye*(kx(ii)*wy(ii,jj,kk)-ky(jj)*wx(ii,jj,kk))
        xix(ii,jj,kk) = ignore_me
        xiy(ii,jj,kk) = tmp
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

      rux(1:nx:2,:,:) = real(  wx(1:lx,:,:) )
      rux(2:nx:2,:,:) = aimag( wx(1:lx,:,:) )
      ruy(1:nx:2,:,:) = real(  wy(1:lx,:,:) )
      ruy(2:nx:2,:,:) = aimag( wy(1:lx,:,:) )
      ruz(1:nx:2,:,:) = real(  wz(1:lx,:,:) )
      ruz(2:nx:2,:,:) = aimag( wz(1:lx,:,:) )

      call rfftwnd_f77_one_complex_to_real(c2r3d,xix,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,xiy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,xiz,ignore_me)

      xinorm(1:nx:2,:,:) = real( xix(1:lx,:,:) ) * real( xix(1:lx,:,:) )  &
                         + real( xiy(1:lx,:,:) ) * real( xiy(1:lx,:,:) )  &
                         + real( xiz(1:lx,:,:) ) * real( xiz(1:lx,:,:) )

      xinorm(2:nx:2,:,:) = aimag( xix(1:lx,:,:) ) * aimag( xix(1:lx,:,:) )  &
                         + aimag( xiy(1:lx,:,:) ) * aimag( xiy(1:lx,:,:) )  &
                         + aimag( xiz(1:lx,:,:) ) * aimag( xiz(1:lx,:,:) )
      xinorm = sqrt(xinorm)

      ! Change suffix and appending 0 at the end for C output.
      fpath = "ovtr-xinorm-"//str(1:len_trim(str)-3)//"vtk0"
      ll = vtkopenfile (fpath( 1:len_trim(fpath) )) 
      ll = vtkheader (nx)
      ll = vtkscalar (xinorm, nx)
      ll = vtk3dvtr (rux, ruy, ruz, nx)
      ll = vtkclosefile()

    end do
  close(20)

  call destroyplan3d

  deallocate(wx, wy, wz, xix, xiy, xiz, xinorm)
  deallocate(rux, ruy, ruz)

  write(*,*) 'Finished'

end program otovtk
