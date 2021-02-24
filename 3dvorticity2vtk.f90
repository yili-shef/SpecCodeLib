program otovtk
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: ii, jj, kk, nx, ny, nz, lx, lx1, ly, lz, ll, ndel
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  real(sp),    allocatable, dimension(:,:,:) :: g, onorm, rux, ruy, ruz
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  character(80) :: fnm, str, fpath
  real(sp) :: ignore_me, delta_c
  complex(sp) :: ox, oy
  integer :: vtkopenfile, vtk_sp_header, vtk3dvtr, vtkscalar, vtkclosefile


  write(*,*) 
  write(*,'(''>>>>>> Generate 3d vorticity fields for paraview visualization <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./3dvorticity2vtk.x nx filelist ndel'
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
  allocate( rux(lx1,ly,lz), ruy(lx1,ly,lz), ruz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), onorm(nx,ny,nz) )

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
        read(10) wx
      close(10)
      fpath='./out/uy'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10) wy
      close(10)
      fpath='./out/uz'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10) wz
      close(10)
      write(*,*) 'after reading data files'

      wx = wx * g; wy = wy * g; wz = wz * g

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ox = eye*(ky(jj)*wz(ii,jj,kk)-kz(kk)*wy(ii,jj,kk))
        oy = eye*(kz(kk)*wx(ii,jj,kk)-kx(ii)*wz(ii,jj,kk))
        wz(ii,jj,kk) = eye*(kx(ii)*wy(ii,jj,kk)-ky(jj)*wx(ii,jj,kk))
        wx(ii,jj,kk) = ox
        wy(ii,jj,kk) = oy
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

      onorm(1:nx:2,:,:) =  real( wx(1:lx,:,:), sp ) * real( wx(1:lx,:,:), sp )  &
                         + real( wy(1:lx,:,:), sp ) * real( wy(1:lx,:,:), sp )  &
                         + real( wz(1:lx,:,:), sp ) * real( wz(1:lx,:,:), sp )

      onorm(2:nx:2,:,:) =  aimag( wx(1:lx,:,:) ) * aimag( wx(1:lx,:,:) )  &
                         + aimag( wy(1:lx,:,:) ) * aimag( wy(1:lx,:,:) )  &
                         + aimag( wz(1:lx,:,:) ) * aimag( wz(1:lx,:,:) )
      onorm = sqrt(onorm)

      rux(1:nx:2,:,:) =  real( wx(1:lx,:,:), sp )
      rux(2:nx:2,:,:) = aimag( wx(1:lx,:,:)     )
      ruy(1:nx:2,:,:) =  real( wy(1:lx,:,:), sp )
      ruy(2:nx:2,:,:) = aimag( wy(1:lx,:,:)     )
      ruz(1:nx:2,:,:) =  real( wz(1:lx,:,:), sp )
      ruz(2:nx:2,:,:) = aimag( wz(1:lx,:,:)     )

      ! Change suffix and appending 0 at the end for C output.
      fpath = "onorm-"//str(1:len_trim(str)-3)//"vtk0"
      ll = vtkopenfile (fpath( 1:len_trim(fpath) )) 
      ll = vtk_sp_header (sp, nz, ny, nx, 0._sp, 0._sp, 0._sp, 1._sp, 1._sp, 1._sp)
      fpath = "onorm0"
      ll = vtkscalar (onorm, sp, nz, ny, nx, fpath(1:len_trim(fpath)) )
      fpath = "omega0"
      ll = vtk3dvtr (ruz, ruy, rux, sp, nz, ny, nx, fpath(1 : len_trim(fpath)) )
      ll = vtkclosefile()

    end do

  close(20)

  call destroyplan3d

  deallocate(wx, wy, wz, onorm, rux, ruy, ruz, g, kx, ky, kz)

  write(*,*) 'Finished'

end program otovtk
