program phivtk
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ll, ndel
  complex(sp), allocatable, dimension(:,:,:) :: phi
  real(sp),    allocatable, dimension(:,:,:) :: g, rphi 
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  character(80) :: fnm, str, fpath
  real(sp) :: ignore_me, delta_c
  integer :: vtkopenfile, vtk_sp_header, vtkscalar, vtkclosefile


  write(*,*) 
  write(*,'(''>>>>>> Generate 3d vorticity fields for paraview visualization <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./phi-vtk.x nx filelist ndel'
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

  allocate( phi(lx1,ly,lz), rphi(nx,ny,nz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz) )

  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=fnm(1:len_trim(fnm))//'.list')
    do while ( .not. eof(20)) 
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      fpath='./out/phi'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)phi
      close(10)
      write(*,*) 'after reading data files'

      phi = phi * g


      call rfftwnd_f77_one_complex_to_real(c2r3d,phi,ignore_me)

      rphi(1:nx:2,:,:) = real(  phi(1:lx,:,:) )
      rphi(2:nx:2,:,:) = aimag( phi(1:lx,:,:) )

      ! Change suffix and appending 0 at the end for C output.
      fpath = "phi-"//str(1:len_trim(str)-3)//"vtk0"

      ll = vtkopenfile (fpath( 1:len_trim(fpath) )) 
      ll = vtk_sp_header (sp, nz, ny, nx, 0._sp, 0._SP, 0._sp, 1._sp, 1._sp, 1._sp)

      fpath = "phi0"
      ll = vtkscalar (rphi, sp, nz, ny, nx, fpath(1 : len_trim(fpath)) )

      ll = vtkclosefile()

    end do
  close(20)

  call destroyplan3d

  deallocate(phi, rphi)
  deallocate(kx, ky, kz, g)

  write(*,*) 'Finished'

end program phivtk
