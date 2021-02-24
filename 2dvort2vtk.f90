program twodvortvtk
  use mconstant
  use mfftwplan2d
  use mwavenumber
  implicit none

  integer :: ii, jj, nx, ny, lx, lx1, ly, ll, npart
  complex(sp), allocatable, dimension(:,:) :: ux, uy
  real(sp),    allocatable, dimension(:,:) :: wz
  real(sp),    allocatable, dimension(:) :: kx, ky, xco, yco

  character(80) :: fnm, str, fpath
  real(sp) :: ignore_me
  integer :: vtkopenfile, vtk_rg_header, vtkscalar, vtkclosefile


  write(*,*) 
  write(*,'(''>>>>>> Generate 3d vorticity fields for paraview visualization <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dvort2vtk.x nx filelist npart'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: .list file of data files'
          write(*,*) '        npart: 1/npart/npart of whole domain'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! npart
  call getarg(3,fnm)
  read(fnm, '(I20)') npart

  if ( (nx/npart) * npart .ne. nx ) stop "nx/npart is not an integer"

  ! file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny=nx; 
  lx=nx/2; lx1=lx+1
  ly=nx; 


  call fftwplan2de(nx,ny)
  write(*,*) 'after fftwplan3d'

  allocate( wz(nx/npart,ny/npart), ux(lx1,ly), uy(lx1,ly) )
  allocate( kx(lx1), ky(ly), xco(nx/npart), yco(ny/npart) )

  xco = (/ ( (ii-1)*2*pi/nx, ii = 1, nx/npart) /)
  yco = (/ ( (jj-1)*2*pi/ny, jj = 1, ny/npart) /)

  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,lx1,ly)
  write(*,*) 'after wavenumber'

  open(20,file=fnm(1:len_trim(fnm))//'.list')
    do while ( .not. eof(20)) 
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      fpath='./out/ux'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)ux
      close(10)
      fpath='./out/uy'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      write(*,*) 'after reading data files'


      do jj = 1, ly
      do ii = 1, lx1
        ux(ii,jj) = eye*( kx(ii) * uy(ii,jj) - ky(jj) * ux(ii,jj) )
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r2d,ux,ignore_me)

      wz(1:nx/npart:2,:) = real(  ux(1:lx/npart,:) )
      wz(2:nx/npart:2,:) = aimag( ux(1:lx/npart,:) )

      ! Change suffix and appending 0 at the end for C output.
      ll = npart*npart
      write(fpath,*) ll
      fpath = adjustl(fpath)

      fpath = "2dvort-1of"//fpath( 1:len_trim(fpath) )//"-"//str( 1:len_trim(str)-3 )//"vtk0"
      ii = nx/npart;    jj = ny/npart

      ll = vtkopenfile (fpath( 1:len_trim(fpath) )) 
      ll = vtk_rg_header(xco, yco, 0, sp, ii, jj, 1)

      fpath = "vorticity0"
      fpath = adjustl(fpath)
      ll = vtkscalar (wz, sp, ii, jj, 1, fpath( 1:len_trim(fpath) ) )
      ll = vtkclosefile()

    end do
  close(20)

  call destroyplan2d

  deallocate(wz, ux, uy, kx, ky, xco, yco)

  write(*,*) 'Finished'

end program twodvortvtk
