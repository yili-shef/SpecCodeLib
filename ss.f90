#include "constant.f90"
#include "fftw3.f90"
#include "wavenumber.f90"

program otovtk
  use mconstant
  use mfftw3
  use mwavenumber
  implicit none

  integer :: ii, jj, kk, nx, ny, nz, lx, lx1, ly, lz, ll, ndel
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz, ux, uy, uz
  real(sp),    allocatable, dimension(:,:,:) :: g, rux, ruy, ruz, arrsect
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: fnm, str, fpath, str1, str2
  real(sp) :: ignore_me, tmp, delta_c, dx

  integer(8) :: wxc2r, wyc2r, wzc2r
  integer(8) :: uxc2r, uyc2r, uzc2r


  write(*,*) 
  write(*,'(''>>> Generate 3d ss fields <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./ss.x nx filelist ndel '
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

  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
  allocate( g(lx1,ly,lz) )
  allocate( rux(nx,ny,nz), ruy(nx,ny,nz), ruz(nx,ny,nz) )
  allocate( arrsect(nx,ny,nz) )
  write(*,*) 'arrays allocated'

  call dfftwplan3dc2r(wx, nx, wxc2r)
  call dfftwplan3dc2r(wy, nx, wyc2r)
  call dfftwplan3dc2r(wz, nx, wzc2r)
  call dfftwplan3dc2r(ux, nx, uxc2r)
  call dfftwplan3dc2r(uy, nx, uyc2r)
  call dfftwplan3dc2r(uz, nx, uzc2r)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=fnm(1:len_trim(fnm))//'.list')
    do while ( .true. )
      read(20,*, end=900) str1
      write(*,*) str1(1:len_trim(str1))

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

      ux = ux * g; uy = uy * g; uz = uz * g

      wx = 0.5 * eye * (kx * uy + ky * ux) ! s12
      wy = 0.5 * eye * (kx * uz + kz * ux) ! s13
      wz = 0.5 * eye * (ky * uz + kz * uy) ! s23
      ux = eye * kx * ux ! s11
      uy = eye * ky * uy ! s22
      uz = eye * kz * uz ! s33

      call dfftw_execute(wxc2r)
      call dfftw_execute(wyc2r)
      call dfftw_execute(wzc2r)
      call dfftw_execute(uxc2r)
      call dfftw_execute(uyc2r)
      call dfftw_execute(uzc2r)

      rux(1:nx:2,:,:) = real(  ux(1:lx,:,:) )
      rux(2:nx:2,:,:) = aimag( ux(1:lx,:,:) )
      ruy(1:nx:2,:,:) = real(  uy(1:lx,:,:) )
      ruy(2:nx:2,:,:) = aimag( uy(1:lx,:,:) )
      ruz(1:nx:2,:,:) = real(  uz(1:lx,:,:) )
      ruz(2:nx:2,:,:) = aimag( uz(1:lx,:,:) )

      arrsect = rux * rux + ruy * ruy + ruz * ruz 

      rux(1:nx:2,:,:) = real(  wx(1:lx,:,:) )
      rux(2:nx:2,:,:) = aimag( wx(1:lx,:,:) )
      ruy(1:nx:2,:,:) = real(  wy(1:lx,:,:) )
      ruy(2:nx:2,:,:) = aimag( wy(1:lx,:,:) )
      ruz(1:nx:2,:,:) = real(  wz(1:lx,:,:) )
      ruz(2:nx:2,:,:) = aimag( wz(1:lx,:,:) )

      arrsect = arrsect + 2 * (rux * rux + ruy * ruy + ruz * ruz)
      arrsect = sqrt(arrsect)

      fpath = "./out/ss_"//trim(str)//"dx"//trim(str1)
      open(11, file = fpath, form='binary')
        write(11) arrsect
      close(11)

    end do
    900 continue
  close(20)

  deallocate(wx, wy, wz, ux, uy, uz, arrsect)
  deallocate(rux, ruy, ruz, kx, ky, kz, g)

  write(*,*) 'Finished'

end program otovtk
