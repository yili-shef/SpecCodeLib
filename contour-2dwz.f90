program vorticity
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz, iii,jjj,kkk
  real :: ignore_me
  real, allocatable, dimension(:,:,:) :: kx,ky,kz,g
  complex, allocatable, dimension(:,:,:) :: ux,uy,wz
  character(80) :: str,flnm, path
  
  nx=iargc()
  if (nx .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./contour-2dwz.x nx nfile iii idir ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        nfile: data file number'
          write(*,*) '        iii: location of the contour plane is iii * dx'
          write(*,*) '        idir: 1 for (y,z) plane, 2 (x,z), 3 (x,y)'
          write(*,*) '        ndel: filter scale = ndel * dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file number
  call getarg(2,flnm)
  flnm=adjustl(flnm)

  ! iii
  call getarg(3,str)
  read(str, '(I20)') iii

  ! idir
  call getarg(4,str)
  read(str, '(I20)') kkk

  ! ndel
  call getarg(5,str)
  read(str, '(I20)') ndel

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

  write(*,*) 'fftwplan'
  call fftwplan3de(nx,ny,nz)
  
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz))
  allocate(wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz),g(lx1,ly,lz))
  
  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  
  open(15,file='./out/ux'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) ux
  close(15)
  open(15,file='./out/uy'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) uy
  close(15)
  
  
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wz(ii,jj,kk)=eye*(kx(ii)*uy(ii,jj,kk)-ky(jj)*ux(ii,jj,kk))
  end do
  end do
  end do
  
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
  
  open(16,file='wz-xyplane'//str(1:len_trim(str))//'file'//flnm(1:len_trim(flnm))//'.dat')
    write(16,*) 'Variables = "x" "y" "wz" '
    write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,ny
    do jj=1,ny
    do ii=1,nx
      write(16, '(2I5,5e18.3)') ii,jj, wz(ii,jj,kkk)
    end do
    end do
  close(16)
  
  deallocate(kx,ky,kz,k2)
  deallocate(ux,uy,wz)

end program vorticity
