program vorticity
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz, iii,jjj,kkk
  real :: ignore_me
  real, allocatable, dimension(:,:,:) :: uxr,uyr,uzr,kx,ky,kz,k2,uxn,uyn,uzn
  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz
  character(80) :: str,flnm, path
  
  nx=iargc()
  if (nx .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: vorticity-contour.x nx iii jjj kkk nfile'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        iii,jjj,kkk: labels of (y,z),(x,z),(x,y) planes. 0=skip'
          write(*,*) '        nfile: data file number'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! (iii,jjj,kkk)
  call getarg(2,str)
  read(str, '(I20)') iii
  call getarg(3,str)
  read(str, '(I20)') jjj
  call getarg(4,str)
  read(str, '(I20)') kkk

  ! file number
  call getarg(5,flnm)
  flnm=adjustl(flnm)

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

  write(*,*) 'fftwplan'
  call fftwplan3d(nx,ny,nz)
  
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(uxr(nx,ny,nz),uyr(nx,ny,nz),uzr(nx,ny,nz))
  allocate(uxn(nx,ny,nz),uyn(nx,ny,nz),uzn(nx,ny,nz))
  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  
  open(15,file='./out/ux'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) ux
  close(15)
  open(15,file='./out/uy'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) uy
  close(15)
  open(15,file='./out/uz'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) uz
  close(15)
  
  
  
  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  
  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
  uxr(1:nx:2,:,:)= real(wx(1:lx,:,:))
  uxr(2:nx:2,:,:)=aimag(wx(1:lx,:,:))
  uyr(1:nx:2,:,:)= real(wy(1:lx,:,:))
  uyr(2:nx:2,:,:)=aimag(wy(1:lx,:,:))
  uzr(1:nx:2,:,:)= real(wz(1:lx,:,:))
  uzr(2:nx:2,:,:)=aimag(wz(1:lx,:,:))
  
  call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
  uxn(1:nx:2,:,:)= real(ux(1:lx,:,:))
  uxn(2:nx:2,:,:)=aimag(ux(1:lx,:,:))
  uyn(1:nx:2,:,:)= real(uy(1:lx,:,:))
  uyn(2:nx:2,:,:)=aimag(uy(1:lx,:,:))
  uzn(1:nx:2,:,:)= real(uz(1:lx,:,:))
  uzn(2:nx:2,:,:)=aimag(uz(1:lx,:,:))
  
  if (kkk .ne. 0) then
  write(str,'(I7)') kkk
  str=adjustl(str)
  
  open(16,file='wz-xyplane'//str(1:len_trim(str))//'file'//flnm(1:len_trim(flnm))//'.dat')
    write(16,*) 'Variables = "x" "y" "wz" "ux" "uy" '
    write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,ny
    do jj=1,ny
    do ii=1,nx
      write(16, '(2I5,5e18.3)') ii,jj, uzr(ii,jj,kkk),uxn(ii,jj,kkk),uyn(ii,jj,kkk)
    end do
    end do
  close(16)
  
  end if
  
  if (jjj .ne. 0) then
  write(str,'(I7)') jjj
  str=adjustl(str)
  
  open(16,file='wy-xzplane'//str(1:len_trim(str))//'file'//flnm(1:len_trim(flnm))//'.dat')
    write(16,*) 'Variables = "x" "z" "wy" "ux" "uz" '
    write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,nz
    do kk=1,nz
    do ii=1,nx
      write(16, '(2I5,5e18.3)') ii,kk, uyr(ii,jjj,kk),uxn(ii,jjj,kk),uzn(ii,jjj,kk)
    end do
    end do
  close(16)
  
  end if
  
  if (iii .ne. 0) then
  write(str,'(I7)') iii
  str=adjustl(str)
  
  open(16,file='wx-yzplane'//str(1:len_trim(str))//'file'//flnm(1:len_trim(flnm))//'.dat')
    write(16,*) 'Variables = "y" "z" "wx" "uy" "uz" '
    write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') ny,nz
    do kk=1,nz
    do jj=1,ny
      write(16, '(2I5,5e18.3)') jj,kk, uxr(iii,jj,kk),uyn(iii,jj,kk),uzn(iii,jj,kk)
    end do
    end do
  close(16)
  
  end if
  
  deallocate(ux,uy,uz)
  deallocate(uxr,uyr,uzr)
  deallocate(uxn,uyn,uzn)
  deallocate(kx,ky,kz,k2)
  deallocate(wx,wy,wz)

end program vorticity
