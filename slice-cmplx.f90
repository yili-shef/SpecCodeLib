program slice
use mfftwplan3d
implicit none

integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz
real :: ignore_me
real, allocatable, dimension(:,:,:) :: uxr,uyr,uzr
complex, allocatable, dimension(:,:,:) :: ux,uy,uz
character(80) :: str,flnm, path

write(*,*) 'data resolution nx:'
read(*,*) nx

ny=nx; nz=nx
lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

write(*,*) 'fftwplan'
call fftwplan3d(nx,ny,nz)

write(*,*) 'file number :'
read(*,*) ii
write(flnm, '(I6)') ii
flnm=adjustl(flnm)

allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
allocate(uxr(nx,ny,nz),uyr(nx,ny,nz),uzr(nx,ny,nz))
open(15,file='./out/ux'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
  read(15) ux
close(15)
open(15,file='./out/uy'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
  read(15) uy
close(15)
open(15,file='./out/uz'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
  read(15) uz
close(15)


call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
uxr(1:nx:2,:,:)= real(ux(1:lx,:,:))
uxr(2:nx:2,:,:)=aimag(ux(1:lx,:,:))
uyr(1:nx:2,:,:)= real(uy(1:lx,:,:))
uyr(2:nx:2,:,:)=aimag(uy(1:lx,:,:))
uzr(1:nx:2,:,:)= real(uz(1:lx,:,:))
uzr(2:nx:2,:,:)=aimag(uz(1:lx,:,:))

write(*,*) 'input kk='
read(*,*) kk
if (kk .ne. 0) then
write(str,'(I7)') kk
str=adjustl(str)

open(16,file='uxy-slice'//str(1:len_trim(str))//'.dat')
  write(16,*) 'Variables = "x" "y" "ux" "uy"'
  write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,ny
  do jj=1,ny
  do ii=1,nx
    write(16, '(2I5,5e18.3)') ii,jj, uxr(ii,jj,kk),uyr(ii,jj,kk)
  end do
  end do
close(16)

end if

write(*,*) 'input jj='

read(*,*) jj
if (jj .ne. 0) then
write(str,'(I7)') jj
str=adjustl(str)

open(16,file='uxz-slice'//str(1:len_trim(str))//'.dat')
  write(16,*) 'Variables = "x" "z" "ux" "uz"'
  write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,ny
  do kk=1,nz
  do ii=1,nx
    write(16, '(2I5,5e18.3)') ii,kk, uxr(ii,jj,kk),uzr(ii,jj,kk)
  end do
  end do
close(16)

end if

write(*,*) 'input ii='

read(*,*) ii
if (ii .ne. 0) then
write(str,'(I7)') ii
str=adjustl(str)

open(16,file='uyz-slice'//str(1:len_trim(str))//'.dat')
  write(16,*) 'Variables = "y" "z" "uy" "uz"'
  write(16,'(''zone i='',i6,'' j='',i6,'',f=point'')') nx,ny
  do kk=1,nz
  do jj=1,nx
    write(16, '(2I5,5e18.3)') jj,kk, uyr(ii,jj,kk),uzr(ii,jj,kk)
  end do
  end do
close(16)

end if



deallocate(ux,uy,uz)
deallocate(uxr,uyr,uzr)

end program slice
