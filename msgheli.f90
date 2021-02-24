program msgheli
  use mconstant
  use mfftwplan3d
  implicit none

  integer, parameter :: nx=256, neta=80
  integer, parameter :: ny=nx,nz=nx,lx=nx/2,lx1=lx+1,ly=ny,lz=nz

  real, parameter :: alpha=2., rnu=0.0015, eps=0.1, eta=(rnu**3/eps)**.25

  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz,ux0,uy0,uz0
  complex, allocatable, dimension(:,:,:) :: tij
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  complex, allocatable, dimension(:,:,:) :: helidiss,enerdiss
  real, allocatable, dimension(:,:,:) :: k2,g,g0
  real, allocatable, dimension(:) :: kx,ky,kz

  character(80) :: fnm,fpath
  integer :: ifile,ndel,ii,jj,kk
  real :: delta_c, del, ignore_me

  ifile=5

  delta_c=neta*eta
  ndel=floor(log(delta_c*nx/(2.*pi))/log(2.))
  write(*,*) 'ndel = ', ndel

  write(fnm,'(i30)') ifile
  fnm=adjustl(fnm)

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(ux0(lx1,ly,lz),uy0(lx1,ly,lz),uz0(lx1,ly,lz))
  allocate(tij(lx1,ly,lz),g(lx1,ly,lz),g0(lx1,ly,lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))

  allocate(helidiss(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)ux0
  close(10)
  fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)uy0
  close(10)
  fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)uz0
  close(10)
  write(*,*) 'after reading data files'

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  g=exp(-k2*delta_c**2/24.)
  write(*,*) 'after filter'

  ux=ux0*g; uy=uy0*g; uz=uz0*g

  do kk=1,lz
  do ii=1,lx1
    wx(ii,:,kk)=ky(:)*uz(ii,:,kk)
    wz(ii,:,kk)=ky(:)*ux(ii,:,kk)
  end do
  end do
  do jj=1,ly
  do ii=1,lx1
    wx(ii,jj,:)=wx(ii,jj,:)-kz(:)*uy(ii,jj,:)
    wy(ii,jj,:)=kz(:)*ux(ii,jj,:)
  end do
  end do
  do kk=1,lz
  do jj=1,ly
    wy(:,jj,kk)=wy(:,jj,kk)-kx(:)*uz(:,jj,kk)
    wz(:,jj,kk)=kx(:)*uy(:,jj,kk)-wz(:,jj,kk)
  end do
  end do
  wx=eye*wx; wy=eye*wy; wz=eye*wz
  write(*,*) 'after vorticity'

  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    r11(ii,jj,kk)=eye*kx(ii)*wx(ii,jj,kk)
    r22(ii,jj,kk)=eye*ky(jj)*wy(ii,jj,kk)
    r33(ii,jj,kk)=eye*kz(kk)*wz(ii,jj,kk)
    r12(ii,jj,kk)=eye*.5*(kx(ii)*wy(ii,jj,kk)+ky(jj)*wx(ii,jj,kk))
    r13(ii,jj,kk)=eye*.5*(kx(ii)*wz(ii,jj,kk)+kz(kk)*wx(ii,jj,kk))
    r23(ii,jj,kk)=eye*.5*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
  end do
  end do
  end do
  write(*,*) 'after rij'

  call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
  write(*,*) 'after fftw rij'

  call rsgstauij(ux0,ux0,tij,g,nx,ny,nz) 
  helidiss=cmplx( real(r11)*real(tij), aimag(r11)*aimag(tij))
  call rsgstauij(uy0,uy0,tij,g,nx,ny,nz)
  helidiss=helidiss+cmplx(real(r22)*real(tij), aimag(r22)*aimag(tij))
  call rsgstauij(uz0,uz0,tij,g,nx,ny,nz)
  helidiss=helidiss+cmplx(real(r33)*real(tij), aimag(r33)*aimag(tij))
  call rsgstauij(ux0,uy0,tij,g,nx,ny,nz) 
  helidiss=helidiss+2.*cmplx(real(r12)*real(tij), aimag(r12)*aimag(tij))
  call rsgstauij(ux0,uz0,tij,g,nx,ny,nz)
  helidiss=helidiss+2.*cmplx(real(r13)*real(tij), aimag(r13)*aimag(tij))
  call rsgstauij(uy0,uz0,tij,g,nx,ny,nz)
  helidiss=helidiss+2.*cmplx(real(r23)*real(tij), aimag(r23)*aimag(tij))
  write(*,*) 'after helidiss 0'

  write(*,*) -2.*(sum(real(helidiss(1:lx,:,:)))+sum(aimag(helidiss(1:lx,:,:))))/(nx*ny*nz)

  do ii=0,ndel

    del=delta_c/2**ii
    g0=exp(-k2*del**2/24.)
    ux=ux0*g0
    uy=uy0*g0 
    uz=uz0*g0
    
    call rsgstauij(ux,ux,tij,g,nx,ny,nz) 
    helidiss=cmplx( real(r11)*real(tij), aimag(r11)*aimag(tij))
    call rsgstauij(uy,uy,tij,g,nx,ny,nz)
    helidiss=helidiss+cmplx(real(r22)*real(tij), aimag(r22)*aimag(tij))
    call rsgstauij(uz,uz,tij,g,nx,ny,nz)
    helidiss=helidiss+cmplx(real(r33)*real(tij), aimag(r33)*aimag(tij))
    call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
    helidiss=helidiss+2.*cmplx(real(r12)*real(tij), aimag(r12)*aimag(tij))
    call rsgstauij(ux,uz,tij,g,nx,ny,nz)
    helidiss=helidiss+2.*cmplx(real(r13)*real(tij), aimag(r13)*aimag(tij))
    call rsgstauij(uy,uz,tij,g,nx,ny,nz)
    helidiss=helidiss+2.*cmplx(real(r23)*real(tij), aimag(r23)*aimag(tij))
 
    write(*,*) -2.*(sum(real(helidiss(1:lx,:,:)))+sum(aimag(helidiss(1:lx,:,:))))/(nx*ny*nz)

  end do


  deallocate(helidiss)
  allocate(enerdiss(lx1,ly,lz))

  ux=ux0*g; uy=uy0*g; uz=uz0*g

  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
  r11(ii,jj,kk)=eye*kx(ii)*ux(ii,jj,kk)
  r22(ii,jj,kk)=eye*ky(jj)*uy(ii,jj,kk)
  r33(ii,jj,kk)=eye*kz(kk)*uz(ii,jj,kk)
  r12(ii,jj,kk)=eye*.5*(kx(ii)*uy(ii,jj,kk)+ky(jj)*ux(ii,jj,kk))
  r13(ii,jj,kk)=eye*.5*(kx(ii)*uz(ii,jj,kk)+kz(kk)*ux(ii,jj,kk))
  r23(ii,jj,kk)=eye*.5*(ky(jj)*uz(ii,jj,kk)+kz(kk)*uy(ii,jj,kk))
  end do
  end do
  end do
  write(*,*) 'after sij'

  call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
  write(*,*) 'after fftw rij'

  call rsgstauij(ux0,ux0,tij,g,nx,ny,nz) 
  enerdiss=cmplx( real(r11)*real(tij), aimag(r11)*aimag(tij))
  call rsgstauij(uy0,uy0,tij,g,nx,ny,nz)
  enerdiss=enerdiss+cmplx(real(r22)*real(tij), aimag(r22)*aimag(tij))
  call rsgstauij(uz0,uz0,tij,g,nx,ny,nz)
  enerdiss=enerdiss+cmplx(real(r33)*real(tij), aimag(r33)*aimag(tij))
  call rsgstauij(ux0,uy0,tij,g,nx,ny,nz) 
  enerdiss=enerdiss+2.*cmplx(real(r12)*real(tij), aimag(r12)*aimag(tij))
  call rsgstauij(ux0,uz0,tij,g,nx,ny,nz)
  enerdiss=enerdiss+2.*cmplx(real(r13)*real(tij), aimag(r13)*aimag(tij))
  call rsgstauij(uy0,uz0,tij,g,nx,ny,nz)
  enerdiss=enerdiss+2.*cmplx(real(r23)*real(tij), aimag(r23)*aimag(tij))
  write(*,*) 'after enerdiss 0'

  write(*,*) -(sum(real(enerdiss(1:lx,:,:)))+sum(aimag(enerdiss(1:lx,:,:))))/(nx*ny*nz)

  do ii=0,ndel

    del=delta_c/2**ii
    g0=exp(-k2*del**2/24.)
    ux=ux0*g0
    uy=uy0*g0 
    uz=uz0*g0
    
    call rsgstauij(ux,ux,tij,g,nx,ny,nz) 
    enerdiss=cmplx( real(r11)*real(tij), aimag(r11)*aimag(tij))
    call rsgstauij(uy,uy,tij,g,nx,ny,nz)
    enerdiss=enerdiss+cmplx(real(r22)*real(tij), aimag(r22)*aimag(tij))
    call rsgstauij(uz,uz,tij,g,nx,ny,nz)
    enerdiss=enerdiss+cmplx(real(r33)*real(tij), aimag(r33)*aimag(tij))
    call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
    enerdiss=enerdiss+2.*cmplx(real(r12)*real(tij), aimag(r12)*aimag(tij))
    call rsgstauij(ux,uz,tij,g,nx,ny,nz)
    enerdiss=enerdiss+2.*cmplx(real(r13)*real(tij), aimag(r13)*aimag(tij))
    call rsgstauij(uy,uz,tij,g,nx,ny,nz)
    enerdiss=enerdiss+2.*cmplx(real(r23)*real(tij), aimag(r23)*aimag(tij))
 
    write(*,*) -(sum(real(enerdiss(1:lx,:,:)))+sum(aimag(enerdiss(1:lx,:,:))))/(nx*ny*nz)

  end do

  deallocate(wx,wy,wz,kx,ky,kz,k2)
  deallocate(ux0,uy0,uz0,ux,uy,uz)
  deallocate(tij,enerdiss,g,g0)
  deallocate(r11,r12,r13,r22,r23,r33)

  call destroyplan3d

end program msgheli 
