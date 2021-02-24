program msgheli
  use mconstant
  use mfftwplan3d
  implicit none

  integer, parameter :: nx=256, neta=80
  integer, parameter :: ny=nx,nz=nx,lx=nx/2,lx1=lx+1,ly=ny,lz=nz

  real, parameter :: alpha=2., rnu=0.0015, eps=0.1, eta=(rnu**3/eps)**.25

  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz,ux0,uy0,uz0
  complex, allocatable, dimension(:,:,:) :: tij,rij
  complex, allocatable, dimension(:,:,:) :: helidiss,enerdiss
  real, allocatable, dimension(:,:,:) :: k2,g,g0
  real, allocatable, dimension(:) :: kx,ky,kz

  character(80) :: fnm,fpath
  integer :: ifile,ndel,ii,jj,kk,ll
  real :: delta_c, del, ignore_me,rr,ss

  ifile=5

  delta_c=neta*eta
  ndel=floor(log(delta_c*nx/(2.*pi))/log(2.))
  write(*,*) 'ndel = ', ndel

  write(fnm,'(i30)') ifile
  fnm=adjustl(fnm)

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(ux0(lx1,ly,lz),uy0(lx1,ly,lz),uz0(lx1,ly,lz))
  allocate(tij(lx1,ly,lz),g(lx1,ly,lz))
  allocate(rij(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))
  allocate(helidiss(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

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

  do ll=0,ndel

    del=delta_c/2**ll
    g=exp(-k2*del**2/24.)
    ux=ux0*g
    uy=uy0*g 
    uz=uz0*g
    
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
 
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*kx(ii)*wx(ii,jj,kk)
      tij(ii,jj,kk)=eye*kx(ii)*ux(ii,jj,kk)
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))

    rr=sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*ky(jj)*wy(ii,jj,kk)
      tij(ii,jj,kk)=eye*ky(ii)*uy(ii,jj,kk)
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=helidiss+cmplx(real(rij)*real(tij), aimag(rij)*aimag(tij))

    rr=rr+sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=ss+sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*kz(kk)*wz(ii,jj,kk)
      tij(ii,jj,kk)=eye*kz(kk)*uz(ii,jj,kk)
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=helidiss+cmplx(real(rij)*real(tij), aimag(rij)*aimag(tij))

    rr=rr+sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=ss+sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*.5*(kx(ii)*wy(ii,jj,kk)+ky(jj)*wx(ii,jj,kk))
      tij(ii,jj,kk)=eye*.5*(kx(ii)*uy(ii,jj,kk)+ky(jj)*ux(ii,jj,kk))
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=helidiss+2.*cmplx(real(rij)*real(tij), aimag(rij)*aimag(tij))

    rr=rr+2.*sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+2*sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=ss+2.*sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+2*sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*.5*(kx(ii)*wz(ii,jj,kk)+kz(kk)*wx(ii,jj,kk))
      tij(ii,jj,kk)=eye*.5*(kx(ii)*uz(ii,jj,kk)+kz(kk)*ux(ii,jj,kk))
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=helidiss+2.*cmplx(real(rij)*real(tij), aimag(rij)*aimag(tij))

    rr=rr+2.*sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+2*sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=ss+2.*sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+2*sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      rij(ii,jj,kk)=eye*.5*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
      tij(ii,jj,kk)=eye*.5*(ky(jj)*uz(ii,jj,kk)+kz(kk)*uy(ii,jj,kk))
    end do
    end do
    end do
    call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
    helidiss=helidiss+2.*cmplx(real(rij)*real(tij), aimag(rij)*aimag(tij))
 
    rr=rr+2.*sum(real(rij(1:lx,:,:))*real(rij(1:lx,:,:)))+2.*sum(aimag(rij(1:lx,:,:))*aimag(rij(1:lx,:,:)))
    ss=ss+2.*sum(real(tij(1:lx,:,:))*real(tij(1:lx,:,:)))+2.*sum(aimag(tij(1:lx,:,:))*aimag(tij(1:lx,:,:)))

    rr=rr/(nx*ny*nz); ss=ss/(nx*ny*nz)

    write(*,*) (sum(real(helidiss(1:lx,:,:)))+sum(aimag(helidiss(1:lx,:,:))))/(nx*ny*nz)/sqrt(rr)/sqrt(ss), &
               sqrt(rr), sqrt(ss)

  end do

  deallocate(wx,wy,wz,kx,ky,kz,k2)
  deallocate(ux0,uy0,uz0,ux,uy,uz)
  deallocate(tij,g,helidiss)
  deallocate(rij)

  call destroyplan3d

end program msgheli 
