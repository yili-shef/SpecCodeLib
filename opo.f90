program opo
  use mconstant
  use mfftwplan3d
  implicit none


  complex, allocatable, dimension(:,:,:) :: wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33
  real, allocatable, dimension(:,:,:) :: g
  real, allocatable, dimension(:) :: kx,ky,kz

  integer,  parameter :: npnt = 100
  real(sp), parameter :: bnd = 40, binw=2*bnd/npnt
  real(dp), dimension(npnt) :: popo

  real(sp) :: cc(3,3), om(3)

  character(80) :: fnm,str,fpath
  integer :: ndel,ii,jj,kk,ll,mm,nn
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real(sp) :: delta_c, ignore_me, ompom
  real(dp) :: meanopo, rmsopo

  write(*,*) 
  write(*,'(''>>>>>> PDFs of omega pij omega <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./opo.x nx nfile ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        nfile: data file number'
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

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(p11(lx1,ly,lz),p12(lx1,ly,lz),p13(lx1,ly,lz))
  allocate(p22(lx1,ly,lz),p23(lx1,ly,lz),p33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)p11
  close(10)
  fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)p22
  close(10)
  fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)p33
  close(10)
  fpath='./out/p'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)p12
  close(10)
  write(*,*) 'after reading data files'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)
  p11=p11*g
  p22=p22*g 
  p33=p33*g
  p12=p12*g

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wx(ii,jj,kk)=eye*(ky(jj)*p33(ii,jj,kk)-kz(kk)*p22(ii,jj,kk))
    wy(ii,jj,kk)=eye*(kz(kk)*p11(ii,jj,kk)-kx(ii)*p33(ii,jj,kk))
    wz(ii,jj,kk)=eye*(kx(ii)*p22(ii,jj,kk)-ky(jj)*p11(ii,jj,kk))

    p11(ii,jj,kk) = - kx(ii) * kx(ii) * p12(ii,jj,kk)
    p22(ii,jj,kk) = - ky(jj) * ky(jj) * p12(ii,jj,kk)
    p33(ii,jj,kk) = - kz(kk) * kz(kk) * p12(ii,jj,kk)
    p13(ii,jj,kk) = - kx(ii) * kz(kk) * p12(ii,jj,kk)
    p23(ii,jj,kk) = - ky(jj) * kz(kk) * p12(ii,jj,kk)
    p12(ii,jj,kk) = - kx(ii) * ky(jj) * p12(ii,jj,kk)
  end do
  end do
  end do



  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

  g=real(wx)*real(wx)+real(wy)*real(wy)+real(wz)*real(wz)
  g=sqrt(g)

  wx=cmplx(real(wx)/g,aimag(wx))
  wy=cmplx(real(wy)/g,aimag(wy))
  wz=cmplx(real(wz)/g,aimag(wz))

  g=aimag(wx)*aimag(wx)+aimag(wy)*aimag(wy)+aimag(wz)*aimag(wz)
  g=sqrt(g)

  wx=cmplx(real(wx),aimag(wx)/g)
  wy=cmplx(real(wy),aimag(wy)/g)
  wz=cmplx(real(wz),aimag(wz)/g)


  popo=0.
  meanopo = 0.0_dp
  rmsopo = 0.0_dp
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx

    cc(1,1)=real(p11(ii,jj,kk))
    cc(1,2)=real(p12(ii,jj,kk))
    cc(1,3)=real(p13(ii,jj,kk))
    cc(2,1)=cc(1,2)
    cc(2,2)=real(p22(ii,jj,kk))
    cc(2,3)=real(p23(ii,jj,kk))
    cc(3,1)=cc(1,3)
    cc(3,2)=cc(2,3)
    cc(3,3)=real(p33(ii,jj,kk))

    ignore_me = -(cc(1,1) + cc(2,2) + cc(3,3))/3.
    cc(1,1) = cc(1,1) + ignore_me
    cc(2,2) = cc(2,2) + ignore_me
    cc(3,3) = cc(3,3) + ignore_me

    om(1) = real(wx(ii,jj,kk))
    om(2) = real(wy(ii,jj,kk))
    om(3) = real(wz(ii,jj,kk))

    ompom = 0.
    do mm = 1, 3
    do nn = 1, 3
      ompom = ompom + om(mm) * cc(mm,nn) * om(nn)
    end do
    end do

    ll = floor( (ompom + bnd)/binw ) + 1
    if (ll .ge. 1 .and. ll .le. npnt) popo(ll)=popo(ll)+1


    meanopo = meanopo + ompom
    rmsopo = rmsopo + ompom * ompom

    cc(1,1)=aimag(p11(ii,jj,kk))
    cc(1,2)=aimag(p12(ii,jj,kk))
    cc(1,3)=aimag(p13(ii,jj,kk))
    cc(2,1)=cc(1,2)
    cc(2,2)=aimag(p22(ii,jj,kk))
    cc(2,3)=aimag(p23(ii,jj,kk))
    cc(3,1)=cc(1,3)
    cc(3,2)=cc(2,3)
    cc(3,3)=aimag(p33(ii,jj,kk))

    ignore_me = -(cc(1,1) + cc(2,2) + cc(3,3))/3.
    cc(1,1) = cc(1,1) + ignore_me
    cc(2,2) = cc(2,2) + ignore_me
    cc(3,3) = cc(3,3) + ignore_me

    om(1) = aimag(wx(ii,jj,kk))
    om(2) = aimag(wy(ii,jj,kk))
    om(3) = aimag(wz(ii,jj,kk))

    ompom = 0.
    do mm = 1, 3
    do nn = 1, 3
      ompom = ompom + om(mm) * cc(mm,nn) * om(nn)
    end do
    end do

    ll = floor( (ompom + bnd)/binw ) + 1
    if (ll .ge. 1 .and. ll .le. npnt) popo(ll)=popo(ll)+1

    meanopo = meanopo + ompom
    rmsopo = rmsopo + ompom * ompom

  end do
  end do
  end do
  ignore_me=1./(nx*ny*nz)
  write(*,*) 'check popo:', sum(popo)*ignore_me
  popo=popo*ignore_me/binw

  meanopo = meanopo * ignore_me
  rmsopo = sqrt(rmsopo * ignore_me - meanopo * meanopo)

  write(*,*) 'mean opo', meanopo
  write(*,*) 'rms opo', rmsopo

  fnm =fnm(1:len_trim(fnm))//'-'//str(1:len_trim(str))//'dx'

  open(15, file='popo'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,*) (-bnd+(ii-.5)*binw - meanopo)/rmsopo, popo(ii)*rmsopo
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(g)

  call destroyplan3d

  write(*,*) 'done'

end program opo
