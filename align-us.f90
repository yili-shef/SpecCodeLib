program alignos
  use mconstant
  use mfftwplan3d
  implicit none


  complex, allocatable, dimension(:,:,:) :: ux,uy,uz
  complex, allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real, allocatable, dimension(:,:,:) :: g
  real, allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer, parameter :: npnt=100
  real,    parameter :: binw=2./npnt
  real, dimension(npnt) :: palpha, pbeta, pgamma

  character(80) :: fnm,str,fpath
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real :: delta_c, ignore_me
  real(dp) :: alpha, beta, gamm

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between u and eigenvectors of strain rate <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-us.x nx nfile ndel'
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

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)ux
  close(10)
  fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uy
  close(10)
  fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uz
  close(10)
  write(*,*) 'after reading data files'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    s11(ii,jj,kk)=eye*kx(ii)*ux(ii,jj,kk)*g(ii,jj,kk)
    s22(ii,jj,kk)=eye*ky(jj)*uy(ii,jj,kk)*g(ii,jj,kk)
    s33(ii,jj,kk)=eye*kz(kk)*uz(ii,jj,kk)*g(ii,jj,kk)
    s12(ii,jj,kk)=.5*eye*(kx(ii)*uy(ii,jj,kk)+ky(jj)*ux(ii,jj,kk))*g(ii,jj,kk)
    s13(ii,jj,kk)=.5*eye*(kx(ii)*uz(ii,jj,kk)+kz(kk)*ux(ii,jj,kk))*g(ii,jj,kk)
    s23(ii,jj,kk)=.5*eye*(ky(jj)*uz(ii,jj,kk)+kz(kk)*uy(ii,jj,kk))*g(ii,jj,kk)
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)

  g=real(ux)*real(ux)+real(uy)*real(uy)+real(uz)*real(uz)
  g=sqrt(g)

  ux=cmplx(real(ux)/g,aimag(ux))
  uy=cmplx(real(uy)/g,aimag(uy))
  uz=cmplx(real(uz)/g,aimag(uz))
  
  g=aimag(ux)*aimag(ux)+aimag(uy)*aimag(uy)+aimag(uz)*aimag(uz)
  g=sqrt(g)
  ux=cmplx(real(ux),aimag(ux)/g)
  uy=cmplx(real(uy),aimag(uy)/g)
  uz=cmplx(real(uz),aimag(uz)/g)

  palpha=0.
  pbeta=0.
  pgamma=0.
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx

    cc(1,1)=real(s11(ii,jj,kk))
    cc(1,2)=real(s12(ii,jj,kk))
    cc(1,3)=real(s13(ii,jj,kk))
    cc(2,1)=cc(1,2)
    cc(2,2)=real(s22(ii,jj,kk))
    cc(2,3)=real(s23(ii,jj,kk))
    cc(3,1)=cc(1,3)
    cc(3,2)=cc(2,3)
    cc(3,3)=real(s33(ii,jj,kk))

    call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
    do ll=1,3
      evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
    end do
    alpha = real(ux(ii,jj,kk))*evectors(1,3)+real(uy(ii,jj,kk))*evectors(2,3)+real(uz(ii,jj,kk))*evectors(3,3)
    beta  = real(ux(ii,jj,kk))*evectors(1,2)+real(uy(ii,jj,kk))*evectors(2,2)+real(uz(ii,jj,kk))*evectors(3,2)
    gamm  = real(ux(ii,jj,kk))*evectors(1,1)+real(uy(ii,jj,kk))*evectors(2,1)+real(uz(ii,jj,kk))*evectors(3,1)
    
    ll=floor((alpha+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
    ll=floor((beta+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
    ll=floor((gamm+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1

    cc(1,1)=aimag(s11(ii,jj,kk))
    cc(1,2)=aimag(s12(ii,jj,kk))
    cc(1,3)=aimag(s13(ii,jj,kk))
    cc(2,1)=cc(1,2)
    cc(2,2)=aimag(s22(ii,jj,kk))
    cc(2,3)=aimag(s23(ii,jj,kk))
    cc(3,1)=cc(1,3)
    cc(3,2)=cc(2,3)
    cc(3,3)=aimag(s33(ii,jj,kk))

    call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
    do ll=1,3
      evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
    end do
    alpha = aimag(ux(ii,jj,kk))*evectors(1,3)+aimag(uy(ii,jj,kk))*evectors(2,3)+aimag(uz(ii,jj,kk))*evectors(3,3)
    beta  = aimag(ux(ii,jj,kk))*evectors(1,2)+aimag(uy(ii,jj,kk))*evectors(2,2)+aimag(uz(ii,jj,kk))*evectors(3,2)
    gamm  = aimag(ux(ii,jj,kk))*evectors(1,1)+aimag(uy(ii,jj,kk))*evectors(2,1)+aimag(uz(ii,jj,kk))*evectors(3,1)
    
    ll=floor((alpha+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
    ll=floor((beta+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
    ll=floor((gamm+1)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1

  end do
  end do
  end do
  write(*,*) 'check palpha:', sum(palpha)/(nx*ny*nz)
  write(*,*) 'check pbeta:',  sum(pbeta)/(nx*ny*nz)
  write(*,*) 'check pgamma:', sum(pgamma)/(nx*ny*nz)
  palpha=palpha/(nx*ny*nz)/binw
  pbeta=pbeta/(nx*ny*nz)/binw
  pgamma=pgamma/(nx*ny*nz)/binw

  fnm =fnm(1:len_trim(fnm))//'-'//str(1:len_trim(str))//'dx'
  open(15, file='pusalign'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,*) -1+(ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  deallocate(kx,ky,kz)
  deallocate(ux,uy,uz,s11,s12,s13,s22,s23,s33)
  deallocate(g)

  call destroyplan3d

  write(*,*) 'done'

end program alignos
