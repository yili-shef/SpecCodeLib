program alignot
  use mconstant
  use mfftwplan3d
  implicit none


  complex, allocatable, dimension(:,:,:) :: wx,wy,wz,tmp
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  real, allocatable, dimension(:,:,:) :: g, k2
  real, allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  real(sp), dimension(3) :: wvector
  integer :: ierr

  integer, parameter :: npnt=100
  real,    parameter :: binw=1./npnt
  real, dimension(npnt) :: palpha, pbeta, pgamma

  character(80) :: fnm,str,fpath
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real :: delta_c, ignore_me, const, kcut2
  real(dp) :: alpha, beta, gamm

  write(*,*) 
  write(*,'(''>>>> PDFs of the cos angles between omega and eigenvectors of the sgs stress <<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-ot-trun.x nx nfile ndel'
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

  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx
  kcut2 = ( pi / (delta_c/2) )**2

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(tmp(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz), k2(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)tmp
  close(10)
  fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)wz
  close(10)
  fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
    read(10)wy
  close(10)
  write(*,*) 'after reading data files'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)

  call rsgstauij(tmp,tmp,t11,g,nx,ny,nz) 
  call rsgstauij( wz, wz,t22,g,nx,ny,nz)
  call rsgstauij( wy, wy,t33,g,nx,ny,nz)
  t12=-(t11+t22+t33)/3.
  t11=t11+t12
  t22=t22+t12
  t33=t33+t12
  call rsgstauij(tmp, wz,t12,g,nx,ny,nz)
  call rsgstauij(tmp, wy,t13,g,nx,ny,nz)
  call rsgstauij( wz, wy,t23,g,nx,ny,nz)

  call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)
  
  where ( k2 .ge. kcut2 )
    t11 = 0.
    t12 = 0.
    t13 = 0.
    t22 = 0.
    t23 = 0.
    t33 = 0. 
  elsewhere
    t11 = t11 * const
    t12 = t12 * const
    t13 = t13 * const
    t22 = t22 * const
    t23 = t23 * const
    t33 = t33 * const 
  endwhere

  call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)

  ! tmp = ux; wz = uy; wy = uz
  ! staggered use to save memory
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wx(ii,jj,kk)=eye*(ky(jj)* wy(ii,jj,kk)-kz(kk)* wz(ii,jj,kk))*g(ii,jj,kk)
    wy(ii,jj,kk)=eye*(kz(kk)*tmp(ii,jj,kk)-kx(ii)* wy(ii,jj,kk))*g(ii,jj,kk)
    wz(ii,jj,kk)=eye*(kx(ii)* wz(ii,jj,kk)-ky(jj)*tmp(ii,jj,kk))*g(ii,jj,kk)
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  tmp=cmplx(real(wx)*real(wx)+real(wy)*real(wy)+real(wz)*real(wz), &
     aimag(wx)*aimag(wx)+aimag(wy)*aimag(wy)+aimag(wz)*aimag(wz))
  tmp=cmplx(sqrt(real(tmp)),sqrt(aimag(tmp)))

  wx=cmplx(real(wx)/real(tmp),aimag(wx)/aimag(tmp))
  wy=cmplx(real(wy)/real(tmp),aimag(wy)/aimag(tmp))
  wz=cmplx(real(wz)/real(tmp),aimag(wz)/aimag(tmp))

  palpha=0.
  pbeta =0.
  pgamma=0.
  do kk=1,nz
  do jj=1,ny
  do ii=1,nx

    if ( mod(ii,2) .eq. 1 ) then
      ll = (ii+1)/2
      cc(1,1)=real(t11(ll,jj,kk))
      cc(1,2)=real(t12(ll,jj,kk))
      cc(1,3)=real(t13(ll,jj,kk))
      cc(2,2)=real(t22(ll,jj,kk))
      cc(2,3)=real(t23(ll,jj,kk))
      cc(3,3)=real(t33(ll,jj,kk))

      wvector(1) = real(wx(ll,jj,kk))
      wvector(2) = real(wy(ll,jj,kk))
      wvector(3) = real(wz(ll,jj,kk))
    else
      ll = ii/2
      cc(1,1)=aimag(t11(ll,jj,kk))
      cc(1,2)=aimag(t12(ll,jj,kk))
      cc(1,3)=aimag(t13(ll,jj,kk))
      cc(2,2)=aimag(t22(ll,jj,kk))
      cc(2,3)=aimag(t23(ll,jj,kk))
      cc(3,3)=aimag(t33(ll,jj,kk))

      wvector(1) = aimag(wx(ll,jj,kk))
      wvector(2) = aimag(wy(ll,jj,kk))
      wvector(3) = aimag(wz(ll,jj,kk))
    end if
    cc(2,1)=cc(1,2)
    cc(3,1)=cc(1,3)
    cc(3,2)=cc(2,3)

    call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
    do ll=1,3
      evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
    end do
    alpha = wvector(1)*evectors(1,3)+wvector(2)*evectors(2,3)+wvector(3)*evectors(3,3)
    beta  = wvector(1)*evectors(1,2)+wvector(2)*evectors(2,2)+wvector(3)*evectors(3,2)
    gamm  = wvector(1)*evectors(1,1)+wvector(2)*evectors(2,1)+wvector(3)*evectors(3,1)

    alpha = abs(alpha)
    beta = abs(beta)
    gamm = abs(gamm)
    
    ll=floor((alpha)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
    ll=floor((beta)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
    ll=floor((gamm)/binw)+1
    if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1

  end do
  end do
  end do
  ignore_me=1./(nx*ny*nz)
  write(*,*) 'check palpha:', sum(palpha)*ignore_me
  write(*,*) 'check pbeta:',  sum(pbeta)*ignore_me
  write(*,*) 'check pgamma:', sum(pgamma)*ignore_me
  palpha=palpha*ignore_me/binw
  pbeta=pbeta*ignore_me/binw
  pgamma=pgamma*ignore_me/binw

  fnm =fnm(1:len_trim(fnm))//'-'//str(1:len_trim(str))//'dx'

  open(15, file='potalign-trun-'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,*) (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(tmp,t11,t12,t13,t22,t23,t33)
  deallocate(g,k2)

  call destroyplan3d

  write(*,*) 'done'

end program alignot
