program evpresshess
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none


  complex(sp), allocatable, dimension(:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33
  real(sp),    allocatable, dimension(:,:,:) :: kx,ky,kz,g

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer,  parameter :: npnt=100
  real(sp), parameter :: binw=2./npnt, ebnd=9.
  real(sp), dimension(npnt) :: palpha, pbeta, pgamma, pea, peb, peg

  character(80) :: fnm,str,fpath
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real(sp) :: delta_c, ignore_me, erangea, erangeb, erangeg, binwea, binweb, binweg
  real(dp) :: alpha, beta, gamm, rmsealpha,rmsebeta,rmsegamma,meanealpha,meanebeta,meanegamma

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of pressure Hessian <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./evpresshess.x nx nfile ndel'
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
  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
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

  wx = eye * (ky*p33 - kz*p22)
  wy = eye * (kz*p11 - kx*p33)
  wz = eye * (kx*p22 - ky*p11)

  p11 = - kx * kx * p12
  p22 = - ky * ky * p12
  p33 = - kz * kz * p12
  p13 = - kx * kz * p12
  p23 = - ky * kz * p12
  p12 = - kx * ky * p12


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


  palpha=0.
  pbeta=0.
  pgamma=0.
  meanealpha=0.0_dp
  meanebeta=0.0_dp
  meanegamma=0.0_dp
  rmsealpha=0.0_dp
  rmsebeta=0.0_dp
  rmsegamma=0.0_dp
  do kk=1,nz
  do jj=1,ny
  do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          cc(1,1) = real(p11(ll,jj,kk), sp)
          cc(1,2) = real(p12(ll,jj,kk), sp)
          cc(1,3) = real(p13(ll,jj,kk), sp)
          cc(2,2) = real(p22(ll,jj,kk), sp)
          cc(2,3) = real(p23(ll,jj,kk), sp)
          cc(3,3) = real(p33(ll,jj,kk), sp)
          ox = real(wx(ll,jj,kk), sp)
          oy = real(wy(ll,jj,kk), sp)
          oz = real(wz(ll,jj,kk), sp)
      else
          ll = ii/2
          cc(1,1) = aimag(p11(ll,jj,kk))
          cc(1,2) = aimag(p12(ll,jj,kk))
          cc(1,3) = aimag(p13(ll,jj,kk))
          cc(2,2) = aimag(p22(ll,jj,kk))
          cc(2,3) = aimag(p23(ll,jj,kk))
          cc(3,3) = aimag(p33(ll,jj,kk))
          ox = aimag(wx(ll,jj,kk))
          oy = aimag(wy(ll,jj,kk))
          oz = aimag(wz(ll,jj,kk))
      end if
      cc(2,1)=cc(1,2)
      cc(3,1)=cc(1,3)
      cc(3,2)=cc(2,3)

      ignore_me = -(cc(1,1) + cc(2,2) + cc(3,3))/3.
      cc(1,1) = cc(1,1) + ignore_me
      cc(2,2) = cc(2,2) + ignore_me
      cc(3,3) = cc(3,3) + ignore_me

      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
      alpha=ox*evectors(1,3)+oy*evectors(2,3)+oz*evectors(3,3)
      beta =ox*evectors(1,2)+oy*evectors(2,2)+oz*evectors(3,2)
      gamm =ox*evectors(1,1)+oy*evectors(2,1)+oz*evectors(3,1)
    
      ll=floor((alpha+1)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
      ll=floor((beta+1)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
      ll=floor((gamm+1)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1
 
      wx(ii,jj,kk)=cmplx(real(evalues(3)),aimag(wx(ii,jj,kk)))
      wy(ii,jj,kk)=cmplx(real(evalues(2)),aimag(wy(ii,jj,kk)))
      wz(ii,jj,kk)=cmplx(real(evalues(1)),aimag(wz(ii,jj,kk)))
 
      meanealpha=meanealpha+evalues(3)
      meanebeta=meanebeta+evalues(2)
      meanegamma=meanegamma+evalues(1)
 
      rmsealpha=rmsealpha+evalues(3)*evalues(3)
      rmsebeta=rmsebeta+evalues(2)*evalues(2)
      rmsegamma=rmsegamma+evalues(1)*evalues(1)

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

  meanealpha=meanealpha*ignore_me
  meanebeta=meanebeta*ignore_me
  meanegamma=meanegamma*ignore_me

  rmsealpha = sqrt(rmsealpha*ignore_me-meanealpha*meanealpha)
  rmsebeta = sqrt(rmsebeta*ignore_me-meanebeta*meanebeta)
  rmsegamma = sqrt(rmsegamma*ignore_me-meanegamma*meanebeta)


  ! for the alpha eigenvalue, lower bound is set to zero
  erangea = ebnd * real(rmsealpha)
  binwea = erangea/npnt

  ! for beta eigenvalue, lower is minus upper
  erangeb = ebnd * real(rmsebeta)
  binweb = 2.*erangeb/npnt

  ! for gamma eigenvalue, upper is zero, erangeg gives the abs of lower
  erangeg = ebnd * real(rmsegamma)
  binweg = erangeg/npnt

  pea=0.
  peb=0.
  peg=0.
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx
    ll=floor(real(wx(ii,jj,kk))/binwea)+1
    if (ll .ge. 1 .and. ll .le. npnt) pea(ll)=pea(ll)+1
    ll=floor((real(wy(ii,jj,kk))+erangeb)/binweb)+1
    if (ll .ge. 1 .and. ll .le. npnt) peb(ll)=peb(ll)+1
    ll=floor((real(wz(ii,jj,kk))+erangeg)/binweg)+1
    if (ll .ge. 1 .and. ll .le. npnt) peg(ll)=peg(ll)+1

    ll=floor(aimag(wx(ii,jj,kk))/binwea)+1
    if (ll .ge. 1 .and. ll .le. npnt) pea(ll)=pea(ll)+1
    ll=floor((aimag(wy(ii,jj,kk))+erangeb)/binweb)+1
    if (ll .ge. 1 .and. ll .le. npnt) peb(ll)=peb(ll)+1
    ll=floor((aimag(wz(ii,jj,kk))+erangeg)/binweg)+1
    if (ll .ge. 1 .and. ll .le. npnt) peg(ll)=peg(ll)+1

  end do
  end do
  end do
  write(*,*) 'check pea:', sum(pea)*ignore_me
  write(*,*) 'check peb:', sum(peb)*ignore_me
  write(*,*) 'check peg:', sum(peg)*ignore_me
  pea=pea/binwea*ignore_me
  peb=peb/binweb*ignore_me
  peg=peg/binweg*ignore_me

  fnm =fnm(1:len_trim(fnm))//'-'//str(1:len_trim(str))//'dx'

  open(15,file='prmtrpijev'//fnm(1:len_trim(fnm))//'.dat')
    write(15,*) 'mean of alpha', meanealpha
    write(15,*) 'mean of beta' , meanebeta
    write(15,*) 'mean of gamma', meanegamma
    write(15,*) 'rms of alpha:', rmsealpha
    write(15,*) 'rms of beta:' , rmsebeta
    write(15,*) 'rms of gamma:', rmsegamma
  close(15)
  
  open(15, file='popalign'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,*) -1+(ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  open(15, file='ppijevalue'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binwea, pea(ii), &
                                -erangeb+(ii-.5)*binweb, peb(ii), &
                                -erangeg+(ii-.5)*binweg, peg(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(g)

  call destroyplan3d

  write(*,*) 'done'

end program evpresshess
