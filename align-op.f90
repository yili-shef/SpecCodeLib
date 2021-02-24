program alignop
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none


  complex(sp), allocatable, dimension(:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33
  real(sp),    allocatable, dimension(:,:,:) :: kx,ky,kz,g,tmp

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer,  parameter :: npnt=100
  real(sp), parameter :: binw=1./npnt
  real(sp), dimension(npnt) :: palpha, pbeta, pgamma

  character(80) :: fnm,str,fpath,strdel
  integer  :: ndel,ii,jj,kk,ll,nx,ny,nz,lx,lx1,ly,lz, nfile
  real(sp) :: delta_c, ignore_me, ox, oy, oz
  real(dp) :: alpha, beta, gamm

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of pressure Hessian <<<<<<'')')
  write(*,*)
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-op.x nx fnamelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        fnamelist: name of list file'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter parameter
  call getarg(3,strdel)
  read(strdel, '(I20)') ndel
  strdel = adjustl(strdel)

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
  allocate(g(lx1,ly,lz), tmp(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g = exp(-g*delta_c**2/24.)

  open(33, file = fnm(1:len_trim(fnm))//'.list')

  palpha=0.
  pbeta=0.
  pgamma=0.
  nfile = 0
  do while ( .not. eof(33) )

    read(33,*) str
    str = adjustl(str)
    write(*,*) 'reading data', str

    fpath='./out/ux'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p11
    close(10)
    fpath='./out/uy'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p22
    close(10)
    fpath='./out/uz'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p33
    close(10)
    fpath='./out/p'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p12
    close(10)
    write(*,*) 'after reading data files'
 
    p11 = p11 * g
    p22 = p22 * g 
    p33 = p33 * g
    p12 = p12 * g
 
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
 
    tmp=real(wx)*real(wx)+real(wy)*real(wy)+real(wz)*real(wz)
    tmp=sqrt(tmp)
 
    wx=cmplx(real(wx)/tmp,aimag(wx))
    wy=cmplx(real(wy)/tmp,aimag(wy))
    wz=cmplx(real(wz)/tmp,aimag(wz))
 
    tmp=aimag(wx)*aimag(wx)+aimag(wy)*aimag(wy)+aimag(wz)*aimag(wz)
    tmp=sqrt(tmp)
 
    wx=cmplx(real(wx),aimag(wx)/tmp)
    wy=cmplx(real(wy),aimag(wy)/tmp)
    wz=cmplx(real(wz),aimag(wz)/tmp)
 
 
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
        alpha = abs( ox*evectors(1,3)+oy*evectors(2,3)+oz*evectors(3,3) )
        beta  = abs( ox*evectors(1,2)+oy*evectors(2,2)+oz*evectors(3,2) )
        gamm  = abs( ox*evectors(1,1)+oy*evectors(2,1)+oz*evectors(3,1) )
      
        ll=floor(alpha/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
        ll=floor(beta/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
        ll=floor(gamm/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1
  
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(33)

  ignore_me = 1._sp / (nx*ny*nz*nfile)
  write(*,*) 'check palpha:', sum(palpha)*ignore_me
  write(*,*) 'check pbeta: ', sum(pbeta )*ignore_me
  write(*,*) 'check pgamma:', sum(pgamma)*ignore_me

  palpha=palpha*ignore_me/binw
  pbeta=pbeta*ignore_me/binw
  pgamma=pgamma*ignore_me/binw


  fpath =fnm(1:len_trim(fnm))//'-'//strdel(1:len_trim(strdel))//'dx'
  open(15, file='opalign-'//fpath(1:len_trim(fpath))//'.dat')
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(g,tmp)

  call destroyplan3d

  write(*,*) 'done'

end program alignop
