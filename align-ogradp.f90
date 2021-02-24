program alignogradp
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none


  complex(sp), allocatable, dimension(:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33
  real(sp),    allocatable, dimension(:,:,:) :: kx,ky,kz,g,tmp

  integer,  parameter :: npnt=100
  real(sp), parameter :: binw=1._sp/npnt
  real(sp), dimension(npnt) :: pcos

  character(80) :: fnm,str,fpath,strdel
  integer  :: ndel,ii,jj,kk,ll,nx,ny,nz,lx,lx1,ly,lz, nfile
  real(sp) :: delta_c, ignore_me, ox, oy, oz, p1, p2, p3

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and pressure gradient <<<<<<'')')
  write(*,*)
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-ogradp.x nx fnamelist ndel'
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
  allocate(p11(lx1,ly,lz),p22(lx1,ly,lz),p33(lx1,ly,lz))
  allocate(g(lx1,ly,lz), tmp(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g = exp(-g*delta_c**2/24.)

  open(33, file = fnm(1:len_trim(fnm))//'.list')

  pcos=0._sp
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
 
    p11 = p11 * g
    p22 = p22 * g 
    p33 = p33 * g
 
    wx = eye * (ky*p33 - kz*p22)
    wy = eye * (kz*p11 - kx*p33)
    wz = eye * (kx*p22 - ky*p11)

    fpath='./out/p'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p11
    close(10)
    write(*,*) 'after reading data files'
    p11 = p11 * g
 
    p22 = eye * ky * p11
    p33 = eye * kz * p11
    p11 = eye * kx * p11
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
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
 
    tmp=real(p11)*real(p11)+real(p22)*real(p22)+real(p33)*real(p33)
    tmp=sqrt(tmp)
 
    p11=cmplx(real(p11)/tmp,aimag(p11))
    p22=cmplx(real(p22)/tmp,aimag(p22))
    p33=cmplx(real(p33)/tmp,aimag(p33))
 
    tmp=aimag(p11)*aimag(p11)+aimag(p22)*aimag(p22)+aimag(p33)*aimag(p33)
    tmp=sqrt(tmp)
 
    p11=cmplx(real(p11),aimag(p11)/tmp)
    p22=cmplx(real(p22),aimag(p22)/tmp)
    p33=cmplx(real(p33),aimag(p33)/tmp)
 
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
 
        if ( mod(ii, 2) .eq. 1 ) then
            ll = ii/2 + 1
            p1 = real(p11(ll,jj,kk), sp)
            p2 = real(p22(ll,jj,kk), sp)
            p3 = real(p33(ll,jj,kk), sp)
            ox = real(wx(ll,jj,kk), sp)
            oy = real(wy(ll,jj,kk), sp)
            oz = real(wz(ll,jj,kk), sp)
        else
            ll = ii/2
            p1 = aimag(p11(ll,jj,kk))
            p2 = aimag(p22(ll,jj,kk))
            p3 = aimag(p33(ll,jj,kk))
            ox = aimag(wx(ll,jj,kk))
            oy = aimag(wy(ll,jj,kk))
            oz = aimag(wz(ll,jj,kk))
        end if
        ignore_me = abs(p1 * ox + p2 * oy + p3 * oz)
      
        ll=floor(ignore_me/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) pcos(ll)=pcos(ll)+1
  
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(33)

  ignore_me = 1._sp / (nx*ny*nz*nfile)
  write(*,*) 'check pcos:', sum(pcos)*ignore_me

  pcos=pcos*ignore_me/binw


  fpath =fnm(1:len_trim(fnm))//'-'//strdel(1:len_trim(strdel))//'dx'
  open(15, file='ogradpalign-'//fpath(1:len_trim(fpath))//'.dat')
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5)*binw, pcos(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(p11,p22,p33)
  deallocate(g,tmp)

  call destroyplan3d

  write(*,*) 'done'

end program alignogradp
