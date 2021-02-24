program ptcovac
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, dt, time, meancurv, rmscurv, curvp, covalfcurv
  real(sp) :: wxp, wyp, wzp, meanalpha, rmsalpha, alphap, const

  complex(sp) :: a12, a13, a21, a23, a31, a32

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: wxr, wyr, wzr, alpha, sijr, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp
  real(sp),    allocatable, dimension(:)     :: alfp


  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3),idnty(3,3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-cov-alf-curv.x nx filelist ndel nprtcle dt prefix'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*) '                     nprtcle: number of particles'
          write(*,*) '                     dt: time step size'
          write(*,*) '                     prefix: prefix of coordinates data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! number of particles
  call getarg(4,flnm)
  read(flnm, '(I20)') nprtcle

  ! time step size
  call getarg(5,flnm)
  read(flnm, '(F20.10)') dt

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! prefix
  call getarg(6,prefix)

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz))
  allocate(wxr(nx,ny,nz), wyr(nx,ny,nz), wzr(nx,ny,nz))
  allocate(s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz))
  allocate(s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz))
  allocate(alpha(nx,ny,nz), sijr(nx,ny,nz) )
  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate(alfp(nprtcle))
  allocate(g(lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)


  open(25, file = 'cov-alf-curv-'//flnm(1:len_trim(flnm))//'.dat')
  open(30, file = flnm(1:len_trim(flnm))//'.list')

  time = 0.
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) ux
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uy
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uz
    close(15)

    open(15,file='./out/'//prefix(1:len_trim(prefix))//str1(1:len_trim(str1)),form='unformatted')
      read(15) xp, yp, zp
    close(15)
    write(*,*) 'finishing reading data'

    wx = ux * g
    wy = uy * g
    wz = uz * g

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a12 = eye * ky(jj) * wx(ii,jj,kk)
      a13 = eye * kz(kk) * wx(ii,jj,kk)
      a21 = eye * kx(ii) * wy(ii,jj,kk)
      a23 = eye * kz(kk) * wy(ii,jj,kk)
      a31 = eye * kx(ii) * wz(ii,jj,kk) 
      a32 = eye * ky(jj) * wz(ii,jj,kk)

      s11(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
      s22(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
      s33(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)

      s12(ii,jj,kk) = .5 * (a12 + a21)
      s13(ii,jj,kk) = .5 * (a13 + a31)
      s23(ii,jj,kk) = .5 * (a23 + a32)

      wx(ii,jj,kk) = a32 - a23
      wy(ii,jj,kk) = a13 - a31
      wz(ii,jj,kk) = a21 - a12

    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    
    wxr(1:nx:2,:,:)=real(wx(1:lx,:,:),sp); wxr(2:nx:2,:,:)=aimag(wx(1:lx,:,:))
    wyr(1:nx:2,:,:)=real(wy(1:lx,:,:),sp); wyr(2:nx:2,:,:)=aimag(wy(1:lx,:,:))
    wzr(1:nx:2,:,:)=real(wz(1:lx,:,:),sp); wzr(2:nx:2,:,:)=aimag(wz(1:lx,:,:))
 
    sijr(1:nx:2,:,:)=real(s11(1:lx,:,:),sp); sijr(2:nx:2,:,:)=aimag(s11(1:lx,:,:))
    alpha = wxr*sijr*wxr 

    sijr(1:nx:2,:,:)=real(s22(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s22(1:lx,:,:))
    alpha = alpha + wyr*sijr*wyr 

    sijr(1:nx:2,:,:)=real(s33(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s33(1:lx,:,:))
    alpha = alpha + wzr*sijr*wzr 

    sijr(1:nx:2,:,:)=real(s12(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s12(1:lx,:,:))
    alpha = alpha + 2 * wxr*sijr*wyr 

    sijr(1:nx:2,:,:)=real(s13(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s13(1:lx,:,:))
    alpha = alpha + 2 * wxr*sijr*wzr 

    sijr(1:nx:2,:,:)=real(s23(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s23(1:lx,:,:))
    alpha = alpha + 2 * wyr*sijr*wzr

    alpha = alpha/ (wxr*wxr + wyr*wyr + wzr*wzr)

    meanalpha = 0.
    rmsalpha = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(alpha,alphap,lhnode,bg,nx,ny,nz)

      alfp(ip) = alphap

      meanalpha = meanalpha + alphap
      rmsalpha = rmsalpha + alphap * alphap

    end do
    meanalpha = meanalpha / nprtcle
    rmsalpha = sqrt(rmsalpha / nprtcle - meanalpha*meanalpha)


    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx

      a12 = wx(ii,jj,kk) 
      a13 = wy(ii,jj,kk) 
      a23 = wz(ii,jj,kk) 

      wxp = sqrt( real(a12)**2 + real(a13)**2 + real(a23)**2 )
      wyp = sqrt( aimag(a12)**2 + aimag(a13)**2 + aimag(a23)**2 )

      wx(ii,jj,kk) = cmplx( real(a12) / wxp,  aimag(a12)/ wyp )  
      wy(ii,jj,kk) = cmplx( real(a13) / wxp,  aimag(a13)/ wyp )  
      wz(ii,jj,kk) = cmplx( real(a23) / wxp,  aimag(a23)/ wyp )  

      wxr(2*ii-1,jj,kk) = wxr(2*ii-1,jj,kk) / wxp
      wxr(2*ii  ,jj,kk) = wxr(2*ii  ,jj,kk) / wyp

      wyr(2*ii-1,jj,kk) = wyr(2*ii-1,jj,kk) / wxp
      wyr(2*ii  ,jj,kk) = wyr(2*ii  ,jj,kk) / wyp

      wzr(2*ii-1,jj,kk) = wzr(2*ii-1,jj,kk) / wxp
      wzr(2*ii  ,jj,kk) = wzr(2*ii  ,jj,kk) / wyp

    end do
    end do
    end do


    wx = wx * const; wy = wy * const; wz = wz * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      s11(ii,jj,kk) = eye * wx(ii,jj,kk) * ky(jj)
      s12(ii,jj,kk) = eye * wx(ii,jj,kk) * kz(kk)
      s13(ii,jj,kk) = eye * wy(ii,jj,kk) * kx(ii)
      s22(ii,jj,kk) = eye * wy(ii,jj,kk) * kz(kk)
      s23(ii,jj,kk) = eye * wz(ii,jj,kk) * kx(ii)
      s33(ii,jj,kk) = eye * wz(ii,jj,kk) * ky(jj)
      wx(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
      wy(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
      wz(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2
          wxp = wxr(ii,jj,kk) * real( wx (ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( s11(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( s12(ll,jj,kk) )
          wyp = wxr(ii,jj,kk) * real( s13(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( wy (ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( s22(ll,jj,kk) )
          wzp = wxr(ii,jj,kk) * real( s23(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( s33(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( wz (ll,jj,kk) )
      else
          ll = ii / 2
          wxp = wxr(ii,jj,kk) * aimag( wx (ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( s11(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( s12(ll,jj,kk) )
          wyp = wxr(ii,jj,kk) * aimag( s13(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( wy (ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( s22(ll,jj,kk) )
          wzp = wxr(ii,jj,kk) * aimag( s23(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( s33(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( wz (ll,jj,kk) )
      end if

      alpha(ii,jj,kk) = sqrt( wxp*wxp + wyp*wyp + wzp*wzp )
      
    end do
    end do
    end do

    meancurv = 0.
    rmscurv = 0.
    covalfcurv =0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(alpha,curvp,lhnode,bg,nx,ny,nz)

      meancurv = meancurv + curvp
      rmscurv = rmscurv + curvp * curvp
      covalfcurv = covalfcurv + curvp * alfp(ip)

    end do
    meancurv = meancurv/nprtcle
    rmscurv = sqrt( rmscurv/nprtcle - meancurv * meancurv)

    covalfcurv = covalfcurv/nprtcle - meancurv * meanalpha
    covalfcurv = covalfcurv / (rmsalpha * rmscurv)

    write(25, '(20E15.4)') time, covalfcurv, meancurv, rmscurv, meanalpha, rmsalpha
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(wx, wy, wz, wxr, wyr, wzr, kx, ky, kz, g, xp, yp, zp)
  deallocate(s11,s12,s13,s22,s23,s33,sijr,alpha,ux,uy,uz,alfp)

  call destroyplan3d

  write(*,*) 'Finished'


end program ptcovac      
