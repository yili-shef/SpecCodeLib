program ptvort
  use mconstant
  use mfftwplan3d
  implicit none

  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, const, dt, time 
  real(sp) :: meanww, meanalpha, meanalphaww, meandiffww, meandissww
  real(sp) :: wwp, alphap, alphawwp, diffwwp, disswwp

  complex(sp) :: a12, a13, a21, a23, a31, a32

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: ww, alpha, alphaww, diffww, dissww, g, k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-vort-new.x nx filelist ndel nprtcle dt prefix'
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

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz))
  allocate(s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz))
  allocate(s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz))

  allocate(ww(nx,ny,nz), alpha(nx,ny,nz), alphaww(nx,ny,nz))
  allocate(diffww(nx,ny,nz), dissww(nx,ny,nz))
  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle))
  allocate(g(lx1,ly,lz), k2(lx1,ly,lz), kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)


  open(25, file = 'ww-new-'//flnm(1:len_trim(flnm))//'.dat')
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

    ww(1:nx:2,:,:) = real( wx(1:lx,:,:) ) * real( wx(1:lx,:,:) ) &
                   + real( wy(1:lx,:,:) ) * real( wy(1:lx,:,:) ) &
                   + real( wz(1:lx,:,:) ) * real( wz(1:lx,:,:) ) 
    ww(2:nx:2,:,:) = aimag( wx(1:lx,:,:) ) * aimag( wx(1:lx,:,:) ) &
                   + aimag( wy(1:lx,:,:) ) * aimag( wy(1:lx,:,:) ) &
                   + aimag( wz(1:lx,:,:) ) * aimag( wz(1:lx,:,:) ) 

    alphaww(1:nx:2,:,:) = real( wx(1:lx,:,:) ) * real( s11(1:lx,:,:) ) * real( wx(1:lx,:,:) ) &
                        + real( wx(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( wy(1:lx,:,:) ) &
                        + real( wx(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( wz(1:lx,:,:) ) &
                        + real( wy(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( wx(1:lx,:,:) ) &
                        + real( wy(1:lx,:,:) ) * real( s22(1:lx,:,:) ) * real( wy(1:lx,:,:) ) &
                        + real( wy(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( wz(1:lx,:,:) ) &
                        + real( wz(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( wx(1:lx,:,:) ) &
                        + real( wz(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( wy(1:lx,:,:) ) &
                        + real( wz(1:lx,:,:) ) * real( s33(1:lx,:,:) ) * real( wz(1:lx,:,:) ) 

    alphaww(2:nx:2,:,:) = aimag( wx(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) * aimag( wx(1:lx,:,:) ) &
                        + aimag( wx(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( wy(1:lx,:,:) ) &
                        + aimag( wx(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( wz(1:lx,:,:) ) &
                        + aimag( wy(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( wx(1:lx,:,:) ) &
                        + aimag( wy(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) * aimag( wy(1:lx,:,:) ) &
                        + aimag( wy(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( wz(1:lx,:,:) ) &
                        + aimag( wz(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( wx(1:lx,:,:) ) &
                        + aimag( wz(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( wy(1:lx,:,:) ) &
                        + aimag( wz(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) * aimag( wz(1:lx,:,:) ) 

    alpha = alphaww / max(ww, mytiny)

    alphaww = 2 * alphaww 

    s11(1:lx,:,:) = cmplx( ww(1:nx:2,:,:), ww(2:nx:2,:,:) ) * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,s11,ignore_me)
    s11 = - s11 * k2 * rnu
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)

    diffww(1:nx:2,:,:) =  real( s11(1:lx,:,:) )
    diffww(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) 

    wx = wx * const; wy = wy * const; wz = wz * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
        s11(ii,jj,kk) = eye * ky(jj) * wx(ii,jj,kk) ! G12
        s12(ii,jj,kk) = eye * kz(kk) * wx(ii,jj,kk) ! G13
        s13(ii,jj,kk) = eye * kx(ii) * wy(ii,jj,kk) ! G21
        s22(ii,jj,kk) = eye * kz(kk) * wy(ii,jj,kk) ! G23
        s23(ii,jj,kk) = eye * kx(ii) * wz(ii,jj,kk) ! G31
        s33(ii,jj,kk) = eye * ky(jj) * wz(ii,jj,kk) ! G32
        wx(ii,jj,kk)  = eye * kx(ii) * wx(ii,jj,kk) ! G11
        wy(ii,jj,kk)  = eye * ky(jj) * wy(ii,jj,kk) ! G22
        wz(ii,jj,kk)  = eye * kz(kk) * wz(ii,jj,kk) ! G33
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

    dissww(1:nx:2,:,:) = real(  wx(1:lx,:,:) ) * real(  wx(1:lx,:,:) ) &
                       + real(  wy(1:lx,:,:) ) * real(  wy(1:lx,:,:) ) &
                       + real(  wz(1:lx,:,:) ) * real(  wz(1:lx,:,:) ) &
                       + real( s11(1:lx,:,:) ) * real( s11(1:lx,:,:) ) &
                       + real( s12(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &
                       + real( s13(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &
                       + real( s22(1:lx,:,:) ) * real( s22(1:lx,:,:) ) &
                       + real( s23(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &
                       + real( s33(1:lx,:,:) ) * real( s33(1:lx,:,:) ) 

    dissww(2:nx:2,:,:) = aimag(  wx(1:lx,:,:) ) * aimag(  wx(1:lx,:,:) ) &
                       + aimag(  wy(1:lx,:,:) ) * aimag(  wy(1:lx,:,:) ) &
                       + aimag(  wz(1:lx,:,:) ) * aimag(  wz(1:lx,:,:) ) &
                       + aimag( s11(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) &
                       + aimag( s12(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &
                       + aimag( s13(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &
                       + aimag( s22(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) &
                       + aimag( s23(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &
                       + aimag( s33(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) 

    dissww = -2 * rnu * dissww

    meanww = 0.
    meanalpha = 0.
    meanalphaww = 0.
    meandiffww = 0.
    meandissww = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(ww,wwp,lhnode,bg,nx,ny,nz)
      call value(alpha,alphap,lhnode,bg,nx,ny,nz)
      call value(alphaww,alphawwp,lhnode,bg,nx,ny,nz)
      call value(diffww,diffwwp,lhnode,bg,nx,ny,nz)
      call value(dissww,disswwp,lhnode,bg,nx,ny,nz)

      meanww = meanww + wwp
      meanalpha = meanalpha + alphap
      meanalphaww = meanalphaww + alphawwp
      meandiffww = meandiffww + diffwwp
      meandissww = meandissww + disswwp

    end do
    meanww = meanww / nprtcle 
    meanalpha = meanalpha / nprtcle
    meanalphaww = meanalphaww / nprtcle
    meandiffww = meandiffww / nprtcle
    meandissww = meandissww / nprtcle

    write(25, '(10E15.4)') time, meanww, meanalpha, meanalphaww, meandiffww, meandissww
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(ux, uy, uz, wx, wy, wz, kx, ky, kz, g, k2, xp, yp, zp)
  deallocate(s11,s12,s13,s22,s23,s33,ww,alpha,alphaww,diffww,dissww)

  call destroyplan3d

  write(*,*) 'Finished'


end program ptvort      
