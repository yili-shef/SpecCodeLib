program ptvortvis
  use mconstant
  use mfftwplan3d
  implicit none

  !--------------------------------------
  ! Depends on order of Lagrangian interpolation
  integer,  parameter :: norder = 6
  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, const, dt, time 
  real(sp) :: relww, reldiffww, reldissww, meanww, meandiffww, meandissww
  real(sp) :: wwspec, diffwwspec, disswwspec, wwp, diffwwp, disswwp

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: ww, diffww, dissww
  complex(sp), allocatable, dimension(:,:,:) :: wwc, diffwwc, disswwc
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, k2, g

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(norder,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./interp-vis.x nx filelist ndel nprtcle dt prefix'
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

  allocate(wwc(lx1,ly,lz), disswwc(lx1,ly,lz), diffwwc(lx1,ly,lz))

  allocate(ww(nx,ny,nz))
  allocate(diffww(nx,ny,nz), dissww(nx,ny,nz))
  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle))
  allocate(g(lx1,ly,lz), k2(lx1,ly,lz) )
  allocate(kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber3d(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)


  open(25, file = 'interp6-vis-'//flnm(1:len_trim(flnm))//'.dat')
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

    s11 = ux * g
    s22 = uy * g
    s33 = uz * g

    wx = eye * ( ky * s33 - kz * s22 )
    wy = eye * ( kz * s11 - kx * s33 )
    wz = eye * ( kx * s22 - ky * s11 )
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    ww(1:nx:2,:,:) = real( wx(1:lx,:,:) ) * real( wx(1:lx,:,:) ) &
                   + real( wy(1:lx,:,:) ) * real( wy(1:lx,:,:) ) &
                   + real( wz(1:lx,:,:) ) * real( wz(1:lx,:,:) ) 
    ww(2:nx:2,:,:) = aimag( wx(1:lx,:,:) ) * aimag( wx(1:lx,:,:) ) &
                   + aimag( wy(1:lx,:,:) ) * aimag( wy(1:lx,:,:) ) &
                   + aimag( wz(1:lx,:,:) ) * aimag( wz(1:lx,:,:) ) 

    s11(1:lx,:,:) = cmplx( ww(1:nx:2,:,:), ww(2:nx:2,:,:) ) * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,s11,ignore_me)
    wwc = s11
    s11 = - s11 * k2 * rnu
    diffwwc = s11
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)

    diffww(1:nx:2,:,:) =  real( s11(1:lx,:,:) )
    diffww(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) 

    wx = wx * const; wy = wy * const; wz = wz * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
    
    s11 = eye * ky * wx ! G12
    s12 = eye * kz * wx ! G13
    s13 = eye * kx * wy ! G21
    s22 = eye * kz * wy ! G23
    s23 = eye * kx * wz ! G31
    s33 = eye * ky * wz ! G32
    wx  = eye * kx * wx ! G11
    wy  = eye * ky * wy ! G22
    wz  = eye * kz * wz ! G33

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

    disswwc = cmplx( dissww(1:nx:2,:,:), dissww(2:nx:2,:,:) ) * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,disswwc,ignore_me)

    meanww = 0.
    meandiffww = 0.
    meandissww = 0.
    relww = 0.
    reldiffww = 0.
    reldissww = 0.
    !$omp parallel default(shared) &
    !$omp private(ip, y, bg, lhnode, wwp, diffwwp, disswwp, wwspec, diffwwspec, disswwspec) &
    !$omp reduction(+:meanww, meandiffww, meandissww, relww, reldiffww, reldissww)
    !$omp do
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(ww,wwp,lhnode,bg,nx,ny,nz)
      call value(diffww,diffwwp,lhnode,bg,nx,ny,nz)
      call value(dissww,disswwp,lhnode,bg,nx,ny,nz)

      call specint(wwc,wwspec,y(1),y(2),y(3),kx,ky,kz,lx1,ly,lz)
      call specint(diffwwc,diffwwspec,y(1),y(2),y(3),kx,ky,kz,lx1,ly,lz)
      call specint(disswwc,disswwspec,y(1),y(2),y(3),kx,ky,kz,lx1,ly,lz)

      relww = relww + (wwp - wwspec)**2
      reldiffww = reldiffww + (diffwwp - diffwwspec)**2
      reldissww = reldissww + (disswwp - disswwspec)**2

      meanww = meanww + wwspec
      meandiffww = meandiffww + diffwwspec
      meandissww = meandissww + disswwspec

    end do
    !$omp end do
    !$omp end parallel

    meanww = meanww / nprtcle
    meandiffww = meandiffww / nprtcle
    meandissww = meandissww / nprtcle

    relww = sqrt( relww / nprtcle ) / meanww
    reldiffww = sqrt( reldiffww / nprtcle ) / abs(meandiffww)
    reldissww = sqrt( reldissww / nprtcle ) / abs(meandissww)

    write(25, '(10E15.4)') time, meanww, meandiffww, meandissww, relww, reldiffww, reldissww
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(ux, uy, uz, wx, wy, wz, kx, ky, kz, g, k2, xp, yp, zp)
  deallocate(s11,s12,s13,s22,s23,s33,ww,diffww,dissww, wwc, diffwwc, disswwc)

  call destroyplan3d

  write(*,*) 'Finished'


end program ptvortvis
