program lagpolyvis
  use mconstant
  use mfftwplan3d
  use mwavenumber
  use mlagpoly
  implicit none

  !--------------------------------------
  ! Depends on order of Lagrangian interpolation
  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,lx1,lx,ly,lz,ndel, nprtcle, ip, norder
  real(sp) :: ignore_me, delta_c, const, dt, time 
  real(sp) :: relww, reldiffww, reldissww, meanww, meandiffww, meandissww
  real(sp) :: meanwwspec, meandiffwwspec, meandisswwspec
  real(sp) :: wwp, diffwwp, disswwp

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: ww, diffww, dissww
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, k2, g
  real(sp),    allocatable, dimension(:) :: wwspec, diffwwspec, disswwspec 

  character(80) :: str, flnm, prefix, str1, strorder
  
  real(sp), allocatable, dimension(:,:) :: bg
  real(sp) :: y(3),dxyz(3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 7) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./lagpoly-vis.x nx filelist ndel nprtcle dt prefix norder'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*) '                     nprtcle: number of particles'
          write(*,*) '                     dt: time step size'
          write(*,*) '                     prefix: prefix of coordinates data files'
          write(*,*) '                     norder: order of lagrange interpolation'
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

  ! norder
  call getarg(7, strorder)
  read(strorder, '(I20)') norder
  strorder = adjustl(strorder)

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

  allocate(ww(nx,ny,nz))
  allocate(diffww(nx,ny,nz), dissww(nx,ny,nz))
  allocate(g(lx1,ly,lz), k2(lx1,ly,lz) )
  allocate(kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  allocate(bg(norder,3))

  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle))
  allocate(wwspec(nprtcle), diffwwspec(nprtcle), disswwspec(nprtcle))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)

  open(25, file = 'lagpoly'//strorder(1:len_trim(strorder))//'-vis-'//flnm(1:len_trim(flnm))//'.dat')
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
    open(15,file='./out/specint-vis-'//prefix(1:len_trim(prefix))//str1(1:len_trim(str1)),&
            form='unformatted')
      read(15) wwspec, diffwwspec, disswwspec
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
    s11 = - s11 * k2 * rnu
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

    meanwwspec = 0.
    meandiffwwspec = 0.
    meandisswwspec = 0.
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
 
      call pre_interp(y,dxyz,bg,lhnode,norder)
      call value(ww,wwp,lhnode,bg,nx,ny,nz,norder)
      call value(diffww,diffwwp,lhnode,bg,nx,ny,nz,norder)
      call value(dissww,disswwp,lhnode,bg,nx,ny,nz,norder)

      relww = relww + (wwp - wwspec(ip))**2
      reldiffww = reldiffww + (diffwwp - diffwwspec(ip))**2
      reldissww = reldissww + (disswwp - disswwspec(ip))**2

      meanww = meanww + wwp
      meandiffww = meandiffww + diffwwp
      meandissww = meandissww + disswwp

      meanwwspec = meanwwspec + wwspec(ip)
      meandiffwwspec = meandiffwwspec + diffwwspec(ip)
      meandisswwspec = meandisswwspec + disswwspec(ip)

    end do
    !$omp end do
    !$omp end parallel

    meanww = meanww / nprtcle
    meandiffww = meandiffww / nprtcle
    meandissww = meandissww / nprtcle

    meanwwspec = meanwwspec / nprtcle
    meandiffwwspec = meandiffwwspec / nprtcle
    meandisswwspec = meandisswwspec / nprtcle

    relww = sqrt( relww / nprtcle ) / meanwwspec
    reldiffww = sqrt( reldiffww / nprtcle ) / abs(meandiffwwspec)
    reldissww = sqrt( reldissww / nprtcle ) / abs(meandisswwspec)

    write(25, '(20E15.4)') time, meanww, meandiffww, meandissww, &
                            relww, reldiffww, reldissww, &
                            meanwwspec, meandiffwwspec, meandisswwspec
    time = time + dt

  end do
  close(25)
  close(30)

  deallocate(ux, uy, uz, wx, wy, wz, kx, ky, kz, g, k2, xp, yp, zp)
  deallocate(s11,s12,s13,s22,s23,s33,ww,diffww,dissww)
  deallocate(wwspec, diffwwspec, disswwspec)

  call destroyplan3d

  write(*,*) 'Finished'


end program lagpolyvis
