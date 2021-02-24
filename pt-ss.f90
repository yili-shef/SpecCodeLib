program ptss
  use mconstant
  use mfftwplan3d
  use mwavenumber
  use mlagpoly
  implicit none

  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,lx1,lx,ly,lz,ndel, nprtcle, ip, norder
  real(sp) :: ignore_me, delta_c, const, dt, time
  real(sp) :: ssp, sssp, spp, visssp
  real(sp) :: meanss, meansss, meansp, meanvis

  complex(4), allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),    allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  complex(sp), allocatable, dimension(:,:,:) :: p11, p12, p13, p22, p23, p33
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  real(sp),    allocatable, dimension(:,:,:) :: ssr, sssr, sijpij, visss
  real(sp),    allocatable, dimension(:,:,:) :: g, kx, ky, kz, k2

  character(80) :: str, flnm, prefix, str1, strnorder
  
  real(sp), allocatable, dimension(:,:) :: bg
  real(sp) :: y(3),dxyz(3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 7) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-ss.x nx filelist ndel nprtcle dt prefix norder'
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

  !========================== Reading parameters from the command line =============
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

  call getarg(7, strnorder)
  read(strnorder, '(I20)') norder
  strnorder = adjustl(strnorder)

  !================================================================================

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle) )

  allocate(s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz))
  allocate(s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz))
  allocate(p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz))
  allocate(p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )

  allocate(g(lx1,ly,lz), k2(lx1,ly,lz), kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  allocate( ssr(nx,ny,nz), visss(nx,ny,nz), sijpij(nx,ny,nz), sssr(nx,ny,nz) )

  allocate( bg(norder, 3) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)


  open(25, file = 'pt-ss-'//prefix(1:len_trim(prefix))//'.dat')
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

    s12 = eye * ( kx * s22 + ky * s11 ) * .5
    s13 = eye * ( kx * s33 + kz * s11 ) * .5
    s23 = eye * ( ky * s33 + kz * s22 ) * .5
    s11 = eye * kx * s11
    s22 = eye * ky * s22
    s33 = eye * kz * s33

    wx = eye * ( ky * s33 - kz * s22 )
    wy = eye * ( kz * s11 - kx * s33 )
    wz = eye * ( kx * s22 - ky * s11 )

    ! viscous diffusion for sij
    p11 = - rnu * k2 * s11
    p12 = - rnu * k2 * s12
    p13 = - rnu * k2 * s13
    p22 = - rnu * k2 * s22
    p23 = - rnu * k2 * s23
    p33 = - rnu * k2 * s33
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    

    visss(1:nx:2,:,:) = real( s11(1:lx,:,:) ) * real( p11(1:lx,:,:) ) &
                      + real( s22(1:lx,:,:) ) * real( p22(1:lx,:,:) ) &
                      + real( s33(1:lx,:,:) ) * real( p33(1:lx,:,:) ) &
                      + 2. * ( real( s12(1:lx,:,:) ) * real( p12(1:lx,:,:) )  &
                      +        real( s13(1:lx,:,:) ) * real( p13(1:lx,:,:) )  &
                      +        real( s23(1:lx,:,:) ) * real( p23(1:lx,:,:) ) )

    visss(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) * aimag( p11(1:lx,:,:) ) &
                      + aimag( s22(1:lx,:,:) ) * aimag( p22(1:lx,:,:) ) &
                      + aimag( s33(1:lx,:,:) ) * aimag( p33(1:lx,:,:) ) &
                      + 2. * ( aimag( s12(1:lx,:,:) ) * aimag( p12(1:lx,:,:) )  &
                      +        aimag( s13(1:lx,:,:) ) * aimag( p13(1:lx,:,:) )  &
                      +        aimag( s23(1:lx,:,:) ) * aimag( p23(1:lx,:,:) ) )


    ssr(1:nx:2,:,:) = real( s11(1:lx,:,:) ) ** 2 &
                    + real( s22(1:lx,:,:) ) ** 2 &
                    + real( s33(1:lx,:,:) ) ** 2 &
                    + 2 * ( real( s12(1:lx,:,:) ) ** 2 &
                    +       real( s13(1:lx,:,:) ) ** 2 &
                    +       real( s23(1:lx,:,:) ) ** 2 )
    
    ssr(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) ** 2 &
                    + aimag( s22(1:lx,:,:) ) ** 2 &
                    + aimag( s33(1:lx,:,:) ) ** 2 &
                    + 2 * ( aimag( s12(1:lx,:,:) ) ** 2 &
                    +       aimag( s13(1:lx,:,:) ) ** 2 &
                    +       aimag( s23(1:lx,:,:) ) ** 2 )


    p12(1:lx,:,:) = cmplx( ssr(1:nx:2,:,:), ssr(2:nx:2,:,:) ) * const
    p13 = cmplx( real( wx ) * real( wx ) + real( wy ) * real( wy ) + real( wz ) * real( wz ), &
            aimag( wx ) * aimag( wx ) + aimag( wy ) * aimag( wy ) + aimag( wz ) * aimag( wz ) )
    p13 = p13 * const / 2 

    wx = p13 - p12 ! source for pressure Laplacian
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)

    p11 = kx * kx * wx / k2
    p22 = ky * ky * wx / k2 
    p33 = kz * kz * wx / k2
    p12 = kx * ky * wx / k2
    p13 = kx * kz * wx / k2
    p23 = ky * kz * wx / k2

    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)

    wy = ( p11 + p22 + p33 ) / 3
    p11 = p11 - wy; p22 = p22 - wy; p33 = p33 - wy

    sijpij(1:nx:2,:,:) = real( s11(1:lx,:,:) ) * real( p11(1:lx,:,:) ) &
                       + real( s22(1:lx,:,:) ) * real( p22(1:lx,:,:) ) &
                       + real( s33(1:lx,:,:) ) * real( p33(1:lx,:,:) ) &
                       + 2. * ( real( s12(1:lx,:,:) ) * real( p12(1:lx,:,:) )  &
                       +        real( s13(1:lx,:,:) ) * real( p13(1:lx,:,:) )  &
                       +        real( s23(1:lx,:,:) ) * real( p23(1:lx,:,:) ) )

    sijpij(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) * aimag( p11(1:lx,:,:) ) &
                       + aimag( s22(1:lx,:,:) ) * aimag( p22(1:lx,:,:) ) &
                       + aimag( s33(1:lx,:,:) ) * aimag( p33(1:lx,:,:) ) &
                       + 2. * ( aimag( s12(1:lx,:,:) ) * aimag( p12(1:lx,:,:) )  &
                       +        aimag( s13(1:lx,:,:) ) * aimag( p13(1:lx,:,:) )  &
                       +        aimag( s23(1:lx,:,:) ) * aimag( p23(1:lx,:,:) ) )

    sijpij = - sijpij
    

    sssr(1:nx:2,:,:) = real( s11(1:lx,:,:) ) * real( s11(1:lx,:,:) ) * real( s11(1:lx,:,:) ) &
                     + real( s12(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s11(1:lx,:,:) ) &
                     + real( s13(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s11(1:lx,:,:) ) &

                     + real( s11(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &
                     + real( s12(1:lx,:,:) ) * real( s22(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &
                     + real( s13(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &

                     + real( s11(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &
                     + real( s12(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &
                     + real( s13(1:lx,:,:) ) * real( s33(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &

                     + real( s12(1:lx,:,:) ) * real( s11(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &
                     + real( s22(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s12(1:lx,:,:) ) &

                     + real( s12(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s22(1:lx,:,:) ) &
                     + real( s22(1:lx,:,:) ) * real( s22(1:lx,:,:) ) * real( s22(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s22(1:lx,:,:) ) &

                     + real( s12(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &
                     + real( s22(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s33(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &

                     + real( s13(1:lx,:,:) ) * real( s11(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &
                     + real( s33(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s13(1:lx,:,:) ) &

                     + real( s13(1:lx,:,:) ) * real( s12(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s22(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &
                     + real( s33(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s23(1:lx,:,:) ) &

                     + real( s13(1:lx,:,:) ) * real( s13(1:lx,:,:) ) * real( s33(1:lx,:,:) ) &
                     + real( s23(1:lx,:,:) ) * real( s23(1:lx,:,:) ) * real( s33(1:lx,:,:) ) &
                     + real( s33(1:lx,:,:) ) * real( s33(1:lx,:,:) ) * real( s33(1:lx,:,:) ) 

    sssr(2:nx:2,:,:) = aimag( s11(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) &
                     + aimag( s12(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) &
                     + aimag( s13(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) &

                     + aimag( s11(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &
                     + aimag( s12(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &
                     + aimag( s13(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &

                     + aimag( s11(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &
                     + aimag( s12(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &
                     + aimag( s13(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &

                     + aimag( s12(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &
                     + aimag( s22(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) &

                     + aimag( s12(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) &
                     + aimag( s22(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) &

                     + aimag( s12(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &
                     + aimag( s22(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &

                     + aimag( s13(1:lx,:,:) ) * aimag( s11(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &
                     + aimag( s33(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) &

                     + aimag( s13(1:lx,:,:) ) * aimag( s12(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s22(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &
                     + aimag( s33(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) &

                     + aimag( s13(1:lx,:,:) ) * aimag( s13(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) &
                     + aimag( s23(1:lx,:,:) ) * aimag( s23(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) &
                     + aimag( s33(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) * aimag( s33(1:lx,:,:) ) 

    sssr = - 2. * sssr

    meanss = 0.
    meansss = 0.
    meansp = 0.
    meanvis = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode, norder)
      call value(ssr,ssp,lhnode,bg,nx,ny,nz, norder)
      call value(sssr,sssp,lhnode,bg,nx,ny,nz, norder)
      call value(sijpij,spp,lhnode,bg,nx,ny,nz, norder)
      call value(visss,visssp,lhnode,bg,nx,ny,nz, norder)

      meanss = meanss + ssp
      meansss = meansss + sssp
      meansp = meansp + spp
      meanvis = meanvis + visssp

    end do
    meanss = meanss / nprtcle 
    meansss = meansss / nprtcle
    meansp = meansp / nprtcle
    meanvis = meanvis / nprtcle

    write(25, '(20E15.4)') time, meanss, meansss, meansp, meanvis
    time = time + dt

  end do
  close(25)
  close(30)

  deallocate(kx, ky, kz, k2, g, ux, uy, uz, xp, yp, zp)
  deallocate(s11, s12, s13, s22, s23, s33)
  deallocate(p11, p12, p13, p22, p23, p33)
  deallocate(wx, wy, wz, bg)
  deallocate(ssr, sssr, visss, sijpij)

  call destroyplan3d

  write(*,*) 'Finished'

end program ptss 
