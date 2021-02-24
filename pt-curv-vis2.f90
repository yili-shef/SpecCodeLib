program ptcv
  use mconstant
  use mfftwplan3d
  implicit none

  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, dt, time, const 
  real(sp) :: meanvis, visp, c
  real(sp) :: wxp, wyp, wzp, cx, cy, cz, vfx, vfy, vfz, gvx, gvy, gvz

  complex(sp) :: wxc, wyc, wzc

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz, visx, visy, visz
  complex(sp), allocatable, dimension(:,:,:) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
  real(sp),    allocatable, dimension(:,:,:) :: g, k2, vis
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-curv-vis2.x nx filelist ndel nprtcle dt prefix'
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

  allocate( ux (lx1,ly,lz), uy (lx1,ly,lz), uz (lx1,ly,lz) )
  allocate( wx (lx1,ly,lz), wy (lx1,ly,lz), wz (lx1,ly,lz) )
  allocate( visx (lx1,ly,lz), visy (lx1,ly,lz), visz (lx1,ly,lz) )
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz) )
  allocate( g21(lx1,ly,lz), g22(lx1,ly,lz), g23(lx1,ly,lz) )
  allocate( g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  allocate( vis(nx,ny,nz), k2(lx1,ly,lz) )
  allocate( xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate( g(lx1,ly,lz),kx(lx1), ky(ly), kz(lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)

  open(25, file = 'curv-vis2-'//flnm(1:len_trim(flnm))//'.dat')
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
      wxc = eye * ky(jj) * wz(ii,jj,kk) - eye * kz(kk) * wy(ii,jj,kk) 
      wyc = eye * kz(kk) * wx(ii,jj,kk) - eye * kx(ii) * wz(ii,jj,kk)
      wzc = eye * kx(ii) * wy(ii,jj,kk) - eye * ky(jj) * wx(ii,jj,kk)

      wx(ii,jj,kk) = wxc
      wy(ii,jj,kk) = wyc
      wz(ii,jj,kk) = wzc

      g11(ii,jj,kk) = eye * wx(ii,jj,kk) * kx(ii)
      g12(ii,jj,kk) = eye * wx(ii,jj,kk) * ky(jj)
      g13(ii,jj,kk) = eye * wx(ii,jj,kk) * kz(kk)
      g21(ii,jj,kk) = eye * wy(ii,jj,kk) * kx(ii)
      g22(ii,jj,kk) = eye * wy(ii,jj,kk) * ky(jj)
      g23(ii,jj,kk) = eye * wy(ii,jj,kk) * kz(kk)
      g31(ii,jj,kk) = eye * wz(ii,jj,kk) * kx(ii)
      g32(ii,jj,kk) = eye * wz(ii,jj,kk) * ky(jj)
      g33(ii,jj,kk) = eye * wz(ii,jj,kk) * kz(kk)

      visx(ii,jj,kk) = -rnu * k2(ii,jj,kk) * wx(ii,jj,kk)
      visy(ii,jj,kk) = -rnu * k2(ii,jj,kk) * wy(ii,jj,kk)
      visz(ii,jj,kk) = -rnu * k2(ii,jj,kk) * wz(ii,jj,kk)

    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,visx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,visy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,visz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2

          ! omega_i 
          wxp = real( wx(ll,jj,kk) )
          wyp = real( wy(ll,jj,kk) )
          wzp = real( wz(ll,jj,kk) )

          cx =  real( g11(ll,jj,kk) ) * wxp &
             +  real( g12(ll,jj,kk) ) * wyp &
             +  real( g13(ll,jj,kk) ) * wzp 

          cy =  real( g21(ll,jj,kk) ) * wxp &
             +  real( g22(ll,jj,kk) ) * wyp &
             +  real( g23(ll,jj,kk) ) * wzp 

          cz =  real( g31(ll,jj,kk) ) * wxp &
             +  real( g32(ll,jj,kk) ) * wyp &
             +  real( g33(ll,jj,kk) ) * wzp 

          vfx = real( visx(ll,jj,kk) )
          vfy = real( visy(ll,jj,kk) )
          vfz = real( visz(ll,jj,kk) )

          gvx =  real( g11(ll,jj,kk) ) * vfx &
              +  real( g12(ll,jj,kk) ) * vfy &
              +  real( g13(ll,jj,kk) ) * vfz 

          gvy =  real( g21(ll,jj,kk) ) * vfx &
              +  real( g22(ll,jj,kk) ) * vfy &
              +  real( g23(ll,jj,kk) ) * vfz 

          gvz =  real( g31(ll,jj,kk) ) * vfx &
              +  real( g32(ll,jj,kk) ) * vfy &
              +  real( g33(ll,jj,kk) ) * vfz 
          
      else
          ll = ii / 2

          ! omega_i 
          wxp = aimag( wx(ll,jj,kk) )
          wyp = aimag( wy(ll,jj,kk) )
          wzp = aimag( wz(ll,jj,kk) )

          cx =  aimag( g11(ll,jj,kk) ) * wxp &
             +  aimag( g12(ll,jj,kk) ) * wyp &
             +  aimag( g13(ll,jj,kk) ) * wzp 

          cy =  aimag( g21(ll,jj,kk) ) * wxp &
             +  aimag( g22(ll,jj,kk) ) * wyp &
             +  aimag( g23(ll,jj,kk) ) * wzp 

          cz =  aimag( g31(ll,jj,kk) ) * wxp &
             +  aimag( g32(ll,jj,kk) ) * wyp &
             +  aimag( g33(ll,jj,kk) ) * wzp 

          vfx = aimag( visx(ll,jj,kk) )
          vfy = aimag( visy(ll,jj,kk) )
          vfz = aimag( visz(ll,jj,kk) )

          gvx =  aimag( g11(ll,jj,kk) ) * vfx &
              +  aimag( g12(ll,jj,kk) ) * vfy &
              +  aimag( g13(ll,jj,kk) ) * vfz 

          gvy =  aimag( g21(ll,jj,kk) ) * vfx &
              +  aimag( g22(ll,jj,kk) ) * vfy &
              +  aimag( g23(ll,jj,kk) ) * vfz 

          gvz =  aimag( g31(ll,jj,kk) ) * vfx &
              +  aimag( g32(ll,jj,kk) ) * vfy &
              +  aimag( g33(ll,jj,kk) ) * vfz 
          
      end if
      ! omega^2
      delta_c = wxp * wxp + wyp * wyp + wzp * wzp

      ignore_me = cx * wxp + cy * wyp + cz * wzp

      cx = cx - ignore_me * wxp / delta_c 
      cy = cy - ignore_me * wyp / delta_c
      cz = cz - ignore_me * wzp / delta_c

      cx = cx / delta_c
      cy = cy / delta_c
      cz = cz / delta_c

      c = sqrt( cx ** 2 + cy ** 2 + cz ** 2 )

      const = ( cx * gvx + cy * gvy + cz * gvz ) / c &
            - c * (vfx * wxp + vfy * wyp + vfz * wzp) 
      const = const / delta_c

      vis(ii,jj,kk) = const
      
    end do
    end do
    end do

    meanvis = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(vis,visp,lhnode,bg,nx,ny,nz)

      meanvis = meanvis + visp

    end do
    meanvis =  meanvis / nprtcle 

    write(25, '(10E15.4)') time, meanvis
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(wx, wy, wz, kx, ky, kz, g, xp, yp, zp)
  deallocate(g11,g12,g13,g21,g22,g23,g31,g32,g33,vis)
  deallocate(ux, uy, uz, visx, visy, visz)

  call destroyplan3d

  write(*,*) 'Finished'

end program ptcv      
