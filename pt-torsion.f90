program pttorsion
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, dt, time, const 
  real(sp) :: meantsb, meanxigradsn, meantor, torsionp, c2
  real(sp) :: wxrc, wyrc, wzrc, wgx, wgy, wgz, wxp, wyp, wzp, bix, biy, biz

  complex(sp) :: wxc, wyc, wzc

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: wgradwx, wgradwy, wgradwz
  complex(sp), allocatable, dimension(:,:,:) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
  real(sp),    allocatable, dimension(:,:,:) :: g, torsion
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-torsion.x nx filelist ndel nprtcle dt prefix'
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

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( wgradwx(lx1,ly,lz), wgradwy(lx1,ly,lz), wgradwz(lx1,ly,lz))
  allocate( g11(lx1,ly,lz), g12(lx1,ly,lz), g13(lx1,ly,lz) )
  allocate( g21(lx1,ly,lz), g22(lx1,ly,lz), g23(lx1,ly,lz) )
  allocate( g31(lx1,ly,lz), g32(lx1,ly,lz), g33(lx1,ly,lz) )
  allocate( torsion(nx,ny,nz) )
  allocate( xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate( g(lx1,ly,lz),kx(lx1), ky(ly), kz(lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(25, file = 'torsion-'//flnm(1:len_trim(flnm))//'.dat')
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

      ! vorticity
      wxc = eye * ky(jj) * wz(ii,jj,kk) - eye * kz(kk) * wy(ii,jj,kk) 
      wyc = eye * kz(kk) * wx(ii,jj,kk) - eye * kx(ii) * wz(ii,jj,kk)
      wzc = eye * kx(ii) * wy(ii,jj,kk) - eye * ky(jj) * wx(ii,jj,kk)

      wx(ii,jj,kk) = wxc
      wy(ii,jj,kk) = wyc
      wz(ii,jj,kk) = wzc

      ! gradient of vorticity
      g11(ii,jj,kk) = eye * wx(ii,jj,kk) * kx(ii)
      g12(ii,jj,kk) = eye * wx(ii,jj,kk) * ky(jj) 
      g13(ii,jj,kk) = eye * wx(ii,jj,kk) * kz(kk) 
      g21(ii,jj,kk) = eye * wy(ii,jj,kk) * kx(ii)
      g22(ii,jj,kk) = eye * wy(ii,jj,kk) * ky(jj)
      g23(ii,jj,kk) = eye * wy(ii,jj,kk) * kz(kk) 
      g31(ii,jj,kk) = eye * wz(ii,jj,kk) * kx(ii)
      g32(ii,jj,kk) = eye * wz(ii,jj,kk) * ky(jj)
      g33(ii,jj,kk) = eye * wz(ii,jj,kk) * kz(kk)

    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx

        ! omega_i dot grad_i omega_j
        wxp = real( wx(ii,jj,kk) ) * real( g11(ii,jj,kk) ) &
            + real( wy(ii,jj,kk) ) * real( g12(ii,jj,kk) ) &
            + real( wz(ii,jj,kk) ) * real( g13(ii,jj,kk) )
        wyp = real( wx(ii,jj,kk) ) * real( g21(ii,jj,kk) ) &
            + real( wy(ii,jj,kk) ) * real( g22(ii,jj,kk) ) &
            + real( wz(ii,jj,kk) ) * real( g23(ii,jj,kk) )
        wzp = real( wx(ii,jj,kk) ) * real( g31(ii,jj,kk) ) &
            + real( wy(ii,jj,kk) ) * real( g32(ii,jj,kk) ) &
            + real( wz(ii,jj,kk) ) * real( g33(ii,jj,kk) )

        ! omega_i dot grad_i omega_j
        ignore_me = aimag( wx(ii,jj,kk) ) * aimag( g11(ii,jj,kk) ) &
                  + aimag( wy(ii,jj,kk) ) * aimag( g12(ii,jj,kk) ) &
                  + aimag( wz(ii,jj,kk) ) * aimag( g13(ii,jj,kk) )
        delta_c   = aimag( wx(ii,jj,kk) ) * aimag( g21(ii,jj,kk) ) &
                  + aimag( wy(ii,jj,kk) ) * aimag( g22(ii,jj,kk) ) &
                  + aimag( wz(ii,jj,kk) ) * aimag( g23(ii,jj,kk) )
        c2        = aimag( wx(ii,jj,kk) ) * aimag( g31(ii,jj,kk) ) &
                  + aimag( wy(ii,jj,kk) ) * aimag( g32(ii,jj,kk) ) &
                  + aimag( wz(ii,jj,kk) ) * aimag( g33(ii,jj,kk) )

        ! Vector wgradw = omega_i dot grad_i omega_j
        wgradwx(ii,jj,kk) = cmplx(wxp, ignore_me)
        wgradwy(ii,jj,kk) = cmplx(wyp, delta_c)
        wgradwz(ii,jj,kk) = cmplx(wzp, c2)
      
    end do
    end do
    end do

    wgradwx = wgradwx * const; wgradwy = wgradwy * const; wgradwz = wgradwz * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,wgradwx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wgradwy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wgradwz,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1

      ! gradient of wgradw
      g11(ii,jj,kk) = eye * wgradwx(ii,jj,kk) * kx(ii)
      g12(ii,jj,kk) = eye * wgradwx(ii,jj,kk) * ky(jj) 
      g13(ii,jj,kk) = eye * wgradwx(ii,jj,kk) * kz(kk) 
      g21(ii,jj,kk) = eye * wgradwy(ii,jj,kk) * kx(ii)
      g22(ii,jj,kk) = eye * wgradwy(ii,jj,kk) * ky(jj)
      g23(ii,jj,kk) = eye * wgradwy(ii,jj,kk) * kz(kk) 
      g31(ii,jj,kk) = eye * wgradwz(ii,jj,kk) * kx(ii)
      g32(ii,jj,kk) = eye * wgradwz(ii,jj,kk) * ky(jj)
      g33(ii,jj,kk) = eye * wgradwz(ii,jj,kk) * kz(kk)

    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwz,ignore_me)

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

        if ( mod(ii,2) .eq. 1 ) then

            ll = (ii + 1) / 2

            wxrc = real(wx(ll,jj,kk))
            wyrc = real(wy(ll,jj,kk))
            wzrc = real(wz(ll,jj,kk))

            wgx = real(wgradwx(ll,jj,kk))
            wgy = real(wgradwy(ll,jj,kk))
            wgz = real(wgradwz(ll,jj,kk))

            ! vector (omega_j dot grad_j) wgradw
            wxp = wxrc * real(g11(ll,jj,kk)) &
                + wyrc * real(g12(ll,jj,kk)) &
                + wzrc * real(g13(ll,jj,kk))
            wyp = wxrc * real(g21(ll,jj,kk)) &
                + wyrc * real(g22(ll,jj,kk)) &
                + wzrc * real(g23(ll,jj,kk))
            wzp = wxrc * real(g31(ll,jj,kk)) &
                + wyrc * real(g32(ll,jj,kk)) &
                + wzrc * real(g33(ll,jj,kk))

        else
            ll = ii / 2

            wxrc = aimag(wx(ll,jj,kk))
            wyrc = aimag(wy(ll,jj,kk))
            wzrc = aimag(wz(ll,jj,kk))

            wgx = aimag(wgradwx(ll,jj,kk))
            wgy = aimag(wgradwy(ll,jj,kk))
            wgz = aimag(wgradwz(ll,jj,kk))

            wxp = wxrc * aimag(g11(ll,jj,kk)) &
                + wyrc * aimag(g12(ll,jj,kk)) &
                + wzrc * aimag(g13(ll,jj,kk))
            wyp = wxrc * aimag(g21(ll,jj,kk)) &
                + wyrc * aimag(g22(ll,jj,kk)) &
                + wzrc * aimag(g23(ll,jj,kk))
            wzp = wxrc * aimag(g31(ll,jj,kk)) &
                + wyrc * aimag(g32(ll,jj,kk)) &
                + wzrc * aimag(g33(ll,jj,kk))

        end if

        ! omega^2
        delta_c = wxrc**2 + wyrc**2 + wzrc**2

        ! curvature squared
        c2 = ( wgx ** 2 + wgy ** 2 + wgz ** 2 ) / delta_c ** 2 &
           - ( wxrc * wgx + wyrc * wgy + wzrc * wgz ) ** 2 / delta_c ** 3

        ! omega cross wgradw dot (omega dot grad) wgradw
        ignore_me = wxp * ( wyrc * wgz - wzrc * wgy) &
                  + wyp * ( wzrc * wgx - wxrc * wgz) &
                  + wzp * ( wxrc * wgy - wyrc * wgx)
        torsion(ii,jj,kk) = ignore_me / c2 / delta_c**3

        
    end do
    end do
    end do

    meantor = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(torsion,torsionp,lhnode,bg,nx,ny,nz)

      meantor = meantor + torsionp

    end do
    meantor =  meantor / nprtcle 

    g21 = ux * g
    g31 = uy * g
    g32 = uz * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1

        ! strain rate tensor
        g11(ii,jj,kk) = eye * g21(ii,jj,kk) * kx(ii)
        g22(ii,jj,kk) = eye * g31(ii,jj,kk) * ky(jj)
        g33(ii,jj,kk) = eye * g32(ii,jj,kk) * kz(kk)
        g12(ii,jj,kk) = eye * ( g21(ii,jj,kk) * ky(jj) + g31(ii,jj,kk) * kx(ii) ) / 2
        g13(ii,jj,kk) = eye * ( g21(ii,jj,kk) * kz(kk) + g32(ii,jj,kk) * kx(ii) ) / 2
        g23(ii,jj,kk) = eye * ( g31(ii,jj,kk) * kz(kk) + g32(ii,jj,kk) * ky(jj) ) / 2
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
            ll = (ii + 1) / 2

            ! omega_i
            wxrc = real( wx(ll,jj,kk) )
            wyrc = real( wy(ll,jj,kk) )
            wzrc = real( wz(ll,jj,kk) )

            ! sij*omega_j
            wxp = wxrc * real( g11(ll,jj,kk) ) + wyrc * real( g12(ll,jj,kk) ) &
                + wzrc * real( g13(ll,jj,kk) )
            wyp = wxrc * real( g12(ll,jj,kk) ) + wyrc * real( g22(ll,jj,kk) ) &
                + wzrc * real( g23(ll,jj,kk) )
            wzp = wxrc * real( g13(ll,jj,kk) ) + wyrc * real( g23(ll,jj,kk) ) &
                + wzrc * real( g33(ll,jj,kk) )

            ! omega
            delta_c   = sqrt( wxrc * wxrc + wyrc * wyrc + wzrc * wzrc )

            ! wgradw_i
            wgx = real( wgradwx(ll,jj,kk) )
            wgy = real( wgradwy(ll,jj,kk) )
            wgz = real( wgradwz(ll,jj,kk) )

            ! wgradw dot omega
            ignore_me = wgx * wxrc + wgy * wyrc + wgz * wzrc

            ! vector c_i, |c| = curvature
            wgx = ( wgx - wxrc * ignore_me / delta_c**2 ) / delta_c**2
            wgy = ( wgy - wyrc * ignore_me / delta_c**2 ) / delta_c**2
            wgz = ( wgz - wzrc * ignore_me / delta_c**2 ) / delta_c**2

            ! curvature c
            ignore_me = sqrt( wgx * wgx + wgy * wgy + wgz * wgz )

            ! save for next calculation
            ! g21 = s_n * c * omega
            ! g31 = c * omega
            ! g32 = 1 / (c * omega)
            g21(ll,jj,kk) = cmplx( wxp * wgx + wyp * wgy + wzp * wgz, aimag( g21(ll,jj,kk) ) )
            g31(ll,jj,kk) = cmplx( delta_c * ignore_me, aimag( g31(ll,jj,kk) ) )
            g32(ll,jj,kk) = cmplx( 1. / ( delta_c * ignore_me + mytiny ), aimag( g32(ll,jj,kk) ) )

        else
            ll = ii / 2

            ! omega_i
            wxrc = aimag( wx(ll,jj,kk) )
            wyrc = aimag( wy(ll,jj,kk) )
            wzrc = aimag( wz(ll,jj,kk) )

            ! sij*omega_j
            wxp = wxrc * aimag( g11(ll,jj,kk) ) + wyrc * aimag( g12(ll,jj,kk) ) &
                + wzrc * aimag( g13(ll,jj,kk) )
            wyp = wxrc * aimag( g12(ll,jj,kk) ) + wyrc * aimag( g22(ll,jj,kk) ) &
                + wzrc * aimag( g23(ll,jj,kk) )
            wzp = wxrc * aimag( g13(ll,jj,kk) ) + wyrc * aimag( g23(ll,jj,kk) ) &
                + wzrc * aimag( g33(ll,jj,kk) )

            ! omega
            delta_c = sqrt( wxrc * wxrc + wyrc * wyrc + wzrc * wzrc )

            ! wgradw_i
            wgx = aimag( wgradwx(ll,jj,kk) )
            wgy = aimag( wgradwy(ll,jj,kk) )
            wgz = aimag( wgradwz(ll,jj,kk) )

            ! wgradw dot omega
            ignore_me = wgx * wxrc + wgy * wyrc + wgz * wzrc

            ! vector c_i, |c| = curvature
            wgx = ( wgx - wxrc * ignore_me / delta_c**2 ) / delta_c**2
            wgy = ( wgy - wyrc * ignore_me / delta_c**2 ) / delta_c**2
            wgz = ( wgz - wzrc * ignore_me / delta_c**2 ) / delta_c**2

            ! curvature c
            ignore_me = sqrt( wgx * wgx + wgy * wgy + wgz * wgz )

            ! g21 = s_n * c * omega = c_i*sij*omega_j
            ! g31 = c * omega
            ! g32 = 1 / (c * omega)
            g21(ll,jj,kk) = cmplx( real( g21(ll,jj,kk) ), wxp * wgx + wyp * wgy + wzp * wgz )
            g31(ll,jj,kk) = cmplx( real( g31(ll,jj,kk) ), delta_c * ignore_me )
            g32(ll,jj,kk) = cmplx( real( g32(ll,jj,kk) ), 1. / ( delta_c * ignore_me + mytiny ) )

        end if

        ! binormal bi_i = omega_i cross c_i / omega / c
        bix = ( wyrc * wgz - wzrc * wgy ) / ( delta_c * ignore_me )
        biy = ( wzrc * wgx - wxrc * wgz ) / ( delta_c * ignore_me )
        biz = ( wxrc * wgy - wyrc * wgx ) / ( delta_c * ignore_me )

        ! torsion = T * s_b
        torsion(ii,jj,kk) = torsion(ii,jj,kk) * ( bix * wxp + biy * wyp + biz * wzp)  &
                          / delta_c 


    end do
    end do
    end do


    meantsb = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(torsion,torsionp,lhnode,bg,nx,ny,nz)

      meantsb = meantsb + torsionp

    end do
    meantsb =  meantsb / nprtcle 

    ! back up sn * c * omega = c_i s_ij omega_j in real space
    g12 = g21

    ! g21 = s_n * c * omega
    ! g31 = c * omega
    g21 = g21 * const; g31 = g31 * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,g21,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,g31,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
        ! gradient of c_i*s_ij*omega_j
        wgradwx(ii,jj,kk) = eye * kx(ii) * g21(ii,jj,kk)
        wgradwy(ii,jj,kk) = eye * ky(jj) * g21(ii,jj,kk)
        wgradwz(ii,jj,kk) = eye * kz(kk) * g21(ii,jj,kk)

        ! gradient of c * omega
        g11(ii,jj,kk) = eye * kx(ii) * g31(ii,jj,kk)
        g22(ii,jj,kk) = eye * ky(jj) * g31(ii,jj,kk)
        g33(ii,jj,kk) = eye * kz(kk) * g31(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wgradwz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,g33,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
            ll = ( ii + 1 ) / 2

            ! omega_i
            wxrc = real( wx(ll,jj,kk) )
            wyrc = real( wy(ll,jj,kk) )
            wzrc = real( wz(ll,jj,kk) )

            ! gradient of c_i*s_ij*omega_j
            wgx = real( wgradwx(ll,jj,kk) )
            wgy = real( wgradwy(ll,jj,kk) )
            wgz = real( wgradwz(ll,jj,kk) )

            ! gradient of c * omega
            wxp = real( g11(ll,jj,kk) )
            wyp = real( g22(ll,jj,kk) )
            wzp = real( g33(ll,jj,kk) )

            ! 1/(c*omega)
            delta_c = real( g32(ll,jj,kk) ) 
            ! c_i*sij*omega_j
            c2 = real( g12(ll,jj,kk) )
        else 
            ll = ii / 2

            ! omega_i
            wxrc = aimag( wx(ll,jj,kk) )
            wyrc = aimag( wy(ll,jj,kk) )
            wzrc = aimag( wz(ll,jj,kk) )

            ! gradient of c_i*s_ij*omega_j
            wgx = aimag( wgradwx(ll,jj,kk) )
            wgy = aimag( wgradwy(ll,jj,kk) )
            wgz = aimag( wgradwz(ll,jj,kk) )

            ! gradient of c*omega
            wxp = aimag( g11(ll,jj,kk) )
            wyp = aimag( g22(ll,jj,kk) )
            wzp = aimag( g33(ll,jj,kk) )

            ! 1/(c*omega)
            delta_c = aimag( g32(ll,jj,kk) ) 
            ! c_i*sij*omega_j
            c2 = aimag( g12(ll,jj,kk) )
        end if

        ! omega
        ignore_me = sqrt( wxrc * wxrc + wyrc * wyrc + wzrc * wzrc )
        ! xi
        wxrc = wxrc / ignore_me
        wyrc = wyrc / ignore_me
        wzrc = wzrc / ignore_me

        ! xi dot grad s_n
        torsion(ii,jj,kk) = ( wxrc * wgx + wyrc * wgy + wzrc * wgz ) * delta_c &
                          - c2 * delta_c**2 * (wxrc * wxp + wyrc * wyp + wzrc * wzp)
    end do
    end do
    end do

    meanxigradsn = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(torsion,torsionp,lhnode,bg,nx,ny,nz)

      meanxigradsn = meanxigradsn + torsionp

    end do
    meanxigradsn =  meanxigradsn / nprtcle 


    write(25, '(10E15.4)') time, meantor, -meantsb, meanxigradsn
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(wx, wy, wz, kx, ky, kz, g, xp, yp, zp)
  deallocate(g11,g12,g13,g21,g22,g23,g31,g32,g33,torsion)
  deallocate(ux, uy, uz, wgradwx, wgradwy, wgradwz)

  call destroyplan3d

  write(*,*) 'Finished'

end program pttorsion      
