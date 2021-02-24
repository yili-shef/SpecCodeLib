program ptalpha
  use mconstant
  use mfftwplan3d
  implicit none

  real(sp), parameter :: rnu = 0.0025

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip

  real(sp) :: ignore_me, delta_c, dt, time, const
  real(sp) :: wxp, wyp, wzp, sxp, syp, szp, vxp, vyp, vzp, dsxp, dsyp, dszp
  real(sp) :: meanss, meanvs, meandiffs, meanph, meanphloc, meanphnonloc
  real(sp) :: ssp, vsp, diffsp, php, phlocp, phnonlocp

  complex(sp) :: a12, a13, a21, a23, a31, a32

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz, visx, visy, visz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  complex(sp), allocatable, dimension(:,:,:) :: p11, p12, p13, p22, p23, p33
  real(sp),    allocatable, dimension(:,:,:) :: diffs, ss, vs, phloc, ph, phnonloc, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-alpha.x nx filelist ndel nprtcle dt prefix'
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
  allocate( visx(lx1,ly,lz), visy(lx1,ly,lz), visz(lx1,ly,lz) )
  allocate( s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate( s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )
  allocate( ss(nx,ny,nz), vs(nx,ny,nz), diffs(nx,ny,nz) )
  allocate( phloc(nx,ny,nz), ph(nx,ny,nz), phnonloc(nx,ny,nz) )
  allocate( xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate( g(lx1,ly,lz), kx(lx1), ky(ly), kz(lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)


  open(25, file = 'alpha-'//flnm(1:len_trim(flnm))//'.dat')
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

      ignore_me = max( kx(ii) ** 2 + ky(jj) ** 2 + kz(kk) ** 2, mytiny ) 

      visx(ii,jj,kk) = - ignore_me * wx(ii,jj,kk) * rnu
      visy(ii,jj,kk) = - ignore_me * wy(ii,jj,kk) * rnu
      visz(ii,jj,kk) = - ignore_me * wz(ii,jj,kk) * rnu

      ! rnu * laplacian of sij
      p11(ii,jj,kk) = - ignore_me * s11(ii,jj,kk) * rnu
      p12(ii,jj,kk) = - ignore_me * s12(ii,jj,kk) * rnu
      p13(ii,jj,kk) = - ignore_me * s13(ii,jj,kk) * rnu
      p22(ii,jj,kk) = - ignore_me * s22(ii,jj,kk) * rnu
      p23(ii,jj,kk) = - ignore_me * s23(ii,jj,kk) * rnu
      p33(ii,jj,kk) = - ignore_me * s33(ii,jj,kk) * rnu

    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,visx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,visy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,visz,ignore_me)

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

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then 
            ll = (ii + 1) / 2

            wxp = real( wx(ll,jj,kk) ) 
            wyp = real( wy(ll,jj,kk) )
            wzp = real( wz(ll,jj,kk) )

            delta_c = max( sqrt( wxp**2 + wyp**2 + wzp**2 ), mytiny )
            wxp = wxp / delta_c 
            wyp = wyp / delta_c
            wzp = wzp / delta_c

            vxp = real( visx(ll,jj,kk) )
            vyp = real( visy(ll,jj,kk) )
            vzp = real( visz(ll,jj,kk) )

            vxp = vxp / delta_c
            vyp = vyp / delta_c
            vzp = vzp / delta_c

            sxp = wxp * real( s11(ll,jj,kk) ) &
                + wyp * real( s12(ll,jj,kk) ) &
                + wzp * real( s13(ll,jj,kk) ) 
            syp = wxp * real( s12(ll,jj,kk) ) &
                + wyp * real( s22(ll,jj,kk) ) &
                + wzp * real( s23(ll,jj,kk) ) 
            szp = wxp * real( s13(ll,jj,kk) ) &
                + wyp * real( s23(ll,jj,kk) ) &
                + wzp * real( s33(ll,jj,kk) ) 

            dsxp = wxp * real( p11(ll,jj,kk) ) &
                 + wyp * real( p12(ll,jj,kk) ) &
                 + wzp * real( p13(ll,jj,kk) ) 
            dsyp = wxp * real( p12(ll,jj,kk) ) &
                 + wyp * real( p22(ll,jj,kk) ) &
                 + wzp * real( p23(ll,jj,kk) ) 
            dszp = wxp * real( p13(ll,jj,kk) ) &
                 + wyp * real( p23(ll,jj,kk) ) &
                 + wzp * real( p33(ll,jj,kk) ) 

            ignore_me = real( s11(ll,jj,kk) )**2 + real( s22(ll,jj,kk) )**2 &
                      + real( s33(ll,jj,kk) )**2 + 2 * real( s12(ll,jj,kk) )**2 &
                      + 2 * real( s13(ll,jj,kk) )**2 + 2 * real( s23(ll,jj,kk) )**2

            ignore_me = delta_c**2 / 2 - ignore_me

            s11(ll,jj,kk) = cmplx( ignore_me, aimag( s11(ll,jj,kk) ) )

        else
            ll = ii / 2

            wxp = aimag( wx(ll,jj,kk) ) 
            wyp = aimag( wy(ll,jj,kk) )
            wzp = aimag( wz(ll,jj,kk) )

            delta_c = max( sqrt( wxp**2 + wyp**2 + wzp**2 ), mytiny )
            wxp = wxp / delta_c
            wyp = wyp / delta_c
            wzp = wzp / delta_c

            vxp = aimag( visx(ll,jj,kk) )
            vyp = aimag( visy(ll,jj,kk) )
            vzp = aimag( visz(ll,jj,kk) )

            vxp = vxp / delta_c
            vyp = vyp / delta_c
            vzp = vzp / delta_c

            sxp = wxp * aimag( s11(ll,jj,kk) ) &
                + wyp * aimag( s12(ll,jj,kk) ) &
                + wzp * aimag( s13(ll,jj,kk) ) 
            syp = wxp * aimag( s12(ll,jj,kk) ) &
                + wyp * aimag( s22(ll,jj,kk) ) &
                + wzp * aimag( s23(ll,jj,kk) ) 
            szp = wxp * aimag( s13(ll,jj,kk) ) &
                + wyp * aimag( s23(ll,jj,kk) ) &
                + wzp * aimag( s33(ll,jj,kk) ) 

            dsxp = wxp * aimag( p11(ll,jj,kk) ) &
                 + wyp * aimag( p12(ll,jj,kk) ) &
                 + wzp * aimag( p13(ll,jj,kk) ) 
            dsyp = wxp * aimag( p12(ll,jj,kk) ) &
                 + wyp * aimag( p22(ll,jj,kk) ) &
                 + wzp * aimag( p23(ll,jj,kk) ) 
            dszp = wxp * aimag( p13(ll,jj,kk) ) &
                 + wyp * aimag( p23(ll,jj,kk) ) &
                 + wzp * aimag( p33(ll,jj,kk) ) 

            ignore_me = aimag( s11(ll,jj,kk) )**2 + aimag( s22(ll,jj,kk) )**2 &
                      + aimag( s33(ll,jj,kk) )**2 + 2 * aimag( s12(ll,jj,kk) )**2 &
                      + 2 * aimag( s13(ll,jj,kk) )**2 + 2 * aimag( s23(ll,jj,kk) )**2

            ignore_me = delta_c**2 / 2 - ignore_me
            
            s11(ll,jj,kk) = cmplx( real( s11(ll,jj,kk) ), ignore_me ) 

        end if

        ! Local pressure term
        phloc(ii,jj,kk) = -ignore_me / 3.

        ! ignore_me = alpha
        ignore_me = sxp * wxp + syp * wyp + szp * wzp 

        ! self amplification
        ss(ii,jj,kk) = sxp**2 + syp**2 + szp**2 - 2 * ignore_me**2

        ! delta_c = V_alpha 
        delta_c = vxp * wxp + vyp * wyp + vzp * wzp

        ! viscous titlting
        vs(ii,jj,kk) = 2. * (sxp * vxp + syp * vyp + szp * vzp - ignore_me * delta_c)

        ! diffusion of strain rate
        diffs(ii,jj,kk) = wxp * dsxp + wyp * dsyp + wzp * dszp 

    end do
    end do
    end do

    s11 = s11 * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,s11,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
        ignore_me = max( kx(ii) ** 2 + ky(jj) ** 2 + kz(kk) ** 2, mytiny ) 

        ! TODO check this formula
        p11(ii,jj,kk) = kx(ii) * kx(ii) * s11(ii,jj,kk) / ignore_me
        p22(ii,jj,kk) = ky(jj) * ky(jj) * s11(ii,jj,kk) / ignore_me
        p33(ii,jj,kk) = kz(kk) * kz(kk) * s11(ii,jj,kk) / ignore_me

        p12(ii,jj,kk) = kx(ii) * ky(jj) * s11(ii,jj,kk) / ignore_me
        p13(ii,jj,kk) = kx(ii) * kz(kk) * s11(ii,jj,kk) / ignore_me
        p23(ii,jj,kk) = ky(jj) * kz(kk) * s11(ii,jj,kk) / ignore_me

    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then 
            ll = (ii + 1) / 2

            wxp = real( wx(ll,jj,kk) )
            wyp = real( wy(ll,jj,kk) )
            wzp = real( wz(ll,jj,kk) )

            delta_c = max( sqrt( wxp**2 + wyp**2 + wzp**2 ), mytiny )

            wxp = wxp / delta_c
            wyp = wyp / delta_c
            wzp = wzp / delta_c

            ph(ii,jj,kk) = wxp * real( p11(ll,jj,kk) ) * wxp &
                         + wxp * real( p12(ll,jj,kk) ) * wyp &
                         + wxp * real( p13(ll,jj,kk) ) * wzp &
                         + wyp * real( p12(ll,jj,kk) ) * wxp &
                         + wyp * real( p22(ll,jj,kk) ) * wyp &
                         + wyp * real( p23(ll,jj,kk) ) * wzp &
                         + wzp * real( p13(ll,jj,kk) ) * wxp &
                         + wzp * real( p23(ll,jj,kk) ) * wyp &
                         + wzp * real( p33(ll,jj,kk) ) * wzp
        else
            ll = ii / 2

            wxp = aimag( wx(ll,jj,kk) )
            wyp = aimag( wy(ll,jj,kk) )
            wzp = aimag( wz(ll,jj,kk) )

            delta_c = max( sqrt( wxp**2 + wyp**2 + wzp**2 ), mytiny )

            wxp = wxp / delta_c
            wyp = wyp / delta_c
            wzp = wzp / delta_c

            ph(ii,jj,kk) = wxp * aimag( p11(ll,jj,kk) ) * wxp &
                         + wxp * aimag( p12(ll,jj,kk) ) * wyp &
                         + wxp * aimag( p13(ll,jj,kk) ) * wzp &
                         + wyp * aimag( p12(ll,jj,kk) ) * wxp &
                         + wyp * aimag( p22(ll,jj,kk) ) * wyp &
                         + wyp * aimag( p23(ll,jj,kk) ) * wzp &
                         + wzp * aimag( p13(ll,jj,kk) ) * wxp &
                         + wzp * aimag( p23(ll,jj,kk) ) * wyp &
                         + wzp * aimag( p33(ll,jj,kk) ) * wzp
        end if

        ph(ii,jj,kk) = - ph(ii,jj,kk)

    
    end do
    end do
    end do

    phnonloc = ph - phloc

    meanss = 0.
    meanvs = 0.
    meandiffs = 0.
    meanphloc = 0.
    meanphnonloc = 0.
    meanph = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(ss,ssp,lhnode,bg,nx,ny,nz)
      call value(vs,vsp,lhnode,bg,nx,ny,nz)
      call value(diffs,diffsp,lhnode,bg,nx,ny,nz)
      call value(phloc,phlocp,lhnode,bg,nx,ny,nz)
      call value(phnonloc,phnonlocp,lhnode,bg,nx,ny,nz)
      call value(ph,php,lhnode,bg,nx,ny,nz)

      meanss = meanss + ssp
      meanvs = meanvs + vsp
      meandiffs = meandiffs + diffsp
      meanphloc = meanphloc + phlocp
      meanphnonloc = meanphnonloc + phnonlocp
      meanph = meanph + php

    end do
    meanss = meanss / nprtcle 
    meanvs = meanvs / nprtcle
    meandiffs = meandiffs / nprtcle
    meanphloc = meanphloc / nprtcle
    meanphnonloc = meanphnonloc / nprtcle
    meanph = meanph / nprtcle

    write(25, '(20E15.4)') time, meanss, meanvs, meandiffs, meanphloc, meanphnonloc, meanph
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate( ux, uy, uz, wx, wy, wz, kx, ky, kz, g, xp, yp, zp )
  deallocate( s11,s12,s13,s22,s23,s33,visx,visy,visz )
  deallocate( p11,p12,p13,p22,p23,p33,diffs,ss,vs,ph,phloc,phnonloc)

  call destroyplan3d

  write(*,*) 'Finished'


end program ptalpha      
