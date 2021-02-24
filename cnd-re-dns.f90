program cndredns
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  real(sp), parameter :: eps=0.1
  real(sp), parameter :: dt=0.001

  integer,  parameter :: ndt = 1, ndir = 8

  integer,  parameter :: nx=256,ny=nx,nz=nx
  integer,  parameter :: lx=nx/2,ly=ny,lz=nz,lx1=lx+1
  real(sp), parameter :: meshx=2._sp*pi/nx, meshy=2._sp*pi/ny, meshz=2._sp*pi/nz

  complex(sp), allocatable, dimension (:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension (:,:,:) :: uxn,uyn,uzn, uxnp, uynp, uznp, g
  real(sp),    allocatable, dimension (:,:,:) :: durdt_dns,dusdt_dns,dutdt_dns,durdt,dusdt,dutdt
  real(sp),    allocatable, dimension (:)     :: kx,ky,kz

  character(1) :: nullchr
  character(80) :: path,path1,dflnm1a,dflnm1b,dflnm1c,dflnm2a,dflnm2b,dflnm2c,cndel

  integer :: ii, jj, kk, ndim, nfile, ll, mm, nsets
  integer :: ndel, nxp, nyp, nzp, fnum1, fnum2, fskip

  real(sp) :: wxp, wyp, wzp, xpe, ype, zpe, uxi, uyi, uzi, uxe, uye, uze
  real(sp) :: del_uxp, del_uyp, del_uzp, del_urp, del_usp, del_utp
  real(sp) :: del_x, del_y, del_z, del_xp, del_yp, del_zp, xp, yp, zp
  real(sp) :: dl, ignore_me, alpha, del_us, del_ur, del_ut

  real(sp), dimension(3) :: rhat,shat,that, shatn, thatn

  integer,  parameter :: npntx=100
  real(sp), parameter :: width=14._sp
  real(8),  dimension(npntx) :: cndur, cndus, cndut, pdfur, pdfus, pdfut

  real(8) :: tt

  real(sp) :: mdusdt_dns0, rdusdt_dns0
  real(sp) :: maxmx, minmx, bnwx
  real(sp) :: aa, bb

  real(sp) :: ro, omega
  
  ll = iargc()
  if ( ll .ne. 5 ) then
      write(*,*) ' Wrong number of arguments. Usage:'
      write(*,*) 
      write(*,*) ' ./correldns.x ndel fnum1 fnum2 fskip ro'
      write(*,*)
      write(*,*) ' ro > 100 for non-rotating turbulence'
      write(*,*)
      write(*,*) ' Stoppted'
      stop
  end if

  call getarg(1, cndel)
  read(cndel, '(i20)') ndel
  cndel = adjustl(cndel)

  call getarg(2, path)
  read(path, '(i20)') fnum1

  call getarg(3, path)
  read(path, '(i20)') fnum2

  call getarg(4, path)
  read(path, '(i20)') fskip

  call getarg(5, path)
  read(path, '(f20.6)') ro

  if ( ro .ge. 100. ) then
      omega = 0.
  else
      omega = eps**(1./3.)*1.5**(2./3.)/ro ! 1.5 is the mean forcing wavenumber. see spectral_rot_dns_ft
  end if

  write(*,*) 'nx = ', nx, 'ny = ', ny, 'nz = ', nz
  write(*,*) 'ndel = ', ndel, 'fnum1 = ', fnum1, 'fnum2 = ', fnum2
  write(*,*) 'omega = ', omega

  call fftwplan3d(nx,ny,nz)
  nullchr=char(0)

  write(*,*) 'allocating'
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz),g(lx1,ly,lz))
  allocate(uxn(nx,ny,nz),uyn(nx,ny,nz),uzn(nx,ny,nz))
  allocate(uxnp(nx,ny,nz),uynp(nx,ny,nz),uznp(nx,ny,nz))
  allocate(durdt(nx,ny,nz),dusdt(nx,ny,nz),dutdt(nx,ny,nz))
  allocate(durdt_dns(nx,ny,nz),dusdt_dns(nx,ny,nz),dutdt_dns(nx,ny,nz))

  write(*,*) 'wavenumber'
  call wavenumber(kx,ky,kz,g,lx1,ly,lz)

  write(*,*) 'define filter and the displacement'
  dl=ndel*(2._sp*pi/real(nx,sp))
  g=exp(-g*dl**2/24._sp)


  tt = 1./(nx*ny*nz)
  nsets = 0

  cndur = 0.
  cndus = 0.
  cndut = 0.
  pdfur = 0.
  pdfus = 0.
  pdfut = 0.
  do nfile = fnum1, fnum2, fskip

    write(path,'(I5)') nfile
    write(path1,'(I5)') nfile+1
    path = adjustl(path)
    path1 = adjustl(path1)

    dflnm1a='./out/ux'//path(1:len_trim(path))//'.dat'//char(0)
    dflnm1b='./out/uy'//path(1:len_trim(path))//'.dat'//char(0)
    dflnm1c='./out/uz'//path(1:len_trim(path))//'.dat'//char(0)

    dflnm2a='./out/ux'//path1(1:len_trim(path1))//'.dat'//char(0)
    dflnm2b='./out/uy'//path1(1:len_trim(path1))//'.dat'//char(0)
    dflnm2c='./out/uz'//path1(1:len_trim(path1))//'.dat'//char(0)
 
    write(*,*) 'reading data files: ', path(1:len_trim(path))//'.dat'
    open(16,file=dflnm1a,form='unformatted')
      read(16) ux
    close(16)
    open(16,file=dflnm1b,form='unformatted')
      read(16) uy
    close(16)
    open(16,file=dflnm1c,form='unformatted')
      read(16) uz
    close(16)
 
    ux=g*ux
    uy=g*uy
    uz=g*uz
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
 
    uxn(1:nx:2,:,:) =  real(ux(1:lx,:,:))
    uxn(2:nx:2,:,:) = aimag(ux(1:lx,:,:))
    uyn(1:nx:2,:,:) =  real(uy(1:lx,:,:))
    uyn(2:nx:2,:,:) = aimag(uy(1:lx,:,:))
    uzn(1:nx:2,:,:) =  real(uz(1:lx,:,:))
    uzn(2:nx:2,:,:) = aimag(uz(1:lx,:,:))
 
    write(*,*) 'reading data files: ', path1(1:len_trim(path1))//'.dat'
    open(16,file=dflnm2a,form='unformatted')
      read(16) ux
    close(16)
    open(16,file=dflnm2b,form='unformatted')
      read(16) uy
    close(16)
    open(16,file=dflnm2c,form='unformatted')
      read(16) uz
    close(16)
    ux=g*ux
    uy=g*uy
    uz=g*uz
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
 
    uxnp(1:nx:2,:,:)= real(ux(1:lx,:,:))
    uxnp(2:nx:2,:,:)=aimag(ux(1:lx,:,:))
    uynp(1:nx:2,:,:)= real(uy(1:lx,:,:))
    uynp(2:nx:2,:,:)=aimag(uy(1:lx,:,:))
    uznp(1:nx:2,:,:)= real(uz(1:lx,:,:))
    uznp(2:nx:2,:,:)=aimag(uz(1:lx,:,:))
    
 
    do ndim = 1, ndir
 
      alpha = (ndim - 1) * pi / ndir
 
      rhat(1) = cos(alpha); rhat(2) = sin(alpha); rhat(3) = 0.
      shat(1) = cos(alpha + pi/2); shat(2) = sin(alpha + pi/2); shat(3) = 0.
      that(1) = 0.; that(2) = 0.; that(3) = 1.
 
      del_x = dl * rhat(1); del_y = dl * rhat(2); del_z = dl * rhat(3)
      
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
        xpe = (ii-1)*meshx + del_x
        ype = (jj-1)*meshy + del_y
        zpe = (kk-1)*meshz + del_z
        wxp = modulo(xpe, meshx) / meshx
        wyp = modulo(ype, meshy) / meshy
        wzp = modulo(zpe, meshz) / meshz
        nxp = 1 + modulo(floor(xpe/meshx),nx)
        nyp = 1 + modulo(floor(ype/meshy),ny)
        nzp = 1 + modulo(floor(zpe/meshz),nz)
  
        uxi = uxn(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uxn(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uxn(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uxn(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uxn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uxn(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uxn(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uxn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uyi = uyn(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uyn(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uyn(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uyn(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uyn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uyn(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uyn(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uyn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uzi = uzn(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uzn(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uzn(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uzn(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uzn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uzn(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uzn(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uzn(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
  
        xpe=uxi-uxn(ii,jj,kk)
        ype=uyi-uyn(ii,jj,kk)
        zpe=uzi-uzn(ii,jj,kk)
  
        ! Initial velocity increments
        del_ur=xpe*rhat(1)+ype*rhat(2)+zpe*rhat(3)
        del_us=xpe*shat(1)+ype*shat(2)+zpe*shat(3)
        del_ut=xpe*that(1)+ype*that(2)+zpe*that(3)
  
        ! Model rates of changes
        durdt(ii,jj,kk)=(del_ut*del_ut+del_us*del_us-del_ur*del_ur/3.)/dl + 2*omega*del_us
        dusdt(ii,jj,kk)=-2.*del_ur*del_us/dl-2.*omega*del_ur
        dutdt(ii,jj,kk)=-2.*del_ur*del_ut/dl
  
        ! dl at new time
        wxp=xpe*ndt*dt+del_x
        wyp=ype*ndt*dt+del_y
        wzp=zpe*ndt*dt+del_z
        xpe = sqrt(wxp*wxp + wyp*wyp + wzp*wzp)
  
        ! direction of dl at new time
        del_xp=wxp/xpe
        del_yp=wyp/xpe
        del_zp=wzp/xpe
  
        ! New locations 
        xp = (ii-1) * meshx + uxn(ii,jj,kk) * ndt * dt
        yp = (jj-1) * meshy + uyn(ii,jj,kk) * ndt * dt
        zp = (kk-1) * meshz + uzn(ii,jj,kk) * ndt * dt
  
        wxp = modulo(xp, meshx) / meshx
        wyp = modulo(yp, meshy) / meshy
        wzp = modulo(zp, meshz) / meshz
  
        nxp=1+modulo(floor(xp/meshx),nx)
        nyp=1+modulo(floor(yp/meshy),ny)
        nzp=1+modulo(floor(zp/meshz),nz)
  
        uxi = uxnp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uxnp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uxnp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uxnp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uxnp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uxnp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uxnp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uxnp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uyi = uynp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uynp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uynp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uynp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uynp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uynp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uynp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uynp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uzi = uznp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
             +uznp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
             +uznp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
             +uznp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
             +uznp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
             +uznp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
             +uznp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
             +uznp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
  
        ! New locations of the other ends of the line elements
        xpe=xp+del_xp*dl
        ype=yp+del_yp*dl
        zpe=zp+del_zp*dl
    
        ! Intepolation 
        wxp = modulo(xpe,meshx)/meshx
        wyp = modulo(ype,meshy)/meshy
        wzp = modulo(zpe,meshz)/meshz
        nxp = 1+modulo(floor(xpe/meshx),nx)
        nyp = 1+modulo(floor(ype/meshy),ny)
        nzp = 1+modulo(floor(zpe/meshz),nz)
    
        uxe = uxnp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
            + uxnp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
            + uxnp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
            + uxnp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
            + uxnp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
            + uxnp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
            + uxnp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
            + uxnp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uye = uynp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
            + uynp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
            + uynp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
            + uynp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
            + uynp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
            + uynp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
            + uynp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
            + uynp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
        uze = uznp(nxp,nyp,nzp)*(1._sp-wxp)*(1._sp-wyp)*(1._sp-wzp) &
            + uznp(modulo(nxp,nx)+1,nyp,nzp)*wxp*(1._sp-wyp)*(1._sp-wzp) &
            + uznp(nxp,modulo(nyp,ny)+1,nzp)*(1._sp-wxp)*wyp*(1._sp-wzp) &
            + uznp(nxp,nyp,modulo(nzp,nz)+1)*(1._sp-wxp)*(1._sp-wyp)*wzp &
            + uznp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,nzp)*wxp*wyp*(1._sp-wzp) &
            + uznp(modulo(nxp,nx)+1,nyp,modulo(nzp,nz)+1)*wxp*(1._sp-wyp)*wzp &
            + uznp(nxp,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*(1._sp-wxp)*wyp*wzp &
            + uznp(modulo(nxp,nx)+1,modulo(nyp,ny)+1,modulo(nzp,nz)+1)*wxp*wyp*wzp
  
        del_uxp=uxe-uxi
        del_uyp=uye-uyi
        del_uzp=uze-uzi
  
        shatn(1)=-del_yp; shatn(2)=del_xp; shatn(3)=0.
        shatn(:)=shatn(:)/sqrt(shatn(1)*shatn(1)+shatn(2)*shatn(2)+shatn(3)*shatn(3))
  
        thatn(1)=del_yp*shatn(3)-del_zp*shatn(2)
        thatn(2)=del_zp*shatn(1)-del_xp*shatn(3)
        thatn(3)=del_xp*shatn(2)-del_yp*shatn(1)
  
        del_urp=del_uxp*del_xp+del_uyp*del_yp+del_uzp*del_zp
        del_usp=del_uxp*shatn(1)+del_uyp*shatn(2)+del_uzp*shatn(3)
        del_utp=del_uxp*thatn(1)+del_uyp*thatn(2)+del_uzp*thatn(3)
  
        durdt_dns(ii,jj,kk)=(del_urp-del_ur)/(ndt*dt)
        dusdt_dns(ii,jj,kk)=(del_usp-del_us)/(ndt*dt)
        dutdt_dns(ii,jj,kk)=(del_utp-del_ut)/(ndt*dt)
      end do 
      end do
      end do
 
      ! Estimate the ranges for plotting
      if ( nsets .eq. 0 ) then
          mdusdt_dns0 = sum(dusdt_dns)*tt
          rdusdt_dns0 = sqrt(sum((dusdt_dns-mdusdt_dns0)**2)*tt)
          
          maxmx = mdusdt_dns0 + width*rdusdt_dns0
          minmx = mdusdt_dns0 - width*rdusdt_dns0
          bnwx = (maxmx - minmx)/real(npntx,sp)
      end if
 
      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
 
        aa = durdt_dns(ii,jj,kk)
        bb = durdt(ii,jj,kk)
 
        ll = floor( ( aa - minmx ) / bnwx ) + 1
        if ( ll .ge. 1 .and. ll .le. npntx) then
            cndur(ll) = cndur(ll) + bb
            pdfur(ll) = pdfur(ll) + 1
        end if
 
        aa = dusdt_dns(ii,jj,kk)
        bb = dusdt(ii,jj,kk)
 
        ll = floor( ( aa - minmx ) / bnwx ) + 1
        if ( ll .ge. 1 .and. ll .le. npntx) then
            cndus(ll) = cndus(ll) + bb
            pdfus(ll) = pdfus(ll) + 1
        end if
 
        aa = dutdt_dns(ii,jj,kk)
        bb = dutdt(ii,jj,kk)
 
        ll = floor( ( aa - minmx ) / bnwx ) + 1
        if ( ll .ge. 1 .and. ll .le. npntx) then
            cndut(ll) = cndut(ll) + bb
            pdfut(ll) = pdfut(ll) + 1
        end if

 
      end do
      end do
      end do
      nsets = nsets + 1
    end do
  end do

  tt = tt / nsets

  cndur = cndur / (pdfur + mytiny)
  cndus = cndus / (pdfus + mytiny)
  cndut = cndut / (pdfut + mytiny)

  open(15,file='rrrsrt'//cndel(1:len_trim(cndel))//'dx-cnd.dat')
    ! Format for gnuplot
    write(15,'(''# variables = "x", "cnddurdt", "cnddusdt, "cnddutdt"'')')
    do ii=1,npntx
      write(15,'(15e18.3)') minmx+(ii-.5)*bnwx, cndur(ii), cndus(ii), cndut(ii)
    end do
  close(15)

  deallocate(ux,uy,uz,kx,ky,kz,g)
  deallocate(durdt_dns,dusdt_dns,dutdt_dns,durdt,dusdt,dutdt)
  deallocate(uxn,uyn,uzn,uxnp,uynp,uznp)

  write(*,*) 'destroy plans'
  call destroyplan3d

  write(*,*) 'Finished'
  
end program cndredns
