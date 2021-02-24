program tracking
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer, parameter :: norder = 6

  integer :: nx, ny, nz, lx, ly, lz, lx1, nprtcle, forward, ndel
  integer :: nstep, nfstart
  real(sp) :: dt

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension(:,:,:) :: uxr,uyr,uzr
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(sp),    allocatable, dimension(:)    :: xp,yp,zp
  real(sp),    allocatable, dimension(:)    :: xps,yps,zps

  real(sp) :: y(3), dxyz(3), bg(norder,3)
  integer  :: ix(norder), iy(norder), iz(norder), lhnode(3) 
  integer  :: ii, jj, kk, ll, ip, ix1, iy1, iz1
  real(sp) :: uxp, uyp, uzp, ignore_me, const, delta_c
  
  character(80) :: str, str1, prefix, strnfstart

  write(*,*) 
  write(*,'(''>>> Particle tracking starting from given initial coordinates and velocity <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 8) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./particletracking.x nx nfstart nstep ndel nprtcle dt forward prefix'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        nfstart: the starting file number'
          write(*,*) '        nstep: # of steps of tracking (generating nstep+1 data files)'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        nprtcle: number of particles'
          write(*,*) '        dt: time step size'
          write(*,*) '        forward: 1 is forward tracking, 0 backward'
          write(*,*) '        prefix: prefix for the coordinates data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(3,strnfstart)
  read(strnfstart, '(I20)') nstep

  ! starting file number string
  call getarg(2,strnfstart)
  read(strnfstart, '(I20)') nfstart
  strnfstart = adjustl(strnfstart)

  ! filter scale
  call getarg(4,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! number of particle
  call getarg(5,str1)
  read(str1,'(I20)') nprtcle

  ! time step size
  call getarg(6,str1)
  read(str1,'(F20.10)') dt

  ! flag for forward or backward 
  call getarg(7, str1) 
  read(str1, '(I20)') forward

  if (forward .eq. 0) dt = -dt

  ! prefix for coordinates data file
  call getarg(8, prefix)
  prefix = adjustl(prefix)

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( uxr(nx,ny,nz), uyr(nx,ny,nz), uzr(nx,ny,nz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate( xps(nprtcle), yps(nprtcle), zps(nprtcle) )
  allocate( g(lx1,ly,lz) )
  write(*,*) 'Arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'After wavenumber'

  g = exp(-g*delta_c**2 / 24)

  open(16,file='./out/ux'//strnfstart(1:len_trim(strnfstart))//'.dat',form='unformatted')
    read(16) ux
  close(16)
  open(16,file='./out/uy'//strnfstart(1:len_trim(strnfstart))//'.dat',form='unformatted')
    read(16) uy
  close(16)
  open(16,file='./out/uz'//strnfstart(1:len_trim(strnfstart))//'.dat',form='unformatted')
    read(16) uz
  close(16)

  write(*,*) 'partiletracking.x: data file: ', strnfstart(1:len_trim(strnfstart))

  ux = ux * g
  uy = uy * g
  uz = uz * g

  call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
 
  uxr(1:nx:2,:,:)=real(ux(1:lx,:,:),sp); uxr(2:nx:2,:,:)=aimag(ux(1:lx,:,:))
  uyr(1:nx:2,:,:)=real(uy(1:lx,:,:),sp); uyr(2:nx:2,:,:)=aimag(uy(1:lx,:,:))
  uzr(1:nx:2,:,:)=real(uz(1:lx,:,:),sp); uzr(2:nx:2,:,:)=aimag(uz(1:lx,:,:))

  str1 = './out/'//prefix(1:len_trim(prefix))//strnfstart(1:len_trim(strnfstart))//'.dat'
  open(16, file = str1(1:len_trim(str1)), form = 'unformatted')
    read(16) xp, yp, zp
  close(16)


  nstep = nstep -1
  do while ( nstep .ge. 0 )

    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)

      uxp = 0._sp; uyp = 0._sp; uzp = 0._sp
      do  ii=1,norder
        ix1 = ix(ii)
        do  jj=1,norder
          iy1 = iy(jj)
          do  kk=1,norder
            iz1 = iz(kk)

            uxp = uxp + uxr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)
            uyp = uyp + uyr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)
            uzp = uzp + uzr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)

          end do
        end do
      end do
 
      xps(ip)=xp(ip)+dt*uxp
      yps(ip)=yp(ip)+dt*uyp
      zps(ip)=zp(ip)+dt*uzp
      
    end do
 
    if (forward .eq. 1) then 
        nfstart = nfstart + 1
    else
        nfstart = nfstart - 1
    end if

    write(str1,'(I8)') nfstart
    str1 = adjustl(str1)
 
    open(16,file='./out/ux'//str1(1:len_trim(str1))//'.dat',form='unformatted')
      read(16) ux
    close(16)
    open(16,file='./out/uy'//str1(1:len_trim(str1))//'.dat',form='unformatted')
      read(16) uy
    close(16)
    open(16,file='./out/uz'//str1(1:len_trim(str1))//'.dat',form='unformatted')
      read(16) uz
    close(16)
 
    write(*,*) 'particle tracking.x: data file: ', str1(1:len_trim(str1))
    nstep = nstep - 1

    ux = ux * g
    uy = uy * g
    uz = uz * g

    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
    
    uxr(1:nx:2,:,:)=real(ux(1:lx,:,:),sp); uxr(2:nx:2,:,:)=aimag(ux(1:lx,:,:))
    uyr(1:nx:2,:,:)=real(uy(1:lx,:,:),sp); uyr(2:nx:2,:,:)=aimag(uy(1:lx,:,:))
    uzr(1:nx:2,:,:)=real(uz(1:lx,:,:),sp); uzr(2:nx:2,:,:)=aimag(uz(1:lx,:,:))
 
    do ip = 1, nprtcle

      y(1)=xps(ip)
      y(2)=yps(ip)
      y(3)=zps(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)

      uxp = 0._sp; uyp = 0._sp; uzp = 0._sp
      do  ii=1,norder
        ix1 = ix(ii)
        do  jj=1,norder
          iy1 = iy(jj)
          do  kk=1,norder
            iz1 = iz(kk)

            uxp = uxp + uxr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)
            uyp = uyp + uyr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)
            uzp = uzp + uzr(ix1,iy1,iz1) * bg(ii,1) * bg(jj,2) * bg(kk,3)

          end do
        end do
      end do
 
      xp(ip)=.5*(xp(ip)+xps(ip))+.5*dt*uxp
      yp(ip)=.5*(yp(ip)+yps(ip))+.5*dt*uyp
      zp(ip)=.5*(zp(ip)+zps(ip))+.5*dt*uzp
 
    end do

    str1 = './out/'//prefix(1:len_trim(prefix))//str1(1:len_trim(str1))//'.dat'
    open(16, file = str1(1:len_trim(str1)), form = 'unformatted')
      write(16) xp, yp, zp
    close(16)

  end do

  call destroyplan3d

  deallocate(kx, ky, kz, g, ux, uy, uz, uxr, uyr, uzr, xp, yp, zp) 
  deallocate(xps, yps, zps)

  write(*,*) 'particletracking.x: finished'
  stop

end program tracking      
