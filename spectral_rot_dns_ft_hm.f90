program spectral

  use mconstant
  use mwavenumber
  use msymmetrize
  use mheliwave
  use minitialize
  use minput
  use moutput
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, nstep, itout, ieout, idp
  real(sp) :: dt, time, eps, fk1, fk2, akp, u0, om, rnu
  logical  :: new, forced
  ! parameters should be input from 'parameter_rot_dns.d' initially.

  integer  :: lx, ly, lz, lx1
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am,efp,efm
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  real(sp),    allocatable, dimension (:,:,:) ::  k2,kx,ky,kz
  real(sp),    allocatable, dimension (:)     ::  eks

  integer  :: ii, jj, kk, ll, istep, jstep, ifile
  real(sp) :: s, dt_h, ek, kcut, kcut2, fk12, fk22

  open(90,file='parameter_rot_dns_ft.d',status='unknown')
    read(90,*) forced
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) nstep
    read(90,*) itout
    read(90,*) ieout
    read(90,*) dt
    read(90,*) rnu
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) eps
    read(90,*) fk1
    read(90,*) fk2
    read(90,*) akp
    read(90,*) u0
    read(90,*) om
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'nstep', nstep
  write(*,*) 'itout', itout
  write(*,*) 'ieout', ieout
  write(*,*) 'dt', dt
  write(*,*) 'rnu', rnu
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'om', om
  write(*,*) 'iseed', iseed

  !if (forced) then
  !  om=eps**(1./3.)*(.5*(fk1+fk2))**(2./3.)/ro   ! ro is defined by forcing scales
  !else
  !  om=u0*akp/(2._sp*ro*pi)  ! ro is defined by initial integral scales
  !end if

  !if (ro .gt. 1000) om = 0.
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz
  
  lx1=lx+1
  kcut = 2.*lx/3
  kcut2 = kcut * kcut
  fk12 = fk1 * fk1
  fk22 = fk2 * fk2


  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(efp(lx1,ly,lz),efm(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  efp = exp(-k2*rnu*dt+2._sp*eye*dt*kz*om/(sqrt(k2)+mytiny))  ! positive helical wave
  efm = exp(-k2*rnu*dt-2._sp*eye*dt*kz*om/(sqrt(k2)+mytiny))  ! negative helical wave

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'


  if (new) then
    ! generate initial condition -------
    ! The velocity fields have been dealiased in initialize
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'after initialize'
  else
    ! reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  endif
  
  ox=vx
  oy=vy
  oz=vz
  ap=vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am=vx*hpx+vy*hpy+vz*hpz

  dt_h=.5_sp*dt
  ifile = idp
  istep=0
  jstep=0 

!============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (jstep .le. nstep)

    ! write out the total energy in k space

    if (mod(jstep,ieout).eq.0) then

      !call skewness(vx,kx,s,nx,ny,nz)

      wx = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      wx(1,:,:)=0.5_sp*wx(1,:,:)
      ek = sum(real(wx))
      write( *,'(i6,6f12.6)') jstep,time,ek,s
      write(25,'(i6,6f12.6)') jstep,time,ek,s
    end if
 
    ! output the field data and spectra, including the initial data when jstep=0.
    if (mod(jstep,itout).eq.0) then

      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      
      ! wx = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where wx=.5*q_ij(k)
      wx = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      wx(1,:,:)=0.5_sp*wx(1,:,:)
      write(24,*)time
      eks = 0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + real(wx(ii,jj,kk))
            end if
      end do
      end do
      end do
      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
 
    end if

    ! =========== calculate vorticity =============================
    wx = eye*(ky*vz - kz*vy)
    wy = eye*(kz*vx - kx*vz)
    wz = eye*(kx*vy - ky*vx)
  
    !=================== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========================== forcing term =================================  
    if (forced) call force_rot(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fk12,fk22,eps)
 
    !call symmetrize(wx,lx1,ly,lz)
    !call symmetrize(wy,lx1,ly,lz)
    !call symmetrize(wz,lx1,ly,lz)
 
    !======================== projecting ================================= 
    !call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
    vx = (kx*wx + ky*wy + kz*wz)/k2
    wx = wx - kx * vx 
    wy = wy - ky * vx 
    wz = wz - kz * vx 
 

    !======================= time stepping ======================================
    select case (istep)
    case (0)
      ap = ( ap + dt_h*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz)) ) * sqrt(efp)
      am = ( am + dt_h*(wx*hpx+wy*hpy+wz*hpz)                      ) * sqrt(efm)
      
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! only half of first step, set to one to avoid duplicate output

      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      ! ox,oy,oz are the initial velocity components now
      ap = (ox*conjg(hpx)+oy*conjg(hpy)+oz*conjg(hpz))*efp + & 
            dt*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz))*sqrt(efp)
      am = (ox*hpx+oy*hpy+oz*hpz)*efm + dt*(wx*hpx+wy*hpy+wz*hpz)*sqrt(efm)
 
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! the actual first step.

      call input(ox,oy,oz,0,lx1,ly,lz)
      ! after the input, ox,oy,oz store the components of the rhs of the equation
      ! at previous time step: wx,wy,wz.
 
    case (2:)
      vx = 3._sp*wx - efp*ox
      vy = 3._sp*wy - efp*oy
      vz = 3._sp*wz - efp*oz
      ap = (ap + dt_h*(vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)))*efp
      
      vx = 3._sp*wx - efm*ox
      vy = 3._sp*wy - efm*oy
      vz = 3._sp*wz - efm*oz
      am = (am + dt_h*(vx*hpx+vy*hpy+vz*hpz))*efm
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      istep = istep + 1
      jstep = jstep + 1
 
    case default
      write(*,*) "istep wrong. istep = ", istep
      stop
    end select
    vx=ap*hpx+am*conjg(hpx)
    vy=ap*hpy+am*conjg(hpy)
    vz=ap*hpz+am*conjg(hpz)

    ! dealiasing
    where(k2 .ge. kcut2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
    vx(1,1,1) = 0._sp
    vy(1,1,1) = 0._sp
    vz(1,1,1) = 0._sp

  end do
!=================== end of loop ==============================

  close(24)
  close(25)

!  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(ap, am, hpx, hpy, hpz)
  deallocate(kx, ky, kz, k2, efp, efm)

!  destroy plans
  call destroyplan3d

  write(*,*) 'spectral_rot_dns_ft.x finished'
  stop

end program spectral
