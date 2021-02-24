program spectral

  use mconstant
  use mwavenumber
  use minitialize
  use minput
  use moutput
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, nstep, itout, ieout, idp
  real(sp) :: dt, time, eps, fkmax, akp, u0, rnu
  logical  :: new, forced
  ! parameters should be input from 'parameter_dns_ft.d' initially.

  integer  :: lx, ly, lz, lx1
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  real(sp),    allocatable, dimension (:,:,:) ::  k2, expk2
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz, eks

  integer  :: ii, jj, kk, ll, istep, jstep, ifile
  real(sp) :: s, dt_h, ek, kcut, kcut2, fkmax2, kinet, coeforce

  open(90,file='parameter_dns_ft.d',status='unknown')
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
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
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
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz
  
  lx1=lx+1
  kcut = 2.*lx/3
  kcut2 = kcut * kcut
  fkmax2 = fkmax * fkmax


  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(expk2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  expk2 = exp(-k2*rnu*dt)
  expk2(1,1,1) = 1._sp

  if (new) then
    ! ----- generate initial condition -------
    ! The velocity fields have been dealiased in initialize
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'after initialize'
  else
    ! ----- reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  endif
  
  ox=vx
  oy=vy
  oz=vz

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
      write(*,'(i6,6f12.4)') jstep,time,ek,s
      write(25,'(i6,6f12.4)')jstep,time,ek,s
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

    ! ================== coefficients for forcing terms =====================
    where (k2 .lt. fkmax2)
      wx = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
    elsewhere
      wx = 0.
    end where
    wx(1,:,:) = .5_sp * wx(1,:,:)
    kinet = sum( real(wx) )
 
    coeforce = eps / (2._sp * kinet)

    ! =========== calculate vorticity =============================
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      wx(ii,jj,kk) = eye*(ky(jj)*vz(ii,jj,kk) - kz(kk)*vy(ii,jj,kk))
      wy(ii,jj,kk) = eye*(kz(kk)*vx(ii,jj,kk) - kx(ii)*vz(ii,jj,kk))
      wz(ii,jj,kk) = eye*(kx(ii)*vy(ii,jj,kk) - ky(jj)*vx(ii,jj,kk))
    end do
    end do
    end do
  
    !=================== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========================== forcing term =================================  
    ! if (forced) call force(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fkmax,eps)
    if (forced) then
        where ( k2 .lt. fkmax2 )
            wx = wx + coeforce * vx
            wy = wy + coeforce * vy
            wz = wz + coeforce * vz
        endwhere
    end if
 
    !======================== projecting ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)

    !======================= time stepping ======================================
    select case (istep)
    case (0)
      vx = ( vx + dt_h * wx )  * sqrt(expk2)
      vy = ( vy + dt_h * wy )  * sqrt(expk2)
      vz = ( vz + dt_h * wz )  * sqrt(expk2)
      
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! only half of first step, set to one to avoid duplicate output

      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      ! ox,oy,oz are the initial velocity components now
      vx = ox * expk2 + dt * wx * sqrt(expk2)
      vy = oy * expk2 + dt * wy * sqrt(expk2)
      vz = oz * expk2 + dt * wz * sqrt(expk2)
 
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! the actual first step.

      call input(ox,oy,oz,0,lx1,ly,lz)
      ! after the input, ox,oy,oz store the components of the rhs of the equation
      ! at previous time step: wx,wy,wz.
 
    case (2:)
      vx = ( vx + ( 3._sp*wx - expk2*ox ) * dt_h ) * expk2
      vy = ( vy + ( 3._sp*wy - expk2*oy ) * dt_h ) * expk2
      vz = ( vz + ( 3._sp*wz - expk2*oz ) * dt_h ) * expk2
      
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      istep = istep + 1
      jstep = jstep + 1
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

    ! dealiasing
    where(k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere
    vx(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    vy(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    vz(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
  

  end do
!=================== end of loop ==============================

  close(24)
  close(25)

!  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(kx, ky, kz, k2, expk2)

!  destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral
