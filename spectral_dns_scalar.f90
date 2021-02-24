program spectral_dns_scalar

  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  real(sp), parameter :: cphiuz = -0.3 !TODO: find the correct value.

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, rnu, dtw, time, eps, epsphi, fkmax, akp, u0, pr, phirms
  logical  :: new, forced, scalarforced
  ! parameters should be input from 'parameter_scalar_dns.d' initially.

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,phi
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz,wphi
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz,ophi
  real(sp),    allocatable, dimension (:,:,:) ::  k2,k2_e
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz, eks

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: delta, S, dtt, ek, dt, dto, betadx, umax, twrite
  real(sp) :: meanphigrad

  open(90,file='parameter_scalar_dns.d',status='unknown')
    read(90,*) forced
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) timemax
    read(90,*) ieout
    read(90,*) beta
    read(90,*) rnu
    read(90,*) dtw
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) eps
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    ! parameter for scalar advection
    read(90,*) scalarforced
    read(90,*) pr
    read(90,*) epsphi
    read(90,*) phirms

    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'scalarforced', scalarforced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'beta', beta
  write(*,*) 'rnu', rnu
  write(*,*) 'dtw', dtw
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'prandtl number', pr
  write(*,*) 'eps of phi', epsphi
  write(*,*) 'rms of phi', phirms
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  delta = pi / real(lx,sp)
  betadx = beta * delta

  ! Mean scalar gradient to generate scalar dissipation epsphi
  meanphigrad = - epsphi/(u0*phirms*cphiuz)

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(k2_e(lx1,ly,lz),phi(lx1,ly,lz),wphi(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(ophi(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (new) then
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    call initializephi(phi,k2,iseed,lx1,ly,lz,akp,phirms)

    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    call inputphi(phi,idp,lx1,ly,lz)
    write(*,*) 'initial velocity and scalar field readed.'
  endif

  
  ox=vx
  oy=vy
  oz=vz
  ophi=phi

  ! specify umax empirically as three times the rms value, in reference to
  ! Gaussian distribution
  umax=3._sp*sqrt(3._sp)*u0
  
  dto = timemax
  twrite = time
  ifile = idp
  istep = -1

  !============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  open(26,file='./post/scalar_spectrum.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

    istep = istep + 1
 
    ! Output the total energy in k space
    ! For monitoring purpose only. May want to reduce this part of calculation by
    ! increasing ieout.
    if (mod(istep,ieout).eq.0) then
        
      ! calculating skewness is relatively time consuming. May need it only when
      ! debugging.
      !call skewness(vx,kx,s,nx,ny,nz)

      k2_e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)
      ek = sum(k2_e)

      k2_e = phi * conjg(phi)
      k2_e(1,:,:) = .5_sp * k2_e(1,:,:)
      phirms = sum(k2_e)

      write( *,'(I6, 10E12.5)')istep,time,dto,ek,s,phirms
      write(25,'(I6, 10E12.5)')istep,time,dto,ek,s,phirms

      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3._sp*sqrt(2._sp*ek)
    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. smallest) then

      call output(vx,vy,vz,ifile,lx1,ly,lz)
      call outputphi(phi,ifile,lx1,ly,lz)

      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      k2_e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)

      write(24,*) '#', time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ll=floor(sqrt(k2(ii,jj,kk))+.5)
        if (ll .ge. 1 .and. ll .le. lx) then
          eks(ll) = eks(ll) + k2_e(ii,jj,kk)
        end if
      end do
      end do
      end do

      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
      write(24,*)

      ! Scalar spectrum
      k2_e = phi * conjg(phi)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)

      write(26,*) '#', time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + k2_e(ii,jj,kk)
            end if
      end do
      end do
      end do

      do ii = 1, lx
        write(26,*) ii, eks(ii)
      end do
      write(26,*)

 
    end if

    ! =========== Calculate vorticity =============================
    
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk) = eye * ky(jj) * vz(ii,jj,kk) - eye * kz(kk) * vy(ii,jj,kk)
      wy(ii,jj,kk) = eye * kz(kk) * vx(ii,jj,kk) - eye * kx(ii) * vz(ii,jj,kk)
      wz(ii,jj,kk) = eye * kx(ii) * vy(ii,jj,kk) - eye * ky(jj) * vx(ii,jj,kk)
    end do
    end do
    end do
 
    !========== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
    call convecphi_dns(vx,vy,vz,phi,wphi,kx,ky,kz,lx1,ly,lz,nx,ny,nz)
 
    !========== combine convection term, and forcing ============
    if (forced) call force(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fkmax,eps)

    !if (scalarforced) call forcephi(phi,wphi,k2,lx1,ly,lz,fkmax,epsphi)

    ! Forcing scalar with a mean scalar gradient in z direction. 
    if (scalarforced) then 
        wphi = wphi - meanphigrad * vz
    end if
 
    call symmetrize(wx,k2,lx1,ly,lz)
    call symmetrize(wy,k2,lx1,ly,lz)
    call symmetrize(wz,k2,lx1,ly,lz)
    call symmetrize(wphi,k2,lx1,ly,lz)

    !========== projection ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
 
    !========== time stepping ===================================
    ! umax is updated only every ieout steps
    dt=min(betadx/umax,3._sp*dto)
    if (dt .ge. twrite-time) then
            dt=twrite-time
    else if (twrite-time-dt .lt. .2_sp*dt) then
            dt=.5_sp*(twrite-time)
    end if

    select case (istep)
    case (0)
      k2_e=exp(-rnu*k2*dt)
      vx = ( vx + dt*wx ) * k2_e 
      vy = ( vy + dt*wy ) * k2_e 
      vz = ( vz + dt*wz ) * k2_e 

      k2_e = exp(-rnu*k2*dt/pr)
      phi = ( phi  + dt * wphi ) * k2_e
 
      time = time + dt
      dto=dt
      call output(wx,wy,wz,0,lx1,ly,lz)
      call outputphi(phi,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto

      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2

      k2_e=exp(-rnu*k2*dto)

      ! ox,oy,oz are the initial velocity components now
      vx = S*ox*k2_e + ek*vx + delta*wx 
      vy = S*oy*k2_e + ek*vy + delta*wy 
      vz = S*oz*k2_e + ek*vz + delta*wz 

      k2_e = exp(-rnu*k2*dto/pr)
      phi = S * ophi * k2_e + ek * phi + delta * wphi
 
      k2_e=exp(-rnu*k2*dt)
      vx = vx*k2_e
      vy = vy*k2_e
      vz = vz*k2_e

      k2_e = exp(-rnu*k2*dt/pr)
      phi = phi * k2_e

      time = time + dt
      dto=dt+dto
      call input(ox,oy,oz,0,lx1,ly,lz)
      call inputphi(ophi,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt

      k2_e=exp(-rnu*k2*dto)
      vx = vx + S*wx + ek*ox*k2_e
      vy = vy + S*wy + ek*oy*k2_e
      vz = vz + S*wz + ek*oz*k2_e

      k2_e = exp(-rnu*k2*dto/pr)
      phi = phi + S * wphi + ek * ophi * k2_e
 
      k2_e=exp(-rnu*k2*dt)
      vx = vx*k2_e
      vy = vy*k2_e
      vz = vz*k2_e

      k2_e = exp(-rnu*k2*dt/pr)
      phi = phi * k2_e
 
      ox = wx
      oy = wy
      oz = wz

      ophi = wphi
 
      time = time + dt
      dto=dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)
  close(26)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(phi, wphi, ophi)
  deallocate(kx, ky, kz, k2, k2_e, eks)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral_dns_scalar
