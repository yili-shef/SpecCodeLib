program spectral_dns

  use mconstant
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, rnu, dtw, time, eps, fkmax, akp, u0
  logical  :: new, forced
  ! parameters should be input from 'parameter_dns.d' initially.

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  p
  real(sp),    allocatable, dimension (:,:,:) ::  k2,k2_e
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz

  integer  :: lx, ly, lz, lx1
  integer  :: i, ii, iii, ifile, istep
  real(sp) :: delta, S, dtt, ek, dt, dto, betadx, umax, twrite

  open(90,file='parameter_dns.d',status='unknown')
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
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
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
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  delta=pi/real(lx,sp)
  betadx=beta*delta

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(k2_e(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(p(lx1,ly,lz))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (new) then
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  endif
  
  ox=vx
  oy=vy
  oz=vz

  ! specify umax empirically as three times the rms value, in reference to
  ! Gaussian distribution
  umax=3.*sqrt(3.)*u0
  
  dto=timemax
  twrite=time
  ifile = idp
  istep=-1

  !============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

    istep = istep + 1
 
    ! Output the total energy in k space
    ! For monitoring purpose only. May want to reduce this part of calculation by
    ! increasing ieout.
    if (mod(istep,ieout).eq.0) then
      ! calculating skewness is relatively time consuming. May need it only when
      ! debugging.
      ! call skewness(vx,kx,s,nx,ny,nz)
      k2_e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)
      ek = sum(k2_e)
      write( *,*)istep,time,dto,ek,s
      write(25,*)istep,time,dto,ek,s

      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.*sqrt(2.*ek)
    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. smallest) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      
      ! Calculate pressure field.
      call pressure_dns(vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,nx,ny,nz,p)
      call output_p(p,ii,lx1,ly,lz)

      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      k2_e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)
      write(24,*)time
      do i=1,lx
        ek=sum(k2_e,mask=(abs(sqrt(k2)-i).lt.0.5_sp))
        write(24,*)i,ek
      end do
 
    end if

    ! =========== Calculate vorticity =============================
    do iii=1,lz
    do i=1,lx1
      wx(i,:,iii)=ky(:)*vz(i,:,iii)
      wz(i,:,iii)=ky(:)*vx(i,:,iii)
    end do
    end do
    do ii=1,ly
    do i=1,lx1
      wx(i,ii,:)=wx(i,ii,:)-kz(:)*vy(i,ii,:)
      wy(i,ii,:)=kz(:)*vx(i,ii,:)
    end do
    end do
    do iii=1,lz
    do ii=1,ly
      wy(:,ii,iii)=wy(:,ii,iii)-kx(:)*vz(:,ii,iii)
      wz(:,ii,iii)=kx(:)*vy(:,ii,iii)-wz(:,ii,iii)
    end do
    end do
    wx=eye*wx; wy=eye*wy; wz=eye*wz
 
    !========== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========== combine convection term, and forcing ============
    if (forced) call force(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fkmax,eps)
 
    call symmetrize(wx,k2,lx1,ly,lz)
    call symmetrize(wy,k2,lx1,ly,lz)
    call symmetrize(wz,k2,lx1,ly,lz)
 
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
      vx = vx * (1._sp-rnu*dt*k2) + dt*wx 
      vy = vy * (1._sp-rnu*dt*k2) + dt*wy 
      vz = vz * (1._sp-rnu*dt*k2) + dt*wz 
 
      time = time + dt
      dto=dt
      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto

      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2

      vx = S*ox + ek*vx + delta*(wx - rnu*k2*vx) 
      vy = S*oy + ek*vy + delta*(wy - rnu*k2*vy) 
      vz = S*oz + ek*vz + delta*(wz - rnu*k2*vz) 
 
      time = time + dt
      dto=dt+dto
      call input(ox,oy,oz,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      k2_e=exp(-rnu*k2*dto)

      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt

      vx = vx + S*wx + ek*ox*k2_e
      vy = vy + S*wy + ek*oy*k2_e
      vz = vz + S*wz + ek*oz*k2_e
 
      k2_e=exp(-rnu*k2*dt)
      ! k2_e=k2_e**dtt ! Is this faster than evaluating exp?

      vx = vx*k2_e
      vy = vy*k2_e
      vz = vz*k2_e
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      dto=dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(kx, ky, kz, k2, k2_e, p)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral_dns
