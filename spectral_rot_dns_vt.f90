program spectral

  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, dtw, rnu, time, eps, fk1, fk2, akp, u0, hrel, ro
  logical  :: new, forced
  ! parameters should be input from 'parameter_rot_dns.d' initially.

  integer  :: lx, ly, lz, lx1
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am,ef
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  real(sp),    allocatable, dimension (:,:,:) ::  k2
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz

  integer  :: i, ii, iii, istep, ifile, q
  real(sp) :: delta, om, s, dtt, dt, ek, dto, betadx, umax, twrite

  open(90,file='parameter_rot_dns.d',status='unknown')
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
    read(90,*) fk1
    read(90,*) fk2
    read(90,*) akp
    read(90,*) u0
    read(90,*) iseed
    read(90,*) q
    read(90,*) ro
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
  write(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  write(*,*) 'q', q
  write(*,*) 'ro', ro

  om=u0*akp/(2._sp*ro*pi)
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  delta=pi/real(lx,sp)
  betadx=beta*delta

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(ef(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'


  if (new) then
! generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    write(*,*) 'after initialize'
  else
! rreading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    print*, 'initial velocity field readed.'
  endif
  
  ox=vx
  oy=vy
  oz=vz
  ap=vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am=vx*hpx+vy*hpy+vz*hpz

  ! specify umax empirically as three times the rms value, in reference to
  ! gaussian distribution
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
  
    !write out the total energy in k space
 
    if (mod(istep,ieout).eq.0) then
      call skewness(vx,kx,s,lx1,ly,lz,nx,ny,nz)
      ef = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      ef(1,:,:)=0.5_sp*ef(1,:,:)
      ek = sum(real(ef))
      write(*,'(i6,6e12.4)')istep,time,ek,s
      write(25,'(i6,6e12.4)')istep,time,ek,s
 
      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.*sqrt(2.*ek)
    end if
 
    ! output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. smallest) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! ef = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where ef=.5*q_ij(k)
      ef= vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      ef(1,:,:)=0.5_sp*ef(1,:,:)
      write(24,*)time
      do i=1,lx
        ek=sum(real(ef),mask=(abs(sqrt(k2)-i).lt.0.5_sp))
        write(24,*)i,ek
      end do
    end if
  
      ! =========== calculate vorticity =============================
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
  
 
    !=================== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========================== forcing term =================================  
    if (forced) call force_rot(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fk1,fk2,eps)
 
    call symmetrize(wx,k2,lx1,ly,lz)
    call symmetrize(wy,k2,lx1,ly,lz)
    call symmetrize(wz,k2,lx1,ly,lz)
 
    !======================== projecting ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
 
    !========== time stepping ===================================
    ! umax is updated only every ieout steps
    dt=min(betadx/umax,3._sp*dto)
    if (dt .ge. twrite-time) then
            dt=twrite-time
    else if (twrite-time-dt .lt. .2_sp*dt) then
            dt=.5_sp*(twrite-time)
    end if
 
    !======================= time stepping ======================================
    select case (istep)
    case (0)
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)**q*dt*rnu+2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! positive helical wave
      enddo
      enddo 
      ap=ap*ef+dt*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz))
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dt*rnu-2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! negative helical wave
      enddo
      enddo 
      am=am*ef+dt*(wx*hpx+wy*hpy+wz*hpz)
      
      time = time + dt
      dto=dt
      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto
 
      s= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dto*rnu+2._sp*eye*dto*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! positive helical wave
      end do
      end do 
 
      ! ox,oy,oz are the initial velocity components now
      ap = s*(ox*conjg(hpx)+oy*conjg(hpy)+oz*conjg(hpz))*ef + ek*ap + &
           delta*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz))
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dt*rnu+2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! positive helical wave
      end do
      end do 
      ap = ap*ef
 
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dto*rnu-2._sp*eye*dto*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! negative helical wave
      end do
      end do 
 
      ! ox,oy,oz are the initial velocity components now
      am = s*(ox*hpx+oy*hpy+oz*hpz)*ef + ek*am + delta*(wx*hpx+wy*hpy+wz*hpz)
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dt*rnu-2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! negative helical wave
      end do
      end do 
      am = am*ef
 
      time = time + dt
      dto=dt+dto
      call input(ox,oy,oz,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      s=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dto*rnu+2._sp*eye*dto*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! positive helical wave
      end do
      end do 
 
      vx = s*wx + ek*ox*ef
      vy = s*wy + ek*oy*ef
      vz = s*wz + ek*oz*ef
      ap = ap + vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dt*rnu+2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! positive helical wave
      end do
      end do 
      ap=ap*ef
 
      
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dto*rnu-2._sp*eye*dto*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! negative helical wave
      end do
      end do 
 
      vx = s*wx + ek*ox*ef
      vy = s*wy + ek*oy*ef
      vz = s*wz + ek*oz*ef
      am = am + vx*hpx+vy*hpy+vz*hpz
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = exp(-k2(i,ii,:)*dt*rnu-2._sp*eye*dt*kz*om/(sqrt(k2(i,ii,:))+smallest))  ! negative helical wave
      end do
      end do 
      am=am*ef
 
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      dto=dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select
    vx=ap*hpx+am*conjg(hpx)
    vy=ap*hpy+am*conjg(hpy)
    vz=ap*hpz+am*conjg(hpz)

  end do
!=================== end of loop ==============================

  close(24)
  close(25)

!  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(ap, am, hpx, hpy, hpz)
  deallocate(kx, ky, kz, k2, ef)

!  destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral
