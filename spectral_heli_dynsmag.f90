program spectral

  use mconstant
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, forcetype, iseed, nstep, itout, ieout, idp, iupdate
  real(sp) :: dt, rnu, time, eps, eta, fkmax, akp, u0, hrel, delta_ratio
  logical  :: new, update
  ! parameters should be input from 'parameter_heli_dynsmag.d' initially.

  real(sp) :: delta, delta_test,k_test,cs2
  integer  :: nxb, nyb, nzb, lx, ly, lz, lxb, lyb, lzb, lx1, lxb1
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz,hmx,hmy,hmz
  complex(sp), allocatable, dimension (:,:,:) ::  fx,fy,fz
  complex(sp), allocatable, dimension (:,:,:) ::  t11,t12,t13,t22,t23,t33
  real(sp),    allocatable, dimension (:,:,:) ::  kx,ky,kz,k2,g_test
  real(sp),    allocatable, dimension (:,:,:) ::  k2_e,tmp

  integer  :: i, ii, iii, istep, jstep
  real(sp) :: s, dt_h, ek, hek, lm, mm

  open(90,file='parameter_heli_dynsmag.d',status='unknown')
    read(90,*) forcetype
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
    read(90,*) eta
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    read(90,*) hrel
    read(90,*) iupdate
    read(90,*) delta_ratio
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forcetype', forcetype
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
  write(*,*) 'eta', eta
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'hrel', hrel
  write(*,*) 'iupdate', iupdate
  write(*,*) 'delta_ratio', delta_ratio
  write(*,*) 'iseed', iseed
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  nxb=nx*3/2
  nyb=ny*3/2
  nzb=nz*3/2

  lx=nx/2
  ly=ny
  lz=nz

  lxb=lx*3/2
  lyb=ly*3/2
  lzb=lz*3/2
  
  lx1=lx+1
  lxb1=lxb+1

  delta=pi/real(lx,sp)
  delta_test=delta_ratio*delta
  k_test=lx/delta_ratio

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(k2_e(lx1,ly,lz),tmp(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(fx(lx1,ly,lz),fy(lx1,ly,lz),fz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(hmx(lx1,ly,lz),hmy(lx1,ly,lz),hmz(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g_test(lx1,ly,lz))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  k2_e = exp(-k2*dt*rnu)
  k2_e(1,1,1) = 1.0_sp
  write(*,*) 'after wavenumber'

  where(k2 .le. k_test*k_test)
          g_test=1._sp
  elsewhere
          g_test=0._sp
  endwhere
  write(*,*) 'after g_test'

  call heliwave(hpx,hpy,hpz,hmx,hmy,hmz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'


  if (new) then
! generate initial condition -------
    call initialize_heli(hpx,hpy,hpz,hmx,hmy,hmz,vx,vy,vz,k2,lx1,ly,lz,iseed,hrel,akp,u0)
    write(*,*) 'after initialize_heli'
  else
! rreading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'after input'
  endif
  
  ox=vx
  oy=vy
  oz=vz

  dt_h=.5_sp*dt
  ii = idp
  istep=-1
  jstep=-1 
  ! istep counts the actual steps of time advancing.
  ! because the step size for the first  two steps are half of dt,
  ! for the output purpose, the first two steps are counted as one
  ! step and the corresponding number of steps, with step size being dt, 
  ! is counted by jstep. therefore jstep counts the number of steps on each
  ! of which dt has been advanced.
  ! before the first step, istep=jstep=0; when
  ! istep=1 before the next actual time step, jstep is arbitrarily set (since 
  ! this step is not counted in jstep) to be itout+1 to avoid output 
  ! on this step (requiring itout>1). then before further next step, istep=2,
  ! jstep=1.

!============== begin the loop ====================================

  open(23,file='./post/helicity_spectrum.data',status='unknown') 
  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  open(26,file='./post/lm.data',status='unknown')
  write(25,*) '# jstep,time,ener,heli,s,cs2'
  write(26,*) '# jstep,time,lm,mm'
  
  write(*,*) 'starting the loop'
  do while (jstep .le. nstep)

  istep = istep + 1
  jstep = jstep + 1

! calculate vorticity
  wx = eye * (ky*vz - kz*vy)
  wy = eye * (kz*vx - kx*vz)
  wz = eye * (kx*vy - ky*vx)

! write out the total energy in k space

  if (mod(jstep,ieout).eq.0) then
    call skewness(vx,kx,s,lx1,ly,lz,nx,ny,nz)
    tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
    tmp(1,:,:)=0.5_sp*tmp(1,:,:)
    ek = sum(tmp)
    tmp = (vx*conjg(wx)+vy*conjg(wy)+vz*conjg(wz)) &
        +(conjg(vx)*wx+conjg(vy)*wy+conjg(vz)*wz)
    tmp(1,:,:) = 0.5_sp*tmp(1,:,:)
    hek = sum(tmp)
    write(*,'(i8,6e12.4)')jstep,time,ek,hek,s,cs2
    write(25,'(i8,6e12.4)')jstep,time,ek,hek,s,cs2
    write(26,'(i8,6e12.4)')jstep,time,lm,mm
  end if

  ! output the field data and spectra, including the initial data when jstep=0.
  if (mod(jstep,itout).eq.0) then
    call output(vx,vy,vz,ii,lx1,ly,lz)
    ii = ii + 1
    
    ! tmp = u_i w_i
    tmp = (vx*conjg(wx)+vy*conjg(wy)+vz*conjg(wz)) &
         +(conjg(vx)*wx+conjg(vy)*wy+conjg(vz)*wz)
    tmp(1,:,:) = 0.5_sp*tmp(1,:,:)
    write(23,*) time          
    do i=1,lx
      hek=sum(tmp(:,:,:),mask=(abs(sqrt(k2)-i-0.5_sp*oneless).lt.0.5_sp))
      write(23,*) i,hek
    end do  
    
    ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
    ! at kx=0, where tmp=.5*q_ij(k)
    tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
    tmp(1,:,:)=0.5_sp*tmp(1,:,:)
    write(24,*)time
    do i=1,lx
      ek=sum(tmp(:,:,:),mask=(abs(sqrt(k2)-i-0.5_sp*oneless).lt.0.5_sp))
      write(24,*)i,ek
    end do

  end if

!========================== forcing term =================================  
  if (forcetype .ne. 0) call force_heli(fx,fy,fz,vx,vy,vz,wx,wy,wz,k2,kx,ky,kz,lx1,ly,lz,   &
                                            forcetype,fkmax,eps,eta)
!========================== calculate sgs stress ==============================
  if (mod(jstep,iupdate)==0) then 
          ! notice that at the first step jstep=0, so cs2 is calculated.
          update=.true.
  else
          update=.false.
  end if
  call dynsmag (vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,g_test,  &
                lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,mm,lm,cs2,  &
                update)
!=================== calculate convection term: lamb vector =======================
  call convec(vx,vy,vz,wx,wy,wz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
  !  (wx,wy,wz) stores lamb vector on return.

!========== combine convection term,divergence of sgs stress and forcing ============
  wx=wx-eye*(kx*t11+ky*t12+kz*t13)
  wy=wy-eye*(kx*t12+ky*t22+kz*t23)
  wz=wz-eye*(kx*t13+ky*t23+kz*t33)
  if (forcetype .ne. 0) then
    wx=wx+fx
    wy=wy+fy
    wz=wz+fz
  end if

  call symmetrize(wx,k2,lx1,ly,lz)
  call symmetrize(wy,k2,lx1,ly,lz)
  call symmetrize(wz,k2,lx1,ly,lz)

!======================== projecting ================================= 
  tmp = (kx*real(wx,sp) + ky*real(wy,sp) + kz*real(wz,sp))/k2
  wx = cmplx(real(wx,sp) - kx*tmp, aimag(wx))
  wy = cmplx(real(wy,sp) - ky*tmp, aimag(wy))
  wz = cmplx(real(wz,sp) - kz*tmp, aimag(wz))
  tmp = (kx*aimag(wx) + ky*aimag(wy) + kz*aimag(wz))/k2
  wx = cmplx(real(wx,sp), aimag(wx) - kx*tmp)
  wy = cmplx(real(wy,sp), aimag(wy) - ky*tmp)
  wz = cmplx(real(wz,sp), aimag(wz) - kz*tmp)

!======================= time stepping ======================================
  select case (istep)
  case (0)
    vx = vx * (1._sp-rnu*dt_h*k2) + dt_h*wx 
    vy = vy * (1._sp-rnu*dt_h*k2) + dt_h*wy 
    vz = vz * (1._sp-rnu*dt_h*k2) + dt_h*wz 

    time = time + dt_h
    jstep=itout ! this should work for itout that is greater than 1.
    call output(wx,wy,wz,0,lx1,ly,lz)

  case (1)
    vx = ox + dt*wx - rnu*dt*k2*vx 
    vy = oy + dt*wy - rnu*dt*k2*vy 
    vz = oz + dt*wz - rnu*dt*k2*vz 

    time = time + dt_h
    jstep = 0
    call input(ox,oy,oz,0,lx1,ly,lz)

  case (2:)
    vx = vx + dt_h*(3._sp*wx - k2_e*ox) 
    vy = vy + dt_h*(3._sp*wy - k2_e*oy) 
    vz = vz + dt_h*(3._sp*wz - k2_e*oz)

    vx = vx*k2_e
    vy = vy*k2_e
    vz = vz*k2_e

    ox = wx
    oy = wy
    oz = wz

    time = time + dt

  case default
    write(*,*) "istep wrong. istep = ", istep
  end select

  end do
!=================== end of loop ==============================

  close(23)
  close(24)
  close(25)
  close(26)

!  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz )
  deallocate(fx, fy, fz)
  deallocate(kx, ky, kz, k2, k2_e)
  deallocate(tmp)
  deallocate(t11,t12,t13,t22,t23,t33)

!  destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral
