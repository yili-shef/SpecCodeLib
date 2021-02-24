PROGRAM spectral_dns

  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx, ny, nz, iseed, nstep, itout, ieout, idp
  REAL(DP) :: dt, rnu, time, eps, fkmax, akp, u0
  LOGICAL  :: new, forced
  ! Parameters should be input from 'parameter_dns.d' initially.

  REAL(DP) :: delta
  INTEGER  :: nxb, nyb, nzb, lx, ly, lz, lxb, lyb, lzb, lx1, lxb1
  
  COMPLEX(DP), ALLOCATABLE, DIMENSION (:,:,:) ::  vx,vy,vz
  COMPLEX(DP), ALLOCATABLE, DIMENSION (:,:,:) ::  wx,wy,wz
  COMPLEX(DP), ALLOCATABLE, DIMENSION (:,:,:) ::  ox,oy,oz
  COMPLEX(DP), ALLOCATABLE, DIMENSION (:,:,:) ::  fx,fy,fz
  REAL(DP),    ALLOCATABLE, DIMENSION (:,:,:) ::  kx,ky,kz,k2
  REAL(DP),    ALLOCATABLE, DIMENSION (:,:,:) ::  k2_e,tmp

  INTEGER  :: i, ii, iii, istep, jstep
  REAL(DP) :: S, dt_h, ek, hek

  OPEN(90,file='parameter_dns.d',status='unknown')
    READ(90,*) forced
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) nstep
    READ(90,*) itout
    READ(90,*) ieout
    READ(90,*) dt
    READ(90,*) rnu
    READ(90,*) time
    READ(90,*) new
    READ(90,*) idp
    READ(90,*) eps
    READ(90,*) fkmax
    READ(90,*) akp
    READ(90,*) u0
    READ(90,*) iseed
  CLOSE(90)
  
  WRITE(*,*) 'forced', forced
  WRITE(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  WRITE(*,*) 'nstep', nstep
  WRITE(*,*) 'itout', itout
  WRITE(*,*) 'ieout', ieout
  WRITE(*,*) 'dt', dt
  WRITE(*,*) 'rnu', rnu
  WRITE(*,*) 'time', time
  WRITE(*,*) 'new', new
  WRITE(*,*) 'idp', idp
  WRITE(*,*) 'eps', eps
  WRITE(*,*) 'fkmax', fkmax
  WRITE(*,*) 'akp', akp
  WRITE(*,*) 'u0', u0
  WRITE(*,*) 'iseed', iseed
  
  CALL fftwplan3d(nx,ny,nz)
  WRITE(*,*) 'after fftwplan3d'
  
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

  delta=pi/REAL(lx,DP)

  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(k2_e(lx1,ly,lz),tmp(lx1,ly,lz))
  ALLOCATE(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  ALLOCATE(fx(lx1,ly,lz),fy(lx1,ly,lz),fz(lx1,ly,lz))
  WRITE(*,*) 'after allocate'

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  k2_e = EXP(-k2*dt*rnu)
  k2_e(1,1,1) = 1.0_DP
  WRITE(*,*) 'after wavenumber'

  IF (new) THEN
! Generate initial condition -------
    CALL initialize (vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    WRITE(*,*) 'after initialize'
  ELSE
! Rreading from DISK ------
    CALL input(vx,vy,vz,idp,lx1,ly,lz)
    PRINT*, 'initial velocity field readed.'
  ENDIF
  
  ox=vx
  oy=vy
  oz=vz

  dt_h=.5_DP*dt
  ii = idp
  istep=-1
  jstep=-1 
  ! istep counts the actual steps of time advancing.
  ! Because the step size for the first  two steps are half of dt,
  ! for the output purpose, the first two steps are counted as one
  ! step and the corresponding number of steps, with step size being dt, 
  ! is counted by jstep. Therefore jstep counts the number of steps on each
  ! of which dt has been advanced.
  ! Before the first step, istep=jstep=0; when
  ! istep=1 before the next actual time step, jstep is arbitrarily set (since 
  ! this step is not counted in jstep) to be itout+1 to avoid output 
  ! on this step (requiring itout>1). Then before further next step, istep=2,
  ! jstep=1.

!============== Begin the loop ====================================

  OPEN(24,file='./post/spectrum.data',status='unknown')
  OPEN(25,file='./post/ener_time.data',status='unknown')
  
  WRITE(*,*) 'starting the loop'
  DO WHILE (jstep .LE. nstep)

  istep = istep + 1
  jstep = jstep + 1

  IF(istep .EQ. 2) jstep=istep-1

! calculate vorticity
  wx = eye * (ky*vz - kz*vy)
  wy = eye * (kz*vx - kx*vz)
  wz = eye * (kx*vy - ky*vx)

! WRITE out the total energy in K space

  IF (istep .NE. 1 .AND. MOD(jstep,ieout).EQ.0) THEN
    CALL skewness(vx,kx,S,lx1,ly,lz,nx,ny,nz)
    tmp = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
    tmp(1,:,:)=0.5_DP*tmp(1,:,:)
    ek = SUM(tmp)
    WRITE(*,*)jstep,time,ek,S
    WRITE(25,*)jstep,time,ek,S
  END IF

  ! Output the field data and spectra, including the initial data when jstep=0.
  IF (istep.NE.1 .AND. MOD(jstep,itout).EQ.0) THEN
    CALL output(vx,vy,vz,ii,lx1,ly,lz)
    ii = ii + 1
    
    ! tmp = Q_ij(k), Q_ij(k) is the energy density spectral tensor, except
    ! at kx=0, where tmp=.5*Q_ij(k)
    tmp = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
    tmp(1,:,:)=0.5_DP*tmp(1,:,:)
    WRITE(24,*)time
    DO i=1,lx
      ek=SUM(tmp(:,:,:),mask=(ABS(SQRT(k2)-i).LT.0.5_DP))
      WRITE(24,*)i,ek
    END DO

  END IF

!========================== Forcing term =================================  
  IF (forced) CALL force(fx,fy,fz,vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,fkmax,eps)
!=================== Calculate convection term: Lamb vector =======================
  CALL convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)

!========== Combine convection term,divergence of SGS stress and forcing ============
  IF (forced) THEN
    wx=wx+fx
    wy=wy+fy
    wz=wz+fz
  END IF

  CALL symmetrize(wx,k2,lx1,ly,lz)
  CALL symmetrize(wy,k2,lx1,ly,lz)
  CALL symmetrize(wz,k2,lx1,ly,lz)

!======================== Projecting ================================= 
  tmp = (kx*REAL(wx,DP) + ky*REAL(wy,DP) + kz*REAL(wz,DP))/k2
  wx = CMPLX(REAL(wx,DP) - kx*tmp, AIMAG(wx))
  wy = CMPLX(REAL(wy,DP) - ky*tmp, AIMAG(wy))
  wz = CMPLX(REAL(wz,DP) - kz*tmp, AIMAG(wz))
  tmp = (kx*AIMAG(wx) + ky*AIMAG(wy) + kz*AIMAG(wz))/k2
  wx = CMPLX(REAL(wx,DP), AIMAG(wx) - kx*tmp)
  wy = CMPLX(REAL(wy,DP), AIMAG(wy) - ky*tmp)
  wz = CMPLX(REAL(wz,DP), AIMAG(wz) - kz*tmp)

!======================= Time stepping ======================================
  SELECT CASE (istep)
  CASE (0)
    vx = vx * (1._DP-rnu*dt_h*k2) + dt_h*wx 
    vy = vy * (1._DP-rnu*dt_h*k2) + dt_h*wy 
    vz = vz * (1._DP-rnu*dt_h*k2) + dt_h*wz 

    time = time + dt_h
    CALL output(wx,wy,wz,0,lx1,ly,lz)

  CASE (1)
    vx = ox + dt*wx - rnu*dt*k2*vx 
    vy = oy + dt*wy - rnu*dt*k2*vy 
    vz = oz + dt*wz - rnu*dt*k2*vz 

    time = time + dt_h
    CALL input(ox,oy,oz,0,lx1,ly,lz)

  CASE (2:)
    vx = vx + dt_h*(3._DP*wx - k2_e*ox) 
    vy = vy + dt_h*(3._DP*wy - k2_e*oy) 
    vz = vz + dt_h*(3._DP*wz - k2_e*oz)

    vx = vx*k2_e
    vy = vy*k2_e
    vz = vz*k2_e

    ox = wx
    oy = wy
    oz = wz

    time = time + dt

  CASE DEFAULT
    WRITE(*,*) "istep wrong. istep = ", istep
  END SELECT

  END DO
!=================== End of loop ==============================

  CLOSE(24)
  CLOSE(25)

!  Deallocate arrays
  DEALLOCATE(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  DEALLOCATE(fx, fy, fz)
  DEALLOCATE(kx, ky, kz, k2, k2_e)
  DEALLOCATE(tmp)

!  Destroy plans
  CALL destroyplan3d

  WRITE(*,*)'finished '
  STOP

END PROGRAM spectral_dns
