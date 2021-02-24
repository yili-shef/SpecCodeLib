PROGRAM spectral_les

  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: model, nx, ny, nz, iseed, ieout, idp
  REAL(SP) :: timemax, beta, rnu, dtw, time, eps, fkmax, akp, u0
  LOGICAL  :: new, forced
  ! Parameters should be input from 'parameter_les.d' initially.

  REAL(SP) :: delta
  INTEGER  :: lx, ly, lz, lx1, lxb, lyb, lzb, lxb1, nxb, nyb, nzb
  
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  vx,vy,vz,p
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  ox,oy,oz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  fx,fy,fz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  t11,t12,t13,t22,t23,t33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) ::  kx,ky,kz,k2
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) ::  k2_e,tmp

  INTEGER  :: i, ii, iii, istep
  REAL(SP) :: S, dtt, ek, hek, dt, dto, betadx, umax, twrite

  OPEN(90,file='parameter_les.d',status='unknown')
    READ(90,*) model
    READ(90,*) forced
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) timemax
    READ(90,*) ieout
    READ(90,*) beta
    READ(90,*) rnu
    READ(90,*) dtw
    READ(90,*) time
    READ(90,*) new
    READ(90,*) idp
    READ(90,*) eps
    READ(90,*) fkmax
    READ(90,*) akp
    READ(90,*) u0
    READ(90,*) iseed
  CLOSE(90)
  
  WRITE(*,*) 'model', model
  WRITE(*,*) 'forced', forced
  WRITE(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  WRITE(*,*) 'timemax', timemax
  WRITE(*,*) 'ieout', ieout
  WRITE(*,*) 'beta', beta
  WRITE(*,*) 'rnu', rnu
  WRITE(*,*) 'dtw', dtw
  WRITE(*,*) 'time', time
  WRITE(*,*) 'new', new
  WRITE(*,*) 'idp', idp
  WRITE(*,*) 'eps', eps
  WRITE(*,*) 'fkmax', fkmax
  WRITE(*,*) 'akp', akp
  WRITE(*,*) 'u0', u0
  WRITE(*,*) 'iseed', iseed
  WRITE(*,*) 'smallest', smallest
  WRITE(*,*) 'oneless', oneless
  
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

  delta=pi/REAL(lx,SP)
  betadx=beta*delta

  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(k2_e(lx1,ly,lz),tmp(lx1,ly,lz),p(lx1,ly,lz))
  ALLOCATE(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  ALLOCATE(fx(lx1,ly,lz),fy(lx1,ly,lz),fz(lx1,ly,lz))
  ALLOCATE(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  ALLOCATE(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  WRITE(*,*) 'after allocate'

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  WRITE(*,*) 'after wavenumber'

  IF (new) THEN
    ! Generate initial condition 
    CALL initialize (vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    WRITE(*,*) 'after initialize'
  ELSE
    ! Reading from disk
    CALL input(vx,vy,vz,idp,lx1,ly,lz)
    WRITE(*,*) 'initial velocity field readed.'
  ENDIF
  
  ox=vx
  oy=vy
  oz=vz

  dto=timemax
  twrite=time
  ii = idp
  istep=-1

!-------------------- Begin the loop ---------------------------

  OPEN(24,file='./post/spectrum.data',status='unknown')
  OPEN(25,file='./post/ener_time.data',status='unknown')
  
  WRITE(*,*) 'starting the loop'
  DO WHILE (time .LE. timemax)

  istep = istep + 1

  ! Calculate vorticity
  wx = eye * (ky*vz - kz*vy)
  wy = eye * (kz*vx - kx*vz)
  wz = eye * (kx*vy - ky*vx)

  ! Calculate pressure field.
  CALL pressure(vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,p)

  ! Write out the total energy in K space
  IF (MOD(istep,ieout).EQ.0) THEN
    CALL skewness(vx,kx,S,lx1,ly,lz,nx,ny,nz)
    tmp = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
    tmp(1,:,:)=0.5_SP*tmp(1,:,:)
    ek = SUM(tmp)
    WRITE( *,*)istep,time,dto,ek,S
    WRITE(25,*)istep,time,dto,ek,S
  END IF

  ! Output the field data and spectra, including the initial data when jstep=0.
  IF (ABS(time-twrite) .LE. smallest) THEN
    CALL output(vx,vy,vz,ii,lx1,ly,lz)
    CALL output_p(p,ii,lx1,ly,lz)
    ii = ii + 1
    twrite=twrite+dtw
    
    ! tmp = Q_ij(k), Q_ij(k) is the energy density spectral tensor, except
    ! at kx=0, where tmp=.5*Q_ij(k)
    tmp = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
    tmp(1,:,:)=0.5_SP*tmp(1,:,:)
    WRITE(24,*)time
    DO i=1,lx
      ek=SUM(tmp(:,:,:),mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
      WRITE(24,*)i,ek
    END DO

  END IF

  ! Subgrid scale model
!  CALL smag(vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz, &
!            lxb1,lyb,lzb,nxb,nyb,nzb,delta)

  CALL sgsm_p(p,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,delta)

  ! Convection term, finding umax for selecting dt
  CALL convec(vx,vy,vz,wx,wy,wz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,umax)

  ! Combining subgrid-scale stress and the convection term
  wx=wx-eye*(kx*t11+ky*t12+kz*t13)
  wy=wy-eye*(kx*t12+ky*t22+kz*t23)
  wz=wz-eye*(kx*t13+ky*t23+kz*t33)

  ! Optional force term
  IF (forced) THEN
    CALL force(fx,fy,fz,vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,fkmax,eps)
    wx=wx+fx
    wy=wy+fy
    wz=wz+fz
  END IF

  ! Nullifying Nyquist wavenumbers etc
  CALL symmetrize(wx,k2,lx1,ly,lz)
  CALL symmetrize(wy,k2,lx1,ly,lz)
  CALL symmetrize(wz,k2,lx1,ly,lz)

  ! Projecting wx wy wz to the k-normal plane
  tmp = (kx*REAL(wx,SP) + ky*REAL(wy,SP) + kz*REAL(wz,SP))/k2
  wx = CMPLX(REAL(wx,SP) - kx*tmp, AIMAG(wx))
  wy = CMPLX(REAL(wy,SP) - ky*tmp, AIMAG(wy))
  wz = CMPLX(REAL(wz,SP) - kz*tmp, AIMAG(wz))
  tmp = (kx*AIMAG(wx) + ky*AIMAG(wy) + kz*AIMAG(wz))/k2
  wx = CMPLX(REAL(wx,SP), AIMAG(wx) - kx*tmp)
  wy = CMPLX(REAL(wy,SP), AIMAG(wy) - ky*tmp)
  wz = CMPLX(REAL(wz,SP), AIMAG(wz) - kz*tmp)

  ! Selecting dt
  dt=MIN(betadx/umax,3._SP*dto)
  IF (dt .GE. twrite-time) THEN
          dt=twrite-time
  ELSE IF (twrite-time-dt .LT. dt/5._SP) THEN
          dt=.5_SP*(twrite-time)
  END IF

  ! Stepping forward
  SELECT CASE (istep)
  CASE (0)
    vx = vx * (1._SP-rnu*dt*k2) + dt*wx 
    vy = vy * (1._SP-rnu*dt*k2) + dt*wy 
    vz = vz * (1._SP-rnu*dt*k2) + dt*wz 

    time = time + dt
    dto=dt
    CALL output(wx,wy,wz,0,lx1,ly,lz)

  CASE (1)
    dtt=dt/dto
    vx = .5*(1.+dtt*dtt)*ox + .5*(1.-dtt*dtt)*vx + .5*dto*(1.+dtt)**2*(wx - rnu*k2*vx) 
    vy = .5*(1.+dtt*dtt)*oy + .5*(1.-dtt*dtt)*vy + .5*dto*(1.+dtt)**2*(wy - rnu*k2*vy) 
    vz = .5*(1.+dtt*dtt)*oz + .5*(1.-dtt*dtt)*vz + .5*dto*(1.+dtt)**2*(wz - rnu*k2*vz) 

    time = time + dt
    dto=dt+dto
    CALL input(ox,oy,oz,0,lx1,ly,lz)

  CASE (2:)
    dtt=dt/dto
    k2_e=EXP(-rnu*k2*dto)
    k2_e(1,1,1)=1.
    vx = vx + (1.+.5*dtt)*dt*wx-.5*dtt*dt*ox*k2_e
    vy = vy + (1.+.5*dtt)*dt*wy-.5*dtt*dt*oy*k2_e
    vz = vz + (1.+.5*dtt)*dt*wz-.5*dtt*dt*oz*k2_e

    k2_e=EXP(-rnu*k2*dt)
    k2_e(1,1,1)=1.
    vx = vx*k2_e
    vy = vy*k2_e
    vz = vz*k2_e

    ox = wx
    oy = wy
    oz = wz

    time = time + dt
    dto=dt

  CASE DEFAULT
    WRITE(*,*) 'istep wrong. istep = ', istep
    STOP
  END SELECT

  END DO
  WRITE(*,*) 'loop ended.'
!-------------------- End of loop ----------------------

  CLOSE(24)
  CLOSE(25)

!  Deallocate arrays
  DEALLOCATE(vx, vy, vz, wx, wy, wz, ox, oy, oz, p)
  DEALLOCATE(fx, fy, fz, t11, t12, t13, t22, t23, t33)
  DEALLOCATE(kx, ky, kz, k2, k2_e)
  DEALLOCATE(tmp)

!  Destroy plans
  CALL destroyplan3d

  WRITE(*,*) 'finished '
  STOP

END PROGRAM spectral_les
