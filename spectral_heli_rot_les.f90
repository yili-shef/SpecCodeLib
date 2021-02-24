PROGRAM spectral

  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx, ny, nz, model, forcetype, iseed, nstep, itout, ieout, idp, iupdate
  REAL(SP) :: dt, rnu, time, eps, eta, fkmax, akp, u0, hrel, Ro, delta_ratio
  LOGICAL  :: new
  ! Parameters should be input from 'parameter_heli_rot_les.d' initially.

  REAL(SP) :: delta, Om, delta_test,k_test,c1,c2
  INTEGER  :: nxb, nyb, nzb, lx, ly, lz, lxb, lyb, lzb, lx1, lxb1
  LOGICAL  :: update
  
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  vx,vy,vz,ap,am,efp,efm
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  ox,oy,oz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  hpx,hpy,hpz,hmx,hmy,hmz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  fx,fy,fz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  t11,t12,t13,t22,t23,t33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) ::  kx,ky,kz,k2,tmp,G_test

  INTEGER  :: i, ii, iii, istep, jstep
  REAL(SP) :: S, dt_h, ek, hek, LM, LP, MM, PP, PM

  OPEN(90,file='parameter_heli_rot_les.d',status='unknown')
    READ(90,*) model
    READ(90,*) forcetype
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) iseed
    READ(90,*) nstep
    READ(90,*) itout
    READ(90,*) ieout
    READ(90,*) dt
    READ(90,*) rnu
    READ(90,*) time
    READ(90,*) new
    READ(90,*) idp
    READ(90,*) eps
    READ(90,*) eta
    READ(90,*) fkmax
    READ(90,*) akp
    READ(90,*) u0
    READ(90,*) hrel
    READ(90,*) Ro
    READ(90,*) iupdate
    READ(90,*) delta_ratio
  CLOSE(90)
  
  WRITE(*,*) 'model', model
  WRITE(*,*) 'forcetype', forcetype
  WRITE(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  WRITE(*,*) 'iseed', iseed
  WRITE(*,*) 'nstep', nstep
  WRITE(*,*) 'itout', itout
  WRITE(*,*) 'ieout', ieout
  WRITE(*,*) 'dt', dt
  WRITE(*,*) 'rnu', rnu
  WRITE(*,*) 'time', time
  WRITE(*,*) 'new', new
  WRITE(*,*) 'idp', idp
  WRITE(*,*) 'eps', eps
  WRITE(*,*) 'eta', eta
  WRITE(*,*) 'fkmax', fkmax
  WRITE(*,*) 'akp', akp
  WRITE(*,*) 'u0', u0
  WRITE(*,*) 'hrel', hrel
  WRITE(*,*) 'Ro', Ro
  WRITE(*,*) 'iupdate', iupdate
  WRITE(*,*) 'delta_ratio', delta_ratio

  Om=u0*akp/(2._SP*Ro*Pi)
  
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
  delta_test=delta_ratio*delta
  k_test=lx/delta_ratio

  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(efp(lx1,ly,lz),efm(lx1,ly,lz),tmp(lx1,ly,lz))
  ALLOCATE(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  ALLOCATE(fx(lx1,ly,lz),fy(lx1,ly,lz),fz(lx1,ly,lz))
  ALLOCATE(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  ALLOCATE(hmx(lx1,ly,lz),hmy(lx1,ly,lz),hmz(lx1,ly,lz))
  ALLOCATE(ap(lx1,ly,lz),am(lx1,ly,lz),G_test(lx1,ly,lz))
  ALLOCATE(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  ALLOCATE(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  WRITE(*,*) 'after allocate'

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  efp = EXP(-k2*dt*rnu+2._SP*eye*dt*kz*Om/(SQRT(k2)+smallest))
  efm = EXP(-k2*dt*rnu-2._SP*eye*dt*kz*Om/(SQRT(k2)+smallest))
  efp(1,1,1) = 1.0_SP
  efm(1,1,1) = 1.0_SP
  WRITE(*,*) 'after wavenumber'

  WHERE(k2 .LE. k_test*k_test)
          G_test=1._SP
  ELSEWHERE
          G_test=0._SP
  ENDWHERE
  WRITE(*,*) 'after G_test'

  CALL heliwave(hpx,hpy,hpz,hmx,hmy,hmz,kx,ky,kz,k2,lx1,ly,lz)
  WRITE(*,*) 'after heliwave'


  IF (new) THEN
! Generate initial condition -------
    CALL initialize_heli(hpx,hpy,hpz,hmx,hmy,hmz,vx,vy,vz,k2,lx1,ly,lz,iseed,hrel,akp,u0)
    WRITE(*,*) 'after initialize_heli'
  ELSE
! Rreading from DISK ------
    CALL input(vx,vy,vz,idp,lx1,ly,lz)
    PRINT*, 'initial velocity field readed.'
  ENDIF
  
  ox=vx
  oy=vy
  oz=vz

  dt_h=.5_SP*dt
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

  OPEN(23,file='./post/helicity_spectrum.data',status='unknown') 
  OPEN(24,file='./post/spectrum.data',status='unknown')
  OPEN(25,file='./post/ener_time.data',status='unknown')
  
  WRITE(*,*) 'starting the loop'
  DO WHILE (jstep .LE. nstep)

  istep = istep + 1
  jstep = jstep + 1

! calculate vorticity
  wx = eye * (ky*vz - kz*vy)
  wy = eye * (kz*vx - kx*vz)
  wz = eye * (kx*vy - ky*vx)

! WRITE out the total energy in K space

  IF (MOD(jstep,ieout).EQ.0) THEN
    CALL skewness(vx,kx,S,lx1,ly,lz,nx,ny,nz)
    tmp = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
    tmp(1,:,:)=0.5_SP*tmp(1,:,:)
    ek = SUM(tmp)
    tmp = (vx*CONJG(wx)+vy*CONJG(wy)+vz*CONJG(wz)) &
        +(CONJG(vx)*wx+CONJG(vy)*wy+CONJG(vz)*wz)
    tmp(1,:,:) = 0.5_SP*tmp(1,:,:)
    hek = SUM(tmp)
    WRITE(*,'(I6,6E12.4)')jstep,time,ek,hek,S,c1,c2
    WRITE(25,'(I6,6E12.4)')jstep,time,ek,hek,S,c1,c2
  END IF

  ! Output the field data and spectra, including the initial data when jstep=0.
  IF (MOD(jstep,itout).EQ.0) THEN
    CALL output(vx,vy,vz,ii,lx1,ly,lz)
    ii = ii + 1
    
    ! tmp = u_i w_i
    tmp = (vx*CONJG(wx)+vy*CONJG(wy)+vz*CONJG(wz)) &
         +(CONJG(vx)*wx+CONJG(vy)*wy+CONJG(vz)*wz)
    tmp(1,:,:) = 0.5_SP*tmp(1,:,:)
    WRITE(23,*) time          
    DO i=1,lx
      hek=SUM(tmp(:,:,:),mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
      WRITE(23,*) i,hek
    END DO  
    
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

!========================== Forcing term =================================  
  IF (forcetype .NE. 0) CALL force_heli(fx,fy,fz,vx,vy,vz, wx,wy,wz,k2,kx,ky,kz,lx1,ly,lz,   &
                                            forcetype,fkmax,eps,eta)
!========================== Calculate SGS stress ==============================
  SELECT CASE (model)
  CASE (1) ! Smagorinsky model
    CALL smag(vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz, &
              lxb1,lyb,lzb,nxb,nyb,nzb,delta)
  CASE (2) ! Local helicity model
    CALL localheli(vx,vy,vz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                   lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
  CASE (3) ! Nonlocal helicity model, Not applicable in decaying turbulence.
    CALL nonlocalheli(vx,vy,vz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)

  CASE (4) 
    IF (MOD(jstep,iupdate)==0) THEN 
          ! Notice that at the first step jstep=0, so c1 and c2 are calculated.
          update=.true.
    ELSE
          update=.false.
    END IF
    CALL dynheli_a (vx,vy,vz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                    G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,&
                    PP,LM,LP,c1,c2,update)
  CASE (5) 
    IF (MOD(jstep,iupdate)==0) THEN 
          ! Notice that at the first step jstep=0, so c1 and c2 are calculated.
          update=.true.
    ELSE
          update=.false.
    END IF
    CALL dynheli_b (vx,vy,vz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                    G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,&
                    PP,LM,LP,c1,c2,update)
  CASE DEFAULT
    WRITE(*,*) "model # wrong, model = ", model
    STOP
  END SELECT

!=================== Calculate convection term: Lamb vector =======================
  CALL convec(vx,vy,vz,wx,wy,wz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
  !  (wx,wy,wz) stores Lamb vector on return.

!========== Combine convection term,divergence of SGS stress and forcing ============
  wx=wx-eye*(kx*t11+ky*t12+kz*t13)
  wy=wy-eye*(kx*t12+ky*t22+kz*t23)
  wz=wz-eye*(kx*t13+ky*t23+kz*t33)
  IF (forcetype .NE. 0) THEN
    wx=wx+fx
    wy=wy+fy
    wz=wz+fz
  END IF

  CALL symmetrize(wx,k2,lx1,ly,lz)
  CALL symmetrize(wy,k2,lx1,ly,lz)
  CALL symmetrize(wz,k2,lx1,ly,lz)

!======================== Projecting ================================= 
  tmp = (kx*REAL(wx,SP) + ky*REAL(wy,SP) + kz*REAL(wz,SP))/k2
  wx = CMPLX(REAL(wx,SP) - kx*tmp, AIMAG(wx))
  wy = CMPLX(REAL(wy,SP) - ky*tmp, AIMAG(wy))
  wz = CMPLX(REAL(wz,SP) - kz*tmp, AIMAG(wz))
  tmp = (kx*AIMAG(wx) + ky*AIMAG(wy) + kz*AIMAG(wz))/k2
  wx = CMPLX(REAL(wx,SP), AIMAG(wx) - kx*tmp)
  wy = CMPLX(REAL(wy,SP), AIMAG(wy) - ky*tmp)
  wz = CMPLX(REAL(wz,SP), AIMAG(wz) - kz*tmp)

!======================= Time stepping ======================================
  ap=vx*CONJG(hpx)+vy*CONJG(hpy)+vz*CONJG(hpz)
  am=vx*CONJG(hmx)+vy*CONJG(hmy)+vz*CONJG(hmz)
  SELECT CASE (istep)
  CASE (0)
    ap=ap*SQRT(efp)+dt_h*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))
    am=am*SQRT(efm)+dt_h*(wx*CONJG(hmx)+wy*CONJG(hmy)+wz*CONJG(hmz))
    
    time = time + dt_h
    jstep=itout ! This should work for itout that is greater than 1.
    CALL output(wx,wy,wz,0,lx1,ly,lz)

  CASE (1)
    ! ox,oy,oz are the initial velocity components now
    ap = (ox*CONJG(hpx)+oy*CONJG(hpy)+oz*CONJG(hpz))*efp + dt*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))*SQRT(efp)
    am = (ox*CONJG(hmx)+oy*CONJG(hmy)+oz*CONJG(hmz))*efm + dt*(wx*CONJG(hmx)+wy*CONJG(hmy)+wz*CONJG(hmz))*SQRT(efm)

    time = time + dt_h
    jstep = 0
    CALL input(ox,oy,oz,0,lx1,ly,lz)
    ! after the input, ox,oy,oz store the components of the rhs of the equation
    ! at previous time step: wx,wy,wz.

  CASE (2:)
    fx = 3._SP*wx - efp*ox
    fy = 3._SP*wy - efp*oy
    fz = 3._SP*wz - efp*oz
    ap = (ap + dt_h*(fx*CONJG(hpx)+fy*CONJG(hpy)+fz*CONJG(hpz)))*efp
    
    fx = 3._SP*wx - efm*ox
    fy = 3._SP*wy - efm*oy
    fz = 3._SP*wz - efm*oz
    am = (am + dt_h*(fx*CONJG(hmx)+fy*CONJG(hmy)+fz*CONJG(hmz)))*efm

    ox = wx
    oy = wy
    oz = wz

    time = time + dt

  CASE DEFAULT
    WRITE(*,*) "istep wrong. istep = ", istep
  END SELECT
  vx=ap*hpx+am*hmx
  vy=ap*hpy+am*hmy
  vz=ap*hpz+am*hmz

  END DO
!=================== End of loop ==============================

  CLOSE(23)
  CLOSE(24)
  CLOSE(25)

!  Deallocate arrays
  DEALLOCATE(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  DEALLOCATE(fx, fy, fz, ap, am, hpx, hpy, hpz, hmx, hmy, hmz)
  DEALLOCATE(kx, ky, kz, k2, efp, efm, G_test)
  DEALLOCATE(tmp)
  DEALLOCATE(t11,t12,t13,t22,t23,t33)

!  Destroy plans
  CALL destroyplan3d

  WRITE(*,*)'finished '
  STOP

END PROGRAM spectral
