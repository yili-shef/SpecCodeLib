program spectral
  use mconstant
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, nstep, itout, ieout, idp, q
  real(sp) :: dt, rnukq, time, eps, fk1, fk2, akp, u0, ro
  logical  :: new, forced
  ! parameters should be input from 'parameter_rot_dns.d' initially.

  real(sp) :: delta, om
  integer  :: lx, ly, lz, lx1, lx2
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am,efp
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  real(sp),    allocatable, dimension (:,:,:) ::  k2
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz, ek

  INTEGER  :: i, ii, iii, ll, istep, jstep, ifile
  REAL(SP) :: S, dt_h, rnu, ktrun, ktrun2

  OPEN(90,file='parameter_rot_dns_ft.d',status='unknown')
    READ(90,*) forced
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) nstep
    READ(90,*) itout
    READ(90,*) ieout
    READ(90,*) dt
    read(90,*) q
    READ(90,*) rnukq ! rnukq=rnu_q * k_max^2q  
    READ(90,*) time
    READ(90,*) new
    READ(90,*) idp
    READ(90,*) eps
    READ(90,*) fk1
    read(90,*) fk2
    READ(90,*) akp
    READ(90,*) u0
    READ(90,*) Ro
    READ(90,*) iseed
  CLOSE(90)
  
  WRITE(*,*) 'forced', forced
  WRITE(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  WRITE(*,*) 'nstep', nstep
  WRITE(*,*) 'itout', itout
  WRITE(*,*) 'ieout', ieout
  WRITE(*,*) 'dt', dt
  WRITE(*,*) 'rnukq', rnukq
  WRITE(*,*) 'time', time
  WRITE(*,*) 'new', new
  WRITE(*,*) 'idp', idp
  WRITE(*,*) 'eps', eps
  WRITE(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  WRITE(*,*) 'akp', akp
  WRITE(*,*) 'u0', u0
  WRITE(*,*) 'Ro', Ro
  write(*,*) 'q', q
  WRITE(*,*) 'iseed', iseed

  if (forced) then
    Om=eps**(1./3.)*(.5*(fk1+fk2))**(2./3.)/Ro   ! Ro is defined by forcing scales
  else
    Om=u0*akp/(2._SP*Ro*Pi)  ! Ro is defined by initial integral scales
  end if
  !Om=0. ! Turn off rotation
  
  call fftwplan3d(nx,ny,nz)
  WRITE(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz
  
  lx1=lx+1
  ktrun=real(ny,sp)/3.
  ktrun2=ktrun*ktrun

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz),ek(lx))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(efp(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  do iii=1,lz
  do ii=1,ly
  do i=1,lx1
    efp(i,ii,iii) = exp(-(k2(i,ii,iii)/ktrun2)**q*dt*rnukq &
                   +2._sp*eye*dt*kz(iii)*Om/(sqrt(k2(i,ii,iii))+smallest))  ! positive helical wave
  enddo
  enddo 
  enddo
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'


  IF (new) THEN
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)

    call dealiasing(vx,k2,lx1,ly,lz)
    call dealiasing(vy,k2,lx1,ly,lz)
    call dealiasing(vz,k2,lx1,ly,lz)

    write(*,*) 'after initialize'
  else
    ! Reading from DISK ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  ENDIF
  
  ox=vx
  oy=vy
  oz=vz
  ap=vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am=vx*hpx+vy*hpy+vz*hpz

  dt_h=.5_SP*dt
  ifile = idp
  istep=0
  jstep=0 

!============== Begin the loop ====================================

  OPEN(24,file='./post/spectrum.data',status='unknown')
  OPEN(25,file='./post/ener_time.data',status='unknown')
  
  WRITE(*,*) 'starting the loop'
  DO WHILE (jstep .LE. nstep)

    ! WRITE out the total energy in K space

    IF (MOD(jstep,ieout).EQ.0) THEN
      ! CALL skewness(vx,kx,S,nx,ny,nz)
      wx = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
      wx(1,:,:)=0.5_SP*wx(1,:,:)
      ktrun = SUM(wx) ! Total TKE
      WRITE(*,'(I6,6f12.4)') jstep,time,ktrun,S
      WRITE(25,'(I6,6f12.4)')jstep,time,ktrun,S
    END IF
 
    ! Output the field data and spectra, including the initial data when jstep=0.
    IF (MOD(jstep,itout).EQ.0) THEN
      CALL output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      
      ! wx = Q_ij(k), Q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where wx=.5*Q_ij(k)
      wx = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
      wx(1,:,:)=0.5_SP*wx(1,:,:)
      ek = 0.
      do iii = 1, lz
      do ii = 1, ly
      do i = 1, lx1
        ll = floor(sqrt(k2(i,ii,iii)) + .5)
        if ( ll .ge. 1 .and. ll .le. lx ) then
          ek(ll) = ek(ll) + real(wx(i,ii,iii))
        end if
      end do
      end do
      end do

      WRITE(24,*)time
      do i = 1, lx
        write(24,*)i,ek(i)
      end do
 
    END IF

    ! =========== Calculate vorticity =============================
    do iii=1,lz
    do ii=1,ly
    do i=1,lx1
      wx(i,ii,iii) = eye*(ky(ii)*vz(i,ii,iii) - kz(iii)*vy(i,ii,iii))
      wy(i,ii,iii) = eye*(kz(iii)*vx(i,ii,iii) - kx(i)*vz(i,ii,iii))
      wz(i,ii,iii) = eye*(kx(i)*vy(i,ii,iii) - ky(ii)*vx(i,ii,iii))
    end do
    end do
    end do
  
    !=================== Calculate convection term: Lamb vector =======================
    CALL convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    call symmetrize(wx,k2,lx1,ly,lz)
    call symmetrize(wy,k2,lx1,ly,lz)
    call symmetrize(wz,k2,lx1,ly,lz)

    call dealiasing(wx,k2,lx1,ly,lz)
    call dealiasing(wy,k2,lx1,ly,lz)
    call dealiasing(wz,k2,lx1,ly,lz)
 
    !======================== Projecting ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)

    !======================= Time stepping ======================================
    SELECT CASE (istep)
    CASE (0)
      ap=ap*SQRT(efp)+dt_h*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))
      am=am*SQRT(conjg(efp))+dt_h*(wx*hpx+wy*hpy+wz*hpz)
      
      time = time + dt_h
      istep=istep+1
      jstep=1 ! only half of first step, set to one to avoid duplicate output

      CALL output(wx,wy,wz,0,lx1,ly,lz)
 
    CASE (1)
      ! ox,oy,oz are the initial velocity components now
      ap = (ox*CONJG(hpx)+oy*CONJG(hpy)+oz*CONJG(hpz))*efp + dt*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))*SQRT(efp)
      am = (ox*hpx+oy*hpy+oz*hpz)*conjg(efp) + dt*(wx*hpx+wy*hpy+wz*hpz)*SQRT(conjg(efp))
 
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! the actual first step.

      CALL input(ox,oy,oz,0,lx1,ly,lz)
      ! after the input, ox,oy,oz store the components of the rhs of the equation
      ! at previous time step: wx,wy,wz.
 
    CASE (2:)
      vx = 3._SP*wx - efp*ox
      vy = 3._SP*wy - efp*oy
      vz = 3._SP*wz - efp*oz
      ap = (ap + dt_h*(vx*CONJG(hpx)+vy*CONJG(hpy)+vz*CONJG(hpz)))*efp
      
      vx = 3._SP*wx - conjg(efp)*ox
      vy = 3._SP*wy - conjg(efp)*oy
      vz = 3._SP*wz - conjg(efp)*oz
      am = (am + dt_h*(vx*hpx+vy*hpy+vz*hpz))*conjg(efp)
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      istep = istep + 1
      jstep = jstep + 1
 
    CASE DEFAULT
      WRITE(*,*) "istep wrong. istep = ", istep
    END SELECT
    vx=ap*hpx+am*conjg(hpx)
    vy=ap*hpy+am*conjg(hpy)
    vz=ap*hpz+am*conjg(hpz)

    !========================== Forcing term =================================  
    if (forced) call force_rescale(vx,vy,vz,k2,lx1,ly,lz)

  END DO
!=================== End of loop ==============================

  CLOSE(24)
  CLOSE(25)

!  Deallocate arrays
  DEALLOCATE(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  DEALLOCATE(ap, am, hpx, hpy, hpz)
  DEALLOCATE(kx, ky, kz, k2, efp)

!  Destroy plans
  CALL destroyplan3d

  WRITE(*,*)'finished '
  STOP

end program spectral
