PROGRAM spectral

  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx, ny, nz, iseed, ieout, idp
  REAL(SP) :: timemax, beta, dtw, rnu, time, eps, fk1, fk2, akp, u0, hrel, Ro
  LOGICAL  :: new, forced
  ! Parameters should be input from 'parameter_rot_dns.d' initially.

  INTEGER  :: lx, ly, lz, lx1
  
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  vx,vy,vz,ap,am,ef
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  ox,oy,oz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) ::  hpx,hpy,hpz
  real(sp),    allocatable, dimension (:,:,:) ::  k2,tmp
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     ::  kx,ky,kz

  INTEGER  :: i, ii, iii, istep, ifile, q
  REAL(SP) :: delta, Om, S, dtt, dt, ek, dto, betadx, umax, twrite

  OPEN(90,file='parameter_rot_dns.d',status='unknown')
    READ(90,*) forced
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
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
    READ(90,*) Ro
  CLOSE(90)
  
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
  WRITE(*,*) 'Ro', Ro

  Om=u0*akp/(2._SP*Ro*Pi)
  
  CALL fftwplan3d(nx,ny,nz)
  WRITE(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  delta=pi/REAL(lx,SP)
  betadx=beta*delta

  ALLOCATE(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  ALLOCATE(efp(lx1,ly,lz),efm(lx1,ly,lz))
  ALLOCATE(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  ALLOCATE(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  ALLOCATE(ap(lx1,ly,lz),am(lx1,ly,lz))
  WRITE(*,*) 'after allocate'

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  do i=1,lx1
  do ii=1,ly
    efp(i,ii,:) = EXP(-(k2(i,ii,:)/nx/nx)**q*dt*rnuqk &
                 +2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
    emp(i,ii,:) = EXP(-(k2(i,ii,:)/nx/nx)**q*dt*rnuqk &
                 -2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
  enddo
  enddo 
  WRITE(*,*) 'after wavenumber'

  CALL heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  WRITE(*,*) 'after heliwave'


  IF (new) THEN
! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    write(*,*) 'after initialize'
  ELSE
! Reading from DISK ------
    CALL input(vx,vy,vz,idp,lx1,ly,lz)
    PRINT*, 'initial velocity field readed.'
  ENDIF
  
  ox=vx
  oy=vy
  oz=vz
  ap=vx*CONJG(hpx)+vy*CONJG(hpy)+vz*CONJG(hpz)
  am=vx*hpx+vy*hpy+vz*hpz

  ! specify umax empirically as three times the rms value, in reference to
  ! Gaussian distribution
  umax=3.*sqrt(3.)*u0
  
  dto=timemax
  twrite=time
  ifile = idp
  istep=-1

!============== Begin the loop ====================================

  OPEN(24,file='./post/spectrum.data',status='unknown')
  OPEN(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

    istep = istep + 1
  
    !WRITE out the total energy in K space
 
    IF (MOD(istep,ieout).EQ.0) THEN
      CALL skewness(vx,kx,S,lx1,ly,lz,nx,ny,nz)
      ef = vx*CONJG(vx) + vy*CONJG(vy) + vz*CONJG(vz)
      ef(1,:,:)=0.5_SP*ef(1,:,:)
      ek = SUM(real(ef))
      WRITE(*,'(I6,6E12.4)')istep,time,ek,S
      WRITE(25,'(I6,6E12.4)')istep,time,ek,S
 
      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.*sqrt(2.*ek)
    END IF
 
    ! Output the field data and spectra, including the initial data when time=0.
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
  
 
    !=================== Calculate convection term: Lamb vector =======================
    CALL convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========================== Forcing term =================================  
    if (forced) call force_rot(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fk1,fk2,eps)
 
    CALL symmetrize(wx,k2,lx1,ly,lz)
    CALL symmetrize(wy,k2,lx1,ly,lz)
    CALL symmetrize(wz,k2,lx1,ly,lz)
 
    !======================== Projecting ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
 
    !======================= Time stepping ======================================
    SELECT CASE (istep)
    CASE (0)
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)**q*dt*rnu+2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
      enddo
      enddo 
      ap=ap*ef+dt*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dt*rnu-2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
      enddo
      enddo 
      am=am*ef+dt*(wx*hpx+wy*hpy+wz*hpz)
      
      time = time + dt
      dto=dt
      CALL output(wx,wy,wz,0,lx1,ly,lz)
 
    CASE (1)
      dtt=dt/dto
 
      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dto*rnu+2._SP*eye*dto*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
      end do
      end do 
 
      ! ox,oy,oz are the initial velocity components now
      ap = S*(ox*CONJG(hpx)+oy*CONJG(hpy)+oz*CONJG(hpz))*ef + ek*ap + &
           delta*(wx*CONJG(hpx)+wy*CONJG(hpy)+wz*CONJG(hpz))
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dt*rnu+2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
      end do
      end do 
      ap = ap*ef
 
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dto*rnu-2._SP*eye*dto*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
      end do
      end do 
 
      ! ox,oy,oz are the initial velocity components now
      am = S*(ox*hpx+oy*hpy+oz*hpz)*ef + ek*am + delta*(wx*hpx+wy*hpy+wz*hpz)
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dt*rnu-2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
      end do
      end do 
      am = am*ef
 
      time = time + dt
      dto=dt+dto
      CALL input(ox,oy,oz,0,lx1,ly,lz)
 
    CASE (2:)
      dtt=dt/dto
      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dto*rnu+2._SP*eye*dto*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
      end do
      end do 
 
      vx = S*wx + ek*ox*ef
      vy = S*wy + ek*oy*ef
      vz = S*wz + ek*oz*ef
      ap = ap + vx*CONJG(hpx)+vy*CONJG(hpy)+vz*CONJG(hpz)
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dt*rnu+2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Positive helical wave
      end do
      end do 
      ap=ap*ef
 
      
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dto*rnu-2._SP*eye*dto*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
      end do
      end do 
 
      vx = S*wx + ek*ox*ef
      vy = S*wy + ek*oy*ef
      vz = S*wz + ek*oz*ef
      am = am + vx*hpx+vy*hpy+vz*hpz
 
      do i=1,lx1
      do ii=1,ly
      ef(i,ii,:) = EXP(-k2(i,ii,:)*dt*rnu-2._SP*eye*dt*kz*Om/(SQRT(k2(i,ii,:))+smallest))  ! Negative helical wave
      end do
      end do 
      am=am*ef
 
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      dto=dt
 
    CASE DEFAULT
      WRITE(*,*) "istep wrong. istep = ", istep
    END SELECT
    vx=ap*hpx+am*conjg(hpx)
    vy=ap*hpy+am*conjg(hpy)
    vz=ap*hpz+am*conjg(hpz)

  END DO
!=================== End of loop ==============================

  CLOSE(24)
  CLOSE(25)

!  Deallocate arrays
  DEALLOCATE(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  DEALLOCATE(ap, am, hpx, hpy, hpz)
  DEALLOCATE(kx, ky, kz, k2, ef)

!  Destroy plans
  CALL destroyplan3d

  WRITE(*,*)'finished '
  STOP

END PROGRAM spectral
