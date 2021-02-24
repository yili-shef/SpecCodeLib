module mfftw3plan
  use mconstant
  use mfftw3param

  integer*8 :: c2rvx, c2rvy, c2rvz, c2rwx, c2rwy, c2rwz, r2cwx, r2cwy, r2cwz 
  integer*8 :: r2cvx, r2cvy, r2cvz, c2rwfi

contains

  subroutine destroyplans

    call dfftw_destroy_plan(c2rvx)
    call dfftw_destroy_plan(c2rvy)
    call dfftw_destroy_plan(c2rvz)

    call dfftw_destroy_plan(c2rwx)
    call dfftw_destroy_plan(c2rwy)
    call dfftw_destroy_plan(c2rwz)

    call dfftw_destroy_plan(r2cwx)
    call dfftw_destroy_plan(r2cwy)
    call dfftw_destroy_plan(r2cwz)

    call dfftw_destroy_plan(r2cvx)
    call dfftw_destroy_plan(r2cvy)
    call dfftw_destroy_plan(r2cvz)

    call dfftw_destroy_plan(c2rwfi)
  end subroutine destroyplans

end module mfftw3plan

program spectral

  use mconstant
  use mwavenumber
  use mheliwave
  use moutput
  use minput
  use mfftw3param
  use mfftw3plan
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, dtw, rnu, time, eps, fk1, fk2, akp, u0, lmbd
  logical  :: new, forced
  ! parameters should be input from 'parameter_rot_dns.d' initially.

  integer  :: lx, ly, lz, lx1
  
  integer, parameter :: numth = 4

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  oap,oam
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  complex(sp), allocatable, dimension (:,:,:) ::  fi,ofi,wfi,ffi
  real(sp),    allocatable, dimension (:,:,:) ::  k2,kx,ky,kz
  real(sp),    allocatable, dimension (:)     ::  ekk

  integer  :: istep, ifile
  real(sp) :: delta, om, s, dtt, dt, ek, dto, betadx, umax, fcoe
  real(sp) :: twrite, kinet, kcut, kcut2, fk12, fk22, const

  real(sp), parameter :: pr = 1.0_sp

  open(90,file='parameter_rot_dns_scalar.d',status='unknown')
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
    read(90,*) om
    read(90,*) lmbd
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
  write(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'om', om
  write(*,*) 'mean scalar gradient', lmbd
  write(*,*) 'iseed', iseed

  lx=nx/2
  ly=ny
  lz=nz
  lx1=lx+1

  const = 1._sp/(nx*ny*nz)

  kcut = 2._sp*lx/3
  delta = pi / real(kcut,sp)
  betadx = beta*delta

  kcut2 = kcut * kcut

  fk12 = fk1 * fk1
  fk22 = fk2 * fk2

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(oap(lx1,ly,lz),oam(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz),ffi(lx1,ly,lz))
  allocate(fi(lx1,ly,lz),ofi(lx1,ly,lz),wfi(lx1,ly,lz))
  allocate(ekk(lx))
  write(*,*) 'after allocate'

  call dfftw_plan_dft_c2r_3d(c2rvx,nx,ny,nz,vx,vx,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_3d(c2rvy,nx,ny,nz,vy,vy,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_3d(c2rvz,nx,ny,nz,vz,vz,FFTW_MEASURE)

  call dfftw_plan_dft_c2r_3d(c2rwx,nx,ny,nz,wx,wx,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_3d(c2rwy,nx,ny,nz,wy,wy,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_3d(c2rwz,nx,ny,nz,wz,wz,FFTW_MEASURE)
 
  call dfftw_plan_dft_r2c_3d(r2cwx,nx,ny,nz,wx,wx,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_3d(r2cwy,nx,ny,nz,wy,wy,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_3d(r2cwz,nx,ny,nz,wz,wz,FFTW_MEASURE)

  call dfftw_plan_dft_r2c_3d(r2cvx,nx,ny,nz,vx,vx,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_3d(r2cvy,nx,ny,nz,vy,vy,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_3d(r2cvz,nx,ny,nz,vz,vz,FFTW_MEASURE)

  call dfftw_plan_dft_c2r_3d(c2rwfi,nx,ny,nz,wfi,wfi,FFTW_MEASURE)
 
  write(*,*) 'after creating plans'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'

  if (new) then
    ! generate initial condition -------
    call initrotphi(vx,vy,vz,fi,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'after initialize'
  else
    ! reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    call inputphi(fi,idp,lx1,ly,lz)
    print*, 'initial velocity field read.'
  endif
  
  ap = vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am = vx*hpx+vy*hpy+vz*hpz
  oap = ap
  oam = am
  ofi = fi

  ! specify umax empirically as three times the rms value, in reference to
  ! gaussian distribution
  umax = 3._sp*sqrt(3._sp)*u0
  
  dto = timemax
  twrite = time
  ifile = idp
  istep = -1
  s = 0._sp

!============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

    istep = istep + 1
  
    !write out the total energy in k space
 
    if (mod(istep,ieout).eq.0) then

      !call skewness(vz,kz,s,lx1,ly,lz,nx,ny,nz)
      wx = ap*conjg(ap) + am*conjg(am) 
      wx(1,:,:) = 0.5_sp*wx(1,:,:)
      ek = sum(real(wx))
      write( *,'(i6,6e15.6)')istep, time, dto, ek, s
      write(25,'(i6,6e15.6)')istep, time, dto, ek, s
 
      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.5_sp*sqrt(2._sp*ek)
    end if
 
    ! output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. myeps) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      call outputphi(fi,ifile,lx1,ly,lz)
      ifile = ifile + 1
      twrite=twrite+dtw
    end if
  
    ! =========== coefficient for forcing ==============
    if (forced) then 

      where (k2 .lt. fk22 .and. k2 .ge. fk12)
        wx = ap*conjg(ap) + am*conjg(am) 
      elsewhere
        wx = 0._sp
      end where
      wx(1,:,:) = .5_sp * wx(1,:,:)
      kinet=sum(real(wx,sp))
  
      fcoe = eps/(2._sp*kinet+mytiny)
    end if

    ffi = - lmbd * vz

    ! =========== calculate vorticity =============================

    wx = eye * ( ky * vz - kz * vy )
    wy = eye * ( kz * vx - kx * vz )
    wz = eye * ( kx * vy - ky * vx )

    !=================== calculate convection term: lamb vector =======================

    call dfftw_execute_dft_c2r(c2rvx, vx, vx)
    call dfftw_execute_dft_c2r(c2rvy, vy, vy)
    call dfftw_execute_dft_c2r(c2rvz, vz, vz)
    
    call dfftw_execute_dft_c2r(c2rwx, wx, wx)
    call dfftw_execute_dft_c2r(c2rwy, wy, wy)
    call dfftw_execute_dft_c2r(c2rwz, wz, wz)
 
    wfi =  cmplx( real(vy)* real(wz) -  real(vz)* real(wy), &
                  aimag(vy)*aimag(wz) - aimag(vz)*aimag(wy) )

    wz = cmplx(   real (vz)*real (wx) - real (vx)*real (wz)               &
               ,  aimag(vz)*aimag(wx) - aimag(vx)*aimag(wz)                       &
               )                        
    wx = cmplx(   real (vx)*real (wy) - real (vy)*real (wx)               &
               ,  aimag(vx)*aimag(wy) - aimag(vy)*aimag(wx)                       &
               )
 
    wy = wz*const
    wz = wx*const
    wx = wfi*const
 
    call dfftw_execute_dft_r2c(r2cwx, wx, wx)
    call dfftw_execute_dft_r2c(r2cwy, wy, wy)
    call dfftw_execute_dft_r2c(r2cwz, wz, wz)
 
 
    !======================== projecting ================================= 
    wfi = wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz)
    wy = wx*hpx+wy*hpy+wz*hpz
    wx = wfi
    if (forced) then
      where (k2 .lt. fk22 .and. k2 .ge. fk12)
        wx = wx + fcoe * ap
        wy = wy + fcoe * am
      endwhere
    end if

    !=============== advection term for scalar ====================

    wfi = fi
    call dfftw_execute_dft_c2r(c2rwfi, wfi, wfi) 
    ! input will be destroyed by fftw even for out of place fftw.

    vx = vx * wfi * const
    vy = vy * wfi * const
    vz = vz * wfi * const
    call dfftw_execute_dft_r2c(r2cvx, vx, vx)
    call dfftw_execute_dft_r2c(r2cvy, vy, vy)
    call dfftw_execute_dft_r2c(r2cvz, vz, vz)

    wfi = - eye * ( kx * vx + ky * vy + kz * vz )

    !================== forcing by a mean gradient ================
    wfi = wfi + ffi
 
    !====================== time step size ===================================
    ! umax is updated only every ieout steps
    dt=min(min(beta/(2*om+mytiny), betadx/umax),3._sp*dto)
    if (dt .ge. twrite-time) then
            dt = twrite - time
    else if (twrite-time-dt .lt. .2_sp*dt) then
            dt = .5_sp * (twrite - time)
    end if
 
    !======================= time stepping ======================================
    select case (istep)
    case (0)
      dt = dt / 2

      open(88,file = "./out/oap.dat",status="unknown",form="unformatted")
        write(88) wx
      close(88)
      open(88,file = "./out/oam.dat",status="unknown",form="unformatted")
        write(88) wy
      close(88)
      open(88,file = "./out/ophi.dat",status="unknown",form="unformatted")
        write(88) wfi
      close(88)

      vz = exp(-k2*dt*rnu+2._sp*eye*dt*kz*om/(sqrt(k2)+mytiny))  ! positive helical wave
      ap = ( ap + dt * wx ) * vz
 
      vz = conjg(vz)                                             ! negative helical wave
      am = ( am + dt * wy ) * vz
      
      vz = exp(-rnu * k2 * dt / pr)
      fi = ( fi + dt * wfi ) * vz

      time = time + dt
      dto = dt
 
    case (1)
      dt = dt / 2
      dtt=dt/dto
 
      s= .5_sp*(1.+dtt*dtt)
      ek=.5_sp*(1.-dtt*dtt)
      delta=.5_sp*dto*(1.+dtt)**2
 
      vz = -k2*rnu + 2._sp*eye*kz*om/(sqrt(k2)+mytiny)  ! positive helical wave
      vx = exp(vz * dto)
      vy = exp(vz * dt)

      ap = s * oap * vx + ek * ap + delta * wx
      ap = ap * vy

      am = s * oam * conjg(vx) + ek * am + delta * wy
      am = am * conjg(vy)

      vz = -rnu * k2 / pr
      vx = exp(vz * dto )
      vy = exp(vz *  dt )
      fi = s * ofi * vx + ek * fi + delta * wfi
      fi = fi * vy
      ! ==========================================
 
      time = time + dt
      dto = dt + dto

      open(88,file = "./out/oap.dat",status="unknown",form="unformatted")
        read(88) oap
      close(88)
      open(88,file = "./out/oam.dat",status="unknown",form="unformatted")
        read(88) oam
      close(88)
      open(88,file = "./out/ophi.dat",status="unknown",form="unformatted")
        read(88) ofi
      close(88)
 
    case (2:)
      dtt=dt/dto
      s=(1._sp+.5_sp*dtt)*dt
      ek=-.5_sp*dtt*dt
 

      vz = -k2*rnu + 2._sp*eye*kz*om/(sqrt(k2)+mytiny)  ! positive helical wave
      vx = exp(vz * dto)
      vy = exp(vz * dt)

      ! =================== ap ===================
      ap = ap + s * wx + ek * oap * vx
      ap = ap * vy
 
      ! =================== am ===================
      am = am + s * wy + ek * oam * conjg(vx)
      am = am * conjg(vy)
      ! ==========================================

      vz = -k2*rnu/pr
      vx = exp(vz * dto)
      vy = exp(vz * dt)
      fi = fi + s * wfi + ek * ofi * vx
      fi = fi * vy
 
      oap = wx
      oam = wy
      ofi = wfi

 
      time = time + dt
      dto = dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
      stop
    end select
    ! dealiasing
    where(k2 .ge. kcut2)
        ap = 0._sp
        am = 0._sp
        fi = 0._sp
    endwhere
    ap(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    am(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    
    vx = ap*hpx + am*conjg(hpx)
    vy = ap*hpy + am*conjg(hpy)
    vz = ap*hpz + am*conjg(hpz)
  
  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  !  destroy plans
  call destroyplans

  !  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, oap, oam)
  deallocate(ap, am, hpx, hpy, hpz, fi, ofi, wfi, ffi)
  deallocate(kx, ky, kz, k2, ekk)

  write(*,*) 'spectral_rot_dns_vt_hm.x finished '
  stop

end program spectral

