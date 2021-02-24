program spectral

  use mconstant
  use mfftwplan3d
  use mwavenumber
  use mheliwave
  use moutput
  use minput
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, dtw, rnu, time, eps, fk1, fk2, akp, u0, lmbdx, lmbdy, lmbdz, pr
  logical  :: newvel, velforced, newphi, phiforced
  ! parameters should be input from 'parameter_rot_dns_scalar.d' initially.

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
  real(sp) :: twrite, kinet, kcut, kcut2, fk12, fk22, const, ignore_me

  write(*,*) 'Starting spectral_rot_scalar_dns_fftw2.x'

  open(90,file='parameter_rot_dns_scalar.d',status='unknown')
    read(90,*) velforced
    read(90,*) phiforced
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) timemax
    read(90,*) ieout
    read(90,*) beta
    read(90,*) rnu
    read(90,*) dtw
    read(90,*) time
    read(90,*) newvel
    read(90,*) newphi
    read(90,*) idp
    read(90,*) eps
    read(90,*) fk1
    read(90,*) fk2
    read(90,*) akp
    read(90,*) u0
    read(90,*) om
    read(90,*) lmbdx
    read(90,*) lmbdy
    read(90,*) lmbdz
    read(90,*) pr
    read(90,*) iseed
  close(90)
  
  write(*,*) 'velforced', velforced
  write(*,*) 'phiforced', phiforced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'beta', beta
  write(*,*) 'rnu', rnu
  write(*,*) 'dtw', dtw
  write(*,*) 'time', time
  write(*,*) 'newvel', newvel
  write(*,*) 'newphi', newphi
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'om', om
  write(*,*) 'mean scalar gradient (lmbdx, lmbdy, lmbdz)', lmbdx, lmbdy, lmbdz
  write(*,*) 'Prandtl number ', pr
  write(*,*) 'iseed', iseed

  lx=nx/2
  ly=ny
  lz=nz
  lx1=lx+1

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'FFTW plans created'

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
  write(*,*) 'Arrays allocated'


  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'Wavenumber defined'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'Heliwave defined'

  if (newvel) then
    ! generate initial condition -------
    call initvel(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'Velocity initialized'
  else 
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'Initial velocity field read.'
  end if

  if (newphi) then
    call initphi(fi,kx,ky,kz,k2,iseed,lx1,ly,lz,kcut)
    write(*,*) 'Scalar field initialized.'
  else
    call input(fi,'phi',idp,lx1,ly,lz)
    write(*,*) 'Initial phi fields read.'
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

      wx = fi*conjg(fi)
      wx(1,:,:) = 0.5_sp*wx(1,:,:)
      s = sum(real(wx)) ! scalar variance / 2
      write( *,'(i6,6e15.6)')istep, time, dto, ek, s
      write(25,'(i6,6e15.6)')istep, time, dto, ek, s
 
      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.5_sp*sqrt(2._sp*ek)
    end if
 
    ! output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. myeps) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      call output(fi,'phi',ifile,lx1,ly,lz)
      ifile = ifile + 1
      twrite=twrite+dtw
    end if
  
    ! =========== coefficient for forcing ==============
    if (velforced) then 

      where (k2 .lt. fk22 .and. k2 .ge. fk12)
        wx = ap*conjg(ap) + am*conjg(am) 
      elsewhere
        wx = 0._sp
      end where
      wx(1,:,:) = .5_sp * wx(1,:,:)
      kinet=sum(real(wx,sp))
  
      fcoe = eps/(2._sp*kinet+mytiny)
    end if

    ! ============= forcing term for scalar ===============
    if (phiforced) ffi = - lmbdx * vx - lmbdy * vy - lmbdz * vz

    ! =========== calculate vorticity =============================

    wx = eye * ( ky * vz - kz * vy )
    wy = eye * ( kz * vx - kx * vz )
    wz = eye * ( kx * vy - ky * vx )

    !=================== calculate convection term: lamb vector =======================

    call rfftwnd_f77_one_complex_to_real(c2r3d,vx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,vy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,vz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    
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
 
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
 
    !======================== projecting ================================= 
    wfi = wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz)
    wy = wx*hpx+wy*hpy+wz*hpz
    wx = wfi
    if (velforced) then
      where (k2 .lt. fk22 .and. k2 .ge. fk12)
        wx = wx + fcoe * ap
        wy = wy + fcoe * am
      endwhere
    end if

    !=============== advection term for scalar ====================
    wfi = fi
    call rfftwnd_f77_one_complex_to_real(c2r3d,wfi,ignore_me)

    vx = cmplx( real(vx,sp) * real(wfi,sp), aimag(vx) * aimag(wfi) ) * const
    vy = cmplx( real(vy,sp) * real(wfi,sp), aimag(vy) * aimag(wfi) ) * const
    vz = cmplx( real(vz,sp) * real(wfi,sp), aimag(vz) * aimag(wfi) ) * const

    call rfftwnd_f77_one_real_to_complex(r2c3d,vx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,vy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,vz,ignore_me)

    wfi = - eye * ( kx * vx + ky * vy + kz * vz )

    !================== forcing by a mean gradient ================
    if (phiforced) wfi = wfi + ffi
 
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
      vx = exp(vz * dto)
      vy = exp(vz *  dt)

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

      ! =================== fi ===================
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
    am(1,1,1) = (0.0_sp,0.0_sp) 
    fi(1,1,1) = (0.0_sp,0.0_sp)
    
    vx = ap*hpx + am*conjg(hpx)
    vy = ap*hpy + am*conjg(hpy)
    vz = ap*hpz + am*conjg(hpz)
  
  end do
  !=================== end of loop ==============================

  close(25)

  !  destroy plans
  call destroyplan3d

  !  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, oap, oam)
  deallocate(ap, am, hpx, hpy, hpz, fi, ofi, wfi, ffi)
  deallocate(kx, ky, kz, k2, ekk)

  write(*,*) 'spectral_rot_scalar_dns_fftw2.x finished '
  stop

end program spectral

