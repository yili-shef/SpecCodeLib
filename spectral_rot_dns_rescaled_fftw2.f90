program spectral

  use mconstant
  use mwavenumber
  use mheliwave
  use moutput
  use minput
  use mfftwplan3d
  implicit none

  integer  :: iseed, ieout, idp, nx, ny, nz
  real(sp) :: timemax, beta, dtw, rnu, time, fp1, fp2, akp, u0 
  logical  :: new
  ! parameters should be input from 'parameter_rot_dns_rescaled.d' initially.

  integer  :: lx, ly, lz, lx1
  
  integer, parameter :: numth = 4

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  oap,oam
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  real(sp),    allocatable, dimension (:,:,:) ::  k2,kx,ky,kz
  real(sp),    allocatable, dimension (:)     ::  ekk

  integer  :: istep, ifile
  real(sp) :: delta, om, s, dtt, dt, ek, dto, betadx, umax
  real(sp) :: twrite, kcut, kcut2, const, ignore_me

  open(90,file='parameter_rot_dns_rescaled.d',status='unknown')
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
    read(90,*) fp1
    read(90,*) fp2
    read(90,*) akp
    read(90,*) u0
    read(90,*) om
    read(90,*) iseed
  close(90)
  
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'beta', beta
  write(*,*) 'rnu', rnu
  write(*,*) 'dtw', dtw
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'fp1', fp1
  write(*,*) 'fp2', fp2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'om', om
  write(*,*) 'iseed', iseed

  lx=nx/2
  ly=ny
  lz=nz
  lx1=lx+1

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after creating plans'

  const = 1._sp/(nx*ny*nz)

  kcut = 2._sp*lx/3
  delta = pi / real(kcut,sp)
  betadx = beta*delta

  kcut2 = kcut * kcut

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(oap(lx1,ly,lz),oam(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))
  allocate(ekk(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'

  if (new) then
    ! generate initial condition -------
    call initvel(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'after initialize'
  else
    ! reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    print*, 'initial velocity field read.'
  endif
  
  ap = vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am = vx*hpx+vy*hpy+vz*hpz
  oap = ap
  oam = am

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
      write( *,'(i6,6e15.6)')istep, time, dto, ek, s
      write(25,'(i6,6e15.6)')istep, time, dto, ek, s
 
      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.5_sp*sqrt(2._sp*ek)
    end if
 
    ! output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. myeps) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      twrite=twrite+dtw
    end if
  
  
    ! =========== calculate vorticity =============================
    !$omp parallel workshare num_threads(numth)
    kx = sqrt(k2)
    wx = kx * (ap * hpx - am * conjg(hpx))
    wy = kx * (ap * hpy - am * conjg(hpy))
    wz = kx * (ap * hpz - am * conjg(hpz))
    !$omp end parallel workshare 

    !=================== calculate convection term: lamb vector =======================
    call rfftwnd_f77_one_complex_to_real(c2r3d,vx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,vy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,vz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    
 
    !$omp parallel workshare num_threads(numth)
    kx =  real(vy)* real(wz) -  real(vz)* real(wy)           
    ky = aimag(vy)*aimag(wz) - aimag(vz)*aimag(wy)

    wz = cmplx(   real (vz)*real (wx) - real (vx)*real (wz)               &
               ,  aimag(vz)*aimag(wx) - aimag(vx)*aimag(wz)                       &
               )                        
    wx = cmplx(   real (vx)*real (wy) - real (vy)*real (wx)               &
               ,  aimag(vx)*aimag(wy) - aimag(vy)*aimag(wx)                       &
               )
 
    wy = wz*const
    wz = wx*const
    wx = cmplx(kx,ky)*const
    !$omp end parallel workshare 
 
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
 
    !======================== projecting ================================= 
    vx = wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz)
    vy = wx*hpx+wy*hpy+wz*hpz
 
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
        write(88) vx
      close(88)
      open(88,file = "./out/oam.dat",status="unknown",form="unformatted")
        write(88) vy
      close(88)

      vz = exp(-k2*dt*rnu+2._sp*eye*dt*kz*om/(sqrt(k2)+mytiny))  ! positive helical wave
      ap = ( ap + dt * vx ) * vz
 
      vz = conjg(vz)                                             ! negative helical wave
      am = ( am + dt * vy ) * vz
      
      time = time + dt
      dto = dt
 
    case (1)
      dtt=dt/dto
 
      s= .5_sp*(1._sp+dtt*dtt)
      ek=.5_sp*(1._sp-dtt*dtt)
      delta=.5_sp*dto*(1._sp+dtt)**2
 
      vz = -k2*rnu + 2._sp*eye*kz*om/(sqrt(k2)+mytiny)  ! positive helical wave
      wx = exp(vz * dto)
      wy = exp(vz * dt)

      ap = s * oap * wx + ek * ap + delta * vx
      ap = ap * wy

      am = s * oam * conjg(wx) + ek * am + delta * vy
      am = am * conjg(wy)
      ! ==========================================
 
      time = time + dt
      dto = dt + dto

      open(88,file = "./out/oap.dat",status="unknown",form="unformatted")
        read(88) oap
      close(88)
      open(88,file = "./out/oam.dat",status="unknown",form="unformatted")
        read(88) oam
      close(88)
 
    case (2:)
      dtt=dt/dto
      s=(1._sp+.5_sp*dtt)*dt
      ek=-.5_sp*dtt*dt
 
      !$omp parallel workshare num_threads(numth)

      vz = -k2*rnu + 2._sp*eye*kz*om/(sqrt(k2)+mytiny)  ! positive helical wave
      wx = exp(vz * dto)
      wy = exp(vz * dt)

      ! =================== ap ===================
      ap = ap + s * vx + ek * oap * wx
      ap = ap * wy
 
      ! =================== am ===================
      am = am + s * vy + ek * oam * conjg(wx)
      am = am * conjg(wy)
      ! ==========================================
 
      oap = vx
      oam = vy

      !$omp end parallel workshare 
 
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
    endwhere
    ap(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    am(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode

    call force_rescale(ap,am,k2,lx1,ly,lz,fp1,fp2)
    
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
  deallocate(ap, am, hpx, hpy, hpz)
  deallocate(kx, ky, kz, k2, ekk)

  write(*,*) 'spectral_rot_dns_rescaled_fftw2.x finished '
  stop

end program spectral

