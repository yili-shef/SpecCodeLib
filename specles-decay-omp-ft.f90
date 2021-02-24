program spectralles

  use mconstant
  use mwavenumber
  use mfftwplan3d
  use minitialize
  use moutput
  use minput
  implicit none

  integer, parameter :: numth = 4
  real(sp), parameter  :: c1=0.026 !cut-off filter

  integer  :: nx, ny, nz, iseed, ieout, idp, itout
  real(sp) :: timemax, dt, rnu, time, akp, u0
  logical  :: new
  ! parameters to read from input file

  complex(sp), allocatable, dimension (:,:,:) :: vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) :: ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) :: t11,t12,t13,t22,t23,t33
  real(sp),    allocatable, dimension (:,:,:) :: k2,expk2,tmp
  real(sp),    allocatable, dimension (:,:,:) :: kx,ky,kz 
  real(sp),    allocatable, dimension (:)     :: eks

  integer  :: lx, ly, lz, lx1 
  ! integer  :: OMP_get_thread_num
  integer  :: ii, jj, kk, ifile, istep, jstep, ll, lx2, nles
  real(sp) :: ek, dt_h, kcut2, s, delta, const, ignore_me

  integer, save :: TID

  !$omp threadprivate(TID)

  open(90,file='parameter_decaysmagomp_ft.d',status='unknown')
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) nles
    read(90,*) timemax
    read(90,*) ieout
    read(90,*) itout
    read(90,*) dt
    read(90,*) rnu
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) akp
    read(90,*) u0
    read(90,*) iseed
  close(90)
  
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz, 'nles', nles
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'itout', itout
  write(*,*) 'dt', dt
  write(*,*) 'rnu', rnu
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1
  lx2 = lx * lx

  const = 1._sp/(nx*ny*nz)
  delta = 2*pi/nles
  kcut2 = (nles/2.)**2

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(expk2(lx1,ly,lz),tmp(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  expk2 = exp(-k2*rnu*dt)
  expk2(1,1,1) = 1._sp

  if (new) then
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)

    where(k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere

    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)

    where(k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere
    write(*,*) 'initial velocity field readed.'
  endif

  
  ox=vx
  oy=vy
  oz=vz

  dt_h=.5_sp*dt
  ifile = idp
  istep = 0
  jstep = 0

  !============== begin the loop ====================================
  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

 
    ! Output the total energy in k space
    ! For monitoring purpose only. May want to reduce this part of calculation by
    ! increasing ieout.
    if (mod(jstep,ieout).eq.0) then

      ! calculating skewness is relatively time consuming. May need it only when
      ! debugging.
      ! call skewness(vx,kx,s,nx,ny,nz)

      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      tmp(1,:,:)=0.5_sp*tmp(1,:,:)
      ek = sum(tmp)
      write( *,'(I6,20E12.4)')istep,time,ek,s
      write(25,'(I6,20E12.4)')istep,time,ek,s

    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    ! if (abs(time-twrite) .le. myeps) then
    if (jstep .ge. 0 .and. mod(jstep,itout).eq.0) then

      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      tmp(1,:,:)=0.5_sp*tmp(1,:,:)

      write(24,*)time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + tmp(ii,jj,kk)
            end if
      end do
      end do
      end do
      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
 
    end if

    !$omp parallel num_threads(numth)
    !$omp workshare
    wx = eye * ky * vz - eye * kz * vy
    wy = eye * kz * vx - eye * kx * vz
    wz = eye * kx * vy - eye * ky * vx
    !$omp end workshare
    !$omp end parallel
 
    !========== calculate convection term: lamb vector =======================

    !$omp parallel workshare num_threads(numth)
    t11=vx
    t22=vy
    t33=vz
    !$omp end parallel workshare
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
 
    !$omp parallel workshare num_threads(numth)
    t12=cmplx(  real (t22)*real (wz)-real (t33)*real (wy)               &
             ,  aimag(t22)*aimag(wz)-aimag(t33)*aimag(wy)                       &
             )
    t13=cmplx(   real(t33)*real (wx)-real (t11)*real (wz)               &
             ,  aimag(t33)*aimag(wx)-aimag(t11)*aimag(wz)                       &
             )
    t23=cmplx(   real(t11)*real (wy)-real (t22)*real (wx)               &
             ,  aimag(t11)*aimag(wy)-aimag(t22)*aimag(wx)                       &
             )
 
    wx=t12*const
    wy=t13*const
    wz=t23*const
    !$omp end parallel workshare 
 
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
 
    !==================== Smagorinsky model ==========================

    !$omp parallel num_threads(numth)
    !$omp workshare
    t11 = eye * kx * vx
    t22 = eye * ky * vy
    t33 = eye * kz * vz
    t12 = eye * ( kx * vy + ky * vx ) * .5
    t13 = eye * ( kx * vz + kz * vx ) * .5
    t23 = eye * ( ky * vz + kz * vy ) * .5
    !$omp end workshare
    !$omp end parallel

    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)

    !$omp parallel num_threads(numth)
    !$omp workshare
    tmp = real(t11)**2 + real(t22)**2 + real(t33)**2 + 2.*(real(t12)**2 + real(t13)**2 + real(t23)**2)   
    tmp = sqrt(2. * tmp)

    t11 = cmplx( tmp * real(t11), aimag(t11) )
    t22 = cmplx( tmp * real(t22), aimag(t22) )
    t33 = cmplx( tmp * real(t33), aimag(t33) )
    t12 = cmplx( tmp * real(t12), aimag(t12) )
    t13 = cmplx( tmp * real(t13), aimag(t13) )
    t23 = cmplx( tmp * real(t23), aimag(t23) )

    tmp = aimag(t11)**2 + aimag(t22)**2 + aimag(t33)**2 + 2.*(aimag(t12)**2 + aimag(t13)**2 + aimag(t23)**2)   
    tmp = sqrt(2. * tmp)

    t11 = -2. * c1 * delta * delta *  cmplx( real(t11), tmp * aimag(t11) ) * const
    t22 = -2. * c1 * delta * delta *  cmplx( real(t22), tmp * aimag(t22) ) * const
    t33 = -2. * c1 * delta * delta *  cmplx( real(t33), tmp * aimag(t33) ) * const
    t12 = -2. * c1 * delta * delta *  cmplx( real(t12), tmp * aimag(t12) ) * const
    t13 = -2. * c1 * delta * delta *  cmplx( real(t13), tmp * aimag(t13) ) * const
    t23 = -2. * c1 * delta * delta *  cmplx( real(t23), tmp * aimag(t23) ) * const

    !$omp end workshare
    !$omp end parallel

    call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)

    !$omp parallel num_threads(numth)
    !$omp workshare
    wx = wx - eye * ( kx*t11 + ky*t12 + kz*t13 )
    wy = wy - eye * ( kx*t12 + ky*t22 + kz*t23 )
    wz = wz - eye * ( kx*t13 + ky*t23 + kz*t33 )
    !$omp end workshare
    !$omp end parallel

    ! ================ projection ================================= 

    !$omp parallel num_threads(numth)
    !$omp workshare
    tmp = (kx*real(wx) + ky*real(wy) + kz*real(wz)) / k2
    wx = cmplx(real(wx) - kx*tmp, aimag(wx))
    wy = cmplx(real(wy) - ky*tmp, aimag(wy))
    wz = cmplx(real(wz) - kz*tmp, aimag(wz))
    tmp = (kx*aimag(wx) + ky*aimag(wy) + kz*aimag(wz)) / k2
    wx = cmplx(real(wx), aimag(wx) - kx*tmp)
    wy = cmplx(real(wy), aimag(wy) - ky*tmp)
    wz = cmplx(real(wz), aimag(wz) - kz*tmp)
    !$omp end workshare
    !$omp end parallel
 
    ! ============== time stepping ===================================
    select case (istep)
    case (0)
      vx = ( vx + dt_h * wx )  * sqrt(expk2)
      vy = ( vy + dt_h * wy )  * sqrt(expk2)
      vz = ( vz + dt_h * wz )  * sqrt(expk2)
      
      time = time + dt_h
      istep = istep + 1
      jstep = -1 ! only half of first step, set to -1 to avoid duplicate output

      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      ! ox,oy,oz are the initial velocity components now

      !$omp parallel num_threads(numth)
      !$omp workshare
      vx = ox * expk2 + dt * wx * sqrt(expk2)
      vy = oy * expk2 + dt * wy * sqrt(expk2)
      vz = oz * expk2 + dt * wz * sqrt(expk2)
      !$omp end workshare
      !$omp end parallel
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! the actual first step.

      call input(ox,oy,oz,0,lx1,ly,lz)
      ! after the input, ox,oy,oz store the components of the rhs of the equation
      ! at previous time step: wx,wy,wz.
 
    case (2:)
      !$omp parallel num_threads(numth)
      !$omp workshare
      vx = ( vx + ( 3._sp*wx - expk2*ox ) * dt_h ) * expk2
      vy = ( vy + ( 3._sp*wy - expk2*oy ) * dt_h ) * expk2
      vz = ( vz + ( 3._sp*wz - expk2*oz ) * dt_h ) * expk2
      
      ox = wx
      oy = wy
      oz = wz
      !$omp end workshare
      !$omp end parallel

      time = time + dt
      istep = istep + 1
      jstep = jstep + 1
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

    ! ===== Dealiasing ===========================
    vx(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    vy(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode
    vz(1,1,1) = (0.0_sp,0.0_sp) ! zero kx=ky=kz=0 mode

    !!$omp parallel num_threads(numth)
    !!$omp workshare
    where(k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere
    !!$omp end workshare
    !!$omp end parallel
    
  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(kx, ky, kz, k2, expk2, eks, tmp)
  deallocate(t11,t12,t13,t22,t23,t33)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectralles
