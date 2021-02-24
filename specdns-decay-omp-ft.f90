program spectral_dns

  use mconstant
  use mwavenumber
  use mfftwplan3d
  use minitialize
  use moutput
  use minput
  implicit none

  integer, parameter :: numth = 2

  integer  :: nx, ny, nz, iseed, ieout, idp, itout
  real(sp) :: timemax, dt, rnu, time, akp, u0
  logical  :: new
  ! parameters to read from input file

  complex(sp), allocatable, dimension (:,:,:) :: vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) :: ox,oy,oz
  real(sp),    allocatable, dimension (:,:,:) :: k2,expk2,tmp
  real(sp),    allocatable, dimension (:,:,:) :: kx,ky,kz 
  real(sp),    allocatable, dimension (:)     :: eks

  integer  :: lx, ly, lz, lx1 
  ! integer  :: OMP_get_thread_num
  integer  :: ii, jj, kk, ifile, istep, jstep, ll, lx2
  real(sp) :: ek, dt_h, kcut2, s

  integer, save :: TID

  !$omp threadprivate(TID)

  open(90,file='parameter_decaydnsomp_ft.d',status='unknown')
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
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
  
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
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
  kcut2 = (2.*lx/3)**2

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(expk2(lx1,ly,lz),tmp(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
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
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    where(k2 .ge. kcut2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
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
    if ( jstep .ge. 0 .and. mod(jstep,itout) .eq. 0 ) then

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
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========== projection ================================= 
    !$omp parallel num_threads(numth)
    !$omp workshare
    tmp = (kx*real(wx,sp) + ky*real(wy,sp) + kz*real(wz,sp))/k2
    wx = cmplx(real(wx,sp) - kx*tmp, aimag(wx))
    wy = cmplx(real(wy,sp) - ky*tmp, aimag(wy))
    wz = cmplx(real(wz,sp) - kz*tmp, aimag(wz))

    tmp = (kx*aimag(wx) + ky*aimag(wy) + kz*aimag(wz))/k2
    wx = cmplx(real(wx,sp), aimag(wx) - kx*tmp)
    wy = cmplx(real(wy,sp), aimag(wy) - ky*tmp)
    wz = cmplx(real(wz,sp), aimag(wz) - kz*tmp)
    !$omp end workshare
    !$omp end parallel
 
    !========== time stepping ===================================
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
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
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

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral_dns
