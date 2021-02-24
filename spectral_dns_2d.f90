program spectral_dns_2d
  use mconstant
  use mwavenumber
  use minput
  use moutput
  use minitialize
  use mfftwplan2d
  use msymmetrize
  implicit none

  integer  :: nx, ny, iseed, ieout, itout, idp, nhyper, nhypo
  real(sp) :: timemax, sigma, rnu, dt, time, eps, fkmin, fkmax, akp, u0
  logical  :: new, forced
  ! parameters should be input from 'parameter_dns_2d.d' initially.

  complex(sp), allocatable, dimension (:,:) ::  vx,vy
  complex(sp), allocatable, dimension (:,:) ::  wz,oz
  real(sp),    allocatable, dimension (:,:) ::  kx,ky,k2,k2exp
  real(sp),    allocatable, dimension (:)   ::  eks

  integer  :: lx, ly, lx1
  integer  :: ii, jj, ifile, istep, jstep, ll
  real(sp) :: ignore_me, const, S, ek, dt_h, fkmin2, fkmax2
  real(sp) :: wzrms

  open(90,file='parameter_dns_2d.d',status='unknown')
    read(90,*) forced
    read(90,*) nx
    read(90,*) ny
    read(90,*) timemax
    read(90,*) ieout
    read(90,*) itout
    read(90,*) sigma
    read(90,*) rnu
    read(90,*) nhyper
    read(90,*) nhypo
    read(90,*) dt
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) eps
    read(90,*) fkmin
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'nx', nx, 'ny', ny
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'itout', itout
  write(*,*) 'sigma', sigma
  write(*,*) 'rnu', rnu
  write(*,*) 'nhyper', nhyper
  write(*,*) 'nhypo', nhypo
  write(*,*) 'dt', dt
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fkmin', fkmin
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  
  call fftwplan2d(nx,ny)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lx1=lx+1

  const = 1. / ( nx * ny )
  dt_h = dt / 2.
  fkmin2 = fkmin * fkmin
  fkmax2 = fkmax * fkmax

  allocate(kx(lx1,ly),ky(lx1,ly),k2(lx1,ly))
  allocate(k2exp(lx1,ly))
  allocate(vx(lx1,ly),vy(lx1,ly))
  allocate(wz(lx1,ly),oz(lx1,ly))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,k2,lx1,ly)
  write(*,*) 'after wavenumber'

  ! Integration factor with hypervisocity and hypoviscocity
  k2exp = exp( -((sqrt(k2)*sigma)**nhyper + (rnu/sqrt(k2))**nhypo) * dt )

  if (new) then
    ! Generate initial condition -------
    call initialize(wz,k2,iseed,lx1,ly,akp,wzrms,2.0*ly/3)
    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(wz,idp,lx1,ly)
    write(*,*) 'initial velocity field readed.'
  endif

  oz=wz

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
    if ( mod(istep,ieout) .eq. 0 ) then
      ! calculating skewness is relatively time consuming. May need it only when
      ! debugging.
      !call skewness(vx,kx,s,nx,ny,nz)
      vx = wz * conjg(wz) / k2
      vx(1,:)=0.5_sp*vx(1,:)
      ek = sum(real(vx))
      write( *,'(I6,20E12.4)')istep,time,ek,s
      write(25,'(I6,20E12.4)')istep,time,ek,s

    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if ( mod(jstep, itout) .eq. 0 ) then

      call output(wz,ifile,lx1,ly)
      ifile = ifile + 1
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      vx = wz * conjg(wz) / k2
      vx(1,:)=0.5_sp*vx(1,:)

      write(24,*)time
      eks=0.
      do jj = 1, ly
      do ii = 1, lx1
          ll=floor(sqrt(k2(ii,jj))+.5)
          if (ll .ge. 1 .and. ll .le. lx) eks(ll) = eks(ll) + real( vx(ii,jj) )
      end do
      end do
      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
 
    end if

    ! Velocity components
    vx =   eye * ky * wz / k2
    vy = - eye * kx * wz / k2

    
    ! =========== advection term =============================
    call rfftwnd_f77_one_complex_to_real(c2r2d,vx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,vy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,wz,ignore_me)

    vx = vx * wz
    vy = vy * wz
    call rfftwnd_f77_one_real_to_complex(r2c2d,vx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c2d,vy,ignore_me)
    vx = vx * const
    vy = vy * const

    vx = eye * kx * vx
    vy = vx + eye * ky * vy  ! vy is the advection term
 
    !========== combine advection term and forcing ============
    if (forced) then 
        where ( k2 .ge. fkmin2 .and. k2 .le. fkmax2 )
            vy = -vy + eps / (conjg(wz) + mytiny) ! forcing is added
        endwhere
    end if 
 
    ! Remove Nyquist and dealiasing
    call symmetrize(vy,k2,2.*lx/3,lx1,ly)

    !========== time stepping ===================================

    select case (istep)
    case (0)
        wz = ( wz + dt_h*vy ) * k2exp 
       
        time = time + dt_h
        istep = istep + 1
        jstep = 1 ! only half of first step, set to one to avoid duplicate output
       
        call output(vy,0,lx1,ly)
    case (1)
        wz = oz * k2exp + dt * vy * sqrt(k2exp)

        time = time + dt_h
        istep = istep + 1
        jstep = 1 ! the actual first step.

        call input(oz,0,lx1,ly)
    case (2:)
        wz = ( wz + dt_h * ( 3. * vy - oz * k2exp ) ) * k2exp
        oz = vy

        time = time + dt
        istep = istep + 1
        jstep = jstep + 1

    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, wz, oz)
  deallocate(kx, ky, k2, k2exp, eks)

  ! Destroy plans
  call destroyplan2d

  write(*,*)'finished '
  stop

end program spectral_dns_2d
