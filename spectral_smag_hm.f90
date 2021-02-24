program spectral_dns

  use mconstant
  use mwavenumber
  use mfftwplan3d
  use msymmetrize
  use minitialize
  use minput
  use moutput
  implicit none

  real(sp), parameter  :: c1=0.026 !cut-off filter

  integer  :: nx, ny, nz, nles, iseed, ieout, idp
  real(sp) :: timemax, beta, rnu, dtw, time, eps, fkmax, akp, u0
  logical  :: new, forced
  ! parameters should be input from 'parameter_les.d' initially.

  complex(sp), allocatable, dimension (:,:,:) :: vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) :: ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) :: t11,t12,t13,t22,t23,t33
  real(sp),    allocatable, dimension (:,:,:) :: k2,k2e,kx,ky,kz
  real(sp),    allocatable, dimension (:)     :: eks

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: delta, S, dtt, ek, dt, dto, betadx, umax, twrite
  real(sp) :: kcut, kcut2, ignore_me, const
  real(sp) :: coeforce, kinet, fkmax2

  open(90,file='parameter_les.d',status='unknown')
    read(90,*) forced
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) nles
    read(90,*) timemax
    read(90,*) ieout
    read(90,*) beta
    read(90,*) rnu
    read(90,*) dtw
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) eps
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz, 'nles', nles
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'beta', beta
  write(*,*) 'rnu', rnu
  write(*,*) 'dtw', dtw
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  const = 1._sp/(nx*ny*nz)
  delta = 2._sp*pi/nles
  kcut = nles/2._sp
  kcut2 = kcut * kcut

  betadx = beta * (2*pi/nx) 
  fkmax2 = fkmax * fkmax

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(k2e(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (new) then
    ! Generate initial condition -------
    ! The initial fields have been dealiased in initialize
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
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

  ! specify umax empirically as three times the rms value, in reference to
  ! Gaussian distribution
  umax=3.*sqrt(3.)*u0
  
  dto = timemax
  twrite = time
  ifile = idp
  istep = -1

  !============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (time .le. timemax)

    istep = istep + 1
 
    ! Output the total energy in k space
    ! For monitoring purpose only. May want to reduce this part of calculation by
    ! increasing ieout.
    if (mod(istep,ieout).eq.0) then
      ! calculating skewness is relatively time consuming. May need it only when
      ! debugging.
      !call skewness(vx,kx,s,nx,ny,nz)
      k2e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2e(1,:,:)=0.5_sp*k2e(1,:,:)
      ek = sum(k2e)
      write( *,'(I6,20E12.4)')istep,time,dto,ek,s
      write(25,'(I6,20E12.4)')istep,time,dto,ek,s

      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.5*sqrt(2.*ek)
    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. myeps) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      k2e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2e(1,:,:)=0.5_sp*k2e(1,:,:)

      write(24,*)time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + k2e(ii,jj,kk)
            end if
      end do
      end do
      end do
      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
 
    end if

    ! =========== Calculate vorticity =============================
    wx = eye * ky * vz - eye * kz * vy
    wy = eye * kz * vx - eye * kx * vz
    wz = eye * kx * vy - eye * ky * vx
 
    !========== calculate convection term: lamb vector =======================
    t11=vx
    t22=vy
    t33=vz
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d, wx,ignore_me) 
    call rfftwnd_f77_one_complex_to_real(c2r3d, wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d, wz,ignore_me)
 
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
 
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)
 
    !==================== Smagorinsky model ==========================
    t11 = eye * kx * vx
    t22 = eye * ky * vy
    t33 = eye * kz * vz
    t12 = eye * ( kx * vy + ky * vx ) * .5
    t13 = eye * ( kx * vz + kz * vx ) * .5
    t23 = eye * ( ky * vz + kz * vy ) * .5

    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)

    k2e = real(t11)**2 + real(t22)**2 + real(t33)**2 + 2.*(real(t12)**2 + real(t13)**2 + real(t23)**2)   
    k2e = sqrt(2. * k2e)

    t11 = cmplx( k2e * real(t11), aimag(t11) )
    t22 = cmplx( k2e * real(t22), aimag(t22) )
    t33 = cmplx( k2e * real(t33), aimag(t33) )
    t12 = cmplx( k2e * real(t12), aimag(t12) )
    t13 = cmplx( k2e * real(t13), aimag(t13) )
    t23 = cmplx( k2e * real(t23), aimag(t23) )

    k2e = aimag(t11)**2 + aimag(t22)**2 + aimag(t33)**2 + 2.*(aimag(t12)**2 + aimag(t13)**2 + aimag(t23)**2)   
    k2e = sqrt(2. * k2e)

    t11 = -2. * c1 * delta * delta *  cmplx( real(t11), k2e * aimag(t11) ) * const
    t22 = -2. * c1 * delta * delta *  cmplx( real(t22), k2e * aimag(t22) ) * const
    t33 = -2. * c1 * delta * delta *  cmplx( real(t33), k2e * aimag(t33) ) * const
    t12 = -2. * c1 * delta * delta *  cmplx( real(t12), k2e * aimag(t12) ) * const
    t13 = -2. * c1 * delta * delta *  cmplx( real(t13), k2e * aimag(t13) ) * const
    t23 = -2. * c1 * delta * delta *  cmplx( real(t23), k2e * aimag(t23) ) * const

    call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)

    wx = wx - eye * ( kx*t11 + ky*t12 + kz*t13 )
    wy = wy - eye * ( kx*t12 + ky*t22 + kz*t23 )
    wz = wz - eye * ( kx*t13 + ky*t23 + kz*t33 )

    !========== combine convection, sgs stress and forcing ============
    if (forced) then
        where (k2 .lt. fkmax2)
          k2e = real(vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz))
        elsewhere
          k2e = 0.
        end where
        k2e(1,:,:) = .5_sp * k2e(1,:,:)
        kinet = sum( k2e )
  
        coeforce = eps / (2._sp * kinet)

        where ( k2 .lt. fkmax2 )
            wx = wx + coeforce * vx
            wy = wy + coeforce * vy
            wz = wz + coeforce * vz
        endwhere
    end if

    !========== projection ================================= 
    t11 = (kx*wx + ky*wy + kz*wz)/k2
    wx = wx - kx * t11
    wy = wy - ky * t11
    wz = wz - kz * t11

    !========== time stepping ===================================
    ! umax is updated only every ieout steps
    dt=min(betadx/umax,3._sp*dto)
    if (dt .ge. twrite-time) then
            dt=twrite-time
    else if (twrite-time-dt .lt. .2_sp*dt) then
            dt=.5_sp*(twrite-time)
    end if

    select case (istep)
    case (0)
      k2e=exp(-rnu*k2*dt)
      vx = ( vx + dt*wx ) * k2e 
      vy = ( vy + dt*wy ) * k2e 
      vz = ( vz + dt*wz ) * k2e 
 
      time = time + dt
      dto=dt
      call output(wx,wy,wz,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto

      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2
      k2e=exp(-rnu*k2*dto)

      ! ox,oy,oz are the initial velocity components now
      vx = S*ox*k2e + ek*vx + delta*wx 
      vy = S*oy*k2e + ek*vy + delta*wy 
      vz = S*oz*k2e + ek*vz + delta*wz 
 
      k2e=exp(-rnu*k2*dt)
      vx = vx*k2e
      vy = vy*k2e
      vz = vz*k2e

      time = time + dt
      dto=dt+dto
      call input(ox,oy,oz,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      k2e=exp(-rnu*k2*dto)

      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt

      vx = vx + S*wx + ek*ox*k2e
      vy = vy + S*wy + ek*oy*k2e
      vz = vz + S*wz + ek*oz*k2e
 
      k2e=exp(-rnu*k2*dt)
      ! k2e=k2e**dtt ! Is this faster than evaluating exp?

      vx = vx*k2e
      vy = vy*k2e
      vz = vz*k2e
 
      ox = wx
      oy = wy
      oz = wz
 
      time = time + dt
      dto=dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

    ! dealiasing
    where(k2 .ge. kcut2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
    vx(1,1,1) = 0._sp
    vy(1,1,1) = 0._sp
    vz(1,1,1) = 0._sp

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(kx, ky, kz, k2, k2e, eks)
  deallocate(t11,t12,t13,t22,t23,t33)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral_dns
