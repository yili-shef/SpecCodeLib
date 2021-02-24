program dnsdeform

  use mconstant
  use mwavenumber
  use mfftwplan3d
  use msymmetrize
  use minitialize
  use minput
  use moutput
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, itout, idp
  real(sp) :: timemax, rnu, dt, time, eps, fkmax, akp, u0
  logical  :: new, forced
  ! parameters should be input from 'parameter_dns_deform.d' initially.

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  dux,duy,duz
  complex(sp), allocatable, dimension (:,:,:) ::  kux,kuy,kuz
  real(sp),    allocatable, dimension (:,:,:) ::  k2,k2e
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz, eks

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: kcut, kcut2, S, ek

  open(90,file='parameter_dns_deform.d',status='unknown')
    read(90,*) forced
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
    read(90,*) eps
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'timemax', timemax
  write(*,*) 'ieout', ieout
  write(*,*) 'itout', itout
  write(*,*) 'dt', dt
  write(*,*) 'rnu', rnu
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fkmax', fkmax
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

  kcut = 2.*lx/3
  kcut2 = kcut * kcut

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(k2e(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kux(lx1,ly,lz),kuy(lx1,ly,lz),kuz(lx1,ly,lz))
  allocate(dux(lx1,ly,lz),duy(lx1,ly,lz),duz(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  k2e = exp(-rnu*k2*dt/2)

  if (new) then
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)

    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  endif

  
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
      dux = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      dux(1,:,:)=0.5_sp*dux(1,:,:)
      ek = sum(real(dux))
      write( *,'(I6,20E12.4)')istep,time,ek,s
      write(25,'(I6,20E12.4)')istep,time,ek,s

    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (mod(istep,itout) .eq. 0) then
      call output(vx,vy,vz,ifile,lx1,ly,lz)
      ifile = ifile + 1
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      dux = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      dux(1,:,:)=0.5_sp*dux(1,:,:)

      write(24,*)time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
          ll=floor(sqrt(k2(ii,jj,kk))+.5)
          if (ll .ge. 1 .and. ll .le. lx) then
              eks(ll) = eks(ll) + real(dux(ii,jj,kk))
          end if
      end do
      end do
      end do
      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
 
    end if

    ! == Classical RK4 method ==

    ! -- First substep: k1u --

    call ratesdudt(vx,vy,vz,wx,wy,wz)
    kux = dt * wx 
    kuy = dt * wy
    kuz = dt * wz

    ! increases in velocity
    dux = kux
    duy = kuy
    duz = kuz

    ! -- Second substep: k2u --

    ! un + k1u/2 to calculate k2u
    kux = k2e * (vx + kux/2)
    kuy = k2e * (vy + kuy/2)
    kuz = k2e * (vz + kuz/2)

    call ratesdudt(kux,kuy,kuz,wx,wy,wz)
    kux = dt * wx / k2e
    kuy = dt * wy / k2e
    kuz = dt * wz / k2e

    ! increases in velocity
    dux = dux + 2 * kux
    duy = duy + 2 * kuy
    duz = duz + 2 * kuz

    ! -- Third substep: k3u -- 

    ! un + k2u/2 to calculate k3u
    kux = k2e * (vx + kux/2) 
    kuy = k2e * (vy + kuy/2)
    kuz = k2e * (vz + kuz/2)

    call ratesdudt(kux,kuy,kuz,wx,wy,wz)
    kux = dt * wx / k2e
    kuy = dt * wy / k2e
    kuz = dt * wz / k2e

    ! increases in velocity
    dux = dux + 2 * kux
    duy = duy + 2 * kuy
    duz = duz + 2 * kuz

    ! -- Fourth substep: k4u --

    ! un + k3u to calculate k4u
    kux = k2e * k2e * (vx + kux) 
    kuy = k2e * k2e * (vy + kuy)
    kuz = k2e * k2e * (vz + kuz)

    call ratesdudt(kux,kuy,kuz,wx,wy,wz)
    kux = dt * wx / k2e / k2e
    kuy = dt * wy / k2e / k2e
    kuz = dt * wz / k2e / k2e

    ! increases in velocity
    dux = dux + kux
    duy = duy + kuy
    duz = duz + kuz

    ! -- Velocity at next time step --
    vx = k2e * k2e * ( vx + dux / 6 )
    vy = k2e * k2e * ( vy + duy / 6 )
    vz = k2e * k2e * ( vz + duz / 6 )

    ! == End of RK4 stepping ==

    ! Dealiasing
    where (k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere

    time = time + dt

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, kux, kuy, kuz, dux, duy, duz)
  deallocate(kx, ky, kz, k2, k2e, eks)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

contains

  subroutine ratesdudt(lvx,lvy,lvz,lwx,lwy,lwz)
    
    complex(sp), dimension(:,:,:) :: lvx,lvy,lvz,lwx,lwy,lwz
    integer :: ii, jj, kk

    ! =========== Calculate vorticity =============================
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      lwx(ii,jj,kk) = eye * ky(jj) * lvz(ii,jj,kk) - eye * kz(kk) * lvy(ii,jj,kk)
      lwy(ii,jj,kk) = eye * kz(kk) * lvx(ii,jj,kk) - eye * kx(ii) * lvz(ii,jj,kk)
      lwz(ii,jj,kk) = eye * kx(ii) * lvy(ii,jj,kk) - eye * ky(jj) * lvx(ii,jj,kk)
    end do
    end do
    end do
 
    !========== calculate conlvection term: lamb lvector =======================
    call convec_dns(lvx,lvy,lvz,lwx,lwy,lwz,lx1,ly,lz,nx,ny,nz)
 
    !========== combine conlvection term, and forcing ============
    if (forced) call force(lwx,lwy,lwz,lvx,lvy,lvz,k2,lx1,ly,lz,fkmax,eps)
 
    call symmetrize(lwx,lx1,ly,lz)
    call symmetrize(lwy,lx1,ly,lz)
    call symmetrize(lwz,lx1,ly,lz)

    !========== projection ================================= 
    call projection(lwx,lwy,lwz,kx,ky,kz,lx1,ly,lz)
 
  end subroutine ratesdudt

end program dnsdeform
