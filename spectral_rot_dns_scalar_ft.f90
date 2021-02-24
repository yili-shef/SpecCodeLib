program spectral

  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, nstep, itout, ieout, idp
  real(sp) :: dt, rnu, time, eps, fk1, fk2, akp, u0, om, pr, phirms, epsphi 
  logical  :: new, forced, scalarforced
  ! parameters should be input from 'parameter_rot_dns_scalar_ft.d' initially.

  integer  :: lx, ly, lz, lx1
  
  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,ap,am,efp,efm
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz
  complex(sp), allocatable, dimension (:,:,:) ::  hpx,hpy,hpz
  complex(sp), allocatable, dimension (:,:,:) ::  phi, ophi, wphi
  real(sp),    allocatable, dimension (:,:,:) ::  k2,tmp,efphi,kx,ky,kz

  integer  :: ii, jj, kk, istep, jstep, ifile
  real(sp) :: s, dt_h, ek 

  open(90,file='parameter_rot_dns_scalar_ft.d',status='unknown')
    read(90,*) forced
    read(90,*) nx
    read(90,*) ny
    read(90,*) nz
    read(90,*) nstep
    read(90,*) itout
    read(90,*) ieout
    read(90,*) dt
    read(90,*) rnu 
    read(90,*) time
    read(90,*) new
    read(90,*) idp
    read(90,*) eps
    read(90,*) fk1
    read(90,*) fk2
    read(90,*) akp
    read(90,*) u0
    read(90,*) om
    ! parameter for scalar advection
    read(90,*) scalarforced
    read(90,*) pr
    read(90,*) epsphi
    read(90,*) phirms

    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'nx', nx, 'ny', ny, 'nz', nz
  write(*,*) 'nstep', nstep
  write(*,*) 'itout', itout
  write(*,*) 'ieout', ieout
  write(*,*) 'dt', dt
  write(*,*) 'rnu', rnu
  write(*,*) 'time', time
  write(*,*) 'new', new
  write(*,*) 'idp', idp
  write(*,*) 'eps', eps
  write(*,*) 'fk1', fk1
  write(*,*) 'fk2', fk2
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'om', om
  write(*,*) 'scalarforced', scalarforced
  write(*,*) 'pr', pr
  write(*,*) 'epsphi', epsphi
  write(*,*) 'phirms', phirms
  write(*,*) 'iseed', iseed

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz
  
  lx1=lx+1

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(efp(lx1,ly,lz),efm(lx1,ly,lz),tmp(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(ox(lx1,ly,lz),oy(lx1,ly,lz),oz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))

  allocate(phi(lx1,ly,lz), wphi(lx1,ly,lz), ophi(lx1,ly,lz))
  allocate(efphi(lx1,ly,lz))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1
    efp(ii,jj,kk) = exp(-k2(ii,jj,kk)*rnu*dt &
                   +2._sp*eye*dt*kz(kk)*om/(sqrt(k2(ii,jj,kk))+mytiny))  ! positive helical wave
    efm(ii,jj,kk) = exp(-k2(ii,jj,kk)*rnu*dt & 
                   -2._sp*eye*dt*kz(kk)*om/(sqrt(k2(ii,jj,kk))+mytiny))  ! negative helical wave
    efphi(ii,jj,kk) = exp(-k2(ii,jj,kk)*rnu*dt/pr)
  enddo
  enddo 
  enddo
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after heliwave'


  if (new) then
    ! generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    call initializephi(phi,k2,iseed,lx1,ly,lz,akp,phirms)

    call dealiasing(vx,k2,lx1,ly,lz)
    call dealiasing(vy,k2,lx1,ly,lz)
    call dealiasing(vz,k2,lx1,ly,lz)

    call dealiasing(phi,k2,lx1,ly,lz)

    write(*,*) 'after initialize'
  else
    ! reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    call inputphi(phi,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field readed.'
  endif
  
  ox=vx
  oy=vy
  oz=vz
  ophi = phi
  ap=vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)
  am=vx*hpx+vy*hpy+vz*hpz

  dt_h=.5_sp*dt
  ifile = idp
  istep=0
  jstep=0 

!============== begin the loop ====================================

  open(24,file='./post/spectrum.data',status='unknown')
  open(25,file='./post/ener_time.data',status='unknown')
  open(26,file='./post/scalar_spectrum.data',status='unknown')
  
  write(*,*) 'starting the loop'
  do while (jstep .le. nstep)

    ! write out the total energy in k space

    if (mod(jstep,ieout).eq.0) then

      !call skewness(vx,kx,s,nx,ny,nz)

      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      tmp(1,:,:)=0.5_sp*tmp(1,:,:)
      ek = sum(tmp)

      tmp = phi * conjg(phi)
      tmp(1,:,:) = .5_sp * tmp(1,:,:)
      phirms = sum(tmp)

      write( *,'(i6,6f12.4)') jstep,time,ek,s, phirms
      write(25,'(i6,6f12.4)') jstep,time,ek,s, phirms
    end if
 
    ! output the field data and spectra, including the initial data when jstep=0.
    if (mod(jstep,itout).eq.0) then

      call output(vx,vy,vz,ifile,lx1,ly,lz)
      call outputphi(phi,ifile,lx1,ly,lz)

      ifile = ifile + 1
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      tmp(1,:,:)=0.5_sp*tmp(1,:,:)

      write(24,*) '#', time
      do ii=1,lx
        ek=sum(tmp,mask=(abs(sqrt(k2)-ii).lt.0.5_sp))
        write(24,*)ii,ek
      end do
 
      ! Scalar spectrum
      tmp = phi * conjg(phi)
      tmp(1,:,:)=0.5_sp*tmp(1,:,:)

      write(26,*) '#', time
      do ii=1,lx
        ek=sum(tmp,mask=(abs(sqrt(k2)-ii).lt.0.5_sp))
        write(26,*)ii,ek
      end do

    end if

    ! =========== calculate vorticity =============================
    wx = eye*(ky*vz - kz*vy)
    wy = eye*(kz*vx - kx*vz)
    wz = eye*(kx*vy - ky*vx)
  
    !=================== calculate convection term: lamb vector =======================
    call convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
    call convecphi_dns(vx,vy,vz,phi,wphi,kx,ky,kz,lx1,ly,lz,nx,ny,nz)
 
    !========================== forcing term =================================  
    if (forced) call force_rot(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fk1,fk2,eps)
    if (scalarforced) call forcephi_rot(phi,wphi,k2,lx1,ly,lz,fk1,fk2,epsphi)
 
    call dealiasing(wx,k2,lx1,ly,lz)
    call dealiasing(wy,k2,lx1,ly,lz)
    call dealiasing(wz,k2,lx1,ly,lz)
    call dealiasing(wphi,k2,lx1,ly,lz)
 
    !======================== projecting ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)

    !======================= time stepping ======================================
    select case (istep)
    case (0)
      ap = ( ap + dt_h*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz)) ) * sqrt(efp)
      am = ( am + dt_h*(wx*hpx+wy*hpy+wz*hpz)                      ) * sqrt(efm)
      
      phi = ( phi + dt_h * wphi ) * sqrt(efphi)

      time = time + dt_h
      istep=istep+1
      jstep=1 ! only half of first step, set to one to avoid duplicate output

      call output(wx,wy,wz,0,lx1,ly,lz)
      call outputphi(wphi,0,lx1,ly,lz)
 
    case (1)

      ! The first two steps constitute the Euler 2 method. Note the handling of the integration
      ! factor.

      ! ox,oy,oz are the initial velocity components now
      ap = (ox*conjg(hpx)+oy*conjg(hpy)+oz*conjg(hpz))*efp + & 
            dt*(wx*conjg(hpx)+wy*conjg(hpy)+wz*conjg(hpz))*sqrt(efp)
      am = (ox*hpx+oy*hpy+oz*hpz)*efm + dt*(wx*hpx+wy*hpy+wz*hpz)*sqrt(efm)

      phi = ophi * efphi + dt * wphi * sqrt(efphi)
 
      time = time + dt_h
      istep = istep + 1
      jstep = 1 ! the actual first step.

      call input(ox,oy,oz,0,lx1,ly,lz)
      call inputphi(ophi,0,lx1,ly,lz)
      ! after the input, ox,oy,oz store the components of the rhs of the equation
      ! at previous time step: wx,wy,wz.
 
    case (2:)

      ! AB2 method

      vx = 3._sp*wx - efp*ox
      vy = 3._sp*wy - efp*oy
      vz = 3._sp*wz - efp*oz
      ap = (ap + dt_h*(vx*conjg(hpx)+vy*conjg(hpy)+vz*conjg(hpz)))*efp
      
      vx = 3._sp*wx - efm*ox
      vy = 3._sp*wy - efm*oy
      vz = 3._sp*wz - efm*oz
      am = (am + dt_h*(vx*hpx+vy*hpy+vz*hpz))*efm
 
      ox = wx
      oy = wy
      oz = wz

      phi = ( phi + dt_h * ( 3._sp * wphi - ophi * efphi ) ) * efphi
      ophi = wphi
 
      time = time + dt
      istep = istep + 1
      jstep = jstep + 1
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

    vx=ap*hpx+am*conjg(hpx)
    vy=ap*hpy+am*conjg(hpy)
    vz=ap*hpz+am*conjg(hpz)

  end do
!=================== end of loop ==============================

  close(24)
  close(25)

!  deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, ox, oy, oz)
  deallocate(ap, am, hpx, hpy, hpz)
  deallocate(kx, ky, kz, k2, efp, efm)
  deallocate(phi, ophi, wphi, efphi)
  deallocate(tmp)

!  destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program spectral
