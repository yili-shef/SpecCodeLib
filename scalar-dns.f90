program scalardns

  use mconstant
  use mfftwplan3d
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, rnu, dtw, time, eps, epsphi, fkmax, akp, u0, pr, phirms
  logical  :: new, forced, scalarforced
  ! parameters should be input from 'parameter_scalar_dns.d' initially.

  complex(sp), allocatable, dimension (:,:,:) ::  va,vb,phi
  complex(sp), allocatable, dimension (:,:,:) ::  wa,wb,wphi
  complex(sp), allocatable, dimension (:,:,:) ::  oa,ob,ophi
  real(sp),    allocatable, dimension (:,:,:) ::  k2,k2_e
  real(sp),    allocatable, dimension (:)     ::  kx,ky,kz,eks

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: delta, S, dtt, ek, dt, dto, betadx, umax, twrite, const

  open(90,file='parameter_scalar_dns.d',status='unknown')
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
    read(90,*) fkmax
    read(90,*) akp
    read(90,*) u0
    ! parameter for scalar advection
    read(90,*) scalarforced
    read(90,*) pr
    read(90,*) epsphi
    read(90,*) phirms

    read(90,*) iseed
  close(90)
  
  write(*,*) 'forced', forced
  write(*,*) 'scalarforced', scalarforced
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
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'prandtl number', pr
  write(*,*) 'eps of phi', epsphi
  write(*,*) 'rms of phi', phirms
  write(*,*) 'iseed', iseed
  write(*,*) 'smallest', smallest
  write(*,*) 'oneless', oneless
  
  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'
  
  lx=nx/2
  ly=ny
  lz=nz

  lx1=lx+1

  delta = pi / real(lx,sp)
  betadx = beta * delta
  const = 1./(nx*ny*nz)

  allocate(kx(lx1),ky(ly),kz(lz), k2(lx1,ly,lz))
  allocate(k2_e(lx1,ly,lz),phi(lx1,ly,lz),wphi(lx1,ly,lz))
  allocate(va(lx1,ly,lz),vb(lx1,ly,lz))
  allocate(wa(lx1,ly,lz),wb(lx1,ly,lz))
  allocate(oa(lx1,ly,lz),ob(lx1,ly,lz))
  allocate(ophi(lx1,ly,lz))
  allocate(eks(lx))
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (new) then
    ! Generate initial condition -------
    call velinit(va,vb,k2,iseed,lx1,ly,lz,akp,u0)
    call phiinit(phi,k2,iseed,lx1,ly,lz,akp,phirms)

    write(*,*) 'after initialize'
  else
    ! Reading from disk ------
    call input(wa,wb,wphi,idp,lx1,ly,lz) !wa,wb,wphi are used as intermediate storage.
    call xyztocrayafull(wa,wb,wphi,va,vb,kx,ky,kz,lx1,ly,lz)
    call inputphi(phi,idp,lx1,ly,lz)
    write(*,*) 'initial velocity and scalar field readed.'
  endif

  
  oa=va
  ob=vb
  ophi=phi

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
  open(26,file='./post/scalar_spectrum.data',status='unknown')
  
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

      k2_e = va*conjg(va) + vb*conjg(vb)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)
      ek = sum(k2_e)

      k2_e = phi * conjg(phi)
      k2_e(1,:,:) = .5_sp * k2_e(1,:,:)
      phirms = sum(k2_e)

      write( *,'(I6, 10E12.5)')istep,time,dto,ek,s,phirms
      write(25,'(I6, 10E12.5)')istep,time,dto,ek,s,phirms

      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.*sqrt(2.*ek)
    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. smallest) then

      call crayatoxyzfull(va,vb,kx,ky,kz,wa,wb,wphi,lx1,ly,lz)
      call output(wa,wb,wphi,ifile,lx1,ly,lz)
      call outputphi(phi,ifile,lx1,ly,lz)

      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      k2_e = va*conjg(va) + vb*conjg(vb)
      k2_e(1,:,:)=0.5_sp*k2_e(1,:,:)

      write(24,*) '#', time
      eks=0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + k2_e(ii,jj,kk)
            end if
      end do
      end do
      end do

      do ii = 1, lx
        write(24,*) ii, eks(ii)
      end do
      write(24,*)

      ! Scalar spectrum
      k2_e = phi * conjg(phi)
      k2_e(1,:,:) = 0.5_sp * k2_e(1,:,:)

      write(26,*) '#', time
      eks = 0.
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    eks(ll) = eks(ll) + k2_e(ii,jj,kk)
            end if
      end do
      end do
      end do

      do ii = 1, lx
        write(26,*) ii, eks(ii)
      end do
      write(26,*)

 
    end if

    !========== calculate convection term: lamb vector =======================
    call convecdnsvel2cphi(va,vb,phi,kx,ky,kz,wa,wb,wphi,lx1,ly,lz,const)
    call setzero(wphi, k2, lx1, ly, lz)
    call setzero(wa, k2, lx1, ly, lz)
    call setzero(wb, k2, lx1, ly, lz)

 
    !========== combine convection term, and forcing ============
    if (forced) call velforce2c(va,vb,wa,wb,k2,lx1,ly,lz,fkmax,eps)
    if (scalarforced) call phiforce(phi,wphi,k2,lx1,ly,lz,fkmax,epsphi)
 
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
      k2_e=exp(-rnu*k2*dt)
      va = va * k2_e + dt*wa 
      vb = vb * k2_e + dt*wb 

      k2_e = exp(-rnu*k2*dt/pr)
      phi = phi * k2_e + dt * wphi
 
      time = time + dt
      dto=dt
      call outputvel2c(wa,wb,0,lx1,ly,lz)
      call outputphi(wphi,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto

      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2

      k2_e=exp(-rnu*k2*dto)

      ! ox,oy,oz are the initial velocity components now
      va = S*oa*k2_e + ek*va + delta*wa 
      vb = S*ob*k2_e + ek*vb + delta*wb 

      k2_e = exp(-rnu*k2*dto/pr)
      phi = S * ophi * k2_e + ek * phi + delta * wphi
 
      k2_e=exp(-rnu*k2*dt)
      va = va*k2_e
      vb = vb*k2_e

      k2_e = exp(-rnu*k2*dt/pr)
      phi = phi * k2_e

      time = time + dt
      dto=dt+dto
      call inputvel2c(oa,ob,0,lx1,ly,lz)
      call inputphi(ophi,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt

      k2_e=exp(-rnu*k2*dto)
      va = va + S*wa + ek*oa*k2_e
      vb = vb + S*wb + ek*ob*k2_e

      k2_e = exp(-rnu*k2*dto/pr)
      phi = phi + S * wphi + ek * ophi * k2_e
 
      k2_e=exp(-rnu*k2*dt)
      va = va*k2_e
      vb = vb*k2_e

      k2_e = exp(-rnu*k2*dt/pr)
      phi = phi * k2_e
 
      oa = wa
      ob = wb
      ophi = wphi
 
      time = time + dt
      dto=dt
 
    case default
      write(*,*) "istep wrong. istep = ", istep
    end select

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)
  close(26)

  ! Deallocate arrays
  deallocate(va, vb, wa, wb, oa, ob)
  deallocate(phi, wphi, ophi)
  deallocate(kx, ky, kz, k2, k2_e, eks)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

end program scalardns
