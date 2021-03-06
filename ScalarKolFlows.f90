#include "fftw3param.f90"
#include "constant_dp.f90"
#include "fftw3plan.f90"
#include "symmetrize.f90"
#include "wavenumber.f90"
#include "initialize.f90"
#include "input.f90"
#include "output.f90"

program ScalarKolFlows

  use mconstant
  use mfftw3plan
  use mwavenumber
  use minitialize
  use minput
  use moutput
  implicit none

  integer  :: nx, ny, nz, iseed, ieout, idp
  real(sp) :: timemax, beta, rnu, dtw, time, eps, fkmax, akp, u0, mgrd, fa
  logical  :: new, forced
  ! parameters should be input from 'parameter_scalarkol.d' initially.

  complex(sp), allocatable, dimension (:,:,:) ::  vx,vy,vz,fi
  complex(sp), allocatable, dimension (:,:,:) ::  vxt,vyt,vzt
  complex(sp), allocatable, dimension (:,:,:) ::  wx,wy,wz,wf
  complex(sp), allocatable, dimension (:,:,:) ::  ox,oy,oz,of
  real(sp),    allocatable, dimension (:,:,:) ::  k2,k2e,kx,ky,kz,tmp
  real(sp),    allocatable, dimension (:)     ::  eks, ekf

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: delta, S, dtt, ek, dt, dto, betadx, umax, twrite
  real(sp) :: kcut, kcut2, const

  integer(dp) :: IFTvx, IFTvy, IFTvz, FTvx, FTvy, FTvz, IFTwx, IFTwy, IFTwz, FTwx, FTwy, FTwz, IFTwf

  open(90,file='parameter_scalarkol.d',status='unknown')
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
    read(90,*) iseed
    read(90,*) fa
    read(90,*) mgrd
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
  write(*,*) 'fkmax', fkmax
  write(*,*) 'akp', akp
  write(*,*) 'u0', u0
  write(*,*) 'iseed', iseed
  write(*,*) 'fa', fa
  write(*,*) 'mgrd', mgrd
  
  lx=nx/2
  ly=ny
  lz=nz

  const = 1._sp/(nx*ny*nz)

  lx1=lx+1

  kcut = 2._sp*lx/3 
  !kcut = lx
  kcut2 = kcut * kcut
  delta = pi / real(kcut,sp)
  betadx = beta * delta

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  allocate(k2e(lx1,ly,lz),tmp(lx1,ly,lz))
  allocate( vx(lx1,ly,lz), vy(lx1,ly,lz), vz(lx1,ly,lz), fi(lx1,ly,lz))
  allocate(vxt(lx1,ly,lz),vyt(lx1,ly,lz),vzt(lx1,ly,lz))
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz), wf(lx1,ly,lz))
  allocate( ox(lx1,ly,lz), oy(lx1,ly,lz), oz(lx1,ly,lz), of(lx1,ly,lz))
  allocate(eks(lx),ekf(lx))
  write(*,*) 'after allocate'

  call dfftwplan3dc2r(vxt,nx,IFTvx)
  call dfftwplan3dc2r(vyt,nx,IFTvy)
  call dfftwplan3dc2r(vzt,nx,IFTvz)
 
  call dfftwplan3dr2c(vxt,nx,FTvx)
  call dfftwplan3dr2c(vyt,nx,FTvy)
  call dfftwplan3dr2c(vzt,nx,FTvz)
 
  call dfftwplan3dc2r(wx,nx,IFTwx)
  call dfftwplan3dc2r(wy,nx,IFTwy)
  call dfftwplan3dc2r(wz,nx,IFTwz)
 
  call dfftwplan3dr2c(wx,nx,FTwx)
  call dfftwplan3dr2c(wy,nx,FTwy)
  call dfftwplan3dr2c(wz,nx,FTwz)
 
  call dfftwplan3dc2r(wf,nx,IFTwf)
  write(*,*) 'after creating plans'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (new) then
    call initialize(vx,vy,vz,fi,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    write(*,*) 'after initialize'
  else
    call input(vx,vy,vz,fi,idp,lx1,ly,lz)
    write(*,*) 'initial fields readed.'
  endif

  ox=vx; oy=vy; oz=vz; of=fi

  ! specify umax empirically as three times the rms value, in reference to
  ! Gaussian distribution
  umax=3._sp*sqrt(3._sp)*u0
  
  dto = timemax; twrite = time
  ifile = idp; istep = -1

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
      k2e = real( vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz) )
      k2e(1,:,:)=0.5_sp*k2e(1,:,:)
      ek = sum(k2e)
      k2e = real( fi*conjg(fi))
      s = sum(k2)
      write( *,'(I6,20E12.4)')istep,time,dto,ek,s
      write(25,'(I6,20E12.4)')istep,time,dto,ek,s

      ! update umax
      ! umax is approximated as three times the rms velocity.
      umax=3.5_sp*sqrt(2._sp*ek)
    end if
 
    ! Output the field data and spectra, including the initial data when time=0.
    if (abs(time-twrite) .le. myeps) then
      call output(vx,vy,vz,fi,ifile,lx1,ly,lz)

      ifile = ifile + 1
      twrite=twrite+dtw
      
      ! tmp = q_ij(k), q_ij(k) is the energy density spectral tensor, except
      ! at kx=0, where tmp=.5*q_ij(k)
      k2e = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      k2e(1,:,:)=0.5_sp*k2e(1,:,:)

      eks=0._sp
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

      ! scalar spectrum
      k2e = real(fi * conjg(fi))
      k2e(1,:,:)=0.5_sp*k2e(1,:,:)

      ekf=0._sp
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ll=floor(sqrt(k2(ii,jj,kk))+.5)
        if (ll .ge. 1 .and. ll .le. lx) then
          ekf(ll) = ekf(ll) + k2e(ii,jj,kk)
        end if
      end do
      end do
      end do

      write(24,*)time
      do ii = 1, lx
        write(24,*) ii, eks(ii), ekf(ii)
      end do
 
    end if

    !======= Advection term: lamb vector ====================
    wx = eye * ky * vz - eye * kz * vy
    wy = eye * kz * vx - eye * kx * vz
    wz = eye * kx * vy - eye * ky * vx
 
    vxt = vx; vyt = vy; vzt = vz
    call dfftw_execute(IFTvx)
    call dfftw_execute(IFTvy)
    call dfftw_execute(IFTvz)

    call dfftw_execute(IFTwx)
    call dfftw_execute(IFTwy)
    call dfftw_execute(IFTwz)

    tmp = real (vyt)*real (wz)-real (vzt)*real (wy) 
    k2e = aimag(vyt)*aimag(wz)-aimag(vzt)*aimag(wy)    

    wz=cmplx(   real (vzt)*real (wx)-real (vxt)*real (wz)               &
             ,  aimag(vzt)*aimag(wx)-aimag(vxt)*aimag(wz)               &
             )
    wx=cmplx(   real (vxt)*real (wy)-real (vyt)*real (wx)               &
             ,  aimag(vxt)*aimag(wy)-aimag(vyt)*aimag(wx)               &
             )
 
    wy=wz*const
    wz=wx*const
    wx=cmplx(tmp,k2e)*const
 
    call dfftw_execute(FTwx)
    call dfftw_execute(FTwy)
    call dfftw_execute(FTwz)

    ! ====== Advection term of scalar equation ==============
    wf = fi
    call dfftw_execute(IFTwf)

    vxt = cmplx(real(vxt) * real(wf), aimag(vxt) * aimag(wf)) 
    vyt = cmplx(real(vyt) * real(wf), aimag(vyt) * aimag(wf)) 
    vzt = cmplx(real(vzt) * real(wf), aimag(vzt) * aimag(wf)) 

    call dfftw_execute(FTvx)
    call dfftw_execute(FTvy)
    call dfftw_execute(FTvz)

    wf = - eye * (kx * vxt + ky * vyt + kz * vzt)

    !======= combine advection term with forcing ============
    ! Fy = fa * cos(kf * x); kf = 1
    wy(2,1,1) = wy(2,1,1) + fa/2
    wf = wf - mgrd * vy ! mean gradient is in the mean velocity direction.


    !========== projection ==========

    vxt = (kx*wx + ky*wy + kz*wz)/k2
    wx = wx - kx * vxt
    wy = wy - ky * vxt
    wz = wz - kz * vxt
 
    !========== time stepping ===============
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

      fi = ( fi + dt * wf) * k2e
 
      time = time + dt
      dto=dt
      call output(wx,wy,wz,wf,0,lx1,ly,lz)
 
    case (1)
      dtt=dt/dto

      S= .5*(1.+dtt*dtt)
      ek=.5*(1.-dtt*dtt)
      delta=.5*dto*(1.+dtt)**2

      tmp = -rnu * k2
      k2e = S * exp(tmp*dto)

      ! ox,oy,oz are the initial velocity components now
      vx = ox*k2e + ek*vx + delta*wx 
      vy = oy*k2e + ek*vy + delta*wy 
      vz = oz*k2e + ek*vz + delta*wz 

      fi = of * k2e + ek * fi + delta * wf ! Pr = 1 assumed
 
      k2e = exp(tmp * dt)

      vx = vx*k2e
      vy = vy*k2e
      vz = vz*k2e

      fi = fi * k2e

      time = time + dt
      dto = dt + dto
      call input(ox,oy,oz,of,0,lx1,ly,lz)
 
    case (2:)
      dtt=dt/dto
      S=(1.+.5*dtt)*dt
      ek=-.5*dtt*dt

      tmp = -rnu * k2
      k2e = ek * exp(tmp * dto)

      vx = vx + S*wx + ox*k2e
      vy = vy + S*wy + oy*k2e
      vz = vz + S*wz + oz*k2e

      fi = fi + S * wf + of * k2e
 
      k2e = exp(tmp * dt)

      vx = vx*k2e
      vy = vy*k2e
      vz = vz*k2e

      fi = fi * k2e
 
      ox = wx
      oy = wy
      oz = wz

      of = wf
 
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
        fi = 0._sp
    endwhere
    vx(1,1,1) = 0._sp
    vy(1,1,1) = 0._sp
    vz(1,1,1) = 0._sp
    fi(1,1,1) = 0._sp

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, fi, wx, wy, wz, wf, ox, oy, oz, of)
  deallocate(kx, ky, kz, k2, k2e, eks, vxt, vyt, vzt, tmp)

  write(*,*) 'ScalarKolFlows.x finished '
  stop

end program ScalarKolFlows
