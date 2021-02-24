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
  logical  :: new, newbij, forced
  ! parameters should be input from 'parameter_dns_deform.d' initially.

  complex(sp), allocatable, dimension (:,:,:) :: vx,vy,vz
  complex(sp), allocatable, dimension (:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension (:,:,:) :: dux,duy,duz
  complex(sp), allocatable, dimension (:,:,:) :: kux,kuy,kuz

  complex(sp), allocatable, dimension (:,:,:) :: b11, b12, b13
  complex(sp), allocatable, dimension (:,:,:) :: b21, b22, b23
  complex(sp), allocatable, dimension (:,:,:) :: b31, b32, b33
  complex(sp), allocatable, dimension (:,:,:) :: wb11, wb12, wb13
  complex(sp), allocatable, dimension (:,:,:) :: wb21, wb22, wb23
  complex(sp), allocatable, dimension (:,:,:) :: wb31, wb32, wb33
  complex(sp), allocatable, dimension (:,:,:) :: mb11, mb12, mb13
  complex(sp), allocatable, dimension (:,:,:) :: mb21, mb22, mb23
  complex(sp), allocatable, dimension (:,:,:) :: mb31, mb32, mb33
  complex(sp), allocatable, dimension (:,:,:) :: db11, db12, db13
  complex(sp), allocatable, dimension (:,:,:) :: db21, db22, db23
  complex(sp), allocatable, dimension (:,:,:) :: db31, db32, db33

  complex(sp), allocatable, dimension (:,:,:) :: lbt,tmp ! Two arrays used as temporary storage

  real(sp),    allocatable, dimension (:,:,:) :: k2,k2e
  real(sp),    allocatable, dimension (:)     :: kx,ky,kz, eks

  integer  :: lx, ly, lz, lx1
  integer  :: ii, jj, kk, ifile, istep, ll
  real(sp) :: kcut, kcut2, S, ek, const, ignore_me

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
    read(90,*) newbij
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
  write(*,*) 'newbij', newbij
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

  const = 1./(nx*ny*nz)
  kcut = 2.*lx/3
  kcut2 = kcut * kcut

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kux(lx1,ly,lz),kuy(lx1,ly,lz),kuz(lx1,ly,lz))
  allocate(dux(lx1,ly,lz),duy(lx1,ly,lz),duz(lx1,ly,lz))
  allocate(eks(lx),k2e(lx1,ly,lz))
  allocate( b11(lx1,ly,lz),  b12(lx1,ly,lz),  b13(lx1,ly,lz))
  allocate( b21(lx1,ly,lz),  b22(lx1,ly,lz),  b23(lx1,ly,lz)) 
  allocate( b31(lx1,ly,lz),  b32(lx1,ly,lz),  b33(lx1,ly,lz))
  allocate(mb11(lx1,ly,lz), mb12(lx1,ly,lz), mb13(lx1,ly,lz))
  allocate(mb21(lx1,ly,lz), mb22(lx1,ly,lz), mb23(lx1,ly,lz)) 
  allocate(mb31(lx1,ly,lz), mb32(lx1,ly,lz), mb33(lx1,ly,lz))
  allocate(db11(lx1,ly,lz), db12(lx1,ly,lz), db13(lx1,ly,lz))
  allocate(db21(lx1,ly,lz), db22(lx1,ly,lz), db23(lx1,ly,lz)) 
  allocate(db31(lx1,ly,lz), db32(lx1,ly,lz), db33(lx1,ly,lz))
  allocate(wb11(lx1,ly,lz), wb12(lx1,ly,lz), wb13(lx1,ly,lz))
  allocate(wb21(lx1,ly,lz), wb22(lx1,ly,lz), wb23(lx1,ly,lz)) 
  allocate(wb31(lx1,ly,lz), wb32(lx1,ly,lz), wb33(lx1,ly,lz))

  allocate(lbt(lx1,ly,lz),tmp(lx1,ly,lz))

  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  k2e = exp(-rnu*k2*dt/2)

  if (new) then
    ! Generate initial condition -------
    call initialize(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)

    write(*,*) 'after initialization of velocity '
  else
    ! Reading from disk ------
    call input(vx,vy,vz,idp,lx1,ly,lz)
    write(*,*) 'initial velocity field read.'
  endif

  if (newbij) then
    ! initialize bij
    b11 = (1,1); b12 = 0; b13 = 0
    b21 = 0; b22 = (1,1); b23 = 0
    b31 = 0; b32 = 0; b33 = (1,1)

    b11(3,4,5) = (4., 4.)
    b22(3,4,5) = (.5, .5)
    b33(3,4,5) = (.5, .5)

    write(*,*) 'after initializing bij'
  else
    call input(b11,b12,b13,b21,b22,b23,b31,b32,b33,idp,lx1,ly,lz)
    write(*,*) 'initial bij read.' 
  end if

  
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
      call output(b11,b12,b13,b21,b22,b23,b31,b32,b33,ifile,lx1,ly,lz)
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

    kux = vx; kuy = vy; kuz = vz
    mb11 = b11; mb12 = b12; mb13 = b13
    mb21 = b21; mb22 = b22; mb23 = b23
    mb31 = b31; mb32 = b32; mb33 = b33

    ! Calulate the rate of change of velocity 
    call ratesdudt
    ! kui are the same on return

    ! Calculate the stretching term in the equation of bij
    call dbdtstretch
    ! data in mbij are destroyed in dbdtstretch
    ! mbij returns as its Fourier transform (real values.
    ! wbij is the stretching term in the equation of bij

    ! Calculate the convection term in the equation of bij
    call renormalize ! renormalize the volume of bij
    call dbdtconvec(mb11,wb11)
    call dbdtconvec(mb12,wb12)
    call dbdtconvec(mb13,wb13)
    call dbdtconvec(mb21,wb21)
    call dbdtconvec(mb22,wb22)
    call dbdtconvec(mb23,wb23)
    call dbdtconvec(mb31,wb31)
    call dbdtconvec(mb32,wb32)
    call dbdtconvec(mb33,wb33)

    kux = dt * wx; kuy = dt * wy; kuz = dt * wz

    mb11 = dt * wb11; mb12 = dt * wb12; mb13 = dt * wb13
    mb21 = dt * wb21; mb22 = dt * wb22; mb23 = dt * wb23
    mb31 = dt * wb31; mb32 = dt * wb32; mb33 = dt * wb33

    ! increases in velocity
    dux = kux; duy = kuy; duz = kuz

    ! increase in bij
    db11 = mb11; db12 = mb12; db13 = mb13
    db21 = mb21; db22 = mb22; db23 = mb23
    db31 = mb31; db32 = mb32; db33 = mb33

    ! -- Second substep: k2u --

    ! un + k1u/2 to calculate k2u
    kux = k2e * (vx + kux/2)
    kuy = k2e * (vy + kuy/2)
    kuz = k2e * (vz + kuz/2)

    mb11 = b11 + mb11 / 2; mb12 = b12 + mb12 / 2; mb13 = b13 + mb13 / 2
    mb21 = b21 + mb21 / 2; mb22 = b22 + mb22 / 2; mb23 = b23 + mb23 / 2
    mb31 = b31 + mb31 / 2; mb32 = b32 + mb32 / 2; mb33 = b33 + mb33 / 2

    call ratesdudt
    ! kux kuy kuz are not changed on return

    ! Calculate the stretching term in the equation of bij
    call dbdtstretch
    ! mbij returns as its Fourier transform (real values.
    ! wbij is the stretching term in the equation of bij

    ! Calculate the convection term in the equation of bij
    call renormalize ! renormalize the volume of bij
    call dbdtconvec(mb11,wb11)
    call dbdtconvec(mb12,wb12)
    call dbdtconvec(mb13,wb13)
    call dbdtconvec(mb21,wb21)
    call dbdtconvec(mb22,wb22)
    call dbdtconvec(mb23,wb23)
    call dbdtconvec(mb31,wb31)
    call dbdtconvec(mb32,wb32)
    call dbdtconvec(mb33,wb33)

    kux = dt * wx / k2e; kuy = dt * wy / k2e; kuz = dt * wz / k2e

    mb11 = dt * wb11; mb12 = dt * wb12; mb13 = dt * wb13
    mb21 = dt * wb21; mb22 = dt * wb22; mb23 = dt * wb23
    mb31 = dt * wb31; mb32 = dt * wb32; mb33 = dt * wb33

    ! increases in velocity
    dux = dux + 2 * kux; duy = duy + 2 * kuy; duz = duz + 2 * kuz

    ! increase in bij
    db11 = db11 + 2 * mb11; db12 = db12 + 2 * mb12; db13 = db13 + 2 * mb13
    db21 = db21 + 2 * mb21; db22 = db22 + 2 * mb22; db23 = db23 + 2 * mb23
    db31 = db31 + 2 * mb31; db32 = db32 + 2 * mb32; db33 = db33 + 2 * mb33

    ! -- Third substep: k3u -- 

    ! un + k2u/2 to calculate k3u
    kux = k2e * (vx + kux/2); kuy = k2e * (vy + kuy/2); kuz = k2e * (vz + kuz/2)

    mb11 = b11 + mb11 / 2; mb12 = b12 + mb12 / 2; mb13 = b13 + mb13 / 2
    mb21 = b21 + mb21 / 2; mb22 = b22 + mb22 / 2; mb23 = b23 + mb23 / 2
    mb31 = b31 + mb31 / 2; mb32 = b32 + mb32 / 2; mb33 = b33 + mb33 / 2

    call ratesdudt
    ! kux kuy kuz are not changed on return

    ! Calculate the stretching term in the equation of bij
    call dbdtstretch
    ! mbij returns as its Fourier transform (real values.
    ! wbij is the stretching term in the equation of bij

    ! Calculate the convection term in the equation of bij
    call renormalize ! renormalize the volume of bij
    call dbdtconvec(mb11,wb11)
    call dbdtconvec(mb12,wb12)
    call dbdtconvec(mb13,wb13)
    call dbdtconvec(mb21,wb21)
    call dbdtconvec(mb22,wb22)
    call dbdtconvec(mb23,wb23)
    call dbdtconvec(mb31,wb31)
    call dbdtconvec(mb32,wb32)
    call dbdtconvec(mb33,wb33)

    kux = dt * wx / k2e; kuy = dt * wy / k2e; kuz = dt * wz / k2e

    mb11 = dt * wb11; mb12 = dt * wb12; mb13 = dt * wb13
    mb21 = dt * wb21; mb22 = dt * wb22; mb23 = dt * wb23
    mb31 = dt * wb31; mb32 = dt * wb32; mb33 = dt * wb33

    ! increases in velocity
    dux = dux + 2 * kux; duy = duy + 2 * kuy; duz = duz + 2 * kuz

    ! increase in bij
    db11 = db11 + 2 * mb11; db12 = db12 + 2 * mb12; db13 = db13 + 2 * mb13
    db21 = db21 + 2 * mb21; db22 = db22 + 2 * mb22; db23 = db23 + 2 * mb23
    db31 = db31 + 2 * mb31; db32 = db32 + 2 * mb32; db33 = db33 + 2 * mb33

    ! -- Fourth substep: k4u --

    ! un + k3u to calculate k4u
    kux = k2e * k2e * (vx + kux) 
    kuy = k2e * k2e * (vy + kuy)
    kuz = k2e * k2e * (vz + kuz)

    mb11 = b11 + mb11; mb12 = b12 + mb12; mb13 = b13 + mb13
    mb21 = b21 + mb21; mb22 = b22 + mb22; mb23 = b23 + mb23
    mb31 = b31 + mb31; mb32 = b32 + mb32; mb33 = b33 + mb33

    call ratesdudt
    ! kux kuy kuz are the not changed on return

    ! Calculate the stretching term in the equation of bij
    call dbdtstretch    
    ! mbij returns as its Fourier transform (real values.
    ! wbij is the stretching term in the equation of bij

    ! Calculate the convection term in the equation of bij
    call renormalize ! renormalize the volume of bij
    call dbdtconvec(mb11,wb11)
    call dbdtconvec(mb12,wb12)
    call dbdtconvec(mb13,wb13)
    call dbdtconvec(mb21,wb21)
    call dbdtconvec(mb22,wb22)
    call dbdtconvec(mb23,wb23)
    call dbdtconvec(mb31,wb31)
    call dbdtconvec(mb32,wb32)
    call dbdtconvec(mb33,wb33)

    kux = dt * wx / k2e / k2e
    kuy = dt * wy / k2e / k2e
    kuz = dt * wz / k2e / k2e

    mb11 = dt * wb11; mb12 = dt * wb12; mb13 = dt * wb13
    mb21 = dt * wb21; mb22 = dt * wb22; mb23 = dt * wb23
    mb31 = dt * wb31; mb32 = dt * wb32; mb33 = dt * wb33

    ! increases in velocity
    dux = dux + kux; duy = duy + kuy; duz = duz + kuz

    ! increase in bij
    db11 = db11 + mb11; db12 = db12 + mb12; db13 = db13 + mb13
    db21 = db21 + mb21; db22 = db22 + mb22; db23 = db23 + mb23
    db31 = db31 + mb31; db32 = db32 + mb32; db33 = db33 + mb33

    ! -- Velocity at next time step --
    vx = k2e * k2e * ( vx + dux / 6 )
    vy = k2e * k2e * ( vy + duy / 6 )
    vz = k2e * k2e * ( vz + duz / 6 )

    b11 = b11 + db11 / 6; b12 = b12 + db12 / 6; b13 = b13 + db13 / 6
    b21 = b21 + db21 / 6; b22 = b22 + db22 / 6; b23 = b23 + db23 / 6
    b31 = b31 + db31 / 6; b32 = b32 + db32 / 6; b33 = b33 + db33 / 6

    ! == End of RK4 stepping ==

    ! Dealiasing
    where (k2 .ge. kcut2)
        vx = 0. ; vy = 0. ; vz = 0.
    endwhere

    time = time + dt

  end do
  !=================== end of loop ==============================

  close(24)
  close(25)

  ! Deallocate arrays
  deallocate(vx, vy, vz, wx, wy, wz, kux, kuy, kuz, dux, duy, duz)
  deallocate(kx, ky, kz, k2, k2e, eks)

  deallocate(b11,b12,b13,b21,b22,b23,wb11,wb12,wb13,wb21,wb22,wb23)
  deallocate(wb31,wb32,wb33,db11,db12,db13,db21,db22,db23,db31,db32,db33)
  deallocate(mb11,mb12,mb13,mb21,mb22,mb23,mb31,mb32,mb33,lbt,tmp)

  ! Destroy plans
  call destroyplan3d

  write(*,*)'finished '
  stop

contains

  subroutine ratesdudt
    
    integer :: ii, jj, kk

    ! =========== Calculate vorticity =============================
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk) = eye * ky(jj) * kuz(ii,jj,kk) - eye * kz(kk) * kuy(ii,jj,kk)
      wy(ii,jj,kk) = eye * kz(kk) * kux(ii,jj,kk) - eye * kx(ii) * kuz(ii,jj,kk)
      wz(ii,jj,kk) = eye * kx(ii) * kuy(ii,jj,kk) - eye * ky(jj) * kux(ii,jj,kk)
    end do
    end do
    end do
 
    !========== calculate conkuection term: lamb kuector =======================
    call convec_dns(kux,kuy,kuz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
 
    !========== combine conkuection term, and forcing ============
    if (forced) call force(wx,wy,wz,kux,kuy,kuz,k2,lx1,ly,lz,fkmax,eps)
 
    call symmetrize(wx,lx1,ly,lz)
    call symmetrize(wy,lx1,ly,lz)
    call symmetrize(wz,lx1,ly,lz)

    !========== projection ================================= 
    call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
 
  end subroutine ratesdudt

  subroutine dbdtstretch

    integer :: ii, jj, kk

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wb11(ii,jj,kk) = eye * kx(ii) * kux(ii,jj,kk)
      wb12(ii,jj,kk) = eye * ky(jj) * kux(ii,jj,kk)
      wb13(ii,jj,kk) = eye * kz(kk) * kux(ii,jj,kk)
      wb21(ii,jj,kk) = eye * kx(ii) * kuy(ii,jj,kk)
      wb22(ii,jj,kk) = eye * ky(jj) * kuy(ii,jj,kk)
      wb23(ii,jj,kk) = eye * kz(kk) * kuy(ii,jj,kk)
      wb31(ii,jj,kk) = eye * kx(ii) * kuz(ii,jj,kk)
      wb32(ii,jj,kk) = eye * ky(jj) * kuz(ii,jj,kk)
      wb33(ii,jj,kk) = eye * kz(kk) * kuz(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wb11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wb33,ignore_me)
    
    lbt  = cmplx( real(wb11) * real(mb11) &
                + real(wb12) * real(mb21) &
                + real(wb13) * real(mb31),  &
                  aimag(wb11) * aimag(mb11) &
                + aimag(wb12) * aimag(mb21) &
                + aimag(wb13) * aimag(mb31) )  !wb11
    tmp  = cmplx( real(wb11) * real(mb12) &
                + real(wb12) * real(mb22) &
                + real(wb13) * real(mb32),  &
                  aimag(wb11) * aimag(mb12) &
                + aimag(wb12) * aimag(mb22) &
                + aimag(wb13) * aimag(mb32) )  !wb12
    wb13  = cmplx( real(wb11) * real(mb13) &
                 + real(wb12) * real(mb23) &
                 + real(wb13) * real(mb33),  &
                   aimag(wb11) * aimag(mb13) &
                 + aimag(wb12) * aimag(mb23) &
                 + aimag(wb13) * aimag(mb33) )  !wb13
    wb11 = lbt; wb12 = tmp

    lbt  = cmplx( real(wb21) * real(mb11) &
                + real(wb22) * real(mb21) &
                + real(wb23) * real(mb31),  &
                  aimag(wb21) * aimag(mb11) &
                + aimag(wb22) * aimag(mb21) &
                + aimag(wb23) * aimag(mb31) )  !wb21
    tmp  = cmplx( real(wb21) * real(mb12) &
                + real(wb22) * real(mb22) &
                + real(wb23) * real(mb32),  &
                  aimag(wb21) * aimag(mb12) &
                + aimag(wb22) * aimag(mb22) &
                + aimag(wb23) * aimag(mb32) )  !wb22
    wb23  = cmplx( real(wb21) * real(mb13) &
                 + real(wb22) * real(mb23) &
                 + real(wb23) * real(mb33),  &
                   aimag(wb21) * aimag(mb13) &
                 + aimag(wb22) * aimag(mb23) &
                 + aimag(wb23) * aimag(mb33) )  !wb23
    wb21 = lbt; wb22 = tmp

    lbt  = cmplx( real(wb31) * real(mb11) &
                + real(wb32) * real(mb21) &
                + real(wb33) * real(mb31),  &
                  aimag(wb31) * aimag(mb11) &
                + aimag(wb32) * aimag(mb21) &
                + aimag(wb33) * aimag(mb31) )  !wb31
    tmp  = cmplx( real(wb31) * real(mb12) &
                + real(wb32) * real(mb22) &
                + real(wb33) * real(mb32),  &
                  aimag(wb31) * aimag(mb12) &
                + aimag(wb32) * aimag(mb22) &
                + aimag(wb33) * aimag(mb32) )  !wb32
    wb33  = cmplx( real(wb31) * real(mb13) &
                 + real(wb32) * real(mb23) &
                 + real(wb33) * real(mb33),  &
                   aimag(wb31) * aimag(mb13) &
                 + aimag(wb32) * aimag(mb23) &
                 + aimag(wb33) * aimag(mb33) )  !wb23
    wb31 = lbt; wb32 = tmp



  end subroutine dbdtstretch

  subroutine dbdtconvec(lbij, lwbij)

    ! lbij is real when entering this subroutine

    complex(sp), dimension(:,:,:) :: lbij, lwbij
    complex(sp), dimension(lx1,ly,lz) :: tmp2

!    lbt = lbij
!
!    tmp = kux
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp = const * cmplx( real(tmp) * real(lbt), aimag(tmp) * aimag(lbt)) 
!    call rfftwnd_f77_one_real_to_complex(r2c3d,tmp,ignore_me)
!    do ii = 1, lx1
!      tmp(ii,:,:) = eye * kx(ii) * tmp(ii,:,:)
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    lwbij = lwbij - tmp 
!
!    tmp = kuy
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp = const * cmplx( real(tmp) * real(lbt), aimag(tmp) * aimag(lbt)) 
!    call rfftwnd_f77_one_real_to_complex(r2c3d,tmp,ignore_me)
!    do jj = 1, ly
!      tmp(:,jj,:) = eye * ky(jj) * tmp(:,jj,:)
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    lwbij = lwbij - tmp
!
!    tmp = kuz
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp = const * cmplx( real(tmp) * real(lbt), aimag(tmp) * aimag(lbt)) 
!    call rfftwnd_f77_one_real_to_complex(r2c3d,tmp,ignore_me)
!    do kk = 1, lz
!      tmp(:,:,kk) = eye * kz(kk) * tmp(:,:,kk)
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    lwbij = lwbij - tmp

!    lbt = lbij
!    call rfftwnd_f77_one_real_to_complex(r2c3d,lbt,ignore_me)
!    lbt = lbt * const
!
!    ! zero ky=-ly/2 modes
!    lbt(:,ly/2+1,:) = (0.0_SP,0.0_SP)
!    ! zero kz=-lz/2 modes
!    lbt(:,:,lz/2+1) = (0.0_SP,0.0_SP)
!    ! zero kx=lx modes
!    lbt(lx1,:,:)=(0.0_SP,0.0_SP)
!
!    do ii = 1, lx1
!      tmp(ii,:,:) = eye * kx(ii) * lbt(ii,:,:) 
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp2 = kux 
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp2,ignore_me)
!
!    lwbij = lwbij - cmplx( real(tmp) * real(tmp2), aimag(tmp) * aimag(tmp2) )
!
!    do jj = 1, ly
!      tmp(:,jj,:) = eye * ky(jj) * lbt(:,jj,:)
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp2 = kuy
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp2,ignore_me)
!
!    lwbij = lwbij - cmplx( real(tmp) * real(tmp2), aimag(tmp) * aimag(tmp2) )
!
!    do kk = 1, lz
!      tmp(:,:,kk) = eye * kz(kk) * lbt(:,:,kk)
!    end do
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
!    tmp2 = kuz 
!    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp2,ignore_me)
!
!    lwbij = lwbij - cmplx( real(tmp) * real(tmp2), aimag(tmp) * aimag(tmp2) )

  end subroutine dbdtconvec

  subroutine renormalize

    integer :: ii, jj, kk, ll
    real(sp) :: vecx1, vecy1, vecz1
    real(sp) :: vecx2, vecy2, vecz2
    real(sp) :: vecx3, vecy3, vecz3
    real(sp) :: vol

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx 

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2

          vecx1 = real(mb11(ll,jj,kk))
          vecy1 = real(mb21(ll,jj,kk))
          vecz1 = real(mb31(ll,jj,kk))
          vecx2 = real(mb12(ll,jj,kk))
          vecy2 = real(mb22(ll,jj,kk))
          vecz2 = real(mb32(ll,jj,kk))
          vecx3 = real(mb13(ll,jj,kk))
          vecy3 = real(mb23(ll,jj,kk))
          vecz3 = real(mb33(ll,jj,kk))
      else
          ll = ii/2

          vecx1 = aimag(mb11(ll,jj,kk))
          vecy1 = aimag(mb21(ll,jj,kk))
          vecz1 = aimag(mb31(ll,jj,kk))
          vecx2 = aimag(mb12(ll,jj,kk))
          vecy2 = aimag(mb22(ll,jj,kk))
          vecz2 = aimag(mb32(ll,jj,kk))
          vecx3 = aimag(mb13(ll,jj,kk))
          vecy3 = aimag(mb23(ll,jj,kk))
          vecz3 = aimag(mb33(ll,jj,kk))
      end if

      vol = vecx1 * (vecy2 * vecz3 - vecz2 * vecy3) &
          + vecy1 * (vecz2 * vecx3 - vecx2 * vecz3) &
          + vecz1 * (vecx2 * vecy3 - vecy2 * vecx3)

      vol = vol**(1./3)
      vecx1 = vecx1 / vol; vecy1 = vecy1 / vol; vecz1 = vecz1 / vol
      vecx2 = vecx2 / vol; vecy2 = vecy2 / vol; vecz2 = vecz2 / vol
      vecx3 = vecx3 / vol; vecy3 = vecy3 / vol; vecz3 = vecz3 / vol

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2

          mb11(ll,jj,kk) = cmplx( vecx1, aimag(mb11(ll,jj,kk)) )
          mb21(ll,jj,kk) = cmplx( vecy1, aimag(mb21(ll,jj,kk)) ) 
          mb31(ll,jj,kk) = cmplx( vecz1, aimag(mb31(ll,jj,kk)) ) 
          mb12(ll,jj,kk) = cmplx( vecx2, aimag(mb12(ll,jj,kk)) )
          mb22(ll,jj,kk) = cmplx( vecy2, aimag(mb22(ll,jj,kk)) ) 
          mb32(ll,jj,kk) = cmplx( vecz2, aimag(mb32(ll,jj,kk)) ) 
          mb13(ll,jj,kk) = cmplx( vecx3, aimag(mb13(ll,jj,kk)) )
          mb23(ll,jj,kk) = cmplx( vecy3, aimag(mb23(ll,jj,kk)) ) 
          mb33(ll,jj,kk) = cmplx( vecz3, aimag(mb33(ll,jj,kk)) ) 

      else
          ll = ii/2

          mb11(ll,jj,kk) = cmplx( real(mb11(ll,jj,kk)), vecx1 )
          mb21(ll,jj,kk) = cmplx( real(mb21(ll,jj,kk)), vecy1 ) 
          mb31(ll,jj,kk) = cmplx( real(mb31(ll,jj,kk)), vecz1 ) 
          mb12(ll,jj,kk) = cmplx( real(mb12(ll,jj,kk)), vecx2 )
          mb22(ll,jj,kk) = cmplx( real(mb22(ll,jj,kk)), vecy2 ) 
          mb32(ll,jj,kk) = cmplx( real(mb32(ll,jj,kk)), vecz2 ) 
          mb13(ll,jj,kk) = cmplx( real(mb13(ll,jj,kk)), vecx3 )
          mb23(ll,jj,kk) = cmplx( real(mb23(ll,jj,kk)), vecy3 ) 
          mb33(ll,jj,kk) = cmplx( real(mb33(ll,jj,kk)), vecz3 ) 

      end if

    end do
    end do
    end do

  end subroutine renormalize

!  subroutine dbdtconvec_det
!    
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb11,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb12,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb13,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb21,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb22,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb23,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb31,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb32,ignore_me)
!    call rfftwnd_f77_one_real_to_complex(r2c3d,mb33,ignore_me)
!    mb11 = mb11 * const; mb12 = mb12 * const; mb13 = mb13 * const
!    mb21 = mb21 * const; mb22 = mb22 * const; mb23 = mb23 * const
!    mb31 = mb31 * const; mb32 = mb32 * const; mb33 = mb33 * const
!    
!    ! projection to make the determinant constant
!
!    do kk = 1, lz
!    do jj = 1, ly
!    do ii = 1, nx + 2
!
!      if ( mod(ii,2) .eq. 1) then
!          ll = (ii + 1) / 2
!
!          vec1x = real(mb11(ll,jj,kk))
!          vec1y = real(mb21(ll,jj,kk))
!          vec1z = real(mb31(ll,jj,kk))
!          vec2x = real(mb12(ll,jj,kk))
!          vec2y = real(mb22(ll,jj,kk))
!          vec2z = real(mb32(ll,jj,kk))
!          vec3x = real(mb13(ll,jj,kk))
!          vec3y = real(mb23(ll,jj,kk))
!          vec3z = real(mb33(ll,jj,kk))
!      else
!          ll = ii/2
!          vec1x = aimag(mb11(ll,jj,kk))
!          vec1y = aimag(mb21(ll,jj,kk))
!          vec1z = aimag(mb31(ll,jj,kk))
!          vec2x = aimag(mb12(ll,jj,kk))
!          vec2y = aimag(mb22(ll,jj,kk))
!          vec2z = aimag(mb32(ll,jj,kk))
!          vec3x = aimag(mb13(ll,jj,kk))
!          vec3y = aimag(mb23(ll,jj,kk))
!          vec3z = aimag(mb33(ll,jj,kk))
!      end if
!
!      cross23x = vec2y * vec3z - vec2z * vec3y
!      cross23y = vec2z * vec3x - vec2x * vec3z
!      cross23z = vec2x * vec3y - vec2y * vec3x
!      normc23 = sqrt(cross23x**2 + cross23y**2 + cross23z**2)
!      cross23x = cross23x / normc23
!      cross23y = cross23y / normc23
!      cross23z = cross23z / normc23
!
!      normc23 = vec1x * cross23x + vec1y * cross23y + vec1z * cross23z
!
!      vec1x = vec1x - normc23 * cross23x
!      vec1y = vec1y - normc23 * cross23y
!      vec1z = vec1z - normc23 * cross23z
!
!      if ( mod(ii,2) .eq. 1) then
!          ll = (ii + 1) / 2
!
!          mb11(ll,jj,kk) = cmplx( vec1x, aimag(mb11(ll,jj,kk)) )
!          mb21(ll,jj,kk) = cmplx( vec1y, aimag(mb11(ll,jj,kk)) ) 
!          mb31(ll,jj,kk) = cmplx( vec1z, aimag(mb11(ll,jj,kk)) ) 
!
!      else
!          ll = ii/2
!
!          mb11(ll,jj,kk) = cmplx( aimag(mb11(ll,jj,kk)), vec1x )
!          mb21(ll,jj,kk) = cmplx( aimag(mb11(ll,jj,kk)), vec1y ) 
!          mb31(ll,jj,kk) = cmplx( aimag(mb11(ll,jj,kk)), vec1z ) 
!
!      end if
!
!    end do
!    end do
!    end do
!
!  end subroutine dbdtconvec_det

end program dnsdeform












