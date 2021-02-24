module mconstant
  integer,     parameter :: ppp=precision(1.d0)
  integer,     parameter :: sp=selected_real_kind(ppp)
  complex(sp), parameter :: eye=(0._sp,1._sp)
  real(sp),    parameter :: mytiny=tiny(1.0_sp)
  real(sp),    parameter :: pi = 3.1415926536_sp
end module mconstant

module mfftwheader
  integer, parameter :: fftw_r2hc=0, fftw_hc2r=1, fftw_dht=2
  integer, parameter :: fftw_redft00=3, fftw_redft01=4, fftw_redft10=5
  integer, parameter :: fftw_redft11=6, fftw_rodft00=7, fftw_rodft01=8
  integer, parameter :: fftw_rodft10=9, fftw_rodft11=10, fftw_forward=-1
  integer, parameter :: fftw_backward=+1, fftw_measure=0, fftw_destroy_input=1
  integer, parameter :: fftw_unaligned=2, fftw_conserve_memory=4, fftw_exhaustive=8
  integer, parameter :: fftw_preserve_input=16, fftw_patient=32, fftw_estimate=64
  integer, parameter :: fftw_estimate_patient=128, fftw_believe_pcost=256
  integer, parameter :: fftw_no_dft_r2hc=512, fftw_no_nonthreaded=1024 
  integer, parameter :: fftw_no_buffering=2048, fftw_no_indirect_op=4096
  integer, parameter :: fftw_allow_large_generic=8192, fftw_no_rank_splits=16384
  integer, parameter :: fftw_no_vrank_splits=32768, fftw_no_vrecurse=65536
  integer, parameter :: fftw_no_simd=131072, fftw_no_slow=262144
  integer, parameter :: fftw_no_fixed_radix_large_n=524288
  integer, parameter :: fftw_allow_pruning=1048576
end module mfftwheader

subroutine wavenumber2d(kx, ky, k2, lx1, ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly

  real(sp), dimension(lx1), intent(out) :: kx
  real(sp), dimension(ly),  intent(out) :: ky
  real(sp), dimension(lx1,ly), intent(out) :: k2

  integer :: ii, jj

  do ii=1,lx1
    kx(ii)=real(ii-1)
  end do

  do jj=1,ly
    ! ky=0, 1, ..., ly/2-1, -ly/2, -ly/2+1, ..., -1
    ! Ref. fftw_manual.pdf page 22
    ky(jj) = real(mod(jj-1+ly/2,ly)-ly/2,sp)
  end do


  do jj=1,ly
  do ii=1,lx1
    k2(ii,jj)=kx(ii)*kx(ii)+ky(jj)*ky(jj)
  end do
  end do
  k2(1,1)=mytiny


  return
end subroutine wavenumber2d

subroutine output2d(mb,tmprt,idump,lx1,ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly
  complex(sp), dimension(lx1,ly), intent(in) :: mb, tmprt

  integer      :: idump
  character*50 :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)

  fpath='./out/mb'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  write(10)mb
  close(10)

  fpath='./out/tmprt'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  write(10)tmprt
  close(10)

  return
end subroutine output2d

subroutine input2d(mb,tmprt,idump,lx1,ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly
  complex(sp), dimension(lx1,ly), intent(out) :: mb, tmprt

  integer      :: idump
  character*50 :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)

  fpath='./out/mb'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)mb
  close(10)

  fpath='./out/tmprt'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)tmprt
  close(10)

  return
end subroutine input2d

subroutine initbt (mb, tmprt, k2, lx1, ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1, ly
  complex(sp), dimension(lx1,ly), intent(out) :: mb, tmprt
  real(sp),    dimension(lx1,ly), intent(in)  :: k2

  real(sp), parameter :: bam = 0.005, fkmax2 = 50. * 50.

  integer :: ranseed

  real :: ran2

  ranseed = -197

  ! temperature fluctuation is initially zero.
  tmprt = 0. 

  ! initializing magnetic field
  where ( k2 .le. fkmax2 ) 
      mb = bam * exp( eye * ran2(ranseed) * 2 * pi )
  elsewhere
      mb = 0.
  endwhere

end subroutine initbt

subroutine setzeros2d (c, k2, lx1, ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1, ly
  complex(sp), dimension(lx1,ly), intent(inout) :: c
  real(sp),    dimension(lx1,ly), intent(in)    :: k2

  real(sp) :: kcut2

  c(:,ly/2+1) = 0.
  c(lx1,:) = 0.
  c(1,1) = 0.

  kcut2 = ( (lx1 - 1)*2. /3.)**2 

  where ( k2 .ge. kcut2 ) c = 0.

end subroutine setzeros2d

subroutine hermitianize2d(c,lx1,ly)
  use mconstant
  implicit none

  integer, intent(in) :: lx1, ly
  complex(sp), dimension(lx1,ly), intent(inout) :: c

  integer :: iy

  ! imposing the constraint on the line kx=0.
  do iy=ly/2+2,ly
  c(1,iy) = CONJG(c(1,ly+2-iy))
  end do

  return
end subroutine hermitianize2d


real function ran2(idum)
!  Random generator from Numerical Recipes, p.271
implicit none
integer idum,ia,im,iq,ir,ntab,ndiv
real am,eps,rnmx
parameter(ia = 16807, im = 2147483647, am = 1./im, iq = 127773, ir = 2836,  &
          ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7, rnmx = 1. - eps)
integer j,k
integer, dimension(ntab), save :: iv=(/(0,j=1,ntab)/)
integer, save :: iy = 0


if (idum <= 0 .or. iy == 0) then
  idum = max(-idum,1)
  do j = ntab + 8, 1, -1 
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k
    if (idum < 0) idum = idum + im
    if (j <= ntab) iv(j) = idum
  end do
  iy = iv(1)
end if

k = idum/iq
idum = ia*(idum - k*iq) - ir*k
if (idum < 0) idum = idum + im
j = 1 + iy/ndiv
iy = iv(j)
iv(j) = idum
ran2 = min(am*iy,rnmx)

End function


program medv
  use mconstant
  use mfftwheader
  implicit none

  ! FFTW plan
  integer(sp) :: IFTb, IFTdxb, IFTdyb, IFTdxt, IFTdyt
  integer(sp) :: FTb, FTdxb, FTdxt

  ! Physical parameters
  real(sp), parameter :: sigma = -1.  ! sigma = 1 or -1 only.
  real(sp), parameter :: kappa = 0.  ! kappa = 0 set the scalar nonlinear term off

  integer :: nx, ny, lx, lx1, ly, ifile, iscreen, istep
  integer :: ii, jj
  real(sp) :: time, starttime, timemax, dtw, dt, dt_h, twrite
  real(sp) :: const, etot, emb, etmprt, enstrophy

  logical :: newrun

  complex(sp), allocatable, dimension(:,:) :: tmprt, mb
  complex(sp), allocatable, dimension(:,:) :: dxb, dyb
  complex(sp), allocatable, dimension(:,:) :: dxt, dyt
  complex(sp), allocatable, dimension(:,:) :: fb, ft
  complex(sp), allocatable, dimension(:,:) :: fbo, fto
  real(sp),    allocatable, dimension(:,:) :: k2
  real(sp),    allocatable, dimension(:)   :: kx, ky

  open(90,file='parameter_medv.d',status='unknown')
    read(90,*) newrun
    read(90,*) nx
    read(90,*) ny
    read(90,*) dt
    read(90,*) starttime
    read(90,*) timemax
    read(90,*) dtw
    read(90,*) ifile
    read(90,*) iscreen
  close(90)


  lx = nx/2
  lx1 = lx + 1
  ly = ny
  const = 1._sp / (nx * ny)
  dt_h = dt/ 2.

  allocate ( tmprt(lx1,ly), mb(lx1,ly) )
  allocate ( dxb(lx1,ly), dyb(lx1,ly), dxt(lx1,ly), dyt(lx1,ly) )
  allocate ( fb(lx1,ly), ft(lx1,ly), fbo(lx1,ly), fto(lx1,ly) )
  allocate ( kx(lx1), ky(ly), k2(lx1,ly) )

  call dfftw_plan_dft_c2r_2d(IFTb,nx,ny,mb,mb,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_2d(IFTdxb,nx,ny,dxb,dxb,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_2d(IFTdyb,nx,ny,dyb,dyb,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_2d(IFTdxt,nx,ny,dxt,dxt,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_2d(IFTdyt,nx,ny,dyt,dyt,FFTW_MEASURE)

  call dfftw_plan_dft_r2c_2d(FTb,nx,ny,mb,mb,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_2d(FTdxb,nx,ny,dxb,dxb,FFTW_MEASURE)
  call dfftw_plan_dft_r2c_2d(FTdxt,nx,ny,dxt,dxt,FFTW_MEASURE)
  Write(*,*) 'FFTW plan done.'


  call wavenumber2d(kx, ky, k2, lx1, ly)
  Write(*,*) 'wavenumber2d done'

  if (newrun) then

      call initbt(mb, tmprt, k2, lx1, ly)
      call setzeros2d(mb, k2, lx1, ly)
      call setzeros2d(tmprt, k2, lx1, ly)
      call hermitianize2d(mb, lx1, ly)
      call hermitianize2d(tmprt, lx1, ly)
      write(*,*) 'Initialization done.'

  else

      call input2d(mb, tmprt, ifile, lx1, ly)
      write(*,*) 'Reading initial field done.'

  end if


  fbo = mb
  fto = tmprt

  istep = -1
  time = starttime
  twrite = time

  write(*,*) 'Starting loop...'
  do while (time .le. timemax)

    istep = istep + 1

    if ( abs(time - twrite) .le. dt_h ) then

        call output2d(fb,ft,ifile,lx1,ly)

        twrite = twrite + dtw
        ifile = ifile + 1
    end if
    !========= RHS of the equation for B ===========
 
    do jj = 1, ly
    do ii = 1, lx1
      ! linear terms on the RHS of the equations
      fb(ii,jj) =  - eye * ky(jj) * tmprt(ii,jj) / ( 1 + k2(ii,jj) )
 
      ! gradient of B
      dxb(ii,jj) = eye * kx(ii) * mb(ii,jj)
      dyb(ii,jj) = eye * ky(jj) * mb(ii,jj)
 
      ! gradient of Laplacian of B. dxt and dyt are used temporarily. 
      dxt(ii,jj) = - eye * kx(ii) * mb(ii,jj) * k2(ii,jj)
      dyt(ii,jj) = - eye * ky(jj) * mb(ii,jj) * k2(ii,jj)
    end do
    end do

    ! prepare output on the screen: energy
    if ( mod(istep, iscreen) .eq. 0 ) then 

        ft = mb * conjg(mb) + dxb * conjg(dxb) + dyb * conjg(dyb) 
        ft(1,:) = .5 * ft(1,:)

        emb = 2.*sum( real(ft) )

        ft = tmprt * conjg(tmprt) / sigma
        ft(1,:) = .5 * ft(1,:)
        
        etmprt = 2. * sum ( real(ft) )

        etot = emb + etmprt

    end if
 
    call dfftw_execute_dft_c2r(IFTdxb, dxb, dxb) ! dxB in real space
    call dfftw_execute_dft_c2r(IFTdyb, dyb, dyb) ! dyB in real space
    call dfftw_execute_dft_c2r(IFTdxt, dxt, dxt) ! dx LaplacianB in real space
    call dfftw_execute_dft_c2r(IFTdyt, dyt, dyt) ! dy LaplacianB in real space
 
    ! {B, Laplacian B} in Real space. dxt as temporary storage
    dxt = cmplx( real(dxb) * real(dyt) - real(dyb) * real(dxt), &
                 aimag(dxb) * aimag(dyt) - aimag(dyb) * aimag(dxt) )
 
    dxt = dxt * const
    call dfftw_execute_dft_r2c(FTdxt, dxt, dxt) ! {B, Laplacian B} in Fourier space
    fb = fb + dxt / ( 1 + k2 )
 
    call dfftw_execute_dft_c2r(IFTb, mb, mb) ! B in Real space
    dxt = cmplx( real(mb) * real(mb), aimag(mb) * aimag(mb) )
 
    dxt = dxt * const
    call dfftw_execute_dft_r2c(FTdxt, dxt, dxt) ! dxt in Fourier space
    do jj = 1, ly
      fb(:,jj) = fb(:,jj) - .5 * eye * kappa * ky(jj) * dxt(:,jj) / ( 1 + k2(:,jj) )
    end do
 
    mb = mb * const
    call dfftw_execute_dft_r2c(FTb, mb, mb) ! B back in Fourier space 
 
    !======== RHS of equation for temperature ===============
    do jj = 1, ly
    do ii = 1, lx1
 
      ! gradient of temperature
      dxt(ii,jj) = eye * kx(ii) * tmprt(ii,jj)
      dyt(ii,jj) = eye * ky(jj) * tmprt(ii,jj)
    end do
    end do
 
    call dfftw_execute_dft_c2r(IFTdxt, dxt, dxt) ! dxt in Real space
    call dfftw_execute_dft_c2r(IFTdyt, dyt, dyt) ! dyt in Real space

    ! screen output 
    if ( mod(istep, iscreen) .eq. 0 ) then

        ! ft used as temporary storage
        ft = conjg(mb) * tmprt
        ft(1,:) = .5 * ft(1,:)

        enstrophy = sum( real(ft) )

        enstrophy = enstrophy + sum( real(dxb(1:lx,:)) * real(dxt(1:lx,:)) ) * const &
                              + sum( aimag(dxb(1:lx,:)) * aimag(dxt(1:lx,:)) ) * const

        write(*,'(10E12.3)') time, emb, etmprt, etot, enstrophy

    end if

    ! The linear term
    do jj = 1, ly
      ft(:,jj) = - eye * sigma * ky(jj) * mb(:,jj)
    end do
 
    ! Nonlinear term
    dxb = cmplx( real(dxb) * real(dyt) - real(dyb) * real(dxt), &
                 aimag(dxb) * aimag(dyt) - aimag(dyb) * aimag(dxt) )
 
    dxb = dxb * const
    call dfftw_execute_dft_r2c(FTdxb, dxb, dxb) ! {B, T} in Fourier space
    ft = ft - dxb

    call setzeros2d(fb, k2, lx1, ly)
    call setzeros2d(ft, k2, lx1, ly)
    call hermitianize2d(fb, lx1, ly)
    call hermitianize2d(ft, lx1, ly)
 
    !=============== Advance by 2nd order Adams-Bashforth ============

    select case (istep)

    case (0) 
        tmprt = tmprt  + dt_h * ft 
        mb = mb + dt_h * fb

        time = time + dt_h
        call output2d(fb,ft,0,lx1,ly)

    case (1)
        tmprt = fto + dt * ft 
        mb = fbo + dt * fb

        time = time + dt_h
        call input2d(fbo,fto,0,lx1,ly)

    case (2:)
        tmprt = tmprt + dt_h * (3. * ft - fto)
        mb = mb + dt_h * (3. * fb - fbo)

        fto = ft; fbo = fb

        time = time + dt

    case default
        write(*,*) "wrong istep. istep = ", istep

    end select

  end do

  call dfftw_destroy_plan(IFTb)
  call dfftw_destroy_plan(IFTdxb)
  call dfftw_destroy_plan(IFTdyb)
  call dfftw_destroy_plan(IFTdxt)
  call dfftw_destroy_plan(IFTdyt)

  call dfftw_destroy_plan(FTb)
  call dfftw_destroy_plan(FTdxb)
  call dfftw_destroy_plan(FTdxt)

  deallocate(mb, tmprt, fb, ft, fbo, fto, dxb, dyb, dxt, dyt)
  deallocate(kx, ky, k2)

  write(*,*) 'Finished.'

end program medv
