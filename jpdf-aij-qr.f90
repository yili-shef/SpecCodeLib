program jpdfaijqr
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,mm,lx1,lx,ly,lz,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(sp), dimension(3,3) :: gg

  integer, parameter :: npnt = 160
  real(sp), parameter :: boundq = 8._sp, boundr = 8._sp
  real(sp), parameter :: binwq = 2*boundq/npnt, binwr = 2*boundr/npnt
  real(dp), dimension(npnt, npnt) :: pqr

  real(dp) :: const
  real(sp) :: delta_c, ignore_me, rmsrr, rmsrr0
  character(80) :: str, flnm, str1
  
  nx=iargc()
  write(*,*) '======= Joint PDF of Q and R for the velocity gradient ========'
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./jpdf-aij-qr.x nx filelist ndel '
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz))
  allocate(a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz))
  allocate(a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz))
  allocate(g(lx1,ly,lz), kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24._sp)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pqr = 0._dp
  rmsrr = 0._sp
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) ux
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uy
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uz
    close(15)

    a11 = ux * g
    a22 = uy * g
    a33 = uz * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a12(ii,jj,kk) = eye * ky(jj) * a11(ii,jj,kk)
      a13(ii,jj,kk) = eye * kz(kk) * a11(ii,jj,kk)
      a21(ii,jj,kk) = eye * kx(ii) * a22(ii,jj,kk)
      a23(ii,jj,kk) = eye * kz(kk) * a22(ii,jj,kk)
      a31(ii,jj,kk) = eye * kx(ii) * a33(ii,jj,kk) 
      a32(ii,jj,kk) = eye * ky(jj) * a33(ii,jj,kk)

      a11(ii,jj,kk) = eye * kx(ii) * a11(ii,jj,kk)
      a22(ii,jj,kk) = eye * ky(jj) * a22(ii,jj,kk)
      a33(ii,jj,kk) = eye * kz(kk) * a33(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

    if ( nfile .eq. 0 ) then
      rmsrr0 = sum( a11(1:lx,:,:) * conjg( a11(1:lx,:,:) ) ) &
             + sum( a12(1:lx,:,:) * conjg( a12(1:lx,:,:) ) ) &
             + sum( a13(1:lx,:,:) * conjg( a13(1:lx,:,:) ) ) &
             + sum( a21(1:lx,:,:) * conjg( a21(1:lx,:,:) ) ) &
             + sum( a22(1:lx,:,:) * conjg( a22(1:lx,:,:) ) ) &
             + sum( a23(1:lx,:,:) * conjg( a23(1:lx,:,:) ) ) &
             + sum( a31(1:lx,:,:) * conjg( a31(1:lx,:,:) ) ) &
             + sum( a32(1:lx,:,:) * conjg( a32(1:lx,:,:) ) ) &
             + sum( a33(1:lx,:,:) * conjg( a33(1:lx,:,:) ) ) 
      rmsrr0 = sqrt( rmsrr0 * const )         

      write(*,*) 'Estimated rms of rr: ', rmsrr0
    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii + 1 )/2
        gg(1,1) = real( a11(ll,jj,kk), kind=sp )
        gg(1,2) = real( a12(ll,jj,kk), kind=sp )
        gg(1,3) = real( a13(ll,jj,kk), kind=sp )
        gg(2,1) = real( a21(ll,jj,kk), kind=sp )
        gg(2,2) = real( a22(ll,jj,kk), kind=sp )
        gg(2,3) = real( a23(ll,jj,kk), kind=sp )
        gg(3,1) = real( a31(ll,jj,kk), kind=sp )
        gg(3,2) = real( a32(ll,jj,kk), kind=sp )
        gg(3,3) = real( a33(ll,jj,kk), kind=sp )
      else
        ll = ii / 2
        gg(1,1) = aimag( a11(ll,jj,kk) )
        gg(1,2) = aimag( a12(ll,jj,kk) )
        gg(1,3) = aimag( a13(ll,jj,kk) )
        gg(2,1) = aimag( a21(ll,jj,kk) )
        gg(2,2) = aimag( a22(ll,jj,kk) )
        gg(2,3) = aimag( a23(ll,jj,kk) )
        gg(3,1) = aimag( a31(ll,jj,kk) )
        gg(3,2) = aimag( a32(ll,jj,kk) )
        gg(3,3) = aimag( a33(ll,jj,kk) )
      end if

      rmsrr = rmsrr + sum( gg * gg )

      ignore_me = -.5_sp * sum(  gg * transpose(gg) ) / rmsrr0**2
      delta_c = -sum ( matmul(gg,gg) * transpose(gg) ) / 3._sp /rmsrr0**3

      ll = floor( (ignore_me + boundq) / binwq ) + 1
      mm = floor( (delta_c + boundr) / binwr ) + 1

      if ( ll .ge. 1 .and. ll .le. npnt .and. mm .ge. 1 .and. mm .le. npnt) then
        pqr(mm,ll) = pqr(mm,ll) + 1
      end if

    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  rmsrr = rmsrr / nfile * const
  rmsrr = sqrt(rmsrr)
  write(*,*) 'rmsrr is: ', rmsrr

  pqr = pqr / nfile * const
  write(*,*) 'Check normalizationof pqr:', sum(pqr)

  pqr = pqr / binwq / binwr

  open(15,file='jpdf-qr-aij-'//str(1:len_trim(str))//'dx-'//trim(flnm)//'.dat', form = 'binary')
    write(15) pqr
    write(15) ((-boundr + (ii-.5_sp)*binwr)*rmsrr0**3/rmsrr**3, ii = 1, npnt)
    write(15) ((-boundq + (jj-.5_sp)*binwq)*rmsrr0**2/rmsrr**2, jj = 1, npnt)
  close(15)

  deallocate(ux, uy, uz, kx, ky, kz, g)
  deallocate(a11, a12, a13, a21, a22, a23, a31, a32, a33)

  call destroyplan3d

  write(*,*) 'Finished'

end program jpdfaijqr
