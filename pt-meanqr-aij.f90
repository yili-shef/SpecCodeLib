program ptqraij
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, const, dt, time 
  real(sp) :: meanqaij, meanraij
  real(sp) :: qaijp, raijp

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
  real(sp),    allocatable, dimension(:,:,:) :: qaij, raij, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-qr-aij.x nx filelist ndel nprtcle dt prefix'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*) '                     nprtcle: number of particles'
          write(*,*) '                     dt: time step size'
          write(*,*) '                     prefix: prefix of coordinates data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! number of particles
  call getarg(4,flnm)
  read(flnm, '(I20)') nprtcle

  ! time step size
  call getarg(5,flnm)
  read(flnm, '(F20.10)') dt

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! prefix
  call getarg(6,prefix)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(xp(nprtcle), yp(nprtcle), zp(nprtcle))

  allocate(a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz))
  allocate(a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz))
  allocate(a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz))
  allocate(qaij(nx,ny,nz), raij(nx,ny,nz))

  allocate(g(lx1,ly,lz), kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)


  open(25, file = 'pt-meanqr-aij-'//flnm(1:len_trim(flnm))//'.dat')
  open(30, file = flnm(1:len_trim(flnm))//'.list')

  time = 0.
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

    open(15,file='./out/'//prefix(1:len_trim(prefix))//str1(1:len_trim(str1)),form='unformatted')
      read(15) xp, yp, zp
    close(15)
    write(*,*) 'finishing reading data'

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

    qaij(1:nx:2,:,:) = real( a11(1:lx,:,:) ) * real( a11(1:lx,:,:) ) &
                     + real( a12(1:lx,:,:) ) * real( a21(1:lx,:,:) ) &
                     + real( a13(1:lx,:,:) ) * real( a31(1:lx,:,:) ) &
                     + real( a21(1:lx,:,:) ) * real( a12(1:lx,:,:) ) &
                     + real( a22(1:lx,:,:) ) * real( a22(1:lx,:,:) ) &
                     + real( a23(1:lx,:,:) ) * real( a32(1:lx,:,:) ) &
                     + real( a31(1:lx,:,:) ) * real( a13(1:lx,:,:) ) &
                     + real( a32(1:lx,:,:) ) * real( a23(1:lx,:,:) ) &
                     + real( a33(1:lx,:,:) ) * real( a33(1:lx,:,:) ) 

    qaij(2:nx:2,:,:) = aimag( a11(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) &
                     + aimag( a12(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) &
                     + aimag( a13(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) &
                     + aimag( a21(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) &
                     + aimag( a22(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) &
                     + aimag( a23(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) &
                     + aimag( a31(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) &
                     + aimag( a32(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) &
                     + aimag( a33(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) 

    qaij = -.5 * qaij

    raij(1:nx:2,:,:) = real( a11(1:lx,:,:) ) * real( a11(1:lx,:,:) ) * real( a11(1:lx,:,:) ) &
                     + real( a12(1:lx,:,:) ) * real( a21(1:lx,:,:) ) * real( a11(1:lx,:,:) ) &
                     + real( a13(1:lx,:,:) ) * real( a31(1:lx,:,:) ) * real( a11(1:lx,:,:) ) &

                     + real( a11(1:lx,:,:) ) * real( a12(1:lx,:,:) ) * real( a21(1:lx,:,:) ) &
                     + real( a12(1:lx,:,:) ) * real( a22(1:lx,:,:) ) * real( a21(1:lx,:,:) ) &
                     + real( a13(1:lx,:,:) ) * real( a32(1:lx,:,:) ) * real( a21(1:lx,:,:) ) &

                     + real( a11(1:lx,:,:) ) * real( a13(1:lx,:,:) ) * real( a31(1:lx,:,:) ) &
                     + real( a12(1:lx,:,:) ) * real( a23(1:lx,:,:) ) * real( a31(1:lx,:,:) ) &
                     + real( a13(1:lx,:,:) ) * real( a33(1:lx,:,:) ) * real( a31(1:lx,:,:) ) &

                     + real( a21(1:lx,:,:) ) * real( a11(1:lx,:,:) ) * real( a12(1:lx,:,:) ) &
                     + real( a22(1:lx,:,:) ) * real( a21(1:lx,:,:) ) * real( a12(1:lx,:,:) ) &
                     + real( a23(1:lx,:,:) ) * real( a31(1:lx,:,:) ) * real( a12(1:lx,:,:) ) &

                     + real( a21(1:lx,:,:) ) * real( a12(1:lx,:,:) ) * real( a22(1:lx,:,:) ) &
                     + real( a22(1:lx,:,:) ) * real( a22(1:lx,:,:) ) * real( a22(1:lx,:,:) ) &
                     + real( a23(1:lx,:,:) ) * real( a32(1:lx,:,:) ) * real( a22(1:lx,:,:) ) &

                     + real( a21(1:lx,:,:) ) * real( a13(1:lx,:,:) ) * real( a32(1:lx,:,:) ) &
                     + real( a22(1:lx,:,:) ) * real( a23(1:lx,:,:) ) * real( a32(1:lx,:,:) ) &
                     + real( a23(1:lx,:,:) ) * real( a33(1:lx,:,:) ) * real( a32(1:lx,:,:) ) &

                     + real( a31(1:lx,:,:) ) * real( a11(1:lx,:,:) ) * real( a13(1:lx,:,:) ) &
                     + real( a32(1:lx,:,:) ) * real( a21(1:lx,:,:) ) * real( a13(1:lx,:,:) ) &
                     + real( a33(1:lx,:,:) ) * real( a31(1:lx,:,:) ) * real( a13(1:lx,:,:) ) &

                     + real( a31(1:lx,:,:) ) * real( a12(1:lx,:,:) ) * real( a23(1:lx,:,:) ) &
                     + real( a32(1:lx,:,:) ) * real( a22(1:lx,:,:) ) * real( a23(1:lx,:,:) ) &
                     + real( a33(1:lx,:,:) ) * real( a32(1:lx,:,:) ) * real( a23(1:lx,:,:) ) &

                     + real( a31(1:lx,:,:) ) * real( a13(1:lx,:,:) ) * real( a33(1:lx,:,:) ) &
                     + real( a32(1:lx,:,:) ) * real( a23(1:lx,:,:) ) * real( a33(1:lx,:,:) ) &
                     + real( a33(1:lx,:,:) ) * real( a33(1:lx,:,:) ) * real( a33(1:lx,:,:) ) 

    raij(2:nx:2,:,:) = aimag( a11(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) &
                     + aimag( a12(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) &
                     + aimag( a13(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) &

                     + aimag( a11(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) &
                     + aimag( a12(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) &
                     + aimag( a13(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) &

                     + aimag( a11(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) &
                     + aimag( a12(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) &
                     + aimag( a13(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) &

                     + aimag( a21(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) &
                     + aimag( a22(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) &
                     + aimag( a23(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) &

                     + aimag( a21(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) &
                     + aimag( a22(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) &
                     + aimag( a23(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) &

                     + aimag( a21(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) &
                     + aimag( a22(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) &
                     + aimag( a23(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) &

                     + aimag( a31(1:lx,:,:) ) * aimag( a11(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) &
                     + aimag( a32(1:lx,:,:) ) * aimag( a21(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) &
                     + aimag( a33(1:lx,:,:) ) * aimag( a31(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) &

                     + aimag( a31(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) &
                     + aimag( a32(1:lx,:,:) ) * aimag( a22(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) &
                     + aimag( a33(1:lx,:,:) ) * aimag( a32(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) &

                     + aimag( a31(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) &
                     + aimag( a32(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) &
                     + aimag( a33(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) * aimag( a33(1:lx,:,:) ) 

    raij = - raij / 3.


    meanqaij = 0.
    meanraij = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(qaij,qaijp,lhnode,bg,nx,ny,nz)
      call value(raij,raijp,lhnode,bg,nx,ny,nz)

      meanqaij = meanqaij + qaijp
      meanraij = meanraij + raijp

    end do
    meanqaij = meanqaij / nprtcle
    meanraij = meanraij / nprtcle

    write(25, '(10E15.4)') time, meanqaij, meanraij
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(ux, uy, uz, kx, ky, kz, g, xp, yp, zp)
  deallocate(a11, a12, a13, a21, a22, a23, a31, a32, a33, qaij, raij)

  call destroyplan3d

  write(*,*) 'Finished'


end program ptqraij      
