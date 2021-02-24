program ptn
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel, nprtcle, ip
  real(sp) :: ignore_me, delta_c, dt, time, const 
  real(sp) :: meannonalign 
  real(sp) :: nonalign, wxp, wyp, wzp

  complex(sp) :: wxc, wyc, wzc

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz
  real(4),     allocatable, dimension(:)     :: xp, yp, zp

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
  real(sp),    allocatable, dimension(:,:,:) :: g, wwr
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  real(sp) :: y(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  nx=iargc()
  if (nx .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pt-nonalign.x nx filelist ndel nprtcle dt prefix'
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
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  dxyz(1)=2.*pi/real(nx)
  dxyz(2)=2.*pi/real(ny)
  dxyz(3)=2.*pi/real(nz)

  allocate( ux (lx1,ly,lz), uy (lx1,ly,lz), uz (lx1,ly,lz) )
  allocate( wx (lx1,ly,lz), wy (lx1,ly,lz), wz (lx1,ly,lz) )
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s21(lx1,ly,lz), s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( s31(lx1,ly,lz), s32(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate( wwr(nx,ny,nz) )
  allocate( xp(nprtcle), yp(nprtcle), zp(nprtcle) )
  allocate( g(lx1,ly,lz),kx(lx1), ky(ly), kz(lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(25, file = 'nonalign-'//flnm(1:len_trim(flnm))//'.dat')
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

    wx = ux * g
    wy = uy * g
    wz = uz * g

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wxc = eye * ky(jj) * wz(ii,jj,kk) - eye * kz(kk) * wy(ii,jj,kk) 
      wyc = eye * kz(kk) * wx(ii,jj,kk) - eye * kx(ii) * wz(ii,jj,kk)
      wzc = eye * kx(ii) * wy(ii,jj,kk) - eye * ky(jj) * wx(ii,jj,kk)

      wx(ii,jj,kk) = wxc
      wy(ii,jj,kk) = wyc
      wz(ii,jj,kk) = wzc

      s11(ii,jj,kk) = eye * wx(ii,jj,kk) * kx(ii)
      s12(ii,jj,kk) = eye * wx(ii,jj,kk) * ky(jj)
      s13(ii,jj,kk) = eye * wx(ii,jj,kk) * kz(kk)
      s21(ii,jj,kk) = eye * wy(ii,jj,kk) * kx(ii)
      s22(ii,jj,kk) = eye * wy(ii,jj,kk) * ky(jj)
      s23(ii,jj,kk) = eye * wy(ii,jj,kk) * kz(kk)
      s31(ii,jj,kk) = eye * wz(ii,jj,kk) * kx(ii)
      s32(ii,jj,kk) = eye * wz(ii,jj,kk) * ky(jj)
      s33(ii,jj,kk) = eye * wz(ii,jj,kk) * kz(kk)

    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2

          delta_c = real(wx(ll,jj,kk))**2 + real(wy(ll,jj,kk))**2 + real(wz(ll,jj,kk))**2

          nonalign =  real( s11(ll,jj,kk) )**2 &
                   +  real( s12(ll,jj,kk) )**2 &
                   +  real( s13(ll,jj,kk) )**2 &
                   +  real( s21(ll,jj,kk) )**2 &
                   +  real( s22(ll,jj,kk) )**2 &
                   +  real( s23(ll,jj,kk) )**2 &
                   +  real( s31(ll,jj,kk) )**2 &
                   +  real( s32(ll,jj,kk) )**2 &
                   +  real( s33(ll,jj,kk) )**2

          nonalign = nonalign / (delta_c + mytiny)

          wxp = real(wx(ll,jj,kk)) * real(s11(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s12(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s13(ll,jj,kk))
          wyp = real(wx(ll,jj,kk)) * real(s21(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s22(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s23(ll,jj,kk))
          wzp = real(wx(ll,jj,kk)) * real(s31(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s32(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s33(ll,jj,kk))

          nonalign = nonalign - (wxp*wxp+wyp*wyp+wzp*wzp)/(delta_c**2+mytiny)
          nonalign = nonalign + ( wxp * real( wx(ll,jj,kk) ) &
                                + wyp * real( wy(ll,jj,kk) ) &
                                + wzp * real( wz(ll,jj,kk) ) ) ** 2 / (delta_c**3+mytiny)

          wxp = real(wx(ll,jj,kk)) * real(s11(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s21(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s31(ll,jj,kk))
          wyp = real(wx(ll,jj,kk)) * real(s12(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s22(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s32(ll,jj,kk))
          wzp = real(wx(ll,jj,kk)) * real(s13(ll,jj,kk)) &
              + real(wy(ll,jj,kk)) * real(s23(ll,jj,kk)) &
              + real(wz(ll,jj,kk)) * real(s33(ll,jj,kk))

          nonalign = nonalign - (wxp*wxp+wyp*wyp+wzp*wzp)/(delta_c**2+mytiny)
      else
          ll = ii / 2

          delta_c = aimag(wx(ll,jj,kk))**2 + aimag(wy(ll,jj,kk))**2 + aimag(wz(ll,jj,kk))**2

          nonalign =  aimag( s11(ll,jj,kk) )**2 &
                   +  aimag( s12(ll,jj,kk) )**2 &
                   +  aimag( s13(ll,jj,kk) )**2 &
                   +  aimag( s21(ll,jj,kk) )**2 &
                   +  aimag( s22(ll,jj,kk) )**2 &
                   +  aimag( s23(ll,jj,kk) )**2 &
                   +  aimag( s31(ll,jj,kk) )**2 &
                   +  aimag( s32(ll,jj,kk) )**2 &
                   +  aimag( s33(ll,jj,kk) )**2

          nonalign = nonalign / (delta_c + mytiny)

          wxp = aimag(wx(ll,jj,kk)) * aimag(s11(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s12(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s13(ll,jj,kk))
          wyp = aimag(wx(ll,jj,kk)) * aimag(s21(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s22(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s23(ll,jj,kk))
          wzp = aimag(wx(ll,jj,kk)) * aimag(s31(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s32(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s33(ll,jj,kk))

          nonalign = nonalign - (wxp*wxp+wyp*wyp+wzp*wzp)/(delta_c**2+mytiny)
          nonalign = nonalign + ( wxp * aimag( wx(ll,jj,kk) ) &
                                + wyp * aimag( wy(ll,jj,kk) ) &
                                + wzp * aimag( wz(ll,jj,kk) ) ) ** 2 / (delta_c**2+mytiny)

          wxp = aimag(wx(ll,jj,kk)) * aimag(s11(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s21(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s31(ll,jj,kk))
          wyp = aimag(wx(ll,jj,kk)) * aimag(s12(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s22(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s32(ll,jj,kk))
          wzp = aimag(wx(ll,jj,kk)) * aimag(s13(ll,jj,kk)) &
              + aimag(wy(ll,jj,kk)) * aimag(s23(ll,jj,kk)) &
              + aimag(wz(ll,jj,kk)) * aimag(s33(ll,jj,kk))

          nonalign = nonalign - (wxp*wxp+wyp*wyp+wzp*wzp)/(delta_c**2+mytiny)
      end if

      wwr(ii,jj,kk) = nonalign
      
    end do
    end do
    end do

    meannonalign = 0.
    do ip=1, nprtcle
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)
 
      call pre_interp(y,dxyz,bg,lhnode)
      call value(wwr,nonalign,lhnode,bg,nx,ny,nz)

      meannonalign = meannonalign + nonalign

    end do
    meannonalign =  meannonalign / nprtcle 

    write(25, '(10E15.4)') time, meannonalign
    time = time + dt

    nfile = nfile + 1
  end do
  close(25)
  close(30)

  deallocate(wx, wy, wz, kx, ky, kz, g, xp, yp, zp)
  deallocate(s11,s12,s13,s21,s22,s23,s31,s32,s33,wwr)
  deallocate(ux, uy, uz)

  call destroyplan3d

  write(*,*) 'Finished'

end program ptn      
