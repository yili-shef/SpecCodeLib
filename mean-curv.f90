program meancurvbalance
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,nfile,ndel
  real(sp) :: ignore_me, delta_c, meancurv
  real(sp) :: ww, wxp, wyp, wzp, meanalpha, const

  complex(sp) :: a12, a13, a21, a23, a31, a32

  complex(4),  allocatable, dimension(:,:,:) :: ux, uy, uz

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: wxr, wyr, wzr, alpha, sijr, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  character(80) :: str, flnm, prefix, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./mean-curv.x nx filelist ndel'
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

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
  const = 1. / (nx*ny*nz)
  delta_c=ndel*2*pi/nx

  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz))
  allocate(wxr(nx,ny,nz), wyr(nx,ny,nz), wzr(nx,ny,nz))
  allocate(s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz))
  allocate(s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz))
  allocate(alpha(nx,ny,nz), sijr(nx,ny,nz) )
  allocate(g (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  ww = 0.
  meanalpha = 0.
  meancurv = 0.
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

    wx = ux * g
    wy = uy * g
    wz = uz * g

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a12 = eye * ky(jj) * wx(ii,jj,kk)
      a13 = eye * kz(kk) * wx(ii,jj,kk)
      a21 = eye * kx(ii) * wy(ii,jj,kk)
      a23 = eye * kz(kk) * wy(ii,jj,kk)
      a31 = eye * kx(ii) * wz(ii,jj,kk) 
      a32 = eye * ky(jj) * wz(ii,jj,kk)

      s11(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
      s22(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
      s33(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)

      s12(ii,jj,kk) = .5 * (a12 + a21)
      s13(ii,jj,kk) = .5 * (a13 + a31)
      s23(ii,jj,kk) = .5 * (a23 + a32)

      wx(ii,jj,kk) = a32 - a23
      wy(ii,jj,kk) = a13 - a31
      wz(ii,jj,kk) = a21 - a12

    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    wxr(1:nx:2,:,:)=real(wx(1:lx,:,:)); wxr(2:nx:2,:,:)=aimag(wx(1:lx,:,:))
    wyr(1:nx:2,:,:)=real(wy(1:lx,:,:)); wyr(2:nx:2,:,:)=aimag(wy(1:lx,:,:))
    wzr(1:nx:2,:,:)=real(wz(1:lx,:,:)); wzr(2:nx:2,:,:)=aimag(wz(1:lx,:,:))
 
    sijr(1:nx:2,:,:)=real(s11(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s11(1:lx,:,:))
    alpha = wxr*sijr*wxr 

    sijr(1:nx:2,:,:)=real(s22(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s22(1:lx,:,:))
    alpha = alpha + wyr*sijr*wyr 

    sijr(1:nx:2,:,:)=real(s33(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s33(1:lx,:,:))
    alpha = alpha + wzr*sijr*wzr 

    sijr(1:nx:2,:,:)=real(s12(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s12(1:lx,:,:))
    alpha = alpha + 2 * wxr*sijr*wyr 

    sijr(1:nx:2,:,:)=real(s13(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s13(1:lx,:,:))
    alpha = alpha + 2 * wxr*sijr*wzr 

    sijr(1:nx:2,:,:)=real(s23(1:lx,:,:)); sijr(2:nx:2,:,:)=aimag(s23(1:lx,:,:))
    alpha = alpha + 2 * wyr*sijr*wzr

    alpha = alpha/ (wxr*wxr + wyr*wyr + wzr*wzr)

    ! Accumulate the means
    ww = ww + sum(wxr * wxr + wyr * wyr + wzr * wzr) 
    meanalpha = meanalpha + sum(alpha)

    wx = wx * const; wy = wy * const; wz = wz * const
    call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      s11(ii,jj,kk) = eye * wx(ii,jj,kk) * ky(jj)
      s12(ii,jj,kk) = eye * wx(ii,jj,kk) * kz(kk)
      s13(ii,jj,kk) = eye * wy(ii,jj,kk) * kx(ii)
      s22(ii,jj,kk) = eye * wy(ii,jj,kk) * kz(kk)
      s23(ii,jj,kk) = eye * wz(ii,jj,kk) * kx(ii)
      s33(ii,jj,kk) = eye * wz(ii,jj,kk) * ky(jj)
      wx (ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
      wy (ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
      wz (ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2
          wxp = wxr(ii,jj,kk) * real( wx (ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( s11(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( s12(ll,jj,kk) )
          wyp = wxr(ii,jj,kk) * real( s13(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( wy (ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( s22(ll,jj,kk) )
          wzp = wxr(ii,jj,kk) * real( s23(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * real( s33(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * real( wz (ll,jj,kk) )
      else
          ll = ii / 2
          wxp = wxr(ii,jj,kk) * aimag( wx (ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( s11(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( s12(ll,jj,kk) )
          wyp = wxr(ii,jj,kk) * aimag( s13(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( wy (ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( s22(ll,jj,kk) )
          wzp = wxr(ii,jj,kk) * aimag( s23(ll,jj,kk) ) &
              + wyr(ii,jj,kk) * aimag( s33(ll,jj,kk) ) &
              + wzr(ii,jj,kk) * aimag( wz (ll,jj,kk) )
      end if
      delta_c = wxr(ii,jj,kk) * wxp &
              + wyr(ii,jj,kk) * wyp &
              + wzr(ii,jj,kk) * wzp
      ignore_me =  wxr(ii,jj,kk) ** 2 &
                +  wyr(ii,jj,kk) ** 2 &
                +  wzr(ii,jj,kk) ** 2

      wxp = wxp - wxr(ll,jj,kk) * delta_c / (ignore_me + mytiny)
      wyp = wyp - wyr(ll,jj,kk) * delta_c / (ignore_me + mytiny)
      wzp = wzp - wzr(ll,jj,kk) * delta_c / (ignore_me + mytiny)

      ! curvature
      alpha(ii,jj,kk) = sqrt( wxp*wxp + wyp*wyp + wzp*wzp ) / (ignore_me + mytiny)
      
    end do
    end do
    end do

    ! accumulate the mean of curvature
    meancurv = meancurv + sum(alpha)

    nfile = nfile + 1
  end do
  close(30)

  ww = ww / nfile * const
  meanalpha = meanalpha / nfile * const
  meancurv = meancurv / nfile * const

  open(25, file = 'mean-curv-'//flnm(1:len_trim(flnm))//'.dat')
    write(25, *) '# ww <alpha> <c>'
    write(25, '(10E15.4)') ww, meanalpha, meancurv 
  close(25)

  deallocate(ux, uy, uz, wx, wy, wz, wxr, wyr, wzr, kx, ky, kz, g)
  deallocate(s11,s12,s13,s22,s23,s33,sijr,alpha)

  call destroyplan3d

  write(*,*) 'Finished'


end program meancurvbalance      
