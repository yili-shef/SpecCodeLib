program spcndss
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer  :: ndel,nx,ny,nz,lx,lx1,ly,lz, nfile

  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33
  real(sp),    allocatable, dimension(:,:,:) :: kx,ky,kz,g, csijpij, ss

  integer,  parameter :: npnt = 60
  real(sp), parameter :: bnd = 45._sp, binw = bnd/npnt
  real(sp), dimension(npnt) :: pdfss, sijpijcndss

  integer :: ii,jj,kk,ll

  character(80) :: fnm,str,fpath,strdel
  real(sp) :: delta_c, ignore_me

  write(*,*) 
  write(*,'(''>>>>>> Correlation between aij and pij <<<<<<'')')
  write(*,*)
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./spcndss.x nx fnamelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        fnamelist: name of list file'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter parameter
  call getarg(3,strdel)
  read(strdel, '(I20)') ndel
  strdel = adjustl(strdel)

  ! file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p12(lx1,ly,lz), p13(lx1,ly,lz) )
  allocate( p22(lx1,ly,lz), p23(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate(   g(lx1,ly,lz), csijpij(nx,ny,nz), ss(nx,ny,nz) )                                  
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g = exp(-g*delta_c**2/24.)

  open(33, file = fnm(1:len_trim(fnm))//'.list')

  pdfss = 0._sp
  sijpijcndss = 0._sp
  nfile = 0
  do while ( .not. eof(33) )

    read(33,*) str
    str = adjustl(str)
    write(*,*) 'reading data', str

    fpath='./out/ux'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p11
    close(10)
    fpath='./out/uy'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p22
    close(10)
    fpath='./out/uz'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p33
    close(10)
    fpath='./out/p'//str(1:len_trim(str))
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)p12
    close(10)
    write(*,*) 'after reading data files'
 
    p11 = p11 * g
    p22 = p22 * g 
    p33 = p33 * g
    p12 = p12 * g
 
    s11 = eye * kx * p11
    s22 = eye * ky * p22
    s33 = eye * kz * p33
    s12 = eye * ( kx * p22 + ky * p11 ) * 0.5_sp
    s13 = eye * ( kx * p33 + kz * p11 ) * 0.5_sp
    s23 = eye * ( ky * p33 + kz * p22 ) * 0.5_sp
 
    p11 = - kx * kx * p12
    p22 = - ky * ky * p12
    p33 = - kz * kz * p12
    p13 = - kx * kz * p12
    p23 = - ky * kz * p12
    p12 = - kx * ky * p12
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
 
    csijpij(1:nx:2,:,:) = real(s11(1:lx,:,:)) * real(p11(1:lx,:,:)) &
                        + real(s22(1:lx,:,:)) * real(p22(1:lx,:,:)) &
                        + real(s33(1:lx,:,:)) * real(p33(1:lx,:,:)) &
                        + 2 * ( real(s12(1:lx,:,:)) * real(p12(1:lx,:,:)) &
                        +       real(s13(1:lx,:,:)) * real(p13(1:lx,:,:)) &
                        +       real(s23(1:lx,:,:)) * real(p23(1:lx,:,:)) ) 

    csijpij(2:nx:2,:,:) = aimag(s11(1:lx,:,:)) * aimag(p11(1:lx,:,:)) &
                        + aimag(s22(1:lx,:,:)) * aimag(p22(1:lx,:,:)) &
                        + aimag(s33(1:lx,:,:)) * aimag(p33(1:lx,:,:)) &
                        + 2 * ( aimag(s12(1:lx,:,:)) * aimag(p12(1:lx,:,:)) &
                        +       aimag(s13(1:lx,:,:)) * aimag(p13(1:lx,:,:)) &
                        +       aimag(s23(1:lx,:,:)) * aimag(p23(1:lx,:,:)) ) 
 

    ss(1:nx:2,:,:) = real(s11(1:lx,:,:)) * real(s11(1:lx,:,:)) &
                   + real(s22(1:lx,:,:)) * real(s22(1:lx,:,:)) &
                   + real(s33(1:lx,:,:)) * real(s33(1:lx,:,:)) &
                   + 2 * ( real(s12(1:lx,:,:)) * real(s12(1:lx,:,:)) &
                   +       real(s13(1:lx,:,:)) * real(s13(1:lx,:,:)) &
                   +       real(s23(1:lx,:,:)) * real(s23(1:lx,:,:)) ) 

    ss(2:nx:2,:,:) = aimag(s11(1:lx,:,:)) * aimag(s11(1:lx,:,:)) &
                   + aimag(s22(1:lx,:,:)) * aimag(s22(1:lx,:,:)) &
                   + aimag(s33(1:lx,:,:)) * aimag(s33(1:lx,:,:)) &
                   + 2 * ( aimag(s12(1:lx,:,:)) * aimag(s12(1:lx,:,:)) &
                   +       aimag(s13(1:lx,:,:)) * aimag(s13(1:lx,:,:)) &
                   +       aimag(s23(1:lx,:,:)) * aimag(s23(1:lx,:,:)) ) 
 
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      ll = floor( ss(ii,jj,kk) / binw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) then
          pdfss(ll) = pdfss(ll) + 1
          sijpijcndss(ll) = sijpijcndss(ll) + csijpij(ii,jj,kk)
      end if

    end do
    end do
    end do
 

    nfile = nfile + 1
  end do
  close(33)

  ignore_me = 1._sp / (nx*ny*nz*nfile)

  sijpijcndss = sijpijcndss / (pdfss + mytiny)
  pdfss = pdfss * ignore_me

  write(*,*) 'check normalization of pdfss: ', sum(pdfss)
  pdfss = pdfss / binw

  open(25, file = 'spcndss-'//fnm(1:len_trim(fnm))//'-'//strdel(1:len_trim(strdel))//'dx.dat')
  do ll = 1, npnt
    write(25, '(15E15.3)') (ll-.5_sp)*binw, sijpijcndss(ll), pdfss(ll)
  end do
  close(25)

  deallocate(kx,ky,kz)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(s11,s12,s13,s22,s23,s33)
  deallocate(g, csijpij, ss)

  call destroyplan3d

  write(*,*) 'done'

end program spcndss
