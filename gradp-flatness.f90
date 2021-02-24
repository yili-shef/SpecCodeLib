program gradp
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile, ndel
  real(sp) :: ignore_me, const

  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, g

  character(80) :: str, flnm, prf, str1
  real(sp) :: afntmp, dpdx, dpdy, dpdz, mp4, mp2, delta_c
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./gradp-flatness.x nx filelist ndel'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: ndel * dx = filter scale'
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

  ! ndel
  call getarg(3,prf)
  read(prf, '(I20)') ndel
  prf = adjustl(prf)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate(   g(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  g = exp(-g*delta_c**2/24.)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  mp4 = 0._sp
  mp2 = 0._sp
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p33
    close(15)
    write(*,*) 'finishing reading ux uy uz and p'

    p33 = p33 * g

    ! dp/dx, dp/dy, dp/dz
    p11 = p33 * eye * kx ; p22 = p33 * eye * ky ; p33 = p33 * eye * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          dpdx = real( p11(ll,jj,kk),sp )
          dpdy = real( p22(ll,jj,kk),sp )
          dpdz = real( p33(ll,jj,kk),sp )

      else
          ll = ii/2
          dpdx = aimag( p11(ll,jj,kk) )
          dpdy = aimag( p22(ll,jj,kk) )
          dpdz = aimag( p33(ll,jj,kk) )

      endif
      afntmp = dpdx * dpdx + dpdy * dpdy + dpdz * dpdz
      mp4 = mp4 + afntmp * afntmp
      mp2 = mp2 + afntmp

    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  mp4 = mp4 * const / nfile
  mp2 = mp2 * const / nfile

  open(28, file = 'flatness-gradp-'//trim(prf)//'dx-'//trim(flnm)//'.dat')
    write(28, '(15E15.3)') mp4/mp2*2
  close(28)

  deallocate(kx,ky,kz,p11,p22,p33, g)
  
  call destroyplan3d

  write(*,*) 'gradp-flatness.x done.'

end program gradp
