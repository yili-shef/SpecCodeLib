program scalargradpdf
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10._sp, binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfx,pdfy,pdfz,pdfall
  
  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,ndel,nfile
  real(sp) :: ignore_me, const, delta_c, gx, gy, gz
  real(dp) :: rmsx, rmsy, rmsz, rmsx0, rmsall

  complex(sp), allocatable, dimension(:,:,:) :: phi, gphix, gphiy, gphiz
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./scalargradpdf.x nx ndel filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     ndel: filter scale = ndel * dx'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  const = 1./(nx*ny*nz)

  allocate(   phi(lx1,ly,lz),     g(lx1,ly,lz) )
  allocate( gphix(lx1,ly,lz), gphiy(lx1,ly,lz), gphiz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  delta_c=ndel*2*pi/nx

  ! file list 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  g=exp(-g*delta_c**2/24.)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfx = 0.d0; pdfy = 0.d0; pdfz = 0.d0
  rmsx = 0.d0; rmsy = 0.d0; rmsz = 0.d0
  pdfall = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/phi'//str1(1:len_trim(str1)),form='unformatted')
      read(15) phi
    close(15)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
        gphix(ii,jj,kk) = eye * kx(ii) * phi(ii,jj,kk)
        gphiy(ii,jj,kk) = eye * ky(jj) * phi(ii,jj,kk)
        gphiz(ii,jj,kk) = eye * kz(kk) * phi(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphix,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphiy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphiz,ignore_me)
    
    
    if ( nfile .eq. 0 ) rmsx0 = sqrt( sum( gphix(1:lx,:,:) * conjg(gphix(1:lx,:,:))) * const )
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
        if ( mod(ii,2) .eq. 1) then
            ll = (ii + 1)/2
            gx = real(gphix(ll,jj,kk))
            gy = real(gphiy(ll,jj,kk))
            gz = real(gphiz(ll,jj,kk))
        else
            ll = ii / 2
            gx = aimag(gphix(ll,jj,kk))
            gy = aimag(gphiy(ll,jj,kk))
            gz = aimag(gphiz(ll,jj,kk))
        end if

        rmsx = rmsx + gx * gx
        rmsy = rmsy + gy * gy
        rmsz = rmsz + gz * gz
    
        ignore_me = gx / rmsx0

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfx(ll) = pdfx(ll) + 1
          pdfall(ll) = pdfall(ll) + 1
        end if
 
        ignore_me = gy / rmsx0

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfy(ll) = pdfy(ll) + 1
          pdfall(ll) = pdfall(ll) + 1
        end if
 
        ignore_me = gz / rmsx0

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfz(ll) = pdfz(ll) + 1
          pdfall(ll) = pdfall(ll) + 1
        end if
 
    end do
    end do
    end do


    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsall = rmsx + rmsy + rmsz

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  rmsz = sqrt( rmsz * ignore_me )

  pdfx = pdfx / binw * ignore_me
  pdfy = pdfy / binw * ignore_me
  pdfz = pdfz / binw * ignore_me

  rmsall = sqrt( rmsall * ignore_me / 3 )
  pdfall = pdfall / binw * ignore_me / 3
 
  
  path='scalargradpdf-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)

    write(15,'(''# Title = "'',4E12.4, ''"'')') rmsx, rmsy, rmsz, rmsall
    write(15,'(''# variables = "xx", "pdfgx","yy","pdfgy","zz", "pdfgz", "all", "pdfall"'')')

    do lz=1,npnt

      ignore_me = -bnd + (lz - .5) * binw

      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, pdfx(lz) * rmsx / rmsx0, &
      ignore_me * rmsx0 / rmsy, rmsy / rmsx0 * pdfy(lz), ignore_me * rmsx0 / rmsz, &
      rmsz / rmsx0 * pdfz(lz), ignore_me * rmsx0 / rmsall, &
      rmsall / rmsx0 * pdfall(lz)

    end do

  close(15)
  
  
  deallocate(phi, gphix, gphiy, gphiz, g, kx, ky, kz)
  
  call destroyplan3d

  write(*,*) 'scalargradpdf.x done.'

end program scalargradpdf
