program scalarincpdf
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfx,pdfy,pdfz,pdfall
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ndel, nfile
  real(sp) :: ignore_me
  real(dp) :: rmsx, rmsy, rmsz, rmsx0, rmsall

  complex(sp), allocatable, dimension(:,:,:) :: phi
  real(sp), allocatable, dimension(:,:,:) :: phir, dphir
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./scalarincpdf.x nx ndel filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     ndel: length of displacement (l=ndel*dx)'
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

  allocate( phi(lx1,ly,lz) )
  allocate( phir(nx,ny,nz), dphir(nx,ny,nz) )

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  ! file list 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  write(*,*) 'fftwplan'
  call fftwplan3de(nx,ny,nz)


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfx = 0.d0; pdfy = 0.d0; pdfz = 0.d0
  rmsx = 0.d0; rmsy = 0.d0; rmsz = 0.d0
  pdfall = 0.d0; rmsall = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/phi'//str1(1:len_trim(str1)),form='unformatted')
      read(15) phi
    close(15)
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,phi,ignore_me)
    
    phir(1:nx:2,:,:) =  real( phi(1:lx,:,:) )
    phir(2:nx:2,:,:) = aimag( phi(1:lx,:,:) )
 
    
    dphir = cshift(phir, ndel, 1) - phir
 
    if ( nfile .eq. 0 ) rmsx0 = sqrt( sum(dphir * dphir) / (nx*ny*nz) )
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
      ignore_me = dphir(ii,jj,kk)

      rmsx = rmsx + ignore_me * ignore_me
      rmsall = rmsall + ignore_me * ignore_me

      ignore_me = ignore_me / rmsx0

      ly = 1 + floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) then
          pdfx(ly) = pdfx(ly) + 1
          pdfall(ly) = pdfall(ly) + 1
      end if
 
    end do
    end do
    end do

    dphir = cshift(phir, ndel, 2) - phir
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
      ignore_me = dphir(ii,jj,kk)

      rmsy = rmsy + ignore_me * ignore_me
      rmsall = rmsall + ignore_me * ignore_me

      ignore_me = ignore_me / rmsx0

      ly = 1 + floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) then
          pdfy(ly) = pdfy(ly) + 1
          pdfall(ly) = pdfall(ly) + 1
      end if
 
    end do
    end do
    end do

    dphir = cshift(phir, ndel, 3) - phir
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
      ignore_me = dphir(ii,jj,kk)

      rmsz = rmsz + ignore_me * ignore_me
      rmsall = rmsall + ignore_me * ignore_me

      ignore_me = ignore_me / rmsx0

      ly = 1 + floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) then
          pdfz(ly) = pdfz(ly) + 1
          pdfall(ly) = pdfall(ly) + 1
      end if
 
    end do
    end do
    end do


    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  rmsz = sqrt( rmsz * ignore_me )
  pdfx = pdfx / binw * ignore_me
  pdfy = pdfy / binw * ignore_me
  pdfz = pdfz / binw * ignore_me

  rmsall = sqrt( rmsall * ignore_me / 3 )
  pdfall = pdfall / binw * ignore_me / 3
 
  
  path='scalarpdf-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)

    write(15,'(''# Title = "'',4E12.4, ''"'')') rmsx, rmsy, rmsz, rmsall
    write(15,'(''# variables = "xx", "pdfdx","yy","pdfdy","zz", "pdfdz", "all", "pdfall"'')')

    do lz=1,npnt

      ignore_me = -bnd + (lz - .5) * binw

      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, pdfx(lz) * rmsx / rmsx0, &
      ignore_me * rmsx0 / rmsy, rmsy / rmsx0 * pdfy(lz), ignore_me * rmsx0 / rmsz, &
      rmsz / rmsx0 * pdfz(lz), ignore_me * rmsx0 / rmsall, &
      rmsall / rmsx0 * pdfall(lz)

    end do

  close(15)
  
  
  deallocate(phi, phir, dphir)
  
  call destroyplan3d

  write(*,*) 'scalarincpdf.x done.'
end program scalarincpdf
