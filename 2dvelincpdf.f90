program twodvelincpdf
  use mconstant
  use mfftwplan2d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfx,pdfy
  
  integer :: nx,ny,ii,jj,lx1,lx,ly,ndel, nfile
  real(sp) :: ignore_me
  real(dp) :: rmsx, rmsy, rmsx0, rmsy0

  complex(sp), allocatable, dimension(:,:) :: ux,uy
  real(sp), allocatable, dimension(:,:) :: uxr, uyr
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dvelincpdf.x nx ndel filelist'
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

  ny=nx
  lx=nx/2;lx1=nx/2+1;ly=ny

  allocate( ux(lx1,ly), uy(lx1,ly) )
  allocate( uxr(nx,ny), uyr(nx,ny) )

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  ! file number 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  write(*,*) 'fftwplan'
  call fftwplan2de(nx,ny)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfx = 0.d0; pdfy = 0.d0
  rmsx = 0.d0; rmsy = 0.d0
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
    
    call rfftwnd_f77_one_complex_to_real(c2r2d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,uy,ignore_me)
    
    uxr(1:nx:2,:) =  real( ux(1:lx,:) )
    uxr(2:nx:2,:) = aimag( ux(1:lx,:) )
    uyr(1:nx:2,:) =  real( uy(1:lx,:) )
    uyr(2:nx:2,:) = aimag( uy(1:lx,:) )
 
    if ( nfile .eq. 0 ) then

      rmsx0 = 0.
      rmsy0 = 0.
      do jj = 1, ny
      do ii = 1, nx
        ignore_me = uxr( modulo(ii+ndel-1, nx) + 1, jj ) - uxr(ii,jj)
        rmsx0 = rmsx0 + ignore_me * ignore_me
        ignore_me = uyr( modulo(ii+ndel-1, nx) + 1, jj ) - uyr(ii,jj)
        rmsy0 = rmsy0 + ignore_me * ignore_me
      end do
      end do
      ignore_me = 1./(nx*ny)
      rmsx0 = sqrt( rmsx0 * ignore_me )
      rmsy0 = sqrt( rmsy0 * ignore_me )

    end if
 
    do jj=1,ny
    do ii=1,nx
      
      ! increment along x direction
      ignore_me = uxr( modulo(ii+ndel-1, nx) + 1, jj ) - uxr(ii, jj)

      rmsx = rmsx + ignore_me * ignore_me
      ignore_me = ignore_me / rmsx0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfx(ly)=pdfx(ly)+1
 
      ignore_me = uyr( modulo(ii+ndel-1, nx) + 1, jj ) - uyr(ii, jj)

      rmsy = rmsy + ignore_me * ignore_me
      ignore_me = ignore_me / rmsy0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfy(ly)=pdfy(ly)+1

      ! increment along y direction
      ignore_me = uyr( ii, modulo(jj+ndel-1,ny) + 1 ) - uyr(ii, jj)
 
      rmsx = rmsx + ignore_me * ignore_me
      ignore_me = ignore_me / rmsx0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfx(ly)=pdfx(ly)+1
 
      ignore_me = -( uxr( ii, modulo(jj+ndel-1,ny) + 1 ) - uxr(ii, jj) )

      rmsy = rmsy + ignore_me * ignore_me
      ignore_me = ignore_me / rmsy0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfy(ly)=pdfy(ly)+1

    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (2.*nx*ny) / nfile

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  pdfx=pdfx*ignore_me
  pdfy=pdfy*ignore_me

  write(*,*) 'Check normalization: pdfx', sum(pdfx)
  write(*,*) 'Check normalization: pdfy', sum(pdfy)

  pdfx=pdfx/binw
  pdfy=pdfy/binw 

  path='2dvelincpdf-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  open(15,file=path)
    write(15,'(''# Title = "'',2E12.4, ''"'')') rmsx, rmsy
    write(15,'(''# variables = "dux", "pdfdux","duy","pdfduy"'')') 
    do ii=1,npnt
      ignore_me = -bnd + (ii - .5) * binw
      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, pdfx(ii) * rmsx / rmsx0, &
      ignore_me * rmsy0 / rmsy, rmsy / rmsy0 * pdfy(ii)
    end do
  close(15)
  
  
  deallocate(ux,uy,uxr,uyr)
  
  call destroyplan2d

  write(*,*) 'done.'
end program twodvelincpdf
