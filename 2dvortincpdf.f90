program twodvortincpdf
  use mconstant
  use mwavenumber
  use mfftwplan2d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfdom
  
  integer :: nx,ny,ii,jj,lx1,lx,ly,ndel, nfile
  real(sp) :: ignore_me
  real(dp) :: rmsdom, rmsdom0

  complex(sp), allocatable, dimension(:,:) :: ux,uy
  real(sp), allocatable, dimension(:,:) :: omr
  real(sp), allocatable, dimension(:)   :: kx, ky
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dvortincpdf.x nx ndel filelist'
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
  allocate( omr(nx,ny) )
  allocate( kx(lx1), ky(ly) )

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  ! file number 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  write(*,*) 'fftwplan'
  call fftwplan2de(nx,ny)
  call wavenumber(kx,ky,lx1,ly)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfdom = 0.d0
  rmsdom = 0.d0
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

    do jj = 1, ly
    do ii = 1, lx1
        ux(ii,jj) = eye * ( kx(ii) * uy(ii,jj) - ky(jj) * ux(ii,jj) )
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r2d,ux,ignore_me)
    
    omr(1:nx:2,:) =  real( ux(1:lx,:) )
    omr(2:nx:2,:) = aimag( ux(1:lx,:) )
 
    if ( nfile .eq. 0 ) then

      rmsdom0 = 0.
      do jj = 1, ny
      do ii = 1, nx
        ignore_me = omr( modulo(ii+ndel-1, nx) + 1, jj ) - omr(ii,jj)
        rmsdom0 = rmsdom0 + ignore_me * ignore_me
      end do
      end do
      ignore_me = 1./(nx*ny)
      rmsdom0 = sqrt( rmsdom0 * ignore_me )

    end if
 
    do jj=1,ny
    do ii=1,nx
      
      ! increment along x direction
      ignore_me = omr( modulo(ii+ndel-1, nx) + 1, jj ) - omr(ii, jj)

      rmsdom = rmsdom + ignore_me * ignore_me
      ignore_me = ignore_me / rmsdom0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfdom(ly)=pdfdom(ly)+1
 
      ! increment along y direction
      ignore_me = omr( ii, modulo(jj+ndel-1,ny) + 1 ) - omr(ii, jj)
 
      rmsdom = rmsdom + ignore_me * ignore_me
      ignore_me = ignore_me / rmsdom0

      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfdom(ly)=pdfdom(ly)+1
 
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (2.*nx*ny) / nfile

  rmsdom = sqrt( rmsdom * ignore_me )
  pdfdom = pdfdom * ignore_me

  write(*,*) 'Check normalization: pdfdom', sum(pdfdom)

  pdfdom = pdfdom / binw

  path='2dvortincpdf-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  open(15,file=path)
    write(15,'(''# Title = "'',1E12.4, ''"'')') rmsdom
    write(15,'(''# variables = "dom", "pdfdom"'')') 
    do ii=1,npnt
      ignore_me = -bnd + (ii - .5) * binw
      write(15,'(10E12.4)') ignore_me * rmsdom0 / rmsdom, pdfdom(ii) * rmsdom / rmsdom0
    end do
  close(15)
  
  
  deallocate(ux,uy,omr,kx,ky)
  
  call destroyplan2d

  write(*,*) 'done.'
end program twodvortincpdf
