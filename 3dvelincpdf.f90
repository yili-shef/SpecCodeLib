program vincpdf
  use mconstant
  use mfftw3
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10._sp, binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfx,pdfy,pdfz
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ndel, nfile
  real(sp) :: ignore_me
  real(dp) :: rmsx, rmsy, rmsz, rmsx0, rmsy0, rmsz0

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension(:,:,:) :: uxr, uyr, uzr
  character(80) :: str, flnm, path, str1
  
  integer(8) :: uxc2r, uyc2r, uzc2r
  
  nx=iargc()
  if (nx .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./3dvelincpdf.x nx ndel filelist idir'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     ndel: length of displacement (l=ndel*dx)'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     idir: direction of displacement(1=x,2=y,3=z,0=quit)'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(uxr(nx,ny,nz), uyr(nx,ny,nz), uzr(nx,ny,nz) )

  ! direction of displacement
  call getarg(4,str)
  read(str,'(I20)') lx1

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  ! file number 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  !write(*,*) 'fftwplan'
  !call fftwplan3de(nx,ny,nz)

  ! fftw3
  call dfftwplan3dc2r(ux, nx, uxc2r)
  call dfftwplan3dc2r(uy, nx, uyc2r)
  call dfftwplan3dc2r(uz, nx, uzc2r)
  write(*,*) 'fftwplan'


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfx = 0._dp; pdfy = 0._dp; pdfz = 0._dp
  rmsx = 0._dp; rmsy = 0._dp; rmsz = 0._dp
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
    
    ! fftw2
    !call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    !call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    !call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)

    ! fftw3
    call dfftw_execute(uxc2r)
    call dfftw_execute(uyc2r)
    call dfftw_execute(uzc2r)
    
    uxr(1:nx:2,:,:) =  real( ux(1:lx,:,:) )
    uxr(2:nx:2,:,:) = aimag( ux(1:lx,:,:) )
    uyr(1:nx:2,:,:) =  real( uy(1:lx,:,:) )
    uyr(2:nx:2,:,:) = aimag( uy(1:lx,:,:) )
    uzr(1:nx:2,:,:) =  real( uz(1:lx,:,:) )
    uzr(2:nx:2,:,:) = aimag( uz(1:lx,:,:) )
 
 
    
    uxr = cshift(uxr, ndel, lx1) - uxr
    uyr = cshift(uyr, ndel, lx1) - uyr
    uzr = cshift(uzr, ndel, lx1) - uzr
 
    if ( nfile .eq. 0 ) then
      rmsx0 = sqrt( sum(uxr * uxr) / (nx*ny*nz) )
      rmsy0 = sqrt( sum(uyr * uyr) / (nx*ny*nz) )
      rmsz0 = sqrt( sum(uzr * uzr) / (nx*ny*nz) )
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
      ignore_me = uxr(ii,jj,kk)
      rmsx = rmsx + ignore_me * ignore_me
      ignore_me = ignore_me / rmsx0
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfx(ly)=pdfx(ly)+1
 
      ignore_me = uyr(ii,jj,kk)
      rmsy = rmsy + ignore_me * ignore_me
      ignore_me = ignore_me / rmsy0
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfy(ly)=pdfy(ly)+1
 
      ignore_me = uzr(ii,jj,kk)
      rmsz = rmsz + ignore_me * ignore_me
      ignore_me = ignore_me / rmsz0
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfz(ly)=pdfz(ly)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1._sp / (nx*ny*nz) / nfile

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  rmsz = sqrt( rmsz * ignore_me )
  pdfx=pdfx/binw*ignore_me
  pdfy=pdfy/binw*ignore_me
  pdfz=pdfz/binw*ignore_me
 
  
  if (lx1 .eq. 1) then 
    path='pdfvinc-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  else if (lx1 .eq. 2) then
    path='pdfvinc-'//str(1:len_trim(str))//'dy-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  else if (lx1 .eq. 3) then
    path='pdfvinc-'//str(1:len_trim(str))//'dz-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  end if
  
  open(15,file=path)
    write(15,'(''# Title = "'',3E12.4, ''"'')') rmsx, rmsy, rmsz
    write(15,*) '# variables = "xx", "pdfdux","yy","pdfduy","zz", "pdfduz"'
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5_sp) * binw
      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, pdfx(lz) * rmsx / rmsx0, &
      ignore_me * rmsy0 / rmsy, rmsy / rmsy0 * pdfy(lz), ignore_me * rmsz0 / rmsz, &
      rmsz / rmsz0 * pdfz(lz)
    end do
  close(15)
  
  
  deallocate(ux,uy,uz,uxr,uyr,uzr)
  

  write(*,*) '3dvelincpdf.x done.'
end program vincpdf
