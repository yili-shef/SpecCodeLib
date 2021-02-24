program gradpdfaii
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=160
  real(sp), parameter ::  bnd = 16., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfa11,pdfa22,pdfa33, pdfaii
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,nfile
  real(sp) :: ignore_me
  real(dp) :: rmsa11, rmsa22, rmsa33, rmsa110, rmsaii

  complex(sp), allocatable, dimension(:,:,:) :: a11, a22, a33
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfaii-s(d)p.x nx filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
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

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

  allocate(a11(lx1,ly,lz),a22(lx1,ly,lz),a33(lx1,ly,lz))
  allocate(kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfa11 = 0.d0; pdfa22 = 0.d0; pdfa33 = 0.d0
  rmsa11 = 0.d0; rmsa22 = 0.d0; rmsa33 = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a33
    close(15)


    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a11(ii,jj,kk) = eye * kx(ii) * a11(ii,jj,kk)
      a22(ii,jj,kk) = eye * ky(jj) * a22(ii,jj,kk)
      a33(ii,jj,kk) = eye * kz(kk) * a33(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)
    
    if ( nfile .eq. 0 ) then
      rmsa110 = sum( real( a11(1:lx,:,:) * conjg(a11(1:lx,:,:)) ) )
      rmsa110 = sqrt( rmsa110 / (nx * ny * nz) )
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a11(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a11(ii/2,jj,kk) )
      endif
      rmsa11 = rmsa11 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa110
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfa11(ly)=pdfa11(ly)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a22(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a22(ii/2,jj,kk) )
      endif
      rmsa22 = rmsa22 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa110
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfa22(ly)=pdfa22(ly)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a33(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a33(ii/2,jj,kk) )
      endif
      rmsa33 = rmsa33 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa110
      ly=1+floor((ignore_me + bnd) / binw)
      if (ly .ge. 1 .and. ly .le. npnt) pdfa33(ly)=pdfa33(ly)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsaii = rmsa11 + rmsa22 + rmsa33
  rmsaii = sqrt( rmsaii * ignore_me / 3. )

  rmsa11 = sqrt( rmsa11 * ignore_me )
  rmsa22 = sqrt( rmsa22 * ignore_me )
  rmsa33 = sqrt( rmsa33 * ignore_me )

  pdfa11=pdfa11*ignore_me
  pdfa22=pdfa22*ignore_me
  pdfa33=pdfa33*ignore_me

  write(*,*) 'check pdfa11: ', sum(pdfa11)
  write(*,*) 'check pdfa22: ', sum(pdfa22)
  write(*,*) 'check pdfa33: ', sum(pdfa33)

  pdfa11=pdfa11/binw
  pdfa22=pdfa22/binw
  pdfa33=pdfa33/binw

  pdfaii=(pdfa11 + pdfa22 + pdfa33)/3.
  
  path='pdfaii-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  open(15,file=path)
    write(15,'(''#Title = "'',4E12.3, ''"'')') rmsa11, rmsa22, rmsa33, rmsaii
    write(15,'(''#variables = "a11", "pdfa11","a22","pdfa22","a33", "pdfa33", "aii", "pdfaii"'')')
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsa110 / rmsa11, pdfa11(lz) * rmsa11 / rmsa110, &
      ignore_me * rmsa110 / rmsa22, rmsa22 / rmsa110 * pdfa22(lz), ignore_me * rmsa110 / rmsa33, &
      rmsa33 / rmsa110 * pdfa33(lz), ignore_me * rmsa110 / rmsaii, rmsaii / rmsa110 * pdfaii(lz)
    end do
  close(15)
  
  
  deallocate(a11,a22,a33)
  deallocate(kx,ky,kz)
  
  call destroyplan3d

  write(*,*) 'done.'
end program gradpdfaii
