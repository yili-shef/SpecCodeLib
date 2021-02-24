program gradpdf
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfa11,pdfa12,pdfa13
  real(dp), dimension(npnt) :: pdfa21,pdfa22,pdfa23
  real(dp), dimension(npnt) :: pdfa31,pdfa32,pdfa33
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ndel, nfile
  real(sp) :: ignore_me
  real(dp) :: rmsx, rmsy, rmsz, rmsx0, rmsy0, rmsz0

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  real(sp),    allocatable, dimension(:,:,:) :: k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./gradpdf.x nx ndel filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ny=nx; nz=nx
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz

  allocate(a11(lx1,ly,lz),a12(lx1,ly,lz),a13(lx1,ly,lz))
  allocate(a21(lx1,ly,lz),a22(lx1,ly,lz),a23(lx1,ly,lz))
  allocate(a31(lx1,ly,lz),a32(lx1,ly,lz),a33(lx1,ly,lz))
  allocate(k2 (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

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


  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfx = 0.d0; pdfy = 0.d0; pdfz = 0.d0
  rmsx = 0.d0; rmsy = 0.d0; rmsz = 0.d0
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
      a12(ii,jj,kk) = eye * ky(jj) * a11(ii,jj,kk) 
      a13(ii,jj,kk) = eye * kz(kk) * a11(ii,jj,kk)
      a21(ii,jj,kk) = eye * kx(ii) * a22(ii,jj,kk)
      a23(ii,jj,kk) = eye * kz(kk) * a22(ii,jj,kk)
      a31(ii,jj,kk) = eye * kx(ii) * a33(ii,jj,kk)
      a32(ii,jj,kk) = eye * ky(jj) * a33(ii,jj,kk)
      a11(ii,jj,kk) = eye * kx(ii) * a11(ii,jj,kk)
      a22(ii,jj,kk) = eye * ky(jj) * a22(ii,jj,kk)
      a33(ii,jj,kk) = eye * kz(kk) * a33(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)
    
    
    if ( nfile .eq. 0 ) then
      rmsa110 = sum( real( a11(1:lx,:,:) * conjg(a11(1:lx,:,:)) ) )
      rmsa110 = sqrt( rmsa110 / (nx * ny * nz) )

      rmsa120 = sum( real( a12(1:lx,:,:) * conjg(a12(1:lx,:,:)) ) )
      rmswz0 = rmsa120
      rmsa120 = sqrt( rmsa120 / (nx * ny * nz) )

      rmsa130 = sum( real( a13(1:lx,:,:) * conjg(a13(1:lx,:,:)) ) )
      rmswy0 = rmsa130
      rmsa130 = sqrt( rmsa130 / (nx * ny * nz) )

      rmsa210 = sum( real( a21(1:lx,:,:) * conjg(a21(1:lx,:,:)) ) )
      rmswz0 = rmsz0 + rmsa210
      rmsa210 = sqrt( rmsa210 / (nx * ny * nz) )

      rmsa220 = sum( real( a22(1:lx,:,:) * conjg(a22(1:lx,:,:)) ) )
      rmsa220 = sqrt( rmsa220 / (nx * ny * nz) )

      rmsa230 = sum( real( a23(1:lx,:,:) * conjg(a23(1:lx,:,:)) ) )
      rmswx0 = rmsa230
      rmsa230 = sqrt( rmsa230 / (nx * ny * nz) )

      rmsa310 = sum( real( a31(1:lx,:,:) * conjg(a31(1:lx,:,:)) ) )
      rmswy0 = rmswy0 + rmsa310
      rmsa310 = sqrt( rmsa310 / (nx * ny * nz) )

      rmsa320 = sum( real( a32(1:lx,:,:) * conjg(a32(1:lx,:,:)) ) )
      rmswx0 = rmswx0 + rmsa320
      rmsa320 = sqrt( rmsa320 / (nx * ny * nz) )

      rmsa330 = sum( real( a33(1:lx,:,:) * conjg(a33(1:lx,:,:)) ) )
      rmsa330 = sqrt( rmsa330 / (nx * ny * nz) )

      rmswx0 = rmswx0 - 2 * sum( real( a32(1:lx,:,:) ) * real ( a23(1:lx,:,:) ) ) &
                  - 2 * sum( aimag( a32(1:lx,:,:) ) * aimag( a23(1:lx,:,:) ) )                    

      rmswy0 = rmswy0 - 2 * sum( real( a31(1:lx,:,:) ) * real ( a13(1:lx,:,:) ) ) &
                  - 2 * sum( aimag( a31(1:lx,:,:) ) * aimag( a13(1:lx,:,:) ) )                    

      rmswz0 = rmswz0 - 2 * sum( real( a21(1:lx,:,:) ) * real ( a12(1:lx,:,:) ) ) &
                  - 2 * sum( aimag( a21(1:lx,:,:) ) * aimag( a12(1:lx,:,:) ) )                    
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

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  rmsz = sqrt( rmsz * ignore_me )
  pdfx=pdfx/binw*ignore_me
  pdfy=pdfy/binw*ignore_me
  pdfz=pdfz/binw*ignore_me
 
  
  if (lx1 .eq. 1) then 
    path='pdf-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  else if (lx1 .eq. 2) then
    path='pdf-'//str(1:len_trim(str))//'dy-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  else if (lx1 .eq. 3) then
    path='pdf-'//str(1:len_trim(str))//'dz-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  end if
  
  open(15,file=path)
    write(15,'(''Title = "'',3E12.4, ''"'')') rmsx, rmsy, rmsz
    write(15,*) 'variables = "xx", "pdfdux","yy","pdfduy","zz", "pdfduz"'
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, pdfx(lz) * rmsx / rmsx0, &
      ignore_me * rmsy0 / rmsy, rmsy / rmsy0 * pdfy(lz), ignore_me * rmsz0 / rmsz, &
      rmsz / rmsz0 * pdfz(lz)
    end do
  close(15)
  
  
  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33)
  deallocate(kx,ky,kz,k2)
  
  call destroyplan3d

  write(*,*) 'done.'
end program gradpdf
