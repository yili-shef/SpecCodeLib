program gradpdfaij
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=160
  real(sp), parameter ::  bnd = 16., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfr11,pdfr12,pdfr13
  real(dp), dimension(npnt) :: pdfr22,pdfr23,pdfr33
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me
  real(dp) :: rmsr11, rmsr12, rmsr13, rmsr22, rmsr23, rmsr33, rmsr110, rmsaij

  complex(sp), allocatable, dimension(:,:,:) :: r11, r12, r13
  complex(sp), allocatable, dimension(:,:,:) :: r22, r23, r33
  real(sp),    allocatable, dimension(:,:,:) :: k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfrij.x nx filelist'
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

  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(k2 (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfr11 = 0.d0; pdfr12 = 0.d0; pdfr13 = 0.d0
  pdfr22 = 0.d0; pdfr23 = 0.d0; pdfr33 = 0.d0
  rmsr11 = 0.d0; rmsr12 = 0.d0; rmsr13 = 0.d0
  rmsr22 = 0.d0; rmsr23 = 0.d0; rmsr33 = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) r12
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) r13
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) r23
    close(15)
    write(*,*) 'finishing reading data'


    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      r11(ii,jj,kk)=eye*(ky(jj)*r13(ii,jj,kk)-kz(kk)*r12(ii,jj,kk))
      r22(ii,jj,kk)=eye*(kz(kk)*r12(ii,jj,kk)-kx(ii)*r13(ii,jj,kk))
      r33(ii,jj,kk)=eye*(kx(ii)*r13(ii,jj,kk)-ky(jj)*r12(ii,jj,kk))

      r12(ii,jj,kk)=.5*eye*(kx(ii)*r22(ii,jj,kk)+ky(jj)*r11(ii,jj,kk))
      r13(ii,jj,kk)=.5*eye*(kx(ii)*r33(ii,jj,kk)+kz(kk)*r11(ii,jj,kk))
      r23(ii,jj,kk)=.5*eye*(ky(jj)*r33(ii,jj,kk)+kz(kk)*r22(ii,jj,kk))
      r11(ii,jj,kk)=eye*kx(ii)*r11(ii,jj,kk)
      r22(ii,jj,kk)=eye*ky(jj)*r22(ii,jj,kk)
      r33(ii,jj,kk)=eye*kz(kk)*r33(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
    
    if ( nfile .eq. 0 ) then
      rmsr110 = sum( real( r11(1:lx,:,:) * conjg(r11(1:lx,:,:)) ) )
      rmsr110 = sqrt( rmsr110 / (nx * ny * nz) )
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r11(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r11(ii/2,jj,kk) )
      endif
      rmsr11 = rmsr11 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr11(ll)=pdfr11(ll)+1
 

      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r12(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r12(ii/2,jj,kk) )
      endif
      rmsr12 = rmsr12 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr12(ll)=pdfr12(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r13(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r13(ii/2,jj,kk) )
      endif
      rmsr13 = rmsr13 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr13(ll)=pdfr13(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r22(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r22(ii/2,jj,kk) )
      endif
      rmsr22 = rmsr22 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr22(ll)=pdfr22(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r23(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r23(ii/2,jj,kk) )
      endif
      rmsr23 = rmsr23 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr23(ll)=pdfr23(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( r33(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( r33(ii/2,jj,kk) )
      endif
      rmsr33 = rmsr33 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsr110
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfr33(ll)=pdfr33(ll)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsaij = rmsr11 + rmsr12 + rmsr13 + rmsr22 + rmsr23 + rmsr33
  rmsaij = sqrt( rmsaij * ignore_me / 6. )

  rmsr11 = sqrt( rmsr11 * ignore_me )
  rmsr12 = sqrt( rmsr12 * ignore_me )
  rmsr13 = sqrt( rmsr13 * ignore_me )
  rmsr22 = sqrt( rmsr22 * ignore_me )
  rmsr23 = sqrt( rmsr23 * ignore_me )
  rmsr33 = sqrt( rmsr33 * ignore_me )

  pdfr11=pdfr11*ignore_me
  pdfr12=pdfr12*ignore_me
  pdfr13=pdfr13*ignore_me
  pdfr22=pdfr22*ignore_me
  pdfr23=pdfr23*ignore_me
  pdfr33=pdfr33*ignore_me

  write(*,*) 'check pdfr11:', sum(pdfr11)
  write(*,*) 'check pdfr12:', sum(pdfr12)
  write(*,*) 'check pdfr13:', sum(pdfr13)
  write(*,*) 'check pdfr22:', sum(pdfr22)
  write(*,*) 'check pdfr23:', sum(pdfr23)
  write(*,*) 'check pdfr33:', sum(pdfr33)

  pdfr11=pdfr11/binw
  pdfr12=pdfr12/binw
  pdfr13=pdfr13/binw
  pdfr22=pdfr22/binw
  pdfr23=pdfr23/binw
  pdfr33=pdfr33/binw

  path='pdfrij-'//flnm(1:len_trim(flnm))//'.dat'//char(0)

  open(15,file=path)
    write(15,'(''variables = "aij", "pdfaij"'')')

    write(15,'(''zone T = " r11 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr11, npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr11, pdfr11(lz) * rmsr11 / rmsr110
    end do

    write(15,'(''zone T = " r12 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr12,  npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr12, pdfr12(lz) * rmsr12 / rmsr110
    end do

    write(15,'(''zone T = " r13 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr13,  npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr13, pdfr13(lz) * rmsr13 / rmsr110
    end do

    write(15,'(''zone T = " r22 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr22,  npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr22, pdfr22(lz) * rmsr22 / rmsr110
    end do

    write(15,'(''zone T = " r23 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr23,  npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr23, pdfr23(lz) * rmsr23 / rmsr110
    end do

    write(15,'(''zone T = " r33 '',1E12.3, ''", I='', I4, '' F = point'')') rmsr33,  npnt
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsr110 / rmsr33, pdfr33(lz) * rmsr33 / rmsr110
    end do

  close(15)
  
  
  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(kx,ky,kz,k2)
  
  call destroyplan3d

  write(*,*) 'done.'
end program gradpdfaij
