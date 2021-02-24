program gradpdfaij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=80
  real(sp), parameter ::  bnd = 16., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfa12,pdfa13,pdfa21
  real(dp), dimension(npnt) :: pdfa23,pdfa31,pdfa32
  real(dp), dimension(npnt) :: pdfaij
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me
  real(dp) :: rmsa12, rmsa13, rmsa21, rmsa23, rmsa31, rmsa32, rmsa120, rmsaij

  complex(sp), allocatable, dimension(:,:,:) :: a12, a13, a21
  complex(sp), allocatable, dimension(:,:,:) :: a23, a31, a32
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfaij-s(d)p.x nx filelist'
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

  allocate(a12(lx1,ly,lz),a13(lx1,ly,lz),a21(lx1,ly,lz))
  allocate(a23(lx1,ly,lz),a31(lx1,ly,lz),a32(lx1,ly,lz))
  allocate(kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfa12 = 0.d0; pdfa13 = 0.d0; pdfa21 = 0.d0
  pdfa23 = 0.d0; pdfa31 = 0.d0; pdfa32 = 0.d0
  rmsa12 = 0.d0; rmsa13 = 0.d0; rmsa21 = 0.d0
  rmsa23 = 0.d0; rmsa31 = 0.d0; rmsa32 = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a12
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a21
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a31
    close(15)
    write(*,*) 'finishing reading data'


    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a13(ii,jj,kk) = eye * kz(kk) * a12(ii,jj,kk) ! a12 = ux
      a12(ii,jj,kk) = eye * ky(jj) * a12(ii,jj,kk)
      a23(ii,jj,kk) = eye * kz(kk) * a21(ii,jj,kk) ! a21 = uy
      a21(ii,jj,kk) = eye * kx(ii) * a21(ii,jj,kk)
      a32(ii,jj,kk) = eye * ky(jj) * a31(ii,jj,kk) ! a31 = uz
      a31(ii,jj,kk) = eye * kx(ii) * a31(ii,jj,kk)
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    
    if ( nfile .eq. 0 ) then
      rmsa120 = sum( real( a12(1:lx,:,:) * conjg(a12(1:lx,:,:)) ) )
      rmsa120 = sqrt( rmsa120 / (nx * ny * nz) )
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a12(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a12(ii/2,jj,kk) )
      endif
      rmsa12 = rmsa12 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa12(ll)=pdfa12(ll)+1
 

      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a13(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a13(ii/2,jj,kk) )
      endif
      rmsa13 = rmsa13 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa13(ll)=pdfa13(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a21(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a21(ii/2,jj,kk) )
      endif
      rmsa21 = rmsa21 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa21(ll)=pdfa21(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a23(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a23(ii/2,jj,kk) )
      endif
      rmsa23 = rmsa23 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa23(ll)=pdfa23(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a31(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a31(ii/2,jj,kk) )
      endif
      rmsa31 = rmsa31 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa31(ll)=pdfa31(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ignore_me = real( a32(ii/2+1,jj,kk) )
      else
        ignore_me = aimag( a32(ii/2,jj,kk) )
      endif
      rmsa32 = rmsa32 + ignore_me * ignore_me
      ignore_me = ignore_me / rmsa120
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfa32(ll)=pdfa32(ll)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsaij = rmsa12 + rmsa13 + rmsa21 + rmsa23 + rmsa31 + rmsa32
  rmsaij = sqrt( rmsaij * ignore_me / 6. )

  rmsa12 = sqrt( rmsa12 * ignore_me )
  rmsa13 = sqrt( rmsa13 * ignore_me )
  rmsa21 = sqrt( rmsa21 * ignore_me )
  rmsa23 = sqrt( rmsa23 * ignore_me )
  rmsa31 = sqrt( rmsa31 * ignore_me )
  rmsa32 = sqrt( rmsa32 * ignore_me )

  pdfa12=pdfa12*ignore_me
  pdfa13=pdfa13*ignore_me
  pdfa21=pdfa21*ignore_me
  pdfa23=pdfa23*ignore_me
  pdfa31=pdfa31*ignore_me
  pdfa32=pdfa32*ignore_me

  write(*,*) 'check pdfa12:', sum(pdfa12)
  write(*,*) 'check pdfa13:', sum(pdfa13)
  write(*,*) 'check pdfa21:', sum(pdfa21)
  write(*,*) 'check pdfa23:', sum(pdfa23)
  write(*,*) 'check pdfa31:', sum(pdfa31)
  write(*,*) 'check pdfa32:', sum(pdfa32)

  pdfa12=pdfa12/binw
  pdfa13=pdfa13/binw
  pdfa21=pdfa21/binw
  pdfa23=pdfa23/binw
  pdfa31=pdfa31/binw
  pdfa32=pdfa32/binw

  pdfaij=(pdfa12 + pdfa13 + pdfa21 + pdfa23 + pdfa31 + pdfa32)/6.
  
  path='pdfaij-'//flnm(1:len_trim(flnm))//'.dat'//char(0)

  open(15,file=path)
    write(15,'(''#variables = "a12", "a13", "a21", "a23", "a31", "a32", "aij"'')')

    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmsa120 / rmsa12, pdfa12(lz) * rmsa12 / rmsa120, &
                            ignore_me * rmsa120 / rmsa13, pdfa13(lz) * rmsa13 / rmsa120, &
                            ignore_me * rmsa120 / rmsa21, pdfa21(lz) * rmsa21 / rmsa120, &
                            ignore_me * rmsa120 / rmsa23, pdfa23(lz) * rmsa23 / rmsa120, &
                            ignore_me * rmsa120 / rmsa31, pdfa31(lz) * rmsa31 / rmsa120, &
                            ignore_me * rmsa120 / rmsa32, pdfa32(lz) * rmsa32 / rmsa120, &
                            ignore_me * rmsa120 / rmsaij, pdfaij(lz) * rmsaij / rmsa120
    end do

  close(15)
  
  
  deallocate(a12,a13,a21,a23,a31,a32)
  deallocate(kx,ky,kz)
  
  call destroyplan3d

  write(*,*) 'done.'
end program gradpdfaij
