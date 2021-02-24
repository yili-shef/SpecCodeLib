program alignic
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 50
  real(sp), parameter :: bnd = 1._sp, binw = bnd/npnt
  real(sp), parameter :: bndphi = pi/2._sp, binwphi = bndphi/npnt
  real(dp), dimension(npnt,npnt) :: pdfic2d
  real(dp), dimension(npnt) :: pdfalfa, pdfbeta, pdfgmma


  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------
  real(dp), dimension(evsize,evsize)  :: cc, aa, bb, evtrcc
  real(dp), dimension(evsize) :: evcc, angm
  real(dp) :: costh, cos1, cos2, phi, angma, angmb, angmg
  real(dp), parameter :: safety = 8._dp
  integer :: pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: a11,a22,a33,a12,a13,a23,a21,a31,a32
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-ic.x nx filelist prefix'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
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

  ! prefix
  call getarg(3,prf)
  prf = adjustl(prf)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
  allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
  allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'align-ic2d-'//trim(flnm)//'.dat')
  open(26, file = 'align-ic1d-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b33
    close(15)
    write(*,*) 'finishing reading ux uy uz'

    a11 = eye * kx * b11; a12 = eye * ky * b11; a13 = eye * kz * b11
    a21 = eye * kx * b22; a22 = eye * ky * b22; a23 = eye * kz * b22
    a31 = eye * kx * b33; a32 = eye * ky * b33; a33 = eye * kz * b33 

    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)


    open(15,file='./out/b11'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    open(15,file='./out/b12'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b12
    close(15)
    open(15,file='./out/b13'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b13
    close(15)
    open(15,file='./out/b21'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b21
    close(15)
    open(15,file='./out/b22'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b22
    close(15)
    open(15,file='./out/b23'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b23
    close(15)
    open(15,file='./out/b31'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b31
    close(15)
    open(15,file='./out/b32'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b32
    close(15)
    open(15,file='./out/b33'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b33
    close(15)
    write(*,*) 'finishing reading bij'


    pdfic2d = 0._dp; pdfalfa = 0._dp; pdfbeta = 0._dp; pdfgmma = 0._dp
    pointcnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          aa(1,1) = real( a11(ll,jj,kk),sp )
          aa(1,2) = real( a12(ll,jj,kk),sp )
          aa(1,3) = real( a13(ll,jj,kk),sp )
          aa(2,1) = real( a21(ll,jj,kk),sp )
          aa(2,2) = real( a22(ll,jj,kk),sp )
          aa(2,3) = real( a23(ll,jj,kk),sp )
          aa(3,1) = real( a31(ll,jj,kk),sp )
          aa(3,2) = real( a32(ll,jj,kk),sp )
          aa(3,3) = real( a33(ll,jj,kk),sp )

          cc(1,1) = real( b11(ll,jj,kk),sp )
          cc(1,2) = real( b12(ll,jj,kk),sp )
          cc(1,3) = real( b13(ll,jj,kk),sp )
          cc(2,1) = real( b21(ll,jj,kk),sp )
          cc(2,2) = real( b22(ll,jj,kk),sp )
          cc(2,3) = real( b23(ll,jj,kk),sp )
          cc(3,1) = real( b31(ll,jj,kk),sp )
          cc(3,2) = real( b32(ll,jj,kk),sp )
          cc(3,3) = real( b33(ll,jj,kk),sp )
      else
          ll = ii/2

          aa(1,1) = aimag( a11(ll,jj,kk) )
          aa(1,2) = aimag( a12(ll,jj,kk) )
          aa(1,3) = aimag( a13(ll,jj,kk) )
          aa(2,1) = aimag( a21(ll,jj,kk) )
          aa(2,2) = aimag( a22(ll,jj,kk) )
          aa(2,3) = aimag( a23(ll,jj,kk) )
          aa(3,1) = aimag( a31(ll,jj,kk) )
          aa(3,2) = aimag( a32(ll,jj,kk) )
          aa(3,3) = aimag( a33(ll,jj,kk) )

          cc(1,1) = aimag( b11(ll,jj,kk) )
          cc(1,2) = aimag( b12(ll,jj,kk) )
          cc(1,3) = aimag( b13(ll,jj,kk) )
          cc(2,1) = aimag( b21(ll,jj,kk) )
          cc(2,2) = aimag( b22(ll,jj,kk) )
          cc(2,3) = aimag( b23(ll,jj,kk) )
          cc(3,1) = aimag( b31(ll,jj,kk) )
          cc(3,2) = aimag( b32(ll,jj,kk) )
          cc(3,3) = aimag( b33(ll,jj,kk) )
      endif
      cc = matmul( cc, transpose(cc) )

      ! Angular momentum
      bb = matmul(cc, transpose(aa))
      angm(1) = bb(2,3) - bb(3,2)
      angm(2) = bb(3,1) - bb(1,3)
      angm(3) = bb(1,2) - bb(2,1)
      angm = angm / sqrt( dot_product(angm, angm) )

      ! The interface of subroutine dsyevr
      ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
      ! lwork, iwork, liwork, info)

      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info
      ! The eigenvalues are in ascending order: smallest first.

      if (      abs(evcc(1) - evcc(2)) .le. safety * dlamch('e') * abs(evcc(2)) &
           .or. abs(evcc(2) - evcc(3)) .le. safety * dlamch('e') * abs(evcc(3)) &
           .or. abs(evcc(3) - evcc(1)) .le. safety * dlamch('e') * abs(evcc(3)) ) then 
          cycle
      else
          pointcnt = pointcnt + 1
      end if

      angma = abs( dot_product(angm, evtrcc(:,3)) )
      angmb = abs( dot_product(angm, evtrcc(:,2)) )
      angmg = abs( dot_product(angm, evtrcc(:,1)) )

      ll = floor(angma/binw) + 1
      pdfalfa(ll) = pdfalfa(ll) + 1

      ll = floor(angmb/binw) + 1
      pdfbeta(ll) = pdfbeta(ll) + 1

      ll = floor(angmg/binw) + 1
      pdfgmma(ll) = pdfgmma(ll) + 1

      costh = abs( dot_product(evtrcc(:,3), angm) )
      cos1  = dot_product(evtrcc(:,2), angm)
      cos2  = dot_product(evtrcc(:,1), angm)
      phi = acos( abs( cos1 / sqrt(cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      pdfic2d(ll,mm)=pdfic2d(ll,mm)+1

    end do
    end do
    end do
    pdfic2d = pdfic2d * const / binw / binwphi
    pdfalfa = pdfalfa * const / binw
    pdfbeta = pdfbeta * const / binw
    pdfgmma = pdfgmma * const / binw

    write(*,*) 'point counted is: ', pointcnt * const

    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
    do mm = 1, npnt
      write(27, '(15E15.3)') (ll - 0.5_sp) * binw, (mm-0.5_sp) * binwphi, pdfic2d(ll,mm)
    end do
    write(27, *)
    end do
    write(27, *)
    write(27, *)

    write(26, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(26, '(15E15.3)') (ll - 0.5_sp) * binw, pdfalfa(ll), pdfbeta(ll), pdfgmma(ll)
    end do
    write(26, *)
    write(26, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(26)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,a11,a22,a33,a12,a13,a23,a21,a31,a32)
  
  call destroyplan3d

  write(*,*) 'align-ic.x done.'

end program alignic
