program alignsc
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 50
  real(sp), parameter :: bnd = 1._sp, binw = bnd/npnt
  real(sp), parameter :: bndphi = pi/2._sp, binwphi = bndphi/npnt
  real(dp), dimension(npnt,npnt) :: pdfsc_alfa, pdfsc_beta, pdfsc_gamma


  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  real(dp), dimension(evsize,evsize)  :: cc, sij, evtrcc, evtrss
  real(dp), dimension(evsize) :: evcc, evss
  real(dp) :: costh, cos1, cos2, phi
  ! ----------------------------------------------------------------------------------

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s22,s33,s12,s13,s23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-sc.x nx filelist prefix'
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
  allocate( s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate( s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'align-sijcijevtr-2d-'//trim(flnm)//'.dat')

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

    s11 = eye * kx * b11
    s22 = eye * ky * b22
    s33 = eye * kz * b33
    s12 = eye * (kx * b22 + ky * b11)*0.5_sp
    s13 = eye * (kx * b33 + kz * b11)*0.5_sp
    s23 = eye * (ky * b33 + kz * b22)*0.5_sp
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)


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


    pdfsc_alfa = 0._dp; pdfsc_beta = 0._dp; pdfsc_gamma = 0._dp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          sij(1,1) = real( s11(ll,jj,kk),sp )
          sij(1,2) = real( s12(ll,jj,kk),sp )
          sij(1,3) = real( s13(ll,jj,kk),sp )
          sij(2,2) = real( s22(ll,jj,kk),sp )
          sij(2,3) = real( s23(ll,jj,kk),sp )
          sij(3,3) = real( s33(ll,jj,kk),sp )

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
          sij(1,1) = aimag( s11(ll,jj,kk) )
          sij(1,2) = aimag( s12(ll,jj,kk) )
          sij(1,3) = aimag( s13(ll,jj,kk) )
          sij(2,2) = aimag( s22(ll,jj,kk) )
          sij(2,3) = aimag( s23(ll,jj,kk) )
          sij(3,3) = aimag( s33(ll,jj,kk) )

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
      sij(2,1) = sij(1,2); sij(3,1) = sij(1,3); sij(3,2) = sij(2,3)

      ! The interface of subroutine dsyevr
      ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
      ! lwork, iwork, liwork, info)

      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info
     
      call dsyevr("V", "A", "U", evsize, sij, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evss, evtrss, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev sij. info = ', info
      ! The eigenvalues are in ascending order: smallest first.

      costh = abs( dot_product(evtrcc(:,3), evtrss(:,3)) )
      cos1  = dot_product(evtrcc(:,3), evtrss(:,2))
      cos2  = dot_product(evtrcc(:,3), evtrss(:,1))
      phi = acos( abs( cos1 / sqrt(cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      if (ll .ge. 1 .and. ll .le. npnt .and.  &
          mm .ge. 1 .and. mm .le. npnt) then
          pdfsc_alfa(ll,mm)=pdfsc_alfa(ll,mm)+1
      end if

      costh = abs( dot_product(evtrcc(:,2), evtrss(:,3)) )
      cos1  = dot_product(evtrcc(:,2), evtrss(:,2))
      cos2  = dot_product(evtrcc(:,2), evtrss(:,1))
      phi = acos( abs( cos1 / sqrt(cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      if (ll .ge. 1 .and. ll .le. npnt .and.  &
          mm .ge. 1 .and. mm .le. npnt) then
          pdfsc_beta(ll,mm)=pdfsc_beta(ll,mm)+1
      end if

      costh = abs( dot_product(evtrcc(:,1), evtrss(:,3)) )
      cos1  = dot_product(evtrcc(:,1), evtrss(:,2))
      cos2  = dot_product(evtrcc(:,1), evtrss(:,1))
      phi = acos( abs( cos1 / sqrt(cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      if (ll .ge. 1 .and. ll .le. npnt .and.  &
          mm .ge. 1 .and. mm .le. npnt) then
          pdfsc_gamma(ll,mm)=pdfsc_gamma(ll,mm)+1
      end if

    end do
    end do
    end do
    pdfsc_alfa = pdfsc_alfa * const / binw / binwphi
    pdfsc_beta = pdfsc_beta * const / binw / binwphi
    pdfsc_gamma = pdfsc_gamma * const / binw / binwphi

    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
    do mm = 1, npnt
      write(27, '(15E15.3)') (ll - 0.5_sp) * binw, (mm-0.5_sp) * binwphi, pdfsc_alfa(ll,mm), &
                             pdfsc_beta(ll,mm), pdfsc_gamma(ll,mm)
    end do
    write(27, *)
    end do
    write(27, *)
    write(27, *)


    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,s11,s22,s33,s12,s13,s23)
  
  call destroyplan3d

  write(*,*) 'align-sc.x done.'

end program alignsc
