program alignpijcevtr
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
  ! ----------------------------------------------------------------------------------

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(dp) :: costh, cos1, cos2, phi

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-pijcevtr.x nx filelist prefix'
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
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'align-pijcevtr-2d-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    write(*,*) 'finishing reading p'

    p11 = - b11 * kx * kx
    p22 = - b11 * ky * ky
    p33 = - b11 * kz * kz
    p12 = - b11 * kx * ky
    p13 = - b11 * kx * kz
    p23 = - b11 * ky * kz
    b11 = - (p11 + p22 + p33)/3._sp
    p11 = p11 + b11
    p22 = p22 + b11
    p33 = p33 + b11

    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)


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
          sij(1,1) = real( p11(ll,jj,kk),sp )
          sij(1,2) = real( p12(ll,jj,kk),sp )
          sij(1,3) = real( p13(ll,jj,kk),sp )
          sij(2,2) = real( p22(ll,jj,kk),sp )
          sij(2,3) = real( p23(ll,jj,kk),sp )
          sij(3,3) = real( p33(ll,jj,kk),sp )

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
          sij(1,1) = aimag( p11(ll,jj,kk) )
          sij(1,2) = aimag( p12(ll,jj,kk) )
          sij(1,3) = aimag( p13(ll,jj,kk) )
          sij(2,2) = aimag( p22(ll,jj,kk) )
          sij(2,3) = aimag( p23(ll,jj,kk) )
          sij(3,3) = aimag( p33(ll,jj,kk) )

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
      phi = acos( abs( cos1 / (cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      if (ll .ge. 1 .and. ll .le. npnt .and.  &
          mm .ge. 1 .and. mm .le. npnt) then
          pdfsc_alfa(ll,mm)=pdfsc_alfa(ll,mm)+1
      end if

      costh = abs( dot_product(evtrcc(:,2), evtrss(:,3)) )
      cos1  = dot_product(evtrcc(:,2), evtrss(:,2))
      cos2  = dot_product(evtrcc(:,2), evtrss(:,1))
      phi = acos( abs( cos1 / (cos1*cos1 + cos2*cos2) ) )

      ll = floor(costh/binw) + 1
      mm = floor(phi/binwphi) + 1
      if (ll .ge. 1 .and. ll .le. npnt .and.  &
          mm .ge. 1 .and. mm .le. npnt) then
          pdfsc_beta(ll,mm)=pdfsc_beta(ll,mm)+1
      end if

      costh = abs( dot_product(evtrcc(:,1), evtrss(:,3)) )
      cos1  = dot_product(evtrcc(:,1), evtrss(:,2))
      cos2  = dot_product(evtrcc(:,1), evtrss(:,1))
      phi = acos( abs( cos1 / (cos1*cos1 + cos2*cos2) ) )

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

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,p12,p13,p23)
  
  call destroyplan3d

  write(*,*) 'align-pijcevtr.x done.'

end program alignpijcevtr
