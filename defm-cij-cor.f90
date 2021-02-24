program defmcijcor
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile,npnt
  real(sp) :: ignore_me, const

  integer, parameter :: naverage = 20
  real(sp), dimension(naverage) :: xx, yy, zz

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------
  external :: ran1
  real(sp) :: ran1

  real(dp), dimension(evsize,evsize)  :: cc, cc0, evtrcc
  real(dp), dimension(evsize) :: evcc

  real(dp), dimension(evsize, naverage) :: evcc0
  real(dp), dimension(evsize, evsize, naverage) :: evtrcc0

  real(dp) :: costh, cos1, cos2, phi, mlogal
  real(dp), parameter :: safety = 8._dp
  integer :: pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: a11,a22,a33,a12,a13,a23,a21,a31,a32
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  real(sp),    allocatable, dimension(:)     :: malcosth, malphi, corelal
  integer,     allocatable, dimension(:)     :: pntcnt

  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-cij-cor.x nx filelist prefix'
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

  npnt = floor( sqrt( real(nx-1)**2 + real(ny-1)**2 + real(nz-1)**2 ) + 0.5_sp ) + 1
  allocate(malcosth(npnt), malphi(npnt), corelal(npnt), pntcnt(npnt))

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

  open(27, file = 'defm-cij-cor-mlogal-'//trim(flnm)//'.dat')
  open(26, file = 'defm-cij-cor-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))


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

    ll = -100
    do mm = 1, naverage

      ii = floor( ran1(ll) * nx ) + 1
      jj = floor( ran1(ll) * ny ) + 1
      kk = floor( ran1(ll) * nz ) + 1

      xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

      cc0(1,1) = real( b11(ii/2+1,jj,kk),sp )
      cc0(1,2) = real( b12(ii/2+1,jj,kk),sp )
      cc0(1,3) = real( b13(ii/2+1,jj,kk),sp )
      cc0(2,1) = real( b21(ii/2+1,jj,kk),sp )
      cc0(2,2) = real( b22(ii/2+1,jj,kk),sp )
      cc0(2,3) = real( b23(ii/2+1,jj,kk),sp )
      cc0(3,1) = real( b31(ii/2+1,jj,kk),sp )
      cc0(3,2) = real( b32(ii/2+1,jj,kk),sp )
      cc0(3,3) = real( b33(ii/2+1,jj,kk),sp )
      cc0 = matmul(cc0, transpose(cc0))

      ! The eigenvalues are in ascending order: smallest first.
      call dsyevr("V", "A", "U", evsize, cc0, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc0(:,mm), evtrcc0(:,:,mm), evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

    end do


    malcosth = 0._dp; malphi = 0._dp; corelal = 0._dp; mlogal = 0._dp
    pntcnt = 0; pointcnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
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

      do ll = 1, naverage

        mm = floor( sqrt( (ii-xx(ll))**2 + (jj-yy(ll))**2 + (kk-zz(ll))**2 ) + 0.5_sp )  + 1
 
        costh = abs( dot_product( evtrcc(:,3), evtrcc0(:,3,ll) ) )
        cos1  = dot_product( evtrcc(:,3), evtrcc0(:,2,ll) )
        cos2  = dot_product( evtrcc(:,3), evtrcc0(:,1,ll) )
        phi = acos( abs(cos1) / sqrt( cos1 * cos1 + cos2 * cos2 ) )
 
        malcosth(mm) = malcosth(mm) + costh
        malphi(mm)   = malphi(mm)   + phi

        corelal(mm) = corelal(mm) + log(evcc(3)) * log(evcc0(3,ll))

        pntcnt(mm) = pntcnt(mm) + 1
 
      end do
      mlogal = mlogal + log(evcc(3))

    end do
    end do
    end do

    write(*,*) 'point counted is: ', pointcnt * const


    mlogal = mlogal * const

    malcosth = malcosth / (pntcnt + mytiny)
    malphi   = malphi   / (pntcnt + mytiny)
    corelal  = corelal  / (pntcnt + mytiny)

    write(26, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(26, '(I6, 15E15.3)') ll-1, malcosth(ll), malphi(ll), corelal(ll)
    end do
    write(26, *)
    write(26, *)

    write(27, '(I6, 15E15.3)') nfile, mlogal

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(26)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,a11,a22,a33,a12,a13,a23,a21,a31,a32)
  deallocate(malcosth, malphi, corelal, pntcnt)
  
  call destroyplan3d

  write(*,*) 'defm-cij-cor.x done.'

end program defmcijcor
