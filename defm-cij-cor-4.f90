program defmcijcor
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile,npnt
  real(sp) :: const

  integer, parameter :: naverage = 20
  real(sp), dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll
  integer :: iiloc, jjloc, kkloc

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

  real(dp) :: mshape, shapecc0, shapecc
  real(dp), parameter :: safety = 8._dp
  integer :: pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: a11,a22,a33,a12,a13,a23,a21,a31,a32
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  real(sp),    allocatable, dimension(:)     :: corelshape
  integer,     allocatable, dimension(:)     :: pntcnt

  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-cij-cor-2.x nx filelist prefix'
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

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate(corelshape(npnt), pntcnt(npnt))

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
  allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
  allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(26, file = 'defm-shape-cor-'//trim(flnm)//'.dat')
  open(27, file = 'defm-mshape-'//trim(flnm)//'.dat')

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

    mm = -100
    do ll = 1, naverage

      ii = floor( ran1(mm) * nx ) + 1
      jj = floor( ran1(mm) * ny ) + 1
      kk = floor( ran1(mm) * nz ) + 1

      xx(ll) = ii; yy(ll) = jj; zz(ll) = kk

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
                  nfound, evcc0(:,ll), evtrcc0(:,:,ll), evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      evcc0(:,ll) = log(evcc0(:,ll))
      shapecc0 = - 3 * sqrt(6._sp) * evcc0(1,ll) * evcc0(2,ll) * evcc0(3,ll)
      evcc0(1,ll) = shapecc0 / ( sum(evcc0(:,ll) * evcc0(:,ll)) )**1.5_sp
         
    end do


    corelshape = 0._dp; mshape = 0._dp
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

         
      evcc = log(evcc)
      shapecc = - 3 * sqrt(6._sp) * evcc(1) * evcc(2) * evcc(3)
      shapecc = shapecc / ( sum(evcc * evcc) )**1.5_sp

      do ll = 1, naverage 

        do kkloc = -1, 1 !checking the images of the flow field
        do jjloc = -1, 1
        do iiloc = -1, 1

          xxll = xx(ll) + iiloc * nx
          yyll = yy(ll) + jjloc * ny
          zzll = zz(ll) + kkloc * nz

          mm = floor( sqrt( (ii-xxll)**2 + (jj-yyll)**2 + (kk-zzll)**2 ) + 0.5_sp )  + 1
         
          if ( mm .ge. 1 .and. mm .le. npnt) then
         
            corelshape(mm) = corelshape(mm) + shapecc * evcc0(1,ll)
            pntcnt(mm) = pntcnt(mm) + 1
         
          end if

        end do
        end do
        end do

 
      end do

      mshape = mshape + shapecc

    end do
    end do
    end do

    write(*,*) 'point counted is: ', pointcnt * const

    mshape = mshape * const
    corelshape  = corelshape  / (pntcnt + mytiny)

    write(26, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(26, '(I6, 15E15.3)') ll-1, corelshape(ll)
    end do
    write(26, *)
    write(26, *)

    write(27, '(I6, 15E15.3)') nfile, mshape

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(26)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,a11,a22,a33,a12,a13,a23,a21,a31,a32)
  deallocate(corelshape, pntcnt)
  
  write(*,*) 'defm-cij-cor-2.x done.'

end program defmcijcor
