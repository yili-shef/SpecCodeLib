program defmcijcor
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,nfile,npnt
  real(sp) :: const

  integer, parameter :: naverage = 200
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll
  integer  :: iiloc, jjloc, kkloc

  integer, parameter :: npdf = 50
  real(sp), parameter :: binw = 1._sp/npdf

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------
  real(sp) :: ran1
  external ran1

  real(dp), dimension(evsize,evsize)  :: cc, evtrcc
  real(dp), dimension(evsize) :: evcc

  real(dp), dimension(evsize, naverage) :: evcc0
  real(dp), dimension(evsize, evsize, naverage) :: evtrcc0

  real(dp) :: costh
  real(dp), parameter :: safety = 8._dp
  integer :: pointcnt, noffset

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  integer,     allocatable, dimension(:)     :: pntcnt
  real(dp), allocatable, dimension(:,:) :: pcthal, pcthbe, pcthgm

  character(80) :: str, flnm, prf, str1, str4
  
  if (iargc() .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-cij-alignpdf.x nx filelist prefix noffset'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     noffset: offset in file number in the list'
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

  call getarg(4,str4)
  read(str4, '(I20)') noffset

  write(*,*) 'Running defm-cij-alignpdf.x '//trim(str)//' '//trim(flnm)//' '//trim(prf)//' '//trim(str4)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), pcthal(npdf,npnt), pcthbe(npdf,npnt), pcthgm(npdf,npnt))
  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )

  ! Will calculate the pdf for the "noffset" data file in each list only.
  ! For this we need to generate the file names instead of reading them from the lists. 
  ! Will use the smartrunner(.py) to run the code, so the interface is kept the same.

  str1 = flnm( len_trim(flnm) - 3 : len_trim(flnm) ) !last 4 letters are the starting data file number 
  read(str1, '(I20)') ll
  ll = ll + noffset
  write(str1, '(I20)') ll 
  str1 = adjustl(str1) ! The number of the last data file in this series
  flnm = trim(flnm)//'-'//trim(str1)//'.dat'


  open(15,file='./out/b11'//trim(flnm),form='unformatted')
    read(15) b11
  close(15)
  open(15,file='./out/b12'//trim(flnm),form='unformatted')
    read(15) b12
  close(15)
  open(15,file='./out/b13'//trim(flnm),form='unformatted')
    read(15) b13
  close(15)
  open(15,file='./out/b21'//trim(flnm),form='unformatted')
    read(15) b21
  close(15)
  open(15,file='./out/b22'//trim(flnm),form='unformatted')
    read(15) b22
  close(15)
  open(15,file='./out/b23'//trim(flnm),form='unformatted')
    read(15) b23
  close(15)
  open(15,file='./out/b31'//trim(flnm),form='unformatted')
    read(15) b31
  close(15)
  open(15,file='./out/b32'//trim(flnm),form='unformatted')
    read(15) b32
  close(15)
  open(15,file='./out/b33'//trim(flnm),form='unformatted')
    read(15) b33
  close(15)

  ll = -100
  do mm = 1, naverage

    ii = floor( ran1(ll) * nx ) + 1
    jj = floor( ran1(ll) * ny ) + 1
    kk = floor( ran1(ll) * nz ) + 1

    xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

    cc(1,1) = real( b11(ii/2+1,jj,kk),sp )
    cc(1,2) = real( b12(ii/2+1,jj,kk),sp )
    cc(1,3) = real( b13(ii/2+1,jj,kk),sp )
    cc(2,1) = real( b21(ii/2+1,jj,kk),sp )
    cc(2,2) = real( b22(ii/2+1,jj,kk),sp )
    cc(2,3) = real( b23(ii/2+1,jj,kk),sp )
    cc(3,1) = real( b31(ii/2+1,jj,kk),sp )
    cc(3,2) = real( b32(ii/2+1,jj,kk),sp )
    cc(3,3) = real( b33(ii/2+1,jj,kk),sp )
    cc = matmul(cc, transpose(cc))

    ! The eigenvalues are in ascending order: smallest first.
    call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                nfound, evcc0(:,mm), evtrcc0(:,:,mm), evsize, isuppz, work, lwork, iwork, liwork, info)
    if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

  end do


  pcthal = 0._dp; pcthbe = 0._dp; pcthgm = 0._dp
  pntcnt = 0
  pointcnt = 0
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


    ! loop over different groups
    do ll = 1, naverage

      do kkloc = -1, 1 ! checking the images of the flow field
      do jjloc = -1, 1
      do iiloc = -1, 1

        xxll = xx(ll) + iiloc * nx
        yyll = yy(ll) + jjloc * ny
        zzll = zz(ll) + kkloc * nz

        mm = floor( sqrt( (ii-xxll)**2 + (jj-yyll)**2 + (kk-zzll)**2 ) + 0.5_sp )  + 1
 
        if ( mm .ge. 1 .and. mm .le. npnt) then
  
          costh = abs( dot_product( evtrcc(:,3), evtrcc0(:,3,ll) ) )
          pp = floor( costh / binw ) + 1
          pcthal(pp,mm) = pcthal(pp,mm) + 1

          costh = abs( dot_product( evtrcc(:,2), evtrcc0(:,2,ll) ) )
          pp = floor( costh / binw ) + 1
          pcthbe(pp,mm) = pcthbe(pp,mm) + 1

          costh = abs( dot_product( evtrcc(:,1), evtrcc0(:,1,ll) ) )
          pp = floor( costh / binw ) + 1
          pcthgm(pp,mm) = pcthgm(pp,mm) + 1

          pntcnt(mm) = pntcnt(mm) + 1
 
        end if

      end do
      end do
      end do
 
    end do

  end do
  end do
  end do
  write(*,*) 'defm-cij-alignpdf.x: point counted is: ', pointcnt * const

  open(26, file = 'defm-cij-alignpdf-'//trim(flnm))

    do mm = 1, npnt

      pcthal(:,mm) = pcthal(:,mm) / (pntcnt(mm) + mytiny) / binw
      pcthbe(:,mm) = pcthbe(:,mm) / (pntcnt(mm) + mytiny) / binw
      pcthgm(:,mm) = pcthgm(:,mm) / (pntcnt(mm) + mytiny) / binw

      write(26, '( "#index = ", I6)') mm
      do pp = 1, npdf
        write(26, '(15E15.3)') (pp-0.5_sp)*binw, pcthal(pp,mm), pcthbe(pp,mm), pcthgm(pp,mm)
      end do
      write(26, *)
      write(26, *)

    end do

  close(26)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(pcthal, pcthbe, pcthgm, pntcnt)
  
  write(*,*) 'defm-cij-alignpdf.x done.'

end program defmcijcor
