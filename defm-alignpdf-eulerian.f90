program defmalignpdfeulerian
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,npnt,rr
  real(sp) :: const

  integer, parameter :: naverage = 200
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll, dx, xpp, ypp, zpp, xxrr, yyrr, zzrr
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
  integer :: pointcnt, noffset, forward, nlst

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  real(sp), allocatable, dimension(:) :: xp, yp, zp
  integer,  allocatable, dimension(:) :: pntcnt
  real(dp), allocatable, dimension(:,:) :: pcthal, pcthbe, pcthgm

  character(80) :: str, flnm, prf, str1, strll, mark, xyzlst
  
  if (iargc() .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-alignpdf-eulerian.x nx filelist prefix noffset forward'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     noffset: offset in file number in the list'
          write(*,*) '                     fw: 1 for forward, 0 for backward'
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

  call getarg(4,str)
  read(str, '(I20)') noffset

  call getarg(5,mark)
  read(mark, '(I20)') forward


  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  dx = 2*pi/nx

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), pcthal(npdf,npnt), pcthbe(npdf,npnt), pcthgm(npdf,npnt))
  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  allocate( xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz) )

  !-------------------------------------------------------------------------------------
  ! Will calculate the pdf for the "noffset" data file in each list only.
  ! For this we need to generate the file names instead of reading them from the lists. 
  ! Will use the smartrunner(.py) to run the code, so the interface is kept the same.

  xyzlst = flnm( len_trim(flnm) - 3 : len_trim(flnm) ) ! extract the starting file number
  read(xyzlst, "(I20)") nlst

  if (forward .eq. 0) then 
      ll = nlst - noffset
      nlst = nlst - 150 
      mark = 'bw'
  else
      ll = nlst + noffset
      nlst = nlst + 150
      mark = 'fw'
  end if

  write(str1, "(I20)") nlst ! get the last file number
  str1 = adjustl(str1)

  write(strll, "(I20)") ll  ! The file number in the list
  strll = adjustl(strll)

  xyzlst = trim(mark)//trim(xyzlst)//"to"//trim(str1)//'-'//trim(strll)//'.dat'
  ! file name for coordinates: xyzfwxxxxtoyyyy-zzzz.dat
  flnm = trim(flnm)//'-'//trim(strll)//'.dat'
  ! file name for bij: b11fwfromxxxx-zzzz.dat etc
  !--------------------------------------------------------------

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
  open(15,file='./out/xyz'//trim(xyzlst),form='unformatted')
    read(15) xp, yp, zp
  close(15)


  ll = -100
  do mm = 1, naverage

    ii = floor( ran1(ll) * nx ) + 1
    jj = floor( ran1(ll) * ny ) + 1
    kk = floor( ran1(ll) * nz ) + 1

    xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

    if ( mod(ii, 2) .eq. 1 ) then
      cc(1,1) = real( b11(ii/2+1,jj,kk),sp )
      cc(1,2) = real( b12(ii/2+1,jj,kk),sp )
      cc(1,3) = real( b13(ii/2+1,jj,kk),sp )
      cc(2,1) = real( b21(ii/2+1,jj,kk),sp )
      cc(2,2) = real( b22(ii/2+1,jj,kk),sp )
      cc(2,3) = real( b23(ii/2+1,jj,kk),sp )
      cc(3,1) = real( b31(ii/2+1,jj,kk),sp )
      cc(3,2) = real( b32(ii/2+1,jj,kk),sp )
      cc(3,3) = real( b33(ii/2+1,jj,kk),sp )
    else
      cc(1,1) = aimag( b11(ii/2,jj,kk) )
      cc(1,2) = aimag( b12(ii/2,jj,kk) )
      cc(1,3) = aimag( b13(ii/2,jj,kk) )
      cc(2,1) = aimag( b21(ii/2,jj,kk) )
      cc(2,2) = aimag( b22(ii/2,jj,kk) )
      cc(2,3) = aimag( b23(ii/2,jj,kk) )
      cc(3,1) = aimag( b31(ii/2,jj,kk) )
      cc(3,2) = aimag( b32(ii/2,jj,kk) )
      cc(3,3) = aimag( b33(ii/2,jj,kk) )
    end if
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

    pp = (kk-1) * nx * ny + (jj-1) * nx + ii
    xpp = xp(pp); ypp = yp(pp); zpp = zp(pp)
    xpp = modulo(xpp, 2*pi)
    ypp = modulo(ypp, 2*pi)
    zpp = modulo(zpp, 2*pi)

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

      rr = ( zz(ll) - 1 ) * nx * ny + ( yy(ll) - 1 ) * nx + xx(ll)
      xxrr = xp(rr); yyrr = yp(rr); zzrr = zp(rr)

      xxrr = modulo(xxrr, 2*pi)
      yyrr = modulo(yyrr, 2*pi)
      zzrr = modulo(zzrr, 2*pi)

      do kkloc = -1, 1 ! checking the images of the flow field
      do jjloc = -1, 1
      do iiloc = -1, 1

        xxll = xxrr + iiloc * 2 * pi
        yyll = yyrr + jjloc * 2 * pi
        zzll = zzrr + kkloc * 2 * pi

        mm = floor( sqrt( (xpp-xxll)**2 + (ypp-yyll)**2 + (zpp-zzll)**2 ) / dx + 0.5_sp )  + 1
 
        if ( mm .ge. 1 .and. mm .le. npnt) then
  
          costh = abs( dot_product( evtrcc(:,3), evtrcc0(:,3,ll) ) )
          rr = floor( costh / binw ) + 1
          pcthal(rr,mm) = pcthal(rr,mm) + 1

          costh = abs( dot_product( evtrcc(:,2), evtrcc0(:,2,ll) ) )
          rr = floor( costh / binw ) + 1
          pcthbe(rr,mm) = pcthbe(rr,mm) + 1

          costh = abs( dot_product( evtrcc(:,1), evtrcc0(:,1,ll) ) )
          rr = floor( costh / binw ) + 1
          pcthgm(rr,mm) = pcthgm(rr,mm) + 1

          pntcnt(mm) = pntcnt(mm) + 1
 
        end if

      end do
      end do
      end do
 
    end do

  end do
  end do
  end do
  write(*,*) 'defm-cij-alignpdf-eulerian.x: point counted is: ', pointcnt * const

  open(26, file = 'defm-alignpdf-eulerian-'//trim(flnm))

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
  deallocate(xp, yp, zp)
  
  write(*,*) 'defm-alignpdf-eulerian.x done.'

end program defmalignpdfeulerian
