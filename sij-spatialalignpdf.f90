program sijspatialalign
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,nfile,npnt
  real(sp) :: const, ignore_me

  integer, parameter :: naverage = 300
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

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp), allocatable, dimension(:,:,:) :: kx, ky, kz
  integer,  allocatable, dimension(:)     :: pntcnt
  real(dp), allocatable, dimension(:,:) :: pcthal, pcthbe, pcthgm

  character(80) :: str, flnm, fpath
  
  if (iargc() .ne. 2) then
    write(*,*)
    write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
    write(*,*) 
    write(*,*) ' Usage: ./sij-spatialalignpdf.x nx filelist'
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
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), pcthal(npdf,npnt), pcthbe(npdf,npnt), pcthgm(npdf,npnt))
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( s33(lx1,ly,lz) )
  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,lx1,ly,lz)

  open(20, file = flnm(1 : len_trim(flnm))//'.list')

  pcthal = 0._dp; pcthbe = 0._dp; pcthgm = 0._dp
  pntcnt = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) "sij-spatialalignpdf.x: ", fpath(1 : len_trim(fpath))

    open(15,file='./out/ux'//trim(fpath),form='unformatted')
      read(15) s11
    close(15)
    open(15,file='./out/uy'//trim(fpath),form='unformatted')
      read(15) s22
    close(15)
    open(15,file='./out/uz'//trim(fpath),form='unformatted')
      read(15) s33
    close(15)
  
    s12=.5_sp*eye*(kx*s22+ky*s11)
    s13=.5_sp*eye*(kx*s33+kz*s11)
    s23=.5_sp*eye*(ky*s33+kz*s22)
    s11=eye*kx*s11
    s22=eye*ky*s22
    s33=eye*kz*s33
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
    ll = -100
    do mm = 1, naverage
 
      ii = floor( ran1(ll) * nx ) + 1
      jj = floor( ran1(ll) * ny ) + 1
      kk = floor( ran1(ll) * nz ) + 1
 
      xx(mm) = ii; yy(mm) = jj; zz(mm) = kk
 
      if ( mod(ii, 2) .eq. 1 ) then
        cc(1,1) = real( s11(ii/2+1,jj,kk),sp )
        cc(1,2) = real( s12(ii/2+1,jj,kk),sp )
        cc(1,3) = real( s13(ii/2+1,jj,kk),sp )
        cc(2,2) = real( s22(ii/2+1,jj,kk),sp )
        cc(2,3) = real( s23(ii/2+1,jj,kk),sp )
        cc(3,3) = real( s33(ii/2+1,jj,kk),sp )
      else 
        cc(1,1) = aimag( s11(ii/2,jj,kk) )
        cc(1,2) = aimag( s12(ii/2,jj,kk) )
        cc(1,3) = aimag( s13(ii/2,jj,kk) )
        cc(2,2) = aimag( s22(ii/2,jj,kk) )
        cc(2,3) = aimag( s23(ii/2,jj,kk) )
        cc(3,3) = aimag( s33(ii/2,jj,kk) )
      end if
      cc(2,1) = cc(1,2); cc(3,1) = cc(1,3); cc(3,2) = cc(2,3)
 
      ! The eigenvalues are in ascending order: smallest first.
      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc0(:,mm), evtrcc0(:,:,mm), evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info
 
    end do
 
 
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          cc(1,1) = real( s11(ll,jj,kk),sp )
          cc(1,2) = real( s12(ll,jj,kk),sp )
          cc(1,3) = real( s13(ll,jj,kk),sp )
          cc(2,1) = real( s12(ll,jj,kk),sp )
          cc(2,2) = real( s22(ll,jj,kk),sp )
          cc(2,3) = real( s23(ll,jj,kk),sp )
          cc(3,1) = real( s13(ll,jj,kk),sp )
          cc(3,2) = real( s23(ll,jj,kk),sp )
          cc(3,3) = real( s33(ll,jj,kk),sp )
      else
          ll = ii/2
 
          cc(1,1) = aimag( s11(ll,jj,kk) )
          cc(1,2) = aimag( s12(ll,jj,kk) )
          cc(1,3) = aimag( s13(ll,jj,kk) )
          cc(2,1) = aimag( s12(ll,jj,kk) )
          cc(2,2) = aimag( s22(ll,jj,kk) )
          cc(2,3) = aimag( s23(ll,jj,kk) )
          cc(3,1) = aimag( s13(ll,jj,kk) )
          cc(3,2) = aimag( s23(ll,jj,kk) )
          cc(3,3) = aimag( s33(ll,jj,kk) )
      endif
 
 
      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info
      ! The eigenvalues are in ascending order: smallest first.
 
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

  end do
  close(20)

  open(26, file = 'sij-spatialalignpdf-'//trim(flnm)//'.dat')

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

  call destroyplan3d

  deallocate(s11,s12,s13,s22,s23,s33,kx,ky,kz)
  deallocate(pcthal, pcthbe, pcthgm, pntcnt)
  
  write(*,*) 'sij-spatialalignpdf.x done.'

end program sijspatialalign
