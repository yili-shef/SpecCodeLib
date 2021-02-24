program align3dcijsij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer  :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nn
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt=25
  real(sp), parameter :: bwcth=1._sp/npnt, bwphi=pi/2._sp/npnt, bwzeta=pi/2._sp/npnt
  real(sp), dimension(npnt,npnt,npnt) :: p3d

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s22,s33,s12,s13,s23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf
  real(dp) :: ctheta, cos1, cos2, cosphi, phi, cos3, cos4, coszeta, zeta

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), evcc(evsize), evss(evsize), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  real(dp), dimension(evsize,evsize)  :: cc, sij, evtrcc, evtrss
  ! ----------------------------------------------------------------------------------

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3dcs.x nx file# prefix'
          write(*,*) '                       nx: resolution of data'
          write(*,*) '                       file#: number of data file'
          write(*,*) '                       prefix: prefix for datafiles for bij'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file number 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! prefix
  call getarg(3,prf)
  prf = adjustl(prf)

  ny = nx; nz = nx
  lx = nx/2; lx1 = nx/2 + 1; ly = ny; lz = nz
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

  open(15,file='./out/ux'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b21
  close(15)
  open(15,file='./out/uy'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b31
  close(15)
  open(15,file='./out/uz'//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b32
  close(15)
  write(*,*) 'finishing reading velocity and pressure'

  s11 = eye * kx * b21
  s22 = eye * ky * b31
  s33 = eye * kz * b32
  s12 = eye * ( ky * b21 + kx * b31 ) * .5_sp
  s13 = eye * ( kz * b21 + kx * b32 ) * .5_sp
  s23 = eye * ( ky * b32 + kz * b31 ) * .5_sp

  call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 

  open(15,file='./out/b11'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b11
  close(15)
  open(15,file='./out/b12'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b12
  close(15)
  open(15,file='./out/b13'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b13
  close(15)
  open(15,file='./out/b21'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b21
  close(15)
  open(15,file='./out/b22'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b22
  close(15)
  open(15,file='./out/b23'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b23
  close(15)
  open(15,file='./out/b31'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b31
  close(15)
  open(15,file='./out/b32'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b32
  close(15)
  open(15,file='./out/b33'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat',form='unformatted')
    read(15) b33
  close(15)
  write(*,*) 'finishing reading bij'


  p3d = 0._sp
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
    !call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
    !lwork, iwork, liwork, info)

    call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
    if ( .not. (info .eq. 0 ) ) write(*,*) 'Something wrong with ev cc. info = ', info

    call dsyevr("V", "A", "U", evsize, sij, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                nfound, evss, evtrss, evsize, isuppz, work, lwork, iwork, liwork, info)
    if ( .not. (info .eq. 0 ) ) write(*,*) 'Something wrong with ev ss. info = ', info
    ! The eigenvalues are in ascending order: smallest first.

    ctheta = abs( sum( evtrcc(:,3) * evtrss(:,3) ) )
    cos1 = sum( evtrss(:,2) * evtrcc(:,3) )
    cos2 = sum( evtrss(:,1) * evtrcc(:,3) )
    cosphi = cos1 / sqrt( cos1 * cos1 + cos2 * cos2 )
    phi = acos( abs( cosphi ) )
    cos3 = sum( evtrcc(:,2) * evtrss(:,1) )
    cos4 = sum( evtrcc(:,1) * evtrss(:,1) )
    coszeta = cos4 / sqrt( cos3 * cos3 + cos4 * cos4 )
    zeta = acos( abs( coszeta ) )

    ll = floor(ctheta/bwcth)+1
    mm = floor(phi/bwphi)+1
    nn = floor(zeta/bwzeta)+1

    p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1
  
  end do
  end do
  end do

  p3d = p3d * const / bwcth / bwphi / bwzeta

  open(15, file = 'p3d-align3dcijsij-'//prf(1:len_trim(prf))//flnm(1:len_trim(flnm))//'.dat')
    write(15,'(''#zone t="3d align pdf", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') npnt,npnt,npnt
    ! Note the order of looping must be the same as that in 3dscalarpdf2vtk. Any change must be
    ! applied to both.
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
      p3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,s11,s22,s33,s12,s13,s23)
  call destroyplan3d

  write(*,*) 'align-3dcs.x done.'

end program align3dcijsij
