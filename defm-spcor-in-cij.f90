program defmsijincij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  ! integer,  parameter :: npnt = 50
  ! real(sp), parameter :: bnd = 1._sp, binw = bnd/npnt
  ! real(sp), parameter :: bndphi = pi/2._sp, binwphi = bndphi/npnt

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize,evsize)  :: cc, sij, evtrcc, evtrss
  real(dp), dimension(evsize) :: evcc, evss

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
          write(*,*) ' Usage: ./defm-spcor-in-cij.x nx filelist prefix'
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

  open(27, file = 'defm-saa-in-cij-'//trim(flnm)//'.dat')
  write(27,'("# msaa, msaa1, msaa2, msaa3, vsaa, vsaa1, vsaa2, vsaa3")')
  open(26, file = 'defm-sbb-in-cij-'//trim(flnm)//'.dat')
  write(26,'("# msbb, msbb1, msbb2, msbb3, vsbb, vsbb1, vsbb2, vsbb3")')
  open(25, file = 'defm-sgg-in-cij-'//trim(flnm)//'.dat')
  write(25,'("# msgg, msgg1, msgg2, msgg3, vsgg, vsgg1, vsgg2, vsgg3")')
  open(24, file = 'defm-sab-in-cij-'//trim(flnm)//'.dat')
  write(24,'("# msab, msab1, msab2, msab3, vsab, vsab1, vsab2, vsab3")')
  open(23, file = 'defm-sag-in-cij-'//trim(flnm)//'.dat')
  write(23,'("# msag, msag1, msag2, msag3, vsag, vsag1, vsag2, vsag3")')
  open(22, file = 'defm-sbg-in-cij-'//trim(flnm)//'.dat')
  write(22,'("# msbg, msbg1, msbg2, msbg3, vsbg, vsbg1, vsbg2, vsbg3")')

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

    msaa = 0._dp; vsaa = 0._dp
    msbb = 0._dp; vsbb = 0._dp
    msgg = 0._dp; vsgg = 0._dp
    msab = 0._dp; vsab = 0._dp
    msag = 0._dp; vsag = 0._dp
    msbg = 0._dp; vsbg = 0._dp
    msaa1 = 0._dp; msaa2 = 0._dp; msaa3 = 0._dp
    vsaa1 = 0._dp; vsaa2 = 0._dp; vsaa3 = 0._dp
    msbb1 = 0._dp; msbb2 = 0._dp; msbb3 = 0._dp
    vsbb1 = 0._dp; vsbb2 = 0._dp; vsbb3 = 0._dp
    msgg1 = 0._dp; msgg2 = 0._dp; msgg3 = 0._dp
    vsgg1 = 0._dp; vsgg2 = 0._dp; vsgg3 = 0._dp
    msab1 = 0._dp; msab2 = 0._dp; msab3 = 0._dp
    vsab1 = 0._dp; vsab2 = 0._dp; vsab3 = 0._dp
    msag1 = 0._dp; msag2 = 0._dp; msag3 = 0._dp
    vsag1 = 0._dp; vsag2 = 0._dp; vsag3 = 0._dp
    msbg1 = 0._dp; msbg2 = 0._dp; msbg3 = 0._dp
    vsbg1 = 0._dp; vsbg2 = 0._dp; vsbg3 = 0._dp
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
     
      saa = dot_product(evtrcc(:,3), matmul(sij, evtrcc(:,3)))
      sbb = dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,2)))
      sgg = dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,1)))
      sab = dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,3)))
      sag = dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,3)))
      sbg = dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,2)))

      call dsyevr("V", "A", "U", evsize, sij, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evss, evtrss, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev sij. info = ', info
      ! The eigenvalues are in ascending order: smallest first.

      saa1 = evss(3) * (dot_product(evtrcc(:,3), evtrss(:,3)))**2 
      saa2 = evss(2) * (dot_product(evtrcc(:,3), evtrss(:,2)))**2 
      saa3 = evss(1) * (dot_product(evtrcc(:,3), evtrss(:,1)))**2 

      sbb1 = evss(3) * (dot_product(evtrcc(:,2), evtrss(:,3)))**2 
      sbb2 = evss(2) * (dot_product(evtrcc(:,2), evtrss(:,2)))**2 
      sbb3 = evss(1) * (dot_product(evtrcc(:,2), evtrss(:,1)))**2 

      sgg1 = evss(3) * (dot_product(evtrcc(:,1), evtrss(:,3)))**2 
      sgg2 = evss(2) * (dot_product(evtrcc(:,1), evtrss(:,2)))**2 
      sgg3 = evss(1) * (dot_product(evtrcc(:,1), evtrss(:,1)))**2 

      sab1 = evss(3) * dot_product(evtrcc(:,3), evtrss(:,3)) * dot_product(evtrcc(:,2), evtrss(:,3))
      sab2 = evss(2) * dot_product(evtrcc(:,3), evtrss(:,2)) * dot_product(evtrcc(:,2), evtrss(:,2))
      sab3 = evss(1) * dot_product(evtrcc(:,3), evtrss(:,1)) * dot_product(evtrcc(:,2), evtrss(:,1))

      sag1 = evss(3) * dot_product(evtrcc(:,3), evtrss(:,3)) * dot_product(evtrcc(:,1), evtrss(:,3))
      sag2 = evss(2) * dot_product(evtrcc(:,3), evtrss(:,2)) * dot_product(evtrcc(:,1), evtrss(:,2))
      sag3 = evss(1) * dot_product(evtrcc(:,3), evtrss(:,1)) * dot_product(evtrcc(:,1), evtrss(:,1))

      sbg1 = evss(3) * dot_product(evtrcc(:,2), evtrss(:,3)) * dot_product(evtrcc(:,1), evtrss(:,3))
      sbg2 = evss(2) * dot_product(evtrcc(:,2), evtrss(:,2)) * dot_product(evtrcc(:,1), evtrss(:,2))
      sbg3 = evss(1) * dot_product(evtrcc(:,2), evtrss(:,1)) * dot_product(evtrcc(:,1), evtrss(:,1))

      msaa = msaa + saa; vsaa = vsaa + saa * saa
      msbb = msbb + sbb; vsbb = vsbb + sbb * sbb
      msgg = msgg + sgg; vsgg = vsgg + sgg * sgg
      msab = msab + sab; vsab = vsab + sab * sab
      msag = msag + sag; vsag = vsag + sag * sag
      msbg = msbg + sbg; vsbg = vsbg + sbg * sbg

      msaa1 = msaa1 + saa1; vsaa1 = vsaa1 + saa1 * saa1
      msaa2 = msaa2 + saa2; vsaa2 = vsaa2 + saa2 * saa2
      msaa3 = msaa3 + saa3; vsaa3 = vsaa3 + saa3 * saa3

      msbb1 = msbb1 + sbb1; vsbb1 = vsbb1 + sbb1 * sbb1
      msbb2 = msbb2 + sbb2; vsbb2 = vsbb2 + sbb2 * sbb2
      msbb3 = msbb3 + sbb3; vsbb3 = vsbb3 + sbb3 * sbb3

      msgg1 = msgg1 + sgg1; vsgg1 = vsgg1 + sgg1 * sgg1
      msgg2 = msgg2 + sgg2; vsgg2 = vsgg2 + sgg2 * sgg2
      msgg3 = msgg3 + sgg3; vsgg3 = vsgg3 + sgg3 * sgg3

      msab1 = msab1 + sab1; vsab1 = vsab1 + sab1 * sab1
      msab2 = msab2 + sab2; vsab2 = vsab2 + sab2 * sab2
      msab3 = msab3 + sab3; vsab3 = vsab3 + sab3 * sab3

      msag1 = msag1 + sag1; vsag1 = vsag1 + sag1 * sag1
      msag2 = msag2 + sag2; vsag2 = vsag2 + sag2 * sag2
      msag3 = msag3 + sag3; vsag3 = vsag3 + sag3 * sag3

      msbg1 = msbg1 + sbg1; vsbg1 = vsbg1 + sbg1 * sbg1
      msbg2 = msbg2 + sbg2; vsbg2 = vsbg2 + sbg2 * sbg2
      msbg3 = msbg3 + sbg3; vsbg3 = vsbg3 + sbg3 * sbg3

    end do
    end do
    end do
    msaa = msaa * const; vsaa = vsaa * const
    msbb = msbb * const; vsbb = vsbb * const
    msgg = msgg * const; vsgg = vsgg * const
    msag = msag * const; vsag = vsag * const
    msbg = msbg * const; vsbg = vsbg * const
    msab = msab * const; vsab = vsab * const

    msaa1 = msaa1 * const; vsaa1 = vsaa1 * const
    msbb1 = msbb1 * const; vsbb1 = vsbb1 * const
    msgg1 = msgg1 * const; vsgg1 = vsgg1 * const
    msag1 = msag1 * const; vsag1 = vsag1 * const
    msbg1 = msbg1 * const; vsbg1 = vsbg1 * const
    msab1 = msab1 * const; vsab1 = vsab1 * const

    msaa2 = msaa2 * const; vsaa2 = vsaa2 * const
    msbb2 = msbb2 * const; vsbb2 = vsbb2 * const
    msgg2 = msgg2 * const; vsgg2 = vsgg2 * const
    msag2 = msag2 * const; vsag2 = vsag2 * const
    msbg2 = msbg2 * const; vsbg2 = vsbg2 * const
    msab2 = msab2 * const; vsab2 = vsab2 * const

    msaa3 = msaa3 * const; vsaa3 = vsaa3 * const
    msbb3 = msbb3 * const; vsbb3 = vsbb3 * const
    msgg3 = msgg3 * const; vsgg3 = vsgg3 * const
    msag3 = msag3 * const; vsag3 = vsag3 * const
    msbg3 = msbg3 * const; vsbg3 = vsbg3 * const
    msab3 = msab3 * const; vsab3 = vsab3 * const

    write(27, '(I4, 15E13.4)') nfile, msaa, msaa1, msaa2, msaa3, vsaa, vsaa1, vsaa2, vsaa3
    write(26, '(I4, 15E13.4)') nfile, msbb, msbb1, msbb2, msbb3, vsbb, vsbb1, vsbb2, vsbb3
    write(25, '(I4, 15E13.4)') nfile, msgg, msgg1, msgg2, msgg3, vsgg, vsgg1, vsgg2, vsgg3
    write(24, '(I4, 15E13.4)') nfile, msab, msab1, msab2, msab3, vsab, vsab1, vsab2, vsab3
    write(23, '(I4, 15E13.4)') nfile, msag, msag1, msag2, msag3, vsag, vsag1, vsag2, vsag3
    write(22, '(I4, 15E13.4)') nfile, msbg, msbg1, msbg2, msbg3, vsbg, vsbg1, vsbg2, vsbg3

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(26)
  close(25)
  close(24)
  close(23)
  close(22)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,s11,s22,s33,s12,s13,s23)
  
  call destroyplan3d

  write(*,*) 'defm-sij-in-cij.x done.'

end program defmsijincij
