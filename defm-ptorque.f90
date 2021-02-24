program defmptorque
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 120
  real(sp), parameter :: bnd = 15._sp, binw = 2._sp*bnd/npnt
  real(dp), dimension(npnt) :: pdfpab, pdfpag, pdfpbg

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  real(dp), dimension(evsize,evsize)  :: cc, pij, evtrcc
  real(dp), dimension(evsize) :: evcc
  ! ----------------------------------------------------------------------------------
  real(dp) :: mta, mtb, mtg, vta, vtb, vtg, pab, pag, pbg, torque_a, torque_b, torque_g
  real(dp) :: mpab, mpag, mpbg, vpab, vpag, vpbg
  real(dp), parameter :: ugly = 8._dp
  integer :: pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-ptorque.x nx filelist prefix'
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

  open(27, file = 'defm-ptorque-mean-torque-'//trim(flnm)//'.dat')
  open(26, file = 'defm-ptorque-mean-pij-'//trim(flnm)//'.dat')
  open(25, file = 'defm-ptoruqe-pdfpij-'//trim(flnm)//'.dat')

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


    mta = 0._dp; mtb = 0._dp; mtg = 0._dp
    vta = 0._dp; vtb = 0._dp; vtg = 0._dp
    mpab = 0._dp; vpab = 0._dp
    mpag = 0._dp; vpag = 0._dp
    mpbg = 0._dp; vpbg = 0._dp
    pdfpab = 0._dp; pdfpag = 0._dp; pdfpbg = 0._dp
    pointcnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          pij(1,1) = real( p11(ll,jj,kk),sp )
          pij(1,2) = real( p12(ll,jj,kk),sp )
          pij(1,3) = real( p13(ll,jj,kk),sp )
          pij(2,2) = real( p22(ll,jj,kk),sp )
          pij(2,3) = real( p23(ll,jj,kk),sp )
          pij(3,3) = real( p33(ll,jj,kk),sp )

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
          pij(1,1) = aimag( p11(ll,jj,kk) )
          pij(1,2) = aimag( p12(ll,jj,kk) )
          pij(1,3) = aimag( p13(ll,jj,kk) )
          pij(2,2) = aimag( p22(ll,jj,kk) )
          pij(2,3) = aimag( p23(ll,jj,kk) )
          pij(3,3) = aimag( p33(ll,jj,kk) )

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
      pij(2,1) = pij(1,2); pij(3,1) = pij(1,3); pij(3,2) = pij(2,3)

      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      if (      abs(evcc(1) - evcc(2)) .le. ugly * dlamch('e') * abs(evcc(2)) &
           .or. abs(evcc(2) - evcc(3)) .le. ugly * dlamch('e') * abs(evcc(3)) &
           .or. abs(evcc(3) - evcc(1)) .le. ugly * dlamch('e') * abs(evcc(3)) ) then 
          cycle
      else
          pointcnt = pointcnt + 1
      end if

      ! The eigenvalues are in ascending order: smallest first.

      ! paa = dot_product(evtrcc(3), matmul(pij, evtrcc(3)))
      ! pbb = dot_product(evtrcc(2), matmul(pij, evtrcc(2)))
      ! pgg = dot_product(evtrcc(1), matmul(pij, evtrcc(1)))

      pab = dot_product(evtrcc(:,2), matmul(pij, evtrcc(:,3)))
      pag = dot_product(evtrcc(:,1), matmul(pij, evtrcc(:,3)))
      pbg = dot_product(evtrcc(:,1), matmul(pij, evtrcc(:,2)))

      mpab = mpab + pab
      mpag = mpag + pag
      mpbg = mpbg + pbg

      vpab = vpab + pab * pab
      vpag = vpag + pag * pag
      vpbg = vpbg + pbg * pbg

      ll = floor( (pab + bnd) / binw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfpab(ll) = pdfpab(ll) + 1
      ll = floor( (pag + bnd) / binw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfpag(ll) = pdfpag(ll) + 1
      ll = floor( (pbg + bnd) / binw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfpbg(ll) = pdfpbg(ll) + 1
     

      torque_a = (evcc(1) - evcc(2)) * pbg
      torque_b = (evcc(3) - evcc(1)) * pag
      torque_g = (evcc(2) - evcc(3)) * pab

      mta = mta + torque_a
      mtb = mtb + torque_b
      mtg = mtg + torque_g

      vta = vta + torque_a * torque_a
      vtb = vtb + torque_b * torque_b
      vtg = vtg + torque_g * torque_g
     
    end do
    end do
    end do
    mta = mta * const
    mtb = mtb * const
    mtg = mtg * const
    vta = vta * const
    vtb = vtb * const
    vtg = vtg * const

    mpab = mpab * const
    mpag = mpag * const
    mpbg = mpbg * const
    vpab = vpab * const
    vpag = vpag * const
    vpbg = vpbg * const

    pdfpab = pdfpab * const
    pdfpag = pdfpag * const
    pdfpbg = pdfpbg * const

    write(*,*) 'normalization pdfpab: ', sum(pdfpab)
    write(*,*) 'normalization pdfpag: ', sum(pdfpag)
    write(*,*) 'normalization pdfpbg: ', sum(pdfpbg)
    pdfpab = pdfpab / binw
    pdfpag = pdfpag / binw
    pdfpbg = pdfpbg / binw

    write(*,*) 'Points counted: ', pointcnt * const

    write(27, '(I4,15E15.3)') nfile,  mta,  mtb,  mtg,  vta,  vtb,  vtg
    write(26, '(I4,15E15.3)') nfile, mpab, mpag, mpbg, vpab, vpag, vpbg

    write(25, '("# index = ", I6)') nfile
    do ll = 1, npnt
      write(25, '(15E13.4)') -bnd + (ll -.5_sp) * binw, pdfpab(ll), pdfpag(ll), pdfpbg(ll)
    end do
    write(25,*)
    write(25,*)


    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(26)
  close(25)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,p12,p13,p23)
  
  call destroyplan3d

  write(*,*) 'defm-ptorque.x done.'

end program defmptorque
