program condpress
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, nn, ndel, nfile

  integer,  parameter :: npdf = 60
  real(sp), parameter :: bnd = 5., binw = 2.*bnd/npdf
  real(dp), dimension(npdf, npdf) :: pdfqr, cndpq, cndpr
  real(sp), dimension(3,3) :: aij, pij, aij2

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  complex(sp), allocatable, dimension(:,:,:) :: p11, p12, p13
  complex(sp), allocatable, dimension(:,:,:) :: p22, p23, p33
  real(sp),    allocatable, dimension(:,:,:) :: g, k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(dp) :: sijsij, sijsij0
  real(sp) :: tmp, ignore_me, const, delta_c, q, r, pq, pr
  complex(sp) :: tmpc
  character(90) :: str, strndel, fnm


  write(*,*) 
  write(*,'(''>>> Conditional pressure Hessian on Q R plane <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cndpress.x nx dnsfilelist ndel filtertype '
          write(*,*) '        nx: resolution'
          write(*,*) '        dnsfilelist: dns data file list'
          write(*,*) '        ndel: ndel * dx = filter scale'
          write(*,*) '        filtertype: 0, Gaussian; 1, cutoff'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if 

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,strndel)
  read(strndel, '(I20)') ndel
  strndel = adjustl(strndel)

  ! filter type
  call getarg(4,str)
  read(str, '(I20)') mm

  fnm = fnm(1:len_trim(fnm))//'.list'

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1. / (nx * ny * nz)

  delta_c = ndel * 2 * pi / nx

  allocate( a11(lx1,ly,lz), a21(lx1,ly,lz), a31(lx1,ly,lz) )
  allocate( a12(lx1,ly,lz), a22(lx1,ly,lz), a32(lx1,ly,lz) )
  allocate( a13(lx1,ly,lz), a23(lx1,ly,lz), a33(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )
  allocate( g(lx1,ly,lz), kx(lx1), ky(ly), kz(lz), k2(lx1,ly,lz) )
  write(*,*) 'after allocation'

  call fftwplan3de(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call definefilter(mm,nx,delta_c, k2, g)
  write(*,*) 'after definefilter'


  sijsij = 0.d0
  pdfqr = 0.d0
  cndpq = 0.d0
  cndpr = 0.d0
  open(25, file = fnm(1:len_trim(fnm)) )

    nfile = 0
    do while ( .not. eof(25) )
      read(25,*) str
      write(*,*) str( 1:len_trim(str) )

      open(10,file='./out/ux'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a13
      close(10)
      open(10,file='./out/uy'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a23
      close(10)
      open(10,file='./out/uz'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)a33
      close(10)

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        a11(ii,jj,kk) = eye * a13(ii,jj,kk) * kx(ii)
        a12(ii,jj,kk) = eye * a13(ii,jj,kk) * ky(jj)
        a13(ii,jj,kk) = eye * a13(ii,jj,kk) * kz(kk)
        a21(ii,jj,kk) = eye * a23(ii,jj,kk) * kx(ii)
        a22(ii,jj,kk) = eye * a23(ii,jj,kk) * ky(jj)
        a23(ii,jj,kk) = eye * a23(ii,jj,kk) * kz(kk)
        a31(ii,jj,kk) = eye * a33(ii,jj,kk) * kx(ii)
        a32(ii,jj,kk) = eye * a33(ii,jj,kk) * ky(jj)
        a33(ii,jj,kk) = eye * a33(ii,jj,kk) * kz(kk)
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

      ! - Tr AA
      p23 = - cmplx( &
            real(a11)*real(a11) + real(a22)*real(a22) + real(a33)*real(a33) &
            + 2 * ( real(a12)*real(a21) + real(a13)*real(a31) + &
            real(a23)*real(a32) ), &
            aimag(a11)*aimag(a11) + aimag(a22)*aimag(a22) + aimag(a33)*aimag(a33) &
            + 2 * ( aimag(a12)*aimag(a21) + aimag(a13)*aimag(a31) + &
            aimag(a23)*aimag(a32) ) )

      call rfftwnd_f77_one_real_to_complex(r2c3d,p23,ignore_me)

      p23 = - p23 * const / k2       ! Pressure in Fourier space
      p23 = p23 * g                  ! Filtered pressure

      ! Pressure Hessian
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        p11(ii,jj,kk) = - kx(ii) * kx(ii) * p23(ii,jj,kk)
        p22(ii,jj,kk) = - ky(jj) * ky(jj) * p23(ii,jj,kk)
        p33(ii,jj,kk) = - kz(kk) * kz(kk) * p23(ii,jj,kk)
        p12(ii,jj,kk) = - kx(ii) * ky(jj) * p23(ii,jj,kk)
        p13(ii,jj,kk) = - kx(ii) * kz(kk) * p23(ii,jj,kk)
        p23(ii,jj,kk) = - ky(jj) * kz(kk) * p23(ii,jj,kk)

        tmpc = p11(ii,jj,kk) + p22(ii,jj,kk) + p33(ii,jj,kk)
        tmpc = - tmpc / 3.

        p11(ii,jj,kk) = p11(ii,jj,kk) + tmpc
        p22(ii,jj,kk) = p22(ii,jj,kk) + tmpc
        p33(ii,jj,kk) = p33(ii,jj,kk) + tmpc
      end do
      end do
      end do
                                                                                                          
      call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)


      if ( nfile .eq. 0 ) then

        ! Estimate the normalization factor

        sijsij0 = 0.
        do kk = 1, lz
        do jj = 1, ly
        do ii = 1, lx
          sijsij0 = sijsij0 + real(a11(ii,jj,kk))**2 + aimag(a11(ii,jj,kk))**2
          sijsij0 = sijsij0 + real(a22(ii,jj,kk))**2 + aimag(a22(ii,jj,kk))**2
          sijsij0 = sijsij0 + real(a33(ii,jj,kk))**2 + aimag(a33(ii,jj,kk))**2

          tmp = .5 * real( ( a12(ii,jj,kk) + a21(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp
          tmp = .5 * aimag( ( a12(ii,jj,kk) + a21(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp

          tmp = .5 * real( ( a13(ii,jj,kk) + a31(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp
          tmp = .5 * aimag( ( a13(ii,jj,kk) + a31(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp

          tmp = .5 * real( ( a23(ii,jj,kk) + a32(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp
          tmp = .5 * aimag( ( a23(ii,jj,kk) + a32(ii,jj,kk) ) )
          sijsij0 = sijsij0 + 2. * tmp * tmp
        end do
        end do
        end do
        sijsij0 = const * sijsij0

      end if

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, nx 

        if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1) / 2

          aij(1,1) = real( a11(ll,jj,kk) ); aij(1,2) = real( a12(ll,jj,kk) ); aij(1,3) = real( a13(ll,jj,kk) )
          aij(2,1) = real( a21(ll,jj,kk) ); aij(2,2) = real( a22(ll,jj,kk) ); aij(2,3) = real( a23(ll,jj,kk) )
          aij(3,1) = real( a31(ll,jj,kk) ); aij(3,2) = real( a32(ll,jj,kk) ); aij(3,3) = real( a33(ll,jj,kk) )
 
          pij(1,1) = real( p11(ll,jj,kk) )
          pij(1,2) = real( p12(ll,jj,kk) )
          pij(1,3) = real( p13(ll,jj,kk) )
          pij(2,2) = real( p22(ll,jj,kk) )
          pij(2,3) = real( p23(ll,jj,kk) )
          pij(3,3) = real( p33(ll,jj,kk) )
        else 
          ll = ii / 2

          aij(1,1) = aimag( a11(ll,jj,kk) ); aij(1,2) = aimag( a12(ll,jj,kk) ); aij(1,3) = aimag( a13(ll,jj,kk) )
          aij(2,1) = aimag( a21(ll,jj,kk) ); aij(2,2) = aimag( a22(ll,jj,kk) ); aij(2,3) = aimag( a23(ll,jj,kk) )
          aij(3,1) = aimag( a31(ll,jj,kk) ); aij(3,2) = aimag( a32(ll,jj,kk) ); aij(3,3) = aimag( a33(ll,jj,kk) )
 
          pij(1,1) = aimag( p11(ll,jj,kk) )
          pij(1,2) = aimag( p12(ll,jj,kk) )
          pij(1,3) = aimag( p13(ll,jj,kk) )
          pij(2,2) = aimag( p22(ll,jj,kk) )
          pij(2,3) = aimag( p23(ll,jj,kk) )
          pij(3,3) = aimag( p33(ll,jj,kk) )
        end if
        pij(2,1) = pij(1,2); pij(3,1) = pij(1,3); pij(3,2) = pij(2,3)


        aij2 = matmul(aij, pij)
        pq = aij2(1,1) + aij2(2,2) + aij2(3,3)
        pq = pq / sijsij0**1.5

        aij2 = matmul(aij, aij)
        q = aij2(1,1) + aij2(2,2) + aij2(3,3)
        q = -.5 * q / sijsij0 

        pij = matmul(aij2, pij)
        pr = pij(1,1) + pij(2,2) + pij(3,3)
        pr = pr / sijsij0**2

        aij2 = matmul(aij2, aij)
        r = aij2(1,1) + aij2(2,2) + aij2(3,3)
        r = - r / sijsij0**1.5 / 3


        ll = floor( (q + bnd) / binw ) + 1
        mm = floor( (r + bnd) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npdf .and. mm .ge. 1 .and. mm .le. npdf ) then 
          pdfqr(ll,mm) = pdfqr(ll,mm) + 1
          cndpq(ll,mm) = cndpq(ll,mm) + pq
          cndpr(ll,mm) = cndpr(ll,mm) + pr
        end if

        pij = .5 * (aij + transpose(aij))
        pij = matmul(pij, pij)
        sijsij = sijsij + pij(1,1) + pij(2,2) + pij(3,3)

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(25)

  ignore_me = const / nfile

  sijsij = sijsij * ignore_me

  cndpq = cndpq / pdfqr
  cndpr = cndpr / pdfqr

  pdfqr = pdfqr * ignore_me 
  write(*,*) 'check pdfqr: ', sum(pdfqr)
  pdfqr = pdfqr / binw / binw


  ! output
  str='cnd-p-hessian-dns-'//strndel(1:len_trim(strndel))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  open(15, file = str( 1:len_trim(str) ) )
    tmp = sijsij0 / sijsij
    do jj = 1, npdf
    do ii = 1, npdf
      write(15,'(8e18.5)') (-bnd + (ii-.5)*binw ) * tmp,  (-bnd + (ii-.5)*binw ) * tmp, &
                           cndpq(ii,jj) * tmp**3 * pdfqr(ii,jj) / tmp**2, &
                           cndpr(ii,jj) * tmp**4 * pdfqr(ii,jj) / tmp**2
    end do
    end do
  close(15) 

  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(g,kx,ky,kz,k2)

  call destroyplan3d

  write(*,*) 'Done'

end program condpress
