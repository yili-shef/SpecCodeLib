program condpress
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, nn, ndel, nfile

  integer,  parameter :: npdf = 60
  real(sp), parameter :: bnd = 5., binw = 2.*bnd/npdf, nu = 0.0015

  real(dp), dimension(npdf) :: pdfaat, cndahp, cndahv, cnda2at
  real(dp) :: sijsij, sijsij0, dtmp, ahp, ahv, aat, a2at

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  
  complex(sp), allocatable, dimension(:,:,:) :: va11, va12, va13
  complex(sp), allocatable, dimension(:,:,:) :: va21, va22, va23
  complex(sp), allocatable, dimension(:,:,:) :: va31, va32, va33

  complex(sp), allocatable, dimension(:,:,:) :: p11, p12, p13
  complex(sp), allocatable, dimension(:,:,:) :: p22, p23, p33
  real(sp),    allocatable, dimension(:,:,:) :: g, k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(sp), dimension(3,3) :: aij, pij, vaij
  real(sp) :: tmp, ignore_me, const, delta_c
  complex(sp) :: tmpc
  character(90) :: str, strndel, fnm


  write(*,*) 
  write(*,'(''>>> Conditional terms in the Eq for Tr(A A^T) <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cnd-aat-dns.x nx dnsfilelist ndel filtertype '
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

  fnm = fnm(1:len_trim(fnm))

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1. / (nx * ny * nz)

  delta_c = ndel * 2 * pi / nx

  allocate( a11(lx1,ly,lz), a21(lx1,ly,lz), a31(lx1,ly,lz) )
  allocate( a12(lx1,ly,lz), a22(lx1,ly,lz), a32(lx1,ly,lz) )
  allocate( a13(lx1,ly,lz), a23(lx1,ly,lz), a33(lx1,ly,lz) )
  allocate( va11(lx1,ly,lz), va21(lx1,ly,lz), va31(lx1,ly,lz) )
  allocate( va12(lx1,ly,lz), va22(lx1,ly,lz), va32(lx1,ly,lz) )
  allocate( va13(lx1,ly,lz), va23(lx1,ly,lz), va33(lx1,ly,lz) )
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
  pdfaat = 0.d0
  cndahp = 0.d0
  cndahv = 0.d0
  cnda2at = 0.d0
  open(25, file = fnm(1:len_trim(fnm))//'.list' )

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

      a13 = a13 * g
      a23 = a23 * g
      a33 = a33 * g

      open(10,file='./out/p'//str( 1:len_trim(str) ),status='unknown',form='unformatted')
        read(10)p23
      close(10)

      p23 = p23 * g

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
        
        va11(ii,jj,kk) = -nu * a11(ii,jj,kk) * k2(ii,jj,kk)
        va12(ii,jj,kk) = -nu * a12(ii,jj,kk) * k2(ii,jj,kk)
        va13(ii,jj,kk) = -nu * a13(ii,jj,kk) * k2(ii,jj,kk)
        va21(ii,jj,kk) = -nu * a21(ii,jj,kk) * k2(ii,jj,kk)
        va22(ii,jj,kk) = -nu * a22(ii,jj,kk) * k2(ii,jj,kk)
        va23(ii,jj,kk) = -nu * a23(ii,jj,kk) * k2(ii,jj,kk)
        va31(ii,jj,kk) = -nu * a31(ii,jj,kk) * k2(ii,jj,kk)
        va32(ii,jj,kk) = -nu * a32(ii,jj,kk) * k2(ii,jj,kk)
        va33(ii,jj,kk) = -nu * a33(ii,jj,kk) * k2(ii,jj,kk)
        
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

      call rfftwnd_f77_one_complex_to_real(c2r3d,va11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va21,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va31,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va32,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,va33,ignore_me)

      ! Pressure Hessian  d^2p /dxi dxj Trace free
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

        ! Remove the trace
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
 
          vaij(1,1) = real( va11(ll,jj,kk) ); vaij(1,2) = real( va12(ll,jj,kk) ); vaij(1,3) = real( va13(ll,jj,kk) )
          vaij(2,1) = real( va21(ll,jj,kk) ); vaij(2,2) = real( va22(ll,jj,kk) ); vaij(2,3) = real( va23(ll,jj,kk) )
          vaij(3,1) = real( va31(ll,jj,kk) ); vaij(3,2) = real( va32(ll,jj,kk) ); vaij(3,3) = real( va33(ll,jj,kk) )
 
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
 
          vaij(1,1) = aimag( va11(ll,jj,kk) ); vaij(1,2) = aimag( va12(ll,jj,kk) ); vaij(1,3) = aimag( va13(ll,jj,kk) )
          vaij(2,1) = aimag( va21(ll,jj,kk) ); vaij(2,2) = aimag( va22(ll,jj,kk) ); vaij(2,3) = aimag( va23(ll,jj,kk) )
          vaij(3,1) = aimag( va31(ll,jj,kk) ); vaij(3,2) = aimag( va32(ll,jj,kk) ); vaij(3,3) = aimag( va33(ll,jj,kk) )
 
          pij(1,1) = aimag( p11(ll,jj,kk) )
          pij(1,2) = aimag( p12(ll,jj,kk) )
          pij(1,3) = aimag( p13(ll,jj,kk) )
          pij(2,2) = aimag( p22(ll,jj,kk) )
          pij(2,3) = aimag( p23(ll,jj,kk) )
          pij(3,3) = aimag( p33(ll,jj,kk) )
        end if
        pij(2,1) = pij(1,2); pij(3,1) = pij(1,3); pij(3,2) = pij(2,3)


        aat = sum(aij * aij) / sijsij0

        ahp = - 2 * sum(aij * pij) / sijsij0**1.5 ! -2 * Tr(A Hp)

        pij = vaij + transpose(vaij)
        ahv = sum( aij * pij ) / sijsij0**1.5  ! Tr(A Hv) + Tr(A Hv^T)

        pij = matmul(aij, aij)
        a2at = -2 * sum( pij * aij ) / sijsij0**1.5  ! -2 * Tr(A^2 A^T)

        ll = floor( (aat) / binw ) + 1  ! aat is positive definite
        if ( ll .ge. 1 .and. ll .le. npdf ) then 
          pdfaat(ll) = pdfaat(ll) + 1
          cndahp(ll) = cndahp(ll) + ahp
          cndahv(ll) = cndahv(ll) + ahv
          cnda2at(ll) = cnda2at(ll) + a2at
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

  cndahp = cndahp / pdfaat
  cndahv = cndahv / pdfaat
  cnda2at = cnda2at / pdfaat

  pdfaat = pdfaat * ignore_me 
  write(*,*) 'check pdfaat: ', sum(pdfaat)
  pdfaat = pdfaat / binw / binw

  ! output
  str='cnd-aph-dns-'//strndel(1:len_trim(strndel))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  open(15, file = str( 1:len_trim(str) ) )
    write(15,'(a)') 'Variables = "aat" "cndahp" "cndahv" "cnda2at" "pdfaat"'
    write(15,'(a, I4, a)') 'zone i=', npdf, ' f=point'
    dtmp = sijsij0 / sijsij
    do ii = 1, npdf
      write(15,'(8e18.5)') (ii-.5)*binw * dtmp, & 
                           cndahp(ii) * dtmp**3 * pdfaat(ii) / dtmp, &
                           cndahv(ii) * dtmp**3 * pdfaat(ii) / dtmp, &
                           cnda2at(ii) * dtmp**3 * pdfaat(ii) / dtmp, & 
                           pdfaat(ii) / dtmp 
    end do
  close(15) 

  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33)
  deallocate(va11,va12,va13,va21,va22,va23,va31,va32,va33)
  deallocate(p11,p12,p13,p22,p23,p33)
  deallocate(g,kx,ky,kz,k2)

  call destroyplan3d

  write(*,*) 'Done'

end program condpress
