program condpress
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, nn, ndel, nfile

  real(sp), parameter :: beta = -.1, tau = 0.1

  integer,  parameter :: npdf = 60
  real(sp), parameter :: bnd = 5., binw = 2.*bnd/npdf
  real(dp), dimension(npdf, npdf) :: pdfqr, cndpq, cndpr
  real(sp), dimension(3,3) :: aij, aij2, aij3, aijt, hpij

  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  real(sp),    allocatable, dimension(:,:,:) :: g, k2
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  real(dp) :: sijsij, sijsij0
  real(sp) :: tmp, ignore_me, const, delta_c, q, r, pq, pr, traa, traat
  character(90) :: str, strndel, fnm

  real(sp) :: discrim, alpha


  write(*,*) 
  write(*,'(''>>> Conditional pressure Hessian on Q R plane: model <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cndpress-cml.x nx dnsfilelist ndel filtertype '
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
 
        else 
          ll = ii / 2

          aij(1,1) = aimag( a11(ll,jj,kk) ); aij(1,2) = aimag( a12(ll,jj,kk) ); aij(1,3) = aimag( a13(ll,jj,kk) )
          aij(2,1) = aimag( a21(ll,jj,kk) ); aij(2,2) = aimag( a22(ll,jj,kk) ); aij(2,3) = aimag( a23(ll,jj,kk) )
          aij(3,1) = aimag( a31(ll,jj,kk) ); aij(3,2) = aimag( a32(ll,jj,kk) ); aij(3,3) = aimag( a33(ll,jj,kk) )
 
        end if


        ! Q and R
        aij2 = matmul(aij, aij)
        q = aij2(1,1) + aij2(2,2) + aij2(3,3)
        q = -.5 * q / sijsij0 

        aij3 = matmul(aij2, aij)
        r = aij3(1,1) + aij3(2,2) + aij3(3,3)
        r = - r / sijsij0**1.5 / 3

        ! Projection pressure Hessian
        ! discrim = ( (27./4.) * r*r + q*q*q ) / ( sum(aij * aij) / sijsij0 )**3
        ! discrim =  .25 * sum( (aij + transpose(aij)) * (aij + transpose(aij)) ) / sum(aij * aij) 
        ! alpha = beta * discrim
        alpha = beta

        aijt = transpose(aij)
        traa = aij2(1,1) + aij2(2,2) + aij2(3,3)

        hpij = & ! -alpha * (matmul(aijt, aij) - matmul(aij, aijt)) & 
               + tau * (traa / 3) * (aij + aijt)  &
               + alpha * tau * ( matmul( matmul(aijt, aijt), aij ) - matmul( matmul(aijt, aij), aijt )  &
                              + matmul( matmul(aijt, aij), aij ) - matmul( matmul(aij, aijt), aij ) )

        aij2 = matmul(aij, hpij) 
        pq = ( aij2(1,1) + aij2(2,2) + aij2(3,3) ) / sijsij0**1.5

        aij3 = matmul(aij, aij2)
        pr = ( aij3(1,1) + aij3(2,2) + aij3(3,3) ) / sijsij0**2


        ! PDF and conditional averages
        ll = floor( (q + bnd) / binw ) + 1
        mm = floor( (r + bnd) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npdf .and. mm .ge. 1 .and. mm .le. npdf ) then 
          pdfqr(ll,mm) = pdfqr(ll,mm) + 1
          cndpq(ll,mm) = cndpq(ll,mm) + pq
          cndpr(ll,mm) = cndpr(ll,mm) + pr
        end if

        aij2 = .5 * (aij + transpose(aij))
        sijsij = sijsij + sum(aij2 * aij2)

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(25)

  ignore_me = const / nfile

  sijsij = sijsij * ignore_me

  cndpq = cndpq / (pdfqr + tiny)
  cndpr = cndpr / (pdfqr + tiny)

  pdfqr = pdfqr * ignore_me 
  write(*,*) 'check pdfqr: ', sum(pdfqr)
  pdfqr = pdfqr / binw / binw


  ! output
  str='p-hessian-cml-'//strndel(1:len_trim(strndel))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  open(15, file = str( 1:len_trim(str) ) )
    write(15,'(a)') 'Variables = "q" "r" "cndpq" "cndpr" "pdfqr"'
    write(15,'(a, I4, a, I4, a)') 'zone i=', npdf, ' j=', npdf, ' f=point'
    tmp = sijsij0 / sijsij
    do ii = 1, npdf
    do jj = 1, npdf
      write(15,'(8e18.5)') (-bnd + (jj-.5)*binw ) * tmp,  (-bnd + (ii-.5)*binw ) * tmp**1.5, &
                           cndpr(ii,jj) * tmp**4 * pdfqr(ii,jj) / tmp**2.5, &
                           cndpq(ii,jj) * tmp**3 * pdfqr(ii,jj) / tmp**2.5, &
                           pdfqr(ii,jj) / tmp**2
    end do
    end do
  close(15) 

  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33)
  deallocate(g,kx,ky,kz,k2)

  call destroyplan3d

  write(*,*) 'Done'

end program condpress
