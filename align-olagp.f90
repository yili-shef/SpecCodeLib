program alignolagp
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 80
  real(sp), parameter :: evbnd = 1._sp, evbinw = evbnd/npnt
  real(sp), dimension(npnt) :: pdfal, pdfbe, pdfgm

  integer, parameter :: matz = 5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: evectors
  integer :: ierr

  integer, parameter, dimension(3,3) :: delij = (/1,0,0,0,1,0,0,0,1/)
  real(sp), dimension(3,3) :: aij

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33,p, wx, wy, wz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(sp) :: akk, meanakk, trc, ox, oy, oz
  real(dp) :: alpha, beta, gmma
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-olagp.x nx filelist prefix'
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
  allocate(  wx(lx1,ly,lz),  wy(lx1,ly,lz),  wz(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p12(lx1,ly,lz), p13(lx1,ly,lz) )
  allocate( p22(lx1,ly,lz), p23(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(28, file = 'align-olagp-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p33
    close(15)
    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p
    close(15)
    write(*,*) 'finishing reading ux uy uz and p'

    wx = eye * ( ky * p33 - kz * p22 )
    wy = eye * ( kz * p11 - kx * p33 )
    wz = eye * ( kx * p22 - ky * p11 )

    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

    ! dp/dx, dp/dy, dp/dz
    p11 = p * eye * kx ; p22 = p * eye * ky ; p33 = p * eye * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

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

    ! Trace of Cauchy-Green tensor
    p12 = cmplx( real(b11) * real(b11), aimag(b11) * aimag(b11) ) + &
          cmplx( real(b12) * real(b12), aimag(b12) * aimag(b12) ) + &
          cmplx( real(b13) * real(b13), aimag(b13) * aimag(b13) ) + &
          cmplx( real(b21) * real(b21), aimag(b21) * aimag(b21) ) + &
          cmplx( real(b22) * real(b22), aimag(b22) * aimag(b22) ) + &
          cmplx( real(b23) * real(b23), aimag(b23) * aimag(b23) ) + &
          cmplx( real(b31) * real(b31), aimag(b31) * aimag(b31) ) + &
          cmplx( real(b32) * real(b32), aimag(b32) * aimag(b32) ) + &
          cmplx( real(b33) * real(b33), aimag(b33) * aimag(b33) )

    trc = ( sum(real(p12(1:lx,:,:))) + sum(aimag(p12(1:lx,:,:))) ) * const


    ! dp/dX
    p12  = cmplx( real(p11,sp) * real(b11,sp), aimag(p11) * aimag(b11) ) + &
           cmplx( real(p22,sp) * real(b21,sp), aimag(p22) * aimag(b21) ) + &
           cmplx( real(p33,sp) * real(b31,sp), aimag(p33) * aimag(b31) )
    ! dp/dY
    p13  = cmplx( real(p11,sp) * real(b12,sp), aimag(p11) * aimag(b12) ) + &
           cmplx( real(p22,sp) * real(b22,sp), aimag(p22) * aimag(b22) ) + &
           cmplx( real(p33,sp) * real(b32,sp), aimag(p33) * aimag(b32) )
    ! dp/dZ
    p23  = cmplx( real(p11,sp) * real(b13,sp), aimag(p11) * aimag(b13) ) + &
           cmplx( real(p22,sp) * real(b23,sp), aimag(p22) * aimag(b23) ) + &
           cmplx( real(p33,sp) * real(b33,sp), aimag(p33) * aimag(b33) )

    !** p11 p22 p33 are released

    ! p12 = dp/dX, p13 = dp/dY, p23 = dp/dZ
    call rfftwnd_f77_one_real_to_complex(r2c3d,p12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,p13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,p23,ignore_me)

    ! p11 = d2p/dxdX, p22 = d2p/dydX, p33 = d2p/dzdX, p12 = dp/dX
    ! Note the '* const' part.
    p11 = p12 * eye * kx * const; p22 = p12 * eye * ky * const; p33 = p12 * eye * kz * const
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    ! d2p/dXdX
    b11 = cmplx( real(p11,sp) * real(b11,sp), aimag(p11) * aimag(b11) ) + &
          cmplx( real(p22,sp) * real(b21,sp), aimag(p22) * aimag(b21) ) + &
          cmplx( real(p33,sp) * real(b31,sp), aimag(p33) * aimag(b31) )
    !******* b11, b21, b31 are released, and b11 is re-used. ****

    ! b21 = d2p/dXdY afterwards. 
    b21 = cmplx( real(p11,sp) * real(b12,sp), aimag(p11) * aimag(b12) ) + &
          cmplx( real(p22,sp) * real(b22,sp), aimag(p22) * aimag(b22) ) + &
          cmplx( real(p33,sp) * real(b32,sp), aimag(p33) * aimag(b32) )
    !******* b21 re-used ****

    ! b31 = d2p/dXdZ afterwards. 
    b31 = cmplx( real(p11,sp) * real(b13,sp), aimag(p11) * aimag(b13) ) + &
          cmplx( real(p22,sp) * real(b23,sp), aimag(p22) * aimag(b23) ) + &
          cmplx( real(p33,sp) * real(b33,sp), aimag(p33) * aimag(b33) )

    ! p11 = d2p/dxdY, p22 = d2p/dydY, p33 = d2p/dzdY, p13 = dp/dY
    p11 = p13 * eye * kx * const; p22 = p13 * eye * ky * const; p33 = p13 * eye * kz * const
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    ! b22 = d2p/dYdY afterwards.
    b22 = cmplx( real(p11,sp) * real(b12,sp), aimag(p11) * aimag(b12) ) + &
          cmplx( real(p22,sp) * real(b22,sp), aimag(p22) * aimag(b22) ) + &
          cmplx( real(p33,sp) * real(b32,sp), aimag(p33) * aimag(b32) )
    !******* b12, b22, b32 are released, and b22 is re-used. ****

    ! b32 = d2p/dYdZ afterwards. 
    b32 = cmplx( real(p11,sp) * real(b13,sp), aimag(p11) * aimag(b13) ) + &
          cmplx( real(p22,sp) * real(b23,sp), aimag(p22) * aimag(b23) ) + &
          cmplx( real(p33,sp) * real(b33,sp), aimag(p33) * aimag(b33) )
    !******* b12, b22, b32 are released, and b32 is re-used ****

    ! p11 = d2p/dxdZ, p22 = d2p/dydZ, p33 = d2p/dzdZ, p23 = dp/dZ
    p11 = p23 * eye * kx * const; p22 = p23 * eye * ky * const; p33 = p23 * eye * kz * const
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    ! b33 = d2p/dZdZ afterwards
    b33 = cmplx( real(p11,sp) * real(b13,sp), aimag(p11) * aimag(b13) ) + &
          cmplx( real(p22,sp) * real(b23,sp), aimag(p22) * aimag(b23) ) + &
          cmplx( real(p33,sp) * real(b33,sp), aimag(p33) * aimag(b33) )

    ! Now:
    ! b11 = p11, b21 = p12, b31 = p13, b22 = p22, b32 = p23, b33 = p33
    ! and p11, p12, ..., p33 are all released. 

    ! Measure of Anisotropy:
    meanakk = 0._sp
    pdfal = 0._sp
    pdfbe = 0._sp
    pdfgm = 0._sp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          aij(1,1) = real( b11(ll,jj,kk),sp )
          aij(1,2) = real( b21(ll,jj,kk),sp )
          aij(1,3) = real( b31(ll,jj,kk),sp )
          aij(2,2) = real( b22(ll,jj,kk),sp )
          aij(2,3) = real( b32(ll,jj,kk),sp )
          aij(3,3) = real( b33(ll,jj,kk),sp )

          ox = real( wx(ll,jj,kk), sp )
          oy = real( wy(ll,jj,kk), sp )
          oz = real( wz(ll,jj,kk), sp )
      else
          ll = ii/2
          aij(1,1) = aimag( b11(ll,jj,kk) )
          aij(1,2) = aimag( b21(ll,jj,kk) )
          aij(1,3) = aimag( b31(ll,jj,kk) )
          aij(2,2) = aimag( b22(ll,jj,kk) )
          aij(2,3) = aimag( b32(ll,jj,kk) )
          aij(3,3) = aimag( b33(ll,jj,kk) )

          ox = aimag( wx(ll,jj,kk) )
          oy = aimag( wy(ll,jj,kk) )
          oz = aimag( wz(ll,jj,kk) )

      endif
      aij(2,1) = aij(1,2); aij(3,1) = aij(1,3); aij(3,2) = aij(2,3)
      ignore_me = sqrt(ox * ox + oy * oy + oz * oz )
      ox = ox / ignore_me; oy = oy / ignore_me; oz = oz / ignore_me
      
      akk = aij(1,1) + aij(2,2) + aij(3,3)
      aij = aij  - delij / 3._sp * akk

      meanakk = meanakk + akk

      call rs(3,3,aij,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
      alpha = abs( ox*evectors(1,3)+oy*evectors(2,3)+oz*evectors(3,3) )
      beta  = abs( ox*evectors(1,2)+oy*evectors(2,2)+oz*evectors(3,2) )
      gmma  = abs( ox*evectors(1,1)+oy*evectors(2,1)+oz*evectors(3,1) )
      
      ll=floor(alpha/evbinw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfal(ll)=pdfal(ll)+1
      ll=floor(beta/evbinw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfbe(ll)=pdfbe(ll)+1
      ll=floor(gmma/evbinw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfgm(ll)=pdfgm(ll)+1

    end do
    end do
    end do
    meanakk = meanakk * const

    pdfal = pdfal * const 
    pdfbe = pdfbe * const
    pdfgm = pdfgm * const
    write(*,*) 'normalization of pdfal', sum(pdfal)
    write(*,*) 'normalization of pdfbe', sum(pdfbe)
    write(*,*) 'normalization of pdfgm', sum(pdfgm)
    pdfal = pdfal / evbinw
    pdfbe = pdfbe / evbinw
    pdfgm = pdfgm / evbinw

    write(28, '( "#index = ", I6, "pKK=", E15.3, "TrC=", E15.3)') nfile, meanakk, trc
    do ll = 1, npnt
      write(28, '(15E15.3)') (ll - 0.5_sp) * evbinw, pdfal(ll), pdfbe(ll), pdfgm(ll)
    end do
    write(28, *)
    write(28, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p12,p13,p22,p23,p33,p,wx,wy,wz)
  
  call destroyplan3d

  write(*,*) 'align-olagp.x done.'

end program alignolagp
