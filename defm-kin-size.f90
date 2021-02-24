program defmkinsize
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 50
  real(sp), parameter :: bnd = 10._sp, binw = bnd/npnt

  real(sp), dimension(npnt) :: pdfcc, cndkin, cndmom
  real(sp), dimension(3,3)  :: cc, aa, bb

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: a11,a12,a13
  complex(sp), allocatable, dimension(:,:,:) :: a21,a22,a23
  complex(sp), allocatable, dimension(:,:,:) :: a31,a32,a33, trcc
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(sp) :: kin, ix, iy, iz, angularmom, meantrcc, meanaijaij, meankin, meanmom, ccii
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-kin-size.x nx filelist prefix'
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
  allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
  allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
  allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )
  allocate( trcc(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'kin-cnd-size-'//trim(flnm)//'.dat')
  write(27,'("# cc pdfcc cndkin cndmom")')
  open(28, file = 'mean-kin-mom-'//trim(flnm)//'.dat')
  write(28,'("# nfile meantrcc meankin meanmom meanaijaij")')

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

    a11 = eye * kx * b11; a12 = eye * ky * b11; a13 = eye * kz * b11
    a21 = eye * kx * b22; a22 = eye * ky * b22; a23 = eye * kz * b22
    a31 = eye * kx * b33; a32 = eye * ky * b33; a33 = eye * kz * b33 

    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

    ! mean aijaij
    trcc = cmplx( real(a11,sp) * real(a11,sp), aimag(a11) * aimag(a11) ) + &
           cmplx( real(a12,sp) * real(a12,sp), aimag(a12) * aimag(a12) ) + &
           cmplx( real(a13,sp) * real(a13,sp), aimag(a13) * aimag(a13) ) + &
           cmplx( real(a21,sp) * real(a21,sp), aimag(a21) * aimag(a21) ) + &
           cmplx( real(a22,sp) * real(a22,sp), aimag(a22) * aimag(a22) ) + &
           cmplx( real(a23,sp) * real(a23,sp), aimag(a23) * aimag(a23) ) + &
           cmplx( real(a31,sp) * real(a31,sp), aimag(a31) * aimag(a31) ) + &
           cmplx( real(a32,sp) * real(a32,sp), aimag(a32) * aimag(a32) ) + &
           cmplx( real(a33,sp) * real(a33,sp), aimag(a33) * aimag(a33) )
    meanaijaij = ( sum(real(trcc(1:lx,:,:),sp)) + sum(aimag(trcc(1:lx,:,:))) ) * const

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
    trcc = cmplx( real(b11,sp) * real(b11,sp), aimag(b11) * aimag(b11) ) + &
           cmplx( real(b12,sp) * real(b12,sp), aimag(b12) * aimag(b12) ) + &
           cmplx( real(b13,sp) * real(b13,sp), aimag(b13) * aimag(b13) ) + &
           cmplx( real(b21,sp) * real(b21,sp), aimag(b21) * aimag(b21) ) + &
           cmplx( real(b22,sp) * real(b22,sp), aimag(b22) * aimag(b22) ) + &
           cmplx( real(b23,sp) * real(b23,sp), aimag(b23) * aimag(b23) ) + &
           cmplx( real(b31,sp) * real(b31,sp), aimag(b31) * aimag(b31) ) + &
           cmplx( real(b32,sp) * real(b32,sp), aimag(b32) * aimag(b32) ) + &
           cmplx( real(b33,sp) * real(b33,sp), aimag(b33) * aimag(b33) )

    meantrcc = ( sum(real(trcc(1:lx,:,:),sp)) + sum(aimag(trcc(1:lx,:,:))) ) * const

    pdfcc = 0._sp
    cndkin = 0._sp
    cndmom = 0._sp
    meankin = 0._sp
    meanmom = 0._sp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          aa(1,1) = real( a11(ll,jj,kk),sp )
          aa(1,2) = real( a12(ll,jj,kk),sp )
          aa(1,3) = real( a13(ll,jj,kk),sp )
          aa(2,1) = real( a21(ll,jj,kk),sp )
          aa(2,2) = real( a22(ll,jj,kk),sp )
          aa(2,3) = real( a23(ll,jj,kk),sp )
          aa(3,1) = real( a31(ll,jj,kk),sp )
          aa(3,2) = real( a32(ll,jj,kk),sp )
          aa(3,3) = real( a33(ll,jj,kk),sp )

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

          aa(1,1) = aimag( a11(ll,jj,kk) )
          aa(1,2) = aimag( a12(ll,jj,kk) )
          aa(1,3) = aimag( a13(ll,jj,kk) )
          aa(2,1) = aimag( a21(ll,jj,kk) )
          aa(2,2) = aimag( a22(ll,jj,kk) )
          aa(2,3) = aimag( a23(ll,jj,kk) )
          aa(3,1) = aimag( a31(ll,jj,kk) )
          aa(3,2) = aimag( a32(ll,jj,kk) )
          aa(3,3) = aimag( a33(ll,jj,kk) )

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
      ! Kinetic energy
      bb = matmul(aa, cc)
      kin = sum(bb * bb) * .5_sp

      meankin = meankin + kin
      
      ! Angular momentum
      bb = matmul(bb, transpose(cc))
      ix = bb(2,3) - bb(3,2)
      iy = bb(3,1) - bb(1,3)
      iz = bb(1,2) - bb(2,1)
      angularmom = .5_sp * ( ix * ix + iy * iy + iz * iz)

      meanmom = meanmom + angularmom

      ! Shape of the volume
      cc = matmul( cc, transpose(cc) )
      ccii = cc(1,1) + cc(2,2) + cc(3,3)
      ccii = ccii / meantrcc

      ll = floor(ccii/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) then 
          cndkin(ll) = cndkin(ll) + kin
          cndmom(ll) = cndmom(ll) + angularmom
          pdfcc(ll) = pdfcc(ll)+1
      end if

    end do
    end do
    end do
    cndkin = cndkin / (pdfcc + mytiny)
    cndmom = cndmom / (pdfcc + mytiny)
    pdfcc = pdfcc * const / binw

    cndkin = cndkin/meanaijaij/meantrcc
    cndmom = cndmom/meanaijaij/meantrcc**2

    meankin = meankin * const
    meanmom = meanmom * const

    write(28, '(I6, 15E15.3)') nfile, meantrcc, meankin, meanmom, meanaijaij

    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(27, '(15E15.3)') (ll - 0.5_sp) * binw, pdfcc(ll), cndkin(ll), cndmom(ll)
    end do
    write(27, *)
    write(27, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33,trcc)
  
  call destroyplan3d

  write(*,*) 'defm-kin-size.x done.'

end program defmkinsize
