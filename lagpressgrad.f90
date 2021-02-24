program lagpg
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 80
  real(sp), parameter :: bnd = 45._sp, binw = bnd/npnt
  real(sp), dimension(npnt) :: pdfafn, cndpkpk

  integer, parameter, dimension(3,3) :: delij = (/1,0,0,0,1,0,0,0,1/)

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p12,p13,p22,p23,p33,p, afn
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(sp) :: meanpX2, meanpY2, meanpZ2, meanafn, afntmp, pkpk, meanpkpk, dpdX, dpdY, dpdZ, trc
  real(sp) :: meaneulerpx
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./lagpressgrad.x nx filelist prefix'
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
  allocate( p11(lx1,ly,lz), p12(lx1,ly,lz), p13(lx1,ly,lz) )
  allocate( p22(lx1,ly,lz), p23(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p(lx1,ly,lz), afn(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'pkpk-'//trim(flnm)//'.dat')
  write(27,*) '# file pkpk pX2 pY2 pZ2 AA sqrt(trc) meaneulerpkpk'
  open(28, file = 'pkpkcndaa-'//trim(flnm)//'.dat')

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


    b11 = eye * kx * p11
    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
    afn = cmplx( real(b11,sp) * real(b11,sp), aimag(b11) * aimag(b11) )

    b11 = eye * ky * p22
    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
    afn = afn + cmplx( real(b11,sp) * real(b11,sp), aimag(b11) * aimag(b11) )

    b11 = eye * kz * p33
    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
    afn = afn + cmplx( real(b11,sp) * real(b11,sp), aimag(b11) * aimag(b11) )

    b12 = eye * ky * p11
    b21 = eye * kx * p22
    call rfftwnd_f77_one_complex_to_real(c2r3d,b12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b21,ignore_me)
    afn = afn + cmplx( real(b12,sp) * real(b12,sp), aimag(b12) * aimag(b12) ) &
              + cmplx( real(b21,sp) * real(b21,sp), aimag(b21) * aimag(b21) ) 

    b13 = eye * kz * p11
    b31 = eye * kx * p33
    call rfftwnd_f77_one_complex_to_real(c2r3d,b13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b31,ignore_me)
    afn = afn + cmplx( real(b13,sp) * real(b13,sp), aimag(b13) * aimag(b13) ) &
              + cmplx( real(b31,sp) * real(b31,sp), aimag(b31) * aimag(b31) )

    b23 = eye * kz * p22
    b32 = eye * ky * p33
    call rfftwnd_f77_one_complex_to_real(c2r3d,b23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b32,ignore_me)
    afn = afn + cmplx( real(b23,sp) * real(b23,sp), aimag(b23) * aimag(b23) ) &
              + cmplx( real(b32,sp) * real(b32,sp), aimag(b32) * aimag(b32) )


    ! dp/dx, dp/dy, dp/dz
    p11 = p * eye * kx ; p22 = p * eye * ky ; p33 = p * eye * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    p = cmplx( real(p11) * real(p11), aimag(p11) * aimag(p11) ) + &
        cmplx( real(p22) * real(p22), aimag(p22) * aimag(p22) ) + &
        cmplx( real(p33) * real(p33), aimag(p33) * aimag(p33) )

    meaneulerpx = sqrt( sum( real(p(1:lx,:,:)) + aimag(p(1:lx,:,:)) ) * const )


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
    trc = sqrt(trc)


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

    p12 = p12 / trc
    p13 = p13 / trc
    p23 = p23 / trc 

    meanafn = 0._sp; meanpkpk = 0._sp
    meanpX2 = 0._sp; meanpY2 = 0._sp; meanpZ2 = 0._sp
    pdfafn = 0._sp; cndpkpk = 0._sp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          dpdX = real( p12(ll,jj,kk),sp )
          dpdY = real( p13(ll,jj,kk),sp )
          dpdZ = real( p23(ll,jj,kk),sp )

          afntmp = real( afn(ll,jj,kk) )
      else
          ll = ii/2
          dpdX = aimag( p12(ll,jj,kk) )
          dpdY = aimag( p13(ll,jj,kk) )
          dpdZ = aimag( p23(ll,jj,kk) )

          afntmp = aimag( afn(ll,jj,kk) ) 
      endif
      
      meanpX2 = meanpX2 + dpdX * dpdX
      meanpY2 = meanpY2 + dpdY * dpdY
      meanpZ2 = meanpZ2 + dpdZ * dpdZ

      meanafn = meanafn + afntmp

      pkpk = dpdX * dpdX + dpdY * dpdY + dpdZ * dpdZ
      meanpkpk = meanpkpk + pkpk

      ll = floor( afntmp / binw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) then 
          pdfafn(ll) = pdfafn(ll) + 1
          cndpkpk(ll) = cndpkpk(ll) + pkpk
      end if

    end do
    end do
    end do
    meanpkpk = meanpkpk * const
    meanpX2 = meanpX2 * const
    meanpY2 = meanpY2 * const
    meanpZ2 = meanpZ2 * const
    meanafn = meanafn * const

    write(27, '(I5, 15E15.3)') nfile, meanpkpk, meanpX2, meanpY2, meanpZ2, meanafn, trc, meaneulerpx

    cndpkpk = cndpkpk / (pdfafn + mytiny)

    pdfafn = pdfafn * const 
    write(*,*) 'normalization of pdfafn', sum(pdfafn)
    pdfafn = pdfafn / binw

    write(28, *) '#', nfile
    do ll = 1, npnt
      write(28, '(15E15.3)') (ll - 0.5_sp) * binw, cndpkpk(ll), pdfafn(ll)
    end do
    write(28, *)
    write(28, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p12,p13,p22,p23,p33,p,afn)
  
  call destroyplan3d

  write(*,*) 'lagpressgrad.x done.'

end program lagpg
