program alignlagpc
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 50
  real(sp), parameter :: bnd = 1._sp, binw = bnd/npnt
  real(sp), dimension(npnt) :: pdfal, pdfbe, pdfgm
  real(sp), dimension(3,3)  :: cc

  integer, parameter :: matz = 5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: evectors
  integer :: ierr

  integer, parameter, dimension(3,3) :: delij = (/1,0,0,0,1,0,0,0,1/)

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(sp) :: afntmp, dpdx, dpdy, dpdz
  real(dp) :: alpha, beta, gmma
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-lagpc.x nx filelist prefix'
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

  open(28, file = 'align-lagpc-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) p12
    close(15)
    write(*,*) 'finishing reading p'

    ! dp/dx, dp/dy, dp/dz
    p11 = p12 * eye * kx ; p22 = p12 * eye * ky ; p33 = p12 * eye * kz
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


    pdfal = 0._sp; pdfbe = 0._sp; pdfgm = 0._sp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          dpdx = real( p12(ll,jj,kk),sp )
          dpdy = real( p13(ll,jj,kk),sp )
          dpdz = real( p23(ll,jj,kk),sp )

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
          dpdx = aimag( p12(ll,jj,kk) )
          dpdy = aimag( p13(ll,jj,kk) )
          dpdz = aimag( p23(ll,jj,kk) )

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
      afntmp = sqrt(dpdx * dpdx + dpdy * dpdy + dpdz * dpdz)
      dpdx = dpdx / afntmp; dpdy = dpdy / afntmp; dpdz = dpdz / afntmp
      
      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
      alpha = abs( dpdx*evectors(1,3)+dpdy*evectors(2,3)+dpdz*evectors(3,3) )
      beta  = abs( dpdx*evectors(1,2)+dpdy*evectors(2,2)+dpdz*evectors(3,2) )
      gmma  = abs( dpdx*evectors(1,1)+dpdy*evectors(2,1)+dpdz*evectors(3,1) )

      ll=floor(alpha/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfal(ll)=pdfal(ll)+1
      ll=floor(beta/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfbe(ll)=pdfbe(ll)+1
      ll=floor(gmma/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt) pdfgm(ll)=pdfgm(ll)+1

    end do
    end do
    end do
    pdfal = pdfal * const / binw
    pdfbe = pdfbe * const / binw
    pdfgm = pdfgm * const / binw


    write(28, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(28, '(15E15.3)') (ll - 0.5_sp) * binw, pdfal(ll), pdfbe(ll), pdfgm(ll)
    end do
    write(28, *)
    write(28, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,p12,p13,p23)
  
  call destroyplan3d

  write(*,*) 'align-lagpc.x done.'

end program alignlagpc
