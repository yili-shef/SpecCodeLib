program alignpc
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

  real(sp), parameter :: alss = 1.0_sp

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33, ss
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(sp) :: afntmp, dpdx, dpdy, dpdz, meanss, sstmp
  real(dp) :: alpha, beta, gmma
  integer  :: cnt
  
  if (iargc() .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-pc-cndss.x nx filelist prefix meanss'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     meanss: <ss> as a unit for conditioning'
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

  ! meanss
  call getarg(4,str)
  read(str, '(E15.3)') meanss

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( ss(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open (28, file = 'align-pc-cndss-'//trim(flnm)//'.dat')
  write(28,*)'# meanss = ', meanss
  write(28,*) 

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
    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) ss
    close(15)
    write(*,*) 'finishing reading ux uy uz and p'

    ! dp/dx, dp/dy, dp/dz
    p11 = ss * eye * kx ; p22 = ss * eye * ky ; p33 = ss * eye * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    b12 = eye * (kx * b22 + ky * b11)/2
    b13 = eye * (kx * b33 + kz * b11)/2
    b23 = eye * (ky * b33 + kz * b22)/2
    b11 = eye * kx * b11
    b22 = eye * ky * b22
    b33 = eye * kz * b33
    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b23,ignore_me)

    ss = cmplx( real(b11, sp) ** 2, aimag(b11) ** 2 ) + &
         cmplx( real(b22, sp) ** 2, aimag(b22) ** 2 ) + &
         cmplx( real(b33, sp) ** 2, aimag(b33) ** 2 ) + &
         2 * cmplx( real(b12, sp) ** 2, aimag(b12) ** 2 ) + &
         2 * cmplx( real(b13, sp) ** 2, aimag(b13) ** 2 ) + &
         2 * cmplx( real(b23, sp) ** 2, aimag(b23) ** 2 )

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


    pdfal = 0._sp; pdfbe = 0._sp; pdfgm = 0._sp
    cnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          dpdx = real( p11(ll,jj,kk),sp )
          dpdy = real( p22(ll,jj,kk),sp )
          dpdz = real( p33(ll,jj,kk),sp )

          cc(1,1) = real( b11(ll,jj,kk),sp )
          cc(1,2) = real( b12(ll,jj,kk),sp )
          cc(1,3) = real( b13(ll,jj,kk),sp )
          cc(2,1) = real( b21(ll,jj,kk),sp )
          cc(2,2) = real( b22(ll,jj,kk),sp )
          cc(2,3) = real( b23(ll,jj,kk),sp )
          cc(3,1) = real( b31(ll,jj,kk),sp )
          cc(3,2) = real( b32(ll,jj,kk),sp )
          cc(3,3) = real( b33(ll,jj,kk),sp )

          sstmp = real( ss(ll,jj,kk), sp )
      else
          ll = ii/2
          dpdx = aimag( p11(ll,jj,kk) )
          dpdy = aimag( p22(ll,jj,kk) )
          dpdz = aimag( p33(ll,jj,kk) )

          cc(1,1) = aimag( b11(ll,jj,kk) )
          cc(1,2) = aimag( b12(ll,jj,kk) )
          cc(1,3) = aimag( b13(ll,jj,kk) )
          cc(2,1) = aimag( b21(ll,jj,kk) )
          cc(2,2) = aimag( b22(ll,jj,kk) )
          cc(2,3) = aimag( b23(ll,jj,kk) )
          cc(3,1) = aimag( b31(ll,jj,kk) )
          cc(3,2) = aimag( b32(ll,jj,kk) )
          cc(3,3) = aimag( b33(ll,jj,kk) )
          
          sstmp = aimag( ss(ll,jj,kk) )
      endif

      if ( sstmp .ge. alss * meanss ) then
          cnt = cnt + 1

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
      end if

    end do
    end do
    end do
    pdfal = pdfal / cnt / binw
    pdfbe = pdfbe / cnt / binw
    pdfgm = pdfgm / cnt / binw


    write(28, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(28, '(15E15.3)') (ll - 0.5_sp) * binw, pdfal(ll), pdfbe(ll), pdfgm(ll)
    end do
    write(28, *)
    write(28, *)

    nfile = nfile + 1
  end do
  close(30)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,ss)
  
  call destroyplan3d

  write(*,*) 'align-pc-cndss.x done.'

end program alignpc
