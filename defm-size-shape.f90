program defmsizeshape
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: const

  integer,  parameter :: npnt = 30
  real(sp), parameter :: bnd = 1._sp, binw = 2._sp* bnd/npnt

  real(dp), dimension(npnt) :: pdfsstar, cndsize
  real(dp), dimension(3,3)  :: cc

  integer, parameter :: matz = 5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: evectors
  integer :: ierr

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: trcc
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(dp) :: meantrcc, trcctmp, mcndsize
  real(dp) :: sstar
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-size-shape.x nx filelist prefix'
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
  allocate( trcc(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'size-cnd-shape-'//trim(flnm)//'.dat')
  write(27, '("# sstar  pdfsstar  cndsize   ")') 

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

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

    pdfsstar = 0._dp
    cndsize = 0._dp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

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

      ! size of the volume
      cc = matmul( cc, transpose(cc) )
      trcctmp = cc(1,1) + cc(2,2) + cc(3,3)

      ! shape of the volume
      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)

      evalues = log(evalues)

      sstar = sum(evalues)/3._sp ! removing the trace in log eigenvalues.
      evalues = evalues - sstar

      sstar = - 3 * sqrt(6._sp) * evalues(1) * evalues(2) * evalues(3)
      sstar = sstar / ( sum(evalues * evalues) )**1.5_sp

      ll = floor((sstar+1)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) then 
          cndsize(ll) = cndsize(ll) + trcctmp
          pdfsstar(ll) = pdfsstar(ll) + 1
      end if

    end do
    end do
    end do
    mcndsize = sum(cndsize) * const

    cndsize = cndsize / (pdfsstar + mytiny)
    pdfsstar = pdfsstar * const / binw

    cndsize = cndsize/mcndsize

    write(27, '( "# index = ", I6)') nfile
    write(27, '( "# meantrcc = ", E18.5)') meantrcc
    do ll = 1, npnt
      write(27, '(15E15.3)') -1 + (ll - 0.5_sp) * binw, pdfsstar(ll), cndsize(ll)
    end do
    write(27, *)
    write(27, *)

    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33)
  
  call destroyplan3d

  write(*,*) 'defm-size-shape.x done.'

end program defmsizeshape
