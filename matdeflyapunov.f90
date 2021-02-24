program matdeflyapunov
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 80
  real(sp), parameter :: evbnd = 10._sp, evbinw = 2._sp*evbnd/npnt
  real(sp), dimension(npnt) :: pdfal, pdfbe, pdfgm

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: evectors
  integer :: ierr

  integer, parameter, dimension(3,3) :: delij = (/1,0,0,0,1,0,0,0,1/)
  real(sp), dimension(3,3) :: dij, cij

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1
  real(dp) :: al, be, ga, trc, logtrc, meanlogal, meanlogbe, meanlogga, meanlogtrc, meantrc
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./matdeflyapunov.x nx filelist prefix'
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

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'meanlogev-'//trim(flnm)//'.dat')
  open(28, file = 'pdflogev-'//trim(flnm)//'.dat')

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


    ! === PDF and mean of log eigenvalues ===
    meanlogal = 0._sp; meanlogbe = 0._sp; meanlogga = 0._sp; meanlogtrc = 0._sp
    pdfal = 0._sp; pdfbe = 0._sp; pdfgm = 0._sp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          dij(1,1) = real( b11(ll,jj,kk),sp )
          dij(1,2) = real( b12(ll,jj,kk),sp )
          dij(1,3) = real( b13(ll,jj,kk),sp )
          dij(2,1) = real( b21(ll,jj,kk),sp )
          dij(2,2) = real( b22(ll,jj,kk),sp )
          dij(2,3) = real( b23(ll,jj,kk),sp )
          dij(3,1) = real( b31(ll,jj,kk),sp )
          dij(3,2) = real( b32(ll,jj,kk),sp )
          dij(3,3) = real( b33(ll,jj,kk),sp )
      else
          ll = ii/2

          dij(1,1) = aimag( b11(ll,jj,kk) )
          dij(1,2) = aimag( b12(ll,jj,kk) )
          dij(1,3) = aimag( b13(ll,jj,kk) )
          dij(2,1) = aimag( b21(ll,jj,kk) )
          dij(2,2) = aimag( b22(ll,jj,kk) )
          dij(2,3) = aimag( b23(ll,jj,kk) )
          dij(3,1) = aimag( b31(ll,jj,kk) )
          dij(3,2) = aimag( b32(ll,jj,kk) )
          dij(3,3) = aimag( b33(ll,jj,kk) )
      endif
      cij = matmul( dij, transpose(dij) )

      trc = ( cij(1,1) + cij(2,2) + cij(3,3) )
      meantrc = meantrc + trc
      logtrc = log( cij(1,1) + cij(2,2) + cij(3,3) )
      meanlogtrc = meanlogtrc + logtrc
      
      call rs(3,3,cij,evalues,matz,evectors,fv1,fv2,ierr)
      al = log(max(evalues(3), mytiny))
      be = log(max(evalues(2), mytiny))
      ga = log(max(evalues(1), mytiny))

      meanlogal = meanlogal + al
      meanlogbe = meanlogbe + be
      meanlogga = meanlogga + ga

      al = al + evbnd
      ll = floor( al / evbinw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfal(ll) = pdfal(ll) + 1
      be = be + evbnd
      ll = floor( be / evbinw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfbe(ll) = pdfbe(ll) + 1
      ga = ga + evbnd
      ll = floor( ga / evbinw ) + 1
      if ( ll .ge. 1 .and. ll .le. npnt ) pdfgm(ll) = pdfgm(ll) + 1

    end do
    end do
    end do

    meanlogtrc = meanlogtrc * const
    meantrc = meantrc * const

    meanlogal = meanlogal * const
    meanlogbe = meanlogbe * const
    meanlogga = meanlogga * const

    pdfal = pdfal * const 
    pdfbe = pdfbe * const
    pdfgm = pdfgm * const
    write(*,*) 'normalization of pdfal', sum(pdfal)
    write(*,*) 'normalization of pdfbe', sum(pdfbe)
    write(*,*) 'normalization of pdfgm', sum(pdfgm)
    pdfal = pdfal / evbinw
    pdfbe = pdfbe / evbinw
    pdfgm = pdfgm / evbinw


    write(27,'(I6, 15E15.3)') nfile, meanlogal, meanlogbe, meanlogga, meanlogtrc, meantrc

    write(28, *) '#', nfile
    do ll = 1, npnt
      write(28, '(15E15.3)') -evbnd + (ll - 0.5_sp) * evbinw, pdfal(ll), pdfbe(ll), pdfgm(ll)
    end do
    write(28, *)
    write(28, *)

    nfile = nfile + 1

  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33)
  
  call destroyplan3d

  write(*,*) 'matdeflyapunov.x done.'

end program matdeflyapunov
