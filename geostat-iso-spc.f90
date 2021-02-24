program geostatisospc
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer  :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nn,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: matz = 5 
  real(dp), dimension(3) :: evcc, fv1,fv2
  real(dp), dimension(3,3) :: evtrcc
  integer :: ierr

  real(sp), dimension(3,3)  :: cc, pp

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: flnm, prf, str1, str
  real(dp) :: mal, mbe, mgm, malsize, mbesize, mgmsize, alsp, besp, gmsp
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./geostat-iso-spc.x nx filelist prefix'
          write(*,*) '                       nx: resolution of data'
          write(*,*) '                       filelist: list file of data '
          write(*,*) '                       prefix: prefix for datafiles for bij'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! list file 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! prefix
  call getarg(3,prf)
  prf = adjustl(prf)

  ny = nx; nz = nx
  lx = nx/2; lx1 = nx/2 + 1; ly = ny; lz = nz
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

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  open(27, file = 'meanspc-iso-'//flnm(1:len_trim(flnm))//'.dat')

  nfile = 0
  do while ( .not. eof(30) )

    read(30,*) str1
    write(*,*) str1(1:len_trim(str1))

    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
 
    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b21
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b31
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b32
    close(15)
    write(*,*) 'finishing reading velocity and pressure'

    ! Pii
    p11 = - b11 * kx * kx
    p22 = - b11 * ky * ky
    p33 = - b11 * kz * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)

    b11 = (p11 + p22 + p33)/3
 
    p11 = eye * kx * b21
    p22 = eye * ky * b31
    p33 = eye * kz * b32
    p12 = eye * ( ky * b21 + kx * b31 ) * .5_sp
    p13 = eye * ( kz * b21 + kx * b32 ) * .5_sp
    p23 = eye * ( ky * b32 + kz * b31 ) * .5_sp

    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
    
    p11 = cmplx( real(b11,sp) * real(p11,sp), aimag(b11) * aimag(p11) )
    p12 = cmplx( real(b11,sp) * real(p12,sp), aimag(b11) * aimag(p12) )
    p13 = cmplx( real(b11,sp) * real(p13,sp), aimag(b11) * aimag(p13) )
    p22 = cmplx( real(b11,sp) * real(p22,sp), aimag(b11) * aimag(p22) )
    p23 = cmplx( real(b11,sp) * real(p23,sp), aimag(b11) * aimag(p23) )
    p33 = cmplx( real(b11,sp) * real(p33,sp), aimag(b11) * aimag(p33) )
 
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
 
    mal = 0._dp
    mbe = 0._dp
    mgm = 0._dp
    malsize = 0._dp
    mbesize = 0._dp
    mgmsize = 0._dp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
 
          pp(1,1) = real( p11(ll,jj,kk),sp )
          pp(1,2) = real( p12(ll,jj,kk),sp )
          pp(1,3) = real( p13(ll,jj,kk),sp )
          pp(2,2) = real( p22(ll,jj,kk),sp )
          pp(2,3) = real( p23(ll,jj,kk),sp )
          pp(3,3) = real( p33(ll,jj,kk),sp )
 
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
 
          pp(1,1) = aimag( p11(ll,jj,kk) )
          pp(1,2) = aimag( p12(ll,jj,kk) )
          pp(1,3) = aimag( p13(ll,jj,kk) )
          pp(2,2) = aimag( p22(ll,jj,kk) )
          pp(2,3) = aimag( p23(ll,jj,kk) )
          pp(3,3) = aimag( p33(ll,jj,kk) )
 
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
      pp(2,1) = pp(1,2); pp(3,1) = pp(1,3); pp(3,2) = pp(2,3)
      cc = matmul( cc, transpose(cc) )
      
      call rs(3,3,cc,evcc,matz,evtrcc,fv1,fv2,ierr)
      ! The eigenvalues are ordered smallest first. 
      do ll=1,3
        evtrcc(:,ll)=evtrcc(:,ll)/sqrt(sum(evtrcc(:,ll)**2))
      end do
 
      alsp = sum(evtrcc(:,3) * matmul(pp,evtrcc(:,3)))
      besp = sum(evtrcc(:,2) * matmul(pp,evtrcc(:,2)))
      gmsp = sum(evtrcc(:,1) * matmul(pp,evtrcc(:,1)))
 
      mal = mal + alsp
      mbe = mbe + besp
      mgm = mgm + gmsp
      
      malsize = malsize + evcc(3) * alsp
      mbesize = mbesize + evcc(2) * besp
      mgmsize = mgmsize + evcc(1) * gmsp
 
    end do
    end do
    end do

    mal = mal * const
    mbe = mbe * const
    mgm = mgm * const
    malsize = malsize * const
    mbesize = mbesize * const
    mgmsize = mgmsize * const

    write(27,'(I4, 15E15.3)') nfile, mal, mbe, mgm, malsize, mbesize, mgmsize

    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,p12,p13,p23)

  call destroyplan3d

  write(*,*) 'geostat-iso-spc.x done.'

end program geostatisospc
