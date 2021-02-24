program align2dspc
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer  :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nn,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt=25
  real(sp), parameter :: bwcth=1._sp/npnt, bwphi=pi/2._sp/npnt
  real(sp), dimension(npnt,npnt) :: p2dcc1insp, p2dcc2insp, p2dcc3insp

  integer,  parameter :: matz = 5 
  real(dp), dimension(3) :: evcc, evpp, fv1,fv2
  real(dp), dimension(3,3) :: evtrcc, evtrpp
  integer :: ierr

  real(sp), dimension(3,3)  :: cc, pp

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23,sp11,sp12,sp13
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: flnm, prf, str1, str
  real(dp) :: malpp, mbepp, mgapp
  real(dp) :: mcc1pp1size, mcc1pp2size, mcc1pp3size
  real(dp) :: mcc2pp1size, mcc2pp2size, mcc2pp3size
  real(dp) :: mcc3pp1size, mcc3pp2size, mcc3pp3size
  real(dp) :: mcc1pp1, mcc1pp2, mcc1pp3
  real(dp) :: mcc2pp1, mcc2pp2, mcc2pp3
  real(dp) :: mcc3pp1, mcc3pp2, mcc3pp3
  real(dp) :: cth, phi, c1, c2
  
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3dspc.x nx filelist prefix'
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
  allocate( sp11(lx1,ly,lz), sp12(lx1,ly,lz), sp13(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  open(25, file = 'mevpp-align2dspc-'//flnm(1:len_trim(flnm))//'.dat')
  open(26, file = 'mccippj-align2dspc-'//flnm(1:len_trim(flnm))//'.dat')
  open(27, file = 'p2d-align2dspc-'//flnm(1:len_trim(flnm))//'.dat')
  open(28, file = 'mccippjsize-align2dspc-'//flnm(1:len_trim(flnm))//'.dat')

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

    ! Pij
    p11 = - b11 * kx * kx
    p22 = - b11 * ky * ky
    p33 = - b11 * kz * kz
    p12 = - b11 * kx * ky
    p13 = - b11 * kx * kz
    p23 = - b11 * ky * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
 
    b11 = eye * kx * b21
    b22 = eye * ky * b31
    b33 = eye * kz * b32
    b12 = eye * ( ky * b21 + kx * b31 ) * .5_sp
    b13 = eye * ( kz * b21 + kx * b32 ) * .5_sp
    b23 = eye * ( ky * b32 + kz * b31 ) * .5_sp

    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,b23,ignore_me)
 
    sp11 = cmplx( real(p11,sp) * real(b11,sp), aimag(p11) * aimag(b11) ) + &
           cmplx( real(p12,sp) * real(b12,sp), aimag(p12) * aimag(b12) ) + &
           cmplx( real(p13,sp) * real(b13,sp), aimag(p13) * aimag(b13) ) 
    sp12 = cmplx( real(p11,sp) * real(b12,sp), aimag(p11) * aimag(b12) ) + &
           cmplx( real(p12,sp) * real(b22,sp), aimag(p12) * aimag(b22) ) + &
           cmplx( real(p13,sp) * real(b23,sp), aimag(p13) * aimag(b23) ) 
    sp13 = cmplx( real(p11,sp) * real(b13,sp), aimag(p11) * aimag(b13) ) + &
           cmplx( real(p12,sp) * real(b23,sp), aimag(p12) * aimag(b23) ) + &
           cmplx( real(p13,sp) * real(b33,sp), aimag(p13) * aimag(b33) ) 
    ! sp21
    b21  = cmplx( real(p12,sp) * real(b11,sp), aimag(p12) * aimag(b11) ) + &
           cmplx( real(p22,sp) * real(b12,sp), aimag(p22) * aimag(b12) ) + &
           cmplx( real(p23,sp) * real(b13,sp), aimag(p23) * aimag(b13) ) 
    ! sp22
    b31  = cmplx( real(p12,sp) * real(b12,sp), aimag(p12) * aimag(b12) ) + &
           cmplx( real(p22,sp) * real(b22,sp), aimag(p22) * aimag(b22) ) + &
           cmplx( real(p23,sp) * real(b23,sp), aimag(p23) * aimag(b23) ) 
    ! sp23
    b32  = cmplx( real(p12,sp) * real(b13,sp), aimag(p12) * aimag(b13) ) + &
           cmplx( real(p22,sp) * real(b23,sp), aimag(p22) * aimag(b23) ) + &
           cmplx( real(p23,sp) * real(b33,sp), aimag(p23) * aimag(b33) ) 
 
    ! sp31
    p11  = cmplx( real(p13,sp) * real(b11,sp), aimag(p13) * aimag(b11) ) + &
           cmplx( real(p23,sp) * real(b12,sp), aimag(p23) * aimag(b12) ) + &
           cmplx( real(p33,sp) * real(b13,sp), aimag(p33) * aimag(b13) ) 
    ! sp32
    p12  = cmplx( real(p13,sp) * real(b12,sp), aimag(p13) * aimag(b12) ) + &
           cmplx( real(p23,sp) * real(b22,sp), aimag(p23) * aimag(b22) ) + &
           cmplx( real(p33,sp) * real(b23,sp), aimag(p33) * aimag(b23) ) 
    ! sp33
    p22  = cmplx( real(p13,sp) * real(b13,sp), aimag(p13) * aimag(b13) ) + &
           cmplx( real(p23,sp) * real(b23,sp), aimag(p23) * aimag(b23) ) + &
           cmplx( real(p33,sp) * real(b33,sp), aimag(p33) * aimag(b33) ) 

    p13 = b21; p23 = b31; p33 = b32
 
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
 
    mcc1pp1 = 0._dp; mcc1pp2 = 0._dp; mcc1pp3 = 0._sp
    mcc2pp1 = 0._dp; mcc2pp2 = 0._dp; mcc2pp3 = 0._sp
    mcc3pp1 = 0._dp; mcc3pp2 = 0._dp; mcc3pp3 = 0._sp

    mcc1pp1size = 0._dp; mcc1pp2size = 0._dp; mcc1pp3size = 0._sp
    mcc2pp1size = 0._dp; mcc2pp2size = 0._dp; mcc2pp3size = 0._sp
    mcc3pp1size = 0._dp; mcc3pp2size = 0._dp; mcc3pp3size = 0._sp

    malpp = 0._dp
    mbepp = 0._dp
    mgapp = 0._dp

    p2dcc1insp = 0._dp
    p2dcc2insp = 0._dp
    p2dcc3insp = 0._dp
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
 
          pp(1,1) = -real( sp11(ll,jj,kk),sp )
          pp(1,2) = -real( sp12(ll,jj,kk),sp )
          pp(1,3) = -real( sp13(ll,jj,kk),sp )
          pp(2,1) = -real(  p13(ll,jj,kk),sp )
          pp(2,2) = -real(  p23(ll,jj,kk),sp )
          pp(2,3) = -real(  p33(ll,jj,kk),sp )
          pp(3,1) = -real(  p11(ll,jj,kk),sp )
          pp(3,2) = -real(  p12(ll,jj,kk),sp )
          pp(3,3) = -real(  p22(ll,jj,kk),sp )
 
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
 
          pp(1,1) = -aimag( sp11(ll,jj,kk) )
          pp(1,2) = -aimag( sp12(ll,jj,kk) )
          pp(1,3) = -aimag( sp13(ll,jj,kk) )
          pp(2,1) = -aimag(  p13(ll,jj,kk) )
          pp(2,2) = -aimag(  p23(ll,jj,kk) )
          pp(2,3) = -aimag(  p33(ll,jj,kk) )
          pp(3,1) = -aimag(  p11(ll,jj,kk) )
          pp(3,2) = -aimag(  p12(ll,jj,kk) )
          pp(3,3) = -aimag(  p22(ll,jj,kk) )
 
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
      pp = .5_sp * (pp + transpose(pp))
      cc = matmul( cc, transpose(cc) )
      
      call rs(3,3,cc,evcc,matz,evtrcc,fv1,fv2,ierr)
      ! The eigenvalues are ordered smallest first. 
      do ll=1,3
        evtrcc(:,ll)=evtrcc(:,ll)/sqrt(sum(evtrcc(:,ll)**2))
      end do

      call rs(3,3,pp,evpp,matz,evtrpp,fv1,fv2,ierr)
      ! The eigenvalues are ordered smallest first. 
      do ll=1,3
        evtrpp(:,ll)=evtrpp(:,ll)/sqrt(sum(evtrpp(:,ll)**2))
      end do

      mcc1pp1 = mcc1pp1 + ( sum( evtrcc(:,3) * evtrpp(:,3) ) )**2
      mcc1pp2 = mcc1pp2 + ( sum( evtrcc(:,3) * evtrpp(:,2) ) )**2
      mcc1pp3 = mcc1pp3 + ( sum( evtrcc(:,3) * evtrpp(:,1) ) )**2

      mcc2pp1 = mcc2pp1 + ( sum( evtrcc(:,2) * evtrpp(:,3) ) )**2
      mcc2pp2 = mcc2pp2 + ( sum( evtrcc(:,2) * evtrpp(:,2) ) )**2
      mcc2pp3 = mcc2pp3 + ( sum( evtrcc(:,2) * evtrpp(:,1) ) )**2

      mcc3pp1 = mcc3pp1 + ( sum( evtrcc(:,1) * evtrpp(:,3) ) )**2
      mcc3pp2 = mcc3pp2 + ( sum( evtrcc(:,1) * evtrpp(:,2) ) )**2
      mcc3pp3 = mcc3pp3 + ( sum( evtrcc(:,1) * evtrpp(:,1) ) )**2

      mcc1pp1size = mcc1pp1size + evpp(3) * ( sum( evtrcc(:,3) * evtrpp(:,3) ) )**2
      mcc1pp2size = mcc1pp2size + evpp(2) * ( sum( evtrcc(:,3) * evtrpp(:,2) ) )**2
      mcc1pp3size = mcc1pp3size + evpp(1) * ( sum( evtrcc(:,3) * evtrpp(:,1) ) )**2

      mcc2pp1size = mcc2pp1size + evpp(3) * ( sum( evtrcc(:,2) * evtrpp(:,3) ) )**2
      mcc2pp2size = mcc2pp2size + evpp(2) * ( sum( evtrcc(:,2) * evtrpp(:,2) ) )**2
      mcc2pp3size = mcc2pp3size + evpp(1) * ( sum( evtrcc(:,2) * evtrpp(:,1) ) )**2

      mcc3pp1size = mcc3pp1size + evpp(3) * ( sum( evtrcc(:,1) * evtrpp(:,3) ) )**2
      mcc3pp2size = mcc3pp2size + evpp(2) * ( sum( evtrcc(:,1) * evtrpp(:,2) ) )**2
      mcc3pp3size = mcc3pp3size + evpp(1) * ( sum( evtrcc(:,1) * evtrpp(:,1) ) )**2

      malpp = malpp + evpp(3)
      mbepp = mbepp + evpp(2)
      mgapp = mgapp + evpp(1)

      ! Vector alpha of cc in the frame of pp
      cth = abs( sum( evtrcc(:,3) * evtrpp(:,3) ) )
      c1 = sum( evtrcc(:,3) * evtrpp(:,2) )
      c2 = sum( evtrcc(:,3) * evtrpp(:,1) )
      c1 = c1 / sqrt(c1 * c1 + c2 * c2)
      phi = acos(abs(c1))

      ll = floor(cth/bwcth) + 1
      mm = floor(phi/bwphi) + 1
      p2dcc1insp(ll,mm) = p2dcc1insp(ll,mm) + 1
 
      ! Vector beta of cc in the frame of pp
      cth = abs( sum( evtrcc(:,2) * evtrpp(:,3) ) )
      c1 = sum( evtrcc(:,2) * evtrpp(:,2) )
      c2 = sum( evtrcc(:,2) * evtrpp(:,1) )
      c1 = c1 / sqrt(c1 * c1 + c2 * c2)
      phi = acos(abs(c1))

      ll = floor(cth/bwcth) + 1
      mm = floor(phi/bwphi) + 1
      p2dcc2insp(ll,mm) = p2dcc2insp(ll,mm) + 1
 
 
      ! Vector gamma of cc in the frame of pp
      cth = abs( sum( evtrcc(:,1) * evtrpp(:,3) ) )
      c1 = sum( evtrcc(:,1) * evtrpp(:,2) )
      c2 = sum( evtrcc(:,1) * evtrpp(:,1) )
      c1 = c1 / sqrt(c1 * c1 + c2 * c2)
      phi = acos(abs(c1))

      ll = floor(cth/bwcth) + 1
      mm = floor(phi/bwphi) + 1
      p2dcc3insp(ll,mm) = p2dcc3insp(ll,mm) + 1
 
    end do
    end do
    end do
    p2dcc1insp = p2dcc1insp * const / bwcth / bwphi
    p2dcc2insp = p2dcc2insp * const / bwcth / bwphi
    p2dcc3insp = p2dcc3insp * const / bwcth / bwphi

    mcc1pp1 = mcc1pp1 * const
    mcc1pp2 = mcc1pp2 * const
    mcc1pp3 = mcc1pp3 * const

    mcc2pp1 = mcc2pp1 * const
    mcc2pp2 = mcc2pp2 * const
    mcc2pp3 = mcc2pp3 * const
    
    mcc3pp1 = mcc3pp1 * const
    mcc3pp2 = mcc3pp2 * const
    mcc3pp3 = mcc3pp3 * const

    mcc1pp1size = mcc1pp1size * const
    mcc1pp2size = mcc1pp2size * const
    mcc1pp3size = mcc1pp3size * const

    mcc2pp1size = mcc2pp1size * const
    mcc2pp2size = mcc2pp2size * const
    mcc2pp3size = mcc2pp3size * const
    
    mcc3pp1size = mcc3pp1size * const
    mcc3pp2size = mcc3pp2size * const
    mcc3pp3size = mcc3pp3size * const

    malpp = malpp * const
    mbepp = mbepp * const
    mgapp = mgapp * const

    write(25, '(I3, 15E13.3)') nfile, malpp, mbepp, mgapp

    write(26, '(I3, 15E13.3)') nfile, mcc1pp1, mcc1pp2, mcc1pp3, mcc2pp1, mcc2pp2, mcc2pp3, &
                                      mcc3pp1, mcc3pp2, mcc3pp3
    
    write(28, '(I3, 15E13.3)') nfile, mcc1pp1size, mcc1pp2size, mcc1pp3size, &
                                      mcc2pp1size, mcc2pp2size, mcc2pp3size, &
                                      mcc3pp1size, mcc3pp2size, mcc3pp3size
    
    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
    do ii = 1, npnt
      write(27, '(15E15.3)') (ll - 0.5_sp) * bwcth, (ii-.5_sp)*bwphi, p2dcc1insp(ll,ii), &
                              p2dcc2insp(ll,ii), p2dcc3insp(ll,ii)
    end do
    write(27,*) 
    end do
    write(27, *)
    write(27, *)

    nfile = nfile + 1
  end do
  close(30)
  close(25)
  close(26)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,p11,p22,p33,p12,p13,p23,sp11,sp12,sp13)

  call destroyplan3d

  write(*,*) 'align-2dspc.x done.'

end program align2dspc
