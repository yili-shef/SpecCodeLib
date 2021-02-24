program align3dsp
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer  :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nn,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt=25
  real(sp), parameter :: bwcth=1._sp/npnt, bwphi=pi/2._sp/npnt, bwzeta=pi/2._sp/npnt
  real(sp), dimension(npnt,npnt,npnt) :: p3d
  real(sp), dimension(npnt,npnt) :: p2dpp3inss, p2dpp2inss, p2dpp1inss
  real(sp), dimension(npnt) :: p1dpp3ss2

  real(sp), dimension(3,3)  :: ss, pp

  integer,  parameter :: matz = 5 
  real(dp), dimension(3) :: evss, evpp, fv1,fv2
  real(dp), dimension(3,3) :: evtrss, evtrpp
  integer :: ierr

  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p33,p12,p13,p23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, str1
  real(dp) :: ctheta, cos1, cos2, cosphi, phi, cos3, cos4, coszeta, zeta
  
  if (iargc() .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3dsp.x nx filelist'
          write(*,*) '                       nx: resolution of data'
          write(*,*) '                       filelist: list of data file'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! list file name 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ny = nx; nz = nx
  lx = nx/2; lx1 = nx/2 + 1; ly = ny; lz = nz
  const = 1._sp/(nx*ny*nz)

  allocate( S11(lx1,ly,lz), S12(lx1,ly,lz), S13(lx1,ly,lz) )
  allocate( S22(lx1,ly,lz), S23(lx1,ly,lz), S33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(20,file = flnm(1:len_trim(flnm))//'.list')

    p3d = 0._sp
    p2dpp3inss = 0._sp
    p2dpp2inss = 0._sp
    p2dpp1inss = 0._sp
    p1dpp3ss2 = 0._sp
    nfile = 0

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

    open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
      read(15) S11
    close(15)
    write(*,*) 'finishing reading p'
 
    ! dp/dx, dp/dy, dp/dz
    p11 = - s11 * kx * kx
    p22 = - s11 * ky * ky
    p33 = - s11 * kz * kz
    p12 = - s11 * kx * ky
    p13 = - s11 * kx * kz
    p23 = - s11 * ky * kz
    call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
 
    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) S11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) S22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) S33
    close(15)
    write(*,*) 'finishing reading ux uy uz'
 
    s12 = .5_sp * eye * (kx * s22 + ky * s11)
    s13 = .5_sp * eye * (kx * s33 + kz * s11)
    s23 = .5_sp * eye * (ky * s33 + kz * s22)
    s11 = eye * kx * s11
    s22 = eye * ky * s22
    s33 = eye * kz * s33
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 
 
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          pp(1,1) = real( p11(ll,jj,kk),sp )
          pp(2,2) = real( p22(ll,jj,kk),sp )
          pp(3,3) = real( p33(ll,jj,kk),sp )
          pp(1,2) = real( p12(ll,jj,kk),sp )
          pp(1,3) = real( p13(ll,jj,kk),sp )
          pp(2,3) = real( p23(ll,jj,kk),sp )
 
          ss(1,1) = real( s11(ll,jj,kk),sp )
          ss(1,2) = real( s12(ll,jj,kk),sp )
          ss(1,3) = real( s13(ll,jj,kk),sp )
          ss(2,2) = real( s22(ll,jj,kk),sp )
          ss(2,3) = real( s23(ll,jj,kk),sp )
          ss(3,3) = real( s33(ll,jj,kk),sp )
      else
          ll = ii/2
          pp(1,1) = aimag( p11(ll,jj,kk) )
          pp(2,2) = aimag( p22(ll,jj,kk) )
          pp(3,3) = aimag( p33(ll,jj,kk) )
          pp(1,2) = aimag( p12(ll,jj,kk) )
          pp(1,3) = aimag( p13(ll,jj,kk) )
          pp(2,3) = aimag( p23(ll,jj,kk) )
 
          ss(1,1) = aimag( s11(ll,jj,kk) )
          ss(1,2) = aimag( s12(ll,jj,kk) )
          ss(1,3) = aimag( s13(ll,jj,kk) )
          ss(2,2) = aimag( s22(ll,jj,kk) )
          ss(2,3) = aimag( s23(ll,jj,kk) )
          ss(3,3) = aimag( s33(ll,jj,kk) )
      endif
      pp(2,1) = pp(1,2); pp(3,1) = pp(1,3); pp(3,2) = pp(2,3)
      ss(2,1) = ss(1,2); ss(3,1) = ss(1,3); ss(3,2) = ss(2,3)
      
      call rs(3,3,ss,evss,matz,evtrss,fv1,fv2,ierr)
      ! The eigenvalues are ordered smallest first. 
      do ll=1,3
        evtrss(:,ll)=evtrss(:,ll)/sqrt(sum(evtrss(:,ll)**2))
      end do
      call rs(3,3,pp,evpp,matz,evtrpp,fv1,fv2,ierr)
      do ll=1,3
        evtrpp(:,ll)=evtrpp(:,ll)/sqrt(sum(evtrpp(:,ll)**2))
      end do
 
      ctheta = abs( sum( evtrss(:,3) * evtrpp(:,3) ) )
      cos1 = sum( evtrpp(:,3) * evtrss(:,2) )
      cos2 = sum( evtrpp(:,3) * evtrss(:,1) )
      cosphi = cos1 / sqrt( cos1 * cos1 + cos2 * cos2 )
      phi = acos( abs( cosphi ) )
      cos3 = sum( evtrss(:,1) * evtrpp(:,2) )
      cos4 = sum( evtrss(:,1) * evtrpp(:,1) )
      coszeta = cos4 / sqrt( cos3 * cos3 + cos4 * cos4 )
      zeta = acos( abs( coszeta ) )
 
      ll=floor(ctheta/bwcth)+1
      mm=floor(phi/bwphi)+1
      nn=floor(zeta/bwzeta)+1
 
      p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1
      p2dpp3inss(ll,mm) = p2dpp3inss(ll,mm) + 1

      ! p2dpp2inss
      ctheta = abs( sum( evtrss(:,3) * evtrpp(:,2) ) )
      cos1 = sum( evtrpp(:,2) * evtrss(:,2) )
      cos2 = sum( evtrpp(:,2) * evtrss(:,1) )
      cosphi = cos1 / sqrt( cos1 * cos1 + cos2 * cos2 )
      phi = acos( abs( cosphi ) )

      ll=floor(ctheta/bwcth)+1
      mm=floor(phi/bwphi)+1
      p2dpp2inss(ll,mm) = p2dpp2inss(ll,mm) + 1
 
      ! p2dpp1inss
      ctheta = abs( sum( evtrss(:,3) * evtrpp(:,1) ) )
      cos1 = sum( evtrpp(:,1) * evtrss(:,2) )
      cos2 = sum( evtrpp(:,1) * evtrss(:,1) )
      cosphi = cos1 / sqrt( cos1 * cos1 + cos2 * cos2 )
      phi = acos( abs( cosphi ) )

      ll=floor(ctheta/bwcth)+1
      mm=floor(phi/bwphi)+1
      p2dpp1inss(ll,mm) = p2dpp1inss(ll,mm) + 1
 
      cos1 = abs( sum( evtrss(:,2) * evtrpp(:,3) ) )
      ll = floor(cos1 / bwcth) + 1
      p1dpp3ss2(ll) = p1dpp3ss2(ll) + 1
    
    end do
    end do
    end do

    nfile = nfile + 1
  end do


  p3d = p3d * const / bwcth / bwphi / bwzeta / nfile
  p2dpp3inss = p2dpp3inss * const / bwcth / bwphi / nfile
  p2dpp2inss = p2dpp2inss * const / bwcth / bwphi / nfile
  p2dpp1inss = p2dpp1inss * const / bwcth / bwphi / nfile
  p1dpp3ss2 = p1dpp3ss2 * const / bwcth / nfile

  open(15, file = 'align1dpdf-pp3ss2-3dsp-'//flnm(1:len_trim(flnm))//'.dat')
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, p1dpp3ss2(ii)
    end do
  close(15)

  ! pp3 correspondings to the largest eigenvalue
  open(15, file = 'align2dpdf-pp3inss-3dsp-'//flnm(1:len_trim(flnm))//'.dat')
    do ii=1,npnt
    do jj=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, p2dpp3inss(ii,jj)
    end do
    write(15,*) 
    end do
  close(15)

  ! pp2 correspondings to the intermediate eigenvalue
  open(15, file = 'align2dpdf-pp2inss-3dsp-'//flnm(1:len_trim(flnm))//'.dat')
    do ii=1,npnt
    do jj=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, p2dpp2inss(ii,jj)
    end do
    write(15,*) 
    end do
  close(15)

  ! pp1 correspondings to the smallest eigenvalue
  open(15, file = 'align2dpdf-pp1inss-3dsp-'//flnm(1:len_trim(flnm))//'.dat')
    do ii=1,npnt
    do jj=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, p2dpp1inss(ii,jj)
    end do
    write(15,*) 
    end do
  close(15)

  open(15, file = 'align3dpdf-3dsp-'//flnm(1:len_trim(flnm))//'.dat')
    write(15,'(''#zone t="3d align pdf", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') npnt,npnt,npnt
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
      p3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,s11,s12,s13,s22,s23,s33,p11,p22,p33,p12,p13,p23)
  call destroyplan3d

  write(*,*) 'align-3dsp.x done.'

end program align3dsp
