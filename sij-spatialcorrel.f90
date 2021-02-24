program sijspatialcorrel
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,nfile,npnt
  real(sp) :: const, ignore_me

  integer, parameter :: naverage = 300
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll
  integer  :: iiloc, jjloc, kkloc

  integer, parameter :: evsize = 3

  real(sp) :: ran1
  external ran1

  real(dp), dimension(evsize,evsize)  :: cc
  real(dp), dimension(evsize, evsize, naverage) :: evtrcc0
  real(dp) :: costh

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  real(dp), allocatable, dimension(:) :: r_s11, r_s12, r_s13, r_s22, r_s23, r_s33
  integer,  allocatable, dimension(:) :: pntcnt

  character(80) :: str, flnm, fpath
  
  if (iargc() .ne. 2) then
    write(*,*)
    write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
    write(*,*) 
    write(*,*) ' Usage: ./sij-spatialcorrel.x nx filelist'
    write(*,*) '                     nx: resolution of data'
    write(*,*) '                     filelist: data file list: *.list'
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

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt) )
  allocate( r_s11(npnt), r_s12(npnt), r_s13(npnt) )
  allocate( r_s22(npnt), r_s23(npnt), r_s33(npnt) )
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( s33(lx1,ly,lz) )
  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,lx1,ly,lz)

  open(20, file = flnm(1 : len_trim(flnm))//'.list')

  r_s11 = 0._dp; r_s12 = 0._dp; r_s13 = 0._dp
  r_s22 = 0._dp; r_s23 = 0._dp; r_s33 = 0._dp
  pntcnt = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) "sij-spatialcorel.x: ", trim(fpath)

    open(15,file='./out/ux'//trim(fpath),form='unformatted')
      read(15) s11
    close(15)
    open(15,file='./out/uy'//trim(fpath),form='unformatted')
      read(15) s22
    close(15)
    open(15,file='./out/uz'//trim(fpath),form='unformatted')
      read(15) s33
    close(15)
  
    s12 = .5_sp * eye * (kx * s22 + ky * s11)
    s13 = .5_sp * eye * (kx * s33 + kz * s11)
    s23 = .5_sp * eye * (ky * s33 + kz * s22)
    s11 = eye * kx * s11
    s22 = eye * ky * s22
    s33 = eye * kz * s33
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
    ll = -100
    do mm = 1, naverage
 
      ii = floor( ran1(ll) * nx ) + 1
      jj = floor( ran1(ll) * ny ) + 1
      kk = floor( ran1(ll) * nz ) + 1
 
      xx(mm) = ii; yy(mm) = jj; zz(mm) = kk
 
      if ( mod(ii, 2) .eq. 1 ) then
        cc(1,1) = real( s11(ii/2+1,jj,kk),sp )
        cc(1,2) = real( s12(ii/2+1,jj,kk),sp )
        cc(1,3) = real( s13(ii/2+1,jj,kk),sp )
        cc(2,2) = real( s22(ii/2+1,jj,kk),sp )
        cc(2,3) = real( s23(ii/2+1,jj,kk),sp )
        cc(3,3) = real( s33(ii/2+1,jj,kk),sp )
      else
        cc(1,1) = aimag( s11(ii/2,jj,kk) )
        cc(1,2) = aimag( s12(ii/2,jj,kk) )
        cc(1,3) = aimag( s13(ii/2,jj,kk) )
        cc(2,2) = aimag( s22(ii/2,jj,kk) )
        cc(2,3) = aimag( s23(ii/2,jj,kk) )
        cc(3,3) = aimag( s33(ii/2,jj,kk) )
      end if
      cc(2,1) = cc(1,2); cc(3,1) = cc(1,3); cc(3,2) = cc(2,3)
      evtrcc0(:,:,mm) = cc
 
    end do
 
 
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          cc(1,1) = real( s11(ll,jj,kk),sp )
          cc(1,2) = real( s12(ll,jj,kk),sp )
          cc(1,3) = real( s13(ll,jj,kk),sp )
          cc(2,1) = real( s12(ll,jj,kk),sp )
          cc(2,2) = real( s22(ll,jj,kk),sp )
          cc(2,3) = real( s23(ll,jj,kk),sp )
          cc(3,1) = real( s13(ll,jj,kk),sp )
          cc(3,2) = real( s23(ll,jj,kk),sp )
          cc(3,3) = real( s33(ll,jj,kk),sp )
      else
          ll = ii/2
 
          cc(1,1) = aimag( s11(ll,jj,kk) )
          cc(1,2) = aimag( s12(ll,jj,kk) )
          cc(1,3) = aimag( s13(ll,jj,kk) )
          cc(2,1) = aimag( s12(ll,jj,kk) )
          cc(2,2) = aimag( s22(ll,jj,kk) )
          cc(2,3) = aimag( s23(ll,jj,kk) )
          cc(3,1) = aimag( s13(ll,jj,kk) )
          cc(3,2) = aimag( s23(ll,jj,kk) )
          cc(3,3) = aimag( s33(ll,jj,kk) )
      endif
 
 
      ! loop over different groups
      do ll = 1, naverage
 
        do kkloc = -1, 1 ! checking the images of the flow field
        do jjloc = -1, 1
        do iiloc = -1, 1
 
          xxll = xx(ll) + iiloc * nx
          yyll = yy(ll) + jjloc * ny
          zzll = zz(ll) + kkloc * nz
 
          mm = floor( sqrt( (ii-xxll)**2 + (jj-yyll)**2 + (kk-zzll)**2 ) + 0.5_sp )  + 1
  
          if ( mm .ge. 1 .and. mm .le. npnt) then

            costh = cc(1,1) * evtrcc0(1,1,ll)
            r_s11(mm) = r_s11(mm) + costh

            costh = cc(1,2) * evtrcc0(1,2,ll)
            r_s12(mm) = r_s12(mm) + costh
  
            costh = cc(1,3) * evtrcc0(1,3,ll)
            r_s13(mm) = r_s13(mm) + costh
  
            costh = cc(2,2) * evtrcc0(2,2,ll)
            r_s22(mm) = r_s22(mm) + costh
  
            costh = cc(2,3) * evtrcc0(2,3,ll)
            r_s23(mm) = r_s23(mm) + costh
  
            costh = cc(3,3) * evtrcc0(3,3,ll)
            r_s33(mm) = r_s33(mm) + costh

            pntcnt(mm) = pntcnt(mm) + 1
          end if
 
        end do
        end do
        end do
  
      end do
 
    end do
    end do
    end do

  end do
  close(20)

  r_s11 = r_s11 / (pntcnt + mytiny)
  r_s12 = r_s12 / (pntcnt + mytiny)
  r_s13 = r_s13 / (pntcnt + mytiny)
  r_s22 = r_s22 / (pntcnt + mytiny)
  r_s23 = r_s23 / (pntcnt + mytiny)
  r_s33 = r_s33 / (pntcnt + mytiny)

  open(26, file = 'sij-spatialcorrel-'//trim(flnm)//'.dat')

    do mm = 1, npnt
      write(26, '(15E15.3)') real(mm), r_s11(mm), r_s12(mm), r_s13(mm), &
                                       r_s22(mm), r_s23(mm), r_s33(mm)
    end do

  close(26)

  call destroyplan3d

  deallocate(s11,s12,s13,s22,s23,s33,kx,ky,kz)
  deallocate(r_s11,r_s12,r_s13,r_s22,r_s23,r_s33,pntcnt)
  
  write(*,*) 'sij-spatialcorrel.x'

end program sijspatialcorrel
