program linecijcor
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,npnt

  integer, parameter :: naverage = 200
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll
  integer  :: iiloc, jjloc, kkloc

  integer, parameter :: npdf = 50, evsize = 3
  real(sp), parameter :: binw = 1._sp/npdf

  real(sp) :: ran1
  external ran1

  real(dp), dimension(evsize,evsize)  :: cc
  real(dp), dimension(evsize) :: line
  real(dp), dimension(naverage, evsize) :: cline

  real(sp), parameter, dimension(evsize) :: iniline = (/1,0,0/)

  real(dp) :: costh
  integer :: noffset

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  integer,     allocatable, dimension(:)     :: pntcnt
  real(dp),    allocatable, dimension(:,:)   :: pcthline

  character(80) :: str, flnm, prf, str1, str4
  
  if (iargc() .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./line-cij-alignpdf.x nx filelist prefix noffset'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     noffset: offset in file number in the list'
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

  call getarg(4,str4)
  read(str4, '(I20)') noffset

  write(*,*) 'Running line-cij-alignpdf.x '//trim(str)//' '//trim(flnm)//' '//trim(prf)//' '//trim(str4)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), pcthline(npdf,npnt) )
  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )

  ! Will calculate the pdf for the "noffset" data file in each list only.
  ! For this we need to generate the file names instead of reading them from the lists. 
  ! Will use the smartrunner(.py) to run the code, so the interface is kept the same.

  str1 = flnm( len_trim(flnm) - 3 : len_trim(flnm) ) !last 4 letters are the starting data file number 
  read(str1, '(I20)') ll
  ll = ll + noffset
  write(str1, '(I20)') ll 
  str1 = adjustl(str1) ! The number of the last data file in this series
  flnm = trim(flnm)//'-'//trim(str1)//'.dat'

  open(15,file='./out/b11'//trim(flnm),form='unformatted')
    read(15) b11
  close(15)
  open(15,file='./out/b12'//trim(flnm),form='unformatted')
    read(15) b12
  close(15)
  open(15,file='./out/b13'//trim(flnm),form='unformatted')
    read(15) b13
  close(15)
  open(15,file='./out/b21'//trim(flnm),form='unformatted')
    read(15) b21
  close(15)
  open(15,file='./out/b22'//trim(flnm),form='unformatted')
    read(15) b22
  close(15)
  open(15,file='./out/b23'//trim(flnm),form='unformatted')
    read(15) b23
  close(15)
  open(15,file='./out/b31'//trim(flnm),form='unformatted')
    read(15) b31
  close(15)
  open(15,file='./out/b32'//trim(flnm),form='unformatted')
    read(15) b32
  close(15)
  open(15,file='./out/b33'//trim(flnm),form='unformatted')
    read(15) b33
  close(15)

  ll = -100
  do mm = 1, naverage

    ii = floor( ran1(ll) * nx ) + 1
    jj = floor( ran1(ll) * ny ) + 1
    kk = floor( ran1(ll) * nz ) + 1

    xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

    cc(1,1) = real( b11(ii/2+1,jj,kk),sp )
    cc(1,2) = real( b12(ii/2+1,jj,kk),sp )
    cc(1,3) = real( b13(ii/2+1,jj,kk),sp )
    cc(2,1) = real( b21(ii/2+1,jj,kk),sp )
    cc(2,2) = real( b22(ii/2+1,jj,kk),sp )
    cc(2,3) = real( b23(ii/2+1,jj,kk),sp )
    cc(3,1) = real( b31(ii/2+1,jj,kk),sp )
    cc(3,2) = real( b32(ii/2+1,jj,kk),sp )
    cc(3,3) = real( b33(ii/2+1,jj,kk),sp )

    cline(mm,:) = matmul(cc,iniline)
    cline(mm,:) = cline(mm,:)/sqrt( sum( cline(mm,:) * cline(mm,:) ) )

  end do


  pcthline = 0._dp; pntcnt = 0
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
    line = matmul(cc, iniline)
    line = line / sqrt( sum(line * line) )


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
  
          costh = abs( dot_product( line, cline(ll,:) ) )
          pp = floor( costh / binw ) + 1
          pcthline(pp,mm) = pcthline(pp,mm) + 1

          pntcnt(mm) = pntcnt(mm) + 1
 
        end if

      end do
      end do
      end do
 
    end do

  end do
  end do
  end do

  open(26, file = 'line-cij-alignpdf-'//trim(flnm))

    do mm = 1, npnt

      pcthline(:,mm) = pcthline(:,mm) / (pntcnt(mm) + mytiny) 
      write(*,*) 'Normalization:', sum(pcthline(:,mm))
      pcthline(:,mm) = pcthline(:,mm) / binw

      write(26, '( "#index = ", I6)') mm
      do pp = 1, npdf
        write(26, '(15E15.3)') (pp-0.5_sp)*binw, pcthline(pp,mm)
      end do
      write(26, *)
      write(26, *)

    end do

  close(26)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(pcthline, pntcnt)
  
  write(*,*) 'line-cij-alignpdf.x done.'

end program linecijcor
