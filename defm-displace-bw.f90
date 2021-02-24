program defmcijcorbw
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,ll,mm,pp,nfile
  integer :: iiloc, jjloc, kkloc, iimap, jjmap, kkmap, iic, jjc, kkc

  integer, parameter :: ngroups = 512
  integer :: npntgrp
  integer, allocatable, dimension(:, :) :: xx, yy, zz
  real(sp), allocatable, dimension(:) :: xp, yp, zp

  real(sp) :: mxv, myv, mzv, xc, yc, zc, xv, yv, zv, rr, xpll, ypll, zpll
  real(sp) :: rad1, rad2
  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 7) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-displace-bw.x nx filelist prefix inif npntgrp rad1 rad2'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     inif: initial coordinats data file'
          write(*,*) '                     npntgrp: number of points in each group'
          write(*,*) '                     rad1: inner radius of the ring'
          write(*,*) '                     rad2: outer radius of the ring'
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

  ! initial coordinates
  call getarg(4,str)

  ! # of points in each group
  call getarg(5,str1)
  read(str1, '(I20)') npntgrp

  call getarg(6,str1)
  read(str1, '(F20.8)') rad1

  call getarg(7,str1)
  read(str1, '(F20.8)') rad2

  ny=nx; nz=nx

  allocate(xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz))
  allocate(xx(npntgrp, ngroups), yy(npntgrp, ngroups), zz(npntgrp, ngroups))

  open(30, file = './out/'//trim(str)//'.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  open(27, file = "defm-displace-"//trim(str)//"-"//trim(flnm)//'.dat' )

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    if ( mod(nfile, 10) .eq. 0) write(*,*) 'data file ', str1(1:len_trim(str1))

    open(25, file = './out/xyz'//prf(1:len_trim(prf))//str1(1:len_trim(str1)), form = 'unformatted')
      read(25) xp, yp, zp
    close(25)

    mxv = 0._sp; myv = 0._sp; mzv = 0._sp
    do ll = 1, ngroups

      iic = xx(1,ll); jjc = yy(1,ll); kkc = zz(1,ll)

      xc = 0._sp; yc = 0._sp; zc = 0._sp
      xv = 0._sp; yv = 0._sp; zv = 0._sp
      do mm = 2, npntgrp
 
        ii = xx(mm, ll); jj = yy(mm,ll); kk = zz(mm,ll)
        pp = (kk-1) * nx * ny + (jj-1) * nx + ii
        xpll = xp(pp); ypll = yp(pp); zpll = zp(pp)
 

        do kkloc = -1, 1
        do jjloc = -1, 1
        do iiloc = -1, 1
 
          iimap = ii + iiloc * nx
          jjmap = jj + jjloc * ny
          kkmap = kk + kkloc * nz
          pp = floor( sqrt( real(iimap-iic)**2 + (jjmap-jjc)**2 + (kkmap-kkc)**2 ) + 0.5d0 ) + 1
  
          if ( pp .gt. rad1 .and. pp .le. rad2 ) then
              xpll = xpll + iiloc * 2 * pi
              ypll = ypll + jjloc * 2 * pi
              zpll = zpll + kkloc * 2 * pi
          end if
        end do
        end do
        end do

        xc = xc + xpll; yc = yc + ypll; zc = zc + zpll
 
        xv = xv + xpll * xpll
        yv = yv + ypll * ypll
        zv = zv + zpll * zpll
      end do
      xc = xc / (npntgrp-1); yc = yc / (npntgrp-1); zc = zc / (npntgrp-1)
      xv = xv / (npntgrp-1); yv = yv / (npntgrp-1); zv = zv / (npntgrp-1)
 
      xv = xv - xc * xc
      yv = yv - yc * yc
      zv = zv - zc * zc

      mxv = mxv + xv; myv = myv + yv; mzv = mzv + zv
    end do
    mxv = mxv / ngroups; myv = myv / ngroups; mzv = mzv / ngroups
    rr =  mxv + myv + mzv 

    write(27,'(I6, 15E14.4)') nfile, mxv, myv, mzv, rr
    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(xp, yp, zp, xx, yy, zz)
  
  write(*,*) 'defm-displace-bw.x done.'

end program defmcijcorbw
