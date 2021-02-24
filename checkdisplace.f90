program defmcijcorbw
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,ll,mm,pp,nfile
  integer :: iiloc, jjloc, kkloc, iimap, jjmap, kkmap, iic, jjc, kkc

  integer, parameter :: ngroups = 512
  integer, allocatable, dimension(:, :) :: xx, yy, zz
  integer :: npntgrp

  real(sp) :: mxv, myv, mzv, xc, yc, zc, xv, yv, zv, rr, xpll, ypll, zpll, dx, dy, dz
  real(sp) :: rad1, rad2
  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./checkdisplace.x nx npntgrp inif rad1 rad2'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     npntgrp: number of points in each group'
          write(*,*) '                     inif: initial coordinats data file'
          write(*,*) '                     rad1: inner radius of the ring'
          write(*,*) '                     rad2: outer radius of the ring'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! # of points in each group
  call getarg(2,str)
  read(str, '(I20)') npntgrp

  ! initial coordinates
  call getarg(3,str)

  call getarg(4,str1)
  read(str1, '(F20.8)') rad1

  call getarg(5,str1)
  read(str1, '(F20.8)') rad2

  ny=nx; nz=nx
  dx = 2 * pi / nx
  dy = 2 * pi / ny
  dz = 2 * pi / nz

  allocate(xx(npntgrp, ngroups), yy(npntgrp, ngroups), zz(npntgrp, ngroups))

  open(30, file = './out/'//trim(str)//'.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  mxv = 0._sp; myv = 0._sp; mzv = 0._sp
  do ll = 1, ngroups

    iic = xx(1,ll); jjc = yy(1,ll); kkc = zz(1,ll)

    xc = 0._sp; yc = 0._sp; zc = 0._sp
    xv = 0._sp; yv = 0._sp; zv = 0._sp
    do mm = 2, npntgrp

      ii = xx(mm, ll); jj = yy(mm,ll); kk = zz(mm,ll)
      
      ! coordinates
      xpll = (ii-1)*dx; ypll = (jj-1)*dy; zpll = (kk-1)*dz

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
    write(*,*) 'xv', xv
  end do
  mxv = mxv / ngroups; myv = myv / ngroups; mzv = mzv / ngroups
  rr = sqrt( mxv + myv + mzv )
  mxv = sqrt(mxv); myv = sqrt(myv); mzv = sqrt(mzv)
  
  write(*,*) "mxv, myv, mzv, rr"
  write(*,*) mxv, myv, mzv, rr

  deallocate(xx,yy,zz)

  write(*,*) 'checkdisplace.x done.'

end program defmcijcorbw
