program cmvortmean
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall

  real(dp) :: meanwx, meanwy, meanwz
  real(dp) :: rmswx, rmswy, rmswz
  real(dp) :: skewwx, skewwy, skewwz
  real(dp) :: flatwx, flatwy, flatwz

  integer :: ii, jj, kk, ll, nfiles, nstart, nreal, ifirst
  real(sp) :: tmp, ro, wx, wy, wz
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cmmodel-vortmean.x nfiles nstart nreal ifirst'
          write(*,*) '        nfiles: number of files, each 100,000 points'
          write(*,*) '        nstart: the index of the first file'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,str)
  read(str, '(I20)') nfiles
  call getarg(2,str)
  read(str, '(I20)') nstart
  call getarg(3,strreal)
  read(strreal, '(I20)') nreal
  call getarg(4,str)
  read(str, '(I20)') ifirst

  str = 'vortmean-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  write(15, *)'meanwx, meanwy, meanwz'

  str = 'vortrms-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  write(16, *)'rmswx, rmswy, rmswz'

  str = 'vortskew-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(17, file = str(1:len_trim(str)) )
  write(17, *)'skewwx, skewwy, skewez'

  str = 'vortflat-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(18, file = str(1:len_trim(str)) )
  write(18, *)'flatwx, flatwy, flatwz'

  do ll = nstart, nfiles + nstart -1
    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-'//str(1:len_trim(str))//'-'

    meanwx = 0.d0; rmswx = 0.d0
    meanwy = 0.d0; rmswy = 0.d0
    meanwz = 0.d0; rmswz = 0.d0
    skewwx = 0.d0; flatwx = 0.d0
    skewwy = 0.d0; flatwy = 0.d0
    skewwz = 0.d0; flatwz = 0.d0

    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) aijall
      close(10)

      do ii = 1, npnt
        wx = aijall(3,2,ii) - aijall(2,3,ii)
        wy = aijall(1,3,ii) - aijall(3,1,ii)
        wz = aijall(2,1,ii) - aijall(1,2,ii)

        meanwx = meanwx + wx
        meanwy = meanwy + wy
        meanwz = meanwz + wz

        rmswx = rmswx + wx * wx
        rmswy = rmswy + wy * wy
        rmswz = rmswz + wz * wz

        skewwx = skewwx + wx * wx * wx
        skewwy = skewwy + wy * wy * wy
        skewwz = skewwz + wz * wz * wz

        flatwx = flatwx + wx * wx * wx * wx
        flatwy = flatwy + wy * wy * wy * wy
        flatwz = flatwz + wz * wz * wz * wz
      end do
    
    end do

    tmp = 1. / nreal / npnt
 
    meanwx = meanwx * tmp
    meanwy = meanwy * tmp
    meanwz = meanwz * tmp
 
    rmswx = rmswx * tmp 
    flatwx = flatwx * tmp 
    skewwx = skewwx * tmp 
    flatwx = flatwx + 6.*meanwx**2*rmswx - 4*meanwx*skewwx &
                   - 3.*meanwx**4
    skewwx = skewwx - 3*rmswx*meanwx + 2*meanwx**3
    rmswx = sqrt( rmswx - meanwx * meanwx )
    skewwx = skewwx / rmswx**3
    flatwx = flatwx / rmswx**4
 
 
    rmswy = rmswy * tmp 
    flatwy = flatwy * tmp 
    skewwy = skewwy * tmp 
    flatwy = flatwy + 6.*meanwy**2*rmswy - 4*meanwy*skewwy &
                   - 3.*meanwy**4
    skewwy = skewwy - 3*rmswy*meanwy + 2*meanwy**3
    rmswy = sqrt( rmswy - meanwy * meanwy )
    skewwy = skewwy / rmswy**3
    flatwy = flatwy / rmswy**4
 
 
    rmswz = rmswz * tmp 
    flatwz = flatwz * tmp 
    skewwz = skewwz * tmp 
    flatwz = flatwz + 6.*meanwz**2*rmswz - 4*meanwz*skewwz &
                   - 3.*meanwz**4
    skewwz = skewwz - 3*rmswz*meanwz + 2*meanwz**3
    rmswz = sqrt( rmswz - meanwz * meanwz )
    skewwz = skewwz / rmswz**3
    flatwz = flatwz / rmswz**4
 
 
 
    write(15, '(15E13.4)') meanwx, meanwy, meanwz
    write(16, '(15E13.4)') rmswx, rmswy, rmswz
    write(17, '(15E13.4)') skewwx, skewwy, skewwz
    write(18, '(15E13.4)') flatwx, flatwy, flatwz

  end do
  close(15)
  close(16)
  close(17)
  close(18)

  write(*,*) 'finished'

end program cmvortmean 
