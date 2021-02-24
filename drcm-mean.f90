program rcmmean
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(dp), parameter :: dt = 1.e-8

  real(dp), dimension(three,three,npnt) :: aijall

  real(dp) :: meanaiiperp, meanaiipara, rmsaiiperp, rmsaiipara
  real(dp) :: skewaiiperp, skewaiipara, flataiiperp, flataiipara
  real(dp) :: meana12, meana13, meana21, meana23, meana31, meana32
  real(dp) :: rmsa12, rmsa13, rmsa21, rmsa23, rmsa31, rmsa32
  real(dp) :: skewa12, skewa13, skewa21, skewa23, skewa31, skewa32
  real(dp) :: flata12, flata13, flata21, flata23, flata31, flata32

  integer :: ii, jj, kk, ll, nfiles, nstart, nreal, ifirst
  real(dp) :: tmp, ro
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./drcm-mean.x nfiles nstart ro nreal ifirst'
          write(*,*) '        nfiles: number of files, each 100,000 points'
          write(*,*) '        nstart: the index of the first file'
          write(*,*) '        ro: Rossby number'
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
  call getarg(3,strro)
  read(strro, '(F10.6)') ro
  strro = adjustl(strro)
  call getarg(4,strreal)
  read(strreal, '(I20)') nreal
  call getarg(5,str)
  read(str, '(I20)') ifirst

  str = 'mean-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  write(15, *)'meanaiiperp meanaiipara meana12 meana13 meana21 meana23 meana31 meana32'

  str = 'rms-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  write(16, *)'rmsaiiperp rmsaiipara rmsa12 rmsa13 rmsa21 rmsa23 rmsa31 rmsa32'

  str = 'skew-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(17, file = str(1:len_trim(str)) )
  write(17, *)'skewaiiperp skewaiipara skewa12 skewa13 skewa21 skewa23 skewa31 skewa32'

  str = 'flat-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(18, file = str(1:len_trim(str)) )
  write(18, *)'flataiiperp flataiipara flata12 flata13 flata21 flata23 flata31 flata32'

  do ll = nstart, nfiles + nstart -1
    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-ro='//strro(1:len_trim(strro))//'-'//str(1:len_trim(str))//'-'

    meanaiipara = 0.d0; rmsaiipara = 0.d0
    meanaiiperp = 0.d0; rmsaiiperp = 0.d0
    skewaiipara = 0.d0; flataiipara = 0.d0
    skewaiiperp = 0.d0; flataiiperp = 0.d0
    
    meana12 = 0.d0; rmsa12 = 0.d0
    meana13 = 0.d0; rmsa13 = 0.d0
    meana21 = 0.d0; rmsa21 = 0.d0
    meana23 = 0.d0; rmsa23 = 0.d0
    meana31 = 0.d0; rmsa31 = 0.d0
    meana32 = 0.d0; rmsa32 = 0.d0
    skewa12 = 0.d0; flata12 = 0.d0
    skewa13 = 0.d0; flata13 = 0.d0
    skewa21 = 0.d0; flata21 = 0.d0
    skewa23 = 0.d0; flata23 = 0.d0
    skewa31 = 0.d0; flata31 = 0.d0
    skewa32 = 0.d0; flata32 = 0.d0

    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'unformatted')
        read(10) aijall
      close(10)
    
      meanaiiperp = meanaiiperp + sum(aijall(1,1,:)) + sum(aijall(2,2,:))
      rmsaiiperp = rmsaiiperp + sum(aijall(1,1,:)**2) + sum(aijall(2,2,:)**2)
      skewaiiperp = skewaiiperp + sum(aijall(1,1,:)**3) + sum(aijall(2,2,:)**3)
      flataiiperp = flataiiperp + sum(aijall(1,1,:)**4) + sum(aijall(2,2,:)**4)
  
      meanaiipara = meanaiipara + sum(aijall(3,3,:)) 
      rmsaiipara = rmsaiipara + sum(aijall(3,3,:)**2) 
      skewaiipara = skewaiipara + sum(aijall(3,3,:)**3)
      flataiipara = flataiipara + sum(aijall(3,3,:)**4)
  
      meana12 = meana12 + sum(aijall(1,2,:)) 
      rmsa12 = rmsa12 + sum(aijall(1,2,:)**2) 
      skewa12 = skewa12 + sum(aijall(1,2,:)**3) 
      flata12 = flata12 + sum(aijall(1,2,:)**4) 
  
      meana13 = meana13 + sum(aijall(1,3,:)) 
      rmsa13 = rmsa13 + sum(aijall(1,3,:)**2) 
      skewa13 = skewa13 + sum(aijall(1,3,:)**3) 
      flata13 = flata13 + sum(aijall(1,3,:)**4) 

      meana21 = meana21 + sum(aijall(2,1,:)) 
      rmsa21 = rmsa21 + sum(aijall(2,1,:)**2) 
      skewa21 = skewa21 + sum(aijall(2,1,:)**3) 
      flata21 = flata21 + sum(aijall(2,1,:)**4) 

      meana23 = meana23 + sum(aijall(2,3,:)) 
      rmsa23 = rmsa23 + sum(aijall(2,3,:)**2) 
      skewa23 = skewa23 + sum(aijall(2,3,:)**3) 
      flata23 = flata23 + sum(aijall(2,3,:)**4) 

      meana31 = meana31 + sum(aijall(3,1,:)) 
      rmsa31 = rmsa31 + sum(aijall(3,1,:)**2) 
      skewa31 = skewa31 + sum(aijall(3,1,:)**3) 
      flata31 = flata31 + sum(aijall(3,1,:)**4) 

      meana32 = meana32 + sum(aijall(3,2,:)) 
      rmsa32 = rmsa32 + sum(aijall(3,2,:)**2) 
      skewa32 = skewa32 + sum(aijall(3,2,:)**3) 
      flata32 = flata32 + sum(aijall(3,2,:)**4) 

    end do

    tmp = 1. / nreal / npnt
 
    meanaiiperp = meanaiiperp * tmp / 2
    meanaiipara = meanaiipara * tmp 
    meana12 = meana12 * tmp
    meana13 = meana13 * tmp
    meana21 = meana21 * tmp
    meana23 = meana23 * tmp
    meana31 = meana31 * tmp
    meana32 = meana32 * tmp
 
    rmsaiiperp = rmsaiiperp * tmp / 2
    flataiiperp = flataiiperp * tmp / 2
    skewaiiperp = skewaiiperp * tmp / 2
    flataiiperp = flataiiperp + 6.*meanaiiperp**2*rmsaiiperp - 4*meanaiiperp*skewaiiperp &
                   - 3.*meanaiiperp**4
    skewaiiperp = skewaiiperp - 3*rmsaiiperp*meanaiiperp + 2*meanaiiperp**3
    rmsaiiperp = sqrt( rmsaiiperp - meanaiiperp * meanaiiperp )
    skewaiiperp = skewaiiperp / rmsaiiperp**3
    flataiiperp = flataiiperp / rmsaiiperp**4
 
 
    rmsaiipara = rmsaiipara * tmp 
    flataiipara = flataiipara * tmp 
    skewaiipara = skewaiipara * tmp 
    flataiipara = flataiipara + 6.*meanaiipara**2*rmsaiipara - 4*meanaiipara*skewaiipara &
                   - 3.*meanaiipara**4
    skewaiipara = skewaiipara - 3*rmsaiipara*meanaiipara + 2*meanaiipara**3
    rmsaiipara = sqrt( rmsaiipara - meanaiipara * meanaiipara )
    skewaiipara = skewaiipara / rmsaiipara**3
    flataiipara = flataiipara / rmsaiipara**4
    
    rmsa12 = rmsa12 * tmp
    flata12 = flata12 * tmp 
    skewa12 = skewa12 * tmp 
    flata12 = flata12 + 6.*meana12**2*rmsa12 - 4*meana12*skewa12 &
                   - 3.*meana12**4
    skewa12 = skewa12 - 3*rmsa12*meana12 + 2*meana12**3
    rmsa12 = sqrt( rmsa12 - meana12 * meana12 )
    skewa12 = skewa12 / rmsa12**3
    flata12 = flata12 / rmsa12**4
 
    rmsa13 = rmsa13 * tmp
    flata13 = flata13 * tmp 
    skewa13 = skewa13 * tmp 
    flata13 = flata13 + 6.*meana13**2*rmsa13 - 4*meana13*skewa13 &
                   - 3.*meana13**4
    skewa13 = skewa13 - 3*rmsa13*meana13 + 2*meana13**3
    rmsa13 = sqrt( rmsa13 - meana13 * meana13 )
    skewa13 = skewa13 / rmsa13**3
    flata13 = flata13 / rmsa13**4
 
    rmsa21 = rmsa21 * tmp
    flata21 = flata21 * tmp 
    skewa21 = skewa21 * tmp 
    flata21 = flata21 + 6.*meana21**2*rmsa21 - 4*meana21*skewa21 &
                   - 3.*meana21**4
    skewa21 = skewa21 - 3*rmsa21*meana21 + 2*meana21**3
    rmsa21 = sqrt( rmsa21 - meana21 * meana21 )
    skewa21 = skewa21 / rmsa21**3
    flata21 = flata21 / rmsa21**4
 
    rmsa23 = rmsa23 * tmp
    flata23 = flata23 * tmp 
    skewa23 = skewa23 * tmp 
    flata23 = flata23 + 6.*meana23**2*rmsa23 - 4*meana23*skewa23 &
                   - 3.*meana23**4
    skewa23 = skewa23 - 3*rmsa23*meana23 + 2*meana23**3
    rmsa23 = sqrt( rmsa23 - meana23 * meana23 )
    skewa23 = skewa23 / rmsa23**3
    flata23 = flata23 / rmsa23**4
 
    rmsa31 = rmsa31 * tmp
    flata31 = flata31 * tmp 
    skewa31 = skewa31 * tmp 
    flata31 = flata31 + 6.*meana31**2*rmsa31 - 4*meana31*skewa31 &
                   - 3.*meana31**4
    skewa31 = skewa31 - 3*rmsa31*meana31 + 2*meana31**3
    rmsa31 = sqrt( rmsa31 - meana31 * meana31 )
    skewa31 = skewa31 / rmsa31**3
    flata31 = flata31 / rmsa31**4
 
    rmsa32 = rmsa32 * tmp
    flata32 = flata32 * tmp 
    skewa32 = skewa32 * tmp 
    flata32 = flata32 + 6.*meana32**2*rmsa32 - 4*meana32*skewa32 &
                   - 3.*meana32**4
    skewa32 = skewa32 - 3*rmsa32*meana32 + 2*meana32**3
    rmsa32 = sqrt( rmsa32 - meana32 * meana32 )
    skewa32 = skewa32 / rmsa32**3
    flata32 = flata32 / rmsa32**4
 
 
    write(15, '(15E13.4)') meanaiiperp, meanaiipara, meana12, meana13, meana21, &
                           meana23, meana31, meana32 
 
    write(16, '(15E13.4)') rmsaiiperp, rmsaiipara, rmsa12, rmsa13, rmsa21, &
                           rmsa23, rmsa31, rmsa32 
 
    write(17, '(15E13.4)') skewaiiperp, skewaiipara, skewa12, skewa13, skewa21, &
                           skewa23, skewa31, skewa32 
 
    write(18, '(15E13.4)') flataiiperp, flataiipara, flata12, flata13, flata21, &
                           flata23, flata31, flata32 

  end do
  close(15)
  close(16)
  close(17)
  close(18)

  write(*,*) 'finished'

end program rcmmean 
