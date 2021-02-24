program rcmmeanavr
  use mconstant
  implicit none

  integer, parameter :: three = 3

  real(sp), allocatable, dimension(:,:,:) :: aijall

  real(dp) :: meanaiiperp, meanaiipara, rmsaiiperp, rmsaiipara
  real(dp) :: skewaiiperp, skewaiipara, flataiiperp, flataiipara
  real(dp) :: meanaiiall, rmsaiiall, skewaiiall, flataiiall
  real(dp) :: meana12, meana13, meana21, meana23, meana31, meana32
  real(dp) :: rmsa12, rmsa13, rmsa21, rmsa23, rmsa31, rmsa32
  real(dp) :: skewa12, skewa13, skewa21, skewa23, skewa31, skewa32
  real(dp) :: flata12, flata13, flata21, flata23, flata31, flata32

  integer :: ll, npnt
  real(sp) :: tmp
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-mean-avr-s(d)p.x filelist npnt '
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        npnt: number of realizations'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,strro)
  strro = adjustl(strro)
  call getarg(2,strreal)
  read(strreal, '(I20)') npnt

  allocate(aijall(three,three,npnt))

  meanaiipara = 0.d0; rmsaiipara = 0.d0
  skewaiipara = 0.d0; flataiipara = 0.d0
  meanaiiperp = 0.d0; rmsaiiperp = 0.d0
  skewaiiperp = 0.d0; flataiiperp = 0.d0
  meanaiiall  = 0.d0; rmsaiiall = 0.d0
  skewaiiall  = 0.d0; flataiiall = 0.d0
  
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

  ll = 0
  open( 20, file = strro(1:len_trim(strro))//'.list' )
  do while ( .not. eof(20))
    read(20, *) str
    write(*,*) str
    str = adjustl(str)

      open(10, file = str(1:len_trim(str)), form = 'binary')
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
  
      meanaiiall = meanaiiall + sum(aijall(1,1,:)) + sum(aijall(2,2,:)) + sum(aijall(3,3,:))
      rmsaiiall = rmsaiiall + sum(aijall(1,1,:)**2) + sum(aijall(2,2,:)**2) + sum(aijall(3,3,:)**2)
      skewaiiall = skewaiiall + sum(aijall(1,1,:)**3) + sum(aijall(2,2,:)**3) + sum(aijall(3,3,:)**3)
      flataiiall = flataiiall + sum(aijall(1,1,:)**4) + sum(aijall(2,2,:)**4) + sum(aijall(3,3,:)**4)
  
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

    ll = ll + 1
  end do
  close(20)

  tmp = 1. / npnt / ll
 
  meanaiiall  = meanaiiall  * tmp / 3
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

  rmsaiiall  = rmsaiiall  * tmp / 3
  flataiiall = flataiiall * tmp / 3
  skewaiiall = skewaiiall * tmp / 3
  flataiiall = flataiiall + 6.*meanaiiall**2*rmsaiiall - 4*meanaiiall*skewaiiall &
                 - 3.*meanaiiall**4
  skewaiiall = skewaiiall - 3*rmsaiiall*meanaiiall + 2*meanaiiall**3
  rmsaiiall = sqrt( rmsaiiall - meanaiiall * meanaiiall )
  skewaiiall = skewaiiall / rmsaiiall**3
  flataiiall = flataiiall / rmsaiiall**4
  
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
 
  str = 'meanstats-'//strro(1:len_trim(strro))//'.dat'
  open(15, file = str(1:len_trim(str)) )
    write(15, '(''# mean rms skew flat '')')
    write(15, '('' aiiall     '', 10E14.4)') meanaiiall, rmsaiiall, skewaiiall, flataiiall
    write(15, '('' aiiperp    '', 10E14.4)') meanaiiperp, rmsaiiperp, skewaiiperp, flataiiperp
    write(15, '('' aiipara    '', 10E14.4)') meanaiipara, rmsaiipara, skewaiipara, flataiipara
    write(15, '('' a12        '', 10E14.4)') meana12, rmsa12, skewa12, flata12
    write(15, '('' a13        '', 10E14.4)') meana13, rmsa13, skewa13, flata13
    write(15, '('' a21        '', 10E14.4)') meana21, rmsa21, skewa21, flata21
    write(15, '('' a23        '', 10E14.4)') meana23, rmsa23, skewa23, flata23
    write(15, '('' a31        '', 10E14.4)') meana31, rmsa31, skewa31, flata31
    write(15, '('' a32        '', 10E14.4)') meana32, rmsa32, skewa32, flata32
  close(15) 

  deallocate(aijall)
  write(*,*) 'finished'

end program rcmmeanavr
