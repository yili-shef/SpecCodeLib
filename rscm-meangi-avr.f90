program rscmmeangi
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall

  real(dp) :: meang1, meang3, rmsg1, rmsg3
  real(dp) :: skewg1, skewg3, flatg1, flatg3
  real(dp) :: meang2, meanperp
  real(dp) :: rmsg2, rmsperp
  real(dp) :: skewg2, skewperp
  real(dp) :: flatg2, flatperp

  integer :: kk, ll, nreal, ifirst
  real(sp) :: tmp
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rscm-meangi-s(d)p.x filelist nreal ifirst'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,strro)
  strro = adjustl(strro)
  call getarg(2,strreal)
  read(strreal, '(I20)') nreal
  call getarg(3,str)
  read(str, '(I20)') ifirst

  meang1 = 0.d0; rmsg1 = 0.d0
  meang2 = 0.d0; rmsg2 = 0.d0
  meang3 = 0.d0; rmsg3 = 0.d0
  skewg1 = 0.d0; flatg1 = 0.d0
  skewg2 = 0.d0; flatg2 = 0.d0
  skewg3 = 0.d0; flatg3 = 0.d0
  ll = 0
  open( 20, file = strro(1:len_trim(strro))//'.list' )
  do while ( .not. eof(20))
    read(20, *) str
    write(*,*) str
    str = adjustl(str)

    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) giall
      close(10)
    
      meang1 = meang1 + sum(giall(1,:))
      rmsg1 = rmsg1 + sum(giall(1,:)**2) 
      skewg1 = skewg1 + sum(giall(1,:)**3)
      flatg1 = flatg1 + sum(giall(1,:)**4)
  
      meang2 = meang2 + sum(giall(2,:)) 
      rmsg2 = rmsg2 + sum(giall(2,:)**2) 
      skewg2 = skewg2 + sum(giall(2,:)**3) 
      flatg2 = flatg2 + sum(giall(2,:)**4) 
  
      meang3 = meang3 + sum(giall(3,:)) 
      rmsg3 = rmsg3 + sum(giall(3,:)**2) 
      skewg3 = skewg3 + sum(giall(3,:)**3)
      flatg3 = flatg3 + sum(giall(3,:)**4)
  
    end do

    ll = ll + 1

  end do
  close(20)

  tmp = 1. / nreal / npnt / ll

  meanperp = meang1 + meang2
  rmsperp = rmsg1 + rmsg2
  skewperp = skewg1 + skewg2
  flatperp= flatg1 + flatg2

  meanperp = meanperp * tmp / 2
  rmsperp = rmsperp * tmp / 2
  flatperp = flatperp * tmp / 2
  skewperp = skewperp * tmp / 2
  flatperp = flatperp + 6.*meanperp**2*rmsperp - 4*meanperp*skewperp &
                 - 3.*meanperp**4
  skewperp = skewperp - 3*rmsperp*meanperp + 2*meanperp**3
  rmsperp = sqrt( rmsperp - meanperp * meanperp )
  skewperp = skewperp / rmsperp**3
  flatperp = flatperp / rmsperp**4

  meang1 = meang1 * tmp
  meang2 = meang2 * tmp
  meang3 = meang3 * tmp 

  rmsg1 = rmsg1 * tmp
  flatg1 = flatg1 * tmp
  skewg1 = skewg1 * tmp
  flatg1 = flatg1 + 6.*meang1**2*rmsg1 - 4*meang1*skewg1 &
                 - 3.*meang1**4
  skewg1 = skewg1 - 3*rmsg1*meang1 + 2*meang1**3
  rmsg1 = sqrt( rmsg1 - meang1 * meang1 )
  skewg1 = skewg1 / rmsg1**3
  flatg1 = flatg1 / rmsg1**4

  rmsg2 = rmsg2 * tmp
  flatg2 = flatg2 * tmp 
  skewg2 = skewg2 * tmp 
  flatg2 = flatg2 + 6.*meang2**2*rmsg2 - 4*meang2*skewg2 &
                 - 3.*meang2**4
  skewg2 = skewg2 - 3*rmsg2*meang2 + 2*meang2**3
  rmsg2 = sqrt( rmsg2 - meang2 * meang2 )
  skewg2 = skewg2 / rmsg2**3
  flatg2 = flatg2 / rmsg2**4

  rmsg3 = rmsg3 * tmp 
  flatg3 = flatg3 * tmp 
  skewg3 = skewg3 * tmp 
  flatg3 = flatg3 + 6.*meang3**2*rmsg3 - 4*meang3*skewg3 &
                 - 3.*meang3**4
  skewg3 = skewg3 - 3*rmsg3*meang3 + 2*meang3**3
  rmsg3 = sqrt( rmsg3 - meang3 * meang3 )
  skewg3 = skewg3 / rmsg3**3
  flatg3 = flatg3 / rmsg3**4
  
  open(15, file = 'meanstats-'//strro(1:len_trim(strro))//'.dat')
    write(15, *)  '# statistics  g1 g2 g3 gperp'
    write(15, '(''mean     '', 15E13.4)') meang1, meang2, meang3, meanperp
    write(15, '(''rms      '', 15E13.4)') rmsg1,  rmsg2,  rmsg3, rmsperp
    write(15, '(''skewness '', 15E13.4)') skewg1, skewg2, skewg3, skewperp
    write(15, '(''flatness '', 15E13.4)') flatg1, flatg2, flatg3, flatperp
  close(15)

  write(*,*) 'finished'

end program rscmmeangi 
