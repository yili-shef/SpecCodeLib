program rscmmeangi
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001, bnd = 25.

  real(sp), dimension(three,npnt) :: giall

  real(dp) :: meang1, meang3, rmsg1, rmsg3
  real(dp) :: skewg1, skewg3, flatg1, flatg3
  real(dp) :: meang2
  real(dp) :: rmsg2
  real(dp) :: skewg2
  real(dp) :: flatg2

  integer :: ii, kk, ll, nreal, ifirst
  integer(8) :: pntcount
  real(dp) :: tmp
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

  str = 'mean-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  write(15, *)'meang1 meang2 meang3'

  str = 'rms-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  write(16, *)'rmsg1 rmsg2 rmsg3'

  str = 'skew-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(17, file = str(1:len_trim(str)) )
  write(17, *)'skewg1 skewg2 skewg3'

  str = 'flat-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(18, file = str(1:len_trim(str)) )
  write(18, *)'flatg1 flatg2 flatg3'

  open( 20, file = strro(1:len_trim(strro))//'.list' )
  do while ( .not. eof(20))
    read(20, *) str
    write(*,*) str
    str = adjustl(str)

    meang1 = 0.d0; rmsg1 = 0.d0
    meang2 = 0.d0; rmsg2 = 0.d0
    meang3 = 0.d0; rmsg3 = 0.d0
    skewg1 = 0.d0; flatg1 = 0.d0
    skewg2 = 0.d0; flatg2 = 0.d0
    skewg3 = 0.d0; flatg3 = 0.d0
    pntcount = 0
    
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) giall
      close(10)
    
      do ii = 1, npnt

        if (abs(giall(1,ii)) .ge. bnd) cycle
        if (abs(giall(2,ii)) .ge. bnd) cycle
        if (abs(giall(3,ii)) .ge. bnd) cycle

        pntcount = pntcount + 1

        meang1 = meang1 + (giall(1,ii))
        rmsg1 = rmsg1 + (giall(1,ii)**2) 
        skewg1 = skewg1 + (giall(1,ii)**3)
        flatg1 = flatg1 + (giall(1,ii)**4)
      
        meang2 = meang2 + (giall(2,ii)) 
        rmsg2 = rmsg2 + (giall(2,ii)**2) 
        skewg2 = skewg2 + (giall(2,ii)**3) 
        flatg2 = flatg2 + (giall(2,ii)**4) 
      
        meang3 = meang3 + (giall(3,ii)) 
        rmsg3 = rmsg3 + (giall(3,ii)**2) 
        skewg3 = skewg3 + (giall(3,ii)**3)
        flatg3 = flatg3 + (giall(3,ii)**4)
      end do
  
    end do

    tmp = 1. / pntcount
 
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
    
    write(15, '(15E13.4)') meang1, meang2, meang3
    write(16, '(15E13.4)') rmsg1,  rmsg2,  rmsg3
    write(17, '(15E13.4)') skewg1, skewg2, skewg3
    write(18, '(15E13.4)') flatg1, flatg2, flatg3

  end do
  close(20)
  close(15)
  close(16)
  close(17)
  close(18)

  write(*,*) 'finished'

end program rscmmeangi 
