program rcmpost
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 15., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: paiiperp, paiipara
  real(dp), dimension(npdf) :: pa12, pa13, pa21, pa23, pa31, pa32

  real(dp) :: meanaiiperp0, meanaiipara0, rmsaiiperp0, rmsaiipara0
  real(dp) :: meana120, meana130, meana210, meana230, meana310, meana320
  real(dp) :: rmsa120, rmsa130, rmsa210, rmsa230, rmsa310, rmsa320

  real(dp) :: meanaiiperp, meanaiipara, rmsaiiperp, rmsaiipara
  real(dp) :: skewaiiperp, skewaiipara, flataiiperp, flataiipara
  real(dp) :: meana12, meana13, meana21, meana23, meana31, meana32
  real(dp) :: rmsa12, rmsa13, rmsa21, rmsa23, rmsa31, rmsa32
  real(dp) :: skewa12, skewa13, skewa21, skewa23, skewa31, skewa32
  real(dp) :: flata12, flata13, flata21, flata23, flata31, flata32

  integer :: ii, jj, ll, nfiles, nstart
  real(sp) :: tmp, ro
  character(80) :: str, strro

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-post.x nfiles nstart ro'
          write(*,*) '        nfiles: number of files, each 100,000 points'
          write(*,*) '        nstart: the index of the first file, used for output '
          write(*,*) '        ro: Rossby number'
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

  paiipara = 0.d0
  paiiperp = 0.d0
  pa12 = 0.d0
  pa13 = 0.d0
  pa21 = 0.d0
  pa23 = 0.d0
  pa31 = 0.d0
  pa32 = 0.d0
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

  do ll = nstart, nfiles + nstart -1

    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-ro='//strro(1:len_trim(strro))//'-'//str(1:len_trim(str))//'.data'
 
    open(10, file = str(1:len_trim(str)), form = 'unformatted')
      read(10) aijall
    close(10)

    if ( ll .eq. nstart ) then
      meanaiiperp0 = sum(aijall(1,1,:)) + sum(aijall(2,2,:))
      meanaiiperp0 = meanaiiperp0 / npnt / 2
      rmsaiiperp0 = sum(aijall(1,1,:)**2) + sum(aijall(2,2,:)**2)
      rmsaiiperp0 = rmsaiiperp0 / npnt / 2
      rmsaiiperp0 = sqrt( rmsaiiperp0 - meanaiiperp0**2)

      meanaiipara0 = sum(aijall(3,3,:)) 
      meanaiipara0 = meanaiipara0 / npnt 
      rmsaiipara0 = sum(aijall(3,3,:)**2) 
      rmsaiipara0 = rmsaiipara0 / npnt 
      rmsaiipara0 = sqrt( rmsaiipara0 - meanaiipara0**2)

      meana120 = sum(aijall(1,2,:)) 
      meana120 = meana120 / npnt 
      rmsa120 = sum(aijall(1,2,:)**2) 
      rmsa120 = rmsa120 / npnt 
      rmsa120 = sqrt( rmsa120 - meana120**2)

      meana130 = sum(aijall(1,3,:)) 
      meana130 = meana130 / npnt 
      rmsa130 = sum(aijall(1,3,:)**2) 
      rmsa130 = rmsa130 / npnt 
      rmsa130 = sqrt( rmsa130 - meana130**2)

      meana210 = sum(aijall(2,1,:)) 
      meana210 = meana210 / npnt 
      rmsa210 = sum(aijall(2,1,:)**2) 
      rmsa210 = rmsa210 / npnt 
      rmsa210 = sqrt( rmsa210 - meana210**2)

      meana230 = sum(aijall(2,3,:)) 
      meana230 = meana230 / npnt 
      rmsa230 = sum(aijall(2,3,:)**2) 
      rmsa230 = rmsa230 / npnt 
      rmsa230 = sqrt( rmsa230 - meana230**2)

      meana310 = sum(aijall(3,1,:)) 
      meana310 = meana310 / npnt 
      rmsa310 = sum(aijall(3,1,:)**2) 
      rmsa310 = rmsa310 / npnt 
      rmsa310 = sqrt( rmsa310 - meana310**2)

      meana320 = sum(aijall(3,2,:)) 
      meana320 = meana320 / npnt 
      rmsa320 = sum(aijall(3,2,:)**2) 
      rmsa320 = rmsa320 / npnt 
      rmsa320 = sqrt( rmsa320 - meana320**2)

      write(*,*) 'meanaiiperp0: ', meanaiiperp0, 'rmsaiiperp0: ', rmsaiiperp0
      write(*,*) 'meanaiipara0: ', meanaiipara0, 'rmsaiipara0: ', rmsaiipara0
      write(*,*) 'meana120: ', meana120, 'rmsa120: ', rmsa120
      write(*,*) 'meana130: ', meana130, 'rmsa130: ', rmsa130
      write(*,*) 'meana210: ', meana210, 'rmsa210: ', rmsa210
      write(*,*) 'meana230: ', meana230, 'rmsa230: ', rmsa230
      write(*,*) 'meana310: ', meana310, 'rmsa310: ', rmsa310
      write(*,*) 'meana320: ', meana320, 'rmsa320: ', rmsa320

    end if

    do ii = 1, npnt
      aij = aijall(:,:,ii)

      meanaiiperp = meanaiiperp + aij(1,1) + aij(2,2)
      rmsaiiperp = rmsaiiperp + aij(1,1) * aij(1,1) + aij(2,2) * aij(2,2)
      skewaiiperp = skewaiiperp + aij(1,1)**3 + aij(2,2)**3
      flataiiperp = flataiiperp + aij(1,1)**4 + aij(2,2)**4

      meanaiipara = meanaiipara + aij(3,3) 
      rmsaiipara = rmsaiiperp + aij(3,3) * aij(3,3)
      skewaiipara = skewaiipara + aij(3,3)**3
      flataiipara = flataiipara + aij(3,3)**4

      meana12 = meana12 + aij(1,2)
      rmsa12 = rmsa12 + aij(1,2) * aij(1,2)
      skewa12 = skewa12 + aij(1,2)**3
      flata12 = flata12 + aij(1,2)**4

      meana13 = meana13 + aij(1,3)
      rmsa13 = rmsa13 + aij(1,3) * aij(1,3)
      skewa13 = skewa13 + aij(1,3)**3
      flata13 = flata13 + aij(1,3)**4

      meana21 = meana21 + aij(2,1)
      rmsa21 = rmsa21 + aij(2,1) * aij(2,1)
      skewa21 = skewa21 + aij(2,1)**3
      flata21 = flata21 + aij(2,1)**4

      meana23 = meana23 + aij(2,3)
      rmsa23 = rmsa23 + aij(2,3) * aij(2,3)
      skewa23 = skewa23 + aij(2,3)**3
      flata23 = flata23 + aij(2,3)**4

      meana31 = meana31 + aij(3,1)
      rmsa31 = rmsa31 + aij(3,1) * aij(3,1)
      skewa31 = skewa31 + aij(3,1)**3
      flata31 = flata31 + aij(3,1)**4

      meana32 = meana32 + aij(3,2)
      rmsa32 = rmsa32 + aij(3,2) * aij(3,2)
      skewa32 = skewa32 + aij(3,2)**3
      flata32 = flata32 + aij(3,2)**4

      aij(1,1) = (aij(1,1) - meanaiiperp0)/rmsaiiperp0
      jj = floor( ( aij(1,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) paiiperp(jj) = paiiperp(jj) + 1

      aij(2,2) = (aij(2,2) - meanaiiperp0)/rmsaiiperp0
      jj = floor( ( aij(2,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) paiiperp(jj) = paiiperp(jj) + 1


      aij(3,3) = ( aij(3,3) - meanaiipara0 ) / rmsaiipara0
      jj = floor( ( aij(3,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) paiipara(jj) = paiipara(jj) + 1

      aij(1,2) = ( aij(1,2) - meana120 ) / rmsa120
      jj = floor( ( aij(1,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa12(jj) = pa12(jj) + 1

      aij(1,3) = ( aij(1,3) - meana130 ) / rmsa130
      jj = floor( ( aij(1,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa13(jj) = pa13(jj) + 1

      aij(2,1) = ( aij(2,1) - meana210 ) / rmsa210
      jj = floor( ( aij(2,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa21(jj) = pa21(jj) + 1

      aij(2,3) = ( aij(2,3) - meana230 ) / rmsa230
      jj = floor( ( aij(2,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa23(jj) = pa23(jj) + 1

      aij(3,1) = ( aij(3,1) - meana310 ) / rmsa310
      jj = floor( ( aij(3,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa31(jj) = pa31(jj) + 1

      aij(3,2) = ( aij(3,2) - meana320 ) / rmsa320
      jj = floor( ( aij(3,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npdf ) pa32(jj) = pa32(jj) + 1

    end do

  end do

  tmp = 1. / nfiles / npnt

  paiiperp = paiiperp * tmp / 2 
  paiipara = paiipara * tmp 
  pa12 = pa12 * tmp
  pa13 = pa13 * tmp
  pa21 = pa21 * tmp
  pa23 = pa23 * tmp
  pa31 = pa31 * tmp 
  pa32 = pa32 * tmp

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
  flataiiperp = flataiiperp + 6.*rmsaiiperp**2 - 8*meanaiiperp*skewaiiperp + &
                meanaiiperp**4
  skewaiiperp = skewaiiperp - 3*rmsaiiperp*meanaiiperp + 2*meanaiiperp**3
  rmsaiiperp = sqrt( rmsaiiperp - meanaiiperp * meanaiiperp )
  skewaiiperp = skewaiiperp / rmsaiiperp**3
  flataiiperp = flataiiperp / rmsaiiperp**4


  rmsaiipara = rmsaiipara * tmp 
  flataiipara = flataiipara * tmp 
  skewaiipara = skewaiipara * tmp 
  flataiipara = flataiipara + 6.*rmsaiipara**2 - 8*meanaiipara*skewaiipara + &
                meanaiipara**4
  skewaiipara = skewaiipara - 3*rmsaiipara*meanaiipara + 2*meanaiipara**3
  rmsaiipara = sqrt( rmsaiipara - meanaiipara * meanaiipara )
  skewaiipara = skewaiipara / rmsaiipara**3
  flataiipara = flataiipara / rmsaiipara**4
  
  rmsa12 = rmsa12 * tmp
  flata12 = flata12 * tmp 
  skewa12 = skewa12 * tmp 
  flata12 = flata12 + 6.*rmsa12**2 - 8*meana12*skewa12 + &
                meana12**4
  skewa12 = skewa12 - 3*rmsa12*meana12 + 2*meana12**3
  rmsa12 = sqrt( rmsa12 - meana12 * meana12 )
  skewa12 = skewa12 / rmsa12**3
  flata12 = flata12 / rmsa12**4

  rmsa13 = rmsa13 * tmp
  flata13 = flata13 * tmp 
  skewa13 = skewa13 * tmp 
  flata13 = flata13 + 6.*rmsa13**2 - 8*meana13*skewa13 + &
                meana13**4
  skewa13 = skewa13 - 3*rmsa13*meana13 + 2*meana13**3
  rmsa13 = sqrt( rmsa13 - meana13 * meana13 )
  skewa13 = skewa13 / rmsa13**3
  flata13 = flata13 / rmsa13**4

  rmsa21 = rmsa21 * tmp
  flata21 = flata21 * tmp 
  skewa21 = skewa21 * tmp 
  flata21 = flata21 + 6.*rmsa21**2 - 8*meana21*skewa21 + &
                meana21**4
  skewa21 = skewa21 - 3*rmsa21*meana21 + 2*meana21**3
  rmsa21 = sqrt( rmsa21 - meana21 * meana21 )
  skewa21 = skewa21 / rmsa21**3
  flata21 = flata21 / rmsa21**4

  rmsa23 = rmsa23 * tmp
  flata23 = flata23 * tmp 
  skewa23 = skewa23 * tmp 
  flata23 = flata23 + 6.*rmsa23**2 - 8*meana23*skewa23 + &
                meana23**4
  skewa23 = skewa23 - 3*rmsa23*meana23 + 2*meana23**3
  rmsa23 = sqrt( rmsa23 - meana23 * meana23 )
  skewa23 = skewa23 / rmsa23**3
  flata23 = flata23 / rmsa23**4

  rmsa31 = rmsa31 * tmp
  flata31 = flata31 * tmp 
  skewa31 = skewa31 * tmp 
  flata31 = flata31 + 6.*rmsa31**2 - 8*meana31*skewa31 + &
                meana31**4
  skewa31 = skewa31 - 3*rmsa31*meana31 + 2*meana31**3
  rmsa31 = sqrt( rmsa31 - meana31 * meana31 )
  skewa31 = skewa31 / rmsa31**3
  flata31 = flata31 / rmsa31**4

  rmsa32 = rmsa32 * tmp
  flata32 = flata32 * tmp 
  skewa32 = skewa32 * tmp 
  flata32 = flata32 + 6.*rmsa32**2 - 8*meana32*skewa32 + &
                meana32**4
  skewa32 = skewa32 - 3*rmsa32*meana32 + 2*meana32**3
  rmsa32 = sqrt( rmsa32 - meana32 * meana32 )
  skewa32 = skewa32 / rmsa32**3
  flata32 = flata32 / rmsa32**4

  write(*,*) 'check normalization paiipara: ', sum(paiipara)
  write(*,*) 'check normalization paiiperp: ', sum(paiiperp)
  write(*,*) 'check normalization pa12: ', sum(pa12)
  write(*,*) 'check normalization pa13: ', sum(pa13)
  write(*,*) 'check normalization pa21: ', sum(pa21)
  write(*,*) 'check normalization pa23: ', sum(pa23)
  write(*,*) 'check normalization pa31: ', sum(pa31)
  write(*,*) 'check normalization pa32: ', sum(pa32)

  paiiperp = paiiperp / binw
  paiipara = paiipara / binw
  pa12 = pa12 / binw
  pa13 = pa13 / binw
  pa21 = pa21 / binw
  pa23 = pa23 / binw
  pa31 = pa31 / binw
  pa32 = pa32 / binw

  open(15, file = 'mean-ro='//strro(1:len_trim(strro))//'.dat')
  write(15, *)'meanaiiperp meanaiipara meana12 meana13 meana21 meana23 meana31 meana32'
  write(15, '(15E13.4)') meanaiiperp, meanaiipara, meana12, meana13, meana21, &
                         meana23, meana31, meana32 
  close(15) 

  open(15, file = 'rms-ro='//strro(1:len_trim(strro))//'.dat')
  write(15, *)'rmsaiiperp rmsaiipara rmsa12 rmsa13 rmsa21 rmsa23 rmsa31 rmsa32'
  write(15, '(15E13.4)') rmsaiiperp, rmsaiipara, rmsa12, rmsa13, rmsa21, &
                         rmsa23, rmsa31, rmsa32 
  close(15) 

  open(15, file = 'skew-ro='//strro(1:len_trim(strro))//'.dat')
  write(15, *)'skewaiiperp skewaiipara skewa12 skewa13 skewa21 skewa23 skewa31 skewa32'
  write(15, '(15E13.4)') skewaiiperp, skewaiipara, skewa12, skewa13, skewa21, &
                         skewa23, skewa31, skewa32 
  close(15) 

  open(15, file = 'flat-ro='//strro(1:len_trim(strro))//'.dat')
  write(15, *)'flataiiperp flataiipara flata12 flata13 flata21 flata23 flata31 flata32'
  write(15, '(15E13.4)') flataiiperp, flataiipara, flata12, flata13, flata21, &
                         flata23, flata31, flata32 
  close(15) 

  open(15, file = 'pdfcmaii-ro='//strro(1:len_trim(strro))//'.dat')
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(15,'(15E15.4)') tmp*rmsaiiperp0/rmsaiiperp+meanaiiperp0-meanaiiperp, & 
                          paiiperp(ii)*rmsaiiperp/rmsaiiperp0, &
                          tmp*rmsaiipara0/rmsaiipara+meanaiipara0-meanaiipara, &
                          paiipara(ii)*rmsaiipara/rmsaiipara0
  end do
  close(15)

  open(15, file = 'pdfcmaij-ro='//strro(1:len_trim(strro))//'.dat')
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(15,'(15E15.4)') tmp*rmsa120/rmsa12+meana120-meana12, pa12(ii)*rmsa12/rmsa120, &
                          tmp*rmsa130/rmsa13+meana130-meana13, pa13(ii)*rmsa13/rmsa130, &
                          tmp*rmsa210/rmsa21+meana210-meana21, pa21(ii)*rmsa21/rmsa210, &
                          tmp*rmsa230/rmsa23+meana230-meana23, pa23(ii)*rmsa23/rmsa230, &
                          tmp*rmsa310/rmsa31+meana310-meana31, pa31(ii)*rmsa31/rmsa310, &
                          tmp*rmsa320/rmsa32+meana320-meana32, pa32(ii)*rmsa32/rmsa320
  end do
  close(15)

  write(*,*) 'finished'

end program rcmpost 
