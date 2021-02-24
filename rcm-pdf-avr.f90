program rcmpdfavr
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

  real(dp) :: meanaiiperp, meanaiipara, rmsaiiperp, rmsaiipara
  real(dp) :: meana12, meana13, meana21, meana23, meana31, meana32
  real(dp) :: rmsa12, rmsa13, rmsa21, rmsa23, rmsa31, rmsa32

  real(dp) :: meanaiiperp0, meanaiipara0, rmsaiiperp0, rmsaiipara0
  real(dp) :: meana120, meana130, meana210, meana230, meana310, meana320
  real(dp) :: rmsa120, rmsa130, rmsa210, rmsa230, rmsa310, rmsa320

  integer :: ii, jj, kk, ll, nreal, ifirst
  real(sp) :: tmp
  character(80) :: str, strro, strreal, strout

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-pdf-avr-s(d)p.x filelist nreal ifirst'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Note: To be used with rcm-mean-s(d)p.x'
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

  open( 22, file = strro(1:len_trim(strro))//'.list' )

  ll = 0
  paiipara = 0.d0
  paiiperp = 0.d0
  pa12 = 0.d0
  pa13 = 0.d0
  pa21 = 0.d0
  pa23 = 0.d0
  pa31 = 0.d0
  pa32 = 0.d0
  meanaiipara = 0.d0
  meanaiiperp = 0.d0
  meana12 = 0.d0
  meana13 = 0.d0
  meana21 = 0.d0
  meana23 = 0.d0
  meana31 = 0.d0
  meana32 = 0.d0
  rmsaiipara = 0.d0
  rmsaiiperp = 0.d0
  rmsa12 = 0.d0
  rmsa13 = 0.d0
  rmsa21 = 0.d0
  rmsa23 = 0.d0
  rmsa31 = 0.d0
  rmsa32 = 0.d0
  do while ( .not. eof(22))

    read(22,*) str
    write(*, *) str
    str = adjustl(str)
 
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strout, '(I20)') kk
      strout = adjustl(strout)
      strout = str(1:len_trim(str))//strout(1:len_trim(strout))//'.data'
 
      open(10, file = strout(1:len_trim(strout)), form = 'binary')
        read(10) aijall
      close(10)

      if (kk .eq. ifirst .and. ll .eq. 0 ) then
          meanaiiperp0 = sum( aijall(1,1,:) + aijall(2,2,:) )
          meanaiiperp0 = meanaiiperp0 / npnt / 2
          meanaiipara0 = sum( aijall(3,3,:) ) / npnt

          meana120 = sum( aijall(1,2,:) ) / npnt
          !meana130 = sum( aijall(1,3,:) ) / npnt
          !meana210 = sum( aijall(2,1,:) ) / npnt
          !meana230 = sum( aijall(2,3,:) ) / npnt
          !meana310 = sum( aijall(3,1,:) ) / npnt
          !meana320 = sum( aijall(3,2,:) ) / npnt

          rmsaiiperp0 = sum( (aijall(1,1,:) - meanaiiperp0)**2 + (aijall(2,2,:) - meanaiiperp0)**2 )
          rmsaiiperp0 = sqrt(rmsaiiperp0 / npnt / 2)
          rmsaiipara0 = sqrt( sum( (aijall(3,3,:) - meanaiipara0)**2 ) / npnt )
          rmsa120 = sqrt( sum( (aijall(1,2,:) - meana120)**2 ) / npnt )
          !rmsa130 = sqrt( sum( (aijall(1,3,:) - meana130)**2 ) / npnt )
          !rmsa210 = sqrt( sum( (aijall(2,1,:) - meana210)**2 ) / npnt )
          !rmsa230 = sqrt( sum( (aijall(2,3,:) - meana230)**2 ) / npnt )
          !rmsa310 = sqrt( sum( (aijall(3,1,:) - meana310)**2 ) / npnt )
          !rmsa320 = sqrt( sum( (aijall(3,2,:) - meana320)**2 ) / npnt )
      end if
    
      do ii = 1, npnt
        aij = aijall(:,:,ii)
     
        meanaiiperp = meanaiiperp + aij(1,1) + aij(2,2)
        meanaiipara = meanaiipara + aij(3,3)
        meana12 = meana12 + aij(1,2)
        meana13 = meana13 + aij(1,3)
        meana21 = meana21 + aij(2,1)
        meana23 = meana23 + aij(2,3)
        meana31 = meana31 + aij(3,1)
        meana32 = meana32 + aij(3,2)

        rmsaiiperp = rmsaiiperp + aij(1,1) * aij(1,1) + aij(2,2) * aij(2,2)
        rmsaiipara = rmsaiipara + aij(3,3) * aij(3,3)
        rmsa12 = rmsa12 + aij(1,2) * aij(1,2)
        rmsa13 = rmsa13 + aij(1,3) * aij(1,3)
        rmsa21 = rmsa21 + aij(2,1) * aij(2,1)
        rmsa23 = rmsa23 + aij(2,3) * aij(2,3)
        rmsa31 = rmsa31 + aij(3,1) * aij(3,1)
        rmsa32 = rmsa32 + aij(3,2) * aij(3,2)

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
     
        aij(1,3) = ( aij(1,3) - meana120 ) / rmsa120
        jj = floor( ( aij(1,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa13(jj) = pa13(jj) + 1
     
        aij(2,1) = ( aij(2,1) - meana120 ) / rmsa120
        jj = floor( ( aij(2,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa21(jj) = pa21(jj) + 1
     
        aij(2,3) = ( aij(2,3) - meana120 ) / rmsa120
        jj = floor( ( aij(2,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa23(jj) = pa23(jj) + 1
     
        aij(3,1) = ( aij(3,1) - meana120 ) / rmsa120
        jj = floor( ( aij(3,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa31(jj) = pa31(jj) + 1
     
        aij(3,2) = ( aij(3,2) - meana120 ) / rmsa120
        jj = floor( ( aij(3,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa32(jj) = pa32(jj) + 1
     
      end do

    end do

    ll = ll + 1
  end do
  close(22)

  tmp = 1. / nreal / npnt / ll
 
  meanaiiperp = meanaiiperp * tmp / 2 
  meanaiipara = meanaiipara * tmp 
  meana12 = meana12 * tmp
  meana13 = meana13 * tmp
  meana21 = meana21 * tmp
  meana23 = meana23 * tmp
  meana31 = meana31 * tmp 
  meana32 = meana32 * tmp
 
  rmsaiiperp = rmsaiiperp * tmp / 2 
  rmsaiipara = rmsaiipara * tmp 
  rmsa12 = rmsa12 * tmp
  rmsa13 = rmsa13 * tmp
  rmsa21 = rmsa21 * tmp
  rmsa23 = rmsa23 * tmp
  rmsa31 = rmsa31 * tmp 
  rmsa32 = rmsa32 * tmp

  rmsaiiperp = sqrt( rmsaiiperp - meanaiiperp*meanaiiperp)
  rmsaiipara = sqrt( rmsaiipara - meanaiipara*meanaiipara)
  rmsa12 = sqrt(rmsa12 - meana12 * meana12) 
  rmsa13 = sqrt(rmsa13 - meana13 * meana13)
  rmsa21 = sqrt(rmsa21 - meana21 * meana21)
  rmsa23 = sqrt(rmsa23 - meana23 * meana23)
  rmsa31 = sqrt(rmsa31 - meana31 * meana31)
  rmsa32 = sqrt(rmsa32 - meana32 * meana32)

  paiiperp = paiiperp * tmp / 2 
  paiipara = paiipara * tmp 
  pa12 = pa12 * tmp
  pa13 = pa13 * tmp
  pa21 = pa21 * tmp
  pa23 = pa23 * tmp
  pa31 = pa31 * tmp 
  pa32 = pa32 * tmp

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

  str = 'pdfaiiavr-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

    write(20,'(''# x pdfaiiperp  x pdfaiipara '')') 
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') (tmp * rmsaiiperp0 + meanaiiperp0 - meanaiiperp) &
                            / rmsaiiperp, paiiperp(ii), &
                            (tmp * rmsaiipara0 + meanaiipara0 - meanaiipara) &
                            / rmsaiipara, paiipara(ii)
    end do
  close(20)
 
  str = 'pdfaijavr-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(21, file = str(1:len_trim(str)) )

    write(21,'(''# x p12 p13 p21 p23 p31 p32 '')')
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(21,'(15E15.4)') (tmp * rmsa120 + meana120 - meana12) / rmsa12, pa12(ii), &
                            (tmp * rmsa120 + meana120 - meana13) / rmsa13, pa13(ii), &
                            (tmp * rmsa120 + meana120 - meana21) / rmsa21, pa21(ii), &
                            (tmp * rmsa120 + meana120 - meana23) / rmsa23, pa23(ii), &
                            (tmp * rmsa120 + meana120 - meana31) / rmsa31, pa31(ii), &
                            (tmp * rmsa120 + meana120 - meana32) / rmsa32, pa32(ii)
    end do
  close(21)


  write(*,*) 'finished'

end program rcmpdfavr
