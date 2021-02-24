program rscmpdf
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

  integer :: ii, jj, kk, ll, nreal, ifirst
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
          write(*,*) ' Usage: ./rcm-pdf-s(d)p.x filelist nreal ifirst'
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

  str = 'mean-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  read(15,*)

  str = 'rms-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  read(16,*)

  str = 'pdfaii-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  str = 'pdfaij-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(21, file = str(1:len_trim(str)) )

  open( 22, file = strro(1:len_trim(strro))//'.list' )
  do while ( .not. eof(22))

    read(15, '(15E13.4)') meanaiiperp, meanaiipara, meana12, meana13, meana21, &
                          meana23, meana31, meana32 
 
    read(16, '(15E13.4)') rmsaiiperp, rmsaiipara, rmsa12, rmsa13, rmsa21, &
                          rmsa23, rmsa31, rmsa32 
 
    read(22,*) str
    write(*, *) str
    str = adjustl(str)
 
    paiipara = 0.d0
    paiiperp = 0.d0
    pa12 = 0.d0
    pa13 = 0.d0
    pa21 = 0.d0
    pa23 = 0.d0
    pa31 = 0.d0
    pa32 = 0.d0
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) aijall
      close(10)
    
      do ii = 1, npnt
        aij = aijall(:,:,ii)
     
        aij(1,1) = (aij(1,1) - meanaiiperp)/rmsaiiperp
        jj = floor( ( aij(1,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paiiperp(jj) = paiiperp(jj) + 1
     
        aij(2,2) = (aij(2,2) - meanaiiperp)/rmsaiiperp
        jj = floor( ( aij(2,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paiiperp(jj) = paiiperp(jj) + 1
     
     
        aij(3,3) = ( aij(3,3) - meanaiipara ) / rmsaiipara
        jj = floor( ( aij(3,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paiipara(jj) = paiipara(jj) + 1
     
        aij(1,2) = ( aij(1,2) - meana12 ) / rmsa12
        jj = floor( ( aij(1,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa12(jj) = pa12(jj) + 1
     
        aij(1,3) = ( aij(1,3) - meana13 ) / rmsa13
        jj = floor( ( aij(1,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa13(jj) = pa13(jj) + 1
     
        aij(2,1) = ( aij(2,1) - meana21 ) / rmsa21
        jj = floor( ( aij(2,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa21(jj) = pa21(jj) + 1
     
        aij(2,3) = ( aij(2,3) - meana23 ) / rmsa23
        jj = floor( ( aij(2,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa23(jj) = pa23(jj) + 1
     
        aij(3,1) = ( aij(3,1) - meana31 ) / rmsa31
        jj = floor( ( aij(3,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa31(jj) = pa31(jj) + 1
     
        aij(3,2) = ( aij(3,2) - meana32 ) / rmsa32
        jj = floor( ( aij(3,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pa32(jj) = pa32(jj) + 1
     
      end do

    end do

    tmp = 1. / nreal / npnt
 
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

    write(20,'(''zone T = " step '',I4, ''", I='', I4, '' F = point'')') ll, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') tmp, paiiperp(ii), paiipara(ii)
    end do
 
    write(21,'(''zone T = " step '',I4, ''", I='', I4, '' F = point'')') ll, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(21,'(15E15.4)') tmp, pa12(ii), pa13(ii), pa21(ii), pa23(ii), pa31(ii), pa32(ii)
    end do

  end do
  close(15)
  close(16)
  close(20)
  close(21)
  close(22)

  write(*,*) 'finished'

end program rscmpdf 
