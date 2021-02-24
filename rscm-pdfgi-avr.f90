program rscmpdfgi
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three) :: gi

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 15., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: pg1, pg2, pg3, pgperp

  real(dp) :: meang1, meang2, rmsg1, rmsg2
  real(dp) :: meang3, rmsg3, meangperp, rmsgperp

  integer :: ii, jj, kk, ll, nreal, ifirst
  real(sp) :: tmp, meang10, meang20, meang30, rmsg10, rmsg20, rmsg30
  real(sp) :: meangperp0, rmsgperp0
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rscm-pdfgi-avr-s(d)p.x filelist nreal ifirst'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
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

  ll = 0
  pg1 = 0.d0; meang1 = 0.d0; rmsg1 = 0.d0
  pg2 = 0.d0; meang2 = 0.d0; rmsg2 = 0.d0
  pg3 = 0.d0; meang3 = 0.d0; rmsg3 = 0.d0
  pgperp = 0.d0; meangperp = 0.d0; rmsgperp = 0.d0
  open(22, file = strro(1:len_trim(strro))//'.list')
  do while ( .not. eof(22) )

    read(22, *) str
    write(*, *) str
    str = adjustl(str)
 
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) giall
      close(10)

      if ( ll .eq. 0 .and. kk .eq. ifirst ) then
          meang10 = sum(giall(1,:)) / npnt
          meang20 = sum(giall(2,:)) / npnt
          meang30 = sum(giall(3,:)) / npnt
          meangperp0 = .5 * (meang10 + meang20)

          rmsg10 = sum( ( giall(1,:) - meang10 )**2 ) / npnt
          rmsg10 = sqrt(rmsg10)
          rmsg20 = sum( ( giall(2,:) - meang20 )**2 ) / npnt
          rmsg20 = sqrt(rmsg20)
          rmsg30 = sum( ( giall(3,:) - meang30 )**2 ) / npnt
          rmsg30 = sqrt(rmsg30)

          rmsgperp0 = sum( (giall(1,:) - meangperp0)**2 ) + sum( (giall(2,:) - meangperp0)**2 )
          rmsgperp0 = sqrt( rmsgperp0 / npnt / 2)
          write(*,*) 'estimate mean gi and gperp:', meang10, meang20, meang30, meangperp0
          write(*,*) 'estimate mean rms and rmsperp:', rmsg10, rmsg20, rmsg30, rmsgperp0
      end if
    
      do ii = 1, npnt
        gi = giall(:,ii)

        meang1 = meang1 + gi(1)
        meang2 = meang2 + gi(2)
        meang3 = meang3 + gi(3)
        meangperp = meangperp + gi(1) + gi(2)
        rmsg1 = rmsg1 + gi(1) * gi(1)
        rmsg2 = rmsg2 + gi(2) * gi(2)
        rmsg3 = rmsg3 + gi(3) * gi(3)
        rmsgperp = rmsgperp + gi(1) * gi(1) + gi(2) * gi(2)
     
        tmp = (gi(1) - meang10)/rmsg10
        jj = floor( ( tmp + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg1(jj) = pg1(jj) + 1
     
        tmp = (gi(2) - meang20)/rmsg20
        jj = floor( ( tmp + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg2(jj) = pg2(jj) + 1
     
     
        tmp = ( gi(3) - meang30 ) / rmsg30
        jj = floor( ( tmp + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg3(jj) = pg3(jj) + 1

        tmp = ( gi(1) - meangperp0 ) / rmsgperp0
        jj = floor( ( tmp + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pgperp(jj) = pgperp(jj) + 1
     
        tmp = ( gi(2) - meangperp0 ) / rmsgperp0
        jj = floor( ( tmp + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pgperp(jj) = pgperp(jj) + 1

      end do

    end do
    ll = ll + 1
 
  end do
  close(22)


  tmp = 1. / nreal / npnt / ll

  meang1 = meang1 * tmp
  meang2 = meang2 * tmp
  meang3 = meang3 * tmp
  rmsg1 = sqrt( rmsg1 * tmp - meang1 * meang1)
  rmsg2 = sqrt( rmsg2 * tmp - meang2 * meang2)
  rmsg3 = sqrt( rmsg3 * tmp - meang3 * meang3)

  meangperp = meangperp * tmp / 2
  rmsgperp = sqrt( rmsgperp * tmp / 2 - meangperp * meangperp)

 
  pg1 = pg1 * tmp
  pg2 = pg2 * tmp 
  pg3 = pg3 * tmp 

  pgperp = pgperp * tmp / 2

  write(*,*) 'check normalization pg1: ', sum(pg1)
  write(*,*) 'check normalization pg2: ', sum(pg2)
  write(*,*) 'check normalization pg3: ', sum(pg3)

  pg1 = pg1 / binw
  pg2 = pg2 / binw
  pg3 = pg3 / binw

  pgperp = pgperp / binw

  str = 'pdfavrgi-'//strro(1:len_trim(strro))//'.dat'
  open(20, file = str(1:len_trim(str)) )
  write(20, *) '#', 'rmsg1 = ', rmsg1, 'meang1 = ', meang1
  write(20, *) '#', 'rmsg2 = ', rmsg2, 'meang2 = ', meang2
  write(20, *) '#', 'rmsg3 = ', rmsg3, 'meang3 = ', meang3
  write(20, *) '#', 'rmsgperp =', rmsgperp, 'meangperp = ', meangperp
  write(20,'(''# g1 pdfg1 g2 pdfg2 g3 pdfg3 gperp pdfgperp '')')
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') (tmp*rmsg10 + meang10-meang1)/rmsg1, pg1(ii)*rmsg1/rmsg10, &
                          (tmp*rmsg20 + meang20-meang2)/rmsg2, pg2(ii)*rmsg2/rmsg20, &
                          (tmp*rmsg30 + meang30-meang3)/rmsg3, pg3(ii)*rmsg3/rmsg30, &
                          (tmp*rmsgperp0 + meangperp0 - meangperp) / rmsgperp, &
                          pgperp(ii) * rmsgperp / rmsgperp0
  end do
  close(20)

  write(*,*) 'finished'

end program rscmpdfgi
