program rscmtrccndgi
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001
  real(sp), parameter :: gmma = 0.1
  real(sp), parameter :: timetheta = 0.4

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three) :: gi, gjaji
  real(sp), dimension(three,three) :: aij, cinv

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 15., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: pg1, pg2, pg3
  real(dp), dimension(npdf) :: cndpg1, cndpg2, cndpg3
  real(dp), dimension(npdf) :: trccndg1, trccndg2, trccndg3

  integer :: ii, jj, kk, ll, mm, nreal, ifirst
  real(sp) :: tmp, trcinv
  character(80) :: str, straij, strro, strreal, str1, str2

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rscm-trcgjajicndgi-s(d)p.x gifilelist aijfilelist nreal ifirst'
          write(*,*) '        gifilelist: gifilelist.list is the list of gi data files'
          write(*,*) '        aijfilelist: aijfilelist.list is the list of aij data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,strro)
  strro = adjustl(strro)
  call getarg(2, straij)
  straij = adjustl(straij)

  call getarg(3,strreal)
  read(strreal, '(I20)') nreal
  call getarg(4,str)
  read(str, '(I20)') ifirst


  open(22, file = strro(1:len_trim(strro))//'.list')
  open(23, file = straij(1:len_trim(straij))//'.list')

  ll = 0
  pg1 = 0.d0; pg2 = 0.d0; pg3 = 0.d0
  cndpg1 = 0.d0 ; cndpg2 = 0.d0 ; cndpg3 = 0.d0
  trccndg1 = 0.d0; trccndg2 = 0.d0; trccndg3 = 0.d0

  do while ( (.not. eof(22) ) .and. (.not. eof(23)) )
 
    read(22,*) str
    write(*,*) str
    str = adjustl(str)

    read(23,*) straij
    write(*,*) straij
    straij = adjustl(straij)
 
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(str2, '(I20)') kk

      str2 = adjustl(str2)

      str1 = straij(1:len_trim(straij))//str2(1:len_trim(str2))//'.data'
      str2 = str(1:len_trim(str))//str2(1:len_trim(str2))//'.data'
 
      open(10, file = str1(1:len_trim(str1)), form = 'binary')
        read(10) aijall
      close(10)
      open(10, file = str2(1:len_trim(str2)), form = 'binary')
        read(10) giall
      close(10)
    
      do ii = 1, npnt
        gi = giall(:,ii)
        aij = aijall(:,:,ii)

        ! Matrix exponential 
        cinv = -gmma * aij
        tmp = sqrt( sum ( cinv * cinv ) )
        call taylorch(cinv, tmp, three)
   
        cinv = matmul(transpose(cinv), cinv)
        trcinv = cinv(1,1) + cinv(2,2) + cinv(3,3)

        trcinv = - trcinv/(3*timetheta)

        do jj = 1, three
        do mm = 1, three
          aij(mm,jj) = - gi(mm)*aij(mm,jj)
        end do
        gjaji(jj) = sum(aij(:,jj))
        end do

        jj = floor( ( gi(1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) then
            pg1(jj) = pg1(jj) + 1
            cndpg1(jj) = cndpg1(jj) + gjaji(1)
            trccndg1(jj) = trccndg1(jj) + trcinv
        end if
     
        jj = floor( ( gi(2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) then 
            pg2(jj) = pg2(jj) + 1
            cndpg2(jj) = cndpg2(jj) + gjaji(2)
            trccndg2(jj) = trccndg2(jj) + trcinv
        end if
     
        jj = floor( ( gi(3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) then 
            pg3(jj) = pg3(jj) + 1
            cndpg3(jj) = cndpg3(jj) + gjaji(3)
            trccndg3(jj) = trccndg3(jj) + trcinv
        end if

      end do
    end do
    ll = ll + 1
  end do
  close(22)

  tmp = 1. / nreal / npnt / ll
 
  cndpg1 = cndpg1 / (pg1 + mytiny)
  cndpg2 = cndpg2 / (pg2 + mytiny) 
  cndpg3 = cndpg3 / (pg3 + mytiny) 

  trccndg1 = trccndg1 / (pg1 + mytiny)
  trccndg2 = trccndg2 / (pg2 + mytiny) 
  trccndg3 = trccndg3 / (pg3 + mytiny) 

  pg1 = pg1 * tmp
  pg2 = pg2 * tmp 
  pg3 = pg3 * tmp 

  write(*,*) 'check normalization pg1: ', sum(pg1)
  write(*,*) 'check normalization pg2: ', sum(pg2)
  write(*,*) 'check normalization pg3: ', sum(pg3)

  pg1 = pg1 / binw
  pg2 = pg2 / binw
  pg3 = pg3 / binw
 
  str = 'cndgjaj1-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pg1 cndpg1 trccndg1'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg1(ii), cndpg1(ii), trccndg1(ii)
  end do
  close(20)
 
  str = 'cndgjaj2-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pg2 cndpg2 trccndg2'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg2(ii), cndpg2(ii), trccndg2(ii)
  end do
  close(20)
 
  str = 'cndgjaj3-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pg3 cndpg3 trccndg3'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg3(ii), cndpg3(ii), trccndg3(ii)
  end do
  close(20)

  write(*,*) 'finished'

end program rscmtrccndgi
