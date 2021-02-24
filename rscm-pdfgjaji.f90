program rscmpdfgjaji
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three) :: gi, gjaji
  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer, parameter :: npdf = 100 
  real(sp), parameter :: bnd = 25., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: pg1, pg2, pg3
  real(dp), dimension(npdf) :: pg3a31, pg2a21, pg1a11
  real(dp), dimension(npdf) :: pg3a32, pg2a22, pg1a12
  real(dp), dimension(npdf) :: pg3a33, pg2a23, pg1a13

  integer :: ii, jj, kk, ll, mm, nreal, ifirst
  real(sp) :: tmp
  character(80) :: str, straij, strro, strreal, str1, str2

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rscm-pdfgjaji-s(d)p.x gifilelist aijfilelist nreal ifirst'
          write(*,*) '        gifilelist: gifilelist.list is the list of gi data files'
          write(*,*) '        aijfilelist: aijfilelist.list is the list of aij data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Note: To be used with rscm-mean-s(d)p.x'
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
  pg1 = 0.d0 ; pg2 = 0.d0 ; pg3 = 0.d0
  pg3a31 = 0.d0; pg2a21 = 0.d0; pg1a11 = 0.d0
  pg3a32 = 0.d0; pg2a22 = 0.d0; pg1a12 = 0.d0;
  pg3a33 = 0.d0; pg2a23 = 0.d0; pg1a13 = 0.d0;

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

        gjaji = - matmul(transpose(aij), gi)

        do jj = 1, three
        do mm = 1, three
          aij(mm,jj) = - gi(mm)*aij(mm,jj)
        end do
        end do

        jj = floor( ( gjaji(1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg1(jj) = pg1(jj) + 1
     
        jj = floor( ( gjaji(2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg2(jj) = pg2(jj) + 1
     
        jj = floor( ( gjaji(3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg3(jj) = pg3(jj) + 1

        jj = floor( ( aij(1,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg1a11(jj) = pg1a11(jj) + 1
     
        jj = floor( ( aij(2,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg2a21(jj) = pg2a21(jj) + 1
     
        jj = floor( ( aij(3,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg3a31(jj) = pg3a31(jj) + 1

        jj = floor( ( aij(1,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg1a12(jj) = pg1a12(jj) + 1
     
        jj = floor( ( aij(2,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg2a22(jj) = pg2a22(jj) + 1
     
        jj = floor( ( aij(3,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg3a32(jj) = pg3a32(jj) + 1

        jj = floor( ( aij(1,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg1a13(jj) = pg1a13(jj) + 1
     
        jj = floor( ( aij(2,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg2a23(jj) = pg2a23(jj) + 1
     
        jj = floor( ( aij(3,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pg3a33(jj) = pg3a33(jj) + 1
      end do

    end do
    ll = ll + 1
  end do
  close(22)

  tmp = 1. / nreal / npnt / ll
 
  pg1 = pg1 * tmp
  pg2 = pg2 * tmp 
  pg3 = pg3 * tmp 

  pg1a11 = pg1a11 * tmp
  pg2a21 = pg2a21 * tmp 
  pg3a31 = pg3a31 * tmp 

  pg1a12 = pg1a12 * tmp
  pg2a22 = pg2a22 * tmp 
  pg3a32 = pg3a32 * tmp 

  pg1a13 = pg1a13 * tmp
  pg2a23 = pg2a23 * tmp 
  pg3a33 = pg3a33 * tmp 

  write(*,*) 'check normalization pg1: ', sum(pg1)
  write(*,*) 'check normalization pg2: ', sum(pg2)
  write(*,*) 'check normalization pg3: ', sum(pg3)

  pg1 = pg1 / binw
  pg2 = pg2 / binw
  pg3 = pg3 / binw

  pg1a11 = pg1a11 / binw
  pg2a21 = pg2a21 / binw 
  pg3a31 = pg3a31 / binw 

  pg1a12 = pg1a12 / binw
  pg2a22 = pg2a22 / binw 
  pg3a32 = pg3a32 / binw 

  pg1a13 = pg1a13 / binw
  pg2a23 = pg2a23 / binw 
  pg3a33 = pg3a33 / binw 

  str = 'pdfgjaji-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pdfgjaj1 pdfgjaj2 pdfgjaj3'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg1(ii), pg2(ii), pg3(ii)
  end do
  close(20)
 
  str = 'pdfgjaj1-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pdfg1a11 pdfg2a21 pdfg3a31'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg1a11(ii), pg2a21(ii), pg3a31(ii)
  end do
  close(20)
 
  str = 'pdfgjaj2-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pdfg1a12 pdfg2a22 pdfg3a32'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg1a12(ii), pg2a22(ii), pg3a32(ii)
  end do
  close(20)
 
  str = 'pdfgjaj3-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,*) '# x pdfg1a13 pdfg2a23 pdfg3a33'
  do ii = 1, npdf
    tmp = -bnd + (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pg1a13(ii), pg2a23(ii), pg3a33(ii)
  end do
  close(20)
 


  write(*,*) 'finished'

end program rscmpdfgjaji
