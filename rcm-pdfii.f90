program rcmpdf
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 15., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: paii

  real(dp) :: meanaiiperp, meanaiipara, rmsaiiperp, rmsaiipara
  real(dp) :: meanaii, rmsaii

  integer :: ii, jj, kk, ll, nfiles, nstart, nreal, ifirst
  real(sp) :: tmp, ro
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-pdfii.x nfiles nstart ro nreal ifirst'
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


  str = 'mean-ro-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  read(15,*)

  str = 'rms-ro-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  read(16,*)

  str = 'pdfcmaii-all-ro-'//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  do ll = nstart, nfiles + nstart -1

    read(15, '(15E13.4)') meanaiiperp, meanaiipara
    read(16, '(15E13.4)') rmsaiiperp, rmsaiipara
    meanaii = ( meanaiipara + meanaiiperp * 2 ) / 3
    rmsaii = sqrt( ( rmsaiipara**2 + rmsaiiperp**2 * 2 ) / 3 )
 
    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-ro='//strro(1:len_trim(strro))//'-'//str(1:len_trim(str))//'-'
 
    paii = 0.d0
    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'unformatted')
        read(10) aijall
      close(10)
    
      do ii = 1, npnt
        aij = aijall(:,:,ii)
     
        aij(1,1) = (aij(1,1) - meanaii)/rmsaii
        jj = floor( ( aij(1,1) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paii(jj) = paii(jj) + 1
     
        aij(2,2) = (aij(2,2) - meanaii)/rmsaii
        jj = floor( ( aij(2,2) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paii(jj) = paii(jj) + 1
     
     
        aij(3,3) = ( aij(3,3) - meanaii ) / rmsaii
        jj = floor( ( aij(3,3) + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) paii(jj) = paii(jj) + 1
     
      end do

    end do

    tmp = 1. / nreal / npnt / 3
 
    paii= paii * tmp 

    write(*,*) 'check normalization paii: ', sum(paii)

    paii = paii / binw

    write(20,'(''zone T = " step '',I4, ''", I='', I4, '' F = point'')') ll, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') tmp, paii(ii)
    end do
 
  end do
  close(15)
  close(16)
  close(20)

  write(*,*) 'finished'

end program rcmpdf 
