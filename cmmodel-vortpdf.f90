program cmvortpdf
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 15., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: pwx, pwy
  real(dp), dimension(npdf) :: pwz

  real(dp) :: meanwx, meanwy, rmswx, rmswy
  real(dp) :: meanwz
  real(dp) :: rmswz

  integer :: ii, jj, kk, ll, nfiles, nstart, nreal, ifirst
  real(sp) :: tmp, wx, wy, wz
  character(80) :: str, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cmmodel-vortpdf.x nfiles nstart nreal ifirst'
          write(*,*) '        nfiles: number of files, each 100,000 points'
          write(*,*) '        nstart: the index of the first file'
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
  call getarg(3,strreal)
  read(strreal, '(I20)') nreal
  call getarg(4,str)
  read(str, '(I20)') ifirst


  str = 'vortmean-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  read(15,*)

  str = 'vortrms-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(16, file = str(1:len_trim(str)) )
  read(16,*)

  str = 'pdfcmvort-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(20, file = str(1:len_trim(str)) )

  do ll = nstart, nfiles + nstart -1

    read(15, '(15E13.4)') meanwx, meanwy, meanwz
    read(16, '(15E13.4)') rmswx, rmswy, rmswz
 
    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-'//str(1:len_trim(str))//'-'
 
    pwx = 0.d0
    pwy = 0.d0
    pwz = 0.d0
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

        wx = aij(3,2) - aij(2,3)
        wy = aij(1,3) - aij(3,1)
        wz = aij(2,1) - aij(1,2)
     
        wx = (wx - meanwx)/rmswx
        jj = floor( ( wx + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pwx(jj) = pwx(jj) + 1
     
        wy = ( wy - meanwy ) / rmswy
        jj = floor( ( wy + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pwy(jj) = pwy(jj) + 1
     
        wz = ( wz - meanwz ) / rmswz
        jj = floor( ( wz + bnd ) / binw ) + 1
        if ( jj .ge. 1 .and. jj .le. npdf ) pwz(jj) = pwz(jj) + 1
     
      end do

    end do

    tmp = 1. / nreal / npnt
 
    pwx = pwx * tmp 
    pwy = pwy * tmp 
    pwz = pwz * tmp

    write(*,*) 'check normalization pwy: ', sum(pwy)
    write(*,*) 'check normalization pwx: ', sum(pwx)
    write(*,*) 'check normalization pwz: ', sum(pwz)

    pwx = pwx / binw
    pwy = pwy / binw
    pwz = pwz / binw

    write(20,'(''zone T = " step '',I4, ''", I='', I4, '' F = point'')') ll, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') tmp, pwx(ii), pwy(ii), pwz(ii)
    end do
 
  end do
  close(15)
  close(16)
  close(20)

  write(*,*) 'finished'

end program cmvortpdf 
