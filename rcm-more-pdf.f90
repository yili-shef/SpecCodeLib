program rcmmore
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: sij, aij

  integer, parameter :: npdf = 60 
  real(sp), parameter :: bnd = 20., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: poos, plnaijaij, paijaij

  real(dp) :: rmsoos, meanoos, rmslnaijaij, meanlnaijaij
  real(dp) :: meanaijaij, rmsaijaij

  integer :: ii, jj, kk, ll, nfiles, nstart, nreal, ifirst
  real(sp) :: tmp, ro, wx, wy, wz
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-more-pdf.x nfiles nstart ro nreal ifirst'
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

  str = 'more-means-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(15, file = str(1:len_trim(str)) )
  read(15, *)

  str = 'more-pdf-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(21, file = str(1:len_trim(str)) )


  do ll = nstart, nfiles + nstart -1

    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-ro='//strro(1:len_trim(strro))//'-'//str(1:len_trim(str))//'-'

    ! mean and rms parameters for normalization.
    read(15, '(15E13.4)') rmsoos, meanoos, rmsoos, meanlnaijaij, rmslnaijaij, &
                          meanaijaij, rmsaijaij

    poos = 0.d0
    plnaijaij = 0.d0
    paijaij = 0.d0
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

        sij = .5 * ( aij + transpose(aij) )
        wx = aij(3,2) - aij(2,3)
        wy = aij(1,3) - aij(3,1)
        wz = aij(2,1) - aij(1,2)

        ! vortex stretching
        tmp =  wx * wx * sij(1,1) + wy * wy * sij(2,2)  &
             + wz * wz * sij(3,3) + 2. * ( wx * wy * sij(1,2) &
             + wx * wz * sij(1,3) + wy * wz * sij(2,3) )
        tmp = (tmp - meanoos) / rmsoos

        jj = floor((tmp + bnd)/binw) + 1
        if (jj .ge. 1 .and. jj .le. npdf) poos(jj) = poos(jj) + 1

        ! pseudodissipation lnaijaij
        sij = matmul(aij, transpose(aij))
        tmp = sij(1,1) + sij(2,2) + sij(3,3) 
        tmp = tmp / meanaijaij

        jj = floor((tmp)/binw) + 1
        if (jj .ge. 1 .and. jj .le. npdf) paijaij(jj) = paijaij(jj) + 1

        tmp = log( sij(1,1) + sij(2,2) + sij(3,3) )
        tmp = (tmp - meanlnaijaij) / rmslnaijaij

        jj = floor((tmp + bnd)/binw) + 1
        if (jj .ge. 1 .and. jj .le. npdf) plnaijaij(jj) = plnaijaij(jj) + 1

      end do
  
    end do

    tmp = 1. / nreal / npnt
 
    poos = poos * tmp
    write(*,*) 'check normalization poos: ', sum(poos)
    poos = poos / binw

    plnaijaij = plnaijaij * tmp
    write(*,*) 'check normalization plnaijaij: ', sum(plnaijaij)
    plnaijaij = plnaijaij / binw

    paijaij = paijaij * tmp
    write(*,*) 'check normalization paijaij: ', sum(paijaij)
    paijaij = paijaij / binw
 
    write(21,'(''zone T = " step '',I4, ''", I='', I4, '' F = point'')') ll, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(21,'(15E15.4)') tmp, poos(ii), plnaijaij(ii), paijaij(ii)
    end do

  end do
  close(15)
  close(21)

  write(*,*) 'finished'

end program rcmmore 
