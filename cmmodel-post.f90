program cmmodelpost
  use mconstant
  implicit none

  integer, parameter :: nstep = 1000000, three = 3
  real(sp), parameter :: dt = 0.001

  real(sp), dimension(three,three,nstep) :: aijall
  real(sp), dimension(three,three) :: aij

  integer, parameter :: npnt = 50 
  real(sp), parameter :: bnd = 10, binw = 2.*bnd/npnt
  real(dp), dimension(npnt) :: paii, paij

  real(dp) :: meanaii, meanaij, rmsaii, rmsaij

  integer :: ii, jj, kk, ll, nini, nruns, nstart
  real(sp) :: tmp
  character(80) :: str

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cmmodel-post.x nini nruns nstart'
          write(*,*) '        nini: 0 - initialized with Gaussian, 1 - read in '
          write(*,*) '        nruns: number of runs, each 100,000 steps'
          write(*,*) '        nstart: the index of the first run, used for output '
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,str)
  read(str, '(I20)') nini
  call getarg(2,str)
  read(str, '(I20)') nruns
  call getarg(3,str)
  read(str, '(I20)') nstart

  paii = 0.d0
  paij = 0.d0
  meanaii = 0.d0
  rmsaii = 0.d0
  meanaij = 0.d0
  rmsaij = 0.d0
  do ll = nstart, nruns + nstart -1

    write(str, '(I20)') ll
    str = adjustl(str)
 
    open(10, file = 'cm-aij-time-'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(10) aijall
    close(10)

    do ii = 1, nstep
      aij = aijall(:,:,ii)

      jj =  floor( ( aij(1,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paii(jj) = paii(jj) + 1

      jj =  floor( ( aij(2,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paii(jj) = paii(jj) + 1

      jj =  floor( ( aij(3,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paii(jj) = paii(jj) + 1

      jj = floor( ( aij(1,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1

      jj = floor( ( aij(1,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1

      jj = floor( ( aij(2,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1
      
      jj = floor( ( aij(2,3) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1
      
      jj = floor( ( aij(3,1) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1

      jj = floor( ( aij(3,2) + bnd ) / binw ) + 1
      if ( jj .ge. 1 .and. jj .le. npnt ) paij(jj) = paij(jj) + 1

      meanaii = meanaii + aij(1,1) + aij(2,2) + aij(3,3)
      rmsaii = rmsaii + aij(1,1) * aij(1,1) + aij(2,2) * aij(2,2) + aij(3,3) * aij(3,3)

      meanaij = meanaij + aij(1,2) + aij(1,3) + aij(2,1) + aij(2,3) + &
                          aij(3,1) + aij(3,2)

      rmsaij = rmsaij + aij(1,2)**2 + aij(1,3)**2 + aij(2,1)**2 + aij(2,3)**2 &
                      + aij(3,1)**2 + aij(3,2)**2
    end do

  end do

  paii = paii / nruns / nstep / 3 
  paij = paij / nruns / nstep / 6 

  meanaii = meanaii / nruns / nstep / 3
  meanaij = meanaij / nruns / nstep / 6

  rmsaii = rmsaii / nruns / nstep / 3
  rmsaii = sqrt( rmsaii - meanaii * meanaii )

  rmsaij = rmsaij / nruns / nstep / 6
  rmsaij = sqrt( rmsaij - meanaij * meanaij )

  write(*,*) 'mean aii, mean aij: ', meanaii, meanaij
  write(*,*) 'rms aii, aij: ', rmsaii, rmsaij
  write(*,*) 'check normalization: ', sum(paii)
  write(*,*) 'check normalization: ', sum(paij)

  paii = paii / binw
  paij = paij / binw

  open(15, file = 'pdfcm.dat')
  do ii = 1, npnt
    tmp = -bnd + (ii-.5) * binw
    write(15,'(15E15.4)') (tmp-meanaii)/rmsaii, paii(ii)*rmsaii, (tmp-meanaij)/rmsaij, paij(ii)*rmsaij
  end do
  close(15)

  write(*,*) 'finished'

end program cmmodelpost      
