program rcmmore
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: sij, aij

  real(dp) :: meana12, meana21
  real(dp) :: rmsa12, rmsa21
  real(dp) :: skewa12, skewa21
  real(dp) :: c1, c2, c3, meano3, rmso3
  real(dp) :: oos, rmsoos
  real(dp) :: lnaijaij, rmslnaijaij
  real(dp) :: aijaij, rmsaijaij

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
          write(*,*) ' Usage: ./rcm-more-means.x nfiles nstart ro nreal ifirst'
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
  write(15, *)'ro c3'

  do ll = nstart, nfiles + nstart -1
    write(*,*) 'file ', ll
    write(str, '(I20)') ll
    str = adjustl(str)
    str = 'cm-aij-ro='//strro(1:len_trim(strro))//'-'//str(1:len_trim(str))//'-'


    meana12 = 0.d0; rmsa12 = 0.d0
    meano3 = 0.d0; rmso3 = 0.d0
    c1 = 0.d0; c2 = 0.d0
    meana21 = 0.d0; rmsa21 = 0.d0
    skewa12 = 0.d0; skewa21 = 0.d0
    oos = 0.d0; rmsoos = 0.d0
    lnaijaij = 0.d0; rmslnaijaij = 0.d0
    aijaij = 0.d0; rmsaijaij = 0.d0

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

        ! correlation <a12^2 a21> 
        meano3 = meano3 + wz
        rmso3 = rmso3 + wz * wz
        c1 = c1 + ( aij(1,2) * aij(1,2) * aij(2,1) )
        c2 = c2 + ( aij(1,2) * aij(2,1) )
        meana12 = meana12 + aij(1,2)
        meana21 = meana21 + aij(2,1)
        rmsa12 = rmsa12 +  aij(1,2) * aij(1,2) 
        rmsa21 = rmsa21 +  aij(2,1) * aij(2,1) 
        skewa12 = skewa12 + aij(1,2)**3 
        skewa21 = skewa21 + aij(2,1)**3 

        ! vortex stretching
        tmp =   wx * wx * sij(1,1) + wy * wy * sij(2,2)  &
              + wz * wz * sij(3,3) + 2. * ( wx * wy * sij(1,2) &
              + wx * wz * sij(1,3) + wy * wz * sij(2,3) )
        oos = oos + tmp
        rmsoos = rmsoos + tmp * tmp

        ! log of pseudodissipation lnaijaij
        sij = matmul(aij, transpose(aij))
        tmp = sij(1,1) + sij(2,2) + sij(3,3) 

        aijaij = aijaij + tmp
        rmsaijaij = rmsaijaij + tmp * tmp

        tmp = log( sij(1,1) + sij(2,2) + sij(3,3) )

        lnaijaij = lnaijaij + tmp
        rmslnaijaij = rmslnaijaij + tmp * tmp


      end do
  
    end do

    tmp = 1. / nreal / npnt
 
    meana12 = meana12 * tmp
    meana21 = meana21 * tmp
    c1 = c1 * tmp
    c2 = c2 * tmp
    rmsa12 = rmsa12 * tmp

    meano3 = meano3 * tmp
    rmso3 = rmso3 * tmp
    rmso3 = sqrt( rmso3 - meano3**2 )

    ! correlation <a12^2 a21>
    c3 = c1 - 2 * c2 * meana12 + 2 * meana12**2 * meana21 - rmsa12 * meana21
    ! c3 = c3 / rmso3**3

    oos = oos * tmp
    rmsoos = sqrt(rmsoos * tmp - oos * oos)

    lnaijaij = lnaijaij * tmp
    rmslnaijaij = sqrt(rmslnaijaij * tmp - lnaijaij * lnaijaij)
 
    aijaij = aijaij * tmp
    rmsaijaij = sqrt(rmsaijaij * tmp - aijaij * aijaij)

    write(15, '(15E13.4)') c3, oos, rmsoos, lnaijaij, rmslnaijaij, aijaij, rmsaijaij

  end do
  close(15)

  write(*,*) 'finished'

end program rcmmore 
