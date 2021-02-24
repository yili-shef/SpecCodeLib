program cml
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 400000, nine = three * three

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 2.e-4, dtwrite = .5, gam = 0.08
  ! Integral time scale is 1.

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: dwij, aij, id, expaij, a2ij

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: twrite, tmax, time, tra2, trinvcij
  real(sp) :: gasdev
  character(80) :: str0, str, strireal


  id = 0.; id(1,1) = 1.; id(2,2) = 1.; id(3,3) = 1.

  write(*,*)
  write(*,'('' >>> Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cml-12012009.x nini tmax idum ireal'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        idum: seed for random number'
          write(*,*) '        ireal: the ith realization'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,str)
  read(str, '(I20)') nini
  call getarg(2,str)
  read(str, '(F10.6)') tmax
  call getarg(3, str)
  read(str, '(I20)') idum
  call getarg(4, strireal)
  read(strireal, '(I20)') ireal
  strireal = adjustl(strireal)

  write(*,*) 'tmax = ', tmax

  id = 0.; id(1,1) = 1.; id(2,2) = 1.; id(3,3) = 1.

  if ( nini .eq. 0 ) then

    idum = - idum
    write(*,*) 'idum = ', idum
    do ll = 1, npnt

      do ii = 1, three
      do jj = 1, three
        aij(ii,jj)=gasdev(idum)
      end do
      end do
      aall(:,:,ll) = af * id * sum(aij * id) + bf * aij + cf * transpose(aij)

    end do

  else

    write(str0, '(I20)') nini
    str0 = adjustl(str0)

    str = 'cm-aij-'//str0(1:len_trim(str0)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .lt. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then
      write(str0, '(I20)') ll
      str0 = adjustl(str0)

      str = 'cm-aij-'//str0(1:len_trim(str0)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite
    end if

    do ii = 1, npnt
      
      ! diffusion coefficient
      do kk = 1, three
      do jj = 1, three
        aij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do

      ! force term
      dwij = af * id * sum(aij * id) + bf * aij + cf * transpose(aij)

      aij = aall(:,:,ii)
 
      expaij = -gam * aij
      tra2 = sqrt( sum(expaij * expaij) )
      call taylorch(expaij, tra2, three)
      expaij = matmul( transpose(expaij), expaij )
      trinvcij = expaij(1,1) + expaij(2,2) + expaij(3,3)
  
      a2ij = matmul(aij, aij) 
      tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3) 

      a2ij = (-a2ij + (tra2 / trinvcij) * expaij - (trinvcij/3) * aij) * dt + (trinvcij) * dwij

      aij = aij + a2ij

      ! order 2 weak convergent predictor-corrector method
      ! call drift(aij, driftaij0)
      ! aij1 = aij + driftaij0 * dt + dwij

      ! call drift(aij1, driftaij1)
      ! aij1 = aij + .5 * (driftaij1 + driftaij0) * dt + dwij

      ! call drift(aij1, driftaij1)
      ! aij = aij + .5 * (driftaij1 + driftaij0) * dt + dwij

      ! removing numerical error in trace
      tra2 = - (aij(1,1) + aij(2,2) + aij(3,3)) / 3.
      aij(1,1) = aij(1,1) + tra2
      aij(2,2) = aij(2,2) + tra2
      aij(3,3) = aij(3,3) + tra2

      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cml      

