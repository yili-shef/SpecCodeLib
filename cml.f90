program cml
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100000

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 1.e-4, gmma = 0.08, dtwrite = .5 
  ! Integral time scale is 1.

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: aij, dwij, id, dbij, cinvij

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: anorm, trcinvij, traa, trace
  real(sp) :: twrite, tmax, time
  real(sp) :: gasdev
  character(80) :: str, strireal


  write(*,*)
  write(*,'('' >>> Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cml.x nini tmax idum ireal'
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
        dbij(ii,jj)=gasdev(idum)
      end do
      end do
      aall(:,:,ll) = af * id * sum(dbij * id) + bf * dbij + cf * transpose(dbij)

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'cm-aij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'unformatted')
      read(10) aall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then
      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'cm-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'unformatted')
        write(10) aall
      close(10)
      ll = ll + 1
      twrite = twrite + dtwrite
    end if

    do ii = 1, npnt

      aij = aall(:,:,ii)
 
      do kk = 1, three
      do jj = 1, three
        dbij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do

      ! White noise force
      dwij = af * id * sum(dbij * id) + bf * dbij + cf * transpose(dbij)

      ! nonlinear term
      dbij = matmul(aij, aij)
      traa = dbij(1,1) + dbij(2,2) + dbij(3,3)
 
      ! Matrix exponential: taking into account changes of aij in time tau due to
      ! the nonlinear term

      cinvij = gmma * aij - gmma*gmma/2 * matmul(aij,aij)
      cinvij = - transpose(cinvij)

      trace = cinvij(1,1) + cinvij(2,2) + cinvij(3,3) ! trace of the matrix
      cinvij = cinvij - trace * id / 3.        ! remove the trace

      anorm = sqrt( sum ( matmul(cinvij, transpose(cinvij)) * id ) )
      call taylorch(cinvij, anorm, three) ! cinvij has to be traceless
      cinvij = cinvij * exp(trace/3.)

      ! Cauchy-Green tensor
      cinvij = matmul(cinvij, transpose(cinvij))

      trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)

      ! Simple Euler forward for daij
      dbij = dt * ( -dbij + traa / trcinvij * cinvij - (trcinvij / 3) * aij ) + dwij
 
      aij = aij + dbij

      ! remove the trace 
      traa = -( aij(1,1) + aij(2,2) + aij(3,3) )/3.

      aij(1,1) = aij(1,1) + traa
      aij(2,2) = aij(2,2) + traa
      aij(3,3) = aij(3,3) + traa

      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cml      
