program cmlcouple
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100000, nine = three * three

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 1.e-4, dtwrite = .5 
  ! Integral time scale is 1.

  real(sp), dimension(three, three, npnt) :: aall, gall
  real(sp), dimension(three, three) :: dwij, aij, gij, id
  real(sp), dimension(three, three) :: gbaraij, gbargij
  real(sp), dimension(three, three) :: driftaij0, driftgij0
  real(sp), dimension(three, three) :: driftaij1, driftgij1

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: twrite, tmax, time, trace
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
          write(*,*) ' Usage: ./cml-coupleca.x nini tmax idum ireal'
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
      gall(:,:,ll) = id

    end do

  else

    write(str0, '(I20)') nini
    str0 = adjustl(str0)

    str = 'cm-aij-'//str0(1:len_trim(str0)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

    str = 'cm-gij-'//str0(1:len_trim(str0)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) gall
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

      str = 'cm-gij-'//str0(1:len_trim(str0)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) gall
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
      gij = gall(:,:,ii)
 
      ! order 2 weak convergent predictor-corrector method
      call drift(aij, gij, driftaij0, driftgij0)
      gbaraij = aij + driftaij0 * dt + dwij
      gbargij = gij + driftgij0 * dt

      call drift(gbaraij, gbargij, driftaij1, driftgij1)
      gbaraij = aij + .5 * (driftaij1 + driftaij0) * dt + dwij
      gbargij = gij + .5 * (driftgij1 + driftgij0) * dt

      call drift(gbaraij, gbargij, driftaij1, driftgij1)
      aij = aij + .5 * (driftaij1 + driftaij0) * dt + dwij
      gij = gij + .5 * (driftgij1 + driftgij0) * dt

      ! removing numerical error in trace
      trace = - (aij(1,1) + aij(2,2) + aij(3,3)) / 3.
      aij(1,1) = aij(1,1) + trace
      aij(2,2) = aij(2,2) + trace
      aij(3,3) = aij(3,3) + trace

      gall(:,:,ii) = gij
      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cmlcouple      

subroutine drift(aij, gij, driftaij, driftgij)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(sp), parameter :: gam = 0.08
  real(sp), dimension(three, three) :: aij, gij, driftaij, driftgij

  real(sp), dimension(three, three) :: a2ij
  real(sp) :: tra2, trcinv

  ! drift term for the aij equation

  a2ij = matmul(aij, aij) 
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)

  driftgij = matmul( transpose(gij), gij )
  trcinv = driftgij(1,1) + driftgij(2,2) + driftgij(3,3)

  driftaij = - a2ij + (tra2 / trcinv) * driftgij - trcinv / 3. * aij

  ! drift term for the gij equation

  a2ij = gam * aij
  driftgij = matmul( transpose(a2ij), a2ij )
  tra2 = sqrt( driftgij(1,1) + driftgij(2,2) + driftgij(3,3) )
  call taylorch(a2ij, tra2, three)

  driftgij = - gam * aij
  call taylorch(driftgij, tra2, three)
  driftgij = matmul( transpose(driftgij), driftgij )
  tra2 = driftgij(1,1) + driftgij(2,2) + driftgij(3,3)

  ! use exp(-gam*aij) as damping
  driftgij = - matmul(gij, aij) + (tra2/3) * (gij - matmul(gij, matmul(a2ij,gij) ) ) 

end subroutine drift

