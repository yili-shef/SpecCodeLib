program cmmodel
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100000

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 1.e-4, dtwrite = .5 
  ! Integral time scale is 1.

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: aij, dwij, id
  real(sp), dimension(three, three) :: aijbar, driftaij0, driftaij1

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: twrite, tmax, time, trace
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
          write(*,*) ' Usage: ./cmmodel-a.x nini tmax idum ireal'
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

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'cm-aij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
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
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite
    end if

    do ii = 1, npnt

      do kk = 1, three
      do jj = 1, three
        aij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do

      ! White noise force
      dwij = af * id * sum(aij * id) + bf * aij + cf * transpose(aij)

      aij = aall(:,:,ii)

      ! order 2 weak convergent method
      call drift(aij, driftaij0)
      aijbar = aij + driftaij0 * dt + dwij
      call drift(aijbar, driftaij1)
      aijbar = aij + .5 * (driftaij1 + driftaij0) * dt + dwij
      call drift(aijbar, driftaij1)
      aij = aij + .5 * (driftaij1 + driftaij0) * dt + dwij

      ! removing numerical error in trace
      trace = - (aij(1,1) + aij(2,2) + aij(3,3)) / 3.
      aij(1,1) = aij(1,1) + trace
      aij(2,2) = aij(2,2) + trace
      aij(3,3) = aij(3,3) + trace

      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cmmodel      

subroutine drift(aij, driftaij)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(sp), parameter :: gam = 0.01
  real(sp), dimension(three, three), intent(in) :: aij
  real(sp), dimension(three, three), intent(out) ::  driftaij

  real(sp), dimension(three, three) :: a2ij
  real(sp) :: tra2, trcinv

  ! drift term for the aij equation

  driftaij = - gam * transpose(aij)

  a2ij = matmul( transpose(driftaij), driftaij )
  trcinv = sqrt(a2ij(1,1) + a2ij(2,2) + a2ij(3,3))
  call taylorch(driftaij, trcinv, three)

  driftaij = matmul(driftaij, transpose(driftaij))
  trcinv = driftaij(1,1) + driftaij(2,2) + driftaij(3,3)

  a2ij = matmul(aij, aij) 
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)

  driftaij = - a2ij + (tra2 / trcinv) * driftaij - trcinv / 3. * aij

end subroutine drift
