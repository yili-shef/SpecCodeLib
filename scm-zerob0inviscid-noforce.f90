! This version includes stratification

module mprmtr
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100000

  ! Gaussian random noise forcing term.
  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  ! controls
  real(sp), parameter :: dt = 2.e-4, dtwrite = 0.02 
  ! Integral time scale is 1.

  ! parameters for models
  real(sp), parameter :: gmma = 0.1
  real(sp), parameter :: timetheta = 0.4

end module mprmtr

program scm
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: aij, id 
  real(sp), dimension(three, three) :: aij1, driftaij0, driftaij1

  real(sp), dimension(three, npnt) :: gall
  real(sp), dimension(three) :: gi, gi1, driftgi0, driftgi1

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: ro, twrite, tmax, time, beta
  real(sp) :: gasdev
  character(80) :: str, strireal


  write(*,*)
  write(*,'('' >>> Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./scm-zerob0.x nini tmax beta idum ireal'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        beta: buoyancy force '
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

  call getarg(3,str)
  read(str, '(F10.6)') beta ! Buoyancy force beta g gamma in Claude 2001 notation.

  call getarg(4, str)
  read(str, '(I20)') idum

  call getarg(5, strireal)
  read(strireal, '(I20)') ireal
  strireal = adjustl(strireal)

  write(*,*) 'tmax = ', tmax
  write(*,*) 'dtwrite = ', dtwrite
  write(*,*) 'Buoyancy force = ', beta

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

      do ii = 1, three
      gall(ii,ll) = 0.
      end do

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'scm-aij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'scm-gi-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) gall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'scm-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'scm-gi-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) gall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time', time, 'aii', sum(aall(1,1,:) + aall(2,2,:)+aall(3,3,:))/npnt

    end if

    do ii = 1, npnt

      aij = aall(:,:,ii)
      gi = gall(:,ii)

      ! order 2 weak convergent predictor-corrector method
      call drift(aij, gi, beta, driftaij0, driftgi0)
      aij1 = aij + driftaij0 * dt 
      gi1 = gi + driftgi0 * dt

      call drift(aij1, gi1, beta, driftaij1, driftgi1)
      aij1 = aij + .5 * (driftaij1 + driftaij0) * dt 
      gi1 = gi + .5 * (driftgi1 + driftgi0) * dt 

      call drift(aij1, gi1, beta, driftaij1, driftgi1)
      aij = aij + .5 * (driftaij1 + driftaij0) * dt
      gi = gi + .5 * (driftgi1 + driftgi0) * dt 

      ! remove the trace 
      ro = -( aij(1,1) + aij(2,2) + aij(3,3) )/3.

      aij(1,1) = aij(1,1) + ro
      aij(2,2) = aij(2,2) + ro
      aij(3,3) = aij(3,3) + ro

      aall(:,:,ii) = aij
      gall(:,ii) = gi
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program scm      

subroutine drift(aij, gi, beta, driftaij, driftgi)
  use mconstant
  use mprmtr
  implicit none

  real(sp) :: beta
  real(sp), dimension(three, three) :: aij, driftaij
  real(sp), dimension(three) :: gi, driftgi

  real(sp), dimension(three, three) :: a2ij, cinvij
  real(sp) :: tra2, trcinvij

  ! drift term for the aij equation

  ! Matrix exponential 
  cinvij = -gmma * aij
  tra2 = sqrt( sum ( cinvij * cinvij ) )
  call taylorch(cinvij, tra2, three)

  cinvij = matmul(transpose(cinvij), cinvij)
  trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)

  a2ij = matmul(aij, aij)
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)
 
  driftaij =  -a2ij - ((-tra2 + gi(3)) / trcinvij) * cinvij 
  driftaij(3,:) = driftaij(3,:) + gi ! stratification is in the vertical direction

  driftgi = - beta * aij(3,:) - matmul(transpose(aij), gi) 

end subroutine drift
