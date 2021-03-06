!!!!!!!!!!!!!!!!! Based on rcm-12012009.f90 !!!!!!!!!!!
!!! Simply removing the forcing term !!!
program decayrcm
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100000

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 1.e-4, dtwrite = .2 
  ! Integral time scale is 1.

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: aij, id 
  real(sp), dimension(three, three) :: aij1, driftaij0, driftaij1

  integer :: idum, ii, jj, ll, nini, ireal

  real(sp) :: ro, omega, twrite, tmax, time
  real(sp) :: gasdev
  character(80) :: str, strro, strireal


  write(*,*)
  write(*,'('' >>> Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-decay-27102009-s(d)p.x nini tmax ro idum ireal'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        ro: rossby number'
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
  call getarg(3,strro)
  read(strro, '(F10.6)') ro
  call getarg(4, str)
  read(str, '(I20)') idum
  call getarg(5, strireal)
  read(strireal, '(I20)') ireal
  strireal = adjustl(strireal)

  write(*,*) 'tmax = ', tmax
  write(*,*) 'Rossby # = ', ro
  omega = 1./ro

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
    str = 'rcm-decay-'//strro(1:len_trim(strro))//'-'//str(1:len_trim(str)) 
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
      str = 'rcm-decay-'//strro(1:len_trim(strro))//'-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)
      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time', time, 'aii', sum(aall(1,1,:)+aall(2,2,:)+aall(3,3,:))/npnt
    end if

    do ii = 1, npnt

      aij = aall(:,:,ii)

      ! order 2 weak convergent predictor-corrector method
      call drift(aij, omega, driftaij0)
      aij1 = aij + driftaij0 * dt 

      call drift(aij1, omega, driftaij1)
      aij1 = aij + .5 * (driftaij1 + driftaij0) * dt 

      call drift(aij1, omega, driftaij1)
      aij = aij + .5 * (driftaij1 + driftaij0) * dt 

      ! remove the trace 
      ro = -( aij(1,1) + aij(2,2) + aij(3,3) )/3.

      aij(1,1) = aij(1,1) + ro
      aij(2,2) = aij(2,2) + ro
      aij(3,3) = aij(3,3) + ro

      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program decayrcm      

subroutine drift(aij, omega, driftaij)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(sp), parameter :: gmma = 0.08

  real(sp) :: omega
  real(sp), dimension(three, three) :: aij, driftaij

  real(sp), dimension(three, three) :: a2ij, coij, cinvij
  real(sp) :: tra2, trcinvij, o3

  ! drift term for the aij equation

  o3 = aij(2,1) - aij(1,2) ! a(1,2) = d u1 / dx2
  call coriolis(aij, coij)

  ! Matrix exponential 
  cinvij = -gmma * aij
  tra2 = sqrt( sum ( cinvij * cinvij ) )
  call taylorch(cinvij, tra2, three)

  cinvij = matmul(transpose(cinvij), cinvij)
  trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)

  a2ij = matmul(aij, aij)
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)
 
  driftaij =  -a2ij - omega*coij - ((2.*o3*omega-tra2) / trcinvij) * cinvij - (trcinvij / 3) * aij 

end subroutine drift

subroutine coriolis(aij, coij) ! Nondimensionalized, so that omega = 1
  use mconstant
  implicit none

  real(sp), dimension(3,3), intent(in) :: aij
  real(sp), dimension(3,3), intent(out) :: coij

  coij(1,1) = -2. * aij(2,1)
  coij(1,2) = -2. * aij(2,2)
  coij(1,3) = -2. * aij(2,3)
  coij(2,1) =  2. * aij(1,1)
  coij(2,2) =  2. * aij(1,2)
  coij(2,3) =  2. * aij(1,3)
  coij(3,1) =  0.
  coij(3,2) =  0.
  coij(3,3) =  0.


end subroutine coriolis
