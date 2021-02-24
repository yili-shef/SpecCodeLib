! This version includes stratification

module mprmtr
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 1000

  ! Gaussian random noise forcing term.
  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  ! controls
  real(sp), parameter :: dt = 4.e-4, dtwrite = 0.5 
  ! Integral time scale is 1.

  ! parameters for models
  !real(sp), parameter :: gmma = 0.08
  real(sp), parameter :: gmma = 0.1

end module mprmtr

program cm
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three, npnt) :: aall
  real(sp), dimension(three, three) :: aij, dwij, id 
  real(sp), dimension(three, three) :: aij1, driftaij0, driftaij1

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: twrite, tmax, time, tmp
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
          write(*,*) ' Usage: ./cm-s(d)p.x nini tmax idum ireal'
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
    str = 'cmaij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .le. tmax + dt/2)

    if ( abs(time - twrite) .le. .5 * dt ) then

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'cmaij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time', time, 'aii', sum(aall(1,1,:) + aall(2,2,:)+aall(3,3,:))/npnt

    end if

    do ii = 1, npnt

      do kk = 1, three
      do jj = 1, three
        aij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do

      ! White noise force for velocity gradient
      dwij = af * id * sum(aij * id) + bf * aij + cf * transpose(aij)

      !!!! Test TODO !!!!
      dwij = 0.

      aij = aall(:,:,ii)

      ! order 2 weak convergent predictor-corrector method
      call drift(aij, driftaij0)
      aij1 = aij + driftaij0 * dt + dwij

      call drift(aij1, driftaij1)
      aij1 = aij + .5 * (driftaij1 + driftaij0) * dt + dwij

      call drift(aij1, driftaij1)
      aij = aij + .5 * (driftaij1 + driftaij0) * dt + dwij

      ! remove the trace 
      tmp = -( aij(1,1) + aij(2,2) + aij(3,3) )/3.

      aij(1,1) = aij(1,1) + tmp
      aij(2,2) = aij(2,2) + tmp
      aij(3,3) = aij(3,3) + tmp

      aall(:,:,ii) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cm      

subroutine drift(aij, driftaij)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three) :: aij, driftaij
  real(sp), dimension(three, three) :: a2ij, cinvij, delij
  real(sp) :: tra2, trcinvij

  delij = 0; delij(1,1) = 1; delij(2,2) = 1; delij(3,3) = 1
  !========== drift term for the aij equation

  ! Matrix exponential 
  !cinvij = -gmma * aij
  !tra2 = sqrt( sum ( cinvij * cinvij ) )
  !call taylorch(cinvij, tra2, three)

  !cinvij = matmul(transpose(cinvij), cinvij)
  !trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)

  a2ij = matmul(aij, aij)
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)
 
  !driftaij =  -a2ij + tra2 / trcinvij * cinvij - (trcinvij / 3) * aij 

  !driftaij = -a2ij -  gmma * tra2 * ( aij + transpose(aij) ) / 3 +  tra2 * delij / 3  - aij

  ! testing .....
  driftaij = -a2ij +  tra2 * delij / 3  - (trcinvij / 3) * aij

end subroutine drift
