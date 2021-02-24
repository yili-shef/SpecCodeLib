program uo
  use mconstant
  implicit none

  integer, parameter :: three = 3

  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  real(sp), parameter :: dt = 1.e-4, dtwrite = .1 
  ! Integral time scale is 1.

  real(sp), allocatable, dimension(:, :, :) :: aall, ball
  real(sp), dimension(three, three) :: aij, bij, dwij, id, aijo, bijo
  real(sp), dimension(three, three) :: aijbar, driftaij0, driftaij1

  integer :: idum, ii, jj, kk, ll, nini, npnt

  real(sp) :: twrite, tmax, time, trace
  real(sp) :: gasdev
  character(80) :: str, strpnt


  write(*,*)
  write(*,'('' >>> Ornstein-Uhlenbeck process for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./uo.x nini tmax idum npnt'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        idum: seed for random number'
          write(*,*) '        npnt: number of realizations'
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
  call getarg(4, strpnt)
  read(strpnt, '(I20)') npnt
  strpnt = adjustl(strpnt)

  allocate( aall(three,three,npnt) )
  allocate( ball(three,three,npnt) )

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
      ball(:,:,ll) = id

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'uo-aij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//trim(strpnt)//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'uo-bij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//trim(strpnt)//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) ball
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then
      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'uo-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//trim(strpnt)//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'uo-bij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//trim(strpnt)//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) ball
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

      aijo = aall(:,:,ii)
      bijo = ball(:,:,ii)

      ! order 2 weak convergent method
      call drift(aijo, driftaij0)
      aijbar = aijo + driftaij0 * dt + dwij
      call drift(aijbar, driftaij1)
      aijbar = aijo + .5 * (driftaij1 + driftaij0) * dt + dwij
      call drift(aijbar, driftaij1)
      aij = aijo + .5 * (driftaij1 + driftaij0) * dt + dwij

      ! removing numerical error in trace
      trace = - (aij(1,1) + aij(2,2) + aij(3,3)) / 3.
      aij(1,1) = aij(1,1) + trace
      aij(2,2) = aij(2,2) + trace
      aij(3,3) = aij(3,3) + trace

      ! solving for bij
      call calbij(aijo, aij, bijo, bij, dt)

      aall(:,:,ii) = aij
      ball(:,:,ii) = bij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'uo.x finished'
  deallocate(aall)
  
end program uo      

subroutine drift(aij, driftaij)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(sp), dimension(three, three), intent(in) :: aij
  real(sp), dimension(three, three), intent(out) ::  driftaij

  ! drift term for the aij equation

  driftaij = - aij

end subroutine drift


subroutine calbij(aijo, aij, bijo, bij, dt)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(sp), intent(in) :: dt
  real(sp), dimension(three, three), intent(in) :: aij, aijo, bijo
  real(sp), dimension(three, three), intent(out) ::  bij

  real(sp), dimension(three,three) :: b1, b2, b3

  b1 = bijo +  .5_sp * dt * matmul( aijo, bijo )
  b2 = bijo + .25_sp * dt * matmul( aijo + aij, b1 )
  b3 = bijo +  .5_sp * dt * matmul( aijo + aij, b2 )
  bij = bijo + (dt/6._sp) * ( matmul(aijo, bijo) + matmul(aijo + aij, b1 + b2) + matmul(aij, b3) )

end subroutine calbij
