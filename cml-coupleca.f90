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
  real(sp), dimension(three, three) :: dwij, dbij, id
  real(sp), dimension(three, three, 2) :: gaij

  integer :: idum, ii, jj, kk, ll, nini, ireal

  real(sp) :: twrite, tmax, time, trace
  real(sp) :: gasdev
  character(80) :: str, strireal

  ! stuff for rk5
  real(sp) :: h1, hmin, hnext, eps
  integer :: nok, nbad
  external :: dcadtcml


  h1 = dt/10.
  hmin = 1.e-20_dp
  eps = 1.e-3_dp

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
        dbij(ii,jj)=gasdev(idum)
      end do
      end do
      aall(:,:,ll) = af * id * sum(dbij * id) + bf * dbij + cf * transpose(dbij)
      gall(:,:,ll) = id

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'cm-aij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
      read(10) gall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .lt. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then
      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'cm-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'-'//strireal(1:len_trim(strireal))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
        write(10) gall
      close(10)
      ll = ll + 1
      twrite = twrite + dtwrite
    end if

    do ii = 1, npnt

      gaij(:,:,1) = gall(:,:,ii)
      gaij(:,:,2) = aall(:,:,ii)
 
      call rkcaller(dcadtcml,gaij,2*nine,time,time+dt,eps,h1,hmin,hnext,nok,nbad)

      do kk = 1, three
      do jj = 1, three
        dbij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do

      ! White noise force
      dwij = af * id * sum(dbij * id) + bf * dbij + cf * transpose(dbij)

      trace = dwij(1,1) + dwij(2,2) + dwij(3,3)
      dwij(1,1) = dwij(1,1) - trace / 3
      dwij(2,2) = dwij(2,2) - trace / 3
      dwij(3,3) = dwij(3,3) - trace / 3

      gaij(:,:,2) = gaij(:,:,2) + dwij


      gall(:,:,ii) = gaij(:,:,1)
      aall(:,:,ii) = gaij(:,:,2)
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program cmlcouple      
