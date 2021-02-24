
! The stochastic model for vorticity in Zybin etal 
module mprmtr
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100

  ! controls
  real(sp), parameter :: dt = 1.e-5_sp, dtwrite = 0.001_sp
  ! Integral time scale is 1.

  integer, dimension(three, three), parameter :: delij = (/(/1,0,0/),(/0,1,0/),(/0,0,1/)/)

end module mprmtr

program vortmod
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, npnt) :: omall, doall

  real(sp), dimension(three) :: omi, doi, omi1, doi1
  real(sp), dimension(three) :: uomi, romi, phidoi
  real(sp), dimension(three,three) :: dwij
  real(sp), dimension(three, three, three) :: bijm, tmpval1, tmpval2

  integer :: idum, ii, jj, kk, ll, mm, nn, rr, nini 

  real(sp) :: twrite, tmax, time, gasdev
  character(80) :: str


  write(*,*)
  write(*,'('' >>> Vorticity Model <<< '')')
  write(*,*) 
  ll = iarg()
  if (iarg() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./vortmod.x nini tmax idum'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        idum: seed for random number'
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


  write(*,*) 'tmax = ', tmax
  write(*,*) 'dtwrite = ', dtwrite

  if ( nini .eq. 0 ) then

    idum = - idum
    write(*,*) 'idum = ', idum
    do ll = 1, npnt

      do ii = 1, three
        omi(ii) = gasdev(idum)
        doi(ii) = gasdev(idum)
      end do
      omall(:,ll) = omi
      doall(:,ll) = doi

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'vort-omi-'//str(1:len_trim(str))//'.data' 
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) omall
    close(10)

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'vort-doi-'//str(1:len_trim(str))//'.data' 
    open(10, file = str(1:len_trim(str)), form = 'binary')
    read(10) doall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini

  open(15, file = 'omiseries.dat')
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'vort-omi-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) omall
      close(10)

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'vort-doi-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) doall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time',  time, 'max omi', maxval(abs(omall))
      write(15,'(10E12.3)') time, omall(1,1), omall(2,1), omall(3,1)

    end if

    do rr = 1, npnt

      ! White noise
      do kk = 1, three
      do jj = 1, three
        dwij(jj,kk)=gasdev(idum) * sqrt(2._sp*dt)
      end do
      end do

      omi = omall(:,rr)
      doi = doall(:,rr)

      ! order 2 weak convergent predictor-corrector method

      call diffusion(omi, bijm)

      ! 
      omi1 = omi + doi * dt        ! Drift term for omi. The diffusion term for omi is zero
      doi1 = doi 
      do nn = 1, three
      do mm = 1, three
        doi1 = doi1 + bijm(:,mm,nn) * dwij(mm,nn) !Diffusion term for doi
      end do
      end do

      ! Calculate phi
      ! ---- For our problem the general formulas can be simplied to the following ----

      uomi = omi
      romi = omi + doi * dt

      call diffusion(uomi, tmpval1)
      tmpval1 = (three * three - 1) * 2 * ( tmpval1 - bijm ) / sqrt(dt)

      call diffusion(romi, tmpval2)

      phidoi = 0
      do nn = 1, three
      do mm = 1, three
        phidoi = phidoi + ( 2 * tmpval2(:,mm,nn) + 2 * bijm(:,mm,nn) + tmpval1(:,mm,nn) ) &
                          * dwij(mm,nn) 
      end do
      end do
      phidoi = phidoi / 4
      ! ------------------ For the full formula see end of the program --------------------
      ! -----------------------------------------------------------------------------------

      omi1 = omi + (doi1 + doi ) * dt / 2
      doi1 = doi + phidoi

      omi  = omi + (doi1 + doi ) * dt / 2
      doi  = doi + + phidoi

      omall(:,rr) = omi
      doall(:,rr) = doi
 
    end do
 
    time = time + dt
  end do
  close(15)

  write(*,*) 'finished'
  
end program vortmod      

subroutine diffusion(omi, bijm)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three,three,three), intent(out) :: bijm
  real(sp), dimension(three), intent(in) :: omi

  integer :: ii, jj, mm

  
  do mm = 1, three
  do jj = 1, three
  do ii = 1, three

   ! The following should be correct now.
    bijm(ii,jj,mm) = sqrt(2._sp)/ 2 * ( delij(ii,jj) * omi(mm)  &
                                      + delij(ii,mm) * omi(jj) )
  end do
  end do
  end do
  
  
end subroutine diffusion 
