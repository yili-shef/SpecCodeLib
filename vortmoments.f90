program vortmoments
  use mconstant
  implicit none
  
  integer,  parameter :: three = 3, npnt = 100
  real(sp), parameter :: dtw = 0.001_sp

  real(sp), dimension(three, npnt) :: omall
  real(sp), dimension(npnt) :: ommag

  real(sp) :: time, om2, om4, om6
  character(80) :: str, filelst


  write(*,*) 
  write(*,'('' >>> Moments of vorticity in Zybin model <<< '')')
  write(*,*)
  if ( iarg() .ne. 1 ) then
      write(*,*) 
      write(*,*) ' >>> Wrong number of arguments <<<'
      write(*,*)
      write(*,*) ' Usage: ./vortmoments.x filelist'
      write(*,*) 
      stop 'Stopped'
  end if

  call getarg(1,filelst)
  filelst = adjustl(filelst)

  open(30, file = filelst(1:len_trim(filelst))//'.list')

  time = 0.
  open(15, file = 'vmtimeseries'//filelst(1:len_trim(filelst))//'.dat')

  do while ( .not. eof(30))

    read(30,*) str
    str = adjustl(str)
    write(*,*) 'reading data file:', str(1:len_trim(str))

    open(20, file = str(1:len_trim(str)), form = 'binary')
      read(20) omall
    close(20)

    ommag = omall(1,:) * omall(1,:) + omall(2,:) * omall(2,:) + omall(3,:) * omall(3,:)

    om2 = sum(ommag)/npnt

    omall(1,:) = ommag * ommag
    om4 = sum( omall(1,:) ) / npnt

    omall(1,:) = omall(1,:) * ommag
    om6 = sum( omall(1,:) ) / npnt

    ! output
    write(15,*)  time, om2, om4, om6

    time = time + dtw

  end do

  close(15)
  close(30)

  write(*,*) 'Finished'

end program vortmoments
