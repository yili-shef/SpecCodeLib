program uo
  use mconstant
  implicit none

  integer, parameter :: three = 3

  real(sp), parameter :: dt = 1.e-4, dtwrite = .1 
  ! Integral time scale is 1.

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: lwork = 26 * three, liwork = 10 * three
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*three), iignore, nfound, info
  real(dp) :: evbij(three), evtrbij(three,three), mbijeval, mbijevbe, mbijevgm
  ! ----------------------------------------------------------------------------------

  real(sp), allocatable, dimension(:, :, :) :: ball
  real(sp), dimension(three, three) :: bij

  integer :: idum, ii, jj, kk, ll, nini, npnt, nlast

  real(sp) :: twrite, tmax, time, trace
  character(80) :: str, strpnt


  write(*,*)
  write(*,'('' >>> Ornstein-Uhlenbeck process for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./postuo.x nini nlast npnt'
          write(*,*) '        nini:  first file number'
          write(*,*) '        nlast: last file number'
          write(*,*) '        npnt: number of realizations'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,str)
  read(str, '(I20)') nini
  call getarg(2,str)
  read(str, '(I20)') nlast
  call getarg(3, strpnt)
  read(strpnt, '(I20)') npnt
  strpnt = adjustl(strpnt)

  allocate( ball(three,three,npnt) )

  open(30, file = 'bijev.dat')

  do ii = nini, nlast

    write(str, '(I20)') ii
    str = adjustl(str)
    str = 'uo-bij-'//str(1:len_trim(str)) 
    str = str(1:len_trim(str))//'-'//trim(strpnt)//'.data'
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) ball
    close(10)

    mbijeval = 0._dp
    mbijevbe = 0._dp
    mbijevgm = 0._dp
    do jj = 1, npnt

      bij = ball(:,:,jj)
      bij = matmul(bij, transpose(bij))

      ! The eigenvalues are in ascending order: smallest first.
      call dsyevr("V", "A", "U", three, bij, three, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evbij, evtrbij, three, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      mbijeval = mbijeval + log(evbij(3))
      mbijevbe = mbijevbe + log(evbij(2))
      mbijevgm = mbijevgm + log(evbij(1))

    end do
    mbijeval = mbijeval / npnt
    mbijevbe = mbijevbe / npnt
    mbijevgm = mbijevgm / npnt

    write(30, *) ii, mbijeval, mbijevbe, mbijevgm

  end do
  close(30)

  deallocate(ball)
  write(*,*) 'postuo.x finished'
end 

