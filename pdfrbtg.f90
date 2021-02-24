program pdfrbetatgamma
  use mconstant
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, nthresh
  integer(8) :: numpnt

  real, allocatable, dimension(:,:,:) :: rbeta, tgamma
  complex, allocatable, dimension (:,:,:) :: dissrate

  integer, parameter :: npnt = 160, nmean = 40
  real(dp), dimension (npnt) :: prbtg

  real(dp) :: mrbtg, mdissrate, meandiss

  real :: binw, bound

  real :: btmp, gtmp, ss
  character(80) :: strdiss, strthresh, strrij, strtij, strtmp, strtmp1, strtmp2


  write(*,*) 
  write(*,'(''>>> Statistics of beta of r and gamma of tau <<<'')')
  write(*,*)

  ll = iargc()
  if ( ll .ne. 5 ) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfrbetatgamma.x nx nthresh revfilelist tevfilelist dissratefilelist'
          write(*,*) '        nx: resolution'
          write(*,*) '        nthresh: threshold of dissrate = nthresh * mdissrate'
          write(*,*) '        revflielist: rij ev data file list'
          write(*,*) '        tevfilelist: tau eigenvalue data file list'
          write(*,*) '        dissratefilelist: dissrate data file list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,strdiss)
  read(strdiss, '(I20)') nx

  ! nthresh
  call getarg(2,strthresh)
  read(strthresh, '(I20)') nthresh
  strthresh = adjustl(strthresh)

  ! file list for rij ev data
  call getarg(3,strrij) 
  strrij = adjustl(strrij)

  ! file list for tauij ev data
  call getarg(4, strtij)
  strtij = adjustl(strtij)

  ! file list for dissipation data
  call getarg(5,strdiss)
  strdiss = adjustl(strdiss)

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1

  allocate(rbeta(nx,ny,nz), tgamma(nx,ny,nz))
  allocate(dissrate(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  ! strdiss is the list of dissrate data
  open(21, file = strdiss( 1 : len_trim( strdiss ) )//'.list')

    read(21,*) strtmp

    open(10, file = './out/'//strtmp(1 : len_trim(strtmp)), form = 'unformatted')
      read(10) dissrate
    close(10)

    ! Estimate the mean dissipation using the first data file
    mdissrate = sum( real( dissrate(1:lx,:,:) ) ) + sum( aimag( &
                dissrate(1:lx,:,:) ) ) 
    mdissrate = mdissrate / (nx * ny * nz)

  close(21)

  bound = nmean * mdissrate
  binw = 2. * bound / npnt

  write(*,*) 'Estimated mean dissipation: ', mdissrate

  open(20, file = strdiss(1:len_trim(strdiss))//'.list')
  open(21, file = strrij(1:len_trim(strrij))//'.list')
  open(22, file = strtij(1:len_trim(strtij))//'.list')
  

  numpnt = 0
  prbtg = 0.d0
  meandiss = 0.d0
  mrbtg = 0.d0
  do while ( .not. eof(20) )
    read(20,*) strtmp
    read(21,*) strtmp1
    read(22,*) strtmp2

    open(10, file = './out/'//strtmp(1 : len_trim(strtmp)), form = 'unformatted')
      read(10) dissrate
    close(10)

    open(10, file = './out/'//strtmp1(1 : len_trim(strtmp1)), form = 'unformatted')
      read(10) 
      read(10) rbeta
    close(10)

    open(10, file = './out/'//strtmp2(1 : len_trim(strtmp2)), form = 'unformatted')
      read(10) 
      read(10) 
      read(10) tgamma
    close(10)

    write(*,*) 'after reading data files'

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      
      if ( mod(ii, 2) .eq. 1) then
              ss = real( dissrate( (ii+1)/2, jj, kk ) )
      else 
              ss = aimag( dissrate( ii/2, jj,kk ) )
      end if

      if ( ss .ge. nthresh * mdissrate) then  
      ! conditioned on higher dissipation

        meandiss = meandiss + ss

        numpnt = numpnt + 1

        btmp = rbeta(ii,jj,kk)
        gtmp = tgamma(ii,jj,kk)
       
        ss = 3. * btmp * gtmp

        mrbtg = mrbtg + ss
       
        ll = floor( (ss + bound) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) prbtg(ll) = prbtg(ll) + 1

      end if

    end do
    end do
    end do

  end do

  write(*,*) 'Number of points: ', numpnt

  meandiss = meandiss / numpnt
  mrbtg = mrbtg / numpnt

  write(*,*) 'Mean sgs diss: ', meandiss
  write(*,*) 'eigenvalue approximation: ', mrbtg

  prbtg = prbtg / numpnt 
  write(*,*) 'check prbtg: ', sum(prbtg)
  prbtg = prbtg / binw

  strtmp2 = 'prbtg-'//strdiss(1:len_trim(strdiss))//strthresh(1:len_trim(strthresh))//'.dat'
  open(10, file = strtmp2(1:len_trim(strtmp2)) )
    do ii = 1, npnt
      write(10,*) - bound + (ii - .5) * binw, prbtg(ii)
    end do
  close(10)

  deallocate(rbeta,tgamma)
  deallocate(dissrate)

  write(*,*) 'finished'

end program pdfrbetatgamma 
