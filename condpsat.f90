program condpsat
  use mconstant
  implicit none

  !real, parameter :: mdissrate = 0.2014974 ! delta = 8dx

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, nthresh
  integer(8) :: numpnt

  real, allocatable, dimension(:,:,:) :: alpha, beta, gmma
  complex, allocatable, dimension(:,:,:) :: dissrate

  integer, parameter :: npnt = 100
  real,    parameter :: binw = 2. / npnt
  real(dp), dimension (npnt) :: psatcond

  real :: atmp, btmp, gtmp, ss, mdissrate
  character(80) :: str, str1, str2, str3, str5

  write(*,*) 
  write(*,'(''>>> Conditional PDF of s asterik <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./condpsat.x nx nthresh evfilelist dissratefilelist'
          write(*,*) '        nx: resolution'
          write(*,*) '        nthresh: threshold of dissrate = nthresh * mdissrate'
          write(*,*) '        evfilelist: eigenvalue data file list'
          write(*,*) '        dissratefilelist: dissrate data file list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if 

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! nthresh
  call getarg(2,str5)
  read(str5, '(I20)') nthresh
  str5 = adjustl(str5)

  ! file list for ev data
  call getarg(3,str1) 
  str1 = adjustl(str1)

  ! file list for dissipation data
  call getarg(4,str)
  str = adjustl(str)

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1

  allocate(alpha(nx,ny,nz), beta(nx,ny,nz), gmma(nx,ny,nz))
  allocate(dissrate(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  ! str is the list of dissrate data
  open(21, file = str (1:len_trim(str ))//'.list')
    read(21,*) str3
    open(10, file = './out/'//str3(1 : len_trim(str3)), form = 'unformatted')
      read(10) dissrate
    close(10)
    mdissrate = sum( real( dissrate(1:lx,:,:) ) ) + sum( aimag( &
                dissrate(1:lx,:,:) ) ) 
    mdissrate = mdissrate / (nx * ny * nz)
  close(21)

  write(*,*) 'Estimated mean dissipation: ', mdissrate

  ! str1 is the list of eigenvalues
  ! str  is the list of dissipation data
  open(20, file = str1(1:len_trim(str1))//'.list')
  open(21, file = str (1:len_trim(str ))//'.list')

  numpnt = 0
  psatcond = 0.d0
  do while ( .not. eof(20) )

    read(20,*) str2
    read(21,*) str3
    write(*,*) 'data files: ', str2(1:len_trim(str2)),' ', str3(1:len_trim(str3))

    open(10, file = './out/'//str3(1 : len_trim(str3)), form = 'unformatted')
      read(10) dissrate
    close(10)

    open(10, file = './out/'//str2(1 : len_trim(str2)), form = 'unformatted')
      read(10) alpha
      read(10) beta
      read(10) gmma
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

        numpnt = numpnt + 1

        atmp = alpha(ii,jj,kk)
        btmp = beta(ii,jj,kk)
        gtmp = gmma(ii,jj,kk)
       
        ss = sqrt( 2 * ( atmp * atmp + btmp * btmp + gtmp * gtmp ) )
        ss = -12 * sqrt(3.) * atmp * btmp * gtmp / ss**3  ! s^* defined by Lund and Rogers
       
        ll = floor( (ss + 1) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) psatcond(ll) = psatcond(ll) + 1

      end if

    end do
    end do
    end do

  end do

  write(*,*) 'Number of points: ', numpnt

  psatcond = psatcond / numpnt 
  write(*,*) 'check psatcond: ', sum(psatcond)
  psatcond = psatcond / binw

  str2 = 'condpsat-'//str(1:len_trim(str))//str5(1:len_trim(str5))//'.dat'
  open(10, file = str2(1:len_trim(str2)))
    do ii = 1, npnt
      write(10,*) -1. + (ii - .5) * binw, psatcond(ii)
    end do
  close(10)

  deallocate(alpha,beta,gmma)
  deallocate(dissrate)

  write(*,*) 'finished'
end program condpsat      
