program hdcndbeta

  use mconstant
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll
  integer(8) :: numpnt

  real, parameter :: mdissrate = 0.243 ! see helidiss.dat

  real, allocatable, dimension(:,:,:) :: rbeta
  complex, allocatable, dimension (:,:,:) :: dissrate

  integer, parameter :: npnt = 160, nmean = 40
  real(dp), dimension (npnt) :: cndpdfhd

  real(dp) :: meandiss

  real :: binw, bound

  real :: btmp, ss
  character(80) :: strdiss, strcnd, strrij, strtmp, strtmp1, strtmp2


  nx = 256

  strdiss = 'helidiss8dx'
  strrij = 'rijev8dx'
  strcnd = 'gt0'


  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1

  allocate( rbeta(nx, ny, nz) )
  allocate( dissrate(lx1, ly, lz) )
  write(*,*) 'arrays allocated'


  bound = nmean * mdissrate
  binw = 2. * bound / npnt

  open(20, file = strdiss(1:len_trim(strdiss))//'.list')
  open(21, file = strrij(1:len_trim(strrij))//'.list')
  
    numpnt = 0
    meandiss = 0.d0
    cndpdfhd = 0.d0

    do while ( .not. eof(20) )

      read(20,*) strtmp
      read(21,*) strtmp1
 
      open(10, file = './out/'//strtmp(1 : len_trim(strtmp)), form = 'unformatted')
        read(10) dissrate
      close(10)
 
      open(10, file = './out/'//strtmp1(1 : len_trim(strtmp1)), form = 'unformatted')
        read(10) 
        read(10) rbeta
      close(10)
 
 
      write(*,*) 'after reading data files'
 
      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        
        btmp = rbeta(ii,jj,kk)

        if ( btmp .ge. 0.) then  
        ! conditioned on rbeta > 0
 
          if ( mod(ii, 2) .eq. 1) then
                ss = real( dissrate( (ii+1)/2, jj, kk ) )
          else 
                ss = aimag( dissrate( ii/2, jj,kk ) )
          end if
 
 
          numpnt = numpnt + 1
 
          meandiss = meandiss + ss
         
          ll = floor( (ss + bound) / binw ) + 1
          if (ll .ge. 1 .and. ll .le. npnt) cndpdfhd(ll) = cndpdfhd(ll) + 1
 
        end if
 
      end do
      end do
      end do
 
    end do

  close(20)
  close(21)

  write(*,*) 'Number of points: ', numpnt

  meandiss = meandiss / numpnt

  write(*,*) 'Mean sgs diss: ', meandiss

  cndpdfhd = cndpdfhd / numpnt 
  write(*,*) 'check cndpdfhd: ', sum(cndpdfhd)
  cndpdfhd = cndpdfhd / binw

  strtmp2 = 'cndpdfhd-'//strdiss(1:len_trim(strdiss))//'-'//strcnd(1:len_trim(strcnd))//'.dat'
  open(10, file = strtmp2(1:len_trim(strtmp2)) )
    do ii = 1, npnt
      write(10,*) - bound + (ii - .5) * binw, cndpdfhd(ii)
    end do
  close(10)

  deallocate(rbeta)
  deallocate(dissrate)

  write(*,*) 'finished'

end program hdcndbeta      
