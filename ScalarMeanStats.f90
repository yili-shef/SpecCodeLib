! Will reading the energy spectrum and calculate the mean energy
! dissipation rate with given viscosity
program scalarmeanstats
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, ly, lz, lx1, ii, jj, kk, ll
  real(dp) :: meangphix, rmsgphix, skwgphix, fltgphix
  real(dp) :: meangphiy, rmsgphiy, skwgphiy, fltgphiy
  real(dp) :: meangphiz, rmsgphiz, skwgphiz, fltgphiz
  real(dp) :: meanphi,  rmsphi,  skwphi,  fltphi
  real(dp) :: phiux, phiuy, phiuz, rmsux, rmsuy, rmsuz

  real(sp) :: ignore_me, const

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, phi
  complex(sp), allocatable, dimension(:,:,:) :: gphix, gphiy, gphiz
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, fnm


  write(*,*) 
  write(*,'(''       >>>>>> Mean parameters of scalar field <<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 2) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./ScalarMeanStats.x nx filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: list of data files, *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! file list string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny = nx; nz = nx;
  lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz

  const = 1./(nx*ny*nz)

  allocate( phi(lx1,ly,lz), gphix(lx1,ly,lz), gphiy(lx1,ly,lz), gphiz(lx1,ly,lz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open( 14, file = 'scalarmeans-'//fnm(1:len_trim(fnm))//'.dat' )
  write(14, '(''# nfile mean rms skewness flatness'')')
  open( 15, file = 'scalargradmeansx-'//fnm(1:len_trim(fnm))//'.dat' )
  write(15,'(''# nfile meangphix rmsgphix skewnessgphix flatnessgphix'')')
  open( 16, file = 'scalargradmeansy-'//fnm(1:len_trim(fnm))//'.dat' )
  write(16,'(''# nfile meangphiy rmsgphiy skewnessgphiy flatnessgphiy'')')
  open( 17, file = 'scalargradmeansz-'//fnm(1:len_trim(fnm))//'.dat' )
  write(17,'(''# nfile meangphiz rmsgphiz skewnessgphiz flatnessgphiz'')')
  open( 18, file = 'meanscalarflux-'//fnm(1:len_trim(fnm))//'.dat' )
  write(18, '(''# nfile phiux phiuy phiuz rmsux rmsuy rmsuz'')')

  open(20, file = fnm(1:len_trim(fnm))//'.list')

  do while ( .not. eof(20) )

    read(20,*) str
    write(*,*) 'data file: ', str(1:len_trim(str))

    open(10, file = './out/ux'//str(1:len_trim(str)), form = 'unformatted')
      read(10) ux
    close(10)
    open(10, file = './out/uy'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uy
    close(10)
    open(10, file = './out/uz'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uz
    close(10)
    open(10, file = './out/phi'//str(1:len_trim(str)), form = 'unformatted')
      read(10) phi
    close(10)
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
        gphix(ii,jj,kk) = eye * kx(ii) * phi(ii,jj,kk)
        gphiy(ii,jj,kk) = eye * ky(jj) * phi(ii,jj,kk)
        gphiz(ii,jj,kk) = eye * kz(kk) * phi(ii,jj,kk)
    end do
    end do
    end do

    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,phi,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphix,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphiy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gphiz,ignore_me)

    meanphi = sum( real(phi(1:lx,:,:)) ) + sum( aimag(phi(1:lx,:,:)) )
    meanphi = meanphi * const
    rmsphi =  sum((real(phi(1:lx,:,:)) - meanphi)**2) + sum((aimag(phi(1:lx,:,:)) - meanphi)**2)
    rmsphi = sqrt(rmsphi * const)
    skwphi =  sum((real(phi(1:lx,:,:)) - meanphi)**3) + sum((aimag(phi(1:lx,:,:)) - meanphi)**3)
    skwphi = skwphi * const
    skwphi = skwphi / rmsphi**3
    fltphi =  sum((real(phi(1:lx,:,:)) - meanphi)**4) + sum((aimag(phi(1:lx,:,:)) - meanphi)**4)
    fltphi = fltphi * const
    fltphi = fltphi / rmsphi**4

    rmsux = sum( real(ux(1:lx,:,:))**2 ) + sum( aimag(ux(1:lx,:,:))**2 )
    rmsux = sqrt( rmsux * const )
    phiux = sum(  real(ux(1:lx,:,:)) *  real(phi(1:lx,:,:)) ) + &
            sum( aimag(ux(1:lx,:,:)) * aimag(phi(1:lx,:,:)) )
    phiux = phiux * const / (rmsphi * rmsux) 

    rmsuy = sum( real(uy(1:lx,:,:))**2 ) + sum( aimag(uy(1:lx,:,:))**2 )
    rmsuy = sqrt( rmsuy * const )
    phiuy = sum(  real(uy(1:lx,:,:)) *  real(phi(1:lx,:,:)) ) + &
            sum( aimag(uy(1:lx,:,:)) * aimag(phi(1:lx,:,:)) )
    phiuy = phiuy * const / (rmsphi * rmsuy) 

    rmsuz = sum( real(uz(1:lx,:,:))**2 ) + sum( aimag(uz(1:lx,:,:))**2 )
    rmsuz = sqrt( rmsuz * const )
    phiuz = sum(  real(uz(1:lx,:,:)) *  real(phi(1:lx,:,:)) ) + &
            sum( aimag(uz(1:lx,:,:)) * aimag(phi(1:lx,:,:)) )
    phiuz = phiuz * const / (rmsphi * rmsuz) 

    meangphix = sum( real(gphix(1:lx,:,:)) ) + sum( aimag(gphix(1:lx,:,:)) )
    meangphix = meangphix * const
    rmsgphix =  sum((real(gphix(1:lx,:,:)) - meangphix)**2) + sum((aimag(gphix(1:lx,:,:)) - meangphix)**2)
    rmsgphix = sqrt(rmsgphix * const)
    skwgphix =  sum((real(gphix(1:lx,:,:)) - meangphix)**3) + sum((aimag(gphix(1:lx,:,:)) - meangphix)**3)
    skwgphix = skwgphix * const
    skwgphix = skwgphix / rmsgphix**3
    fltgphix =  sum((real(gphix(1:lx,:,:)) - meangphix)**4) + sum((aimag(gphix(1:lx,:,:)) - meangphix)**4)
    fltgphix = fltgphix * const
    fltgphix = fltgphix / rmsgphix**4

    meangphiy = sum( real(gphiy(1:lx,:,:)) ) + sum( aimag(gphiy(1:lx,:,:)) )
    meangphiy = meangphiy * const
    rmsgphiy =  sum((real(gphiy(1:lx,:,:)) - meangphiy)**2) + sum((aimag(gphiy(1:lx,:,:)) - meangphiy)**2)
    rmsgphiy = sqrt(rmsgphiy * const)
    skwgphiy =  sum((real(gphiy(1:lx,:,:)) - meangphiy)**3) + sum((aimag(gphiy(1:lx,:,:)) - meangphiy)**3)
    skwgphiy = skwgphiy * const
    skwgphiy = skwgphiy / rmsgphiy**3
    fltgphiy =  sum((real(gphiy(1:lx,:,:)) - meangphiy)**4) + sum((aimag(gphiy(1:lx,:,:)) - meangphiy)**4)
    fltgphiy = fltgphiy * const
    fltgphiy = fltgphiy / rmsgphiy**4

    meangphiz = sum( real(gphiz(1:lx,:,:)) ) + sum( aimag(gphiz(1:lx,:,:)) )
    meangphiz = meangphiz * const
    rmsgphiz =  sum((real(gphiz(1:lx,:,:)) - meangphiz)**2) + sum((aimag(gphiz(1:lx,:,:)) - meangphiz)**2)
    rmsgphiz = sqrt(rmsgphiz * const)
    skwgphiz =  sum((real(gphiz(1:lx,:,:)) - meangphiz)**3) + sum((aimag(gphiz(1:lx,:,:)) - meangphiz)**3)
    skwgphiz = skwgphiz * const
    skwgphiz = skwgphiz / rmsgphiz**3
    fltgphiz =  sum((real(gphiz(1:lx,:,:)) - meangphiz)**4) + sum((aimag(gphiz(1:lx,:,:)) - meangphiz)**4)
    fltgphiz = fltgphiz * const
    fltgphiz = fltgphiz / rmsgphiz**4


    str = str(1: len_trim(str) - 4) ! removing .dat part
    read(str,'(I5)') ll
    write(14, '(I5, 20E12.3)') ll, meanphi, rmsphi, skwphi, fltphi
    write(15, '(I5, 20E12.3)') ll, meangphix, rmsgphix, skwgphix, fltgphix
    write(16, '(I5, 20E12.3)') ll, meangphiy, rmsgphiy, skwgphiy, fltgphiy
    write(17, '(I5, 20E12.3)') ll, meangphiz, rmsgphiz, skwgphiz, fltgphiz
    write(18, '(I5, 20E12.3)') ll, phiux, phiuy, phiuz, rmsux, rmsuy, rmsuz

  end do


  close(14)
  close(15)
  close(16)
  close(17)
  close(18)


  deallocate(ux,uy,uz,phi,gphix,gphiy,gphiz,kx,ky,kz)

  call destroyplan3d
  write(*,*) 'done'

end program scalarmeanstats 
