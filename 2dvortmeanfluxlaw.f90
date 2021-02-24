program twodvortmeanfluxlaw
  use mconstant
  use mwavenumber
  use mfftwplan2d
  implicit none
  
  integer :: nx,ny,ii,jj,kk,lx1,lx,ly,ndr, nfile
  real(sp) :: ignore_me

  complex(sp), allocatable, dimension(:,:) :: ux, uy, om
  real(sp), allocatable, dimension(:,:) :: uxr, uyr, omr, dur, domr
  real(dp), allocatable, dimension(:)  :: meanflux, kx, ky
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dvortmeanfluxlaw.x nx filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ny=nx
  lx=nx/2;lx1=nx/2+1;ly=ny

  ndr = 0
  ii = nx/2
  do while (ii .gt. 1)
    ii = ii / 2
    ndr = ndr + 1
  end do
  ! ndr is log2(nx) - 1
  write(*,*) 'ndr = ', ndr


  allocate( ux(lx1,ly), uy(lx1,ly), om(lx1,ly) )
  allocate( uxr(nx,ny), uyr(nx,ny), dur(nx,ny), omr(nx,ny), domr(nx,ny) )
  allocate( meanflux(ndr), kx(lx1), ky(ly) )

  ! file number 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  write(*,*) 'fftwplan'
  call fftwplan2de(nx,ny)
  call wavenumber(kx,ky,lx1,ly)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  meanflux = 0.d0
  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) ux
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uy
    close(15)

    do jj = 1, ly
    do ii = 1, lx1
      om(ii,jj) = eye * ( kx(ii) * uy(ii,jj) - ky(jj) * ux(ii,jj) )
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r2d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,om,ignore_me)
    
    uxr(1:nx:2,:) =  real( ux(1:lx,:) )
    uxr(2:nx:2,:) = aimag( ux(1:lx,:) )
    uyr(1:nx:2,:) =  real( uy(1:lx,:) )
    uyr(2:nx:2,:) = aimag( uy(1:lx,:) )

    omr(1:nx:2,:) =  real( om(1:lx,:) )
    omr(2:nx:2,:) = aimag( om(1:lx,:) )
 
    jj = 1
    do kk = 1, ndr

      jj = jj * 2

      dur  = cshift(uxr, jj, 1) - uxr
      domr = cshift(omr, jj, 1) - omr
      meanflux(kk) = meanflux(kk) + sum( abs(dur) * domr**2 )

      dur  = cshift(uyr, jj, 2) - uyr
      domr = cshift(omr, jj, 2) - omr
      meanflux(kk) = meanflux(kk) + sum( abs(dur) * domr**2 )

    end do

    nfile = nfile + 1

  end do
  close(30)

  ignore_me = 1. / (nx*ny) / nfile / 2
  meanflux = meanflux * ignore_me
  
  path='2dvortmeanfluxlaw-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(''# variables = "dr", "meanflux"'')') 
    ignore_me = 2*pi/nx
    jj = 1
    do ii=1,ndr
      jj = jj * 2
      write(15,'(10E12.4)') jj*ignore_me, meanflux(ii)
    end do
  close(15)
  
  
  deallocate(ux,uy,om,uxr,uyr,omr,dur,domr)
  deallocate(meanflux, kx, ky)
  
  call destroyplan2d

  write(*,*) 'done.'
end program twodvortmeanfluxlaw
