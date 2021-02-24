program twodvelincskew
  use mconstant
  use mfftwplan2d
  implicit none
  
  integer :: nx,ny,ii,jj,kk,lx1,lx,ly,ndr, nfile
  real(sp) :: ignore_me

  complex(sp), allocatable, dimension(:,:) :: ux,uy
  real(sp), allocatable, dimension(:,:) :: uxr, uyr, dur, tmp
  real(dp), allocatable, dimension(:)  :: mdu2, mdv2, mdu3, mdv3
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dvelincskewness.x nx filelist'
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


  allocate( ux(lx1,ly), uy(lx1,ly) )
  allocate( uxr(nx,ny), uyr(nx,ny), dur(nx,ny), tmp(nx,ny) )
  allocate( mdu2(ndr), mdv2(ndr), mdu3(ndr), mdv3(ndr) )

  ! file number 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  write(*,*) 'fftwplan'
  call fftwplan2de(nx,ny)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  mdu2 = 0.d0; mdv2 = 0.d0
  mdu3 = 0.d0; mdv3 = 0.d0
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
    
    call rfftwnd_f77_one_complex_to_real(c2r2d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,uy,ignore_me)
    
    uxr(1:nx:2,:) =  real( ux(1:lx,:) )
    uxr(2:nx:2,:) = aimag( ux(1:lx,:) )
    uyr(1:nx:2,:) =  real( uy(1:lx,:) )
    uyr(2:nx:2,:) = aimag( uy(1:lx,:) )
 
    
    jj = 1
    do kk = 1, ndr

      jj = jj * 2

      dur = cshift(uxr, jj, 1) - uxr
      tmp = dur * dur
      mdu2(kk) = mdu2(kk) + sum(tmp)
      tmp = tmp * dur
      mdu3(kk) = mdu3(kk) + sum(tmp)

      dur = cshift(uyr, jj, 2) - uyr
      tmp = dur * dur
      mdu2(kk) = mdu2(kk) + sum(tmp)
      tmp = tmp * dur
      mdu3(kk) = mdu3(kk) + sum(tmp)

      dur = cshift(uyr, jj, 1) - uyr
      tmp = dur * dur
      mdv2(kk) = mdv2(kk) + sum(tmp)
      tmp = tmp * dur
      mdv3(kk) = mdv3(kk) + sum(tmp)

      dur = - cshift(uxr, jj, 2) + uxr
      tmp = dur * dur
      mdv2(kk) = mdv2(kk) + sum(tmp)
      tmp = tmp * dur
      mdv3(kk) = mdv3(kk) + sum(tmp)

    end do

    nfile = nfile + 1

  end do
  close(30)

  ignore_me = 1. / (nx*ny) / nfile / 2

  mdu2 = mdu2 * ignore_me
  mdu3 = mdu3 * ignore_me
  mdv2 = mdv2 * ignore_me
  mdv3 = mdv3 * ignore_me
  
  path='2dvelincskewness-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(''# variables = "dr", "du2","du3","dv2", "dv3"'')') 
    ignore_me = 2*pi/nx
    jj = 1
    do ii=1,ndr
      jj = jj * 2
      write(15,'(10E12.4)') jj*ignore_me, mdu2(ii), mdu3(ii), mdv2(ii), mdv3(ii)
    end do
  close(15)
  
  
  deallocate(ux,uy,uxr,uyr,dur,tmp)
  deallocate(mdu2, mdu3, mdv2, mdv3)
  
  call destroyplan2d

  write(*,*) 'done.'
end program twodvelincskew
