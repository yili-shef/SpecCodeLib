program pdfoo
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=160
  real(sp), parameter ::  bnd = 6., binw = bnd / npnt
  real(dp), dimension(npnt) :: pdfo
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile,ndel
  real(sp) :: ignore_me, rmso0, delta_c
  real(dp) :: rmso, tmp

  complex(sp) :: a12, a13, a21, a23, a31, a32 

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfoo.x nx filelist ndel'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*) ' Output: PDFs of vorticity magnitude |w|.'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  delta_c=ndel*2*pi/nx

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate( g(lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfo = 0._dp
  rmso = 0._dp
  nfile = 0

  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) wx
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) wy
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) wz
    close(15)
    write(*,*) 'finishing reading data'

    wx = wx * g
    wy = wy * g
    wz = wz * g

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      a12 = eye * ky(jj) * wx(ii,jj,kk)
      a13 = eye * kz(kk) * wx(ii,jj,kk)
      a21 = eye * kx(ii) * wy(ii,jj,kk)
      a23 = eye * kz(kk) * wy(ii,jj,kk)
      a31 = eye * kx(ii) * wz(ii,jj,kk) 
      a32 = eye * ky(jj) * wz(ii,jj,kk)

      wx(ii,jj,kk) = a32 - a23
      wy(ii,jj,kk) = a13 - a31
      wz(ii,jj,kk) = a21 - a12
    end do
    end do
    end do
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
    
    ! squared magnitude of vorticity.
    wz = cmplx( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz), &
                aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) )

    if ( nfile .eq. 0 ) then
      rmso0 = sum( real( wz(1:lx,:,:) ) ) + sum( aimag( wz(1:lx,:,:) ) ) 
      rmso0 = sqrt( rmso0 / (nx * ny * nz) )
      write(*,*) 'Estimate of |w|: ', rmso0
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
        ll = ( ii + 1 ) / 2
        ignore_me = sqrt( real( wz(ll,jj,kk) ) )
      else
        ll = ii / 2
        ignore_me = sqrt( aimag( wz(ll,jj,kk) ) )
      endif

      rmso = rmso + real(ignore_me * ignore_me, dp)
      ignore_me = ignore_me / rmso0

      ll = 1 + floor((ignore_me / binw))
      if (ll .ge. 1 .and. ll .le. npnt) pdfo(ll)=pdfo(ll)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  write(*,*) 'nfile = ', nfile

  tmp = 1._dp / (nx*ny*nz*nfile)

  rmso = sqrt( rmso * tmp )
  write(*,*) 'Real rms of |w|:', rmso

  pdfo = pdfo * tmp

  write(*,*) 'check pdfo:', sum(pdfo)

  pdfo=pdfo/binw

  path='pdfo-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)

  open(15,file=path)
    do lz=1,npnt
      ignore_me = (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmso0 / rmso, pdfo(lz) * rmso / rmso0
    end do
  close(15)
  
  
  deallocate(wx,wy,wz)
  deallocate(kx,ky,kz,g)
  
  call destroyplan3d

  write(*,*) 'done.'
end program pdfoo
