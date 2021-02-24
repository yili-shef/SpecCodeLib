program vorticitypdf
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none
  
  integer, parameter :: npnt=160
  real(sp), parameter ::  bnd = 16., binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfwx,pdfwy,pdfwz, pdfw
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,nfile,ndel
  real(sp) :: ignore_me, delta_c
  real(dp) :: rmswx, rmswy, rmswz, rmswz0, rmsw

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
          write(*,*) ' Usage: ./pdfvort.x nx filelist ndel'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: delta_c = ndel * dx'
          write(*,*) ' Output: PDFs of vorticity components.'
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
  lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
  delta_c=ndel*2*pi/nx

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(g (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)


  open(30, file = flnm(1:len_trim(flnm))//'.list')

  pdfwx = 0.d0; pdfwy = 0.d0; pdfwz = 0.d0
  rmswx = 0.d0; rmswy = 0.d0; rmswz = 0.d0
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
    
    if ( nfile .eq. 0 ) then
      rmswz0 = sum( real( wz(1:lx,:,:) * conjg(wz(1:lx,:,:)) ) )
      rmswz0 = sqrt( rmswz0 / (nx * ny * nz) )
      write(*,*) 'Estimate rms of |wz|:', rmswz0 
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii, 2) .eq. 1 ) then
        ll = (ii + 1)/2
        ignore_me = real( wx(ll,jj,kk) )
      else
        ignore_me = aimag( wx(ii/2,jj,kk) )
      endif
      rmswx = rmswx + ignore_me * ignore_me
      ignore_me = ignore_me / rmswz0
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfwx(ll)=pdfwx(ll)+1
 

      if ( mod(ii, 2) .eq. 1 ) then
        ll = (ii + 1)/2
        ignore_me = real( wy(ll,jj,kk) )
      else
        ignore_me = aimag( wy(ii/2,jj,kk) )
      endif
      rmswy = rmswy + ignore_me * ignore_me
      ignore_me = ignore_me / rmswz0
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfwy(ll)=pdfwy(ll)+1
 
      if ( mod(ii, 2) .eq. 1 ) then
        ll = (ii + 1)/2
        ignore_me = real( wz(ll,jj,kk) )
      else
        ignore_me = aimag( wz(ii/2,jj,kk) )
      endif
      rmswz = rmswz + ignore_me * ignore_me
      ignore_me = ignore_me / rmswz0
      ll=1+floor((ignore_me + bnd) / binw)
      if (ll .ge. 1 .and. ll .le. npnt) pdfwz(ll)=pdfwz(ll)+1
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile

  rmsw = rmswx + rmswy + rmswz
  rmsw = sqrt( rmsw * ignore_me / 3. )

  rmswx = sqrt( rmswx * ignore_me )
  rmswy = sqrt( rmswy * ignore_me )
  rmswz = sqrt( rmswz * ignore_me )

  write(*,*) ' rms wx:', rmswx
  write(*,*) ' rms wy:', rmswy
  write(*,*) ' rms wz:', rmswz
  write(*,*) ' rms w :', rmsw

  pdfwx=pdfwx*ignore_me
  pdfwy=pdfwy*ignore_me
  pdfwz=pdfwz*ignore_me

  write(*,*) 'check pdfwx:', sum(pdfwx)
  write(*,*) 'check pdfwy:', sum(pdfwy)
  write(*,*) 'check pdfwz:', sum(pdfwz)

  pdfwx=pdfwx/binw
  pdfwy=pdfwy/binw
  pdfwz=pdfwz/binw

  pdfw=(pdfwx + pdfwy + pdfwz)/3.
  
  path='pdfvort-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)

  open(15,file=path)
    write(15,'(''# wx pdfwx wy pdfwy wz pdfwz wi pdfwi'')')
    do lz=1,npnt
      ignore_me = -bnd + (lz - .5) * binw
      write(15,'(20E12.4)') ignore_me * rmswz0 / rmswx, pdfwx(lz) * rmswx / rmswz0, &
                            ignore_me * rmswz0 / rmswy, pdfwy(lz) * rmswy / rmswz0, &
                            ignore_me * rmswz0 / rmswz, pdfwz(lz) * rmswz / rmswz0, &
                            ignore_me * rmswz0 / rmsw, pdfw(lz) * rmsw / rmswz0
    end do

  close(15)
  
  
  deallocate(wx,wy,wz)
  deallocate(kx,ky,kz,g)
  
  call destroyplan3d

  write(*,*) 'done.'
end program vorticitypdf
