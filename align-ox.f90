program alignox
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex, allocatable, dimension(:,:,:) :: wx,wy,wz,xix,xiy,xiz
  real, allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt=160
  real,    parameter :: binw=1./npnt, binwphi=pi/2/npnt
  real(dp), dimension(npnt) :: pox

  real :: delta_c, ignore_me, tmp, tmp1, tmp2, tmp3
  real(dp) :: alpha
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of Rij <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-ox.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: dns data file list'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! filter parameter
  call getarg(3,str)
  read(str, '(I20)') ndel
  ! file list string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(xix(lx1,ly,lz),xiy(lx1,ly,lz),xiz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pox=0._dp
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xix
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xiy
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)xiz
    close(10)
    write(*,*) 'after reading data files'
 
    xix=xix*g
    xiy=xiy*g 
    xiz=xiz*g
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk)=eye*(ky(jj)*xiz(ii,jj,kk)-kz(kk)*xiy(ii,jj,kk))
      wy(ii,jj,kk)=eye*(kz(kk)*xix(ii,jj,kk)-kx(ii)*xiz(ii,jj,kk))
      wz(ii,jj,kk)=eye*(kx(ii)*xiy(ii,jj,kk)-ky(jj)*xix(ii,jj,kk))

      xix(ii,jj,kk) = eye * ( ky(jj) * wz(ii,jj,kk) - kz(kk) * wy(ii,jj,kk) )
      xiy(ii,jj,kk) = eye * ( kz(kk) * wx(ii,jj,kk) - kx(ii) * wz(ii,jj,kk) )
      xiz(ii,jj,kk) = eye * ( kx(ii) * wy(ii,jj,kk) - ky(jj) * wx(ii,jj,kk) )
    end do
    end do
    end do

 
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,xix,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiz,ignore_me)
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
 
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii+1)/2
 
        delta_c = real( wx(ll, jj, kk) )
        ignore_me = real( wy(ll, jj, kk) )
        tmp = real( wz(ll, jj, kk) )

        tmp1 = real( xix(ll,jj,kk) )
        tmp2 = real( xiy(ll,jj,kk) )
        tmp3 = real( xiz(ll,jj,kk) )
  
      else
        ll = ii/2
 
        delta_c = aimag( wx(ll, jj, kk) )
        ignore_me = aimag( wy(ll, jj, kk) )
        tmp = aimag( wz(ll, jj, kk) )
  
        tmp1 = aimag( xix(ll,jj,kk) )
        tmp2 = aimag( xiy(ll,jj,kk) )
        tmp3 = aimag( xiz(ll,jj,kk) )
      end if
 
 
      alpha = delta_c*tmp1+ignore_me*tmp2+tmp*tmp3
 
      delta_c = sqrt( delta_c * delta_c + ignore_me * ignore_me + tmp * tmp )
      ignore_me = sqrt( tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3 )
 
      alpha = abs(alpha) / delta_c / ignore_me
      
      ll=floor((alpha)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt )  then 
 
        pox(ll)=pox(ll)+1
      end if

    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  alpha = 1._dp/(nx*ny*nz)/nfile

  write(*,*) 'check pox:', sum(pox)*alpha

  pox=pox*alpha/binw

  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='poxalign'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, pox(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,g,kx,ky,kz)
  deallocate(xix,xiy,xiz)

  call destroyplan3d

  write(*,*) 'done'

end program alignox
