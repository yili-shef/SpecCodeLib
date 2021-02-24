program helipdf
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx1,lx,ly,lz,i,ndel,ii,jj,kk
  real,    allocatable, dimension(:)     :: kx,ky,kz
  real,    allocatable, dimension(:,:,:) :: g,heli, relheli
  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz
  character(80) :: str,str1

  integer, parameter :: npnt=100
  real, dimension(npnt) :: pheli, prelheli

  real :: mh, rmsh, const, ignore_me, delta_c


  write(*,*) 
  write(*,*) ' >>>>>> Calculate the PDF of resolved helicity <<<<<<'
  write(*,*) 

  i=iargc()
  if (i .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./helipdf.x nx nfile ndel'
          write(*,*) '         nx: resolution of data'
          write(*,*) '         nfile: data file number'
          write(*,*) '         ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  call getarg(3,str1)
  read(str1, '(I20)') ndel
  str1=adjustl(str1)
  ! file number string
  call getarg(2,str)
  str = adjustl(str)

  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx
  ny=nx; nz=nx

  const=1./(nx*ny*nz)

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(g(lx1,ly,lz),heli(nx,ny,nz),relheli(nx,ny,nz))
  write(*,*) 'arrays allocated'

  open(15,file='./out/ux'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) ux
  close(15)
  open(15,file='./out/uy'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uy
  close(15)
  open(15,file='./out/uz'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uz
  close(15)
  write(*,*) 'after reading data'


  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (ndel .ne. 0) then
    delta_c=ndel*2.*pi/nx
    g=exp(-g*delta_c**2/24.)
    ux=ux*g
    uy=uy*g 
    uz=uz*g
  end if

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wx(ii,jj,kk)=eye*(ky(jj)*uz(ii,jj,kk)-kz(kk)*uy(ii,jj,kk))
    wy(ii,jj,kk)=eye*(kz(kk)*ux(ii,jj,kk)-kx(ii)*uz(ii,jj,kk))
    wz(ii,jj,kk)=eye*(kx(ii)*uy(ii,jj,kk)-ky(jj)*ux(ii,jj,kk))
  end do
  end do
  end do
 
  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
  
  heli(1:nx:2,:,:)=real(ux(1:lx,:,:))*real(wx(1:lx,:,:)) &
                  +real(uy(1:lx,:,:))*real(wy(1:lx,:,:)) &
                  +real(uz(1:lx,:,:))*real(wz(1:lx,:,:))
  heli(2:nx:2,:,:)=aimag(ux(1:lx,:,:))*aimag(wx(1:lx,:,:))+ &
                   aimag(uy(1:lx,:,:))*aimag(wy(1:lx,:,:))+ &
                   aimag(uz(1:lx,:,:))*aimag(wz(1:lx,:,:))


  ! Relative helicity
  g =  sqrt( real(ux) * real(ux) + real(uy) * real(uy) + real(uz) * real(uz) )
  g = g * sqrt( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz) )
  relheli(1:nx:2,:,:) = heli(1:nx:2,:,:) / g(1:lx,:,:)

  g =  sqrt( aimag(ux) * aimag(ux) + aimag(uy) * aimag(uy) + aimag(uz) * aimag(uz) )
  g = g * sqrt( aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) )
  relheli(2:nx:2,:,:) = heli(2:nx:2,:,:) / g(1:lx,:,:)

  delta_c = 2.0/npnt
  prelheli = 0.
  do kk=1,nz
  do jj=1,ny
  do ii=1,nx
    i=floor((relheli(ii,jj,kk)+1)/delta_c)+1
    if( i .ge. 1 .and. i .le. npnt) prelheli(i)=prelheli(i)+1
  end do
  end do
  end do
  write(*,*) 'check prelheli:', sum(prelheli)/(nx*ny*nz)
  prelheli = prelheli / (nx*ny*nz) / delta_c


  ! true helicity
  mh=sum(heli)*const
  rmsh=sqrt(sum((heli-mh)**2)*const)

  heli=(heli-mh)/rmsh

  ignore_me=maxval(abs(heli))
  write(*,*) 'max heli', ignore_me
  const=2.*ignore_me/npnt

  pheli=0.
  do kk=1,nz
  do jj=1,ny
  do ii=1,nx
    i=floor((heli(ii,jj,kk)+ignore_me)/const)+1
    if( i .ge. 1 .and. i .le. npnt) pheli(i)=pheli(i)+1
  end do
  end do
  end do
  write(*,*) 'check pheli:', sum(pheli)/(nx*ny*nz)
  pheli=pheli/(nx*ny*nz)/const


  open(24,file='helipdf'//str(1:len_trim(str))//'-'//str1(1:len_trim(str1))//'dx.dat')
  write(24,*) '# title=" mean=', mh, '  rms =', rmsh, '"'
  do i=1,npnt
    write(24,*) -ignore_me+(i-.5)*const, pheli(i), -1+(i-0.5)*delta_c, prelheli(i)
  end do
  close(24)


  deallocate(ux,uy,uz,wx,wy,wz,kx,ky,kz,g,heli,relheli)

  call destroyplan3d

  write(*,*) 'done'

end program helipdf
