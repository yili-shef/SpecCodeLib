program meanhk
  use mconstant
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, wx, wy, wz
  real(sp), allocatable, dimension(:,:,:) :: q
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  character(80) :: str, fnm



  write(*,*)
  write(*,'('' >>> Calculate the mean helicity for a given list of data files <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meanek.x nx filelist'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: the list of data files, *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  
  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny = nx; nz = nx;
  lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz

  allocate( ux(lx1, ly, lz), uy(lx1, ly, lz), uz(lx1, ly, lz) )
  allocate( wx(lx1, ly, lz), wy(lx1, ly, lz), wz(lx1, ly, lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( q(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,q,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(20, file = fnm(1:len_trim(fnm))//'.list')

  do while ( .not. eof(20) )

    read(20,*) str
 
    open(10, file = './out/ux'//str(1:len_trim(str)), form = 'unformatted')
      read(10) ux
    close(10)
    open(10, file = './out/uy'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uy
    close(10)
    open(10, file = './out/uz'//str(1:len_trim(str)), form = 'unformatted')
      read(10) uz
    close(10)

    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk) = eye * ( ky(jj) * uz(ii,jj,kk) - kz(kk) * uy(ii,jj,kk) )
      wy(ii,jj,kk) = eye * ( kz(kk) * ux(ii,jj,kk) - kx(ii) * uz(ii,jj,kk) )
      wz(ii,jj,kk) = eye * ( kx(ii) * uy(ii,jj,kk) - ky(jj) * ux(ii,jj,kk) )
    end do
    end do
    end do
 
    q = 2. * real(conjg(ux)*wx + conjg(uy)*wy + conjg(uz)*wz)
    q(1,:,:)=0.5_sp*q(1,:,:)
    write(*,*) ' data file ', str(1:len_trim(str)), sum(q)

  end do

  close(20)

  deallocate ( ux, uy, uz, wx, wy, wz, q, kx, ky, kz )

  write(*,*) 'finished'

end program meanhk      
