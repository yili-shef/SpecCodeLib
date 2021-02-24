program velc2r
  use mconstant
  use mfftwplan3d
  implicit none

  
  integer :: nx, ny, nz, lx, ly, lz, lx1, ll, ndel
  complex(sp), allocatable, dimension(:,:,:)  :: ux
  real(sp), allocatable, dimension(:,:,:)  :: rux, g
  real(sp), allocatable, dimension(:) :: kx,ky,kz

  character(80) :: str, fnm, datafn
  real(sp) :: ignore_me, delta_c


  write(*,*) 
  write(*,'(''>>>>>> Filter and convert complex data files to REAL BINARY files <<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: velc2r.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: .list file containing the list of data files'
          write(*,*) '        ndel: filter scale delta = ndel * dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  endif

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter parameter
  call getarg(3,str)
  read(str, '(I20)') ndel
  str = adjustl(str)

  ! list file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny = nx; nz = nx
  lx=nx/2; ly=nx; lz=nx; lx1=lx+1
  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( ux(lx1,ly,lz) )
  allocate( rux(nx,ny,nz) )
  allocate( g(lx1,ly,lz) )
  write(*,*) 'after allocate'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  do while ( .not. eof(20) )

    read(20,*) datafn
    datafn = adjustl(datafn)
    write(*,*) datafn( 1 : len_trim(datafn) )

    open(10,file='./out/'//datafn(1 : len_trim(datafn)),form='unformatted')
      read(10) ux
    close(10)
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  
    rux(1:nx:2,:,:)= real(ux(1:lx,:,:))
    rux(2:nx:2,:,:)=aimag(ux(1:lx,:,:))

    open(10, file = './out/r'//datafn( 1:len_trim(datafn) ), form = 'binary')
      write(10) rux
    close(10)

  end do
  close(20)
   
  call destroyplan3d

  deallocate(kx, ky, kz, ux, rux, g)

  write(*,*) 'finished'

end program velc2r      
