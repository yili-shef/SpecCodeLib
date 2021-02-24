program meanek
  use mconstant
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  real(sp), allocatable, dimension(:,:,:) :: q

  character(80) :: str, fnm



  write(*,*)
  write(*,'('' >>> Calculate the mean TKE for a given list of data files <<< '')')
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
  allocate( q(lx1,ly,lz) )

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
 
    q = real(ux*conjg(ux) + uy*conjg(uy) + uz*conjg(uz))
    q(1,:,:)=0.5_sp*q(1,:,:)
    write(*,*) ' data file ', str(1:len_trim(str)), sum(q)

  end do

  close(20)

  deallocate ( ux, uy, uz, q )

  write(*,*) 'finished'

end program meanek      
