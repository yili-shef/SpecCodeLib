program r2c
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, ly, lz, lx1, ll

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension(:,:,:) :: uxr,uyr,uzr

  real(sp) :: ignore_me
  character(80) :: str, fnm, str1


  write(*,*) 
  write(*,'(''>>> Complex data to real data <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 2) then
      write(*,*)
      write(*,*) ' >>>>>> Wrong number of arguments <<<<<<' 
      write(*,*) 
      write(*,*) ' Usage: ./r2c.x nx filelist' 
      write(*,*) '                nx: resolution'
      write(*,*) '                filelist: filelist.list gives the list of data files'
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

  ny=nx; nz=nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1
  ignore_me = 1._sp/(nx*ny*nz)

  call fftwplan3de(nx,ny,nz)

  write(*,*) 'Allocate array'
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(uxr(nx,ny,nz),uyr(nx,ny,nz),uzr(nx,ny,nz))

  open( 20, file=fnm(1:len_trim(fnm))//'.list' )

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      write(*,*) 'reading data'
      open(16,file='./out/uxr'//str1(1:len_trim(str1)),form='unformatted')
        read(16) uxr
      close(16)
      open(16,file='./out/uyr'//str1(1:len_trim(str1)),form='unformatted')
        read(16) uyr
      close(16)
      open(16,file='./out/uzr'//str1(1:len_trim(str1)),form='unformatted')
        read(16) uzr
      close(16)
  
      ux(1:lx,:,:) = cmplx( uxr(1:nx:2,:,:), uxr(2:nx:2,:,:) ) * ignore_me
      uy(1:lx,:,:) = cmplx( uyr(1:nx:2,:,:), uyr(2:nx:2,:,:) ) * ignore_me
      uz(1:lx,:,:) = cmplx( uzr(1:nx:2,:,:), uzr(2:nx:2,:,:) ) * ignore_me
  
      call rfftwnd_f77_one_real_to_complex(r2c3d,ux,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,uy,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,uz,ignore_me)
  
      open(12,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
        write(12)ux
      close(12)
      open(12,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
        write(12)uy
      close(12)
      open(12,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
        write(12)uz
      close(12)

    end do

  close(20)

  deallocate(ux,uy,uz,uxr,uyr,uzr)

  call destroyplan3d
  write(*,*) 'Finished'

end program r2c
