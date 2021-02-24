program c2r
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, ly, lz, lx1, ndel, nfilter

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:,:,:) :: uxr,uyr,uzr

  real(sp) :: delta_c, ignore_me, kcut, kcut2
  character(80) :: str, fnm, str1


  write(*,*) 
  write(*,'(''>>> Complex data to real data <<<'')')
  write(*,*)
  if (iargc() .ne. 4) then
      write(*,*)
      write(*,*) ' >>>>>> Wrong number of arguments <<<<<<' 
      write(*,*) 
      write(*,*) ' Usage: ./c2r.x nx filelist ndel filter' 
      write(*,*) '                nx: resolution'
      write(*,*) '                filelist: filelist.list gives the list of data files'
      write(*,*) '                ndel: filter scale = ndel * dx'
      write(*,*) '                filter: 1 for Gaussian, 0 for cutoff '
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

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel

  ! filter type
  call getarg(4,str)
  read(str, '(I20)') nfilter

  ny=nx; nz=nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1

  delta_c=ndel*(2._sp*pi/real(nx,sp))
  kcut = (nx/2._sp)/ndel
  kcut2 = kcut * kcut

  call fftwplan3de(nx,ny,nz)

  write(*,*) 'Allocate array'
  allocate(G(lx1,ly,lz),ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(uxr(nx,ny,nz),uyr(nx,ny,nz),uzr(nx,ny,nz))

  write(*,*) 'wavenumber'
  call wavenumber(g,lx1,ly,lz)

  write(*,*) 'define filter'

  if ( nfilter .eq. 1 ) g=exp(-g*delta_c**2/24._sp)  ! Gaussian filter

  open( 20, file=fnm(1:len_trim(fnm))//'.list' )

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      write(*,*) 'reading data'
      open(16,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
        read(16) ux
      close(16)
      open(16,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
        read(16) uy
      close(16)
      open(16,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
        read(16) uz
      close(16)
  
      if ( nfilter .eq. 1) then 
          ux=g*ux; uy=g*uy; uz=g*uz ! Gaussian filter
      else
          where( g .ge. kcut2 ) ! cutoff filter 
              ux = 0._sp; uy = 0._sp; uz = 0._sp
          endwhere
      endif
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
  
      uxr(1:nx:2,:,:) =  real(ux(1:lx,:,:))
      uxr(2:nx:2,:,:) = aimag(ux(1:lx,:,:))
      uyr(1:nx:2,:,:) =  real(uy(1:lx,:,:))
      uyr(2:nx:2,:,:) = aimag(uy(1:lx,:,:))
      uzr(1:nx:2,:,:) =  real(uz(1:lx,:,:))
      uzr(2:nx:2,:,:) = aimag(uz(1:lx,:,:))
  
  
      open(12,file='./out/uxr'//str1(1:len_trim(str1)),form='unformatted')
        write(12)uxr
      close(12)
      open(12,file='./out/uyr'//str1(1:len_trim(str1)),form='unformatted')
        write(12)uyr
      close(12)
      open(12,file='./out/uzr'//str1(1:len_trim(str1)),form='unformatted')
        write(12)uzr
      close(12)

    end do

  close(20)

  deallocate(G,ux,uy,uz,uxr,uyr,uzr)

  call destroyplan3d
  write(*,*) 'Finished'

end program c2r
