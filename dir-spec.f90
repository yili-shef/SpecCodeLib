program dirspec
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly,lz, i
  real(sp), allocatable, dimension(:,:,:) :: k2
  real(sp), allocatable, dimension(:)     :: kx,ky,kz
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  character(80) :: str

  i=iargc()
  if (i .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./dir-spec.x nx nfile'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     nfile: data file number'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! file number string
  call getarg(2,str)
  str = adjustl(str)


  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  allocate(kx(lx1),ky(ly),kz(lz),k2(lx1,ly,lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))

  open(15,file='./out/ux'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) ux
  close(15)
  open(15,file='./out/uy'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uy
  close(15)
  open(15,file='./out/uz'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uz
  close(15)


  write(*,*) 'wavenumber'
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  open(24,file='dir-spec'//str(1:len_trim(str))//'.dat')
  ux = ux*conjg(ux) 
  uy = uy*conjg(uy) 
  uz = uz*conjg(uz)
  ux(1,:,:)=0.5*ux(1,:,:)
  uy(1,:,:)=0.5*uy(1,:,:)
  uz(1,:,:)=0.5*uz(1,:,:)
  do i=1,lx
    write(24,*)i,sum(real(ux),mask=(abs(sqrt(k2)-i).lt.0.5)),sum(real(uy),mask=(abs(sqrt(k2)-i).lt.0.5)), &
                 sum(real(uz),mask=(abs(sqrt(k2)-i).lt.0.5))
  end do
  close(24)

  deallocate(kx,ky,kz,k2,ux,uy,uz)

  write(*,*) 'done.'
end program dirspec
