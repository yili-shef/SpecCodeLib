program highpass
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly,lz,kcut
  real(sp),    allocatable, dimension(:,:,:) :: k2
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp) :: kcut2
  character(80) :: str, str1, fln

  write(*,*)
  write(*,*) ' >>>>>> High pass filtering <<<<<<'
  write(*,*)
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./highpass.x nx kcut filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        kcut: threshold wavenumber, lower is removed'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! threshold wavenumber
  call getarg(2,str1)
  read(str1, '(I20)') kcut
  str1 = adjustl(str1)

  ! file number string
  call getarg(3,fln)
  fln = adjustl(fln)

  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx
  kcut2 = kcut * kcut

  allocate(k2(lx1,ly,lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))

  call wavenumber(k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(30, file = fln(1:len_trim(fln))//'.list')

  do while( .not. eof(30) )
    read(30,*) str
    write(*,*) 'reading data', str(1:len_trim(str))

    open(15,file='./out/ux'//str(1:len_trim(str)),form='unformatted')
      read(15) ux
    close(15)
    open(15,file='./out/uy'//str(1:len_trim(str)),form='unformatted')
      read(15) uy
    close(15)
    open(15,file='./out/uz'//str(1:len_trim(str)),form='unformatted')
      read(15) uz
    close(15)
 
    where (k2 .le. kcut2 + mytiny)
            ux = (0.,0); uy = (0.,0.); uz = (0.,0.)
    end where
 
    str = '-kc'//str1(1:len_trim(str1))//'-'//str(1:len_trim(str))
 
    open(15,file='./out/ux'//str(1:len_trim(str)),form='unformatted')
      write(15) ux
    close(15)
    open(15,file='./out/uy'//str(1:len_trim(str)),form='unformatted')
      write(15) uy
    close(15)
    open(15,file='./out/uz'//str(1:len_trim(str)),form='unformatted')
      write(15) uz
    close(15)

  end do
  close(30)


  deallocate(k2,ux,uy,uz)

  write(*,*) 'highpass.x done'

end program highpass      
