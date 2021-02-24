program twodspec
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly, nfile, ii,jj,ll
  real(sp) :: kxii,kyjj, ektmp
  real(sp), allocatable, dimension(:) :: kx, ky
  real(sp), allocatable, dimension(:) :: ek
  complex(sp), allocatable, dimension(:,:) :: ux,uy
  real(sp),    allocatable, dimension(:,:) :: k2
  character(80) :: str, str1


  write(*,*)
  write(*,*) ' >>>>>> Energy spectrum <<<<<<'
  write(*,*)
  ii=iargc()
  if (ii .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./2dspec_s(d)p.x nx filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: data file list *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! file list string
  call getarg(2,str)
  str = adjustl(str)

  lx=nx/2;   lx1=lx+1;   ly=nx

  allocate(kx(lx1),ky(ly))
  allocate(ek(lx), k2(lx1,ly))
  allocate(ux(lx1,ly),uy(lx1,ly))
  write(*,*) 'arrays allocated'
  
  call wavenumber(kx,ky,k2,lx1,ly)
  write(*,*) 'after wavenumber'

  open(20, file=str(1:len_trim(str))//'.list')

    ek=0._sp
    nfile = 0
    do while ( .not. eof(20))
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) ux
      close(15)
      open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uy
      close(15)
  
      do jj=1,ly
    
        kyjj = ky(jj)
        kyjj = kyjj * kyjj
        do ii=1,lx1
    
          kxii = kx(ii)
          kxii = kxii * kxii
    
          ektmp = real( ux(ii,jj) * conjg(ux(ii,jj)) + uy(ii,jj) * conjg(uy(ii,jj)) ) 

          if ( ii .eq. 1 ) ektmp = .5 * ektmp

          ll=floor(sqrt(kxii+kyjj)+.5)
          if (ll .ge. 1 .and. ll .le. lx) ek(ll)=ek(ll)+ektmp
    
        end do
      end do

      nfile = nfile + 1
    end do
  
  open(24,file='2dspec-'//str(1:len_trim(str))//'.dat')
  do ii=1,lx
    write(24,*) ii,ek(ii)/nfile
  end do
  close(24)

  deallocate(kx,ky,ux,uy,ek,k2)

  write(*,*) 'finished'
end program twodspec
