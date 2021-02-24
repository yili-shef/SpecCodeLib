program spec
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly,lz, nfile, ii,jj,kk,ll
  real(sp) :: kxii,kyjj,kzkk, ektmp
  real(sp), allocatable, dimension(:) :: kx, ky, kz
  real(sp), allocatable, dimension(:) :: ek
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp),    allocatable, dimension(:,:,:) :: k2
  character(80) :: str, str1


  write(*,*)
  write(*,*) ' >>>>>> Energy spectrum <<<<<<'
  write(*,*)
  ii=iargc()
  if (ii .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./spec_s(d)p.x nx filelist'
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

  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(ek(lx), k2(lx1,ly,lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  write(*,*) 'arrays allocated'
  
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
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
      open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) uz
      close(15)
  
      do kk=1,lz
      
        kzkk=kz(kk)
        do jj=1,ly
      
          kyjj=ky(jj)
          do ii=1,lx1
      
            kxii=kx(ii)
      
            ektmp = real(ux(ii,jj,kk)*conjg(ux(ii,jj,kk)) &
                    + uy(ii,jj,kk)*conjg(uy(ii,jj,kk)) + &
                       uz(ii,jj,kk)*conjg(uz(ii,jj,kk)))

            if ( ii .eq. 1 ) ektmp = .5 * ektmp

            ll=floor(sqrt(kxii*kxii+kyjj*kyjj+kzkk*kzkk)+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    ek(ll)=ek(ll)+ektmp
            end if
      
          end do
        end do
      end do

      nfile = nfile + 1
    end do
  
  open(24,file='spec-'//str(1:len_trim(str))//'.dat')
  do ii=1,lx
    write(24,*) ii,ek(ii)/nfile
  end do
  close(24)

  deallocate(kx,ky,kz,ux,uy,uz,ek, k2)

  write(*,*) 'finished'
end program spec
 
