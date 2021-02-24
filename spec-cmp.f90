#include "constant_dp.f90"
#include "wavenumber.f90"

program spec
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly,lz, nfile, ii,jj,kk,ll
  real(sp) :: kxii,kyjj,kzkk, ektmpx, ektmpy, ektmpz
  real(sp), allocatable, dimension(:) :: kx, ky, kz
  real(sp), allocatable, dimension(:) :: ekx, eky, ekz
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
          write(*,*) ' Usage: ./spec-cmp.x nx filelist'
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
  allocate(ekx(lx), eky(lx), ekz(lx), k2(lx1,ly,lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  write(*,*) 'arrays allocated'
  
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(20, file=str(1:len_trim(str))//'.list')

    ekx=0._sp; eky = 0._sp; ekz = 0._sp
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
      
            ektmpx = real(ux(ii,jj,kk)*conjg(ux(ii,jj,kk)))
            ektmpy = real(uy(ii,jj,kk)*conjg(uy(ii,jj,kk)))
            ektmpz = real(uz(ii,jj,kk)*conjg(uz(ii,jj,kk)))

            if ( ii .eq. 1 ) then
              ektmpx = .5 * ektmpx
              ektmpy = .5 * ektmpy
              ektmpz = .5 * ektmpz
            end if

            ll=floor(sqrt(kxii*kxii+kyjj*kyjj+kzkk*kzkk)+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
                    ekx(ll)=ekx(ll)+ektmpx
                    eky(ll)=eky(ll)+ektmpy
                    ekz(ll)=ekz(ll)+ektmpz
            end if
      
          end do
        end do
      end do

      nfile = nfile + 1
    end do
  
  open(24,file='spec-cmp-'//str(1:len_trim(str))//'.dat')
  do ii=1,lx
    write(24,'(I6, 10E12.4)') ii, ekx(ii)/nfile, eky(ii)/nfile, ekz(ii)/nfile
  end do
  close(24)

  deallocate(kx,ky,kz,ux,uy,uz,ekx,eky,ekz,k2)

  write(*,*) 'spec-cmp.x finished'
end program spec
 
