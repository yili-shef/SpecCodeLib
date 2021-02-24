program helispec
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx,lx1,lx,ly,lz, nfile, ii,jj,kk,ll
  real(sp)    :: kxii,kyjj,kzkk, hktmp
  real(sp), allocatable, dimension(:) :: kx,ky,kz
  real(sp), allocatable, dimension(:) :: hk
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  complex(sp) :: wx, wy, wz
  character(80) :: str,str1

  write(*,*)
  write(*,*) ' >>>>>> Helicity spectrum <<<<<<'
  write(*,*)

  ii=iargc()
  if (ii .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./helispec-s(d)p.x nx filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list'
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

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(hk(lx))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,lx1,ly,lz) 
  write(*,*) 'after wavenumber'

  open(20, file=str(1:len_trim(str))//'.list')

    hk=0.
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

            wx=eye*(kyjj*uz(ii,jj,kk)-kzkk*uy(ii,jj,kk))
            wy=eye*(kzkk*ux(ii,jj,kk)-kxii*uz(ii,jj,kk))
            wz=eye*(kxii*uy(ii,jj,kk)-kyjj*ux(ii,jj,kk))
          
            hktmp = 2.*(wx*conjg(ux(ii,jj,kk))+wy*conjg(uy(ii,jj,kk))+wz*conjg(uz(ii,jj,kk)))
            if ( ii .eq. 1 ) hktmp = .5_sp * hktmp

            ll=floor(sqrt(kxii*kxii+kyjj*kyjj+kzkk*kzkk)+.5)
            if (ll .ge. 1 .and. ll .le. lx) then 
                    hk(ll) = hk(ll) + hktmp
            end if
      
          end do
        end do
      end do

      nfile = nfile + 1
    end do

  close(20)
  
  open(24,file='helispec-'//str(1:len_trim(str))//'.dat')
  do ii=1,lx
    write(24,*) ii,hk(ii)/nfile
  end do
  close(24)

  deallocate(kx,ky,kz,ux,uy,uz,hk)

end program helispec
 
