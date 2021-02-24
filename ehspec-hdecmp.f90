program ehspechdecmp
  use mconstant
  implicit none

  integer :: nx,lx1,lx,ly,lz,nfile,ii,jj,kk,ll
  real(sp) :: kxii, kyjj, kzkk 

  real(sp), allocatable, dimension(:) :: kx, ky, kz
  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz

  real(dp) :: sqrtk2, k2d, ehktmp1, ehktmp2
  real(dp), allocatable, dimension(:) ::  ek, hk
  complex(dp) :: hpx, hpy, hpz, hmx, hmy, hmz, ap, am
  !!!! Experience showed that double precisionis needed to obtain "correct" results.
  !!!! Physically, this is because the energy and helcity in 
  !!!! negative and positive helical modes tend to be 
  !!!! an order of magnitude larger than the total field, see Eyink et al.
                                                                                                             
  character(80) :: str, str1

  write(*,*)
  write(*,*) ' >>>>>> Energy and Helicity spectra (using helical decomposition) <<<<<<'
  write(*,*)

  ii=iargc()
  if (ii .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./ehspec-hdecmp.x nx filelist'
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

  allocate(kx(lx1), ky(ly), kz(lz) )
  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz))
  allocate(ek(lx), hk(lx))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,real(ux),lx1,ly,lz) !!! data in real(ux) not used
  write(*,*) 'after wavenumber'

  open(20, file=str(1:len_trim(str))//'.list')

    hk=0._sp; ek=0._sp
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
      
            sqrtk2=max( sqrt(kxii * kxii + kyjj * kyjj + kzkk * kzkk), mytiny_sp )
         
            if (jj .eq. 1 .and. kk .eq. 1) then  
              ! on the kx axis, another definition for hp and hp is used.
              k2d=max(sqrt(kxii*kxii+kyjj*kyjj), mytiny_sp)
              
              hpx=(sqrt(2.0_dp)/2.0_dp)*(-kyjj/k2d-eye*kxii*kzkk/sqrtk2/k2d)
              hpy=(sqrt(2.0_dp)/2.0_dp)*(kxii/k2d-eye*kyjj*kzkk/sqrtk2/k2d)
              hpz=(sqrt(2.0_dp)/2.0_dp)*eye*(kyjj*kyjj+kxii*kxii)/sqrtk2/k2d
            else
              k2d = max( sqrt(kzkk*kzkk+kyjj*kyjj), mytiny_sp )
              
              hpx=(sqrt(2.0_dp)/2.0_dp)*eye*(kyjj*kyjj+kzkk*kzkk)/sqrtk2/k2d
              hpy=(sqrt(2.0_dp)/2.0_dp)*(-kzkk/k2d-eye*kxii*kyjj/sqrtk2/k2d)
              hpz=(sqrt(2.0_dp)/2.0_dp)*(kyjj/k2d-eye*kzkk*kxii/sqrtk2/k2d)
            endif
            hmx=conjg(hpx)
            hmy=conjg(hpy)
            hmz=conjg(hpz)
         
            ap = ux(ii,jj,kk) * hmx + uy(ii,jj,kk) * hmy + uz(ii,jj,kk) * hmz
            am = ux(ii,jj,kk) * hpx + uy(ii,jj,kk) * hpy + uz(ii,jj,kk) * hpz
         
            ll=floor(sqrtk2+.5) !!! ll-.5 <= k < ll+.5
            if (ll .ge. 1 .and. ll .le. lx) then  
                 if ( ii .eq. 1 ) then 
                     ehktmp1 =.5*ap*conjg(ap) !!! kx = 0 plane
                     ehktmp2 =.5*am*conjg(am) 
                 else 
                     ehktmp1 = ap * conjg(ap) !!! kx /=0 
                     ehktmp2 = am * conjg(am)
                 end if
                 ek(ll) = ek(ll) + ehktmp1 + ehktmp2
                 hk(ll) = hk(ll) + 2. * sqrtk2 * (ehktmp1 - ehktmp2)
            end if
      
          end do
        end do
      end do

      nfile = nfile + 1
    end do

  close(20)

  open(24,file='ehspec-hdecmp-'//str(1:len_trim(str))//'.dat')
  do ii=1,lx
    write(24,*) ii, ek(ii)/nfile, hk(ii)/nfile
  end do
  close(24)

  deallocate(kx,ky,kz,ux,uy,uz,ek,hk)

  write(*,*)
  write(*,*) 'finished'

end program ehspechdecmp      
