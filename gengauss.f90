! Generate Gaussian field with given energy and helicity spectra
program generategaussian
  use mconstant
  use mwavenumber
  implicit none

  integer :: nx, lx, lx1, ly, lz, ii, jj, kk, ll, idum
  real(sp)    :: kxii, kyjj, kzkk, k,k2d
  complex(sp) :: hpx,hpy,hpz,hmx,hmy,hmz

  real(sp), allocatable, dimension(:) :: fpk, fmk, ek, hk, cntep, cntem,cnthp,cnthm
  real(sp), allocatable, dimension(:) :: kx,ky,kz
  real(sp), allocatable, dimension(:,:,:) :: k2
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz,ap,am

  character(80) :: str, str1, str2

  real :: ran1
  integer :: iseed

  write(*,*) 
  write(*,'(''>>>>>> Generate a Gaussian velocity field with read-in energy and helicity spectra <<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 5) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./gengauss.x nx filename1 filename2 filename3 iseed '
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filename1: data file for energy spectrum'
          write(*,*) '        filename2: data file for helicity spectrum'
          write(*,*) '        filename3: output file number'
          write(*,*) '        iseed: minus iseed is the seed of ran1'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if



  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! iseed
  call getarg(5,str)
  read(str, '(I20)') iseed
  ! file name for energy spectrum
  call getarg(2,str)
  str = adjustl(str)
  ! file name for helicity spectrum
  call getarg(3,str1)
  str1 = adjustl(str1)
  ! file number for output
  call getarg(4,str2)
  str2=adjustl(str2)

  idum = -iseed

  lx=nx/2
  lx1=lx+1
  ly=nx
  lz=nx

  allocate(fpk(lx),fmk(lx),ek(lx),hk(lx),cntep(lx),cntem(lx),cnthp(lx),cnthm(lx))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))
  allocate(ap(lx1,ly,lz),am(lx1,ly,lz))

  open(15,file=str(1:len_trim(str))) ! data for energy spectrum
  open(16,file=str1(1:len_trim(str1))) ! data for helicity spectrum
  do ii=1,lx
    read(15,*) jj, ek(ii)  ! data in jj is not needed
    read(16,*) jj, hk(ii)
  end do
  close(15)
  close(16)

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  do kk=1,lz
  do jj=1,ly
  do ii=1,lx1

    k = ran1(idum)
    k2d = ran1(idum)
    k2d = 2._sp*pi*k2d
    am(ii,jj,kk)=gasdev(k,k2d)

    k = ran1(idum)
    k2d = ran1(idum)
    k2d = 2._sp*pi*k2d
    ap(ii,jj,kk)=gasdev(k,k2d)

  end do
  end do
  end do

  call setzero(ap, k2, lx1, ly, lz)
  call setzero(am, k2, lx1, ly, lz)
  call antihermit(ap, lx1, ly, lz)
  call antihermit(am, lx1, ly, lz)

  cntep=0.
  cntem=0.
  cnthp=0.
  cnthm=0.
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1

    if ( ii .eq. 1 ) then ! kx = 0 plane
        kxii=.5*ap(ii,jj,kk)*conjg(ap(ii,jj,kk))
        k2d =.5*am(ii,jj,kk)*conjg(am(ii,jj,kk))
    else
        kxii=ap(ii,jj,kk)*conjg(ap(ii,jj,kk))
        k2d =am(ii,jj,kk)*conjg(am(ii,jj,kk))
    end if

    k=sqrt(k2(ii,jj,kk))
    ll=floor(k+.5)   !!! When 0.5 =< k < 1.5
    if (ll .ge. 1 .and. ll .le. lx) then
        cntep(ll)=cntep(ll)+kxii
        cntem(ll)=cntem(ll)+k2d
        cnthp(ll)=cnthp(ll)+k*kxii
        cnthm(ll)=cnthm(ll)+k*k2d
    end if

  end do
  end do
  end do

  ! To match the spectra, we follow exactly the way the spectra were calculated.
  ! Formula: u(k)= fpk * ap *hp + fmk * am * hm
  ! where fpk and fmk were fixed by the energy and helicity spectra.
  do ii=1,lx
    fpk(ii)=sqrt(ek(ii)*cnthm(ii)+.5*hk(ii)*cntem(ii))/sqrt(cntep(ii)*cnthm(ii)+cnthp(ii)*cntem(ii))
    fmk(ii)=sqrt(max(ek(ii)*cnthp(ii)-.5*hk(ii)*cntep(ii),0.))/sqrt(cntem(ii)*cnthp(ii)+cnthm(ii)*cntep(ii))
  end do

  do kk=1,lz
    kzkk=kz(kk)

    do jj=1,ly
      kyjj=ky(jj)

      do ii=1,lx1
        kxii=kx(ii)
 
        k=sqrt(k2(ii,jj,kk))
        ll=floor(k+.5)
        if (ll .ge. 1 .and. ll .le. lx) then
            am(ii,jj,kk)=am(ii,jj,kk)*fmk(ll)
            ap(ii,jj,kk)=ap(ii,jj,kk)*fpk(ll)
        end if
 
        if (jj .eq. 1 .and. kk .eq. 1) then  ! on the kx axis
            k2d=max(sqrt(kxii*kxii+kyjj*kyjj), mytiny_sp)
            
            hpx=(sqrt(2._sp)/2._sp)*(-kyjj/k2d-eye*kxii*kzkk/k/k2d)
            hpy=(sqrt(2._sp)/2._sp)*(kxii/k2d-eye*kyjj*kzkk/k/k2d)
            hpz=(sqrt(2._sp)/2._sp)*eye*(kyjj*kyjj+kxii*kxii)/k/k2d
        else
            k2d=sqrt(kzkk*kzkk+kyjj*kyjj)
            
            hpx=(sqrt(2._sp)/2._sp)*eye*(kyjj*kyjj+kzkk*kzkk)/k/k2d
            hpy=(sqrt(2._sp)/2._sp)*(-kzkk/k2d-eye*kxii*kyjj/k/k2d)
            hpz=(sqrt(2._sp)/2._sp)*(kyjj/k2d-eye*kzkk*kxii/k/k2d)
        end if
        hmx=conjg(hpx)
        hmy=conjg(hpy)
        hmz=conjg(hpz)
 
        ux(ii,jj,kk)=ap(ii,jj,kk)*hpx+am(ii,jj,kk)*hmx
        uy(ii,jj,kk)=ap(ii,jj,kk)*hpy+am(ii,jj,kk)*hmy
        uz(ii,jj,kk)=ap(ii,jj,kk)*hpz+am(ii,jj,kk)*hmz
      end do
    end do
  end do

  call hermitianize(ux,lx1,ly,lz)
  call hermitianize(uy,lx1,ly,lz)
  call hermitianize(uz,lx1,ly,lz)

  open(10,file='./out/ux'//str2(1:len_trim(str2))//'.dat',form='unformatted')
    write(10) ux
  close(10)
  open(10,file='./out/uy'//str2(1:len_trim(str2))//'.dat',form='unformatted')
    write(10) uy
  close(10)
  open(10,file='./out/uz'//str2(1:len_trim(str2))//'.dat',form='unformatted')
    write(10) uz
  close(10)

  deallocate(ux,uy,uz,am,ap,kx,ky,kz,k2,fpk,fmk,ek,hk,cnthp,cnthm,cntep,cntem)

  write(*,*) 'Done. The output data files are u[xyz]'//str2(1:len_trim(str2))//'.dat.'

contains

  complex(sp) function gasdev(t1,t2)
    real(sp) :: t1, t2

    gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
    ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
  end function gasdev

end program generategaussian      
