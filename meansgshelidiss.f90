program meansgshelidiss
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ndel,numfilter,ii,ll,jj,kk,nfile

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz,wa,wb
  complex(sp), allocatable, dimension(:,:,:) :: tij, rij, helidiss
  real(sp),    allocatable, dimension(:,:,:) :: k2, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz, meanhelidiss

  real(sp) :: delta_c, ignore_me
  character(80) :: str, str1, filelist

  write(*,*) 
  write(*,'(''>>> Mean sgs helicity dissipation <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgshelidiss.x nx filelist ndel numfilter'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: filelist.list = the list of data files'
          write(*,*) '        ndel: maximum filter scale delta = ndel * dx'
          write(*,*) '        numfilter: number of scales, each halved from a larger one'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if 

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(2,filelist) 
  filelist=adjustl(filelist)

  call getarg(3,str)
  read(str, '(I20)') ndel

  call getarg(4,str)
  read(str, '(I20)') numfilter

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wa(lx1,ly,lz),wb(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1), ky(ly), kz(lz))
  allocate(rij(lx1,ly,lz),tij(lx1,ly,lz))
  allocate(g(lx1,ly,lz), helidiss(lx1,ly,lz))
  allocate(meanhelidiss(numfilter))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  meanhelidiss = 0._dp
  nfile = 0
  open(30, file = filelist(1:len_trim(filelist))//'.list')

  do while ( .not. eof(30) )

    read(30, *) str1
    str1 = adjustl(str1)
    write(*, *) 'reading ', str1(1:len_trim(str1))

    open(10,file='./out/ux'//str1(1:len_trim(str1)),status='unknown',form='unformatted')
      read(10)ux
    close(10)
    open(10,file='./out/uy'//str1(1:len_trim(str1)),status='unknown',form='unformatted')
      read(10)uy
    close(10)
    open(10,file='./out/uz'//str1(1:len_trim(str1)),status='unknown',form='unformatted')
      read(10)uz
    close(10)
    write(*,*) 'after reading data files'
 
    do ll=1,numfilter
 
      delta_c=ndel*2*pi/nx/ 2**(ll-1)
 
      ! Gaussian filter
      g=exp(-k2*delta_c**2/24.)
 
 
      ! wa = wx; wb = wy
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
      wa(ii,jj,kk) = eye * (ky(jj) * uz(ii,jj,kk) - kz(kk) * uy(ii,jj,kk)) * g(ii,jj,kk)
      wb(ii,jj,kk) = eye * (kz(kk) * ux(ii,jj,kk) - kx(ii) * uz(ii,jj,kk)) * g(ii,jj,kk)
      end do
      end do
      end do
 
      
      ! rxx
      do ii = 1, lx1
      rij(ii,:,:) = eye * kx(ii) * wa(ii,:,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      ! txx
      call rsgstauij(ux,ux,tij,g,nx,ny,nz) 
 
      helidiss=cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))
 
      ! tyy
      call rsgstauij(uy,uy,tij,g,nx,ny,nz) 
 
      ! ryy
      do jj = 1, ly
      rij(:,jj,:) = eye * ky(jj) * wb(:,jj,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      helidiss=helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))
 
      ! txy
      call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
 
      ! rxy
      do jj = 1, ly
      do ii = 1, lx1
        rij(ii,jj,:) = .5 * eye * (kx(ii) * wb(ii,jj,:) + ky(jj) * wa(ii,jj,:))
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      helidiss=helidiss + 2.*cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))
 
 
      ! wa = wx; wb = wz
      do jj = 1, ly
      do ii = 1, lx1
      wb(ii,jj,:) = eye * (kx(ii) * uy(ii,jj,:) - ky(jj) * ux(ii,jj,:)) * g(ii,jj,:)
      end do
      end do
 
      ! rzz
      do kk = 1, lz
      rij(:,:,kk) = eye * kz(kk) * wb(:,:,kk)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      ! tzz
      call rsgstauij(uz,uz,tij,g,nx,ny,nz) 
      
      helidiss=helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))
 
      ! rxz
      do kk = 1, lz
      do ii = 1, lx1
      rij(ii,:,kk) = .5 * eye * (kx(ii) * wb(ii,:,kk) + kz(kk) * wa(ii,:,kk))
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      ! txz
      call rsgstauij(ux,uz,tij,g,nx,ny,nz) 
 
      helidiss=helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))
 
      ! wa = wy; wb=wz
      do kk = 1, lz
      do ii = 1, lx1
        wa(ii,:,kk) = eye * (kz(kk) * ux(ii,:,kk) - kx(ii) * uz(ii,:,kk)) * g(ii,:,kk)
      end do
      end do
 
      ! ryz
      do kk = 1, lz
      do jj = 1, ly
      rij(:,jj,kk) = .5 * eye * ( kz(kk) * wa(:,jj,kk) + ky(jj) * wb(:,jj,kk) )
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
 
      ! tyz
      call rsgstauij(uy,uz,tij,g,nx,ny,nz) 
 
      helidiss=helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij))

      helidiss = -2._sp * helidiss
 
      meanhelidiss(ll) = meanhelidiss(ll) &
          +(sum(real(helidiss(1:lx,:,:)))+sum(aimag(helidiss(1:lx,:,:))))/(nx*ny*nz)
 
    end do

    nfile = nfile + 1
  end do
  close(30)

  meanhelidiss = meanhelidiss / nfile

  open(15, file = 'meansgshd-'//filelist(1:len_trim(filelist))//'.dat')
    do ll = 1, numfilter
      delta_c=ndel*2*pi/nx/ 2**(ll-1)
      write(15,*) ndel/2**(ll-1), delta_c, meanhelidiss(ll)
    end do
  close(15)

  call destroyplan3d

  deallocate(ux,uy,uz,wa,wb)
  deallocate(tij,rij,helidiss,k2,kx,ky,kz,g, meanhelidiss)

  write(*,*) 'finished'
end program meansgshelidiss      
