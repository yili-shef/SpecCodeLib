program meansgsenerdiss
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ndel,numfilter,ii,jj,kk,ll, nfile, ifilter

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  complex(sp), allocatable, dimension(:,:,:) :: tij, sij, enerdiss
  real(sp),    allocatable, dimension(:,:,:) :: k2, g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz 
  real(dp),    allocatable, dimension(:)     :: meanenerdiss

  real(sp) :: delta_c, ignore_me, kc2
  character(80) :: str, str1, filelist, strfltr

  write(*,*) 
  write(*,'(''>>> Mean sgs energy dissipation <<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgsenerdiss.x nx filelist ndel numfilter ifilter'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        ndel: maximum filter scale delta = ndel * dx'
          write(*,*) '        numfilter: number of scales, each halved from a larger one'
          write(*,*) '        ifilter: = 1 for Gaussian; = 0 for cutoff'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if 

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(2,filelist) 
  filelist = adjustl(filelist)

  call getarg(3,str)
  read(str, '(I20)') ndel

  call getarg(4,str)
  read(str, '(I20)') numfilter

  call getarg(5,strfltr)
  read(strfltr, '(I20)') ifilter
  strfltr = adjustl(strfltr)

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1),ky(ly),kz(lz))
  allocate(sij(lx1,ly,lz),tij(lx1,ly,lz))
  allocate(g(lx1,ly,lz), enerdiss(lx1,ly,lz))
  allocate(meanenerdiss(numfilter))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  if ( ifilter .eq. 1 ) then
      strfltr = 'gau'
  else if ( ifilter .eq. 0 ) then
      strfltr = 'cut'
  else
      stop 'Wrong filter type. Stopping '
  end if


  meanenerdiss = 0._dp
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
 
    do ll = 1, numfilter
 
      delta_c = ndel*2*pi/nx / 2**(ll-1)
      kc2 = (pi/delta_c)**2

      if ( ifilter .eq. 1 ) then
          ! Gaussian filter
          g=exp(-k2*delta_c**2/24.)
      else if ( ifilter .eq. 0 ) then
          ! Cutoff 
          where ( k2 .ge. kc2 )
              g = 0
          elsewhere
              g = 1
          endwhere
      end if
 
      ! sxx
      do ii = 1, lx1
      sij(ii,:,:) = eye * kx(ii) * ux(ii,:,:) * g(ii,:,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      ! txx
      call rsgstauij(ux,ux,tij,g,nx,ny,nz) 
 
      enerdiss = cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      ! tyy
      call rsgstauij(uy,uy,tij,g,nx,ny,nz) 
 
      ! syy
      do jj = 1, ly
      sij(:,jj,:) = eye * ky(jj) * uy(:,jj,:) * g(:,jj,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      enerdiss = enerdiss + cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      ! txy
      call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
 
      ! sxy
      do jj = 1, ly
      do ii = 1, lx1
      sij(ii,jj,:) = .5 * eye * (kx(ii) * uy(ii,jj,:) + ky(jj) * ux(ii,jj,:)) * g(ii,jj,:)
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      enerdiss = enerdiss + 2.*cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      ! szz
      do kk = 1, lz
      sij(:,:,kk) = eye * kz(kk) * uz(:,:,kk) * g(:,:,kk)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      ! tzz
      call rsgstauij(uz,uz,tij,g,nx,ny,nz) 
      
      enerdiss = enerdiss + cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      ! sxz
      do kk = 1, lz
      do ii = 1, lx1
        sij(ii,:,kk) = .5 * eye * (kx(ii) * uz(ii,:,kk) + kz(kk) * ux(ii,:,kk)) * g(ii,:,kk)
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      ! txz
      call rsgstauij(ux,uz,tij,g,nx,ny,nz) 
 
      enerdiss = enerdiss + 2. * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      ! syz
      do kk = 1, lz
      do jj = 1, ly
      sij(:,jj,kk) = .5 * eye * ( kz(kk) * uy(:,jj,kk) + ky(jj) * uz(:,jj,kk) ) * g(:,jj,kk)
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,sij,ignore_me)
 
      ! tyz
      call rsgstauij(uy,uz,tij,g,nx,ny,nz) 
 
      enerdiss = enerdiss + 2. * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )
 
      meanenerdiss(ll) = meanenerdiss(ll) &
          -(sum(real(enerdiss(1:lx,:,:)))+sum(aimag(enerdiss(1:lx,:,:))))/(nx*ny*nz)

    end do

    nfile = nfile + 1

  end do
  close(30)

  meanenerdiss = meanenerdiss / nfile

  open(15, file = 'meansgsed-'//strfltr(1:len_trim(strfltr))//'-'//filelist(1:len_trim(filelist))//'.dat')
    do ll = 1, numfilter
      delta_c=ndel*2*pi/nx/ 2**(ll-1)
      write(15,*) ndel/2**(ll-1), delta_c, meanenerdiss(ll)
    end do
  close(15)

  call destroyplan3d

  deallocate(ux,uy,uz, meanenerdiss)
  deallocate(tij,sij,enerdiss,kx,ky,kz,g,k2)

  write(*,*) 'meansgsenerdiss.x finished'
end program meansgsenerdiss      
