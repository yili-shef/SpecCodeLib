program meansgsmhdiss
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ndel,numfilter,ii,ll,nfile

  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wa,wb
  complex, allocatable, dimension(:,:,:) :: hpx, hpy, hpz, ap
  complex, allocatable, dimension(:,:,:) :: tij, rij
  real,    allocatable, dimension(:,:,:) :: k2, g, kx, ky, kz
  real,    allocatable, dimension(:)     :: helidiss

  real :: delta_c, ignore_me
  character(80) :: str,str1,fnm

  write(*,*) 
  write(*,'(''>>> Mean sgs helicity dissipation by negative helical modes<<< '')')
  write(*,*)

  ll = iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgsmhdiss.x nx filelist ndel numfilter'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: list of data files'
          write(*,*) '        ndel: maximum filter scale delta = ndel * dx'
          write(*,*) '        numfilter: number of scales, each halved from a larger one'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if 

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  call getarg(2,str1) 
  str1=adjustl(str1)
  fnm = str1(1:len_trim(str1))//'.list'
  str1 = 'mhdiss-'//str1(1:len_trim(str1))//'.dat'

  call getarg(3,str)
  read(str, '(I20)') ndel

  call getarg(4,str)
  read(str, '(I20)') numfilter

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wa(lx1,ly,lz),wb(lx1,ly,lz),ap(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  allocate(rij(lx1,ly,lz),tij(lx1,ly,lz))
  allocate(g(lx1,ly,lz),helidiss(numfilter))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx(:,1,1),ky(1,:,1),kz(1,1,:),k2,lx1,ly,lz)
  hpx=conjg(hpx)
  hpy=conjg(hpy)
  hpz=conjg(hpz)
  ! Make it negative helical wave
  write(*,*) 'after heliwave'

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  open(20,file=fnm(1:len_trim(fnm)))

    helidiss=0.
    nfile=0
    do while ( .not. eof(20))
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      open(10,file='./out/ux'//str(1:len_trim(str)),status='unknown',form='unformatted')
        read(10)ux
      close(10)
      open(10,file='./out/uy'//str(1:len_trim(str)),status='unknown',form='unformatted')
        read(10)uy
      close(10)
      open(10,file='./out/uz'//str(1:len_trim(str)),status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'
  
      ap = ux * conjg(hpx) + uy * conjg(hpy) + uz * conjg(hpz)
      ap = ap * sqrt(k2) 
  
      do ii=1,numfilter
  
        delta_c=ndel*2*pi/nx/ 2**(ii-1)
  
        ! Gaussian filter
        g=exp(-k2*delta_c**2/24.)
  
  
        ! wa = wxp; wb = wyp
        wa = ap * hpx * g
        wb = ap * hpy * g
  
        
        ! rxx
        rij = eye * kx * wa
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        ! txx
        call rsgstauij(ux,ux,tij,g,nx,ny,nz) 
  
        helidiss(ii) = sum( real(rij)*real(tij) ) + sum( aimag(rij)*aimag(tij) )
  
        ! tyy
        call rsgstauij(uy,uy,tij,g,nx,ny,nz) 
  
        ! ryy
        rij = eye * ky * wb
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        helidiss(ii) = helidiss(ii) + sum( real(rij)*real(tij) )+ sum( aimag(rij)*aimag(tij) )
  
        ! txy
        call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
  
        ! rxy
        rij = .5 * eye * (kx * wb + ky * wa)
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        helidiss(ii) = helidiss(ii) + 2.*sum( real(rij)*real(tij)) +2. *  sum( aimag(rij)*aimag(tij) )
  
  
        ! wa = wx; wb = wz
        wb  = ap * hpz * g
  
        ! rzz
        rij = eye * kz * wb 
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        ! tzz
        call rsgstauij(uz,uz,tij,g,nx,ny,nz) 
        
        helidiss(ii) = helidiss(ii) + sum( real(rij)*real(tij) ) + sum( aimag(rij)*aimag(tij) )
  
        ! rxz
        rij = .5 * eye * (kx * wb + kz * wa)
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        ! txz
        call rsgstauij(ux,uz,tij,g,nx,ny,nz) 
  
        helidiss(ii) = helidiss(ii) + 2. * sum( real(rij)*real(tij)) + 2. * sum( aimag(rij)*aimag(tij))
  
        ! wa = wy; wb=wz
        wa = ap * hpy * g
  
        ! ryz
        rij = .5 * eye * ( kz * wa + ky * wb )
        call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
  
        ! tyz
        call rsgstauij(uy,uz,tij,g,nx,ny,nz) 
  
        helidiss(ii) = helidiss(ii) + 2. * sum( real(rij)*real(tij)) +2. *  sum( aimag(rij)*aimag(tij))
  
      end do

      nfile = nfile + 1
    end do

    helidiss = -2. * helidiss / (nx*ny*nz*nfile)

  close(20)

  open(20, file = str1)
    do ii = 1, numfilter
      write(20,*) ndel/2**(ii-1), ndel*2*pi/nx/ 2**(ii-1), helidiss(ii)
    end do
  close(20)


  call destroyplan3d

  deallocate(ux,uy,uz,wa,wb)
  deallocate(tij,rij,helidiss,k2,kx,ky,kz,g)
  deallocate(hpx,hpy,hpz,ap)

  write(*,*) 'finished'

end program meansgsmhdiss 
