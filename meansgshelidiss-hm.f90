program meansgshelidiss
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ndel,numfilter,ll,nfile

  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz,fux,fuy,fuz
  complex(sp), allocatable, dimension(:,:,:) :: tij, rij
  real(sp),    allocatable, dimension(:,:,:) :: k2, g, kx, ky, kz
  real(dp),    allocatable, dimension(:)     :: meanhelidiss

  real(sp) :: delta_c, ignore_me
  real(dp) :: const
  character(80) :: str,str1, filelist

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
          write(*,*) '        filelist: filelist.list is the lis tof data files'
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
  const = 1._dp/(nx*ny*nz)

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(fux(lx1,ly,lz),fuy(lx1,ly,lz),fuz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz), kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))
  allocate(rij(lx1,ly,lz), tij(lx1,ly,lz))
  allocate(g(lx1,ly,lz), meanhelidiss(numfilter))
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
 
    wx = ux; wy = uy; wz = uz
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
 
    do ll=1,numfilter
 
      delta_c=ndel*2*pi/nx/ 2**(ll-1)
 
      ! Gaussian filter
      g=exp(-k2*delta_c**2/24.)
 
      fux = wx * g; fuy = wy * g; fuz = wz * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,fux,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,fuy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,fuz,ignore_me)
 
      rij = eye * kx * eye * (ky * wz - kz * wy ) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = ux * ux * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fux * fux
 
      ignore_me = sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
      rij = eye * ky * eye * ( kz * wx - kx * wz) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = uy * uy * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fuy * fuy
 
      ignore_me = sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
      rij = eye * kz * eye * ( kx * wy - ky * wx ) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = uz * uz * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fuz * fuz
 
      ignore_me = sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
      rij = .5 * eye * ( kx * eye * ( kz * wx - kx * wz) + ky * eye * (ky * wz - kz * wy) ) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = ux * uy * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fux * fuy
 
      ignore_me = 2 * sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + 2 * sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
      rij = .5 * eye * ( kx * eye * (kx * wy - ky * wx) + kz * eye * (ky * wz - kz * wy) ) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = ux * uz * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fux * fuz
 
      ignore_me = 2 * sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + 2 * sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
      rij = .5 * eye * ( ky * eye * ( kx * wy - ky * wx ) + kz * eye * ( kz * wx - kx * wz) ) * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)
      tij = uy * uz * const
      call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
      tij = tij * g
      call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)
      tij = tij - fuy * fuz
 
      ignore_me = 2 * sum( real(rij(1:lx,:,:))*real(tij(1:lx,:,:)) )
      meanhelidiss(ll) = meanhelidiss(ll) + ignore_me + 2 * sum( aimag(rij(1:lx,:,:))*aimag(tij(1:lx,:,:)) )
 
    end do

    nfile = nfile + 1

  end do
  close(30)

  meanhelidiss = -2 * meanhelidiss / nfile / (nx * ny * nz)

  open(15, file = 'meansgshd-'//filelist(1:len_trim(filelist))//'.dat')
    do ll = 1, numfilter
      delta_c=ndel*2*pi/nx/ 2**(ll-1)
      write(15,*) ndel/2**(ll-1), delta_c, meanhelidiss(ll)
    end do
  close(15)


  call destroyplan3d

  deallocate(ux, uy, uz, wx, wy, wz, fux, fuy, fuz)
  deallocate(tij,rij,k2,kx,ky,kz,g)

  write(*,*) 'finished'

end program meansgshelidiss      
