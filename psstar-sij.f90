#include "constant.f90"
#include "fftw3plan.f90"
#include "wavenumber.f90"

program psstartij
  use mconstant
  use mfftw3plan
  use mwavenumber
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s22,s33
  real(sp),    allocatable, dimension(:,:,:) :: g, kx, ky, kz

  integer(8) :: t11c2r, t22c2r, t33c2r, t12c2r, t23c2r, t13c2r

  real(sp), dimension(3,3) :: tij

  integer,  parameter :: npdf = 100
  real(sp), parameter :: binw = 2._sp/npdf
  real(dp), dimension(npdf) :: psstar

  real(sp) :: delta_c, invii, inviii
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> psstar-sij.x for DNS<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./psstar-sij.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  str='-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  delta_c=ndel*2*pi/nx

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  allocate(s11(lx1,ly,lz),s22(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call dfftwplan3dc2r(t11,nx,t11c2r)
  call dfftwplan3dc2r(t22,nx,t22c2r)
  call dfftwplan3dc2r(t33,nx,t33c2r)
  call dfftwplan3dc2r(t12,nx,t12c2r)
  call dfftwplan3dc2r(t13,nx,t13c2r)
  call dfftwplan3dc2r(t23,nx,t23c2r)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24._sp)

  nfile=0
  psstar = 0._dp
  open(20,file=fnm(1:len_trim(fnm)))

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s11
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s22
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s33
      close(10)

  
      t11 = eye * kx * s11 * g
      t22 = eye * ky * s22 * g
      t33 = eye * kz * s33 * g
      t12 = eye * (kx * s22 + ky * s11) * g
      t13 = eye * (kx * s33 + kz * s11) * g
      t23 = eye * (ky * s33 + kz * s22) * g
      
      call dfftw_execute(t11c2r)
      call dfftw_execute(t12c2r)
      call dfftw_execute(t13c2r)
      call dfftw_execute(t22c2r)
      call dfftw_execute(t23c2r)
      call dfftw_execute(t33c2r)
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx

        if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1)/2
          tij(1,1) = real(t11(ll,jj,kk),sp)
          tij(1,2) = real(t12(ll,jj,kk),sp)
          tij(1,3) = real(t13(ll,jj,kk),sp)
          tij(2,2) = real(t22(ll,jj,kk),sp)
          tij(2,3) = real(t23(ll,jj,kk),sp)
          tij(3,3) = real(t33(ll,jj,kk),sp)
        else
          ll = ii/2
          tij(1,1) = aimag(t11(ll,jj,kk))
          tij(1,2) = aimag(t12(ll,jj,kk))
          tij(1,3) = aimag(t13(ll,jj,kk))
          tij(2,2) = aimag(t22(ll,jj,kk))
          tij(2,3) = aimag(t23(ll,jj,kk))
          tij(3,3) = aimag(t33(ll,jj,kk))
        end if
        tij(2,1) = tij(1,2)
        tij(3,1) = tij(1,3)
        tij(3,2) = tij(2,3)

        invii = -0.5_sp * sum(tij * tij)
        tij = matmul(matmul(tij,tij),tij)
        inviii = (1._sp/3._sp) * (tij(1,1) + tij(2,2) + tij(3,3))

        invii = - 3 * sqrt(6._sp) * inviii / (-2 * invii)**1.5
        ll = floor((invii+1)/binw) + 1
   
        if ( ll .ge. 1 .and. ll .le. npdf) psstar(ll) = psstar(ll) + 1

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  psstar = psstar / (nfile * nx * ny * nz)
  write(*,*) 'check normalization psstar:', sum(psstar)
  psstar = psstar / binw

  open(15,file='psstar-sij'//str(1:len_trim(str)))
  do ii = 1, npdf
    write(15,'(14E12.4)') -1 + (ii-0.5)*binw, psstar(ii)
  end do
  close(15)


  deallocate(s11,s22,s33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz)


  write(*,*) 'psstar-sij.x finished.'

end program psstartij 
