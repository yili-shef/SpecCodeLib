#include "constant.f90"
#include "fftw3plan.f90"
#include "wavenumber.f90"
#include "rsgstauii.f90"
#include "rsgstauij.f90"

program msgsed
  use mconstant
  use mwavenumber
  use mfftw3plan
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, ndel, nfile, npnts

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33, vx, vy, vz
  complex(sp), allocatable, dimension(:,:,:) :: t11, t12, t13, t22, t23, t33, pie
  real(sp),    allocatable, dimension(:,:,:) :: g, kx, ky, kz, k2

  real(dp), allocatable, dimension(:) :: menerdiss
  real(sp) :: delta_c, const

  character(80) :: str, str1, dnslist, fpath
  integer(8) :: s11c2r, s12c2r, s13c2r, s22c2r, s23c2r, s33c2r
  integer(8) :: s12r2c, s13r2c

  write(*,*) 
  write(*,'(''>>>>>> mean SGS energy dissipation as function of delta <<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 2) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./meansgsed.x nx datalist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        datalist: list for dns data files'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! list of dns data files
  call getarg(2,dnslist)
  dnslist = adjustl(dnslist)

  str = '-'//dnslist(1:len_trim(dnslist))//'.dat'
  dnslist = dnslist( 1:len_trim(dnslist) )//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1
  const=1._sp/(nx*ny*nz)
  npnts = nx/2

  allocate(menerdiss(npnts))
  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),k2(lx1,ly,lz),pie(lx1,ly,lz))
  allocate(vx(lx1,ly,lz),vy(lx1,ly,lz),vz(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call dfftwplan3dc2r(s11, nx, s11c2r)
  call dfftwplan3dc2r(s12, nx, s12c2r)
  call dfftwplan3dc2r(s13, nx, s13c2r)
  call dfftwplan3dc2r(s22, nx, s22c2r)
  call dfftwplan3dc2r(s23, nx, s23c2r)
  call dfftwplan3dc2r(s33, nx, s33c2r)

  call dfftwplan3dr2c(s12, nx, s12r2c)
  call dfftwplan3dr2c(s13, nx, s13r2c)
  write(*,*) 'after fftwplan3d'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(20,file=dnslist(1:len_trim(dnslist)))
    menerdiss = 0._dp
    nfile = 0
    do while ( .not. eof(20) )

      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)vx
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)vy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)vz
      close(10)

      do ndel = 1, npnts

        delta_c = ndel*2*pi/nx

        ! Gaussian filter
        g=exp(-k2*delta_c**2/24._sp)
       
        s12 = vx
        call rsgstauii(s12,s12r2c,s12c2r,t11,g,nx,ny,nz) 
        s12 = vy
        call rsgstauii(s12,s12r2c,s12c2r,t22,g,nx,ny,nz) 
        s12 = vz
        call rsgstauii(s12,s12r2c,s12c2r,t33,g,nx,ny,nz) 
        s12 = vx; s13 = vy
        call rsgstauij(s12,s13,s12r2c,s12c2r,s13r2c,s13c2r,t12,g,nx,ny,nz) 
        s12 = vx; s13 = vz
        call rsgstauij(s12,s13,s12r2c,s12c2r,s13r2c,s13c2r,t13,g,nx,ny,nz) 
        s12 = vy; s13 = vz
        call rsgstauij(s12,s13,s12r2c,s12c2r,s13r2c,s13c2r,t23,g,nx,ny,nz) 
 
        s12 = -(t11+t22+t33) / 3._sp
 
        t11 = t11 + s12
        t22 = t22 + s12
        t33 = t33 + s12
    
        s11=vx*g
        s22=vy*g 
        s33=vz*g
        s12 = .5_sp * eye * ( kx * s22 + ky * s11 )
        s13 = .5_sp * eye * ( kx * s33 + kz * s11 )
        s23 = .5_sp * eye * ( ky * s33 + kz * s22 )
        s11 = eye * kx * s11
        s22 = eye * ky * s22
        s33 = eye * kz * s33
 
        call dfftw_execute(s11c2r)
        call dfftw_execute(s12c2r)
        call dfftw_execute(s13c2r)
        call dfftw_execute(s22c2r)
        call dfftw_execute(s23c2r)
        call dfftw_execute(s33c2r)
 
        pie = cmplx(  real(t11) *  real(s11)           +  real(t22) *  real(s22) & 
                   +  real(t33) *  real(s33) + 2._sp * (  real(t12) *  real(s12) & 
                   +  real(t13) *  real(s13)           +  real(t23) *  real(s23) ), &
                     aimag(t11) * aimag(s11)           + aimag(t22) * aimag(s22) & 
                   + aimag(t33) * aimag(s33) + 2._sp * ( aimag(t12) * aimag(s12) & 
                   + aimag(t13) * aimag(s13)           + aimag(t23) * aimag(s23) ) )
        pie = - pie
 
        menerdiss(ndel) = menerdiss(ndel) + sum( real(pie(1:lx,:,:))) &
                                          + sum(aimag(pie(1:lx,:,:))) 
 
      end do

      nfile = nfile + 1
    end do
  close(20) 

  menerdiss = menerdiss / nfile / (nx*ny*nz)

  open( 15, file ='msgsed-dns-vs-delta'//str(1:len_trim(str)) )
    do ii = 1, npnts
      write(15,'(20E18.5)') ii*2*pi/nx, menerdiss(ii)
    end do
  close(15)

  deallocate(kx,ky,kz,k2,g,vx,vy,vz,pie, menerdiss)
  deallocate(s11,s12,s13,s22,s23,s33,t11,t12,t13,t22,t23,t33)

  write(*,*) 'meansgsed.x finished'

end program msgsed
