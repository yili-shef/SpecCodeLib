#include "constant.f90"
#include "wavenumber.f90"
#include "fftw3plan.f90"

#include "rsgstauii.f90"
#include "rsgstauij.f90"

program sgsedpdf
  use mconstant
  use mfftw3plan
  use mwavenumber
  implicit none

  integer :: lx1,lx,ly,lz,nx,ny,nz,ii,jj,kk,ll,nfile,ndel

  integer(8) :: sijc2r, sijr2c, tmpc2r, tmpr2c

  integer,  parameter :: npnt = 150
  real(sp), parameter :: bound = 25., binw = 2.*bound/npnt
  real(sp), dimension(npnt) :: pdfed

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, tmp 
  complex(sp), allocatable, dimension(:,:,:) :: sij, tij, enerdiss
  real(sp),    allocatable, dimension(:,:,:) :: g, kx, ky, kz

  real(sp) :: meaned, rmsed
  real(sp) :: meaned0, rmsed0, ignore_me, const, delta_c

  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> PDF of SGS energy dissipation<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfsgsed-s(d)p.x nx filelist ndel'
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

  ! some parameters 
  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  ! filter scale
  delta_c=ndel*2*pi/nx

  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
  allocate( g(lx1,ly,lz), enerdiss(lx1,ly,lz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( sij(lx1,ly,lz), tij(lx1,ly,lz), tmp(lx1,ly,lz) )
  write(*,*) 'allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  call dfftwplan3dc2r(sij, nx, sijc2r)
  call dfftwplan3dr2c(sij, nx, sijr2c)
  call dfftwplan3dc2r(tmp, nx, tmpc2r)
  call dfftwplan3dr2c(tmp, nx, tmpr2c)

  nfile = 0
  meaned = 0.0_dp
  rmsed = 0.0_dp
  pdfed = 0.0_dp
  open( 20, file = fnm( 1:len_trim(fnm) )//'.list' )
    do while ( .not. eof(20) )
      read(20,*) str1
      write(*,*) str1( 1:len_trim(str1) )

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)ux
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'

      sij = ux
      call rsgstauii( sij, sijr2c, sijc2r, tij, g, nx, ny, nz )
      sij = eye * kx * ux * g
      call dfftw_execute(sijc2r)
      enerdiss = cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      sij = uy
      call rsgstauii( sij, sijr2c, sijc2r, tij, g, nx, ny, nz )
      sij = eye * ky * uy * g
      call dfftw_execute(sijc2r)
      enerdiss = enerdiss + cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      sij = ux; tmp = uy
      call rsgstauij( sij, tmp, sijr2c, sijc2r, tmpr2c, tmpc2r, tij, g, nx, ny, nz )
      sij = .5_sp * eye * ( kx * uy + ky * ux ) * g
      call dfftw_execute(sijc2r)
      enerdiss = enerdiss + 2._sp * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      sij = ux; tmp = uz
      call rsgstauij( sij, tmp, sijr2c, sijc2r, tmpr2c, tmpc2r, tij, g, nx, ny, nz )
      sij = .5_sp * eye * ( kx * uz + kz * ux ) * g
      call dfftw_execute(sijc2r)
      enerdiss = enerdiss + 2._sp * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      sij = uy; tmp = uz
      call rsgstauij( sij, tmp, sijr2c, sijc2r, tmpr2c, tmpc2r, tij, g, nx, ny, nz )
      sij = .5_sp * eye * ( kz * uy + ky * uz ) * g
      call dfftw_execute(sijc2r)
      enerdiss = enerdiss + 2._sp * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      enerdiss = - enerdiss

      if (nfile .eq. 0) then 
        meaned0 = sum(real(enerdiss(1:lx,:,:))) + sum(aimag(enerdiss(1:lx,:,:)))
        meaned0 = meaned0 * const
        rmsed0 = sum(real(enerdiss(1:lx,:,:)**2)) + sum(aimag(enerdiss(1:lx,:,:)**2)) 
        rmsed0 = rmsed0 * const
        rmsed0 = sqrt(rmsed0 - meaned0**2)
        write(*,*) 'estimate meaned0 = ', meaned0
        write(*,*) 'estimate rmsed0  = ', rmsed0
      end if

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx

        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1)/2
          ignore_me = real(enerdiss(ll,jj,kk))
        else
          ll = ii / 2
          ignore_me = aimag(enerdiss(ll,jj,kk))
        end if

        meaned = meaned + ignore_me
        rmsed = rmsed + ignore_me * ignore_me 

        ignore_me = (ignore_me - meaned0) / rmsed0
        ll = floor( ( ignore_me + bound ) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npnt ) then
          pdfed(ll) = pdfed(ll) + 1
        end if

      end do
      end do
      end do

      nfile = nfile + 1

    end do
  close(20)
  pdfed = pdfed * const / nfile
  write(*,*) 'check pdfed: ', sum(pdfed)
  pdfed = pdfed / binw

  meaned = meaned * const / nfile
  rmsed = rmsed * const / nfile

  rmsed = sqrt( rmsed - meaned * meaned )

  open(15, file = 'pdfed'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat')
    write(15,*) '# Calculated mean enerdiss: ', meaned
    write(15,*) '# Calculated rms of enerdiss:', rmsed

    write(15,"('# Title = meaned is', E12.4, 'rmsed is', E12.4)") meaned, rmsed
    write(15,*) '# variables = "(pi_H-<pi_H>)/rms", "PDF"'
    do ii = 1, npnt
      write(15,*) ( (-bound + (ii-.5) * binw) * rmsed0 + meaned0 - meaned ) / rmsed, & 
                  pdfed(ii) * rmsed / rmsed0
    end do
  close(15)

  deallocate( kx, ky, kz, g, ux, uy, uz, tmp, sij, tij, enerdiss )

  write(*,*) 'pdfsgsed.x finished'

end program sgsedpdf
