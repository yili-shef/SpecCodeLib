#include "constant.f90"
#include "wavenumber.f90"
#include "fftw3plan.f90"
!-------------------------
! subroutines for 
! eigenvalue calculations
#include "pythag.f90"
#include "rs.f90"
#include "tql1.f90"
#include "tql2.f90"
#include "tred1.f90"
#include "tred2.f90"
!-------------------------

program alignos
  use mconstant
  use mwavenumber
  use mfftw3plan
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  integer(8) :: wxc2r, wyc2r, wzc2r, s11c2r, s12c2r, s13c2r, s22c2r, s23c2r, s33c2r

  complex(sp), allocatable, dimension(:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  !---- declaration needed to call subroutine rs ---
  ! rs calculates the eigenvectors and eigenvalues for input array
  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer, parameter :: npnt=80
  real(sp), parameter :: binw=1._sp/npnt, binwphi=pi/2/npnt
  real(dp), dimension(npnt) :: palpha, pbeta, pgamma
  real(dp), dimension(npnt, npnt) :: p2d

  real(sp) :: delta_c, tmp1, tmp2
  real(dp) :: alpha, beta, gamm
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between xi and eigenvectors of sij <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-os.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: dns data file list'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! filter parameter
  call getarg(3,str)
  read(str, '(I20)') ndel
  ! file list string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  ! filter length scale
  delta_c=ndel*2*pi/nx

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call dfftwplan3dc2r(wx, nx, wxc2r)
  call dfftwplan3dc2r(wy, nx, wyc2r)
  call dfftwplan3dc2r(wz, nx, wzc2r)
        
  call dfftwplan3dc2r(s11, nx, s11c2r)
  call dfftwplan3dc2r(s12, nx, s12c2r)
  call dfftwplan3dc2r(s13, nx, s13c2r)
  call dfftwplan3dc2r(s22, nx, s22c2r)
  call dfftwplan3dc2r(s23, nx, s23c2r)
  call dfftwplan3dc2r(s33, nx, s33c2r)
  write(*,*) 'after fftwplan3d'

  ! g is used to store k2 temporarily
  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  palpha=0._dp
  pbeta=0._dp
  pgamma=0._dp
  p2d = 0._dp
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s11 ! s11 is used to store ux
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s22 ! s22 is used to store uy
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s33 ! s33 is used to store uz
    close(10)
    write(*,*) 'after reading data files'
 
    s11=s11*g
    s22=s22*g 
    s33=s33*g
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      ! vorticity in Fourier space
      wx(ii,jj,kk) = eye * ( ky(jj) * s33(ii,jj,kk) - kz(kk) * s22(ii,jj,kk) )
      wy(ii,jj,kk) = eye * ( kz(kk) * s11(ii,jj,kk) - kx(ii) * s33(ii,jj,kk) )
      wz(ii,jj,kk) = eye * ( kx(ii) * s22(ii,jj,kk) - ky(jj) * s11(ii,jj,kk) )

      ! strain rate tensor is Fourier space
      s12(ii,jj,kk)=.5_sp*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
      s13(ii,jj,kk)=.5_sp*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
      s23(ii,jj,kk)=.5_sp*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))
      s11(ii,jj,kk)=eye*kx(ii)*s11(ii,jj,kk)
      s22(ii,jj,kk)=eye*ky(jj)*s22(ii,jj,kk)
      s33(ii,jj,kk)=eye*kz(kk)*s33(ii,jj,kk)
    end do
    end do
    end do

 
    call dfftw_execute(wxc2r)
    call dfftw_execute(wyc2r)
    call dfftw_execute(wzc2r)

    call dfftw_execute(s11c2r)
    call dfftw_execute(s12c2r)
    call dfftw_execute(s13c2r)
    call dfftw_execute(s22c2r)
    call dfftw_execute(s23c2r)
    call dfftw_execute(s33c2r)

    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
 
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii+1)/2
 
        delta_c   = real( wx(ll, jj, kk), sp )
        tmp2 = real( wy(ll, jj, kk), sp )
        tmp1       = real( wz(ll, jj, kk), sp)
     
        cc(1,1)=real(s11(ll,jj,kk), sp)
        cc(1,2)=real(s12(ll,jj,kk), sp)
        cc(1,3)=real(s13(ll,jj,kk), sp)
        cc(2,2)=real(s22(ll,jj,kk), sp)
        cc(2,3)=real(s23(ll,jj,kk), sp)
        cc(3,3)=real(s33(ll,jj,kk), sp)
 
      else
        ll = ii/2

        delta_c   = aimag( wx(ll, jj, kk) )
        tmp2      = aimag( wy(ll, jj, kk) )
        tmp1      = aimag( wz(ll, jj, kk) )
  
        cc(1,1)=aimag(s11(ll,jj,kk))
        cc(1,2)=aimag(s12(ll,jj,kk))
        cc(1,3)=aimag(s13(ll,jj,kk))
        cc(2,2)=aimag(s22(ll,jj,kk))
        cc(2,3)=aimag(s23(ll,jj,kk))
        cc(3,3)=aimag(s33(ll,jj,kk))

      end if
      cc(2,1)=cc(1,2)
      cc(3,1)=cc(1,3)
      cc(3,2)=cc(2,3)

      ! rs calculate the eigenvalues and eigenvectors for the input array
      ! In this case, cc is the input array. eigenvalues are returned in array
      ! 'evalues'. eigenvectors are returned in array 'evectors'. Each column
      ! of evectors corresponds to an eigenvector.
      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
 
      alpha = delta_c*evectors(1,3)+tmp2*evectors(2,3)+tmp1*evectors(3,3)
      beta  = delta_c*evectors(1,2)+tmp2*evectors(2,2)+tmp1*evectors(3,2)
      gamm  = delta_c*evectors(1,1)+tmp2*evectors(2,1)+tmp1*evectors(3,1)
 
      delta_c = sqrt( delta_c * delta_c + tmp2 * tmp2 + tmp1 * tmp1 )
 
      alpha = abs(alpha) / delta_c
      beta  = abs(beta) / delta_c
      gamm  = abs(gamm) / delta_c
      
      ll=floor((alpha)/binw)+1
      mm=floor((beta)/binw)+1
      nn=floor((gamm)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt .and. &
          mm .ge. 1 .and. mm .le. npnt .and. &
          nn .ge. 1 .and. nn .le. npnt)  then 
 
        palpha(ll)=palpha(ll)+1
        pbeta(mm)=pbeta(mm)+1
        pgamma(nn)=pgamma(nn)+1
      end if

      beta = abs( beta ) / sqrt( beta * beta + gamm * gamm )
      beta = acos( beta )
      
      mm = floor( beta / binwphi ) + 1
      if ( mm .ge. 1 .and. mm .le. npnt .and.  &
           ll .ge. 1 .and. ll .le. npnt ) then
           p2d( ll, mm ) = p2d( ll, mm ) + 1
      end if
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  alpha = 1._dp/(nx*ny*nz)/nfile

  write(*,*) 'check palpha:', sum(palpha)*alpha
  write(*,*) 'check pbeta:',  sum(pbeta)*alpha
  write(*,*) 'check pgamma:', sum(pgamma)*alpha
  write(*,*) 'check p2d:   ', sum(p2d)*alpha

  palpha = palpha * alpha/binw
  pbeta  = pbeta  * alpha/binw
  pgamma = pgamma * alpha/binw

  p2d = p2d * alpha / ( binw * binwphi )


  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='posalign-1d'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  open(15, file = 'posalign-2d'//fnm(1 : len_trim( fnm )))
    write(15,'(''zone t="2d align pdf", i='',I4,'', j='', I4, '', f=point'')') npnt,npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.6)') (ii-.5)*binw, (jj-.5)*binwphi, p2d(ii,jj)
    end do
    end do
  close(15)

  deallocate(wx,wy,wz,g,kx,ky,kz)
  deallocate(s11,s12,s13,s22,s23,s33)

  write(*,*) 'align-os.x done'

end program alignos
