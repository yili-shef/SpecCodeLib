program sgsedpdf2d
  use mconstant
  use mfftwplan2d
  use mwavenumber
  use msgstij
  implicit none

  integer :: lx1,lx,ly,nx,ny,ii,jj,ll,nfile,ndel

  integer,  parameter :: npnt = 150
  real(sp), parameter :: bound = 25., binw = 2.*bound/npnt
  real(dp), dimension(npnt) :: pdfed

  complex(sp), allocatable, dimension(:,:) :: ux, uy
  complex(sp), allocatable, dimension(:,:) :: sij, tij, enerdiss
  real(sp),    allocatable, dimension(:,:) :: g
  real(sp),    allocatable, dimension(:) :: kx, ky

  real(dp) :: meaned, rmsed
  real(sp) :: meaned0, rmsed0, ignore_me, const, delta_c

  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> PDF of SGS energy dissipation in 2D <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfsgsed2d-s(d)p.x nx filelist ndel'
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
  ny=nx
  lx=nx/2; ly=ny; lx1=lx+1
  const=1./(nx*ny)

  ! filter scale
  delta_c=ndel*2*pi/nx

  call fftwplan2de(nx,ny)
  write(*,*) 'after fftwplan2de'

  allocate( kx(lx1), ky(ly) )
  allocate( g(lx1,ly), enerdiss(lx1,ly) )
  allocate( ux(lx1,ly), uy(lx1,ly) )
  allocate( sij(lx1,ly), tij(lx1,ly) )
  write(*,*) 'allocated'

  call wavenumber(kx,ky,g,lx1,ly)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

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
      write(*,*) 'after reading data files'

      call sgstij( ux, ux, tij, g, nx, ny )
      do ii = 1, lx1
      sij(ii,:) = eye * kx(ii) * ux(ii,:) * g(ii,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r2d,sij,ignore_me)

      enerdiss = cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      call sgstij( uy, uy, tij, g, nx, ny )
      do jj = 1, ly
      sij(:,jj) = eye * ky(jj) * uy(:,jj) * g(:,jj)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r2d,sij,ignore_me)

      enerdiss = enerdiss + cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      call sgstij( ux, uy, tij, g, nx, ny )
      do jj = 1, ly
      do ii = 1, lx1
        sij(ii,jj) = .5 * eye * ( kx(ii) * uy(ii,jj) + ky(jj) * ux(ii,jj) ) * g(ii,jj)
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r2d,sij,ignore_me)

      enerdiss = enerdiss + 2. * cmplx( real(sij)*real(tij), aimag(sij)*aimag(tij) )

      enerdiss = - enerdiss

      if (nfile .eq. 0) then 
        meaned0 = sum(real(enerdiss(1:lx,:))) + sum(aimag(enerdiss(1:lx,:)))
        meaned0 = meaned0 * const
        rmsed0 = sum(real(enerdiss(1:lx,:)**2)) + sum(aimag(enerdiss(1:lx,:)**2)) 
        rmsed0 = rmsed0 * const
        rmsed0 = sqrt(rmsed0 - meaned0**2)
        write(*,*) 'estimate meaned0 = ', meaned0
        write(*,*) 'estimate rmsed0  = ', rmsed0
      end if

      do jj = 1, ny
      do ii = 1, nx

        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1)/2
          ignore_me = real(enerdiss(ll,jj))
        else
          ll = ii / 2
          ignore_me = aimag(enerdiss(ll,jj))
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

      nfile = nfile + 1
    end do
  close(20)
  pdfed = pdfed * const / nfile
  write(*,*) 'check pdfed: ', sum(pdfed)
  pdfed = pdfed / binw

  meaned = meaned * const / nfile
  rmsed = rmsed * const / nfile

  rmsed = sqrt( rmsed - meaned * meaned )

  write(*,*) 'Calculated mean enerdiss: ', meaned
  write(*,*) 'Calculated rms of enerdiss:', rmsed

  open(15, file = '2dpdfed'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat')
    write(15,*) 'Title = " meaned is', meaned, 'rmsed is', rmsed, '"'
    write(15,*) 'variables = "(pi_H-<pi_H>)/rms", "PDF"'
    do ii = 1, npnt
      write(15,*) ( (-bound + (ii-.5) * binw) * rmsed0 + meaned0 - meaned ) / rmsed, & 
                  pdfed(ii) * rmsed / rmsed0
    end do
  close(15)

  deallocate( kx, ky, g, ux, uy, sij, tij, enerdiss )

  call destroyplan2d

  write(*,*) 'finished'

end program sgsedpdf2d
