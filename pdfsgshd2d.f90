program sgshdpdf2d
  use mconstant
  use mfftwplan2d
  use mwavenumber
  use msgstij
  implicit none

  integer :: lx1,lx,ly,nx,ny,ii,jj,ll,nfile,ndel

  integer,  parameter :: npnt = 150
  real(sp), parameter :: bound = 25., binw = 2.*bound/npnt
  real(dp), dimension(npnt) :: pdfhd

  complex(sp), allocatable, dimension(:,:) :: ux, uy, wz
  complex(sp), allocatable, dimension(:,:) :: gradw, tij, enstdiss
  real(sp),    allocatable, dimension(:,:) :: g
  real(sp),    allocatable, dimension(:) :: kx, ky

  real(dp) :: meanhd, rmshd
  real(sp) :: meanhd0, rmshd0, ignore_me, const, delta_c

  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> PDF of SGS energy dissipation in 2D <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfsgshd2d-s(d)p.x nx filelist ndel'
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
  allocate( g(lx1,ly), enstdiss(lx1,ly) )
  allocate( ux(lx1,ly), uy(lx1,ly), wz(lx1,ly) )
  allocate( gradw(lx1,ly), tij(lx1,ly) )
  write(*,*) 'allocated'

  call wavenumber(kx,ky,g,lx1,ly)
  write(*,*) 'after wavenumber'

  ! Gaussian filter in 2D TODO numerical factor need change
  g=exp(-g*delta_c**2/24.)

  nfile = 0
  meanhd = 0.0_dp
  rmshd = 0.0_dp
  pdfhd = 0.0_dp
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

      do jj = 1, ly
      do ii = 1, lx1
        wz(ii,jj) = eye * ( kx(ii) * uy(ii,jj) - ky(jj) * ux(ii,jj) )
        gradw(ii,jj) = eye * kx(ii) * wz(ii,jj) * g(ii,jj)
      end do
      end do
      call sgstij( ux, wz, tij, g, nx, ny )
      call rfftwnd_f77_one_complex_to_real(c2r2d,gradw,ignore_me)
      
      enstdiss = cmplx( real(gradw)*real(tij), aimag(gradw)*aimag(tij) )

      do jj = 1, ly
      do ii = 1, lx1
        gradw(ii,jj) = eye * ky(jj) * wz(ii,jj) * g(ii,jj)
      end do
      end do
      call sgstij( uy, wz, tij, g, nx, ny )
      call rfftwnd_f77_one_complex_to_real(c2r2d,gradw,ignore_me)

      enstdiss = enstdiss + cmplx( real(gradw)*real(tij), aimag(gradw)*aimag(tij) )
      enstdiss = - enstdiss

      if (nfile .eq. 0) then 
        meanhd0 = sum(real(enstdiss(1:lx,:))) + sum(aimag(enstdiss(1:lx,:)))
        meanhd0 = meanhd0 * const
        rmshd0 = sum(real(enstdiss(1:lx,:)**2)) + sum(aimag(enstdiss(1:lx,:)**2)) 
        rmshd0 = rmshd0 * const
        rmshd0 = sqrt(rmshd0 - meanhd0**2)
        write(*,*) 'estimate meanhd0 = ', meanhd0
        write(*,*) 'estimate rmshd0  = ', rmshd0
      end if

      do jj = 1, ny
      do ii = 1, nx

        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1)/2
          ignore_me = real(enstdiss(ll,jj))
        else
          ll = ii / 2
          ignore_me = aimag(enstdiss(ll,jj))
        end if

        meanhd = meanhd + ignore_me
        rmshd = rmshd + ignore_me * ignore_me 

        ignore_me = (ignore_me - meanhd0) / rmshd0
        ll = floor( ( ignore_me + bound ) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npnt ) then
          pdfhd(ll) = pdfhd(ll) + 1
        end if

      end do
      end do

      nfile = nfile + 1
    end do
  close(20)
  pdfhd = pdfhd * const / nfile
  write(*,*) 'check pdfhd: ', sum(pdfhd)
  pdfhd = pdfhd / binw

  meanhd = meanhd * const / nfile
  rmshd = rmshd * const / nfile
  rmshd = sqrt( rmshd - meanhd * meanhd )

  write(*,*) 'Calculated mean enerdiss: ', meanhd
  write(*,*) 'Calculated rms of enerdiss:', rmshd

  open(15, file = '2dpdfhd'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat')
    write(15,*) 'Title = " meanhd is', meanhd, 'rmshd is', rmshd, '"'
    write(15,*) 'variables = "(pi_H-<pi_H>)/rms", "PDF"'
    do ii = 1, npnt
      write(15,*) ( (-bound + (ii-.5) * binw) * rmshd0 + meanhd0 - meanhd ) / rmshd, & 
                  pdfhd(ii) * rmshd / rmshd0
    end do
  close(15)

  deallocate( kx, ky, g, ux, uy, gradw, tij, enstdiss )

  call destroyplan2d

  write(*,*) 'finished'

end program sgshdpdf2d
