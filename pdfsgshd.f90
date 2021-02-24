program sgshdpdf
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  integer :: lx1,lx,ly,lz,nx,ny,nz,ii,jj,kk,ll,nfile,ndel

  integer,  parameter :: npnt = 150
  real(sp), parameter :: bound = 25., binw = 2.*bound/npnt
  real(dp), dimension(npnt) :: pdfhd

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: rij, tij, helidiss
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:) :: kx, ky, kz

  real(dp) :: meanhd, rmshd
  real(sp) :: meanhd0, rmshd0, ignore_me, const, delta_c

  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> PDF of SGS helicity dissipation<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfsgshd.x nx filelist ndel'
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

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3de'

  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), helidiss(lx1,ly,lz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( rij(lx1,ly,lz), tij(lx1,ly,lz) )
  write(*,*) 'allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
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
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk) = eye * ( ky(jj) * uz(ii,jj,kk) - kz(kk) * uy(ii,jj,kk) )
        wy(ii,jj,kk) = eye * ( kz(kk) * ux(ii,jj,kk) - kx(ii) * uz(ii,jj,kk) )
        wz(ii,jj,kk) = eye * ( kx(ii) * uy(ii,jj,kk) - ky(jj) * ux(ii,jj,kk) )
      end do
      end do
      end do

      wx = wx * g
      wy = wy * g
      wz = wz * g

      call rsgstauij( ux, ux, tij, g, nx, ny, nz )
      do ii = 1, lx1
      rij(ii,:,:) = eye * kx(ii) * wx(ii,:,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      call rsgstauij( uy, uy, tij, g, nx, ny, nz )
      do jj = 1, ly
      rij(:,jj,:) = eye * ky(jj) * wy(:,jj,:)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      call rsgstauij( uz, uz, tij, g, nx, ny, nz )
      do kk = 1, lz
      rij(:,:,kk) = eye * kz(kk) * wz(:,:,kk)
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = helidiss + cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      call rsgstauij( ux, uy, tij, g, nx, ny, nz )
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        rij(ii,jj,kk) = .5 * eye * ( kx(ii) * wy(ii,jj,kk) + ky(jj) * wx(ii,jj,kk) )
      end do
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      call rsgstauij( ux, uz, tij, g, nx, ny, nz )
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        rij(ii,jj,kk) = .5 * eye * ( kx(ii) * wz(ii,jj,kk) + kz(kk) * wx(ii,jj,kk) )
      end do
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )

      call rsgstauij( uy, uz, tij, g, nx, ny, nz )
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        rij(ii,jj,kk) = .5 * eye * ( ky(jj) * wz(ii,jj,kk) + kz(kk) * wy(ii,jj,kk) )
      end do
      end do
      end do
      call rfftwnd_f77_one_complex_to_real(c2r3d,rij,ignore_me)

      helidiss = helidiss + 2. * cmplx( real(rij)*real(tij), aimag(rij)*aimag(tij) )


      helidiss = -2. * helidiss

      if (nfile .eq. 0) then 
        meanhd0 = sum(real(helidiss(1:lx,:,:))) + sum(aimag(helidiss(1:lx,:,:)))
        meanhd0 = meanhd0 * const
        rmshd0 = sum(real(helidiss(1:lx,:,:)**2)) + sum(aimag(helidiss(1:lx,:,:)**2)) 
        rmshd0 = rmshd0 * const
        rmshd0 = sqrt(rmshd0 - meanhd0**2)
        write(*,*) 'estimate meanhd0 = ', meanhd0
        write(*,*) 'estimate rmshd0  = ', rmshd0
      end if

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx

        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1)/2
          ignore_me = real(helidiss(ll,jj,kk))
        else
          ll = ii / 2
          ignore_me = aimag(helidiss(ll,jj,kk))
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

  write(*,*) 'Calculated mean helidiss: ', meanhd
  write(*,*) 'Calculated rms of helidiss: ', rmshd

  open(15, file = 'pdfhd'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat')
    write(15,*) '# Title = " meanhd is', meanhd, 'rmshd is', rmshd, '"'
    write(15,*) '# variables = "(pi_H-<pi_H>)/rms", "PDF"'
    do ii = 1, npnt
      write(15,*) ( (-bound + (ii-.5) * binw) * rmshd0 + meanhd0 - meanhd ) / rmshd, & 
                  pdfhd(ii) * rmshd / rmshd0
    end do
  close(15)

  deallocate( kx, ky, kz, g, ux, uy, uz, wx, wy, wz, rij, tij, helidiss )

  call destroyplan3d

  write(*,*) 'pdfsgshd.x finished'

end program sgshdpdf
