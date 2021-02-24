program jpdfsno
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: ii, jj, kk, nx, ny, nz, lx, lx1, ly, lz, ll, mm, ndel, nfile
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: sn11, sn12, sn13
  complex(sp), allocatable, dimension(:,:,:) :: sn22, sn23, sn33

  real(sp),    allocatable, dimension(:,:,:) :: g, ooo
  real(sp),    allocatable, dimension(:) :: kx, ky, kz

  character(80) :: fnm, str, str1, fpath
  real(sp) :: ignore_me, tmp, delta_c, nxnorm, nwnorm

  integer,  parameter :: npnto=100, npntx = 100
  real(sp), parameter :: boundo = 15., boundx = 10.
  real(sp), parameter :: bwo = boundo/npnto, bwx = boundx/npntx

  real(dp), dimension(npnto, npntx) :: jox
  real(dp), dimension(npnto) :: pdfo
  real(dp), dimension(npntx) :: pdfx

  real(dp) :: rmso0, rmsx0, rmso, rmsx


  write(*,*) 
  write(*,'(''>>>>>> JPDF of omega and the symmetric part of grad of omega/|omega| <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./jpdf-sno-s(d)p.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: .list file of data files'
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
  ! list file name string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( sn11(lx1,ly,lz), sn12(lx1,ly,lz), sn13(lx1,ly,lz) )
  allocate( sn22(lx1,ly,lz), sn23(lx1,ly,lz), sn33(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz), ooo(nx,ny,nz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  nfile = 0
  jox = 0._dp
  rmso = 0._dp
  rmsx = 0._dp
  pdfo = 0._dp; pdfx = 0._dp
  open(20,file=fnm(1:len_trim(fnm))//'.list')
    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wx
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)wz
      close(10)
      write(*,*) 'after reading data files'

      wx = wx * g; wy = wy * g; wz = wz * g

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ignore_me = eye*(ky(jj)*wz(ii,jj,kk)-kz(kk)*wy(ii,jj,kk))
        tmp       = eye*(kz(kk)*wx(ii,jj,kk)-kx(ii)*wz(ii,jj,kk))
        wz(ii,jj,kk) = eye*(kx(ii)*wy(ii,jj,kk)-ky(jj)*wx(ii,jj,kk))
        wx(ii,jj,kk) = ignore_me
        wy(ii,jj,kk) = tmp
        
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

      delta_c = 1./(nx*ny*nz)
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx
        ignore_me = sqrt( real(wx(ii,jj,kk))**2 +  real(wy(ii,jj,kk))**2 +  real(wz(ii,jj,kk))**2) 
        tmp       = sqrt(aimag(wx(ii,jj,kk))**2 + aimag(wy(ii,jj,kk))**2 + aimag(wz(ii,jj,kk))**2) 

        ignore_me = max(ignore_me, mytiny_sp)
        tmp = max(tmp, mytiny_sp)

        wx(ii,jj,kk) = cmplx( real(wx(ii,jj,kk))/ignore_me, aimag(wx(ii,jj,kk))/tmp ) *delta_c
        wy(ii,jj,kk) = cmplx( real(wy(ii,jj,kk))/ignore_me, aimag(wy(ii,jj,kk))/tmp ) *delta_c
        wz(ii,jj,kk) = cmplx( real(wz(ii,jj,kk))/ignore_me, aimag(wz(ii,jj,kk))/tmp ) *delta_c

        ooo(2*ii-1,jj,kk) = ignore_me
        ooo(2*ii,  jj,kk) = tmp

      end do
      end do
      end do
      call rfftwnd_f77_one_real_to_complex(r2c3d,wx,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,wy,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,wz,ignore_me)

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        sn11(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
        sn22(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
        sn33(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)
        sn12(ii,jj,kk) = eye * ( ky(jj) * wx(ii,jj,kk) + kx(ii) * wy(ii,jj,kk) ) / 2
        sn13(ii,jj,kk) = eye * ( kz(kk) * wx(ii,jj,kk) + kx(ii) * wz(ii,jj,kk) ) / 2
        sn23(ii,jj,kk) = eye * ( kz(kk) * wy(ii,jj,kk) + ky(jj) * wz(ii,jj,kk) ) / 2
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,sn11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,sn12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,sn13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,sn22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,sn23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,sn33,ignore_me)

      if ( nfile .eq. 0 ) then
          rmsx0 = 2.* sum( sn11(1:lx,:,:) * conjg( sn11(1:lx,:,:) ) + &
                       sn22(1:lx,:,:) * conjg( sn22(1:lx,:,:) ) + &
                       sn33(1:lx,:,:) * conjg( sn33(1:lx,:,:) ) + &
                       sn12(1:lx,:,:) * conjg( sn12(1:lx,:,:) ) + &
                       sn13(1:lx,:,:) * conjg( sn13(1:lx,:,:) ) + &
                       sn12(1:lx,:,:) * conjg( sn12(1:lx,:,:) ) + &
                       sn23(1:lx,:,:) * conjg( sn23(1:lx,:,:) ) + &
                       sn13(1:lx,:,:) * conjg( sn13(1:lx,:,:) ) + &
                       sn23(1:lx,:,:) * conjg( sn23(1:lx,:,:) ) )

          rmsx0 = sqrt( rmsx0 / (nx*ny*nz) )

          rmso0 = sum( ooo ) / (nx*ny*nz)
          rmso0 = sqrt(rmso0)
  
          write(*,*) 'Estimated rmso0 rmsx0: ', rmso0, rmsx0
      end if


      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1) / 2
          
          nxnorm = 2.*real( sn11(ll,jj,kk) ) * real( sn11(ll,jj,kk) ) &
                 + real( sn22(ll,jj,kk) ) * real( sn22(ll,jj,kk) ) &
                 + real( sn33(ll,jj,kk) ) * real( sn33(ll,jj,kk) ) &
                 + real( sn12(ll,jj,kk) ) * real( sn12(ll,jj,kk) ) &
                 + real( sn13(ll,jj,kk) ) * real( sn13(ll,jj,kk) ) &
                 + real( sn12(ll,jj,kk) ) * real( sn12(ll,jj,kk) ) &
                 + real( sn23(ll,jj,kk) ) * real( sn23(ll,jj,kk) ) &
                 + real( sn13(ll,jj,kk) ) * real( sn13(ll,jj,kk) ) &
                 + real( sn23(ll,jj,kk) ) * real( sn23(ll,jj,kk) ) 

        else
          ll = ii/2

          nxnorm = 2.*aimag( sn11(ll,jj,kk) ) * aimag( sn11(ll,jj,kk) ) &
                 + aimag( sn22(ll,jj,kk) ) * aimag( sn22(ll,jj,kk) ) &
                 + aimag( sn33(ll,jj,kk) ) * aimag( sn33(ll,jj,kk) ) &
                 + aimag( sn12(ll,jj,kk) ) * aimag( sn12(ll,jj,kk) ) &
                 + aimag( sn13(ll,jj,kk) ) * aimag( sn13(ll,jj,kk) ) &
                 + aimag( sn12(ll,jj,kk) ) * aimag( sn12(ll,jj,kk) ) &
                 + aimag( sn23(ll,jj,kk) ) * aimag( sn23(ll,jj,kk) ) &
                 + aimag( sn13(ll,jj,kk) ) * aimag( sn13(ll,jj,kk) ) &
                 + aimag( sn23(ll,jj,kk) ) * aimag( sn23(ll,jj,kk) ) 

        end if
        nwnorm = ooo(ii,jj,kk) 
        nwnorm = nwnorm * nwnorm

        rmso = rmso + nwnorm
        rmsx = rmsx + nxnorm
        nwnorm = sqrt( nwnorm )/ rmso0 
        nxnorm = sqrt( nxnorm )/ rmsx0 

        ll = floor( nwnorm / bwo) + 1
        mm = floor( nxnorm / bwx) + 1
        if ( ll .ge. 1 .and. ll .le. npnto .and. mm .ge. 1 .and. &
             mm .le. npntx ) then
          jox(ll,mm) = jox(ll,mm) + 1
          pdfo(ll) = pdfo(ll) + 1
          pdfx(mm) = pdfx(mm) + 1
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20) 
  rmso = rmso / (nfile * nx * ny * nz)
  rmso = sqrt(rmso)
  rmsx = rmsx / (nfile * nx * ny * nz)
  rmsx = sqrt(rmsx)

  write(*,*) 'rmso, rmsx: ', rmso, rmsx

  jox = jox / (nfile * nx * ny * nz)
  write(*,*) 'Check normalization of jox: ', sum(jox)
  jox = jox / bwo / bwx

  pdfo = pdfo / (nfile * nx * ny * nz) / bwo
  pdfx = pdfx / (nfile * nx * ny * nz) / bwx

  open(15, file ='jpdfsno-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat' )
    write(15,*) '# Zone T= "(omega,|grad n|)", i=', npnto, ', j=', npntx, ', F=point'
    do ii = 1, npnto
    do jj = 1, npntx
      write(15,'(15E15.4)') (ii-.5)*bwo*rmso0/rmso, (jj-.5)*bwx*rmsx0/rmsx, &
                  jox(ii,jj)*rmso*rmsx/rmso0/rmsx0, pdfo(ii)*pdfx(jj)*rmso*rmsx/rmso0/rmsx0
    end do
    write(15,*)
    end do
  close(15)

  open(15, file ='pdfsno-1d-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat' )
    write(15,*) '# 1D PDFs for omega '
    do ii = 1, npnto
      write(15,*) (ii-.5)*bwo*rmso0/rmso, pdfo(ii)*rmso/rmso0
    end do
    write(15,*)
    write(15,*) '# 1D PDFs for symmetric part of grad n '
    do ii = 1, npntx
      write(15,*) (ii-.5)*bwx*rmsx0/rmsx, pdfx(ii)*rmsx/rmsx0
    end do
  close(15)

  deallocate(kx,ky,kz,wx,wy,wz,g,sn11,sn12,sn13,sn22,sn23,sn33,ooo)

  call destroyplan3d

  write(*,*) 'Finished'


end program jpdfsno 
