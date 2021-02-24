program jpdfpihoo
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, ndel, nfile

  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: r11, r12, r13, r22, r23, r33
  complex(sp), allocatable, dimension(:,:,:) :: t11, t12, t13, t22, t23, t33
  complex(sp), allocatable, dimension(:,:,:) :: oo, pih
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  integer, parameter :: npnts=160, npnth = 900
  real(sp), parameter :: bounds = 8., boundh = 45.
  real(sp), parameter :: bws = bounds/npnts, bwh = 2.*boundh/npnth

  real(dp), dimension(npnts, npnth) :: jpihoo
  real(dp) :: rmso

  real(sp) :: mhelidiss, delta_c, ns, npih
  real(sp) :: const, ignore_me, rmso0

  character(80) :: str, str1, dnslist, fpath

  write(*,*) 
  write(*,'(''>>>>>> Joint PDF of helicity dissipations and magnitude of strain rate tensor <<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 4) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./jpdf-pih-oo.x nx datalist ndel normhelidiss'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        datalist: list for dns data files'
          write(*,*) '        ndel: delta_c = ndel * dx'
          write(*,*) '        normhelidiss: normalization factor for helidiss, usually its mean'
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

  ! ndel
  call getarg(3,str1)
  read(str1, '(I20)') ndel
  str1 = adjustl(str1)

  call getarg(4,str)
  read(str, '(F15.6)') mhelidiss

  str = '-'//str1(1:len_trim(str1))//'dx-'//dnslist(1:len_trim(dnslist))//'.dat'
  dnslist = dnslist( 1:len_trim(dnslist) )//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),oo(lx1,ly,lz),pih(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=dnslist(1:len_trim(dnslist)))

    jpihoo = 0._dp
    nfile = 0
    rmso = 0._dp
    do while ( .not. eof(20) )

      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r11
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r22
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r33
      close(10)
      write(*,*) 'after reading data files'

      call rsgstauij(r11,r11,t11,g,nx,ny,nz) 
      call rsgstauij(r11,r22,t12,g,nx,ny,nz)
      call rsgstauij(r11,r33,t13,g,nx,ny,nz)
      call rsgstauij(r22,r22,t22,g,nx,ny,nz)
      call rsgstauij(r22,r33,t23,g,nx,ny,nz)
      call rsgstauij(r33,r33,t33,g,nx,ny,nz)

      wx = -(t11+t22+t33) / 3.

      t11 = t11 + wx
      t22 = t22 + wx
      t33 = t33 + wx
  
      r11=r11*g
      r22=r22*g 
      r33=r33*g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk)=eye*(ky(jj)*r33(ii,jj,kk)-kz(kk)*r22(ii,jj,kk))
        wy(ii,jj,kk)=eye*(kz(kk)*r11(ii,jj,kk)-kx(ii)*r33(ii,jj,kk))
        wz(ii,jj,kk)=eye*(kx(ii)*r22(ii,jj,kk)-ky(jj)*r11(ii,jj,kk))

        r11(ii,jj,kk) = eye * kx(ii) * wx(ii,jj,kk)
        r22(ii,jj,kk) = eye * ky(jj) * wy(ii,jj,kk)
        r33(ii,jj,kk) = eye * kz(kk) * wz(ii,jj,kk)
        r12(ii,jj,kk) = .5 * eye * ( kx(ii) * wy(ii,jj,kk) + ky(jj) * wx(ii,jj,kk) )
        r13(ii,jj,kk) = .5 * eye * ( kx(ii) * wz(ii,jj,kk) + kz(kk) * wx(ii,jj,kk) )
        r23(ii,jj,kk) = .5 * eye * ( ky(jj) * wz(ii,jj,kk) + kz(kk) * wy(ii,jj,kk) )
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
 
      call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)

      ! squared magnitude of vorticity.
      oo = cmplx( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz), &
                  aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) )
  
      if ( nfile .eq. 0 ) then
        rmso0 = sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )
        rmso0 = sqrt( rmso0 / (nx*ny*nz) )
 
        write(*,*) 'Estimated rms of omega: ', rmso0
      end if
      rmso = rmso + sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )

      pih = cmplx( real(t11) * real(r11) + real(t22) * real(r22) & 
                  + real(t33) * real(r33) + 2. * ( real(t12) * real(r12) & 
                  + real(t13) * real(r13) + real(t23) * real(r23) ), &
                    aimag(t11) * aimag(r11) + aimag(t22) * aimag(r22) & 
                  + aimag(t33) * aimag(r33) + 2. * ( aimag(t12) * aimag(r12) & 
                  + aimag(t13) * aimag(r13) + aimag(t23) * aimag(r23) ) )

      pih = - 2. * pih

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1) / 2
          ns = sqrt( real( oo(ll,jj,kk) ) ) / rmso0
          npih = real( pih(ll,jj,kk) ) / mhelidiss
        else
          ll = ii/2
          ns = sqrt( aimag( oo(ll,jj,kk) ) ) / rmso0
          npih = aimag( pih(ll,jj,kk) ) / mhelidiss
        end if

        ll = floor( ns / bws ) + 1
        mm = floor( (npih + boundh) / bwh) + 1
        if ( ll .ge. 1 .and. ll .le. npnts .and. mm .ge. 1 .and. &
             mm .le. npnth ) then
          jpihoo(ll,mm) = jpihoo(ll,mm) + 1
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20) 

  jpihoo = jpihoo / nfile * const

  rmso = rmso / nfile / (nx*ny*nz)
  rmso = sqrt(rmso)

  write(*,*) 'rms of oo: ', rmso

  write(*,*) 'Check normalization of jpihoo: ', sum(jpihoo) 

  jpihoo = jpihoo / bws / bwh

  open(15, file ='jpihoo'//str(1:len_trim(str)) )
    write(15,*) 'Zone T= "(|w|,pih)", i=', npnts, ', j=', npnth, ', F=point'
    do jj = 1, npnth
    do ii = 1, npnts
      write(15,*) ((ii-.5)*bws)*rmso0/rmso, -boundh+(jj-.5)*bwh, jpihoo(ii,jj)
    end do
    end do
  close(15)

  deallocate(kx,ky,kz,wx,wy,wz,g,oo,pih)
  deallocate(r11,r12,r13,r22,r23,r33,t11,t12,t13,t22,t23,t33)

  call destroyplan3d

  write(*,*) 'Finished'

end program jpdfpihoo
