program msgsedcndssoo
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, ndel, nfile

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  complex(sp), allocatable, dimension(:,:,:) :: t11, t12, t13, t22, t23, t33
  complex(sp), allocatable, dimension(:,:,:) :: oo, ss, pie
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  integer, parameter :: npnts = 80
  real(sp), parameter :: bounds = 4., bws = bounds/npnts

  real(dp), dimension(npnts, npnts) :: jpdfssoo, mpiecndssoo
  real(dp) :: meanoo, meanss

  real(sp) :: menerdiss, delta_c, nss, noo, npie
  real(sp) :: const, ignore_me, meanoo0, meanss0

  character(80) :: str, str1, dnslist, fpath

  write(*,*) 
  write(*,'(''>>>>>> SGS energy dissipation conditioned on vorticity magnitude <<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 4) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./meansgsed-cnd-ssoo.x nx datalist ndel normenerdiss'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        datalist: list for dns data files'
          write(*,*) '        ndel: delta_c = ndel * dx'
          write(*,*) '        normenerdiss: normalization factor for enerdiss, usually its mean'
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
  read(str, '(F15.6)') menerdiss

  str = '-'//str1(1:len_trim(str1))//'dx-'//dnslist(1:len_trim(dnslist))//'.dat'
  dnslist = dnslist( 1:len_trim(dnslist) )//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),oo(lx1,ly,lz),ss(lx1,ly,lz),pie(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=dnslist(1:len_trim(dnslist)))

    mpiecndssoo = 0._dp
    jpdfssoo = 0._dp
    nfile = 0
    meanoo = 0._dp
    meanss = 0._dp
    do while ( .not. eof(20) )

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
      write(*,*) 'after reading data files'

      call rsgstauij(s11,s11,t11,g,nx,ny,nz) 
      call rsgstauij(s11,s22,t12,g,nx,ny,nz)
      call rsgstauij(s11,s33,t13,g,nx,ny,nz)
      call rsgstauij(s22,s22,t22,g,nx,ny,nz)
      call rsgstauij(s22,s33,t23,g,nx,ny,nz)
      call rsgstauij(s33,s33,t33,g,nx,ny,nz)

      s12 = -(t11+t22+t33) / 3.

      t11 = t11 + s12
      t22 = t22 + s12
      t33 = t33 + s12
  
      s11 = s11 * g
      s22 = s22 * g 
      s33 = s33 * g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        s12(ii,jj,kk)=eye*(ky(jj)*s33(ii,jj,kk)-kz(kk)*s22(ii,jj,kk))
        s13(ii,jj,kk)=eye*(kz(kk)*s11(ii,jj,kk)-kx(ii)*s33(ii,jj,kk))
        s23(ii,jj,kk)=eye*(kx(ii)*s22(ii,jj,kk)-ky(jj)*s11(ii,jj,kk))
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 
      ! squared magnitude of vorticity.
      oo = cmplx( real(s12) * real(s12) + real(s13) * real(s13) + real(s23) * real(s23), &
                  aimag(s12) * aimag(s12) + aimag(s13) * aimag(s13) + aimag(s23) * aimag(s23) )
  
      if ( nfile .eq. 0 ) then
        meanoo0 = sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )
        meanoo0 = sqrt( meanoo0 / (nx*ny*nz) )
 
        write(*,*) 'Estimated rms of omega: ', meanoo0
      end if
      meanoo = meanoo + sum( real( oo(1:lx,:,:) ) ) + sum( aimag( oo(1:lx,:,:) ) )

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        s12(ii,jj,kk) = .5 * eye * ( kx(ii) * s22(ii,jj,kk) + ky(jj) * s11(ii,jj,kk) )
        s13(ii,jj,kk) = .5 * eye * ( kx(ii) * s33(ii,jj,kk) + kz(kk) * s11(ii,jj,kk) )
        s23(ii,jj,kk) = .5 * eye * ( ky(jj) * s33(ii,jj,kk) + kz(kk) * s22(ii,jj,kk) )
        s11(ii,jj,kk) = eye * kx(ii) * s11(ii,jj,kk)
        s22(ii,jj,kk) = eye * ky(jj) * s22(ii,jj,kk)
        s33(ii,jj,kk) = eye * kz(kk) * s33(ii,jj,kk)
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)

      ss = cmplx( real(s11) * real(s11) + real(s22) * real(s22) & 
                + real(s33) * real(s33) + 2. * ( real(s12) * real(s12) & 
                + real(s13) * real(s13) + real(s23) * real(s23) ), &
                  aimag(s11) * aimag(s11) + aimag(s22) * aimag(s22) & 
                + aimag(s33) * aimag(s33) + 2. * ( aimag(s12) * aimag(s12) & 
                + aimag(s13) * aimag(s13) + aimag(s23) * aimag(s23) ) )
      ss = 2. * ss

      if ( nfile .eq. 0 ) then
        meanss0 = sum( real( ss(1:lx,:,:) ) ) + sum( aimag( ss(1:lx,:,:) ) )
        meanss0 = sqrt( meanss0 / (nx*ny*nz) )
     
        write(*,*) 'Estimated mean of sijsij: ', meanss0
      end if
      meanss = meanss + sum( real( ss(1:lx,:,:) ) ) + sum( aimag( ss(1:lx,:,:) ) )

      pie = cmplx( real(t11) * real(s11) + real(t22) * real(s22) & 
                 + real(t33) * real(s33) + 2. * ( real(t12) * real(s12) & 
                 + real(t13) * real(s13) + real(t23) * real(s23) ), &
                   aimag(t11) * aimag(s11) + aimag(t22) * aimag(s22) & 
                 + aimag(t33) * aimag(s33) + 2. * ( aimag(t12) * aimag(s12) & 
                 + aimag(t13) * aimag(s13) + aimag(t23) * aimag(s23) ) )

      pie = - pie

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1) / 2
          noo = sqrt( real( oo(ll,jj,kk) ) ) / meanoo0
          nss = sqrt( real( ss(ll,jj,kk) ) ) / meanss0
          npie = real( pie(ll,jj,kk) ) / menerdiss
        else
          ll = ii/2
          noo = sqrt( aimag( oo(ll,jj,kk) ) ) / meanoo0
          nss = sqrt( aimag( ss(ll,jj,kk) ) ) / meanss0
          npie = aimag( pie(ll,jj,kk) ) / menerdiss
        end if

        mm = floor( nss / bws ) + 1
        ll = floor( noo / bws ) + 1
        if ( ll .ge. 1 .and. ll .le. npnts .and. mm .ge. 1 .and. mm .le. npnts  ) then
          jpdfssoo(ll,mm) = jpdfssoo(ll,mm) + 1
          mpiecndssoo(ll,mm) = mpiecndssoo(ll,mm) + npie
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20) 

  mpiecndssoo = mpiecndssoo / (jpdfssoo + mytiny)
  jpdfssoo = jpdfssoo / nfile * const

  meanoo = meanoo / nfile * const
  meanoo = sqrt(meanoo)
  meanss = meanss / nfile * const
  meanss = sqrt(meanss)

  write(*,*) 'rms of oo: ', meanoo
  write(*,*) 'rms of ss: ', meanss
  write(*,*) 'Check normalization of jpdfssoo: ', sum(jpdfssoo) 

  jpdfssoo = jpdfssoo / bws / bws

  open(15, file ='msgsed-cndssoo'//str(1:len_trim(str)) )
    write(15,*) '# oo   ss   jpdfssoo    meanedcndoo'
    do jj = 1, npnts

      do ii = 1, npnts
        write(15,'(20E18.5)') (ii-.5)*bws*meanoo0/meanoo, (jj-.5)*bws*meanss0/meanss, &
                              jpdfssoo(ii,jj)*meanoo*meanss/meanoo0/meanss0, &
                              mpiecndssoo(ii,jj)
      end do
      write(15,*)

    end do
  close(15)

  deallocate(kx,ky,kz,g,oo,ss,pie)
  deallocate(s11,s12,s13,s22,s23,s33,t11,t12,t13,t22,t23,t33)

  call destroyplan3d

  write(*,*) 'Finished'

end program msgsedcndssoo
