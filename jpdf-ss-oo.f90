program jpssoo
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex(sp), allocatable, dimension(:,:,:) :: s11,s22,s33,s12,s13,s23,wx,wy,wz
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt = 80
  real(sp), parameter :: bndss = 10., bndoo = 10.
  real(sp), parameter :: binwss = bndss/npnt, binwoo = bndoo/npnt

  real(dp), dimension(npnt, npnt) :: jpdfssoo

  real(dp) :: const, meanoo, meanss
  real(sp) :: delta_c, ignore_me, meanoo0, ss, oo
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>>>>> JPDF of gg and aa <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./jpdf-ss-oo.x nx filelist ndel'
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

  str='-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx
  const = 1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( kx(lx1), ky(ly), kz(lz), g(lx1,ly,lz) )
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate(  wx(lx1,ly,lz),  wy(lx1,ly,lz),  wz(lx1,ly,lz) ) 
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  jpdfssoo = 0._dp
  meanoo = 0._dp
  meanss = 0._dp
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s11
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s22
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)s33
    close(10)
    write(*,*) 'after reading data files'

    s11 = s11 * g
    s22 = s22 * g
    s33 = s33 * g
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1

      ! vorticity 
      wx(ii,jj,kk) = eye * ( ky(jj) * s33(ii,jj,kk) - kz(kk) * s22(ii,jj,kk) )
      wy(ii,jj,kk) = eye * ( kz(kk) * s11(ii,jj,kk) - kx(ii) * s33(ii,jj,kk) )
      wz(ii,jj,kk) = eye * ( kx(ii) * s22(ii,jj,kk) - ky(jj) * s11(ii,jj,kk) )

      ! strain rate tensor
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
    call rfftwnd_f77_one_complex_to_real(c2r3d, wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d, wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d, wz,ignore_me)


    if ( nfile .eq. 0 ) then
        meanoo0 = sum( wx(1:lx,:,:) * conjg(wx(1:lx,:,:)) + &
                       wy(1:lx,:,:) * conjg(wy(1:lx,:,:)) + &
                       wz(1:lx,:,:) * conjg(wz(1:lx,:,:)) )
        meanoo0 = meanoo0 / (nx * ny * nz)
                  
        write(*,*) 'Estimated mean of oo: ', meanoo0
    end if

    
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii + 1 )/2

        ss = real( s11(ll,jj,kk) ) * real( s11(ll,jj,kk) ) &
           + real( s22(ll,jj,kk) ) * real( s22(ll,jj,kk) ) &
           + real( s33(ll,jj,kk) ) * real( s33(ll,jj,kk) ) &
           + 2. * ( real( s12(ll,jj,kk) ) * real( s12(ll,jj,kk) ) &
                  + real( s13(ll,jj,kk) ) * real( s13(ll,jj,kk) ) &
                  + real( s23(ll,jj,kk) ) * real( s23(ll,jj,kk) ) )

        oo = real( wx(ll,jj,kk) ) * real( wx(ll,jj,kk) )  &
           + real( wy(ll,jj,kk) ) * real( wy(ll,jj,kk) )  &
           + real( wz(ll,jj,kk) ) * real( wz(ll,jj,kk) )
      else
        ll = ii / 2

        ss = real( s11(ll,jj,kk) ) * real( s11(ll,jj,kk) ) &
           + real( s22(ll,jj,kk) ) * real( s22(ll,jj,kk) ) &
           + real( s33(ll,jj,kk) ) * real( s33(ll,jj,kk) ) &
           + 2. * ( real( s12(ll,jj,kk) ) * real( s12(ll,jj,kk) ) &
                  + real( s13(ll,jj,kk) ) * real( s13(ll,jj,kk) ) &
                  + real( s23(ll,jj,kk) ) * real( s23(ll,jj,kk) ) )

        oo = aimag( wx(ll,jj,kk) ) * aimag( wx(ll,jj,kk) )  &
           + aimag( wy(ll,jj,kk) ) * aimag( wy(ll,jj,kk) )  &
           + aimag( wz(ll,jj,kk) ) * aimag( wz(ll,jj,kk) )
      end if
      meanoo = meanoo + oo
      meanss = meanss + ss

      oo = oo / meanoo0
      ss = ss / meanoo0

      ll = floor( ss / binwss ) + 1
      mm = floor( oo / binwoo ) + 1

      if ( ll .ge. 1 .and. ll .le. npnt .and. mm .ge. 1 .and. mm .le. npnt) then
        jpdfssoo(mm,ll) = jpdfssoo(mm,ll) + 1
      end if

    end do
    end do
    end do
    nfile = nfile + 1
  end do
  close(20)

  meanoo = meanoo / nfile * const
  write(*,*) 'meanoo is: ', meanoo

  meanss = meanss / nfile * const
  write(*,*) 'meanss is: ', meanss

  jpdfssoo = jpdfssoo / nfile * const
  write(*,*) 'Check normalizationof jpdfssoo:', sum(jpdfssoo)

  jpdfssoo = jpdfssoo / binwss / binwoo

  open(15,file='jpdf-ss-oo'//str(1:len_trim(str)))
      write(15,*) '# Zone T= "(r,q)", i=', npnt, ', j=', npnt, ', F=point'
      do ii=1,npnt
      do jj=1,npnt
          write(15,'(15E15.5)') ((ii-.5)*binwoo)*meanoo0/meanoo, & 
                                ((jj-.5)*binwss)*meanoo0/meanoo, &
              jpdfssoo(ii,jj)*meanoo*meanoo/meanoo0/meanoo0
      end do
      write(15,*)
      end do
  close(15)

  deallocate(kx,ky,kz,g,s11,s12,s13,s22,s23,s33,wx,wy,wz)

  call destroyplan3d
  write(*,*) 'Finished'

end program jpssoo
