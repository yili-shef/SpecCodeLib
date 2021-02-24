program jpdfro
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: ii, jj, kk, nx, ny, nz, lx, lx1, ly, lz, ll, mm, nn, ndel, nfile
  complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
  complex(sp), allocatable, dimension(:,:,:) :: r11, r12, r13, r22, r23, r33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  character(80) :: fnm, str, str1, fpath
  real(sp) :: ignore_me, tmp, delta_c, nrnorm, nwnorm

  integer,  parameter :: npnto=100, npntr = 200
  real(sp), parameter :: boundo = 15., boundr = 30.
  real(sp), parameter :: bwo = boundo/npnto, bwr = boundr/npntr

  real(dp), dimension(npnto, npntr) :: jor

  real(dp) :: rmso0, rmsr0, rmso, rmsr


  write(*,*) 
  write(*,'(''>>>>>> JPDF of omega and |rij| <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./jpdf-ro-s(d)p.x nx filelist ndel'
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
  allocate( r11(lx1,ly,lz), r12(lx1,ly,lz), r13(lx1,ly,lz) )
  allocate( r22(lx1,ly,lz), r23(lx1,ly,lz), r33(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  nfile = 0
  jor = 0._dp
  rmso = 0._dp
  rmsr = 0._dp
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
        r12(ii,jj,kk) = .5 * eye*(ky(jj)*wx(ii,jj,kk)+kx(kk)*wy(ii,jj,kk))
        r13(ii,jj,kk) = .5 * eye*(kz(kk)*wx(ii,jj,kk)+kx(ii)*wz(ii,jj,kk))
        r23(ii,jj,kk) = .5 * eye*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
        r11(ii,jj,kk) = eye*kx(ii)*wx(ii,jj,kk)
        r22(ii,jj,kk) = eye*ky(jj)*wy(ii,jj,kk)
        r33(ii,jj,kk) = eye*kz(kk)*wz(ii,jj,kk)
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

      if ( nfile .eq. 0 ) then
          rmso0 = sum( wx(1:lx,:,:)*conjg( wx(1:lx,:,:) ) + wy(1:lx,:,:)*conjg(wy(1:lx,:,:)) &
                     + wz(1:lx,:,:)*conjg( wz(1:lx,:,:) ) ) / (nx*ny*nz)
          rmso0 = sqrt(rmso0)

          rmsr0 = sum( r11(1:lx,:,:) * conjg( r11(1:lx,:,:) ) + &
                       r22(1:lx,:,:) * conjg( r22(1:lx,:,:) ) + &
                       r33(1:lx,:,:) * conjg( r33(1:lx,:,:) ) + &
                   2.*(r12(1:lx,:,:) * conjg( r12(1:lx,:,:) ) + &
                       r13(1:lx,:,:) * conjg( r13(1:lx,:,:) ) + &
                       r23(1:lx,:,:) * conjg( r23(1:lx,:,:) ) ) ) 
  
          rmsr0 = sqrt( rmsr0 / (nx*ny*nz) )
          write(*,*) 'Estimated rmso0 rmsr0: ', rmso0, rmsr0
      end if


      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1) / 2
          nwnorm = real( wx(ll,jj,kk) ) * real( wx(ll,jj,kk) ) &
                 + real( wy(ll,jj,kk) ) * real( wy(ll,jj,kk) ) &
                 + real( wz(ll,jj,kk) ) * real( wz(ll,jj,kk) )

          nrnorm = real( r11(ll,jj,kk) ) * real( r11(ll,jj,kk) ) &
                 + real( r22(ll,jj,kk) ) * real( r22(ll,jj,kk) ) &
                 + real( r33(ll,jj,kk) ) * real( r33(ll,jj,kk) ) &
                 + 2. * ( real(r12(ll,jj,kk)) ** 2  & 
                 + real(r13(ll,jj,kk)) **2 + real(r23(ll,jj,kk)) **2 ) 
        else
          ll = ii/2
          nwnorm = aimag( wx(ll,jj,kk) ) * aimag( wx(ll,jj,kk) ) &
                 + aimag( wy(ll,jj,kk) ) * aimag( wy(ll,jj,kk) ) &
                 + aimag( wz(ll,jj,kk) ) * aimag( wz(ll,jj,kk) )

          nrnorm = aimag( r11(ll,jj,kk) ) * aimag( r11(ll,jj,kk) ) &
                 + aimag( r22(ll,jj,kk) ) * aimag( r22(ll,jj,kk) ) &
                 + aimag( r33(ll,jj,kk) ) * aimag( r33(ll,jj,kk) ) &
                 + 2. * ( aimag(r12(ll,jj,kk)) ** 2  & 
                 + aimag(r13(ll,jj,kk)) ** 2 + aimag(r23(ll,jj,kk)) ** 2 ) 
        end if

        rmso = rmso + nwnorm
        rmsr = rmsr + nrnorm
        nwnorm = sqrt( nwnorm )/ rmso0 
        nrnorm = sqrt( nrnorm )/ rmsr0 

        ll = floor( nwnorm / bwo) + 1
        mm = floor( nrnorm / bwr) + 1
        if ( ll .ge. 1 .and. ll .le. npnto .and. mm .ge. 1 .and. &
             mm .le. npntr ) then
          jor(ll,mm) = jor(ll,mm) + 1
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20) 
  rmso = rmso / (nfile * nx * ny * nz)
  rmso = sqrt(rmso)
  rmsr = rmsr / (nfile * nx * ny * nz)
  rmsr = sqrt(rmsr)

  write(*,*) 'rmso, rmsr: ', rmso, rmsr

  jor = jor / (nfile * nx * ny * nz)
  write(*,*) 'Check normalization of jor: ', sum(jor)
  jor = jor / bwo / bwr

  open(15, file ='jpdfro-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat' )
    write(15,*) '# Zone T= "(omega,|rij|)", i=', npnto, ', j=', npntr, ', F=point'
    do ii = 1, npnto
    do jj = 1, npntr
      write(15,*) (ii-.5)*bwo*rmso0/rmso, (jj-.5)*bwr*rmsr0/rmsr, &
                  jor(ii,jj)*rmso*rmsr/rmso0/rmsr0
    end do
    write(15,*)
    end do
  close(15)

  deallocate(kx,ky,kz,wx,wy,wz,g,r11,r22,r33,r12,r13,r23)

  call destroyplan3d

  write(*,*) 'Finished'


end program jpdfro 
