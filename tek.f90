program tek
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: lx1, lx, ly, lz, nx, ny, nz, ii, jj, kk, ll, nfile
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz
  real(sp),    allocatable, dimension(:,:,:) :: k2, atmp
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz
  real(sp),    allocatable, dimension(:)     :: tekx, teky, tekz, tektotal

  character(80) :: str, fnm

  write(*,*) 
  write(*,'(''>>> Energy transfer spectra <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./tek.x nx filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
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

  ny = nx; nz = nx
  lx = nx / 2; lx1 = lx + 1
  ly = ny; lz = nz

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz), k2(lx1,ly,lz) )
  allocate( tekx(lx), teky(lx), tekz(lx), tektotal(lx) )
  allocate( atmp(lx1,ly,lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(20, file=fnm(1:len_trim(fnm))//'.list')

    tektotal = 0.
    tekx = 0.
    teky = 0.
    tekz = 0.
    nfile = 0
    do while ( .not. eof(20))
      read(20,*) str
      write(*,*) str(1:len_trim(str))

      open(15,file='./out/ux'//str(1:len_trim(str)),form='unformatted')
      read(15) ux
      close(15)
      open(15,file='./out/uy'//str(1:len_trim(str)),form='unformatted')
      read(15) uy
      close(15)
      open(15,file='./out/uz'//str(1:len_trim(str)),form='unformatted')
      read(15) uz
      close(15)
  
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk)=eye*(ky(jj)*uz(ii,jj,kk)-kz(kk)*uy(ii,jj,kk))
        wy(ii,jj,kk)=eye*(kz(kk)*ux(ii,jj,kk)-kx(ii)*uz(ii,jj,kk))
        wz(ii,jj,kk)=eye*(kx(ii)*uy(ii,jj,kk)-ky(jj)*ux(ii,jj,kk))
      end do
      end do
      end do

      call convec_dns(ux,uy,uz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
      call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)
  
      ! note the definition of t is the nonlinear transfer for the kinetic energy of
      ! the kth mode |u(vector k)|^2/2.

      atmp = 2.*real(wx*conjg(ux)+wy*conjg(uy)+wz*conjg(uz))
      atmp(1,:,:)=.5*atmp(1,:,:)

      ux = 2. * real( wx * conjg(ux) )
      ux(1,:,:) = .5 * ux(1,:,:)

      uy = 2. * real( wy * conjg(uy) )
      uy(1,:,:) = .5 * uy(1,:,:)

      uz = 2. * real( wz * conjg(uz) )
      uz(1,:,:) = .5 * uz(1,:,:)

      do kk=1,lz
        do jj=1,ly
          do ii=1,lx1
      
            ll=floor(sqrt(k2(ii,jj,kk))+.5)
            if (ll .ge. 1 .and. ll .le. lx) then 
                    tekx(ll) = tekx(ll) + real(ux(ii,jj,kk))
                    teky(ll) = teky(ll) + real(uy(ii,jj,kk))
                    tekz(ll) = tekz(ll) + real(uz(ii,jj,kk))
                    tektotal(ll) = tektotal(ll) + atmp(ii,jj,kk)
            end if
      
          end do
        end do
      end do

      nfile = nfile + 1
 
    end do
  close(20)

  tekx = tekx / nfile
  teky = teky / nfile
  tekz = tekz / nfile
  tektotal  = tektotal  / nfile

  str='tek-'//fnm(1:len_trim(fnm))//'.dat'
  open(15, file = str)
    do ii = 1, lx
      write(15,'(I6, 10E15.4)') ii, tekx(ii), teky(ii), tekz(ii), tektotal(ii)
    end do
  close(15)
  
  deallocate( ux, uy, uz, wx, wy, wz, kx, ky, kz, k2 )
  deallocate( atmp, tekx, teky, tekz, tektotal )

  call destroyplan3d

  write(*,*) 'done.'

end program tek
