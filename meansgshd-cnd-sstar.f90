program msgshdcndsstar
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, ii, jj, kk, ll, mm, ndel, nfile

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(sp), allocatable, dimension(:,:,:) :: r11, r12, r13, r22, r23, r33
  complex(sp), allocatable, dimension(:,:,:) :: t11, t12, t13, t22, t23, t33
  complex(sp), allocatable, dimension(:,:,:) :: pih
  real(dp),    allocatable, dimension(:,:,:) :: sstar
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx, ky, kz

  integer, parameter :: npnts=100
  real(sp), parameter :: bounds = 1.
  real(sp), parameter :: bws = 2.*bounds/npnts

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evrij
  integer :: ierr

  real(dp), dimension(npnts) :: mpihcndsstar, pdfsstar
  real(dp) :: alpha, beta, gamm, ns

  real(sp) :: mhelidiss, delta_c, npih
  real(sp) :: const, ignore_me

  character(80) :: str, str1, dnslist, fpath

  write(*,*) 
  write(*,'(''>>>>>> Mean helicity dissipations conditioned on the geometry of strain rate s*<<<<<<'')')
  write(*,*) 

  ii=iargc()
  if (ii .ne. 4) then
          write(*,*) 
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*)
          write(*,*) ' Usage: ./meansgshd-cnd-sstar.x nx datalist ndel normhelidiss'
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
  allocate(ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),sstar(nx,ny,nz),pih(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=dnslist(1:len_trim(dnslist)))

    mpihcndsstar = 0._dp
    pdfsstar = 0._dp
    nfile = 0
    do while ( .not. eof(20) )

      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

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

      call rsgstauij(ux,ux,t11,g,nx,ny,nz) 
      call rsgstauij(ux,uy,t12,g,nx,ny,nz)
      call rsgstauij(ux,uz,t13,g,nx,ny,nz)
      call rsgstauij(uy,uy,t22,g,nx,ny,nz)
      call rsgstauij(uy,uz,t23,g,nx,ny,nz)
      call rsgstauij(uz,uz,t33,g,nx,ny,nz)

      r11 = -(t11+t22+t33) / 3.

      t11 = t11 + r11
      t22 = t22 + r11
      t33 = t33 + r11
  
      ux=ux*g
      uy=uy*g 
      uz=uz*g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ! rij is now temporarily sij
        r11(ii,jj,kk) = eye * kx(ii) * ux(ii,jj,kk)
        r22(ii,jj,kk) = eye * ky(jj) * uy(ii,jj,kk)
        r33(ii,jj,kk) = eye * kz(kk) * uz(ii,jj,kk)
        r12(ii,jj,kk) = .5 * eye * ( kx(ii) * uy(ii,jj,kk) + ky(jj) * ux(ii,jj,kk) )
        r13(ii,jj,kk) = .5 * eye * ( kx(ii) * uz(ii,jj,kk) + kz(kk) * ux(ii,jj,kk) )
        r23(ii,jj,kk) = .5 * eye * ( ky(jj) * uz(ii,jj,kk) + kz(kk) * uy(ii,jj,kk) )
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii,2) .eq. 1 ) then 
          ll = (ii + 1)/2
          cc(1,1)=real(r11(ll,jj,kk))
          cc(1,2)=real(r12(ll,jj,kk))
          cc(1,3)=real(r13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=real(r22(ll,jj,kk))
          cc(2,3)=real(r23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=real(r33(ll,jj,kk))
        else
          ll = ii / 2
          cc(1,1)=aimag(r11(ll,jj,kk))
          cc(1,2)=aimag(r12(ll,jj,kk))
          cc(1,3)=aimag(r13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=aimag(r22(ll,jj,kk))
          cc(2,3)=aimag(r23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=aimag(r33(ll,jj,kk))
        end if
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)

        alpha = evalues(3); beta = evalues(2); gamm = evalues(1)

        sstar(ii,jj,kk) = - 3 * sqrt(6.) * alpha * beta * gamm &
                     / (alpha*alpha + beta*beta + gamm*gamm)**1.5
    
      end do
      end do
      end do

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ! Vorticity vector
        r11(ii,jj,kk)=eye*(ky(jj)*uz(ii,jj,kk)-kz(kk)*uy(ii,jj,kk))
        r22(ii,jj,kk)=eye*(kz(kk)*ux(ii,jj,kk)-kx(ii)*uz(ii,jj,kk))
        r33(ii,jj,kk)=eye*(kx(ii)*uy(ii,jj,kk)-ky(jj)*ux(ii,jj,kk))
        ! Rij, the symmetric part of the vorticity gradient
        r12(ii,jj,kk)=.5*eye*(kx(ii)*r22(ii,jj,kk)+ky(jj)*r11(ii,jj,kk))
        r13(ii,jj,kk)=.5*eye*(kx(ii)*r33(ii,jj,kk)+kz(kk)*r11(ii,jj,kk))
        r23(ii,jj,kk)=.5*eye*(ky(jj)*r33(ii,jj,kk)+kz(kk)*r22(ii,jj,kk))
        r11(ii,jj,kk)=eye*kx(ii)*r11(ii,jj,kk)
        r22(ii,jj,kk)=eye*ky(jj)*r22(ii,jj,kk)
        r33(ii,jj,kk)=eye*kz(kk)*r33(ii,jj,kk)
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)

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

        ns = sstar(ii,jj,kk)
        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii + 1) / 2
          npih = real( pih(ll,jj,kk) ) / mhelidiss
        else
          ll = ii/2
          npih = aimag( pih(ll,jj,kk) ) / mhelidiss
        end if

        ll = floor( (ns + bounds) / bws ) + 1
        if ( ll .ge. 1 .and. ll .le. npnts ) then
          pdfsstar(ll) = pdfsstar(ll) + 1
          mpihcndsstar(ll) = mpihcndsstar(ll) + npih
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20) 

  mpihcndsstar = mpihcndsstar / (pdfsstar + tiny)

  pdfsstar = pdfsstar / nfile * const

  write(*,*) 'Check normalization of pdfsstar: ', sum(pdfsstar) 

  pdfsstar = pdfsstar / bws

  open(15, file ='mpihcndsstar'//str(1:len_trim(str)) )
    do ii = 1, npnts
      write(15,*) -bounds+(ii-.5)*bws, pdfsstar(ii), mpihcndsstar(ii), pdfsstar(ii)*mpihcndsstar(ii)
    end do
  close(15)

  deallocate(kx,ky,kz,ux,uy,uz,g,sstar,pih)
  deallocate(r11,r12,r13,r22,r23,r33,t11,t12,t13,t22,t23,t33)

  call destroyplan3d

  write(*,*) 'Finished'

end program msgshdcndsstar
