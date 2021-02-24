program pdfevcndoo
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,ndel
  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz
  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, wx, wy, wz

  integer, parameter :: npnt=100
  real(sp), parameter :: binw=2./npnt, bndoo = 10., bwoo = 2.*bndoo/npnt
  real(dp), dimension(npnt,npnt) :: jpro
  real(sp), dimension(3) :: oi

  integer :: nfile
  real(sp) :: atmp, btmp, gtmp, rstar, oo, moo0, delta_c, ignore_me
  real(dp) :: const, moo
  character(80) :: fnm, str, str1, list2, fpath

  write(*,*) 
  write(*,'(''>>>>>> Joint PDFs of the eigenvalues of given data file<<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfev-cndoo.x nx evdatafilelist ndel dnsdatalist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        evdatafilelist: the list of the eigenvalue data files'
          write(*,*) '        ndel:  filter scale = ndel * dx'
          write(*,*) '        dnsdatalist: dnsdatalist.list will be the list of dns data files'
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

  call getarg(3,list2)
  read(list2, '(I20)') ndel
  list2 = adjustl(list2)

  str='-'//list2(1:len_trim(list2))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  call getarg(4,list2)
  
  ny=nx; nz=nx
  const=1.d0/(nx*ny*nz)

  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( alpha(nx,ny,nz), beta(nx,ny,nz), gmma(nx,ny,nz) )
  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( g(lx1,ly,lz) )
  write(*,*) 'arrays allocated'
  
  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  ! Start calculation
  open(21, file = list2(1:len_trim(list2))//'.list')
  open(20, file = fnm(1:len_trim(fnm))//'.list')

    moo = 0.d0
    jpro = 0.d0

    nfile = 0
    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      open(15,file='./out/'//str1(1:len_trim(str1)),form='unformatted')
        read(15) alpha
        read(15) beta
        read(15) gmma
      close(15)

      read(21,*) str1
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


      ux = ux * g
      uy = uy * g
      uz = uz * g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk) = eye * ( ky(jj) * uz(ii,jj,kk) - kz(kk) * uy(ii,jj,kk) )
        wy(ii,jj,kk) = eye * ( kz(kk) * ux(ii,jj,kk) - kx(ii) * uz(ii,jj,kk) )
        wz(ii,jj,kk) = eye * ( kx(ii) * uy(ii,jj,kk) - ky(jj) * ux(ii,jj,kk) )
      end do
      end do
      end do

      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

      if (nfile .eq. 0) then
        moo0 = sum( wx(1:lx,:,:) * conjg(wx(1:lx,:,:)) )
        moo0 = moo0 + sum( wy(1:lx,:,:) * conjg(wy(1:lx,:,:)) )
        moo0 = moo0 + sum( wz(1:lx,:,:) * conjg(wz(1:lx,:,:)) )
        moo0 = moo0 * const

        write(*,*) 'Estimate of mean oo = ', moo0
      end if
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
        if ( mod(ii,2) .eq. 1 ) then 
          ll = (ii + 1)/2
          oi(1) = real(wx(ll,jj,kk))
          oi(2) = real(wy(ll,jj,kk))
          oi(3) = real(wz(ll,jj,kk))
        else 
          ll = ii/2
          oi(1) = aimag(wx(ll,jj,kk))
          oi(2) = aimag(wy(ll,jj,kk))
          oi(3) = aimag(wz(ll,jj,kk))
        end if

        oo = sum (oi * oi)
        moo = moo + oo

        atmp = alpha(ii,jj,kk)
        btmp = beta(ii,jj,kk)
        gtmp = gmma(ii,jj,kk)

        rstar =  atmp * atmp + btmp * btmp + gtmp * gtmp 
        rstar = -3 * sqrt(6.) * atmp * btmp * gtmp / rstar**1.5  ! s^* defined by Lund and Rogers

        oo = oo / moo0
  
        mm = floor( (oo + bndoo) / bwoo ) + 1
        nn = floor( (rstar + 1) / binw ) + 1
        if ( mm .ge. 1 .and. mm .le. npnt .and. nn .ge. 1 .and. nn .le. npnt ) then
          jpro(mm,nn) = jpro(mm,nn) + 1
        end if

      end do
      end do
      end do
 
      nfile = nfile + 1

    end do

  close(20)
  close(21) 

  write(*,*) 'finished loops'

  jpro = jpro / nfile * const

  moo = moo *const /nfile
  write(*,*) 'mean oo is ', moo

  write(*,*) 'check jpro: ', sum(jpro)
  jpro = jpro / binw / bwoo

  open(15,file='jpro-cndoo'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(oo,rstar)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (-bndoo + (ii-.5)*bwoo)*moo0/moo, -1+(jj-.5)*binw, &
      jpro(ii,jj)* moo/moo0
    end do
    end do
  close(15)
 
  deallocate(alpha,beta,gmma,ux,uy,uz,wx,wy,wz,kx,ky,kz,g)

  write(*,*) 'finished'

end program pdfevcndoo
