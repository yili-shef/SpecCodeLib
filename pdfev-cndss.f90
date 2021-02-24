program pdfevcndss
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,ndel
  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz
  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp), dimension(3,3) :: sij

  integer, parameter :: npnt=100
  real(sp), parameter :: binw=2./npnt
  real(dp), dimension(npnt) :: prstar
  real(dp), dimension(npnt,npnt) :: jprs

  integer :: nfile
  integer(8) :: numpnt
  real(sp) :: atmp, btmp, gtmp, rstar, sstar, ss, mss0, delta_c, ignore_me
  real(dp) :: const, mss
  character(80) :: fnm,str,str1, list2, fpath

  write(*,*) 
  write(*,'(''>>>>>> Joint PDFs of the eigenvalues of given data file<<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfev-cndss.x nx evdatafilelist ndel dnsdatalist'
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
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'
  
  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  ! Start calculation
  open(21, file = list2(1:len_trim(list2))//'.list')
  open(20, file = fnm(1:len_trim(fnm))//'.list')

    mss = 0.d0
    prstar = 0.d0
    jprs = 0.d0

    nfile = 0
    numpnt = 0
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


      s11 = s11 * g
      s22 = s22 * g
      s33 = s33 * g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        s12(ii,jj,kk) = eye * ( kx(ii) * s22(ii,jj,kk) + ky(jj) * s11(ii,jj,kk) )
        s13(ii,jj,kk) = eye * ( kx(ii) * s33(ii,jj,kk) + kz(kk) * s11(ii,jj,kk) )
        s23(ii,jj,kk) = eye * ( ky(jj) * s33(ii,jj,kk) + kz(kk) * s22(ii,jj,kk) )

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

      if (nfile .eq. 0) then
        mss0 = sum( s11(1:lx,:,:) * conjg(s11(1:lx,:,:)) )
        mss0 = mss0 + sum( s22(1:lx,:,:) * conjg(s22(1:lx,:,:)) )
        mss0 = mss0 + sum( s33(1:lx,:,:) * conjg(s33(1:lx,:,:)) )
        mss0 = mss0 + sum( s12(1:lx,:,:) * conjg(s12(1:lx,:,:)) )
        mss0 = mss0 + sum( s13(1:lx,:,:) * conjg(s13(1:lx,:,:)) )
        mss0 = mss0 + sum( s23(1:lx,:,:) * conjg(s23(1:lx,:,:)) )
        mss0 = mss0 * const

        write(*,*) 'Estimate of mean ss = ', mss0
      end if
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
        if ( mod(ii,2) .eq. 1 ) then 
          ll = (ii + 1)/2
          sij(1,1) = real(s11(ll,jj,kk))
          sij(1,2) = real(s12(ll,jj,kk))
          sij(1,3) = real(s13(ll,jj,kk))
          sij(2,2) = real(s22(ll,jj,kk))
          sij(2,3) = real(s23(ll,jj,kk))
          sij(3,3) = real(s33(ll,jj,kk))
        else 
          ll = ii/2
          sij(1,1) = aimag(s11(ll,jj,kk))
          sij(1,2) = aimag(s12(ll,jj,kk))
          sij(1,3) = aimag(s13(ll,jj,kk))
          sij(2,2) = aimag(s22(ll,jj,kk))
          sij(2,3) = aimag(s23(ll,jj,kk))
          sij(3,3) = aimag(s33(ll,jj,kk))
        end if
        sij(2,1) = sij(1,2)
        sij(3,1) = sij(1,3)
        sij(3,2) = sij(2,3)

        ss = sum (sij * sij)
        mss = mss + ss

        sstar = -sqrt(6.) * sum( matmul(sij, sij) * sij) / ss**1.5

        atmp = alpha(ii,jj,kk)
        btmp = beta(ii,jj,kk)
        gtmp = gmma(ii,jj,kk)

        rstar =  atmp * atmp + btmp * btmp + gtmp * gtmp 
        rstar = -3 * sqrt(6.) * atmp * btmp * gtmp / rstar**1.5  ! s^* defined by Lund and Rogers
  
        if ( ss .ge. 2. * mss0 ) then
          numpnt = numpnt + 1
          ll = floor( (rstar + 1) / binw ) + 1
          if (ll .ge. 1 .and. ll .le. npnt) prstar(ll) = prstar(ll) + 1
        end if

        mm = floor( (sstar + 1) / binw ) + 1
        nn = floor( (rstar + 1) / binw ) + 1
        if ( mm .ge. 1 .and. mm .le. npnt .and. nn .ge. 1 .and. nn .le. npnt ) then
          jprs(mm,nn) = jprs(mm,nn) + 1
        end if

      end do
      end do
      end do
 
      nfile = nfile + 1

    end do

  close(20)
  close(21) 

  write(*,*) 'finished loops'

  prstar = prstar / numpnt
  jprs = jprs / nfile * const

  write(*,*) 'check prstar: ', sum(prstar)
  write(*,*) 'check jprs: ', sum(jprs)
  prstar = prstar / binw
  jprs = jprs / binw / binw

  open(15,file='prstar-cndss'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,*) -1.+(ii-.5)*binw, prstar(ii)
    end do
  close(15)

  open(15,file='jprs-cndss'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(sstar,rstar)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) -1 + (ii-.5)*binw, -1+(jj-.5)*binw, jprs(ii,jj)
    end do
    end do
  close(15)
 
  deallocate(alpha,beta,gmma,s11,s12,s13,s22,s23,s33,kx,ky,kz,g)

  write(*,*) 'finished'

end program pdfevcndss
