program alignscalargradsij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: ndel,ii,jj,kk,ll, nfile
  integer :: nx,ny,nz,lx,lx1,ly,lz

  complex(sp), allocatable, dimension(:,:,:) :: wx,wy,wz
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real(sp), allocatable, dimension(:,:,:) :: g
  real(sp), allocatable, dimension(:) :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer,  parameter :: npnt=100
  real(sp), parameter :: binw=1./npnt
  real(sp), dimension(npnt) :: palpha, pbeta, pgamma

  character(80) :: fnm,str,fpath, str1
  real(sp) :: delta_c, ignore_me, gphix, gphiy, gphiz
  real(dp) :: alpha, beta, gamm

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of strain rate <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-scalargradsij.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: data file list'
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
  ! file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(30, file = fnm( 1:len_trim(fnm) )//'.list' )

    palpha=0.
    pbeta=0.
    pgamma=0.
    nfile = 0
    do while ( .not. eof(30) )
      read(30, *) str1
      write(*, *) str1(1:len_trim(str1))
 
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
 
      fpath='./out/phi'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s12
      close(10)
      write(*,*) 'after reading data files'
  
      s12 = s12 * g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk)=eye*kx(ii)*s12(ii,jj,kk)
        wy(ii,jj,kk)=eye*ky(jj)*s12(ii,jj,kk)
        wz(ii,jj,kk)=eye*kz(kk)*s12(ii,jj,kk)
      end do
      end do
      end do
 
      s11=s11*g
      s22=s22*g 
      s33=s33*g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        s12(ii,jj,kk)=.5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
        s13(ii,jj,kk)=.5*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
        s23(ii,jj,kk)=.5*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))
        s11(ii,jj,kk)=eye*kx(ii)*s11(ii,jj,kk)
        s22(ii,jj,kk)=eye*ky(jj)*s22(ii,jj,kk)
        s33(ii,jj,kk)=eye*kz(kk)*s33(ii,jj,kk)
      end do
      end do
      end do
  
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
  
        if ( mod(ii,2) .eq. 1 ) then
            ll = (ii + 1) / 2
            cc(1,1)=real(s11(ll,jj,kk))
            cc(1,2)=real(s12(ll,jj,kk))
            cc(1,3)=real(s13(ll,jj,kk))
            cc(2,1)=cc(1,2)
            cc(2,2)=real(s22(ll,jj,kk))
            cc(2,3)=real(s23(ll,jj,kk))
            cc(3,1)=cc(1,3)
            cc(3,2)=cc(2,3)
            cc(3,3)=real(s33(ll,jj,kk))
 
            gphix = real(wx(ll,jj,kk))
            gphiy = real(wy(ll,jj,kk))
            gphiz = real(wz(ll,jj,kk))
 
        else
            ll = ii / 2
            cc(1,1)=aimag(s11(ll,jj,kk))
            cc(1,2)=aimag(s12(ll,jj,kk))
            cc(1,3)=aimag(s13(ll,jj,kk))
            cc(2,1)=cc(1,2)
            cc(2,2)=aimag(s22(ll,jj,kk))
            cc(2,3)=aimag(s23(ll,jj,kk))
            cc(3,1)=cc(1,3)
            cc(3,2)=cc(2,3)
            cc(3,3)=aimag(s33(ll,jj,kk))
 
            gphix = aimag(wx(ll,jj,kk))
            gphiy = aimag(wy(ll,jj,kk))
            gphiz = aimag(wz(ll,jj,kk))
        end if
        ignore_me = sqrt(gphix * gphix + gphiy * gphiy + gphiz * gphiz)
        gphix = gphix / ignore_me
        gphiy = gphiy / ignore_me
        gphiz = gphiz / ignore_me
  
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        do ll=1,3
          evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
        end do
        alpha = abs( gphix*evectors(1,3)+gphiy*evectors(2,3)+gphiz*evectors(3,3))
        beta  = abs( gphix*evectors(1,2)+gphiy*evectors(2,2)+gphiz*evectors(3,2))
        gamm  = abs( gphix*evectors(1,1)+gphiy*evectors(2,1)+gphiz*evectors(3,1))
        
        ll=floor((alpha)/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) palpha(ll)=palpha(ll)+1
        ll=floor((beta)/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll)=pbeta(ll)+1
        ll=floor((gamm)/binw)+1
        if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll)=pgamma(ll)+1
      end do
      end do
      end do
 
      nfile = nfile + 1
    end do
  close(30)

  ignore_me=1./(nx*ny*nz*nfile)
  write(*,*) 'check palpha:', sum(palpha)*ignore_me
  write(*,*) 'check pbeta:',  sum(pbeta)*ignore_me
  write(*,*) 'check pgamma:', sum(pgamma)*ignore_me
  palpha=palpha*ignore_me/binw
  pbeta =pbeta *ignore_me/binw
  pgamma=pgamma*ignore_me/binw
 
  fnm =fnm(1:len_trim(fnm))//'-'//str(1:len_trim(str))//'dx'
  open(15, file='pdf-scalargrad-sij-align-'//fnm(1:len_trim(fnm))//'.dat')
    do ii=1,npnt
      write(15,*) (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  deallocate(wx,wy,wz,kx,ky,kz)
  deallocate(s11,s12,s13,s22,s23,s33)
  deallocate(g)

  call destroyplan3d

  write(*,*) 'done'

end program alignscalargradsij
