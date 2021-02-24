program alignxr
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex(sp), allocatable, dimension(:,:,:) :: xix,xiy,xiz
  complex(sp), allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer, parameter :: npnt=160
  real(sp), parameter :: binw=1._sp/npnt, binwphi=pi/2/npnt
  real(dp), dimension(npnt) :: palpha, pbeta, pgamma
  real(dp), dimension(npnt, npnt) :: p2d

  real(sp) :: delta_c, ignore_me, tmp
  real(dp) :: alpha, beta, gamm
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between xi and eigenvectors of Rij <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-xr.x nx filelist ndel'
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


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(xix(lx1,ly,lz),xiy(lx1,ly,lz),xiz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  palpha=0._dp
  pbeta=0._dp
  pgamma=0._dp
  p2d = 0._dp
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r12
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r13
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r23
    close(10)
    write(*,*) 'after reading data files'
 
    r12=r12*g
    r13=r13*g 
    r23=r23*g
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      r11(ii,jj,kk)=eye*(ky(jj)*r23(ii,jj,kk)-kz(kk)*r13(ii,jj,kk))
      r22(ii,jj,kk)=eye*(kz(kk)*r12(ii,jj,kk)-kx(ii)*r23(ii,jj,kk))
      r33(ii,jj,kk)=eye*(kx(ii)*r13(ii,jj,kk)-ky(jj)*r12(ii,jj,kk))

      xix(ii,jj,kk) = eye * ( ky(jj) * r33(ii,jj,kk) - kz(kk) * r22(ii,jj,kk) )
      xiy(ii,jj,kk) = eye * ( kz(kk) * r11(ii,jj,kk) - kx(ii) * r33(ii,jj,kk) )
      xiz(ii,jj,kk) = eye * ( kx(ii) * r22(ii,jj,kk) - ky(jj) * r11(ii,jj,kk) )

      r12(ii,jj,kk)=.5*eye*(kx(ii)*r22(ii,jj,kk)+ky(jj)*r11(ii,jj,kk))
      r13(ii,jj,kk)=.5*eye*(kx(ii)*r33(ii,jj,kk)+kz(kk)*r11(ii,jj,kk))
      r23(ii,jj,kk)=.5*eye*(ky(jj)*r33(ii,jj,kk)+kz(kk)*r22(ii,jj,kk))
      r11(ii,jj,kk)=eye*kx(ii)*r11(ii,jj,kk)
      r22(ii,jj,kk)=eye*ky(jj)*r22(ii,jj,kk)
      r33(ii,jj,kk)=eye*kz(kk)*r33(ii,jj,kk)
    end do
    end do
    end do

 
    call rfftwnd_f77_one_complex_to_real(c2r3d,xix,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,xiz,ignore_me)
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
 
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii+1)/2
 
        delta_c = real( xix(ll, jj, kk), sp )
        ignore_me = real( xiy(ll, jj, kk), sp )
        tmp = real( xiz(ll, jj, kk), sp)
     
        cc(1,1)=real(r11(ll,jj,kk), sp)
        cc(1,2)=real(r12(ll,jj,kk), sp)
        cc(1,3)=real(r13(ll,jj,kk), sp)
        cc(2,2)=real(r22(ll,jj,kk), sp)
        cc(2,3)=real(r23(ll,jj,kk), sp)
        cc(3,3)=real(r33(ll,jj,kk), sp)
 
      else
        ll = ii/2

        delta_c = aimag( xix(ll, jj, kk) )
        ignore_me = aimag( xiy(ll, jj, kk) )
        tmp = aimag( xiz(ll, jj, kk) )
  
        cc(1,1)=aimag(r11(ll,jj,kk))
        cc(1,2)=aimag(r12(ll,jj,kk))
        cc(1,3)=aimag(r13(ll,jj,kk))
        cc(2,2)=aimag(r22(ll,jj,kk))
        cc(2,3)=aimag(r23(ll,jj,kk))
        cc(3,3)=aimag(r33(ll,jj,kk))

      end if
      cc(2,1)=cc(1,2)
      cc(3,1)=cc(1,3)
      cc(3,2)=cc(2,3)

      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
 
      alpha = delta_c*evectors(1,3)+ignore_me*evectors(2,3)+tmp*evectors(3,3)
      beta  = delta_c*evectors(1,2)+ignore_me*evectors(2,2)+tmp*evectors(3,2)
      gamm  = delta_c*evectors(1,1)+ignore_me*evectors(2,1)+tmp*evectors(3,1)
 
      delta_c = sqrt( delta_c * delta_c + ignore_me * ignore_me + tmp * tmp )
 
      alpha = abs(alpha) / delta_c
      beta = abs(beta) / delta_c
      gamm = abs(gamm) / delta_c
      
      ll=floor((alpha)/binw)+1
      mm=floor((beta)/binw)+1
      nn=floor((gamm)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt .and. &
          mm .ge. 1 .and. mm .le. npnt .and. &
          nn .ge. 1 .and. nn .le. npnt)  then 
 
        palpha(ll)=palpha(ll)+1
        pbeta(mm)=pbeta(mm)+1
        pgamma(nn)=pgamma(nn)+1
      end if

      beta = abs( beta ) / sqrt( beta * beta + gamm * gamm )
      beta = acos( beta )
      
      mm = floor( beta / binwphi ) + 1
      if ( mm .ge. 1 .and. mm .le. npnt .and.  &
           ll .ge. 1 .and. ll .le. npnt ) then
           p2d( ll, mm ) = p2d( ll, mm ) + 1
      end if
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  alpha = 1._dp/(nx*ny*nz)/nfile

  write(*,*) 'check palpha:', sum(palpha)*alpha
  write(*,*) 'check pbeta:',  sum(pbeta)*alpha
  write(*,*) 'check pgamma:', sum(pgamma)*alpha
  write(*,*) 'check p2d:   ', sum(p2d)*alpha

  palpha=palpha*alpha/binw
  pbeta=pbeta*alpha/binw
  pgamma=pgamma*alpha/binw

  p2d = p2d * alpha / ( binw * binwphi )


  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='pxralign-1d'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  open(15, file = 'pxralign-2d'//fnm(1 : len_trim( fnm )))
    write(15,'(''zone t="2d align pdf", i='',I4,'', j='', I4, '', f=point'')') npnt,npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.6)') (ii-.5)*binw, (jj-.5)*binwphi, p2d(ii,jj)
    end do
    end do
  close(15)

  deallocate(xix,xiy,xiz,g,kx,ky,kz)
  deallocate(r11,r12,r13,r22,r23,r33)

  call destroyplan3d

  write(*,*) 'done'

end program alignxr
