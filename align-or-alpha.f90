program alignor
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex, allocatable, dimension(:,:,:) :: wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real, allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer, parameter :: npnt=160
  real,    parameter :: binw=pi/(2.*npnt)
  real(dp), dimension(npnt) :: palpha

  real :: delta_c, ignore_me, tmp
  real(dp) :: alpha, beta, gamm
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of Rij (in helical turb) <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-or-alpha.x nx filelist ndel'
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

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
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
  nfile = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r11
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r22
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r33
    close(10)
    write(*,*) 'after reading data files'
 
    r11=r11*g
    r22=r22*g 
    r33=r33*g
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk)=eye*(ky(jj)*r33(ii,jj,kk)-kz(kk)*r22(ii,jj,kk))
      wy(ii,jj,kk)=eye*(kz(kk)*r11(ii,jj,kk)-kx(ii)*r33(ii,jj,kk))
      wz(ii,jj,kk)=eye*(kx(ii)*r22(ii,jj,kk)-ky(jj)*r11(ii,jj,kk))
      r11(ii,jj,kk)=eye*kx(ii)*wx(ii,jj,kk)
      r22(ii,jj,kk)=eye*ky(jj)*wy(ii,jj,kk)
      r33(ii,jj,kk)=eye*kz(kk)*wz(ii,jj,kk)
      r12(ii,jj,kk)=.5*eye*(kx(ii)*wy(ii,jj,kk)+ky(jj)*wx(ii,jj,kk))
      r13(ii,jj,kk)=.5*eye*(kx(ii)*wz(ii,jj,kk)+kz(kk)*wx(ii,jj,kk))
      r23(ii,jj,kk)=.5*eye*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
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
 
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx
 
      delta_c = real( wx(ii, jj, kk) )
      ignore_me = real( wy(ii, jj, kk) )
      tmp = real( wz(ii, jj, kk) )
 
      cc(1,1)=real(r11(ii,jj,kk))
      cc(1,2)=real(r12(ii,jj,kk))
      cc(1,3)=real(r13(ii,jj,kk))
      cc(2,1)=cc(1,2)
      cc(2,2)=real(r22(ii,jj,kk))
      cc(2,3)=real(r23(ii,jj,kk))
      cc(3,1)=cc(1,3)
      cc(3,2)=cc(2,3)
      cc(3,3)=real(r33(ii,jj,kk))
 
      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
 
      alpha = delta_c*evectors(1,3)+ignore_me*evectors(2,3)+tmp*evectors(3,3)
      beta  = delta_c*evectors(1,2)+ignore_me*evectors(2,2)+tmp*evectors(3,2)
      gamm  = delta_c*evectors(1,1)+ignore_me*evectors(2,1)+tmp*evectors(3,1)
 
      delta_c = sqrt( delta_c * delta_c + ignore_me * ignore_me + tmp * tmp )
 
      alpha = abs(alpha) / delta_c
      alpha = acos(alpha)
      
      ll=floor((alpha)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt ) then
        palpha(ll)=palpha(ll)+1
      end if
 
      delta_c = aimag( wx(ii, jj, kk) )
      ignore_me = aimag( wy(ii, jj, kk) )
      tmp = aimag( wz(ii, jj, kk) )
 
      cc(1,1)=aimag(r11(ii,jj,kk))
      cc(1,2)=aimag(r12(ii,jj,kk))
      cc(1,3)=aimag(r13(ii,jj,kk))
      cc(2,1)=cc(1,2)
      cc(2,2)=aimag(r22(ii,jj,kk))
      cc(2,3)=aimag(r23(ii,jj,kk))
      cc(3,1)=cc(1,3)
      cc(3,2)=cc(2,3)
      cc(3,3)=aimag(r33(ii,jj,kk))
 
      call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
      do ll=1,3
        evectors(:,ll)=evectors(:,ll)/sqrt(sum(evectors(:,ll)**2))
      end do
      
      alpha = delta_c*evectors(1,3)+ignore_me*evectors(2,3)+tmp*evectors(3,3)
      beta  = delta_c*evectors(1,2)+ignore_me*evectors(2,2)+tmp*evectors(3,2)
      gamm  = delta_c*evectors(1,1)+ignore_me*evectors(2,1)+tmp*evectors(3,1)
 
      delta_c = sqrt( delta_c * delta_c + ignore_me * ignore_me + tmp * tmp )
 
      alpha = abs(alpha) / delta_c
      alpha = acos(alpha)
      
      ll=floor((alpha)/binw)+1
      if (ll .ge. 1 .and. ll .le. npnt ) then
        palpha(ll)=palpha(ll)+1
      end if
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  alpha = 1._dp/(nx*ny*nz)/nfile

  write(*,*) 'check palpha:', sum(palpha)*alpha

  palpha=palpha*alpha/binw

  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='poralign-1d-alpha-'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, palpha(ii)
    end do
  close(15)


  deallocate(wx,wy,wz,g,kx,ky,kz)
  deallocate(r11,r12,r13,r22,r23,r33)

  call destroyplan3d

  write(*,*) 'done'

end program alignor
