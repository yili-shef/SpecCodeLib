program alignor
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx, ny, nz, lx, lx1, ly, lz, nthresh
  integer :: ndel, ii, jj, kk, ll, mm, nn

  integer(8) :: numpnt


  complex, allocatable, dimension(:,:,:) :: wx, wy, wz, dissrate
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real,    allocatable, dimension(:,:,:) :: g, gcut
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer,  parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  integer, parameter :: npnt=160
  real,    parameter :: binw=1./npnt, binwphi=pi/2/npnt
  real(dp), dimension(npnt) :: palpha, pbeta, pgamma
  real(dp), dimension(npnt, npnt) :: p2d

  real :: mdissrate

  real :: delta_c, ignore_me, tmp, kcut2
  real(dp) :: alpha, beta, gamm
  character(80) :: fnm,str,fpath, disslist, thresh

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of Rij (in helical turb) <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-or-cnd.x nx filelist ndel disslist nthresh'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: dns data file list'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*) '        disslist: list of dissipation data files'
          write(*,*) '        nthresh: threshold dissipation = nthresh * mdissrate'
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

  ! file name list for the dissipation data files
  call getarg(4, disslist)
  disslist = adjustl( disslist )

  ! the threshold value parameter
  call getarg(5, thresh)
  read(thresh, '(I20)') nthresh
  thresh = adjustl( thresh )


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx
  kcut2 = ( pi / (delta_c/2) )**2

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(g(lx1,ly,lz), gcut(lx1,ly,lz), dissrate(lx1, ly, lz) )
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  where (g .ge. kcut2)
    gcut =0.
  elsewhere
    gcut = 1.
  endwhere
  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')
  open(21, file = disslist( 1 : len_trim( disslist ) )//'.list' )

    palpha=0._dp
    pbeta=0._dp
    pgamma=0._dp
    p2d = 0._dp

    numpnt = 0

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

      read(21,*) fpath
      write(*,*) fpath( 1 : len_trim(fpath) )

      open(10, file = './out/'//fpath(1 : len_trim(fpath) ), form = 'unformatted')
        read(10) dissrate
      close(10)
      write(*,*) 'after reading data files'

      if ( numpnt .eq. 0 ) then
        mdissrate = sum( real( dissrate(1:lx,:,:) ) ) + sum( aimag( &
                    dissrate(1:lx,:,:) ) ) 
        mdissrate = mdissrate / (nx * ny * nz)
        write(*,*) 'estimated dissipation: ', mdissrate
        mdissrate = 0.26158
        write(*,*) 'mdissrate =', mdissrate
      end if
  
      r11=r11*g*gcut
      r22=r22*g*gcut
      r33=r33*g*gcut
  
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
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx

        if ( mod(ii, 2) .eq. 1 ) then
          ll = (ii+1)/2
          ignore_me = real( dissrate(ll,jj,kk) )
        else
          ll = ii/2
          ignore_me = aimag( dissrate(ll,jj,kk) )
        end if

        if ( ignore_me .ge. nthresh * mdissrate ) then

          numpnt = numpnt + 1
  
          if ( mod(ii, 2) .eq. 1 ) then
            delta_c = real( wx(ll, jj, kk) )
            ignore_me = real( wy(ll, jj, kk) )
            tmp = real( wz(ll, jj, kk) )
      
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
            delta_c = aimag( wx(ll, jj, kk) )
            ignore_me = aimag( wy(ll, jj, kk) )
            tmp = aimag( wz(ll, jj, kk) )
      
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
    

        end if
    

      end do
      end do
      end do
 
    end do

  close(20)
  close(21)

  alpha = 1._dp/numpnt

  write(*,*) 'check palpha:', sum(palpha)*alpha
  write(*,*) 'check pbeta:',  sum(pbeta)*alpha
  write(*,*) 'check pgamma:', sum(pgamma)*alpha
  write(*,*) 'check p3d:   ', sum(p2d)*alpha

  palpha=palpha*alpha/binw
  pbeta=pbeta*alpha/binw
  pgamma=pgamma*alpha/binw

  p2d = p2d * alpha / ( binw * binwphi )


  fnm = '-'//disslist(1:len_trim(disslist))//'-cnd'//thresh(1:len_trim(thresh))//'-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='poralign-1d'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, palpha(ii),pbeta(ii),pgamma(ii)
    end do
  close(15)

  open(15, file = 'poralign-2d'//fnm(1 : len_trim( fnm )))
    write(15,'(''zone t="2d align pdf", i='',I4,'', j='', I4, '', f=point'')') npnt,npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.6)') (ii-.5)*binw, (jj-.5)*binwphi, p2d(ii,jj)
    end do
    end do
  close(15)

  deallocate(wx,wy,wz,g,kx,ky,kz, dissrate)
  deallocate(r11,r12,r13,r22,r23,r33, gcut)

  call destroyplan3d

  write(*,*) 'done'

end program alignor
