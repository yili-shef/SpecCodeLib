program alignxr
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
  real,    parameter :: binw=1./npnt, binwphi=pi/2/npnt
  real(dp), dimension(npnt) :: palpha1, pbeta1, pgamma1
  real(dp), dimension(npnt) :: palpha2, pbeta2, pgamma2
  real(dp), dimension(npnt) :: palpha3, pbeta3, pgamma3
  real(dp), dimension(npnt) :: palpha4, pbeta4, pgamma4

  real :: delta_c, ignore_me, tmp, q, r
  real(dp) :: alpha, beta, gamm
  integer(8) :: num1, num2, num3, num4
  character(80) :: fnm,str,fpath

  write(*,*) 
  write(*,'(''>>>>>> PDFs of the cos angles between omega and eigenvectors of Rij <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-xr-cnd-gij-qr.x nx filelist ndel'
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

  palpha4=0._dp
  pbeta4=0._dp
  pgamma4=0._dp
  palpha3=0._dp
  pbeta3=0._dp
  pgamma3=0._dp
  palpha2=0._dp
  pbeta2=0._dp
  pgamma2=0._dp
  palpha1=0._dp
  pbeta1=0._dp
  pgamma1=0._dp
  num1 = 0; num2 = 0; num3 = 0; num4 = 0
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

      wx(ii,jj,kk) = eye * ( ky(jj) * r33(ii,jj,kk) - kz(kk) * r22(ii,jj,kk) )
      wy(ii,jj,kk) = eye * ( kz(kk) * r11(ii,jj,kk) - kx(ii) * r33(ii,jj,kk) )
      wz(ii,jj,kk) = eye * ( kx(ii) * r22(ii,jj,kk) - ky(jj) * r11(ii,jj,kk) )

      r12(ii,jj,kk)=.5*eye*(kx(ii)*r22(ii,jj,kk)+ky(jj)*r11(ii,jj,kk))
      r13(ii,jj,kk)=.5*eye*(kx(ii)*r33(ii,jj,kk)+kz(kk)*r11(ii,jj,kk))
      r23(ii,jj,kk)=.5*eye*(ky(jj)*r33(ii,jj,kk)+kz(kk)*r22(ii,jj,kk))
      r11(ii,jj,kk)=eye*kx(ii)*r11(ii,jj,kk)
      r22(ii,jj,kk)=eye*ky(jj)*r22(ii,jj,kk)
      r33(ii,jj,kk)=eye*kz(kk)*r33(ii,jj,kk)
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
 
      if ( mod(ii,2) .eq. 1 ) then
        ll = (ii+1)/2

        delta_c = real( wx(ll, jj, kk) )
        ignore_me = real( wy(ll, jj, kk) )
        tmp = real( wz(ll, jj, kk) )
       
        cc(1,1)=real(r11(ll,jj,kk))
        cc(1,2)=real(r12(ll,jj,kk)) - .5 * tmp
        cc(1,3)=real(r13(ll,jj,kk)) + .5 * ignore_me
        cc(2,1)=cc(1,2) + tmp
        cc(2,2)=real(r22(ll,jj,kk))
        cc(2,3)=real(r23(ll,jj,kk)) - .5 * delta_c
        cc(3,1)=cc(1,3) - ignore_me
        cc(3,2)=cc(2,3) + delta_c
        cc(3,3)=real(r33(ll,jj,kk))

      else
        ll = ii / 2
 
        delta_c = aimag( wx(ll, jj, kk) )
        ignore_me = aimag( wy(ll, jj, kk) )
        tmp = aimag( wz(ll, jj, kk) )
       
        cc(1,1)=aimag(r11(ll,jj,kk))
        cc(1,2)=aimag(r12(ll,jj,kk)) - .5 * tmp
        cc(1,3)=aimag(r13(ll,jj,kk)) + .5 * ignore_me
        cc(2,1)=cc(1,2) + tmp
        cc(2,2)=aimag(r22(ll,jj,kk))
        cc(2,3)=aimag(r23(ll,jj,kk)) - .5 * delta_c
        cc(3,1)=cc(1,3) - ignore_me
        cc(3,2)=cc(2,3) + delta_c
        cc(3,3)=aimag(r33(ll,jj,kk))

      end if
      q = -.5 * sum( cc * transpose(cc) )
      r = -(1./3.) * sum( matmul(cc,cc) * transpose(cc) )

      cc(1,2) = cc(1,2) + .5 * tmp
      cc(1,3) = cc(1,3) - .5 * ignore_me
      cc(2,1) = cc(1,2)
      cc(2,3) = cc(2,3) + .5 * delta_c
      cc(3,1) = cc(1,3)
      cc(3,2) = cc(2,3)
 
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
 
        if ( q .gt. - (27 * r * r / 4)**(1./3.) ) then
          if ( r .ge. 0.) then
            palpha1(ll)=palpha1(ll)+1
            pbeta1(mm)=pbeta1(mm)+1
            pgamma1(nn)=pgamma1(nn)+1
            num1 = num1 + 1
          else
            palpha2(ll)=palpha2(ll)+1
            pbeta2(mm)=pbeta2(mm)+1
            pgamma2(nn)=pgamma2(nn)+1
            num2 = num2 + 1
          end if
        else
          if ( r .lt. 0. ) then
            palpha3(ll)=palpha3(ll)+1
            pbeta3(mm)=pbeta3(mm)+1
            pgamma3(nn)=pgamma3(nn)+1
            num3 = num3 + 1
          else
            palpha4(ll)=palpha4(ll)+1
            pbeta4(mm)=pbeta4(mm)+1
            pgamma4(nn)=pgamma4(nn)+1
            num4 = num4 + 1
          end if
        end if
            
      end if

    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  palpha1 = palpha1 / num1
  pbeta1 = pbeta1 / num1
  pgamma1 = pgamma1 / num1

  palpha2 = palpha2 / num2
  pbeta2 = pbeta2 / num2
  pgamma2 = pgamma2 / num2

  palpha3 = palpha3 / num3
  pbeta3 = pbeta3 / num3
  pgamma3 = pgamma3 / num3

  palpha4 = palpha4 / num4
  pbeta4 = pbeta4 / num4
  pgamma4 = pgamma4 / num4

  write(*,*) 'check palpha1:', sum(palpha1)
  write(*,*) 'check pbeta1:',  sum(pbeta1)
  write(*,*) 'check pgamma1:', sum(pgamma1)

  palpha1=palpha1/binw
  pbeta1=pbeta1/binw
  pgamma1=pgamma1/binw

  palpha2=palpha2/binw
  pbeta2=pbeta2/binw
  pgamma2=pgamma2/binw

  palpha3=palpha3/binw
  pbeta3=pbeta3/binw
  pgamma3=pgamma3/binw

  palpha4=palpha4/binw
  pbeta4=pbeta4/binw
  pgamma4=pgamma4/binw

  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='pxralign-cnd-gij-qr'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (ii-.5)*binw, palpha1(ii),pbeta1(ii),pgamma1(ii), &
       palpha2(ii),pbeta2(ii),pgamma2(ii), &
       palpha3(ii),pbeta3(ii),pgamma3(ii), &
       palpha4(ii),pbeta4(ii),pgamma4(ii)

 
    end do
  close(15)

  deallocate(wx,wy,wz,g,kx,ky,kz)
  deallocate(r11,r12,r13,r22,r23,r33)

  call destroyplan3d

  write(*,*) 'done'

end program alignxr
