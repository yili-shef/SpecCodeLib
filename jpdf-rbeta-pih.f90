program jpdfrbetapih
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel,nfilter
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real,    allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,bb,evrij
  integer :: ierr

  real, parameter :: rbound = 2., pbound = 45.
  integer, parameter :: npntr = 40, npntp = 125
  real, parameter :: rbinw = 2.*rbound/npntr, pbinw = 2.*pbound/npntp
  real(dp), dimension(npntp,npntr) :: jpdf

  real(dp) :: const, rmsr, meanhd, rmshd, rmsr0, meanhd0
  real :: delta_c, ignore_me, pih, rbeta
  character(80) :: fnm, str, str1, fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of between the eigenvalues of rij and sgs helicity dissipation <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./jpdf-rbeta-pih.x nx filelist ndel nfilter normdissrate'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        nfileter: 0 for cutoff 1 for Gaussian'
          write(*,*) '        normdissrate: normalization of dissrate, normally the mean.'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter type
  call getarg(4,str)
  read(str, '(I20)') nfilter

  ! normalization factor
  call getarg(5,str)
  read(str, '(F8.6)') meanhd0

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  str='-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  if (nfilter .eq. 0) then
    ! Cutoff filter
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      if ( g(ii,jj,kk) .gt. (lx/ndel)**2 ) then
        g(ii,jj,kk) = 0.
      else
        g(ii,jj,kk) = 1.
      end if
    end do
    end do
    end do
  else if (nfilter .eq. 1) then
    ! Gaussian filter
    g=exp(-g*delta_c**2/24.)
  else
    write(*,*) 'Wrong filter type. Stopped'
    stop
  end if


  jpdf = 0.d0
  rmsr = 0.d0
  meanhd = 0.d0
  rmshd = 0.d0
  nfile = 0
  open(20,file=fnm(1:len_trim(fnm)))

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r12
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r13
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r23
      close(10)
      write(*,*) 'after reading data files'

      call rsgstauij(r12,r12,t11,g,nx,ny,nz) 
      call rsgstauij(r12,r13,t12,g,nx,ny,nz)
      call rsgstauij(r12,r23,t13,g,nx,ny,nz)
      call rsgstauij(r13,r13,t22,g,nx,ny,nz)
      call rsgstauij(r13,r23,t23,g,nx,ny,nz)
      call rsgstauij(r23,r23,t33,g,nx,ny,nz)

      r11 = -(t11+t22+t33) / 3.

      t11 = t11 + r11
      t22 = t22 + r11
      t33 = t33 + r11
  
      r12=r12*g
      r13=r13*g 
      r23=r23*g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ! Vorticity vector
        r11(ii,jj,kk)=eye*(ky(jj)*r23(ii,jj,kk)-kz(kk)*r13(ii,jj,kk))
        r22(ii,jj,kk)=eye*(kz(kk)*r12(ii,jj,kk)-kx(ii)*r23(ii,jj,kk))
        r33(ii,jj,kk)=eye*(kx(ii)*r13(ii,jj,kk)-ky(jj)*r12(ii,jj,kk))
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
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx

        if ( mod(ii,2) .eq. 1 ) then
          ll = (ii+1)/2
          cc(1,1)=real(r11(ll,jj,kk))
          cc(1,2)=real(r12(ll,jj,kk))
          cc(1,3)=real(r13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=real(r22(ll,jj,kk))
          cc(2,3)=real(r23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=real(r33(ll,jj,kk))
 
          bb(1,1)=-real(t11(ll,jj,kk))
          bb(1,2)=-real(t12(ll,jj,kk))
          bb(1,3)=-real(t13(ll,jj,kk))
          bb(2,1)=bb(1,2)
          bb(2,2)=-real(t22(ll,jj,kk))
          bb(2,3)=-real(t23(ll,jj,kk))
          bb(3,1)=bb(1,3)
          bb(3,2)=bb(2,3)
          bb(3,3)=-real(t33(ll,jj,kk))
        else
          ll = ii/2
          cc(1,1)=aimag(r11(ll,jj,kk))
          cc(1,2)=aimag(r12(ll,jj,kk))
          cc(1,3)=aimag(r13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=aimag(r22(ll,jj,kk))
          cc(2,3)=aimag(r23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=aimag(r33(ll,jj,kk))
 
          bb(1,1)=-aimag(t11(ll,jj,kk))
          bb(1,2)=-aimag(t12(ll,jj,kk))
          bb(1,3)=-aimag(t13(ll,jj,kk))
          bb(2,1)=bb(1,2)
          bb(2,2)=-aimag(t22(ll,jj,kk))
          bb(2,3)=-aimag(t23(ll,jj,kk))
          bb(3,1)=bb(1,3)
          bb(3,2)=bb(2,3)
          bb(3,3)=-aimag(t33(ll,jj,kk))
        end if

        ! SGS helicity dissipation
        pih = 2. * sum (bb * cc)

        ! eigenvalues of Rij
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do

  
        if ( mod(ii,2) .eq. 1) then
          ll = (ii+1)/2

          ! t11 will be the SGS helicity dissipation rate
          t11(ll,jj,kk) = cmplx( pih, aimag(t11(ll,jj,kk)) )

          ! r22 is the R_beta eigenvalue
          r22(ll,jj,kk) = cmplx( evalues(2), aimag(r22(ll,jj,kk)) )

          ! t22 is the normalization factor
          pih = evalues(1)* evalues(1) + evalues(2)*evalues(2) + evalues(3)*evalues(3)
          t22(ll,jj,kk) = cmplx( pih, aimag(t22(ll,jj,kk)) )
        else
          ll = ii/2
          
          t11(ll,jj,kk) = cmplx( real(t11(ll,jj,kk)), pih )

          r22(ll,jj,kk) = cmplx( real(r22(ll,jj,kk)), evalues(2) )

          pih = evalues(1)* evalues(1) + evalues(2)*evalues(2) + evalues(3)*evalues(3)
          t22(ll,jj,kk) = cmplx( real(t22(ll,jj,kk)), pih )
        end if

      end do
      end do
      end do

      if (nfile .eq. 0) then

        rmsr0 = sum( real(t22(1:lx,:,:)) ) + sum( aimag(t22(1:lx,:,:)) )
        rmsr0 = sqrt( rmsr0 * const )
      end if

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
        if ( mod(ii, 2) .eq. 1 ) then
          ll = (ii + 1) / 2
          pih = real( t11(ll,jj,kk) )
          rbeta = real( r22(ll,jj,kk) )
          rmsr = rmsr + real( t22(ll,jj,kk) )
          rmshd = rmshd + pih * pih
        else
          ll = ii / 2
          pih = aimag( t11(ll,jj,kk) )
          rbeta = aimag( r22(ll,jj,kk) )
          rmsr = rmsr + aimag( t22(ll,jj,kk) )
        end if
        meanhd = meanhd + pih
        rmshd = rmshd + pih * pih

        rbeta = rbeta / rmsr0
        pih = pih / meanhd0

        nn = floor( (pih + pbound) / pbinw ) + 1
        mm = floor( (rbeta + rbound) / rbinw ) + 1
        if ( mm .ge. 1 .and. mm .le. npntr .and. &
             nn .ge. 1 .and. nn .le. npntp ) then
             jpdf(nn, mm) = jpdf(nn, mm) + 1
        end if

      end do
      end do
      end do

      nfile = nfile + 1
    end do
  close(20)

  meanhd = meanhd * const / nfile
  write(*,*) 'meanhd = ', meanhd

  rmshd = rmshd * const / nfile
  rmshd = sqrt( rmshd - meanhd * meanhd )

  rmsr = rmsr * const / nfile
  rmsr = sqrt(rmsr)

  jpdf = jpdf * const / nfile
  write(*,*) 'check normalization of jpdf: ', sum(jpdf)

  jpdf = jpdf / rbinw / pbinw

  open(15, file = 'jpdf-beta-pih'//str( 1:len_trim(str) ) )
    write(15,*) 'Zone T= "(Pi_H,beta)", i=', npntp, ', j=', npntr, ', F=point'
    do jj=1,npntr
    do ii=1,npntp
      write(15,*) ( -pbound + (ii-.5) * pbinw ), &
                  ( -rbound + (jj-.5) * rbinw ) * rmsr0 / rmsr, jpdf(ii,jj)
    end do
    end do
  close(15)


  call destroyplan3d

  deallocate(t11, t12, t13, t22, t23, t33, r11, r12, r13, r22, r23, r33, g, kx, ky, kz)

  write(*,*) 'finished'
  
end program jpdfrbetapih
