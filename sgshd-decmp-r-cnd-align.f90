program sgshddecmpcnd
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalrij, evaltij, fv1, fv2
  real(dp), dimension(3,3) :: cc, bb, evrij, evtauij
  integer :: ierr

!  real(sp), parameter :: mdissrate0 = 0.26158 ! Delta = 16 dx
!  real(sp), parameter :: mdissrate0 = 0.2375 ! Delta =  8 dx

  integer, parameter :: npnt = 50
  real(sp), parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)

  real(dp), dimension(npnt,npnt,npnt) :: jpdf3d, cndhdral, cndhdrbe, cndhdrga

  real(dp) :: const, mdissrate, mpihral, mpihrbe, mpihrga

  real :: delta_c, ignore_me, pih, mdissrate0
  real :: ctheta, cos1, cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  real :: cthtalral, cthtberal, cthtgaral, cthtalrbe, cthtberbe, cthtgarbe
  real :: cthtalrga, cthtberga, cthtgarga, pihral, pihrbe, pihrga
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Mean SGS helicity dissipation conditioned on alignment, decomposed into r (al,be,ga) <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sgshd-decmp-r-cnd-align.x nx filelist ndel normdissrate'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        normdissrate: normalization of dissipation, normally the mean' 
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! normalization factor for dissipation
  call getarg(4, str)
  read(str, '(F10.5)') mdissrate0

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

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  jpdf3d = 0.d0; cndhdral = 0.d0; cndhdrbe = 0.d0; cndhdrga = 0.d0
  nfile=0
  mdissrate = 0._dp
  mpihral = 0._dp; mpihrbe = 0._dp; mpihrga = 0._dp
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
  
        if ( mod(ii,2) .eq. 1) then
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

        pih = 2 * sum(cc * bb)
        mdissrate = mdissrate + pih
  
        call rs(3,3,cc,evalrij,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
        call rs(3,3,bb,evaltij,matz,evtauij,fv1,fv2,ierr)
        do ll=1,3
          evtauij(:,ll)=evtauij(:,ll)/sqrt(sum(evtauij(:,ll)**2))
        end do

        cthtalral = abs( sum( evrij(:,3) * evtauij(:,3) ) )
        cthtberal = abs( sum( evrij(:,3) * evtauij(:,2) ) )
        cthtgaral = abs( sum( evrij(:,3) * evtauij(:,1) ) )

        cthtalrbe = abs( sum( evrij(:,2) * evtauij(:,3) ) )
        cthtberbe = abs( sum( evrij(:,2) * evtauij(:,2) ) )
        cthtgarbe = abs( sum( evrij(:,2) * evtauij(:,1) ) )

        cthtalrga = abs( sum( evrij(:,1) * evtauij(:,3) ) )
        cthtberga = abs( sum( evrij(:,1) * evtauij(:,2) ) )
        cthtgarga = abs( sum( evrij(:,1) * evtauij(:,1) ) )
  
        pihral = evalrij(3) * evaltij(3) * cthtalral**2 + &
                 evalrij(3) * evaltij(2) * cthtberal**2 + &
                 evalrij(3) * evaltij(1) * cthtgaral**2
        pihral = 2 * pihral

        pihrbe = evalrij(2) * evaltij(3) * cthtalrbe**2 + &
                 evalrij(2) * evaltij(2) * cthtberbe**2 + &
                 evalrij(2) * evaltij(1) * cthtgarbe**2
        pihrbe = 2 * pihrbe

        pihrga = evalrij(1) * evaltij(3) * cthtalrga**2 + &
                 evalrij(1) * evaltij(2) * cthtberga**2 + &
                 evalrij(1) * evaltij(1) * cthtgarga**2
        pihrga = 2 * pihrga

        mpihral = mpihral + pihral
        mpihrbe = mpihrbe + pihrbe
        mpihrga = mpihrga + pihrga

        pih = pih / mdissrate0
        pihral = pihral / mdissrate0
        pihrbe = pihrbe / mdissrate0
        pihrga = pihrga / mdissrate0

        ctheta=abs(sum(evrij(:,3)*evtauij(:,3))) 
        
        cos1 = sum(evtauij(:,3)*evrij(:,2))
        cos2 = sum(evtauij(:,3)*evrij(:,1))
        cosphi= cos1/sqrt(cos1*cos1+cos2*cos2)
        phi=acos(abs(cosphi))

        cos3 = sum(evrij(:,1)*evtauij(:,2))
        cos4 = sum(evrij(:,1)*evtauij(:,1))
        coszeta = cos4/sqrt(cos3*cos3+cos4*cos4)
        zeta = acos(abs(coszeta))
  
        ll=floor(ctheta/bwcth)+1
        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1

        if ( ll .ge. 1 .and. ll .le. npnt .and. &
             mm .ge. 1 .and. mm .le. npnt .and. &
             nn .ge. 1 .and. nn .le. npnt ) then

          jpdf3d(ll,mm,nn) = jpdf3d(ll,mm,nn) + 1
          cndhdral(ll,mm,nn) = cndhdral(ll,mm,nn) + pihral
          cndhdrbe(ll,mm,nn) = cndhdrbe(ll,mm,nn) + pihrbe
          cndhdrga(ll,mm,nn) = cndhdrga(ll,mm,nn) + pihrga
 
        end if

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  mdissrate = mdissrate * const / nfile
  mpihral = mpihral * const / nfile
  mpihrbe = mpihrbe * const / nfile
  mpihrga = mpihrga * const / nfile

  cndhdral = cndhdral / (jpdf3d + tiny) 
  cndhdrbe = cndhdrbe / (jpdf3d + tiny)
  cndhdrga = cndhdrga / (jpdf3d + tiny)

  jpdf3d = jpdf3d * const / nfile

  write(*,*) 'check jpdf3d: ', sum(jpdf3d)
  write(*,*) 'check mean dissipation rate: ', mdissrate
  write(*,*) 'check decomp. ral: ', mpihral, 'rbe: ', mpihrbe, 'rga: ', mpihrga

  jpdf3d = jpdf3d / bwcth / bwphi / bwzeta

  open(15,file='cndsgshd-decmp-r-al'//str(1:len_trim(str)))
    write(15,'(''zone t = "decmp cnd hd", i='', i6, '',j='', i6, '', k='', i6, '',f=point'')') npnt, npnt, npnt
    do kk = 1, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, (kk-.5)*bwzeta, &
                            cndhdral(ii,jj,kk), jpdf3d(ii,jj,kk) * cndhdral(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  open(15,file='cndsgshd-decmp-r-be'//str(1:len_trim(str)))
    write(15,'(''zone t = "decmp cnd hd", i='', i6, '',j='', i6, '', k='', i6, '',f=point'')') npnt, npnt, npnt
    do kk = 1, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, (kk-.5)*bwzeta, &
                            cndhdrbe(ii,jj,kk), jpdf3d(ii,jj,kk) * cndhdrbe(ii,jj,kk) 
    end do
    end do
    end do
  close(15)

  open(15,file='cndsgshd-decmp-r-ga'//str(1:len_trim(str)))
    write(15,'(''zone t = "decmp cnd hd", i='', i6, '',j='', i6, '', k='', i6, '',f=point'')') npnt, npnt, npnt
    do kk = 1, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, (kk-.5)*bwzeta, &
                            cndhdrga(ii,jj,kk), jpdf3d(ii,jj,kk) * cndhdrga(ii,jj,kk)
    end do
    end do
    end do
  close(15)


  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz)

  call destroyplan3d

  write(*,*) 'finished'

end program sgshddecmpcnd 
