program sgshdcnd
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,bb,evrij,evtauij
  integer :: ierr

!  real(sp), parameter :: mdissrate0 = 0.3

  integer, parameter :: npnt = 40
  real(sp), parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)

  real(dp), dimension(npnt) :: pcth,pphi,pzeta, cndcth,cndphi,cndzeta
  real(dp), dimension(npnt) :: pcthbeta, pcthgamma, cndcthbeta,cndcthgamma

  real(dp) :: const, mdissrate

  real :: delta_c, ignore_me, pih, mdissrate0
  real :: ctheta, cos1, cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Mean SGS helicity dissipation conditioned on alignment <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgshd-cnd-align.x nx filelist ndel mdissrate0'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        mdissrate0: normalization factor of the dissipation'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! mdissrate0
  call getarg(4,str)
  read(str, '(F15.6)') mdissrate0

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

  pcth = 0.d0; pzeta = 0.d0; pphi = 0.d0
  pcthbeta = 0.d0; pcthgamma = 0.d0
  cndcth = 0.d0; cndphi = 0.d0; cndzeta = 0.d0
  cndcthbeta = 0.d0; cndcthgamma = 0.d0
  nfile=0
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
  
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
        call rs(3,3,bb,evalues,matz,evtauij,fv1,fv2,ierr)
        do ll=1,3
          evtauij(:,ll)=evtauij(:,ll)/sqrt(sum(evtauij(:,ll)**2))
        end do
  
        ctheta=abs(sum(evrij(:,3)*evtauij(:,3))) 
        
        cos1 = sum(evtauij(:,3)*evrij(:,2))
        cos2 = sum(evtauij(:,3)*evrij(:,1))
        cosphi= cos1/sqrt(cos1*cos1+cos2*cos2)
        phi=acos(abs(cosphi))

        cos3 = sum(evrij(:,1)*evtauij(:,2))
        cos4 = sum(evrij(:,1)*evtauij(:,1))
        coszeta = cos4/sqrt(cos3*cos3+cos4*cos4)
        zeta = acos(abs(coszeta))
  
        pih = pih / mdissrate0
  
        ll=floor(ctheta/bwcth)+1
        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1

        pcth(ll) = pcth(ll) + 1
        cndcth(ll) = cndcth(ll) + pih

        pphi(mm) = pphi(mm) + 1
        cndphi(mm) = cndphi(mm) + pih

        pzeta(nn) = pzeta(nn) + 1
        cndzeta(nn) = cndzeta(nn) + pih

        ctheta=abs(sum(evrij(:,2)*evtauij(:,2))) 
        ll=floor(ctheta/bwcth)+1
        pcthbeta(ll) = pcthbeta(ll) + 1
        cndcthbeta(ll) = cndcthbeta(ll) + pih

        ctheta=abs(sum(evrij(:,1)*evtauij(:,1))) 
        mm=floor(ctheta/bwcth)+1
        pcthgamma(mm) = pcthgamma(mm) + 1
        cndcthgamma(mm) = cndcthgamma(mm) + pih

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  cndcth = cndcth / (pcth + tiny) 
  cndphi = cndphi / (pphi + tiny)
  cndzeta = cndzeta / (pzeta + tiny)
  cndcthbeta = cndcthbeta / (pcthbeta + tiny)
  cndcthgamma = cndcthgamma / (pcthgamma + tiny)

  pcth =pcth *const/nfile
  pphi =pphi *const/nfile
  pzeta = pzeta*const/nfile
  pcthbeta = pcthbeta * const / nfile
  pcthgamma = pcthgamma * const / nfile

  write(*,*) 'check pcth: ', sum(pcth)
  write(*,*) 'check pphi:', sum(pphi)
  write(*,*) 'check pzeta:', sum(pzeta)
  write(*,*) 'check pcthbeta:    ', sum( pcthbeta)
  write(*,*) 'check pcthgamma:    ', sum( pcthgamma)
  write(*,*) 'input normalization mdissrate0: ', mdissrate0

  pcth = pcth/bwcth
  pzeta=pzeta/bwzeta
  pphi=pphi/bwphi
  pcthbeta = pcthbeta / bwcth 
  pcthgamma = pcthgamma / bwcth

  open(15,file='cndsgshd-align-cth'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, pcth(ii),cndcth(ii),pcthbeta(ii), &
                            cndcthbeta(ii),pcthgamma(ii),cndcthgamma(ii)
    end do
  close(15)

  open(15,file='cndsgshd-align-phi'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5)*bwphi, pphi(ii), cndphi(ii)
    end do
  close(15)

  open(15,file='cndsgshd-align-zeta'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5)*bwzeta, pzeta(ii), cndzeta(ii)
    end do
  close(15)

  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz)

  call destroyplan3d

  write(*,*) 'finished'

end program sgshdcnd 
