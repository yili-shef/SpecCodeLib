program align3dst
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  integer :: dnsnx, dnsny, dnsnz, dnslx, dnslx1, dnsly, dnslz
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  complex, allocatable, dimension(:,:,:) :: wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: dnsux,dnsuy,dnsuz
  real,    allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: bb, cc,evsij,evtauij
  integer :: ierr

  integer, parameter :: npnt=30
  real, parameter :: bwcth=1./npnt, bwphi=pi/2./npnt, bwzeta=pi/2./npnt
  real, dimension(npnt,npnt) :: pcthphi,pcthzeta,pphizeta
  real, dimension(npnt,npnt,npnt) :: p3d
  real, dimension(npnt) :: ptha, pthb, pthg

  real :: const, delta_c, ignore_me
  real :: ctheta,cos1,cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of the angles between eigenvectors of sij and sgs stress minor Leonard stress <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3dst-sub-horiuti.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') dnsnx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  str='-'//str(1:len_trim(str))//'dx-st-sub-horiuti-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  dnsny = dnsnx; dnsnz = dnsnx
  dnslx = dnsnx/2; dnsly = dnsny; dnslz = dnsnz
  dnslx1 = dnslx + 1
  
  nx = 2 * dnsnx / ndel ! coarse-graining

  ny = nx; nz = nx
  lx=nx/2; ly=ny; lz=nz; lx1=lx+1

  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/dnsnx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( dnsux(dnslx1, dnsly, dnslz) )
  allocate( dnsuy(dnslx1, dnsly, dnslz) )
  allocate( dnsuz(dnslx1, dnsly, dnslz) )

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  pcthphi = 0.; pcthzeta = 0.; pphizeta = 0.
  ptha = 0. ; pthb = 0. ; pthg = 0. ;
  p3d = 0.
  nfile=0
  open(20,file=fnm(1:len_trim(fnm)))

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)dnsux
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)dnsuy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)dnsuz
      close(10)
      write(*,*) 'after reading data files'

      ! padding
      call unpadd ( dnsux, s11, lx1, ly, lz, dnslx1, dnsly, dnslz )
      call unpadd ( dnsuy, s22, lx1, ly, lz, dnslx1, dnsly, dnslz )
      call unpadd ( dnsuz, s33, lx1, ly, lz, dnslx1, dnsly, dnslz )

      call rsgstauij(s11,s11,t11,g,nx,ny,nz) 
      call rsgstauij(s11,s22,t12,g,nx,ny,nz)
      call rsgstauij(s11,s33,t13,g,nx,ny,nz)
      call rsgstauij(s22,s22,t22,g,nx,ny,nz)
      call rsgstauij(s22,s33,t23,g,nx,ny,nz)
      call rsgstauij(s33,s33,t33,g,nx,ny,nz)
  
      s11=s11*g
      s22=s22*g 
      s33=s33*g

      s12 = -(t11 + t22 + t33)/3.
      t11 = t11 + s12
      t22 = t22 + s12
      t33 = t33 + s12

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk) = eye * ( ky(jj) * s33(ii,jj,kk) - kz(kk) * s22(ii,jj,kk) )
        wy(ii,jj,kk) = eye * ( kz(kk) * s11(ii,jj,kk) - kx(ii) * s33(ii,jj,kk) )
        wz(ii,jj,kk) = eye * ( kx(ii) * s22(ii,jj,kk) - ky(jj) * s11(ii,jj,kk) )
        s12(ii,jj,kk)=.5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
        s13(ii,jj,kk)=.5*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
        s23(ii,jj,kk)=.5*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))
        s11(ii,jj,kk)=eye*kx(ii)*s11(ii,jj,kk)
        s22(ii,jj,kk)=eye*ky(jj)*s22(ii,jj,kk)
        s33(ii,jj,kk)=eye*kz(kk)*s33(ii,jj,kk)
      end do
      end do
      end do
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
  
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

          bb(1,1) = 0.; bb(2,2) = 0.; bb(3,3) = 0.

          bb(1,2) = -.5 * real( wz(ll,jj,kk) )
          bb(1,3) =  .5 * real( wy(ll,jj,kk) )
          bb(2,3) = -.5 * real( wx(ll,jj,kk) )
          bb(2,1) = - bb(1,2)
          bb(3,1) = - bb(1,3)
          bb(3,2) = - bb(2,3)

          bb = - delta_c**2 / 12 * ( matmul(cc,bb) + transpose( matmul(cc,bb) ) )
         
          bb(1,1)=-real(t11(ll,jj,kk)) + bb(1,1)
          bb(1,2)=-real(t12(ll,jj,kk)) + bb(1,2)
          bb(1,3)=-real(t13(ll,jj,kk)) + bb(1,3)
          bb(2,1)=bb(1,2)
          bb(2,2)=-real(t22(ll,jj,kk)) + bb(2,2)
          bb(2,3)=-real(t23(ll,jj,kk)) + bb(2,3) 
          bb(3,1)=bb(1,3)
          bb(3,2)=bb(2,3)
          bb(3,3)=-real(t33(ll,jj,kk)) + bb(3,3)

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
         
          bb(1,1) = 0.; bb(2,2) = 0.; bb(3,3) = 0.

          bb(1,2) = -.5 * aimag( wz(ll,jj,kk) )
          bb(1,3) =  .5 * aimag( wy(ll,jj,kk) )
          bb(2,3) = -.5 * aimag( wx(ll,jj,kk) )
          bb(2,1) = - bb(1,2)
          bb(3,1) = - bb(1,3)
          bb(3,2) = - bb(2,3)

          bb = - delta_c**2 / 12 * ( matmul(cc,bb) + transpose( matmul(cc,bb) ) )
         
          bb(1,1)=-aimag(t11(ll,jj,kk)) + bb(1,1)
          bb(1,2)=-aimag(t12(ll,jj,kk)) + bb(1,2)
          bb(1,3)=-aimag(t13(ll,jj,kk)) + bb(1,3)
          bb(2,1)=bb(1,2)                        
          bb(2,2)=-aimag(t22(ll,jj,kk)) + bb(2,2)
          bb(2,3)=-aimag(t23(ll,jj,kk)) + bb(2,3)
          bb(3,1)=bb(1,3)                        
          bb(3,2)=bb(2,3)                        
          bb(3,3)=-aimag(t33(ll,jj,kk)) + bb(3,3)
        end if

  
        call rs(3,3,cc,evalues,matz,evsij,fv1,fv2,ierr)
        do ll=1,3
          evsij(:,ll)=evsij(:,ll)/sqrt(sum(evsij(:,ll)**2))
        end do
  
        call rs(3,3,bb,evalues,matz,evtauij,fv1,fv2,ierr)
        do ll=1,3
          evtauij(:,ll)=evtauij(:,ll)/sqrt(sum(evtauij(:,ll)**2))
        end do
  
        ctheta=abs(sum(evsij(:,3)*evtauij(:,3))) 
        cos1 = sum(evtauij(:,3)*evsij(:,2))
        cos2 = sum(evtauij(:,3)*evsij(:,1))
        cosphi= cos1/sqrt(cos1*cos1+cos2*cos2)
        phi=acos(abs(cosphi))
        cos3 = sum(evsij(:,1)*evtauij(:,2))
        cos4 = sum(evsij(:,1)*evtauij(:,1))
        coszeta = cos4/sqrt(cos3*cos3+cos4*cos4)
        zeta = acos(abs(coszeta))
  
  
        ll=floor(ctheta/bwcth)+1
        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1
  
        pcthphi(ll,mm) =pcthphi(ll,mm) +1
        pcthzeta(ll,nn)=pcthzeta(ll,nn)+1
        pphizeta(mm,nn)=pphizeta(mm,nn)+1
  
        p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1

        ptha(ll) = ptha(ll) + 1

        ctheta = abs( sum( evsij(:,2) * evtauij(:,2) ) )
        ll = floor( ctheta / bwcth ) + 1
        pthb(ll) = pthb(ll) + 1

        ctheta = abs( sum( evsij(:,1) * evtauij(:,1) ) )
        ll = floor( ctheta / bwcth ) + 1
        pthg(ll) = pthg(ll) + 1
  
      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  pcthphi =pcthphi *const/nfile
  pcthzeta=pcthzeta*const/nfile
  pphizeta=pphizeta*const/nfile

  p3d = p3d * const / nfile

  ptha = ptha * const / nfile
  pthb = pthb * const / nfile
  pthg = pthg * const / nfile

  write(*,*) 'check pcthphi: ', sum(pcthphi)
  write(*,*) 'check pcthzeta:', sum(pcthzeta)
  write(*,*) 'check pphizeta:', sum(pphizeta)
  write(*,*) 'check p3d:     ', sum(p3d)
  write(*,*) 'check ptha:    ', sum( ptha )

  ptha = ptha / bwcth
  pthb = pthb / bwcth 
  pthg = pthg / bwcth

  pcthphi = pcthphi/bwcth/bwphi
  pcthzeta=pcthzeta/bwcth/bwzeta
  pphizeta=pphizeta/bwphi/bwzeta

  p3d = p3d / bwcth / bwphi / bwzeta

  open(15,file='jpalign-1d'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,*) (ii-.5)*bwcth, ptha(ii), pthb(ii), pthg(ii)
    end do
  close(15)

  open(15, file = 'jpalign-3d'//str(1:len_trim(str)))
    write(15,'(''zone t="3d align pdf", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') npnt,npnt,npnt
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
      p3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)


  open(15,file='jpalign-2d'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(costheta,phi)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5)*bwcth, (jj-.5)*bwphi, pcthphi(ii,jj)
    end do
    end do

    write(15,*) 'Zone T= "(costheta,zeta)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5)*bwcth, (jj-.5)*bwzeta, pcthzeta(ii,jj)
    end do
    end do

    write(15,*) 'Zone T= "(phi,zeta)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5)*bwphi, (jj-.5)*bwzeta, pphizeta(ii,jj)
    end do
    end do
  close(15)

  deallocate(s11,s12,s13,s22,s23,s33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz,wx,wy,wz)
  deallocate(dnsux, dnsuy, dnsuz)

  call destroyplan3d

  write(*,*) 'Finished.'

end program align3dst 
