program align3drprm
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  
  complex, allocatable, dimension(:,:,:) :: hpx,hpy,hpz
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real,    allocatable, dimension(:,:,:) :: g,kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evrij,evtauij
  integer :: ierr

  integer, parameter :: npnt=50
  real, parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)
  real, dimension(npnt,npnt) :: pcthphi,pcthzeta,pphizeta
  real, dimension(npnt,npnt,npnt) :: p3d
  real, dimension(npnt) :: ptha, pthb, pthg

  real :: const, delta_c, ignore_me
  real :: ctheta,cos1,cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of the angles between eigenvectors of rij(p) and rij(m) <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3drprm.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
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

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  str='-'//str(1:len_trim(str))//'dx-rprm-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  allocate(hpx(lx1,ly,lz),hpy(lx1,ly,lz),hpz(lx1,ly,lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  call heliwave(hpx,hpy,hpz,kx(:,1,1),ky(1,:,1),kz(1,1,:),g,lx1,ly,lz)
  write(*,*) 'after heliwave'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  pcthphi = 0.; pcthzeta = 0.; pphizeta = 0.
  ptha = 0.; pthb = 0.; pthg = 0.
  p3d = 0.
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

      r12 = r12 * g
      r13 = r13 * g
      r23 = r23 * g
  
      ! r11 is the ap component
      r11 = r12 * conjg(hpx) + r13 * conjg(hpy) + r23 * conjg(hpz)

      ! t11 is the am component
      t11 = r12 * hpx + r13 * hpy + r23 * hpz

      r12=r11*hpx
      r13=r11*hpy 
      r23=r11*hpz
      r11=eye*(ky*r23-kz*r13)
      r22=eye*(kz*r12-kx*r23)
      r33=eye*(kx*r13-ky*r12)
      r12=.5*eye*(kx*r22+ky*r11)
      r13=.5*eye*(kx*r33+kz*r11)
      r23=.5*eye*(ky*r33+kz*r22)
      r11=eye*kx*r11
      r22=eye*ky*r22
      r33=eye*kz*r33
  
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
  

      t12=t11*conjg(hpx)
      t13=t11*conjg(hpy) 
      t23=t11*conjg(hpz)
      t11=eye*(ky*t23-kz*t13)
      t22=eye*(kz*t12-kx*t23)
      t33=eye*(kx*t13-ky*t12)
      t12=.5*eye*(kx*t22+ky*t11)
      t13=.5*eye*(kx*t33+kz*t11)
      t23=.5*eye*(ky*t33+kz*t22)
      t11=eye*kx*t11
      t22=eye*ky*t22
      t33=eye*kz*t33
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)
  

  
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx
  
        cc(1,1)=real(r11(ii,jj,kk))
        cc(1,2)=real(r12(ii,jj,kk))
        cc(1,3)=real(r13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=real(r22(ii,jj,kk))
        cc(2,3)=real(r23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=real(r33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
        cc(1,1)=real(t11(ii,jj,kk))
        cc(1,2)=real(t12(ii,jj,kk))
        cc(1,3)=real(t13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=real(t22(ii,jj,kk))
        cc(2,3)=real(t23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=-real(t33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evtauij,fv1,fv2,ierr)
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
  
  
        ll=floor(ctheta/bwcth)+1
        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1
  
        pcthphi(ll,mm)=pcthphi(ll,mm)+1
        pcthzeta(ll,nn)=pcthzeta(ll,nn)+1
        pphizeta(mm,nn)=pphizeta(mm,nn)+1

        p3d(ll,mm,nn)=p3d(ll,mm,nn)+1

        ptha(ll) = ptha(ll) + 1

        ctheta = abs( sum( evrij(:,2) * evtauij(:,2) ) ) 
        ll = floor( ctheta / bwcth ) + 1
        pthb(ll) = pthb(ll) + 1

        ctheta = abs( sum( evrij(:,1) * evtauij(:,1) ) )
        ll = floor( ctheta / bwcth ) + 1
        pthg(ll) = pthg(ll) + 1
  
  
        cc(1,1)=aimag(r11(ii,jj,kk))
        cc(1,2)=aimag(r12(ii,jj,kk))
        cc(1,3)=aimag(r13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=aimag(r22(ii,jj,kk))
        cc(2,3)=aimag(r23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=aimag(r33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
        cc(1,1)=aimag(t11(ii,jj,kk))
        cc(1,2)=aimag(t12(ii,jj,kk))
        cc(1,3)=aimag(t13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=aimag(t22(ii,jj,kk))
        cc(2,3)=aimag(t23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=aimag(t33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evtauij,fv1,fv2,ierr)
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
  
  
        ll=floor(ctheta/bwcth)+1
        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1
  
        pcthphi(ll,mm)=pcthphi(ll,mm)+1
        pcthzeta(ll,nn)=pcthzeta(ll,nn)+1
        pphizeta(mm,nn)=pphizeta(mm,nn)+1

        p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1
  
        ptha(ll) = ptha(ll) + 1

        ctheta = abs( sum( evrij(:,2) * evtauij(:,2) ) ) 
        ll = floor( ctheta / bwcth ) + 1
        pthb(ll) = pthb(ll) + 1

        ctheta = abs( sum( evrij(:,1) * evtauij(:,1) ) )
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
  write(*,*) 'check pthb:    ', sum( pthb )
  write(*,*) 'check pthg:    ', sum( pthg )

  pcthphi = pcthphi/bwcth/bwphi
  pcthzeta=pcthzeta/bwcth/bwzeta
  pphizeta=pphizeta/bwphi/bwzeta

  p3d = p3d / bwcth / bwphi / bwzeta

  ptha = ptha / bwcth
  pthb = pthb / bwcth
  pthg = pthg / bwcth

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

  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz,hpx,hpy,hpz)

  call destroyplan3d

  write(*,*) 'finished'

end program align3drprm 
