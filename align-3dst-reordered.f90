program align3dst
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel,mx,my,mz
  
  complex, allocatable, dimension(:,:,:) :: wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real,    allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1, fv2 
  real(dp), dimension(3,3) :: cc,evsij,evtauij
  real(dp), dimension(3) :: cosines
  integer :: ierr

  integer, parameter :: npnt=40
  real, parameter :: bwcth=1./npnt, bwphi=pi/2./npnt, bwzeta=pi/2./npnt
  real, dimension(npnt,npnt) :: pcthphi,pcthzeta,pphizeta
  real, dimension(npnt,npnt,npnt) :: p3d
  real, dimension(npnt) :: pthz,pthp,pthm

  real :: const, delta_c, ignore_me
  real :: ctheta,cos1,cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of the angles between reordered eigenvectors of sij and sgs stress <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3dst-reordered.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files'
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

  str='-'//str(1:len_trim(str))//'dx-reorderst-'//fnm(1:len_trim(fnm))//'.dat'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
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
  pthz= 0.; pthp=0.; pthm=0.
  p3d = 0.
  nfile=0
  open(20,file=fnm(1:len_trim(fnm))//'.list')

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

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
      write(*,*) 'after reading data'

  
      call rsgstauij(s11,s11,t11,g,nx,ny,nz) 
      call rsgstauij(s11,s22,t12,g,nx,ny,nz)
      call rsgstauij(s11,s33,t13,g,nx,ny,nz)
      call rsgstauij(s22,s22,t22,g,nx,ny,nz)
      call rsgstauij(s22,s33,t23,g,nx,ny,nz)
      call rsgstauij(s33,s33,t33,g,nx,ny,nz)
      write(*,*) 'after sgs stress'

      wx=-(t11+t22+t33)/3.
      t11=t11+wx
      t22=t22+wx
      t33=t33+wx
      write(*,*) 'trace removed'
  
      s11=s11*g
      s22=s22*g 
      s33=s33*g

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        wx(ii,jj,kk)=eye*(ky(jj)*s33(ii,jj,kk)-kz(kk)*s22(ii,jj,kk))
        wy(ii,jj,kk)=eye*(kz(kk)*s11(ii,jj,kk)-kx(ii)*s33(ii,jj,kk))
        wz(ii,jj,kk)=eye*(kx(ii)*s22(ii,jj,kk)-ky(jj)*s11(ii,jj,kk))
        s12(ii,jj,kk)=.5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
        s13(ii,jj,kk)=.5*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
        s23(ii,jj,kk)=.5*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))
        s11(ii,jj,kk)=eye*kx(ii)*s11(ii,jj,kk)
        s22(ii,jj,kk)=eye*ky(jj)*s22(ii,jj,kk)
        s33(ii,jj,kk)=eye*kz(kk)*s33(ii,jj,kk)
      end do
      end do
      end do
      write(*,*) 'vorticity and strain rate'
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
      write(*,*) 'fftw vorticity'

      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
      write(*,*) 'fftw sij'
  
  
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx
  
        cc(1,1)=real(s11(ii,jj,kk))
        cc(1,2)=real(s12(ii,jj,kk))
        cc(1,3)=real(s13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=real(s22(ii,jj,kk))
        cc(2,3)=real(s23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=real(s33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evsij,fv1,fv2,ierr)
        do ndel=1,3
          evsij(:,ndel)=evsij(:,ndel)/sqrt(sum(evsij(:,ndel)**2))
        end do
        cosines(3) = real(wx(ii,jj,kk))*evsij(1,3)+real(wy(ii,jj,kk))*evsij(2,3)+real(wz(ii,jj,kk))*evsij(3,3)
        cosines(2) = real(wx(ii,jj,kk))*evsij(1,2)+real(wy(ii,jj,kk))*evsij(2,2)+real(wz(ii,jj,kk))*evsij(3,2)
        cosines(1) = real(wx(ii,jj,kk))*evsij(1,1)+real(wy(ii,jj,kk))*evsij(2,1)+real(wz(ii,jj,kk))*evsij(3,1)
        cosines = abs(cosines)

        ! call reorders
        ! reordering 
        ll = maxloc(cosines,dim=1)
        select case (ll)
        case (1)
          if (evalues(2) .ge. evalues(3)) then
            mm = 2
            nn = 3
          else 
            mm = 3
            nn = 2
          end if
        case (2)
          if (evalues(1) .ge. evalues(3)) then
            mm = 1
            nn = 3
          else
            mm = 3
            nn = 1
          end if
        case (3)
          if (evalues(1) .ge. evalues(2)) then
            mm = 1
            nn = 2
          else
            mm = 2
            nn = 1
          end if
        end select

  
        cc(1,1)=-real(t11(ii,jj,kk))
        cc(1,2)=-real(t12(ii,jj,kk))
        cc(1,3)=-real(t13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=-real(t22(ii,jj,kk))
        cc(2,3)=-real(t23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=-real(t33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evtauij,fv1,fv2,ierr)
        do ndel=1,3
          evtauij(:,ndel)=evtauij(:,ndel)/sqrt(sum(evtauij(:,ndel)**2))
        end do
        cosines(3) = real(wx(ii,jj,kk))*evtauij(1,3)+real(wy(ii,jj,kk))*evtauij(2,3)+real(wz(ii,jj,kk))*evtauij(3,3)
        cosines(2) = real(wx(ii,jj,kk))*evtauij(1,2)+real(wy(ii,jj,kk))*evtauij(2,2)+real(wz(ii,jj,kk))*evtauij(3,2)
        cosines(1) = real(wx(ii,jj,kk))*evtauij(1,1)+real(wy(ii,jj,kk))*evtauij(2,1)+real(wz(ii,jj,kk))*evtauij(3,1)
        cosines = abs(cosines)

!        call reordert

       ! reordering 
        mx = maxloc(cosines,dim=1)
        select case (mx)
        case (1)
          if (evalues(2) .ge. evalues(3)) then
            my = 2
            mz = 3
          else 
            my = 3
            mz = 2
          end if
        case (2)
          if (evalues(1) .ge. evalues(3)) then
            my = 1
            mz = 3
          else
            my = 3
            mz = 1
          end if
        case (3)
          if (evalues(1) .ge. evalues(2)) then
            my = 1
            mz = 2
          else
            my = 2
            mz = 1
          end if
        end select


!        call pdf

        ctheta=abs(sum(evsij(:,ll)*evtauij(:,mx))) 
        cos1 = sum(evtauij(:,mx)*evsij(:,mm))
        cos2 = sum(evtauij(:,mx)*evsij(:,nn))
        cosphi= cos1/sqrt(cos1*cos1+cos2*cos2)
        phi=acos(abs(cosphi))

        ! Tao etal used the third eigen vector
        !cos3 = sum(evsij(:,nn)*evtauij(:,my))
        !cos4 = sum(evsij(:,nn)*evtauij(:,mz))
        !coszeta = cos4/sqrt(cos3*cos3+cos4*cos4)

        ! Horiuti used the second 
        cos3 = sum(evsij(:,mm)*evtauij(:,my))
        cos4 = sum(evsij(:,mm)*evtauij(:,mz))
        coszeta = cos3/sqrt(cos3*cos3+cos4*cos4)
        zeta = acos(abs(coszeta))
  
  
        cos1 = abs(sum(evsij(:,mm)*evtauij(:,my)))
        ll=floor(cos1/bwcth)+1
        pthp(ll)=pthp(ll)+1
        cos1 = abs(sum(evsij(:,nn)*evtauij(:,mz)))
        ll=floor(cos1/bwcth)+1
        pthm(ll)=pthm(ll)+1

        ll=floor(ctheta/bwcth)+1
        pthz(ll)=pthz(ll)+1

        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1

  
        pcthphi(ll,mm)=pcthphi(ll,mm)+1
        pcthzeta(ll,nn)=pcthzeta(ll,nn)+1
        pphizeta(mm,nn)=pphizeta(mm,nn)+1

        p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1

  
        cc(1,1)=aimag(s11(ii,jj,kk))
        cc(1,2)=aimag(s12(ii,jj,kk))
        cc(1,3)=aimag(s13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=aimag(s22(ii,jj,kk))
        cc(2,3)=aimag(s23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=aimag(s33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evsij,fv1,fv2,ierr)
        do ndel=1,3
          evsij(:,ndel)=evsij(:,ndel)/sqrt(sum(evsij(:,ndel)**2))
        end do
        cosines(3) = aimag(wx(ii,jj,kk))*evsij(1,3)+aimag(wy(ii,jj,kk))*evsij(2,3)+aimag(wz(ii,jj,kk))*evsij(3,3)
        cosines(2) = aimag(wx(ii,jj,kk))*evsij(1,2)+aimag(wy(ii,jj,kk))*evsij(2,2)+aimag(wz(ii,jj,kk))*evsij(3,2)
        cosines(1) = aimag(wx(ii,jj,kk))*evsij(1,1)+aimag(wy(ii,jj,kk))*evsij(2,1)+aimag(wz(ii,jj,kk))*evsij(3,1)
        cosines = abs(cosines)

        !call reorders
        ! reordering 
        ll = maxloc(cosines,dim=1)
        select case (ll)
        case (1)
          if (evalues(2) .ge. evalues(3)) then
            mm = 2
            nn = 3
          else 
            mm = 3
            nn = 2
          end if
        case (2)
          if (evalues(1) .ge. evalues(3)) then
            mm = 1
            nn = 3
          else
            mm = 3
            nn = 1
          end if
        case (3)
          if (evalues(1) .ge. evalues(2)) then
            mm = 1
            nn = 2
          else
            mm = 2
            nn = 1
          end if
        end select

  
        cc(1,1)=-aimag(t11(ii,jj,kk))
        cc(1,2)=-aimag(t12(ii,jj,kk))
        cc(1,3)=-aimag(t13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=-aimag(t22(ii,jj,kk))
        cc(2,3)=-aimag(t23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=-aimag(t33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evtauij,fv1,fv2,ierr)
        do ndel=1,3
          evtauij(:,ndel)=evtauij(:,ndel)/sqrt(sum(evtauij(:,ndel)**2))
        end do
        cosines(3) = aimag(wx(ii,jj,kk))*evtauij(1,3)+aimag(wy(ii,jj,kk))*evtauij(2,3)+aimag(wz(ii,jj,kk))*evtauij(3,3)
        cosines(2) = aimag(wx(ii,jj,kk))*evtauij(1,2)+aimag(wy(ii,jj,kk))*evtauij(2,2)+aimag(wz(ii,jj,kk))*evtauij(3,2)
        cosines(1) = aimag(wx(ii,jj,kk))*evtauij(1,1)+aimag(wy(ii,jj,kk))*evtauij(2,1)+aimag(wz(ii,jj,kk))*evtauij(3,1)
        cosines = abs(cosines)

!        call reordert

        ! reordering 
        mx = maxloc(cosines,dim=1)
        select case (mx)
        case (1)
          if (evalues(2) .ge. evalues(3)) then
            my = 2
            mz = 3
          else 
            my = 3
            mz = 2
           end if
         case (2)
          if (evalues(1) .ge. evalues(3)) then
            my = 1
            mz = 3
          else
            my = 3
            mz = 1
          end if
        case (3)
          if (evalues(1) .ge. evalues(2)) then
            my = 1
            mz = 2
          else
            my = 2
            mz = 1
          end if
        end select

  
!        call pdf

        ctheta=abs(sum(evsij(:,ll)*evtauij(:,mx))) 
        cos1 = sum(evtauij(:,mx)*evsij(:,mm))
        cos2 = sum(evtauij(:,mx)*evsij(:,nn))
        cosphi= cos1/sqrt(cos1*cos1+cos2*cos2)
        phi=acos(abs(cosphi))
        ! Tao use the third eigven vector
        !cos3 = sum(evsij(:,nn)*evtauij(:,my))
        !cos4 = sum(evsij(:,nn)*evtauij(:,mz))
        !coszeta = cos4/sqrt(cos3*cos3+cos4*cos4)
        ! Horiuti used the second
        cos3 = sum(evsij(:,mm)*evtauij(:,my))
        cos4 = sum(evsij(:,mm)*evtauij(:,mz))
        coszeta = cos3/sqrt(cos3*cos3+cos4*cos4)

        zeta = acos(abs(coszeta))
 
 
        cos1 = abs(sum(evsij(:,mm)*evtauij(:,my)))
        ll=floor(cos1/bwcth)+1
        pthp(ll)=pthp(ll)+1
        cos1 = abs(sum(evsij(:,nn)*evtauij(:,mz)))
        ll=floor(cos1/bwcth)+1
        pthm(ll)=pthm(ll)+1

        ll=floor(ctheta/bwcth)+1
        pthz(ll)=pthz(ll)+1

        mm=floor(phi/bwphi)+1
        nn=floor(zeta/bwzeta)+1
  
        pcthphi(ll,mm)=pcthphi(ll,mm)+1
        pcthzeta(ll,nn)=pcthzeta(ll,nn)+1
        pphizeta(mm,nn)=pphizeta(mm,nn)+1

        p3d(ll,mm,nn) = p3d(ll,mm,nn) + 1

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)


  pthz=pthz*const/nfile
  pthp=pthp*const/nfile
  pthm=pthm*const/nfile
  write(*,*) 'check pthz: ', sum(pthz)
  write(*,*) 'check pthp: ', sum(pthp)
  write(*,*) 'check pthm: ', sum(pthm)
  pthz=pthz/bwcth
  pthp=pthp/bwcth
  pthm=pthm/bwcth

  open(15,file='jpalign-1d'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,*) (ii-.5)*bwcth, pthz(ii), pthp(ii), pthm(ii)
    end do
  close(15)


  pcthphi =pcthphi *const/nfile
  pcthzeta=pcthzeta*const/nfile
  pphizeta=pphizeta*const/nfile


  write(*,*) 'check pcthphi: ', sum(pcthphi)
  write(*,*) 'check pcthzeta:', sum(pcthzeta)
  write(*,*) 'check pphizeta:', sum(pphizeta)

  pcthphi = pcthphi/bwcth/bwphi
  pcthzeta=pcthzeta/bwcth/bwzeta
  pphizeta=pphizeta/bwphi/bwzeta

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

  p3d = p3d * const / nfile

  write(*,*) 'check p3d: ', sum(p3d)

  p3d = p3d / bwcth / bwphi / bwzeta

  open(15, file = 'jpalign-3d'//str(1:len_trim(str)))
    write(15,*) 'zone t="3d align pdf", i=',npnt,', j=', npnt, ', k=', npnt, ', f=points'
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
      p3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  deallocate(wx,wy,wz,s11,s12,s13,s22,s23,s33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz)

  call destroyplan3d

  write(*,*) 'done'

end program align3dst 
!
!
!
!
!
!
!
!
!
 
 
!
!
!
!
!
!
!
!
!
!
!
