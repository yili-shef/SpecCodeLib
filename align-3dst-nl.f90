program align3drtnl
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  real, parameter :: cnl=1./12.
  
  complex, allocatable, dimension(:,:,:) :: ux,uy,uz,wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real,    allocatable, dimension(:,:,:) :: k2,g,kx,ky,kz

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
  write(*,'(''>>> Joint PDFs of the angles between eigenvectors of rij and nonlinear model stress <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3drt-nl.x nx filelist ndel'
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

  str='-'//str(1:len_trim(str))//'dx-stnl-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(k2(lx1,ly,lz),kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-k2*delta_c**2/24.)

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
        read(10)ux
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'

      ux=ux*g
      uy=uy*g 
      uz=uz*g
  
      ! s11 used as temporary storage for the velocity gradient
      s11=eye*ky*ux; s12=eye*kz*ux
      s13=eye*kx*uy; s22=eye*kz*uy
      s23=eye*kx*uz; s33=eye*ky*uz
      wx=eye*kx*ux; wy=eye*ky*uy; wz=eye*kz*uz
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

      t11 = cmplx(real(wx)*real(wx) + real(s11)*real(s11) + real(s12)*real(s12), &
            aimag(wx)*aimag(wx) + aimag(s11)*aimag(s11) + aimag(s12)*aimag(s12))
      t12 = cmplx(real(wx)*real(s13) + real(s11)*real(wy) + real(s12)*real(s22), &
            aimag(wx)*aimag(s13) + aimag(s11)*aimag(wy) + aimag(s12)*aimag(s22))
      t13 = cmplx(real(wx)*real(s23) + real(s11)*real(s33) + real(s12)*real(wz), &
            aimag(wx)*aimag(s23) + aimag(s11)*aimag(s33) + aimag(s12)*aimag(wz))
      t22 = cmplx(real(s13)*real(s13) + real(wy)*real(wy) + real(s22)*real(s22), &
            aimag(s13)*aimag(s13) + aimag(wy)*aimag(wy) + aimag(s22)*aimag(s22))
      t23 = cmplx(real(s13)*real(s23) + real(wy)*real(s33) + real(s22)*real(wz), &
            aimag(s13)*aimag(s23) + aimag(wy)*aimag(s33) + aimag(s22)*aimag(wz))
      t33 = cmplx(real(s23)*real(s23) + real(s33)*real(s33) + real(wz)*real(wz), &
            aimag(s23)*aimag(s23) + aimag(s33)*aimag(s33) + aimag(wz)*aimag(wz))


      wx = -(t11+t22+t33) / 3.

      t11 = t11 + wx
      t22 = t22 + wx
      t33 = t33 + wx

      t11 = cnl*delta_c**2*t11
      t12 = cnl*delta_c**2*t12
      t13 = cnl*delta_c**2*t13
      t22 = cnl*delta_c**2*t22
      t23 = cnl*delta_c**2*t23
      t33 = cnl*delta_c**2*t33
  
      s11=eye*kx*ux
      s22=eye*ky*uy
      s33=eye*kz*uz
      s12=.5*eye*(kx*uy+ky*ux)
      s13=.5*eye*(kx*uz+kz*ux)
      s23=.5*eye*(ky*uz+kz*uy)
  
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
  
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
  
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
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
  
  
        cc(1,1)=aimag(s11(ii,jj,kk))
        cc(1,2)=aimag(s12(ii,jj,kk))
        cc(1,3)=aimag(s13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=aimag(s22(ii,jj,kk))
        cc(2,3)=aimag(s23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=aimag(s33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evrij,fv1,fv2,ierr)
        do ll=1,3
          evrij(:,ll)=evrij(:,ll)/sqrt(sum(evrij(:,ll)**2))
        end do
  
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

  deallocate(ux,uy,uz,wx,wy,wz,s11,s12,s13,s22,s23,s33)
  deallocate(t11,t12,t13,t22,t23,t33)

end program align3drtnl 
