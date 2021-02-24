program align3drtsim
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  complex, allocatable, dimension(:,:,:) :: ux, uy, uz
  real,    allocatable, dimension(:,:,:) :: g, gcut
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,bb,evrij,evtauij
  integer :: ierr

  integer, parameter :: npnt=40
  real, parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)
  real(dp), dimension(npnt,npnt) :: pcthphi,pcthzeta,pphizeta
  real(dp), dimension(npnt,npnt,npnt) :: p3d
  real(dp), dimension(npnt) :: ptha, pthb, pthg

  real(dp) :: const
  real :: delta_c, ignore_me, kcut2
  real :: ctheta,cos1,cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of the angles between eigenvectors of rij and Leonard stress <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3drt-sub.x nx filelist ndel'
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

  str='-'//str(1:len_trim(str))//'dx-rtsub-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx
  kcut2 = ( pi / (delta_c/2) )**2

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),gcut(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  where ( g .ge. kcut2 )
    gcut = 0.
  elsewhere
    gcut = 1.
  endwhere
  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  pcthphi = 0.d0; pcthzeta = 0.d0; pphizeta = 0.d0
  ptha = 0.d0; pthb = 0.d0; pthg = 0.d0
  p3d = 0.d0
  nfile=0
  open(20,file=fnm(1:len_trim(fnm)))

    do while ( .not. eof(20) ) 
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

      ux = ux * gcut
      uy = uy * gcut
      uz = uz * gcut

      call rsgstauij(ux,ux,t11,g,nx,ny,nz) 
      call rsgstauij(ux,uy,t12,g,nx,ny,nz)
      call rsgstauij(ux,uz,t13,g,nx,ny,nz)
      call rsgstauij(uy,uy,t22,g,nx,ny,nz)
      call rsgstauij(uy,uz,t23,g,nx,ny,nz)
      call rsgstauij(uz,uz,t33,g,nx,ny,nz)

      ux=ux*g
      uy=uy*g 
      uz=uz*g

      call rsgstauij(ux,ux,r11,g,nx,ny,nz) 
      call rsgstauij(ux,uy,r12,g,nx,ny,nz)
      call rsgstauij(ux,uz,r13,g,nx,ny,nz)
      call rsgstauij(uy,uy,r22,g,nx,ny,nz)
      call rsgstauij(uy,uz,r23,g,nx,ny,nz)
      call rsgstauij(uz,uz,r33,g,nx,ny,nz)

      t11 = t11 - r11
      t12 = t12 - r12
      t13 = t13 - r13
      t22 = t22 - r22
      t23 = t23 - r23
      t33 = t33 - r33
  
      r11 = -(t11+t22+t33) / 3.

      t11 = t11 + r11
      t22 = t22 + r11
      t33 = t33 + r11

      call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)
      t11 = t11 * gcut * const
      t12 = t12 * gcut * const
      t13 = t13 * gcut * const
      t22 = t22 * gcut * const
      t23 = t23 * gcut * const
      t33 = t33 * gcut * const 
      call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        r11(ii,jj,kk)=eye*(ky(jj)*uz(ii,jj,kk)-kz(kk)*uy(ii,jj,kk))
        r22(ii,jj,kk)=eye*(kz(kk)*ux(ii,jj,kk)-kx(ii)*uz(ii,jj,kk))
        r33(ii,jj,kk)=eye*(kx(ii)*uy(ii,jj,kk)-ky(jj)*ux(ii,jj,kk))

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

          ll = (ii + 1) / 2

          cc(1,1)=real(r11(ll,jj,kk))
          cc(1,2)=real(r12(ll,jj,kk))
          cc(1,3)=real(r13(ll,jj,kk))
          cc(2,2)=real(r22(ll,jj,kk))
          cc(2,3)=real(r23(ll,jj,kk))
          cc(3,3)=real(r33(ll,jj,kk))

          bb(1,1)=-real(t11(ll,jj,kk))
          bb(1,2)=-real(t12(ll,jj,kk))
          bb(1,3)=-real(t13(ll,jj,kk))
          bb(2,2)=-real(t22(ll,jj,kk))
          bb(2,3)=-real(t23(ll,jj,kk))
          bb(3,3)=-real(t33(ll,jj,kk))
        else
          ll = ii / 2

          cc(1,1)=aimag(r11(ll,jj,kk))
          cc(1,2)=aimag(r12(ll,jj,kk))
          cc(1,3)=aimag(r13(ll,jj,kk))
          cc(2,2)=aimag(r22(ll,jj,kk))
          cc(2,3)=aimag(r23(ll,jj,kk))
          cc(3,3)=aimag(r33(ll,jj,kk))

          bb(1,1)=-aimag(t11(ll,jj,kk))
          bb(1,2)=-aimag(t12(ll,jj,kk))
          bb(1,3)=-aimag(t13(ll,jj,kk))
          bb(2,2)=-aimag(t22(ll,jj,kk))
          bb(2,3)=-aimag(t23(ll,jj,kk))
          bb(3,3)=-aimag(t33(ll,jj,kk))
        end if
        cc(2,1)=cc(1,2)
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
  
        bb(2,1)=bb(1,2)
        bb(3,1)=bb(1,3)
        bb(3,2)=bb(2,3)

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
      write(15,'(15E15.4)') (ii-.5)*bwcth, ptha(ii), pthb(ii), pthg(ii)
    end do
  close(15)

  open(15, file = 'jpalign-3d'//str(1:len_trim(str)))
    write(15,'(''zone t="3d align pdf", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') npnt,npnt,npnt
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
      p3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  open(15,file='jpalign-2d'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(costheta,phi)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5)*bwcth, (jj-.5)*bwphi, pcthphi(ii,jj)
    end do
    end do

    write(15,*) 'Zone T= "(costheta,zeta)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5)*bwcth, (jj-.5)*bwzeta, pcthzeta(ii,jj)
    end do
    end do

    write(15,*) 'Zone T= "(phi,zeta)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.4)') (ii-.5)*bwphi, (jj-.5)*bwzeta, pphizeta(ii,jj)
    end do
    end do
  close(15)

  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,gcut,kx,ky,kz,ux,uy,uz)

  call destroyplan3d

  write(*,*) 'finished'

end program align3drtsim 
