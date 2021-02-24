program align3drtcndpss
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nthresh,ndel,nfile
  integer(8) :: numpnt
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  complex, allocatable, dimension(:,:,:) :: wx, wy, wz, ss
  real,    allocatable, dimension(:,:,:) :: g
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc, bb, evrij, evtauij
  integer :: ierr

  integer, parameter :: npnt=40
  real, parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)
  real(dp), dimension(npnt,npnt) :: pcthphi,pcthzeta,pphizeta
  real(dp), dimension(npnt,npnt,npnt) :: p3d
  real(dp), dimension(npnt) :: pcth, pphi, pzeta

  real(dp) :: rmss

  real :: delta_c, ignore_me, const, rmss0
  real :: ctheta, cos1, cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Joint PDFs of the angles between eigenvectors of rij and sgs stress <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./align-3drt-cndpss.x nx filelist ndel nthresh'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        nthresh: defining the threshold value for ss'
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

  ! threshold value for condition
  call getarg(4,fpath)
  read(fpath, '(I20)') nthresh
  fpath=adjustl(fpath)

  str='-'//str(1:len_trim(str))//'dx-rt-cndpss'//fpath(1:len_trim(fpath))//'-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny = nx; nz = nx
  lx = nx/2; ly = ny; lz= nz; lx1 =lx + 1
  const = 1. / (nx * ny * nz)

  delta_c = ndel * 2 * pi / nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),wx(lx1,ly,lz), wy(lx1,ly,lz), wz(lx1,ly,lz) )
  allocate(ss(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1:len_trim(fnm)) )

    pcthphi = 0.d0; pcthzeta = 0.d0; pphizeta = 0.d0
    pcth = 0.d0; pphi = 0.d0; pzeta = 0.d0
    p3d = 0.d0

    rmss = 0.0_dp

    numpnt = 0
    nfile = 0
    do while ( .not. eof(20) ) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))


      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r11
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r22
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r33
      close(10)


      write(*,*) 'after reading data files'

      call rsgstauij(r11,r11,t11,g,nx,ny,nz) 
      call rsgstauij(r11,r22,t12,g,nx,ny,nz)
      call rsgstauij(r11,r33,t13,g,nx,ny,nz)
      call rsgstauij(r22,r22,t22,g,nx,ny,nz)
      call rsgstauij(r22,r33,t23,g,nx,ny,nz)
      call rsgstauij(r33,r33,t33,g,nx,ny,nz)

      r12 = -(t11+t22+t33) / 3.

      t11 = t11 + r12
      t22 = t22 + r12
      t33 = t33 + r12
  
      ! Filtered velocity
      r11=r11*g
      r22=r22*g 
      r33=r33*g

      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ! Vorticity
        wx(ii,jj,kk)=eye*(ky(jj)*r33(ii,jj,kk)-kz(kk)*r22(ii,jj,kk))
        wy(ii,jj,kk)=eye*(kz(kk)*r11(ii,jj,kk)-kx(ii)*r33(ii,jj,kk))
        wz(ii,jj,kk)=eye*(kx(ii)*r22(ii,jj,kk)-ky(jj)*r11(ii,jj,kk))

      ! rij is now temporarily sij
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
  
      ! squared magnitude of sij defined as : 2 sij sij.
      ss = cmplx( real(r11) * real(r11) + real(r22) * real(r22) + real(r33) * real(r33) &
                 + 2 * ( real(r12) * real(r12) + real(r13) * real(r13) + real(r23) * real(r23) ), &
                  aimag(r11) * aimag(r11) + aimag(r22) * aimag(r22) + aimag(r33) * aimag(r33) &
                 + 2 * ( aimag(r12) * aimag(r12) + aimag(r13) * aimag(r13) + aimag(r23) * aimag(r23) ) &
                 )
      ss = 2.*ss
    
      if ( numpnt .eq. 0 ) then
        rmss0 = sum( real( ss(1:lx,:,:) ) ) + sum( aimag( ss(1:lx,:,:) ) )
        rmss0 = sqrt( rmss0 / (nx*ny*nz) )
    
        write(*,*) 'Estimated rms of sij: ', rmss0
      end if
      rmss = rmss + sum( real( ss(1:lx,:,:) ) ) + sum( aimag( ss(1:lx,:,:) ) )
      nfile = nfile + 1

  
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        ! rij 
        r11(ii,jj,kk)=eye*kx(ii)*wx(ii,jj,kk)
        r22(ii,jj,kk)=eye*ky(jj)*wy(ii,jj,kk)
        r33(ii,jj,kk)=eye*kz(kk)*wz(ii,jj,kk)
        r12(ii,jj,kk)=.5*eye*(kx(ii)*wy(ii,jj,kk)+ky(jj)*wx(ii,jj,kk))
        r13(ii,jj,kk)=.5*eye*(kx(ii)*wz(ii,jj,kk)+kz(kk)*wx(ii,jj,kk))
        r23(ii,jj,kk)=.5*eye*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
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
          ll = (ii + 1)/2
          ignore_me = sqrt( real( ss(ll,jj,kk) ) )
        else
          ll = ii / 2
          ignore_me = sqrt( aimag( ss(ll,jj,kk) ) )
        end if

        if ( ignore_me .ge. nthresh * rmss0 ) then

          numpnt = numpnt + 1
  
          if ( mod(ii,2) .eq. 1 ) then
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
 
          pcth(ll) = pcth(ll) + 1
          pphi(mm) = pphi(mm) + 1
          pzeta(nn) = pzeta(nn) + 1

        end if

      end do
      end do
      end do

    end do
  close(20)

  rmss = rmss * const / nfile
  rmss = sqrt( rmss )
  write(*,*) 'Mean rms s :', rmss

  write(*,*) 'numpnt = ', numpnt
  write(*,*) 'percentage = ', numpnt * const / nfile

  pcthphi  = pcthphi  / numpnt
  pcthzeta = pcthzeta / numpnt
  pphizeta = pphizeta / numpnt
  

  p3d = p3d / numpnt

  pcth = pcth / numpnt
  pphi = pphi / numpnt
  pzeta = pzeta / numpnt

  write(*,*) 'check pcthphi: ', sum(pcthphi)
  write(*,*) 'check pcthzeta:', sum(pcthzeta)
  write(*,*) 'check pphizeta:', sum(pphizeta)
  write(*,*) 'check p3d:     ', sum(p3d)
  write(*,*) 'check pcth:    ', sum( pcth )
  write(*,*) 'check pphi:    ', sum( pphi )
  write(*,*) 'check pzeta:    ', sum( pzeta )

  pcthphi = pcthphi/bwcth/bwphi
  pcthzeta=pcthzeta/bwcth/bwzeta
  pphizeta=pphizeta/bwphi/bwzeta

  p3d = p3d / bwcth / bwphi / bwzeta

  pcth = pcth / bwcth
  pphi = pphi / bwphi
  pzeta = pzeta / bwzeta

  open(15,file='jpalign-1d'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, pcth(ii), (ii-.5)*bwphi, pphi(ii), (ii-.5)*bwzeta, pzeta(ii)
    end do
  close(15)

  open(15, file = 'jpalign-3d'//str(1:len_trim(str)))
    write(15,'(''zone t="3d align pdf", i='',I4,'', j='', I4, '', k='', I4, '', f=point'')') npnt,npnt,npnt
    do kk=1,npnt
    do jj=1,npnt
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5) * bwcth, (jj-.5) * bwphi, (kk-.5) * bwzeta, &
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
  deallocate(g,kx,ky,kz,wx,wy,wz,ss)

  call destroyplan3d

  write(*,*) 'finished'

end program align3drtcndpss 
