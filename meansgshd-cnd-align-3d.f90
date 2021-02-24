program sgshdcnd3d
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  
  complex, allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real(sp),    allocatable, dimension(:,:,:) :: g, gcut
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,bb,evrij,evtauij
  integer :: ierr

  real(sp), parameter :: mdissrate0 = 0.26158 ! Delta = 16 dx
!  real(sp), parameter :: mdissrate0 = 0.2375 ! Delta =  8 dx

  integer, parameter :: npnt = 40
  real(sp), parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)

  real(dp), dimension(npnt,npnt) :: jpdfcthphi, cndhdcthphi
  real(dp), dimension(npnt,npnt) :: jpdfcthzeta, cndhdcthzeta
  real(dp), dimension(npnt,npnt) :: jpdfphizeta, cndhdphizeta
  real(dp), dimension(npnt,npnt,npnt) :: jpdf3d, cndhd3d

  real(dp) :: const, mdissrate

  real :: delta_c, ignore_me, pih, kcut2
  real :: ctheta, cos1, cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Mean SGS helicity dissipation conditioned on alignment <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgshd-cnd-align-3d.x nx filelist ndel'
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

  str='-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny=nx; nz=nx
  lx=nx/2; ly=ny;lz=nz;lx1=lx+1
  const=1./(nx*ny*nz)

  delta_c=ndel*2*pi/nx
  kcut2 = ( pi / (delta_c/2) )**2

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

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

  jpdf3d = 0.d0; cndhd3d = 0.d0
  jpdfcthphi = 0.d0; cndhdcthphi = 0.d0
  jpdfcthzeta = 0.d0; cndhdcthzeta = 0.d0
  jpdfphizeta = 0.d0; cndhdphizeta = 0.d0
  nfile=0
  mdissrate = 0._dp
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

      r12 = r12 * gcut
      r13 = r13 * gcut
      r23 = r23 * gcut
  
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

        if ( ll .ge. 1 .and. ll .le. npnt .and. &
             mm .ge. 1 .and. mm .le. npnt .and. &
             nn .ge. 1 .and. nn .le. npnt ) then

          jpdf3d(ll,mm,nn) = jpdf3d(ll,mm,nn) + 1
          cndhd3d(ll,mm,nn) = cndhd3d(ll,mm,nn) + pih
 
          jpdfcthphi(ll,mm) = jpdfcthphi(ll,mm) + 1
          cndhdcthphi(ll,mm) = cndhdcthphi(ll,mm) + pih

          jpdfcthzeta(ll,nn) = jpdfcthzeta(ll,nn) + 1
          cndhdcthzeta(ll,nn) = cndhdcthzeta(ll,nn) + pih

          jpdfphizeta(mm,nn) = jpdfphizeta(mm,nn) + 1
          cndhdphizeta(mm,nn) = cndhdphizeta(mm,nn) + pih

        end if

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  mdissrate = mdissrate * const / nfile

  cndhd3d = cndhd3d / (jpdf3d + tiny) 
  cndhdcthphi = cndhdcthphi / (jpdfcthphi + tiny)
  cndhdcthzeta = cndhdcthzeta / (jpdfcthzeta + tiny)
  cndhdphizeta = cndhdphizeta / (jpdfphizeta + tiny)

  jpdf3d = jpdf3d * const / nfile
  jpdfcthphi = jpdfcthphi * const / nfile
  jpdfcthzeta = jpdfcthzeta * const / nfile
  jpdfphizeta = jpdfphizeta * const / nfile

  write(*,*) 'check jpdf3d: ', sum(jpdf3d)
  write(*,*) 'check jpdfcthphi:', sum(jpdfcthphi)
  write(*,*) 'check jpdfcthzeta:', sum(jpdfcthzeta)
  write(*,*) 'check jpdfphizeta:', sum(jpdfphizeta)
  write(*,*) 'check mean dissipation rate: ', mdissrate

  jpdf3d = jpdf3d / bwcth / bwphi / bwzeta
  jpdfcthphi = jpdfcthphi / bwcth / bwphi
  jpdfcthzeta = jpdfcthzeta / bwcth / bwzeta
  jpdfphizeta = jpdfphizeta / bwphi / bwzeta

  open(15,file='cndsgshd-align-3d'//str(1:len_trim(str)))
    write(15,'(''zone t = "cnd hd", i='', i6, '',j='', i6, ''k='', i6, '',f=point'')') npnt, npnt, npnt
    do kk = 1, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, (kk-.5)*bwzeta, &
                            jpdf3d(ii,jj,kk), cndhd3d(ii,jj,kk), &
                            jpdf3d(ii,jj,kk) * cndhd3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  open(15,file='cndsgshd-align-2d'//str(1:len_trim(str)))
    write(15,'(''zone t = "(cth, phi)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, jpdfcthphi(ii,jj), &
                            cndhdcthphi(ii,jj), jpdfcthphi(ii,jj) * &
                            cndhdcthphi(ii,jj)
    end do
    end do

    write(15,'(''zone t = "(cth, zeta)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwzeta, jpdfcthzeta(ii,jj), &
                            cndhdcthzeta(ii,jj), jpdfcthzeta(ii,jj) * &
                            cndhdcthzeta(ii,jj)
    end do
    end do

    write(15,'(''zone t = "(phi, zeta)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do jj = 1, npnt
    do ii = 1, npnt
      write(15,'(15E15.5)') (ii-.5)*bwphi, (jj-.5)*bwzeta, jpdfphizeta(ii,jj), &
                            cndhdphizeta(ii,jj), jpdfphizeta(ii,jj) * &
                            cndhdphizeta(ii,jj)
    end do
    end do
  close(15)


  deallocate(r11,r12,r13,r22,r23,r33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,gcut,kx,ky,kz)

  call destroyplan3d

  write(*,*) 'finished'

end program sgshdcnd3d 
