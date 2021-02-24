program sgsedcnd3d
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel
  
  complex(sp), allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,bb,evsij,evtauij
  integer :: ierr

  integer, parameter :: npnt = 160
  real(sp), parameter :: bwcth=1./npnt, bwphi=pi/(2.*npnt), bwzeta=pi/(2.*npnt)

  integer, parameter :: npnt3d = 40
  real(sp), parameter :: bwcth3d=1./npnt3d, bwphi3d=pi/(2.*npnt3d), bwzeta3d=pi/(2.*npnt3d)

  real(dp), dimension(npnt,npnt) :: jpdfcthphi, cndedcthphi
  real(dp), dimension(npnt,npnt) :: jpdfcthzeta, cndedcthzeta
  real(dp), dimension(npnt,npnt) :: jpdfphizeta, cndedphizeta
  real(dp), dimension(npnt3d,npnt3d,npnt3d) :: jpdf3d, cnded3d

  real(dp) :: const, mdissrate

  real(sp) :: delta_c, ignore_me, pie, mdissrate0
  real(sp) :: ctheta, cos1, cos2,cosphi,phi,cos3,cos4,coszeta,zeta
  character(80) :: fnm,str,str1,fpath

  write(*,*) 
  write(*,'(''>>> Mean SGS energy dissipation conditioned on alignment <<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./meansgsed-cnd-align-st-3d.x nx filelist ndel normdissrate'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files, *.list'
          write(*,*) '        ndel: filter scale Delta=ndel*dx'
          write(*,*) '        normdissrate: dissrate for normalization'
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

  ! dissipation rate normalization 
  call getarg(4, str)
  read(str, '(F15.6)') mdissrate0

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
  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  !where ( g .ge. kcut2 )
  !  gcut = 0.
  !elsewhere
  !  gcut = 1.
  !endwhere

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  jpdf3d = 0.d0; cnded3d = 0.d0
  jpdfcthphi = 0.d0; cndedcthphi = 0.d0
  jpdfcthzeta = 0.d0; cndedcthzeta = 0.d0
  jpdfphizeta = 0.d0; cndedphizeta = 0.d0
  nfile=0
  mdissrate = 0._dp
  open(20,file=fnm(1:len_trim(fnm)))

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
      write(*,*) 'after reading data files'

      call rsgstauij(s11,s11,t11,g,nx,ny,nz) 
      call rsgstauij(s11,s22,t12,g,nx,ny,nz)
      call rsgstauij(s11,s33,t13,g,nx,ny,nz)
      call rsgstauij(s22,s22,t22,g,nx,ny,nz)
      call rsgstauij(s22,s33,t23,g,nx,ny,nz)
      call rsgstauij(s33,s33,t33,g,nx,ny,nz)

      s12 = -(t11+t22+t33) / 3.

      t11 = t11 + s12
      t22 = t22 + s12
      t33 = t33 + s12
  
      s11=s11*g
      s22=s22*g 
      s33=s33*g
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        s12(ii,jj,kk)=.5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
        s13(ii,jj,kk)=.5*eye*(kz(kk)*s11(ii,jj,kk)+kx(ii)*s33(ii,jj,kk))
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
  
  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
  
        if ( mod(ii,2) .eq. 1) then
          ll = (ii + 1)/2
          cc(1,1)=real(s11(ll,jj,kk))
          cc(1,2)=real(s12(ll,jj,kk))
          cc(1,3)=real(s13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=real(s22(ll,jj,kk))
          cc(2,3)=real(s23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=real(s33(ll,jj,kk))
   
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
          cc(1,1)=aimag(s11(ll,jj,kk))
          cc(1,2)=aimag(s12(ll,jj,kk))
          cc(1,3)=aimag(s13(ll,jj,kk))
          cc(2,1)=cc(1,2)
          cc(2,2)=aimag(s22(ll,jj,kk))
          cc(2,3)=aimag(s23(ll,jj,kk))
          cc(3,1)=cc(1,3)
          cc(3,2)=cc(2,3)
          cc(3,3)=aimag(s33(ll,jj,kk))
   
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

        pie = sum(cc * bb)
        mdissrate = mdissrate + pie
  
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
  
        pie = pie / mdissrate0
  
        ll = floor(ctheta/bwcth)+1
        mm = floor(phi/bwphi)+1
        nn = floor(zeta/bwzeta)+1

        if ( ll .ge. 1 .and. ll .le. npnt .and. &
             mm .ge. 1 .and. mm .le. npnt .and. &
             nn .ge. 1 .and. nn .le. npnt ) then

          jpdfcthphi(ll,mm) = jpdfcthphi(ll,mm) + 1
          cndedcthphi(ll,mm) = cndedcthphi(ll,mm) + pie

          jpdfcthzeta(ll,nn) = jpdfcthzeta(ll,nn) + 1
          cndedcthzeta(ll,nn) = cndedcthzeta(ll,nn) + pie

          jpdfphizeta(mm,nn) = jpdfphizeta(mm,nn) + 1
          cndedphizeta(mm,nn) = cndedphizeta(mm,nn) + pie

        end if

        ll = floor(ctheta/bwcth3d)+1
        mm = floor(phi/bwphi3d)+1
        nn = floor(zeta/bwzeta3d)+1

        if ( ll .ge. 1 .and. ll .le. npnt3d .and. &
             mm .ge. 1 .and. mm .le. npnt3d .and. &
             nn .ge. 1 .and. nn .le. npnt3d ) then

          jpdf3d(ll,mm,nn) = jpdf3d(ll,mm,nn) + 1
          cnded3d(ll,mm,nn) = cnded3d(ll,mm,nn) + pie

        end if
 

      end do
      end do
      end do

      nfile=nfile+1
    end do
  close(20)

  mdissrate = mdissrate * const / nfile

  cnded3d = cnded3d / (jpdf3d + mytiny) 
  cndedcthphi = cndedcthphi / (jpdfcthphi + mytiny)
  cndedcthzeta = cndedcthzeta / (jpdfcthzeta + mytiny)
  cndedphizeta = cndedphizeta / (jpdfphizeta + mytiny)

  jpdf3d = jpdf3d * const / nfile
  jpdfcthphi = jpdfcthphi * const / nfile
  jpdfcthzeta = jpdfcthzeta * const / nfile
  jpdfphizeta = jpdfphizeta * const / nfile

  write(*,*) 'check jpdf3d: ', sum(jpdf3d)
  write(*,*) 'check jpdfcthphi:', sum(jpdfcthphi)
  write(*,*) 'check jpdfcthzeta:', sum(jpdfcthzeta)
  write(*,*) 'check jpdfphizeta:', sum(jpdfphizeta)
  write(*,*) 'check mean dissipation rate: ', mdissrate

  jpdf3d = jpdf3d / bwcth3d / bwphi3d / bwzeta3d
  jpdfcthphi = jpdfcthphi / bwcth / bwphi
  jpdfcthzeta = jpdfcthzeta / bwcth / bwzeta
  jpdfphizeta = jpdfphizeta / bwphi / bwzeta

  open(15,file='cndsgsed-align-st-3d'//str(1:len_trim(str)))
    write(15,'(''zone t = "cnded", i='', i6, '',j='', i6, ''k='', i6, '',f=point'')') & 
      npnt3d, npnt3d, npnt3d

    do kk = 1, npnt3d
    do jj = 1, npnt3d
    do ii = 1, npnt3d
      write(15,'(15E15.5)') (ii-.5)*bwcth3d, (jj-.5)*bwphi3d, (kk-.5)*bwzeta3d, &
                            jpdf3d(ii,jj,kk), cnded3d(ii,jj,kk), &
                            jpdf3d(ii,jj,kk) * cnded3d(ii,jj,kk)
    end do
    end do
    end do
  close(15)

  open(15,file='cndsgsed-align-st-2d-cthphi'//str(1:len_trim(str)))
    write(15,'(''# zone t = "(cth, phi)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do ii = 1, npnt
      do jj = 1, npnt
        write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwphi, jpdfcthphi(ii,jj), &
                              cndedcthphi(ii,jj)
      end do
      write(15,*)
    end do
  close(15)

  open(15,file='cndsgsed-align-st-2d-cthzeta'//str(1:len_trim(str)))
    write(15,'(''# zone t = "(cth, zeta)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do ii = 1, npnt
      do jj = 1, npnt
        write(15,'(15E15.5)') (ii-.5)*bwcth, (jj-.5)*bwzeta, jpdfcthzeta(ii,jj), &
                              cndedcthzeta(ii,jj)
      end do
      write(15,*)
    end do
  close(15)

  open(15,file='cndsgsed-align-st-2d-phizeta'//str(1:len_trim(str)))
    write(15,'(''# zone t = "(phi, zeta)", i='', i6, '',j='', i6, '',f=point'')') npnt, npnt
    do ii = 1, npnt
      do jj = 1, npnt
        write(15,'(15E15.5)') (ii-.5)*bwphi, (jj-.5)*bwzeta, jpdfphizeta(ii,jj), &
                              cndedphizeta(ii,jj)
      end do
      write(15,*)
    end do
  close(15)


  deallocate(s11,s12,s13,s22,s23,s33)
  deallocate(t11,t12,t13,t22,t23,t33)
  deallocate(g,kx,ky,kz)

  call destroyplan3d

  write(*,*) 'finished'

end program sgsedcnd3d 
