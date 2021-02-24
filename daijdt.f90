PROGRAM daijdt
  USE mTecIOInterface
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx,ny,nz
  INTEGER  :: lx,ly,lz,lx1
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz,wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: t11,t12,t13,t22,t23,t33
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: A11,A12,A13,A21,A22,A23,A31,A32,A33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,tmp,G
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: Q,duxdt,duydt,duzdt
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: del_ux,del_uy,del_uz
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: del_duxdt,del_duydt,del_duzdt
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: dulldt_dns,dutldt_dns,dulldt,dutldt

  INTEGER :: VIsDouble,Debug
  CHARACTER(1) :: NULLCHR

  CHARACTER(80) :: path
  INTEGER  :: ndel,i,iii
  REAL(SP) :: aa,eps,eta,dx,dx0,ignore_me
  REAL(SP) :: rnu,delta_c
  REAL(SP), PARAMETER :: alpha=40._SP



  VIsDouble=0
  Debug=1
  NULLCHR=CHAR(0)

  path='./'
  sufix='_dns.data'
  OPEN(90,FILE=path(1:LEN_TRIM(path))//'parameter_dns.d',STATUS='unknown')
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) rnu
  CLOSE(90)
  CALL fftwplan3d(nx,ny,nz)
  
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lx1 = lx+1 

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  ALLOCATE(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  ALLOCATE(A11(lx1,ly,lz),A12(lx1,ly,lz),A13(lx1,ly,lz))
  ALLOCATE(A21(lx1,ly,lz),A22(lx1,ly,lz),A23(lx1,ly,lz))
  ALLOCATE(A31(lx1,ly,lz),A32(lx1,ly,lz),A33(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(tmp(lx1,ly,lz),G(lx1,ly,lz),Q(nx,ny,nz))
  ALLOCATE(duxdt(nx,ny,nz),duydt(nx,ny,nz),duzdt(nx,ny,nz))
  ALLOCATE(del_ux(nx,ny,nz),del_uy(nx,ny,nz),del_uz(nx,ny,nz))
  ALLOCATE(del_duxdt(nx,ny,nz),del_duydt(nx,ny,nz),del_duzdt(nx,ny,nz))
  ALLOCATE(dulldt_dns(nx,ny,nz),dutldt_dns(nx,ny,nz),dulldt(nx,ny,nz),dutldt(nx,ny,nz))


  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  OPEN(16,FILE='vel0012.dat',FORM='unformatted')
    READ(16) ux
    READ(16) uy
    READ(16) uz
  CLOSE(16)

  tmp=2._SP*(.5_SP*(ux*CONJG(ux)+uy*CONJG(uy)+uz*CONJG(uz)))
  tmp(1,:,:)=.5_SP*tmp(1,:,:) 
  eps=0._SP
  DO i=1,lx1-1
     aa=SUM(tmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
     eps=eps+2._SP*rnu*REAL(i,SP)*REAL(i,SP)*aa
  END DO

  eta=(rnu*rnu*rnu/eps)**0.25_SP
  delta_c=alpha*eta

  dx=delta_c
  ndel=INT(dx*nx/(2._SP*pi))
  dx=ndel*(2._SP*pi/REAL(nx,SP))
  dx0=dx

  G=EXP(-k2*delta_c**2/24._SP)

  CALL sgsstress(ux,uy,uz,t11,t22,t33,t12,t13,t23,G,lx1,ly,lz,lx,nx,ny,nz)
  ux=G*ux
  uy=G*uy
  uz=G*uz

  wx = eye * (ky*uz - kz*uy)
  wy = eye * (kz*ux - kx*uz)
  wz = eye * (kx*uy - ky*ux)

  CALL convec_dns(ux,uy,uz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
  wx=wx-eye*(kx*t11+ky*t12+kz*t13)-rnu*k2*ux
  wy=wy-eye*(kx*t12+ky*t22+kz*t23)-rnu*k2*uy
  wz=wz-eye*(kx*t13+ky*t23+kz*t33)-rnu*k2*uz
  
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wx,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wz,ignore_me)

  A11=eye*kx*ux; A12=eye*kx*uy; A13=eye*kx*uz
  A21=eye*ky*ux; A22=eye*ky*uy; A23=eye*ky*uz
  A31=eye*kz*ux; A32=eye*kz*uy; A33=eye*kz*uz

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A21,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A31,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A32,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A33,ignore_me)
  
  Q(1:nx:2,:,:)=-(1._SP/2._SP)*(REAL(A11(1:lx,:,:))*REAL(A11(1:lx,:,:))+REAL(A22(1:lx,:,:))*REAL(A22(1:lx,:,:)) &
                               +REAL(A33(1:lx,:,:))*REAL(A33(1:lx,:,:))+REAL(A12(1:lx,:,:))*REAL(A21(1:lx,:,:)) &
                               +REAL(A13(1:lx,:,:))*REAL(A31(1:lx,:,:))+REAL(A21(1:lx,:,:))*REAL(A12(1:lx,:,:)) &
                               +REAL(A23(1:lx,:,:))*REAL(A32(1:lx,:,:))+REAL(A31(1:lx,:,:))*REAL(A13(1:lx,:,:)) &
                               +REAL(A32(1:lx,:,:))*REAL(A23(1:lx,:,:)))

  Q(2:nx:2,:,:)=-(1._SP/2._SP)*(AIMAG(A11(1:lx,:,:))*AIMAG(A11(1:lx,:,:))+AIMAG(A22(1:lx,:,:))*AIMAG(A22(1:lx,:,:)) &
                               +AIMAG(A33(1:lx,:,:))*AIMAG(A33(1:lx,:,:))+AIMAG(A12(1:lx,:,:))*AIMAG(A21(1:lx,:,:)) &
                               +AIMAG(A13(1:lx,:,:))*AIMAG(A31(1:lx,:,:))+AIMAG(A21(1:lx,:,:))*AIMAG(A12(1:lx,:,:)) &
                               +AIMAG(A23(1:lx,:,:))*AIMAG(A32(1:lx,:,:))+AIMAG(A31(1:lx,:,:))*AIMAG(A13(1:lx,:,:)) &
                               +AIMAG(A32(1:lx,:,:))*AIMAG(A23(1:lx,:,:)))

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ux,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)

  duxdt(1:nx:2,:,:)=REAL(wx(1:lx,:,:))+REAL(ux(1:lx,:,:))*REAL(A11(1:lx,:,:)) &
                   +REAL(uy(1:lx,:,:))*REAL(A21(1:lx,:,:))+REAL(uz(1:lx,:,:))*REAL(A31(1:lx,:,:))
  duxdt(2:nx:2,:,:)=AIMAG(wx(1:lx,:,:))+AIMAG(ux(1:lx,:,:))*AIMAG(A11(1:lx,:,:)) &
                   +AIMAG(uy(1:lx,:,:))*AIMAG(A21(1:lx,:,:))+AIMAG(uz(1:lx,:,:))*AIMAG(A31(1:lx,:,:))
  duydt(1:nx:2,:,:)=REAL(wy(1:lx,:,:))+REAL(ux(1:lx,:,:))*REAL(A12(1:lx,:,:)) &
                   +REAL(uy(1:lx,:,:))*REAL(A22(1:lx,:,:))+REAL(uz(1:lx,:,:))*REAL(A32(1:lx,:,:))
  duydt(2:nx:2,:,:)=AIMAG(wy(1:lx,:,:))+AIMAG(ux(1:lx,:,:))*AIMAG(A12(1:lx,:,:)) &
                   +AIMAG(uy(1:lx,:,:))*AIMAG(A22(1:lx,:,:))+AIMAG(uz(1:lx,:,:))*AIMAG(A32(1:lx,:,:))
  duzdt(1:nx:2,:,:)=REAL(wz(1:lx,:,:))+REAL(ux(1:lx,:,:))*REAL(A13(1:lx,:,:)) &
                   +REAL(uy(1:lx,:,:))*REAL(A23(1:lx,:,:))+REAL(uz(1:lx,:,:))*REAL(A33(1:lx,:,:))
  duzdt(2:nx:2,:,:)=AIMAG(wz(1:lx,:,:))+AIMAG(ux(1:lx,:,:))*AIMAG(A13(1:lx,:,:)) &
                   +AIMAG(uy(1:lx,:,:))*AIMAG(A23(1:lx,:,:))+AIMAG(uz(1:lx,:,:))*AIMAG(A33(1:lx,:,:))


  wx=CSHIFT(ux,ndel,1)-ux
  wy=CSHIFT(uy,ndel,1)-uy
  wz=CSHIFT(uz,ndel,1)-uz

  del_ux(1:nx:2,:,:)=REAL(wx(1:lx,:,:))
  del_ux(2:nx:2,:,:)=AIMAG(wx(1:lx,:,:))
  del_uy(1:nx:2,:,:)=REAL(wy(1:lx,:,:))
  del_uy(2:nx:2,:,:)=AIMAG(wy(1:lx,:,:))
  del_uz(1:nx:2,:,:)=REAL(wz(1:lx,:,:))
  del_uz(2:nx:2,:,:)=AIMAG(wz(1:lx,:,:))
  del_ut=SQRT(del_uy*del_uy+del_uz*del_uz)

  del_duxdt=CSHIFT(duxdt,ndel,1)-duxdt
  del_duydt=CSHIFT(duydt,ndel,1)-duydt
  del_duzdt=CSHIFT(duzdt,ndel,1)-duzdt

  dulldt_dns=del_duxdt*dx0/dx+(del_ut*del_ut-del_ux*del_ux)*(dx0/dx)**2/dx0
  dutldt_dns=(del_duydt*del_uy+del_duzdt*del_uz)/(del_ut+smallest)*(dx0/dx)-2._SP*del_ux*del_ut*(dx0/dx)**2/dx0

  dulldt=-(2._SP/3._SP)*Q*dx0+(del_ut*del_ut-del_ux*del_ux)*(dx0/dx)**2/dx0
  dutldt=-2._SP*del_ux*del_ut*(dx0/dx)**2/dx0

  i = TecIni('Velocity increments'//NULLCHR, 'dulldt_dns dulldt dutldt_dns dutldt'//NULLCHR, 'du_x.plt'//NULLCHR, &
             '.'//NULLCHR, Debug, VIsDouble)
  i = TecZne('Zone'//NULLCHR,nx,ny,nz,'BLOCK'//NULLCHR,NULLCHR)

  iii=nx*ny*nz
  i = TecDat(iii,dulldt_dns,0)
  i = TecDat(iii,dulldt,0)
  i = TecDat(iii,dutldt_dns,0)
  i = TecDat(iii,dutldt,0)

  i = TecEnd()

  
  DEALLOCATE(ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33)
  DEALLOCATE(A11,A12,A13,A21,A22,A23,A31,A32,A33)
  DEALLOCATE(kx,ky,kz,k2,tmp,G,Q,duxdt,duydt,duzdt)
  DEALLOCATE(del_ux,del_uy,del_uz,del_duxdt,del_duydt,del_duzdt)
  DEALLOCATE(dulldt_dns,dutldt_dns,dulldt,dutldt)
  CALL destroyplan3d
  
END PROGRAM daijdt      
