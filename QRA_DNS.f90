PROGRAM materialdeform
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nx=128,ny=128,nz=128, nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=211
  REAL(SP), PARAMETER :: dt=0.017

  INTEGER,  PARAMETER :: lx=nx/2,ly=ny,lz=nz,lx1=lx+1

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: ux,uy,uz
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: A11,A12,A13
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: A21,A22,A23
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: A31,A32,A33

  REAL(SP),    DIMENSION(lx1,ly,lz) :: kx,ky,kz,k2
  
  REAL(SP),    DIMENSION(nx,ny,nz)  :: uxr,uyr,uzr
  REAL(SP),    DIMENSION(nx,ny,nz)  :: rA11,rA12,rA13
  REAL(SP),    DIMENSION(nx,ny,nz)  :: rA21,rA22,rA23
  REAL(SP),    DIMENSIOn(nx,ny,nz)  :: rA31,rA32,rA33
  REAL(SP),    DIMENSION(nprtcl)    :: xp,yp,zp
  REAL(SP),    DIMENSION(nprtcl)    :: xps,yps,zps
  REAL(SP),    DIMENSIOn(nprtcl)    :: Q,R
  REAL(SP),    DIMENSION(3,3,nprtcl):: Ap

  REAL(SP) :: y(3),dxyz(3),bg(6,3)
  INTEGER  :: lhnode(3),idnty(3,3)
  INTEGER  :: i,ifile,id1,id2,id3,id4,ip,idum
  REAL(SP) :: uxp,uyp,uzp
  REAL(SP) :: ignore_me,ran1
  

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
  dxyz(1)=2.*pi/REAL(nx)
  dxyz(2)=2.*pi/REAL(ny)
  dxyz(3)=2.*pi/REAL(nz)

  CALL fftwplan3d(nx,ny,nz)
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  idum=-7
  DO i=1,nprtcl
  xp(i)=2.*pi*ran1(idum)
  yp(i)=2.*pi*ran1(idum)
  zp(i)=2.*pi*ran1(idum)
  END DO

  ifile=if1
  WRITE(*,*) ifile
  id1=MOD(ifile,10)+48
  id2=MOD(INT(ifile/10),10)+48
  id3=MOD(INT(ifile/100),10)+48
  id4=MOD(INT(ifile/1000),10)+48
  OPEN(16,FILE='./out/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='unformatted')
    READ(16) ux
    READ(16) uy
    READ(16) uz
  CLOSE(16)

  A11=eye*kx*ux; A12=eye*kx*uy; A13=eye*kx*uz
  A21=eye*ky*ux; A22=eye*ky*uy; A23=eye*ky*uz
  A31=eye*kz*ux; A32=eye*kz*uy; A33=eye*kz*uz

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ux,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A21,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A31,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A32,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,A33,ignore_me)
  
  uxr(1:nx:2,:,:)=REAL(ux(1:lx,:,:)); uxr(2:nx:2,:,:)=AIMAG(ux(1:lx,:,:))
  uyr(1:nx:2,:,:)=REAL(uy(1:lx,:,:)); uyr(2:nx:2,:,:)=AIMAG(uy(1:lx,:,:))
  uzr(1:nx:2,:,:)=REAL(uz(1:lx,:,:)); uzr(2:nx:2,:,:)=AIMAG(uz(1:lx,:,:))

  rA11(1:nx:2,:,:)=REAL(A11(1:lx,:,:)); rA11(2:nx:2,:,:)=AIMAG(A11(1:lx,:,:))
  rA12(1:nx:2,:,:)=REAL(A12(1:lx,:,:)); rA12(2:nx:2,:,:)=AIMAG(A12(1:lx,:,:))
  rA13(1:nx:2,:,:)=REAL(A13(1:lx,:,:)); rA13(2:nx:2,:,:)=AIMAG(A13(1:lx,:,:))
  rA21(1:nx:2,:,:)=REAL(A21(1:lx,:,:)); rA21(2:nx:2,:,:)=AIMAG(A21(1:lx,:,:))
  rA22(1:nx:2,:,:)=REAL(A22(1:lx,:,:)); rA22(2:nx:2,:,:)=AIMAG(A22(1:lx,:,:))
  rA23(1:nx:2,:,:)=REAL(A23(1:lx,:,:)); rA23(2:nx:2,:,:)=AIMAG(A23(1:lx,:,:))
  rA31(1:nx:2,:,:)=REAL(A31(1:lx,:,:)); rA31(2:nx:2,:,:)=AIMAG(A31(1:lx,:,:))
  rA32(1:nx:2,:,:)=REAL(A32(1:lx,:,:)); rA32(2:nx:2,:,:)=AIMAG(A32(1:lx,:,:))
  rA33(1:nx:2,:,:)=REAL(A33(1:lx,:,:)); rA33(2:nx:2,:,:)=AIMAG(A33(1:lx,:,:))

  DO ifile = if1+1,if2
    DO ip=1, nprtcl
      y(1)=xp(ip)
      y(2)=yp(ip)
      y(3)=zp(ip)

      CALL pre_interp(y,dxyz,bg,lhnode)
      CALL value(uxr,uxp,lhnode,bg,nx,ny,nz)
      CALL value(uyr,uyp,lhnode,bg,nx,ny,nz)
      CALL value(uzr,uzp,lhnode,bg,nx,ny,nz)

      CALL value(rA11,Ap(1,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA12,Ap(1,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA13,Ap(1,3,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA21,Ap(2,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA22,Ap(2,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA23,Ap(2,3,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA31,Ap(3,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA32,Ap(3,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA33,Ap(3,3,ip),lhnode,bg,nx,ny,nz)

      xps(ip)=xp(ip)+dt*uxp
      yps(ip)=yp(ip)+dt*uyp
      zps(ip)=zp(ip)+dt*uzp
      
      Q(ip)=-.5*SUM(MATMUL(Ap(:,:,ip),Ap(:,:,ip)), MASK=idnty .EQ. 1)
      R(ip)=-(1./3.)*SUM(MATMUL(MATMUL(Ap(:,:,ip),Ap(:,:,ip)),Ap(:,:,ip)), MASK=idnty .EQ. 1)
    END DO
    WRITE(20) Q,R,Ap

    WRITE(*,*) ifile
    id1=MOD(ifile,10)+48
    id2=MOD(INT((ifile)/10),10)+48
    id3=MOD(INT((ifile)/100),10)+48
    id4=MOD(INT((ifile)/1000),10)+48
    OPEN(16,FILE='./out/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='unformatted')
      READ(16) ux
      READ(16) uy
      READ(16) uz
    CLOSE(16)

    A11=eye*kx*ux; A12=eye*kx*uy; A13=eye*kx*uz
    A21=eye*ky*ux; A22=eye*ky*uy; A23=eye*ky*uz
    A31=eye*kz*ux; A32=eye*kz*uy; A33=eye*kz*uz

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,ux,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uy,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A21,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A23,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A31,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A32,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,A33,ignore_me)
    
    uxr(1:nx:2,:,:)=REAL(ux(1:lx,:,:)); uxr(2:nx:2,:,:)=AIMAG(ux(1:lx,:,:))
    uyr(1:nx:2,:,:)=REAL(uy(1:lx,:,:)); uyr(2:nx:2,:,:)=AIMAG(uy(1:lx,:,:))
    uzr(1:nx:2,:,:)=REAL(uz(1:lx,:,:)); uzr(2:nx:2,:,:)=AIMAG(uz(1:lx,:,:))

    rA11(1:nx:2,:,:)=REAL(A11(1:lx,:,:)); rA11(2:nx:2,:,:)=AIMAG(A11(1:lx,:,:))
    rA12(1:nx:2,:,:)=REAL(A12(1:lx,:,:)); rA12(2:nx:2,:,:)=AIMAG(A12(1:lx,:,:))
    rA13(1:nx:2,:,:)=REAL(A13(1:lx,:,:)); rA13(2:nx:2,:,:)=AIMAG(A13(1:lx,:,:))
    rA21(1:nx:2,:,:)=REAL(A21(1:lx,:,:)); rA21(2:nx:2,:,:)=AIMAG(A21(1:lx,:,:))
    rA22(1:nx:2,:,:)=REAL(A22(1:lx,:,:)); rA22(2:nx:2,:,:)=AIMAG(A22(1:lx,:,:))
    rA23(1:nx:2,:,:)=REAL(A23(1:lx,:,:)); rA23(2:nx:2,:,:)=AIMAG(A23(1:lx,:,:))
    rA31(1:nx:2,:,:)=REAL(A31(1:lx,:,:)); rA31(2:nx:2,:,:)=AIMAG(A31(1:lx,:,:))
    rA32(1:nx:2,:,:)=REAL(A32(1:lx,:,:)); rA32(2:nx:2,:,:)=AIMAG(A32(1:lx,:,:))
    rA33(1:nx:2,:,:)=REAL(A33(1:lx,:,:)); rA33(2:nx:2,:,:)=AIMAG(A33(1:lx,:,:))

    DO ip=1,nprtcl
      y(1)=xps(ip)
      y(2)=yps(ip)
      y(3)=zps(ip)

      CALL pre_interp(y,dxyz,bg,lhnode)
      CALL value(uxr,uxp,lhnode,bg,nx,ny,nz)
      CALL value(uyr,uyp,lhnode,bg,nx,ny,nz)
      CALL value(uzr,uzp,lhnode,bg,nx,ny,nz)

      CALL value(rA11,Ap(1,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA12,Ap(1,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA13,Ap(1,3,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA21,Ap(2,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA22,Ap(2,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA23,Ap(2,3,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA31,Ap(3,1,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA32,Ap(3,2,ip),lhnode,bg,nx,ny,nz)
      CALL value(rA33,Ap(3,3,ip),lhnode,bg,nx,ny,nz)

      xp(ip)=.5*(xp(ip)+xps(ip))+.5*dt*uxp
      yp(ip)=.5*(yp(ip)+yps(ip))+.5*dt*uyp
      zp(ip)=.5*(zp(ip)+zps(ip))+.5*dt*uzp

    END DO
  END DO
  DO ip=1, nprtcl
    y(1)=xp(ip)
    y(2)=yp(ip)
    y(3)=zp(ip)

    CALL pre_interp(y,dxyz,bg,lhnode)
    CALL value(uxr,uxp,lhnode,bg,nx,ny,nz)
    CALL value(uyr,uyp,lhnode,bg,nx,ny,nz)
    CALL value(uzr,uzp,lhnode,bg,nx,ny,nz)

    CALL value(rA11,Ap(1,1,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA12,Ap(1,2,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA13,Ap(1,3,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA21,Ap(2,1,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA22,Ap(2,2,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA23,Ap(2,3,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA31,Ap(3,1,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA32,Ap(3,2,ip),lhnode,bg,nx,ny,nz)
    CALL value(rA33,Ap(3,3,ip),lhnode,bg,nx,ny,nz)

    Q(ip)=-.5*SUM(MATMUL(Ap(:,:,ip),Ap(:,:,ip)), MASK=idnty .EQ. 1)
    R(ip)=-(1./3.)*SUM(MATMUL(MATMUL(Ap(:,:,ip),Ap(:,:,ip)),Ap(:,:,ip)), MASK=idnty .EQ. 1)
  END DO
  WRITE(20) Q,R,Ap

  CLOSE(20)

  CALL destroyplan3d
  WRITE(*,*) 'finished'
  STOP

END PROGRAM materialdeform      
