PROGRAM materialdeform
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nx=128,ny=128,nz=128
  INTEGER,  PARAMETER :: if1=1,if2=211
  INTEGER,  PARAMETER :: lx=nx/2,ly=ny,lz=nz,lx1=lx+1

  REAL(DP), PARAMETER :: rnu=0.0032, eps=0.1

  COMPLEX(DP), DIMENSION(lx1,ly,lz) :: ux,uy,uz
  COMPLEX(DP), DIMENSION(lx1,ly,lz) :: S11,S12,S13,S22,S23,S33
  REAL(DP),    DIMENSION(lx1,ly,lz) :: G,k2,kx,ky,kz
  REAL(DP),    DIMENSION(nx, ny,nz) :: S

  REAL(DP) :: eta,delta,ignore_me,ms,mms
  INTEGER  :: ifile,id1,id2,id3,id4

  eta=(rnu**3/eps)**.25
  delta=20.*eta

  CALL fftwplan3d(nx,ny,nz)
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  mms=0.
  DO ifile=if1,if2

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

  S11=eye*kx*ux; S22=eye*ky*uy; S33=eye*kz*uz
  S12=.5*eye*(ky*ux+kx*uy)
  S13=.5*eye*(kz*ux+kx*uz)
  S23=.5*eye*(ky*uz+kz*uy)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,S33,ignore_me)

  S(1:nx:2,:,:)=REAL(S11(1:lx,:,:))**2+REAL(S22(1:lx,:,:))**2+REAL(S33(1:lx,:,:))**2    &
               +2.*(REAL(S12(1:lx,:,:))**2+REAL(S13(1:lx,:,:))**2+REAL(S23(1:lx,:,:))**2)
  
  S(2:nx:2,:,:)=AIMAG(S11(1:lx,:,:))**2+AIMAG(S22(1:lx,:,:))**2+AIMAG(S33(1:lx,:,:))**2 &
               +2.*(AIMAG(S12(1:lx,:,:))**2+AIMAG(S13(1:lx,:,:))**2+AIMAG(S23(1:lx,:,:))**2)
  
  S=SQRT(2.*S)
  ms=SUM(S)/(nx*ny*nz)
  ms=1./ms
  WRITE(15,*) ms
  mms=mms+ms
  END DO
  mms=mms/(if2-if1+1)
  WRITE(15,*) 'mean', mms

END PROGRAM materialdeform
