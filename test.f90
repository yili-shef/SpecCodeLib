PROGRAM transform
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, PARAMETER :: nx=64,ny=64,nz=64
  INTEGER, PARAMETER :: lx=nx/2,ly=ny,lz=nz,lx1=lx+1
  
  REAL(SP), DIMENSION(nx,ny,nz) :: uxr,uyr,uzr,pr
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: ux,uy,uz,p

  INTEGER, PARAMETER :: if1=1,if2=1707

  INTEGER :: ifile,id1,id2,id3,id4
  REAL(SP) :: ignore_me
  
  CALL fftwplan3d(nx,ny,nz)

  DO ifile = if1, if2
  
  id1=MOD(ifile,10)+48
  id2=MOD(INT(ifile/10),10)+48
  id3=MOD(INT(ifile/100),10)+48
  id4=MOD(INT(ifile/1000),10)+48
  OPEN(16,FILE='./out/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='unformatted')
    READ(16) ux
    READ(16) uy
    READ(16) uz
  CLOSE(16)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ux,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)
  
  uxr(1:nx:2,:,:)= REAL(ux(1:lx,:,:))
  uxr(2:nx:2,:,:)=AIMAG(ux(1:lx,:,:))

  uyr(1:nx:2,:,:)= REAL(uy(1:lx,:,:))
  uyr(2:nx:2,:,:)=AIMAG(uy(1:lx,:,:))

  uzr(1:nx:2,:,:)= REAL(uz(1:lx,:,:))
  uzr(2:nx:2,:,:)=AIMAG(uz(1:lx,:,:))

  OPEN(16,FILE='./data/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='binary',ACCESS='sequential')
    WRITE(16) uxr
    WRITE(16) uyr
    WRITE(16) uzr
  CLOSE(16)

  OPEN(16,FILE='./out/p'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='unformatted')
    READ(16) p
  CLOSE(16)


  CALL rfftwnd_f77_one_complex_to_real(C2R3D,p,ignore_me)
  pr(1:nx:2,:,:)= REAL(p(1:lx,:,:))
  pr(2:nx:2,:,:)=AIMAG(p(1:lx,:,:))
  OPEN(16,FILE='./data/p'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat',FORM='binary',ACCESS='sequential')
    WRITE(16) pr
  CLOSE(16)

  END DO

  CALL destroyplan3d

END PROGRAM transform      
