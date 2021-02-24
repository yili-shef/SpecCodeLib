PROGRAM vort_rot
  USE mTecIOInterface
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, PARAMETER :: nd=4
  INTEGER  :: nx,ny,nz
  INTEGER  :: lx,ly,lz,lx1

  COMPLEX(SP), ALLOCATABLE, DIMENSION(:,:,:) :: ux,uy,uz
  REAL(SP),    ALLOCATABLE, DIMENSION(:,:,:) :: enstr,wz,kx,ky,kz,k2

  INTEGER  :: i,ii,iii
  REAL(SP) :: ignore_me

  INTEGER      :: VIsDouble, Debug
  CHARACTER(1) :: NULLCHR
  CHARACTER(80):: path
  REAL(SP),    ALLOCATABLE, DIMENSION(:,:,:) :: X,Y,Z

  VIsDouble=0
  Debug=1
  NULLCHR=CHAR(0)
  path='./'

  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli_rot_dns.d',status='unknown')
    READ(90,*) 
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) 
  CLOSE(90)
  CALL fftwplan3d(nx,ny,nz)
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lx1 = lx+1

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz))
  ALLOCATE(k2(lx1,ly,lz),enstr(nx,ny,nz),wz(nx,ny,nz))
  ALLOCATE(X(nx/nd,ny/nd,nz/nd),Y(nx/nd,ny/nd,nz/nd),Z(nx/nd,ny/nd,nz/nd))
  
  DO i=1,nx/nd
  X(i,:,:)=i
  END DO
  DO i=1,ny/nd
  Y(:,i,:)=i
  END DO
  DO i=1,nz/nd
  Z(:,:,i)=i
  END DO

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  OPEN(16,FILE=path(1:LEN_TRIM(path))//'out/vel0051.dat',FORM='unformatted')
    READ(16) ux
    READ(16) uy
    READ(16) uz
  CLOSE(16)
  
  CALL enstrophy(ux,uy,uz,kx,ky,kz,enstr,lx1,ly,lz,lx,nx,ny,nz)

  uz=eye*(kx*uy-ky*ux)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)
  wz(1:nx:2,:,:)=REAL(uz(1:lx,:,:))
  wz(2:nx:2,:,:)=AIMAG(uz(1:lx,:,:))

  i = TecIni('3D Enstrophy field'//NULLCHR, 'X Y Z Enstr'//NULLCHR, 'enstr.plt'//NULLCHR, &
             '.'//NULLCHR, Debug, VIsDouble)
  i = TecZne('Zone'//NULLCHR,nx/nd,ny/nd,nz/nd,'BLOCK'//NULLCHR,NULLCHR)

  iii=nx*ny*nz/(nd*nd*nd)
  i = TecDat(iii,X,0)
  i = TecDat(iii,Y,0)
  i = TecDat(iii,Z,0)
  i = TecDat(iii,Enstr(1:nx/nd,1:ny/nd,1:nz/nd),0)

  i = TecEnd()
   
  i = TecIni('3D field of wz'//NULLCHR, 'X Y Z wz'//NULLCHR, 'wz.plt'//NULLCHR, &
             '.'//NULLCHR, Debug, VIsDouble)
  i = TecZne('Zone'//NULLCHR,nx/nd,ny/nd,nz/nd,'BLOCK'//NULLCHR,NULLCHR)

  iii=nx*ny*nz/(nd*nd*nd)
  i = TecDat(iii,X,0)
  i = TecDat(iii,Y,0)
  i = TecDat(iii,Z,0)
  i = TecDat(iii,wz(1:nx/nd,1:ny/nd,1:nz/nd),0)

  i = TecEnd()

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,wz,enstr,X,Y,Z)

  CALL destroyplan3d

END PROGRAM vort_rot      
