PROGRAM vort_rot
  IMPLICIT NONE

  INTEGER, PARAMETER :: nx=16,ny=16,nz=16
  INTEGER  :: i,ii,iii

  INTEGER      :: VIsDouble, Debug
  CHARACTER(1) :: NULLCHR
  REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: X,Y,Z,enstr

  INTEGER :: TecIni, TecZne, TecDat,TecEnd

  VIsDouble=0
  Debug=1
  NULLCHR=CHAR(0)
  ALLOCATE(X(nx,ny,nz),Y(nx,ny,nz),Z(nx,ny,nz),Enstr(nx,ny,nz))
  
  DO i=1,nx
  X(i,:,:)=i
  END DO
  DO i=1,ny
  Y(:,i,:)=i
  END DO
  DO i=1,nz
  Z(:,:,i)=i
  END DO
  Enstr=SQRT(X*X+Y*Y)


  i = TecIni('test tecio'//NULLCHR, 'X Y Z F'//NULLCHR, 'f.plt'//NULLCHR, &
             '.'//NULLCHR, Debug, VIsDouble)
  i = TecZne('Zone'//NULLCHR,nx,ny,nz,'BLOCK'//NULLCHR,NULLCHR)

  iii=nx*ny*nz
  i = TecDat(iii,X,0)
  i = TecDat(iii,Y,0)
  i = TecDat(iii,Z,0)
  i = TecDat(iii,Enstr,0)

  i = TecEnd()
   
  DEALLOCATE(enstr,X,Y,Z)

END PROGRAM vort_rot      
