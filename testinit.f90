PROGRAM TestInit
  USE mprmtr
  USE mfftwplan
  IMPLICIT NONE

  COMPLEX(SP), ALLOCATABLE, DIMENSION(lx1,ly,lz) :: ux,uy,uz
  REAL(SP), ALLOCATABLE, DIMENSION(lx1,ly,lz) :: kx,ky,kz,k2,atmp
  REAL(SP) :: tmp, pi, S
  INTEGER :: iseed, i, i1, i2

  WRITE(*,*) 'Random seed and file number:'
  READ(*,*) iseed, i
  OPEN(90,FILE='parameter_les.d',STATUS='old')
  READ(90,*) 
  READ(90,*) nx
  READ(90,*) ny
  READ(90,*) nz
  CLOSE(90)
  
  lx=nx/2
  ly=ny
  lz=nz
  lx1=lx+1
  WRITE(*,*) 'lx,ly,lz,lx1', lx,ly,lz,lx1

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz),atmp(lx1,ly,lz))
  
  CALL setfftwplans

  pi=4._SP*ATAN(1._SP)
  CALL wavenumber(kx,ky,kz,k2,atmp,tmp,tmp)
  
  OPEN(15, FILE='inivel.dat',FORM='unformatted')
  READ(15) ux
  READ(15) uy
  READ(15) uz
  CLOSE(15)

  CALL Skewness (ux,S)
  WRITE(*,*) 'Skewness before rumbling:', S
  CALL InitVel(kx,ky,kz,SQRT(k2),ux,uy,uz,pi,iseed)
  CALL Skewness (ux,S)
  WRITE(*,*) 'Skewness after rumbling:', S

  i2=int(i/10)
  i1=mod(i,10)
  OPEN(16,FILE='inivel'//char(i2+48)//char(i1+48)//'.dat', FORM='unformatted')
  WRITE(16) ux
  WRITE(16) uy
  WRITE(16) uz
  CLOSE(16)

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,atmp)
  CALL destroyplans

END PROGRAM TestInit      
