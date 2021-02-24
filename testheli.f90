PROGRAM testheli
  USE mprmtr
  USE mfftwplan
  IMPLICIT NONE

  COMPLEX(SP), ALLOCATABLE, DIMENSION(lx1,ly,lz) :: ux,uy,uz,wx,wy,wz
  COMPLEX,     ALLOCATABLE, DIMENSION(lx1,ly,lz) :: uxs,uys,uzs
  REAL(SP),    ALLOCATABLE, DIMENSION(lx1,ly,lz) :: kx,ky,kz,k2,atmp
  REAL(SP),    ALLOCATABLE, DIMENSION(lx1,ly,lz) :: ModRe, ModIm, handed, phi
  REAL(SP),    ALLOCATABLE, DIMENSION(:)     :: hk, hk1

  REAL(SP) :: tmp, S
  INTEGER  :: nk, ni, nnn, i, ii, iii, nksum

  Pi=4._SP*ATAN(1._SP)
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

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz),wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz),atmp(lx1,ly,lz))
  ALLOCATE(ModRe(lx1,ly,lz),ModIm(lx1,ly,lz),handed(lx1,ly,lz),phi(lx1,ly,lz))
  ALLOCATE(hk(lx),hk1(lx))
  ALLOCATE(uxs(lx1,ly,lz),uys(lx1,ly,lz),uzs(lx1,ly,lz))

  CALL setfftwplans
  CALL wavenumber(kx,ky,kz,k2,atmp,tmp,tmp)

  OPEN(15, FILE='inivel.dat',FORM='unformatted')
    READ(15) uxs
    READ(15) uys
    READ(15) uzs
  CLOSE(15)
  ux=CMPLX(REAL(uxs,SP),REAL(AIMAG(uxs),SP))
  uy=CMPLX(REAL(uys,SP),REAL(AIMAG(uys),SP))
  uz=CMPLX(REAL(uzs,SP),REAL(AIMAG(uzs),SP))

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  
  atmp = (ux*CONJG(wx)+uy*CONJG(wy)+uz*CONJG(wz)) &
         +(CONJG(ux)*wx+CONJG(uy)*wy+CONJG(uz)*wz)
  atmp(1,:,:) = 0.5_SP*atmp(1,:,:)
  DO i=1,lx
    hk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
  END DO
  
  atmp=kx*REAL(uy,SP)*AIMAG(uz) &
      -kx*REAL(uz,SP)*AIMAG(uy) &
      +ky*REAL(uz,SP)*AIMAG(ux) &
      -ky*REAL(ux,SP)*AIMAG(uz) &
      +kz*REAL(ux,SP)*AIMAG(uy) &
      -kz*REAL(uy,SP)*AIMAG(ux)
  handed=SIGN(1._SP,atmp)
  ModRe=SQRT(REAL(ux,SP)**2+REAL(uy,SP)**2+REAL(uz,SP)**2)
  ModIm=SQRT(AIMAG(ux)**2+AIMAG(uy)**2+AIMAG(uz)**2)
  atmp=(REAL(ux,SP)*AIMAG(ux)+REAL(uy,SP)*AIMAG(uy) &
      +REAL(uz,SP)*AIMAG(uz))/(ModRe*ModIm+smallest)
  phi=handed*ACOS(atmp)+Pi*(1._SP-handed)
  atmp=2._SP*2._SP*SQRT(k2)*ModRe*ModIm*SIN(phi)
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx
   hk1(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
  END DO

  OPEN(16, FILE='testh.dat')
  DO i=1,lx
  WRITE(*,*) i, hk(i), hk1(i)
  WRITE(16,*) i, hk(i), hk1(i)
  END DO
  CLOSE(16)

  
  DEALLOCATE(ux,uy,uz,wx,wy,wz,kx,ky,kz,atmp,modre,modim,handed,phi,hk,hk1,uxs,uys,uzs)
  
  CALL destroyplans 
END PROGRAM testheli
