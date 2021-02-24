PROGRAM InitVel
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,lx,ly,lz,lx1
  COMPLEX(SP), ALLOCATABLE, DIMENSION(:,:,:) :: ux,uy,uz,wx,wy,wz
  REAL(SP),    ALLOCATABLE, DIMENSION(:,:,:) :: kx,ky,kz,k2,atmp
  REAL(SP),    ALLOCATABLE, DIMENSION(:)     :: hkinit,ekinit,htmp,etmp
  
  COMPLEX(SP), ALLOCATABLE, DIMENSION(:)     :: uxsh, uysh, uzsh
  REAL(SP),    ALLOCATABLE, DIMENSION(:)     :: kxsh, kysh, kzsh
  REAL(SP),    ALLOCATABLE, DIMENSION(:)     :: modk
  REAL(SP),    ALLOCATABLE, DIMENSION(:)     :: handed

  REAL(SP) :: tmp,tmp1, S
  INTEGER  :: nk, ni, nnn, i, ii, iii, nksum

  OPEN(90,FILE='parameter_heli.d',STATUS='old')
    READ(90,*) 
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
  ALLOCATE(hkinit(lx),ekinit(lx),htmp(lx),etmp(lx))

  CALL fftwplan3d(nx,ny,nz)
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  OPEN(15, FILE='inivel.dat',FORM='unformatted')
    READ(15) ux
    READ(15) uy
    READ(15) uz
  CLOSE(15)

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  
  atmp = (ux*CONJG(wx)+uy*CONJG(wy)+uz*CONJG(wz)) &
         +(CONJG(ux)*wx+CONJG(uy)*wy+CONJG(uz)*wz)
  atmp(1,:,:) = 0.5_SP*atmp(1,:,:)
  DO i=1,lx
    hkinit(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
  END DO
  atmp = ux*CONJG(ux) + uy*CONJG(uy) + uz*CONJG(uz)
  atmp(1,:,:)=0.5_SP*atmp(1,:,:)
  DO i=1,lx
    ekinit(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
  END DO

  OPEN(16,FILE='initspec.dat')
  DO i=1,lx
    WRITE(16,*) i,ekinit(i),hkinit(i)
  END DO
  CLOSE(16)
  
  OPEN(15, FILE='inivel1.dat',FORM='unformatted')
    READ(15) ux
    READ(15) uy
    READ(15) uz
  CLOSE(15)

  CALL Skewness (ux,kx,S)
  WRITE(*,*) 'Skewness before matching:', S
  
  atmp = ux*CONJG(ux) + uy*CONJG(uy) + uz*CONJG(uz)
  atmp(1,:,:)=0.5_SP*atmp(1,:,:)
  DO i=1,lx
    tmp=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
    WHERE(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP)
            ux=ux*SQRT(ekinit(i)/(tmp+smallest))
            uy=uy*SQRT(ekinit(i)/(tmp+smallest))
            uz=uz*SQRT(ekinit(i)/(tmp+smallest))
    ENDWHERE     
  END DO
  
  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  
  atmp = (ux*CONJG(wx)+uy*CONJG(wy)+uz*CONJG(wz)) &
         +(CONJG(ux)*wx+CONJG(uy)*wy+CONJG(uz)*wz)
  atmp(1,:,:) = 0.5_SP*atmp(1,:,:)
  DO nnn=1, lx
    nk=COUNT(ABS(SQRT(k2)-nnn-0.5_SP*oneless).LT.0.5_SP)
    WRITE(*,*) 'nnn, nk', nnn, nk
    ALLOCATE(uxsh(nk),uysh(nk),uzsh(nk),kxsh(nk),kysh(nk),kzsh(nk))
    ALLOCATE(handed(nk),modk(nk))
    
    ni=1
    DO iii=1,lz
      DO ii=1,ly
        DO i=1,lx1
          IF (ABS(SQRT(k2(i,ii,iii))-nnn-0.5_SP*oneless).LT.0.5_SP) THEN 
                  kxsh(ni)=kx(i,ii,iii)
                  kysh(ni)=ky(i,ii,iii)
                  kzsh(ni)=kz(i,ii,iii)
                  uxsh(ni)=ux(i,ii,iii)
                  uysh(ni)=uy(i,ii,iii)
                  uzsh(ni)=uz(i,ii,iii)
                  tmp=kxsh(ni)*REAL(uysh(ni),SP)*AIMAG(uzsh(ni)) &
                     -kxsh(ni)*REAL(uzsh(ni),SP)*AIMAG(uysh(ni)) &
                     +kysh(ni)*REAL(uzsh(ni),SP)*AIMAG(uxsh(ni)) &
                     -kysh(ni)*REAL(uxsh(ni),SP)*AIMAG(uzsh(ni)) &
                     +kzsh(ni)*REAL(uxsh(ni),SP)*AIMAG(uysh(ni)) &
                     -kzsh(ni)*REAL(uysh(ni),SP)*AIMAG(uxsh(ni))
                  handed(ni)=SIGN(1._SP,tmp)
                  modk(ni)=SQRT(k2(i,ii,iii))
                  ni=ni+1
          END IF
        END DO
      END DO
    END DO
    CALL MatchHeli (uxsh,uysh,uzsh,kxsh,kysh,kzsh,handed,modk,hkinit(nnn),nk)
    ni=1
    DO iii=1,lz
      DO ii=1,ly
        DO i=1,lx1
          IF (ABS(SQRT(k2(i,ii,iii))-nnn-0.5_SP*oneless).LT.0.5_SP) THEN 
                  ux(i,ii,iii)=uxsh(ni)
                  uy(i,ii,iii)=uysh(ni)
                  uz(i,ii,iii)=uzsh(ni)
                  ni=ni+1
          END IF
        END DO
      END DO
    END DO
    
    DEALLOCATE(uxsh,uysh,uzsh,kxsh,kysh,kzsh,handed,modk)
  END DO

  CALL symmetrize(ux,k2)
  CALL symmetrize(uy,k2)
  CALL symmetrize(uz,k2)

  CALL Skewness (ux,kx,S)
  WRITE(*,*) 'Skewness after matching:', S

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  
  atmp = (ux*CONJG(wx)+uy*CONJG(wy)+uz*CONJG(wz)) &
       +(CONJG(ux)*wx+CONJG(uy)*wy+CONJG(uz)*wz)
  atmp(1,:,:)=0.5_SP*atmp(1,:,:)
  DO i=1,lx
    htmp(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).lt.0.5_SP))
  END DO
  atmp = ux*CONJG(ux) + uy*CONJG(uy) + uz*CONJG(uz)
  atmp(1,:,:)=0.5_SP*atmp(1,:,:)
  DO i=1,lx
    etmp(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).lt.0.5_SP))
  END DO

  OPEN(16,FILE='scaledspec.dat')
  DO i=1,lx
    WRITE(16,*) i, etmp(i), htmp(i)
  END DO
  CLOSE(16)

  OPEN(16,FILE='scaledvel.dat',FORM='unformatted')
    WRITE(16) ux
    WRITE(16) uy
    WRITE(16) uz
  CLOSE(16)
  
  DEALLOCATE(ux,uy,uz,wx,wy,wz,kx,ky,kz,k2,hkinit,ekinit,htmp,etmp,atmp)
  CALL destroyplan3d

END PROGRAM InitVel      

