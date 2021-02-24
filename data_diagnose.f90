PROGRAM diagnose
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx,ny,nz,model,iupdate,itout
  INTEGER  :: nxb,nyb,nzb,lx,ly,lz,lxb,lyb,lzb,lxb1,lx1
  REAL(SP) :: rnu,eps,eta,delta_ratio
  LOGICAL  :: update

  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz,wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: t11,t12,t13,t22,t23,t33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,tmp,G_test
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: sgsed,sgshd
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: visced,vischd
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: SF

  INTEGER  :: ifile, startfno, numfile, i, ii, id1, id2, id3, id4, jstep
  REAL(SP) :: asgsed,asgshd,avisced,avischd
  REAL(SP) :: tt,S,k_test,delta,delta_test,MM,PM,PP,LP,LM,dyn_c1,dyn_c2,cs2
  CHARACTER(80) :: path, sufix

  path='./'
  sufix='_nl.data'
  startfno=1
  numfile=51
  model=7
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli_nl.d',status='unknown')
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) itout
    READ(90,*) 
    READ(90,*)
    READ(90,*) rnu
    READ(90,*) 
    READ(90,*)
    READ(90,*)
    READ(90,*) eps
    READ(90,*) eta
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) delta_ratio
  CLOSE(90)

  CALL fftwplan3d(nx,ny,nz)

  nxb = nx*3/2 ;  nyb  = ny*3/2 ;  nzb = nz*3/2
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lxb = lx*3/2 ;  lyb  = ly*3/2 ;  lzb = lz*3/2
  lx1 = lx+1   ;  lxb1 = lxb+1

  delta=pi/REAL(lx,SP)
  delta_test=delta_ratio*delta
  k_test=lx/delta_ratio
  

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  ALLOCATE(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(sgsed(nx,ny,nz),sgshd(nx,ny,nz))
  ALLOCATE(visced(nx,ny,nz),vischd(nx,ny,nz))
  ALLOCATE(tmp(lx1,ly,lz),G_test(lx1,ly,lz),SF(lx))

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  tmp(2:lx1,:,:)=2._SP
  tmp(1,:,:)=1._SP
  WHERE(k2 .GE. lx*lx) tmp=0._SP
  DO i=1,lx
    SF(i)=SUM(tmp,mask=(ABS(SQRT(k2)-i).LT..5_SP))
    tt=(4._SP/3._SP)*Pi*((i+.5_SP)**3-(i-.5_SP)**3)
    SF(i)=tt/(SF(i)+smallest)
  END DO
  
  WHERE(k2 .LE. k_test*k_test)
          G_test=1._SP
  ELSEWHERE
          G_test=0._SP
  ENDWHERE

  OPEN(15,FILE=path(1:LEN_TRIM(path))//'check1.dat',STATUS='unknown')
  DO ifile=startfno,startfno+numfile-1
    jstep=(ifile-1)*itout
    WRITE(*,*) ifile
    id1=MOD(ifile,10)+48
    id2=MOD(INT(ifile/10),10)+48
    id3=MOD(INT(ifile/100),10)+48
    id4=MOD(INT(ifile/1000),10)+48
    OPEN(16,FILE=path(1:LEN_TRIM(path))//'out/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat', &
         FORM='unformatted')
      READ(16) ux
      READ(16) uy
      READ(16) uz
    CLOSE(16)

    wx=eye*(ky*uz-kz*uy)
    wy=eye*(kz*ux-kx*uz)
    wz=eye*(kx*uy-ky*ux)

    update=.true.
    SELECT CASE (model)
    CASE (1) ! Smagorinsky model
      CALL smag(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz, &
                lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (2) ! Local helicity model
      CALL localheli(ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                     lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (3) ! Nonlocal helicity model
      CALL nonlocalheli(ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                        lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)
    CASE (4)
      CALL dynheli_a (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,PP,LM,LP,   &
                      dyn_c1,dyn_c2,update)
    CASE (5)
      CALL dynheli_b (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,PP,LM,LP,   &
                      dyn_c1,dyn_c2,update)
    CASE (6)
      CALL dynsmag (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,  &
                    lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,LM,cs2,  &
                    update)
    CASE (7)
      CALL nonlinear (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,    &
                      lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,PP,LM,LP,dyn_c1, &
                      dyn_c2,update)
    CASE DEFAULT
      WRITE(*,*) "model # wrong, model = ", model
      STOP
    END SELECT

    CALL sgsenerdiss (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,sgsed,lx1,lx,ly,lz,nx,ny,nz)
    CALL sgshelidiss (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,sgshd,lx1,lx,ly,lz,nx,ny,nz)
    asgsed=SUM(sgsed)/REAL(nx*ny*nz,SP)
    asgshd=SUM(sgshd)/REAL(nx*ny*nz,SP)
    CALL viscenerdiss (ux,uy,uz,kx,ky,kz,rnu,visced,lx1,lx,ly,lz,nx,ny,nz)
    CALL vischelidiss (ux,uy,uz,wx,wy,wz,kx,ky,kz,rnu,vischd,lx1,lx,ly,lz,nx,ny,nz)
    avisced=SUM(visced)/REAL(nx*ny*nz,SP)
    avischd=SUM(vischd)/REAL(nx*ny*nz,SP)

    WRITE(15,'(I4,20E15.4)') ifile,eps,asgsed,avisced,asgsed+avisced,eta,asgshd,avischd,asgshd+avischd
    
  END DO
  CLOSE(15)

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2)
  DEALLOCATE(wx,wy,wz,t11,t12,t13,t22,t23,t33)
  DEALLOCATE(tmp,G_test,SF)
  CALL destroyplan3d



END PROGRAM diagnose
