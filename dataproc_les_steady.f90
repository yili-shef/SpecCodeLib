PROGRAM dataproc
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx,ny,nz,itout,ieout,model
  INTEGER  :: nxb,nyb,nzb,lx,ly,lz,lxb,lyb,lzb,lxb1,lx1
  REAL(SP) :: dt,rnu,time
  REAL(SP) :: delta,eps,eta,delta_ratio
  LOGICAL  :: update

  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz,wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: t11,t12,t13,t22,t23,t33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,tmp,G_test
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: TEk,Ek,DEk,PiEk,SGSDEk,SGSDEk_s
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: THk,Hk,DHk,PiHk,SGSDHk,SGSDHk_s
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: aTEk,aEk,aDEk,aPiEk,aSGSDEk,aSGSDEk_s
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: aTHk,aHk,aDHk,aPiHk,aSGSDHk,aSGSDHk_s
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: SF
  
  
  INTEGER  :: ifile, startfno, numfile, i, ii, id1, id2, id3, id4
  REAL(SP) :: tt, S, avrsgsed,avrsgshd,k_test,delta_test,MM,PM,PP,LP,LM,dyn_c1,dyn_c2,cs2,NM,NN,LN
  CHARACTER(80) :: path, sufix

  path='./'
  sufix='_smag.data'
  startfno=13
  numfile=27
  model=1
  ! 1: Smag
  ! 2: local helical 
  ! 3: nonlocal helical
  ! 4: dyn_a 
  ! 5: dyn_b
  ! 6: dynsmag
  ! 7: nonlinear
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli_smag.d',status='unknown')
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) itout
    READ(90,*) ieout
    READ(90,*) dt
    READ(90,*) rnu
    READ(90,*) time
    READ(90,*)
    READ(90,*)
    READ(90,*) eps
    READ(90,*) eta
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
!    READ(90,*) delta_ratio
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
  ALLOCATE(tmp(lx1,ly,lz),G_test(lx1,ly,lz),SF(lx))
  ALLOCATE(TEk(lx),Ek(lx),DEk(lx),PiEk(lx),SGSDEk(lx),SGSDEk_s(lx))
  ALLOCATE(THk(lx),Hk(lx),DHk(lx),PiHk(lx),SGSDHk(lx),SGSDHk_s(lx))
  ALLOCATE(aTEk(lx),aEk(lx),aDEk(lx),aPiEk(lx),aSGSDEk(lx),aSGSDEk_s(lx))
  ALLOCATE(aTHk(lx),aHk(lx),aDHk(lx),aPiHk(lx),aSGSDHk(lx),aSGSDHk_s(lx))
  
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  tmp(2:lx1,:,:)=2._SP
  tmp(1,:,:)=1._SP
  WHERE(k2 .GE. lx*lx) tmp=0._SP
  DO i=1,lx
    SF(i)=SUM(tmp,mask=(ABS(SQRT(k2)-i).LT..5_SP))
    tt=(4._SP/3._SP)*Pi*((i+.5)**3-(i-0.5)**3)
    SF(i)=tt/(SF(i)+smallest)

    WRITE(*,*) SF(i)
  END DO
  STOP
  
  WHERE(k2 .LE. k_test*k_test)
          G_test=1._SP
  ELSEWHERE
          G_test=0._SP
  ENDWHERE

  aTEk=0._SP
  aEk=0._SP
  aDEk=0._SP
  aPiEk=0._SP
  aSGSDEk=0._SP
  aTHk=0._SP
  aHk=0._SP
  aDHk=0._SP
  aPiHk=0._SP
  aSGSDHk=0._SP
  aSGSDEk_s=0._SP
  aSGSDHk_s=0._SP
  DO ifile=startfno,startfno+numfile-1
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
      update=.true.
      CALL dynheli_a (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,PP,LM,LP,   &
                      dyn_c1,dyn_c2,update)
    CASE (5)
      update=.true.
      CALL dynheli_b (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,PP,LM,LP,   &
                      dyn_c1,dyn_c2,update)
    CASE (6)
      update=.true.
      CALL dynsmag (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,  &
                    lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,LM,cs2,  &
                    update)
    CASE (7)
      update=.true.
      CALL nonlinear (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,    &
                      lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,NM,NN,LM,LN,dyn_c1, &
                      dyn_c2,update)
    CASE DEFAULT
      WRITE(*,*) "model # wrong, model = ", model
      STOP
    END SELECT

    CALL energy_budget_les(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,k2,rnu,TEk,PiEk,Ek,DEk,SGSDEk, &
                           lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
    CALL helicity_budget_les(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,k2,rnu,THk,PiHk,Hk,DHk,SGSDHk, &
                             lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)

    aTEk=aTEk+TEk
    aPiEk=aPiEk+PiEk
    aEk=aEk+Ek
    aDEk=aDEk+DEk
    aSGSDEk=aSGSDEk+SGSDEk

    aTHk=aTHk+THk
    aPiHk=aPiHk+PiHk
    aHk=aHk+Hk
    aDHk=aDHk+DHk
    aSGSDHk=aSGSDHk+SGSDHk

    IF (model .NE. 1 .AND. model .NE. 6) THEN
                             
      CALL smag(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz, &
                lxb1,lyb,lzb,nxb,nyb,nzb,delta)
      IF (model .EQ. 4 .OR. model .EQ. 5) THEN
              tt=dyn_c1/c1
      ELSE
              tt=1._SP
      END IF
      t11=t11*tt
      t12=t12*tt
      t13=t13*tt
      t22=t22*tt
      t23=t23*tt
      t33=t33*tt
 
      CALL energy_budget_les(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,k2,rnu,TEk,PiEk,Ek,DEk, &
                             SGSDEk_s,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
      CALL helicity_budget_les(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,k2,rnu,THk,PiHk,Hk,DHk, &
                               SGSDHk_s,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
       
      aSGSDEk_s=aSGSDEk_s+SGSDEk_s
      aSGSDHk_s=aSGSDHk_s+SGSDHk_s
    END IF

  END DO
  aTEk=aTEk/REAL(numfile,SP)*SF
  aPiEk=aPiEk/REAL(numfile,SP)*SF
  aEk=aEk/REAL(numfile,SP)*SF
  aDEk=aDEk/REAL(numfile,SP)*SF
  aSGSDEk=aSGSDEk/REAL(numfile,SP)*SF
  aTHk=aTHk/REAL(numfile,SP)*SF
  aPiHk=aPiHk/REAL(numfile,SP)*SF
  aHk=aHk/REAL(numfile,SP)*SF
  aSGSDHk=aSGSDHk/REAL(numfile,SP)*SF

  aSGSDEk_s=aSGSDEk_s/REAL(numfile,SP)*SF
  aSGSDHk_s=aSGSDHk_s/REAL(numfile,SP)*SF

  OPEN(20,FILE='EnerBudget'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,'('' Variables = "`k","E(`k)","D_E(`k)","T_E(`k)","`P_E(`k)","D_E_s_g_s(`k)","D^S_E_s_g_s(`k)","D^H_E_s_g_s(`k)"'')')
  WRITE(20,*) 'zone I=',lx,',F=Point'
  DO ii=1,lx
    WRITE(20,'(I4,20E14.6)') ii,aEk(ii),aDEk(ii),aTEk(ii),aPiEk(ii),aSGSDEk(ii),aSGSDEk_s(ii),aSGSDEk(ii)-aSGSDEk_s(ii)
  END DO
  CLOSE(20)

  OPEN(20,FILE='HeliBudget'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,'('' Variables = "`k","H(`k)","D_H(`k)","T_H(`k)","`P_H(`k)","D_H_s_g_s(`k)","D^S_H_s_g_s(`k)","D^H_H_s_g_s(`k)"'')')
  WRITE(20,*) 'zone I=',lx,',F=Point'
  DO ii=1,lx
    WRITE(20,'(I4,20E14.6)') ii,aHk(ii),aDHk(ii),aTHk(ii),aPiHk(ii),aSGSDHk(ii),aSGSDHk_s(ii),aSGSDHk(ii)-aSGSDHk_s(ii)
  END DO
  CLOSE(20)

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,TEk,Ek,DEk,THk,Hk,DHk)
  DEALLOCATE(wx,wy,wz,t11,t12,t13,t22,t23,t33)
  DEALLOCATE(PiEk,PiHk,SGSDEk,SGSDHk,tmp,SF,G_test)
  DEALLOCATE(aTEk,aEk,aDEk,aTHk,aHk,aDHk,aPiEk,aPiHk,aSGSDEk,aSGSDHk)
  DEALLOCATE(SGSDEk_s,SGSDHk_s,aSGSDEk_s,aSGSDHk_s)
  CALL destroyplan3d

END PROGRAM dataproc
