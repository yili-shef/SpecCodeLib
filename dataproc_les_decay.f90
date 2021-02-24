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
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: TEk,Ek,DEk,PiEk,SGSDEk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: THk,Hk,DHk,PiHk,SGSDHk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: SF
  
  INTEGER :: ifile, startfno, numfile, i, id1, id2, id3, id4
  REAL(SP) :: tt, S, k_test,delta_test,MM,PM,LM,PP,LP,dyn_c1,dyn_c2,cs2,NM,NN,LN
  CHARACTER(80) :: path, sufix

  path='./'
  sufix='_dynsmag.data'
  startfno=1
  numfile=51
  model=6
  ! 1: Smag
  ! 2: local helical 
  ! 3: nonlocal helical
  ! 4: dyn_a 
  ! 5: dyn_b
  ! 6: dynsmag
  ! 7: nonlinear
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli_dynsmag.d',status='unknown')
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
  ALLOCATE(tmp(lx1,ly,lz),G_test(lx1,ly,lz),SF(lx))
  ALLOCATE(TEk(lx),Ek(lx),DEk(lx),PiEk(lx),SGSDEk(lx))
  ALLOCATE(THk(lx),Hk(lx),DHk(lx),PiHk(lx),SGSDHk(lx))
  
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  tmp(2:lx1,:,:)=2._SP
  tmp(1,:,:)=1._SP
  WHERE(k2 .GE. lx*lx) tmp=0._SP
  DO i=1,lx
    SF(i)=SUM(tmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT..5_SP))
    tt=(4._SP/3._SP)*Pi*((i+1._SP)**3-i**3)
    SF(i)=tt/(SF(i)+smallest)
  END DO

  WHERE(k2 .LE. k_test*k_test)
          G_test=1._SP
  ELSEWHERE
          G_test=0._SP
  ENDWHERE

  OPEN(18,FILE='GlobalE'//sufix(1:LEN_TRIM(sufix)))
  WRITE(18,*) 'Variables = "t", "k","`e","T_E","D_E_s_g_s","S"'
  WRITE(18,*) 'zone I=',numfile,',F=Point'
  OPEN(19,FILE='GlobalH'//sufix(1:LEN_TRIM(sufix)))
  WRITE(19,*) 'variables = "t", "h","`h","T_H","D_H_s_g_s"' 
  WRITE(19,*) 'zone I=',numfile,',F=Point'
  OPEN(20,FILE=    'Ek'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "E(`k)"'
  OPEN(21,FILE=    'DEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(21,*) 'Variables = "`k", "D_E(`k)"'
  OPEN(22,FILE= 'TEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(22,*) 'Variables = "`k", "T_E(`k)"'
  OPEN(23,FILE=    'Hk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(23,*) 'Variables = "`k", "H(`k)"'
  OPEN(24,FILE=   'DHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(24,*) 'Variables = "`k", "D_H(`k)"'
  OPEN(25,FILE='THk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(25,*) 'Variables = "`k", "T_H(`k)"'
  OPEN(26,FILE='PiEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(26,*) 'Variables = "`k", "`P_E(`k)"'
  OPEN(27,FILE='PiHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(27,*) 'Variables = "`k", "`P_H(`k)"'
  OPEN(28,FILE='SGSDEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(28,*) 'Variables = "`k", "D_E_s_g_s(`k)"'
  OPEN(29,FILE='SGSDHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(29,*) 'Variables = "`k", "D_H_s_g_s(`k)"'

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


    tt = (ifile-1)*itout*dt 

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
    CALL skewness(ux,kx,S,lx1,ly,lz,nx,ny,nz)
    WRITE(18,'(20E14.6)') tt,SUM(Ek),SUM(DEk),SUM(TEk),SUM(SGSDEk), S
    WRITE(19,'(20E14.6)') tt,SUM(Hk),SUM(DHk),SUM(THk),SUM(SGSDHk)
    
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(21,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(22,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(23,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(24,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(25,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(26,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(27,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(28,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(29,*) ' Zone T="',tt,'",I=',lx,',F=point'
    DO i=1,lx
      WRITE(20,*) i,Ek(i)*SF(i)
      WRITE(21,*) i,DEk(i)*SF(i)
      WRITE(22,*) i,TEk(i)*SF(i)
      WRITE(23,*) i,Hk(i)*SF(i)
      WRITE(24,*) i,DHk(i)*SF(i)
      WRITE(25,*) i,THk(i)*SF(i)
      WRITE(26,*) i,PiEk(i)*SF(i)
      WRITE(27,*) i,PiHk(i)*SF(i)
      WRITE(28,*) i,SGSDEk(i)*SF(i)
      WRITE(29,*) i,SGSDHk(i)*SF(i)
    END DO 
  END DO
  CLOSE(29)
  CLOSE(28)
  CLOSE(27)
  CLOSE(26)
  CLOSE(25)
  CLOSE(24)
  CLOSE(23)
  CLOSE(22)
  CLOSE(21)
  CLOSE(20)
  CLOSE(19)
  CLOSE(18)

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,TEk,Ek,DEk,THk,Hk,DHk)
  DEALLOCATE(wx,wy,wz,t11,t12,t13,t22,t23,t33)
  DEALLOCATE(PiEk,PiHk,SGSDEk,SGSDHk,tmp,SF,G_test)
  CALL destroyplan3d

END PROGRAM dataproc
