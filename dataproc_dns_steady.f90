PROGRAM dataproc
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER  :: nx,ny,nz
  INTEGER  :: lx,ly,lz,lx1
  REAL(SP) :: rnu,delta_c

  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz,wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: t11,t12,t13,t22,t23,t33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,tmp,G
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: TEk,Ek,DEk,PiEk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: FTEk,FEk,FDEk,FPiEk,SGSDEk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: aTEk,aEk,aDEk,aPiEk,aSGSDEk,aFTEk,aFPiEk,aFEk,aFDEk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: SF
  
  
  INTEGER  :: ifile, startfno, numfile, i, ii, id1, id2, id3, id4
  REAL(SP) :: tt, S
  CHARACTER(80) :: path, sufix


  
  INTEGER, PARAMETER :: k_c=32
  delta_c=pi/REAL(k_c,SP)

  path='./'
  sufix='_dns.data'
  startfno=10
  numfile=32
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_dns.d',status='unknown')
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) 
    READ(90,*) rnu
  CLOSE(90)
  CALL fftwplan3d(nx,ny,nz)
  
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lx1 = lx+1 

  
  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  ALLOCATE(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(tmp(lx1,ly,lz),G(lx1,ly,lz),SF(lx))
  ALLOCATE(TEk(lx),Ek(lx),DEk(lx),PiEk(lx))
  ALLOCATE(FTEk(lx),FEk(lx),FDEk(lx),FPiEk(lx),SGSDEk(lx))
  ALLOCATE(aTEk(lx),aEk(lx),aDEk(lx),aPiEk(lx),aSGSDEk(lx),aFTEk(lx),aFPiEk(lx))
  ALLOCATE(aFEk(lx),aFDEk(lx))
  
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  tmp(2:lx1,:,:)=2._SP
  tmp(1,:,:)=1._SP
  WHERE(k2 .GE. lx*lx) tmp=0._SP
  DO i=1,lx
    SF(i)=SUM(tmp,mask=(ABS(SQRT(k2)-i).LT..5_SP))
    tt=(4._SP/3._SP)*Pi*((i+.5_SP)**3-(i-.5_SP)**3)
    SF(i)=tt/(SF(i)+smallest)
  END DO
  
  WHERE(k2 .LE. k_c*k_c)
          G=1._SP
  ELSEWHERE
          G=0._SP
  ENDWHERE

  aTEk=0._SP;  aEk=0._SP;  aDEk=0._SP;  aPiEk=0._SP;  aSGSDEk=0._SP
  aFTEk=0._SP; aFPiEk=0._SP; aFEk=0._SP; aFDEk=0._SP
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

    CALL energy_budget(ux,uy,uz,kx,ky,kz,k2,rnu,TEk,PiEk,Ek,DEk,lx1,ly,lz,nx,ny,nz)

    aTEk=aTEk+TEk;       aPiEk=aPiEk+PiEk;     aEk=aEk+Ek;       aDEk=aDEk+DEk

    CALL sgsstress(ux,uy,uz,t11,t22,t33,t12,t13,t23,G,lx1,ly,lz,lx,nx,ny,nz)
!    ux=G*ux
!    uy=G*uy
!    uz=G*uz
    WHERE(k2 .GT. k_c*k_c)
            ux=0._SP
            uy=0._SP
            uz=0._SP
    END WHERE
    CALL energy_budget_fdns(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,k2,rnu,FTEk,FPiEk,FEk,FDEk,SGSDEk, &
                          lx1,ly,lz,nx,ny,nz)

    aFEk=aFEk+FEk; aFDEk=aFDEk+FDEk; aFTEk=aFTEk+FTEk; aFPiEk=aFPiEk+FPiEk; aSGSDEk=aSGSDEk+SGSDEk

  END DO
  aTEk=aTEk/REAL(numfile,SP)*SF;    aPiEk=aPiEk/REAL(numfile,SP)*SF
  aEk=aEk/REAL(numfile,SP)*SF;      aDEk=aDEk/REAL(numfile,SP)*SF
  aFTEk=aFTEk/REAL(numfile,SP)*SF;  aFPiEk=aFPiEk/REAL(numfile,SP)*SF
  aFEk=aFEk/REAL(numfile,SP)*SF;    aFDEk=aFDEk/REAL(numfile,SP)*SF
  aSGSDEk=aSGSDEk/REAL(numfile,SP)*SF

  OPEN(20,FILE='EnerBudget'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,'('' Variables = "k","E(k)","D_E(k)","T_E(k)","P_E(k)"'')')
  WRITE(20,*) 'zone I=',lx,',F=Point'
  DO ii=1,lx
    WRITE(20,'(I4,20E14.6)')ii,aEk(ii),aDEk(ii),aTEk(ii),aPiEk(ii)
  END DO
  CLOSE(20)

  OPEN(20,FILE='EnerBudget_f'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,'('' Variables = "k","FE(k)","FD_E(k)","FT_E(k)","FP_E(k)","D_E_s_g_s(k)"'')')
  WRITE(20,*) 'zone I=',lx,',F=Point'
  DO ii=1,lx
    WRITE(20,'(I4,20E14.6)')ii,aFEk(ii),aFDEk(ii),aFTEk(ii),aFPiEk(ii),aSGSDEk(ii)
  END DO
  CLOSE(20)


  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,TEk,Ek,DEk,PiEk)
  DEALLOCATE(wx,wy,wz,t11,t12,t13,t22,t23,t33,FTEk,FEk,FDEk,FPiEk)
  DEALLOCATE(SGSDEk,tmp,SF,G,aFEk,aFDEk,aFTEk,aFPiEk,)
  DEALLOCATE(aTEk,aEk,aDEk,aPiEk,aSGSDEk)
  CALL destroyplan3d

END PROGRAM dataproc
