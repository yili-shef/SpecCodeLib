PROGRAM dataproc
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, PARAMETER :: k_c=8
  INTEGER  :: nx,ny,nz,itout,ieout,model
  INTEGER  :: nxb,nyb,nzb,lx,ly,lz,lxb,lyb,lzb,lxb1,lx1
  REAL(SP) :: dt,rnu
  REAL(SP) :: delta,eps,eta,delta_c

  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,G
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: TEk,Ek,DEk,PiEk,SGSDEk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: THk,Hk,DHk,PiHk,SGSDHk
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: E11,E22,E33,V11,V22,V33
  REAL(SP),    ALLOCATABLE, DIMENSION (:)     :: SF
  
  INTEGER :: ifile, startfno, numfile, i, id1, id2, id3, id4, ii
  REAL(SP) :: tt, S, ignore_me
  CHARACTER(80) :: path, sufix

  path='./'
  sufix='_dns.data'
  startfno=1
  numfile=101
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli_rot.d',status='unknown')
    READ(90,*) 
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
    READ(90,*) 
    READ(90,*) 
    READ(90,*) itout
    READ(90,*) ieout
    READ(90,*) dt
    READ(90,*) rnu
  CLOSE(90)
  CALL fftwplan3d(nx,ny,nz)
  
  nxb = nx*3/2 ;  nyb  = ny*3/2 ;  nzb = nz*3/2
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lxb = lx*3/2 ;  lyb  = ly*3/2 ;  lzb = lz*3/2
  lx1 = lx+1   ;  lxb1 = lxb+1

  delta=pi/REAL(lx,SP)
  delta_c=pi/REAL(k_c,SP)
  
  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz),G(lx1,ly,lz))
  ALLOCATE(tau11(lx1,ly,lz),tau12(lx1,ly,lz),tau13(lx1,ly,lz))
  ALLOCATE(tau22(lx1,ly,lz),tau23(lx1,ly,lz),tau33(lx1,ly,lz))
  ALLOCATE(TEk(lx),Ek(lx),DEk(lx),PiEk(lx),SGSDEk(lx))
  ALLOCATE(THk(lx),Hk(lx),DHk(lx),PiHk(lx),SGSDHk(lx))
  ALLOCATE(E11(lx),E22(lx),E33(lx),V11(lx),V22(lx),V33(lx))
  ALLOCATE(SF(lx))
  
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  G(2:lx1,:,:)=2._SP
  G(1,:,:)=1._SP
  WHERE(k2 .GE. lx*lx) G=0._SP
  DO i=1,lx
    SF(i)=SUM(G,mask=(ABS(SQRT(k2)-i).LT..5_SP))
    tt=(4._SP/3._SP)*Pi*((i+.5_SP)**3-(i-.5_SP)**3)
    SF(i)=tt/(SF(i)+smallest)
  END DO

  WHERE (k2 .GT. k_c*k_c) 
    G=0._SP
  ELSEWHERE
    G=1._SP
  END WHERE
  

  OPEN(19,FILE='Global'//sufix(1:LEN_TRIM(sufix)))
  WRITE(19,*) 'variables = "t", "k", "`e","T_E","h","`h","T_H", "S"' 
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
  OPEN(30,FILE='FGlobal'//sufix(1:LEN_TRIM(sufix)))
  WRITE(30,*) 'variables = "t", "k_f","`e_f","D_E_s_g_s","h_f","`h_f","D_H_s_g_s"'
  OPEN(31,FILE='Ekaa'//sufix(1:LEN_TRIM(sufix)))
  WRITE(31,*) 'Variables = "`k", "E11(`k)", "E22(`k)", "E33(`k)"'
  OPEN(32,FILE='Vkaa'//sufix(1:LEN_TRIM(sufix)))
  WRITE(32,*) 'Variables = "`k", "V11(`k)", "V22(`k)", "V33(`k)"'

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

    CALL energy_budget(ux,uy,uz,kx,ky,kz,k2,rnu,TEk,PiEk,Ek,DEk,lx1,ly,lz,nx,ny,nz)
    CALL helicity_budget(ux,uy,uz,kx,ky,kz,k2,rnu,THk,PiHk,Hk,DHk,lx1,ly,lz,nx,ny,nz)
    CALL enercompspec(ux,uy,uz,kx,ky,kz,k2,E11,E22,E33,V11,V22,V33,lx1,ly,lz,nx,ny,nz)
    CALL skewness(ux,kx,S,lx1,ly,lz,nx,ny,nz)
    CALL sgsstress(ux,uy,uz,tau11,tau22,tau33,tau12,tau13,tau23,G,lx1,ly,lz,lx,nx,ny,nz)
    

    CALL sgsenerdiss_spec(ux,uy,uz,tau11,tau22,tau33,tau12,tau13,tau23,kx,ky,kz,k2, &
                          SGSDEk,G,lx1,ly,lz,lx)
    CALL sgshelidiss_spec(ux,uy,uz,tau11,tau22,tau33,tau12,tau13,tau23,kx,ky,kz,k2, &
                          SGSDHk,G,lx1,ly,lz,lx)
    WRITE(19,'(27E14.6)') tt, SUM(Ek),SUM(DEk),SUM(TEk),SUM(Hk),SUM(DHk),SUM(THk),S
    WRITE(30,'(27E14.6)') tt, SUM(Ek(1:k_c)),SUM(DEK(1:k_c)),SUM(SGSDEk),SUM(Hk(1:k_c)),SUM(DHk(1:k_c)),SUM(SGSDHk)
    
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
    WRITE(31,*) ' Zone T="',tt,'",I=',lx,',F=point'
    WRITE(32,*) ' Zone T="',tt,'",I=',lx,',F=point'
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
      WRITE(31,*) i,E11(i)*SF(i),E22(i)*SF(i),E33(i)*SF(i)
      WRITE(32,*) i,V11(i)*SF(i),V22(i)*SF(i),V33(i)*SF(i)
    END DO 
  END DO
  CLOSE(32)
  CLOSE(31)
  CLOSE(30)
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

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2,TEk,Ek,DEk,THk,Hk,DHk)
  DEALLOCATE(tau11,tau12,tau13,tau22,tau23,tau33)
  DEALLOCATE(PiEk,PiHk,SGSDEk,SGSDHk)
  DEALLOCATE(E11,E22,E33,V11,V22,V33)
  DEALLOCATE(G,SF)
  CALL destroyplan3d

END PROGRAM dataproc
