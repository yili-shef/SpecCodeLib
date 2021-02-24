PROGRAM ensemble_avr
  USE mconstant
  IMPLICIT NONE

  INTEGER, PARAMETER :: nensmble=10
  CHARACTER(80) :: path,sufix

  INTEGER :: startfno,numfile,irun,ifile,id1,id2,model,nx,ny,nz,lx,ly,lz,lx1
  INTEGER :: i,ieout,itout,itt
  REAL(SP) :: time,dt,rnu,tt,tk,tde,tte,tsgsde,th,tdh,tth,tsgsdh,ts
  REAL(SP), ALLOCATABLE, DIMENSION(:) :: sk,sde,ste,ssgsde,sh,sdh,sth,ssgsdh,ss
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: Ek,DEk,TEk,Hk,DHk,THk,PiEk,PiHk,SGSDEk,SGSDHk

  path='/home/yili/Workspace/RawData/LES_DynSmag_Unst_Heli_64/'
  sufix='_dynsmag.data'
  startfno=1
  numfile=51
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
  CLOSE(90)
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lx1 = lx + 1



  OPEN(19,FILE='../GlobalE'//sufix(1:LEN_TRIM(sufix)))
  WRITE(19,*) 'Variables = "t", "k", "`e","T_E","D_E_s_g_s","S"'
  WRITE(19,*) 'zone I=',numfile,',F=Point'
  OPEN(20,FILE='../GlobalH'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "t", "h","`h","T_H","D_H_s_g_s"' 
  WRITE(20,*) 'zone I=',numfile,',F=Point'
  
  ALLOCATE(sk(numfile),sde(numfile),ste(numfile),ssgsde(numfile))
  ALLOCATE(sh(numfile),sdh(numfile),sth(numfile),ssgsdh(numfile))
  ALLOCATE(ss(numfile))
  sk=0._SP; sde=0._SP; ste=0._SP; ssgsde=0._SP; sh=0._SP; sdh=0._SP; sth=0._SP; ssgsdh=0._SP; ss=0._SP

  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='GlobalE'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    READ(77,*)
    OPEN(57,FILE='GlobalH'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(57,*)
    READ(57,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,'(20E14.6)') tt,tk,tde,tte,tsgsde,ts
      READ(57,'(20E14.6)') tt,th,tdh,tth,tsgsdh
      sk(ifile)=sk(ifile)+tk
      sde(ifile)=sde(ifile)+tde
      ste(ifile)=ste(ifile)+tte
      ssgsde(ifile)=ssgsde(ifile)+tsgsde
      sh(ifile)=sh(ifile)+th
      sdh(ifile)=sdh(ifile)+tdh
      sth(ifile)=sth(ifile)+tth
      ssgsdh(ifile)=ssgsdh(ifile)+tsgsdh
      ss(ifile)=ss(ifile)+ts
    END DO
    CLOSE(77)
    CLOSE(57)
  END DO
  sk=sk/REAL(nensmble)
  sde=sde/REAL(nensmble)
  ste=ste/REAL(nensmble)
  ssgsde=ssgsde/REAL(nensmble)
  sh=sh/REAL(nensmble)
  sdh=sdh/REAL(nensmble)
  sth=sth/REAL(nensmble)
  ssgsdh=ssgsdh/REAL(nensmble)
  ss=ss/REAL(nensmble)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(19,'(20E14.6)') tt,sk(ifile),sde(ifile),ste(ifile),ssgsde(ifile),ss(ifile)
    WRITE(20,'(20E14.6)') tt,sh(ifile),sdh(ifile),sth(ifile),ssgsdh(ifile)
  END DO
  DEALLOCATE(sk,sde,ste,ssgsde,sh,sdh,sth,ssgsdh,ss)
  CLOSE(19)
  CLOSE(20)


  OPEN(20,FILE=    '../Ek'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "E(`k)"'
  ALLOCATE(Ek(lx,numfile))
  Ek=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='Ek'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        Ek(i,ifile)=Ek(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  Ek=Ek/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,Ek(i,ifile)
    END DO
  END DO
  DEALLOCATE(Ek)
  CLOSE(20)


  OPEN(20,FILE='../DEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "D_E(`k)"'
  ALLOCATE(DEk(lx,numfile))
  DEk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='DEk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        DEk(i,ifile)=DEk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  DEk=DEk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,DEk(i,ifile)
    END DO
  END DO
  DEALLOCATE(DEk)
  CLOSE(20)

  OPEN(20,FILE='../TEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "T_E(`k)"'
  ALLOCATE(TEk(lx,numfile))
  TEk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='TEk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        TEk(i,ifile)=TEk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  TEk=TEk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,TEk(i,ifile)
    END DO
  END DO
  DEALLOCATE(TEk)
  CLOSE(20)

  OPEN(20,FILE='../Hk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "H(`k)"'
  ALLOCATE(Hk(lx,numfile))
  Hk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='Hk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        Hk(i,ifile)=Hk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  Hk=Hk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,Hk(i,ifile)
    END DO
  END DO
  DEALLOCATE(Hk)
  CLOSE(20)

  OPEN(20,FILE='../DHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "D_H(`k)"'
  ALLOCATE(DHk(lx,numfile))
  DHk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='DHk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        DHk(i,ifile)=DHk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  DHk=DHk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,DHk(i,ifile)
    END DO
  END DO
  DEALLOCATE(DHk)
  CLOSE(20)


  OPEN(20,FILE='../THk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "T_H(`k)"'
  ALLOCATE(THk(lx,numfile))
  THk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='THk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        THk(i,ifile)=THk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  THk=THk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,THk(i,ifile)
    END DO
  END DO
  DEALLOCATE(THk)
  CLOSE(20)


  OPEN(20,FILE='../PiEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "`P_E(`k)"'
  ALLOCATE(PiEk(lx,numfile))
  PiEk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='PiEk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        PiEk(i,ifile)=PiEk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  PiEk=PiEk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,PiEk(i,ifile)
    END DO
  END DO
  DEALLOCATE(PiEk)
  CLOSE(20)


  OPEN(20,FILE='../PiHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "`P_H(`k)"'
  ALLOCATE(PiHk(lx,numfile))
  PiHk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='PiHk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        PiHk(i,ifile)=PiHk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  PiHk=PiHk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,PiHk(i,ifile)
    END DO
  END DO
  DEALLOCATE(PiHk)
  CLOSE(20)


  OPEN(20,FILE='../SGSDEk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "D_E_s_g_s(`k)"'
  ALLOCATE(SGSDEk(lx,numfile))
  SGSDEk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='SGSDEk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        SGSDEk(i,ifile)=SGSDEk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  SGSDEk=SGSDEk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,SGSDEk(i,ifile)
    END DO
  END DO
  DEALLOCATE(SGSDEk)
  CLOSE(20)


  OPEN(20,FILE='../SGSDHk'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables = "`k", "D_H_s_g_s(`k)"'
  ALLOCATE(SGSDHk(lx,numfile))
  SGSDHk=0._SP
  DO irun=1, nensmble
    id1=MOD(irun,10)+48
    id2=INT(irun/10)+48
    OPEN(77,FILE='SGSDHk'//CHAR(id2)//CHAR(id1)//sufix(1:LEN_TRIM(sufix)))
    READ(77,*)
    DO ifile=startfno,startfno+numfile-1
      READ(77,*)
      DO i=1,lx
        READ(77,*) itt, tt
        SGSDHk(i,ifile)=SGSDHk(i,ifile)+tt
      END DO
    END DO
    CLOSE(77)
  END DO    
  SGSDHk=SGSDHk/REAL(nensmble,SP)
  DO ifile=startfno,startfno+numfile-1
    tt = (ifile-1)*itout*dt 
    WRITE(20,*) ' Zone T="',tt,'",I=',lx,',F=point' 
    DO i=1,lx
      WRITE(20,*) i,SGSDHk(i,ifile)
    END DO
  END DO
  DEALLOCATE(SGSDHk)
  CLOSE(20)

  WRITE(*,*) 'ensmbleavr done'
END PROGRAM ensemble_avr      
