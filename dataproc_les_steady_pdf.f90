PROGRAM dataproc
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, PARAMETER :: npnt=800
  REAL(SP), DIMENSION(npnt) :: xpdf,kdypdf
  REAL(SP), DIMENSION(npnt) :: ypdf,hdypdf

  
  INTEGER  :: nx,ny,nz,itout,ieout,model
  INTEGER  :: nxb,nyb,nzb,lx,ly,lz,lxb,lyb,lzb,lxb1,lx1
  REAL(SP) :: dt,rnu,time
  REAL(SP) :: delta,eps,eta
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz,wx,wy,wz
  COMPLEX(SP), ALLOCATABLE, DIMENSION (:,:,:) :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: sgsdiss
  REAL(SP),    ALLOCATABLE, DIMENSION (:,:,:) :: kx,ky,kz,k2,tmp
  
  INTEGER  :: ifile, startfno, numfile, i, ii, id1, id2, id3, id4
  REAL(SP) :: tt, S, avrsgsed,avrsgshd,sdvsgsed,sdvsgshd
  CHARACTER(80) :: path, sufix

  path='./'
  sufix='_smag.data'
  startfno=20
  numfile=32
  OPEN(90,file=path(1:LEN_TRIM(path))//'parameter_heli.d',status='unknown')
    READ(90,*) model
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
    READ(90,*) time
    READ(90,*)
    READ(90,*)
    READ(90,*) eps
    READ(90,*) eta
  CLOSE(90)
  CALL fftwplan3d(nx,ny,nz)
  
  nxb = nx*3/2 ;  nyb  = ny*3/2 ;  nzb = nz*3/2
  lx  = nx/2   ;  ly   = ny     ;  lz  = nz
  lxb = lx*3/2 ;  lyb  = ly*3/2 ;  lzb = lz*3/2
  lx1 = lx+1   ;  lxb1 = lxb+1

  delta=pi/REAL(lx,SP)
  
  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  ALLOCATE(tau11(lx1,ly,lz),tau12(lx1,ly,lz),tau13(lx1,ly,lz))
  ALLOCATE(tau22(lx1,ly,lz),tau23(lx1,ly,lz),tau33(lx1,ly,lz))
  ALLOCATE(kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))
  ALLOCATE(sgsdiss(nx,ny,nz))
  
  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  avrsgsed=0._SP
  avrsgshd=0._SP
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
      CALL smag(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,lx1,ly,lz, &
                lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (2) ! Local helicity model
      CALL localheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                     lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (3) ! Nonlocal helicity model
      CALL nonlocalheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                        lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)
    CASE DEFAULT
      WRITE(*,*) "model # wrong, model = ", model
      STOP
    END SELECT

    CALL pdf_sgsenerdiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    avrsgsed=avrsgsed+SUM(sgsdiss)/REAL(nx*ny*nz,SP)
    CALL pdf_sgshelidiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    avrsgshd=avrsgshd+SUM(sgsdiss)/REAL(nx*ny*nz,SP)

  END DO
  avrsgsed=avrsgsed/REAL(numfile,SP)
  avrsgshd=avrsgshd/REAL(numfile,SP)

  sdvsgsed=0._SP
  sdvsgshd=0._SP
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
      CALL smag(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,lx1,ly,lz, &
                lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (2) ! Local helicity model
      CALL localheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                     lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (3) ! Nonlocal helicity model
      CALL nonlocalheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                        lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)
    CASE DEFAULT
      WRITE(*,*) "model # wrong, model = ", model
      STOP
    END SELECT

    CALL pdf_sgsenerdiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    sdvsgsed=sdvsgsed+SQRT(SUM((sgsdiss-avrsgsed)**2)/REAL(nx*ny*nz,SP))
    CALL pdf_sgshelidiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    sdvsgshd=sdvsgshd+SQRT(SUM((sgsdiss-avrsgshd)**2)/REAL(nx*ny*nz,SP))

  END DO
  sdvsgsed=sdvsgsed/REAL(numfile,SP)
  sdvsgshd=sdvsgshd/REAL(numfile,SP)

  kdypdf=0._SP
  hdypdf=0._SP
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
      CALL smag(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,lx1,ly,lz, &
                lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (2) ! Local helicity model
      CALL localheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                     lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
    CASE (3) ! Nonlocal helicity model
      CALL nonlocalheli(ux,uy,uz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                        lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)
    CASE DEFAULT
      WRITE(*,*) "model # wrong, model = ", model
      STOP
    END SELECT

    CALL pdf_sgsenerdiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    CALL pdfsampling(sgsdiss,avrsgsed,sdvsgsed,40._SP,npnt,ypdf,xpdf,nx,ny,nz)
    kdypdf=kdypdf+ypdf
    CALL pdf_sgshelidiss(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,sgsdiss,lx1,lx,ly,lz,nx,ny,nz)
    CALL pdfsampling(sgsdiss,avrsgshd,sdvsgshd,40._SP,npnt,ypdf,xpdf,nx,ny,nz)
    hdypdf=hdypdf+ypdf

  END DO
  kdypdf=kdypdf/REAL(numfile,SP)
  hdypdf=hdypdf/REAL(numfile,SP)

  
  OPEN(20,FILE='AvrSDV'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'ener avr and sdv:', avrsgsed, sdvsgsed
  WRITE(20,*) 'heli avr and sdv:', avrsgshd, sdvsgshd
  CLOSE(20)

  OPEN(20,FILE='EnerDissPdf'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables="(`P_E-`m_E)/`s_E", "PDF"'
  WRITE(20,*) 'zone I=',npnt,',F=Point'
  DO ii=1,npnt
  WRITE(20,*) xpdf(ii),kdypdf(ii)
  END DO
  CLOSE(20)

  OPEN(20,FILE='HeliDissPdf'//sufix(1:LEN_TRIM(sufix)))
  WRITE(20,*) 'Variables="(`P_H-`m_H/`s_H", "PDF"'
  WRITE(20,*) 'zone I=',npnt,',F=Point'
  DO ii=1,npnt
  WRITE(20,*) xpdf(ii),hdypdf(ii)
  END DO
  CLOSE(20)

  DEALLOCATE(ux,uy,uz,kx,ky,kz,k2)
  DEALLOCATE(wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33)
  CALL destroyplan3d

END PROGRAM dataproc
