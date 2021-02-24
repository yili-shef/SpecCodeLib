PROGRAM alignment
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024,npt=400
  INTEGER,  PARAMETER :: if1=1,if2=20 !40 for filtered, 20 for DNS
  REAL(DP), PARAMETER :: dt=0.017
  REAL(DP), PARAMETER :: mnt=0.,mxt=1.!pi/2.
  
  REAL(DP), DIMENSION(3,3,nprtcl) :: EVecDNS, EVecRE
  REAL(DP), DIMENSION(3,3)        :: EVecDNSn, EVecREn
  INTEGER,  DIMENSION(nprtcl)     :: marker

  REAL(DP), DIMENSION(npt) :: xt
  REAL(DP), DIMENSION(npt) :: rpdfa,rpdfb,rpdfg
  INTEGER,  DIMENSION(npt) ::  pdfa, pdfb, pdfg

  REAL(DP) :: time,a1,a2,a3,normt,bwt,ran1,gasdev

  INTEGER :: i,ii,iii,ifile,ip,nn,ns,idum
  
  idum=-7
  
  bwt=(mxt-mnt)/npt
  xt=(/(mnt+(i-.5_DP)*bwt,i=1,npt)/)
  

  OPEN(20,FILE='marker.dat',FORM='formatted')
  READ(20,*) marker
  CLOSE(20)
  ns=SUM(marker)
  WRITE(*,*) 'ns=',ns

  OPEN(21,FILE='./xBdata/evtr_re.dat',FORM='unformatted')
  OPEN(20,FILE='./xBdata/evtr_dns.dat',FORM='unformatted')

  OPEN(14,FILE='align_pdfa.dat')
  OPEN(13,FILE='align_pdfb.dat')
  OPEN(12,FILE='align_pdfg.dat')
  WRITE(14,'(''VARIABLES = "xt","pdfa"'')')
  WRITE(13,'(''VARIABLES = "xt","pdfb"'')')
  WRITE(12,'(''VARIABLES = "xt","pdfg"'')')


  time=0.
  DO ifile=if1,if2

    READ(20) EVecDNS
    READ(21) EVecRE

    pdfa=0; pdfb=0; pdfg=0
    DO i=1,nprtcl
    
      IF (marker(i) .EQ. 0) CYCLE

      EVecDNSn=EVecDNS(:,:,i)
      EVecREn =EVecRE (:,:,i)
     
      DO ii=1,3
      EVecDNSn(:,ii)=EVecDNSn(:,ii)/SQRT(SUM(EVecDNSn(:,ii)**2))
      EVecREn(:,ii)=EVecREn(:,ii)/SQRT(SUM(EVecREn(:,ii)**2))
      END DO

!      a1=ACOS(ABS(DOT_PRODUCT(EVecDNSn(:,3),EVecREn(:,3))))
!      a2=ACOS(ABS(DOT_PRODUCT(EVecDNSn(:,2),EVecREn(:,2))))
!      a3=ACOS(ABS(DOT_PRODUCT(EVecDNSn(:,1),EVecREn(:,1))))
      a1=(ABS(DOT_PRODUCT(EVecDNSn(:,3),EVecREn(:,3))))
      a2=(ABS(DOT_PRODUCT(EVecDNSn(:,2),EVecREn(:,2))))
      a3=(ABS(DOT_PRODUCT(EVecDNSn(:,1),EVecREn(:,1))))
      normt=(a1-mnt)/bwt
      IF (normt .GE. 0._DP .AND. normt .LT. npt) THEN
              nn=FLOOR(normt)+1
              pdfa(nn)=pdfa(nn)+1
      END IF
      normt=(a2-mnt)/bwt
      IF (normt .GE. 0._DP .AND. normt .LT. npt) THEN
              nn=FLOOR(normt)+1
              pdfb(nn)=pdfb(nn)+1
      END IF
      normt=(a3-mnt)/bwt
      IF (normt .GE. 0._DP .AND. normt .LT. npt) THEN
              nn=FLOOR(normt)+1
              pdfg(nn)=pdfg(nn)+1
      END IF
    
    END DO
    rpdfa=DBLE(pdfa)/ns
    rpdfb=DBLE(pdfb)/ns
    rpdfg=DBLE(pdfg)/ns
    WRITE(*,*) 'check ca cb cg', SUM(rpdfa),SUM(rpdfb),SUM(rpdfg)
    rpdfa=rpdfa/bwt
    rpdfb=rpdfb/bwt
    rpdfg=rpdfg/bwt
    
    WRITE(14,*) 'ZONE T="t=',time,'"'
    DO i=1,npt
!      WRITE(14,'(20E12.3)') xt(i)/pi*180,rpdfa(i)/SIN(xt(i))
      WRITE(14,'(20E12.3)') xt(i),rpdfa(i)
    END DO
    WRITE(13,*) 'ZONE T="t=',time,'"'
    DO i=1,npt
!      WRITE(13,'(20E12.3)') xt(i)/pi*180,rpdfb(i)/SIN(xt(i))
      WRITE(13,'(20E12.3)') xt(i),rpdfb(i)
    END DO
    WRITE(12,*) 'ZONE T="t=',time,'"'
    DO i=1,npt
!      WRITE(12,'(20E12.3)') xt(i)/pi*180,rpdfg(i)/SIN(xt(i))
      WRITE(12,'(20E12.3)') xt(i),rpdfg(i)
    END DO

    time=time+dt
  END DO
  CLOSE(21)
  CLOSE(20)
  CLOSE(14)
  CLOSE(13)
  CLOSE(12)

  WRITE(*,*) 'finished'

END PROGRAM alignment
