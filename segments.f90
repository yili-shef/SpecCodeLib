program seg
      if(ifile==150) then
              open(80,file='helidissA_2.dat',status='unknown')
                  write(80,*) 'VARIABLES="X","Y","helidissA"'
                  write(80,*) 'ZONE I=',2*lx,',J=',ly,',F=POINT'
                  do ii=1,ny
                  do i=1,nx
                  write(80,'(2I6,E18.5)') i, ii, helidissA(i,ii,50)
                  END do
                  END do
              close(80)
      END if

  call rfftwnd_f77_one_complex_to_real(C2R3D,vx,ignore_me)

  open(30,file='vx.dat',status='unknown')
  write(30,*) 'VARIABLES="X","Y","vx"'
  write(30,*) 'ZONE I=',2*lx,',J=',ly,',F=POINT'
  do ii=1,ly
  do i=1,lx
  write(30,'(2I6,E18.5)') 2*i-1, ii, REAL(DP)(vx(i,ii,50))
  write(30,'(2I6,E18.5)') 2*i,   ii, aimag(vx(i,ii,50))
  END do
  END do
  close(30)
  vx=vx/REAL(DP)(nx*ny*nz)
  call rfftwnd_f77_one_real_to_complex(R2C3D,vx,ignore_me)
  
  call symmpart(vx,vy,vz,tau11,tau12,tau13,tau22,tau23,tau33,lx1,ly,lz,kx,ky,kz)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,tau33,ignore_me)
  
  open(30,file='S11.dat',status='unknown')
  write(30,*) 'VARIABLES="X","Y","S11"'
  write(30,*) 'ZONE I=',2*lx,',J=',ly,',F=POINT'
  do ii=1,ly
  do i=1,lx
  write(30,'(2I6,E18.5)') 2*i-1, ii, REAL(DP)(tau11(i,ii,50))
  write(30,'(2I6,E18.5)') 2*i,   ii, aimag(tau11(i,ii,50))
  END do
  END do
  close(30)
  
  S = cmplx(                                                             &
             sqrt(                                                       &
                  2*(REAL(DP)(tau11)**2+REAL(DP)(tau22)**2+REAL(DP)(tau33)**2        &
                  +2.*(REAL(DP)(tau12)**2+REAL(DP)(tau13)**2+REAL(DP)(tau23)**2))    &
                 )                                                       &
           ,                                                             &
             sqrt(                                                       &
                  2*(aimag(tau11)**2+aimag(tau22)**2+aimag(tau33)**2     &
                  +2.*(aimag(tau12)**2+aimag(tau13)**2+aimag(tau23)**2)) &
                 )                                                       &
           )

  tau11=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau11),aimag(S)*aimag(tau11))
  tau12=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau12),aimag(S)*aimag(tau12))
  tau13=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau13),aimag(S)*aimag(tau13))
  tau22=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau22),aimag(S)*aimag(tau22))
  tau23=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau23),aimag(S)*aimag(tau23))
  tau33=-2.*c1*delta**2*cmplx(REAL(DP)(S)*REAL(DP)(tau33),aimag(S)*aimag(tau33))
  
  open(30,file='tau1.dat',status='unknown')
  write(30,*) 'VARIABLES="X","Y","tau"'
  write(30,*) 'ZONE I=',2*lx,',J=',ly,',F=POINT'
  do ii=1,ly
  do i=1,lx
  write(30,'(2I6,E18.5)') 2*i-1, ii, REAL(DP)(tau11(i,ii,50))
  write(30,'(2I6,E18.5)') 2*i,   ii, aimag(tau11(i,ii,50))
  END do
  END do
  close(30)
  open(57,file='heli_lambda6_10.dat',status='unknown')
  write(57,*) "# time, lambda6, lambda7, lambda8, lambda9, lambda10"
  open(58,file='heli_enerdiss_model3.dat',status='unknown')
  write(58,*) "# time, enerdiss3A, enerdiss3B, sumAB, eps, SGSeps"
  open(59,file='heli_enerdiss_model4.dat',status='unknown')
  write(59,*) "# time, enerdiss4A, enerdiss4B, sumAB, eps, SGSeps"
  open(60,file='heli_enerdiss_model5.dat',status='unknown')
  write(60,*) "# time, enerdiss5A, enerdiss5B, sumAB, eps, SGSeps"
  open(61,file='heli_helidiss_model3.dat',status='unknown')
  write(61,*) "# time, helidiss3A, helidiss3B, sumAB, eta, SGSeta" 
  open(62,file='heli_helidiss_model4.dat',status='unknown')
  write(62,*) "# time, helidiss4A, helidiss4B, sumAB, eta, SGSeta" 
  open(63,file='heli_helidiss_model5.dat',status='unknown')
  write(63,*) "# time, helidiss5A, helidiss5B, sumAB, eta, SGSeta" 
  open(64,file='heli_lambda11.dat',status='unknown')
  write(64,*) "# time, lambda11"
  open(65,file='heli_SROP.dat', status='unknown')
  write(65,*) "# time, avrS2, avrR2, avrSR, avrOmega2, avrP2, avrOmegaP"
  open(66,file='heli_dnsc.dat', status='unknown')
  write(66,*) "# time, dnsc1_model1,dnsc2_model1, dnsc1_model2, dnsbeta_model2"
  open(67,file='heli_enerdiss_dns.dat', status='unknown')
  write(67,*) "# time, dnsenerdiss1A, dnsenerdiss1B, dnsenerdiss2A, dnsenerdiss2B"
  open(68,file='heli_helidiss_dns.dat', status='unknown')
  write(68,*) "# time, dnshelidiss1A, dnshelidiss1B, dnshelidiss2A, dnshelidiss2B"

    write(58,'(6E15.5)') ifile*200*dt, enerdiss3A,enerdiss3B,enerdiss3A+enerdiss3B,eps,avrenerflux
    write(59,'(6E15.5)') ifile*200*dt, enerdiss4A,enerdiss4B,enerdiss4A+enerdiss4B,eps,avrenerflux
    write(60,'(6E15.5)') ifile*200*dt, enerdiss5A,enerdiss5B,enerdiss5A+enerdiss5B,eps,avrenerflux
    write(61,'(6E15.5)') ifile*200*dt, helidiss3A,helidiss3B,helidiss3A+helidiss3B,eta,avrheliflux
    write(62,'(6E15.5)') ifile*200*dt, helidiss4A,helidiss4B,helidiss4A+helidiss4B,eta,avrheliflux
    write(63,'(6E15.5)') ifile*200*dt, helidiss5A,helidiss5B,helidiss5A+helidiss5B,eta,avrheliflux
    write(65,'(7E15.5)') ifile*200*dt, avrS2, avrR2, avrSR, avrOmega2, avrP2, avrOmegaP

  close(57)
  close(58)
  close(59)
  close(60)
  close(61)
  close(62)
  close(63)
  close(64)
  close(65)
  close(66)
  close(67)
  close(68)

  if (id==0 .and. ifile==52) then
  open(81, file='enerflux.dat', status='unknown')
  open(82, file='heliflux.dat', status='unknown')
  write(81,*) 'VARIABLES="X","Y","enerflux"'
  write(81,*) 'ZONE I=', 2*lx, ',J=', ly, ',F=POINT'
  write(82,*) 'VARIABLES="X","Y","heliflux"'
  write(82,*) 'ZONE I=', 2*lx, ',J=', ly, ',F=POINT'
  do ii=1,ly
  do i=1,2*lx
  write(81,'(2I6,E15.5)') i, ii,enerflux(i,ii,8)
  write(82,'(2I6,E15.5)') i, ii,heliflux(i,ii,8)
  END do
  END do
  close(81)
  close(82)
  END if

--Lambda6
  avrsth=sum((SR+OmegaP)*S*SR/(sqrt(.5*(S**2+Omega**2))*sqrt(.5*(R**2+P**2))))
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda6=(tmp/REAL(DP)(lx2*ly*nz))*sqrt(avrS2+avrOmega2)*sqrt(avrR2+avrP2) &
          /((avrSR+avrOmegaP)*sqrt(avrS2)*avrSR)


--Lambda7
  avrsth=sum((SR+OmegaP)*S*.5*R**2/(sqrt(.5*(S**2+Omega**2))*sqrt(.5*(R**2+P**2))))
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda7=(tmp/REAL(DP)(lx2*ly*nz))*sqrt(avrS2+avrOmega2)*sqrt(avrR2+avrP2) &
          /((avrSR+avrOmegaP)*sqrt(avrS2)*avrR2)

  if (id==0 .and. ifile==55) then
  open(80, file='SR_RR.dat', status='unknown')
  write(80,*) 'VARIABLES="X","Y","SR+OmegaP","RR"'
  write(80,*) 'ZONE I=', 2*lx, ',J=', ly, ',F=POINT'
  do ii=1,ly
  do i=1,lx*2
  write(80,'(2I6,2E15.5)') i, ii, SR(i,ii,8)+OmegaP(i,ii,8), 0.5*R(i,ii,8)**2
  END do
  END do
  close(80)
  END if

--Lambda8
  avrsth=sum((SR+OmegaP)*S*SR/(.5*(S**2+Omega**2)))
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda8=(tmp/REAL(DP)(lx2*ly*nz))*(avrS2+avrOmega2) &
          /((avrSR+avrOmegaP)*sqrt(avrS2)*avrSR)

--Lambda9
  avrsth=sum((SR+OmegaP)*S*.5*R**2/(.5*(S**2+Omega**2)))
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda9=(tmp/REAL(DP)(lx2*ly*nz))*(avrS2+avrOmega2) &
          /((avrSR+avrOmegaP)*sqrt(avrS2)*avrR2)

--Lambda10
  avrsth=sum(2.*SR**2/S)
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda10=(tmp/REAL(DP)(lx2*ly*nz))*sqrt(avrS2)/avrSR**2

--Lambda11
  avrsth=sum(2.*SR*.5*R**2/S)
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
  lambda11=(tmp/REAL(DP)(lx2*ly*nz))*sqrt(avrS2)/(avrSR*avrR2)

  if (id==0) then
    write(51,'(6E15.5)') ifile*200*dt, lambda1,lambda2,lambda3,lambda4,lambda5
    write(57,'(6E15.5)') ifile*200*dt, lambda6,lambda7,lambda8,lambda9,lambda10
    write(64,'(6E15.5)') ifile*200*dt, lambda11
  END if 

  lambda1=1.
  lambda2=1.
  lambda3=1.
  lambda4=1.
  lambda5=1.
  lambda6=1.
  lambda7=1.
  lambda8=1.
  lambda9=1.
  lambda10=1.
  lambda11=1.


  A11=2.*sqrt(2.)*delta**2*lambda1*avrS2**1.5
  A12=-delta**3*lambda3*sqrt(avrS2)*avrSR
  A21=4.*delta**2*lambda3*sqrt(avrS2)*avrSR
  A22=-2.*delta**3*lambda5*sqrt(avrS2)*avrR2

  dnsc1_model1=(eps*A22-eta*A12)/(A22*A11-A12*A21)
  dnsc2_model1=(A11*eta-A21*eps)/(A22*A11-A12*A21)

  A12=-delta**3*lambda2*avrSR**2/sqrt(avrR2)
  A22=-2.*delta**3*lambda4*sqrt(avrR2)*avrSR

  dnsc1_model2=(eps*A22-eta*A12)/(A22*A11-A12*A21)
  dnsbeta_model2=(A11*eta-A21*eps)/(A22*A11-A12*A21)

  if(id==0) then
  write(66,'(6E15.5)') ifile*200*dt, dnsc1_model1, dnsc2_model1, dnsc1_model2, dnsbeta_model2
  END if
 Linearized version
  c1=1./(lambda1*2.*sqrt(2.)*pi**2*ck**1.5*I13**1.5)
  c2=(eta*delta/eps)*(lambda3*ch-sqrt(2.)*lambda1*ck)/  &
     (2.*sqrt(2.)*lambda1*lambda5*ck**2.5*pi**4*I13**.5*I73)
  beta = -(sqrt(2.)*lambda1*ck-ch*lambda3)  &
         /(sqrt(2.)*pi**3*lambda1*lambda4*ck**1.5*ch*I13*I73**.5)

 Full version
  det = (4.*sqrt(2.)*lambda1*lambda5*ck**3*pi**6*I13**2*I73 &
     -lambda3**2*ck*I13**3*pi**4*(eta*delta/eps)**2)
  c1=(2.*lambda5*ck**1.5*pi**4*I13**1.5*I73+ &
     .5*lambda3*ck**.5*ch*I13**1.5*pi**2*(eta*delta/eps))/det 
  c2=(eta*delta/eps)*(2.*lambda3*ck**0.5*ch*I13**1.5*pi**2- &
     2.*sqrt(2.)*lambda1*ck**1.5*pi**2*I13**1.5)/det 

  enerdiss1A=2.*sqrt(2.)*c1*delta**2*lambda1*avrS2**1.5
  enerdiss2A=enerdiss1A
  enerdiss3A=enerdiss1A
  enerdiss4A=enerdiss1A
  enerdiss5A=enerdiss1A

  helidiss1A=4.*c1*delta**2*lambda3*sqrt(avrS2)*avrSR
  helidiss2A=helidiss1A
  helidiss3A=helidiss1A
  helidiss4A=helidiss1A
  helidiss5A=helidiss1A

  enerdiss1B=-c2*delta**3*lambda3*sqrt(avrS2)*avrSR
  enerdiss2B=-beta*delta**3*lambda2*(avrSR**2/sqrt(avrR2))
  enerdiss3B=-beta*delta**3*lambda6*((avrSR+avrOmegaP)*sqrt(avrS2)*avrSR) &
             /(sqrt(avrS2+avrOmega2)*sqrt(avrR2+avrP2))
  enerdiss4B=-beta*delta**4*lambda8*((avrSR+avrOmegaP)*sqrt(avrS2)*avrSR) &
             /(avrS2+avrOmega2)
  enerdiss5B=-beta*delta**4*lambda10**avrSR**2/sqrt(avrS2)

  helidiss1B=-2.*c2*delta**3*lambda5*sqrt(avrS2)*avrR2
  helidiss2B=-2.*beta*delta**3*lambda4*sqrt(avrR2)*avrSR
  helidiss3B=-2.*beta*delta**3*lambda7*((avrSR+avrOmegaP)*sqrt(avrS2)*avrR2) &
             /(sqrt(avrS2+avrOmega2)*sqrt(avrR2+avrP2))
  helidiss4B=-2.*beta*delta**4*lambda9*((avrSR+avrOmegaP)*sqrt(avrS2)*avrR2) &
             /(avrS2+avrOmega2)
  helidiss5B=-2.*beta*delta**4*lambda11*avrSR*avrR2/sqrt(avrS2)

  tmp=4.*dnsc1_model2*delta**2*lambda3*sqrt(avrS2)*avrSR-2.*dnsbeta_model2*delta**3*lambda4*sqrt(avrR2)*avrSR

  dnsenerdiss1A=2.*sqrt(2.)*dnsc1_model1*delta**2*lambda1*avrS2**1.5
  dnsenerdiss2A=2.*sqrt(2.)*dnsc1_model2*delta**2*lambda1*avrS2**1.5
  dnsenerdiss1B=-dnsc2_model1*delta**3*lambda3*sqrt(avrS2)*avrSR
  dnsenerdiss2B=-dnsbeta_model2*delta**3*lambda2*(avrSR**2/sqrt(avrR2))
  dnshelidiss1A=4.*dnsc1_model1*delta**2*lambda3*sqrt(avrS2)*avrSR
  dnshelidiss2A=4.*dnsc1_model2*delta**2*lambda3*sqrt(avrS2)*avrSR
  dnshelidiss1B=-2.*dnsc2_model1*delta**3*lambda5*sqrt(avrS2)*avrR2
  dnshelidiss2B=-2.*dnsbeta_model2*delta**3*lambda4*sqrt(avrR2)*avrSR

  if(id==0) then
    write(67,'(6E15.5)') ifile*200*dt, dnsenerdiss1A, dnsenerdiss1B, dnsenerdiss2A, dnsenerdiss2B
    write(68,'(6E15.5)') ifile*200*dt, dnshelidiss1A, dnshelidiss1B, dnshelidiss2A, dnshelidiss2B
  END if

  pi = 2.*asin(1.)
  delta=pi/REAL(DP)(lx)

  open(52,file='heli_enerdiss_model2.dat',status='unknown')
  write(52,*) "# time, enerdiss2A, enerdiss2B, sumAB, eps, SGSeps"
  open(53,file='heli_helidiss_model2.dat',status='unknown')
  write(53,*) "# time, helidiss2A, helidiss2B, sumAB, eta, SGSeta"
  open(54,file='compenenerspec.dat', status='unknown')
  open(55,file='compenhelispec.dat', status='unknown')

  istartno = 50
  numfile=70

!---Loop reading data files
  do ifile=istartno, istartno+numfile

  if(id==0)  write(*,*) "file number:", ifile
  i4d=int(ifile/1000)
  i5d=int((ifile-i4d*1000)/100)
  i6d=int((ifile-i4d*1000-i5d*100)/10)
  i7d=mod(ifile,10)

  fnm1='./out/vel'//char(i4d+48)           &
       //char(i5d+48)//char(i6d+48)      &
       //char(i7d+48)//'.dat'

  open(30, file=fnm1, status='unknown',form='unformatted')
  read(30) ux
  read(30) uy
  read(30) uz
  close(30)


  call symmpart(ux,uy,uz,S11,S12,S13,S22,S23,S33,lx1,ly,lz,kx,ky,kz)
  call curl(ux,uy,uz,curlx,curly,curlz,lx1,ly,lz,kx,ky,kz)
  call symmpart(curlx,curly,curlz,R11,R12,R13,R22,R23,R33,lx1,ly,lz,kx,ky,kz)

  call mpifft3DCR(S11,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(S22,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(S33,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(S12,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(S13,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(S23,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R11,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R22,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R33,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R12,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R13,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(R23,lx,ly,lz,nproc,id,icouples)

!--Calculate S, R, and SR fileds.
  call CalSR

!--Calculate the lambda's, the modeled dissipation terms.
  call ModelDiss

  eps = enerdiss2A
  eta = helidiss2A

!--Calculate compensated energy and helicity spectra
  call spec

  if(id==0) then
    write(52,'(6E15.5)') ifile*itout*dt, enerdiss2A,enerdiss2B,enerdiss2A+enerdiss2B,eps
    write(53,'(6E15.5)') ifile*itout*dt, helidiss2A,helidiss2B,helidiss2A+helidiss2B,eta
  END if

  END do

  if(id==0) then
  close(52)
  close(53)
  close(54)
  close(55)
  END if

  if(nproc.gt.1)then
    deallocate( icouples )
  endif
  deallocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  deallocate(curlx(lx1,ly,lz),curly(lx1,ly,lz),curlz(lx1,ly,lz))
  deallocate(S11(lx1,ly,lz),S22(lx1,ly,lz),S33(lx1,ly,lz))
  deallocate(S12(lx1,ly,lz),S13(lx1,ly,lz),S23(lx1,ly,lz))
  deallocate(R11(lx1,ly,lz),R22(lx1,ly,lz),R33(lx1,ly,lz))
  deallocate(R12(lx1,ly,lz),R13(lx1,ly,lz),R23(lx1,ly,lz))
  deallocate(G(lx1,ly,lz),kx(lx1,ly,lz),ky(lx1,ly,lz),kz(lx1,ly,lz),k2(lx1,ly,lz))

  deallocate(S(lx2,ly,lz),R(lx2,ly,lz),SR(lx2,ly,lz),Omega(lx2,ly,lz),P(lx2,ly,lz),OmegaP(lx2,ly,lz))
  deallocate(tau11(lx2,ly,lz),tau12(lx2,ly,lz),tau13(lx2,ly,lz))
  deallocate(tau22(lx2,ly,lz),tau23(lx2,ly,lz),tau33(lx2,ly,lz))
  deallocate(enerflux(lx2,ly,lz),heliflux(lx2,ly,lz)) 
  deallocate(tmparr(lx1,ly,lz) )

  call destroy_plan1D
  call destroy_plan2D

  call mpi_finalize(ierr)

CONTAINS

  SUBROUTINE CalSR
!--Calculate instantaneous S, R and SR fields.
    do i=1,lx
      S(2*i-1,:,:) = sqrt(2*(REAL(DP)(S11(i,:,:))**2+REAL(DP)(S22(i,:,:))**2 &
                     +REAL(DP)(S33(i,:,:))**2+2.*(REAL(DP)(S12(i,:,:))**2   &
                     +REAL(DP)(S13(i,:,:))**2+REAL(DP)(S23(i,:,:))**2)))
      S(2*i,:,:) = sqrt(2*(aimag(S11(i,:,:))**2+aimag(S22(i,:,:))**2 &
                     +aimag(S33(i,:,:))**2+2.*(aimag(S12(i,:,:))**2 &
                     +aimag(S13(i,:,:))**2+aimag(S23(i,:,:))**2)))
      R(2*i-1,:,:) = sqrt(2*(REAL(DP)(R11(i,:,:))**2+REAL(DP)(R22(i,:,:))**2 &
                     +REAL(DP)(R33(i,:,:))**2+2.*(REAL(DP)(R12(i,:,:))**2     &
                     +REAL(DP)(R13(i,:,:))**2+REAL(DP)(R23(i,:,:))**2)))
      R(2*i,:,:) = sqrt(2*(aimag(R11(i,:,:))**2+aimag(R22(i,:,:))**2 &
                     +aimag(R33(i,:,:))**2+2.*(aimag(R12(i,:,:))**2     &
                     +aimag(R13(i,:,:))**2+aimag(R23(i,:,:))**2)))
      SR(2*i-1,:,:) = REAL(DP)(S11(i,:,:))*REAL(DP)(R11(i,:,:)) &
                     +REAL(DP)(S22(i,:,:))*REAL(DP)(R22(i,:,:))  &
                     +REAL(DP)(S33(i,:,:))*REAL(DP)(R33(i,:,:))  &
                     +2*(REAL(DP)(S12(i,:,:))*REAL(DP)(R12(i,:,:)) &
                     +REAL(DP)(S13(i,:,:))*REAL(DP)(R13(i,:,:))    &
                     +REAL(DP)(S23(i,:,:))*REAL(DP)(R23(i,:,:)))
      SR(2*i,:,:) = aimag(S11(i,:,:))*aimag(R11(i,:,:)) &
                     +aimag(S22(i,:,:))*aimag(R22(i,:,:))  &
                     +aimag(S33(i,:,:))*aimag(R33(i,:,:))  &
                     +2*(aimag(S12(i,:,:))*aimag(R12(i,:,:)) &
                     +aimag(S13(i,:,:))*aimag(R13(i,:,:))    &
                     +aimag(S23(i,:,:))*aimag(R23(i,:,:)))
    END do

  END SUBROUTINE CalSR

  SUBROUTINE ModelDiss
!--Calculate modeled dissipations and lambda's.
    avrsth=sum((0.5*S**2)**1.5)
    tmp=0.
    call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
    tmp = tmp/REAL(DP)(lx2*ly*nz)
 
    enerdiss2A=2.*c1*sqrt(2.)*delta**2*tmp
 
    avrsth=sum(S*SR/sqrt(2.))
    tmp=0.
    call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
    tmp = tmp/REAL(DP)(lx2*ly*nz)
 
    helidiss2A=4.*sqrt(2.)*c1*delta**2*tmp
 
    avrsth=sum(SR**2/(R/sqrt(2.)))
    tmp=0.
    call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
    tmp = tmp/REAL(DP)(lx2*ly*nz)
 
    enerdiss2B=-(sqrt(2.)/2.)*beta*delta**3*tmp
 
    avrsth=sum(SR*R/sqrt(2.))
    tmp=0.
    call mpi_allreduce(avrsth,tmp,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
    tmp=tmp/REAL(DP)(lx2*ly*nz)
 
    helidiss2B=-sqrt(2.)*beta*delta**3*tmp

  END SUBROUTINE ModelDiss

  SUBROUTINE spec
    tmparr = ux*conjg(ux) + uy*conjg(uy) + uz*conjg(uz)
    tmparr(1,:,:)=0.5*tmparr(1,:,:)
    if(id==0) write(54,*) dt*ifile*itout
    do i=1,nek
      e_t = 0.0
      ek=sum(tmparr(lx1,ly,lz),mask=(abs(sqrt(k2)-i-0.499999).lt.0.5))
      call mpi_allreduce(ek,e_t,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
      if (id==0) then
        write(54,*) i, e_t
      END if
    END do

    tmparr = (ux*conjg(curlx)+uy*conjg(curly)+uz*conjg(curlz)) &
         +(conjg(ux)*curlx+conjg(uy)*curly+conjg(uz)*curlz)
    tmparr(1,:,:) = 0.5*tmparr(1,:,:)
    if(id==0) write(55,*) dt*ifile*itout
    do i=1,nek
      he_t = 0.0
      hek=sum(tmparr(lx1,ly,lz),mask=(abs(sqrt(k2)-i-0.499999).lt.0.5))
!    MPI_allreduce broadcasts the results to all nodes, while MPI_reduce
!    only returns the results to the root.
      call mpi_allreduce(hek,he_t,1,MPI_REAL,MPI_SUM,nallgrp,ierr)
      if (id==0) then
        write(55,*) i, he_t
      END if
    END do
  END SUBROUTINE spec
SUBROUTINE sgsstress (vx,vy,vz,tau11,tau22,tau33,tau12,tau13,tau23, G,  &
                      lx,ly,lz,nproc,id,icouples)
 
  INTEGER, DIMENSION(nproc-1) :: icouples
  INTEGER :: lx,ly,lz,nproc,id, i
  COMPLEX(DP), DIMENSION(lx+1,ly,lz) :: vx, vy, vz, vxt, vyt, vzt, arrt
  REAL(DP), DIMENSION(lx+1,ly,lz) :: G
  REAL(DP), DIMENSION(2*lx,ly,lz) :: tau11,tau12,tau13,tau22,tau23,tau33

  vxt=vx 
  vyt=vy
  vzt=vz

  call mpifft3DCR(vx,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(vy,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(vz,lx,ly,lz,nproc,id,icouples)

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vx(i,:,:))*REAL(DP)(vx(i,:,:)),aimag(vx(i,:,:))*aimag(vx(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau11(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau11(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vy(i,:,:))*REAL(DP)(vy(i,:,:)),aimag(vy(i,:,:))*aimag(vy(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau22(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau22(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vz(i,:,:))*REAL(DP)(vz(i,:,:)),aimag(vz(i,:,:))*aimag(vz(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau33(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau33(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vx(i,:,:))*REAL(DP)(vy(i,:,:)),aimag(vx(i,:,:))*aimag(vy(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau12(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau12(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vx(i,:,:))*REAL(DP)(vz(i,:,:)),aimag(vx(i,:,:))*aimag(vz(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau13(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau13(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  do i = 1, lx
  arrt(i,:,:) = cmplx(REAL(DP)(vy(i,:,:))*REAL(DP)(vz(i,:,:)),aimag(vy(i,:,:))*aimag(vz(i,:,:)))
  END do
  call mpifft3DRC(arrt,lx,ly,lz,nproc,id,icouples)
  arrt=G*arrt
  call mpifft3DCR(arrt,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  tau23(2*i-1,:,:)=REAL(DP)(arrt(i,:,:))
  tau23(2*i,:,:) = aimag(arrt(i,:,:))
  END do

  vx=G*vxt
  vy=G*vyt
  vz=G*vzt

  call mpifft3DCR(vx,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(vy,lx,ly,lz,nproc,id,icouples)
  call mpifft3DCR(vz,lx,ly,lz,nproc,id,icouples)

  do i=1,lx
  tau11(2*i-1,:,:)=tau11(2*i-1,:,:)-(REAL(DP)(vx(i,:,:)))**2
  tau11(2*i,:,:)=tau11(2*i,:,:)-(aimag(vx(i,:,:)))**2
  tau22(2*i-1,:,:)=tau22(2*i-1,:,:)-(REAL(DP)(vy(i,:,:)))**2
  tau22(2*i,:,:)=tau22(2*i,:,:)-(aimag(vy(i,:,:)))**2
  tau33(2*i-1,:,:)=tau33(2*i-1,:,:)-(REAL(DP)(vz(i,:,:)))**2
  tau33(2*i,:,:)=tau33(2*i,:,:)-(aimag(vz(i,:,:)))**2
  tau12(2*i-1,:,:)=tau12(2*i-1,:,:)-REAL(DP)(vx(i,:,:))*REAL(DP)(vy(i,:,:))
  tau12(2*i,:,:)=tau12(2*i,:,:)-aimag(vx(i,:,:))*aimag(vy(i,:,:))
  tau13(2*i-1,:,:)=tau13(2*i-1,:,:)-REAL(DP)(vx(i,:,:))*REAL(DP)(vz(i,:,:))
  tau13(2*i,:,:)=tau13(2*i,:,:)-aimag(vx(i,:,:))*aimag(vz(i,:,:))
  tau23(2*i-1,:,:)=tau23(2*i-1,:,:)-REAL(DP)(vy(i,:,:))*REAL(DP)(vz(i,:,:))
  tau23(2*i,:,:)=tau23(2*i,:,:)-aimag(vy(i,:,:))*aimag(vz(i,:,:))
  END do

  vx=vxt
  vy=vyt
  vz=vzt
END SUBROUTINE sgsstress

SUBROUTINE gradskew (ux,uy,uz,lx,ly,lz,kx,ky,kz,nproc,id,icouples,gradskewx)
  IMPLICIT NONE
  include 'mpif.h'

  INTEGER :: i, lx,ly,lz, ierr, nproc, id
  INTEGER, DIMENSION(nproc-1) :: icouples
  COMPLEX(DP), DIMENSION(lx+1,ly,lz) :: ux, uy, uz, dudx 
  REAL(DP), DIMENSION(lx+1,ly,lz) :: kx,ky,kz
  REAL(DP), DIMENSION(lx*2,ly,lz) :: arrt
  REAL(DP) :: gradskewx, avrsth, tmp
  COMPLEX(DP), PARAMETER :: eye=(0.,1.)

  dudx=eye*ux
  call mpifft3DCR(dudx,lx,ly,lz,nproc,id,icouples)
  do i = 1, lx
  arrt(2*i-1,:,:)=REAL(DP)(dudx(i,:,:))
  arrt(2*i,:,:)=aimag(dudx(i,:,:))
  END do
  avrsth = sum(arrt(1:lx*2,:,:)**3)
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1, mpi_real, mpi_sum,mpi_comm_world,ierr)
  gradskewx=tmp/REAL(DP)(2*lx*ly*lz*nproc)
  avrsth=sum(arrt(1:lx,:,:)**2)
  tmp=0.
  call mpi_allreduce(avrsth,tmp,1, mpi_real, mpi_sum,mpi_comm_world,ierr)
  gradskewx=gradskewx/(tmp/REAL(DP)(2*lx*ly*lz**nproc))**1.5

END SUBROUTINE gradskew
!  if (istep==1) then 
!    open(77,file='test1.dat')
!    do i=1,ny
!    write(77,'(50E12.5)') (REAL(DP)(vx(2,i,iii)),aimag(vx(2,i,iii)), iii=1,nz)
!    END do
!   close(77)
!   open(78,file='test2.dat')
!   do i=1,nyb
!   write(78,'(50E12.5)') (REAL(DP)(vxb(2,i,iii)),aimag(vxb(2,i,iii)),iii=1,nzb)
!   END do
!    close(78)
!  END if
END program seg
