! The subroutine is used to calculate the energy budget from a filtered DNS
! velocity field. The difference between this and the subroutine
! energy_budget_les is that here subroutine convec_dns, other than convec, is
! invoked. The idea is that, the input velocity fields ux,uy,uz have already
! been filtered; from another perspective, the large scale velocity fields have
! been padded with zeros. Therefore, if the cutoff wavenumber is no greater than
! two thirds of the kmax of DNS, no additional padding is needed. This way, the
! definition of too large arrays is avoided.
SUBROUTINE energy_budget_fdns(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,k2,nu,TEk,PiEk,Ek,DEk,SGSDEk, &
                          lx1,ly,lz,nx,ny,nz)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: lx1,ly,lz,nx,ny,nz
  REAL(SP), INTENT(IN) :: nu
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,k2
  REAL(SP),    DIMENSION(lx1-1),     INTENT(OUT) :: TEk,Ek,DEk,PiEk,SGSDEk

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: wx,wy,wz
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: s11,s12,s13,s22,s23,s33
  REAL(SP),    DIMENSION(lx1,ly,lz) :: atmp
  INTEGER :: i

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)

  s11=eye*kx*ux
  s22=eye*ky*uy
  s33=eye*kz*uz
  s12=.5_SP*eye*(ky*ux+kx*uy)
  s13=.5_SP*eye*(kz*ux+kx*uz)
  s23=.5_SP*eye*(ky*uz+kz*uy)
 
  atmp=2._SP*(.5_SP*(ux*CONJG(ux)+uy*CONJG(uy)+uz*CONJG(uz)))
  atmp(1,:,:)=.5_SP*atmp(1,:,:) 
  DO i=1,lx1-1
     Ek(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
     DEk(i)=2._SP*nu*REAL(i,SP)*REAL(i,SP)*Ek(i)
  END DO

  atmp=-2._SP*REAL(tau11*CONJG(s11)+tau22*CONJG(s22)+tau33*CONJG(s33)+              &
                  2._SP*(tau12*CONJG(s12)+tau13*CONJG(s13)+tau23*CONJG(s23)))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    SGSDEk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT..5_SP))
  END DO
  
  CALL convec_dns(ux,uy,uz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)

  ! Note the definition of T is the nonlinear transfer for the kinetic energy of
  ! the kth mode |u(vector k)|^2/2.
  atmp=2._SP*REAL(wx*CONJG(ux)+wy*CONJG(uy)+wz*CONJG(uz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)

  DO i=lx1-1,1,-1
     TEk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
     PiEk(i)=SUM(TEk(i:lx1-1))
  END DO
END SUBROUTINE energy_budget_fdns
