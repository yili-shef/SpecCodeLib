USE mconstant
IMPLICIT NONE

INTEGER,  PARAMETER :: lx1=33,ly=64,lz=64,iseed=-39,ire=10
REAL(SP), PARAMETER :: hrel=1._SP, akp=3._SP, u0=0.4_SP
REAL(SP),    DIMENSION(lx1,ly,lz) :: kx,ky,kz,k2,tmp,atmp
COMPLEX(SP), DIMENSION(lx1,ly,lz) :: hpx,hpy,hpz,hmx,hmy,hmz
COMPLEX(SP), DIMENSION(lx1,ly,lz) :: ux,uy,uz

INTEGER  :: i
REAL(SP) :: sk

CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
CALL heliwave(hpx,hpy,hpz,hmx,hmy,hmz,kx,ky,kz,k2,lx1,ly,lz)

atmp = u0 * SQRT(8._SP*SQRT(2._SP/pi)/(3._SP*pi*akp**5)) * SQRT(k2) * EXP (-k2/(akp*akp)) 
tmp=2._SP*atmp*atmp
tmp(1,:,:)=.5_SP*tmp(1,:,:)
WRITE(*,*) SUM(tmp), SUM(tmp)*3._SP/2._SP, SUM(tmp)/(u0*u0)
sk=0._SP
DO i=1,ire
CALL initialize_heli(hpx,hpy,hpz,hmx,hmy,hmz,ux,uy,uz,k2,lx1,ly,lz,iseed,hrel,akp,u0)
tmp=ux*CONJG(ux)+uy*CONJG(uy)+uz*CONJG(uz)
tmp(1,:,:)=.5_SP*tmp(1,:,:)
sk=sk+SUM(tmp)
END DO
WRITE(*,*) 'mean k:', sk/REAL(ire,SP), sk/REAL(ire,SP)*2._SP/3._SP,sk/REAL(ire,SP)/(3._SP/2._SP*u0*u0)
END
