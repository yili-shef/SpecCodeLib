!==================================================================
SUBROUTINE symmpart (vx,vy,vz,sv11,sv12,sv13,sv22,sv23,sv33,lx,ly,lz,kx,ky,kz)
  USE mconstant, ONLY : SP
  IMPLICIT NONE

  INTEGER :: lx,ly,lz
  REAL(SP), DIMENSION(lx,ly,lz) :: kx,ky,kz
  COMPLEX(SP), DIMENSION(lx,ly,lz) :: vx,vy,vz,sv11,sv12,sv13,sv22,sv23,sv33
  COMPLEX(SP), PARAMETER :: eye=(0.,1.)

  sv11=eye*kx*vx
  sv22=eye*ky*vy
  sv33=eye*kz*vz
  sv12=.5_SP*eye*(ky*vx+kx*vy)
  sv13=.5_SP*eye*(kz*vx+kx*vz)
  sv23=.5_SP*eye*(ky*vz+kz*vy)

END SUBROUTINE symmpart
