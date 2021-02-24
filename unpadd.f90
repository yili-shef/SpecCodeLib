SUBROUTINE unpadd (ub,u,lx1,ly,lz,lxb1,lyb,lzb)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz,lxb1,lyb,lzb
  COMPLEX(SP), DIMENSION(lxb1,lyb,lzb), INTENT(IN)  :: ub
  COMPLEX(SP), DIMENSION(lx1,ly,lz),    INTENT(OUT) :: u

  u(1:lx1,1:ly/2,1:lz/2)=ub(1:lx1,1:ly/2,1:lz/2)
  u(1:lx1,1:ly/2,lz/2+1:lz)=ub(1:lx1,1:ly/2,lzb-lz/2+1:lzb)
  u(1:lx1,ly/2+1:ly,1:lz/2)=ub(1:lx1,lyb-ly/2+1:lyb,1:lz/2)
  u(1:lx1,ly/2+1:ly,lz/2+1:lz)=ub(1:lx1,lyb-ly/2+1:lyb,lzb-lz/2+1:lzb)
END SUBROUTINE unpadd
