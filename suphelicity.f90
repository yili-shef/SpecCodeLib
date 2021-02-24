!===================================================================
SUBROUTINE suphelicity(vx,vy,vz,kx,ky,kz,k2)

   use mprmtr
   IMPLICIT NONE
                      
   COMPLEX(SP),DIMENSION(lx1,ly,lz) :: vx,vy,vz
   REAL(SP),DIMENSION(lx1,ly,lz)    :: kx,ky,kz,tmp,k2
   REAL(SP),DIMENSION(lx1,ly,lz)    :: tmp1
   INTEGER :: is

   do is = 1, 2
       where( abs(sqrt(k2)-is-0.4999).lt.0.5 )

        tmp = sqrt((ky*REAL(SP)(vz)-kz*REAL(SP)(vy))**2+ &
             (kx*REAL(SP)(vz)-kz*REAL(SP)(vx))**2+ &
             (kx*REAL(SP)(vy)-ky*REAL(SP)(vx))**2)

        tmp1 =  sqrt(aimag(vx)**2+aimag(vy)**2+aimag(vz)**2)

        vx = cmplx( REAL(SP)(vx),  (ky*REAL(SP)(vz)-kz*REAL(SP)(vy))*tmp1/tmp)
        vy = cmplx( REAL(SP)(vy), -(kx*REAL(SP)(vz)-kz*REAL(SP)(vx))*tmp1/tmp)
        vz = cmplx( REAL(SP)(vz),  (kx*REAL(SP)(vy)-ky*REAL(SP)(vx))*tmp1/tmp)

      END where
   END do

  RETURN
END SUBROUTINE suphelicity  
