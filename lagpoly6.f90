!*********************************************************

      SUBROUTINE pre_interp(yp,dxyz,bg,lhnode)
      USE mconstant

      REAL(SP) :: yp(3),dxyz(3),bg(6,3)
      INTEGER  :: lhnode(3)

      REAL(SP) :: z,z2,z3,z4,z5
      INTEGER  :: ix,node

      do ix=1,3
          node=int( floor( yp(ix)/dxyz(ix) ) )
          lhnode(ix)=node
 
          z= ( yp(ix)-float(node)*dxyz(ix) )/dxyz(ix)
          z2=z*z
          z3=z2*z
          z4=z3*z
          z5=z4*z
          bg(1,ix)=( 6.0*z - 5.0*z2 - 5.0*z3 + 5.0*z4 - z5 )/120.0
          bg(2,ix)=( - 12.0*z + 16.0*z2 - z3 - 4.0*z4 + z5 )/24.0
          bg(3,ix)=( 12.0 - 4.0*z - 15.0*z2 + 5.0*z3 + 3.0*z4 &
                        - z5 )/12.0
          bg(4,ix)=( 12.0*z + 8.0*z2 - 7.0*z3 - 2.0*z4 + z5 )/12.0
          bg(5,ix)=( - 6.0*z - z2 + 7.0*z3 + z4 - z5 )/24.0
          bg(6,ix)=( 4.0*z - 5.0*z3 + z5 )/120.0
      end do
      return
      end
      
!****************************************************

      subroutine value(u,up,lhnode,bg,n1,n2,n3)
      USE mconstant
      IMPLICIT NONE
      INTEGER,  INTENT(IN)  :: lhnode(3),n1,n2,n3
      REAL(SP), INTENT(IN)  :: u(n1,n2,n3),bg(6,3)
      REAL(SP), INTENT(OUT) :: up

      INTEGER :: ix(6),iy(6),iz(6),ndx,ndy,ndz,ix1,iy1,iz1,jx,jy,jz
!
      up=0.0
!
      jx=lhnode(1)-2
      ix(1)=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1

      if(ix(1).le.(n1-5) )then
      ix(2) = ix(1) + 1
      ix(3) = ix(2) + 1
      ix(4) = ix(3) + 1
      ix(5) = ix(4) + 1
      ix(6) = ix(5) + 1
      else
      ix(2) = mod(ix(1),n1) + 1
      ix(3) = mod(ix(2),n1) + 1
      ix(4) = mod(ix(3),n1) + 1
      ix(5) = mod(ix(4),n1) + 1
      ix(6) = mod(ix(5),n1) + 1
      endif

      jy=lhnode(2)-2
      iy(1)=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1

      if(iy(1).le.(n2-5) )then
      iy(2) = iy(1) + 1
      iy(3) = iy(2) + 1
      iy(4) = iy(3) + 1
      iy(5) = iy(4) + 1
      iy(6) = iy(5) + 1
      else
      iy(2) = mod(iy(1),n2) + 1
      iy(3) = mod(iy(2),n2) + 1
      iy(4) = mod(iy(3),n2) + 1
      iy(5) = mod(iy(4),n2) + 1
      iy(6) = mod(iy(5),n2) + 1
      endif

      jz=lhnode(3)-2
      iz(1)=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1

      if(iz(1).le.(n3-5) ) then
      iz(2) = iz(1) + 1
      iz(3) = iz(2) + 1
      iz(4) = iz(3) + 1
      iz(5) = iz(4) + 1
      iz(6) = iz(5) + 1
      else
      iz(2) = mod(iz(1),n3) + 1
      iz(3) = mod(iz(2),n3) + 1
      iz(4) = mod(iz(3),n3) + 1
      iz(5) = mod(iz(4),n3) + 1
      iz(6) = mod(iz(5),n3) + 1
      endif
 
      do  ndx=1,6
        ix1 = ix(ndx)
        do  ndy=1,6
          iy1 = iy(ndy)
          do  ndz=1,6
            iz1 = iz(ndz)
            up=up+u(ix1,iy1,iz1)*bg(ndx,1)* &
                         bg(ndy,2)*bg(ndz,3)
          end do
        end do
      end do
      return
      end

      subroutine ixiyiz(lhnode,n1,n2,n3,ix,iy,iz)
      USE mconstant
      IMPLICIT NONE
      INTEGER,  INTENT(IN)  :: lhnode(3),n1,n2,n3

      INTEGER :: ix(6),iy(6),iz(6),jx,jy,jz

      jx=lhnode(1)-2
      ix(1)=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1

      if(ix(1).le.(n1-5) )then
      ix(2) = ix(1) + 1
      ix(3) = ix(2) + 1
      ix(4) = ix(3) + 1
      ix(5) = ix(4) + 1
      ix(6) = ix(5) + 1
      else
      ix(2) = mod(ix(1),n1) + 1
      ix(3) = mod(ix(2),n1) + 1
      ix(4) = mod(ix(3),n1) + 1
      ix(5) = mod(ix(4),n1) + 1
      ix(6) = mod(ix(5),n1) + 1
      endif

      jy=lhnode(2)-2
      iy(1)=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1

      if(iy(1).le.(n2-5) )then
      iy(2) = iy(1) + 1
      iy(3) = iy(2) + 1
      iy(4) = iy(3) + 1
      iy(5) = iy(4) + 1
      iy(6) = iy(5) + 1
      else
      iy(2) = mod(iy(1),n2) + 1
      iy(3) = mod(iy(2),n2) + 1
      iy(4) = mod(iy(3),n2) + 1
      iy(5) = mod(iy(4),n2) + 1
      iy(6) = mod(iy(5),n2) + 1
      endif

      jz=lhnode(3)-2
      iz(1)=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1

      if(iz(1).le.(n3-5) ) then
      iz(2) = iz(1) + 1
      iz(3) = iz(2) + 1
      iz(4) = iz(3) + 1
      iz(5) = iz(4) + 1
      iz(6) = iz(5) + 1
      else
      iz(2) = mod(iz(1),n3) + 1
      iz(3) = mod(iz(2),n3) + 1
      iz(4) = mod(iz(3),n3) + 1
      iz(5) = mod(iz(4),n3) + 1
      iz(6) = mod(iz(5),n3) + 1
      endif

      end subroutine ixiyiz
 
