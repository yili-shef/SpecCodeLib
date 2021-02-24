      implicit real*4(a-h,o-z)
      include 'params.inc'
      parameter(n1pp=n1+2)
      dimension dxyz(3),yp(3),bg(6,3),  &
                lhnode(3),u(n1pp,n2,n3),up(nn1,nn2,nn3), &
                x(nn1),y(nn2),z(nn3)
      pi=3.14159265358979323846
      pi2=pi*2.0
      hx=pi2
      hy=pi2
      hz=pi2
      dxyz(1)=hx/real(n1)
      dxyz(2)=hy/real(n2)
      dxyz(3)=hz/real(n3)
      dx=hx/real(nn1)
      dy=hy/real(nn2)
      dz=hz/real(nn3)

      do i=1,nn1
      x(i)=real(i-1)*dx
      end do

      do j=1,nn2
      y(j)=real(j-1)*dy
      end do

      do k=1,nn3
      z(k)=real(k-1)*dz
      end do
      open(2,file='vort00002.dat',form='formatted')
      u=0.0
      do k=1,n3
        do j=1,n2
          do i=1,n1
            read(2,100)a,b,c,u(i,j,k)
          end do
        end do
      end do

      do i=1,nn1
      do j=1,nn2
      do k=1,nn3
        yp(1)=x(i)
        yp(2)=y(j)
        yp(3)=z(k)

        call pre_interp(yp,dxyz,bg,lhnode)
        call value(u,up(i,j,k),lhnode,bg)
      end do
      end do
!      write(*,*)u(i,10,10),up(i,10,10)
      end do

      open(3,file='vort-200.dat',form='formatted')
      write(3,*)'VARIABLES="X","Y","Z","VORT"'
      write(3,*)'ZONE',' i=',nn1,' j=',nn2,' k=',nn3
      do k=1,nn3
        do j=1,nn2
          do i=1,nn1
            write(3,100)x(i),y(j),z(k),up(i,j,k)
          end do
        end do
      end do
100   format(4(e12.6,1x))
      end 
      
!*********************************************************

      subroutine pre_interp(yp,dxyz,bg,lhnode)
      implicit real*4(a-h,o-z)
      include 'params.inc'
      dimension yp(3),dxyz(3),bg(6,3),lhnode(3)        
      do ix=1,3
          node=int( yp(ix)/dxyz(ix) )
          if(yp(ix).lt.0.0)node=node-1
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

      subroutine value(u,up,lhnode,bg)
!
      implicit real*4(a-h,o-z)
      include 'params.inc'
      parameter ( n1pp=n1+2)
      real up
      dimension u(n1pp,n2,n3)
      dimension bg(6,3), lhnode(3)
      dimension ix(6),iy(6),iz(6)
!
      up=0.0
!
      mn1=n1
      mn2=n2
      mn3=n3
              
      jx=lhnode(1)-2
      ix(1)=( mod(jx,mn1) + mn1*(1-isign(1,jx))/2) + 1

      if(ix(1).le.(mn1-5) )then
      ix(2) = ix(1) + 1
      ix(3) = ix(2) + 1
      ix(4) = ix(3) + 1
      ix(5) = ix(4) + 1
      ix(6) = ix(5) + 1
      else
      ix(2) = mod(ix(1),mn1) + 1
      ix(3) = mod(ix(2),mn1) + 1
      ix(4) = mod(ix(3),mn1) + 1
      ix(5) = mod(ix(4),mn1) + 1
      ix(6) = mod(ix(5),mn1) + 1
      endif

      jy=lhnode(2)-2
      iy(1)=( mod(jy,mn2) + mn2*(1-isign(1,jy))/2) + 1

      if(iy(1).le.(mn2-5) )then
      iy(2) = iy(1) + 1
      iy(3) = iy(2) + 1
      iy(4) = iy(3) + 1
      iy(5) = iy(4) + 1
      iy(6) = iy(5) + 1
      else
      iy(2) = mod(iy(1),mn2) + 1
      iy(3) = mod(iy(2),mn2) + 1
      iy(4) = mod(iy(3),mn2) + 1
      iy(5) = mod(iy(4),mn2) + 1
      iy(6) = mod(iy(5),mn2) + 1
      endif

      jz=lhnode(3)-2
      iz(1)=( mod(jz,mn3) + mn3*(1-isign(1,jz))/2) + 1

      if(iz(1).le.(mn3-5) ) then
      iz(2) = iz(1) + 1
      iz(3) = iz(2) + 1
      iz(4) = iz(3) + 1
      iz(5) = iz(4) + 1
      iz(6) = iz(5) + 1
      else
      iz(2) = mod(iz(1),mn3) + 1
      iz(3) = mod(iz(2),mn3) + 1
      iz(4) = mod(iz(3),mn3) + 1
      iz(5) = mod(iz(4),mn3) + 1
      iz(6) = mod(iz(5),mn3) + 1
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
