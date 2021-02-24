subroutine pre_interp(yp, dxyz, bg, lhnode, norder)
    use mconstant
    implicit none

    integer, intent(in) :: norder
    real(sp), intent(in) :: yp(3), dxyz(3)
    real(sp), intent(out) :: bg(norder, 3)
    integer, intent(out) :: lhnode(3)

    select case(norder)
    case (2)
        call pre_interp2(yp, dxyz, bg, lhnode)
    case (4)
        call pre_interp4(yp, dxyz, bg, lhnode)
    case (6)
        call pre_interp6(yp, dxyz, bg, lhnode)
    case (8)
        call pre_interp8(yp, dxyz, bg, lhnode)
    end select
end subroutine pre_interp

subroutine value(u, up, lhnode, bg, n1, n2, n3,norder)
    use mconstant
    implicit none

    integer,  intent(in)  :: lhnode(3),n1,n2,n3, norder
    real(sp), intent(in)  :: u(n1,n2,n3),bg(norder,3)
    real(sp), intent(out) :: up
 
    select case(norder)
    case (2)
        call value2(u, up, lhnode, bg, n1, n2, n3)
    case (4)
        call value4(u, up, lhnode, bg, n1, n2, n3)
    case (6)
        call value6(u, up, lhnode, bg, n1, n2, n3)
    case (8)
        call value8(u, up, lhnode, bg, n1, n2, n3)
    end select

end subroutine value

subroutine pre_interp2(yp, dxyz, bg, lhnode)
    use mconstant
    implicit none

    real(sp), intent(in)  :: yp(3), dxyz(3) 
    real(sp), intent(out) :: bg(2,3)
    integer,  intent(out) :: lhnode(3)

    real(sp) :: z
    integer :: ix, node

    do ix = 1, 3
        node = int( floor( yp(ix)/dxyz(ix) ) )
        lhnode(ix)=node
       
        z= ( yp(ix)-float(node)*dxyz(ix) )/dxyz(ix)

        bg(1,ix) = 1 - z
        bg(2,ix) = z
    end do 

end subroutine pre_interp2

subroutine value2(u, up, lhnode, bg, n1, n2, n3)
    use mconstant
    implicit none

    integer,  intent(in)  :: lhnode(3),n1,n2,n3
    real(sp), intent(in)  :: u(n1,n2,n3),bg(2,3)
    real(sp), intent(out) :: up
 
    integer :: ix(2),iy(2),iz(2),ndx,ndy,ndz,ix1,iy1,iz1,jx,jy,jz

    jx=lhnode(1)-1
    ix(1)=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1
    if( ix(1) .le. n1 - 1) then
        ix(2) = ix(1) + 1
    else
        ix(2) = 1
    end if

    jy=lhnode(2)-1
    iy(1)=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1
    if( iy(1) .le. n2 - 1) then
        iy(2) = iy(1) + 1
    else
        iy(2) = 1
    end if

    jz=lhnode(3)-1
    iz(1)=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1
    if( iz(1) .le. n3 - 1) then
        iz(2) = iz(1) + 1
    else
        iz(2) = 1
    end if

    up=0.0
    do ndx=1,2
      ix1 = ix(ndx)
      do ndy=1,2
        iy1 = iy(ndy)
        do ndz=1,2
          iz1 = iz(ndz)
          up = up + u(ix1,iy1,iz1)*bg(ndx,1)* &
                       bg(ndy,2)*bg(ndz,3)
        end do
      end do
    end do

end subroutine value2

subroutine pre_interp4(yp, dxyz, bg, lhnode)

  use mconstant
  implicit none

  real(sp) :: yp(3),dxyz(3),bg(4,3)
  integer  :: lhnode(3)

  real(sp) :: z,z2,z3
  integer  :: ix,node

  do ix=1,3
      node=int( floor( yp(ix)/dxyz(ix) ) )
      lhnode(ix)=node

      z= ( yp(ix)-float(node)*dxyz(ix) )/dxyz(ix)
      z2=z*z
      z3=z2*z
      bg(1,ix)=( 2*z - 3*z2 +z3 )/6
      bg(2,ix)=( 2 - z - 2*z2 + z3 )/2
      bg(3,ix)=-( -2*z - z2 + z3 )/2 
      bg(4,ix)=( z3 - z)/6
  end do

  return

end subroutine pre_interp4

subroutine value4(u, up, lhnode, bg, n1, n2, n3)
  use mconstant
  implicit none

  integer,  intent(in)  :: lhnode(3),n1,n2,n3
  real(sp), intent(in)  :: u(n1,n2,n3),bg(4,3)
  real(sp), intent(out) :: up

  integer :: ix(4),iy(4),iz(4),ndx,ndy,ndz,ix1,iy1,iz1,jx,jy,jz

  up=0.0

  jx=lhnode(1)-1
  ix(1)=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1

  if(ix(1).le.(n1-3) )then
  ix(2) = ix(1) + 1
  ix(3) = ix(2) + 1
  ix(4) = ix(3) + 1
  else
  ix(2) = mod(ix(1),n1) + 1
  ix(3) = mod(ix(2),n1) + 1
  ix(4) = mod(ix(3),n1) + 1
  endif

  jy=lhnode(2)-1
  iy(1)=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1

  if(iy(1).le.(n2-3) )then
  iy(2) = iy(1) + 1
  iy(3) = iy(2) + 1
  iy(4) = iy(3) + 1
  else
  iy(2) = mod(iy(1),n2) + 1
  iy(3) = mod(iy(2),n2) + 1
  iy(4) = mod(iy(3),n2) + 1
  endif


  jz=lhnode(3)-1
  iz(1)=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1

  if(iz(1).le.(n3-3) )then
  iz(2) = iz(1) + 1
  iz(3) = iz(2) + 1
  iz(4) = iz(3) + 1
  else
  iz(2) = mod(iz(1),n3) + 1
  iz(3) = mod(iz(2),n3) + 1
  iz(4) = mod(iz(3),n3) + 1
  endif

  do  ndx=1,4
    ix1 = ix(ndx)
    do  ndy=1,4
      iy1 = iy(ndy)
      do  ndz=1,4
        iz1 = iz(ndz)
        up=up+u(ix1,iy1,iz1)*bg(ndx,1)* &
                     bg(ndy,2)*bg(ndz,3)
      end do
    end do
  end do

end subroutine value4

subroutine pre_interp6(yp,dxyz,bg,lhnode)

  use mconstant
  implicit none

  real(sp) :: yp(3),dxyz(3),bg(6,3)
  integer  :: lhnode(3)

  real(sp) :: z,z2,z3,z4,z5
  integer  :: ix,node

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

end subroutine pre_interp6
  
subroutine value6(u,up,lhnode,bg,n1,n2,n3)
  use mconstant
  implicit none

  integer,  intent(in)  :: lhnode(3),n1,n2,n3
  real(sp), intent(in)  :: u(n1,n2,n3),bg(6,3)
  real(sp), intent(out) :: up

  integer :: ix(6),iy(6),iz(6),ndx,ndy,ndz,ix1,iy1,iz1,jx,jy,jz

  up=0.0

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
end subroutine value6

subroutine pre_interp8(yp,dxyz,bg,lhnode)
  use mconstant

  real(sp) :: yp(3),dxyz(3),bg(8,3)
  integer  :: lhnode(3)

  real(sp) :: z,z2,z3,z4,z5,z6,z7
  integer  :: ix,node

  do ix=1,3
      node=int( floor(yp(ix)/dxyz(ix) ))
      lhnode(ix)=node

      z= ( yp(ix)-float(node)*dxyz(ix) )/dxyz(ix)
      z2=z*z
      z3=z2*z
      z4=z3*z
      z5=z4*z
      z6=z5*z
      z7=z6*z
      ! lagrangian interpolants, for evenly distributed points.
      bg(1,ix)=-z*(z6-7.*z5+7.*z4+35.*z3-56.*z2-28.*z+48)/5040.
      bg(2,ix)=z*(z6-6.*z5-2.*z4+60.*z3-71.*z2-54.*z+72.)/720.
      bg(3,ix)=-z*(z6-5.*z5-9.*z4+65.*z3-16.*z2-180.*z+144)/240.
      bg(4,ix)=(z7-4.*z6-14.*z5+56.*z4+49.*z3-196.*z2-36.*z+144.)/144.
      bg(5,ix)=-z*(z6-3.*z5-17.*z4+39.*z3+88.*z2-108.*z-144.)/144.
      bg(6,ix)=z*(z6-2.*z5-18.*z4+20.*z3+89.*z2-18.*z-72)/240.
      bg(7,ix)=-z*(z6-z5-17.*z4+5.*z3+64.*z2-4.*z-48.)/720.
      bg(8,ix)=z*(z6-14.*z4+49.*z2-36)/5040.
  end do
  return
end subroutine pre_interp8
          
!***!****************************************************

subroutine value8(u,up,lhnode,bg,n1,n2,n3)
  use mconstant
  implicit none

  integer,  intent(in)  :: lhnode(3),n1,n2,n3
  real(sp), intent(in)  :: u(n1,n2,n3),bg(8,3)
  real(sp), intent(out) :: up

  integer :: ix(8),iy(8),iz(8),ndx,ndy,ndz,ix1,iy1,iz1,jx,jy,jz

  up=0.0

  jx=lhnode(1)-4
  ix(1)=( mod(jx,n1) + n1*(1-isign(1,jx))/2) + 1

  if(ix(1).le.(n1-7) )then
  ix(2) = ix(1) + 1
  ix(3) = ix(2) + 1
  ix(4) = ix(3) + 1
  ix(5) = ix(4) + 1
  ix(6) = ix(5) + 1
  ix(7) = ix(6) + 1
  ix(8) = ix(7) + 1
  else
  ix(2) = mod(ix(1),n1) + 1
  ix(3) = mod(ix(2),n1) + 1
  ix(4) = mod(ix(3),n1) + 1
  ix(5) = mod(ix(4),n1) + 1
  ix(6) = mod(ix(5),n1) + 1
  ix(7) = mod(ix(6),n1) + 1
  ix(8) = mod(ix(7),n1) + 1
  endif

  jy=lhnode(2)-4
  iy(1)=( mod(jy,n2) + n2*(1-isign(1,jy))/2) + 1

  if(iy(1).le.(n2-7) )then
  iy(2) = iy(1) + 1
  iy(3) = iy(2) + 1
  iy(4) = iy(3) + 1
  iy(5) = iy(4) + 1
  iy(6) = iy(5) + 1
  iy(7) = iy(6) + 1
  iy(8) = iy(7) + 1
  else
  iy(2) = mod(iy(1),n2) + 1
  iy(3) = mod(iy(2),n2) + 1
  iy(4) = mod(iy(3),n2) + 1
  iy(5) = mod(iy(4),n2) + 1
  iy(6) = mod(iy(5),n2) + 1
  iy(7) = mod(iy(6),n2) + 1
  iy(8) = mod(iy(7),n2) + 1
  endif


  jz=lhnode(3)-4
  iz(1)=( mod(jz,n3) + n3*(1-isign(1,jz))/2) + 1

  if(iz(1).le.(n3-7) )then
  iz(2) = iz(1) + 1
  iz(3) = iz(2) + 1
  iz(4) = iz(3) + 1
  iz(5) = iz(4) + 1
  iz(6) = iz(5) + 1
  iz(7) = iz(6) + 1
  iz(8) = iz(7) + 1
  else
  iz(2) = mod(iz(1),n3) + 1
  iz(3) = mod(iz(2),n3) + 1
  iz(4) = mod(iz(3),n3) + 1
  iz(5) = mod(iz(4),n3) + 1
  iz(6) = mod(iz(5),n3) + 1
  iz(7) = mod(iz(6),n3) + 1
  iz(8) = mod(iz(7),n3) + 1
  endif


  do  ndx=1,8
  ix1 = ix(ndx)
    do  ndy=1,8
    iy1 = iy(ndy)
      do  ndz=1,8
        iz1 = iz(ndz)
        up=up+u(ix1,iy1,iz1)*bg(ndx,1)*bg(ndy,2)*bg(ndz,3)
      end do
    end do
  end do

  return
end subroutine value8

