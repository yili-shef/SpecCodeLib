subroutine padd (u,ub,lx1,ly,lz,lxb1,lyb,lzb)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz,lxb1,lyb,lzb
  complex(sp), dimension(lx1,ly,lz),    intent(in)  :: u
  complex(sp), dimension(lxb1,lyb,lzb), intent(out) :: ub

  ub = (0._sp,0._sp)
  ub(1:lx1,1:ly/2,1:lz/2)=u(1:lx1,1:ly/2,1:lz/2)
  ub(1:lx1,1:ly/2,lzb-lz/2+1:lzb)=u(1:lx1,1:ly/2,lz/2+1:lz)
  ub(1:lx1,lyb-ly/2+1:lyb,1:lz/2)=u(1:lx1,ly/2+1:ly,1:lz/2)
  ub(1:lx1,lyb-ly/2+1:lyb,lzb-lz/2+1:lzb)=u(1:lx1,ly/2+1:ly,lz/2+1:lz)
end subroutine padd
