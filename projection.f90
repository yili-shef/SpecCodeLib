! Project the vx vy vz on the plane normal to k.
subroutine projection(vx,vy,vz,kx,ky,kz,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp), dimension(lx1), intent(in) :: kx
  real(sp), dimension(ly),  intent(in) :: ky
  real(sp), dimension(lz),  intent(in) :: kz
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz

  complex(sp) :: tmp
  real(sp) :: kxii,kyjj,kzkk
  integer :: ii,jj,kk

  do kk=1,lz
    kzkk=kz(kk)
    do jj=1,ly
      kyjj=ky(jj)
      do ii=1,lx1
        kxii=kx(ii)

        tmp=kxii*vx(ii,jj,kk)+kyjj*vy(ii,jj,kk)+kzkk*vz(ii,jj,kk)
        tmp=tmp/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
        vx(ii,jj,kk)=vx(ii,jj,kk)-kxii*tmp
        vy(ii,jj,kk)=vy(ii,jj,kk)-kyjj*tmp
        vz(ii,jj,kk)=vz(ii,jj,kk)-kzkk*tmp
      end do
    end do
  end do
end subroutine projection
