! Combine the right hand side of the Navier-Stokes equation with the forcing
! term. When only the forcing term is needed, set wx wy wz =0 upon input.
subroutine velforce2c(va,vb,wa,wb,k2,lx1,ly,lz,fkmax,eps)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fkmax,eps
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wa,wb
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: va,vb
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2
  
  real(sp) :: kinet, a, tmp 
  integer :: ii, jj, kk

  kinet = 0.
  a = fkmax*fkmax
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    if ( k2(ii,jj,kk) .lt. a) then 
        tmp = va(ii,jj,kk) * conjg(va(ii,jj,kk)) + &
              vb(ii,jj,kk) * conjg(vb(ii,jj,kk))
        if ( ii .eq. 1) then
            kinet = kinet + tmp/2.
        else
            kinet = kinet + tmp
        end if 
    end if
  end do
  end do
  end do

  tmp=eps/(2._sp*kinet)
  where ( k2 .lt. a) 
      wa = wa + tmp * va
      wb = wb + tmp * vb
  endwhere

end subroutine velforce2c
