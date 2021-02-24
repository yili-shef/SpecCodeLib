module msgstij
  use mconstant
  use mfftwplan2d
  use mfftwplan3d

  interface sgstij
      module procedure rsgstauij3d, rsgstauij2d
  end interface sgstij

contains

  subroutine rsgstauij3d (vi,vj,tij,g,nx,ny,nz) 
    implicit none
   
    integer, intent(in) :: nx,ny,nz
    complex(sp), dimension(nx/2+1,ny,nz), intent(in)  :: vi, vj
    complex(sp), dimension(nx/2+1,ny,nz), intent(out) :: tij
    real(sp),    dimension(nx/2+1,ny,nz), intent(in)  :: g
  
    complex(sp), dimension(nx/2+1,ny,nz) :: ui, uj
    real(sp) :: const, ignore_me
  
    const=1._sp/real(nx*ny*nz,sp)
    ui=vi 
    uj=vj
  
    call rfftwnd_f77_one_complex_to_real(c2r3d,ui,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uj,ignore_me)
  
    tij=cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))
    call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
    tij=g*const*tij
    call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)

  
    ui=g*vi
    uj=g*vj
  
    call rfftwnd_f77_one_complex_to_real(c2r3d,ui,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uj,ignore_me)
  
    tij=tij-cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))

  end subroutine rsgstauij3d

  subroutine rsgstauij2d (vi,vj,tij,g,nx,ny) 
    implicit none
   
    integer, intent(in) :: nx,ny
    complex(sp), dimension(nx/2+1,ny), intent(in)  :: vi, vj
    complex(sp), dimension(nx/2+1,ny), intent(out) :: tij
    real(sp),    dimension(nx/2+1,ny), intent(in)  :: g
  
    complex(sp), dimension(nx/2+1,ny) :: ui, uj
    real(sp) :: const, ignore_me
  
    const=1._sp/real(nx*ny,sp)
    ui=vi 
    uj=vj
  
    call rfftwnd_f77_one_complex_to_real(c2r2d,ui,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,uj,ignore_me)
  
    tij=cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))
    call rfftwnd_f77_one_real_to_complex(r2c2d,tij,ignore_me)
    tij=g*const*tij
    call rfftwnd_f77_one_complex_to_real(c2r2d,tij,ignore_me)

  
    ui=g*vi
    uj=g*vj
  
    call rfftwnd_f77_one_complex_to_real(c2r2d,ui,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,uj,ignore_me)
  
    tij=tij-cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))

  end subroutine rsgstauij2d

end module msgstij
