MODULE mfftwplan3d
  IMPLICIT NONE
  include 'fftw_f77.i'     

  INTEGER :: R2C3D, C2R3D, R2C3DB, C2R3DB
CONTAINS
  SUBROUTINE fftwplan3d(ix,iy,iz)

    INTEGER :: ix,iy,iz,ixb,iyb,izb
    
    ixb=3*ix/2
    iyb=3*iy/2
    izb=3*iz/2
    call rfftw3d_f77_create_plan(R2C3D,ix,iy,iz,FFTW_REAL_TO_COMPLEX, &
    FFTW_MEASURE+FFTW_IN_PLACE)
    call rfftw3d_f77_create_plan(C2R3D,ix,iy,iz,FFTW_COMPLEX_TO_REAL, &
    FFTW_MEASURE+FFTW_IN_PLACE)
    call rfftw3d_f77_create_plan(R2C3DB,ixb,iyb,izb,FFTW_REAL_TO_COMPLEX, &
    FFTW_MEASURE+FFTW_IN_PLACE)
    call rfftw3d_f77_create_plan(C2R3DB,ixb,iyb,izb,FFTW_COMPLEX_TO_REAL, &
    FFTW_MEASURE+FFTW_IN_PLACE)
  END SUBROUTINE fftwplan3d

  SUBROUTINE destroyplan3d
    call rfftwnd_f77_destroy_plan(R2C3D)
    call rfftwnd_f77_destroy_plan(C2R3D)
    call rfftwnd_f77_destroy_plan(R2C3DB)
    call rfftwnd_f77_destroy_plan(C2R3DB)
  END SUBROUTINE destroyplan3d
END MODULE mfftwplan3d
