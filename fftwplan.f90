MODULE mfftwplan3d
  IMPLICIT NONE
  include 'fftw_f77.i'

  INTEGER*8 :: R2C3D, C2R3D
CONTAINS
  SUBROUTINE fftwplan3d(ix,iy,iz)

    INTEGER :: ix,iy,iz
    
    call rfftw3d_f77_create_plan(R2C3D,ix,iy,iz,FFTW_REAL_TO_COMPLEX, &
    FFTW_MEASURE+FFTW_IN_PLACE)
    call rfftw3d_f77_create_plan(C2R3D,ix,iy,iz,FFTW_COMPLEX_TO_REAL, &
    FFTW_MEASURE+FFTW_IN_PLACE)
  END SUBROUTINE fftwplan3d

  SUBROUTINE fftwplan3de(ix,iy,iz)

    INTEGER :: ix,iy,iz
    
    call rfftw3d_f77_create_plan(R2C3D,ix,iy,iz,FFTW_REAL_TO_COMPLEX, &
    FFTW_ESTIMATE+FFTW_IN_PLACE)
    call rfftw3d_f77_create_plan(C2R3D,ix,iy,iz,FFTW_COMPLEX_TO_REAL, &
    FFTW_ESTIMATE+FFTW_IN_PLACE)
  END SUBROUTINE fftwplan3de

  SUBROUTINE destroyplan3d
    call rfftwnd_f77_destroy_plan(R2C3D)
    call rfftwnd_f77_destroy_plan(C2R3D)
  END SUBROUTINE destroyplan3d
END MODULE mfftwplan3d

MODULE mfftwplan2d
  IMPLICIT NONE
  include 'fftw_f77.i'

  INTEGER*8 :: R2C2D, C2R2D
CONTAINS
  SUBROUTINE fftwplan2d(ix,iy)

    INTEGER :: ix,iy
    
    call rfftw2d_f77_create_plan(R2C2D,ix,iy,FFTW_REAL_TO_COMPLEX, &
    FFTW_MEASURE+FFTW_IN_PLACE)
    call rfftw2d_f77_create_plan(C2R2D,ix,iy,FFTW_COMPLEX_TO_REAL, &
    FFTW_MEASURE+FFTW_IN_PLACE)
  END SUBROUTINE fftwplan2d

  SUBROUTINE fftwplan2de(ix,iy)

    INTEGER :: ix,iy
    
    call rfftw2d_f77_create_plan(R2C2D,ix,iy,FFTW_REAL_TO_COMPLEX, &
    FFTW_ESTIMATE+FFTW_IN_PLACE)
    call rfftw2d_f77_create_plan(C2R2D,ix,iy,FFTW_COMPLEX_TO_REAL, &
    FFTW_ESTIMATE+FFTW_IN_PLACE)
  END SUBROUTINE fftwplan2de

  SUBROUTINE destroyplan2d
    call rfftwnd_f77_destroy_plan(R2C2D)
    call rfftwnd_f77_destroy_plan(C2R2D)
  END SUBROUTINE destroyplan2d
END MODULE mfftwplan2d

MODULE mfftwplan1d
  IMPLICIT NONE
  include 'fftw_f77.i'

  INTEGER*8 :: R2C1D, C2R1D

CONTAINS

  SUBROUTINE fftwplan1d(ix)

    INTEGER :: ix

    CALL rfftw_f77_create_plan(R2C1D,ix,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
    CALL rfftw_f77_create_plan(C2R1D,ix,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE)

  END SUBROUTINE fftwplan1d
  
  SUBROUTINE destroyplan1d
    CALL rfftw_f77_destroy_plan(R2C1D)
    CALL rfftw_f77_destroy_plan(C2R1D)
  END SUBROUTINE destroyplan1d

END MODULE mfftwplan1d
