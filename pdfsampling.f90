SUBROUTINE pdfsampling (rawdata,translation,scaling,width,npnt,ypdf,xpdf,nx,ny,nz)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: nx,ny,nz,npnt
  REAL(SP), INTENT(IN) :: translation, scaling, width
  
  REAL(SP), DIMENSION(nx,ny,nz), INTENT(IN)  :: rawdata
  REAL(SP), DIMENSION(npnt),     INTENT(OUT) :: xpdf,ypdf
  
  REAL(SP), DIMENSION(nx,ny,nz) :: normeddata
  REAL(SP) :: tmp, binwidth
  INTEGER  :: ii

  binwidth=2._SP*width/REAL(npnt,SP)
  normeddata=(rawdata-translation)/scaling
  xpdf=(/(-width+(ii-.5_SP)*binwidth, ii=1,npnt)/)
  DO ii=1,npnt
  tmp=COUNT(normeddata .GE. xpdf(ii)-.5_SP*binwidth .AND. normeddata .LT. xpdf(ii)+.5_SP*binwidth)
  ypdf(ii)=tmp/REAL(nx*ny*nz,SP)/binwidth
  END DO
  
  
END SUBROUTINE pdfsampling
