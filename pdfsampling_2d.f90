SUBROUTINE pdfsampling_2d (rawdatax,translationx,scalingx,minvaluex,maxvaluex,npntx,xpdf, &
                           rawdatay,translationy,scalingy,minvaluey,maxvaluey,npnty,ypdf, &
                           nrawdata,zpdf,check)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: npntx,npnty,nrawdata
  REAL(SP), INTENT(IN) :: translationx,scalingx,minvaluex,maxvaluex
  REAL(SP), INTENT(IN) :: translationy,scalingy,minvaluey,maxvaluey
  REAL(SP), INTENT(OUT) :: check
  
  REAL(SP), DIMENSION(nrawdata), INTENT(IN)  :: rawdatax
  REAL(SP), DIMENSION(nrawdata), INTENT(IN)  :: rawdatay
  REAL(SP), DIMENSION(npntx),     INTENT(OUT) :: xpdf
  REAL(SP), DIMENSION(npnty),     INTENT(OUT) :: ypdf
  REAL(SP), DIMENSION(npntx,npnty), INTENT(OUT) :: zpdf
  
  REAL(SP), DIMENSION(nrawdata) :: normeddatax
  REAL(SP), DIMENSION(nrawdata) :: normeddatay
  REAL(SP) :: binwidthx,binwidthy,ttx,tty
  INTEGER  :: i,ii,iii,iiii

  binwidthx=(maxvaluex-minvaluex)/REAL(npntx,SP)
  xpdf=(/(minvaluex+(ii-.5_SP)*binwidthx, ii=1,npntx)/)
  binwidthy=(maxvaluey-minvaluey)/REAL(npnty,SP)
  ypdf=(/(minvaluey+(ii-.5_SP)*binwidthy, ii=1,npnty)/)

  normeddatax=(rawdatax-translationx)/(scalingx+TINY)-minvaluex
  normeddatay=(rawdatay-translationy)/(scalingy+TINY)-minvaluey

  ttx=maxvaluex-minvaluex
  tty=maxvaluey-minvaluey
  zpdf=0._SP
  DO i=1,nrawdata
    IF (normeddatax(i) .LT. 0._SP) CYCLE
    IF (normeddatax(i) .GT. ttx) CYCLE
    IF (normeddatay(i) .LT. 0._SP) CYCLE
    IF (normeddatay(i) .GT. tty) CYCLE
    ii=FLOOR(normeddatax(i)/binwidthx)+1
    iii=FLOOR(normeddatay(i)/binwidthy)+1
    zpdf(ii,iii)=zpdf(ii,iii)+1._SP
  END DO
  zpdf=zpdf/REAL(nrawdata,SP)
  zpdf=zpdf/binwidthx/binwidthy
  check=SUM(zpdf)*binwidthx*binwidthy
  
END SUBROUTINE pdfsampling_2d
