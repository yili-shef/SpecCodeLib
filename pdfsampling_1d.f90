SUBROUTINE pdfsampling_1d (rawdata,translation,scaling,minvalue,maxvalue,npnt,ypdf,xpdf,nrawdata)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: npnt,nrawdata
  REAL(SP), INTENT(IN) :: translation,scaling,minvalue,maxvalue
  
  REAL(SP), DIMENSION(nrawdata), INTENT(IN)  :: rawdata
  REAL(SP), DIMENSION(npnt),     INTENT(OUT) :: xpdf,ypdf
  
  REAL(SP), DIMENSION(nrawdata) :: normeddata
  REAL(SP) :: binwidth
  INTEGER  :: ii,i

  binwidth=(maxvalue-minvalue)/REAL(npnt,SP)
  xpdf=(/(minvalue+(ii-.5_SP)*binwidth, ii=1,npnt)/)
  normeddata=(rawdata-translation)/(scaling+TINY)

  normeddata=normeddata-minvalue
  WHERE(normeddata .LE. 0._SP)
          normeddata=.5_SP*binwidth
  ELSE WHERE(normeddata .GE. maxvalue-minvalue)
          normeddata=maxvalue-minvalue-0.5_SP*binwidth
  END WHERE

  ypdf=0._SP
  DO i=1,nrawdata
    ii=FLOOR(normeddata(i)/binwidth)+1
    ypdf(ii)=ypdf(ii)+1._SP
  END DO
  ypdf=ypdf/REAL(nrawdata,SP)/binwidth
  
END SUBROUTINE pdfsampling_1d
