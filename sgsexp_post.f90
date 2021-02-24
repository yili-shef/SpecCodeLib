program sgsexp_post
  use mconstant
  use mfftwplan3d
  implicit none

  integer, parameter :: nx=128,ny=128,nz=128, lx=nx/2, lx1=lx+1, ly=ny, lz=nz

  complex(sp), dimension (lx1,ly,lz) :: ux,uy,uz  
  real(sp), dimension (lx1,ly,lz) :: kx,ky,kz,k2

  integer :: i, ifile, startfno, numfile, id1, id2, id3, id4

  real(sp), dimension(lx1,ly,lz) :: atmp
  real(sp), dimension(lx) :: Ek, aEk

  startfno=6
  numfile=13

  CALL wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  aEk=0._SP
  DO ifile=startfno,startfno+numfile-1
    WRITE(*,*) ifile
    id1=MOD(ifile,10)+48
    id2=MOD(INT(ifile/10),10)+48
    id3=MOD(INT(ifile/100),10)+48
    id4=MOD(INT(ifile/1000),10)+48
    OPEN(16,FILE='./out/vel'//CHAR(id4)//CHAR(id3)//CHAR(id2)//CHAR(id1)//'.dat', &
         FORM='unformatted')
      READ(16) ux
      READ(16) uy
      READ(16) uz
    CLOSE(16)

    atmp=real(2._SP*(.5_SP*(ux*CONJG(ux)+uy*CONJG(uy)+uz*CONJG(uz))))
    atmp(1,:,:)=.5_SP*atmp(1,:,:) 
    DO i=1,lx1-1
       Ek(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT.0.5_SP))
    END DO

    aEk=aEk+Ek
  END DO
  aEk=aEk/13.

  open(16, file='aEk.dat')
  do i=1,lx1-1
    write(16,*) i, aEk(i)
  end do
  close(16)
  
end program sgsexp_post      
