SUBROUTINE outputphi(phi,idump,lx1,ly,lz)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN) :: phi

  INTEGER      :: idump
  CHARACTER*50 :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)

  fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  write(10)phi
  close(10)

  RETURN
END SUBROUTINE outputphi

