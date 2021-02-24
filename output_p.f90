SUBROUTINE output_p(p,idump,lx1,ly,lz)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN) :: p

  INTEGER       :: idump
  CHARACTER(60) :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)
  fpath='./out/p'//fnm(1:len_trim(fnm))//'.dat'

  open(10,file=fnm,status='unknown',form='unformatted')
  write(10)p
  close(10)

  RETURN
END SUBROUTINE output_p

