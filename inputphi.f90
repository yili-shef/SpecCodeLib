subroutine inputphi(phi,idump,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: phi 

  integer      :: idump
  character*50 :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)

  fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  read(10)phi
  close(10)


  return

end subroutine inputphi
