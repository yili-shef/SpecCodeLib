subroutine outputvel2c(va,vb,idump,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  complex(sp), dimension(lx1,ly,lz), intent(in) :: va, vb

  integer       :: idump
  character(50) :: fnm,fpath

  write(fnm,'(i30)') idump
  fnm=adjustl(fnm)

  fpath='./out/va'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  write(10)va
  close(10)

  fpath='./out/vb'//fnm(1:len_trim(fnm))//'.dat'
  open(10,file=fpath,status='unknown',form='unformatted')
  write(10)vb
  close(10)

  return

end subroutine outputvel2c
