module minput

  interface input
      module procedure input3d, input3dwithphi, input3dphi
      module procedure input2d, input2dwithphi, input2dphi
      module procedure inputbij
      module procedure inputbyname3d
  end interface input

contains

  subroutine inputbyname3d(vx, strvx, iin, lx1, ly, lz)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly,lz
    character(*), intent(in) :: strvx
    complex(sp), dimension(lx1,ly,lz), intent(out) :: vx

    integer      :: iin
    character*50 :: fnm,fpath
  

    write(fnm,'(i30)') iin
    fnm=adjustl(fnm)
    fpath = adjustl(strvx)
  
    fpath='./out/'//fpath(1:len_trim(fpath))//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
      read(10)vx
    close(10)

  end subroutine inputbyname3d

  subroutine input2d(wz,idump,lx1,ly)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly
    complex(sp), dimension(lx1,ly), intent(out) :: wz
  
    integer      :: idump
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/wz'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)wz
    close(10)

    return

  end subroutine input2d

  subroutine input2dphi(idump,phi,lx1,ly)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly
    complex(sp), dimension(lx1,ly), intent(out) :: phi
  
    integer      :: idump
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)phi
    close(10)

    return

  end subroutine input2dphi

  subroutine input2dwithphi(wz,phi,idump,lx1,ly)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly
    complex(sp), dimension(lx1,ly), intent(out) :: wz, phi
  
    integer      :: idump
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/wz'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)wz
    close(10)

    fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)phi
    close(10)

    return

  end subroutine input2dwithphi

  subroutine input3d(ux,uy,uz,idump,lx1,ly,lz)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly,lz
    complex(sp), dimension(lx1,ly,lz), intent(out) :: ux,uy,uz
  
    integer      :: idump
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)ux
    close(10)
  
    fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uy
    close(10)
  
    fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uz
    close(10)
  
    return

  end subroutine input3d

  subroutine input3dwithphi(ux,uy,uz,phi,idump,lx1,ly,lz)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly,lz
    complex(sp), dimension(lx1,ly,lz), intent(out) :: ux,uy,uz,phi
  
    integer      :: idump
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/ux'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)ux
    close(10)
  
    fpath='./out/uy'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uy
    close(10)
  
    fpath='./out/uz'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)uz
    close(10)
  
    fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)phi
    close(10)
  
    return

  end subroutine input3dwithphi

  subroutine input3dphi(phi,idump,lx1,ly,lz)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly,lz, idump
    complex(sp), dimension(lx1,ly,lz), intent(out) :: phi
  
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/phi'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)phi
    close(10)
  
    return

  end subroutine input3dphi

  subroutine inputbij(b11,b12,b13,b21,b22,b23,b31,b32,b33,idump,lx1,ly,lz)
    use mconstant
    implicit none
  
    integer, intent(in) :: lx1,ly,lz, idump
    complex(sp), dimension(lx1,ly,lz), intent(out) :: b11, b12, b13
    complex(sp), dimension(lx1,ly,lz), intent(out) :: b21, b22, b23
    complex(sp), dimension(lx1,ly,lz), intent(out) :: b31, b32, b33
  
    character*50 :: fnm,fpath
  
    write(fnm,'(i30)') idump
    fnm=adjustl(fnm)
  
    fpath='./out/b11'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b11
    close(10)
    fpath='./out/b12'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b12
    close(10)
    fpath='./out/b13'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b13
    close(10)
    fpath='./out/b21'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b21
    close(10)
    fpath='./out/b22'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b22
    close(10)
    fpath='./out/b23'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b23
    close(10)
    fpath='./out/b31'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b31
    close(10)
    fpath='./out/b32'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b32
    close(10)
    fpath='./out/b33'//fnm(1:len_trim(fnm))//'.dat'
    open(10,file=fpath,status='unknown',form='unformatted')
    read(10)b33
    close(10)
  
    return
  end subroutine inputbij

end module minput
