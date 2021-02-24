program rmfiles
  
  integer :: filen1, filen2, ll
  character(80) :: prefix, str1, str2

  if (iargc() .ne. 3) stop 'Usage: rmfiles.x filenum1 filenum2 prefix'

  call getarg(1, str1)
  str1 = adjustl(str1)
  read(str1, '(I20)') filen1

  call getarg(2, str2)
  str2 = adjustl(str2)
  read(str2, '(I20)') filen2

  call getarg(3, prefix)
  prefix = adjustl(prefix)

  
  do ll = filen1, filen2
    write(str1, '(I20)') ll
    str1 = adjustl(str1)

    open(unit=24, file = prefix(1:len_trim(prefix))//str1(1:len_trim(str1))//'.dat')
    close(unit=24, status="delete")

  end do

end program rmfiles
