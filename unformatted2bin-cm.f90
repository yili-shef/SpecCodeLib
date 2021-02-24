implicit none
integer, parameter :: npnt = 100000
double precision :: aijall(3,3,npnt)
integer :: ll
character(90) :: str

ll = iarg()
call getarg(1,str)
str = adjustl(str)
str = str(1:len_trim(str))//'.data'

open(16, file = str(1:len_trim(str)), form = 'unformatted')
  read(16) aijall
close(16)

open(16, file = 'bin-'//str(1:len_trim(str)), form = 'binary')
  write(16) aijall
close(16)

end 
