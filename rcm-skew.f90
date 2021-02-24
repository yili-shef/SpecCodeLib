program rcmskew
  use mconstant
  implicit none

  real(sp) :: skewaiiperp, skewaiipara, skew

  str = 'skew-ro='//strro(1:len_trim(strro))//'-'//strreal(1:len_trim(strreal))//'reals.dat'
  open(17, file = str(1:len_trim(str)) )
  read(17, '(15E13.4)') skewaiiperp, skewaiipara

end program rcmskew
