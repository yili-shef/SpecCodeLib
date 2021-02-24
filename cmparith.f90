module mcmparith
  implicit none

contains

  complex function add(a,b)
    implicit none
  
    complex, intent(in) :: a, b
  
    add=cmplx(real(a)+real(b), aimag(a)+aimag(b))

  end function add

  complex function add(a,b)
    implicit none

    complex, dimension(:) :: add
    complex, dimension(:), intent(in) :: a, b

    add=cmplx(real(a)+real(b), aimag(a)+aimag(b))
  end function add

  complex function minus(a,b)
    implicit none
  
    complex, intent(in) :: a, b

    minus = cmplx(real(a)-real(b), aimag(a)-aimag(b))

  end function minus

  complex function prod(a,b)
    implicit none

    complex, intent(in) :: a, b

    prod = cmplx(real(a)*real(b),aimag(a)*aimag(b))

  end function prod

  complex function div(a,b)
    implicit none
  
    complex, intent(in) :: a, b
  
    div = cmplx(real(a)/real(b), aimag(a)/aimag(b))

  end function div

end module mcmparith
