module mconstant
  implicit none

  integer,     parameter :: sp=selected_real_kind( precision(1.d0) )
  integer,     parameter :: dp=selected_real_kind( precision(1.d0) )


  real(dp),    parameter :: myeps    = 2._dp * epsilon(1._dp)
  real(dp),    parameter :: myeps_dp = 2._dp * epsilon(1._dp)
  real(sp),    parameter :: myeps_sp = 2._sp * epsilon(1._sp)
  real(dp),    parameter :: mytiny    = tiny(1.0_dp)
  real(dp),    parameter :: mytiny_dp = tiny(1.0_dp)
  real(sp),    parameter :: mytiny_sp = tiny(1.0_sp)


  complex(dp), parameter :: eye = (0._dp,1._dp)
  real(dp),    parameter :: pi = 3.1415926536_dp
  
  real(sp),    parameter :: smallest=2._sp*epsilon(1._sp)
  real(sp),    parameter :: oneless=1._sp-smallest

end module

