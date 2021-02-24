module miderivs
  interface iderivs
    subroutine derivs(x,y,dydx,ny)
      use mconstant
      implicit none
      integer,  intent(in) :: ny
      real(dp), intent(in) :: x
      real(dp), dimension(ny), intent(in) :: y
      real(dp), dimension(ny), intent(out) :: dydx
    end subroutine derivs
  end interface iderivs
end module miderivs

subroutine rkcaller(derivs, ystart,ny,x1,x2,eps,h1,hmin,hnext,nok,nbad)
  use mconstant
  implicit none
  !use miderivs
  interface iderivs
    subroutine derivs(x,y,dydx,ny)
      use mconstant
      implicit none
      integer,  intent(in) :: ny
      real(dp), intent(in) :: x
      real(dp), dimension(ny), intent(in) :: y
      real(dp), dimension(ny), intent(out) :: dydx
    end subroutine derivs
  end interface iderivs

  integer,  intent(in) :: ny
  real(dp), dimension(ny), intent(inout) :: ystart
  real(dp), intent(in) :: x1,x2,eps,h1,hmin
  real(dp), intent(out) :: hnext
  integer, intent(inout) :: nok,nbad

  integer, parameter :: maxstp=10000
  integer :: nstp
  real(dp) :: h, hdid, x
  real(dp), dimension(ny) :: dydx, y, yscal

  x=x1
  h=sign(h1,x2-x1)
  y=ystart
  do nstp=1,maxstp
    call derivs(x,y,dydx,ny)
    yscal(:)=abs(y(:))+abs(h*dydx(:))+tiny
    if ((x+h-x2)*(x+h-x1)>0.0_dp) h=x2-x
    call rkqs(derivs, y,dydx,x,ny,h,eps,yscal,hdid,hnext)
    if (hdid==h) then
            nok=nok+1
    else
            nbad=nbad+1
    end if
    if ((x-x2)*(x2-x1)>=0.0_dp) then
            ystart=y
            return
    end if
    if (abs(hnext)<hmin) then
            write(*,*) 'stepsize smaller than hmin in rkcaller. stop.'
            stop
    end if
    h=hnext
  end do
  write(*,*) 'too many steps in rkcaller. stop.'
  stop
end subroutine rkcaller


subroutine rkqs(derivs, y,dydx,x,ny,htry,eps,yscal,hdid,hnext)
  use mconstant
  implicit none
  !use miderivs
  interface iderivs
    subroutine derivs(x,y,dydx,ny)
      use mconstant
      implicit none
      integer,  intent(in) :: ny
      real(dp), intent(in) :: x
      real(dp), dimension(ny), intent(in) :: y
      real(dp), dimension(ny), intent(out) :: dydx
    end subroutine derivs
  end interface iderivs
  
  integer,  intent(in) :: ny
  real(dp), dimension(ny), intent(inout) :: y
  real(dp), dimension(ny), intent(in) :: dydx,yscal
  real(dp), intent(inout) :: x
  real(dp), intent(in) :: htry,eps
  real(dp), intent(out) :: hdid,hnext

  real(dp) :: errmax, h, htemp, xnew
  real(dp), dimension(ny) :: yerr,ytemp
  real(dp), parameter :: safety=0.9_dp, pgrow=-0.2_dp, pshrnk=-0.25_dp, &
    errcon=1.89e-4_dp

  h=htry
  do
    call rkck(derivs, y,dydx,x,ny,h,ytemp,yerr)
    errmax=maxval(abs(yerr(:)/yscal(:)))/eps
    if (errmax<=1.0_dp) exit
    htemp=safety*h*(errmax**pshrnk)
    h=sign(max(abs(htemp),0.1_dp*abs(h)),h)
    xnew=x+h
    if (xnew==x) then
            write(*,*) 'stepsize underflow in rkqs. stop'
            stop
    end if
  end do
  if (errmax > errcon) then
          hnext=safety*h*(errmax**pgrow)
  else
          hnext=5.0_dp*h
  end if
  hdid=h
  x=x+h
  y=ytemp
end subroutine rkqs


! adopted from numerical recipe. using the fifth
! order cash-karp runge-kutta method to advance 
! the solution.

subroutine rkck(derivs, y,dydx,x,ny,h,yout,yerr)
  use mconstant
  implicit none
  !use miderivs
  interface iderivs
    subroutine derivs(x,y,dydx,ny)
      use mconstant
      implicit none
      integer,  intent(in) :: ny
      real(dp), intent(in) :: x
      real(dp), dimension(ny), intent(in) :: y
      real(dp), dimension(ny), intent(out) :: dydx
    end subroutine derivs
  end interface iderivs

  integer,  intent(in) :: ny
  real(dp), dimension(ny), intent(in) :: y, dydx
  real(dp), intent(in) :: x, h
  real(dp), dimension(ny), intent(out) :: yout, yerr

  real(dp), dimension(ny) :: ak2,ak3,ak4,ak5,ak6,ytemp
  real(dp), parameter :: a2=0.2_dp, a3=0.3_dp, a4=0.6_dp, a5=1.0_dp, &
    a6=0.875_dp, b21=0.2_dp, b31=3.0_dp/40.0_dp, b32=9.0_dp/40.0_dp, &
    b41=0.3_dp, b42=-0.9_dp, b43=1.2_dp, b51=-11.0_dp/54.0_dp, &
    b52=2.5_dp, b53=-70.0_dp/27.0_dp, b54=35.0_dp/27.0_dp, &
    b61=1631.0_dp/55296.0_dp, b62=175.0_dp/512.0_dp, &
    b63=575.0_dp/13824.0_dp, b64=44275.0_dp/110592.0_dp, &
    b65=253.0_dp/4096.0_dp, c1=37.0_dp/378.0_dp, &
    c3=250.0_dp/621.0_dp, c4=125.0_dp/594.0_dp,  &
    c6=512.0_dp/1771.0_dp, dc1=c1-2825.0_dp/27648.0_dp, &
    dc3=c3-18575.0_dp/48384.0_dp, dc4=c4-13525.0_dp/55296.0_dp, &
    dc5=-277.0_dp/14336.0_dp, dc6=c6-0.25_dp

  ytemp=y+b21*h*dydx                                   !first step
  call derivs(x+a2*h,ytemp,ak2,ny)                     !second step
  ytemp=y+h*(b31*dydx+b32*ak2)
  call derivs(x+a3*h,ytemp,ak3,ny)                     !third step
  ytemp=y+h*(b41*dydx+b42*ak2+b43*ak3)
  call derivs(x+a4*h,ytemp,ak4,ny)                     !fourth step
  ytemp=y+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  call derivs(x+a5*h,ytemp,ak5,ny)                     !fifth step
  ytemp=y+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
  call derivs(x+a6*h,ytemp,ak6,ny)                     !sixth step
  yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
  yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
end subroutine rkck
