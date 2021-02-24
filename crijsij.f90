program checkrs
  use mconstant
  implicit none

  integer, parameter :: nx=256, lx=nx/2,neta=80
  real,    parameter :: rnu=0.0015,eps=0.1,eta=(rnu**3/eps)**.25
  real,    parameter :: alpha=2.

  real :: delta,g,rs,hk
  integer :: ll,ndel,kk,ii

  delta=neta*eta
  ndel=floor(log(delta*nx/(2.*pi))/log(alpha))
  do ll=0,ndel
    open(15,file='helispec5.dat')

      rs=0.
      do ii=1,lx
        read(15,*) kk,hk
        g=exp(-kk**2*delta**2/24)
        rs=rs+.5*kk**2*hk*g*g
      end do

      write(*,*) rs
      delta=delta/alpha

    close(15)
  end do
end program checkrs
