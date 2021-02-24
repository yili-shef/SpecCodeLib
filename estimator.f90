! Estimate the viscosity, Reynolds number, time step size, eddy turn-over time
! scale for the DNS given resolution kmax, energy injection rate eps, and
! alpha, defined as kmax times eta.

implicit none

real, parameter :: ck=1.5, pi=3.1415926, ch=1.0
integer, parameter :: q=8 ! for hyperviscosity
real :: alpha, eps, beta, eps_h
integer :: kmax,i

real :: k, re, rnu, dt, tau_l, urms, dx, eta, lambda, tt, h
real*8 :: rnuq

write(*,*) 'input kmax please:'
read(*,*) kmax
write(*,*) 'kmax = ', kmax

alpha=1.5
beta=.15  ! cfl number based on maximum u.
eps=0.1
eps_h=0.3

tt=0.      
do i=1,kmax
      tt=tt+i**(-5./3.)
end do ! tt is a constant depends only on kmax

k=ck*eps**(2./3.)*tt   ! Estimate e, energy, from Kolmogorov spectrum
urms=(k*2./3.)**.5
h=ch*eps_h*eps**(-1./3.)*tt

eta = alpha/real(kmax) ! Kolmogorov length scale.
rnu = (eta**4*eps)**(1./3.) ! viscosity
rnuq= (dble(eta)**(6*q-2)*dble(eps))**(1./3.) ! hyperviscosity

lambda = (rnu*urms**2/eps)**.5 ! lambda

re = 15**(.5)*(2./3.)*ck*tt/(alpha/kmax)**(2./3.) ! Re_lambda, depends only on alpha and kmax.

dx = pi/real(kmax)

dt = beta*dx/urms/3. ! umax is estimated as 3urms

tau_l = k/eps ! Eddy turn-over time scale

write(*,*) 'eps', eps
write(*,*) 'alpha', alpha
write(*,*) 'Reynolds no.', re
write(*,*) 'TKE', k
write(*,*) 'viscosity', rnu
write(*,*) 'time step', dt
write(*,*) 'tau_L', tau_l
write(*,*) 'u rms', urms
write(*,*) 'mean helicity', h
write(*,*) 'eta = l_K', eta
write(*,*) 'tau_K', (eta**2/eps)**(1./3.)
write(*,*) 'rnuq', rnuq

end  
