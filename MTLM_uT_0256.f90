MODULE global
real(8), parameter :: pi = 3.141592653589793238462643d0, two_Pi = 2.d0 * pi

!========  parameter for velocity field   =========
real(8), parameter :: urms = 1.000d00
real(8), parameter :: Ck = 1.50d+00
real(8), parameter :: cl = 9.811747074127197d-01   
real(8), parameter :: l = 2.073500000000000d+00
real(8), parameter :: p0 = 4.d0
real(8), parameter :: Ldom = 6.283185308d0
real(8), parameter :: phi = pi


!========  parameter for scalar field   =========
real(8), parameter :: theta_rms = 1.000d00
real(8), parameter :: qdiss = 2.893658066071859d-01
real(8), parameter :: gamma_diff = 1.955936002032939d-03
real(8), parameter :: Cq =  6.800000d-01
real(8), parameter :: clq = 2.562338113784790d-01
real(8), parameter :: etaq = 1.116071428717157d-02
real(8), parameter :: lq = 1.430000000000000d+00

integer, parameter :: nt = 256
!----------------------------------------------

integer, parameter :: nx = nt, ny = nt, nz = nt
integer, parameter :: ld = nx + 2
integer, parameter :: lh = ld/2


complex(8), dimension(nt/2 + 1,nt,nt) :: vx,vy,vz,theta


real(8) E(lh-2),Et(lh-2),eta,diss,rnu
real(8) hx,hy,hz,CFL,max_k2,resc_factor,epsilon,epsl,tiny_real
real(8) urms0,thetarms0,Int0kmax,Int0kmaxscalar
integer nxh,nyh,nzh,nyhp,nzhp,lh2,lhm,nth,nbytes
logical output
END MODULE




PROGRAM MTLM_uT
use global
implicit none

!....... definitions in fftw3.f ......................................
integer,parameter::FFTW_R2HC = 0   , FFTW_HC2R = 1    , FFTW_DHT = 2
integer,parameter::FFTW_REDFT00 = 3, FFTW_REDFT01 = 4 , FFTW_REDFT10 = 5 
integer,parameter::FFTW_REDFT11 = 6, FFTW_RODFT00 = 7 , FFTW_RODFT01 = 8 
integer,parameter::FFTW_RODFT10 = 9, FFTW_RODFT11 = 10, FFTW_FORWARD = -1

integer,parameter::FFTW_BACKWARD = +1          , FFTW_MEASURE = 0  
integer,parameter::FFTW_DESTROY_INPUT = 1      , FFTW_UNALIGNED = 2
integer,parameter::FFTW_CONSERVE_MEMORY = 4    , FFTW_EXHAUSTIVE = 8
integer,parameter::FFTW_PRESERVE_INPUT = 16    , FFTW_PATIENT = 32
integer,parameter::FFTW_ESTIMATE = 64          , FFTW_ESTIMATE_PATIENT = 128
integer,parameter::FFTW_BELIEVE_PCOST = 256    , FFTW_DFT_R2HC_ICKY = 512
integer,parameter::FFTW_NONTHREADED_ICKY = 1024, FFTW_NO_BUFFERING = 2048
integer,parameter::FFTW_NO_INDIRECT_OP = 4096  , FFTW_ALLOW_LARGE_GENERIC = 8192
integer,parameter::FFTW_NO_RANK_SPLITS = 16384 , FFTW_NO_VRANK_SPLITS = 32768
integer,parameter::FFTW_NO_VRECURSE = 65536    , FFTW_NO_SIMD = 131072
!......................................................................

!...for reduced velocity field
complex(8), dimension(nt/4 + 1 ,nt/2 ,nt/2 ) :: vx128,vy128,vz128
complex(8), dimension(nt/8 + 1 ,nt/4 ,nt/4 ) :: vx64 ,vy64 ,vz64
complex(8), dimension(nt/16 + 1,nt/8 ,nt/8 ) :: vx32 ,vy32 ,vz32
complex(8), dimension(nt/32 + 1,nt/16,nt/16) :: vx16 ,vy16 ,vz16
complex(8), dimension(nt/64 + 1,nt/32,nt/32) :: vx08 ,vy08 ,vz08

!...for reduced scalar field
complex(8), dimension(nt/4 + 1 ,nt/2 ,nt/2 ) :: theta128
complex(8), dimension(nt/8 + 1 ,nt/4 ,nt/4 ) :: theta64
complex(8), dimension(nt/16 + 1,nt/8 ,nt/8 ) :: theta32
complex(8), dimension(nt/32 + 1,nt/16,nt/16) :: theta16
complex(8), dimension(nt/64 + 1,nt/32,nt/32) :: theta08


real(8) plane(nx,ny,1),kyc(ny),kzc(nz),kxc(1)
integer ikx,iky,ikz

real(8), dimension(lh-2) :: Eloc,Etloc

integer(8) FTvx,FTvy,FTvz,IFTvx,IFTvy,IFTvz
integer(8) FTvx128,FTvy128,FTvz128,IFTvx128,IFTvy128,IFTvz128
integer(8) FTvx64,FTvy64,FTvz64,IFTvx64,IFTvy64,IFTvz64
integer(8) FTvx32,FTvy32,FTvz32,IFTvx32,IFTvy32,IFTvz32
integer(8) FTvx16,FTvy16,FTvz16,IFTvx16,IFTvy16,IFTvz16
integer(8) FTvx08,FTvy08,FTvz08,IFTvx08,IFTvy08,IFTvz08

integer(8) FTth,IFTth
integer(8) FTth128,IFTth128
integer(8) FTth64,IFTth64
integer(8) FTth32,IFTth32
integer(8) FTth16,IFTth16
integer(8) FTth08,IFTth08


complex(8) cn1,cn2,cn3,cp1,cp2,cp3
complex(8) a,b,c,t,vxx,vyy,vzz,cn

real(8) d,mag,ratio,wksp,rn,wn,t2,fact,G,const,cns,ur,u1,u2,v1,v2,w1,w2,tadv
real(8) hb,spec3d,k2r,kx,k22,norm_factor,s2,upl,scalar_spec3d,theta_rmsl
real(8) Rex,Imx,Rey,Imy,Rez,Imz,ss,tn,time_scale,total_adv,next_total
real(8) thetapl

real ran2

integer num(nt+1),ind,jx,jy,jz,idum,ans,rit,ky,kz,r,rs,nadv,ad
integer kyky,kzkz,jpy,jpz,bl,nblocks,jpx,nnr
integer jxr,jyr,jzr,jx2,jy2,jz2,nbytesre,i,j,k,ii,ipos
integer n,nit,kxi
logical ext,last_adv
character(24) cmd
character(4) gen
character(2) par2
character par,conf



!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
nxh = nx/2
nyh = ny/2
nzh = nz/2
nyhp = nyh + 1
nzhp = nzh + 1
lh2 = lh*lh
lhm = lh - 1
nth = nt/2

max_k2 = real((nt/2 + 1)**2,kind=8)

tiny_real =tiny(d)
output = .false.

eta = Ldom/(nt*phi)
diss = urms*urms*urms/l
rnu = (eta*eta*eta*eta * diss)**(1.d0/3.d0)

write(*,300) diss,rnu,nx,ny,nz


rit = 0

idum = -7993
t2 = dble(ran2(idum))    !...initialize random generator

if (mod(ld,2) /= 0) stop ' ld must be even'

hx = two_Pi/nx
hy = two_Pi/ny
hz = two_Pi/nz

write(*,"(/,' initialising FFT...')")
nnr = nt


print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx,nnr,nnr,nnr,vx,vx,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy,nnr,nnr,nnr,vy,vy,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz,nnr,nnr,nnr,vz,vz,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth,nnr,nnr,nnr,theta,theta,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx,nnr,nnr,nnr,vx,vx,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy,nnr,nnr,nnr,vy,vy,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz,nnr,nnr,nnr,vz,vz,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTth,nnr,nnr,nnr,theta,theta,FFTW_MEASURE)


nnr = nt/2
print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx128,nnr,nnr,nnr,vx128,vx128,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy128,nnr,nnr,nnr,vy128,vy128,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz128,nnr,nnr,nnr,vz128,vz128,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth128,nnr,nnr,nnr,theta128,theta128,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx128,nnr,nnr,nnr,vx128,vx128,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy128,nnr,nnr,nnr,vy128,vy128,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz128,nnr,nnr,nnr,vz128,vz128,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTth128,nnr,nnr,nnr,theta128,theta128,FFTW_MEASURE)

nnr = nt/4
print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx64,nnr,nnr,nnr,vx64,vx64,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy64,nnr,nnr,nnr,vy64,vy64,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz64,nnr,nnr,nnr,vz64,vz64,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth64,nnr,nnr,nnr,theta64,theta64,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx64,nnr,nnr,nnr,vx64,vx64,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy64,nnr,nnr,nnr,vy64,vy64,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz64,nnr,nnr,nnr,vz64,vz64,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTth64,nnr,nnr,nnr,theta64,theta64,FFTW_MEASURE)

nnr = nt/8
print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx32,nnr,nnr,nnr,vx32,vx32,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy32,nnr,nnr,nnr,vy32,vy32,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz32,nnr,nnr,nnr,vz32,vz32,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth32,nnr,nnr,nnr,theta32,theta32,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx32,nnr,nnr,nnr,vx32,vx32,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy32,nnr,nnr,nnr,vy32,vy32,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz32,nnr,nnr,nnr,vz32,vz32,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTth32,nnr,nnr,nnr,theta32,theta32,FFTW_MEASURE)

nnr = nt/16
print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx16,nnr,nnr,nnr,vx16,vx16,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy16,nnr,nnr,nnr,vy16,vy16,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz16,nnr,nnr,nnr,vz16,vz16,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth16,nnr,nnr,nnr,theta16,theta16,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx16,nnr,nnr,nnr,vx16,vx16,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy16,nnr,nnr,nnr,vy16,vy16,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz16,nnr,nnr,nnr,vz16,vz16,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTth16,nnr,nnr,nnr,theta16,theta16,FFTW_MEASURE)

nnr = nt/32
print *,nnr
call dfftw_plan_dft_c2r_3d(IFTvx08,nnr,nnr,nnr,vx08,vx08,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvy08,nnr,nnr,nnr,vy08,vy08,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTvz08,nnr,nnr,nnr,vz08,vz08,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(IFTth08,nnr,nnr,nnr,theta08,theta08,FFTW_MEASURE)

call dfftw_plan_dft_r2c_3d(FTvx08,nnr,nnr,nnr,vx08,vx08,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvy08,nnr,nnr,nnr,vy08,vy08,FFTW_MEASURE)
call dfftw_plan_dft_r2c_3d(FTvz08,nnr,nnr,nnr,vz08,vz08,FFTW_MEASURE)
call dfftw_plan_dft_c2r_3d(FTth08,nnr,nnr,nnr,theta08,theta08,FFTW_MEASURE)

999 continue



s2 = 0.5d0 * spec3d(1.0d0)
do jx=2,lh-2
  s2 = s2 + 0.5d0*(spec3d(dble(jx-1)) + spec3d(dble(jx)))
end do
Int0kmax = s2
urms0 = sqrt(2.d0*Int0kmax/3.d0)


s2 = 0.5d0 * scalar_spec3d(1.0d0)
do jx=2,lh-2
  s2 = s2 + 0.5d0*(scalar_spec3d(dble(jx-1)) + scalar_spec3d(dble(jx)))
end do
Int0kmaxscalar = s2
thetarms0 = sqrt(2.d0*Int0kmaxscalar)




! calculate num(k)
! this is the number of wavenumbers in shell k
num = 0
do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky   
do jx = 1,lhm
   kx = real(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if ((k22 < 1.d-20) .or. (k22 >= max_k2)) cycle
   ind = int(sqrt(k22) + 0.5d0) 
   if (kx < 1.d-20) then
      num(ind) = num(ind) + 1
   else
      num(ind) = num(ind) + 2
   end if
end do
end do
end do

open(unit=99, file = 'num.cnt')
write(99,"(I4,I7)") (ind,num(ind),ind=1,nt+1)
close(99)

!................... end of initialisation .....................


!...Calculate random velocity field.  abs(v(k))=sqrt(2E(k)/N(k))
do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz   
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky   
do jx = 1,nth
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if(k22 < max_k2) then
       wn = sqrt(k22)
       ind = int(wn + 0.5d0)
       if (ind < 1) ind = 1

       t2 = dble(ran2(idum))
       if (t2 <= 1.d-10) t2 = 1.d-10
       t2 = two_Pi*t2
       a = cmplx(dcos(t2),dsin(t2),kind=8)

       t2 = dble(ran2(idum))
       if (t2 <= 1.d-10) t2 = 1.d-10
       t2 = two_Pi*t2
       b = cmplx(dcos(t2),dsin(t2),kind=8)

       t2 = dble(ran2(idum))
       if (t2 <= 1.d-10) t2 = 1.d-10
       t2 = two_Pi*t2
       c = cmplx(dcos(t2),dsin(t2),kind=8)

       d = sqrt(2.d0*spec3d(wn)/dble(num(ind)))
       t = (a*kx + b*ky + c*kz)/k22
       vxx = a - t*kx
       vyy = b - t*ky
       vzz = c - t*kz
       mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))

       vx(jx,jy,jz) = vxx/mag*d     !u_i k_i=0
       vy(jx,jy,jz) = vyy/mag*d
       vz(jx,jy,jz) = vzz/mag*d
   end if
end do
end do
end do


!...first spectrum, before re-scaling

call spectrum(8,'spec001.dat')

Eloc = 0.d0
!...re-scale spectrum to improve fitting at low wave numbers

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1, nth

   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if (k22 < max_k2) then
     wn = sqrt(k22)
     ind = int(wn + 0.5d0)
     if (ind <= 1) then
       ind = 1
       wn = 1.0d0
     end if  
     if (ind < nth) then
       if (E(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(spec3d(wn)/E(ind))
       end if    
       vx(jx,jy,jz) = vx(jx,jy,jz)*ratio
       vy(jx,jy,jz) = vy(jx,jy,jz)*ratio
       vz(jx,jy,jz) = vz(jx,jy,jz)*ratio

!...compute contribution to the new energy spectrum
       Rex =  real(vx(jx,jy,jz))
       Imx = aimag(vx(jx,jy,jz))
       Rey =  real(vy(jx,jy,jz))
       Imy = aimag(vy(jx,jy,jz))
       Rez =  real(vz(jx,jy,jz))
       Imz = aimag(vz(jx,jy,jz))
       ss = Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy + Rez*Rez + Imz*Imz
       if (kx < 1.d-20) ss = ss * 0.5d0
       if (ss < tiny_real) ss = 0.d0
       Eloc(ind) = Eloc(ind) +  ss
    else
      vx(jx,jy,jz) = (0.d0,0.d0)
      vy(jx,jy,jz) = (0.d0,0.d0)
      vz(jx,jy,jz) = (0.d0,0.d0)
    end if
  else
      vx(jx,jy,jz) = (0.d0,0.d0)
      vy(jx,jy,jz) = (0.d0,0.d0)
      vz(jx,jy,jz) = (0.d0,0.d0)
  end if
end do
end do
end do


!...second rescale to match the kinetic energy
s2 = 0.5d0 * Eloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
end do
upl = sqrt(s2*2.d0/3.d0)
resc_factor = urms/upl
print *,'resc= ',resc_factor

vx = vx * resc_factor
vy = vy * resc_factor
vz = vz * resc_factor

output = .true.
call spectrum(1,'spec002.dat')
output = .false.

epsilon = epsl
open(unit=21, file='random_vxvyvz0256.bin', form='unformatted')
write(21) vx
write(21) vy
write(21) vz
close(unit=21)

!============================================================================

!...generates the random scalar field
do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz   
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky   
do jx = 1,nth
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if(k22 < max_k2) then
       wn = sqrt(k22)
       ind = int(wn + 0.5d0)
       if (ind < 1) ind = 1

       t2 = dble(ran2(idum))
       if (t2 <= 1.d-10) t2 = 1.d-10
       t2 = two_Pi*t2
       a = cmplx(dcos(t2),dsin(t2),kind=8)

       d = sqrt(2.d0*scalar_spec3d(wn)/dble(num(ind)))
       theta(jx,jy,jz) = a*d
   end if
end do
end do
end do


!...first scalar spectrum, before re-scaling

call scalar_spectrum(8,'scalar_spec001.dat')

Eloc = 0.d0
!...re-scale spectrum to improve fitting at low wave numbers

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1, nth

   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if (k22 < max_k2) then
     wn = sqrt(k22)
     ind = int(wn + 0.5d0)
     if (ind <= 1) then
       ind = 1
       wn = 1.0d0
     end if  
     if (ind < nth) then
       if (E(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(scalar_spec3d(wn)/Et(ind))
       end if    
       theta(jx,jy,jz) = theta(jx,jy,jz)*ratio

!...compute contribution to the new energy spectrum
       Rex =  real(theta(jx,jy,jz))
       Imx = aimag(theta(jx,jy,jz))
       ss = Rex*Rex + Imx*Imx
       if (kx < 1.d-20) ss = ss * 0.5d0
       if (ss < tiny_real) ss = 0.d0
       Eloc(ind) = Eloc(ind) +  ss
    else
      theta(jx,jy,jz) = (0.d0,0.d0)
    end if
  else
      theta(jx,jy,jz) = (0.d0,0.d0)
  end if
end do
end do
end do


!...second rescale to match the scalar variance
s2 = 0.5d0 * Eloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
end do

theta_rmsl = sqrt(2.d0*s2)
resc_factor = theta_rms/theta_rmsl

print *,'scalar resc= ',resc_factor

theta = theta * resc_factor

output = .true.
call scalar_spectrum(1,'scalar_spec002.dat')
output = .false.

epsilon = epsl
open(unit=22, file='random_theta0256.bin', form='unformatted')
write(22) theta
close(unit=22)

open(unit=35, file='control.out')


!...apply the MTLM procedure

call apply_MTLM(vx08,vy08,vz08,theta08,nt/32 + 2,                         &
                IFTvx08,IFTvy08,IFTvz08,IFTth08,FTvx08,FTvy08,FTvz08,FTth08)

call apply_MTLM(vx16,vy16,vz16,theta16,nt/16 + 2,                         &
                IFTvx16,IFTvy16,IFTvz16,IFTth16,FTvx16,FTvy16,FTvz16,FTth16)

call apply_MTLM(vx32,vy32,vz32,theta32,nt/8 + 2,                          &
                IFTvx32,IFTvy32,IFTvz32,IFTth32,FTvx32,FTvy32,FTvz32,FTth32)

call apply_MTLM(vx64,vy64,vz64,theta64,nt/4 + 2,                          &
                IFTvx64,IFTvy64,IFTvz64,IFTth64,FTvx64,FTvy64,FTvz64,FTth64)

call apply_MTLM(vx128,vy128,vz128,theta128,nt/2 + 2,IFTvx128,IFTvy128,    &
                IFTvz128,IFTth128,FTvx128,FTvy128,FTvz128,FTth128)



!...the last scale is direct (no filtered fields):

write(*,"('using epsilon value =',ES25.15)") epsilon



write(*,"('52/',$)") 
call dfftw_execute(IFTvx)
call dfftw_execute(IFTvy)
call dfftw_execute(IFTvz)
call dfftw_execute(IFTth)
write(*,"('53/',$)") 



!.............    compute urms   .................
ur = 0.0d0
do jz = 1,nz
do jy = 1,ny
do jx = 1,nxh
   u1 =  abs(real(vx(jx,jy,jz)))
   u2 = abs(aimag(vx(jx,jy,jz)))

   v1 =  abs(real(vy(jx,jy,jz)))
   v2 = abs(aimag(vy(jx,jy,jz)))

   w1 =  abs(real(vz(jx,jy,jz)))
   w2 = abs(aimag(vz(jx,jy,jz)))

   ur = ur + u1*u1 + v1*v1 + w1*w1 + u2*u2 + v2*v2 + w2*w2
end do
end do
end do


ur = sqrt(ur/(3.0d0*nx*ny*nz))
tn = hx/ur
time_scale = (hx*hx/epsilon)**(1.d0/3.d0)

write(*,"(/,'    urms,n = ',ES12.5)") ur
write(*,"('        tn = ',ES12.5)") tn
write(*,"('time_scale = ',ES12.5)") time_scale
write(911,"(/,'    urms,n = ',ES24.14)") ur
write(911,"('        tn = ',ES24.14)") tn
write(911,"('time_scale = ',ES24.14)") time_scale


total_adv = 0.d0
last_adv = .false.

nadv = nint(time_scale/tn)
tadv = time_scale/nadv
CFL = tadv/tn
if (CFL > 1.05) then
  nadv = nadv + 1
  tadv = time_scale/nadv
  CFL = tadv/tn
end if  
write(*,"('       CFL = ',ES12.5)") CFL
write(911,"('       CFL = ',ES12.5)") CFL  




for_time_scale: DO ad = 1,nadv
   next_total = total_adv + tadv

!... apply MTLM on disk
   call MTLM(vx,vy,vz,theta,ld,tadv,urms0,thetarms0)

!... Back to Fourier space. First 2-D FFT by x-y planes. 
   norm_factor = 1.d0/dble(nx*ny*nz)

   write(*,"('54/',$)") 
   !... back to Fourier space
   call dfftw_execute(FTvx)
   call dfftw_execute(FTvy)
   call dfftw_execute(FTvz)
   write(*,"('55/',$)") 

   vx = vx*norm_factor
   vy = vy*norm_factor
   vz = vz*norm_factor

   vx(lh,1:ny,1:nz) = (0.d0,0.d0)      !Oddballs
   vy(lh,1:ny,1:nz) = (0.d0,0.d0)
   vz(lh,1:ny,1:nz) = (0.d0,0.d0)
   vx(1:lh,nyhp,1:nz) = (0.d0,0.d0)
   vy(1:lh,nyhp,1:nz) = (0.d0,0.d0)
   vz(1:lh,nyhp,1:nz) = (0.d0,0.d0)
   vx(1:lh,1:ny,nzhp) = (0.d0,0.d0)
   vy(1:lh,1:ny,nzhp) = (0.d0,0.d0)
   vz(1:lh,1:ny,nzhp) = (0.d0,0.d0)

   vx(1,1,1) = (0.d0,0.d0)   !<u_i>=0
   vy(1,1,1) = (0.d0,0.d0)
   vz(1,1,1) = (0.d0,0.d0)


!!!!!   call spectrum(1,'specbbbb.dat') !xxxxxxx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.............. project on to divergence-free part ......................


   write(*,"('Main projection')")

   do jz = 1,nz
     if (jz == nzhp) cycle
     kz = jz - 1
     if (jz > nzh) kz = kz - nz
     kzkz = kz * kz
   do jy = 1,ny
     if (jy == nyhp) cycle
     ky = jy - 1
     if (jy > nyh) ky = ky - ny
     kyky = ky * ky
   do jx = 1,lhm
     kx = dble(jx - 1)
     k22 = kx*kx + kyky + kzkz
     if (k22 < 1.d-20) cycle

     if (k22 < max_k2) then
       cn1 = vx(jx,jy,jz)
       cn2 = vy(jx,jy,jz)
       cn3 = vz(jx,jy,jz)
       t = (cn1 * kx + cn2 * ky + cn3 * kz)/k22
       vxx = cn1 - t*kx
       vyy = cn2 - t*ky
       vzz = cn3 - t*kz
       mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))
       d = sqrt(cn1 * conjg(cn1) + cn2 * conjg(cn2) + cn3 * conjg(cn3))
       vx(jx,jy,jz) = vxx/mag*d    
       vy(jx,jy,jz) = vyy/mag*d
       vz(jx,jy,jz) = vzz/mag*d
     else
       vx(jx,jy,jz) = (0.d0,0.d0)
       vy(jx,jy,jz) = (0.d0,0.d0)
       vz(jx,jy,jz) = (0.d0,0.d0)
     end if
   end do
   end do
   end do

   total_adv = total_adv + tadv
   write(*,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv, total_adv   
   write(911,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv,total_adv 

   if (ad /= nadv) then    
!.........................................................................
!................     IFT for the next advancement .......................
      write(*,"('sIFT/',$)") 
      call dfftw_execute(IFTvx)
      call dfftw_execute(IFTvy)
      call dfftw_execute(IFTvz)
      write(*,"('eIFT',$)") 
   end if
END DO for_time_scale


!...transform the scalar field back to Fourier space
write(*,"(/,'scalar FT')") 
call dfftw_execute(FTth)

theta = theta * norm_factor

theta(lh,1:ny,1:nz) = (0.d0,0.d0)
theta(1:lh,nyhp,1:nz) = (0.d0,0.d0)
theta(1:lh,1:ny,nzhp) = (0.d0,0.d0)

theta(1,1,1) = (0.d0,0.d0)



!...calculate new spectrum for rescaling 
write(gen,"(I4.4)") nx
call spectrum(19,'afterspec'//gen//'.dat')

call scalar_spectrum(21,'safterspec'//gen//'.dat') 



Eloc = 0.d0
Etloc = 0.d0

write(*,"('Final Rescaling (1/2)')")

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1,lhm
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if (k22 < max_k2) then
     wn = sqrt(k22)
     ind = int(wn + 0.5d0)
     if (ind == 1) wn = 1.0d0
      if (ind < lh - 1) then
!...for velocity field
       if (E(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(spec3d(wn)/E(ind))
       end if    
       vx(jx,jy,jz) = vx(jx,jy,jz)*ratio
       vy(jx,jy,jz) = vy(jx,jy,jz)*ratio
       vz(jx,jy,jz) = vz(jx,jy,jz)*ratio

!...for scalar field
       if (Et(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(spec3d(wn)/Et(ind))
       end if    
       theta(jx,jy,jz) = theta(jx,jy,jz)*ratio


!....compute contribution to the new energy spectrum
        Rex =  real(vx(jx,jy,jz))
        Imx = aimag(vx(jx,jy,jz))
        Rey =  real(vy(jx,jy,jz))
        Imy = aimag(vy(jx,jy,jz))
        Rez =  real(vz(jx,jy,jz))
        Imz = aimag(vz(jx,jy,jz))
        ss = Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy + Rez*Rez + Imz*Imz
        if (kx < 1.d-20) ss = ss * 0.5d0
        if (ss < tiny_real) ss = 0.d0
        Eloc(ind) = Eloc(ind) +  ss

!...compute contribution to the new scalar spectrum
        Rez =  real(theta(jx,jy,jz))
        Imz = aimag(theta(jx,jy,jz))
        ss = Rez*Rez + Imz*Imz
        if (kx < 1.d-20) ss = ss * 0.5d0
        if (ss < tiny_real) ss = 0.d0
        Etloc(ind) = Etloc(ind) +  ss
     else
       vx(jx,jy,jz) = (0.d0,0.d0)
       vy(jx,jy,jz) = (0.d0,0.d0)
       vz(jx,jy,jz) = (0.d0,0.d0)
       theta(jx,jy,jz) = (0.d0,0.d0)
     end if
   else
       vx(jx,jy,jz) = (0.d0,0.d0)
       vy(jx,jy,jz) = (0.d0,0.d0)
       vz(jx,jy,jz) = (0.d0,0.d0)
       theta(jx,jy,jz) = (0.d0,0.d0)
   end if
end do
end do
end do



!...second rescale to match the kinetic energy
write(*,"('Final energy spectrum rescaling (2/2)')")
s2 = 0.5d0 * Eloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
end do
upl = sqrt(s2*2.d0/3.d0)
resc_factor = urms0/upl
print *,'resc= ',resc_factor

vx = vx * resc_factor
vy = vy * resc_factor
vz = vz * resc_factor



write(*,"('final scalar spectrum rescaling (2/2)')")
s2 = 0.5d0 * Etloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Etloc(jx-1) + Etloc(jx))
end do
thetapl = sqrt(s2*2.d0)
resc_factor = thetarms0/thetapl
print *,'scalar resc= ',resc_factor

theta = theta * resc_factor


!................  final spectrum  ....................

close(35)


output = .true.
call spectrum(7,'spec_final.dat')
call scalar_spectrum(1,'scalar_specfinal.dat')
output = .false.


open(unit=21, file='MTLM_vxvyvz0256.bin', form='unformatted')
write(21) vx
write(21) vy
write(21) vz
close(unit=21)

open(unit=21, file='MTLM_theta0256.bin', form='unformatted')
write(21) theta
close(unit=21)


WRITE(*,"(//,' MTLM0256 HAS FINISHED NORMALLY')")


300 format(6X,'+',61('-'),'+',/,6X,'|',61X,'|',/,6X,                          &
 '|  Using the following main parameters in physical space:',5X,'|',/,        &
  6X,'|',61X,'|',/,6X,                                                        &
 '|     turbulent dissipation        : ',ES11.5,14X,'|',/,6X,'|',61X,'|',/,6X,&
 '|     kinematic viscosity          : ',ES11.5,14X,'|',/,6X,'|',61X,'|',/,6X,&
 '|     mesh',25X,': ',I4.4,'x',I3.3,'x',I3.3,13X,'|',/,6X,'|',61X,'|',/,6X,  &
 '|',61X,'|',/,6X,'+',61('-'),'+')

END PROGRAM









!======================================================================
subroutine spectrum(iu,name)
! calculates the energy spectrum and related Fourier-space quantities
!======================================================================
use global
implicit none
integer iu,kw,jx,jy,jz,ind,ky,kz,kyky,kzkz,n,nblocks,jp,jpy,jpz,jxs,jxe,bl
integer jpx
real(8) dm,emax,s2,s3,spec3d,kx,k22,upl
real(8) Rex,Imx,Rey,Imy,Rez,Imz,ss
character(*) name


write(*,"(/,'calculating energy spectrum...',/)")

E = 0.0d0

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1,nth
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   ind = int(sqrt(k22) + 0.5d0)
   if (ind < 1) ind = 1
   if (ind < nth) then
     Rex =  real(vx(jx,jy,jz))
     Imx = aimag(vx(jx,jy,jz))
     Rey =  real(vy(jx,jy,jz))
     Imy = aimag(vy(jx,jy,jz))
     Rez =  real(vz(jx,jy,jz))
     Imz = aimag(vz(jx,jy,jz))
     ss = (Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy + Rez*Rez + Imz*Imz)         
     if (kx < 1.d-20) ss = ss * 0.5d0
     if (ss < tiny_real) ss = 0.d0
     E(ind) = E(ind) +  ss
   end if   
end do
end do
end do


s2 = 0.5d0 * E(1)
do jx=2,lh-2
  s2 = s2 + 0.5d0*(E(jx-1) + E(jx))
end do
upl = sqrt(s2*2.d0/3.d0)
resc_factor = urms/upl

if (iu > 0) then
  open(unit = iu, file = name, status = 'unknown')
  write(iu,"('#      k      E(k) discrete   E(k) analytical')")
  do jx = 1,nth - 1
    write(iu,"(f11.5,2(es16.6))") real(jx), E(jx), spec3d(dble(jx))    

  end do

  s3 = 0.5d0 * E(1)
  do jx=2,lh-2
    s3 = s3 + 0.5d0*(E(jx-1)*((jx-1)*(jx-1)) + E(jx)*(jx*jx))
  end do

  epsl = 2.d0*rnu*s3
  
  write(iu,"(//,'# (Fourier) dissipation    =',ES13.6)") epsl
  write(iu,"('# (Fourier) kinetic energy =',ES13.6)") s2
  write(iu,"('# (Fourier) rms velocity   =',ES13.6)") upl

  write(iu,"('# Length scale/Domain size =',ES13.6)") upl*upl*upl/epsl/Ldom
  write(iu,"('# Eff. Reynolds_{Taylor}   =',ES13.6)") &
                                         sqrt(15.d0/(rnu*epsl))*upl*upl
  write(iu,"('# (Fourier) Viscosity      =',ES14.7)") rnu
  write(iu,"('# Effective phi            =',ES14.7)") &
                                         (Ldom/nx)/((rnu*rnu*rnu/epsl)**0.25d0)
  if(output) then
  write(*,"(//,'# (Fourier) dissipation    =',ES13.6)") epsl
  write(*,"('# (Fourier) kinetic energy =',ES13.6)") s2
  write(*,"('# (Fourier) rms velocity   =',ES13.6)") upl
  write(*,"('# Length scale/Domain size =',ES13.6)") upl*upl*upl/epsl/Ldom
  write(*,"('# Eff. Reynolds_{Taylor}   =',ES13.6)") &
                                        sqrt(15.d0/(rnu*epsl))*upl*upl
  write(*,"('# (Fourier) Viscosity      =',ES14.7)") rnu
  write(*,"('# Effective phi            =',ES14.7)") &
                                        (Ldom/nx)/((rnu*rnu*rnu/epsl)**0.25d0)
  output = .false.
  end if
close(iu)
end if
END subroutine spectrum









!======================================================================
subroutine scalar_spectrum(iu,name)
! calculates the scalar spectrum and related Fourier-space quantities
!======================================================================
use global
implicit none
integer iu,kw,jx,jy,jz,ind,ky,kz,kyky,kzkz,n,nblocks,jp,jpy,jpz,jxs,jxe,bl
integer jpx
real(8) dm,emax,s2,s3,scalar_spec3d,kx,k22,theta_rmsl,scaldissl
real(8) Re,Im,ss
character(*) name


write(*,"(/,'calculating scalar spectrum...',/)")

Et = 0.0d0

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1,nth
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   ind = int(sqrt(k22) + 0.5d0)
   if (ind < 1) ind = 1
   if (ind < nth) then
     Re =  real(theta(jx,jy,jz))
     Im = aimag(theta(jx,jy,jz))
     ss = Re*Re + Im*Im
     if (kx < 1.d-20) ss = ss * 0.5d0
     if (ss < tiny_real) ss = 0.d0
     Et(ind) = Et(ind) +  ss
   end if   
end do
end do
end do


s2 = 0.5d0 * Et(1)
do jx=2,lh-2
  s2 = s2 + 0.5d0*(Et(jx-1) + Et(jx))
end do
theta_rmsl = sqrt(2.d0*s2)
resc_factor = theta_rms/theta_rmsl

if (iu > 0) then
  open(unit = iu, file = name, status = 'unknown')
  write(iu,"('#      k      E(k) discrete   E(k) analytical')")
  do jx = 1,nth - 1
    write(iu,"(f11.5,2(es16.6))") real(jx), Et(jx), scalar_spec3d(dble(jx))    
  end do

  s3 = 0.5d0 * E(1)
  do jx=2,lh-2
    s3 = s3 + 0.5d0*(Et(jx-1)*((jx-1)*(jx-1)) + Et(jx)*(jx*jx))
  end do

  scaldissl = 2.d0 * gamma_diff * s3
  
  write(iu,"(//,'# (Fourier) scalar dissipation    =',ES13.6)") scaldissl
  write(iu,"('# (Fourier) rms scalar =',ES13.6)") theta_rmsl
  write(iu,"('# Effective phi            =',ES14.7)") &
                                      (Ldom/nx)/((gamma_diff**3/diss)**0.25d0)
  if(output) then
  write(*,"(//,'# (Fourier) scalar dissipation    =',ES13.6)") scaldissl
  write(*,"('# (Fourier) rms scalar =',ES13.6)") theta_rmsl
  write(*,"('# Effective phi            =',ES14.7)") &
                                      (Ldom/nx)/((gamma_diff**3/diss)**0.25d0)
  output = .false.
  end if
close(iu)
end if
END subroutine scalar_spectrum












real(8) function spec3d(rk)
!****** 3-D spectrum function ********
use global
implicit none
real(8) rk,kl
kl = rk * l
spec3d = Ck * diss**(2.d0/3.d0) * rk**(-5.d0/3.d0)       &
            * (kl/sqrt(kl*kl + cl))**(5.d0/3.d0 + p0)    &
            * dexp(-1.5d0*Ck*(rk*eta)**4.d0/3.d0)
end function spec3d







real(8) function scalar_spec3d(rk)
!****** 3-D spectrum function for passive scalar********
use global
implicit none
real(8) rk,kl
kl = rk * lq
scalar_spec3d = Cq * qdiss * diss**(-1.d0/3.d0) * rk**(-5.d0/3.d0)   &
                   * (kl/sqrt(kl*kl + clq))**(5.d0/3.d0 + p0)        &
                   * dexp(-1.5d0*Cq*(rk*etaq)**4.d0/3.d0)
end function scalar_spec3d






!============================================================================





Subroutine apply_MTLM(vxred,vyred,vzred,thetared,ldr,IFTvxred,IFTvyred,   &
                      IFTvzred,IFTthred,FTvxred,FTvyred,FTvzred,FTthred)
!...Apply the MMLM with:
!        reduced arrays in RAM memory.
!        main field is disk
!        3-D FFT for reduced arrays.        

use global
implicit none
integer ldr
complex(8), dimension(ldr/2,ldr-2,ldr-2) :: vxred,vyred,vzred,thetared
real(8), dimension(lh-2) :: Eloc,Etloc
complex(8) t,vxx,vyy,vzz,vxa,vya,vza
integer(8) IFTvxred,IFTvyred,IFTvzred,FTvxred,FTvyred,FTvzred,IFTthred,FTthred
real(8) hxr,hyr,hzr,cns,wn,spec3d,ratio,mag,d,k2r,kx,k22,ur,tn,time_scale
real(8) u1,u2,v1,v2,w1,w2,total_adv,tadv,next_total,Int0kmaxr,urb,s2,thetarb
real(8) Rex,Imx,Rey,Imy,Rez,Imz,ss,dm,upl,thetapl,scalar_spec3d
integer lhr,nxr,nyr,nzr,nxrh,nyrh,nzrh,dny,dnz,jxr,jyr,jzr,ind,jx,jy,jz,r
integer ky,kz,kyky,kzkz,Sup_shell,lhrm,n,jp,bl,nblocks,jpx,jpy,jpz
integer ipos,i,j,k,ad,nadv,i1,i2
logical ext,last_adv
character(2) par2
character(4) suf
character conf

lhr = ldr/2
nxr = ldr - 2
nyr = (ldr - 2)
nzr = (ldr - 2)
nxrh = nxr/2
nyrh = nyr/2
nzrh = nzr/2
dny = ny - nyr
dnz = nz - nzr
Sup_shell = (ldr - 2)/2
lhrm = lhr-1

vxred = (0.d0,0.d0)
vyred = (0.d0,0.d0)
vzred = (0.d0,0.d0)
thetared = (0.d0,0.d0)


!......apply cut-off spectral filter 

lhrm = lhr-1

write(*,"('filter 1')")

do jzr = 1,nzrh
   kz = jzr - 1
   kzkz = kz * kz
loop_jyr: do jyr = 1,nyrh
   ky = jyr - 1
   kyky = ky * ky
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxred(jxr,jyr,jzr) = vx(jxr,jyr,jzr)
         vyred(jxr,jyr,jzr) = vy(jxr,jyr,jzr)
         vzred(jxr,jyr,jzr) = vz(jxr,jyr,jzr)
         thetared(jxr,jyr,jzr) = theta(jxr,jyr,jzr)
      end if
   end do
end do loop_jyr
end do



write(*,"('filter 2')")

do jzr = 1,nzrh
   kz = jzr - 1
   kzkz = kz * kz
loop_jyr2: do jyr = nyrh+1,nyr
   jy = jyr + dny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxred(jxr,jyr,jzr) = vx(jxr,jy,jzr)
         vyred(jxr,jyr,jzr) = vy(jxr,jy,jzr)
         vzred(jxr,jyr,jzr) = vz(jxr,jy,jzr)
         thetared(jxr,jyr,jzr) = theta(jxr,jy,jzr)
      end if
   end do
end do loop_jyr2
end do


write(*,"('filter 3')")


do jzr = nzrh+1,nzr
   jz = jzr + dnz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz   
loop_jyr3: do jyr = 1,nyrh
   ky = jyr - 1
   kyky = ky * ky
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxred(jxr,jyr,jzr) = vx(jxr,jyr,jz)
         vyred(jxr,jyr,jzr) = vy(jxr,jyr,jz)
         vzred(jxr,jyr,jzr) = vz(jxr,jyr,jz)
         thetared(jxr,jyr,jzr) = theta(jxr,jyr,jz)
      end if
  end do
end do loop_jyr3
end do



write(*,"('filter 4')") 

do jzr = nzrh+1,nzr
   jz = jzr + dnz  
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
loop_jyr4: do jyr = nyrh+1,nyr
   jy = jyr + dny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky 
   do jxr = 1,lhrm
     kx = dble(jxr - 1)
     k22 = kx*kx + kyky + kzkz
     if (k22 < 1.d-20) cycle
     ind = int(sqrt(k22) + 0.5d0)
     if (ind < Sup_shell) then
        vxred(jxr,jyr,jzr) = vx(jxr,jy,jz)
        vyred(jxr,jyr,jzr) = vy(jxr,jy,jz)
        vzred(jxr,jyr,jzr) = vz(jxr,jy,jz)
        thetared(jxr,jyr,jzr) = theta(jxr,jy,jz)
     end if
  end do
end do loop_jyr4
end do


!...compute rms velocity in the band, for rescaling
s2 = 0.5d0 * spec3d(1.0d0)
do jx=2,lhr-1
  s2 = s2 + 0.5d0*(spec3d(dble(jx-1)) + spec3d(dble(jx)))
end do
Int0kmaxr = s2
urb = sqrt(Int0kmaxr/Int0kmax) * urms0


!...compute rms scalar in the band, for rescaling
s2 = 0.5d0 * scalar_spec3d(1.0d0)
do jx=2,lhr-1
  s2 = s2 + 0.5d0*(scalar_spec3d(dble(jx-1)) + scalar_spec3d(dble(jx)))
end do
Int0kmaxr = s2
thetarb = sqrt(Int0kmaxr/Int0kmaxscalar) * thetarms0

print *,'thetarb ',thetarb


!...first IFT
write(*,"('52/',$)") 
call dfftw_execute(IFTvxred)
call dfftw_execute(IFTvyred)
call dfftw_execute(IFTvzred)
call dfftw_execute(IFTthred)
write(*,"('53/',$)") 


hxr = two_Pi/nxr

ur = 0.0d0
do k = 1,nzr
do j = 1,nyr
do i = 1,nxr/2
   i2 = 2*i
   i1 = i2-1
   u1 =  abs(real(vxred(i,j,k)))
   u2 = abs(aimag(vxred(i,j,k)))
   
   v1 =  abs(real(vyred(i,j,k)))
   v2 = abs(aimag(vyred(i,j,k)))

   w1 =  abs(real(vzred(i,j,k)))
   w2 = abs(aimag(vzred(i,j,k)))
   
   ur = ur + u1*u1 + v1*v1 + w1*w1 + u2*u2 + v2*v2 + w2*w2
end do
end do
end do
ur = sqrt(ur/(3.0d0*nxr*nyr*nzr))
tn = hxr/ur
time_scale = (hxr*hxr/epsilon)**(1.d0/3.d0)

write(*,"(/,'    urms,n = ',ES12.5)") ur
write(*,"('        tn = ',ES12.5)") tn
write(*,"('time_scale = ',ES12.5)") time_scale
write(911,"(//,60('='),/,'  LEVEL = ',I5)") nxr
write(911,"(/,'    urms,n = ',ES12.5)") ur
write(911,"('        tn = ',ES12.5)") tn
write(911,"('time_scale = ',ES12.5)") time_scale


total_adv = 0.d0
last_adv = .false.

nadv = nint(time_scale/tn)
tadv = time_scale/nadv
CFL = tadv/tn
if (CFL > 1.05) then
  nadv = nadv + 1
  tadv = time_scale/nadv
  CFL = tadv/tn
end if  
write(*,"('       CFL = ',ES12.5)") CFL
write(911,"('       CFL = ',ES12.5)") CFL  

DO ad = 1, nadv
   next_total = total_adv + tadv

   call MTLM(vxred,vyred,vzred,thetared,ldr,tadv,urb,thetarb)

   cns = dble(nxr*nyr*nzr)

   vxred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)
   vyred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)
   vzred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)

   write(*,"('54/',$)") 
   !... back to Fourier space
   call dfftw_execute(FTvxred)
   call dfftw_execute(FTvyred)
   call dfftw_execute(FTvzred)
   write(*,"('55/',$)") 

   vxred = vxred/cns
   vyred = vyred/cns
   vzred = vzred/cns


!....velocity field
   vxred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)      !Oddballs
   vyred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)
   vzred(lhr,1:nyr,1:nzr) = (0.d0,0.d0)
   vxred(1:lhr,nyr/2+1,1:nzr) = (0.d0,0.d0)
   vyred(1:lhr,nyr/2+1,1:nzr) = (0.d0,0.d0)
   vzred(1:lhr,nyr/2+1,1:nzr) = (0.d0,0.d0)
   vxred(1:lhr,1:nyr,nzr/2+1) = (0.d0,0.d0)
   vyred(1:lhr,1:nyr,nzr/2+1) = (0.d0,0.d0)
   vzred(1:lhr,1:nyr,nzr/2+1) = (0.d0,0.d0)

   vxred(1,1,1) = (0.d0,0.d0)   !<u_i>=0
   vyred(1,1,1) = (0.d0,0.d0)
   vzred(1,1,1) = (0.d0,0.d0)


!....perform projection in Fourier space
   do jzr = 1,nzrh
      kz = jzr - 1
      kzkz = kz * kz
   do jyr = 1,nyrh
      ky = jyr - 1
      kyky = ky * ky
   do jxr = 1,lhr-1
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxa = vxred(jxr,jyr,jzr)
         vya = vyred(jxr,jyr,jzr)
         vza = vzred(jxr,jyr,jzr)
         t = (vxa*kx + vya*ky + vza*kz)/k22
         vxx = vxa - t*kx
         vyy = vya - t*ky
         vzz = vza - t*kz
         mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))
         d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) + vza * conjg(vza))/mag
         vxred(jxr,jyr,jzr) = vxx*d    
         vyred(jxr,jyr,jzr) = vyy*d
         vzred(jxr,jyr,jzr) = vzz*d
      end if
   end do
   end do
   end do


   do jzr = 1,nzrh
      kz = jzr - 1
      kzkz = kz * kz
   do jyr = nyrh+1,nyr
      jy = jyr + dny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky
   do jxr = 1,lhr-1
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxa = vxred(jxr,jyr,jzr)
         vya = vyred(jxr,jyr,jzr)
         vza = vzred(jxr,jyr,jzr)
         t = (vxa*kx + vya*ky + vza*kz)/k22
         vxx = vxa - t*kx
         vyy = vya - t*ky
         vzz = vza - t*kz
         mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))
         d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) + vza * conjg(vza))/mag
         vxred(jxr,jyr,jzr) = vxx*d    
         vyred(jxr,jyr,jzr) = vyy*d
         vzred(jxr,jyr,jzr) = vzz*d
      end if   
   end do
   end do
   end do


   do jzr = nzrh+1,nzr
      jz = jzr + dnz
      if (jz == nzhp) cycle
      kz = jz - 1
      if (jz > nzh) kz = kz - nz
      kzkz = kz * kz  
   do jyr = 1,nyrh
      ky = jyr - 1
      kyky = ky * ky
   do jxr = 1,lhr-1
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxa = vxred(jxr,jyr,jzr)
         vya = vyred(jxr,jyr,jzr)
         vza = vzred(jxr,jyr,jzr)
         t = (vxa*kx + vya*ky + vza*kz)/k22
         vxx = vxa - t*kx
         vyy = vya - t*ky
         vzz = vza - t*kz
         mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))
         d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) + vza * conjg(vza))/mag
         vxred(jxr,jyr,jzr) = vxx*d    
         vyred(jxr,jyr,jzr) = vyy*d
         vzred(jxr,jyr,jzr) = vzz*d
      end if
   end do
   end do
   end do


   do jzr = nzrh+1,nzr
      jz = jzr + dnz   
      if (jz == nzhp) cycle
      kz = jz - 1
      if (jz > nzh) kz = kz - nz
      kzkz = kz * kz
   do jyr = nyrh+1,nyr
      jy = jyr + dny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky
   do jxr = 1,lhr-1
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vxa = vxred(jxr,jyr,jzr)
         vya = vyred(jxr,jyr,jzr)
         vza = vzred(jxr,jyr,jzr)
         t = (vxa*kx + vya*ky + vza*kz)/k22
         vxx = vxa - t*kx
         vyy = vya - t*ky
         vzz = vza - t*kz
         mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) + vzz*conjg(vzz))
         d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) + vza * conjg(vza))/mag
         vxred(jxr,jyr,jzr) = vxx*d    
         vyred(jxr,jyr,jzr) = vyy*d
         vzred(jxr,jyr,jzr) = vzz*d
      end if
   end do
   end do
   end do

   vxred(1,1,1) = 0.d0   !Mean velocity=0
   vyred(1,1,1) = 0.d0
   vzred(1,1,1) = 0.d0
   write(*,"('57/',$)") 
 
   total_adv = total_adv + tadv
   write(*,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv, total_adv   
   write(911,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv,total_adv
   write(911,"(/,30('+-'),/)")
   if (ad /= nadv) then
!...IFT for the next advancement
      write(*,"('sIFT/',$)") 
      call dfftw_execute(IFTvxred)
      call dfftw_execute(IFTvyred)
      call dfftw_execute(IFTvzred)
      write(*,"('eIFT',$)") 
   end if
END DO


!...transform the scalar field back to Fourier space
write(*,"(/,'scalar FT')") 
call dfftw_execute(FTthred)

thetared = thetared/cns

thetared(lhr,1:nyr,1:nzr) = (0.d0,0.d0)
thetared(1:lhr,nyr/2+1,1:nzr) = (0.d0,0.d0)
thetared(1:lhr,1:nyr,nzr/2+1) = (0.d0,0.d0)

thetared(1,1,1) = (0.d0,0.d0)


!...put back the modified shell into the whole field

write(*,"('write-back 1')")
do jzr = 1,nzrh
   kz = jzr - 1
   kzkz = kz * kz
loop_jyr5: do jyr = 1,nyrh
   ky = jyr - 1
   kyky = ky * ky
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
         vx(jxr,jyr,jzr) = vxred(jxr,jyr,jzr)
         vy(jxr,jyr,jzr) = vyred(jxr,jyr,jzr)
         vz(jxr,jyr,jzr) = vzred(jxr,jyr,jzr)         
         theta(jxr,jyr,jzr) = thetared(jxr,jyr,jzr)
      end if
   end do  
end do loop_jyr5
end do


!----------------------------------------------------------------

write(*,"('write-back 2')")

do jzr = 1,nzrh
   kz = jzr - 1
   kzkz = kz * kz
loop_jyr6: do jyr = nyrh+1,nyr
   jy = jyr + dny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky 
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
        vx(jxr,jy,jzr) = vxred(jxr,jyr,jzr)
        vy(jxr,jy,jzr) = vyred(jxr,jyr,jzr)
        vz(jxr,jy,jzr) = vzred(jxr,jyr,jzr)
        theta(jxr,jy,jzr) = thetared(jxr,jyr,jzr)
      end if
   end do
end do loop_jyr6
end do     


!----------------------------------------------------------------

write(*,"('write-back 3')")

do jzr = nzrh+1,nzr
   jz = jzr + dnz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz  
loop_jyr7: do jyr = 1,nyrh
   ky = jyr - 1
   kyky = ky * ky         
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
        vx(jxr,jyr,jz) = vxred(jxr,jyr,jzr)
        vy(jxr,jyr,jz) = vyred(jxr,jyr,jzr)
        vz(jxr,jyr,jz) = vzred(jxr,jyr,jzr)         
        theta(jxr,jyr,jz) = thetared(jxr,jyr,jzr)         
      end if
   end do
end do loop_jyr7
end do     

!----------------------------------------------------------------

write(*,"('write-back 4')")

do jzr = nzrh+1,nzr
   jz = jzr + dnz   
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
loop_jyr8: do jyr = nyrh+1,nyr
   jy = jyr + dny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
   do jxr = 1,lhrm
      kx = dble(jxr - 1)
      k22 = kx*kx + kyky + kzkz
      if (k22 < 1.d-20) cycle
      ind = int(sqrt(k22) + 0.5d0)
      if (ind < Sup_shell) then
        vx(jxr,jy,jz) = vxred(jxr,jyr,jzr)
        vy(jxr,jy,jz) = vyred(jxr,jyr,jzr)
        vz(jxr,jy,jz) = vzred(jxr,jyr,jzr)  
        theta(jxr,jy,jz) = thetared(jxr,jyr,jzr)         
      end if
   end do
end do loop_jyr8
end do



!...calculate new spectra for rescaling (velocity and scalar fields)

write(suf,"(I4.4)") nxr
call spectrum(19,'afterspec'//suf//'.dat') 

call scalar_spectrum(21,'safterspec'//suf//'.dat') 


Eloc = 0.0d0
Etloc = 0.0d0

!...re-scale the spectra
write(*,"('Rescaling the spectra (1/2)')")

do jz = 1,nz
   if (jz == nzhp) cycle
   kz = jz - 1
   if (jz > nzh) kz = kz - nz
   kzkz = kz * kz
do jy = 1,ny
   if (jy == nyhp) cycle
   ky = jy - 1
   if (jy > nyh) ky = ky - ny
   kyky = ky * ky
do jx = 1,lhm
   kx = dble(jx - 1)
   k22 = kx*kx + kyky + kzkz
   if (k22 < 1.d-20) cycle
   if (k22 < max_k2) then
     wn = sqrt(k22)
     ind = int(wn + 0.5d0)
     if (ind <= 1) then
       ind = 1
       wn = 1.0d0
     end if  
     if (ind < nth) then
!...for velocity field
       if (E(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(spec3d(wn)/E(ind))
       end if    
       vx(jx,jy,jz) = vx(jx,jy,jz)*ratio
       vy(jx,jy,jz) = vy(jx,jy,jz)*ratio
       vz(jx,jy,jz) = vz(jx,jy,jz)*ratio

!...for scalar field
       if (Et(ind) < 2.d-38) then
          ratio = 0.0d0
       else  
          ratio = sqrt(scalar_spec3d(wn)/Et(ind))
       end if    
       theta(jx,jy,jz) = theta(jx,jy,jz)*ratio

!...compute contribution to the new energy spectrum
        Rex =  real(vx(jx,jy,jz))
        Imx = aimag(vx(jx,jy,jz))
        Rey =  real(vy(jx,jy,jz))
        Imy = aimag(vy(jx,jy,jz))
        Rez =  real(vz(jx,jy,jz))
        Imz = aimag(vz(jx,jy,jz))
        ss = Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy + Rez*Rez + Imz*Imz
        if (kx < 1.d-20) ss = ss * 0.5d0
        if (ss < tiny_real) ss = 0.d0
        Eloc(ind) = Eloc(ind) +  ss

!...compute contribution to the new scalar spectrum
        Rez =  real(theta(jx,jy,jz))
        Imz = aimag(theta(jx,jy,jz))
        ss = Rez*Rez + Imz*Imz
        if (kx < 1.d-20) ss = ss * 0.5d0
        if (ss < tiny_real) ss = 0.d0
        Etloc(ind) = Etloc(ind) +  ss
     else
       vx(jx,jy,jz) = (0.d0,0.d0)
       vy(jx,jy,jz) = (0.d0,0.d0)
       vz(jx,jy,jz) = (0.d0,0.d0)
       theta(jx,jy,jz) = (0.d0,0.d0)
     end if
   else
       vx(jx,jy,jz) = (0.d0,0.d0)
       vy(jx,jy,jz) = (0.d0,0.d0)
       vz(jx,jy,jz) = (0.d0,0.d0)
       theta(jx,jy,jz) = (0.d0,0.d0)
   end if
end do
end do
end do


!...second rescale to match the kinetic energy
write(*,"('Rescaling the energy spectrum (2/2)')")

s2 = 0.5d0 * Eloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
end do
upl = sqrt(s2*2.d0/3.d0)
resc_factor = urms0/upl
print *,'resc= ',resc_factor

vx = vx * resc_factor
vy = vy * resc_factor
vz = vz * resc_factor


write(*,"('Rescaling the scalar spectrum (2/2)')")

s2 = 0.5d0 * Etloc(1)
do jx = 2,lh-2
  s2 = s2 + 0.5d0*(Etloc(jx-1) + Etloc(jx))
end do
thetapl = sqrt(s2*2.d0)
resc_factor = thetarms0/thetapl
print *,'scalar resc= ',resc_factor

theta = theta * resc_factor

End subroutine apply_MTLM





!============================================================================


SUBROUTINE MTLM(vxred,vyred,vzred,thetared,ldr,tadv,urmsb,thetarmsb)
!...field distortion using 'MTLM' displacement. Velocities are recovered
!...by gridding on the irregular distorted mesh.
!...Full RAM mapping.
use global
implicit none

integer ldr
complex(8), dimension(ldr/2,ldr-2,ldr-2) :: vxred,vyred,vzred,thetared
real(8), dimension(ldr-2,ldr-2,ldr-2) :: ugrid,vgrid,wgrid,Swf,thetagrid

integer, dimension(ldr-2) :: is,ie
integer, dimension(ldr-2) :: js,je
integer, dimension(ldr-2) :: ks,ke

real(8) dt,maxd,wf,hxr,hyr,hzr,vel1,vel2,f1,f2,f3,R,minwf,summat,tadv,urmsb,ur
real(8) ua1,va1,wa1,ua2,va2,wa2,ir,jr,kr,dirl2,dirh2,djrl2,djrh2,dkrl2,dkrh2
integer il,ih,jl,jh,kl,kh,ia,ja,ka,ia1,ia2,i,j,k,nxr,nyr,nzr,i1,i2
real(8) swfs,Swfu,Swfv,Swfw,R2,Ra,u1,u2,v1,v2,w1,w2,rfac,thetarmsb
real(8) thetaa1,thetaa2,rfacth,Swftheta,thetamean,th,thetar
real(8) dx,dy,dz,dist2,ip,jp,kp,nxr8,nyr8,nzr8,norm,S11,S12,S13,S22,S23,S33
integer in,jn,kn,isearch,jsearch,ksearch,inh,ini
integer voids



write(*,"('90/',$)") 

nxr = ldr - 2
nyr = nxr
nzr = nxr
hxr = two_Pi/nxr
hyr = two_Pi/nyr
hzr = two_Pi/nzr

nxr8 = dble(nxr) 
nyr8 = dble(nyr) 
nzr8 = dble(nzr) 

dt = tadv

!...define sphere for neighbors searching
R = 1.01d0
minwf = 1.d0/R

!...apply 'MMLM' displacements
dt = dt/hxr
ugrid = 0.d0
vgrid = 0.d0
wgrid = 0.d0
thetagrid = 0.d0
Swf = 0.d0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


do ka = 1,nzr
if((ka/4)*4 == ka) write(*,"('*',$)")
do ja = 1,nyr
do ia = 1,nxr/2
   ia2 = 2*ia
   ia1 = ia2 - 1

   ua1 = real(vxred(ia,ja,ka))
   va1 = real(vyred(ia,ja,ka))
   wa1 = real(vzred(ia,ja,ka))      
   thetaa1 = real(thetared(ia,ja,ka))

   ua2 = aimag(vxred(ia,ja,ka))
   va2 = aimag(vyred(ia,ja,ka))
   wa2 = aimag(vzred(ia,ja,ka))
   thetaa2 = aimag(thetared(ia,ja,ka))

!****(1) 
   ir = dble(ia1) + dt*ua1
   jr = dble(ja) + dt*va1
   kr = dble(ka) + dt*wa1
   
   il = floor(ir)
   ih = il + 1
   jl = floor(jr)
   jh = jl + 1
   kl = floor(kr)
   kh = kl + 1

   dirl2 = ir - dble(il)
   dirh2 = ir - dble(ih)   
   djrl2 = jr - dble(jl)                   
   djrh2 = jr - dble(jh)
   dkrl2 = kr - dble(kl)                  
   dkrh2 = kr - dble(kh)

   dirl2 = dirl2 * dirl2
   dirh2 = dirh2 * dirh2
   djrl2 = djrl2 * djrl2
   djrh2 = djrh2 * djrh2
   dkrl2 = dkrl2 * dkrl2
   dkrh2 = dkrh2 * dkrh2

   if (il > nxr) then
     il = il - nxr
   else if (il < 1) then
     il = il + nxr
   end if  

   if (ih > nxr) then
     ih = ih - nxr
   else if (ih < 1) then
     ih = ih + nxr
   end if  
     
   if (jl > nyr) then
     jl = jl - nyr
   else if (jl < 1) then
     jl = jl + nyr
   end if  

   if (jh > nyr) then
     jh = jh - nyr
   else if (jh < 1) then
     jh = jh + nyr
   end if  

   if (kl > nzr) then
     kl = kl - nzr
   else if (kl < 1) then
     kl = kl + nzr
   end if  

   if (kh > nzr) then
     kh = kh - nzr
   else if (kh < 1) then
     kh = kh + nzr
   end if  
   
!print *,il,ih,jl,jh,kl,kh
       
!...000
   wf = 1.d0/sqrt(dirl2 + djrl2 + dkrl2)
   if (wf > minwf) then
   ugrid(il,jl,kl) = ugrid(il,jl,kl) + wf*ua1
   vgrid(il,jl,kl) = vgrid(il,jl,kl) + wf*va1
   wgrid(il,jl,kl) = wgrid(il,jl,kl) + wf*wa1  
   thetagrid(il,jl,kl) = thetagrid(il,jl,kl) + wf*thetaa1  
   Swf(il,jl,kl) = Swf(il,jl,kl) + wf
   end if 
!...001
   wf = 1.d0/sqrt(dirl2 + djrl2 + dkrh2)
   if (wf > minwf) then
   ugrid(il,jl,kh) = ugrid(il,jl,kh) + wf*ua1
   vgrid(il,jl,kh) = vgrid(il,jl,kh) + wf*va1
   wgrid(il,jl,kh) = wgrid(il,jl,kh) + wf*wa1   
   thetagrid(il,jl,kh) = thetagrid(il,jl,kh) + wf*thetaa1
   Swf(il,jl,kh) = Swf(il,jl,kh) + wf
   end if
   
!...010
   wf = 1.d0/sqrt(dirl2 + djrh2 + dkrl2)
   if (wf > minwf) then
   ugrid(il,jh,kl) = ugrid(il,jh,kl) + wf*ua1
   vgrid(il,jh,kl) = vgrid(il,jh,kl) + wf*va1
   wgrid(il,jh,kl) = wgrid(il,jh,kl) + wf*wa1   
   thetagrid(il,jh,kl) = thetagrid(il,jh,kl) + wf*thetaa1   
   Swf(il,jh,kl) = Swf(il,jh,kl) + wf
   end if

!...011
   wf = 1.d0/sqrt(dirl2 + djrh2 + dkrh2)
   if (wf > minwf) then
   ugrid(il,jh,kh) = ugrid(il,jh,kh) + wf*ua1
   vgrid(il,jh,kh) = vgrid(il,jh,kh) + wf*va1
   wgrid(il,jh,kh) = wgrid(il,jh,kh) + wf*wa1   
   thetagrid(il,jh,kh) = thetagrid(il,jh,kh) + wf*thetaa1   
   Swf(il,jh,kh) = Swf(il,jh,kh) + wf
   end if

!...100
   wf = 1.d0/sqrt(dirh2 + djrl2 + dkrl2)
   if (wf > minwf) then
   ugrid(ih,jl,kl) = ugrid(ih,jl,kl) + wf*ua1
   vgrid(ih,jl,kl) = vgrid(ih,jl,kl) + wf*va1
   wgrid(ih,jl,kl) = wgrid(ih,jl,kl) + wf*wa1  
   thetagrid(ih,jl,kl) = thetagrid(ih,jl,kl) + wf*thetaa1  
   Swf(ih,jl,kl) = Swf(ih,jl,kl) + wf
   end if

!...101
   wf = 1.d0/sqrt(dirh2 + djrl2 + dkrh2)
   if (wf > minwf) then
   ugrid(ih,jl,kh) = ugrid(ih,jl,kh) + wf*ua1
   vgrid(ih,jl,kh) = vgrid(ih,jl,kh) + wf*va1
   wgrid(ih,jl,kh) = wgrid(ih,jl,kh) + wf*wa1   
   thetagrid(ih,jl,kh) = thetagrid(ih,jl,kh) + wf*thetaa1   
   Swf(ih,jl,kh) = Swf(ih,jl,kh) + wf
   end if

!...110
   wf = 1.d0/sqrt(dirh2 + djrh2 + dkrl2)
   if (wf > minwf) then
   ugrid(ih,jh,kl) = ugrid(ih,jh,kl) + wf*ua1
   vgrid(ih,jh,kl) = vgrid(ih,jh,kl) + wf*va1
   wgrid(ih,jh,kl) = wgrid(ih,jh,kl) + wf*wa1   
   thetagrid(ih,jh,kl) = thetagrid(ih,jh,kl) + wf*thetaa1   
   Swf(ih,jh,kl) = Swf(ih,jh,kl) + wf
   end if

!...111
   wf = 1.d0/sqrt(dirh2 + djrh2 + dkrh2)
   if (wf > minwf) then
   ugrid(ih,jh,kh) = ugrid(ih,jh,kh) + wf*ua1
   vgrid(ih,jh,kh) = vgrid(ih,jh,kh) + wf*va1
   wgrid(ih,jh,kh) = wgrid(ih,jh,kh) + wf*wa1   
   thetagrid(ih,jh,kh) = thetagrid(ih,jh,kh) + wf*thetaa1   
   Swf(ih,jh,kh) = Swf(ih,jh,kh) + wf
   end if


!****(2) 
   ir = dble(ia2) + dt*ua2
   jr = dble(ja) + dt*va2
   kr = dble(ka) + dt*wa2

   il = floor(ir)
   ih = il + 1
   jl = floor(jr)
   jh = jl + 1
   kl = floor(kr)
   kh = kl + 1

   dirl2 = ir - dble(il)
   dirh2 = ir - dble(ih)   
   djrl2 = jr - dble(jl)                   
   djrh2 = jr - dble(jh)
   dkrl2 = kr - dble(kl)                   
   dkrh2 = kr - dble(kh)

   dirl2 = dirl2 * dirl2
   dirh2 = dirh2 * dirh2
   djrl2 = djrl2 * djrl2
   djrh2 = djrh2 * djrh2
   dkrl2 = dkrl2 * dkrl2
   dkrh2 = dkrh2 * dkrh2

!!!!print *,il,ih,jl,jh,kl,kh

   if (il > nxr) then
     il = il - nxr
   else if (il < 1) then
     il = il + nxr
   end if  

   if (ih > nxr) then
     ih = ih - nxr
   else if (ih < 1) then
     ih = ih + nxr
   end if  
     
   if (jl > nyr) then
     jl = jl - nyr
   else if (jl < 1) then
     jl = jl + nyr
   end if  

   if (jh > nyr) then
     jh = jh - nyr
   else if (jh < 1) then
     jh = jh + nyr
   end if  

   if (kl > nzr) then
     kl = kl - nzr
   else if (kl < 1) then
     kl = kl + nzr
   end if  

   if (kh > nzr) then
     kh = kh - nzr
   else if (kh < 1) then
     kh = kh + nzr
   end if  


!...000
   wf = 1.d0/sqrt(dirl2 + djrl2 + dkrl2)
   if (wf > minwf) then
   ugrid(il,jl,kl) = ugrid(il,jl,kl) + wf*ua2
   vgrid(il,jl,kl) = vgrid(il,jl,kl) + wf*va2
   wgrid(il,jl,kl) = wgrid(il,jl,kl) + wf*wa2   
   thetagrid(il,jl,kl) = thetagrid(il,jl,kl) + wf*thetaa2   
   Swf(il,jl,kl) = Swf(il,jl,kl) + wf
   end if

!...001
   wf = 1.d0/sqrt(dirl2 + djrl2 + dkrh2)
   if (wf > minwf) then
   ugrid(il,jl,kh) = ugrid(il,jl,kh) + wf*ua2
   vgrid(il,jl,kh) = vgrid(il,jl,kh) + wf*va2
   wgrid(il,jl,kh) = wgrid(il,jl,kh) + wf*wa2   
   thetagrid(il,jl,kh) = thetagrid(il,jl,kh) + wf*thetaa2   
   Swf(il,jl,kh) = Swf(il,jl,kh) + wf
   end if

!...010
   wf = 1.d0/sqrt(dirl2 + djrh2 + dkrl2)
   if (wf > minwf) then
   ugrid(il,jh,kl) = ugrid(il,jh,kl) + wf*ua2
   vgrid(il,jh,kl) = vgrid(il,jh,kl) + wf*va2
   wgrid(il,jh,kl) = wgrid(il,jh,kl) + wf*wa2   
   thetagrid(il,jh,kl) = thetagrid(il,jh,kl) + wf*thetaa2   
   Swf(il,jh,kl) = Swf(il,jh,kl) + wf
   end if

!...011
   wf = 1.d0/sqrt(dirl2 + djrh2 + dkrh2)
   if (wf > minwf) then
   ugrid(il,jh,kh) = ugrid(il,jh,kh) + wf*ua2
   vgrid(il,jh,kh) = vgrid(il,jh,kh) + wf*va2
   wgrid(il,jh,kh) = wgrid(il,jh,kh) + wf*wa2   
   thetagrid(il,jh,kh) = thetagrid(il,jh,kh) + wf*thetaa2   
   Swf(il,jh,kh) = Swf(il,jh,kh) + wf
   end if

!...100
   wf = 1.d0/sqrt(dirh2 + djrl2 + dkrl2)
   if (wf > minwf) then
   ugrid(ih,jl,kl) = ugrid(ih,jl,kl) + wf*ua2
   vgrid(ih,jl,kl) = vgrid(ih,jl,kl) + wf*va2
   wgrid(ih,jl,kl) = wgrid(ih,jl,kl) + wf*wa2   
   thetagrid(ih,jl,kl) = thetagrid(ih,jl,kl) + wf*thetaa2   
   Swf(ih,jl,kl) = Swf(ih,jl,kl) + wf
   end if

!...101
   wf = 1.d0/sqrt(dirh2 + djrl2 + dkrh2)
   if (wf > minwf) then   
   ugrid(ih,jl,kh) = ugrid(ih,jl,kh) + wf*ua2
   vgrid(ih,jl,kh) = vgrid(ih,jl,kh) + wf*va2
   wgrid(ih,jl,kh) = wgrid(ih,jl,kh) + wf*wa2   
   thetagrid(ih,jl,kh) = thetagrid(ih,jl,kh) + wf*thetaa2   
   Swf(ih,jl,kh) = Swf(ih,jl,kh) + wf
   end if

!...110
   wf = 1.d0/sqrt(dirh2 + djrh2 + dkrl2)
   if (wf > minwf) then   
   ugrid(ih,jh,kl) = ugrid(ih,jh,kl) + wf*ua2
   vgrid(ih,jh,kl) = vgrid(ih,jh,kl) + wf*va2
   wgrid(ih,jh,kl) = wgrid(ih,jh,kl) + wf*wa2   
   thetagrid(ih,jh,kl) = thetagrid(ih,jh,kl) + wf*thetaa2   
   Swf(ih,jh,kl) = Swf(ih,jh,kl) + wf
   end if

!...111
   wf = 1.d0/sqrt(dirh2 + djrh2 + dkrh2)
   if (wf > minwf) then
   ugrid(ih,jh,kh) = ugrid(ih,jh,kh) + wf*ua2
   vgrid(ih,jh,kh) = vgrid(ih,jh,kh) + wf*va2
   wgrid(ih,jh,kh) = wgrid(ih,jh,kh) + wf*wa2   
   thetagrid(ih,jh,kh) = thetagrid(ih,jh,kh) + wf*thetaa2   
   Swf(ih,jh,kh) = Swf(ih,jh,kh) + wf
   end if
end do
end do
end do




voids = 0

do k = 1,nzr
do j = 1,nyr
do i = 1,nxr
  if (Swf(i,j,k) /= 0.d0) then
     ugrid(i,j,k) = ugrid(i,j,k)/Swf(i,j,k)
     vgrid(i,j,k) = vgrid(i,j,k)/Swf(i,j,k)
     wgrid(i,j,k) = wgrid(i,j,k)/Swf(i,j,k)
     thetagrid(i,j,k) = thetagrid(i,j,k)/Swf(i,j,k) 
  else 
     voids = voids + 1
  end if
end do
end do
end do
  

!#########  second search (only for isolated nodes)   #########

maxd = max(CFL,1.d0)     !......because hx = hy = hz

DO

write(*,"(' % voids : ',ES10.3)") real(voids)/(nxr*nyr*nzr) * 100.
write(911,"(' % voids : ',ES10.3)") real(voids)/(nxr*nyr*nzr) * 100.d0
R = R + 1.0d0
R2 = R*R
Ra = R + maxd

if (voids == 0) exit

voids = 0

!...calculate limits fpr searching in the mesh 
if (Ra < 0.5d0*nyr) then
  do i = 1,nxr
    is(i) = ceiling(dble(i) - Ra)
    ie(i) = floor(dble(i) + Ra)
  end do
  js(1:nyr) = is(1:nyr)
  je(1:nyr) = ie(1:nyr)
else
  is = 1
  ie = nxr
  js = 1
  je = nyr
  print *,'  Searching in the whole mesh...'
end if
ks = js
ke = je

write(*,"(/)")     
 
do k = 1,nzr
   if((k/4)*4 == k) write(*,"('#',$)")     
do j = 1,nyr
do i = 1,nxr
   if (Swf(i,j,k) /= 0.d0) cycle

   Swfs = 0.d0
   Swfu = 0.d0
   Swfv = 0.d0
   Swfw = 0.d0
   Swftheta = 0.d0

   do ksearch = ks(k),ke(k)
      kn = ksearch
      if (kn < 1) then 
        kn = kn + nzr
      else if (kn > nzr) then
        kn = kn - nzr
      end if   

   do jsearch = js(j),je(j)
      jn = jsearch
      if (jn < 1) then 
        jn = jn + nyr
      else if (jn > nyr) then
        jn = jn - nyr
      end if   

   do isearch = is(i),ie(i)
      in = isearch
      if (in < 1) then 
        in = in + nxr
      else if (in > nxr) then
        in = in - nxr
      end if   

      inh = in/2
      if (inh*2 == in) then
        ua1 = aimag(vxred(inh,jn,kn))
        va1 = aimag(vyred(inh,jn,kn))
        wa1 = aimag(vzred(inh,jn,kn))
        thetaa1 = aimag(thetared(inh,jn,kn))
      else
        ini = (in+1)/2
        ua1 = real(vxred(ini,jn,kn))
        va1 = real(vyred(ini,jn,kn))
        wa1 = real(vzred(ini,jn,kn))
        thetaa1 = real(thetared(ini,jn,kn))
      end if
       
      ip = dble(in) + ua1 * dt
      jp = dble(jn) + va1 * dt
      kp = dble(kn) + wa1 * dt   

      if (ip < 1.d0) then
      ip = ip + nxr8
      else if (ip > nxr8) then
      ip = ip - nxr8
      end if

      if (jp < 1.d0) then
      jp = jp + nyr8 
      else if (jp > nyr8) then
      jp = jp - nyr8
      end if   

      if (kp < 1.d0) then
      kp = kp + nzr8 
      else if (kp > nzr8) then
      kp = kp - nzr8
      end if

      dx = dabs(ip - dble(i))
      dx = dmin1(dx,nxr8 - dx)
      dy = dabs(jp - dble(j))
      dy = dmin1(dy,nyr8 - dy)
      dz = dabs(kp - dble(k))
      dz = dmin1(dz,nzr8 - dz)
      dist2 = dx*dx + dy*dy + dz*dz  
  
      if (dist2 <= R2) then
        wf = 1.d0/sqrt(dist2)
        Swfu = Swfu + wf * ua1
        Swfv = Swfv + wf * va1
        Swfw = Swfw + wf * wa1
        Swftheta = Swftheta + wf * thetaa1
        Swfs = Swfs + wf
      end if
   end do
   end do
   end do
   if (Swfs /= 0.d0) then 
      ugrid(i,j,k) = Swfu/Swfs
      vgrid(i,j,k) = Swfv/Swfs
      wgrid(i,j,k) = Swfw/Swfs
      thetagrid(i,j,k) = Swftheta/Swfs
   else
      voids = voids + 1
   end if   
end do
end do
end do

END DO


!...rescaling of the Urms
ur = 0.0d0
do k = 1,nzr
do j = 1,nyr
do i = 1,nxr
   u1 =  ugrid(i,j,k)
   v1 =  vgrid(i,j,k)
   w1 =  wgrid(i,j,k)  
   ur = ur + u1*u1 + v1*v1 + w1*w1
end do
end do
end do
ur = sqrt(ur/(3.0d0*nxr*nyr*nzr))

rfac = urmsb/ur
write(911,"('urmsb = ',ES14.6)") urmsb
write(911,"('   ur = ',ES14.6)") ur
write(911,"(' rfac = ',ES14.6)") rfac



!...rescaling of the scalar variance
thetamean = 0.0d0
do k = 1,nzr
do j = 1,nyr
do i = 1,nxr
   thetamean = thetamean + thetagrid(i,j,k)
end do
end do
end do
thetamean = thetamean/(nxr*nyr*nzr)

thetagrid(:,:,:) = thetagrid(:,:,:) - thetamean


thetar = 0.0d0
do k = 1,nzr
do j = 1,nyr
do i = 1,nxr
   th =  thetagrid(i,j,k)
   thetar = thetar + th*th
end do
end do
end do
thetar = sqrt(thetar/(nxr*nyr*nzr))

rfacth = thetarmsb/thetar
write(911,"('thetamean = ',ES14.6,'  before')") thetamean
write(911,"('thetarmsb = ',ES14.6)") thetarmsb
write(911,"('   thetar = ',ES14.6)") thetar
write(911,"('   rfacth = ',ES14.6)") rfacth


do k = 1,nzr
do j = 1,nyr
do i = 1,nxr/2
   i2 = 2*i
   i1 = i2 - 1
   vxred(i,j,k) = cmplx(rfac*ugrid(i1,j,k),rfac*ugrid(i2,j,k),kind=8)
   vyred(i,j,k) = cmplx(rfac*vgrid(i1,j,k),rfac*vgrid(i2,j,k),kind=8)
   vzred(i,j,k) = cmplx(rfac*wgrid(i1,j,k),rfac*wgrid(i2,j,k),kind=8)
   
   thetared(i,j,k) = cmplx(rfacth*thetagrid(i1,j,k),      &
                           rfacth*thetagrid(i2,j,k),kind=8)
end do
end do
end do


!!!!if(k > 0) stop

write(*,"(' ')")
!10 format(' *** FATAL: node ',3I3,' is alone!')
12 format(' *** Using the nearest two neighbors at node ',3I3)
!100 format(3I4)
!110 format(10X,3I4)

END SUBROUTINE MTLM



!============================================================================


real function ran2_true(idum)
implicit none
integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
real AM,EPS,RNMX
parameter(IM1 = 2147483563, IM2 = 2147483399, AM = 1./IM1, IMM1 = IM1 - 1, &
          IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774, IR1 = 12211, &
          IR2 = 3791, NTAB = 32, NDIV = 1+IMM1/NTAB, EPS=1.2e-7, RNMX=1. - EPS) 
integer idum2,j,k,iv(NTAB),iy
save iv,iy,idum2
data idum2/123456789/, iv/NTAB*0/, iy/0/

if (idum <= 0) then
   idum = max(-idum,1)
   idum2 = idum
   do j = NTAB+8,1,-1
      k = idum/IQ1
      idum = IA1*(idum - k*IQ1) - k*IR1
      if (idum < 0) idum = idum + IM1
      if (j <= NTAB) iv(j) = idum
   end do
   iy = iv(1)
end if

k = idum/IQ1
idum = IA1 * (idum - k*IQ1) - k*IR1
if (idum < 0) idum = idum + IM1
k = idum2/IQ2
idum2 = IA2 * (idum2 - k*IQ2) - k*IR2
j = 1 + iy/NDIV
iy = iv(j) - idum2
iv(j) = idum
if (iy < 1) iy = iy + IMM1
ran2_true = min(AM*iy,RNMX)
END function






real function ran2(idum)
!  Random generator from Numerical Recipes, p.271
implicit none
integer idum,ia,im,iq,ir,ntab,ndiv
real am,eps,rnmx
parameter(ia = 16807, im = 2147483647, am = 1./im, iq = 127773, ir = 2836,  &
          ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7, rnmx = 1. - eps)
integer j,k
integer, dimension(ntab), save :: iv=(/(0,j=1,ntab)/)
integer, save :: iy = 0


if (idum <= 0 .or. iy == 0) then
  idum = max(-idum,1)
  do j = ntab + 8, 1, -1 
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k
    if (idum < 0) idum = idum + im
    if (j <= ntab) iv(j) = idum
  end do
  iy = iv(1)
end if

k = idum/iq
idum = ia*(idum - k*iq) - ir*k
if (idum < 0) idum = idum + im
j = 1 + iy/ndiv
iy = iv(j)
iv(j) = idum
ran2 = min(am*iy,rnmx)
End function

