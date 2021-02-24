module minitialize

  interface initialize
      module procedure initialize3dvellowmem
      module procedure initialize3dvelhighmem
      module procedure initialize2dvelhighmem
      module procedure initialize3dvellowmemdealiased
      module procedure initialize3dvelhighmemdealiased
      module procedure initialize3dvelhighmemdealiasedwithphi
  end interface initialize 

contains

  subroutine initialize3dvelhighmemdealiasedwithphi(vx,vy,vz,phi,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    use mconstant
    use msymmetrize
    implicit none
  
    integer,  intent(in) :: lx1,ly,lz,iseed
    real(sp), intent(in) :: akp, u0, kcut
    complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz,phi
    real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2,kx,ky,kz
  
    real(sp), dimension(lx1,ly,lz) :: tmp
  
    integer  :: idum,ii,jj,kk
    real(sp) :: tmp1,tmp2, kcut2, kxii, kyjj, kzkk
    complex(sp) :: tmpz

  
    kcut2 = kcut * kcut

    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1

      kzkk=kz(ii,jj,kk)
      kyjj=ky(ii,jj,kk)
      kxii=kx(ii,jj,kk)

      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      phi(ii,jj,kk) = gasdev(tmp1,tmp2)

      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vx(ii,jj,kk)=gasdev(tmp1,tmp2)
      
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vy(ii,jj,kk)=gasdev(tmp1,tmp2)
  
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vz(ii,jj,kk)=gasdev(tmp1,tmp2)
  
      tmpz=kxii*vx(ii,jj,kk)+kyjj*vy(ii,jj,kk)+kzkk*vz(ii,jj,kk)
      tmpz=tmpz/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
      vx(ii,jj,kk)=vx(ii,jj,kk)-kxii*tmpz
      vy(ii,jj,kk)=vy(ii,jj,kk)-kyjj*tmpz
      vz(ii,jj,kk)=vz(ii,jj,kk)-kzkk*tmpz
    end do
    end do
    end do

    call symmetrize(vx,lx1,ly,lz)
    call symmetrize(vy,lx1,ly,lz)
    call symmetrize(vz,lx1,ly,lz) 
    call symmetrize(phi,lx1,ly,lz) 

    where(k2 .ge. kcut2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
        phi = 0._sp
    endwhere
  
    ! different initial conditions can be imposed here
  
    tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 
  
    vx = vx * tmp
    vy = vy * tmp
    vz = vz * tmp

    phi = phi * tmp
  
    tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
    tmp(1,:,:)=.5_sp*tmp(1,:,:)
    tmp1=sum(tmp)
    tmp2=(3._sp/2._sp)*u0*u0
    vx=sqrt(tmp2/tmp1)*vx
    vy=sqrt(tmp2/tmp1)*vy
    vz=sqrt(tmp2/tmp1)*vz
    ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

  end subroutine initialize3dvelhighmemdealiasedwithphi


  subroutine initialize3dvelhighmem(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    use mconstant
    use msymmetrize
    implicit none
  
    integer,  intent(in) :: lx1,ly,lz,iseed
    real(sp), intent(in) :: akp, u0
    complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
    real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2, kx, ky, kz
  
    real(sp), dimension(lx1,ly,lz) :: tmp
  
    integer  :: idum,ii,jj,kk,lx2
    real(sp) :: tmp1,tmp2

    lx2 = (lx1-1)*(lx1-1)
  
    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vx(ii,jj,kk)=gasdev(tmp1,tmp2)
      
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vy(ii,jj,kk)=gasdev(tmp1,tmp2)
  
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vz(ii,jj,kk)=gasdev(tmp1,tmp2)
  
    end do
    end do
    end do
  
    !call projection(vx,vy,vz,kx,ky,kz,lx1,ly,lz)
    tmp = (kx*real(vy) + ky*real(vz) + kz*real(vx))/k2
    vy = cmplx(real(vy) - kx*tmp, aimag(vy))
    vz = cmplx(real(vz) - ky*tmp, aimag(vz))
    vx = cmplx(real(vx) - kz*tmp, aimag(vx))
    tmp = (kx*aimag(vy) + ky*aimag(vz) + kz*aimag(vx))/k2
    vy = cmplx(real(vy), aimag(vy) - kx*tmp)
    vz = cmplx(real(vz), aimag(vz) - ky*tmp)
    vx = cmplx(real(vx), aimag(vx) - kz*tmp)
  
    call symmetrize(vx,lx1,ly,lz)
    call symmetrize(vy,lx1,ly,lz)
    call symmetrize(vz,lx1,ly,lz) 

    where(k2 .ge. lx2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
  
    ! different initial conditions can be imposed here
  
    tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 
  
    ! Another initial condition
    !where ( k2 .le. akp*akp)
    !    tmp = u0
    !else
    !    tmp = 0.
    !endwhere
  
    vx = vx * tmp
    vy = vy * tmp
    vz = vz * tmp
  
    tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
    tmp(1,:,:)=.5_sp*tmp(1,:,:)
    tmp1=sum(tmp)
    tmp2=(3._sp/2._sp)*u0*u0
    vx=sqrt(tmp2/tmp1)*vx
    vy=sqrt(tmp2/tmp1)*vy
    vz=sqrt(tmp2/tmp1)*vz
    ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

  end subroutine initialize3dvelhighmem

  subroutine initialize3dvellowmem(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0)
    use mconstant
    use msymmetrize
    implicit none
  
    integer,  intent(in) :: lx1,ly,lz,iseed
    real(sp), intent(in) :: akp, u0
    complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
    real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2
    real(sp),    dimension(lx1),       intent(in)    :: kx
    real(sp),    dimension(ly),        intent(in)    :: ky
    real(sp),    dimension(lz),        intent(in)    :: kz
  
    real(sp), dimension(lx1,ly,lz) :: tmp
  
    integer  :: idum,ii,jj,kk,lx2
    real(sp) :: tmp1,tmp2,kxii,kyjj,kzkk
    complex(sp) :: vxt, vyt, vzt

    lx2 = (lx1-1)*(lx1-1)
  
    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do kk=1,lz
    kzkk = kz(kk)
    do jj=1,ly
    kyjj = ky(jj)
    do ii=1,lx1
    kxii = kx(ii)
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vxt=gasdev(tmp1,tmp2)
      
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vyt=gasdev(tmp1,tmp2)
  
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vzt=gasdev(tmp1,tmp2)
  
      tmp1=kxii*vxt+kyjj*vyt+kzkk*vzt
      tmp1=tmp1/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
      vx(ii,jj,kk)=vxt-kxii*tmp1
      vy(ii,jj,kk)=vyt-kyjj*tmp1
      vz(ii,jj,kk)=vzt-kzkk*tmp1
    end do
    end do
    end do
  
    !call projection(vx,vy,vz,kx,ky,kz,lx1,ly,lz)
  
    call symmetrize(vx,lx1,ly,lz)
    call symmetrize(vy,lx1,ly,lz)
    call symmetrize(vz,lx1,ly,lz) 

    where(k2 .ge. lx2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
  
    ! different initial conditions can be imposed here
  
    tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 
  
    ! Another initial condition
    !where ( k2 .le. akp*akp)
    !    tmp = u0
    !else
    !    tmp = 0.
    !endwhere
  
    vx = vx * tmp
    vy = vy * tmp
    vz = vz * tmp
  
    tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
    tmp(1,:,:)=.5_sp*tmp(1,:,:)
    tmp1=sum(tmp)
    tmp2=(3._sp/2._sp)*u0*u0
    vx=sqrt(tmp2/tmp1)*vx
    vy=sqrt(tmp2/tmp1)*vy
    vz=sqrt(tmp2/tmp1)*vz
    ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

  end subroutine initialize3dvellowmem

  subroutine initialize3dvelhighmemdealiased(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    use mconstant
    use msymmetrize
    implicit none
  
    integer,  intent(in) :: lx1,ly,lz,iseed
    real(sp), intent(in) :: akp, u0, kcut
    complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
    real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2,kx,ky,kz
  
    real(sp), dimension(lx1,ly,lz) :: tmp
  
    integer  :: idum,ii,jj,kk
    real(sp) :: tmp1,tmp2, kcut2, kxii, kyjj, kzkk
    complex(sp) :: tmpz

  
    kcut2 = kcut * kcut

    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do kk=1,lz
    do jj=1,ly
    do ii=1,lx1

      kzkk=kz(ii,jj,kk)
      kyjj=ky(ii,jj,kk)
      kxii=kx(ii,jj,kk)

      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vx(ii,jj,kk)=gasdev(tmp1,tmp2)
      
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vy(ii,jj,kk)=gasdev(tmp1,tmp2)
  
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vz(ii,jj,kk)=gasdev(tmp1,tmp2)
  
      tmpz=kxii*vx(ii,jj,kk)+kyjj*vy(ii,jj,kk)+kzkk*vz(ii,jj,kk)
      tmpz=tmpz/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
      vx(ii,jj,kk)=vx(ii,jj,kk)-kxii*tmpz
      vy(ii,jj,kk)=vy(ii,jj,kk)-kyjj*tmpz
      vz(ii,jj,kk)=vz(ii,jj,kk)-kzkk*tmpz
    end do
    end do
    end do

    call symmetrize(vx,lx1,ly,lz)
    call symmetrize(vy,lx1,ly,lz)
    call symmetrize(vz,lx1,ly,lz) 

    where(k2 .ge. kcut2)
        vx = 0.
        vy = 0.
        vz = 0.
    endwhere
  
    ! different initial conditions can be imposed here
  
    tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 
  
    ! Another initial condition
    !where ( k2 .le. akp*akp)
    !    tmp = u0
    !else
    !    tmp = 0.
    !endwhere
  
    vx = vx * tmp
    vy = vy * tmp
    vz = vz * tmp
  
    tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
    tmp(1,:,:)=.5_sp*tmp(1,:,:)
    tmp1=sum(tmp)
    tmp2=(3._sp/2._sp)*u0*u0
    vx=sqrt(tmp2/tmp1)*vx
    vy=sqrt(tmp2/tmp1)*vy
    vz=sqrt(tmp2/tmp1)*vz
    ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

  end subroutine initialize3dvelhighmemdealiased

  subroutine initialize3dvellowmemdealiased(vx,vy,vz,kx,ky,kz,k2,iseed,lx1,ly,lz,akp,u0,kcut)
    use mconstant
    use msymmetrize
    implicit none
  
    integer,  intent(in) :: lx1,ly,lz,iseed
    real(sp), intent(in) :: akp, u0, kcut
    complex(sp), dimension(lx1,ly,lz), intent(inout) :: vx,vy,vz
    real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2
    real(sp),    dimension(lx1),       intent(in)    :: kx
    real(sp),    dimension(ly),        intent(in)    :: ky
    real(sp),    dimension(lz),        intent(in)    :: kz
  
    real(sp), dimension(lx1,ly,lz) :: tmp
  
    integer  :: idum,ii,jj,kk
    real(sp) :: tmp1,tmp2, kcut2, kxii, kyjj, kzkk
    complex(sp) :: vxt, vyt, vzt

  
    kcut2 = kcut * kcut

    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do kk=1,lz
    kzkk = kz(kk)
    do jj=1,ly
    kyjj = ky(jj)
    do ii=1,lx1
    kxii = kx(ii)
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vxt=gasdev(tmp1,tmp2)
      
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vyt=gasdev(tmp1,tmp2)
  
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      vzt=gasdev(tmp1,tmp2)
  
      tmp1=kxii*vxt+kyjj*vyt+kzkk*vzt
      tmp1=tmp1/(kxii*kxii+kyjj*kyjj+kzkk*kzkk+mytiny)
      vx(ii,jj,kk)=vxt-kxii*tmp1
      vy(ii,jj,kk)=vyt-kyjj*tmp1
      vz(ii,jj,kk)=vzt-kzkk*tmp1
    end do
    end do
    end do
  
    !call projection(vx,vy,vz,kx,ky,kz,lx1,ly,lz)
  
    call symmetrize(vx,lx1,ly,lz)
    call symmetrize(vy,lx1,ly,lz)
    call symmetrize(vz,lx1,ly,lz) 

    where(k2 .ge. kcut2)
        vx = 0._sp
        vy = 0._sp
        vz = 0._sp
    endwhere
  
    ! different initial conditions can be imposed here
  
    tmp = u0 * sqrt(8._sp*sqrt(2._sp/pi)/(3._sp*pi*akp**5)) * sqrt(k2) * exp (-k2/(akp*akp)) 
  
    ! Another initial condition
    !where ( k2 .le. akp*akp)
    !    tmp = u0
    !else
    !    tmp = 0.
    !endwhere
  
    vx = vx * tmp
    vy = vy * tmp
    vz = vz * tmp
  
    tmp=vx*conjg(vx)+vy*conjg(vy)+vz*conjg(vz)
    tmp(1,:,:)=.5_sp*tmp(1,:,:)
    tmp1=sum(tmp)
    tmp2=(3._sp/2._sp)*u0*u0
    vx=sqrt(tmp2/tmp1)*vx
    vy=sqrt(tmp2/tmp1)*vy
    vz=sqrt(tmp2/tmp1)*vz
    ! the rescaling is to make the mean kinetic energy equal 3/2 u0^2.

  end subroutine initialize3dvellowmemdealiased

  subroutine initialize2dvelhighmem(wz, k2, iseed, lx1, ly, wzrms, kcut)
    use mconstant
    use msymmetrize

    integer,  intent(in) :: lx1,ly,iseed
    real(sp), intent(in) :: wzrms, kcut
    complex(sp), dimension(lx1,ly), intent(inout) :: wz
    real(sp),    dimension(lx1,ly), intent(in)    :: k2
  
    real(sp), dimension(lx1,ly) :: tmp
  
    integer  :: idum,ii,jj
    real(sp) :: tmp1, tmp2, const, kcut2
  
    
    const = 1./(ly*(lx1-1)*2)
    kcut2 = kcut * kcut 
  
    ! this is used to initialize ran1 with a given seed iseed.
    idum=iseed
    do jj=1,ly
    do ii=1,lx1
      tmp1 = ran1(idum)
      tmp2 = ran1(idum)
      tmp2 = 2._sp*pi*tmp2
      wz(ii,jj)=gasdev(tmp1,tmp2)
    end do
    end do
  
    call symmetrize(wz,lx1,ly)

    where(k2 .ge. kcut2) wz = 0.
  
    tmp = wz*conjg(wz)
    tmp(1,:)=.5_sp*tmp(1,:)
    tmp1=sum(tmp)
    tmp2= wzrms * wzrms  / 2._sp
    wz=sqrt(tmp2/tmp1)*wz
    ! TODO: Rescaling to match rms wz. Mean wz is assumed to be ZERO.
  end subroutine initialize2dvelhighmem
  
  
  complex(sp) function gasdev(t1,t2)
    use mconstant
    real(sp) :: t1, t2
  
    gasdev=sqrt(-2.0_sp*log(t1))*cmplx(cos(t2),sin(t2))/sqrt(2.0_sp)
    ! devided by sqrt(2) to set the variance of the complex velocity to be 1.
  end function gasdev
  
  !  Random generator from Numerical Recipes, p.271
  real(sp) function ran1(idum)
    USE mconstant
    implicit none

    integer :: idum,ia,im,iq,ir,ntab,ndiv
    real(sp) :: am,eps,rnmx
    parameter(ia = 16807, im = 2147483647, am = 1._SP/im, iq = 127773, ir = 2836,  &
              ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7_SP, rnmx = 1._SP - eps)
    integer :: j,k
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
    ran1 = min(am*iy,rnmx)
  
  end function ran1
  
end module minitialize  
