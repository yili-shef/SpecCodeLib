PROGRAM init
  USE mprmtr
  USE mfftwplan
  IMPLICIT NONE

  INTEGER, PARAMETER :: nxlb=128,nylb=128,nzlb=128
  COMPLEX(SP), ALLOCATABLE, DIMENSION(lx1,ly,lz) :: ux,uy,uz
  REAL(SP),    ALLOCATABLE, DIMENSION(lx1,ly,lz) :: uxlb,uylb,uzlb
  REAL(SP) :: ignore_me

  OPEN(90,FILE='parameter.d',STATUS='old')
    READ(90,*)
    READ(90,*) nx
    READ(90,*) ny
    READ(90,*) nz
  CLOSE(90)

  lx=nx/2
  ly=ny
  lz=nz
  lx1=lx+1

  ALLOCATE(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  ALLOCATE(uxlb(nxlb,nylb,nylb),uylb(nxlb,nylb,nzlb),uzlb(nxlb,nylb,nzlb))
  
  OPEN(16, FILE='vel0001.dat', FORM='unformatted')
  READ(16) ux
  READ(16) uy
  READ(16) uz
  CLOSE(16)
  CALL setfftwplans
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ux,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,uz,ignore_me)

  uxlb(1:nxlb:2,:,:)=REAL(ux(1:nxlb/2,1:nylb,1:nzlb))
  uxlb(2:nxlb:2,:,:)=AIMAG(ux(1:nxlb/2,1:nylb,1:nzlb))
  uylb(1:nxlb:2,:,:)=REAL(uy(1:nxlb/2,1:nylb,1:nzlb))
  uylb(2:nxlb:2,:,:)=AIMAG(uy(1:nxlb/2,1:nylb,1:nzlb))
  uzlb(1:nxlb:2,:,:)=REAL(uz(1:nxlb/2,1:nylb,1:nzlb))
  uzlb(2:nxlb:2,:,:)=AIMAG(uz(1:nxlb/2,1:nylb,1:nzlb))
  
  OPEN(16, FILE='LB_init_DNS.dat', FORM='unformatted')
  WRITE(16) uxlb
  WRITE(16) uylb
  WRITE(16) uzlb 
  CLOSE(16)

  DEALLOCATE(ux,uy,uz,uxlb,uylb,uzlb)
  CALL destroyplans
END PROGRAM init
