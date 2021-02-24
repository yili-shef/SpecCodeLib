subroutine sgpadm( ideg,m,H,ns )
      use mconstant
      implicit none

      integer, intent(in) :: ideg,m
      real(dp), dimension(m,m), intent(inout) :: H
      integer, intent(out) :: ns

      real(dp) :: Ainf, c
      real(dp), dimension(m,m) :: D, N, X, cX

      integer, dimension(m) :: ipiv
      integer :: i, iflag
      logical :: iisodd

     
      Ainf=0.
      do i = 1, m
        Ainf = max(Ainf, sum(abs(H(i,:))))
      end do

      ns=max(floor(log(Ainf)/log(2.))+1,0)
      H=H/2**ns
      
      D=0.; N=0.; X=0.
      do i=1,m
        D(i,i)=1; N(i,i)=1; X(i,i)=1
      end do

      c=1.
      iisodd=.true.
      do i=1,ideg
        c=c*(ideg-i+1)/(i*(2*ideg-i+1))
        X=matmul(H,X)
        cX=c*X
        N=N+cX
        if(iisodd) then
                D=D-cX
        else
                D=D+cX
        end if
        iisodd= .not. iisodd
      end do
      
      call SGESV( m,m,D,m,ipiv,N,m,iflag )
      if ( iflag.ne.0 ) stop 'Problem in SGESV (within SGPADM)'
      do i=1,ns
        N=matmul(N,N)
      end do
      H=N

END
!--------------
