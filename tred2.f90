      subroutine tred2(nm,n,a,d,e,z) 
!                                                                       
      integer i,j,k,l,n,ii,nm,jp1 
      double precision a(nm,n),d(n),e(n),z(nm,n) 
      double precision f,g,h,hh,scale 
!                                                                       
!     this subroutine is a translation of the algol procedure tred2,    
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine reduces a real symmetric matrix to a              
!     symmetric tridiagonal matrix using and accumulating               
!     orthogonal similarity transformations.                            
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        a contains the real symmetric input matrix.  only the          
!          lower triangle of the matrix need be supplied.               
!                                                                       
!     on output                                                         
!                                                                       
!        d contains the diagonal elements of the tridiagonal matrix.    
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is set to zero.      
!                                                                       
!        z contains the orthogonal transformation matrix                
!          produced in the reduction.                                   
!                                                                       
!        a and z may coincide.  if distinct, a is unaltered.            
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      do 100 i = 1, n 
!                                                                       
         do 80 j = i, n 
   80    z(j,i) = a(j,i) 
!                                                                       
         d(i) = a(n,i) 
  100 continue 
!                                                                       
      if (n .eq. 1) go to 510 
!     .......... for i=n step -1 until 2 do -- ..........               
      do 300 ii = 2, n 
         i = n + 2 - ii 
         l = i - 1 
         h = 0.0d0 
         scale = 0.0d0 
         if (l .lt. 2) go to 130 
!     .......... scale row (algol tol then not needed) ..........       
         do 120 k = 1, l 
  120    scale = scale + dabs(d(k)) 
!                                                                       
         if (scale .ne. 0.0d0) go to 140 
  130    e(i) = d(l) 
!                                                                       
         do 135 j = 1, l 
            d(j) = z(l,j) 
            z(i,j) = 0.0d0 
            z(j,i) = 0.0d0 
  135    continue 
!                                                                       
         go to 290 
!                                                                       
  140    do 150 k = 1, l 
            d(k) = d(k) / scale 
            h = h + d(k) * d(k) 
  150    continue 
!                                                                       
         f = d(l) 
         g = -dsign(dsqrt(h),f) 
         e(i) = scale * g 
         h = h - f * g 
         d(l) = f - g 
!     .......... form a*u ..........                                    
         do 170 j = 1, l 
  170    e(j) = 0.0d0 
!                                                                       
         do 240 j = 1, l 
            f = d(j) 
            z(j,i) = f 
            g = e(j) + z(j,j) * f 
            jp1 = j + 1 
            if (l .lt. jp1) go to 220 
!                                                                       
            do 200 k = jp1, l 
               g = g + z(k,j) * d(k) 
               e(k) = e(k) + z(k,j) * f 
  200       continue 
!                                                                       
  220       e(j) = g 
  240    continue 
!     .......... form p ..........                                      
         f = 0.0d0 
!                                                                       
         do 245 j = 1, l 
            e(j) = e(j) / h 
            f = f + e(j) * d(j) 
  245    continue 
!                                                                       
         hh = f / (h + h) 
!     .......... form q ..........                                      
         do 250 j = 1, l 
  250    e(j) = e(j) - hh * d(j) 
!     .......... form reduced a ..........                              
         do 280 j = 1, l 
            f = d(j) 
            g = e(j) 
!                                                                       
            do 260 k = j, l 
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k) 
!                                                                       
            d(j) = z(l,j) 
            z(i,j) = 0.0d0 
  280    continue 
!                                                                       
  290    d(i) = h 
  300 continue 
!     .......... accumulation of transformation matrices ..........     
      do 500 i = 2, n 
         l = i - 1 
         z(n,l) = z(l,l) 
         z(l,l) = 1.0d0 
         h = d(i) 
         if (h .eq. 0.0d0) go to 380 
!                                                                       
         do 330 k = 1, l 
  330    d(k) = z(k,i) / h 
!                                                                       
         do 360 j = 1, l 
            g = 0.0d0 
!                                                                       
            do 340 k = 1, l 
  340       g = g + z(k,i) * z(k,j) 
!                                                                       
            do 360 k = 1, l 
               z(k,j) = z(k,j) - g * d(k) 
  360    continue 
!                                                                       
  380    do 400 k = 1, l 
  400    z(k,i) = 0.0d0 
!                                                                       
  500 continue 
!                                                                       
  510 do 520 i = 1, n 
         d(i) = z(n,i) 
         z(n,i) = 0.0d0 
  520 continue 
!                                                                       
      z(n,n) = 1.0d0 
      e(1) = 0.0d0 
      return 
      END                                           
