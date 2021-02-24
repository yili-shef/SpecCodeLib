      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr) 
!                                                                       
      integer n,nm,ierr,matz 
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     of a real symmetric matrix.                                       
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrix  a.                              
!                                                                       
!        a  contains the real symmetric matrix.                         
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired.  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output                                                         
!                                                                       
!        w  contains the eigenvalues in ascending order.                
!                                                                       
!        z  contains the eigenvectors if matz is not zero.              
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for tqlrat   
!           and tql2.  the normal completion code is zero.              
!                                                                       
!        fv1  and  fv2  are temporary storage arrays.                   
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n .le. nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 if (matz .ne. 0) go to 20 
!     .......... find eigenvalues only ..........                       
      call  tred1(nm,n,a,w,fv1,fv2) 
!  tqlrat encounters catastrophic underflow on the Vax                  
!     call  tqlrat(n,w,fv2,ierr)                                        
      call  tql1(n,w,fv1,ierr) 
      go to 50 
!     .......... find both eigenvalues and eigenvectors ..........      
   20 call  tred2(nm,n,a,w,fv1,z) 
      call  tql2(nm,n,w,fv1,z,ierr) 
   50 return 
      END                                           
