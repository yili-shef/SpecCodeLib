C
C Simple example fortran program to write a
C binary datafile for tecplot.  This example
C does the following:
C
C   1.  Open a datafile called "t.plt"
C   2.  Assign values for X,Y, and P
C   3.  Write out a zone dimensioned 4x5
C   4.  Close the datafile.
C
C
      program test
      USE TecIOInterface

      character*1 NULLCHR
      Integer*4   Debug,III,NPts,NElm

      DIMENSION X1(4),P1(4), XP1(2,4)
      Dimension X(4,5), Y(4,5), P(4,5)
      Integer*4 TecIni,TecDat,TecZne,TecNod,TecFil
      Integer*4 VIsDouble
      DIMENSION X3(4,5,3), Y3(4,5,3), Z3(4,5,3), P3(4,5,3)

      NULLCHR = CHAR(0)
      Debug   = 1
      VIsDouble = 0
      IMax    = 4
      JMax    = 5
      KMax    = 1
C
C... Open the file and write the tecplot datafile 
C... header information.
C
      I = TecIni('SIMPLE DATASET'//NULLCHR,
     &           'X Y P'//NULLCHR,
     &           't.plt'//NULLCHR,
     &           '.'//NULLCHR,
     &           Debug,
     &           VIsDouble)

      I = TecIni('SIMPLE 3D DATASET'//NULLCHR,
     &           'X3 Y3 Z3 P3'//NULLCHR,
     &           't3.plt'//NULLCHR,
     &           '.'//NULLCHR,
     &           Debug,
     &           VIsDouble)

      III=1
      I=TecFil(III)
      
      Do 10 I = 1,4
      Do 10 J = 1,5
        X(I,J) = I
        Y(I,J) = J
        P(I,J) = I*J
   10 Continue
C
C... Write the zone header information.
C
      I = TecZne('Simple Zone'//NULLCHR, 
     &           IMax,
     &           JMax,
     &           KMax,
     &           'BLOCK'//NULLCHR,
     &           NULLCHR)
C
C... Write out the field data.
C
      III = IMax*JMax
      I   = TecDat(III,X,0)
      I   = TecDat(III,Y,0)
      I   = TecDat(III,P,0)

      III=2
      I=TecFil(III)
      
      KMax=3

      DO I=1,4
      DO J=1,5
      DO K=1,3
        X3(I,J,K)=I
        Y3(I,J,K)=J
        Z3(I,J,K)=K
        P3(I,J,K)=I*J*K
      END DO
      END DO
      END DO

      I = TecZne('Simple Zone'//NULLCHR, 
     &           IMax,
     &           JMax,
     &           KMax,
     &           'BLOCK'//NULLCHR,
     &           NULLCHR)
      
      III = IMax*JMax*KMax
      I   = TecDat(III,X3,0)
      I   = TecDat(III,Y3,0)
      I   = TecDat(III,Z3,0)
      I   = TecDat(III,P3,0)


      I = TecEnd()

      III=1
      I=TecFil(III)
      I=TecEnd()


      I = TecIni('SIMPLE 1D DATASET'//NULLCHR,
     &           'X1 P1'//NULLCHR,
     &           't1.plt'//NULLCHR,
     &           '.'//NULLCHR,
     &           Debug,
     &           VIsDouble)

      JMax=1
      KMax=1

      I = TecZne('Simple Zone'//NULLCHR, 
     &           IMax,
     &           JMax,
     &           KMax,
     &           'BLOCK'//NULLCHR,
     &           NULLCHR)

      DO I=1,IMax
      X1(I)=I
      P1(I)=I*I
      END DO
      III=IMax
      I =TECDat(III,X1,0)
      I =TecDat(III,P1,0)

      I=TecEnd()

      I = TecIni('SIMPLE 1D DATASET'//NULLCHR,
     &           'X1 P1'//NULLCHR,
     &           't1p.plt'//NULLCHR,
     &           '.'//NULLCHR,
     &           Debug,
     &           VIsDouble)

      JMax=1
      KMax=1

      I = TecZne('Simple Zone'//NULLCHR, 
     &           IMax,
     &           JMax,
     &           KMax,
     &           'POINT'//NULLCHR,
     &           NULLCHR)

      DO I=1,IMAX
      XP1(1,I)=I
      XP1(2,I)=I*I
      END DO
      III=IMax*2
      I =TecDat(III,XP1,0)

      I=TecEnd()
      Stop
      
      End
