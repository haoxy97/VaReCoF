      real*8 function cc_000_000_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms (not used)
c     x = cartesian coordinates 
c     e = energy correction in au
c     iparm(1) = number of the first carbon atom in the active cc bond
c     iparm(2) = number of the second carbon atom in the active cc bond
c     rparm (not used)
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_000(rcc,cc_000_000_corr)
      return 
      end

      subroutine cc_000_000(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c     
c         umbrella angle B      = 90.0 
c         deviation from planar =  0.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -2.2220713048D-01 , -1.5660823904D-01 /
      data fpp(  3), fpp(  4) / -7.7714913348D-02 , -4.5352107565D-02 /
      data fpp(  5), fpp(  6) / -1.7206656392D-02 , -3.7062668653D-03 /
      data fpp(  7), fpp(  8) /  2.1667238535D-03 ,  3.0793714513D-03 /
      data fpp(  9), fpp( 10) /  1.5657903413D-03 , -7.9253281646D-04 /
      data fpp( 11), fpp( 12) / -2.6406590755D-03 , -3.6898308817D-03 /
      data fpp( 13), fpp( 14) / -4.1850173978D-03 , -4.4550995272D-03 /
      data fpp( 15), fpp( 16) / -4.6945844935D-03 , -4.8465624988D-03 /
      data fpp( 17), fpp( 18) / -4.8841655114D-03 , -4.3629722163D-03 /
      data fpp( 19), fpp( 20) / -3.3664456233D-03 , -2.2530179832D-03 /
      data fpp( 21), fpp( 22) / -1.5510824439D-03 , -1.0826522410D-03 /
      data fpp( 23), fpp( 24) / -6.2230859202D-04 , -9.6748103430D-05 /
      data fpp( 25), fpp( 26) / -2.9898994256D-05 , -1.0728965517D-05 /
      data fpp( 27), fpp( 28) / -4.4374521074D-06 ,  1.1525226054D-05 /
 
      data f(  1), f(  2) / -6.7230600000D-02 , -5.2468100000D-02 /
      data f(  3), f(  4) / -4.3881300000D-02 , -3.8713300000D-02 /
      data f(  5), f(  6) / -3.5387500000D-02 , -3.2847600000D-02 /
      data f(  7), f(  8) / -3.0506800000D-02 , -2.8112400000D-02 /
      data f(  9), f( 10) / -2.5611000000D-02 , -2.3052600000D-02 /
      data f( 11), f( 12) / -2.0522500000D-02 , -1.8092700000D-02 /
      data f( 13), f( 14) / -1.5806800000D-02 , -1.3686800000D-02 /
      data f( 15), f( 16) / -1.1744800000D-02 , -9.9900000000D-03 /
      data f( 17), f( 18) / -8.4283000000D-03 , -5.8766000000D-03 /
      data f( 19), f( 20) / -4.0103000000D-03 , -2.4217000000D-03 /
      data f( 21), f( 22) / -1.4135000000D-03 , -8.0280000000D-04 /
      data f( 23), f( 24) / -4.6310000000D-04 , -2.0120000000D-04 /
      data f( 25), f( 26) / -1.1250000000D-04 , -3.4300000000D-05 /
      data f( 27), f( 28) / -7.2000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p05_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p05(rcc,cc_000_p05_corr)
      return 
      end

      subroutine cc_000_p05(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 95.0 
c         deviation from planar = +5.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) /  6.5739333455D-03 , -2.9210366691D-02 /
      data fpp(  3), fpp(  4) / -1.0946746658D-01 , -3.4069766983D-02 /
      data fpp(  5), fpp(  6) / -1.9573465485D-02 , -4.4413710752D-03 /
      data fpp(  7), fpp(  8) /  8.5894978602D-04 ,  2.5505719311D-03 /
      data fpp(  9), fpp( 10) /  1.4487624897D-03 , -3.8062188970D-04 /
      data fpp( 11), fpp( 12) / -1.8312749308D-03 , -2.6442783869D-03 /
      data fpp( 13), fpp( 14) / -3.0866115215D-03 , -3.4442755271D-03 /
      data fpp( 15), fpp( 16) / -3.8512863701D-03 , -4.1455789924D-03 /
      data fpp( 17), fpp( 18) / -4.3763976601D-03 , -4.1142675234D-03 /
      data fpp( 19), fpp( 20) / -3.2902822465D-03 , -2.2803698941D-03 /
      data fpp( 21), fpp( 22) / -1.5898381772D-03 , -1.0970773973D-03 /
      data fpp( 23), fpp( 24) / -6.1225223384D-04 , -9.0904599869D-05 /
      data fpp( 25), fpp( 26) / -2.9129366690D-05 , -1.1009599995D-05 /
      data fpp( 27), fpp( 28) / -4.3817555567D-06 ,  1.1675377778D-05 /
 
      data f(  1), f(  2) / -5.9051800000D-02 , -4.9328900000D-02 /
      data f(  3), f(  4) / -4.1070900000D-02 , -3.6153900000D-02 /
      data f(  5), f(  6) / -3.3005700000D-02 , -3.0636200000D-02 /
      data f(  7), f(  8) / -2.8509900000D-02 , -2.6373300000D-02 /
      data f(  9), f( 10) / -2.4153300000D-02 , -2.1880200000D-02 /
      data f( 11), f( 12) / -1.9619800000D-02 , -1.7428400000D-02 /
      data f( 13), f( 14) / -1.5340300000D-02 , -1.3375100000D-02 /
      data f( 15), f( 16) / -1.1548000000D-02 , -9.8742000000D-03 /
      data f( 17), f( 18) / -8.3658000000D-03 , -5.8641000000D-03 /
      data f( 19), f( 20) / -4.0057000000D-03 , -2.4084000000D-03 /
      data f( 21), f( 22) / -1.3945000000D-03 , -7.8630000000D-04 /
      data f( 23), f( 24) / -4.5270000000D-04 , -1.9820000000D-04 /
      data f( 25), f( 26) / -1.1120000000D-04 , -3.3100000000D-05 /
      data f( 27), f( 28) / -6.7000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p10_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p10(rcc,cc_000_p10_corr)
      return 
      end

      subroutine cc_000_p10(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 100.0 
c         deviation from planar = +10.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -3.2504905980D-01 , -1.8934688041D-01 /
      data fpp(  3), fpp(  4) /  5.0076581418D-02 , -7.3524445268D-02 /
      data fpp(  5), fpp(  6) / -7.8138003453D-03 , -8.3803533507D-03 /
      data fpp(  7), fpp(  8) /  1.1521374805D-04 ,  1.0794983585D-03 /
      data fpp(  9), fpp( 10) /  9.2179281797D-04 , -2.9666963039D-04 /
      data fpp( 11), fpp( 12) / -1.2801142964D-03 , -1.8728731840D-03 /
      data fpp( 13), fpp( 14) / -2.2983929677D-03 , -2.7485549454D-03 /
      data fpp( 15), fpp( 16) / -3.2373872507D-03 , -3.6968960516D-03 /
      data fpp( 17), fpp( 18) / -4.0250285428D-03 , -3.9789663457D-03 /
      data fpp( 19), fpp( 20) / -3.2716060744D-03 , -2.3202450556D-03 /
      data fpp( 21), fpp( 22) / -1.6098137034D-03 , -1.0901001310D-03 /
      data fpp( 23), fpp( 24) / -5.9858577265D-04 , -8.6792616555D-05 /
      data fpp( 25), fpp( 26) / -2.9243761131D-05 , -1.0922408329D-05 /
      data fpp( 27), fpp( 28) / -4.3294648158D-06 ,  1.1751732408D-05 /
 
      data f(  1), f(  2) / -5.6810700000D-02 , -4.4735200000D-02 /
      data f(  3), f(  4) / -3.9542100000D-02 , -3.4766100000D-02 /
      data f(  5), f(  6) / -3.1669000000D-02 , -2.9326300000D-02 /
      data f(  7), f(  8) / -2.7258400000D-02 , -2.5236100000D-02 /
      data f(  9), f( 10) / -2.3178100000D-02 , -2.1090300000D-02 /
      data f( 11), f( 12) / -1.9012800000D-02 , -1.6983900000D-02 /
      data f( 13), f( 14) / -1.5028800000D-02 , -1.3165800000D-02 /
      data f( 15), f( 16) / -1.1413000000D-02 , -9.7895000000D-03 /
      data f( 17), f( 18) / -8.3130000000D-03 , -5.8374000000D-03 /
      data f( 19), f( 20) / -3.9808000000D-03 , -2.3801000000D-03 /
      data f( 21), f( 22) / -1.3695000000D-03 , -7.6930000000D-04 /
      data f( 23), f( 24) / -4.4280000000D-04 , -1.9440000000D-04 /
      data f( 25), f( 26) / -1.0850000000D-04 , -3.1400000000D-05 /
      data f( 27), f( 28) / -6.1000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p15_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p15(rcc,cc_000_p15_corr)
      return 
      end

      subroutine cc_000_p15(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 105.0 
c         deviation from planar = +15.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -3.1898407863D-01 , -1.8465934274D-01 /
      data fpp(  3), fpp(  4) /  5.3386449589D-02 , -7.1591455617D-02 /
      data fpp(  5), fpp(  6) / -7.6206271229D-03 , -9.2560358917D-03 /
      data fpp(  7), fpp(  8) / -1.1352293102D-03 ,  7.1953132494D-05 /
      data fpp(  9), fpp( 10) /  4.4241678022D-04 , -2.8162025337D-04 /
      data fpp( 11), fpp( 12) / -9.5093576675D-04 , -1.4346366796D-03 /
      data fpp( 13), fpp( 14) / -1.8305175147D-03 , -2.3132932614D-03 /
      data fpp( 15), fpp( 16) / -2.8663094396D-03 , -3.3814689801D-03 /
      data fpp( 17), fpp( 18) / -3.8128146400D-03 , -3.9183215900D-03 /
      data fpp( 19), fpp( 20) / -3.3026490001D-03 , -2.3572063276D-03 /
      data fpp( 21), fpp( 22) / -1.6109256896D-03 , -1.0750909141D-03 /
      data fpp( 23), fpp( 24) / -5.8311065417D-04 , -8.5522580451D-05 /
      data fpp( 25), fpp( 26) / -2.8799024022D-05 , -1.0641637710D-05 /
      data fpp( 27), fpp( 28) / -4.1951916201D-06 ,  1.1617595810D-05 /
 
      data f(  1), f(  2) / -5.5495100000D-02 , -4.3799500000D-02 /
      data f(  3), f(  4) / -3.8798800000D-02 , -3.4082800000D-02 /
      data f(  5), f(  6) / -3.0970800000D-02 , -2.8601000000D-02 /
      data f(  7), f(  8) / -2.6536400000D-02 , -2.4563300000D-02 /
      data f(  9), f( 10) / -2.2592900000D-02 , -2.0612100000D-02 /
      data f( 11), f( 12) / -1.8642200000D-02 , -1.6709100000D-02 /
      data f( 13), f( 14) / -1.4832800000D-02 , -1.3030300000D-02 /
      data f( 15), f( 16) / -1.1320800000D-02 , -9.7257000000D-03 /
      data f( 17), f( 18) / -8.2653000000D-03 , -5.7991000000D-03 /
      data f( 19), f( 20) / -3.9406000000D-03 , -2.3417000000D-03 /
      data f( 21), f( 22) / -1.3404000000D-03 , -7.5060000000D-04 /
      data f( 23), f( 24) / -4.3140000000D-04 , -1.8840000000D-04 /
      data f( 25), f( 26) / -1.0440000000D-04 , -2.9600000000D-05 /
      data f( 27), f( 28) / -5.7000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p20_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p20(rcc,cc_000_p20_corr)
      return 
      end

      subroutine cc_000_p20(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 110.0 
c         deviation from planar = +20.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -3.1248493424D-01 , -1.7999763153D-01 /
      data fpp(  3), fpp(  4) /  5.6140460342D-02 , -7.0499209841D-02 /
      data fpp(  5), fpp(  6) / -8.0086209789D-03 , -1.0191306244D-02 /
      data fpp(  7), fpp(  8) / -2.0911540467D-03 , -5.8407756952D-04 /
      data fpp(  9), fpp( 10) /  1.9746432481D-04 , -2.8077972973D-04 /
      data fpp( 11), fpp( 12) / -7.8434540590D-04 , -1.2168386467D-03 /
      data fpp( 13), fpp( 14) / -1.5483000074D-03 , -2.0249613237D-03 /
      data fpp( 15), fpp( 16) / -2.5618546978D-03 , -3.1326198853D-03 /
      data fpp( 17), fpp( 18) / -3.6276657612D-03 , -3.8894427737D-03 /
      data fpp( 19), fpp( 20) / -3.3558131441D-03 , -2.3887184624D-03 /
      data fpp( 21), fpp( 22) / -1.6069130062D-03 , -1.0572295127D-03 /
      data fpp( 23), fpp( 24) / -5.7216894316D-04 , -8.5278414208D-05 /
      data fpp( 25), fpp( 26) / -2.8117400014D-05 , -1.0108592854D-05 /
      data fpp( 27), fpp( 28) / -4.0264238102D-06 ,  1.1213711905D-05 /
 
      data f(  1), f(  2) / -5.4694500000D-02 , -4.3294500000D-02 /
      data f(  3), f(  4) / -3.8403400000D-02 , -3.3685200000D-02 /
      data f(  5), f(  6) / -3.0526100000D-02 , -2.8118500000D-02 /
      data f(  7), f(  8) / -2.6050000000D-02 , -2.4109100000D-02 /
      data f(  9), f( 10) / -2.2196400000D-02 , -2.0284200000D-02 /
      data f( 11), f( 12) / -1.8383400000D-02 , -1.6513500000D-02 /
      data f( 13), f( 14) / -1.4691600000D-02 , -1.2932600000D-02 /
      data f( 15), f( 16) / -1.1255000000D-02 , -9.6801000000D-03 /
      data f( 17), f( 18) / -8.2300000000D-03 , -5.7655000000D-03 /
      data f( 19), f( 20) / -3.9021000000D-03 , -2.3054000000D-03 /
      data f( 21), f( 22) / -1.3136000000D-03 , -7.3320000000D-04 /
      data f( 23), f( 24) / -4.1980000000D-04 , -1.8140000000D-04 /
      data f( 25), f( 26) / -9.9900000000D-05 , -2.8300000000D-05 /
      data f( 27), f( 28) / -5.6000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p25_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p25(rcc,cc_000_p25_corr)
      return 
      end

      subroutine cc_000_p25(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 115.0 
c         deviation from planar = +25.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -2.9594122457D-01 , -1.7169755086D-01 /
      data fpp(  3), fpp(  4) /  4.8921428019D-02 , -6.7773161212D-02 /
      data fpp(  5), fpp(  6) / -9.3987831717D-03 , -1.0866706101D-02 /
      data fpp(  7), fpp(  8) / -2.8143924226D-03 , -9.1572420830D-04 /
      data fpp(  9), fpp( 10) /  1.2289255786D-05 , -3.1843281484D-04 /
      data fpp( 11), fpp( 12) / -7.9355799643D-04 , -1.1123351994D-03 /
      data fpp( 13), fpp( 14) / -1.3721012058D-03 , -1.7242599774D-03 /
      data fpp( 15), fpp( 16) / -2.2458588848D-03 , -2.8223044835D-03 /
      data fpp( 17), fpp( 18) / -3.3999231813D-03 , -3.8503282143D-03 /
      data fpp( 19), fpp( 20) / -3.4062639615D-03 , -2.4321871672D-03 /
      data fpp( 21), fpp( 22) / -1.6057873698D-03 , -1.0446633535D-03 /
      data fpp( 23), fpp( 24) / -5.6595921608D-04 , -8.6190674985D-05 /
      data fpp( 25), fpp( 26) / -2.7678083977D-05 , -9.4704105763D-06 /
      data fpp( 27), fpp( 28) / -3.9132420941D-06 ,  1.0652621047D-05 /
 
      data f(  1), f(  2) / -5.4383100000D-02 , -4.3115000000D-02 /
      data f(  3), f(  4) / -3.8072300000D-02 , -3.3321500000D-02 /
      data f(  5), f(  6) / -3.0114500000D-02 , -2.7682400000D-02 /
      data f(  7), f(  8) / -2.5621500000D-02 , -2.3714200000D-02 /
      data f(  9), f( 10) / -2.1850000000D-02 , -1.9993700000D-02 /
      data f( 11), f( 12) / -1.8151100000D-02 , -1.6339200000D-02 /
      data f( 13), f( 14) / -1.4571400000D-02 , -1.2859100000D-02 /
      data f( 15), f( 16) / -1.1216900000D-02 , -9.6649000000D-03 /
      data f( 17), f( 18) / -8.2258000000D-03 , -5.7599000000D-03 /
      data f( 19), f( 20) / -3.8862000000D-03 , -2.2847000000D-03 /
      data f( 21), f( 22) / -1.2974000000D-03 , -7.2260000000D-04 /
      data f( 23), f( 24) / -4.1240000000D-04 , -1.7640000000D-04 /
      data f( 25), f( 26) / -9.6800000000D-05 , -2.8000000000D-05 /
      data f( 27), f( 28) / -5.7000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p30_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p30(rcc,cc_000_p30_corr)
      return 
      end

      subroutine cc_000_p30(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 120.0 
c         deviation from planar = +30.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -6.0587651229D-01 , -3.0874197541D-01 /
      data fpp(  3), fpp(  4) /  2.2660441394D-01 , -7.1790680352D-02 /
      data fpp(  5), fpp(  6) / -9.1316925340D-03 , -1.1817549512D-02 /
      data fpp(  7), fpp(  8) / -3.1631094173D-03 , -1.1950128185D-03 /
      data fpp(  9), fpp( 10) / -1.8683930851D-04 , -5.1762994743D-04 /
      data fpp( 11), fpp( 12) / -9.2264090176D-04 , -1.0868064455D-03 /
      data fpp( 13), fpp( 14) / -1.1801333162D-03 , -1.3776602897D-03 /
      data fpp( 15), fpp( 16) / -1.8142255250D-03 , -2.4054376104D-03 /
      data fpp( 17), fpp( 18) / -3.0690240333D-03 , -3.7689590948D-03 /
      data fpp( 19), fpp( 20) / -3.4701395874D-03 , -2.4893302096D-03 /
      data fpp( 21), fpp( 22) / -1.6205395742D-03 , -1.0389114936D-03 /
      data fpp( 23), fpp( 24) / -5.6461445154D-04 , -8.8300898589D-05 /
      data fpp( 25), fpp( 26) / -2.7781954099D-05 , -9.0036884073D-06 /
      data fpp( 27), fpp( 28) / -3.8664025762D-06 ,  1.0310701288D-05 /
 
      data f(  1), f(  2) / -4.7745800000D-02 , -3.6786800000D-02 /
      data f(  3), f(  4) / -3.6589400000D-02 , -3.2886100000D-02 /
      data f(  5), f(  6) / -2.9647400000D-02 , -2.7209600000D-02 /
      data f(  7), f(  8) / -2.5168900000D-02 , -2.3299300000D-02 /
      data f(  9), f( 10) / -2.1483900000D-02 , -1.9684900000D-02 /
      data f( 11), f( 12) / -1.7907100000D-02 , -1.6164600000D-02 /
      data f( 13), f( 14) / -1.4465100000D-02 , -1.2813500000D-02 /
      data f( 15), f( 16) / -1.1218600000D-02 , -9.6973000000D-03 /
      data f( 17), f( 18) / -8.2727000000D-03 , -5.8016000000D-03 /
      data f( 19), f( 20) / -3.9069000000D-03 , -2.2884000000D-03 /
      data f( 21), f( 22) / -1.2969000000D-03 , -7.2250000000D-04 /
      data f( 23), f( 24) / -4.1230000000D-04 , -1.7550000000D-04 /
      data f( 25), f( 26) / -9.6300000000D-05 , -2.8900000000D-05 /
      data f( 27), f( 28) / -6.4000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p35_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p35(rcc,cc_000_p35_corr)
      return 
      end

      subroutine cc_000_p35(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 125.0 
c         deviation from planar = +35.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -1.7007466119D-01 , -1.1378067761D-01 /
      data fpp(  3), fpp(  4) / -6.2102628358D-02 ,  5.7311910457D-03 /
      data fpp(  5), fpp(  6) / -3.0632135824D-02 , -6.6376477480D-03 /
      data fpp(  7), fpp(  8) / -4.9022731835D-03 , -1.0382595181D-03 /
      data fpp(  9), fpp( 10) / -6.1968874402D-04 , -8.1798550578D-04 /
      data fpp( 11), fpp( 12) / -1.1333692329D-03 , -1.1285375628D-03 /
      data fpp( 13), fpp( 14) / -9.6748051590D-04 , -9.7154037358D-04 /
      data fpp( 15), fpp( 16) / -1.2663579898D-03 , -1.8380276674D-03 /
      data fpp( 17), fpp( 18) / -2.5865313408D-03 , -3.6288921439D-03 /
      data fpp( 19), fpp( 20) / -3.5454000836D-03 , -2.5806459839D-03 /
      data fpp( 21), fpp( 22) / -1.6552159808D-03 , -1.0416900930D-03 /
      data fpp( 23), fpp( 24) / -5.6922364723D-04 , -9.2884011798D-05 /
      data fpp( 25), fpp( 26) / -2.8240305575D-05 , -8.6370773768D-06 /
      data fpp( 27), fpp( 28) / -3.9162050274D-06 ,  9.9701025137D-06 /
 
      data f(  1), f(  2) / -5.2821300000D-02 , -4.2151600000D-02 /
      data f(  3), f(  4) / -3.6063900000D-02 , -3.2352600000D-02 /
      data f(  5), f(  6) / -2.9106700000D-02 , -2.6683700000D-02 /
      data f(  7), f(  8) / -2.4674600000D-02 , -2.2847400000D-02 /
      data f(  9), f( 10) / -2.1084700000D-02 , -1.9350900000D-02 /
      data f( 11), f( 12) / -1.7650600000D-02 , -1.5993500000D-02 /
      data f( 13), f( 14) / -1.4380500000D-02 , -1.2807300000D-02 /
      data f( 15), f( 16) / -1.1274900000D-02 , -9.7950000000D-03 /
      data f( 17), f( 18) / -8.3898000000D-03 , -5.9076000000D-03 /
      data f( 19), f( 20) / -3.9760000000D-03 , -2.3218000000D-03 /
      data f( 21), f( 22) / -1.3144000000D-03 , -7.3380000000D-04 /
      data f( 23), f( 24) / -4.1950000000D-04 , -1.7780000000D-04 /
      data f( 25), f( 26) / -9.7600000000D-05 , -3.0400000000D-05 /
      data f( 27), f( 28) / -6.9000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p40_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p40(rcc,cc_000_p40_corr)
      return 
      end

      subroutine cc_000_p40(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 130.0 
c         deviation from planar = +40.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -1.8382573307D-01 , -1.2139853386D-01 /
      data fpp(  3), fpp(  4) / -6.0540131485D-02 ,  5.8090598001D-03 /
      data fpp(  5), fpp(  6) / -3.0826107716D-02 , -6.5696289364D-03 /
      data fpp(  7), fpp(  8) / -4.9953765383D-03 , -1.3938649103D-03 /
      data fpp(  9), fpp( 10) / -1.1591638207D-03 , -1.3794798071D-03 /
      data fpp( 11), fpp( 12) / -1.5279169508D-03 , -1.2688523897D-03 /
      data fpp( 13), fpp( 14) / -8.5167349046D-04 , -5.1445364848D-04 /
      data fpp( 15), fpp( 16) / -5.8551191563D-04 , -1.0284986890D-03 /
      data fpp( 17), fpp( 18) / -1.8404933284D-03 , -3.3542706703D-03 /
      data fpp( 19), fpp( 20) / -3.6424239903D-03 , -2.7294570986D-03 /
      data fpp( 21), fpp( 22) / -1.7213476153D-03 , -1.0627524401D-03 /
      data fpp( 23), fpp( 24) / -5.8444262430D-04 , -9.9095907063D-05 /
      data fpp( 25), fpp( 26) / -2.8373747450D-05 , -8.5808041178D-06 /
      data fpp( 27), fpp( 28) / -3.9814879738D-06 ,  9.8012439869D-06 /
 
      data f(  1), f(  2) / -5.2384300000D-02 , -4.1459500000D-02 /
      data f(  3), f(  4) / -3.5401100000D-02 , -3.1727700000D-02 /
      data f(  5), f(  6) / -2.8508500000D-02 , -2.6116400000D-02 /
      data f(  7), f(  8) / -2.4138300000D-02 , -2.2346500000D-02 /
      data f(  9), f( 10) / -2.0632900000D-02 , -1.8968700000D-02 /
      data f( 11), f( 12) / -1.7359200000D-02 , -1.5808100000D-02 /
      data f( 13), f( 14) / -1.4306700000D-02 , -1.2839900000D-02 /
      data f( 15), f( 16) / -1.1396400000D-02 , -9.9788000000D-03 /
      data f( 17), f( 18) / -8.6048000000D-03 , -6.1072000000D-03 /
      data f( 19), f( 20) / -4.1136000000D-03 , -2.3935000000D-03 /
      data f( 21), f( 22) / -1.3518000000D-03 , -7.5500000000D-04 /
      data f( 23), f( 24) / -4.3140000000D-04 , -1.8150000000D-04 /
      data f( 25), f( 26) / -9.9800000000D-05 , -3.1900000000D-05 /
      data f( 27), f( 28) / -7.3000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p45_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p45(rcc,cc_000_p45_corr)
      return 
      end

      subroutine cc_000_p45(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 135.0 
c         deviation from planar =  45.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -1.9250488781D-01 , -1.2431772437D-01 /
      data fpp(  3), fpp(  4) / -5.1449214688D-02 ,  5.0345831284D-03 /
      data fpp(  5), fpp(  6) / -3.0069117825D-02 , -6.4981118275D-03 /
      data fpp(  7), fpp(  8) / -5.2584348649D-03 , -2.0931487129D-03 /
      data fpp(  9), fpp( 10) / -2.0139702834D-03 , -2.1359701533D-03 /
      data fpp( 11), fpp( 12) / -2.1171491033D-03 , -1.5604334334D-03 /
      data fpp( 13), fpp( 14) / -7.7611716314D-04 , -6.0097914072D-05 /
      data fpp( 15), fpp( 16) /  3.1150881942D-04 ,  1.1906263638D-04 /
      data fpp( 17), fpp( 18) / -6.3775936492D-04 , -2.8750032234D-03 /
      data fpp( 19), fpp( 20) / -3.7884777414D-03 , -2.9746775521D-03 /
      data fpp( 21), fpp( 22) / -1.8448120500D-03 , -1.1156742479D-03 /
      data fpp( 23), fpp( 24) / -6.0929095841D-04 , -1.0609000082D-04 /
      data fpp( 25), fpp( 26) / -2.8949038315D-05 , -8.4578846469D-06 /
      data fpp( 27), fpp( 28) / -4.0743589674D-06 ,  9.6566794837D-06 /
 
      data f(  1), f(  2) / -5.1071000000D-02 , -4.0371500000D-02 /
      data f(  3), f(  4) / -3.4613500000D-02 , -3.1022700000D-02 /
      data f(  5), f(  6) / -2.7841100000D-02 , -2.5471100000D-02 /
      data f(  7), f(  8) / -2.3509900000D-02 , -2.1746200000D-02 /
      data f(  9), f( 10) / -2.0086800000D-02 , -1.8509300000D-02 /
      data f( 11), f( 12) / -1.7016300000D-02 , -1.5604400000D-02 /
      data f( 13), f( 14) / -1.4253400000D-02 , -1.2933900000D-02 /
      data f( 15), f( 16) / -1.1619100000D-02 , -1.0295600000D-02 /
      data f( 17), f( 18) / -8.9711000000D-03 , -6.4482000000D-03 /
      data f( 19), f( 20) / -4.3500000000D-03 , -2.5153000000D-03 /
      data f( 21), f( 22) / -1.4111000000D-03 , -7.8480000000D-04 /
      data f( 23), f( 24) / -4.4670000000D-04 , -1.8580000000D-04 /
      data f( 25), f( 26) / -1.0200000000D-04 , -3.3300000000D-05 /
      data f( 27), f( 28) / -7.6000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
      real*8 function cc_000_p50_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(2)

      na=iparm(1)
      nb=iparm(2)
      
      rcc = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rcc.le.2.0) rcc = 2.0
      if(rcc.ge.20.0) rcc = 20.0

      call cc_000_p50(rcc,cc_000_p50_corr)
      return 
      end

      subroutine cc_000_p50(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rcc = xi on the ch3+ch3 mep
c
c         umbrella angle A      = 90.0 
c         deviation from planar =  0.0
c
c         umbrella angle B      = 140.0 
c         deviation from planar = +50.0
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -5.8878979719D-01 , -2.9992040563D-01 /
      data fpp(  3), fpp(  4) /  2.1698141970D-01 , -6.4230273159D-02 /
      data fpp(  5), fpp(  6) / -1.1015327060D-02 , -1.1648418601D-02 /
      data fpp(  7), fpp(  8) / -4.5509985367D-03 , -3.3725872522D-03 /
      data fpp(  9), fpp( 10) / -2.8836524544D-03 , -3.0028029302D-03 /
      data fpp( 11), fpp( 12) / -2.7451358250D-03 , -1.9616537700D-03 /
      data fpp( 13), fpp( 14) / -8.6824909495D-04 ,  3.1965014982D-04 /
      data fpp( 15), fpp( 16) /  1.2846484957D-03 ,  1.6217558675D-03 /
      data fpp( 17), fpp( 18) /  1.1383280344D-03 , -2.0396120371D-03 /
      data fpp( 19), fpp( 20) / -4.0161298861D-03 , -3.3490427802D-03 /
      data fpp( 21), fpp( 22) / -2.0420989930D-03 , -1.2049612477D-03 /
      data fpp( 23), fpp( 24) / -6.4285601604D-04 , -1.1375132802D-04 /
      data fpp( 25), fpp( 26) / -3.0138671894D-05 , -8.4083203083D-06 /
      data fpp( 27), fpp( 28) / -4.2131510426D-06 ,  9.7030755213D-06 /
 
      data f(  1), f(  2) / -4.4441800000D-02 , -3.3826600000D-02 /
      data f(  3), f(  4) / -3.3688000000D-02 , -3.0190900000D-02 /
      data f(  5), f(  6) / -2.7033500000D-02 , -2.4675700000D-02 /
      data f(  7), f(  8) / -2.2732300000D-02 , -2.1010400000D-02 /
      data f(  9), f( 10) / -1.9428000000D-02 , -1.7965000000D-02 /
      data f( 11), f( 12) / -1.6619600000D-02 , -1.5380500000D-02 /
      data f( 13), f( 14) / -1.4217800000D-02 , -1.3089200000D-02 /
      data f( 15), f( 16) / -1.1949300000D-02 , -1.0762200000D-02 /
      data f( 17), f( 18) / -9.5157000000D-03 , -6.9644000000D-03 /
      data f( 19), f( 20) / -4.7074000000D-03 , -2.6961000000D-03 /
      data f( 21), f( 22) / -1.4954000000D-03 , -8.2480000000D-04 /
      data f( 23), f( 24) / -4.6690000000D-04 , -1.9190000000D-04 /
      data f( 25), f( 26) / -1.0490000000D-04 , -3.4700000000D-05 /
      data f( 27), f( 28) / -7.9000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  2.0000000000D+00 ,  2.2000000000D+00 /
      data x(  3), x(  4) /  2.4000000000D+00 ,  2.6000000000D+00 /
      data x(  5), x(  6) /  2.8000000000D+00 ,  3.0000000000D+00 /
      data x(  7), x(  8) /  3.2000000000D+00 ,  3.4000000000D+00 /
      data x(  9), x( 10) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x( 11), x( 12) /  4.0000000000D+00 ,  4.2000000000D+00 /
      data x( 13), x( 14) /  4.4000000000D+00 ,  4.6000000000D+00 /
      data x( 15), x( 16) /  4.8000000000D+00 ,  5.0000000000D+00 /
      data x( 17), x( 18) /  5.2000000000D+00 ,  5.6000000000D+00 /
      data x( 19), x( 20) /  6.0000000000D+00 ,  6.5000000000D+00 /
      data x( 21), x( 22) /  7.0000000000D+00 ,  7.5000000000D+00 /
      data x( 23), x( 24) /  8.0000000000D+00 ,  9.0000000000D+00 /
      data x( 25), x( 26) /  1.0000000000D+01 ,  1.2000000000D+01 /
      data x( 27), x( 28) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  5), del(  6) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  7), del(  8) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  9), del( 10) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 11), del( 12) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 13), del( 14) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 15), del( 16) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del( 17), del( 18) /  4.0000000000D-01 ,  4.0000000000D-01 /
      data del( 19), del( 20) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 21), del( 22) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del( 23), del( 24) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del( 25), del( 26) /  2.0000000000D+00 ,  3.0000000000D+00 /
      data del( 27) /  5.0000000000D+00 /
      data npts / 28/
 

      IF(XI .LE. X(1)) THEN
        ii=1
      ELSE IF( XI .GE. X(npts)) THEN
        ii=npts-1
      ELSE
        call hunt(x,npts,xi,ii)
      ENDIF
      
 20      FI = FPP(II)   * (X(II+1)-XI)**3 / (6.0*DEL(II)) + 
     X     FPP(II+1) * (XI-X(II))**3   / (6.0*DEL(II)) +
     X     ((F(II+1)/DEL(II))-(FPP(II+1)*DEL(II)/6.0)) * (XI-X(II)) + 
     X     ((F(II)  /DEL(II))-(FPP(II)  *DEL(II)/6.0)) * (X(II+1)-XI)
      FI = FI + 0.00 
      RETURN
      END  
