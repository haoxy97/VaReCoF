      real*8 function cc_caspt2_corr(natoms,x,rparm,iparm)
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

      call cc_caspt2(rcc,cc_caspt2_corr)
      return 
      end

      subroutine cc_caspt2(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(caspt2/cc-pvdz)
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
 
      data fpp(  1), fpp(  2) / -1.8894066309D-01 , -1.3812117383D-01 /
      data fpp(  3), fpp(  4) / -8.1339641605D-02 , -4.5425259753D-02 /
      data fpp(  5), fpp(  6) / -2.2304319384D-02 , -8.8174627108D-03 /
      data fpp(  7), fpp(  8) / -1.4658297727D-03 ,  2.3057818016D-03 /
      data fpp(  9), fpp( 10) /  3.4927025663D-03 ,  3.1334079333D-03 /
      data fpp( 11), fpp( 12) /  2.3636657006D-03 ,  1.7669292642D-03 /
      data fpp( 13), fpp( 14) /  1.3386172426D-03 ,  8.4360176549D-04 /
      data fpp( 15), fpp( 16) /  2.6697569547D-04 , -3.2150454739D-04 /
      data fpp( 17), fpp( 18) / -8.4095750593D-04 , -1.1863752085D-03 /
      data fpp( 19), fpp( 20) / -9.0104166002D-04 , -5.5194985713D-04 /
      data fpp( 21), fpp( 22) / -4.5755891146D-04 , -3.9461449703D-04 /
      data fpp( 23), fpp( 24) / -2.0558310042D-04 ,  4.4856549774D-05 /
      data fpp( 25), fpp( 26) /  8.1569013226D-06 , -9.5489788548D-06 /
      data fpp( 27), fpp( 28) / -1.1080046989D-06 ,  4.4270023495D-06 /
 
      data f(  1), f(  2) / -2.7202900000D-02 , -1.5181600000D-02 /
      data f(  3), f(  4) / -8.6454000000D-03 , -5.5019000000D-03 /
      data f(  5), f(  6) / -4.2607000000D-03 , -3.9759000000D-03 /
      data f(  7), f(  8) / -4.0847000000D-03 , -4.2760000000D-03 /
      data f(  9), f( 10) / -4.3923000000D-03 , -4.3792000000D-03 /
      data f( 11), f( 12) / -4.2435000000D-03 , -4.0121000000D-03 /
      data f( 13), f( 14) / -3.7089000000D-03 , -3.3526000000D-03 /
      data f( 15), f( 16) / -2.9631000000D-03 , -2.5630000000D-03 /
      data f( 17), f( 18) / -2.1753000000D-03 , -1.5031000000D-03 /
      data f( 19), f( 20) / -1.0039000000D-03 , -5.7760000000D-04 /
      data f( 21), f( 22) / -2.9990000000D-04 , -1.3790000000D-04 /
      data f( 23), f( 24) / -6.9300000000D-05 , -6.0300000000D-05 /
      data f( 25), f( 26) / -5.4300000000D-05 , -1.7400000000D-05 /
      data f( 27), f( 28) / -3.3000000000D-06 ,  0.0000000000D+00 /
 
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
