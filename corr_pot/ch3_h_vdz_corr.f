c     implicit real*8 (a-h,o-z)
c     dimension x(3,5)

c     do i=1,3
c       do j=1,5
c         x(i,j) = 0.0d0
c       enddo
c     enddo

c     x(3,5) = 4.00d0
c     call func(5,x,e)
c     write(6,2) e
c2    format(' correction = ',1f10.6)
c     stop
c     end

c     subroutine func(na,x,e)
      real*8 function ch3_h_vdz_corr(na,x,rpar,ipar)
c
c     na = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,5)

c     if(na.ne.5) then
c       write(6,*) "ch3_h_vdz_corr: wrong number or atoms"
c       stop
c     endif
      
      rch = dsqrt( (x(1,1)-x(1,na))**2 + (x(2,1)-x(2,na))**2 + 
     x             (x(3,1)-x(3,na))**2 )

      if(rch.le.2.0) rch = 2.0
      if(rch.ge.20.) rch = 20.

      call dz_cas(rch,ch3_h_vdz_corr)
      return 
      end

      subroutine dz_cas(xi,fi)
c
c     this subroutine returns the difference:
c
c            e(cas+1+2+qc/aug-cc-pvtz) - e(cas/cc-pvdz)
c
c     for Rch = xi on the ch3+h mep
c     
c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension x( 28), f( 28),fpp( 28),del( 27)
 
      data fpp(  1), fpp(  2) / -2.0498832818D-02 , -1.3702334363D-02 /
      data fpp(  3), fpp(  4) / -4.7318297288D-03 , -3.3703467216D-03 /
      data fpp(  5), fpp(  6) / -2.9067833849D-03 , -2.7025197389D-03 /
      data fpp(  7), fpp(  8) / -2.3031376597D-03 , -1.9149296224D-03 /
      data fpp(  9), fpp( 10) / -1.7971438507D-03 , -1.9664949749D-03 /
      data fpp( 11), fpp( 12) / -2.3068762497D-03 , -2.6810000263D-03 /
      data fpp( 13), fpp( 14) / -2.9291236450D-03 , -3.0175053938D-03 /
      data fpp( 15), fpp( 16) / -2.9858547798D-03 , -2.8440754870D-03 /
      data fpp( 17), fpp( 18) / -2.6478432721D-03 , -2.1319324403D-03 /
      data fpp( 19), fpp( 20) / -1.5369269667D-03 , -9.0591696762D-04 /
      data fpp( 21), fpp( 22) / -5.0820516282D-04 , -2.8686238110D-04 /
      data fpp( 23), fpp( 24) / -1.3714531276D-04 , -2.4732871156D-05 /
      data fpp( 25), fpp( 26) / -1.2923202613D-05 , -2.1139565829D-06 /
      data fpp( 27), fpp( 28) / -1.5046763149D-06 ,  3.3153381575D-06 /
 
      data f(  1), f(  2) / -2.6909000000D-02 , -2.4108500000D-02 /
      data f(  3), f(  4) / -2.1841600000D-02 , -1.9814700000D-02 /
      data f(  5), f(  6) / -1.7928600000D-02 , -1.6160500000D-02 /
      data f(  7), f(  8) / -1.4499200000D-02 , -1.2930100000D-02 /
      data f(  9), f( 10) / -1.1439400000D-02 , -1.0022500000D-02 /
      data f( 11), f( 12) / -8.6854000000D-03 , -7.4408000000D-03 /
      data f( 13), f( 14) / -6.3026000000D-03 , -5.2805000000D-03 /
      data f( 15), f( 16) / -4.3783000000D-03 , -3.5948000000D-03 /
      data f( 17), f( 18) / -2.9247000000D-03 , -1.8911000000D-03 /
      data f( 19), f( 20) / -1.1965000000D-03 , -6.6760000000D-04 /
      data f( 21), f( 22) / -3.7490000000D-04 , -2.1660000000D-04 /
      data f( 23), f( 24) / -1.3300000000D-04 , -6.2400000000D-05 /
      data f( 25), f( 26) / -3.3300000000D-05 , -1.0600000000D-05 /
      data f( 27), f( 28) / -2.3000000000D-06 ,  0.0000000000D+00 /
 
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

