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

      if(rch.le.3.6) rch = 3.6
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
      dimension x( 14), f( 14),fpp( 14),del( 13)
 
      data fpp(  1), fpp(  2) / -1.6756776845D-03 , -1.9881208215D-03 /
      data fpp(  3), fpp(  4) / -2.3418390295D-03 , -3.0600023889D-03 /
      data fpp(  5), fpp(  6) / -2.8709514151D-03 , -2.2849919509D-03 /
      data fpp(  7), fpp(  8) / -1.4914807815D-03 , -4.6826168015D-04 /
      data fpp(  9), fpp( 10) / -1.1427249793D-04 ,  1.2483338691D-06 /
      data fpp( 11), fpp( 12) / -6.3708375447D-06 ,  3.8501630982D-07 /
      data fpp( 13), fpp( 14) / -9.6842276945D-07 ,  9.9921138473D-07 /
 
      data f(  1), f(  2) / -1.1439500000D-02 , -1.0022600000D-02 /
      data f(  3), f(  4) / -8.6855000000D-03 , -5.7766000000D-03 /
      data f(  5), f(  6) / -3.5949000000D-03 , -2.1144000000D-03 /
      data f(  7), f(  8) / -1.1965000000D-03 , -3.7490000000D-04 /
      data f(  9), f( 10) / -1.3310000000D-04 , -3.3300000000D-05 /
      data f( 11), f( 12) / -1.0600000000D-05 , -3.8000000000D-06 /
      data f( 13), f( 14) / -2.3000000000D-06 ,  0.0000000000D+00 /
 
      data x(  1), x(  2) /  3.6000000000D+00 ,  3.8000000000D+00 /
      data x(  3), x(  4) /  4.0000000000D+00 ,  4.5000000000D+00 /
      data x(  5), x(  6) /  5.0000000000D+00 ,  5.5000000000D+00 /
      data x(  7), x(  8) /  6.0000000000D+00 ,  7.0000000000D+00 /
      data x(  9), x( 10) /  8.0000000000D+00 ,  1.0000000000D+01 /
      data x( 11), x( 12) /  1.2000000000D+01 ,  1.4000000000D+01 /
      data x( 13), x( 14) /  1.5000000000D+01 ,  2.0000000000D+01 /
 
 
      data del(  1), del(  2) /  2.0000000000D-01 ,  2.0000000000D-01 /
      data del(  3), del(  4) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del(  5), del(  6) /  5.0000000000D-01 ,  5.0000000000D-01 /
      data del(  7), del(  8) /  1.0000000000D+00 ,  1.0000000000D+00 /
      data del(  9), del( 10) /  2.0000000000D+00 ,  2.0000000000D+00 /
      data del( 11), del( 12) /  2.0000000000D+00 ,  1.0000000000D+00 /
      data del( 13) /  5.0000000000D+00 /
      data npts / 14/
 

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

