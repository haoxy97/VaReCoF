
      real*8 function c2h5_h_special_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms
c     x = cartesian coordinates 
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,8),r1(3),r2(3),r3(3),r4(3),r5(3),cptmp1(3),
     $ cptmp2(3),xx(3)
      dimension iparm(3)

      na=iparm(1)
      nb=iparm(2)
      
      rch = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      if(rch.le.2.0) rch = 2.0
      if(rch.ge.20.0) rch = 20.0

      call ch_000(rch,c2h5_h_special_corr)
c     write (6,*) 'rch test',rch,c2h5_h_special_corr,iparm(1),iparm(2),
c    $ iparm(3)

c
c insert special function to restrict to certain sides
c

      r12 = 0.0d0
      r22 = 0.0d0
      r32 = 0.0d0
      r42 = 0.0d0
      r52 = 0.0d0
      do idim = 1 , 3
         xx(idim) = 0.5*(x(idim,1)+x(idim,2))
         r1(idim) = x(idim,nb) - xx(idim) 
         r2(idim) = xx(idim) - x(idim,2)
         r3(idim) = x(idim,5) - x(idim,2)
         r4(idim) = x(idim,nb) - x(idim,2)
         r5(idim) = xx(idim) - x(idim,5)
         r12 = r12 + r1(idim)**2
         r22 = r22 + r2(idim)**2
         r32 = r32 + r3(idim)**2
         r42 = r42 + r4(idim)**2
         r52 = r52 + r5(idim)**2
      enddo
      r1abs = sqrt(r12)
      r2abs = sqrt(r22)
      r3abs = sqrt(r32)
      r4abs = sqrt(r42)
      r5abs = sqrt(r52)

      cthe1 = (r12 + r22 - r42)/(2.0d0*r1abs*r2abs)
      if (dabs(cthe1).gt.1.0d0) then
          if(dabs(cthe1)-1.d0.gt.1.e-6)
     x     write (6,*) 'error in the1, cthe1 = ',cthe1
         the1 = 0.0d0
         if (cthe1.lt.1.0d0) the1 = 2.0d0*dasin(1.0d0)
      else
         the1 = dacos(cthe1)
      endif
 
      cthe1p = (r22 + r42 - r12)/(2.0d0*r2abs*r4abs)
      if (dabs(cthe1p).gt.1.0d0) then
          if(dabs(cthe1p)-1.d0.gt.1.e-6)
     x     write (6,*) 'error in the1, cthe1p = ',cthe1p
         the1p = 0.0d0
         if (cthe1p.lt.1.0d0) the1p = 2.0d0*dasin(1.0d0)
      else
         the1p = dacos(cthe1p)
      endif
 
      cthe2 = (r22 + r32 - r52)/(2.0d0*r2abs*r3abs)
      if (dabs(cthe2).gt.1.0d0) then
          if(dabs(cthe2)-1.d0.gt.1.e-6)
     x     write (6,*) 'error in the2, cthe2 = ',cthe2
         the2 = 0.0d0
         if (cthe2.lt.1.0d0) the2 = 2.0d0*dasin(1.0d0)
      else
         the2 = dacos(cthe2)
      endif
 
      sum = 0.0d0
      call cross(r1,r2,cptmp1)
      call cross(r2,r3,cptmp2)
      do idim = 1 , 3
         sum = sum - cptmp1(idim)*cptmp2(idim)
      enddo
      ctau1 = sum/(r22*r1abs*r3abs*sin(the1)*sin(the2))
      if (dabs(ctau1).gt.1.0d0) then
         write (6,*) 'error in tau1, ctau1 = ',ctau1
         write (6,*) sum,r22,r1abs,r3abs,the1,the2
         tau1 = 0.0d0
         if (ctau1.lt.1.0d0) tau1 = 2.0d0*dasin(1.0d0)
      else
         tau1 = dacos(ctau1)
      endif

      pi2 = asin(1.0d0)
      the1pdeg = the1p*90.0d0/pi2
      tau1deg = tau1*90.0d0/pi2
      if (iparm(3).eq.1) then
         if ((the1pdeg.gt.rparm).or.(tau1.lt.pi2)) 
     $    c2h5_h_special_corr=1.0d20
      endif
      if (iparm(3).eq.2) then
         if ((the1pdeg.gt.90.0d0).or.(tau1.gt.pi2)) 
     $    c2h5_h_special_corr=1.0d20
      endif
      if (iparm(3).eq.3) then
         if (((tau1.gt.pi2).and.(the1pdeg.lt.rparm)).or.
     $    ((tau1.lt.pi2).and.(the1pdeg.lt.90.0d0))) 
     $    c2h5_h_special_corr=1.0d20
      endif

      return 
      end

