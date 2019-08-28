
c***********************************************************************

      real*8 function ch3_o(r, rpar, ipar)

c     
c     this is a general subroutine for evaluating the potential 
c     energy for the transitional modes.
c     
c     the potential should be output in au in vtot
c     the cartesian coordinates of the atoms of each fragment are 
c     stored in r
c     the bond length or some other suitable approximation to a 
c     distinguished reaction coordinate should be passed back in rmepi
c     
c     
c     the most important current utility of this routine is its 
c     evaluation of a particular set of valence bending and 
c     torsional angles for the incipient bond (qq).
c     it also generates the atom-atom interfragment distances and the 
c     center-of-mass to center-of-mass separation. 
c     it is also currently setup to consider
c     (i) valence stretching and bending potentials and/or 
c     lennard-jones and/or hard sphere interactions.
c     (ii) ion-molecule potentials including ion-induced dipole, 
c     ion-dipole, hard sphere, and/or lennard-jones interactions.
c     
c     the valence stretching and bending potential works in terms  of
c     either 
c     (i) an exponential extrapolation of an equilibrium set 
c     of force constants 
c     or
c     (ii) an exponential interpolation of force constants along the 
c     bond stretching path.
c     
c     for the ion-dipole component some work needs to be done on 
c     defining the axis for the dipole.  once this is done it should 
c     also be easy to include the ion-quadrupole interaction.
c     

      implicit double precision (a-h,o-z)

      include 'param2.fi'
      include 'commonpot.fi'
c     
c     vectors for primary coordinates 
c     

      dimension r(nfrag,natommx,ndim)
      dimension cm(nfrag,ndim),rij(natommx,natommx)
      dimension qq(5)

c     
c     vectors from mcpot calling routine
c     

      dimension ift(nfrag),natom(nfrag)
c      dimension dmss(nfrag,natommx),dmsst(nfrag)

c     
c     intermediate use vectors
c     

      dimension r1(ndim),r2p(ndim),r2d(ndim),r2q(ndim)
      dimension nr(nfrag)
      dimension rbvec(ndim),rax1v(ndim),rbx1v(ndim),rax2v(ndim),
     $     rbx2v(ndim),rx1x2v(ndim),ray1v(ndim),rby1v(ndim),
     $     ray2v(ndim),rby2v(ndim),ry1y2v(ndim)
      dimension cptmp1(ndim),cptmp2(ndim),cptmp3(ndim)

c     
c     specific potential parameter vectors
c     

      dimension fijrb(5,5),qqerb(5)
      dimension rpathv(nrpotm),vrr(nrpotm)

      include 'data.fi'

c     
c     for generality/ease of modification purposes i have set this up to 
c     begin by evaluating the interfragment atom-atom distances
c     and a particular set of valence stretching, bending and torsional 
c     angles.
c     for efficient use one might wish to only evaluate the components 
c     that are required for the particular potential being used.  for a 
c     variety of reasons this seems to me to be best left for the specific 
c     applications.
c     
      ift(1)=2
      ift(2)=0
      iel = 1
      ieff = 1
      ich = 1
      ictf1(1, 1, 1) = 1

      pi2 = asin(1.0d0)
      pi = 2.0d0*pi2

c     
c     start by evaluating bonding separation
c     

      rij(1,1) = dsqrt((r(1,1,1)-r(2,1,1))**2+(r(1,1,2)-r(2,1,2))**2+ 
     $     (r(1,1,3)-r(2,1,3))**2)
      rbond = rij(1,1)
      rmepi = rbond

c     
c     now proceed to the valence bending/torsional coordinates.
c     

c     
c     start by evaluating fragment's number of angular degrees of freedom
c     
      do 500 ifrag = 1, nfrag
         if (ift(ifrag).eq.0) nr(ifrag) = 0
         if (ift(ifrag).eq.1) then
            nr(ifrag) = 2
         endif
         if (ift(ifrag).ge.2) then
            nr(ifrag) = 3
         endif
 500  continue
c     
c     determine distances needed for bending angles
c     
      rbond2 = 0.0d0
      rax12 = 0.0d0
      rbx12 = 0.0d0
      rby12 = 0.0d0
      ray12 = 0.0d0
      do 1500 idim = 1 , ndim
         rbvec(idim) = r(2,1,idim) - r(1,1,idim)
         rbond2 = rbond2 + rbvec(idim)**2
         rax1v(idim) = r(1,2,idim) - r(1,1,idim)
         rax12 = rax12 + rax1v(idim)**2
         rbx1v(idim) = r(2,1,idim) - r(1,2,idim)
         rbx12 = rbx12 + rbx1v(idim)**2
         if (nr(2).lt.2) go to 1500
         ray1v(idim) = r(1,1,idim) - r(2,2,idim)
         ray12 = ray12 + ray1v(idim)**2
         rby1v(idim) = r(2,2,idim) - r(2,1,idim)
         rby12 = rby12 + rby1v(idim)**2
 1500 continue
      rbond0 = sqrt(rbond2)
      rax1 = sqrt(rax12)
      rbx1 = sqrt(rbx12)
      ray1 = sqrt(ray12)
      rby1 = sqrt(rby12)
c     
c     calculate the bending angles.
c     

      cthe1 = (rax12 + rbond2 - rbx12)/(2.0d0*rax1*rbond0)
      if (dabs(cthe1).le.1.0d0) go to 2030
      if(dabs(cthe1)-1.d0.gt.1.e-6) 
     x     write (6,*) 'error in the1, cthe1 = ',cthe1
      the1 = 0.0d0
      if (cthe1.lt.1.0d0) the1 = 2.0d0*dasin(1.0d0)
      go to 2050
 2030 continue
      the1 = dacos(cthe1)
 2050 continue
      
      if (nr(2).lt.2) go to 2200

      cthe2 = (rby12 + rbond2 - ray12)/(2.0d0*rby1*rbond0)
      if (dabs(cthe2).le.1.0d0) go to 2130
      if(dabs(cthe2)-1.d0.gt.1.e-6) 
     x     write (6,*) 'error in the2, cthe2 = ',cthe2
      the2 = 0.0d0
      if (cthe2.lt.1.0d0) the2 = 2.0d0*dasin(1.0d0)
      go to 2150
 2130 continue
      the2 = dacos(cthe2)
 2150 continue

 2200 continue

c     
c     calculate the torsional angles if ivb.ne.0
c     
c     
c     evaluate the number of degrees of freedom
c     
      ntdof = nr(1) + nr(2) - 1
      if (nangles(iel,ieff,ich).eq.0) go to 520
      if (ntdof.ne.nangles(iel,ieff,ich)) go to 9380
 520  continue     
c     
c     note: if nr(i) is equal to 3 then we must  
c     specify whether the `bonding' atom is a central or 
c     terminal atom in order to determine what type of torsional angle 
c     to use. 
c     ict = 1 implies a terminal atom
c     ict = 2 implies a central atom
c     
      if (nr(1).eq.3) then
         if ((ictf1(iel,ieff,ich).lt.1).or.(ictf1(iel,ieff,ich).gt.2))
     $        go to 9010
      endif
      if (nr(2).eq.3) then
         if ((ictf2(iel,ieff,ich).lt.1).or.(ictf2(iel,ieff,ich).gt.2)) 
     $        go to 9020
      endif
c     
c     determine distances needed for torsional angles
c     
      rax22 = 0.0d0
      rbx22 = 0.0d0
      rby22 = 0.0d0
      ray22 = 0.0d0
      ry1y22 = 0.0d0
      rx1x22 = 0.0d0
      do idim = 1 , ndim
         if (nr(1).ge.3) then
            rax2v(idim) = r(1,3,idim) - r(1,1,idim)
            rax22 = rax22 + rax2v(idim)**2
            rbx2v(idim) = r(2,1,idim) - r(1,3,idim)
            rbx22 = rbx22 + rbx2v(idim)**2
            rx1x2v(idim) = r(1,3,idim) - r(1,2,idim)
            rx1x22 = rx1x22 + rx1x2v(idim)**2
         endif
         if (nr(2).ge.3) then
            ray2v(idim) = r(1,1,idim) - r(2,3,idim)
            ray22 = ray22 + ray2v(idim)**2
            rby2v(idim) = r(2,3,idim) - r(2,1,idim)
            rby22 = rby22 + rby2v(idim)**2
            ry1y2v(idim) = r(2,3,idim) - r(2,2,idim)
            ry1y22 = ry1y22 + ry1y2v(idim)**2
         endif
      enddo
      rax2 = sqrt(rax22)
      rbx2 = sqrt(rbx22)
      ray2 = sqrt(ray22)
      rby2 = sqrt(rby22)
      ry1y2 = sqrt(ry1y22)
      rx1x2 = sqrt(rx1x22)
c     
c     start with tau 1 which corresponds to the torsion angle 
c     for y1,21,11,x1.
c     tau1b corresponds to the torsion angle for perp to by1,rb;
c     y1,21,11,perp.  
c     

      if (nr(2).ge.2) then

         sum = 0.0d0
         call cross(rby1v,rbvec,cptmp1)
         call cross(rbvec,rax1v,cptmp2)
         do 2500 idim = 1 , ndim
            sum = sum - cptmp1(idim)*cptmp2(idim)
 2500    continue
         ctau1 = sum/(rbond2*rby1*rax1*sin(the1)*sin(the2))
         if (dabs(ctau1).le.1.0d0) go to 2530
         write (6,*) 'error in tau1, ctau1 = ',ctau1
         tau1 = 0.0d0
         if (ctau1.lt.1.0d0) tau1 = 2.0d0*dasin(1.0d0)
         go to 2550
 2530    continue
         tau1 = dacos(ctau1)
 2550    continue
         
c     
c     still need to determine whether it is +/-
c     
         
         sum = 0.0d0
         call cross(rby1v,rbvec,cptmp1)
         call cross(rbvec,cptmp2,cptmp3)

         do 2570 idim = 1 , ndim
            sum = sum - cptmp1(idim)*cptmp3(idim)
 2570    continue
         ctau1b = sum/(rby1*rbond2*rbond0*rax1*sin(the1)*sin(the2))
         if (dabs(ctau1b).le.1.0d0) go to 2580
         write (6,*) 'error in tau1, ctau1b = ',ctau1b
         tau1b = 0.0d0
         if (ctau1b.lt.1.0d0) tau1b = 2.0d0*dasin(1.0d0)
         go to 2590
 2580    continue
         tau1b = dacos(ctau1b)
 2590    continue
         
         if (tau1b.gt.pi2) tau1 = -tau1
         
      endif

      if (nr(2).eq.3) then

c     
c     now evaluate tau 2 (torsion angle for y2,21,y1,11).
c     tau2b corresponds to torsion angle for perp to by1,rb;
c     y2,21,y1,perp
c     

         if (ictf2(iel,ieff,ich).eq.2) then 

c     
c     first evaluate internal y2,21,y1 bending angle
c     

            ctbint = (rby12 + rby22-ry1y22)/(2.0d0*rby1*rby2)
            if (dabs(ctbint).le.1.0d0) go to 2630
            write (6,*) 'error in tbint, ctbint = ',ctbint
            tbint = 0.0d0
            if (ctbint.lt.1.0d0) tbint = 2.0d0*dasin(1.0d0)
            go to 2650
 2630       continue
            tbint = dacos(ctbint)
 2650       continue
            
            sum = 0.0d0
            call cross(rby2v,rby1v,cptmp1)
            call cross(rby1v,rbvec,cptmp2)
            do 2660 idim = 1 , ndim
               sum = sum + cptmp1(idim)*cptmp2(idim)
 2660       continue
            ctau2 = sum/(rby2*rby12*rbond0*sin(the2)*sin(tbint))
            if (dabs(ctau2).le.1.0d0) go to 2670
            write (6,*) 'error in tau2, ctau2 = ',ctau2
            tau2 = 0.0d0
            if (ctau2.lt.1.0d0) tau2 = 2.0d0*dasin(1.0d0)
            go to 2680
 2670       continue
            tau2 = dacos(ctau2)
 2680       continue
            
c     
c     still need to determine whether it is +/-
c     
            
            sum = 0.0d0
            call cross(rby2v,rby1v,cptmp1)
            call cross(rby1v,cptmp2,cptmp3)
            do 2685 idim = 1 , ndim
               sum = sum + cptmp1(idim)*cptmp3(idim)
 2685       continue
            ctau2b = sum/(rby2*rby12*rby1*rbond0*sin(the2)*
     $           sin(tbint))
            if (dabs(ctau2b).le.1.0d0) go to 2690
            write (6,*) 'error in tau2b, ctau2b = ',ctau2b
            tau2b = 0.0d0
            if (ctau2b.lt.1.0d0) tau2b = 2.0d0*dasin(1.0d0)
            go to 2695
 2690       continue
            tau2b = dacos(ctau2b)
 2695       continue
            
            if (tau2b.lt.pi2) tau2 = -tau2
            
         endif

         if (ictf2(iel,ieff,ich).eq.1) then 

c     
c     first evaluate internal y2,21,y1 bending angle
c     

            ctbint = (rby12 + ry1y22 - rby22)/(2.0d0*rby1*ry1y2)
            if (dabs(ctbint).le.1.0d0) go to 2930
            write (6,*) 'error in tbint, ctbint = ',ctbint
            tbint = 0.0d0
            if (ctbint.lt.1.0d0) tbint = 2.0d0*dasin(1.0d0)
            go to 2950
 2930       continue
            tbint = dacos(ctbint)
 2950       continue
            
            sum = 0.0d0
            call cross(ry1y2v,rby1v,cptmp1)
            call cross(rby1v,rbvec,cptmp2)
            do 2960 idim = 1 , ndim
               sum = sum + cptmp1(idim)*cptmp2(idim)
 2960       continue
            ctau2 = sum/(ry1y2*rby12*rbond0*sin(the2)*sin(tbint))
            if (dabs(ctau2).le.1.0d0) go to 2970
            write (6,*) 'error in tau2, ctau2 = ',ctau2
            tau2 = 0.0d0
            if (ctau2.lt.1.0d0) tau2 = 2.0d0*dasin(1.0d0)
            go to 2980
 2970       continue
            tau2 = dacos(ctau2)
 2980       continue
            
c     
c     still need to determine whether it is +/-
c     
            
            sum = 0.0d0
            call cross(ry1y2v,rby1v,cptmp1)
            call cross(rby1v,cptmp2,cptmp3)
            do 2985 idim = 1 , ndim
               sum = sum + cptmp1(idim)*cptmp3(idim)
 2985       continue
            ctau2b = sum/(ry1y2*rby12*rby1*rbond0*sin(the2)*
     $           sin(tbint))
            if (dabs(ctau2b).le.1.0d0) go to 2990
            write (6,*) 'error in tau2b, ctau2b = ',ctau2b
            tau2b = 0.0d0
            if (ctau2b.lt.1.0d0) tau2b = 2.0d0*dasin(1.0d0)
            go to 2995
 2990       continue
            tau2b = dacos(ctau2b)
 2995       continue
            
            if (tau2b.gt.pi2) tau2 = -tau2
            
         endif

      endif

      if (nr(1).eq.3) then

c     
c     now evaluate tau 3 (torsion angle for x2,11,x1,21.
c     tau3b corresponds to torsion angle for perp to x1,21;
c     x2,21,x1,perp
c     

         if (ictf1(iel,ieff,ich).eq.2) then 

c     
c     consider the central atom case first
c     
c     first evaluate internal x2,11,x1 bending angle
c     

            ctaint = (rax12 + rax22 - rx1x22)/(2.0d0*rax1*rax2)
            if (dabs(ctaint).le.1.0d0) go to 2710
            write (6,*) 'error in taint, ctaint = ',ctaint
            taint = 0.0d0
            if (ctaint.lt.1.0d0) taint = 2.0d0*dasin(1.0d0)
            go to 2720
 2710       continue
            taint = dacos(ctaint)
 2720       continue
            
            sum = 0.0d0
            call cross(rbvec,rax1v,cptmp1)
            call cross(rax1v,rax2v,cptmp2)
            do 2730 idim = 1 , ndim
               sum = sum - cptmp1(idim)*cptmp2(idim)
 2730       continue
            ctau3 = sum/(rbond0*rax12*rax2*sin(the1)*sin(taint))
            if (dabs(ctau3).le.1.0d0) go to 2740
            write (6,*) 'error in tau3, ctau3 = ',ctau3
            tau3 = 0.0d0
            if (ctau3.lt.1.0d0) tau3 = 2.0d0*dasin(1.0d0)
            go to 2750
 2740       continue
            tau3 = dacos(ctau3)
 2750       continue
            
c     
c     still need to determine whether it is +/-
c     
            
            sum = 0.0d0
            call cross(rbvec,rax1v,cptmp1)
            call cross(rax1v,cptmp2,cptmp3)
            do 2760 idim = 1 , ndim
               sum = sum - cptmp1(idim)*cptmp3(idim)
 2760       continue
            ctau3b = sum/(rbond0*rax12*rax1*rax2*sin(the1)*
     $           sin(taint))
            if (dabs(ctau3b).le.1.0d0) go to 2770
            write (6,*) 'error in tau3b, ctau3b = ',ctau3b
            tau3b = 0.0d0
            if (ctau3b.lt.1.0d0) tau3b = 2.0d0*dasin(1.0d0)
            go to 2780
 2770       continue
            tau3b = dacos(ctau3b)
 2780       continue
            
            if (tau3b.lt.pi2) tau3 = -tau3

         endif

         if (ictf1(iel,ieff,ich).eq.1) then 

c     
c     now consider the terminal atom case.
c     
c     first evaluate internal x2,11,x1 bending angle
c     

            ctaint = (rax12 + rx1x22 - rax22)/(2.0d0*rax1*rx1x2)
            if (dabs(ctaint).le.1.0d0) go to 2810
            write (6,*) 'error in taint, ctaint = ',ctaint
            taint = 0.0d0
            if (ctaint.lt.1.0d0) taint = 2.0d0*dasin(1.0d0)
            go to 2820
 2810       continue
            taint = dacos(ctaint)
 2820       continue
            
            sum = 0.0d0
            call cross(rbvec,rax1v,cptmp1)
            call cross(rax1v,rx1x2v,cptmp2)
            do 2830 idim = 1 , ndim
               sum = sum - cptmp1(idim)*cptmp2(idim)
 2830       continue
            ctau3 = sum/(rbond0*rax12*rx1x2*sin(the1)*sin(taint))
            if (dabs(ctau3).le.1.0d0) go to 2840
            write (6,*) 'error in tau3, ctau3 = ',ctau3
            tau3 = 0.0d0
            if (ctau3.lt.1.0d0) tau3 = 2.0d0*dasin(1.0d0)
            go to 2850
 2840       continue
            tau3 = dacos(ctau3)
 2850       continue
            
c     
c     still need to determine whether it is +/-
c     
            
            sum = 0.0d0
            call cross(rbvec,rax1v,cptmp1)
            call cross(rax1v,cptmp2,cptmp3)
            do 2860 idim = 1 , ndim
               sum = sum - cptmp1(idim)*cptmp3(idim)
 2860       continue
            ctau3b = sum/(rbond0*rax12*rax1*rx1x2*sin(the1)*
     $           sin(taint))
            if (dabs(ctau3b).le.1.0d0) go to 2870
            write (6,*) 'error in tau3b, ctau3b = ',ctau3b
            tau3b = 0.0d0
            if (ctau3b.lt.1.0d0) tau3b = 2.0d0*dasin(1.0d0)
            go to 2880
 2870       continue
            tau3b = dacos(ctau3b)
 2880       continue
            
            if (tau3b.lt.pi2) tau3 = -tau3

         endif

      endif

c     
c     finished generating the valence bending and torsional angles. 
c     now need to put them in the order that one is likely to 
c     generate them in during a z-matrix calculation.
c     

      qq(1) = the1
      if (nr(1).eq.2) then
         if (nr(2).eq.2) then
            qq(2) = the2
            qq(3) = tau1
         endif
         if (nr(2).eq.3) then
            qq(2) = the2
            qq(3) = tau1
            qq(4) = tau2
         endif
      endif
      if (nr(1).eq.3) then 
         qq(2) = tau3
         if (nr(2).eq.2) then
            qq(3) = the2
            qq(4) = tau1
         endif
         if (nr(2).eq.3) then
            qq(3) = the2
            qq(4) = tau1
            qq(5) = tau2
         endif
      endif

c ***************************************************************
c
c     add up all the contributions
c
c ***************************************************************

      the = the1*180.0d0/pi
      if (the.gt.90.0d0) then
         the=180.0d0-the
c         vtot=1.0d20
c         return
      endif
      tau = tau3*180.0d0/pi
      if (rbond0.lt.3.6d0) rbond0=3.6d0

      call ch3po_1ap(rbond0,the,tau,vtot1ap)
      ch3_o = vtot1ap*ckctoerg/cautoerg
c     call ch3po_1app(rbond0,the,tau,vtot1app)
c     vtot1app = vtot1app*ckctoerg/cautoerg
c     call ch3po_2app(rbond0,the,tau,vtot2app)
c     vtot2app = vtot2app*ckctoerg/cautoerg
c     write (6,9001) rbond0*cautoang,the1*180.0d0/pi,
c    $ tau3*180.0d0/pi,vtot1ap*cautoicm,vtot1app*cautoicm,
c    $ vtot2app*cautoicm
c9001 format (1x,'vtot test',6g12.5)
c     vtot = min(vtot1ap,vtot1app,vtot2app)
c     if (iel.eq.1) vtot=vtot1ap
c     if (iel.eq.2) vtot=vtot1app
c     if (iel.eq.3) vtot=vtot2app

      return

c
c     error messages
c

 9010 continue
      write (6,*) 'error in potcalc.f.  must designate atom ',
     $     'type in potential input.'
      go to 9900

 9020 continue
      write (6,*) 'error in potcalc.f.  must designate atom ',
     $     'type in potential input.'
      go to 9900

 9100 continue
      write (6,*) 'error in potcalc.f. incorrect stretching',
     $     ' or bending parameters in potential input.'
      go to 9900

 9200 continue
      write (6,*) 'error in potcalc.f.  incorrect atom-atoma',
     $     ' inverse 6 or anisotropy parameters in potential input.'
      go to 9900

 9300 continue
      write (6,*) 'error in potcalc.f.  incorrect ion-',
     $     'molecule parameters in potential input.'
      go to 9900

 9350 write (6,*) 'error in potcalc.f.  distance coordinate out',
     $     ' of range of interpolation points.'
      go to 9900

 9380 write (6,*) 'error in potcalc.f.  input number of angles',
     $     ' nanglesd not equal to calculated degrees of freedom'
      go to 9900 

 9400 continue

 9410 continue

 9900 continue
      stop
      end



      subroutine ch3po_1ap(r,theta,phi,energy)
      implicit real*8 (a-h,o-z)
      data degrad / 57.29577951308232d 00 /
      call ch3po_1ap_s(r,theta,es)
      call ch3po_1ap_d(r,theta,ed)
      energy = es + 0.5d0 * cos(3.0d0*phi/degrad)*ed 
      return
      end

      subroutine ch3po_1ap_s(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     1(2a') surface
c     average of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -2.91250000d+01 , -3.63120000d+01 /
      data f( 1, 3),f( 1, 4) / -4.06010000d+01 , -4.20220000d+01 /
      data f( 1, 5),f( 1, 6) / -4.06010000d+01 , -3.63120000d+01 /
      data f( 1, 7),f( 1, 8) / -2.91250000d+01 , -1.91900000d+01 /
      data f( 1, 9),f( 1,10) / -7.37700000d+00 ,  5.24000000d+00 /
      data f( 1,11),f( 1,12) /  1.91780000d+01 ,  3.97640000d+01 /
      data f( 1,13),f( 1,14) /  3.97640000d+01 ,  1.91780000d+01 /
      data f( 1,15) /           5.24000000d+00 /
      data f( 2, 1),f( 2, 2) / -2.69710000d+01 , -3.35520000d+01 /
      data f( 2, 3),f( 2, 4) / -3.74660000d+01 , -3.87620000d+01 /
      data f( 2, 5),f( 2, 6) / -3.74660000d+01 , -3.35520000d+01 /
      data f( 2, 7),f( 2, 8) / -2.69710000d+01 , -1.77890000d+01 /
      data f( 2, 9),f( 2,10) / -6.58100000d+00 ,  4.99500000d+00 /
      data f( 2,11),f( 2,12) /  1.72900000d+01 ,  3.27840000d+01 /
      data f( 2,13),f( 2,14) /  3.27840000d+01 ,  1.72900000d+01 /
      data f( 2,15) /           4.99500000d+00 /
      data f( 3, 1),f( 3, 2) / -2.47770000d+01 , -3.07620000d+01 /
      data f( 3, 3),f( 3, 4) / -3.43170000d+01 , -3.54930000d+01 /
      data f( 3, 5),f( 3, 6) / -3.43170000d+01 , -3.07620000d+01 /
      data f( 3, 7),f( 3, 8) / -2.47770000d+01 , -1.63730000d+01 /
      data f( 3, 9),f( 3,10) / -5.92500000d+00 ,  5.30000000d+00 /
      data f( 3,11),f( 3,12) /  1.63570000d+01 ,  2.81240000d+01 /
      data f( 3,13),f( 3,14) /  2.81240000d+01 ,  1.63570000d+01 /
      data f( 3,15) /           5.30000000d+00 /
      data f( 4, 1),f( 4, 2) / -2.25910000d+01 , -2.80020000d+01 /
      data f( 4, 3),f( 4, 4) / -3.12120000d+01 , -3.22740000d+01 /
      data f( 4, 5),f( 4, 6) / -3.12120000d+01 , -2.80020000d+01 /
      data f( 4, 7),f( 4, 8) / -2.25910000d+01 , -1.49660000d+01 /
      data f( 4, 9),f( 4,10) / -5.37600000d+00 ,  5.22500000d+00 /
      data f( 4,11),f( 4,12) /  1.53420000d+01 ,  2.51380000d+01 /
      data f( 4,13),f( 4,14) /  2.51380000d+01 ,  1.53420000d+01 /
      data f( 4,15) /           5.22500000d+00 /
      data f( 5, 1),f( 5, 2) / -2.04580000d+01 , -2.53190000d+01 /
      data f( 5, 3),f( 5, 4) / -2.82020000d+01 , -2.91560000d+01 /
      data f( 5, 5),f( 5, 6) / -2.82020000d+01 , -2.53190000d+01 /
      data f( 5, 7),f( 5, 8) / -2.04580000d+01 , -1.35960000d+01 /
      data f( 5, 9),f( 5,10) / -4.90700000d+00 ,  4.86500000d+00 /
      data f( 5,11),f( 5,12) /  1.42930000d+01 ,  2.23060000d+01 /
      data f( 5,13),f( 5,14) /  2.23060000d+01 ,  1.42930000d+01 /
      data f( 5,15) /           4.86500000d+00 /
      data f( 6, 1),f( 6, 2) / -1.64810000d+01 , -2.03340000d+01 /
      data f( 6, 3),f( 6, 4) / -2.26200000d+01 , -2.33770000d+01 /
      data f( 6, 5),f( 6, 6) / -2.26200000d+01 , -2.03340000d+01 /
      data f( 6, 7),f( 6, 8) / -1.64810000d+01 , -1.10450000d+01 /
      data f( 6, 9),f( 6,10) / -4.13600000d+00 ,  3.74600000d+00 /
      data f( 6,11),f( 6,12) /  1.14370000d+01 ,  1.73410000d+01 /
      data f( 6,13),f( 6,14) /  1.73410000d+01 ,  1.14370000d+01 /
      data f( 6,15) /           3.74600000d+00 /
      data f( 7, 1),f( 7, 2) / -1.30370000d+01 , -1.60180000d+01 /
      data f( 7, 3),f( 7, 4) / -1.77910000d+01 , -1.83780000d+01 /
      data f( 7, 5),f( 7, 6) / -1.77910000d+01 , -1.60180000d+01 /
      data f( 7, 7),f( 7, 8) / -1.30370000d+01 , -8.84400000d+00 /
      data f( 7, 9),f( 7,10) / -3.52800000d+00 ,  2.54900000d+00 /
      data f( 7,11),f( 7,12) /  8.50300000d+00 ,  1.29680000d+01 /
      data f( 7,13),f( 7,14) /  1.29680000d+01 ,  8.50300000d+00 /
      data f( 7,15) /           2.54900000d+00 /
      data f( 8, 1),f( 8, 2) / -1.01940000d+01 , -1.24510000d+01 /
      data f( 8, 3),f( 8, 4) / -1.37960000d+01 , -1.42410000d+01 /
      data f( 8, 5),f( 8, 6) / -1.37960000d+01 , -1.24510000d+01 /
      data f( 8, 7),f( 8, 8) / -1.01940000d+01 , -7.03300000d+00 /
      data f( 8, 9),f( 8,10) / -3.04100000d+00 ,  1.51700000d+00 /
      data f( 8,11),f( 8,12) /  5.99600000d+00 ,  9.36100000d+00 /
      data f( 8,13),f( 8,14) /  9.36100000d+00 ,  5.99600000d+00 /
      data f( 8,15) /           1.51700000d+00 /
      data f( 9, 1),f( 9, 2) / -7.92600000d+00 , -9.60400000d+00 /
      data f( 9, 3),f( 9, 4) / -1.06060000d+01 , -1.09390000d+01 /
      data f( 9, 5),f( 9, 6) / -1.06060000d+01 , -9.60400000d+00 /
      data f( 9, 7),f( 9, 8) / -7.92600000d+00 , -5.58500000d+00 /
      data f( 9, 9),f( 9,10) / -2.63900000d+00 ,  7.17000000d-01 /
      data f( 9,11),f( 9,12) /  4.02100000d+00 ,  6.51800000d+00 /
      data f( 9,13),f( 9,14) /  6.51800000d+00 ,  4.02100000d+00 /
      data f( 9,15) /           7.17000000d-01 /
      data f(10, 1),f(10, 2) / -6.15500000d+00 , -7.38500000d+00 /
      data f(10, 3),f(10, 4) / -8.12300000d+00 , -8.36900000d+00 /
      data f(10, 5),f(10, 6) / -8.12300000d+00 , -7.38500000d+00 /
      data f(10, 7),f(10, 8) / -6.15500000d+00 , -4.44300000d+00 /
      data f(10, 9),f(10,10) / -2.29600000d+00 ,  1.43000000d-01 /
      data f(10,11),f(10,12) /  2.54400000d+00 ,  4.36500000d+00 /
      data f(10,13),f(10,14) /  4.36500000d+00 ,  2.54400000d+00 /
      data f(10,15) /           1.43000000d-01 /
      data f(11, 1),f(11, 2) / -4.78400000d+00 , -5.67800000d+00 /
      data f(11, 3),f(11, 4) / -6.21900000d+00 , -6.40000000d+00 /
      data f(11, 5),f(11, 6) / -6.21900000d+00 , -5.67800000d+00 /
      data f(11, 7),f(11, 8) / -4.78400000d+00 , -3.54500000d+00 /
      data f(11, 9),f(11,10) / -1.99500000d+00 , -2.40000000d-01 /
      data f(11,11),f(11,12) /  1.48400000d+00 ,  2.79100000d+00 /
      data f(11,13),f(11,14) /  2.79100000d+00 ,  1.48400000d+00 /
      data f(11,15) /          -2.40000000d-01 /
      data f(12, 1),f(12, 2) / -3.29700000d+00 , -3.84700000d+00 /
      data f(12, 3),f(12, 4) / -4.18400000d+00 , -4.29800000d+00 /
      data f(12, 5),f(12, 6) / -4.18400000d+00 , -3.84700000d+00 /
      data f(12, 7),f(12, 8) / -3.29700000d+00 , -2.54200000d+00 /
      data f(12, 9),f(12,10) / -1.60400000d+00 , -5.47000000d-01 /
      data f(12,11),f(12,12) /  4.83000000d-01 ,  1.25500000d+00 /
      data f(12,13),f(12,14) /  1.25500000d+00 ,  4.83000000d-01 /
      data f(12,15) /          -5.47000000d-01 /
      data f(13, 1),f(13, 2) / -1.81200000d+00 , -2.05600000d+00 /
      data f(13, 3),f(13, 4) / -2.20900000d+00 , -2.26200000d+00 /
      data f(13, 5),f(13, 6) / -2.20900000d+00 , -2.05600000d+00 /
      data f(13, 7),f(13, 8) / -1.81200000d+00 , -1.48700000d+00 /
      data f(13, 9),f(13,10) / -1.09100000d+00 , -6.52000000d-01 /
      data f(13,11),f(13,12) / -2.34000000d-01 ,  6.90000000d-02 /
      data f(13,13),f(13,14) /  6.90000000d-02 , -2.34000000d-01 /
      data f(13,15) /          -6.52000000d-01 /
      data f(14, 1),f(14, 2) / -1.03800000d+00 , -1.14800000d+00 /
      data f(14, 3),f(14, 4) / -1.21800000d+00 , -1.24200000d+00 /
      data f(14, 5),f(14, 6) / -1.21800000d+00 , -1.14800000d+00 /
      data f(14, 7),f(14, 8) / -1.03800000d+00 , -8.95000000d-01 /
      data f(14, 9),f(14,10) / -7.28000000d-01 , -5.51000000d-01 /
      data f(14,11),f(14,12) / -3.90000000d-01 , -2.80000000d-01 /
      data f(14,13),f(14,14) / -2.80000000d-01 , -3.90000000d-01 /
      data f(14,15) /          -5.51000000d-01 /
      data f(15, 1),f(15, 2) / -6.25000000d-01 , -6.77000000d-01 /
      data f(15, 3),f(15, 4) / -7.10000000d-01 , -7.22000000d-01 /
      data f(15, 5),f(15, 6) / -7.10000000d-01 , -6.77000000d-01 /
      data f(15, 7),f(15, 8) / -6.25000000d-01 , -5.60000000d-01 /
      data f(15, 9),f(15,10) / -4.87000000d-01 , -4.12000000d-01 /
      data f(15,11),f(15,12) / -3.50000000d-01 , -3.10000000d-01 /
      data f(15,13),f(15,14) / -3.10000000d-01 , -3.50000000d-01 /
      data f(15,15) /          -4.12000000d-01 /
      data f(16, 1),f(16, 2) / -3.96000000d-01 , -4.20000000d-01 /
      data f(16, 3),f(16, 4) / -4.37000000d-01 , -4.43000000d-01 /
      data f(16, 5),f(16, 6) / -4.37000000d-01 , -4.20000000d-01 /
      data f(16, 7),f(16, 8) / -3.96000000d-01 , -3.65000000d-01 /
      data f(16, 9),f(16,10) / -3.31000000d-01 , -2.96000000d-01 /
      data f(16,11),f(16,12) / -2.70000000d-01 , -2.54000000d-01 /
      data f(16,13),f(16,14) / -2.54000000d-01 , -2.70000000d-01 /
      data f(16,15) /          -2.96000000d-01 /
      data f(17, 1),f(17, 2) / -2.60000000d-01 , -2.73000000d-01 /
      data f(17, 3),f(17, 4) / -2.82000000d-01 , -2.85000000d-01 /
      data f(17, 5),f(17, 6) / -2.82000000d-01 , -2.73000000d-01 /
      data f(17, 7),f(17, 8) / -2.60000000d-01 , -2.45000000d-01 /
      data f(17, 9),f(17,10) / -2.29000000d-01 , -2.14000000d-01 /
      data f(17,11),f(17,12) / -2.02000000d-01 , -1.94000000d-01 /
      data f(17,13),f(17,14) / -1.94000000d-01 , -2.02000000d-01 /
      data f(17,15) /          -2.14000000d-01 /
      data f(18, 1),f(18, 2) / -1.26000000d-01 , -1.30000000d-01 /
      data f(18, 3),f(18, 4) / -1.33000000d-01 , -1.34000000d-01 /
      data f(18, 5),f(18, 6) / -1.33000000d-01 , -1.30000000d-01 /
      data f(18, 7),f(18, 8) / -1.26000000d-01 , -1.21000000d-01 /
      data f(18, 9),f(18,10) / -1.17000000d-01 , -1.13000000d-01 /
      data f(18,11),f(18,12) / -1.12000000d-01 , -1.11000000d-01 /
      data f(18,13),f(18,14) / -1.11000000d-01 , -1.12000000d-01 /
      data f(18,15) /          -1.13000000d-01 /
      data f(19, 1),f(19, 2) / -6.70000000d-02 , -7.00000000d-02 /
      data f(19, 3),f(19, 4) / -7.10000000d-02 , -7.20000000d-02 /
      data f(19, 5),f(19, 6) / -7.10000000d-02 , -7.00000000d-02 /
      data f(19, 7),f(19, 8) / -6.70000000d-02 , -6.50000000d-02 /
      data f(19, 9),f(19,10) / -6.30000000d-02 , -6.10000000d-02 /
      data f(19,11),f(19,12) / -6.10000000d-02 , -6.10000000d-02 /
      data f(19,13),f(19,14) / -6.10000000d-02 , -6.10000000d-02 /
      data f(19,15) /          -6.10000000d-02 /
      data f(20, 1),f(20, 2) / -3.70000000d-02 , -3.80000000d-02 /
      data f(20, 3),f(20, 4) / -4.00000000d-02 , -4.00000000d-02 /
      data f(20, 5),f(20, 6) / -4.00000000d-02 , -3.80000000d-02 /
      data f(20, 7),f(20, 8) / -3.70000000d-02 , -3.60000000d-02 /
      data f(20, 9),f(20,10) / -3.50000000d-02 , -3.40000000d-02 /
      data f(20,11),f(20,12) / -3.40000000d-02 , -3.40000000d-02 /
      data f(20,13),f(20,14) / -3.40000000d-02 , -3.40000000d-02 /
      data f(20,15) /          -3.40000000d-02 /
      data f(21, 1),f(21, 2) / -1.90000000d-02 , -2.10000000d-02 /
      data f(21, 3),f(21, 4) / -2.10000000d-02 , -2.10000000d-02 /
      data f(21, 5),f(21, 6) / -2.10000000d-02 , -2.10000000d-02 /
      data f(21, 7),f(21, 8) / -1.90000000d-02 , -1.90000000d-02 /
      data f(21, 9),f(21,10) / -1.90000000d-02 , -1.90000000d-02 /
      data f(21,11),f(21,12) / -1.90000000d-02 , -1.90000000d-02 /
      data f(21,13),f(21,14) / -1.90000000d-02 , -1.90000000d-02 /
      data f(21,15) /          -1.90000000d-02 /
      data f(22, 1),f(22, 2) / -9.00000000d-03 , -9.00000000d-03 /
      data f(22, 3),f(22, 4) / -1.00000000d-02 , -1.00000000d-02 /
      data f(22, 5),f(22, 6) / -1.00000000d-02 , -9.00000000d-03 /
      data f(22, 7),f(22, 8) / -9.00000000d-03 , -9.00000000d-03 /
      data f(22, 9),f(22,10) / -9.00000000d-03 , -9.00000000d-03 /
      data f(22,11),f(22,12) / -9.00000000d-03 , -9.00000000d-03 /
      data f(22,13),f(22,14) / -9.00000000d-03 , -9.00000000d-03 /
      data f(22,15) /          -9.00000000d-03 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 8.79698035d+00, 2.92842407d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 8.95048848d+00, 2.89715185d-02/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 7.17935222d+00, 2.87096851d-02/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 6.77499155d+00, 2.82697410d-02/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 7.17935222d+00, 2.87313509d-02/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 8.95048848d+00, 2.88848554d-02/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 8.79698035d+00, 2.96092277d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 3.92440207d+00, 1.75582339d-02/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-1.72795962d+01, 1.28378365d-02/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.43800953d+02,-2.06695801d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.94276473d+02, 1.49100484d-01/
      data fpp( 1,12,1),fpp( 1,12,2)/ 3.02014327d+02,-1.76852355d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.02014327d+02,-1.61573176d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.94276473d+02, 8.79837645d-02/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.43800953d+02, 2.08518118d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 4.00603931d+00, 2.71550764d-02/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 3.09902303d+00, 2.66798471d-02/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 1.44129556d+00, 2.61455351d-02/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 9.50016908d-01, 2.58180123d-02/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 1.44129556d+00, 2.61024156d-02/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 3.09902303d+00, 2.68523251d-02/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 4.00603931d+00, 2.65082839d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 1.45119586d+00, 2.31745392d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-1.40408075d+01, 2.35355942d-03/
      data fpp( 2,10,1),fpp( 2,10,2)/ 6.33980932d+01,-1.05087768d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.05347053d+02, 8.28215479d-02/
      data fpp( 2,12,1),fpp( 2,12,2)/ 2.21171346d+02,-1.28837415d-01/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.21171346d+02,-1.19718530d-01/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.05347053d+02, 4.63460085d-02/
      data fpp( 2,15,1),fpp( 2,15,2)/ 6.33980932d+01, 1.26274496d-01/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-8.21137582d-01, 2.48052456d-02/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-3.34658062d+00, 2.43095088d-02/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-4.54453446d+00, 2.37567190d-02/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-5.17505918d+00, 2.34036150d-02/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-4.54453446d+00, 2.37488210d-02/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-3.34658062d+00, 2.43411011d-02/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-8.21137582d-01, 2.46867745d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-7.29185517d-01, 2.20518008d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.05571737d+01, 9.74602213d-03/
      data fpp( 3,10,1),fpp( 3,10,2)/-6.73933263d+01,-1.44158894d-02/
      data fpp( 3,11,1),fpp( 3,11,2)/-4.26646863d+01, 3.78375353d-02/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.05300288d+02,-9.43342519d-02/
      data fpp( 3,13,1),fpp( 3,13,2)/ 2.05300288d+02,-8.89260118d-02/
      data fpp( 3,14,1),fpp( 3,14,2)/-4.26646863d+01, 1.62045748d-02/
      data fpp( 3,15,1),fpp( 3,15,2)/-6.73933263d+01, 6.67077126d-02/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-5.52148898d+00, 2.25325120d-02/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-7.71270057d+00, 2.20249760d-02/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-9.66315773d+00, 2.14275839d-02/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.02497802d+01, 2.11446882d-02/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-9.66315773d+00, 2.14336631d-02/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-7.71270057d+00, 2.20006594d-02/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-5.52148898d+00, 2.26236994d-02/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-3.93445380d+00, 2.03445429d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-7.93049753d+00, 1.38981290d-02/
      data fpp( 4,10,1),fpp( 4,10,2)/-2.18247878d+01,-1.52770590d-02/
      data fpp( 4,11,1),fpp( 4,11,2)/ 1.61116920d+01, 1.81701071d-02/
      data fpp( 4,12,1),fpp( 4,12,2)/-3.79724980d+01,-7.66633693d-02/
      data fpp( 4,13,1),fpp( 4,13,2)/-3.79724980d+01,-7.29749458d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.61116920d+01, 3.41641307d-03/
      data fpp( 4,15,1),fpp( 4,15,2)/-2.18247878d+01, 4.00492935d-02/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-8.89290649d+00, 2.02622675d-02/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.20026171d+01, 1.97954650d-02/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-1.38028346d+01, 1.92358723d-02/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-1.44258200d+01, 1.90010456d-02/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-1.38028346d+01, 1.92399453d-02/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-1.20026171d+01, 1.97791731d-02/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-8.89290649d+00, 2.03233624d-02/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-5.73299930d+00, 1.89873773d-02/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-5.72083614d+00, 1.33471283d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/-1.63075223d+01,-7.39589063d-03/
      data fpp( 5,11,1),fpp( 5,11,2)/-4.21820815d+01,-4.40356583d-03/
      data fpp( 5,12,1),fpp( 5,12,2)/ 3.89897042d+01,-5.98898461d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 3.89897042d+01,-5.85186789d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/-4.21820815d+01,-9.88823461d-03/
      data fpp( 5,15,1),fpp( 5,15,2)/-1.63075223d+01, 1.31716173d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-1.39105360d+01, 1.60430770d-02/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.72857984d+01, 1.56838461d-02/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.94599173d+01, 1.52415388d-02/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-2.01476498d+01, 1.50899987d-02/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.94599173d+01, 1.52384663d-02/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.72857984d+01, 1.56961363d-02/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.39105360d+01, 1.59969887d-02/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-9.18377520d+00, 1.52959088d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-3.92224282d+00, 1.11993759d-02/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.50392617d-02,-1.71341250d-03/
      data fpp( 6,11,1),fpp( 6,11,2)/ 4.79039857d+00,-1.58057259d-02/
      data fpp( 6,12,1),fpp( 6,12,2)/ 6.86713631d+00,-4.22836838d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/ 6.86713631d+00,-4.23660857d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/ 4.79039857d+00,-1.54761184d-02/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.50392617d-02,-2.94944081d-03/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-1.54149494d+01, 1.22981852d-02/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.92041893d+01, 1.20836295d-02/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-2.13074962d+01, 1.18472967d-02/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-2.19835807d+01, 1.16871836d-02/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-2.13074962d+01, 1.18439687d-02/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.92041893d+01, 1.20969415d-02/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-1.54149494d+01, 1.22482655d-02/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.00318999d+01, 1.16299967d-02/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-3.04019260d+00, 8.61174769d-03/
      data fpp( 7,10,1),fpp( 7,10,2)/ 4.66767932d+00,-4.16987461d-04/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.13204872d+01,-1.43237978d-02/
      data fpp( 7,12,1),fpp( 7,12,2)/ 2.23417505d+01,-3.16278212d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 2.23417505d+01,-3.19046376d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.13204872d+01,-1.32165321d-02/
      data fpp( 7,15,1),fpp( 7,15,2)/ 4.66767932d+00,-4.56923394d-03/
      data fpp( 8, 1,1),fpp( 8, 1,2)/-1.45796664d+01, 9.24124705d-03/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.82474445d+01, 9.11750591d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-2.04100980d+01, 9.00872932d-03/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-2.12180273d+01, 8.84757680d-03/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-2.04100980d+01, 9.00096348d-03/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.82474445d+01, 9.14856928d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.45796664d+01, 9.12475941d-03/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-9.18862522d+00, 8.59239308d-03/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-2.06698680d+00, 6.36566829d-03/
      data fpp( 8,10,1),fpp( 8,10,2)/ 6.09432199d+00,-9.50662167d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/ 1.39776525d+01,-1.07254034d-02/
      data fpp( 8,12,1),fpp( 8,12,2)/ 1.86658617d+01,-2.38433201d-02/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.86658617d+01,-2.40573380d-02/
      data fpp( 8,14,1),fpp( 8,14,2)/ 1.39776525d+01,-9.86933201d-03/
      data fpp( 8,15,1),fpp( 8,15,2)/ 6.09432199d+00,-3.30533399d-03/
      data fpp( 9, 1,1),fpp( 9, 1,2)/-1.25163851d+01, 6.82877865d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-1.58060327d+01, 6.76244269d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.78021119d+01, 6.68145057d-03/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-1.83943103d+01, 6.65175502d-03/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.78021119d+01, 6.67152935d-03/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-1.58060327d+01, 6.80212759d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-1.25163851d+01, 6.67996030d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-7.66359922d+00, 6.25803122d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-1.44186020d+00, 4.58791481d-03/
      data fpp( 9,10,1),fpp( 9,10,2)/ 5.75503273d+00,-9.69045164d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.25689026d+01,-7.66915300d-03/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.75948029d+01,-1.77336975d-02/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.75948029d+01,-1.78743309d-02/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.25689026d+01,-7.10661975d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 5.75503273d+00,-2.11919012d-03/
      data fpp(10, 1,1),fpp(10, 1,2)/-9.90479319d+00, 4.91996613d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.27284245d+01, 4.92006774d-03/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.44314546d+01, 4.91976292d-03/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.50047317d+01, 4.92088058d-03/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.44314546d+01, 4.91671475d-03/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.27284245d+01, 4.93226041d-03/
      data fpp(10, 7,1),fpp(10, 7,2)/-9.90479319d+00, 4.87424363d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/-6.05697791d+00, 4.49076509d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.01557239d+00, 3.26269601d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 4.78554708d+00,-2.15491312d-05/
      data fpp(10,11,1),fpp(10,11,2)/ 1.04467369d+01,-5.45649949d-03/
      data fpp(10,12,1),fpp(10,12,2)/ 1.44549268d+01,-1.29524529d-02/
      data fpp(10,13,1),fpp(10,13,2)/ 1.44549268d+01,-1.30443915d-02/
      data fpp(10,14,1),fpp(10,14,2)/ 1.04467369d+01,-5.08874529d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 4.78554708d+00,-1.40062735d-03/
      data fpp(11, 1,1),fpp(11, 1,2)/-7.86444214d+00, 3.46179858d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/-1.00802691d+01, 3.52640285d-03/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.13220699d+01, 3.61259003d-03/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.17367628d+01, 3.62323702d-03/
      data fpp(11, 5,1),fpp(11, 5,2)/-1.13220699d+01, 3.61446189d-03/
      data fpp(11, 6,1),fpp(11, 6,2)/-1.00802691d+01, 3.51891543d-03/
      data fpp(11, 7,1),fpp(11, 7,2)/-7.86444214d+00, 3.48987641d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/-4.70848914d+00, 3.22157895d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/-7.95850244d-01, 2.28380781d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 3.75277896d+00,-5.68101677d-05/
      data fpp(11,11,1),fpp(11,11,2)/ 8.19414969d+00,-3.91656713d-03/
      data fpp(11,12,1),fpp(11,12,2)/ 1.14354898d+01,-9.29692130d-03/
      data fpp(11,13,1),fpp(11,13,2)/ 1.14354898d+01,-9.36095255d-03/
      data fpp(11,14,1),fpp(11,14,2)/ 8.19414969d+00,-3.66044213d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 3.75277896d+00,-1.01727894d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/-5.14866408d+00, 2.03114221d-03/
      data fpp(12, 2,1),fpp(12, 2,2)/-6.54681982d+00, 2.12771558d-03/
      data fpp(12, 3,1),fpp(12, 3,2)/-7.37213061d+00, 2.23799546d-03/
      data fpp(12, 4,1),fpp(12, 4,2)/-7.64096952d+00, 2.30030257d-03/
      data fpp(12, 5,1),fpp(12, 5,2)/-7.37213061d+00, 2.24079427d-03/
      data fpp(12, 6,1),fpp(12, 6,2)/-6.54681982d+00, 2.11652035d-03/
      data fpp(12, 7,1),fpp(12, 7,2)/-5.14866408d+00, 2.07312434d-03/
      data fpp(12, 8,1),fpp(12, 8,2)/-3.20038426d+00, 1.89098228d-03/
      data fpp(12, 9,1),fpp(12, 9,2)/-7.03450927d-01, 1.34294655d-03/
      data fpp(12,10,1),fpp(12,10,2)/ 2.13370542d+00,-1.22768481d-04/
      data fpp(12,11,1),fpp(12,11,2)/ 4.98834309d+00,-2.47187263d-03/
      data fpp(12,12,1),fpp(12,12,2)/ 7.24508284d+00,-5.46974101d-03/
      data fpp(12,13,1),fpp(12,13,2)/ 7.24508284d+00,-5.51484065d-03/
      data fpp(12,14,1),fpp(12,14,2)/ 4.98834309d+00,-2.29147410d-03/
      data fpp(12,15,1),fpp(12,15,2)/ 2.13370542d+00,-7.99262949d-04/
      data fpp(13, 1,1),fpp(13, 1,2)/-2.64560965d+00, 8.19986046d-04/
      data fpp(13, 2,1),fpp(13, 2,2)/-3.25801509d+00, 9.10027908d-04/
      data fpp(13, 3,1),fpp(13, 3,2)/-3.61594011d+00, 9.99902322d-04/
      data fpp(13, 4,1),fpp(13, 4,2)/-3.72283984d+00, 1.09036281d-03/
      data fpp(13, 5,1),fpp(13, 5,2)/-3.61594011d+00, 9.98646458d-04/
      data fpp(13, 6,1),fpp(13, 6,2)/-3.25801509d+00, 9.15051365d-04/
      data fpp(13, 7,1),fpp(13, 7,2)/-2.64560965d+00, 8.01148084d-04/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.73367689d+00, 7.40356299d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/-5.99446886d-01, 4.97426719d-04/
      data fpp(13,10,1),fpp(13,10,2)/ 6.80475285d-01,-1.50063175d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 1.95281229d+00,-1.15717402d-03/
      data fpp(13,12,1),fpp(13,12,2)/ 2.93044104d+00,-2.12124075d-03/
      data fpp(13,13,1),fpp(13,13,2)/ 2.93044104d+00,-2.14769074d-03/
      data fpp(13,14,1),fpp(13,14,2)/ 1.95281229d+00,-1.05137408d-03/
      data fpp(13,15,1),fpp(13,15,2)/ 6.80475285d-01,-5.46812962d-04/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.33289731d+00, 3.41327640d-04/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.61311983d+00, 3.97344721d-04/
      data fpp(14, 3,1),fpp(14, 3,2)/-1.78010897d+00, 4.69293477d-04/
      data fpp(14, 4,1),fpp(14, 4,2)/-1.85167110d+00, 4.85481372d-04/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.78010897d+00, 4.68781037d-04/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.61311983d+00, 3.99394481d-04/
      data fpp(14, 7,1),fpp(14, 7,2)/-1.33289731d+00, 3.33641038d-04/
      data fpp(14, 8,1),fpp(14, 8,2)/-9.76908169d-01, 2.46041367d-04/
      data fpp(14, 9,1),fpp(14, 9,2)/-4.98761529d-01, 1.22193495d-04/
      data fpp(14,10,1),fpp(14,10,2)/ 8.83934432d-02,-1.34815346d-04/
      data fpp(14,11,1),fpp(14,11,2)/ 6.64407730d-01,-5.42932112d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 1.12115300d+00,-7.53456205d-04/
      data fpp(14,13,1),fpp(14,13,2)/ 1.12115300d+00,-7.68165328d-04/
      data fpp(14,14,1),fpp(14,14,2)/ 6.64407730d-01,-4.84095621d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 8.83934432d-02,-3.55452190d-04/
      data fpp(15, 1,1),fpp(15, 1,2)/-6.86801125d-01, 1.68856279d-04/
      data fpp(15, 2,1),fpp(15, 2,2)/-7.77505610d-01, 1.92287442d-04/
      data fpp(15, 3,1),fpp(15, 3,2)/-8.55624027d-01, 2.01993953d-04/
      data fpp(15, 4,1),fpp(15, 4,2)/-8.70475772d-01, 2.59736746d-04/
      data fpp(15, 5,1),fpp(15, 5,2)/-8.55624027d-01, 1.99059064d-04/
      data fpp(15, 6,1),fpp(15, 6,2)/-7.77505610d-01, 2.04026998d-04/
      data fpp(15, 7,1),fpp(15, 7,2)/-6.86801125d-01, 1.24832944d-04/
      data fpp(15, 8,1),fpp(15, 8,2)/-5.26690430d-01, 7.66412250d-05/
      data fpp(15, 9,1),fpp(15, 9,2)/-3.33506996d-01, 4.86021558d-05/
      data fpp(15,10,1),fpp(15,10,2)/-1.22049057d-01,-1.51049848d-04/
      data fpp(15,11,1),fpp(15,11,2)/ 9.35567870d-02,-2.24402763d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 2.40946950d-01,-2.71339101d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 2.40946950d-01,-2.73781315d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 9.35567870d-02,-2.14633910d-04/
      data fpp(15,15,1),fpp(15,15,2)/-1.22049057d-01,-1.87683045d-04/
      data fpp(16, 1,1),fpp(16, 1,2)/-3.35898195d-01, 3.11538836d-05/
      data fpp(16, 2,1),fpp(16, 2,2)/-4.12857736d-01, 6.76922328d-05/
      data fpp(16, 3,1),fpp(16, 3,2)/-4.37394924d-01, 1.18077185d-04/
      data fpp(16, 4,1),fpp(16, 4,2)/-4.50425818d-01, 1.19999026d-04/
      data fpp(16, 5,1),fpp(16, 5,2)/-4.37394924d-01, 1.21926710d-04/
      data fpp(16, 6,1),fpp(16, 6,2)/-4.12857736d-01, 5.22941340d-05/
      data fpp(16, 7,1),fpp(16, 7,2)/-3.35898195d-01, 8.88967539d-05/
      data fpp(16, 8,1),fpp(16, 8,2)/-2.76330110d-01, 1.21188502d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/-2.07210485d-01, 4.26278453d-05/
      data fpp(16,10,1),fpp(16,10,2)/-1.52197214d-01,-1.22630231d-04/
      data fpp(16,11,1),fpp(16,11,2)/-7.86348774d-02,-9.21069195d-05/
      data fpp(16,12,1),fpp(16,12,2)/-2.09408037d-02,-1.08942091d-04/
      data fpp(16,13,1),fpp(16,13,2)/-2.09408037d-02,-1.07120268d-04/
      data fpp(16,14,1),fpp(16,14,2)/-7.86348774d-02,-9.93942091d-05/
      data fpp(16,15,1),fpp(16,15,2)/-1.52197214d-01,-9.53028955d-05/
      data fpp(17, 1,1),fpp(17, 1,2)/-2.01606094d-01, 2.08858970d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/-2.11063448d-01, 3.82282061d-05/
      data fpp(17, 3,1),fpp(17, 3,2)/-2.26796278d-01, 6.62012787d-05/
      data fpp(17, 4,1),fpp(17, 4,2)/-2.31820958d-01, 5.69666793d-05/
      data fpp(17, 5,1),fpp(17, 5,2)/-2.26796278d-01, 6.59320042d-05/
      data fpp(17, 6,1),fpp(17, 6,2)/-2.11063448d-01, 3.93053039d-05/
      data fpp(17, 7,1),fpp(17, 7,2)/-2.01606094d-01, 1.68467801d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/-1.67989131d-01, 1.33075757d-05/
      data fpp(17, 9,1),fpp(17, 9,2)/-1.33651064d-01,-1.00770829d-05/
      data fpp(17,10,1),fpp(17,10,2)/-8.51620873d-02,-3.29992440d-05/
      data fpp(17,11,1),fpp(17,11,2)/-6.70172772d-02,-3.79259409d-05/
      data fpp(17,12,1),fpp(17,12,2)/-6.11837356d-02,-5.52969923d-05/
      data fpp(17,13,1),fpp(17,13,2)/-6.11837356d-02,-5.51460527d-05/
      data fpp(17,14,1),fpp(17,14,2)/-6.70172772d-02,-3.85296992d-05/
      data fpp(17,15,1),fpp(17,15,2)/-8.51620873d-02,-3.07351504d-05/
      data fpp(18, 1,1),fpp(18, 1,2)/-5.52326213d-02, 4.59145683d-07/
      data fpp(18, 2,1),fpp(18, 2,2)/-6.63807889d-02, 9.08170863d-06/
      data fpp(18, 3,1),fpp(18, 3,2)/-6.69137045d-02, 2.32140198d-05/
      data fpp(18, 4,1),fpp(18, 4,2)/-6.93242186d-02, 1.80622123d-05/
      data fpp(18, 5,1),fpp(18, 5,2)/-6.69137045d-02, 2.45371312d-05/
      data fpp(18, 6,1),fpp(18, 6,2)/-6.63807889d-02, 3.78926292d-06/
      data fpp(18, 7,1),fpp(18, 7,2)/-5.52326213d-02, 2.03058171d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/-5.38675534d-02,-2.50125314d-05/
      data fpp(18, 9,1),fpp(18, 9,2)/-4.74415657d-02, 1.97443084d-05/
      data fpp(18,10,1),fpp(18,10,2)/-4.64151313d-02,-5.39647024d-05/
      data fpp(18,11,1),fpp(18,11,2)/-3.56307296d-02, 1.61145011d-05/
      data fpp(18,12,1),fpp(18,12,2)/-2.79783913d-02,-1.04933021d-05/
      data fpp(18,13,1),fpp(18,13,2)/-2.79783913d-02,-6.57734427d-06/
      data fpp(18,14,1),fpp(18,14,2)/-3.56307296d-02, 4.50669790d-07/
      data fpp(18,15,1),fpp(18,15,2)/-4.64151313d-02, 4.77466510d-06/
      data fpp(19, 1,1),fpp(19, 1,2)/-2.74634211d-02, 3.77301008d-05/
      data fpp(19, 2,1),fpp(19, 2,2)/-2.14133968d-02, 2.45397983d-05/
      data fpp(19, 3,1),fpp(19, 3,2)/-2.75489039d-02,-1.58892942d-05/
      data fpp(19, 4,1),fpp(19, 4,2)/-2.48821682d-02, 3.90173784d-05/
      data fpp(19, 5,1),fpp(19, 5,2)/-2.75489039d-02,-2.01802193d-05/
      data fpp(19, 6,1),fpp(19, 6,2)/-2.14133968d-02, 4.17034987d-05/
      data fpp(19, 7,1),fpp(19, 7,2)/-2.74634211d-02,-2.66337757d-05/
      data fpp(19, 8,1),fpp(19, 8,2)/-2.45406558d-02, 4.83160406d-06/
      data fpp(19, 9,1),fpp(19, 9,2)/-2.45826731d-02, 7.30735945d-06/
      data fpp(19,10,1),fpp(19,10,2)/-2.31773875d-02,-3.40610419d-05/
      data fpp(19,11,1),fpp(19,11,2)/-2.44598044d-02, 8.93680801d-06/
      data fpp(19,12,1),fpp(19,12,2)/-2.49026992d-02,-1.68619019d-06/
      data fpp(19,13,1),fpp(19,13,2)/-2.49026992d-02, 5.90166567d-07/
      data fpp(19,14,1),fpp(19,14,2)/-2.44598044d-02,-1.68619019d-07/
      data fpp(19,15,1),fpp(19,15,2)/-2.31773875d-02, 8.43095096d-08/
      data fpp(20, 1,1),fpp(20, 1,2)/-8.91369432d-03,-3.73000394d-05/
      data fpp(20, 2,1),fpp(20, 2,2)/-1.59656239d-02,-1.53999212d-05/
      data fpp(20, 3,1),fpp(20, 3,2)/-8.89067971d-03, 3.88997242d-05/
      data fpp(20, 4,1),fpp(20, 4,2)/-1.11471086d-02,-2.01989756d-05/
      data fpp(20, 5,1),fpp(20, 5,2)/-8.89067971d-03, 4.18961783d-05/
      data fpp(20, 6,1),fpp(20, 6,2)/-1.59656239d-02,-2.73857376d-05/
      data fpp(20, 7,1),fpp(20, 7,2)/-8.91369432d-03, 7.64677203d-06/
      data fpp(20, 8,1),fpp(20, 8,2)/-9.96982344d-03,-3.20135054d-06/
      data fpp(20, 9,1),fpp(20, 9,2)/-1.02277419d-02, 5.15863011d-06/
      data fpp(20,10,1),fpp(20,10,2)/-1.08753188d-02,-1.74331699d-05/
      data fpp(20,11,1),fpp(20,11,2)/-1.05300527d-02, 4.57404953d-06/
      data fpp(20,12,1),fpp(20,12,2)/-1.04108117d-02,-8.63028213d-07/
      data fpp(20,13,1),fpp(20,13,2)/-1.04108117d-02, 3.02059875d-07/
      data fpp(20,14,1),fpp(20,14,2)/-1.05300527d-02,-8.63028213d-08/
      data fpp(20,15,1),fpp(20,15,2)/-1.08753188d-02, 4.31514107d-08/
      data fpp(21, 1,1),fpp(21, 1,2)/-8.88180162d-03, 3.90496277d-05/
      data fpp(21, 2,1),fpp(21, 2,2)/-4.72410745d-03, 2.19007446d-05/
      data fpp(21, 3,1),fpp(21, 3,2)/-8.88837723d-03,-6.65260601d-06/
      data fpp(21, 4,1),fpp(21, 4,2)/-8.52939756d-03, 4.70967945d-06/
      data fpp(21, 5,1),fpp(21, 5,2)/-8.88837723d-03,-1.21861118d-05/
      data fpp(21, 6,1),fpp(21, 6,2)/-4.72410745d-03, 4.40347678d-05/
      data fpp(21, 7,1),fpp(21, 7,2)/-8.88180162d-03,-4.39529593d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/-7.58005045d-03, 1.17770694d-05/
      data fpp(21, 9,1),fpp(21, 9,2)/-6.50635947d-03,-3.15531846d-06/
      data fpp(21,10,1),fpp(21,10,2)/-5.32133750d-03, 8.44204410d-07/
      data fpp(21,11,1),fpp(21,11,2)/-5.41998496d-03,-2.21499177d-07/
      data fpp(21,12,1),fpp(21,12,2)/-5.45405379d-03, 4.17922975d-08/
      data fpp(21,13,1),fpp(21,13,2)/-5.45405379d-03,-1.46273041d-08/
      data fpp(21,14,1),fpp(21,14,2)/-5.41998496d-03, 4.17922975d-09/
      data fpp(21,15,1),fpp(21,15,2)/-5.32133750d-03,-2.08961488d-09/
      data fpp(22, 1,1),fpp(22, 1,2)/-3.55909919d-03,-2.84189600d-05/
      data fpp(22, 2,1),fpp(22, 2,2)/ 4.86205372d-03,-1.31620800d-05/
      data fpp(22, 3,1),fpp(22, 3,2)/-3.55581139d-03, 2.10672800d-05/
      data fpp(22, 4,1),fpp(22, 4,2)/-2.73530122d-03,-1.11070400d-05/
      data fpp(22, 5,1),fpp(22, 5,2)/-3.55581139d-03, 2.33608799d-05/
      data fpp(22, 6,1),fpp(22, 6,2)/ 4.86205372d-03,-2.23364796d-05/
      data fpp(22, 7,1),fpp(22, 7,2)/-3.55909919d-03, 5.98503836d-06/
      data fpp(22, 8,1),fpp(22, 8,2)/-1.70997478d-03,-1.60367387d-06/
      data fpp(22, 9,1),fpp(22, 9,2)/ 2.53179735d-04, 4.29657123d-07/
      data fpp(22,10,1),fpp(22,10,2)/ 2.16066875d-03,-1.14954621d-07/
      data fpp(22,11,1),fpp(22,11,2)/ 2.20999248d-03, 3.01613609d-08/
      data fpp(22,12,1),fpp(22,12,2)/ 2.22702689d-03,-5.69082282d-09/
      data fpp(22,13,1),fpp(22,13,2)/ 2.22702689d-03, 1.99178799d-09/
      data fpp(22,14,1),fpp(22,14,2)/ 2.20999248d-03,-5.69082282d-10/
      data fpp(22,15,1),fpp(22,15,2)/ 2.16066875d-03, 2.84541141d-10/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-5.00855998d-02,-2.33965321d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.81930643d-02,-7.36918976d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/ 4.98069757d-02,-1.09852178d-01/
      data fpppp( 1, 7),fpppp( 1, 8)/ 2.74123071d-01,-1.26978431d+00/
      data fpppp( 1, 9),fpppp( 1,10)/ 3.82512898d+00,-3.09365874d+00/
      data fpppp( 1,11),fpppp( 1,12)/ 1.91320420d+00,-1.12341805d+00/
      data fpppp( 1,13),fpppp( 1,14)/-8.18483556d-01, 6.93466227d-01/
      data fpppp( 1,15) /             1.48035866d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/-2.55815536d-02,-9.70142913d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.93445983d-02, 2.30996542d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/ 3.03689782d-02,-5.37989487d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.39784145d-01,-7.13049214d-01/
      data fpppp( 2, 9),fpppp( 2,10)/ 1.93618312d+00,-1.45582901d+00/
      data fpppp( 2,11),fpppp( 2,12)/ 1.75773647d+00,-1.14259690d+00/
      data fpppp( 2,13),fpppp( 2,14)/-9.25806340d-01, 8.90574247d-01/
      data fpppp( 2,15) /             1.79602934d+00 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 2.00118591d-02, 1.50021590d-02/
      data fpppp( 3, 3),fpppp( 3, 4)/-3.71143456d-04, 2.05281620d-02/
      data fpppp( 3, 5),fpppp( 3, 6)/-6.07853812d-03, 3.78317376d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-6.55990609d-02, 7.85550478d-02/
      data fpppp( 3, 9),fpppp( 3,10)/-8.43817548d-01, 4.76225280d-01/
      data fpppp( 3,11),fpppp( 3,12)/ 3.83280399d+00,-2.41326117d+00/
      data fpppp( 3,13),fpppp( 3,14)/-2.11556770d+00, 2.64203010d+00/
      data fpppp( 3,15) /             4.94162734d+00 /
      data fpppp( 4, 1),fpppp( 4, 2)/-8.20751092d-03, 1.17604918d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 1.79485798d-02, 8.85971365d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.70072614d-02, 4.94132282d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/-2.23272871d-02, 4.81172411d-02/
      data fpppp( 4, 9),fpppp( 4,10)/-5.05126413d-01, 1.37849361d+00/
      data fpppp( 4,11),fpppp( 4,12)/-1.89900184d+00, 6.96273553d-01/
      data fpppp( 4,13),fpppp( 4,14)/ 4.83205960d-01,-1.04673146d+00/
      data fpppp( 4,15) /            -1.81752029d+00 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 1.43116285d-02, 1.33067563d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/ 1.10309323d-02, 1.32034398d-02/
      data fpppp( 5, 5),fpppp( 5, 6)/ 1.09135590d-02, 1.37762492d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/ 1.25510304d-02,-6.09685755d-02/
      data fpppp( 5, 9),fpppp( 5,10)/ 4.24586296d-02,-7.44796900d-01/
      data fpppp( 5,11),fpppp( 5,12)/ 2.01945659d+00,-9.10248741d-01/
      data fpppp( 5,13),fpppp( 5,14)/-7.14135641d-01, 1.23500419d+00/
      data fpppp( 5,15) /             2.19689960d+00 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 9.36364789d-03, 1.16021504d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.62963592d-02, 1.23955962d-02/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.66491568d-02, 1.01909604d-02/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.46556107d-02, 1.22765036d-02/
      data fpppp( 6, 9),fpppp( 6,10)/-3.16753316d-02, 3.31650928d-02/
      data fpppp( 6,11),fpppp( 6,12)/-4.70909830d-02,-8.52316663d-03/
      data fpppp( 6,13),fpppp( 6,14)/-1.31871408d-02,-2.84350861d-02/
      data fpppp( 6,15) /            -3.67945204d-02 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 1.93785968d-02, 1.69950082d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/ 1.37973487d-02, 1.34489375d-02/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.35370482d-02, 1.80362104d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.54740885d-02, 1.56960136d-02/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.82613246d-02,-4.57714350d-02/
      data fpppp( 7,11),fpppp( 7,12)/ 1.01520575d-01,-9.82035409d-02/
      data fpppp( 7,13),fpppp( 7,14)/-8.67875632d-02, 5.58566639d-02/
      data fpppp( 7,15) /             1.25468231d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.63149943d-02, 1.55317581d-02/
      data fpppp( 8, 3),fpppp( 8, 4)/ 1.18654509d-02, 1.82898920d-02/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.19264926d-02, 1.52875916d-02/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.72306187d-02, 1.91857151d-02/
      data fpppp( 8, 9),fpppp( 8,10)/ 9.86235713d-03, 3.74507828d-03/
      data fpppp( 8,11),fpppp( 8,12)/-4.15213646d-02,-2.93669058d-02/
      data fpppp( 8,13),fpppp( 8,14)/-3.17848737d-02,-3.18494930d-02/
      data fpppp( 8,15) /            -3.25244401d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 1.20318967d-02, 1.25370178d-02/
      data fpppp( 9, 3),fpppp( 9, 4)/ 1.54341430d-02, 9.95925383d-03/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.57926483d-02, 1.11029966d-02/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.74094762d-02, 1.30473939d-02/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.25381354d-02,-4.69069994d-03/
      data fpppp( 9,11),fpppp( 9,12)/-1.67567177d-02,-3.55606091d-02/
      data fpppp( 9,13),fpppp( 9,14)/-3.57168212d-02,-1.61318692d-02/
      data fpppp( 9,15) /            -7.03388155d-03 /
      data fpppp(10, 1),fpppp(10, 2)/ 1.11057996d-02, 1.12234092d-02/
      data fpppp(10, 3),fpppp(10, 4)/ 1.12366420d-02, 1.16151937d-02/
      data fpppp(10, 5),fpppp(10, 6)/ 1.10958449d-02, 1.17865976d-02/
      data fpppp(10, 7),fpppp(10, 8)/ 8.99384292d-03, 1.36890670d-02/
      data fpppp(10, 9),fpppp(10,10)/ 7.86530368d-03, 4.32554797d-04/
      data fpppp(10,11),fpppp(10,12)/-1.79913003d-02,-2.76473493d-02/
      data fpppp(10,13),fpppp(10,14)/-2.83079995d-02,-1.53486992d-02/
      data fpppp(10,15) /            -9.47719913d-03 /
      data fpppp(11, 1),fpppp(11, 2)/ 1.11416500d-02, 9.87585548d-03/
      data fpppp(11, 3),fpppp(11, 4)/ 7.79650339d-03, 8.56460299d-03/
      data fpppp(11, 5),fpppp(11, 6)/ 7.70823138d-03, 1.02289435d-02/
      data fpppp(11, 7),fpppp(11, 8)/ 9.81756981d-03, 6.90833637d-03/
      data fpppp(11, 9),fpppp(11,10)/ 7.95023871d-03,-5.49872808d-04/
      data fpppp(11,11),fpppp(11,12)/-1.21862559d-02,-2.27069420d-02/
      data fpppp(11,13),fpppp(11,14)/-2.30262488d-02,-1.09090287d-02/
      data fpppp(11,15) /            -5.33947453d-03 /
      data fpppp(12, 1),fpppp(12, 2)/ 5.89696483d-03, 5.71888014d-03/
      data fpppp(12, 3),fpppp(12, 4)/ 5.59821137d-03, 5.27658690d-03/
      data fpppp(12, 5),fpppp(12, 6)/ 5.55611081d-03, 5.88728237d-03/
      data fpppp(12, 7),fpppp(12, 8)/ 5.26545647d-03, 6.05833708d-03/
      data fpppp(12, 9),fpppp(12,10)/ 3.42040542d-03, 6.73422249d-04/
      data fpppp(12,11),fpppp(12,12)/-5.06521484d-03,-1.62864382d-02/
      data fpppp(12,13),fpppp(12,14)/-1.63102705d-02,-4.96988582d-03/
      data fpppp(12,15) /             3.15938413d-04 /
      data fpppp(13, 1),fpppp(13, 2)/ 2.60322934d-03, 2.49705670d-03/
      data fpppp(13, 3),fpppp(13, 4)/ 2.67736900d-03, 1.85498400d-03/
      data fpppp(13, 5),fpppp(13, 6)/ 2.73066368d-03, 2.28387795d-03/
      data fpppp(13, 7),fpppp(13, 8)/ 3.40264964d-03, 2.07716291d-03/
      data fpppp(13, 9),fpppp(13,10)/ 1.62653362d-03, 1.58232375d-04/
      data fpppp(13,11),fpppp(13,12)/-2.71457271d-03,-6.98243749d-03/
      data fpppp(13,13),fpppp(13,14)/-7.02426351d-03,-2.54726863d-03/
      data fpppp(13,15) /            -4.69157955d-04 /
      data fpppp(14, 1),fpppp(14, 2)/ 1.27006397d-03, 1.21300075d-03/
      data fpppp(14, 3),fpppp(14, 4)/ 6.71935754d-04, 1.82487696d-03/
      data fpppp(14, 5),fpppp(14, 6)/ 6.16011898d-04, 1.43669617d-03/
      data fpppp(14, 7),fpppp(14, 8)/ 4.31206124d-04, 1.38447634d-03/
      data fpppp(14, 9),fpppp(14,10)/ 1.36033871d-03,-2.85331187d-04/
      data fpppp(14,11),fpppp(14,12)/-8.87455139d-04,-3.32098905d-03/
      data fpppp(14,13),fpppp(14,14)/-3.29566348d-03,-9.88757393d-04/
      data fpppp(14,15) /             9.45522655d-05 /
      data fpppp(15, 1),fpppp(15, 2)/-3.34996108d-04, 3.39621148d-05/
      data fpppp(15, 3),fpppp(15, 4)/ 9.54311690d-04,-5.52084665d-05/
      data fpppp(15, 5),fpppp(15, 6)/ 1.04873148d-03,-3.43717041d-04/
      data fpppp(15, 7),fpppp(15, 8)/ 1.08130072d-03, 1.82886710d-04/
      data fpppp(15, 9),fpppp(15,10)/ 1.71516801d-04, 2.27516399d-04/
      data fpppp(15,11),fpppp(15,12)/-8.32708081d-04,-9.89624938d-04/
      data fpppp(15,13),fpppp(15,14)/-1.03647604d-03,-6.45303661d-04/
      data fpppp(15,15) /            -4.75250177d-04 /
      data fpppp(16, 1),fpppp(16, 2)/ 9.04375630d-04, 5.82240468d-04/
      data fpppp(16, 3),fpppp(16, 4)/-8.79963705d-05, 4.60122656d-04/
      data fpppp(16, 5),fpppp(16, 6)/-1.88786969d-04, 9.85402860d-04/
      data fpppp(16, 7),fpppp(16, 8)/-6.07483342d-04, 4.01043232d-04/
      data fpppp(16, 9),fpppp(16,10)/-4.23597231d-04, 4.46964462d-04/
      data fpppp(16,11),fpppp(16,12)/-2.51316696d-04,-3.93793434d-04/
      data fpppp(16,13),fpppp(16,14)/-4.23783562d-04,-1.31356187d-04/
      data fpppp(16,15) /            -2.88744668d-06 /
      data fpppp(17, 1),fpppp(17, 2)/-2.24080601d-04,-7.97756047d-05/
      data fpppp(17, 3),fpppp(17, 4)/ 1.66654461d-04, 5.56467849d-05/
      data fpppp(17, 5),fpppp(17, 6)/ 2.13719958d-04,-2.68037589d-04/
      data fpppp(17, 7),fpppp(17, 8)/ 4.81901840d-04,-2.09993231d-04/
      data fpppp(17, 9),fpppp(17,10)/ 4.01337297d-04,-5.46301353d-04/
      data fpppp(17,11),fpppp(17,12)/-3.67818828d-05,-4.52472202d-05/
      data fpppp(17,13),fpppp(17,14)/-2.08736461d-05,-1.34276179d-04/
      data fpppp(17,15) /            -1.80697742d-04 /
      data fpppp(18, 1),fpppp(18, 2)/ 2.20709322d-04, 1.26895916d-04/
      data fpppp(18, 3),fpppp(18, 4)/-9.13778733d-05, 1.25959678d-04/
      data fpppp(18, 5),fpppp(18, 6)/-1.23199157d-04, 2.54181052d-04/
      data fpppp(18, 7),fpppp(18, 8)/-2.56609940d-04, 1.85272729d-04/
      data fpppp(18, 9),fpppp(18,10)/-1.80825788d-04, 2.14057230d-04/
      data fpppp(18,11),fpppp(18,12)/-8.99250936d-05,-4.22806617d-05/
      data fpppp(18,13),fpppp(18,14)/-5.77656170d-05,-2.79852724d-05/
      data fpppp(18,15) /            -1.82170996d-05 /
      data fpppp(19, 1),fpppp(19, 2)/-3.12570447d-04,-1.60180535d-04/
      data fpppp(19, 3),fpppp(19, 4)/ 2.22160701d-04,-2.00327698d-04/
      data fpppp(19, 5),fpppp(19, 6)/ 2.59141807d-04,-3.08104958d-04/
      data fpppp(19, 7),fpppp(19, 8)/ 2.42146140d-04,-1.22112228d-04/
      data fpppp(19, 9),fpppp(19,10)/ 6.84158145d-05,-6.47128539d-05/
      data fpppp(19,11),fpppp(19,12)/ 2.91734465d-05,-1.60960259d-06/
      data fpppp(19,13),fpppp(19,14)/ 3.52892858d-06, 8.61932182d-06/
      data fpppp(19,15) /             1.23651136d-05 /
      data fpppp(20, 1),fpppp(20, 2)/ 3.55765892d-04, 1.81439369d-04/
      data fpppp(20, 3),fpppp(20, 4)/-2.33910937d-04, 1.94321995d-04/
      data fpppp(20, 5),fpppp(20, 6)/-2.72605582d-04, 3.36217948d-04/
      data fpppp(20, 7),fpppp(20, 8)/-2.24653779d-04, 7.59136432d-05/
      data fpppp(20, 9),fpppp(20,10)/-3.11081515d-05, 2.51394541d-05/
      data fpppp(20,11),fpppp(20,12)/-9.87908458d-06, 8.15372495d-07/
      data fpppp(20,13),fpppp(20,14)/-1.08380244d-06,-2.28238485d-06/
      data fpppp(20,15) /            -3.34816995d-06 /
      data fpppp(21, 1),fpppp(21, 2)/-2.02214279d-04,-1.02134628d-04/
      data fpppp(21, 3),fpppp(21, 4)/ 1.11434956d-04,-7.22102277d-05/
      data fpppp(21, 5),fpppp(21, 6)/ 1.34328394d-04,-1.93708383d-04/
      data fpppp(21, 7),fpppp(21, 8)/ 1.41187300d-04,-4.34740972d-05/
      data fpppp(21, 9),fpppp(21,10)/ 1.90254762d-05,-2.59479481d-05/
      data fpppp(21,11),fpppp(21,12)/ 7.74615049d-06,-1.16193619d-06/
      data fpppp(21,13),fpppp(21,14)/ 6.34798257d-07, 5.59212694d-07/
      data fpppp(21,15) /             1.00306862d-06 /
      data fpppp(22, 1),fpppp(22, 2)/-4.09818305d-04,-2.07081797d-04/
      data fpppp(22, 3),fpppp(22, 4)/ 2.27804411d-04,-1.49833331d-04/
      data fpppp(22, 5),fpppp(22, 6)/ 2.73067695d-04,-3.88134932d-04/
      data fpppp(22, 7),fpppp(22, 8)/ 2.69130951d-04,-7.21722342d-05/
      data fpppp(22, 9),fpppp(22,10)/ 2.63997915d-05,-3.67668617d-05/
      data fpppp(22,11),fpppp(22,12)/ 9.17773809d-06,-1.88144952d-06/
      data fpppp(22,13),fpppp(22,14)/ 5.44447035d-07,-5.25848108d-07/
      data fpppp(22,15) /            -3.78413428d-07 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /

      iprint=0

      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
      subroutine ch3po_1ap_d(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     2a' surface
c     difference of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  1.09800000d+00 ,  1.49500000d+00 /
      data f( 1, 3),f( 1, 4) /  9.66000000d-01 ,  0.00000000d+00 /
      data f( 1, 5),f( 1, 6) / -9.66000000d-01 , -1.49500000d+00 /
      data f( 1, 7),f( 1, 8) / -1.09800000d+00 ,  5.59000000d-01 /
      data f( 1, 9),f( 1,10) /  2.78900000d+00 ,  5.14100000d+00 /
      data f( 1,11),f( 1,12) /  1.10870000d+01 ,  3.48490000d+01 /
      data f( 1,13),f( 1,14) /  3.48490000d+01 ,  1.10870000d+01 /
      data f( 1,15) /           5.14100000d+00 /
      data f( 2, 1),f( 2, 2) /  6.49000000d-01 ,  1.14400000d+00 /
      data f( 2, 3),f( 2, 4) /  7.70000000d-01 ,  0.00000000d+00 /
      data f( 2, 5),f( 2, 6) / -7.70000000d-01 , -1.14400000d+00 /
      data f( 2, 7),f( 2, 8) / -6.49000000d-01 ,  1.16200000d+00 /
      data f( 2, 9),f( 2,10) /  4.12000000d+00 ,  6.52100000d+00 /
      data f( 2,11),f( 2,12) /  1.15040000d+01 ,  2.74960000d+01 /
      data f( 2,13),f( 2,14) /  2.74960000d+01 ,  1.15040000d+01 /
      data f( 2,15) /           6.52100000d+00 /
      data f( 3, 1),f( 3, 2) /  3.30000000d-01 ,  8.69000000d-01 /
      data f( 3, 3),f( 3, 4) /  6.11000000d-01 ,  0.00000000d+00 /
      data f( 3, 5),f( 3, 6) / -6.11000000d-01 , -8.69000000d-01 /
      data f( 3, 7),f( 3, 8) / -3.30000000d-01 ,  1.51300000d+00 /
      data f( 3, 9),f( 3,10) /  4.83000000d+00 ,  8.63500000d+00 /
      data f( 3,11),f( 3,12) /  1.33960000d+01 ,  2.40700000d+01 /
      data f( 3,13),f( 3,14) /  2.40700000d+01 ,  1.33960000d+01 /
      data f( 3,15) /           8.63500000d+00 /
      data f( 4, 1),f( 4, 2) /  1.05000000d-01 ,  6.54000000d-01 /
      data f( 4, 3),f( 4, 4) /  4.83000000d-01 ,  0.00000000d+00 /
      data f( 4, 5),f( 4, 6) / -4.83000000d-01 , -6.54000000d-01 /
      data f( 4, 7),f( 4, 8) / -1.05000000d-01 ,  1.68500000d+00 /
      data f( 4, 9),f( 4,10) /  5.09600000d+00 ,  9.69200000d+00 /
      data f( 4,11),f( 4,12) /  1.47080000d+01 ,  2.33060000d+01 /
      data f( 4,13),f( 4,14) /  2.33060000d+01 ,  1.47080000d+01 /
      data f( 4,15) /           9.69200000d+00 /
      data f( 5, 1),f( 5, 2) / -4.80000000d-02 ,  4.87000000d-01 /
      data f( 5, 3),f( 5, 4) /  3.80000000d-01 ,  0.00000000d+00 /
      data f( 5, 5),f( 5, 6) / -3.80000000d-01 , -4.87000000d-01 /
      data f( 5, 7),f( 5, 8) /  4.80000000d-02 ,  1.73700000d+00 /
      data f( 5, 9),f( 5,10) /  5.06500000d+00 ,  9.94300000d+00 /
      data f( 5,11),f( 5,12) /  1.55570000d+01 ,  2.22070000d+01 /
      data f( 5,13),f( 5,14) /  2.22070000d+01 ,  1.55570000d+01 /
      data f( 5,15) /           9.94300000d+00 /
      data f( 6, 1),f( 6, 2) / -2.10000000d-01 ,  2.59000000d-01 /
      data f( 6, 3),f( 6, 4) /  2.31000000d-01 ,  0.00000000d+00 /
      data f( 6, 5),f( 6, 6) / -2.31000000d-01 , -2.59000000d-01 /
      data f( 6, 7),f( 6, 8) /  2.10000000d-01 ,  1.63000000d+00 /
      data f( 6, 9),f( 6,10) /  4.52700000d+00 ,  9.10000000d+00 /
      data f( 6,11),f( 6,12) /  1.46330000d+01 ,  1.96590000d+01 /
      data f( 6,13),f( 6,14) /  1.96590000d+01 ,  1.46330000d+01 /
      data f( 6,15) /           9.10000000d+00 /
      data f( 7, 1),f( 7, 2) / -2.60000000d-01 ,  1.27000000d-01 /
      data f( 7, 3),f( 7, 4) /  1.38000000d-01 ,  0.00000000d+00 /
      data f( 7, 5),f( 7, 6) / -1.38000000d-01 , -1.27000000d-01 /
      data f( 7, 7),f( 7, 8) /  2.60000000d-01 ,  1.40000000d+00 /
      data f( 7, 9),f( 7,10) /  3.75600000d+00 ,  7.57000000d+00 /
      data f( 7,11),f( 7,12) /  1.22720000d+01 ,  1.63110000d+01 /
      data f( 7,13),f( 7,14) /  1.63110000d+01 ,  1.22720000d+01 /
      data f( 7,15) /           7.57000000d+00 /
      data f( 8, 1),f( 8, 2) / -2.55000000d-01 ,  5.40000000d-02 /
      data f( 8, 3),f( 8, 4) /  7.80000000d-02 ,  0.00000000d+00 /
      data f( 8, 5),f( 8, 6) / -7.80000000d-02 , -5.40000000d-02 /
      data f( 8, 7),f( 8, 8) /  2.55000000d-01 ,  1.14500000d+00 /
      data f( 8, 9),f( 8,10) /  2.98600000d+00 ,  5.98500000d+00 /
      data f( 8,11),f( 8,12) /  9.72000000d+00 ,  1.29370000d+01 /
      data f( 8,13),f( 8,14) /  1.29370000d+01 ,  9.72000000d+00 /
      data f( 8,15) /           5.98500000d+00 /
      data f( 9, 1),f( 9, 2) / -2.27000000d-01 ,  1.20000000d-02 /
      data f( 9, 3),f( 9, 4) /  4.20000000d-02 ,  0.00000000d+00 /
      data f( 9, 5),f( 9, 6) / -4.20000000d-02 , -1.20000000d-02 /
      data f( 9, 7),f( 9, 8) /  2.27000000d-01 ,  9.06000000d-01 /
      data f( 9, 9),f( 9,10) /  2.30700000d+00 ,  4.58100000d+00 /
      data f( 9,11),f( 9,12) /  7.42700000d+00 ,  9.90000000d+00 /
      data f( 9,13),f( 9,14) /  9.90000000d+00 ,  7.42700000d+00 /
      data f( 9,15) /           4.58100000d+00 /
      data f(10, 1),f(10, 2) / -1.91000000d-01 , -1.00000000d-02 /
      data f(10, 3),f(10, 4) /  1.80000000d-02 ,  0.00000000d+00 /
      data f(10, 5),f(10, 6) / -1.80000000d-02 ,  1.00000000d-02 /
      data f(10, 7),f(10, 8) /  1.91000000d-01 ,  7.02000000d-01 /
      data f(10, 9),f(10,10) /  1.74700000d+00 ,  3.43200000d+00 /
      data f(10,11),f(10,12) /  5.53600000d+00 ,  7.37100000d+00 /
      data f(10,13),f(10,14) /  7.37100000d+00 ,  5.53600000d+00 /
      data f(10,15) /           3.43200000d+00 /
      data f(11, 1),f(11, 2) / -1.59000000d-01 , -2.30000000d-02 /
      data f(11, 3),f(11, 4) /  4.00000000d-03 ,  0.00000000d+00 /
      data f(11, 5),f(11, 6) / -4.00000000d-03 ,  2.30000000d-02 /
      data f(11, 7),f(11, 8) /  1.59000000d-01 ,  5.38000000d-01 /
      data f(11, 9),f(11,10) /  1.30400000d+00 ,  2.52900000d+00 /
      data f(11,11),f(11,12) /  4.04900000d+00 ,  5.37300000d+00 /
      data f(11,13),f(11,14) /  5.37300000d+00 ,  4.04900000d+00 /
      data f(11,15) /           2.52900000d+00 /
      data f(12, 1),f(12, 2) / -1.24000000d-01 , -3.40000000d-02 /
      data f(12, 3),f(12, 4) / -8.00000000d-03 ,  0.00000000d+00 /
      data f(12, 5),f(12, 6) /  8.00000000d-03 ,  3.40000000d-02 /
      data f(12, 7),f(12, 8) /  1.24000000d-01 ,  3.57000000d-01 /
      data f(12, 9),f(12,10) /  8.24000000d-01 ,  1.56200000d+00 /
      data f(12,11),f(12,12) /  2.46500000d+00 ,  3.24300000d+00 /
      data f(12,13),f(12,14) /  3.24300000d+00 ,  2.46500000d+00 /
      data f(12,15) /           1.56200000d+00 /
      data f(13, 1),f(13, 2) / -8.50000000d-02 , -3.80000000d-02 /
      data f(13, 3),f(13, 4) / -1.40000000d-02 ,  0.00000000d+00 /
      data f(13, 5),f(13, 6) /  1.40000000d-02 ,  3.80000000d-02 /
      data f(13, 7),f(13, 8) /  8.50000000d-02 ,  1.84000000d-01 /
      data f(13, 9),f(13,10) /  3.72000000d-01 ,  6.65000000d-01 /
      data f(13,11),f(13,12) /  1.01300000d+00 ,  1.30000000d+00 /
      data f(13,13),f(13,14) /  1.30000000d+00 ,  1.01300000d+00 /
      data f(13,15) /           6.65000000d-01 /
      data f(14, 1),f(14, 2) / -5.80000000d-02 , -3.00000000d-02 /
      data f(14, 3),f(14, 4) / -1.20000000d-02 ,  0.00000000d+00 /
      data f(14, 5),f(14, 6) /  1.20000000d-02 ,  3.00000000d-02 /
      data f(14, 7),f(14, 8) /  5.80000000d-02 ,  1.02000000d-01 /
      data f(14, 9),f(14,10) /  1.75000000d-01 ,  2.79000000d-01 /
      data f(14,11),f(14,12) /  3.97000000d-01 ,  4.87000000d-01 /
      data f(14,13),f(14,14) /  4.87000000d-01 ,  3.97000000d-01 /
      data f(14,15) /           2.79000000d-01 /
      data f(15, 1),f(15, 2) / -3.70000000d-02 , -2.20000000d-02 /
      data f(15, 3),f(15, 4) / -1.00000000d-02 ,  0.00000000d+00 /
      data f(15, 5),f(15, 6) /  1.00000000d-02 ,  2.20000000d-02 /
      data f(15, 7),f(15, 8) /  3.70000000d-02 ,  6.00000000d-02 /
      data f(15, 9),f(15,10) /  9.30000000d-02 ,  1.33000000d-01 /
      data f(15,11),f(15,12) /  1.71000000d-01 ,  1.93000000d-01 /
      data f(15,13),f(15,14) /  1.93000000d-01 ,  1.71000000d-01 /
      data f(15,15) /           1.33000000d-01 /
      data f(16, 1),f(16, 2) / -2.40000000d-02 , -1.50000000d-02 /
      data f(16, 3),f(16, 4) / -8.00000000d-03 ,  0.00000000d+00 /
      data f(16, 5),f(16, 6) /  8.00000000d-03 ,  1.50000000d-02 /
      data f(16, 7),f(16, 8) /  2.40000000d-02 ,  3.60000000d-02 /
      data f(16, 9),f(16,10) /  5.30000000d-02 ,  7.30000000d-02 /
      data f(16,11),f(16,12) /  9.00000000d-02 ,  9.80000000d-02 /
      data f(16,13),f(16,14) /  9.80000000d-02 ,  9.00000000d-02 /
      data f(16,15) /           7.30000000d-02 /
      data f(17, 1),f(17, 2) / -1.60000000d-02 , -1.20000000d-02 /
      data f(17, 3),f(17, 4) / -6.00000000d-03 ,  0.00000000d+00 /
      data f(17, 5),f(17, 6) /  6.00000000d-03 ,  1.20000000d-02 /
      data f(17, 7),f(17, 8) /  1.60000000d-02 ,  2.20000000d-02 /
      data f(17, 9),f(17,10) /  3.10000000d-02 ,  4.10000000d-02 /
      data f(17,11),f(17,12) /  5.10000000d-02 ,  5.60000000d-02 /
      data f(17,13),f(17,14) /  5.60000000d-02 ,  5.10000000d-02 /
      data f(17,15) /           4.10000000d-02 /
      data f(18, 1),f(18, 2) / -1.10000000d-02 , -7.00000000d-03 /
      data f(18, 3),f(18, 4) / -4.00000000d-03 ,  0.00000000d+00 /
      data f(18, 5),f(18, 6) /  4.00000000d-03 ,  7.00000000d-03 /
      data f(18, 7),f(18, 8) /  1.10000000d-02 ,  1.40000000d-02 /
      data f(18, 9),f(18,10) /  1.60000000d-02 ,  1.90000000d-02 /
      data f(18,11),f(18,12) /  2.10000000d-02 ,  2.30000000d-02 /
      data f(18,13),f(18,14) /  2.30000000d-02 ,  2.10000000d-02 /
      data f(18,15) /           1.90000000d-02 /
      data f(19, 1),f(19, 2) / -8.00000000d-03 , -5.00000000d-03 /
      data f(19, 3),f(19, 4) / -2.00000000d-03 ,  0.00000000d+00 /
      data f(19, 5),f(19, 6) /  2.00000000d-03 ,  5.00000000d-03 /
      data f(19, 7),f(19, 8) /  8.00000000d-03 ,  1.10000000d-02 /
      data f(19, 9),f(19,10) /  1.40000000d-02 ,  1.70000000d-02 /
      data f(19,11),f(19,12) /  1.80000000d-02 ,  1.90000000d-02 /
      data f(19,13),f(19,14) /  1.90000000d-02 ,  1.80000000d-02 /
      data f(19,15) /           1.70000000d-02 /
      data f(20, 1),f(20, 2) / -5.00000000d-03 , -3.00000000d-03 /
      data f(20, 3),f(20, 4) / -1.00000000d-03 ,  0.00000000d+00 /
      data f(20, 5),f(20, 6) /  1.00000000d-03 ,  3.00000000d-03 /
      data f(20, 7),f(20, 8) /  5.00000000d-03 ,  7.00000000d-03 /
      data f(20, 9),f(20,10) /  9.00000000d-03 ,  1.20000000d-02 /
      data f(20,11),f(20,12) /  1.20000000d-02 ,  1.20000000d-02 /
      data f(20,13),f(20,14) /  1.20000000d-02 ,  1.20000000d-02 /
      data f(20,15) /           1.20000000d-02 /
      data f(21, 1),f(21, 2) / -3.00000000d-03 , -1.00000000d-03 /
      data f(21, 3),f(21, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(21, 5),f(21, 6) /  0.00000000d+00 ,  1.00000000d-03 /
      data f(21, 7),f(21, 8) /  3.00000000d-03 ,  4.00000000d-03 /
      data f(21, 9),f(21,10) /  5.00000000d-03 ,  5.00000000d-03 /
      data f(21,11),f(21,12) /  5.00000000d-03 ,  4.00000000d-03 /
      data f(21,13),f(21,14) /  4.00000000d-03 ,  5.00000000d-03 /
      data f(21,15) /           5.00000000d-03 /
      data f(22, 1),f(22, 2) / -1.00000000d-03 , -1.00000000d-03 /
      data f(22, 3),f(22, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(22, 5),f(22, 6) /  0.00000000d+00 ,  1.00000000d-03 /
      data f(22, 7),f(22, 8) /  1.00000000d-03 ,  2.00000000d-03 /
      data f(22, 9),f(22,10) /  2.00000000d-03 ,  1.00000000d-03 /
      data f(22,11),f(22,12) /  0.00000000d+00 , -1.00000000d-03 /
      data f(22,13),f(22,14) / -1.00000000d-03 ,  0.00000000d+00 /
      data f(22,15) /           1.00000000d-03 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 1.65435391d+01,-1.41285142d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 9.18489988d+00,-9.30297166d-03/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 4.30060283d+00,-4.21959920d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 0.00000000d+00,-3.86315259d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-4.30060283d+00, 4.37412531d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-9.18489988d+00, 8.76213029d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-1.65435391d+01, 1.61373535d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-3.24592728d+01, 2.28845562d-03/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-7.97011910d+01, 9.08882399d-03/
      data fpp( 1,10,1),fpp( 1,10,2)/ 2.44116937d+02,-3.13237516d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 3.43137527d+02, 3.31846182d-01/
      data fpp( 1,12,1),fpp( 1,12,2)/ 5.29961286d+02,-2.27100978d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 5.29961286d+02,-1.97480158d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 3.43137527d+02, 2.13362902d-01/
      data fpp( 1,15,1),fpp( 1,15,2)/ 2.44116937d+02, 4.12988549d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 1.31129218d+01,-1.33903179d-02/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 7.63020024d+00,-8.74936429d-03/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.69879433d+00,-3.75222498d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 0.00000000d+00,-1.73580549d-06/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-3.69879433d+00, 3.75916820d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-7.63020024d+00, 8.72506302d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-1.31129218d+01, 1.34805797d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-2.52814545d+01, 1.63126180d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-6.22976181d+01,-9.91105185d-03/
      data fpp( 2,10,1),fpp( 2,10,2)/ 9.01661267d+01,-1.00884106d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.67224947d+02, 2.05184694d-01/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.71177429d+02,-1.50110367d-01/
      data fpp( 2,13,1),fpp( 2,13,2)/ 3.71177429d+02,-1.32021247d-01/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.67224947d+02, 1.32828213d-01/
      data fpp( 2,15,1),fpp( 2,15,2)/ 9.01661267d+01, 2.61248393d-01/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 9.00477365d+00,-1.23742550d-02/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 5.89429917d+00,-8.04148999d-03/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 3.10421983d+00,-3.27978505d-03/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 0.00000000d+00,-1.93698117d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-3.10421983d+00, 3.35726430d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-5.89429917d+00, 7.77031262d-03/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-9.00477365d+00, 1.33814852d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-1.76149093d+01, 1.69437465d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-4.37083368d+01, 7.28352870d-03/
      data fpp( 3,10,1),fpp( 3,10,2)/-1.64381444d+02,-1.67978613d-02/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.27037314d+02, 1.17267917d-01/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.41528999d+02,-9.74938050d-02/
      data fpp( 3,13,1),fpp( 3,13,2)/ 3.41528999d+02,-8.63725432d-02/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.27037314d+02, 7.27828695d-02/
      data fpp( 3,15,1),fpp( 3,15,2)/-1.64381444d+02, 1.50021065d-01/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 7.26798360d+00,-1.12419702d-02/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 4.79260307d+00,-7.27605950d-03/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 2.48432634d+00,-2.85379174d-03/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 0.00000000d+00,-2.87735391d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-2.48432634d+00, 2.96888590d-03/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-4.79260307d+00, 6.87322996d-03/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-7.26798360d+00, 1.27381943d-02/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-1.16589084d+01, 1.66339929d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-2.92690349d+01, 1.79858340d-02/
      data fpp( 4,10,1),fpp( 4,10,2)/-6.68403525d+01,-1.74773291d-02/
      data fpp( 4,11,1),fpp( 4,11,2)/-7.07569256d+00, 7.71234824d-02/
      data fpp( 4,12,1),fpp( 4,12,2)/-1.40093423d+02,-7.60966005d-02/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.40093423d+02,-6.82119398d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/-7.07569256d+00, 4.55848400d-02/
      data fpp( 4,15,1),fpp( 4,15,2)/-6.68403525d+01, 1.00792580d-01/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 5.12329195d+00,-1.00719679d-02/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 3.73528856d+00,-6.49606429d-03/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.95847481d+00,-2.46377497d-03/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 0.00000000d+00,-2.88358165d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-1.95847481d+00, 2.57911824d-03/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-3.73528856d+00, 6.09236286d-03/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-5.12329195d+00, 1.15714303d-02/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-7.74945712d+00, 1.68619159d-02/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-1.74155237d+01, 1.93209061d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/-5.18571463d+01,-1.14554045d-03/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.22459916d+02, 2.94212557d-02/
      data fpp( 5,12,1),fpp( 5,12,2)/ 1.78446949d+01,-5.43794822d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 1.78446949d+01,-5.10721812d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/-1.22459916d+02, 1.61920518d-02/
      data fpp( 5,15,1),fpp( 5,15,2)/-5.18571463d+01, 4.84639741d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 2.59613234d+00,-7.87392423d-03/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 2.29783277d+00,-5.04215154d-03/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.43241241d+00,-1.77746960d-03/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 0.00000000d+00,-2.79700663d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.43241241d+00, 1.88934986d-03/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-2.29783277d+00, 4.65057061d-03/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-2.59613234d+00, 9.32836768d-03/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-2.57217445d+00, 1.50959587d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-4.51891143d+00, 1.89077976d-02/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.27583848d+01, 9.83285083d-03/
      data fpp( 6,11,1),fpp( 6,11,2)/-2.23824050d+01,-6.39200962d-04/
      data fpp( 6,12,1),fpp( 6,12,2)/-3.59873729d+01,-3.76960470d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.59873729d+01,-3.73722586d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/-2.23824050d+01,-1.93435470d-03/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.27583848d+01, 1.46896773d-02/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 1.29217867d+00,-5.99901099d-03/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.47338034d+00,-3.82197802d-03/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 7.11875558d-01,-1.27307694d-03/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 0.00000000d+00,-2.57142256d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-7.11875558d-01, 1.37593384d-03/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.47338034d+00, 3.46197886d-03/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-1.29217867d+00, 7.33615072d-03/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-4.11845074d-01, 1.23734183d-02/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 5.41169413d-01, 1.61301763d-02/
      data fpp( 7,10,1),fpp( 7,10,2)/-1.59314315d-01, 1.05858767d-02/
      data fpp( 7,11,1),fpp( 7,11,2)/-3.56046364d+00,-5.19368299d-03/
      data fpp( 7,12,1),fpp( 7,12,2)/ 6.10479678d+00,-2.95911447d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 6.10479678d+00,-2.97997243d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/-3.56046364d+00,-4.35936447d-03/
      data fpp( 7,15,1),fpp( 7,15,2)/-1.59314315d-01, 7.45718224d-03/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 4.85152972d-01,-4.64811653d-03/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 6.58645867d-01,-2.91376693d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 6.70085360d-01,-7.96815742d-04/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 0.00000000d+00,-1.89701025d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-6.70085360d-01, 8.72696152d-04/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-6.58645867d-01, 2.64818550d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-4.85152972d-01, 5.63456186d-03/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 4.69554747d-01, 9.67356705d-03/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 2.50423377d+00, 1.27311699d-02/
      data fpp( 8,10,1),fpp( 8,10,2)/ 5.14564211d+00, 8.88175321d-03/
      data fpp( 8,11,1),fpp( 8,11,2)/ 7.97425959d+00,-4.09818277d-03/
      data fpp( 8,12,1),fpp( 8,12,2)/ 7.66818579d+00,-2.35690221d-02/
      data fpp( 8,13,1),fpp( 8,13,2)/ 7.66818579d+00,-2.37538423d-02/
      data fpp( 8,14,1),fpp( 8,14,2)/ 7.97425959d+00,-3.35890221d-03/
      data fpp( 8,15,1),fpp( 8,15,2)/ 5.14564211d+00, 6.10945111d-03/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 2.17209440d-01,-3.43437704d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 5.42036191d-01,-2.14124592d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 2.07783002d-01,-5.40639278d-04/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 0.00000000d+00,-1.61969669d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-2.07783002d-01, 6.05427146d-04/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-5.42036191d-01, 1.91448838d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-2.17209440d-01, 4.27661932d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 9.33626086d-01, 7.37903435d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 3.09189549d+00, 9.52724328d-03/
      data fpp( 9,10,1),fpp( 9,10,2)/ 6.72674589d+00, 6.89199253d-03/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.05134253d+01,-2.77521339d-03/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.37724600d+01,-1.81711390d-02/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.37724600d+01,-1.82889764d-02/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.05134253d+01,-2.30386390d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 6.72674589d+00, 5.12443195d-03/
      data fpp(10, 1,1),fpp(10, 1,2)/-1.53990733d-01,-2.57605029d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.73209368d-01,-1.57789941d-03/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.98782633d-01,-2.92352051d-04/
      data fpp(10, 4,1),fpp(10, 4,2)/ 0.00000000d+00,-1.26923818d-05/
      data fpp(10, 5,1),fpp(10, 5,2)/-2.98782633d-01, 3.43121578d-04/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.73209368d-01, 1.40020607d-03/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.53990733d-01, 3.23605415d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.04594091d+00, 5.45557734d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/ 2.97818425d+00, 6.98163648d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 6.19737434d+00, 5.01787673d-03/
      data fpp(10,11,1),fpp(10,11,2)/ 1.02720394d+01,-1.91314340d-03/
      data fpp(10,12,1),fpp(10,12,2)/ 1.34419740d+01,-1.35053031d-02/
      data fpp(10,13,1),fpp(10,13,2)/ 1.34419740d+01,-1.35775189d-02/
      data fpp(10,14,1),fpp(10,14,2)/ 1.02720394d+01,-1.62428031d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 6.19737434d+00, 3.93464016d-03/
      data fpp(11, 1,1),fpp(11, 1,2)/-2.01246508d-01,-1.85150665d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/ 1.15126336d-01,-1.12698671d-03/
      data fpp(11, 3,1),fpp(11, 3,2)/ 9.70864679d-02,-1.80546529d-04/
      data fpp(11, 4,1),fpp(11, 4,2)/ 0.00000000d+00,-1.08271762d-05/
      data fpp(11, 5,1),fpp(11, 5,2)/-9.70864679d-02, 2.23855234d-04/
      data fpp(11, 6,1),fpp(11, 6,2)/-1.15126336d-01, 9.75406239d-04/
      data fpp(11, 7,1),fpp(11, 7,2)/ 2.01246508d-01, 2.41451981d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/ 8.82610275d-01, 3.94651452d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/ 2.54536751d+00, 5.01942212d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 5.38375676d+00, 3.51579700d-03/
      data fpp(11,11,1),fpp(11,11,2)/ 8.99841730d+00,-1.38261010d-03/
      data fpp(11,12,1),fpp(11,12,2)/ 1.21096439d+01,-9.74535658d-03/
      data fpp(11,13,1),fpp(11,13,2)/ 1.21096439d+01,-9.79262520d-03/
      data fpp(11,14,1),fpp(11,14,2)/ 8.99841730d+00,-1.19353566d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 5.38375676d+00, 2.80676783d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/-9.31844835d-02,-1.08900418d-03/
      data fpp(12, 2,1),fpp(12, 2,2)/ 6.74393005d-02,-6.61991632d-04/
      data fpp(12, 3,1),fpp(12, 3,2)/ 7.71900185d-02,-1.03029289d-04/
      data fpp(12, 4,1),fpp(12, 4,2)/ 0.00000000d+00,-5.89121404d-06/
      data fpp(12, 5,1),fpp(12, 5,2)/-7.71900185d-02, 1.26594145d-04/
      data fpp(12, 6,1),fpp(12, 6,2)/-6.74393005d-02, 5.79514635d-04/
      data fpp(12, 7,1),fpp(12, 7,2)/ 9.31844835d-02, 1.39534731d-03/
      data fpp(12, 8,1),fpp(12, 8,2)/ 6.94005142d-01, 2.41909611d-03/
      data fpp(12, 9,1),fpp(12, 9,2)/ 1.82998546d+00, 2.96826825d-03/
      data fpp(12,10,1),fpp(12,10,2)/ 3.75589458d+00, 1.96783088d-03/
      data fpp(12,11,1),fpp(12,11,2)/ 6.25724941d+00,-9.39591765d-04/
      data fpp(12,12,1),fpp(12,12,2)/ 8.47320428d+00,-5.70946382d-03/
      data fpp(12,13,1),fpp(12,13,2)/ 8.47320428d+00,-5.74181266d-03/
      data fpp(12,14,1),fpp(12,14,2)/ 6.25724941d+00,-8.10196382d-04/
      data fpp(12,15,1),fpp(12,15,2)/ 3.75589458d+00, 1.48259819d-03/
      data fpp(13, 1,1),fpp(13, 1,2)/-4.50617477d-02,-3.58739517d-04/
      data fpp(13, 2,1),fpp(13, 2,2)/ 5.91184365d-02,-2.32520965d-04/
      data fpp(13, 3,1),fpp(13, 3,2)/ 3.07400601d-02,-9.11766221d-05/
      data fpp(13, 4,1),fpp(13, 4,2)/ 0.00000000d+00,-2.77254663d-06/
      data fpp(13, 5,1),fpp(13, 5,2)/-3.07400601d-02, 1.02266809d-04/
      data fpp(13, 6,1),fpp(13, 6,2)/-5.91184365d-02, 1.93705312d-04/
      data fpp(13, 7,1),fpp(13, 7,2)/ 4.50617477d-02, 5.02911942d-04/
      data fpp(13, 8,1),fpp(13, 8,2)/ 3.37617379d-01, 9.14646919d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/ 9.68826032d-01, 1.17850038d-03/
      data fpp(13,10,1),fpp(13,10,2)/ 1.90288329d+00, 6.71351553d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 3.08975150d+00,-5.63906596d-04/
      data fpp(13,12,1),fpp(13,12,2)/ 4.18795995d+00,-2.07572517d-03/
      data fpp(13,13,1),fpp(13,13,2)/ 4.18795995d+00,-2.10087119d-03/
      data fpp(13,14,1),fpp(13,14,2)/ 3.08975150d+00,-4.63322517d-04/
      data fpp(13,15,1),fpp(13,15,2)/ 1.90288329d+00, 2.94161259d-04/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.45685257d-02,-1.40773000d-04/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.59130467d-02,-9.84539995d-05/
      data fpp(14, 3,1),fpp(14, 3,2)/-8.15025894d-03,-6.54110017d-05/
      data fpp(14, 4,1),fpp(14, 4,2)/ 0.00000000d+00, 9.80064932d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/ 8.15025894d-03, 6.50189758d-05/
      data fpp(14, 6,1),fpp(14, 6,2)/ 1.59130467d-02, 9.98260904d-05/
      data fpp(14, 7,1),fpp(14, 7,2)/ 1.45685257d-02, 1.35676663d-04/
      data fpp(14, 8,1),fpp(14, 8,2)/ 1.39525341d-01, 3.17467259d-04/
      data fpp(14, 9,1),fpp(14, 9,2)/ 4.14710415d-01, 3.34454301d-04/
      data fpp(14,10,1),fpp(14,10,2)/ 8.96572271d-01, 2.04715538d-04/
      data fpp(14,11,1),fpp(14,11,2)/ 1.44774461d+00,-3.13316453d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 1.89495594d+00,-6.31449726d-04/
      data fpp(14,13,1),fpp(14,13,2)/ 1.89495594d+00,-6.48992596d-04/
      data fpp(14,14,1),fpp(14,14,2)/ 1.44774461d+00,-2.43144973d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 8.96572271d-01,-5.84275137d-05/
      data fpp(15, 1,1),fpp(15, 1,2)/-4.06641496d-02,-4.03487639d-05/
      data fpp(15, 2,1),fpp(15, 2,2)/ 4.53375036d-03,-2.93024721d-05/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.86097566d-03,-2.24413475d-05/
      data fpp(15, 4,1),fpp(15, 4,2)/ 0.00000000d+00,-9.32137679d-07/
      data fpp(15, 5,1),fpp(15, 5,2)/-1.86097566d-03, 2.61698983d-05/
      data fpp(15, 6,1),fpp(15, 6,2)/-4.53375036d-03, 1.62525446d-05/
      data fpp(15, 7,1),fpp(15, 7,2)/ 4.06641496d-02, 8.88199233d-05/
      data fpp(15, 8,1),fpp(15, 8,2)/ 6.42812586d-02, 1.08467762d-04/
      data fpp(15, 9,1),fpp(15, 9,2)/ 1.32332307d-01, 7.73090275d-05/
      data fpp(15,10,1),fpp(15,10,2)/ 2.70827629d-01, 2.29612775d-06/
      data fpp(15,11,1),fpp(15,11,2)/ 4.79270081d-01,-2.06493538d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 6.88216309d-01,-1.36321974d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 6.88216309d-01,-1.47787309d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 4.79270081d-01,-1.60632197d-04/
      data fpp(15,15,1),fpp(15,15,2)/ 2.70827629d-01,-1.69683901d-04/
      data fpp(16, 1,1),fpp(16, 1,2)/-1.47748758d-02,-4.84808401d-05/
      data fpp(16, 2,1),fpp(16, 2,2)/-2.62219547d-02,-2.30383198d-05/
      data fpp(16, 3,1),fpp(16, 3,2)/ 7.06356311d-04, 2.06341194d-05/
      data fpp(16, 4,1),fpp(16, 4,2)/ 0.00000000d+00, 5.01842092d-07/
      data fpp(16, 5,1),fpp(16, 5,2)/-7.06356311d-04,-2.26414878d-05/
      data fpp(16, 6,1),fpp(16, 6,2)/ 2.62219547d-02, 3.00641091d-05/
      data fpp(16, 7,1),fpp(16, 7,2)/ 1.47748758d-02, 2.23850513d-05/
      data fpp(16, 8,1),fpp(16, 8,2)/ 3.53496252d-02, 6.03956856d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/ 6.39603572d-02, 3.60322063d-05/
      data fpp(16,10,1),fpp(16,10,2)/ 8.41172118d-02,-2.45245108d-05/
      data fpp(16,11,1),fpp(16,11,2)/ 1.15175070d-01,-1.17934163d-04/
      data fpp(16,12,1),fpp(16,12,2)/ 1.28178827d-01,-4.37388372d-05/
      data fpp(16,13,1),fpp(16,13,2)/ 1.28178827d-01,-4.98164070d-05/
      data fpp(16,14,1),fpp(16,14,2)/ 1.15175070d-01,-9.36238837d-05/
      data fpp(16,15,1),fpp(16,15,2)/ 8.41172118d-02,-1.15688058d-04/
      data fpp(17, 1,1),fpp(17, 1,2)/-2.02363473d-02, 3.92826325d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/ 4.35406856d-03, 2.14347350d-05/
      data fpp(17, 3,1),fpp(17, 3,2)/-4.68640090d-03,-5.02157235d-06/
      data fpp(17, 4,1),fpp(17, 4,2)/ 0.00000000d+00,-1.34844555d-06/
      data fpp(17, 5,1),fpp(17, 5,2)/ 4.68640090d-03, 1.04153545d-05/
      data fpp(17, 6,1),fpp(17, 6,2)/-4.35406856d-03,-4.03129726d-05/
      data fpp(17, 7,1),fpp(17, 7,2)/ 2.02363473d-02, 3.08365359d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/ 3.43202405d-02, 3.69668291d-05/
      data fpp(17, 9,1),fpp(17, 9,2)/ 4.38262643d-02, 1.29614763d-06/
      data fpp(17,10,1),fpp(17,10,2)/ 6.47035234d-02, 1.78485804d-05/
      data fpp(17,11,1),fpp(17,11,2)/ 6.80296383d-02,-7.26904691d-05/
      data fpp(17,12,1),fpp(17,12,2)/ 7.10683846d-02,-2.70867039d-05/
      data fpp(17,13,1),fpp(17,13,2)/ 7.10683846d-02,-3.23946536d-05/
      data fpp(17,14,1),fpp(17,14,2)/ 6.80296383d-02,-5.14586704d-05/
      data fpp(17,15,1),fpp(17,15,2)/ 6.47035234d-02,-6.17706648d-05/
      data fpp(18, 1,1),fpp(18, 1,2)/ 2.09647980d-03,-2.88799896d-05/
      data fpp(18, 2,1),fpp(18, 2,2)/-5.95122832d-03,-1.22400208d-05/
      data fpp(18, 3,1),fpp(18, 3,2)/ 1.70602455d-03, 1.78400729d-05/
      data fpp(18, 4,1),fpp(18, 4,2)/ 0.00000000d+00, 8.79729285d-07/
      data fpp(18, 5,1),fpp(18, 5,2)/-1.70602455d-03,-2.13589900d-05/
      data fpp(18, 6,1),fpp(18, 6,2)/ 5.95122832d-03, 2.45562308d-05/
      data fpp(18, 7,1),fpp(18, 7,2)/-2.09647980d-03,-1.68659332d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/-6.35534050d-04,-1.70924979d-05/
      data fpp(18, 9,1),fpp(18, 9,2)/ 1.05410286d-02, 2.52359247d-05/
      data fpp(18,10,1),fpp(18,10,2)/ 1.58308238d-02,-2.38512011d-05/
      data fpp(18,11,1),fpp(18,11,2)/ 2.63235500d-02, 1.01688795d-05/
      data fpp(18,12,1),fpp(18,12,2)/ 2.87054328d-02,-1.68243169d-05/
      data fpp(18,13,1),fpp(18,13,2)/ 2.87054328d-02,-1.46114891d-05/
      data fpp(18,14,1),fpp(18,14,2)/ 2.63235500d-02, 1.31756831d-06/
      data fpp(18,15,1),fpp(18,15,2)/ 1.58308238d-02, 9.34121584d-06/
      data fpp(19, 1,1),fpp(19, 1,2)/-1.49571906d-04, 9.23895952d-06/
      data fpp(19, 2,1),fpp(19, 2,2)/ 1.45084470d-03, 1.52208096d-06/
      data fpp(19, 3,1),fpp(19, 3,2)/-2.13769730d-03,-1.53272834d-05/
      data fpp(19, 4,1),fpp(19, 4,2)/ 0.00000000d+00,-2.12947502d-07/
      data fpp(19, 5,1),fpp(19, 5,2)/ 2.13769730d-03, 1.61790734d-05/
      data fpp(19, 6,1),fpp(19, 6,2)/-1.45084470d-03,-4.50334598d-06/
      data fpp(19, 7,1),fpp(19, 7,2)/ 1.49571906d-04, 1.83431057d-06/
      data fpp(19, 8,1),fpp(19, 8,2)/-1.77810428d-03,-2.83389628d-06/
      data fpp(19, 9,1),fpp(19, 9,2)/-7.99037881d-03, 9.50127455d-06/
      data fpp(19,10,1),fpp(19,10,2)/-8.02681876d-03,-3.51712019d-05/
      data fpp(19,11,1),fpp(19,11,2)/-1.13238381d-02, 1.11835332d-05/
      data fpp(19,12,1),fpp(19,12,2)/-1.18901160d-02,-9.56293079d-06/
      data fpp(19,13,1),fpp(19,13,2)/-1.18901160d-02,-6.90297422d-06/
      data fpp(19,14,1),fpp(19,14,2)/-1.13238381d-02, 5.43706921d-07/
      data fpp(19,15,1),fpp(19,15,2)/-8.02681876d-03, 4.72814654d-06/
      data fpp(20, 1,1),fpp(20, 1,2)/-1.49819218d-03, 9.24011810d-06/
      data fpp(20, 2,1),fpp(20, 2,2)/ 1.47849503d-04, 1.51976379d-06/
      data fpp(20, 3,1),fpp(20, 3,2)/ 8.44764658d-04,-1.53191733d-05/
      data fpp(20, 4,1),fpp(20, 4,2)/ 0.00000000d+00,-2.43070666d-07/
      data fpp(20, 5,1),fpp(20, 5,2)/-8.44764658d-04, 1.62914559d-05/
      data fpp(20, 6,1),fpp(20, 6,2)/-1.47849503d-04,-4.92275312d-06/
      data fpp(20, 7,1),fpp(20, 7,2)/ 1.49819218d-03, 3.39955655d-06/
      data fpp(20, 8,1),fpp(20, 8,2)/ 1.74795115d-03,-8.67547308d-06/
      data fpp(20, 9,1),fpp(20, 9,2)/ 3.42048660d-03, 3.13023358d-05/
      data fpp(20,10,1),fpp(20,10,2)/-1.72354880d-03,-5.65338700d-05/
      data fpp(20,11,1),fpp(20,11,2)/ 9.71802574d-04, 1.48331441d-05/
      data fpp(20,12,1),fpp(20,12,2)/ 8.55031236d-04,-2.79870643d-06/
      data fpp(20,13,1),fpp(20,13,2)/ 8.55031236d-04, 9.79547252d-07/
      data fpp(20,14,1),fpp(20,14,2)/ 9.71802574d-04,-2.79870643d-07/
      data fpp(20,15,1),fpp(20,15,2)/-1.72354880d-03, 1.39935322d-07/
      data fpp(21, 1,1),fpp(21, 1,2)/ 1.42340623d-04,-1.04086817d-05/
      data fpp(21, 2,1),fpp(21, 2,2)/-2.04224272d-03,-9.18263663d-06/
      data fpp(21, 3,1),fpp(21, 3,2)/-1.24136133d-03,-1.28607718d-05/
      data fpp(21, 4,1),fpp(21, 4,2)/ 0.00000000d+00, 6.25723848d-07/
      data fpp(21, 5,1),fpp(21, 5,2)/ 1.24136133d-03, 1.03578764d-05/
      data fpp(21, 6,1),fpp(21, 6,2)/ 2.04224272d-03, 1.79427705d-05/
      data fpp(21, 7,1),fpp(21, 7,2)/-1.42340623d-04,-2.21289584d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/ 7.86299671d-04, 1.05730631d-05/
      data fpp(21, 9,1),fpp(21, 9,2)/ 3.08432399d-04,-2.01632941d-05/
      data fpp(21,10,1),fpp(21,10,2)/ 2.92101394d-03, 1.00801131d-05/
      data fpp(21,11,1),fpp(21,11,2)/ 1.43662784d-03,-2.01571584d-05/
      data fpp(21,12,1),fpp(21,12,2)/ 2.46999108d-03, 1.05485205d-05/
      data fpp(21,13,1),fpp(21,13,2)/ 2.46999108d-03, 8.43301784d-06/
      data fpp(21,14,1),fpp(21,14,2)/ 1.43662784d-03,-1.16951480d-05/
      data fpp(21,15,1),fpp(21,15,2)/ 2.92101394d-03,-2.16524260d-05/
      data fpp(22, 1,1),fpp(22, 1,2)/ 9.28829689d-04, 2.88871143d-05/
      data fpp(22, 2,1),fpp(22, 2,2)/-3.97887864d-03, 1.22257714d-05/
      data fpp(22, 3,1),fpp(22, 3,2)/-1.87931933d-03,-1.77902000d-05/
      data fpp(22, 4,1),fpp(22, 4,2)/ 0.00000000d+00,-1.06497156d-06/
      data fpp(22, 5,1),fpp(22, 5,2)/ 1.87931933d-03, 2.20500862d-05/
      data fpp(22, 6,1),fpp(22, 6,2)/ 3.97887864d-03,-2.71353733d-05/
      data fpp(22, 7,1),fpp(22, 7,2)/-9.28829689d-04, 2.64914069d-05/
      data fpp(22, 8,1),fpp(22, 8,2)/ 1.10685016d-03,-1.88302543d-05/
      data fpp(22, 9,1),fpp(22, 9,2)/ 1.34578380d-03,-1.11703896d-05/
      data fpp(22,10,1),fpp(22,10,2)/ 8.03949303d-03, 3.51181286d-06/
      data fpp(22,11,1),fpp(22,11,2)/ 5.28168608d-03,-2.87686179d-06/
      data fpp(22,12,1),fpp(22,12,2)/ 7.26500446d-03, 7.99563430d-06/
      data fpp(22,13,1),fpp(22,13,2)/ 7.26500446d-03, 7.45152800d-06/
      data fpp(22,14,1),fpp(22,14,2)/ 5.28168608d-03,-7.00436570d-07/
      data fpp(22,15,1),fpp(22,15,2)/ 8.03949303d-03,-4.64978171d-06/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 4.35890129d-02, 2.48651980d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/ 5.41072495d-03,-1.14864448d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/ 4.05350542d-02,-1.85675425d-01/
      data fpppp( 1, 7),fpppp( 1, 8)/ 5.53706116d-01,-2.54257471d+00/
      data fpppp( 1, 9),fpppp( 1,10)/ 7.73702164d+00,-6.14190910d+00/
      data fpppp( 1,11),fpppp( 1,12)/ 3.34276250d+00,-1.96095077d+00/
      data fpppp( 1,13),fpppp( 1,14)/-1.39324170d+00, 1.07192621d+00/
      data fpppp( 1,15) /             2.37372699d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 2.85108557d-02, 1.58918411d-02/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.00072062d-03,-5.93802949d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/ 2.27513973d-02,-9.90242539d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 2.80266678d-01,-1.42319112d+00/
      data fpppp( 2, 9),fpppp( 2,10)/ 3.92163996d+00,-2.89457420d+00/
      data fpppp( 2,11),fpppp( 2,12)/ 3.13236137d+00,-2.02125154d+00/
      data fpppp( 2,13),fpppp( 2,14)/-1.62100052d+00, 1.53135727d+00/
      data fpppp( 2,15) /             3.10919118d+00 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 9.08748099d-03, 4.12760419d-03/
      data fpppp( 3, 3),fpppp( 3, 4)/-6.37419006d-03, 2.52072675d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-3.70871696d-03, 3.11625704d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-1.40165272d-01, 1.99518849d-01/
      data fpppp( 3, 9),fpppp( 3,10)/-1.70690763d+00, 9.53330929d-01/
      data fpppp( 3,11),fpppp( 3,12)/ 7.37461813d+00,-4.57847250d+00/
      data fpppp( 3,13),fpppp( 3,14)/-4.00888091d+00, 5.09625177d+00/
      data fpppp( 3,15) /             9.49720478d+00 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 4.87236775d-03, 2.13144735d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/-3.37192850d-03, 7.93289724d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.98769607d-04, 8.97460878d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/-4.61234334d-02, 6.05864687d-02/
      data fpppp( 4, 9),fpppp( 4,10)/-9.89374542d-01, 2.69924023d+00/
      data fpppp( 4,11),fpppp( 4,12)/-3.96742772d+00, 1.60352719d+00/
      data fpppp( 4,13),fpppp( 4,14)/ 1.16366421d+00,-2.20797577d+00/
      data fpppp( 4,15) /            -3.89870456d+00 /
      data fpppp( 5, 1),fpppp( 5, 2)/-5.94637005d-03,-3.91455742d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.72402248d-03,-8.90156478d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/ 2.08008507d-03, 2.66833835d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/ 1.05751837d-02,-1.19258780d-01/
      data fpppp( 5, 9),fpppp( 5,10)/ 4.40658498d-02,-1.54353798d+00/
      data fpppp( 5,11),fpppp( 5,12)/ 3.96041722d+00,-1.64368805d+00/
      data fpppp( 5,13),fpppp( 5,14)/-1.25828278d+00, 2.41879615d+00/
      data fpppp( 5,15) /             4.23754105d+00 /
      data fpppp( 6, 1),fpppp( 6, 2)/-5.89124857d-03,-5.23370173d-03/
      data fpppp( 6, 3),fpppp( 6, 4)/-7.20119218d-03, 1.89478408d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/ 7.12540082d-03, 5.49897151d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/ 4.90596084d-03,-5.78736716d-03/
      data fpppp( 6, 9),fpppp( 6,10)/-9.99981840d-02, 2.82159164d-02/
      data fpppp( 6,11),fpppp( 6,12)/-9.59382872d-02, 1.16680371d-01/
      data fpppp( 6,13),fpppp( 6,14)/ 1.06077068d-01,-5.35250763d-02/
      data fpppp( 6,15) /            -1.30833624d-01 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.89581939d-02,-1.02115192d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/ 3.24188370d-03, 2.21737866d-04/
      data fpppp( 7, 5),fpppp( 7, 6)/-4.12883516d-03, 1.33158493d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/ 7.42782486d-03,-1.07923306d-03/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.24996075d-03,-1.03130503d-01/
      data fpppp( 7,11),fpppp( 7,12)/ 2.49232115d-01,-1.09813371d-01/
      data fpppp( 7,13),fpppp( 7,14)/-8.51337578d-02, 1.50513663d-01/
      data fpppp( 7,15) /             2.67063690d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 3.12120717d-03,-7.14587401d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/-9.98606171d-03,-2.32656950d-04/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.09166895d-02,-2.54260991d-03/
      data fpppp( 8, 7),fpppp( 8, 8)/ 8.97695425d-03, 1.35076824d-02/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.79059473d-03, 1.57336972d-02/
      data fpppp( 8,11),fpppp( 8,12)/-5.34928341d-02, 1.01561617d-02/
      data fpppp( 8,13),fpppp( 8,14)/ 5.46014605d-03,-3.47087716d-02/
      data fpppp( 8,15) /            -5.47065373d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.40998791d-02,-7.28364258d-03/
      data fpppp( 9, 3),fpppp( 9, 4)/ 3.68965300d-03, 1.13241812d-04/
      data fpppp( 9, 5),fpppp( 9, 6)/-4.14262025d-03, 8.86902795d-03/
      data fpppp( 9, 7),fpppp( 9, 8)/ 8.21130488d-03, 7.84627905d-03/
      data fpppp( 9, 9),fpppp( 9,10)/ 2.08496119d-02,-2.64986757d-03/
      data fpppp( 9,11),fpppp( 9,12)/-1.14040290d-03,-2.44471959d-02/
      data fpppp( 9,13),fpppp( 9,14)/-2.38592544d-02,-3.49216899d-03/
      data fpppp( 9,15) /             6.16925527d-03 /
      data fpppp(10, 1),fpppp(10, 2)/-3.36030225d-05,-1.52701785d-03/
      data fpppp(10, 3),fpppp(10, 4)/-5.95593577d-03,-1.10592884d-04/
      data fpppp(10, 5),fpppp(10, 6)/ 6.39830730d-03,-2.12825166d-05/
      data fpppp(10, 7),fpppp(10, 8)/ 5.78443297d-03, 1.07685552d-02/
      data fpppp(10, 9),fpppp(10,10)/ 1.35589360d-02, 1.22125061d-02/
      data fpppp(10,11),fpppp(10,12)/-1.10804645d-02,-2.21744704d-02/
      data fpppp(10,13),fpppp(10,14)/-2.30343961d-02,-7.64076172d-03/
      data fpppp(10,15) /            -6.86379197d-04 /
      data fpppp(11, 1),fpppp(11, 2)/-5.83317394d-03,-3.47335577d-03/
      data fpppp(11, 3),fpppp(11, 4)/-3.38165754d-04, 8.32228025d-05/
      data fpppp(11, 5),fpppp(11, 6)/ 5.27454378d-06, 4.63847500d-03/
      data fpppp(11, 7),fpppp(11, 8)/ 1.50558821d-03, 1.12386275d-02/
      data fpppp(11, 9),fpppp(11,10)/ 1.24235100d-02, 9.60525269d-03/
      data fpppp(11,11),fpppp(11,12)/-4.26824273d-03,-2.27383180d-02/
      data fpppp(11,13),fpppp(11,14)/-2.29877228d-02,-3.27062369d-03/
      data fpppp(11,15) /             5.86418129d-03 /
      data fpppp(12, 1),fpppp(12, 2)/-2.15455192d-03,-1.49573475d-03/
      data fpppp(12, 3),fpppp(12, 4)/-9.14893069d-04,-6.11371642d-05/
      data fpppp(12, 5),fpppp(12, 6)/ 1.15944173d-03, 6.39814447d-04/
      data fpppp(12, 7),fpppp(12, 8)/ 5.33368445d-03, 4.43726022d-03/
      data fpppp(12, 9),fpppp(12,10)/ 9.02685396d-03, 6.85105257d-03/
      data fpppp(12,11),fpppp(12,12)/-1.90432175d-03,-1.63577633d-02/
      data fpppp(12,13),fpppp(12,14)/-1.64531953d-02,-1.52259360d-03/
      data fpppp(12,15) /             5.41957200d-03 /
      data fpppp(13, 1),fpppp(13, 2)/-2.57927602d-03,-1.42214231d-03/
      data fpppp(13, 3),fpppp(13, 4)/ 3.14331629d-04, 2.31147741d-05/
      data fpppp(13, 5),fpppp(13, 6)/-4.06790725d-04, 1.74574915d-03/
      data fpppp(13, 7),fpppp(13, 8)/ 1.37730777d-03, 4.04754660d-03/
      data fpppp(13, 9),fpppp(13,10)/ 2.75168711d-03, 3.11662107d-03/
      data fpppp(13,11),fpppp(13,12)/-4.95142050d-05,-8.23814975d-03/
      data fpppp(13,13),fpppp(13,14)/-8.20704715d-03,-1.73924581d-04/
      data fpppp(13,15) /             3.58315998d-03 /
      data fpppp(14, 1),fpppp(14, 2)/ 1.76500427d-04, 9.46151647d-05/
      data fpppp(14, 3),fpppp(14, 4)/-8.52255503d-06,-3.72766749d-05/
      data fpppp(14, 5),fpppp(14, 6)/ 1.57629255d-04,-6.16488614d-04/
      data fpppp(14, 7),fpppp(14, 8)/ 1.76188667d-03, 1.14702209d-03/
      data fpppp(14, 9),fpppp(14,10)/ 2.66372055d-03, 5.98702571d-04/
      data fpppp(14,11),fpppp(14,12)/-8.99902068d-04,-3.23675458d-03/
      data fpppp(14,13),fpppp(14,14)/-3.25612515d-03,-8.22419764d-04/
      data fpppp(14,15) /             3.08143931d-04 /
      data fpppp(15, 1),fpppp(15, 2)/-9.47740426d-04,-5.14288865d-04/
      data fpppp(15, 3),fpppp(15, 4)/ 1.32655404d-04, 3.23751906d-05/
      data fpppp(15, 5),fpppp(15, 6)/-2.62156167d-04, 9.67541534d-04/
      data fpppp(15, 7),fpppp(15, 8)/-7.35769485d-04, 6.80688941d-04/
      data fpppp(15, 9),fpppp(15,10)/ 6.79050085d-04, 8.29767166d-04/
      data fpppp(15,11),fpppp(15,12)/ 1.98709002d-04,-1.59437658d-03/
      data fpppp(15,13),fpppp(15,14)/-1.58461162d-03, 1.59649172d-04/
      data fpppp(15,15) /             9.76241530d-04 /
      data fpppp(16, 1),fpppp(16, 2)/ 1.00910978d-03, 4.53243295d-04/
      data fpppp(16, 3),fpppp(16, 4)/-5.19559556d-04,-3.30851128d-05/
      data fpppp(16, 5),fpppp(16, 6)/ 6.51900007d-04,-9.16434874d-04/
      data fpppp(16, 7),fpppp(16, 8)/ 7.11316089d-04,-7.51977528d-06/
      data fpppp(16, 9),fpppp(16,10)/-1.99078038d-04, 2.96599286d-04/
      data fpppp(16,11),fpppp(16,12)/-3.33258881d-04,-4.68098825d-05/
      data fpppp(16,13),fpppp(16,14)/-8.30536024d-05,-1.88284001d-04/
      data fpppp(16,15) /            -2.47056513d-04 /
      data fpppp(17, 1),fpppp(17, 2)/-7.87205510d-04,-3.81670655d-04/
      data fpppp(17, 3),fpppp(17, 4)/ 2.96035008d-04, 2.11428428d-05/
      data fpppp(17, 5),fpppp(17, 6)/-3.80606380d-04, 6.77670454d-04/
      data fpppp(17, 7),fpppp(17, 8)/-3.12222315d-04,-5.91725543d-05/
      data fpppp(17, 9),fpppp(17,10)/ 2.74240369d-04,-3.55514797d-04/
      data fpppp(17,11),fpppp(17,12)/ 9.47501647d-05,-4.07279766d-05/
      data fpppp(17,13),fpppp(17,14)/-1.63535418d-05,-2.74757483d-06/
      data fpppp(17,15) /             1.01017257d-05 /
      data fpppp(18, 1),fpppp(18, 2)/ 3.94924237d-04, 1.82665124d-04/
      data fpppp(18, 3),fpppp(18, 4)/-1.83287075d-04,-1.13134713d-05/
      data fpppp(18, 5),fpppp(18, 6)/ 2.28540960d-04,-3.41053723d-04/
      data fpppp(18, 7),fpppp(18, 8)/ 1.93376272d-04, 1.38067866d-04/
      data fpppp(18, 9),fpppp(18,10)/-1.62710720d-04, 1.59568966d-04/
      data fpppp(18,11),fpppp(18,12)/-1.63389289d-04, 7.33759474d-06/
      data fpppp(18,13),fpppp(18,14)/-1.17746268d-05,-8.69404024d-05/
      data fpppp(18,15) /            -1.27114357d-04 /
      data fpppp(19, 1),fpppp(19, 2)/-1.54802010d-04,-6.43686969d-05/
      data fpppp(19, 3),fpppp(19, 4)/ 1.00939281d-04, 4.18593281d-06/
      data fpppp(19, 5),fpppp(19, 6)/-1.17683012d-04, 1.22971756d-04/
      data fpppp(19, 7),fpppp(19, 8)/-6.28664963d-05,-8.31913386d-05/
      data fpppp(19, 9),fpppp(19,10)/ 1.38555949d-04,-1.00482384d-04/
      data fpppp(19,11),fpppp(19,12)/ 6.77388202d-05,-6.62840812d-06/
      data fpppp(19,13),fpppp(19,14)/ 3.00415093d-06, 2.92085840d-05/
      data fpppp(19,15) /             4.40060019d-05 /
      data fpppp(20, 1),fpppp(20, 2)/-4.42293981d-06,-7.77685044d-06/
      data fpppp(20, 3),fpppp(20, 4)/-2.14172500d-05, 9.45061805d-07/
      data fpppp(20, 5),fpppp(20, 6)/ 1.76370028d-05, 2.10077157d-05/
      data fpppp(20, 7),fpppp(20, 8)/-4.47202741d-05, 7.40964179d-05/
      data fpppp(20, 9),fpppp(20,10)/-1.66298809d-04, 1.82104566d-04/
      data fpppp(20,11),fpppp(20,12)/-9.17562489d-05, 1.61930673d-05/
      data fpppp(20,13),fpppp(20,14)/ 8.02062737d-07,-3.01922307d-05/
      data fpppp(20,15) /            -4.87605022d-05 /
      data fpppp(21, 1),fpppp(21, 2)/ 5.45940190d-05, 3.12755991d-05/
      data fpppp(21, 3),fpppp(21, 4)/-5.68532081d-07,-2.57267396d-06/
      data fpppp(21, 5),fpppp(21, 6)/ 1.08592279d-05,-6.72930345d-05/
      data fpppp(21, 7),fpppp(21, 8)/ 7.91850269d-05,-6.26536550d-05/
      data fpppp(21, 9),fpppp(21,10)/ 8.70391393d-05,-1.00075973d-04/
      data fpppp(21,11),fpppp(21,12)/ 6.74466949d-05,-1.86458458d-05/
      data fpppp(21,13),fpppp(21,14)/-8.78670719d-06, 2.80101404d-05/
      data fpppp(21,15) /             4.78111063d-05 /
      data fpppp(22, 1),fpppp(22, 2)/ 1.39714438d-04, 7.53393049d-05/
      data fpppp(22, 3),fpppp(22, 4)/-2.06355994d-05,-6.01130551d-06/
      data fpppp(22, 5),fpppp(22, 6)/ 4.46808215d-05,-1.59497582d-04/
      data fpppp(22, 7),fpppp(22, 8)/ 1.72873448d-04,-1.15392920d-04/
      data fpppp(22, 9),fpppp(22,10)/ 1.80893457d-04,-2.20894373d-04/
      data fpppp(22,11),fpppp(22,12)/ 1.35593064d-04,-3.70103622d-05/
      data fpppp(22,13),fpppp(22,14)/-1.62649966d-05, 5.26116013d-05/
      data fpppp(22,15) /             9.02861111d-05 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /

      iprint=0

      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
      subroutine ch3po_1app(r,theta,phi,energy)
      implicit real*8 (a-h,o-z)
      data degrad / 57.29577951308232d 00 /
      call ch3po_1app_s(r,theta,es)
      call ch3po_1app_d(r,theta,ed)
      energy = es + 0.5d0 * cos(3.0d0*phi/degrad)*ed 
      return
      end

      subroutine ch3po_1app_s(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     1(2a") surface
c     average of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -2.82190000d+01 , -3.62010000d+01 /
      data f( 1, 3),f( 1, 4) / -4.09180000d+01 , -4.24710000d+01 /
      data f( 1, 5),f( 1, 6) / -4.09180000d+01 , -3.62010000d+01 /
      data f( 1, 7),f( 1, 8) / -2.82190000d+01 , -1.67670000d+01 /
      data f( 1, 9),f( 1,10) / -1.49400000d+00 ,  1.77970000d+01 /
      data f( 1,11),f( 1,12) /  4.18920000d+01 ,  6.59770000d+01 /
      data f( 1,13),f( 1,14) /  6.59770000d+01 ,  4.18920000d+01 /
      data f( 1,15) /           1.77970000d+01 /
      data f( 2, 1),f( 2, 2) / -2.64310000d+01 , -3.35960000d+01 /
      data f( 2, 3),f( 2, 4) / -3.78340000d+01 , -3.92290000d+01 /
      data f( 2, 5),f( 2, 6) / -3.78340000d+01 , -3.35960000d+01 /
      data f( 2, 7),f( 2, 8) / -2.64310000d+01 , -1.61650000d+01 /
      data f( 2, 9),f( 2,10) / -2.53600000d+00 ,  1.45600000d+01 /
      data f( 2,11),f( 2,12) /  3.49150000d+01 ,  5.41660000d+01 /
      data f( 2,13),f( 2,14) /  5.41660000d+01 ,  3.49150000d+01 /
      data f( 2,15) /           1.45600000d+01 /
      data f( 3, 1),f( 3, 2) / -2.45170000d+01 , -3.09280000d+01 /
      data f( 3, 3),f( 3, 4) / -3.47240000d+01 , -3.59730000d+01 /
      data f( 3, 5),f( 3, 6) / -3.47240000d+01 , -3.09280000d+01 /
      data f( 3, 7),f( 3, 8) / -2.45170000d+01 , -1.53450000d+01 /
      data f( 3, 9),f( 3,10) / -3.22200000d+00 ,  1.19020000d+01 /
      data f( 3,11),f( 3,12) /  2.91690000d+01 ,  4.49470000d+01 /
      data f( 3,13),f( 3,14) /  4.49470000d+01 ,  2.91690000d+01 /
      data f( 3,15) /           1.19020000d+01 /
      data f( 4, 1),f( 4, 2) / -2.25400000d+01 , -2.82580000d+01 /
      data f( 4, 3),f( 4, 4) / -3.16460000d+01 , -3.27620000d+01 /
      data f( 4, 5),f( 4, 6) / -3.16460000d+01 , -2.82580000d+01 /
      data f( 4, 7),f( 4, 8) / -2.25400000d+01 , -1.43750000d+01 /
      data f( 4, 9),f( 4,10) / -3.62500000d+00 ,  9.67200000d+00 /
      data f( 4,11),f( 4,12) /  2.45420000d+01 ,  3.75910000d+01 /
      data f( 4,13),f( 4,14) /  3.75910000d+01 ,  2.45420000d+01 /
      data f( 4,15) /           9.67200000d+00 /
      data f( 5, 1),f( 5, 2) / -2.05580000d+01 , -2.56390000d+01 /
      data f( 5, 3),f( 5, 4) / -2.86540000d+01 , -2.91760000d+01 /
      data f( 5, 5),f( 5, 6) / -2.86540000d+01 , -2.56390000d+01 /
      data f( 5, 7),f( 5, 8) / -2.05580000d+01 , -1.33160000d+01 /
      data f( 5, 9),f( 5,10) / -3.81700000d+00 ,  7.83200000d+00 /
      data f( 5,11),f( 5,12) /  2.07450000d+01 ,  3.16890000d+01 /
      data f( 5,13),f( 5,14) /  3.16890000d+01 ,  2.07450000d+01 /
      data f( 5,15) /           7.83200000d+00 /
      data f( 6, 1),f( 6, 2) / -1.67530000d+01 , -2.07180000d+01 /
      data f( 6, 3),f( 6, 4) / -2.30750000d+01 , -2.38510000d+01 /
      data f( 6, 5),f( 6, 6) / -2.30750000d+01 , -2.07180000d+01 /
      data f( 6, 7),f( 6, 8) / -1.67530000d+01 , -1.11260000d+01 /
      data f( 6, 9),f( 6,10) / -3.80200000d+00 ,  5.07800000d+00 /
      data f( 6,11),f( 6,12) /  1.47320000d+01 ,  2.30190000d+01 /
      data f( 6,13),f( 6,14) /  2.30190000d+01 ,  1.47320000d+01 /
      data f( 6,15) /           5.07800000d+00 /
      data f( 7, 1),f( 7, 2) / -1.33640000d+01 , -1.64050000d+01 /
      data f( 7, 3),f( 7, 4) / -1.82150000d+01 , -1.88100000d+01 /
      data f( 7, 5),f( 7, 6) / -1.82150000d+01 , -1.64050000d+01 /
      data f( 7, 7),f( 7, 8) / -1.33640000d+01 , -9.07000000d+00 /
      data f( 7, 9),f( 7,10) / -3.52000000d+00 ,  3.15100000d+00 /
      data f( 7,11),f( 7,12) /  1.03450000d+01 ,  1.65320000d+01 /
      data f( 7,13),f( 7,14) /  1.65320000d+01 ,  1.03450000d+01 /
      data f( 7,15) /           3.15100000d+00 /
      data f( 8, 1),f( 8, 2) / -1.05140000d+01 , -1.28030000d+01 /
      data f( 8, 3),f( 8, 4) / -1.41680000d+01 , -1.46160000d+01 /
      data f( 8, 5),f( 8, 6) / -1.41680000d+01 , -1.28030000d+01 /
      data f( 8, 7),f( 8, 8) / -1.05140000d+01 , -7.29400000d+00 /
      data f( 8, 9),f( 8,10) / -3.15800000d+00 ,  1.78300000d+00 /
      data f( 8,11),f( 8,12) /  7.09100000d+00 ,  1.16880000d+01 /
      data f( 8,13),f( 8,14) /  1.16880000d+01 ,  7.09100000d+00 /
      data f( 8,15) /           1.78300000d+00 /
      data f( 9, 1),f( 9, 2) / -8.20900000d+00 , -9.90200000d+00 /
      data f( 9, 3),f( 9, 4) / -1.09150000d+01 , -1.12490000d+01 /
      data f( 9, 5),f( 9, 6) / -1.09150000d+01 , -9.90200000d+00 /
      data f( 9, 7),f( 9, 8) / -8.20900000d+00 , -5.83200000d+00 /
      data f( 9, 9),f( 9,10) / -2.79300000d+00 ,  8.21000000d-01 /
      data f( 9,11),f( 9,12) /  4.69600000d+00 ,  8.07600000d+00 /
      data f( 9,13),f( 9,14) /  8.07600000d+00 ,  4.69600000d+00 /
      data f( 9,15) /           8.21000000d-01 /
      data f(10, 1),f(10, 2) / -6.39000000d+00 , -7.62700000d+00 /
      data f(10, 3),f(10, 4) / -8.36900000d+00 , -8.61600000d+00 /
      data f(10, 5),f(10, 6) / -8.36900000d+00 , -7.62700000d+00 /
      data f(10, 7),f(10, 8) / -6.39000000d+00 , -4.65900000d+00 /
      data f(10, 9),f(10,10) / -2.45100000d+00 ,  1.65000000d-01 /
      data f(10,11),f(10,12) /  2.96600000d+00 ,  5.43000000d+00 /
      data f(10,13),f(10,14) /  5.43000000d+00 ,  2.96600000d+00 /
      data f(10,15) /           1.65000000d-01 /
      data f(11, 1),f(11, 2) / -4.97600000d+00 , -5.87100000d+00 /
      data f(11, 3),f(11, 4) / -6.41100000d+00 , -6.59100000d+00 /
      data f(11, 5),f(11, 6) / -6.41100000d+00 , -5.87100000d+00 /
      data f(11, 7),f(11, 8) / -4.97600000d+00 , -3.72800000d+00 /
      data f(11, 9),f(11,10) / -2.13800000d+00 , -2.60000000d-01 /
      data f(11,11),f(11,12) /  1.74700000d+00 ,  3.52900000d+00 /
      data f(11,13),f(11,14) /  3.52900000d+00 ,  1.74700000d+00 /
      data f(11,15) /          -2.60000000d-01 /
      data f(12, 1),f(12, 2) / -3.43300000d+00 , -3.97800000d+00 /
      data f(12, 3),f(12, 4) / -4.31200000d+00 , -4.42400000d+00 /
      data f(12, 5),f(12, 6) / -4.31200000d+00 , -3.97800000d+00 /
      data f(12, 7),f(12, 8) / -3.43300000d+00 , -2.68000000d+00 /
      data f(12, 9),f(12,10) / -1.72400000d+00 , -5.99000000d-01 /
      data f(12,11),f(12,12) /  6.00000000d-01 ,  1.68200000d+00 /
      data f(12,13),f(12,14) /  1.68200000d+00 ,  6.00000000d-01 /
      data f(12,15) /          -5.99000000d-01 /
      data f(13, 1),f(13, 2) / -1.88700000d+00 , -2.12300000d+00 /
      data f(13, 3),f(13, 4) / -2.27100000d+00 , -2.32100000d+00 /
      data f(13, 5),f(13, 6) / -2.27100000d+00 , -2.12300000d+00 /
      data f(13, 7),f(13, 8) / -1.88700000d+00 , -1.57200000d+00 /
      data f(13, 9),f(13,10) / -1.18100000d+00 , -7.26000000d-01 /
      data f(13,11),f(13,12) / -2.45000000d-01 ,  1.99000000d-01 /
      data f(13,13),f(13,14) /  1.99000000d-01 , -2.45000000d-01 /
      data f(13,15) /          -7.26000000d-01 /
      data f(14, 1),f(14, 2) / -1.08200000d+00 , -1.18300000d+00 /
      data f(14, 3),f(14, 4) / -1.24800000d+00 , -1.26900000d+00 /
      data f(14, 5),f(14, 6) / -1.24800000d+00 , -1.18300000d+00 /
      data f(14, 7),f(14, 8) / -1.08200000d+00 , -9.52000000d-01 /
      data f(14, 9),f(14,10) / -8.00000000d-01 , -6.33000000d-01 /
      data f(14,11),f(14,12) / -4.69000000d-01 , -3.31000000d-01 /
      data f(14,13),f(14,14) / -3.31000000d-01 , -4.69000000d-01 /
      data f(14,15) /          -6.33000000d-01 /
      data f(15, 1),f(15, 2) / -6.53000000d-01 , -6.96000000d-01 /
      data f(15, 3),f(15, 4) / -7.24000000d-01 , -7.34000000d-01 /
      data f(15, 5),f(15, 6) / -7.24000000d-01 , -6.96000000d-01 /
      data f(15, 7),f(15, 8) / -6.53000000d-01 , -6.01000000d-01 /
      data f(15, 9),f(15,10) / -5.46000000d-01 , -4.95000000d-01 /
      data f(15,11),f(15,12) / -4.60000000d-01 , -4.42000000d-01 /
      data f(15,13),f(15,14) / -4.42000000d-01 , -4.60000000d-01 /
      data f(15,15) /          -4.95000000d-01 /
      data f(16, 1),f(16, 2) / -4.15000000d-01 , -4.32000000d-01 /
      data f(16, 3),f(16, 4) / -4.44000000d-01 , -4.48000000d-01 /
      data f(16, 5),f(16, 6) / -4.44000000d-01 , -4.32000000d-01 /
      data f(16, 7),f(16, 8) / -4.15000000d-01 , -3.97000000d-01 /
      data f(16, 9),f(16,10) / -3.82000000d-01 , -3.76000000d-01 /
      data f(16,11),f(16,12) / -3.82000000d-01 , -3.97000000d-01 /
      data f(16,13),f(16,14) / -3.97000000d-01 , -3.82000000d-01 /
      data f(16,15) /          -3.76000000d-01 /
      data f(17, 1),f(17, 2) / -2.75000000d-01 , -2.80000000d-01 /
      data f(17, 3),f(17, 4) / -2.85000000d-01 , -2.87000000d-01 /
      data f(17, 5),f(17, 6) / -2.85000000d-01 , -2.80000000d-01 /
      data f(17, 7),f(17, 8) / -2.75000000d-01 , -2.72000000d-01 /
      data f(17, 9),f(17,10) / -2.74000000d-01 , -2.84000000d-01 /
      data f(17,11),f(17,12) / -3.01000000d-01 , -3.20000000d-01 /
      data f(17,13),f(17,14) / -3.20000000d-01 , -3.01000000d-01 /
      data f(17,15) /          -2.84000000d-01 /
      data f(18, 1),f(18, 2) / -1.37000000d-01 , -1.35000000d-01 /
      data f(18, 3),f(18, 4) / -1.34000000d-01 , -1.35000000d-01 /
      data f(18, 5),f(18, 6) / -1.34000000d-01 , -1.35000000d-01 /
      data f(18, 7),f(18, 8) / -1.37000000d-01 , -1.41000000d-01 /
      data f(18, 9),f(18,10) / -1.49000000d-01 , -1.62000000d-01 /
      data f(18,11),f(18,12) / -1.77000000d-01 , -1.92000000d-01 /
      data f(18,13),f(18,14) / -1.92000000d-01 , -1.77000000d-01 /
      data f(18,15) /          -1.62000000d-01 /
      data f(19, 1),f(19, 2) / -7.50000000d-02 , -7.30000000d-02 /
      data f(19, 3),f(19, 4) / -7.20000000d-02 , -7.20000000d-02 /
      data f(19, 5),f(19, 6) / -7.20000000d-02 , -7.30000000d-02 /
      data f(19, 7),f(19, 8) / -7.50000000d-02 , -7.90000000d-02 /
      data f(19, 9),f(19,10) / -8.60000000d-02 , -9.50000000d-02 /
      data f(19,11),f(19,12) / -1.06000000d-01 , -1.13000000d-01 /
      data f(19,13),f(19,14) / -1.13000000d-01 , -1.06000000d-01 /
      data f(19,15) /          -9.50000000d-02 /
      data f(20, 1),f(20, 2) / -4.10000000d-02 , -4.00000000d-02 /
      data f(20, 3),f(20, 4) / -4.00000000d-02 , -3.90000000d-02 /
      data f(20, 5),f(20, 6) / -4.00000000d-02 , -4.00000000d-02 /
      data f(20, 7),f(20, 8) / -4.10000000d-02 , -4.50000000d-02 /
      data f(20, 9),f(20,10) / -5.10000000d-02 , -5.80000000d-02 /
      data f(20,11),f(20,12) / -6.40000000d-02 , -6.90000000d-02 /
      data f(20,13),f(20,14) / -6.90000000d-02 , -6.40000000d-02 /
      data f(20,15) /          -5.80000000d-02 /
      data f(21, 1),f(21, 2) / -2.20000000d-02 , -2.10000000d-02 /
      data f(21, 3),f(21, 4) / -2.10000000d-02 , -2.10000000d-02 /
      data f(21, 5),f(21, 6) / -2.10000000d-02 , -2.10000000d-02 /
      data f(21, 7),f(21, 8) / -2.20000000d-02 , -2.50000000d-02 /
      data f(21, 9),f(21,10) / -2.90000000d-02 , -3.40000000d-02 /
      data f(21,11),f(21,12) / -3.90000000d-02 , -4.20000000d-02 /
      data f(21,13),f(21,14) / -4.20000000d-02 , -3.90000000d-02 /
      data f(21,15) /          -3.40000000d-02 /
      data f(22, 1),f(22, 2) / -1.10000000d-02 , -1.00000000d-02 /
      data f(22, 3),f(22, 4) / -9.00000000d-03 , -1.00000000d-02 /
      data f(22, 5),f(22, 6) / -9.00000000d-03 , -1.00000000d-02 /
      data f(22, 7),f(22, 8) / -1.10000000d-02 , -1.30000000d-02 /
      data f(22, 9),f(22,10) / -1.60000000d-02 , -1.90000000d-02 /
      data f(22,11),f(22,12) / -2.40000000d-02 , -2.60000000d-02 /
      data f(22,13),f(22,14) / -2.60000000d-02 , -2.40000000d-02 /
      data f(22,15) /          -1.90000000d-02 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 1.88894732d+01, 3.36540304d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.23741203d+01, 3.26619391d-02/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 8.39847284d+00, 3.15982131d-02/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.29897735d+00, 3.07852085d-02/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 8.39847284d+00, 3.16209530d-02/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.23741203d+01, 3.25709796d-02/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.88894732d+01, 3.39951286d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.85781547d+01, 3.96485060d-02/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 4.29015482d+01, 3.66708474d-02/
      data fpp( 1,10,1),fpp( 1,10,2)/ 7.25127124d+01, 5.47481044d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.35003235d+02, 3.25767350d-02/
      data fpp( 1,12,1),fpp( 1,12,2)/ 3.30930264d+02,-1.85655044d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.30930264d+02,-1.81873234d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.35003235d+02, 1.74494956d-02/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.25127124d+01, 1.11475252d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 1.26210537d+01, 3.01084115d-02/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 6.35175936d+00, 2.92731770d-02/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 2.60305432d+00, 2.84188805d-02/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 9.40204530d+00, 2.76313008d-02/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 2.60305432d+00, 2.84559161d-02/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 6.35175936d+00, 2.91250348d-02/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.26210537d+01, 3.06639448d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 2.18436906d+01, 3.42791859d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 3.55969037d+01, 3.39993114d-02/
      data fpp( 2,10,1),fpp( 2,10,2)/ 5.88745752d+01, 3.77435683d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.21693530d+02, 1.05664152d-02/
      data fpp( 2,12,1),fpp( 2,12,2)/ 2.61539473d+02,-1.46249229d-01/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.61539473d+02,-1.44065520d-01/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.21693530d+02, 1.83157707d-03/
      data fpp( 2,15,1),fpp( 2,15,2)/ 5.88745752d+01, 7.04992115d-02/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 6.22631221d+00, 2.68329864d-02/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 1.88422278d-02, 2.61440272d-02/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.21069010d+00, 2.54909048d-02/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-3.25071586d+01, 2.47123536d-02/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-3.21069010d+00, 2.55396808d-02/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 1.88422278d-02, 2.59489232d-02/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 6.22631221d+00, 2.75646263d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 1.48470830d+01, 2.94525716d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 2.83108371d+01, 3.16850873d-02/
      data fpp( 3,10,1),fpp( 3,10,2)/ 3.93889868d+01, 2.38670793d-02/
      data fpp( 3,11,1),fpp( 3,11,2)/ 1.16822646d+02, 1.42659553d-03/
      data fpp( 3,12,1),fpp( 3,12,2)/ 1.78111845d+02,-1.18913461d-01/
      data fpp( 3,13,1),fpp( 3,13,2)/ 1.78111845d+02,-1.17312914d-01/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.16822646d+02,-4.97559614d-03/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.93889868d+01, 4.78752981d-02/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 2.73697501d-01, 2.38813039d-02/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-5.22712827d+00, 2.32973921d-02/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-8.96029390d+00, 2.27291275d-02/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 9.36265890d+01, 2.21060979d-02/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-8.96029390d+00, 2.27664810d-02/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-5.22712827d+00, 2.31479781d-02/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 2.73697501d-01, 2.44416065d-02/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 8.76797739d+00, 2.59055957d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.09597479d+01, 2.70360105d-02/
      data fpp( 4,10,1),fpp( 4,10,2)/ 4.03694778d+01, 1.87703621d-02/
      data fpp( 4,11,1),fpp( 4,11,2)/ 8.24158868d+01,-7.73745886d-03/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.43813148d+02,-9.70805266d-02/
      data fpp( 4,13,1),fpp( 4,13,2)/ 1.43813148d+02,-9.63596907d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/ 8.24158868d+01,-1.06208027d-02/
      data fpp( 4,15,1),fpp( 4,15,2)/ 4.03694778d+01, 2.95829013d-02/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-4.32110221d+00, 1.75494468d-02/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-9.71032913d+00, 1.83411064d-02/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-1.25481343d+01, 3.30461278d-02/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-1.16999197d+02,-9.45617391d-04/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-1.25481343d+01, 3.33763418d-02/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-9.71032913d+00, 1.70202502d-02/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-4.32110221d+00, 2.25026575d-02/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 3.48100743d+00, 2.26291199d-02/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.44501711d+01, 2.24008630d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/ 3.31331021d+01, 1.67674280d-02/
      data fpp( 5,11,1),fpp( 5,11,2)/ 5.15138069d+01,-1.36305752d-02/
      data fpp( 5,12,1),fpp( 5,12,2)/ 1.19035563d+02,-8.03851273d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 1.19035563d+02,-8.03493304d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/ 5.15138069d+01,-1.37737627d-02/
      data fpp( 5,15,1),fpp( 5,15,2)/ 3.31331021d+01, 1.73043814d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-1.10235421d+01, 1.63575553d-02/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.58054485d+01, 1.60648894d-02/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.86254501d+01, 1.58628872d-02/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 2.71342972d+01, 1.53435620d-02/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.86254501d+01, 1.58828649d-02/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.58054485d+01, 1.59849784d-02/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.10235421d+01, 1.66572217d-02/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-4.02701098d+00, 1.71061350d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 6.01961265d+00, 1.67382383d-02/
      data fpp( 6,10,1),fpp( 6,10,2)/ 1.93159547d+01, 9.30091174d-03/
      data fpp( 6,11,1),fpp( 6,11,2)/ 4.14006359d+01,-7.50188526d-03/
      data fpp( 6,12,1),fpp( 6,12,2)/ 4.10867365d+01,-6.13133707d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/ 4.10867365d+01,-6.09189453d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/ 4.14006359d+01,-9.07958707d-03/
      data fpp( 6,15,1),fpp( 6,15,2)/ 1.93159547d+01, 1.52172935d-02/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-1.39847294d+01, 1.24796153d-02/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.82678770d+01, 1.22907694d-02/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-2.08000652d+01, 1.22173070d-02/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-3.41379917d+01, 1.17400027d-02/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-2.08000652d+01, 1.22226822d-02/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.82678770d+01, 1.22692683d-02/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-1.39847294d+01, 1.25602444d-02/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-7.47296350d+00, 1.26697541d-02/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.52137829d+00, 1.21207392d-02/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.36530789d+01, 6.10728897d-03/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.67836496d+01,-5.16989513d-03/
      data fpp( 7,12,1),fpp( 7,12,2)/ 4.40674908d+01,-4.58477085d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 4.40674908d+01,-4.54819270d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/ 2.67836496d+01,-6.63302085d-03/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.36530789d+01, 1.15940104d-02/
      data fpp( 8, 1,1),fpp( 8, 1,2)/-1.38875405d+01, 9.32091752d-03/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.77730437d+01, 9.21816495d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-2.01242891d+01, 9.24642266d-03/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-1.76323303d+01, 8.81614440d-03/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-2.01242891d+01, 9.24899974d-03/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.77730437d+01, 9.20785663d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.38875405d+01, 9.35957374d-03/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-8.08113500d+00, 9.21384841d-03/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.05125792d-01, 8.74503263d-03/
      data fpp( 8,10,1),fpp( 8,10,2)/ 9.92172955d+00, 4.10602107d-03/
      data fpp( 8,11,1),fpp( 8,11,2)/ 2.14147658d+01,-3.14911691d-03/
      data fpp( 8,12,1),fpp( 8,12,2)/ 2.90933002d+01,-3.41695534d-02/
      data fpp( 8,13,1),fpp( 8,13,2)/ 2.90933002d+01,-3.38267813d-02/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.14147658d+01,-4.52020534d-03/
      data fpp( 8,15,1),fpp( 8,15,2)/ 9.92172955d+00, 9.24760267d-03/
      data fpp( 9, 1,1),fpp( 9, 1,2)/-1.22151087d+01, 6.81695831d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-1.57899483d+01, 6.78608338d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.78027783d+01, 6.83870817d-03/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-1.93826870d+01, 6.59908394d-03/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.78027783d+01, 6.84495606d-03/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-1.57899483d+01, 6.76109181d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-1.22151087d+01, 6.91067669d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-7.30249650d+00, 6.63620145d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-6.50875119d-01, 6.26451753d-03/
      data fpp( 9,10,1),fpp( 9,10,2)/ 7.56000288d+00, 2.80572842d-03/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.64072873d+01,-1.82743122d-03/
      data fpp( 9,12,1),fpp( 9,12,2)/ 2.43593083d+01,-2.51960035d-02/
      data fpp( 9,13,1),fpp( 9,13,2)/ 2.43593083d+01,-2.48982738d-02/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.64072873d+01,-3.01835035d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 7.56000288d+00, 7.27167518d-03/
      data fpp(10, 1,1),fpp(10, 1,2)/-1.01520246d+01, 4.95068431d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.29671630d+01, 4.94863138d-03/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.47145978d+01, 4.95479017d-03/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.49369217d+01, 4.93220795d-03/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.47145978d+01, 4.95637804d-03/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.29671630d+01, 4.94227988d-03/
      data fpp(10, 7,1),fpp(10, 7,2)/-1.01520246d+01, 4.97450242d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/-6.05887902d+00, 4.79971044d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/-7.41373733d-01, 4.44665583d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 5.73825893d+00, 1.89366626d-03/
      data fpp(10,11,1),fpp(10,11,2)/ 1.27060851d+01,-9.21320849d-04/
      data fpp(10,12,1),fpp(10,12,2)/ 1.83694667d+01,-1.84283829d-02/
      data fpp(10,13,1),fpp(10,13,2)/ 1.83694667d+01,-1.81741910d-02/
      data fpp(10,14,1),fpp(10,14,2)/ 1.27060851d+01,-1.93808829d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 5.73825893d+00, 5.70654414d-03/
      data fpp(11, 1,1),fpp(11, 1,2)/-7.92679270d+00, 3.50225582d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/-1.01913998d+01, 3.54548836d-03/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.15388306d+01, 3.61579073d-03/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.20696264d+01, 3.59134871d-03/
      data fpp(11, 5,1),fpp(11, 5,2)/-1.15388306d+01, 3.61881444d-03/
      data fpp(11, 6,1),fpp(11, 6,2)/-1.01913998d+01, 3.53339354d-03/
      data fpp(11, 7,1),fpp(11, 7,2)/-7.92679270d+00, 3.54761141d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/-4.76198743d+00, 3.45616084d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/-7.33629950d-01, 3.14774523d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 4.13696141d+00, 1.23285823d-03/
      data fpp(11,11,1),fpp(11,11,2)/ 9.41837238d+00,-3.39178149d-04/
      data fpp(11,12,1),fpp(11,12,2)/ 1.39128249d+01,-1.33761456d-02/
      data fpp(11,13,1),fpp(11,13,2)/ 1.39128249d+01,-1.31619740d-02/
      data fpp(11,14,1),fpp(11,14,2)/ 9.41837238d+00,-1.19586456d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 4.13696141d+00, 4.44543228d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/-5.34267456d+00, 2.00361634d-03/
      data fpp(12, 2,1),fpp(12, 2,2)/-6.78389216d+00, 2.10276732d-03/
      data fpp(12, 3,1),fpp(12, 3,2)/-7.59416619d+00, 2.24531438d-03/
      data fpp(12, 4,1),fpp(12, 4,2)/-7.84329763d+00, 2.23597514d-03/
      data fpp(12, 5,1),fpp(12, 5,2)/-7.59416619d+00, 2.25078505d-03/
      data fpp(12, 6,1),fpp(12, 6,2)/-6.78389216d+00, 2.08088467d-03/
      data fpp(12, 7,1),fpp(12, 7,2)/-5.34267456d+00, 2.08567627d-03/
      data fpp(12, 8,1),fpp(12, 8,2)/-3.32078923d+00, 2.05641027d-03/
      data fpp(12, 9,1),fpp(12, 9,2)/-7.60317678d-01, 1.86868266d-03/
      data fpp(12,10,1),fpp(12,10,2)/ 2.28462268d+00, 6.08859075d-04/
      data fpp(12,11,1),fpp(12,11,2)/ 5.56803535d+00, 1.35881035d-04/
      data fpp(12,12,1),fpp(12,12,2)/ 8.34427251d+00,-8.17238321d-03/
      data fpp(12,13,1),fpp(12,13,2)/ 8.34427251d+00,-8.01079088d-03/
      data fpp(12,14,1),fpp(12,14,2)/ 5.56803535d+00,-5.10488321d-04/
      data fpp(12,15,1),fpp(12,15,2)/ 2.28462268d+00, 3.03274416d-03/
      data fpp(13, 1,1),fpp(13, 1,2)/-2.76336577d+00, 7.83105362d-04/
      data fpp(13, 2,1),fpp(13, 2,2)/-3.37670522d+00, 8.73789275d-04/
      data fpp(13, 3,1),fpp(13, 3,2)/-3.75136984d+00, 1.00173754d-03/
      data fpp(13, 4,1),fpp(13, 4,2)/-3.86767174d+00, 9.99260575d-04/
      data fpp(13, 5,1),fpp(13, 5,2)/-3.75136984d+00, 1.00122016d-03/
      data fpp(13, 6,1),fpp(13, 6,2)/-3.37670522d+00, 8.75858776d-04/
      data fpp(13, 7,1),fpp(13, 7,2)/-2.76336577d+00, 7.75344732d-04/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.84428201d+00, 7.62762294d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/-6.54805460d-01, 7.33606092d-04/
      data fpp(13,10,1),fpp(13,10,2)/ 7.19030564d-01, 1.42813337d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 2.13126344d+00, 2.55140560d-04/
      data fpp(13,12,1),fpp(13,12,2)/ 3.23863302d+00,-3.38337558d-03/
      data fpp(13,13,1),fpp(13,13,2)/ 3.23863302d+00,-3.29744355d-03/
      data fpp(13,14,1),fpp(13,14,2)/ 2.13126344d+00,-8.85875577d-05/
      data fpp(13,15,1),fpp(13,15,2)/ 7.19030564d-01, 1.43179378d-03/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.38786234d+00, 2.84897083d-04/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.66928696d+00, 3.50205834d-04/
      data fpp(14, 3,1),fpp(14, 3,2)/-1.83235446d+00, 4.74279581d-04/
      data fpp(14, 4,1),fpp(14, 4,2)/-1.91001540d+00, 3.92675843d-04/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.83235446d+00, 4.75017047d-04/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.66928696d+00, 3.47255968d-04/
      data fpp(14, 7,1),fpp(14, 7,2)/-1.38786234d+00, 2.95959083d-04/
      data fpp(14, 8,1),fpp(14, 8,2)/-1.01408271d+00, 2.08907702d-04/
      data fpp(14, 9,1),fpp(14, 9,2)/-5.08460483d-01, 1.88410108d-04/
      data fpp(14,10,1),fpp(14,10,2)/ 1.19255059d-01,-6.25481349d-05/
      data fpp(14,11,1),fpp(14,11,2)/ 8.10910889d-01,-1.18217569d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 1.57319540d+00,-1.02458159d-03/
      data fpp(14,13,1),fpp(14,13,2)/ 1.57319540d+00,-1.00714644d-03/
      data fpp(14,14,1),fpp(14,14,2)/ 8.10910889d-01,-1.87958159d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 1.19255059d-01, 1.98979080d-04/
      data fpp(15, 1,1),fpp(15, 1,2)/-7.09184868d-01, 1.19976652d-04/
      data fpp(15, 2,1),fpp(15, 2,2)/-8.18146959d-01, 1.50046697d-04/
      data fpp(15, 3,1),fpp(15, 3,2)/-8.95212310d-01, 1.79836561d-04/
      data fpp(15, 4,1),fpp(15, 4,2)/-9.00266662d-01, 2.10607060d-04/
      data fpp(15, 5,1),fpp(15, 5,2)/-8.95212310d-01, 1.77735198d-04/
      data fpp(15, 6,1),fpp(15, 6,2)/-8.18146959d-01, 1.58452148d-04/
      data fpp(15, 7,1),fpp(15, 7,2)/-7.09184868d-01, 8.84562105d-05/
      data fpp(15, 8,1),fpp(15, 8,2)/-5.55387127d-01, 2.77230103d-05/
      data fpp(15, 9,1),fpp(15, 9,2)/-3.59352610d-01,-1.93482515d-05/
      data fpp(15,10,1),fpp(15,10,2)/-1.16050801d-01,-1.90330004d-04/
      data fpp(15,11,1),fpp(15,11,2)/ 2.17093005d-01,-1.79331732d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 5.24585385d-01,-1.12343070d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 5.24585385d-01,-1.13304926d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 2.17093005d-01,-1.75484307d-04/
      data fpp(15,15,1),fpp(15,15,2)/-1.16050801d-01,-2.04757847d-04/
      data fpp(16, 1,1),fpp(16, 1,2)/-3.59398190d-01, 2.13233100d-05/
      data fpp(16, 2,1),fpp(16, 2,2)/-4.10125207d-01, 4.73533800d-05/
      data fpp(16, 3,1),fpp(16, 3,2)/-4.42796298d-01, 8.92631698d-05/
      data fpp(16, 4,1),fpp(16, 4,2)/-4.64917953d-01, 7.55939406d-05/
      data fpp(16, 5,1),fpp(16, 5,2)/-4.42796298d-01, 8.83610679d-05/
      data fpp(16, 6,1),fpp(16, 6,2)/-4.10125207d-01, 5.09617880d-05/
      data fpp(16, 7,1),fpp(16, 7,2)/-3.59398190d-01, 7.79178028d-06/
      data fpp(16, 8,1),fpp(16, 8,2)/-2.92368776d-01,-2.21289091d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/-2.14129078d-01,-9.92761439d-05/
      data fpp(16,10,1),fpp(16,10,2)/-1.11051856d-01,-1.20766515d-04/
      data fpp(16,11,1),fpp(16,11,2)/-2.32829091d-02,-1.37657796d-04/
      data fpp(16,12,1),fpp(16,12,2)/ 7.24630606d-02, 1.31397697d-04/
      data fpp(16,13,1),fpp(16,13,2)/ 7.24630606d-02, 1.24635806d-04/
      data fpp(16,14,1),fpp(16,14,2)/-2.32829091d-02,-1.10610230d-04/
      data fpp(16,15,1),fpp(16,15,2)/-1.11051856d-01,-2.22194885d-04/
      data fpp(17, 1,1),fpp(17, 1,2)/-2.05222374d-01,-2.93240410d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/-2.29352212d-01,-1.35191805d-06/
      data fpp(17, 3,1),fpp(17, 3,2)/-2.37602497d-01, 3.47317132d-05/
      data fpp(17, 4,1),fpp(17, 4,2)/-2.40061527d-01, 4.24250654d-05/
      data fpp(17, 5,1),fpp(17, 5,2)/-2.37602497d-01, 3.55680252d-05/
      data fpp(17, 6,1),fpp(17, 6,2)/-2.29352212d-01,-4.69716614d-06/
      data fpp(17, 7,1),fpp(17, 7,2)/-2.05222374d-01,-1.67793606d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/-1.71137767d-01,-4.81853914d-05/
      data fpp(17, 9,1),fpp(17, 9,2)/-1.28131079d-01,-9.04790736d-05/
      data fpp(17,10,1),fpp(17,10,2)/-8.77417735d-02,-6.98983141d-05/
      data fpp(17,11,1),fpp(17,11,2)/-5.19613685d-02,-4.99276701d-05/
      data fpp(17,12,1),fpp(17,12,2)/-4.64376279d-02, 1.49608994d-04/
      data fpp(17,13,1),fpp(17,13,2)/-4.64376279d-02, 1.46136852d-04/
      data fpp(17,14,1),fpp(17,14,2)/-5.19613685d-02,-3.60391006d-05/
      data fpp(17,15,1),fpp(17,15,2)/-8.77417735d-02,-1.21980450d-04/
      data fpp(18, 1,1),fpp(18, 1,2)/-5.66337838d-02,-3.12395281d-06/
      data fpp(18, 2,1),fpp(18, 2,2)/-6.08807609d-02,-3.75209439d-06/
      data fpp(18, 3,1),fpp(18, 3,2)/-6.77943584d-02,-4.18676696d-05/
      data fpp(18, 4,1),fpp(18, 4,2)/-6.73564416d-02, 5.12227729d-05/
      data fpp(18, 5,1),fpp(18, 5,2)/-6.77943584d-02,-4.30234221d-05/
      data fpp(18, 6,1),fpp(18, 6,2)/-6.08807609d-02, 8.70915519d-07/
      data fpp(18, 7,1),fpp(18, 7,2)/-5.66337838d-02,-2.04602400d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/-5.44023111d-02,-3.90299557d-05/
      data fpp(18, 9,1),fpp(18, 9,2)/-5.45422247d-02,-6.34199374d-05/
      data fpp(18,10,1),fpp(18,10,2)/-5.32487513d-02,-7.29029491d-06/
      data fpp(18,11,1),fpp(18,11,2)/-6.04744399d-02,-2.74188830d-05/
      data fpp(18,12,1),fpp(18,12,2)/-5.29186467d-02, 1.16965827d-04/
      data fpp(18,13,1),fpp(18,13,2)/-5.29186467d-02, 1.12811961d-04/
      data fpp(18,14,1),fpp(18,14,2)/-6.04744399d-02,-1.08034173d-05/
      data fpp(18,15,1),fpp(18,15,2)/-5.32487513d-02,-6.95982913d-05/
      data fpp(19, 1,1),fpp(19, 1,2)/-2.42424909d-02,-1.06727931d-05/
      data fpp(19, 2,1),fpp(19, 2,2)/-2.51247445d-02,-8.65441378d-06/
      data fpp(19, 3,1),fpp(19, 3,2)/-2.52200689d-02,-1.47095518d-05/
      data fpp(19, 4,1),fpp(19, 4,2)/-2.45127064d-02, 7.49262083d-06/
      data fpp(19, 5,1),fpp(19, 5,2)/-2.52200689d-02,-1.52609315d-05/
      data fpp(19, 6,1),fpp(19, 6,2)/-2.51247445d-02,-6.44889466d-06/
      data fpp(19, 7,1),fpp(19, 7,2)/-2.42424909d-02,-1.89434898d-05/
      data fpp(19, 8,1),fpp(19, 8,2)/-2.52529888d-02,-3.77771461d-05/
      data fpp(19, 9,1),fpp(19, 9,2)/-2.57000222d-02,-9.94792572d-06/
      data fpp(19,10,1),fpp(19,10,2)/-2.92632213d-02,-4.24311510d-05/
      data fpp(19,11,1),fpp(19,11,2)/-2.41408718d-02, 5.96725297d-05/
      data fpp(19,12,1),fpp(19,12,2)/-3.58877854d-02, 4.37410321d-05/
      data fpp(19,13,1),fpp(19,13,2)/-3.58877854d-02, 4.89406388d-05/
      data fpp(19,14,1),fpp(19,14,2)/-2.41408718d-02, 3.88741032d-05/
      data fpp(19,15,1),fpp(19,15,2)/-2.92632213d-02, 3.55629484d-05/
      data fpp(20, 1,1),fpp(20, 1,2)/-1.43962524d-02,-2.71188541d-05/
      data fpp(20, 2,1),fpp(20, 2,2)/-1.26202611d-02,-1.57622918d-05/
      data fpp(20, 3,1),fpp(20, 3,2)/-1.13253661d-02, 3.01680213d-05/
      data fpp(20, 4,1),fpp(20, 4,2)/-1.45927329d-02,-4.49097936d-05/
      data fpp(20, 5,1),fpp(20, 5,2)/-1.13253661d-02, 2.94711529d-05/
      data fpp(20, 6,1),fpp(20, 6,2)/-1.26202611d-02,-1.29748179d-05/
      data fpp(20, 7,1),fpp(20, 7,2)/-1.43962524d-02,-3.75718813d-05/
      data fpp(20, 8,1),fpp(20, 8,2)/-1.25857338d-02,-1.67376570d-05/
      data fpp(20, 9,1),fpp(20, 9,2)/-1.06576863d-02,-1.54774908d-05/
      data fpp(20,10,1),fpp(20,10,2)/-9.69836350d-03, 1.86476200d-05/
      data fpp(20,11,1),fpp(20,11,2)/-1.69620730d-02, 8.87010582d-07/
      data fpp(20,12,1),fpp(20,12,2)/-1.35302116d-02, 3.78043376d-05/
      data fpp(20,13,1),fpp(20,13,2)/-1.35302116d-02, 3.61434818d-05/
      data fpp(20,14,1),fpp(20,14,2)/-1.69620730d-02, 7.53043376d-06/
      data fpp(20,15,1),fpp(20,15,2)/-9.69836350d-03,-6.26521688d-06/
      data fpp(21, 1,1),fpp(21, 1,2)/-8.17249930d-03,-1.95606569d-05/
      data fpp(21, 2,1),fpp(21, 2,2)/-8.39421112d-03,-1.08786862d-05/
      data fpp(21, 3,1),fpp(21, 3,2)/-7.47846684d-03, 3.07540172d-06/
      data fpp(21, 4,1),fpp(21, 4,2)/-7.11636203d-03,-1.42292066d-06/
      data fpp(21, 5,1),fpp(21, 5,2)/-7.47846684d-03, 2.61628092d-06/
      data fpp(21, 6,1),fpp(21, 6,2)/-8.39421112d-03,-9.04220303d-06/
      data fpp(21, 7,1),fpp(21, 7,2)/-8.17249930d-03,-2.64474688d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/-8.40407606d-03,-5.16792174d-06/
      data fpp(21, 9,1),fpp(21, 9,2)/-9.66923248d-03,-1.28808442d-05/
      data fpp(21,10,1),fpp(21,10,2)/-9.94332472d-03,-3.30870127d-06/
      data fpp(21,11,1),fpp(21,11,2)/-1.00108363d-02, 2.61156493d-05/
      data fpp(21,12,1),fpp(21,12,2)/-1.19913681d-02, 1.88461039d-05/
      data fpp(21,13,1),fpp(21,13,2)/-1.19913681d-02, 2.04038636d-05/
      data fpp(21,14,1),fpp(21,14,2)/-1.00108363d-02, 1.98846104d-05/
      data fpp(21,15,1),fpp(21,15,2)/-9.94332472d-03, 2.00576948d-05/
      data fpp(22, 1,1),fpp(22, 1,2)/-9.13750349d-04, 1.64210941d-05/
      data fpp(22, 2,1),fpp(22, 2,2)/-1.80289444d-03, 7.15781180d-06/
      data fpp(22, 3,1),fpp(22, 3,2)/-7.60766582d-04,-4.50523413d-05/
      data fpp(22, 4,1),fpp(22, 4,2)/ 1.05818101d-03, 5.30515534d-05/
      data fpp(22, 5,1),fpp(22, 5,2)/-7.60766582d-04,-4.71538723d-05/
      data fpp(22, 6,1),fpp(22, 6,2)/-1.80289444d-03, 1.55639357d-05/
      data fpp(22, 7,1),fpp(22, 7,2)/-9.13750349d-04,-1.51018704d-05/
      data fpp(22, 8,1),fpp(22, 8,2)/-1.79796197d-03,-1.51564543d-05/
      data fpp(22, 9,1),fpp(22, 9,2)/-4.66538376d-03, 1.57276874d-05/
      data fpp(22,10,1),fpp(22,10,2)/-4.52833764d-03,-4.77542952d-05/
      data fpp(22,11,1),fpp(22,11,2)/-2.99458185d-03, 5.52894933d-05/
      data fpp(22,12,1),fpp(22,12,2)/-4.50431595d-03, 6.59632202d-06/
      data fpp(22,13,1),fpp(22,13,2)/-4.50431595d-03, 1.25662873d-05/
      data fpp(22,14,1),fpp(22,14,2)/-2.99458185d-03, 3.14096322d-05/
      data fpp(22,15,1),fpp(22,15,2)/-4.52833764d-03, 4.17951839d-05/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 5.28451126d-02, 4.37720040d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/-7.55508060d-02, 1.91000339d-01/
      data fpppp( 1, 5),fpppp( 1, 6)/-7.65110909d-02, 4.76131436d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 3.84408392d-02,-1.09767790d-02/
      data fpppp( 1, 9),fpppp( 1,10)/ 2.83548991d-01,-2.05952936d-01/
      data fpppp( 1,11),fpppp( 1,12)/ 2.51302426d+00,-1.83995377d+00/
      data fpppp( 1,13),fpppp( 1,14)/-1.61446167d+00, 1.61105585d+00/
      data fpppp( 1,15) /             3.17642859d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/-3.53418368d-02,-1.42407851d-02/
      data fpppp( 2, 3),fpppp( 2, 4)/ 2.43540332d-01,-3.27058779d-01/
      data fpppp( 2, 5),fpppp( 2, 6)/ 2.48815866d-01,-3.53429245d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 4.37911860d-02, 3.73787380d-02/
      data fpppp( 2, 9),fpppp( 2,10)/ 7.85284343d-02, 2.19975029d-01/
      data fpppp( 2,11),fpppp( 2,12)/ 1.41404843d+00,-1.25454943d+00/
      data fpppp( 2,13),fpppp( 2,14)/-1.13875422d+00, 9.50867594d-01/
      data fpppp( 2,15) /             1.95690316d+00 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 2.50638533d-01, 1.68958539d-01/
      data fpppp( 3, 3),fpppp( 3, 4)/-7.47796432d-01, 1.25821102d+00/
      data fpppp( 3, 5),fpppp( 3, 6)/-7.69471430d-01, 2.55658532d-01/
      data fpppp( 3, 7),fpppp( 3, 8)/-7.44864382d-02, 1.87085270d-01/
      data fpppp( 3, 9),fpppp( 3,10)/-3.83275645d-01, 1.20288104d+00/
      data fpppp( 3,11),fpppp( 3,12)/-4.46917960d-01,-3.83876816d-01/
      data fpppp( 3,13),fpppp( 3,14)/-4.63586540d-01,-1.28079060d-01/
      data fpppp( 3,15) /             7.23517020d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-7.72966714d-01,-4.92084534d-01/
      data fpppp( 4, 3),fpppp( 4, 4)/ 2.84736446d+00,-4.51817039d+00/
      data fpppp( 4, 5),fpppp( 4, 6)/ 2.91489116d+00,-7.62191347d-01/
      data fpppp( 4, 7),fpppp( 4, 8)/ 2.39933836d-01,-1.79367503d-02/
      data fpppp( 4, 9),fpppp( 4,10)/ 5.36626051d-02, 2.36363887d-01/
      data fpppp( 4,11),fpppp( 4,12)/ 3.59082600d-01,-5.11643155d-01/
      data fpppp( 4,13),fpppp( 4,14)/-4.86529671d-01, 2.58628663d-01/
      data fpppp( 4,15) /             6.13066148d-01 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 8.13123543d-01, 5.33589155d-01/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.79439486d+00, 4.54719482d+00/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.86025686d+00, 7.97037163d-01/
      data fpppp( 5, 7),fpppp( 5, 8)/-1.74806487d-01, 4.69617467d-02/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.76982743d-01,-2.92066681d-01/
      data fpppp( 5,11),fpppp( 5,12)/ 9.73150406d-01,-6.52071856d-01/
      data fpppp( 5,13),fpppp( 5,14)/-5.56012324d-01, 5.88912278d-01/
      data fpppp( 5,15) /             1.14882630d+00 /
      data fpppp( 6, 1),fpppp( 6, 2)/-3.32113869d-01,-2.09272007d-01/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.28691618d+00,-2.02360777d+00/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.31634523d+00,-3.26988208d-01/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.09321883d-01, 2.25781617d-02/
      data fpppp( 6, 9),fpppp( 6,10)/-1.66289796d-02, 2.38920864d-01/
      data fpppp( 6,11),fpppp( 6,12)/-4.11754135d-01, 6.41808460d-02/
      data fpppp( 6,13),fpppp( 6,14)/ 2.27515109d-02,-2.46036795d-01/
      data fpppp( 6,15) /            -3.82519161d-01 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 1.12424654d-01, 7.88134277d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/-3.22620802d-01, 5.63325479d-01/
      data fpppp( 7, 5),fpppp( 7, 6)/-3.30129929d-01, 1.08849938d-01/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.12260893d-04, 2.57161995d-02/
      data fpppp( 7, 9),fpppp( 7,10)/ 4.63020195d-02,-2.26827463d-02/
      data fpppp( 7,11),fpppp( 7,12)/ 1.04361166d-01,-1.45565685d-01/
      data fpppp( 7,13),fpppp( 7,14)/-1.33998765d-01, 5.80934873d-02/
      data fpppp( 7,15) /             1.50821050d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/-9.69840729d-03,-7.54384499d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/ 1.04771409d-01,-1.27738995d-01/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.07149513d-01,-1.02668007d-02/
      data fpppp( 8, 7),fpppp( 8, 8)/ 2.59731535d-02, 2.16283241d-02/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.76897735d-02, 3.06633496d-02/
      data fpppp( 8,11),fpppp( 8,12)/-5.23723183d-02,-5.00441838d-02/
      data fpppp( 8,13),fpppp( 8,14)/-5.40373229d-02,-3.63997618d-02/
      data fpppp( 8,15) /            -2.92337371d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 2.45797567d-02, 2.02825470d-02/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.19893634d-02, 5.36501784d-02/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.30223015d-02, 2.44142991d-02/
      data fpppp( 9, 7),fpppp( 9, 8)/ 9.08568610d-03, 1.95093127d-02/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.72176123d-02, 5.17563565d-03/
      data fpppp( 9,11),fpppp( 9,12)/ 2.64229234d-04,-5.99483577d-02/
      data fpppp( 9,13),fpppp( 9,14)/-5.88476710d-02,-4.13851776d-03/
      data fpppp( 9,15) /             2.16859369d-02 /
      data fpppp(10, 1),fpppp(10, 2)/ 7.03216228d-03, 8.81863271d-03/
      data fpppp(10, 3),fpppp(10, 4)/ 2.17555182d-02,-4.33404908d-03/
      data fpppp(10, 5),fpppp(10, 6)/ 2.22595421d-02, 6.80253717d-03/
      data fpppp(10, 7),fpppp(10, 8)/ 1.45925205d-02, 1.15078183d-02/
      data fpppp(10, 9),fpppp(10,10)/ 1.28377860d-02, 6.86868032d-03/
      data fpppp(10,11),fpppp(10,12)/-1.10208974d-02,-4.10517630d-02/
      data fpppp(10,13),fpppp(10,14)/-4.12357111d-02,-1.02851049d-02/
      data fpppp(10,15) /             4.10945850d-03 /
      data fpppp(11, 1),fpppp(11, 2)/ 9.96842054d-03, 9.58926936d-03/
      data fpppp(11, 3),fpppp(11, 4)/ 6.70507543d-03, 1.25885311d-02/
      data fpppp(11, 5),fpppp(11, 6)/ 6.63629536d-03, 9.86438963d-03/
      data fpppp(11, 7),fpppp(11, 8)/ 8.93671954d-03, 8.40062546d-03/
      data fpppp(11, 9),fpppp(11,10)/ 9.27391086d-03, 5.03776398d-03/
      data fpppp(11,11),fpppp(11,12)/-4.77579027d-03,-3.31521092d-02/
      data fpppp(11,13),fpppp(11,14)/-3.29893531d-02,-5.42681456d-03/
      data fpppp(11,15) /             7.47910506d-03 /
      data fpppp(12, 1),fpppp(12, 2)/ 7.01737387d-03, 6.28957929d-03/
      data fpppp(12, 3),fpppp(12, 4)/ 5.68092342d-03, 4.65528204d-03/
      data fpppp(12, 5),fpppp(12, 6)/ 5.59372160d-03, 6.63838656d-03/
      data fpppp(12, 7),fpppp(12, 8)/ 5.70934659d-03, 5.36429122d-03/
      data fpppp(12, 9),fpppp(12,10)/ 5.14866133d-03, 3.10919226d-03/
      data fpppp(12,11),fpppp(12,12)/-3.27709191d-03,-2.04313555d-02/
      data fpppp(12,13),fpppp(12,14)/-2.03545023d-02,-3.58450435d-03/
      data fpppp(12,15) /             4.26198888d-03 /
      data fpppp(13, 1),fpppp(13, 2)/ 2.21666695d-03, 2.33315336d-03/
      data fpppp(13, 3),fpppp(13, 4)/ 2.77120938d-03, 2.08377179d-03/
      data fpppp(13, 5),fpppp(13, 6)/ 2.84993211d-03, 2.01826241d-03/
      data fpppp(13, 7),fpppp(13, 8)/ 3.39750801d-03, 2.73636436d-03/
      data fpppp(13, 9),fpppp(13,10)/ 1.88060224d-03, 8.02794820d-04/
      data fpppp(13,11),fpppp(13,12)/-2.78797046d-03,-7.94271046d-03/
      data fpppp(13,13),fpppp(13,14)/-7.99897090d-03,-2.56292870d-03/
      data fpppp(13,15) /            -4.11117853d-05 /
      data fpppp(14, 1),fpppp(14, 2)/ 1.45122309d-03, 1.30727772d-03/
      data fpppp(14, 3),fpppp(14, 4)/ 4.21092450d-04, 2.13274685d-03/
      data fpppp(14, 5),fpppp(14, 6)/ 3.67232413d-04, 1.52271787d-03/
      data fpppp(14, 7),fpppp(14, 8)/ 6.43322544d-04, 1.44529252d-03/
      data fpppp(14, 9),fpppp(14,10)/ 1.48606380d-03,-6.39491407d-05/
      data fpppp(14,11),fpppp(14,12)/ 2.60615005d-03,-6.12293030d-03/
      data fpppp(14,13),fpppp(14,14)/-5.80281938d-03, 1.32570637d-03/
      data fpppp(14,15) /             4.73771464d-03 /
      data fpppp(15, 1),fpppp(15, 2)/-2.26509547d-05, 1.99919036d-04/
      data fpppp(15, 3),fpppp(15, 4)/ 1.13677931d-03,-4.26376401d-04/
      data fpppp(15, 5),fpppp(15, 6)/ 1.17524860d-03, 4.60418530d-05/
      data fpppp(15, 7),fpppp(15, 8)/ 5.54388481d-04, 4.26543147d-04/
      data fpppp(15, 9),fpppp(15,10)/ 2.73645541d-04, 1.31491220d-03/
      data fpppp(15,11),fpppp(15,12)/-1.42774558d-04,-2.28289949d-03/
      data fpppp(15,13),fpppp(15,14)/-2.30468565d-03,-5.56299145d-05/
      data fpppp(15,15) /             9.88119787d-04 /
      data fpppp(16, 1),fpppp(16, 2)/ 2.29790887d-04, 2.32225819d-04/
      data fpppp(16, 3),fpppp(16, 4)/-7.53385707d-05, 7.02094663d-04/
      data fpppp(16, 5),fpppp(16, 6)/-7.84415540d-05, 2.44637752d-04/
      data fpppp(16, 7),fpppp(16, 8)/ 1.83246137d-04, 5.21435989d-07/
      data fpppp(16, 9),fpppp(16,10)/ 4.87285248d-04,-4.59411063d-04/
      data fpppp(16,11),fpppp(16,12)/ 4.31862561d-04,-7.89417836d-04/
      data fpppp(16,13),fpppp(16,14)/-7.20056864d-04, 1.54418673d-04/
      data fpppp(16,15) /             5.81003517d-04 /
      data fpppp(17, 1),fpppp(17, 2)/ 2.55631181d-04, 1.66890141d-04/
      data fpppp(17, 3),fpppp(17, 4)/ 2.95813966d-05, 6.22596193d-05/
      data fpppp(17, 5),fpppp(17, 6)/ 1.64637121d-05, 2.19360879d-04/
      data fpppp(17, 7),fpppp(17, 8)/ 5.88659140d-05, 1.42461596d-04/
      data fpppp(17, 9),fpppp(17,10)/-9.33874258d-05, 7.40451372d-05/
      data fpppp(17,11),fpppp(17,12)/-4.79327141d-04, 2.78635652d-05/
      data fpppp(17,13),fpppp(17,14)/-9.63934345d-06,-3.29315507d-04/
      data fpppp(17,15) /            -4.88498492d-04 /
      data fpppp(18, 1),fpppp(18, 2)/-1.16820246d-04,-4.67208141d-05/
      data fpppp(18, 3),fpppp(18, 4)/ 1.43706279d-04,-8.70134425d-05/
      data fpppp(18, 5),fpppp(18, 6)/ 1.51797469d-04,-7.90855741d-05/
      data fpppp(18, 7),fpppp(18, 8)/ 4.54760408d-06,-6.00351039d-05/
      data fpppp(18, 9),fpppp(18,10)/ 9.33096271d-05,-2.27200178d-04/
      data fpppp(18,11),fpppp(18,12)/ 3.04341359d-04,-1.03276346d-04/
      data fpppp(18,13),fpppp(18,14)/-6.90154382d-05, 1.67297727d-04/
      data fpppp(18,15) /             2.86713446d-04 /
      data fpppp(19, 1),fpppp(19, 2)/ 9.15623215d-06, 4.98025692d-06/
      data fpppp(19, 3),fpppp(19, 4)/ 1.81384907d-05,-2.93730071d-05/
      data fpppp(19, 5),fpppp(19, 6)/ 1.44700384d-05, 1.96540663d-05/
      data fpppp(19, 7),fpppp(19, 8)/-4.58705532d-05, 5.02630618d-05/
      data fpppp(19, 9),fpppp(19,10)/-1.21373830d-04, 2.48262322d-04/
      data fpppp(19,11),fpppp(19,12)/-3.50542542d-04, 1.41752056d-04/
      data fpppp(19,13),fpppp(19,14)/ 1.02422514d-04,-1.93224376d-04/
      data fpppp(19,15) /            -3.41680802d-04 /
      data fpppp(20, 1),fpppp(20, 2)/ 2.68119933d-05, 1.35664346d-05/
      data fpppp(20, 3),fpppp(20, 4)/-1.09943511d-04, 1.52471897d-04/
      data fpppp(20, 5),fpppp(20, 6)/-1.07860060d-04, 5.23263050d-06/
      data fpppp(20, 7),fpppp(20, 8)/ 5.80637588d-05,-2.22970659d-05/
      data fpppp(20, 9),fpppp(20,10)/ 3.81762338d-05,-1.88531348d-04/
      data fpppp(20,11),fpppp(20,12)/ 2.22567218d-04,-6.00032739d-05/
      data fpppp(20,13),fpppp(20,14)/-3.42296285d-05, 1.19472637d-04/
      data fpppp(20,15) /             1.98073333d-04 /
      data fpppp(21, 1),fpppp(21, 2)/ 2.76213438d-05, 1.27029064d-05/
      data fpppp(21, 3),fpppp(21, 4)/-1.01856039d-05,-5.17885888d-06/
      data fpppp(21, 5),fpppp(21, 6)/-1.25515374d-05, 2.21666403d-05/
      data fpppp(21, 7),fpppp(21, 8)/-7.86765825d-06,-1.78933215d-05/
      data fpppp(21, 9),fpppp(21,10)/ 1.74261644d-05, 7.65251468d-06/
      data fpppp(21,11),fpppp(21,12)/-3.56413835d-05, 2.01318046d-05/
      data fpppp(21,13),fpppp(21,14)/ 1.68412325d-05,-2.24790950d-05/
      data fpppp(21,15) /            -4.17060671d-05 /
      data fpppp(22, 1),fpppp(22, 2)/ 3.32439544d-05, 1.45392942d-05/
      data fpppp(22, 3),fpppp(22, 4)/ 2.44751862d-05,-6.58308549d-05/
      data fpppp(22, 5),fpppp(22, 6)/ 2.05745219d-05, 3.01419513d-05/
      data fpppp(22, 7),fpppp(22, 8)/-2.52660098d-05,-3.54792548d-05/
      data fpppp(22, 9),fpppp(22,10)/ 4.81904189d-05, 2.29856536d-05/
      data fpppp(22,11),fpppp(22,12)/-5.63304532d-05, 1.97267666d-05/
      data fpppp(22,13),fpppp(22,14)/ 1.42769497d-05,-3.45311856d-05/
      data fpppp(22,15) /            -5.87615999d-05 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /

      iprint=0

      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
      subroutine ch3po_1app_d(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     1(2a") surface
c     difference of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -4.51400000d+00 , -2.45100000d+00 /
      data f( 1, 3),f( 1, 4) / -1.08100000d+00 ,  0.00000000d+00 /
      data f( 1, 5),f( 1, 6) /  1.08100000d+00 ,  2.45100000d+00 /
      data f( 1, 7),f( 1, 8) /  4.51400000d+00 ,  7.97100000d+00 /
      data f( 1, 9),f( 1,10) /  1.41330000d+01 ,  2.47030000d+01 /
      data f( 1,11),f( 1,12) /  4.38600000d+01 ,  6.70360000d+01 /
      data f( 1,13),f( 1,14) /  6.70360000d+01 ,  4.38600000d+01 /
      data f( 1,15) /           2.47030000d+01 /
      data f( 2, 1),f( 2, 2) / -3.69700000d+00 , -1.98300000d+00 /
      data f( 2, 3),f( 2, 4) / -8.70000000d-01 ,  0.00000000d+00 /
      data f( 2, 5),f( 2, 6) /  8.70000000d-01 ,  1.98300000d+00 /
      data f( 2, 7),f( 2, 8) /  3.69700000d+00 ,  6.64600000d+00 /
      data f( 2, 9),f( 2,10) /  1.19700000d+01 ,  2.11820000d+01 /
      data f( 2,11),f( 2,12) /  3.63790000d+01 ,  5.33610000d+01 /
      data f( 2,13),f( 2,14) /  5.33610000d+01 ,  3.63790000d+01 /
      data f( 2,15) /           2.11820000d+01 /
      data f( 3, 1),f( 3, 2) / -3.01800000d+00 , -1.60000000d+00 /
      data f( 3, 3),f( 3, 4) / -6.98000000d-01 ,  0.00000000d+00 /
      data f( 3, 5),f( 3, 6) /  6.98000000d-01 ,  1.60000000d+00 /
      data f( 3, 7),f( 3, 8) /  3.01800000d+00 ,  5.53100000d+00 /
      data f( 3, 9),f( 3,10) /  1.01460000d+01 ,  1.82490000d+01 /
      data f( 3,11),f( 3,12) /  3.04980000d+01 ,  4.36000000d+01 /
      data f( 3,13),f( 3,14) /  4.36000000d+01 ,  3.04980000d+01 /
      data f( 3,15) /           1.82490000d+01 /
      data f( 4, 1),f( 4, 2) / -2.45800000d+00 , -1.28600000d+00 /
      data f( 4, 3),f( 4, 4) / -5.57000000d-01 ,  0.00000000d+00 /
      data f( 4, 5),f( 4, 6) /  5.57000000d-01 ,  1.28600000d+00 /
      data f( 4, 7),f( 4, 8) /  2.45800000d+00 ,  4.59800000d+00 /
      data f( 4, 9),f( 4,10) /  8.62000000d+00 ,  1.56980000d+01 /
      data f( 4,11),f( 4,12) /  2.60820000d+01 ,  3.64030000d+01 /
      data f( 4,13),f( 4,14) /  3.64030000d+01 ,  2.60820000d+01 /
      data f( 4,15) /           1.56980000d+01 /
      data f( 5, 1),f( 5, 2) / -1.99700000d+00 , -1.02800000d+00 /
      data f( 5, 3),f( 5, 4) / -4.43000000d-01 ,  0.00000000d+00 /
      data f( 5, 5),f( 5, 6) /  4.43000000d-01 ,  1.02800000d+00 /
      data f( 5, 7),f( 5, 8) /  1.99700000d+00 ,  3.82300000d+00 /
      data f( 5, 9),f( 5,10) /  7.33800000d+00 ,  1.35380000d+01 /
      data f( 5,11),f( 5,12) /  2.26330000d+01 ,  3.10660000d+01 /
      data f( 5,13),f( 5,14) /  3.10660000d+01 ,  2.26330000d+01 /
      data f( 5,15) /           1.35380000d+01 /
      data f( 6, 1),f( 6, 2) / -1.30500000d+00 , -6.48000000d-01 /
      data f( 6, 3),f( 6, 4) / -2.75000000d-01 ,  0.00000000d+00 /
      data f( 6, 5),f( 6, 6) /  2.75000000d-01 ,  6.48000000d-01 /
      data f( 6, 7),f( 6, 8) /  1.30500000d+00 ,  2.64200000d+00 /
      data f( 6, 9),f( 6,10) /  5.34100000d+00 ,  1.01760000d+01 /
      data f( 6,11),f( 6,12) /  1.71360000d+01 ,  2.39700000d+01 /
      data f( 6,13),f( 6,14) /  2.39700000d+01 ,  1.71360000d+01 /
      data f( 6,15) /           1.01760000d+01 /
      data f( 7, 1),f( 7, 2) / -8.41000000d-01 , -3.95000000d-01 /
      data f( 7, 3),f( 7, 4) / -1.66000000d-01 ,  0.00000000d+00 /
      data f( 7, 5),f( 7, 6) /  1.66000000d-01 ,  3.95000000d-01 /
      data f( 7, 7),f( 7, 8) /  8.41000000d-01 ,  1.82100000d+00 /
      data f( 7, 9),f( 7,10) /  3.88600000d+00 ,  7.63100000d+00 /
      data f( 7,11),f( 7,12) /  1.30040000d+01 ,  1.83430000d+01 /
      data f( 7,13),f( 7,14) /  1.83430000d+01 ,  1.30040000d+01 /
      data f( 7,15) /           7.63100000d+00 /
      data f( 8, 1),f( 8, 2) / -5.25000000d-01 , -2.29000000d-01 /
      data f( 8, 3),f( 8, 4) / -9.50000000d-02 ,  0.00000000d+00 /
      data f( 8, 5),f( 8, 6) /  9.50000000d-02 ,  2.29000000d-01 /
      data f( 8, 7),f( 8, 8) /  5.25000000d-01 ,  1.24200000d+00 /
      data f( 8, 9),f( 8,10) /  2.80300000d+00 ,  5.64400000d+00 /
      data f( 8,11),f( 8,12) /  9.72000000d+00 ,  1.38330000d+01 /
      data f( 8,13),f( 8,14) /  1.38330000d+01 ,  9.72000000d+00 /
      data f( 8,15) /           5.64400000d+00 /
      data f( 9, 1),f( 9, 2) / -3.11000000d-01 , -1.21000000d-01 /
      data f( 9, 3),f( 9, 4) / -5.00000000d-02 ,  0.00000000d+00 /
      data f( 9, 5),f( 9, 6) /  5.00000000d-02 ,  1.21000000d-01 /
      data f( 9, 7),f( 9, 8) /  3.11000000d-01 ,  8.27000000d-01 /
      data f( 9, 9),f( 9,10) /  1.98400000d+00 ,  4.09200000d+00 /
      data f( 9,11),f( 9,12) /  7.11300000d+00 ,  1.02010000d+01 /
      data f( 9,13),f( 9,14) /  1.02010000d+01 ,  7.11300000d+00 /
      data f( 9,15) /           4.09200000d+00 /
      data f(10, 1),f(10, 2) / -1.69000000d-01 , -5.20000000d-02 /
      data f(10, 3),f(10, 4) / -2.20000000d-02 ,  0.00000000d+00 /
      data f(10, 5),f(10, 6) /  2.20000000d-02 ,  5.20000000d-02 /
      data f(10, 7),f(10, 8) /  1.69000000d-01 ,  5.31000000d-01 /
      data f(10, 9),f(10,10) /  1.37000000d+00 ,  2.90200000d+00 /
      data f(10,11),f(10,12) /  5.09100000d+00 ,  7.35600000d+00 /
      data f(10,13),f(10,14) /  7.35600000d+00 ,  5.09100000d+00 /
      data f(10,15) /           2.90200000d+00 /
      data f(11, 1),f(11, 2) / -7.60000000d-02 , -1.20000000d-02 /
      data f(11, 3),f(11, 4) / -5.00000000d-03 ,  0.00000000d+00 /
      data f(11, 5),f(11, 6) /  5.00000000d-03 ,  1.20000000d-02 /
      data f(11, 7),f(11, 8) /  7.60000000d-02 ,  3.21000000d-01 /
      data f(11, 9),f(11,10) /  9.16000000d-01 ,  2.00500000d+00 /
      data f(11,11),f(11,12) /  3.56200000d+00 ,  5.19000000d+00 /
      data f(11,13),f(11,14) /  5.19000000d+00 ,  3.56200000d+00 /
      data f(11,15) /           2.00500000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  1.60000000d-02 /
      data f(12, 3),f(12, 4) /  7.00000000d-03 ,  0.00000000d+00 /
      data f(12, 5),f(12, 6) / -7.00000000d-03 , -1.60000000d-02 /
      data f(12, 7),f(12, 8) /  0.00000000d+00 ,  1.22000000d-01 /
      data f(12, 9),f(12,10) /  4.53000000d-01 ,  1.08200000d+00 /
      data f(12,11),f(12,12) /  1.98000000d+00 ,  2.94400000d+00 /
      data f(12,13),f(12,14) /  2.94400000d+00 ,  1.98000000d+00 /
      data f(12,15) /           1.08200000d+00 /
      data f(13, 1),f(13, 2) /  4.00000000d-02 ,  2.90000000d-02 /
      data f(13, 3),f(13, 4) /  1.30000000d-02 ,  0.00000000d+00 /
      data f(13, 5),f(13, 6) / -1.30000000d-02 , -2.90000000d-02 /
      data f(13, 7),f(13, 8) / -4.00000000d-02 , -2.40000000d-02 /
      data f(13, 9),f(13,10) /  6.60000000d-02 ,  2.69000000d-01 /
      data f(13,11),f(13,12) /  5.74000000d-01 ,  9.14000000d-01 /
      data f(13,13),f(13,14) /  9.14000000d-01 ,  5.74000000d-01 /
      data f(13,15) /           2.69000000d-01 /
      data f(14, 1),f(14, 2) /  4.60000000d-02 ,  2.90000000d-02 /
      data f(14, 3),f(14, 4) /  1.30000000d-02 ,  0.00000000d+00 /
      data f(14, 5),f(14, 6) / -1.30000000d-02 , -2.90000000d-02 /
      data f(14, 7),f(14, 8) / -4.60000000d-02 , -6.20000000d-02 /
      data f(14, 9),f(14,10) / -6.40000000d-02 , -4.10000000d-02 /
      data f(14,11),f(14,12) /  1.00000000d-03 ,  2.50000000d-02 /
      data f(14,13),f(14,14) /  2.50000000d-02 ,  1.00000000d-03 /
      data f(14,15) /          -4.10000000d-02 /
      data f(15, 1),f(15, 2) /  4.30000000d-02 ,  2.40000000d-02 /
      data f(15, 3),f(15, 4) /  1.00000000d-02 ,  0.00000000d+00 /
      data f(15, 5),f(15, 6) / -1.00000000d-02 , -2.40000000d-02 /
      data f(15, 7),f(15, 8) / -4.30000000d-02 , -6.50000000d-02 /
      data f(15, 9),f(15,10) / -9.00000000d-02 , -1.24000000d-01 /
      data f(15,11),f(15,12) / -1.75000000d-01 , -2.42000000d-01 /
      data f(15,13),f(15,14) / -2.42000000d-01 , -1.75000000d-01 /
      data f(15,15) /          -1.24000000d-01 /
      data f(16, 1),f(16, 2) /  3.50000000d-02 ,  2.00000000d-02 /
      data f(16, 3),f(16, 4) /  8.00000000d-03 ,  0.00000000d+00 /
      data f(16, 5),f(16, 6) / -8.00000000d-03 , -2.00000000d-02 /
      data f(16, 7),f(16, 8) / -3.50000000d-02 , -5.80000000d-02 /
      data f(16, 9),f(16,10) / -8.80000000d-02 , -1.33000000d-01 /
      data f(16,11),f(16,12) / -1.93000000d-01 , -2.53000000d-01 /
      data f(16,13),f(16,14) / -2.53000000d-01 , -1.93000000d-01 /
      data f(16,15) /          -1.33000000d-01 /
      data f(17, 1),f(17, 2) /  2.80000000d-02 ,  1.50000000d-02 /
      data f(17, 3),f(17, 4) /  6.00000000d-03 ,  0.00000000d+00 /
      data f(17, 5),f(17, 6) / -6.00000000d-03 , -1.50000000d-02 /
      data f(17, 7),f(17, 8) / -2.80000000d-02 , -4.80000000d-02 /
      data f(17, 9),f(17,10) / -7.80000000d-02 , -1.19000000d-01 /
      data f(17,11),f(17,12) / -1.69000000d-01 , -2.10000000d-01 /
      data f(17,13),f(17,14) / -2.10000000d-01 , -1.69000000d-01 /
      data f(17,15) /          -1.19000000d-01 /
      data f(18, 1),f(18, 2) /  1.70000000d-02 ,  9.00000000d-03 /
      data f(18, 3),f(18, 4) /  4.00000000d-03 ,  0.00000000d+00 /
      data f(18, 5),f(18, 6) / -4.00000000d-03 , -9.00000000d-03 /
      data f(18, 7),f(18, 8) / -1.70000000d-02 , -2.90000000d-02 /
      data f(18, 9),f(18,10) / -4.60000000d-02 , -7.00000000d-02 /
      data f(18,11),f(18,12) / -9.50000000d-02 , -1.15000000d-01 /
      data f(18,13),f(18,14) / -1.15000000d-01 , -9.50000000d-02 /
      data f(18,15) /          -7.00000000d-02 /
      data f(19, 1),f(19, 2) /  9.00000000d-03 ,  5.00000000d-03 /
      data f(19, 3),f(19, 4) /  2.00000000d-03 ,  0.00000000d+00 /
      data f(19, 5),f(19, 6) / -2.00000000d-03 , -5.00000000d-03 /
      data f(19, 7),f(19, 8) / -9.00000000d-03 , -1.60000000d-02 /
      data f(19, 9),f(19,10) / -2.50000000d-02 , -3.60000000d-02 /
      data f(19,11),f(19,12) / -4.70000000d-02 , -5.60000000d-02 /
      data f(19,13),f(19,14) / -5.60000000d-02 , -4.70000000d-02 /
      data f(19,15) /          -3.60000000d-02 /
      data f(20, 1),f(20, 2) /  5.00000000d-03 ,  2.00000000d-03 /
      data f(20, 3),f(20, 4) /  1.00000000d-03 ,  0.00000000d+00 /
      data f(20, 5),f(20, 6) / -1.00000000d-03 , -2.00000000d-03 /
      data f(20, 7),f(20, 8) / -5.00000000d-03 , -9.00000000d-03 /
      data f(20, 9),f(20,10) / -1.50000000d-02 , -2.10000000d-02 /
      data f(20,11),f(20,12) / -2.80000000d-02 , -3.40000000d-02 /
      data f(20,13),f(20,14) / -3.40000000d-02 , -2.80000000d-02 /
      data f(20,15) /          -2.10000000d-02 /
      data f(21, 1),f(21, 2) /  3.00000000d-03 ,  2.00000000d-03 /
      data f(21, 3),f(21, 4) /  1.00000000d-03 ,  0.00000000d+00 /
      data f(21, 5),f(21, 6) / -1.00000000d-03 , -2.00000000d-03 /
      data f(21, 7),f(21, 8) / -3.00000000d-03 , -6.00000000d-03 /
      data f(21, 9),f(21,10) / -1.00000000d-02 , -1.50000000d-02 /
      data f(21,11),f(21,12) / -2.00000000d-02 , -2.50000000d-02 /
      data f(21,13),f(21,14) / -2.50000000d-02 , -2.00000000d-02 /
      data f(21,15) /          -1.50000000d-02 /
      data f(22, 1),f(22, 2) /  3.00000000d-03 ,  1.00000000d-03 /
      data f(22, 3),f(22, 4) /  1.00000000d-03 ,  0.00000000d+00 /
      data f(22, 5),f(22, 6) / -1.00000000d-03 , -1.00000000d-03 /
      data f(22, 7),f(22, 8) / -3.00000000d-03 , -5.00000000d-03 /
      data f(22, 9),f(22,10) / -8.00000000d-03 , -1.30000000d-02 /
      data f(22,11),f(22,12) / -1.70000000d-02 , -2.10000000d-02 /
      data f(22,13),f(22,14) / -2.10000000d-02 , -1.70000000d-02 /
      data f(22,15) /          -1.30000000d-02 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/-1.57047320d+01,-1.09250170d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-1.00937062d+01,-7.01996594d-03/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-4.68083045d+00,-2.57511920d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 0.00000000d+00,-1.95572449d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 4.68083045d+00, 2.65334818d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.00937062d+01, 6.74616451d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.57047320d+01, 1.19419938d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.37778188d+01, 2.91258604d-02/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 3.80638421d+01, 3.38545645d-02/
      data fpp( 1,10,1),fpp( 1,10,2)/ 7.84587974d+01, 9.99358814d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.74947946d+02, 8.16219098d-02/
      data fpp( 1,12,1),fpp( 1,12,2)/ 5.24039606d+02,-1.85283521d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 5.24039606d+02,-1.80240393d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.74947946d+02, 6.14493979d-02/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.84587974d+01, 1.75582801d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/-1.37905361d+01,-9.54465502d-03/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-8.51258764d+00,-6.10068997d-03/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-3.93833910d+00,-2.11258511d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 0.00000000d+00,-2.89695784d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.93833910d+00, 2.22846343d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 8.51258764d+00, 5.69511587d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.37905361d+01, 1.10510731d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 2.10443624d+01, 2.42005918d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 3.37723158d+01, 3.46465599d-02/
      data fpp( 2,10,1),fpp( 2,10,2)/ 6.06824052d+01, 7.04931687d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.57104108d+02, 4.24807651d-02/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.96120787d+02,-1.33316229d-01/
      data fpp( 2,13,1),fpp( 2,13,2)/ 3.96120787d+02,-1.30751695d-01/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.57104108d+02, 3.22226271d-02/
      data fpp( 2,15,1),fpp( 2,15,2)/ 6.06824052d+01, 1.08961186d-01/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-1.19331238d+01,-8.23698505d-03/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-6.85594327d+00,-5.24602991d-03/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-2.96581315d+00,-1.73889532d-03/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 0.00000000d+00,-3.83888276d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 2.96581315d+00, 1.89245063d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 6.85594327d+00, 4.70858632d-03/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.19331238d+01, 1.02332041d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 1.80447317d+01, 2.00585974d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 3.02468947d+01, 3.56524065d-02/
      data fpp( 3,10,1),fpp( 3,10,2)/ 3.16115817d+01, 4.66117767d-02/
      data fpp( 3,11,1),fpp( 3,11,2)/ 1.56635621d+02, 2.66604868d-02/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.39877244d+02,-1.02073724d-01/
      data fpp( 3,13,1),fpp( 3,13,2)/ 2.39877244d+02,-1.00169072d-01/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.56635621d+02, 1.90418776d-02/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.16115817d+01, 7.51815612d-02/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-9.87696866d+00,-7.10052680d-03/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-5.46363927d+00,-4.50894639d-03/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-2.79840828d+00,-1.44368763d-03/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 0.00000000d+00,-3.63031000d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 2.79840828d+00, 1.58890003d-03/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 5.46363927d+00, 4.00070299d-03/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 9.87696866d+00, 8.98828800d-03/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 1.59767107d+01, 1.81261450d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.40401052d+01, 3.14271320d-02/
      data fpp( 4,10,1),fpp( 4,10,2)/ 4.20712680d+01, 3.95253269d-02/
      data fpp( 4,11,1),fpp( 4,11,2)/ 9.53534085d+01, 8.83156026d-03/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.82770235d+02,-7.86315680d-02/
      data fpp( 4,13,1),fpp( 4,13,2)/ 1.82770235d+02,-7.81510762d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/ 9.53534085d+01, 6.90959320d-03/
      data fpp( 4,15,1),fpp( 4,15,2)/ 4.20712680d+01, 4.67327034d-02/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-7.95900155d+00,-6.22036795d-03/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-4.88949964d+00,-3.91926410d-03/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-2.04055371d+00,-1.14257564d-03/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 0.00000000d+00,-3.04333266d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 2.04055371d+00, 1.26430895d-03/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 4.88949964d+00, 3.49319753d-03/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 7.95900155d+00, 7.80290093d-03/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.28484254d+01, 1.67151987d-02/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.99926843d+01, 2.66763041d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/ 3.47033462d+01, 3.76795849d-02/
      data fpp( 5,11,1),fpp( 5,11,2)/ 4.21507450d+01,-3.69464355d-03/
      data fpp( 5,12,1),fpp( 5,12,2)/ 1.45041814d+02,-6.26210107d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 1.45041814d+02,-6.32796463d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/ 4.21507450d+01,-1.06010107d-03/
      data fpp( 5,15,1),fpp( 5,15,2)/ 3.47033462d+01, 2.78000505d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-5.68451101d+00,-4.66496755d-03/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-2.99968143d+00,-2.91006491d-03/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.47913473d+00,-7.34772827d-04/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 0.00000000d+00,-3.08437840d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.47913473d+00, 8.58147964d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 2.99968143d+00, 2.47825193d-03/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 5.68451101d+00, 6.26884432d-03/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 8.81636833d+00, 1.32463708d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.30518944d+01, 2.24656725d-02/
      data fpp( 6,10,1),fpp( 6,10,2)/ 1.85543273d+01, 2.50509393d-02/
      data fpp( 6,11,1),fpp( 6,11,2)/ 3.60210606d+01, 4.83057039d-03/
      data fpp( 6,12,1),fpp( 6,12,2)/ 1.01894400d+01,-5.19332208d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/ 1.01894400d+01,-5.16356227d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/ 3.60210606d+01, 3.64017792d-03/
      data fpp( 6,15,1),fpp( 6,15,2)/ 1.85543273d+01, 2.95149110d-02/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-3.50295441d+00,-3.67410430d-03/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-2.16177462d+00,-2.24179139d-03/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-8.92907373d-01,-3.78730134d-04/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 0.00000000d+00,-2.32880721d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 8.92907373d-01, 4.71882423d-04/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 2.16177462d+00, 1.91575838d-03/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 3.50295441d+00, 4.88508405d-03/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 5.88610125d+00, 1.05839054d-02/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 9.09973807d+00, 1.78792943d-02/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.36293447d+01, 1.86989173d-02/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.85150126d+01, 5.00503656d-03/
      data fpp( 7,12,1),fpp( 7,12,2)/ 3.45504259d+01,-4.07590635d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 3.45504259d+01,-4.03953278d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.85150126d+01, 3.55009365d-03/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.36293447d+01, 2.41549532d-02/
      data fpp( 8, 1,1),fpp( 8, 1,2)/-2.50367133d+00,-2.81689179d-03/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.40322009d+00,-1.68621642d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-6.49235780d-01,-1.58242540d-04/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 0.00000000d+00,-2.08134214d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 6.49235780d-01, 2.41496226d-04/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 1.40322009d+00, 1.39482852d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 2.50367133d+00, 3.89918971d-03/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 3.93922668d+00, 8.26841266d-03/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 6.34915331d+00, 1.36671597d-02/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.06282939d+01, 1.38629487d-02/
      data fpp( 8,11,1),fpp( 8,11,2)/ 1.71188891d+01, 4.98104565d-03/
      data fpp( 8,12,1),fpp( 8,12,2)/ 1.91588565d+01,-3.15671313d-02/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.91588565d+01,-3.11791291d-02/
      data fpp( 8,14,1),fpp( 8,14,2)/ 1.71188891d+01, 3.42903687d-03/
      data fpp( 8,15,1),fpp( 8,15,2)/ 1.06282939d+01, 1.96829816d-02/
      data fpp( 9, 1,1),fpp( 9, 1,2)/-1.78236025d+00,-2.13970137d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-9.25345029d-01,-1.25059726d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-4.10149508d-01, 2.09040567d-06/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 0.00000000d+00,-1.77643639d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 4.10149508d-01, 6.89670500d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 9.25345029d-01, 1.00189616d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.78236025d+00, 3.06344829d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 2.95699203d+00, 6.30431066d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 5.10364869d+00, 1.01793091d-02/
      data fpp( 9,10,1),fpp( 9,10,2)/ 9.10747967d+00, 1.00384531d-02/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.45594312d+01, 4.44687865d-03/
      data fpp( 9,12,1),fpp( 9,12,2)/ 2.05141480d+01,-2.38059677d-02/
      data fpp( 9,13,1),fpp( 9,13,2)/ 2.05141480d+01,-2.34455363d-02/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.45594312d+01, 3.00515323d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 9.10747967d+00, 1.54449234d-02/
      data fpp(10, 1,1),fpp(10, 1,2)/-1.16688765d+00,-1.63219037d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/-7.45399795d-01,-9.25619251d-04/
      data fpp(10, 3,1),fpp(10, 3,2)/-2.60166189d-01, 1.14667379d-04/
      data fpp(10, 4,1),fpp(10, 4,2)/ 0.00000000d+00,-1.30502646d-05/
      data fpp(10, 5,1),fpp(10, 5,2)/ 2.60166189d-01,-6.24663206d-05/
      data fpp(10, 6,1),fpp(10, 6,2)/ 7.45399795d-01, 7.42915547d-04/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.16688765d+00, 2.31080413d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/ 2.08280519d+00, 4.71386792d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/ 3.98625192d+00, 7.45372417d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 7.24178740d+00, 7.05123539d-03/
      data fpp(10,11,1),fpp(10,11,2)/ 1.23933862d+01, 3.76133428d-03/
      data fpp(10,12,1),fpp(10,12,2)/ 1.68345516d+01,-1.75365725d-02/
      data fpp(10,13,1),fpp(10,13,2)/ 1.68345516d+01,-1.72209496d-02/
      data fpp(10,14,1),fpp(10,14,2)/ 1.23933862d+01, 2.49884275d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 7.24178740d+00, 1.17855786d-02/
      data fpp(11, 1,1),fpp(11, 1,2)/-9.00089155d-01,-1.09914757d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/-4.43055793d-01,-6.11704853d-04/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.99185738d-01, 1.25966985d-04/
      data fpp(11, 4,1),fpp(11, 4,2)/ 0.00000000d+00,-1.21630888d-05/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.99185738d-01,-7.73146304d-05/
      data fpp(11, 6,1),fpp(11, 6,2)/ 4.43055793d-01, 4.41421610d-04/
      data fpp(11, 7,1),fpp(11, 7,2)/ 9.00089155d-01, 1.73162819d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.61178721d+00, 3.49206563d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/ 2.95134363d+00, 5.30010928d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 5.87537072d+00, 4.94749725d-03/
      data fpp(11,11,1),fpp(11,11,2)/ 9.81702397d+00, 2.98990171d-03/
      data fpp(11,12,1),fpp(11,12,2)/ 1.39976456d+01,-1.26471041d-02/
      data fpp(11,13,1),fpp(11,13,2)/ 1.39976456d+01,-1.23936386d-02/
      data fpp(11,14,1),fpp(11,14,2)/ 9.81702397d+00, 1.97603959d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 5.87537072d+00, 8.74948020d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/-4.55111053d-01,-5.08449124d-04/
      data fpp(12, 2,1),fpp(12, 2,2)/-1.59547495d-01,-2.73101751d-04/
      data fpp(12, 3,1),fpp(12, 3,2)/-6.26034144d-02, 1.00856130d-04/
      data fpp(12, 4,1),fpp(12, 4,2)/ 0.00000000d+00,-1.03227680d-05/
      data fpp(12, 5,1),fpp(12, 5,2)/ 6.26034144d-02,-5.95650579d-05/
      data fpp(12, 6,1),fpp(12, 6,2)/ 1.59547495d-01, 1.28583000d-04/
      data fpp(12, 7,1),fpp(12, 7,2)/ 4.55111053d-01, 1.04523306d-03/
      data fpp(12, 8,1),fpp(12, 8,2)/ 9.72172521d-01, 2.05048476d-03/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.03801996d+00, 3.29282789d-03/
      data fpp(12,10,1),fpp(12,10,2)/ 3.75423935d+00, 2.65820366d-03/
      data fpp(12,11,1),fpp(12,11,2)/ 6.44766264d+00, 2.21435745d-03/
      data fpp(12,12,1),fpp(12,12,2)/ 8.98481366d+00,-7.55563348d-03/
      data fpp(12,13,1),fpp(12,13,2)/ 8.98481366d+00,-7.36027828d-03/
      data fpp(12,14,1),fpp(12,14,2)/ 6.44766264d+00, 1.43293665d-03/
      data fpp(12,15,1),fpp(12,15,2)/ 3.75423935d+00, 5.58853167d-03/
      data fpp(13, 1,1),fpp(13, 1,2)/-8.35911375d-02,-1.25706222d-04/
      data fpp(13, 2,1),fpp(13, 2,2)/-3.16145396d-02,-5.85875568d-05/
      data fpp(13, 3,1),fpp(13, 3,2)/-1.61576312d-02, 6.00564487d-05/
      data fpp(13, 4,1),fpp(13, 4,2)/ 0.00000000d+00,-1.63823805d-06/
      data fpp(13, 5,1),fpp(13, 5,2)/ 1.61576312d-02,-5.35034965d-05/
      data fpp(13, 6,1),fpp(13, 6,2)/ 3.16145396d-02, 3.56522241d-05/
      data fpp(13, 7,1),fpp(13, 7,2)/ 8.35911375d-02, 2.10894600d-04/
      data fpp(13, 8,1),fpp(13, 8,2)/ 3.77975608d-01, 7.40769375d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/ 9.39529945d-01, 1.26602790d-03/
      data fpp(13,10,1),fpp(13,10,2)/ 1.86921166d+00, 9.75119025d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 3.01326518d+00, 9.53495998d-04/
      data fpp(13,12,1),fpp(13,12,2)/ 3.97000893d+00,-2.68910302d-03/
      data fpp(13,13,1),fpp(13,13,2)/ 3.97000893d+00,-2.60943894d-03/
      data fpp(13,14,1),fpp(13,14,2)/ 3.01326518d+00, 6.34839698d-04/
      data fpp(13,15,1),fpp(13,15,2)/ 1.86921166d+00, 2.17008015d-03/
      data fpp(14, 1,1),fpp(14, 1,2)/-2.65243969d-02,-8.09994773d-06/
      data fpp(14, 2,1),fpp(14, 2,2)/-2.59943464d-02, 6.19989546d-06/
      data fpp(14, 3,1),fpp(14, 3,2)/-1.67660610d-02, 4.33003659d-05/
      data fpp(14, 4,1),fpp(14, 4,2)/ 0.00000000d+00, 5.98641037d-07/
      data fpp(14, 5,1),fpp(14, 5,2)/ 1.67660610d-02,-4.56949300d-05/
      data fpp(14, 6,1),fpp(14, 6,2)/ 2.59943464d-02, 2.18107906d-06/
      data fpp(14, 7,1),fpp(14, 7,2)/ 2.65243969d-02,-2.30293862d-05/
      data fpp(14, 8,1),fpp(14, 8,2)/ 1.07925045d-01, 1.49936466d-04/
      data fpp(14, 9,1),fpp(14, 9,2)/ 3.71860258d-01, 2.63283523d-04/
      data fpp(14,10,1),fpp(14,10,2)/ 8.40914015d-01, 2.96929441d-04/
      data fpp(14,11,1),fpp(14,11,2)/ 1.49127665d+00,-3.11001289d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 2.51915060d+00,-1.32924285d-04/
      data fpp(14,13,1),fpp(14,13,2)/ 2.51915060d+00,-1.65726500d-04/
      data fpp(14,14,1),fpp(14,14,2)/ 1.49127665d+00,-1.79792429d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 8.40914015d-01,-1.95103786d-04/
      data fpp(15, 1,1),fpp(15, 1,2)/-2.63112750d-02, 6.11789915d-05/
      data fpp(15, 2,1),fpp(15, 2,2)/ 1.55919253d-02, 4.76420170d-05/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.12218752d-02, 4.82529405d-05/
      data fpp(15, 4,1),fpp(15, 4,2)/ 0.00000000d+00,-6.53778832d-07/
      data fpp(15, 5,1),fpp(15, 5,2)/-1.12218752d-02,-4.56378251d-05/
      data fpp(15, 6,1),fpp(15, 6,2)/-1.55919253d-02,-5.67949207d-05/
      data fpp(15, 7,1),fpp(15, 7,2)/ 2.63112750d-02,-2.71824922d-05/
      data fpp(15, 8,1),fpp(15, 8,2)/ 3.03242120d-02,-1.44751105d-05/
      data fpp(15, 9,1),fpp(15, 9,2)/ 6.90290217d-02,-9.49170659d-05/
      data fpp(15,10,1),fpp(15,10,2)/ 2.15132279d-01,-1.45856626d-04/
      data fpp(15,11,1),fpp(15,11,2)/ 5.49628207d-01,-3.41656430d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 8.81388667d-01, 5.52482345d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 8.81388667d-01, 5.23381179d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 5.49628207d-01,-2.25251765d-04/
      data fpp(15,15,1),fpp(15,15,2)/ 2.15132279d-01,-5.82374117d-04/
      data fpp(16, 1,1),fpp(16, 1,2)/ 1.17694970d-02, 2.18570370d-05/
      data fpp(16, 2,1),fpp(16, 2,2)/-1.23733547d-02, 2.62859260d-05/
      data fpp(16, 3,1),fpp(16, 3,2)/-4.12143961d-03, 5.29992589d-05/
      data fpp(16, 4,1),fpp(16, 4,2)/ 0.00000000d+00, 1.71703847d-06/
      data fpp(16, 5,1),fpp(16, 5,2)/ 4.12143961d-03,-5.98674128d-05/
      data fpp(16, 6,1),fpp(16, 6,2)/ 1.23733547d-02,-2.24738740d-06/
      data fpp(16, 7,1),fpp(16, 7,2)/-1.17694970d-02,-1.11143038d-04/
      data fpp(16, 8,1),fpp(16, 8,2)/ 1.07781071d-02,-3.31804621d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.40236550d-02,-1.76135114d-04/
      data fpp(16,10,1),fpp(16,10,2)/ 7.45568678d-02,-1.62279082d-04/
      data fpp(16,11,1),fpp(16,11,2)/ 1.02210517d-01,-7.47485577d-05/
      data fpp(16,12,1),fpp(16,12,2)/ 9.92947306d-02, 4.61273313d-04/
      data fpp(16,13,1),fpp(16,13,2)/ 9.92947306d-02, 4.53554341d-04/
      data fpp(16,14,1),fpp(16,14,2)/ 1.02210517d-01,-4.38726687d-05/
      data fpp(16,15,1),fpp(16,15,2)/ 7.45568678d-02,-2.78063666d-04/
      data fpp(17, 1,1),fpp(17, 1,2)/ 3.23328721d-03, 5.07548812d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/ 9.90149364d-03, 3.84902376d-05/
      data fpp(17, 3,1),fpp(17, 3,2)/ 5.26388329d-03, 3.52841684d-05/
      data fpp(17, 4,1),fpp(17, 4,2)/ 0.00000000d+00, 3.73088901d-07/
      data fpp(17, 5,1),fpp(17, 5,2)/-5.26388329d-03,-3.67765240d-05/
      data fpp(17, 6,1),fpp(17, 6,2)/-9.90149364d-03,-3.32669930d-05/
      data fpp(17, 7,1),fpp(17, 7,2)/-3.23328721d-03,-7.01555041d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/-1.43664043d-03,-1.06110991d-04/
      data fpp(17, 9,1),fpp(17, 9,2)/ 2.68763584d-02,-1.05400533d-04/
      data fpp(17,10,1),fpp(17,10,2)/ 3.86402497d-02,-1.32286877d-04/
      data fpp(17,11,1),fpp(17,11,2)/ 4.95297233d-02, 9.45480421d-05/
      data fpp(17,12,1),fpp(17,12,2)/ 1.74324103d-02, 2.94094709d-04/
      data fpp(17,13,1),fpp(17,13,2)/ 1.74324103d-02, 3.00441852d-04/
      data fpp(17,14,1),fpp(17,14,2)/ 4.95297233d-02, 6.91594709d-05/
      data fpp(17,15,1),fpp(17,15,2)/ 3.86402497d-02,-3.70797355d-05/
      data fpp(18, 1,1),fpp(18, 1,2)/ 2.41538989d-03, 4.96213813d-05/
      data fpp(18, 2,1),fpp(18, 2,2)/ 4.82196457d-04, 3.07572373d-05/
      data fpp(18, 3,1),fpp(18, 3,2)/-1.73093005d-03, 7.34966938d-06/
      data fpp(18, 4,1),fpp(18, 4,2)/ 0.00000000d+00,-1.55914831d-07/
      data fpp(18, 5,1),fpp(18, 5,2)/ 1.73093005d-03,-6.72601005d-06/
      data fpp(18, 6,1),fpp(18, 6,2)/-4.82196457d-04,-3.29400450d-05/
      data fpp(18, 7,1),fpp(18, 7,2)/-2.41538989d-03,-4.15138101d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/-7.07913226d-03,-4.10047146d-05/
      data fpp(18, 9,1),fpp(18, 9,2)/-2.06409027d-02,-9.44673316d-05/
      data fpp(18,10,1),fpp(18,10,2)/-2.71991829d-02,-1.12595896d-06/
      data fpp(18,11,1),fpp(18,11,2)/-4.36944286d-02, 3.89711674d-05/
      data fpp(18,12,1),fpp(18,12,2)/-4.79445961d-02, 1.45241289d-04/
      data fpp(18,13,1),fpp(18,13,2)/-4.79445961d-02, 1.44790549d-04/
      data fpp(18,14,1),fpp(18,14,2)/-4.36944286d-02, 4.07741289d-05/
      data fpp(18,15,1),fpp(18,15,2)/-2.71991829d-02,-7.88706446d-06/
      data fpp(19, 1,1),fpp(19, 1,2)/ 5.10515323d-03, 1.03614926d-05/
      data fpp(19, 2,1),fpp(19, 2,2)/ 1.69720537d-04, 9.27701479d-06/
      data fpp(19, 3,1),fpp(19, 3,2)/ 1.65983692d-03, 1.25304482d-05/
      data fpp(19, 4,1),fpp(19, 4,2)/ 0.00000000d+00, 6.01192257d-07/
      data fpp(19, 5,1),fpp(19, 5,2)/-1.65983692d-03,-1.49352173d-05/
      data fpp(19, 6,1),fpp(19, 6,2)/-1.69720537d-04,-8.60323193d-07/
      data fpp(19, 7,1),fpp(19, 7,2)/-5.10515323d-03,-4.16234900d-05/
      data fpp(19, 8,1),fpp(19, 8,2)/-6.24683053d-03,-1.26457170d-05/
      data fpp(19, 9,1),fpp(19, 9,2)/-1.03127477d-02,-2.77936422d-05/
      data fpp(19,10,1),fpp(19,10,2)/-1.98435180d-02, 3.82028573d-06/
      data fpp(19,11,1),fpp(19,11,2)/-3.07520088d-02, 1.25124993d-05/
      data fpp(19,12,1),fpp(19,12,2)/-4.16540258d-02, 6.61297171d-05/
      data fpp(19,13,1),fpp(19,13,2)/-4.16540258d-02, 6.53545990d-05/
      data fpp(19,14,1),fpp(19,14,2)/-3.07520088d-02, 1.56129717d-05/
      data fpp(19,15,1),fpp(19,15,2)/-1.98435180d-02,-7.80648586d-06/
      data fpp(20, 1,1),fpp(20, 1,2)/ 1.16399721d-03, 3.92609678d-05/
      data fpp(20, 2,1),fpp(20, 2,2)/ 4.83892139d-03, 2.14780644d-05/
      data fpp(20, 3,1),fpp(20, 3,2)/ 1.09158237d-03,-5.17322557d-06/
      data fpp(20, 4,1),fpp(20, 4,2)/ 0.00000000d+00,-7.85162174d-07/
      data fpp(20, 5,1),fpp(20, 5,2)/-1.09158237d-03, 8.31387426d-06/
      data fpp(20, 6,1),fpp(20, 6,2)/-4.83892139d-03,-3.24703349d-05/
      data fpp(20, 7,1),fpp(20, 7,2)/-1.16399721d-03, 1.56746527d-06/
      data fpp(20, 8,1),fpp(20, 8,2)/-3.93354563d-03,-3.37995262d-05/
      data fpp(20, 9,1),fpp(20, 9,2)/-4.10810638d-03, 1.36306396d-05/
      data fpp(20,10,1),fpp(20,10,2)/-7.42674516d-03,-2.07230320d-05/
      data fpp(20,11,1),fpp(20,11,2)/-7.29753609d-03, 9.26148861d-06/
      data fpp(20,12,1),fpp(20,12,2)/-7.43930075d-03, 4.36770776d-05/
      data fpp(20,13,1),fpp(20,13,2)/-7.43930075d-03, 4.43380228d-05/
      data fpp(20,14,1),fpp(20,14,2)/-7.29753609d-03, 6.61770776d-06/
      data fpp(20,15,1),fpp(20,15,2)/-7.42674516d-03,-1.08088539d-05/
      data fpp(21, 1,1),fpp(21, 1,2)/ 2.23885794d-03,-2.31662812d-08/
      data fpp(21, 2,1),fpp(21, 2,2)/-1.52540611d-03, 4.63325623d-08/
      data fpp(21, 3,1),fpp(21, 3,2)/-2.61663907d-05,-1.62163968d-07/
      data fpp(21, 4,1),fpp(21, 4,2)/ 0.00000000d+00, 6.02323310d-07/
      data fpp(21, 5,1),fpp(21, 5,2)/ 2.61663907d-05,-2.24712927d-06/
      data fpp(21, 6,1),fpp(21, 6,2)/ 1.52540611d-03, 8.38619378d-06/
      data fpp(21, 7,1),fpp(21, 7,2)/-2.23885794d-03,-3.12976459d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/-2.01898696d-03,-3.19561034d-06/
      data fpp(21, 9,1),fpp(21, 9,2)/-3.25482675d-03,-1.59199128d-05/
      data fpp(21,10,1),fpp(21,10,2)/-4.44950138d-03, 6.87526143d-06/
      data fpp(21,11,1),fpp(21,11,2)/-6.05784683d-03,-1.15811329d-05/
      data fpp(21,12,1),fpp(21,12,2)/-6.58877121d-03, 3.94492704d-05/
      data fpp(21,13,1),fpp(21,13,2)/-6.58877121d-03, 3.74427554d-05/
      data fpp(21,14,1),fpp(21,14,2)/-6.05784683d-03,-3.55507296d-06/
      data fpp(21,15,1),fpp(21,15,2)/-4.44950138d-03,-2.32224635d-05/
      data fpp(22, 1,1),fpp(22, 1,2)/ 1.88057103d-03, 4.85064311d-05/
      data fpp(22, 2,1),fpp(22, 2,2)/-4.73729694d-03, 2.29871379d-05/
      data fpp(22, 3,1),fpp(22, 3,2)/-9.86916805d-04,-2.04549826d-05/
      data fpp(22, 4,1),fpp(22, 4,2)/ 0.00000000d+00,-1.16720747d-06/
      data fpp(22, 5,1),fpp(22, 5,2)/ 9.86916805d-04, 2.51238125d-05/
      data fpp(22, 6,1),fpp(22, 6,2)/ 4.73729694d-03,-3.93280424d-05/
      data fpp(22, 7,1),fpp(22, 7,2)/-1.88057103d-03, 1.21883573d-05/
      data fpp(22, 8,1),fpp(22, 8,2)/ 9.49348185d-06,-9.42538659d-06/
      data fpp(22, 9,1),fpp(22, 9,2)/-8.72586625d-04,-3.44868109d-05/
      data fpp(22,10,1),fpp(22,10,2)/ 1.22475069d-03, 2.73726302d-05/
      data fpp(22,11,1),fpp(22,11,2)/ 1.52892342d-03,-1.50037099d-05/
      data fpp(22,12,1),fpp(22,12,2)/ 3.79438561d-03, 3.26422094d-05/
      data fpp(22,13,1),fpp(22,13,2)/ 3.79438561d-03, 2.95752267d-05/
      data fpp(22,14,1),fpp(22,14,2)/ 1.52892342d-03,-2.73577906d-06/
      data fpp(22,15,1),fpp(22,15,2)/ 1.22475069d-03,-1.86321105d-05/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 2.91199352d-03,-1.09058451d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/-1.04386592d-02,-1.07749555d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.47486414d-02,-1.39943533d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 5.31177753d-02,-5.07530850d-02/
      data fpppp( 1, 9),fpppp( 1,10)/ 5.22670751d-01,-4.73394001d-01/
      data fpppp( 1,11),fpppp( 1,12)/ 4.73655684d+00,-3.31668264d+00/
      data fpppp( 1,13),fpppp( 1,14)/-2.89098030d+00, 3.03374749d+00/
      data fpppp( 1,15) /             5.91214108d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/-7.92766191d-03,-6.61148124d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/-7.84840574d-03,-1.49462120d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 8.44625422d-03, 4.51901155d-03/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.56996922d-02, 5.12348936d-02/
      data fpppp( 2, 9),fpppp( 2,10)/ 1.07808360d-01, 3.68459827d-01/
      data fpppp( 2,11),fpppp( 2,12)/ 2.58904915d+00,-2.16895789d+00/
      data fpppp( 2,13),fpppp( 2,14)/-1.95815128d+00, 1.74582271d+00/
      data fpppp( 2,15) /             3.53055900d+00 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.48165709d-02,-1.12330405d-02/
      data fpppp( 3, 3),fpppp( 3, 4)/-1.14742930d-02, 1.67119459d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/ 4.78951463d-03, 3.46297648d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-7.20855479d-02, 3.15778068d-01/
      data fpppp( 3, 9),fpppp( 3,10)/-8.25593419d-01, 2.33634704d+00/
      data fpppp( 3,11),fpppp( 3,12)/-1.10023363d+00,-4.42357464d-01/
      data fpppp( 3,13),fpppp( 3,14)/-6.20059500d-01,-3.89425486d-01/
      data fpppp( 3,15) /            -3.29183489d-01 /
      data fpppp( 4, 1),fpppp( 4, 2)/-3.55330925d-02,-1.90022809d-02/
      data fpppp( 4, 3),fpppp( 4, 4)/ 6.65631223d-03, 3.67669832d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-8.12699156d-03, 2.41496586d-02/
      data fpppp( 4, 7),fpppp( 4, 8)/ 1.64142612d-02, 1.13780566d-02/
      data fpppp( 4, 9),fpppp( 4,10)/ 5.58926606d-02, 3.63117397d-01/
      data fpppp( 4,11),fpppp( 4,12)/ 6.06696418d-01,-7.41821888d-01/
      data fpppp( 4,13),fpppp( 4,14)/-7.00387351d-01, 4.40958273d-01/
      data fpppp( 4,15) /             9.84635440d-01 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 3.16930891d-03,-1.19857222d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.16083787d-02,-8.71446446d-04/
      data fpppp( 5, 5),fpppp( 5, 6)/ 1.50941645d-02,-1.10016780d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/ 4.21459063d-02,-4.83866286d-02/
      data fpppp( 5, 9),fpppp( 5,10)/ 2.86690709d-01,-6.44392024d-01/
      data fpppp( 5,11),fpppp( 5,12)/ 1.85508160d+00,-1.04931416d+00/
      data fpppp( 5,13),fpppp( 5,14)/-8.66330382d-01, 1.12314648d+00/
      data fpppp( 5,15) /             2.10036469d+00 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.24789766d-02,-1.24279506d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/ 2.33380708d-03, 6.08003785d-04/
      data fpppp( 6, 5),fpppp( 6, 6)/-4.76582222d-03, 2.09400036d-02/
      data fpppp( 6, 7),fpppp( 6, 8)/-9.13722020d-03, 4.24305420d-02/
      data fpppp( 6, 9),fpppp( 6,10)/-9.43648228d-02, 4.11043156d-01/
      data fpppp( 6,11),fpppp( 6,12)/-8.31949775d-01, 3.18854710d-01/
      data fpppp( 6,13),fpppp( 6,14)/ 2.34359376d-01,-4.93968441d-01/
      data fpppp( 6,15) /            -8.56386845d-01 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 2.06602359d-03,-2.28477124d-04/
      data fpppp( 7, 3),fpppp( 7, 4)/-5.49086801d-03,-3.65643284d-04/
      data fpppp( 7, 5),fpppp( 7, 6)/ 6.95344115d-03,-4.89052886d-03/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.69474272d-02,-3.81157676d-04/
      data fpppp( 7, 9),fpppp( 7,10)/ 3.44066030d-02,-5.82870658d-02/
      data fpppp( 7,11),fpppp( 7,12)/ 2.20105335d-01,-1.53149550d-01/
      data fpppp( 7,13),fpppp( 7,14)/-1.31666416d-01, 1.34172801d-01/
      data fpppp( 7,15) /             2.63959937d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/-5.83272697d-03,-3.56292234d-03/
      data fpppp( 8, 3),fpppp( 8, 4)/-7.03599874d-04, 9.24101550d-05/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.33959254d-04, 4.85666451d-03/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.02739891d-03, 1.11399860d-02/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.28749339d-02, 4.95131164d-02/
      data fpppp( 8,11),fpppp( 8,12)/-7.82401260d-02,-3.59027236d-03/
      data fpppp( 8,13),fpppp( 8,14)/-1.13081445d-02,-4.73686372d-02/
      data fpppp( 8,15) /            -6.62549664d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/-5.73842960d-03,-3.51320574d-03/
      data fpppp( 9, 3),fpppp( 9, 4)/-7.17929694d-04, 8.21636975d-05/
      data fpppp( 9, 5),fpppp( 9, 6)/ 3.89274904d-04, 4.66349750d-03/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.46591732d-03, 8.52982632d-03/
      data fpppp( 9, 9),fpppp( 9,10)/ 2.27362704d-02, 1.19555513d-02/
      data fpppp( 9,11),fpppp( 9,12)/ 1.63287563d-02,-4.71046602d-02/
      data fpppp( 9,13),fpppp( 9,14)/-4.54919008d-02, 9.87771841d-03/
      data fpppp( 9,15) /             3.61469432d-02 /
      data fpppp(10, 1),fpppp(10, 2)/ 3.33453165d-03, 1.01957270d-03/
      data fpppp(10, 3),fpppp(10, 4)/-3.58807727d-03,-1.71308674d-04/
      data fpppp(10, 5),fpppp(10, 6)/ 4.27331196d-03,-3.41789414d-03/
      data fpppp(10, 7),fpppp(10, 8)/ 5.57351941d-03, 1.07895979d-02/
      data fpppp(10, 9),fpppp(10,10)/ 1.05198403d-02, 2.82563661d-02/
      data fpppp(10,11),fpppp(10,12)/-9.78150522d-03,-3.17563497d-02/
      data fpppp(10,13),fpppp(10,14)/-3.30751603d-02,-4.50626270d-03/
      data fpppp(10,15) /             8.47420669d-03 /
      data fpppp(11, 1),fpppp(11, 2)/-3.76992628d-03,-2.22462648d-03/
      data fpppp(11, 3),fpppp(11, 4)/-1.21366246d-04, 2.90324703d-05/
      data fpppp(11, 5),fpppp(11, 6)/ 5.23636466d-06, 2.63108106d-03/
      data fpppp(11, 7),fpppp(11, 8)/ 2.26023783d-03, 3.60784896d-03/
      data fpppp(11, 9),fpppp(11,10)/ 2.09798686d-02, 7.54091658d-03/
      data fpppp(11,11),fpppp(11,12)/ 9.91403483d-03,-3.28589542d-02/
      data fpppp(11,13),fpppp(11,14)/-3.17988032d-02, 5.67343107d-03/
      data fpppp(11,15) /             2.34431807d-02 /
      data fpppp(12, 1),fpppp(12, 2)/-3.58098364d-03,-2.08219324d-03/
      data fpppp(12, 3),fpppp(12, 4)/-7.41202640d-06, 5.14013524d-05/
      data fpppp(12, 5),fpppp(12, 6)/-1.98193383d-04, 2.80181217d-03/
      data fpppp(12, 7),fpppp(12, 8)/ 9.08113311d-04, 6.85560922d-03/
      data fpppp(12, 9),fpppp(12,10)/ 4.59660818d-03, 1.37802748d-02/
      data fpppp(12,11),fpppp(12,12)/-1.08547282d-03,-1.88147195d-02/
      data fpppp(12,13),fpppp(12,14)/-1.91276356d-02, 1.66191566d-04/
      data fpppp(12,15) /             9.08653331d-03 /
      data fpppp(13, 1),fpppp(13, 2)/-7.21004892d-04,-3.97989147d-04/
      data fpppp(13, 3),fpppp(13, 4)/ 1.21780108d-04,-4.70879224d-05/
      data fpppp(13, 5),fpppp(13, 6)/ 6.65715812d-05,-2.61241767d-04/
      data fpppp(13, 7),fpppp(13, 8)/ 3.16957686d-03, 2.12740672d-03/
      data fpppp(13, 9),fpppp(13,10)/ 4.35098821d-03, 2.55628314d-03/
      data fpppp(13,11),fpppp(13,12)/-1.71381258d-03,-6.93961839d-03/
      data fpppp(13,13),fpppp(13,14)/-7.02655128d-03,-1.36608100d-03/
      data fpppp(13,15) /             1.25228969d-03 /
      data fpppp(14, 1),fpppp(14, 2)/ 1.01724115d-04, 8.07080081d-05/
      data fpppp(14, 3),fpppp(14, 4)/ 9.73379524d-05,-1.77932832d-05/
      data fpppp(14, 5),fpppp(14, 6)/-2.61648194d-05,-3.29813974d-04/
      data fpppp(14, 7),fpppp(14, 8)/ 8.23526614d-04, 1.88794337d-03/
      data fpppp(14, 9),fpppp(14,10)/ 2.57677382d-03, 1.12073965d-04/
      data fpppp(14,11),fpppp(14,12)/ 7.85346322d-03,-8.87524843d-03/
      data fpppp(14,13),fpppp(14,14)/-8.13720469d-03, 4.90128828d-03/
      data fpppp(14,15) /             1.11827300d-02 /
      data fpppp(15, 1),fpppp(15, 2)/-8.45747600d-04,-4.85130822d-04/
      data fpppp(15, 3),fpppp(15, 4)/ 9.87586011d-06, 3.45178804d-05/
      data fpppp(15, 5),fpppp(15, 6)/-1.47947382d-04, 9.68381148d-04/
      data fpppp(15, 7),fpppp(15, 8)/-9.49182182d-04, 5.54931779d-04/
      data fpppp(15, 9),fpppp(15,10)/ 8.10967426d-04, 2.64510539d-03/
      data fpppp(15,11),fpppp(15,12)/-8.78287664d-05,-2.45791838d-03/
      data fpppp(15,13),fpppp(15,14)/-2.53514428d-03, 2.21074842d-04/
      data fpppp(15,15) /             1.48671686d-03 /
      data fpppp(16, 1),fpppp(16, 2)/ 6.74844836d-04, 3.52658178d-04/
      data fpppp(16, 3),fpppp(16, 4)/-1.41791541d-04,-3.33205462d-05/
      data fpppp(16, 5),fpppp(16, 6)/ 2.75073726d-04,-8.19145825d-04/
      data fpppp(16, 7),fpppp(16, 8)/ 1.05782357d-03,-6.10721100d-04/
      data fpppp(16, 9),fpppp(16,10)/ 8.26937462d-04,-4.59768852d-04/
      data fpppp(16,11),fpppp(16,12)/-3.60635845d-04, 6.81460449d-05/
      data fpppp(16,13),fpppp(16,14)/ 6.33533926d-05,-3.41465236d-04/
      data fpppp(16,15) /            -5.31658636d-04 /
      data fpppp(17, 1),fpppp(17, 2)/-2.16356575d-04,-1.20052230d-04/
      data fpppp(17, 3),fpppp(17, 4)/ 1.82164886d-05, 9.60989939d-06/
      data fpppp(17, 5),fpppp(17, 6)/-5.66560862d-05, 2.54590822d-04/
      data fpppp(17, 7),fpppp(17, 8)/-2.83358193d-04, 5.86548373d-04/
      data fpppp(17, 9),fpppp(17,10)/-4.71854175d-04, 3.07921877d-04/
      data fpppp(17,11),fpppp(17,12)/-8.12298392d-04, 3.62064492d-04/
      data fpppp(17,13),fpppp(17,14)/ 2.82875111d-04,-4.95540870d-04/
      data fpppp(17,15) /            -8.79918830d-04 /
      data fpppp(18, 1),fpppp(18, 2)/-4.19278227d-05,-9.02213955d-06/
      data fpppp(18, 3),fpppp(18, 4)/ 6.12203964d-05, 7.83947415d-07/
      data fpppp(18, 5),fpppp(18, 6)/-6.43561861d-05, 1.99974033d-05/
      data fpppp(18, 7),fpppp(18, 8)/ 1.16255711d-06,-1.88480568d-04/
      data fpppp(18, 9),fpppp(18,10)/ 2.18878032d-04,-2.66822153d-04/
      data fpppp(18,11),fpppp(18,12)/ 2.52192654d-04,-7.24377259d-06/
      data fpppp(18,13),fpppp(18,14)/ 2.31400157d-05, 1.30657501d-04/
      data fpppp(18,15) /             1.88934671d-04 /
      data fpppp(19, 1),fpppp(19, 2)/ 1.55309532d-04, 7.36574548d-05/
      data fpppp(19, 3),fpppp(19, 4)/-6.44064072d-05,-5.02902413d-06/
      data fpppp(19, 5),fpppp(19, 6)/ 8.45225037d-05,-1.44063793d-04/
      data fpppp(19, 7),fpppp(19, 8)/ 1.06199722d-04,-5.31097727d-05/
      data fpppp(19, 9),fpppp(19,10)/-6.92150257d-05, 2.07869400d-06/
      data fpppp(19,11),fpppp(19,12)/-2.17629866d-05, 8.53616856d-05/
      data fpppp(19,13),fpppp(19,14)/ 8.18569453d-05,-7.74402566d-06/
      data fpppp(19,15) /            -5.04924095d-05 /
      data fpppp(20, 1),fpppp(20, 2)/-1.70354124d-04,-8.35200464d-05/
      data fpppp(20, 3),fpppp(20, 4)/ 5.90985165d-05, 6.47137997d-06/
      data fpppp(20, 5),fpppp(20, 6)/-8.49840364d-05, 1.74119366d-04/
      data fpppp(20, 7),fpppp(20, 8)/-1.66157635d-04, 1.03842818d-04/
      data fpppp(20, 9),fpppp(20,10)/-9.35143755d-05, 8.15700025d-05/
      data fpppp(20,11),fpppp(20,12)/-2.58947631d-05, 5.75062563d-06/
      data fpppp(20,13),fpppp(20,14)/-5.15553937d-08,-2.68603899d-06/
      data fpppp(20,15) /            -5.46271295d-06 /
      data fpppp(21, 1),fpppp(21, 2)/ 1.17005307d-04, 5.86260414d-05/
      data fpppp(21, 3),fpppp(21, 4)/-3.56992463d-05,-4.21345628d-06/
      data fpppp(21, 5),fpppp(21, 6)/ 5.25530714d-05,-1.17614429d-04/
      data fpppp(21, 7),fpppp(21, 8)/ 1.02094420d-04,-5.17151471d-05/
      data fpppp(21, 9),fpppp(21,10)/ 1.74235232d-05,-1.55090365d-05/
      data fpppp(21,11),fpppp(21,12)/ 1.97923738d-05, 9.84805332d-07/
      data fpppp(21,13),fpppp(21,14)/ 3.07712854d-06, 1.14230810d-05/
      data fpppp(21,15) /             1.58758116d-05 /
      data fpppp(22, 1),fpppp(22, 2)/ 2.29215284d-04, 1.15251104d-04/
      data fpppp(22, 3),fpppp(22, 4)/-6.81248129d-05,-8.55965243d-06/
      data fpppp(22, 5),fpppp(22, 6)/ 1.02363423d-04,-2.35086238d-04/
      data fpppp(22, 7),fpppp(22, 8)/ 2.15886642d-04,-1.17984382d-04/
      data fpppp(22, 9),fpppp(22,10)/ 8.97222080d-05,-6.21394047d-05/
      data fpppp(22,11),fpppp(22,12)/ 5.12455355d-05,-2.51653692d-05/
      data fpppp(22,13),fpppp(22,14)/-1.80905260d-05, 2.29461629d-05/
      data fpppp(22,15) /             4.39832425d-05 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /

      iprint=0
      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
      subroutine ch3po_2app(r,theta,phi,energy)
      implicit real*8 (a-h,o-z)
      data degrad / 57.29577951308232d 00 /
      call ch3po_2app_s(r,theta,es)
      call ch3po_2app_d(r,theta,ed)
      energy = es + 0.5d0 * cos(3.0d0*phi/degrad)*ed 
      return
      end

      subroutine ch3po_2app_s(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     2(2a") surface
c     average of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  5.26480000d+01 ,  4.78330000d+01 /
      data f( 1, 3),f( 1, 4) /  4.50910000d+01 ,  4.42300000d+01 /
      data f( 1, 5),f( 1, 6) /  4.50910000d+01 ,  4.78330000d+01 /
      data f( 1, 7),f( 1, 8) /  5.26480000d+01 ,  5.92970000d+01 /
      data f( 1, 9),f( 1,10) /  6.42360000d+01 ,  6.64340000d+01 /
      data f( 1,11),f( 1,12) /  7.14220000d+01 ,  8.44210000d+01 /
      data f( 1,13),f( 1,14) /  8.44210000d+01 ,  7.14220000d+01 /
      data f( 1,15) /           6.64340000d+01 /
      data f( 2, 1),f( 2, 2) /  4.58640000d+01 ,  4.16430000d+01 /
      data f( 2, 3),f( 2, 4) /  3.92720000d+01 ,  3.85300000d+01 /
      data f( 2, 5),f( 2, 6) /  3.92720000d+01 ,  4.16430000d+01 /
      data f( 2, 7),f( 2, 8) /  4.58640000d+01 ,  5.19240000d+01 /
      data f( 2, 9),f( 2,10) /  5.79670000d+01 ,  5.92950000d+01 /
      data f( 2,11),f( 2,12) /  6.28020000d+01 ,  7.15940000d+01 /
      data f( 2,13),f( 2,14) /  7.15940000d+01 ,  6.28020000d+01 /
      data f( 2,15) /           5.92950000d+01 /
      data f( 3, 1),f( 3, 2) /  3.98030000d+01 ,  3.61400000d+01 /
      data f( 3, 3),f( 3, 4) /  3.41050000d+01 ,  3.34730000d+01 /
      data f( 3, 5),f( 3, 6) /  3.41050000d+01 ,  3.61400000d+01 /
      data f( 3, 7),f( 3, 8) /  3.98030000d+01 ,  4.52130000d+01 /
      data f( 3, 9),f( 3,10) /  5.14560000d+01 ,  5.35500000d+01 /
      data f( 3,11),f( 3,12) /  5.57820000d+01 ,  6.15610000d+01 /
      data f( 3,13),f( 3,14) /  6.15610000d+01 ,  5.57820000d+01 /
      data f( 3,15) /           5.35500000d+01 /
      data f( 4, 1),f( 4, 2) /  3.44310000d+01 ,  3.12820000d+01 /
      data f( 4, 3),f( 4, 4) /  2.95510000d+01 ,  2.90160000d+01 /
      data f( 4, 5),f( 4, 6) /  2.95510000d+01 ,  3.12820000d+01 /
      data f( 4, 7),f( 4, 8) /  3.44310000d+01 ,  3.91810000d+01 /
      data f( 4, 9),f( 4,10) /  4.50950000d+01 ,  4.94030000d+01 /
      data f( 4,11),f( 4,12) /  4.96730000d+01 ,  5.38230000d+01 /
      data f( 4,13),f( 4,14) /  5.38230000d+01 ,  4.96730000d+01 /
      data f( 4,15) /           4.94030000d+01 /
      data f( 5, 1),f( 5, 2) /  2.97000000d+01 ,  2.70170000d+01 /
      data f( 5, 3),f( 5, 4) /  2.55560000d+01 ,  6.31140000d+01 /
      data f( 5, 5),f( 5, 6) /  2.55560000d+01 ,  2.70170000d+01 /
      data f( 5, 7),f( 5, 8) /  2.97000000d+01 ,  3.38120000d+01 /
      data f( 5, 9),f( 5,10) /  3.91660000d+01 ,  4.43400000d+01 /
      data f( 5,11),f( 5,12) /  4.58040000d+01 ,  4.75660000d+01 /
      data f( 5,13),f( 5,14) /  4.75660000d+01 ,  4.58040000d+01 /
      data f( 5,15) /           4.43400000d+01 /
      data f( 6, 1),f( 6, 2) /  2.19390000d+01 ,  2.00420000d+01 /
      data f( 6, 3),f( 6, 4) /  1.90290000d+01 ,  1.87240000d+01 /
      data f( 6, 5),f( 6, 6) /  1.90290000d+01 ,  2.00420000d+01 /
      data f( 6, 7),f( 6, 8) /  2.19390000d+01 ,  2.49220000d+01 /
      data f( 6, 9),f( 6,10) /  2.89900000d+01 ,  3.35740000d+01 /
      data f( 6,11),f( 6,12) /  3.73710000d+01 ,  3.96510000d+01 /
      data f( 6,13),f( 6,14) /  3.96510000d+01 ,  3.73710000d+01 /
      data f( 6,15) /           3.35740000d+01 /
      data f( 7, 1),f( 7, 2) /  1.60790000d+01 ,  1.47800000d+01 /
      data f( 7, 3),f( 7, 4) /  1.41010000d+01 ,  1.39020000d+01 /
      data f( 7, 5),f( 7, 6) /  1.41010000d+01 ,  1.47800000d+01 /
      data f( 7, 7),f( 7, 8) /  1.60790000d+01 ,  1.81660000d+01 /
      data f( 7, 9),f( 7,10) /  2.10810000d+01 ,  2.45030000d+01 /
      data f( 7,11),f( 7,12) /  2.75880000d+01 ,  2.93770000d+01 /
      data f( 7,13),f( 7,14) /  2.93770000d+01 ,  2.75880000d+01 /
      data f( 7,15) /           2.45030000d+01 /
      data f( 8, 1),f( 8, 2) /  1.17090000d+01 ,  1.08480000d+01 /
      data f( 8, 3),f( 8, 4) /  1.04090000d+01 ,  1.02850000d+01 /
      data f( 8, 5),f( 8, 6) /  1.04090000d+01 ,  1.08480000d+01 /
      data f( 8, 7),f( 8, 8) /  1.17090000d+01 ,  1.31230000d+01 /
      data f( 8, 9),f( 8,10) /  1.51300000d+01 ,  1.75210000d+01 /
      data f( 8,11),f( 8,12) /  1.97130000d+01 ,  2.09540000d+01 /
      data f( 8,13),f( 8,14) /  2.09540000d+01 ,  1.97130000d+01 /
      data f( 8,15) /           1.75210000d+01 /
      data f( 9, 1),f( 9, 2) /  8.47600000d+00 ,  7.92700000d+00 /
      data f( 9, 3),f( 9, 4) /  7.65600000d+00 ,  7.58000000d+00 /
      data f( 9, 5),f( 9, 6) /  7.65600000d+00 ,  7.92700000d+00 /
      data f( 9, 7),f( 9, 8) /  8.47600000d+00 ,  9.40400000d+00 /
      data f( 9, 9),f( 9,10) /  1.07390000d+01 ,  1.23380000d+01 /
      data f( 9,11),f( 9,12) /  1.37960000d+01 ,  1.45720000d+01 /
      data f( 9,13),f( 9,14) /  1.45720000d+01 ,  1.37960000d+01 /
      data f( 9,15) /           1.23380000d+01 /
      data f(10, 1),f(10, 2) /  6.10400000d+00 ,  5.76800000d+00 /
      data f(10, 3),f(10, 4) /  5.61000000d+00 ,  5.56600000d+00 /
      data f(10, 5),f(10, 6) /  5.61000000d+00 ,  5.76800000d+00 /
      data f(10, 7),f(10, 8) /  6.10400000d+00 ,  6.68800000d+00 /
      data f(10, 9),f(10,10) /  7.54400000d+00 ,  8.56900000d+00 /
      data f(10,11),f(10,12) /  9.48200000d+00 ,  9.91200000d+00 /
      data f(10,13),f(10,14) /  9.91200000d+00 ,  9.48200000d+00 /
      data f(10,15) /           8.56900000d+00 /
      data f(11, 1),f(11, 2) /  4.37300000d+00 ,  4.18100000d+00 /
      data f(11, 3),f(11, 4) /  4.09500000d+00 ,  4.07200000d+00 /
      data f(11, 5),f(11, 6) /  4.09500000d+00 ,  4.18100000d+00 /
      data f(11, 7),f(11, 8) /  4.37300000d+00 ,  4.72000000d+00 /
      data f(11, 9),f(11,10) /  5.24200000d+00 ,  5.86700000d+00 /
      data f(11,11),f(11,12) /  6.39900000d+00 ,  6.58800000d+00 /
      data f(11,13),f(11,14) /  6.58800000d+00 ,  6.39900000d+00 /
      data f(11,15) /           5.86700000d+00 /
      data f(12, 1),f(12, 2) /  2.62500000d+00 ,  2.56300000d+00 /
      data f(12, 3),f(12, 4) /  2.53800000d+00 ,  2.53300000d+00 /
      data f(12, 5),f(12, 6) /  2.53800000d+00 ,  2.56300000d+00 /
      data f(12, 7),f(12, 8) /  2.62500000d+00 ,  2.75400000d+00 /
      data f(12, 9),f(12,10) /  2.96600000d+00 ,  3.21900000d+00 /
      data f(12,11),f(12,12) /  3.40200000d+00 ,  3.37800000d+00 /
      data f(12,13),f(12,14) /  3.37800000d+00 ,  3.40200000d+00 /
      data f(12,15) /           3.21900000d+00 /
      data f(13, 1),f(13, 2) /  1.09000000d+00 ,  1.11400000d+00 /
      data f(13, 3),f(13, 4) /  1.13100000d+00 ,  1.13700000d+00 /
      data f(13, 5),f(13, 6) /  1.13100000d+00 ,  1.11400000d+00 /
      data f(13, 7),f(13, 8) /  1.09000000d+00 ,  1.06500000d+00 /
      data f(13, 9),f(13,10) /  1.04800000d+00 ,  1.03400000d+00 /
      data f(13,11),f(13,12) /  9.84000000d-01 ,  8.55000000d-01 /
      data f(13,13),f(13,14) /  8.55000000d-01 ,  9.84000000d-01 /
      data f(13,15) /           1.03400000d+00 /
      data f(14, 1),f(14, 2) /  4.23000000d-01 ,  4.66000000d-01 /
      data f(14, 3),f(14, 4) /  4.92000000d-01 ,  5.01000000d-01 /
      data f(14, 5),f(14, 6) /  4.92000000d-01 ,  4.66000000d-01 /
      data f(14, 7),f(14, 8) /  4.23000000d-01 ,  3.65000000d-01 /
      data f(14, 9),f(14,10) /  2.95000000d-01 ,  2.17000000d-01 /
      data f(14,11),f(14,12) /  1.30000000d-01 ,  3.90000000d-02 /
      data f(14,13),f(14,14) /  3.90000000d-02 ,  1.30000000d-01 /
      data f(14,15) /           2.17000000d-01 /
      data f(15, 1),f(15, 2) /  1.40000000d-01 ,  1.78000000d-01 /
      data f(15, 3),f(15, 4) /  2.02000000d-01 ,  2.10000000d-01 /
      data f(15, 5),f(15, 6) /  2.02000000d-01 ,  1.78000000d-01 /
      data f(15, 7),f(15, 8) /  1.40000000d-01 ,  8.60000000d-02 /
      data f(15, 9),f(15,10) /  2.40000000d-02 , -4.30000000d-02 /
      data f(15,11),f(15,12) / -1.04000000d-01 , -1.51000000d-01 /
      data f(15,13),f(15,14) / -1.51000000d-01 , -1.04000000d-01 /
      data f(15,15) /          -4.30000000d-02 /
      data f(16, 1),f(16, 2) /  2.60000000d-02 ,  5.50000000d-02 /
      data f(16, 3),f(16, 4) /  7.20000000d-02 ,  7.80000000d-02 /
      data f(16, 5),f(16, 6) /  7.20000000d-02 ,  5.50000000d-02 /
      data f(16, 7),f(16, 8) /  2.60000000d-02 , -1.40000000d-02 /
      data f(16, 9),f(16,10) / -5.90000000d-02 , -1.03000000d-01 /
      data f(16,11),f(16,12) / -1.37000000d-01 , -1.59000000d-01 /
      data f(16,13),f(16,14) / -1.59000000d-01 , -1.37000000d-01 /
      data f(16,15) /          -1.03000000d-01 /
      data f(17, 1),f(17, 2) / -1.50000000d-02 ,  5.00000000d-03 /
      data f(17, 3),f(17, 4) /  1.70000000d-02 ,  2.10000000d-02 /
      data f(17, 5),f(17, 6) /  1.70000000d-02 ,  5.00000000d-03 /
      data f(17, 7),f(17, 8) / -1.50000000d-02 , -4.20000000d-02 /
      data f(17, 9),f(17,10) / -7.40000000d-02 , -1.01000000d-01 /
      data f(17,11),f(17,12) / -1.20000000d-01 , -1.32000000d-01 /
      data f(17,13),f(17,14) / -1.32000000d-01 , -1.20000000d-01 /
      data f(17,15) /          -1.01000000d-01 /
      data f(18, 1),f(18, 2) / -2.70000000d-02 , -1.80000000d-02 /
      data f(18, 3),f(18, 4) / -1.30000000d-02 , -1.10000000d-02 /
      data f(18, 5),f(18, 6) / -1.30000000d-02 , -1.80000000d-02 /
      data f(18, 7),f(18, 8) / -2.70000000d-02 , -4.00000000d-02 /
      data f(18, 9),f(18,10) / -5.20000000d-02 , -6.40000000d-02 /
      data f(18,11),f(18,12) / -7.20000000d-02 , -7.80000000d-02 /
      data f(18,13),f(18,14) / -7.80000000d-02 , -7.20000000d-02 /
      data f(18,15) /          -6.40000000d-02 /
      data f(19, 1),f(19, 2) / -1.80000000d-02 , -1.50000000d-02 /
      data f(19, 3),f(19, 4) / -1.20000000d-02 , -1.20000000d-02 /
      data f(19, 5),f(19, 6) / -1.20000000d-02 , -1.50000000d-02 /
      data f(19, 7),f(19, 8) / -1.80000000d-02 , -2.40000000d-02 /
      data f(19, 9),f(19,10) / -2.90000000d-02 , -3.50000000d-02 /
      data f(19,11),f(19,12) / -3.90000000d-02 , -4.20000000d-02 /
      data f(19,13),f(19,14) / -4.20000000d-02 , -3.90000000d-02 /
      data f(19,15) /          -3.50000000d-02 /
      data f(20, 1),f(20, 2) / -9.00000000d-03 , -7.00000000d-03 /
      data f(20, 3),f(20, 4) / -6.00000000d-03 , -5.00000000d-03 /
      data f(20, 5),f(20, 6) / -6.00000000d-03 , -7.00000000d-03 /
      data f(20, 7),f(20, 8) / -9.00000000d-03 , -1.20000000d-02 /
      data f(20, 9),f(20,10) / -1.40000000d-02 , -1.80000000d-02 /
      data f(20,11),f(20,12) / -2.00000000d-02 , -2.10000000d-02 /
      data f(20,13),f(20,14) / -2.10000000d-02 , -2.00000000d-02 /
      data f(20,15) /          -1.80000000d-02 /
      data f(21, 1),f(21, 2) / -1.00000000d-03 ,  1.00000000d-03 /
      data f(21, 3),f(21, 4) /  1.00000000d-03 ,  2.00000000d-03 /
      data f(21, 5),f(21, 6) /  1.00000000d-03 ,  1.00000000d-03 /
      data f(21, 7),f(21, 8) / -1.00000000d-03 , -4.00000000d-03 /
      data f(21, 9),f(21,10) / -5.00000000d-03 , -7.00000000d-03 /
      data f(21,11),f(21,12) / -9.00000000d-03 , -1.10000000d-02 /
      data f(21,13),f(21,14) / -1.10000000d-02 , -9.00000000d-03 /
      data f(21,15) /          -7.00000000d-03 /
      data f(22, 1),f(22, 2) /  4.00000000d-03 ,  6.00000000d-03 /
      data f(22, 3),f(22, 4) /  7.00000000d-03 ,  7.00000000d-03 /
      data f(22, 5),f(22, 6) /  7.00000000d-03 ,  6.00000000d-03 /
      data f(22, 7),f(22, 8) /  4.00000000d-03 ,  2.00000000d-03 /
      data f(22, 9),f(22,10) /  1.00000000d-03 , -2.00000000d-03 /
      data f(22,11),f(22,12) / -4.00000000d-03 , -4.00000000d-03 /
      data f(22,13),f(22,14) / -4.00000000d-03 , -4.00000000d-03 /
      data f(22,15) /          -2.00000000d-03 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 7.57500125d+01, 2.26749604d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 7.29367888d+01, 2.06800791d-02/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 6.91633416d+01, 1.89847231d-02/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-2.52512424d+02, 1.62410284d-02/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 6.91633416d+01, 1.93711634d-02/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 7.29367888d+01, 1.91343179d-02/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 7.57500125d+01, 2.84715651d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 6.46238150d+01,-2.29805784d-02/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-6.30961178d+01,-3.91492515d-02/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.33568404d+02, 1.51175845d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 2.16100668d+02, 1.46078913d-01/
      data fpp( 1,12,1),fpp( 1,12,2)/ 3.32760273d+02,-1.18773238d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.32760273d+02,-1.06689742d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.16100668d+02, 9.77449262d-02/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.33568404d+02, 1.96370037d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 7.21999750d+01, 2.07071441d-02/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 6.86264223d+01, 1.85057118d-02/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 6.50733167d+01, 1.62700088d-02/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 7.06524848d+02, 1.41542529d-02/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 6.50733167d+01, 1.61529794d-02/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 6.86264223d+01, 1.89738294d-02/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 7.21999750d+01, 1.89517032d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 6.59523700d+01, 1.55593580d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-2.48077644d+01,-8.22091351d-02/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.10263192d+02, 3.03771825d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.85598664d+02, 9.14404051d-02/
      data fpp( 2,12,1),fpp( 2,12,2)/ 2.72479453d+02,-7.90388028d-02/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.72479453d+02,-7.23637940d-02/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.85598664d+02, 6.47403697d-02/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.10263192d+02, 1.30502315d-01/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 6.92500873d+01, 1.85203702d-02/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 6.47575219d+01, 1.62992596d-02/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 6.17433915d+01, 1.39625915d-02/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-2.18778697d+03, 1.20303744d-02/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 6.17433915d+01, 1.37559107d-02/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 6.47575219d+01, 1.71259826d-02/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 6.92500873d+01, 1.54201588d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 6.87667050d+01, 2.60133821d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 1.71271755d+01,-6.94936874d-02/
      data fpp( 3,10,1),fpp( 3,10,2)/ 2.61778829d+02, 3.02136735d-03/
      data fpp( 3,11,1),fpp( 3,11,2)/ 1.50467717d+00, 6.56882180d-02/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.53721914d+02,-5.29542392d-02/
      data fpp( 3,13,1),fpp( 3,13,2)/ 2.53721914d+02,-4.73513913d-02/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.50467717d+00, 4.32768261d-02/
      data fpp( 3,15,1),fpp( 3,15,2)/ 2.61778829d+02, 8.70640870d-02/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 6.41996756d+01, 1.63855550d-02/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 5.93434900d+01, 1.42088900d-02/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 5.57531172d+01, 1.18588851d-02/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 8.40462302d+03, 1.01155697d-02/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 5.57531172d+01, 1.18788361d-02/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 5.93434900d+01, 1.41290858d-02/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 6.41996756d+01, 1.66848208d-02/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 6.63808099d+01, 1.51916312d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 4.62990624d+01,-7.61134547d-03/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.98578507d+02,-8.11062493d-02/
      data fpp( 4,11,1),fpp( 4,11,2)/ 3.54982628d+02, 8.97563426d-02/
      data fpp( 4,12,1),fpp( 4,12,2)/ 8.96328899d+01,-4.51191213d-02/
      data fpp( 4,13,1),fpp( 4,13,2)/ 8.96328899d+01,-3.40208076d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.54982628d+02, 4.53630879d-02/
      data fpp( 4,15,1),fpp( 4,15,2)/-1.98578507d+02, 8.53684561d-02/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 5.85512103d+01,-2.72591821d-01/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 5.36685182d+01,-1.74096358d-01/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 5.06441397d+01, 1.04229725d+00/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-8.29770512d+03,-1.65395265d+00/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 5.06441397d+01, 1.06655335d+00/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 5.36685182d+01,-2.71120765d-01/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 5.85512103d+01, 9.12497060d-02/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 6.35100555d+01,-8.13805900d-03/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 5.68765748d+01, 1.58225300d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/-1.70647990d+01,-6.59520610d-02/
      data fpp( 5,11,1),fpp( 5,11,2)/-7.74351878d+01, 2.53857140d-02/
      data fpp( 5,12,1),fpp( 5,12,2)/ 2.76346526d+02,-1.77107951d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 2.76346526d+02,-1.24204717d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/-7.74351878d+01, 4.22442049d-03/
      data fpp( 5,15,1),fpp( 5,15,2)/-1.70647990d+01, 1.34027898d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 4.73965312d+01, 1.05875601d-02/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 4.25727003d+01, 8.86487988d-03/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 3.96410222d+01, 6.99292042d-03/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 3.80290384d+03, 5.64343844d-03/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 3.96410222d+01, 7.03332584d-03/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 4.25727003d+01, 8.70325822d-03/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 4.73965312d+01, 1.11936413d-02/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 5.34794285d+01, 1.16821766d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 5.85207444d+01, 7.17765229d-03/
      data fpp( 6,10,1),fpp( 6,10,2)/ 5.44836507d+01,-9.43278578d-03/
      data fpp( 6,11,1),fpp( 6,11,2)/-4.94357505d+01,-1.66665092d-02/
      data fpp( 6,12,1),fpp( 6,12,2)/-1.84006024d+02,-1.49211775d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.84006024d+02,-1.53032129d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/-4.94357505d+01,-1.51383678d-02/
      data fpp( 6,15,1),fpp( 6,15,2)/ 5.44836507d+01,-1.51633161d-02/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 3.70126649d+01, 7.59280128d-03/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 3.29906805d+01, 6.21439744d-03/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 3.06417713d+01, 4.74960896d-03/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-9.78710243d+02, 3.58716674d-03/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 3.06417713d+01, 4.78172409d-03/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 3.29906805d+01, 6.08593690d-03/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 3.70126649d+01, 8.07452832d-03/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 4.26722306d+01, 8.89594981d-03/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 4.90904478d+01, 6.02167242d-03/
      data fpp( 7,10,1),fpp( 7,10,2)/ 5.33801963d+01,-2.56263951d-03/
      data fpp( 7,11,1),fpp( 7,11,2)/ 7.26781899d+01,-1.59911144d-02/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.05827568d+02,-1.12329029d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.05827568d+02,-1.19757340d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/ 7.26781899d+01,-1.30197903d-02/
      data fpp( 7,15,1),fpp( 7,15,2)/ 5.33801963d+01,-1.37051049d-02/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 2.80528091d+01, 5.28739630d-03/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 2.49645778d+01, 4.22520740d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 2.31918926d+01, 3.13177408d-03/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 2.92687131d+02, 2.14769626d-03/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 2.31918926d+01, 3.15744086d-03/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 2.49645778d+01, 4.12254029d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 2.80528091d+01, 5.67239800d-03/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 3.27816492d+01, 6.36786773d-03/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 3.88174646d+01, 4.43613107d-03/
      data fpp( 8,10,1),fpp( 8,10,2)/ 4.53455640d+01,-1.07239202d-03/
      data fpp( 8,11,1),fpp( 8,11,2)/ 4.49229910d+01,-1.20865630d-02/
      data fpp( 8,12,1),fpp( 8,12,2)/ 3.83457514d+01,-7.64135604d-03/
      data fpp( 8,13,1),fpp( 8,13,2)/ 3.83457514d+01,-8.26265039d-03/
      data fpp( 8,14,1),fpp( 8,14,2)/ 4.49229910d+01,-9.60138560d-03/
      data fpp( 8,15,1),fpp( 8,15,2)/ 4.53455640d+01,-1.03918072d-02/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 2.13260987d+01, 3.60200924d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.88010082d+01, 2.79598152d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.74406583d+01, 1.89406468d-03/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-5.52382823d+01, 1.32775975d-03/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.74406583d+01, 1.91489631d-03/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.88010082d+01, 2.71265502d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.13260987d+01, 3.91448362d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 2.48011725d+01, 4.36941049d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 2.96396938d+01, 3.02787441d-03/
      data fpp( 9,10,1),fpp( 9,10,2)/ 3.50875475d+01,-6.40908119d-04/
      data fpp( 9,11,1),fpp( 9,11,2)/ 4.13298462d+01,-8.92424193d-03/
      data fpp( 9,12,1),fpp( 9,12,2)/ 4.69394262d+01,-4.58212416d-03/
      data fpp( 9,13,1),fpp( 9,13,2)/ 4.69394262d+01,-5.07150654d-03/
      data fpp( 9,14,1),fpp( 9,14,2)/ 4.13298462d+01,-6.96671242d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 3.50875475d+01,-7.98164379d-03/
      data fpp(10, 1,1),fpp(10, 1,2)/ 1.57927960d+01, 2.40898368d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.41313895d+01, 1.80203264d-03/
      data fpp(10, 3,1),fpp(10, 3,2)/ 1.30954741d+01, 1.06288575d-03/
      data fpp(10, 4,1),fpp(10, 4,2)/ 3.19159978d+01, 7.86424361d-04/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.30954741d+01, 1.07141681d-03/
      data fpp(10, 6,1),fpp(10, 6,2)/ 1.41313895d+01, 1.76790841d-03/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.57927960d+01, 2.53694955d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.84636607d+01, 2.96429340d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/ 2.20237602d+01, 1.92587686d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 2.64042460d+01,-5.27800828d-04/
      data fpp(10,11,1),fpp(10,11,2)/ 3.02076241d+01,-6.53467355d-03/
      data fpp(10,12,1),fpp(10,12,2)/ 3.21965439d+01,-2.31350499d-03/
      data fpp(10,13,1),fpp(10,13,2)/ 3.21965439d+01,-2.69214825d-03/
      data fpp(10,14,1),fpp(10,14,2)/ 3.02076241d+01,-5.02010050d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 2.64042460d+01,-6.20744975d-03/
      data fpp(11, 1,1),fpp(11, 1,2)/ 1.16527174d+01, 1.48229374d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/ 1.04734340d+01, 1.07541252d-03/
      data fpp(11, 3,1),fpp(11, 3,2)/ 9.82744530d+00, 5.76056183d-04/
      data fpp(11, 4,1),fpp(11, 4,2)/ 5.57429100d+00, 4.00362749d-04/
      data fpp(11, 5,1),fpp(11, 5,2)/ 9.82744530d+00, 5.82492822d-04/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.04734340d+01, 1.04966596d-03/
      data fpp(11, 7,1),fpp(11, 7,2)/ 1.16527174d+01, 1.57884333d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.35441847d+01, 1.93496073d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/ 1.62152654d+01, 1.18131375d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 1.93454685d+01,-4.80215732d-04/
      data fpp(11,11,1),fpp(11,11,2)/ 2.24896574d+01,-4.84045082d-03/
      data fpp(11,12,1),fpp(11,12,2)/ 2.46743984d+01,-7.37980977d-04/
      data fpp(11,13,1),fpp(11,13,2)/ 2.46743984d+01,-1.03583166d-03/
      data fpp(11,14,1),fpp(11,14,2)/ 2.24896574d+01,-3.64904810d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 1.93454685d+01,-4.94797595d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/ 7.19574455d+00, 5.39217173d-04/
      data fpp(12, 2,1),fpp(12, 2,2)/ 6.50096050d+00, 3.71565654d-04/
      data fpp(12, 3,1),fpp(12, 3,2)/ 6.21153295d+00, 1.94520211d-04/
      data fpp(12, 4,1),fpp(12, 4,2)/ 6.94169813d+00, 5.03535014d-05/
      data fpp(12, 5,1),fpp(12, 5,2)/ 6.21153295d+00, 2.04065783d-04/
      data fpp(12, 6,1),fpp(12, 6,2)/ 6.50096050d+00, 3.33383365d-04/
      data fpp(12, 7,1),fpp(12, 7,2)/ 7.19574455d+00, 6.82400755d-04/
      data fpp(12, 8,1),fpp(12, 8,2)/ 8.27694389d+00, 9.57013613d-04/
      data fpp(12, 9,1),fpp(12, 9,2)/ 9.73327513d+00, 4.69544791d-04/
      data fpp(12,10,1),fpp(12,10,2)/ 1.15789409d+01,-3.75192779d-04/
      data fpp(12,11,1),fpp(12,11,2)/ 1.33960594d+01,-3.16877368d-03/
      data fpp(12,12,1),fpp(12,12,2)/ 1.46876427d+01, 6.30287486d-04/
      data fpp(12,13,1),fpp(12,13,2)/ 1.46876427d+01, 4.13524380d-04/
      data fpp(12,14,1),fpp(12,14,2)/ 1.33960594d+01,-2.30172125d-03/
      data fpp(12,15,1),fpp(12,15,2)/ 1.15789409d+01,-3.62663937d-03/
      data fpp(13, 1,1),fpp(13, 1,2)/ 3.06198698d+00,-3.11065748d-05/
      data fpp(13, 2,1),fpp(13, 2,2)/ 2.85686603d+00,-6.77868505d-05/
      data fpp(13, 3,1),fpp(13, 3,2)/ 2.73862740d+00,-1.17746023d-04/
      data fpp(13, 4,1),fpp(13, 4,2)/ 2.49799139d+00,-1.21229056d-04/
      data fpp(13, 5,1),fpp(13, 5,2)/ 2.73862740d+00,-1.17337751d-04/
      data fpp(13, 6,1),fpp(13, 6,2)/ 2.85686603d+00,-6.94199378d-05/
      data fpp(13, 7,1),fpp(13, 7,2)/ 3.06198698d+00,-2.49824972d-05/
      data fpp(13, 8,1),fpp(13, 8,2)/ 3.49126874d+00, 1.09349927d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/ 4.13236034d+00, 6.75827912d-05/
      data fpp(13,10,1),fpp(13,10,2)/ 4.82010810d+00,-1.99681091d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 5.48681547d+00,-1.42885843d-03/
      data fpp(13,12,1),fpp(13,12,2)/ 6.04290422d+00, 1.17511480d-03/
      data fpp(13,13,1),fpp(13,13,2)/ 6.04290422d+00, 1.05908482d-03/
      data fpp(13,14,1),fpp(13,14,2)/ 5.48681547d+00,-9.64738520d-04/
      data fpp(13,15,1),fpp(13,15,2)/ 4.82010810d+00,-1.94013074d-03/
      data fpp(14, 1,1),fpp(14, 1,2)/ 1.38830752d+00,-1.69318907d-04/
      data fpp(14, 2,1),fpp(14, 2,2)/ 1.29557540d+00,-1.71362187d-04/
      data fpp(14, 3,1),fpp(14, 3,2)/ 1.26595747d+00,-1.65232347d-04/
      data fpp(14, 4,1),fpp(14, 4,2)/ 1.30633632d+00,-1.87708425d-04/
      data fpp(14, 5,1),fpp(14, 5,2)/ 1.26595747d+00,-1.63933954d-04/
      data fpp(14, 6,1),fpp(14, 6,2)/ 1.29557540d+00,-1.76555759d-04/
      data fpp(14, 7,1),fpp(14, 7,2)/ 1.38830752d+00,-1.49843008d-04/
      data fpp(14, 8,1),fpp(14, 8,2)/ 1.49398114d+00,-1.24072207d-04/
      data fpp(14, 9,1),fpp(14, 9,2)/ 1.69728353d+00,-7.38681631d-05/
      data fpp(14,10,1),fpp(14,10,2)/ 1.97262672d+00,-6.04551405d-05/
      data fpp(14,11,1),fpp(14,11,2)/ 2.19267871d+00,-2.24311275d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 2.10874039d+00, 7.17700241d-04/
      data fpp(14,13,1),fpp(14,13,2)/ 2.10874039d+00, 6.89054916d-04/
      data fpp(14,14,1),fpp(14,14,2)/ 2.19267871d+00,-1.09729976d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 1.97262672d+00,-4.90135012d-04/
      data fpp(15, 1,1),fpp(15, 1,2)/ 6.00782944d-01,-1.20946161d-04/
      data fpp(15, 2,1),fpp(15, 2,2)/ 6.00832392d-01,-1.38107678d-04/
      data fpp(15, 3,1),fpp(15, 3,2)/ 5.73542737d-01,-1.66623125d-04/
      data fpp(15, 4,1),fpp(15, 4,2)/ 5.56663341d-01,-1.55399820d-04/
      data fpp(15, 5,1),fpp(15, 5,2)/ 5.73542737d-01,-1.71777593d-04/
      data fpp(15, 6,1),fpp(15, 6,2)/ 6.00832392d-01,-1.17489808d-04/
      data fpp(15, 7,1),fpp(15, 7,2)/ 6.00782944d-01,-1.98263177d-04/
      data fpp(15, 8,1),fpp(15, 8,2)/ 6.36806689d-01,-4.94574848d-05/
      data fpp(15, 9,1),fpp(15, 9,2)/ 6.46505540d-01,-8.39068840d-05/
      data fpp(15,10,1),fpp(15,10,2)/ 6.57385029d-01, 8.50850206d-05/
      data fpp(15,11,1),fpp(15,11,2)/ 6.22469696d-01, 1.03566802d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 5.46134224d-01, 3.40647773d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 5.46134224d-01, 3.36273279d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 6.22469696d-01, 1.21064777d-04/
      data fpp(15,15,1),fpp(15,15,2)/ 6.57385029d-01, 1.94676113d-05/
      data fpp(16, 1,1),fpp(16, 1,2)/ 2.64560707d-01,-1.28884793d-04/
      data fpp(16, 2,1),fpp(16, 2,2)/ 2.61095037d-01,-1.22230414d-04/
      data fpp(16, 3,1),fpp(16, 3,2)/ 2.79871585d-01,-1.02193551d-04/
      data fpp(16, 4,1),fpp(16, 4,2)/ 2.83010319d-01,-1.28995383d-04/
      data fpp(16, 5,1),fpp(16, 5,2)/ 2.79871585d-01,-1.01824917d-04/
      data fpp(16, 6,1),fpp(16, 6,2)/ 2.61095037d-01,-1.23704949d-04/
      data fpp(16, 7,1),fpp(16, 7,2)/ 2.64560707d-01,-1.23355288d-04/
      data fpp(16, 8,1),fpp(16, 8,2)/ 2.54792102d-01,-4.28738992d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.28694310d-01,-5.14911500d-06/
      data fpp(16,10,1),fpp(16,10,2)/ 1.97833167d-01, 1.23470359d-04/
      data fpp(16,11,1),fpp(16,11,2)/ 1.41442509d-01, 1.11267678d-04/
      data fpp(16,12,1),fpp(16,12,2)/ 7.47227159d-02, 1.51458929d-04/
      data fpp(16,13,1),fpp(16,13,2)/ 7.47227159d-02, 1.49989375d-04/
      data fpp(16,14,1),fpp(16,14,2)/ 1.41442509d-01, 1.17145893d-04/
      data fpp(16,15,1),fpp(16,15,2)/ 1.97833167d-01, 1.01427054d-04/
      data fpp(17, 1,1),fpp(17, 1,2)/ 9.29742286d-02,-7.99886952d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/ 1.06787460d-01,-8.00226096d-05/
      data fpp(17, 3,1),fpp(17, 3,2)/ 1.06970923d-01,-7.99208664d-05/
      data fpp(17, 4,1),fpp(17, 4,2)/ 1.11295381d-01,-8.02939247d-05/
      data fpp(17, 5,1),fpp(17, 5,2)/ 1.06970923d-01,-7.89034348d-05/
      data fpp(17, 6,1),fpp(17, 6,2)/ 1.06787460d-01,-8.40923362d-05/
      data fpp(17, 7,1),fpp(17, 7,2)/ 9.29742286d-02,-6.47272204d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/ 7.20249035d-02,-7.69987823d-05/
      data fpp(17, 9,1),fpp(17, 9,2)/ 7.07172192d-02, 7.27223494d-05/
      data fpp(17,10,1),fpp(17,10,2)/ 3.92823021d-02, 8.61093847d-05/
      data fpp(17,11,1),fpp(17,11,2)/ 1.17602689d-02, 6.28401119d-05/
      data fpp(17,12,1),fpp(17,12,2)/-5.02508755d-03, 8.25301676d-05/
      data fpp(17,13,1),fpp(17,13,2)/-5.02508755d-03, 8.09894414d-05/
      data fpp(17,14,1),fpp(17,14,2)/ 1.17602689d-02, 6.90030168d-05/
      data fpp(17,15,1),fpp(17,15,2)/ 3.92823021d-02, 6.29984916d-05/
      data fpp(18, 1,1),fpp(18, 1,2)/ 8.79696086d-03,-4.88965034d-05/
      data fpp(18, 2,1),fpp(18, 2,2)/ 1.10901004d-02,-4.22069931d-05/
      data fpp(18, 3,1),fpp(18, 3,2)/ 1.91514377d-02,-2.22755241d-05/
      data fpp(18, 4,1),fpp(18, 4,2)/ 1.66086962d-02,-4.86909106d-05/
      data fpp(18, 5,1),fpp(18, 5,2)/ 1.91514377d-02,-2.29608334d-05/
      data fpp(18, 6,1),fpp(18, 6,2)/ 1.10901004d-02,-3.94657557d-05/
      data fpp(18, 7,1),fpp(18, 7,2)/ 8.79696086d-03,-5.91761439d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/ 4.52923857d-03, 3.61703312d-05/
      data fpp(18, 9,1),fpp(18, 9,2)/-1.44988129d-02,-2.55051810d-05/
      data fpp(18,10,1),fpp(18,10,2)/-1.87634898d-02, 6.58503928d-05/
      data fpp(18,11,1),fpp(18,11,2)/-2.20020610d-02, 2.10360981d-06/
      data fpp(18,12,1),fpp(18,12,2)/-2.22860953d-02, 4.57351680d-05/
      data fpp(18,13,1),fpp(18,13,2)/-2.22860953d-02, 4.17426912d-05/
      data fpp(18,14,1),fpp(18,14,2)/-2.20020610d-02, 1.80735168d-05/
      data fpp(18,15,1),fpp(18,15,2)/-1.87634898d-02, 5.96324160d-06/
      data fpp(19, 1,1),fpp(19, 1,2)/-2.16207198d-03, 2.66007071d-05/
      data fpp(19, 2,1),fpp(19, 2,2)/ 4.85213804d-03, 6.79858577d-06/
      data fpp(19, 3,1),fpp(19, 3,2)/ 2.42332598d-03,-5.37950502d-05/
      data fpp(19, 4,1),fpp(19, 4,2)/ 8.26983402d-03, 2.83816150d-05/
      data fpp(19, 5,1),fpp(19, 5,2)/ 2.42332598d-03,-5.97314099d-05/
      data fpp(19, 6,1),fpp(19, 6,2)/ 4.85213804d-03, 3.05440245d-05/
      data fpp(19, 7,1),fpp(19, 7,2)/-2.16207198d-03,-6.24446880d-05/
      data fpp(19, 8,1),fpp(19, 8,2)/-6.14185776d-03, 3.92347277d-05/
      data fpp(19, 9,1),fpp(19, 9,2)/-6.72196768d-03,-3.44942226d-05/
      data fpp(19,10,1),fpp(19,10,2)/-1.22283429d-02, 3.87421629d-05/
      data fpp(19,11,1),fpp(19,11,2)/-1.37520249d-02,-4.74428874d-07/
      data fpp(19,12,1),fpp(19,12,2)/-1.38305311d-02, 2.31555526d-05/
      data fpp(19,13,1),fpp(19,13,2)/-1.38305311d-02, 2.07705566d-05/
      data fpp(19,14,1),fpp(19,14,2)/-1.37520249d-02, 9.06555526d-06/
      data fpp(19,15,1),fpp(19,15,2)/-1.22283429d-02, 2.96722237d-06/
      data fpp(20, 1,1),fpp(20, 1,2)/-1.48672928d-04,-1.82233276d-05/
      data fpp(20, 2,1),fpp(20, 2,2)/-4.98652549d-04,-1.35533448d-05/
      data fpp(20, 3,1),fpp(20, 3,2)/ 1.15525839d-03, 1.24367068d-05/
      data fpp(20, 4,1),fpp(20, 4,2)/-1.68803224d-03,-3.61934823d-05/
      data fpp(20, 5,1),fpp(20, 5,2)/ 1.15525839d-03, 1.23372223d-05/
      data fpp(20, 6,1),fpp(20, 6,2)/-4.98652549d-04,-1.31554069d-05/
      data fpp(20, 7,1),fpp(20, 7,2)/-1.48672928d-04,-1.97155945d-05/
      data fpp(20, 8,1),fpp(20, 8,2)/-3.96180753d-03, 3.20177851d-05/
      data fpp(20, 9,1),fpp(20, 9,2)/-6.61331639d-03,-4.83555457d-05/
      data fpp(20,10,1),fpp(20,10,2)/-4.32313844d-03, 4.14043977d-05/
      data fpp(20,11,1),fpp(20,11,2)/-6.98983945d-03, 2.73795507d-06/
      data fpp(20,12,1),fpp(20,12,2)/-1.23917801d-02, 7.64378206d-06/
      data fpp(20,13,1),fpp(20,13,2)/-1.23917801d-02, 5.69967628d-06/
      data fpp(20,14,1),fpp(20,14,2)/-6.98983945d-03, 1.05143782d-05/
      data fpp(20,15,1),fpp(20,15,2)/-4.32313844d-03, 1.22428109d-05/
      data fpp(21, 1,1),fpp(21, 1,2)/-3.24323631d-03,-4.66443818d-05/
      data fpp(21, 2,1),fpp(21, 2,2)/-2.85752784d-03,-2.67112363d-05/
      data fpp(21, 3,1),fpp(21, 3,2)/-1.04435954d-03, 3.34893271d-05/
      data fpp(21, 4,1),fpp(21, 4,2)/-1.51770508d-03,-4.72460721d-05/
      data fpp(21, 5,1),fpp(21, 5,2)/-1.04435954d-03, 3.54949612d-05/
      data fpp(21, 6,1),fpp(21, 6,2)/-2.85752784d-03,-3.47337727d-05/
      data fpp(21, 7,1),fpp(21, 7,2)/-3.24323631d-03,-1.65598704d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/-2.01091214d-03, 4.09732542d-05/
      data fpp(21, 9,1),fpp(21, 9,2)/-2.82476674d-03,-2.73331465d-05/
      data fpp(21,10,1),fpp(21,10,2)/-6.47910330d-03, 8.35933192d-06/
      data fpp(21,11,1),fpp(21,11,2)/-6.28861730d-03,-6.10418115d-06/
      data fpp(21,12,1),fpp(21,12,2)/-2.60234855d-03, 1.60573927d-05/
      data fpp(21,13,1),fpp(21,13,2)/-2.60234855d-03, 1.48799126d-05/
      data fpp(21,14,1),fpp(21,14,2)/-6.28861730d-03,-1.39426073d-06/
      data fpp(21,15,1),fpp(21,15,2)/-6.47910330d-03,-9.30286963d-06/
      data fpp(22, 1,1),fpp(22, 1,2)/-4.87838185d-03,-1.06617032d-05/
      data fpp(22, 2,1),fpp(22, 2,2)/-6.07123608d-03,-8.67659357d-06/
      data fpp(22, 3,1),fpp(22, 3,2)/-2.97782023d-03,-1.46319225d-05/
      data fpp(22, 4,1),fpp(22, 4,2)/-4.24114746d-03, 7.20428364d-06/
      data fpp(22, 5,1),fpp(22, 5,2)/-2.97782023d-03,-1.41852120d-05/
      data fpp(22, 6,1),fpp(22, 6,2)/-6.07123608d-03,-1.04634355d-05/
      data fpp(22, 7,1),fpp(22, 7,2)/-4.87838185d-03,-3.96104582d-06/
      data fpp(22, 8,1),fpp(22, 8,2)/ 5.45606768d-06, 2.63076188d-05/
      data fpp(22, 9,1),fpp(22, 9,2)/-8.76166278d-05,-4.12694295d-05/
      data fpp(22,10,1),fpp(22,10,2)/-5.76044835d-03, 1.87700990d-05/
      data fpp(22,11,1),fpp(22,11,2)/-3.85569135d-03, 2.61890334d-05/
      data fpp(22,12,1),fpp(22,12,2)/ 4.80117427d-03,-3.52623272d-06/
      data fpp(22,13,1),fpp(22,13,2)/ 4.80117427d-03,-2.51581855d-06/
      data fpp(22,14,1),fpp(22,14,2)/-3.85569135d-03, 2.21473767d-05/
      data fpp(22,15,1),fpp(22,15,2)/-5.76044835d-03, 3.39263116d-05/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 2.37644281d+00, 1.55714956d+00/
      data fpppp( 1, 3),fpppp( 1, 4)/-8.66265447d+00, 1.40193292d+01/
      data fpppp( 1, 5),fpppp( 1, 6)/-8.81357051d+00, 2.16081372d+00/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.12702213d-01,-3.44798784d+00/
      data fpppp( 1, 9),fpppp( 1,10)/ 6.68362502d+00,-3.82344494d+00/
      data fpppp( 1,11),fpppp( 1,12)/ 1.76221929d+00,-1.17779173d+00/
      data fpppp( 1,13),fpppp( 1,14)/-8.47522614d-01, 4.41142824d-01/
      data fpppp( 1,15) /             1.13059179d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/-4.87275616d+00,-3.15355806d+00/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.74882152d+01,-2.80990247d+01/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.79336999d+01,-4.93549654d+00/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.80951312d+00,-2.89182542d+00/
      data fpppp( 2, 9),fpppp( 2,10)/ 4.68703681d+00,-2.30645639d+00/
      data fpppp( 2,11),fpppp( 2,12)/ 9.54659707d-01,-8.19463386d-01/
      data fpppp( 2,13),fpppp( 2,14)/-6.25363376d-01, 1.78259667d-01/
      data fpppp( 2,15) /             6.05043757d-01 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 1.69832073d+01, 1.10378318d+01/
      data fpppp( 3, 3),fpppp( 3, 4)/-6.10458283d+01, 9.83545078d+01/
      data fpppp( 3, 5),fpppp( 3, 6)/-6.24285597d+01, 1.65687575d+01/
      data fpppp( 3, 7),fpppp( 3, 8)/-3.75776414d+00,-1.83625781d+00/
      data fpppp( 3, 9),fpppp( 3,10)/ 8.03342654d+00,-1.25199774d+01/
      data fpppp( 3,11),fpppp( 3,12)/ 1.17509346d+01,-3.73427789d+00/
      data fpppp( 3,13),fpppp( 3,14)/-2.23915077d+00, 5.77042619d+00/
      data fpppp( 3,15) /             9.90692934d+00 /
      data fpppp( 4, 1),fpppp( 4, 2)/-6.30288977d+01,-4.09281195d+01/
      data fpppp( 4, 3),fpppp( 4, 4)/ 2.26817325d+02,-3.65193562d+02/
      data fpppp( 4, 5),fpppp( 4, 6)/ 2.32092536d+02,-6.20289635d+01/
      data fpppp( 4, 7),fpppp( 4, 8)/ 1.60992672d+01,-2.52860841d+00/
      data fpppp( 4, 9),fpppp( 4,10)/-7.32060649d+00, 1.83232850d+01/
      data fpppp( 4,11),fpppp( 4,12)/-1.80662112d+01, 4.80690757d+00/
      data fpppp( 4,13),fpppp( 4,14)/ 2.57287505d+00,-9.13008117d+00/
      data fpppp( 4,15) /            -1.51872027d+01 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 6.30268469d+01, 4.09457194d+01/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.26698226d+02, 3.65127691d+02/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.32010627d+02, 6.21953241d+01/
      data fpppp( 5, 7),fpppp( 5, 8)/-1.66591707d+01, 4.44592790d+00/
      data fpppp( 5, 9),fpppp( 5,10)/-1.82008045d+00,-1.20407970d+00/
      data fpppp( 5,11),fpppp( 5,12)/ 7.45065834d+00,-3.74942748d+00/
      data fpppp( 5,13),fpppp( 5,14)/-3.09049814d+00, 4.81494098d+00/
      data fpppp( 5,15) /             8.67986039d+00 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.83959719d+01,-1.84373385d+01/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.02258855d+02,-1.64626412d+02/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.04655253d+02,-2.80229318d+01/
      data fpppp( 6, 7),fpppp( 6, 8)/ 7.55000323d+00,-2.10153710d+00/
      data fpppp( 6, 9),fpppp( 6,10)/ 7.93650284d-01,-1.61776861d+00/
      data fpppp( 6,11),fpppp( 6,12)/-3.15514288d-01, 1.04077345d+00/
      data fpppp( 6,13),fpppp( 6,14)/ 1.07254497d+00,-4.42600372d-01/
      data fpppp( 6,15) /            -1.14119580d+00 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 7.63091971d+00, 4.96187645d+00/
      data fpppp( 7, 3),fpppp( 7, 4)/-2.73780410d+01, 4.41301012d+01/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.80201223d+01, 7.53020156d+00/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.00029946d+00, 5.69251141d-01/
      data fpppp( 7, 9),fpppp( 7,10)/-2.31186015d-01, 2.27784800d-01/
      data fpppp( 7,11),fpppp( 7,12)/ 2.20541515d-01,-2.78867785d-01/
      data fpppp( 7,13),fpppp( 7,14)/-2.68148747d-01, 1.77665365d-01/
      data fpppp( 7,15) /             3.88570361d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.02430385d+00,-1.31097348d+00/
      data fpppp( 8, 3),fpppp( 8, 4)/ 7.34713053d+00,-1.18014732d+01/
      data fpppp( 8, 5),fpppp( 8, 6)/ 7.51933358d+00,-1.99978571d+00/
      data fpppp( 8, 7),fpppp( 8, 8)/ 5.58742000d-01,-1.36745759d-01/
      data fpppp( 8, 9),fpppp( 8,10)/ 6.66595521d-02,-1.00355406d-01/
      data fpppp( 8,11),fpppp( 8,12)/-8.22782792d-02, 6.01885360d-02/
      data fpppp( 8,13),fpppp( 8,14)/ 5.78907170d-02,-7.30870032d-02/
      data fpppp( 8,15) /            -1.34822691d-01 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 5.59298383d-01, 3.66012086d-01/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.95346228d+00, 3.16872160d+00/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.99995125d+00, 5.51967961d-01/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.38036148d-01, 5.71756221d-02/
      data fpppp( 9, 9),fpppp( 9,10)/-8.85949084d-03, 1.48222860d-02/
      data fpppp( 9,11),fpppp( 9,12)/-2.76295003d-03,-4.17336138d-02/
      data fpppp( 9,13),fpppp( 9,14)/-4.17050818d-02,-2.87707796d-03/
      data fpppp( 9,15) /             1.52502657d-02 /
      data fpppp(10, 1),fpppp(10, 2)/-1.39035456d-01,-8.77833137d-02/
      data fpppp(10, 3),fpppp(10, 4)/ 5.27698178d-01,-8.31623052d-01/
      data fpppp(10, 5),fpppp(10, 6)/ 5.40331184d-01,-1.38315338d-01/
      data fpppp(10, 7),fpppp(10, 8)/ 5.04596352d-02,-2.95570714d-03/
      data fpppp(10, 9),fpppp(10,10)/ 1.47172787d-02,-6.69022961d-03/
      data fpppp(10,11),fpppp(10,12)/-2.25828214d-02,-1.18459863d-02/
      data fpppp(10,13),fpppp(10,14)/-1.28382228d-02,-1.86138756d-02/
      data fpppp(10,15) /            -2.15737766d-02 /
      data fpppp(11, 1),fpppp(11, 2)/ 3.67937863d-02, 2.52204818d-02/
      data fpppp(11, 3),fpppp(11, 4)/-1.05678024d-01, 1.81061674d-01/
      data fpppp(11, 5),fpppp(11, 6)/-1.08190157d-01, 3.52690155d-02/
      data fpppp(11, 7),fpppp(11, 8)/-8.88215050d-04, 1.10148703d-02/
      data fpppp(11, 9),fpppp(11,10)/ 3.60554249d-03, 2.11030318d-03/
      data fpppp(11,11),fpppp(11,12)/-1.12076135d-02,-1.48467146d-02/
      data fpppp(11,13),fpppp(11,14)/-1.53982812d-02,-9.00134717d-03/
      data fpppp(11,15) /            -6.16319573d-03 /
      data fpppp(12, 1),fpppp(12, 2)/-1.58427397d-04, 1.92824596d-04/
      data fpppp(12, 3),fpppp(12, 4)/ 2.37085185d-02,-3.38513343d-02/
      data fpppp(12, 5),fpppp(12, 6)/ 2.40769966d-02,-1.28108779d-03/
      data fpppp(12, 7),fpppp(12, 8)/ 5.36874404d-03, 2.99102931d-03/
      data fpppp(12, 9),fpppp(12,10)/ 5.17505269d-03,-3.31170281d-04/
      data fpppp(12,11),fpppp(12,12)/-5.56320306d-03,-8.94813115d-03/
      data fpppp(12,13),fpppp(12,14)/-9.12150455d-03,-4.86970945d-03/
      data fpppp(12,15) /            -2.93177131d-03 /
      data fpppp(13, 1),fpppp(13, 2)/ 2.46598866d-03, 1.86008672d-03/
      data fpppp(13, 3),fpppp(13, 4)/-4.69339582d-03, 9.56965377d-03/
      data fpppp(13, 5),fpppp(13, 6)/-4.70889821d-03, 1.92209626d-03/
      data fpppp(13, 7),fpppp(13, 8)/ 2.23345290d-03, 2.59374022d-03/
      data fpppp(13, 9),fpppp(13,10)/ 1.00176260d-04,-1.95074644d-04/
      data fpppp(13,11),fpppp(13,12)/-5.82301908d-04,-4.11283463d-03/
      data fpppp(13,13),fpppp(13,14)/-4.05300767d-03,-8.21609758d-04/
      data fpppp(13,15) /             7.02329791d-04 /
      data fpppp(14, 1),fpppp(14, 2)/ 6.65146756d-04, 4.25480557d-04/
      data fpppp(14, 3),fpppp(14, 4)/ 1.41978263d-03,-1.90480424d-03/
      data fpppp(14, 5),fpppp(14, 6)/ 1.35397219d-03, 6.88722314d-04/
      data fpppp(14, 7),fpppp(14, 8)/-3.22009830d-04, 1.37580709d-03/
      data fpppp(14, 9),fpppp(14,10)/ 6.76507374d-04, 2.40611176d-04/
      data fpppp(14,11),fpppp(14,12)/-4.95642368d-03, 1.34566484d-03/
      data fpppp(14,13),fpppp(14,14)/ 9.59366918d-04,-3.41123200d-03/
      data fpppp(14,15) /            -5.55385761d-03 /
      data fpppp(15, 1),fpppp(15, 2)/-6.48533807d-04,-2.78092694d-04/
      data fpppp(15, 3),fpppp(15, 4)/ 1.20558410d-04, 4.20474549d-04/
      data fpppp(15, 5),fpppp(15, 6)/ 2.23070974d-04,-6.88142949d-04/
      data fpppp(15, 7),fpppp(15, 8)/ 8.89154649d-04,-7.04084044d-04/
      data fpppp(15, 9),fpppp(15,10)/ 3.47687846d-04,-6.15829008d-04/
      data fpppp(15,11),fpppp(15,12)/-6.32061183d-04, 6.58865434d-04/
      data fpppp(15,13),fpppp(15,14)/ 6.29498444d-04,-5.14593222d-04/
      data fpppp(15,15) /            -1.05633386d-03 /
      data fpppp(16, 1),fpppp(16, 2)/ 5.77515397d-04, 2.69836371d-04/
      data fpppp(16, 3),fpppp(16, 4)/-3.22327804d-04, 8.12060436d-05/
      data fpppp(16, 5),fpppp(16, 6)/-3.79144523d-04, 4.97103247d-04/
      data fpppp(16, 7),fpppp(16, 8)/-2.74735387d-04,-1.92218201d-04/
      data fpppp(16, 9),fpppp(16,10)/ 6.38569984d-05,-3.49010881d-04/
      data fpppp(16,11),fpppp(16,12)/-1.99584391d-04, 5.27600369d-04/
      data fpppp(16,13),fpppp(16,14)/ 5.18584876d-04,-1.63522417d-04/
      data fpppp(16,15) /            -4.84243285d-04 /
      data fpppp(17, 1),fpppp(17, 2)/-2.97494260d-04,-1.69319828d-04/
      data fpppp(17, 3),fpppp(17, 4)/ 1.56987439d-04,-2.10170219d-04/
      data fpppp(17, 5),fpppp(17, 6)/ 1.64758469d-04,-2.00403948d-04/
      data fpppp(17, 7),fpppp(17, 8)/-1.80928812d-04, 4.95953602d-04/
      data fpppp(17, 9),fpppp(17,10)/-6.24387149d-04, 1.93961016d-04/
      data fpppp(17,11),fpppp(17,12)/ 8.33161250d-05, 1.16975093d-04/
      data fpppp(17,13),fpppp(17,14)/ 1.10977352d-04, 1.07307089d-04/
      data fpppp(17,15) /             1.03994902d-04 /
      data fpppp(18, 1),fpppp(18, 2)/ 2.03525059d-04, 9.34413447d-05/
      data fpppp(18, 3),fpppp(18, 4)/-2.31198574d-04, 1.95108221d-04/
      data fpppp(18, 5),fpppp(18, 6)/-2.44105326d-04, 1.45068353d-04/
      data fpppp(18, 7),fpppp(18, 8)/ 9.92377797d-06,-3.03238430d-04/
      data fpppp(18, 9),fpppp(18,10)/ 3.17410192d-04,-8.05998630d-05/
      data fpppp(18,11),fpppp(18,12)/ 6.65556006d-05,-8.35032661d-06/
      data fpppp(18,13),fpppp(18,14)/ 2.94209604d-07, 3.19774557d-05/
      data fpppp(18,15) /             4.90681802d-05 /
      data fpppp(19, 1),fpppp(19, 2)/-2.50175865d-04,-1.37305776d-04/
      data fpppp(19, 3),fpppp(19, 4)/ 2.32817644d-04,-2.97445594d-04/
      data fpppp(19, 5),fpppp(19, 6)/ 2.55383767d-04,-2.27570267d-04/
      data fpppp(19, 7),fpppp(19, 8)/ 8.83159765d-05, 5.63718153d-05/
      data fpppp(19, 9),fpppp(19,10)/-1.09822686d-04, 8.73430081d-05/
      data fpppp(19,11),fpppp(19,12)/-5.87746628d-07, 1.71851938d-06/
      data fpppp(19,13),fpppp(19,14)/-2.50649702d-06, 1.63123190d-05/
      data fpppp(19,15) /             2.39677621d-05 /
      data fpppp(20, 1),fpppp(20, 2)/ 7.52509286d-05, 3.96367020d-05/
      data fpppp(20, 3),fpppp(20, 4)/-1.13564303d-04, 1.44788417d-04/
      data fpppp(20, 5),fpppp(20, 6)/-1.24394489d-04, 8.29574452d-05/
      data fpppp(20, 7),fpppp(20, 8)/-8.72018582d-05, 1.60631345d-05/
      data fpppp(20, 9),fpppp(20,10)/ 9.26468640d-05,-9.01493810d-05/
      data fpppp(20,11),fpppp(20,12)/-2.94620783d-05, 4.38833178d-05/
      data fpppp(20,13),fpppp(20,14)/ 4.51393044d-05,-3.44860247d-05/
      data fpppp(20,15) /            -7.13095819d-05 /
      data fpppp(21, 1),fpppp(21, 2)/ 4.75900334d-05, 2.19232020d-05/
      data fpppp(21, 3),fpppp(21, 4)/-4.96352509d-05, 3.94269714d-05/
      data fpppp(21, 5),fpppp(21, 6)/-5.12711705d-05, 2.84668802d-05/
      data fpppp(21, 7),fpppp(21, 8)/ 2.30512401d-05,-2.35898825d-05/
      data fpppp(21, 9),fpppp(21,10)/-5.14624369d-05, 5.90107130d-05/
      data fpppp(21,11),fpppp(21,12)/ 4.61089388d-05,-3.36995036d-05/
      data fpppp(21,13),fpppp(21,14)/-3.25440211d-05, 4.14870086d-05/
      data fpppp(21,15) /             7.63429512d-05 /
      data fpppp(22, 1),fpppp(22, 2)/ 1.20811335d-04, 5.98256956d-05/
      data fpppp(22, 3),fpppp(22, 4)/-1.02937913d-04, 9.05213706d-05/
      data fpppp(22, 5),fpppp(22, 6)/-1.07548302d-04, 7.82672515d-05/
      data fpppp(22, 7),fpppp(22, 8)/ 5.16555004d-05,-6.34302320d-05/
      data fpppp(22, 9),fpppp(22,10)/-9.65492089d-05, 1.14841526d-04/
      data fpppp(22,11),fpppp(22,12)/ 9.18384277d-05,-7.70687194d-05/
      data fpppp(22,13),fpppp(22,14)/-7.44190246d-05, 8.12396486d-05/
      data fpppp(22,15) /             1.54586948d-04 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /


      iprint=0
      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
      subroutine ch3po_2app_d(xi,yi,fi)
      implicit real*8 (a-h,o-z)
c
c     ch3+o
c     cas+1+2+qc/aug-cc-pvdz
c     2(2a") surface
c     difference of eclipsed and staggered energies
c
      dimension fpp(22,15,2),f(22,15),fpppp(22,15)
      dimension delx(21),dely(14),x(22),y(15)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -7.14300000d+00 , -2.15200000d+00 /
      data f( 1, 3),f( 1, 4) / -2.61000000d-01 ,  0.00000000d+00 /
      data f( 1, 5),f( 1, 6) /  2.61000000d-01 ,  2.15200000d+00 /
      data f( 1, 7),f( 1, 8) /  7.14300000d+00 ,  1.55450000d+01 /
      data f( 1, 9),f( 1,10) /  2.10130000d+01 ,  2.23340000d+01 /
      data f( 1,11),f( 1,12) /  3.13460000d+01 ,  5.88410000d+01 /
      data f( 1,13),f( 1,14) /  5.88410000d+01 ,  3.13460000d+01 /
      data f( 1,15) /           2.23340000d+01 /
      data f( 2, 1),f( 2, 2) / -6.37700000d+00 , -1.90600000d+00 /
      data f( 2, 3),f( 2, 4) / -2.33000000d-01 ,  0.00000000d+00 /
      data f( 2, 5),f( 2, 6) /  2.33000000d-01 ,  1.90600000d+00 /
      data f( 2, 7),f( 2, 8) /  6.37700000d+00 ,  1.42790000d+01 /
      data f( 2, 9),f( 2,10) /  2.25270000d+01 ,  2.24960000d+01 /
      data f( 2,11),f( 2,12) /  2.87700000d+01 ,  4.79790000d+01 /
      data f( 2,13),f( 2,14) /  4.79790000d+01 ,  2.87700000d+01 /
      data f( 2,15) /           2.24960000d+01 /
      data f( 3, 1),f( 3, 2) / -5.65100000d+00 , -1.68000000d+00 /
      data f( 3, 3),f( 3, 4) / -2.06000000d-01 ,  0.00000000d+00 /
      data f( 3, 5),f( 3, 6) /  2.06000000d-01 ,  1.68000000d+00 /
      data f( 3, 7),f( 3, 8) /  5.65100000d+00 ,  1.28990000d+01 /
      data f( 3, 9),f( 3,10) /  2.21250000d+01 ,  2.40470000d+01 /
      data f( 3,11),f( 3,12) /  2.79870000d+01 ,  4.11670000d+01 /
      data f( 3,13),f( 3,14) /  4.11670000d+01 ,  2.79870000d+01 /
      data f( 3,15) /           2.40470000d+01 /
      data f( 4, 1),f( 4, 2) / -4.97500000d+00 , -1.47300000d+00 /
      data f( 4, 3),f( 4, 4) / -1.82000000d-01 ,  0.00000000d+00 /
      data f( 4, 5),f( 4, 6) /  1.82000000d-01 ,  1.47300000d+00 /
      data f( 4, 7),f( 4, 8) /  4.97500000d+00 ,  1.14990000d+01 /
      data f( 4, 9),f( 4,10) /  2.06200000d+01 ,  2.73810000d+01 /
      data f( 4,11),f( 4,12) /  2.75980000d+01 ,  3.74540000d+01 /
      data f( 4,13),f( 4,14) /  3.74540000d+01 ,  2.75980000d+01 /
      data f( 4,15) /           2.73810000d+01 /
      data f( 5, 1),f( 5, 2) / -4.35100000d+00 , -1.28600000d+00 /
      data f( 5, 3),f( 5, 4) / -1.60000000d-01 ,  0.00000000d+00 /
      data f( 5, 5),f( 5, 6) /  1.60000000d-01 ,  1.28600000d+00 /
      data f( 5, 7),f( 5, 8) /  4.35100000d+00 ,  1.01390000d+01 /
      data f( 5, 9),f( 5,10) /  1.86490000d+01 ,  2.75270000d+01 /
      data f( 5,11),f( 5,12) /  3.03080000d+01 ,  3.52890000d+01 /
      data f( 5,13),f( 5,14) /  3.52890000d+01 ,  3.03080000d+01 /
      data f( 5,15) /           2.75270000d+01 /
      data f( 6, 1),f( 6, 2) / -3.27600000d+00 , -9.66000000d-01 /
      data f( 6, 3),f( 6, 4) / -1.22000000d-01 ,  0.00000000d+00 /
      data f( 6, 5),f( 6, 6) /  1.22000000d-01 ,  9.66000000d-01 /
      data f( 6, 7),f( 6, 8) /  3.27600000d+00 ,  7.69400000d+00 /
      data f( 6, 9),f( 6,10) /  1.44740000d+01 ,  2.28220000d+01 /
      data f( 6,11),f( 6,12) /  3.05470000d+01 ,  3.63600000d+01 /
      data f( 6,13),f( 6,14) /  3.63600000d+01 ,  3.05470000d+01 /
      data f( 6,15) /           2.28220000d+01 /
      data f( 7, 1),f( 7, 2) / -2.42100000d+00 , -7.14000000d-01 /
      data f( 7, 3),f( 7, 4) / -9.10000000d-02 ,  0.00000000d+00 /
      data f( 7, 5),f( 7, 6) /  9.10000000d-02 ,  7.14000000d-01 /
      data f( 7, 7),f( 7, 8) /  2.42100000d+00 ,  5.70300000d+00 /
      data f( 7, 9),f( 7,10) /  1.07950000d+01 ,  1.72940000d+01 /
      data f( 7,11),f( 7,12) /  2.37880000d+01 ,  2.84410000d+01 /
      data f( 7,13),f( 7,14) /  2.84410000d+01 ,  2.37880000d+01 /
      data f( 7,15) /           1.72940000d+01 /
      data f( 8, 1),f( 8, 2) / -1.76300000d+00 , -5.18000000d-01 /
      data f( 8, 3),f( 8, 4) / -6.50000000d-02 ,  0.00000000d+00 /
      data f( 8, 5),f( 8, 6) /  6.50000000d-02 ,  5.18000000d-01 /
      data f( 8, 7),f( 8, 8) /  1.76300000d+00 ,  4.15700000d+00 /
      data f( 8, 9),f( 8,10) /  7.87200000d+00 ,  1.26340000d+01 /
      data f( 8,11),f( 8,12) /  1.74540000d+01 ,  2.08630000d+01 /
      data f( 8,13),f( 8,14) /  2.08630000d+01 ,  1.74540000d+01 /
      data f( 8,15) /           1.26340000d+01 /
      data f( 9, 1),f( 9, 2) / -1.26700000d+00 , -3.69000000d-01 /
      data f( 9, 3),f( 9, 4) / -4.60000000d-02 ,  0.00000000d+00 /
      data f( 9, 5),f( 9, 6) /  4.60000000d-02 ,  3.69000000d-01 /
      data f( 9, 7),f( 9, 8) /  1.26700000d+00 ,  2.99200000d+00 /
      data f( 9, 9),f( 9,10) /  5.65300000d+00 ,  9.03900000d+00 /
      data f( 9,11),f( 9,12) /  1.24450000d+01 ,  1.48070000d+01 /
      data f( 9,13),f( 9,14) /  1.48070000d+01 ,  1.24450000d+01 /
      data f( 9,15) /           9.03900000d+00 /
      data f(10, 1),f(10, 2) / -9.01000000d-01 , -2.61000000d-01 /
      data f(10, 3),f(10, 4) / -3.20000000d-02 ,  0.00000000d+00 /
      data f(10, 5),f(10, 6) /  3.20000000d-02 ,  2.61000000d-01 /
      data f(10, 7),f(10, 8) /  9.01000000d-01 ,  2.13000000d+00 /
      data f(10, 9),f(10,10) /  4.01000000d+00 ,  6.37300000d+00 /
      data f(10,11),f(10,12) /  8.70900000d+00 ,  1.02710000d+01 /
      data f(10,13),f(10,14) /  1.02710000d+01 ,  8.70900000d+00 /
      data f(10,15) /           6.37300000d+00 /
      data f(11, 1),f(11, 2) / -6.35000000d-01 , -1.85000000d-01 /
      data f(11, 3),f(11, 4) / -2.30000000d-02 ,  0.00000000d+00 /
      data f(11, 5),f(11, 6) /  2.30000000d-02 ,  1.85000000d-01 /
      data f(11, 7),f(11, 8) /  6.35000000d-01 ,  1.50000000d+00 /
      data f(11, 9),f(11,10) /  2.81200000d+00 ,  4.43500000d+00 /
      data f(11,11),f(11,12) /  6.00100000d+00 ,  6.99100000d+00 /
      data f(11,13),f(11,14) /  6.99100000d+00 ,  6.00100000d+00 /
      data f(11,15) /           4.43500000d+00 /
      data f(12, 1),f(12, 2) / -3.72000000d-01 , -1.12000000d-01 /
      data f(12, 3),f(12, 4) / -1.60000000d-02 ,  0.00000000d+00 /
      data f(12, 5),f(12, 6) /  1.60000000d-02 ,  1.12000000d-01 /
      data f(12, 7),f(12, 8) /  3.72000000d-01 ,  8.68000000d-01 /
      data f(12, 9),f(12,10) /  1.60900000d+00 ,  2.50600000d+00 /
      data f(12,11),f(12,12) /  3.32900000d+00 ,  3.78300000d+00 /
      data f(12,13),f(12,14) /  3.78300000d+00 ,  3.32900000d+00 /
      data f(12,15) /           2.50600000d+00 /
      data f(13, 1),f(13, 2) / -1.52000000d-01 , -4.80000000d-02 /
      data f(13, 3),f(13, 4) / -7.00000000d-03 ,  0.00000000d+00 /
      data f(13, 5),f(13, 6) /  7.00000000d-03 ,  4.80000000d-02 /
      data f(13, 7),f(13, 8) /  1.52000000d-01 ,  3.29000000d-01 /
      data f(13, 9),f(13,10) /  5.80000000d-01 ,  8.71000000d-01 /
      data f(13,11),f(13,12) /  1.11300000d+00 ,  1.19700000d+00 /
      data f(13,13),f(13,14) /  1.19700000d+00 ,  1.11300000d+00 /
      data f(13,15) /           8.71000000d-01 /
      data f(14, 1),f(14, 2) / -5.60000000d-02 , -1.90000000d-02 /
      data f(14, 3),f(14, 4) / -2.00000000d-03 ,  0.00000000d+00 /
      data f(14, 5),f(14, 6) /  2.00000000d-03 ,  1.90000000d-02 /
      data f(14, 7),f(14, 8) /  5.60000000d-02 ,  1.14000000d-01 /
      data f(14, 9),f(14,10) /  1.82000000d-01 ,  2.51000000d-01 /
      data f(14,11),f(14,12) /  3.05000000d-01 ,  3.38000000d-01 /
      data f(14,13),f(14,14) /  3.38000000d-01 ,  3.05000000d-01 /
      data f(14,15) /           2.51000000d-01 /
      data f(15, 1),f(15, 2) / -1.50000000d-02 , -5.00000000d-03 /
      data f(15, 3),f(15, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(15, 5),f(15, 6) /  0.00000000d+00 ,  5.00000000d-03 /
      data f(15, 7),f(15, 8) /  1.50000000d-02 ,  2.90000000d-02 /
      data f(15, 9),f(15,10) /  4.40000000d-02 ,  5.80000000d-02 /
      data f(15,11),f(15,12) /  7.50000000d-02 ,  9.70000000d-02 /
      data f(15,13),f(15,14) /  9.70000000d-02 ,  7.50000000d-02 /
      data f(15,15) /           5.80000000d-02 /
      data f(16, 1),f(16, 2) /  1.00000000d-03 ,  0.00000000d+00 /
      data f(16, 3),f(16, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(16, 5),f(16, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(16, 7),f(16, 8) / -1.00000000d-03 , -2.00000000d-03 /
      data f(16, 9),f(16,10) / -1.00000000d-03 ,  5.00000000d-03 /
      data f(16,11),f(16,12) /  1.80000000d-02 ,  3.20000000d-02 /
      data f(16,13),f(16,14) /  3.20000000d-02 ,  1.80000000d-02 /
      data f(16,15) /           5.00000000d-03 /
      data f(17, 1),f(17, 2) /  4.00000000d-03 ,  1.00000000d-03 /
      data f(17, 3),f(17, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(17, 5),f(17, 6) /  0.00000000d+00 , -1.00000000d-03 /
      data f(17, 7),f(17, 8) / -4.00000000d-03 , -9.00000000d-03 /
      data f(17, 9),f(17,10) / -1.30000000d-02 , -9.00000000d-03 /
      data f(17,11),f(17,12) /  2.00000000d-03 ,  9.00000000d-03 /
      data f(17,13),f(17,14) /  9.00000000d-03 ,  2.00000000d-03 /
      data f(17,15) /          -9.00000000d-03 /
      data f(18, 1),f(18, 2) /  2.00000000d-03 ,  0.00000000d+00 /
      data f(18, 3),f(18, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(18, 5),f(18, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(18, 7),f(18, 8) / -2.00000000d-03 , -3.00000000d-03 /
      data f(18, 9),f(18,10) / -4.00000000d-03 , -5.00000000d-03 /
      data f(18,11),f(18,12) / -3.00000000d-03 , -2.00000000d-03 /
      data f(18,13),f(18,14) / -2.00000000d-03 , -3.00000000d-03 /
      data f(18,15) /          -5.00000000d-03 /
      data f(19, 1),f(19, 2) / -1.00000000d-03 ,  0.00000000d+00 /
      data f(19, 3),f(19, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(19, 5),f(19, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(19, 7),f(19, 8) /  1.00000000d-03 ,  1.00000000d-03 /
      data f(19, 9),f(19,10) /  3.00000000d-03 ,  3.00000000d-03 /
      data f(19,11),f(19,12) /  4.00000000d-03 ,  4.00000000d-03 /
      data f(19,13),f(19,14) /  4.00000000d-03 ,  4.00000000d-03 /
      data f(19,15) /           3.00000000d-03 /
      data f(20, 1),f(20, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(20, 3),f(20, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(20, 5),f(20, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(20, 7),f(20, 8) /  0.00000000d+00 ,  2.00000000d-03 /
      data f(20, 9),f(20,10) /  3.00000000d-03 ,  3.00000000d-03 /
      data f(20,11),f(20,12) /  4.00000000d-03 ,  3.00000000d-03 /
      data f(20,13),f(20,14) /  3.00000000d-03 ,  4.00000000d-03 /
      data f(20,15) /           3.00000000d-03 /
      data f(21, 1),f(21, 2) /  0.00000000d+00 , -1.00000000d-03 /
      data f(21, 3),f(21, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(21, 5),f(21, 6) /  0.00000000d+00 ,  1.00000000d-03 /
      data f(21, 7),f(21, 8) /  0.00000000d+00 ,  1.00000000d-03 /
      data f(21, 9),f(21,10) /  1.00000000d-03 ,  1.00000000d-03 /
      data f(21,11),f(21,12) / -1.00000000d-03 , -2.00000000d-03 /
      data f(21,13),f(21,14) / -2.00000000d-03 , -1.00000000d-03 /
      data f(21,15) /           1.00000000d-03 /
      data f(22, 1),f(22, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(22, 3),f(22, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(22, 5),f(22, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(22, 7),f(22, 8) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(22, 9),f(22,10) / -1.00000000d-03 , -1.00000000d-03 /
      data f(22,11),f(22,12) / -3.00000000d-03 , -4.00000000d-03 /
      data f(22,13),f(22,14) / -4.00000000d-03 , -3.00000000d-03 /
      data f(22,15) /          -1.00000000d-03 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/-3.03902475d+00,-4.57551287d-02/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-2.08824485d+00,-3.08897426d-02/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 8.64664222d-02,-1.66859008d-02/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 0.00000000d+00,-1.66654019d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-8.64664222d-02, 1.73525169d-02/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.08824485d+00, 2.85565864d-02/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 3.03902475d+00, 5.44211377d-02/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-2.06689518d+01,-4.15811371d-02/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-2.72456716d+02,-6.41365894d-02/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.28362498d+02, 4.93074948d-02/
      data fpp( 1,11,1),fpp( 1,11,2)/ 2.93348972d+02, 3.28366610d-01/
      data fpp( 1,12,1),fpp( 1,12,2)/ 5.06928452d+02,-2.53793936d-01/
      data fpp( 1,13,1),fpp( 1,13,2)/ 5.06928452d+02,-2.27651497d-01/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.93348972d+02, 2.23796856d-01/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.28362498d+02, 4.41444072d-01/
      data fpp( 2, 1,1),fpp( 2, 1,2)/-3.92195049d+00,-4.15960643d-02/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-2.02351030d+00,-2.79078713d-02/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-7.29328444d-02,-1.46524504d-02/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 0.00000000d+00, 1.17672746d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 7.29328444d-02, 1.41817594d-02/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 2.02351030d+00, 2.95552898d-02/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 3.92195049d+00, 3.54770815d-02/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-1.16620964d+01, 3.43963841d-02/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-1.92486569d+02,-1.52302618d-01/
      data fpp( 2,10,1),fpp( 2,10,2)/ 8.11750046d+01, 7.80740879d-02/
      data fpp( 2,11,1),fpp( 2,11,2)/ 2.31002056d+02, 2.18306267d-01/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.91343096d+02,-1.75199154d-01/
      data fpp( 2,13,1),fpp( 2,13,2)/ 3.91343096d+02,-1.59825671d-01/
      data fpp( 2,14,1),fpp( 2,14,2)/ 2.31002056d+02, 1.56812335d-01/
      data fpp( 2,15,1),fpp( 2,15,2)/ 8.11750046d+01, 3.08676333d-01/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-5.27317328d+00,-3.72813595d-02/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-1.81771396d+00,-2.49272810d-02/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.94735045d-01,-1.28295165d-02/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 0.00000000d+00, 1.65346984d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 3.94735045d-01, 1.21681286d-02/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 1.81771396d+00, 2.72421388d-02/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 5.27317328d+00, 2.86833163d-02/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-1.08266245d+00, 5.46445958d-02/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.07197009d+02,-1.28581700d-01/
      data fpp( 3,10,1),fpp( 3,10,2)/ 3.80337484d+02, 2.14422031d-02/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.41557196d+02, 1.63892887d-01/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.57699166d+02,-1.22613752d-01/
      data fpp( 3,13,1),fpp( 3,13,2)/ 3.57699166d+02,-1.09505187d-01/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.41557196d+02, 1.11458625d-01/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.80337484d+02, 2.18070688d-01/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-4.98535639d+00,-3.31338809d-02/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-2.10563385d+00,-2.21022382d-02/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.48126977d-01,-1.11171663d-02/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 0.00000000d+00, 3.09032639d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.48126977d-01, 1.09935532d-02/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 2.10563385d+00, 2.25348839d-02/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 4.98535639d+00, 3.15269112d-02/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 3.99274623d+00, 3.26774714d-02/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-4.05253958d+01,-6.41679681d-03/
      data fpp( 4,10,1),fpp( 4,10,2)/-5.32724941d+02,-1.48610284d-01/
      data fpp( 4,11,1),fpp( 4,11,2)/ 5.71626727d+02, 2.08217933d-01/
      data fpp( 4,12,1),fpp( 4,12,2)/ 3.72602420d+01,-1.05921450d-01/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.72602420d+01,-8.20246176d-02/
      data fpp( 4,14,1),fpp( 4,14,2)/ 5.71626727d+02, 1.12630605d-01/
      data fpp( 4,15,1),fpp( 4,15,2)/-5.32724941d+02, 2.09842197d-01/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-5.98540117d+00,-2.91170915d-02/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.75975062d+00,-1.93958171d-02/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-2.12757048d-01,-9.63964018d-03/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 0.00000000d+00,-5.62218857d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 2.12757048d-01, 9.66212893d-03/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 1.75975062d+00, 1.93171065d-02/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 5.98540117d+00, 2.94094453d-02/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 9.11167753d+00, 2.64251125d-02/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-1.03014080d+01, 2.82101047d-02/
      data fpp( 5,10,1),fpp( 5,10,2)/-1.62237722d+02,-1.17185531d-01/
      data fpp( 5,11,1),fpp( 5,11,2)/-2.85549712d+02, 7.47120206d-02/
      data fpp( 5,12,1),fpp( 5,12,2)/ 4.22059866d+02,-4.96625511d-02/
      data fpp( 5,13,1),fpp( 5,13,2)/ 4.22059866d+02,-3.77983571d-02/
      data fpp( 5,14,1),fpp( 5,14,2)/-2.85549712d+02, 2.72552449d-02/
      data fpp( 5,15,1),fpp( 5,15,2)/-1.62237722d+02, 6.07773776d-02/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-5.50111830d+00,-2.20918586d-02/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.76793122d+00,-1.46762828d-02/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.87665367d-01,-7.16301034d-03/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 0.00000000d+00, 8.32410916d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.87665367d-01, 7.12971390d-03/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.76793122d+00, 1.47928203d-02/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 5.50111830d+00, 2.16590049d-02/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 1.19185943d+01, 2.50511599d-02/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.62169220d+01, 1.98563553d-02/
      data fpp( 6,10,1),fpp( 6,10,2)/ 3.52563604d+00,-1.03965812d-02/
      data fpp( 6,11,1),fpp( 6,11,2)/-2.06314227d+02,-1.56500307d-02/
      data fpp( 6,12,1),fpp( 6,12,2)/-4.74659720d+02,-4.17232961d-02/
      data fpp( 6,13,1),fpp( 6,13,2)/-4.74659720d+02,-4.13950964d-02/
      data fpp( 6,14,1),fpp( 6,14,2)/-2.06314227d+02,-1.69628296d-02/
      data fpp( 6,15,1),fpp( 6,15,2)/ 3.52563604d+00,-5.47358520d-03/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-5.01012565d+00,-1.63524768d-02/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.36852451d+00,-1.08550465d-02/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-8.65814853d-02,-5.26733742d-03/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 0.00000000d+00, 4.39614082d-06/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 8.65814853d-02, 5.24975286d-03/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 1.36852451d+00, 1.09165924d-02/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 5.01012565d+00, 1.61238775d-02/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.13139453d+01, 1.90878978d-02/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.98337202d+01, 1.61245315d-02/
      data fpp( 7,10,1),fpp( 7,10,2)/ 2.46851778d+01, 8.33976348d-04/
      data fpp( 7,11,1),fpp( 7,11,2)/ 6.11066193d+01,-1.97604369d-02/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.28079013d+02,-3.22522289d-02/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.28079013d+02,-3.29530949d-02/
      data fpp( 7,14,1),fpp( 7,14,2)/ 6.11066193d+01,-1.69569729d-02/
      data fpp( 7,15,1),fpp( 7,15,2)/ 2.46851778d+01,-9.67901356d-03/
      data fpp( 8, 1,1),fpp( 8, 1,2)/-4.00837912d+00,-1.19539942d-02/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.15797073d+00,-7.93201161d-03/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-2.16008692d-01,-3.83795938d-03/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 0.00000000d+00, 3.84911662d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 2.16008692d-01, 3.82256291d-03/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 1.15797073d+00, 7.98589924d-03/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 4.00837912d+00, 1.17538401d-02/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 9.57562463d+00, 1.39387402d-02/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 1.78481973d+01, 1.17511989d-02/
      data fpp( 8,10,1),fpp( 8,10,2)/ 2.79336529d+01, 1.87646399d-03/
      data fpp( 8,11,1),fpp( 8,11,2)/ 2.56377497d+01,-1.57770549d-02/
      data fpp( 8,12,1),fpp( 8,12,2)/ 1.34936664d+01,-2.34282444d-02/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.34936664d+01,-2.40967395d-02/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.56377497d+01,-1.31030744d-02/
      data fpp( 8,15,1),fpp( 8,15,2)/ 2.79336529d+01,-8.15096278d-03/
      data fpp( 9, 1,1),fpp( 9, 1,2)/-3.25635788d+00,-8.72210860d-03/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-1.04959256d+00,-5.76578281d-03/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-9.93837456d-02,-2.71476017d-03/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 0.00000000d+00, 4.82347844d-06/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 9.93837456d-02, 2.69546625d-03/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.04959256d+00, 5.83331151d-03/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 3.25635788d+00, 8.47128772d-03/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 7.53355622d+00, 9.90153763d-03/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.43734906d+01, 8.08256178d-03/
      data fpp( 9,10,1),fpp( 9,10,2)/ 2.33302106d+01, 1.26821526d-03/
      data fpp( 9,11,1),fpp( 9,11,2)/ 3.50923819d+01,-1.19554228d-02/
      data fpp( 9,12,1),fpp( 9,12,2)/ 4.62463210d+01,-1.60865240d-02/
      data fpp( 9,13,1),fpp( 9,13,2)/ 4.62463210d+01,-1.66227166d-02/
      data fpp( 9,14,1),fpp( 9,14,2)/ 3.50923819d+01,-9.81065240d-03/
      data fpp( 9,15,1),fpp( 9,15,2)/ 2.33302106d+01,-6.77467380d-03/
      data fpp(10, 1,1),fpp(10, 1,2)/-2.46618936d+00,-6.24359211d-03/
      data fpp(10, 2,1),fpp(10, 2,2)/-7.93659023d-01,-4.12281578d-03/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.36456325d-01,-1.92514478d-03/
      data fpp(10, 4,1),fpp(10, 4,2)/ 0.00000000d+00, 3.39488748d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.36456325d-01, 1.91156523d-03/
      data fpp(10, 6,1),fpp(10, 6,2)/ 7.93659023d-01, 4.17034420d-03/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.46618936d+00, 6.06705796d-03/
      data fpp(10, 8,1),fpp(10, 8,2)/ 5.74015050d+00, 6.90142395d-03/
      data fpp(10, 9,1),fpp(10, 9,2)/ 1.10578404d+01, 5.38724625d-03/
      data fpp(10,10,1),fpp(10,10,2)/ 1.80955048d+01, 5.29591071d-04/
      data fpp(10,11,1),fpp(10,11,2)/ 2.49427229d+01,-9.12561053d-03/
      data fpp(10,12,1),fpp(10,12,2)/ 2.95210497d+01,-1.04671490d-02/
      data fpp(10,13,1),fpp(10,13,2)/ 2.95210497d+01,-1.08957479d-02/
      data fpp(10,14,1),fpp(10,14,2)/ 2.49427229d+01,-7.41121490d-03/
      data fpp(10,15,1),fpp(10,15,2)/ 1.80955048d+01,-5.89939255d-03/
      data fpp(11, 1,1),fpp(11, 1,2)/-1.87888466d+00,-4.36619624d-03/
      data fpp(11, 2,1),fpp(11, 2,2)/-5.75771345d-01,-2.88760751d-03/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.04790953d-01,-1.36337370d-03/
      data fpp(11, 4,1),fpp(11, 4,2)/ 0.00000000d+00, 1.10231421d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.04790953d-01, 1.35896444d-03/
      data fpp(11, 6,1),fpp(11, 6,2)/ 5.75771345d-01, 2.90303991d-03/
      data fpp(11, 7,1),fpp(11, 7,2)/ 1.87888466d+00, 4.30887590d-03/
      data fpp(11, 8,1),fpp(11, 8,2)/ 4.30584178d+00, 4.76145647d-03/
      data fpp(11, 9,1),fpp(11, 9,2)/ 8.14514773d+00, 3.46529821d-03/
      data fpp(11,10,1),fpp(11,10,2)/ 1.34877704d+01, 3.73506768d-05/
      data fpp(11,11,1),fpp(11,11,2)/ 1.93367267d+01,-7.03470092d-03/
      data fpp(11,12,1),fpp(11,12,2)/ 2.40694801d+01,-6.45854700d-03/
      data fpp(11,13,1),fpp(11,13,2)/ 2.40694801d+01,-6.80700855d-03/
      data fpp(11,14,1),fpp(11,14,2)/ 1.93367267d+01,-5.64085470d-03/
      data fpp(11,15,1),fpp(11,15,2)/ 1.34877704d+01,-5.18957265d-03/
      data fpp(12, 1,1),fpp(12, 1,2)/-1.15959155d+00,-2.47846317d-03/
      data fpp(12, 2,1),fpp(12, 2,2)/-2.84989500d-01,-1.64307366d-03/
      data fpp(12, 3,1),fpp(12, 3,2)/ 6.94072802d-03,-7.89242191d-04/
      data fpp(12, 4,1),fpp(12, 4,2)/ 0.00000000d+00, 4.24222296d-08/
      data fpp(12, 5,1),fpp(12, 5,2)/-6.94072802d-03, 7.89072502d-04/
      data fpp(12, 6,1),fpp(12, 6,2)/ 2.84989500d-01, 1.64366757d-03/
      data fpp(12, 7,1),fpp(12, 7,2)/ 1.15959155d+00, 2.47625721d-03/
      data fpp(12, 8,1),fpp(12, 8,2)/ 2.68709372d+00, 2.61130357d-03/
      data fpp(12, 9,1),fpp(12, 9,2)/ 5.07761395d+00, 1.77852850d-03/
      data fpp(12,10,1),fpp(12,10,2)/ 8.17709548d+00,-3.65417558d-04/
      data fpp(12,11,1),fpp(12,11,2)/ 1.15824290d+01,-4.75685826d-03/
      data fpp(12,12,1),fpp(12,12,2)/ 1.42210330d+01,-2.74714938d-03/
      data fpp(12,13,1),fpp(12,13,2)/ 1.42210330d+01,-3.00012272d-03/
      data fpp(12,14,1),fpp(12,14,2)/ 1.15824290d+01,-3.74496494d-03/
      data fpp(12,15,1),fpp(12,15,2)/ 8.17709548d+00,-4.16001753d-03/
      data fpp(13, 1,1),fpp(13, 1,2)/-4.01976234d-01,-9.22070897d-04/
      data fpp(13, 2,1),fpp(13, 2,2)/-1.26570794d-01,-6.25858206d-04/
      data fpp(13, 3,1),fpp(13, 3,2)/-2.33357577d-02,-3.54496281d-04/
      data fpp(13, 4,1),fpp(13, 4,2)/ 0.00000000d+00, 3.84332797d-06/
      data fpp(13, 5,1),fpp(13, 5,2)/ 2.33357577d-02, 3.39122969d-04/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.26570794d-01, 6.79664797d-04/
      data fpp(13, 7,1),fpp(13, 7,2)/ 4.01976234d-01, 7.22217843d-04/
      data fpp(13, 8,1),fpp(13, 8,2)/ 1.16179501d+00, 8.11463832d-04/
      data fpp(13, 9,1),fpp(13, 9,2)/ 2.28854672d+00, 4.71926830d-04/
      data fpp(13,10,1),fpp(13,10,2)/ 3.66063223d+00,-2.99171152d-04/
      data fpp(13,11,1),fpp(13,11,2)/ 5.03019111d+00,-2.21524222d-03/
      data fpp(13,12,1),fpp(13,12,2)/ 6.30700620d+00,-3.19859958d-04/
      data fpp(13,13,1),fpp(13,13,2)/ 6.30700620d+00,-4.52799015d-04/
      data fpp(13,14,1),fpp(13,14,2)/ 5.03019111d+00,-1.68348600d-03/
      data fpp(13,15,1),fpp(13,15,2)/ 3.66063223d+00,-2.29325700d-03/
      data fpp(14, 1,1),fpp(14, 1,2)/-2.08503510d-01,-2.53851257d-04/
      data fpp(14, 2,1),fpp(14, 2,2)/-4.87273260d-02,-1.92297485d-04/
      data fpp(14, 3,1),fpp(14, 3,2)/-9.59769739d-03,-1.76958802d-04/
      data fpp(14, 4,1),fpp(14, 4,2)/ 0.00000000d+00, 1.32694608d-07/
      data fpp(14, 5,1),fpp(14, 5,2)/ 9.59769739d-03, 1.76428024d-04/
      data fpp(14, 6,1),fpp(14, 6,2)/ 4.87273260d-02, 1.94155210d-04/
      data fpp(14, 7,1),fpp(14, 7,2)/ 2.08503510d-01, 2.46951138d-04/
      data fpp(14, 8,1),fpp(14, 8,2)/ 4.41726226d-01, 7.80402390d-05/
      data fpp(14, 9,1),fpp(14, 9,2)/ 9.12199153d-01, 4.08879062d-05/
      data fpp(14,10,1),fpp(14,10,2)/ 1.54037560d+00,-1.81591864d-04/
      data fpp(14,11,1),fpp(14,11,2)/ 2.08880655d+00,-2.14520452d-04/
      data fpp(14,12,1),fpp(14,12,2)/ 1.99894216d+00,-2.20326330d-04/
      data fpp(14,13,1),fpp(14,13,2)/ 1.99894216d+00,-2.21760785d-04/
      data fpp(14,14,1),fpp(14,14,2)/ 2.08880655d+00,-2.08782633d-04/
      data fpp(14,15,1),fpp(14,15,2)/ 1.54037560d+00,-2.03108684d-04/
      data fpp(15, 1,1),fpp(15, 1,2)/-8.40097265d-02,-5.19276612d-05/
      data fpp(15, 2,1),fpp(15, 2,2)/-3.85199025d-02,-4.61446775d-05/
      data fpp(15, 3,1),fpp(15, 3,2)/-1.02734528d-02,-6.34936287d-05/
      data fpp(15, 4,1),fpp(15, 4,2)/ 0.00000000d+00, 1.19192262d-07/
      data fpp(15, 5,1),fpp(15, 5,2)/ 1.02734528d-02, 6.30168596d-05/
      data fpp(15, 6,1),fpp(15, 6,2)/ 3.85199025d-02, 4.78133692d-05/
      data fpp(15, 7,1),fpp(15, 7,2)/ 8.40097265d-02, 4.57296636d-05/
      data fpp(15, 8,1),fpp(15, 8,2)/ 1.91300085d-01, 9.26797642d-06/
      data fpp(15, 9,1),fpp(15, 9,2)/ 3.02656665d-01,-2.28015693d-05/
      data fpp(15,10,1),fpp(15,10,2)/ 4.25865360d-01, 2.19383007d-05/
      data fpp(15,11,1),fpp(15,11,2)/ 4.86582685d-01, 1.15048367d-04/
      data fpp(15,12,1),fpp(15,12,2)/ 5.29225169d-01,-1.82131767d-04/
      data fpp(15,13,1),fpp(15,13,2)/ 5.29225169d-01,-1.71128881d-04/
      data fpp(15,14,1),fpp(15,14,2)/ 4.86582685d-01, 7.10368233d-05/
      data fpp(15,15,1),fpp(15,15,2)/ 4.25865360d-01, 1.86981588d-04/
      data fpp(16, 1,1),fpp(16, 1,2)/-5.54575842d-02, 1.96295118d-05/
      data fpp(16, 2,1),fpp(16, 2,2)/-1.31930639d-02, 1.07409763d-05/
      data fpp(16, 3,1),fpp(16, 3,2)/ 2.69150856d-03,-2.59341711d-06/
      data fpp(16, 4,1),fpp(16, 4,2)/ 0.00000000d+00,-3.67307871d-07/
      data fpp(16, 5,1),fpp(16, 5,2)/-2.69150856d-03, 4.06264859d-06/
      data fpp(16, 6,1),fpp(16, 6,2)/ 1.31930639d-02,-1.58832865d-05/
      data fpp(16, 7,1),fpp(16, 7,2)/ 5.54575842d-02,-5.29502568d-07/
      data fpp(16, 8,1),fpp(16, 8,2)/ 8.90734347d-02, 1.80012968d-05/
      data fpp(16, 9,1),fpp(16, 9,2)/ 1.09174188d-01, 4.85243154d-05/
      data fpp(16,10,1),fpp(16,10,2)/ 1.16162956d-01, 8.79014414d-05/
      data fpp(16,11,1),fpp(16,11,2)/ 1.16862709d-01, 1.98699188d-05/
      data fpp(16,12,1),fpp(16,12,2)/ 1.08157168d-01,-1.07381117d-04/
      data fpp(16,13,1),fpp(16,13,2)/ 1.08157168d-01,-1.07791609d-04/
      data fpp(16,14,1),fpp(16,14,2)/ 1.16862709d-01, 2.15118883d-05/
      data fpp(16,15,1),fpp(16,15,2)/ 1.16162956d-01, 8.17440558d-05/
      data fpp(17, 1,1),fpp(17, 1,2)/-6.15993656d-03, 3.00045607d-05/
      data fpp(17, 2,1),fpp(17, 2,2)/-4.70784168d-03, 1.99908786d-05/
      data fpp(17, 3,1),fpp(17, 3,2)/-4.92581445d-04, 1.00319248d-05/
      data fpp(17, 4,1),fpp(17, 4,2)/ 0.00000000d+00,-1.18577771d-07/
      data fpp(17, 5,1),fpp(17, 5,2)/ 4.92581445d-04,-9.55761370d-06/
      data fpp(17, 6,1),fpp(17, 6,2)/ 4.70784168d-03,-2.16509674d-05/
      data fpp(17, 7,1),fpp(17, 7,2)/ 6.15993656d-03,-2.38385166d-05/
      data fpp(17, 8,1),fpp(17, 8,2)/ 2.84061763d-02,-2.99496627d-06/
      data fpp(17, 9,1),fpp(17, 9,2)/ 5.26465845d-02, 9.58183817d-05/
      data fpp(17,10,1),fpp(17,10,2)/ 4.54828160d-02, 9.97214396d-05/
      data fpp(17,11,1),fpp(17,11,2)/ 2.99664811d-02,-7.47041401d-05/
      data fpp(17,12,1),fpp(17,12,2)/ 4.61461604d-02,-4.09048792d-05/
      data fpp(17,13,1),fpp(17,13,2)/ 4.61461604d-02,-4.99332923d-05/
      data fpp(17,14,1),fpp(17,14,2)/ 2.99664811d-02,-3.85904879d-05/
      data fpp(17,15,1),fpp(17,15,2)/ 4.54828160d-02,-3.57047560d-05/
      data fpp(18, 1,1),fpp(18, 1,2)/-1.79139822d-03, 3.92786819d-05/
      data fpp(18, 2,1),fpp(18, 2,2)/ 2.72005703d-03, 2.14426363d-05/
      data fpp(18, 3,1),fpp(18, 3,2)/ 1.31990056d-04,-5.04922699d-06/
      data fpp(18, 4,1),fpp(18, 4,2)/ 0.00000000d+00,-1.24572833d-06/
      data fpp(18, 5,1),fpp(18, 5,2)/-1.31990056d-04, 1.00321403d-05/
      data fpp(18, 6,1),fpp(18, 6,2)/-2.72005703d-03,-3.88828329d-05/
      data fpp(18, 7,1),fpp(18, 7,2)/ 1.79139822d-03, 2.54991912d-05/
      data fpp(18, 8,1),fpp(18, 8,2)/-9.75524612d-03,-3.11393180d-06/
      data fpp(18, 9,1),fpp(18, 9,2)/-1.45268472d-02,-1.30434640d-05/
      data fpp(18,10,1),fpp(18,10,2)/-2.52992582d-03, 5.52877877d-05/
      data fpp(18,11,1),fpp(18,11,2)/ 1.36692026d-02,-2.81076869d-05/
      data fpp(18,12,1),fpp(18,12,2)/ 1.74829349d-02,-2.85704021d-06/
      data fpp(18,13,1),fpp(18,13,2)/ 1.74829349d-02,-7.37503593d-06/
      data fpp(18,14,1),fpp(18,14,2)/ 1.36692026d-02,-1.00357040d-05/
      data fpp(18,15,1),fpp(18,15,2)/-2.52992582d-03,-1.24821480d-05/
      data fpp(19, 1,1),fpp(19, 1,2)/ 7.32552942d-03,-1.96551697d-05/
      data fpp(19, 2,1),fpp(19, 2,2)/-1.72386419d-04,-1.06896606d-05/
      data fpp(19, 3,1),fpp(19, 3,2)/-3.53787778d-05, 2.41381202d-06/
      data fpp(19, 4,1),fpp(19, 4,2)/ 0.00000000d+00, 1.03441249d-06/
      data fpp(19, 5,1),fpp(19, 5,2)/ 3.53787778d-05,-6.55146199d-06/
      data fpp(19, 6,1),fpp(19, 6,2)/ 1.72386419d-04, 2.51714355d-05/
      data fpp(19, 7,1),fpp(19, 7,2)/-7.32552942d-03,-3.41342799d-05/
      data fpp(19, 8,1),fpp(19, 8,2)/-1.38519176d-03, 5.13656842d-05/
      data fpp(19, 9,1),fpp(19, 9,2)/-6.53919559d-03,-5.13284569d-05/
      data fpp(19,10,1),fpp(19,10,2)/-1.13631127d-02, 3.39481434d-05/
      data fpp(19,11,1),fpp(19,11,2)/-1.26432914d-02,-2.44641168d-05/
      data fpp(19,12,1),fpp(19,12,2)/-1.40779001d-02, 3.90832393d-06/
      data fpp(19,13,1),fpp(19,13,2)/-1.40779001d-02, 5.07086624d-07/
      data fpp(19,14,1),fpp(19,14,2)/-1.26432914d-02,-1.08591676d-05/
      data fpp(19,15,1),fpp(19,15,2)/-1.13631127d-02,-1.70704162d-05/
      data fpp(20, 1,1),fpp(20, 1,2)/-3.51071946d-03, 2.77738972d-08/
      data fpp(20, 2,1),fpp(20, 2,2)/-2.03051135d-03,-5.55477944d-08/
      data fpp(20, 3,1),fpp(20, 3,2)/ 9.52505557d-06, 1.94417280d-07/
      data fpp(20, 4,1),fpp(20, 4,2)/ 0.00000000d+00,-7.22121327d-07/
      data fpp(20, 5,1),fpp(20, 5,2)/-9.52505557d-06, 2.69406803d-06/
      data fpp(20, 6,1),fpp(20, 6,2)/ 2.03051135d-03,-1.00541508d-05/
      data fpp(20, 7,1),fpp(20, 7,2)/ 3.51071946d-03, 3.75225351d-05/
      data fpp(20, 8,1),fpp(20, 8,2)/-2.70398683d-03,-2.00359897d-05/
      data fpp(20, 9,1),fpp(20, 9,2)/-1.31637042d-03,-1.73785762d-05/
      data fpp(20,10,1),fpp(20,10,2)/-1.76235130d-05, 2.95502946d-05/
      data fpp(20,11,1),fpp(20,11,2)/-5.09603693d-03,-4.08226020d-05/
      data fpp(20,12,1),fpp(20,12,2)/-3.17133459d-03, 1.37401136d-05/
      data fpp(20,13,1),fpp(20,13,2)/-3.17133459d-03, 9.19096024d-06/
      data fpp(20,14,1),fpp(20,14,2)/-5.09603693d-03,-2.26259886d-05/
      data fpp(20,15,1),fpp(20,15,2)/-1.76235130d-05,-3.86870057d-05/
      data fpp(21, 1,1),fpp(21, 1,2)/ 7.17348417d-04, 4.85347394d-05/
      data fpp(21, 2,1),fpp(21, 2,2)/ 2.29443181d-03, 2.29305211d-05/
      data fpp(21, 3,1),fpp(21, 3,2)/-2.72144445d-06,-2.02568239d-05/
      data fpp(21, 4,1),fpp(21, 4,2)/ 0.00000000d+00,-1.90322560d-06/
      data fpp(21, 5,1),fpp(21, 5,2)/ 2.72144445d-06, 2.78697263d-05/
      data fpp(21, 6,1),fpp(21, 6,2)/-2.29443181d-03,-4.95756795d-05/
      data fpp(21, 7,1),fpp(21, 7,2)/-7.17348417d-04, 5.04329916d-05/
      data fpp(21, 8,1),fpp(21, 8,2)/ 2.01139095d-04,-3.21562870d-05/
      data fpp(21, 9,1),fpp(21, 9,2)/-1.95322738d-04, 1.81921564d-05/
      data fpp(21,10,1),fpp(21,10,2)/-5.66393282d-04,-4.06123385d-05/
      data fpp(21,11,1),fpp(21,11,2)/ 3.02743912d-03, 2.42571977d-05/
      data fpp(21,12,1),fpp(21,12,2)/ 2.76323846d-03, 3.58354760d-06/
      data fpp(21,13,1),fpp(21,13,2)/ 2.76323846d-03, 7.12075834d-06/
      data fpp(21,14,1),fpp(21,14,2)/ 3.02743912d-03, 1.01083548d-05/
      data fpp(21,15,1),fpp(21,15,2)/-5.66393282d-04, 1.24458226d-05/
      data fpp(22, 1,1),fpp(22, 1,2)/ 6.41325792d-04, 4.90389031d-09/
      data fpp(22, 2,1),fpp(22, 2,2)/ 4.85278409d-03,-9.80778062d-09/
      data fpp(22, 3,1),fpp(22, 3,2)/ 1.36072222d-06, 3.43272322d-08/
      data fpp(22, 4,1),fpp(22, 4,2)/ 0.00000000d+00,-1.27501148d-07/
      data fpp(22, 5,1),fpp(22, 5,2)/-1.36072222d-06, 4.75677360d-07/
      data fpp(22, 6,1),fpp(22, 6,2)/-4.85278409d-03,-1.77520829d-06/
      data fpp(22, 7,1),fpp(22, 7,2)/-6.41325792d-04, 6.62515581d-06/
      data fpp(22, 8,1),fpp(22, 8,2)/ 1.89943045d-03,-2.47254150d-05/
      data fpp(22, 9,1),fpp(22, 9,2)/ 2.09766137d-03, 3.22765040d-05/
      data fpp(22,10,1),fpp(22,10,2)/ 2.28319664d-03,-4.43806011d-05/
      data fpp(22,11,1),fpp(22,11,2)/ 1.09862804d-02, 2.52459003d-05/
      data fpp(22,12,1),fpp(22,12,2)/ 1.01183808d-02, 3.39699995d-06/
      data fpp(22,13,1),fpp(22,13,2)/ 1.01183808d-02, 7.18605002d-06/
      data fpp(22,14,1),fpp(22,14,2)/ 1.09862804d-02, 1.00897000d-05/
      data fpp(22,15,1),fpp(22,15,2)/ 2.28319664d-03, 1.24551500d-05/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 4.60505349d-02, 1.43190528d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/-2.98908636d-02,-3.04262601d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.51595904d-01,-4.40286694d-01/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.53611499d+00,-7.18369865d+00/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.35138924d+01,-7.71545228d+00/
      data fpppp( 1,11),fpppp( 1,12)/ 3.19795241d+00,-2.16077700d+00/
      data fpppp( 1,13),fpppp( 1,14)/-1.52402961d+00, 6.50962836d-01/
      data fpppp( 1,15) /             1.83575861d+00 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 1.91644177d-02, 1.83091968d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/-2.33598610d-02,-2.10501523d-02/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.07560470d-01,-2.96533052d-01/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.07544350d+00,-5.05419018d+00/
      data fpppp( 2, 9),fpppp( 2,10)/ 9.22689169d+00,-4.58421384d+00/
      data fpppp( 2,11),fpppp( 2,12)/ 1.67989234d+00,-1.50451622d+00/
      data fpppp( 2,13),fpppp( 2,14)/-1.13662871d+00, 2.08342306d-01/
      data fpppp( 2,15) /             9.34098781d-01 /
      data fpppp( 3, 1),fpppp( 3, 2)/-2.93239164d-02,-2.24113097d-02/
      data fpppp( 3, 3),fpppp( 3, 4)/-2.97966880d-03,-2.73646475d-02/
      data fpppp( 3, 5),fpppp( 3, 6)/ 1.12438259d-01,-3.60693755d-01/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.45228558d+00,-6.03712629d+00/
      data fpppp( 3, 9),fpppp( 3,10)/ 1.67107089d+01,-2.51867790d+01/
      data fpppp( 3,11),fpppp( 3,12)/ 2.34706569d+01,-7.42678619d+00/
      data fpppp( 3,13),fpppp( 3,14)/-4.43266074d+00, 1.14941551d+01/
      data fpppp( 3,15) /             1.97251027d+01 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.20449136d-03,-5.51420198d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/-3.10716400d-02, 2.12379678d-02/
      data fpppp( 4, 5),fpppp( 4, 6)/-5.38802312d-02, 3.02845751d-01/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.10216983d+00, 3.87349363d+00/
      data fpppp( 4, 9),fpppp( 4,10)/-1.70033366d+01, 3.72789685d+01/
      data fpppp( 4,11),fpppp( 4,12)/-3.63194648d+01, 9.67580149d+00/
      data fpppp( 4,13),fpppp( 4,14)/ 5.16332249d+00,-1.82695488d+01/
      data fpppp( 4,15) /            -3.04082165d+01 /
      data fpppp( 5, 1),fpppp( 5, 2)/-4.02716166d-02,-2.67048854d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.36282606d-02, 1.16373666d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/ 8.97331401d-03, 4.29971986d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/-2.02426895d-02,-2.79888922d-02/
      data fpppp( 5, 9),fpppp( 5,10)/-1.22016346d+00,-3.04275099d+00/
      data fpppp( 5,11),fpppp( 5,12)/ 1.51086268d+01,-7.53646213d+00/
      data fpppp( 5,13),fpppp( 5,14)/-6.17321438d+00, 9.65563580d+00/
      data fpppp( 5,15) /             1.74059653d+01 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.93348114d-02,-2.11244289d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.53427467d-02,-1.06061347d-03/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.95852005d-02, 6.27584039d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/ 8.44867116d-02,-1.83165351d-01/
      data fpppp( 6, 9),fpppp( 6,10)/ 5.21025793d-01,-2.92031463d+00/
      data fpppp( 6,11),fpppp( 6,12)/-6.68681870d-01, 2.08470429d+00/
      data fpppp( 6,13),fpppp( 6,14)/ 2.13059286d+00,-8.52236152d-01/
      data fpppp( 6,15) /            -2.23198607d+00 /
      data fpppp( 7, 1),fpppp( 7, 2)/-3.52667862d-02,-2.35421021d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/-1.21442919d-02, 3.97577295d-04/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.05539828d-02, 2.91081842d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.45927669d-02, 7.22538576d-02/
      data fpppp( 7, 9),fpppp( 7,10)/-1.70650880d-01, 3.90250620d-01/
      data fpppp( 7,11),fpppp( 7,12)/ 5.03847433d-01,-5.72583194d-01/
      data fpppp( 7,13),fpppp( 7,14)/-5.43345958d-01, 3.86898489d-01/
      data fpppp( 7,15) /             8.28809160d-01 /
      data fpppp( 8, 1),fpppp( 8, 2)/-3.07329339d-02,-1.94373828d-02/
      data fpppp( 8, 3),fpppp( 8, 4)/-6.02431599d-03,-2.25540908d-05/
      data fpppp( 8, 5),fpppp( 8, 6)/ 6.11453235d-03, 1.91216255d-02/
      data fpppp( 8, 7),fpppp( 8, 8)/ 3.19057466d-02, 1.62656156d-02/
      data fpppp( 8, 9),fpppp( 8,10)/ 6.53514216d-02,-1.68898327d-01/
      data fpppp( 8,11),fpppp( 8,12)/-1.32639642d-01, 1.08566086d-01/
      data fpppp( 8,13),fpppp( 8,14)/ 1.04944061d-01,-1.18151543d-01/
      data fpppp( 8,15) /            -2.23228698d-01 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.67910249d-02,-1.22292740d-02/
      data fpppp( 9, 3),fpppp( 9, 4)/-9.68526930d-03,-7.91530290d-05/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.00018814d-02, 1.11211316d-02/
      data fpppp( 9, 7),fpppp( 9, 8)/ 2.09069824d-02, 2.94769202d-02/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.49494973d-02, 3.77322310d-02/
      data fpppp( 9,11),fpppp( 9,12)/ 2.44865415d-03,-8.40207773d-02/
      data fpppp( 9,13),fpppp( 9,14)/-8.37801685d-02, 1.48621911d-03/
      data fpppp( 9,15) /             4.13413623d-02 /
      data fpppp(10, 1),fpppp(10, 2)/-1.51083343d-02,-1.01347861d-02/
      data fpppp(10, 3),fpppp(10, 4)/-5.27217989d-03,-2.12767101d-05/
      data fpppp(10, 5),fpppp(10, 6)/ 5.35728673d-03, 9.83691216d-03/
      data fpppp(10, 7),fpppp(10, 8)/ 1.62147232d-02, 2.13900425d-02/
      data fpppp(10, 9),fpppp(10,10)/ 2.08488342d-02,-1.58691543d-03/
      data fpppp(10,11),fpppp(10,12)/-2.59279459d-02,-3.08347750d-02/
      data fpppp(10,13),fpppp(10,14)/-3.18815081d-02,-2.17410136d-02/
      data fpppp(10,15) /            -1.72879117d-02 /
      data fpppp(11, 1),fpppp(11, 2)/-1.29437303d-02,-8.39539677d-03/
      data fpppp(11, 3),fpppp(11, 4)/-3.40265802d-03, 3.46625235d-05/
      data fpppp(11, 5),fpppp(11, 6)/ 3.26400792d-03, 8.88067210d-03/
      data fpppp(11, 7),fpppp(11, 8)/ 1.11412791d-02, 1.39848399d-02/
      data fpppp(11, 9),fpppp(11,10)/ 1.76602910d-02, 5.57299992d-03/
      data fpppp(11,11),fpppp(11,12)/-9.57227224d-03,-3.42560847d-02/
      data fpppp(11,13),fpppp(11,14)/-3.44282125d-02,-8.88376091d-03/
      data fpppp(11,15) /             2.99108242d-03 /
      data fpppp(12, 1),fpppp(12, 2)/-8.67272470d-03,-5.81072276d-03/
      data fpppp(12, 3),fpppp(12, 4)/-3.04469380d-03, 5.72405897d-05/
      data fpppp(12, 5),fpppp(12, 6)/ 2.81573144d-03, 6.61209101d-03/
      data fpppp(12, 7),fpppp(12, 8)/ 5.69621403d-03, 9.77705996d-03/
      data fpppp(12, 9),fpppp(12,10)/ 6.97662931d-03, 4.85410097d-03/
      data fpppp(12,11),fpppp(12,12)/-8.04191236d-03,-1.86902230d-02/
      data fpppp(12,13),fpppp(12,14)/-1.90664953d-02,-6.53682343d-03/
      data fpppp(12,15) /            -7.89982514d-04 /
      data fpppp(13, 1),fpppp(13, 2)/-2.63695114d-03,-1.73663239d-03/
      data fpppp(13, 3),fpppp(13, 4)/-7.46743580d-04,-7.03499841d-05/
      data fpppp(13, 5),fpppp(13, 6)/ 1.02814352d-03, 7.51732615d-04/
      data fpppp(13, 7),fpppp(13, 8)/ 6.29515032d-03, 3.13246635d-03/
      data fpppp(13, 9),fpppp(13,10)/ 3.19096030d-03,-1.17627992d-03/
      data fpppp(13,11),fpppp(13,12)/ 1.36256169d-03,-9.83859379d-03/
      data fpppp(13,13),fpppp(13,14)/-9.46995230d-03,-1.12004293d-04/
      data fpppp(13,15) /             4.35334251d-03 /
      data fpppp(14, 1),fpppp(14, 2)/-2.09561393d-03,-1.25046128d-03/
      data fpppp(14, 3),fpppp(14, 4)/-1.41334249d-04, 4.38824084d-05/
      data fpppp(14, 5),fpppp(14, 6)/-3.41953845d-05, 1.86481500d-03/
      data fpppp(14, 7),fpppp(14, 8)/-1.86271310d-04, 3.28706216d-03/
      data fpppp(14, 9),fpppp(14,10)/ 1.27303533d-03, 1.08300789d-03/
      data fpppp(14,11),fpppp(14,12)/-1.03897969d-02, 2.17845905d-03/
      data fpppp(14,13),fpppp(14,14)/ 1.35545315d-03,-7.09777331d-03/
      data fpppp(14,15) /            -1.12620806d-02 /
      data fpppp(15, 1),fpppp(15, 2)/-1.71698196d-04,-1.59312382d-04/
      data fpppp(15, 3),fpppp(15, 4)/-2.25654732d-04,-1.64485051d-05/
      data fpppp(15, 5),fpppp(15, 6)/ 2.91448753d-04,-7.09666891d-05/
      data fpppp(15, 7),fpppp(15, 8)/ 1.02702046d-03,-3.29083080d-04/
      data fpppp(15, 9),fpppp(15,10)/ 5.33285156d-04,-1.09293061d-03/
      data fpppp(15,11),fpppp(15,12)/ 8.89550349d-05,-3.47379962d-04/
      data fpppp(15,13),fpppp(15,14)/-2.81612147d-04,-1.74116227d-04/
      data fpppp(15,15) /            -1.06413381d-04 /
      data fpppp(16, 1),fpppp(16, 2)/-3.46387445d-04,-2.54700877d-04/
      data fpppp(16, 3),fpppp(16, 4)/-2.17605913d-04, 1.05596649d-05/
      data fpppp(16, 5),fpppp(16, 6)/ 1.75367253d-04, 4.02536186d-04/
      data fpppp(16, 7),fpppp(16, 8)/-2.02715130d-04,-1.10595854d-04/
      data fpppp(16, 9),fpppp(16,10)/-1.65807307d-04,-1.28940004d-05/
      data fpppp(16,11),fpppp(16,12)/-1.59957629d-04, 8.84069134d-05/
      data fpppp(16,13),fpppp(16,14)/ 7.59242987d-05,-1.10027171d-04/
      data fpppp(16,15) /            -2.00133220d-04 /
      data fpppp(17, 1),fpppp(17, 2)/ 8.88769953d-05, 3.48578537d-05/
      data fpppp(17, 3),fpppp(17, 4)/-6.25184881d-05,-8.14462903d-06/
      data fpppp(17, 5),fpppp(17, 6)/ 9.50970042d-05,-1.48882660d-04/
      data fpppp(17, 7),fpppp(17, 8)/ 3.34643715d-04, 5.79564919d-05/
      data fpppp(17, 9),fpppp(17,10)/-4.46819572d-04,-1.54928807d-04/
      data fpppp(17,11),fpppp(17,12)/ 5.65380817d-04,-2.04833603d-04/
      data fpppp(17,13),fpppp(17,14)/-1.53579979d-04, 3.60366319d-04/
      data fpppp(17,15) /             6.13875558d-04 /
      data fpppp(18, 1),fpppp(18, 2)/-1.62293122d-04,-7.95114044d-05/
      data fpppp(18, 3),fpppp(18, 4)/ 5.43674073d-05, 9.40639012d-06/
      data fpppp(18, 5),fpppp(18, 6)/-9.19929678d-05, 2.11200866d-04/
      data fpppp(18, 7),fpppp(18, 8)/-3.26839164d-04, 1.32669814d-04/
      data fpppp(18, 9),fpppp(18,10)/ 2.02662501d-04, 6.27915323d-05/
      data fpppp(18,11),fpppp(18,12)/-2.01696209d-04, 8.69537060d-07/
      data fpppp(18,13),fpppp(18,14)/-1.61724767d-05,-1.33528154d-04/
      data fpppp(18,15) /            -1.92838673d-04 /
      data fpppp(19, 1),fpppp(19, 2)/ 1.51047657d-04, 8.16834370d-05/
      data fpppp(19, 3),fpppp(19, 4)/-1.96859964d-05,-9.03718318d-06/
      data fpppp(19, 5),fpppp(19, 6)/ 5.58347291d-05,-2.08204001d-04/
      data fpppp(19, 7),fpppp(19, 8)/ 3.18885868d-04,-2.61044261d-04/
      data fpppp(19, 9),fpppp(19,10)/ 5.96306856d-05, 4.23267238d-05/
      data fpppp(19,11),fpppp(19,12)/-1.63132811d-05, 1.36606059d-05/
      data fpppp(19,13),fpppp(19,14)/ 1.02130828d-05,-2.52318895d-06/
      data fpppp(19,15) /            -9.38612193d-06 /
      data fpppp(20, 1),fpppp(20, 2)/ 2.98032648d-05, 9.37611424d-06/
      data fpppp(20, 3),fpppp(20, 4)/-3.37180241d-05, 2.52229474d-06/
      data fpppp(20, 5),fpppp(20, 6)/ 2.36288452d-05, 2.59360121d-05/
      data fpppp(20, 7),fpppp(20, 8)/-1.60962591d-04, 1.56219489d-04/
      data fpppp(20, 9),fpppp(20,10)/-7.77600111d-06,-1.30447655d-04/
      data fpppp(20,11),fpppp(20,12)/ 1.46937002d-04,-3.71134075d-05/
      data fpppp(20,13),fpppp(20,14)/-1.98693483d-05, 7.79607649d-05/
      data fpppp(20,15) /             1.28213234d-04 /
      data fpppp(21, 1),fpppp(21, 2)/-9.73165848d-05,-4.50761573d-05/
      data fpppp(21, 3),fpppp(21, 4)/ 4.51670145d-05, 2.40058152d-06/
      data fpppp(21, 5),fpppp(21, 6)/-5.47693406d-05, 7.86842986d-05/
      data fpppp(21, 7),fpppp(21, 8)/-2.75136544d-05,-8.14543413d-06/
      data fpppp(21, 9),fpppp(21,10)/-1.88015698d-05, 8.48751908d-05/
      data fpppp(21,11),fpppp(21,12)/-8.28050163d-05, 1.48628903d-05/
      data fpppp(21,13),fpppp(21,14)/ 4.73985723d-06,-4.23128840d-05/
      data fpppp(21,15) /            -6.69703054d-05 /
      data fpppp(22, 1),fpppp(22, 2)/-2.22757732d-04,-1.04629872d-04/
      data fpppp(22, 3),fpppp(22, 4)/ 9.75043198d-05, 5.61635149d-06/
      data fpppp(22, 5),fpppp(22, 6)/-1.19969726d-04, 1.83258793d-04/
      data fpppp(22, 7),fpppp(22, 8)/-6.92925450d-05,-6.33073632d-06/
      data fpppp(22, 9),fpppp(22,10)/-4.59360294d-05, 1.89313115d-04/
      data fpppp(22,11),fpppp(22,12)/-2.00263520d-04, 3.74819564d-05/
      data fpppp(22,13),fpppp(22,14)/ 1.37228809d-05,-1.05227218d-04/
      data fpppp(22,15) /            -1.67073017d-04 /
 
      data x( 1), x( 2) /  3.60000000d+00 ,  3.70000000d+00 /
      data x( 3), x( 4) /  3.80000000d+00 ,  3.90000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.40000000d+00 ,  4.60000000d+00 /
      data x( 9), x(10) /  4.80000000d+00 ,  5.00000000d+00 /
      data x(11), x(12) /  5.20000000d+00 ,  5.50000000d+00 /
      data x(13), x(14) /  6.00000000d+00 ,  6.50000000d+00 /
      data x(15), x(16) /  7.00000000d+00 ,  7.50000000d+00 /
      data x(17), x(18) /  8.00000000d+00 ,  9.00000000d+00 /
      data x(19), x(20) /  1.00000000d+01 ,  1.10000000d+01 /
      data x(21), x(22) /  1.20000000d+01 ,  1.30000000d+01 /
 
      data y( 1), y( 2) / -3.00000000d+01 , -2.00000000d+01 /
      data y( 3), y( 4) / -1.00000000d+01 ,  0.00000000d+00 /
      data y( 5), y( 6) /  1.00000000d+01 ,  2.00000000d+01 /
      data y( 7), y( 8) /  3.00000000d+01 ,  4.00000000d+01 /
      data y( 9), y(10) /  5.00000000d+01 ,  6.00000000d+01 /
      data y(11), y(12) /  7.00000000d+01 ,  8.00000000d+01 /
      data y(13), y(14) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(15) /         1.20000000d+02 /
 
      data delx( 1), delx( 2) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 3), delx( 4) /  1.00000000d-01 ,  1.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 9), delx(10) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx(11), delx(12) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx(13), delx(14) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(15), delx(16) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx(17), delx(18) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(19), delx(20) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(21) /            1.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  22 , 15 /

      iprint=0

      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        do 10 i=1,nptx
          if(xi .gt. x(i))go to 10
          ix=i-1
          go to 15
 10     continue
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
   30 c(1,1)=fpppp(ix+1,iy+1)
      c(1,2)=fpppp(ix+1,iy  )
      c(2,1)=fpppp(ix  ,iy+1)
      c(2,2)=fpppp(ix,  iy  )
      c(1,3)=fpp(ix+1,iy+1,1)
      c(1,4)=fpp(ix+1,iy  ,1)
      c(2,3)=fpp(ix,  iy+1,1)
      c(2,4)=fpp(ix,  iy,  1)
      c(3,1)=fpp(ix+1,iy+1,2)
      c(3,2)=fpp(ix+1,iy  ,2)
      c(4,1)=fpp(ix,  iy+1,2)
      c(4,2)=fpp(ix,  iy,  2)
      c(3,3)=f(ix+1,iy+1)
      c(3,4)=f(ix+1,iy  )
      c(4,3)=f(ix,  iy+1)
      c(4,4)=f(ix,  iy  )
      px(1)=((xi-xix)**3/(6.0*delxi))-(xi-xix)*delxi/6.0
      px(2)=(xi-xixp1)*delxi/6.0-((xi-xixp1)**3/(6.0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0*delyi))-(yi-yiy)*delyi/6.0
      py(2)=(yi-yiyp1)*delyi/6.0-((yi-yiyp1)**3/(6.0*delyi))
      py(3)=(yi-yiy)/delyi
      py(4)=(yiyp1-yi)/delyi

      fi=0.0d 00
      do l=1,4
        sum=0.0d 00
        do i=1,4
          sum=sum + c(i,l) * px(i)
        enddo
        fi = fi + sum * py(l)
      enddo
      return
      end  
