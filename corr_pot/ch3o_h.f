c     ch3o_h and ch3o_h_rlx (with relaxation) potentials 

      real*8 function ch3o_h_rlx(r, rpar, ipar)

c     
c     this is a general subroutine for evaluating the potential 
c     energy for the transitional modes.
c
c the potential should be output in au in vtot
c the cartesian coordinates of the atoms of each fragment are 
c stored in r
c the bond length or some other suitable approximation to a 
c distinguished reaction coordinate should be passed back in rmepi
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
c the valence stretching and bending potential works in terms  of
c either 
c (i) an exponential extrapolation of an equilibrium set 
c of force constants 
c or
c (ii) an exponential interpolation of force constants along the 
c bond stretching path.
c
c for the ion-dipole component some work needs to be done on 
c defining the axis for the dipole.  once this is done it should 
c also be easy to include the ion-quadrupole interaction.
c     

      implicit double precision (a-h,o-z)

      include 'param2.fi'
      include 'commonpot.fi'

c
c vectors for primary coordinates 
c

      dimension r(nfrag,natommx,ndim)
      dimension cm(nfrag,ndim),rij(natommx,natommx)
      dimension qq(5),rmov(ndim),rmhv(ndim),regrlx(9),egrlxr(9)

c
c vectors from mcpot calling routine
c

      dimension ift(nfrag),natom(nfrag),ineutf(nchmx)
      dimension dmss(nfrag,natommx),dmsst(nfrag)

c
c intermediate use vectors
c

      dimension r1(ndim),r2p(ndim),r2d(ndim),r2q(ndim)
      dimension nr(nfrag)
      dimension rbvec(ndim),rax1v(ndim),rbx1v(ndim),rax2v(ndim),
     $          rbx2v(ndim),rx1x2v(ndim),ray1v(ndim),rby1v(ndim),
     $          ray2v(ndim),rby2v(ndim),ry1y2v(ndim)
      dimension cptmp1(ndim),cptmp2(ndim),cptmp3(ndim)

c
c specific potential parameter vectors
c

      dimension fijrb(5,5),qqerb(5)
      dimension rpathv(nrpotm),vrr(nrpotm)

      include 'data.fi'
      data dmss(1,1), dmss(1,2) /15.99491 , 12./
c
c for generality/ease of modification purposes i have set this up to 
c begin by evaluating the interfragment atom-atom distances
c and a particular set of valence stretching, bending and torsional 
c angles.
c for efficient use one might wish to only evaluate the components 
c that are required for the particular potential being used.  for a 
c variety of reasons this seems to me to be best left for the specific 
c applications.
c      
      pi2 = asin(1.0d0)
      pi = 2.0d0*pi2

c
c now evaluate atom-atom separations
c

      rij(1,1) = dsqrt((r(1,1,1)-r(2,1,1))**2+(r(1,1,2)-r(2,1,2))**2+ 
     $                 (r(1,1,3)-r(2,1,3))**2)
      rbond = rij(1,1)
      rmepi = rbond

c
c now proceed to the valence bending/torsional coordinates.
c

c
c determine distances needed for bending angles
c
      rbond2 = 0.0d0
      rax12 = 0.0d0
      rbx12 = 0.0d0
      rby12 = 0.0d0
      ray12 = 0.0d0
      rmo2 = 0.0d0
      rmh2 = 0.0d0
      dm12t = dmss(1,1)+dmss(1,2)
      do 1500 idim = 1 , ndim
         rbvec(idim) = r(2,1,idim) - r(1,1,idim)
         rbond2 = rbond2 + rbvec(idim)**2
         rax1v(idim) = r(1,2,idim) - r(1,1,idim)
         rax12 = rax12 + rax1v(idim)**2
         rbx1v(idim) = r(2,1,idim) - r(1,2,idim)
         rbx12 = rbx12 + rbx1v(idim)**2
         rmov(idim)=(dmss(1,2)*r(1,2,idim)+dmss(1,1)*r(1,1,idim))/
     $    dm12t - r(1,1,idim)
         rmo2 = rmo2 + rmov(idim)**2
         rmhv(idim)=(dmss(1,2)*r(1,2,idim)+dmss(1,1)*r(1,1,idim))/
     $    dm12t - r(2,1,idim)
         rmh2 = rmh2 + rmhv(idim)**2
 1500 continue
      rbond = sqrt(rbond2)
      rax1 = sqrt(rax12)
      rbx1 = sqrt(rbx12)
      ray1 = sqrt(ray12)
      rby1 = sqrt(rby12)
      rmo = sqrt(rmo2)
      rmh = sqrt(rmh2)
c     write (6,*) 'rmo test',rmov,rmo
c     write (6,*) 'rmh test',rmhv,rmh
c     write (6,*) 'rbvec test',rbvec,rbond
c     write (6,*) 'r11 test',(r(1,1,idim),idim=1,ndim)
c     write (6,*) 'r12 test',(r(1,2,idim),idim=1,ndim)
c     write (6,*) 'r21 test',(r(2,1,idim),idim=1,ndim)
c     write (6,*) 'mss test',dmss(1,1),dmss(1,2),dm12t

c     
c     calculate the bending angles.
c     

      cthep = (rmo2+rmh2-rbond2)/(2.0d0*rmo*rmh)
      if (dabs(cthep).le.1.0d0) go to 1630
      if(dabs(cthep)-1.d0.gt.1.d-10) 
     x   write (6,*) 'error in thep, cthep = ',cthep
      thep = 0.0d0
      if (cthep.lt.1.0d0) thep = 2.0d0*dasin(1.0d0)
      go to 1650
 1630 continue
      thep = dacos(cthep)
 1650 continue

      cthe1 = (rax12 + rbond2 - rbx12)/(2.0d0*rax1*rbond)
      if (dabs(cthe1).le.1.0d0) go to 2030
      if(dabs(cthe1)-1.d0.gt.1.e-6) 
     x   write (6,*) 'error in the1, cthe1 = ',cthe1
      the1 = 0.0d0
      if (cthe1.lt.1.0d0) the1 = 2.0d0*dasin(1.0d0)
      go to 2050
 2030 continue
c      write (6,*) 'cthe1 test',cthe1
      the1 = dacos(cthe1)
 2050 continue
      
c     
c     calculate the torsional angles if ivb.ne.0
c     
c
c
c determine distances needed for torsional angles
c
      rax22 = 0.0d0
      rbx22 = 0.0d0
      rby22 = 0.0d0
      ray22 = 0.0d0
      ry1y22 = 0.0d0
      rx1x22 = 0.0d0
      do idim = 1 , ndim
         rax2v(idim) = r(1,3,idim) - r(1,1,idim)
         rax22 = rax22 + rax2v(idim)**2
         rbx2v(idim) = r(2,1,idim) - r(1,3,idim)
         rbx22 = rbx22 + rbx2v(idim)**2
         rx1x2v(idim) = r(1,3,idim) - r(1,2,idim)
         rx1x22 = rx1x22 + rx1x2v(idim)**2
      enddo
      rax2 = sqrt(rax22)
      rbx2 = sqrt(rbx22)
      ray2 = sqrt(ray22)
      rby2 = sqrt(rby22)
      ry1y2 = sqrt(ry1y22)
      rx1x2 = sqrt(rx1x22)

c
c     now evaluate tau 3 (torsion angle for x2,11,x1,21.
c     tau3b corresponds to torsion angle for perp to x1,21;
c     x2,21,x1,perp
c

c
c     consider the terminal atom case.
c
c     first evaluate internal x2,11,x1 bending angle
c

      ctaint = (rax12 + rx1x22 - rax22)/(2.0d0*rax1*rx1x2)
      if (dabs(ctaint).le.1.0d0) go to 2810
      write (6,*) 'error in taint, ctaint = ',ctaint
      taint = 0.0d0
      if (ctaint.lt.1.0d0) taint = 2.0d0*dasin(1.0d0)
      go to 2820
 2810 continue
c      write (6,*) 'ctaint test',ctaint
      taint = dacos(ctaint)
 2820 continue
 
      sum = 0.0d0
      call cross(rbvec,rax1v,cptmp1)
      call cross(rax1v,rx1x2v,cptmp2)
      do 2830 idim = 1 , ndim
         sum = sum - cptmp1(idim)*cptmp2(idim)
 2830 continue
      ctau3 = sum/(rbond*rax12*rx1x2*sin(the1)*sin(taint))
      if (dabs(ctau3).le.1.0d0) go to 2840
      write (6,*) 'error in tau3, ctau3 = ',ctau3
      tau3 = 0.0d0
      if (ctau3.lt.1.0d0) tau3 = 2.0d0*dasin(1.0d0)
      go to 2850
 2840 continue
c      write (6,*) 'ctau3 test',ctau3
      tau3 = dacos(ctau3)
 2850 continue
         
c     
c     still need to determine whether it is +/-
c     
         
      sum = 0.0d0
      call cross(rbvec,rax1v,cptmp1)
      call cross(rax1v,cptmp2,cptmp3)
      do 2860 idim = 1 , ndim
         sum = sum - cptmp1(idim)*cptmp3(idim)
 2860 continue
      ctau3b = sum/(rbond*rax12*rax1*rx1x2*sin(the1)*
     $ sin(taint))
      if (dabs(ctau3b).le.1.0d0) go to 2870
      write (6,*) 'error in tau3b, ctau3b = ',ctau3b
      tau3b = 0.0d0
      if (ctau3b.lt.1.0d0) tau3b = 2.0d0*dasin(1.0d0)
      go to 2880
 2870 continue
c      write (6,*) 'ctau3b test',ctau3b
      tau3b = dacos(ctau3b)
 2880 continue
            
      if (tau3b.lt.pi2) tau3 = -tau3

      thetap = 180.0d0 - thep*180.0d0/pi
      phi = tau3*180.0d0/pi

      if (rmh.lt.3.2d0) then
c         write (6,*) 'error rmh too small',rmh,theta,phi
         rmh=3.2d0
      endif
      if (rmh.gt.15.0d0) then
c         write (6,*) 'error rmh too large',rmh,theta,phi
         rmh=15.0d0
      endif
c     write (6,*) 'entering hpch3o',rmh,thetap,phi,vtot
      call hpch3o(rmh,thetap,phi,vtot)
c     write (6,*) 'exiting hpch3o',rmh,thetap,phi,vtot
c add in correction from geometry optimization
      regrlx(9) = 3.8/0.529
      regrlx(8) = 3.6/0.529
      regrlx(7) = 3.4/0.529
      regrlx(6) = 3.2/0.529
      regrlx(5) = 3.0/0.529
      regrlx(4) = 2.8/0.529
      regrlx(3) = 2.6/0.529
      regrlx(2) = 2.4/0.529
      regrlx(1) = 2.2/0.529
      egrlxr(9) = .0013602
      egrlxr(8) = .0013404
      egrlxr(7) = .0012822
      egrlxr(6) = .0011844
      egrlxr(5) = .0009623
      egrlxr(4) = .0014367
      egrlxr(3) = .0014656
      egrlxr(2) = .0015279
      egrlxr(1) = .0016330
      if (rbond.gt.7.183) then
         egrlx = .0013602+(.0013602-.0013602)*(rbond-7.183)/
     $    (0.2/0.529)
         if (egrlx.gt..001398) egrlx=.001398
      else
         if (rbond.lt.4.1588) then
            egrlx=0.0016330
         else
            call spline(regrlx,egrlxr,9,1.d32,1.d32,y2)
            call splint(regrlx,egrlxr,y2,9,rbond,egrlx)
         endif
      endif
      egrlx = .001398 - egrlx
c     write (6,*) 'egrlx test',egrlx,rbond
      ch3o_h_rlx = vtot + egrlx

      return
      end

c**************************************************

      real*8 function ch3o_h(r, rpar, ipar)

c     
c     this is a general subroutine for evaluating the potential 
c     energy for the transitional modes.
c
c the potential should be output in au in vtot
c the cartesian coordinates of the atoms of each fragment are 
c stored in r
c the bond length or some other suitable approximation to a 
c distinguished reaction coordinate should be passed back in rmepi
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
c the valence stretching and bending potential works in terms  of
c either 
c (i) an exponential extrapolation of an equilibrium set 
c of force constants 
c or
c (ii) an exponential interpolation of force constants along the 
c bond stretching path.
c
c for the ion-dipole component some work needs to be done on 
c defining the axis for the dipole.  once this is done it should 
c also be easy to include the ion-quadrupole interaction.
c     

      implicit double precision (a-h,o-z)

      include 'param2.fi'
      include 'commonpot.fi'

c
c vectors for primary coordinates 
c

      dimension r(nfrag,natommx,ndim)
      dimension cm(nfrag,ndim),rij(natommx,natommx)
      dimension qq(5),rmov(ndim),rmhv(ndim)

c
c vectors from mcpot calling routine
c

      dimension ift(nfrag),natom(nfrag),ineutf(nchmx)
      dimension dmss(nfrag,natommx),dmsst(nfrag)

c
c intermediate use vectors
c

      dimension r1(ndim),r2p(ndim),r2d(ndim),r2q(ndim)
      dimension nr(nfrag)
      dimension rbvec(ndim),rax1v(ndim),rbx1v(ndim),rax2v(ndim),
     $          rbx2v(ndim),rx1x2v(ndim),ray1v(ndim),rby1v(ndim),
     $          ray2v(ndim),rby2v(ndim),ry1y2v(ndim)
      dimension cptmp1(ndim),cptmp2(ndim),cptmp3(ndim)

c
c specific potential parameter vectors
c

      dimension fijrb(5,5),qqerb(5)
      dimension rpathv(nrpotm),vrr(nrpotm)

      include 'data.fi'

      data dmss(1,1), dmss(1,2) /15.99491 , 12./
c
c for generality/ease of modification purposes i have set this up to 
c begin by evaluating the interfragment atom-atom distances
c and a particular set of valence stretching, bending and torsional 
c angles.
c for efficient use one might wish to only evaluate the components 
c that are required for the particular potential being used.  for a 
c variety of reasons this seems to me to be best left for the specific 
c applications.
c      
      pi2 = asin(1.0d0)
      pi = 2.0d0*pi2

c
c now evaluate atom-atom separations
c

      rij(1,1) = dsqrt((r(1,1,1)-r(2,1,1))**2+(r(1,1,2)-r(2,1,2))**2+ 
     $                 (r(1,1,3)-r(2,1,3))**2)
      rbond = rij(1,1)
      rmepi = rbond

c
c now proceed to the valence bending/torsional coordinates.
c

c
c determine distances needed for bending angles
c
      rbond2 = 0.0d0
      rax12 = 0.0d0
      rbx12 = 0.0d0
      rby12 = 0.0d0
      ray12 = 0.0d0
      rmo2 = 0.0d0
      rmh2 = 0.0d0
      dm12t = dmss(1,1)+dmss(1,2)
      do 1500 idim = 1 , ndim
         rbvec(idim) = r(2,1,idim) - r(1,1,idim)
         rbond2 = rbond2 + rbvec(idim)**2
         rax1v(idim) = r(1,2,idim) - r(1,1,idim)
         rax12 = rax12 + rax1v(idim)**2
         rbx1v(idim) = r(2,1,idim) - r(1,2,idim)
         rbx12 = rbx12 + rbx1v(idim)**2
         rmov(idim)=(dmss(1,2)*r(1,2,idim)+dmss(1,1)*r(1,1,idim))/
     $    dm12t - r(1,1,idim)
         rmo2 = rmo2 + rmov(idim)**2
         rmhv(idim)=(dmss(1,2)*r(1,2,idim)+dmss(1,1)*r(1,1,idim))/
     $    dm12t - r(2,1,idim)
         rmh2 = rmh2 + rmhv(idim)**2
 1500 continue
      rbond = sqrt(rbond2)
      rax1 = sqrt(rax12)
      rbx1 = sqrt(rbx12)
      ray1 = sqrt(ray12)
      rby1 = sqrt(rby12)
      rmo = sqrt(rmo2)
      rmh = sqrt(rmh2)
c     write (6,*) 'rmo test',rmov,rmo
c     write (6,*) 'rmh test',rmhv,rmh
c     write (6,*) 'rbvec test',rbvec,rbond
c     write (6,*) 'r11 test',(r(1,1,idim),idim=1,ndim)
c     write (6,*) 'r12 test',(r(1,2,idim),idim=1,ndim)
c     write (6,*) 'r21 test',(r(2,1,idim),idim=1,ndim)
c     write (6,*) 'mss test',dmss(1,1),dmss(1,2),dm12t

c     
c     calculate the bending angles.
c     

      cthep = (rmo2+rmh2-rbond2)/(2.0d0*rmo*rmh)
      if (dabs(cthep).le.1.0d0) go to 1630
      if(dabs(cthep)-1.d0.gt.1.d-10) 
     x   write (6,*) 'error in thep, cthep = ',cthep
      thep = 0.0d0
      if (cthep.lt.1.0d0) thep = 2.0d0*dasin(1.0d0)
      go to 1650
 1630 continue
      thep = dacos(cthep)
 1650 continue

      cthe1 = (rax12 + rbond2 - rbx12)/(2.0d0*rax1*rbond)
      if (dabs(cthe1).le.1.0d0) go to 2030
      if(dabs(cthe1)-1.d0.gt.1.e-6) 
     x   write (6,*) 'error in the1, cthe1 = ',cthe1
      the1 = 0.0d0
      if (cthe1.lt.1.0d0) the1 = 2.0d0*dasin(1.0d0)
      go to 2050
 2030 continue
c      write (6,*) 'cthe1 test',cthe1
      the1 = dacos(cthe1)
 2050 continue
      
c     
c     calculate the torsional angles if ivb.ne.0
c     
c
c
c determine distances needed for torsional angles
c
      rax22 = 0.0d0
      rbx22 = 0.0d0
      rby22 = 0.0d0
      ray22 = 0.0d0
      ry1y22 = 0.0d0
      rx1x22 = 0.0d0
      do idim = 1 , ndim
         rax2v(idim) = r(1,3,idim) - r(1,1,idim)
         rax22 = rax22 + rax2v(idim)**2
         rbx2v(idim) = r(2,1,idim) - r(1,3,idim)
         rbx22 = rbx22 + rbx2v(idim)**2
         rx1x2v(idim) = r(1,3,idim) - r(1,2,idim)
         rx1x22 = rx1x22 + rx1x2v(idim)**2
      enddo
      rax2 = sqrt(rax22)
      rbx2 = sqrt(rbx22)
      ray2 = sqrt(ray22)
      rby2 = sqrt(rby22)
      ry1y2 = sqrt(ry1y22)
      rx1x2 = sqrt(rx1x22)

c
c     now evaluate tau 3 (torsion angle for x2,11,x1,21.
c     tau3b corresponds to torsion angle for perp to x1,21;
c     x2,21,x1,perp
c

c
c     consider the terminal atom case.
c
c     first evaluate internal x2,11,x1 bending angle
c

      ctaint = (rax12 + rx1x22 - rax22)/(2.0d0*rax1*rx1x2)
      if (dabs(ctaint).le.1.0d0) go to 2810
      write (6,*) 'error in taint, ctaint = ',ctaint
      taint = 0.0d0
      if (ctaint.lt.1.0d0) taint = 2.0d0*dasin(1.0d0)
      go to 2820
 2810 continue
c      write (6,*) 'ctaint test',ctaint
      taint = dacos(ctaint)
 2820 continue
 
      sum = 0.0d0
      call cross(rbvec,rax1v,cptmp1)
      call cross(rax1v,rx1x2v,cptmp2)
      do 2830 idim = 1 , ndim
         sum = sum - cptmp1(idim)*cptmp2(idim)
 2830 continue
      ctau3 = sum/(rbond*rax12*rx1x2*sin(the1)*sin(taint))
      if (dabs(ctau3).le.1.0d0) go to 2840
      write (6,*) 'error in tau3, ctau3 = ',ctau3
      tau3 = 0.0d0
      if (ctau3.lt.1.0d0) tau3 = 2.0d0*dasin(1.0d0)
      go to 2850
 2840 continue
c      write (6,*) 'ctau3 test',ctau3
      tau3 = dacos(ctau3)
 2850 continue
         
c     
c     still need to determine whether it is +/-
c     
         
      sum = 0.0d0
      call cross(rbvec,rax1v,cptmp1)
      call cross(rax1v,cptmp2,cptmp3)
      do 2860 idim = 1 , ndim
         sum = sum - cptmp1(idim)*cptmp3(idim)
 2860 continue
      ctau3b = sum/(rbond*rax12*rax1*rx1x2*sin(the1)*
     $ sin(taint))
      if (dabs(ctau3b).le.1.0d0) go to 2870
      write (6,*) 'error in tau3b, ctau3b = ',ctau3b
      tau3b = 0.0d0
      if (ctau3b.lt.1.0d0) tau3b = 2.0d0*dasin(1.0d0)
      go to 2880
 2870 continue
c      write (6,*) 'ctau3b test',ctau3b
      tau3b = dacos(ctau3b)
 2880 continue
            
      if (tau3b.lt.pi2) tau3 = -tau3

      thetap = 180.0d0-thep*180.0d0/pi
      phi = tau3*180.0d0/pi

      if (rmh.lt.3.2d0) then
c         write (6,*) 'error rmh too small',rmh,theta,phi
         rmh=3.2d0
      endif
      if (rmh.gt.15.0d0) then
c         write (6,*) 'error rmh too large',rmh,theta,phi
         rmh=15.0d0
      endif
      call hpch3o(rmh,thetap,phi,vtot)
      ch3o_h = vtot
c     write (6,*) 'pot test',rmh,thetap,phi,vtot*cautoicm/350.

      return
      end

c                hpch3o input:
c
c                rmh   = r(m-ha)   (au)
c                        m is the center of mass of the co bond
c
c                theta = polar angle  (degrees)
c                        angle between m-ha axis and the
c                        co bond. theta=0 correspnds to
c                        approach along the c end of the co axis.
c
c                phi   = azimuthal angle (degrees)
c                        dihedral angle between the plane
c                        defined by the co bond and
c                        the ha atom and the plane defined
c                        by the co bond and one of the
c                        inactive hydrogens
c
c
c         ranges:
c
c                       3.2 < rmh   < 15.0 au
c                       0.0 < theta < 180.0
c                        phi is periodic 
c
c     output: energy in au relative to the h+ch3o asymptote
c
c     energies calculated at the cas+1+2/aug-pvdz level
c
c     the geometry of the methoxy fragment is fixed at the
c      following (c3v) geometry:
c
c       r(co) = 2.60 au
c       r(ch) = 2.10 au
c       hco angles = 110.0
c       hhco dihedral angles = 120.0
c

      subroutine hpch3o(rmh,theta,phi,energy)
      implicit real*8 (a-h,o-z)
      data degrad / 57.29577951308232d 00/
c
c     fcs is based on points having cs symmetry.
c     by convention if theta>0, fcs returns the
c     energy for phi=0 (eclipsed), if theta<0, fcs
c     returns the energy for phi=180 (staggered)
c
c
      call fcs_ci_1ap(rmh,-theta,es)
      call fcs_ci_1ap(rmh,theta,ee)

      p = phi/degrad
      sum = 0.5d0*(ee+es)
      diff =0.5d0*(ee-es)
      energy = 115.2465333d0 + sum + diff*cos( 3.0d0*p) 

      return
      end
      subroutine fcs_ci_1ap(xi,yi,fi)
c
c     cas+1+2/aug-cc-pvdz
c     1a' state
c
      implicit real*8 (a-h,o-z)
      dimension fpp(15,42,2),f(15,42),fpppp(15,42)
      dimension delx(14),dely(41),x(15),y(42)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -1.15348328d+02 , -1.15351647d+02 /
      data f( 1, 3),f( 1, 4) / -1.15346749d+02 , -1.15338611d+02 /
      data f( 1, 5),f( 1, 6) / -1.15338629d+02 , -1.15346920d+02 /
      data f( 1, 7),f( 1, 8) / -1.15352260d+02 , -1.15349851d+02 /
      data f( 1, 9),f( 1,10) / -1.15338976d+02 , -1.15321422d+02 /
      data f( 1,11),f( 1,12) / -1.15299617d+02 , -1.15275576d+02 /
      data f( 1,13),f( 1,14) / -1.15250587d+02 , -1.15225233d+02 /
      data f( 1,15),f( 1,16) / -1.15199327d+02 , -1.15171967d+02 /
      data f( 1,17),f( 1,18) / -1.15148837d+02 , -1.15200728d+02 /
      data f( 1,19),f( 1,20) / -1.15200593d+02 , -1.15178523d+02 /
      data f( 1,21),f( 1,22) / -1.15142677d+02 , -1.15140934d+02 /
      data f( 1,23),f( 1,24) / -1.15143708d+02 , -1.14866162d+02 /
      data f( 1,25),f( 1,26) / -1.12849565d+02 , -1.14495165d+02 /
      data f( 1,27),f( 1,28) / -1.15123574d+02 , -1.15174983d+02 /
      data f( 1,29),f( 1,30) / -1.15192000d+02 , -1.15226218d+02 /
      data f( 1,31),f( 1,32) / -1.15259596d+02 , -1.15289739d+02 /
      data f( 1,33),f( 1,34) / -1.15315644d+02 , -1.15335839d+02 /
      data f( 1,35),f( 1,36) / -1.15348328d+02 , -1.15351647d+02 /
      data f( 1,37),f( 1,38) / -1.15346749d+02 , -1.15338611d+02 /
      data f( 1,39),f( 1,40) / -1.15338629d+02 , -1.15346920d+02 /
      data f( 1,41),f( 1,42) / -1.15352260d+02 , -1.15349851d+02 /
      data f( 2, 1),f( 2, 2) / -1.15327047d+02 , -1.15328502d+02 /
      data f( 2, 3),f( 2, 4) / -1.15322315d+02 , -1.15313252d+02 /
      data f( 2, 5),f( 2, 6) / -1.15313270d+02 , -1.15322498d+02 /
      data f( 2, 7),f( 2, 8) / -1.15329136d+02 , -1.15328558d+02 /
      data f( 2, 9),f( 2,10) / -1.15320261d+02 , -1.15305973d+02 /
      data f( 2,11),f( 2,12) / -1.15287936d+02 , -1.15268016d+02 /
      data f( 2,13),f( 2,14) / -1.15247464d+02 , -1.15226957d+02 /
      data f( 2,15),f( 2,16) / -1.15206617d+02 , -1.15186383d+02 /
      data f( 2,17),f( 2,18) / -1.15170556d+02 , -1.15185900d+02 /
      data f( 2,19),f( 2,20) / -1.15207081d+02 , -1.15197721d+02 /
      data f( 2,21),f( 2,22) / -1.15174217d+02 , -1.15174233d+02 /
      data f( 2,23),f( 2,24) / -1.15179906d+02 , -1.15006574d+02 /
      data f( 2,25),f( 2,26) / -1.14332296d+02 , -1.14813765d+02 /
      data f( 2,27),f( 2,28) / -1.15159189d+02 , -1.15192473d+02 /
      data f( 2,29),f( 2,30) / -1.15200765d+02 , -1.15225952d+02 /
      data f( 2,31),f( 2,32) / -1.15253356d+02 , -1.15278688d+02 /
      data f( 2,33),f( 2,34) / -1.15300484d+02 , -1.15317232d+02 /
      data f( 2,35),f( 2,36) / -1.15327047d+02 , -1.15328502d+02 /
      data f( 2,37),f( 2,38) / -1.15322315d+02 , -1.15313252d+02 /
      data f( 2,39),f( 2,40) / -1.15313270d+02 , -1.15322498d+02 /
      data f( 2,41),f( 2,42) / -1.15329136d+02 , -1.15328558d+02 /
      data f( 3, 1),f( 3, 2) / -1.15307838d+02 , -1.15306815d+02 /
      data f( 3, 3),f( 3, 4) / -1.15298031d+02 , -1.15285907d+02 /
      data f( 3, 5),f( 3, 6) / -1.15285925d+02 , -1.15298223d+02 /
      data f( 3, 7),f( 3, 8) / -1.15307442d+02 , -1.15309274d+02 /
      data f( 3, 9),f( 3,10) / -1.15303731d+02 , -1.15292650d+02 /
      data f( 3,11),f( 3,12) / -1.15278197d+02 , -1.15262148d+02 /
      data f( 3,13),f( 3,14) / -1.15245709d+02 , -1.15229577d+02 /
      data f( 3,15),f( 3,16) / -1.15214000d+02 , -1.15199101d+02 /
      data f( 3,17),f( 3,18) / -1.15187300d+02 , -1.15188883d+02 /
      data f( 3,19),f( 3,20) / -1.15203736d+02 , -1.15201536d+02 /
      data f( 3,21),f( 3,22) / -1.15183459d+02 , -1.15185311d+02 /
      data f( 3,23),f( 3,24) / -1.15200119d+02 , -1.15110906d+02 /
      data f( 3,25),f( 3,26) / -1.14851248d+02 , -1.15020720d+02 /
      data f( 3,27),f( 3,28) / -1.15190420d+02 , -1.15206184d+02 /
      data f( 3,29),f( 3,30) / -1.15209569d+02 , -1.15227625d+02 /
      data f( 3,31),f( 3,32) / -1.15249294d+02 , -1.15269898d+02 /
      data f( 3,33),f( 3,34) / -1.15287646d+02 , -1.15300923d+02 /
      data f( 3,35),f( 3,36) / -1.15307838d+02 , -1.15306815d+02 /
      data f( 3,37),f( 3,38) / -1.15298031d+02 , -1.15285907d+02 /
      data f( 3,39),f( 3,40) / -1.15285925d+02 , -1.15298223d+02 /
      data f( 3,41),f( 3,42) / -1.15307442d+02 , -1.15309274d+02 /
      data f( 4, 1),f( 4, 2) / -1.15291650d+02 , -1.15288381d+02 /
      data f( 4, 3),f( 4, 4) / -1.15276895d+02 , -1.15259600d+02 /
      data f( 4, 5),f( 4, 6) / -1.15259615d+02 , -1.15277082d+02 /
      data f( 4, 7),f( 4, 8) / -1.15288964d+02 , -1.15292950d+02 /
      data f( 4, 9),f( 4,10) / -1.15289839d+02 , -1.15281606d+02 /
      data f( 4,11),f( 4,12) / -1.15270349d+02 , -1.15257745d+02 /
      data f( 4,13),f( 4,14) / -1.15244906d+02 , -1.15232461d+02 /
      data f( 4,15),f( 4,16) / -1.15220660d+02 , -1.15209599d+02 /
      data f( 4,17),f( 4,18) / -1.15200532d+02 , -1.15197688d+02 /
      data f( 4,19),f( 4,20) / -1.15201613d+02 , -1.15198100d+02 /
      data f( 4,21),f( 4,22) / -1.15178061d+02 , -1.15183174d+02 /
      data f( 4,23),f( 4,24) / -1.15209689d+02 , -1.15175016d+02 /
      data f( 4,25),f( 4,26) / -1.15071249d+02 , -1.15136590d+02 /
      data f( 4,27),f( 4,28) / -1.15213175d+02 , -1.15216539d+02 /
      data f( 4,29),f( 4,30) / -1.15217303d+02 , -1.15230198d+02 /
      data f( 4,31),f( 4,32) / -1.15246892d+02 , -1.15263187d+02 /
      data f( 4,33),f( 4,34) / -1.15277222d+02 , -1.15287342d+02 /
      data f( 4,35),f( 4,36) / -1.15291650d+02 , -1.15288381d+02 /
      data f( 4,37),f( 4,38) / -1.15276895d+02 , -1.15259600d+02 /
      data f( 4,39),f( 4,40) / -1.15259615d+02 , -1.15277082d+02 /
      data f( 4,41),f( 4,42) / -1.15288964d+02 , -1.15292950d+02 /
      data f( 5, 1),f( 5, 2) / -1.15278807d+02 , -1.15274104d+02 /
      data f( 5, 3),f( 5, 4) / -1.15261230d+02 , -1.15239959d+02 /
      data f( 5, 5),f( 5, 6) / -1.15239958d+02 , -1.15261385d+02 /
      data f( 5, 7),f( 5, 8) / -1.15274604d+02 , -1.15279925d+02 /
      data f( 5, 9),f( 5,10) / -1.15278653d+02 , -1.15272760d+02 /
      data f( 5,11),f( 5,12) / -1.15264210d+02 , -1.15254534d+02 /
      data f( 5,13),f( 5,14) / -1.15244704d+02 , -1.15235234d+02 /
      data f( 5,15),f( 5,16) / -1.15226327d+02 , -1.15218064d+02 /
      data f( 5,17),f( 5,18) / -1.15211016d+02 , -1.15206920d+02 /
      data f( 5,19),f( 5,20) / -1.15205382d+02 , -1.15198598d+02 /
      data f( 5,21),f( 5,22) / -1.15179669d+02 , -1.15183217d+02 /
      data f( 5,23),f( 5,24) / -1.15213999d+02 , -1.15210237d+02 /
      data f( 5,25),f( 5,26) / -1.15171587d+02 , -1.15197589d+02 /
      data f( 5,27),f( 5,28) / -1.15227781d+02 , -1.15224161d+02 /
      data f( 5,29),f( 5,30) / -1.15223732d+02 , -1.15233000d+02 /
      data f( 5,31),f( 5,32) / -1.15245643d+02 , -1.15258236d+02 /
      data f( 5,33),f( 5,34) / -1.15269058d+02 , -1.15276522d+02 /
      data f( 5,35),f( 5,36) / -1.15278807d+02 , -1.15274104d+02 /
      data f( 5,37),f( 5,38) / -1.15261230d+02 , -1.15239959d+02 /
      data f( 5,39),f( 5,40) / -1.15239958d+02 , -1.15261385d+02 /
      data f( 5,41),f( 5,42) / -1.15274604d+02 , -1.15279925d+02 /
      data f( 6, 1),f( 6, 2) / -1.15269139d+02 , -1.15263986d+02 /
      data f( 6, 3),f( 6, 4) / -1.15251740d+02 , -1.15232961d+02 /
      data f( 6, 5),f( 6, 6) / -1.15232932d+02 , -1.15251842d+02 /
      data f( 6, 7),f( 6, 8) / -1.15264381d+02 , -1.15270053d+02 /
      data f( 6, 9),f( 6,10) / -1.15269972d+02 , -1.15265885d+02 /
      data f( 6,11),f( 6,12) / -1.15259530d+02 , -1.15252244d+02 /
      data f( 6,13),f( 6,14) / -1.15244839d+02 , -1.15237705d+02 /
      data f( 6,15),f( 6,16) / -1.15230992d+02 , -1.15224767d+02 /
      data f( 6,17),f( 6,18) / -1.15219288d+02 , -1.15215192d+02 /
      data f( 6,19),f( 6,20) / -1.15211741d+02 , -1.15204364d+02 /
      data f( 6,21),f( 6,22) / -1.15190165d+02 , -1.15191497d+02 /
      data f( 6,23),f( 6,24) / -1.15217399d+02 , -1.15228091d+02 /
      data f( 6,25),f( 6,26) / -1.15217861d+02 , -1.15228041d+02 /
      data f( 6,27),f( 6,28) / -1.15236226d+02 , -1.15229737d+02 /
      data f( 6,29),f( 6,30) / -1.15228938d+02 , -1.15235663d+02 /
      data f( 6,31),f( 6,32) / -1.15245133d+02 , -1.15254690d+02 /
      data f( 6,33),f( 6,34) / -1.15262858d+02 , -1.15268220d+02 /
      data f( 6,35),f( 6,36) / -1.15269139d+02 , -1.15263986d+02 /
      data f( 6,37),f( 6,38) / -1.15251740d+02 , -1.15232961d+02 /
      data f( 6,39),f( 6,40) / -1.15232932d+02 , -1.15251842d+02 /
      data f( 6,41),f( 6,42) / -1.15264381d+02 , -1.15270053d+02 /
      data f( 7, 1),f( 7, 2) / -1.15259501d+02 , -1.15254993d+02 /
      data f( 7, 3),f( 7, 4) / -1.15245888d+02 , -1.15234467d+02 /
      data f( 7, 5),f( 7, 6) / -1.15234427d+02 , -1.15245926d+02 /
      data f( 7, 7),f( 7, 8) / -1.15255241d+02 , -1.15260128d+02 /
      data f( 7, 9),f( 7,10) / -1.15260855d+02 , -1.15258599d+02 /
      data f( 7,11),f( 7,12) / -1.15254655d+02 , -1.15250027d+02 /
      data f( 7,13),f( 7,14) / -1.15245298d+02 , -1.15240707d+02 /
      data f( 7,15),f( 7,16) / -1.15236332d+02 , -1.15232225d+02 /
      data f( 7,17),f( 7,18) / -1.15228473d+02 , -1.15225135d+02 /
      data f( 7,19),f( 7,20) / -1.15221503d+02 , -1.15215457d+02 /
      data f( 7,21),f( 7,22) / -1.15206678d+02 , -1.15206676d+02 /
      data f( 7,23),f( 7,24) / -1.15223378d+02 , -1.15238790d+02 /
      data f( 7,25),f( 7,26) / -1.15243445d+02 , -1.15245514d+02 /
      data f( 7,27),f( 7,28) / -1.15241970d+02 , -1.15235518d+02 /
      data f( 7,29),f( 7,30) / -1.15234815d+02 , -1.15239053d+02 /
      data f( 7,31),f( 7,32) / -1.15245112d+02 , -1.15251273d+02 /
      data f( 7,33),f( 7,34) / -1.15256468d+02 , -1.15259623d+02 /
      data f( 7,35),f( 7,36) / -1.15259501d+02 , -1.15254993d+02 /
      data f( 7,37),f( 7,38) / -1.15245888d+02 , -1.15234467d+02 /
      data f( 7,39),f( 7,40) / -1.15234427d+02 , -1.15245926d+02 /
      data f( 7,41),f( 7,42) / -1.15255241d+02 , -1.15260128d+02 /
      data f( 8, 1),f( 8, 2) / -1.15251685d+02 , -1.15249148d+02 /
      data f( 8, 3),f( 8, 4) / -1.15244849d+02 , -1.15240407d+02 /
      data f( 8, 5),f( 8, 6) / -1.15240383d+02 , -1.15244852d+02 /
      data f( 8, 7),f( 8, 8) / -1.15249253d+02 , -1.15251985d+02 /
      data f( 8, 9),f( 8,10) / -1.15252743d+02 , -1.15251960d+02 /
      data f( 8,11),f( 8,12) / -1.15250270d+02 , -1.15248199d+02 /
      data f( 8,13),f( 8,14) / -1.15246057d+02 , -1.15243948d+02 /
      data f( 8,15),f( 8,16) / -1.15241873d+02 , -1.15239833d+02 /
      data f( 8,17),f( 8,18) / -1.15237850d+02 , -1.15235873d+02 /
      data f( 8,19),f( 8,20) / -1.15233546d+02 , -1.15230226d+02 /
      data f( 8,21),f( 8,22) / -1.15226241d+02 , -1.15225805d+02 /
      data f( 8,23),f( 8,24) / -1.15232998d+02 , -1.15242646d+02 /
      data f( 8,25),f( 8,26) / -1.15248724d+02 , -1.15248729d+02 /
      data f( 8,27),f( 8,28) / -1.15244647d+02 , -1.15241223d+02 /
      data f( 8,29),f( 8,30) / -1.15240954d+02 , -1.15242902d+02 /
      data f( 8,31),f( 8,32) / -1.15245705d+02 , -1.15248546d+02 /
      data f( 8,33),f( 8,34) / -1.15250874d+02 , -1.15252126d+02 /
      data f( 8,35),f( 8,36) / -1.15251685d+02 , -1.15249148d+02 /
      data f( 8,37),f( 8,38) / -1.15244849d+02 , -1.15240407d+02 /
      data f( 8,39),f( 8,40) / -1.15240383d+02 , -1.15244852d+02 /
      data f( 8,41),f( 8,42) / -1.15249253d+02 , -1.15251985d+02 /
      data f( 9, 1),f( 9, 2) / -1.15248736d+02 , -1.15247549d+02 /
      data f( 9, 3),f( 9, 4) / -1.15245701d+02 , -1.15243931d+02 /
      data f( 9, 5),f( 9, 6) / -1.15243921d+02 , -1.15245703d+02 /
      data f( 9, 7),f( 9, 8) / -1.15247598d+02 , -1.15248876d+02 /
      data f( 9, 9),f( 9,10) / -1.15249306d+02 , -1.15249031d+02 /
      data f( 9,11),f( 9,12) / -1.15248325d+02 , -1.15247434d+02 /
      data f( 9,13),f( 9,14) / -1.15246507d+02 , -1.15245597d+02 /
      data f( 9,15),f( 9,16) / -1.15244683d+02 , -1.15243733d+02 /
      data f( 9,17),f( 9,18) / -1.15242735d+02 , -1.15241670d+02 /
      data f( 9,19),f( 9,20) / -1.15240433d+02 , -1.15238852d+02 /
      data f( 9,21),f( 9,22) / -1.15237114d+02 , -1.15236767d+02 /
      data f( 9,23),f( 9,24) / -1.15239695d+02 , -1.15244027d+02 /
      data f( 9,25),f( 9,26) / -1.15247168d+02 , -1.15247410d+02 /
      data f( 9,27),f( 9,28) / -1.15245707d+02 , -1.15244274d+02 /
      data f( 9,29),f( 9,30) / -1.15244150d+02 , -1.15244993d+02 /
      data f( 9,31),f( 9,32) / -1.15246233d+02 , -1.15247493d+02 /
      data f( 9,33),f( 9,34) / -1.15248510d+02 , -1.15249015d+02 /
      data f( 9,35),f( 9,36) / -1.15248736d+02 , -1.15247549d+02 /
      data f( 9,37),f( 9,38) / -1.15245701d+02 , -1.15243931d+02 /
      data f( 9,39),f( 9,40) / -1.15243921d+02 , -1.15245703d+02 /
      data f( 9,41),f( 9,42) / -1.15247598d+02 , -1.15248876d+02 /
      data f(10, 1),f(10, 2) / -1.15247583d+02 , -1.15247072d+02 /
      data f(10, 3),f(10, 4) / -1.15246302d+02 , -1.15245585d+02 /
      data f(10, 5),f(10, 6) / -1.15245582d+02 , -1.15246307d+02 /
      data f(10, 7),f(10, 8) / -1.15247100d+02 , -1.15247653d+02 /
      data f(10, 9),f(10,10) / -1.15247850d+02 , -1.15247743d+02 /
      data f(10,11),f(10,12) / -1.15247450d+02 , -1.15247079d+02 /
      data f(10,13),f(10,14) / -1.15246697d+02 , -1.15246331d+02 /
      data f(10,15),f(10,16) / -1.15245960d+02 , -1.15245581d+02 /
      data f(10,17),f(10,18) / -1.15245137d+02 , -1.15244622d+02 /
      data f(10,19),f(10,20) / -1.15244037d+02 , -1.15243361d+02 /
      data f(10,21),f(10,22) / -1.15242661d+02 , -1.15242450d+02 /
      data f(10,23),f(10,24) / -1.15243582d+02 , -1.15245338d+02 /
      data f(10,25),f(10,26) / -1.15246726d+02 , -1.15247008d+02 /
      data f(10,27),f(10,28) / -1.15246423d+02 , -1.15245813d+02 /
      data f(10,29),f(10,30) / -1.15245697d+02 , -1.15246015d+02 /
      data f(10,31),f(10,32) / -1.15246526d+02 , -1.15247060d+02 /
      data f(10,33),f(10,34) / -1.15247497d+02 , -1.15247710d+02 /
      data f(10,35),f(10,36) / -1.15247583d+02 , -1.15247072d+02 /
      data f(10,37),f(10,38) / -1.15246302d+02 , -1.15245585d+02 /
      data f(10,39),f(10,40) / -1.15245582d+02 , -1.15246307d+02 /
      data f(10,41),f(10,42) / -1.15247100d+02 , -1.15247653d+02 /
      data f(11, 1),f(11, 2) / -1.15246867d+02 , -1.15246787d+02 /
      data f(11, 3),f(11, 4) / -1.15246660d+02 , -1.15246541d+02 /
      data f(11, 5),f(11, 6) / -1.15246541d+02 , -1.15246663d+02 /
      data f(11, 7),f(11, 8) / -1.15246797d+02 , -1.15246890d+02 /
      data f(11, 9),f(11,10) / -1.15246920d+02 , -1.15246895d+02 /
      data f(11,11),f(11,12) / -1.15246838d+02 , -1.15246774d+02 /
      data f(11,13),f(11,14) / -1.15246719d+02 , -1.15246682d+02 /
      data f(11,15),f(11,16) / -1.15246661d+02 , -1.15246643d+02 /
      data f(11,17),f(11,18) / -1.15246606d+02 , -1.15246539d+02 /
      data f(11,19),f(11,20) / -1.15246464d+02 , -1.15246397d+02 /
      data f(11,21),f(11,22) / -1.15246321d+02 , -1.15246261d+02 /
      data f(11,23),f(11,24) / -1.15246409d+02 , -1.15246673d+02 /
      data f(11,25),f(11,26) / -1.15246922d+02 , -1.15247008d+02 /
      data f(11,27),f(11,28) / -1.15246904d+02 , -1.15246739d+02 /
      data f(11,29),f(11,30) / -1.15246638d+02 , -1.15246626d+02 /
      data f(11,31),f(11,32) / -1.15246675d+02 , -1.15246754d+02 /
      data f(11,33),f(11,34) / -1.15246832d+02 , -1.15246879d+02 /
      data f(11,35),f(11,36) / -1.15246867d+02 , -1.15246787d+02 /
      data f(11,37),f(11,38) / -1.15246660d+02 , -1.15246541d+02 /
      data f(11,39),f(11,40) / -1.15246541d+02 , -1.15246663d+02 /
      data f(11,41),f(11,42) / -1.15246797d+02 , -1.15246890d+02 /
      data f(12, 1),f(12, 2) / -1.15246678d+02 , -1.15246669d+02 /
      data f(12, 3),f(12, 4) / -1.15246651d+02 , -1.15246633d+02 /
      data f(12, 5),f(12, 6) / -1.15246633d+02 , -1.15246652d+02 /
      data f(12, 7),f(12, 8) / -1.15246672d+02 , -1.15246685d+02 /
      data f(12, 9),f(12,10) / -1.15246687d+02 , -1.15246682d+02 /
      data f(12,11),f(12,12) / -1.15246671d+02 , -1.15246662d+02 /
      data f(12,13),f(12,14) / -1.15246657d+02 , -1.15246653d+02 /
      data f(12,15),f(12,16) / -1.15246677d+02 , -1.15246698d+02 /
      data f(12,17),f(12,18) / -1.15246717d+02 , -1.15246729d+02 /
      data f(12,19),f(12,20) / -1.15246738d+02 , -1.15246754d+02 /
      data f(12,21),f(12,22) / -1.15246767d+02 , -1.15246756d+02 /
      data f(12,23),f(12,24) / -1.15246773d+02 , -1.15246821d+02 /
      data f(12,25),f(12,26) / -1.15246869d+02 , -1.15246872d+02 /
      data f(12,27),f(12,28) / -1.15246827d+02 , -1.15246759d+02 /
      data f(12,29),f(12,30) / -1.15246698d+02 , -1.15246649d+02 /
      data f(12,31),f(12,32) / -1.15246641d+02 , -1.15246645d+02 /
      data f(12,33),f(12,34) / -1.15246660d+02 , -1.15246674d+02 /
      data f(12,35),f(12,36) / -1.15246678d+02 , -1.15246669d+02 /
      data f(12,37),f(12,38) / -1.15246651d+02 , -1.15246633d+02 /
      data f(12,39),f(12,40) / -1.15246633d+02 , -1.15246652d+02 /
      data f(12,41),f(12,42) / -1.15246672d+02 , -1.15246685d+02 /
      data f(13, 1),f(13, 2) / -1.15246574d+02 , -1.15246575d+02 /
      data f(13, 3),f(13, 4) / -1.15246567d+02 , -1.15246567d+02 /
      data f(13, 5),f(13, 6) / -1.15246567d+02 , -1.15246567d+02 /
      data f(13, 7),f(13, 8) / -1.15246575d+02 , -1.15246575d+02 /
      data f(13, 9),f(13,10) / -1.15246575d+02 , -1.15246575d+02 /
      data f(13,11),f(13,12) / -1.15246577d+02 , -1.15246571d+02 /
      data f(13,13),f(13,14) / -1.15246575d+02 , -1.15246579d+02 /
      data f(13,15),f(13,16) / -1.15246584d+02 , -1.15246588d+02 /
      data f(13,17),f(13,18) / -1.15246602d+02 , -1.15246608d+02 /
      data f(13,19),f(13,20) / -1.15246615d+02 , -1.15246621d+02 /
      data f(13,21),f(13,22) / -1.15246615d+02 , -1.15246616d+02 /
      data f(13,23),f(13,24) / -1.15246629d+02 , -1.15246634d+02 /
      data f(13,25),f(13,26) / -1.15246634d+02 , -1.15246625d+02 /
      data f(13,27),f(13,28) / -1.15246611d+02 , -1.15246597d+02 /
      data f(13,29),f(13,30) / -1.15246578d+02 , -1.15246569d+02 /
      data f(13,31),f(13,32) / -1.15246564d+02 , -1.15246562d+02 /
      data f(13,33),f(13,34) / -1.15246571d+02 , -1.15246572d+02 /
      data f(13,35),f(13,36) / -1.15246574d+02 , -1.15246575d+02 /
      data f(13,37),f(13,38) / -1.15246567d+02 , -1.15246567d+02 /
      data f(13,39),f(13,40) / -1.15246567d+02 , -1.15246567d+02 /
      data f(13,41),f(13,42) / -1.15246575d+02 , -1.15246575d+02 /
      data f(14, 1),f(14, 2) / -1.15246542d+02 , -1.15246542d+02 /
      data f(14, 3),f(14, 4) / -1.15246543d+02 , -1.15246543d+02 /
      data f(14, 5),f(14, 6) / -1.15246543d+02 , -1.15246543d+02 /
      data f(14, 7),f(14, 8) / -1.15246542d+02 , -1.15246542d+02 /
      data f(14, 9),f(14,10) / -1.15246542d+02 , -1.15246542d+02 /
      data f(14,11),f(14,12) / -1.15246543d+02 , -1.15246544d+02 /
      data f(14,13),f(14,14) / -1.15246545d+02 , -1.15246547d+02 /
      data f(14,15),f(14,16) / -1.15246548d+02 , -1.15246550d+02 /
      data f(14,17),f(14,18) / -1.15246551d+02 , -1.15246554d+02 /
      data f(14,19),f(14,20) / -1.15246557d+02 , -1.15246560d+02 /
      data f(14,21),f(14,22) / -1.15246561d+02 , -1.15246562d+02 /
      data f(14,23),f(14,24) / -1.15246561d+02 , -1.15246561d+02 /
      data f(14,25),f(14,26) / -1.15246560d+02 , -1.15246557d+02 /
      data f(14,27),f(14,28) / -1.15246553d+02 , -1.15246549d+02 /
      data f(14,29),f(14,30) / -1.15246546d+02 , -1.15246544d+02 /
      data f(14,31),f(14,32) / -1.15246542d+02 , -1.15246541d+02 /
      data f(14,33),f(14,34) / -1.15246541d+02 , -1.15246541d+02 /
      data f(14,35),f(14,36) / -1.15246542d+02 , -1.15246542d+02 /
      data f(14,37),f(14,38) / -1.15246543d+02 , -1.15246543d+02 /
      data f(14,39),f(14,40) / -1.15246543d+02 , -1.15246543d+02 /
      data f(14,41),f(14,42) / -1.15246542d+02 , -1.15246542d+02 /
      data f(15, 1),f(15, 2) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15, 3),f(15, 4) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15, 5),f(15, 6) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15, 7),f(15, 8) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15, 9),f(15,10) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,11),f(15,12) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,13),f(15,14) / -1.15246533d+02 , -1.15246534d+02 /
      data f(15,15),f(15,16) / -1.15246534d+02 , -1.15246535d+02 /
      data f(15,17),f(15,18) / -1.15246535d+02 , -1.15246536d+02 /
      data f(15,19),f(15,20) / -1.15246537d+02 , -1.15246537d+02 /
      data f(15,21),f(15,22) / -1.15246538d+02 , -1.15246538d+02 /
      data f(15,23),f(15,24) / -1.15246538d+02 , -1.15246538d+02 /
      data f(15,25),f(15,26) / -1.15246538d+02 , -1.15246537d+02 /
      data f(15,27),f(15,28) / -1.15246537d+02 , -1.15246536d+02 /
      data f(15,29),f(15,30) / -1.15246535d+02 , -1.15246534d+02 /
      data f(15,31),f(15,32) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,33),f(15,34) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,35),f(15,36) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,37),f(15,38) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,39),f(15,40) / -1.15246533d+02 , -1.15246533d+02 /
      data f(15,41),f(15,42) / -1.15246533d+02 , -1.15246533d+02 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/-2.85419899d-02, 1.33799608d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 7.72249148d-03, 7.83627832d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 7.07957297d-02, 4.56972589d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 1.27009555d-01,-6.66678188d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 1.26857566d-01,-6.72806730d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 7.05539620d-02, 4.01086759d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 8.12175005d-03, 8.39239695d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-2.69536896d-02, 8.91414461d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-4.36092291d-02, 6.74642460d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/-4.94872039d-02, 4.17355698d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-4.98926344d-02, 2.06294748d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/-4.79887913d-02, 9.91853090d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/-4.45522765d-02,-3.44759843d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/-3.80265481d-02, 2.58258628d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/-2.22021522d-02,-6.67838529d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 3.01465149d-02, 3.28609549d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 1.60334052d-01,-1.50152034d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/-7.39183413d-01, 1.17624782d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 5.17851918d-01,-8.19229359d-05/
      data fpp( 1,20,1),fpp( 1,20,2)/ 5.93352735d-01, 4.67561924d-04/
      data fpp( 1,21,1),fpp( 1,21,2)/ 7.67958574d-01,-9.61770760d-04/
      data fpp( 1,22,1),fpp( 1,22,2)/ 7.88916879d-01, 1.60229332d-03/
      data fpp( 1,23,1),fpp( 1,23,2)/ 5.33248584d-01,-7.90896039d-03/
      data fpp( 1,24,1),fpp( 1,24,2)/ 8.14174287d-01, 4.68527602d-02/
      data fpp( 1,25,1),fpp( 1,25,2)/ 4.02411297d+01,-7.51590205d-02/
      data fpp( 1,26,1),fpp( 1,26,2)/ 3.32409129d+00, 3.40515499d-02/
      data fpp( 1,27,1),fpp( 1,27,2)/ 1.14432977d-02,-1.57671613d-05/
      data fpp( 1,28,1),fpp( 1,28,2)/ 1.05264833d-01, 6.31476720d-04/
      data fpp( 1,29,1),fpp( 1,29,2)/-2.79443142d-02,-4.46577719d-04/
      data fpp( 1,30,1),fpp( 1,30,2)/-7.41572895d-02, 1.22810156d-04/
      data fpp( 1,31,1),fpp( 1,31,2)/-6.74415522d-02, 5.67109353d-06/
      data fpp( 1,32,1),fpp( 1,32,2)/-6.12094757d-02, 4.86594695d-05/
      data fpp( 1,33,1),fpp( 1,33,2)/-5.59358161d-02, 5.39470286d-05/
      data fpp( 1,34,1),fpp( 1,34,2)/-4.69923095d-02, 7.81344162d-05/
      data fpp( 1,35,1),fpp( 1,35,2)/-2.85419899d-02, 9.58933067d-05/
      data fpp( 1,36,1),fpp( 1,36,2)/ 7.72249148d-03, 8.85163569d-05/
      data fpp( 1,37,1),fpp( 1,37,2)/ 7.07957297d-02, 4.29892657d-05/
      data fpp( 1,38,1),fpp( 1,38,2)/ 1.27009555d-01,-6.59894197d-05/
      data fpp( 1,39,1),fpp( 1,39,2)/ 1.26857566d-01,-6.79618736d-05/
      data fpp( 1,40,1),fpp( 1,40,2)/ 7.05539620d-02, 4.28390814d-05/
      data fpp( 1,41,1),fpp( 1,41,2)/ 8.12175005d-03, 7.36835482d-05/
      data fpp( 1,42,1),fpp( 1,42,2)/-2.69536896d-02, 1.27372726d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/-5.07110201d-02, 1.26170551d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-3.50024830d-02, 7.22378978d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-2.97895931d-03, 4.34038577d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 4.61333901d-02,-7.32753284d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 4.61148676d-02,-7.40494435d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-2.90792399d-03, 3.77163179d-05/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-3.43460001d-02, 7.85601718d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-4.91776208d-02, 8.10329948d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-5.39440417d-02, 6.04238490d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-5.27130921d-02, 3.67256091d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-4.82897312d-02, 1.76377146d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/-4.21924174d-02, 5.70953249d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/-3.42979469d-02,-2.57984456d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/-2.28669037d-02, 1.94584576d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/-3.31819555d-03,-1.52835385d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 4.10319701d-02, 5.28763082d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.25731896d-01,-4.60665694d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-4.56843174d-01,-8.04495313d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 2.54536164d-01, 4.32219819d-04/
      data fpp( 2,20,1),fpp( 2,20,2)/ 3.73562030d-01, 1.84084254d-04/
      data fpp( 2,21,1),fpp( 2,21,2)/ 5.19127851d-01,-3.19952834d-04/
      data fpp( 2,22,1),fpp( 2,22,2)/ 5.39046242d-01, 1.62439876d-04/
      data fpp( 2,23,1),fpp( 2,23,2)/ 3.99302832d-01,-6.74564585d-04/
      data fpp( 2,24,1),fpp( 2,24,2)/ 8.70411426d-01, 1.32760765d-02/
      data fpp( 2,25,1),fpp( 2,25,2)/ 2.50425906d+01,-2.23730233d-02/
      data fpp( 2,26,1),fpp( 2,26,2)/ 2.75340991d+00, 6.87125065d-03/
      data fpp( 2,27,1),fpp( 2,27,2)/ 1.01328405d-01, 3.05073267d-03/
      data fpp( 2,28,1),fpp( 2,28,2)/ 9.41128334d-02,-3.45829319d-04/
      data fpp( 2,29,1),fpp( 2,29,2)/-2.65637166d-03,-1.67823389d-04/
      data fpp( 2,30,1),fpp( 2,30,2)/-4.91004210d-02, 3.35087596d-06/
      data fpp( 2,31,1),fpp( 2,31,2)/-5.44493955d-02, 2.14118853d-05/
      data fpp( 2,32,1),fpp( 2,32,2)/-5.62635485d-02, 3.53275827d-05/
      data fpp( 2,33,1),fpp( 2,33,2)/-5.75908677d-02, 4.94377837d-05/
      data fpp( 2,34,1),fpp( 2,34,2)/-5.67353810d-02, 6.98132823d-05/
      data fpp( 2,35,1),fpp( 2,35,2)/-5.07110201d-02, 8.72830872d-05/
      data fpp( 2,36,1),fpp( 2,36,2)/-3.50024830d-02, 8.26543691d-05/
      data fpp( 2,37,1),fpp( 2,37,2)/-2.97895931d-03, 4.06254365d-05/
      data fpp( 2,38,1),fpp( 2,38,2)/ 4.61333901d-02,-7.25781152d-05/
      data fpp( 2,39,1),fpp( 2,39,2)/ 4.61148676d-02,-7.47518725d-05/
      data fpp( 2,40,1),fpp( 2,40,2)/-2.90792399d-03, 4.05364657d-05/
      data fpp( 2,41,1),fpp( 2,41,2)/-3.43460001d-02, 6.79820098d-05/
      data fpp( 2,42,1),fpp( 2,42,2)/-4.91776208d-02, 1.20525495d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-7.92939295d-02, 1.25184740d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-8.63375597d-02, 7.09215194d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-8.14248924d-02, 5.68431821d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.36581154d-02,-9.78822478d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.35670365d-02,-9.87858477d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-8.11222660d-02, 5.11755820d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-8.53427497d-02, 7.87755197d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-7.76858273d-02, 7.70203393d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-6.83196040d-02, 5.55771231d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-5.84254276d-02, 3.29631683d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-4.81284407d-02, 1.49202038d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-3.69815392d-02, 3.09801661d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/-2.34709358d-02,-3.94227020d-06/
      data fpp( 3,14,1),fpp( 3,14,2)/-5.02583696d-03,-5.74293580d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 2.15699344d-02,-6.34998662d-06/
      data fpp( 3,16,1),fpp( 3,16,2)/ 6.04856045d-02,-9.56111773d-06/
      data fpp( 3,17,1),fpp( 3,17,2)/ 8.30783642d-02,-1.41261542d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/-1.05063891d-01,-2.28456712d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/-6.10165745d-02, 2.58876392d-04/
      data fpp( 3,20,1),fpp( 3,20,2)/ 2.19774147d-01, 2.16179143d-04/
      data fpp( 3,21,1),fpp( 3,21,2)/ 5.00050021d-01,-1.71038965d-04/
      data fpp( 3,22,1),fpp( 3,22,2)/ 3.88108154d-01,-1.65047677d-04/
      data fpp( 3,23,1),fpp( 3,23,2)/ 2.67140087d-01, 4.99461994d-04/
      data fpp( 3,24,1),fpp( 3,24,2)/ 1.11606001d+00, 4.40845970d-03/
      data fpp( 3,25,1),fpp( 3,25,2)/ 4.15541773d+00,-7.90662480d-03/
      data fpp( 3,26,1),fpp( 3,26,2)/ 2.40921405d+00, 1.47028149d-03/
      data fpp( 3,27,1),fpp( 3,27,2)/ 2.40873084d-01, 2.01177083d-03/
      data fpp( 3,28,1),fpp( 3,28,2)/ 8.52088332d-02,-2.81174804d-04/
      data fpp( 3,29,1),fpp( 3,29,2)/ 3.25398008d-02,-1.44349611d-04/
      data fpp( 3,30,1),fpp( 3,30,2)/-2.02910264d-02,-2.16927534d-05/
      data fpp( 3,31,1),fpp( 3,31,2)/-4.15058657d-02, 1.43406241d-05/
      data fpp( 3,32,1),fpp( 3,32,2)/-5.28413301d-02, 2.82422569d-05/
      data fpp( 3,33,1),fpp( 3,33,2)/-6.18957129d-02, 4.40623484d-05/
      data fpp( 3,34,1),fpp( 3,34,2)/-7.06761665d-02, 6.37803497d-05/
      data fpp( 3,35,1),fpp( 3,35,2)/-7.92939295d-02, 8.25242530d-05/
      data fpp( 3,36,1),fpp( 3,36,2)/-8.63375597d-02, 8.23486384d-05/
      data fpp( 3,37,1),fpp( 3,37,2)/-8.14248924d-02, 5.37951935d-05/
      data fpp( 3,38,1),fpp( 3,38,2)/-1.36581154d-02,-9.71174124d-05/
      data fpp( 3,39,1),fpp( 3,39,2)/-1.35670365d-02,-9.95563594d-05/
      data fpp( 3,40,1),fpp( 3,40,2)/-8.11222660d-02, 5.42689814d-05/
      data fpp( 3,41,1),fpp( 3,41,2)/-8.53427497d-02, 6.71724339d-05/
      data fpp( 3,42,1),fpp( 3,42,2)/-7.76858273d-02, 1.20339283d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-8.54282620d-02, 1.13039458d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.07537278d-01, 6.85770832d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.43521471d-01, 1.05660209d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.47245929d-01,-1.42665918d-04/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-1.46961721d-01,-1.43908851d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-1.42703012d-01, 1.01187942d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-1.06488001d-01, 7.43170840d-05/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-8.41990700d-02, 7.52257221d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-6.85525422d-02, 5.06360276d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/-5.54351974d-02, 2.95921677d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-4.29215061d-02, 1.24113017d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/-2.97214260d-02, 1.56462535d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.46483100d-02,-4.57580316d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.43025156d-03,-6.91341272d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.54434578d-02,-6.35654596d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/ 4.99506118d-02,-1.21264035d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 6.86046471d-02,-6.47238402d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 3.67873850d-03,-1.02388236d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/-1.93739866d-01, 6.81847829d-05/
      data fpp( 4,20,1),fpp( 4,20,2)/-1.65023616d-01, 2.75893104d-04/
      data fpp( 4,21,1),fpp( 4,21,2)/-3.23162933d-01,-1.80203199d-04/
      data fpp( 4,22,1),fpp( 4,22,2)/-1.09093857d-01,-2.75194456d-04/
      data fpp( 4,23,1),fpp( 4,23,2)/ 1.28736822d-01, 5.74042134d-04/
      data fpp( 4,24,1),fpp( 4,24,2)/ 6.98768532d-01, 1.65032992d-03/
      data fpp( 4,25,1),fpp( 4,25,2)/ 3.17846342d+00,-3.02973381d-03/
      data fpp( 4,26,1),fpp( 4,26,2)/ 1.27231888d+00, 3.22077331d-04/
      data fpp( 4,27,1),fpp( 4,27,2)/ 2.06609260d-01, 1.06685049d-03/
      data fpp( 4,28,1),fpp( 4,28,2)/ 6.84368339d-02,-1.96207289d-04/
      data fpp( 4,29,1),fpp( 4,29,2)/ 3.30571684d-02,-1.26081335d-04/
      data fpp( 4,30,1),fpp( 4,30,2)/-4.61547325d-03,-2.72733714d-05/
      data fpp( 4,31,1),fpp( 4,31,2)/-2.83921417d-02, 7.20482051d-06/
      data fpp( 4,32,1),fpp( 4,32,2)/-4.40861311d-02, 2.24060893d-05/
      data fpp( 4,33,1),fpp( 4,33,2)/-5.69262806d-02, 3.87408222d-05/
      data fpp( 4,34,1),fpp( 4,34,2)/-6.99249530d-02, 5.75306220d-05/
      data fpp( 4,35,1),fpp( 4,35,2)/-8.54282620d-02, 7.99046897d-05/
      data fpp( 4,36,1),fpp( 4,36,2)/-1.07537278d-01, 7.74526192d-05/
      data fpp( 4,37,1),fpp( 4,37,2)/-1.43521471d-01, 1.03292834d-04/
      data fpp( 4,38,1),fpp( 4,38,2)/-1.47245929d-01,-1.42071953d-04/
      data fpp( 4,39,1),fpp( 4,39,2)/-1.46961721d-01,-1.44507057d-04/
      data fpp( 4,40,1),fpp( 4,40,2)/-1.42703012d-01, 1.03589246d-04/
      data fpp( 4,41,1),fpp( 4,41,2)/-1.06488001d-01, 6.53100726d-05/
      data fpp( 4,42,1),fpp( 4,42,2)/-8.41990700d-02, 1.08852464d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-8.06980227d-02, 8.96110036d-05/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.07108327d-01, 6.13359929d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-1.65169224d-01, 1.55269025d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-3.97138170d-01,-1.78526093d-04/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-3.96581078d-01,-1.80199735d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-1.64560687d-01, 1.52603594d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-1.06525245d-01, 6.22833592d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-8.01728926d-02, 7.21249694d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-6.32802270d-02, 4.47667633d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/-4.95937826d-02, 2.61219775d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/-3.65655350d-02, 1.01593267d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-2.29027569d-02, 7.94715641d-07/
      data fpp( 5,13,1),fpp( 5,13,2)/-8.08582429d-03,-4.11618928d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 8.02983071d-03,-5.91795852d-06/
      data fpp( 5,15,1),fpp( 5,15,2)/ 2.55312344d-02,-6.00997663d-06/
      data fpp( 5,16,1),fpp( 5,16,2)/ 4.47519485d-02,-8.66413497d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 5.47480472d-02,-3.22334835d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 2.64039371d-02,-3.95159311d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/-4.79889615d-02, 3.67932078d-05/
      data fpp( 5,20,1),fpp( 5,20,2)/-1.49749682d-01, 2.07169100d-04/
      data fpp( 5,21,1),fpp( 5,21,2)/-2.58298287d-01,-1.36805608d-04/
      data fpp( 5,22,1),fpp( 5,22,2)/-2.78927725d-01,-3.14271227d-04/
      data fpp( 5,23,1),fpp( 5,23,2)/ 6.80762733d-03, 4.18773575d-04/
      data fpp( 5,24,1),fpp( 5,24,2)/ 4.22110861d-01, 7.11828925d-04/
      data fpp( 5,25,1),fpp( 5,25,2)/ 1.08005857d+00,-1.17281528d-03/
      data fpp( 5,26,1),fpp( 5,26,2)/ 7.32295421d-01, 1.00282181d-04/
      data fpp( 5,27,1),fpp( 5,27,2)/ 1.55009876d-01, 5.20298553d-04/
      data fpp( 5,28,1),fpp( 5,28,2)/ 5.08138311d-02,-1.52756392d-04/
      data fpp( 5,29,1),fpp( 5,29,2)/ 3.09515257d-02,-1.00714985d-04/
      data fpp( 5,30,1),fpp( 5,30,2)/ 4.28291942d-03,-2.61856688d-05/
      data fpp( 5,31,1),fpp( 5,31,2)/-1.80255676d-02, 2.92165990d-06/
      data fpp( 5,32,1),fpp( 5,32,2)/-3.50241455d-02, 1.74990292d-05/
      data fpp( 5,33,1),fpp( 5,33,2)/-4.93541647d-02, 3.33962235d-05/
      data fpp( 5,34,1),fpp( 5,34,2)/-6.35940215d-02, 5.03360769d-05/
      data fpp( 5,35,1),fpp( 5,35,2)/-8.06980227d-02, 7.59994688d-05/
      data fpp( 5,36,1),fpp( 5,36,2)/-1.07108327d-01, 6.49820478d-05/
      data fpp( 5,37,1),fpp( 5,37,2)/-1.65169224d-01, 1.54296340d-04/
      data fpp( 5,38,1),fpp( 5,38,2)/-3.97138170d-01,-1.78281407d-04/
      data fpp( 5,39,1),fpp( 5,39,2)/-3.96581078d-01,-1.80447449d-04/
      data fpp( 5,40,1),fpp( 5,40,2)/-1.64560687d-01, 1.53600505d-04/
      data fpp( 5,41,1),fpp( 5,41,2)/-1.06525245d-01, 5.85434270d-05/
      data fpp( 5,42,1),fpp( 5,42,2)/-8.01728926d-02, 8.60877865d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-6.79696472d-02, 8.47900042d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-8.78194150d-02, 5.43439915d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.21991635d-01, 1.23354030d-04/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-1.60771391d-01,-1.55732110d-04/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.61453968d-01,-1.57432184d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.22259242d-01, 1.20608322d-04/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-8.79610179d-02, 5.72828942d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-6.81943596d-02, 6.22261007d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-5.41065497d-02, 3.89927028d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/-4.18396721d-02, 2.21990881d-05/
      data fpp( 6,11,1),fpp( 6,11,2)/-2.96663537d-02, 8.29094492d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/-1.68625465d-02, 4.85132255d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.57339287d-03,-3.09147394d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/ 9.67542560d-03,-4.37923649d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 2.27616047d-02,-4.67558011d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 3.54015944d-02,-6.12644309d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/ 4.42331639d-02,-1.56566475d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 3.46455130d-02,-1.42089667d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/-2.54928814d-03, 3.38585145d-05/
      data fpp( 6,20,1),fpp( 6,20,2)/-2.62526573d-02, 1.14244909d-04/
      data fpp( 6,21,1),fpp( 6,21,2)/ 2.30510815d-02,-8.14641493d-05/
      data fpp( 6,22,1),fpp( 6,22,2)/-1.06402431d-02,-2.58687507d-04/
      data fpp( 6,23,1),fpp( 6,23,2)/-1.94673309d-02, 2.00860338d-04/
      data fpp( 6,24,1),fpp( 6,24,2)/ 2.17883023d-01, 3.67906155d-04/
      data fpp( 6,25,1),fpp( 6,25,2)/ 6.10947285d-01,-4.17182960d-04/
      data fpp( 6,26,1),fpp( 6,26,2)/ 3.80519435d-01, 7.61836830d-05/
      data fpp( 6,27,1),fpp( 6,27,2)/ 9.75462370d-02, 2.32166228d-04/
      data fpp( 6,28,1),fpp( 6,28,2)/ 3.54928417d-02,-1.24354594d-04/
      data fpp( 6,29,1),fpp( 6,29,2)/ 2.67517289d-02,-7.61958535d-05/
      data fpp( 6,30,1),fpp( 6,30,2)/ 8.33379556d-03,-2.23379925d-05/
      data fpp( 6,31,1),fpp( 6,31,2)/-1.03405878d-02, 8.77823454d-07/
      data fpp( 6,32,1),fpp( 6,32,2)/-2.63872871d-02, 1.36546987d-05/
      data fpp( 6,33,1),fpp( 6,33,2)/-4.03020606d-02, 2.78253819d-05/
      data fpp( 6,34,1),fpp( 6,34,2)/-5.34889611d-02, 4.33557739d-05/
      data fpp( 6,35,1),fpp( 6,35,2)/-6.79696472d-02, 6.53615225d-05/
      data fpp( 6,36,1),fpp( 6,36,2)/-8.78194150d-02, 5.95481360d-05/
      data fpp( 6,37,1),fpp( 6,37,2)/-1.21991635d-01, 1.21965933d-04/
      data fpp( 6,38,1),fpp( 6,38,2)/-1.60771391d-01,-1.55383870d-04/
      data fpp( 6,39,1),fpp( 6,39,2)/-1.61453968d-01,-1.57782858d-04/
      data fpp( 6,40,1),fpp( 6,40,2)/-1.22259242d-01, 1.22015885d-04/
      data fpp( 6,41,1),fpp( 6,41,2)/-8.79610179d-02, 5.20033186d-05/
      data fpp( 6,42,1),fpp( 6,42,2)/-6.81943596d-02, 8.19368407d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-4.38958274d-02, 7.22961708d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-4.81763990d-02, 3.89216585d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-4.21217356d-02, 4.78251953d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 4.93416417d-04,-9.12444395d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 3.50610558d-04,-9.22197790d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-4.25920689d-02, 4.47245531d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-4.87931101d-02, 4.42715666d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-4.47535394d-02, 4.39231805d-05/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-3.77080164d-02, 2.96537116d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/-2.92185712d-02, 1.63939732d-05/
      data fpp( 7,11,1),fpp( 7,11,2)/-1.97417975d-02, 6.06239555d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/-9.71300695d-03, 4.14444591d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/ 2.31859089d-04,-1.63017391d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/ 9.41202751d-03,-2.19174894d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.76471616d-02,-2.62283033d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 2.51667196d-02,-3.35492975d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 3.09207555d-02,-5.25145066d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/ 3.12789986d-02,-4.91267601d-07/
      data fpp( 7,19,1),fpp( 7,19,2)/ 2.54836014d-02, 2.48985211d-05/
      data fpp( 7,20,1),fpp( 7,20,2)/ 2.44486453d-02, 4.56711833d-05/
      data fpp( 7,21,1),fpp( 7,21,2)/ 4.41052529d-02,-4.35732544d-05/
      data fpp( 7,22,1),fpp( 7,22,2)/ 3.74559602d-02,-1.55449828d-04/
      data fpp( 7,23,1),fpp( 7,23,2)/ 1.78601802d-03, 1.76354794d-05/
      data fpp( 7,24,1),fpp( 7,24,2)/ 6.43860172d-02, 1.62295911d-04/
      data fpp( 7,25,1),fpp( 7,25,2)/ 1.65223335d-01,-2.13331230d-05/
      data fpp( 7,26,1),fpp( 7,26,2)/ 1.23748269d-01, 7.81725810d-05/
      data fpp( 7,27,1),fpp( 7,27,2)/ 3.30892930d-02, 4.54047989d-05/
      data fpp( 7,28,1),fpp( 7,28,2)/ 1.98679736d-02,-8.53297766d-05/
      data fpp( 7,29,1),fpp( 7,29,2)/ 1.88598866d-02,-4.90136925d-05/
      data fpp( 7,30,1),fpp( 7,30,2)/ 9.72206853d-03,-1.50334534d-05/
      data fpp( 7,31,1),fpp( 7,31,2)/-3.05099560d-03,-1.48493970d-07/
      data fpp( 7,32,1),fpp( 7,32,2)/-1.55496128d-02, 9.51342926d-06/
      data fpp( 7,33,1),fpp( 7,33,2)/-2.67736882d-02, 2.00727769d-05/
      data fpp( 7,34,1),fpp( 7,34,2)/-3.63241152d-02, 3.25894631d-05/
      data fpp( 7,35,1),fpp( 7,35,2)/-4.38958274d-02, 4.61653708d-05/
      data fpp( 7,36,1),fpp( 7,36,2)/-4.81763990d-02, 4.59210536d-05/
      data fpp( 7,37,1),fpp( 7,37,2)/-4.21217356d-02, 4.59584147d-05/
      data fpp( 7,38,1),fpp( 7,38,2)/ 4.93416417d-04,-9.07767124d-05/
      data fpp( 7,39,1),fpp( 7,39,2)/ 3.50610558d-04,-9.26895702d-05/
      data fpp( 7,40,1),fpp( 7,40,2)/-4.25920689d-02, 4.66078458d-05/
      data fpp( 7,41,1),fpp( 7,41,2)/-4.87931101d-02, 3.72081869d-05/
      data fpp( 7,42,1),fpp( 7,42,2)/-4.47535394d-02, 7.02934065d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/-1.67107639d-02, 3.45267888d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.25822743d-02, 1.61334224d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-1.17706533d-03, 6.61752166d-06/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 1.25615020d-02,-3.39695090d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 1.26064268d-02,-3.43182338d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.25943429d-03, 4.98842067d-06/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.29286368d-02, 1.85045511d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-1.74200580d-02, 2.10793750d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.68816178d-02, 1.56299490d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/-1.34887688d-02, 8.84882887d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-8.75123583d-03, 3.40673548d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-3.58964984d-03, 3.54229212d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.53408664d-03,-5.21652330d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/ 6.33145659d-03,-2.89619892d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/ 1.04905202d-02,-2.99868102d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.39855406d-02,-6.52907704d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/ 1.68520842d-02,-5.50501083d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.91198967d-02, 2.56091203d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 2.14220482d-02, 1.12948529d-05/
      data fpp( 8,20,1),fpp( 8,20,2)/ 2.67871295d-02, 1.18516762d-05/
      data fpp( 8,21,1),fpp( 8,21,2)/ 3.60549417d-02,-1.88315576d-05/
      data fpp( 8,22,1),fpp( 8,22,2)/ 3.45770732d-02,-6.24381652d-05/
      data fpp( 8,23,1),fpp( 8,23,2)/ 1.42323409d-02,-3.23618934d-05/
      data fpp( 8,24,1),fpp( 8,24,2)/-1.30506850d-03, 4.45977390d-05/
      data fpp( 8,25,1),fpp( 8,25,2)/ 1.38495693d-03, 6.81349375d-05/
      data fpp( 8,26,1),fpp( 8,26,2)/-2.57492125d-03, 4.72665109d-05/
      data fpp( 8,27,1),fpp( 8,27,2)/ 1.07532027d-03,-1.19629813d-05/
      data fpp( 8,28,1),fpp( 8,28,2)/ 9.47477942d-03,-3.89185859d-05/
      data fpp( 8,29,1),fpp( 8,29,2)/ 1.13733254d-02,-2.16686753d-05/
      data fpp( 8,30,1),fpp( 8,30,2)/ 7.08990337d-03,-7.40271284d-06/
      data fpp( 8,31,1),fpp( 8,31,2)/ 8.87538608d-04, 9.52667185d-09/
      data fpp( 8,32,1),fpp( 8,32,2)/-5.62726665d-03, 5.03660615d-06/
      data fpp( 8,33,1),fpp( 8,33,2)/-1.14925614d-02, 1.06000487d-05/
      data fpp( 8,34,1),fpp( 8,34,2)/-1.56550546d-02, 1.71651989d-05/
      data fpp( 8,35,1),fpp( 8,35,2)/-1.67107639d-02, 2.23011557d-05/
      data fpp( 8,36,1),fpp( 8,36,2)/-1.25822743d-02, 1.94081785d-05/
      data fpp( 8,37,1),fpp( 8,37,2)/-1.17706533d-03, 5.74413053d-06/
      data fpp( 8,38,1),fpp( 8,38,2)/ 1.25615020d-02,-3.37507006d-05/
      data fpp( 8,39,1),fpp( 8,39,2)/ 1.26064268d-02,-3.45379636d-05/
      data fpp( 8,40,1),fpp( 8,40,2)/-1.25943429d-03, 5.86918250d-06/
      data fpp( 8,41,1),fpp( 8,41,2)/-1.29286368d-02, 1.52012336d-05/
      data fpp( 8,42,1),fpp( 8,42,2)/-1.74200580d-02, 3.34118832d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/-6.06911698d-03, 1.41740871d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-3.39850365d-03, 6.20782587d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.45559692d-03, 6.06609443d-07/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 7.24457555d-03,-1.32722636d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 7.26768229d-03,-1.34425138d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.45380607d-03,-6.23900312d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-3.50834259d-03, 6.97207390d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-6.40142864d-03, 9.16409442d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-6.98711249d-03, 7.23954842d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/-5.88075365d-03, 4.22571189d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-3.83965919d-03, 1.68160402d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-1.43559369d-03, 1.77872022d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.04059436d-03,-2.87092112d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/ 3.50374612d-03, 1.64964258d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/ 5.92035746d-03, 4.19106407d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 7.87111807d-03, 4.67077944d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/ 9.48850785d-03, 5.50581810d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.08350144d-02, 1.44059481d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 1.25938056d-02, 3.97103895d-06/
      data fpp( 9,20,1),fpp( 9,20,2)/ 1.58252367d-02, 3.29724940d-06/
      data fpp( 9,21,1),fpp( 9,21,2)/ 2.02373801d-02,-7.72203654d-06/
      data fpp( 9,22,1),fpp( 9,22,2)/ 2.02701472d-02,-2.54175151d-05/
      data fpp( 9,23,1),fpp( 9,23,2)/ 1.14534185d-02,-1.81528364d-05/
      data fpp( 9,24,1),fpp( 9,24,2)/ 2.17456815d-04, 1.37648608d-05/
      data fpp( 9,25,1),fpp( 9,25,2)/-6.71116274d-03, 3.45713934d-05/
      data fpp( 9,26,1),fpp( 9,26,2)/-4.60858377d-03, 2.19195658d-05/
      data fpp( 9,27,1),fpp( 9,27,2)/ 1.42222591d-03,-5.57365659d-06/
      data fpp( 9,28,1),fpp( 9,28,2)/ 5.92650871d-03,-1.58249395d-05/
      data fpp( 9,29,1),fpp( 9,29,2)/ 6.28361160d-03,-9.65458557d-06/
      data fpp( 9,30,1),fpp( 9,30,2)/ 4.12711797d-03,-3.55871825d-06/
      data fpp( 9,31,1),fpp( 9,31,2)/ 1.02964117d-03, 1.54585827d-08/
      data fpp( 9,32,1),fpp( 9,32,2)/-2.13172055d-03, 2.35688392d-06/
      data fpp( 9,33,1),fpp( 9,33,2)/-4.76886626d-03, 5.08300573d-06/
      data fpp( 9,34,1),fpp( 9,34,2)/-6.31966645d-03, 8.04909315d-06/
      data fpp( 9,35,1),fpp( 9,35,2)/-6.06911698d-03, 9.76662168d-06/
      data fpp( 9,36,1),fpp( 9,36,2)/-3.39850365d-03, 7.38842013d-06/
      data fpp( 9,37,1),fpp( 9,37,2)/ 1.45559692d-03, 2.91697808d-07/
      data fpp( 9,38,1),fpp( 9,38,2)/ 7.24457555d-03,-1.31932114d-05/
      data fpp( 9,39,1),fpp( 9,39,2)/ 7.26768229d-03,-1.35222148d-05/
      data fpp( 9,40,1),fpp( 9,40,2)/ 1.45380607d-03, 2.57711685d-07/
      data fpp( 9,41,1),fpp( 9,41,2)/-3.50834259d-03, 5.77136809d-06/
      data fpp( 9,42,1),fpp( 9,42,2)/-6.40142864d-03, 1.36468160d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/-2.09276818d-03, 5.77434962d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/-7.54111039d-04, 2.46230077d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 1.37867765d-03,-6.55526831d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 3.34739578d-03,-5.33809004d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 3.36364406d-03,-5.42995355d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 1.38901000d-03,-3.40098611d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/-7.86792810d-04, 2.71634799d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/-2.21662748d-03, 3.86870664d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/-2.68513224d-03, 3.14482546d-06/
      data fpp(10,10,1),fpp(10,10,2)/-2.36501663d-03, 1.82799151d-06/
      data fpp(10,11,1),fpp(10,11,2)/-1.55332739d-03, 6.97208508d-07/
      data fpp(10,12,1),fpp(10,12,2)/-5.22375410d-04, 5.11744564d-08/
      data fpp(10,13,1),fpp(10,13,2)/ 5.45935935d-04,-2.71906336d-07/
      data fpp(10,14,1),fpp(10,14,2)/ 1.59195893d-03, 1.48450886d-07/
      data fpp(10,15,1),fpp(10,15,2)/ 2.59844990d-03,-1.05897208d-07/
      data fpp(10,16,1),fpp(10,16,2)/ 3.77078715d-03, 8.33137945d-07/
      data fpp(10,17,1),fpp(10,17,2)/ 4.80988441d-03, 6.73345424d-07/
      data fpp(10,18,1),fpp(10,18,2)/ 5.79604548d-03, 6.49480355d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 6.97792952d-03, 1.00673315d-06/
      data fpp(10,20,1),fpp(10,20,2)/ 8.72472370d-03, 7.71587030d-07/
      data fpp(10,21,1),fpp(10,21,2)/ 1.08075379d-02,-2.68308128d-06/
      data fpp(10,22,1),fpp(10,22,2)/ 1.10119382d-02,-1.01625497d-05/
      data fpp(10,23,1),fpp(10,23,2)/ 7.38198519d-03,-7.88453932d-06/
      data fpp(10,24,1),fpp(10,24,2)/ 2.12484124d-03, 4.21870698d-06/
      data fpp(10,25,1),fpp(10,25,2)/-1.27630596d-03, 1.30837114d-05/
      data fpp(10,26,1),fpp(10,26,2)/-9.98743659d-04, 9.86644745d-06/
      data fpp(10,27,1),fpp(10,27,2)/ 1.50617607d-03,-5.59501199d-07/
      data fpp(10,28,1),fpp(10,28,2)/ 3.12398575d-03,-6.13444266d-06/
      data fpp(10,29,1),fpp(10,29,2)/ 3.07542816d-03,-4.54272817d-06/
      data fpp(10,30,1),fpp(10,30,2)/ 2.05282474d-03,-1.72864466d-06/
      data fpp(10,31,1),fpp(10,31,2)/ 6.62696708d-04,-1.46693182d-07/
      data fpp(10,32,1),fpp(10,32,2)/-7.28251164d-04, 9.41417391d-07/
      data fpp(10,33,1),fpp(10,33,2)/-1.84397359d-03, 2.23102362d-06/
      data fpp(10,34,1),fpp(10,34,2)/-2.40067961d-03, 3.52648814d-06/
      data fpp(10,35,1),fpp(10,35,2)/-2.09276818d-03, 4.13502384d-06/
      data fpp(10,36,1),fpp(10,36,2)/-7.54111039d-04, 2.90141649d-06/
      data fpp(10,37,1),fpp(10,37,2)/ 1.37867765d-03,-1.82689814d-07/
      data fpp(10,38,1),fpp(10,38,2)/ 3.34739578d-03,-5.30865724d-06/
      data fpp(10,39,1),fpp(10,39,2)/ 3.36364406d-03,-5.45968338d-06/
      data fpp(10,40,1),fpp(10,40,2)/ 1.38901000d-03,-2.20585244d-07/
      data fpp(10,41,1),fpp(10,41,2)/-7.86792810d-04, 2.26802435d-06/
      data fpp(10,42,1),fpp(10,42,2)/-2.21662748d-03, 5.54248782d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/-2.40936957d-04, 1.01305309d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/-4.76150572d-05, 4.40893815d-07/
      data fpp(11, 3,1),fpp(11, 3,2)/ 2.03168580d-04,-1.06283532d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/ 4.43924874d-04,-8.84380402d-07/
      data fpp(11, 5,1),fpp(11, 5,2)/ 4.46026689d-04,-9.11544619d-07/
      data fpp(11, 6,1),fpp(11, 6,2)/ 2.04266962d-04,-6.39714814d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-5.24502744d-05, 4.59430544d-07/
      data fpp(11, 8,1),fpp(11, 8,2)/-2.55803244d-04, 6.86249307d-07/
      data fpp(11, 9,1),fpp(11, 9,2)/-3.46647019d-04, 5.69572227d-07/
      data fpp(11,10,1),fpp(11,10,2)/-3.34973275d-04, 3.11461784d-07/
      data fpp(11,11,1),fpp(11,11,2)/-2.58388225d-04, 6.85806350d-08/
      data fpp(11,12,1),fpp(11,12,2)/-1.47476925d-04,-1.05784326d-07/
      data fpp(11,13,1),fpp(11,13,2)/-7.10498410d-06,-1.91443333d-07/
      data fpp(11,14,1),fpp(11,14,2)/ 1.76050142d-04,-1.90442345d-07/
      data fpp(11,15,1),fpp(11,15,2)/ 3.76271579d-04,-4.27872896d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 5.60879526d-04, 2.05591502d-07/
      data fpp(11,17,1),fpp(11,17,2)/ 8.23492837d-04, 3.60421280d-07/
      data fpp(11,18,1),fpp(11,18,2)/ 1.12355635d-03, 9.27233720d-08/
      data fpp(11,19,1),fpp(11,19,2)/ 1.45770865d-03,-1.91314768d-07/
      data fpp(11,20,1),fpp(11,20,2)/ 1.80161054d-03, 1.80535699d-07/
      data fpp(11,21,1),fpp(11,21,2)/ 2.06689637d-03, 2.11719716d-08/
      data fpp(11,22,1),fpp(11,22,2)/ 2.16331188d-03,-1.54578376d-06/
      data fpp(11,23,1),fpp(11,23,2)/ 1.80933520d-03,-1.39364136d-06/
      data fpp(11,24,1),fpp(11,24,2)/ 1.23634788d-03, 1.00349192d-07/
      data fpp(11,25,1),fpp(11,25,2)/ 7.03899261d-04, 1.88024459d-06/
      data fpp(11,26,1),fpp(11,26,2)/ 4.72322864d-04, 2.22467245d-06/
      data fpp(11,27,1),fpp(11,27,2)/ 4.69158831d-04, 5.85065614d-07/
      data fpp(11,28,1),fpp(11,28,2)/ 5.70788402d-04,-8.80934905d-07/
      data fpp(11,29,1),fpp(11,29,2)/ 5.41509709d-04,-9.61325998d-07/
      data fpp(11,30,1),fpp(11,30,2)/ 3.68766792d-04,-5.95761107d-07/
      data fpp(11,31,1),fpp(11,31,2)/ 1.11889289d-04,-2.91629577d-07/
      data fpp(11,32,1),fpp(11,32,2)/-1.03386235d-04,-7.72058708d-09/
      data fpp(11,33,1),fpp(11,33,2)/-2.56846103d-04, 3.16511923d-07/
      data fpp(11,34,1),fpp(11,34,2)/-3.12727940d-04, 6.19672894d-07/
      data fpp(11,35,1),fpp(11,35,2)/-2.40936957d-04, 7.50796507d-07/
      data fpp(11,36,1),fpp(11,36,2)/-4.76150572d-05, 5.11141080d-07/
      data fpp(11,37,1),fpp(11,37,2)/ 2.03168580d-04,-2.93608281d-08/
      data fpp(11,38,1),fpp(11,38,2)/ 4.43924874d-04,-8.79697768d-07/
      data fpp(11,39,1),fpp(11,39,2)/ 4.46026689d-04,-9.16226282d-07/
      data fpp(11,40,1),fpp(11,40,2)/ 2.04266962d-04,-4.52467704d-08/
      data fpp(11,41,1),fpp(11,41,2)/-5.24502744d-05, 3.89213362d-07/
      data fpp(11,42,1),fpp(11,42,2)/-2.55803244d-04, 9.48393320d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/-9.70839892d-05, 1.81541541d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/-6.10287324d-05, 8.59169166d-08/
      data fpp(12, 3,1),fpp(12, 3,2)/ 7.04802703d-06, 8.79079081d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/ 6.15047210d-05,-1.33080079d-07/
      data fpp(12, 5,1),fpp(12, 5,2)/ 6.08491892d-05,-1.39155158d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.32215147d-06,-8.90889036d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/-6.66060927d-05, 7.87907210d-08/
      data fpp(12, 8,1),fpp(12, 8,2)/-1.07559546d-04, 1.13746006d-07/
      data fpp(12, 9,1),fpp(12, 9,2)/-1.18679680d-04, 9.02252592d-08/
      data fpp(12,10,1),fpp(12,10,2)/-1.04490266d-04, 4.13529570d-08/
      data fpp(12,11,1),fpp(12,11,2)/-7.65197083d-05, 8.36291408d-09/
      data fpp(12,12,1),fpp(12,12,2)/-3.85168888d-05,-1.16804613d-07/
      data fpp(12,13,1),fpp(12,13,2)/-1.47159988d-05, 1.46855537d-07/
      data fpp(12,14,1),fpp(12,14,2)/-1.79594994d-05,-4.88617537d-07/
      data fpp(12,15,1),fpp(12,15,2)/ 6.46378889d-06, 1.63614607d-07/
      data fpp(12,16,1),fpp(12,16,2)/ 2.70947478d-05,-3.98408948d-08/
      data fpp(12,17,1),fpp(12,17,2)/ 4.29442404d-05, 1.15748968d-07/
      data fpp(12,18,1),fpp(12,18,2)/ 7.65291203d-05, 2.68450186d-08/
      data fpp(12,19,1),fpp(12,19,2)/ 1.09835880d-04,-1.03129044d-07/
      data fpp(12,20,1),fpp(12,20,2)/ 1.44634134d-04, 7.67115656d-09/
      data fpp(12,21,1),fpp(12,21,2)/ 2.07676647d-04, 2.64444418d-07/
      data fpp(12,22,1),fpp(12,22,2)/ 2.35014297d-04,-2.46668831d-07/
      data fpp(12,23,1),fpp(12,23,2)/ 1.57474025d-04,-3.95875849d-07/
      data fpp(12,24,1),fpp(12,24,2)/ 5.29672363d-05,-2.98277741d-08/
      data fpp(12,25,1),fpp(12,25,2)/-3.86910812d-05, 5.27186944d-07/
      data fpp(12,26,1),fpp(12,26,2)/-7.63477972d-05, 5.55079996d-07/
      data fpp(12,27,1),fpp(12,27,2)/-3.18113944d-05, 2.22493068d-07/
      data fpp(12,28,1),fpp(12,28,2)/ 2.22606432d-05,-1.55052272d-07/
      data fpp(12,29,1),fpp(12,29,2)/ 4.51330001d-05, 2.57160158d-08/
      data fpp(12,30,1),fpp(12,30,2)/ 7.90809280d-06,-6.37811794d-07/
      data fpp(12,31,1),fpp(12,31,2)/-8.65386598d-06, 2.35311574d-08/
      data fpp(12,32,1),fpp(12,32,2)/-4.62038967d-05,-1.76312836d-07/
      data fpp(12,33,1),fpp(12,33,2)/-8.48419983d-05, 3.37201854d-08/
      data fpp(12,34,1),fpp(12,34,2)/-1.02608630d-04, 1.07432094d-07/
      data fpp(12,35,1),fpp(12,35,2)/-9.70839892d-05, 1.42551439d-07/
      data fpp(12,36,1),fpp(12,36,2)/-6.10287324d-05, 9.63621485d-08/
      data fpp(12,37,1),fpp(12,37,2)/ 7.04802703d-06, 5.99996780d-09/
      data fpp(12,38,1),fpp(12,38,2)/ 6.15047210d-05,-1.32362019d-07/
      data fpp(12,39,1),fpp(12,39,2)/ 6.08491892d-05,-1.39913926d-07/
      data fpp(12,40,1),fpp(12,40,2)/ 4.32215147d-06,-5.79240447d-09/
      data fpp(12,41,1),fpp(12,41,2)/-6.66060927d-05, 6.70835434d-08/
      data fpp(12,42,1),fpp(12,42,2)/-1.07559546d-04, 1.57458227d-07/
      data fpp(13, 1,1),fpp(13, 1,2)/-1.79554205d-07, 2.57281878d-07/
      data fpp(13, 2,1),fpp(13, 2,2)/-5.50627421d-06, 1.14436241d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/-2.34283710d-05,-1.57026843d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-3.11766000d-05, 2.16711316d-08/
      data fpp(13, 5,1),fpp(13, 5,2)/-3.11609118d-05, 2.25000257d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/-2.29999353d-05,-1.78342415d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/-4.95658475d-06, 1.86869635d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/ 8.80258779d-07,-5.31361239d-08/
      data fpp(13, 9,1),fpp(13, 9,2)/ 2.71254889d-06, 1.96748630d-08/
      data fpp(13,10,1),fpp(13,10,2)/ 1.40743686d-06,-6.75633309d-08/
      data fpp(13,11,1),fpp(13,11,2)/-1.09676275d-06, 1.78578459d-07/
      data fpp(13,12,1),fpp(13,12,2)/-1.26108711d-05,-2.02750506d-07/
      data fpp(13,13,1),fpp(13,13,2)/-1.46995115d-05, 5.04235632d-08/
      data fpp(13,14,1),fpp(13,14,2)/-9.39657272d-06,-2.29437475d-08/
      data fpp(13,15,1),fpp(13,15,2)/-2.18271560d-05, 5.33514237d-08/
      data fpp(13,16,1),fpp(13,16,2)/-3.24740067d-05,-1.96461949d-07/
      data fpp(13,17,1),fpp(13,17,2)/-3.38791398d-05, 1.68496370d-07/
      data fpp(13,18,1),fpp(13,18,2)/-4.31655364d-05,-4.55235348d-08/
      data fpp(13,19,1),fpp(13,19,2)/-5.14119636d-05,-1.64022321d-08/
      data fpp(13,20,1),fpp(13,20,2)/-6.33076717d-05, 1.89132463d-07/
      data fpp(13,21,1),fpp(13,21,2)/-9.01781272d-05,-5.01276191d-08/
      data fpp(13,22,1),fpp(13,22,2)/-9.27488311d-05,-1.30183374d-07/
      data fpp(13,23,1),fpp(13,23,2)/-6.81896741d-05, 1.37355480d-07/
      data fpp(13,24,1),fpp(13,24,2)/-5.37756498d-05, 1.87614530d-08/
      data fpp(13,25,1),fpp(13,25,2)/-4.65763869d-05, 1.11598707d-07/
      data fpp(13,26,1),fpp(13,26,2)/-4.37180405d-05, 6.88437157d-08/
      data fpp(13,27,1),fpp(13,27,2)/-4.80952323d-05,-6.29735742d-08/
      data fpp(13,28,1),fpp(13,28,2)/-4.66261307d-05, 1.59050579d-07/
      data fpp(13,29,1),fpp(13,29,2)/-4.42038549d-05,-2.13228744d-07/
      data fpp(13,30,1),fpp(13,30,2)/-2.13576742d-05, 3.86439679d-09/
      data fpp(13,31,1),fpp(13,31,2)/-1.73830467d-05, 5.77115334d-09/
      data fpp(13,32,1),fpp(13,32,2)/-1.02451926d-05,-2.06949012d-07/
      data fpp(13,33,1),fpp(13,33,2)/ 1.34904643d-06, 1.74024894d-07/
      data fpp(13,34,1),fpp(13,34,2)/ 1.73985927d-06,-3.91505653d-08/
      data fpp(13,35,1),fpp(13,35,2)/-1.79554205d-07,-3.54226284d-08/
      data fpp(13,36,1),fpp(13,36,2)/-5.50627421d-06, 1.92841082d-07/
      data fpp(13,37,1),fpp(13,37,2)/-2.34283710d-05,-1.77941698d-07/
      data fpp(13,38,1),fpp(13,38,2)/-3.11766000d-05, 2.69257130d-08/
      data fpp(13,39,1),fpp(13,39,2)/-3.11609118d-05, 1.71937114d-08/
      data fpp(13,40,1),fpp(13,40,2)/-2.29999353d-05,-1.57013691d-07/
      data fpp(13,41,1),fpp(13,41,2)/-4.95658475d-06, 1.06861054d-07/
      data fpp(13,42,1),fpp(13,42,2)/ 8.80258779d-07, 2.45569473d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.12477940d-05,-4.73828121d-09/
      data fpp(14, 2,1),fpp(14, 2,2)/-9.04617078d-06,-1.52343907d-09/
      data fpp(14, 3,1),fpp(14, 3,2)/-3.18454296d-06, 4.83203654d-09/
      data fpp(14, 4,1),fpp(14, 4,2)/ 3.51678992d-07, 1.95292407d-10/
      data fpp(14, 5,1),fpp(14, 5,2)/ 3.44458240d-07, 2.99810431d-09/
      data fpp(14, 6,1),fpp(14, 6,2)/-3.37241018d-06,-3.79208732d-10/
      data fpp(14, 7,1),fpp(14, 7,2)/-9.41756830d-06,-1.48126809d-09/
      data fpp(14, 8,1),fpp(14, 8,2)/-1.20614895d-05, 3.04279863d-10/
      data fpp(14, 9,1),fpp(14, 9,2)/-1.27705157d-05,-5.73584937d-09/
      data fpp(14,10,1),fpp(14,10,2)/-1.18394809d-05,-1.36088213d-09/
      data fpp(14,11,1),fpp(14,11,2)/-1.01432407d-05,-6.82062432d-09/
      data fpp(14,12,1),fpp(14,12,2)/-6.43962699d-06,-1.35662020d-09/
      data fpp(14,13,1),fpp(14,13,2)/-5.23595529d-06, 2.47103305d-10/
      data fpp(14,14,1),fpp(14,14,2)/-5.80420979d-06, 3.68206433d-10/
      data fpp(14,15,1),fpp(14,15,2)/-4.80516499d-06,-1.71993325d-09/
      data fpp(14,16,1),fpp(14,16,2)/-3.99872113d-06,-5.48847467d-09/
      data fpp(14,17,1),fpp(14,17,2)/-4.92768141d-06,-6.32616987d-09/
      data fpp(14,18,1),fpp(14,18,2)/-3.31697463d-06,-5.20685006d-09/
      data fpp(14,19,1),fpp(14,19,2)/-1.68802539d-06, 3.15356917d-09/
      data fpp(14,20,1),fpp(14,20,2)/ 1.46552998d-07, 1.05925747d-08/
      data fpp(14,21,1),fpp(14,21,2)/ 5.43586185d-06, 1.44761305d-08/
      data fpp(14,22,1),fpp(14,22,2)/ 7.13102747d-06,-7.24678593d-10/
      data fpp(14,23,1),fpp(14,23,2)/-6.53287892d-08,-6.04189241d-10/
      data fpp(14,24,1),fpp(14,24,2)/-7.96463702d-06, 9.14143470d-09/
      data fpp(14,25,1),fpp(14,25,2)/-1.44033712d-05, 1.80384479d-08/
      data fpp(14,26,1),fpp(14,26,2)/-1.60800409d-05, 8.70477147d-09/
      data fpp(14,27,1),fpp(14,27,2)/-1.16076764d-05, 1.14246129d-09/
      data fpp(14,28,1),fpp(14,28,2)/-6.45612053d-06,-7.27462083d-09/
      data fpp(14,29,1),fpp(14,29,2)/-2.11758035d-06,-8.04397864d-09/
      data fpp(14,30,1),fpp(14,30,2)/-3.32739599d-06,-8.54946763d-09/
      data fpp(14,31,1),fpp(14,31,2)/-3.86394699d-06,-5.75815320d-09/
      data fpp(14,32,1),fpp(14,32,2)/-6.11533284d-06,-1.04179208d-08/
      data fpp(14,33,1),fpp(14,33,2)/-1.01041874d-05,-5.70162051d-10/
      data fpp(14,34,1),fpp(14,34,2)/-1.14508072d-05,-5.30143114d-09/
      data fpp(14,35,1),fpp(14,35,2)/-1.12477940d-05, 3.77589116d-09/
      data fpp(14,36,1),fpp(14,36,2)/-9.04617078d-06,-3.80213270d-09/
      data fpp(14,37,1),fpp(14,37,2)/-3.18454296d-06, 5.43264179d-09/
      data fpp(14,38,1),fpp(14,38,2)/ 3.51678992d-07, 7.15655343d-11/
      data fpp(14,39,1),fpp(14,39,2)/ 3.44458240d-07, 3.06898401d-09/
      data fpp(14,40,1),fpp(14,40,2)/-3.37241018d-06,-5.57033510d-10/
      data fpp(14,41,1),fpp(14,41,2)/-9.41756830d-06,-8.40848236d-10/
      data fpp(14,42,1),fpp(14,42,2)/-1.20614895d-05,-2.07957584d-09/
      data fpp(15, 1,1),fpp(15, 1,2)/ 1.18456827d-05,-2.95305337d-11/
      data fpp(15, 2,1),fpp(15, 2,2)/ 7.35808538d-06, 5.90590814d-11/
      data fpp(15, 3,1),fpp(15, 3,2)/ 8.53405719d-06,-2.06708444d-10/
      data fpp(15, 4,1),fpp(15, 4,2)/ 1.97880337d-06, 7.67773860d-10/
      data fpp(15, 5,1),fpp(15, 5,2)/ 2.19241374d-06, 8.00031719d-10/
      data fpp(15, 6,1),fpp(15, 6,2)/ 8.80799080d-06,-3.35736252d-10/
      data fpp(15, 7,1),fpp(15, 7,2)/ 8.19628414d-06, 5.42913108d-10/
      data fpp(15, 8,1),fpp(15, 8,2)/ 1.31514590d-05,-1.83591686d-09/
      data fpp(15, 9,1),fpp(15, 9,2)/ 1.43266864d-05, 8.00756707d-10/
      data fpp(15,10,1),fpp(15,10,2)/ 1.20933119d-05,-1.36711133d-09/
      data fpp(15,11,1),fpp(15,11,2)/ 7.34197750d-06,-1.33231107d-09/
      data fpp(15,12,1),fpp(15,12,2)/ 1.00726706d-05, 6.96354177d-10/
      data fpp(15,13,1),fpp(15,13,2)/ 5.48619192d-06,-1.45310684d-09/
      data fpp(15,14,1),fpp(15,14,2)/ 1.41174774d-06,-8.83927526d-10/
      data fpp(15,15,1),fpp(15,15,2)/ 4.06865392d-06,-1.01118572d-09/
      data fpp(15,16,1),fpp(15,16,2)/ 6.34507484d-06,-1.07133215d-09/
      data fpp(15,17,1),fpp(15,17,2)/-3.88302162d-07,-7.03489159d-10/
      data fpp(15,18,1),fpp(15,18,2)/-2.19972698d-06,-2.11471423d-09/
      data fpp(15,19,1),fpp(15,19,2)/-4.66527302d-06, 3.16234584d-09/
      data fpp(15,20,1),fpp(15,20,2)/-4.45006221d-06, 1.46533051d-09/
      data fpp(15,21,1),fpp(15,21,2)/ 4.19921192d-06, 2.97633023d-09/
      data fpp(15,22,1),fpp(15,22,2)/-2.70870888d-07, 8.38344857d-10/
      data fpp(15,23,1),fpp(15,23,2)/-5.88912132d-06,-1.98272732d-09/
      data fpp(15,24,1),fpp(15,24,2)/ 4.66588992d-06, 1.09256331d-09/
      data fpp(15,25,1),fpp(15,25,2)/ 1.91288284d-05, 3.61247259d-09/
      data fpp(15,26,1),fpp(15,26,2)/ 2.72121633d-05, 2.45754532d-09/
      data fpp(15,27,1),fpp(15,27,2)/ 2.37224096d-05, 4.55734178d-09/
      data fpp(15,28,1),fpp(15,28,2)/ 1.37044889d-05,-2.68691574d-09/
      data fpp(15,29,1),fpp(15,29,2)/ 1.26945045d-05, 1.90317569d-10/
      data fpp(15,30,1),fpp(15,30,2)/ 6.32976941d-06,-4.07435802d-09/
      data fpp(15,31,1),fpp(15,31,2)/ 8.23518776d-06,-1.89288613d-09/
      data fpp(15,32,1),fpp(15,32,2)/ 1.17812378d-05,-3.54097328d-10/
      data fpp(15,33,1),fpp(15,33,2)/ 8.61459369d-06,-2.69072586d-09/
      data fpp(15,34,1),fpp(15,34,2)/ 1.19761179d-05,-8.82996512d-10/
      data fpp(15,35,1),fpp(15,35,2)/ 1.18456827d-05, 2.22717646d-10/
      data fpp(15,36,1),fpp(15,36,2)/ 7.35808538d-06,-7.87292481d-12/
      data fpp(15,37,1),fpp(15,37,2)/ 8.53405719d-06,-1.91225278d-10/
      data fpp(15,38,1),fpp(15,38,2)/ 1.97880337d-06, 7.72775900d-10/
      data fpp(15,39,1),fpp(15,39,2)/ 2.19241374d-06, 7.77286386d-10/
      data fpp(15,40,1),fpp(15,40,2)/ 8.80799080d-06,-2.09269014d-10/
      data fpp(15,41,1),fpp(15,41,2)/ 8.19628414d-06, 5.97902353d-11/
      data fpp(15,42,1),fpp(15,42,2)/ 1.31514590d-05,-2.98955748d-11/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 6.05166096d-04, 2.67293903d-04/
      data fpppp( 1, 3),fpppp( 1, 4)/-6.58163039d-05,-4.15593460d-04/
      data fpppp( 1, 5),fpppp( 1, 6)/-4.09006057d-04,-8.84333304d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/ 3.95022914d-04, 1.49748010d-04/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.11179053d-04, 5.21896611d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 8.41496453d-06, 5.27068917d-05/
      data fpppp( 1,13),fpppp( 1,14)/-1.27282229d-04, 6.41774844d-04/
      data fpppp( 1,15),fpppp( 1,16)/-1.88189710d-03, 9.07726981d-03/
      data fpppp( 1,17),fpppp( 1,18)/-2.97568500d-02, 4.81678299d-02/
      data fpppp( 1,19),fpppp( 1,20)/-3.35213018d-02, 1.50253066d-02/
      data fpppp( 1,21),fpppp( 1,22)/-2.06336233d-02, 4.94644159d-02/
      data fpppp( 1,23),fpppp( 1,24)/-2.71488096d-01, 1.06868361d+00/
      data fpppp( 1,25),fpppp( 1,26)/-1.65448455d+00, 9.68614971d-01/
      data fpppp( 1,27),fpppp( 1,28)/-2.03711910d-01, 5.06208402d-02/
      data fpppp( 1,29),fpppp( 1,30)/-1.23932921d-02, 4.17209836d-03/
      data fpppp( 1,31),fpppp( 1,32)/-1.11937861d-03, 2.76396437d-04/
      data fpppp( 1,33),fpppp( 1,34)/-4.37121518d-05, 1.18642993d-04/
      data fpppp( 1,35),fpppp( 1,36)/ 1.39548956d-04, 3.92010892d-04/
      data fpppp( 1,37),fpppp( 1,38)/-9.90671193d-05,-4.07307188d-04/
      data fpppp( 1,39),fpppp( 1,40)/-4.17239466d-04,-5.56054188d-05/
      data fpppp( 1,41),fpppp( 1,42)/ 2.71944678d-04, 6.09233046d-04/
      data fpppp( 2, 1),fpppp( 2, 2)/ 1.78047179d-04, 1.17878449d-04/
      data fpppp( 2, 3),fpppp( 2, 4)/ 3.29338214d-04,-4.09901757d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/-4.08612157d-04, 3.30664637d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.41036535d-04, 1.01576548d-04/
      data fpppp( 2, 9),fpppp( 2,10)/ 5.65692579d-05, 3.19886532d-05/
      data fpppp( 2,11),fpppp( 2,12)/ 7.02080608d-06, 4.03653012d-05/
      data fpppp( 2,13),fpppp( 2,14)/-6.06526142d-05, 4.14439520d-04/
      data fpppp( 2,15),fpppp( 2,16)/-1.11004557d-03, 5.51383020d-03/
      data fpppp( 2,17),fpppp( 2,18)/-1.85242896d-02, 2.85468285d-02/
      data fpppp( 2,19),fpppp( 2,20)/-1.80257600d-02, 8.01500319d-03/
      data fpppp( 2,21),fpppp( 2,22)/-1.24418553d-02, 2.92498656d-02/
      data fpppp( 2,23),fpppp( 2,24)/-1.59597639d-01, 6.45791812d-01/
      data fpppp( 2,25),fpppp( 2,26)/-1.00150537d+00, 5.72548072d-01/
      data fpppp( 2,27),fpppp( 2,28)/-1.10460964d-01, 2.79877388d-02/
      data fpppp( 2,29),fpppp( 2,30)/-6.86320964d-03, 2.48460910d-03/
      data fpppp( 2,31),fpppp( 2,32)/-6.09522285d-04, 1.65569322d-04/
      data fpppp( 2,33),fpppp( 2,34)/-2.35449728d-05, 5.95789251d-05/
      data fpppp( 2,35),fpppp( 2,36)/ 9.53617193d-05, 1.40024779d-04/
      data fpppp( 2,37),fpppp( 2,38)/ 3.23438354d-04,-4.08448647d-04/
      data fpppp( 2,39),fpppp( 2,40)/-4.10021557d-04, 3.36214815d-04/
      data fpppp( 2,41),fpppp( 2,42)/ 1.20245223d-04, 1.79191620d-04/
      data fpppp( 3, 1),fpppp( 3, 2)/-3.20867835d-04,-1.75316558d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 1.10837230d-03,-6.44710966d-04/
      data fpppp( 3, 5),fpppp( 3, 6)/-6.51690382d-04, 1.14351809d-03/
      data fpppp( 3, 7),fpppp( 3, 8)/-1.22297227d-04, 5.83151770d-05/
      data fpppp( 3, 9),fpppp( 3,10)/-8.40542479d-06, 6.98370903d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 4.63922075d-06, 2.54542827d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 3.53657598d-05, 1.29152404d-04/
      data fpppp( 3,15),fpppp( 3,16)/-6.29350193d-05, 8.61781595d-04/
      data fpppp( 3,17),fpppp( 3,18)/-4.36356598d-03, 3.94838144d-03/
      data fpppp( 3,19),fpppp( 3,20)/ 2.50141456d-03, 2.50564608d-04/
      data fpppp( 3,21),fpppp( 3,22)/-3.53456384d-03, 3.91004986d-04/
      data fpppp( 3,23),fpppp( 3,24)/ 8.23269709d-04, 5.45091957d-02/
      data fpppp( 3,25),fpppp( 3,26)/-8.74337845d-02, 8.09225828d-03/
      data fpppp( 3,27),fpppp( 3,28)/ 2.97365142d-02,-6.27771227d-03/
      data fpppp( 3,29),fpppp( 3,30)/ 1.55404792d-03, 5.18128898d-05/
      data fpppp( 3,31),fpppp( 3,32)/ 1.35659799d-04,-1.68959206d-06/
      data fpppp( 3,33),fpppp( 3,34)/ 7.96346194d-06,-1.37285006d-05/
      data fpppp( 3,35),fpppp( 3,36)/ 5.67119770d-05,-1.18671440d-04/
      data fpppp( 3,37),fpppp( 3,38)/ 1.13535163d-03,-6.51488478d-04/
      data fpppp( 3,39),fpppp( 3,40)/-6.44847507d-04, 1.11601586d-03/
      data fpppp( 3,41),fpppp( 3,42)/-1.91311992d-05,-3.26846710d-04/
      data fpppp( 4, 1),fpppp( 4, 2)/-5.68496508d-04,-2.01960487d-04/
      data fpppp( 4, 3),fpppp( 4, 4)/ 5.43827891d-04,-3.77669790d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/-4.26161714d-05, 5.78227345d-04/
      data fpppp( 4, 7),fpppp( 4, 8)/-3.52915147d-04,-2.13151481d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-3.71030016d-05,-1.20745833d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 5.71362952d-06, 1.95362623d-05/
      data fpppp( 4,13),fpppp( 4,14)/ 2.85234775d-05, 4.66965600d-05/
      data fpppp( 4,15),fpppp( 4,16)/ 2.07689635d-05, 1.98644499d-05/
      data fpppp( 4,17),fpppp( 4,18)/-4.51413878d-04,-3.22900558d-03/
      data fpppp( 4,19),fpppp( 4,20)/ 5.41787443d-03,-4.87440090d-03/
      data fpppp( 4,21),fpppp( 4,22)/ 2.86839517d-03, 1.78723059d-03/
      data fpppp( 4,23),fpppp( 4,24)/-8.61240544d-03, 5.25944531d-02/
      data fpppp( 4,25),fpppp( 4,26)/-8.71856160d-02, 3.29976449d-02/
      data fpppp( 4,27),fpppp( 4,28)/ 5.62113162d-03, 1.70060378d-04/
      data fpppp( 4,29),fpppp( 4,30)/-1.33807495d-04, 2.27591040d-04/
      data fpppp( 4,31),fpppp( 4,32)/ 5.72017263d-05, 2.85627926d-05/
      data fpppp( 4,33),fpppp( 4,34)/-2.22498631d-07,-3.71841742d-05/
      data fpppp( 4,35),fpppp( 4,36)/-1.31899622d-06,-3.53882292d-04/
      data fpppp( 4,37),fpppp( 4,38)/ 5.84337598d-04,-4.78840043d-05/
      data fpppp( 4,39),fpppp( 4,40)/-3.25199490d-05, 5.37884061d-04/
      data fpppp( 4,41),fpppp( 4,42)/-2.01638236d-04,-5.66895875d-04/
      data fpppp( 5, 1),fpppp( 5, 2)/ 9.02377935d-04, 9.08754717d-05/
      data fpppp( 5, 3),fpppp( 5, 4)/-3.16491540d-03, 2.13430316d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/ 2.14697300d-03,-3.24593361d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/ 3.97664432d-04,-2.45709444d-04/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.75921215d-05,-1.70323110d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 1.10453122d-05, 1.09228991d-05/
      data fpppp( 5,13),fpppp( 5,14)/ 1.45123552d-05, 8.95102509d-06/
      data fpppp( 5,15),fpppp( 5,16)/ 3.28284654d-05,-3.71062634d-05/
      data fpppp( 5,17),fpppp( 5,18)/-4.37880329d-04,-5.11784955d-04/
      data fpppp( 5,19),fpppp( 5,20)/-2.77907159d-04,-1.86557022d-05/
      data fpppp( 5,21),fpppp( 5,22)/-5.47431480d-05, 3.12057389d-03/
      data fpppp( 5,23),fpppp( 5,24)/-8.50952770d-04, 8.05731010d-03/
      data fpppp( 5,25),fpppp( 5,26)/-1.68196190d-02,-1.12148603d-03/
      data fpppp( 5,27),fpppp( 5,28)/ 7.53421954d-03,-6.30022110d-04/
      data fpppp( 5,29),fpppp( 5,30)/ 4.58932552d-05, 3.80710374d-05/
      data fpppp( 5,31),fpppp( 5,32)/ 6.34297476d-05, 2.68045263d-05/
      data fpppp( 5,33),fpppp( 5,34)/-1.05343380d-05, 2.07425748d-05/
      data fpppp( 5,35),fpppp( 5,36)/-2.44284630d-04, 3.98017784d-04/
      data fpppp( 5,37),fpppp( 5,38)/-3.24682209d-03, 2.15478758d-03/
      data fpppp( 5,39),fpppp( 5,40)/ 2.12647307d-03,-3.16390289d-03/
      data fpppp( 5,41),fpppp( 5,42)/ 9.00415114d-05, 9.02751525d-04/
      data fpppp( 6, 1),fpppp( 6, 2)/-2.49896905d-04,-1.24178051d-04/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.12738005d-04, 2.98677870d-04/
      data fpppp( 6, 5),fpppp( 6, 6)/ 3.13489430d-04,-1.06131482d-04/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.82753593d-04,-3.47480977d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-1.89849149d-05, 1.43181380d-06/
      data fpppp( 6,11),fpppp( 6,12)/ 7.64410668d-06, 5.82109018d-06/
      data fpppp( 6,13),fpppp( 6,14)/-1.80767934d-06,-1.01048486d-06/
      data fpppp( 6,15),fpppp( 6,16)/-3.90874635d-06,-1.01258869d-05/
      data fpppp( 6,17),fpppp( 6,18)/-1.84092923d-04,-3.58655646d-04/
      data fpppp( 6,19),fpppp( 6,20)/-3.77135085d-05, 1.31899560d-03/
      data fpppp( 6,21),fpppp( 6,22)/-8.57842423d-04,-7.04525654d-05/
      data fpppp( 6,23),fpppp( 6,24)/ 2.61951471d-03, 4.36304020d-03/
      data fpppp( 6,25),fpppp( 6,26)/-1.07288410d-02, 1.14279703d-03/
      data fpppp( 6,27),fpppp( 6,28)/ 3.00493195d-03, 9.26633374d-05/
      data fpppp( 6,29),fpppp( 6,30)/-1.76848357d-04, 3.41208622d-05/
      data fpppp( 6,31),fpppp( 6,32)/ 2.49779066d-05, 2.36285569d-05/
      data fpppp( 6,33),fpppp( 6,34)/ 8.42340891d-06,-1.36498123d-05/
      data fpppp( 6,35),fpppp( 6,36)/-3.14512928d-05,-1.82689916d-04/
      data fpppp( 6,37),fpppp( 6,38)/-9.71361551d-05, 2.94782335d-04/
      data fpppp( 6,39),fpppp( 6,40)/ 3.17375110d-04,-1.21654491d-04/
      data fpppp( 6,41),fpppp( 6,42)/-1.24547237d-04,-2.52050513d-04/
      data fpppp( 7, 1),fpppp( 7, 2)/-1.19417910d-04, 2.43877910d-05/
      data fpppp( 7, 3),fpppp( 7, 4)/ 6.41980839d-04,-3.98681829d-04/
      data fpppp( 7, 5),fpppp( 7, 6)/-4.05541581d-04, 6.58336550d-04/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.33063253d-05, 4.93254664d-05/
      data fpppp( 7, 9),fpppp( 7,10)/ 6.36160257d-06, 1.18634476d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 5.42432436d-06,-4.39739065d-07/
      data fpppp( 7,13),fpppp( 7,14)/-8.70083750d-06,-1.06387675d-05/
      data fpppp( 7,15),fpppp( 7,16)/-5.44615582d-06,-1.05111675d-05/
      data fpppp( 7,17),fpppp( 7,18)/-5.84405072d-05,-7.94743677d-05/
      data fpppp( 7,19),fpppp( 7,20)/ 7.11956307d-06, 3.36622572d-04/
      data fpppp( 7,21),fpppp( 7,22)/-1.12116021d-04,-5.21400843d-04/
      data fpppp( 7,23),fpppp( 7,24)/ 1.41191935d-03, 7.69919938d-04/
      data fpppp( 7,25),fpppp( 7,26)/-2.19735998d-03,-5.19223080d-04/
      data fpppp( 7,27),fpppp( 7,28)/ 1.32321773d-03,-1.27388434d-04/
      data fpppp( 7,29),fpppp( 7,30)/-8.08700452d-05,-3.69152532d-05/
      data fpppp( 7,31),fpppp( 7,32)/ 1.04162965d-05, 1.17168803d-05/
      data fpppp( 7,33),fpppp( 7,34)/ 1.91886956d-05, 1.19472372d-05/
      data fpppp( 7,35),fpppp( 7,36)/ 5.17452454d-05,-2.14597791d-05/
      data fpppp( 7,37),fpppp( 7,38)/ 6.54207965d-04,-4.01742760d-04/
      data fpppp( 7,39),fpppp( 7,40)/-4.02472351d-04, 6.46043035d-04/
      data fpppp( 7,41),fpppp( 7,42)/ 2.27985040d-05,-1.22800336d-04/
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.25363627d-04, 6.64415528d-05/
      data fpppp( 8, 3),fpppp( 8, 4)/ 4.54733296d-05,-1.08333372d-04/
      data fpppp( 8, 5),fpppp( 8, 6)/-1.09219697d-04, 3.86855205d-05/
      data fpppp( 8, 7),fpppp( 8, 8)/ 8.62771288d-05, 4.68728461d-05/
      data fpppp( 8, 9),fpppp( 8,10)/ 2.80231682d-05, 1.22990093d-05/
      data fpppp( 8,11),fpppp( 8,12)/ 3.46183125d-06,-7.03151593d-07/
      data fpppp( 8,13),fpppp( 8,14)/-2.92019558d-06,-7.19805775d-06/
      data fpppp( 8,15),fpppp( 8,16)/-6.58595135d-06,-6.30073614d-06/
      data fpppp( 8,17),fpppp( 8,18)/-5.91970863d-06,-5.94429046d-06/
      data fpppp( 8,19),fpppp( 8,20)/ 3.17572064d-05, 6.26912505d-05/
      data fpppp( 8,21),fpppp( 8,22)/-4.83583497d-05,-1.86472972d-04/
      data fpppp( 8,23),fpppp( 8,24)/ 3.92066535d-05, 3.18085733d-04/
      data fpppp( 8,25),fpppp( 8,26)/-2.17903496d-04, 1.54534036d-04/
      data fpppp( 8,27),fpppp( 8,28)/ 5.63745334d-05,-9.50791131d-05/
      data fpppp( 8,29),fpppp( 8,30)/-6.61128685d-05,-1.13874983d-05/
      data fpppp( 8,31),fpppp( 8,32)/-3.47370044d-06, 6.53587053d-06/
      data fpppp( 8,33),fpppp( 8,34)/ 1.63008499d-05, 3.04288222d-05/
      data fpppp( 8,35),fpppp( 8,36)/ 4.83908939d-05, 8.70595345d-05/
      data fpppp( 8,37),fpppp( 8,38)/ 3.99741361d-05,-1.06954580d-04/
      data fpppp( 8,39),fpppp( 8,40)/-1.10606478d-04, 4.42486182d-05/
      data fpppp( 8,41),fpppp( 8,42)/ 6.54115186d-05, 1.24772189d-04/
      data fpppp( 9, 1),fpppp( 9, 2)/ 3.59637265d-05, 1.85493482d-05/
      data fpppp( 9, 3),fpppp( 9, 4)/ 2.08481156d-05,-4.58491269d-05/
      data fpppp( 9, 5),fpppp( 9, 6)/-4.61994351d-05, 1.93690895d-05/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.98267301d-05, 2.54677476d-05/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.67464106d-05, 9.06917212d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 3.06103728d-06, 4.64941974d-07/
      data fpppp( 9,13),fpppp( 9,14)/-5.93452959d-07, 1.12669298d-06/
      data fpppp( 9,15),fpppp( 9,16)/-6.70574397d-06,-2.25476158d-06/
      data fpppp( 9,17),fpppp( 9,18)/-4.27745876d-06, 3.11160493d-06/
      data fpppp( 9,19),fpppp( 9,20)/ 1.65681103d-05, 1.89743557d-05/
      data fpppp( 9,21),fpppp( 9,22)/-2.16227972d-05,-7.64915824d-05/
      data fpppp( 9,23),fpppp( 9,24)/-2.77916443d-05, 4.25041805d-05/
      data fpppp( 9,25),fpppp( 9,26)/ 1.16215449d-04, 3.45059340d-05/
      data fpppp( 9,27),fpppp( 9,28)/-1.85453423d-05,-5.19161782d-05/
      data fpppp( 9,29),fpppp( 9,30)/-2.26207390d-05,-8.41665676d-06/
      data fpppp( 9,31),fpppp( 9,32)/-1.71624297d-07, 5.27005892d-06/
      data fpppp( 9,33),fpppp( 9,34)/ 1.05443490d-05, 1.77332759d-05/
      data fpppp( 9,35),fpppp( 9,36)/ 2.66035275d-05, 2.10564454d-05/
      data fpppp( 9,37),fpppp( 9,38)/ 2.01799258d-05,-4.56834648d-05/
      data fpppp( 9,39),fpppp( 9,40)/-4.63623263d-05, 2.00151123d-05/
      data fpppp( 9,41),fpppp( 9,42)/ 1.74055300d-05, 3.45065251d-05/
      data fpppp(10, 1),fpppp(10, 2)/ 1.76526102d-05, 7.68276807d-06/
      data fpppp(10, 3),fpppp(10, 4)/-7.35789735d-07,-1.45838429d-05/
      data fpppp(10, 5),fpppp(10, 6)/-1.46983964d-05,-1.60742762d-06/
      data fpppp(10, 7),fpppp(10, 8)/ 9.05798152d-06, 1.01335899d-05/
      data fpppp(10, 9),fpppp(10,10)/ 8.08745300d-06, 4.83382063d-06/
      data fpppp(10,11),fpppp(10,12)/ 2.07168230d-06, 3.52146730d-08/
      data fpppp(10,13),fpppp(10,14)/ 2.90208110d-08,-1.48859878d-06/
      data fpppp(10,15),fpppp(10,16)/ 3.55345227d-06,-2.77443307d-06/
      data fpppp(10,17),fpppp(10,18)/-4.50119054d-07, 1.39873717d-06/
      data fpppp(10,19),fpppp(10,20)/ 6.59854952d-06, 6.10167261d-06/
      data fpppp(10,21),fpppp(10,22)/-1.08440414d-05,-2.99371320d-05/
      data fpppp(10,23),fpppp(10,24)/-2.26183143d-05, 2.27789318d-05/
      data fpppp(10,25),fpppp(10,26)/ 4.28623921d-05, 2.64940702d-05/
      data fpppp(10,27),fpppp(10,28)/-1.51972273d-05,-1.89317641d-05/
      data fpppp(10,29),fpppp(10,30)/-9.05775193d-06,-3.27997852d-06/
      data fpppp(10,31),fpppp(10,32)/ 1.26189410d-07, 2.72603049d-06/
      data fpppp(10,33),fpppp(10,34)/ 5.48321545d-06, 8.88209195d-06/
      data fpppp(10,35),fpppp(10,36)/ 1.08654637d-05, 9.50079646d-06/
      data fpppp(10,37),fpppp(10,38)/-1.22075674d-06,-1.44620032d-05/
      data fpppp(10,39),fpppp(10,40)/-1.48214318d-05,-1.11289420d-06/
      data fpppp(10,41),fpppp(10,42)/ 7.20288330d-06, 1.70594494d-05/
      data fpppp(11, 1),fpppp(11, 2)/ 1.28937624d-06, 4.94881278d-07/
      data fpppp(11, 3),fpppp(11, 4)/ 1.78802915d-07,-1.81173353d-06/
      data fpppp(11, 5),fpppp(11, 6)/-1.84536247d-06, 1.27003826d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 4.39896595d-07, 1.31526583d-06/
      data fpppp(11, 9),fpppp(11,10)/ 1.04959175d-06, 6.37418309d-07/
      data fpppp(11,11),fpppp(11,12)/ 2.95413391d-07, 2.40503088d-07/
      data fpppp(11,13),fpppp(11,14)/ 5.10212757d-07, 2.85636960d-07/
      data fpppp(11,15),fpppp(11,16)/-6.28781916d-07, 1.29268134d-06/
      data fpppp(11,17),fpppp(11,18)/ 1.38378323d-07, 4.00817595d-07/
      data fpppp(11,19),fpppp(11,20)/ 3.03678315d-07,-1.03055523d-06/
      data fpppp(11,21),fpppp(11,22)/-8.98420903d-07,-3.30180204d-06/
      data fpppp(11,23),fpppp(11,24)/-2.52341220d-06, 2.54812989d-07/
      data fpppp(11,25),fpppp(11,26)/ 3.93648195d-06, 2.05159264d-06/
      data fpppp(11,27),fpppp(11,28)/ 1.56188933d-06,-2.01153368d-06/
      data fpppp(11,29),fpppp(11,30)/-1.37025048d-06,-1.11531783d-06/
      data fpppp(11,31),fpppp(11,32)/ 7.83446735d-07, 4.77649568d-07/
      data fpppp(11,33),fpppp(11,34)/ 1.01489436d-06, 1.31745489d-06/
      data fpppp(11,35),fpppp(11,36)/ 1.37565525d-06, 4.71779072d-07/
      data fpppp(11,37),fpppp(11,38)/ 1.84932729d-07,-1.81315058d-06/
      data fpppp(11,39),fpppp(11,40)/-1.84417622d-06, 1.22720423d-07/
      data fpppp(11,41),fpppp(11,42)/ 4.55843954d-07, 1.25575979d-06/
      data fpppp(12, 1),fpppp(12, 2)/ 7.68594651d-07, 3.36287138d-07/
      data fpppp(12, 3),fpppp(12, 4)/-1.92453043d-07,-3.83678892d-07/
      data fpppp(12, 5),fpppp(12, 6)/-3.96270601d-07,-2.26974914d-07/
      data fpppp(12, 7),fpppp(12, 8)/ 4.40097868d-07, 2.65070911d-07/
      data fpppp(12, 9),fpppp(12,10)/ 2.89617608d-07, 9.50315053d-08/
      data fpppp(12,11),fpppp(12,12)/ 1.57125055d-07,-1.21596049d-07/
      data fpppp(12,13),fpppp(12,14)/-5.22856639d-07, 5.90359175d-07/
      data fpppp(12,15),fpppp(12,16)/-1.78572733d-07,-1.03607998d-07/
      data fpppp(12,17),fpppp(12,18)/ 3.06116749d-07,-5.67357677d-08/
      data fpppp(12,19),fpppp(12,20)/-9.58609026d-08, 5.29669060d-07/
      data fpppp(12,21),fpppp(12,22)/-3.28159802d-07,-7.61565763d-07/
      data fpppp(12,23),fpppp(12,24)/-2.46831618d-07, 1.30901202d-07/
      data fpppp(12,25),fpppp(12,26)/ 4.94135086d-07, 1.13265454d-06/
      data fpppp(12,27),fpppp(12,28)/-9.31661383d-08,-1.87851902d-07/
      data fpppp(12,29),fpppp(12,30)/-1.02740709d-06, 6.91644395d-07/
      data fpppp(12,31),fpppp(12,32)/-4.99393581d-07, 4.66456146d-08/
      data fpppp(12,33),fpppp(12,34)/ 2.47526865d-07, 2.15535136d-07/
      data fpppp(12,35),fpppp(12,36)/ 2.87808925d-07, 4.65066137d-07/
      data fpppp(12,37),fpppp(12,38)/-2.26783317d-07,-3.75136792d-07/
      data fpppp(12,39),fpppp(12,40)/-4.04731763d-07,-1.93292139d-07/
      data fpppp(12,41),fpppp(12,42)/ 3.13827932d-07, 7.36467882d-07/
      data fpppp(13, 1),fpppp(13, 2)/-3.42326912d-07,-1.48592371d-07/
      data fpppp(13, 3),fpppp(13, 4)/ 1.80973788d-07, 3.51292882d-08/
      data fpppp(13, 5),fpppp(13, 6)/ 3.68074335d-08, 1.98084769d-07/
      data fpppp(13, 7),fpppp(13, 8)/-2.36204063d-07, 1.43410629d-08/
      data fpppp(13, 9),fpppp(13,10)/-6.14333920d-08, 4.31483800d-08/
      data fpppp(13,11),fpppp(13,12)/-1.83105381d-07, 1.48678622d-07/
      data fpppp(13,13),fpppp(13,14)/ 1.53918969d-07,-3.20859745d-07/
      data fpppp(13,15),fpppp(13,16)/ 6.55086911d-08, 1.65848934d-07/
      data fpppp(13,17),fpppp(13,18)/-1.74401371d-07, 5.88807360d-08/
      data fpppp(13,19),fpppp(13,20)/ 1.27659582d-09,-2.82943975d-07/
      data fpppp(13,21),fpppp(13,22)/ 2.32014465d-07, 2.12981697d-07/
      data fpppp(13,23),fpppp(13,24)/-1.91248574d-07,-5.66953641d-08/
      data fpppp(13,25),fpppp(13,26)/-1.48556519d-08,-1.44337017d-07/
      data fpppp(13,27),fpppp(13,28)/ 1.58071424d-07,-1.37171076d-07/
      data fpppp(13,29),fpppp(13,30)/ 4.47803325d-07,-4.28607922d-07/
      data fpppp(13,31),fpppp(13,32)/ 1.34335167d-07, 8.10608561d-08/
      data fpppp(13,33),fpppp(13,34)/-1.91195496d-07, 1.15155536d-08/
      data fpppp(13,35),fpppp(13,36)/ 6.51970118d-09,-2.42032752d-07/
      data fpppp(13,37),fpppp(13,38)/ 2.05888701d-07, 2.89100177d-08/
      data fpppp(13,39),fpppp(13,40)/ 4.30077889d-08, 1.73321176d-07/
      data fpppp(13,41),fpppp(13,42)/-1.43350052d-07,-3.32311392d-07/
      data fpppp(14, 1),fpppp(14, 2)/ 9.40722539d-08, 4.13638407d-08/
      data fpppp(14, 3),fpppp(14, 4)/-3.99273391d-08,-2.11788363d-08/
      data fpppp(14, 5),fpppp(14, 6)/-2.26947913d-08,-4.42690628d-08/
      data fpppp(14, 7),fpppp(14, 8)/ 6.00736607d-08, 8.04863721d-09/
      data fpppp(14, 9),fpppp(14,10)/ 2.38254857d-08,-4.94692038d-09/
      data fpppp(14,11),fpppp(14,12)/ 4.18745248d-08,-4.21087708d-08/
      data fpppp(14,13),fpppp(14,14)/-2.34359635d-08, 2.95370533d-08/
      data fpppp(14,15),fpppp(14,16)/-6.74291189d-10,-3.83959449d-08/
      data fpppp(14,17),fpppp(14,18)/ 5.01338219d-08,-9.75931934d-09/
      data fpppp(14,19),fpppp(14,20)/-1.00019961d-08, 6.21050510d-08/
      data fpppp(14,21),fpppp(14,22)/-3.11343799d-08,-7.09011667d-08/
      data fpppp(14,23),fpppp(14,24)/ 5.03941570d-09, 8.56638576d-09/
      data fpppp(14,25),fpppp(14,26)/ 4.83294876d-08, 8.38395268d-08/
      data fpppp(14,27),fpppp(14,28)/-1.47455420d-08, 1.58941296d-08/
      data fpppp(14,29),fpppp(14,30)/-9.76119220d-08, 4.16522103d-08/
      data fpppp(14,31),fpppp(14,32)/-2.86010407d-08,-3.01381385d-08/
      data fpppp(14,33),fpppp(14,34)/ 4.49054707d-08, 9.05034081d-09/
      data fpppp(14,35),fpppp(14,36)/ 1.18711514d-08, 6.33816505d-08/
      data fpppp(14,37),fpppp(14,38)/-4.57974763d-08,-1.97160972d-08/
      data fpppp(14,39),fpppp(14,40)/-2.41479399d-08,-3.84756494d-08/
      data fpppp(14,41),fpppp(14,42)/ 3.83531559d-08, 8.91372429d-08/
      data fpppp(15, 1),fpppp(15, 2)/ 1.79803088d-07, 7.81967904d-08/
      data fpppp(15, 3),fpppp(15, 4)/-1.52776104d-07, 6.90340890d-08/
      data fpppp(15, 5),fpppp(15, 6)/ 6.91475550d-08,-1.62427194d-07/
      data fpppp(15, 7),fpppp(15, 8)/ 1.46924197d-07,-9.12567009d-08/
      data fpppp(15, 9),fpppp(15,10)/-8.69424115d-09,-7.84824536d-08/
      data fpppp(15,11),fpppp(15,12)/ 1.71546464d-07,-1.58781748d-07/
      data fpppp(15,13),fpppp(15,14)/ 2.45502192d-08, 9.13029427d-08/
      data fpppp(15,15),fpppp(15,16)/ 1.41190319d-08,-1.70608185d-07/
      data fpppp(15,17),fpppp(15,18)/ 1.27725834d-07,-4.49780204d-08/
      data fpppp(15,19),fpppp(15,20)/ 1.29389728d-08, 1.54067539d-07/
      data fpppp(15,21),fpppp(15,22)/-1.23165327d-07,-3.40672548d-08/
      data fpppp(15,23),fpppp(15,24)/ 2.47741641d-07, 1.34963942d-08/
      data fpppp(15,25),fpppp(15,26)/-6.72515816d-08,-1.27266285d-07/
      data fpppp(15,27),fpppp(15,28)/-1.18068591d-07, 2.07850622d-07/
      data fpppp(15,29),fpppp(15,30)/-1.72857711d-07, 1.62295185d-07/
      data fpppp(15,31),fpppp(15,32)/ 1.98861756d-08,-1.43401985d-07/
      data fpppp(15,33),fpppp(15,34)/ 1.50960109d-07,-6.87483495d-08/
      data fpppp(15,35),fpppp(15,36)/-8.54842739d-08, 1.49255721d-07/
      data fpppp(15,37),fpppp(15,38)/-1.71724464d-07, 7.37685999d-08/
      data fpppp(15,39),fpppp(15,40)/ 6.44182026d-08,-1.43520101d-07/
      data fpppp(15,41),fpppp(15,42)/ 7.60251764d-08, 1.73432288d-07/
 
      data x( 1), x( 2) /  3.20000000d+00 ,  3.40000000d+00 /
      data x( 3), x( 4) /  3.60000000d+00 ,  3.80000000d+00 /
      data x( 5), x( 6) /  4.00000000d+00 ,  4.20000000d+00 /
      data x( 7), x( 8) /  4.50000000d+00 ,  5.00000000d+00 /
      data x( 9), x(10) /  5.50000000d+00 ,  6.00000000d+00 /
      data x(11), x(12) /  7.00000000d+00 ,  8.00000000d+00 /
      data x(13), x(14) /  1.00000000d+01 ,  1.20000000d+01 /
      data x(15) /         1.50000000d+01 /
 
      data y( 1), y( 2) / -2.20000000d+02 , -2.10000000d+02 /
      data y( 3), y( 4) / -2.00000000d+02 , -1.90000000d+02 /
      data y( 5), y( 6) / -1.70000000d+02 , -1.60000000d+02 /
      data y( 7), y( 8) / -1.50000000d+02 , -1.40000000d+02 /
      data y( 9), y(10) / -1.30000000d+02 , -1.20000000d+02 /
      data y(11), y(12) / -1.10000000d+02 , -1.00000000d+02 /
      data y(13), y(14) / -9.00000000d+01 , -8.00000000d+01 /
      data y(15), y(16) / -7.00000000d+01 , -6.00000000d+01 /
      data y(17), y(18) / -5.00000000d+01 , -4.00000000d+01 /
      data y(19), y(20) / -3.00000000d+01 , -2.00000000d+01 /
      data y(21), y(22) / -1.00000000d+01 ,  1.00000000d+01 /
      data y(23), y(24) /  2.00000000d+01 ,  3.00000000d+01 /
      data y(25), y(26) /  4.00000000d+01 ,  5.00000000d+01 /
      data y(27), y(28) /  6.00000000d+01 ,  7.00000000d+01 /
      data y(29), y(30) /  8.00000000d+01 ,  9.00000000d+01 /
      data y(31), y(32) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(33), y(34) /  1.20000000d+02 ,  1.30000000d+02 /
      data y(35), y(36) /  1.40000000d+02 ,  1.50000000d+02 /
      data y(37), y(38) /  1.60000000d+02 ,  1.70000000d+02 /
      data y(39), y(40) /  1.90000000d+02 ,  2.00000000d+02 /
      data y(41), y(42) /  2.10000000d+02 ,  2.20000000d+02 /
 
      data delx( 1), delx( 2) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 3), delx( 4) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  3.00000000d-01 /
      data delx( 7), delx( 8) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 9), delx(10) /  5.00000000d-01 ,  1.00000000d+00 /
      data delx(11), delx(12) /  1.00000000d+00 ,  2.00000000d+00 /
      data delx(13), delx(14) /  2.00000000d+00 ,  3.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(15), dely(16) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(17), dely(18) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(19), dely(20) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(21), dely(22) /  2.00000000d+01 ,  1.00000000d+01 /
      data dely(23), dely(24) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(25), dely(26) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(27), dely(28) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(29), dely(30) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(31), dely(32) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(33), dely(34) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(35), dely(36) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(37), dely(38) /  1.00000000d+01 ,  2.00000000d+01 /
      data dely(39), dely(40) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(41) /            1.00000000d+01 /
      data nptx,npty /  15 , 42 /

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
c     if(iprint .gt. 2) then
c       write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
c    x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
c     endif
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
c     if(iprint .gt. 2) then
c       write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
c    x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
c     endif
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
      px(1)=((xi-xix)**3/(6.0d0*delxi))-(xi-xix)*delxi/6.0d0
      px(2)=(xi-xixp1)*delxi/6.0d0-((xi-xixp1)**3/(6.0d0*delxi))
      px(3)=(xi-xix)/delxi
      px(4)=(xixp1-xi)/delxi
      py(1)=((yi-yiy)**3/(6.0d0*delyi))-(yi-yiy)*delyi/6.0d0
      py(2)=(yi-yiyp1)*delyi/6.0d0-((yi-yiyp1)**3/(6.0d0*delyi))
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
