c***********************************************************************

      real*8 function ch2oh_h(r, rpar, ipar)

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
      data dmss(1,1), dmss(1,2) /12., 15.99491/

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

      cthe2 = (rax12 + rbx12 - rbond2)/(2.0d0*rax1*rbx1)
      if (dabs(cthe2).le.1.0d0) go to 2130
      if(dabs(cthe2)-1.d0.gt.1.e-6) 
     x   write (6,*) 'error in the2, cthe2 = ',cthe2
      the2 = 0.0d0
      if (cthe2.lt.1.0d0) the2 = 2.0d0*dasin(1.0d0)
      go to 2150
 2130 continue
c      write (6,*) 'cthe2 test',cthe2
      the2 = dacos(cthe2)
 2150 continue


      
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

      if (rmh.lt.3.8d0) then
c         write (6,*) 'error rmh too small',rmh,thetap,phi
         rmh=3.8d0
      endif
      if (rmh.gt.12.0d0) then
c         write (6,*) 'error rmh too large',rmh,thetap,phi
         rmh=12.0d0
      endif
      call hpch2oh(rmh,thetap,phi,vtot)
      ch2oh_h = vtot
c     write (6,*) 'pot test',rmh,thetap,phi,vtot*cautoicm/350.
c      if (phi.lt.0.0d0) then
c         vtot = 1.0d20
c      endif

      return
      end

c
c     test program to call h2coh + h 3d surface
c
c
c  ***  coordinate system  ***
c
c               
c                               r(cha)  = 2.10 au
c                               r(chb)  = 2.05 au
c       ha         hd           r(ohd)  = 1.85 au
c        \        /             r(c-o)  = 2.60 au
c         c===m==o              ha-c-o = 120.0 deg
c        /     \                hb-c-o = 115.0 deg
c       hb      \               hd-o-c = 105.0 deg 
c                \              habco  = 145.0 deg
c                 \             hacohd =  30.0 deg  
c                  hc           m is the center-of-mass of
c                               the co bond
c
c
c  ***   input   ***
c
c     rm :     r(m-hd) au
c     alpha :  hd-m-o angle  (degrees)
c     phi :    hd-m-o-ha dihedral angle
c
c     ranges :     3.8 < rd < 12.0
c                    0 < alpha < +180
c                  no limits on phi
c
c
c  ***   output   ***
c
c     energy - potential energy relative to h+h2coh in au
c
c

      subroutine hpch2oh(rm,alpha,phi,energy)
      implicit real*8 (a-h,o-z)
c
c     (2e,2o)-cas+1+2+dv/aug-cc-pvdz potential
c
      data ezero  / -115.2893987d0 /
      data degrad / 57.29577951308232d 00 /
c
c     find location in spline grids
c
      call fnd_grd_ch2oh_h(rm,alpha,ix,iy,delxi,delyi,
     $     xix,xixp1,yiy,yiyp1)
c
      call c0_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c0)
      call c1_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c1)
      call c2_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c2)
      call c3_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c3)
      call c4_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c4)
      call c5_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c5)
      call s1_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,s1)
      call s2_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,s2)
      call s3_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,s3)
      call s4_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,s4)
      call s5_spl_ch2oh_h(rm,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,s5)
c
      p = phi / degrad
      energy = c0 + c1*cos(1.0d0*p) + c2*cos(2.0d0*p) + c3*cos(3.0d0*p)
     x            + c4*cos(4.0d0*p) + c5*cos(5.0d0*p)
     x            + s1*sin(1.0d0*p) + s2*sin(2.0d0*p) + s3*sin(3.0d0*p)
     x            + s4*sin(4.0d0*p) + s5*sin(5.0d0*p)
     x            -  ezero

c      write(6,113) rm,alpha,phi,c0,c1,c2,c3,c4,c5,s1,s2,s3,s4,s5
c113   format(' coordinates:  ',3f12.6,/,
c     x       ' coefficients: ',1f12.6,/,15x,5f12.6,/,15x,5f12.6)

      return
      end
      subroutine fnd_grd_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,
     $     xixp1,yiy,yiyp1)
      implicit real*8 (a-h,o-z)
      dimension delx(12),dely(18),x(13),y(19)
 
 
      data x( 1), x( 2) /  3.80000000d+00 ,  4.00000000d+00 /
      data x( 3), x( 4) /  4.20000000d+00 ,  4.50000000d+00 /
      data x( 5), x( 6) /  5.00000000d+00 ,  5.50000000d+00 /
      data x( 7), x( 8) /  6.00000000d+00 ,  6.50000000d+00 /
      data x( 9), x(10) /  7.00000000d+00 ,  8.00000000d+00 /
      data x(11), x(12) /  9.00000000d+00 ,  1.00000000d+01 /
      data x(13) /         1.20000000d+01 /
 
      data y( 1), y( 2) /  0.00000000d+00 ,  1.00000000d+01 /
      data y( 3), y( 4) /  2.00000000d+01 ,  3.00000000d+01 /
      data y( 5), y( 6) /  4.00000000d+01 ,  5.00000000d+01 /
      data y( 7), y( 8) /  6.00000000d+01 ,  7.00000000d+01 /
      data y( 9), y(10) /  8.00000000d+01 ,  9.00000000d+01 /
      data y(11), y(12) /  1.00000000d+02 ,  1.10000000d+02 /
      data y(13), y(14) /  1.20000000d+02 ,  1.30000000d+02 /
      data y(15), y(16) /  1.40000000d+02 ,  1.50000000d+02 /
      data y(17), y(18) /  1.60000000d+02 ,  1.70000000d+02 /
      data y(19) /         1.80000000d+02 /
 
      data delx( 1), delx( 2) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 3), delx( 4) /  3.00000000d-01 ,  5.00000000d-01 /
      data delx( 5), delx( 6) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 7), delx( 8) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 9), delx(10) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx(11), delx(12) /  1.00000000d+00 ,  2.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(15), dely(16) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(17), dely(18) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  13 , 19 /
 

      if(xi .le. x(1)) then
        ix=1
      else if( xi .ge. x(nptx)) then
        ix=nptx-1
      else
        call hunt(x,nptx,xi,ix)
      endif
 15   xix=x(ix)
      xixp1=x(ix+1)
      delxi = delx(ix)
c      if(iprint .gt. 2) then
c        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' ix=',ix,
c     x       '  xix,xixp1=',xix,xixp1,'  delxi=',delxi
c      endif
c
      if(yi .le. y(1)) then
        iy=1
      else if( yi .ge. y(npty)) then
        iy=npty-1
      else
        call hunt(y,npty,yi,iy)
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
c      if(iprint .gt. 2) then
c        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
c     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
c      endif
c
      return
      end  
      subroutine c0_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -1.15244931d+02 , -1.15247255d+02 /
      data f( 1, 3),f( 1, 4) / -1.15252888d+02 , -1.15259315d+02 /
      data f( 1, 5),f( 1, 6) / -1.15264943d+02 , -1.15269283d+02 /
      data f( 1, 7),f( 1, 8) / -1.15271940d+02 , -1.15273589d+02 /
      data f( 1, 9),f( 1,10) / -1.15276513d+02 , -1.15282155d+02 /
      data f( 1,11),f( 1,12) / -1.15290058d+02 , -1.15297735d+02 /
      data f( 1,13),f( 1,14) / -1.15300761d+02 , -1.15292710d+02 /
      data f( 1,15),f( 1,16) / -1.15273925d+02 , -1.15259089d+02 /
      data f( 1,17),f( 1,18) / -1.15253771d+02 , -1.15232780d+02 /
      data f( 1,19) /          -1.15229320d+02 /
      data f( 2, 1),f( 2, 2) / -1.15256194d+02 , -1.15257920d+02 /
      data f( 2, 3),f( 2, 4) / -1.15262161d+02 , -1.15267117d+02 /
      data f( 2, 5),f( 2, 6) / -1.15271495d+02 , -1.15274683d+02 /
      data f( 2, 7),f( 2, 8) / -1.15276475d+02 , -1.15277674d+02 /
      data f( 2, 9),f( 2,10) / -1.15279995d+02 , -1.15284504d+02 /
      data f( 2,11),f( 2,12) / -1.15290961d+02 , -1.15297557d+02 /
      data f( 2,13),f( 2,14) / -1.15300725d+02 , -1.15294970d+02 /
      data f( 2,15),f( 2,16) / -1.15280841d+02 , -1.15267846d+02 /
      data f( 2,17),f( 2,18) / -1.15256707d+02 , -1.15235908d+02 /
      data f( 2,19) /          -1.15232784d+02 /
      data f( 3, 1),f( 3, 2) / -1.15264918d+02 , -1.15266184d+02 /
      data f( 3, 3),f( 3, 4) / -1.15269326d+02 , -1.15273046d+02 /
      data f( 3, 5),f( 3, 6) / -1.15276333d+02 , -1.15278639d+02 /
      data f( 3, 7),f( 3, 8) / -1.15279887d+02 , -1.15280811d+02 /
      data f( 3, 9),f( 3,10) / -1.15282657d+02 , -1.15286223d+02 /
      data f( 3,11),f( 3,12) / -1.15291445d+02 , -1.15297045d+02 /
      data f( 3,13),f( 3,14) / -1.15300269d+02 , -1.15296728d+02 /
      data f( 3,15),f( 3,16) / -1.15285794d+02 , -1.15273610d+02 /
      data f( 3,17),f( 3,18) / -1.15260877d+02 , -1.15244006d+02 /
      data f( 3,19) /          -1.15238982d+02 /
      data f( 4, 1),f( 4, 2) / -1.15274229d+02 , -1.15275012d+02 /
      data f( 4, 3),f( 4, 4) / -1.15276973d+02 , -1.15279317d+02 /
      data f( 4, 5),f( 4, 6) / -1.15281374d+02 , -1.15282771d+02 /
      data f( 4, 7),f( 4, 8) / -1.15283529d+02 , -1.15284176d+02 /
      data f( 4, 9),f( 4,10) / -1.15285472d+02 , -1.15287947d+02 /
      data f( 4,11),f( 4,12) / -1.15291667d+02 , -1.15295925d+02 /
      data f( 4,13),f( 4,14) / -1.15298941d+02 , -1.15297785d+02 /
      data f( 4,15),f( 4,16) / -1.15290818d+02 , -1.15280771d+02 /
      data f( 4,17),f( 4,18) / -1.15269081d+02 , -1.15256320d+02 /
      data f( 4,19) /          -1.15251009d+02 /
      data f( 5, 1),f( 5, 2) / -1.15282983d+02 , -1.15283328d+02 /
      data f( 5, 3),f( 5, 4) / -1.15284198d+02 , -1.15285237d+02 /
      data f( 5, 5),f( 5, 6) / -1.15286132d+02 , -1.15286716d+02 /
      data f( 5, 7),f( 5, 8) / -1.15287039d+02 , -1.15287372d+02 /
      data f( 5, 9),f( 5,10) / -1.15288053d+02 , -1.15289375d+02 /
      data f( 5,11),f( 5,12) / -1.15291409d+02 , -1.15293930d+02 /
      data f( 5,13),f( 5,14) / -1.15296088d+02 , -1.15296437d+02 /
      data f( 5,15),f( 5,16) / -1.15294048d+02 , -1.15287741d+02 /
      data f( 5,17),f( 5,18) / -1.15279620d+02 , -1.15271782d+02 /
      data f( 5,19) /          -1.15268256d+02 /
      data f( 6, 1),f( 6, 2) / -1.15286984d+02 , -1.15287136d+02 /
      data f( 6, 3),f( 6, 4) / -1.15287504d+02 , -1.15287942d+02 /
      data f( 6, 5),f( 6, 6) / -1.15288323d+02 , -1.15288551d+02 /
      data f( 6, 7),f( 6, 8) / -1.15288655d+02 , -1.15288817d+02 /
      data f( 6, 9),f( 6,10) / -1.15289123d+02 , -1.15289856d+02 /
      data f( 6,11),f( 6,12) / -1.15291023d+02 , -1.15292460d+02 /
      data f( 6,13),f( 6,14) / -1.15293778d+02 , -1.15294209d+02 /
      data f( 6,15),f( 6,16) / -1.15292885d+02 , -1.15289561d+02 /
      data f( 6,17),f( 6,18) / -1.15285063d+02 , -1.15280855d+02 /
      data f( 6,19) /          -1.15278983d+02 /
      data f( 7, 1),f( 7, 2) / -1.15288685d+02 , -1.15288751d+02 /
      data f( 7, 3),f( 7, 4) / -1.15288917d+02 , -1.15289107d+02 /
      data f( 7, 5),f( 7, 6) / -1.15289249d+02 , -1.15289317d+02 /
      data f( 7, 7),f( 7, 8) / -1.15289345d+02 , -1.15289395d+02 /
      data f( 7, 9),f( 7,10) / -1.15289563d+02 , -1.15289927d+02 /
      data f( 7,11),f( 7,12) / -1.15290528d+02 , -1.15291313d+02 /
      data f( 7,13),f( 7,14) / -1.15292076d+02 , -1.15292405d+02 /
      data f( 7,15),f( 7,16) / -1.15291846d+02 , -1.15290218d+02 /
      data f( 7,17),f( 7,18) / -1.15287945d+02 , -1.15285779d+02 /
      data f( 7,19) /          -1.15284845d+02 /
      data f( 8, 1),f( 8, 2) / -1.15289334d+02 , -1.15289362d+02 /
      data f( 8, 3),f( 8, 4) / -1.15289434d+02 , -1.15289512d+02 /
      data f( 8, 5),f( 8, 6) / -1.15289566d+02 , -1.15289583d+02 /
      data f( 8, 7),f( 8, 8) / -1.15289581d+02 , -1.15289589d+02 /
      data f( 8, 9),f( 8,10) / -1.15289666d+02 , -1.15289856d+02 /
      data f( 8,11),f( 8,12) / -1.15290188d+02 , -1.15290626d+02 /
      data f( 8,13),f( 8,14) / -1.15291073d+02 , -1.15291302d+02 /
      data f( 8,15),f( 8,16) / -1.15291069d+02 , -1.15290300d+02 /
      data f( 8,17),f( 8,18) / -1.15289197d+02 , -1.15288172d+02 /
      data f( 8,19) /          -1.15287740d+02 /
      data f( 9, 1),f( 9, 2) / -1.15289543d+02 , -1.15289556d+02 /
      data f( 9, 3),f( 9, 4) / -1.15289586d+02 , -1.15289619d+02 /
      data f( 9, 5),f( 9, 6) / -1.15289639d+02 , -1.15289640d+02 /
      data f( 9, 7),f( 9, 8) / -1.15289628d+02 , -1.15289622d+02 /
      data f( 9, 9),f( 9,10) / -1.15289653d+02 , -1.15289750d+02 /
      data f( 9,11),f( 9,12) / -1.15289937d+02 , -1.15290192d+02 /
      data f( 9,13),f( 9,14) / -1.15290456d+02 , -1.15290614d+02 /
      data f( 9,15),f( 9,16) / -1.15290541d+02 , -1.15290196d+02 /
      data f( 9,17),f( 9,18) / -1.15289682d+02 , -1.15289208d+02 /
      data f( 9,19) /          -1.15289012d+02 /
      data f(10, 1),f(10, 2) / -1.15289561d+02 , -1.15289564d+02 /
      data f(10, 3),f(10, 4) / -1.15289569d+02 , -1.15289574d+02 /
      data f(10, 5),f(10, 6) / -1.15289575d+02 , -1.15289571d+02 /
      data f(10, 7),f(10, 8) / -1.15289563d+02 , -1.15289555d+02 /
      data f(10, 9),f(10,10) / -1.15289559d+02 , -1.15289585d+02 /
      data f(10,11),f(10,12) / -1.15289644d+02 , -1.15289730d+02 /
      data f(10,13),f(10,14) / -1.15289827d+02 , -1.15289903d+02 /
      data f(10,15),f(10,16) / -1.15289921d+02 , -1.15289864d+02 /
      data f(10,17),f(10,18) / -1.15289758d+02 , -1.15289657d+02 /
      data f(10,19) /          -1.15289615d+02 /
      data f(11, 1),f(11, 2) / -1.15289497d+02 , -1.15289497d+02 /
      data f(11, 3),f(11, 4) / -1.15289498d+02 , -1.15289498d+02 /
      data f(11, 5),f(11, 6) / -1.15289497d+02 , -1.15289494d+02 /
      data f(11, 7),f(11, 8) / -1.15289489d+02 , -1.15289486d+02 /
      data f(11, 9),f(11,10) / -1.15289487d+02 , -1.15289497d+02 /
      data f(11,11),f(11,12) / -1.15289517d+02 , -1.15289546d+02 /
      data f(11,13),f(11,14) / -1.15289579d+02 , -1.15289612d+02 /
      data f(11,15),f(11,16) / -1.15289631d+02 , -1.15289627d+02 /
      data f(11,17),f(11,18) / -1.15289604d+02 , -1.15289578d+02 /
      data f(11,19) /          -1.15289567d+02 /
      data f(12, 1),f(12, 2) / -1.15289452d+02 , -1.15289452d+02 /
      data f(12, 3),f(12, 4) / -1.15289452d+02 , -1.15289452d+02 /
      data f(12, 5),f(12, 6) / -1.15289452d+02 , -1.15289450d+02 /
      data f(12, 7),f(12, 8) / -1.15289447d+02 , -1.15289446d+02 /
      data f(12, 9),f(12,10) / -1.15289448d+02 , -1.15289453d+02 /
      data f(12,11),f(12,12) / -1.15289461d+02 , -1.15289472d+02 /
      data f(12,13),f(12,14) / -1.15289483d+02 , -1.15289495d+02 /
      data f(12,15),f(12,16) / -1.15289504d+02 , -1.15289508d+02 /
      data f(12,17),f(12,18) / -1.15289504d+02 , -1.15289498d+02 /
      data f(12,19) /          -1.15289495d+02 /
      data f(13, 1),f(13, 2) / -1.15289416d+02 , -1.15289416d+02 /
      data f(13, 3),f(13, 4) / -1.15289416d+02 , -1.15289417d+02 /
      data f(13, 5),f(13, 6) / -1.15289417d+02 , -1.15289416d+02 /
      data f(13, 7),f(13, 8) / -1.15289415d+02 , -1.15289415d+02 /
      data f(13, 9),f(13,10) / -1.15289415d+02 , -1.15289417d+02 /
      data f(13,11),f(13,12) / -1.15289419d+02 , -1.15289421d+02 /
      data f(13,13),f(13,14) / -1.15289424d+02 , -1.15289426d+02 /
      data f(13,15),f(13,16) / -1.15289428d+02 , -1.15289430d+02 /
      data f(13,17),f(13,18) / -1.15289431d+02 , -1.15289431d+02 /
      data f(13,19) /          -1.15289432d+02 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 7.47320601d-02,-5.79824742d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 7.06890031d-02,-3.36050517d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 6.24165010d-02,-6.13731919d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 5.68958081d-02, 1.05143284d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 5.40547977d-02, 1.20200055d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 4.62292836d-02, 1.86856496d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 3.51286360d-02, 1.42173962d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.86301742d-02,-1.50752344d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.45744629d-02,-3.04164587d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.94307723d-02,-2.63389307d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.36434154d-02, 1.12181666d-07/
      data fpp( 1,12,1),fpp( 1,12,2)/ 1.14602828d-02, 3.94502041d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.22395528d-02, 1.21147002d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/ 5.80875100d-03, 1.40581788d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 6.30293363d-02,-3.94341531d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 1.18713704d-01,-2.19785175d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/-3.70433409d-02, 3.47494855d-04/
      data fpp( 1,18,1),fpp( 1,18,2)/-2.23799556d-01,-2.29814244d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/-9.66752358d-02,-4.80097878d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 6.34894513d-02,-4.30092673d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 6.00427081d-02,-2.54314653d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 5.27412837d-02,-6.16487139d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 4.70119552d-02, 7.19095089d-06/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 4.32461190d-02, 1.20810678d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 3.65385756d-02, 1.58847778d-05/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.83591566d-02, 8.13982100d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 2.38282230d-02,-1.28640618d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 2.05682170d-02,-2.40035738d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.58570269d-02,-2.24016429d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.05753121d-02,-3.26985439d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 8.39800586d-03, 2.71410605d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.03008943d-02, 1.00385612d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.14324980d-02, 1.06696490d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 5.03120418d-02,-2.47315734d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 8.13761644d-02,-7.58101968d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/-2.68004610d-02, 2.16612361d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.34333745d-01,-2.11039246d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/-6.64709570d-02,-4.32955377d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 5.21601349d-02,-3.16636502d-05/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 4.92901646d-02,-1.89126996d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 4.28183640d-02,-5.24555128d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 3.60063711d-02, 5.21490476d-06/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 3.00607265d-02, 1.03659322d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 2.42164138d-02, 1.21813662d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.98847377d-02, 4.38860277d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 1.82569336d-02,-1.02957773d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 1.61526691d-02,-1.85254935d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.16411201d-02,-1.88022487d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/ 6.90533622d-03,-5.62551157d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 5.04769376d-03, 1.86242950d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 9.55686999d-03, 7.36883315d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.37612571d-02, 9.25223790d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.01724967d-02,-1.97847405d-07/
      data fpp( 3,16,1),fpp( 3,16,2)/ 4.73163883d-03,-1.67309894d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-4.08548149d-02, 1.00061805d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.56345357d-02,-1.35236230d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/-4.75409363d-02,-2.69936885d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 3.54732495d-02,-1.97029105d-05/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 3.35376461d-02,-1.18341791d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 2.88112640d-02,-3.64037323d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 2.34707928d-02, 3.41567199d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.87001658d-02, 7.19768526d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 1.50529036d-02, 7.39358695d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 1.32114365d-02, 1.56796692d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 1.26247392d-02,-7.00545465d-06/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 1.09789582d-02,-1.24861483d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/ 7.59158176d-03,-1.37899521d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/ 3.53200452d-03,-7.05404340d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.04235022d-03, 9.72612568d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 4.20983717d-03, 4.26695407d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.85074778d-02, 6.99157115d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.62503165d-02, 2.63276132d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/ 2.89770943d-02, 9.57383563d-06/
      data fpp( 4,17,1),fpp( 4,17,2)/ 2.41163571d-02, 3.39570443d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 2.63073777d-02,-8.11420126d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 2.07837591d-02,-1.56388994d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 1.75335206d-02,-8.80355354d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.66414338d-02,-5.26289293d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.45929367d-02,-1.64487475d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 1.20496404d-02, 1.70239192d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 9.57103350d-03, 3.47530708d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 7.90086033d-03, 3.05637976d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 7.23256041d-03,-4.08261307d-08/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 6.54267434d-03,-3.49307524d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 5.83173230d-03,-6.86687289d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/ 3.41026633d-03,-7.49943318d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/-3.73616200d-04,-5.85539438d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/-3.28413695d-03, 1.70101071d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/-3.85360092d-03, 2.08313516d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 1.15131690d-03, 2.35135831d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 2.13354893d-02, 4.93943161d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 2.35943150d-02, 1.39891524d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 2.25645462d-02, 3.48907436d-06/
      data fpp( 5,18,1),fpp( 5,18,2)/ 2.79076698d-02,-4.49254498d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 2.91685328d-02,-8.25072751d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 8.46466813d-03,-3.61342201d-06/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 8.08861871d-03,-2.17315597d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 6.87298931d-03,-6.53954100d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 5.49064556d-03, 5.88972371d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 4.62370017d-03, 1.71806462d-06/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 3.98365512d-03, 1.71876917d-06/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 3.31432181d-03,-1.15314129d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 3.22856341d-03,-5.86204019d-07/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.95811262d-03,-5.14204264d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 1.49535293d-03,-4.46562544d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 1.03446028d-03,-3.03545562d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/-5.05802418d-04, 4.07447905d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.82743347d-03, 8.54566400d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/-1.99274539d-03, 1.86298961d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/-6.16027348d-03, 2.22347516d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/ 2.45645702d-04, 1.24310975d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 7.92945805d-03,-1.51914164d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/ 1.53979430d-02,-2.37545310d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 1.90221097d-02,-4.36227345d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 3.80780690d-03,-1.76180876d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 3.63609135d-03,-9.96382487d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 3.34710609d-03,-2.52661297d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 2.94777735d-03, 5.67027675d-07/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 2.29416583d-03, 8.64550601d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 1.82051920d-03, 4.14769925d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 1.73415235d-03,-1.23630298d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.35107200d-03,-1.24024874d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.45581724d-03,-1.99537476d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/ 4.48321949d-04,-2.53825223d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-1.14822492d-03,-2.07161632d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/-2.44465338d-03,-2.15282484d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/-3.42866520d-03, 4.25274626d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/-3.35633537d-03, 9.24429745d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 3.29604684d-04, 1.20500639d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 3.33510219d-03, 6.69544678d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 7.18162159d-03,-1.31851056d-07/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.00765582d-02,-1.25880426d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 1.15030286d-02,-2.34359787d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 1.55210427d-03,-8.17664427d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.46301591d-03,-4.44671147d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 1.24258632d-03,-4.36509889d-08/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 9.58245028d-04, 2.59275103d-07/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 8.15636500d-04, 4.46550578d-07/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 7.34268097d-04, 1.74522586d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 6.45068797d-04,-4.64092075d-09/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 5.83148592d-04,-7.55958905d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 3.06618442d-04,-1.11152346d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.19359273d-04,-1.57794726d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-1.61560587d-04,-1.09668750d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-7.55584057d-04,-3.95302738d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.23390574d-03, 2.13789845d-06/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.40591315d-03, 4.92370893d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/-1.44614525d-03, 5.88726584d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.13945544d-04, 3.68722771d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 2.46405559d-03,-5.96176693d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/ 5.03982404d-03,-5.98252095d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 6.17377611d-03,-1.10537395d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 5.43776023d-04,-3.11658623d-07/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 5.19845027d-04,-1.66682753d-07/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 4.42548638d-04,-4.16103650d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 3.71242533d-04, 1.53124215d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 2.99288166d-04, 2.09113509d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 2.58408415d-04, 1.50421750d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.21572462d-04,-3.08005059d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 1.80333634d-04,-3.87219729d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.01708996d-04,-6.40320576d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/-8.57590424d-05,-1.01149797d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-3.41532728d-04,-7.13687562d-07/
      data fpp( 9,12,1),fpp( 9,12,2)/-6.05010392d-04,-2.13751783d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/-8.99711836d-04, 1.02869470d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-9.80012026d-04, 2.45897300d-06/
      data fpp( 9,15,1),fpp( 9,15,2)/-5.21023678d-04, 2.99541331d-06/
      data fpp( 9,16,1),fpp( 9,16,2)/ 2.73115635d-04, 1.87937375d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/ 1.37015607d-03,-3.72908318d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 2.33214560d-03,-2.78774048d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 2.75386700d-03,-5.15612976d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/-7.38020335d-06,-4.11076265d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.10430356d-05,-1.77847477d-08/
      data fpp(10, 3,1),fpp(10, 3,2)/-2.29390721d-05,-7.75338569d-09/
      data fpp(10, 4,1),fpp(10, 4,2)/-3.88501146d-05, 4.87982922d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-4.56827487d-05, 5.25602201d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-4.43592943d-05, 4.09608313d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-3.32517857d-05, 2.35964550d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/-3.45751979d-05,-1.35346654d-07/
      data fpp(10, 9,1),fpp(10, 9,2)/-5.04362104d-05,-2.02209839d-07/
      data fpp(10,10,1),fpp(10,10,2)/-8.44025098d-05,-3.75813992d-07/
      data fpp(10,11,1),fpp(10,11,2)/-1.48621524d-04,-2.74534191d-07/
      data fpp(10,12,1),fpp(10,12,2)/-2.43176797d-04,-1.46049243d-07/
      data fpp(10,13,1),fpp(10,13,2)/-3.13911621d-04, 1.98731164d-07/
      data fpp(10,14,1),fpp(10,14,2)/-3.47007345d-04, 6.11124588d-07/
      data fpp(10,15,1),fpp(10,15,2)/-3.29856340d-04, 8.36770484d-07/
      data fpp(10,16,1),fpp(10,16,2)/-1.82319677d-04, 5.41793478d-07/
      data fpp(10,17,1),fpp(10,17,2)/ 2.15040106d-05,-6.39443984d-08/
      data fpp(10,18,1),fpp(10,18,2)/ 2.21651177d-04,-5.86015886d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 2.97510943d-04,-1.13199206d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/-2.22552096d-05,-2.90686072d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-2.56728854d-05,-1.18627863d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-2.67923494d-05, 1.65197511d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/-2.98420751d-05, 5.78378340d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-3.25571714d-05, 2.03451165d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-3.29712383d-05, 3.28357527d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-3.45653195d-05,-3.16881251d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/-3.00328426d-05,-2.60832543d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/-3.19641552d-05,-1.03978857d-07/
      data fpp(11,10,1),fpp(11,10,2)/-3.86309184d-05,-9.80013151d-08/
      data fpp(11,11,1),fpp(11,11,2)/-5.99811779d-05,-1.04015880d-07/
      data fpp(11,12,1),fpp(11,12,2)/-9.02824219d-05,-2.59351611d-08/
      data fpp(11,13,1),fpp(11,13,2)/-1.30641679d-04,-3.22434741d-08/
      data fpp(11,14,1),fpp(11,14,2)/-1.51958593d-04, 1.54909058d-07/
      data fpp(11,15,1),fpp(11,15,2)/-1.39550963d-04, 2.52607242d-07/
      data fpp(11,16,1),fpp(11,16,2)/-1.13836927d-04, 2.14661977d-07/
      data fpp(11,17,1),fpp(11,17,2)/-7.61721081d-05, 2.87448519d-08/
      data fpp(11,18,1),fpp(11,18,2)/-5.07503097d-05,-1.49641387d-07/
      data fpp(11,19,1),fpp(11,19,2)/-3.79107746d-05,-3.30179306d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/-1.75989581d-05, 2.91723900d-10/
      data fpp(12, 2,1),fpp(12, 2,2)/-1.82654229d-05,-5.83448797d-10/
      data fpp(12, 3,1),fpp(12, 3,2)/-1.98915301d-05, 2.04206848d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/-2.17815850d-05,-7.58482331d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/-2.20885657d-05, 2.82972285d-08/
      data fpp(12, 6,1),fpp(12, 6,2)/-2.17557523d-05, 1.43959127d-08/
      data fpp(12, 7,1),fpp(12, 7,2)/-2.04869361d-05,-2.58808763d-08/
      data fpp(12, 8,1),fpp(12, 8,2)/-1.92934314d-05,-3.08724096d-08/
      data fpp(12, 9,1),fpp(12, 9,2)/-1.97071689d-05,-3.06294851d-08/
      data fpp(12,10,1),fpp(12,10,2)/-2.50738163d-05,-2.66096517d-08/
      data fpp(12,11,1),fpp(12,11,2)/-3.74537644d-05,-4.29319072d-08/
      data fpp(12,12,1),fpp(12,12,2)/-5.56935156d-05, 1.83372844d-08/
      data fpp(12,13,1),fpp(12,13,2)/-7.55216643d-05,-3.04172294d-08/
      data fpp(12,14,1),fpp(12,14,2)/-8.91582813d-05, 4.33316340d-08/
      data fpp(12,15,1),fpp(12,15,2)/-8.99398073d-05, 3.70906945d-08/
      data fpp(12,16,1),fpp(12,16,2)/-7.03326145d-05, 1.08305588d-07/
      data fpp(12,17,1),fpp(12,17,2)/-4.08155784d-05, 9.68695650d-09/
      data fpp(12,18,1),fpp(12,18,2)/-1.26499380d-05,-2.70534166d-08/
      data fpp(12,19,1),fpp(12,19,2)/-1.86784507d-06,-8.14732911d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/-1.70755209d-05, 8.55739950d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/-1.33672885d-05, 2.88519962d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/-1.09292349d-05,-2.00981995d-08/
      data fpp(13, 4,1),fpp(13, 4,2)/-5.23420748d-06, 1.75076002d-08/
      data fpp(13, 5,1),fpp(13, 5,2)/ 4.42828800d-08, 1.00678010d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/ 7.52876180d-07, 2.22119772d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 7.43468069d-07,-1.89525909d-08/
      data fpp(13, 8,1),fpp(13, 8,2)/-6.03284200d-07, 1.35891635d-08/
      data fpp(13, 9,1),fpp(13, 9,2)/ 7.60358452d-06,-3.54040632d-08/
      data fpp(13,10,1),fpp(13,10,2)/ 1.65369082d-05, 8.02708972d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 3.73518823d-05, 3.29570634d-09/
      data fpp(13,12,1),fpp(13,12,2)/ 6.67217578d-05,-2.12099133d-08/
      data fpp(13,13,1),fpp(13,13,2)/ 9.23858322d-05, 2.15439493d-08/
      data fpp(13,14,1),fpp(13,14,2)/ 9.59541407d-05,-4.96588436d-09/
      data fpp(13,15,1),fpp(13,15,2)/ 7.25949037d-05,-1.68041205d-09/
      data fpp(13,16,1),fpp(13,16,2)/ 2.79163073d-05, 1.16875337d-08/
      data fpp(13,17,1),fpp(13,17,2)/-2.99672107d-05, 1.49302792d-08/
      data fpp(13,18,1),fpp(13,18,2)/-7.61750310d-05,-1.14086516d-08/
      data fpp(13,19,1),fpp(13,19,2)/-9.69410774d-05,-2.92956751d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-1.10236113d-04,-4.60362151d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/ 4.06142659d-05, 4.86877068d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/-7.45841490d-05,-4.94213264d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/ 7.57614418d-05, 2.25067068d-05/
      data fpppp( 1, 9),fpppp( 1,10)/-1.92232357d-05,-1.08925295d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 2.41733819d-05, 1.30452459d-04/
      data fpppp( 1,13),fpppp( 1,14)/-3.68239056d-04, 9.09899451d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 5.47724477d-04,-3.19297044d-03/
      data fpppp( 1,17),fpppp( 1,18)/-4.62327413d-04, 3.18232984d-03/
      data fpppp( 1,19) /             6.56584018d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-9.14608757d-05,-4.12542229d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 2.51968969d-05, 3.47923827d-05/
      data fpppp( 2, 5),fpppp( 2, 6)/-4.65568873d-05,-2.50672599d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/ 5.85133814d-05, 9.92286748d-06/
      data fpppp( 2, 9),fpppp( 2,10)/-2.19492039d-05,-9.19709170d-06/
      data fpppp( 2,11),fpppp( 2,12)/ 2.45060859d-05, 9.74372628d-05/
      data fpppp( 2,13),fpppp( 2,14)/-1.69443458d-04, 5.34059484d-04/
      data fpppp( 2,15),fpppp( 2,16)/ 2.98081928d-04,-2.19531246d-03/
      data fpppp( 2,17),fpppp( 2,18)/ 1.28723033d-04, 1.71902083d-03/
      data fpppp( 2,19) /             3.51895796d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-6.80588342d-05,-3.71699927d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 6.28994122d-07, 1.42424723d-05/
      data fpppp( 3, 5),fpppp( 3, 6)/-5.61798826d-06, 1.43094014d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/ 3.91385772d-05,-8.63139112d-06/
      data fpppp( 3, 9),fpppp( 3,10)/-3.32006359d-05,-3.00313992d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 3.17591069d-05, 4.86551968d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 1.55629227d-04,-8.94594555d-05/
      data fpppp( 3,15),fpppp( 3,16)/-2.65380252d-04,-7.60145385d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 2.09722604d-03,-1.50421051d-03/
      data fpppp( 3,19) /            -3.26027335d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-4.93423698d-05,-2.85724067d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/-3.81471958d-06, 6.98593299d-06/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.00616447d-05, 2.01693705d-05/
      data fpppp( 4, 7),fpppp( 4, 8)/ 1.76085885d-05,-1.53175442d-05/
      data fpppp( 4, 9),fpppp( 4,10)/-1.98834329d-05,-9.64444929d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 1.81291821d-05, 3.13230965d-05/
      data fpppp( 4,13),fpppp( 4,14)/ 1.96006907d-04,-1.47541506d-04/
      data fpppp( 4,15),fpppp( 4,16)/ 8.71002221d-07,-1.56906158d-04/
      data fpppp( 4,17),fpppp( 4,18)/ 1.71502730d-04,-1.05999294d-04/
      data fpppp( 4,19) /            -2.10383915d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.81141389d-05,-1.16962575d-05/
      data fpppp( 5, 3),fpppp( 5, 4)/-4.48545230d-06,-4.98800600d-08/
      data fpppp( 5, 5),fpppp( 5, 6)/ 8.56633286d-06, 1.42905728d-05/
      data fpppp( 5, 7),fpppp( 5, 8)/-5.61622916d-06, 6.87917504d-06/
      data fpppp( 5, 9),fpppp( 5,10)/-2.31638297d-05,-1.68552916d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 8.84000271d-06, 3.38969876d-05/
      data fpppp( 5,13),fpppp( 5,14)/-3.96454661d-06, 3.16424107d-04/
      data fpppp( 5,15),fpppp( 5,16)/-3.50976610d-04, 1.19615378d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 1.05814786d-04,-5.28471380d-05/
      data fpppp( 5,19) /            -1.39361870d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/-1.52224842d-05,-8.19974394d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/-2.35333925d-06, 7.61023995d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 2.83628108d-06,-5.34134358d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.67717977d-05,-2.67313524d-05/
      data fpppp( 6, 9),fpppp( 6,10)/ 1.90720677d-05,-1.09545160d-06/
      data fpppp( 6,11),fpppp( 6,12)/-1.45782391d-05,-5.35379492d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 4.91113175d-05,-1.21712327d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 1.97605019d-04,-3.43009121d-05/
      data fpppp( 6,17),fpppp( 6,18)/ 1.62722192d-05,-4.37076093d-05/
      data fpppp( 6,19) /            -7.21008779d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.12036248d-06,-1.41589010d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-2.52258945d-07,-4.19568337d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.77802560d-06, 7.88147393d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.00671339d-05, 1.45842518d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-1.90003381d-05,-5.31733088d-06/
      data fpppp( 7,11),fpppp( 7,12)/ 4.92656644d-06, 3.61817000d-06/
      data fpppp( 7,13),fpppp( 7,14)/-6.54247910d-07, 6.23793204d-05/
      data fpppp( 7,15),fpppp( 7,16)/-3.20464208d-05, 2.49798099d-05/
      data fpppp( 7,17),fpppp( 7,18)/-1.74115048d-05,-1.24287556d-05/
      data fpppp( 7,19) /            -2.09814532d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.07532010d-06,-1.13818717d-06/
      data fpppp( 8, 3),fpppp( 8, 4)/-1.25240483d-06, 2.31310454d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 5.03952340d-07,-6.54506397d-07/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.64421944d-06,-4.28562573d-06/
      data fpppp( 8, 9),fpppp( 8,10)/ 2.62168679d-06,-8.44862573d-07/
      data fpppp( 8,11),fpppp( 8,12)/-4.86187803d-06, 1.50615815d-06/
      data fpppp( 8,13),fpppp( 8,14)/ 5.77935254d-06,-6.24471178d-06/
      data fpppp( 8,15),fpppp( 8,16)/ 2.71060132d-05,-1.59967324d-07/
      data fpppp( 8,17),fpppp( 8,18)/ 8.93501087d-06,-1.60405715d-05/
      data fpppp( 8,19) /            -3.12817079d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.09368174d-06,-6.00711948d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/ 2.94605862d-07,-2.18294389d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/ 5.39675916d-07,-7.59322853d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/ 6.68110091d-09,-2.14964637d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-1.38997110d-06,-7.55755038d-07/
      data fpppp( 9,11),fpppp( 9,12)/ 3.14652469d-07,-9.65093563d-07/
      data fpppp( 9,13),fpppp( 9,14)/ 1.67229497d-06, 7.13998889d-06/
      data fpppp( 9,15),fpppp( 9,16)/ 2.12506177d-06, 4.46882191d-06/
      data fpppp( 9,17),fpppp( 9,18)/-1.82628238d-06,-5.26674610d-06/
      data fpppp( 9,19) /            -9.52282142d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-1.30058416d-07,-7.12432636d-08/
      data fpppp(10, 3),fpppp(10, 4)/-7.89607793d-08, 1.46186026d-07/
      data fpppp(10, 5),fpppp(10, 6)/ 3.89211756d-08, 1.87494578d-07/
      data fpppp(10, 7),fpppp(10, 8)/-2.01856234d-07,-1.25924882d-07/
      data fpppp(10, 9),fpppp(10,10)/-1.66700258d-07,-2.93591297d-07/
      data fpppp(10,11),fpppp(10,12)/-4.74097415d-07, 3.69805403d-07/
      data fpppp(10,13),fpppp(10,14)/ 4.24102689d-07, 1.92129909d-07/
      data fpppp(10,15),fpppp(10,16)/ 1.82218142d-06, 3.42283842d-07/
      data fpppp(10,17),fpppp(10,18)/ 1.85904719d-07,-1.30649398d-06/
      data fpppp(10,19) /            -2.41717279d-06 /
      data fpppp(11, 1),fpppp(11, 2)/ 6.28435378d-08, 2.78287386d-08/
      data fpppp(11, 3),fpppp(11, 4)/-3.62658006d-08, 1.41876890d-09/
      data fpppp(11, 5),fpppp(11, 6)/ 5.06684957d-08,-6.60309853d-08/
      data fpppp(11, 7),fpppp(11, 8)/ 1.42654589d-07,-1.36993887d-07/
      data fpppp(11, 9),fpppp(11,10)/ 1.74935948d-08,-2.17107548d-07/
      data fpppp(11,11),fpppp(11,12)/-3.00731772d-08,-1.99658806d-07/
      data fpppp(11,13),fpppp(11,14)/ 2.25227642d-07, 4.41288749d-07/
      data fpppp(11,15),fpppp(11,16)/ 3.30900636d-08, 2.24735355d-07/
      data fpppp(11,17),fpppp(11,18)/-2.14984515d-07,-9.93785427d-08/
      data fpppp(11,19) /            -1.42437103d-07 /
      data fpppp(12, 1),fpppp(12, 2)/-1.73878639d-08,-7.92743517d-09/
      data fpppp(12, 3),fpppp(12, 4)/-8.48093428d-09, 2.60143114d-08/
      data fpppp(12, 5),fpppp(12, 6)/-5.91865003d-10, 1.47407957d-08/
      data fpppp(12, 7),fpppp(12, 8)/-2.21114777d-09,-1.04149015d-08/
      data fpppp(12, 9),fpppp(12,10)/-5.25637736d-08,-7.65045929d-08/
      data fpppp(12,11),fpppp(12,12)/-6.22158995d-08,-2.62199973d-08/
      data fpppp(12,13),fpppp(12,14)/ 7.17920413d-08, 1.10543730d-07/
      data fpppp(12,15),fpppp(12,16)/ 2.57338497d-07, 8.34254105d-08/
      data fpppp(12,17),fpppp(12,18)/ 3.55046575d-09,-1.78711024d-07/
      data fpppp(12,19) /            -3.31719214d-07 /
      data fpppp(13, 1),fpppp(13, 2)/-5.49406572d-08,-1.87671036d-08/
      data fpppp(13, 3),fpppp(13, 4)/ 5.37983489d-08,-1.00786327d-09/
      data fpppp(13, 5),fpppp(13, 6)/-7.47591300d-08, 2.58505529d-08/
      data fpppp(13, 7),fpppp(13, 8)/-7.17231692d-08, 1.80801474d-07/
      data fpppp(13, 9),fpppp(13,10)/-7.82654609d-08, 1.75847675d-07/
      data fpppp(13,11),fpppp(13,12)/ 8.77737878d-08,-1.36487385d-08/
      data fpppp(13,13),fpppp(13,14)/-2.55526911d-07,-2.89989571d-07/
      data fpppp(13,15),fpppp(13,16)/-2.00167541d-07,-1.88501822d-07/
      data fpppp(13,17),fpppp(13,18)/ 1.61879530d-07, 2.41525579d-07/
      data fpppp(13,19) /             3.98524574d-07 /
 
      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine c1_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  4.27420000d-03 /
      data f( 1, 3),f( 1, 4) /  7.85310000d-03 ,  9.99630000d-03 /
      data f( 1, 5),f( 1, 6) /  1.01891000d-02 ,  8.41390000d-03 /
      data f( 1, 7),f( 1, 8) /  6.21330000d-03 ,  4.81850000d-03 /
      data f( 1, 9),f( 1,10) /  3.68420000d-03 ,  2.60100000d-03 /
      data f( 1,11),f( 1,12) /  1.83390000d-03 ,  2.00440000d-03 /
      data f( 1,13),f( 1,14) /  3.91510000d-03 ,  1.01059000d-02 /
      data f( 1,15),f( 1,16) /  1.87949000d-02 ,  2.85159000d-02 /
      data f( 1,17),f( 1,18) /  2.53894000d-02 ,  2.02553000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  2.97980000d-03 /
      data f( 2, 3),f( 2, 4) /  5.51020000d-03 ,  6.99440000d-03 /
      data f( 2, 5),f( 2, 6) /  7.04010000d-03 ,  5.89810000d-03 /
      data f( 2, 7),f( 2, 8) /  4.61110000d-03 ,  3.71830000d-03 /
      data f( 2, 9),f( 2,10) /  2.84300000d-03 ,  1.96140000d-03 /
      data f( 2,11),f( 2,12) /  1.39560000d-03 ,  1.47270000d-03 /
      data f( 2,13),f( 2,14) /  2.78790000d-03 ,  7.67880000d-03 /
      data f( 2,15),f( 2,16) /  1.58071000d-02 ,  2.28440000d-02 /
      data f( 2,17),f( 2,18) /  2.14924000d-02 ,  1.99694000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  2.06790000d-03 /
      data f( 3, 3),f( 3, 4) /  3.86060000d-03 ,  4.92550000d-03 /
      data f( 3, 5),f( 3, 6) /  4.96810000d-03 ,  4.25360000d-03 /
      data f( 3, 7),f( 3, 8) /  3.46390000d-03 ,  2.84240000d-03 /
      data f( 3, 9),f( 3,10) /  2.16380000d-03 ,  1.48780000d-03 /
      data f( 3,11),f( 3,12) /  1.06420000d-03 ,  1.09530000d-03 /
      data f( 3,13),f( 3,14) /  1.97470000d-03 ,  5.42290000d-03 /
      data f( 3,15),f( 3,16) /  1.20532000d-02 ,  1.81055000d-02 /
      data f( 3,17),f( 3,18) /  1.87068000d-02 ,  1.63668000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  1.17870000d-03 /
      data f( 4, 3),f( 4, 4) /  2.24300000d-03 ,  2.91490000d-03 /
      data f( 4, 5),f( 4, 6) /  2.99790000d-03 ,  2.65680000d-03 /
      data f( 4, 7),f( 4, 8) /  2.23760000d-03 ,  1.84410000d-03 /
      data f( 4, 9),f( 4,10) /  1.39040000d-03 ,  9.53600000d-04 /
      data f( 4,11),f( 4,12) /  7.04200000d-04 ,  7.28300000d-04 /
      data f( 4,13),f( 4,14) /  1.21700000d-03 ,  3.06600000d-03 /
      data f( 4,15),f( 4,16) /  7.29230000d-03 ,  1.19103000d-02 /
      data f( 4,17),f( 4,18) /  1.39388000d-02 ,  1.14510000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  4.25100000d-04 /
      data f( 5, 3),f( 5, 4) /  8.47500000d-04 ,  1.15660000d-03 /
      data f( 5, 5),f( 5, 6) /  1.24500000d-03 ,  1.15000000d-03 /
      data f( 5, 7),f( 5, 8) /  9.88300000d-04 ,  8.00500000d-04 /
      data f( 5, 9),f( 5,10) /  5.95900000d-04 ,  4.04800000d-04 /
      data f( 5,11),f( 5,12) /  3.34300000d-04 ,  3.35400000d-04 /
      data f( 5,13),f( 5,14) /  6.60400000d-04 ,  1.39980000d-03 /
      data f( 5,15),f( 5,16) /  4.81560000d-03 ,  6.64070000d-03 /
      data f( 5,17),f( 5,18) /  7.77890000d-03 ,  6.12070000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.15400000d-04 /
      data f( 6, 3),f( 6, 4) /  2.69400000d-04 ,  4.01900000d-04 /
      data f( 6, 5),f( 6, 6) /  4.46100000d-04 ,  4.21300000d-04 /
      data f( 6, 7),f( 6, 8) /  3.69100000d-04 ,  2.68200000d-04 /
      data f( 6, 9),f( 6,10) /  2.28700000d-04 ,  2.33000000d-04 /
      data f( 6,11),f( 6,12) /  2.87400000d-04 ,  3.35100000d-04 /
      data f( 6,13),f( 6,14) /  5.12500000d-04 ,  8.98900000d-04 /
      data f( 6,15),f( 6,16) /  1.83020000d-03 ,  3.33570000d-03 /
      data f( 6,17),f( 6,18) /  4.26790000d-03 ,  3.16490000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.02000000d-05 /
      data f( 7, 3),f( 7, 4) /  3.92000000d-05 ,  7.70000000d-05 /
      data f( 7, 5),f( 7, 6) /  1.00100000d-04 ,  9.93000000d-05 /
      data f( 7, 7),f( 7, 8) /  7.36000000d-05 ,  4.81000000d-05 /
      data f( 7, 9),f( 7,10) /  2.99000000d-05 ,  2.13000000d-05 /
      data f( 7,11),f( 7,12) /  3.88000000d-05 ,  9.62000000d-05 /
      data f( 7,13),f( 7,14) /  2.32700000d-04 ,  5.24900000d-04 /
      data f( 7,15),f( 7,16) /  1.01560000d-03 ,  1.86030000d-03 /
      data f( 7,17),f( 7,18) /  2.24340000d-03 ,  1.60210000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -2.26000000d-05 /
      data f( 8, 3),f( 8, 4) / -3.08000000d-05 , -2.61000000d-05 /
      data f( 8, 5),f( 8, 6) / -1.94000000d-05 , -1.99000000d-05 /
      data f( 8, 7),f( 8, 8) / -2.82000000d-05 , -3.29000000d-05 /
      data f( 8, 9),f( 8,10) / -3.02000000d-05 , -2.39000000d-05 /
      data f( 8,11),f( 8,12) / -1.80000000d-06 ,  3.39000000d-05 /
      data f( 8,13),f( 8,14) /  1.17100000d-04 ,  2.86300000d-04 /
      data f( 8,15),f( 8,16) /  5.92300000d-04 ,  9.71300000d-04 /
      data f( 8,17),f( 8,18) /  1.12850000d-03 ,  7.87800000d-04 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -2.49000000d-05 /
      data f( 9, 3),f( 9, 4) / -4.14000000d-05 , -4.81000000d-05 /
      data f( 9, 5),f( 9, 6) / -4.92000000d-05 , -5.24000000d-05 /
      data f( 9, 7),f( 9, 8) / -5.45000000d-05 , -5.13000000d-05 /
      data f( 9, 9),f( 9,10) / -4.16000000d-05 , -3.69000000d-05 /
      data f( 9,11),f( 9,12) / -1.49000000d-05 ,  7.40000000d-06 /
      data f( 9,13),f( 9,14) /  5.00000000d-05 ,  1.37500000d-04 /
      data f( 9,15),f( 9,16) /  2.89200000d-04 ,  4.70200000d-04 /
      data f( 9,17),f( 9,18) /  5.40100000d-04 ,  3.71100000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -9.90000000d-06 /
      data f(10, 3),f(10, 4) / -2.03000000d-05 , -2.94000000d-05 /
      data f(10, 5),f(10, 6) / -3.53000000d-05 , -3.79000000d-05 /
      data f(10, 7),f(10, 8) / -3.75000000d-05 , -3.33000000d-05 /
      data f(10, 9),f(10,10) / -2.55000000d-05 , -2.01000000d-05 /
      data f(10,11),f(10,12) / -1.25000000d-05 , -7.40000000d-06 /
      data f(10,13),f(10,14) /  9.00000000d-07 ,  1.64000000d-05 /
      data f(10,15),f(10,16) /  4.55000000d-05 ,  8.28000000d-05 /
      data f(10,17),f(10,18) /  1.02700000d-04 ,  7.33000000d-05 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -3.00000000d-06 /
      data f(11, 3),f(11, 4) / -7.50000000d-06 , -1.30000000d-05 /
      data f(11, 5),f(11, 6) / -1.71000000d-05 , -1.84000000d-05 /
      data f(11, 7),f(11, 8) / -1.70000000d-05 , -1.39000000d-05 /
      data f(11, 9),f(11,10) / -9.90000000d-06 , -6.20000000d-06 /
      data f(11,11),f(11,12) / -5.30000000d-06 , -5.50000000d-06 /
      data f(11,13),f(11,14) / -7.50000000d-06 , -4.60000000d-06 /
      data f(11,15),f(11,16) / -8.00000000d-07 ,  4.90000000d-06 /
      data f(11,17),f(11,18) /  1.01000000d-05 ,  9.00000000d-06 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 , -1.90000000d-06 /
      data f(12, 3),f(12, 4) / -4.40000000d-06 , -7.10000000d-06 /
      data f(12, 5),f(12, 6) / -8.90000000d-06 , -9.00000000d-06 /
      data f(12, 7),f(12, 8) / -7.50000000d-06 , -5.60000000d-06 /
      data f(12, 9),f(12,10) / -3.60000000d-06 , -2.10000000d-06 /
      data f(12,11),f(12,12) / -1.80000000d-06 , -2.30000000d-06 /
      data f(12,13),f(12,14) / -4.00000000d-06 , -3.80000000d-06 /
      data f(12,15),f(12,16) / -2.70000000d-06 , -1.40000000d-06 /
      data f(12,17),f(12,18) / -4.00000000d-07 ,  0.00000000d+00 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 , -1.20000000d-06 /
      data f(13, 3),f(13, 4) / -2.50000000d-06 , -3.70000000d-06 /
      data f(13, 5),f(13, 6) / -4.50000000d-06 , -4.30000000d-06 /
      data f(13, 7),f(13, 8) / -3.40000000d-06 , -2.20000000d-06 /
      data f(13, 9),f(13,10) / -1.10000000d-06 , -6.00000000d-07 /
      data f(13,11),f(13,12) / -6.00000000d-07 , -1.00000000d-06 /
      data f(13,13),f(13,14) / -1.30000000d-06 , -1.10000000d-06 /
      data f(13,15),f(13,16) / -3.00000000d-07 ,  5.00000000d-07 /
      data f(13,17),f(13,18) /  7.00000000d-07 ,  5.00000000d-07 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 3.67677378d-07/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.22423313d-02,-6.78635476d-06/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 2.23068856d-02,-1.49402584d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.06710382d-02,-1.95946118d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 3.67333975d-02,-2.37052943d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 3.02522020d-02,-3.66421096d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.53439191d-02, 1.28381382d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 6.77796682d-03, 6.59658330d-07/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 4.71129766d-03, 1.53228521d-07/
      data fpp( 1,10,1),fpp( 1,10,2)/ 5.64016753d-03, 1.79342759d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 3.40786883d-03, 1.16390611d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 4.94175857d-03, 7.90632788d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/ 9.33432194d-03, 6.11476273d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-3.25506689d-03, 4.30916277d-06/
      data fpp( 1,15,1),fpp( 1,15,2)/-4.54524608d-02, 7.15077216d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 3.11642228d-02,-2.28420049d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 5.54834329d-02, 7.13224748d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.55265706d-01,-1.77325850d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-2.69291075d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 1.43518725d-06/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 9.65512320d-03,-4.41637450d-06/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 1.75125859d-02,-1.07336892d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.36414950d-02,-1.54208685d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 2.74584908d-02,-1.38928367d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 2.22999531d-02,-2.69784826d-07/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.16343046d-02, 6.27197597d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 5.66685207d-03,-1.16611904d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 4.06111897d-03,-5.57499815d-07/
      data fpp( 2,10,1),fpp( 2,10,2)/ 4.25652208d-03, 3.01811830d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 2.64947663d-03, 7.43302663d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.75098287d-03, 5.82377519d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/ 7.77621327d-03, 4.35578726d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 3.21413378d-03, 3.44867344d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/-1.92700070d-02, 1.27391900d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 2.68211258d-02,-1.50927494d-04/
      data fpp( 2,17,1),fpp( 2,17,2)/ 3.34961342d-02, 8.76607869d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-9.15219452d-02,-2.09999653d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-3.54446173d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 1.75586773d-06/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 6.51217594d-03,-2.71573545d-06/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 1.16377709d-02,-7.40492591d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 1.47129817d-02,-1.13325609d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.49826395d-02,-8.60283049d-06/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 1.12429855d-02, 3.17882850d-07/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 6.36886258d-03, 2.81929909d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 4.19962488d-03,-1.50307920d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 3.34422645d-03,-2.32982290d-07/
      data fpp( 3,10,1),fpp( 3,10,2)/ 2.23374416d-03, 2.59100836d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 2.02922463d-03, 5.01294885d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.19930996d-03, 4.63919625d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/ 6.66082499d-03, 2.73282661d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.60785318d-02, 4.01757392d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 7.61748895d-03, 2.89477707d-06/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.56127399d-03,-8.64348475d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-2.27579696d-02, 1.57846128d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 2.38484869d-02,-1.53181604d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-2.44666198d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 1.63550331d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 3.76599806d-03,-1.14300663d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 6.65237299d-03,-3.92747681d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 8.04573086d-03,-6.69108615d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 7.60554120d-03,-4.64217860d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 5.65341282d-03,-1.86199465d-07/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 3.98092169d-03, 7.00976458d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 3.26001567d-03,-1.07570637d-06/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.50516586d-03,-1.01509974d-08/
      data fpp( 4,10,1),fpp( 4,10,2)/ 1.46317141d-03, 2.13031035d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/ 6.09600127d-04, 2.73290958d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.08311543d-04, 3.34805133d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.41977451d-03, 1.17508851d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.27254716d-02, 3.12664083d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 4.54517082d-02, 5.82148176d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/ 3.77516695d-02,-3.10503353d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.42224760d-02,-3.69901405d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 1.40596739d-02,-9.19671027d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-1.32933449d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 1.09082894d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.52310063d-03,-5.06578891d-08/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 2.94174387d-03,-1.05019739d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 3.65067221d-03,-2.54655256d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.41108448d-03,-2.00559238d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.87208764d-03,-4.35077929d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 2.50853304d-03,-2.56095907d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.93377493d-03,-1.06538442d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.84493337d-03,-3.25750326d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/ 2.17440500d-03, 2.21953975d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 2.35414481d-03,-1.31640866d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/ 2.98421708d-03, 7.34209489d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 2.00982656d-03,-8.61797088d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 3.91857197d-03, 5.19937886d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/-1.90207596d-02,-3.87731837d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/-4.04506825d-04, 7.65694615d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 1.10252587d-02,-3.30686009d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 9.40475126d-03,-4.31665426d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-6.20152287d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 9.97812484d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 7.95199412d-04, 3.64375031d-07/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.19825152d-03,-1.39312610d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.43798029d-03,-1.09712459d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.64612088d-03,-7.70189027d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.53263660d-03, 3.78806992d-08/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 1.10734613d-03,-1.02533377d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 1.27608461d-03, 1.14145438d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 3.70300641d-04, 1.43516254d-07/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.11279141d-03, 9.12480606d-07/
      data fpp( 6,11,1),fpp( 6,11,2)/-2.27417938d-03,-7.87438677d-07/
      data fpp( 6,12,1),fpp( 6,12,2)/-2.62277988d-03, 1.83527410d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.65028074d-03, 1.22834226d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/-4.32559413d-04, 5.79135685d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 1.84225301d-02, 8.30023035d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.10167578d-02,-4.54027824d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/ 5.25008923d-03,-2.45371174d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 5.30932101d-03,-1.94232522d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-2.14838739d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 2.95181759d-07/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 2.04101720d-04, 1.73636482d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 6.14850054d-04, 1.38272314d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 9.12606618d-04,-1.98725739d-07/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 8.74032005d-04,-2.25369356d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 7.58165949d-04,-3.33796835d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 8.30882429d-04, 6.65566946d-08/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 4.54686641d-04, 7.95700562d-08/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 7.15464063d-04, 5.31630806d-08/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.31916064d-03, 2.83777622d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.90177269d-03, 3.77726433d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.78050244d-03, 5.99316646d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.42569642d-03, 1.97100698d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/ 8.57265687d-04, 8.58655420d-07/
      data fpp( 7,15,1),fpp( 7,15,2)/-2.57016091d-03, 6.50437134d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 2.47875661d-04,-5.63614077d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 3.65038438d-03,-1.16558083d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/ 2.78996468d-03,-9.20462621d-06/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-9.17368689d-06/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 1.62933480d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.25993709d-04, 1.36133041d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 1.87148264d-04, 1.56534357d-07/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 2.34793236d-04, 1.17295302d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 2.93751102d-04,-8.34524780d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 3.01899603d-04,-1.09919618d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 2.17924151d-04, 5.51309509d-08/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 2.43568827d-04, 1.05395815d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 9.66431085d-05,-3.27142093d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/-1.67851168d-04, 2.41461023d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/-3.40911380d-04, 1.48701184d-08/
      data fpp( 8,12,1),fpp( 8,12,2)/-2.60829876d-04, 5.15058504d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.11704937d-04, 7.74895867d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.53096667d-04, 1.54535803d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/ 1.24931353d-03, 1.25167202d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.06533956d-03,-2.17204611d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 1.97877323d-03,-5.87148759d-06/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.49482026d-03,-4.21600355d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-4.09049823d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 7.16080813d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 2.39234450d-05, 8.07838374d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 6.21568899d-05, 1.09256569d-07/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 9.46204388d-05, 7.01898858d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.03763587d-04,-5.40161126d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.15035640d-04, 1.98745644d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.09420967d-04, 4.05178551d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 7.34380497d-05, 1.36054015d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 6.67635033d-05,-1.94733916d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/ 1.25044028d-04, 3.42881650d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.21872830d-04,-1.38792684d-07/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.22017065d-04, 2.30289085d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.85123326d-04, 4.35636344d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/ 2.85547646d-04, 7.21165539d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/ 4.57706803d-04, 5.31701499d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 8.00366085d-04,-1.08997154d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/ 1.07052268d-03,-2.83781536d-06/
      data fpp( 9,18,1),fpp( 9,18,2)/ 7.73154283d-04,-1.89276704d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-1.71711648d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00,-2.32240170d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.71671895d-05,-4.55196597d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/-2.62448015d-05, 1.14318809d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/-2.50579343d-05, 3.68244424d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.71663120d-05, 3.32703493d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.90567205d-05, 2.80941602d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-1.96249771d-05, 3.43530099d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/-1.32985626d-05, 6.24938001d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.52120643d-05,-6.83282102d-08/
      data fpp(10,10,1),fpp(10,10,2)/-3.44065005d-05, 6.68190406d-08/
      data fpp(10,11,1),fpp(10,11,2)/-2.35627995d-05,-6.69479521d-08/
      data fpp(10,12,1),fpp(10,12,2)/-6.43625682d-06, 5.09727680d-08/
      data fpp(10,13,1),fpp(10,13,2)/ 1.10824902d-05, 5.50568803d-08/
      data fpp(10,14,1),fpp(10,14,2)/ 7.58087285d-05, 1.60799711d-07/
      data fpp(10,15,1),fpp(10,15,2)/ 1.77222827d-04, 1.17744276d-07/
      data fpp(10,16,1),fpp(10,16,2)/ 2.55031963d-04,-1.39776816d-07/
      data fpp(10,17,1),fpp(10,17,2)/ 2.35445339d-04,-6.02637011d-07/
      data fpp(10,18,1),fpp(10,18,2)/ 1.46727023d-04,-4.07675140d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-4.00662430d-07/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-2.08614175d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-3.85468697d-06,-1.32771650d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-6.97768381d-06,-1.60299225d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/-8.18870149d-06, 1.73968551d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/-9.29833894d-06, 3.04425021d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-8.80875777d-06, 2.88331366d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-9.92105864d-06, 1.62249516d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/-1.18437993d-05, 8.26705693d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-8.91524624d-06, 4.70682067d-09/
      data fpp(11,10,1),fpp(11,10,2)/-4.81802618d-06,-4.50943396d-08/
      data fpp(11,11,1),fpp(11,11,2)/ 1.17836829d-06, 7.67053782d-09/
      data fpp(11,12,1),fpp(11,12,2)/ 3.92796232d-06,-5.15878116d-08/
      data fpp(11,13,1),fpp(11,13,2)/ 1.47467131d-05, 9.06807088d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 1.18174399d-05,-1.71350234d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 1.78018875d-05, 3.18593847d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 3.65060623d-05, 3.69748451d-09/
      data fpp(11,17,1),fpp(11,17,2)/ 5.64959635d-05,-7.66493228d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 4.09376255d-05,-7.51001935d-08/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-9.69499033d-08/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-1.03074829d-08/
      data fpp(12, 2,1),fpp(12, 2,2)/-2.21406261d-06,-5.38503412d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/-4.04446324d-06,-4.15238059d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/-5.18725970d-06, 9.99455648d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/-5.64033221d-06, 1.81741547d-08/
      data fpp(12, 6,1),fpp(12, 6,2)/-6.30824845d-06, 1.93088248d-08/
      data fpp(12, 7,1),fpp(12, 7,2)/-6.69078827d-06, 5.90546136d-10/
      data fpp(12, 8,1),fpp(12, 8,2)/-5.92624014d-06, 2.32899066d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/-4.92695075d-06,-3.90650877d-09/
      data fpp(12,10,1),fpp(12,10,2)/-5.12139476d-06,-1.67029556d-08/
      data fpp(12,11,1),fpp(12,11,2)/-3.35067366d-06,-1.28166897d-09/
      data fpp(12,12,1),fpp(12,12,2)/-1.47559246d-06,-2.61703686d-08/
      data fpp(12,13,1),fpp(12,13,2)/ 1.33065738d-06, 3.39631432d-08/
      data fpp(12,14,1),fpp(12,14,2)/ 7.72151203d-06, 4.31779583d-09/
      data fpp(12,15,1),fpp(12,15,2)/ 1.79696225d-05, 2.76567349d-09/
      data fpp(12,16,1),fpp(12,16,2)/ 2.85437875d-05,-3.38048980d-09/
      data fpp(12,17,1),fpp(12,17,2)/ 3.11708073d-05,-7.24371428d-09/
      data fpp(12,18,1),fpp(12,18,2)/ 2.13224749d-05,-3.64465306d-09/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-2.17767347d-09/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-2.98229723d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/ 6.31953130d-06,-1.03540554d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/ 9.17223162d-06, 1.12391937d-09/
      data fpp(13, 4,1),fpp(13, 4,2)/ 7.05612985d-06, 2.53972804d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/ 3.57016611d-06, 1.27171685d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/ 2.17912422d-06, 6.59159809d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 2.68289414d-06, 2.91643919d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/ 3.90062007d-06,-2.57354830d-10/
      data fpp(13, 9,1),fpp(13, 9,2)/ 4.08847538d-06,-7.88701987d-09/
      data fpp(13,10,1),fpp(13,10,2)/ 7.72319738d-06,-4.19456571d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 7.62836829d-07,-5.33471730d-09/
      data fpp(13,12,1),fpp(13,12,2)/-5.18720377d-06, 1.53343491d-09/
      data fpp(13,13,1),fpp(13,13,2)/-1.78153287d-05, 5.20097768d-09/
      data fpp(13,14,1),fpp(13,14,2)/-2.74232560d-05, 7.66265439d-09/
      data fpp(13,15,1),fpp(13,15,2)/-5.35098112d-05, 1.48404763d-10/
      data fpp(13,16,1),fpp(13,16,2)/-8.21343938d-05,-8.25627344d-09/
      data fpp(13,17,1),fpp(13,17,2)/-8.86104037d-05,-3.12331100d-09/
      data fpp(13,18,1),fpp(13,18,2)/-5.66862374d-05,-3.25048257d-09/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-1.87475871d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-2.74873255d-05,-1.99061561d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/-2.35546622d-05, 1.21006955d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/-1.62955718d-04,-1.12891109d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.08894909d-04, 5.78513050d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.96568604d-05,-7.67464039d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 6.76586400d-05, 3.20831508d-05/
      data fpppp( 1,13),fpppp( 1,14)/-2.44708254d-05,-9.53116981d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 2.06045845d-03,-1.59872160d-04/
      data fpppp( 1,17),fpppp( 1,18)/-4.55881822d-03, 4.29104408d-03/
      data fpppp( 1,19) /             9.35553257d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-1.92985822d-05,-1.67147929d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/-2.17018788d-05,-1.90901867d-07/
      data fpppp( 2, 5),fpppp( 2, 6)/-1.16249318d-04,-7.33438277d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/ 7.91979750d-05, 3.84436905d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 2.87304269d-05,-4.52972259d-05/
      data fpppp( 2,11),fpppp( 2,12)/ 4.43115644d-05, 3.05640689d-05/
      data fpppp( 2,13),fpppp( 2,14)/ 8.85560997d-06,-5.81225102d-04/
      data fpppp( 2,15),fpppp( 2,16)/ 1.24072112d-03,-2.67142953d-04/
      data fpppp( 2,17),fpppp( 2,18)/-2.53711678d-03, 2.51402479d-03/
      data fpppp( 2,19) /             5.47341910d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-7.34275534d-06,-1.36358532d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-2.13086890d-05,-2.41524417d-05/
      data fpppp( 3, 5),fpppp( 3, 6)/-5.04147278d-05,-1.47473499d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/ 4.13359878d-05, 1.16965141d-05/
      data fpppp( 3, 9),fpppp( 3,10)/-9.29168893d-06, 1.01652109d-05/
      data fpppp( 3,11),fpppp( 3,12)/ 2.29886108d-05,-1.96433628d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 1.93070622d-04,-3.95267621d-04/
      data fpppp( 3,15),fpppp( 3,16)/ 3.15274884d-04,-7.21542245d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 1.47511238d-03,-9.23365255d-04/
      data fpppp( 3,19) /            -2.00894796d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.75076926d-06,-8.61927719d-06/
      data fpppp( 4, 3),fpppp( 4, 4)/-1.55495102d-05,-1.87637054d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.94085200d-05, 5.68146271d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/ 1.34609033d-05,-2.42996899d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-5.77765438d-06, 8.31190755d-06/
      data fpppp( 4,11),fpppp( 4,12)/-1.61645853d-05, 7.74833953d-05/
      data fpppp( 4,13),fpppp( 4,14)/-6.50039026d-05, 5.42186259d-04/
      data fpppp( 4,15),fpppp( 4,16)/-6.98508758d-04,-1.73727748d-04/
      data fpppp( 4,17),fpppp( 4,18)/ 4.43670457d-04,-1.98970588d-04/
      data fpppp( 4,19) /            -4.81600419d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 4.95050653d-06,-9.29584459d-07/
      data fpppp( 5, 3),fpppp( 5, 4)/-7.49961203d-06,-1.16548616d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.79190590d-06, 4.85793874d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/-6.11331467d-06, 6.92310908d-06/
      data fpppp( 5, 9),fpppp( 5,10)/ 7.57587175d-06,-1.21278050d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 3.19514395d-05,-8.86580051d-05/
      data fpppp( 5,13),fpppp( 5,14)/ 2.26412813d-04,-6.44005091d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 8.58722935d-04,-2.97551590d-04/
      data fpppp( 5,17),fpppp( 5,18)/-9.97058073d-05,-8.66415589d-05/
      data fpppp( 5,19) /            -2.07825863d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/-6.21146706d-06,-3.91796458d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.64551301d-06, 7.00616744d-07/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.05224525d-06,-7.78912747d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.55003835d-05,-1.85706698d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-5.68905082d-06, 6.68838789d-06/
      data fpppp( 6,11),fpppp( 6,12)/-1.76225541d-06, 4.91278812d-05/
      data fpppp( 6,13),fpppp( 6,14)/-1.15483291d-04, 4.27518614d-04/
      data fpppp( 6,15),fpppp( 6,16)/-5.36349074d-04, 1.42225969d-04/
      data fpppp( 6,17),fpppp( 6,18)/ 6.57914252d-05,-5.58376485d-05/
      data fpppp( 6,19) /            -1.64553999d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 5.26461935d-06, 2.06292749d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-1.11753239d-06,-4.37230420d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/-1.57312142d-06, 6.02730328d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.12211395d-05, 1.19225187d-05/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.74945744d-06, 1.65480124d-06/
      data fpppp( 7,11),fpppp( 7,12)/-9.63373471d-06,-5.35279998d-06/
      data fpppp( 7,13),fpppp( 7,14)/ 1.70327885d-05,-7.55958368d-05/
      data fpppp( 7,15),fpppp( 7,16)/ 1.13810807d-04,-4.91960002d-06/
      data fpppp( 7,17),fpppp( 7,18)/-5.90640778d-05,-1.45997945d-05/
      data fpppp( 7,19) /             1.69055705d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/-1.15783163d-06,-6.56102752d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/-1.08106580d-07, 2.77954081d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/-3.24936085d-07,-2.02677166d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 2.90458556d-06,-3.01436288d-06/
      data fpppp( 8, 9),fpppp( 8,10)/-1.20135775d-06, 7.65680408d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 3.62468005d-06,-7.58977046d-08/
      data fpppp( 8,13),fpppp( 8,14)/ 8.21516926d-07, 9.73042984d-06/
      data fpppp( 8,15),fpppp( 8,16)/-1.85832086d-06,-1.31085959d-05/
      data fpppp( 8,17),fpppp( 8,18)/ 1.37162433d-07,-1.12832525d-05/
      data fpppp( 8,19) /            -1.56561894d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 3.48856345d-07, 1.33185217d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-2.29972255d-08,-3.87390066d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.73333444d-07,-1.78209441d-07/
      data fpppp( 9, 7),fpppp( 9, 8)/-4.73699186d-07, 2.50911475d-07/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.22855556d-06,-1.26782945d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 1.55658824d-07, 8.44120166d-07/
      data fpppp( 9,13),fpppp( 9,14)/ 2.45582067d-07, 4.12635099d-07/
      data fpppp( 9,15),fpppp( 9,16)/ 2.40796775d-06, 1.85501394d-07/
      data fpppp( 9,17),fpppp( 9,18)/-7.50013446d-06,-4.23646326d-06/
      data fpppp( 9,19) /            -4.10116550d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 5.92589733d-08, 8.06713461d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 1.03430294d-07, 1.21476230d-07/
      data fpppp(10, 5),fpppp(10, 6)/-1.87049904d-07, 3.98015389d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 1.07172856d-07,-5.48126867d-08/
      data fpppp(10, 9),fpppp(10,10)/-3.82317084d-07, 5.47224950d-07/
      data fpppp(10,11),fpppp(10,12)/-4.29448128d-09,-1.53076522d-07/
      data fpppp(10,13),fpppp(10,14)/ 6.40132830d-07, 4.24994680d-07/
      data fpppp(10,15),fpppp(10,16)/-1.38839917d-07,-1.28593279d-06/
      data fpppp(10,17),fpppp(10,18)/-5.61174518d-07,-6.17270614d-07/
      data fpppp(10,19) /            -4.50265460d-07 /
      data fpppp(11, 1),fpppp(11, 2)/-2.72405077d-09, 3.79302450d-09/
      data fpppp(11, 3),fpppp(11, 4)/ 3.14533603d-08,-1.48877157d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 3.41803160d-08,-2.58804309d-08/
      data fpppp(11, 7),fpppp(11, 8)/-2.67715153d-08, 8.43401044d-08/
      data fpppp(11, 9),fpppp(11,10)/-1.95112770d-08, 6.38250228d-08/
      data fpppp(11,11),fpppp(11,12)/-1.21838349d-07, 2.28720347d-07/
      data fpppp(11,13),fpppp(11,14)/-3.08893634d-07, 1.81972746d-07/
      data fpppp(11,15),fpppp(11,16)/ 1.15825903d-07, 1.17907270d-07/
      data fpppp(11,17),fpppp(11,18)/-5.10311400d-07,-2.09556016d-07/
      data fpppp(11,19) /            -1.74221790d-07 /
      data fpppp(12, 1),fpppp(12, 2)/ 7.85240073d-10, 3.86053526d-09/
      data fpppp(12, 3),fpppp(12, 4)/ 6.79233738d-09, 1.02263653d-08/
      data fpppp(12, 5),fpppp(12, 6)/-6.31436112d-09, 2.14045572d-09/
      data fpppp(12, 7),fpppp(12, 8)/ 1.48751228d-08, 7.18433048d-09/
      data fpppp(12, 9),fpppp(12,10)/-2.95279698d-08, 3.93035448d-08/
      data fpppp(12,11),fpppp(12,12)/-9.77630238d-09, 6.06327013d-09/
      data fpppp(12,13),fpppp(12,14)/ 4.13933408d-08, 4.34396548d-08/
      data fpppp(12,15),fpppp(12,16)/ 1.62833891d-08,-8.90099365d-08/
      data fpppp(12,17),fpppp(12,18)/-1.37072359d-07,-1.11221757d-07/
      data fpppp(12,19) /            -1.06489162d-07 /
      data fpppp(13, 1),fpppp(13, 2)/-2.16677200d-08,-3.06300678d-08/
      data fpppp(13, 3),fpppp(13, 4)/-6.38218681d-08,-1.22105847d-08/
      data fpppp(13, 5),fpppp(13, 6)/ 3.04724882d-08, 1.60159436d-08/
      data fpppp(13, 7),fpppp(13, 8)/ 1.91524452d-08,-4.97883630d-08/
      data fpppp(13, 9),fpppp(13,10)/ 1.18208770d-07,-2.16234713d-07/
      data fpppp(13,11),fpppp(13,12)/ 1.11025130d-07,-1.67246608d-07/
      data fpppp(13,13),fpppp(13,14)/ 1.57276243d-07,-2.80646507d-07/
      data fpppp(13,15),fpppp(13,16)/-2.34078880d-08, 2.21996422d-07/
      data fpppp(13,17),fpppp(13,18)/ 4.64336557d-07, 2.24667913d-07/
      data fpppp(13,19) /             1.22716063d-07 /
 
 
      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine c2_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  1.69790000d-03 /
      data f( 1, 3),f( 1, 4) /  5.32640000d-03 ,  8.16340000d-03 /
      data f( 1, 5),f( 1, 6) /  8.51060000d-03 ,  6.16580000d-03 /
      data f( 1, 7),f( 1, 8) /  2.66720000d-03 ,  2.97400000d-04 /
      data f( 1, 9),f( 1,10) /  1.07290000d-03 ,  6.16710000d-03 /
      data f( 1,11),f( 1,12) /  1.50358000d-02 ,  2.68592000d-02 /
      data f( 1,13),f( 1,14) /  4.21482000d-02 ,  6.24326000d-02 /
      data f( 1,15),f( 1,16) /  7.98035000d-02 ,  7.70153000d-02 /
      data f( 1,17),f( 1,18) /  4.09759000d-02 ,  7.08720000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  1.26470000d-03 /
      data f( 2, 3),f( 2, 4) /  3.99950000d-03 ,  6.13750000d-03 /
      data f( 2, 5),f( 2, 6) /  6.32270000d-03 ,  4.51030000d-03 /
      data f( 2, 7),f( 2, 8) /  1.88390000d-03 ,  1.07900000d-04 /
      data f( 2, 9),f( 2,10) /  7.65900000d-04 ,  4.91590000d-03 /
      data f( 2,11),f( 2,12) /  1.24375000d-02 ,  2.25667000d-02 /
      data f( 2,13),f( 2,14) /  3.53902000d-02 ,  5.22787000d-02 /
      data f( 2,15),f( 2,16) /  6.62825000d-02 ,  6.29564000d-02 /
      data f( 2,17),f( 2,18) /  3.68843000d-02 ,  9.21190000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  9.39100000d-04 /
      data f( 3, 3),f( 3, 4) /  2.98790000d-03 ,  4.59630000d-03 /
      data f( 3, 5),f( 3, 6) /  4.70570000d-03 ,  3.31270000d-03 /
      data f( 3, 7),f( 3, 8) /  1.32510000d-03 , -1.16000000d-05 /
      data f( 3, 9),f( 3,10) /  5.40200000d-04 ,  3.91230000d-03 /
      data f( 3,11),f( 3,12) /  1.02012000d-02 ,  1.88006000d-02 /
      data f( 3,13),f( 3,14) /  2.94307000d-02 ,  4.27591000d-02 /
      data f( 3,15),f( 3,16) /  5.41668000d-02 ,  5.22793000d-02 /
      data f( 3,17),f( 3,18) /  3.31468000d-02 ,  9.84260000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  5.95100000d-04 /
      data f( 4, 3),f( 4, 4) /  1.90720000d-03 ,  2.94700000d-03 /
      data f( 4, 5),f( 4, 6) /  3.00250000d-03 ,  2.07520000d-03 /
      data f( 4, 7),f( 4, 8) /  7.62100000d-04 , -1.10800000d-04 /
      data f( 4, 9),f( 4,10) /  3.02700000d-04 ,  2.73460000d-03 /
      data f( 4,11),f( 4,12) /  7.44230000d-03 ,  1.40685000d-02 /
      data f( 4,13),f( 4,14) /  2.20605000d-02 ,  3.11189000d-02 /
      data f( 4,15),f( 4,16) /  3.88309000d-02 ,  3.80658000d-02 /
      data f( 4,17),f( 4,18) /  2.56476000d-02 ,  8.28490000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  2.70000000d-04 /
      data f( 5, 3),f( 5, 4) /  8.73600000d-04 ,  1.35710000d-03 /
      data f( 5, 5),f( 5, 6) /  1.37350000d-03 ,  9.13600000d-04 /
      data f( 5, 7),f( 5, 8) /  2.64800000d-04 , -1.61200000d-04 /
      data f( 5, 9),f( 5,10) /  7.97000000d-05 ,  1.43090000d-03 /
      data f( 5,11),f( 5,12) /  4.18700000d-03 ,  8.26010000d-03 /
      data f( 5,13),f( 5,14) /  1.32301000d-02 ,  1.82824000d-02 /
      data f( 5,15),f( 5,16) /  2.12521000d-02 ,  2.11980000d-02 /
      data f( 5,17),f( 5,18) /  1.49201000d-02 ,  4.98920000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.22000000d-04 /
      data f( 6, 3),f( 6, 4) /  3.89100000d-04 ,  5.89600000d-04 /
      data f( 6, 5),f( 6, 6) /  5.90100000d-04 ,  3.73900000d-04 /
      data f( 6, 7),f( 6, 8) /  6.19000000d-05 , -1.45800000d-04 /
      data f( 6, 9),f( 6,10) / -2.86000000d-05 ,  6.75400000d-04 /
      data f( 6,11),f( 6,12) /  2.05410000d-03 ,  4.34490000d-03 /
      data f( 6,13),f( 6,14) /  7.28410000d-03 ,  1.03106000d-02 /
      data f( 6,15),f( 6,16) /  1.24153000d-02 ,  1.21381000d-02 /
      data f( 6,17),f( 6,18) /  8.39700000d-03 ,  2.93260000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  5.07000000d-05 /
      data f( 7, 3),f( 7, 4) /  1.63600000d-04 ,  2.52300000d-04 /
      data f( 7, 5),f( 7, 6) /  2.46300000d-04 ,  1.37000000d-04 /
      data f( 7, 7),f( 7, 8) / -1.67000000d-05 , -1.23100000d-04 /
      data f( 7, 9),f( 7,10) / -3.47000000d-05 ,  3.38700000d-04 /
      data f( 7,11),f( 7,12) /  1.12390000d-03 ,  2.37830000d-03 /
      data f( 7,13),f( 7,14) /  4.00940000d-03 ,  5.69390000d-03 /
      data f( 7,15),f( 7,16) /  6.78400000d-03 ,  6.46750000d-03 /
      data f( 7,17),f( 7,18) /  4.28150000d-03 ,  1.36900000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  2.10000000d-05 /
      data f( 8, 3),f( 8, 4) /  6.76000000d-05 ,  1.02200000d-04 /
      data f( 8, 5),f( 8, 6) /  9.32000000d-05 ,  3.74000000d-05 /
      data f( 8, 7),f( 8, 8) / -4.20000000d-05 , -8.56000000d-05 /
      data f( 8, 9),f( 8,10) / -4.00000000d-05 ,  1.50100000d-04 /
      data f( 8,11),f( 8,12) /  5.41400000d-04 ,  1.19250000d-03 /
      data f( 8,13),f( 8,14) /  2.03320000d-03 ,  2.91990000d-03 /
      data f( 8,15),f( 8,16) /  3.48780000d-03 ,  3.27560000d-03 /
      data f( 8,17),f( 8,18) /  2.12620000d-03 ,  6.70400000d-04 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  8.50000000d-06 /
      data f( 9, 3),f( 9, 4) /  2.68000000d-05 ,  3.83000000d-05 /
      data f( 9, 5),f( 9, 6) /  2.82000000d-05 , -3.40000000d-06 /
      data f( 9, 7),f( 9, 8) / -3.84000000d-05 , -5.68000000d-05 /
      data f( 9, 9),f( 9,10) / -3.26000000d-05 ,  6.55000000d-05 /
      data f( 9,11),f( 9,12) /  2.50300000d-04 ,  5.62200000d-04 /
      data f( 9,13),f( 9,14) /  9.85300000d-04 ,  1.42900000d-03 /
      data f( 9,15),f( 9,16) /  1.70430000d-03 ,  1.58540000d-03 /
      data f( 9,17),f( 9,18) /  1.01740000d-03 ,  3.18100000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.30000000d-06 /
      data f(10, 3),f(10, 4) /  3.50000000d-06 ,  3.00000000d-06 /
      data f(10, 5),f(10, 6) / -2.50000000d-06 , -1.18000000d-05 /
      data f(10, 7),f(10, 8) / -2.08000000d-05 , -2.33000000d-05 /
      data f(10, 9),f(10,10) / -1.54000000d-05 ,  8.60000000d-06 /
      data f(10,11),f(10,12) /  4.90000000d-05 ,  1.17200000d-04 /
      data f(10,13),f(10,14) /  2.07700000d-04 ,  3.00200000d-04 /
      data f(10,15),f(10,16) /  3.54300000d-04 ,  3.25500000d-04 /
      data f(10,17),f(10,18) /  2.07400000d-04 ,  6.46000000d-05 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(11, 3),f(11, 4) /  3.00000000d-07 , -7.00000000d-07 /
      data f(11, 5),f(11, 6) / -2.90000000d-06 , -5.60000000d-06 /
      data f(11, 7),f(11, 8) / -7.70000000d-06 , -7.90000000d-06 /
      data f(11, 9),f(11,10) / -5.60000000d-06 , -4.00000000d-07 /
      data f(11,11),f(11,12) /  9.10000000d-06 ,  2.30000000d-05 /
      data f(11,13),f(11,14) /  4.24000000d-05 ,  5.63000000d-05 /
      data f(11,15),f(11,16) /  6.17000000d-05 ,  5.35000000d-05 /
      data f(11,17),f(11,18) /  3.35000000d-05 ,  1.05000000d-05 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) / -4.00000000d-07 , -1.00000000d-06 /
      data f(12, 5),f(12, 6) / -1.80000000d-06 , -2.30000000d-06 /
      data f(12, 7),f(12, 8) / -2.40000000d-06 , -2.10000000d-06 /
      data f(12, 9),f(12,10) / -1.40000000d-06 , -3.00000000d-07 /
      data f(12,11),f(12,12) /  1.50000000d-06 ,  4.20000000d-06 /
      data f(12,13),f(12,14) /  8.40000000d-06 ,  1.08000000d-05 /
      data f(12,15),f(12,16) /  1.10000000d-05 ,  8.20000000d-06 /
      data f(12,17),f(12,18) /  4.40000000d-06 ,  1.20000000d-06 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(13, 3),f(13, 4) / -4.00000000d-07 , -8.00000000d-07 /
      data f(13, 5),f(13, 6) / -1.20000000d-06 , -1.30000000d-06 /
      data f(13, 7),f(13, 8) / -1.10000000d-06 , -7.00000000d-07 /
      data f(13, 9),f(13,10) / -3.00000000d-07 ,  0.00000000d+00 /
      data f(13,11),f(13,12) /  2.00000000d-07 ,  1.00000000d-07 /
      data f(13,13),f(13,14) /  1.00000000d-07 ,  2.00000000d-07 /
      data f(13,15),f(13,16) /  1.00000000d-07 , -1.00000000d-07 /
      data f(13,17),f(13,18) / -2.00000000d-07 , -1.00000000d-07 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 4.62232598d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 3.33498168d-03, 1.99134804d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 9.62414447d-03,-1.00411814d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 1.48874984d-02,-2.72387549d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 1.81614408d-02,-3.03917991d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.47924145d-02,-1.27140488d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 7.24044437d-03, 1.20199941d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.32188777d-03, 3.23620723d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.60210973d-03, 4.72497168d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/ 7.69301830d-03, 3.77610606d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/ 9.95891028d-03, 2.81760407d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 1.39491154d-02, 2.68167766d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.92836356d-02, 7.24928528d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.17565509d-04,-1.70641878d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 3.28576842d-02,-1.79046102d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/ 1.32584173d-01,-4.76297406d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 3.44377118d-02, 8.91637246d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-4.20880397d-02, 2.48684507d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 5.24188246d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 3.51882544d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 2.71089379d-03, 1.50644911d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 7.93035392d-03,-7.24021886d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 1.22075033d-02,-2.19116157d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 1.44520470d-02,-2.22813185d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.16070997d-02,-8.81911048d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 5.68775413d-03, 8.71776039d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 1.77651017d-03, 2.49720689d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 2.06785196d-03, 3.74339639d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/ 6.30510627d-03, 3.48120756d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/ 9.13389373d-03, 2.56137337d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/ 1.31806263d-02, 1.91889897d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.96776573d-02, 5.92883074d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.39970833d-02,-1.24422193d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/ 3.50487030d-02,-1.82601430d-04/
      data fpp( 2,16,1),fpp( 2,16,2)/ 9.22089402d-02,-2.96946060d-04/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.61252192d-02, 5.62566999d-06/
      data fpp( 2,18,1),fpp( 2,18,2)/-3.46739205d-02, 1.78425380d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 3.88302810d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 2.64869825d-05/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 1.96144316d-03, 1.13190350d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 5.94943984d-03,-5.18112242d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 8.98748857d-03,-1.70185453d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 9.66537116d-03,-1.66846963d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 7.46418690d-03,-6.38666934d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 3.68353913d-03, 6.55537369d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 1.07207153d-03, 1.92191746d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 1.32148241d-03, 2.98779279d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/ 4.22655664d-03, 3.04871136d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/ 7.80551480d-03, 2.31816176d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 1.22883795d-02, 1.54164161d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 2.17807351d-02, 3.69947179d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 3.89391014d-02,-1.49728785d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.77425037d-02,-1.46247567d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 5.85006652d-03,-2.11224446d-04/
      data fpp( 3,17,1),fpp( 3,17,2)/-4.58235887d-02,-4.35546491d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/-4.33162781d-02, 1.35141043d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 3.10686479d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 1.70113482d-05/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 1.28126029d-03, 7.27330353d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 3.99496458d-03,-3.08456237d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 6.07003593d-03,-1.12730541d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 6.30073144d-03,-1.08812214d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 4.64131055d-03,-4.17006049d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 2.27636681d-03, 4.41346332d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 5.78754769d-04, 1.29282072d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 9.53157328d-04, 2.10577078d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/ 3.55474037d-03, 2.39449616d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/ 7.59568817d-03, 1.97104459d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.13883176d-02, 1.23232549d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/ 1.88824449d-02, 1.29445344d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.68182731d-02,-1.17392645d-07/
      data fpp( 4,15,1),fpp( 4,15,2)/ 4.00025190d-02,-9.32589638d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/ 3.91704848d-02,-1.35472752d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.57984828d-02,-6.40360283d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 5.86873988d-04, 9.49468652d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 2.28916567d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 7.85839585d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 6.80701192d-04, 3.36520830d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 2.06804943d-03,-1.30322907d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.99779187d-03,-5.35829204d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.07043668d-03,-5.28960277d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.29089411d-03,-2.06129687d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.09030274d-03, 2.20079024d-06/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 2.63141819d-04, 6.62613589d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 3.05007106d-04, 1.13086662d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/ 1.90809685d-03, 1.47571993d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/ 3.23928898d-03, 1.39565365d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/ 6.06675592d-03, 8.43665450d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 9.38613523d-03, 6.11084545d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 1.63500653d-02,-2.79420363d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 4.08912370d-02,-1.92987002d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.48572087d-02,-7.62911628d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 1.94470082d-02,-4.89646485d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 7.32297009d-03, 5.29697567d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 1.33587622d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 3.54505908d-06/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 2.46334947d-04, 1.49688183d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 9.11237706d-04,-8.26586412d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.67639659d-03,-2.18653618d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.71192184d-03,-2.42726886d-06/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.12071302d-03,-1.10638839d-06/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 4.28022248d-04, 1.10482242d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-5.21220447d-05, 2.94509872d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 5.79614248d-04, 6.60878271d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 1.96967224d-03, 5.82777043d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 6.38475592d-03, 1.05621356d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/ 9.78145868d-03, 6.64968733d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/ 1.27986142d-02, 1.74311512d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/ 1.45342658d-02,-8.38414780d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 6.24053313d-03,-2.35145239d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/ 8.79028043d-03,-4.04717565d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 7.31908442d-03,-2.24324502d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.40354352d-04, 2.68035572d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 6.71262214d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 1.48437776d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.74759022d-04, 6.25244487d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 5.02999748d-04,-2.53355703d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 6.21421776d-04,-1.06382167d-06/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 6.32275978d-04,-1.17335760d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 4.93453805d-04,-4.40747923d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 1.80808273d-04, 2.72349293d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.20546360d-04, 2.18935075d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-1.70664099d-04, 2.65824770d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/ 2.64414201d-04, 4.27765845d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/ 8.64873531d-05, 4.93911850d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.57380937d-03, 4.11786754d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/ 3.53060815d-03, 1.19141134d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/ 6.03527148d-03,-5.67951288d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.10786305d-02,-1.41373598d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 1.13248696d-02,-2.21670479d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/ 9.05905411d-03,-9.36444865d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/ 5.07044732d-03, 1.60348425d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 3.78350788d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 6.32185378d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 5.30289666d-05, 2.55629244d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 1.84763300d-04,-1.18702355d-07/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 3.30716307d-04,-5.00819826d-07/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 3.35774254d-04,-4.94018341d-07/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 2.00671758d-04,-3.31106809d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.27944659d-04, 4.02445576d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-7.48633948d-05, 8.69324505d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 1.22242147d-04, 1.47225640d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/ 5.27070959d-04, 1.91164988d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/ 1.61409467d-03, 2.95314408d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/ 2.66250386d-03, 1.86377380d-06/
      data fpp( 8,13,1),fpp( 8,13,2)/ 4.24295323d-03, 9.67760730d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/ 5.54944825d-03,-2.97481672d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/ 5.48734488d-03,-8.19649386d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/ 5.39904126d-03,-1.10452078d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/ 3.48949913d-03,-3.85467481d-06/
      data fpp( 8,18,1),fpp( 8,18,2)/ 6.18565085d-04, 8.07990709d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 1.86590465d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 2.64382912d-07/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 2.59251118d-05, 9.72341764d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 8.27470507d-05,-6.53196173d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 1.24512995d-04,-2.43955707d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.39027005d-04,-2.54857554d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.15059163d-04,-2.66140771d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.01309076d-06, 1.57313862d-07/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-2.98927809d-05, 3.93358627d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-1.35044909d-05, 8.25251629d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/ 1.23301964d-04, 7.39634858d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/ 4.50733965d-04, 1.41820894d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.10817520d-03, 1.21352939d-06/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.77677892d-03, 3.99673501d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/ 2.56133550d-03,-1.57622340d-06/
      data fpp( 9,15,1),fpp( 9,15,2)/ 3.27678998d-03,-4.19877992d-06/
      data fpp( 9,16,1),fpp( 9,16,2)/ 3.11976539d-03,-5.28065693d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/ 2.09894939d-03,-1.62459237d-06/
      data fpp( 9,18,1),fpp( 9,18,2)/ 7.66492343d-04, 3.90102639d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 8.89248680d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 4.47876957d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/ 2.51018136d-06, 9.42460864d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/ 9.17719769d-06,-2.84861303d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.61028604d-05,-5.74800876d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.08318578d-05,-4.15935193d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-6.31336752d-06,-4.14583526d-09/
      data fpp(10, 7,1),fpp(10, 7,2)/-4.61160175d-06, 7.61768603d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/-1.74899601d-05, 8.94383940d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-6.20760103d-06, 1.90069564d-07/
      data fpp(10,10,1),fpp(10,10,2)/ 4.03586294d-05, 1.16283352d-07/
      data fpp(10,11,1),fpp(10,11,2)/ 1.26150770d-04, 3.28797029d-07/
      data fpp(10,12,1),fpp(10,12,2)/ 2.37822483d-04, 2.36528532d-07/
      data fpp(10,13,1),fpp(10,13,2)/ 4.57386627d-04, 6.30888445d-08/
      data fpp(10,14,1),fpp(10,14,2)/ 6.59269373d-04,-3.68883909d-07/
      data fpp(10,15,1),fpp(10,15,2)/ 7.27957610d-04,-8.91553207d-07/
      data fpp(10,16,1),fpp(10,16,2)/ 6.64183206d-04,-1.03890326d-06/
      data fpp(10,17,1),fpp(10,17,2)/ 4.04002277d-04,-3.10833737d-07/
      data fpp(10,18,1),fpp(10,18,2)/ 9.78404294d-05, 8.00238210d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 1.80188089d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 8.70660004d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 6.34162801d-07,-4.13200077d-10/
      data fpp(11, 3,1),fpp(11, 3,2)/ 1.14415850d-06,-1.30537997d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/ 6.75563052d-07,-1.33716010d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/-5.54436268d-07,-5.45979628d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/-2.20569276d-06, 5.21078613d-09/
      data fpp(11, 7,1),fpp(11, 7,2)/-9.56668375d-06, 2.06166518d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/-8.74737893d-06, 2.63226068d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/-6.06510499d-06, 2.40929209d-08/
      data fpp(11,10,1),fpp(11,10,2)/ 2.66351857d-06, 5.13057097d-08/
      data fpp(11,11,1),fpp(11,11,2)/ 1.30629553d-05, 2.86842403d-08/
      data fpp(11,12,1),fpp(11,12,2)/ 4.53348729d-05, 9.79573290d-08/
      data fpp(11,13,1),fpp(11,13,2)/ 6.74745718d-05,-9.05135564d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 1.10987007d-04,-6.59031033d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 1.55779576d-04,-1.55874030d-07/
      data fpp(11,16,1),fpp(11,16,2)/ 1.50901788d-04,-1.26600775d-07/
      data fpp(11,17,1),fpp(11,17,2)/ 1.01641506d-04,-4.57228682d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 3.85459396d-05, 1.29492248d-07/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 2.77753876d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-5.85052283d-09/
      data fpp(12, 2,1),fpp(12, 2,2)/ 3.53167440d-07,-4.29895434d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.24616830d-06,-9.53659812d-10/
      data fpp(12, 4,1),fpp(12, 4,2)/ 1.59488739d-06,-3.88640641d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/ 3.85887254d-07, 4.49928546d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/-2.26386145d-06, 3.88926456d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/-3.92166325d-06, 3.94365631d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/-5.12052421d-06, 4.33611020d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/-3.13197900d-06, 2.71190290d-09/
      data fpp(12,10,1),fpp(12,10,2)/ 3.58729629d-06, 8.81627821d-09/
      data fpp(12,11,1),fpp(12,11,2)/ 1.53974089d-05, 4.02298426d-09/
      data fpp(12,12,1),fpp(12,12,2)/ 3.32380254d-05, 2.90917847d-08/
      data fpp(12,13,1),fpp(12,13,2)/ 6.05150856d-05,-3.03901232d-08/
      data fpp(12,14,1),fpp(12,14,2)/ 8.71825986d-05,-1.55312919d-08/
      data fpp(12,15,1),fpp(12,15,2)/ 1.00324085d-04,-3.94847091d-08/
      data fpp(12,16,1),fpp(12,16,2)/ 9.24096424d-05,-6.52987180d-09/
      data fpp(12,17,1),fpp(12,17,2)/ 5.82316988d-05, 5.60419625d-09/
      data fpp(12,18,1),fpp(12,18,2)/ 1.67758121d-05, 2.01130868d-08/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 3.39434566d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-2.96941522d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/-9.26583720d-07,-2.06116956d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/-2.21058415d-06,-7.85906547d-10/
      data fpp(13, 4,1),fpp(13, 4,2)/-3.92244369d-06,-7.95204254d-10/
      data fpp(13, 5,1),fpp(13, 5,2)/-3.28044363d-06, 3.96672356d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/-5.05569276d-07, 2.92831001d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 2.59833163d-06, 2.32003641d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/ 4.43526211d-06,-2.08455633d-10/
      data fpp(13, 9,1),fpp(13, 9,2)/ 1.47848950d-06,-1.48621388d-09/
      data fpp(13,10,1),fpp(13,10,2)/-1.19436481d-05, 1.53311136d-10/
      data fpp(13,11,1),fpp(13,11,2)/-3.18737045d-05,-5.12703067d-09/
      data fpp(13,12,1),fpp(13,12,2)/-7.21315127d-05, 2.35481153d-09/
      data fpp(13,13,1),fpp(13,13,2)/-1.25732543d-04, 1.70778453d-09/
      data fpp(13,14,1),fpp(13,14,2)/-1.96441299d-04,-3.18594967d-09/
      data fpp(13,15,1),fpp(13,15,2)/-2.43112042d-04,-9.63985863d-10/
      data fpp(13,16,1),fpp(13,16,2)/-2.29229821d-04, 1.04189312d-09/
      data fpp(13,17,1),fpp(13,17,2)/-1.45115849d-04, 2.79641339d-09/
      data fpp(13,18,1),fpp(13,18,2)/-4.36504060d-05,-2.27546683d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-1.88622666d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 6.75134648d-05, 3.31983040d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/-2.30558140d-05,-2.52358154d-06/
      data fpppp( 1, 5),fpppp( 1, 6)/-8.62145491d-05,-5.11963461d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/ 4.00233076d-05, 4.91079253d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 7.54717045d-05,-6.23535474d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 4.44149028d-06, 1.48046377d-04/
      data fpppp( 1,13),fpppp( 1,14)/-5.15968095d-04, 4.51790584d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 1.81117709d-03,-3.67131674d-03/
      data fpppp( 1,17),fpppp( 1,18)/ 1.00171291d-03, 9.61707681d-04/
      data fpppp( 1,19) /             2.26828385d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 5.82768260d-05, 2.77208807d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/-1.86463685d-05,-9.67405380d-06/
      data fpppp( 2, 5),fpppp( 2, 6)/-6.46137520d-05,-3.72404043d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/ 2.91114784d-05, 4.12805857d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 5.79213233d-05,-3.62111282d-05/
      data fpppp( 2,11),fpppp( 2,12)/ 2.41517932d-06, 9.96271161d-05/
      data fpppp( 2,13),fpppp( 2,14)/-2.53905733d-04, 1.85339511d-04/
      data fpppp( 2,15),fpppp( 2,16)/ 1.11647932d-03,-2.48473974d-03/
      data fpppp( 2,17),fpppp( 2,18)/ 8.27842172d-04, 6.90445929d-04/
      data fpppp( 2,19) /             1.53875773d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 4.94392439d-05, 2.14481478d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-1.36386235d-05,-2.38905315d-05/
      data fpppp( 3, 5),fpppp( 3, 6)/-3.24092183d-05,-1.92166066d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.45078341d-05, 3.13360806d-05/
      data fpppp( 3, 9),fpppp( 3,10)/ 3.18005518d-05, 8.01513346d-07/
      data fpppp( 3,11),fpppp( 3,12)/ 5.42643107d-06, 3.17271516d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 1.68234420d-04,-2.44704184d-04/
      data fpppp( 3,15),fpppp( 3,16)/-2.90715527d-04,-4.34184072d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 8.40578735d-04, 3.22727081d-04/
      data fpppp( 3,19) /             3.17050988d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 3.47695448d-05, 1.48557699d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/-8.24598364d-06,-2.01898122d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/-2.16573182d-05,-6.58789960d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/ 5.67754604d-06, 2.39176175d-05/
      data fpppp( 4, 9),fpppp( 4,10)/ 2.29728600d-05, 1.78217713d-05/
      data fpppp( 4,11),fpppp( 4,12)/-7.89805927d-06,-1.12863503d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 2.34502469d-04,-3.10379190d-04/
      data fpppp( 4,15),fpppp( 4,16)/ 1.21919356d-04,-4.18275040d-04/
      data fpppp( 4,17),fpppp( 4,18)/ 1.98782735d-04, 1.12767690d-04/
      data fpppp( 4,19) /             2.27630594d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 1.84284241d-05, 7.62761991d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-6.54008108d-06,-8.92364323d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/-9.19120385d-06,-5.44278434d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 5.69941324d-06, 5.05095865d-06/
      data fpppp( 5, 9),fpppp( 5,10)/ 2.62383244d-05,-1.63307890d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 2.27709747d-05, 1.50233791d-05/
      data fpppp( 5,13),fpppp( 5,14)/-5.33497493d-05, 4.17048662d-04/
      data fpppp( 5,15),fpppp( 5,16)/-5.60210401d-04,-1.07190561d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 4.05162920d-05, 4.58236321d-05/
      data fpppp( 5,19) /             6.42532601d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 7.66753377d-06, 3.58820064d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 3.09373244d-06,-9.94776301d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/-7.08069853d-06, 6.66513503d-07/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.67427310d-06, 1.87833678d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-6.74636300d-06, 5.37013860d-05/
      data fpppp( 6,11),fpppp( 6,12)/-2.65576395d-05,-8.57368320d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 3.80795352d-05,-2.20634687d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 2.42696152d-04,-9.95411225d-05/
      data fpppp( 6,17),fpppp( 6,18)/-8.57882603d-05, 8.33995978d-05/
      data fpppp( 6,19) /             2.08177457d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 4.96517133d-06, 1.94011659d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-3.51673536d-06,-4.62297087d-07/
      data fpppp( 7, 5),fpppp( 7, 6)/-1.08814586d-06,-4.16570192d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 7.32155198d-06,-9.97748891d-06/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.87314909d-05,-2.13711493d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 2.99727975d-05, 1.39489078d-06/
      data fpppp( 7,13),fpppp( 7,14)/-7.38375416d-06, 6.10119984d-05/
      data fpppp( 7,15),fpppp( 7,16)/-8.43424985d-05,-1.14692003d-05/
      data fpppp( 7,17),fpppp( 7,18)/-2.05039728d-05,-9.88238820d-06/
      data fpppp( 7,19) /            -4.87690560d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.47551569d-06, 6.99863493d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 4.47352351d-07,-1.63615250d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/-2.35644596d-06, 2.65230970d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/-4.51026902d-06, 7.58390908d-06/
      data fpppp( 8, 9),fpppp( 8,10)/-1.83055155d-06, 1.22016933d-05/
      data fpppp( 8,11),fpppp( 8,12)/-6.04452745d-06, 9.65954520d-06/
      data fpppp( 8,13),fpppp( 8,14)/-6.71242358d-07,-2.34118368d-05/
      data fpppp( 8,15),fpppp( 8,16)/ 1.22026857d-05,-2.69709208d-05/
      data fpppp( 8,17),fpppp( 8,18)/-1.35933132d-05, 2.36606593d-05/
      data fpppp( 8,19) /             5.40928133d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.58123805d-07, 3.29713635d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-2.23168716d-07,-3.40398433d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-5.03536499d-08,-1.76709808d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.71405216d-06,-1.00698548d-07/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.52639172d-06, 1.22022153d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 5.03025494d-06,-1.54068751d-06/
      data fpppp( 9,13),fpppp( 9,14)/ 1.80224464d-06, 1.28888044d-06/
      data fpppp( 9,15),fpppp( 9,16)/-1.11038922d-05,-9.22205618d-06/
      data fpppp( 9,17),fpppp( 9,18)/-3.83536742d-06, 5.86506330d-06/
      data fpppp( 9,19) /             1.43329963d-05 /
      data fpppp(10, 1),fpppp(10, 2)/ 8.39290335d-08, 3.48107544d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 2.62380478d-08,-1.24244163d-07/
      data fpppp(10, 5),fpppp(10, 6)/-2.61061314d-07, 4.56036054d-07/
      data fpppp(10, 7),fpppp(10, 8)/-4.32263434d-07, 3.98210239d-07/
      data fpppp(10, 9),fpppp(10,10)/ 2.89065517d-07, 5.62559978d-07/
      data fpppp(10,11),fpppp(10,12)/-1.85750825d-07, 1.73321766d-06/
      data fpppp(10,13),fpppp(10,14)/-2.73573930d-07,-1.69980586d-06/
      data fpppp(10,15),fpppp(10,16)/-9.18873175d-07,-2.57245990d-06/
      data fpppp(10,17),fpppp(10,18)/-5.75678739d-07, 2.11631977d-06/
      data fpppp(10,19) /             4.60968474d-06 /
      data fpppp(11, 1),fpppp(11, 2)/ 7.24478342d-09,-1.12609880d-09/
      data fpppp(11, 3),fpppp(11, 4)/-1.01904142d-08,-1.68277135d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 3.18170359d-08,-1.35715860d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 1.68462335d-07,-4.73157325d-08/
      data fpppp(11, 9),fpppp(11,10)/ 1.32578742d-07,-1.20218259d-07/
      data fpppp(11,11),fpppp(11,12)/ 4.48543081d-07,-3.61605208d-07/
      data fpppp(11,13),fpppp(11,14)/ 3.89944618d-07, 8.41909235d-08/
      data fpppp(11,15),fpppp(11,16)/-6.49900274d-07,-4.64811287d-07/
      data fpppp(11,17),fpppp(11,18)/-1.53804183d-07, 2.49910936d-07/
      data fpppp(11,19) /             6.27138051d-07 /
      data fpppp(12, 1),fpppp(12, 2)/ 1.62356694d-08, 5.40596759d-09/
      data fpppp(12, 3),fpppp(12, 4)/-5.46953457d-09,-1.61847355d-08/
      data fpppp(12, 5),fpppp(12, 6)/-2.32546771d-08, 2.27585299d-08/
      data fpppp(12, 7),fpppp(12, 8)/-8.26262842d-09, 3.78284341d-08/
      data fpppp(12, 9),fpppp(12,10)/ 4.81932625d-08, 5.32423205d-08/
      data fpppp(12,11),fpppp(12,12)/ 4.42876978d-08, 1.31437117d-07/
      data fpppp(12,13),fpppp(12,14)/-3.84953794d-09,-1.52611802d-07/
      data fpppp(12,15),fpppp(12,16)/-1.97264860d-07,-3.21684467d-07/
      data fpppp(12,17),fpppp(12,18)/-9.18073531d-08, 2.52237296d-07/
      data fpppp(12,19) /             5.63662648d-07 /
      data fpppp(13, 1),fpppp(13, 2)/-4.33636989d-09,-6.40913419d-10/
      data fpppp(13, 3),fpppp(13, 4)/-1.45449790d-08, 3.31492826d-08/
      data fpppp(13, 5),fpppp(13, 6)/ 2.31794252d-08, 2.10547342d-09/
      data fpppp(13, 7),fpppp(13, 8)/-1.18597259d-08,-3.06847949d-08/
      data fpppp(13, 9),fpppp(13,10)/-1.53023280d-07, 1.48560117d-08/
      data fpppp(13,11),fpppp(13,12)/-2.96875888d-07,-4.70175729d-08/
      data fpppp(13,13),fpppp(13,14)/-3.15647133d-07, 2.83142525d-07/
      data fpppp(13,15),fpppp(13,16)/ 6.25357837d-07, 8.48603980d-07/
      data fpppp(13,17),fpppp(13,18)/ 1.94131283d-07,-5.84040821d-07/
      data fpppp(13,19) /            -1.32687024d-06 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine c3_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -2.06000000d-05 /
      data f( 1, 3),f( 1, 4) / -3.47000000d-05 ,  1.58100000d-04 /
      data f( 1, 5),f( 1, 6) /  5.41100000d-04 ,  7.56500000d-04 /
      data f( 1, 7),f( 1, 8) /  5.43100000d-04 ,  3.86000000d-05 /
      data f( 1, 9),f( 1,10) / -2.83000000d-04 , -2.67000000d-04 /
      data f( 1,11),f( 1,12) /  2.55200000d-04 ,  1.47100000d-03 /
      data f( 1,13),f( 1,14) /  4.63500000d-03 ,  1.33574000d-02 /
      data f( 1,15),f( 1,16) /  2.59861000d-02 ,  2.81951000d-02 /
      data f( 1,17),f( 1,18) /  1.35233000d-02 ,  3.26570000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.02000000d-05 /
      data f( 2, 3),f( 2, 4) /  9.00000000d-07 ,  1.60400000d-04 /
      data f( 2, 5),f( 2, 6) /  4.35200000d-04 ,  5.71000000d-04 /
      data f( 2, 7),f( 2, 8) /  3.82100000d-04 ,  5.80000000d-06 /
      data f( 2, 9),f( 2,10) / -2.22300000d-04 , -2.22900000d-04 /
      data f( 2,11),f( 2,12) /  1.94800000d-04 ,  1.08060000d-03 /
      data f( 2,13),f( 2,14) /  3.24170000d-03 ,  9.36110000d-03 /
      data f( 2,15),f( 2,16) /  1.78202000d-02 ,  1.94053000d-02 /
      data f( 2,17),f( 2,18) /  1.16608000d-02 ,  3.91660000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.70000000d-06 /
      data f( 3, 3),f( 3, 4) /  2.19000000d-05 ,  1.48800000d-04 /
      data f( 3, 5),f( 3, 6) /  3.41700000d-04 ,  4.18100000d-04 /
      data f( 3, 7),f( 3, 8) /  2.60300000d-04 , -1.36000000d-05 /
      data f( 3, 9),f( 3,10) / -1.78100000d-04 , -1.71600000d-04 /
      data f( 3,11),f( 3,12) /  1.49700000d-04 ,  7.94300000d-04 /
      data f( 3,13),f( 3,14) /  2.21970000d-03 ,  6.19540000d-03 /
      data f( 3,15),f( 3,16) /  1.20891000d-02 ,  1.38515000d-02 /
      data f( 3,17),f( 3,18) /  9.48960000d-03 ,  3.20150000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  2.80000000d-06 /
      data f( 4, 3),f( 4, 4) /  3.16000000d-05 ,  1.19600000d-04 /
      data f( 4, 5),f( 4, 6) /  2.28900000d-04 ,  2.54500000d-04 /
      data f( 4, 7),f( 4, 8) /  1.43500000d-04 , -2.56000000d-05 /
      data f( 4, 9),f( 4,10) / -1.32400000d-04 , -1.17500000d-04 /
      data f( 4,11),f( 4,12) /  1.03300000d-04 ,  5.13500000d-04 /
      data f( 4,13),f( 4,14) /  1.26540000d-03 ,  3.14300000d-03 /
      data f( 4,15),f( 4,16) /  6.32990000d-03 ,  7.81580000d-03 /
      data f( 4,17),f( 4,18) /  5.80570000d-03 ,  1.78240000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  3.10000000d-06 /
      data f( 5, 3),f( 5, 4) /  2.43000000d-05 ,  6.88000000d-05 /
      data f( 5, 5),f( 5, 6) /  1.10500000d-04 ,  1.09200000d-04 /
      data f( 5, 7),f( 5, 8) /  5.38000000d-05 , -2.30000000d-05 /
      data f( 5, 9),f( 5,10) / -8.60000000d-05 , -6.88000000d-05 /
      data f( 5,11),f( 5,12) /  5.60000000d-05 ,  2.39600000d-04 /
      data f( 5,13),f( 5,14) /  5.76900000d-04 ,  1.11800000d-03 /
      data f( 5,15),f( 5,16) /  2.80690000d-03 ,  3.18870000d-03 /
      data f( 5,17),f( 5,18) /  2.32040000d-03 ,  6.28200000d-04 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.68000000d-05 /
      data f( 6, 3),f( 6, 4) /  3.06000000d-05 ,  3.62000000d-05 /
      data f( 6, 5),f( 6, 6) /  4.82000000d-05 ,  4.66000000d-05 /
      data f( 6, 7),f( 6, 8) /  2.22000000d-05 ,  4.00000000d-07 /
      data f( 6, 9),f( 6,10) / -7.66000000d-05 ,  1.43400000d-04 /
      data f( 6,11),f( 6,12) /  2.92700000d-04 ,  4.19000000d-04 /
      data f( 6,13),f( 6,14) /  5.77400000d-04 ,  7.77200000d-04 /
      data f( 6,15),f( 6,16) /  1.14200000d-03 ,  1.49190000d-03 /
      data f( 6,17),f( 6,18) /  1.18870000d-03 ,  3.69900000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.00000000d-06 /
      data f( 7, 3),f( 7, 4) /  6.70000000d-06 ,  1.73000000d-05 /
      data f( 7, 5),f( 7, 6) /  2.51000000d-05 ,  2.13000000d-05 /
      data f( 7, 7),f( 7, 8) /  1.13000000d-05 , -1.28000000d-05 /
      data f( 7, 9),f( 7,10) / -1.64000000d-05 , -2.53000000d-05 /
      data f( 7,11),f( 7,12) /  5.10000000d-06 ,  6.43000000d-05 /
      data f( 7,13),f( 7,14) /  1.49400000d-04 ,  2.69700000d-04 /
      data f( 7,15),f( 7,16) /  3.90700000d-04 ,  5.32800000d-04 /
      data f( 7,17),f( 7,18) /  3.59000000d-04 ,  7.28000000d-05 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  4.00000000d-07 /
      data f( 8, 3),f( 8, 4) /  3.30000000d-06 ,  8.40000000d-06 /
      data f( 8, 5),f( 8, 6) /  1.17000000d-05 ,  1.00000000d-05 /
      data f( 8, 7),f( 8, 8) / -4.90000000d-06 , -9.00000000d-06 /
      data f( 8, 9),f( 8,10) / -1.22000000d-05 , -2.01000000d-05 /
      data f( 8,11),f( 8,12) / -7.20000000d-06 ,  1.90000000d-05 /
      data f( 8,13),f( 8,14) /  6.26000000d-05 ,  1.22200000d-04 /
      data f( 8,15),f( 8,16) /  1.93900000d-04 ,  2.20600000d-04 /
      data f( 8,17),f( 8,18) /  1.35500000d-04 ,  2.44000000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f( 9, 3),f( 9, 4) /  1.50000000d-06 ,  4.10000000d-06 /
      data f( 9, 5),f( 9, 6) /  6.40000000d-06 ,  5.50000000d-06 /
      data f( 9, 7),f( 9, 8) / -2.00000000d-06 , -5.50000000d-06 /
      data f( 9, 9),f( 9,10) / -9.60000000d-06 , -1.92000000d-05 /
      data f( 9,11),f( 9,12) / -1.36000000d-05 , -3.30000000d-06 /
      data f( 9,13),f( 9,14) /  1.61000000d-05 ,  4.46000000d-05 /
      data f( 9,15),f( 9,16) /  7.49000000d-05 ,  8.49000000d-05 /
      data f( 9,17),f( 9,18) /  5.17000000d-05 ,  9.90000000d-06 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(10, 3),f(10, 4) /  3.00000000d-07 ,  9.00000000d-07 /
      data f(10, 5),f(10, 6) /  1.30000000d-06 ,  2.00000000d-07 /
      data f(10, 7),f(10, 8) / -1.30000000d-06 , -2.70000000d-06 /
      data f(10, 9),f(10,10) / -4.70000000d-06 , -9.80000000d-06 /
      data f(10,11),f(10,12) / -1.32000000d-05 , -1.44000000d-05 /
      data f(10,13),f(10,14) / -1.17000000d-05 , -7.30000000d-06 /
      data f(10,15),f(10,16) / -2.10000000d-06 ,  1.80000000d-06 /
      data f(10,17),f(10,18) /  2.50000000d-06 ,  7.00000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(11, 3),f(11, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(11, 5),f(11, 6) /  0.00000000d+00 , -5.00000000d-07 /
      data f(11, 7),f(11, 8) / -1.00000000d-06 , -1.40000000d-06 /
      data f(11, 9),f(11,10) / -1.80000000d-06 , -2.80000000d-06 /
      data f(11,11),f(11,12) / -5.70000000d-06 , -8.30000000d-06 /
      data f(11,13),f(11,14) / -1.02000000d-05 , -9.60000000d-06 /
      data f(11,15),f(11,16) / -8.90000000d-06 , -6.50000000d-06 /
      data f(11,17),f(11,18) / -2.80000000d-06 , -5.00000000d-07 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) / -1.00000000d-07 , -1.00000000d-07 /
      data f(12, 5),f(12, 6) / -2.00000000d-07 , -4.00000000d-07 /
      data f(12, 7),f(12, 8) / -6.00000000d-07 , -7.00000000d-07 /
      data f(12, 9),f(12,10) / -8.00000000d-07 , -1.00000000d-06 /
      data f(12,11),f(12,12) / -1.80000000d-06 , -2.90000000d-06 /
      data f(12,13),f(12,14) / -4.40000000d-06 , -4.70000000d-06 /
      data f(12,15),f(12,16) / -4.80000000d-06 , -3.60000000d-06 /
      data f(12,17),f(12,18) / -1.50000000d-06 , -2.00000000d-07 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 5),f(13, 6) / -1.00000000d-07 , -2.00000000d-07 /
      data f(13, 7),f(13, 8) / -2.00000000d-07 , -3.00000000d-07 /
      data f(13, 9),f(13,10) / -3.00000000d-07 , -4.00000000d-07 /
      data f(13,11),f(13,12) / -8.00000000d-07 , -1.20000000d-06 /
      data f(13,13),f(13,14) / -1.70000000d-06 , -1.80000000d-06 /
      data f(13,15),f(13,16) / -1.40000000d-06 , -8.00000000d-07 /
      data f(13,17),f(13,18) / -3.00000000d-07 , -1.00000000d-07 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-1.88749388d-06/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.49726817d-06,-3.80122459d-08/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-4.29907670d-04, 2.42954286d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-5.06951959d-04, 2.73384080d-06/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 2.63690595d-04,-1.95290608d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 7.73986014d-04,-4.97821650d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.07658400d-03,-3.86222791d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 4.28673959d-04, 2.96112816d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-5.33242319d-04, 2.99171527d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/ 5.97501472d-04, 5.32801077d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 5.14131131d-04, 6.06824165d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 3.20551979d-03, 1.20150226d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.07215983d-02, 6.27636679d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.95720111d-02, 7.04343059d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.93429628d-02,-1.10122891d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/ 1.21217099d-01,-2.55124740d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/-1.04850818d-02, 1.17773853d-04/
      data fpp( 1,18,1),fpp( 1,18,2)/-5.85170797d-02, 4.88813277d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 1.06214836d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-1.01622318d-06/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-3.83516792d-05, 1.29446360d-07/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-3.62613231d-04, 1.77643774d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-3.54596083d-04, 1.66880268d-06/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.06618810d-04,-1.53364846d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 7.91313687d-04,-3.87420884d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 9.60546277d-04,-2.45151619d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 3.31080653d-04, 2.43627362d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-4.07301077d-04, 1.59842173d-06/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.75854198d-04, 4.82003945d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 2.89237738d-04, 4.21942047d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 2.46017471d-03, 6.38827868d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/ 9.08973198d-03, 4.67454648d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.99789778d-02, 4.41278621d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/ 6.34166458d-02,-8.28749133d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 8.64943738d-02,-1.25068209d-04/
      data fpp( 2,17,1),fpp( 2,17,2)/-5.65005062d-03, 2.33717485d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-3.60398406d-02, 3.15992147d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 7.98873926d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-4.95723441d-07/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-1.33090551d-04, 1.90446882d-07/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.09639406d-04, 1.25193591d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.59663711d-04, 9.99809463d-07/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 3.69834164d-04,-1.29117377d-06/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 9.50759238d-04,-2.82511440d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 9.61230888d-04,-1.46036863d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 2.57003428d-04, 1.70058891d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-3.12553374d-04, 1.22201300d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/-2.20918266d-04, 3.67135911d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 6.23917917d-04, 2.98055056d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.56878137d-03, 3.80443865d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/ 8.61447380d-03, 2.86496948d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.51020779d-02, 3.46147820d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.22104539d-02,-5.20288228d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.82054061d-02,-7.43774908d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.32197157d-02,-1.79192140d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/-2.22355803d-03, 3.04823469d-05/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 8.11858266d-05/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-5.06158420d-08/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-8.07970426d-05, 2.17231684d-07/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.79459827d-04, 7.41689106d-07/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.80569091d-05, 3.68011891d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 3.92806914d-04,-9.35736670d-07/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 6.86593417d-04,-1.64706521d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 5.48866189d-04,-6.72002489d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 6.26014722d-05, 8.49075167d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-5.99547008d-05, 1.01370182d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/-9.04175247d-04, 2.39811754d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/-8.55884882d-04, 1.74782800d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/-2.92721033d-04, 1.97457045d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.80526601d-03, 1.08558902d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.60837553d-02, 2.21438688d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 3.95173896d-02,-2.08733654d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/ 3.46523971d-02,-4.07104072d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.93590860d-02,-2.60450058d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 8.34175382d-03, 2.40984302d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 6.41052849d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 1.39903673d-07/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.65604867d-04, 1.59192654d-07/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.96855091d-04, 3.09325712d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 1.02380336d-04, 1.50449694d-09/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 1.91517378d-04,-4.83343700d-07/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.89245524d-04,-6.48129696d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.86089664d-04,-1.70137514d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.87873232d-04, 4.46797532d-08/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-3.35012933d-04, 8.19418502d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/ 2.03071175d-03, 1.48964624d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 3.08528087d-03,-3.22003465d-07/
      data fpp( 5,12,1),fpp( 5,12,2)/ 4.05383848d-03, 3.32636762d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 4.30246448d-03,-3.76146701d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 6.96673627d-03, 2.39475004d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 3.40807652d-05,-2.31605346d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 8.56668567d-03,-9.73136185d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 9.69155408d-03,-1.29200180d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 3.70372258d-03, 1.19774337d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 2.88502831d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 1.04437013d-08/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-2.60022426d-04,-6.88740254d-09/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-2.81560535d-04,-1.62894091d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 4.53355666d-05, 1.66463767d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.87523573d-04,-1.18960977d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.41224487d-04,-5.06619859d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 1.01175157d-04, 7.77440414d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-3.14894401d-04,-2.44714179d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 5.12006432d-04, 5.69912677d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-3.29467175d-03,-2.52936527d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-4.66923861d-03, 1.76334310d-07/
      data fpp( 6,12,1),fpp( 6,12,2)/-5.04343291d-03, 4.44028028d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/-4.47912391d-03,-2.64464233d-08/
      data fpp( 6,14,1),fpp( 6,14,2)/-3.52990038d-03, 2.14575766d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 4.94068729d-03, 1.34341576d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.40806026d-03,-8.41342072d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/-1.63890235d-03,-6.87573288d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.65504416d-03, 4.98035225d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 1.38883239d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 4.74854216d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.66484838d-04, 4.20291567d-08/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 2.04587049d-04, 6.63979515d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 4.50773981d-05,-1.36209627d-08/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-8.11671235d-07,-1.79914101d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 4.10565272d-05, 3.72773660d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-9.39902922d-05,-3.41195363d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.93304373d-04, 4.81504087d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-4.93812794d-04,-3.54820985d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/ 2.00637525d-03, 6.19779854d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 3.00847355d-03, 2.33701571d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 3.30149314d-03, 1.73413864d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/ 3.33003118d-03, 6.26642973d-07/
      data fpp( 7,14,1),fpp( 7,14,2)/ 3.15206527d-03,-5.67985756d-07/
      data fpp( 7,15,1),fpp( 7,15,2)/ 2.12957006d-03, 1.68730005d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 3.50587331d-03,-4.91521445d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 4.11205531d-03,-9.80442264d-07/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.98525407d-03, 2.09298350d-06/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 5.41250825d-06/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 3.00967788d-08/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-4.11169258d-05, 2.08064425d-08/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-4.47876603d-05, 3.66774513d-08/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 1.43548410d-05,-3.55162476d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 4.85231117d-05,-2.61246083d-09/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 3.05494041d-05,-2.54033909d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.47586012d-04, 2.26748097d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-5.03230903d-05,-4.95847913d-09/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 1.19244743d-04,-1.52914181d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/-5.57229234d-04, 3.34615201d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/-7.57455596d-04, 6.24533757d-08/
      data fpp( 8,12,1),fpp( 8,12,2)/-7.36939662d-04, 2.13571296d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-6.52200809d-04, 1.27261440d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/-4.38360684d-04, 2.37382945d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/-1.50967545d-04,-3.50793220d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/ 9.40464919d-05,-1.53421006d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/-2.60518878d-04,-2.20366521d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/-3.17172098d-04, 8.55676149d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 1.99966193d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 9.75518945d-09/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 7.58286529d-06, 9.48962111d-09/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.29635925d-05, 1.82863261d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 7.90323786d-06,-4.63492561d-09/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.11922440d-06,-1.77466237d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-5.41435316d-08,-1.16378580d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-3.79537547d-05, 8.72609424d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 7.87988167d-07, 7.33480992d-09/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-2.15661771d-05,-1.52600182d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/ 1.19341690d-04, 2.73065919d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.62948833d-04,-2.76634920d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.98265505d-04, 1.19588049d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/ 2.45972055d-04, 9.53112942d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/ 2.78977471d-04, 4.51667739d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/ 3.41500118d-04,-1.67978390d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 3.53940719d-04,-5.91253215d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/ 2.82820207d-04,-5.90087498d-08/
      data fpp( 9,18,1),fpp( 9,18,2)/ 9.70343272d-05, 3.11288214d-07/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 7.27855893d-07/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 7.86548841d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.59013297d-06, 5.26902318d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/-2.09694728d-06, 1.05841887d-09/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.51286593d-06, 2.49730132d-09/
      data fpp(10, 5,1),fpp(10, 5,2)/ 5.38077093d-06,-2.30476242d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/ 7.08772855d-06,-3.06804639d-10/
      data fpp(10, 7,1),fpp(10, 7,2)/ 9.46825831d-06, 2.74842728d-10/
      data fpp(10, 8,1),fpp(10, 8,2)/-2.40241937d-06, 5.20743373d-09/
      data fpp(10, 9,1),fpp(10, 9,2)/ 3.27616008d-06,-5.71045776d-08/
      data fpp(10,10,1),fpp(10,10,2)/-3.38104522d-05, 3.72108768d-08/
      data fpp(10,11,1),fpp(10,11,2)/-3.09187019d-05, 1.02610704d-08/
      data fpp(10,12,1),fpp(10,12,2)/-2.53266845d-05, 5.37448415d-08/
      data fpp(10,13,1),fpp(10,13,2)/-2.06157613d-05, 8.75956349d-09/
      data fpp(10,14,1),fpp(10,14,2)/ 2.04793040d-06, 1.32169045d-08/
      data fpp(10,15,1),fpp(10,15,2)/ 1.69834192d-05,-1.36271816d-08/
      data fpp(10,16,1),fpp(10,16,2)/ 2.09545963d-05,-3.67081781d-08/
      data fpp(10,17,1),fpp(10,17,2)/-7.80118218d-06,-3.15401059d-08/
      data fpp(10,18,1),fpp(10,18,2)/-1.37169326d-05, 1.28686017d-08/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 4.60656992d-08/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 3.98017153d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 5.77666570d-07, 2.03965694d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/ 8.24196654d-07,-1.38799297d-10/
      data fpp(11, 4,1),fpp(11, 4,2)/ 4.45298441d-07,-1.48445975d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.57691859d-07,-5.92336169d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/-6.96770671d-07, 1.17790651d-09/
      data fpp(11, 7,1),fpp(11, 7,2)/-2.31927850d-06, 1.21173566d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/-1.78310692d-07,-2.48491341d-11/
      data fpp(11, 9,1),fpp(11, 9,2)/-3.53846318d-06,-1.11233912d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 1.50011899d-06,-3.15257944d-08/
      data fpp(11,11,1),fpp(11,11,2)/ 3.32597417d-06, 1.32155167d-08/
      data fpp(11,12,1),fpp(11,12,2)/ 6.24123276d-06,-3.33627225d-09/
      data fpp(11,13,1),fpp(11,13,2)/ 1.22909898d-05, 4.21295724d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 1.04308078d-05,-1.51820172d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 1.17662055d-05, 2.45984963d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 1.10408957d-05, 1.87880319d-08/
      data fpp(11,17,1),fpp(11,17,2)/ 1.17845216d-05,-2.17506240d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 5.83340333d-06,-1.57855360d-08/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-2.31072320d-08/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-2.81442099d-09/
      data fpp(12, 2,1),fpp(12, 2,2)/-1.20533314d-07,-1.37115802d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.60669287d-10, 2.29905306d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/ 3.05940312d-07,-1.82505421d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/ 5.88461628d-07,-9.98836232d-10/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.99354134d-07,-1.79600868d-10/
      data fpp(12, 7,1),fpp(12, 7,2)/ 4.08855700d-07, 1.71723970d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/-4.84337862d-07,-6.89357945d-10/
      data fpp(12, 9,1),fpp(12, 9,2)/-5.22307364d-07, 1.04019208d-09/
      data fpp(12,10,1),fpp(12,10,2)/-3.39002380d-06,-9.47141037d-09/
      data fpp(12,11,1),fpp(12,11,2)/-3.98519483d-06, 8.45449390d-10/
      data fpp(12,12,1),fpp(12,12,2)/-3.83824655d-06,-1.19103872d-08/
      data fpp(12,13,1),fpp(12,13,2)/-2.74819796d-06, 2.27960994d-08/
      data fpp(12,14,1),fpp(12,14,2)/-5.71161558d-07,-7.27401035d-09/
      data fpp(12,15,1),fpp(12,15,2)/ 1.35175890d-06, 1.82999420d-08/
      data fpp(12,16,1),fpp(12,16,2)/ 2.08182086d-06, 1.20742424d-08/
      data fpp(12,17,1),fpp(12,17,2)/ 2.63095675d-07,-1.25969114d-08/
      data fpp(12,18,1),fpp(12,18,2)/-6.16680665d-07,-9.68659674d-09/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-1.46567016d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 6.02167865d-11/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.27233343d-07,-1.20433573d-10/
      data fpp(13, 3,1),fpp(13, 3,2)/ 3.74196654d-08, 4.21517505d-10/
      data fpp(13, 4,1),fpp(13, 4,2)/-3.90470156d-07,-1.56563645d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/-1.09423081d-06,-1.58971713d-10/
      data fpp(13, 6,1),fpp(13, 6,2)/-1.14967707d-06, 2.20152330d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/-6.66927850d-07,-2.64712148d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/ 4.21689308d-08, 2.38696263d-09/
      data fpp(13, 9,1),fpp(13, 9,2)/ 1.08615368d-06,-9.00729023d-10/
      data fpp(13,10,1),fpp(13,10,2)/ 4.92001190d-06,-4.78404653d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 9.25974174d-08, 2.03691515d-09/
      data fpp(13,12,1),fpp(13,12,2)/-5.25587672d-06,-3.36361408d-09/
      data fpp(13,13,1),fpp(13,13,2)/-1.12509010d-05, 5.41754118d-09/
      data fpp(13,14,1),fpp(13,14,2)/-1.38519192d-05, 5.69344937d-09/
      data fpp(13,15,1),fpp(13,15,2)/-1.71383795d-05, 1.80866136d-09/
      data fpp(13,16,1),fpp(13,16,2)/-1.62659104d-05,-9.28094797d-10/
      data fpp(13,17,1),fpp(13,17,2)/-8.78154784d-06,-4.09628217d-09/
      data fpp(13,18,1),fpp(13,18,2)/-1.81665967d-06,-6.86776523d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 8.43388261d-10/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-1.22993652d-05,-4.13359285d-06/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.85960423d-06, 1.39568149d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/-7.82565341d-06, 1.72497063d-06/
      data fpppp( 1, 7),fpppp( 1, 8)/-1.15360748d-05,-1.26111537d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.31403154d-05,-3.43905040d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 2.15748527d-05, 1.14576633d-04/
      data fpppp( 1,13),fpppp( 1,14)/-1.90399994d-04, 7.27083403d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 3.37298714d-04,-3.15008720d-03/
      data fpppp( 1,17),fpppp( 1,18)/ 1.84847108d-03, 7.76413847d-04/
      data fpppp( 1,19) /             1.43881819d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-9.07031505d-06,-2.80043753d-06/
      data fpppp( 2, 3),fpppp( 2, 4)/ 3.11747282d-06, 1.02672683d-05/
      data fpppp( 2, 5),fpppp( 2, 6)/-4.99468115d-06,-8.79744605d-07/
      data fpppp( 2, 7),fpppp( 2, 8)/-1.04140777d-05,-5.38583754d-06/
      data fpppp( 2, 9),fpppp( 2,10)/ 2.54224614d-05,-1.70117878d-05/
      data fpppp( 2,11),fpppp( 2,12)/ 1.44383858d-05, 8.27114505d-05/
      data fpppp( 2,13),fpppp( 2,14)/-7.77669703d-05, 4.83937742d-04/
      data fpppp( 2,15),fpppp( 2,16)/ 9.49213395d-05,-2.08521951d-03/
      data fpppp( 2,17),fpppp( 2,18)/ 1.33262755d-03, 4.59987390d-04/
      data fpppp( 2,19) /             8.13200723d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-4.03501824d-06,-6.33369648d-07/
      data fpppp( 3, 3),fpppp( 3, 4)/ 3.96099866d-06, 4.38084796d-06/
      data fpppp( 3, 5),fpppp( 3, 6)/ 1.28694027d-06,-6.44297709d-06/
      data fpppp( 3, 7),fpppp( 3, 8)/-9.74223731d-06, 2.52997968d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 7.70255810d-06, 6.33130258d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 1.21642960d-05, 1.10131495d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 1.89832845d-04,-1.43829832d-04/
      data fpppp( 3,15),fpppp( 3,16)/-1.77267194d-04,-4.13906826d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 7.87690059d-04,-1.91576643d-04/
      data fpppp( 3,19) /            -5.47739464d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.87919649d-06,-3.36208181d-07/
      data fpppp( 4, 3),fpppp( 4, 4)/ 3.15208469d-06, 3.33181157d-06/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.51167671d-06,-4.30974391d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/-7.14017154d-06, 1.19581808d-05/
      data fpppp( 4, 9),fpppp( 4,10)/-1.88700390d-05, 2.02221127d-05/
      data fpppp( 4,11),fpppp( 4,12)/-8.46775717d-06, 4.45413250d-05/
      data fpppp( 4,13),fpppp( 4,14)/ 4.23918490d-05, 2.76721415d-04/
      data fpppp( 4,15),fpppp( 4,16)/-4.79968806d-04,-5.47638030d-05/
      data fpppp( 4,17),fpppp( 4,18)/ 7.33249100d-05, 1.80228920d-05/
      data fpppp( 4,19) /             1.51182257d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.59739589d-06,-1.00844085d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.43011936d-06, 3.18541960d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/ 7.05148794d-07,-5.49054860d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 9.20400523d-06,-2.50291066d-05/
      data fpppp( 5, 9),fpppp( 5,10)/ 5.94322370d-05,-3.93831908d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 1.94311925d-05,-4.35022699d-05/
      data fpppp( 5,13),fpppp( 5,14)/ 1.11381990d-04,-2.57086942d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 3.41150139d-04,-1.79597989d-04/
      data fpppp( 5,17),fpppp( 5,18)/-6.72223715d-05, 2.17254812d-05/
      data fpppp( 5,19) /             1.17366981d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.58927238d-06, 1.77698694d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 5.61183893d-06,-3.31829005d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.42116439d-06, 5.69372207d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.89787385d-05, 4.76600183d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-9.70831111d-05, 6.26576854d-05/
      data fpppp( 6,11),fpppp( 6,12)/-7.62095119d-06, 2.78484727d-05/
      data fpppp( 6,13),fpppp( 6,14)/-4.74627418d-05, 1.85097367d-04/
      data fpppp( 6,15),fpppp( 6,16)/-2.41644877d-04, 6.12892580d-05/
      data fpppp( 6,17),fpppp( 6,18)/ 2.56277111d-05, 1.80491446d-05/
      data fpppp( 6,19) /             2.44686937d-06 /
      data fpppp( 7, 1),fpppp( 7, 2)/-7.63652992d-07,-9.39588152d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/-3.18095204d-06, 1.80668462d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 2.77144845d-06,-7.62704236d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.71218199d-05,-3.55197482d-05/
      data fpppp( 7, 9),fpppp( 7,10)/ 6.64924631d-05,-3.92117918d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 4.69319863d-07,-5.21021048d-06/
      data fpppp( 7,13),fpppp( 7,14)/ 4.50262886d-06,-2.51905421d-05/
      data fpppp( 7,15),fpppp( 7,16)/ 4.55877822d-05,-1.32326794d-05/
      data fpppp( 7,17),fpppp( 7,18)/-3.88643400d-05, 4.71104537d-06/
      data fpppp( 7,19) /             2.85129890d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.59010874d-07, 2.98023103d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 8.95668196d-07,-1.11901738d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/-1.94651508d-06, 4.76944336d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/-9.03063945d-06, 1.24563718d-05/
      data fpppp( 8, 9),fpppp( 8,10)/-1.87462318d-05, 1.17660470d-05/
      data fpppp( 8,11),fpppp( 8,12)/ 2.56900807d-07, 4.50887564d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 1.79292406d-06, 1.23492507d-07/
      data fpppp( 8,15),fpppp( 8,16)/ 2.12628678d-06,-1.11713858d-05/
      data fpppp( 8,17),fpppp( 8,18)/ 6.58449186d-06, 2.70814736d-06/
      data fpppp( 8,19) /             5.01243776d-06 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 5.67387806d-08,-1.47628304d-08/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.29815745d-07,-9.24390964d-08/
      data fpppp( 9, 5),fpppp( 9, 6)/ 3.96152601d-07,-1.15553258d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/ 2.02240312d-06,-2.33559864d-06/
      data fpppp( 9, 9),fpppp( 9,10)/ 3.65423696d-06,-2.48562728d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 4.50228761d-07, 1.87283932d-07/
      data fpppp( 9,13),fpppp( 9,14)/-4.55971793d-07, 7.54535152d-07/
      data fpppp( 9,15),fpppp( 9,16)/-7.91134902d-07,-5.94918289d-07/
      data fpppp( 9,17),fpppp( 9,18)/-1.84285876d-06, 1.08643126d-06/
      data fpppp( 9,19) /             2.82222686d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-1.69609135d-08, 5.75520895d-09/
      data fpppp(10, 3),fpppp(10, 4)/ 5.89391966d-08, 5.48565610d-09/
      data fpppp(10, 5),fpppp(10, 6)/-6.53963131d-08, 1.26442753d-07/
      data fpppp(10, 7),fpppp(10, 8)/-3.99960369d-07, 6.18326279d-07/
      data fpppp(10, 9),fpppp(10,10)/-1.02038932d-06, 8.97319493d-07/
      data fpppp(10,11),fpppp(10,12)/-1.70186901d-07,-5.45558653d-08/
      data fpppp(10,13),fpppp(10,14)/ 3.35544711d-07,-2.10456871d-07/
      data fpppp(10,15),fpppp(10,16)/ 4.25905976d-08,-6.17764222d-07/
      data fpppp(10,17),fpppp(10,18)/ 4.64848959d-07, 1.28770066d-07/
      data fpppp(10,19) /             1.98031765d-07 /
      data fpppp(11, 1),fpppp(11, 2)/-8.68821572d-10,-2.31061525d-09/
      data fpppp(11, 3),fpppp(11, 4)/-9.75690662d-09, 3.81254395d-09/
      data fpppp(11, 5),fpppp(11, 6)/-1.57712674d-11,-3.77608159d-08/
      data fpppp(11, 7),fpppp(11, 8)/ 1.04976317d-07,-1.56335912d-07/
      data fpppp(11, 9),fpppp(11,10)/ 1.90300115d-07,-1.00940470d-07/
      data fpppp(11,11),fpppp(11,12)/ 2.06981449d-08, 8.35120942d-08/
      data fpppp(11,13),fpppp(11,14)/-1.66676614d-07, 1.08598016d-07/
      data fpppp(11,15),fpppp(11,16)/-7.59806693d-08, 7.16822140d-08/
      data fpppp(11,17),fpppp(11,18)/-1.22612047d-07, 1.70813227d-08/
      data fpppp(11,19) /             6.13496551d-08 /
      data fpppp(12, 1),fpppp(12, 2)/ 2.98054598d-09, 2.39855971d-09/
      data fpppp(12, 3),fpppp(12, 4)/ 1.89885300d-09, 1.11116786d-09/
      data fpppp(12, 5),fpppp(12, 6)/-7.73902401d-09, 7.54719956d-09/
      data fpppp(12, 7),fpppp(12, 8)/-2.25332307d-08, 3.44240154d-08/
      data fpppp(12, 9),fpppp(12,10)/-6.38493873d-08, 5.11887178d-08/
      data fpppp(12,11),fpppp(12,12)/-4.55275999d-09, 1.15494814d-08/
      data fpppp(12,13),fpppp(12,14)/ 1.49408528d-08,-6.09362393d-09/
      data fpppp(12,15),fpppp(12,16)/-5.81331364d-09,-4.22246321d-08/
      data fpppp(12,17),fpppp(12,18)/ 2.17846143d-08, 1.14231055d-08/
      data fpppp(12,19) /             2.23103840d-08 /
      data fpppp(13, 1),fpppp(13, 2)/ 1.61885071d-08, 6.06816002d-09/
      data fpppp(13, 3),fpppp(13, 4)/-1.09479661d-08,-3.82886555d-09/
      data fpppp(13, 5),fpppp(13, 6)/ 9.71117806d-09, 3.88301761d-09/
      data fpppp(13, 7),fpppp(13, 8)/ 7.04847970d-09,-1.84960826d-08/
      data fpppp(13, 9),fpppp(13,10)/ 8.70291288d-08,-1.62228025d-07/
      data fpppp(13,11),fpppp(13,12)/ 4.22066080d-08,-3.78619868d-08/
      data fpppp(13,13),fpppp(13,14)/ 7.04483299d-08,-4.02909673d-08/
      data fpppp(13,15),fpppp(13,16)/ 4.95890176d-08, 9.14706524d-08/
      data fpppp(13,17),fpppp(13,18)/-1.87580131d-08,-4.76070654d-08/
      data fpppp(13,19) /            -9.97074357d-08 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine c4_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -1.04000000d-05 /
      data f( 1, 3),f( 1, 4) / -1.15500000d-04 , -3.02800000d-04 /
      data f( 1, 5),f( 1, 6) / -3.88700000d-04 , -2.83100000d-04 /
      data f( 1, 7),f( 1, 8) / -1.76900000d-04 , -1.98300000d-04 /
      data f( 1, 9),f( 1,10) / -1.36800000d-04 ,  6.65000000d-05 /
      data f( 1,11),f( 1,12) /  4.51200000d-04 ,  9.92200000d-04 /
      data f( 1,13),f( 1,14) /  1.82520000d-03 ,  3.02440000d-03 /
      data f( 1,15),f( 1,16) /  3.19030000d-03 ,  3.86320000d-03 /
      data f( 1,17),f( 1,18) /  8.14500000d-04 ,  2.05500000d-04 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.08000000d-05 /
      data f( 2, 3),f( 2, 4) / -9.17000000d-05 , -2.19000000d-04 /
      data f( 2, 5),f( 2, 6) / -2.54600000d-04 , -1.78500000d-04 /
      data f( 2, 7),f( 2, 8) / -1.37900000d-04 , -1.66300000d-04 /
      data f( 2, 9),f( 2,10) / -1.14500000d-04 ,  3.77000000d-05 /
      data f( 2,11),f( 2,12) /  3.50600000d-04 ,  7.92700000d-04 /
      data f( 2,13),f( 2,14) /  1.43860000d-03 ,  2.57440000d-03 /
      data f( 2,15),f( 2,16) /  3.31760000d-03 ,  3.18160000d-03 /
      data f( 2,17),f( 2,18) /  3.63700000d-04 , -2.58200000d-04 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -6.30000000d-06 /
      data f( 3, 3),f( 3, 4) / -6.62000000d-05 , -1.56500000d-04 /
      data f( 3, 5),f( 3, 6) / -1.77500000d-04 , -1.31400000d-04 /
      data f( 3, 7),f( 3, 8) / -1.18200000d-04 , -1.37100000d-04 /
      data f( 3, 9),f( 3,10) / -9.58000000d-05 ,  1.99000000d-05 /
      data f( 3,11),f( 3,12) /  2.65700000d-04 ,  6.24200000d-04 /
      data f( 3,13),f( 3,14) /  1.11450000d-03 ,  1.96250000d-03 /
      data f( 3,15),f( 3,16) /  2.77120000d-03 ,  2.64000000d-03 /
      data f( 3,17),f( 3,18) /  6.68300000d-04 , -4.35000000d-05 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -3.40000000d-06 /
      data f( 4, 3),f( 4, 4) / -3.83000000d-05 , -9.27000000d-05 /
      data f( 4, 5),f( 4, 6) / -1.10100000d-04 , -9.39000000d-05 /
      data f( 4, 7),f( 4, 8) / -9.12000000d-05 , -9.72000000d-05 /
      data f( 4, 9),f( 4,10) / -7.22000000d-05 ,  3.50000000d-06 /
      data f( 4,11),f( 4,12) /  1.67600000d-04 ,  4.22500000d-04 /
      data f( 4,13),f( 4,14) /  7.48400000d-04 ,  1.20520000d-03 /
      data f( 4,15),f( 4,16) /  1.72650000d-03 ,  1.64480000d-03 /
      data f( 4,17),f( 4,18) /  7.05200000d-04 ,  6.57000000d-05 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -1.40000000d-06 /
      data f( 5, 3),f( 5, 4) / -1.48000000d-05 , -3.80000000d-05 /
      data f( 5, 5),f( 5, 6) / -5.06000000d-05 , -5.00000000d-05 /
      data f( 5, 7),f( 5, 8) / -4.71000000d-05 , -4.43000000d-05 /
      data f( 5, 9),f( 5,10) / -4.32000000d-05 , -7.50000000d-06 /
      data f( 5,11),f( 5,12) /  7.17000000d-05 ,  1.96000000d-04 /
      data f( 5,13),f( 5,14) /  3.71000000d-04 ,  5.39600000d-04 /
      data f( 5,15),f( 5,16) /  4.83500000d-04 ,  5.37300000d-04 /
      data f( 5,17),f( 5,18) /  3.41100000d-04 ,  4.78000000d-05 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -4.20000000d-06 /
      data f( 6, 3),f( 6, 4) / -3.65000000d-05 , -6.93000000d-05 /
      data f( 6, 5),f( 6, 6) / -5.25000000d-05 , -2.74000000d-05 /
      data f( 6, 7),f( 6, 8) / -3.02000000d-05 , -9.00000000d-06 /
      data f( 6, 9),f( 6,10) / -5.17000000d-05 , -4.66000000d-05 /
      data f( 6,11),f( 6,12) / -7.70000000d-05 , -4.81000000d-05 /
      data f( 6,13),f( 6,14) /  3.69000000d-05 ,  1.27500000d-04 /
      data f( 6,15),f( 6,16) /  1.98000000d-04 ,  2.17800000d-04 /
      data f( 6,17),f( 6,18) /  1.46300000d-04 ,  1.29300000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f( 7, 3),f( 7, 4) / -1.70000000d-06 , -4.90000000d-06 /
      data f( 7, 5),f( 7, 6) / -7.00000000d-06 , -5.80000000d-06 /
      data f( 7, 7),f( 7, 8) /  3.70000000d-06 ,  1.74000000d-05 /
      data f( 7, 9),f( 7,10) /  9.00000000d-06 ,  1.10000000d-06 /
      data f( 7,11),f( 7,12) /  1.48000000d-05 ,  4.31000000d-05 /
      data f( 7,13),f( 7,14) /  8.06000000d-05 ,  1.15600000d-04 /
      data f( 7,15),f( 7,16) /  1.27700000d-04 ,  1.21500000d-04 /
      data f( 7,17),f( 7,18) /  6.54000000d-05 ,  6.30000000d-06 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f( 8, 3),f( 8, 4) / -7.00000000d-07 , -1.90000000d-06 /
      data f( 8, 5),f( 8, 6) / -2.30000000d-06 ,  7.00000000d-07 /
      data f( 8, 7),f( 8, 8) /  1.67000000d-05 ,  1.78000000d-05 /
      data f( 8, 9),f( 8,10) /  1.25000000d-05 ,  4.90000000d-06 /
      data f( 8,11),f( 8,12) /  8.20000000d-06 ,  2.24000000d-05 /
      data f( 8,13),f( 8,14) /  3.84000000d-05 ,  5.36000000d-05 /
      data f( 8,15),f( 8,16) /  5.95000000d-05 ,  4.89000000d-05 /
      data f( 8,17),f( 8,18) /  2.07000000d-05 ,  1.49000000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 9, 3),f( 9, 4) / -3.00000000d-07 , -6.00000000d-07 /
      data f( 9, 5),f( 9, 6) /  0.00000000d+00 ,  2.50000000d-06 /
      data f( 9, 7),f( 9, 8) /  1.17000000d-05 ,  1.45000000d-05 /
      data f( 9, 9),f( 9,10) /  1.15000000d-05 ,  7.80000000d-06 /
      data f( 9,11),f( 9,12) /  5.60000000d-06 ,  9.70000000d-06 /
      data f( 9,13),f( 9,14) /  1.81000000d-05 ,  2.54000000d-05 /
      data f( 9,15),f( 9,16) /  2.65000000d-05 ,  1.91000000d-05 /
      data f( 9,17),f( 9,18) /  7.00000000d-06 ,  6.00000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(10, 3),f(10, 4) /  2.00000000d-07 ,  6.00000000d-07 /
      data f(10, 5),f(10, 6) /  1.40000000d-06 ,  3.50000000d-06 /
      data f(10, 7),f(10, 8) /  6.30000000d-06 ,  7.60000000d-06 /
      data f(10, 9),f(10,10) /  7.10000000d-06 ,  6.40000000d-06 /
      data f(10,11),f(10,12) /  3.00000000d-06 ,  3.40000000d-06 /
      data f(10,13),f(10,14) /  4.60000000d-06 ,  5.50000000d-06 /
      data f(10,15),f(10,16) /  5.10000000d-06 ,  3.00000000d-06 /
      data f(10,17),f(10,18) /  9.00000000d-07 ,  1.00000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(11, 3),f(11, 4) /  1.00000000d-07 ,  5.00000000d-07 /
      data f(11, 5),f(11, 6) /  1.00000000d-06 ,  1.80000000d-06 /
      data f(11, 7),f(11, 8) /  2.50000000d-06 ,  2.90000000d-06 /
      data f(11, 9),f(11,10) /  3.00000000d-06 ,  2.50000000d-06 /
      data f(11,11),f(11,12) /  2.10000000d-06 ,  1.70000000d-06 /
      data f(11,13),f(11,14) /  2.10000000d-06 ,  1.90000000d-06 /
      data f(11,15),f(11,16) /  1.50000000d-06 ,  6.00000000d-07 /
      data f(11,17),f(11,18) /  1.00000000d-07 ,  2.00000000d-07 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) /  1.00000000d-07 ,  2.00000000d-07 /
      data f(12, 5),f(12, 6) /  4.00000000d-07 ,  6.00000000d-07 /
      data f(12, 7),f(12, 8) /  8.00000000d-07 ,  1.00000000d-06 /
      data f(12, 9),f(12,10) /  1.10000000d-06 ,  1.20000000d-06 /
      data f(12,11),f(12,12) /  1.10000000d-06 ,  9.00000000d-07 /
      data f(12,13),f(12,14) /  9.00000000d-07 ,  9.00000000d-07 /
      data f(12,15),f(12,16) /  8.00000000d-07 ,  4.00000000d-07 /
      data f(12,17),f(12,18) /  1.00000000d-07 , -1.00000000d-07 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(13, 5),f(13, 6) /  2.00000000d-07 ,  3.00000000d-07 /
      data f(13, 7),f(13, 8) /  3.00000000d-07 ,  4.00000000d-07 /
      data f(13, 9),f(13,10) /  4.00000000d-07 ,  4.00000000d-07 /
      data f(13,11),f(13,12) /  5.00000000d-07 ,  4.00000000d-07 /
      data f(13,13),f(13,14) /  4.00000000d-07 ,  4.00000000d-07 /
      data f(13,15),f(13,16) /  3.00000000d-07 ,  2.00000000d-07 /
      data f(13,17),f(13,18) /  1.00000000d-07 ,  0.00000000d+00 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-1.14525620d-06/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 2.59468357d-04,-8.00487590d-07/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 1.81820681d-04,-1.33479343d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-6.54643534d-04, 1.20766133d-06/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-2.07753956d-03, 2.58814813d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-2.25771926d-03,-7.02538423d-08/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-8.47786985d-04,-2.27113276d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-8.39010892d-05, 1.49878488d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-1.20392497d-04, 1.24999324d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/ 3.81278134d-04, 2.00924218d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 3.77185341d-04, 1.59703806d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 8.34646952d-04, 9.80605574d-07/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.52609493d-03, 1.20005396d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-9.11347362d-03,-2.70107641d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/-2.82350680d-02, 3.40445169d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 7.98331856d-03,-7.87473035d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/ 3.84028167d-02, 5.76486971d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/ 3.27885917d-02,-5.46548488d-06/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-1.15767576d-05/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-9.93129069d-07/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 1.46563286d-04,-5.90741863d-07/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 7.32872088d-05,-8.49903481d-07/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-5.16498646d-04, 1.20635579d-06/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-1.45992088d-03, 1.52648034d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-1.50363292d-03,-6.10277132d-07/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-5.20783174d-04,-1.21537181d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-7.30549644d-05, 1.33176437d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-8.17864341d-05, 7.00314332d-07/
      data fpp( 2,10,1),fpp( 2,10,2)/ 2.98443731d-04, 1.89097830d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 4.27415032d-04, 1.37777245d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 8.17420383d-04, 3.49931886d-07/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.56988156d-03, 9.45050000d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/-4.52383847d-03,-8.75793190d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/-1.77850782d-02, 2.02522761d-06/
      data fpp( 2,16,1),fpp( 2,16,2)/ 4.71164859d-03,-5.20949785d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/ 2.18236523d-02, 4.54406865d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/ 1.92422451d-02, 2.09223242d-06/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-1.00361621d-06/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-8.07563190d-07/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-1.10721501d-04,-4.56873620d-07/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-2.19969517d-04,-5.80942331d-07/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-4.74361882d-04, 9.56642944d-07/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-6.32776911d-04, 9.12370557d-07/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-3.52749078d-04,-5.80125171d-07/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 3.59196799d-05,-5.65869874d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-4.38790533d-05, 9.17604668d-07/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-9.24617663d-05, 5.07451202d-07/
      data fpp( 3,10,1),fpp( 3,10,2)/ 7.49469403d-05, 1.51659052d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 2.68154531d-04, 1.23218670d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 5.45671518d-04, 3.16662681d-07/
      data fpp( 3,13,1),fpp( 3,13,2)/ 1.56937881d-03, 5.40916258d-06/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.92382750d-03,-4.91312993d-07/
      data fpp( 3,15,1),fpp( 3,15,2)/-1.67961915d-03,-5.80191061d-06/
      data fpp( 3,16,1),fpp( 3,16,2)/-5.82991293d-03,-3.26950446d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.23874260d-02, 2.61520889d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/-7.99757205d-03, 3.68068888d-06/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 4.44315556d-06/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-4.58474527d-07/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 1.46961448d-05,-2.68050945d-07/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-5.62641715d-06,-3.59321692d-07/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-7.11279641d-05, 5.35337712d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-1.34129709d-04, 4.37970844d-07/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-3.17477958d-05,-2.71221087d-07/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 5.74565163d-05,-1.63086494d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-6.50331795d-05, 4.01567063d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 6.60635103d-05, 4.16818241d-07/
      data fpp( 4,10,1),fpp( 4,10,2)/ 2.37881045d-04, 9.73159975d-07/
      data fpp( 4,11,1),fpp( 4,11,2)/ 7.71208210d-04, 9.94541861d-07/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.03948135d-03, 4.96672580d-07/
      data fpp( 4,13,1),fpp( 4,13,2)/ 1.72548291d-03, 1.27876782d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.97313397d-03, 2.24225614d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.44878264d-03,-6.37779239d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/ 4.10527737d-03,-1.29110866d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/-1.25768159d-03, 6.54813869d-06/
      data fpp( 4,18,1),fpp( 4,18,2)/-3.59589869d-04, 4.72453180d-06/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 8.98173410d-06/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-1.51027775d-07/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-4.85947630d-05,-1.01944449d-07/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-4.02013755d-04,-1.61194428d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-7.26973386d-04, 1.58722162d-07/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-4.59118785d-04, 1.62305781d-07/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-1.33157607d-04,-1.59452841d-08/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-2.27012660d-04, 3.94753560d-08/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-9.19663936d-05,-1.47956140d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-4.03926173d-04, 4.50349204d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/-4.14187508d-04, 4.22559325d-07/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.00635899d-03, 4.69413496d-07/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.02174324d-03, 4.05786690d-07/
      data fpp( 5,13,1),fpp( 5,13,2)/-8.76772613d-04, 9.49439742d-07/
      data fpp( 5,14,1),fpp( 5,14,2)/-1.50725219d-04,-4.58754566d-06/
      data fpp( 5,15,1),fpp( 5,15,2)/ 5.12766703d-03, 3.91874290d-06/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.58906017d-03,-4.49342593d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 1.24263669d-03,-9.45039172d-07/
      data fpp( 5,18,1),fpp( 5,18,2)/ 1.15163081d-03, 2.44758262d-06/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 5.88470869d-06/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-5.75476545d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 6.44829070d-05,-2.44046910d-07/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 5.28881438d-04,-1.34335814d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 9.15021508d-04, 7.51390166d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 4.97004849d-04, 1.04775150d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 5.31782224d-05,-6.72490765d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 1.97794124d-04, 9.11187909d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 1.04987541d-05,-1.53226087d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 6.49641182d-04, 1.38385558d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 7.44468985d-04,-1.13516145d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 1.98702775d-03, 1.02679022d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/ 2.62509160d-03, 5.86000569d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/ 2.82080754d-03,-4.79249681d-09/
      data fpp( 6,14,1),fpp( 6,14,2)/ 2.71376690d-03,-2.30830582d-07/
      data fpp( 6,15,1),fpp( 6,15,2)/ 2.05492234d-05,-2.77885176d-07/
      data fpp( 6,16,1),fpp( 6,16,2)/ 4.50481948d-04,-1.69962872d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/ 3.50334824d-04, 1.59840004d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.86133339d-03,-1.42397144d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-2.64051428d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-1.52975992d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-4.37368651d-05,-1.24048016d-08/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-3.57511996d-04,-2.50831943d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-6.36312645d-04, 1.67375787d-08/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-3.91300611d-04, 2.41328796d-08/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.03555283d-04, 8.47309031d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-1.56163836d-04, 1.34943508d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.63628623d-04,-3.72504935d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-5.33838555d-04, 2.90762328d-08/
      data fpp( 7,10,1),fpp( 7,10,2)/-4.80488433d-04, 2.86200004d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/-1.16975201d-03, 1.22123751d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.43142317d-03, 1.01304993d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/-1.33925754d-03, 2.46562782d-08/
      data fpp( 7,14,1),fpp( 7,14,2)/-1.09954238d-03,-3.49930106d-07/
      data fpp( 7,15,1),fpp( 7,15,2)/-4.50639270d-05, 1.06414433d-09/
      data fpp( 7,16,1),fpp( 7,16,2)/-3.41879634d-05,-7.52326472d-07/
      data fpp( 7,17,1),fpp( 7,17,2)/ 8.96240111d-05, 1.42417424d-08/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.38570275d-03, 5.15359502d-07/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 1.09232025d-06/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-5.07797524d-09/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.20645535d-05,-2.84404952d-09/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 8.99665459d-05,-1.35458267d-08/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 1.56629072d-04, 2.10273563d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 8.89975942d-05,-2.25635983d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.35708990d-06, 2.73227037d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-7.47387789d-05,-2.90344550d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 2.00157364d-05,-5.84883659d-09/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 1.12913038d-04,-7.02601035d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.23884748d-04, 1.48889251d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/ 3.30380305d-04, 1.28703101d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/ 4.15001073d-04,-9.70165463d-09/
      data fpp( 8,13,1),fpp( 8,13,2)/ 4.74622638d-04, 1.81035175d-08/
      data fpp( 8,14,1),fpp( 8,14,2)/ 4.82002630d-04,-1.10712416d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/ 2.10106485d-04,-1.33253856d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.55069905d-04,-3.46272162d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/ 1.59969131d-04, 4.62342505d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/-5.23077593d-04,-1.59097859d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-3.71951071d-07/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-6.51951977d-09/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-2.12134870d-06,-1.96096046d-09/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.67541875d-05,-3.63663837d-09/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-3.10036447d-05, 1.65075140d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-2.22897659d-05,-8.39341745d-09/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-3.81635726d-06, 1.31066156d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.31189518d-05,-1.13871206d-07/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-5.23432289d-06,-5.95813321d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-2.58135972d-05, 4.19653437d-09/
      data fpp( 9,10,1),fpp( 9,10,2)/-3.66505596d-05, 7.95194640d-10/
      data fpp( 9,11,1),fpp( 9,11,2)/-5.57692053d-05, 8.26226871d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/-3.65811242d-05, 4.67140571d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/-3.36330085d-05,-1.14789153d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-1.72681358d-05,-6.67983957d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/ 4.94379886d-05,-9.33275017d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 4.11083421d-05,-6.98915975d-08/
      data fpp( 9,17,1),fpp( 9,17,2)/ 1.44994631d-05, 9.08938916d-08/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.57007627d-04, 4.83160310d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 6.38419845d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00,-7.62130606d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/ 3.31769370d-07,-2.75738789d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/ 3.47928962d-06, 6.65085761d-09/
      data fpp(10, 4,1),fpp(10, 4,2)/ 6.29639788d-06, 1.53957462d-10/
      data fpp(10, 5,1),fpp(10, 5,2)/ 3.17050075d-06, 1.67333125d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-3.47238327d-06, 1.09127924d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-4.38746598d-06,-1.83844820d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 3.89510049d-06,-2.73748645d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/ 6.58427248d-06, 1.98839399d-08/
      data fpp(10,10,1),fpp(10,10,2)/ 4.80930468d-06,-6.41608953d-08/
      data fpp(10,11,1),fpp(10,11,2)/ 1.77174634d-05, 7.47596413d-08/
      data fpp(10,12,1),fpp(10,12,2)/ 1.68428361d-05,-6.87766973d-09/
      data fpp(10,13,1),fpp(10,13,2)/ 2.61877065d-05, 7.51037643d-10/
      data fpp(10,14,1),fpp(10,14,2)/ 2.98030927d-05,-1.41264808d-08/
      data fpp(10,15,1),fpp(10,15,2)/ 1.42327918d-05,-2.22451142d-08/
      data fpp(10,16,1),fpp(10,16,2)/ 1.01400211d-05, 1.10693784d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 4.31704510d-06, 1.78173629d-08/
      data fpp(10,18,1),fpp(10,18,2)/-4.08840832d-05, 5.62361060d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 1.68819470d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-1.04172873d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-4.05728782d-07,-4.16542541d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-7.62970952d-07, 9.07898894d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.98194681d-06,-2.15053034d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-1.19223704d-06, 5.52313242d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.50589033d-06,-1.94199934d-09/
      data fpp(11, 7,1),fpp(11, 7,2)/ 4.03091210d-06,-3.75513504d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 2.85392092d-06,-1.03746049d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/ 1.27650724d-06,-1.00950230d-08/
      data fpp(11,10,1),fpp(11,10,2)/ 2.41334087d-06, 5.41755251d-09/
      data fpp(11,11,1),fpp(11,11,2)/-4.90064826d-06,-5.57518705d-09/
      data fpp(11,12,1),fpp(11,12,2)/-3.19022002d-06, 1.68831957d-08/
      data fpp(11,13,1),fpp(11,13,2)/-5.11781750d-06,-1.39575957d-08/
      data fpp(11,14,1),fpp(11,14,2)/-4.14423492d-06, 2.94718700d-09/
      data fpp(11,15,1),fpp(11,15,2)/ 4.30844271d-07,-9.83115232d-09/
      data fpp(11,16,1),fpp(11,16,2)/ 5.31573385d-07, 6.37742227d-09/
      data fpp(11,17,1),fpp(11,17,2)/ 3.23565534d-08, 8.32146323d-09/
      data fpp(11,18,1),fpp(11,18,2)/ 1.01287061d-05,-3.66327521d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-1.16683624d-08/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 1.89228884d-09/
      data fpp(12, 2,1),fpp(12, 2,2)/ 9.11457563d-08, 1.21542232d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.72594190d-07,-7.53978112d-10/
      data fpp(12, 4,1),fpp(12, 4,2)/ 4.31389362d-07, 1.80049013d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/ 3.98447408d-07,-4.47982403d-10/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.48821933d-07,-8.56051532d-12/
      data fpp(12, 7,1),fpp(12, 7,2)/ 8.63817580d-07, 4.82224465d-10/
      data fpp(12, 8,1),fpp(12, 8,2)/ 1.48921582d-06,-1.92033734d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 1.50969855d-06, 1.19912491d-09/
      data fpp(12,10,1),fpp(12,10,2)/ 1.13733183d-06,-2.87616229d-09/
      data fpp(12,11,1),fpp(12,11,2)/ 1.28512965d-06,-1.69447575d-09/
      data fpp(12,12,1),fpp(12,12,2)/ 1.31804400d-06, 3.65406530d-09/
      data fpp(12,13,1),fpp(12,13,2)/ 2.08356350d-06,-9.21785446d-10/
      data fpp(12,14,1),fpp(12,14,2)/ 2.37384698d-06, 3.30764852d-11/
      data fpp(12,15,1),fpp(12,15,2)/ 1.44383115d-06,-5.21052049d-09/
      data fpp(12,16,1),fpp(12,16,2)/ 9.33685323d-07, 2.80900549d-09/
      data fpp(12,17,1),fpp(12,17,2)/ 3.53528689d-07,-2.55014790d-11/
      data fpp(12,18,1),fpp(12,18,2)/-2.03074122d-06, 3.29300042d-09/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 4.85349979d-09/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-8.99618502d-10/
      data fpp(13, 2,1),fpp(13, 2,2)/ 5.29427122d-07,-2.00762995d-10/
      data fpp(13, 3,1),fpp(13, 3,2)/-2.86297095d-07, 1.70267048d-09/
      data fpp(13, 4,1),fpp(13, 4,2)/ 4.46805319d-07,-6.09918941d-10/
      data fpp(13, 5,1),fpp(13, 5,2)/ 9.00776296d-07, 7.37005279d-10/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.05058903d-06,-2.33810218d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/-2.56908790d-07, 2.61540343d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.09460791d-06,-2.12351154d-09/
      data fpp(13, 9,1),fpp(13, 9,2)/-5.17349276d-07,-1.21357287d-10/
      data fpp(13,10,1),fpp(13,10,2)/-1.91866591d-06, 2.60894068d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 6.94935174d-07,-4.31440545d-09/
      data fpp(13,12,1),fpp(13,12,2)/-7.09022002d-07, 2.64868110d-09/
      data fpp(13,13,1),fpp(13,13,2)/-8.41781750d-07,-2.80318954d-10/
      data fpp(13,14,1),fpp(13,14,2)/-2.79942349d-06,-1.52740528d-09/
      data fpp(13,15,1),fpp(13,15,2)/-3.19691557d-06, 3.89940090d-10/
      data fpp(13,16,1),fpp(13,16,2)/-2.76684266d-06,-3.23550756d-11/
      data fpp(13,17,1),fpp(13,17,2)/-1.07676434d-06,-2.60519787d-10/
      data fpp(13,18,1),fpp(13,18,2)/ 2.07787061d-06, 1.07443422d-09/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 1.96278289d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 8.03941231d-07,-3.28735331d-06/
      data fpppp( 1, 3),fpppp( 1, 4)/-7.88148996d-06,-1.07156792d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.55582983d-05, 2.30454658d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/-1.23334432d-05,-1.24744754d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.42087066d-05,-1.20706286d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 3.72800227d-06, 2.48518837d-05/
      data fpppp( 1,13),fpppp( 1,14)/-8.90963549d-05,-3.48327456d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 9.73484628d-04,-2.25212194d-04/
      data fpppp( 1,17),fpppp( 1,18)/-4.20569162d-04,-2.54534545d-04/
      data fpppp( 1,19) /            -1.91754667d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 7.16201438d-07,-2.09417547d-06/
      data fpppp( 2, 3),fpppp( 2, 4)/-5.52986135d-06,-6.77696581d-06/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.14195417d-05, 1.50814112d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/-1.01514798d-05,-6.58278396d-06/
      data fpppp( 2, 9),fpppp( 2,10)/ 9.09503489d-06,-6.45965749d-06/
      data fpppp( 2,11),fpppp( 2,12)/ 1.66806320d-06, 1.54494477d-05/
      data fpppp( 2,13),fpppp( 2,14)/-4.17185041d-05,-2.59346304d-04/
      data fpppp( 2,15),fpppp( 2,16)/ 6.49052538d-04,-1.91385858d-04/
      data fpppp( 2,17),fpppp( 2,18)/-2.06592490d-04,-1.63848841d-04/
      data fpppp( 2,19) /            -1.37662415d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 1.34706438d-06, 2.82432455d-07/
      data fpppp( 3, 3),fpppp( 3, 4)/-2.38838512d-06, 5.62447096d-07/
      data fpppp( 3, 5),fpppp( 3, 6)/ 5.89723688d-06, 2.15517711d-06/
      data fpppp( 3, 7),fpppp( 3, 8)/-7.99948983d-06, 1.73473276d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 2.93352000d-06,-5.09327576d-07/
      data fpppp( 3,11),fpppp( 3,12)/ 6.51723332d-07, 2.96099806d-06/
      data fpppp( 3,13),fpppp( 3,14)/ 3.22757029d-05,-1.12219326d-04/
      data fpppp( 3,15),fpppp( 3,16)/ 5.91278812d-05,-9.71030259d-05/
      data fpppp( 3,17),fpppp( 3,18)/ 1.84851066d-04, 1.45407824d-05/
      data fpppp( 3,19) /            -2.65511080d-05 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.51692554d-07,-3.43970528d-07/
      data fpppp( 4, 3),fpppp( 4, 4)/-4.73547737d-07,-4.72577627d-07/
      data fpppp( 4, 5),fpppp( 4, 6)/ 2.51384638d-06, 3.40211598d-07/
      data fpppp( 4, 7),fpppp( 4, 8)/-4.66534883d-06, 5.61954324d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-2.59764100d-06, 7.21427146d-06/
      data fpppp( 4,11),fpppp( 4,12)/-4.56886698d-06,-4.84204489d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 4.90007517d-05,-9.74619921d-05/
      data fpppp( 4,15),fpppp( 4,16)/ 1.14527073d-04,-1.69795537d-04/
      data fpppp( 4,17),fpppp( 4,18)/ 1.43487852d-04,-2.84928316d-05/
      data fpppp( 4,19) /            -6.18266375d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/-6.58764692d-06,-2.63510485d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.16138742d-06, 8.98821621d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/ 7.77376465d-07,-8.61132743d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 8.47895933d-06,-1.15704307d-05/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.09824008d-05,-1.42572657d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 1.11320531d-05, 4.33628752d-06/
      data fpppp( 5,13),fpppp( 5,14)/-1.88559107d-05, 1.05951962d-04/
      data fpppp( 5,15),fpppp( 5,16)/-1.31811244d-04, 1.22730685d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 3.42499733d-05,-1.39479057d-05/
      data fpppp( 5,19) /            -4.20958466d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 9.03678772d-06, 3.48737497d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.00864983d-06,-1.22174819d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.88125858d-07, 1.22213873d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.31908716d-05, 2.06274229d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-1.97325522d-05, 2.56439083d-05/
      data fpppp( 6,11),fpppp( 6,12)/-1.39792231d-05,-5.99671076d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 1.14251913d-05,-5.78694492d-05/
      data fpppp( 6,15),fpppp( 6,16)/ 6.48819832d-05,-1.42694596d-05/
      data fpppp( 6,17),fpppp( 6,18)/-3.96089358d-05, 4.60139374d-05/
      data fpppp( 6,19) /             9.99332823d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/-5.92843878d-06,-2.34452535d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-8.95755759d-07, 8.02601728d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 2.20447626d-07,-6.34381019d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 4.73356027d-06,-9.88180488d-06/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.30289505d-05,-1.68203938d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 9.69580270d-06, 3.69272858d-06/
      data fpppp( 7,13),fpppp( 7,14)/-3.23651027d-06, 1.81062848d-05/
      data fpppp( 7,15),fpppp( 7,16)/-2.03028312d-05, 4.88890521d-07/
      data fpppp( 7,17),fpppp( 7,18)/ 2.51234298d-05,-3.06466040d-05/
      data fpppp( 7,19) /            -6.34439026d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.47794370d-06, 5.60773870d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.29207159d-07,-2.15197045d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.21034366d-07,-4.95559361d-07/
      data fpppp( 8, 7),fpppp( 8, 8)/ 2.67958278d-06,-1.34599522d-07/
      data fpppp( 8, 9),fpppp( 8,10)/-2.25261751d-06, 4.22953408d-06/
      data fpppp( 8,11),fpppp( 8,12)/-2.93408802d-06, 1.94330729d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 6.56812899d-07,-5.95607675d-06/
      data fpppp( 8,15),fpppp( 8,16)/ 6.41092591d-06,-6.76052935d-07/
      data fpppp( 8,17),fpppp( 8,18)/-1.21105658d-05, 1.38415593d-05/
      data fpppp( 8,19) /             2.91117878d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-2.63204046d-07,-1.06834046d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-6.01491774d-08, 3.70433653d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-4.37852802d-08, 3.90279263d-07/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.00961775d-06, 3.30876709d-07/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.52550937d-07,-3.56541747d-07/
      data fpppp( 9,11),fpppp( 9,12)/ 7.76715057d-07,-4.51914881d-07/
      data fpppp( 9,13),fpppp( 9,14)/ 5.65465448d-08, 1.03073412d-06/
      data fpppp( 9,15),fpppp( 9,16)/-1.15900794d-06,-8.96848645d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 3.64964857d-06,-3.55472309d-06/
      data fpppp( 9,19) /            -7.40170363d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 6.13281072d-08, 2.47395690d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 8.65866915d-09,-7.91989648d-08/
      data fpppp(10, 5),fpppp(10, 6)/-4.84431334d-08, 6.19522854d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 1.44302070d-07,-8.73016157d-08/
      data fpppp(10, 9),fpppp(10,10)/-1.30699276d-07, 3.42250333d-07/
      data fpppp(10,11),fpppp(10,12)/-3.57314468d-07, 2.60040376d-07/
      data fpppp(10,13),fpppp(10,14)/-6.96771717d-08,-3.25100746d-07/
      data fpppp(10,15),fpppp(10,16)/ 2.18938932d-07, 1.37996835d-07/
      data fpppp(10,17),fpppp(10,18)/-8.74738598d-07, 9.98268422d-07/
      data fpppp(10,19) /             2.04677760d-06 /
      data fpppp(11, 1),fpppp(11, 2)/ 7.72644765d-09, 4.20610902d-09/
      data fpppp(11, 3),fpppp(11, 4)/-2.16416871d-08, 3.06566179d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 1.95363530d-08, 5.70302602d-09/
      data fpppp(11, 7),fpppp(11, 8)/-5.27347935d-08,-1.68846284d-08/
      data fpppp(11, 9),fpppp(11,10)/ 9.62479569d-08,-2.05252360d-07/
      data fpppp(11,11),fpppp(11,12)/ 2.17712119d-07,-1.24131072d-07/
      data fpppp(11,13),fpppp(11,14)/ 6.05306241d-08, 5.60793794d-08/
      data fpppp(11,15),fpppp(11,16)/-6.87583456d-08,-4.95070014d-08/
      data fpppp(11,17),fpppp(11,18)/ 2.30789594d-07,-2.37917394d-07/
      data fpppp(11,19) /            -4.92623358d-07 /
      data fpppp(12, 1),fpppp(12, 2)/-1.60055362d-09,-8.30693634d-10/
      data fpppp(12, 3),fpppp(12, 4)/ 4.34148881d-09,-5.89445737d-09/
      data fpppp(12, 5),fpppp(12, 6)/ 1.73211311d-09, 3.96499370d-09/
      data fpppp(12, 7),fpppp(12, 8)/ 4.28517935d-09,-8.48155582d-09/
      data fpppp(12, 9),fpppp(12,10)/-6.65388599d-09, 1.15261320d-08/
      data fpppp(12,11),fpppp(12,12)/-8.24076896d-09, 1.45439354d-08/
      data fpppp(12,13),fpppp(12,14)/-5.97866373d-09,-1.91434413d-08/
      data fpppp(12,15),fpppp(12,16)/ 9.33446981d-09, 6.99776298d-09/
      data fpppp(12,17),fpppp(12,18)/-4.15261704d-08, 5.08601220d-08/
      data fpppp(12,19) /             1.02986350d-07 /
      data fpppp(13, 1),fpppp(13, 2)/-4.01781524d-08,-1.78777948d-08/
      data fpppp(13, 3),fpppp(13, 4)/ 3.09802511d-08,-1.31136119d-08/
      data fpppp(13, 5),fpppp(13, 6)/ 4.72631015d-09,-2.40411231d-08/
      data fpppp(13, 7),fpppp(13, 8)/ 3.99954862d-09, 3.62308510d-08/
      data fpppp(13, 9),fpppp(13,10)/-6.40254875d-08, 1.01156583d-07/
      data fpppp(13,11),fpppp(13,12)/-9.97057808d-08, 5.66130446d-08/
      data fpppp(13,13),fpppp(13,14)/-5.04745518d-08, 3.57922432d-08/
      data fpppp(13,15),fpppp(13,16)/ 9.14558791d-10, 1.02034212d-08/
      data fpppp(13,17),fpppp(13,18)/ 3.38720806d-08,-5.78183454d-08/
      data fpppp(13,19) /            -1.16549033d-07 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine c5_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -3.10000000d-06 /
      data f( 1, 3),f( 1, 4) / -5.05000000d-05 , -1.69600000d-04 /
      data f( 1, 5),f( 1, 6) / -2.51600000d-04 , -1.92900000d-04 /
      data f( 1, 7),f( 1, 8) / -9.69000000d-05 , -6.20000000d-05 /
      data f( 1, 9),f( 1,10) /  2.00000000d-07 ,  5.54000000d-05 /
      data f( 1,11),f( 1,12) /  3.36700000d-04 ,  8.64700000d-04 /
      data f( 1,13),f( 1,14) /  2.08730000d-03 ,  5.12840000d-03 /
      data f( 1,15),f( 1,16) /  9.53990000d-03 ,  9.43350000d-03 /
      data f( 1,17),f( 1,18) /  2.12300000d-03 , -2.25400000d-04 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -3.40000000d-06 /
      data f( 2, 3),f( 2, 4) / -3.74000000d-05 , -1.19400000d-04 /
      data f( 2, 5),f( 2, 6) / -1.70000000d-04 , -1.28500000d-04 /
      data f( 2, 7),f( 2, 8) / -7.27000000d-05 , -5.01000000d-05 /
      data f( 2, 9),f( 2,10) /  3.90000000d-06 ,  6.19000000d-05 /
      data f( 2,11),f( 2,12) /  2.76100000d-04 ,  6.86500000d-04 /
      data f( 2,13),f( 2,14) /  1.55690000d-03 ,  3.73260000d-03 /
      data f( 2,15),f( 2,16) /  6.14470000d-03 ,  5.76330000d-03 /
      data f( 2,17),f( 2,18) /  1.92160000d-03 , -3.31900000d-04 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.40000000d-06 /
      data f( 3, 3),f( 3, 4) / -2.52000000d-05 , -8.33000000d-05 /
      data f( 3, 5),f( 3, 6) / -1.19100000d-04 , -9.39000000d-05 /
      data f( 3, 7),f( 3, 8) / -5.98000000d-05 , -4.13000000d-05 /
      data f( 3, 9),f( 3,10) /  4.70000000d-06 ,  5.87000000d-05 /
      data f( 3,11),f( 3,12) /  2.23800000d-04 ,  5.40100000d-04 /
      data f( 3,13),f( 3,14) /  1.14920000d-03 ,  2.58980000d-03 /
      data f( 3,15),f( 3,16) /  4.26220000d-03 ,  4.14080000d-03 /
      data f( 3,17),f( 3,18) /  1.81510000d-03 , -2.83000000d-05 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -6.00000000d-07 /
      data f( 4, 3),f( 4, 4) / -1.36000000d-05 , -4.76000000d-05 /
      data f( 4, 5),f( 4, 6) / -7.18000000d-05 , -6.23000000d-05 /
      data f( 4, 7),f( 4, 8) / -4.58000000d-05 , -3.23000000d-05 /
      data f( 4, 9),f( 4,10) /  3.00000000d-06 ,  5.34000000d-05 /
      data f( 4,11),f( 4,12) /  1.63900000d-04 ,  3.73400000d-04 /
      data f( 4,13),f( 4,14) /  7.34800000d-04 ,  1.44150000d-03 /
      data f( 4,15),f( 4,16) /  2.37630000d-03 ,  2.40120000d-03 /
      data f( 4,17),f( 4,18) /  1.19060000d-03 ,  8.50000000d-05 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -3.00000000d-07 /
      data f( 5, 3),f( 5, 4) / -5.10000000d-06 , -1.92000000d-05 /
      data f( 5, 5),f( 5, 6) / -3.21000000d-05 , -3.33000000d-05 /
      data f( 5, 7),f( 5, 8) / -2.99000000d-05 , -2.24000000d-05 /
      data f( 5, 9),f( 5,10) / -2.50000000d-06 ,  3.46000000d-05 /
      data f( 5,11),f( 5,12) /  1.01700000d-04 ,  2.08900000d-04 /
      data f( 5,13),f( 5,14) /  3.66200000d-04 ,  5.98400000d-04 /
      data f( 5,15),f( 5,16) /  7.98300000d-04 ,  8.37500000d-04 /
      data f( 5,17),f( 5,18) /  4.46000000d-04 ,  3.86000000d-05 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -4.20000000d-06 /
      data f( 6, 3),f( 6, 4) / -1.20000000d-05 , -1.69000000d-05 /
      data f( 6, 5),f( 6, 6) / -1.85000000d-05 , -1.98000000d-05 /
      data f( 6, 7),f( 6, 8) / -2.92000000d-05 , -1.78000000d-05 /
      data f( 6, 9),f( 6,10) / -2.30000000d-05 ,  1.76300000d-04 /
      data f( 6,11),f( 6,12) /  2.64300000d-04 ,  3.39000000d-04 /
      data f( 6,13),f( 6,14) /  4.21900000d-04 ,  5.11700000d-04 /
      data f( 6,15),f( 6,16) /  6.14300000d-04 ,  6.35800000d-04 /
      data f( 6,17),f( 6,18) /  4.92300000d-04 ,  4.38200000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 7, 3),f( 7, 4) / -4.00000000d-07 , -2.00000000d-06 /
      data f( 7, 5),f( 7, 6) / -5.00000000d-06 , -9.10000000d-06 /
      data f( 7, 7),f( 7, 8) / -9.30000000d-06 ,  3.80000000d-06 /
      data f( 7, 9),f( 7,10) /  7.00000000d-07 ,  1.21000000d-05 /
      data f( 7,11),f( 7,12) /  3.38000000d-05 ,  6.51000000d-05 /
      data f( 7,13),f( 7,14) /  9.86000000d-05 ,  1.31600000d-04 /
      data f( 7,15),f( 7,16) /  1.55200000d-04 ,  1.36500000d-04 /
      data f( 7,17),f( 7,18) /  5.47000000d-05 ,  2.60000000d-06 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  1.00000000d-07 /
      data f( 8, 3),f( 8, 4) /  1.00000000d-07 , -5.00000000d-07 /
      data f( 8, 5),f( 8, 6) / -2.40000000d-06 , -4.20000000d-06 /
      data f( 8, 7),f( 8, 8) /  5.80000000d-06 ,  3.70000000d-06 /
      data f( 8, 9),f( 8,10) /  1.80000000d-06 ,  4.80000000d-06 /
      data f( 8,11),f( 8,12) /  1.51000000d-05 ,  3.60000000d-05 /
      data f( 8,13),f( 8,14) /  5.11000000d-05 ,  6.34000000d-05 /
      data f( 8,15),f( 8,16) /  6.32000000d-05 ,  4.43000000d-05 /
      data f( 8,17),f( 8,18) /  1.40000000d-05 , -5.20000000d-06 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  1.00000000d-07 /
      data f( 9, 3),f( 9, 4) / -1.00000000d-07 , -1.00000000d-07 /
      data f( 9, 5),f( 9, 6) /  2.00000000d-07 ,  5.20000000d-06 /
      data f( 9, 7),f( 9, 8) /  4.50000000d-06 ,  2.40000000d-06 /
      data f( 9, 9),f( 9,10) /  1.20000000d-06 ,  9.00000000d-07 /
      data f( 9,11),f( 9,12) /  4.70000000d-06 ,  1.32000000d-05 /
      data f( 9,13),f( 9,14) /  2.40000000d-05 ,  3.07000000d-05 /
      data f( 9,15),f( 9,16) /  2.76000000d-05 ,  1.58000000d-05 /
      data f( 9,17),f( 9,18) /  4.00000000d-06 ,  1.00000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(10, 3),f(10, 4) /  0.00000000d+00 ,  3.00000000d-07 /
      data f(10, 5),f(10, 6) /  1.30000000d-06 ,  2.20000000d-06 /
      data f(10, 7),f(10, 8) /  2.10000000d-06 ,  1.20000000d-06 /
      data f(10, 9),f(10,10) / -3.00000000d-07 , -2.00000000d-06 /
      data f(10,11),f(10,12) / -2.90000000d-06 , -3.00000000d-07 /
      data f(10,13),f(10,14) /  3.50000000d-06 ,  5.70000000d-06 /
      data f(10,15),f(10,16) /  4.60000000d-06 ,  2.10000000d-06 /
      data f(10,17),f(10,18) /  4.00000000d-07 ,  2.00000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(11, 3),f(11, 4) /  1.00000000d-07 ,  2.00000000d-07 /
      data f(11, 5),f(11, 6) /  7.00000000d-07 ,  1.00000000d-06 /
      data f(11, 7),f(11, 8) /  9.00000000d-07 ,  4.00000000d-07 /
      data f(11, 9),f(11,10) / -5.00000000d-07 , -1.50000000d-06 /
      data f(11,11),f(11,12) / -2.50000000d-06 , -2.30000000d-06 /
      data f(11,13),f(11,14) / -4.00000000d-07 ,  6.00000000d-07 /
      data f(11,15),f(11,16) /  6.00000000d-07 ,  1.00000000d-07 /
      data f(11,17),f(11,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(12, 5),f(12, 6) /  2.00000000d-07 ,  2.00000000d-07 /
      data f(12, 7),f(12, 8) /  0.00000000d+00 , -2.00000000d-07 /
      data f(12, 9),f(12,10) / -5.00000000d-07 , -8.00000000d-07 /
      data f(12,11),f(12,12) / -1.10000000d-06 , -1.30000000d-06 /
      data f(12,13),f(12,14) / -1.20000000d-06 , -4.00000000d-07 /
      data f(12,15),f(12,16) / -4.00000000d-07 , -1.00000000d-07 /
      data f(12,17),f(12,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(13, 5),f(13, 6) /  1.00000000d-07 ,  1.00000000d-07 /
      data f(13, 7),f(13, 8) /  0.00000000d+00 , -1.00000000d-07 /
      data f(13, 9),f(13,10) / -2.00000000d-07 , -3.00000000d-07 /
      data f(13,11),f(13,12) / -4.00000000d-07 , -4.00000000d-07 /
      data f(13,13),f(13,14) / -3.00000000d-07 , -3.00000000d-07 /
      data f(13,15),f(13,16) / -2.00000000d-07 , -1.00000000d-07 /
      data f(13,17),f(13,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-2.20356278d-07/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.25209737d-04,-3.40287444d-07/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 2.79921367d-05,-1.07649395d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-4.43293848d-04, 3.44263231d-07/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.08274933d-03, 1.92544102d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-1.13395319d-03, 3.95972681d-07/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-4.55459415d-04,-1.27133175d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-9.72287590d-05, 1.02335431d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-1.03903310d-04,-1.18408548d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/-4.01249917d-04, 3.29298762d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 2.11829883d-04, 1.57813499d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 9.08030454d-04, 5.19647241d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.50381405d-03, 1.93119754d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/ 5.45223759d-03, 2.66656261d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 5.84123663d-02,-4.37504799d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 8.50352332d-02,-1.22737706d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 8.87383354d-03, 1.02455306d-04/
      data fpp( 1,18,1),fpp( 1,18,2)/ 2.20491057d-02, 1.06424841d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 9.40275794d-06/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-1.70162942d-07/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 7.09376688d-05,-2.29674116d-07/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-8.91284484d-06,-7.47140594d-07/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-3.53483732d-04, 3.38236492d-07/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-7.88715619d-04, 1.27819463d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-7.80236485d-04, 7.49850061d-08/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-2.98581171d-04,-7.20134650d-07/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-7.48996250d-05, 8.13553593d-07/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-6.76933795d-05,-6.50079723d-07/
      data fpp( 2,10,1),fpp( 2,10,2)/-3.29285880d-04, 2.02676530d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.30554519d-04, 1.91501853d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 7.22653377d-04, 2.08516058d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.94715762d-03, 1.73443392d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 5.97852482d-03, 6.85548277d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/ 4.00211960d-02,-3.05822702d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 5.53993194d-02,-5.21364018d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/ 4.05983292d-03, 3.15098774d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/ 1.20544315d-02, 2.13888922d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 3.80585539d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.32661360d-07/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-6.39604122d-05,-1.68677280d-07/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-1.27340757d-04,-5.36629519d-07/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-2.57771224d-04, 2.57195355d-07/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-3.67388189d-04, 8.45848100d-07/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-2.15100873d-04, 1.94122462d-08/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-4.52159028d-05,-3.89497085d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-6.81727412d-05, 6.02576092d-07/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-6.03231719d-05,-3.70807283d-07/
      data fpp( 3,10,1),fpp( 3,10,2)/ 2.63393438d-04, 1.36065304d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 5.10952040d-04, 1.59419512d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 9.71356038d-04, 1.33456647d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/ 3.11255547d-03, 1.06355390d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 8.58366315d-03, 6.01327748d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 8.40784964d-03,-2.07806489d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 5.22489264d-04,-3.05186817d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.08781652d-02, 1.05973758d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/-8.75183151d-03, 1.70671783d-05/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 3.34359108d-05/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-5.35596448d-08/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 1.92429281d-05,-9.28807104d-08/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.62555789d-05,-3.18917514d-07/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.35106766d-04, 1.08550765d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-1.86228956d-04, 4.72714455d-07/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-1.16172766d-04, 2.25914159d-08/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-6.89287702d-06,-1.43080119d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-2.82444592d-06, 3.69729058d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 5.28728259d-05,-2.78361142d-08/
      data fpp( 4,10,1),fpp( 4,10,2)/-6.91787541d-04, 6.47615399d-07/
      data fpp( 4,11,1),fpp( 4,11,2)/-5.53543145d-04, 1.04337452d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/-1.92955711d-04, 1.11888652d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 8.03376682d-04, 3.59507940d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/ 5.12877297d-03, 5.21879590d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/ 7.81637052d-03,-1.07842630d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/ 7.60215619d-03,-1.66757439d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 2.57066208d-03, 3.35723875d-06/
      data fpp( 4,18,1),fpp( 4,18,2)/-1.67018260d-03, 9.54678893d-06/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 1.96916055d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-2.78769355d-09/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-4.80011226d-05,-3.34246129d-08/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-1.31577693d-04,-1.33513855d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-1.59395614d-04, 9.48003231d-09/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-1.22834429d-04, 1.67593726d-07/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-6.71866238d-05, 2.21450653d-08/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-1.29213252d-04, 1.98260133d-08/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-7.24581283d-05, 1.44550881d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-1.96999140d-04, 1.45970461d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/ 1.81648407d-03, 3.03567275d-07/
      data fpp( 5,11,1),fpp( 5,11,2)/ 2.36796684d-03, 4.39760439d-07/
      data fpp( 5,12,1),fpp( 5,12,2)/ 2.75464465d-03, 3.43390969d-07/
      data fpp( 5,13,1),fpp( 5,13,2)/ 3.29126133d-03, 1.19267568d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 4.13532862d-03,-6.20093707d-07/
      data fpp( 5,15,1),fpp( 5,15,2)/ 7.50690455d-03,-6.50300858d-07/
      data fpp( 5,16,1),fpp( 5,16,2)/ 7.41480663d-03,-6.42070286d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 5.41038047d-03, 4.91112309d-07/
      data fpp( 5,18,1),fpp( 5,18,2)/ 4.95008322d-03, 3.50225363d-06/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 7.62787319d-06/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-9.80011141d-08/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 7.19615625d-05,-4.19977717d-08/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.72966351d-04, 4.99922010d-08/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.46289223d-04, 1.60289677d-08/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 5.11666697d-05, 8.38919281d-08/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.29192614d-05,-3.33596680d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 1.58945885d-04, 7.64494793d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 1.65456959d-04,-1.47638249d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 3.75123733d-04, 4.14503517d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-2.72214873d-03,-2.83375818d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-3.52312421d-03, 5.11997561d-07/
      data fpp( 6,12,1),fpp( 6,12,2)/-3.75522290d-03,-1.22320616d-08/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.78522202d-03, 2.89306853d-08/
      data fpp( 6,14,1),fpp( 6,14,2)/-3.51648746d-03, 3.10509320d-07/
      data fpp( 6,15,1),fpp( 6,15,2)/-4.38798872d-03,-5.02967967d-07/
      data fpp( 6,16,1),fpp( 6,16,2)/-4.57338270d-03,-3.16463745d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/-5.23058397d-03, 3.26151778d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/-7.42615027d-03,-4.51743365d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-8.23778317d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 3.82520749d-09/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-4.54451272d-05,-3.65041498d-09/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-1.16287711d-04,-1.32235476d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-1.23361279d-04,-1.54553948d-08/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-8.42322503d-05,-8.95487341d-09/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-5.16904220d-05,-1.47251116d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-4.57702863d-05, 3.01855320d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.81369709d-04,-3.94696168d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-2.42695791d-04, 3.04929352d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.73051084d-03, 4.49787602d-08/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.29013000d-03, 1.33155607d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 2.57024696d-03,-1.60118889d-09/
      data fpp( 7,13,1),fpp( 7,13,2)/ 2.75362674d-03, 5.24914839d-09/
      data fpp( 7,14,1),fpp( 7,14,2)/ 2.88902121d-03,-4.93954047d-08/
      data fpp( 7,15,1),fpp( 7,15,2)/ 3.44265034d-03,-3.71667530d-07/
      data fpp( 7,16,1),fpp( 7,16,2)/ 3.73632418d-03,-1.00193448d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 3.89835539d-03, 5.93405436d-07/
      data fpp( 7,18,1),fpp( 7,18,2)/ 4.70971786d-03, 4.10312733d-07/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 7.35343634d-07/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 3.64469291d-09/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.14189462d-05,-2.89385822d-10/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 2.57844914d-05,-8.48714962d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 2.55558921d-05,-1.76201568d-09/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 2.41623316d-05,-6.24647877d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 5.46424266d-05, 2.57621166d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-9.10647393d-05,-2.60019878d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 3.92218750d-05, 5.64583438d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 5.32594307d-05, 4.61865023d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/-4.34294636d-04, 5.27956468d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-5.54195808d-04, 1.80630910d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-6.50564936d-04,-1.39319288d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-6.10084925d-04, 2.86462433d-08/
      data fpp( 8,14,1),fpp( 8,14,2)/-5.53997387d-04,-1.43265685d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/-5.72212649d-04,-2.05583504d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/-6.01514033d-04,-1.56400298d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/-8.37237605d-04, 1.47184696d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/-1.14552115d-03, 2.33661516d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 3.82169242d-07/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-6.79652177d-09/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-2.63065768d-06,-5.40695646d-09/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-3.65025486d-06, 1.04243476d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-5.26228953d-06,-2.42904340d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.24170761d-05, 1.04737388d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-5.88792842d-05,-1.12659119d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.64292435d-05, 3.89908829d-09/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-4.31779157d-06, 1.30627660d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-1.11419319d-05,-2.15015217d-09/
      data fpp( 9,10,1),fpp( 9,10,2)/ 8.82677040d-05, 4.95378427d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.25853228d-04, 4.99987813d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.83212785d-04, 3.24670322d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.76312963d-04,-4.18669101d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.78968336d-04,-1.10999392d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/ 1.99800253d-04,-1.02135522d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 1.98531949d-04,-2.45851987d-09/
      data fpp( 9,17,1),fpp( 9,17,2)/ 1.87395026d-04, 1.11969602d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.86766757d-04, 2.85801139d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 1.70994307d-09/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 1.66024342d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/ 9.82499943d-07, 2.67951316d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/ 1.05851891d-06,-3.78296075d-10/
      data fpp(10, 4,1),fpp(10, 4,2)/ 6.08922551d-07, 1.08336711d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/ 5.70062631d-07,-9.56388470d-10/
      data fpp(10, 6,1),fpp(10, 6,2)/ 1.85166394d-05,-1.30081173d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-2.55536095d-06,-7.01114251d-09/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.74243718d-06,-6.94731271d-09/
      data fpp(10, 9,1),fpp(10, 9,2)/ 4.99608037d-06,-1.19960665d-09/
      data fpp(10,10,1),fpp(10,10,2)/-1.82557940d-05,-2.54260698d-10/
      data fpp(10,11,1),fpp(10,11,2)/-2.12617793d-05, 5.02166494d-08/
      data fpp(10,12,1),fpp(10,12,2)/-3.17558860d-05, 9.38766294d-09/
      data fpp(10,13,1),fpp(10,13,2)/-2.16964267d-05,-1.57673012d-08/
      data fpp(10,14,1),fpp(10,14,2)/-1.75063153d-05,-4.23184582d-08/
      data fpp(10,15,1),fpp(10,15,2)/-2.40944338d-05,-1.29588662d-08/
      data fpp(10,16,1),fpp(10,16,2)/-3.50388315d-05, 1.01539229d-08/
      data fpp(10,17,1),fpp(10,17,2)/-4.51662745d-05, 2.03431746d-08/
      data fpp(10,18,1),fpp(10,18,2)/-5.05396956d-05,-1.52662132d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-1.42366893d-08/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 6.49469083d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/-9.93420904d-08, 4.01061834d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-5.83820765d-07,-4.53716418d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.73400671d-07, 8.13803837d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-6.31743766d-08,-4.01498931d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/-4.38727354d-06,-4.07808112d-09/
      data fpp(11, 7,1),fpp(11, 7,2)/ 9.92200249d-07,-3.67268622d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/-2.51957152d-07,-5.23117400d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.04238957d-06, 5.97382229d-10/
      data fpp(11,10,1),fpp(11,10,2)/ 5.15547212d-06,-3.15835491d-09/
      data fpp(11,11,1),fpp(11,11,2)/ 7.19388928d-06, 1.20360374d-08/
      data fpp(11,12,1),fpp(11,12,2)/ 1.28107595d-05, 2.70142052d-08/
      data fpp(11,13,1),fpp(11,13,2)/ 1.00727439d-05,-1.80928582d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 1.04569251d-05,-8.64277226d-09/
      data fpp(11,15,1),fpp(11,15,2)/ 1.05774826d-05,-7.33605274d-09/
      data fpp(11,16,1),fpp(11,16,2)/ 1.18233767d-05, 7.98698321d-09/
      data fpp(11,17,1),fpp(11,17,2)/ 1.24700722d-05,-6.11880094d-10/
      data fpp(11,18,1),fpp(11,18,2)/ 1.35920252d-05, 4.60537170d-10/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-1.23026858d-09/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-9.13927264d-10/
      data fpp(12, 2,1),fpp(12, 2,2)/ 1.48684181d-08,-1.72145473d-10/
      data fpp(12, 3,1),fpp(12, 3,2)/ 7.67641530d-08, 1.60250915d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/ 8.46801343d-08,-2.37891147d-10/
      data fpp(12, 5,1),fpp(12, 5,2)/ 2.82634875d-07,-6.50944567d-10/
      data fpp(12, 6,1),fpp(12, 6,2)/ 1.43245471d-06,-3.15833059d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/ 3.86559950d-07, 1.28426691d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/ 4.65391430d-07,-1.97873705d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 3.73477914d-07, 6.30681304d-10/
      data fpp(12,10,1),fpp(12,10,2)/-1.16609442d-06,-5.43988161d-10/
      data fpp(12,11,1),fpp(12,11,2)/-1.51377786d-06, 1.54527134d-09/
      data fpp(12,12,1),fpp(12,12,2)/-1.48715190d-06, 3.62902795d-10/
      data fpp(12,13,1),fpp(12,13,2)/ 5.45122520d-09, 1.50031175d-08/
      data fpp(12,14,1),fpp(12,14,2)/ 2.78614982d-07,-1.83753727d-08/
      data fpp(12,15,1),fpp(12,15,2)/-2.15496518d-07, 1.04983733d-08/
      data fpp(12,16,1),fpp(12,16,2)/-1.45467534d-06,-5.61812069d-09/
      data fpp(12,17,1),fpp(12,17,2)/-2.31401444d-06,-2.58905830d-11/
      data fpp(12,18,1),fpp(12,18,2)/-2.62840503d-06,-2.78316976d-10/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 1.13915849d-09/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-8.34558325d-10/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.94934209d-07,-3.30883350d-10/
      data fpp(13, 3,1),fpp(13, 3,2)/ 3.61617924d-07, 2.15809172d-09/
      data fpp(13, 4,1),fpp(13, 4,2)/ 1.32659933d-07,-2.30148355d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/ 5.33682562d-07, 1.04784247d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.46272646d-07,-1.88988634d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 1.04422002d-06, 5.11702880d-10/
      data fpp(13, 8,1),fpp(13, 8,2)/ 6.79804285d-07,-1.56925183d-10/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.49238957d-07, 1.15997851d-10/
      data fpp(13,10,1),fpp(13,10,2)/-4.29452788d-07,-3.07066223d-10/
      data fpp(13,11,1),fpp(13,11,2)/-2.20561107d-06, 1.11226704d-09/
      data fpp(13,12,1),fpp(13,12,2)/-3.59392405d-06, 1.85799806d-09/
      data fpp(13,13,1),fpp(13,13,2)/-1.30272561d-06,-2.54425930d-09/
      data fpp(13,14,1),fpp(13,14,2)/-2.91430749d-06, 2.31903913d-09/
      data fpp(13,15,1),fpp(13,15,2)/-1.34225174d-06,-7.31897226d-10/
      data fpp(13,16,1),fpp(13,16,2)/-9.47662328d-07, 6.08549772d-10/
      data fpp(13,17,1),fpp(13,17,2)/ 7.07007222d-07,-1.70230186d-09/
      data fpp(13,18,1),fpp(13,18,2)/ 1.08920252d-06, 2.00657675d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 8.99671163d-10/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-7.88565105d-07,-2.06286897d-06/
      data fpppp( 1, 3),fpppp( 1, 4)/-4.30559927d-06,-3.15883702d-06/
      data fpppp( 1, 5),fpppp( 1, 6)/ 6.85077735d-06, 1.10508256d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/-7.27222213d-06,-1.17772398d-06/
      data fpppp( 1, 9),fpppp( 1,10)/-9.91119439d-06, 2.33821782d-05/
      data fpppp( 1,11),fpppp( 1,12)/-2.89919341d-05, 9.75728044d-05/
      data fpppp( 1,13),fpppp( 1,14)/-2.47324302d-04, 8.52882802d-04/
      data fpppp( 1,15),fpppp( 1,16)/-1.03504598d-04,-2.01910012d-03/
      data fpppp( 1,17),fpppp( 1,18)/ 2.01284908d-03,-6.72095881d-04/
      data fpppp( 1,19) /            -1.43792822d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-4.40954312d-07,-1.36309303d-06/
      data fpppp( 2, 3),fpppp( 2, 4)/-3.15396452d-06,-1.90427128d-06/
      data fpppp( 2, 5),fpppp( 2, 6)/ 5.33138965d-06, 7.20137398d-06/
      data fpppp( 2, 7),fpppp( 2, 8)/-5.74631472d-06, 3.05458783d-07/
      data fpppp( 2, 9),fpppp( 2,10)/-8.46403842d-06, 1.74227701d-05/
      data fpppp( 2,11),fpppp( 2,12)/-1.79410680d-05, 6.22770094d-05/
      data fpppp( 2,13),fpppp( 2,14)/-1.33222646d-04, 5.19025353d-04/
      data fpppp( 2,15),fpppp( 2,16)/-8.22005253d-05,-1.31009612d-03/
      data fpppp( 2,17),fpppp( 2,18)/ 1.31952842d-03,-4.07972463d-04/
      data fpppp( 2,19) /            -8.90580368d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 6.48619265d-07, 7.27672470d-08/
      data fpppp( 3, 3),fpppp( 3, 4)/-9.04884231d-07,-4.76237599d-07/
      data fpppp( 3, 5),fpppp( 3, 6)/ 4.05864467d-06,-4.40841631d-08/
      data fpppp( 3, 7),fpppp( 3, 8)/-2.82644879d-06,-2.20629202d-07/
      data fpppp( 3, 9),fpppp( 3,10)/ 5.55735007d-06,-3.05674863d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 2.10016392d-06, 7.42681680d-06/
      data fpppp( 3,13),fpppp( 3,14)/ 6.90402950d-05,-8.37935021d-05/
      data fpppp( 3,15),fpppp( 3,16)/-7.26815576d-05,-8.80530795d-05/
      data fpppp( 3,17),fpppp( 3,18)/ 2.13976230d-04, 4.37674501d-05/
      data fpppp( 3,19) /             8.48383868d-06 /
      data fpppp( 4, 1),fpppp( 4, 2)/-3.44805230d-07,-3.80407697d-07/
      data fpppp( 4, 3),fpppp( 4, 4)/-1.41805009d-06, 1.05144726d-06/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.27600090d-06, 1.11525187d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/-3.38358635d-06, 6.10640603d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-1.79443073d-05, 1.76493651d-05/
      data fpppp( 4,11),fpppp( 4,12)/ 3.21132831d-07,-5.59331416d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 6.01968215d-05,-3.54501383d-05/
      data fpppp( 4,15),fpppp( 4,16)/-1.66641918d-05,-7.20018077d-05/
      data fpppp( 4,17),fpppp( 4,18)/ 1.56346354d-05, 5.69022324d-05/
      data fpppp( 4,19) /             1.11418071d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.23133659d-06,-4.31272178d-07/
      data fpppp( 5, 3),fpppp( 5, 4)/ 8.21898445d-07, 4.89197340d-07/
      data fpppp( 5, 5),fpppp( 5, 6)/ 1.08405862d-06,-3.68023468d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 6.57641411d-06,-1.54985167d-05/
      data fpppp( 5, 9),fpppp( 5,10)/ 4.45398845d-05,-3.43795683d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 5.25836278d-06, 3.45781978d-06/
      data fpppp( 5,13),fpppp( 5,14)/-1.00933100d-05, 5.53624567d-05/
      data fpppp( 5,15),fpppp( 5,16)/-5.97059986d-05,-2.43588934d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 4.24018781d-05,-5.26008851d-05/
      data fpppp( 5,19) /            -1.01385495d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.76650268d-06, 4.72794259d-07/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.91508616d-06,-4.73364597d-07/
      data fpppp( 6, 5),fpppp( 6, 6)/-2.98181013d-07, 5.07859737d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/-8.95976657d-06, 2.23895360d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-6.84090356d-05, 5.28302523d-05/
      data fpppp( 6,11),fpppp( 6,12)/-5.13415510d-06, 1.83897557d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 9.90422752d-06,-2.35318653d-05/
      data fpppp( 6,15),fpppp( 6,16)/ 1.58090843d-05, 1.46196516d-06/
      data fpppp( 6,17),fpppp( 6,18)/-4.99653820d-05, 1.06097660d-04/
      data fpppp( 6,19) /             2.02877734d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.10572032d-06,-3.33812483d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/ 9.17122872d-07, 4.91461900d-07/
      data fpppp( 7, 5),fpppp( 7, 6)/-1.10814664d-07,-4.43435256d-07/
      data fpppp( 7, 7),fpppp( 7, 8)/ 2.87254128d-07,-9.19675473d-06/
      data fpppp( 7, 9),fpppp( 7,10)/ 4.09561652d-05,-3.25559433d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 4.45235978d-06,-2.02362838d-06/
      data fpppp( 7,13),fpppp( 7,14)/-2.16207702d-06, 7.79281848d-06/
      data fpppp( 7,15),fpppp( 7,16)/-3.91511755d-06,-7.72966571d-06/
      data fpppp( 7,17),fpppp( 7,18)/ 2.69352225d-05,-6.10513489d-05/
      data fpppp( 7,19) /            -1.13994646d-04 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 2.08214512d-07, 2.27838103d-08/
      data fpppp( 8, 3),fpppp( 8, 4)/-1.22553818d-07,-4.08217204d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.68552496d-06,-4.42146332d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 5.42909267d-06,-7.35280533d-07/
      data fpppp( 8, 9),fpppp( 8,10)/-9.46291405d-06, 8.49143940d-06/
      data fpppp( 8,11),fpppp( 8,12)/-2.44366985d-06, 2.69516263d-06/
      data fpppp( 8,13),fpppp( 8,14)/-1.26032282d-07,-1.25458193d-06/
      data fpppp( 8,15),fpppp( 8,16)/ 6.86192006d-07,-2.15535344d-06/
      data fpppp( 8,17),fpppp( 8,18)/-4.45010946d-06, 1.56021926d-05/
      data fpppp( 8,19) /             2.92696210d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 2.66054644d-08, 3.91908463d-08/
      data fpppp( 9, 3),fpppp( 9, 4)/-8.67052191d-08, 2.72083781d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.33419502d-06, 2.70625101d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/-2.18456488d-06, 2.68674717d-07/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.94523969d-06,-1.67560692d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 1.04774124d-06,-1.32891604d-06/
      data fpppp( 9,13),fpppp( 9,14)/ 4.12360189d-07, 2.52786977d-07/
      data fpppp( 9,15),fpppp( 9,16)/-3.32915503d-07,-2.47138158d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 7.29350908d-07,-2.03974614d-06/
      data fpppp( 9,19) /            -3.73867572d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-7.88848206d-09,-1.90347785d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 2.96387373d-08,-1.31057090d-07/
      data fpppp(10, 5),fpppp(10, 6)/ 5.19233808d-07,-8.66751940d-07/
      data fpppp(10, 7),fpppp(10, 8)/ 6.06659318d-07,-3.76974217d-08/
      data fpppp(10, 9),fpppp(10,10)/-5.18518928d-07, 5.21442076d-07/
      data fpppp(10,11),fpppp(10,12)/-3.52496025d-07, 4.39254731d-07/
      data fpppp(10,13),fpppp(10,14)/-1.71308932d-07,-1.06179878d-07/
      data fpppp(10,15),fpppp(10,16)/-5.06653469d-08, 4.74645150d-08/
      data fpppp(10,17),fpppp(10,18)/-9.01754267d-08, 5.98478495d-07/
      data fpppp(10,19) /             1.05104846d-06 /
      data fpppp(11, 1),fpppp(11, 2)/-1.67854428d-08,-3.58391888d-09/
      data fpppp(11, 3),fpppp(11, 4)/ 8.01292333d-09, 2.52261516d-08/
      data fpppp(11, 5),fpppp(11, 6)/-1.26929158d-07, 2.16430952d-07/
      data fpppp(11, 7),fpppp(11, 8)/-1.56580274d-07, 1.24722743d-08/
      data fpppp(11, 9),fpppp(11,10)/ 1.33914676d-07,-1.28833333d-07/
      data fpppp(11,11),fpppp(11,12)/ 1.31851984d-07,-1.83867419d-07/
      data fpppp(11,13),fpppp(11,14)/ 1.02324545d-07,-3.80989510d-08/
      data fpppp(11,15),fpppp(11,16)/ 3.42538363d-08,-3.13961963d-08/
      data fpppp(11,17),fpppp(11,18)/ 5.53790315d-08,-1.61604483d-07/
      data fpppp(11,19) /            -2.91799786d-07 /
      data fpppp(12, 1),fpppp(12, 2)/ 1.57186037d-09, 2.87240171d-10/
      data fpppp(12, 3),fpppp(12, 4)/ 1.00817955d-10,-3.92929721d-09/
      data fpppp(12, 5),fpppp(12, 6)/ 2.70186965d-08,-4.70335832d-08/
      data fpppp(12, 7),fpppp(12, 8)/ 2.93727609d-08,-2.97388616d-09/
      data fpppp(12, 9),fpppp(12,10)/-2.77219161d-08, 2.70020211d-08/
      data fpppp(12,11),fpppp(12,12)/-8.77283411d-09, 3.05478788d-08/
      data fpppp(12,13),fpppp(12,14)/-2.54600513d-08,-1.87403546d-09/
      data fpppp(12,15),fpppp(12,16)/-1.30803223d-08, 9.49128493d-09/
      data fpppp(12,17),fpppp(12,18)/-2.09443398d-09, 3.15833617d-08/
      data fpppp(12,19) /             5.23287246d-08 /
      data fpppp(13, 1),fpppp(13, 2)/ 2.60220480d-08, 1.32404235d-08/
      data fpppp(13, 3),fpppp(13, 4)/-2.18945617d-08, 2.12072158d-08/
      data fpppp(13, 5),fpppp(13, 6)/-2.51354642d-08, 3.20286882d-08/
      data fpppp(13, 7),fpppp(13, 8)/-2.58578509d-08,-4.33907171d-09/
      data fpppp(13, 9),fpppp(13,10)/ 1.53364877d-08,-2.40771143d-08/
      data fpppp(13,11),fpppp(13,12)/-8.78469773d-09, 8.24866234d-08/
      data fpppp(13,13),fpppp(13,14)/-1.00391111d-07, 8.49110010d-08/
      data fpppp(13,15),fpppp(13,16)/-4.82346355d-08, 3.73795608d-08/
      data fpppp(13,17),fpppp(13,18)/-2.56787993d-08,-1.10128191d-08/
      data fpppp(13,19) /            -1.85537932d-08 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine s1_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -7.04510000d-03 /
      data f( 1, 3),f( 1, 4) / -1.11215000d-02 , -1.19545000d-02 /
      data f( 1, 5),f( 1, 6) / -1.05244000d-02 , -7.87560000d-03 /
      data f( 1, 7),f( 1, 8) / -5.34150000d-03 , -3.41930000d-03 /
      data f( 1, 9),f( 1,10) / -1.02060000d-03 ,  2.37430000d-03 /
      data f( 1,11),f( 1,12) /  6.70970000d-03 ,  1.22378000d-02 /
      data f( 1,13),f( 1,14) /  2.02449000d-02 ,  3.29711000d-02 /
      data f( 1,15),f( 1,16) /  5.17336000d-02 ,  6.28771000d-02 /
      data f( 1,17),f( 1,18) /  6.33255000d-02 ,  5.78454000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -5.24480000d-03 /
      data f( 2, 3),f( 2, 4) / -8.36960000d-03 , -8.99710000d-03 /
      data f( 2, 5),f( 2, 6) / -7.83070000d-03 , -5.86600000d-03 /
      data f( 2, 7),f( 2, 8) / -4.09180000d-03 , -2.65160000d-03 /
      data f( 2, 9),f( 2,10) / -7.44800000d-04 ,  2.00060000d-03 /
      data f( 2,11),f( 2,12) /  5.56910000d-03 ,  1.01552000d-02 /
      data f( 2,13),f( 2,14) /  1.67313000d-02 ,  2.67461000d-02 /
      data f( 2,15),f( 2,16) /  4.03712000d-02 ,  5.07003000d-02 /
      data f( 2,17),f( 2,18) /  5.54038000d-02 ,  5.15567000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -3.89500000d-03 /
      data f( 3, 3),f( 3, 4) / -6.27760000d-03 , -6.77480000d-03 /
      data f( 3, 5),f( 3, 6) / -5.88370000d-03 , -4.43590000d-03 /
      data f( 3, 7),f( 3, 8) / -3.14550000d-03 , -2.02920000d-03 /
      data f( 3, 9),f( 3,10) / -5.25800000d-04 ,  1.66180000d-03 /
      data f( 3,11),f( 3,12) /  4.57510000d-03 ,  8.34600000d-03 /
      data f( 3,13),f( 3,14) /  1.36759000d-02 ,  2.15601000d-02 /
      data f( 3,15),f( 3,16) /  3.19002000d-02 ,  4.07846000d-02 /
      data f( 3,17),f( 3,18) /  4.58290000d-02 ,  4.01972000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -2.47450000d-03 /
      data f( 4, 3),f( 4, 4) / -4.03710000d-03 , -4.39410000d-03 /
      data f( 4, 5),f( 4, 6) / -3.83670000d-03 , -2.92310000d-03 /
      data f( 4, 7),f( 4, 8) / -2.08800000d-03 , -1.31650000d-03 /
      data f( 4, 9),f( 4,10) / -2.87500000d-04 ,  1.23500000d-03 /
      data f( 4,11),f( 4,12) /  3.34140000d-03 ,  6.11950000d-03 /
      data f( 4,13),f( 4,14) /  9.95870000d-03 ,  1.54506000d-02 /
      data f( 4,15),f( 4,16) /  2.24788000d-02 ,  2.90737000d-02 /
      data f( 4,17),f( 4,18) /  3.25312000d-02 ,  2.67331000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -1.13040000d-03 /
      data f( 5, 3),f( 5, 4) / -1.86860000d-03 , -2.05950000d-03 /
      data f( 5, 5),f( 5, 6) / -1.81600000d-03 , -1.39130000d-03 /
      data f( 5, 7),f( 5, 8) / -9.79400000d-04 , -5.78500000d-04 /
      data f( 5, 9),f( 5,10) / -6.94000000d-05 ,  7.19900000d-04 /
      data f( 5,11),f( 5,12) /  1.88110000d-03 ,  3.51360000d-03 /
      data f( 5,13),f( 5,14) /  5.69220000d-03 ,  8.70470000d-03 /
      data f( 5,15),f( 5,16) /  1.20044000d-02 ,  1.59199000d-02 /
      data f( 5,17),f( 5,18) /  1.76395000d-02 ,  1.34219000d-02 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -4.99000000d-04 /
      data f( 6, 3),f( 6, 4) / -8.47600000d-04 , -9.55400000d-04 /
      data f( 6, 5),f( 6, 6) / -8.27500000d-04 , -6.11900000d-04 /
      data f( 6, 7),f( 6, 8) / -4.20900000d-04 , -2.08100000d-04 /
      data f( 6, 9),f( 6,10) / -2.70000000d-06 ,  4.05300000d-04 /
      data f( 6,11),f( 6,12) /  9.95200000d-04 ,  1.89420000d-03 /
      data f( 6,13),f( 6,14) /  3.12750000d-03 ,  4.80420000d-03 /
      data f( 6,15),f( 6,16) /  6.87770000d-03 ,  8.80580000d-03 /
      data f( 6,17),f( 6,18) /  9.32110000d-03 ,  6.47860000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -2.14900000d-04 /
      data f( 7, 3),f( 7, 4) / -3.54500000d-04 , -3.85900000d-04 /
      data f( 7, 5),f( 7, 6) / -3.30000000d-04 , -2.37000000d-04 /
      data f( 7, 7),f( 7, 8) / -1.38300000d-04 , -5.97000000d-05 /
      data f( 7, 9),f( 7,10) /  4.32000000d-05 ,  2.18500000d-04 /
      data f( 7,11),f( 7,12) /  5.08000000d-04 ,  9.58200000d-04 /
      data f( 7,13),f( 7,14) /  1.60530000d-03 ,  2.48760000d-03 /
      data f( 7,15),f( 7,16) /  3.58800000d-03 ,  4.53470000d-03 /
      data f( 7,17),f( 7,18) /  4.67060000d-03 ,  3.15850000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -8.99000000d-05 /
      data f( 8, 3),f( 8, 4) / -1.46400000d-04 , -1.54700000d-04 /
      data f( 8, 5),f( 8, 6) / -1.24400000d-04 , -7.73000000d-05 /
      data f( 8, 7),f( 8, 8) / -3.16000000d-05 ,  2.90000000d-06 /
      data f( 8, 9),f( 8,10) /  4.41000000d-05 ,  1.17500000d-04 /
      data f( 8,11),f( 8,12) /  2.46100000d-04 ,  4.58700000d-04 /
      data f( 8,13),f( 8,14) /  7.79000000d-04 ,  1.22300000d-03 /
      data f( 8,15),f( 8,16) /  1.76570000d-03 ,  2.22860000d-03 /
      data f( 8,17),f( 8,18) /  2.24830000d-03 ,  1.48270000d-03 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -3.75000000d-05 /
      data f( 9, 3),f( 9, 4) / -5.89000000d-05 , -5.83000000d-05 /
      data f( 9, 5),f( 9, 6) / -4.05000000d-05 , -1.54000000d-05 /
      data f( 9, 7),f( 9, 8) /  7.50000000d-06 ,  2.34000000d-05 /
      data f( 9, 9),f( 9,10) /  3.74000000d-05 ,  6.47000000d-05 /
      data f( 9,11),f( 9,12) /  1.13900000d-04 ,  2.06000000d-04 /
      data f( 9,13),f( 9,14) /  3.52600000d-04 ,  5.63700000d-04 /
      data f( 9,15),f( 9,16) /  8.25700000d-04 ,  1.04780000d-03 /
      data f( 9,17),f( 9,18) /  1.05010000d-03 ,  6.79700000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -8.80000000d-06 /
      data f(10, 3),f(10, 4) / -1.15000000d-05 , -7.30000000d-06 /
      data f(10, 5),f(10, 6) /  1.30000000d-06 ,  1.03000000d-05 /
      data f(10, 7),f(10, 8) /  1.78000000d-05 ,  2.22000000d-05 /
      data f(10, 9),f(10,10) /  2.30000000d-05 ,  2.26000000d-05 /
      data f(10,11),f(10,12) /  2.33000000d-05 ,  3.20000000d-05 /
      data f(10,13),f(10,14) /  5.34000000d-05 ,  9.20000000d-05 /
      data f(10,15),f(10,16) /  1.46700000d-04 ,  2.00600000d-04 /
      data f(10,17),f(10,18) /  2.11600000d-04 ,  1.40600000d-04 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -3.00000000d-06 /
      data f(11, 3),f(11, 4) / -3.40000000d-06 , -5.00000000d-07 /
      data f(11, 5),f(11, 6) /  3.70000000d-06 ,  7.40000000d-06 /
      data f(11, 7),f(11, 8) /  9.80000000d-06 ,  1.10000000d-05 /
      data f(11, 9),f(11,10) /  1.10000000d-05 ,  9.60000000d-06 /
      data f(11,11),f(11,12) /  7.10000000d-06 ,  4.30000000d-06 /
      data f(11,13),f(11,14) /  3.90000000d-06 ,  6.30000000d-06 /
      data f(11,15),f(11,16) /  1.46000000d-05 ,  2.77000000d-05 /
      data f(11,17),f(11,18) /  3.63000000d-05 ,  2.72000000d-05 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 , -8.00000000d-07 /
      data f(12, 3),f(12, 4) / -5.00000000d-07 ,  8.00000000d-07 /
      data f(12, 5),f(12, 6) /  2.40000000d-06 ,  3.50000000d-06 /
      data f(12, 7),f(12, 8) /  4.00000000d-06 ,  4.00000000d-06 /
      data f(12, 9),f(12,10) /  3.80000000d-06 ,  3.30000000d-06 /
      data f(12,11),f(12,12) /  2.20000000d-06 ,  4.00000000d-07 /
      data f(12,13),f(12,14) / -1.30000000d-06 , -2.90000000d-06 /
      data f(12,15),f(12,16) / -2.90000000d-06 ,  0.00000000d+00 /
      data f(12,17),f(12,18) /  3.90000000d-06 ,  4.50000000d-06 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(13, 3),f(13, 4) /  6.00000000d-07 ,  1.20000000d-06 /
      data f(13, 5),f(13, 6) /  1.70000000d-06 ,  1.70000000d-06 /
      data f(13, 7),f(13, 8) /  1.40000000d-06 ,  9.00000000d-07 /
      data f(13, 9),f(13,10) /  4.00000000d-07 , -1.00000000d-07 /
      data f(13,11),f(13,12) / -8.00000000d-07 , -1.60000000d-06 /
      data f(13,13),f(13,14) / -2.50000000d-06 , -3.20000000d-06 /
      data f(13,15),f(13,16) / -3.60000000d-06 , -3.20000000d-06 /
      data f(13,17),f(13,18) / -2.10000000d-06 , -1.00000000d-06 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 2.74647468d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-1.39690233d-02, 2.86375064d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-2.03227786d-02, 3.61072274d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-2.31461703d-02, 2.15375838d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-2.45318374d-02, 1.35284372d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-1.95207431d-02,-2.52933274d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-9.87489160d-03,-1.02931062d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-4.21121366d-03, 6.98775774d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-1.60957331d-03, 1.09320753d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/ 7.05365611d-04, 9.05594108d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 3.86831343d-03, 9.27416039d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 7.13552355d-03, 2.54094174d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.14473372d-02, 3.78281702d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.91593521d-02, 1.06423902d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 9.61695907d-02,-1.01345778d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/ 6.85703966d-02,-1.58180792d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/-8.73611434d-02, 9.23629439d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-2.67098572d-01,-5.66980984d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-9.66357008d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 1.78468519d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-1.13463105d-02, 2.03602962d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-1.66140857d-02, 2.79119634d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-1.85673023d-02, 1.78298503d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-1.89745395d-02, 8.40263549d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-1.48030139d-02,-3.54239223d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-7.73493108d-03,-5.66306655d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-3.65307268d-03, 6.15465844d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-1.41342481d-03, 9.04043278d-06/
      data fpp( 2,10,1),fpp( 2,10,2)/ 8.41911634d-04, 7.99961042d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 3.66008743d-03, 8.34712552d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 6.81281004d-03, 1.96678875d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.13154685d-02, 3.23813245d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 2.59748672d-02, 5.71288146d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/ 7.33449616d-02,-4.42785829d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 5.70568496d-02,-7.77744829d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/-4.44364989d-02, 1.78405146d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.45150570d-01,-5.06623576d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-8.53922212d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 1.17182359d-05/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-8.22073475d-03, 1.44755283d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-1.22058785d-02, 2.11236510d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.28496206d-02, 1.41538678d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.15750047d-02, 5.55887786d-06/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-8.19220145d-03,-2.98737923d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-4.69538409d-03,-3.05336093d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-2.97149563d-03, 4.75482294d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.25672746d-03, 7.26006915d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.16198785d-03, 7.25690046d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 3.48133686d-03, 7.25432903d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 6.62323630d-03, 1.51817834d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 1.20207889d-02, 2.55585372d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.27911792d-02, 3.58420678d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 4.41605631d-02,-2.15728082d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 4.23672049d-02,-3.68928349d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/ 1.71421390d-02,-6.12558521d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 8.70808511d-02,-3.58655757d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-5.78045122d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 6.39402289d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-5.31334386d-03, 8.69495421d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-8.07101441d-03, 1.35401603d-05/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-8.30639660d-03, 9.48040474d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-7.00029131d-03, 3.40222078d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-4.98065258d-03,-1.71728785d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-3.32209899d-03,-1.24306939d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-2.38629943d-03, 2.87356542d-06/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-8.81958593d-04, 5.19880770d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/ 9.92099406d-04, 5.94120379d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/ 3.10881884d-03, 6.07037713d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/ 5.86733897d-03, 1.00792877d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/ 1.01137247d-02, 1.72784721d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.80128246d-02, 1.99688239d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.29081485d-02,-4.97576773d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/ 3.15814172d-02,-2.60637530d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.34438692d-02,-7.90132204d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 4.48408760d-02,-2.13219366d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-3.24209317d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 2.46610139d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-2.62645881d-03, 3.73179722d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-4.42522678d-03, 6.13870974d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-4.90735855d-03, 4.55136383d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-4.03726497d-03, 1.71983495d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-2.89539087d-03,-5.58703640d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-2.24565277d-03,-2.53020395d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-1.37694443d-03, 9.10785219d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-7.21296025d-04, 3.10187952d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/ 8.37689191d-04, 3.49369670d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 2.26377760d-03, 5.23733367d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/ 3.76897353d-03, 3.83496863d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 6.71560766d-03, 1.21887918d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 1.11626537d-02,-2.55613593d-06/
      data fpp( 5,15,1),fpp( 5,15,2)/ 2.56679870d-02, 1.52677519d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 2.62639420d-02,-2.15668716d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 2.52055350d-02,-6.07542653d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 2.33558862d-02,-9.16480670d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-1.24911466d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 6.24259099d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.28562090d-03, 1.45548180d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.76807848d-03, 2.57781370d-06/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-1.59616920d-03, 2.68126342d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.62344882d-03, 8.39132639d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.49538392d-03,-7.75793973d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-8.97689912d-04, 7.88043254d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-9.28322852d-04,-1.06837904d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.33542694d-04, 3.04147291d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 4.69143832d-04, 1.05848740d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 1.62167075d-03, 3.63857750d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/ 2.73276692d-03, 2.93320260d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/ 3.86704467d-03, 4.68661212d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/ 5.62616042d-03, 4.92434894d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 2.76470361d-03,-5.76007864d-07/
      data fpp( 6,16,1),fpp( 6,16,2)/ 8.31561467d-03,-1.13443175d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.34931906d-02,-3.88147222d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 1.45651793d-02,-3.48647937d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-3.98921032d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 4.42385288d-07/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-5.66257576d-04, 7.16229425d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-1.17205929d-03, 1.21069701d-06/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-1.53836466d-03, 9.32982523d-07/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-1.25293974d-03, 2.95372896d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-8.31073430d-04, 1.11525893d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-7.85187580d-04,-3.99476468d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-2.37764163d-04, 2.80379978d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-3.12074753d-04, 7.35956554d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/ 3.52935483d-04, 1.11979381d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/ 8.18339410d-04, 1.63686822d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.70155878d-03, 1.97473329d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/ 2.83621367d-03, 2.27819860d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/ 4.34630457d-03, 3.02447232d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 7.36119857d-03,-1.29008789d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 8.70559927d-03,-7.08612077d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 8.85130266d-03,-1.90134290d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/ 5.34019659d-03,-1.57401631d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-1.68099184d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 1.94664879d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-2.67748792d-04, 3.16670242d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-3.83684366d-04, 5.42654154d-07/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-3.69572155d-04, 4.04713143d-07/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-3.70392222d-04, 1.54493276d-07/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-3.45122356d-04,-1.46862466d-08/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.83159770d-04,-1.79748290d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-1.79820496d-04, 6.16794049d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 3.47563160d-05, 3.35030670d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.78314235d-04, 5.30197916d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/ 5.12171613d-04, 8.56177667d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/ 9.36997971d-04, 1.08509142d-06/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.48970064d-03, 1.26545666d-06/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.23662129d-03, 1.27508193d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/ 3.00810210d-03,-4.43784385d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/ 4.02198824d-03,-4.28794439d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 4.57839880d-03,-8.99643805d-06/
      data fpp( 8,18,1),fpp( 8,18,2)/ 3.53723432d-03,-6.84430341d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-6.65234829d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 1.05880727d-07/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-1.05147257d-04, 1.53238546d-07/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.87603247d-04, 2.47165089d-07/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-2.18546719d-04, 1.78101098d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.86291372d-04, 7.24305177d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-1.35637146d-04,-2.98231692d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-1.04573341d-04,-8.51378408d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-5.33538525d-05,-4.96254674d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-9.35051144d-06, 1.69639710d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/ 9.06075755d-05, 1.69066626d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/ 2.45774138d-04, 4.68093787d-07/
      data fpp( 9,12,1),fpp( 9,12,2)/ 4.73649340d-04, 5.32558225d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/ 8.02583775d-04, 6.71673312d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/ 1.23441028d-03, 6.50748528d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/ 1.78159304d-03,-2.20667423d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 2.21364777d-03,-2.16207884d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/ 2.21350215d-03,-4.31901724d-06/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.45806614d-03,-2.92385222d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-2.54357389d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 5.41407416d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-7.28383258d-06, 5.87185169d-08/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.09480746d-05, 7.69851909d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.03737653d-05, 4.73407194d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.19297738d-05,-2.34806865d-09/
      data fpp(10, 6,1),fpp(10, 6,2)/-9.12738341d-06,-1.39484448d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-2.10009182d-06,-3.18581521d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/-3.22819440d-06,-4.46189470d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/ 4.67337630d-06,-5.66606011d-09/
      data fpp(10,10,1),fpp(10,10,2)/ 2.00201559d-05,-4.71681260d-09/
      data fpp(10,11,1),fpp(10,11,2)/ 4.93917800d-05, 9.05333105d-08/
      data fpp(10,12,1),fpp(10,12,2)/ 9.89529934d-05, 1.22583571d-07/
      data fpp(10,13,1),fpp(10,13,2)/ 1.68998356d-04, 1.81132407d-07/
      data fpp(10,14,1),fpp(10,14,2)/ 2.59858517d-04, 1.84886801d-07/
      data fpp(10,15,1),fpp(10,15,2)/ 3.57169819d-04, 4.53203875d-08/
      data fpp(10,16,1),fpp(10,16,2)/ 4.34462568d-04,-4.14168351d-07/
      data fpp(10,17,1),fpp(10,17,2)/ 4.17694142d-04,-9.62646982d-07/
      data fpp(10,18,1),fpp(10,18,2)/ 2.58584412d-04,-6.55243719d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-5.92378140d-07/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 2.01188434d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-3.11741248d-06, 2.37623132d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-4.40445404d-06, 4.08319039d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/-5.15821966d-06, 1.09100714d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/-2.38953322d-06,-6.47218942d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/ 5.46679845d-07,-1.50213137d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/ 3.17370837d-06,-1.14425558d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/ 6.26663010d-06,-1.12084633d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/ 5.05700624d-06,-1.57235912d-08/
      data fpp(11,10,1),fpp(11,10,2)/ 3.91180109d-06,-9.89717187d-09/
      data fpp(11,11,1),fpp(11,11,2)/ 3.05874210d-06,-1.06877213d-08/
      data fpp(11,12,1),fpp(11,12,2)/ 8.33868594d-06, 3.46480571d-08/
      data fpp(11,13,1),fpp(11,13,2)/ 1.96228011d-05, 1.60954928d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 4.21556535d-05, 6.89699716d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 7.11276792d-05, 6.20246209d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 9.43019557d-05,-2.90684551d-08/
      data fpp(11,17,1),fpp(11,17,2)/ 9.49212784d-05,-2.15750801d-07/
      data fpp(11,18,1),fpp(11,18,2)/ 6.17962073d-05,-1.69928343d-07/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-1.90535829d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 1.22177356d-08/
      data fpp(12, 2,1),fpp(12, 2,2)/-1.84651750d-06, 1.05645288d-08/
      data fpp(12, 3,1),fpp(12, 3,2)/-2.63410919d-06, 1.15241491d-08/
      data fpp(12, 4,1),fpp(12, 4,2)/-1.99335607d-06, 3.33887467d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/-7.12093357d-07,-6.87964779d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/ 9.40664031d-07,-5.82028350d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/ 2.60525833d-06,-5.83921822d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/ 3.36167398d-06,-8.22843627d-10/
      data fpp(12, 9,1),fpp(12, 9,2)/ 3.89859875d-06,-2.86940727d-09/
      data fpp(12,10,1),fpp(12,10,2)/ 4.53263978d-06,-5.69952727d-09/
      data fpp(12,11,1),fpp(12,11,2)/ 6.17325158d-06,-1.03324836d-08/
      data fpp(12,12,1),fpp(12,12,2)/ 1.04922628d-05, 5.02946179d-09/
      data fpp(12,13,1),fpp(12,13,2)/ 1.83104398d-05,-3.78536352d-09/
      data fpp(12,14,1),fpp(12,14,2)/ 3.05188693d-05, 1.61119923d-08/
      data fpp(12,15,1),fpp(12,15,2)/ 4.59194642d-05, 3.53373943d-08/
      data fpp(12,16,1),fpp(12,16,2)/ 5.95296089d-05, 1.65384304d-08/
      data fpp(12,17,1),fpp(12,17,2)/ 6.00207443d-05,-4.14911159d-08/
      data fpp(12,18,1),fpp(12,18,2)/ 3.84307585d-05,-4.85739669d-08/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-7.02130166d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 2.09431162d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/ 1.99825875d-06, 1.81137676d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/ 3.05455460d-06, 2.66018136d-09/
      data fpp(13, 4,1),fpp(13, 4,2)/ 5.25917803d-06,-4.52102179d-10/
      data fpp(13, 5,1),fpp(13, 5,2)/ 6.18104668d-06,-6.85177264d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/ 5.90466798d-06,-2.14080726d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 4.09737084d-06,-2.58499832d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/ 3.13166301d-06, 4.80800527d-10/
      data fpp(13, 9,1),fpp(13, 9,2)/ 2.27570062d-06, 6.61796210d-10/
      data fpp(13,10,1),fpp(13,10,2)/-1.75381989d-06,-3.12798537d-09/
      data fpp(13,11,1),fpp(13,11,2)/-9.84912579d-06,-1.49854746d-10/
      data fpp(13,12,1),fpp(13,12,2)/-2.69461314d-05,-2.27259565d-09/
      data fpp(13,13,1),fpp(13,13,2)/-5.09427199d-05, 3.24023735d-09/
      data fpp(13,14,1),fpp(13,14,2)/-8.54844346d-05, 1.31164625d-09/
      data fpp(13,15,1),fpp(13,15,2)/-1.21872232d-04, 9.51317766d-09/
      data fpp(13,16,1),fpp(13,16,2)/-1.47439804d-04, 8.63564310d-09/
      data fpp(13,17,1),fpp(13,17,2)/-1.39322872d-04,-2.05575007d-09/
      data fpp(13,18,1),fpp(13,18,2)/-8.63403793d-05,-4.12642838d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-2.29367858d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 1.17085297d-04, 7.59855404d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/ 3.58886271d-05,-7.71823769d-06/
      data fpppp( 1, 5),fpppp( 1, 6)/ 8.12477994d-05, 6.65327280d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/-6.92932845d-05,-2.82900013d-05/
      data fpppp( 1, 9),fpppp( 1,10)/-1.26896593d-06, 1.61637792d-05/
      data fpppp( 1,11),fpppp( 1,12)/-1.25056172d-05, 4.01144278d-05/
      data fpppp( 1,13),fpppp( 1,14)/-8.52758834d-05, 1.10500118d-03/
      data fpppp( 1,15),fpppp( 1,16)/-1.37683544d-03,-1.27422539d-03/
      data fpppp( 1,17),fpppp( 1,18)/-1.22620377d-03, 4.75068712d-03/
      data fpppp( 1,19) /             9.03361536d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 8.86381536d-05, 6.03592809d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 3.46368370d-05,-3.31079496d-08/
      data fpppp( 2, 5),fpppp( 2, 6)/ 5.82543576d-05, 4.17414443d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/-5.14267032d-05,-1.52080950d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 1.72645129d-06, 9.24360422d-06/
      data fpppp( 2,11),fpppp( 2,12)/-4.93050712d-06, 3.05512332d-05/
      data fpppp( 2,13),fpppp( 2,14)/-3.62782759d-05, 7.23966285d-04/
      data fpppp( 2,15),fpppp( 2,16)/-8.96945125d-04,-9.55678166d-04/
      data fpppp( 2,17),fpppp( 2,18)/-3.92656408d-04, 2.57306046d-03/
      data fpppp( 2,19) /             4.85229301d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 5.17539092d-05, 4.14436945d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 3.66067705d-05, 1.26133289d-05/
      data fpppp( 3, 5),fpppp( 3, 6)/ 2.80413868d-05, 1.71236825d-06/
      data fpppp( 3, 7),fpppp( 3, 8)/-2.80500135d-05, 4.11195103d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 1.10549927d-05,-6.09509337d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 7.36340285d-06, 2.59945076d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 2.39977560d-05, 2.00384730d-04/
      data fpppp( 3,15),fpppp( 3,16)/-1.89597054d-04,-8.31761043d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 2.11073877d-03,-1.90136734d-03/
      data fpppp( 3,19) /            -3.92644319d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 2.63089420d-05, 2.47200141d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/ 2.81513998d-05, 1.40116881d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/ 8.29109591d-06,-4.36406459d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.24999464d-05, 1.09986084d-05/
      data fpppp( 4, 9),fpppp( 4,10)/ 2.61798999d-06, 7.12461121d-07/
      data fpppp( 4,11),fpppp( 4,12)/ 9.09185145d-06, 1.42817498d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 7.44673842d-05,-8.01348598d-05/
      data fpppp( 4,15),fpppp( 4,16)/ 6.58454919d-05, 4.34295828d-05/
      data fpppp( 4,17),fpppp( 4,18)/-4.82128233d-05,-4.78505007d-04/
      data fpppp( 4,19) /            -8.12040112d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 3.45744241d-06, 8.13693344d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/ 1.36562744d-05, 1.62361407d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/ 2.53268462d-06,-1.00600488d-05/
      data fpppp( 5, 7),fpppp( 5, 8)/ 8.17935129d-06,-9.51914157d-06/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.71136184d-05,-4.73512326d-06/
      data fpppp( 5,11),fpppp( 5,12)/-6.14693352d-06, 3.40693079d-05/
      data fpppp( 5,13),fpppp( 5,14)/-4.36440055d-05, 2.30531431d-04/
      data fpppp( 5,15),fpppp( 5,16)/-2.74984489d-04, 3.48438359d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 3.63474233d-05,-2.27708042d-04/
      data fpppp( 5,19) /            -4.15889493d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 9.92718943d-06, 7.21644993d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 9.39681022d-06,-5.54167874d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 8.18569844d-07, 1.15880710d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.89931072d-05, 2.66847406d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-2.21959460d-05, 1.85231787d-05/
      data fpppp( 6,11),fpppp( 6,12)/-2.88122228d-06,-9.48413385d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 4.22086515d-05,-1.21860191d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 1.67997761d-04,-4.53887786d-05/
      data fpppp( 6,17),fpppp( 6,18)/-8.84275592d-06,-1.65575428d-04/
      data fpppp( 6,19) /            -2.67085615d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/-3.32695663d-06,-1.13220226d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/ 1.40718943d-06, 8.85424269d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 2.27965771d-06,-9.78639046d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.43070766d-05,-1.73496621d-05/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.77875316d-05,-9.44121452d-06/
      data fpppp( 7,11),fpppp( 7,12)/ 8.00094797d-06, 2.50634900d-06/
      data fpppp( 7,13),fpppp( 7,14)/-2.94021207d-06, 3.17806593d-05/
      data fpppp( 7,15),fpppp( 7,16)/-3.38942392d-05, 3.56669952d-06/
      data fpppp( 7,17),fpppp( 7,18)/-5.22943976d-05,-1.37976758d-05/
      data fpppp( 7,19) /            -2.26033089d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.80031876d-06, 1.38906761d-06/
      data fpppp( 8, 3),fpppp( 8, 4)/ 1.75220382d-06,-5.95015716d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/-2.68077723d-07, 3.23272267d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/-4.46124979d-06, 5.09487774d-06/
      data fpppp( 8, 9),fpppp( 8,10)/-3.24400887d-06, 3.62002418d-06/
      data fpppp( 8,11),fpppp( 8,12)/ 1.81879656d-07, 1.11059600d-06/
      data fpppp( 8,13),fpppp( 8,14)/ 3.04831487d-06,-1.65077657d-06/
      data fpppp( 8,15),fpppp( 8,16)/ 5.02840100d-06,-3.91850734d-06/
      data fpppp( 8,17),fpppp( 8,18)/-1.68029068d-05,-2.47243679d-05/
      data fpppp( 8,19) /            -3.40638115d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-5.98482503d-08, 2.24009480d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/ 5.25286349d-07, 7.65596234d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/ 2.04257867d-07,-4.78695024d-07/
      data fpppp( 9, 7),fpppp( 9, 8)/ 5.35097008d-07,-4.52352001d-07/
      data fpppp( 9, 9),fpppp( 9,10)/ 8.41342146d-07, 4.44268168d-07/
      data fpppp( 9,11),fpppp( 9,12)/ 6.94093705d-07, 1.14187542d-06/
      data fpppp( 9,13),fpppp( 9,14)/ 8.01958523d-07, 1.82381470d-06/
      data fpppp( 9,15),fpppp( 9,16)/-1.17584176d-06,-4.02812993d-06/
      data fpppp( 9,17),fpppp( 9,18)/-8.64365916d-06,-6.71465698d-06/
      data fpppp( 9,19) /            -6.65552092d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 3.43552854d-08, 2.74979276d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 7.28284355d-08,-6.44985862d-08/
      data fpppp(10, 5),fpppp(10, 6)/ 5.73468409d-08, 9.66151534d-08/
      data fpppp(10, 7),fpppp(10, 8)/-1.90313382d-07, 1.75314723d-07/
      data fpppp(10, 9),fpppp(10,10)/ 3.08348869d-08, 1.48058261d-07/
      data fpppp(10,11),fpppp(10,12)/ 2.18422747d-07, 1.89626105d-07/
      data fpppp(10,13),fpppp(10,14)/ 2.52121782d-07, 5.07746523d-08/
      data fpppp(10,15),fpppp(10,16)/-6.81518800d-08,-9.79280324d-07/
      data fpppp(10,17),fpppp(10,18)/-1.65839736d-06,-9.27608447d-07/
      data fpppp(10,19) /            -5.99649806d-07 /
      data fpppp(11, 1),fpppp(11, 2)/ 2.87356324d-08, 2.33817620d-08/
      data fpppp(11, 3),fpppp(11, 4)/-1.24404254d-08, 5.83764965d-08/
      data fpppp(11, 5),fpppp(11, 6)/-9.71843752d-09,-9.45114932d-09/
      data fpppp(11, 7),fpppp(11, 8)/ 2.89719630d-08,-7.84831104d-08/
      data fpppp(11, 9),fpppp(11,10)/ 2.68077428d-08,-2.48827375d-08/
      data fpppp(11,11),fpppp(11,12)/ 9.02519764d-08, 3.18550020d-08/
      data fpppp(11,13),fpppp(11,14)/ 1.42578292d-07, 7.27560698d-08/
      data fpppp(11,15),fpppp(11,16)/-4.72521794d-08,-2.31612302d-07/
      data fpppp(11,17),fpppp(11,18)/-3.79595842d-07,-2.74667961d-07/
      data fpppp(11,19) /            -2.42000486d-07 /
      data fpppp(12, 1),fpppp(12, 2)/ 7.42692709d-09, 9.52554047d-09/
      data fpppp(12, 3),fpppp(12, 4)/ 1.80064600d-08, 4.14930794d-09/
      data fpppp(12, 5),fpppp(12, 6)/ 3.82688356d-09, 2.83283841d-09/
      data fpppp(12, 7),fpppp(12, 8)/-1.44480228d-08, 4.68534477d-10/
      data fpppp(12, 9),fpppp(12,10)/-5.95567891d-10, 7.74071242d-09/
      data fpppp(12,11),fpppp(12,12)/ 3.00269644d-08, 3.28553961d-08/
      data fpppp(12,13),fpppp(12,14)/ 4.85013960d-08, 3.65541719d-08/
      data fpppp(12,15),fpppp(12,16)/-3.18816198d-09,-1.31228534d-07/
      data fpppp(12,17),fpppp(12,18)/-2.59038256d-07,-1.57485715d-07/
      data fpppp(12,19) /            -1.21465249d-07 /
      data fpppp(13, 1),fpppp(13, 2)/-2.81304100d-08,-1.38038774d-08/
      data fpppp(13, 3),fpppp(13, 4)/ 2.68281451d-08,-2.46090472d-08/
      data fpppp(13, 5),fpppp(13, 6)/-5.35724407d-09,-2.58568168d-08/
      data fpppp(13, 7),fpppp(13, 8)/ 1.69294042d-08, 8.63455913d-09/
      data fpppp(13, 9),fpppp(13,10)/-4.48829143d-08,-1.95163894d-08/
      data fpppp(13,11),fpppp(13,12)/-1.20998651d-07,-3.65909894d-08/
      data fpppp(13,13),fpppp(13,14)/-1.46612364d-07,-9.66713187d-09/
      data fpppp(13,15),fpppp(13,16)/ 7.45159303d-08, 3.60816916d-07/
      data fpppp(13,17),fpppp(13,18)/ 5.03286684d-07, 3.17969984d-07/
      data fpppp(13,19) /             2.26306561d-07 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine s2_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -2.06700000d-03 /
      data f( 1, 3),f( 1, 4) / -6.56750000d-03 , -1.03002000d-02 /
      data f( 1, 5),f( 1, 6) / -1.10815000d-02 , -8.60080000d-03 /
      data f( 1, 7),f( 1, 8) / -4.85100000d-03 , -2.19750000d-03 /
      data f( 1, 9),f( 1,10) / -1.75460000d-03 , -4.64450000d-03 /
      data f( 1,11),f( 1,12) / -1.10298000d-02 , -2.04443000d-02 /
      data f( 1,13),f( 1,14) / -3.30256000d-02 , -4.94916000d-02 /
      data f( 1,15),f( 1,16) / -6.24118000d-02 , -5.65858000d-02 /
      data f( 1,17),f( 1,18) / -3.18334000d-02 , -7.15840000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.57850000d-03 /
      data f( 2, 3),f( 2, 4) / -5.01290000d-03 , -7.78020000d-03 /
      data f( 2, 5),f( 2, 6) / -8.20280000d-03 , -6.29410000d-03 /
      data f( 2, 7),f( 2, 8) / -3.56460000d-03 , -1.57530000d-03 /
      data f( 2, 9),f( 2,10) / -1.26340000d-03 , -3.67300000d-03 /
      data f( 2,11),f( 2,12) / -9.13270000d-03 , -1.72610000d-02 /
      data f( 2,13),f( 2,14) / -2.79270000d-02 , -4.16124000d-02 /
      data f( 2,15),f( 2,16) / -5.18015000d-02 , -4.76043000d-02 /
      data f( 2,17),f( 2,18) / -2.99477000d-02 , -6.67890000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.18520000d-03 /
      data f( 3, 3),f( 3, 4) / -3.77890000d-03 , -5.84950000d-03 /
      data f( 3, 5),f( 3, 6) / -6.11110000d-03 , -4.66410000d-03 /
      data f( 3, 7),f( 3, 8) / -2.62950000d-03 , -1.10800000d-03 /
      data f( 3, 9),f( 3,10) / -8.99800000d-04 , -2.90440000d-03 /
      data f( 3,11),f( 3,12) / -7.50860000d-03 , -1.44527000d-02 /
      data f( 3,13),f( 3,14) / -2.33979000d-02 , -3.43769000d-02 /
      data f( 3,15),f( 3,16) / -4.28937000d-02 , -4.06029000d-02 /
      data f( 3,17),f( 3,18) / -2.67788000d-02 , -7.30480000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -7.56900000d-04 /
      data f( 4, 3),f( 4, 4) / -2.42590000d-03 , -3.76140000d-03 /
      data f( 4, 5),f( 4, 6) / -3.91650000d-03 , -2.97450000d-03 /
      data f( 4, 7),f( 4, 8) / -1.63840000d-03 , -6.15700000d-04 /
      data f( 4, 9),f( 4,10) / -5.19100000d-04 , -2.01590000d-03 /
      data f( 4,11),f( 4,12) / -5.50590000d-03 , -1.08938000d-02 /
      data f( 4,13),f( 4,14) / -1.77288000d-02 , -2.54767000d-02 /
      data f( 4,15),f( 4,16) / -3.16250000d-02 , -3.06970000d-02 /
      data f( 4,17),f( 4,18) / -2.07634000d-02 , -6.43650000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -3.41900000d-04 /
      data f( 5, 3),f( 5, 4) / -1.10540000d-03 , -1.72460000d-03 /
      data f( 5, 5),f( 5, 6) / -1.79550000d-03 , -1.34240000d-03 /
      data f( 5, 7),f( 5, 8) / -6.80200000d-04 , -1.64400000d-04 /
      data f( 5, 9),f( 5,10) / -1.67200000d-04 , -1.04180000d-03 /
      data f( 5,11),f( 5,12) / -3.12490000d-03 , -6.46370000d-03 /
      data f( 5,13),f( 5,14) / -1.07778000d-02 , -1.53221000d-02 /
      data f( 5,15),f( 5,16) / -1.81674000d-02 , -1.79242000d-02 /
      data f( 5,17),f( 5,18) / -1.23077000d-02 , -3.97990000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -1.35300000d-04 /
      data f( 6, 3),f( 6, 4) / -4.78100000d-04 , -7.65600000d-04 /
      data f( 6, 5),f( 6, 6) / -7.84900000d-04 , -5.60300000d-04 /
      data f( 6, 7),f( 6, 8) / -2.42200000d-04 ,  3.32000000d-05 /
      data f( 6, 9),f( 6,10) / -5.80000000d-06 , -4.94200000d-04 /
      data f( 6,11),f( 6,12) / -1.59710000d-03 , -3.48980000d-03 /
      data f( 6,13),f( 6,14) / -6.05290000d-03 , -8.79310000d-03 /
      data f( 6,15),f( 6,16) / -1.06959000d-02 , -1.03351000d-02 /
      data f( 6,17),f( 6,18) / -6.93780000d-03 , -2.28600000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -6.06000000d-05 /
      data f( 7, 3),f( 7, 4) / -1.95800000d-04 , -3.04800000d-04 /
      data f( 7, 5),f( 7, 6) / -3.09800000d-04 , -2.05200000d-04 /
      data f( 7, 7),f( 7, 8) / -4.15000000d-05 ,  6.85000000d-05 /
      data f( 7, 9),f( 7,10) /  3.78000000d-05 , -2.32800000d-04 /
      data f( 7,11),f( 7,12) / -8.46900000d-04 , -1.88080000d-03 /
      data f( 7,13),f( 7,14) / -3.29330000d-03 , -4.80730000d-03 /
      data f( 7,15),f( 7,16) / -5.78590000d-03 , -5.45870000d-03 /
      data f( 7,17),f( 7,18) / -3.51880000d-03 , -1.10700000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -2.35000000d-05 /
      data f( 8, 3),f( 8, 4) / -7.62000000d-05 , -1.17300000d-04 /
      data f( 8, 5),f( 8, 6) / -1.13000000d-04 , -5.83000000d-05 /
      data f( 8, 7),f( 8, 8) /  1.64000000d-05 ,  6.66000000d-05 /
      data f( 8, 9),f( 8,10) /  4.52000000d-05 , -9.64000000d-05 /
      data f( 8,11),f( 8,12) / -4.08500000d-04 , -9.44100000d-04 /
      data f( 8,13),f( 8,14) / -1.67420000d-03 , -2.46260000d-03 /
      data f( 8,15),f( 8,16) / -2.96080000d-03 , -2.74870000d-03 /
      data f( 8,17),f( 8,18) / -1.75060000d-03 , -5.20900000d-04 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -8.20000000d-06 /
      data f( 9, 3),f( 9, 4) / -2.68000000d-05 , -3.94000000d-05 /
      data f( 9, 5),f( 9, 6) / -3.07000000d-05 , -1.00000000d-07 /
      data f( 9, 7),f( 9, 8) /  3.41000000d-05 ,  5.39000000d-05 /
      data f( 9, 9),f( 9,10) /  3.86000000d-05 , -3.64000000d-05 /
      data f( 9,11),f( 9,12) / -1.88900000d-04 , -4.50900000d-04 /
      data f( 9,13),f( 9,14) / -8.17300000d-04 , -1.20970000d-03 /
      data f( 9,15),f( 9,16) / -1.44920000d-03 , -1.33390000d-03 /
      data f( 9,17),f( 9,18) / -8.41100000d-04 , -2.59000000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -5.00000000d-07 /
      data f(10, 3),f(10, 4) / -1.10000000d-06 ,  1.00000000d-06 /
      data f(10, 5),f(10, 6) /  6.90000000d-06 ,  1.61000000d-05 /
      data f(10, 7),f(10, 8) /  2.37000000d-05 ,  2.54000000d-05 /
      data f(10, 9),f(10,10) /  1.83000000d-05 , -1.50000000d-06 /
      data f(10,11),f(10,12) / -3.74000000d-05 , -9.63000000d-05 /
      data f(10,13),f(10,14) / -1.77500000d-04 , -2.62500000d-04 /
      data f(10,15),f(10,16) / -3.11600000d-04 , -2.84300000d-04 /
      data f(10,17),f(10,18) / -1.78600000d-04 , -5.48000000d-05 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(11, 3),f(11, 4) /  9.00000000d-07 ,  2.90000000d-06 /
      data f(11, 5),f(11, 6) /  6.00000000d-06 ,  8.60000000d-06 /
      data f(11, 7),f(11, 8) /  9.60000000d-06 ,  8.70000000d-06 /
      data f(11, 9),f(11,10) /  5.80000000d-06 ,  8.00000000d-07 /
      data f(11,11),f(11,12) / -7.30000000d-06 , -1.97000000d-05 /
      data f(11,13),f(11,14) / -3.71000000d-05 , -5.19000000d-05 /
      data f(11,15),f(11,16) / -5.91000000d-05 , -5.26000000d-05 /
      data f(11,17),f(11,18) / -3.31000000d-05 , -1.02000000d-05 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  3.00000000d-07 /
      data f(12, 3),f(12, 4) /  1.00000000d-06 ,  2.40000000d-06 /
      data f(12, 5),f(12, 6) /  3.60000000d-06 ,  4.00000000d-06 /
      data f(12, 7),f(12, 8) /  3.40000000d-06 ,  2.30000000d-06 /
      data f(12, 9),f(12,10) /  9.00000000d-07 , -5.00000000d-07 /
      data f(12,11),f(12,12) / -2.20000000d-06 , -4.40000000d-06 /
      data f(12,13),f(12,14) / -7.80000000d-06 , -1.05000000d-05 /
      data f(12,15),f(12,16) / -1.15000000d-05 , -9.80000000d-06 /
      data f(12,17),f(12,18) / -6.00000000d-06 , -2.00000000d-06 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(13, 3),f(13, 4) /  8.00000000d-07 ,  1.60000000d-06 /
      data f(13, 5),f(13, 6) /  2.20000000d-06 ,  2.20000000d-06 /
      data f(13, 7),f(13, 8) /  1.60000000d-06 ,  8.00000000d-07 /
      data f(13, 9),f(13,10) /  1.00000000d-07 , -3.00000000d-07 /
      data f(13,11),f(13,12) / -4.00000000d-07 , -3.00000000d-07 /
      data f(13,13),f(13,14) / -1.00000000d-07 ,  1.00000000d-07 /
      data f(13,15),f(13,16) /  1.00000000d-07 ,  1.00000000d-07 /
      data f(13,17),f(13,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-5.60713925d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-2.58641785d-03,-2.48882150d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-9.20129558d-03, 9.61425243d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-1.80713532d-02, 3.24992053d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-2.56080105d-02, 3.74729265d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-2.26130465d-02, 1.33290888d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-1.15235665d-02,-1.46432817d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-4.78685812d-03,-2.05339620d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-4.02153014d-03,-3.58568704d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/-6.35780214d-03,-3.60065564d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-7.69312376d-03,-2.98409041d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/-9.93927903d-03,-2.63818271d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/-1.37028799d-02,-5.46397873d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-8.09174803d-03, 1.18589763d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/-5.44257045d-02, 2.19951882d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/-8.23374806d-02, 2.33105496d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 4.56538363d-02,-1.67898642d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-6.82353930d-02,-1.70590039d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-3.51845981d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-4.35837102d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-2.35230716d-03,-1.89695796d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-7.99955169d-03, 8.10802870d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-1.48430793d-02, 2.65634648d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-1.99844076d-02, 2.63201120d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-1.72615498d-02, 8.03408715d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-8.94922420d-03,-9.20846062d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-3.91321232d-03,-1.56122447d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-3.23093971d-03,-2.89865607d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-5.15418143d-03,-3.17315126d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-6.88132390d-03,-2.70933888d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/-9.37958480d-03,-2.00109323d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/-1.40098830d-02,-4.51248821d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/-1.50400039d-02, 1.93464606d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/-4.39603767d-02, 1.77517040d-04/
      data fpp( 2,16,1),fpp( 2,16,2)/-5.50665388d-02, 1.33763381d-04/
      data fpp( 2,17,1),fpp( 2,17,2)/ 3.10797560d-02, 9.49934358d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-3.51147855d-02,-1.77005125d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-3.82366938d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-3.32575408d-05/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-2.28435353d-03,-1.43719184d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-6.89049767d-03, 6.23521431d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.09513297d-02, 2.08170611d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.25043590d-02, 1.90365411d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-9.84575434d-03, 5.55277429d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-5.37453673d-03,-5.99163832d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-2.79529259d-03,-1.23722210d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-2.19471100d-03,-2.33174776d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-3.46047215d-03,-2.71258687d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-5.73158063d-03,-2.41550477d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-8.79238178d-03,-1.66479406d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/-1.56825880d-02,-2.93191901d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/-2.83032362d-02, 1.18967008d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/-2.51077887d-02, 1.29464387d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 5.58863586d-03, 1.18701752d-04/
      data fpp( 3,17,1),fpp( 3,17,2)/ 2.25071396d-02, 8.77266052d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 4.28845350d-02,-1.30614173d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-2.95421914d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-2.14969005d-05/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.59395014d-03,-9.28119891d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-4.89863999d-03, 3.89569617d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-7.46351492d-03, 1.37084142d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-7.85919834d-03, 1.20946469d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-6.03311901d-03, 3.73899809d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-3.55539475d-03,-3.40463927d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-1.98354982d-03,-8.92444100d-06/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-1.51033686d-03,-1.64635967d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/-2.65563854d-03,-2.08251721d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-5.20384863d-03,-1.98277150d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/-8.00900420d-03,-1.37379681d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.33547847d-02,-1.20464128d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/-2.58325433d-02, 7.14961911d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/-2.65337866d-02, 7.94239363d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/-2.16644270d-02, 9.97326356d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/-1.16069693d-02, 6.19815212d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 9.38073780d-04,-8.40607203d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-1.99162640d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-9.84708332d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-7.00747429d-04,-4.27183335d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-2.61805344d-03, 1.63841673d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-4.18675445d-03, 6.37616642d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-4.22794992d-03, 5.75491758d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-3.20016656d-03, 2.04416327d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-2.04521474d-03,-1.38557065d-06/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-8.36265011d-04,-5.28588068d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-6.32495454d-04,-8.58690662d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/-1.58727339d-03,-1.26744928d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/-2.87273599d-03,-1.32251221d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-5.12935748d-03,-9.76701884d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/-7.79513611d-03,-6.22480257d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/-1.26515196d-02, 2.08542291d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/-2.77928098d-02, 2.47478861d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/-2.37158151d-02, 6.54642264d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/-1.40411821d-02, 3.57932083d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/-4.50615708d-03,-4.59590595d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-1.12830970d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-4.70442566d-06/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-6.04660140d-04,-2.07214867d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.26594627d-03, 5.43020350d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-1.65666726d-03, 3.21806727d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.87860196d-03, 2.67671056d-06/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.56621474d-03, 7.09090483d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-7.48546267d-04, 9.69275085d-08/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-7.60190132d-04,-3.65880052d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-5.31681327d-04,-4.32572544d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.23126792d-03,-6.00229771d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-3.78200739d-03,-8.53508370d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/-6.42236587d-03,-7.24536749d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-8.89107083d-03,-2.70744636d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/-1.05757781d-02, 7.44915291d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/-5.96137424d-03, 2.31548347d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/-7.88111250d-03, 3.57475082d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/-6.28750230d-03, 1.60451324d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.21824547d-03,-2.46580378d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-5.93609811d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-1.75137369d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-4.62120090d-05,-7.51252612d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-5.98161493d-04, 2.80384143d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-1.14337650d-03, 1.20171604d-06/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-1.10964224d-03, 1.15275169d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-7.82974480d-04, 7.63277196d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-6.55800186d-04,-6.59860476d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.81744607d-05,-1.34583529d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-6.79792366d-05,-2.39879836d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/-3.56454953d-04,-3.45297128d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-6.61634458d-04,-4.39931651d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.93877906d-03,-4.13776267d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-3.80778058d-03,-1.76563281d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/-6.08216798d-03, 5.11029391d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/-9.83769323d-03, 1.34484572d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/-9.86453488d-03, 1.94438775d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/-7.63040871d-03, 5.53803299d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/-2.97846103d-03,-1.32820094d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-3.06979953d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-6.99351962d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.12891823d-04,-2.93296076d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-2.46207759d-04, 1.20536268d-07/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-3.29026752d-04, 5.07151006d-07/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-3.62029075d-04, 5.74859710d-07/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-2.98687341d-04, 2.17410154d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-5.54529898d-05,-2.44500326d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-5.99120252d-05,-7.09408848d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-6.52017260d-05,-1.21386428d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/-3.42912274d-04,-1.64713403d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-1.05465478d-03,-2.42759959d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-1.95771791d-03,-2.05246760d-06/
      data fpp( 8,13,1),fpp( 8,13,2)/-3.24980684d-03,-1.03253000d-06/
      data fpp( 8,14,1),fpp( 8,14,2)/-4.48194998d-03, 2.68458761d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/-4.72545286d-03, 7.70617956d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/-4.65434799d-03, 9.10869414d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/-2.81006288d-03, 3.01904389d-06/
      data fpp( 8,18,1),fpp( 8,18,2)/-1.09751041d-03,-7.28886968d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.63915652d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-2.68610232d-07/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-2.54206973d-05,-1.02779535d-07/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.01807472d-04, 5.57283741d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-1.70916496d-04, 2.39866039d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.90241458d-04, 2.62807470d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-1.51076157d-04, 2.29040822d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-8.71878550d-05,-1.38423799d-07/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-1.37743870d-06,-3.33208888d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-7.21385942d-06,-6.34740649d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.05495950d-04,-7.09828516d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/-3.70946419d-04,-1.17594529d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-8.74349300d-04,-1.15639033d-06/
      data fpp( 9,13,1),fpp( 9,13,2)/-1.48579205d-03,-4.62493388d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/-2.19323209d-03, 1.44636388d-06/
      data fpp( 9,15,1),fpp( 9,15,2)/-2.78449535d-03, 3.85103785d-06/
      data fpp( 9,16,1),fpp( 9,16,2)/-2.60287318d-03, 4.43748470d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.73813978d-03, 1.04902335d-06/
      data fpp( 9,18,1),fpp( 9,18,2)/-4.12297321d-04,-3.27557810d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-7.33271095d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00,-2.83411274d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-4.69199651d-06,-2.31774527d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.00737047d-05, 3.16121084d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.51371348d-05, 3.78693115d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.02610887d-05, 4.49106456d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/ 1.37214053d-06,-1.95118939d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.44900598d-05,-6.28630702d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.54883287d-05,-8.30358254d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/ 1.16424413d-05,-1.32993628d-07/
      data fpp(10,10,1),fpp(10,10,2)/-2.26560134d-05,-1.46989662d-07/
      data fpp(10,11,1),fpp(10,11,2)/-8.60333527d-05,-2.45047725d-07/
      data fpp(10,12,1),fpp(10,12,2)/-1.88893145d-04,-2.52819438d-07/
      data fpp(10,13,1),fpp(10,13,2)/-3.61720438d-04,-8.16745218d-08/
      data fpp(10,14,1),fpp(10,14,2)/-5.30928736d-04, 3.51517525d-07/
      data fpp(10,15,1),fpp(10,15,2)/-5.97387511d-04, 8.29604420d-07/
      data fpp(10,16,1),fpp(10,16,2)/-5.44206474d-04, 9.14064795d-07/
      data fpp(10,17,1),fpp(10,17,2)/-3.19549213d-04, 2.18136401d-07/
      data fpp(10,18,1),fpp(10,18,2)/-1.31952829d-04,-7.00610400d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-1.55569480d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 6.89641486d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 9.88683292d-07, 9.20717029d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-9.77092958d-08, 1.02749040d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/ 4.65035461d-07, 1.56932137d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/ 2.85812824d-07,-7.04775895d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/ 3.38759460d-06,-1.75021779d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/ 7.02761585d-06,-1.89435294d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.02241240d-05,-2.07237047d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/ 7.44409441d-06,-1.81616520d-08/
      data fpp(11,10,1),fpp(11,10,2)/ 5.20003536d-07,-3.26296872d-08/
      data fpp(11,11,1),fpp(11,11,2)/-1.33201704d-05,-3.73195992d-08/
      data fpp(11,12,1),fpp(11,12,2)/-3.80781197d-05,-7.60919159d-08/
      data fpp(11,13,1),fpp(11,13,2)/-6.37262005d-05, 4.16872627d-08/
      data fpp(11,14,1),fpp(11,14,2)/-1.02652964d-04, 6.53428652d-08/
      data fpp(11,15,1),fpp(11,15,2)/-1.36554602d-04, 1.52941276d-07/
      data fpp(11,16,1),fpp(11,16,2)/-1.27700928d-04, 1.44892029d-07/
      data fpp(11,17,1),fpp(11,17,2)/-8.56633650d-05, 4.74906076d-08/
      data fpp(11,18,1),fpp(11,18,2)/-1.74913607d-05,-1.30854459d-07/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-2.86072770d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 1.52522969d-09/
      data fpp(12, 2,1),fpp(12, 2,2)/-4.62736658d-07, 2.94954062d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/-9.35458141d-07, 1.06766078d-08/
      data fpp(12, 4,1),fpp(12, 4,2)/-1.12300709d-06,-3.65597199d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/ 1.17837435d-07,-8.05271990d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/ 2.47748108d-06,-1.21331484d-08/
      data fpp(12, 7,1),fpp(12, 7,2)/ 4.79947683d-06,-3.41468641d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/ 5.41517519d-06,-4.20810593d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 4.18118112d-06, 2.24711013d-09/
      data fpp(12,10,1),fpp(12,10,2)/-1.02400071d-06,-4.78033461d-09/
      data fpp(12,11,1),fpp(12,11,2)/-1.06859659d-05,-1.12577170d-09/
      data fpp(12,12,1),fpp(12,12,2)/-2.65943761d-05,-2.07165786d-08/
      data fpp(12,13,1),fpp(12,13,2)/-4.99747599d-05, 1.19920860d-08/
      data fpp(12,14,1),fpp(12,14,2)/-7.36594072d-05, 1.47482345d-08/
      data fpp(12,15,1),fpp(12,15,2)/-8.57940795d-05, 3.10149758d-08/
      data fpp(12,16,1),fpp(12,16,2)/-7.83898144d-05, 2.31918621d-08/
      data fpp(12,17,1),fpp(12,17,2)/-4.81973270d-05, 2.21757558d-09/
      data fpp(12,18,1),fpp(12,18,2)/-1.64817279d-05,-2.00621645d-08/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-4.19689178d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 6.07437229d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/-1.56131671d-07, 3.85125542d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/ 2.25522907d-06, 2.52060602d-09/
      data fpp(13, 4,1),fpp(13, 4,2)/ 3.43650355d-06,-1.93367952d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/ 4.60358128d-06,-6.78588795d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.97375946d-06,-6.92276867d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/-2.01223842d-06,-1.52303738d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/-4.40758760d-06, 1.01491819d-09/
      data fpp(13, 9,1),fpp(13, 9,2)/-2.76559056d-06, 3.46336463d-09/
      data fpp(13,10,1),fpp(13,10,2)/ 7.01200035d-06, 3.13162330d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 2.61179830d-05, 2.01014217d-09/
      data fpp(13,12,1),fpp(13,12,2)/ 5.90721880d-05, 8.27808001d-10/
      data fpp(13,13,1),fpp(13,13,2)/ 1.05437380d-04, 6.78625821d-10/
      data fpp(13,14,1),fpp(13,14,2)/ 1.64004704d-04,-3.54231129d-09/
      data fpp(13,15,1),fpp(13,15,2)/ 2.00259540d-04, 1.49061932d-09/
      data fpp(13,16,1),fpp(13,16,2)/ 1.85469907d-04,-2.42016600d-09/
      data fpp(13,17,1),fpp(13,17,2)/ 1.15123664d-04, 2.19004469d-09/
      data fpp(13,18,1),fpp(13,18,2)/ 3.65908639d-05,-3.40012770d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-8.29993615d-10/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-5.79012733d-05,-4.05168496d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/-2.17389212d-05,-7.83825961d-06/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.33095983d-04, 1.07351599d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/-7.68314104d-05,-6.11922607d-05/
      data fpppp( 1, 9),fpppp( 1,10)/-3.66823687d-05, 2.18257365d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 9.43644578d-06,-1.14221538d-04/
      data fpppp( 1,13),fpppp( 1,14)/ 3.56402971d-04,-7.48906379d-04/
      data fpppp( 1,15),fpppp( 1,16)/-4.77482758d-04, 3.76416823d-03/
      data fpppp( 1,17),fpppp( 1,18)/-5.22500460d-03, 2.62301740d-03/
      data fpppp( 1,19) /             5.66041231d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-5.37217361d-05,-3.33777353d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/-1.04635651d-05, 3.45501228d-06/
      data fpppp( 2, 5),fpppp( 2, 6)/ 9.87754697d-05, 7.32942823d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/-5.65845357d-05,-4.35349618d-05/
      data fpppp( 2, 9),fpppp( 2,10)/-3.04999733d-05, 9.20399590d-06/
      data fpppp( 2,11),fpppp( 2,12)/ 5.44994381d-06,-7.72708763d-05/
      data fpppp( 2,13),fpppp( 2,14)/ 1.75711322d-04,-4.09563771d-04/
      data fpppp( 2,15),fpppp( 2,16)/-2.10871350d-04, 2.32190181d-03/
      data fpppp( 2,17),fpppp( 2,18)/-3.24158847d-03, 1.50400187d-03/
      data fpppp( 2,19) /             3.30414060d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-5.14387250d-05,-2.41183233d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 8.60458161d-06, 2.24187251d-05/
      data fpppp( 3, 5),fpppp( 3, 6)/ 5.21886799d-05, 2.15245910d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-2.95302651d-05,-1.69219383d-05/
      data fpppp( 3, 9),fpppp( 3,10)/-2.15017352d-05,-9.05168552d-06/
      data fpppp( 3,11),fpppp( 3,12)/-2.61236205d-06,-2.78804267d-05/
      data fpppp( 3,13),fpppp( 3,14)/-1.15630234d-04, 1.46574841d-04/
      data fpppp( 3,15),fpppp( 3,16)/ 4.78296616d-04,-4.09702684d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 3.33838869d-04,-7.18119290d-04/
      data fpppp( 3,19) /            -1.25707753d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-4.12973453d-05,-1.77385929d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/ 9.60733463d-06, 2.36981493d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/ 2.57515588d-05, 6.60138079d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.30583865d-05,-8.72059428d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-1.79771543d-05,-1.64816671d-05/
      data fpppp( 4,11),fpppp( 4,12)/-2.70682074d-07, 2.14766702d-06/
      data fpppp( 4,13),fpppp( 4,14)/-1.60757483d-04, 2.12963579d-04/
      data fpppp( 4,15),fpppp( 4,16)/ 1.54940892d-05, 5.92962324d-05/
      data fpppp( 4,17),fpppp( 4,18)/ 5.86068720d-05,-1.44468602d-04/
      data fpppp( 4,19) /            -2.89719472d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.78152572d-05,-1.21695142d-05/
      data fpppp( 5, 3),fpppp( 5, 4)/ 3.49979930d-06, 1.90866162d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/ 1.18040686d-05,-2.16416058d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 4.48268109d-06,-1.25266889d-05/
      data fpppp( 5, 9),fpppp( 5,10)/-1.46867361d-05, 1.76078391d-06/
      data fpppp( 5,11),fpppp( 5,12)/-1.21974800d-05,-1.12403966d-05/
      data fpppp( 5,13),fpppp( 5,14)/ 3.26096378d-05,-2.50634449d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 3.52833760d-04,-7.60350002d-06/
      data fpppp( 5,17),fpppp( 5,18)/ 1.34385410d-05,-5.45271427d-05/
      data fpppp( 5,19) /            -9.70620478d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/-3.61345533d-06,-1.01569134d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 4.27866147d-06, 1.34953418d-07/
      data fpppp( 6, 5),fpppp( 6, 6)/ 5.30870277d-06, 1.06895504d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/-1.77500292d-05, 1.05518262d-05/
      data fpppp( 6, 9),fpppp( 6,10)/-1.00481154d-05,-2.60450880d-05/
      data fpppp( 6,11),fpppp( 6,12)/ 3.15929436d-06, 8.03077014d-06/
      data fpppp( 6,13),fpppp( 6,14)/-2.49831638d-05, 1.38941747d-04/
      data fpppp( 6,15),fpppp( 6,16)/-1.52837156d-04, 8.03583528d-05/
      data fpppp( 6,17),fpppp( 6,18)/ 4.22046528d-05,-4.06381668d-05/
      data fpppp( 6,19) /            -1.10712666d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.02814580d-05,-4.85864739d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-6.28200973d-07, 7.77552016d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 4.26307583d-06,-7.25181305d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.27745683d-05,-1.32193745d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-1.14290043d-06, 3.47071980d-06/
      data fpppp( 7,11),fpppp( 7,12)/-1.37422061d-05,-6.81980079d-06/
      data fpppp( 7,13),fpppp( 7,14)/ 5.50999355d-06,-3.95433257d-05/
      data fpppp( 7,15),fpppp( 7,16)/ 6.37950385d-05, 8.08418699d-06/
      data fpppp( 7,17),fpppp( 7,18)/ 3.95262830d-05,-2.11200288d-05/
      data fpppp( 7,19) /            -5.54553669d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/-8.85438375d-07,-2.60267682d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 7.01062396d-07, 4.85834625d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.44599291d-07, 3.91641167d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/-5.21668897d-06, 2.08874103d-06/
      data fpppp( 8, 9),fpppp( 8,10)/-3.18811509d-06,-5.68153154d-06/
      data fpppp( 8,11),fpppp( 8,12)/-1.27676239d-07,-5.28700095d-06/
      data fpppp( 8,13),fpppp( 8,14)/-2.06586803d-06, 1.71472206d-05/
      data fpppp( 8,15),fpppp( 8,16)/-7.20459840d-06, 3.05476376d-05/
      data fpppp( 8,17),fpppp( 8,18)/-8.59513756d-06,-4.07104592d-06/
      data fpppp( 8,19) /            -1.20232019d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.08750856d-06,-5.18841757d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/ 1.04910939d-07, 5.35863018d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/ 7.38680769d-07, 1.88296639d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/ 6.69380608d-07,-1.38102523d-06/
      data fpppp( 9, 9),fpppp( 9,10)/-6.44089920d-07,-1.58935527d-06/
      data fpppp( 9,11),fpppp( 9,12)/-3.02859172d-06,-5.73422535d-07/
      data fpppp( 9,13),fpppp( 9,14)/-1.16011011d-06,-5.45974825d-07/
      data fpppp( 9,15),fpppp( 9,16)/ 1.03146163d-05, 5.66063600d-06/
      data fpppp( 9,17),fpppp( 9,18)/ 8.02951285d-06,-1.01121435d-05/
      data fpppp( 9,19) /            -2.23936474d-05 /
      data fpppp(10, 1),fpppp(10, 2)/-2.20412521d-08, 3.23135935d-09/
      data fpppp(10, 3),fpppp(10, 4)/-3.22668849d-08, 1.44932866d-07/
      data fpppp(10, 5),fpppp(10, 6)/ 4.89039873d-08, 6.48821792d-08/
      data fpppp(10, 7),fpppp(10, 8)/-2.19351304d-07, 8.53440168d-08/
      data fpppp(10, 9),fpppp(10,10)/-4.12674141d-07,-2.61801490d-07/
      data fpppp(10,11),fpppp(10,12)/-2.84852973d-07,-9.67733813d-07/
      data fpppp(10,13),fpppp(10,14)/-4.22618080d-08, 1.35392074d-06/
      data fpppp(10,15),fpppp(10,16)/ 7.91550234d-07, 2.65826705d-06/
      data fpppp(10,17),fpppp(10,18)/-1.13604505d-06,-3.37739461d-07/
      data fpppp(10,19) /            -8.51610354d-07 /
      data fpppp(11, 1),fpppp(11, 2)/-5.45158061d-08,-2.77049287d-08/
      data fpppp(11, 3),fpppp(11, 4)/ 4.08309680d-08,-3.66707027d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 6.13337990d-08,-1.18042287d-08/
      data fpppp(11, 7),fpppp(11, 8)/ 1.81774843d-08,-8.75164923d-08/
      data fpppp(11, 9),fpppp(11,10)/-2.67037840d-08,-5.43120462d-08/
      data fpppp(11,11),fpppp(11,12)/-1.71013012d-07, 8.32975678d-08/
      data fpppp(11,13),fpppp(11,14)/-2.15585145d-07,-1.76779628d-08/
      data fpppp(11,15),fpppp(11,16)/ 5.87804522d-07, 2.31778631d-07/
      data fpppp(11,17),fpppp(11,18)/ 4.76114275d-07,-5.68169254d-07/
      data fpppp(11,19) /            -1.24427588d-06 /
      data fpppp(12, 1),fpppp(12, 2)/-3.54788788d-09, 8.93083916d-10/
      data fpppp(12, 3),fpppp(12, 4)/-6.23537242d-10, 1.87114169d-08/
      data fpppp(12, 5),fpppp(12, 6)/ 1.14814783d-08, 2.49061704d-09/
      data fpppp(12, 7),fpppp(12, 8)/-2.37028201d-08,-1.00571798d-08/
      data fpppp(12, 9),fpppp(12,10)/-4.70500070d-08,-4.00140575d-08/
      data fpppp(12,11),fpppp(12,12)/-6.03007670d-08,-9.35695689d-08/
      data fpppp(12,13),fpppp(12,14)/-1.37393803d-08, 1.30271285d-07/
      data fpppp(12,15),fpppp(12,16)/ 1.85652736d-07, 2.99454020d-07/
      data fpppp(12,17),fpppp(12,18)/-1.61754789d-08,-1.43365399d-07/
      data fpppp(12,19) /            -3.24395199d-07 /
      data fpppp(13, 1),fpppp(13, 2)/ 6.08388361d-08, 3.12986738d-08/
      data fpppp(13, 3),fpppp(13, 4)/-3.19839865d-08, 2.28320964d-08/
      data fpppp(13, 5),fpppp(13, 6)/-6.01962036d-08,-9.86125562d-09/
      data fpppp(13, 7),fpppp(13, 8)/ 1.82706629d-08, 3.22175255d-08/
      data fpppp(13, 9),fpppp(13,10)/ 9.51000081d-08, 7.55180746d-08/
      data fpppp(13,11),fpppp(13,12)/ 1.62531195d-07, 1.05250491d-07/
      data fpppp(13,13),fpppp(13,14)/ 2.21126053d-07,-2.57626800d-07/
      data fpppp(13,15),fpppp(13,16)/-5.29368102d-07,-6.87568916d-07/
      data fpppp(13,17),fpppp(13,18)/-5.37529012d-08, 4.11387168d-07/
      data fpppp(13,19) /             9.24720365d-07 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine s3_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -1.92100000d-04 /
      data f( 1, 3),f( 1, 4) / -1.18590000d-03 , -2.73780000d-03 /
      data f( 1, 5),f( 1, 6) / -3.86470000d-03 , -3.65360000d-03 /
      data f( 1, 7),f( 1, 8) / -2.53520000d-03 , -1.39260000d-03 /
      data f( 1, 9),f( 1,10) / -4.21300000d-04 ,  1.78000000d-04 /
      data f( 1,11),f( 1,12) /  7.98300000d-04 ,  1.74330000d-03 /
      data f( 1,13),f( 1,14) /  3.58320000d-03 ,  7.22500000d-03 /
      data f( 1,15),f( 1,16) /  1.52897000d-02 ,  1.58519000d-02 /
      data f( 1,17),f( 1,18) /  6.85040000d-03 ,  8.80600000d-04 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.50700000d-04 /
      data f( 2, 3),f( 2, 4) / -9.19400000d-04 , -2.05680000d-03 /
      data f( 2, 5),f( 2, 6) / -2.78750000d-03 , -2.61730000d-03 /
      data f( 2, 7),f( 2, 8) / -1.87950000d-03 , -1.06130000d-03 /
      data f( 2, 9),f( 2,10) / -3.14000000d-04 ,  1.61200000d-04 /
      data f( 2,11),f( 2,12) /  5.35700000d-04 ,  1.18070000d-03 /
      data f( 2,13),f( 2,14) /  2.55000000d-03 ,  4.99530000d-03 /
      data f( 2,15),f( 2,16) /  8.97760000d-03 ,  9.55600000d-03 /
      data f( 2,17),f( 2,18) /  5.85820000d-03 ,  2.48310000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.13000000d-04 /
      data f( 3, 3),f( 3, 4) / -6.97300000d-04 , -1.54140000d-03 /
      data f( 3, 5),f( 3, 6) / -2.05180000d-03 , -1.92890000d-03 /
      data f( 3, 7),f( 3, 8) / -1.41340000d-03 , -8.02500000d-04 /
      data f( 3, 9),f( 3,10) / -2.27700000d-04 ,  1.27200000d-04 /
      data f( 3,11),f( 3,12) /  3.49600000d-04 ,  7.49600000d-04 /
      data f( 3,13),f( 3,14) /  1.70960000d-03 ,  3.40260000d-03 /
      data f( 3,15),f( 3,16) /  5.65090000d-03 ,  6.05830000d-03 /
      data f( 3,17),f( 3,18) /  4.37710000d-03 ,  1.83820000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -7.22000000d-05 /
      data f( 4, 3),f( 4, 4) / -4.47000000d-04 , -9.86400000d-04 /
      data f( 4, 5),f( 4, 6) / -1.30850000d-03 , -1.24120000d-03 /
      data f( 4, 7),f( 4, 8) / -9.18900000d-04 , -5.13500000d-04 /
      data f( 4, 9),f( 4,10) / -1.30300000d-04 ,  9.84000000d-05 /
      data f( 4,11),f( 4,12) /  1.79200000d-04 ,  3.25700000d-04 /
      data f( 4,13),f( 4,14) /  8.19600000d-04 ,  1.78870000d-03 /
      data f( 4,15),f( 4,16) /  2.91170000d-03 ,  3.19640000d-03 /
      data f( 4,17),f( 4,18) /  2.33540000d-03 ,  8.53300000d-04 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -3.20000000d-05 /
      data f( 5, 3),f( 5, 4) / -1.99700000d-04 , -4.46600000d-04 /
      data f( 5, 5),f( 5, 6) / -5.99000000d-04 , -5.73900000d-04 /
      data f( 5, 7),f( 5, 8) / -4.22100000d-04 , -2.18500000d-04 /
      data f( 5, 9),f( 5,10) / -3.60000000d-05 ,  6.85000000d-05 /
      data f( 5,11),f( 5,12) /  6.06000000d-05 ,  7.23000000d-05 /
      data f( 5,13),f( 5,14) /  1.13800000d-04 ,  4.29100000d-04 /
      data f( 5,15),f( 5,16) / -1.34400000d-04 ,  5.75400000d-04 /
      data f( 5,17),f( 5,18) /  7.38900000d-04 ,  2.38200000d-04 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -1.77000000d-05 /
      data f( 6, 3),f( 6, 4) / -1.07400000d-04 , -2.17800000d-04 /
      data f( 6, 5),f( 6, 6) / -2.67300000d-04 , -2.47000000d-04 /
      data f( 6, 7),f( 6, 8) / -1.89100000d-04 , -7.31000000d-05 /
      data f( 6, 9),f( 6,10) /  5.18000000d-05 ,  8.06000000d-05 /
      data f( 6,11),f( 6,12) /  4.32000000d-05 ,  2.31000000d-05 /
      data f( 6,13),f( 6,14) /  1.50000000d-05 ,  1.12900000d-04 /
      data f( 6,15),f( 6,16) /  3.09800000d-04 ,  4.42600000d-04 /
      data f( 6,17),f( 6,18) /  4.12300000d-04 ,  1.60700000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -4.90000000d-06 /
      data f( 7, 3),f( 7, 4) / -3.12000000d-05 , -7.06000000d-05 /
      data f( 7, 5),f( 7, 6) / -9.62000000d-05 , -9.22000000d-05 /
      data f( 7, 7),f( 7, 8) / -5.50000000d-05 ,  1.20000000d-06 /
      data f( 7, 9),f( 7,10) /  2.50000000d-05 ,  2.78000000d-05 /
      data f( 7,11),f( 7,12) /  4.60000000d-06 , -3.29000000d-05 /
      data f( 7,13),f( 7,14) / -7.02000000d-05 , -7.25000000d-05 /
      data f( 7,15),f( 7,16) /  3.12000000d-05 ,  6.04000000d-05 /
      data f( 7,17),f( 7,18) /  5.53000000d-05 ,  1.33000000d-05 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -1.90000000d-06 /
      data f( 8, 3),f( 8, 4) / -1.17000000d-05 , -2.60000000d-05 /
      data f( 8, 5),f( 8, 6) / -3.44000000d-05 , -2.86000000d-05 /
      data f( 8, 7),f( 8, 8) /  2.60000000d-06 ,  1.83000000d-05 /
      data f( 8, 9),f( 8,10) /  2.41000000d-05 ,  2.00000000d-05 /
      data f( 8,11),f( 8,12) / -3.50000000d-06 , -2.63000000d-05 /
      data f( 8,13),f( 8,14) / -5.53000000d-05 , -6.26000000d-05 /
      data f( 8,15),f( 8,16) / -3.90000000d-05 , -4.30000000d-06 /
      data f( 8,17),f( 8,18) /  8.90000000d-06 ,  1.56000000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -7.00000000d-07 /
      data f( 9, 3),f( 9, 4) / -3.90000000d-06 , -8.10000000d-06 /
      data f( 9, 5),f( 9, 6) / -9.40000000d-06 , -1.20000000d-06 /
      data f( 9, 7),f( 9, 8) /  1.42000000d-05 ,  2.28000000d-05 /
      data f( 9, 9),f( 9,10) /  2.21000000d-05 ,  1.70000000d-05 /
      data f( 9,11),f( 9,12) / -6.60000000d-06 , -2.72000000d-05 /
      data f( 9,13),f( 9,14) / -4.30000000d-05 , -4.72000000d-05 /
      data f( 9,15),f( 9,16) / -3.44000000d-05 , -1.37000000d-05 /
      data f( 9,17),f( 9,18) / -7.00000000d-07 ,  7.00000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(10, 3),f(10, 4) /  7.00000000d-07 ,  1.90000000d-06 /
      data f(10, 5),f(10, 6) /  5.00000000d-06 ,  9.40000000d-06 /
      data f(10, 7),f(10, 8) /  1.40000000d-05 ,  1.64000000d-05 /
      data f(10, 9),f(10,10) /  1.39000000d-05 ,  8.00000000d-06 /
      data f(10,11),f(10,12) / -3.70000000d-06 , -1.43000000d-05 /
      data f(10,13),f(10,14) / -2.24000000d-05 , -2.58000000d-05 /
      data f(10,15),f(10,16) / -2.33000000d-05 , -1.44000000d-05 /
      data f(10,17),f(10,18) / -5.30000000d-06 , -6.00000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(11, 3),f(11, 4) /  7.00000000d-07 ,  2.00000000d-06 /
      data f(11, 5),f(11, 6) /  3.70000000d-06 ,  5.00000000d-06 /
      data f(11, 7),f(11, 8) /  6.10000000d-06 ,  6.70000000d-06 /
      data f(11, 9),f(11,10) /  6.20000000d-06 ,  4.30000000d-06 /
      data f(11,11),f(11,12) /  1.10000000d-06 , -3.30000000d-06 /
      data f(11,13),f(11,14) / -6.20000000d-06 , -1.11000000d-05 /
      data f(11,15),f(11,16) / -1.27000000d-05 , -9.70000000d-06 /
      data f(11,17),f(11,18) / -4.20000000d-06 , -6.00000000d-07 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) /  3.00000000d-07 ,  1.00000000d-06 /
      data f(12, 5),f(12, 6) /  1.70000000d-06 ,  2.10000000d-06 /
      data f(12, 7),f(12, 8) /  2.20000000d-06 ,  2.40000000d-06 /
      data f(12, 9),f(12,10) /  2.50000000d-06 ,  2.40000000d-06 /
      data f(12,11),f(12,12) /  1.60000000d-06 ,  1.00000000d-07 /
      data f(12,13),f(12,14) / -1.50000000d-06 , -3.40000000d-06 /
      data f(12,15),f(12,16) / -5.00000000d-06 , -4.40000000d-06 /
      data f(12,17),f(12,18) / -2.10000000d-06 , -4.00000000d-07 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  2.00000000d-07 ,  5.00000000d-07 /
      data f(13, 5),f(13, 6) /  9.00000000d-07 ,  1.00000000d-06 /
      data f(13, 7),f(13, 8) /  1.00000000d-06 ,  9.00000000d-07 /
      data f(13, 9),f(13,10) /  7.00000000d-07 ,  6.00000000d-07 /
      data f(13,11),f(13,12) /  3.00000000d-07 , -1.00000000d-07 /
      data f(13,13),f(13,14) / -6.00000000d-07 , -1.00000000d-06 /
      data f(13,15),f(13,16) / -1.20000000d-06 , -1.00000000d-06 /
      data f(13,17),f(13,18) / -3.00000000d-07 , -1.00000000d-07 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-1.07325140d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-2.17420196d-06,-7.45797207d-06/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-1.12786018d-03,-7.53759774d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-5.17788163d-03, 4.12236303d-06/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.16431632d-02, 1.65481456d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-1.20996348d-02, 9.96505455d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-6.41043164d-03,-1.97036382d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-2.22335063d-03,-6.31599281d-07/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-5.98598920d-04,-5.78123906d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/-9.96176974d-04, 1.43655551d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 2.31460937d-03, 1.29501702d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 3.55854656d-03, 1.28653764d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 4.76765166d-03, 9.37477281d-07/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.05984361d-02, 9.14987145d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.12015146d-01,-1.01558335d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/ 1.01878978d-01,-1.35415374d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/-2.44035658d-02, 6.93978315d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.01475231d-01, 3.97260481d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 7.70499759d-05/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-8.89970902d-06/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-7.17230246d-05,-5.72658196d-06/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-1.08342251d-03,-5.27396314d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-4.17623674d-03, 4.70043450d-06/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-8.72760227d-03, 1.08742251d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-8.91980189d-03, 5.85666502d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-4.84713671d-03,-2.44885183d-07/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-1.83051302d-03,-5.31242845d-08/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-5.45230730d-04,-3.79661768d-06/
      data fpp( 2,10,1),fpp( 2,10,2)/-5.42217480d-04,-1.08640500d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.90113840d-03, 2.10023767d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.28969260d-03, 8.91545430d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/ 4.71669668d-03, 5.69594511d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.61604135d-02, 3.28607653d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/ 7.63724232d-02,-4.49190062d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 7.15236150d-02,-5.74187406d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/-1.29263685d-02, 1.80219686d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-6.15205387d-02, 4.69286611d-06/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 1.67265669d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-7.00377846d-06/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-2.65933699d-04,-4.36144308d-06/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-1.19844980d-03,-3.82844921d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-2.95717142d-03, 4.08723991d-06/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-4.67142777d-03, 7.50148958d-06/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-4.40615768d-03, 3.90480178d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-2.64102151d-03, 4.35303305d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-1.32959728d-03, 7.79849998d-08/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-3.70478158d-04,-2.91324330d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/ 5.85046895d-04,-1.61901178d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 1.55583704d-03, 1.43929043d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.00768305d-03, 6.51785006d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/ 5.28556162d-03, 6.08930931d-06/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.03099099d-02, 1.31049127d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.03051617d-02,-2.51909600d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 3.17565619d-02,-2.27950726d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/ 2.77403967d-03,-8.94474969d-06/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.04473856d-02, 7.11207134d-06/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 2.25384643d-05/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-4.52175248d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.15738985d-04,-2.79449505d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-8.06218990d-04,-2.45626734d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.89860410d-03, 2.74356441d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-2.62683927d-03, 4.52000970d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-2.35960649d-03, 2.54039679d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-1.60850382d-03, 6.18403140d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-9.61000393d-04,-2.80093485d-08/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-5.38252321d-04,-1.83836575d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.08677995d-04,-1.88852767d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/ 7.96450928d-04, 5.18476421d-07/
      data fpp( 4,12,1),fpp( 4,12,2)/ 2.63126143d-03, 3.75662198d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.94366349d-03, 5.29903564d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/ 6.53669129d-03, 3.55923545d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/-1.87548770d-03,-1.03019775d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/ 5.43905046d-03,-1.26493256d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.13674468d-02,-7.84272002d-06/
      data fpp( 4,18,1),fpp( 4,18,2)/ 5.01907388d-03, 6.75420572d-06/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 1.85538971d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-1.97685345d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.37275027d-04,-1.24729311d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-7.77829350d-04,-1.17597412d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-1.39496403d-03, 1.19918957d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-1.49525768d-03, 2.04921582d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-1.29836463d-03, 1.25394714d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-1.12497486d-03, 5.36995614d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-6.07040377d-04,-2.93929598d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 3.11894322d-04,-6.27277221d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/ 4.31141448d-04,-1.87696152d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 4.87454805d-04, 1.39112329d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/ 6.49753597d-04,-2.51153163d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 2.86973986d-03, 1.04430032d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 4.82224191d-03,-2.28324812d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 2.42800636d-02, 2.81589218d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 1.51131014d-02,-1.34052058d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 5.31174656d-03,-7.31609843d-06/
      data fpp( 5,18,1),fpp( 5,18,2)/ 2.30413223d-03, 2.81759955d-06/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 1.17957002d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-1.25189147d-06/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 4.32390941d-05,-6.82217065d-07/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.97536390d-04,-3.39240274d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.44602276d-05, 7.97178161d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-4.59330001d-04, 8.04527630d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-6.16534990d-04, 1.72711318d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-2.22796742d-04, 7.60627098d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-2.01238098d-04, 2.70780290d-07/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-8.65324966d-04,-1.30974826d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-6.07887797d-04,-7.97787253d-07/
      data fpp( 6,11,1),fpp( 6,11,2)/-3.17470150d-04, 5.28897271d-07/
      data fpp( 6,12,1),fpp( 6,12,2)/-3.29475817d-04,-2.79801830d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/-8.54622912d-04, 1.31031005d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/-7.84058948d-04, 1.39856163d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.14775669d-02,-9.64556569d-07/
      data fpp( 6,16,1),fpp( 6,16,2)/-6.17465608d-03,-1.38633535d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/-2.13683301d-03,-3.27610202d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.33320282d-03, 1.21274343d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 3.87912828d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-3.05603208d-07/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-7.16813490d-05,-1.96793583d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-3.98716211d-04,-1.91222459d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-6.21276878d-04, 1.75683419d-07/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-5.21822312d-04, 3.16488784d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-3.65895408d-04, 3.34361446d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-3.57438173d-04, 3.38065432d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-2.94407230d-04,-5.46623175d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 3.99005541d-04,-9.55727333d-08/
      data fpp( 7,10,1),fpp( 7,10,2)/ 4.42809741d-04,-3.31085892d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.73625793d-04,-1.40083698d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 5.04949672d-04, 3.34206837d-08/
      data fpp( 7,13,1),fpp( 7,13,2)/ 8.75151794d-04, 1.84009632d-08/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.45319388d-03, 1.99297546d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 4.28300386d-03,-1.63030282d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 3.59992293d-03, 5.82358068d-08/
      data fpp( 7,17,1),fpp( 7,17,2)/ 2.50598547d-03,-6.60640410d-07/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.35107903d-03, 3.70325831d-07/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 9.01337084d-07/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-1.16736137d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 8.28630211d-06,-7.15277270d-08/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 3.65284537d-05,-7.11529556d-08/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 8.24728474d-06, 8.61395495d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-7.65807515d-05, 8.05947578d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.08683377d-04, 4.43481419d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.83450566d-04,-3.30520435d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 6.06701914d-06,-5.13996791d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.09097198d-04,-5.78808485d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/-8.33511679d-05,-3.11076927d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/-4.50330235d-05, 1.38188556d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-1.87922872d-04,-1.99677297d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-2.43584265d-04, 2.88520632d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/-3.41516559d-04, 3.47594768d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.52848582d-04, 1.75100297d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/-6.05035650d-04,-3.81995956d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/-4.32708889d-04, 6.28835266d-08/
      data fpp( 8,18,1),fpp( 8,18,2)/-4.78313312d-04,-2.59538150d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-3.62730925d-07/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-4.05391110d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-4.66385941d-06,-2.39217781d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-2.81976036d-05,-1.37737767d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-5.25122609d-05, 1.90168847d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-5.50546821d-05, 1.11706238d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-6.81710830d-05, 1.04158164d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-1.27595613d-05,-9.63388941d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-3.22608463d-05,-1.26802588d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.09832501d-05, 4.55492452d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/ 5.79493020d-06,-3.19394393d-07/
      data fpp( 9,11,1),fpp( 9,11,2)/ 2.65063007d-05, 1.22028327d-07/
      data fpp( 9,12,1),fpp( 9,12,2)/ 6.67418141d-05, 1.12810849d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/ 3.67852650d-05, 1.20847333d-07/
      data fpp( 9,14,1),fpp( 9,14,2)/ 4.48723601d-05, 2.01329582d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/ 1.23590462d-04, 9.38343385d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 1.47419668d-04,-1.02666936d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/ 1.08050080d-04,-1.45166594d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.49374214d-04,-1.26666874d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 6.98333437d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 3.30014842d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/ 2.48427165d-07, 6.39970316d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/ 3.28584071d-07, 1.10103895d-09/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.38685977d-06, 2.51961410d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.01455780d-05, 1.21143969d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-6.34506229d-06, 4.34627127d-09/
      data fpp(10, 7,1),fpp(10, 7,2)/-1.03960329d-05,-1.74994820d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.34902929d-06,-6.63483433d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-3.60115144d-06,-1.11071449d-08/
      data fpp(10,10,1),fpp(10,10,2)/ 6.29079333d-06,-9.32230772d-08/
      data fpp(10,11,1),fpp(10,11,2)/-2.40239033d-06, 3.59994538d-08/
      data fpp(10,12,1),fpp(10,12,2)/-1.80640066d-05, 1.52252619d-08/
      data fpp(10,13,1),fpp(10,13,2)/-1.25636625d-05, 5.30994984d-08/
      data fpp(10,14,1),fpp(10,14,2)/-2.02588006d-05, 5.43767443d-08/
      data fpp(10,15,1),fpp(10,15,2)/-3.29470955d-05, 8.33935243d-08/
      data fpp(10,16,1),fpp(10,16,2)/-3.11411788d-05,-3.95084157d-09/
      data fpp(10,17,1),fpp(10,17,2)/-2.01957961d-05,-5.55901580d-08/
      data fpp(10,18,1),fpp(10,18,2)/-3.79659869d-05,-3.76885263d-08/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-3.96557369d-08/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-1.70441256d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/-5.29849254d-07, 2.40882512d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-7.16732650d-07, 1.00691121d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.34030006d-06, 5.31472653d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.43699422d-06,-7.32801821d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/ 3.55133218d-06,-2.65370821d-12/
      data fpp(11, 7,1),fpp(11, 7,2)/ 8.14369287d-06,-4.66136696d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 7.06472913d-06,-1.13518784d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/ 6.42135564d-06,-1.59311193d-08/
      data fpp(11,10,1),fpp(11,10,2)/ 8.41896493d-07,-8.92364451d-09/
      data fpp(11,11,1),fpp(11,11,2)/-5.49673939d-06,-2.63743027d-08/
      data fpp(11,12,1),fpp(11,12,2)/-5.88578774d-06, 4.24208553d-08/
      data fpp(11,13,1),fpp(11,13,2)/-1.29306151d-05,-5.33091185d-08/
      data fpp(11,14,1),fpp(11,14,2)/-4.03715774d-06, 5.08156187d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 5.19791987d-06, 4.80466436d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 9.54504704d-06, 3.29978069d-08/
      data fpp(11,17,1),fpp(11,17,2)/ 6.93310424d-06,-3.00378711d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 1.02897334d-05,-2.68463225d-08/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-4.25768387d-08/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 2.21403987d-09/
      data fpp(12, 2,1),fpp(12, 2,2)/ 7.09698508d-08, 2.57192026d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.38346530d-07, 5.49827910d-09/
      data fpp(12, 4,1),fpp(12, 4,2)/ 1.48060012d-07,-5.65036662d-10/
      data fpp(12, 5,1),fpp(12, 5,2)/ 1.97601156d-07,-3.23813245d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/ 1.13973356d-06,-4.48243352d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/ 1.82126143d-06, 3.16786654d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/ 2.79205417d-06,-2.18903265d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 1.91572887d-06,-4.11735960d-10/
      data fpp(12,10,1),fpp(12,10,2)/ 1.14162070d-06,-8.16402351d-09/
      data fpp(12,11,1),fpp(12,11,2)/-1.41065212d-06,-8.93216999d-09/
      data fpp(12,12,1),fpp(12,12,2)/-3.99284245d-06, 1.89270346d-09/
      data fpp(12,13,1),fpp(12,13,2)/-4.71387697d-06,-4.63864387d-09/
      data fpp(12,14,1),fpp(12,14,2)/-5.59256845d-06,-1.33812799d-09/
      data fpp(12,15,1),fpp(12,15,2)/-5.24458397d-06, 2.79911558d-08/
      data fpp(12,16,1),fpp(12,16,2)/-3.43900941d-06, 2.13735046d-08/
      data fpp(12,17,1),fpp(12,17,2)/-1.53662085d-06,-1.14851743d-08/
      data fpp(12,18,1),fpp(12,18,2)/-1.99294668d-06,-1.14328073d-08/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-2.07835963d-08/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 2.89487799d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/ 6.52015075d-07, 2.21024402d-09/
      data fpp(13, 3,1),fpp(13, 3,2)/ 9.93326735d-07, 2.64145932d-10/
      data fpp(13, 4,1),fpp(13, 4,2)/ 2.47596999d-06, 2.73317225d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/ 3.48869942d-06,-5.19683495d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.85513322d-06, 5.41675314d-11/
      data fpp(13, 7,1),fpp(13, 7,2)/ 3.64369287d-07,-1.01983518d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.25852709d-06,-1.97482681d-09/
      data fpp(13, 9,1),fpp(13, 9,2)/-5.57864436d-07, 2.91914243d-09/
      data fpp(13,10,1),fpp(13,10,2)/-8.45810351d-07,-3.70174292d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 3.53032606d-06,-1.12170756d-10/
      data fpp(13,12,1),fpp(13,12,2)/ 4.42142123d-06,-1.84957406d-09/
      data fpp(13,13,1),fpp(13,13,2)/ 7.85693849d-06, 1.51046699d-09/
      data fpp(13,14,1),fpp(13,14,2)/-7.03715774d-07, 1.80770611d-09/
      data fpp(13,15,1),fpp(13,15,2)/-4.26520801d-06, 3.25870858d-09/
      data fpp(13,16,1),fpp(13,16,2)/-5.25549530d-06, 9.15745955d-09/
      data fpp(13,17,1),fpp(13,17,2)/-2.45668958d-06,-9.88854680d-09/
      data fpp(13,18,1),fpp(13,18,2)/ 6.83973340d-07, 3.96727658d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 2.30163617d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 6.86472767d-06,-1.14183343d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/-2.86020968d-05,-4.96334073d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/ 8.22201222d-05, 8.12815125d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/-3.86056874d-05,-1.69860895d-05/
      data fpppp( 1, 9),fpppp( 1,10)/-4.71897126d-05, 8.44051540d-05/
      data fpppp( 1,11),fpppp( 1,12)/-6.79290393d-05, 6.33000533d-05/
      data fpppp( 1,13),fpppp( 1,14)/-1.87361099d-04, 1.56344510d-03/
      data fpppp( 1,15),fpppp( 1,16)/-1.53126382d-03,-1.53156245d-03/
      data fpppp( 1,17),fpppp( 1,18)/ 6.88731004d-04, 1.72929118d-03/
      data fpppp( 1,19) /             3.10691801d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 1.98166076d-06,-9.33984934d-06/
      data fpppp( 2, 3),fpppp( 2, 4)/-2.10208508d-05,-3.14436325d-05/
      data fpppp( 2, 5),fpppp( 2, 6)/ 5.92823024d-05, 5.58643779d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/-2.68479263d-05,-1.18351622d-05/
      data fpppp( 2, 9),fpppp( 2,10)/-2.96919087d-05, 5.36666547d-05/
      data fpppp( 2,11),fpppp( 2,12)/-3.85541524d-05, 3.72618544d-05/
      data fpppp( 2,13),fpppp( 2,14)/-1.08186272d-04, 9.96485998d-04/
      data fpppp( 2,15),fpppp( 2,16)/-9.51660147d-04,-1.09349449d-03/
      data fpppp( 2,17),fpppp( 2,18)/ 5.49567579d-04, 1.04657296d-03/
      data fpppp( 2,19) /             1.87102311d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-5.32252800d-06,-6.15995376d-06/
      data fpppp( 3, 3),fpppp( 3, 4)/-1.00326011d-05,-3.28197287d-06/
      data fpppp( 3, 5),fpppp( 3, 6)/ 2.58284094d-05, 1.87399212d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-1.07961296d-05,-2.77811855d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 7.70297045d-07,-5.18713706d-07/
      data fpppp( 3,11),fpppp( 3,12)/ 2.22046352d-06, 2.05002113d-05/
      data fpppp( 3,13),fpppp( 3,14)/-3.46593556d-05, 2.82925395d-04/
      data fpppp( 3,15),fpppp( 3,16)/-1.98788017d-04,-6.00404422d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 7.74370362d-04,-2.97724938d-04/
      data fpppp( 3,19) /            -6.70714500d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-7.70400461d-06,-5.29093971d-06/
      data fpppp( 4, 3),fpppp( 4, 4)/-5.61669774d-06, 3.64342450d-06/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.28919960d-05, 4.51666852d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.92647689d-06,-3.02671504d-06/
      data fpppp( 4, 9),fpppp( 4,10)/ 5.48015492d-07, 1.24422832d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 2.30083471d-05,-3.74967219d-05/
      data fpppp( 4,13),fpppp( 4,14)/ 9.56340345d-05,-2.68201872d-04/
      data fpppp( 4,15),fpppp( 4,16)/ 3.16861045d-04,-5.56392769d-05/
      data fpppp( 4,17),fpppp( 4,18)/-1.77472449d-04, 2.89229205d-05/
      data fpppp( 4,19) /             1.41538707d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.03571309d-05,-4.91809578d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.67243725d-07, 6.99224910d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/ 3.20870926d-06,-1.99588399d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 3.36462984d-06, 9.21004720d-06/
      data fpppp( 5, 9),fpppp( 5,10)/-1.61448056d-05, 7.38792077d-06/
      data fpppp( 5,11),fpppp( 5,12)/-1.71829037d-05, 6.77028199d-05/
      data fpppp( 5,13),fpppp( 5,14)/-1.30167128d-04, 4.36916641d-04/
      data fpppp( 5,15),fpppp( 5,16)/-5.67180254d-04, 1.14317337d-04/
      data fpppp( 5,17),fpppp( 5,18)/ 7.18473501d-05, 5.91769382d-06/
      data fpppp( 5,19) /            -5.33091997d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 5.45403892d-06, 1.39230144d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/-4.35975256d-06,-4.19569872d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 3.69970347d-06, 8.39199926d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/-4.21110627d-06,-1.38783504d-05/
      data fpppp( 6, 9),fpppp( 6,10)/ 1.85857773d-05,-5.17331676d-06/
      data fpppp( 6,11),fpppp( 6,12)/ 4.08631847d-06,-2.93173561d-05/
      data fpppp( 6,13),fpppp( 6,14)/ 8.23946201d-05,-2.64518461d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 3.29834910d-04,-9.50360541d-05/
      data fpppp( 6,17),fpppp( 6,18)/-2.55959561d-05, 3.36830545d-06/
      data fpppp( 6,19) /             4.38970917d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/-6.14719401d-06,-2.56277152d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/ 1.07706933d-06, 4.52294588d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.52061146d-07,-1.74285023d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.02884035d-06, 1.31326341d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-1.26787862d-05,-1.39400339d-06/
      data fpppp( 7,11),fpppp( 7,12)/ 5.47551087d-06, 3.52242953d-06/
      data fpppp( 7,13),fpppp( 7,14)/-1.12325344d-05, 5.38781057d-05/
      data fpppp( 7,15),fpppp( 7,16)/-6.91738142d-05, 1.20436958d-05/
      data fpppp( 7,17),fpppp( 7,18)/-3.65236088d-06,-1.09239131d-06/
      data fpppp( 7,19) /            -3.74842926d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 9.56628614d-07, 2.15001652d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/-6.19284255d-07,-1.12926386d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.74352764d-06,-2.68132207d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 6.42188684d-06,-7.14913880d-06/
      data fpppp( 8, 9),fpppp( 8,10)/ 3.89376019d-06, 2.87128210d-08/
      data fpppp( 8,11),fpppp( 8,12)/-3.25428461d-06, 2.11594606d-06/
      data fpppp( 8,13),fpppp( 8,14)/ 2.42076683d-08,-4.74903080d-06/
      data fpppp( 8,15),fpppp( 8,16)/ 6.16793184d-06, 1.62600065d-06/
      data fpppp( 8,17),fpppp( 8,18)/-5.20110461d-06, 6.10254671d-06/
      data fpppp( 8,19) /             1.22259819d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-3.83401561d-07,-1.61072858d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.04500095d-07, 5.32218457d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-7.18039570d-07, 1.70550104d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.99228921d-06, 1.76888741d-06/
      data fpppp( 9, 9),fpppp( 9,10)/-1.31853753d-06, 5.99317746d-07/
      data fpppp( 9,11),fpppp( 9,12)/ 4.75247972d-07,-1.32886106d-06/
      data fpppp( 9,13),fpppp( 9,14)/ 6.28672499d-07, 1.09678972d-06/
      data fpppp( 9,15),fpppp( 9,16)/-7.77970959d-07,-1.27823965d-06/
      data fpppp( 9,17),fpppp( 9,18)/ 2.09900195d-06,-2.27614484d-06/
      data fpppp( 9,19) /            -4.43632351d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 2.05616259d-08,-1.36247499d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 2.38411580d-08,-1.89475927d-07/
      data fpppp(10, 5),fpppp(10, 6)/ 3.11466083d-07,-3.02834365d-07/
      data fpppp(10, 7),fpppp(10, 8)/ 4.28782196d-07,-4.64532451d-07/
      data fpppp(10, 9),fpppp(10,10)/ 4.27633035d-07,-3.55472159d-07/
      data fpppp(10,11),fpppp(10,12)/-1.20852105d-07, 4.20774620d-07/
      data fpppp(10,13),fpppp(10,14)/-2.92528750d-07,-4.23885576d-08/
      data fpppp(10,15),fpppp(10,16)/ 1.62493572d-07, 2.62066967d-07/
      data fpppp(10,17),fpppp(10,18)/-6.62393485d-07, 6.64572563d-07/
      data fpppp(10,19) /             1.34827390d-06 /
      data fpppp(11, 1),fpppp(11, 2)/ 8.09213559d-09, 9.69770201d-09/
      data fpppp(11, 3),fpppp(11, 4)/-2.63049922d-08, 6.93212257d-08/
      data fpppp(11, 5),fpppp(11, 6)/-4.69282092d-08, 7.86142318d-08/
      data fpppp(11, 7),fpppp(11, 8)/-1.18847354d-07, 5.64957203d-08/
      data fpppp(11, 9),fpppp(11,10)/-8.10001125d-08,-2.86604095d-08/
      data fpppp(11,11),fpppp(11,12)/ 1.50091146d-07,-2.14728924d-07/
      data fpppp(11,13),fpppp(11,14)/ 3.09477807d-07,-6.68852162d-08/
      data fpppp(11,15),fpppp(11,16)/-2.14397299d-08,-1.40632890d-07/
      data fpppp(11,17),fpppp(11,18)/ 1.66427092d-07,-1.66961161d-07/
      data fpppp(11,19) /            -3.17364203d-07 /
      data fpppp(12, 1),fpppp(12, 2)/ 6.07912231d-10,-2.42219101d-10/
      data fpppp(12, 3),fpppp(12, 4)/ 1.45373878d-10,-3.79906823d-09/
      data fpppp(12, 5),fpppp(12, 6)/ 1.74405587d-08,-1.24076909d-08/
      data fpppp(12, 7),fpppp(12, 8)/ 1.65539320d-08,-3.64521441d-08/
      data fpppp(12, 9),fpppp(12,10)/ 1.84275614d-08,-3.11250737d-08/
      data fpppp(12,11),fpppp(12,12)/-6.17145676d-10, 3.17986061d-08/
      data fpppp(12,13),fpppp(12,14)/-1.49079298d-08, 1.83736956d-08/
      data fpppp(12,15),fpppp(12,16)/ 1.50137049d-08, 9.02689011d-09/
      data fpppp(12,17),fpppp(12,18)/-4.53124255d-08, 3.06999484d-08/
      data fpppp(12,19) /             6.94689825d-08 /
      data fpppp(13, 1),fpppp(13, 2)/-1.66039046d-08,-5.15399348d-09/
      data fpppp(13, 3),fpppp(13, 4)/ 1.85776737d-08,-6.76805222d-10/
      data fpppp(13, 5),fpppp(13, 6)/-4.40652826d-08, 1.81601977d-08/
      data fpppp(13, 7),fpppp(13, 8)/-2.00073719d-08, 5.39413434d-08/
      data fpppp(13, 9),fpppp(13,10)/-5.63444603d-08, 1.12119984d-07/
      data fpppp(13,11),fpppp(13,12)/-1.12290535d-07, 1.27939683d-07/
      data fpppp(13,13),fpppp(13,14)/-2.46802872d-07, 1.39501513d-07/
      data fpppp(13,15),fpppp(13,16)/-1.12534573d-08, 5.97846140d-08/
      data fpppp(13,17),fpppp(13,18)/-5.39418508d-10,-3.71155081d-08/
      data fpppp(13,19) /            -8.04767243d-08 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine s4_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -2.42000000d-05 /
      data f( 1, 3),f( 1, 4) / -2.61500000d-04 , -7.76400000d-04 /
      data f( 1, 5),f( 1, 6) / -1.23740000d-03 , -1.19270000d-03 /
      data f( 1, 7),f( 1, 8) / -7.66500000d-04 , -3.34900000d-04 /
      data f( 1, 9),f( 1,10) / -1.31600000d-04 , -4.81500000d-04 /
      data f( 1,11),f( 1,12) / -1.62450000d-03 , -3.69790000d-03 /
      data f( 1,13),f( 1,14) / -7.04870000d-03 , -1.18911000d-02 /
      data f( 1,15),f( 1,16) / -1.36459000d-02 , -1.06579000d-02 /
      data f( 1,17),f( 1,18) / -2.60490000d-03 ,  1.01910000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.80000000d-05 /
      data f( 2, 3),f( 2, 4) / -1.94200000d-04 , -5.61200000d-04 /
      data f( 2, 5),f( 2, 6) / -8.55000000d-04 , -8.13800000d-04 /
      data f( 2, 7),f( 2, 8) / -5.34500000d-04 , -2.41100000d-04 /
      data f( 2, 9),f( 2,10) / -1.12100000d-04 , -3.97400000d-04 /
      data f( 2,11),f( 2,12) / -1.36890000d-03 , -3.10480000d-03 /
      data f( 2,13),f( 2,14) / -5.79000000d-03 , -9.86760000d-03 /
      data f( 2,15),f( 2,16) / -1.17906000d-02 , -9.03090000d-03 /
      data f( 2,17),f( 2,18) / -2.75910000d-03 ,  1.53080000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.29000000d-05 /
      data f( 3, 3),f( 3, 4) / -1.41100000d-04 , -4.04400000d-04 /
      data f( 3, 5),f( 3, 6) / -6.04000000d-04 , -5.73700000d-04 /
      data f( 3, 7),f( 3, 8) / -3.82900000d-04 , -1.75700000d-04 /
      data f( 3, 9),f( 3,10) / -9.12000000d-05 , -3.31200000d-04 /
      data f( 3,11),f( 3,12) / -1.13240000d-03 , -2.57700000d-03 /
      data f( 3,13),f( 3,14) / -4.70950000d-03 , -7.81280000d-03 /
      data f( 3,15),f( 3,16) / -9.78140000d-03 , -8.01170000d-03 /
      data f( 3,17),f( 3,18) / -3.04310000d-03 ,  5.99800000d-04 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -7.60000000d-06 /
      data f( 4, 3),f( 4, 4) / -8.49000000d-05 , -2.45100000d-04 /
      data f( 4, 5),f( 4, 6) / -3.65800000d-04 , -3.54700000d-04 /
      data f( 4, 7),f( 4, 8) / -2.39800000d-04 , -1.09200000d-04 /
      data f( 4, 9),f( 4,10) / -6.00000000d-05 , -2.39900000d-04 /
      data f( 4,11),f( 4,12) / -8.24100000d-04 , -1.90070000d-03 /
      data f( 4,13),f( 4,14) / -3.42410000d-03 , -5.34880000d-03 /
      data f( 4,15),f( 4,16) / -6.79810000d-03 , -5.85180000d-03 /
      data f( 4,17),f( 4,18) / -2.60230000d-03 , -6.96000000d-05 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -3.20000000d-06 /
      data f( 5, 3),f( 5, 4) / -3.60000000d-05 , -1.05400000d-04 /
      data f( 5, 5),f( 5, 6) / -1.62200000d-04 , -1.63900000d-04 /
      data f( 5, 7),f( 5, 8) / -1.12700000d-04 , -4.55000000d-05 /
      data f( 5, 9),f( 5,10) / -1.97000000d-05 , -1.19400000d-04 /
      data f( 5,11),f( 5,12) / -4.52800000d-04 , -1.05770000d-03 /
      data f( 5,13),f( 5,14) / -1.93400000d-03 , -2.87700000d-03 /
      data f( 5,15),f( 5,16) / -3.09950000d-03 , -2.80140000d-03 /
      data f( 5,17),f( 5,18) / -1.36650000d-03 , -1.29900000d-04 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -1.27000000d-05 /
      data f( 6, 3),f( 6, 4) / -3.87000000d-05 , -6.49000000d-05 /
      data f( 6, 5),f( 6, 6) / -7.98000000d-05 , -7.69000000d-05 /
      data f( 6, 7),f( 6, 8) / -6.22000000d-05 , -2.36000000d-05 /
      data f( 6, 9),f( 6,10) /  7.61000000d-05 , -5.24000000d-05 /
      data f( 6,11),f( 6,12) / -1.64800000d-04 , -4.56200000d-04 /
      data f( 6,13),f( 6,14) / -9.07100000d-04 , -1.41760000d-03 /
      data f( 6,15),f( 6,16) / -1.71080000d-03 , -1.47240000d-03 /
      data f( 6,17),f( 6,18) / -7.58800000d-04 , -3.77700000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -5.00000000d-07 /
      data f( 7, 3),f( 7, 4) / -5.50000000d-06 , -1.79000000d-05 /
      data f( 7, 5),f( 7, 6) / -3.03000000d-05 , -3.20000000d-05 /
      data f( 7, 7),f( 7, 8) / -1.89000000d-05 ,  1.36000000d-05 /
      data f( 7, 9),f( 7,10) /  1.68000000d-05 , -3.07000000d-05 /
      data f( 7,11),f( 7,12) / -1.15900000d-04 , -2.68900000d-04 /
      data f( 7,13),f( 7,14) / -4.85500000d-04 , -7.10800000d-04 /
      data f( 7,15),f( 7,16) / -7.90900000d-04 , -6.14600000d-04 /
      data f( 7,17),f( 7,18) / -2.37500000d-04 , -2.17000000d-05 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f( 8, 3),f( 8, 4) / -2.40000000d-06 , -7.90000000d-06 /
      data f( 8, 5),f( 8, 6) / -1.28000000d-05 , -1.15000000d-05 /
      data f( 8, 7),f( 8, 8) /  1.00000000d-05 ,  1.72000000d-05 /
      data f( 8, 9),f( 8,10) /  1.33000000d-05 , -1.31000000d-05 /
      data f( 8,11),f( 8,12) / -5.59000000d-05 , -1.37100000d-04 /
      data f( 8,13),f( 8,14) / -2.29800000d-04 , -3.23800000d-04 /
      data f( 8,15),f( 8,16) / -3.47100000d-04 , -2.46700000d-04 /
      data f( 8,17),f( 8,18) / -8.82000000d-05 ,  1.82000000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 9, 3),f( 9, 4) / -1.00000000d-06 , -2.80000000d-06 /
      data f( 9, 5),f( 9, 6) / -2.30000000d-06 ,  5.60000000d-06 /
      data f( 9, 7),f( 9, 8) /  1.02000000d-05 ,  1.40000000d-05 /
      data f( 9, 9),f( 9,10) /  1.08000000d-05 , -6.30000000d-06 /
      data f( 9,11),f( 9,12) / -2.77000000d-05 , -6.39000000d-05 /
      data f( 9,13),f( 9,14) / -1.09200000d-04 , -1.47000000d-04 /
      data f( 9,15),f( 9,16) / -1.48400000d-04 , -9.94000000d-05 /
      data f( 9,17),f( 9,18) / -3.38000000d-05 , -2.80000000d-06 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(10, 3),f(10, 4) /  2.00000000d-07 ,  1.00000000d-06 /
      data f(10, 5),f(10, 6) /  2.10000000d-06 ,  4.10000000d-06 /
      data f(10, 7),f(10, 8) /  6.20000000d-06 ,  6.80000000d-06 /
      data f(10, 9),f(10,10) /  5.10000000d-06 , -4.00000000d-07 /
      data f(10,11),f(10,12) / -8.20000000d-06 , -1.81000000d-05 /
      data f(10,13),f(10,14) / -2.73000000d-05 , -3.20000000d-05 /
      data f(10,15),f(10,16) / -2.83000000d-05 , -1.66000000d-05 /
      data f(10,17),f(10,18) / -5.10000000d-06 , -5.00000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(11, 3),f(11, 4) /  2.00000000d-07 ,  7.00000000d-07 /
      data f(11, 5),f(11, 6) /  1.50000000d-06 ,  2.20000000d-06 /
      data f(11, 7),f(11, 8) /  2.50000000d-06 ,  2.10000000d-06 /
      data f(11, 9),f(11,10) /  1.00000000d-06 , -6.00000000d-07 /
      data f(11,11),f(11,12) / -3.00000000d-06 , -5.50000000d-06 /
      data f(11,13),f(11,14) / -8.60000000d-06 , -8.10000000d-06 /
      data f(11,15),f(11,16) / -6.30000000d-06 , -3.30000000d-06 /
      data f(11,17),f(11,18) / -9.00000000d-07 ,  0.00000000d+00 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) /  1.00000000d-07 ,  3.00000000d-07 /
      data f(12, 5),f(12, 6) /  6.00000000d-07 ,  7.00000000d-07 /
      data f(12, 7),f(12, 8) /  6.00000000d-07 ,  4.00000000d-07 /
      data f(12, 9),f(12,10) /  0.00000000d+00 , -6.00000000d-07 /
      data f(12,11),f(12,12) / -1.20000000d-06 , -1.80000000d-06 /
      data f(12,13),f(12,14) / -2.60000000d-06 , -2.20000000d-06 /
      data f(12,15),f(12,16) / -1.60000000d-06 , -7.00000000d-07 /
      data f(12,17),f(12,18) / -2.00000000d-07 ,  0.00000000d+00 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(13, 5),f(13, 6) /  3.00000000d-07 ,  3.00000000d-07 /
      data f(13, 7),f(13, 8) /  3.00000000d-07 ,  1.00000000d-07 /
      data f(13, 9),f(13,10) /  0.00000000d+00 , -2.00000000d-07 /
      data f(13,11),f(13,12) / -3.00000000d-07 , -4.00000000d-07 /
      data f(13,13),f(13,14) / -5.00000000d-07 , -5.00000000d-07 /
      data f(13,15),f(13,16) / -3.00000000d-07 , -2.00000000d-07 /
      data f(13,17),f(13,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-1.62202124d-06/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-2.66844781d-05,-1.85895752d-06/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-3.96245303d-04,-3.72814868d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-1.83998831d-03, 1.15552257d-07/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-4.48537205d-03, 6.49993966d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-4.79611409d-03, 4.22668911d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-2.75097022d-03,-5.16696111d-07/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-9.52856178d-04,-1.83590467d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 8.00239299d-05,-5.83768521d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/-7.19255064d-04,-8.00535448d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/-3.45066751d-04,-9.72689687d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/-1.69911852d-03,-8.91105803d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/-4.47732006d-03,-3.12728710d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/ 8.10903989d-03, 4.45065420d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.51049835d-03, 3.85027028d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/-3.36154541d-02, 8.60506467d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/-1.47209533d-02,-7.88052895d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-7.27061646d-02,-3.65694887d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-5.35027556d-05/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-1.36965847d-06/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-2.25596153d-05,-1.35468306d-06/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-3.45080822d-04,-2.70360930d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-1.46802338d-03, 7.21120243d-07/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-3.35454161d-03, 4.21112832d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-3.54005753d-03, 2.53436646d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-2.04691670d-03,-6.25941556d-08/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-7.19144787d-04,-1.43798984d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 8.38071162d-06,-4.04944650d-06/
      data fpp( 2,10,1),fpp( 2,10,2)/-4.88275586d-04,-7.22222415d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/-4.99223641d-04,-8.23365690d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/-1.66012010d-03,-5.70714827d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/-4.38235988d-03,-2.58957500d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.60113450d-03, 2.57461484d-05/
      data fpp( 2,15,1),fpp( 2,15,2)/ 3.81436044d-03, 5.21871564d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/-1.88278061d-02, 4.64672261d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/-5.66852194d-03,-2.73300609d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-4.12395994d-02,-5.60609826d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-9.76680087d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.04285634d-06/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-4.80770608d-05,-9.77287314d-07/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.53431410d-04,-1.96599440d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.04791819d-03, 7.35264914d-07/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.80646152d-03, 2.84693474d-06/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-1.86365579d-03, 1.67099611d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-1.12136296d-03, 9.90807990d-08/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-4.30564675d-04,-1.08331931d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 9.64532236d-05,-3.12780356d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/-1.26425907d-05,-5.87546646d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/-5.23038685d-04,-7.04233060d-06/
      data fpp( 3,12,1),fpp( 3,12,2)/-1.45540107d-03,-4.55921115d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/-4.72324042d-03,-1.59948248d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/-9.81857790d-03, 1.02905104d-05/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.17059899d-04, 4.29147831d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.77566785d-02, 4.23483571d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/ 1.79250411d-02,-2.03742115d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 2.12595623d-02,-4.03935110d-05/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-7.26137445d-05/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-6.21454750d-07/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 1.86299463d-05,-5.84090499d-07/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.55174754d-04,-1.22418325d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-5.88257132d-04, 5.06823513d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-9.62100522d-04, 1.56688920d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-8.37775682d-04, 1.13361967d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-5.17512315d-04, 1.26632099d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-1.92021225d-04,-6.98148071d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-3.37097886d-04,-2.21803982d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.65674307d-04,-4.17569266d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.02038862d-03,-5.33718953d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/-1.73524970d-03,-4.01954920d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/-3.69095867d-03,-5.39261365d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/-9.55216332d-03, 1.51200381d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/-5.63310662d-03, 2.78685984d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/-4.56372420d-03, 3.07496025d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.81554433d-03,-1.26750084d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 5.10119186d-03,-2.30575690d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-4.28807155d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-2.48951974d-07/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.37169592d-04,-2.50096053d-07/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-3.65781942d-04,-5.26663815d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-5.08026266d-04, 1.60751314d-07/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-4.79001418d-04, 6.39658560d-07/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-3.81724345d-04, 5.86614447d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-3.44742812d-04, 1.87883654d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-2.58393276d-04,-3.78149062d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 7.40041302d-04,-1.15928741d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/-2.22256664d-04,-2.51470132d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 1.58266800d-04,-2.80390733d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/-3.93960320d-04,-2.55966935d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/-1.00858801d-03,-3.24141525d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/-2.77873062d-03, 1.15233304d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/-1.27298947d-02, 3.78093798d-07/
      data fpp( 5,16,1),fpp( 5,16,2)/-9.23648964d-03, 1.82002944d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/-4.53756652d-03,-4.97127158d-06/
      data fpp( 5,18,1),fpp( 5,18,2)/-3.75075136d-03,-1.02132081d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-2.05778959d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-2.63718511d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 1.96448420d-04,-1.33562977d-07/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 3.79902523d-04,-2.95800978d-11/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 2.39562196d-04, 1.21681298d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-3.06938079d-05, 1.91304390d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.26526938d-04, 1.81101142d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 5.80835627d-05,-2.07708959d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 2.22394329d-04, 2.08373469d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-1.29106732d-03,-4.46122981d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-2.29299037d-04, 2.06918456d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-1.61187858d-03,-2.84950844d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/-2.48490902d-03,-1.41115079d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.39148930d-03,-1.07588839d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/-3.63051419d-03, 2.13870437d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 1.11508561d-03, 5.55907093d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.96082753d-04, 7.52101192d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.26032175d-03,-7.13111859d-06/
      data fpp( 6,18,1),fpp( 6,18,2)/ 5.40181358d-03, 1.05346246d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 2.71326877d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-1.95804444d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.27824090d-04,-3.78391112d-08/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-2.92228150d-04,-9.90631108d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-2.94222519d-04,-9.90844552d-09/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-1.87823351d-04, 1.38696893d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.22567902d-04, 9.71208740d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-6.03914390d-05, 3.60819611d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-2.63984041d-04,-3.76399319d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 7.01827991d-04,-6.13222335d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/ 5.22528133d-05,-2.12711341d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 5.50847512d-04,-7.97932300d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 3.92796407d-04,-6.63559459d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/ 4.73452239d-05,-3.63829863d-07/
      data fpp( 7,14,1),fpp( 7,14,2)/-7.61612614d-04, 1.59687891d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/-2.98164770d-03, 2.68831421d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/-2.85664138d-03, 3.03386423d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/-2.57732048d-03,-2.77577114d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/-3.36530297d-03,-1.60877967d-06/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-2.43511016d-06/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-1.45199088d-08/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 3.16479384d-05,-1.69601824d-08/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 6.66100757d-05,-4.96393615d-08/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 4.93278779d-05, 2.35176284d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 1.39872114d-05,-8.43115203d-09/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 3.11985455d-05, 3.82206980d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.62117807d-04,-3.08396767d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 2.71418352d-05,-6.61991204d-09/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.77044639d-04,-3.31123585d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/-7.81122158d-05,-1.88857483d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-3.25111468d-04,-5.77333422d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-4.18276608d-04, 2.42194352d-08/
      data fpp( 8,13,1),fpp( 8,13,2)/-7.79491591d-04,-2.09544319d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/-9.98235352d-04, 7.35957841d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.14894794d-04, 1.50771295d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/-5.27117249d-04, 6.55190342d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/ 1.20960159d-04,-6.42474323d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/ 4.72998315d-04,-1.21129305d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.98835347d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-1.23898016d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-5.96766399d-06,-9.22039676d-09/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.50121530d-05,-1.07286114d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-2.06889931d-05, 4.13484217d-09/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-3.61254949d-05, 1.32189243d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-8.38262802d-05,-8.88918129d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.00626657d-05, 2.53780090d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-7.78329974d-06,-6.06202231d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 3.03505654d-05,-2.02897117d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/ 9.96049792d-07, 3.82086899d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/-1.36016385d-05,-2.07937643d-07/
      data fpp( 9,12,1),fpp( 9,12,2)/-1.26089976d-04,-9.44581179d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/-1.71778858d-04, 3.97701145d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-2.90245979d-04, 3.85377660d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/-4.41173118d-04, 6.02719246d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/-3.29289627d-04, 2.27745357d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.84120158d-04,-5.17700673d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.17097156d-05,-2.32942665d-07/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-2.42528668d-07/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00,-6.40032340d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.47902276d-06,-1.19935320d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.13142123d-06, 1.11977362d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/-9.96959587d-07,-1.59159158d-09/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.78287891d-06, 1.31686301d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/ 2.16795680d-05, 2.91707111d-09/
      data fpp(10, 7,1),fpp(10, 7,2)/-5.52909363d-06,-1.88369145d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 4.97898162d-06,-1.75694129d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-6.72937665d-06,-4.88854338d-08/
      data fpp(10,10,1),fpp(10,10,2)/-1.01320415d-05,-1.48888519d-08/
      data fpp(10,11,1),fpp(10,11,2)/-1.80393502d-05,-2.95591586d-08/
      data fpp(10,12,1),fpp(10,12,2)/-1.61917669d-05, 7.12548641d-09/
      data fpp(10,13,1),fpp(10,13,2)/-5.07176295d-05, 4.30572130d-08/
      data fpp(10,14,1),fpp(10,14,2)/-6.17443858d-05, 9.06456617d-08/
      data fpp(10,15,1),fpp(10,15,2)/-3.28332502d-05, 9.83601404d-08/
      data fpp(10,16,1),fpp(10,16,2)/-1.93724942d-05,-4.08622320d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 1.12803944d-05,-9.40152476d-08/
      data fpp(10,18,1),fpp(10,18,2)/-5.82830417d-06,-3.38527864d-08/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-1.65736068d-08/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-3.82941359d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/-5.48427041d-07,-3.41172815d-10/
      data fpp(11, 3,1),fpp(11, 3,2)/-7.13531903d-07, 5.19410485d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/ 7.68314703d-08, 3.56475340d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-1.00602077d-06,-1.45311846d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/-5.29199157d-06,-3.75227957d-09/
      data fpp(11, 7,1),fpp(11, 7,2)/ 3.85370885d-06,-7.53776325d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 2.86737326d-06,-8.09666745d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/ 6.16694122d-06,-2.07556697d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 2.93211618d-06,-1.36010647d-08/
      data fpp(11,11,1),fpp(11,11,2)/-4.09604608d-08, 8.47982562d-09/
      data fpp(11,12,1),fpp(11,12,2)/-8.34295609d-06,-2.63182378d-08/
      data fpp(11,13,1),fpp(11,13,2)/-4.55062382d-06, 6.07931256d-08/
      data fpp(11,14,1),fpp(11,14,2)/-9.37647743d-06,-8.54264594d-10/
      data fpp(11,15,1),fpp(11,15,2)/-1.60938815d-05, 2.06239328d-08/
      data fpp(11,16,1),fpp(11,16,2)/-1.02203963d-05,-9.64146652d-09/
      data fpp(11,17,1),fpp(11,17,2)/-8.00141957d-06,-1.80580667d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 8.03501098d-07,-8.12626666d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-3.43686667d-09/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 9.59002930d-10/
      data fpp(12, 2,1),fpp(12, 2,2)/ 1.14685408d-07, 1.08199414d-09/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.22706381d-07, 7.13020510d-10/
      data fpp(12, 4,1),fpp(12, 4,2)/ 8.96337059d-08, 2.06592382d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/ 4.41204153d-07,-2.97671578d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/ 1.88839831d-06,-2.15906068d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/ 9.14258230d-07,-3.87041487d-10/
      data fpp(12, 8,1),fpp(12, 8,2)/ 1.55152535d-06,-2.29277337d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 6.61611755d-07,-2.44186503d-09/
      data fpp(12,10,1),fpp(12,10,2)/-3.96423237d-07, 6.02334854d-11/
      data fpp(12,11,1),fpp(12,11,2)/-2.19680791d-06, 2.20093109d-09/
      data fpp(12,12,1),fpp(12,12,2)/-3.83640878d-06,-8.86395783d-09/
      data fpp(12,13,1),fpp(12,13,2)/-7.27987524d-06, 2.12549002d-08/
      data fpp(12,14,1),fpp(12,14,2)/-8.74970451d-06,-4.15564314d-09/
      data fpp(12,15,1),fpp(12,15,2)/-6.59122370d-06, 7.36767233d-09/
      data fpp(12,16,1),fpp(12,16,2)/-3.94592075d-06,-7.31504619d-09/
      data fpp(12,17,1),fpp(12,17,2)/-2.74716086d-07,-2.10748756d-09/
      data fpp(12,18,1),fpp(12,18,2)/-3.85700220d-07,-2.25500355d-09/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-8.72498223d-10/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-1.01136873d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/ 2.30157296d-07, 2.27374565d-11/
      data fpp(13, 3,1),fpp(13, 3,2)/ 1.38646810d-07, 9.20418902d-10/
      data fpp(13, 4,1),fpp(13, 4,2)/ 5.92683147d-07, 2.29558693d-09/
      data fpp(13, 5,1),fpp(13, 5,2)/ 1.42939792d-06,-4.10276664d-09/
      data fpp(13, 6,1),fpp(13, 6,2)/ 8.80800843d-07, 2.11547963d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 5.80370885d-07,-4.35915186d-09/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.43826267d-06, 3.32112783d-09/
      data fpp(13, 9,1),fpp(13, 9,2)/-2.06830588d-06,-2.92535944d-09/
      data fpp(13,10,1),fpp(13,10,2)/ 3.23211618d-07, 2.38030994d-09/
      data fpp(13,11,1),fpp(13,11,2)/ 2.56090395d-06,-5.95880332d-10/
      data fpp(13,12,1),fpp(13,12,2)/ 6.68070439d-06, 3.21138519d-12/
      data fpp(13,13,1),fpp(13,13,2)/ 9.26493762d-06, 5.83034791d-10/
      data fpp(13,14,1),fpp(13,14,2)/ 1.57873523d-05, 3.66464945d-09/
      data fpp(13,15,1),fpp(13,15,2)/ 1.56706118d-05,-3.24163259d-09/
      data fpp(13,16,1),fpp(13,16,2)/ 9.89796037d-06, 3.30188090d-09/
      data fpp(13,17,1),fpp(13,17,2)/ 3.02485804d-06,-3.96589101d-09/
      data fpp(13,18,1),fpp(13,18,2)/ 7.55350110d-07, 5.61683146d-10/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 1.71915843d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 4.11900654d-06,-3.89818677d-06/
      data fpppp( 1, 3),fpppp( 1, 4)/-9.09884028d-06,-2.41573831d-05/
      data fpppp( 1, 5),fpppp( 1, 6)/ 3.36299287d-05, 2.97161705d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/-1.11414559d-05, 2.78633583d-08/
      data fpppp( 1, 9),fpppp( 1,10)/-3.48840335d-05, 2.95787245d-05/
      data fpppp( 1,11),fpppp( 1,12)/-1.30228262d-05,-8.11818246d-05/
      data fpppp( 1,13),fpppp( 1,14)/ 2.52301139d-04,-6.14903987d-06/
      data fpppp( 1,15),fpppp( 1,16)/-1.01879907d-03, 1.64970066d-03/
      data fpppp( 1,17),fpppp( 1,18)/-1.97877637d-03, 1.65262212d-03/
      data fpppp( 1,19) /             3.20977046d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 2.10780499d-06,-3.20526260d-06/
      data fpppp( 2, 3),fpppp( 2, 4)/-7.28445008d-06,-1.56822179d-05/
      data fpppp( 2, 5),fpppp( 2, 6)/ 2.41987812d-05, 2.09472313d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/-7.26830166d-06,-1.79615921d-06/
      data fpppp( 2, 9),fpppp( 2,10)/-2.15618467d-05, 1.45926381d-05/
      data fpppp( 2,11),fpppp( 2,12)/-7.66621108d-06,-5.29246982d-05/
      data fpppp( 2,13),fpppp( 2,14)/ 1.25684405d-04, 7.25311281d-05/
      data fpppp( 2,15),fpppp( 2,16)/-6.42025024d-04, 1.00424542d-03/
      data fpppp( 2,17),fpppp( 2,18)/-1.22686961d-03, 9.79411340d-04/
      data fpppp( 2,19) /             1.91786487d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.38432277d-06,-2.31257031d-06/
      data fpppp( 3, 3),fpppp( 3, 4)/-4.80203325d-06,-1.82724233d-06/
      data fpppp( 3, 5),fpppp( 3, 6)/ 8.26760886d-06, 1.08377511d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-3.64938758d-06, 6.70127109d-07/
      data fpppp( 3, 9),fpppp( 3,10)/-8.85794428d-06,-3.40517280d-06/
      data fpppp( 3,11),fpppp( 3,12)/-1.59938136d-06,-1.55152791d-05/
      data fpppp( 3,13),fpppp( 3,14)/-7.64681203d-05, 2.11737873d-04/
      data fpppp( 3,15),fpppp( 3,16)/ 1.43375147d-04,-3.46999613d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 2.08347946d-04,-2.96422653d-04/
      data fpppp( 3,19) /            -4.98302349d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-1.37756741d-06,-1.68104393d-06/
      data fpppp( 4, 3),fpppp( 4, 4)/-3.44433564d-06,-9.82742210d-08/
      data fpppp( 4, 5),fpppp( 4, 6)/ 7.39177186d-06, 4.21280534d-07/
      data fpppp( 4, 7),fpppp( 4, 8)/ 2.67941758d-06,-1.08252874d-05/
      data fpppp( 4, 9),fpppp( 4,10)/ 1.23876667d-05,-1.97353651d-05/
      data fpppp( 4,11),fpppp( 4,12)/ 4.98551984d-06, 8.18447991d-06/
      data fpppp( 4,13),fpppp( 4,14)/-1.12174313d-04, 2.06183031d-04/
      data fpppp( 4,15),fpppp( 4,16)/-1.25742128d-04, 1.25805024d-04/
      data fpppp( 4,17),fpppp( 4,18)/-5.88848005d-05,-7.58830824d-05/
      data fpppp( 4,19) /            -1.40793234d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.67718768d-06,-9.45123123d-07/
      data fpppp( 5, 3),fpppp( 5, 4)/ 9.71114634d-07, 2.24274621d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/ 3.34050865d-07, 5.16183790d-07/
      data fpppp( 5, 7),fpppp( 5, 8)/-6.01651839d-06, 2.65119699d-05/
      data fpppp( 5, 9),fpppp( 5,10)/-4.53062587d-05, 3.70691122d-05/
      data fpppp( 5,11),fpppp( 5,12)/-2.24009044d-05,-3.43052947d-06/
      data fpppp( 5,13),fpppp( 5,14)/ 3.23789882d-05,-1.95416319d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 2.58424998d-04,-3.16095168d-05/
      data fpppp( 5,17),fpppp( 5,18)/-5.96558500d-05, 3.55064393d-05/
      data fpppp( 5,19) /             9.54662650d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 2.79185459d-06, 2.42463530d-07/
      data fpppp( 6, 3),fpppp( 6, 4)/-4.54136777d-06,-1.50465823d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 2.76506004d-06, 9.09790470d-07/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.04223960d-05,-4.38173584d-05/
      data fpppp( 6, 9),fpppp( 6,10)/ 6.41806925d-05,-5.83916151d-05/
      data fpppp( 6,11),fpppp( 6,12)/ 2.27248985d-05,-1.93503294d-06/
      data fpppp( 6,13),fpppp( 6,14)/-1.69977570d-05, 1.09979385d-04/
      data fpppp( 6,15),fpppp( 6,16)/-1.23842300d-04, 4.55136555d-05/
      data fpppp( 6,17),fpppp( 6,18)/ 6.07821895d-05,-1.04007243d-04/
      data fpppp( 6,19) /            -2.17351542d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/-2.27652015d-06,-5.24152031d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/ 2.17833006d-06, 1.55541325d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/-1.89637089d-06, 3.56144720d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.25341571d-05, 3.06290372d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-3.98177138d-05, 3.17185854d-05/
      data fpppp( 7,11),fpppp( 7,12)/-1.81664353d-05, 1.54840768d-06/
      data fpppp( 7,13),fpppp( 7,14)/ 7.28799820d-07,-3.22740062d-05/
      data fpppp( 7,15),fpppp( 7,16)/ 4.37025900d-05,-1.83386876d-06/
      data fpppp( 7,17),fpppp( 7,18)/-2.71082408d-05, 4.62286281d-05/
      data fpppp( 7,19) /             9.13908567d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 5.76365269d-07, 5.78661039d-08/
      data fpppp( 8, 3),fpppp( 8, 4)/-6.08977755d-07,-7.56615186d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/ 2.55193038d-06,-6.29798629d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.00083536d-05,-1.07808684d-05/
      data fpppp( 8, 9),fpppp( 8,10)/ 9.50835321d-06,-9.06541054d-06/
      data fpppp( 8,11),fpppp( 8,12)/ 5.99738840d-06,-5.69409627d-06/
      data fpppp( 8,13),fpppp( 8,14)/ 6.96006030d-07, 1.14583456d-05/
      data fpppp( 8,15),fpppp( 8,16)/-1.04043292d-05, 1.24251907d-05/
      data fpppp( 8,17),fpppp( 8,18)/-5.67844163d-06,-7.47377925d-06/
      data fpppp( 8,19) /            -1.39286296d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-9.92712063d-08,-2.26518184d-08/
      data fpppp( 9, 3),fpppp( 9, 4)/ 5.26897706d-09, 2.03634846d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.40538806d-06, 3.48206038d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/-3.42746959d-06, 2.32372330d-06/
      data fpppp( 9, 9),fpppp( 9,10)/-1.90863379d-06, 1.26150902d-06/
      data fpppp( 9,11),fpppp( 9,12)/-2.25199266d-06, 1.87302265d-06/
      data fpppp( 9,13),fpppp( 9,14)/-1.23213059d-06,-1.31119464d-06/
      data fpppp( 9,15),fpppp( 9,16)/ 4.52930813d-06,-1.03740018d-06/
      data fpppp( 9,17),fpppp( 9,18)/ 1.61745132d-06,-2.39278083d-06/
      data fpppp( 9,19) /            -4.49870336d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 2.00065676d-08,-5.72876364d-09/
      data fpppp(10, 3),fpppp(10, 4)/-4.66889701d-08,-3.43621135d-08/
      data fpppp(10, 5),fpppp(10, 6)/ 5.38630583d-07,-1.09314919d-06/
      data fpppp(10, 7),fpppp(10, 8)/ 1.00764512d-06,-6.74427097d-07/
      data fpppp(10, 9),fpppp(10,10)/ 3.57077253d-07,-2.55540310d-07/
      data fpppp(10,11),fpppp(10,12)/ 3.94805351d-07,-7.38387567d-07/
      data fpppp(10,13),fpppp(10,14)/ 3.76338157d-07, 6.42981316d-07/
      data fpppp(10,15),fpppp(10,16)/-5.51989911d-07, 6.37955560d-07/
      data fpppp(10,17),fpppp(10,18)/-9.68304382d-07, 3.69566745d-07/
      data fpppp(10,19) /             8.66257566d-07 /
      data fpppp(11, 1),fpppp(11, 2)/-1.26647479d-09, 2.58969389d-09/
      data fpppp(11, 3),fpppp(11, 4)/ 1.39070300d-08,-8.89719981d-10/
      data fpppp(11, 5),fpppp(11, 6)/-1.22741087d-07, 2.99666953d-07/
      data fpppp(11, 7),fpppp(11, 8)/-2.70026452d-07, 1.72516696d-07/
      data fpppp(11, 9),fpppp(11,10)/-1.62886117d-07, 8.69641909d-08/
      data fpppp(11,11),fpppp(11,12)/-1.69265743d-07, 2.70363641d-07/
      data fpppp(11,13),fpppp(11,14)/-1.86529150d-07,-4.13381943d-08/
      data fpppp(11,15),fpppp(11,16)/ 2.38388898d-07,-1.56764039d-07/
      data fpppp(11,17),fpppp(11,18)/ 1.69396746d-07,-1.25666305d-07/
      data fpppp(11,19) /            -2.43236831d-07 /
      data fpppp(12, 1),fpppp(12, 2)/-1.63562820d-09,-1.24009247d-09/
      data fpppp(12, 3),fpppp(12, 4)/ 1.96131909d-10,-2.01005399d-09/
      data fpppp(12, 5),fpppp(12, 6)/ 3.09226713d-08,-5.59432086d-08/
      data fpppp(12, 7),fpppp(12, 8)/ 4.75701087d-08,-3.76527940d-08/
      data fpppp(12, 9),fpppp(12,10)/ 1.14102246d-08,-1.80753882d-08/
      data fpppp(12,11),fpppp(12,12)/ 1.63503475d-08,-3.76789741d-08/
      data fpppp(12,13),fpppp(12,14)/ 2.61336142d-08, 5.15627477d-08/
      data fpppp(12,15),fpppp(12,16)/-1.46859994d-08, 3.63905781d-08/
      data fpppp(12,17),fpppp(12,18)/-6.93222103d-08, 1.39669356d-08/
      data fpppp(12,19) /             4.32555289d-08 /
      data fpppp(13, 1),fpppp(13, 2)/-1.17806851d-08,-3.43295543d-09/
      data fpppp(13, 3),fpppp(13, 4)/ 6.21243987d-09, 1.13160054d-08/
      data fpppp(13, 5),fpppp(13, 6)/-2.85157550d-08, 1.96283033d-08/
      data fpppp(13, 7),fpppp(13, 8)/-3.51074308d-08, 1.77092037d-08/
      data fpppp(13, 9),fpppp(13,10)/ 4.75860372d-08,-2.67597105d-08/
      data fpppp(13,11),fpppp(13,12)/ 5.02232953d-08,-6.12069845d-08/
      data fpppp(13,13),fpppp(13,14)/ 1.02470610d-07,-1.12384572d-07/
      data fpppp(13,15),fpppp(13,16)/-5.12816269d-08,-2.18435850d-08/
      data fpppp(13,17),fpppp(13,18)/ 7.26289156d-08, 7.54358633d-09/
      data fpppp(13,19) /            -1.19537916d-08 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
      subroutine s5_spl_ch2oh_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(13,19,2),f(13,19),fpppp(13,19)
      dimension delx(12),dely(18),x(13),y(19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -1.20000000d-06 /
      data f( 1, 3),f( 1, 4) / -2.09000000d-05 , -1.14200000d-04 /
      data f( 1, 5),f( 1, 6) / -2.67300000d-04 , -3.14700000d-04 /
      data f( 1, 7),f( 1, 8) / -2.17100000d-04 , -6.34000000d-05 /
      data f( 1, 9),f( 1,10) /  2.29000000d-05 ,  5.02000000d-05 /
      data f( 1,11),f( 1,12) /  5.86000000d-05 ,  5.67000000d-05 /
      data f( 1,13),f( 1,14) / -9.43000000d-05 , -9.61900000d-04 /
      data f( 1,15),f( 1,16) / -1.08330000d-03 , -4.70900000d-04 /
      data f( 1,17),f( 1,18) / -5.80400000d-04 ,  1.19500000d-04 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  3.00000000d-07 /
      data f( 2, 3),f( 2, 4) / -1.52000000d-05 , -8.25000000d-05 /
      data f( 2, 5),f( 2, 6) / -1.85100000d-04 , -2.16200000d-04 /
      data f( 2, 7),f( 2, 8) / -1.45500000d-04 , -3.86000000d-05 /
      data f( 2, 9),f( 2,10) /  2.32000000d-05 ,  5.24000000d-05 /
      data f( 2,11),f( 2,12) /  2.38000000d-05 ,  4.00000000d-06 /
      data f( 2,13),f( 2,14) / -6.29000000d-05 , -6.40100000d-04 /
      data f( 2,15),f( 2,16) / -1.07920000d-03 , -7.84000000d-04 /
      data f( 2,17),f( 2,18) / -7.03800000d-04 ,  8.91000000d-05 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -6.00000000d-07 /
      data f( 3, 3),f( 3, 4) / -1.17000000d-05 , -5.91000000d-05 /
      data f( 3, 5),f( 3, 6) / -1.27900000d-04 , -1.46700000d-04 /
      data f( 3, 7),f( 3, 8) / -9.63000000d-05 , -2.37000000d-05 /
      data f( 3, 9),f( 3,10) /  2.19000000d-05 ,  4.22000000d-05 /
      data f( 3,11),f( 3,12) /  4.50000000d-06 , -3.18000000d-05 /
      data f( 3,13),f( 3,14) / -5.86000000d-05 , -3.78300000d-04 /
      data f( 3,15),f( 3,16) / -7.89200000d-04 , -7.59700000d-04 /
      data f( 3,17),f( 3,18) / -5.67800000d-04 ,  1.94000000d-05 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -4.00000000d-07 /
      data f( 4, 3),f( 4, 4) / -7.10000000d-06 , -3.46000000d-05 /
      data f( 4, 5),f( 4, 6) / -7.24000000d-05 , -8.37000000d-05 /
      data f( 4, 7),f( 4, 8) / -5.33000000d-05 , -1.25000000d-05 /
      data f( 4, 9),f( 4,10) /  1.97000000d-05 ,  3.27000000d-05 /
      data f( 4,11),f( 4,12) / -6.30000000d-06 , -5.74000000d-05 /
      data f( 4,13),f( 4,14) / -7.54000000d-05 , -1.67400000d-04 /
      data f( 4,15),f( 4,16) / -3.92400000d-04 , -4.70100000d-04 /
      data f( 4,17),f( 4,18) / -3.23500000d-04 , -2.06000000d-05 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -2.00000000d-07 /
      data f( 5, 3),f( 5, 4) / -2.90000000d-06 , -1.36000000d-05 /
      data f( 5, 5),f( 5, 6) / -2.85000000d-05 , -3.48000000d-05 /
      data f( 5, 7),f( 5, 8) / -2.53000000d-05 , -6.90000000d-06 /
      data f( 5, 9),f( 5,10) /  1.80000000d-05 ,  2.46000000d-05 /
      data f( 5,11),f( 5,12) / -4.40000000d-06 , -2.66000000d-05 /
      data f( 5,13),f( 5,14) / -8.88000000d-05 , -1.01900000d-04 /
      data f( 5,15),f( 5,16) / -5.94100000d-04 , -3.94200000d-04 /
      data f( 5,17),f( 5,18) / -1.25400000d-04 , -8.00000000d-06 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.40000000d-06 /
      data f( 6, 3),f( 6, 4) /  1.38000000d-05 ,  2.41000000d-05 /
      data f( 6, 5),f( 6, 6) /  5.90000000d-06 , -1.35000000d-05 /
      data f( 6, 7),f( 6, 8) / -9.90000000d-06 , -5.40000000d-06 /
      data f( 6, 9),f( 6,10) /  9.13000000d-05 ,  2.73000000d-05 /
      data f( 6,11),f( 6,12) /  1.30000000d-06 , -1.90000000d-06 /
      data f( 6,13),f( 6,14) / -3.17000000d-05 , -3.74000000d-05 /
      data f( 6,15),f( 6,16) / -3.57000000d-05 , -3.34000000d-05 /
      data f( 6,17),f( 6,18) /  3.61000000d-05 ,  1.08600000d-04 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 7, 3),f( 7, 4) / -2.00000000d-07 , -2.00000000d-06 /
      data f( 7, 5),f( 7, 6) / -5.70000000d-06 , -9.30000000d-06 /
      data f( 7, 7),f( 7, 8) / -1.15000000d-05 , -5.20000000d-06 /
      data f( 7, 9),f( 7,10) /  1.94000000d-05 ,  9.10000000d-06 /
      data f( 7,11),f( 7,12) /  7.00000000d-07 , -1.22000000d-05 /
      data f( 7,13),f( 7,14) / -2.89000000d-05 , -4.34000000d-05 /
      data f( 7,15),f( 7,16) / -1.62000000d-05 , -2.91000000d-05 /
      data f( 7,17),f( 7,18) / -1.24000000d-05 , -7.00000000d-07 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 8, 3),f( 8, 4) / -1.00000000d-07 , -9.00000000d-07 /
      data f( 8, 5),f( 8, 6) / -3.00000000d-06 , -5.60000000d-06 /
      data f( 8, 7),f( 8, 8) / -1.02000000d-05 ,  1.50000000d-06 /
      data f( 8, 9),f( 8,10) /  1.14000000d-05 ,  6.40000000d-06 /
      data f( 8,11),f( 8,12) /  6.00000000d-07 , -4.90000000d-06 /
      data f( 8,13),f( 8,14) / -1.29000000d-05 , -1.84000000d-05 /
      data f( 8,15),f( 8,16) / -1.86000000d-05 , -1.21000000d-05 /
      data f( 8,17),f( 8,18) / -3.70000000d-06 ,  2.26000000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f( 9, 3),f( 9, 4) / -1.00000000d-07 , -6.00000000d-07 /
      data f( 9, 5),f( 9, 6) / -2.30000000d-06 , -4.60000000d-06 /
      data f( 9, 7),f( 9, 8) / -1.30000000d-06 ,  3.10000000d-06 /
      data f( 9, 9),f( 9,10) /  6.10000000d-06 ,  6.10000000d-06 /
      data f( 9,11),f( 9,12) /  5.00000000d-07 , -3.30000000d-06 /
      data f( 9,13),f( 9,14) / -6.50000000d-06 , -8.00000000d-06 /
      data f( 9,15),f( 9,16) / -6.70000000d-06 , -3.80000000d-06 /
      data f( 9,17),f( 9,18) / -1.00000000d-06 ,  0.00000000d+00 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(10, 3),f(10, 4) / -1.00000000d-07 , -2.00000000d-07 /
      data f(10, 5),f(10, 6) / -1.00000000d-07 ,  5.00000000d-07 /
      data f(10, 7),f(10, 8) /  1.40000000d-06 ,  2.10000000d-06 /
      data f(10, 9),f(10,10) /  2.40000000d-06 ,  3.10000000d-06 /
      data f(10,11),f(10,12) /  5.00000000d-07 , -6.00000000d-07 /
      data f(10,13),f(10,14) / -1.70000000d-06 , -1.90000000d-06 /
      data f(10,15),f(10,16) / -1.50000000d-06 , -6.00000000d-07 /
      data f(10,17),f(10,18) / -1.00000000d-07 ,  0.00000000d+00 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(11, 3),f(11, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(11, 5),f(11, 6) /  2.00000000d-07 ,  5.00000000d-07 /
      data f(11, 7),f(11, 8) /  8.00000000d-07 ,  9.00000000d-07 /
      data f(11, 9),f(11,10) /  9.00000000d-07 ,  7.00000000d-07 /
      data f(11,11),f(11,12) /  6.00000000d-07 ,  2.00000000d-07 /
      data f(11,13),f(11,14) /  5.00000000d-07 , -1.00000000d-07 /
      data f(11,15),f(11,16) / -3.00000000d-07 , -2.00000000d-07 /
      data f(11,17),f(11,18) /  0.00000000d+00 ,  2.00000000d-07 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 3),f(12, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(12, 5),f(12, 6) /  1.00000000d-07 ,  1.00000000d-07 /
      data f(12, 7),f(12, 8) /  2.00000000d-07 ,  1.00000000d-07 /
      data f(12, 9),f(12,10) /  1.00000000d-07 ,  2.00000000d-07 /
      data f(12,11),f(12,12) /  2.00000000d-07 ,  0.00000000d+00 /
      data f(12,13),f(12,14) /  2.00000000d-07 ,  1.00000000d-07 /
      data f(12,15),f(12,16) /  0.00000000d+00 , -1.00000000d-07 /
      data f(12,17),f(12,18) /  0.00000000d+00 , -1.00000000d-07 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 3),f(13, 4) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 5),f(13, 6) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 7),f(13, 8) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13, 9),f(13,10) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,11),f(13,12) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,13),f(13,14) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,15),f(13,16) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,17),f(13,18) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(13,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 3.57164367d-07/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-1.23266005d-04,-1.67328735d-07/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-8.90155249d-05,-7.97849429d-07/
      data fpp( 1, 4,1),fpp( 1, 4,2)/-2.57203320d-04,-1.05727355d-06/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-8.06368411d-04, 1.43894363d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-8.77863772d-04, 1.64349904d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-6.87442313d-04, 6.87060230d-07/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-3.30711058d-04,-1.02573996d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-5.77909579d-05,-6.28100405d-07/
      data fpp( 1,10,1),fpp( 1,10,2)/-6.17714258d-04,-1.85842480d-09/
      data fpp( 1,11,1),fpp( 1,11,2)/ 5.11179430d-04,-4.98465896d-07/
      data fpp( 1,12,1),fpp( 1,12,2)/ 4.54111603d-04, 1.37772201d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/-9.87380919d-04,-1.39584221d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-7.84010228d-04, 1.14599666d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.37058096d-02, 1.28905559d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/ 1.28981566d-02,-1.89941902d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/ 1.13781856d-02, 1.97722051d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-2.42503899d-03,-1.15306300d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-2.28136850d-05/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 1.91875547d-07/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-7.17537045d-05,-1.37751094d-07/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-6.63975217d-05,-5.88871171d-07/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-2.21521931d-04,-6.14764220d-07/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-6.41120321d-04, 9.29928053d-07/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-7.19272455d-04, 1.18505201d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-5.61115374d-04, 4.37863912d-07/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-2.50506456d-04,-7.64507658d-07/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-6.72752270d-05,-8.58332821d-08/
      data fpp( 2,10,1),fpp( 2,10,2)/-3.58571484d-04,-8.48159214d-07/
      data fpp( 2,11,1),fpp( 2,11,2)/ 3.89569712d-04, 1.04701383d-08/
      data fpp( 2,12,1),fpp( 2,12,2)/ 4.41276793d-04, 1.33427866d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/-6.87738161d-04,-8.17358478d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/-1.34797954d-03, 7.42060467d-07/
      data fpp( 2,15,1),fpp( 2,15,2)/ 7.15688073d-03, 1.34913429d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/ 8.18240103d-03,-1.06494321d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/ 6.89462871d-03, 1.62063856d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.25713631d-03,-1.14141102d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-2.34699449d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 1.43472443d-07/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 5.02808229d-05,-8.59448863d-08/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.46056116d-05,-4.29692898d-07/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-1.01708957d-04,-3.73283522d-07/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-3.79150307d-04, 6.38826984d-07/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-5.95046407d-04, 8.17975584d-07/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-4.28096189d-04, 2.41270678d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-1.52263118d-04,-4.51058296d-07/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 8.68918659d-05,-5.70374925d-08/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.92000194d-04,-8.38791733d-07/
      data fpp( 3,11,1),fpp( 3,11,2)/ 2.55541721d-04,-6.77955736d-08/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.15781224d-04, 1.19397403d-06/
      data fpp( 3,13,1),fpp( 3,13,2)/-3.26666435d-04,-4.13810054d-06/
      data fpp( 3,14,1),fpp( 3,14,2)/-2.82407159d-03,-2.21557188d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 5.51667455d-04, 7.52838804d-06/
      data fpp( 3,16,1),fpp( 3,16,2)/ 4.98223926d-03,-1.47398030d-06/
      data fpp( 3,17,1),fpp( 3,17,2)/-4.67004813d-05, 8.11153316d-06/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.55858423d-03,-7.25415233d-06/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-1.54909238d-05/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 7.48146553d-08/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.64336066d-05,-4.86293107d-08/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-8.10870241d-05,-2.58297413d-07/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-2.19955523d-04,-1.66181039d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-3.28752097d-04, 3.05021569d-07/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-2.86997008d-04, 5.36094764d-07/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-2.52269119d-04, 5.25993747d-08/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-6.87853023d-05,-1.22492263d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-2.61456068d-04,-7.86303235d-08/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.42863232d-05,-7.14986443d-07/
      data fpp( 4,11,1),fpp( 4,11,2)/ 9.84811220d-05,-1.81423904d-07/
      data fpp( 4,12,1),fpp( 4,12,2)/ 5.26544725d-04, 7.14682060d-07/
      data fpp( 4,13,1),fpp( 4,13,2)/-2.61977476d-06,-6.91304336d-07/
      data fpp( 4,14,1),fpp( 4,14,2)/-1.80777499d-03,-2.38946472d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/-9.15681200d-03, 2.26916320d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/-5.18573156d-03, 2.15081193d-06/
      data fpp( 4,17,1),fpp( 4,17,2)/-1.75408420d-03, 2.58558910d-06/
      data fpp( 4,18,1),fpp( 4,18,2)/-5.38565446d-05,-3.11516831d-06/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-7.06291584d-06/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 2.71800289d-08/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.92190475d-05,-1.93600577d-08/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.61515110d-04,-9.97397979d-08/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.88883047d-04,-6.16807507d-08/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 1.13096895d-04, 9.44628006d-08/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-7.09817311d-05, 1.99829548d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.61188950d-05, 5.42190059d-08/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-2.12916190d-06, 1.17294428d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 8.31724299d-04,-1.33396718d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/ 1.16116118d-04,-6.81707555d-07/
      data fpp( 5,11,1),fpp( 5,11,2)/ 9.13537690d-06, 7.24226940d-07/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.11211854d-04,-1.80720020d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 5.54783140d-04, 4.10457387d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 6.15322919d-04,-1.16650953d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 8.25799793d-03, 1.38098073d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.84259744d-03,-2.04813392d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/ 6.23489734d-04,-1.48327164d-06/
      data fpp( 5,18,1),fpp( 5,18,2)/ 1.13959041d-03,-1.10277953d-06/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-6.69610234d-07/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 2.51683873d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-2.68425832d-05, 8.86322544d-08/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-2.64973416d-04, 5.37871096d-08/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-5.34776666d-04,-4.29780693d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-3.51635482d-04,-4.46643385d-08/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-9.14760677d-05, 5.36438047d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.14606461d-04,-7.21087848d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-2.10980501d-05, 2.40191335d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-1.26544113d-03,-3.35456554d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.90978148d-04, 1.37434881d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-4.38226297d-05, 1.37170299d-07/
      data fpp( 6,12,1),fpp( 6,12,2)/-2.28097308d-04,-5.55030007d-07/
      data fpp( 6,13,1),fpp( 6,13,2)/-5.24512787d-04, 4.86949730d-07/
      data fpp( 6,14,1),fpp( 6,14,2)/-6.77516688d-04, 5.32310863d-08/
      data fpp( 6,15,1),fpp( 6,15,2)/-5.63277972d-03,-2.55874076d-07/
      data fpp( 6,16,1),fpp( 6,16,2)/-3.34705821d-03, 1.00626522d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/-1.61827473d-03, 2.62813211d-07/
      data fpp( 6,18,1),fpp( 6,18,2)/-2.00850509d-03,-1.87751806d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-3.61874097d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 1.18090277d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.61512854d-05,-1.61805538d-09/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 1.61578555d-04,-1.73368062d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 3.19023617d-04,-2.50347199d-08/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 1.89445032d-04, 3.47568576d-09/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 2.64860018d-05, 1.71319769d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 3.43069488d-05, 1.19964068d-08/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 5.53213624d-05, 4.44882396d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 7.45240211d-04,-6.93525990d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.46196475d-04, 2.35221565d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.49551417d-05,-1.33360268d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.83601086d-04, 2.82195074d-08/
      data fpp( 7,13,1),fpp( 7,13,2)/ 2.40068007d-04,-2.07517762d-07/
      data fpp( 7,14,1),fpp( 7,14,2)/ 4.02743834d-04, 9.33851539d-07/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.33952094d-03,-1.02588840d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 9.89635407d-04, 7.63702044d-07/
      data fpp( 7,17,1),fpp( 7,17,2)/ 8.09609200d-04,-2.52919781d-07/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.47282994d-03,-5.20229197d-08/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-1.98988540d-07/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 5.45784469d-09/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-4.16255822d-06,-1.91568939d-09/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-4.29408030d-05,-3.79508715d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-8.85178010d-05,-2.49039620d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-6.29446472d-05, 2.54109353d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-2.64679397d-05,-1.06739779d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 4.69786657d-05, 2.81548181d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-4.41873995d-05,-4.14529441d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.81919718d-04,-2.23736404d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/-2.18077532d-05, 4.23985613d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-3.99793719d-06, 6.14215921d-09/
      data fpp( 8,12,1),fpp( 8,12,2)/-8.39070368d-05,-4.89671981d-08/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.18959242d-04, 3.97266332d-08/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.89458648d-04, 4.00606653d-08/
      data fpp( 8,15,1),fpp( 8,15,2)/-2.50904034d-04, 1.18030706d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/-3.06683415d-04,-1.10183488d-07/
      data fpp( 8,17,1),fpp( 8,17,2)/-2.47362066d-04, 4.36703247d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/-7.00414663d-04,-5.62629499d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.12018525d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 2.10767568d-09/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 4.98947509d-07,-1.21535136d-09/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 7.78465725d-06,-3.24627023d-09/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 1.58475873d-05,-9.79956773d-09/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.43335566d-05,-2.95554589d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.45857570d-05, 9.20214032d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-3.98216115d-05,-2.53015375d-09/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-9.71764531d-07,-1.59007881d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 4.72386597d-05,-1.78666937d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.36546269d-06,-9.26324372d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.03660706d-06, 5.23964425d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.52270609d-05,-8.95333279d-09/
      data fpp( 9,13,1),fpp( 9,13,2)/ 5.36895915d-06, 1.94168887d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/ 4.69075774d-06, 3.32857782d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/ 7.29519687d-06, 1.54399987d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 2.82982538d-05, 9.54227161d-10/
      data fpp( 9,17,1),fpp( 9,17,2)/ 3.58390636d-05,-2.52569073d-08/
      data fpp( 9,18,1),fpp( 9,18,2)/ 2.27228714d-04,-7.92659791d-09/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-3.03670104d-09/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 2.80627759d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.55634178d-08, 1.38744481d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/-1.88357025d-06,-2.35605684d-09/
      data fpp(10, 4,1),fpp(10, 4,2)/-4.48386139d-06, 2.03678256d-09/
      data fpp(10, 5,1),fpp(10, 5,2)/-6.72834613d-06, 6.20892660d-09/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.19233012d-05, 3.12751102d-09/
      data fpp(10, 7,1),fpp(10, 7,2)/ 5.37550167d-06,-7.18970698d-10/
      data fpp(10, 8,1),fpp(10, 8,2)/-1.91006675d-07,-1.22516282d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-9.35612020d-06, 2.57254836d-08/
      data fpp(10,10,1),fpp(10,10,2)/ 6.00264664d-07,-6.66503062d-08/
      data fpp(10,11,1),fpp(10,11,2)/ 8.91474055d-08, 4.28757414d-08/
      data fpp(10,12,1),fpp(10,12,2)/-6.72766417d-06,-1.48526592d-08/
      data fpp(10,13,1),fpp(10,13,2)/-4.62725667d-06, 1.65348954d-08/
      data fpp(10,14,1),fpp(10,14,2)/-7.54294926d-06, 2.71307762d-09/
      data fpp(10,15,1),fpp(10,15,2)/-8.03357381d-06, 8.61279412d-09/
      data fpp(10,16,1),fpp(10,16,2)/-1.19530538d-05,-7.16425409d-09/
      data fpp(10,17,1),fpp(10,17,2)/-1.08361579d-05,-3.95577774d-09/
      data fpp(10,18,1),fpp(10,18,2)/-6.02788108d-05,-1.01263493d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 2.00631747d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 4.70519190d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 1.63306163d-07, 2.58961619d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/ 3.49623750d-07,-3.06365668d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/ 8.87858260d-07, 3.66501051d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.17982793d-06, 4.03614621d-10/
      data fpp(11, 6,1),fpp(11, 6,2)/ 2.50744768d-06, 7.20531002d-10/
      data fpp(11, 7,1),fpp(11, 7,2)/-1.48039518d-06,-3.28573863d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 5.35791230d-07, 4.22423520d-10/
      data fpp(11, 9,1),fpp(11, 9,2)/ 3.38582111d-06,-4.40395545d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 2.56440404d-06, 5.19339828d-09/
      data fpp(11,11,1),fpp(11,11,2)/-7.93196686d-07,-1.03696377d-08/
      data fpp(11,12,1),fpp(11,12,2)/ 2.83595835d-07, 1.82851524d-08/
      data fpp(11,13,1),fpp(11,13,2)/-2.45993245d-06,-2.07709720d-08/
      data fpp(11,14,1),fpp(11,14,2)/-3.18960722d-07, 1.07987357d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 8.39098371d-07, 1.57602940d-09/
      data fpp(11,16,1),fpp(11,16,2)/ 2.71396154d-06, 8.97146760d-10/
      data fpp(11,17,1),fpp(11,17,2)/ 2.70556786d-06, 8.35383565d-10/
      data fpp(11,18,1),fpp(11,18,2)/ 1.50865292d-05,-4.23868102d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-7.88065949d-09/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-9.20960053d-11/
      data fpp(12, 2,1),fpp(12, 2,2)/-3.76612325d-08, 1.84192011d-10/
      data fpp(12, 3,1),fpp(12, 3,2)/-1.14924750d-07,-6.44672037d-10/
      data fpp(12, 4,1),fpp(12, 4,2)/-2.67571652d-07, 2.39449614d-09/
      data fpp(12, 5,1),fpp(12, 5,2)/-3.90965586d-07,-2.93331251d-09/
      data fpp(12, 6,1),fpp(12, 6,2)/-5.06489535d-07, 3.33875390d-09/
      data fpp(12, 7,1),fpp(12, 7,2)/ 5.46079035d-07,-4.42170310d-09/
      data fpp(12, 8,1),fpp(12, 8,2)/ 4.47841754d-07, 2.34805850d-09/
      data fpp(12, 9,1),fpp(12, 9,2)/ 1.28357790d-08, 1.02946912d-09/
      data fpp(12,10,1),fpp(12,10,2)/ 5.42119193d-07,-4.65934969d-10/
      data fpp(12,11,1),fpp(12,11,2)/ 8.36393371d-08,-5.16572924d-09/
      data fpp(12,12,1),fpp(12,12,2)/-4.06719167d-07, 9.12885193d-09/
      data fpp(12,13,1),fpp(12,13,2)/-5.33013509d-07,-7.34967849d-09/
      data fpp(12,14,1),fpp(12,14,2)/-7.81207856d-07, 2.26986204d-09/
      data fpp(12,15,1),fpp(12,15,2)/-7.22819674d-07,-1.72976966d-09/
      data fpp(12,16,1),fpp(12,16,2)/-7.02792307d-07, 4.64921661d-09/
      data fpp(12,17,1),fpp(12,17,2)/-5.86113572d-07,-4.86709678d-09/
      data fpp(12,18,1),fpp(12,18,2)/-3.06730583d-06, 2.81917051d-09/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 5.59041475d-09/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 0.00000000d+00/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.68669384d-07, 0.00000000d+00/
      data fpp(13, 3,1),fpp(13, 3,2)/ 1.69962375d-07, 0.00000000d+00/
      data fpp(13, 4,1),fpp(13, 4,2)/ 3.58785826d-07, 0.00000000d+00/
      data fpp(13, 5,1),fpp(13, 5,2)/ 7.32982793d-07, 0.00000000d+00/
      data fpp(13, 6,1),fpp(13, 6,2)/ 1.31574477d-06, 0.00000000d+00/
      data fpp(13, 7,1),fpp(13, 7,2)/ 6.01960482d-07, 0.00000000d+00/
      data fpp(13, 8,1),fpp(13, 8,2)/ 6.38579123d-07, 0.00000000d+00/
      data fpp(13, 9,1),fpp(13, 9,2)/ 5.18582111d-07, 0.00000000d+00/
      data fpp(13,10,1),fpp(13,10,2)/-1.70855960d-06, 0.00000000d+00/
      data fpp(13,11,1),fpp(13,11,2)/ 1.04568033d-06, 0.00000000d+00/
      data fpp(13,12,1),fpp(13,12,2)/ 1.67835958d-06, 0.00000000d+00/
      data fpp(13,13,1),fpp(13,13,2)/ 3.42900675d-06, 0.00000000d+00/
      data fpp(13,14,1),fpp(13,14,2)/ 1.75310393d-06, 0.00000000d+00/
      data fpp(13,15,1),fpp(13,15,2)/ 8.48909837d-07, 0.00000000d+00/
      data fpp(13,16,1),fpp(13,16,2)/ 6.01396154d-07, 0.00000000d+00/
      data fpp(13,17,1),fpp(13,17,2)/ 4.05556786d-07, 0.00000000d+00/
      data fpp(13,18,1),fpp(13,18,2)/ 2.70865292d-06, 0.00000000d+00/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 0.00000000d+00/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 5.24692724d-06, 1.43073528d-06/
      data fpppp( 1, 3),fpppp( 1, 4)/-1.51887926d-06,-7.50151478d-06/
      data fpppp( 1, 5),fpppp( 1, 6)/ 8.66630066d-06, 1.49649590d-06/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.06272498d-06, 4.23119190d-06/
      data fpppp( 1, 9),fpppp( 1,10)/-2.30161619d-05, 3.78628517d-05/
      data fpppp( 1,11),fpppp( 1,12)/-2.71062258d-05,-5.95639534d-07/
      data fpppp( 1,13),fpppp( 1,14)/-5.35766979d-05, 3.13594224d-04/
      data fpppp( 1,15),fpppp( 1,16)/-3.43613248d-04, 1.43010395d-04/
      data fpppp( 1,17),fpppp( 1,18)/-2.71167411d-04, 2.04664029d-04/
      data fpppp( 1,19) /             4.26207113d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 3.19174799d-06, 6.81610214d-07/
      data fpppp( 2, 3),fpppp( 2, 4)/-1.29159561d-06,-5.14406328d-06/
      data fpppp( 2, 5),fpppp( 2, 6)/ 5.99940987d-06, 1.63319909d-06/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.64634668d-06, 9.28524440d-07/
      data fpppp( 2, 9),fpppp( 2,10)/-1.30031058d-05, 2.26122496d-05/
      data fpppp( 2,11),fpppp( 2,12)/-1.50796454d-05,-4.07971479d-06/
      data fpppp( 2,13),fpppp( 2,14)/-3.94448175d-05, 1.89985399d-04/
      data fpppp( 2,15),fpppp( 2,16)/-1.70590680d-04, 4.36169213d-05/
      data fpppp( 2,17),fpppp( 2,18)/-1.42674562d-04, 1.15241766d-04/
      data fpppp( 2,19) /             2.46241577d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-4.81061730d-07,-8.22891105d-07/
      data fpppp( 3, 3),fpppp( 3, 4)/-7.84735902d-07,-2.07652672d-06/
      data fpppp( 3, 5),fpppp( 3, 6)/ 2.32359137d-08, 5.67629806d-06/
      data fpppp( 3, 7),fpppp( 3, 8)/ 2.42350874d-07,-1.12730317d-07/
      data fpppp( 3, 9),fpppp( 3,10)/-1.99211484d-06, 3.83903097d-08/
      data fpppp( 3,11),fpppp( 3,12)/-6.55454439d-07, 2.38530599d-06/
      data fpppp( 3,13),fpppp( 3,14)/-5.10469993d-05, 9.05052410d-05/
      data fpppp( 3,15),fpppp( 3,16)/ 4.14146876d-05,-1.92874026d-04/
      data fpppp( 3,17),fpppp( 3,18)/ 1.62510723d-04,-5.91153973d-05/
      data fpppp( 3,19) /            -1.15881269d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.61168367d-07,-4.04352181d-07/
      data fpppp( 4, 3),fpppp( 4, 4)/-1.01461156d-06, 9.89354091d-09/
      data fpppp( 4, 5),fpppp( 4, 6)/ 2.77935286d-06,-2.09420519d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/ 5.17583583d-06,-9.68378245d-06/
      data fpppp( 4, 9),fpppp( 4,10)/ 1.09900190d-05,-7.88586280d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 1.24892942d-05,-2.31535447d-05/
      data fpppp( 4,13),fpppp( 4,14)/ 2.26911982d-05,-1.44170691d-04/
      data fpppp( 4,15),fpppp( 4,16)/ 2.21358658d-04,-6.20568949d-05/
      data fpppp( 4,17),fpppp( 4,18)/-5.49706306d-06,-1.98400351d-05/
      data fpppp( 4,19) /            -1.39250632d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 2.72745629d-06, 9.97500689d-07/
      data fpppp( 5, 3),fpppp( 5, 4)/ 6.67161865d-07,-4.56183568d-06/
      data fpppp( 5, 5),fpppp( 5, 6)/-6.09064528d-07, 6.50054539d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/-9.12236191d-06, 2.36679813d-05/
      data fpppp( 5, 9),fpppp( 5,10)/-3.44234721d-05, 2.10582085d-05/
      data fpppp( 5,11),fpppp( 5,12)/-1.32917154d-05, 3.13066637d-05/
      data fpppp( 5,13),fpppp( 5,14)/-6.47544058d-05, 1.91383647d-04/
      data fpppp( 5,15),fpppp( 5,16)/-2.45852067d-04, 6.85400900d-05/
      data fpppp( 5,17),fpppp( 5,18)/ 4.34692731d-05,-1.83046792d-05/
      data fpppp( 5,19) /            -6.95920213d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/-4.11733350d-06,-1.69629715d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.77477287d-06, 6.89504364d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.37126434d-06,-7.75900723d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.26673762d-05,-3.59121691d-05/
      data fpppp( 6, 9),fpppp( 6,10)/ 5.07102110d-05,-2.78003116d-05/
      data fpppp( 6,11),fpppp( 6,12)/ 4.85258771d-06,-1.14958510d-05/
      data fpppp( 6,13),fpppp( 6,14)/ 3.44023685d-05,-1.17508928d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 1.47497796d-04,-3.80231857d-05/
      data fpppp( 6,17),fpppp( 6,18)/-2.88213351d-05, 2.61676963d-05/
      data fpppp( 6,19) /             6.80746761d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 2.59678918d-06, 1.02986500d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/ 1.04030987d-06,-4.47003693d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/-3.81580933d-07, 3.99353391d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/-5.34575606d-06, 1.81810983d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-2.72443712d-05, 1.34586312d-05/
      data fpppp( 7,11),fpppp( 7,12)/ 1.47799032d-06,-1.37735584d-06/
      data fpppp( 7,13),fpppp( 7,14)/-2.69930837d-06, 1.85471237d-05/
      data fpppp( 7,15),fpppp( 7,16)/-2.50431098d-05, 4.42555747d-06/
      data fpppp( 7,17),fpppp( 7,18)/ 1.75324393d-05,-2.39604981d-05/
      data fpppp( 7,19) /            -4.98534872d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/-6.65255320d-07,-2.64298626d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/-3.54491372d-07, 1.27433892d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/-4.73855209d-07, 1.27529513d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/-2.40913145d-06,-1.51552957d-06/
      data fpppp( 8, 9),fpppp( 8,10)/ 5.67727453d-06,-3.32291158d-06/
      data fpppp( 8,11),fpppp( 8,12)/-9.23757118d-07, 1.15480512d-06/
      data fpppp( 8,13),fpppp( 8,14)/-1.00404968d-06, 7.34561508d-07/
      data fpppp( 8,15),fpppp( 8,16)/-1.39095511d-06, 5.16921917d-06/
      data fpppp( 8,17),fpppp( 8,18)/-1.23798777d-05, 1.36078550d-05/
      data fpppp( 8,19) /             2.71564935d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 1.36256437d-07, 5.12808326d-08/
      data fpppp( 9, 3),fpppp( 9, 4)/ 6.58259672d-08,-2.67951483d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/ 4.31362320d-07,-1.35152393d-06/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.69515926d-06, 1.66319834d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-1.79880396d-06, 1.22002320d-06/
      data fpppp( 9,11),fpppp( 9,12)/-2.09173028d-08,-4.29050944d-07/
      data fpppp( 9,13),fpppp( 9,14)/ 2.94207750d-07,-1.96986037d-07/
      data fpppp( 9,15),fpppp( 9,16)/ 6.90694830d-07,-1.46187622d-06/
      data fpppp( 9,17),fpppp( 9,18)/ 4.34907520d-06,-4.90349416d-06/
      data fpppp( 9,19) /            -9.85220044d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-3.17480103d-08,-1.44804641d-08/
      data fpppp(10, 3),fpppp(10, 4)/-2.14767384d-08, 5.64503592d-08/
      data fpppp(10, 5),fpppp(10, 6)/-1.82976314d-07, 4.98426679d-07/
      data fpppp(10, 7),fpppp(10, 8)/-4.61104930d-07,-2.59256306d-08/
      data fpppp(10, 9),fpppp(10,10)/ 3.48891142d-07,-2.22349033d-07/
      data fpppp(10,11),fpppp(10,12)/-8.75451389d-08, 1.94187929d-07/
      data fpppp(10,13),fpppp(10,14)/-1.54173433d-07, 1.21539798d-07/
      data fpppp(10,15),fpppp(10,16)/-1.86481679d-07, 4.18655590d-07/
      data fpppp(10,17),fpppp(10,18)/-1.18595812d-06, 1.29160396d-06/
      data fpppp(10,19) /             2.60283010d-06 /
      data fpppp(11, 1),fpppp(11, 2)/-2.06252105d-09,-1.76272511d-09/
      data fpppp(11, 3),fpppp(11, 4)/ 1.04941070d-08,-1.90986874d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 5.11247522d-08,-1.23261317d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 1.22992758d-07,-8.46796063d-09/
      data fpppp(11, 9),fpppp(11,10)/-3.90903076d-08,-5.54576257d-08/
      data fpppp(11,11),fpppp(11,12)/ 1.08749791d-07,-1.13477944d-07/
      data fpppp(11,13),fpppp(11,14)/ 1.15942738d-07,-5.72230054d-08/
      data fpppp(11,15),fpppp(11,16)/ 5.39745254d-08,-1.15666852d-07/
      data fpppp(11,17),fpppp(11,18)/ 2.95697473d-07,-3.23761740d-07/
      data fpppp(11,19) /            -6.48699941d-07 /
      data fpppp(12, 1),fpppp(12, 2)/-1.79118833d-10,-1.14208894d-10/
      data fpppp(12, 3),fpppp(12, 4)/-1.74018269d-09, 2.55193659d-09/
      data fpppp(12, 5),fpppp(12, 6)/-6.71238557d-09, 2.47698047d-08/
      data fpppp(12, 7),fpppp(12, 8)/-2.22812821d-08,-4.69302753d-09/
      data fpppp(12, 9),fpppp(12,10)/ 2.08472706d-08,-2.08386915d-08/
      data fpppp(12,11),fpppp(12,12)/ 3.24169923d-09, 5.95917569d-09/
      data fpppp(12,13),fpppp(12,14)/-5.23455226d-09, 7.66503313d-09/
      data fpppp(12,15),fpppp(12,16)/-7.03062860d-09, 1.81558324d-08/
      data fpppp(12,17),fpppp(12,18)/-5.97936189d-08, 6.51463835d-08/
      data fpppp(12,19) /             1.32117971d-07 /
      data fpppp(13, 1),fpppp(13, 2)/ 1.60986248d-08, 8.16397377d-09/
      data fpppp(13, 3),fpppp(13, 4)/-6.31645127d-09, 2.11333287d-09/
      data fpppp(13, 5),fpppp(13, 6)/ 8.98553075d-09,-2.55415554d-08/
      data fpppp(13, 7),fpppp(13, 8)/ 1.53879152d-08, 9.01407025d-09/
      data fpppp(13, 9),fpppp(13,10)/-6.08411354d-08, 1.07921790d-07/
      data fpppp(13,11),fpppp(13,12)/-7.19631247d-08, 5.26370687d-08/
      data fpppp(13,13),fpppp(13,14)/-7.15070749d-08, 2.77982311d-08/
      data fpppp(13,15),fpppp(13,16)/ 6.61667473d-09,-1.48641056d-08/
      data fpppp(13,17),fpppp(13,18)/ 5.59402066d-08,-5.89605909d-08/
      data fpppp(13,19) /            -1.20802786d-07 /
 

      c(1,1)=fpppp(ix+1,iy+1)
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
