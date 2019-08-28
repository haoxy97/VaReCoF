c***********************************************************************

      real*8 function ch3_h_hirst(r, rpar, ipar)

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
      common /phico2/ r0,ap,bp,cp,phich4   
      common /vinfo/radmake,sqrt3

c
c vectors for primary coordinates 
c

      dimension r(nfrag,natommx,ndim)
      dimension cm(nfrag,ndim),rij(natommx,natommx)
      dimension qq(5)

c
c vectors from mcpot calling routine
c

      dimension natom(nfrag)
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
      do 1500 idim = 1 , ndim
         rbvec(idim) = r(2,1,idim) - r(1,1,idim)
         rbond2 = rbond2 + rbvec(idim)**2
         rax1v(idim) = r(1,2,idim) - r(1,1,idim)
         rax12 = rax12 + rax1v(idim)**2
         rbx1v(idim) = r(2,1,idim) - r(1,2,idim)
         rbx12 = rbx12 + rbx1v(idim)**2
 1500 continue
      rbond = sqrt(rbond2)
      rax1 = sqrt(rax12)
      rbx1 = sqrt(rbx12)
      ray1 = sqrt(ray12)
      rby1 = sqrt(rby12)
c     
c     calculate the bending angles.
c     

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

      theta = the1*180.0d0/pi
      phi = tau3*180.0d0/pi
      alpha = 90.0d0
      rch = rbond

c     if (theta.lt.90.0d0) then 
c        vtot = 1.0d20
c        rmepi = 0.0d0
c     else

c         call hpch3(rch,theta,phi,alpha,vtot)
      call pathinput
      alf=0.0d0
      rchang = rch*0.529177d0
      vtot = vhinder(rchang,theta,phi,alf) / 627.5
          
c         vtotmin = 1.0d20
c         do ialp = 1 , 41
c	      alpha = 70.0d0+0.5d0*float(ialp-1)	
c            call hpch3(rch,theta,phi,alpha,vtot)
c	      if (vtot.lt.vtotmin) vtotmin = vtot
c	write (6,*) 'pot test',vtot*219474.,alpha,rch,theta,phi
c	   enddo
c	   vtot = vtotmin

c        ttmp = 1000.0d0*dkbltz/cautoerg
c 300k	   qalpref = 0.2431d0
c   qalpref = 0.3859d0
c   betatmp = 1.0d0/ttmp
c        qalp = 0.0d0
c        vtotmin = 1.0d20
c        do ialp = 1 , 121
c      alpha = 60.0d0+0.5d0*float(ialp-1)	
c           call hpch3(rch,theta,phi,alpha,vtot)
c      if (vtot.lt.vtotmin) vtotmin = vtot
c      qalp = qalp + 0.5d0*exp(-vtot*betatmp)*pi/180.0d0
c	write (6,*) 'qalp sum',qalp,vtot,alpha,betatmp
c   enddo
c   qalp = qalp*exp(vtotmin*betatmp)
c   vtot = vtotmin
c	write (6,*) 'qalp test',vtot*219474.,qalp,rch,theta,phi,betatmp
c        rmepi = qalp/qalpref
c         rmepi = 1.0d0
c        write (6,*) 'qalp result',qalp,qalp/qalpref,vtot*cautoicm
c     endif
c      write (6,*) 'pot test',energy*219474.,rch,theta,phi,alpha
      
      ch3_h_hirst = vtot
      return
      end


c to test hase/hirst potential for ch4 restricted to c3v symmetry for
c the ch3 fragment
c
c the code uses rch for distance in angstroms of attacking h to c
c               theta for the angle of rch with c3v axis of ch3 where
c                     theta = 0 has ch3 splayed away from attacking h
c               phi for the azimuthal angle of rch with phi = 0 
c                     meaning the attacking h, the c, and one h on the
c                     ch3 are in the same plane
c               alf for the splay angle of the ch3
c               beta for the h-c-h bond angle in ch3 (all three such
c                     angles are common in c3v symmetry
c
c     implicit real*8 (a-h,o-z)
c
c set up potential
c
c
c loop over rch,theta,phi for value of potential
c
c2     read(5,*,end=3)rau,tnow,pnow,alf
c     rnow = rau * 0.529177d0
c     write(6,1) rau,tnow,pnow,alf,v
c1     format(4f8.2,1f10.3)
c      go to 2
c3     stop
c      end

c input for hase/hirst (+struan's fit to hirst's data) for ch3+h
c
      subroutine pathinput
c***********************************************************************
c to read in and initialise data, and to set up potential functions.
c see input file for the definition of variables.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /struan1/a3,a2,aa1,a0,de,re
      common /phico1/ a1,b1,c1,a4,b4,c4,fpch4,gpch4   
      common /phico2/ r0,ap,bp,cp,phich4   
      common /betainfo/bch4,ab,alfb,r0b
      common /dswitch/ad,r0d
      common /tswitch/at,r0t
      common /vinfo/radmake,sqrt3
      common /tconst/fch3,gch3,hch3,fch43,gch43,hch43,fd0,hd0,gn40
c
c title info
c
c      write(6,*) ' '
c      write(6,*) ' ch3+h hase hirst potential restricted to c3v',
c     *   ' geo for ch3'
c      write(6,*) ' references: jacs 109,2916(1987) + jpc 88,1339 (1984)'
c      write(6,*) ' '
c
c input hase/hirst potential parameters        
c
      a3 = -1.017234d-3 
      a2 = 7.7738886d-2 
      aa1 = 7.703640d-2 
      a0 = 1.686690d0  
      de = 103.432d0 
      re = 1.0015d0
c      write(6,*) ' parameters for rch dependent mo along rxn path'
c      write(6,*) ' a3,a2,aa1,a0,de,re'
c      write(6,*) a3,a2,aa1
c      write(6,*) a0,de,re

      phich4 = 109.47d0  
      r0 = 1.086d0
      ap = 0.47418183d0
      bp = -2.6703827d-4
      cp = 0.0d0
      fpch4 = 0.5938d0
      a1 = 0.50453752d0
      b1 = 0.41910920d0 
      c1 = 0.69892603d0
      gpch4 = -0.09964d0 
      a4 = 0.38387689d0
      b4 = -0.16991527d0
      c4 = 0.97118670d0
c      write(6,109) phich4,r0,ap,bp,cp,fpch4,a1,b1,c1,gpch4,a4,b4,c4   
c 109  format(' equilib. splay angle with switch info :'/
c     *       ' phich4 =',f7.2,' r0 =',g15.9,' ap =',g15.9,' bp =',g15.9/
c     *       ' cp=',g15.9/
c     *       ' fp (in mdynes/rad**2) with switch info :'/
c     *       ' fpch4 =',f9.5,' a1 =',g15.9,' b1 =',g15.9,' c1 =',g15.9/
c     *       ' gp (in mdynes/rad**2) with switch info :'/
c     *       ' gpch4 =',f9.5,' a4 =',g15.9,' b4 =',g15.9,' c4 =',g15.9)

      bch4 = 109.47d0
      ab   = .77029473d0
      alfb = -.0046616956
      r0b  = 4.3805318
c      write(6,113) bch4,ab,alfb,r0b
c 113  format(' equilib. ch3 bond angle with switch info :'/
c     *       ' bch4 =',f7.2,' ab =',g15.9,' alfb =',g15.9,' r0b=',g15.9)

      fd0 = .0436d0
      hd0 = .0854d0
      gn40 = .2242d0
      ad  = .14191474d0
      r0d = -.30684503d0
c      write(6,115) fd0,hd0,gn40,ad,r0d   
c 115  format(' fd,hd,gn4 (in mdynes/rad**n) with switch info :'/
c     *       ' fd =',f9.5,' hd =',f9.5,' gn4 =',f9.5,' ad =',g15.9/
c     *       ' r0d =',g15.9)

      fch3 = .4543d0
      gch3 = -.1232d0
      hch3 = -.0101d0
      fch43 =  .5938d0 - fch3
      gch43 = -.0903d0 - gch3
      hch43 =  .0120d0 - hch3
      at  = 1.0147402d-7
      r0t = -12.362798
c      write(6,117) fch3,gch3,hch3,fch43,gch43,hch43,at,r0t
c117   format(' fch3,gch3,hch3 (in mdynes/rad**n) =',3f9.5/
c     *       ' f(g,h)ch4-ch3  (in mdynes/rad**n) =',3f9.5/
c     *       ' associated switch info: at =',g15.9,' r0t =',g15.9)
c      write(6,*) ' '

c
c set constants and convert to atomic units           
c
      aeudya = 0.2293676              
      fpch4  = fpch4*aeudya            
      gpch4  = gpch4*aeudya            
      fch3   = fch3*aeudya
      gch3   = gch3*aeudya
      hch3   = hch3*aeudya
      fch43  = fch43*aeudya
      gch43  = gch43*aeudya
      hch43  = hch43*aeudya
      fd0    = fd0*aeudya
      hd0    = hd0*aeudya
      gn40   = gn40*aeudya
      pi      = acos(-1.d0)
      radmake = pi/180.d0
      sqrt3   = sqrt(3.d0)
      return
      end
c to obtain the hindering potential for given angles and rmep
c
      double precision function vhinder(rch,theta,phi,alf)
      implicit real*8 (a-h,o-z)
      common /struan1/a3,a2,aa1,a0,de,re
      common /phico1/ a1,b1,c1,a4,b4,c4,fpch4,gpch4   
      common /phico2/ r0,ap,bp,cp,phich4   
      common /betainfo/bch4,ab,alfb,r0b
      common /dswitch/ad,r0d
      common /tswitch/at,r0t
      common /vinfo/radmake,sqrt3
      common /tconst/fch3,gch3,hch3,fch43,gch43,hch43,fd0,hd0,gn40
      dimension pphi(3)
c
c the code uses rch for distance in angstroms of attacking h to c
c               theta for the angle of rch with c3v axis of ch3 where
c                     theta = 0 has ch3 splayed away from attacking h
c               phi for the azimuthal angle of rch with phi = 0 
c                     meaning the attacking h, the c, and one h on the
c                     ch3 are in the same plane
c               alf for the splay angle of the ch3
c               beta for the h-c-h bond angle in ch3 (all three such
c                     angles are common in c3v symmetry
c
c               delv is the energy in kcal along the mep
c
c set up angles and trig functions of them
c
      alf0 = 0.0d0    
      if(rch.le.10.0d0) alf0 = (phich4 - 90.d0)*
     x      (1.0d0 - dtanh(ap*(rch-r0)*dexp(bp*(rch-cp)**3)))
      splay0 = (alf0+90.d0)*radmake
      splay = (alf+90.d0)*radmake
      beta0 = 120.d0 + (bch4 - 120.d0)*(1.d0 - tanh(ab*
     x        (rch-r0)*exp(alfb*(rch-r0b)**3)))
      beta0 = beta0*radmake
      beta = acos(.5d0*(3.d0*cos(splay)**2 - 1.d0))
      cphi = cos(phi*radmake)
      sphi = sin(phi*radmake)
c
c compute the three phi's wrt the attacking atom (e.g., phi4)
c
      anow = cos(splay)*cos(theta*radmake)
      bnow = sin(splay)*sin(theta*radmake)
      pphi(1) = acos(anow+bnow*cphi)
      pphi(2) = acos(anow-.5d0*bnow*(cphi+sqrt3*sphi))
      pphi(3) = acos(anow-.5d0*bnow*(cphi-sqrt3*sphi))
c
c make fp and gp
c
      s1= 1.0d0 - dtanh(a1*(rch-r0)*(rch-b1)**c1)       
      fp=s1*fpch4 
      s4= 1.0d0 - dtanh(a4*(rch-r0)*(rch-b4)**c4)       
      gp=s4*gpch4 
c
c make fd and hd and gn4
c
      s3 = 1.d0 - tanh(ad*(rch-r0)*(rch-r0d)**2)
      fd = (1.d0 - s3)*fd0
      hd = (1.d0 - s3)*hd0
      gn4 = s3*gn40
c
c make ft,gt, and ht
c
      s2 = 1.d0 - tanh(at*(rch-r0)*(rch-r0t)**6)
      ft = fch3 + fch43*s2
      gt = gch3 + gch43*s2
      ht = hch3 + hch43*s2
c
c compute reaction path energy
c
      bep  = a3*(rch-re)**3 + a2*(rch-re)**2 + aa1*(rch-re) + a0
      delv = de*(1.d0 - exp(-bep*(rch-re)))**2 - de
c
c compute potential change from rxn path from all the components
c
      v = 0.d0     
      do i=1,3  
         ang = (pphi(i)-splay0)           
         v = v + 0.5d0*fp*ang**2 + gp*(ang**3 + ang**4)  
      enddo
      cc4 = (beta-beta0)*((beta-beta0)**2
     x                        + (pphi(1)-splay0)*(pphi(2)-splay0)
     x                        + (pphi(1)-splay0)*(pphi(3)-splay0)
     x                        + (pphi(2)-splay0)*(pphi(3)-splay0))
      v = v + gn4*cc4
      v = v + 3.d0*(fd*(splay-splay0)**2 + hd*(splay-splay0)**4)
      v = v + 3.d0*(.5d0*ft*(beta-beta0)**2 + gt*(beta-beta0)**3
     x              + ht*(beta-beta0)**4)
      yy = v*627.5095d0
      vhinder = yy + delv
      return
      end
