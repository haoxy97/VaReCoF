c***********************************************************************

      real*8 function c2h3_o(r, rpar, ipar)

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
      dimension qq(5)

c
c vectors from mcpot calling routine
c

      dimension ift(nfrag),natom(nfrag)

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

      ift(1) = 2
      ift(2) = 0
      iel = 1
      ieff = 1
      ich = 1
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
      ro = rbond
      if (ro.lt.3.8d0) ro=3.8d0
      if (ro.gt.15.0d0) ro=15.d0

      call spl_c2h3po_1(ro,theta,phi,eroot1)
      c2h3_o = eroot1
c     call spl_c2h3po_2(ro,theta,phi,eroot2)
c     call spl_c2h3po_3(ro,theta,phi,eroot3)
c     vtot = eroot2
c     vtot = min(eroot1,eroot2,eroot3)
c     vmin = min(eroot1,eroot2,eroot3)
c     vmax = max(eroot1,eroot2,eroot3)
c     vtot = eroot1+eroot2+eroot3-vmin-vmax
c     thedvd=150.d0
c     if (abs(phi).gt.90.0d0) then
c        if (theta.gt.thedvd) then
c           vtot = 0.0d0
c        endif
c     else
c        vtot = 0.0d0
c     endif
c     if (abs(phi).gt.90.0d0) then
c        if (theta.lt.thedvd) then
c           vtot = 0.0d0
c        endif
c     endif

      return
      end

c     implicit real*8 (a-h,o-z)
c
c     test program to call h2cch + o 3d surfaces
c
c
c  ***  coordinate system  ***
c
c                 o
c                /                  r(cha) = 2.05 au
c               /                   r(chb) = 2.05 au
c       ha     /                    r(chc) = 2.05 au
c        \    /                     r(cc)  = 2.50 au
c         ca=m=cb                   ha-c-c = 120.0 deg
c        /      \                   hb-c-c = 135.0 deg
c       hc       hb                 hc-c-c = 120.0 deg
c
c       m is the midpoint of the cc bond
c
c  ***   input   ***
c
c     ro :     r(m-o) au
c     alpha :  o-m-ca angle  (degrees)
c     phi :    dihedral angle between omca and hacacb planes (degrees)
c
c     ranges :     3.8 < rd < 15.0
c                    0 < alpha < +180
c                  no limits on phi
c
c
c  ***   output   ***
c
c     energy - potential energy relative to o+c2h3 in au
c
c1     read (5,*,end=100) ro,alpha,phi

c      call spl_c2h3po_1(ro,alpha,phi,eroot1)
c      call spl_c2h3po_2(ro,alpha,phi,eroot2)
c      call spl_c2h3po_3(ro,alpha,phi,eroot3)

c      write( 6,10) ro,alpha,phi,eroot1,eroot2,eroot3
c10    format(3f8.2,3f15.8)
c      go to 1
c100   stop
c      end

      subroutine spl_c2h3po_1(ro,alpha,phi,energy)
      implicit real*8 (a-h,o-z)

      data degrad / 57.29577951308232d 00 /
c
c     find location in spline grids
c
      call fnd_grd11(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1)
c
      call a_1(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,a1)
      call b_1(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,b1)
      call c_1(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,c1)
      call bsc(ro,ebsc)
c
      p = phi / degrad
      eraw = a1 + b1*cos(p) + c1*cos(2.0d0* p)
      energy = eraw+ebsc

c      write(6,1) a1,b1,c1,eraw,energy
c1     format(' a1,b1,c1 :',5f20.10)

      return
      end
      subroutine spl_c2h3po_2(ro,alpha,phi,energy)
      implicit real*8 (a-h,o-z)

      data degrad / 57.29577951308232d 00 /
c
c     find location in spline grids
c
      call fnd_grd11(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1)
c
      call a_2(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,a2)
      call b_2(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,b2)
      call c_2(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,c2)
      call bsc(ro,ebsc)
c
      p = phi / degrad
      eraw = a2 + b2*cos(p) + c2*cos(2.0d0* p)
      energy = eraw+ebsc

c      write(6,1) a2,b2,c2,eraw,energy
c1     format(' a2,b2,c2 :',5f20.10)

      return
      end
      subroutine spl_c2h3po_3(ro,alpha,phi,energy)
      implicit real*8 (a-h,o-z)

      data degrad / 57.29577951308232d 00 /
c
c     find location in spline grids
c
      call fnd_grd11(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1)
c
      call a_3(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,a3)
      call b_3(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,b3)
      call c_3(ro,alpha,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1,c3)
      call bsc(ro,ebsc)
c
      p = phi / degrad
      eraw = a3 + b3*cos(p) + c3*cos(2.0d0* p)
      energy = eraw+ebsc

c      write(6,1) a3,b3,c3,eraw,energy
c1     format(' a3,b3,c3 :',5f20.10)

      return
      end
      subroutine fnd_grd11(xi,yi,ix,iy,delxi,delyi,xix,xixp1,yiy,yiyp1)
      implicit real*8 (a-h,o-z)
      dimension delx(15),dely(18),x(16),y(19)
      data x( 1), x( 2) /  3.80000000d+00 ,  4.00000000d+00 /
      data x( 3), x( 4) /  4.20000000d+00 ,  4.40000000d+00 /
      data x( 5), x( 6) /  4.60000000d+00 ,  4.80000000d+00 /
      data x( 7), x( 8) /  5.00000000d+00 ,  5.20000000d+00 /
      data x( 9), x(10) /  5.60000000d+00 ,  6.00000000d+00 /
      data x(11), x(12) /  6.50000000d+00 ,  7.00000000d+00 /
      data x(13), x(14) /  8.00000000d+00 ,  1.00000000d+01 /
      data x(15), x(16) /  1.20000000d+01 ,  1.50000000d+01 /
 
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
      data delx( 3), delx( 4) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 5), delx( 6) /  2.00000000d-01 ,  2.00000000d-01 /
      data delx( 7), delx( 8) /  2.00000000d-01 ,  4.00000000d-01 /
      data delx( 9), delx(10) /  4.00000000d-01 ,  5.00000000d-01 /
      data delx(11), delx(12) /  5.00000000d-01 ,  1.00000000d+00 /
      data delx(13), delx(14) /  2.00000000d+00 ,  2.00000000d+00 /
      data delx(15) /            3.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(15), dely(16) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(17), dely(18) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  16 , 19 /

      iprint = 0

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
        call hunt(y,npty,yi,iy)
      endif
 25   yiy=y(iy)
      yiyp1=y(iy+1)
      delyi = dely(iy)
      if(iprint .gt. 2) then
        write(6,'(a,i3,a,2f10.5,a,1f10.5)') ' iy=',iy,
     x       '  yiy,yiyp1=',yiy,yiyp1,'  delyi=',delyi
      endif
c
      return
      end  
      subroutine bsc(xi,fi)
      implicit real*8 (a-h,o-z)
      dimension x( 17), f( 17),fpp( 17),del( 16)
 
      data fpp(  1), fpp(  2) / -5.9061089047d-02 , -3.2235321906d-02 /
      data fpp(  3), fpp(  4) /  4.8637376670d-02 , -5.9654184775d-02 /
      data fpp(  5), fpp(  6) /  7.1434362431d-02 , -5.3433264950d-02 /
      data fpp(  7), fpp(  8) /  1.0283697369d-02 , -4.3815245258d-03 /
      data fpp(  9), fpp( 10) / -8.0977510696d-04 , -1.3006250463d-03 /
      data fpp( 11), fpp( 12) / -4.3292974773d-04 ,  1.6674403724d-04 /
      data fpp( 13), fpp( 14) /  2.7632762160d-05 , -1.2703050990d-06 /
      data fpp( 15), fpp( 16) /  8.4984582355d-06 ,  7.5200928100d-07 /
      data fpp( 17) / -6.6975046405d-06 /
 
      data f(  1), f(  2) / -5.3762000000d-03 , -4.1474000000d-03 /
      data f(  3), f(  4) / -3.8477000000d-03 , -2.8636000000d-03 /
      data f(  5), f(  6) / -2.6698000000d-03 , -1.3250000000d-03 /
      data f(  7), f(  8) / -8.6030000000d-04 , -5.0680000000d-04 /
      data f(  9), f( 10) / -3.4800000000d-05 ,  1.9930000000d-04 /
      data f( 11), f( 12) /  2.5180000000d-04 ,  1.8490000000d-04 /
      data f( 13), f( 14) /  1.0300000000d-04 ,  4.9200000000d-05 /
      data f( 15), f( 16) /  1.6100000000d-05 ,  8.8000000000d-06 /
      data f( 17) /  0.0000000000d+00 /
 
      data x(  1), x(  2) /  3.8000000000d+00 ,  4.0000000000d+00 /
      data x(  3), x(  4) /  4.2000000000d+00 ,  4.4000000000d+00 /
      data x(  5), x(  6) /  4.6000000000d+00 ,  4.8000000000d+00 /
      data x(  7), x(  8) /  5.0000000000d+00 ,  5.2000000000d+00 /
      data x(  9), x( 10) /  5.6000000000d+00 ,  6.0000000000d+00 /
      data x( 11), x( 12) /  6.5000000000d+00 ,  7.0000000000d+00 /
      data x( 13), x( 14) /  8.0000000000d+00 ,  1.0000000000d+01 /
      data x( 15), x( 16) /  1.2000000000d+01 ,  1.5000000000d+01 /
      data x( 17) /  2.0000000000d+01 /
 
      data del(  1), del(  2) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del(  3), del(  4) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del(  5), del(  6) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del(  7), del(  8) /  2.0000000000d-01 ,  4.0000000000d-01 /
      data del(  9), del( 10) /  4.0000000000d-01 ,  5.0000000000d-01 /
      data del( 11), del( 12) /  5.0000000000d-01 ,  1.0000000000d+00 /
      data del( 13), del( 14) /  2.0000000000d+00 ,  2.0000000000d+00 /
      data del( 15), del( 16) /  3.0000000000d+00 ,  5.0000000000d+00 /
      data npts / 17/

      if(xi .le. x(1)) then
        ii=1
      else if( xi .ge. x(npts)) then
        ii=npts-1
      else
        call hunt(x,npts,xi,ii)
      endif
      
 20   fi = fpp(ii)   * (x(ii+1)-xi)**3 / (6.0*del(ii)) + 
     x     fpp(ii+1) * (xi-x(ii))**3   / (6.0*del(ii)) +
     x     ((f(ii+1)/del(ii))-(fpp(ii+1)*del(ii)/6.0)) * (xi-x(ii)) + 
     x     ((f(ii)  /del(ii))-(fpp(ii)  *del(ii)/6.0)) * (x(ii+1)-xi)
      fi = fi + 0.00 
      return
      end  
      subroutine a_1(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  1.42883200d-01 ,  1.04462600d-01 /
      data f( 1, 3),f( 1, 4) /  1.19813800d-01 ,  3.47116100d-01 /
      data f( 1, 5),f( 1, 6) /  5.23517800d-01 ,  2.30874800d-01 /
      data f( 1, 7),f( 1, 8) /  5.87290000d-02 ,  2.47935000d-02 /
      data f( 1, 9),f( 1,10) /  1.52898000d-02 ,  2.05270000d-03 /
      data f( 1,11),f( 1,12) / -1.12785000d-02 , -2.57864000d-02 /
      data f( 1,13),f( 1,14) / -3.67742000d-02 , -2.20902000d-02 /
      data f( 1,15),f( 1,16) /  1.14020000d-01 ,  4.37261800d-01 /
      data f( 1,17),f( 1,18) /  2.24625600d-01 ,  6.11170000d-03 /
      data f( 1,19) /          -2.49770000d-02 /
      data f( 2, 1),f( 2, 2) /  1.10585800d-01 ,  8.07888000d-02 /
      data f( 2, 3),f( 2, 4) /  7.67459000d-02 ,  1.81088000d-01 /
      data f( 2, 5),f( 2, 6) /  2.61219400d-01 ,  1.32339800d-01 /
      data f( 2, 7),f( 2, 8) /  4.28552000d-02 ,  2.38134000d-02 /
      data f( 2, 9),f( 2,10) /  1.45311000d-02 ,  2.75850000d-03 /
      data f( 2,11),f( 2,12) / -8.52450000d-03 , -2.10420000d-02 /
      data f( 2,13),f( 2,14) / -3.04733000d-02 , -2.27377000d-02 /
      data f( 2,15),f( 2,16) /  4.66821000d-02 ,  1.85054800d-01 /
      data f( 2,17),f( 2,18) /  9.63199000d-02 , -1.14220000d-02 /
      data f( 2,19) /          -2.91244000d-02 /
      data f( 3, 1),f( 3, 2) /  8.95927000d-02 ,  6.76304000d-02 /
      data f( 3, 3),f( 3, 4) /  5.50201000d-02 ,  9.97929000d-02 /
      data f( 3, 5),f( 3, 6) /  1.40338200d-01 ,  7.92567000d-02 /
      data f( 3, 7),f( 3, 8) /  3.09124000d-02 ,  1.92843000d-02 /
      data f( 3, 9),f( 3,10) /  1.06763000d-02 ,  1.58270000d-03 /
      data f( 3,11),f( 3,12) / -6.96040000d-03 , -1.69275000d-02 /
      data f( 3,13),f( 3,14) / -2.48683000d-02 , -2.18466000d-02 /
      data f( 3,15),f( 3,16) /  1.79065000d-02 ,  7.34006000d-02 /
      data f( 3,17),f( 3,18) /  3.64915000d-02 , -1.40367000d-02 /
      data f( 3,19) /          -2.28927000d-02 /
      data f( 4, 1),f( 4, 2) /  7.30549000d-02 ,  5.91388000d-02 /
      data f( 4, 3),f( 4, 4) /  4.44138000d-02 ,  6.10004000d-02 /
      data f( 4, 5),f( 4, 6) /  7.66055000d-02 ,  4.62141000d-02 /
      data f( 4, 7),f( 4, 8) /  2.32085000d-02 ,  1.47680000d-02 /
      data f( 4, 9),f( 4,10) /  7.54820000d-03 ,  5.45000000d-04 /
      data f( 4,11),f( 4,12) / -5.90350000d-03 , -1.38810000d-02 /
      data f( 4,13),f( 4,14) / -2.01824000d-02 , -1.96613000d-02 /
      data f( 4,15),f( 4,16) /  7.81500000d-04 ,  2.50785000d-02 /
      data f( 4,17),f( 4,18) /  1.22147000d-02 , -8.88160000d-03 /
      data f( 4,19) /          -1.16546000d-02 /
      data f( 5, 1),f( 5, 2) /  5.71195000d-02 ,  5.09738000d-02 /
      data f( 5, 3),f( 5, 4) /  3.80186000d-02 ,  4.12498000d-02 /
      data f( 5, 5),f( 5, 6) /  4.33311000d-02 ,  2.73174000d-02 /
      data f( 5, 7),f( 5, 8) /  1.60078000d-02 ,  1.09192000d-02 /
      data f( 5, 9),f( 5,10) /  5.11790000d-03 , -2.01800000d-04 /
      data f( 5,11),f( 5,12) / -5.02180000d-03 , -1.09013000d-02 /
      data f( 5,13),f( 5,14) / -1.61320000d-02 , -1.67155000d-02 /
      data f( 5,15),f( 5,16) / -6.25670000d-03 ,  5.78940000d-03 /
      data f( 5,17),f( 5,18) /  4.91040000d-03 , -1.30480000d-03 /
      data f( 5,19) /          -1.83680000d-03 /
      data f( 6, 1),f( 6, 2) /  4.28484000d-02 ,  4.06428000d-02 /
      data f( 6, 3),f( 6, 4) /  3.26828000d-02 ,  2.96215000d-02 /
      data f( 6, 5),f( 6, 6) /  2.54077000d-02 ,  1.69777000d-02 /
      data f( 6, 7),f( 6, 8) /  1.18589000d-02 ,  7.80870000d-03 /
      data f( 6, 9),f( 6,10) /  3.25480000d-03 , -7.28000000d-04 /
      data f( 6,11),f( 6,12) / -4.35870000d-03 , -8.80950000d-03 /
      data f( 6,13),f( 6,14) / -1.27660000d-02 , -1.35185000d-02 /
      data f( 6,15),f( 6,16) / -8.05940000d-03 , -7.16300000d-04 /
      data f( 6,17),f( 6,18) /  3.32130000d-03 ,  5.07750000d-03 /
      data f( 6,19) /           5.11150000d-03 /
      data f( 7, 1),f( 7, 2) /  3.07606000d-02 ,  2.99690000d-02 /
      data f( 7, 3),f( 7, 4) /  2.60791000d-02 ,  2.17143000d-02 /
      data f( 7, 5),f( 7, 6) /  1.59484000d-02 ,  1.13991000d-02 /
      data f( 7, 7),f( 7, 8) /  8.79000000d-03 ,  5.34190000d-03 /
      data f( 7, 9),f( 7,10) /  1.86110000d-03 , -1.07680000d-03 /
      data f( 7,11),f( 7,12) / -3.94050000d-03 , -7.15440000d-03 /
      data f( 7,13),f( 7,14) / -1.00402000d-02 , -1.08806000d-02 /
      data f( 7,15),f( 7,16) / -7.86850000d-03 , -2.39320000d-03 /
      data f( 7,17),f( 7,18) /  3.08700000d-03 ,  7.62150000d-03 /
      data f( 7,19) /           8.28940000d-03 /
      data f( 8, 1),f( 8, 2) /  2.16739000d-02 ,  2.12792000d-02 /
      data f( 8, 3),f( 8, 4) /  1.90955000d-02 ,  1.51555000d-02 /
      data f( 8, 5),f( 8, 6) /  1.08394000d-02 ,  7.97200000d-03 /
      data f( 8, 7),f( 8, 8) /  6.11620000d-03 ,  3.43880000d-03 /
      data f( 8, 9),f( 8,10) /  8.53000000d-04 , -1.28310000d-03 /
      data f( 8,11),f( 8,12) / -3.53880000d-03 , -5.85210000d-03 /
      data f( 8,13),f( 8,14) / -7.77600000d-03 , -8.48260000d-03 /
      data f( 8,15),f( 8,16) / -6.54400000d-03 , -2.01270000d-03 /
      data f( 8,17),f( 8,18) /  3.43580000d-03 ,  7.74140000d-03 /
      data f( 8,19) /           8.66860000d-03 /
      data f( 9, 1),f( 9, 2) /  9.84440000d-03 ,  9.74170000d-03 /
      data f( 9, 3),f( 9, 4) /  9.10300000d-03 ,  7.14040000d-03 /
      data f( 9, 5),f( 9, 6) /  5.21030000d-03 ,  3.77390000d-03 /
      data f( 9, 7),f( 9, 8) /  2.45560000d-03 ,  1.01740000d-03 /
      data f( 9, 9),f( 9,10) / -3.12200000d-04 , -1.43220000d-03 /
      data f( 9,11),f( 9,12) / -2.80360000d-03 , -4.00840000d-03 /
      data f( 9,13),f( 9,14) / -5.00520000d-03 , -5.25580000d-03 /
      data f( 9,15),f( 9,16) / -3.88050000d-03 , -1.05740000d-03 /
      data f( 9,17),f( 9,18) /  2.57420000d-03 ,  4.99650000d-03 /
      data f( 9,19) /           5.86410000d-03 /
      data f(10, 1),f(10, 2) /  3.76330000d-03 ,  3.75470000d-03 /
      data f(10, 3),f(10, 4) /  3.52070000d-03 ,  2.76170000d-03 /
      data f(10, 5),f(10, 6) /  1.97670000d-03 ,  1.22020000d-03 /
      data f(10, 7),f(10, 8) /  4.71900000d-04 , -2.21200000d-04 /
      data f(10, 9),f(10,10) / -7.48300000d-04 , -1.52200000d-03 /
      data f(10,11),f(10,12) / -2.17100000d-03 , -2.81580000d-03 /
      data f(10,13),f(10,14) / -3.37410000d-03 , -3.51120000d-03 /
      data f(10,15),f(10,16) / -2.73930000d-03 , -1.07970000d-03 /
      data f(10,17),f(10,18) /  8.07500000d-04 ,  2.20410000d-03 /
      data f(10,19) /           2.65560000d-03 /
      data f(11, 1),f(11, 2) /  5.01900000d-04 ,  5.37500000d-04 /
      data f(11, 3),f(11, 4) /  3.60700000d-04 ,  1.62300000d-04 /
      data f(11, 5),f(11, 6) / -5.95000000d-05 , -2.96800000d-04 /
      data f(11, 7),f(11, 8) / -5.51000000d-04 , -7.93900000d-04 /
      data f(11, 9),f(11,10) / -1.08330000d-03 , -1.32410000d-03 /
      data f(11,11),f(11,12) / -1.55420000d-03 , -1.86160000d-03 /
      data f(11,13),f(11,14) / -2.16830000d-03 , -2.30170000d-03 /
      data f(11,15),f(11,16) / -2.04720000d-03 , -1.39430000d-03 /
      data f(11,17),f(11,18) / -6.73700000d-04 , -1.02500000d-04 /
      data f(11,19) /           1.03500000d-04 /
      data f(12, 1),f(12, 2) / -6.54400000d-04 , -6.70800000d-04 /
      data f(12, 3),f(12, 4) / -6.94700000d-04 , -7.18300000d-04 /
      data f(12, 5),f(12, 6) / -7.30500000d-04 , -7.43900000d-04 /
      data f(12, 7),f(12, 8) / -8.27700000d-04 , -9.63900000d-04 /
      data f(12, 9),f(12,10) / -1.04050000d-03 , -1.06010000d-03 /
      data f(12,11),f(12,12) / -1.11850000d-03 , -1.26860000d-03 /
      data f(12,13),f(12,14) / -1.45630000d-03 , -1.58720000d-03 /
      data f(12,15),f(12,16) / -1.56550000d-03 , -1.35420000d-03 /
      data f(12,17),f(12,18) / -1.20210000d-03 , -1.15130000d-03 /
      data f(12,19) /          -1.13330000d-03 /
      data f(13, 1),f(13, 2) / -7.67200000d-04 , -7.73400000d-04 /
      data f(13, 3),f(13, 4) / -7.79000000d-04 , -8.19000000d-04 /
      data f(13, 5),f(13, 6) / -8.21500000d-04 , -8.13900000d-04 /
      data f(13, 7),f(13, 8) / -7.91200000d-04 , -7.41400000d-04 /
      data f(13, 9),f(13,10) / -6.80200000d-04 , -6.32100000d-04 /
      data f(13,11),f(13,12) / -6.23900000d-04 , -6.77700000d-04 /
      data f(13,13),f(13,14) / -7.67300000d-04 , -8.78500000d-04 /
      data f(13,15),f(13,16) / -9.96100000d-04 , -1.09460000d-03 /
      data f(13,17),f(13,18) / -1.14680000d-03 , -1.16210000d-03 /
      data f(13,19) /          -1.16500000d-03 /
      data f(14, 1),f(14, 2) / -3.46800000d-04 , -3.52300000d-04 /
      data f(14, 3),f(14, 4) / -3.61000000d-04 , -3.63300000d-04 /
      data f(14, 5),f(14, 6) / -3.54200000d-04 , -3.34200000d-04 /
      data f(14, 7),f(14, 8) / -3.06600000d-04 , -2.74500000d-04 /
      data f(14, 9),f(14,10) / -2.43800000d-04 , -2.23500000d-04 /
      data f(14,11),f(14,12) / -2.20200000d-04 , -2.35300000d-04 /
      data f(14,13),f(14,14) / -2.63700000d-04 , -3.01800000d-04 /
      data f(14,15),f(14,16) / -3.43300000d-04 , -3.80500000d-04 /
      data f(14,17),f(14,18) / -4.05900000d-04 , -4.19200000d-04 /
      data f(14,19) /          -4.23100000d-04 /
      data f(15, 1),f(15, 2) / -1.53900000d-04 , -1.53500000d-04 /
      data f(15, 3),f(15, 4) / -1.51500000d-04 , -1.47300000d-04 /
      data f(15, 5),f(15, 6) / -1.40700000d-04 , -1.31800000d-04 /
      data f(15, 7),f(15, 8) / -1.20100000d-04 , -1.06500000d-04 /
      data f(15, 9),f(15,10) / -9.25000000d-05 , -8.45000000d-05 /
      data f(15,11),f(15,12) / -8.49000000d-05 , -9.19000000d-05 /
      data f(15,13),f(15,14) / -1.03900000d-04 , -1.18100000d-04 /
      data f(15,15),f(15,16) / -1.32700000d-04 , -1.46500000d-04 /
      data f(15,17),f(15,18) / -1.58000000d-04 , -1.66100000d-04 /
      data f(15,19) /          -1.68900000d-04 /
      data f(16, 1),f(16, 2) / -4.23000000d-05 , -4.24000000d-05 /
      data f(16, 3),f(16, 4) / -4.20000000d-05 , -4.16000000d-05 /
      data f(16, 5),f(16, 6) / -4.06000000d-05 , -3.85000000d-05 /
      data f(16, 7),f(16, 8) / -3.47000000d-05 , -2.97000000d-05 /
      data f(16, 9),f(16,10) / -2.46000000d-05 , -2.21000000d-05 /
      data f(16,11),f(16,12) / -2.27000000d-05 , -2.52000000d-05 /
      data f(16,13),f(16,14) / -2.88000000d-05 , -3.29000000d-05 /
      data f(16,15),f(16,16) / -3.64000000d-05 , -3.95000000d-05 /
      data f(16,17),f(16,18) / -4.19000000d-05 , -4.30000000d-05 /
      data f(16,19) /          -4.40000000d-05 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 4.52147131d-01,-9.21023032d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 4.07750158d-01, 2.91614065d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 7.86722588d-01, 2.98087477d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.15765555d+00, 5.01952839d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 5.58513284d+00,-8.04272213d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.75132553d+00, 3.52625368d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 9.68698761d-02, 1.16753939d-03/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-1.74907051d-01, 9.62067379d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-1.69011353d-01,-8.64583461d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/-9.56616920d-02, 2.56226463d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-4.65507865d-02,-2.16782393d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/-6.89449587d-03,-9.51168931d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/-1.24369225d-02, 2.70930996d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/ 4.49301182d-02, 4.66095703d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.60959830d+00, 5.15025819d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.40435290d+00,-9.83923246d-03/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.52180234d+00, 2.05399166d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/ 5.50328642d-01, 1.27060381d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 3.96263369d-01, 4.10910509d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 2.85978238d-01,-4.99755738d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 2.65584684d-01, 1.19516477d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 5.38342324d-01, 1.56693583d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.15118391d+00, 1.15840195d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.64945933d+00,-3.48293861d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.17681145d+00, 1.27525426d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 8.56902479d-02, 7.45621590d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-9.44508977d-02,-3.11726167d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-8.53247931d-02,-3.53611236d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-5.07816160d-02, 2.31991110d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-3.02759269d-02,-2.80593203d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/-1.15485083d-02, 1.49681702d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/-1.61586550d-02, 1.53358640d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/ 3.77547635d-02, 4.01611271d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.01856091d+00, 1.94124828d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.59378919d+00,-4.02943037d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.73847782d+00, 5.50017216d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 3.75727716d-01, 6.88941510d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 2.54540762d-01, 2.09658675d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 9.95849170d-02,-3.52200941d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 1.07221108d-01, 2.43398810d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.61223116d-01, 8.15961416d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 9.47558821d-01, 1.54800453d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.02960985d+00,-1.68881323d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 3.59213682d-01, 5.02844465d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.50019132d-01, 4.41667368d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 2.03606421d-02,-6.65419369d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 4.58955260d-02, 5.70637957d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.65481559d-02, 1.45804186d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.08305057d-02,-3.09980539d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-4.13964711d-02, 2.39717970d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/-2.73134576d-02, 5.66888660d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 3.48408277d-02, 4.07022739d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 1.00503075d-01, 5.19104178d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.30341032d+00,-1.53897945d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 7.95881377d-01, 9.26216211d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.84610495d-01, 3.51346965d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 1.42438584d-01, 1.00232252d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-1.60229061d-02,-3.07252676d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 5.55088563d-03,-5.21716471d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 8.46902113d-02, 4.67405265d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 4.33970808d-01, 6.12465876d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 8.04376285d-01,-7.71281615d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 3.92408824d-01, 2.64089874d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-4.99317774d-02, 1.58070120d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 1.49283293d-02,-2.24643554d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 1.07476893d-02, 5.02930133d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/ 5.30399250d-03, 1.53431501d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-2.48205011d-03,-3.31199018d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.69343926d-02, 2.53964569d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.24525146d-02, 3.21000741d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 1.70119258d-02, 2.55553247d-04/
      data fpp( 4,15,1),fpp( 4,15,2)/ 3.27016795d-01, 1.40988939d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 6.92384513d-01,-5.88257004d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.10736670d-01,-1.76089219d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 5.13003034d-02, 1.64742692d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/-7.33350978d-02, 4.58036154d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 5.48667075d-02,-2.82621846d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-8.04346503d-02,-9.89593076d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 3.16810387d-02, 2.69889077d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 1.72842946d-01,-9.41299865d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.21630015d-01,-3.01231082d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 1.93036021d-01, 1.28637326d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.25187977d-01, 6.89277760d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 2.00510408d-02,-3.10884305d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.57837169d-02, 1.26639458d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/ 5.87087412d-03, 9.32864706d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/-5.52129380d-03,-1.99965341d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-3.63610994d-02, 7.08748934d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/-1.82014838d-02, 3.05745767d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 1.11864690d-02, 1.49446204d-04/
      data fpp( 5,15,1),fpp( 5,15,2)/ 1.04449747d-01, 3.41786084d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 2.82001626d-01,-1.90922637d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 1.07046943d-01,-4.59940592d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/-2.65567090d-02, 5.47268741d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/-6.21431928d-02, 1.68078563d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 4.62010760d-02,-1.57211827d-04/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-8.71228463d-03,-7.12703459d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-5.25043659d-02, 9.70292105d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 9.30024073d-02,-2.29244962d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 2.11753655d-01,-7.44812257d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.18997093d-01, 6.78773992d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 6.94986826d-03, 1.64362909d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 1.56125076d-02,-1.03359155d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.11974433d-02, 9.47803296d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 4.30251101d-03, 6.68978366d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-8.22277468d-03,-1.51111676d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/-4.67499486d-03, 4.54888670d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.74015500d-02, 2.65736208d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/-2.40778017d-02, 8.13966301d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/ 4.05092187d-02, 2.05358588d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/ 9.71189829d-02,-5.05000652d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.83555566d-02,-1.68655978d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.24248467d-01,-1.89215435d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/-1.08517131d-01,-1.07802283d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 8.78239886d-02,-5.53312065d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 6.38637888d-02,-3.47545869d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-1.18485749d-02, 8.45155429d-06/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.33124248d-02,-2.75456302d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 1.00970365d-01, 1.76649666d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 4.51406089d-02, 2.98817638d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 9.01254966d-03,-2.07800216d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.40539287d-02, 2.89832279d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 9.83651009d-03, 7.22473049d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/ 3.52908184d-03, 7.76755272d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.67739251d-03,-5.87975157d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.04439211d-02, 1.73025102d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-8.22231604d-03, 1.86447475d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.25973779d-03, 4.64147590d-05/
      data fpp( 7,15,1),fpp( 7,15,2)/ 3.25533786d-02, 2.68462165d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 5.38424421d-02,-6.00762504d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 2.27508300d-02,-2.52171633d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/-5.21944214d-02,-4.06475096d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/-6.93482833d-02,-6.68842452d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 5.26679695d-02,-1.86760207d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 5.08571295d-02,-1.69719586d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 4.29136655d-02,-2.07761449d-05/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 5.60078937d-02,-5.30146184d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 3.69098850d-02, 1.94159922d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 2.31654717d-02, 1.45594929d-05/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.62649331d-02,-1.69579638d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 1.27267777d-02, 3.97636233d-06/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 7.29651640d-03, 6.54851450d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/ 2.95616162d-03,-3.18842033d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-9.61795364d-04,-9.70833176d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-6.46932063d-03, 3.61575304d-06/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.89491858d-02, 9.87182103d-06/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.69461495d-02, 2.99349628d-05/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.82733081d-04, 2.91003276d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/-3.87875139d-03, 9.22572661d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/-2.18938767d-02,-1.09712341d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/-3.05888471d-02,-3.39147903d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/-3.38947357d-02,-5.60736049d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 3.59803471d-02, 1.55577645d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 3.45754671d-02,-3.43355290d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 2.62345409d-02,-1.99815648d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 1.66638565d-02, 3.92581226d-06/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.08689124d-02, 6.22831580d-06/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 7.53703032d-03, 7.82924537d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 9.96142589d-03,-2.27401395d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 6.72270269d-03, 1.11913126d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 5.10469577d-03, 4.31348889d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/-7.51775768d-04,-5.79708683d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-5.10810164d-04, 3.79085843d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-3.90382754d-03, 6.29653123d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/-4.95128458d-03, 6.17052908d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-8.63642043d-03, 1.94602305d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/-1.36847401d-02, 1.35425487d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/-7.99871689d-03, 1.32375746d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/-4.16378481d-03,-1.79828470d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/ 5.93750197d-03,-1.38641866d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 2.74959888d-03,-1.98424067d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 1.89756419d-02, 3.63779302d-07/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.89847520d-02,-1.49755860d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 1.75306709d-02,-7.89754488d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.37016801d-02, 1.58773814d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 9.44571538d-03,-1.34076594d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/ 8.35140700d-03, 1.75892502d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/ 6.77311336d-03,-1.98162348d-07/
      data fpp(10, 8,1),fpp(10, 8,2)/ 4.73741160d-03, 3.92875689d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/-3.74049480d-04,-5.55686522d-06/
      data fpp(10,10,1),fpp(10,10,2)/ 2.27469146d-03, 3.50270398d-06/
      data fpp(10,11,1),fpp(10,11,2)/-8.42463979d-04,-9.71950699d-07/
      data fpp(10,12,1),fpp(10,12,2)/-2.33161922d-03, 6.37098816d-07/
      data fpp(10,13,1),fpp(10,13,2)/-3.98442587d-03, 3.61355543d-06/
      data fpp(10,14,1),fpp(10,14,2)/-4.09066878d-03, 1.01806794d-05/
      data fpp(10,15,1),fpp(10,15,2)/-1.66455669d-03, 1.02037268d-05/
      data fpp(10,16,1),fpp(10,16,2)/-7.86381043d-04, 2.26641344d-06/
      data fpp(10,17,1),fpp(10,17,2)/ 4.60776596d-03,-5.61338054d-06/
      data fpp(10,18,1),fpp(10,18,2)/ 5.05758923d-03,-9.24889127d-06/
      data fpp(10,19,1),fpp(10,19,2)/ 7.74634023d-03,-1.40970544d-05/
      data fpp(11, 1,1),fpp(11, 1,2)/ 7.06281136d-03,-3.95028543d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/ 6.39171915d-03,-2.28742915d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 7.53095192d-03, 3.56002008d-07/
      data fpp(11, 4,1),fpp(11, 4,2)/ 6.31826633d-03,-4.32578888d-07/
      data fpp(11, 5,1),fpp(11, 5,2)/ 5.43949472d-03,-2.96864568d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/ 4.10831056d-03,-3.78675285d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 2.60905120d-03, 5.30387597d-07/
      data fpp(11, 8,1),fpp(11, 8,2)/ 9.80356099d-04,-1.06487510d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/ 2.30582151d-03, 9.39112811d-07/
      data fpp(11,10,1),fpp(11,10,2)/-1.43868625d-04, 2.24423858d-07/
      data fpp(11,11,1),fpp(11,11,2)/-7.33281546d-04,-1.19480824d-06/
      data fpp(11,12,1),fpp(11,12,2)/-1.36030877d-03,-8.31908824d-08/
      data fpp(11,13,1),fpp(11,13,2)/-1.68883920d-03, 1.56957177d-06/
      data fpp(11,14,1),fpp(11,14,2)/-1.67445605d-03, 4.20290379d-06/
      data fpp(11,15,1),fpp(11,15,2)/-6.85403864d-04, 4.89281307d-06/
      data fpp(11,16,1),fpp(11,16,2)/ 2.34854527d-03, 1.29843918d-07/
      data fpp(11,17,1),fpp(11,17,2)/ 4.19527040d-03,-1.35018875d-06/
      data fpp(11,18,1),fpp(11,18,2)/ 5.45627720d-03,-3.69308893d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 4.91809608d-03,-5.78945554d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 3.29551262d-03,-1.54871879d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/ 3.66197140d-03,-7.12562419d-08/
      data fpp(12, 3,1),fpp(12, 3,2)/ 2.85592139d-03,-1.01031535d-08/
      data fpp(12, 4,1),fpp(12, 4,2)/ 2.27645457d-03, 1.29668856d-07/
      data fpp(12, 5,1),fpp(12, 5,2)/ 1.56110573d-03, 1.75427730d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 8.92950749d-04,-9.03379777d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 6.99481854d-04,-7.85908621d-07/
      data fpp(12, 8,1),fpp(12, 8,2)/ 1.00596401d-03, 9.03014261d-07/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.17963440d-04, 7.49851576d-07/
      data fpp(12,10,1),fpp(12,10,2)/-1.12816954d-04,-4.82420566d-07/
      data fpp(12,11,1),fpp(12,11,2)/-5.70809839d-04,-1.14816931d-06/
      data fpp(12,12,1),fpp(12,12,2)/-8.95945693d-04,-4.26902180d-07/
      data fpp(12,13,1),fpp(12,13,2)/-1.11141732d-03, 5.99778034d-07/
      data fpp(12,14,1),fpp(12,14,2)/-1.09150703d-03, 1.43579005d-06/
      data fpp(12,15,1),fpp(12,15,2)/-6.43427854d-04, 2.81306178d-06/
      data fpp(12,16,1),fpp(12,16,2)/-9.50000295d-05,-1.31203718d-06/
      data fpp(12,17,1),fpp(12,17,2)/ 1.47835243d-03,-1.11691307d-06/
      data fpp(12,18,1),fpp(12,18,2)/ 3.30450195d-03,-2.98310552d-07/
      data fpp(12,19,1),fpp(12,19,2)/ 4.14847545d-03, 3.42155276d-07/
      data fpp(13, 1,1),fpp(13, 1,2)/-2.19143553d-04, 2.99804376d-07/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.97773779d-04, 1.18391248d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/-1.74240136d-04,-7.37369368d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-2.54968810d-05, 7.67086224d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/ 1.02935450d-04,-8.09755291d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/ 2.12192471d-04, 1.62815892d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/ 1.36428841d-04, 3.35711961d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.33070067d-04, 1.20336264d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.58601074d-04,-1.33057016d-07/
      data fpp(13,10,1),fpp(13,10,2)/-1.89614826d-04,-3.74108201d-07/
      data fpp(13,11,1),fpp(13,11,2)/-1.81729711d-04,-7.64510181d-07/
      data fpp(13,12,1),fpp(13,12,2)/-2.02608534d-04,-2.87851074d-07/
      data fpp(13,13,1),fpp(13,13,2)/-2.31328433d-04,-2.32085521d-07/
      data fpp(13,14,1),fpp(13,14,2)/-2.10050886d-04,-7.98068403d-08/
      data fpp(13,15,1),fpp(13,15,2)/-9.10145067d-05, 1.67312883d-07/
      data fpp(13,16,1),fpp(13,16,2)/ 1.87127454d-04, 5.56555310d-07/
      data fpp(13,17,1),fpp(13,17,2)/ 1.39907515d-04, 3.84465878d-07/
      data fpp(13,18,1),fpp(13,18,2)/-1.20844463d-04, 1.19581178d-07/
      data fpp(13,19,1),fpp(13,19,2)/-2.53074380d-04,-1.18790589d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/-2.13256521d-05,-1.26615197d-07/
      data fpp(14, 2,1),fpp(14, 2,2)/ 1.78563537d-06,-3.47696068d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-2.53402875d-05, 7.36936238d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-7.60866432d-05, 1.23995112d-07/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.15409216d-04, 1.14325930d-07/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.53502788d-04, 7.27011697d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/-1.41627449d-04, 5.08693917d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/-7.09218013d-05,-6.17873628d-09/
      data fpp(14, 9,1),fpp(14, 9,2)/-5.94784982d-05,-1.10154447d-07/
      data fpp(14,10,1),fpp(14,10,2)/-4.58470439d-05,-1.77203478d-07/
      data fpp(14,11,1),fpp(14,11,2)/-4.76559473d-05,-2.01031643d-07/
      data fpp(14,12,1),fpp(14,12,2)/-5.33015500d-05,-1.22669950d-07/
      data fpp(14,13,1),fpp(14,13,2)/-6.19060414d-05,-1.06288556d-07/
      data fpp(14,14,1),fpp(14,14,2)/-8.51438272d-05,-3.41758275d-08/
      data fpp(14,15,1),fpp(14,15,2)/-1.34242553d-04, 3.89918654d-08/
      data fpp(14,16,1),fpp(14,16,2)/-2.21532348d-04, 1.36208366d-07/
      data fpp(14,17,1),fpp(14,17,2)/-2.13448760d-04, 1.24174671d-07/
      data fpp(14,18,1),fpp(14,18,2)/-1.42967587d-04, 9.30929512d-08/
      data fpp(14,19,1),fpp(14,19,2)/-1.07064584d-04, 6.74535244d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/-3.68038382d-05, 1.01015076d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/-4.28187629d-05, 1.57969848d-08/
      data fpp(15, 3,1),fpp(15, 3,2)/-3.71487139d-05, 2.27105532d-08/
      data fpp(15, 4,1),fpp(15, 4,2)/-2.97065461d-05, 2.53608024d-08/
      data fpp(15, 5,1),fpp(15, 5,2)/-2.19985879d-05, 1.98462374d-08/
      data fpp(15, 6,1),fpp(15, 6,2)/-1.41313188d-05, 3.32542482d-08/
      data fpp(15, 7,1),fpp(15, 7,2)/-1.70690457d-05, 1.51367697d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-3.15927274d-05, 2.01986728d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/-3.11349332d-05,-7.19314611d-08/
      data fpp(15,10,1),fpp(15,10,2)/-3.13969981d-05,-9.24728285d-08/
      data fpp(15,11,1),fpp(15,11,2)/-3.02464998d-05,-6.21772247d-08/
      data fpp(15,12,1),fpp(15,12,2)/-3.26852655d-05,-5.48182725d-08/
      data fpp(15,13,1),fpp(15,13,2)/-3.67474020d-05,-1.85496854d-08/
      data fpp(15,14,1),fpp(15,14,2)/-3.88738054d-05,-2.98298584d-09/
      data fpp(15,15,1),fpp(15,15,2)/-3.53152817d-05, 6.48162880d-09/
      data fpp(15,16,1),fpp(15,16,2)/-2.11480609d-05, 2.50564706d-08/
      data fpp(15,17,1),fpp(15,17,2)/-2.56124767d-05, 3.12924887d-08/
      data fpp(15,18,1),fpp(15,18,2)/-4.19851897d-05, 5.37735747d-08/
      data fpp(15,19,1),fpp(15,19,2)/-5.02172827d-05, 7.16132127d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 1.83965619d-05, 9.52127574d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/ 1.68054529d-05, 5.95744852d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/ 4.22257122d-06,-3.35106983d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/ 4.21291592d-06, 7.44683082d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/ 3.50143682d-06, 9.56374657d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/ 9.23958795d-06, 2.02981829d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/ 2.17484514d-05, 1.12435219d-08/
      data fpp(16, 8,1),fpp(16, 8,2)/ 3.57902923d-05, 6.72772961d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/ 3.74021095d-05,-3.21544403d-08/
      data fpp(16,10,1),fpp(16,10,2)/ 3.78213562d-05,-3.41099683d-08/
      data fpp(16,11,1),fpp(16,11,2)/ 3.87589642d-05,-1.74056863d-08/
      data fpp(16,12,1),fpp(16,12,2)/ 4.55519185d-05,-1.02672865d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 5.40287010d-05,-7.52516778d-09/
      data fpp(16,14,1),fpp(16,14,2)/ 5.94419027d-05, 1.03679576d-08/
      data fpp(16,15,1),fpp(16,15,2)/ 6.08126408d-05, 2.05333732d-09/
      data fpp(16,16,1),fpp(16,16,2)/ 5.55151019d-05, 5.41869309d-09/
      data fpp(16,17,1),fpp(16,17,2)/ 5.71740955d-05, 1.82718903d-08/
      data fpp(16,18,1),fpp(16,18,2)/ 6.42290234d-05,-5.06254378d-10/
      data fpp(16,19,1),fpp(16,19,2)/ 6.78339985d-05,-1.02468728d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-1.09982278d-02, 3.32571508d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.30975316d-02, 2.38017902d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/-1.14912033d-01, 6.01692639d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 4.99607640d-03, 2.80715386d-03/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.35665672d-04,-5.02578718d-04/
      data fpppp( 1,11),fpppp( 1,12)/ 1.20323846d-04,-5.45993557d-04/
      data fpppp( 1,13),fpppp( 1,14)/-6.48272658d-04, 6.91365223d-03/
      data fpppp( 1,15),fpppp( 1,16)/ 6.34317320d-02,-1.26835394d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 4.32715353d-02, 8.41386521d-03/
      data fpppp( 1,19) /             3.21175093d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/-7.08727984d-03, 2.03044054d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.65545894d-02, 1.21562386d-02/
      data fpppp( 2, 5),fpppp( 2, 6)/-7.20535134d-02, 3.78024170d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 3.73544618d-03, 1.91460151d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-3.78171918d-05,-2.38308386d-04/
      data fpppp( 2,11),fpppp( 2,12)/ 1.48801447d-04,-4.63593623d-04/
      data fpppp( 2,13),fpppp( 2,14)/ 3.05319121d-04, 2.75373105d-03/
      data fpppp( 2,15),fpppp( 2,16)/ 4.42933202d-02,-8.42616833d-02/
      data fpppp( 2,17),fpppp( 2,18)/ 2.69210334d-02, 6.13122558d-03/
      data fpppp( 2,19) /             2.30478533d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.70952734d-03, 9.06716481d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 1.01287898d-02,-8.66580914d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-1.17226340d-02, 1.04095137d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-2.24332385d-03, 3.33594528d-03/
      data fpppp( 3, 9),fpppp( 3,10)/-1.78885481d-03, 5.26538709d-04/
      data fpppp( 3,11),fpppp( 3,12)/-1.99177521d-04, 7.89331513d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 2.56238365d-03,-7.44419144d-03/
      data fpppp( 3,15),fpppp( 3,16)/ 2.74248598d-02,-3.40205476d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 6.03115888d-03, 3.67139593d-03/
      data fpppp( 3,19) /             1.34291957d-02 /
      data fpppp( 4, 1),fpppp( 4, 2)/-1.46559086d-03, 4.06632980d-04/
      data fpppp( 4, 3),fpppp( 4, 4)/ 3.29299097d-03, 2.62987942d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.25450159d-02, 6.07807985d-04/
      data fpppp( 4, 7),fpppp( 4, 8)/ 8.29139547d-03,-3.34134737d-03/
      data fpppp( 4, 9),fpppp( 4,10)/ 9.31549206d-04,-4.60632864d-04/
      data fpppp( 4,11),fpppp( 4,12)/ 7.70441500d-04,-9.88984016d-04/
      data fpppp( 4,13),fpppp( 4,14)/ 2.57293560d-04, 3.49089064d-03/
      data fpppp( 4,15),fpppp( 4,16)/ 2.61156958d-03,-1.06153980d-02/
      data fpppp( 4,17),fpppp( 4,18)/ 1.02908881d-03, 1.83173135d-03/
      data fpppp( 4,19) /             5.73204369d-03 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 4.52427569d-03, 2.74137658d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/-6.44759208d-04, 1.58043336d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/-5.21946455d-03, 2.65456105d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/-1.75402259d-03, 2.12419574d-03/
      data fpppp( 5, 9),fpppp( 5,10)/-6.90583608d-04, 2.99407562d-04/
      data fpppp( 5,11),fpppp( 5,12)/-5.95806153d-04, 9.16958786d-04/
      data fpppp( 5,13),fpppp( 5,14)/-1.32063716d-04, 2.84996308d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 2.82459798d-03,-6.52607211d-03/
      data fpppp( 5,17),fpppp( 5,18)/ 2.12929674d-03, 4.89946991d-04/
      data fpppp( 5,19) /             1.79194542d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/-1.51463358d-03,-2.00645966d-04/
      data fpppp( 6, 3),fpppp( 6, 4)/ 2.98449420d-03,-3.79399568d-04/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.07222746d-03,-2.21591954d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/ 2.00342453d-03,-7.48947092d-04/
      data fpppp( 6, 9),fpppp( 6,10)/ 2.07701615d-04,-2.30651440d-04/
      data fpppp( 6,11),fpppp( 6,12)/ 3.77082940d-04,-3.13296390d-04/
      data fpppp( 6,13),fpppp( 6,14)/-1.00357479d-04, 1.07774452d-03/
      data fpppp( 6,15),fpppp( 6,16)/ 6.51757294d-05,-1.81708280d-03/
      data fpppp( 6,17),fpppp( 6,18)/-9.19235951d-04, 1.66359074d-03/
      data fpppp( 6,19) /             3.76499462d-03 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.99515116d-03,-6.14773136d-04/
      data fpppp( 7, 3),fpppp( 7, 4)/ 1.34911387d-03, 1.27071944d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.68217518d-03, 8.48719520d-04/
      data fpppp( 7, 7),fpppp( 7, 8)/ 4.69398911d-04,-2.56148869d-04/
      data fpppp( 7, 9),fpppp( 7,10)/-3.31290380d-07, 1.32073452d-04/
      data fpppp( 7,11),fpppp( 7,12)/-2.60618182d-04, 2.94221820d-04/
      data fpppp( 7,13),fpppp( 7,14)/-5.56939727d-05, 3.64180996d-04/
      data fpppp( 7,15),fpppp( 7,16)/-9.23347942d-05,-5.95116457d-04/
      data fpppp( 7,17),fpppp( 7,18)/-6.70039915d-04, 6.44057758d-04/
      data fpppp( 7,19) /             1.56129225d-03 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.88576981d-04,-1.50231078d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/ 5.21543855d-04,-6.73682812d-04/
      data fpppp( 8, 5),fpppp( 8, 6)/ 2.41653185d-04, 2.82857959d-05/
      data fpppp( 8, 7),fpppp( 8, 8)/ 5.58361105d-05,-4.98872453d-05/
      data fpppp( 8, 9),fpppp( 8,10)/ 3.01865210d-05,-5.46444990d-06/
      data fpppp( 8,11),fpppp( 8,12)/ 1.70151466d-05,-1.57970234d-04/
      data fpppp( 8,13),fpppp( 8,14)/ 1.96525396d-04, 2.40842738d-04/
      data fpppp( 8,15),fpppp( 8,16)/-3.04273542d-04,-1.91314651d-04/
      data fpppp( 8,17),fpppp( 8,18)/ 1.80385724d-04, 2.89810530d-05/
      data fpppp( 8,19) /             2.70349664d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.27140019d-04,-6.79271110d-05/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.73143106d-05, 6.33988678d-05/
      data fpppp( 9, 5),fpppp( 9, 6)/-9.73674925d-06, 1.23331853d-04/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.38214005d-04, 8.97370404d-05/
      data fpppp( 9, 9),fpppp( 9,10)/-1.23491180d-04, 1.49919800d-04/
      data fpppp( 9,11),fpppp( 9,12)/-1.10341793d-04, 7.34083946d-05/
      data fpppp( 9,13),fpppp( 9,14)/-4.25581656d-05,-6.14364603d-05/
      data fpppp( 9,15),fpppp( 9,16)/ 2.06512980d-04,-1.20554894d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 1.64641132d-04,-1.62028351d-04/
      data fpppp( 9,19) /            -3.13879118d-04 /
      data fpppp(10, 1),fpppp(10, 2)/-6.25522702d-06,-1.31509059d-05/
      data fpppp(10, 3),fpppp(10, 4)/-2.89326158d-05,-1.36132172d-05/
      data fpppp(10, 5),fpppp(10, 6)/ 5.77670491d-05,-2.77555974d-05/
      data fpppp(10, 7),fpppp(10, 8)/ 2.42162251d-05,-9.65537906d-05/
      data fpppp(10, 9),fpppp(10,10)/ 1.77453378d-04,-1.47647603d-04/
      data fpppp(10,11),fpppp(10,12)/ 6.71832495d-05,-2.34053841d-05/
      data fpppp(10,13),fpppp(10,14)/ 1.66192027d-05, 4.97223974d-05/
      data fpppp(10,15),fpppp(10,16)/-6.35674927d-05, 1.11671387d-04/
      data fpppp(10,17),fpppp(10,18)/-1.12159775d-04, 4.03082879d-05/
      data fpppp(10,19) /             8.52622864d-05 /
      data fpppp(11, 1),fpppp(11, 2)/ 5.64531938d-05, 2.46482282d-05/
      data fpppp(11, 3),fpppp(11, 4)/-4.64266081d-05, 1.99431030d-05/
      data fpppp(11, 5),fpppp(11, 6)/-1.33109646d-05, 6.15600197d-06/
      data fpppp(11, 7),fpppp(11, 8)/-2.13975558d-05, 7.16680773d-05/
      data fpppp(11, 9),fpppp(11,10)/-8.80251229d-05, 5.39230817d-05/
      data fpppp(11,11),fpppp(11,12)/-1.60505709d-05, 8.02234359d-06/
      data fpppp(11,13),fpppp(11,14)/ 1.87100420d-06, 5.06845471d-06/
      data fpppp(11,15),fpppp(11,16)/ 3.63353187d-05,-2.77159128d-05/
      data fpppp(11,17),fpppp(11,18)/ 3.29489268d-06,-2.06067581d-05/
      data fpppp(11,19) /            -2.88191357d-05 /
      data fpppp(12, 1),fpppp(12, 2)/-2.49300694d-05,-1.32969641d-05/
      data fpppp(12, 3),fpppp(12, 4)/ 7.76739893d-06,-4.17764030d-06/
      data fpppp(12, 5),fpppp(12, 6)/ 7.90240885d-07, 3.84830853d-06/
      data fpppp(12, 7),fpppp(12, 8)/ 1.22976901d-05,-2.30420060d-05/
      data fpppp(12, 9),fpppp(12,10)/ 1.42013710d-05,-6.33026749d-06/
      data fpppp(12,11),fpppp(12,12)/ 3.48694948d-06, 3.53891380d-07/
      data fpppp(12,13),fpppp(12,14)/ 1.67733855d-06, 7.05966969d-06/
      data fpppp(12,15),fpppp(12,16)/-4.22588427d-06, 1.58647863d-05/
      data fpppp(12,17),fpppp(12,18)/ 2.26221716d-06,-9.74583076d-06/
      data fpppp(12,19) /            -2.22094562d-05 /
      data fpppp(13, 1),fpppp(13, 2)/ 3.74670633d-06, 2.11058845d-06/
      data fpppp(13, 3),fpppp(13, 4)/-5.92280963d-08,-3.61099292d-07/
      data fpppp(13, 5),fpppp(13, 6)/ 2.84969832d-07,-1.92929866d-06/
      data fpppp(13, 7),fpppp(13, 8)/-3.66901428d-06, 4.98123914d-06/
      data fpppp(13, 9),fpppp(13,10)/-1.61786820d-06, 1.16126891d-06/
      data fpppp(13,11),fpppp(13,12)/-6.93275375d-07,-1.14003731d-07/
      data fpppp(13,13),fpppp(13,14)/ 6.78825814d-07, 3.98547157d-07/
      data fpppp(13,15),fpppp(13,16)/ 3.59251552d-06,-5.22227430d-06/
      data fpppp(13,17),fpppp(13,18)/-2.22513231d-06, 1.31088117d-06/
      data fpppp(13,19) /             4.69293137d-06 /
      data fpppp(14, 1),fpppp(14, 2)/-7.78742409d-07,-4.81967041d-07/
      data fpppp(14, 3),fpppp(14, 4)/-3.07622047d-07, 2.95229255d-07/
      data fpppp(14, 5),fpppp(14, 6)/-1.87867972d-07, 5.29982624d-07/
      data fpppp(14, 7),fpppp(14, 8)/ 1.06607219d-06,-1.26445290d-06/
      data fpppp(14, 9),fpppp(14,10)/ 4.35998743d-07,-3.48253000d-07/
      data fpppp(14,11),fpppp(14,12)/ 3.05917918d-08,-4.31612729d-09/
      data fpppp(14,13),fpppp(14,14)/-1.90860600d-07,-1.10239140d-07/
      data fpppp(14,15),fpppp(14,16)/-9.19839235d-07, 1.49813190d-06/
      data fpppp(14,17),fpppp(14,18)/ 6.49714706d-07,-3.53135667d-07/
      data fpppp(14,19) /            -1.31186227d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 2.12967720d-07, 1.22870880d-07/
      data fpppp(15, 3),fpppp(15, 4)/-3.35280941d-09,-3.13252275d-09/
      data fpppp(15, 5),fpppp(15, 6)/ 3.18303287d-08,-1.14630135d-07/
      data fpppp(15, 7),fpppp(15, 8)/-2.21609555d-07, 3.05911065d-07/
      data fpppp(15, 9),fpppp(15,10)/-1.03146147d-07, 6.34819776d-08/
      data fpppp(15,11),fpppp(15,12)/-6.60279725d-08,-1.47259288d-08/
      data fpppp(15,13),fpppp(15,14)/ 2.75294430d-08, 2.07521458d-08/
      data fpppp(15,15),fpppp(15,16)/ 2.30557596d-07,-3.06460705d-07/
      data fpppp(15,17),fpppp(15,18)/-1.22612974d-07, 8.24147719d-08/
      data fpppp(15,19) /             2.81391085d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-3.27834319d-07,-1.45384523d-07/
      data fpppp(16, 3),fpppp(16, 4)/ 2.49866051d-07,-9.96860993d-08/
      data fpppp(16, 5),fpppp(16, 6)/ 1.06768917d-07, 5.95882435d-08/
      data fpppp(16, 7),fpppp(16, 8)/ 6.11208479d-08,-2.12092990d-07/
      data fpppp(16, 9),fpppp(16,10)/ 4.14496911d-08,-2.52600007d-08/
      data fpppp(16,11),fpppp(16,12)/ 9.06919879d-08, 1.38128269d-08/
      data fpppp(16,13),fpppp(16,14)/-4.49136016d-08,-1.79732721d-08/
      data fpppp(16,15),fpppp(16,16)/-1.25741121d-07, 1.20841130d-07/
      data fpppp(16,17),fpppp(16,18)/ 5.97685577d-08,-3.61593042d-08/
      data fpppp(16,19) /            -1.22128512d-07 /
 

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
      subroutine b_1(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  2.55380000d-03 /
      data f( 1, 3),f( 1, 4) /  1.02180000d-03 ,  5.22100000d-04 /
      data f( 1, 5),f( 1, 6) / -4.23490000d-03 ,  7.66670000d-03 /
      data f( 1, 7),f( 1, 8) /  6.87840000d-03 ,  5.59710000d-03 /
      data f( 1, 9),f( 1,10) /  4.89200000d-04 , -6.40730000d-03 /
      data f( 1,11),f( 1,12) / -1.61423000d-02 , -2.84728000d-02 /
      data f( 1,13),f( 1,14) / -5.10470000d-02 , -1.15348400d-01 /
      data f( 1,15),f( 1,16) / -4.05872100d-01 , -1.04510770d+00 /
      data f( 1,17),f( 1,18) / -6.01674800d-01 , -1.20596400d-01 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  3.91590000d-03 /
      data f( 2, 3),f( 2, 4) /  1.39280000d-03 ,  1.34140000d-03 /
      data f( 2, 5),f( 2, 6) /  1.23200000d-03 ,  1.09534000d-02 /
      data f( 2, 7),f( 2, 8) /  8.66440000d-03 ,  6.12250000d-03 /
      data f( 2, 9),f( 2,10) /  2.48300000d-04 , -4.93840000d-03 /
      data f( 2,11),f( 2,12) / -1.26531000d-02 , -2.30508000d-02 /
      data f( 2,13),f( 2,14) / -4.12780000d-02 , -8.75666000d-02 /
      data f( 2,15),f( 2,16) / -2.44304900d-01 , -5.22319700d-01 /
      data f( 2,17),f( 2,18) / -3.29948900d-01 , -8.02873000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  5.35460000d-03 /
      data f( 3, 3),f( 3, 4) /  2.72030000d-03 ,  2.58670000d-03 /
      data f( 3, 5),f( 3, 6) /  1.55025000d-02 ,  1.21227000d-02 /
      data f( 3, 7),f( 3, 8) /  9.44450000d-03 ,  5.40330000d-03 /
      data f( 3, 9),f( 3,10) /  3.70600000d-04 , -3.74730000d-03 /
      data f( 3,11),f( 3,12) / -9.77200000d-03 , -1.82215000d-02 /
      data f( 3,13),f( 3,14) / -3.24195000d-02 , -6.49374000d-02 /
      data f( 3,15),f( 3,16) / -1.60641300d-01 , -2.74147500d-01 /
      data f( 3,17),f( 3,18) / -1.86924100d-01 , -5.59724000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  7.24930000d-03 /
      data f( 4, 3),f( 4, 4) /  4.37780000d-03 ,  4.58040000d-03 /
      data f( 4, 5),f( 4, 6) /  1.51287000d-02 ,  1.18915000d-02 /
      data f( 4, 7),f( 4, 8) /  1.20213000d-02 ,  4.21280000d-03 /
      data f( 4, 9),f( 4,10) /  2.85800000d-04 , -2.81420000d-03 /
      data f( 4,11),f( 4,12) / -7.44680000d-03 , -1.39008000d-02 /
      data f( 4,13),f( 4,14) / -2.52782000d-02 , -4.74183000d-02 /
      data f( 4,15),f( 4,16) / -1.02708100d-01 , -1.53153600d-01 /
      data f( 4,17),f( 4,18) / -1.14446000d-01 , -4.38047000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  9.52680000d-03 /
      data f( 5, 3),f( 5, 4) /  6.31690000d-03 ,  7.90660000d-03 /
      data f( 5, 5),f( 5, 6) /  1.42907000d-02 ,  1.12230000d-02 /
      data f( 5, 7),f( 5, 8) /  8.06700000d-03 ,  3.15880000d-03 /
      data f( 5, 9),f( 5,10) /  1.99000000d-04 , -2.10400000d-03 /
      data f( 5,11),f( 5,12) / -5.60660000d-03 , -1.07498000d-02 /
      data f( 5,13),f( 5,14) / -1.94911000d-02 , -3.43589000d-02 /
      data f( 5,15),f( 5,16) / -6.76284000d-02 , -9.28140000d-02 /
      data f( 5,17),f( 5,18) / -7.79207000d-02 , -3.79301000d-02 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  8.31490000d-03 /
      data f( 6, 3),f( 6, 4) /  9.42670000d-03 ,  9.99230000d-03 /
      data f( 6, 5),f( 6, 6) /  1.30944000d-02 ,  1.03831000d-02 /
      data f( 6, 7),f( 6, 8) /  6.44230000d-03 ,  2.31100000d-03 /
      data f( 6, 9),f( 6,10) /  1.21400000d-04 , -1.57190000d-03 /
      data f( 6,11),f( 6,12) / -4.06690000d-03 , -8.17060000d-03 /
      data f( 6,13),f( 6,14) / -1.43830000d-02 , -2.64030000d-02 /
      data f( 6,15),f( 6,16) / -4.56617000d-02 , -6.12518000d-02 /
      data f( 6,17),f( 6,18) / -5.68195000d-02 , -3.44445000d-02 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  5.63900000d-03 /
      data f( 7, 3),f( 7, 4) /  1.00865000d-02 ,  1.29548000d-02 /
      data f( 7, 5),f( 7, 6) /  1.17296000d-02 ,  8.89010000d-03 /
      data f( 7, 7),f( 7, 8) /  4.75340000d-03 ,  1.66420000d-03 /
      data f( 7, 9),f( 7,10) /  6.35000000d-05 , -1.17800000d-03 /
      data f( 7,11),f( 7,12) / -3.04400000d-03 , -6.13740000d-03 /
      data f( 7,13),f( 7,14) / -1.10095000d-02 , -1.97256000d-02 /
      data f( 7,15),f( 7,16) / -3.25713000d-02 , -4.34286000d-02 /
      data f( 7,17),f( 7,18) / -4.24163000d-02 , -2.85254000d-02 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  4.05680000d-03 /
      data f( 8, 3),f( 8, 4) /  7.86880000d-03 ,  1.04770000d-02 /
      data f( 8, 5),f( 8, 6) /  9.73010000d-03 ,  6.87700000d-03 /
      data f( 8, 7),f( 8, 8) /  3.54280000d-03 ,  1.19050000d-03 /
      data f( 8, 9),f( 8,10) /  2.25000000d-05 , -8.87100000d-04 /
      data f( 8,11),f( 8,12) / -2.25680000d-03 , -4.56280000d-03 /
      data f( 8,13),f( 8,14) / -8.64500000d-03 , -1.51169000d-02 /
      data f( 8,15),f( 8,16) / -2.41297000d-02 , -3.19130000d-02 /
      data f( 8,17),f( 8,18) / -3.32117000d-02 , -2.26668000d-02 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  2.10250000d-03 /
      data f( 9, 3),f( 9, 4) /  4.08210000d-03 ,  5.62330000d-03 /
      data f( 9, 5),f( 9, 6) /  5.30070000d-03 ,  3.64990000d-03 /
      data f( 9, 7),f( 9, 8) /  1.89490000d-03 ,  5.60200000d-04 /
      data f( 9, 9),f( 9,10) / -2.32000000d-05 , -6.49000000d-04 /
      data f( 9,11),f( 9,12) / -1.23270000d-03 , -2.45520000d-03 /
      data f( 9,13),f( 9,14) / -4.77450000d-03 , -8.68200000d-03 /
      data f( 9,15),f( 9,16) / -1.43061000d-02 , -1.94014000d-02 /
      data f( 9,17),f( 9,18) / -1.97048000d-02 , -1.30301000d-02 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.08530000d-03 /
      data f(10, 3),f(10, 4) /  2.23600000d-03 ,  2.81640000d-03 /
      data f(10, 5),f(10, 6) /  2.65080000d-03 ,  1.84600000d-03 /
      data f(10, 7),f(10, 8) /  1.00860000d-03 ,  3.77300000d-04 /
      data f(10, 9),f(10,10) / -1.14400000d-04 , -3.92800000d-04 /
      data f(10,11),f(10,12) / -6.74100000d-04 , -1.27980000d-03 /
      data f(10,13),f(10,14) / -2.50020000d-03 , -4.68190000d-03 /
      data f(10,15),f(10,16) / -7.90280000d-03 , -1.07805000d-02 /
      data f(10,17),f(10,18) / -1.08714000d-02 , -6.86390000d-03 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  5.41900000d-04 /
      data f(11, 3),f(11, 4) /  9.48700000d-04 ,  1.17990000d-03 /
      data f(11, 5),f(11, 6) /  1.11160000d-03 ,  8.01300000d-04 /
      data f(11, 7),f(11, 8) /  4.66100000d-04 ,  1.31500000d-04 /
      data f(11, 9),f(11,10) / -1.63600000d-04 , -2.18900000d-04 /
      data f(11,11),f(11,12) / -3.30800000d-04 , -5.46600000d-04 /
      data f(11,13),f(11,14) / -1.03240000d-03 , -1.98230000d-03 /
      data f(11,15),f(11,16) / -3.44240000d-03 , -4.80310000d-03 /
      data f(11,17),f(11,18) / -4.85180000d-03 , -3.15340000d-03 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  2.26500000d-04 /
      data f(12, 3),f(12, 4) /  4.04100000d-04 ,  5.02600000d-04 /
      data f(12, 5),f(12, 6) /  4.80300000d-04 ,  3.62800000d-04 /
      data f(12, 7),f(12, 8) /  1.10100000d-04 , -8.40000000d-05 /
      data f(12, 9),f(12,10) / -8.50000000d-05 , -1.21900000d-04 /
      data f(12,11),f(12,12) / -1.82400000d-04 , -2.38200000d-04 /
      data f(12,13),f(12,14) / -3.90100000d-04 , -7.53300000d-04 /
      data f(12,15),f(12,16) / -1.35930000d-03 , -2.02480000d-03 /
      data f(12,17),f(12,18) / -1.99780000d-03 , -1.18640000d-03 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  2.06000000d-05 /
      data f(13, 3),f(13, 4) /  2.97000000d-05 ,  1.42000000d-05 /
      data f(13, 5),f(13, 6) / -1.42000000d-05 , -2.60000000d-05 /
      data f(13, 7),f(13, 8) / -2.61000000d-05 , -2.44000000d-05 /
      data f(13, 9),f(13,10) / -2.59000000d-05 , -3.22000000d-05 /
      data f(13,11),f(13,12) / -3.26000000d-05 ,  3.50000000d-06 /
      data f(13,13),f(13,14) /  3.70000000d-05 ,  6.21000000d-05 /
      data f(13,15),f(13,16) /  7.72000000d-05 ,  7.26000000d-05 /
      data f(13,17),f(13,18) /  3.60000000d-05 ,  4.10000000d-06 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -1.90000000d-06 /
      data f(14, 3),f(14, 4) / -2.00000000d-06 , -8.00000000d-07 /
      data f(14, 5),f(14, 6) /  1.30000000d-06 ,  2.20000000d-06 /
      data f(14, 7),f(14, 8) /  1.30000000d-06 , -3.00000000d-07 /
      data f(14, 9),f(14,10) / -1.50000000d-06 ,  2.10000000d-06 /
      data f(14,11),f(14,12) /  1.63000000d-05 ,  4.06000000d-05 /
      data f(14,13),f(14,14) /  6.87000000d-05 ,  9.79000000d-05 /
      data f(14,15),f(14,16) /  1.21900000d-04 ,  1.28000000d-04 /
      data f(14,17),f(14,18) /  1.05800000d-04 ,  5.90000000d-05 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -8.00000000d-07 /
      data f(15, 3),f(15, 4) / -4.00000000d-07 ,  6.00000000d-07 /
      data f(15, 5),f(15, 6) /  1.90000000d-06 ,  3.10000000d-06 /
      data f(15, 7),f(15, 8) /  3.90000000d-06 ,  3.50000000d-06 /
      data f(15, 9),f(15,10) /  4.30000000d-06 ,  5.30000000d-06 /
      data f(15,11),f(15,12) /  1.13000000d-05 ,  2.04000000d-05 /
      data f(15,13),f(15,14) /  3.14000000d-05 ,  4.08000000d-05 /
      data f(15,15),f(15,16) /  4.56000000d-05 ,  4.34000000d-05 /
      data f(15,17),f(15,18) /  3.37000000d-05 ,  1.85000000d-05 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 ,  7.00000000d-07 /
      data f(16, 3),f(16, 4) /  1.30000000d-06 ,  2.40000000d-06 /
      data f(16, 5),f(16, 6) /  4.10000000d-06 ,  5.30000000d-06 /
      data f(16, 7),f(16, 8) /  5.60000000d-06 ,  5.20000000d-06 /
      data f(16, 9),f(16,10) /  4.00000000d-06 ,  4.20000000d-06 /
      data f(16,11),f(16,12) /  6.70000000d-06 ,  1.08000000d-05 /
      data f(16,13),f(16,14) /  1.49000000d-05 ,  1.86000000d-05 /
      data f(16,15),f(16,16) /  2.01000000d-05 ,  1.79000000d-05 /
      data f(16,17),f(16,18) /  1.31000000d-05 ,  6.00000000d-06 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-8.29501314d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-8.56676676d-03,-5.90357372d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 3.96148084d-02, 7.39450802d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 8.59910825d-04,-1.74806584d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 7.63362435d-01, 3.69843254d-04/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-7.12423827d-02,-3.05050434d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-7.80055631d-02, 8.89644827d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-5.04199309d-02,-8.03874966d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.24798493d-02, 2.98950382d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/-7.45232699d-03,-3.88865187d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-1.64543423d-02,-1.77534292d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/-1.57543857d-02,-4.58297646d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/-3.55728346d-03,-4.13549513d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/-1.31301616d-01,-8.03604185d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/-3.19850792d+00,-9.94537175d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/-1.04759738d+01, 1.96623772d-02/
      data fpp( 1,17,1),fpp( 1,17,2)/-4.65015630d+00,-3.74402693d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/-4.98874887d-01,-2.42753945d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-8.17473528d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-1.46126136d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 3.90853352d-03,-7.91317271d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 2.38328833d-02, 7.63130448d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 1.41101783d-02,-7.78184521d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.05940130d-01, 2.31480764d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-5.21652347d-02,-2.58256602d-04/
      data fpp( 2, 7,1),fpp( 2, 7,2)/-5.95613737d-02, 8.09216447d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-3.11701383d-02,-8.06039768d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.07953013d-02, 4.15562627d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-6.92034601d-03,-4.43710738d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.53088155d-02,-1.57519675d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/-1.71487286d-02,-5.36010562d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/-2.08379331d-02,-2.39613808d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/-1.25966768d-01,-6.71627713d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/-2.05441415d+00,-3.70085734d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/-7.01611235d+00, 8.19846707d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/-3.25998991d+00,-8.69874942d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-3.94165226d-01,-1.28151930d-03/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-4.16650585d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.90775168d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 4.42263268d-03,-6.79086638d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 8.52865854d-03,-1.69241767d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 6.59937578d-03, 2.85647371d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-6.66582954d-01,-3.42701306d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-3.77066787d-02, 1.07421853d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.65366058d-01,-4.48901050d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-1.15895160d-02,-9.64143271d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.11810547d-02, 2.39658359d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-6.53628896d-03,-3.13339107d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.35253959d-02,-1.30381930d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-4.55569984d-03,-6.20013172d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/-4.96659842d-02,-8.38665382d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.37721314d-01,-7.01726530d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/-2.69375474d-01,-9.00387341d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/-2.65194679d+00, 3.23513790d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.61504908d+00, 3.61175887d-06/
      data fpp( 3,18,1),fpp( 3,18,2)/-3.23594210d-01,-6.25886931d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-1.99882203d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-2.35703080d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 4.68009358d-02,-9.61158398d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-8.44751742d-03, 1.29184394d-05/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 7.17523185d-02, 2.28888082d-04/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.63746685d-01,-3.07728768d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-7.08305050d-03, 1.74896990d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-3.32397858d-01,-1.89839193d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 6.83320241d-03, 1.08161783d-04/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.86391741d-03,-9.91793892d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/-5.63449814d-03,-1.88700274d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.39746010d-02,-6.55795148d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/-4.09184720d-02,-6.41821667d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/-3.80781301d-02,-3.21173818d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/-8.96629783d-02,-4.53110306d-04/
      data fpp( 4,15,1),fpp( 4,15,2)/-7.27643954d-01,-1.44423393d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/-1.45284550d+00, 1.32146188d-03/
      data fpp( 4,17,1),fpp( 4,17,2)/-8.61818775d-01, 2.07761879d-04/
      data fpp( 4,18,1),fpp( 4,18,2)/-1.33537933d-01,-2.36487394d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-8.72008303d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-2.97845200d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.34206376d-01,-1.37136600d-04/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 6.75014111d-02, 8.21895994d-05/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-9.37336500d-02, 9.63542022d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-5.80337875d-02,-1.79942408d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 4.43880702d-04, 5.63074306d-05/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.84560375d-01,-5.05853141d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 4.73170640d-03, 4.09018257d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-5.74614957d-04, 3.88201111d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/-4.36071848d-03,-1.70218702d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/-3.32620005d-03,-7.77053037d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/-7.22541201d-03,-5.03320083d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/-1.15149545d-03,-6.78743628d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/-1.72581773d-01,-2.90108247d-04/
      data fpp( 5,15,1),fpp( 5,15,2)/-2.48073710d-01, 6.31184224d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/-6.34816217d-01, 5.22668557d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/-3.30595821d-01, 2.50941350d-04/
      data fpp( 5,18,1),fpp( 5,18,2)/-8.62190585d-02,-2.05959572d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-2.92187521d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-1.38805838d-04/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-3.33854327d-02,-7.16193246d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-8.59531271d-02,-6.90286401d-06/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.17107281d-01, 6.64587806d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.46434647d-02,-1.06742258d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-2.04024723d-02, 1.17062533d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-5.64036423d-02,-1.38527548d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 5.16997200d-03, 3.22747658d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 8.14542415d-04, 1.25569171d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-3.63762796d-03,-7.51953259d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-1.77955988d-02,-1.92795613d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/-1.59498799d-02,-1.18842221d-05/
      data fpp( 6,13,1),fpp( 6,13,2)/-5.91658881d-02,-5.97055504d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/ 1.44650719d-02,-9.77495764d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/-2.47011207d-01, 1.63818561d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/-3.24499632d-01, 2.52338152d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/-1.29412942d-01, 1.75609536d-04/
      data fpp( 6,18,1),fpp( 6,18,2)/ 1.20064167d-01, 1.21785704d-04/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 6.14176480d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-6.51818415d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 4.81481066d-02,-1.49546317d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-9.11889027d-02,-5.15328907d-06/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-2.43175475d-01,-5.91842120d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-2.58150713d-02,-3.71986279d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.67989915d-02,-2.27943368d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 3.14241941d-02, 1.70652099d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 4.73840561d-03, 1.73834970d-05/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 2.71445297d-04, 2.71080197d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/-1.81876969d-03,-6.67470490d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-3.01140479d-03,-1.34819824d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.08750683d-02,-1.30413657d-05/
      data fpp( 7,13,1),fpp( 7,13,2)/-2.23749520d-02,-4.10745549d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/-7.70535141d-02,-5.33004146d-05/
      data fpp( 7,15,1),fpp( 7,15,2)/-9.53264615d-02, 6.50021330d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/-1.28035253d-01, 1.46603561d-04/
      data fpp( 7,17,1),fpp( 7,17,2)/-1.56452412d-01, 1.19261541d-04/
      data fpp( 7,18,1),fpp( 7,18,2)/-2.90126082d-02, 1.49066274d-04/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 1.62543363d-04/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 8.00343301d-06/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 4.84800629d-03,-4.17086603d-06/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 1.90837378d-02,-6.00796890d-06/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 3.95496198d-02,-4.40252584d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-6.58817966d-03,-1.91969976d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 9.58343813d-03,-5.55875105d-06/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 2.45186591d-03, 1.25660019d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 1.84140557d-03, 1.42087436d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 6.34676399d-04, 1.65702366d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/-4.53729327d-03,-5.33283827d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-5.51378206d-03,-7.93167057d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-9.33984685d-03,-1.91184795d-05/
      data fpp( 8,13,1),fpp( 8,13,2)/-2.68430370d-03,-2.21664116d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.65560156d-02,-3.55978742d-05/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.90029467d-02, 1.21039083d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/-1.09499355d-01, 6.09522411d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/-2.45674101d-02, 1.33163127d-04/
      data fpp( 8,18,1),fpp( 8,18,2)/-1.30887340d-02, 1.17011249d-04/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 1.26105875d-04/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 2.57741789d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 6.76067784d-03,-2.53183577d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.26694881d-02, 1.75925202d-07/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 6.76012812d-03,-2.44758650d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.65320746d-02,-1.41004651d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 9.61543134d-03, 1.18572526d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 5.93105523d-03, 3.10556403d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 3.99783049d-03, 1.16100186d-05/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-6.78501845d-04,-4.46763852d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/ 1.63251464d-03, 3.71653544d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-2.58920144d-03,-7.87250325d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-5.60292529d-03,-1.05545224d-05/
      data fpp( 9,13,1),fpp( 9,13,2)/-1.29533629d-02,-1.57174070d-05/
      data fpp( 9,14,1),fpp( 9,14,2)/-1.61489463d-02,-2.18678495d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/-1.00629290d-02, 1.92804973d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/-1.96930730d-03, 5.28246296d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/-3.19078135d-02, 7.60226766d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/-2.42462440d-02, 6.17706638d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 5.82186681d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 6.83466227d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/ 3.25053237d-03, 1.00667547d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 3.01081000d-03,-6.93736414d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.01648677d-02,-7.47521891d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 7.19113125d-03,-7.92176023d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 5.32483652d-03, 8.10259811d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.38391318d-03, 2.72472098d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/-1.05522752d-03, 6.56856262d-07/
      data fpp(10, 9,1),fpp(10, 9,2)/ 3.73080980d-04, 3.02385397d-06/
      data fpp(10,10,1),fpp(10,10,2)/-1.31401530d-03, 4.57278579d-08/
      data fpp(10,11,1),fpp(10,11,2)/-1.58566219d-03,-3.38076540d-06/
      data fpp(10,12,1),fpp(10,12,2)/-3.20595198d-03,-5.98666625d-06/
      data fpp(10,13,1),fpp(10,13,2)/-5.35974480d-03,-9.55456959d-06/
      data fpp(10,14,1),fpp(10,14,2)/-1.01531992d-02,-1.34730554d-05/
      data fpp(10,15,1),fpp(10,15,2)/-1.90065871d-02, 1.09479108d-06/
      data fpp(10,16,1),fpp(10,16,2)/-2.85246654d-02, 2.96858911d-05/
      data fpp(10,17,1),fpp(10,17,2)/-2.30575857d-02, 4.73696447d-05/
      data fpp(10,18,1),fpp(10,18,2)/-2.00700401d-02, 2.67395301d-05/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 1.70562350d-05/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-8.97407394d-07/
      data fpp(11, 2,1),fpp(11, 2,2)/ 3.63941193d-04,-1.44818521d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 3.51329354d-03,-1.41585176d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 2.92937384d-03,-3.42440775d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 3.44246783d-03,-2.85651725d-06/
      data fpp(11, 6,1),fpp(11, 6,2)/ 2.18244347d-03, 3.30476740d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 2.42068386d-04, 4.06102877d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.88354695d-04,-4.56917891d-07/
      data fpp(11, 9,1),fpp(11, 9,2)/ 7.54909949d-04, 4.15706128d-06/
      data fpp(11,10,1),fpp(11,10,2)/-8.79566199d-05,-1.78332721d-06/
      data fpp(11,11,1),fpp(11,11,2)/-7.39054965d-04,-4.19752426d-07/
      data fpp(11,12,1),fpp(11,12,2)/-1.64143263d-03,-2.77166308d-06/
      data fpp(11,13,1),fpp(11,13,2)/-3.34402842d-03,-4.69359523d-06/
      data fpp(11,14,1),fpp(11,14,2)/-5.74192574d-03,-6.29995598d-06/
      data fpp(11,15,1),fpp(11,15,2)/-8.57534311d-03,-7.18580857d-07/
      data fpp(11,16,1),fpp(11,16,2)/-1.09051587d-02, 1.51382794d-05/
      data fpp(11,17,1),fpp(11,17,2)/-1.19980407d-02, 1.88854632d-05/
      data fpp(11,18,1),fpp(11,18,2)/-4.28486056d-03, 1.41458676d-05/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 1.18310662d-05/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00,-1.71057715d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/ 7.65702857d-04,-5.20884570d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/ 7.60815836d-04,-6.79404005d-07/
      data fpp(12, 4,1),fpp(12, 4,2)/ 1.13843695d-03,-1.50749941d-06/
      data fpp(12, 5,1),fpp(12, 5,2)/ 8.28597435d-04,-5.38598356d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.94189594d-04,-2.05010717d-06/
      data fpp(12, 7,1),fpp(12, 7,2)/ 1.12381328d-03, 6.27027027d-07/
      data fpp(12, 8,1),fpp(12, 8,2)/ 1.02900874d-03, 3.05799906d-06/
      data fpp(12, 9,1),fpp(12, 9,2)/-3.25520775d-04,-1.27302326d-06/
      data fpp(12,10,1),fpp(12,10,2)/-1.79758217d-04,-1.19906018d-07/
      data fpp(12,11,1),fpp(12,11,2)/-1.35717951d-04, 3.36647332d-07/
      data fpp(12,12,1),fpp(12,12,2)/-4.23517490d-04,-9.44683310d-07/
      data fpp(12,13,1),fpp(12,13,2)/-1.07614151d-03,-2.32391409d-06/
      data fpp(12,14,1),fpp(12,14,2)/-2.17349781d-03,-2.43766032d-06/
      data fpp(12,15,1),fpp(12,15,2)/-3.74724043d-03,-2.49344461d-06/
      data fpp(12,16,1),fpp(12,16,2)/-4.63309984d-03, 8.84143876d-06/
      data fpp(12,17,1),fpp(12,17,2)/-4.92465138d-03, 8.67768956d-06/
      data fpp(12,18,1),fpp(12,18,2)/-4.63451769d-03, 3.51180298d-06/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-2.24901491d-07/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 8.67215799d-09/
      data fpp(13, 2,1),fpp(13, 2,2)/ 7.03208327d-05,-1.00344316d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 2.49705721d-04,-2.97294894d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/ 3.17202240d-04,-1.86476108d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/ 4.01573779d-04, 2.69199325d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/ 3.55409482d-04, 1.05678806d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/-3.76740367d-05, 1.00854487d-08/
      data fpp(13, 8,1),fpp(13, 8,2)/-2.37603578d-04,-3.80206014d-08/
      data fpp(13, 9,1),fpp(13, 9,2)/ 1.05073498d-05,-5.00030432d-08/
      data fpp(13,10,1),fpp(13,10,2)/-4.25470404d-05,-4.99672257d-08/
      data fpp(13,11,1),fpp(13,11,2)/-1.05318664d-04, 6.03871946d-07/
      data fpp(13,12,1),fpp(13,12,2)/-1.59331214d-04,-1.75520559d-07/
      data fpp(13,13,1),fpp(13,13,2)/-2.44561250d-04,-5.77897115d-08/
      data fpp(13,14,1),fpp(13,14,2)/-4.64143711d-04,-9.73205954d-08/
      data fpp(13,15,1),fpp(13,15,2)/-8.48807163d-04,-1.52927907d-07/
      data fpp(13,16,1),fpp(13,16,2)/-1.40332114d-03,-4.72967778d-07/
      data fpp(13,17,1),fpp(13,17,2)/-1.27222550d-03, 1.24799017d-07/
      data fpp(13,18,1),fpp(13,18,2)/-4.15016642d-04, 2.55771709d-07/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 5.20114145d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00, 2.27401357d-08/
      data fpp(14, 2,1),fpp(14, 2,2)/-9.86392654d-06, 1.85197285d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-5.38750798d-05, 1.11809502d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-7.81251932d-05, 1.47564707d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.12270055d-04,-1.62068332d-08/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.04623244d-04,-2.19291381d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/ 8.15469580d-07,-4.07661435d-09/
      data fpp(14, 8,1),fpp(14, 8,2)/ 5.56563635d-05,-3.76440446d-09/
      data fpp(14, 9,1),fpp(14, 9,2)/-9.46166212d-06, 4.31342322d-08/
      data fpp(14,10,1),fpp(14,10,2)/-1.29770518d-07, 1.19227476d-07/
      data fpp(14,11,1),fpp(14,11,2)/ 7.76496891d-06, 1.15955865d-07/
      data fpp(14,12,1),fpp(14,12,2)/ 2.03023872d-05, 2.29490635d-08/
      data fpp(14,13,1),fpp(14,13,2)/ 3.80045061d-05, 2.02478808d-08/
      data fpp(14,14,1),fpp(14,14,2)/ 8.66800372d-05,-3.79405869d-08/
      data fpp(14,15,1),fpp(14,15,2)/ 1.77591702d-04,-1.80485533d-07/
      data fpp(14,16,1),fpp(14,16,2)/ 3.17413338d-04,-3.14117280d-07/
      data fpp(14,17,1),fpp(14,17,2)/ 2.82302186d-04,-2.61045348d-07/
      data fpp(14,18,1),fpp(14,18,2)/ 7.31587733d-05,-1.17701329d-07/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00,-1.49335395d-10/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00, 1.78476097d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/ 4.53487347d-06, 1.23047806d-08/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.57445986d-05, 4.93326803d-09/
      data fpp(15, 4,1),fpp(15, 4,2)/ 1.98985329d-05, 3.96214733d-09/
      data fpp(15, 5,1),fpp(15, 5,2)/ 2.51564416d-05,-2.78185733d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/ 2.21334944d-05, 1.16528199d-09/
      data fpp(15, 7,1),fpp(15, 7,2)/-2.78784158d-06,-2.58792706d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-1.54718754d-05, 3.03518005d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/-5.60701350d-07,-2.35279314d-08/
      data fpp(15,10,1),fpp(15,10,2)/-3.58387753d-06, 7.57599252d-08/
      data fpp(15,11,1),fpp(15,11,2)/-6.59121117d-06, 2.04882307d-08/
      data fpp(15,12,1),fpp(15,12,2)/-7.82833481d-06, 2.82871521d-08/
      data fpp(15,13,1),fpp(15,13,2)/-1.09567745d-05,-1.96368390d-08/
      data fpp(15,14,1),fpp(15,14,2)/-2.19264373d-05,-4.57397959d-08/
      data fpp(15,15,1),fpp(15,15,2)/-4.30596441d-05,-7.34039773d-08/
      data fpp(15,16,1),fpp(15,16,2)/-7.63322141d-05,-8.06442947d-08/
      data fpp(15,17,1),fpp(15,17,2)/-6.98332455d-05,-5.40188437d-08/
      data fpp(15,18,1),fpp(15,18,2)/-2.07184509d-05,-3.32803304d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00,-1.08598348d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00,-6.93996668d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/-8.64029388d-06,-1.12006664d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/-1.70319422d-05, 5.42023325d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/-1.44449807d-05, 9.43913363d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/-8.14143506d-06,-7.17676777d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/-3.46281865d-06,-1.07320626d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/ 7.28249222d-06,-3.89498200d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/ 1.18020091d-05,-1.56880095d-08/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.17677925d-06, 1.86470198d-08/
      data fpp(16,10,1),fpp(16,10,2)/ 8.09943876d-06, 2.50999301d-08/
      data fpp(16,11,1),fpp(16,11,2)/ 1.87273913d-05, 1.89532597d-08/
      data fpp(16,12,1),fpp(16,12,2)/ 2.63595245d-05,-4.91296888d-09/
      data fpp(16,13,1),fpp(16,13,2)/ 3.74862444d-05, 6.98615818d-10/
      data fpp(16,14,1),fpp(16,14,2)/ 5.76014329d-05,-2.18814944d-08/
      data fpp(16,15,1),fpp(16,15,2)/ 8.44376792d-05,-4.51726382d-08/
      data fpp(16,16,1),fpp(16,16,2)/ 1.10431821d-04,-1.94279526d-08/
      data fpp(16,17,1),fpp(16,17,2)/ 1.02942694d-04,-3.31155512d-08/
      data fpp(16,18,1),fpp(16,18,2)/ 5.24556540d-05, 1.38901575d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00, 4.35549213d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 1.04049075d-03, 2.49516504d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/-7.61625041d-03, 2.27536482d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/-3.53228972d-02, 2.27115001d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/-5.85260507d-03, 2.75984891d-03/
      data fpppp( 1, 9),fpppp( 1,10)/-2.46794169d-03, 9.42000467d-04/
      data fpppp( 1,11),fpppp( 1,12)/-4.42505144d-05,-1.82880099d-04/
      data fpppp( 1,13),fpppp( 1,14)/ 1.46559965d-03,-1.40760046d-02/
      data fpppp( 1,15),fpppp( 1,16)/-1.21529300d-01, 2.47577628d-01/
      data fpppp( 1,17),fpppp( 1,18)/-8.25842073d-02,-1.77129663d-02/
      data fpppp( 1,19) /            -6.57083187d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 2.48834190d-04, 8.96063526d-04/
      data fpppp( 2, 3),fpppp( 2, 4)/-2.87213932d-03, 8.81367047d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/-1.42893832d-02, 9.34774340d-03/
      data fpppp( 2, 7),fpppp( 2, 8)/-2.05903689d-03, 1.03564664d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-1.26909742d-03, 4.59877831d-04/
      data fpppp( 2,11),fpppp( 2,12)/-1.07832272d-05,-2.38315443d-05/
      data fpppp( 2,13),fpppp( 2,14)/-4.84807405d-06,-6.04315396d-03/
      data fpppp( 2,15),fpppp( 2,16)/-8.52216490d-02, 1.64934701d-01/
      data fpppp( 2,17),fpppp( 2,18)/-5.14479173d-02,-1.25608975d-02/
      data fpppp( 2,19) /            -4.66080600d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 7.43218502d-04,-1.38156117d-03/
      data fpppp( 3, 3),fpppp( 3, 4)/ 4.76402978d-03,-1.80366765d-02/
      data fpppp( 3, 5),fpppp( 3, 6)/ 2.71074933d-02,-1.22697804d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-3.57658411d-03, 3.77441817d-03/
      data fpppp( 3, 9),fpppp( 3,10)/-8.79246463d-04,-3.25405890d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 1.94230340d-04, 1.83860875d-04/
      data fpppp( 3,13),fpppp( 3,14)/-4.17447267d-03, 1.39373271d-02/
      data fpppp( 3,15),fpppp( 3,16)/-5.41907656d-02, 6.77707059d-02/
      data fpppp( 3,17),fpppp( 3,18)/-1.17239166d-02,-5.60160997d-03/
      data fpppp( 3,19) /            -2.39412831d-02 /
      data fpppp( 4, 1),fpppp( 4, 2)/-3.26354025d-03,-1.28435473d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 2.27799584d-03, 2.99268711d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-2.76739884d-03,-4.99911951d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/ 1.34947726d-02,-9.10721859d-03/
      data fpppp( 4, 9),fpppp( 4,10)/ 2.34208108d-03,-5.32853562d-04/
      data fpppp( 4,11),fpppp( 4,12)/-2.01168073d-04, 2.21299767d-04/
      data fpppp( 4,13),fpppp( 4,14)/ 1.10302179d-03,-7.89889831d-03/
      data fpppp( 4,15),fpppp( 4,16)/-4.69119619d-03, 2.14304489d-02/
      data fpppp( 4,17),fpppp( 4,18)/-2.05690352d-03,-4.96758775d-03/
      data fpppp( 4,19) /            -1.37573200d-02 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 9.72881129d-03, 4.59694251d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/-7.96173157d-03, 5.47341289d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.11582458d-03, 4.35655378d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/-7.77206095d-03, 4.89498023d-03/
      data fpppp( 5, 9),fpppp( 5,10)/-1.33651914d-03, 5.42309392d-04/
      data fpppp( 5,11),fpppp( 5,12)/-5.43481113d-04, 1.33559124d-03/
      data fpppp( 5,13),fpppp( 5,14)/-4.20049612d-03, 4.81614156d-03/
      data fpppp( 5,15),fpppp( 5,16)/-9.30776962d-03, 1.37399026d-02/
      data fpppp( 5,17),fpppp( 5,18)/-4.19406664d-03,-5.54254095d-04/
      data fpppp( 5,19) /            -3.07837921d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.47448256d-03,-1.12271002d-03/
      data fpppp( 6, 3),fpppp( 6, 4)/ 5.81438693d-03,-6.79715151d-03/
      data fpppp( 6, 5),fpppp( 6, 6)/ 3.04276563d-03,-1.32883821d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/ 2.21527325d-03,-1.67776774d-03/
      data fpppp( 6, 9),fpppp( 6,10)/ 5.40055060d-04,-4.88256950d-04/
      data fpppp( 6,11),fpppp( 6,12)/ 8.30624713d-04,-1.87402052d-03/
      data fpppp( 6,13),fpppp( 6,14)/ 3.96175375d-03,-6.96217637d-03/
      data fpppp( 6,15),fpppp( 6,16)/ 3.78051738d-03, 2.87937807d-03/
      data fpppp( 6,17),fpppp( 6,18)/ 1.05647730d-03,-3.84186222d-03/
      data fpppp( 6,19) /            -7.86150495d-03 /
      data fpppp( 7, 1),fpppp( 7, 2)/-3.84363736d-03,-1.43398980d-03/
      data fpppp( 7, 3),fpppp( 7, 4)/-1.66951037d-03, 7.35305750d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/-5.58190103d-03, 2.47388716d-03/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.96122125d-03, 8.76459398d-04/
      data fpppp( 7, 9),fpppp( 7,10)/-2.11486650d-04, 1.12091923d-04/
      data fpppp( 7,11),fpppp( 7,12)/-1.83026249d-04, 2.19751368d-04/
      data fpppp( 7,13),fpppp( 7,14)/-9.14152436d-04, 8.46137677d-04/
      data fpppp( 7,15),fpppp( 7,16)/-2.86061396d-04,-5.68042740d-04/
      data fpppp( 7,17),fpppp( 7,18)/ 2.81573031d-03,-1.34346072d-03/
      data fpppp( 7,19) /            -3.34751916d-03 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.84126547d-04,-2.34698481d-05/
      data fpppp( 8, 3),fpppp( 8, 4)/ 4.73016356d-04,-1.49478654d-03/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.50990890d-03,-8.06284032d-04/
      data fpppp( 8, 7),fpppp( 8, 8)/ 3.17035824d-04,-7.05925527d-05/
      data fpppp( 8, 9),fpppp( 8,10)/-7.04417437d-05, 1.14445098d-04/
      data fpppp( 8,11),fpppp( 8,12)/-1.35609795d-04, 2.57019522d-04/
      data fpppp( 8,13),fpppp( 8,14)/-2.63571817d-04,-4.34367552d-04/
      data fpppp( 8,15),fpppp( 8,16)/-3.13471133d-04, 2.40528344d-03/
      data fpppp( 8,17),fpppp( 8,18)/-1.78196137d-03, 3.15365916d-04/
      data fpppp( 8,19) /             6.17101178d-04 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.68421835d-05, 4.00856559d-05/
      data fpppp( 9, 3),fpppp( 9, 4)/-2.88296864d-04, 4.04011589d-04/
      data fpppp( 9, 5),fpppp( 9, 6)/-3.86871106d-04, 1.42157449d-04/
      data fpppp( 9, 7),fpppp( 9, 8)/ 1.21773390d-05,-8.57977230d-05/
      data fpppp( 9, 9),fpppp( 9,10)/ 1.66427098d-04,-1.60669738d-04/
      data fpppp( 9,11),fpppp( 9,12)/ 8.42879001d-05,-1.04002329d-04/
      data fpppp( 9,13),fpppp( 9,14)/ 7.15185913d-05, 6.72192128d-05/
      data fpppp( 9,15),fpppp( 9,16)/ 2.16500600d-04,-8.12765344d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 7.52633097d-04, 5.82375047d-05/
      data fpppp( 9,19) /             9.49735050d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-1.28373050d-04,-6.56422386d-05/
      data fpppp(10, 3),fpppp(10, 4)/ 1.81526720d-04,-2.16837837d-04/
      data fpppp(10, 5),fpppp(10, 6)/ 7.81569797d-05,-2.93435800d-05/
      data fpppp(10, 7),fpppp(10, 8)/-2.52603763d-05, 1.00492044d-04/
      data fpppp(10, 9),fpppp(10,10)/-8.46608465d-05, 5.12270551d-05/
      data fpppp(10,11),fpppp(10,12)/-3.53204101d-05, 9.13601106d-06/
      data fpppp(10,13),fpppp(10,14)/-3.32338157d-05,-3.45804449d-05/
      data fpppp(10,15),fpppp(10,16)/-7.20404133d-05, 2.82860675d-04/
      data fpppp(10,17),fpppp(10,18)/-1.60292807d-04, 2.09538506d-04/
      data fpppp(10,19) /             3.47088450d-04 /
      data fpppp(11, 1),fpppp(11, 2)/ 8.73822290d-05, 3.91715406d-05/
      data fpppp(11, 3),fpppp(11, 4)/-7.69437220d-05, 4.46070248d-05/
      data fpppp(11, 5),fpppp(11, 6)/-3.56635560d-05,-8.33990154d-06/
      data fpppp(11, 7),fpppp(11, 8)/ 2.82021183d-05, 8.73111203d-06/
      data fpppp(11, 9),fpppp(11,10)/-2.59104298d-05, 1.03452977d-05/
      data fpppp(11,11),fpppp(11,12)/-3.96466771d-06,-9.56338622d-06/
      data fpppp(11,13),fpppp(11,14)/-5.79487473d-06,-8.97520666d-06/
      data fpppp(11,15),fpppp(11,16)/ 1.55644984d-05,-2.30666796d-05/
      data fpppp(11,17),fpppp(11,18)/ 1.50918232d-04,-5.22425131d-05/
      data fpppp(11,19) /            -1.47647356d-04 /
      data fpppp(12, 1),fpppp(12, 2)/-1.81538023d-05,-9.87205197d-06/
      data fpppp(12, 3),fpppp(12, 4)/ 1.14066175d-05,-1.28039302d-05/
      data fpppp(12, 5),fpppp(12, 6)/-1.43853392d-06, 1.70839661d-05/
      data fpppp(12, 7),fpppp(12, 8)/-9.05543853d-06,-2.43279054d-05/
      data fpppp(12, 9),fpppp(12,10)/ 3.07835614d-05,-8.78881568d-06/
      data fpppp(12,11),fpppp(12,12)/-1.73163626d-06,-4.19502755d-06/
      data fpppp(12,13),fpppp(12,14)/-3.37772259d-06,-8.97801827d-06/
      data fpppp(12,15),fpppp(12,16)/ 1.07066160d-05, 7.42454698d-06/
      data fpppp(12,17),fpppp(12,18)/-4.74633154d-06, 4.64618926d-05/
      data fpppp(12,19) /             7.95618017d-05 /
      data fpppp(13, 1),fpppp(13, 2)/ 3.13247533d-06, 1.42601947d-06/
      data fpppp(13, 3),fpppp(13, 4)/-2.29270989d-06, 1.03151798d-06/
      data fpppp(13, 5),fpppp(13, 6)/-8.20860874d-07,-5.58022464d-06/
      data fpppp(13, 7),fpppp(13, 8)/ 2.32660609d-06, 7.86303894d-06/
      data fpppp(13, 9),fpppp(13,10)/-6.89633366d-06, 1.65237660d-06/
      data fpppp(13,11),fpppp(13,12)/-2.96206779d-07, 5.79949795d-08/
      data fpppp(13,13),fpppp(13,14)/-1.80882231d-06,-8.83851295d-07/
      data fpppp(13,15),fpppp(13,16)/-4.56063189d-06, 8.93534733d-06/
      data fpppp(13,17),fpppp(13,18)/ 9.95581962d-06,-5.19183285d-06/
      data fpppp(13,19) /            -1.57200211d-05 /
      data fpppp(14, 1),fpppp(14, 2)/-8.38612511d-07,-4.25357110d-07/
      data fpppp(14, 3),fpppp(14, 4)/ 4.91207347d-07,-3.53809890d-07/
      data fpppp(14, 5),fpppp(14, 6)/ 3.30347302d-07, 1.53992106d-06/
      data fpppp(14, 7),fpppp(14, 8)/-6.22517353d-07,-2.08572084d-06/
      data fpppp(14, 9),fpppp(14,10)/ 1.76786554d-06,-5.18746294d-07/
      data fpppp(14,11),fpppp(14,12)/ 2.20890503d-07,-8.62549841d-08/
      data fpppp(14,13),fpppp(14,14)/ 4.34011468d-07, 2.08613846d-07/
      data fpppp(14,15),fpppp(14,16)/ 1.26570115d-06,-2.33682012d-06/
      data fpppp(14,17),fpppp(14,18)/-2.41438801d-06, 1.55243651d-06/
      data fpppp(14,19) /             4.36372034d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 1.93533959d-07, 8.77904825d-08/
      data fpppp(15, 3),fpppp(15, 4)/-1.44204789d-07, 6.56812226d-08/
      data fpppp(15, 5),fpppp(15, 6)/-5.22816344d-08,-3.53406034d-07/
      data fpppp(15, 7),fpppp(15, 8)/ 1.52002435d-07, 4.79634425d-07/
      data fpppp(15, 9),fpppp(15,10)/-4.14827659d-07, 1.03615193d-07/
      data fpppp(15,11),fpppp(15,12)/ 1.31743767d-09,-2.67234327d-09/
      data fpppp(15,13),fpppp(15,14)/-1.04107031d-07,-5.13729171d-08/
      data fpppp(15,15),fpppp(15,16)/-3.00213941d-07, 5.23866889d-07/
      data fpppp(15,17),fpppp(15,18)/ 5.91038702d-07,-3.31072135d-07/
      data fpppp(15,19) /            -9.70530787d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-9.76707145d-08,-1.17984855d-08/
      data fpppp(16, 3),fpppp(16, 4)/ 1.59783392d-07, 3.13814991d-08/
      data fpppp(16, 5),fpppp(16, 6)/-6.23143363d-08, 1.20380092d-07/
      data fpppp(16, 7),fpppp(16, 8)/-5.52043648d-08,-2.73110269d-07/
      data fpppp(16, 9),fpppp(16,10)/ 2.98960632d-07, 1.01411057d-08/
      data fpppp(16,11),fpppp(16,12)/-5.72074735d-08, 3.89396311d-08/
      data fpppp(16,13),fpppp(16,14)/ 1.11124147d-07, 5.58719027d-08/
      data fpppp(16,15),fpppp(16,16)/ 6.86517055d-08,-3.81004971d-07/
      data fpppp(16,17),fpppp(16,18)/-5.53627978d-07, 1.56421040d-08/
      data fpppp(16,19) /             3.72942732d-07 /
 

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
      subroutine c_1(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -2.46857000d-02 /
      data f( 1, 3),f( 1, 4) /  2.19303000d-02 ,  2.81869900d-01 /
      data f( 1, 5),f( 1, 6) /  4.84253400d-01 ,  2.04763400d-01 /
      data f( 1, 7),f( 1, 8) /  3.90093000d-02 ,  1.25088000d-02 /
      data f( 1, 9),f( 1,10) /  6.24860000d-03 , -3.56330000d-03 /
      data f( 1,11),f( 1,12) / -9.28140000d-03 , -1.24940000d-02 /
      data f( 1,13),f( 1,14) / -1.14764000d-02 ,  1.19359000d-02 /
      data f( 1,15),f( 1,16) /  1.52292800d-01 ,  4.76470800d-01 /
      data f( 1,17),f( 1,18) /  2.61414300d-01 ,  3.61979000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.98574000d-02 /
      data f( 2, 3),f( 2, 4) / -1.43700000d-03 ,  1.26049700d-01 /
      data f( 2, 5),f( 2, 6) /  2.24840100d-01 ,  1.06601500d-01 /
      data f( 2, 7),f( 2, 8) /  2.13707000d-02 ,  5.90730000d-03 /
      data f( 2, 9),f( 2,10) /  1.14490000d-03 , -5.32460000d-03 /
      data f( 2,11),f( 2,12) / -9.81020000d-03 , -1.34986000d-02 /
      data f( 2,13),f( 2,14) / -1.33695000d-02 ,  1.79090000d-03 /
      data f( 2,15),f( 2,16) /  7.58179000d-02 ,  2.16679100d-01 /
      data f( 2,17),f( 2,18) /  1.28549400d-01 ,  1.88759000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.45369000d-02 /
      data f( 3, 3),f( 3, 4) / -1.01752000d-02 ,  5.19879000d-02 /
      data f( 3, 5),f( 3, 6) /  1.06780100d-01 ,  5.48364000d-02 /
      data f( 3, 7),f( 3, 8) /  1.16545000d-02 ,  5.12250000d-03 /
      data f( 3, 9),f( 3,10) /  8.95200000d-04 , -4.28540000d-03 /
      data f( 3,11),f( 3,12) / -8.51270000d-03 , -1.28508000d-02 /
      data f( 3,13),f( 3,14) / -1.39872000d-02 , -5.07140000d-03 /
      data f( 3,15),f( 3,16) /  3.84365000d-02 ,  9.57951000d-02 /
      data f( 3,17),f( 3,18) /  5.93919000d-02 ,  8.56850000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -8.43790000d-03 /
      data f( 4, 3),f( 4, 4) / -1.02851000d-02 ,  1.99550000d-02 /
      data f( 4, 5),f( 4, 6) /  4.68230000d-02 ,  2.46769000d-02 /
      data f( 4, 7),f( 4, 8) /  8.11240000d-03 ,  4.66480000d-03 /
      data f( 4, 9),f( 4,10) /  1.08530000d-03 , -3.05810000d-03 /
      data f( 4,11),f( 4,12) / -6.69090000d-03 , -1.12235000d-02 /
      data f( 4,13),f( 4,14) / -1.33313000d-02 , -8.82730000d-03 /
      data f( 4,15),f( 4,16) /  1.44215000d-02 ,  3.90889000d-02 /
      data f( 4,17),f( 4,18) /  2.53238000d-02 ,  3.25440000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -2.01770000d-03 /
      data f( 5, 3),f( 5, 4) / -5.94230000d-03 ,  7.98820000d-03 /
      data f( 5, 5),f( 5, 6) /  1.92998000d-02 ,  1.06932000d-02 /
      data f( 5, 7),f( 5, 8) /  5.12020000d-03 ,  4.16680000d-03 /
      data f( 5, 9),f( 5,10) /  1.18710000d-03 , -2.06600000d-03 /
      data f( 5,11),f( 5,12) / -5.02590000d-03 , -8.79430000d-03 /
      data f( 5,13),f( 5,14) / -1.15570000d-02 , -9.85560000d-03 /
      data f( 5,15),f( 5,16) /  1.98210000d-03 ,  1.34786000d-02 /
      data f( 5,17),f( 5,18) /  9.98690000d-03 ,  1.40080000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  5.05000000d-04 /
      data f( 6, 3),f( 6, 4) / -3.93500000d-04 ,  4.49060000d-03 /
      data f( 6, 5),f( 6, 6) /  7.56130000d-03 ,  5.07610000d-03 /
      data f( 6, 7),f( 6, 8) /  4.50200000d-03 ,  3.62780000d-03 /
      data f( 6, 9),f( 6,10) /  1.17460000d-03 , -1.33050000d-03 /
      data f( 6,11),f( 6,12) / -3.71140000d-03 , -6.84620000d-03 /
      data f( 6,13),f( 6,14) / -9.36810000d-03 , -8.96820000d-03 /
      data f( 6,15),f( 6,16) / -3.23840000d-03 ,  2.62990000d-03 /
      data f( 6,17),f( 6,18) /  3.10580000d-03 ,  1.29280000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  9.41000000d-04 /
      data f( 7, 3),f( 7, 4) /  2.10430000d-03 ,  3.66050000d-03 /
      data f( 7, 5),f( 7, 6) /  3.39410000d-03 ,  3.35050000d-03 /
      data f( 7, 7),f( 7, 8) /  4.14860000d-03 ,  3.04080000d-03 /
      data f( 7, 9),f( 7,10) /  1.08000000d-03 , -8.13400000d-04 /
      data f( 7,11),f( 7,12) / -2.82260000d-03 , -5.18430000d-03 /
      data f( 7,13),f( 7,14) / -7.21790000d-03 , -7.54000000d-03 /
      data f( 7,15),f( 7,16) / -4.89520000d-03 , -1.36660000d-03 /
      data f( 7,17),f( 7,18) /  3.69500000d-04 ,  9.88500000d-04 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  8.49800000d-04 /
      data f( 8, 3),f( 8, 4) /  2.30500000d-03 ,  2.68570000d-03 /
      data f( 8, 5),f( 8, 6) /  2.40440000d-03 ,  2.83880000d-03 /
      data f( 8, 7),f( 8, 8) /  3.45720000d-03 ,  2.45720000d-03 /
      data f( 8, 9),f( 8,10) /  9.43500000d-04 , -4.64800000d-04 /
      data f( 8,11),f( 8,12) / -2.12490000d-03 , -3.85730000d-03 /
      data f( 8,13),f( 8,14) / -5.24010000d-03 , -5.73170000d-03 /
      data f( 8,15),f( 8,16) / -4.40850000d-03 , -1.88270000d-03 /
      data f( 8,17),f( 8,18) /  3.65000000d-05 ,  6.59300000d-04 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  5.74700000d-04 /
      data f( 9, 3),f( 9, 4) /  1.68630000d-03 ,  1.89460000d-03 /
      data f( 9, 5),f( 9, 6) /  1.99300000d-03 ,  2.22050000d-03 /
      data f( 9, 7),f( 9, 8) /  2.12550000d-03 ,  1.49000000d-03 /
      data f( 9, 9),f( 9,10) /  6.49000000d-04 , -1.36300000d-04 /
      data f( 9,11),f( 9,12) / -1.19910000d-03 , -2.08070000d-03 /
      data f( 9,13),f( 9,14) / -2.78780000d-03 , -2.96810000d-03 /
      data f( 9,15),f( 9,16) / -2.10710000d-03 , -7.58700000d-04 /
      data f( 9,17),f( 9,18) /  4.09100000d-04 ,  2.59300000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  3.01400000d-04 /
      data f(10, 3),f(10, 4) /  8.75600000d-04 ,  1.12170000d-03 /
      data f(10, 5),f(10, 6) /  1.28730000d-03 ,  1.31020000d-03 /
      data f(10, 7),f(10, 8) /  1.11370000d-03 ,  7.60900000d-04 /
      data f(10, 9),f(10,10) /  4.23500000d-04 , -2.07400000d-04 /
      data f(10,11),f(10,12) / -6.83800000d-04 , -1.11380000d-03 /
      data f(10,13),f(10,14) / -1.46220000d-03 , -1.49810000d-03 /
      data f(10,15),f(10,16) / -9.50200000d-04 , -1.07100000d-04 /
      data f(10,17),f(10,18) /  3.50100000d-04 ,  2.26300000d-04 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  1.31300000d-04 /
      data f(11, 3),f(11, 4) /  2.17000000d-04 ,  3.76900000d-04 /
      data f(11, 5),f(11, 6) /  4.84000000d-04 ,  4.89000000d-04 /
      data f(11, 7),f(11, 8) /  3.99700000d-04 ,  2.38100000d-04 /
      data f(11, 9),f(11,10) / -1.01000000d-05 , -2.02500000d-04 /
      data f(11,11),f(11,12) / -3.44100000d-04 , -5.16700000d-04 /
      data f(11,13),f(11,14) / -6.67000000d-04 , -6.68900000d-04 /
      data f(11,15),f(11,16) / -3.94400000d-04 ,  2.73000000d-05 /
      data f(11,17),f(11,18) /  2.05200000d-04 ,  1.24900000d-04 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 , -2.40000000d-06 /
      data f(12, 3),f(12, 4) /  1.58000000d-05 ,  4.66000000d-05 /
      data f(12, 5),f(12, 6) /  8.32000000d-05 ,  1.20800000d-04 /
      data f(12, 7),f(12, 8) /  5.60000000d-05 , -8.23000000d-05 /
      data f(12, 9),f(12,10) / -1.59500000d-04 , -1.60700000d-04 /
      data f(12,11),f(12,12) / -1.72200000d-04 , -2.40700000d-04 /
      data f(12,13),f(12,14) / -3.16700000d-04 , -3.21500000d-04 /
      data f(12,15),f(12,16) / -1.95800000d-04 ,  4.72000000d-05 /
      data f(12,17),f(12,18) /  1.23200000d-04 ,  4.44000000d-05 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 , -3.01000000d-05 /
      data f(13, 3),f(13, 4) / -9.45000000d-05 , -1.26300000d-04 /
      data f(13, 5),f(13, 6) / -1.62300000d-04 , -1.91700000d-04 /
      data f(13, 7),f(13, 8) / -2.02700000d-04 , -1.81900000d-04 /
      data f(13, 9),f(13,10) / -1.39000000d-04 , -9.19000000d-05 /
      data f(13,11),f(13,12) / -6.50000000d-05 , -8.08000000d-05 /
      data f(13,13),f(13,14) / -1.11900000d-04 , -1.42200000d-04 /
      data f(13,15),f(13,16) / -1.58600000d-04 , -1.47100000d-04 /
      data f(13,17),f(13,18) / -9.52000000d-05 , -2.96000000d-05 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -1.60000000d-05 /
      data f(14, 3),f(14, 4) / -4.92000000d-05 , -7.82000000d-05 /
      data f(14, 5),f(14, 6) / -9.17000000d-05 , -9.03000000d-05 /
      data f(14, 7),f(14, 8) / -8.04000000d-05 , -6.54000000d-05 /
      data f(14, 9),f(14,10) / -4.92000000d-05 , -3.36000000d-05 /
      data f(14,11),f(14,12) / -2.50000000d-05 , -2.59000000d-05 /
      data f(14,13),f(14,14) / -3.25000000d-05 , -4.14000000d-05 /
      data f(14,15),f(14,16) / -4.77000000d-05 , -4.45000000d-05 /
      data f(14,17),f(14,18) / -2.88000000d-05 , -9.20000000d-06 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -4.80000000d-06 /
      data f(15, 3),f(15, 4) / -1.54000000d-05 , -2.62000000d-05 /
      data f(15, 5),f(15, 6) / -3.36000000d-05 , -3.67000000d-05 /
      data f(15, 7),f(15, 8) / -3.56000000d-05 , -3.15000000d-05 /
      data f(15, 9),f(15,10) / -2.51000000d-05 , -1.96000000d-05 /
      data f(15,11),f(15,12) / -1.64000000d-05 , -1.65000000d-05 /
      data f(15,13),f(15,14) / -1.86000000d-05 , -2.01000000d-05 /
      data f(15,15),f(15,16) / -1.94000000d-05 , -1.55000000d-05 /
      data f(15,17),f(15,18) / -9.00000000d-06 , -2.70000000d-06 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 , -1.70000000d-06 /
      data f(16, 3),f(16, 4) / -5.80000000d-06 , -1.11000000d-05 /
      data f(16, 5),f(16, 6) / -1.57000000d-05 , -1.83000000d-05 /
      data f(16, 7),f(16, 8) / -1.79000000d-05 , -1.49000000d-05 /
      data f(16, 9),f(16,10) / -1.11000000d-05 , -7.80000000d-06 /
      data f(16,11),f(16,12) / -6.50000000d-06 , -7.40000000d-06 /
      data f(16,13),f(16,14) / -9.60000000d-06 , -1.02000000d-05 /
      data f(16,15),f(16,16) / -9.00000000d-06 , -6.50000000d-06 /
      data f(16,17),f(16,18) / -3.70000000d-06 , -5.00000000d-07 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-5.87631011d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 4.42521685d-03, 4.73875023d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 5.14281368d-01, 2.97023292d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.02127376d+00, 4.44609296d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 5.55736608d+00,-8.20203610d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.75875501d+00, 3.45112512d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 2.48169664d-01, 1.22168963d-03/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.77452683d-01, 1.73323787d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.27813194d-01,-7.66011401d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.33161298d-01, 7.59701819d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/ 7.76978654d-02, 1.83484127d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 5.71971415d-02, 9.66167302d-07/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.19261548d-02, 2.31598918d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/ 8.66337311d-02, 4.16320160d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.59499097d+00, 5.11979644d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.30189168d+00,-9.86623992d-03/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.29672912d+00, 1.99109325d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/ 2.25399940d-01, 1.29227293d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 4.18092504d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-2.58793827d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 1.37495663d-02, 2.50151653d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.68659763d-01, 1.55485521d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.07580749d+00, 7.44054923d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.64928534d+00,-3.57425518d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.20180998d+00, 1.20087524d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.85255671d-01, 7.51222226d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 1.55827134d-01,-1.97201451d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.29133611d-01,-3.02816462d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/ 7.43349031d-02, 3.84207298d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/ 4.66767692d-02,-4.36727289d-06/
      data fpp( 2,12,1),fpp( 2,12,2)/ 4.31807170d-02, 2.68803618d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 3.18926904d-02, 1.25895826d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/ 8.17500377d-02, 3.71414335d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.02838555d+00, 1.92044283d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.55078915d+00,-4.04313367d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.61549676d+00, 5.12637833d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 1.76360119d-01, 6.99954333d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 2.13540083d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.67642078d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 1.44065179d-02, 1.24186156d-04/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.05444578d-01, 8.04813454d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 9.39256288d-01, 1.24644027d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.04848755d+00,-1.74564356d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 3.93525062d-01, 4.53776221d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.99167650d-01, 4.56246677d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-2.82562182d-02,-7.97689288d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.62476396d-02, 1.11103835d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/-1.04259108d-02, 1.81267754d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/ 9.54005796d-03,-1.64201400d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 1.79399903d-02, 4.09057847d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 3.18130837d-02, 4.48990012d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 7.87711179d-02, 3.82630211d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 1.55491819d-01, 5.00106156d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.33110673d+00,-1.55201284d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 7.97393839d-01, 8.22371865d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.21349583d-01, 3.57852090d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 1.02164846d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-1.68229682d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 4.53993620d-02, 2.42483646d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 1.03806925d-01, 4.66678224d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 4.71502357d-01, 3.42767394d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 8.72199440d-01,-8.06111182d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 4.64929769d-01, 2.49321987d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-5.58112715d-02, 1.43719232d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 6.26273890d-03,-3.71849146d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 1.82694711d-03,-2.89357322d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/-4.41625974d-03, 1.49252075d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-6.19200101d-03,-2.61712569d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/ 3.19843216d-02, 3.57718201d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.18949748d-02, 2.85719764d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/ 6.91254907d-02, 2.46648274d-04/
      data fpp( 4,15,1),fpp( 4,15,2)/ 3.54607172d-01, 1.09522927d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 7.51453943d-01,-5.99623980d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.58337885d-01,-1.69770053d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 8.72365496d-02, 1.69274002d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 4.68780999d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-2.01003625d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.47823966d-01,-5.04397508d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 4.72327204d-02, 2.88348628d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 1.84649283d-01,-3.16487606d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.27799687d-01,-3.18887586d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 1.73125864d-01, 1.12107103d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.06562436d-01, 5.24751741d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/-2.83973744d-03,-4.48317994d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-4.30514885d-03, 5.27402351d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/-7.18905022d-03, 7.33170536d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/-8.29205393d-03,-1.70088450d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-2.55922767d-02, 1.21936745d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/ 8.36701699d-03, 2.85761469d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 5.38669191d-02, 1.41347738d-04/
      data fpp( 5,15,1),fpp( 5,15,2)/ 1.62419494d-01, 1.42109017d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.27462502d-01,-2.18663345d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 1.78934623d-01,-3.88495226d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 4.87792190d-02, 6.83974350d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 1.96377782d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-7.97515636d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-3.87284986d-02,-2.63238728d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.11837807d-01, 1.00837055d-04/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 6.02805122d-02,-3.00683466d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.84306812d-01,-8.93676685d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 9.75567770d-02, 5.41850206d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.43384728d-02,-1.27064138d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-1.05378913d-03,-2.13653655d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-1.75135171d-03, 3.42787578d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-5.31753939d-03, 4.53986239d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-1.32147833d-02,-1.41353253d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/-1.78021490d-03, 6.76743893d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.17304278d-03, 2.38395696d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/ 2.76183270d-03, 7.31822826d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/ 7.85498539d-02, 3.22529997d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.52936048d-01,-7.77734825d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 9.42936228d-02,-1.56753701d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/-2.05134255d-02, 3.14096289d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 3.43235186d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 1.98691169d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.02670399d-02,-7.16823380d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-5.75314924d-02, 1.42183818d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-2.56463313d-02,-3.25827039d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 7.06680632d-02, 6.75643394d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 2.03720283d-02, 1.89249682d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-9.48854494d-03,-3.19543067d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.45106019d-04,-5.46174136d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-1.00444431d-03, 2.62127214d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/-4.30079222d-03,-9.79347186d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/-2.70381303d-03,-5.65188339d-06/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.02168637d-02, 2.43688075d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-1.47984587d-03, 1.55903604d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.62057501d-02, 3.78916777d-05/
      data fpp( 7,15,1),fpp( 7,15,2)/ 5.79360909d-02, 1.08569289d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 8.86233050d-02,-2.82913933d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/ 6.56108856d-02,-5.24135565d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/ 3.82948293d-03,-1.77691841d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-2.01319079d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 2.21332374d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 7.16658037d-04, 7.49352512d-06/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-2.60122345d-03,-1.57833379d-05/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 2.05998130d-02,-8.83017347d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 9.64593459d-03, 1.13840318d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 3.04010992d-03, 6.23604630d-06/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.59265255d-03,-2.52882170d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 2.14421321d-03,-2.18717840d-06/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-5.15871065d-04, 3.21493056d-06/
      data fpp( 8,10,1),fpp( 8,10,2)/-2.75429174d-03,-4.34854386d-06/
      data fpp( 8,11,1),fpp( 8,11,2)/-4.63496462d-03,-9.28755122d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-7.58733023d-03, 3.72556435d-06/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.67675737d-02, 7.00249773d-06/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.05698330d-02, 2.17364447d-05/
      data fpp( 8,15,1),fpp( 8,15,2)/ 1.12307826d-02, 1.49397234d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.46307317d-02,-9.33933822d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 3.75783466d-03,-1.39783705d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.46049374d-03,-1.25311799d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.28229101d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 1.88048751d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-4.92704178d-04, 7.30124989d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-1.68433346d-03,-1.57958746d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-5.53252322d-03, 1.68424860d-06/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-5.47183538d-03, 2.46488024d-06/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-4.11509390d-03,-3.79776955d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.88256483d-03,-6.62380202d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 1.13991338d-03,-2.13702235d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.24358535d-03, 2.84189142d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/-3.41297867d-03,-5.88854331d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-2.35319963d-03, 4.06228184d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-5.03207744d-03, 5.11415952d-07/
      data fpp( 9,13,1),fpp( 9,13,2)/-5.33110587d-03, 4.36205435d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-8.38087606d-03, 1.36483666d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/-1.28603933d-02, 3.52247911d-06/
      data fpp( 9,16,1),fpp( 9,16,2)/-7.34634771d-03, 1.50571694d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/-5.13144680d-03,-2.03813469d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/ 3.39377731d-03, 9.63670534d-07/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 9.95666473d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 8.35177033d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.32165867d-03, 3.49845934d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.13855730d-03,-5.97760769d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/ 2.21277993d-03, 7.25971429d-07/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.20515695d-03,-1.75627802d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 2.47026568d-03,-2.26285933d-06/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.87333815d-03,-2.35628465d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/ 2.22488328d-03, 2.30999795d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.87097032d-03,-5.95970714d-06/
      data fpp(10,10,1),fpp(10,10,2)/ 1.42120641d-03, 3.91883062d-06/
      data fpp(10,11,1),fpp(10,11,2)/-1.34598686d-03,-4.45615339d-07/
      data fpp(10,12,1),fpp(10,12,2)/-2.64810999d-03, 6.47630735d-07/
      data fpp(10,13,1),fpp(10,13,2)/-4.15925279d-03, 2.75109240d-06/
      data fpp(10,14,1),fpp(10,14,2)/-4.41666278d-03, 7.09799967d-06/
      data fpp(10,15,1),fpp(10,15,2)/-2.70795955d-03, 3.88490893d-06/
      data fpp(10,16,1),fpp(10,16,2)/-2.96034090d-03,-4.92563538d-06/
      data fpp(10,17,1),fpp(10,17,2)/ 5.82952536d-04,-7.33636740d-06/
      data fpp(10,18,1),fpp(10,18,2)/-1.27310299d-03,-5.88895029d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 3.54194751d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-1.54446731d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/-2.47207882d-04,-6.75065388d-07/
      data fpp(11, 3,1),fpp(11, 3,2)/ 2.16326048d-03, 1.50872886d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 1.77181083d-03,-9.07850048d-07/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.93070329d-03,-1.04532867d-06/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.99931867d-03,-1.03683528d-06/
      data fpp(11, 7,1),fpp(11, 7,2)/ 1.36793081d-03,-4.65330207d-07/
      data fpp(11, 8,1),fpp(11, 8,2)/ 4.04289494d-04,-1.43984389d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/ 2.09922488d-03, 1.02870577d-06/
      data fpp(11,10,1),fpp(11,10,2)/-1.35360138d-04, 6.73020800d-07/
      data fpp(11,11,1),fpp(11,11,2)/-5.78087594d-04,-6.72788973d-07/
      data fpp(11,12,1),fpp(11,12,2)/-1.11774209d-03, 1.58135091d-07/
      data fpp(11,13,1),fpp(11,13,2)/-1.44500526d-03, 1.37824861d-06/
      data fpp(11,14,1),fpp(11,14,2)/-1.59451314d-03, 3.23287047d-06/
      data fpp(11,15,1),fpp(11,15,2)/-1.33083102d-03, 2.27426952d-06/
      data fpp(11,16,1),fpp(11,16,2)/ 2.11905413d-04,-3.49794853d-06/
      data fpp(11,17,1),fpp(11,17,2)/ 2.98928311d-04,-2.91047540d-06/
      data fpp(11,18,1),fpp(11,18,2)/ 4.24548907d-04,-3.52149887d-07/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 1.64307494d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 2.89435996d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/ 5.40772855d-04, 1.99128008d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.86000795d-04, 1.50051973d-07/
      data fpp(12, 4,1),fpp(12, 4,2)/ 6.47976755d-04,-4.33358996d-08/
      data fpp(12, 5,1),fpp(12, 5,2)/ 7.32029908d-04, 3.71291626d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.04459624d-04,-1.38183060d-06/
      data fpp(12, 7,1),fpp(12, 7,2)/ 5.42138598d-04,-9.87969214d-07/
      data fpp(12, 8,1),fpp(12, 8,2)/ 1.01555875d-03, 9.23707461d-07/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.94870783d-04, 9.59139372d-07/
      data fpp(12,10,1),fpp(12,10,2)/ 5.83414360d-06,-2.00264948d-07/
      data fpp(12,11,1),fpp(12,11,2)/-3.68862764d-04,-7.76079578d-07/
      data fpp(12,12,1),fpp(12,12,2)/-5.87321659d-04,-1.15416738d-07/
      data fpp(12,13,1),fpp(12,13,2)/-7.38326164d-04, 7.87746532d-07/
      data fpp(12,14,1),fpp(12,14,2)/-7.68484671d-04, 1.23643061d-06/
      data fpp(12,15,1),fpp(12,15,2)/-5.41516363d-04, 2.09653102d-06/
      data fpp(12,16,1),fpp(12,16,2)/-6.35280749d-04,-2.58455471d-06/
      data fpp(12,17,1),fpp(12,17,2)/-2.69065778d-04,-1.77831219d-06/
      data fpp(12,18,1),fpp(12,18,2)/ 7.65073600d-05, 4.09803484d-07/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 2.20309826d-06/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-9.63699308d-07/
      data fpp(13, 2,1),fpp(13, 2,2)/-6.05146241d-05,-4.39601384d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 1.12967377d-04, 6.64104845d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/ 9.63643199d-05,-2.60817996d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/ 1.75158633d-04, 1.27167139d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/ 3.30361791d-04, 1.48149440d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/ 2.61818799d-04, 3.84235103d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.62098717d-06, 2.22910150d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.84247922d-05, 5.01242990d-08/
      data fpp(13,10,1),fpp(13,10,2)/-3.86223618d-05,-1.71407346d-07/
      data fpp(13,11,1),fpp(13,11,2)/-2.39679109d-05,-5.76494916d-07/
      data fpp(13,12,1),fpp(13,12,2)/-3.17639782d-05,-8.46129900d-08/
      data fpp(13,13,1),fpp(13,13,2)/-3.73188769d-05,-3.05312386d-09/
      data fpp(13,14,1),fpp(13,14,2)/ 9.71058062d-06, 1.44825485d-07/
      data fpp(13,15,1),fpp(13,15,2)/ 1.29964600d-04, 2.57751182d-07/
      data fpp(13,16,1),fpp(13,16,2)/ 3.95289541d-04, 4.98169786d-07/
      data fpp(13,17,1),fpp(13,17,2)/ 3.31333179d-04, 1.73569673d-07/
      data fpp(13,18,1),fpp(13,18,2)/ 8.02034664d-05,-3.70448478d-07/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-8.51775761d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00,-3.82935489d-07/
      data fpp(14, 2,1),fpp(14, 2,2)/ 1.54074448d-05,-1.78129022d-07/
      data fpp(14, 3,1),fpp(14, 3,2)/-3.30525287d-05, 6.34515753d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-2.22313372d-05, 1.76322720d-07/
      data fpp(14, 5,1),fpp(14, 5,2)/-4.90908534d-05, 1.61257544d-07/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.03715186d-04, 7.26471046d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/-9.69756967d-05, 5.81540379d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/-2.93664119d-05, 7.36743750d-10/
      data fpp(14, 9,1),fpp(14, 9,2)/-1.89610150d-05, 1.08989871d-08/
      data fpp(14,10,1),fpp(14,10,2)/-5.99998649d-06,-8.03326921d-08/
      data fpp(14,11,1),fpp(14,11,2)/-5.26488528d-06,-1.09568219d-07/
      data fpp(14,12,1),fpp(14,12,2)/-8.39723569d-06,-5.13944327d-08/
      data fpp(14,13,1),fpp(14,13,2)/-1.41802872d-05,-2.68540503d-08/
      data fpp(14,14,1),fpp(14,14,2)/-3.15894065d-05, 2.08106337d-08/
      data fpp(14,15,1),fpp(14,15,2)/-6.43856193d-05, 9.96115153d-08/
      data fpp(14,16,1),fpp(14,16,2)/-1.31428249d-04, 1.50743305d-07/
      data fpp(14,17,1),fpp(14,17,2)/-1.04666648d-04, 4.74152640d-08/
      data fpp(14,18,1),fpp(14,18,2)/-2.62640792d-05,-1.06404361d-07/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00,-2.45797819d-07/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00,-1.13452642d-07/
      data fpp(15, 2,1),fpp(15, 2,2)/-5.46515508d-06,-5.90947154d-08/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.99273785d-06, 1.83150401d-09/
      data fpp(15, 4,1),fpp(15, 4,2)/-1.58897108d-06, 3.97686994d-08/
      data fpp(15, 5,1),fpp(15, 5,2)/ 2.45478064d-06, 4.30936985d-08/
      data fpp(15, 6,1),fpp(15, 6,2)/ 1.27989512d-05, 4.58565067d-08/
      data fpp(15, 7,1),fpp(15, 7,2)/ 9.83398745d-06, 2.54802746d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-4.81336527d-06, 3.22223950d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/-4.28114773d-06,-1.63698545d-08/
      data fpp(15,10,1),fpp(15,10,2)/-3.82769225d-06,-2.07429769d-08/
      data fpp(15,11,1),fpp(15,11,2)/-2.07254800d-06,-3.86582377d-08/
      data fpp(15,12,1),fpp(15,12,2)/-2.89707900d-06,-2.26240721d-08/
      data fpp(15,13,1),fpp(15,13,2)/-4.20997445d-06, 9.15452627d-09/
      data fpp(15,14,1),fpp(15,14,2)/-2.60295478d-06, 2.20059671d-08/
      data fpp(15,15,1),fpp(15,15,2)/ 3.67787680d-06, 3.48216055d-08/
      data fpp(15,16,1),fpp(15,16,2)/ 2.00234535d-05, 3.07076109d-08/
      data fpp(15,17,1),fpp(15,17,2)/ 1.74334130d-05,-1.65204909d-09/
      data fpp(15,18,1),fpp(15,18,2)/ 4.00285058d-06,-3.60994145d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00,-6.99502927d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00,-3.63556834d-08/
      data fpp(16, 2,1),fpp(16, 2,2)/-1.18777960d-06,-2.32886331d-08/
      data fpp(16, 3,1),fpp(16, 3,2)/-1.20074404d-05,-1.44897840d-08/
      data fpp(16, 4,1),fpp(16, 4,2)/-2.18158716d-05, 9.24776910d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/-2.16220332d-05, 1.94987076d-08/
      data fpp(16, 6,1),fpp(16, 6,2)/-1.48530470d-05, 3.27574005d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/-1.12949373d-06, 2.94716902d-08/
      data fpp(16, 8,1),fpp(16, 8,2)/ 1.27888255d-05, 5.35583855d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/ 1.21445024d-05,-2.89504444d-09/
      data fpp(16,10,1),fpp(16,10,2)/ 1.06256318d-05,-2.37756608d-08/
      data fpp(16,11,1),fpp(16,11,2)/ 8.41841686d-06,-2.20023124d-08/
      data fpp(16,12,1),fpp(16,12,2)/ 1.19217538d-05,-2.02150896d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 1.55867729d-05, 2.48626708d-08/
      data fpp(16,14,1),fpp(16,14,2)/ 1.50361202d-05, 1.67644065d-08/
      data fpp(16,15,1),fpp(16,15,2)/ 9.29749017d-06, 1.60797032d-08/
      data fpp(16,16,1),fpp(16,16,2)/-2.12601244d-06,-3.08321942d-09/
      data fpp(16,17,1),fpp(16,17,2)/-4.60027792d-06, 1.42531745d-08/
      data fpp(16,18,1),fpp(16,18,2)/-8.66782432d-07,-2.99294784d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00,-5.65352608d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-9.43480145d-03, 4.19842491d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.29669579d-02, 2.37619177d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/-1.16268632d-01, 6.12304080d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 8.62854393d-03,-3.35248191d-03/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.60332464d-05, 4.67604496d-04/
      data fpppp( 1,11),fpppp( 1,12)/ 4.34856539d-04,-1.09268109d-04/
      data fpppp( 1,13),fpppp( 1,14)/-2.83999865d-04, 6.04398135d-03/
      data fpppp( 1,15),fpppp( 1,16)/ 6.33270545d-02,-1.27439592d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 4.37075172d-02, 8.63952571d-03/
      data fpppp( 1,19) /             3.24901344d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/-6.25410897d-03, 2.52149891d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.66377512d-02, 1.20617481d-02/
      data fpppp( 2, 5),fpppp( 2, 6)/-7.29049359d-02, 3.83008027d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 5.55698807d-03,-1.30120857d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-1.88052869d-04, 3.67108903d-04/
      data fpppp( 2,11),fpppp( 2,12)/ 3.48051715d-04,-3.09590852d-04/
      data fpppp( 2,13),fpppp( 2,14)/ 4.22793220d-04, 2.28714042d-03/
      data fpppp( 2,15),fpppp( 2,16)/ 4.42353351d-02,-8.46823958d-02/
      data fpppp( 2,17),fpppp( 2,18)/ 2.70324890d-02, 6.32178479d-03/
      data fpppp( 2,19) /             2.34469632d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.20157810d-03, 3.79260279d-04/
      data fpppp( 3, 3),fpppp( 3, 4)/ 1.02824295d-02,-8.94255923d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-1.19870193d-02, 1.10390107d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/-4.53271882d-03, 5.10787720d-03/
      data fpppp( 3, 9),fpppp( 3,10)/-1.53284319d-03, 6.52284555d-04/
      data fpppp( 3,11),fpppp( 3,12)/-2.27640631d-04,-4.35684218d-04/
      data fpppp( 3,13),fpppp( 3,14)/ 2.29876716d-03,-6.77428797d-03/
      data fpppp( 3,15),fpppp( 3,16)/ 2.65841448d-02,-3.36286387d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 5.37074217d-03, 3.60578792d-03/
      data fpppp( 3,19) /             1.34877865d-02 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.70048158d-03,-1.34384126d-04/
      data fpppp( 4, 3),fpppp( 4, 4)/ 4.01851018d-03, 2.61761551d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.25088732d-02,-1.06012796d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/ 9.94110289d-03,-3.73538055d-03/
      data fpppp( 4, 9),fpppp( 4,10)/ 1.00983120d-03,-4.12389153d-04/
      data fpppp( 4,11),fpppp( 4,12)/ 9.07773347d-04,-8.21580404d-04/
      data fpppp( 4,13),fpppp( 4,14)/ 8.26081072d-05, 2.73033974d-03/
      data fpppp( 4,15),fpppp( 4,16)/ 3.89110286d-03,-1.16128458d-02/
      data fpppp( 4,17),fpppp( 4,18)/ 1.16251057d-03, 2.28368690d-03/
      data fpppp( 4,19) /             6.73462895d-03 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 7.17555373d-03, 3.94572762d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.38562509d-03, 2.13836531d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/-5.82380562d-03, 3.28740349d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/-2.03918461d-03, 2.29901018d-03/
      data fpppp( 5, 9),fpppp( 5,10)/-6.80650373d-04, 3.38481916d-04/
      data fpppp( 5,11),fpppp( 5,12)/-5.66423434d-04, 9.55378679d-04/
      data fpppp( 5,13),fpppp( 5,14)/-1.79520298d-04, 4.55139023d-04/
      data fpppp( 5,15),fpppp( 5,16)/ 2.14212454d-03,-5.63421114d-03/
      data fpppp( 5,17),fpppp( 5,18)/ 1.58046673d-03, 4.14692715d-04/
      data fpppp( 5,19) /             1.64333352d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.89783285d-03,-8.27927352d-04/
      data fpppp( 6, 3),fpppp( 6, 4)/ 4.14669366d-03,-1.04518965d-03/
      data fpppp( 6, 5),fpppp( 6, 6)/-2.85145621d-03,-1.95565666d-04/
      data fpppp( 6, 7),fpppp( 6, 8)/ 2.12500601d-03,-7.93662350d-04/
      data fpppp( 6, 9),fpppp( 6,10)/ 2.10708619d-04,-2.21289634d-04/
      data fpppp( 6,11),fpppp( 6,12)/ 4.14586544d-04,-2.77147808d-04/
      data fpppp( 6,13),fpppp( 6,14)/-7.56390859d-05, 1.01936635d-03/
      data fpppp( 6,15),fpppp( 6,16)/ 1.89362419d-04,-1.86092564d-03/
      data fpppp( 6,17),fpppp( 6,18)/-7.27377032d-04, 1.40055640d-03/
      data fpppp( 6,19) /             3.24437988d-03 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.51126020d-03,-4.10342497d-04/
      data fpppp( 7, 3),fpppp( 7, 4)/ 9.32785434d-04, 1.42817758d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.77974175d-03, 8.94163658d-04/
      data fpppp( 7, 7),fpppp( 7, 8)/ 4.29214828d-04,-2.58782242d-04/
      data fpppp( 7, 9),fpppp( 7,10)/-6.25249078d-06, 1.37571628d-04/
      data fpppp( 7,11),fpppp( 7,12)/-2.50434396d-04, 3.17564164d-04/
      data fpppp( 7,13),fpppp( 7,14)/-4.48181474d-05, 3.98623112d-04/
      data fpppp( 7,15),fpppp( 7,16)/-1.06989609d-04,-6.33252276d-04/
      data fpppp( 7,17),fpppp( 7,18)/-5.81979299d-04, 6.35030474d-04/
      data fpppp( 7,19) /             1.51897259d-03 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.95948767d-04,-1.40207799d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/ 6.14707594d-04,-7.27487502d-04/
      data fpppp( 8, 5),fpppp( 8, 6)/ 2.45947530d-04, 4.58060600d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 4.52320837d-05,-6.55678587d-05/
      data fpppp( 8, 9),fpppp( 8,10)/ 2.43406547d-05,-6.49494421d-06/
      data fpppp( 8,11),fpppp( 8,12)/ 2.31039902d-05,-1.50222581d-04/
      data fpppp( 8,13),fpppp( 8,14)/ 2.04113659d-04, 2.56446997d-04/
      data fpppp( 8,15),fpppp( 8,16)/-2.93729156d-04,-1.85570363d-04/
      data fpppp( 8,17),fpppp( 8,18)/ 1.79639834d-04,-1.84556062d-05/
      data fpppp( 8,19) /            -5.56065789d-05 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.68838263d-06, 2.80818868d-06/
      data fpppp( 9, 3),fpppp( 9, 4)/-6.08566439d-05, 8.12247587d-05/
      data fpppp( 9, 5),fpppp( 9, 6)/-2.95097352d-05, 1.14577401d-04/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.50344835d-04, 8.23833272d-05/
      data fpppp( 9, 9),fpppp( 9,10)/-1.28409069d-04, 1.45638790d-04/
      data fpppp( 9,11),fpppp( 9,12)/-1.11165507d-04, 7.47038254d-05/
      data fpppp( 9,13),fpppp( 9,14)/-4.48588317d-05,-6.03130042d-05/
      data fpppp( 9,15),fpppp( 9,16)/ 2.00326027d-04,-1.41377337d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 1.67234644d-04,-1.48941846d-04/
      data fpppp( 9,19) /            -2.86607344d-04 /
      data fpppp(10, 1),fpppp(10, 2)/-1.98260102d-06,-6.41927995d-06/
      data fpppp(10, 3),fpppp(10, 4)/-2.62588170d-06,-2.76377536d-05/
      data fpppp(10, 5),fpppp(10, 6)/ 4.82661597d-05,-2.90629829d-05/
      data fpppp(10, 7),fpppp(10, 8)/ 1.62635959d-05,-9.90830408d-05/
      data fpppp(10, 9),fpppp(10,10)/ 1.73224643d-04,-1.50533712d-04/
      data fpppp(10,11),fpppp(10,12)/ 6.53480057d-05,-2.29541018d-05/
      data fpppp(10,13),fpppp(10,14)/ 1.39272209d-05, 4.24691866d-05/
      data fpppp(10,15),fpppp(10,16)/-6.58371735d-05, 1.03214432d-04/
      data fpppp(10,17),fpppp(10,18)/-1.19280065d-04, 4.99448925d-05/
      data fpppp(10,19) /             1.07250006d-04 /
      data fpppp(11, 1),fpppp(11, 2)/ 7.70700961d-05, 3.47819799d-05/
      data fpppp(11, 3),fpppp(11, 4)/-5.67374414d-05, 2.40527055d-05/
      data fpppp(11, 5),fpppp(11, 6)/-6.45285435d-06,-3.65791217d-06/
      data fpppp(11, 7),fpppp(11, 8)/-2.09156919d-05, 6.73854721d-05/
      data fpppp(11, 9),fpppp(11,10)/-8.91115939d-05, 5.32896788d-05/
      data fpppp(11,11),fpppp(11,12)/-1.65356671d-05, 7.03736716d-06/
      data fpppp(11,13),fpppp(11,14)/ 1.12967773d-06,-8.90760178d-07/
      data fpppp(11,15),fpppp(11,16)/ 2.72247623d-05,-3.12650297d-05/
      data fpppp(11,17),fpppp(11,18)/ 1.04925443d-05,-8.38928571d-06/
      data fpppp(11,19) /            -9.94557172d-06 /
      data fpppp(12, 1),fpppp(12, 2)/-2.47744007d-05,-1.15634048d-05/
      data fpppp(12, 3),fpppp(12, 4)/ 1.72953248d-05,-8.61301335d-06/
      data fpppp(12, 5),fpppp(12, 6)/-5.51863988d-06, 5.99016666d-06/
      data fpppp(12, 7),fpppp(12, 8)/ 9.47292873d-06,-2.37374111d-05/
      data fpppp(12, 9),fpppp(12,10)/ 1.38302289d-05,-5.68442506d-06/
      data fpppp(12,11),fpppp(12,12)/ 3.76785527d-06,-1.27152672d-08/
      data fpppp(12,13),fpppp(12,14)/ 3.30269209d-07, 5.94239832d-06/
      data fpppp(12,15),fpppp(12,16)/-8.67225362d-06, 9.50265453d-06/
      data fpppp(12,17),fpppp(12,18)/-1.73960308d-06,-3.78275220d-06/
      data fpppp(12,19) /            -8.45421803d-06 /
      data fpppp(13, 1),fpppp(13, 2)/ 6.26250261d-06, 2.97652720d-06/
      data fpppp(13, 3),fpppp(13, 4)/-4.12881391d-06, 2.13362494d-06/
      data fpppp(13, 5),fpppp(13, 6)/ 1.31815638d-06,-2.82171978d-06/
      data fpppp(13, 7),fpppp(13, 8)/-3.45604625d-06, 4.95209712d-06/
      data fpppp(13, 9),fpppp(13,10)/-1.55418335d-06, 1.06101040d-06/
      data fpppp(13,11),fpppp(13,12)/-5.98737040d-07,-1.30933386d-08/
      data fpppp(13,13),fpppp(13,14)/ 7.85580513d-07, 2.58326619d-08/
      data fpppp(13,15),fpppp(13,16)/ 3.50456256d-06,-5.33982765d-06/
      data fpppp(13,17),fpppp(13,18)/-1.90213015d-06, 1.71794722d-06/
      data fpppp(13,19) /             5.28591605d-06 /
      data fpppp(14, 1),fpppp(14, 2)/-1.76975585d-06,-8.39482511d-07/
      data fpppp(14, 3),fpppp(14, 4)/ 1.29564080d-06,-7.86210781d-07/
      data fpppp(14, 5),fpppp(14, 6)/-4.11640142d-07, 7.66882393d-07/
      data fpppp(14, 7),fpppp(14, 8)/ 1.02593984d-06,-1.21845399d-06/
      data fpppp(14, 9),fpppp(14,10)/ 4.15642849d-07,-2.90779507d-07/
      data fpppp(14,11),fpppp(14,12)/ 1.39195424d-08, 3.05424015d-09/
      data fpppp(14,13),fpppp(14,14)/-1.85178566d-07, 4.00959524d-08/
      data fpppp(14,15),fpppp(14,16)/-8.98430854d-07, 1.49884247d-06/
      data fpppp(14,17),fpppp(14,18)/ 5.31314778d-07,-5.25643494d-07/
      data fpppp(14,19) /            -1.55705017d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 3.49599340d-07, 1.67745757d-07/
      data fpppp(15, 3),fpppp(15, 4)/-2.45199487d-07, 1.50676079d-07/
      data fpppp(15, 5),fpppp(15, 6)/ 1.00022809d-07,-1.72742182d-07/
      data fpppp(15, 7),fpppp(15, 8)/-2.07602145d-07, 3.02207426d-07/
      data fpppp(15, 9),fpppp(15,10)/-9.04533420d-08, 5.48802182d-08/
      data fpppp(15,11),fpppp(15,12)/-5.09662040d-08,-5.79591706d-09/
      data fpppp(15,13),fpppp(15,14)/ 4.48480047d-08, 1.59880580d-09/
      data fpppp(15,15),fpppp(15,16)/ 2.29185487d-07,-3.14456050d-07/
      data fpppp(15,17),fpppp(15,18)/-1.07498317d-07, 9.40180033d-08/
      data fpppp(15,19) /             2.97089013d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-2.03578774d-07,-9.46610995d-08/
      data fpppp(16, 3),fpppp(16, 4)/ 4.31030326d-09, 1.38093657d-07/
      data fpppp(16, 5),fpppp(16, 6)/ 4.34512505d-08, 8.26102037d-08/
      data fpppp(16, 7),fpppp(16, 8)/ 4.33819662d-08,-2.44452115d-07/
      data fpppp(16, 9),fpppp(16,10)/ 6.06679556d-08,-5.06925599d-08/
      data fpppp(16,11),fpppp(16,12)/ 1.00801621d-07,-9.88080880d-09/
      data fpppp(16,13),fpppp(16,14)/-5.15774518d-08,-3.67496946d-08/
      data fpppp(16,15),fpppp(16,16)/-1.12702413d-07, 1.46466994d-07/
      data fpppp(16,17),fpppp(16,18)/ 6.37886669d-08,-2.91560040d-08/
      data fpppp(16,19) /            -1.19167434d-07 /
 

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

      subroutine a_2(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  1.72826600d-01 ,  1.60608400d-01 /
      data f( 1, 3),f( 1, 4) /  1.66498600d-01 ,  3.85095000d-01 /
      data f( 1, 5),f( 1, 6) /  5.41946400d-01 ,  2.50928500d-01 /
      data f( 1, 7),f( 1, 8) /  8.11914000d-02 ,  4.64065000d-02 /
      data f( 1, 9),f( 1,10) /  3.71461000d-02 ,  1.77484000d-02 /
      data f( 1,11),f( 1,12) / -2.03810000d-03 , -1.43197000d-02 /
      data f( 1,13),f( 1,14) / -2.12946000d-02 , -1.27342000d-02 /
      data f( 1,15),f( 1,16) /  1.36879400d-01 ,  4.81940100d-01 /
      data f( 1,17),f( 1,18) /  3.04907700d-01 ,  1.11729600d-01 /
      data f( 1,19) /           9.19871000d-02 /
      data f( 2, 1),f( 2, 2) /  1.42676600d-01 ,  1.36255100d-01 /
      data f( 2, 3),f( 2, 4) /  1.22629500d-01 ,  2.16592100d-01 /
      data f( 2, 5),f( 2, 6) /  2.77507300d-01 ,  1.43159200d-01 /
      data f( 2, 7),f( 2, 8) /  5.55943000d-02 ,  3.50893000d-02 /
      data f( 2, 9),f( 2,10) /  2.90853000d-02 ,  1.49394000d-02 /
      data f( 2,11),f( 2,12) / -1.52600000d-03 , -1.27994000d-02 /
      data f( 2,13),f( 2,14) / -1.83722000d-02 , -5.10080000d-03 /
      data f( 2,15),f( 2,16) /  6.59214000d-02 ,  2.27934900d-01 /
      data f( 2,17),f( 2,18) /  1.66111200d-01 ,  8.65477000d-02 /
      data f( 2,19) /           7.77550000d-02 /
      data f( 3, 1),f( 3, 2) /  1.17832500d-01 ,  1.13919400d-01 /
      data f( 3, 3),f( 3, 4) /  9.89015000d-02 ,  1.32165100d-01 /
      data f( 3, 5),f( 3, 6) /  1.46213700d-01 ,  8.03917000d-02 /
      data f( 3, 7),f( 3, 8) /  3.84923000d-02 ,  2.71266000d-02 /
      data f( 3, 9),f( 3,10) /  2.16259000d-02 ,  1.01843000d-02 /
      data f( 3,11),f( 3,12) / -2.51910000d-03 , -1.11547000d-02 /
      data f( 3,13),f( 3,14) / -1.61718000d-02 , -1.36193000d-02 /
      data f( 3,15),f( 3,16) /  2.24331000d-02 ,  1.08762800d-01 /
      data f( 3,17),f( 3,18) /  9.03551000d-02 ,  6.24265000d-02 /
      data f( 3,19) /           7.12055000d-02 /
      data f( 4, 1),f( 4, 2) /  1.00148200d-01 ,  9.48034000d-02 /
      data f( 4, 3),f( 4, 4) /  8.40340000d-02 ,  8.68625000d-02 /
      data f( 4, 5),f( 4, 6) /  8.26272000d-02 ,  4.79447000d-02 /
      data f( 4, 7),f( 4, 8) /  2.45472000d-02 ,  1.93396000d-02 /
      data f( 4, 9),f( 4,10) /  1.52793000d-02 ,  6.42970000d-03 /
      data f( 4,11),f( 4,12) / -3.10030000d-03 , -9.62700000d-03 /
      data f( 4,13),f( 4,14) / -1.39011000d-02 , -1.22922000d-02 /
      data f( 4,15),f( 4,16) /  5.31120000d-03 ,  5.17598000d-02 /
      data f( 4,17),f( 4,18) /  5.12326000d-02 ,  4.84388000d-02 /
      data f( 4,19) /           6.71047000d-02 /
      data f( 5, 1),f( 5, 2) /  1.02115300d-01 ,  8.23717000d-02 /
      data f( 5, 3),f( 5, 4) /  6.85726000d-02 ,  5.72203000d-02 /
      data f( 5, 5),f( 5, 6) /  4.78678000d-02 ,  2.87133000d-02 /
      data f( 5, 7),f( 5, 8) /  1.80750000d-02 ,  1.38470000d-02 /
      data f( 5, 9),f( 5,10) /  1.04883000d-02 ,  3.72320000d-03 /
      data f( 5,11),f( 5,12) / -3.31580000d-03 , -7.72410000d-03 /
      data f( 5,13),f( 5,14) / -1.15452000d-02 , -1.12900000d-02 /
      data f( 5,15),f( 5,16) / -1.65050000d-03 ,  2.41172000d-02 /
      data f( 5,17),f( 5,18) /  3.10547000d-02 ,  3.76402000d-02 /
      data f( 5,19) /           5.59690000d-02 /
      data f( 6, 1),f( 6, 2) /  7.27258000d-02 ,  6.41241000d-02 /
      data f( 6, 3),f( 6, 4) /  4.95203000d-02 ,  3.71768000d-02 /
      data f( 6, 5),f( 6, 6) /  2.85012000d-02 ,  1.82530000d-02 /
      data f( 6, 7),f( 6, 8) /  1.27697000d-02 ,  9.67910000d-03 /
      data f( 6, 9),f( 6,10) /  6.93350000d-03 ,  1.81050000d-03 /
      data f( 6,11),f( 6,12) / -3.25900000d-03 , -6.39400000d-03 /
      data f( 6,13),f( 6,14) / -9.34400000d-03 , -9.58280000d-03 /
      data f( 6,15),f( 6,16) / -4.37770000d-03 ,  1.10610000d-02 /
      data f( 6,17),f( 6,18) /  2.15951000d-02 ,  3.02375000d-02 /
      data f( 6,19) /           4.34810000d-02 /
      data f( 7, 1),f( 7, 2) /  5.01433000d-02 ,  4.46311000d-02 /
      data f( 7, 3),f( 7, 4) /  3.50005000d-02 ,  2.47739000d-02 /
      data f( 7, 5),f( 7, 6) /  1.81659000d-02 ,  1.23475000d-02 /
      data f( 7, 7),f( 7, 8) /  8.92450000d-03 ,  6.53110000d-03 /
      data f( 7, 9),f( 7,10) /  4.34460000d-03 ,  4.95300000d-04 /
      data f( 7,11),f( 7,12) / -2.96240000d-03 , -5.28660000d-03 /
      data f( 7,13),f( 7,14) / -7.48030000d-03 , -7.94780000d-03 /
      data f( 7,15),f( 7,16) / -5.14220000d-03 ,  4.90820000d-03 /
      data f( 7,17),f( 7,18) /  1.83376000d-02 ,  2.88668000d-02 /
      data f( 7,19) /           3.16027000d-02 /
      data f( 8, 1),f( 8, 2) /  3.41864000d-02 ,  3.06826000d-02 /
      data f( 8, 3),f( 8, 4) /  2.42320000d-02 ,  1.75522000d-02 /
      data f( 8, 5),f( 8, 6) /  1.23785000d-02 ,  8.69820000d-03 /
      data f( 8, 7),f( 8, 8) /  6.29810000d-03 ,  4.19010000d-03 /
      data f( 8, 9),f( 8,10) /  2.49150000d-03 , -3.72400000d-04 /
      data f( 8,11),f( 8,12) / -2.67430000d-03 , -4.38560000d-03 /
      data f( 8,13),f( 8,14) / -5.83440000d-03 , -6.29520000d-03 /
      data f( 8,15),f( 8,16) / -4.56200000d-03 ,  2.50040000d-03 /
      data f( 8,17),f( 8,18) /  1.29289000d-02 ,  2.08684000d-02 /
      data f( 8,19) /           2.27399000d-02 /
      data f( 9, 1),f( 9, 2) /  1.47991000d-02 ,  1.34194000d-02 /
      data f( 9, 3),f( 9, 4) /  1.08165000d-02 ,  8.35270000d-03 /
      data f( 9, 5),f( 9, 6) /  5.96930000d-03 ,  4.19050000d-03 /
      data f( 9, 7),f( 9, 8) /  2.63930000d-03 ,  1.26300000d-03 /
      data f( 9, 9),f( 9,10) /  3.01100000d-04 , -1.17600000d-03 /
      data f( 9,11),f( 9,12) / -2.12330000d-03 , -3.06170000d-03 /
      data f( 9,13),f( 9,14) / -3.83080000d-03 , -3.96950000d-03 /
      data f( 9,15),f( 9,16) / -2.75840000d-03 ,  7.91400000d-04 /
      data f( 9,17),f( 9,18) /  6.00390000d-03 ,  9.62440000d-03 /
      data f( 9,19) /           1.05664000d-02 /
      data f(10, 1),f(10, 2) /  5.51860000d-03 ,  4.99700000d-03 /
      data f(10, 3),f(10, 4) /  4.11320000d-03 ,  3.31820000d-03 /
      data f(10, 5),f(10, 6) /  2.36850000d-03 ,  1.48110000d-03 /
      data f(10, 7),f(10, 8) /  6.27700000d-04 , -5.16000000d-05 /
      data f(10, 9),f(10,10) / -6.72700000d-04 , -1.14360000d-03 /
      data f(10,11),f(10,12) / -1.63740000d-03 , -2.16960000d-03 /
      data f(10,13),f(10,14) / -2.61150000d-03 , -2.68460000d-03 /
      data f(10,15),f(10,16) / -2.01100000d-03 , -2.18100000d-04 /
      data f(10,17),f(10,18) /  2.12830000d-03 ,  3.54860000d-03 /
      data f(10,19) /           3.86810000d-03 /
      data f(11, 1),f(11, 2) /  7.85100000d-04 ,  6.38900000d-04 /
      data f(11, 3),f(11, 4) /  5.83400000d-04 ,  4.24000000d-04 /
      data f(11, 5),f(11, 6) /  1.63000000d-04 , -1.34700000d-04 /
      data f(11, 7),f(11, 8) / -4.26300000d-04 , -6.47000000d-04 /
      data f(11, 9),f(11,10) / -8.44400000d-04 , -9.51700000d-04 /
      data f(11,11),f(11,12) / -1.16230000d-03 , -1.43400000d-03 /
      data f(11,13),f(11,14) / -1.67960000d-03 , -1.76380000d-03 /
      data f(11,15),f(11,16) / -1.53320000d-03 , -8.95200000d-04 /
      data f(11,17),f(11,18) / -5.56000000d-05 ,  2.90500000d-04 /
      data f(11,19) /           3.51700000d-04 /
      data f(12, 1),f(12, 2) / -5.14300000d-04 , -4.94300000d-04 /
      data f(12, 3),f(12, 4) / -4.75100000d-04 , -5.16900000d-04 /
      data f(12, 5),f(12, 6) / -5.86900000d-04 , -6.29000000d-04 /
      data f(12, 7),f(12, 8) / -6.69500000d-04 , -7.42800000d-04 /
      data f(12, 9),f(12,10) / -7.51600000d-04 , -7.94100000d-04 /
      data f(12,11),f(12,12) / -8.62300000d-04 , -9.68500000d-04 /
      data f(12,13),f(12,14) / -1.11130000d-03 , -1.19320000d-03 /
      data f(12,15),f(12,16) / -1.15210000d-03 , -9.78900000d-04 /
      data f(12,17),f(12,18) / -7.35500000d-04 , -5.51800000d-04 /
      data f(12,19) /          -4.92300000d-04 /
      data f(13, 1),f(13, 2) / -7.23800000d-04 , -6.79600000d-04 /
      data f(13, 3),f(13, 4) / -6.76000000d-04 , -6.32800000d-04 /
      data f(13, 5),f(13, 6) / -6.03400000d-04 , -5.50200000d-04 /
      data f(13, 7),f(13, 8) / -5.02800000d-04 , -4.73800000d-04 /
      data f(13, 9),f(13,10) / -4.65000000d-04 , -4.75900000d-04 /
      data f(13,11),f(13,12) / -5.01300000d-04 , -5.19800000d-04 /
      data f(13,13),f(13,14) / -5.46400000d-04 , -5.77400000d-04 /
      data f(13,15),f(13,16) / -6.01500000d-04 , -6.01000000d-04 /
      data f(13,17),f(13,18) / -5.79000000d-04 , -5.55600000d-04 /
      data f(13,19) /          -5.47700000d-04 /
      data f(14, 1),f(14, 2) / -2.36300000d-04 , -2.32700000d-04 /
      data f(14, 3),f(14, 4) / -2.22500000d-04 , -2.07300000d-04 /
      data f(14, 5),f(14, 6) / -1.90800000d-04 , -1.77700000d-04 /
      data f(14, 7),f(14, 8) / -1.69000000d-04 , -1.65500000d-04 /
      data f(14, 9),f(14,10) / -1.61500000d-04 , -1.58300000d-04 /
      data f(14,11),f(14,12) / -1.56800000d-04 , -1.58200000d-04 /
      data f(14,13),f(14,14) / -1.61600000d-04 , -1.66700000d-04 /
      data f(14,15),f(14,16) / -1.72400000d-04 , -1.77400000d-04 /
      data f(14,17),f(14,18) / -1.79200000d-04 , -1.80800000d-04 /
      data f(14,19) /          -1.80800000d-04 /
      data f(15, 1),f(15, 2) / -8.37000000d-05 , -8.14000000d-05 /
      data f(15, 3),f(15, 4) / -7.59000000d-05 , -6.94000000d-05 /
      data f(15, 5),f(15, 6) / -6.38000000d-05 , -5.96000000d-05 /
      data f(15, 7),f(15, 8) / -5.67000000d-05 , -5.58000000d-05 /
      data f(15, 9),f(15,10) / -5.44000000d-05 , -5.42000000d-05 /
      data f(15,11),f(15,12) / -5.44000000d-05 , -5.56000000d-05 /
      data f(15,13),f(15,14) / -5.72000000d-05 , -5.92000000d-05 /
      data f(15,15),f(15,16) / -6.10000000d-05 , -6.28000000d-05 /
      data f(15,17),f(15,18) / -6.42000000d-05 , -6.62000000d-05 /
      data f(15,19) /          -6.70000000d-05 /
      data f(16, 1),f(16, 2) / -1.15000000d-05 , -1.17000000d-05 /
      data f(16, 3),f(16, 4) / -1.15000000d-05 , -1.17000000d-05 /
      data f(16, 5),f(16, 6) / -1.22000000d-05 , -1.22000000d-05 /
      data f(16, 7),f(16, 8) / -1.19000000d-05 , -1.13000000d-05 /
      data f(16, 9),f(16,10) / -1.08000000d-05 , -1.08000000d-05 /
      data f(16,11),f(16,12) / -1.12000000d-05 , -1.20000000d-05 /
      data f(16,13),f(16,14) / -1.20000000d-05 , -1.24000000d-05 /
      data f(16,15),f(16,16) / -1.24000000d-05 , -1.20000000d-05 /
      data f(16,17),f(16,18) / -1.18000000d-05 , -1.15000000d-05 /
      data f(16,19) /          -1.24000000d-05 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 5.02580583d-02,-1.61542923d-03/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.24292328d-02,-1.17845542d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 7.84438871d-01, 3.17331540d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.20794133d+00, 1.86955948d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 4.94364365d+00,-7.62583919d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.49812080d+00, 3.44424282d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 3.32379642d-01, 1.12571592d-03/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 1.56937459d-01, 1.50025504d-04/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.10638024d-03,-1.94347935d-04/
      data fpp( 1,10,1),fpp( 1,10,2)/-1.19400574d-01, 1.91282366d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-8.35461819d-02, 9.45069887d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 7.88155313d-03, 5.31378087d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/-3.71184183d-02, 1.13437764d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-1.01298982d+00, 8.33605086d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.32594272d-01, 5.11742788d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.15069407d+00,-9.57649061d-03/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.22935496d+00, 1.86294855d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.82320471d-01, 1.15595441d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 3.31832120d-01, 3.91936979d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 2.04731383d-01,-1.13011246d-03/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 6.63565345d-02,-2.51744084d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 5.05734758d-01, 1.70484279d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.13738483d+00,-1.12335094d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.37054520d+00,-3.23834642d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.11295840d+00, 1.34992277d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.39283217d-01, 6.45647322d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 9.66525822d-02, 9.10819383d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.53222395d-02,-1.39915075d-04/
      data fpp( 2,10,1),fpp( 2,10,2)/-5.44863511d-02,-1.99356385d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-4.16526363d-02, 8.04876290d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/ 5.63689375d-03, 9.50512241d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/-1.95281633d-02, 2.23527881d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/-4.85287859d-01, 2.27035352d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 6.50203956d-01, 2.33337871d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.44429436d+00,-4.10107219d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.58966007d+00, 6.40678051d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-9.44655865d-03, 4.73971985d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 1.74225761d-01, 1.70968201d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-7.32985916d-02,-6.57153254d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 2.47846293d-02,-2.06563492d-04/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.13787098d-01, 8.17119220d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 8.53904343d-01,-1.65023390d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.54600056d+00,-1.30992566d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 8.00315587d-01, 6.12490026d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-1.52475092d-02, 2.95321554d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-4.03727878d-02, 3.82457580d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 2.68146617d-02,-9.64045858d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/ 4.54309789d-02,-9.08141491d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/ 2.43767270d-02, 5.70222454d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-1.17691281d-02, 2.50604333d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 6.93107171d-03, 5.98460216d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 5.31356257d-01, 1.89731481d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 7.87044903d-01, 1.19122206d-03/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.29709348d+00,-1.93798171d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 8.68064743d-01, 2.76460767d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/ 3.79211705d-01, 2.60884638d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 1.23654838d-01, 8.82456681d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 1.16243298d+00,-2.30483058d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 3.17459948d-01,-8.22218833d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-3.18081509d-02, 2.33894592d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 3.15657799d-01,-3.74824835d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 6.01517564d-01,-5.07792658d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 2.33854249d-01, 2.41821114d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 2.95241820d-01, 2.17608202d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 9.11935689d-02,-2.08599204d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 4.43391137d-02,-6.53305199d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/ 2.28374355d-02,-5.17599987d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/ 5.93072832d-03, 4.52105194d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/ 2.38896187d-02, 4.53192220d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 2.34887652d-03, 7.18177918d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/-1.63297168d-01, 6.11769107d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 1.56576430d-01, 6.43144565d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 6.92696729d-01,-9.03043172d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.33120954d-01, 1.50480123d-04/
      data fpp( 4,18,1),fpp( 4,18,2)/ 1.26247375d-02, 1.65126679d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/-3.01540113d-01, 4.76595160d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/-1.62872334d+00, 9.08439652d-05/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-2.91979423d-01, 6.66010696d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-1.75639495d-01,-5.78243468d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.32524463d-01, 8.25199043d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.71994183d-01,-2.09513374d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.46607416d-01, 1.67413591d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-4.47847708d-02, 5.08310107d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.97585123d-02, 1.38803664d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 2.91688835d-02,-5.41944762d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/ 2.04342790d-02,-1.48646158d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 6.75535975d-03, 4.37063225d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-2.75093468d-02,-1.54968285d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/-3.54657777d-03, 5.35129913d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 7.30974159d-02, 4.60228632d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 1.10679376d-01, 3.25453556d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.36179608d-01,-3.80145087d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 2.41141441d-01, 6.53147927d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 4.86543446d-02, 9.77659164d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 2.72706124d-02, 2.48219542d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 6.48970382d-01,-1.40957829d-04/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-2.19272569d-02,-6.33953428d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.95731130d-01, 3.44131998d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.94049350d-01, 6.13605435d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 2.19425704d-01,-5.97813737d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 9.53810850d-02, 8.34089513d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 5.89322631d-02, 1.20395684d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 2.84773820d-02, 1.19947750d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 2.44153523d-02,-3.93186686d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/ 1.44954485d-02, 2.63589936d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 7.89283269d-03, 3.19850712d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/ 2.27768489d-04,-1.45061840d-05/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.13675654d-02, 3.71396647d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/-2.33424954d-02, 2.86195251d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/ 3.58810674d-02, 1.75016235d-04/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.50544838d-01,-1.14668465d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/ 2.10058281d-01,-1.06183748d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 3.02142884d-01, 4.36399642d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/-1.03873369d-02, 1.12124518d-04/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 5.38918120d-02,-7.76353321d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.92878451d-01,-3.87293359d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 7.25899737d-02,-1.45513245d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.37368136d-01, 6.11746337d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 1.04998001d-01,-1.30312103d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 5.50882437d-02, 3.83262075d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 2.80707185d-02, 3.45038014d-06/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.93169597d-02, 9.64827190d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.80547072d-02,-2.96294677d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/ 1.12089270d-02, 9.10159909d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-2.35669052d-03, 1.67190714d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/-6.80672715d-03,-7.96788462d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-1.60816050d-03, 2.29824671d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/ 9.44256575d-03, 1.96100162d-05/
      data fpp( 7,15,1),fpp( 7,15,2)/ 4.02013549d-02, 9.49634680d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 9.71510384d-02, 3.52241117d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/-1.51059566d-01,-3.31199147d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/-3.52425881d-01,-7.67564529d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 1.05733735d-01,-1.27452274d-04/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 1.29302370d-01,-5.64350806d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 8.20884543d-02,-2.98858388d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 7.66039750d-02,-8.29564219d-07/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 3.36581063d-02, 1.94520957d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 4.27672908d-02, 1.33871815d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 2.26959403d-02, 1.66031782d-05/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.16048630d-02,-2.98789424d-06/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 1.53047791d-02, 1.28743988d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 1.37358188d-02,-2.39457009d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/ 7.79384336d-03, 1.29904048d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/ 2.58929398d-04, 5.70408188d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-3.96085988d-03,-3.70732254d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.48697926d-02, 1.15288471d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.17877676d-02, 1.35353437d-05/
      data fpp( 8,15,1),fpp( 8,15,2)/ 5.01851313d-03, 6.59697781d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.26010082d-02, 4.23375440d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/ 7.14999819d-02,-3.33539542d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.13405641d-01,-5.82617274d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 3.97773958d-02,-9.76791363d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 5.48907345d-02,-2.51734869d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 5.60629117d-02,-1.35950262d-05/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 3.84493381d-02, 6.16159155d-06/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 2.69878631d-02,-2.70534004d-06/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.29091270d-02, 9.48376862d-06/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 9.02680725d-03, 1.04626557d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.09250519d-02,-1.28308817d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 1.02359329d-02, 9.49905796d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 6.60768988d-03,-1.31194010d-05/
      data fpp( 9,10,1),fpp( 9,10,2)/ 5.95650640d-03, 1.20665459d-05/
      data fpp( 9,11,1),fpp( 9,11,2)/-5.43442933d-04,-3.35878260d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-2.64280678d-03, 1.90258452d-06/
      data fpp( 9,13,1),fpp( 9,13,2)/-2.89404206d-03, 5.90644453d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-6.08923010d-03, 1.22956374d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/-1.10362168d-02, 2.58990060d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/ 1.18956343d-04, 2.44303386d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/ 6.99483718d-03,-2.38583604d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.42260190d-02,-2.45168970d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 3.60046949d-02,-3.87840515d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 3.01396921d-02,-7.78361364d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/ 2.51898990d-02,-4.31877273d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.13061726d-02, 3.32670455d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.45779413d-02,-3.66004546d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.09112011d-02, 2.03147729d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 8.63308071d-03,-7.27863719d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/ 6.46492955d-03, 2.91997758d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/ 4.22023930d-03,-5.06046600d-07/
      data fpp(10, 9,1),fpp(10, 9,2)/ 5.45592163d-03, 2.59620882d-06/
      data fpp(10,10,1),fpp(10,10,2)/-2.69868946d-04,-8.66788687d-07/
      data fpp(10,11,1),fpp(10,11,2)/-5.26407665d-04,-5.03054072d-07/
      data fpp(10,12,1),fpp(10,12,2)/-1.66041298d-03, 5.75004976d-07/
      data fpp(10,13,1),fpp(10,13,2)/-2.96528920d-03, 3.62103417d-06/
      data fpp(10,14,1),fpp(10,14,2)/-2.88531202d-03, 7.06885835d-06/
      data fpp(10,15,1),fpp(10,15,2)/-4.81145830d-04, 1.29055324d-05/
      data fpp(10,16,1),fpp(10,16,2)/ 3.15441647d-03, 8.46701193d-06/
      data fpp(10,17,1),fpp(10,17,2)/ 1.48731694d-02,-1.35635801d-05/
      data fpp(10,18,1),fpp(10,18,2)/ 2.34977834d-02,-9.77869139d-06/
      data fpp(10,19,1),fpp(10,19,2)/ 2.15238246d-02,-1.33696543d-05/
      data fpp(11, 1,1),fpp(11, 1,2)/ 1.23955210d-02, 2.77952228d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/ 1.25436343d-02, 1.05395545d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 8.92210821d-03,-1.55334406d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 7.50332085d-03,-1.07457920d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 5.48437442d-03,-2.44339136d-07/
      data fpp(11, 6,1),fpp(11, 6,2)/ 4.20226366d-03,-1.50064255d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 3.03821213d-03, 1.21059615d-06/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.76679220d-03,-4.38320364d-07/
      data fpp(11, 9,1),fpp(11, 9,2)/ 1.65730209d-04, 1.94068530d-06/
      data fpp(11,10,1),fpp(11,10,2)/-1.60076913d-04,-1.91842084d-06/
      data fpp(11,11,1),fpp(11,11,2)/-8.44778059d-04,-4.65001937d-07/
      data fpp(11,12,1),fpp(11,12,2)/-1.01686784d-03, 1.12428589d-07/
      data fpp(11,13,1),fpp(11,13,2)/-1.22312523d-03, 1.58128758d-06/
      data fpp(11,14,1),fpp(11,14,2)/-1.18929266d-03, 3.24642109d-06/
      data fpp(11,15,1),fpp(11,15,2)/-3.93701551d-04, 4.32102805d-06/
      data fpp(11,16,1),fpp(11,16,2)/ 2.58353563d-03, 3.91346671d-06/
      data fpp(11,17,1),fpp(11,17,2)/ 4.71512054d-03,-7.87889488d-06/
      data fpp(11,18,1),fpp(11,18,2)/ 8.10676467d-03,-2.00788718d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 1.02658754d-02,-1.18355641d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 2.69662402d-03, 5.60015084d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/ 2.03316387d-03, 5.99698312d-08/
      data fpp(12, 3,1),fpp(12, 3,2)/ 2.31659458d-03,-8.47894409d-07/
      data fpp(12, 4,1),fpp(12, 4,2)/ 2.28797529d-03,-3.28392194d-07/
      data fpp(12, 5,1),fpp(12, 5,2)/ 2.08570122d-03, 4.69463185d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 1.47386465d-03, 1.24539455d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 8.41421916d-04,-8.71621006d-07/
      data fpp(12, 8,1),fpp(12, 8,2)/ 7.02991888d-04, 1.39394457d-06/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.29157527d-04,-8.34157276d-07/
      data fpp(12,10,1),fpp(12,10,2)/ 8.69765964d-05,-7.93154684d-08/
      data fpp(12,11,1),fpp(12,11,2)/-2.96880100d-04,-3.90580851d-07/
      data fpp(12,12,1),fpp(12,12,2)/-7.54515653d-04,-6.38361128d-07/
      data fpp(12,13,1),fpp(12,13,2)/-8.68609868d-04, 7.48025364d-07/
      data fpp(12,14,1),fpp(12,14,2)/-7.62317335d-04, 1.30025967d-06/
      data fpp(12,15,1),fpp(12,15,2)/-2.64847964d-04, 1.43093594d-06/
      data fpp(12,16,1),fpp(12,16,2)/ 7.53041013d-04, 9.01996551d-07/
      data fpp(12,17,1),fpp(12,17,2)/ 2.36234849d-03,-8.26922148d-07/
      data fpp(12,18,1),fpp(12,18,2)/ 2.05435796d-03,-1.17630796d-06/
      data fpp(12,19,1),fpp(12,19,2)/ 1.55027381d-03,-1.91984602d-06/
      data fpp(13, 1,1),fpp(13, 1,2)/ 4.81674485d-05,-1.14039876d-06/
      data fpp(13, 2,1),fpp(13, 2,2)/ 1.15291260d-04,-5.41202473d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 8.57621542d-05, 8.69208655d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-2.01862893d-05,-5.59632148d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/-9.94908678d-05, 5.41319937d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/-1.18325781d-04,-1.77647599d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/-1.24771816d-04,-1.78729541d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/-2.28771766d-04,-2.11434237d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.64337687d-04,-1.87533512d-07/
      data fpp(13,10,1),fpp(13,10,2)/-1.62891333d-04,-2.20431714d-07/
      data fpp(13,11,1),fpp(13,11,2)/-1.20970672d-04, 1.99260367d-07/
      data fpp(13,12,1),fpp(13,12,2)/-1.21819121d-04,-1.62609756d-07/
      data fpp(13,13,1),fpp(13,13,2)/-2.12807779d-04,-3.48213447d-08/
      data fpp(13,14,1),fpp(13,14,2)/-2.70801663d-04, 3.78951346d-08/
      data fpp(13,15,1),fpp(13,15,2)/-2.78205331d-04, 2.97240806d-07/
      data fpp(13,16,1),fpp(13,16,2)/-2.79090855d-04, 2.49141639d-07/
      data fpp(13,17,1),fpp(13,17,2)/-3.46805732d-04,-3.80736448d-09/
      data fpp(13,18,1),fpp(13,18,2)/-1.31656206d-04,-1.49912182d-07/
      data fpp(13,19,1),fpp(13,19,2)/ 1.18408834d-05,-3.26543909d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.33064356d-04, 8.26329023d-08/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.36205714d-04, 6.47341954d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-1.32633753d-04, 5.44303163d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-9.74787758d-05, 1.75445396d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/-7.59780061d-05,-4.66084748d-08/
      data fpp(14, 6,1),fpp(14, 6,2)/-5.96049823d-05,-3.51106405d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/-4.57955110d-05,-7.69489633d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/-9.73064708d-06, 3.09064938d-08/
      data fpp(14, 9,1),fpp(14, 9,2)/-2.61157038d-05,-1.66770119d-08/
      data fpp(14,10,1),fpp(14,10,2)/-3.30143000d-05,-1.21984462d-08/
      data fpp(14,11,1),fpp(14,11,2)/-5.48979355d-05,-3.65292034d-08/
      data fpp(14,12,1),fpp(14,12,2)/-6.09848092d-05,-1.56847402d-08/
      data fpp(14,13,1),fpp(14,13,2)/-4.47717283d-05,-2.07318358d-08/
      data fpp(14,14,1),fpp(14,14,2)/-3.77863440d-05,-3.38791669d-09/
      data fpp(14,15,1),fpp(14,15,2)/-4.11100236d-05,-1.71649746d-09/
      data fpp(14,16,1),fpp(14,16,2)/-3.75479426d-05, 5.22539065d-08/
      data fpp(14,17,1),fpp(14,17,2)/-1.05570490d-05,-1.52991287d-08/
      data fpp(14,18,1),fpp(14,18,2)/-5.86103605d-05, 2.09426082d-08/
      data fpp(14,19,1),fpp(14,19,2)/-9.41095538d-05, 2.75286959d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/-1.82600255d-05, 5.39972070d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/-1.38684034d-05, 3.20055859d-08/
      data fpp(15, 3,1),fpp(15, 3,2)/-1.55771421d-05, 9.98044924d-09/
      data fpp(15, 4,1),fpp(15, 4,2)/-2.12986074d-05,-1.19273829d-08/
      data fpp(15, 5,1),fpp(15, 5,2)/-2.49971078d-05,-1.62709176d-08/
      data fpp(15, 6,1),fpp(15, 6,2)/-2.48542899d-05,-6.98894659d-09/
      data fpp(15, 7,1),fpp(15, 7,2)/-2.42961403d-05,-3.37732960d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-3.02056461d-05, 2.20821307d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/-2.57994983d-05,-2.45552266d-08/
      data fpp(15,10,1),fpp(15,10,2)/-2.53014672d-05, 4.13877585d-09/
      data fpp(15,11,1),fpp(15,11,2)/-2.25875866d-05,-1.59998768d-08/
      data fpp(15,12,1),fpp(15,12,2)/-2.27416415d-05,-1.39268820d-10/
      data fpp(15,13,1),fpp(15,13,2)/-2.87053076d-05,-7.44304796d-09/
      data fpp(15,14,1),fpp(15,14,2)/-3.28529611d-05, 5.91146067d-09/
      data fpp(15,15,1),fpp(15,15,2)/-3.39045743d-05,-4.20279471d-09/
      data fpp(15,16,1),fpp(15,16,2)/-3.42173748d-05, 1.08997182d-08/
      data fpp(15,17,1),fpp(15,17,2)/-3.81660725d-05,-1.53960780d-08/
      data fpp(15,18,1),fpp(15,18,2)/-2.42023522d-05, 1.46845937d-08/
      data fpp(15,19,1),fpp(15,19,2)/-1.50526680d-05, 2.86577031d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 4.51096556d-05, 1.17367899d-08/
      data fpp(16, 2,1),fpp(16, 2,2)/ 3.21984874d-05, 4.52642014d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/ 3.66796425d-05,-5.84247051d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/ 3.65478751d-05,-5.15653812d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/ 4.13756968d-05, 8.46862298d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/ 3.60842878d-05, 1.28204621d-09/
      data fpp(16, 7,1),fpp(16, 7,2)/ 2.90841416d-05, 4.40319217d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/ 2.71392516d-05,-8.94814879d-10/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.53754634d-05,-6.82393265d-09/
      data fpp(16,10,1),fpp(16,10,2)/ 3.11810908d-05,-1.80945452d-09/
      data fpp(16,11,1),fpp(16,11,2)/ 3.82905790d-05,-9.93824929d-09/
      data fpp(16,12,1),fpp(16,12,2)/ 4.29286779d-05, 1.75624517d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 5.12655110d-05,-1.23115574d-08/
      data fpp(16,14,1),fpp(16,14,2)/ 5.84007662d-05, 7.68377777d-09/
      data fpp(16,15,1),fpp(16,15,2)/ 6.14219300d-05, 5.57644629d-09/
      data fpp(16,16,1),fpp(16,16,2)/ 5.83565446d-05,-5.98956292d-09/
      data fpp(16,17,1),fpp(16,17,2)/ 5.41916077d-05, 6.38180540d-09/
      data fpp(16,18,1),fpp(16,18,2)/ 4.16147475d-05,-1.35376587d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 3.55152626d-05,-2.42311707d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 3.35435075d-04, 6.79119655d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.10900866d-02, 7.93802667d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/-9.41102021d-02, 5.76292714d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 3.80018176d-04, 2.68594305d-04/
      data fpppp( 1, 9),fpppp( 1,10)/-2.17729155d-04, 2.60176975d-03/
      data fpppp( 1,11),fpppp( 1,12)/-7.47669020d-04, 3.72330687d-03/
      data fpppp( 1,13),fpppp( 1,14)/-2.23312208d-02, 2.97492907d-02/
      data fpppp( 1,15),fpppp( 1,16)/ 6.66213878d-02,-1.35883900d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 3.65478771d-02, 2.02722116d-02/
      data fpppp( 1,19) /             5.79129581d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 9.88843517d-05, 4.84444790d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.51885084d-02, 5.93782958d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/-6.28492091d-02, 3.60141772d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.82719661d-03, 5.39709451d-04/
      data fpppp( 2, 9),fpppp( 2,10)/-3.08016898d-04, 1.38366327d-03/
      data fpppp( 2,11),fpppp( 2,12)/-2.68097839d-04, 1.75607700d-03/
      data fpppp( 2,13),fpppp( 2,14)/-1.11034854d-02, 1.62221863d-02/
      data fpppp( 2,15),fpppp( 2,16)/ 4.22898309d-02,-8.58655945d-02/
      data fpppp( 2,17),fpppp( 2,18)/ 2.22490653d-02, 1.22009927d-02/
      data fpppp( 2,19) /             3.59137008d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/-2.58460755d-03, 6.92881997d-04/
      data fpppp( 3, 3),fpppp( 3, 4)/ 5.26823445d-03, 5.30106673d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-2.33537630d-02, 1.84711394d-03/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.17726199d-02,-1.51132432d-03/
      data fpppp( 3, 9),fpppp( 3,10)/-1.88558903d-04,-6.48708006d-04/
      data fpppp( 3,11),fpppp( 3,12)/ 4.03156779d-04,-1.86941530d-03/
      data fpppp( 3,13),fpppp( 3,14)/ 1.03652677d-02,-9.24815652d-03/
      data fpppp( 3,15),fpppp( 3,16)/ 1.05031660d-02,-1.75029120d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 3.16384347d-03, 1.25807985d-03/
      data fpppp( 3,19) /             5.80160739d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 3.29619467d-03, 4.25817647d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 9.41339560d-03,-1.07715949d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.26789028d-02, 1.16119424d-02/
      data fpppp( 4, 7),fpppp( 4, 8)/-8.02581371d-03, 4.56516312d-03/
      data fpppp( 4, 9),fpppp( 4,10)/-8.03210992d-04, 1.68847471d-04/
      data fpppp( 4,11),fpppp( 4,12)/ 4.03519364d-04, 3.09010930d-04/
      data fpppp( 4,13),fpppp( 4,14)/-4.00954104d-03, 7.08283509d-03/
      data fpppp( 4,15),fpppp( 4,16)/ 4.80937926d-03,-1.33455501d-02/
      data fpppp( 4,17),fpppp( 4,18)/ 8.31056872d-04, 3.66096152d-04/
      data fpppp( 4,19) /             4.08444050d-03 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.64016470d-02,-1.40533862d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/ 9.39095218d-03,-6.00098079d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/-1.50868328d-03,-3.85567530d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/ 6.97105922d-03,-2.67243337d-03/
      data fpppp( 5, 9),fpppp( 5,10)/ 4.10699555d-04,-5.90633920d-05/
      data fpppp( 5,11),fpppp( 5,12)/-4.71104872d-04, 7.08335643d-04/
      data fpppp( 5,13),fpppp( 5,14)/ 1.13141084d-03,-2.07310551d-03/
      data fpppp( 5,15),fpppp( 5,16)/ 4.81728917d-03,-5.92095479d-03/
      data fpppp( 5,17),fpppp( 5,18)/-3.65773995d-04, 1.53711499d-03/
      data fpppp( 5,19) /             4.48351591d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.93410706d-02, 1.01324634d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/-6.55756255d-03, 2.93737680d-03/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.56845662d-03, 2.37119133d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/-6.60560884d-04, 6.30688654d-04/
      data fpppp( 6, 9),fpppp( 6,10)/-2.78622649d-04, 1.32329493d-04/
      data fpppp( 6,11),fpppp( 6,12)/-5.16580423d-05, 1.05557713d-05/
      data fpppp( 6,13),fpppp( 6,14)/-2.26381226d-04, 8.72193369d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 1.00951732d-03,-1.58385014d-03/
      data fpppp( 6,17),fpppp( 6,18)/ 2.01686357d-03,-4.52933452d-03/
      data fpppp( 6,19) /            -8.17641495d-03 /
      data fpppp( 7, 1),fpppp( 7, 2)/-6.69399206d-03,-3.27710444d-03/
      data fpppp( 7, 3),fpppp( 7, 4)/ 4.24590288d-03,-2.60250873d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/ 3.35234257d-04, 2.09194317d-04/
      data fpppp( 7, 7),fpppp( 7, 8)/ 2.01522422d-04, 8.05419803d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-7.41999690d-05,-1.18753763d-04/
      data fpppp( 7,11),fpppp( 7,12)/ 1.46024780d-04, 8.15895010d-05/
      data fpppp( 7,13),fpppp( 7,14)/ 1.06533413d-04,-1.56593578d-04/
      data fpppp( 7,15),fpppp( 7,16)/ 1.70232467d-03,-5.08125145d-03/
      data fpppp( 7,17),fpppp( 7,18)/ 3.13063861d-04, 6.63965333d-03/
      data fpppp( 7,19) /             1.26998787d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 1.11656085d-03, 6.02577905d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/-1.02310629d-03, 1.24216391d-03/
      data fpppp( 8, 5),fpppp( 8, 6)/-8.22246140d-04, 2.95988562d-04/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.77108279d-04,-1.16962072d-04/
      data fpppp( 8, 9),fpppp( 8,10)/-2.53925731d-05,-4.38485499d-05/
      data fpppp( 8,11),fpppp( 8,12)/ 1.05210463d-04,-1.78085822d-04/
      data fpppp( 8,13),fpppp( 8,14)/ 2.05784222d-04, 1.94406395d-04/
      data fpppp( 8,15),fpppp( 8,16)/-1.59954455d-04, 4.91984283d-04/
      data fpppp( 8,17),fpppp( 8,18)/ 7.10060463d-05,-1.19560738d-03/
      data fpppp( 8,19) /            -2.22061074d-03 /
      data fpppp( 9, 1),fpppp( 9, 2)/-4.20395603d-04,-2.21538301d-04/
      data fpppp( 9, 3),fpppp( 9, 4)/ 1.79403763d-04,-1.26950836d-04/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.71363919d-04, 5.32801350d-05/
      data fpppp( 9, 7),fpppp( 9, 8)/-3.76505949d-05,-5.79195712d-05/
      data fpppp( 9, 9),fpppp( 9,10)/ 9.29814360d-05,-1.35382600d-04/
      data fpppp( 9,11),fpppp( 9,12)/ 9.76230153d-05, 8.92566803d-06/
      data fpppp( 9,13),fpppp( 9,14)/-2.24379728d-05,-9.58109429d-05/
      data fpppp( 9,15),fpppp( 9,16)/ 3.00573823d-04,-1.40354755d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 4.08765797d-06, 1.45322184d-04/
      data fpppp( 9,19) /             2.87473248d-04 /
      data fpppp(10, 1),fpppp(10, 2)/ 4.48320515d-05, 2.05293273d-05/
      data fpppp(10, 3),fpppp(10, 4)/-6.29853616d-05, 6.07418261d-05/
      data fpppp(10, 5),fpppp(10, 6)/ 3.70752246d-06, 7.74527250d-06/
      data fpppp(10, 7),fpppp(10, 8)/-2.80904585d-05, 1.00024216d-04/
      data fpppp(10, 9),fpppp(10,10)/-1.63184051d-04, 1.35023612d-04/
      data fpppp(10,11),fpppp(10,12)/-4.87552866d-05, 7.34953845d-06/
      data fpppp(10,13),fpppp(10,14)/ 9.10487858d-06, 3.93221514d-05/
      data fpppp(10,15),fpppp(10,16)/-2.69421439d-05, 1.42330191d-04/
      data fpppp(10,17),fpppp(10,18)/-5.73871848d-05,-9.84297855d-05/
      data fpppp(10,19) /            -1.84808036d-04 /
      data fpppp(11, 1),fpppp(11, 2)/-9.32598394d-05,-4.60170646d-05/
      data fpppp(11, 3),fpppp(11, 4)/ 5.11497350d-05,-2.64175522d-05/
      data fpppp(11, 5),fpppp(11, 6)/ 1.85109290d-05,-3.41602336d-06/
      data fpppp(11, 7),fpppp(11, 8)/ 2.23671844d-06,-1.19729546d-05/
      data fpppp(11, 9),fpppp(11,10)/ 2.58765763d-05,-1.50180582d-05/
      data fpppp(11,11),fpppp(11,12)/ 1.26620150d-05,-4.87331986d-06/
      data fpppp(11,13),fpppp(11,14)/ 4.78120799d-06, 1.53885631d-07/
      data fpppp(11,15),fpppp(11,16)/ 4.03087619d-05,-3.04901689d-05/
      data fpppp(11,17),fpppp(11,18)/ 3.09127774d-05,-1.75573870d-05/
      data fpppp(11,19) /            -3.46352341d-05 /
      data fpppp(12, 1),fpppp(12, 2)/ 2.14238717d-05, 1.07378002d-05/
      data fpppp(12, 3),fpppp(12, 4)/-7.56162030d-06, 7.85680535d-07/
      data fpppp(12, 5),fpppp(12, 6)/-6.00038835d-06,-1.35787716d-06/
      data fpppp(12, 7),fpppp(12, 8)/ 1.01955271d-05,-9.78346886d-06/
      data fpppp(12, 9),fpppp(12,10)/ 8.81408841d-06,-5.57367899d-06/
      data fpppp(12,11),fpppp(12,12)/-1.01991837d-06, 5.22662105d-06/
      data fpppp(12,13),fpppp(12,14)/ 7.25914408d-07, 5.09292621d-06/
      data fpppp(12,15),fpppp(12,16)/ 2.37299106d-06, 1.66402859d-05/
      data fpppp(12,17),fpppp(12,18)/-3.34490250d-05, 2.11793380d-06/
      data fpppp(12,19) /             1.32116727d-05 /
      data fpppp(13, 1),fpppp(13, 2)/-1.21009223d-06,-8.84074701d-07/
      data fpppp(13, 3),fpppp(13, 4)/-1.05278405d-06, 5.10050680d-07/
      data fpppp(13, 5),fpppp(13, 6)/ 6.11213235d-07, 6.73276307d-07/
      data fpppp(13, 7),fpppp(13, 8)/-2.56098576d-06, 3.71743184d-06/
      data fpppp(13, 9),fpppp(13,10)/-2.20269986d-06, 1.31410410d-06/
      data fpppp(13,11),fpppp(13,12)/-6.25258098d-07,-1.37921837d-06/
      data fpppp(13,13),fpppp(13,14)/ 7.33719115d-07, 4.24028362d-07/
      data fpppp(13,15),fpppp(13,16)/ 6.05580327d-07,-2.45526095d-06/
      data fpppp(13,17),fpppp(13,18)/ 5.20570224d-06,-1.39568385d-06/
      data fpppp(13,19) /            -3.92211301d-06 /
      data fpppp(14, 1),fpppp(14, 2)/-1.45870128d-07,-4.25407415d-09/
      data fpppp(14, 3),fpppp(14, 4)/ 5.65685611d-07,-3.63507415d-07/
      data fpppp(14, 5),fpppp(14, 6)/ 6.90915999d-08,-2.20523742d-07/
      data fpppp(14, 7),fpppp(14, 8)/ 6.59190219d-07,-1.08091357d-06/
      data fpppp(14, 9),fpppp(14,10)/ 5.17468842d-07,-4.19774165d-07/
      data fpppp(14,11),fpppp(14,12)/ 2.62525467d-07, 3.17477995d-07/
      data fpppp(14,13),fpppp(14,14)/-1.94440161d-07,-9.33791520d-08/
      data fpppp(14,15),fpppp(14,16)/-5.05870562d-08, 7.08873006d-07/
      data fpppp(14,17),fpppp(14,18)/-1.37917621d-06, 3.05179511d-07/
      data fpppp(14,19) /             9.11705249d-07 /
      data fpppp(15, 1),fpppp(15, 2)/-8.39283002d-08,-5.69069187d-08/
      data fpppp(15, 3),fpppp(15, 4)/-5.44656826d-08, 3.40060626d-08/
      data fpppp(15, 5),fpppp(15, 6)/ 3.98193188d-08, 3.71957647d-08/
      data fpppp(15, 7),fpppp(15, 8)/-1.63682477d-07, 2.29474818d-07/
      data fpppp(15, 9),fpppp(15,10)/-1.35277579d-07, 7.71484946d-08/
      data fpppp(15,11),fpppp(15,12)/-4.03654254d-08,-8.77629239d-08/
      data fpppp(15,13),fpppp(15,14)/ 4.28404488d-08, 2.53618890d-08/
      data fpppp(15,15),fpppp(15,16)/ 4.14744079d-08,-1.46930761d-07/
      data fpppp(15,17),fpppp(15,18)/ 3.28094810d-07,-9.07033991d-08/
      data fpppp(15,19) /            -2.54123385d-07 /
      data fpppp(16, 1),fpppp(16, 2)/ 3.77469836d-07, 2.06934944d-07/
      data fpppp(16, 3),fpppp(16, 4)/-1.61670212d-07, 1.62970553d-07/
      data fpppp(16, 5),fpppp(16, 6)/-1.92636659d-07, 4.22247032d-10/
      data fpppp(16, 7),fpppp(16, 8)/ 8.84234351d-08,-5.08006106d-08/
      data fpppp(16, 9),fpppp(16,10)/ 1.25645114d-07, 2.38508603d-09/
      data fpppp(16,11),fpppp(16,12)/-5.69538020d-08, 7.71467587d-08/
      data fpppp(16,13),fpppp(16,14)/-2.97091825d-08,-3.04046947d-08/
      data fpppp(16,15),fpppp(16,16)/-9.55175310d-08, 4.72818673d-08/
      data fpppp(16,17),fpppp(16,18)/-1.59583025d-07, 8.63348359d-08/
      data fpppp(16,19) /             2.02886196d-07 /
 

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
      subroutine b_2(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  1.79999000d-02 /
      data f( 1, 3),f( 1, 4) /  1.46198000d-02 ,  1.38823000d-02 /
      data f( 1, 5),f( 1, 6) /  1.99302000d-02 ,  1.05835000d-02 /
      data f( 1, 7),f( 1, 8) /  7.58890000d-03 ,  4.17020000d-03 /
      data f( 1, 9),f( 1,10) / -8.06500000d-04 , -6.35480000d-03 /
      data f( 1,11),f( 1,12) / -1.49988000d-02 , -3.21409000d-02 /
      data f( 1,13),f( 1,14) / -5.65160000d-02 , -1.06297300d-01 /
      data f( 1,15),f( 1,16) / -4.17400200d-01 , -1.08590860d+00 /
      data f( 1,17),f( 1,18) / -6.53825600d-01 , -1.50062700d-01 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  1.74386000d-02 /
      data f( 2, 3),f( 2, 4) /  1.65884000d-02 ,  1.59518000d-02 /
      data f( 2, 5),f( 2, 6) /  1.49906000d-02 ,  1.17824000d-02 /
      data f( 2, 7),f( 2, 8) /  8.42710000d-03 ,  4.58110000d-03 /
      data f( 2, 9),f( 2,10) /  6.78000000d-04 , -5.12710000d-03 /
      data f( 2,11),f( 2,12) / -1.19288000d-02 , -2.56439000d-02 /
      data f( 2,13),f( 2,14) / -4.71031000d-02 , -1.02732300d-01 /
      data f( 2,15),f( 2,16) / -2.56093600d-01 , -5.57562200d-01 /
      data f( 2,17),f( 2,18) / -3.68118000d-01 , -1.14965500d-01 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  7.23660000d-03 /
      data f( 3, 3),f( 3, 4) /  1.89650000d-02 ,  1.76233000d-02 /
      data f( 3, 5),f( 3, 6) /  4.57910000d-03 ,  1.27956000d-02 /
      data f( 3, 7),f( 3, 8) /  8.56260000d-03 ,  4.36080000d-03 /
      data f( 3, 9),f( 3,10) /  1.03550000d-03 , -4.07060000d-03 /
      data f( 3,11),f( 3,12) / -9.32980000d-03 , -2.00560000d-02 /
      data f( 3,13),f( 3,14) / -3.67868000d-02 , -6.67873000d-02 /
      data f( 3,15),f( 3,16) / -1.49717900d-01 , -3.03692100d-01 /
      data f( 3,17),f( 3,18) / -2.11004600d-01 , -7.27466000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  2.43330000d-03 /
      data f( 4, 3),f( 4, 4) /  2.03023000d-02 ,  1.79158000d-02 /
      data f( 4, 5),f( 4, 6) /  1.02480000d-02 ,  1.24992000d-02 /
      data f( 4, 7),f( 4, 8) /  5.44210000d-03 ,  3.83230000d-03 /
      data f( 4, 9),f( 4,10) /  1.06310000d-03 , -3.18160000d-03 /
      data f( 4,11),f( 4,12) / -7.17980000d-03 , -1.52555000d-02 /
      data f( 4,13),f( 4,14) / -2.83630000d-02 , -5.15205000d-02 /
      data f( 4,15),f( 4,16) / -9.71736000d-02 , -1.75637900d-01 /
      data f( 4,17),f( 4,18) / -1.28907100d-01 , -5.04433000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  3.20780000d-03 /
      data f( 5, 3),f( 5, 4) /  1.78254000d-02 ,  1.67442000d-02 /
      data f( 5, 5),f( 5, 6) /  1.20961000d-02 ,  1.14910000d-02 /
      data f( 5, 7),f( 5, 8) /  7.49810000d-03 ,  3.13940000d-03 /
      data f( 5, 9),f( 5,10) /  9.72100000d-04 , -2.45850000d-03 /
      data f( 5,11),f( 5,12) / -5.44720000d-03 , -1.17041000d-02 /
      data f( 5,13),f( 5,14) / -2.16514000d-02 , -3.73691000d-02 /
      data f( 5,15),f( 5,16) / -6.63813000d-02 , -1.08096200d-01 /
      data f( 5,17),f( 5,18) / -8.56091000d-02 , -3.87668000d-02 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  9.85160000d-03 /
      data f( 6, 3),f( 6, 4) /  1.46710000d-02 ,  1.57673000d-02 /
      data f( 6, 5),f( 6, 6) /  1.19071000d-02 ,  1.04514000d-02 /
      data f( 6, 7),f( 6, 8) /  5.96000000d-03 ,  2.47990000d-03 /
      data f( 6, 9),f( 6,10) /  8.28200000d-04 , -1.88670000d-03 /
      data f( 6,11),f( 6,12) / -4.20280000d-03 , -8.84800000d-03 /
      data f( 6,13),f( 6,14) / -1.59454000d-02 , -2.83665000d-02 /
      data f( 6,15),f( 6,16) / -4.55138000d-02 , -7.05952000d-02 /
      data f( 6,17),f( 6,18) / -6.45013000d-02 , -3.39788000d-02 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  7.04990000d-03 /
      data f( 7, 3),f( 7, 4) /  1.15369000d-02 ,  1.15291000d-02 /
      data f( 7, 5),f( 7, 6) /  1.12225000d-02 ,  8.88620000d-03 /
      data f( 7, 7),f( 7, 8) /  4.76620000d-03 ,  1.92110000d-03 /
      data f( 7, 9),f( 7,10) /  6.62800000d-04 , -1.44190000d-03 /
      data f( 7,11),f( 7,12) / -3.10360000d-03 , -6.62220000d-03 /
      data f( 7,13),f( 7,14) / -1.21022000d-02 , -2.11180000d-02 /
      data f( 7,15),f( 7,16) / -3.25488000d-02 , -4.89578000d-02 /
      data f( 7,17),f( 7,18) / -5.46940000d-02 , -3.98831000d-02 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  5.04510000d-03 /
      data f( 8, 3),f( 8, 4) /  8.33070000d-03 ,  9.90990000d-03 /
      data f( 8, 5),f( 8, 6) /  9.46120000d-03 ,  6.81090000d-03 /
      data f( 8, 7),f( 8, 8) /  3.46280000d-03 ,  1.47880000d-03 /
      data f( 8, 9),f( 8,10) /  4.91200000d-04 , -1.10150000d-03 /
      data f( 8,11),f( 8,12) / -2.29530000d-03 , -4.91280000d-03 /
      data f( 8,13),f( 8,14) / -9.41070000d-03 , -1.60634000d-02 /
      data f( 8,15),f( 8,16) / -2.40823000d-02 , -3.53424000d-02 /
      data f( 8,17),f( 8,18) / -4.10604000d-02 , -3.00760000d-02 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  2.55250000d-03 /
      data f( 9, 3),f( 9, 4) /  4.47910000d-03 ,  5.42460000d-03 /
      data f( 9, 5),f( 9, 6) /  5.17280000d-03 ,  3.55510000d-03 /
      data f( 9, 7),f( 9, 8) /  1.75570000d-03 ,  9.04700000d-04 /
      data f( 9, 9),f( 9,10) /  1.48000000d-04 , -5.11900000d-04 /
      data f( 9,11),f( 9,12) / -1.25010000d-03 , -2.64860000d-03 /
      data f( 9,13),f( 9,14) / -5.16930000d-03 , -9.13860000d-03 /
      data f( 9,15),f( 9,16) / -1.42720000d-02 , -2.04603000d-02 /
      data f( 9,17),f( 9,18) / -2.28765000d-02 , -1.61389000d-02 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.28460000d-03 /
      data f(10, 3),f(10, 4) /  2.09160000d-03 ,  2.71810000d-03 /
      data f(10, 5),f(10, 6) /  2.56540000d-03 ,  1.75490000d-03 /
      data f(10, 7),f(10, 8) /  8.79700000d-04 ,  4.13500000d-04 /
      data f(10, 9),f(10,10) / -4.18000000d-05 , -3.05100000d-04 /
      data f(10,11),f(10,12) / -6.81600000d-04 , -1.39510000d-03 /
      data f(10,13),f(10,14) / -2.72270000d-03 , -4.92070000d-03 /
      data f(10,15),f(10,16) / -7.88330000d-03 , -1.09987000d-02 /
      data f(10,17),f(10,18) / -1.20714000d-02 , -8.12730000d-03 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  4.73100000d-04 /
      data f(11, 3),f(11, 4) /  8.94000000d-04 ,  1.13320000d-03 /
      data f(11, 5),f(11, 6) /  1.05640000d-03 ,  7.27300000d-04 /
      data f(11, 7),f(11, 8) /  3.67500000d-04 ,  1.10900000d-04 /
      data f(11, 9),f(11,10) /  5.04000000d-05 , -1.61700000d-04 /
      data f(11,11),f(11,12) / -3.31700000d-04 , -6.15500000d-04 /
      data f(11,13),f(11,14) / -1.16190000d-03 , -2.11440000d-03 /
      data f(11,15),f(11,16) / -3.43880000d-03 , -4.72360000d-03 /
      data f(11,17),f(11,18) / -5.04590000d-03 , -3.08820000d-03 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  2.06100000d-04 /
      data f(12, 3),f(12, 4) /  3.83900000d-04 ,  4.78100000d-04 /
      data f(12, 5),f(12, 6) /  4.41100000d-04 ,  3.02300000d-04 /
      data f(12, 7),f(12, 8) /  1.48800000d-04 ,  1.42700000d-04 /
      data f(12, 9),f(12,10) /  7.01000000d-05 ,  8.20000000d-06 /
      data f(12,11),f(12,12) / -1.17400000d-04 , -2.90100000d-04 /
      data f(12,13),f(12,14) / -4.93200000d-04 , -8.66400000d-04 /
      data f(12,15),f(12,16) / -1.39880000d-03 , -1.86520000d-03 /
      data f(12,17),f(12,18) / -1.89260000d-03 , -1.21240000d-03 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 , -3.49000000d-05 /
      data f(13, 3),f(13, 4) /  1.31000000d-05 ,  3.62000000d-05 /
      data f(13, 5),f(13, 6) /  7.07000000d-05 ,  7.34000000d-05 /
      data f(13, 7),f(13, 8) /  6.78000000d-05 ,  5.73000000d-05 /
      data f(13, 9),f(13,10) /  4.88000000d-05 ,  4.98000000d-05 /
      data f(13,11),f(13,12) /  4.06000000d-05 , -1.93000000d-05 /
      data f(13,13),f(13,14) / -9.11000000d-05 , -1.63800000d-04 /
      data f(13,15),f(13,16) / -2.34400000d-04 , -3.00400000d-04 /
      data f(13,17),f(13,18) / -3.12100000d-04 , -2.09400000d-04 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -1.90000000d-06 /
      data f(14, 3),f(14, 4) / -2.30000000d-06 , -3.00000000d-07 /
      data f(14, 5),f(14, 6) /  2.80000000d-06 ,  6.40000000d-06 /
      data f(14, 7),f(14, 8) /  1.08000000d-05 ,  1.01000000d-05 /
      data f(14, 9),f(14,10) /  1.10000000d-05 ,  1.10000000d-05 /
      data f(14,11),f(14,12) /  5.70000000d-06 , -6.30000000d-06 /
      data f(14,13),f(14,14) / -1.97000000d-05 , -3.07000000d-05 /
      data f(14,15),f(14,16) / -3.66000000d-05 , -3.78000000d-05 /
      data f(14,17),f(14,18) / -3.49000000d-05 , -2.00000000d-05 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -1.50000000d-06 /
      data f(15, 3),f(15, 4) / -2.20000000d-06 , -1.90000000d-06 /
      data f(15, 5),f(15, 6) / -9.00000000d-07 ,  8.00000000d-07 /
      data f(15, 7),f(15, 8) /  2.60000000d-06 ,  1.50000000d-06 /
      data f(15, 9),f(15,10) /  1.90000000d-06 ,  0.00000000d+00 /
      data f(15,11),f(15,12) / -3.20000000d-06 , -8.00000000d-06 /
      data f(15,13),f(15,14) / -1.26000000d-05 , -1.62000000d-05 /
      data f(15,15),f(15,16) / -1.84000000d-05 , -1.93000000d-05 /
      data f(15,17),f(15,18) / -1.81000000d-05 , -1.06000000d-05 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 ,  0.00000000d+00 /
      data f(16, 3),f(16, 4) /  0.00000000d+00 ,  1.00000000d-07 /
      data f(16, 5),f(16, 6) /  9.00000000d-07 ,  1.50000000d-06 /
      data f(16, 7),f(16, 8) /  1.30000000d-06 ,  4.00000000d-07 /
      data f(16, 9),f(16,10) / -1.20000000d-06 , -2.20000000d-06 /
      data f(16,11),f(16,12) / -2.40000000d-06 , -3.20000000d-06 /
      data f(16,13),f(16,14) / -3.40000000d-06 , -3.00000000d-06 /
      data f(16,15),f(16,16) / -3.00000000d-06 , -3.20000000d-06 /
      data f(16,17),f(16,18) / -2.80000000d-06 , -2.60000000d-06 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-4.50254363d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-5.99847153d-01,-2.21343273d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 4.95581139d-02, 5.28274562d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 1.47486347d-02, 1.68589448d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-6.26656235d-01,-3.20061249d-04/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.17865961d-02, 1.87979549d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 2.90656847d-02,-5.07309453d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-2.36456743d-02,-1.05017674d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-4.75335278d-02,-7.41985254d-07/
      data fpp( 1,10,1),fpp( 1,10,2)/-4.36800707d-03,-2.08262916d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-1.24051259d-02,-1.01694848d-04/
      data fpp( 1,12,1),fpp( 1,12,2)/-2.48460699d-02,-8.22803154d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 8.96906545d-02,-3.16389010d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.05350872d+00,-1.42943612d-03/
      data fpp( 1,15,1),fpp( 1,15,2)/-1.43772271d+00,-9.95838761d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/-1.05005973d+01, 1.98186566d-02/
      data fpp( 1,17,1),fpp( 1,17,2)/-4.54113130d+00,-3.28075469d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/ 8.16198347d-01,-2.39484380d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-8.36188210d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-3.60316073d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-2.75328195d-01,-1.98079855d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.84877217d-03, 5.53074918d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/-1.02972695d-02,-1.03341123d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-2.34695030d-01,-3.34470427d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-1.30569223d-03, 9.30228305d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.68311307d-02,-1.25880895d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-1.61986514d-02, 1.16080750d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-2.93129445d-02,-3.72702105d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-4.28898585d-03, 2.33527671d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.16147483d-02,-1.15936858d-04/
      data fpp( 2,12,1),fpp( 2,12,2)/-2.45753602d-02, 2.55906647d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.81686910d-02,-4.51071801d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/ 9.74392570d-01,-2.71503461d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/-1.29934708d+00,-4.32684035d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/-7.01754794d+00, 8.69242688d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/-3.24121740d+00,-9.88099160d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 2.53595806d-01,-9.17532240d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-3.63299188d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 2.18574155d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 2.55054931d-01, 4.88436896d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.75320260d-03,-1.44440914d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-3.32595568d-02,-2.55286035d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 7.44651355d-01, 4.63435054d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-4.44188272d-02,-3.22812183d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-2.01795207d-01, 8.08436756d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-6.23972013d-03, 1.30947994d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-4.26469432d-03,-3.34915954d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-4.15604952d-03, 2.58089016d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.17858810d-02,-7.89300112d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-1.32174893d-02,-3.81088569d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/-6.68554184d-02,-1.28910561d-04/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.09407899d+00,-2.42430898d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/-1.60452397d+00,-2.07717185d-03/
      data fpp( 3,16,1),fpp( 3,16,2)/-2.60065595d+00, 4.28850229d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.78312911d+00,-2.77135309d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/-7.62326570d-01,-4.45731055d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-1.87062447d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 4.91636579d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 6.49134703d-02, 1.93621842d-04/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.44730962d-01,-3.39981946d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-6.35145034d-02,-4.90240592d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-3.31850390d-01, 2.19200182d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-1.74589990d-02,-2.32636670d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 3.01949699d-01, 1.52848497d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-5.07246810d-03,-5.19193189d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-3.11327823d-03,-1.47352217d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/-4.21181608d-03, 2.23302056d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-8.59172778d-03,-5.97956009d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/-4.06646826d-02,-2.77978022d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/-4.46220172d-02,-1.30921190d-04/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.00193404d-01,-5.15174360d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/-3.57267050d-01,-1.01274507d-03/
      data fpp( 4,16,1),fpp( 4,16,2)/-1.45221324d+00, 2.13382570d-03/
      data fpp( 4,17,1),fpp( 4,17,2)/-8.78651180d-01,-1.08517264d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/-1.91629526d-01,-1.86438792d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-9.24623104d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 3.68943294d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 3.21961188d-01, 1.46579412d-04/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.05470497d-02,-2.70672942d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 6.77025702d-02,-5.81564438d-06/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 9.63020601d-03, 7.99215194d-05/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 7.48482310d-03,-7.12904333d-05/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-2.29528588d-01, 1.97221381d-06/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.86959252d-03, 4.14535781d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-1.07219276d-03,-3.63025261d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/-3.88168617d-03, 2.79585263d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.64572079d-02,-4.90175789d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.14887803d-02,-2.79802106d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/-1.14865126d-02,-6.04855788d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/-2.74004623d-01,-7.63014741d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/-2.29207833d-01,-4.31978525d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/-6.67366072d-01, 1.04205357d-03/
      data fpp( 5,17,1),fpp( 5,17,2)/-5.22191175d-01, 1.15884230d-04/
      data fpp( 5,18,1),fpp( 5,18,2)/-6.51753245d-02,-4.42784943d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-4.23300253d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-6.09439780d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-4.72363221d-01,-5.52600440d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 9.17762921d-04,-1.99478460d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-1.78090778d-01,-8.83345719d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-1.22354338d-02, 7.58961337d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-1.71902934d-02,-7.09799628d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 7.70496549d-02, 2.58817176d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 2.60409802d-03, 2.81310923d-05/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-5.32950716d-04,-2.87020869d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/-2.95643924d-03, 2.28852553d-05/
      data fpp( 6,11,1),fpp( 6,11,2)/ 1.19055938d-03,-3.89109344d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/-1.76751962d-02,-6.98751788d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-6.02719324d-02,-8.02709941d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/ 2.35050863d-02, 8.64949435d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/-2.14621617d-01,-2.37898983d-04/
      data fpp( 6,16,1),fpp( 6,16,2)/-3.84427467d-01, 4.66900439d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/-3.61114121d-01, 2.40815228d-04/
      data fpp( 6,18,1),fpp( 6,18,2)/-5.80944176d-01, 3.55546491d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-1.75655825d-04/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-9.54129033d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.50666697d-01,-1.91664193d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/-1.11731014d-02,-6.75670323d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.55465540d-01, 1.97465487d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-3.50284709d-02,-2.93471623d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.75636494d-02,-2.41398993d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-2.70250310d-02, 1.88847596d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 2.81901541d-03, 2.50948609d-05/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-2.10043725d-05,-2.40562033d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/-3.34255686d-03, 2.03459522d-05/
      data fpp( 7,11,1),fpp( 7,11,2)/-1.00850296d-02,-3.07476055d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.23554347d-02,-8.76953005d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-2.68457579d-02,-5.18582743d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/-8.31307227d-02, 4.05462716d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/-9.76806986d-02,-1.09260234d-04/
      data fpp( 7,16,1),fpp( 7,16,2)/-1.74464058d-01, 1.34294310d-04/
      data fpp( 7,17,1),fpp( 7,17,2)/ 2.71572659d-01, 2.12450993d-04/
      data fpp( 7,18,1),fpp( 7,18,2)/ 7.85107027d-01, 2.48727716d-04/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 2.96970142d-04/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-1.79958261d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.07685673d-02,-1.78553478d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 3.29596427d-02,-1.61527828d-05/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-5.09213820d-02,-1.99175210d-05/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-9.15568256d-03,-2.58511332d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 1.09298909d-02,-8.77394609d-06/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 1.46104693d-02, 1.90789176d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 3.59484034d-03, 1.43042758d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-3.13031794d-04,-1.65120207d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/ 6.66666683d-04, 1.54378071d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/-4.48544085d-03,-2.13052077d-05/
      data fpp( 8,12,1),fpp( 8,12,2)/-1.03630648d-02,-1.56389762d-05/
      data fpp( 8,13,1),fpp( 8,13,2)/-5.10003611d-03,-2.89628876d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/-2.00671954d-02, 2.20252650d-06/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.94305886d-02,-6.18192184d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/-1.21016299d-01, 5.06023472d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/-1.51231516d-01, 1.91935830d-04/
      data fpp( 8,18,1),fpp( 8,18,2)/-2.02773932d-01, 1.83798334d-04/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 2.18366833d-04/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-2.79528442d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.38598534d-02,-6.08243116d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 2.73762268d-03,-1.04289910d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 2.82726260d-02,-1.10676050d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.62637831d-02,-1.71385889d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 9.54715198d-03,-2.33203937d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 3.41985766d-03, 1.55647464d-05/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-5.50278717d-04,-3.02294621d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 9.49597568d-04, 2.18503843d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/-3.74872162d-03, 9.07925022d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/-2.92866263d-03,-7.24620843d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-6.03058821d-03,-1.07239588d-05/
      data fpp( 9,13,1),fpp( 9,13,2)/-1.40870127d-02,-1.71899565d-05/
      data fpp( 9,14,1),fpp( 9,14,2)/-1.76480525d-02,-7.43221521d-06/
      data fpp( 9,15,1),fpp( 9,15,2)/-9.96913482d-03,-2.29271827d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/-1.27953224d-02, 3.58469459d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/-2.27155317d-02, 1.05865399d-04/
      data fpp( 9,18,1),fpp( 9,18,2)/ 2.87703336d-03, 8.99194574d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 9.85347713d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00,-7.30841020d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.25540355d-03,-5.65317961d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 1.09936166d-02, 1.26512862d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/ 4.53587793d-03,-1.02373349d-05/
      data fpp(10, 5,1),fpp(10, 5,2)/ 7.13805004d-03,-7.06778906d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 5.46650118d-03,-9.59508853d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.87635008d-03, 7.02382448d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.71502453d-03,-2.59578905d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/ 2.26714152d-03, 4.01333174d-06/
      data fpp(10,10,1),fpp(10,10,2)/-2.67802088d-05,-1.93753790d-06/
      data fpp(10,11,1),fpp(10,11,2)/-1.67615863d-03,-3.05518015d-06/
      data fpp(10,12,1),fpp(10,12,2)/-3.41583234d-03,-6.06174149d-06/
      data fpp(10,13,1),fpp(10,13,2)/-5.85691294d-03,-9.54385388d-06/
      data fpp(10,14,1),fpp(10,14,2)/-1.08493445d-02,-7.98684300d-06/
      data fpp(10,15,1),fpp(10,15,2)/-1.90028721d-02,-4.38477411d-06/
      data fpp(10,16,1),fpp(10,16,2)/-3.10711611d-02, 1.63579395d-05/
      data fpp(10,17,1),fpp(10,17,2)/-3.46113573d-02, 6.15150163d-05/
      data fpp(10,18,1),fpp(10,18,2)/-3.09404512d-02, 3.85899953d-05/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 3.51170023d-05/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 7.91320689d-07/
      data fpp(11, 2,1),fpp(11, 2,2)/ 2.95366448d-03,-5.58641378d-07/
      data fpp(11, 3,1),fpp(11, 3,2)/ 1.11548201d-03,-1.68875518d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 4.21013864d-03,-3.58833791d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 3.29799336d-03,-2.91789318d-06/
      data fpp(11, 6,1),fpp(11, 6,2)/ 2.02647417d-03, 1.21910616d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 8.96453602d-04, 5.88250712d-07/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.73973466d-03, 3.71708653d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.01458753d-03,-3.69059685d-06/
      data fpp(11,10,1),fpp(11,10,2)/ 3.32986046d-04, 1.94930086d-06/
      data fpp(11,11,1),fpp(11,11,2)/-2.80298815d-04,-1.58060660d-06/
      data fpp(11,12,1),fpp(11,12,2)/-1.77313301d-03,-2.45487445d-06/
      data fpp(11,13,1),fpp(11,13,2)/-3.58430322d-03,-4.35589560d-06/
      data fpp(11,14,1),fpp(11,14,2)/-6.00971785d-03,-4.48754315d-06/
      data fpp(11,15,1),fpp(11,15,2)/-8.60735263d-03,-7.93180850d-09/
      data fpp(11,16,1),fpp(11,16,2)/-1.11531623d-02, 6.89527038d-06/
      data fpp(11,17,1),fpp(11,17,2)/-1.27676883d-02, 3.01768503d-05/
      data fpp(11,18,1),fpp(11,18,2)/-1.03256025d-02, 9.19732849d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 8.63835755d-07/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 2.75238002d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/-2.06145536d-06,-2.93476003d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.04445534d-03,-7.99333989d-07/
      data fpp(12, 4,1),fpp(12, 4,2)/ 9.38767507d-04,-1.52518804d-06/
      data fpp(12, 5,1),fpp(12, 5,2)/ 1.11877652d-03,-9.71913849d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 8.90002133d-04,-6.95156563d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 5.81835517d-04, 2.87054010d-06/
      data fpp(12, 8,1),fpp(12, 8,2)/-6.48363173d-04,-1.94300385d-06/
      data fpp(12, 9,1),fpp(12, 9,2)/ 5.12085935d-05, 9.11475281d-07/
      data fpp(12,10,1),fpp(12,10,2)/-6.69163977d-04,-1.06089728d-06/
      data fpp(12,11,1),fpp(12,11,2)/-4.57046107d-04,-4.89886166d-07/
      data fpp(12,12,1),fpp(12,12,2)/-3.92435633d-04, 1.94441943d-07/
      data fpp(12,13,1),fpp(12,13,2)/-1.21627416d-03,-2.11188161d-06/
      data fpp(12,14,1),fpp(12,14,2)/-2.51098412d-03,-1.95291551d-06/
      data fpp(12,15,1),fpp(12,15,2)/-4.27571741d-03, 3.71543660d-07/
      data fpp(12,16,1),fpp(12,16,2)/-6.31698967d-03, 4.42674087d-06/
      data fpp(12,17,1),fpp(12,17,2)/-7.25068954d-03, 8.26149284d-06/
      data fpp(12,18,1),fpp(12,18,2)/-3.67633883d-03, 4.98328776d-06/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 3.73735612d-06/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 1.83628765d-06/
      data fpp(13, 2,1),fpp(13, 2,2)/ 2.87352127d-04, 9.70424703d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 2.05292988d-04,-7.43986460d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/ 2.88428158d-04, 5.11521136d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/ 1.55873772d-04,-6.18098084d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/ 4.33565162d-05, 5.28711994d-08/
      data fpp(13, 7,1),fpp(13, 7,2)/-5.53333517d-05,-9.13867136d-08/
      data fpp(13, 8,1),fpp(13, 8,2)/ 1.81222188d-04, 1.86756551d-08/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.05320163d-05, 1.36684093d-07/
      data fpp(13,10,1),fpp(13,10,2)/ 5.17989067d-05, 4.58797221d-09/
      data fpp(13,11,1),fpp(13,11,2)/-1.12312271d-04,-7.67035982d-07/
      data fpp(13,12,1),fpp(13,12,2)/-2.16126599d-04, 2.15559559d-08/
      data fpp(13,13,1),fpp(13,13,2)/-1.70825904d-04,-3.31878414d-08/
      data fpp(13,14,1),fpp(13,14,2)/-2.22588725d-04, 5.71954098d-08/
      data fpp(13,15,1),fpp(13,15,2)/-3.62771462d-04,-6.95937976d-08/
      data fpp(13,16,1),fpp(13,16,2)/-3.84449819d-04, 4.97179781d-07/
      data fpp(13,17,1),fpp(13,17,2)/-2.20687252d-04, 1.33887467d-06/
      data fpp(13,18,1),fpp(13,18,2)/-2.99782254d-04, 1.01132152d-06/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 1.01783924d-06/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00, 6.88161030d-09/
      data fpp(14, 2,1),fpp(14, 2,2)/-8.85256544d-05, 1.32367794d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-4.88066306d-05, 3.01712721d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-6.37182272d-05, 1.00781323d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.76595733d-05,-4.48380127d-09/
      data fpp(14, 6,1),fpp(14, 6,2)/ 1.11293852d-05, 3.78570728d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/ 3.25822965d-05,-9.89444899d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/-3.40849786d-05, 5.19208869d-08/
      data fpp(14, 9,1),fpp(14, 9,2)/ 1.31917521d-05,-1.27390576d-08/
      data fpp(14,10,1),fpp(14,10,2)/-3.81473193d-06,-5.49646566d-08/
      data fpp(14,11,1),fpp(14,11,2)/ 3.91098666d-05,-8.54023162d-08/
      data fpp(14,12,1),fpp(14,12,2)/ 5.16976121d-05,-5.42607864d-09/
      data fpp(14,13,1),fpp(14,13,2)/ 2.14147938d-05, 2.31066307d-08/
      data fpp(14,14,1),fpp(14,14,2)/ 1.51082327d-05, 5.69995556d-08/
      data fpp(14,15,1),fpp(14,15,2)/ 2.96730891d-05, 5.48951467d-08/
      data fpp(14,16,1),fpp(14,16,2)/ 1.13442934d-05, 5.41985759d-09/
      data fpp(14,17,1),fpp(14,17,2)/-3.82934754d-05, 1.69425423d-07/
      data fpp(14,18,1),fpp(14,18,2)/ 1.26161792d-05, 3.68784506d-08/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00,-1.09392253d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00, 6.21797591d-09/
      data fpp(15, 2,1),fpp(15, 2,2)/ 1.78504901d-05, 7.56404818d-09/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.31835349d-05, 1.15258314d-08/
      data fpp(15, 4,1),fpp(15, 4,2)/ 1.87947509d-05, 6.33262629d-09/
      data fpp(15, 5,1),fpp(15, 5,2)/ 1.10645215d-05, 5.14366347d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/ 4.22594297d-06, 1.50927198d-08/
      data fpp(15, 7,1),fpp(15, 7,2)/-1.79583448d-06,-5.95145428d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/ 1.30177261d-05, 4.89654514d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/ 8.15007911d-07,-4.63472627d-08/
      data fpp(15,10,1),fpp(15,10,2)/ 5.16002096d-06,-1.57640071d-09/
      data fpp(15,11,1),fpp(15,11,2)/-5.12719551d-06,-2.53471345d-08/
      data fpp(15,12,1),fpp(15,12,2)/-1.27138499d-05, 6.96493874d-09/
      data fpp(15,13,1),fpp(15,13,2)/-1.12832708d-05, 9.48737957d-09/
      data fpp(15,14,1),fpp(15,14,2)/-1.57442060d-05, 1.50855430d-08/
      data fpp(15,15,1),fpp(15,15,2)/-2.53208949d-05, 1.41704484d-08/
      data fpp(15,16,1),fpp(15,16,2)/-2.70773547d-05, 6.23266331d-09/
      data fpp(15,17,1),fpp(15,17,2)/-1.67388461d-05, 8.68988983d-08/
      data fpp(15,18,1),fpp(15,18,2)/-2.06824623d-05, 2.41717433d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00, 2.41412833d-09/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00,-1.37618086d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/ 1.15469235d-07, 7.52361718d-10/
      data fpp(16, 3,1),fpp(16, 3,2)/-1.00406960d-05,-1.63326601d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/-1.72370183d-05, 1.17807023d-08/
      data fpp(16, 5,1),fpp(16, 5,2)/-2.02086893d-05,-3.48954331d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/-1.54394001d-05,-9.82252910d-09/
      data fpp(16, 7,1),fpp(16, 7,2)/-8.40208276d-06,-5.22034029d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/-1.28024345d-05,-1.12961097d-08/
      data fpp(16, 9,1),fpp(16, 9,2)/-4.47786110d-06, 8.40477928d-09/
      data fpp(16,10,1),fpp(16,10,2)/-5.12358191d-06, 1.36769926d-08/
      data fpp(16,11,1),fpp(16,11,2)/ 4.50740612d-07,-1.51127498d-08/
      data fpp(16,12,1),fpp(16,12,2)/ 1.28144250d-05, 1.07740067d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 2.23677068d-05, 8.01672291d-09/
      data fpp(16,14,1),fpp(16,14,2)/ 3.67085316d-05,-6.84089837d-09/
      data fpp(16,15,1),fpp(16,15,2)/ 5.66875903d-05,-4.65312944d-09/
      data fpp(16,16,1),fpp(16,16,2)/ 7.49283202d-05, 1.34534161d-08/
      data fpp(16,17,1),fpp(16,17,2)/ 7.47251374d-05,-1.31605351d-08/
      data fpp(16,18,1),fpp(16,18,2)/ 5.64640883d-05, 2.71887243d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00, 4.84056378d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 3.13305865d-02, 1.34857429d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/-1.03184129d-02,-1.32649761d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/ 2.69825940d-02,-1.72745377d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 3.64573236d-03,-9.07818565d-04/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.71495223d-03,-1.92878791d-03/
      data fpppp( 1,11),fpppp( 1,12)/ 2.92804102d-03,-1.00476057d-02/
      data fpppp( 1,13),fpppp( 1,14)/ 4.48810419d-02,-5.85196818d-02/
      data fpppp( 1,15),fpppp( 1,16)/-1.38105284d-01, 2.76642229d-01/
      data fpppp( 1,17),fpppp( 1,18)/-6.71231990d-02,-4.42776126d-02/
      data fpppp( 1,19) /            -1.26178030d-01 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 1.37492093d-02, 6.09329960d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/-4.85209803d-03,-4.28428799d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/ 9.37414684d-03,-5.74507347d-03/
      data fpppp( 2, 7),fpppp( 2, 8)/ 6.90996132d-04,-8.89073643d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 8.59562664d-04,-1.06104819d-03/
      data fpppp( 2,11),fpppp( 2,12)/ 1.44364683d-03,-5.05163011d-03/
      data fpppp( 2,13),fpppp( 2,14)/ 2.27051534d-02,-3.21601938d-02/
      data fpppp( 2,15),fpppp( 2,16)/-8.72621899d-02, 1.74541281d-01/
      data fpppp( 2,17),fpppp( 2,18)/-4.12310489d-02,-2.65081257d-02/
      data fpppp( 2,19) /            -7.76409890d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.29846332d-02,-4.30992238d-03/
      data fpppp( 3, 3),fpppp( 3, 4)/-6.07461140d-04, 2.04978737d-02/
      data fpppp( 3, 5),fpppp( 3, 6)/-3.29389977d-02, 1.72392516d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.88361951d-03,-3.59781759d-03/
      data fpppp( 3, 9),fpppp( 3,10)/ 8.92823175d-04,-8.54579670d-05/
      data fpppp( 3,11),fpppp( 3,12)/-1.01529988d-03, 4.51855089d-03/
      data fpppp( 3,13),fpppp( 3,14)/-2.01912829d-02, 1.78314421d-02/
      data fpppp( 3,15),fpppp( 3,16)/-2.01277693d-02, 3.35384143d-02/
      data fpppp( 3,17),fpppp( 3,18)/-5.20635799d-03,-5.16441124d-04/
      data fpppp( 3,19) /            -8.23643546d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-7.68484700d-03,-4.17541893d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 7.91304859d-03,-1.00251220d-02/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.12142987d-02, 1.31563856d-04/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.14395157d-02, 8.04064718d-03/
      data fpppp( 4, 9),fpppp( 4,10)/-2.18419158d-03, 5.12655478d-04/
      data fpppp( 4,11),fpppp( 4,12)/-6.33127632d-05,-1.92098701d-03/
      data fpppp( 4,13),fpppp( 4,14)/ 9.43419803d-03,-1.48894397d-02/
      data fpppp( 4,15),fpppp( 4,16)/-1.00129916d-02, 2.86922616d-02/
      data fpppp( 4,17),fpppp( 4,18)/-4.64555948d-03,-3.30244831d-03/
      data fpppp( 4,19) /            -1.18681749d-02 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.56393352d-02,-7.76148913d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/ 8.68277213d-03,-4.85541988d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/ 3.82523430d-03,-7.08989843d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/ 1.04422777d-02,-6.57451677d-03/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.79539144d-03,-5.99111453d-04/
      data fpppp( 5,11),fpppp( 5,12)/ 1.50926784d-05, 1.59137770d-03/
      data fpppp( 5,13),fpppp( 5,14)/-6.67857307d-03, 9.37169194d-03/
      data fpppp( 5,15),fpppp( 5,16)/-1.23693007d-02, 1.11282092d-02/
      data fpppp( 5,17),fpppp( 5,18)/ 2.85645205d-03,-3.84356026d-03/
      data fpppp( 5,19) /            -1.09926425d-02 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 2.41630896d-02, 1.20018215d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.54317233d-02, 1.05877002d-02/
      data fpppp( 6, 5),fpppp( 6, 6)/-6.22724433d-03, 4.07266494d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/-4.11172694d-03, 2.25311251d-03/
      data fpppp( 6, 9),fpppp( 6,10)/-6.22212626d-04, 2.78551603d-04/
      data fpppp( 6,11),fpppp( 6,12)/-9.77645577d-05,-1.26825863d-03/
      data fpppp( 6,13),fpppp( 6,14)/ 3.74694024d-03,-6.13707703d-03/
      data fpppp( 6,15),fpppp( 6,16)/ 1.48714457d-03, 4.28774994d-03/
      data fpppp( 6,17),fpppp( 6,18)/-7.05099254d-03, 9.32761617d-03/
      data fpppp( 6,19) /             1.77869817d-02 /
      data fpppp( 7, 1),fpppp( 7, 2)/-8.81952941d-03,-4.55583475d-03/
      data fpppp( 7, 3),fpppp( 7, 4)/ 8.29247868d-03,-8.90537357d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/ 5.90105647d-03,-2.22132236d-03/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.36866079d-03,-8.94995127d-04/
      data fpppp( 7, 9),fpppp( 7,10)/ 2.50275738d-04,-1.34999786d-04/
      data fpppp( 7,11),fpppp( 7,12)/ 8.44681905d-05, 6.54510845d-05/
      data fpppp( 7,13),fpppp( 7,14)/-1.07946761d-03, 1.74474085d-03/
      data fpppp( 7,15),fpppp( 7,16)/-3.39539646d-03, 8.10284197d-03/
      data fpppp( 7,17),fpppp( 7,18)/ 2.35323323d-03,-1.34659159d-02/
      data fpppp( 7,19) /            -2.64080534d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 2.13644451d-03, 1.00413455d-03/
      data fpppp( 8, 3),fpppp( 8, 4)/-2.88317605d-03, 2.87201557d-03/
      data fpppp( 8, 5),fpppp( 8, 6)/-1.06608278d-03, 9.15079772d-05/
      data fpppp( 8, 7),fpppp( 8, 8)/-2.84248838d-04, 1.63714933d-04/
      data fpppp( 8, 9),fpppp( 8,10)/ 5.58545145d-05,-9.38787546d-05/
      data fpppp( 8,11),fpppp( 8,12)/-4.82478566d-05, 2.43339196d-04/
      data fpppp( 8,13),fpppp( 8,14)/-2.56669768d-04,-4.30471401d-04/
      data fpppp( 8,15),fpppp( 8,16)/-8.52186688d-05, 6.38007020d-04/
      data fpppp( 8,17),fpppp( 8,18)/-1.18457975d-03, 2.82067999d-03/
      data fpppp( 8,19) /             5.16084072d-03 /
      data fpppp( 9, 1),fpppp( 9, 2)/-7.95458319d-04,-3.91332255d-04/
      data fpppp( 9, 3),fpppp( 9, 4)/ 8.61862287d-04,-8.56682846d-04/
      data fpppp( 9, 5),fpppp( 9, 6)/ 3.12238322d-04,-7.47377402d-05/
      data fpppp( 9, 7),fpppp( 9, 8)/ 2.20728481d-05, 1.15875824d-04/
      data fpppp( 9, 9),fpppp( 9,10)/-1.57375386d-04, 1.41733991d-04/
      data fpppp( 9,11),fpppp( 9,12)/-7.84578868d-05,-6.32215180d-05/
      data fpppp( 9,13),fpppp( 9,14)/ 3.40740226d-05, 1.96648511d-04/
      data fpppp( 9,15),fpppp( 9,16)/-1.46270615d-04,-2.41872367d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 6.88118777d-04,-3.79836279d-04/
      data fpppp( 9,19) /            -8.76949564d-04 /
      data fpppp(10, 1),fpppp(10, 2)/ 3.06126238d-04, 1.35807036d-04/
      data fpppp(10, 3),fpppp(10, 4)/-3.40385810d-04, 2.53979098d-04/
      data fpppp(10, 5),fpppp(10, 6)/-1.31935933d-04, 1.73413765d-05/
      data fpppp(10, 7),fpppp(10, 8)/ 7.45429244d-06, 3.85709874d-05/
      data fpppp(10, 9),fpppp(10,10)/-5.89316900d-05, 2.63934494d-05/
      data fpppp(10,11),fpppp(10,12)/-7.96950935d-06, 6.68711854d-08/
      data fpppp(10,13),fpppp(10,14)/-3.43823891d-05,-1.56183709d-05/
      data fpppp(10,15),fpppp(10,16)/-9.28098919d-05, 1.51972257d-04/
      data fpppp(10,17),fpppp(10,18)/-3.39357600d-06, 2.94268193d-04/
      data fpppp(10,19) /             4.62493504d-04 /
      data fpppp(11, 1),fpppp(11, 2)/-1.35900487d-04,-6.64481544d-05/
      data fpppp(11, 3),fpppp(11, 4)/ 1.14182288d-04,-9.43106525d-05/
      data fpppp(11, 5),fpppp(11, 6)/ 2.26522074d-05,-1.78606115d-05/
      data fpppp(11, 7),fpppp(11, 8)/ 5.72801559d-05,-9.28619145d-05/
      data fpppp(11, 9),fpppp(11,10)/ 9.83113071d-05,-5.42695683d-05/
      data fpppp(11,11),fpppp(11,12)/ 1.11545984d-06,-2.96523094d-06/
      data fpppp(11,13),fpppp(11,14)/-8.35469761d-06,-4.70643193d-07/
      data fpppp(11,15),fpppp(11,16)/-9.59384466d-08, 3.96390189d-06/
      data fpppp(11,17),fpppp(11,18)/ 4.01173545d-05, 7.89633852d-05/
      data fpppp(11,19) /             1.17040108d-04 /
      data fpppp(12, 1),fpppp(12, 2)/ 3.06961859d-05, 1.40806330d-05/
      data fpppp(12, 3),fpppp(12, 4)/-2.41040230d-05, 1.32031818d-05/
      data fpppp(12, 5),fpppp(12, 6)/-1.15668938d-05, 8.53739009d-06/
      data fpppp(12, 7),fpppp(12, 8)/-2.73462005d-05, 4.55254873d-05/
      data fpppp(12, 9),fpppp(12,10)/-3.89695214d-05, 2.51559381d-05/
      data fpppp(12,11),fpppp(12,12)/-5.70480468d-06,-1.11871631d-05/
      data fpppp(12,13),fpppp(12,14)/-2.85348322d-06,-5.65118966d-06/
      data fpppp(12,15),fpppp(12,16)/-2.74315835d-06, 3.14846513d-08/
      data fpppp(12,17),fpppp(12,18)/ 6.90715640d-05,-5.83470682d-06/
      data fpppp(12,19) /            -3.96134490d-05 /
      data fpppp(13, 1),fpppp(13, 2)/-8.57902378d-06,-4.61640200d-06/
      data fpppp(13, 3),fpppp(13, 4)/ 4.87995575d-06,-4.99176241d-06/
      data fpppp(13, 5),fpppp(13, 6)/ 2.14572050d-06,-2.38889177d-06/
      data fpppp(13, 7),fpppp(13, 8)/ 8.23948984d-06,-1.04543431d-05/
      data fpppp(13, 9),fpppp(13,10)/ 7.87929794d-06,-5.81774100d-06/
      data fpppp(13,11),fpppp(13,12)/ 1.80513999d-06, 2.21499204d-06/
      data fpppp(13,13),fpppp(13,14)/-1.71820685d-06,-1.16597551d-06/
      data fpppp(13,15),fpppp(13,16)/ 1.07691392d-06, 3.96858262d-06/
      data fpppp(13,17),fpppp(13,18)/-5.82478899d-06, 4.75911922d-06/
      data fpppp(13,19) /             9.52094749d-06 /
      data fpppp(14, 1),fpppp(14, 2)/ 2.96172040d-06, 1.58140551d-06/
      data fpppp(14, 3),fpppp(14, 4)/-1.59266176d-06, 1.51140431d-06/
      data fpppp(14, 5),fpppp(14, 6)/-7.94740453d-07, 6.31375780d-07/
      data fpppp(14, 7),fpppp(14, 8)/-2.17092550d-06, 2.76511503d-06/
      data fpppp(14, 9),fpppp(14,10)/-2.05289429d-06, 1.58946923d-06/
      data fpppp(14,11),fpppp(14,12)/-7.09117684d-07,-5.73209680d-07/
      data fpppp(14,13),fpppp(14,14)/ 4.29722574d-07, 2.92894816d-07/
      data fpppp(14,15),fpppp(14,16)/-3.49016788d-07,-8.70446795d-07/
      data fpppp(14,17),fpppp(14,18)/ 1.95226559d-06,-9.05770153d-07/
      data fpppp(14,19) /            -2.14073500d-06 /
      data fpppp(15, 1),fpppp(15, 2)/-5.24744215d-07,-2.81947257d-07/
      data fpppp(15, 3),fpppp(15, 4)/ 3.01486528d-07,-3.07308586d-07/
      data fpppp(15, 5),fpppp(15, 6)/ 1.27261085d-07,-1.48236693d-07/
      data fpppp(15, 7),fpppp(15, 8)/ 5.14693749d-07,-6.60418024d-07/
      data fpppp(15, 9),fpppp(15,10)/ 5.06001625d-07,-3.70724603d-07/
      data fpppp(15,11),fpppp(15,12)/ 9.89630178d-08, 1.36906256d-07/
      data fpppp(15,13),fpppp(15,14)/-1.05554031d-07,-6.81809906d-08/
      data fpppp(15,15),fpppp(15,16)/ 7.13327712d-08, 2.52063651d-07/
      data fpppp(15,17),fpppp(15,18)/-3.53889269d-07, 3.06565934d-07/
      data fpppp(15,19) /             6.05190245d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-2.29732021d-07,-1.13314543d-07/
      data fpppp(16, 3),fpppp(16, 4)/ 6.66921218d-08, 2.41366351d-08/
      data fpppp(16, 5),fpppp(16, 6)/ 9.02404168d-08, 7.93593105d-08/
      data fpppp(16, 7),fpppp(16, 8)/-2.71595975d-07, 3.20764451d-07/
      data fpppp(16, 9),fpppp(16,10)/-2.47966324d-07, 1.32883197d-07/
      data fpppp(16,11),fpppp(16,12)/ 8.96361378d-08,-8.40660381d-08/
      data fpppp(16,13),fpppp(16,14)/ 7.80038668d-08, 5.93031429d-08/
      data fpppp(16,15),fpppp(16,16)/ 2.30776011d-08,-2.55913277d-07/
      data fpppp(16,17),fpppp(16,18)/-1.06059261d-07,-4.03321650d-07/
      data fpppp(16,19) /            -5.72836495d-07 /
 

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
      subroutine c_2(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  4.36860000d-03 /
      data f( 1, 3),f( 1, 4) /  4.70281000d-02 ,  3.00994300d-01 /
      data f( 1, 5),f( 1, 6) /  4.79334100d-01 ,  2.14118700d-01 /
      data f( 1, 7),f( 1, 8) /  5.80496000d-02 ,  2.42705000d-02 /
      data f( 1, 9),f( 1,10) /  1.86612000d-02 ,  4.08940000d-03 /
      data f( 1,11),f( 1,12) / -1.23433000d-02 , -1.96644000d-02 /
      data f( 1,13),f( 1,14) / -1.68034000d-02 ,  6.63000000d-05 /
      data f( 1,15),f( 1,16) /  1.46350900d-01 ,  4.69804200d-01 /
      data f( 1,17),f( 1,18) /  2.59907800d-01 ,  3.21972000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  7.75940000d-03 /
      data f( 2, 3),f( 2, 4) /  2.35161000d-02 ,  1.43606200d-01 /
      data f( 2, 5),f( 2, 6) /  2.22893700d-01 ,  1.08799700d-01 /
      data f( 2, 7),f( 2, 8) /  3.24154000d-02 ,  1.56954000d-02 /
      data f( 2, 9),f( 2,10) /  1.12786000d-02 , -3.95100000d-04 /
      data f( 2,11),f( 2,12) / -1.27908000d-02 , -1.81056000d-02 /
      data f( 2,13),f( 2,14) / -1.60104000d-02 ,  3.76380000d-03 /
      data f( 2,15),f( 2,16) /  7.12416000d-02 ,  2.13100100d-01 /
      data f( 2,17),f( 2,18) /  1.19798100d-01 ,  1.67821000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  6.24130000d-03 /
      data f( 3, 3),f( 3, 4) /  1.27994000d-02 ,  6.48236000d-02 /
      data f( 3, 5),f( 3, 6) /  9.80092000d-02 ,  4.99661000d-02 /
      data f( 3, 7),f( 3, 8) /  1.64757000d-02 ,  1.00771000d-02 /
      data f( 3, 9),f( 3,10) /  7.96060000d-03 , -7.06400000d-04 /
      data f( 3,11),f( 3,12) / -1.05688000d-02 , -1.55087000d-02 /
      data f( 3,13),f( 3,14) / -1.53811000d-02 , -8.27670000d-03 /
      data f( 3,15),f( 3,16) /  2.45452000d-02 ,  9.37435000d-02 /
      data f( 3,17),f( 3,18) /  4.78577000d-02 , -2.12900000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  1.05890000d-03 /
      data f( 4, 3),f( 4, 4) /  4.59020000d-03 ,  2.47497000d-02 /
      data f( 4, 5),f( 4, 6) /  4.19257000d-02 ,  2.23590000d-02 /
      data f( 4, 7),f( 4, 8) /  7.57820000d-03 ,  7.24890000d-03 /
      data f( 4, 9),f( 4,10) /  6.20120000d-03 , -5.46800000d-04 /
      data f( 4,11),f( 4,12) / -8.28650000d-03 , -1.27206000d-02 /
      data f( 4,13),f( 4,14) / -1.40552000d-02 , -9.63960000d-03 /
      data f( 4,15),f( 4,16) /  5.41980000d-03 ,  3.80742000d-02 /
      data f( 4,17),f( 4,18) /  1.46563000d-02 , -9.79450000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -8.33740000d-03 /
      data f( 5, 3),f( 5, 4) / -3.22040000d-03 ,  7.62030000d-03 /
      data f( 5, 5),f( 5, 6) /  1.65701000d-02 ,  9.53300000d-03 /
      data f( 5, 7),f( 5, 8) /  5.98650000d-03 ,  5.77160000d-03 /
      data f( 5, 9),f( 5,10) /  4.79990000d-03 , -4.29900000d-04 /
      data f( 5,11),f( 5,12) / -6.33550000d-03 , -9.59910000d-03 /
      data f( 5,13),f( 5,14) / -1.18926000d-02 , -1.03273000d-02 /
      data f( 5,15),f( 5,16) / -3.00030000d-03 ,  1.26275000d-02 /
      data f( 5,17),f( 5,18) /  1.17450000d-03 , -1.04346000d-02 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -2.80190000d-03 /
      data f( 6, 3),f( 6, 4) / -2.70720000d-03 ,  1.15970000d-03 /
      data f( 6, 5),f( 6, 6) /  5.92410000d-03 ,  4.73670000d-03 /
      data f( 6, 7),f( 6, 8) /  4.64280000d-03 ,  4.61120000d-03 /
      data f( 6, 9),f( 6,10) /  3.65600000d-03 , -3.69700000d-04 /
      data f( 6,11),f( 6,12) / -4.71430000d-03 , -7.22830000d-03 /
      data f( 6,13),f( 6,14) / -9.44820000d-03 , -9.28340000d-03 /
      data f( 6,15),f( 6,16) / -6.25430000d-03 ,  1.52570000d-03 /
      data f( 6,17),f( 6,18) / -1.14380000d-03 , -7.02760000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -1.66840000d-03 /
      data f( 7, 3),f( 7, 4) / -1.48890000d-03 , -3.49200000d-04 /
      data f( 7, 5),f( 7, 6) /  2.58740000d-03 ,  3.28690000d-03 /
      data f( 7, 7),f( 7, 8) /  3.78300000d-03 ,  3.62370000d-03 /
      data f( 7, 9),f( 7,10) /  2.73250000d-03 , -3.48800000d-04 /
      data f( 7,11),f( 7,12) / -3.33760000d-03 , -5.31920000d-03 /
      data f( 7,13),f( 7,14) / -7.18380000d-03 , -7.68220000d-03 /
      data f( 7,15),f( 7,16) / -6.80000000d-03 , -2.53840000d-03 /
      data f( 7,17),f( 7,18) /  1.20020000d-03 ,  1.67350000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -8.77700000d-04 /
      data f( 8, 3),f( 8, 4) / -5.91100000d-04 ,  5.83500000d-04 /
      data f( 8, 5),f( 8, 6) /  2.05310000d-03 ,  2.93790000d-03 /
      data f( 8, 7),f( 8, 8) /  3.30710000d-03 ,  2.78380000d-03 /
      data f( 8, 9),f( 8,10) /  1.99310000d-03 , -3.46900000d-04 /
      data f( 8,11),f( 8,12) / -2.34140000d-03 , -3.85310000d-03 /
      data f( 8,13),f( 8,14) / -5.15630000d-03 , -5.79050000d-03 /
      data f( 8,15),f( 8,16) / -5.61960000d-03 , -2.86730000d-03 /
      data f( 8,17),f( 8,18) /  5.28500000d-04 ,  1.26760000d-03 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -1.69000000d-04 /
      data f( 9, 3),f( 9, 4) /  2.82100000d-04 ,  1.41840000d-03 /
      data f( 9, 5),f( 9, 6) /  2.07700000d-03 ,  2.39830000d-03 /
      data f( 9, 7),f( 9, 8) /  2.14560000d-03 ,  1.50840000d-03 /
      data f( 9, 9),f( 9,10) /  9.46300000d-04 , -3.12400000d-04 /
      data f( 9,11),f( 9,12) / -1.12220000d-03 , -1.96610000d-03 /
      data f( 9,13),f( 9,14) / -2.68850000d-03 , -2.96090000d-03 /
      data f( 9,15),f( 9,16) / -2.55990000d-03 , -1.18390000d-03 /
      data f( 9,17),f( 9,18) /  5.08000000d-04 ,  5.84900000d-04 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -5.20000000d-06 /
      data f(10, 3),f(10, 4) /  4.16500000d-04 ,  1.14670000d-03 /
      data f(10, 5),f(10, 6) /  1.46480000d-03 ,  1.46400000d-03 /
      data f(10, 7),f(10, 8) /  1.17100000d-03 ,  7.93300000d-04 /
      data f(10, 9),f(10,10) /  3.13600000d-04 , -8.58000000d-05 /
      data f(10,11),f(10,12) / -5.17200000d-04 , -9.85500000d-04 /
      data f(10,13),f(10,14) / -1.38530000d-03 , -1.49730000d-03 /
      data f(10,15),f(10,16) / -1.14650000d-03 , -2.83600000d-04 /
      data f(10,17),f(10,18) /  4.98600000d-04 ,  3.45600000d-04 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -4.60000000d-06 /
      data f(11, 3),f(11, 4) /  2.96000000d-04 ,  5.31100000d-04 /
      data f(11, 5),f(11, 6) /  6.21000000d-04 ,  5.95300000d-04 /
      data f(11, 7),f(11, 8) /  4.56700000d-04 ,  3.01700000d-04 /
      data f(11, 9),f(11,10) /  1.22000000d-04 ,  1.98000000d-05 /
      data f(11,11),f(11,12) / -1.72200000d-04 , -4.08500000d-04 /
      data f(11,13),f(11,14) / -6.17900000d-04 , -6.87500000d-04 /
      data f(11,15),f(11,16) / -5.22400000d-04 , -1.44500000d-04 /
      data f(11,17),f(11,18) /  2.15500000d-04 ,  9.60000000d-05 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  2.92000000d-05 /
      data f(12, 3),f(12, 4) /  9.62000000d-05 ,  1.52300000d-04 /
      data f(12, 5),f(12, 6) /  1.86100000d-04 ,  1.89200000d-04 /
      data f(12, 7),f(12, 8) /  1.58300000d-04 ,  8.19000000d-05 /
      data f(12, 9),f(12,10) /  6.09000000d-05 ,  8.50000000d-06 /
      data f(12,11),f(12,12) / -6.09000000d-05 , -1.55600000d-04 /
      data f(12,13),f(12,14) / -2.78100000d-04 , -3.38300000d-04 /
      data f(12,15),f(12,16) / -2.93700000d-04 , -1.69300000d-04 /
      data f(12,17),f(12,18) / -4.75000000d-05 , -1.20000000d-06 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  4.18000000d-05 /
      data f(13, 3),f(13, 4) /  4.05000000d-05 , -2.60000000d-06 /
      data f(13, 5),f(13, 6) / -1.49000000d-05 ,  9.90000000d-06 /
      data f(13, 7),f(13, 8) /  3.86000000d-05 ,  5.15000000d-05 /
      data f(13, 9),f(13,10) /  4.20000000d-05 ,  1.53000000d-05 /
      data f(13,11),f(13,12) / -1.69000000d-05 , -3.62000000d-05 /
      data f(13,13),f(13,14) / -5.87000000d-05 , -8.18000000d-05 /
      data f(13,15),f(13,16) / -9.43000000d-05 , -8.10000000d-05 /
      data f(13,17),f(13,18) / -4.70000000d-05 , -1.32000000d-05 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -3.30000000d-06 /
      data f(14, 3),f(14, 4) / -7.50000000d-06 , -4.10000000d-06 /
      data f(14, 5),f(14, 6) /  7.50000000d-06 ,  1.98000000d-05 /
      data f(14, 7),f(14, 8) /  2.63000000d-05 ,  2.46000000d-05 /
      data f(14, 9),f(14,10) /  2.06000000d-05 ,  1.62000000d-05 /
      data f(14,11),f(14,12) /  1.21000000d-05 ,  9.40000000d-06 /
      data f(14,13),f(14,14) /  5.80000000d-06 ,  1.30000000d-06 /
      data f(14,15),f(14,16) / -3.80000000d-06 , -6.30000000d-06 /
      data f(14,17),f(14,18) / -4.30000000d-06 , -1.90000000d-06 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -6.00000000d-07 /
      data f(15, 3),f(15, 4) / -7.00000000d-07 ,  6.00000000d-07 /
      data f(15, 5),f(15, 6) /  3.40000000d-06 ,  6.20000000d-06 /
      data f(15, 7),f(15, 8) /  8.00000000d-06 ,  7.90000000d-06 /
      data f(15, 9),f(15,10) /  8.00000000d-06 ,  6.20000000d-06 /
      data f(15,11),f(15,12) /  3.50000000d-06 ,  1.70000000d-06 /
      data f(15,13),f(15,14) / -2.00000000d-07 , -1.50000000d-06 /
      data f(15,15),f(15,16) / -2.50000000d-06 , -2.50000000d-06 /
      data f(15,17),f(15,18) / -1.30000000d-06 , -6.00000000d-07 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 , -5.00000000d-07 /
      data f(16, 3),f(16, 4) / -1.10000000d-06 , -2.20000000d-06 /
      data f(16, 5),f(16, 6) / -2.60000000d-06 , -2.40000000d-06 /
      data f(16, 7),f(16, 8) / -1.40000000d-06 ,  1.00000000d-07 /
      data f(16, 9),f(16,10) /  1.90000000d-06 ,  2.60000000d-06 /
      data f(16,11),f(16,12) /  1.90000000d-06 ,  6.00000000d-07 /
      data f(16,13),f(16,14) / -1.70000000d-06 , -2.40000000d-06 /
      data f(16,15),f(16,16) / -2.50000000d-06 , -1.70000000d-06 /
      data f(16,17),f(16,18) / -1.20000000d-06 ,  2.00000000d-07 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-1.19703410d-03/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-1.43159978d-01, 8.24791989d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 5.72670301d-01, 3.16457130d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 2.93809766d+00,-6.23624148d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 4.83772301d+00,-7.45270564d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.54741777d+00, 3.25987299d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 3.03068551d-01, 9.61991670d-04/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 7.94548503d-02, 2.29560326d-04/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 1.63212717d-01,-1.90044972d-04/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.93719577d-01,-7.13043631d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.29766664d-01, 1.06912718d-04/
      data fpp( 1,12,1),fpp( 1,12,2)/ 4.56517092d-02, 1.26175566d-04/
      data fpp( 1,13,1),fpp( 1,13,2)/-2.50435021d-02,-6.88982123d-07/
      data fpp( 1,14,1),fpp( 1,14,2)/-1.01284320d+00, 7.17102362d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 7.50923021d-01, 4.89717353d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.23786225d+00,-9.67567449d-03/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.43186793d+00, 1.80454244d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/-4.36702159d-01, 1.38865273d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 4.37165063d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-7.97045950d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-1.44077543d-01,-9.27111000d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.28696899d-01, 1.64772835d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.01405717d+00,-2.38198299d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.32899148d+00,-3.14309115d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.15451446d+00, 1.20767291d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.53565398d-01, 5.74981521d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 7.11852993d-02, 7.22590073d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.03719566d-01,-1.25825550d-04/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.10665846d-01,-4.37080575d-06/
      data fpp( 2,11,1),fpp( 2,11,2)/ 7.11391720d-02, 9.99887734d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/ 2.88990815d-02, 2.92697123d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/-5.20549571d-03, 2.27532378d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/-4.75443603d-01, 1.21340777d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 6.71216459d-01, 2.14932051d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.50834801d+00,-4.25578083d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.72047664d+00, 7.64172800d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.25875682d-01, 6.16249629d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 1.94486269d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-4.08032767d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-1.68648481d-02,-7.74164669d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 3.18371048d-02, 7.36706634d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 7.96498652d-01,-1.41444070d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.57969608d+00,-1.30124635d-03/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 8.07334403d-01, 4.72707483d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.36844859d-01, 2.83578421d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 7.93239523d-02, 1.84868334d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 3.15990185d-02,-1.00599755d-04/
      data fpp( 3,10,1),fpp( 3,10,2)/-1.04029625d-02,-9.11781508d-06/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.38983520d-02, 6.53470149d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-5.53303541d-03, 4.30797554d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 2.13104850d-02, 6.63839634d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 5.53917612d-01, 1.09992391d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 8.26146144d-01, 1.03669647d-03/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.33087072d+00,-2.07419428d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 9.11620527d-01, 3.55034653d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/ 4.15804887d-01, 4.08001671d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 1.13990066d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-1.06845377d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-3.38108064d-01, 4.74675347d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-7.99203177d-02, 2.36202363d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 6.06253221d-01, 4.81357951d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 6.72374209d-01,-6.07755543d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 3.00122930d-01, 1.78324378d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 2.55385168d-01, 1.81612030d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 3.00338913d-02,-3.76824983d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 3.67435971d-03,-7.39860367d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/ 1.58100353d-03,-8.39135495d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/-6.50076416d-03, 4.80494565d-05/
      data fpp( 4,12,1),fpp( 4,12,2)/ 2.19130601d-02, 1.45295290d-05/
      data fpp( 4,13,1),fpp( 4,13,2)/ 2.44535558d-02, 7.98024275d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/-1.38586843d-01, 1.12727611d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 1.59848964d-01, 5.13734528d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 7.21264122d-01,-1.01051087d-03/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.43891257d-01, 1.63970966d-04/
      data fpp( 4,18,1),fpp( 4,18,2)/ 1.49496134d-01, 2.92653010d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 7.20134995d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 2.09934917d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 7.37212105d-01, 1.38376166d-04/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 3.47634166d-01, 4.38244199d-05/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.20163466d-01, 2.97481547d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.39992084d-01,-2.76271039d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.09338876d-01, 1.16121999d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-6.25155318d-02, 2.12190409d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 3.17548250d-03,-1.10216292d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 7.41854263d-03,-6.22183892d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/-2.32605163d-03,-5.51028035d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/-9.79359139d-03, 4.37115106d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-3.21092049d-02,-1.08157620d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/ 6.38029192d-03, 5.77575374d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 1.01709762d-01, 1.13136125d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 1.40253001d-01, 2.42690013d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/ 3.17462795d-01,-4.84025663d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 2.70754445d-01, 6.85646401d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 4.00205778d-02, 2.00401103d-04/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 4.52452949d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 2.08437031d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-3.70970355d-01, 2.76985938d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-6.20463469d-02, 4.21579216d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.13412916d-01, 3.00017198d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.74097453d-01,-1.08314801d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 6.69765664d-02, 4.61494834d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 3.18769590d-02,-1.06731326d-05/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 4.79917873d-03, 2.81047150d-07/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 5.26146978d-03,-4.58670560d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/-7.81797014d-04,-1.04282325d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-3.79487027d-03, 3.09043490d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/-6.08124036d-03,-1.27385726d-05/
      data fpp( 6,13,1),fpp( 6,13,2)/-7.70472347d-03, 3.76959416d-05/
      data fpp( 6,14,1),fpp( 6,14,2)/-8.51220504d-03, 5.03680617d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 5.40540328d-02, 1.14014834d-04/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.60619699d-01,-1.76042141d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.47615964d-01,-3.68162697d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 2.97486555d-01, 1.30449220d-04/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 2.89703390d-04/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 2.58976159d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 8.63693154d-02, 2.13957683d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 6.31622127d-03,-6.06689033d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 6.89398685d-02, 3.86429878d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 6.00131042d-02,-4.61512623d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 2.47298585d-02, 1.17360614d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 7.59269575d-03,-1.29969833d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 3.56280260d-03, 9.27871941d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 4.59557824d-03,-3.46285044d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/-4.41760314d-04, 6.18014577d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-1.17019275d-02, 1.54579214d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.28208336d-02,-7.57983119d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-2.56139805d-03, 2.18814034d-05/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.59340580d-02, 2.02621759d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/ 4.97758681d-02, 5.28497262d-05/
      data fpp( 7,16,1),fpp( 7,16,2)/ 9.57134077d-02,-1.06611225d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/-1.61873301d-01,-4.15852363d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/-4.35851798d-01,-1.89159325d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-1.15590338d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 1.44698376d-05/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-2.59269066d-02, 1.15153249d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-1.12935382d-02, 9.32686298d-06/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-2.29323905d-02, 4.45722320d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 6.21013043d-03,-9.45575580d-06/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-7.76000524d-04,-1.72220002d-06/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-4.66274202d-03,-1.45914441d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 3.08961087d-03, 6.53797651d-06/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 3.97121726d-03,-2.76044619d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/-3.01161731d-04, 1.09218711d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/-6.47241967d-03, 4.64697736d-06/
      data fpp( 8,12,1),fpp( 8,12,2)/-9.08542517d-03,-5.41780562d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.75846843d-02, 1.00301449d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.16490272d-02, 5.61200999d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/ 5.75749487d-03, 3.60310511d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.68066698d-02, 1.01985945d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/ 4.75222381d-02,-3.82154293d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/ 7.98706359d-02,-1.67388773d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.52310613d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 6.21675872d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.86981201d-03, 4.86748256d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/-3.86749601d-03, 1.15193110d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-4.31651284d-03,-9.83272672d-06/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-7.66819338d-03,-8.50404150d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-4.09692769d-03,-7.00365668d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.32812817d-03,-5.57496914d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 4.11476608d-03, 6.23353325d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.98855910d-03,-1.48531639d-05/
      data fpp( 9,10,1),fpp( 9,10,2)/ 2.27561535d-03, 1.13831222d-05/
      data fpp( 9,11,1),fpp( 9,11,2)/-3.72677724d-03,-3.74532508d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-5.52830770d-03, 1.55217809d-06/
      data fpp( 9,13,1),fpp( 9,13,2)/-5.48524796d-03, 4.82661273d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-8.78744757d-03, 6.14137098d-06/
      data fpp( 9,15,1),fpp( 9,15,2)/-1.59516687d-02, 1.10119034d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/-1.04817132d-02, 8.31101560d-06/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.20213142d-02,-2.53019657d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/-1.68447589d-02,-4.00315265d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 1.60657632d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 5.79707629d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.98609147d-03, 3.58084742d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/-9.41477762d-04, 5.49353402d-06/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.29905815d-03,-7.04498351d-06/
      data fpp(10, 5,1),fpp(10, 5,2)/ 6.08893104d-04,-2.03959998d-06/
      data fpp(10, 6,1),fpp(10, 6,2)/ 2.36246129d-03,-3.93061658d-06/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.35897933d-03, 2.30066281d-07/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.46257479d-03,-2.07164855d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/ 3.60329635d-03, 1.93652791d-06/
      data fpp(10,10,1),fpp(10,10,2)/-1.59754967d-03,-8.56463078d-07/
      data fpp(10,11,1),fpp(10,11,2)/-1.65297138d-03,-4.30675594d-07/
      data fpp(10,12,1),fpp(10,12,2)/-2.79134405d-03, 3.65165454d-07/
      data fpp(10,13,1),fpp(10,13,2)/-4.14682382d-03, 3.08001378d-06/
      data fpp(10,14,1),fpp(10,14,2)/-4.42618258d-03, 4.58277943d-06/
      data fpp(10,15,1),fpp(10,15,2)/-3.68707026d-03, 6.35686849d-06/
      data fpp(10,16,1),fpp(10,16,2)/-4.24606708d-03, 7.15746589d-07/
      data fpp(10,17,1),fpp(10,17,2)/ 9.79268507d-04,-1.40618549d-05/
      data fpp(10,18,1),fpp(10,18,2)/ 4.13589983d-03,-5.80327185d-07/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 4.82716359d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 6.65240067d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/ 7.54479671d-04, 3.26519866d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/-4.40683251d-04,-1.40119532d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 1.50641963d-03,-1.59041738d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 2.05733953d-03,-9.49135148d-07/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.95288150d-03,-1.54904203d-06/
      data fpp(11, 7,1),fpp(11, 7,2)/ 1.73997189d-03, 3.71303248d-07/
      data fpp(11, 8,1),fpp(11, 8,2)/ 1.09751787d-03,-9.20170968d-07/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.80114136d-04, 1.82738062d-06/
      data fpp(11,10,1),fpp(11,10,2)/-3.32913469d-04,-1.73935153d-06/
      data fpp(11,11,1),fpp(11,11,2)/-9.37881232d-04,-2.57974500d-07/
      data fpp(11,12,1),fpp(11,12,2)/-1.09851525d-03, 1.13249531d-07/
      data fpp(11,13,1),fpp(11,13,2)/-1.36163589d-03, 1.41897638d-06/
      data fpp(11,14,1),fpp(11,14,2)/-1.50858465d-03, 2.59884497d-06/
      data fpp(11,15,1),fpp(11,15,2)/-1.38881216d-03, 2.26764375d-06/
      data fpp(11,16,1),fpp(11,16,2)/ 6.12020763d-07, 1.09858003d-06/
      data fpp(11,17,1),fpp(11,17,2)/-4.20715296d-04,-7.73596385d-06/
      data fpp(11,18,1),fpp(11,18,2)/-2.24832221d-04, 1.07527539d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 4.84486231d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 8.49290625d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/-2.35027218d-04, 4.09418749d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/ 8.01010765d-04,-2.18965623d-07/
      data fpp(12, 4,1),fpp(12, 4,2)/ 9.56579625d-04,-1.87556258d-07/
      data fpp(12, 5,1),fpp(12, 5,2)/ 9.75348760d-04,-3.68809346d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 9.28412699d-04,-1.79206356d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 6.62733133d-04,-9.54365228d-07/
      data fpp(12, 8,1),fpp(12, 8,2)/ 6.70553710d-04, 1.26666727d-06/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.49160194d-04,-7.88303843d-07/
      data fpp(12,10,1),fpp(12,10,2)/ 1.23603545d-04, 2.54810293d-09/
      data fpp(12,11,1),fpp(12,11,2)/-2.04303687d-04,-2.41888569d-07/
      data fpp(12,12,1),fpp(12,12,2)/-5.92994945d-04,-5.52993827d-07/
      data fpp(12,13,1),fpp(12,13,2)/-6.69032611d-04, 7.85863875d-07/
      data fpp(12,14,1),fpp(12,14,2)/-5.93878805d-04, 1.14753833d-06/
      data fpp(12,15,1),fpp(12,15,2)/-2.47281117d-04, 9.11982825d-07/
      data fpp(12,16,1),fpp(12,16,2)/ 3.10018995d-04,-7.46962306d-09/
      data fpp(12,17,1),fpp(12,17,2)/ 1.18599268d-03,-1.03810433d-06/
      data fpp(12,18,1),fpp(12,18,2)/ 4.21029057d-04,-3.70113048d-07/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00,-1.87443476d-07/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-4.79657442d-07/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.15818222d-06,-3.59685117d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/-1.19290671d-04,-6.67602091d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-6.74869057d-06, 5.22093482d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/ 5.80839533d-05, 4.27228165d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/ 3.57211513d-05,-5.00614008d-09/
      data fpp(13, 7,1),fpp(13, 7,2)/ 4.41465985d-06,-1.73203604d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.05220066d-04,-2.50179443d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/-3.76235130d-05,-1.70078625d-07/
      data fpp(13,10,1),fpp(13,10,2)/-2.79539013d-05,-1.01506059d-07/
      data fpp(13,11,1),fpp(13,11,2)/ 1.02516783d-05, 2.46102861d-07/
      data fpp(13,12,1),fpp(13,12,2)/ 9.84245927d-06,-1.08905385d-07/
      data fpp(13,13,1),fpp(13,13,2)/-7.32842205d-05,-2.48132193d-09/
      data fpp(13,14,1),fpp(13,14,2)/-1.15471259d-04, 8.28306725d-08/
      data fpp(13,15,1),fpp(13,15,2)/-1.11750571d-04, 3.07158632d-07/
      data fpp(13,16,1),fpp(13,16,2)/-1.02962996d-04, 2.36534800d-07/
      data fpp(13,17,1),fpp(13,17,2)/-1.88620388d-04,-1.12978307d-08/
      data fpp(13,18,1),fpp(13,18,2)/-5.62710602d-05,-2.03343477d-07/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00,-4.11328262d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00,-9.16418139d-08/
      data fpp(14, 2,1),fpp(14, 2,2)/ 1.85381555d-05,-1.37163721d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/ 5.24666290d-05, 9.25073024d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/ 4.40625933d-06, 9.96871623d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/-2.53262398d-05, 7.44048185d-10/
      data fpp(14, 6,1),fpp(14, 6,2)/-1.86198036d-05,-6.06633551d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/-3.96054581d-06,-1.06090628d-07/
      data fpp(14, 8,1),fpp(14, 8,2)/ 3.12333435d-05,-6.97413350d-09/
      data fpp(14, 9,1),fpp(14, 9,2)/ 1.28904421d-05,-4.01283815d-09/
      data fpp(14,10,1),fpp(14,10,2)/ 3.00993138d-06,-9.74513915d-10/
      data fpp(14,11,1),fpp(14,11,2)/-1.71031913d-05, 2.59108938d-08/
      data fpp(14,12,1),fpp(14,12,2)/-2.28299055d-05,-1.86690613d-08/
      data fpp(14,13,1),fpp(14,13,2)/-7.08103294d-06,-5.23464859d-09/
      data fpp(14,14,1),fpp(14,14,2)/-1.49681949d-06,-1.43923443d-08/
      data fpp(14,15,1),fpp(14,15,2)/-3.55772771d-06, 2.68040260d-08/
      data fpp(14,16,1),fpp(14,16,2)/ 1.02948997d-06, 6.31762405d-08/
      data fpp(14,17,1),fpp(14,17,2)/ 3.54148240d-05,-9.50898782d-09/
      data fpp(14,18,1),fpp(14,18,2)/ 1.12486521d-05,-1.14028920d-09/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00,-1.59298554d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00,-3.86940079d-09/
      data fpp(15, 2,1),fpp(15, 2,2)/-2.94439960d-07, 4.73880157d-09/
      data fpp(15, 3,1),fpp(15, 3,2)/-8.37584548d-06, 1.49141945d-08/
      data fpp(15, 4,1),fpp(15, 4,2)/-1.57634673d-06, 1.96044204d-08/
      data fpp(15, 5,1),fpp(15, 5,2)/ 3.47100601d-06,-3.33187627d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/ 3.50806303d-06,-6.27691538d-09/
      data fpp(15, 7,1),fpp(15, 7,2)/ 2.42752338d-06,-3.15604622d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-4.41330770d-06, 1.85187642d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/-7.38255274d-07,-3.05145947d-08/
      data fpp(15,10,1),fpp(15,10,2)/-4.35824191d-07,-1.04603856d-08/
      data fpp(15,11,1),fpp(15,11,2)/ 1.76108704d-06, 1.83561370d-08/
      data fpp(15,12,1),fpp(15,12,2)/ 1.52716263d-06,-8.96416233d-09/
      data fpp(15,13,1),fpp(15,13,2)/-4.14164771d-06, 1.15005124d-08/
      data fpp(15,14,1),fpp(15,14,2)/-7.39146264d-06,-1.03788710d-09/
      data fpp(15,15,1),fpp(15,15,2)/-7.81851785d-06, 1.06510360d-08/
      data fpp(15,16,1),fpp(15,16,2)/-7.50496403d-06, 1.84337429d-08/
      data fpp(15,17,1),fpp(15,17,2)/-1.25889082d-05,-1.23860077d-08/
      data fpp(15,18,1),fpp(15,18,2)/-3.72354840d-06, 1.11028792d-09/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00, 1.94485604d-09/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00, 2.17846437d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/-1.40106372d-05, 6.43071265d-10/
      data fpp(16, 3,1),fpp(16, 3,2)/-1.41249344d-05,-1.07507494d-08/
      data fpp(16, 4,1),fpp(16, 4,2)/-4.24968378d-06, 1.23599264d-08/
      data fpp(16, 5,1),fpp(16, 5,2)/ 5.41413985d-06, 3.31104367d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/ 8.58632563d-06, 1.03958989d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/ 6.58195259d-06, 3.10536075d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/ 5.38879671d-06, 7.18265811d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.40055621d-06,-1.38359932d-08/
      data fpp(16,10,1),fpp(16,10,2)/ 7.04612638d-06,-1.78386854d-08/
      data fpp(16,11,1),fpp(16,11,2)/ 1.30651708d-05, 1.19073482d-09/
      data fpp(16,12,1),fpp(16,12,2)/ 1.70960615d-05,-2.29242539d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 2.35261810d-05, 3.05062806d-08/
      data fpp(16,14,1),fpp(16,14,2)/ 2.78360885d-05,-3.10086868d-09/
      data fpp(16,15,1),fpp(16,15,2)/ 2.71335446d-05, 1.78971941d-08/
      data fpp(16,16,1),fpp(16,16,2)/ 2.10635534d-05,-1.44879077d-08/
      data fpp(16,17,1),fpp(16,17,2)/ 1.54198112d-05, 2.20544367d-08/
      data fpp(16,18,1),fpp(16,18,2)/ 4.14605992d-06,-1.97298391d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00,-3.91350805d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 1.17634104d-03, 7.60488909d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 1.99435180d-02, 1.15968639d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/-9.42790946d-02, 5.41236794d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 5.41738218d-04, 4.95349894d-03/
      data fpppp( 1, 9),fpppp( 1,10)/-1.91343994d-03,-4.94799588d-04/
      data fpppp( 1,11),fpppp( 1,12)/-1.77494807d-03, 6.38486934d-03/
      data fpppp( 1,13),fpppp( 1,14)/-2.29593447d-02, 3.04262403d-02/
      data fpppp( 1,15),fpppp( 1,16)/ 6.63483385d-02,-1.32429214d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 2.57925049d-02, 2.55046474d-02/
      data fpppp( 1,19) /             7.05052403d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 6.66736873d-04, 5.25740892d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.53147466d-02, 6.23875465d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/-6.24953231d-02, 3.43778579d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.39556916d-03, 3.15400306d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-1.11671949d-03,-2.22404285d-04/
      data fpppp( 2,11),fpppp( 2,12)/-7.82040653d-04, 3.18776194d-03/
      data fpppp( 2,13),fpppp( 2,14)/-1.14808763d-02, 1.65677314d-02/
      data fpppp( 2,15),fpppp( 2,16)/ 4.22238407d-02,-8.40348050d-02/
      data fpppp( 2,17),fpppp( 2,18)/ 1.64152038d-02, 1.48651334d-02/
      data fpppp( 2,19) /             4.24579427d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/-5.48953341d-03,-6.17850117d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/ 9.67068152d-03, 4.33663460d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-2.59050672d-02, 5.95008814d-03/
      data fpppp( 3, 7),fpppp( 3, 8)/ 8.21704242d-03,-2.04013952d-03/
      data fpppp( 3, 9),fpppp( 3,10)/ 5.31273985d-04, 2.58420743d-04/
      data fpppp( 3,11),fpppp( 3,12)/ 7.45438536d-04,-2.52853252d-03/
      data fpppp( 3,13),fpppp( 3,14)/ 1.04773838d-02,-9.03518627d-03/
      data fpppp( 3,15),fpppp( 3,16)/ 1.00406456d-02,-1.71776339d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 3.23140402d-03,-3.41909190d-04/
      data fpppp( 3,19) /             2.93687791d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 8.11752773d-03, 5.02001923d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 7.58014397d-03,-9.66144761d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/-6.13750649d-03, 7.90913749d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/-5.84823243d-03, 4.64698132d-03/
      data fpppp( 4, 9),fpppp( 4,10)/-8.00188133d-04, 9.74173174d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 4.01916515d-04, 5.72327726d-04/
      data fpppp( 4,13),fpppp( 4,14)/-4.24362713d-03, 6.46732710d-03/
      data fpppp( 4,15),fpppp( 4,16)/ 6.06289110d-03,-1.49401304d-02/
      data fpppp( 4,17),fpppp( 4,18)/ 3.37034925d-03, 4.37397943d-04/
      data fpppp( 4,19) /             3.57399835d-03 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.46416419d-02,-1.22983631d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/ 6.22769187d-03, 3.11402997d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/-3.84585259d-03,-2.75952923d-03/
      data fpppp( 5, 7),fpppp( 5, 8)/ 6.41189756d-03,-2.63533570d-03/
      data fpppp( 5, 9),fpppp( 5,10)/ 4.42567984d-04, 2.58045002d-05/
      data fpppp( 5,11),fpppp( 5,12)/-4.09162715d-04, 7.19961934d-04/
      data fpppp( 5,13),fpppp( 5,14)/ 1.17762160d-03,-2.02004995d-03/
      data fpppp( 5,15),fpppp( 5,16)/ 3.49540432d-03,-3.64157401d-03/
      data fpppp( 5,17),fpppp( 5,18)/-2.36419688d-03, 2.05683052d-03/
      data fpppp( 5,19) /             5.57967214d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.45683137d-02, 7.52738572d-03/
      data fpppp( 6, 3),fpppp( 6, 4)/-3.88419475d-03, 1.50858493d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/-3.00832320d-03, 1.96345886d-03/
      data fpppp( 6, 7),fpppp( 6, 8)/-5.24235477d-04, 6.14792676d-04/
      data fpppp( 6, 9),fpppp( 6,10)/-2.82530948d-04, 1.24997643d-04/
      data fpppp( 6,11),fpppp( 6,12)/-3.56480111d-05, 6.11965918d-05/
      data fpppp( 6,13),fpppp( 6,14)/-1.69365137d-04, 6.65224049d-04/
      data fpppp( 6,15),fpppp( 6,16)/ 1.31089210d-03,-3.26882673d-03/
      data fpppp( 6,17),fpppp( 6,18)/ 4.59025069d-03,-5.31971644d-03/
      data fpppp( 6,19) /            -1.01528137d-02 /
      data fpppp( 7, 1),fpppp( 7, 2)/-4.50984377d-03,-2.15496777d-03/
      data fpppp( 7, 3),fpppp( 7, 4)/ 3.14437028d-03,-1.86190887d-03/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.02405178d-05, 2.39557923d-04/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.20292763d-04, 6.57072015d-05/
      data fpppp( 7, 9),fpppp( 7,10)/-7.93614423d-05,-1.12468284d-04/
      data fpppp( 7,11),fpppp( 7,12)/ 1.55864858d-04, 9.74845170d-05/
      data fpppp( 7,13),fpppp( 7,14)/ 1.36897574d-04,-1.50913583d-04/
      data fpppp( 7,15),fpppp( 7,16)/ 1.38753799d-03,-4.67349462d-03/
      data fpppp( 7,17),fpppp( 7,18)/-9.05014380d-04, 7.31004480d-03/
      data fpppp( 7,19) /             1.42546529d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 9.98807737d-04, 5.55842686d-04/
      data fpppp( 8, 3),fpppp( 8, 4)/-7.88561983d-04, 1.02207201d-03/
      data fpppp( 8, 5),fpppp( 8, 6)/-8.52843662d-04, 2.21583529d-04/
      data fpppp( 8, 7),fpppp( 8, 8)/ 1.52472912d-04,-1.33129516d-04/
      data fpppp( 8, 9),fpppp( 8,10)/-3.21996379d-05,-4.73110551d-05/
      data fpppp( 8,11),fpppp( 8,12)/ 1.07511122d-04,-1.69238285d-04/
      data fpppp( 8,13),fpppp( 8,14)/ 2.16266797d-04, 1.70266080d-04/
      data fpppp( 8,15),fpppp( 8,16)/-2.09079225d-04, 2.84609995d-04/
      data fpppp( 8,17),fpppp( 8,18)/ 2.50622853d-04,-1.18913164d-03/
      data fpppp( 8,19) /            -2.22723832d-03 /
      data fpppp( 9, 1),fpppp( 9, 2)/-1.93908413d-04,-9.83049993d-05/
      data fpppp( 9, 3),fpppp( 9, 4)/ 1.30701209d-04,-1.07202364d-04/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.23948426d-04, 2.67854331d-05/
      data fpppp( 9, 7),fpppp( 9, 8)/-5.98627483d-05,-6.56395171d-05/
      data fpppp( 9, 9),fpppp( 9,10)/ 8.76501227d-05,-1.40165179d-04/
      data fpppp( 9,11),fpppp( 9,12)/ 9.56436648d-05, 9.64224787d-06/
      data fpppp( 9,13),fpppp( 9,14)/-2.35372448d-05,-1.16208829d-04/
      data fpppp( 9,15),fpppp( 9,16)/ 2.56651272d-04,-1.52345664d-04/
      data fpppp( 9,17),fpppp( 9,18)/-6.78420050d-05, 2.26683056d-04/
      data fpppp( 9,19) /             4.61202003d-04 /
      data fpppp(10, 1),fpppp(10, 2)/ 7.05782800d-05, 3.84225805d-05/
      data fpppp(10, 3),fpppp(10, 4)/-4.24262918d-05, 4.71509410d-05/
      data fpppp(10, 5),fpppp(10, 6)/-1.02455729d-05,-1.54316335d-05/
      data fpppp(10, 7),fpppp(10, 8)/-3.34509024d-05, 9.56598888d-05/
      data fpppp(10, 9),fpppp(10,10)/-1.66961088d-04, 1.31690407d-04/
      data fpppp(10,11),fpppp(10,12)/-5.10750838d-05, 7.63287015d-06/
      data fpppp(10,13),fpppp(10,14)/ 7.51717771d-06, 2.68656788d-05/
      data fpppp(10,15),fpppp(10,16)/-5.38716275d-05, 1.10734282d-04/
      data fpppp(10,17),fpppp(10,18)/-4.20055578d-05,-6.68343071d-05/
      data fpppp(10,19) /            -1.28209082d-04 /
      data fpppp(11, 1),fpppp(11, 2)/-6.58368196d-05,-2.86538065d-05/
      data fpppp(11, 3),fpppp(11, 4)/ 6.34734900d-05,-3.67042053d-05/
      data fpppp(11, 5),fpppp(11, 6)/-4.27647746d-07,-9.07879735d-07/
      data fpppp(11, 7),fpppp(11, 8)/-2.44792842d-06,-1.50730702d-05/
      data fpppp(11, 9),fpppp(11,10)/ 2.46295294d-05,-1.59550866d-05/
      data fpppp(11,11),fpppp(11,12)/ 1.20607113d-05,-5.62773398d-06/
      data fpppp(11,13),fpppp(11,14)/ 4.30102708d-06,-4.60606139d-06/
      data fpppp(11,15),fpppp(11,16)/ 3.01264939d-05,-3.97208133d-05/
      data fpppp(11,17),fpppp(11,18)/ 2.01116696d-05,-3.69324161d-06/
      data fpppp(11,19) /            -3.60175451d-06 /
      data fpppp(12, 1),fpppp(12, 2)/ 3.29665724d-05, 1.52294978d-05/
      data fpppp(12, 3),fpppp(12, 4)/-1.76206513d-05, 2.42496017d-06/
      data fpppp(12, 5),fpppp(12, 6)/-2.87172794d-07,-5.21858074d-06/
      data fpppp(12, 7),fpppp(12, 8)/ 8.03688539d-06,-1.05189522d-05/
      data fpppp(12, 9),fpppp(12,10)/ 8.28607777d-06,-4.87514683d-06/
      data fpppp(12,11),fpppp(12,12)/-9.26525494d-07, 4.93420733d-06/
      data fpppp(12,13),fpppp(12,14)/-5.10883465d-08, 4.34163441d-06/
      data fpppp(12,15),fpppp(12,16)/-1.02881642d-06, 1.24157767d-05/
      data fpppp(12,17),fpppp(12,18)/-2.95138761d-05, 7.18348938d-06/
      data fpppp(12,19) /             2.14159925d-05 /
      data fpppp(13, 1),fpppp(13, 2)/-4.31767522d-06,-1.70685424d-06/
      data fpppp(13, 3),fpppp(13, 4)/ 4.24663382d-06,-1.49921292d-06/
      data fpppp(13, 5),fpppp(13, 6)/-1.11234230d-06, 7.16855357d-07/
      data fpppp(13, 7),fpppp(13, 8)/-2.29170050d-06, 3.75025258d-06/
      data fpppp(13, 9),fpppp(13,10)/-2.07543308d-06, 1.07586323d-06/
      data fpppp(13,11),fpppp(13,12)/-5.15861779d-07,-1.32930404d-06/
      data fpppp(13,13),fpppp(13,14)/ 8.70030311d-07, 3.05561256d-07/
      data fpppp(13,15),fpppp(13,16)/ 6.62188282d-07,-2.65030114d-06/
      data fpppp(13,17),fpppp(13,18)/ 4.27231825d-06,-1.35856869d-06/
      data fpppp(13,19) /            -3.40273953d-06 /
      data fpppp(14, 1),fpppp(14, 2)/ 1.03860710d-06, 3.32078566d-07/
      data fpppp(14, 3),fpppp(14, 4)/-1.44350228d-06, 5.22599985d-07/
      data fpppp(14, 5),fpppp(14, 6)/ 4.52774578d-07,-1.47362172d-07/
      data fpppp(14, 7),fpppp(14, 8)/ 6.13843403d-07,-1.07593355d-06/
      data fpppp(14, 9),fpppp(14,10)/ 4.77683370d-07,-3.27056484d-07/
      data fpppp(14,11),fpppp(14,12)/ 2.16585843d-07, 3.23897629d-07/
      data fpppp(14,13),fpppp(14,14)/-2.23641158d-07,-3.92125434d-08/
      data fpppp(14,15),fpppp(14,16)/-7.82159689d-08, 7.50963973d-07/
      data fpppp(14,17),fpppp(14,18)/-1.13775294d-06, 2.86957450d-07/
      data fpppp(14,19) /             7.64974326d-07 /
      data fpppp(15, 1),fpppp(15, 2)/-2.87117993d-07,-1.12730378d-07/
      data fpppp(15, 3),fpppp(15, 4)/ 2.70821571d-07,-7.77016514d-08/
      data fpppp(15, 5),fpppp(15, 6)/-6.51437264d-08, 3.76588138d-08/
      data fpppp(15, 7),fpppp(15, 8)/-1.52547329d-07, 2.26913014d-07/
      data fpppp(15, 9),fpppp(15,10)/-1.24151718d-07, 6.73365758d-08/
      data fpppp(15,11),fpppp(15,12)/-3.15257766d-08,-8.70836077d-08/
      data fpppp(15,13),fpppp(15,14)/ 5.37670510d-08, 1.71551288d-08/
      data fpppp(15,15),fpppp(15,16)/ 4.69780173d-08,-1.60630656d-07/
      data fpppp(15,17),fpppp(15,18)/ 2.71694729d-07,-8.91900249d-08/
      data fpppp(15,19) /            -2.23443310d-07 /
      data fpppp(16, 1),fpppp(16, 2)/ 1.80914742d-07, 1.33196556d-07/
      data fpppp(16, 3),fpppp(16, 4)/ 1.20079431d-07,-1.41414061d-08/
      data fpppp(16, 5),fpppp(16, 6)/-7.61994259d-08,-7.05591613d-08/
      data fpppp(16, 7),fpppp(16, 8)/ 4.78425427d-08,-7.21379807d-08/
      data fpppp(16, 9),fpppp(16,10)/ 1.33004303d-07,-1.85059145d-09/
      data fpppp(16,11),fpppp(16,12)/-4.31934846d-08, 5.53353131d-08/
      data fpppp(16,13),fpppp(16,14)/-3.41940469d-08,-4.57718452d-08/
      data fpppp(16,15),fpppp(16,16)/-8.34656496d-08, 5.75876014d-08/
      data fpppp(16,17),fpppp(16,18)/-1.21309817d-07, 8.98511214d-08/
      data fpppp(16,19) /             1.89566815d-07 /
 

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
      subroutine a_3(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  2.59787200d-01 ,  2.15296500d-01 /
      data f( 1, 3),f( 1, 4) /  2.43675000d-01 ,  4.59406300d-01 /
      data f( 1, 5),f( 1, 6) /  5.72415200d-01 ,  2.81456600d-01 /
      data f( 1, 7),f( 1, 8) /  1.06082400d-01 ,  5.92888000d-02 /
      data f( 1, 9),f( 1,10) /  5.55679000d-02 ,  3.95559000d-02 /
      data f( 1,11),f( 1,12) /  4.20992000d-02 ,  2.69111000d-02 /
      data f( 1,13),f( 1,14) /  1.84636000d-02 ,  4.98685000d-02 /
      data f( 1,15),f( 1,16) /  1.99933500d-01 ,  5.47771500d-01 /
      data f( 1,17),f( 1,18) /  3.37112300d-01 ,  1.51926600d-01 /
      data f( 1,19) /           1.29267700d-01 /
      data f( 2, 1),f( 2, 2) /  2.17846300d-01 ,  1.79310800d-01 /
      data f( 2, 3),f( 2, 4) /  1.82959800d-01 ,  2.64649900d-01 /
      data f( 2, 5),f( 2, 6) /  2.97097200d-01 ,  1.72977100d-01 /
      data f( 2, 7),f( 2, 8) /  8.10408000d-02 ,  5.37842000d-02 /
      data f( 2, 9),f( 2,10) /  4.48523000d-02 ,  3.30154000d-02 /
      data f( 2,11),f( 2,12) /  3.72400000d-02 ,  3.65541000d-02 /
      data f( 2,13),f( 2,14) /  2.48817000d-02 ,  3.09680000d-02 /
      data f( 2,15),f( 2,16) /  1.24513200d-01 ,  2.80803400d-01 /
      data f( 2,17),f( 2,18) /  2.02153300d-01 ,  1.16592300d-01 /
      data f( 2,19) /           1.82571500d-01 /
      data f( 3, 1),f( 3, 2) /  1.81727900d-01 ,  1.53822000d-01 /
      data f( 3, 3),f( 3, 4) /  1.37702600d-01 ,  1.60514100d-01 /
      data f( 3, 5),f( 3, 6) /  1.65644500d-01 ,  1.20335000d-01 /
      data f( 3, 7),f( 3, 8) /  7.37012000d-02 ,  6.00213000d-02 /
      data f( 3, 9),f( 3,10) /  4.45897000d-02 ,  3.54970000d-02 /
      data f( 3,11),f( 3,12) /  3.89168000d-02 ,  4.28085000d-02 /
      data f( 3,13),f( 3,14) /  3.45472000d-02 ,  3.68300000d-02 /
      data f( 3,15),f( 3,16) /  8.23745000d-02 ,  1.61727800d-01 /
      data f( 3,17),f( 3,18) /  1.39896100d-01 ,  1.03850500d-01 /
      data f( 3,19) /           1.55388200d-01 /
      data f( 4, 1),f( 4, 2) /  1.39450600d-01 ,  1.21123400d-01 /
      data f( 4, 3),f( 4, 4) /  1.02129800d-01 ,  1.04326700d-01 /
      data f( 4, 5),f( 4, 6) /  1.06013500d-01 ,  9.03527000d-02 /
      data f( 4, 7),f( 4, 8) /  5.89081000d-02 ,  5.54943000d-02 /
      data f( 4, 9),f( 4,10) /  3.24253000d-02 ,  2.56667000d-02 /
      data f( 4,11),f( 4,12) /  2.85286000d-02 ,  3.69671000d-02 /
      data f( 4,13),f( 4,14) /  3.59399000d-02 ,  4.30293000d-02 /
      data f( 4,15),f( 4,16) /  6.15174000d-02 ,  1.04400900d-01 /
      data f( 4,17),f( 4,18) /  1.02309100d-01 ,  9.48482000d-02 /
      data f( 4,19) /           8.38441000d-02 /
      data f( 5, 1),f( 5, 2) /  9.11848000d-02 ,  8.98439000d-02 /
      data f( 5, 3),f( 5, 4) /  7.80900000d-02 ,  7.87268000d-02 /
      data f( 5, 5),f( 5, 6) /  7.56200000d-02 ,  7.20340000d-02 /
      data f( 5, 7),f( 5, 8) /  4.76861000d-02 ,  4.08881000d-02 /
      data f( 5, 9),f( 5,10) /  2.32631000d-02 ,  1.83751000d-02 /
      data f( 5,11),f( 5,12) /  2.06682000d-02 ,  2.74436000d-02 /
      data f( 5,13),f( 5,14) /  3.20600000d-02 ,  3.45618000d-02 /
      data f( 5,15),f( 5,16) /  4.79797000d-02 ,  7.87385000d-02 /
      data f( 5,17),f( 5,18) /  8.23894000d-02 ,  8.32975000d-02 /
      data f( 5,19) /           7.79503000d-02 /
      data f( 6, 1),f( 6, 2) /  8.72929000d-02 ,  7.89420000d-02 /
      data f( 6, 3),f( 6, 4) /  6.89145000d-02 ,  5.85416000d-02 /
      data f( 6, 5),f( 6, 6) /  5.86628000d-02 ,  4.39678000d-02 /
      data f( 6, 7),f( 6, 8) /  4.72318000d-02 ,  2.94959000d-02 /
      data f( 6, 9),f( 6,10) /  1.64505000d-02 ,  1.30180000d-02 /
      data f( 6,11),f( 6,12) /  1.48032000d-02 ,  2.00174000d-02 /
      data f( 6,13),f( 6,14) /  3.01613000d-02 ,  3.23345000d-02 /
      data f( 6,15),f( 6,16) /  4.21415000d-02 ,  6.05130000d-02 /
      data f( 6,17),f( 6,18) /  7.28495000d-02 ,  7.68913000d-02 /
      data f( 6,19) /           7.63333000d-02 /
      data f( 7, 1),f( 7, 2) /  7.83178000d-02 ,  8.31032000d-02 /
      data f( 7, 3),f( 7, 4) /  5.82466000d-02 ,  4.78940000d-02 /
      data f( 7, 5),f( 7, 6) /  4.28519000d-02 ,  4.12221000d-02 /
      data f( 7, 7),f( 7, 8) /  3.64454000d-02 ,  2.08890000d-02 /
      data f( 7, 9),f( 7,10) /  1.14454000d-02 ,  9.11470000d-03 /
      data f( 7,11),f( 7,12) /  1.04766000d-02 ,  1.44146000d-02 /
      data f( 7,13),f( 7,14) /  2.20589000d-02 ,  3.45550000d-02 /
      data f( 7,15),f( 7,16) /  3.55161000d-02 ,  5.00380000d-02 /
      data f( 7,17),f( 7,18) /  5.91093000d-02 ,  6.68906000d-02 /
      data f( 7,19) /           7.80235000d-02 /
      data f( 8, 1),f( 8, 2) /  5.65014000d-02 ,  5.98882000d-02 /
      data f( 8, 3),f( 8, 4) /  5.78588000d-02 ,  4.44232000d-02 /
      data f( 8, 5),f( 8, 6) /  4.07746000d-02 ,  4.09213000d-02 /
      data f( 8, 7),f( 8, 8) /  2.57062000d-02 ,  1.45171000d-02 /
      data f( 8, 9),f( 8,10) /  7.81710000d-03 ,  6.29180000d-03 /
      data f( 8,11),f( 8,12) /  7.31950000d-03 ,  1.02475000d-02 /
      data f( 8,13),f( 8,14) /  1.59130000d-02 ,  2.51097000d-02 /
      data f( 8,15),f( 8,16) /  3.11210000d-02 ,  4.08130000d-02 /
      data f( 8,17),f( 8,18) /  4.88164000d-02 ,  5.62489000d-02 /
      data f( 8,19) /           5.92421000d-02 /
      data f( 9, 1),f( 9, 2) /  2.79368000d-02 ,  2.95137000d-02 /
      data f( 9, 3),f( 9, 4) /  3.16920000d-02 ,  3.14530000d-02 /
      data f( 9, 5),f( 9, 6) /  2.70952000d-02 ,  1.96180000d-02 /
      data f( 9, 7),f( 9, 8) /  1.20618000d-02 ,  6.53860000d-03 /
      data f( 9, 9),f( 9,10) /  3.40230000d-03 ,  2.83900000d-03 /
      data f( 9,11),f( 9,12) /  3.40710000d-03 ,  4.96220000d-03 /
      data f( 9,13),f( 9,14) /  7.90910000d-03 ,  1.26201000d-02 /
      data f( 9,15),f( 9,16) /  1.87835000d-02 ,  2.47076000d-02 /
      data f( 9,17),f( 9,18) /  2.86216000d-02 ,  3.06440000d-02 /
      data f( 9,19) /           3.16196000d-02 /
      data f(10, 1),f(10, 2) /  1.27829000d-02 ,  1.35073000d-02 /
      data f(10, 3),f(10, 4) /  1.45438000d-02 ,  1.43818000d-02 /
      data f(10, 5),f(10, 6) /  1.22396000d-02 ,  8.67420000d-03 /
      data f(10, 7),f(10, 8) /  5.09350000d-03 ,  2.52540000d-03 /
      data f(10, 9),f(10,10) /  1.26370000d-03 ,  1.14320000d-03 /
      data f(10,11),f(10,12) /  1.44900000d-03 ,  2.22480000d-03 /
      data f(10,13),f(10,14) /  3.67270000d-03 ,  5.93400000d-03 /
      data f(10,15),f(10,16) /  8.85490000d-03 ,  1.17413000d-02 /
      data f(10,17),f(10,18) /  1.37710000d-02 ,  1.50384000d-02 /
      data f(10,19) /           1.55877000d-02 /
      data f(11, 1),f(11, 2) /  3.97730000d-03 ,  4.26310000d-03 /
      data f(11, 3),f(11, 4) /  4.69490000d-03 ,  4.64780000d-03 /
      data f(11, 5),f(11, 6) /  3.83640000d-03 ,  2.51740000d-03 /
      data f(11, 7),f(11, 8) /  1.21500000d-03 ,  4.03500000d-04 /
      data f(11, 9),f(11,10) /  2.37600000d-04 ,  2.37200000d-04 /
      data f(11,11),f(11,12) /  3.63200000d-04 ,  6.56300000d-04 /
      data f(11,13),f(11,14) /  1.19290000d-03 ,  2.00850000d-03 /
      data f(11,15),f(11,16) /  3.02590000d-03 ,  4.02540000d-03 /
      data f(11,17),f(11,18) /  4.82340000d-03 ,  5.49760000d-03 /
      data f(11,19) /           5.75820000d-03 /
      data f(12, 1),f(12, 2) /  6.03800000d-04 ,  7.19800000d-04 /
      data f(12, 3),f(12, 4) /  9.13800000d-04 ,  9.22300000d-04 /
      data f(12, 5),f(12, 6) /  6.48600000d-04 ,  2.17800000d-04 /
      data f(12, 7),f(12, 8) / -1.15400000d-04 , -1.29400000d-04 /
      data f(12, 9),f(12,10) / -6.88000000d-05 , -3.11000000d-05 /
      data f(12,11),f(12,12) / -1.61000000d-05 ,  5.37000000d-05 /
      data f(12,13),f(12,14) /  2.19300000d-04 ,  4.58700000d-04 /
      data f(12,15),f(12,16) /  7.42000000d-04 ,  9.96400000d-04 /
      data f(12,17),f(12,18) /  1.35580000d-03 ,  1.66190000d-03 /
      data f(12,19) /           1.75660000d-03 /
      data f(13, 1),f(13, 2) / -5.60500000d-04 , -5.63500000d-04 /
      data f(13, 3),f(13, 4) / -5.07700000d-04 , -4.65500000d-04 /
      data f(13, 5),f(13, 6) / -3.95900000d-04 , -3.18600000d-04 /
      data f(13, 7),f(13, 8) / -2.54000000d-04 , -1.96800000d-04 /
      data f(13, 9),f(13,10) / -1.51800000d-04 , -1.26100000d-04 /
      data f(13,11),f(13,12) / -1.16700000d-04 , -1.16800000d-04 /
      data f(13,13),f(13,14) / -1.20700000d-04 , -1.25000000d-04 /
      data f(13,15),f(13,16) / -1.28700000d-04 , -1.30700000d-04 /
      data f(13,17),f(13,18) / -1.16900000d-04 , -1.05100000d-04 /
      data f(13,19) /          -1.03900000d-04 /
      data f(14, 1),f(14, 2) / -1.64600000d-04 , -1.63500000d-04 /
      data f(14, 3),f(14, 4) / -1.57400000d-04 , -1.44400000d-04 /
      data f(14, 5),f(14, 6) / -1.27200000d-04 , -1.10700000d-04 /
      data f(14, 7),f(14, 8) / -9.73000000d-05 , -8.07000000d-05 /
      data f(14, 9),f(14,10) / -6.43000000d-05 , -5.20000000d-05 /
      data f(14,11),f(14,12) / -4.66000000d-05 , -4.97000000d-05 /
      data f(14,13),f(14,14) / -5.71000000d-05 , -6.76000000d-05 /
      data f(14,15),f(14,16) / -8.10000000d-05 , -9.69000000d-05 /
      data f(14,17),f(14,18) / -1.11500000d-04 , -1.18900000d-04 /
      data f(14,19) /          -1.21200000d-04 /
      data f(15, 1),f(15, 2) / -5.68000000d-05 , -5.49000000d-05 /
      data f(15, 3),f(15, 4) / -5.06000000d-05 , -4.51000000d-05 /
      data f(15, 5),f(15, 6) / -3.97000000d-05 , -3.47000000d-05 /
      data f(15, 7),f(15, 8) / -2.93000000d-05 , -2.18000000d-05 /
      data f(15, 9),f(15,10) / -1.43000000d-05 , -9.70000000d-06 /
      data f(15,11),f(15,12) / -8.70000000d-06 , -1.10000000d-05 /
      data f(15,13),f(15,14) / -1.56000000d-05 , -2.14000000d-05 /
      data f(15,15),f(15,16) / -2.77000000d-05 , -3.47000000d-05 /
      data f(15,17),f(15,18) / -4.19000000d-05 , -4.69000000d-05 /
      data f(15,19) /          -4.86000000d-05 /
      data f(16, 1),f(16, 2) / -4.70000000d-06 , -4.80000000d-06 /
      data f(16, 3),f(16, 4) / -4.60000000d-06 , -4.60000000d-06 /
      data f(16, 5),f(16, 6) / -4.50000000d-06 , -3.50000000d-06 /
      data f(16, 7),f(16, 8) / -1.00000000d-06 ,  2.30000000d-06 /
      data f(16, 9),f(16,10) /  5.10000000d-06 ,  6.90000000d-06 /
      data f(16,11),f(16,12) /  7.20000000d-06 ,  6.30000000d-06 /
      data f(16,13),f(16,14) /  5.10000000d-06 ,  3.00000000d-06 /
      data f(16,15),f(16,16) /  1.00000000d-06 , -1.50000000d-06 /
      data f(16,17),f(16,18) / -3.10000000d-06 , -3.80000000d-06 /
      data f(16,19) /          -4.40000000d-06 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 4.54543554d-01,-2.77943597d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 6.81307183d-01, 4.52291194d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 5.25530125d-01, 2.84093082d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.30165699d+00,-5.74846473d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 5.37504000d+00,-6.70488893d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.19410364d+00, 3.15635219d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.03564928d+00, 1.01454418d-03/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 8.36674375d-01, 5.00307087d-04/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 7.77257828d-01,-4.31410531d-04/
      data fpp( 1,10,1),fpp( 1,10,2)/ 7.17105728d-01, 4.87869036d-04/
      data fpp( 1,11,1),fpp( 1,11,2)/ 5.89822407d-01,-4.06747615d-04/
      data fpp( 1,12,1),fpp( 1,12,2)/ 1.14939074d-01, 7.52374235d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.56842362d-01, 5.10233921d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.23113422d+00, 2.74970893d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.13898510d+00, 5.50948851d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.79244645d+00,-1.04465449d-02/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.97678584d+00, 2.76685917d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/ 1.02848371d+00, 9.07518237d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/-3.04470909d+00, 3.35467588d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 1.26670392d-01, 1.23238064d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 3.09988134d-01, 3.01926872d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.96969750d-01, 1.20012445d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.32684103d+00,-4.19958668d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.64199749d+00,-2.47485778d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.45847773d+00, 9.25345770d-04/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 5.14126431d-01, 7.04502695d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 3.32568750d-01, 1.37425450d-04/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 3.47199344d-01,-1.54722494d-04/
      data fpp( 2,10,1),fpp( 2,10,2)/ 3.09146044d-01, 3.07164526d-04/
      data fpp( 2,11,1),fpp( 2,11,2)/ 2.40605186d-01,-1.10245608d-04/
      data fpp( 2,12,1),fpp( 2,12,2)/-4.86631490d-02,-1.60812092d-04/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.05880277d-01, 9.43039754d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 6.16179051d-01, 8.49118190d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 8.18149796d-01, 1.75675726d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.81423459d+00,-4.11144725d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.90064331d+00, 5.92613720d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 5.80120090d-01, 1.32633837d-03/
      data fpp( 2,19,1),fpp( 2,19,2)/-1.75342931d+00, 3.19444482d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-8.78501213d-02,-1.18594460d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-3.46724719d-01, 4.78959190d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.05290874d-01, 6.34200783d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 9.84068896d-01,-2.48845053d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.63676502d+00,-6.99686572d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 3.47595460d-01, 2.11973423d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-4.36855007d-01, 5.35439203d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-4.05694375d-01,-1.85720155d-04/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-5.98105203d-01, 1.02339418d-04/
      data fpp( 3,10,1),fpp( 3,10,2)/-6.00374903d-01, 1.56696484d-04/
      data fpp( 3,11,1),fpp( 3,11,2)/-5.71843153d-01, 2.16246462d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/-4.28576479d-01,-2.14881069d-04/
      data fpp( 3,13,1),fpp( 3,13,2)/-2.93253469d-01, 1.08719629d-04/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.85245717d-02, 4.12648552d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 5.80655714d-01, 8.36388164d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.13449018d+00,-1.72967321d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 3.25910906d-01, 1.12046717d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 3.99109352d-02, 8.32020522d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/-2.01463865d+00, 1.91571124d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-6.99104907d-01,-2.06597869d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-4.55925855d-03,-4.39342627d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 2.34526753d-01, 3.42350919d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 9.29143385d-01,-5.40394147d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 5.84197429d-01,-1.56799260d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 5.50110436d-01,-3.59619543d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 1.15268597d-01, 6.48249434d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-3.24406251d-01,-5.51530194d-04/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.59951469d-01, 3.78559341d-04/
      data fpp( 4,10,1),fpp( 4,10,2)/ 2.45568567d-01, 1.59168303d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/ 2.37017424d-01, 1.35003338d-04/
      data fpp( 4,12,1),fpp( 4,12,2)/-5.14009368d-02,-2.21334182d-04/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.73786399d-01, 1.82391391d-04/
      data fpp( 4,14,1),fpp( 4,14,2)/-6.39682338d-01,-2.12353808d-05/
      data fpp( 4,15,1),fpp( 4,15,2)/ 5.14673481d-02, 5.86472132d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 9.10109706d-01,-8.60929149d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 4.96243063d-01, 1.58726463d-04/
      data fpp( 4,18,1),fpp( 4,18,2)/-1.78838831d-01,-9.61227038d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 3.15786393d+00, 1.31723519d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 1.98599475d+00,-3.12526928d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 5.77826753d-01,-1.43410143d-04/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 5.86552113d-01, 2.61387502d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-1.12517436d-01,-1.58697865d-04/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 4.12070266d-01, 1.48787958d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-7.98497202d-01,-4.65205968d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 5.11445617d-01, 4.66321915d-04/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.91439377d-01,-3.47087693d-04/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 8.62932893d-03, 2.72408858d-04/
      data fpp( 5,10,1),fpp( 5,10,2)/-1.09436622d-03, 2.16722627d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/ 2.94345506d-03, 7.17680917d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/ 8.18652259d-02,-3.98066293d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/ 1.97509066d-01,-4.20815745d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 3.40184780d-01, 8.12569274d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 3.11384894d-01, 3.72019865d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/-2.52539988d-02,-5.28882386d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 3.39211842d-01, 1.17035681d-04/
      data fpp( 5,18,1),fpp( 5,18,2)/ 2.93184388d-01,-1.03828337d-04/
      data fpp( 5,19,1),fpp( 5,19,2)/-7.69272064d-01,-7.70403313d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/-5.88789084d-01,-4.04838546d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 7.49892246d-01, 4.04570916d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-3.51090204d-01,-7.62949821d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 3.33131359d-01, 2.80410219d-04/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-2.17033494d-01,-4.15699894d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.18175337d+00, 4.93417358d-04/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-5.45896067d-01,-4.80429537d-04/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 4.07487425d-02, 1.68306788d-04/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 5.79712157d-02, 8.86323831d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/ 4.89838976d-02, 5.39376792d-05/
      data fpp( 6,11,1),fpp( 6,11,2)/ 5.05187553d-02, 8.67890019d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/ 3.85350332d-02, 1.17086720d-04/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.19069865d-01,-1.81243780d-04/
      data fpp( 6,14,1),fpp( 6,14,2)/ 2.14973219d-01, 1.29646402d-04/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.42081923d-01, 1.20686174d-04/
      data fpp( 6,16,1),fpp( 6,16,2)/ 3.06441289d-01,-9.85210982d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/-2.96120432d-01,-8.87017812d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/-2.22223721d-01,-4.43537768d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 5.60744325d-01,-9.87111161d-06/
      data fpp( 7, 1,1),fpp( 7, 1,2)/-3.93318412d-01,-7.14890058d-04/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.31793074d+00,-3.42399884d-04/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 5.93948704d-01, 3.05969594d-04/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 2.10632002d-01,-1.12384925d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 6.28008708d-01, 5.76143760d-05/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.30441289d-01,-1.44810116d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 1.22323651d-01,-1.88504330d-04/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 6.33606531d-02, 1.21716330d-04/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 3.06108081d-02, 6.84070082d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/ 2.32287757d-02, 3.14296370d-05/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.57415236d-02, 2.74304439d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/ 3.75046414d-02, 1.34145873d-05/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.48215395d-01, 1.41289207d-04/
      data fpp( 7,14,1),fpp( 7,14,2)/-5.32907655d-01,-2.87463414d-04/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.38862799d-01, 3.16464451d-04/
      data fpp( 7,16,1),fpp( 7,16,2)/-3.79361590d-02,-1.64746389d-04/
      data fpp( 7,17,1),fpp( 7,17,2)/ 2.15224888d-01, 1.54851047d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/ 5.65354962d-02, 2.54059701d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/-9.77625236d-01, 8.39870150d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 2.35867731d-01,-7.80903700d-06/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 4.15400708d-01,-2.70679260d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-4.82689610d-01,-2.08891259d-04/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-9.91393647d-02, 1.78260962d-04/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-2.34961339d-01, 8.30674107d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-2.93253218d-01,-2.82812605d-04/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 6.36814627d-02, 1.26475009d-04/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 4.10586450d-02, 1.84725706d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 2.61055517d-02, 6.89807090d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/ 2.01609997d-02, 1.60865932d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/ 2.19401503d-02, 1.98529180d-05/
      data fpp( 8,12,1),fpp( 8,12,2)/ 2.68014011d-02, 1.85197346d-05/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.96832838d-02, 7.03181435d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/ 1.66787399d-01,-8.79203086d-05/
      data fpp( 8,15,1),fpp( 8,15,2)/-7.88242721d-02, 9.02390910d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/ 3.28033466d-02,-5.21940553d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/-4.76841177d-02, 1.72211303d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/-1.00068264d-01,-5.09444658d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 2.79016621d-01,-7.98012671d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 5.41135135d-02, 3.58907965d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.48444952d-02, 6.63440700d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.98924479d-01,-2.63444245d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-3.39704067d-02,-4.62947091d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 3.36996643d-02,-3.56047393d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.68666549d-01, 1.54966617d-06/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 4.15687863d-02, 2.46660746d-05/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 2.38424884d-02, 2.17660355d-05/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.29454408d-02, 3.14837834d-05/
      data fpp( 9,10,1),fpp( 9,10,2)/ 1.01401130d-02, 6.67883077d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.13762873d-02, 9.68489348d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.51772261d-02, 1.38015953d-05/
      data fpp( 9,13,1),fpp( 9,13,2)/ 2.76387010d-02, 1.86167253d-05/
      data fpp( 9,14,1),fpp( 9,14,2)/ 6.12912906d-03, 1.75775037d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/ 3.40176668d-02,-1.78273993d-06/
      data fpp( 9,16,1),fpp( 9,16,2)/ 8.48053957d-03,-2.48045439d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/ 5.01024093d-02,-1.96050843d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.09880794d-01,-1.02711188d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 2.45240056d-02,-2.11844062d-06/
      data fpp(10, 1,1),fpp(10, 1,2)/ 5.05794653d-02, 1.80466121d-05/
      data fpp(10, 2,1),fpp(10, 2,2)/ 6.40250612d-02, 3.48177579d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.51891954d-02,-1.32477153d-05/
      data fpp(10, 4,1),fpp(10, 4,2)/ 8.12334914d-02,-2.24009147d-05/
      data fpp(10, 5,1),fpp(10, 5,2)/ 5.60551823d-02,-1.59606259d-05/
      data fpp(10, 6,1),fpp(10, 6,2)/ 7.06827068d-03, 8.51418314d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/ 2.03971421d-02, 1.16369526d-05/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.22701515d-02, 1.33567711d-05/
      data fpp(10, 9,1),fpp(10, 9,2)/ 7.47018524d-03, 1.33199629d-05/
      data fpp(10,10,1),fpp(10,10,2)/ 5.16604810d-03, 1.83537719d-06/
      data fpp(10,11,1),fpp(10,11,2)/ 5.84095050d-03, 4.91652830d-06/
      data fpp(10,12,1),fpp(10,12,2)/ 8.03594452d-03, 6.69850961d-06/
      data fpp(10,13,1),fpp(10,13,2)/ 1.10431623d-02, 8.61543326d-06/
      data fpp(10,14,1),fpp(10,14,2)/ 2.63273344d-02, 7.64375733d-06/
      data fpp(10,15,1),fpp(10,15,2)/ 3.30873547d-02, 3.85537399d-07/
      data fpp(10,16,1),fpp(10,16,2)/ 5.09907451d-02,-1.12559069d-05/
      data fpp(10,17,1),fpp(10,17,2)/ 4.76819805d-02,-6.76390967d-06/
      data fpp(10,18,1),fpp(10,18,2)/ 3.55188395d-02,-7.42645438d-06/
      data fpp(10,19,1),fpp(10,19,2)/ 5.75348567d-02,-6.61627281d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 1.79057142d-02, 7.62037750d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/ 1.59653836d-02, 1.63724500d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 2.82517135d-02,-5.40935749d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 1.32557563d-02,-8.73381503d-06/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.52328124d-02,-5.51338240d-06/
      data fpp(11, 6,1),fpp(11, 6,2)/ 2.01717860d-02, 3.31344615d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 9.28025929d-03, 5.18400394d-06/
      data fpp(11, 8,1),fpp(11, 8,2)/ 6.22386392d-03, 8.38663963d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/ 2.28258055d-03, 5.43752647d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 2.42013641d-03, 1.52161026d-06/
      data fpp(11,11,1),fpp(11,11,2)/ 2.55534837d-03, 1.49212143d-06/
      data fpp(11,12,1),fpp(11,12,2)/ 3.40681886d-03, 2.53590401d-06/
      data fpp(11,13,1),fpp(11,13,2)/ 5.71045492d-03, 2.97426253d-06/
      data fpp(11,14,1),fpp(11,14,2)/ 6.68929303d-03, 2.30704586d-06/
      data fpp(11,15,1),fpp(11,15,2)/ 1.16333895d-02,-9.44459914d-08/
      data fpp(11,16,1),fpp(11,16,2)/ 1.34562861d-02,-3.00326190d-06/
      data fpp(11,17,1),fpp(11,17,2)/ 1.90385429d-02, 1.74935882d-08/
      data fpp(11,18,1),fpp(11,18,2)/ 2.34163430d-02,-4.49471245d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 1.83043115d-02,-6.85464377d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 8.16807792d-03, 3.36744135d-06/
      data fpp(12, 2,1),fpp(12, 2,2)/ 8.93500427d-03, 8.75117308d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/ 7.43115057d-03,-2.18791058d-06/
      data fpp(12, 4,1),fpp(12, 4,2)/ 9.94748356d-03,-3.25347500d-06/
      data fpp(12, 5,1),fpp(12, 5,2)/ 8.18316824d-03,-1.73018943d-06/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.81738524d-03, 7.48232708d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 3.63622070d-03, 4.59325860d-06/
      data fpp(12, 8,1),fpp(12, 8,2)/ 9.70392828d-04, 3.07329066d-08/
      data fpp(12, 9,1),fpp(12, 9,2)/ 6.72292582d-04,-2.40190223d-07/
      data fpp(12,10,1),fpp(12,10,2)/ 4.58206277d-04,-4.43972015d-07/
      data fpp(12,11,1),fpp(12,11,2)/ 8.93656041d-04, 6.54078282d-07/
      data fpp(12,12,1),fpp(12,12,2)/ 1.51838005d-03, 1.11565889d-06/
      data fpp(12,13,1),fpp(12,13,2)/ 2.26381803d-03, 6.31286165d-07/
      data fpp(12,14,1),fpp(12,14,2)/ 3.93229350d-03, 7.87196453d-07/
      data fpp(12,15,1),fpp(12,15,2)/ 5.46148720d-03,-1.14607198d-06/
      data fpp(12,16,1),fpp(12,16,2)/ 7.66971053d-03, 2.06309146d-06/
      data fpp(12,17,1),fpp(12,17,2)/ 7.68384801d-03,-8.06293854d-07/
      data fpp(12,18,1),fpp(12,18,2)/ 7.73818840d-03,-2.03591604d-06/
      data fpp(12,19,1),fpp(12,19,2)/ 9.11749740d-03,-3.73404198d-06/
      data fpp(13, 1,1),fpp(13, 1,2)/ 3.91091526d-05, 1.25849627d-06/
      data fpp(13, 2,1),fpp(13, 2,2)/ 3.20953572d-05, 6.95007463d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 4.24891532d-04,-5.10526121d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-9.11288059d-05, 5.31097020d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/-1.79310911d-04, 3.01380416d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/-1.61248734d-04,-1.89649186d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/-4.15591740d-04,-3.35412974d-08/
      data fpp(13, 8,1),fpp(13, 8,2)/-3.27104437d-05,-1.20185624d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/ 2.06319809d-05,-2.17716205d-07/
      data fpp(13,10,1),fpp(13,10,2)/ 6.49129665d-05,-1.66949554d-07/
      data fpp(13,11,1),fpp(13,11,2)/-1.06423046d-05,-9.24855786d-08/
      data fpp(13,12,1),fpp(13,12,2)/-5.03495763d-05,-3.31081316d-08/
      data fpp(13,13,1),fpp(13,13,2)/-3.48154575d-06,-3.08189492d-09/
      data fpp(13,14,1),fpp(13,14,2)/-4.61270067d-05, 2.14357113d-08/
      data fpp(13,15,1),fpp(13,15,2)/-1.85563529d-05,-4.66609503d-08/
      data fpp(13,16,1),fpp(13,16,2)/-1.51874627d-04, 2.67208090d-07/
      data fpp(13,17,1),fpp(13,17,2)/ 2.04184518d-04,-7.41714088d-08/
      data fpp(13,18,1),fpp(13,18,2)/ 5.03663271d-04,-9.05224546d-08/
      data fpp(13,19,1),fpp(13,19,2)/ 3.51552072d-04,-1.99738773d-07/
      data fpp(14, 1,1),fpp(14, 1,2)/-1.14616416d-04, 3.26776565d-08/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.13888209d-04, 4.66446869d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-2.00299882d-04, 8.07435958d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-5.53053618d-05, 4.43809300d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/-1.71013874d-05,-6.26731562d-09/
      data fpp(14, 6,1),fpp(14, 6,2)/-3.89641889d-06,-6.13116675d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/ 7.95148706d-05, 6.55139856d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/-1.07150830d-05,-8.74427484d-09/
      data fpp(14, 9,1),fpp(14, 9,2)/-1.77922338d-05,-4.25368862d-08/
      data fpp(14,10,1),fpp(14,10,2)/-2.76920380d-05,-6.71081802d-08/
      data fpp(14,11,1),fpp(14,11,2)/-7.95110660d-06,-1.03030393d-07/
      data fpp(14,12,1),fpp(14,12,2)/ 4.00870442d-06,-3.07702481d-08/
      data fpp(14,13,1),fpp(14,13,2)/-6.06437736d-06,-3.18886146d-08/
      data fpp(14,14,1),fpp(14,14,2)/ 9.43427187d-06,-2.76752935d-08/
      data fpp(14,15,1),fpp(14,15,2)/ 8.57545981d-06,-3.14102114d-08/
      data fpp(14,16,1),fpp(14,16,2)/ 5.27686173d-05, 3.31613914d-09/
      data fpp(14,17,1),fpp(14,17,2)/-2.82775618d-05, 9.61456548d-08/
      data fpp(14,18,1),fpp(14,18,2)/-9.97840163d-05, 4.41012415d-08/
      data fpp(14,19,1),fpp(14,19,2)/-5.78549143d-05, 3.34493793d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/-1.27934905d-05, 3.61343196d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/-1.36425223d-05, 2.37313608d-08/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.10579975d-05, 1.29402372d-08/
      data fpp(15, 4,1),fpp(15, 4,2)/-2.03497468d-05,-3.49230955d-09/
      data fpp(15, 5,1),fpp(15, 5,2)/-2.40835391d-05,-4.97099900d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/-2.10155905d-05,-6.23694450d-10/
      data fpp(15, 7,1),fpp(15, 7,2)/-3.55177427d-05, 3.14657768d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-1.02292242d-05, 7.60587252d-10/
      data fpp(15, 9,1),fpp(15, 9,2)/-5.71304584d-06,-3.45081258d-08/
      data fpp(15,10,1),fpp(15,10,2)/-1.84481459d-06,-3.67280840d-08/
      data fpp(15,11,1),fpp(15,11,2)/-5.85326904d-06,-3.45795381d-08/
      data fpp(15,12,1),fpp(15,12,2)/-8.28524138d-06,-2.29537636d-08/
      data fpp(15,13,1),fpp(15,13,2)/-5.41094482d-06,-1.16054077d-08/
      data fpp(15,14,1),fpp(15,14,2)/-8.41008078d-06,-2.62460579d-09/
      data fpp(15,15,1),fpp(15,15,2)/-7.34548634d-06,-7.89616917d-09/
      data fpp(15,16,1),fpp(15,16,2)/-1.65998427d-05,-7.79071754d-09/
      data fpp(15,17,1),fpp(15,17,2)/ 5.22572882d-06, 2.70590393d-08/
      data fpp(15,18,1),fpp(15,18,2)/ 2.41727937d-05, 3.15545602d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 1.47175849d-05, 4.47227199d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 4.59892452d-05, 7.72126553d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/ 4.62005469d-05, 3.55746893d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/ 2.05399298d-05,-3.95114126d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/ 3.24027305d-05, 2.47096116d-10/
      data fpp(16, 5,1),fpp(16, 5,2)/ 2.76460553d-05, 8.96275680d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/ 1.74495810d-05, 1.79018767d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/ 1.62492285d-05, 9.42973642d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/-1.59253078d-06,-7.62082238d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/-6.16169137d-06,-8.94644689d-09/
      data fpp(16,10,1),fpp(16,10,2)/-6.62259270d-06,-1.65933901d-08/
      data fpp(16,11,1),fpp(16,11,2)/-2.48836548d-06,-1.46799929d-08/
      data fpp(16,12,1),fpp(16,12,2)/-2.22166503d-06, 3.31336154d-09/
      data fpp(16,13,1),fpp(16,13,2)/-5.62059902d-06,-1.65734533d-08/
      data fpp(16,14,1),fpp(16,14,2)/-8.18924533d-06, 8.98045158d-09/
      data fpp(16,15,1),fpp(16,15,2)/-1.53986854d-05,-1.33483530d-08/
      data fpp(16,16,1),fpp(16,16,2)/-1.99129358d-05, 1.44129606d-08/
      data fpp(16,17,1),fpp(16,17,2)/-4.23007216d-05, 9.69651060d-09/
      data fpp(16,18,1),fpp(16,18,2)/-5.73199683d-05, 8.00996970d-10/
      data fpp(16,19,1),fpp(16,19,2)/-5.36220067d-05,-6.90049849d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-3.46098620d-02,-8.54538868d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 4.58389755d-02, 1.10372174d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/-9.24184930d-02, 5.33110874d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 5.23064102d-04, 2.16542272d-03/
      data fpppp( 1, 9),fpppp( 1,10)/-8.11253214d-04, 1.03545695d-03/
      data fpppp( 1,11),fpppp( 1,12)/-7.35844788d-03, 7.54233391d-03/
      data fpppp( 1,13),fpppp( 1,14)/ 2.01963094d-02,-5.03842570d-02/
      data fpppp( 1,15),fpppp( 1,16)/ 1.23354259d-01,-1.58296152d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 6.16830307d-02,-3.63944626d-02/
      data fpppp( 1,19) /            -4.35986199d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/-1.88485048d-02,-3.97758993d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 2.89786970d-02,-1.36381818d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/-6.04063131d-02, 3.30684966d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 2.48243498d-03, 2.76938033d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-1.78865983d-03, 1.22422537d-03/
      data fpppp( 2,11),fpppp( 2,12)/-4.93749508d-03, 5.28210626d-03/
      data fpppp( 2,13),fpppp( 2,14)/ 1.04377757d-02,-2.56878883d-02/
      data fpppp( 2,15),fpppp( 2,16)/ 7.38140958d-02,-1.01921652d-01/
      data fpppp( 2,17),fpppp( 2,18)/ 3.92919469d-02,-1.96620522d-02/
      data fpppp( 2,19) /            -2.14253090d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 1.35517535d-02, 8.90575391d-03/
      data fpppp( 3, 3),fpppp( 3, 4)/-5.21357708d-04, 6.78542268d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-3.41852470d-02, 1.34436243d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.06938952d-02,-7.28253913d-03/
      data fpppp( 3, 9),fpppp( 3,10)/ 5.02197372d-03,-1.39688804d-03/
      data fpppp( 3,11),fpppp( 3,12)/ 2.41366543d-03,-1.37367822d-03/
      data fpppp( 3,13),fpppp( 3,14)/ 2.60442756d-03, 1.54326991d-03/
      data fpppp( 3,15),fpppp( 3,16)/ 6.24367886d-03,-2.70157862d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 2.00746422d-02,-2.19280245d-02/
      data fpppp( 3,19) /            -3.84755212d-02 /
      data fpppp( 4, 1),fpppp( 4, 2)/-1.22841421d-02,-7.31531005d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 1.42178041d-02,-2.22240692d-02/
      data fpppp( 4, 5),fpppp( 4, 6)/ 1.23047174d-02,-8.34326272d-03/
      data fpppp( 4, 7),fpppp( 4, 8)/-2.97695725d-03, 1.99611112d-02/
      data fpppp( 4, 9),fpppp( 4,10)/-1.54255333d-02, 5.81658498d-03/
      data fpppp( 4,11),fpppp( 4,12)/-7.49090108d-03, 7.35498622d-03/
      data fpppp( 4,13),fpppp( 4,14)/-1.19670699d-02, 1.99026646d-02/
      data fpppp( 4,15),fpppp( 4,16)/ 1.77914889d-03,-1.69696999d-02/
      data fpppp( 4,17),fpppp( 4,18)/-1.02508895d-02, 4.23003428d-02/
      data fpppp( 4,19) /             8.17565974d-02 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 3.28959672d-02, 1.92086314d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.47168917d-02, 3.71912409d-02/
      data fpppp( 5, 5),fpppp( 5, 6)/-5.06286367d-02, 6.12139957d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/-4.29967287d-02, 1.29759754d-02/
      data fpppp( 5, 9),fpppp( 5,10)/-6.75401308d-04, 1.10811033d-04/
      data fpppp( 5,11),fpppp( 5,12)/ 1.05784816d-03, 1.50833288d-04/
      data fpppp( 5,13),fpppp( 5,14)/ 5.42142845d-04,-6.97492260d-04/
      data fpppp( 5,15),fpppp( 5,16)/-8.04070979d-03, 1.43899910d-02/
      data fpppp( 5,17),fpppp( 5,18)/-7.45297021d-03,-9.20770790d-03/
      data fpppp( 5,19) /            -1.67019380d-02 /
      data fpppp( 6, 1),fpppp( 6, 2)/-6.26039068d-02,-3.24794557d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/ 4.61419026d-02,-4.49759140d-02/
      data fpppp( 6, 5),fpppp( 6, 6)/ 5.96985686d-02,-7.68812572d-02/
      data fpppp( 6, 7),fpppp( 6, 8)/ 6.02402819d-02,-2.52222155d-02/
      data fpppp( 6, 9),fpppp( 6,10)/ 6.48323988d-03,-2.28333150d-03/
      data fpppp( 6,11),fpppp( 6,12)/ 3.28141665d-03,-1.16534499d-02/
      data fpppp( 6,13),fpppp( 6,14)/ 2.25951124d-02,-2.52281207d-02/
      data fpppp( 6,15),fpppp( 6,16)/ 2.48514768d-02,-2.58430853d-02/
      data fpppp( 6,17),fpppp( 6,18)/ 1.54557684d-02, 4.60751791d-03/
      data fpppp( 6,19) /             8.65844008d-03 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 7.55522685d-02, 3.66239743d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/-5.18586597d-02, 3.30988959d-02/
      data fpppp( 7, 5),fpppp( 7, 6)/-3.24953193d-02, 2.63327789d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.21629002d-02, 3.61514566d-03/
      data fpppp( 7, 9),fpppp( 7,10)/-7.24893264d-04, 8.06496145d-04/
      data fpppp( 7,11),fpppp( 7,12)/-1.90740449d-03, 7.37814403d-03/
      data fpppp( 7,13),fpppp( 7,14)/-2.16683134d-02, 3.17850815d-02/
      data fpppp( 7,15),fpppp( 7,16)/-2.42984024d-02, 1.44943634d-02/
      data fpppp( 7,17),fpppp( 7,18)/-7.88145097d-03,-7.67958582d-03/
      data fpppp( 7,19) /            -1.39284862d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/-3.23820835d-02,-1.47498090d-02/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.67239218d-02,-1.52474445d-02/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.10352308d-03, 7.48515797d-03/
      data fpppp( 8, 7),fpppp( 8, 8)/-8.13056137d-03, 2.26363757d-03/
      data fpppp( 8, 9),fpppp( 8,10)/-4.63805443d-04, 1.32096680d-04/
      data fpppp( 8,11),fpppp( 8,12)/ 3.98840881d-04,-1.54253419d-03/
      data fpppp( 8,13),fpppp( 8,14)/ 5.05253381d-03,-9.41426708d-03/
      data fpppp( 8,15),fpppp( 8,16)/ 9.04158730d-03,-5.31772468d-03/
      data fpppp( 8,17),fpppp( 8,18)/ 7.02406456d-04, 4.19429795d-03/
      data fpppp( 8,19) /             8.40854362d-03 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.97778348d-03, 3.55138050d-03/
      data fpppp( 9, 3),fpppp( 9, 4)/-8.78236538d-03, 6.55958890d-03/
      data fpppp( 9, 5),fpppp( 9, 6)/ 5.77907167d-04,-4.83340872d-03/
      data fpppp( 9, 7),fpppp( 9, 8)/ 3.03184881d-03,-7.31698620d-04/
      data fpppp( 9, 9),fpppp( 9,10)/ 3.04700684d-04,-1.60092260d-06/
      data fpppp( 9,11),fpppp( 9,12)/-5.58068759d-05, 3.78714299d-04/
      data fpppp( 9,13),fpppp( 9,14)/-9.39418158d-04, 1.34069552d-03/
      data fpppp( 9,15),fpppp( 9,16)/-1.45947736d-03, 1.29167401d-03/
      data fpppp( 9,17),fpppp( 9,18)/ 3.22321147d-04,-1.49156772d-03/
      data fpppp( 9,19) /            -3.06416059d-03 /
      data fpppp(10, 1),fpppp(10, 2)/-1.83695108d-03,-8.37774151d-04/
      data fpppp(10, 3),fpppp(10, 4)/ 2.05115999d-03,-1.67405611d-03/
      data fpppp(10, 5),fpppp(10, 6)/-2.28291872d-04, 1.15870745d-03/
      data fpppp(10, 7),fpppp(10, 8)/-6.67590934d-04, 2.24304562d-04/
      data fpppp(10, 9),fpppp(10,10)/-3.00058521d-05, 4.54685937d-05/
      data fpppp(10,11),fpppp(10,12)/ 2.68738494d-05,-6.17584941d-05/
      data fpppp(10,13),fpppp(10,14)/ 2.68893554d-04,-2.77198464d-04/
      data fpppp(10,15),fpppp(10,16)/ 3.28451200d-04,-3.68004137d-04/
      data fpppp(10,17),fpppp(10,18)/-1.29163949d-04, 3.53397349d-04/
      data fpppp(10,19) /             7.66324045d-04 /
      data fpppp(11, 1),fpppp(11, 2)/ 5.15893011d-04, 2.25192740d-04/
      data fpppp(11, 3),fpppp(11, 4)/-5.63064348d-04, 3.90127425d-04/
      data fpppp(11, 5),fpppp(11, 6)/ 2.09354487d-05,-2.96154168d-04/
      data fpppp(11, 7),fpppp(11, 8)/ 2.13851200d-04,-8.91427515d-05/
      data fpppp(11, 9),fpppp(11,10)/ 8.96265260d-05,-2.46329984d-05/
      data fpppp(11,11),fpppp(11,12)/ 8.76483358d-06, 3.25491762d-05/
      data fpppp(11,13),fpppp(11,14)/-5.18316044d-05, 9.52893651d-05/
      data fpppp(11,15),fpppp(11,16)/-9.14103541d-05, 8.30800571d-05/
      data fpppp(11,17),fpppp(11,18)/-1.53482627d-05,-9.39544037d-05/
      data fpppp(11,19) /            -1.78224024d-04 /
      data fpppp(12, 1),fpppp(12, 2)/-7.84618004d-05,-3.70191360d-05/
      data fpppp(12, 3),fpppp(12, 4)/ 9.02915408d-05,-8.29358258d-05/
      data fpppp(12, 5),fpppp(12, 6)/-1.53871356d-05, 4.83963073d-05/
      data fpppp(12, 7),fpppp(12, 8)/-4.71209862d-05, 5.10078380d-05/
      data fpppp(12, 9),fpppp(12,10)/-1.48467083d-05, 1.34198317d-05/
      data fpppp(12,11),fpppp(12,12)/ 1.39545482d-07,-2.62155901d-06/
      data fpppp(12,13),fpppp(12,14)/ 1.75895289d-05,-1.23543073d-05/
      data fpppp(12,15),fpppp(12,16)/ 2.34707942d-05,-4.07870919d-05/
      data fpppp(12,17),fpppp(12,18)/ 8.03242315d-06, 1.10695734d-05/
      data fpppp(12,19) /             2.71873993d-05 /
      data fpppp(13, 1),fpppp(13, 2)/ 1.58204589d-05, 6.52591098d-06/
      data fpppp(13, 3),fpppp(13, 4)/-1.79355046d-05, 1.06871167d-05/
      data fpppp(13, 5),fpppp(13, 6)/ 8.57331682d-07,-7.74178650d-06/
      data fpppp(13, 7),fpppp(13, 8)/ 1.37655033d-05,-9.08676868d-06/
      data fpppp(13, 9),fpppp(13,10)/ 2.80923908d-06,-2.69387397d-06/
      data fpppp(13,11),fpppp(13,12)/ 7.76081417d-07, 1.74042827d-06/
      data fpppp(13,13),fpppp(13,14)/-2.54327636d-06, 3.06186769d-06/
      data fpppp(13,15),fpppp(13,16)/-5.49122751d-06, 9.24970668d-06/
      data fpppp(13,17),fpppp(13,18)/-2.14495411d-06,-4.06471375d-06/
      data fpppp(13,19) /            -8.69158801d-06 /
      data fpppp(14, 1),fpppp(14, 2)/-3.74236916d-06,-1.50037958d-06/
      data fpppp(14, 3),fpppp(14, 4)/ 4.51549467d-06,-2.67722744d-06/
      data fpppp(14, 5),fpppp(14, 6)/-2.14017677d-07, 2.03335780d-06/
      data fpppp(14, 7),fpppp(14, 8)/-3.70703426d-06, 2.37630465d-06/
      data fpppp(14, 9),fpppp(14,10)/-8.09016146d-07, 6.90400729d-07/
      data fpppp(14,11),fpppp(14,12)/-1.74142632d-07,-4.60697421d-07/
      data fpppp(14,13),fpppp(14,14)/ 6.94958751d-07,-7.84833720d-07/
      data fpppp(14,15),fpppp(14,16)/ 1.46292845d-06,-2.36376192d-06/
      data fpppp(14,17),fpppp(14,18)/ 4.77759011d-07, 1.02510935d-06/
      data fpppp(14,19) /             2.22793699d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 9.92930672d-07, 4.13781526d-07/
      data fpppp(15, 3),fpppp(15, 4)/-1.11508367d-06, 6.80057314d-07/
      data fpppp(15, 5),fpppp(15, 6)/ 5.52915379d-08,-4.93119014d-07/
      data fpppp(15, 7),fpppp(15, 8)/ 8.62978473d-07,-5.71354638d-07/
      data fpppp(15, 9),fpppp(15,10)/ 1.76099667d-07,-1.71920856d-07/
      data fpppp(15,11),fpppp(15,12)/ 3.89826165d-08, 1.10579317d-07/
      data fpppp(15,13),fpppp(15,14)/-1.62923749d-07, 1.88709731d-07/
      data fpppp(15,15),fpppp(15,16)/-3.48091350d-07, 5.84518621d-07/
      data fpppp(15,17),fpppp(15,18)/-1.25187455d-07,-2.56479198d-07/
      data fpppp(15,19) /            -5.53032180d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-8.37254830d-07,-3.69554632d-07/
      data fpppp(16, 3),fpppp(16, 4)/ 7.63158236d-07,-4.31673244d-07/
      data fpppp(16, 5),fpppp(16, 6)/-3.36338210d-08, 2.39820588d-07/
      data fpppp(16, 7),fpppp(16, 8)/-3.85881222d-07, 3.05219894d-07/
      data fpppp(16, 9),fpppp(16,10)/-3.86424352d-08, 9.58454017d-08/
      data fpppp(16,11),fpppp(16,12)/-6.90314580d-08,-5.17711756d-08/
      data fpppp(16,13),fpppp(16,14)/ 5.61780937d-08,-1.23123938d-07/
      data fpppp(16,15),fpppp(16,16)/ 1.57870034d-07,-3.46644816d-07/
      data fpppp(16,17),fpppp(16,18)/ 1.56297105d-07, 1.63568736d-07/
      data fpppp(16,19) /             3.12460450d-07 /
 

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
      subroutine b_3(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  5.11030000d-03 /
      data f( 1, 3),f( 1, 4) /  1.74350000d-02 ,  2.54305000d-02 /
      data f( 1, 5),f( 1, 6) /  1.54521000d-02 ,  3.12600000d-03 /
      data f( 1, 7),f( 1, 8) /  1.74730000d-03 ,  7.89000000d-05 /
      data f( 1, 9),f( 1,10) / -7.21300000d-04 , -4.11640000d-03 /
      data f( 1,11),f( 1,12) / -1.28760000d-02 , -3.11310000d-03 /
      data f( 1,13),f( 1,14) / -9.92730000d-03 , -8.27330000d-02 /
      data f( 1,15),f( 1,16) / -4.00403200d-01 , -1.08207300d+00 /
      data f( 1,17),f( 1,18) / -5.97521100d-01 , -9.84485000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  6.37870000d-03 /
      data f( 2, 3),f( 2, 4) /  1.75985000d-02 ,  1.69004000d-02 /
      data f( 2, 5),f( 2, 6) /  1.74222000d-02 ,  1.62610000d-03 /
      data f( 2, 7),f( 2, 8) / -3.30600000d-04 , -8.80200000d-04 /
      data f( 2, 9),f( 2,10) / -1.16570000d-03 , -2.71450000d-03 /
      data f( 2,11),f( 2,12) / -9.55950000d-03 ,  1.66020000d-03 /
      data f( 2,13),f( 2,14) / -3.47880000d-03 , -2.66276000d-02 /
      data f( 2,15),f( 2,16) / -2.22147300d-01 , -5.29977200d-01 /
      data f( 2,17),f( 2,18) / -3.32241100d-01 , -5.91463000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  1.91081000d-02 /
      data f( 3, 3),f( 3, 4) /  1.55970000d-02 ,  1.64998000d-02 /
      data f( 3, 5),f( 3, 6) /  1.62928000d-02 , -5.83700000d-04 /
      data f( 3, 7),f( 3, 8) / -1.71010000d-03 , -6.37100000d-04 /
      data f( 3, 9),f( 3,10) / -1.14380000d-03 , -1.66380000d-03 /
      data f( 3,11),f( 3,12) / -6.89000000d-03 , -6.84660000d-03 /
      data f( 3,13),f( 3,14) / -1.02640000d-03 , -2.23210000d-02 /
      data f( 3,15),f( 3,16) / -1.21711400d-01 , -2.65255900d-01 /
      data f( 3,17),f( 3,18) / -1.91987200d-01 , -4.87971000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  1.70076000d-02 /
      data f( 4, 3),f( 4, 4) /  1.29865000d-02 ,  1.46046000d-02 /
      data f( 4, 5),f( 4, 6) /  1.02593000d-02 , -1.71070000d-03 /
      data f( 4, 7),f( 4, 8) / -2.84960000d-03 ,  6.47900000d-04 /
      data f( 4, 9),f( 4,10) / -1.04390000d-03 , -9.15100000d-04 /
      data f( 4,11),f( 4,12) / -4.82570000d-03 , -1.43102000d-02 /
      data f( 4,13),f( 4,14) /  1.30670000d-03 ,  9.02470000d-03 /
      data f( 4,15),f( 4,16) / -6.17690000d-02 , -1.39539700d-01 /
      data f( 4,17),f( 4,18) / -1.11166200d-01 , -4.09810000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  8.82930000d-03 /
      data f( 5, 3),f( 5, 4) /  1.00650000d-02 ,  1.03654000d-02 /
      data f( 5, 5),f( 5, 6) /  5.96870000d-03 , -2.69230000d-03 /
      data f( 5, 7),f( 5, 8) /  1.11330000d-03 ,  9.03000000d-05 /
      data f( 5, 9),f( 5,10) / -8.76100000d-04 , -4.04500000d-04 /
      data f( 5,11),f( 5,12) / -3.27570000d-03 , -1.03024000d-02 /
      data f( 5,13),f( 5,14) / -8.13060000d-03 , -8.15900000d-04 /
      data f( 5,15),f( 5,16) / -2.85184000d-02 , -7.23909000d-02 /
      data f( 5,17),f( 5,18) / -6.38207000d-02 , -2.72086000d-02 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -8.68000000d-04 /
      data f( 6, 3),f( 6, 4) /  5.55490000d-03 ,  7.27020000d-03 /
      data f( 6, 5),f( 6, 6) /  3.60280000d-03 ,  3.77850000d-03 /
      data f( 6, 7),f( 6, 8) /  7.28040000d-03 , -1.39700000d-04 /
      data f( 6, 9),f( 6,10) / -6.90200000d-04 , -7.07000000d-05 /
      data f( 6,11),f( 6,12) / -2.13580000d-03 , -7.25440000d-03 /
      data f( 6,13),f( 6,14) / -1.80392000d-02 , -7.67210000d-03 /
      data f( 6,15),f( 6,16) / -1.07282000d-02 , -3.66097000d-02 /
      data f( 6,17),f( 6,18) / -3.26389000d-02 , -1.44059000d-02 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 , -1.03220000d-03 /
      data f( 7, 3),f( 7, 4) /  2.84320000d-03 ,  4.46830000d-03 /
      data f( 7, 5),f( 7, 6) /  4.69360000d-03 ,  3.49310000d-03 /
      data f( 7, 7),f( 7, 8) /  5.31000000d-05 , -2.23100000d-04 /
      data f( 7, 9),f( 7,10) / -5.02400000d-04 ,  1.40800000d-04 /
      data f( 7,11),f( 7,12) / -1.31160000d-03 , -4.97170000d-03 /
      data f( 7,13),f( 7,14) / -1.28182000d-02 , -2.63628000d-02 /
      data f( 7,15),f( 7,16) / -1.19598000d-02 , -1.99921000d-02 /
      data f( 7,17),f( 7,18) / -1.15173000d-02 ,  2.96500000d-04 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 , -8.84200000d-04 /
      data f( 8, 3),f( 8, 4) /  1.21837000d-02 ,  3.57220000d-03 /
      data f( 8, 5),f( 8, 6) /  3.90060000d-03 , -3.19000000d-05 /
      data f( 8, 7),f( 8, 8) / -1.27500000d-04 , -2.43900000d-04 /
      data f( 8, 9),f( 8,10) / -3.19600000d-04 ,  2.70400000d-04 /
      data f( 8,11),f( 8,12) / -7.23500000d-04 , -3.29360000d-03 /
      data f( 8,13),f( 8,14) / -8.91070000d-03 , -1.85989000d-02 /
      data f( 8,15),f( 8,16) / -1.82623000d-02 , -2.13118000d-02 /
      data f( 8,17),f( 8,18) / -1.75080000d-02 , -1.13945000d-02 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 , -4.93100000d-04 /
      data f( 9, 3),f( 9, 4) / -3.84000000d-04 , -2.92200000d-04 /
      data f( 9, 5),f( 9, 6) / -2.64900000d-04 , -2.44300000d-04 /
      data f( 9, 7),f( 9, 8) / -2.13600000d-04 , -2.03900000d-04 /
      data f( 9, 9),f( 9,10) /  2.52000000d-05 ,  3.84200000d-04 /
      data f( 9,11),f( 9,12) / -2.35000000d-05 , -1.21080000d-03 /
      data f( 9,13),f( 9,14) / -3.87770000d-03 , -8.52520000d-03 /
      data f( 9,15),f( 9,16) / -1.42085000d-02 , -1.74861000d-02 /
      data f( 9,17),f( 9,18) / -1.46956000d-02 , -7.65420000d-03 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 , -2.73800000d-04 /
      data f(10, 3),f(10, 4) / -2.77900000d-04 , -2.70800000d-04 /
      data f(10, 5),f(10, 6) / -2.69500000d-04 , -2.44500000d-04 /
      data f(10, 7),f(10, 8) / -1.87100000d-04 , -1.32300000d-04 /
      data f(10, 9),f(10,10) /  2.73800000d-04 ,  3.94800000d-04 /
      data f(10,11),f(10,12) /  2.89400000d-04 , -1.76100000d-04 /
      data f(10,13),f(10,14) / -1.32810000d-03 , -3.37470000d-03 /
      data f(10,15),f(10,16) / -5.91250000d-03 , -7.44560000d-03 /
      data f(10,17),f(10,18) / -6.40500000d-03 , -3.50230000d-03 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -1.36600000d-04 /
      data f(11, 3),f(11, 4) / -1.75900000d-04 , -1.97900000d-04 /
      data f(11, 5),f(11, 6) / -2.05100000d-04 , -1.77200000d-04 /
      data f(11, 7),f(11, 8) / -1.21700000d-04 ,  5.36000000d-05 /
      data f(11, 9),f(11,10) /  2.00400000d-04 ,  3.34900000d-04 /
      data f(11,11),f(11,12) /  3.99100000d-04 ,  3.43900000d-04 /
      data f(11,13),f(11,14) /  5.30000000d-05 , -5.58300000d-04 /
      data f(11,15),f(11,16) / -1.35930000d-03 , -1.85840000d-03 /
      data f(11,17),f(11,18) / -1.65520000d-03 , -1.07810000d-03 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 , -6.63000000d-05 /
      data f(12, 3),f(12, 4) / -1.01200000d-04 , -1.22000000d-04 /
      data f(12, 5),f(12, 6) / -1.30000000d-04 , -1.16100000d-04 /
      data f(12, 7),f(12, 8) /  4.17000000d-05 ,  8.12000000d-05 /
      data f(12, 9),f(12,10) /  8.22000000d-05 ,  1.56600000d-04 /
      data f(12,11),f(12,12) /  3.13100000d-04 ,  4.39400000d-04 /
      data f(12,13),f(12,14) /  4.57700000d-04 ,  3.77500000d-04 /
      data f(12,15),f(12,16) /  2.19800000d-04 ,  1.27200000d-04 /
      data f(12,17),f(12,18) / -1.21300000d-04 , -2.27500000d-04 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 ,  8.66000000d-05 /
      data f(13, 3),f(13, 4) /  8.08000000d-05 ,  9.36000000d-05 /
      data f(13, 5),f(13, 6) /  7.67000000d-05 ,  4.92000000d-05 /
      data f(13, 7),f(13, 8) /  1.54000000d-05 ,  1.60000000d-06 /
      data f(13, 9),f(13,10) /  1.08000000d-05 ,  3.86000000d-05 /
      data f(13,11),f(13,12) /  8.66000000d-05 ,  1.51100000d-04 /
      data f(13,13),f(13,14) /  2.19200000d-04 ,  2.70600000d-04 /
      data f(13,15),f(13,16) /  2.82400000d-04 ,  2.41300000d-04 /
      data f(13,17),f(13,18) /  1.34900000d-04 ,  3.92000000d-05 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -2.50000000d-06 /
      data f(14, 3),f(14, 4) / -3.20000000d-06 , -2.90000000d-06 /
      data f(14, 5),f(14, 6) / -3.10000000d-06 , -5.50000000d-06 /
      data f(14, 7),f(14, 8) / -1.04000000d-05 , -8.80000000d-06 /
      data f(14, 9),f(14,10) / -6.40000000d-06 , -2.20000000d-06 /
      data f(14,11),f(14,12) /  5.20000000d-06 ,  1.40000000d-05 /
      data f(14,13),f(14,14) /  2.21000000d-05 ,  2.81000000d-05 /
      data f(14,15),f(14,16) /  3.13000000d-05 ,  3.14000000d-05 /
      data f(14,17),f(14,18) /  2.74000000d-05 ,  1.35000000d-05 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -2.10000000d-06 /
      data f(15, 3),f(15, 4) / -3.30000000d-06 , -4.00000000d-06 /
      data f(15, 5),f(15, 6) / -4.30000000d-06 , -5.00000000d-06 /
      data f(15, 7),f(15, 8) / -5.70000000d-06 , -4.90000000d-06 /
      data f(15, 9),f(15,10) / -3.40000000d-06 , -3.50000000d-06 /
      data f(15,11),f(15,12) / -1.40000000d-06 ,  7.00000000d-07 /
      data f(15,13),f(15,14) /  2.30000000d-06 ,  3.20000000d-06 /
      data f(15,15),f(15,16) /  2.80000000d-06 ,  1.70000000d-06 /
      data f(15,17),f(15,18) /  1.30000000d-06 ,  2.00000000d-07 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 , -1.00000000d-07 /
      data f(16, 3),f(16, 4) / -3.00000000d-07 , -4.00000000d-07 /
      data f(16, 5),f(16, 6) / -6.00000000d-07 , -8.00000000d-07 /
      data f(16, 7),f(16, 8) / -6.00000000d-07 , -5.00000000d-07 /
      data f(16, 9),f(16,10) / -8.00000000d-07 , -2.00000000d-07 /
      data f(16,11),f(16,12) /  1.30000000d-06 ,  3.30000000d-06 /
      data f(16,13),f(16,14) /  4.70000000d-06 ,  5.80000000d-06 /
      data f(16,15),f(16,16) /  6.50000000d-06 ,  5.80000000d-06 /
      data f(16,17),f(16,18) /  3.50000000d-06 ,  6.00000000d-07 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 1.91548253d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 9.05385133d-01, 6.42074947d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-9.21423721d-02,-1.55142315d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 4.36352017d-01,-2.61902569d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-4.31922220d-02,-1.53094930d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-5.55142292d-02, 1.82278541d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 2.30625888d-02,-5.69606704d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 3.84656938d-02, 2.81821409d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.09597078d-02,-3.67589307d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/-1.00256819d-02,-1.69172569d-04/
      data fpp( 1,11,1),fpp( 1,11,2)/-1.72666941d-02, 3.58496167d-04/
      data fpp( 1,12,1),fpp( 1,12,2)/-6.93335840d-01,-1.53462101d-04/
      data fpp( 1,13,1),fpp( 1,13,2)/-1.75403648d-01,-7.39273762d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/-3.07648652d+00,-8.48932851d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/-2.85355073d+00,-1.05568648d-02/
      data fpp( 1,16,1),fpp( 1,16,2)/-1.08372028d+01, 2.12364162d-02/
      data fpp( 1,17,1),fpp( 1,17,2)/-4.72672928d+00,-4.41549794d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.37175774d+00,-2.70318245d-03/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00,-8.80921828d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 1.98420005d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 3.63349734d-01, 8.35729901d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-5.58902557d-02,-2.42245965d-04/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.18213465d-01, 1.70336871d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-5.58480560d-02,-3.65907520d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-3.18490417d-02, 3.14219209d-04/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.91748223d-02,-6.06053167d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 2.12486123d-02, 1.26280576d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.24680843d-02, 2.59390863d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-8.74863620d-03,-1.92182403d-04/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.60816118d-02, 4.25018524d-04/
      data fpp( 2,12,1),fpp( 2,12,2)/-3.25500820d-01,-4.24009695d-04/
      data fpp( 2,13,1),fpp( 2,13,2)/-1.42740203d-01, 2.89498257d-04/
      data fpp( 2,14,1),fpp( 2,14,2)/-1.67383197d+00,-1.81457133d-03/
      data fpp( 2,15,1),fpp( 2,15,2)/-1.99572353d+00,-3.37346693d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/-7.29715187d+00, 8.56982706d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/-3.20315893d+00,-5.71881325d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-7.48954519d-01,-1.76077976d-03/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-5.22190962d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-4.85953034d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-6.39634069d-01,-2.47331932d-04/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-9.04660500d-03, 1.18128761d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-8.97808777d-02, 3.96508867d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.98340554d-01,-3.43320308d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 7.64253958d-02, 3.33460346d-04/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-3.50018781d-02,-4.55150762d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 5.68698569d-02,-1.94360413d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-8.87045074d-04, 2.84772415d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-7.65977330d-03,-9.52709245d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.54568588d-02, 7.02344566d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 3.32412115d-03, 1.30509098d-04/
      data fpp( 3,13,1),fpp( 3,13,2)/ 1.46949462d-01,-2.45662850d-04/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.00199439d+00,-7.74745700d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/-8.36555141d-01,-1.34110235d-03/
      data fpp( 3,16,1),fpp( 3,16,2)/-3.08036470d+00, 3.48990910d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.21454998d+00, 3.90257933d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/ 2.46258176d-02,-8.55656838d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00,-2.63121058d-03/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-4.60774196d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-2.92984575d-02,-2.42670607d-04/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 7.26675730d-04, 1.69734625d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-8.32799542d-02,-9.79158919d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.13595272d-01,-1.35875057d-04/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-1.11432542d-01, 1.83934121d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 1.46832690d-01, 5.00045729d-05/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-9.24430398d-02,-1.05768413d-04/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.78009599d-03, 6.17110778d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/-5.91227058d-03,-3.18398985d-05/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.28709529d-02,-1.76715484d-04/
      data fpp( 4,12,1),fpp( 4,12,2)/ 4.68684336d-01, 4.04267834d-04/
      data fpp( 4,13,1),fpp( 4,13,2)/-4.62952645d-01, 6.57281493d-05/
      data fpp( 4,14,1),fpp( 4,14,2)/-2.27828061d+00,-1.14111443d-03/
      data fpp( 4,15,1),fpp( 4,15,2)/-7.32080904d-01,-2.11972426d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/-1.23215435d+00, 1.57038413d-03/
      data fpp( 4,17,1),fpp( 4,17,2)/-8.53576143d-01, 2.99087887d-04/
      data fpp( 4,18,1),fpp( 4,18,2)/ 2.70486249d-01,-2.58033682d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-1.01920516d-03/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-1.38538630d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-1.54842101d-01,-8.38967404d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-4.05100979d-02, 1.85095915d-05/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 7.13006943d-02,-4.62596254d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 5.39446703d-03,-1.15297090d-04/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 3.91114771d-01, 2.51589984d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 2.13031117d-01,-1.43066847d-04/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 3.65123023d-02, 3.09614036d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-4.83388787d-05, 2.26172324d-05/
      data fpp( 5,10,1),fpp( 5,10,2)/-4.40614436d-03,-3.51503331d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.02043297d-02,-8.25839000d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.57351464d-01, 1.16155933d-04/
      data fpp( 5,13,1),fpp( 5,13,2)/-6.06988804d-02, 1.69870168d-04/
      data fpp( 5,14,1),fpp( 5,14,2)/ 9.33183043d-01,-4.87062604d-04/
      data fpp( 5,15,1),fpp( 5,15,2)/-2.38891244d-01,-3.22651753d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/-7.76127919d-01, 8.07469614d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/-3.92470448d-01, 2.39335296d-04/
      data fpp( 5,18,1),fpp( 5,18,2)/-2.13125813d-01,-8.22967989d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-4.74358101d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 1.89398494d-04/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 4.20816862d-01, 7.99000113d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-7.69762841d-02,-7.15445394d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-3.03228232d-02,-7.61778536d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.53531860d-01, 5.32939539d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-3.35166542d-01, 9.35880379d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-6.68327159d-01,-2.28074106d-04/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-4.46616946d-03, 1.63388385d-04/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 1.28259527d-04,-1.33034330d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/-2.98315198d-03,-3.99746528d-05/
      data fpp( 6,11,1),fpp( 6,11,2)/-7.82672829d-03, 1.21260441d-05/
      data fpp( 6,12,1),fpp( 6,12,2)/ 1.67515208d-02,-1.91739524d-04/
      data fpp( 6,13,1),fpp( 6,13,2)/ 6.35053167d-01, 4.14860051d-04/
      data fpp( 6,14,1),fpp( 6,14,2)/-1.00679156d+00,-1.98586679d-04/
      data fpp( 6,15,1),fpp( 6,15,2)/-6.31414121d-01,-4.25905335d-04/
      data fpp( 6,16,1),fpp( 6,16,2)/-3.68473978d-01, 5.32684018d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/-1.09706580d-03, 8.63072645d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 4.36562003d-01,-2.21810756d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-2.27208962d-04/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 1.17332385d-04/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-9.84603453d-02, 5.57212300d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 6.18175234d-01,-4.57613051d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 9.39855986d-02,-7.69400980d-06/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-1.01016908d-01,-7.45065574d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-6.38786036d-02,-4.80513672d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 4.51117519d-01, 6.52861246d-05/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 3.34237552d-03,-2.32651313d-05/
      data fpp( 7, 9,1),fpp( 7, 9,2)/-1.79699229d-04, 2.75884005d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/-2.00624771d-03,-3.17384707d-05/
      data fpp( 7,11,1),fpp( 7,11,2)/-5.84375713d-03,-2.63705178d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/-2.44496190d-02, 4.75854183d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/-2.10073787d-01,-2.43847650d-04/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.31880821d+00, 6.28746056d-04/
      data fpp( 7,15,1),fpp( 7,15,2)/-8.87222713d-02,-5.94280576d-04/
      data fpp( 7,16,1),fpp( 7,16,2)/-6.24516170d-01, 4.02258247d-04/
      data fpp( 7,17,1),fpp( 7,17,2)/-1.11217129d+00,-2.43264127d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/-1.24816720d+00,-1.04612596d-04/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-2.83841202d-04/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 4.62769968d-04/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.98545195d-02, 2.05653063d-04/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-5.87894653d-01,-4.48256221d-04/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-5.97495712d-02, 2.86607820d-04/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-3.20342295d-02,-1.61781059d-04/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 1.04740956d-01, 1.04862417d-04/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-7.91379163d-02,-2.74546096d-05/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 4.86667390d-04, 3.70802108d-06/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.59462613d-04, 1.50645253d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/-1.27685717d-03,-2.40241221d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/-4.21324319d-03,-1.40020369d-05/
      data fpp( 8,12,1),fpp( 8,12,2)/-9.64304474d-03,-1.45397305d-05/
      data fpp( 8,13,1),fpp( 8,13,2)/ 8.21698172d-03,-1.10659041d-04/
      data fpp( 8,14,1),fpp( 8,14,2)/-3.00251276d-01, 2.12909896d-04/
      data fpp( 8,15,1),fpp( 8,15,2)/ 2.25668206d-01,-1.39492541d-04/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.75943658d-01, 1.41894269d-04/
      data fpp( 8,17,1),fpp( 8,17,2)/ 3.82937222d-01,-1.68865339d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/ 5.97096791d-01, 6.42338668d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 7.68110666d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 1.19916700d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-6.76713589d-03, 6.47265994d-06/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 2.82770091d-01,-1.75030978d-06/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 5.45484142d-02,-5.09420817d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 4.98798924d-02,-8.20069520d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-2.58735668d-02, 4.35448625d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.21712396d-02,-1.05378755d-06/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-7.11899282d-05, 2.51970157d-06/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-2.11762547d-04, 4.13898128d-06/
      data fpp( 9,10,1),fpp( 9,10,2)/-6.18804641d-04,-1.12816267d-05/
      data fpp( 9,11,1),fpp( 9,11,2)/-2.29589187d-03,-5.01447453d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/-6.59855628d-03,-1.54364752d-05/
      data fpp( 9,13,1),fpp( 9,13,2)/-2.39390516d-02,-2.20156248d-05/
      data fpp( 9,14,1),fpp( 9,14,2)/ 3.68209743d-02,-1.53370258d-05/
      data fpp( 9,15,1),fpp( 9,15,2)/-7.93848341d-03, 2.12157279d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/ 2.68683622d-02, 7.48161142d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/-3.79585205d-02, 4.36058154d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/-1.50120524d-01, 5.81462417d-06/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00,-3.00963121d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 5.18581867d-06/
      data fpp(10, 2,1),fpp(10, 2,2)/ 7.71524042d-04, 2.88936266d-06/
      data fpp(10, 3,1),fpp(10, 3,2)/-6.79182106d-02,-5.61269323d-07/
      data fpp(10, 4,1),fpp(10, 4,2)/-1.27265855d-02, 2.77146267d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.14515902d-02, 1.02410816d-07/
      data fpp(10, 6,1),fpp(10, 6,2)/ 6.71081116d-03, 9.84642110d-07/
      data fpp(10, 7,1),fpp(10, 7,2)/-5.32454196d-03,-2.09697926d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/ 9.83092323d-04, 7.24727491d-06/
      data fpp(10, 9,1),fpp(10, 9,2)/-2.60098720d-03,-5.81412039d-06/
      data fpp(10,10,1),fpp(10,10,2)/-1.17924268d-04,-1.09679336d-06/
      data fpp(10,11,1),fpp(10,11,2)/-1.11943935d-03,-3.38270617d-06/
      data fpp(10,12,1),fpp(10,12,2)/-3.26648015d-03,-6.97838197d-06/
      data fpp(10,13,1),fpp(10,13,2)/-5.58827533d-03,-9.89376596d-06/
      data fpp(10,14,1),fpp(10,14,2)/-3.16526209d-02,-7.12255420d-06/
      data fpp(10,15,1),fpp(10,15,2)/-3.48317727d-02, 8.91198276d-06/
      data fpp(10,16,1),fpp(10,16,2)/-5.03621064d-02, 3.17566232d-05/
      data fpp(10,17,1),fpp(10,17,2)/-2.56706395d-02, 1.84835245d-05/
      data fpp(10,18,1),fpp(10,18,2)/ 1.88203034d-02, 6.03527871d-06/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00,-6.64863935d-06/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 1.74565015d-06/
      data fpp(11, 2,1),fpp(11, 2,2)/-6.49977842d-04, 1.02769970d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 1.75544854d-02,-1.84489411d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/ 3.28457642d-03, 8.40960668d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/ 3.00541093d-03, 5.70064674d-07/
      data fpp(11, 6,1),fpp(11, 6,2)/-1.83886669d-03,-2.58354762d-07/
      data fpp(11, 7,1),fpp(11, 7,2)/ 2.20595939d-03, 2.11935437d-06/
      data fpp(11, 8,1),fpp(11, 8,2)/-1.16858042d-03,-1.03106273d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/ 3.13363949d-04, 2.94896542d-07/
      data fpp(11,10,1),fpp(11,10,2)/-8.36028922d-04,-8.86523439d-07/
      data fpp(11,11,1),fpp(11,11,2)/-8.87504860d-04,-9.66802786d-07/
      data fpp(11,12,1),fpp(11,12,2)/-1.52282643d-03,-2.41026542d-06/
      data fpp(11,13,1),fpp(11,13,2)/-4.07256752d-03,-3.53413555d-06/
      data fpp(11,14,1),fpp(11,14,2)/-2.42874438d-03,-2.67719240d-06/
      data fpp(11,15,1),fpp(11,15,2)/-7.85803150d-03, 2.86090515d-06/
      data fpp(11,16,1),fpp(11,16,2)/-7.31330661d-03, 9.34757182d-06/
      data fpp(11,17,1),fpp(11,17,2)/-1.19416813d-02, 1.88680759d-06/
      data fpp(11,18,1),fpp(11,18,2)/-1.40328735d-02, 5.53919783d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 6.01640108d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 4.74397398d-07/
      data fpp(12, 2,1),fpp(12, 2,2)/ 2.22787326d-04, 3.39205205d-07/
      data fpp(12, 3,1),fpp(12, 3,2)/-2.95493100d-03, 5.27817830d-08/
      data fpp(12, 4,1),fpp(12, 4,2)/-3.39720199d-04, 2.95667663d-07/
      data fpp(12, 5,1),fpp(12, 5,2)/-3.13253476d-04,-4.67452435d-07/
      data fpp(12, 6,1),fpp(12, 6,2)/ 4.95855616d-04, 2.88814208d-06/
      data fpp(12, 7,1),fpp(12, 7,2)/-1.14729559d-03,-2.45111588d-06/
      data fpp(12, 8,1),fpp(12, 8,2)/-1.07970646d-04,-1.81678566d-07/
      data fpp(12, 9,1),fpp(12, 9,2)/ 2.72331403d-04, 8.67830143d-07/
      data fpp(12,10,1),fpp(12,10,2)/ 6.20439956d-04, 1.11435799d-06/
      data fpp(12,11,1),fpp(12,11,2)/-2.73412152d-05,-3.99262116d-07/
      data fpp(12,12,1),fpp(12,12,2)/-8.30214124d-04,-1.32930953d-06/
      data fpp(12,13,1),fpp(12,13,2)/-1.55505457d-03,-7.63499767d-07/
      data fpp(12,14,1),fpp(12,14,2)/-3.76680162d-03,-1.52669140d-06/
      data fpp(12,15,1),fpp(12,15,2)/-5.11450130d-03, 2.22026538d-06/
      data fpp(12,16,1),fpp(12,16,2)/-6.82306714d-03,-3.44837010d-06/
      data fpp(12,17,1),fpp(12,17,2)/-3.74423518d-03, 2.21921503d-06/
      data fpp(12,18,1),fpp(12,18,2)/-4.55209399d-04, 3.10950999d-06/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 5.36474500d-06/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00,-1.95982063d-06/
      data fpp(13, 2,1),fpp(13, 2,2)/-2.69573058d-04,-1.07235874d-06/
      data fpp(13, 3,1),fpp(13, 3,2)/ 2.83150292d-04, 7.05255576d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-2.40327613d-04,-6.32663566d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/-2.23945038d-04, 4.33986903d-08/
      data fpp(13, 6,1),fpp(13, 6,2)/-3.09533500d-04,-1.76931195d-07/
      data fpp(13, 7,1),fpp(13, 7,2)/ 2.20307075d-04, 2.86326089d-07/
      data fpp(13, 8,1),fpp(13, 8,2)/ 9.94021466d-05, 2.31626841d-07/
      data fpp(13, 9,1),fpp(13, 9,2)/ 1.63238169d-05, 1.67166549d-07/
      data fpp(13,10,1),fpp(13,10,2)/-1.17054084d-05, 2.15706963d-07/
      data fpp(13,11,1),fpp(13,11,2)/ 1.98776075d-04, 1.82005597d-07/
      data fpp(13,12,1),fpp(13,12,2)/ 3.76255587d-04, 4.62706478d-08/
      data fpp(13,13,1),fpp(13,13,2)/ 4.14047477d-04,-1.51088188d-07/
      data fpp(13,14,1),fpp(13,14,2)/ 6.43777047d-04,-4.43917895d-07/
      data fpp(13,15,1),fpp(13,15,2)/ 6.98919641d-04,-4.49240233d-07/
      data fpp(13,16,1),fpp(13,16,2)/ 9.83254717d-04,-9.33121175d-07/
      data fpp(13,17,1),fpp(13,17,2)/ 3.33946196d-04, 2.63724932d-07/
      data fpp(13,18,1),fpp(13,18,2)/-2.24935049d-04, 5.20221448d-07/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 1.04538928d-06/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00, 2.63355043d-08/
      data fpp(14, 2,1),fpp(14, 2,2)/ 1.04975511d-04, 1.73289913d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-4.39853790d-05, 1.23485304d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/ 9.92929373d-05,-6.72311293d-09/
      data fpp(14, 5,1),fpp(14, 5,2)/ 8.86618516d-05,-1.54560787d-08/
      data fpp(14, 6,1),fpp(14, 6,2)/ 1.02722693d-04,-6.34525724d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/-4.70734307d-05, 1.19266368d-07/
      data fpp(14, 8,1),fpp(14, 8,2)/-2.10211171d-05,-2.36129002d-08/
      data fpp(14, 9,1),fpp(14, 9,2)/ 3.26284773d-06, 2.31852327d-08/
      data fpp(14,10,1),fpp(14,10,2)/ 1.76962469d-05, 3.88719694d-08/
      data fpp(14,11,1),fpp(14,11,2)/-2.52576182d-05, 1.33268896d-08/
      data fpp(14,12,1),fpp(14,12,2)/-5.44096983d-05,-8.17952796d-09/
      data fpp(14,13,1),fpp(14,13,2)/-4.47651445d-05,-2.26087778d-08/
      data fpp(14,14,1),fpp(14,14,2)/-9.09803317d-05,-2.73853608d-08/
      data fpp(14,15,1),fpp(14,15,2)/-1.03958275d-04,-3.58497790d-08/
      data fpp(14,16,1),fpp(14,16,2)/-1.95380583d-04,-1.52155232d-08/
      data fpp(14,17,1),fpp(14,17,2)/-5.95709984d-05,-1.49288128d-07/
      data fpp(14,18,1),fpp(14,18,2)/ 6.37598455d-05, 1.83680367d-08/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00, 9.98159817d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00, 1.27738098d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/-1.60789858d-05, 9.45238034d-09/
      data fpp(15, 3,1),fpp(15, 3,2)/ 1.86412236d-05, 3.41666879d-09/
      data fpp(15, 4,1),fpp(15, 4,2)/-1.37441365d-05, 6.88094448d-09/
      data fpp(15, 5,1),fpp(15, 5,2)/-1.28023684d-05,-6.94044671d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/-1.85572723d-05,-3.11915765d-09/
      data fpp(15, 7,1),fpp(15, 7,2)/ 1.37366476d-05, 1.94170773d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/ 6.13232166d-06, 1.54508484d-08/
      data fpp(15, 9,1),fpp(15, 9,2)/ 9.24792130d-07,-3.92204710d-08/
      data fpp(15,10,1),fpp(15,10,2)/ 1.70420906d-07, 4.54310355d-08/
      data fpp(15,11,1),fpp(15,11,2)/ 1.44543976d-05,-1.05036709d-08/
      data fpp(15,12,1),fpp(15,12,2)/ 2.70832063d-05,-3.41635177d-09/
      data fpp(15,13,1),fpp(15,13,2)/ 3.09631012d-05,-5.83092201d-09/
      data fpp(15,14,1),fpp(15,14,2)/ 4.65442797d-05,-1.52599602d-08/
      data fpp(15,15,1),fpp(15,15,2)/ 5.08134597d-05,-1.11292372d-08/
      data fpp(15,16,1),fpp(15,16,2)/ 6.85676163d-05, 1.77769089d-08/
      data fpp(15,17,1),fpp(15,17,2)/ 2.64377979d-05,-1.79783986d-08/
      data fpp(15,18,1),fpp(15,18,2)/-1.15043334d-05, 1.21366853d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00, 2.34316573d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00,-2.81629318d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/-1.54537214d-05,-1.36741364d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/-3.07138261d-05, 2.28594774d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/-1.68815032d-05,-1.77637734d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/-1.27666729d-05,-1.18043839d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/-4.32422097d-06, 6.49813091d-09/
      data fpp(16, 7,1),fpp(16, 7,2)/-1.57065381d-05,-8.12085233d-10/
      data fpp(16, 8,1),fpp(16, 8,2)/-7.39366083d-06,-9.24978997d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/-6.52453892d-06, 1.38112451d-08/
      data fpp(16,10,1),fpp(16,10,2)/-8.86556760d-06, 8.00480947d-09/
      data fpp(16,11,1),fpp(16,11,2)/-2.29429131d-05, 8.16951699d-09/
      data fpp(16,12,1),fpp(16,12,2)/-3.89708889d-05,-1.06828774d-08/
      data fpp(16,13,1),fpp(16,13,2)/-5.19669078d-05,-1.43800726d-09/
      data fpp(16,14,1),fpp(16,14,2)/-6.78607113d-05,-1.56509354d-09/
      data fpp(16,15,1),fpp(16,15,2)/-6.91060156d-05,-1.63016186d-08/
      data fpp(16,16,1),fpp(16,16,2)/-6.58716653d-05,-1.72284321d-08/
      data fpp(16,17,1),fpp(16,17,2)/-2.08453275d-05,-1.07846529d-08/
      data fpp(16,18,1),fpp(16,18,2)/ 9.40788099d-06, 2.43670437d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00, 5.13164782d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-5.03450003d-02,-2.49760691d-02/
      data fpppp( 1, 3),fpppp( 1, 4)/ 3.60745186d-02,-2.77606915d-02/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.44859299d-02,-2.14969393d-03/
      data fpppp( 1, 7),fpppp( 1, 8)/-4.33224639d-04, 9.21697007d-05/
      data fpppp( 1, 9),fpppp( 1,10)/-1.90999962d-03, 6.73906457d-03/
      data fpppp( 1,11),fpppp( 1,12)/-2.36215960d-02, 4.76176314d-02/
      data fpppp( 1,13),fpppp( 1,14)/-9.52088494d-02, 1.28076863d-01/
      data fpppp( 1,15),fpppp( 1,16)/-2.29657483d-01, 2.98157797d-01/
      data fpppp( 1,17),fpppp( 1,18)/-1.17326170d-01, 5.81676190d-03/
      data fpppp( 1,19) /            -2.49337061d-02 /
      data fpppp( 2, 1),fpppp( 2, 2)/-2.11931921d-02,-1.06099762d-02/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.66777135d-02,-1.45002553d-02/
      data fpppp( 2, 5),fpppp( 2, 6)/ 8.43339323d-03,-1.34968549d-03/
      data fpppp( 2, 7),fpppp( 2, 8)/-8.13160304d-04, 4.65322263d-04/
      data fpppp( 2, 9),fpppp( 2,10)/-1.09938783d-03, 3.18605751d-03/
      data fpppp( 2,11),fpppp( 2,12)/-1.08118175d-02, 2.19360385d-02/
      data fpppp( 2,13),fpppp( 2,14)/-4.74015470d-02, 6.48390064d-02/
      data fpppp( 2,15),fpppp( 2,16)/-1.39402466d-01, 1.93998652d-01/
      data fpppp( 2,17),fpppp( 2,18)/-7.28668646d-02,-9.18504842d-04/
      data fpppp( 2,19) /            -2.57741097d-02 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 3.14310500d-02, 1.48754114d-02/
      data fpppp( 3, 3),fpppp( 3, 4)/-1.47194037d-02, 1.32289910d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/ 7.75828308d-03,-9.35649384d-03/
      data fpppp( 3, 7),fpppp( 3, 8)/ 6.49609887d-03,-4.42996110d-03/
      data fpppp( 3, 9),fpppp( 3,10)/ 2.24602732d-03,-1.49509775d-03/
      data fpppp( 3,11),fpppp( 3,12)/ 3.67290225d-03,-1.16018273d-02/
      data fpppp( 3,13),fpppp( 3,14)/ 5.02250686d-02,-8.66132717d-02/
      data fpppp( 3,15),fpppp( 3,16)/ 1.46123501d-02, 6.38482702d-02/
      data fpppp( 3,17),fpppp( 3,18)/-2.34279748d-02,-7.73470598d-03/
      data fpppp( 3,19) /            -2.14612983d-02 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 1.90177364d-03, 1.44326752d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/-4.11542828d-03, 8.17653983d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.17378197d-02, 1.34605566d-02/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.31068239d-02, 9.11428128d-03/
      data fpppp( 4, 9),fpppp( 4,10)/-3.28036927d-03,-2.22773433d-03/
      data fpppp( 4,11),fpppp( 4,12)/ 1.22953276d-02,-1.76427380d-02/
      data fpppp( 4,13),fpppp( 4,14)/-2.65159118d-02, 7.06849262d-02/
      data fpppp( 4,15),fpppp( 4,16)/-5.45321327d-02, 2.46672159d-02/
      data fpppp( 4,17),fpppp( 4,18)/ 8.58236804d-03,-1.42676367d-02/
      data fpppp( 4,19) /            -3.51847395d-02 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 5.52036482d-03, 2.46839978d-03/
      data fpppp( 5, 3),fpppp( 5, 4)/ 7.56482310d-04,-5.64560167d-03/
      data fpppp( 5, 5),fpppp( 5, 6)/ 1.11629032d-02,-1.19084193d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/ 2.64253636d-03, 1.43216414d-03/
      data fpppp( 5, 9),fpppp( 5,10)/ 2.62975056d-05, 3.94815982d-04/
      data fpppp( 5,11),fpppp( 5,12)/-1.69198422d-03,-2.10781603d-03/
      data fpppp( 5,13),fpppp( 5,14)/ 2.47512314d-02,-4.30633493d-02/
      data fpppp( 5,15),fpppp( 5,16)/ 1.75447931d-02, 1.09744337d-02/
      data fpppp( 5,17),fpppp( 5,18)/-6.18887894d-03, 1.52231195d-03/
      data fpppp( 5,19) /             2.12650184d-03 /
      data fpppp( 6, 1),fpppp( 6, 2)/-2.31279993d-02,-1.05634339d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.02651344d-02, 2.16969266d-03/
      data fpppp( 6, 5),fpppp( 6, 6)/-1.07118317d-02, 3.24449043d-04/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.87463026d-02,-1.54883631d-02/
      data fpppp( 6, 9),fpppp( 6,10)/ 3.65115627d-03, 4.21387636d-04/
      data fpppp( 6,11),fpppp( 6,12)/-5.44063670d-03, 2.31064687d-02/
      data fpppp( 6,13),fpppp( 6,14)/-5.13618342d-02, 4.67320856d-02/
      data fpppp( 6,15),fpppp( 6,16)/-1.45331779d-02, 4.65438811d-03/
      data fpppp( 6,17),fpppp( 6,18)/ 2.18183157d-03,-9.16478500d-03/
      data fpppp( 6,19) /            -1.79759559d-02 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 2.69911474d-02, 1.15890058d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/-2.44414151d-02, 1.17271418d-02/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.71592436d-03, 1.30650042d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.08726235d-02, 1.26592138d-02/
      data fpppp( 7, 9),fpppp( 7,10)/-3.10904775d-03,-1.21291255d-04/
      data fpppp( 7,11),fpppp( 7,12)/ 3.47355512d-03,-1.46590304d-02/
      data fpppp( 7,13),fpppp( 7,14)/ 4.51414680d-02,-6.30364716d-02/
      data fpppp( 7,15),fpppp( 7,16)/ 3.08196697d-02,-7.93801218d-03/
      data fpppp( 7,17),fpppp( 7,18)/ 3.82070583d-03, 1.37547414d-02/
      data fpppp( 7,19) /             2.42101149d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.22110362d-02,-9.67599722d-03/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.32588036d-02,-1.52055620d-02/
      data fpppp( 8, 5),fpppp( 8, 6)/ 7.53766006d-03,-8.40148757d-03/
      data fpppp( 8, 7),fpppp( 8, 8)/ 6.82904674d-03,-3.10449200d-03/
      data fpppp( 8, 9),fpppp( 8,10)/ 7.72678444d-04,-1.44976464d-05/
      data fpppp( 8,11),fpppp( 8,12)/-8.23827346d-04, 3.16020210d-03/
      data fpppp( 8,13),fpppp( 8,14)/-1.04195914d-02, 1.89384663d-02/
      data fpppp( 8,15),fpppp( 8,16)/-1.52710094d-02, 7.60692941d-03/
      data fpppp( 8,17),fpppp( 8,18)/ 2.46378515d-04,-8.16248316d-03/
      data fpppp( 8,19) /            -1.62718275d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 1.03255300d-02, 4.51933617d-03/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.06246129d-02, 6.91358138d-03/
      data fpppp( 9, 5),fpppp( 9, 6)/-3.61652327d-03, 3.28741544d-03/
      data fpppp( 9, 7),fpppp( 9, 8)/-2.10524257d-03, 9.16320673d-04/
      data fpppp( 9, 9),fpppp( 9,10)/-2.33928711d-04, 3.40600148d-06/
      data fpppp( 9,11),fpppp( 9,12)/ 1.44101997d-04,-7.37348620d-04/
      data fpppp( 9,13),fpppp( 9,14)/ 2.02302263d-03,-2.66871062d-03/
      data fpppp( 9,15),fpppp( 9,16)/ 2.32065085d-03,-1.83991457d-03/
      data fpppp( 9,17),fpppp( 9,18)/-9.39016269d-04, 2.75587243d-03/
      data fpppp( 9,19) /             5.65247816d-03 /
      data fpppp(10, 1),fpppp(10, 2)/-2.44211707d-03,-1.06645598d-03/
      data fpppp(10, 3),fpppp(10, 4)/ 2.54026548d-03,-1.66172437d-03/
      data fpppp(10, 5),fpppp(10, 6)/ 8.71634218d-04,-8.11568127d-04/
      data fpppp(10, 7),fpppp(10, 8)/ 5.62773020d-04,-3.38944708d-04/
      data fpppp(10, 9),fpppp(10,10)/ 1.99502986d-04,-9.50386878d-05/
      data fpppp(10,11),fpppp(10,12)/-2.84229149d-05, 1.39998804d-04/
      data fpppp(10,13),fpppp(10,14)/-5.42057563d-04, 6.03678426d-04/
      data fpppp(10,15),fpppp(10,16)/-4.99544522d-04, 6.53428752d-04/
      data fpppp(10,17),fpppp(10,18)/ 2.99137554d-04,-6.62010405d-04/
      data fpppp(10,19) /            -1.44977072d-03 /
      data fpppp(11, 1),fpppp(11, 2)/ 6.52808944d-04, 2.86591610d-04/
      data fpppp(11, 3),fpppp(11, 4)/-6.67908919d-04, 4.36581734d-04/
      data fpppp(11, 5),fpppp(11, 6)/-2.38973409d-04, 2.45405172d-04/
      data fpppp(11, 7),fpppp(11, 8)/-2.09301058d-04, 1.46637108d-04/
      data fpppp(11, 9),fpppp(11,10)/-8.58583240d-05, 3.89159535d-05/
      data fpppp(11,11),fpppp(11,12)/-3.93047383d-06,-5.82247962d-05/
      data fpppp(11,13),fpppp(11,14)/ 1.21964487d-04,-1.78019299d-04/
      data fpppp(11,15),fpppp(11,16)/ 1.65726093d-04,-1.26444353d-04/
      data fpppp(11,17),fpppp(11,18)/ 2.96653418d-05, 1.60013938d-04/
      data fpppp(11,19) /             2.97722848d-04 /
      data fpppp(12, 1),fpppp(12, 2)/-1.17064726d-04,-5.17544133d-05/
      data fpppp(12, 3),fpppp(12, 4)/ 1.20052040d-04,-8.08779995d-05/
      data fpppp(12, 5),fpppp(12, 6)/ 4.81353137d-05,-6.47047131d-05/
      data fpppp(12, 7),fpppp(12, 8)/ 6.35479208d-05,-2.85384012d-05/
      data fpppp(12, 9),fpppp(12,10)/ 1.10643103d-05,-1.76504497d-05/
      data fpppp(12,11),fpppp(12,12)/-2.15894896d-07, 9.20852510d-06/
      data fpppp(12,13),fpppp(12,14)/-3.19362579d-05, 2.93221104d-05/
      data fpppp(12,15),fpppp(12,16)/-3.35093416d-05, 8.30632864d-05/
      data fpppp(12,17),fpppp(12,18)/-1.14999358d-05,-2.44519141d-05/
      data fpppp(12,19) /            -6.07213906d-05 /
      data fpppp(13, 1),fpppp(13, 2)/ 2.54394302d-05, 1.17599852d-05/
      data fpppp(13, 3),fpppp(13, 4)/-2.31415864d-05, 1.62342849d-05/
      data fpppp(13, 5),fpppp(13, 6)/-9.40392462d-06, 1.52631513d-05/
      data fpppp(13, 7),fpppp(13, 8)/-1.47229384d-05, 4.58387201d-06/
      data fpppp(13, 9),fpppp(13,10)/-1.34295371d-06, 4.09088910d-06/
      data fpppp(13,11),fpppp(13,12)/-7.09960159d-07,-3.23116680d-06/
      data fpppp(13,13),fpppp(13,14)/ 5.25337006d-06,-6.26605263d-06/
      data fpppp(13,15),fpppp(13,16)/ 9.33562191d-06,-1.73248861d-05/
      data fpppp(13,17),fpppp(13,18)/ 3.94530666d-06, 6.96929611d-06/
      data fpppp(13,19) /             1.52064865d-05 /
      data fpppp(14, 1),fpppp(14, 2)/-7.50403936d-06,-3.53352546d-06/
      data fpppp(14, 3),fpppp(14, 4)/ 6.40195713d-06,-4.53995069d-06/
      data fpppp(14, 5),fpppp(14, 6)/ 2.52328150d-06,-4.07165968d-06/
      data fpppp(14, 7),fpppp(14, 8)/ 3.93193929d-06,-1.10519125d-06/
      data fpppp(14, 9),fpppp(14,10)/ 3.82724788d-07,-1.01674184d-06/
      data fpppp(14,11),fpppp(14,12)/ 2.41006725d-07, 8.80822043d-07/
      data fpppp(14,13),fpppp(14,14)/-1.43649687d-06, 1.51358098d-06/
      data fpppp(14,15),fpppp(14,16)/-2.62359243d-06, 4.27412685d-06/
      data fpppp(14,17),fpppp(14,18)/-8.39001390d-07,-1.66684576d-06/
      data fpppp(14,19) /            -3.71905695d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 1.57715050d-06, 7.27770148d-07/
      data fpppp(15, 3),fpppp(15, 4)/-1.44027938d-06, 1.00701321d-06/
      data fpppp(15, 5),fpppp(15, 6)/-5.88145780d-07, 9.43769588d-07/
      data fpppp(15, 7),fpppp(15, 8)/-9.04003143d-07, 2.78348237d-07/
      data fpppp(15, 9),fpppp(15,10)/-6.55820230d-08, 2.51169354d-07/
      data fpppp(15,11),fpppp(15,12)/-3.67945198d-08,-2.03301349d-07/
      data fpppp(15,13),fpppp(15,14)/ 3.25065086d-07,-3.94881980d-07/
      data fpppp(15,15),fpppp(15,16)/ 5.75742924d-07,-1.09899112d-06/
      data fpppp(15,17),fpppp(15,18)/ 2.27183082d-07, 4.41520015d-07/
      data fpppp(15,19) /             9.73524741d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-2.51307138d-07,-6.95534415d-08/
      data fpppp(16, 3),fpppp(16, 4)/ 5.41137906d-07,-3.49452526d-07/
      data fpppp(16, 5),fpppp(16, 6)/ 2.73622639d-07,-4.85380728d-07/
      data fpppp(16, 7),fpppp(16, 8)/ 4.78414129d-07,-2.46564128d-07/
      data fpppp(16, 9),fpppp(16,10)/ 6.12170655d-08,-1.90913169d-07/
      data fpppp(16,11),fpppp(16,12)/-1.74339820d-09, 8.08489416d-08/
      data fpppp(16,13),fpppp(16,14)/-1.39734953d-07, 3.04223792d-07/
      data fpppp(16,15),fpppp(16,16)/-1.98250259d-07, 7.57556521d-07/
      data fpppp(16,17),fpppp(16,18)/-3.24456579d-07,-3.46117956d-07/
      data fpppp(16,19) /            -6.70736968d-07 /
 

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
      subroutine c_3(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(16,19,2),f(16,19),fpppp(16,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -2.15391000d-02 /
      data f( 1, 3),f( 1, 4) /  6.42785000d-02 ,  3.34638300d-01 /
      data f( 1, 5),f( 1, 6) /  4.93483600d-01 ,  2.21527000d-01 /
      data f( 1, 7),f( 1, 8) /  6.29549000d-02 ,  2.36371000d-02 /
      data f( 1, 9),f( 1,10) /  2.44811000d-02 ,  1.38948000d-02 /
      data f( 1,11),f( 1,12) /  2.36653000d-02 ,  1.41490000d-02 /
      data f( 1,13),f( 1,14) /  4.17270000d-03 ,  2.84104000d-02 /
      data f( 1,15),f( 1,16) /  1.74682500d-01 ,  5.13503000d-01 /
      data f( 1,17),f( 1,18) /  2.84070400d-01 ,  6.34407000d-02 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.86359000d-02 /
      data f( 2, 3),f( 2, 4) /  3.14454000d-02 ,  1.60206700d-01 /
      data f( 2, 5),f( 2, 6) /  2.24520600d-01 ,  1.13661800d-01 /
      data f( 2, 7),f( 2, 8) /  3.86388000d-02 ,  1.97936000d-02 /
      data f( 2, 9),f( 2,10) /  1.53626000d-02 ,  6.11270000d-03 /
      data f( 2,11),f( 2,12) /  1.21216000d-02 ,  1.19911000d-02 /
      data f( 2,13),f( 2,14) / -1.79070000d-03 ,  5.41300000d-04 /
      data f( 2,15),f( 2,16) /  9.21577000d-02 ,  2.41234400d-01 /
      data f( 2,17),f( 2,18) /  1.50355900d-01 ,  3.19696000d-02 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -1.12929000d-02 /
      data f( 3, 3),f( 3, 4) /  1.15140000d-02 ,  7.36349000d-02 /
      data f( 3, 5),f( 3, 6) /  9.65074000d-02 ,  5.93713000d-02 /
      data f( 3, 7),f( 3, 8) /  1.94149000d-02 ,  1.18618000d-02 /
      data f( 3, 9),f( 3,10) / -1.68900000d-04 , -7.25120000d-03 /
      data f( 3,11),f( 3,12) / -2.45850000d-03 ,  2.46740000d-03 /
      data f( 3,13),f( 3,14) / -5.07960000d-03 , -2.75810000d-03 /
      data f( 3,15),f( 3,16) /  4.10606000d-02 ,  1.13378000d-01 /
      data f( 3,17),f( 3,18) /  8.14511000d-02 ,  2.37532000d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 , -6.79920000d-03 /
      data f( 4, 3),f( 4, 4) /  2.40280000d-03 ,  3.03586000d-02 /
      data f( 4, 5),f( 4, 6) /  4.17241000d-02 ,  4.06895000d-02 /
      data f( 4, 7),f( 4, 8) /  1.85608000d-02 ,  2.01978000d-02 /
      data f( 4, 9),f( 4,10) / -2.19500000d-04 , -5.55890000d-03 /
      data f( 4,11),f( 4,12) / -1.99930000d-03 ,  6.34090000d-03 /
      data f( 4,13),f( 4,14) /  3.45570000d-03 ,  5.05500000d-03 /
      data f( 4,15),f( 4,16) /  1.30642000d-02 ,  4.79886000d-02 /
      data f( 4,17),f( 4,18) /  3.90486000d-02 ,  1.96624000d-02 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 , -6.08600000d-04 /
      data f( 5, 3),f( 5, 4) / -2.34610000d-03 ,  6.99340000d-03 /
      data f( 5, 5),f( 5, 6) /  2.47451000d-02 ,  3.48103000d-02 /
      data f( 5, 7),f( 5, 8) /  1.81555000d-02 ,  1.52726000d-02 /
      data f( 5, 9),f( 5,10) / -3.78300000d-04 , -4.28250000d-03 /
      data f( 5,11),f( 5,12) / -1.64890000d-03 ,  4.67130000d-03 /
      data f( 5,13),f( 5,14) /  7.05550000d-03 ,  3.62580000d-03 /
      data f( 5,15),f( 5,16) /  5.17090000d-03 ,  1.65195000d-02 /
      data f( 5,17),f( 5,18) /  1.59423000d-02 ,  9.83480000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 , -5.87660000d-03 /
      data f( 6, 3),f( 6, 4) / -8.33910000d-03 ,  3.07530000d-03 /
      data f( 6, 5),f( 6, 6) /  2.03484000d-02 ,  1.66740000d-02 /
      data f( 6, 7),f( 6, 8) /  2.58843000d-02 ,  1.10630000d-02 /
      data f( 6, 9),f( 6,10) / -5.62000000d-04 , -3.32800000d-03 /
      data f( 6,11),f( 6,12) / -1.38740000d-03 ,  3.30740000d-03 /
      data f( 6,13),f( 6,14) /  1.14383000d-02 ,  8.44460000d-03 /
      data f( 6,15),f( 6,16) /  7.96900000d-03 ,  9.80180000d-03 /
      data f( 6,17),f( 6,18) /  1.22730000d-03 ,  1.55220000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.08031000d-02 /
      data f( 7, 3),f( 7, 4) /  7.91700000d-04 ,  6.76610000d-03 /
      data f( 7, 5),f( 7, 6) /  1.48360000d-02 ,  2.15690000d-02 /
      data f( 7, 7),f( 7, 8) /  2.11942000d-02 ,  7.73540000d-03 /
      data f( 7, 9),f( 7,10) / -7.21200000d-04 , -2.61520000d-03 /
      data f( 7,11),f( 7,12) / -1.19140000d-03 ,  2.28900000d-03 /
      data f( 7,13),f( 7,14) /  8.30330000d-03 ,  1.66793000d-02 /
      data f( 7,15),f( 7,16) /  9.39770000d-03 ,  1.06419000d-02 /
      data f( 7,17),f( 7,18) /  2.79090000d-03 , -4.77060000d-03 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  7.67560000d-03 /
      data f( 8, 3),f( 8, 4) /  1.62391000d-02 ,  1.47395000d-02 /
      data f( 8, 5),f( 8, 6) /  2.07541000d-02 ,  2.70050000d-02 /
      data f( 8, 7),f( 8, 8) /  1.49388000d-02 ,  5.21110000d-03 /
      data f( 8, 9),f( 8,10) / -8.32300000d-04 , -2.08670000d-03 /
      data f( 8,11),f( 8,12) / -1.03980000d-03 ,  1.53390000d-03 /
      data f( 8,13),f( 8,14) /  5.96730000d-03 ,  1.20385000d-02 /
      data f( 8,15),f( 8,16) /  1.17513000d-02 ,  1.12282000d-02 /
      data f( 8,17),f( 8,18) /  6.21040000d-03 ,  1.86810000d-03 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  3.70820000d-03 /
      data f( 9, 3),f( 9, 4) /  1.11748000d-02 ,  1.69942000d-02 /
      data f( 9, 5),f( 9, 6) /  1.74972000d-02 ,  1.29890000d-02 /
      data f( 9, 7),f( 9, 8) /  6.89630000d-03 ,  2.00080000d-03 /
      data f( 9, 9),f( 9,10) / -9.04600000d-04 , -1.39590000d-03 /
      data f( 9,11),f( 9,12) / -8.26500000d-04 ,  5.77000000d-04 /
      data f( 9,13),f( 9,14) /  2.92000000d-03 ,  6.01930000d-03 /
      data f( 9,15),f( 9,16) /  8.83130000d-03 ,  9.21870000d-03 /
      data f( 9,17),f( 9,18) /  6.03340000d-03 ,  1.65690000d-03 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.74160000d-03 /
      data f(10, 3),f(10, 4) /  5.30720000d-03 ,  8.00650000d-03 /
      data f(10, 5),f(10, 6) /  8.07070000d-03 ,  5.76880000d-03 /
      data f(10, 7),f(10, 8) /  2.76010000d-03 ,  3.86700000d-04 /
      data f(10, 9),f(10,10) / -8.52500000d-04 , -9.86300000d-04 /
      data f(10,11),f(10,12) / -6.73000000d-04 ,  7.94000000d-05 /
      data f(10,13),f(10,14) /  1.29560000d-03 ,  2.82770000d-03 /
      data f(10,15),f(10,16) /  4.14240000d-03 ,  4.28220000d-03 /
      data f(10,17),f(10,18) /  2.73410000d-03 ,  7.73400000d-04 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  6.43700000d-04 /
      data f(11, 3),f(11, 4) /  1.95310000d-03 ,  2.84860000d-03 /
      data f(11, 5),f(11, 6) /  2.68490000d-03 ,  1.67160000d-03 /
      data f(11, 7),f(11, 8) /  4.55400000d-04 , -3.85500000d-04 /
      data f(11, 9),f(11,10) / -6.13600000d-04 , -6.52400000d-04 /
      data f(11,11),f(11,12) / -5.15100000d-04 , -1.86200000d-04 /
      data f(11,13),f(11,14) /  3.27700000d-04 ,  9.35000000d-04 /
      data f(11,15),f(11,16) /  1.40580000d-03 ,  1.39770000d-03 /
      data f(11,17),f(11,18) /  8.27000000d-04 ,  2.50300000d-04 /
      data f(11,19) /           0.00000000d+00 /
      data f(12, 1),f(12, 2) /  0.00000000d+00 ,  1.98700000d-04 /
      data f(12, 3),f(12, 4) /  5.85500000d-04 ,  7.65000000d-04 /
      data f(12, 5),f(12, 6) /  5.60900000d-04 ,  1.16300000d-04 /
      data f(12, 7),f(12, 8) / -2.81800000d-04 , -3.73300000d-04 /
      data f(12, 9),f(12,10) / -3.79800000d-04 , -3.76100000d-04 /
      data f(12,11),f(12,12) / -3.52700000d-04 , -2.43100000d-04 /
      data f(12,13),f(12,14) / -4.08000000d-05 ,  1.80900000d-04 /
      data f(12,15),f(12,16) /  3.22700000d-04 ,  2.68300000d-04 /
      data f(12,17),f(12,18) /  1.82900000d-04 ,  7.42000000d-05 /
      data f(12,19) /           0.00000000d+00 /
      data f(13, 1),f(13, 2) /  0.00000000d+00 , -4.39000000d-05 /
      data f(13, 3),f(13, 4) / -7.73000000d-05 , -1.29300000d-04 /
      data f(13, 5),f(13, 6) / -1.44400000d-04 , -1.40900000d-04 /
      data f(13, 7),f(13, 8) / -1.40100000d-04 , -1.39100000d-04 /
      data f(13, 9),f(13,10) / -1.38700000d-04 , -1.36300000d-04 /
      data f(13,11),f(13,12) / -1.22700000d-04 , -9.75000000d-05 /
      data f(13,13),f(13,14) / -6.42000000d-05 , -2.90000000d-05 /
      data f(13,15),f(13,16) / -2.70000000d-06 ,  3.60000000d-06 /
      data f(13,17),f(13,18) /  5.70000000d-06 ,  3.40000000d-06 /
      data f(13,19) /           0.00000000d+00 /
      data f(14, 1),f(14, 2) /  0.00000000d+00 , -6.60000000d-06 /
      data f(14, 3),f(14, 4) / -1.87000000d-05 , -2.46000000d-05 /
      data f(14, 5),f(14, 6) / -2.32000000d-05 , -2.12000000d-05 /
      data f(14, 7),f(14, 8) / -2.43000000d-05 , -2.65000000d-05 /
      data f(14, 9),f(14,10) / -2.72000000d-05 , -2.40000000d-05 /
      data f(14,11),f(14,12) / -1.96000000d-05 , -1.40000000d-05 /
      data f(14,13),f(14,14) / -8.60000000d-06 , -4.90000000d-06 /
      data f(14,15),f(14,16) / -4.00000000d-06 , -5.30000000d-06 /
      data f(14,17),f(14,18) / -6.10000000d-06 , -2.00000000d-06 /
      data f(14,19) /           0.00000000d+00 /
      data f(15, 1),f(15, 2) /  0.00000000d+00 , -1.70000000d-06 /
      data f(15, 3),f(15, 4) / -4.10000000d-06 , -5.90000000d-06 /
      data f(15, 5),f(15, 6) / -7.00000000d-06 , -8.30000000d-06 /
      data f(15, 7),f(15, 8) / -1.00000000d-05 , -1.06000000d-05 /
      data f(15, 9),f(15,10) / -1.09000000d-05 , -1.08000000d-05 /
      data f(15,11),f(15,12) / -1.03000000d-05 , -8.80000000d-06 /
      data f(15,13),f(15,14) / -7.40000000d-06 , -6.10000000d-06 /
      data f(15,15),f(15,16) / -5.20000000d-06 , -4.40000000d-06 /
      data f(15,17),f(15,18) / -3.40000000d-06 , -1.30000000d-06 /
      data f(15,19) /           0.00000000d+00 /
      data f(16, 1),f(16, 2) /  0.00000000d+00 , -8.00000000d-07 /
      data f(16, 3),f(16, 4) / -2.30000000d-06 , -4.50000000d-06 /
      data f(16, 5),f(16, 6) / -6.70000000d-06 , -8.40000000d-06 /
      data f(16, 7),f(16, 8) / -8.70000000d-06 , -8.00000000d-06 /
      data f(16, 9),f(16,10) / -6.60000000d-06 , -5.40000000d-06 /
      data f(16,11),f(16,12) / -4.60000000d-06 , -4.40000000d-06 /
      data f(16,13),f(16,14) / -5.40000000d-06 , -4.90000000d-06 /
      data f(16,15),f(16,16) / -3.90000000d-06 , -3.00000000d-06 /
      data f(16,17),f(16,18) / -1.70000000d-06 ,  1.00000000d-07 /
      data f(16,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00, 4.24037949d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 2.70737540d-01, 8.28915102d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 3.77310683d-01, 2.70170364d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 3.29607268d+00,-5.63197669d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/ 5.18594786d+00,-7.13978297d-03/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.78722991d+00, 3.27421554d-03/
      data fpp( 1, 7,1),fpp( 1, 7,2)/-1.62431276d-01, 8.45990819d-04/
      data fpp( 1, 8,1),fpp( 1, 8,2)/-5.46311009d-01, 4.97079186d-04/
      data fpp( 1, 9,1),fpp( 1, 9,2)/-6.64475137d-01,-4.24599562d-04/
      data fpp( 1,10,1),fpp( 1,10,2)/-6.13751360d-01, 5.15501061d-04/
      data fpp( 1,11,1),fpp( 1,11,2)/-4.89161311d-01,-4.15996683d-04/
      data fpp( 1,12,1),fpp( 1,12,2)/-6.54762718d-01,-8.72232816d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/-1.27370345d-01, 4.23285996d-04/
      data fpp( 1,14,1),fpp( 1,14,2)/ 9.68697065d-01, 3.68418345d-04/
      data fpp( 1,15,1),fpp( 1,15,2)/ 9.87084151d-01, 5.42510463d-03/
      data fpp( 1,16,1),fpp( 1,16,2)/ 5.60972688d+00,-1.05159328d-02/
      data fpp( 1,17,1),fpp( 1,17,2)/ 2.54278191d+00, 2.54344077d-03/
      data fpp( 1,18,1),fpp( 1,18,2)/ 1.05564928d+00, 8.70343781d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 3.40652411d-03/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00, 6.43643846d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 1.55964920d-01, 5.74972309d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 3.17081134d-01, 1.17949892d-03/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.22555463d+00,-5.72167985d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/ 3.58533429d+00,-2.75767098d-03/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.34194269d+00, 1.09248990d-03/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 4.28975527d-02, 5.37859383d-04/
      data fpp( 2, 8,1),fpp( 2, 8,2)/-2.31805482d-01, 1.26740567d-04/
      data fpp( 2, 9,1),fpp( 2, 9,2)/-2.46719727d-01,-1.79969653d-04/
      data fpp( 2,10,1),fpp( 2,10,2)/-2.23032280d-01, 3.04004045d-04/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.53192379d-01,-1.20518526d-04/
      data fpp( 2,12,1),fpp( 2,12,2)/-2.81059563d-01,-1.90293942d-04/
      data fpp( 2,13,1),fpp( 2,13,2)/-2.15681024d-03, 6.26162920d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 5.78193369d-01, 9.06656774d-04/
      data fpp( 2,15,1),fpp( 2,15,2)/ 7.99259197d-01, 1.66782061d-03/
      data fpp( 2,16,1),fpp( 2,16,2)/ 3.70872125d+00,-4.13032123d-03/
      data fpp( 2,17,1),fpp( 2,17,2)/ 1.69053368d+00, 4.56152292d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 5.89258933d-01, 6.55244060d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 2.10787347d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 3.27688382d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-2.28627219d-01, 2.63333236d-04/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.89619780d-01, 6.64966675d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/ 9.80678789d-01,-5.64359935d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 1.61518500d+00,-7.62430934d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 8.81204343d-01, 1.35676721d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 7.54671066d-01, 5.38942246d-04/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 8.60287936d-01,-2.25138655d-04/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 6.89404044d-01, 9.29563756d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/ 6.68610481d-01, 1.50217153d-04/
      data fpp( 3,11,1),fpp( 3,11,2)/ 6.46470826d-01, 1.86750127d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 6.74130971d-01,-2.16925204d-04/
      data fpp( 3,13,1),fpp( 3,13,2)/ 5.37172586d-01, 1.00651803d-04/
      data fpp( 3,14,1),fpp( 3,14,2)/ 4.03984458d-01, 4.06427993d-04/
      data fpp( 3,15,1),fpp( 3,15,2)/ 5.30034059d-01, 7.63468226d-04/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.21721813d+00,-1.75037890d-03/
      data fpp( 3,17,1),fpp( 3,17,2)/ 4.16538385d-01,-1.66106435d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/ 7.55199849d-02, 2.70561470d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 9.71046765d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 1.53727615d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/ 3.31148955d-01, 1.17528771d-04/
      data fpp( 4, 3,1),fpp( 4, 3,2)/ 1.47469745d-01, 3.36229302d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 3.46055212d-01,-3.37217979d-04/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 9.38410706d-01, 1.72246122d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 4.74544941d-01,-4.75686470d-04/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-3.06111815d-01, 6.19875268d-04/
      data fpp( 4, 8,1),fpp( 4, 8,2)/-7.69176264d-01,-5.77872604d-04/
      data fpp( 4, 9,1),fpp( 4, 9,2)/-1.88761449d-01, 3.68357146d-04/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.92979643d-01, 9.11801860d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.76795927d-01, 1.29110779d-04/
      data fpp( 4,12,1),fpp( 4,12,2)/-4.05884321d-01,-2.38725136d-04/
      data fpp( 4,13,1),fpp( 4,13,2)/-3.72903533d-01, 1.52265764d-04/
      data fpp( 4,14,1),fpp( 4,14,2)/-5.27256200d-01,-1.01267921d-04/
      data fpp( 4,15,1),fpp( 4,15,2)/ 5.45709567d-01, 6.37399918d-04/
      data fpp( 4,16,1),fpp( 4,16,2)/ 7.92456240d-01,-8.33419751d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 6.18657783d-01, 6.44150869d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/-2.72498872d-01,-5.10125962d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 1.23063298d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00,-1.27936713d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/-8.41433601d-01,-2.21115734d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/-2.25153760d-01, 1.48649007d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 6.21765365d-01, 9.21355463d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 3.01817172d-01,-1.24591920d-05/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-8.58994105d-01,-5.03488778d-04/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 5.37096194d-01, 4.23214305d-04/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 2.27237118d-01,-3.63054442d-04/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 4.94117505d-02, 2.62923464d-04/
      data fpp( 5,10,1),fpp( 5,10,2)/ 4.09230924d-02, 1.61625853d-05/
      data fpp( 5,11,1),fpp( 5,11,2)/ 4.43928797d-02, 6.46941945d-05/
      data fpp( 5,12,1),fpp( 5,12,2)/ 1.17941314d-01,-5.37433633d-05/
      data fpp( 5,13,1),fpp( 5,13,2)/ 2.14116547d-01,-8.58807412d-05/
      data fpp( 5,14,1),fpp( 5,14,2)/ 3.18695341d-01, 4.84323281d-05/
      data fpp( 5,15,1),fpp( 5,15,2)/ 3.02592672d-01, 1.90639429d-04/
      data fpp( 5,16,1),fpp( 5,16,2)/ 7.01001912d-01,-2.22780043d-04/
      data fpp( 5,17,1),fpp( 5,17,2)/ 3.26048346d-03,-1.50672577d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/ 1.53955505d-01,-4.87689264d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-1.34950368d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00,-6.74298170d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 1.31579545d+00, 2.80266340d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 5.66530295d-01, 1.60169281d-04/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 8.39483300d-02, 1.63910242d-04/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-2.58334395d-01,-4.64288251d-04/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.12286648d+00, 4.36392760d-04/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-6.22157960d-01,-5.08200789d-04/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-3.24322079d-02, 1.54514396d-04/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-1.26205533d-02, 8.19212064d-05/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.89977263d-02, 4.93407789d-05/
      data fpp( 6,11,1),fpp( 6,11,2)/-1.41105924d-02, 3.11167802d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/-2.00259358d-02, 1.03464509d-04/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.66112653d-01,-2.10803714d-04/
      data fpp( 6,14,1),fpp( 6,14,2)/ 1.89674834d-01, 7.22743475d-05/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.52370257d-01, 7.27923240d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.16246113d-01,-2.24939644d-04/
      data fpp( 6,17,1),fpp( 6,17,2)/ 6.26995283d-01, 2.02528250d-04/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.11573147d-01,-5.12093572d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-1.10316821d-04/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-5.53966327d-04/
      data fpp( 7, 2,1),fpp( 7, 2,2)/-1.12959319d+00,-2.52508347d-04/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 2.27602581d-01, 3.15129713d-04/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.83776315d-01,-4.88625065d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 5.64165408d-01, 6.05031287d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.77776816d-01,-5.55527449d-05/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 8.87006453d-02,-2.10307333d-04/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 3.47917138d-02, 1.11742077d-04/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 4.74546290d-03, 6.34710234d-05/
      data fpp( 7,10,1),fpp( 7,10,2)/-1.18718726d-03, 2.81298289d-05/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.22448990d-03, 2.30776608d-05/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.39874290d-02, 2.95552783d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.22664067d-01, 1.17134228d-04/
      data fpp( 7,14,1),fpp( 7,14,2)/-5.65009677d-01,-3.29790439d-04/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.01478356d-01, 2.62571529d-04/
      data fpp( 7,16,1),fpp( 7,16,2)/-3.23163648d-02,-2.08947678d-04/
      data fpp( 7,17,1),fpp( 7,17,2)/-6.94516170d-02, 2.75071825d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/ 5.86307083d-01, 1.16288948d-04/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 2.47263026d-04/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 1.04492797d-04/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 2.31497330d-01, 3.66714052d-05/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-5.29450618d-01,-1.97904418d-04/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-1.76663592d-01, 1.51160268d-04/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-2.83752237d-01, 4.41153462d-05/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-3.30609215d-01,-3.13443653d-04/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 3.25603784d-02, 1.10633265d-04/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 1.37603528d-02, 1.12205913d-05/
      data fpp( 8, 9,1),fpp( 8, 9,2)/ 8.53701751d-04, 6.55423694d-05/
      data fpp( 8,10,1),fpp( 8,10,2)/-3.89852466d-03, 1.39499311d-05/
      data fpp( 8,11,1),fpp( 8,11,2)/-1.44736719d-03, 1.67359063d-05/
      data fpp( 8,12,1),fpp( 8,12,2)/ 3.57121997d-03, 1.07144438d-05/
      data fpp( 8,13,1),fpp( 8,13,2)/-4.69361587d-03, 5.19883186d-05/
      data fpp( 8,14,1),fpp( 8,14,2)/ 1.39038875d-01,-1.20399718d-04/
      data fpp( 8,15,1),fpp( 8,15,2)/-1.14808165d-01, 4.81065543d-05/
      data fpp( 8,16,1),fpp( 8,16,2)/-2.50506540d-02,-8.61804991d-05/
      data fpp( 8,17,1),fpp( 8,17,2)/-7.08038155d-02, 2.69334421d-05/
      data fpp( 8,18,1),fpp( 8,18,2)/-2.89430187d-01, 1.89767308d-05/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 4.56116346d-05/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 9.12679765d-05/
      data fpp( 9, 2,1),fpp( 9, 2,2)/-4.39103913d-02, 3.83280470d-05/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.26084313d-01,-1.90761646d-05/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-7.53511322d-02,-6.08553886d-05/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 3.18275723d-03,-5.64862810d-05/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.47416054d-01,-1.38714873d-05/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 2.55297922d-02, 1.69022303d-05/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 1.02593347d-02, 1.80945662d-05/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 6.87413299d-04, 3.01255050d-05/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.44333239d-03, 6.24941380d-06/
      data fpp( 9,11,1),fpp( 9,11,2)/-1.41393382d-04, 8.51883981d-06/
      data fpp( 9,12,1),fpp( 9,12,2)/ 3.04137561d-03, 9.72122696d-06/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.36750639d-02, 8.96625236d-06/
      data fpp( 9,14,1),fpp( 9,14,2)/-1.22717873d-02,-2.08236383d-07/
      data fpp( 9,15,1),fpp( 9,15,2)/ 7.66531759d-03,-2.53713068d-05/
      data fpp( 9,16,1),fpp( 9,16,2)/-2.80186055d-02,-4.37825363d-05/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.59627451d-02,-1.38605479d-05/
      data fpp( 9,18,1),fpp( 9,18,2)/ 6.93145188d-02, 2.77527280d-05/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 6.60256360d-05/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 4.49450752d-05/
      data fpp(10, 2,1),fpp(10, 2,2)/ 1.91742356d-02, 1.86358497d-05/
      data fpp(10, 3,1),fpp(10, 3,2)/-5.01038308d-03,-1.00484739d-05/
      data fpp(10, 4,1),fpp(10, 4,2)/ 5.64781207d-02,-3.04199542d-05/
      data fpp(10, 5,1),fpp(10, 5,2)/ 3.96612082d-02,-2.63777093d-05/
      data fpp(10, 6,1),fpp(10, 6,2)/-4.21249927d-03,-6.03520876d-06/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.18067029d-02, 8.11054430d-06/
      data fpp(10, 8,1),fpp(10, 8,2)/ 5.05980826d-03, 1.17110315d-05/
      data fpp(10, 9,1),fpp(10, 9,2)/ 1.06164505d-03, 1.30973295d-05/
      data fpp(10,10,1),fpp(10,10,2)/-8.73145774d-04, 2.22365043d-06/
      data fpp(10,11,1),fpp(10,11,2)/-2.29559283d-04, 4.83406876d-06/
      data fpp(10,12,1),fpp(10,12,2)/ 1.48702758d-03, 4.78607452d-06/
      data fpp(10,13,1),fpp(10,13,2)/ 3.35211010d-03, 3.84963315d-06/
      data fpp(10,14,1),fpp(10,14,2)/ 1.60832738d-02,-1.23060710d-06/
      data fpp(10,15,1),fpp(10,15,2)/ 1.78131448d-02,-1.19712047d-05/
      data fpp(10,16,1),fpp(10,16,2)/ 2.73625761d-02,-2.13785740d-05/
      data fpp(10,17,1),fpp(10,17,2)/ 1.75685460d-02,-3.78849931d-06/
      data fpp(10,18,1),fpp(10,18,2)/-1.30391383d-02, 1.17765712d-05/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 2.79202144d-05/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 1.73546612d-05/
      data fpp(11, 2,1),fpp(11, 2,2)/-1.25053530d-03, 6.85367757d-06/
      data fpp(11, 3,1),fpp(11, 3,2)/ 1.26995289d-02,-4.82737148d-06/
      data fpp(11, 4,1),fpp(11, 4,2)/ 2.80107133d-03,-1.23781916d-05/
      data fpp(11, 5,1),fpp(11, 5,2)/ 8.20924484d-03,-9.21186194d-06/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.55053545d-02,-1.75036059d-06/
      data fpp(11, 7,1),fpp(11, 7,2)/ 5.84523579d-03, 4.03930429d-06/
      data fpp(11, 8,1),fpp(11, 8,2)/ 3.46742249d-03, 8.11114342d-06/
      data fpp(11, 9,1),fpp(11, 9,2)/-2.01252831d-04, 2.84122033d-07/
      data fpp(11,10,1),fpp(11,10,2)/ 2.35906985d-05, 2.11036845d-06/
      data fpp(11,11,1),fpp(11,11,2)/ 1.24128123d-04, 1.84040416d-06/
      data fpp(11,12,1),fpp(11,12,2)/ 7.67200216d-04, 2.02401489d-06/
      data fpp(11,13,1),fpp(11,13,2)/ 2.49475250d-03, 1.16353627d-06/
      data fpp(11,14,1),fpp(11,14,2)/ 2.24084423d-03,-1.07415998d-06/
      data fpp(11,15,1),fpp(11,15,2)/ 4.72902482d-03,-5.05689636d-06/
      data fpp(11,16,1),fpp(11,16,2)/ 2.77661030d-03,-7.43225458d-06/
      data fpp(11,17,1),fpp(11,17,2)/ 2.73203046d-03, 1.02991470d-06/
      data fpp(11,18,1),fpp(11,18,2)/ 5.43988290d-03, 2.95259580d-06/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 6.74370210d-06/
      data fpp(12, 1,1),fpp(12, 1,2)/ 0.00000000d+00, 5.77626070d-06/
      data fpp(12, 2,1),fpp(12, 2,2)/ 1.49750556d-03, 1.99847861d-06/
      data fpp(12, 3,1),fpp(12, 3,2)/ 1.88826731d-03,-2.48417513d-06/
      data fpp(12, 4,1),fpp(12, 4,2)/ 6.10079402d-03,-4.49977807d-06/
      data fpp(12, 5,1),fpp(12, 5,2)/ 5.78501249d-03,-2.53271257d-06/
      data fpp(12, 6,1),fpp(12, 6,2)/ 3.19668137d-03, 2.00628362d-07/
      data fpp(12, 7,1),fpp(12, 7,2)/ 2.43235394d-03, 4.52019912d-06/
      data fpp(12, 8,1),fpp(12, 8,2)/-1.03898218d-04, 1.14575146d-07/
      data fpp(12, 9,1),fpp(12, 9,2)/-3.79033730d-04, 1.21500292d-07/
      data fpp(12,10,1),fpp(12,10,2)/-6.03617021d-04, 1.14236856d-08/
      data fpp(12,11,1),fpp(12,11,2)/-1.58953209d-04, 1.01480497d-06/
      data fpp(12,12,1),fpp(12,12,2)/ 4.52971554d-04, 1.10135645d-06/
      data fpp(12,13,1),fpp(12,13,2)/ 1.05447991d-03, 1.41769225d-07/
      data fpp(12,14,1),fpp(12,14,2)/ 2.27974931d-03,-5.04433353d-07/
      data fpp(12,15,1),fpp(12,15,2)/ 2.95475597d-03,-2.91803581d-06/
      data fpp(12,16,1),fpp(12,16,2)/ 3.65338266d-03, 4.04576610d-07/
      data fpp(12,17,1),fpp(12,17,2)/ 1.81533214d-03,-5.60270626d-07/
      data fpp(12,18,1),fpp(12,18,2)/-3.92393280d-04, 4.38505893d-07/
      data fpp(12,19,1),fpp(12,19,2)/ 0.00000000d+00, 8.76247053d-07/
      data fpp(13, 1,1),fpp(13, 1,2)/ 0.00000000d+00, 3.53025521d-07/
      data fpp(13, 2,1),fpp(13, 2,2)/ 1.71509667d-05, 1.90948958d-07/
      data fpp(13, 3,1),fpp(13, 3,2)/ 4.19833602d-04,-4.86821351d-07/
      data fpp(13, 4,1),fpp(13, 4,2)/-6.55177141d-05, 6.40336448d-07/
      data fpp(13, 5,1),fpp(13, 5,2)/-2.03459875d-04, 1.39475558d-07/
      data fpp(13, 6,1),fpp(13, 6,2)/-2.22321344d-04,-8.22386815d-08/
      data fpp(13, 7,1),fpp(13, 7,2)/-5.23079712d-04, 2.74791678d-08/
      data fpp(13, 8,1),fpp(13, 8,2)/-1.63216591d-04,-1.56779898d-08/
      data fpp(13, 9,1),fpp(13, 9,2)/-1.21272396d-04,-7.67208682d-10/
      data fpp(13,10,1),fpp(13,10,2)/-7.77442877d-05, 1.38746825d-07/
      data fpp(13,11,1),fpp(13,11,2)/-1.54004434d-04, 1.17779911d-07/
      data fpp(13,12,1),fpp(13,12,2)/-1.86114770d-04, 8.61335330d-08/
      data fpp(13,13,1),fpp(13,13,2)/-1.29215984d-04, 2.36859573d-08/
      data fpp(13,14,1),fpp(13,14,2)/-1.69870058d-04,-6.68773624d-08/
      data fpp(13,15,1),fpp(13,15,2)/-1.83980329d-04,-2.90176508d-07/
      data fpp(13,16,1),fpp(13,16,2)/-3.83853131d-04, 2.75833939d-08/
      data fpp(13,17,1),fpp(13,17,2)/-1.46011661d-04,-7.21570676d-08/
      data fpp(13,18,1),fpp(13,18,2)/ 1.45638392d-04,-2.95512355d-09/
      data fpp(13,19,1),fpp(13,19,2)/ 0.00000000d+00, 1.79775618d-08/
      data fpp(14, 1,1),fpp(14, 1,2)/ 0.00000000d+00,-1.68395428d-07/
      data fpp(14, 2,1),fpp(14, 2,2)/-1.64556809d-05,-6.22091430d-08/
      data fpp(14, 3,1),fpp(14, 3,2)/-1.27334461d-04, 8.72320006d-08/
      data fpp(14, 4,1),fpp(14, 4,2)/-1.38938661d-05, 8.52811407d-08/
      data fpp(14, 5,1),fpp(14, 5,2)/ 1.55733819d-05, 9.64343647d-09/
      data fpp(14, 6,1),fpp(14, 6,2)/ 1.97733481d-05,-8.78548866d-08/
      data fpp(14, 7,1),fpp(14, 7,2)/ 1.01662167d-04, 3.57761100d-08/
      data fpp(14, 8,1),fpp(14, 8,2)/ 7.89888247d-06,-1.24955349d-09/
      data fpp(14, 9,1),fpp(14, 9,2)/-2.71594850d-06, 5.92221039d-08/
      data fpp(14,10,1),fpp(14,10,2)/-1.59086267d-05,-1.63886217d-09/
      data fpp(14,11,1),fpp(14,11,2)/ 6.13990534d-06, 1.93333448d-08/
      data fpp(14,12,1),fpp(14,12,2)/ 2.03085326d-05,-3.69451689d-09/
      data fpp(14,13,1),fpp(14,13,2)/ 1.40079957d-05,-1.65552772d-08/
      data fpp(14,14,1),fpp(14,14,2)/ 3.55855177d-05,-3.20843743d-08/
      data fpp(14,15,1),fpp(14,15,2)/ 4.88130003d-05,-2.31072257d-08/
      data fpp(14,16,1),fpp(14,16,2)/ 1.05618064d-04,-7.48672301d-09/
      data fpp(14,17,1),fpp(14,17,2)/ 4.42689123d-05, 8.30541177d-08/
      data fpp(14,18,1),fpp(14,18,2)/-3.64185348d-05,-3.07297479d-08/
      data fpp(14,19,1),fpp(14,19,2)/ 0.00000000d+00,-8.61351260d-08/
      data fpp(15, 1,1),fpp(15, 1,2)/ 0.00000000d+00,-1.96154314d-08/
      data fpp(15, 2,1),fpp(15, 2,2)/ 7.17568473d-08,-7.76913723d-09/
      data fpp(15, 3,1),fpp(15, 3,2)/ 2.35042430d-05, 8.69198029d-09/
      data fpp(15, 4,1),fpp(15, 4,2)/-7.90682142d-06, 9.00121606d-09/
      data fpp(15, 5,1),fpp(15, 5,2)/-1.63336529d-05,-2.69684453d-09/
      data fpp(15, 6,1),fpp(15, 6,2)/-1.69720483d-05,-1.02138379d-08/
      data fpp(15, 7,1),fpp(15, 7,2)/-3.58189552d-05, 1.95521963d-08/
      data fpp(15, 8,1),fpp(15, 8,2)/-1.34289387d-05,-1.99494722d-09/
      data fpp(15, 9,1),fpp(15, 9,2)/-1.06638104d-05, 6.42759260d-09/
      data fpp(15,10,1),fpp(15,10,2)/-7.27120548d-06, 2.84576808d-10/
      data fpp(15,11,1),fpp(15,11,2)/-1.12551878d-05, 1.64341002d-08/
      data fpp(15,12,1),fpp(15,12,2)/-1.25693606d-05,-6.02097747d-09/
      data fpp(15,13,1),fpp(15,13,2)/-8.41599899d-06, 1.64980971d-09/
      data fpp(15,14,1),fpp(15,14,2)/-1.04220126d-05,-6.57826137d-09/
      data fpp(15,15,1),fpp(15,15,2)/-1.11216723d-05, 6.63235782d-10/
      data fpp(15,16,1),fpp(15,16,2)/-2.39191243d-05,-2.07468176d-09/
      data fpp(15,17,1),fpp(15,17,2)/-9.31398776d-06, 1.96354912d-08/
      data fpp(15,18,1),fpp(15,18,2)/ 9.18574769d-06,-1.04672832d-08/
      data fpp(15,19,1),fpp(15,19,2)/ 0.00000000d+00,-2.57663584d-08/
      data fpp(16, 1,1),fpp(16, 1,2)/ 0.00000000d+00,-7.33529593d-09/
      data fpp(16, 2,1),fpp(16, 2,2)/ 6.43126443d-06,-6.32940814d-09/
      data fpp(16, 3,1),fpp(16, 3,2)/-6.85783579d-06,-9.34707151d-09/
      data fpp(16, 4,1),fpp(16, 4,2)/ 1.78519821d-05, 1.71769420d-09/
      data fpp(16, 5,1),fpp(16, 5,2)/ 2.80632550d-05, 2.47629473d-09/
      data fpp(16, 6,1),fpp(16, 6,2)/ 3.04245956d-05, 1.83771269d-08/
      data fpp(16, 7,1),fpp(16, 7,2)/ 3.81884062d-05, 8.01519778d-09/
      data fpp(16, 8,1),fpp(16, 8,2)/ 2.53305408d-05, 9.56208201d-09/
      data fpp(16, 9,1),fpp(16, 9,2)/ 2.39233338d-05,-4.26352580d-09/
      data fpp(16,10,1),fpp(16,10,2)/ 2.52431027d-05,-4.50797879d-09/
      data fpp(16,11,1),fpp(16,11,2)/ 2.79240225d-05,-1.70455904d-09/
      data fpp(16,12,1),fpp(16,12,2)/ 2.60921803d-05,-2.46737851d-08/
      data fpp(16,13,1),fpp(16,13,2)/ 1.88479995d-05, 2.83996992d-08/
      data fpp(16,14,1),fpp(16,14,2)/ 1.30163634d-05, 1.07498806d-09/
      data fpp(16,15,1),fpp(16,15,2)/ 6.59690760d-06,-2.69965150d-09/
      data fpp(16,16,1),fpp(16,16,2)/ 9.35170500d-06, 3.72361793d-09/
      data fpp(16,17,1),fpp(16,17,2)/-3.26489748d-08, 1.18051798d-08/
      data fpp(16,18,1),fpp(16,18,2)/-6.10680242d-06,-2.09443371d-08/
      data fpp(16,19,1),fpp(16,19,2)/ 0.00000000d+00,-4.20278315d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-2.90021221d-02,-6.44775294d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 4.49432700d-02,-4.59399539d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/-8.83004981d-02, 4.04804002d-02/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.33223035d-02, 1.77272948d-04/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.91154106d-03, 2.30983706d-03/
      data fpppp( 1,11),fpppp( 1,12)/-6.71891295d-03, 7.15432731d-03/
      data fpppp( 1,13),fpppp( 1,14)/ 1.96812306d-02,-5.17587475d-02/
      data fpppp( 1,15),fpppp( 1,16)/ 1.22692940d-01,-1.62757673d-01/
      data fpppp( 1,17),fpppp( 1,18)/ 6.69624928d-02,-1.03035578d-02/
      data fpppp( 1,19) /             1.40739199d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/-1.60406607d-02,-2.60825947d-03/
      data fpppp( 2, 3),fpppp( 2, 4)/ 2.67827763d-02, 3.18591280d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/-6.09787720d-02, 2.74062219d-02/
      data fpppp( 2, 7),fpppp( 2, 8)/ 8.01467234d-03, 1.99561482d-03/
      data fpppp( 2, 9),fpppp( 2,10)/-4.09804229d-04, 1.95970359d-03/
      data fpppp( 2,11),fpppp( 2,12)/-4.65986287d-03, 4.81732274d-03/
      data fpppp( 2,13),fpppp( 2,14)/ 9.79676814d-03,-2.59175497d-02/
      data fpppp( 2,15),fpppp( 2,16)/ 7.23163696d-02,-1.02044155d-01/
      data fpppp( 2,17),fpppp( 2,18)/ 4.02012735d-02,-3.74616911d-03/
      data fpppp( 2,19) /             5.50435149d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/ 1.28305609d-02, 8.22634888d-03/
      data fpppp( 3, 3),fpppp( 3, 4)/-9.23503384d-04, 5.83638526d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-2.58152054d-02, 1.53152241d-02/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.00115192d-03,-5.39082289d-03/
      data fpppp( 3, 9),fpppp( 3,10)/ 3.97209385d-03,-1.49213273d-03/
      data fpppp( 3,11),fpppp( 3,12)/ 1.91567157d-03,-3.18256561d-03/
      data fpppp( 3,13),fpppp( 3,14)/ 9.37479075d-04,-3.41135260d-04/
      data fpppp( 3,15),fpppp( 3,16)/ 1.59813257d-02,-2.99160996d-02/
      data fpppp( 3,17),fpppp( 3,18)/ 1.44112440d-02,-1.49196019d-04/
      data fpppp( 3,19) /             2.11544498d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/-1.39823785d-02,-5.42194484d-03/
      data fpppp( 4, 3),fpppp( 4, 4)/ 4.78046794d-03, 9.23595367d-03/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.80980809d-02,-2.16905618d-04/
      data fpppp( 4, 7),fpppp( 4, 8)/-4.17559752d-05, 1.94394679d-02/
      data fpppp( 4, 9),fpppp( 4,10)/-1.51073599d-02, 5.91199095d-03/
      data fpppp( 4,11),fpppp( 4,12)/-7.31648923d-03, 8.63763927d-03/
      data fpppp( 4,13),fpppp( 4,14)/-1.15099169d-02, 2.61620208d-02/
      data fpppp( 4,15),fpppp( 4,16)/-1.94990605d-02, 2.26107563d-03/
      data fpppp( 4,17),fpppp( 4,18)/-1.47779498d-02, 1.38092318d-02/
      data fpppp( 4,19) /             2.93603543d-02 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 2.69434189d-02, 1.43860486d-02/
      data fpppp( 5, 3),fpppp( 5, 4)/ 2.97519329d-03,-1.24484647d-02/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.31933735d-02, 5.47701735d-02/
      data fpppp( 5, 7),fpppp( 5, 8)/-4.24732259d-02, 1.27657677d-02/
      data fpppp( 5, 9),fpppp( 5,10)/-6.67822306d-04, 6.57241009d-05/
      data fpppp( 5,11),fpppp( 5,12)/ 1.12243263d-03,-3.50735774d-04/
      data fpppp( 5,13),fpppp( 5,14)/ 1.63811834d-03,-5.69752383d-03/
      data fpppp( 5,15),fpppp( 5,16)/ 1.39110892d-02,-2.50761183d-02/
      data fpppp( 5,17),fpppp( 5,18)/ 2.06243440d-02,-6.51507084d-03/
      data fpppp( 5,19) /            -1.28430922d-02 /
      data fpppp( 6, 1),fpppp( 6, 2)/-4.24656136d-02,-2.36554667d-02/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.31838443d-02,-1.30789191d-02/
      data fpppp( 6, 5),fpppp( 6, 6)/ 4.75497866d-02,-7.37112114d-02/
      data fpppp( 6, 7),fpppp( 6, 8)/ 5.97215400d-02,-2.50899371d-02/
      data fpppp( 6, 9),fpppp( 6,10)/ 6.44336270d-03,-2.25484331d-03/
      data fpppp( 6,11),fpppp( 6,12)/ 3.25186897d-03,-1.14007812d-02/
      data fpppp( 6,13),fpppp( 6,14)/ 2.19409734d-02,-2.22506599d-02/
      data fpppp( 6,15),fpppp( 6,16)/ 1.31917116d-02, 6.12350120d-03/
      data fpppp( 6,17),fpppp( 6,18)/-2.31577484d-02, 1.15484365d-02/
      data fpppp( 6,19) /             2.79724972d-02 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 6.08625246d-02, 3.06348401d-02/
      data fpppp( 7, 3),fpppp( 7, 4)/-3.41945468d-02, 2.20820245d-02/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.86806298d-02, 2.53006157d-02/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.20166517d-02, 3.54280744d-03/
      data fpppp( 7, 9),fpppp( 7,10)/-7.22817256d-04, 7.95277626d-04/
      data fpppp( 7,11),fpppp( 7,12)/-1.89763361d-03, 7.29633252d-03/
      data fpppp( 7,13),fpppp( 7,14)/-2.14728745d-02, 3.08141426d-02/
      data fpppp( 7,15),fpppp( 7,16)/-2.05339890d-02, 3.30484845d-03/
      data fpppp( 7,17),fpppp( 7,18)/ 1.31141633d-02,-1.41878646d-02/
      data fpppp( 7,19) /            -3.08866518d-02 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.92351398d-02,-1.34266837d-02/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.33951581d-02,-1.33298504d-02/
      data fpppp( 8, 5),fpppp( 8, 6)/ 2.33170302d-03, 7.61693829d-03/
      data fpppp( 8, 7),fpppp( 8, 8)/-8.19786190d-03, 2.25633218d-03/
      data fpppp( 8, 9),fpppp( 8,10)/-4.73864331d-04, 1.28390624d-04/
      data fpppp( 8,11),fpppp( 8,12)/ 3.92504868d-04,-1.54436432d-03/
      data fpppp( 8,13),fpppp( 8,14)/ 4.98794701d-03,-9.28758412d-03/
      data fpppp( 8,15),fpppp( 8,16)/ 8.30761756d-03,-3.32661302d-03/
      data fpppp( 8,17),fpppp( 8,18)/-3.13180584d-03, 5.48144380d-03/
      data fpppp( 8,19) /             1.16894242d-02 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.39379195d-03, 3.33627384d-03/
      data fpppp( 9, 3),fpppp( 9, 4)/-7.90458159d-03, 5.99624359d-03/
      data fpppp( 9, 5),fpppp( 9, 6)/ 7.17767309d-04,-4.92534841d-03/
      data fpppp( 9, 7),fpppp( 9, 8)/ 3.01645285d-03,-7.43514735d-04/
      data fpppp( 9, 9),fpppp( 9,10)/ 2.99518256d-04,-8.08774223d-06/
      data fpppp( 9,11),fpppp( 9,12)/-6.12062046d-05, 3.65762360d-04/
      data fpppp( 9,13),fpppp( 9,14)/-9.54788074d-04, 1.25855756d-03/
      data fpppp( 9,15),fpppp( 9,16)/-1.32640481d-03, 7.09800013d-04/
      data fpppp( 9,17),fpppp( 9,18)/ 1.35159177d-03,-1.72288290d-03/
      data fpppp( 9,19) /            -3.73556716d-03 /
      data fpppp(10, 1),fpppp(10, 2)/-1.58120604d-03,-7.18993088d-04/
      data fpppp(10, 3),fpppp(10, 4)/ 1.85564713d-03,-1.56320808d-03/
      data fpppp(10, 5),fpppp(10, 6)/-3.01139795d-04, 1.14435956d-03/
      data fpppp(10, 7),fpppp(10, 8)/-6.82723882d-04, 2.20570157d-04/
      data fpppp(10, 9),fpppp(10,10)/-3.46328579d-05, 4.17636171d-05/
      data fpppp(10,11),fpppp(10,12)/ 2.22810286d-05,-6.65077090d-05/
      data fpppp(10,13),fpppp(10,14)/ 2.52659547d-04,-2.92165607d-04/
      data fpppp(10,15),fpppp(10,16)/ 2.55925321d-04,-2.62362052d-04/
      data fpppp(10,17),fpppp(10,18)/-3.67084804d-04, 4.81882016d-04/
      data fpppp(10,19) /             1.05836610d-03 /
      data fpppp(11, 1),fpppp(11, 2)/ 5.04825290d-04, 2.27349835d-04/
      data fpppp(11, 3),fpppp(11, 4)/-5.02188657d-04, 3.50493482d-04/
      data fpppp(11, 5),fpppp(11, 6)/ 1.86125974d-05,-3.11667704d-04/
      data fpppp(11, 7),fpppp(11, 8)/ 2.10684520d-04,-9.41320528d-05/
      data fpppp(11, 9),fpppp(11,10)/ 8.83919696d-05,-2.58246945d-05/
      data fpppp(11,11),fpppp(11,12)/ 7.44844199d-06, 2.85830066d-05/
      data fpppp(11,13),fpppp(11,14)/-5.67116572d-05, 7.93759891d-05/
      data fpppp(11,15),fpppp(11,16)/-9.62669674d-05, 3.92561740d-05/
      data fpppp(11,17),fpppp(11,18)/ 5.37123525d-05,-8.89596476d-05/
      data fpppp(11,19) /            -1.86737882d-04 /
      data fpppp(12, 1),fpppp(12, 2)/-5.38067965d-05,-2.41588969d-05/
      data fpppp(12, 3),fpppp(12, 4)/ 8.40377552d-05,-8.26862261d-05/
      data fpppp(12, 5),fpppp(12, 6)/-2.49913452d-05, 4.62986320d-05/
      data fpppp(12, 7),fpppp(12, 8)/-5.07629616d-05, 5.04377308d-05/
      data fpppp(12, 9),fpppp(12,10)/-1.53209629d-05, 1.38792540d-05/
      data fpppp(12,11),fpppp(12,12)/-4.12270051d-08,-3.67868886d-06/
      data fpppp(12,13),fpppp(12,14)/ 1.41309981d-05,-1.54196409d-05/
      data fpppp(12,15),fpppp(12,16)/ 1.45318008d-05,-4.12903607d-05/
      data fpppp(12,17),fpppp(12,18)/-1.57099027d-06, 2.53938274d-05/
      data fpppp(12,19) /             5.60028030d-05 /
      data fpppp(13, 1),fpppp(13, 2)/ 1.54091399d-05, 6.21898264d-06/
      data fpppp(13, 3),fpppp(13, 4)/-1.71531703d-05, 9.11166158d-06/
      data fpppp(13, 5),fpppp(13, 6)/ 1.55107338d-06,-8.17111361d-06/
      data fpppp(13, 7),fpppp(13, 8)/ 1.42195671d-05,-9.06986561d-06/
      data fpppp(13, 9),fpppp(13,10)/ 2.98475978d-06,-2.77413877d-06/
      data fpppp(13,11),fpppp(13,12)/ 9.24500082d-07, 1.72512702d-06/
      data fpppp(13,13),fpppp(13,14)/-2.48446083d-06, 2.35954465d-06/
      data fpppp(13,15),fpppp(13,16)/-5.36108952d-06, 7.93906153d-06/
      data fpppp(13,17),fpppp(13,18)/-1.32300246d-07,-4.18134557d-06/
      data fpppp(13,19) /            -9.37962416d-06 /
      data fpppp(14, 1),fpppp(14, 2)/-3.83644741d-06,-1.53464767d-06/
      data fpppp(14, 3),fpppp(14, 4)/ 4.30965210d-06,-2.24479821d-06/
      data fpppp(14, 5),fpppp(14, 6)/-3.68860085d-07, 2.20420164d-06/
      data fpppp(14, 7),fpppp(14, 8)/-3.78661533d-06, 2.40313349d-06/
      data fpppp(14, 9),fpppp(14,10)/-8.37011422d-07, 7.90241366d-07/
      data fpppp(14,11),fpppp(14,12)/-2.09481426d-07,-4.25109947d-07/
      data fpppp(14,13),fpppp(14,14)/ 6.81771363d-07,-6.29291969d-07/
      data fpppp(14,15),fpppp(14,16)/ 1.33439415d-06,-2.09362976d-06/
      data fpppp(14,17),fpppp(14,18)/-4.91280235d-08, 1.12984413d-06/
      data fpppp(14,19) /             2.55611044d-06 /
      data fpppp(15, 1),fpppp(15, 2)/ 9.41951144d-07, 3.81005190d-07/
      data fpppp(15, 3),fpppp(15, 4)/-1.06432814d-06, 5.85694349d-07/
      data fpppp(15, 5),fpppp(15, 6)/ 1.00604723d-07,-5.20807077d-07/
      data fpppp(15, 7),fpppp(15, 8)/ 8.90112893d-07,-5.65429090d-07/
      data fpppp(15, 9),fpppp(15,10)/ 1.94110173d-07,-1.73363004d-07/
      data fpppp(15,11),fpppp(15,12)/ 5.67466061d-08, 1.06565151d-07/
      data fpppp(15,13),fpppp(15,14)/-1.54955146d-07, 1.43692918d-07/
      data fpppp(15,15),fpppp(15,16)/-3.41435299d-07, 4.96180746d-07/
      data fpppp(15,17),fpppp(15,18)/ 8.67622025d-10,-2.65975298d-07/
      data fpppp(15,19) /            -5.98095419d-07 /
      data fpppp(16, 1),fpppp(16, 2)/-7.22798778d-07,-3.00399041d-07/
      data fpppp(16, 3),fpppp(16, 4)/ 7.41173061d-07,-3.84358113d-07/
      data fpppp(16, 5),fpppp(16, 6)/-7.36533122d-08, 2.07975423d-07/
      data fpppp(16, 7),fpppp(16, 8)/-4.34100175d-07, 2.91124719d-07/
      data fpppp(16, 9),fpppp(16,10)/-4.33591972d-08, 4.59306272d-08/
      data fpppp(16,11),fpppp(16,12)/-5.86942643d-08,-8.19192843d-08/
      data fpppp(16,13),fpppp(16,14)/ 6.16310827d-08,-7.98523616d-08/
      data fpppp(16,15),fpppp(16,16)/ 2.22509178d-07,-2.59729155d-07/
      data fpppp(16,17),fpppp(16,18)/ 8.80583611d-08, 1.06107743d-07/
      data fpppp(16,19) /             2.18368020d-07 /
 
 
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
