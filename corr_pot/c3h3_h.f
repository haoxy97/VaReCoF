
      real*8 function c3h3_h(r, rpar, ipar)

c
c  3d potential surface for h2ccch + h modified from lbh
c  for variflex calculations
c
c  ***  coordinate system  ***
c
c               
c                                  r(cha)  = 2.05 au
c                                  r(chb)  = 2.05 au
c       ha                         r(chc)  = 2.05 au
c        \                         r(cacb) = 2.65 au
c         ca===cb==cc--hc          r(cacb) = 2.35 au
c        /      \                  ha-c-c = hb-c-c = 120.0 deg
c       hb       \                 hc-c-c = 180.0 deg 
c                 \                habcc  = 180.0 deg
c                  \               
c                   hd
c
c    cartesian coordinates for the propargyl radical:
c
c              x(au)        y(au)                z(au)
c
c    ca        0.000    0.0000000000d+00   -2.6500000000d+00
c    cb        0.000    0.0000000000d+00    0.0000000000d+00
c    cc        0.000    0.0000000000d+00    2.3500000000d+00
c    ha        0.000    1.7753520778d+00   -3.6750000000d+00
c    hb        0.000   -1.7753520778d+00   -3.6750000000d+00
c    hc        0.000    0.0000000000d+00    4.4000000000d+00
c
c in order of variflex input and in angstroms 
c              x(a)        y(a)                z(a)
c    ca        0.000    0.0000000000d+00   -1.4023000000d+00
c    cb        0.000    0.0000000000d+00    0.0000000000d+00
c    ha        0.000    0.9395000000d+00   -1.9447000000d+00
c    cc        0.000    0.0000000000d+00    1.2436000000d+00
c    hb        0.000   -0.9395000000d+00   -1.9447000000d+00
c    hc        0.000    0.0000000000d+00    2.3284000000d+00
c
c
c  ***   input   ***
c
c     rd :     r(cb-hd) au
c     alpha :  hd-cb-ca angle  (degrees)
c     phi :    dihedral angle between hacacb and hdcacb planes (degrees)
c
c     ranges :     4.0 < rd < 12.0
c                    0 < alpha < +180
c                  no limits on phi
c
c
c  ***   output   ***
c
c     energy - potential energy relative to h+c3h3 in au
c
c
      implicit double precision (a-h,o-z)

      include 'param2.fi'

      include 'commonpot.fi'

      dimension dmsst(nfrag),dmss(nfrag,natommx)

      dimension r(nfrag,natommx,ndim)
      dimension natom(nfrag)
      dimension r11(ndim),r12(ndim),r13(ndim),r23(ndim),r21(ndim)
      dimension tmp1(ndim),tmp2(ndim),tmp3(ndim)

      include 'data.fi'

      pi2 = asin(1.0d0)
      pi = 2.0d0*pi2

c
c     read in the initial parameters
c

c      ichab = ipotpar(iel,ieff,ich,1)
      phidvd = 90.0d0*pi/180.0d0
c      ibsc = ipotpar(iel,ieff,ich,2)
c      dmu12 = dmsst(1)*dmsst(2)/(dmsst(1)+dmsst(2))
c      dmu12p = dmsst(1)*(dmsst(1)+dmsst(2))/(2.0d0*dmsst(1)+dmsst(2))

c
c     for c3h5 + h -- just need to calculate r, theta, phi
c

      r112 = 0.0d0
      r122 = 0.0d0
      r132 = 0.0d0
	r212 = 0.0d0
      r232 = 0.0d0

      do 1100 idim = 1 , ndim

      r11(idim) = r(1,1,idim) - r(2,1,idim)
      r12(idim) = r(1,1,idim) - r(1,2,idim)
      r13(idim) = r(1,1,idim) - r(1,3,idim)
	r21(idim) = r(1,2,idim) - r(2,1,idim)
      r23(idim) = r(1,3,idim) - r(1,2,idim)
      r112 = r112 + r11(idim)**2
      r122 = r122 + r12(idim)**2
      r132 = r132 + r13(idim)**2
      r212 = r212 + r21(idim)**2
      r232 = r232 + r23(idim)**2

1100  continue

      rmepi = sqrt(r112)
      rch = sqrt(r232)
	rchp = sqrt(r132)
	rcc = sqrt(r122)
	rr = sqrt(r212)
c      rhh = sqrt((r(2,1,1)-r(1,3,1))**2 + (r(2,1,2)-r(1,3,2))**2 + 
c     $ (r(2,1,3)-r(1,3,3))**2)
c	rtest = sqrt((r(2,1,1)-r(1,4,1))**2 + (r(2,1,2)-r(1,4,2))**2 + 
c     $ (r(2,1,3)-r(1,4,3))**2)

      cthe = (r122 + r212 - r112)/(2.0d0*rcc*rr)
      cthep = (r232 + r122 - r132)/(2.0d0*rch*rcc)
      if (abs(cthe).gt.1.0d0) then
		write (7,*) 'error in cthe',cthe
	    if (cthe.gt.1.0d0) then 
			the = 0.0d0
		else
			the = pi
		endif
	else
	    the = acos(cthe)
	endif
      if (abs(cthep).gt.1.0d0) then
		write (7,*) 'error in cthep',cthep
	    if (cthep.gt.1.0d0) then 
			thep = 0.0d0
		else
			thep = pi
		endif
	else
	    thep = acos(cthep)
	endif
      call cross(r21,r12,tmp1)
      call cross(r12,r23,tmp2)
      tdot = 0.0d0
      do 1200 idim = 1 , ndim
         tdot = tdot + tmp1(idim)*tmp2(idim)
1200  continue
      rdenom = r122*sqrt(r232*r212)
      ctau = tdot/(rdenom*sin(the)*sin(thep))
      if (abs(ctau).gt.1.0d0) then 
	   write (6,*) 'error in ctau',ctau
         if (ctau.gt.1.0d0) then
	      tau = 0.0d0
	   else
	      tau = pi
	   endif
	else
         tau = acos(ctau)
	endif
c     still need to determine whether it is +/-
c     
         
      sum = 0.0d0
      call cross(r12,tmp2,tmp3)
      do 2570 idim = 1 , ndim
         sum = sum - tmp1(idim)*tmp3(idim)
 2570 continue
      ctaub = sum/(r122*sqrt(r232*r122*r212)*sin(the)*sin(thep))
      if (dabs(ctaub).le.1.0d0) go to 2580
      write (6,*) 'error in tau1, ctau1b = ',ctaub
      taub = 0.0d0
      if (ctaub.lt.1.0d0) taub = 2.0d0*dasin(1.0d0)
      go to 2590
 2580 continue
      taub = dacos(ctaub)
 2590 continue
c	write (6,*) 'tau1b test',ctau1b,tau1b*180./pi,
c     $ rbond*cautoang
         
      if (taub.gt.pi2) tau = -tau
c      write (6,*) 'tau test',tau,taub
      phi = tau*180.0d0/pi
      alpha = the*180.0d0/pi
c      ichtmp = 1
c      if (alpha.gt.90.0d0) ichtmp = 2
c	if (tau.lt.0.0d0) ichtmp = ichtmp + 2

c l ==> ch2 side
c r ==> ch side
c ichtmp = 1 ==> l,+
c ichtmp = 2 ==> r,+
c ichtmp = 3 ==> l,-
c ichtmp = 4 ==> r,-
c ichab = 0 ==> all
c ichab = 1 ==> l,+ and l,-
c ichab = 2 ==> l,+ and r,+
c ichab = 3 ==> l,+ only

c     if ((ichab.eq.1).and.((ichtmp.eq.2).or.(ichtmp.eq.4))) then
c   vtot = 1.0e20
c        go to 5000
cendif
c     if ((ichab.eq.2).and.((ichtmp.eq.3).or.(ichtmp.eq.4))) then
c   vtot = 1.0e20
c        go to 5000
cendif
cif ((ichab.eq.3).and.(ichtmp.ne.1)) then
c   vtot = 1.0e20
c        go to 5000
cendif
	rrold = rr
      if (rr.lt.4.0d0) then
	   rr = 4.00001
         call spl_c3h4(rr,alpha,phi,vtot)
c	   if (vtot.lt.0.0d0) write (6,*) 
c     $    'error in rr',rrold*cautoang,vtot*cautoicm
	   go to 5000
	endif
      if (rr.gt.12.0d0) then
c	   write (6,*) 'error in rr',rr*cautoang
	   rr = 11.9999
	   vtot = 0.0d0
	   go to 5000
	endif
c	if (((rmepi.lt.rtest).and.(phi.gt.90.0d0)).or.
c     $ ((rmepi.gt.rtest).and.(phi.lt.90.0d0))) then
c	  write (6,*) 'pot test',rmepi,rtest,theta,phi,vtot*cautoicm
c      write (6,*) 'rr test',rrold*cautoang,rmepi*cautoang
c	write (6,*) 'r1 geom',(r(1,1,idim)*cautoang,idim=1,3)
c	write (6,*) 'r2 geom',(r(1,2,idim)*cautoang,idim=1,3)
c	write (6,*) 'r3 geom',(r(1,3,idim)*cautoang,idim=1,3)
c	write (6,*) 'r4 geom',(r(1,4,idim)*cautoang,idim=1,3)
c	write (6,*) 'r5 geom',(r(2,1,idim)*cautoang,idim=1,3)
c	endif
      call spl_c3h4(rr,alpha,phi,vtot)
5000  continue
      c3h3_h = vtot
c	write (6,*) 'pot test',rr,rtest,theta,phi,vtot*cautoicm
      return
      end



      subroutine spl_c3h4(rd,alpha,phi,energy)
      implicit real*8 (a-h,o-z)
c
c     (2e,2o)-cas+1+2+q2/cc-pvdz potential
c     including basis set corrections
c
      data ezero  / -1.16460927d0 /
      data degrad / 57.29577951308232d 00 /

      data hay, haz /  1.7753520778d+00,   3.675d+00 /
      data hcz      / -4.4000000000d+00              /
      data caz,ccz  /  2.65d0 ,           -2.35d0    /

      xi = rd*sin(alpha/degrad) * sin(phi/degrad)
      yi = rd*sin(alpha/degrad) * cos(phi/degrad)
      zi = rd*cos(alpha/degrad)

c
c     evaluate h-h repulsive potentials:
c     these were subtracted from the total energies
c     before fitting inorder to reduce repulsive
c     spikes in the potential
c
      rha = dsqrt(xi**2 + (yi-hay)**2 + (zi-haz)**2)
      rhb = dsqrt(xi**2 + (yi+hay)**2 + (zi-haz)**2)
      rhc = dsqrt(xi**2 + yi**2 + (zi-hcz)**2)
      rca = dsqrt(xi**2 + yi**2 + (zi-caz)**2)
      rcc = dsqrt(xi**2 + yi**2 + (zi-ccz)**2)
      call hhspl(rha,eha)
      call hhspl(rhb,ehb)
      call hhspl(rhc,ehc)

c     calculate basis set correction

      call ebasis(rca,rcc,etz)

c
c     find location in spline grids
c
      call fnd_grd_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1)
c
      call c0_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c0)
      call c1_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c1)
      call c2_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c2)
      call c3_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c3)
      call c4_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c4)
      call c5_spl_c3h3_h(rd,alpha,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1,c5)
c
      p = phi / degrad
      energy = c0 + c1*cos( 2.0d0*p)  +  c2*cos(4.0d0*p)
     x            + c3*cos( 6.0d0*p)  +  c4*cos(8.0d0*p)
     x            + c5*cos(10.0d0*p)
     x            -  ezero
c
c     eha, ehb, and ehc are repulsive h-h potentials that were
c     subtracted before fitting to make a less-bumpy potential
c
c     etz is the basis set correction
c
      energy = energy + eha + ehb + ehc
     x                + etz

      return
      end
      subroutine fnd_grd_c3h3_h(xi,yi,ix,iy,delxi,
     $     delyi,xix,xixp1,yiy,yiyp1)
      implicit real*8 (a-h,o-z)
 
      dimension delx(10),dely(18),x(11),y(19)
 
      data x( 1), x( 2) /  4.00000000d+00 ,  4.50000000d+00 /
      data x( 3), x( 4) /  5.00000000d+00 ,  5.50000000d+00 /
      data x( 5), x( 6) /  6.00000000d+00 ,  6.50000000d+00 /
      data x( 7), x( 8) /  7.00000000d+00 ,  8.00000000d+00 /
      data x( 9), x(10) /  9.00000000d+00 ,  1.00000000d+01 /
      data x(11) /         1.20000000d+01 /
 
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
 
      data delx( 1), delx( 2) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 3), delx( 4) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 5), delx( 6) /  5.00000000d-01 ,  5.00000000d-01 /
      data delx( 7), delx( 8) /  1.00000000d+00 ,  1.00000000d+00 /
      data delx( 9), delx(10) /  1.00000000d+00 ,  2.00000000d+00 /
      data dely( 1), dely( 2) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 3), dely( 4) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 5), dely( 6) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 7), dely( 8) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely( 9), dely(10) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(11), dely(12) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(13), dely(14) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(15), dely(16) /  1.00000000d+01 ,  1.00000000d+01 /
      data dely(17), dely(18) /  1.00000000d+01 ,  1.00000000d+01 /
      data nptx,npty /  11 , 19 /

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
        do 20 i=1,npty
          if(yi .gt. y(i))go to 20
          iy=i-1
          go to 25
 20     continue
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
      subroutine c0_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) / -9.04354060d-01 , -1.03082512d+00 /
      data f( 1, 3),f( 1, 4) / -1.19722192d+00 , -1.24644216d+00 /
      data f( 1, 5),f( 1, 6) / -1.21790961d+00 , -1.18752055d+00 /
      data f( 1, 7),f( 1, 8) / -1.17292116d+00 , -1.16261080d+00 /
      data f( 1, 9),f( 1,10) / -1.15586709d+00 , -1.15298125d+00 /
      data f( 1,11),f( 1,12) / -1.15376524d+00 , -1.15869139d+00 /
      data f( 1,13),f( 1,14) / -1.16903711d+00 , -1.18488220d+00 /
      data f( 1,15),f( 1,16) / -1.19980257d+00 , -1.19624662d+00 /
      data f( 1,17),f( 1,18) / -1.15697032d+00 , -1.11903855d+00 /
      data f( 1,19) /          -1.11782120d+00 /
      data f( 2, 1),f( 2, 2) / -1.11223942d+00 , -1.15661557d+00 /
      data f( 2, 3),f( 2, 4) / -1.21644841d+00 , -1.22462222d+00 /
      data f( 2, 5),f( 2, 6) / -1.19725886d+00 , -1.18058063d+00 /
      data f( 2, 7),f( 2, 8) / -1.17145108d+00 , -1.16458054d+00 /
      data f( 2, 9),f( 2,10) / -1.16045017d+00 , -1.15883377d+00 /
      data f( 2,11),f( 2,12) / -1.15928199d+00 , -1.16189285d+00 /
      data f( 2,13),f( 2,14) / -1.16725765d+00 , -1.17538436d+00 /
      data f( 2,15),f( 2,16) / -1.18268522d+00 , -1.17764189d+00 /
      data f( 2,17),f( 2,18) / -1.15259285d+00 , -1.15152071d+00 /
      data f( 2,19) /          -1.21298613d+00 /
      data f( 3, 1),f( 3, 2) / -1.10950219d+00 , -1.14588243d+00 /
      data f( 3, 3),f( 3, 4) / -1.18426288d+00 , -1.19082183d+00 /
      data f( 3, 5),f( 3, 6) / -1.18176388d+00 , -1.17506582d+00 /
      data f( 3, 7),f( 3, 8) / -1.16947768d+00 , -1.16522324d+00 /
      data f( 3, 9),f( 3,10) / -1.16280811d+00 , -1.16189470d+00 /
      data f( 3,11),f( 3,12) / -1.16211824d+00 , -1.16346745d+00 /
      data f( 3,13),f( 3,14) / -1.16615512d+00 , -1.16996836d+00 /
      data f( 3,15),f( 3,16) / -1.17256479d+00 , -1.16613467d+00 /
      data f( 3,17),f( 3,18) / -1.14237271d+00 , -1.12747400d+00 /
      data f( 3,19) /          -1.14416908d+00 /
      data f( 4, 1),f( 4, 2) / -1.11030321d+00 , -1.13287576d+00 /
      data f( 4, 3),f( 4, 4) / -1.16081506d+00 , -1.17306455d+00 /
      data f( 4, 5),f( 4, 6) / -1.17376113d+00 , -1.17102933d+00 /
      data f( 4, 7),f( 4, 8) / -1.16773925d+00 , -1.16527016d+00 /
      data f( 4, 9),f( 4,10) / -1.16391690d+00 , -1.16340932d+00 /
      data f( 4,11),f( 4,12) / -1.16352270d+00 , -1.16421227d+00 /
      data f( 4,13),f( 4,14) / -1.16552945d+00 , -1.16723281d+00 /
      data f( 4,15),f( 4,16) / -1.16792224d+00 , -1.16329043d+00 /
      data f( 4,17),f( 4,18) / -1.14774352d+00 , -1.13367106d+00 /
      data f( 4,19) /          -1.13582305d+00 /
      data f( 5, 1),f( 5, 2) / -1.13792873d+00 , -1.14596048d+00 /
      data f( 5, 3),f( 5, 4) / -1.15959126d+00 , -1.16784497d+00 /
      data f( 5, 5),f( 5, 6) / -1.16960734d+00 , -1.16834773d+00 /
      data f( 5, 7),f( 5, 8) / -1.16649614d+00 , -1.16513024d+00 /
      data f( 5, 9),f( 5,10) / -1.16439325d+00 , -1.16412257d+00 /
      data f( 5,11),f( 5,12) / -1.16418670d+00 , -1.16454280d+00 /
      data f( 5,13),f( 5,14) / -1.16518636d+00 , -1.16593648d+00 /
      data f( 5,15),f( 5,16) / -1.16600779d+00 , -1.16338828d+00 /
      data f( 5,17),f( 5,18) / -1.15520632d+00 , -1.14507081d+00 /
      data f( 5,19) /          -1.14202920d+00 /
      data f( 6, 1),f( 6, 2) / -1.15184980d+00 , -1.15529307d+00 /
      data f( 6, 3),f( 6, 4) / -1.16179244d+00 , -1.16619754d+00 /
      data f( 6, 5),f( 6, 6) / -1.16735207d+00 , -1.16670911d+00 /
      data f( 6, 7),f( 6, 8) / -1.16570372d+00 , -1.16496825d+00 /
      data f( 6, 9),f( 6,10) / -1.16457650d+00 , -1.16443846d+00 /
      data f( 6,11),f( 6,12) / -1.16447928d+00 , -1.16466976d+00 /
      data f( 6,13),f( 6,14) / -1.16499072d+00 , -1.16532945d+00 /
      data f( 6,15),f( 6,16) / -1.16527504d+00 , -1.16395044d+00 /
      data f( 6,17),f( 6,18) / -1.16020059d+00 , -1.15483255d+00 /
      data f( 6,19) /          -1.15223493d+00 /
      data f( 7, 1),f( 7, 2) / -1.15869182d+00 , -1.16024154d+00 /
      data f( 7, 3),f( 7, 4) / -1.16330654d+00 , -1.16549620d+00 /
      data f( 7, 5),f( 7, 6) / -1.16610312d+00 , -1.16576513d+00 /
      data f( 7, 7),f( 7, 8) / -1.16523190d+00 , -1.16484164d+00 /
      data f( 7, 9),f( 7,10) / -1.16463894d+00 , -1.16456851d+00 /
      data f( 7,11),f( 7,12) / -1.16459602d+00 , -1.16470065d+00 /
      data f( 7,13),f( 7,14) / -1.16486686d+00 , -1.16502920d+00 /
      data f( 7,15),f( 7,16) / -1.16498447d+00 , -1.16436650d+00 /
      data f( 7,17),f( 7,18) / -1.16281770d+00 , -1.16062312d+00 /
      data f( 7,19) /          -1.15946880d+00 /
      data f( 8, 1),f( 8, 2) / -1.16346585d+00 , -1.16376891d+00 /
      data f( 8, 3),f( 8, 4) / -1.16440695d+00 , -1.16489759d+00 /
      data f( 8, 5),f( 8, 6) / -1.16503509d+00 , -1.16493932d+00 /
      data f( 8, 7),f( 8, 8) / -1.16479549d+00 , -1.16469177d+00 /
      data f( 8, 9),f( 8,10) / -1.16463963d+00 , -1.16462792d+00 /
      data f( 8,11),f( 8,12) / -1.16464159d+00 , -1.16467447d+00 /
      data f( 8,13),f( 8,14) / -1.16472084d+00 , -1.16476399d+00 /
      data f( 8,15),f( 8,16) / -1.16476001d+00 , -1.16465293d+00 /
      data f( 8,17),f( 8,18) / -1.16445041d+00 , -1.16426214d+00 /
      data f( 8,19) /          -1.16419782d+00 /
      data f( 9, 1),f( 9, 2) / -1.16438573d+00 , -1.16444603d+00 /
      data f( 9, 3),f( 9, 4) / -1.16457982d+00 , -1.16468738d+00 /
      data f( 9, 5),f( 9, 6) / -1.16471881d+00 , -1.16469480d+00 /
      data f( 9, 7),f( 9, 8) / -1.16465871d+00 , -1.16463258d+00 /
      data f( 9, 9),f( 9,10) / -1.16462241d+00 , -1.16462182d+00 /
      data f( 9,11),f( 9,12) / -1.16462768d+00 , -1.16463839d+00 /
      data f( 9,13),f( 9,14) / -1.16465118d+00 , -1.16466195d+00 /
      data f( 9,15),f( 9,16) / -1.16466376d+00 , -1.16465189d+00 /
      data f( 9,17),f( 9,18) / -1.16464142d+00 , -1.16466521d+00 /
      data f( 9,19) /          -1.16468965d+00 /
      data f(10, 1),f(10, 2) / -1.16455246d+00 , -1.16456546d+00 /
      data f(10, 3),f(10, 4) / -1.16459503d+00 , -1.16462222d+00 /
      data f(10, 5),f(10, 6) / -1.16463111d+00 , -1.16462617d+00 /
      data f(10, 7),f(10, 8) / -1.16461953d+00 , -1.16461515d+00 /
      data f(10, 9),f(10,10) / -1.16461362d+00 , -1.16461425d+00 /
      data f(10,11),f(10,12) / -1.16461607d+00 , -1.16461915d+00 /
      data f(10,13),f(10,14) / -1.16462233d+00 , -1.16462469d+00 /
      data f(10,15),f(10,16) / -1.16462524d+00 , -1.16462440d+00 /
      data f(10,17),f(10,18) / -1.16462855d+00 , -1.16464396d+00 /
      data f(10,19) /          -1.16465482d+00 /
      data f(11, 1),f(11, 2) / -1.16460277d+00 , -1.16460317d+00 /
      data f(11, 3),f(11, 4) / -1.16460587d+00 , -1.16460779d+00 /
      data f(11, 5),f(11, 6) / -1.16460904d+00 , -1.16460922d+00 /
      data f(11, 7),f(11, 8) / -1.16460976d+00 , -1.16460992d+00 /
      data f(11, 9),f(11,10) / -1.16461018d+00 , -1.16461019d+00 /
      data f(11,11),f(11,12) / -1.16461034d+00 , -1.16461053d+00 /
      data f(11,13),f(11,14) / -1.16461032d+00 , -1.16461058d+00 /
      data f(11,15),f(11,16) / -1.16461029d+00 , -1.16461022d+00 /
      data f(11,17),f(11,18) / -1.16461104d+00 , -1.16461263d+00 /
      data f(11,19) /          -1.16461393d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 1.67088887d+00,-1.89434644d-03/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.06846265d+00,-5.51125317d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 4.40106213d-01, 1.70330331d-03/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 1.55100334d-01, 7.68505675d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.23325352d-02,-1.12158610d-04/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-5.50904273d-03,-2.08480635d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 4.89976240d-03,-1.29905105d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 8.18949815d-03,-4.36649611d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 1.27555801d-02,-3.80401047d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.60776654d-02,-3.56468200d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.56369336d-02,-3.95624152d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 9.64048474d-03,-5.46331191d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/-3.50907538d-03,-6.70793084d-05/
      data fpp( 1,14,1),fpp( 1,14,2)/-2.19243899d-02,-7.01184734d-06/
      data fpp( 1,15,1),fpp( 1,15,2)/-3.43130033d-02, 1.50609898d-04/
      data fpp( 1,16,1),fpp( 1,16,2)/-2.34597922d-02, 5.13151456d-04/
      data fpp( 1,17,1),fpp( 1,17,2)/ 1.02968390d-01,-5.99947228d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/ 5.53901960d-01,-3.53844365d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 1.50736425d+00,-7.27493018d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 8.98980056d-01,-7.96002926d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 5.75358269d-01,-2.14009247d-04/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 2.17929654d-01, 7.24638516d-04/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 5.77532129d-02, 4.14996983d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-1.85124895d-02,-2.52396249d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/-5.65755453d-03,-4.65197867d-05/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.14663521d-03,-1.44454040d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 5.39556369d-03,-3.12391971d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 8.99843979d-03,-2.50080074d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.13059893d-02,-1.95669732d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.08818127d-02,-2.06012999d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/ 6.61799053d-03,-2.77862271d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/-2.70556924d-03,-3.34901917d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/-1.63444201d-02,-3.96760623d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/-2.74893534d-02, 9.89116166d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/-2.57267757d-02, 3.48972540d-04/
      data fpp( 2,17,1),fpp( 2,17,2)/ 3.56442209d-02,-2.94459176d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/ 2.64723641d-01,-6.09749835d-04/
      data fpp( 2,19,1),fpp( 2,19,2)/ 7.48679267d-01,-1.01879508d-03/
      data fpp( 3, 1,1),fpp( 3, 1,2)/-2.11866936d-01,-3.38262082d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-9.33295608d-02,-5.99163355d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-7.79363473d-02, 4.57914824d-04/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-9.85823852d-02, 1.37547038d-04/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-3.73559867d-02,-7.10889769d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/-6.06337914d-03, 5.21546964d-06/
      data fpp( 3, 7,1),fpp( 3, 7,2)/-1.40662323d-03,-1.63681016d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 2.07720708d-03,-1.97650631d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 4.65402073d-03,-1.49302459d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/ 5.69653752d-03,-1.06171533d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/ 5.16781540d-03,-1.08181407d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.93219316d-03,-1.36504837d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/-1.91496767d-03,-1.48875244d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.06620895d-02, 5.66638130d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/-2.36556630d-02, 6.52305992d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/-4.39733452d-02, 2.75004222d-04/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.05321193d-01,-1.25337087d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/-2.56103642d-01,-3.05450875d-04/
      data fpp( 3,19,1),fpp( 3,19,2)/-5.66513794d-01,-5.48486812d-04/
      data fpp( 4, 1,1),fpp( 4, 1,2)/-1.36430312d-01,-2.54379324d-04/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.47475306d-01,-7.33750519d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.15889304d-01, 2.25874532d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-4.84583123d-02, 1.11265526d-04/
      data fpp( 4, 5,1),fpp( 4, 5,2)/-1.18770838d-02, 2.22379662d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/-5.56860890d-03, 5.48540949d-06/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-2.15942227d-03,-1.06828042d-05/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 5.94328002d-04,-1.20135927d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 2.36507728d-03,-8.21262494d-06/
      data fpp( 4,10,1),fpp( 4,10,2)/ 3.01930065d-03,-5.87670752d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/ 2.80988564d-03,-5.53814498d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.56795684d-03,-6.54211256d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/-1.07920010d-03,-5.95000477d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/-5.33802189d-03, 7.17133162d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/-9.35711440d-03, 3.81004783d-05/
      data fpp( 4,16,1),fpp( 4,16,2)/-6.29136348d-03, 1.59701155d-04/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.14577523d-02,-2.19990995d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 3.38404492d-02,-1.60171757d-04/
      data fpp( 4,19,1),fpp( 4,19,2)/ 6.60714291d-02,-3.10780871d-04/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 1.13800185d-01,-1.62194458d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 5.70374235d-02,-6.31039847d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 8.11708468d-03, 7.86685965d-05/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-8.48916584d-03, 7.10537986d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/-7.51071803d-03, 2.65966090d-05/
      data fpp( 5, 6,1),fpp( 5, 6,2)/-4.17954524d-03, 3.87856541d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/-1.84336767d-03,-6.59207063d-06/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 2.96409147d-05,-6.65168289d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.06423015d-03,-4.53579782d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/ 1.45913988d-03,-3.18372583d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/ 1.36368203d-03,-2.81789885d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/ 7.38939470d-04,-3.06287877d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/-5.50151941d-04,-2.17818607d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/-2.52710296d-03, 5.38202306d-06/
      data fpp( 5,15,1),fpp( 5,15,2)/-4.39027935d-03, 2.13786938d-05/
      data fpp( 5,16,1),fpp( 5,16,2)/-1.47136088d-03, 7.05524017d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 9.28242415d-03, 3.01586996d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/-4.12271414d-03,-7.39741999d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/-4.70242425d-02,-1.59895900d-04/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 1.01363724d-02,-8.06941407d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/ 9.37673164d-03,-3.33021187d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.22144562d-03, 3.05366153d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/-3.31662441d-03, 3.68118574d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/-3.64452407d-03, 1.72501549d-05/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-2.74473013d-03, 2.03692285d-06/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.28366705d-03,-3.65204634d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-1.83211661d-04,-3.62393751d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/ 4.12402118d-04,-2.47540364d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/ 6.80779825d-04,-1.69704795d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/ 6.49466237d-04,-1.46800456d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/ 3.61965277d-04,-1.41053382d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/-2.58992138d-04,-7.18660177d-07/
      data fpp( 6,14,1),fpp( 6,14,2)/-1.09676626d-03, 3.21897453d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/-1.44256820d-03, 1.14311621d-05/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.03336700d-03, 2.72677772d-05/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.06572711d-02, 2.50127292d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 2.19626474d-02,-3.02272941d-05/
      data fpp( 6,19,1),fpp( 6,19,2)/ 2.60356210d-02,-7.03287530d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 1.55515253d-02,-3.85046460d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 1.06745299d-02,-1.62615079d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 3.48705284d-03, 1.26338777d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-9.50496539d-04, 1.82463970d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/-2.06286571d-03, 9.34493437d-06/
      data fpp( 7, 6,1),fpp( 7, 6,2)/-1.51289423d-03, 1.06846554d-06/
      data fpp( 7, 7,1),fpp( 7, 7,2)/-7.16364135d-04,-1.90439654d-06/
      data fpp( 7, 8,1),fpp( 7, 8,2)/-1.45914271d-04,-2.02907938d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.85601378d-04,-1.23288594d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/ 2.77900818d-04,-9.75576856d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 2.58613021d-04,-7.41206634d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.18879421d-04,-6.86796607d-07/
      data fpp( 7,13,1),fpp( 7,13,2)/-1.36599506d-04,-2.06406939d-07/
      data fpp( 7,14,1),fpp( 7,14,2)/-4.48552012d-04, 1.74462436d-06/
      data fpp( 7,15,1),fpp( 7,15,2)/-4.51767855d-04, 5.65210949d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 8.44292881d-04, 1.00413377d-05/
      data fpp( 7,17,1),fpp( 7,17,2)/ 5.14033141d-03, 1.00323399d-05/
      data fpp( 7,18,1),fpp( 7,18,2)/ 1.15802045d-02,-1.14238971d-05/
      data fpp( 7,19,1),fpp( 7,19,2)/ 1.42063984d-02,-2.67523514d-05/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 1.73729780d-03,-8.09252645d-06/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 1.50546438d-03,-3.51194709d-06/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 4.94858677d-04, 2.04151482d-06/
      data fpp( 8, 4,1),fpp( 8, 4,2)/-3.14618180d-04, 4.18988780d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-5.68360831d-04, 2.38733398d-06/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-4.61852249d-04, 2.56976287d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-2.52454069d-04,-5.31639126d-07/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-9.07513571d-05,-5.37019782d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.78651927d-05,-4.15081744d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/ 3.00476332d-05,-2.28453241d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/ 2.68878183d-05,-1.93905292d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/-9.86090137d-06,-1.48525589d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/-7.09054140d-05,-2.13923500d-08/
      data fpp( 8,14,1),fpp( 8,14,2)/-1.17700835d-04, 4.27294989d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/-6.34923345d-05, 1.14001239d-06/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.24577856d-04, 1.19865544d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 8.59430222d-04,-2.08234157d-07/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.93078269d-03,-1.22071881d-06/
      data fpp( 8,19,1),fpp( 8,19,2)/ 2.79531426d-03,-2.34589059d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 6.24183474d-04,-1.71847116d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 4.05112530d-04,-7.62157675d-07/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 9.87524533d-05, 3.57701862d-07/
      data fpp( 9, 4,1),fpp( 9, 4,2)/-1.21430742d-04, 9.05150226d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/-1.74190966d-04, 5.89497234d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/-1.27436777d-04, 6.32608393d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/-7.15995870d-05,-1.17740591d-07/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-3.51603008d-05,-1.89898476d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/-6.68060719d-06,-8.02655047d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/-5.03135084d-06,-6.38395051d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/-9.28429429d-06,-5.13764750d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/-2.00358154d-05,-2.16545949d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/-3.79388384d-05, 1.31948548d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-5.96646488d-05, 9.00751758d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/-6.35228067d-05, 1.64104442d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/-1.77843041d-05, 7.43070568d-08/
      data fpp( 9,17,1),fpp( 9,17,2)/ 7.21477048d-05,-5.45332669d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.12364686d-04, 5.14236198d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 3.54845356d-05, 3.00638190d-07/
      data fpp(10, 1,1),fpp(10, 1,2)/ 2.84868305d-04,-3.55922172d-07/
      data fpp(10, 2,1),fpp(10, 2,2)/ 2.20225494d-04,-1.64255657d-07/
      data fpp(10, 3,1),fpp(10, 3,2)/ 5.60915093d-05, 1.87447993d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/-6.99588516d-05, 2.32076460d-07/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.06355307d-04, 1.50949361d-07/
      data fpp(10, 6,1),fpp(10, 6,2)/-8.37406447d-05,-6.07390566d-09/
      data fpp(10, 7,1),fpp(10, 7,2)/-4.67475826d-05,-2.46537388d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/-1.91674398d-05,-3.09111392d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/-5.99237856d-06,-2.27017045d-08/
      data fpp(10,10,1),fpp(10,10,2)/-1.10222983d-06,-7.88204298d-09/
      data fpp(10,11,1),fpp(10,11,2)/-3.55064114d-06,-1.71701236d-08/
      data fpp(10,12,1),fpp(10,12,2)/-1.10358369d-05, 9.62537233d-10/
      data fpp(10,13,1),fpp(10,13,2)/-2.21992323d-05, 7.31997466d-09/
      data fpp(10,14,1),fpp(10,14,2)/-3.23205702d-05, 1.89575642d-08/
      data fpp(10,15,1),fpp(10,15,2)/-2.87964387d-05, 2.54497687d-08/
      data fpp(10,16,1),fpp(10,16,2)/ 5.25936082d-06,-3.73566390d-08/
      data fpp(10,17,1),fpp(10,17,2)/ 7.52589590d-05,-1.75423213d-07/
      data fpp(10,18,1),fpp(10,18,2)/ 1.65678563d-04, 6.34494893d-08/
      data fpp(10,19,1),fpp(10,19,2)/ 2.22707593d-04, 1.94625255d-07/
      data fpp(11, 1,1),fpp(11, 1,2)/-7.41971653d-04,-5.23460655d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-5.61507747d-04,-2.59078689d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.88280755d-04, 1.79775413d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/ 9.67569258d-05, 7.97703892d-10/
      data fpp(11, 5,1),fpp(11, 5,2)/ 1.76166403d-04, 1.90316432d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/ 1.34475322d-04,-1.27242767d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/ 7.31575413d-05, 1.02654635d-08/
      data fpp(11, 8,1),fpp(11, 8,2)/ 3.06374699d-05,-5.53757732d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/ 1.07439281d-07, 5.88484578d-09/
      data fpp(11,10,1),fpp(11,10,2)/-1.07976351d-05,-3.00180581d-09/
      data fpp(11,11,1),fpp(11,11,2)/-1.09409294d-05,-2.27762253d-09/
      data fpp(11,12,1),fpp(11,12,2)/-1.66458154d-06, 9.71229595d-09/
      data fpp(11,13,1),fpp(11,13,2)/ 1.70321162d-05,-1.25715613d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 3.61790351d-05, 1.23739491d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 2.50157193d-05,-3.92423496d-09/
      data fpp(11,16,1),fpp(11,16,2)/-6.80859304d-05,-9.87700918d-09/
      data fpp(11,17,1),fpp(11,17,2)/-2.74195730d-04,-9.96772830d-09/
      data fpp(11,18,1),fpp(11,18,2)/-5.69973031d-04, 3.54792237d-09/
      data fpp(11,19,1),fpp(11,19,2)/-7.29020046d-04, 1.31760388d-08/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-3.66607890d-03,-8.31363547d-04/
      data fpppp( 1, 3),fpppp( 1, 4)/ 5.43572073d-03,-3.10486220d-04/
      data fpppp( 1, 5),fpppp( 1, 6)/ 2.86060480d-03,-6.76551308d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/ 6.07191897d-05, 6.53038712d-06/
      data fpppp( 1, 9),fpppp( 1,10)/-1.02599667d-05,-4.01303221d-05/
      data fpppp( 1,11),fpppp( 1,12)/-5.49877639d-05,-7.32616521d-05/
      data fpppp( 1,13),fpppp( 1,14)/-8.11523011d-05, 8.19255907d-05/
      data fpppp( 1,15),fpppp( 1,16)/ 1.15052009d-04, 8.52375842d-04/
      data fpppp( 1,17),fpppp( 1,18)/ 3.40994286d-03, 4.97817603d-03/
      data fpppp( 1,19) /             6.82907602d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-2.48903093d-03,-6.57323040d-04/
      data fpppp( 2, 3),fpppp( 2, 4)/ 3.08991340d-03, 1.32799908d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.41353126d-03,-4.39686685d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 4.21707682d-05,-2.31206354d-06/
      data fpppp( 2, 9),fpppp( 2,10)/-1.16856571d-05,-2.86649046d-05/
      data fpppp( 2,11),fpppp( 2,12)/-3.75582864d-05,-5.14806902d-05/
      data fpppp( 2,13),fpppp( 2,14)/-6.01032057d-05, 3.29760443d-05/
      data fpppp( 2,15),fpppp( 2,16)/ 7.78340866d-05, 4.30138271d-04/
      data fpppp( 2,17),fpppp( 2,18)/ 1.77811796d-03, 2.51989528d-03/
      data fpppp( 2,19) /             3.43487329d-03 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.75226149d-03,-9.31900085d-04/
      data fpppp( 3, 3),fpppp( 3, 4)/-7.08787875d-04, 1.60469650d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-7.97651952d-04,-2.10116153d-04/
      data fpppp( 3, 7),fpppp( 3, 8)/ 3.99654670d-05,-2.01212507d-05/
      data fpppp( 3, 9),fpppp( 3,10)/-1.39014633d-05,-1.63307085d-05/
      data fpppp( 3,11),fpppp( 3,12)/-1.50500364d-05,-2.58831538d-05/
      data fpppp( 3,13),fpppp( 3,14)/-3.81096631d-05,-5.56758541d-05/
      data fpppp( 3,15),fpppp( 3,16)/ 6.02597586d-06,-4.07874566d-04/
      data fpppp( 3,17),fpppp( 3,18)/-8.36337665d-04,-1.61285084d-03/
      data fpppp( 3,19) /            -2.28992113d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 5.25482811d-04, 3.63684301d-04/
      data fpppp( 4, 3),fpppp( 4, 4)/ 5.77639670d-04,-5.23543537d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-3.34451343d-04, 4.49837003d-05/
      data fpppp( 4, 7),fpppp( 4, 8)/-1.94407554d-05,-6.54686007d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-1.33518643d-05,-7.03723699d-06/
      data fpppp( 4,11),fpppp( 4,12)/-1.03174906d-05,-1.36436281d-05/
      data fpppp( 4,13),fpppp( 4,14)/-1.94216855d-05,-5.36952062d-06/
      data fpppp( 4,15),fpppp( 4,16)/ 5.52835243d-05, 2.09326030d-04/
      data fpppp( 4,17),fpppp( 4,18)/-1.15857534d-05, 1.15031850d-04/
      data fpppp( 4,19) /             1.42355342d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.50419160d-04, 4.66776848d-05/
      data fpppp( 5, 3),fpppp( 5, 4)/ 4.34253778d-04, 1.55152502d-04/
      data fpppp( 5, 5),fpppp( 5, 6)/ 2.18114801d-07,-1.48614623d-05/
      data fpppp( 5, 7),fpppp( 5, 8)/-4.71978546d-07,-1.10407628d-05/
      data fpppp( 5, 9),fpppp( 5,10)/-5.67013118d-06,-4.65948283d-06/
      data fpppp( 5,11),fpppp( 5,12)/-5.11399236d-06,-6.64163031d-06/
      data fpppp( 5,13),fpppp( 5,14)/-8.18041742d-06,-1.90827679d-06/
      data fpppp( 5,15),fpppp( 5,16)/ 2.26400029d-05, 1.98273957d-04/
      data fpppp( 5,17),fpppp( 5,18)/-3.45643835d-04,-2.65234015d-04/
      data fpppp( 5,19) /            -3.63203512d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/-1.80235357d-04,-8.16558665d-05/
      data fpppp( 6, 3),fpppp( 6, 4)/ 6.31201095d-05, 4.62083880d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/ 4.65656053d-06, 8.82698536d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/-6.28835289d-06,-5.31003556d-06/
      data fpppp( 6, 9),fpppp( 6,10)/-2.76200142d-06,-3.27612305d-06/
      data fpppp( 6,11),fpppp( 6,12)/-2.11498411d-06,-3.63518284d-06/
      data fpppp( 6,13),fpppp( 6,14)/-3.35167190d-06, 4.03286835d-06/
      data fpppp( 6,15),fpppp( 6,16)/ 1.67385290d-05, 9.83172442d-05/
      data fpppp( 6,17),fpppp( 6,18)/ 1.88706288d-05,-7.29114286d-05/
      data fpppp( 6,19) /            -1.61169074d-04 /
      data fpppp( 7, 1),fpppp( 7, 2)/-7.20991462d-05,-2.63243468d-05/
      data fpppp( 7, 3),fpppp( 7, 4)/ 3.87676316d-05, 3.62494833d-05/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.57452473d-05, 5.09966761d-07/
      data fpppp( 7, 7),fpppp( 7, 8)/-2.99159785d-06,-2.10838905d-06/
      data fpppp( 7, 9),fpppp( 7,10)/-2.91089891d-06,-6.00987819d-07/
      data fpppp( 7,11),fpppp( 7,12)/-1.38038404d-06,-1.10422420d-06/
      data fpppp( 7,13),fpppp( 7,14)/-1.14743874d-06, 2.30556435d-06/
      data fpppp( 7,15),fpppp( 7,16)/ 1.04493812d-05, 3.38535057d-05/
      data fpppp( 7,17),fpppp( 7,18)/ 3.41352634d-05,-4.17644832d-05/
      data fpppp( 7,19) /            -9.58980858d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/-1.73861008d-05,-8.18899009d-06/
      data fpppp( 8, 3),fpppp( 8, 4)/ 3.41572360d-06, 6.59382669d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.55302201d-06, 8.09159226d-07/
      data fpppp( 8, 7),fpppp( 8, 8)/-6.16283055d-07,-1.20575503d-06/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.10310297d-07,-7.33886471d-07/
      data fpppp( 8,11),fpppp( 8,12)/-2.39122856d-07,-3.24956391d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 8.12008451d-08, 8.55098527d-07/
      data fpppp( 8,15),fpppp( 8,16)/ 2.55864030d-06, 2.94204166d-06/
      data fpppp( 8,17),fpppp( 8,18)/ 6.48012363d-06,-2.67252976d-06/
      data fpppp( 8,19) /            -8.19925877d-06 /
      data fpppp( 9, 1),fpppp( 9, 2)/-2.58483651d-06,-9.18321292d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/ 1.02077366d-06, 2.00583954d-06/
      data fpppp( 9, 5),fpppp( 9, 6)/ 1.00124651d-06,-3.99608171d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/-2.96423200d-07, 6.17794149d-08/
      data fpppp( 9, 9),fpppp( 9,10)/-4.28270019d-07, 4.14744274d-08/
      data fpppp( 9,11),fpppp( 9,12)/-9.17596788d-08,-6.43503731d-08/
      data fpppp( 9,13),fpppp( 9,14)/-7.99289401d-08, 1.54698892d-07/
      data fpppp( 9,15),fpppp( 9,16)/ 5.33192514d-07, 6.88330686d-07/
      data fpppp( 9,17),fpppp( 9,18)/-6.34904880d-07,-1.13161285d-06/
      data fpppp( 9,19) /            -1.86447156d-06 /
      data fpppp(10, 1),fpppp(10, 2)/-2.34359100d-06,-1.04904913d-06/
      data fpppp(10, 3),fpppp(10, 4)/ 5.70317139d-07, 1.05279800d-06/
      data fpppp(10, 5),fpppp(10, 6)/ 5.97725198d-07, 9.69682594d-08/
      data fpppp(10, 7),fpppp(10, 8)/-1.22894243d-07,-1.70166447d-07/
      data fpppp(10, 9),fpppp(10,10)/-6.07448564d-08,-8.39488804d-08/
      data fpppp(10,11),fpppp(10,12)/-4.37732244d-08,-4.31652898d-08/
      data fpppp(10,13),fpppp(10,14)/-4.25759410d-09, 1.22719114d-07/
      data fpppp(10,15),fpppp(10,16)/ 3.32109308d-07, 3.80743728d-07/
      data fpppp(10,17),fpppp(10,18)/ 3.01543705d-07,-3.61718210d-07/
      data fpppp(10,19) /            -8.58105298d-07 /
      data fpppp(11, 1),fpppp(11, 2)/ 4.69090845d-06, 2.02012367d-06/
      data fpppp(11, 3),fpppp(11, 4)/-1.20561795d-06,-2.48901060d-06/
      data fpppp(11, 5),fpppp(11, 6)/-1.17603181d-06,-7.28956913d-08/
      data fpppp(11, 7),fpppp(11, 8)/ 2.90012576d-07, 4.07079661d-08/
      data fpppp(11, 9),fpppp(11,10)/ 2.66558003d-07, 7.05573970d-08/
      data fpppp(11,11),fpppp(11,12)/ 9.69192099d-08, 1.06944298d-07/
      data fpppp(11,13),fpppp(11,14)/ 4.05245888d-08,-2.42029377d-07/
      data fpppp(11,15),fpppp(11,16)/-8.91021166d-07,-1.11018599d-06/
      data fpppp(11,17),fpppp(11,18)/-1.44872382d-06, 1.52503109d-06/
      data fpppp(11,19) /             3.55241665d-06 /
 
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
      subroutine c1_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  1.33799300d-02 /
      data f( 1, 3),f( 1, 4) /  9.41929000d-03 ,  2.05173400d-02 /
      data f( 1, 5),f( 1, 6) /  5.19237500d-02 ,  5.39359500d-02 /
      data f( 1, 7),f( 1, 8) /  3.02011300d-02 ,  1.14493700d-02 /
      data f( 1, 9),f( 1,10) /  8.76820000d-04 , -2.84298000d-03 /
      data f( 1,11),f( 1,12) / -8.78740000d-04 ,  5.72966000d-03 /
      data f( 1,13),f( 1,14) /  1.53040700d-02 ,  2.44521800d-02 /
      data f( 1,15),f( 1,16) /  3.04592900d-02 ,  3.25292000d-02 /
      data f( 1,17),f( 1,18) /  2.64671000d-02 ,  4.39266000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  7.91092000d-03 /
      data f( 2, 3),f( 2, 4) /  9.58260000d-03 ,  2.39829800d-02 /
      data f( 2, 5),f( 2, 6) /  4.80521500d-02 ,  3.84014300d-02 /
      data f( 2, 7),f( 2, 8) /  1.93596500d-02 ,  6.84353000d-03 /
      data f( 2, 9),f( 2,10) /  4.78110000d-04 , -1.52891000d-03 /
      data f( 2,11),f( 2,12) / -4.05150000d-04 ,  3.25422000d-03 /
      data f( 2,13),f( 2,14) /  9.13188000d-03 ,  1.63589200d-02 /
      data f( 2,15),f( 2,16) /  2.30633600d-02 ,  2.74991200d-02 /
      data f( 2,17),f( 2,18) /  2.58683600d-02 ,  9.37462000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  1.51600700d-02 /
      data f( 3, 3),f( 3, 4) /  2.40974300d-02 ,  3.51065000d-02 /
      data f( 3, 5),f( 3, 6) /  3.76622400d-02 ,  2.44240200d-02 /
      data f( 3, 7),f( 3, 8) /  1.15635800d-02 ,  3.91648000d-03 /
      data f( 3, 9),f( 3,10) /  2.74280000d-04 , -8.08680000d-04 /
      data f( 3,11),f( 3,12) / -1.87740000d-04 ,  1.81555000d-03 /
      data f( 3,13),f( 3,14) /  5.17986000d-03 ,  9.74417000d-03 /
      data f( 3,15),f( 3,16) /  1.47015400d-02 ,  1.84957300d-02 /
      data f( 3,17),f( 3,18) /  1.88633800d-02 ,  1.11920400d-02 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  1.90407500d-02 /
      data f( 4, 3),f( 4, 4) /  3.03041100d-02 ,  3.06220000d-02 /
      data f( 4, 5),f( 4, 6) /  2.42072400d-02 ,  1.42092000d-02 /
      data f( 4, 7),f( 4, 8) /  6.48923000d-03 ,  2.14569000d-03 /
      data f( 4, 9),f( 4,10) /  1.49640000d-04 , -4.18620000d-04 /
      data f( 4,11),f( 4,12) / -8.76900000d-05 ,  9.79590000d-04 /
      data f( 4,13),f( 4,14) /  2.81543000d-03 ,  5.40076000d-03 /
      data f( 4,15),f( 4,16) /  8.33764000d-03 ,  1.05026800d-02 /
      data f( 4,17),f( 4,18) /  9.72187000d-03 ,  4.14620000d-03 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  8.56058000d-03 /
      data f( 5, 3),f( 5, 4) /  1.81523000d-02 ,  1.85725800d-02 /
      data f( 5, 5),f( 5, 6) /  1.36171500d-02 ,  7.67911000d-03 /
      data f( 5, 7),f( 5, 8) /  3.46588000d-03 ,  1.14366000d-03 /
      data f( 5, 9),f( 5,10) /  9.00600000d-05 , -2.10190000d-04 /
      data f( 5,11),f( 5,12) / -4.04200000d-05 ,  5.14780000d-04 /
      data f( 5,13),f( 5,14) /  1.47795000d-03 ,  2.85278000d-03 /
      data f( 5,15),f( 5,16) /  4.42662000d-03 ,  5.51786000d-03 /
      data f( 5,17),f( 5,18) /  4.93350000d-03 ,  2.06956000d-03 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  3.77933000d-03 /
      data f( 6, 3),f( 6, 4) /  8.97707000d-03 ,  9.62623000d-03 /
      data f( 6, 5),f( 6, 6) /  6.99806000d-03 ,  3.93615000d-03 /
      data f( 6, 7),f( 6, 8) /  1.79498000d-03 ,  6.06830000d-04 /
      data f( 6, 9),f( 6,10) /  5.99300000d-05 , -1.01110000d-04 /
      data f( 6,11),f( 6,12) / -1.73800000d-05 ,  2.64640000d-04 /
      data f( 6,13),f( 6,14) /  7.57510000d-04 ,  1.45951000d-03 /
      data f( 6,15),f( 6,16) /  2.25060000d-03 ,  2.76588000d-03 /
      data f( 6,17),f( 6,18) /  2.41709000d-03 ,  1.00587000d-03 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.65394000d-03 /
      data f( 7, 3),f( 7, 4) /  4.15357000d-03 ,  4.62190000d-03 /
      data f( 7, 5),f( 7, 6) /  3.41498000d-03 ,  1.95318000d-03 /
      data f( 7, 7),f( 7, 8) /  9.08390000d-04 ,  3.16580000d-04 /
      data f( 7, 9),f( 7,10) /  3.41400000d-05 , -4.51600000d-05 /
      data f( 7,11),f( 7,12) / -5.94000000d-06 ,  1.34290000d-04 /
      data f( 7,13),f( 7,14) /  3.80840000d-04 ,  7.30910000d-04 /
      data f( 7,15),f( 7,16) /  1.11781000d-03 ,  1.34989000d-03 /
      data f( 7,17),f( 7,18) /  1.14525000d-03 ,  4.59370000d-04 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  3.26500000d-04 /
      data f( 8, 3),f( 8, 4) /  8.53550000d-04 ,  9.93680000d-04 /
      data f( 8, 5),f( 8, 6) /  7.72210000d-04 ,  4.68910000d-04 /
      data f( 8, 7),f( 8, 8) /  2.33120000d-04 ,  9.11900000d-05 /
      data f( 8, 9),f( 8,10) /  2.15900000d-05 , -4.77000000d-06 /
      data f( 8,11),f( 8,12) /  1.42000000d-06 ,  3.35100000d-05 /
      data f( 8,13),f( 8,14) /  9.20900000d-05 ,  1.74920000d-04 /
      data f( 8,15),f( 8,16) /  2.63020000d-04 ,  3.08450000d-04 /
      data f( 8,17),f( 8,18) /  2.48360000d-04 ,  9.33200000d-05 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  6.76800000d-05 /
      data f( 9, 3),f( 9, 4) /  1.76410000d-04 ,  2.10980000d-04 /
      data f( 9, 5),f( 9, 6) /  1.69020000d-04 ,  1.08140000d-04 /
      data f( 9, 7),f( 9, 8) /  5.90000000d-05 ,  2.83700000d-05 /
      data f( 9, 9),f( 9,10) /  8.49000000d-06 ,  9.20000000d-07 /
      data f( 9,11),f( 9,12) /  1.51000000d-06 ,  7.53000000d-06 /
      data f( 9,13),f( 9,14) /  2.03000000d-05 ,  4.03200000d-05 /
      data f( 9,15),f( 9,16) /  5.99600000d-05 ,  6.67200000d-05 /
      data f( 9,17),f( 9,18) /  5.43500000d-05 ,  1.89400000d-05 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  1.41500000d-05 /
      data f(10, 3),f(10, 4) /  3.80400000d-05 ,  4.33000000d-05 /
      data f(10, 5),f(10, 6) /  3.59600000d-05 ,  2.58900000d-05 /
      data f(10, 7),f(10, 8) /  1.45700000d-05 ,  6.92000000d-06 /
      data f(10, 9),f(10,10) /  2.76000000d-06 ,  8.20000000d-07 /
      data f(10,11),f(10,12) /  6.20000000d-07 ,  1.41000000d-06 /
      data f(10,13),f(10,14) /  3.77000000d-06 ,  7.56000000d-06 /
      data f(10,15),f(10,16) /  1.13300000d-05 ,  1.29200000d-05 /
      data f(10,17),f(10,18) /  1.00300000d-05 ,  3.39000000d-06 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 ,  9.70000000d-07 /
      data f(11, 3),f(11, 4) /  7.60000000d-07 ,  1.34000000d-06 /
      data f(11, 5),f(11, 6) /  1.28000000d-06 ,  1.75000000d-06 /
      data f(11, 7),f(11, 8) /  8.20000000d-07 ,  9.00000000d-08 /
      data f(11, 9),f(11,10) / -3.00000000d-08 ,  3.40000000d-07 /
      data f(11,11),f(11,12) /  1.50000000d-07 ,  4.00000000d-08 /
      data f(11,13),f(11,14) /  5.60000000d-07 ,  2.70000000d-07 /
      data f(11,15),f(11,16) /  4.70000000d-07 ,  4.70000000d-07 /
      data f(11,17),f(11,18) /  4.10000000d-07 ,  3.50000000d-07 /
      data f(11,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-4.93563378d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.16032693d-01,-1.81075544d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 1.47574876d-01, 1.77431354d-04/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 1.18729850d-01, 3.74871529d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-4.04882407d-02,-4.58415871d-04/
      data fpp( 1, 6,1),fpp( 1, 6,2)/-2.24929019d-03,-3.04860644d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 1.35316951d-02, 1.33037247d-04/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 8.78592009d-03, 7.16952551d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 1.22449169d-03, 7.09343324d-05/
      data fpp( 1,10,1),fpp( 1,10,2)/-3.41352086d-03, 5.57324154d-05/
      data fpp( 1,11,1),fpp( 1,11,2)/-1.56912182d-03, 4.71784059d-05/
      data fpp( 1,12,1),fpp( 1,12,2)/ 5.85244293d-03, 3.42035611d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.14044294d-02,-6.03205040d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.94440909d-03,-3.56533596d-05/
      data fpp( 1,15,1),fpp( 1,15,2)/-1.53547869d-02,-3.98145114d-05/
      data fpp( 1,16,1),fpp( 1,16,2)/-3.54539489d-02,-4.13205949d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/-4.34640205d-02,-2.82823709d-04/
      data fpp( 1,18,1),fpp( 1,18,2)/ 6.06774251d-03, 2.11875031d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 4.96230384d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-2.49468653d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 4.92455737d-02,-6.75986949d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 5.83458470d-02, 1.45509032d-04/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 4.05620597d-02, 2.49284567d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-2.48689987d-02,-5.62519899d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 5.54006038d-03,-2.23983722d-05/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.20710498d-02, 8.86497874d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 6.75387983d-03, 5.93388225d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 8.15096620d-04, 4.30369228d-05/
      data fpp( 2,10,1),fpp( 2,10,2)/-2.40839828d-03, 3.00174865d-05/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.04647635d-03, 2.47399312d-05/
      data fpp( 2,12,1),fpp( 2,12,2)/ 4.20883413d-03, 2.31593888d-05/
      data fpp( 2,13,1),fpp( 2,13,2)/ 8.89382125d-03, 1.57199135d-05/
      data fpp( 2,14,1),fpp( 2,14,2)/ 5.51066183d-03,-5.07624269d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/-4.59158621d-03,-2.67709427d-05/
      data fpp( 2,16,1),fpp( 2,16,2)/-1.66410222d-02,-2.39607865d-05/
      data fpp( 2,17,1),fpp( 2,17,2)/-2.41045189d-02,-2.41377111d-04/
      data fpp( 2,18,1),fpp( 2,18,2)/-4.52020502d-03, 9.76904318d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 2.77762584d-04/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.37794034d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/-7.77914785d-03,-7.69816311d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.65217846d-02, 7.23579589d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-9.71889691d-02,-8.81476044d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.64752046d-02,-2.26967341d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 1.74596887d-02, 4.83793689d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.12739457d-02, 5.61166655d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 4.48952061d-03, 3.99543692d-05/
      data fpp( 3, 9,1),fpp( 3, 9,2)/ 1.92241831d-04, 2.43598577d-05/
      data fpp( 3,10,1),fpp( 3,10,2)/-1.20504601d-03, 1.61606001d-05/
      data fpp( 3,11,1),fpp( 3,11,2)/-3.93292761d-04, 1.32317418d-05/
      data fpp( 3,12,1),fpp( 3,12,2)/ 2.19470054d-03, 1.38534327d-05/
      data fpp( 3,13,1),fpp( 3,13,2)/ 6.30436561d-03, 1.30157273d-05/
      data fpp( 3,14,1),fpp( 3,14,2)/ 1.04971836d-02, 6.08365812d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 1.05397717d-02,-1.37667598d-05/
      data fpp( 3,16,1),fpp( 3,16,2)/ 6.65859778d-03,-2.08074190d-05/
      data fpp( 3,17,1),fpp( 3,17,2)/-1.38676637d-02,-1.08595964d-04/
      data fpp( 3,18,1),fpp( 3,18,2)/-6.39358824d-02,-2.71481245d-05/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 5.94646227d-06/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-4.93848724d-05/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-9.89722623d-02,-7.11903552d-05/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.11654309d-01,-1.32497107d-04/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-2.63986635d-02,-5.55494179d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.72076571d-02,-4.92642218d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 1.49233449d-02, 3.76095052d-05/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 8.15444734d-03, 3.55104009d-05/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 3.03827775d-03, 2.29346910d-05/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 3.16496057d-04, 1.36002349d-05/
      data fpp( 4,10,1),fpp( 4,10,2)/-6.95497662d-04, 8.33176926d-06/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.96992602d-04, 7.02408805d-06/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.47740371d-03, 7.75287854d-06/
      data fpp( 4,13,1),fpp( 4,13,2)/ 3.99087631d-03, 8.07799778d-06/
      data fpp( 4,14,1),fpp( 4,14,2)/ 7.01276375d-03, 4.90453035d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/ 1.03825793d-02,-6.60311918d-06/
      data fpp( 4,16,1),fpp( 4,16,2)/ 1.42547911d-02,-2.48024536d-05/
      data fpp( 4,17,1),fpp( 4,17,2)/ 2.82984538d-02,-7.09380663d-05/
      data fpp( 4,18,1),fpp( 4,18,2)/ 4.75454947d-02, 2.08631190d-05/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 7.32537905d-05/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 1.06632607d-04/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 5.90077969d-02, 2.17205855d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 4.25352591d-02,-1.31646549d-04/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.12255432d-02,-4.54207885d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 1.64024163d-02,-9.21289681d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 1.12804515d-02, 2.33157757d-05/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 5.33226494d-03, 1.94383939d-05/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.80760840d-03, 1.23912487d-05/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 1.03213941d-04, 7.11381130d-06/
      data fpp( 5,10,1),fpp( 5,10,2)/-3.72083340d-04, 4.35450609d-06/
      data fpp( 5,11,1),fpp( 5,11,2)/-8.54568320d-05, 3.66936435d-06/
      data fpp( 5,12,1),fpp( 5,12,2)/ 8.03284609d-04, 4.09383651d-06/
      data fpp( 5,13,1),fpp( 5,13,2)/ 2.37892916d-03, 4.43348962d-06/
      data fpp( 5,14,1),fpp( 5,14,2)/ 4.54208140d-03, 2.87180500d-06/
      data fpp( 5,15,1),fpp( 5,15,2)/ 6.79903103d-03,-3.98010963d-06/
      data fpp( 5,16,1),fpp( 5,16,2)/ 8.51975774d-03,-1.59073665d-05/
      data fpp( 5,17,1),fpp( 5,17,2)/ 5.14920864d-03,-3.29264244d-05/
      data fpp( 5,18,1),fpp( 5,18,2)/-6.98529634d-03, 1.08382641d-05/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 3.72361679d-05/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 7.10627798d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-2.84845444d-04, 1.97665404d-05/
      data fpp( 6, 3,1),fpp( 6, 3,2)/ 1.29511923d-02,-6.50243415d-05/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.59701709d-02,-3.25839746d-05/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 1.24866779d-02,-1.27956016d-06/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 6.84596886d-03, 1.16778152d-05/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 2.97529290d-03, 9.81269921d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/ 8.96088643d-04, 6.25258792d-06/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-2.25518226d-05, 3.65194911d-06/
      data fpp( 6,10,1),fpp( 6,10,2)/-2.00568980d-04, 2.29121566d-06/
      data fpp( 6,11,1),fpp( 6,11,2)/-4.27000704d-05, 1.86938826d-06/
      data fpp( 6,12,1),fpp( 6,12,2)/ 4.61537851d-04, 2.12863128d-06/
      data fpp( 6,13,1),fpp( 6,13,2)/ 1.30236706d-03, 2.26708660d-06/
      data fpp( 6,14,1),fpp( 6,14,2)/ 2.53195065d-03, 1.35082232d-06/
      data fpp( 6,15,1),fpp( 6,15,2)/ 4.06129657d-03,-2.32497587d-06/
      data fpp( 6,16,1),fpp( 6,16,2)/ 5.25433790d-03,-8.59951884d-06/
      data fpp( 6,17,1),fpp( 6,17,2)/ 5.63175165d-03,-1.51211488d-05/
      data fpp( 6,18,1),fpp( 6,18,2)/ 4.70649069d-03, 5.33831394d-06/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 1.80888930d-05/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 3.60224668d-05/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 5.87222485d-03, 1.08655664d-05/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 1.01014916d-02,-2.87433324d-05/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 9.50225335d-03,-1.77702369d-05/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 6.51511218d-03,-6.90720096d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 3.57543302d-03, 5.24031726d-06/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 1.59000346d-03, 4.75005104d-06/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 5.25957027d-04, 2.93827858d-06/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 9.11533490d-05, 2.05903465d-06/
      data fpp( 7,10,1),fpp( 7,10,2)/-1.00760741d-04, 1.01398283d-06/
      data fpp( 7,11,1),fpp( 7,11,2)/-2.21428866d-05, 9.96234039d-07/
      data fpp( 7,12,1),fpp( 7,12,2)/ 2.25523989d-04, 1.06168102d-06/
      data fpp( 7,13,1),fpp( 7,13,2)/ 6.62082586d-04, 1.13624190d-06/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.28219600d-03, 6.04551394d-07/
      data fpp( 7,15,1),fpp( 7,15,2)/ 1.99330270d-03,-1.34464747d-06/
      data fpp( 7,16,1),fpp( 7,16,2)/ 2.52665064d-03,-4.51516150d-06/
      data fpp( 7,17,1),fpp( 7,17,2)/ 2.19346474d-03,-6.79790652d-06/
      data fpp( 7,18,1),fpp( 7,18,2)/ 5.71893599d-04, 2.83238758d-06/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 9.05895621d-06/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 7.65744496d-06/
      data fpp( 8, 2,1),fpp( 8, 2,2)/ 6.57881741d-05, 2.45101008d-06/
      data fpp( 8, 3,1),fpp( 8, 3,2)/ 1.30180897d-03,-5.42848526d-06/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 1.79079453d-03,-3.95226902d-06/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 1.35166452d-03,-4.58438666d-07/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 7.40736516d-04, 8.76223683d-07/
      data fpp( 8, 7,1),fpp( 8, 7,2)/ 3.29803157d-04, 1.00414393d-06/
      data fpp( 8, 8,1),fpp( 8, 8,2)/ 1.04744598d-04, 7.38800579d-07/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-2.80041356d-05, 3.80453750d-07/
      data fpp( 8,10,1),fpp( 8,10,2)/-2.64932877d-05, 3.33784422d-07/
      data fpp( 8,11,1),fpp( 8,11,2)/-5.34130499d-06, 2.37408563d-07/
      data fpp( 8,12,1),fpp( 8,12,2)/ 5.21791087d-05, 2.70581327d-07/
      data fpp( 8,13,1),fpp( 8,13,2)/ 1.50108709d-04, 2.69666130d-07/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.94696668d-04, 1.05754153d-07/
      data fpp( 8,15,1),fpp( 8,15,2)/ 4.54183608d-04,-3.76482743d-07/
      data fpp( 8,16,1),fpp( 8,16,2)/ 5.36119124d-04,-1.16002318d-06/
      data fpp( 8,17,1),fpp( 8,17,2)/ 4.84469965d-04,-1.31462453d-06/
      data fpp( 8,18,1),fpp( 8,18,2)/ 2.92773861d-04, 7.21521294d-07/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00, 2.13173935d-06/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00, 1.52313650d-06/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 2.76342454d-04, 4.89427002d-07/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 4.28552509d-04,-1.01784451d-06/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 4.07688545d-04,-8.67648974d-07/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 3.15709731d-04,-1.03359596d-07/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 2.02620917d-04, 1.45887358d-07/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 9.76839059d-05, 2.24210164d-07/
      data fpp( 9, 8,1),fpp( 9, 8,2)/ 3.04845793d-05, 6.78719867d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.75631936d-05, 1.49301889d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.46610851d-06, 7.35204555d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/-1.11893423d-07, 4.62162887d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.45595767d-05, 6.74143898d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/ 3.92425765d-05, 8.91261522d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/ 6.73573241d-05, 1.10810014d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/ 1.00342866d-04,-1.56250158d-07/
      data fpp( 9,16,1),fpp( 9,16,2)/ 1.27132862d-04,-1.58880370d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/ 8.59354039d-05,-3.56028362d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/ 7.03095765d-06, 2.00593818d-07/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 5.41853091d-07/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 3.67554314d-07/
      data fpp(10, 2,1),fpp(10, 2,2)/ 6.05820092d-05, 1.24491373d-07/
      data fpp(10, 3,1),fpp(10, 3,2)/ 2.16600998d-04,-2.81119804d-07/
      data fpp(10, 4,1),fpp(10, 4,2)/ 2.68571291d-04,-1.17812155d-07/
      data fpp(10, 5,1),fpp(10, 5,2)/ 2.06276554d-04,-3.63157390d-09/
      data fpp(10, 6,1),fpp(10, 6,2)/ 1.19899817d-04,-3.14615489d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/ 5.76012188d-05, 5.44777697d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 2.15370841d-05, 3.37504703d-08/
      data fpp(10, 9,1),fpp(10, 9,2)/ 1.97136128d-06, 1.99203492d-08/
      data fpp(10,10,1),fpp(10,10,2)/-2.38227830d-06, 1.97681329d-08/
      data fpp(10,11,1),fpp(10,11,2)/-9.11213153d-08, 5.40711905d-09/
      data fpp(10,12,1),fpp(10,12,2)/ 8.74258467d-06, 1.80033909d-08/
      data fpp(10,13,1),fpp(10,13,2)/ 2.44809847d-05, 1.67793175d-08/
      data fpp(10,14,1),fpp(10,14,2)/ 4.69140352d-05, 6.79339050d-10/
      data fpp(10,15,1),fpp(10,15,2)/ 7.10249267d-05,-2.06966737d-08/
      data fpp(10,16,1),fpp(10,16,2)/ 8.29294276d-05,-4.86926442d-08/
      data fpp(10,17,1),fpp(10,17,2)/ 6.99284192d-05,-5.33327497d-08/
      data fpp(10,18,1),fpp(10,18,2)/ 3.20823085d-05, 3.70236428d-08/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 1.00238179d-07/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-2.96790007d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/-1.79097255d-04,-1.54419985d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-5.04889249d-04, 2.06469949d-08/
      data fpp(11, 4,1),fpp(11, 4,2)/-5.69458145d-04,-1.97459811d-08/
      data fpp(11, 5,1),fpp(11, 5,2)/-4.29524527d-04, 1.99369295d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-2.50469908d-04,-2.82017368d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-1.08980609d-04, 8.87001789d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/-2.57485421d-05, 4.72166528d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.69068064d-06, 8.84332099d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 7.45988915d-06,-1.06949492d-08/
      data fpp(11,11,1),fpp(11,11,2)/ 2.29431066d-06, 3.36475945d-10/
      data fpp(11,12,1),fpp(11,12,2)/-1.72025423d-05, 1.41490455d-08/
      data fpp(11,13,1),fpp(11,13,2)/-4.82892424d-05,-1.91326578d-08/
      data fpp(11,14,1),fpp(11,14,2)/-8.70757676d-05, 1.37815856d-08/
      data fpp(11,15,1),fpp(11,15,2)/-1.33646213d-04,-6.59368453d-09/
      data fpp(11,16,1),fpp(11,16,2)/-1.69629714d-04, 5.93152555d-10/
      data fpp(11,17,1),fpp(11,17,2)/-1.34222960d-04, 6.21074312d-10/
      data fpp(11,18,1),fpp(11,18,2)/-5.76724042d-05,-3.07744980d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-5.71127510d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/-9.84720723d-04,-1.04733985d-03/
      data fpppp( 1, 3),fpppp( 1, 4)/ 1.04649547d-04,-2.99449091d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/ 4.05093024d-03,-1.36180757d-03/
      data fpppp( 1, 7),fpppp( 1, 8)/ 4.88221448d-05,-6.50866237d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.25851470d-05, 7.01509865d-05/
      data fpppp( 1,11),fpppp( 1,12)/ 6.57556020d-05, 1.45654857d-06/
      data fpppp( 1,13),fpppp( 1,14)/-1.83756495d-04,-1.07150970d-04/
      data fpppp( 1,15),fpppp( 1,16)/ 2.20098350d-05,-8.88863301d-05/
      data fpppp( 1,17),fpppp( 1,18)/ 1.05888091d-03,-6.94127215d-04/
      data fpppp( 1,19) /            -1.61834238d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/-4.92837054d-04,-4.83909695d-04/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.97578137d-05,-1.20816520d-03/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.95406672d-03,-8.57694626d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 4.40276070d-05,-2.93053646d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 3.58970573d-05, 4.86344338d-05/
      data fpppp( 2,11),fpppp( 2,12)/ 4.46902172d-05, 6.20801104d-06/
      data fpppp( 2,13),fpppp( 2,14)/-1.03741663d-04,-7.53301517d-05/
      data fpppp( 2,15),fpppp( 2,16)/ 1.91695372d-06,-4.91689420d-05/
      data fpppp( 2,17),fpppp( 2,18)/ 4.69915172d-04,-2.07623108d-04/
      data fpppp( 2,19) /            -5.43269275d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-2.22013639d-04, 3.43437863d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-1.17317084d-03, 2.74286671d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-1.31543908d-03,-2.87842656d-04/
      data fpppp( 3, 7),fpppp( 3, 8)/ 5.95715346d-05, 1.36355884d-05/
      data fpppp( 3, 9),fpppp( 3,10)/ 3.51148921d-05, 1.99042992d-05/
      data fpppp( 3,11),fpppp( 3,12)/ 1.78103770d-05, 1.54285956d-05/
      data fpppp( 3,13),fpppp( 3,14)/ 1.17755470d-05,-5.75416080d-05/
      data fpppp( 3,15),fpppp( 3,16)/-3.06229079d-05,-5.53924843d-05/
      data fpppp( 3,17),fpppp( 3,18)/-7.46512407d-04, 1.26892468d-03/
      data fpppp( 3,19) /             2.51105976d-03 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 8.22520663d-04, 7.10715641d-04/
      data fpppp( 4, 3),fpppp( 4, 4)/ 1.51202973d-03,-8.82573073d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-4.80696907d-04, 5.19227365d-05/
      data fpppp( 4, 7),fpppp( 4, 8)/ 3.93083313d-06, 3.15176124d-05/
      data fpppp( 4, 9),fpppp( 4,10)/ 1.36619912d-05, 1.64217013d-05/
      data fpppp( 4,11),fpppp( 4,12)/ 1.12811303d-05, 9.00725288d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 3.03443504d-06, 9.35989764d-06/
      data fpppp( 4,15),fpppp( 4,16)/-1.95983382d-05, 9.91772298d-05/
      data fpppp( 4,17),fpppp( 4,18)/ 2.33176469d-04,-7.19680409d-04/
      data fpppp( 4,19) /            -1.36200697d-03 /
      data fpppp( 5, 1),fpppp( 5, 2)/-1.44546939d-03,-7.86334393d-04/
      data fpppp( 5, 3),fpppp( 5, 4)/ 6.19868765d-05, 2.48156204d-04/
      data fpppp( 5, 5),fpppp( 5, 6)/-6.54163506d-05,-4.42106903d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 3.35273126d-05, 1.57236225d-05/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.27939219d-05, 6.84652069d-06/
      data fpppp( 5,11),fpppp( 5,12)/ 5.53542269d-06, 7.13868457d-06/
      data fpppp( 5,13),fpppp( 5,14)/ 7.12402543d-06,-3.84324519d-07/
      data fpppp( 5,15),fpppp( 5,16)/ 4.11158112d-08,-3.19535137d-05/
      data fpppp( 5,17),fpppp( 5,18)/-1.77703610d-04, 2.16930599d-04/
      data fpppp( 5,19) /             4.57169293d-04 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 3.61261537d-04, 1.57862271d-04/
      data fpppp( 6, 3),fpppp( 6, 4)/-1.81457628d-04,-4.50553131d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/-2.84694115d-05, 2.94999957d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/ 1.66714129d-05, 1.13026547d-05/
      data fpppp( 6, 9),fpppp( 6,10)/ 7.75179551d-06, 2.12756172d-06/
      data fpppp( 6,11),fpppp( 6,12)/ 3.89112164d-06, 3.09009239d-06/
      data fpppp( 6,13),fpppp( 6,14)/ 3.94398636d-06, 4.45922445d-06/
      data fpppp( 6,15),fpppp( 6,16)/-3.79514421d-06,-9.45692244d-06/
      data fpppp( 6,17),fpppp( 6,18)/-7.31482117d-06,-3.94442760d-05/
      data fpppp( 6,19) /            -6.17818577d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 1.31315913d-05,-1.18409853d-05/
      data fpppp( 7, 3),fpppp( 7, 4)/-6.43451346d-05,-2.04887793d-05/
      data fpppp( 7, 5),fpppp( 7, 6)/ 3.02607834d-06, 1.12321865d-05/
      data fpppp( 7, 7),fpppp( 7, 8)/ 9.30015187d-06, 6.85019297d-06/
      data fpppp( 7, 9),fpppp( 7,10)/ 1.05364187d-06, 3.50861484d-06/
      data fpppp( 7,11),fpppp( 7,12)/ 1.14381542d-06, 2.05906477d-06/
      data fpppp( 7,13),fpppp( 7,14)/ 1.95342885d-06, 1.14050897d-06/
      data fpppp( 7,15),fpppp( 7,16)/-1.05586766d-06,-7.58256409d-06/
      data fpppp( 7,17),fpppp( 7,18)/-2.06059066d-05, 1.27030766d-05/
      data fpppp( 7,19) /             3.27742525d-05 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 3.02680824d-05, 1.29161708d-05/
      data fpppp( 8, 3),fpppp( 8, 4)/-1.17188084d-05,-1.08630510d-05/
      data fpppp( 8, 5),fpppp( 8, 6)/-5.15921480d-07, 2.61885675d-06/
      data fpppp( 8, 7),fpppp( 8, 8)/ 2.04017335d-06, 3.72937852d-07/
      data fpppp( 8, 9),fpppp( 8,10)/ 2.00666473d-06,-3.44021858d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 5.47890783d-07, 3.34564584d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 5.38402104d-07, 3.11328505d-07/
      data fpppp( 8,15),fpppp( 8,16)/-8.89777311d-07,-1.40530463d-06/
      data fpppp( 8,17),fpppp( 8,18)/-1.50408472d-06,-9.81173206d-07/
      data fpppp( 8,19) /            -6.35887867d-07 /
      data fpppp( 9, 1),fpppp( 9, 2)/-8.19640611d-07,-1.10585841d-06/
      data fpppp( 9, 3),fpppp( 9, 4)/-2.20486974d-06,-4.59103695d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-2.25606567d-07, 9.49300029d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/ 3.34994745d-07, 8.29352069d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-3.95726571d-07, 3.87079234d-07/
      data fpppp( 9,11),fpppp( 9,12)/ 7.04206639d-08, 1.30273410d-07/
      data fpppp( 9,13),fpppp( 9,14)/ 9.17747816d-09, 3.89215444d-08/
      data fpppp( 9,15),fpppp( 9,16)/ 1.27384025d-07,-9.20190444d-07/
      data fpppp( 9,17),fpppp( 9,18)/-5.25869475d-07, 7.61249056d-07/
      data fpppp( 9,19) /             1.79328257d-06 /
      data fpppp(10, 1),fpppp(10, 2)/ 2.88283771d-06, 1.08714750d-06/
      data fpppp(10, 3),fpppp(10, 4)/-1.50520893d-06,-1.30923358d-06/
      data fpppp(10, 5),fpppp(10, 6)/-1.13758544d-07, 3.19347744d-07/
      data fpppp(10, 7),fpppp(10, 8)/ 2.81055929d-07, 1.30496327d-07/
      data fpppp(10, 9),fpppp(10,10)/ 1.86863471d-07, 3.47747836d-08/
      data fpppp(10,11),fpppp(10,12)/ 7.27251883d-08, 6.68774032d-08/
      data fpppp(10,13),fpppp(10,14)/ 7.40468423d-08, 3.86142542d-08/
      data fpppp(10,15),fpppp(10,16)/-1.27833395d-07,-2.59664114d-07/
      data fpppp(10,17),fpppp(10,18)/-3.27840704d-07, 8.03207886d-08/
      data fpppp(10,19) /             3.52385686d-07 /
      data fpppp(11, 1),fpppp(11, 2)/-5.36779624d-06,-1.82360647d-06/
      data fpppp(11, 3),fpppp(11, 4)/ 3.86053775d-06, 2.05484139d-06/
      data fpppp(11, 5),fpppp(11, 6)/ 1.90247603d-07,-4.68571794d-07/
      data fpppp(11, 7),fpppp(11, 8)/-5.69879606d-07,-7.47343675d-07/
      data fpppp(11, 9),fpppp(11,10)/ 8.80195005d-09,-1.82301624d-07/
      data fpppp(11,11),fpppp(11,12)/-1.38564353d-07,-1.23317436d-07/
      data fpppp(11,13),fpppp(11,14)/-6.35567248d-08,-8.44451778d-08/
      data fpppp(11,15),fpppp(11,16)/-6.56977957d-08, 9.82453081d-07/
      data fpppp(11,17),fpppp(11,18)/ 4.19300750d-07,-1.91028008d-07/
      data fpppp(11,19) /            -7.87877786d-07 /
 
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
      subroutine c2_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 , -5.75790000d-04 /
      data f( 1, 3),f( 1, 4) / -1.30985300d-02 , -1.69531000d-02 /
      data f( 1, 5),f( 1, 6) /  9.90390000d-04 ,  9.02038000d-03 /
      data f( 1, 7),f( 1, 8) /  3.89124000d-03 ,  9.95610000d-04 /
      data f( 1, 9),f( 1,10) /  1.94640000d-04 ,  8.83000000d-06 /
      data f( 1,11),f( 1,12) / -1.09870000d-04 ,  1.64690000d-04 /
      data f( 1,13),f( 1,14) /  1.55339000d-03 ,  3.66488000d-03 /
      data f( 1,15),f( 1,16) /  5.30904000d-03 ,  6.11529000d-03 /
      data f( 1,17),f( 1,18) /  4.39913000d-03 , -2.95721000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 , -1.15502000d-03 /
      data f( 2, 3),f( 2, 4) / -1.02684600d-02 , -1.07027500d-02 /
      data f( 2, 5),f( 2, 6) /  4.67208000d-03 ,  5.89217000d-03 /
      data f( 2, 7),f( 2, 8) /  2.22382000d-03 ,  5.48670000d-04 /
      data f( 2, 9),f( 2,10) /  9.55300000d-05 , -1.74200000d-05 /
      data f( 2,11),f( 2,12) / -5.88500000d-05 ,  2.23300000d-05 /
      data f( 2,13),f( 2,14) /  4.67260000d-04 ,  1.52822000d-03 /
      data f( 2,15),f( 2,16) /  2.97863000d-03 ,  4.39511000d-03 /
      data f( 2,17),f( 2,18) /  4.73483000d-03 ,  1.67178000d-03 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  1.86460000d-04 /
      data f( 3, 3),f( 3, 4) / -2.49639000d-03 ,  4.05650000d-04 /
      data f( 3, 5),f( 3, 6) /  5.05157000d-03 ,  3.43158000d-03 /
      data f( 3, 7),f( 3, 8) /  1.19433000d-03 ,  2.93340000d-04 /
      data f( 3, 9),f( 3,10) /  4.91400000d-05 , -1.77600000d-05 /
      data f( 3,11),f( 3,12) / -3.23500000d-05 , -3.68000000d-06 /
      data f( 3,13),f( 3,14) /  1.33610000d-04 ,  5.09450000d-04 /
      data f( 3,15),f( 3,16) /  1.18075000d-03 ,  1.93170000d-03 /
      data f( 3,17),f( 3,18) /  2.38775000d-03 ,  1.44769000d-03 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  2.25369000d-03 /
      data f( 4, 3),f( 4, 4) /  3.42146000d-03 ,  3.46322000d-03 /
      data f( 4, 5),f( 4, 6) /  3.37714000d-03 ,  1.79856000d-03 /
      data f( 4, 7),f( 4, 8) /  5.94300000d-04 ,  1.48750000d-04 /
      data f( 4, 9),f( 4,10) /  1.73000000d-05 , -1.25700000d-05 /
      data f( 4,11),f( 4,12) / -1.74800000d-05 , -7.46000000d-06 /
      data f( 4,13),f( 4,14) /  3.66100000d-05 ,  1.51070000d-04 /
      data f( 4,15),f( 4,16) /  3.69190000d-04 ,  5.99350000d-04 /
      data f( 4,17),f( 4,18) /  5.23140000d-04 ,  1.04490000d-04 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  3.74350000d-04 /
      data f( 5, 3),f( 5, 4) /  1.70662000d-03 ,  2.24315000d-03 /
      data f( 5, 5),f( 5, 6) /  1.77056000d-03 ,  8.33480000d-04 /
      data f( 5, 7),f( 5, 8) /  2.79220000d-04 ,  7.37400000d-05 /
      data f( 5, 9),f( 5,10) /  7.94000000d-06 , -7.50000000d-06 /
      data f( 5,11),f( 5,12) / -9.18000000d-06 , -5.33000000d-06 /
      data f( 5,13),f( 5,14) /  7.69000000d-06 ,  4.09300000d-05 /
      data f( 5,15),f( 5,16) /  1.07550000d-04 ,  1.64300000d-04 /
      data f( 5,17),f( 5,18) /  1.33930000d-04 ,  2.41900000d-05 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.07920000d-04 /
      data f( 6, 3),f( 6, 4) /  6.96000000d-04 ,  1.03463000d-03 /
      data f( 6, 5),f( 6, 6) /  7.59650000d-04 ,  3.51250000d-04 /
      data f( 6, 7),f( 6, 8) /  1.28910000d-04 ,  3.99300000d-05 /
      data f( 6, 9),f( 6,10) /  5.78000000d-06 , -4.11000000d-06 /
      data f( 6,11),f( 6,12) / -4.76000000d-06 , -4.07000000d-06 /
      data f( 6,13),f( 6,14) /  3.60000000d-07 ,  9.49000000d-06 /
      data f( 6,15),f( 6,16) /  2.54900000d-05 ,  3.96000000d-05 /
      data f( 6,17),f( 6,18) /  3.23300000d-05 ,  8.33000000d-06 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  3.44100000d-05 /
      data f( 7, 3),f( 7, 4) /  2.44880000d-04 ,  3.84830000d-04 /
      data f( 7, 5),f( 7, 6) /  2.85230000d-04 ,  1.37170000d-04 /
      data f( 7, 7),f( 7, 8) /  5.24100000d-05 ,  1.59100000d-05 /
      data f( 7, 9),f( 7,10) / -3.80000000d-06 , -2.20000000d-06 /
      data f( 7,11),f( 7,12) / -2.52000000d-06 , -2.61000000d-06 /
      data f( 7,13),f( 7,14) / -1.72000000d-06 ,  9.10000000d-07 /
      data f( 7,15),f( 7,16) /  5.62000000d-06 ,  9.39000000d-06 /
      data f( 7,17),f( 7,18) /  7.55000000d-06 ,  2.36000000d-06 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  2.87000000d-06 /
      data f( 8, 3),f( 8, 4) /  2.32500000d-05 ,  4.06000000d-05 /
      data f( 8, 5),f( 8, 6) /  3.60000000d-05 ,  2.29100000d-05 /
      data f( 8, 7),f( 8, 8) /  1.06000000d-05 ,  3.74000000d-06 /
      data f( 8, 9),f( 8,10) /  1.85000000d-06 , -8.90000000d-07 /
      data f( 8,11),f( 8,12) / -1.03000000d-06 , -7.50000000d-07 /
      data f( 8,13),f( 8,14) / -7.20000000d-07 , -8.30000000d-07 /
      data f( 8,15),f( 8,16) / -9.40000000d-07 , -1.50000000d-07 /
      data f( 8,17),f( 8,18) / -2.70000000d-07 , -1.00000000d-07 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  2.06000000d-06 /
      data f( 9, 3),f( 9, 4) /  2.01000000d-06 ,  2.62000000d-06 /
      data f( 9, 5),f( 9, 6) /  1.11000000d-06 ,  8.30000000d-07 /
      data f( 9, 7),f( 9, 8) /  7.40000000d-07 ,  2.03000000d-06 /
      data f( 9, 9),f( 9,10) / -1.05000000d-06 , -4.40000000d-07 /
      data f( 9,11),f( 9,12) / -3.80000000d-07 , -5.60000000d-07 /
      data f( 9,13),f( 9,14) / -5.30000000d-07 ,  1.30000000d-06 /
      data f( 9,15),f( 9,16) /  1.58000000d-06 , -8.70000000d-07 /
      data f( 9,17),f( 9,18) /  1.75000000d-06 , -2.70000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  8.10000000d-07 /
      data f(10, 3),f(10, 4) /  2.22000000d-06 , -8.90000000d-07 /
      data f(10, 5),f(10, 6) / -1.19000000d-06 ,  1.17000000d-06 /
      data f(10, 7),f(10, 8) / -9.00000000d-08 , -8.20000000d-07 /
      data f(10, 9),f(10,10) / -6.30000000d-07 ,  2.20000000d-07 /
      data f(10,11),f(10,12) / -7.00000000d-08 , -3.60000000d-07 /
      data f(10,13),f(10,14) / -4.30000000d-07 , -2.80000000d-07 /
      data f(10,15),f(10,16) / -4.60000000d-07 , -6.90000000d-07 /
      data f(10,17),f(10,18) / -5.10000000d-07 , -4.30000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -2.60000000d-07 /
      data f(11, 3),f(11, 4) / -4.90000000d-07 , -2.80000000d-07 /
      data f(11, 5),f(11, 6) / -2.90000000d-07 ,  2.50000000d-07 /
      data f(11, 7),f(11, 8) / -1.60000000d-07 , -5.70000000d-07 /
      data f(11, 9),f(11,10) / -2.70000000d-07 , -5.00000000d-08 /
      data f(11,11),f(11,12) / -9.00000000d-08 , -1.70000000d-07 /
      data f(11,13),f(11,14) / -4.00000000d-08 , -3.00000000d-08 /
      data f(11,15),f(11,16) / -7.00000000d-08 , -8.00000000d-08 /
      data f(11,17),f(11,18) / -6.00000000d-08 , -2.00000000d-08 /
      data f(11,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-3.28164712d-04/
      data fpp( 1, 2,1),fpp( 1, 2,2)/ 1.35806073d-02,-1.14381475d-04/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 4.75114310d-02, 6.88736127d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 6.83668590d-02, 3.58977224d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.84112717d-02,-1.96898909d-04/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.08306686d-03,-1.66191588d-04/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 3.37591637d-03, 7.21174592d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 1.08405767d-03, 1.17323508d-05/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 3.54996985d-04, 6.63273757d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/ 1.82879412d-04,-1.35370109d-06/
      data fpp( 1,11,1),fpp( 1,11,2)/-1.48437279d-04, 2.80866679d-06/
      data fpp( 1,12,1),fpp( 1,12,2)/ 8.29305676d-04, 1.37146339d-05/
      data fpp( 1,13,1),fpp( 1,13,2)/ 5.02076770d-03, 9.18119743d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/ 6.30358064d-03,-7.07202366d-06/
      data fpp( 1,15,1),fpp( 1,15,2)/ 4.62675852d-04,-8.93290278d-06/
      data fpp( 1,16,1),fpp( 1,16,2)/-1.01402451d-02,-7.47096524d-06/
      data fpp( 1,17,1),fpp( 1,17,2)/-2.31399802d-02,-1.12527836d-04/
      data fpp( 1,18,1),fpp( 1,18,2)/-3.42965360d-02, 1.19171510d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 2.54654795d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-2.47309175d-04/
      data fpp( 2, 2,1),fpp( 2, 2,2)/ 5.44698545d-03,-7.68856492d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 1.86508981d-02, 7.73465722d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.48339220d-02, 2.88248360d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-1.27900966d-02,-2.81792814d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 2.56570629d-03,-1.03615058d-05/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.57108727d-03, 2.99324367d-05/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 7.78164657d-04, 1.02237591d-05/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 2.28006029d-04, 2.49312686d-06/
      data fpp( 2,10,1),fpp( 2,10,2)/ 1.08201176d-04, 2.15133441d-07/
      data fpp( 2,11,1),fpp( 2,11,2)/-1.00485442d-04, 9.37539375d-07/
      data fpp( 2,12,1),fpp( 2,12,2)/ 4.90548648d-04, 3.39130906d-06/
      data fpp( 2,13,1),fpp( 2,13,2)/ 3.11486460d-03, 7.32222439d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/ 4.46751871d-03, 4.28159339d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.83468830d-03,-1.08159795d-06/
      data fpp( 2,16,1),fpp( 2,16,2)/-3.63258985d-03,-1.99100158d-06/
      data fpp( 2,17,1),fpp( 2,17,2)/-1.12353996d-02,-5.55599957d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.95156481d-02, 2.00647845d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 5.87770578d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-1.09031797d-04/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 1.07284909d-02,-3.70997062d-05/
      data fpp( 3, 3,1),fpp( 3, 3,2)/-3.50702333d-03, 8.52720216d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-5.11093469d-02, 3.11050198d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-9.68114190d-03,-1.05059301d-04/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 3.67698800d-03, 1.31775840d-05/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 1.65005456d-03, 1.53133649d-05/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 4.01923699d-04, 5.74455621d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.74110282d-06, 1.11581023d-06/
      data fpp( 3,10,1),fpp( 3,10,2)/ 6.15588428d-06, 4.30202882d-07/
      data fpp( 3,11,1),fpp( 3,11,2)/-3.81009543d-05, 3.01978243d-07/
      data fpp( 3,12,1),fpp( 3,12,2)/ 8.99733593d-07, 9.57484147d-07/
      data fpp( 3,13,1),fpp( 3,13,2)/ 5.79293908d-04, 2.38528517d-06/
      data fpp( 3,14,1),fpp( 3,14,2)/ 2.65570451d-03, 3.81437518d-06/
      data fpp( 3,15,1),fpp( 3,15,2)/ 4.97929096d-03, 8.48141090d-08/
      data fpp( 3,16,1),fpp( 3,16,2)/ 6.83308446d-03, 6.25368383d-07/
      data fpp( 3,17,1),fpp( 3,17,2)/ 3.69485852d-03,-2.02802876d-05/
      data fpp( 3,18,1),fpp( 3,18,2)/-4.11479167d-03,-3.27081782d-06/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 2.90575891d-06/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00,-1.12572260d-05/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-3.09429492d-02,-9.26134797d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-4.91240848d-02,-1.68525821d-05/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-1.36164542d-02, 9.11107633d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 2.22058419d-03,-2.72621232d-05/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 2.58802173d-03, 1.03874166d-05/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 1.13573451d-03, 8.17165666d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 2.71900546d-04, 2.44855670d-06/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 1.28158382d-04, 8.80116522d-07/
      data fpp( 4,10,1),fpp( 4,10,2)/-1.04713031d-07, 1.25777208d-07/
      data fpp( 4,11,1),fpp( 4,11,2)/-2.62307413d-05, 1.14374646d-07/
      data fpp( 4,12,1),fpp( 4,12,2)/ 3.93724181d-05, 3.12524209d-07/
      data fpp( 4,13,1),fpp( 4,13,2)/ 2.47559771d-04, 6.78528519d-07/
      data fpp( 4,14,1),fpp( 4,14,2)/ 7.59023263d-04, 1.19676171d-06/
      data fpp( 4,15,1),fpp( 4,15,2)/ 1.91982786d-03, 7.54024624d-07/
      data fpp( 4,16,1),fpp( 4,16,2)/ 3.44569201d-03,-3.49046021d-06/
      data fpp( 4,17,1),fpp( 4,17,2)/ 8.03524550d-03,-5.17438379d-06/
      data fpp( 4,18,1),fpp( 4,18,2)/ 9.11617476d-03, 3.64159537d-06/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 9.45760232d-06/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 2.65459163d-05/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.83256258d-02, 1.07189674d-05/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 1.68188024d-02,-1.19465860d-05/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 2.91180373d-03,-1.06770234d-05/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 2.42720513d-03,-5.89252027d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 2.00148510d-03, 6.37770453d-06/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 6.45807410d-04, 3.35090217d-06/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 1.80394116d-04, 1.14548679d-06/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 2.86275752d-05, 4.47950682d-07/
      data fpp( 5,10,1),fpp( 5,10,2)/-8.61703215d-06, 8.43104842d-08/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.46560806d-05, 4.04073807d-08/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.65494059d-05, 8.58599929d-08/
      data fpp( 5,13,1),fpp( 5,13,2)/ 6.43870066d-05, 1.66352648d-07/
      data fpp( 5,14,1),fpp( 5,14,2)/ 2.65962442d-04, 4.61929416d-07/
      data fpp( 5,15,1),fpp( 5,15,2)/ 5.39477610d-04,-1.12703132d-08/
      data fpp( 5,16,1),fpp( 5,16,2)/ 9.19347503d-04,-1.00904816d-06/
      data fpp( 5,17,1),fpp( 5,17,2)/-4.26240534d-04,-1.17973703d-06/
      data fpp( 5,18,1),fpp( 5,18,2)/-2.04030736d-03, 9.65796295d-07/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 2.44955185d-06/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 1.20418914d-05/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-3.64971390d-03, 4.91321711d-06/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.24984500d-03,-2.88515989d-06/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 2.24643926d-03,-8.33957755d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 2.36667528d-03,-5.73129906d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 9.94437888d-04, 2.62689717d-06/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 2.35515851d-04, 1.22914121d-06/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-4.67701013d-06, 4.58137986d-07/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-6.98686828d-05, 2.28106846d-07/
      data fpp( 6,10,1),fpp( 6,10,2)/-5.74715836d-06, 8.50346311d-08/
      data fpp( 6,11,1),fpp( 6,11,2)/-8.26493635d-06,-1.38453700d-08/
      data fpp( 6,12,1),fpp( 6,12,2)/ 5.94520560d-06, 5.07468489d-08/
      data fpp( 6,13,1),fpp( 6,13,2)/ 1.30522022d-05, 3.52579744d-08/
      data fpp( 6,14,1),fpp( 6,14,2)/ 6.59269700d-05, 9.02212533d-08/
      data fpp( 6,15,1),fpp( 6,15,2)/ 2.32181703d-04, 1.60570122d-08/
      data fpp( 6,16,1),fpp( 6,16,2)/ 3.25317977d-04,-2.67849302d-07/
      data fpp( 6,17,1),fpp( 6,17,2)/ 5.72356632d-04,-2.27459803d-07/
      data fpp( 6,18,1),fpp( 6,18,2)/ 5.91614698d-04, 1.73888515d-07/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00, 4.72105742d-07/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 4.23566340d-06/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 9.03309824d-04, 1.74207319d-06/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 1.60857756d-03,-6.40356176d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.51171923d-03,-3.41184849d-06/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 9.81853737d-04,-8.52498727d-08/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 4.56363350d-04, 8.45247979d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 1.83569187d-04, 5.02257958d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 7.32739245d-05, 4.13201882d-08/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 7.27671559d-05, 3.39861289d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/-3.91433442d-06,-1.22165344d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/-4.60417402d-06, 3.36000880d-08/
      data fpp( 7,12,1),fpp( 7,12,2)/-2.43141649d-06, 1.56499209d-09/
      data fpp( 7,13,1),fpp( 7,13,2)/ 9.40418465d-06, 1.89399436d-08/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.89696782d-05, 2.70752335d-08/
      data fpp( 7,15,1),fpp( 7,15,2)/ 2.43555775d-05,-2.44087748d-09/
      data fpp( 7,16,1),fpp( 7,16,2)/ 4.71405893d-05,-7.37117236d-08/
      data fpp( 7,17,1),fpp( 7,17,2)/-1.95059945d-05,-3.93122283d-08/
      data fpp( 7,18,1),fpp( 7,18,2)/-8.87914266d-05, 2.99606366d-08/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00, 8.92696817d-08/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 3.83519222d-07/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-1.92192523d-04, 1.69061557d-07/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-1.17150193d-04,-9.16544874d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 7.38426698d-05,-3.14199762d-07/
      data fpp( 8, 5,1),fpp( 8, 5,2)/ 6.87611479d-05,-5.10355040d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/ 1.70910048d-05, 8.94177790d-09/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.32548669d-06, 6.20683924d-08/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-2.26326854d-06, 6.97846524d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-3.45071262d-05,-4.30070019d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/-4.43417548d-07, 5.12433551d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/ 4.99023042d-09,-5.96641861d-09/
      data fpp( 8,12,1),fpp( 8,12,2)/-2.03835333d-06,-2.17768068d-09/
      data fpp( 8,13,1),fpp( 8,13,2)/-3.77865505d-06,-3.22858667d-10/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.64748049d-06,-4.93088465d-09/
      data fpp( 8,15,1),fpp( 8,15,2)/ 9.92241588d-06, 2.00463973d-08/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.19924371d-06,-2.12547044d-08/
      data fpp( 8,17,1),fpp( 8,17,2)/ 2.27796676d-05, 1.03724204d-08/
      data fpp( 8,18,1),fpp( 8,18,2)/ 2.74469311d-05,-2.83497726d-09/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-3.23251137d-09/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-4.56453118d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 4.98402693d-05,-2.74093763d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 6.23632088d-05, 2.86828171d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 3.04100869d-05,-4.77218921d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 2.91416716d-05, 3.50047514d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 2.83526303d-05,-1.84971134d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.34327597d-05, 5.03837021d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-1.46085038d-06,-1.00237695d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.39613490d-05, 8.83670773d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/ 5.28004618d-07,-3.18306145d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/-4.55786903d-07, 5.95538083d-09/
      data fpp( 9,12,1),fpp( 9,12,2)/ 5.64829825d-07,-6.39090876d-09/
      data fpp( 9,13,1),fpp( 9,13,2)/ 8.50435539d-07, 3.22082542d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-6.33960013d-06,-1.44421081d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/-9.56524102d-06,-6.74398219d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 9.82435866d-07, 1.20401396d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.25726757d-05,-1.09965760d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/-7.25629765d-06, 4.10616458d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 8.31191771d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 5.14820754d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-9.80855386d-06, 1.74358492d-08/
      data fpp(10, 3,1),fpp(10, 3,2)/-3.60264176d-06,-8.52254722d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.13369826d-05, 5.22660395d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/ 1.02121657d-05, 4.47613142d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/ 4.01847394d-06,-7.17112961d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/ 1.77444807d-06, 2.48838703d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 1.26667008d-06, 3.97581485d-09/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.41826980d-06, 1.44128703d-08/
      data fpp(10,10,1),fpp(10,10,2)/-4.08600924d-07,-2.20272960d-08/
      data fpp(10,11,1),fpp(10,11,2)/-2.21842619d-07, 5.29631374d-09/
      data fpp(10,12,1),fpp(10,12,2)/-1.60965965d-07, 8.42041032d-10/
      data fpp(10,13,1),fpp(10,13,2)/-1.63087108d-07, 4.53552213d-09/
      data fpp(10,14,1),fpp(10,14,2)/ 4.50920026d-07,-5.78412954d-09/
      data fpp(10,15,1),fpp(10,15,2)/ 9.78548204d-07,-1.19900396d-09/
      data fpp(10,16,1),fpp(10,16,2)/ 2.71012827d-07, 7.58014539d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 1.83103513d-06,-4.52157761d-09/
      data fpp(10,18,1),fpp(10,18,2)/ 1.63825953d-06, 4.50616503d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 7.49691748d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00,-3.11096506d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 6.65052693d-06,-1.07806987d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-2.50686791d-05, 9.22324456d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/-3.77709913d-05,-9.41490836d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-3.69573328d-05, 1.52363889d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-2.86317370d-05,-1.85306472d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-9.65472403d-06, 1.88619986d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 5.85541496d-06, 1.09858478d-08/
      data fpp(11, 9,1),fpp(11, 9,2)/-3.44586510d-06,-3.22959087d-09/
      data fpp(11,10,1),fpp(11,10,2)/-1.42319954d-06,-2.86748428d-09/
      data fpp(11,11,1),fpp(11,11,2)/-6.65786903d-08,-9.00472024d-10/
      data fpp(11,12,1),fpp(11,12,2)/-1.14517018d-07, 4.06937237d-09/
      data fpp(11,13,1),fpp(11,13,2)/ 3.49043554d-07,-2.77701747d-09/
      data fpp(11,14,1),fpp(11,14,2)/ 6.93203999d-06,-1.61302499d-10/
      data fpp(11,15,1),fpp(11,15,2)/ 8.55197590d-06, 4.22227465d-10/
      data fpp(11,16,1),fpp(11,16,2)/-9.29256413d-07, 2.72392638d-10/
      data fpp(11,17,1),fpp(11,17,2)/ 8.24823243d-06, 2.88201982d-10/
      data fpp(11,18,1),fpp(11,18,2)/-1.91629765d-07,-2.25200566d-10/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00,-5.87399717d-10/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 6.08782873d-04, 6.14529859d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/ 3.66418167d-04,-2.31164939d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/ 2.42216588d-03,-9.40665957d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.88408609d-04,-2.80509728d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.75631621d-05,-8.78508889d-06/
      data fpppp( 1,11),fpppp( 1,12)/ 8.02524639d-06, 5.52276821d-05/
      data fpppp( 1,13),fpppp( 1,14)/-3.61128307d-05,-8.52953043d-05/
      data fpppp( 1,15),fpppp( 1,16)/-5.01290163d-05, 9.04012818d-08/
      data fpppp( 1,17),fpppp( 1,18)/-9.40414412d-05, 4.86666127d-04/
      data fpppp( 1,19) /             8.74562435d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 2.54893542d-04, 1.84770505d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 1.36613888d-04,-9.86185927d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 1.17970727d-03,-5.53853869d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.14682890d-04,-1.27759062d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 1.09865739d-05,-5.34916307d-06/
      data fpppp( 2,11),fpppp( 2,12)/ 5.07717250d-06, 3.30237155d-05/
      data fpppp( 2,13),fpppp( 2,14)/-1.51751227d-05,-4.86229347d-05/
      data fpppp( 2,15),fpppp( 2,16)/-2.94622102d-05,-3.59508799d-06/
      data fpppp( 2,17),fpppp( 2,18)/-8.42893332d-05, 3.00106093d-04/
      data fpppp( 2,19) /             5.51618757d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-2.50232848d-04,-8.03983756d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-9.26013961d-04, 1.78244566d-03/
      data fpppp( 3, 5),fpppp( 3, 6)/-8.61936948d-04,-1.89023758d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/ 1.44426507d-05, 7.85992813d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 4.78560004d-06,-2.30862097d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 1.31965429d-06, 2.02545539d-06/
      data fpppp( 3,13),fpppp( 3,14)/ 2.29421333d-05,-3.91300320d-06/
      data fpppp( 3,15),fpppp( 3,16)/ 7.54043092d-06,-5.44362979d-05/
      data fpppp( 3,17),fpppp( 3,18)/-8.93164054d-05, 1.31416465d-04/
      data fpppp( 3,19) /             2.79117056d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.26029530d-04, 1.63744156d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/ 9.26240682d-04,-5.00011174d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/-1.06431520d-04,-2.43879878d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/ 7.00322965d-06, 9.73307553d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-2.73002388d-06, 2.11576417d-06/
      data fpppp( 4,11),fpppp( 4,12)/ 3.95191189d-07, 1.80722233d-06/
      data fpppp( 4,13),fpppp( 4,14)/ 9.30971142d-07, 1.26654614d-05/
      data fpppp( 4,15),fpppp( 4,16)/-1.26323506d-05, 5.97675145d-05/
      data fpppp( 4,17),fpppp( 4,18)/-4.26163468d-05,-9.98195817d-05/
      data fpppp( 4,19) /            -1.69931567d-04 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.84963317d-04,-1.73692314d-04/
      data fpppp( 5, 3),fpppp( 5, 4)/-2.10214374d-04, 2.70539287d-04/
      data fpppp( 5, 5),fpppp( 5, 6)/-6.65987674d-05,-6.11503360d-07/
      data fpppp( 5, 7),fpppp( 5, 8)/ 1.32473218d-05, 1.03807979d-06/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.41916428d-06, 1.56579085d-07/
      data fpppp( 5,11),fpppp( 5,12)/-1.73147086d-07, 7.84752648d-07/
      data fpppp( 5,13),fpppp( 5,14)/ 2.00392077d-06,-1.56209436d-06/
      data fpppp( 5,15),fpppp( 5,16)/ 8.56084065d-06,-2.62999847d-05/
      data fpppp( 5,17),fpppp( 5,18)/-6.88837759d-06, 3.77447676d-05/
      data fpppp( 5,19) /             7.51717590d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.10510785d-04, 5.95292619d-05/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.43471347d-05,-5.11328789d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/-1.23785132d-05, 1.10985267d-05/
      data fpppp( 6, 7),fpppp( 6, 8)/ 4.78332795d-06, 8.91912127d-07/
      data fpppp( 6, 9),fpppp( 6,10)/ 2.14909483d-06,-1.72949962d-06/
      data fpppp( 6,11),fpppp( 6,12)/ 7.70545502d-07,-3.49007193d-07/
      data fpppp( 6,13),fpppp( 6,14)/ 1.99294549d-07, 2.29789527d-06/
      data fpppp( 6,15),fpppp( 6,16)/-2.58807772d-06, 3.66730803d-06/
      data fpppp( 6,17),fpppp( 6,18)/-2.84701149d-06,-5.94609746d-06/
      data fpppp( 6,19) /            -1.00209644d-05 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 3.66729884d-06,-1.19418040d-06/
      data fpppp( 7, 3),fpppp( 7, 4)/-1.07731022d-05,-3.84097500d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/ 1.56572304d-07, 3.47719241d-06/
      data fpppp( 7, 7),fpppp( 7, 8)/ 1.09643147d-06, 1.88701576d-06/
      data fpppp( 7, 9),fpppp( 7,10)/-2.05718488d-06, 1.77124046d-06/
      data fpppp( 7,11),fpppp( 7,12)/-4.68277913d-07, 2.73627019d-07/
      data fpppp( 7,13),fpppp( 7,14)/-4.64595474d-08,-2.23995287d-07/
      data fpppp( 7,15),fpppp( 7,16)/ 6.91665047d-07,-1.49871815d-06/
      data fpppp( 7,17),fpppp( 7,18)/-6.26881664d-08, 1.59113993d-06/
      data fpppp( 7,19) /             3.18273999d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/ 4.30458227d-06, 2.43356745d-06/
      data fpppp( 8, 3),fpppp( 8, 4)/ 1.99523912d-06,-3.45749194d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 7.02655242d-08, 3.81112579d-07/
      data fpppp( 8, 7),fpppp( 8, 8)/ 4.00503248d-07,-9.34402989d-07/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.45874416d-06,-9.22119666d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 2.12816450d-07,-7.86512157d-08/
      data fpppp( 8,13),fpppp( 8,14)/ 1.19970924d-07, 8.87537561d-08/
      data fpppp( 8,15),fpppp( 8,16)/-4.24057957d-07, 6.47591620d-07/
      data fpppp( 8,17),fpppp( 8,18)/-3.48092764d-07,-2.70010186d-07/
      data fpppp( 8,19) /            -4.98718168d-07 /
      data fpppp( 9, 1),fpppp( 9, 2)/-3.47750042d-07,-2.80845178d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-7.67909034d-07, 6.83917628d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.26679080d-07,-1.48438869d-07/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.27415205d-07, 6.59675327d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-6.92337536d-07, 3.78342190d-07/
      data fpppp( 9,11),fpppp( 9,12)/-7.40580517d-08, 3.81545118d-08/
      data fpppp( 9,13),fpppp( 9,14)/-1.22660656d-07, 3.94963089d-09/
      data fpppp( 9,15),fpppp( 9,16)/ 3.44725819d-07,-5.56453842d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 4.34922244d-07,-5.09457604d-08/
      data fpppp( 9,19) /            -1.14744025d-07 /
      data fpppp(10, 1),fpppp(10, 2)/ 2.44169419d-07, 1.37710215d-07/
      data fpppp(10, 3),fpppp(10, 4)/ 1.65857681d-07,-2.77118201d-07/
      data fpppp(10, 5),fpppp(10, 6)/-2.12513583d-08, 5.79911462d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 2.62667257d-08,-5.88831764d-08/
      data fpppp(10, 9),fpppp(10,10)/ 7.86362666d-08,-3.39853646d-08/
      data fpppp(10,11),fpppp(10,12)/ 7.93055751d-09,-5.28976438d-09/
      data fpppp(10,13),fpppp(10,14)/ 9.44863219d-09, 4.46293220d-09/
      data fpppp(10,15),fpppp(10,16)/-3.24830983d-08, 5.13596478d-08/
      data fpppp(10,17),fpppp(10,18)/-3.69020317d-08,-8.91939566d-09/
      data fpppp(10,19) /            -1.41494212d-08 /
      data fpppp(11, 1),fpppp(11, 2)/-9.31904958d-07,-4.35014610d-07/
      data fpppp(11, 3),fpppp(11, 4)/ 3.69779418d-07, 9.69105709d-08/
      data fpppp(11, 5),fpppp(11, 6)/ 5.35365389d-08, 1.39659518d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 2.69104150d-08,-4.55313614d-07/
      data fpppp(11, 9),fpppp(11,10)/ 3.05658898d-07,-8.78852391d-08/
      data fpppp(11,11),fpppp(11,12)/ 5.91937592d-09,-2.00658151d-08/
      data fpppp(11,13),fpppp(11,14)/ 1.05033818d-07,-3.29033071d-08/
      data fpppp(11,15),fpppp(11,16)/-2.71204221d-07, 4.51650099d-07/
      data fpppp(11,17),fpppp(11,18)/-4.15872906d-07, 1.54800461d-07/
      data fpppp(11,19) /             3.14560579d-07 /

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
      subroutine c3_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  1.80100000d-04 /
      data f( 1, 3),f( 1, 4) / -3.36133000d-03 , -7.29093000d-03 /
      data f( 1, 5),f( 1, 6) / -1.22053000d-03 ,  2.17421000d-03 /
      data f( 1, 7),f( 1, 8) /  7.14860000d-04 ,  1.47600000d-04 /
      data f( 1, 9),f( 1,10) /  1.78600000d-05 ,  4.10000000d-07 /
      data f( 1,11),f( 1,12) / -2.03000000d-06 ,  3.47000000d-06 /
      data f( 1,13),f( 1,14) /  2.58470000d-04 ,  8.78140000d-04 /
      data f( 1,15),f( 1,16) /  1.46846000d-03 ,  1.86258000d-03 /
      data f( 1,17),f( 1,18) /  1.25701000d-03 , -3.80520000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  4.85430000d-04 /
      data f( 2, 3),f( 2, 4) / -2.45678000d-03 , -4.67064000d-03 /
      data f( 2, 5),f( 2, 6) /  5.48380000d-04 ,  1.18202000d-03 /
      data f( 2, 7),f( 2, 8) /  3.93200000d-04 ,  8.03000000d-05 /
      data f( 2, 9),f( 2,10) /  8.48000000d-06 , -6.10000000d-07 /
      data f( 2,11),f( 2,12) / -1.96000000d-06 , -2.62000000d-06 /
      data f( 2,13),f( 2,14) /  3.10000000d-05 ,  2.11560000d-04 /
      data f( 2,15),f( 2,16) /  5.78200000d-04 ,  1.08659000d-03 /
      data f( 2,17),f( 2,18) /  1.45331000d-03 , -4.15390000d-04 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 ,  8.75300000d-05 /
      data f( 3, 3),f( 3, 4) / -6.36200000d-04 , -3.28520000d-04 /
      data f( 3, 5),f( 3, 6) /  7.86740000d-04 ,  6.39290000d-04 /
      data f( 3, 7),f( 3, 8) /  2.13430000d-04 ,  4.02600000d-05 /
      data f( 3, 9),f( 3,10) /  5.76000000d-06 , -8.60000000d-07 /
      data f( 3,11),f( 3,12) / -1.38000000d-06 , -7.40000000d-07 /
      data f( 3,13),f( 3,14) /  2.42000000d-06 ,  3.25200000d-05 /
      data f( 3,15),f( 3,16) /  1.30490000d-04 ,  2.94310000d-04 /
      data f( 3,17),f( 3,18) /  4.92460000d-04 ,  3.15980000d-04 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  5.81570000d-04 /
      data f( 4, 3),f( 4, 4) /  9.70360000d-04 ,  5.00330000d-04 /
      data f( 4, 5),f( 4, 6) /  5.40550000d-04 ,  3.41790000d-04 /
      data f( 4, 7),f( 4, 8) /  1.02270000d-04 ,  1.52700000d-05 /
      data f( 4, 9),f( 4,10) / -2.17000000d-06 , -7.40000000d-07 /
      data f( 4,11),f( 4,12) / -9.80000000d-07 , -1.29000000d-06 /
      data f( 4,13),f( 4,14) /  1.32000000d-06 ,  4.74000000d-06 /
      data f( 4,15),f( 4,16) /  1.91800000d-05 ,  4.33200000d-05 /
      data f( 4,17),f( 4,18) /  3.35300000d-05 ,  2.35000000d-06 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  2.61000000d-05 /
      data f( 5, 3),f( 5, 4) /  2.05500000d-04 ,  2.91470000d-04 /
      data f( 5, 5),f( 5, 6) /  3.05150000d-04 ,  1.56380000d-04 /
      data f( 5, 7),f( 5, 8) /  4.17500000d-05 ,  4.88000000d-06 /
      data f( 5, 9),f( 5,10) / -1.53000000d-06 , -6.50000000d-07 /
      data f( 5,11),f( 5,12) / -5.70000000d-07 , -9.00000000d-07 /
      data f( 5,13),f( 5,14) / -2.90000000d-07 ,  1.18000000d-06 /
      data f( 5,15),f( 5,16) /  7.70000000d-06 ,  6.52000000d-06 /
      data f( 5,17),f( 5,18) /  4.47000000d-06 ,  3.30000000d-07 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  3.96000000d-06 /
      data f( 6, 3),f( 6, 4) /  5.91400000d-05 ,  1.39380000d-04 /
      data f( 6, 5),f( 6, 6) /  1.32110000d-04 ,  5.89300000d-05 /
      data f( 6, 7),f( 6, 8) /  1.62000000d-05 ,  4.39000000d-06 /
      data f( 6, 9),f( 6,10) /  4.20000000d-07 , -4.10000000d-07 /
      data f( 6,11),f( 6,12) / -3.30000000d-07 , -1.48000000d-06 /
      data f( 6,13),f( 6,14) / -6.80000000d-07 ,  5.00000000d-08 /
      data f( 6,15),f( 6,16) /  5.40000000d-07 ,  1.20000000d-07 /
      data f( 6,17),f( 6,18) /  7.90000000d-07 ,  2.51000000d-06 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.40000000d-06 /
      data f( 7, 3),f( 7, 4) /  1.87400000d-05 ,  5.11000000d-05 /
      data f( 7, 5),f( 7, 6) /  4.44300000d-05 ,  1.52900000d-05 /
      data f( 7, 7),f( 7, 8) /  6.50000000d-07 , -1.44000000d-06 /
      data f( 7, 9),f( 7,10) / -5.52000000d-06 , -1.70000000d-07 /
      data f( 7,11),f( 7,12) / -2.20000000d-07 , -1.19000000d-06 /
      data f( 7,13),f( 7,14) / -1.40000000d-06 , -8.00000000d-07 /
      data f( 7,15),f( 7,16) /  1.70000000d-07 ,  5.40000000d-07 /
      data f( 7,17),f( 7,18) /  7.30000000d-07 ,  1.18000000d-06 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  1.10000000d-07 /
      data f( 8, 3),f( 8, 4) /  1.44000000d-06 ,  3.13000000d-06 /
      data f( 8, 5),f( 8, 6) /  3.54000000d-06 ,  1.98000000d-06 /
      data f( 8, 7),f( 8, 8) / -3.10000000d-07 , -2.40000000d-07 /
      data f( 8, 9),f( 8,10) /  1.66000000d-06 ,  1.00000000d-08 /
      data f( 8,11),f( 8,12) / -7.00000000d-08 , -4.80000000d-07 /
      data f( 8,13),f( 8,14) / -5.20000000d-07 , -6.40000000d-07 /
      data f( 8,15),f( 8,16) / -8.60000000d-07 , -3.20000000d-07 /
      data f( 8,17),f( 8,18) / -4.40000000d-07 , -1.20000000d-07 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  1.72000000d-06 /
      data f( 9, 3),f( 9, 4) / -8.00000000d-08 , -1.13000000d-06 /
      data f( 9, 5),f( 9, 6) / -2.45000000d-06 , -1.93000000d-06 /
      data f( 9, 7),f( 9, 8) / -6.70000000d-07 ,  2.04000000d-06 /
      data f( 9, 9),f( 9,10) / -6.30000000d-07 ,  6.10000000d-07 /
      data f( 9,11),f( 9,12) / -3.90000000d-07 , -5.20000000d-07 /
      data f( 9,13),f( 9,14) / -4.90000000d-07 ,  1.32000000d-06 /
      data f( 9,15),f( 9,16) /  1.61000000d-06 , -7.80000000d-07 /
      data f( 9,17),f( 9,18) /  1.70000000d-06 , -2.60000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  7.50000000d-07 /
      data f(10, 3),f(10, 4) /  2.10000000d-06 , -1.12000000d-06 /
      data f(10, 5),f(10, 6) / -1.17000000d-06 ,  1.53000000d-06 /
      data f(10, 7),f(10, 8) /  5.30000000d-07 , -3.20000000d-07 /
      data f(10, 9),f(10,10) / -4.80000000d-07 ,  7.60000000d-07 /
      data f(10,11),f(10,12) / -1.00000000d-07 , -3.80000000d-07 /
      data f(10,13),f(10,14) / -4.20000000d-07 , -2.50000000d-07 /
      data f(10,15),f(10,16) / -4.30000000d-07 , -6.50000000d-07 /
      data f(10,17),f(10,18) / -4.90000000d-07 , -4.20000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -7.60000000d-07 /
      data f(11, 3),f(11, 4) / -4.70000000d-07 , -9.00000000d-08 /
      data f(11, 5),f(11, 6) / -1.20000000d-07 , -5.60000000d-07 /
      data f(11, 7),f(11, 8) / -2.50000000d-07 , -1.40000000d-07 /
      data f(11, 9),f(11,10) /  3.00000000d-08 , -1.50000000d-07 /
      data f(11,11),f(11,12) / -1.00000000d-08 , -7.00000000d-08 /
      data f(11,13),f(11,14) / -3.80000000d-07 , -1.30000000d-07 /
      data f(11,15),f(11,16) / -2.50000000d-07 , -2.50000000d-07 /
      data f(11,17),f(11,18) / -2.10000000d-07 , -2.40000000d-07 /
      data f(11,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-7.62708762d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-8.42420098d-03,-2.57713475d-05/
      data fpp( 1, 3,1),fpp( 1, 3,2)/ 8.62624817d-03,-4.39355336d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 2.65007612d-02, 1.78223282d-04/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.02344424d-02,-6.89575946d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.60796555d-03,-6.29325037d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 8.51814300d-04, 2.94422093d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 1.55954644d-04,-1.31093362d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 6.90503624d-05, 2.05272517d-06/
      data fpp( 1,10,1),fpp( 1,10,2)/ 4.70940490d-06,-1.62567051d-07/
      data fpp( 1,11,1),fpp( 1,11,2)/ 4.63555984d-06,-5.01856963d-07/
      data fpp( 1,12,1),fpp( 1,12,2)/ 7.09093868d-05, 2.64639490d-06/
      data fpp( 1,13,1),fpp( 1,13,2)/ 1.45878971d-03, 4.88627735d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/ 3.26503642d-03,-3.11304293d-07/
      data fpp( 1,15,1),fpp( 1,15,2)/ 2.22149095d-03,-5.40206018d-06/
      data fpp( 1,16,1),fpp( 1,16,2)/-2.14428717d-03, 1.01475450d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/-1.09971154d-02,-9.51695198d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.71223789d-02, 1.03132134d-04/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 2.14685583d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-7.94295775d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-4.35171805d-03,-2.70898451d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/ 2.78026366d-03,-1.78694422d-05/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 9.54123752d-03, 1.42268614d-04/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-6.26571519d-03,-1.05232214d-04/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.81142890d-03, 3.53744119d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 5.85291400d-04, 5.73484904d-06/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 1.12890712d-04, 2.07836266d-06/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 3.67792753d-05, 4.16500329d-07/
      data fpp( 2,10,1),fpp( 2,10,2)/ 3.02119020d-06, 1.94360246d-08/
      data fpp( 2,11,1),fpp( 2,11,2)/ 2.36888033d-06,-2.98444276d-08/
      data fpp( 2,12,1),fpp( 2,12,2)/ 3.70212264d-05, 1.41341686d-07/
      data fpp( 2,13,1),fpp( 2,13,2)/ 8.40380577d-04, 1.52127768d-06/
      data fpp( 2,14,1),fpp( 2,14,2)/ 2.01064717d-03, 2.58994758d-06/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.71681809d-03,-7.16267998d-07/
      data fpp( 2,16,1),fpp( 2,16,2)/-3.67545669d-04, 8.78012441d-06/
      data fpp( 2,17,1),fpp( 2,17,2)/-5.16412922d-03,-4.29044296d-05/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.05640422d-02, 2.87123942d-05/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 6.51002529d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-2.59809378d-05/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 8.95355317d-03,-9.22932437d-06/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 2.23741719d-03, 1.42226353d-05/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-2.33417913d-02, 1.42233832d-05/
      data fpp( 3, 5,1),fpp( 3, 5,2)/-1.43589682d-03,-2.26613680d-05/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 9.33358848d-04, 6.59488862d-07/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 2.12380101d-04, 3.31881256d-06/
      data fpp( 3, 8,1),fpp( 3, 8,2)/ 4.67225063d-05, 1.22666089d-06/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-5.63274634d-05, 9.47438907d-08/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.68583429d-06, 6.71635508d-08/
      data fpp( 3,11,1),fpp( 3,11,2)/-1.87108115d-06, 2.60190620d-09/
      data fpp( 3,12,1),fpp( 3,12,2)/-2.77142923d-05,-7.97117559d-09/
      data fpp( 3,13,1),fpp( 3,13,2)/-4.69520192d-05, 1.80482796d-07/
      data fpp( 3,14,1),fpp( 3,14,2)/ 3.93334917d-04, 9.02439991d-07/
      data fpp( 3,15,1),fpp( 3,15,2)/ 1.53243668d-03, 2.81957240d-07/
      data fpp( 3,16,1),fpp( 3,16,2)/ 3.22350984d-03, 1.92073105d-06/
      data fpp( 3,17,1),fpp( 3,17,2)/ 3.88203228d-03,-5.90508144d-06/
      data fpp( 3,18,1),fpp( 3,18,2)/-4.42401220d-03,-7.78205304d-07/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 6.47902652d-07/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 3.62609741d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.00559346d-02, 2.85205170d-07/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.68664124d-02,-1.63337181d-05/
      data fpp( 4, 4,1),fpp( 4, 4,2)/-4.92552197d-04, 1.35204672d-05/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 3.80102470d-04,-7.13315077d-06/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 3.40655706d-04, 6.73335852d-07/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 2.11828197d-04, 1.99420736d-06/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 6.14192622d-05, 5.01034713d-07/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 6.34905785d-05, 1.75253789d-07/
      data fpp( 4,10,1),fpp( 4,10,2)/-8.84527349d-07,-6.98498701d-08/
      data fpp( 4,11,1),fpp( 4,11,2)/ 7.95444285d-07, 3.94569114d-09/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.55159429d-05, 4.98671056d-08/
      data fpp( 4,13,1),fpp( 4,13,2)/ 6.94749979d-06,-2.82141134d-08/
      data fpp( 4,14,1),fpp( 4,14,2)/ 4.62531658d-05, 1.11589348d-07/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.27035188d-04, 2.43056722d-07/
      data fpp( 4,16,1),fpp( 4,16,2)/ 4.64466298d-04,-5.01816235d-07/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.68208009d-03,-2.71591783d-07/
      data fpp( 4,18,1),fpp( 4,18,2)/ 3.18009102d-03, 3.04783367d-07/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00, 7.82258317d-07/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 3.86735300d-06/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 6.08194535d-03, 1.79889401d-06/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 8.31415248d-03,-1.86492903d-06/
      data fpp( 5, 4,1),fpp( 5, 4,2)/ 4.06960120d-04, 5.50221277d-08/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 1.74446940d-04,-2.69255948d-06/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 3.94178327d-04, 9.68215778d-07/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 1.55667112d-04, 8.68096364d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 5.80004448d-05, 2.24998764d-07/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 8.04514962d-06, 5.95085784d-08/
      data fpp( 5,10,1),fpp( 5,10,2)/ 1.13227511d-06,-2.56330781d-08/
      data fpp( 5,11,1),fpp( 5,11,2)/-1.07069599d-06,-4.97626620d-09/
      data fpp( 5,12,1),fpp( 5,12,2)/-1.17894791d-05, 2.09381429d-08/
      data fpp( 5,13,1),fpp( 5,13,2)/ 6.92202002d-06,-2.23763052d-08/
      data fpp( 5,14,1),fpp( 5,14,2)/ 2.93241988d-06, 1.20167078d-07/
      data fpp( 5,15,1),fpp( 5,15,2)/-4.46574329d-05,-1.55292007d-07/
      data fpp( 5,16,1),fpp( 5,16,2)/ 5.91849653d-05, 3.90009504d-08/
      data fpp( 5,17,1),fpp( 5,17,2)/-2.93472648d-04,-5.29117943d-08/
      data fpp( 5,18,1),fpp( 5,18,2)/-8.17711886d-04, 4.72462270d-08/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00, 9.25268865d-08/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 8.21591614d-07/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.47192676d-03, 4.16616772d-07/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.54619749d-03, 5.85141299d-07/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 2.27191717d-04,-1.25358197d-06/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 4.18749768d-04,-8.21413427d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.93670985d-04, 5.84635677d-07/
      data fpp( 6, 7,1),fpp( 6, 7,2)/ 4.78335628d-06, 3.09870721d-07/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-5.58210416d-05, 3.10814398d-08/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-6.42311769d-05, 3.62035200d-08/
      data fpp( 6,10,1),fpp( 6,10,2)/-4.45730873d-08, 1.25044802d-08/
      data fpp( 6,11,1),fpp( 6,11,2)/-5.92660335d-07,-3.16214408d-08/
      data fpp( 6,12,1),fpp( 6,12,2)/ 8.36197360d-06, 4.01812831d-08/
      data fpp( 6,13,1),fpp( 6,13,2)/-5.35557986d-06,-1.21036917d-08/
      data fpp( 6,14,1),fpp( 6,14,2)/ 3.37154684d-07, 4.03348364d-09/
      data fpp( 6,15,1),fpp( 6,15,2)/ 5.52745435d-05,-1.84302429d-08/
      data fpp( 6,16,1),fpp( 6,16,2)/ 2.83938410d-05, 1.50874878d-08/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.00930501d-04, 2.34802918d-08/
      data fpp( 6,18,1),fpp( 6,18,2)/ 1.91556524d-04,-4.60086548d-08/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-9.32456726d-08/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00, 1.97594859d-07/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 2.75681688d-04, 1.01410281d-07/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 4.13677465d-04, 3.53164016d-07/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 2.15713010d-04,-6.12866346d-07/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 1.99193986d-04,-2.43498633d-07/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 1.22577732d-04, 2.38660878d-07/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 6.51994632d-05, 1.58855121d-07/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 3.71237215d-05,-1.21081360d-07/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 5.95195582d-05, 2.06070321d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/-9.53982760d-07,-1.37399924d-07/
      data fpp( 7,11,1),fpp( 7,11,2)/ 3.21337329d-07, 1.95293747d-08/
      data fpp( 7,12,1),fpp( 7,12,2)/-7.78415271d-07, 4.08242508d-09/
      data fpp( 7,13,1),fpp( 7,13,2)/ 6.58029941d-06, 9.74092497d-09/
      data fpp( 7,14,1),fpp( 7,14,2)/ 2.43896139d-06, 5.55387503d-09/
      data fpp( 7,15,1),fpp( 7,15,2)/-1.34807412d-05,-9.75642508d-09/
      data fpp( 7,16,1),fpp( 7,16,2)/-9.08032915d-06,-2.52817472d-09/
      data fpp( 7,17,1),fpp( 7,17,2)/-2.33693572d-05, 9.06912396d-09/
      data fpp( 7,18,1),fpp( 7,18,2)/-3.27542093d-05,-1.81483211d-08/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-3.42758394d-08/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00, 2.11891818d-08/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-6.81016841d-05, 1.14216364d-08/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-8.69336526d-05, 6.32427271d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 1.08051099d-05,-1.51187272d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-1.36841235d-07,-2.26493639d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-2.07486889d-05,-1.24838172d-08/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.71500676d-05, 2.87846327d-08/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-6.30064379d-06, 3.89452864d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-3.20830860d-05,-7.47657782d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.08423482d-06, 4.71178263d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-1.08768182d-06,-1.95055270d-08/
      data fpp( 8,12,1),fpp( 8,12,2)/-1.06574098d-06, 1.11042816d-08/
      data fpp( 8,13,1),fpp( 8,13,2)/-3.14310829d-06,-2.71159945d-09/
      data fpp( 8,14,1),fpp( 8,14,2)/ 3.67453850d-06,-5.05788379d-09/
      data fpp( 8,15,1),fpp( 8,15,2)/ 1.10649519d-05, 1.69431346d-08/
      data fpp( 8,16,1),fpp( 8,16,2)/ 2.84406695d-06,-1.71146546d-08/
      data fpp( 8,17,1),fpp( 8,17,2)/ 1.33428209d-05, 1.19154839d-08/
      data fpp( 8,18,1),fpp( 8,18,2)/ 1.06443659d-05,-4.14728113d-09/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-7.32635944d-09/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-7.52214370d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.41250484d-05,-4.05571260d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 2.87371454d-05, 2.62499409d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 3.32655003d-06,-1.94426375d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 1.07533793d-05, 3.53206093d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.68170234d-05,-1.14397995d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 7.00080726d-06, 5.48385888d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-5.44114637d-06,-1.20914556d-07/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 1.19927858d-05, 1.06019635d-07/
      data fpp( 9,10,1),fpp( 9,10,2)/-8.62956532d-07,-6.85639824d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.20938995d-06, 3.38362951d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 5.41379206d-07,-1.45811978d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/ 8.92133761d-07, 3.40884963d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-6.33711540d-06,-1.49727873d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/-9.77906628d-06,-6.53973470d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 1.04061328d-07, 1.15762175d-07/
      data fpp( 9,17,1),fpp( 9,17,2)/-1.01419266d-05,-1.05451355d-07/
      data fpp( 9,18,1),fpp( 9,18,2)/-2.86325417d-06, 3.96432443d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 8.00783779d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 5.17119318d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-3.87850969d-06, 1.79761364d-08/
      data fpp(10, 3,1),fpp(10, 3,2)/-5.81492909d-06,-8.76164775d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.50868999d-06, 5.82897736d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/ 7.43324146d-07, 4.46573830d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-2.29940468d-06,-7.19193058d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-1.49316145d-06, 2.10198401d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 2.25229274d-07,-3.16005462d-09/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.24805716d-06, 3.30203784d-08/
      data fpp(10,10,1),fpp(10,10,2)/-3.32408694d-07,-4.49214589d-08/
      data fpp(10,11,1),fpp(10,11,2)/-8.98779906d-08, 2.06654571d-08/
      data fpp(10,12,1),fpp(10,12,2)/-1.97758412d-08,-2.94036960d-09/
      data fpp(10,13,1),fpp(10,13,2)/-1.85426752d-07, 5.49602129d-09/
      data fpp(10,14,1),fpp(10,14,2)/ 4.93923079d-07,-6.44371555d-09/
      data fpp(10,15,1),fpp(10,15,2)/ 9.91313257d-07,-7.21159093d-10/
      data fpp(10,16,1),fpp(10,16,2)/ 2.79687734d-07, 6.92835192d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 1.24488531d-06,-4.19224859d-09/
      data fpp(10,18,1),fpp(10,18,2)/ 6.88650835d-07, 4.44064246d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 7.42967877d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 1.99103227d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/ 5.21800484d-06, 1.08793546d-08/
      data fpp(11, 3,1),fpp(11, 3,2)/-7.31878546d-06,-4.27741066d-10/
      data fpp(11, 4,1),fpp(11, 4,2)/-4.67434500d-06,-3.76839033d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-9.87166207d-06,-9.09869763d-09/
      data fpp(11, 6,1),fpp(11, 6,2)/-1.50252977d-05, 1.55631808d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-3.79091927d-06,-8.15402576d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 9.39488536d-06, 5.05292221d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.93722142d-06,-8.45766308d-09/
      data fpp(11,10,1),fpp(11,10,2)/-3.86295653d-07, 7.77773009d-09/
      data fpp(11,11,1),fpp(11,11,2)/-1.07006100d-06,-3.45325730d-09/
      data fpp(11,12,1),fpp(11,12,2)/-1.66362079d-07,-5.96470089d-09/
      data fpp(11,13,1),fpp(11,13,2)/-3.97866239d-08, 1.23120609d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 6.57678846d-06,-9.68354259d-09/
      data fpp(11,15,1),fpp(11,15,2)/ 8.30559337d-06, 4.22210948d-09/
      data fpp(11,16,1),fpp(11,16,2)/-6.81093867d-07,-4.89532392d-12/
      data fpp(11,17,1),fpp(11,17,2)/ 8.32630734d-06,-1.80252818d-09/
      data fpp(11,18,1),fpp(11,18,2)/ 1.15674583d-07, 3.01500805d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 5.94249597d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 5.37423250d-04, 1.82404728d-04/
      data fpppp( 1, 3),fpppp( 1, 4)/ 2.61436845d-04,-1.17870827d-03/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.17681324d-03,-5.53888011d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/ 1.62825246d-04,-3.37954790d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 8.89399218d-06,-4.26690264d-07/
      data fpppp( 1,11),fpppp( 1,12)/-3.33120438d-06, 1.77323681d-05/
      data fpppp( 1,13),fpppp( 1,14)/ 1.16981219d-05,-3.94228727d-05/
      data fpppp( 1,15),fpppp( 1,16)/-2.49941611d-05,-5.99344422d-05/
      data fpppp( 1,17),fpppp( 1,18)/-4.49107619d-06, 2.41552631d-04/
      data fpppp( 1,19) /             4.33139097d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 2.48033153d-04, 8.55388380d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 9.88334799d-05,-5.03133228d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 5.59623839d-04,-3.02316317d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 9.14445334d-05,-1.82376078d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 5.28325296d-06,-3.54202891d-07/
      data fpppp( 2,11),fpppp( 2,12)/-1.88009489d-06, 9.99286180d-06/
      data fpppp( 2,13),fpppp( 2,14)/ 8.03106797d-06,-2.01026993d-05/
      data fpppp( 2,15),fpppp( 2,16)/-1.54660105d-05,-2.54653398d-05/
      data fpppp( 2,17),fpppp( 2,18)/-4.54058180d-05, 1.70888845d-04/
      data fpppp( 2,19) /             3.19687753d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.70979682d-04,-6.42636424d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-5.12147097d-04, 9.81067680d-04/
      data fpppp( 3, 5),fpppp( 3, 6)/-5.63017439d-04, 9.88037456d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-1.76116084d-05, 4.96195707d-06/
      data fpppp( 3, 9),fpppp( 3,10)/ 1.52023759d-06,-1.37911137d-06/
      data fpppp( 3,11),fpppp( 3,12)/ 3.01995086d-07,-1.16604672d-06/
      data fpppp( 3,13),fpppp( 3,14)/ 4.75852086d-06, 9.70344305d-06/
      data fpppp( 3,15),fpppp( 3,16)/-1.64340342d-06, 2.99884546d-05/
      data fpppp( 3,17),fpppp( 3,18)/-1.80263458d-04, 1.53191364d-04/
      data fpppp( 3,19) /             3.31301404d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-1.36312931d-04,-2.87880566d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/ 4.46192568d-04,-3.64921933d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/ 8.34228323d-05,-2.34954819d-05/
      data fpppp( 4, 7),fpppp( 4, 8)/ 5.19625062d-06, 1.41559391d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-1.70981121d-06, 1.43686562d-06/
      data fpppp( 4,11),fpppp( 4,12)/-7.43466319d-08,-3.57047479d-07/
      data fpppp( 4,13),fpppp( 4,14)/ 1.05200051d-07, 2.80869382d-06/
      data fpppp( 4,15),fpppp( 4,16)/-2.85139394d-06, 1.19958272d-05/
      data fpppp( 4,17),fpppp( 4,18)/ 1.36790462d-05,-4.98881837d-05/
      data fpppp( 4,19) /            -9.48124283d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 1.14526935d-05,-1.26043082d-05/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.92019754d-04, 1.72319354d-04/
      data fpppp( 5, 5),fpppp( 5, 6)/-3.67769132d-05, 1.92297255d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/ 1.59046692d-06, 1.65832684d-07/
      data fpppp( 5, 9),fpppp( 5,10)/ 6.08884644d-07,-1.88260163d-08/
      data fpppp( 5,11),fpppp( 5,12)/-2.50986373d-07, 5.11822788d-07/
      data fpppp( 5,13),fpppp( 5,14)/-3.04878452d-08,-1.75193736d-06/
      data fpppp( 5,15),fpppp( 5,16)/ 4.42222214d-06,-6.85101614d-06/
      data fpppp( 5,17),fpppp( 5,18)/-4.40815829d-06, 1.41887518d-05/
      data fpppp( 5,19) /             2.81702184d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.16137498d-05, 9.70210349d-06/
      data fpppp( 6, 3),fpppp( 6, 4)/ 3.34371982d-05,-3.25913005d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/ 2.01813485d-06,-4.79448925d-07/
      data fpppp( 6, 7),fpppp( 6, 8)/ 2.07113011d-06,-1.08077664d-07/
      data fpppp( 6, 9),fpppp( 6,10)/ 1.49283630d-06,-1.50746317d-06/
      data fpppp( 6,11),fpppp( 6,12)/ 6.52934911d-07,-5.34113203d-07/
      data fpppp( 6,13),fpppp( 6,14)/ 1.23186659d-07, 1.20598385d-06/
      data fpppp( 6,15),fpppp( 6,16)/-1.99244279d-06, 1.85470181d-06/
      data fpppp( 6,17),fpppp( 6,18)/ 5.38677325d-07,-2.92404938d-06/
      data fpppp( 6,19) /            -5.77343259d-06 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 2.14658473d-07,-5.94407837d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/-6.09818176d-06, 4.82952095d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/-2.33317628d-06, 8.97350422d-07/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.01946335d-07, 1.26858656d-06/
      data fpppp( 7, 9),fpppp( 7,10)/-1.94410520d-06, 1.53567159d-06/
      data fpppp( 7,11),fpppp( 7,12)/-4.93649504d-07, 2.96422064d-07/
      data fpppp( 7,13),fpppp( 7,14)/-1.84530716d-07,-2.48302361d-07/
      data fpppp( 7,15),fpppp( 7,16)/ 4.71038285d-07,-4.16643898d-07/
      data fpppp( 7,17),fpppp( 7,18)/ 7.41708996d-08, 4.14210859d-07/
      data fpppp( 7,19) /             7.97329344d-07 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.13233152d-08, 1.74717784d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.27863511d-06,-2.29501435d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 3.80579489d-07, 1.92502605d-07/
      data fpppp( 8, 7),fpppp( 8, 8)/ 3.02038227d-07,-9.65607358d-07/
      data fpppp( 8, 9),fpppp( 8,10)/ 1.36247924d-06,-9.47323838d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 3.06461861d-07,-1.46892156d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 1.55148273d-07, 5.99999096d-08/
      data fpppp( 8,15),fpppp( 8,16)/-3.60781917d-07, 4.46449860d-07/
      data fpppp( 8,17),fpppp( 8,18)/-3.01839188d-07,-3.09256514d-08/
      data fpppp( 8,19) /            -5.12128531d-08 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 3.50605140d-07, 1.23595995d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-8.15766206d-07, 7.38107285d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.66417458d-07,-1.54228562d-07/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.69459907d-07, 6.74523940d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-7.36082704d-07, 4.52426407d-07/
      data fpppp( 9,11),fpppp( 9,12)/-1.77937596d-07, 9.49025421d-08/
      data fpppp( 9,13),fpppp( 9,14)/-1.40546655d-07, 1.24838542d-08/
      data fpppp( 9,15),fpppp( 9,16)/ 3.17849134d-07,-4.84375681d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 4.11906660d-07,-1.11771341d-07/
      data fpppp( 9,19) /            -2.29746388d-07 /
      data fpppp(10, 1),fpppp(10, 2)/-4.10333203d-08,-6.02961446d-09/
      data fpppp(10, 3),fpppp(10, 4)/ 1.81677196d-07,-1.65076859d-07/
      data fpppp(10, 5),fpppp(10, 6)/-6.70885526d-09, 5.52705015d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 1.65651721d-08,-6.68023398d-08/
      data fpppp(10, 9),fpppp(10,10)/ 5.91435575d-08,-2.64357964d-08/
      data fpppp(10,11),fpppp(10,12)/ 6.21256253d-09,-8.76016691d-09/
      data fpppp(10,13),fpppp(10,14)/ 1.46829215d-08, 7.28525401d-10/
      data fpppp(10,15),fpppp(10,16)/-2.85146023d-08, 4.07889420d-08/
      data fpppp(10,17),fpppp(10,18)/-3.40317797d-08, 4.05225332d-09/
      data fpppp(10,19) /             9.87778498d-09 /
      data fpppp(11, 1),fpppp(11, 2)/-4.80705129d-07,-2.29954114d-07/
      data fpppp(11, 3),fpppp(11, 4)/ 3.35233876d-07,-2.00107545d-07/
      data fpppp(11, 5),fpppp(11, 6)/-5.30914775d-09, 2.23965025d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 9.27298846d-08,-4.77798989d-07/
      data fpppp(11, 9),fpppp(11,10)/ 3.47391386d-07,-1.38784601d-07/
      data fpppp(11,11),fpppp(11,12)/ 7.36655525d-08,-6.06297520d-08/
      data fpppp(11,13),fpppp(11,14)/ 1.22226047d-07,-3.88744591d-08/
      data fpppp(11,15),fpppp(11,16)/-2.59994421d-07, 4.35922614d-07/
      data fpppp(11,17),fpppp(11,18)/-4.04050729d-07, 1.47198265d-07/
      data fpppp(11,19) /             3.00955162d-07 /

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
      subroutine c4_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  1.07440000d-04 /
      data f( 1, 3),f( 1, 4) / -8.99700000d-04 , -2.15331000d-03 /
      data f( 1, 5),f( 1, 6) / -3.70300000d-04 ,  5.83880000d-04 /
      data f( 1, 7),f( 1, 8) /  1.16160000d-04 ,  1.84900000d-05 /
      data f( 1, 9),f( 1,10) /  1.98000000d-06 , -5.00000000d-08 /
      data f( 1,11),f( 1,12) / -5.00000000d-08 , -2.45000000d-06 /
      data f( 1,13),f( 1,14) /  4.07000000d-05 ,  2.03860000d-04 /
      data f( 1,15),f( 1,16) /  3.92360000d-04 ,  5.50840000d-04 /
      data f( 1,17),f( 1,18) /  3.65750000d-04 , -2.90895000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  5.45080000d-04 /
      data f( 2, 3),f( 2, 4) / -3.67900000d-04 , -1.18505000d-03 /
      data f( 2, 5),f( 2, 6) /  2.95810000d-04 ,  2.38130000d-04 /
      data f( 2, 7),f( 2, 8) /  5.56700000d-05 ,  1.08500000d-05 /
      data f( 2, 9),f( 2,10) /  1.05000000d-06 , -2.00000000d-08 /
      data f( 2,11),f( 2,12) / -4.00000000d-08 , -7.70000000d-07 /
      data f( 2,13),f( 2,14) / -1.16000000d-06 ,  2.27400000d-05 /
      data f( 2,15),f( 2,16) /  1.01640000d-04 ,  2.56140000d-04 /
      data f( 2,17),f( 2,18) /  4.38910000d-04 ,  6.70320000d-04 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -7.03000000d-06 /
      data f( 3, 3),f( 3, 4) / -1.87800000d-05 ,  2.11800000d-04 /
      data f( 3, 5),f( 3, 6) /  1.90140000d-04 ,  9.92200000d-05 /
      data f( 3, 7),f( 3, 8) /  2.98400000d-05 ,  5.92000000d-06 /
      data f( 3, 9),f( 3,10) /  1.48000000d-06 , -5.00000000d-08 /
      data f( 3,11),f( 3,12) / -9.00000000d-08 ,  3.00000000d-07 /
      data f( 3,13),f( 3,14) / -1.08000000d-06 , -1.58000000d-06 /
      data f( 3,15),f( 3,16) /  7.14000000d-06 ,  3.64200000d-05 /
      data f( 3,17),f( 3,18) /  1.01510000d-04 ,  6.77300000d-05 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  1.50850000d-04 /
      data f( 4, 3),f( 4, 4) /  3.58080000d-04 ,  1.63750000d-04 /
      data f( 4, 5),f( 4, 6) /  8.55600000d-05 ,  4.90100000d-05 /
      data f( 4, 7),f( 4, 8) /  1.50800000d-05 ,  9.70000000d-07 /
      data f( 4, 9),f( 4,10) / -1.95000000d-06 , -4.00000000d-08 /
      data f( 4,11),f( 4,12) / -3.00000000d-08 , -4.50000000d-07 /
      data f( 4,13),f( 4,14) /  6.40000000d-07 , -2.80000000d-07 /
      data f( 4,15),f( 4,16) / -1.27000000d-06 , -6.50000000d-07 /
      data f( 4,17),f( 4,18) / -1.59000000d-06 , -1.08000000d-06 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  2.48000000d-06 /
      data f( 5, 3),f( 5, 4) /  3.82200000d-05 ,  4.51300000d-05 /
      data f( 5, 5),f( 5, 6) /  4.12500000d-05 ,  2.21000000d-05 /
      data f( 5, 7),f( 5, 8) /  6.13000000d-06 , -5.00000000d-07 /
      data f( 5, 9),f( 5,10) / -1.07000000d-06 ,  3.00000000d-08 /
      data f( 5,11),f( 5,12) / -3.00000000d-08 , -3.60000000d-07 /
      data f( 5,13),f( 5,14) / -1.20000000d-07 ,  3.60000000d-07 /
      data f( 5,15),f( 5,16) /  3.07000000d-06 ,  2.40000000d-07 /
      data f( 5,17),f( 5,18) / -9.00000000d-08 ,  6.00000000d-08 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  5.60000000d-07 /
      data f( 6, 3),f( 6, 4) /  6.16000000d-06 ,  1.55900000d-05 /
      data f( 6, 5),f( 6, 6) /  1.62400000d-05 ,  8.49000000d-06 /
      data f( 6, 7),f( 6, 8) /  2.99000000d-06 ,  1.38000000d-06 /
      data f( 6, 9),f( 6,10) /  2.00000000d-07 , -1.00000000d-08 /
      data f( 6,11),f( 6,12) / -3.00000000d-08 , -8.50000000d-07 /
      data f( 6,13),f( 6,14) / -4.00000000d-07 ,  3.00000000d-08 /
      data f( 6,15),f( 6,16) /  1.40000000d-07 , -3.90000000d-07 /
      data f( 6,17),f( 6,18) /  1.30000000d-07 ,  1.63000000d-06 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  4.00000000d-07 /
      data f( 7, 3),f( 7, 4) /  9.60000000d-07 ,  4.50000000d-06 /
      data f( 7, 5),f( 7, 6) /  5.13000000d-06 ,  9.60000000d-07 /
      data f( 7, 7),f( 7, 8) / -1.61000000d-06 , -1.18000000d-06 /
      data f( 7, 9),f( 7,10) / -3.21000000d-06 ,  0.00000000d+00 /
      data f( 7,11),f( 7,12) /  4.00000000d-08 , -7.40000000d-07 /
      data f( 7,13),f( 7,14) / -8.70000000d-07 , -5.20000000d-07 /
      data f( 7,15),f( 7,16) /  1.00000000d-07 ,  3.20000000d-07 /
      data f( 7,17),f( 7,18) /  4.00000000d-07 ,  7.70000000d-07 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  5.00000000d-08 /
      data f( 8, 3),f( 8, 4) /  8.00000000d-08 ,  2.30000000d-07 /
      data f( 8, 5),f( 8, 6) /  1.35000000d-06 ,  1.19000000d-06 /
      data f( 8, 7),f( 8, 8) / -1.40000000d-07 , -2.30000000d-07 /
      data f( 8, 9),f( 8,10) /  1.13000000d-06 ,  0.00000000d+00 /
      data f( 8,11),f( 8,12) /  3.50000000d-07 , -2.90000000d-07 /
      data f( 8,13),f( 8,14) / -3.40000000d-07 , -4.10000000d-07 /
      data f( 8,15),f( 8,16) / -5.50000000d-07 , -2.00000000d-07 /
      data f( 8,17),f( 8,18) / -2.90000000d-07 , -8.00000000d-08 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  1.14000000d-06 /
      data f( 9, 3),f( 9, 4) / -6.00000000d-08 , -5.60000000d-07 /
      data f( 9, 5),f( 9, 6) / -1.24000000d-06 , -8.40000000d-07 /
      data f( 9, 7),f( 9, 8) / -3.30000000d-07 ,  1.40000000d-06 /
      data f( 9, 9),f( 9,10) / -4.50000000d-07 ,  8.70000000d-07 /
      data f( 9,11),f( 9,12) / -2.90000000d-07 , -3.30000000d-07 /
      data f( 9,13),f( 9,14) / -3.20000000d-07 ,  8.80000000d-07 /
      data f( 9,15),f( 9,16) /  1.07000000d-06 , -5.00000000d-07 /
      data f( 9,17),f( 9,18) /  1.11000000d-06 , -1.60000000d-07 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  5.50000000d-07 /
      data f(10, 3),f(10, 4) /  1.43000000d-06 , -6.40000000d-07 /
      data f(10, 5),f(10, 6) / -5.30000000d-07 ,  1.16000000d-06 /
      data f(10, 7),f(10, 8) /  4.90000000d-07 , -2.20000000d-07 /
      data f(10, 9),f(10,10) / -3.40000000d-07 ,  7.50000000d-07 /
      data f(10,11),f(10,12) / -7.00000000d-08 , -2.60000000d-07 /
      data f(10,13),f(10,14) / -2.70000000d-07 , -1.70000000d-07 /
      data f(10,15),f(10,16) / -2.70000000d-07 , -4.20000000d-07 /
      data f(10,17),f(10,18) / -3.20000000d-07 , -2.80000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -6.00000000d-07 /
      data f(11, 3),f(11, 4) / -3.40000000d-07 , -1.00000000d-07 /
      data f(11, 5),f(11, 6) / -1.70000000d-07 , -1.07000000d-06 /
      data f(11, 7),f(11, 8) / -2.80000000d-07 ,  3.20000000d-07 /
      data f(11, 9),f(11,10) /  2.60000000d-07 , -1.10000000d-07 /
      data f(11,11),f(11,12) /  1.00000000d-07 ,  6.00000000d-08 /
      data f(11,13),f(11,14) / -3.60000000d-07 , -1.10000000d-07 /
      data f(11,15),f(11,16) / -2.20000000d-07 , -2.10000000d-07 /
      data f(11,17),f(11,18) / -1.90000000d-07 , -2.30000000d-07 /
      data f(11,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-2.16957483d-05/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-1.02348376d-02,-7.40810331d-06/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-1.32594339d-03,-1.55466384d-05/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 8.61877031d-03, 5.48064570d-05/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-6.06514054d-03,-2.14819894d-05/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 1.29363109d-03,-1.86082993d-05/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 2.30327967d-04, 1.06011865d-05/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 2.07070955d-05,-1.59344665d-06/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 2.36314459d-05, 6.42200109d-07/
      data fpp( 1,10,1),fpp( 1,10,2)/-6.40319036d-07,-1.06553790d-07/
      data fpp( 1,11,1),fpp( 1,11,2)/-8.53181183d-07,-9.41849491d-08/
      data fpp( 1,12,1),fpp( 1,12,2)/ 1.54093318d-06, 3.39293586d-07/
      data fpp( 1,13,1),fpp( 1,13,2)/ 3.23324323d-04, 1.47001060d-06/
      data fpp( 1,14,1),fpp( 1,14,2)/ 1.13566773d-03, 9.81263998d-07/
      data fpp( 1,15,1),fpp( 1,15,2)/ 1.22126531d-03,-3.87466660d-06/
      data fpp( 1,16,1),fpp( 1,16,2)/-8.47638722d-05, 1.27162024d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/-4.09347691d-03,-6.76043430d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-3.47144694d-02, 7.23245694d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 1.49325065d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-3.12235654d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-5.00524485d-03,-1.23724693d-05/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-1.22363322d-03,-6.77015756d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 2.89345938d-03, 4.52028995d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-3.31403892d-03,-3.61608404d-05/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 8.39937819d-04, 7.12806228d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 1.43984067d-04, 1.61791331d-07/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 1.29458090d-05, 4.83172396d-07/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 1.08171083d-05, 6.71908315d-09/
      data fpp( 2,10,1),fpp( 2,10,2)/-2.39361928d-07, 1.37512710d-08/
      data fpp( 2,11,1),fpp( 2,11,2)/-3.73637635d-07, 1.27583279d-09/
      data fpp( 2,12,1),fpp( 2,12,2)/-7.21866351d-07,-6.14546022d-08/
      data fpp( 2,13,1),fpp( 2,13,2)/ 1.79031354d-04, 2.64942576d-07/
      data fpp( 2,14,1),fpp( 2,14,2)/ 6.59704536d-04, 4.59084299d-07/
      data fpp( 2,15,1),fpp( 2,15,2)/ 7.93149376d-04, 1.19872023d-06/
      data fpp( 2,16,1),fpp( 2,16,2)/ 2.07927744d-04,-7.17965216d-07/
      data fpp( 2,17,1),fpp( 2,17,2)/-1.89864618d-03, 3.36934064d-06/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.84785011d-02,-9.84099732d-06/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00,-1.81091513d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00,-2.13270930d-06/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 6.50181697d-03,-8.17181410d-07/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 1.83615625d-03, 5.11823493d-06/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-9.90644783d-03,-5.11595832d-06/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 7.98576234d-04, 2.11198364d-07/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 3.10777632d-04, 1.15564867d-07/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 2.55757662d-05, 6.18942169d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-7.45033139d-06, 1.36266459d-07/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-3.42598789d-05, 4.79199511d-09/
      data fpp( 3,10,1),fpp( 3,10,2)/ 1.57766748d-07, 1.91655605d-08/
      data fpp( 3,11,1),fpp( 3,11,2)/ 9.07731722d-07, 7.94576294d-09/
      data fpp( 3,12,1),fpp( 3,12,2)/-1.32934678d-05,-2.51486123d-08/
      data fpp( 3,13,1),fpp( 3,13,2)/-3.28897403d-05,-1.35513139d-08/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.12858755d-05, 1.32153868d-07/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.15417185d-04, 3.81358420d-08/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.05257289d-03, 9.48902764d-07/
      data fpp( 3,17,1),fpp( 3,17,2)/ 1.83462164d-03,-1.68514690d-06/
      data fpp( 3,18,1),fpp( 3,18,2)/ 8.26383388d-03,-1.40515172d-07/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 2.10207586d-07/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 4.67245245d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-3.96226304d-03, 1.50529510d-06/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-5.45523180d-03,-7.31083286d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 2.05473194d-03, 3.64443633d-06/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.45893989d-04,-2.98512467d-07/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 4.57516522d-05, 4.80135379d-08/
      data fpp( 4, 7,1),fpp( 4, 7,2)/ 1.93928683d-05, 2.63658316d-07/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 1.63755166d-05, 8.65531993d-08/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 3.35824074d-05, 6.15288872d-08/
      data fpp( 4,10,1),fpp( 4,10,2)/ 5.68294937d-07,-4.28687479d-08/
      data fpp( 4,11,1),fpp( 4,11,2)/-6.17289252d-07,-4.05389558d-09/
      data fpp( 4,12,1),fpp( 4,12,2)/ 1.02157374d-05, 3.32843302d-08/
      data fpp( 4,13,1),fpp( 4,13,2)/-8.11239305d-06,-3.84834252d-08/
      data fpp( 4,14,1),fpp( 4,14,2)/ 3.18966109d-07, 4.93707923d-11/
      data fpp( 4,15,1),fpp( 4,15,2)/ 1.13418829d-05, 3.40859421d-08/
      data fpp( 4,16,1),fpp( 4,16,2)/-3.46193229d-05,-3.97931391d-08/
      data fpp( 4,17,1),fpp( 4,17,2)/ 1.83359626d-04, 3.14866144d-08/
      data fpp( 4,18,1),fpp( 4,18,2)/-1.76611443d-03, 8.46681607d-10/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-6.73340803d-10/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 9.17022322d-07/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 1.99723520d-03, 4.05555357d-07/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 3.26349094d-03,-5.43643749d-07/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-6.15992526d-06, 3.92196399d-08/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 6.43278112d-05,-2.60634810d-07/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 6.54157592d-05, 8.71196011d-08/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 3.62927605d-05, 1.02956406d-07/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 2.54682650d-05, 6.14547760d-08/
      data fpp( 5, 9,1),fpp( 5, 9,2)/ 3.37024920d-06, 1.48244902d-08/
      data fpp( 5,10,1),fpp( 5,10,2)/-9.90946496d-07,-2.05527369d-08/
      data fpp( 5,11,1),fpp( 5,11,2)/ 1.21425286d-07,-2.21354271d-09/
      data fpp( 5,12,1),fpp( 5,12,2)/-7.40948196d-06, 1.32069077d-08/
      data fpp( 5,13,1),fpp( 5,13,2)/ 5.81931252d-06,-1.64140881d-08/
      data fpp( 5,14,1),fpp( 5,14,2)/-5.82998895d-06, 6.68494448d-08/
      data fpp( 5,15,1),fpp( 5,15,2)/-5.47847171d-05,-1.17183691d-07/
      data fpp( 5,16,1),fpp( 5,16,2)/-3.05560292d-06, 6.94853192d-08/
      data fpp( 5,17,1),fpp( 5,17,2)/-5.76601419d-05,-1.07575859d-08/
      data fpp( 5,18,1),fpp( 5,18,2)/ 4.79423822d-04, 2.34502456d-09/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-1.12225123d-08/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 6.83041543d-08/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-5.11877746d-04, 3.87916913d-08/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-6.91531954d-04, 7.89290803d-08/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 1.07827762d-04,-1.24708013d-07/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 5.99947666d-05,-1.06897030d-07/
      data fpp( 6, 6,1),fpp( 6, 6,2)/ 1.17853109d-05, 4.82961330d-08/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-2.51239102d-05, 4.87124982d-08/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-3.78485765d-05,-9.74612559d-09/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-3.77034042d-05, 1.60720042d-08/
      data fpp( 6,10,1),fpp( 6,10,2)/ 7.55491047d-07, 3.65810869d-09/
      data fpp( 6,11,1),fpp( 6,11,2)/ 1.31588109d-07,-1.93044390d-08/
      data fpp( 6,12,1),fpp( 6,12,2)/ 5.50219041d-06, 2.55596472d-08/
      data fpp( 6,13,1),fpp( 6,13,2)/-3.64485703d-06,-6.73414981d-09/
      data fpp( 6,14,1),fpp( 6,14,2)/-2.79010323d-07, 1.76952060d-10/
      data fpp( 6,15,1),fpp( 6,15,2)/ 3.33169853d-05,-1.31736584d-08/
      data fpp( 6,16,1),fpp( 6,16,2)/ 1.03617346d-05, 1.41176816d-08/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.65609420d-05, 1.97029319d-08/
      data fpp( 6,18,1),fpp( 6,18,2)/-1.41260863d-04,-3.41294091d-08/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-7.09852954d-08/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-2.28067416d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 9.25157859d-05,-5.98651678d-09/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 1.47276878d-04, 5.63528087d-08/
      data fpp( 7, 4,1),fpp( 7, 4,2)/ 1.76488760d-05,-4.06247181d-08/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 2.92931226d-05,-6.84539364d-08/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 3.33629972d-05, 2.64404635d-08/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 2.91628805d-05, 5.86920823d-08/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 1.93660411d-05,-8.12087928d-08/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 3.51233678d-05, 1.18543089d-07/
      data fpp( 7,10,1),fpp( 7,10,2)/-8.31017692d-07,-7.85635623d-08/
      data fpp( 7,11,1),fpp( 7,11,2)/ 1.03222228d-06, 5.51116042d-09/
      data fpp( 7,12,1),fpp( 7,12,2)/-1.99279689d-07, 7.31892061d-09/
      data fpp( 7,13,1),fpp( 7,13,2)/ 4.20011559d-06, 4.21315713d-09/
      data fpp( 7,14,1),fpp( 7,14,2)/ 1.66603024d-06, 4.62845087d-09/
      data fpp( 7,15,1),fpp( 7,15,2)/-9.12322411d-06,-6.52696059d-09/
      data fpp( 7,16,1),fpp( 7,16,2)/-6.23133545d-06,-2.52060850d-09/
      data fpp( 7,17,1),fpp( 7,17,2)/-7.38362598d-06, 8.20939460d-09/
      data fpp( 7,18,1),fpp( 7,18,2)/ 2.72996296d-05,-1.29169699d-08/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-2.49415151d-08/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-2.28972498d-09/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-2.17884850d-05, 1.17944996d-09/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-3.89446575d-05,-3.62807484d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 5.99490939d-07, 2.05328494d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-7.23675112d-06,-2.03033228d-08/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-1.42416471d-05,-1.61195581d-08/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-1.09066863d-05, 1.45815552d-08/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-2.75383493d-06, 3.21933373d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-1.95584012d-05,-5.63549045d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.99530755d-06, 4.38262807d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-2.14246089d-06,-3.01502182d-08/
      data fpp( 8,12,1),fpp( 8,12,2)/-7.73256140d-07, 1.73745921d-08/
      data fpp( 8,13,1),fpp( 8,13,2)/-1.95791826d-06,-3.94815033d-09/
      data fpp( 8,14,1),fpp( 8,14,2)/ 2.40141444d-06,-2.78199081d-09/
      data fpp( 8,15,1),fpp( 8,15,2)/ 7.29117969d-06, 1.08761136d-08/
      data fpp( 8,16,1),fpp( 8,16,2)/ 1.87313907d-06,-1.13224634d-08/
      data fpp( 8,17,1),fpp( 8,17,2)/ 6.49040695d-06, 8.01374015d-09/
      data fpp( 8,18,1),fpp( 8,18,2)/-6.04845721d-06,-2.73249719d-09/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-4.88375141d-09/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-5.18328582d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 3.27815394d-06,-2.73342836d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 1.29417520d-05, 2.07699925d-08/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 8.33160279d-07,-1.37456865d-08/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 6.79388187d-06, 2.34127533d-08/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 1.00435914d-05,-1.51053269d-08/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 4.50386482d-06, 4.36085543d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-4.27070133d-06,-8.61288902d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 7.59023717d-06, 8.61070066d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.93021251d-06,-6.80991364d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.83762129d-06, 3.74895388d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 3.52304247d-07,-1.46590187d-08/
      data fpp( 9,13,1),fpp( 9,13,2)/ 5.71557436d-07, 2.41465359d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-4.19168801d-06,-1.05271251d-08/
      data fpp( 9,15,1),fpp( 9,15,2)/-6.42149466d-06,-4.26380355d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 5.87791934d-08, 7.54792673d-08/
      data fpp( 9,17,1),fpp( 9,17,2)/-6.03800183d-06,-6.84790335d-08/
      data fpp( 9,18,1),fpp( 9,18,2)/ 1.51419927d-06, 2.56368667d-08/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 5.17315666d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 3.21312339d-08/
      data fpp(10, 2,1),fpp(10, 2,2)/-1.40413079d-06, 1.12375322d-08/
      data fpp(10, 3,1),fpp(10, 3,2)/-3.04235040d-06,-5.72813627d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 3.27867944d-07, 4.08879186d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.38776375d-07, 2.45296884d-08/
      data fpp(10, 6,1),fpp(10, 6,2)/-1.75271827d-06,-4.42066720d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-1.04877296d-06, 1.06969996d-08/
      data fpp(10, 8,1),fpp(10, 8,2)/ 3.36640267d-07,-9.81326568d-10/
      data fpp(10, 9,1),fpp(10, 9,2)/-6.62547434d-07, 2.86283066d-08/
      data fpp(10,10,1),fpp(10,10,2)/-2.14457497d-07,-4.09318999d-08/
      data fpp(10,11,1),fpp(10,11,2)/-4.80242575d-08, 2.04992932d-08/
      data fpp(10,12,1),fpp(10,12,2)/ 2.40391505d-08,-3.26527268d-09/
      data fpp(10,13,1),fpp(10,13,2)/-1.48311487d-07, 3.36179758d-09/
      data fpp(10,14,1),fpp(10,14,2)/ 3.25337602d-07,-3.58191763d-09/
      data fpp(10,15,1),fpp(10,15,2)/ 6.34798931d-07,-1.03412704d-09/
      data fpp(10,16,1),fpp(10,16,2)/ 1.71744161d-07, 4.71842581d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 6.81600366d-07,-2.83957618d-09/
      data fpp(10,18,1),fpp(10,18,2)/-2.48339853d-07, 3.03987891d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 5.08006055d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 1.70016272d-08/
      data fpp(11, 2,1),fpp(11, 2,2)/ 2.61831539d-06, 9.39674553d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-4.46882480d-06,-2.98860936d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/-3.50183972d-07, 1.35769193d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-4.57061181d-06,-2.10421583d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-9.10864086d-06, 3.30109414d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-2.72061352d-06,-9.60160735d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 6.79542987d-06,-6.00451203d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-1.23747628d-06,-5.98034453d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 6.78478749d-07, 1.13258902d-08/
      data fpp(11,11,1),fpp(11,11,2)/-1.17973787d-06,-4.52321612d-09/
      data fpp(11,12,1),fpp(11,12,2)/ 2.17304247d-08,-8.23302569d-09/
      data fpp(11,13,1),fpp(11,13,2)/-1.25844256d-07, 1.46553189d-08/
      data fpp(11,14,1),fpp(11,14,2)/ 4.35983120d-06,-1.01882498d-08/
      data fpp(11,15,1),fpp(11,15,2)/ 5.40135053d-06, 4.49768018d-09/
      data fpp(11,16,1),fpp(11,16,2)/-4.69622081d-07,-6.02470977d-10/
      data fpp(11,17,1),fpp(11,17,2)/ 5.45919982d-06,-1.48779628d-09/
      data fpp(11,18,1),fpp(11,18,2)/ 4.22919927d-07, 2.95365608d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 5.87317196d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 3.86317029d-04, 1.63836140d-04/
      data fpppp( 1, 3),fpppp( 1, 4)/ 1.06962318d-04,-5.29536239d-04/
      data fpppp( 1, 5),fpppp( 1, 6)/ 5.33465167d-04,-2.81763480d-04/
      data fpppp( 1, 7),fpppp( 1, 8)/ 8.82642678d-05,-2.00726562d-05/
      data fpppp( 1, 9),fpppp( 1,10)/ 4.77907013d-06,-6.75391276d-07/
      data fpppp( 1,11),fpppp( 1,12)/-6.33970863d-07, 3.36769332d-06/
      data fpppp( 1,13),fpppp( 1,14)/ 6.32655411d-06, 7.59691415d-07/
      data fpppp( 1,15),fpppp( 1,16)/-5.29700695d-05, 1.27622981d-04/
      data fpppp( 1,17),fpppp( 1,18)/-6.19682885d-04, 7.54371788d-04/
      data fpppp( 1,19) /             1.52232345d-03 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 1.78024959d-04, 7.65832877d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 4.28532796d-05,-2.27867548d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 2.49141460d-04,-1.47009788d-04/
      data fpppp( 2, 7),fpppp( 2, 8)/ 4.79018615d-05,-1.07027287d-05/
      data fpppp( 2, 9),fpppp( 2,10)/ 2.64362665d-06,-4.07444096d-07/
      data fpppp( 2,11),fpppp( 2,12)/-3.58518599d-07, 1.82868131d-06/
      data fpppp( 2,13),fpppp( 2,14)/ 3.84988031d-06, 8.26995075d-07/
      data fpppp( 2,15),fpppp( 2,16)/-2.79915611d-05, 6.80192611d-05/
      data fpppp( 2,17),fpppp( 2,18)/-3.35366621d-04, 4.05050364d-04/
      data fpppp( 2,19) /             8.18666529d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.71568709d-04,-7.36975989d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-2.03689557d-04, 4.63839224d-04/
      data fpppp( 3, 5),fpppp( 3, 6)/-3.04809650d-04, 8.38300161d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-1.83546104d-05, 4.71897140d-06/
      data fpppp( 3, 9),fpppp( 3,10)/-1.48282238d-07,-4.52210856d-07/
      data fpppp( 3,11),fpppp( 3,12)/-6.29351810d-08,-1.93118288d-07/
      data fpppp( 3,13),fpppp( 3,14)/ 5.11703950d-07, 6.18310733d-07/
      data fpppp( 3,15),fpppp( 3,16)/ 1.53210049d-05,-3.72751713d-05/
      data fpppp( 3,17),fpppp( 3,18)/ 1.36473263d-04,-1.69788069d-04/
      data fpppp( 3,19) /            -3.38903754d-04 /
      data fpppp( 4, 1),fpppp( 4, 2)/-2.66067855d-05,-3.38036440d-06/
      data fpppp( 4, 3),fpppp( 4, 4)/ 1.88285900d-04,-2.09587288d-04/
      data fpppp( 4, 5),fpppp( 4, 6)/ 8.49351488d-05,-2.16315707d-05/
      data fpppp( 4, 7),fpppp( 4, 8)/ 6.01814733d-06,-1.04053265d-06/
      data fpppp( 4, 9),fpppp( 4,10)/-6.42562180d-07, 5.97521171d-07/
      data fpppp( 4,11),fpppp( 4,12)/ 1.62189195d-07,-5.25161299d-07/
      data fpppp( 4,13),fpppp( 4,14)/ 1.88786569d-07, 1.37558440d-06/
      data fpppp( 4,15),fpppp( 4,16)/-5.53563071d-06, 1.73478911d-05/
      data fpppp( 4,17),fpppp( 4,18)/-4.80195243d-05, 4.46830260d-05/
      data fpppp( 4,19) /             9.22227290d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/ 2.45416465d-05, 5.08586634d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-8.87438792d-05, 7.77352541d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/-2.17888214d-05, 5.25604423d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/-1.04801233d-06, 3.39152720d-08/
      data fpppp( 5, 9),fpppp( 5,10)/ 2.35940024d-07, 8.65338378d-08/
      data fpppp( 5,11),fpppp( 5,12)/-2.53661326d-07, 4.09514724d-07/
      data fpppp( 5,13),fpppp( 5,14)/-1.38815465d-07,-1.34693862d-06/
      data fpppp( 5,15),fpppp( 5,16)/ 3.28824435d-06,-5.76500823d-06/
      data fpppp( 5,17),fpppp( 5,18)/ 1.33917694d-05,-1.23007591d-05/
      data fpppp( 5,19) /            -2.51792001d-05 /
      data fpppp( 6, 1),fpppp( 6, 2)/-1.79871330d-06, 6.28324977d-07/
      data fpppp( 6, 3),fpppp( 6, 4)/ 1.92188256d-05,-1.87627921d-05/
      data fpppp( 6, 5),fpppp( 6, 6)/ 5.00077985d-06,-1.26291496d-06/
      data fpppp( 6, 7),fpppp( 6, 8)/ 7.28894042d-07,-2.01587920d-07/
      data fpppp( 6, 9),fpppp( 6,10)/ 8.49647949d-07,-8.98180493d-07/
      data fpppp( 6,11),fpppp( 6,12)/ 3.98106129d-07,-3.34573708d-07/
      data fpppp( 6,13),fpppp( 6,14)/ 6.91297193d-08, 8.08828479d-07/
      data fpppp( 6,15),fpppp( 6,16)/-1.49063470d-06, 1.76063555d-06/
      data fpppp( 6,17),fpppp( 6,18)/-3.80264002d-06, 3.60866379d-06/
      data fpppp( 6,19) /             7.31294493d-06 /
      data fpppp( 7, 1),fpppp( 7, 2)/ 8.35781234d-07, 1.28484736d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/-3.61500180d-06, 3.26817680d-06/
      data fpppp( 7, 5),fpppp( 7, 6)/-9.81370470d-07, 2.02842756d-07/
      data fpppp( 7, 7),fpppp( 7, 8)/-3.26200037d-07, 7.66154033d-07/
      data fpppp( 7, 9),fpppp( 7,10)/-1.20516613d-06, 9.51807737d-07/
      data fpppp( 7,11),fpppp( 7,12)/-3.33007294d-07, 1.94536923d-07/
      data fpppp( 7,13),fpppp( 7,14)/-1.07286565d-07,-1.81399501d-07/
      data fpppp( 7,15),fpppp( 7,16)/ 3.37574431d-07,-3.48029640d-07/
      data fpppp( 7,17),fpppp( 7,18)/ 8.11893380d-07,-7.49411114d-07/
      data fpppp( 7,19) /            -1.53322203d-06 /
      data fpppp( 8, 1),fpppp( 8, 2)/-3.91554201d-07,-1.19282399d-07/
      data fpppp( 8, 3),fpppp( 8, 4)/ 1.14662254d-06,-1.06518850d-06/
      data fpppp( 8, 5),fpppp( 8, 6)/ 2.71308038d-07, 2.98371137d-08/
      data fpppp( 8, 7),fpppp( 8, 8)/ 2.29734918d-07,-6.59703352d-07/
      data fpppp( 8, 9),fpppp( 8,10)/ 9.11633428d-07,-6.85333854d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 2.88213353d-07,-1.37101168d-07/
      data fpppp( 8,13),fpppp( 8,14)/ 1.06959307d-07, 4.19036289d-08/
      data fpppp( 8,15),fpppp( 8,16)/-2.42747870d-07, 3.10619497d-07/
      data fpppp( 8,17),fpppp( 8,18)/-3.97611607d-07, 2.50459010d-07/
      data fpppp( 8,19) /             5.11014850d-07 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 3.10074671d-07, 1.34566658d-07/
      data fpppp( 9, 3),fpppp( 9, 4)/-4.65214659d-07, 4.19960591d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-1.30468906d-07,-6.07456942d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/-1.53914478d-07, 4.82313229d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-5.37208160d-07, 3.83636119d-07/
      data fpppp( 9,11),fpppp( 9,12)/-2.00039306d-07, 1.01332055d-07/
      data fpppp( 9,13),fpppp( 9,14)/-1.03014701d-07, 1.17768302d-08/
      data fpppp( 9,15),fpppp( 9,16)/ 2.07913708d-07,-3.20826834d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 3.20770335d-07,-1.43315579d-07/
      data fpppp( 9,19) /            -2.91492041d-07 /
      data fpppp(10, 1),fpppp(10, 2)/-4.77011305d-08,-1.64709391d-08/
      data fpppp(10, 3),fpppp(10, 4)/ 9.95395575d-08,-8.11810139d-08/
      data fpppp(10, 5),fpppp(10, 6)/-5.02726154d-09, 3.24522054d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 1.42916720d-08,-4.87308180d-08/
      data fpppp(10, 9),fpppp(10,10)/ 3.75555441d-08,-1.46547004d-08/
      data fpppp(10,11),fpppp(10,12)/ 4.16385554d-09,-7.66291172d-09/
      data fpppp(10,13),fpppp(10,14)/ 1.18229486d-08,-8.68899007d-10/
      data fpppp(10,15),fpppp(10,16)/-1.81986182d-08, 2.73124059d-08/
      data fpppp(10,17),fpppp(10,18)/-3.26763467d-08, 1.70051957d-08/
      data fpppp(10,19) /             3.53523682d-08 /
      data fpppp(11, 1),fpppp(11, 2)/-2.86520630d-07,-1.36347139d-07/
      data fpppp(11, 3),fpppp(11, 4)/ 2.49581852d-07,-1.89633409d-07/
      data fpppp(11, 5),fpppp(11, 6)/ 8.60766246d-09, 1.36146686d-07/
      data fpppp(11, 7),fpppp(11, 8)/ 1.02368976d-07,-3.57941630d-07/
      data fpppp(11, 9),fpppp(11,10)/ 2.76460571d-07,-1.50968982d-07/
      data fpppp(11,11),fpppp(11,12)/ 1.00965059d-07,-6.93101599d-08/
      data fpppp(11,13),fpppp(11,14)/ 9.53330017d-08,-3.40268387d-08/
      data fpppp(11,15),fpppp(11,16)/-1.65875014d-07, 2.82777379d-07/
      data fpppp(11,17),fpppp(11,18)/-2.57246829d-07, 8.83038315d-08/
      data fpppp(11,19) /             1.80833101d-07 /

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
      subroutine c5_spl_c3h3_h(xi,yi,ix,iy,delxi,delyi,xix,xixp1,
     x                    yiy,yiyp1,fi)
      implicit real*8 (a-h,o-z)
 
      dimension fpp(11,19,2),f(11,19),fpppp(11,19)
      dimension c(4,4),px(4),py(4)
 
      data f( 1, 1),f( 1, 2) /  0.00000000d+00 ,  7.26400000d-05 /
      data f( 1, 3),f( 1, 4) / -1.92700000d-04 , -3.33560000d-04 /
      data f( 1, 5),f( 1, 6) / -6.48000000d-05 ,  1.12930000d-04 /
      data f( 1, 7),f( 1, 8) /  1.51700000d-05 ,  2.06000000d-06 /
      data f( 1, 9),f( 1,10) /  4.80000000d-07 , -3.00000000d-08 /
      data f( 1,11),f( 1,12) /  4.00000000d-08 , -5.50000000d-07 /
      data f( 1,13),f( 1,14) /  3.76000000d-06 ,  3.03800000d-05 /
      data f( 1,15),f( 1,16) /  6.64700000d-05 ,  1.02380000d-04 /
      data f( 1,17),f( 1,18) /  6.96600000d-05 , -2.96325000d-03 /
      data f( 1,19) /           0.00000000d+00 /
      data f( 2, 1),f( 2, 2) /  0.00000000d+00 ,  3.10600000d-04 /
      data f( 2, 3),f( 2, 4) /  6.91300000d-05 , -2.23960000d-04 /
      data f( 2, 5),f( 2, 6) /  9.63000000d-05 ,  3.99700000d-05 /
      data f( 2, 7),f( 2, 8) /  5.86000000d-06 ,  8.60000000d-07 /
      data f( 2, 9),f( 2,10) /  5.00000000d-08 ,  2.00000000d-08 /
      data f( 2,11),f( 2,12) /  1.00000000d-08 , -1.80000000d-07 /
      data f( 2,13),f( 2,14) / -8.70000000d-07 ,  3.80000000d-07 /
      data f( 2,15),f( 2,16) /  1.04900000d-05 ,  3.83600000d-05 /
      data f( 2,17),f( 2,18) /  8.37200000d-05 , -1.77110000d-04 /
      data f( 2,19) /           0.00000000d+00 /
      data f( 3, 1),f( 3, 2) /  0.00000000d+00 , -6.23900000d-05 /
      data f( 3, 3),f( 3, 4) /  3.63000000d-06 ,  9.14000000d-05 /
      data f( 3, 5),f( 3, 6) /  4.24500000d-05 ,  1.34500000d-05 /
      data f( 3, 7),f( 3, 8) /  2.60000000d-06 ,  5.50000000d-07 /
      data f( 3, 9),f( 3,10) /  4.50000000d-07 , -3.00000000d-08 /
      data f( 3,11),f( 3,12) / -1.00000000d-08 ,  1.60000000d-07 /
      data f( 3,13),f( 3,14) / -2.50000000d-07 , -9.40000000d-07 /
      data f( 3,15),f( 3,16) / -1.40000000d-06 ,  1.54000000d-06 /
      data f( 3,17),f( 3,18) /  1.38000000d-05 ,  9.76000000d-06 /
      data f( 3,19) /           0.00000000d+00 /
      data f( 4, 1),f( 4, 2) /  0.00000000d+00 ,  2.48200000d-05 /
      data f( 4, 3),f( 4, 4) /  7.26100000d-05 ,  3.67300000d-05 /
      data f( 4, 5),f( 4, 6) /  1.41400000d-05 ,  5.43000000d-06 /
      data f( 4, 7),f( 4, 8) /  1.26000000d-06 , -2.50000000d-07 /
      data f( 4, 9),f( 4,10) / -6.90000000d-07 , -2.00000000d-08 /
      data f( 4,11),f( 4,12) / -1.00000000d-08 , -1.60000000d-07 /
      data f( 4,13),f( 4,14) /  3.10000000d-07 , -3.00000000d-08 /
      data f( 4,15),f( 4,16) / -5.60000000d-07 , -9.80000000d-07 /
      data f( 4,17),f( 4,18) / -9.50000000d-07 , -3.70000000d-07 /
      data f( 4,19) /           0.00000000d+00 /
      data f( 5, 1),f( 5, 2) /  0.00000000d+00 ,  1.40000000d-07 /
      data f( 5, 3),f( 5, 4) /  5.16000000d-06 ,  6.95000000d-06 /
      data f( 5, 5),f( 5, 6) /  5.18000000d-06 ,  1.88000000d-06 /
      data f( 5, 7),f( 5, 8) /  4.10000000d-07 , -4.10000000d-07 /
      data f( 5, 9),f( 5,10) / -3.10000000d-07 , -3.00000000d-08 /
      data f( 5,11),f( 5,12) / -3.00000000d-08 , -1.40000000d-07 /
      data f( 5,13),f( 5,14) / -5.00000000d-08 ,  1.50000000d-07 /
      data f( 5,15),f( 5,16) /  1.15000000d-06 ,  1.80000000d-07 /
      data f( 5,17),f( 5,18) / -1.00000000d-08 ,  5.00000000d-08 /
      data f( 5,19) /           0.00000000d+00 /
      data f( 6, 1),f( 6, 2) /  0.00000000d+00 ,  1.70000000d-07 /
      data f( 6, 3),f( 6, 4) /  6.90000000d-07 ,  1.74000000d-06 /
      data f( 6, 5),f( 6, 6) /  1.37000000d-06 ,  6.00000000d-07 /
      data f( 6, 7),f( 6, 8) /  4.80000000d-07 ,  3.90000000d-07 /
      data f( 6, 9),f( 6,10) /  1.00000000d-08 , -1.00000000d-08 /
      data f( 6,11),f( 6,12) / -3.00000000d-08 , -3.00000000d-07 /
      data f( 6,13),f( 6,14) / -1.50000000d-07 ,  2.00000000d-08 /
      data f( 6,15),f( 6,16) /  5.00000000d-08 , -1.40000000d-07 /
      data f( 6,17),f( 6,18) /  2.00000000d-08 ,  6.50000000d-07 /
      data f( 6,19) /           0.00000000d+00 /
      data f( 7, 1),f( 7, 2) /  0.00000000d+00 ,  1.50000000d-07 /
      data f( 7, 3),f( 7, 4) / -1.00000000d-08 ,  2.80000000d-07 /
      data f( 7, 5),f( 7, 6) /  2.20000000d-07 , -3.60000000d-07 /
      data f( 7, 7),f( 7, 8) / -7.00000000d-07 , -3.60000000d-07 /
      data f( 7, 9),f( 7,10) / -1.07000000d-06 , -1.00000000d-08 /
      data f( 7,11),f( 7,12) /  5.00000000d-08 , -2.90000000d-07 /
      data f( 7,13),f( 7,14) / -2.90000000d-07 , -1.90000000d-07 /
      data f( 7,15),f( 7,16) /  3.00000000d-08 ,  1.30000000d-07 /
      data f( 7,17),f( 7,18) /  1.40000000d-07 ,  3.00000000d-07 /
      data f( 7,19) /           0.00000000d+00 /
      data f( 8, 1),f( 8, 2) /  0.00000000d+00 ,  2.00000000d-08 /
      data f( 8, 3),f( 8, 4) /  0.00000000d+00 , -4.00000000d-08 /
      data f( 8, 5),f( 8, 6) /  4.30000000d-07 ,  4.10000000d-07 /
      data f( 8, 7),f( 8, 8) / -8.00000000d-08 , -1.20000000d-07 /
      data f( 8, 9),f( 8,10) /  4.40000000d-07 , -1.00000000d-08 /
      data f( 8,11),f( 8,12) /  2.90000000d-07 , -1.20000000d-07 /
      data f( 8,13),f( 8,14) / -1.40000000d-07 , -1.60000000d-07 /
      data f( 8,15),f( 8,16) / -2.00000000d-07 , -8.00000000d-08 /
      data f( 8,17),f( 8,18) / -1.10000000d-07 , -3.00000000d-08 /
      data f( 8,19) /           0.00000000d+00 /
      data f( 9, 1),f( 9, 2) /  0.00000000d+00 ,  4.50000000d-07 /
      data f( 9, 3),f( 9, 4) /  1.00000000d-08 , -2.30000000d-07 /
      data f( 9, 5),f( 9, 6) / -3.90000000d-07 , -2.70000000d-07 /
      data f( 9, 7),f( 9, 8) / -8.00000000d-08 ,  5.40000000d-07 /
      data f( 9, 9),f( 9,10) / -1.80000000d-07 ,  4.80000000d-07 /
      data f( 9,11),f( 9,12) / -1.20000000d-07 , -1.30000000d-07 /
      data f( 9,13),f( 9,14) / -1.30000000d-07 ,  3.50000000d-07 /
      data f( 9,15),f( 9,16) /  4.20000000d-07 , -1.90000000d-07 /
      data f( 9,17),f( 9,18) /  4.30000000d-07 , -6.00000000d-08 /
      data f( 9,19) /           0.00000000d+00 /
      data f(10, 1),f(10, 2) /  0.00000000d+00 ,  3.10000000d-07 /
      data f(10, 3),f(10, 4) /  5.50000000d-07 , -2.30000000d-07 /
      data f(10, 5),f(10, 6) / -1.50000000d-07 ,  4.70000000d-07 /
      data f(10, 7),f(10, 8) /  3.20000000d-07 , -1.10000000d-07 /
      data f(10, 9),f(10,10) / -1.40000000d-07 ,  3.70000000d-07 /
      data f(10,11),f(10,12) / -3.00000000d-08 , -1.00000000d-07 /
      data f(10,13),f(10,14) / -1.10000000d-07 , -7.00000000d-08 /
      data f(10,15),f(10,16) / -1.00000000d-07 , -1.50000000d-07 /
      data f(10,17),f(10,18) / -1.30000000d-07 , -1.10000000d-07 /
      data f(10,19) /           0.00000000d+00 /
      data f(11, 1),f(11, 2) /  0.00000000d+00 , -2.20000000d-07 /
      data f(11, 3),f(11, 4) / -1.40000000d-07 , -6.00000000d-08 /
      data f(11, 5),f(11, 6) / -1.00000000d-07 , -6.70000000d-07 /
      data f(11, 7),f(11, 8) / -1.40000000d-07 ,  3.50000000d-07 /
      data f(11, 9),f(11,10) /  2.30000000d-07 , -3.00000000d-08 /
      data f(11,11),f(11,12) /  1.10000000d-07 ,  1.50000000d-07 /
      data f(11,13),f(11,14) / -1.50000000d-07 , -5.00000000d-08 /
      data f(11,15),f(11,16) / -9.00000000d-08 , -9.00000000d-08 /
      data f(11,17),f(11,18) / -8.00000000d-08 , -1.00000000d-07 /
      data f(11,19) /           0.00000000d+00 /
 
      data fpp( 1, 1,1),fpp( 1, 1,2)/ 0.00000000d+00,-8.02228250d-06/
      data fpp( 1, 2,1),fpp( 1, 2,2)/-6.42246466d-03,-3.34403500d-06/
      data fpp( 1, 3,1),fpp( 1, 3,2)/-3.01093141d-03, 1.11962251d-06/
      data fpp( 1, 4,1),fpp( 1, 4,2)/ 2.94799188d-03, 6.33434495d-06/
      data fpp( 1, 5,1),fpp( 1, 5,2)/-1.78135111d-03,-1.87980231d-06/
      data fpp( 1, 6,1),fpp( 1, 6,2)/ 2.95746510d-04,-4.27693570d-06/
      data fpp( 1, 7,1),fpp( 1, 7,2)/ 4.04026357d-05, 2.45814512d-06/
      data fpp( 1, 8,1),fpp( 1, 8,2)/ 8.67762700d-06,-4.76644788d-07/
      data fpp( 1, 9,1),fpp( 1, 9,2)/ 1.17073478d-05, 1.40234030d-07/
      data fpp( 1,10,1),fpp( 1,10,2)/-9.93303904d-07,-2.00913331d-08/
      data fpp( 1,11,1),fpp( 1,11,2)/ 1.16148698d-08,-2.50686980d-08/
      data fpp( 1,12,1),fpp( 1,12,2)/ 2.05050121d-06, 8.07661251d-08/
      data fpp( 1,13,1),fpp( 1,13,2)/ 4.16122551d-05,-3.99580225d-09/
      data fpp( 1,14,1),fpp( 1,14,2)/ 2.16777642d-04, 1.27381708d-06/
      data fpp( 1,15,1),fpp( 1,15,2)/ 2.98823261d-04,-4.52307253d-06/
      data fpp( 1,16,1),fpp( 1,16,2)/ 8.76893519d-05, 1.68076731d-05/
      data fpp( 1,17,1),fpp( 1,17,2)/-8.62099607d-04,-6.68254197d-05/
      data fpp( 1,18,1),fpp( 1,18,2)/-1.96520450d-02, 7.04826056d-05/
      data fpp( 1,19,1),fpp( 1,19,2)/ 0.00000000d+00, 1.44664597d-04/
      data fpp( 2, 1,1),fpp( 2, 1,2)/ 0.00000000d+00,-1.08215837d-05/
      data fpp( 2, 2,1),fpp( 2, 2,2)/-3.05567068d-03,-4.92793265d-06/
      data fpp( 2, 3,1),fpp( 2, 3,2)/-1.60057717d-03,-2.59088574d-06/
      data fpp( 2, 4,1),fpp( 2, 4,2)/ 1.17945624d-03, 1.21942756d-05/
      data fpp( 2, 5,1),fpp( 2, 5,2)/-9.40617780d-04,-9.38521662d-06/
      data fpp( 2, 6,1),fpp( 2, 6,2)/ 1.89306981d-04, 2.75119089d-06/
      data fpp( 2, 7,1),fpp( 2, 7,2)/ 2.48347286d-05,-2.86346930d-07/
      data fpp( 2, 8,1),fpp( 2, 8,2)/ 4.36474599d-06, 1.40796834d-07/
      data fpp( 2, 9,1),fpp( 2, 9,2)/ 5.50530436d-06,-2.54404039d-08/
      data fpp( 2,10,1),fpp( 2,10,2)/-4.93392192d-07, 7.76478193d-09/
      data fpp( 2,11,1),fpp( 2,11,2)/ 1.67702603d-08,-4.41872386d-09/
      data fpp( 2,12,1),fpp( 2,12,2)/ 5.78997585d-07,-8.89886502d-10/
      data fpp( 2,13,1),fpp( 2,13,2)/ 2.22554898d-05,-2.20217301d-08/
      data fpp( 2,14,1),fpp( 2,14,2)/ 1.22204716d-04, 2.05376807d-07/
      data fpp( 2,15,1),fpp( 2,15,2)/ 1.82313477d-04,-2.67885498d-07/
      data fpp( 2,16,1),fpp( 2,16,2)/ 9.42212963d-05, 1.93176519d-06/
      data fpp( 2,17,1),fpp( 2,17,2)/-3.96760787d-04,-6.40977524d-06/
      data fpp( 2,18,1),fpp( 2,18,2)/-1.11053101d-02, 5.33593578d-06/
      data fpp( 2,19,1),fpp( 2,19,2)/ 0.00000000d+00, 1.13424321d-05/
      data fpp( 3, 1,1),fpp( 3, 1,2)/ 0.00000000d+00, 2.41183156d-06/
      data fpp( 3, 2,1),fpp( 3, 2,2)/ 3.98234739d-03, 1.16183687d-06/
      data fpp( 3, 3,1),fpp( 3, 3,2)/ 1.55732011d-03, 6.45420948d-07/
      data fpp( 3, 4,1),fpp( 3, 4,2)/-2.72757683d-03,-2.43852066d-06/
      data fpp( 3, 5,1),fpp( 3, 5,2)/ 3.85022229d-04, 9.05461711d-07/
      data fpp( 3, 6,1),fpp( 3, 6,2)/ 6.15855680d-05, 1.36738199d-08/
      data fpp( 3, 7,1),fpp( 3, 7,2)/ 5.45845005d-06, 1.28843009d-07/
      data fpp( 3, 8,1),fpp( 3, 8,2)/-4.77661097d-06,-1.04585649d-09/
      data fpp( 3, 9,1),fpp( 3, 9,2)/-1.38085653d-05,-7.65958317d-09/
      data fpp( 3,10,1),fpp( 3,10,2)/ 5.66872670d-07, 8.88418917d-09/
      data fpp( 3,11,1),fpp( 3,11,2)/ 1.61304089d-07, 2.12282648d-09/
      data fpp( 3,12,1),fpp( 3,12,2)/-5.08649155d-06,-8.37549509d-09/
      data fpp( 3,13,1),fpp( 3,13,2)/-4.63421431d-06,-3.42084612d-09/
      data fpp( 3,14,1),fpp( 3,14,2)/-1.72765062d-05, 5.25887957d-09/
      data fpp( 3,15,1),fpp( 3,15,2)/ 3.00828301d-05,-3.81467215d-09/
      data fpp( 3,16,1),fpp( 3,16,2)/ 1.88225463d-04, 2.13999809d-07/
      data fpp( 3,17,1),fpp( 3,17,2)/ 4.33622754d-04,-2.92984564d-07/
      data fpp( 3,18,1),fpp( 3,18,2)/ 1.69080518d-03,-2.00615532d-08/
      data fpp( 3,19,1),fpp( 3,19,2)/ 0.00000000d+00, 3.00307766d-08/
      data fpp( 4, 1,1),fpp( 4, 1,2)/ 0.00000000d+00, 1.20135584d-06/
      data fpp( 4, 2,1),fpp( 4, 2,2)/-1.82891887d-03, 4.19188321d-07/
      data fpp( 4, 3,1),fpp( 4, 3,2)/-1.40118326d-03,-1.49990912d-06/
      data fpp( 4, 4,1),fpp( 4, 4,2)/ 8.50131081d-04, 5.60248178d-07/
      data fpp( 4, 5,1),fpp( 4, 5,2)/ 1.34888639d-05, 5.63164116d-08/
      data fpp( 4, 6,1),fpp( 4, 6,2)/ 8.35074750d-06, 4.72861753d-08/
      data fpp( 4, 7,1),fpp( 4, 7,2)/-5.88528745d-07, 2.69388872d-08/
      data fpp( 4, 8,1),fpp( 4, 8,2)/ 2.98169788d-06, 4.55827589d-09/
      data fpp( 4, 9,1),fpp( 4, 9,2)/ 1.27689567d-05, 1.90280092d-08/
      data fpp( 4,10,1),fpp( 4,10,2)/-3.34098490d-07,-1.40703128d-08/
      data fpp( 4,11,1),fpp( 4,11,2)/-1.81986615d-07,-2.34675795d-09/
      data fpp( 4,12,1),fpp( 4,12,2)/ 3.92696860d-06, 1.38573446d-08/
      data fpp( 4,13,1),fpp( 4,13,2)/-5.15863256d-06,-1.58826205d-08/
      data fpp( 4,14,1),fpp( 4,14,2)/ 4.21308783d-07, 1.07313740d-09/
      data fpp( 4,15,1),fpp( 4,15,2)/ 2.87520241d-06, 1.90070919d-10/
      data fpp( 4,16,1),fpp( 4,16,2)/-2.39231486d-05, 4.76657893d-09/
      data fpp( 4,17,1),fpp( 4,17,2)/-1.36502305d-05, 7.74361337d-09/
      data fpp( 4,18,1),fpp( 4,18,2)/-3.85910664d-04,-2.74103239d-09/
      data fpp( 4,19,1),fpp( 4,19,2)/ 0.00000000d+00,-9.37948380d-09/
      data fpp( 5, 1,1),fpp( 5, 1,2)/ 0.00000000d+00, 1.26857684d-07/
      data fpp( 5, 2,1),fpp( 5, 2,2)/ 6.47968084d-04, 5.48846317d-08/
      data fpp( 5, 3,1),fpp( 5, 3,2)/ 7.73092927d-04,-5.35962109d-08/
      data fpp( 5, 4,1),fpp( 5, 4,2)/-7.55874948d-05,-3.42997883d-08/
      data fpp( 5, 5,1),fpp( 5, 5,2)/ 2.54223153d-05,-2.28046361d-08/
      data fpp( 5, 6,1),fpp( 5, 6,2)/ 1.22914420d-05, 3.37183328d-08/
      data fpp( 5, 7,1),fpp( 5, 7,2)/ 8.65566493d-06,-2.26869495d-09/
      data fpp( 5, 8,1),fpp( 5, 8,2)/ 8.20981944d-06, 1.43564470d-08/
      data fpp( 5, 9,1),fpp( 5, 9,2)/-7.87261538d-07, 4.29068882d-11/
      data fpp( 5,10,1),fpp( 5,10,2)/ 2.89521290d-07,-3.72807457d-09/
      data fpp( 5,11,1),fpp( 5,11,2)/ 8.66423731d-08,-1.93060862d-09/
      data fpp( 5,12,1),fpp( 5,12,2)/-2.46138286d-06, 4.85050904d-09/
      data fpp( 5,13,1),fpp( 5,13,2)/ 3.18874455d-06,-5.47142755d-09/
      data fpp( 5,14,1),fpp( 5,14,2)/-1.92872892d-06, 2.36352012d-08/
      data fpp( 5,15,1),fpp( 5,15,2)/-2.07036398d-05,-4.10693771d-08/
      data fpp( 5,16,1),fpp( 5,16,2)/-4.21286886d-06, 2.24423073d-08/
      data fpp( 5,17,1),fpp( 5,17,2)/-2.46183230d-06,-1.89985196d-09/
      data fpp( 5,18,1),fpp( 5,18,2)/ 1.06037476d-04, 1.57100559d-10/
      data fpp( 5,19,1),fpp( 5,19,2)/ 0.00000000d+00,-5.32855028d-09/
      data fpp( 6, 1,1),fpp( 6, 1,2)/ 0.00000000d+00, 2.89865312d-09/
      data fpp( 6, 2,1),fpp( 6, 2,2)/-1.69913468d-04, 1.10269376d-09/
      data fpp( 6, 3,1),fpp( 6, 3,2)/-1.79668451d-04, 1.36905718d-08/
      data fpp( 6, 4,1),fpp( 6, 4,2)/ 4.18988980d-05,-2.40649811d-08/
      data fpp( 6, 5,1),fpp( 6, 5,2)/ 8.42187486d-06,-2.63064748d-09/
      data fpp( 6, 6,1),fpp( 6, 6,2)/-3.03651553d-06, 1.05875710d-08/
      data fpp( 6, 7,1),fpp( 6, 7,2)/-1.19541310d-05,-7.19636598d-10/
      data fpp( 6, 8,1),fpp( 6, 8,2)/-1.27809757d-05,-5.90902463d-09/
      data fpp( 6, 9,1),fpp( 6, 9,2)/-1.10599105d-05, 6.95573511d-09/
      data fpp( 6,10,1),fpp( 6,10,2)/-1.03986672d-07,-3.13915814d-10/
      data fpp( 6,11,1),fpp( 6,11,2)/ 3.15417123d-07,-5.70007186d-09/
      data fpp( 6,12,1),fpp( 6,12,2)/ 1.59856283d-06, 8.11420324d-09/
      data fpp( 6,13,1),fpp( 6,13,2)/-1.35634563d-06,-1.55674109d-09/
      data fpp( 6,14,1),fpp( 6,14,2)/-1.46393098d-07,-6.87238879d-10/
      data fpp( 6,15,1),fpp( 6,15,2)/ 1.24993567d-05,-4.09430339d-09/
      data fpp( 6,16,1),fpp( 6,16,2)/ 5.25462398d-06, 3.86445246d-09/
      data fpp( 6,17,1),fpp( 6,17,2)/ 1.65755971d-06, 9.63649357d-09/
      data fpp( 6,18,1),fpp( 6,18,2)/-3.39192409d-05,-1.42104267d-08/
      data fpp( 6,19,1),fpp( 6,19,2)/ 0.00000000d+00,-2.95947866d-08/
      data fpp( 7, 1,1),fpp( 7, 1,2)/ 0.00000000d+00,-1.00053730d-08/
      data fpp( 7, 2,1),fpp( 7, 2,2)/ 3.04857893d-05,-4.48925407d-09/
      data fpp( 7, 3,1),fpp( 7, 3,2)/ 3.60608763d-05, 9.36238924d-09/
      data fpp( 7, 4,1),fpp( 7, 4,2)/-2.00809732d-06,-5.96030289d-09/
      data fpp( 7, 5,1),fpp( 7, 5,2)/ 4.73018527d-06,-6.52117768d-09/
      data fpp( 7, 6,1),fpp( 7, 6,2)/ 7.53462011d-06, 8.45013611d-10/
      data fpp( 7, 7,1),fpp( 7, 7,2)/ 9.16085902d-06, 1.75411232d-08/
      data fpp( 7, 8,1),fpp( 7, 8,2)/ 5.71408317d-06,-3.02095066d-08/
      data fpp( 7, 9,1),fpp( 7, 9,2)/ 1.14269037d-05, 4.02969030d-08/
      data fpp( 7,10,1),fpp( 7,10,2)/-3.53574604d-07,-2.47781054d-08/
      data fpp( 7,11,1),fpp( 7,11,2)/ 5.71689135d-07,-1.18448140d-09/
      data fpp( 7,12,1),fpp( 7,12,2)/ 1.47131544d-07, 5.51603099d-09/
      data fpp( 7,13,1),fpp( 7,13,2)/ 1.27663799d-06,-4.79642577d-10/
      data fpp( 7,14,1),fpp( 7,14,2)/ 5.94301314d-07, 2.40253931d-09/
      data fpp( 7,15,1),fpp( 7,15,2)/-3.37378692d-06,-1.93051468d-09/
      data fpp( 7,16,1),fpp( 7,16,2)/-2.64562707d-06,-1.88048060d-09/
      data fpp( 7,17,1),fpp( 7,17,2)/-2.00840655d-06, 4.05243709d-09/
      data fpp( 7,18,1),fpp( 7,18,2)/ 6.83948738d-06,-5.32926774d-09/
      data fpp( 7,19,1),fpp( 7,19,2)/ 0.00000000d+00,-1.03353661d-08/
      data fpp( 8, 1,1),fpp( 8, 1,2)/ 0.00000000d+00,-1.00456555d-09/
      data fpp( 8, 2,1),fpp( 8, 2,2)/-7.04063376d-06, 4.09131102d-10/
      data fpp( 8, 3,1),fpp( 8, 3,2)/-9.88840353d-06,-3.03195886d-09/
      data fpp( 8, 4,1),fpp( 8, 4,2)/ 6.74842944d-07, 1.05187043d-08/
      data fpp( 8, 5,1),fpp( 8, 5,2)/-3.34149324d-06,-8.44285845d-09/
      data fpp( 8, 6,1),fpp( 8, 6,2)/-4.94560256d-06,-6.14727053d-09/
      data fpp( 8, 7,1),fpp( 8, 7,2)/-3.62551157d-06, 4.83194057d-09/
      data fpp( 8, 8,1),fpp( 8, 8,2)/-3.11761694d-07, 1.38195082d-08/
      data fpp( 8, 9,1),fpp( 8, 9,2)/-6.73075593d-06,-2.41099736d-08/
      data fpp( 8,10,1),fpp( 8,10,2)/ 1.11271715d-06, 2.20203860d-08/
      data fpp( 8,11,1),fpp( 8,11,2)/-1.39277597d-06,-1.89715704d-08/
      data fpp( 8,12,1),fpp( 8,12,2)/-3.40676047d-07, 1.12658956d-08/
      data fpp( 8,13,1),fpp( 8,13,2)/-5.71741152d-07,-2.69201210d-09/
      data fpp( 8,14,1),fpp( 8,14,2)/ 9.90292606d-07,-4.97847240d-10/
      data fpp( 8,15,1),fpp( 8,15,2)/ 2.73168241d-06, 3.48340106d-09/
      data fpp( 8,16,1),fpp( 8,16,2)/ 8.09569216d-07,-3.83575698d-09/
      data fpp( 8,17,1),fpp( 8,17,2)/ 2.25643978d-06, 2.85962688d-09/
      data fpp( 8,18,1),fpp( 8,18,2)/-1.33884169d-06,-1.00275054d-09/
      data fpp( 8,19,1),fpp( 8,19,2)/ 0.00000000d+00,-1.84862473d-09/
      data fpp( 9, 1,1),fpp( 9, 1,2)/ 0.00000000d+00,-1.92379826d-08/
      data fpp( 9, 2,1),fpp( 9, 2,2)/ 1.03674573d-06,-1.00240348d-08/
      data fpp( 9, 3,1),fpp( 9, 3,2)/ 3.49273777d-06, 5.93412188d-09/
      data fpp( 9, 4,1),fpp( 9, 4,2)/ 8.87255409d-08,-1.71245271d-09/
      data fpp( 9, 5,1),fpp( 9, 5,2)/ 2.45578770d-06, 5.71568897d-09/
      data fpp( 9, 6,1),fpp( 9, 6,2)/ 3.54779015d-06,-4.35030317d-09/
      data fpp( 9, 7,1),fpp( 9, 7,2)/ 1.62118726d-06, 1.58855237d-08/
      data fpp( 9, 8,1),fpp( 9, 8,2)/-1.94703640d-06,-3.33917917d-08/
      data fpp( 9, 9,1),fpp( 9, 9,2)/ 2.71611998d-06, 3.72816432d-08/
      data fpp( 9,10,1),fpp( 9,10,2)/-1.15729399d-06,-3.29347811d-08/
      data fpp( 9,11,1),fpp( 9,11,2)/ 1.09941473d-06, 1.88574812d-08/
      data fpp( 9,12,1),fpp( 9,12,2)/ 1.35572644d-07,-7.09514377d-09/
      data fpp( 9,13,1),fpp( 9,13,2)/ 1.70326619d-07, 1.01230939d-08/
      data fpp( 9,14,1),fpp( 9,14,2)/-1.67547174d-06,-4.59723166d-09/
      data fpp( 9,15,1),fpp( 9,15,2)/-2.45294274d-06,-1.63341672d-08/
      data fpp( 9,16,1),fpp( 9,16,2)/ 7.35020643d-09, 2.91339005d-08/
      data fpp( 9,17,1),fpp( 9,17,2)/-2.27735257d-06,-2.64014347d-08/
      data fpp( 9,18,1),fpp( 9,18,2)/ 3.15879393d-07, 9.87183850d-09/
      data fpp( 9,19,1),fpp( 9,19,2)/ 0.00000000d+00, 1.99140808d-08/
      data fpp(10, 1,1),fpp(10, 1,2)/ 0.00000000d+00, 7.42603210d-09/
      data fpp(10, 2,1),fpp(10, 2,2)/-5.26349145d-07, 2.04793579d-09/
      data fpp(10, 3,1),fpp(10, 3,2)/-9.02547554d-07,-1.98177753d-08/
      data fpp(10, 4,1),fpp(10, 4,2)/ 1.10254892d-07, 1.60231653d-08/
      data fpp(10, 5,1),fpp(10, 5,2)/-1.21657539d-07, 7.32511399d-09/
      data fpp(10, 6,1),fpp(10, 6,2)/-7.25558030d-07,-1.29236213d-08/
      data fpp(10, 7,1),fpp(10, 7,2)/-4.59237451d-07,-1.83062887d-09/
      data fpp(10, 8,1),fpp(10, 8,2)/ 2.39907279d-07, 3.44613676d-09/
      data fpp(10, 9,1),fpp(10, 9,2)/-1.73723996d-07, 1.20460818d-08/
      data fpp(10,10,1),fpp(10,10,2)/-8.35412028d-08,-1.92304641d-08/
      data fpp(10,11,1),fpp(10,11,2)/-4.88294556d-09, 1.02757745d-08/
      data fpp(10,12,1),fpp(10,12,2)/ 3.83854712d-08,-2.07263407d-09/
      data fpp(10,13,1),fpp(10,13,2)/-4.95653238d-08, 1.61476173d-09/
      data fpp(10,14,1),fpp(10,14,2)/ 1.31594348d-07,-1.38641284d-09/
      data fpp(10,15,1),fpp(10,15,2)/ 2.40088548d-07,-2.69110372d-10/
      data fpp(10,16,1),fpp(10,16,2)/ 6.10299587d-08, 1.26285433d-09/
      data fpp(10,17,1),fpp(10,17,2)/ 2.52970515d-07,-5.82306934d-10/
      data fpp(10,18,1),fpp(10,18,2)/-4.46758786d-08, 1.06637341d-09/
      data fpp(10,19,1),fpp(10,19,2)/ 0.00000000d+00, 1.71681330d-09/
      data fpp(11, 1,1),fpp(11, 1,2)/ 0.00000000d+00, 5.80969958d-09/
      data fpp(11, 2,1),fpp(11, 2,2)/ 6.85674573d-07, 3.38060085d-09/
      data fpp(11, 3,1),fpp(11, 3,2)/-1.69372622d-06,-1.33210297d-09/
      data fpp(11, 4,1),fpp(11, 4,2)/-1.20127446d-07, 1.94781104d-09/
      data fpp(11, 5,1),fpp(11, 5,2)/-1.50792123d-06,-1.36591412d-08/
      data fpp(11, 6,1),fpp(11, 6,2)/-3.52722099d-06, 2.08887537d-08/
      data fpp(11, 7,1),fpp(11, 7,2)/-1.32288127d-06,-3.89587353d-09/
      data fpp(11, 8,1),fpp(11, 8,2)/ 2.89379636d-06,-7.70525956d-09/
      data fpp(11, 9,1),fpp(11, 9,2)/-4.01888002d-07,-1.88308824d-09/
      data fpp(11,10,1),fpp(11,10,2)/ 5.59270601d-07, 6.83761253d-09/
      data fpp(11,11,1),fpp(11,11,2)/-5.95058527d-07,-1.46736187d-09/
      data fpp(11,12,1),fpp(11,12,2)/ 1.02057264d-07,-6.96816503d-09/
      data fpp(11,13,1),fpp(11,13,2)/-5.64673381d-08, 8.94002201d-09/
      data fpp(11,14,1),fpp(11,14,2)/ 1.73295283d-06,-4.79192300d-09/
      data fpp(11,15,1),fpp(11,15,2)/ 2.08120573d-06, 1.82766997d-09/
      data fpp(11,16,1),fpp(11,16,2)/-2.16764979d-07,-1.18756900d-10/
      data fpp(11,17,1),fpp(11,17,2)/ 2.13476474d-06,-7.52642373d-10/
      data fpp(11,18,1),fpp(11,18,2)/ 1.41087939d-07, 1.32932639d-09/
      data fpp(11,19,1),fpp(11,19,2)/ 0.00000000d+00, 2.63533680d-09/
 
      data fpppp( 1, 1),fpppp( 1, 2)/ 1.77977107d-04, 8.47978794d-05/
      data fpppp( 1, 3),fpppp( 1, 4)/ 7.28712493d-05,-2.23439474d-04/
      data fpppp( 1, 5),fpppp( 1, 6)/ 1.79590669d-04,-8.65367643d-05/
      data fpppp( 1, 7),fpppp( 1, 8)/ 2.66098988d-05,-6.48569903d-06/
      data fpppp( 1, 9),fpppp( 1,10)/ 1.41818110d-06,-1.30847700d-07/
      data fpppp( 1,11),fpppp( 1,12)/-7.24560661d-08, 4.82710018d-07/
      data fpppp( 1,13),fpppp( 1,14)/ 3.92988047d-07, 6.08155577d-06/
      data fpppp( 1,15),fpppp( 1,16)/-3.03063972d-05, 9.75532612d-05/
      data fpppp( 1,17),fpppp( 1,18)/-4.04225951d-04, 4.48941157d-04/
      data fpppp( 1,19) /             9.14980744d-04 /
      data fpppp( 2, 1),fpppp( 2, 2)/ 8.02512412d-05, 3.85369290d-05/
      data fpppp( 2, 3),fpppp( 2, 4)/ 3.62468940d-05,-1.04028111d-04/
      data fpppp( 2, 5),fpppp( 2, 6)/ 8.58591044d-05,-4.44083800d-05/
      data fpppp( 2, 7),fpppp( 2, 8)/ 1.41105947d-05,-3.39386275d-06/
      data fpppp( 2, 9),fpppp( 2,10)/ 7.61488734d-07,-8.04474820d-08/
      data fpppp( 2,11),fpppp( 2,12)/-4.91672656d-08, 2.80240437d-07/
      data fpppp( 2,13),fpppp( 2,14)/ 1.95061411d-07, 3.63587796d-06/
      data fpppp( 2,15),fpppp( 2,16)/-1.71290012d-05, 5.59880702d-05/
      data fpppp( 2,17),fpppp( 2,18)/-2.30996674d-04, 2.54944594d-04/
      data fpppp( 2,19) /             5.20049858d-04 /
      data fpppp( 3, 1),fpppp( 3, 2)/-1.14196878d-04,-5.47775837d-05/
      data fpppp( 3, 3),fpppp( 3, 4)/-5.11352671d-05, 1.47726472d-04/
      data fpppp( 3, 5),fpppp( 3, 6)/-9.59208630d-05, 2.97948364d-05/
      data fpppp( 3, 7),fpppp( 3, 8)/-7.21990987d-06, 1.83832655d-06/
      data fpppp( 3, 9),fpppp( 3,10)/-6.12099137d-08,-1.89043359d-07/
      data fpppp( 3,11),fpppp( 3,12)/-6.94770425d-08, 1.76417905d-07/
      data fpppp( 3,13),fpppp( 3,14)/-2.94190207d-07, 2.14668775d-07/
      data fpppp( 3,15),fpppp( 3,16)/ 3.03561280d-06,-5.71012218d-06/
      data fpppp( 3,17),fpppp( 3,18)/ 2.50401554d-05,-3.37433915d-05/
      data fpppp( 3,19) /            -6.69458457d-05 /
      data fpppp( 4, 1),fpppp( 4, 2)/ 2.98828886d-05, 1.65953720d-05/
      data fpppp( 4, 3),fpppp( 4, 4)/ 3.91348921d-05,-6.37202165d-05/
      data fpppp( 4, 5),fpppp( 4, 6)/ 3.04685805d-05,-8.26385934d-06/
      data fpppp( 4, 7),fpppp( 4, 8)/ 2.35878730d-06,-4.20719704d-07/
      data fpppp( 4, 9),fpppp( 4,10)/-3.02886558d-07, 2.58847095d-07/
      data fpppp( 4,11),fpppp( 4,12)/ 6.28082013d-08,-2.72669300d-07/
      data fpppp( 4,13),fpppp( 4,14)/ 2.36195616d-07, 2.07819387d-07/
      data fpppp( 4,15),fpppp( 4,16)/-1.25503602d-06, 3.05719004d-06/
      data fpppp( 4,17),fpppp( 4,18)/-8.74944798d-06, 8.98860082d-06/
      data fpppp( 4,19) /             1.82853105d-05 /
      data fpppp( 5, 1),fpppp( 5, 2)/-2.07121458d-06,-2.52362756d-06/
      data fpppp( 5, 3),fpppp( 5, 4)/-1.92048696d-05, 2.09147901d-05/
      data fpppp( 5, 5),fpppp( 5, 6)/-7.47287684d-06, 2.12827625d-06/
      data fpppp( 5, 7),fpppp( 5, 8)/-4.70522388d-07,-5.47908034d-08/
      data fpppp( 5, 9),fpppp( 5,10)/ 1.76611472d-07,-4.72232560d-08/
      data fpppp( 5,11),fpppp( 5,12)/-6.44981526d-08, 1.64507088d-07/
      data fpppp( 5,13),fpppp( 5,14)/-1.01641040d-07,-4.03998979d-07/
      data fpppp( 5,15),fpppp( 5,16)/ 8.98190713d-07,-1.07282297d-06/
      data fpppp( 5,17),fpppp( 5,18)/ 2.50871709d-06,-2.55714910d-06/
      data fpppp( 5,19) /            -5.15232779d-06 /
      data fpppp( 6, 1),fpppp( 6, 2)/ 1.21490867d-06, 9.51660311d-07/
      data fpppp( 6, 3),fpppp( 6, 4)/ 4.58795923d-06,-5.42415733d-06/
      data fpppp( 6, 5),fpppp( 6, 6)/ 1.80600775d-06,-4.78755719d-07/
      data fpppp( 6, 7),fpppp( 6, 8)/ 2.61461619d-07,-8.16445091d-08/
      data fpppp( 6, 9),fpppp( 6,10)/ 2.17991004d-07,-2.36227979d-07/
      data fpppp( 6,11),fpppp( 6,12)/ 9.47297076d-08,-9.08663369d-08/
      data fpppp( 6,13),fpppp( 6,14)/ 1.44523899d-08, 2.82948437d-07/
      data fpppp( 6,15),fpppp( 6,16)/-4.60098305d-07, 3.64015836d-07/
      data fpppp( 6,17),fpppp( 6,18)/-7.77104934d-07, 8.25619721d-07/
      data fpppp( 6,19) /             1.64438854d-06 /
      data fpppp( 7, 1),fpppp( 7, 2)/-1.23099989d-07,-1.26453922d-07/
      data fpppp( 7, 3),fpppp( 7, 4)/-8.65726458d-07, 9.70716115d-07/
      data fpppp( 7, 5),fpppp( 7, 6)/-3.28702628d-07, 1.08063531d-07/
      data fpppp( 7, 7),fpppp( 7, 8)/-1.74243252d-07, 2.84528593d-07/
      data fpppp( 7, 9),fpppp( 7,10)/-4.14295335d-07, 3.23054812d-07/
      data fpppp( 7,11),fpppp( 7,12)/-1.15579389d-07, 5.82734653d-08/
      data fpppp( 7,13),fpppp( 7,14)/-2.42706299d-08,-6.99015329d-08/
      data fpppp( 7,15),fpppp( 7,16)/ 1.06731668d-07,-7.52502554d-08/
      data fpppp( 7,17),fpppp( 7,18)/ 1.88812994d-07,-1.87361316d-07/
      data fpppp( 7,19) /            -3.80610609d-07 /
      data fpppp( 8, 1),fpppp( 8, 2)/-2.93012109d-08, 2.52967123d-11/
      data fpppp( 8, 3),fpppp( 8, 4)/ 2.80771864d-07,-3.18451777d-07/
      data fpppp( 8, 5),fpppp( 8, 6)/ 1.18260285d-07,-9.85574980d-09/
      data fpppp( 8, 7),fpppp( 8, 8)/ 9.66147337d-08,-2.56983652d-07/
      data fpppp( 8, 9),fpppp( 8,10)/ 3.47355229d-07,-2.76689224d-07/
      data fpppp( 8,11),fpppp( 8,12)/ 1.38463695d-07,-6.37099729d-08/
      data fpppp( 8,13),fpppp( 8,14)/ 3.93862956d-08, 1.37507224d-08/
      data fpppp( 8,15),fpppp( 8,16)/-8.36278223d-08, 1.00950386d-07/
      data fpppp( 8,17),fpppp( 8,18)/-1.18034698d-07, 6.86592823d-08/
      data fpppp( 8,19) /             1.39444958d-07 /
      data fpppp( 9, 1),fpppp( 9, 2)/ 7.68994862d-08, 3.43634290d-08/
      data fpppp( 9, 3),fpppp( 9, 4)/-1.29198423d-07, 1.30830007d-07/
      data fpppp( 9, 5),fpppp( 9, 6)/-4.78571415d-08,-1.59050230d-08/
      data fpppp( 9, 7),fpppp( 9, 8)/-6.96390874d-08, 1.95964127d-07/
      data fpppp( 9, 9),fpppp( 9,10)/-2.20334620d-07, 1.73180132d-07/
      data fpppp( 9,11),fpppp( 9,12)/-1.04578545d-07, 5.19010021d-08/
      data fpppp( 9,13),fpppp( 9,14)/-4.31096995d-08, 7.70465593d-09/
      data fpppp( 9,15),fpppp( 9,16)/ 7.63907171d-08,-1.19001687d-07/
      data fpppp( 9,17),fpppp( 9,18)/ 1.14916289d-07,-4.79873825d-08/
      data fpppp( 9,19) /            -9.75134402d-08 /
      data fpppp(10, 1),fpppp(10, 2)/-8.91214344d-09,-2.44819340d-09/
      data fpppp(10, 3),fpppp(10, 4)/ 2.77139612d-08,-2.50676002d-08/
      data fpppp(10, 5),fpppp(10, 6)/-2.12645308d-09, 1.12541289d-08/
      data fpppp(10, 7),fpppp(10, 8)/ 9.32320153d-09,-2.25774860d-08/
      data fpppp(10, 9),fpppp(10,10)/ 1.42201820d-08,-4.07439781d-09/
      data fpppp(10,11),fpppp(10,12)/ 1.38593708d-09,-3.59274092d-09/
      data fpppp(10,13),fpppp(10,14)/ 5.11187390d-09,-7.08126676d-10/
      data fpppp(10,15),fpppp(10,16)/-6.63929546d-09, 1.00121411d-08/
      data fpppp(10,17),fpppp(10,18)/-1.11493202d-08, 5.20992293d-09/
      data fpppp(10,19) /             1.08489648d-08 /
      data fpppp(11, 1),fpppp(11, 2)/-9.41298556d-08,-4.40540486d-08/
      data fpppp(11, 3),fpppp(11, 4)/ 8.64415278d-08,-6.45320884d-08/
      data fpppp(11, 5),fpppp(11, 6)/-5.99672780d-09, 5.06286414d-08/
      data fpppp(11, 7),fpppp(11, 8)/ 5.69005299d-08,-1.57490486d-07/
      data fpppp(11, 9),fpppp(11,10)/ 1.22319693d-07,-7.63777085d-08/
      data fpppp(11,11),fpppp(11,12)/ 5.62618769d-08,-3.75831039d-08/
      data fpppp(11,13),fpppp(11,14)/ 4.27321149d-08,-1.64686696d-08/
      data fpppp(11,15),fpppp(11,16)/-6.33274722d-08, 1.11005142d-07/
      data fpppp(11,17),fpppp(11,18)/-1.01723070d-07, 3.51747472d-08/
      data fpppp(11,19) /             7.21794133d-08 /

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
      subroutine ebasis(rca,rcc,etz)
      implicit real*8 (a-h,o-z)
      call ebch2cch(rca,ea)
      call ebhccch2(rcc,ec)
      etz=min(ea,ec)
      return
      end
      subroutine ebch2cch(xi,fi)
c
c     1d spline for h addition to the ch2 side of h2ccch 
c     difference between pvtz and pvdz calculations
c
      implicit real*8 (a-h,o-z)
      dimension x( 17), f( 17),fpp( 17),del( 16)
 
      data fpp(  1), fpp(  2) /  1.1527158116d+00 ,  5.8256837685d-01 /
      data fpp(  3), fpp(  4) /  1.1701068101d-01 , -7.1461110091d-01 /
      data fpp(  5), fpp(  6) / -4.2656627737d-01 , -1.3899561743d-01 /
      data fpp(  7), fpp(  8) /  2.1269990990d-02 , -5.6287110690d-03 /
      data fpp(  9), fpp( 10) /  1.2448532862d-03 , -9.2020432424d-04 /
      data fpp( 11), fpp( 12) / -2.9684806271d-03 ,  7.9412683248d-04 /
      data fpp( 13), fpp( 14) / -2.0802670286d-04 ,  3.7979978945d-05 /
      data fpp( 15), fpp( 16) / -9.9265854060d-06 ,  1.7263626793d-06 /
      data fpp( 17) / -8.6318133965d-07 /
 
      data f(  1), f(  2) / -1.5780000000d+00 , -1.3460000000d+00 /
      data f(  3), f(  4) / -9.6400000000d-01 , -5.6800000000d-01 /
      data f(  5), f(  6) / -3.0400000000d-01 , -7.2000000000d-02 /
      data f(  7), f(  8) / -1.4000000000d-02 , -1.1000000000d-02 /
      data f(  9), f( 10) / -8.0000000000d-03 , -2.0000000000d-03 /
      data f( 11), f( 12) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 13), f( 14) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 15), f( 16) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 17) /  0.0000000000d+00 /
 
      data x(  1), x(  2) /  3.0000000000d+00 ,  3.5000000000d+00 /
      data x(  3), x(  4) /  4.0000000000d+00 ,  4.5000000000d+00 /
      data x(  5), x(  6) /  5.0000000000d+00 ,  6.0000000000d+00 /
      data x(  7), x(  8) /  8.0000000000d+00 ,  9.0000000000d+00 /
      data x(  9), x( 10) /  1.0000000000d+01 ,  1.2000000000d+01 /
      data x( 11), x( 12) /  1.3000000000d+01 ,  1.4000000000d+01 /
      data x( 13), x( 14) /  1.5000000000d+01 ,  1.6000000000d+01 /
      data x( 15), x( 16) /  1.8000000000d+01 ,  2.0000000000d+01 /
      data x( 17) /  2.5000000000d+01 /
 
      data del(  1), del(  2) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  3), del(  4) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  5), del(  6) /  1.0000000000d+00 ,  2.0000000000d+00 /
      data del(  7), del(  8) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del(  9), del( 10) /  2.0000000000d+00 ,  1.0000000000d+00 /
      data del( 11), del( 12) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 13), del( 14) /  1.0000000000d+00 ,  2.0000000000d+00 /
      data del( 15), del( 16) /  2.0000000000d+00 ,  5.0000000000d+00 /
      data npts / 17/
 

      if(xi .le. x(1)) then
        ii=1
      else if( xi .ge. x(npts)) then
        ii=npts-1
      else
        call hunt(x,npts,xi,ii)
      endif
      
 20   fi = fpp(ii)   * (x(ii+1)-xi)**3 / (6.0d0*del(ii)) + 
     x     fpp(ii+1) * (xi-x(ii))**3   / (6.0d0*del(ii)) +
     x   ((f(ii+1)/del(ii))-(fpp(ii+1)*del(ii)/6.0d0)) * (xi-x(ii)) + 
     x   ((f(ii)  /del(ii))-(fpp(ii)  *del(ii)/6.0d0)) * (x(ii+1)-xi)
      fi = fi/627.51d0
      return
      end  
      subroutine ebhccch2(xi,fi)
c
c     1d spline for h addition to the ch side of h2ccch 
c     difference between pvtz and pvdz calculations
c
      implicit real*8 (a-h,o-z)
      dimension x( 18), f( 18),fpp( 18),del( 17)
 
      data fpp(  1), fpp(  2) /  8.3134666005d-01 ,  2.2530667989d-01 /
      data fpp(  3), fpp(  4) / -4.3657337962d-01 , -9.0301316141d-01 /
      data fpp(  5), fpp(  6) / -4.6337397473d-01 , -2.1949093968d-01 /
      data fpp(  7), fpp(  8) / -9.8662266552d-02 ,  1.3528401300d-02 /
      data fpp(  9), fpp( 10) /  2.1541253024d-03 , -4.1449025099d-03 /
      data fpp( 11), fpp( 12) /  3.8576448786d-03 , -5.8560642519d-03 /
      data fpp( 13), fpp( 14) /  1.5666121290d-03 , -4.1038426424d-04 /
      data fpp( 15), fpp( 16) /  7.4924927911d-05 , -1.9582651613d-05 /
      data fpp( 17), fpp( 18) /  3.4056785414d-06 , -1.7028392707d-06 /
 
      data f(  1), f(  2) / -1.9070000000d+00 , -1.4300000000d+00 /
      data f(  3), f(  4) / -8.9900000000d-01 , -4.6900000000d-01 /
      data f(  5), f(  6) / -2.2700000000d-01 , -1.0900000000d-01 /
      data f(  7), f(  8) / -5.1000000000d-02 , -1.1000000000d-02 /
      data f(  9), f( 10) / -1.0000000000d-02 , -6.0000000000d-03 /
      data f( 11), f( 12) / -3.0000000000d-03 ,  0.0000000000d+00 /
      data f( 13), f( 14) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 15), f( 16) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 17), f( 18) /  0.0000000000d+00 ,  0.0000000000d+00 /
 
      data x(  1), x(  2) /  3.0000000000d+00 ,  3.5000000000d+00 /
      data x(  3), x(  4) /  4.0000000000d+00 ,  4.5000000000d+00 /
      data x(  5), x(  6) /  5.0000000000d+00 ,  5.5000000000d+00 /
      data x(  7), x(  8) /  6.0000000000d+00 ,  8.0000000000d+00 /
      data x(  9), x( 10) /  9.0000000000d+00 ,  1.0000000000d+01 /
      data x( 11), x( 12) /  1.2000000000d+01 ,  1.3000000000d+01 /
      data x( 13), x( 14) /  1.4000000000d+01 ,  1.5000000000d+01 /
      data x( 15), x( 16) /  1.6000000000d+01 ,  1.8000000000d+01 /
      data x( 17), x( 18) /  2.0000000000d+01 ,  2.5000000000d+01 /
 
      data del(  1), del(  2) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  3), del(  4) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  5), del(  6) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  7), del(  8) /  2.0000000000d+00 ,  1.0000000000d+00 /
      data del(  9), del( 10) /  1.0000000000d+00 ,  2.0000000000d+00 /
      data del( 11), del( 12) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 13), del( 14) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 15), del( 16) /  2.0000000000d+00 ,  2.0000000000d+00 /
      data del( 17) /  5.0000000000d+00 /
      data npts / 18/
 

      if(xi .le. x(1)) then
        ii=1
      else if( xi .ge. x(npts)) then
        ii=npts-1
      else
        call hunt(x,npts,xi,ii)
      endif
      
 20   fi = fpp(ii)   * (x(ii+1)-xi)**3 / (6.0d0*del(ii)) + 
     x     fpp(ii+1) * (xi-x(ii))**3   / (6.0d0*del(ii)) +
     x   ((f(ii+1)/del(ii))-(fpp(ii+1)*del(ii)/6.0d0)) * (xi-x(ii)) + 
     x   ((f(ii)  /del(ii))-(fpp(ii)  *del(ii)/6.0d0)) * (x(ii+1)-xi)
      fi = fi/627.51d0
      return
      end  
      subroutine hhspl(xi,fi)
c
c     h-h nonbonded potential from:
c     ch-h 3b1 cas+1+2/cc-pvdz
c     xi = hh bond distance in au
c     binding energy returned in au
c
      implicit real*8 (a-h,o-z)
      dimension x( 72), f( 72),fpp( 72),del( 71)
 
      data fpp(  1), fpp(  2) /  1.0589750482d+03 ,  6.8936827155d+02 /
      data fpp(  3), fpp(  4) /  1.9635704958d+02 ,  1.3526138614d+02 /
      data fpp(  5), fpp(  6) /  6.9996117850d+01 ,  4.7066166459d+01 /
      data fpp(  7), fpp(  8) /  3.0917808316d+01 ,  2.1954256276d+01 /
      data fpp(  9), fpp( 10) /  1.5894438581d+01 ,  1.2027173402d+01 /
      data fpp( 11), fpp( 12) /  9.0123558111d+00 ,  7.1215553537d+00 /
      data fpp( 13), fpp( 14) /  5.5933987743d+00 ,  4.4776415492d+00 /
      data fpp( 15), fpp( 16) /  3.6130590289d+00 ,  2.9437463352d+00 /
      data fpp( 17), fpp( 18) /  2.4160676301d+00 ,  1.9964871442d+00 /
      data fpp( 19), fpp( 20) /  1.6584637930d+00 ,  1.3868736839d+00 /
      data fpp( 21), fpp( 22) /  1.1561454714d+00 ,  8.1876274379d-01 /
      data fpp( 23), fpp( 24) /  5.9047355343d-01 ,  4.2831904249d-01 /
      data fpp( 25), fpp( 26) /  3.1265227660d-01 ,  2.2863385110d-01 /
      data fpp( 27), fpp( 28) /  1.6703831901d-01 ,  1.2159887286d-01 /
      data fpp( 29), fpp( 30) /  8.8112189536d-02 ,  6.3582368993d-02 /
      data fpp( 31), fpp( 32) /  4.5824334493d-02 ,  3.3152293036d-02 /
      data fpp( 33), fpp( 34) /  2.4264493365d-02 ,  1.8133733505d-02 /
      data fpp( 35), fpp( 36) /  1.3966572616d-02 ,  1.1163976033d-02 /
      data fpp( 37), fpp( 38) /  9.2675232533d-03 ,  7.9799309539d-03 /
      data fpp( 39), fpp( 40) /  7.0527529310d-03 ,  6.3730573219d-03 /
      data fpp( 41), fpp( 42) /  5.8190177812d-03 ,  5.3548715532d-03 /
      data fpp( 43), fpp( 44) /  4.9454960060d-03 ,  4.5571444228d-03 /
      data fpp( 45), fpp( 46) /  4.1979263026d-03 ,  3.8511503667d-03 /
      data fpp( 47), fpp( 48) /  3.5214722306d-03 ,  3.2009607109d-03 /
      data fpp( 49), fpp( 50) /  2.9046849257d-03 ,  2.6162995864d-03 /
      data fpp( 51), fpp( 52) /  2.1092587780d-03 ,  1.6756653014d-03 /
      data fpp( 53), fpp( 54) /  1.3115800162d-03 ,  1.0070146337d-03 /
      data fpp( 55), fpp( 56) /  7.5336144883d-04 ,  5.5303957095d-04 /
      data fpp( 57), fpp( 58) /  3.9598026736d-04 ,  2.7953935962d-04 /
      data fpp( 59), fpp( 60) /  1.9286229416d-04 ,  1.2501146376d-04 /
      data fpp( 61), fpp( 62) /  3.9022983817d-05 ,  1.0256600974d-05 /
      data fpp( 63), fpp( 64) /  1.0706122885d-06 , -3.7390501275d-06 /
      data fpp( 65), fpp( 66) / -4.3544117786d-06 , -5.4833027583d-06 /
      data fpp( 67), fpp( 68) / -2.6728858360d-06 , -8.0515389792d-07 /
      data fpp( 69), fpp( 70) / -2.2649857235d-07 , -9.5250775167d-10 /
      data fpp( 71), fpp( 72) / -1.8652831783d-08 ,  4.0834749225d-08 /
 
      data f(  1), f(  2) /  8.6995441100d+00 ,  5.3782711800d+00 /
      data f(  3), f(  4) /  3.7290004100d+00 ,  2.7505870800d+00 /
      data f(  5), f(  6) /  2.1085898800d+00 ,  1.6592226900d+00 /
      data f(  7), f(  8) /  1.3303465800d+00 ,  1.0817586600d+00 /
      data f(  9), f( 10) /  8.8926627000d-01 ,  7.3742354000d-01 /
      data f( 11), f( 12) /  6.1600393000d-01 ,  5.1758355000d-01 /
      data f( 13), f( 14) /  4.3711816000d-01 ,  3.7080810000d-01 /
      data f( 15), f( 16) /  3.1579680000d-01 ,  2.6989951000d-01 /
      data f( 17), f( 18) /  2.3142060000d-01 ,  1.9902690000d-01 /
      data f( 19), f( 20) /  1.7165840000d-01 ,  1.4846374000d-01 /
      data f( 21), f( 22) /  1.2875329000d-01 ,  9.7633450000d-02 /
      data f( 23), f( 24) /  7.4883060000d-02 ,  5.8147630000d-02 /
      data f( 25), f( 26) /  4.5772870000d-02 ,  3.6577380000d-02 /
      data f( 27), f( 28) /  2.9705600000d-02 ,  2.4531130000d-02 /
      data f( 29), f( 30) /  2.0592570000d-02 ,  1.7550060000d-02 /
      data f( 31), f( 32) /  1.5154660000d-02 ,  1.3225980000d-02 /
      data f( 33), f( 34) /  1.1635130000d-02 ,  1.0291520000d-02 /
      data f( 35), f( 36) /  9.1325200000d-03 ,  8.1154600000d-03 /
      data f( 37), f( 38) /  7.2115500000d-03 ,  6.4013300000d-03 /
      data f( 39), f( 40) /  5.6715100000d-03 ,  5.0126300000d-03 /
      data f( 41), f( 42) /  4.4176900000d-03 ,  3.8810900000d-03 /
      data f( 43), f( 44) /  3.3981300000d-03 ,  2.9646600000d-03 /
      data f( 45), f( 46) /  2.5768100000d-03 ,  2.2309600000d-03 /
      data f( 47), f( 48) /  1.9236500000d-03 ,  1.6515700000d-03 /
      data f( 49), f( 50) /  1.4115400000d-03 ,  1.2005700000d-03 /
      data f( 51), f( 52) /  8.5470000000d-04 ,  5.9369000000d-04 /
      data f( 53), f( 54) /  4.0017000000d-04 ,  2.5951000000d-04 /
      data f( 55), f( 56) /  1.5947000000d-04 ,  8.9920000000d-05 /
      data f( 57), f( 58) /  4.2780000000d-05 ,  1.1750000000d-05 /
      data f( 59), f( 60) / -7.9000000000d-06 , -1.9710000000d-05 /
      data f( 61), f( 62) / -2.9810000000d-05 , -2.7770000000d-05 /
      data f( 63), f( 64) / -2.2350000000d-05 , -1.6480000000d-05 /
      data f( 65), f( 66) / -1.1370000000d-05 , -7.3700000000d-06 /
      data f( 67), f( 68) / -2.9200000000d-06 , -1.3000000000d-06 /
      data f( 69), f( 70) / -7.0000000000d-07 , -2.1000000000d-07 /
      data f( 71), f( 72) / -5.0000000000d-08 ,  1.0000000000d-08 /
 
      data x(  1), x(  2) /  1.0000000000d-01 ,  1.5000000000d-01 /
      data x(  3), x(  4) /  2.0000000000d-01 ,  2.5000000000d-01 /
      data x(  5), x(  6) /  3.0000000000d-01 ,  3.5000000000d-01 /
      data x(  7), x(  8) /  4.0000000000d-01 ,  4.5000000000d-01 /
      data x(  9), x( 10) /  5.0000000000d-01 ,  5.5000000000d-01 /
      data x( 11), x( 12) /  6.0000000000d-01 ,  6.5000000000d-01 /
      data x( 13), x( 14) /  7.0000000000d-01 ,  7.5000000000d-01 /
      data x( 15), x( 16) /  8.0000000000d-01 ,  8.5000000000d-01 /
      data x( 17), x( 18) /  9.0000000000d-01 ,  9.5000000000d-01 /
      data x( 19), x( 20) /  1.0000000000d+00 ,  1.0500000000d+00 /
      data x( 21), x( 22) /  1.1000000000d+00 ,  1.2000000000d+00 /
      data x( 23), x( 24) /  1.3000000000d+00 ,  1.4000000000d+00 /
      data x( 25), x( 26) /  1.5000000000d+00 ,  1.6000000000d+00 /
      data x( 27), x( 28) /  1.7000000000d+00 ,  1.8000000000d+00 /
      data x( 29), x( 30) /  1.9000000000d+00 ,  2.0000000000d+00 /
      data x( 31), x( 32) /  2.1000000000d+00 ,  2.2000000000d+00 /
      data x( 33), x( 34) /  2.3000000000d+00 ,  2.4000000000d+00 /
      data x( 35), x( 36) /  2.5000000000d+00 ,  2.6000000000d+00 /
      data x( 37), x( 38) /  2.7000000000d+00 ,  2.8000000000d+00 /
      data x( 39), x( 40) /  2.9000000000d+00 ,  3.0000000000d+00 /
      data x( 41), x( 42) /  3.1000000000d+00 ,  3.2000000000d+00 /
      data x( 43), x( 44) /  3.3000000000d+00 ,  3.4000000000d+00 /
      data x( 45), x( 46) /  3.5000000000d+00 ,  3.6000000000d+00 /
      data x( 47), x( 48) /  3.7000000000d+00 ,  3.8000000000d+00 /
      data x( 49), x( 50) /  3.9000000000d+00 ,  4.0000000000d+00 /
      data x( 51), x( 52) /  4.2000000000d+00 ,  4.4000000000d+00 /
      data x( 53), x( 54) /  4.6000000000d+00 ,  4.8000000000d+00 /
      data x( 55), x( 56) /  5.0000000000d+00 ,  5.2000000000d+00 /
      data x( 57), x( 58) /  5.4000000000d+00 ,  5.6000000000d+00 /
      data x( 59), x( 60) /  5.8000000000d+00 ,  6.0000000000d+00 /
      data x( 61), x( 62) /  6.5000000000d+00 ,  7.0000000000d+00 /
      data x( 63), x( 64) /  7.5000000000d+00 ,  8.0000000000d+00 /
      data x( 65), x( 66) /  8.5000000000d+00 ,  9.0000000000d+00 /
      data x( 67), x( 68) /  1.0000000000d+01 ,  1.1000000000d+01 /
      data x( 69), x( 70) /  1.2000000000d+01 ,  1.5000000000d+01 /
      data x( 71), x( 72) /  2.0000000000d+01 ,  3.0000000000d+01 /
 
      data del(  1), del(  2) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del(  3), del(  4) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del(  5), del(  6) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del(  7), del(  8) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del(  9), del( 10) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 11), del( 12) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 13), del( 14) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 15), del( 16) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 17), del( 18) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 19), del( 20) /  5.0000000000d-02 ,  5.0000000000d-02 /
      data del( 21), del( 22) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 23), del( 24) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 25), del( 26) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 27), del( 28) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 29), del( 30) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 31), del( 32) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 33), del( 34) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 35), del( 36) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 37), del( 38) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 39), del( 40) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 41), del( 42) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 43), del( 44) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 45), del( 46) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 47), del( 48) /  1.0000000000d-01 ,  1.0000000000d-01 /
      data del( 49), del( 50) /  1.0000000000d-01 ,  2.0000000000d-01 /
      data del( 51), del( 52) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del( 53), del( 54) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del( 55), del( 56) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del( 57), del( 58) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del( 59), del( 60) /  2.0000000000d-01 ,  5.0000000000d-01 /
      data del( 61), del( 62) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del( 63), del( 64) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del( 65), del( 66) /  5.0000000000d-01 ,  1.0000000000d+00 /
      data del( 67), del( 68) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 69), del( 70) /  3.0000000000d+00 ,  5.0000000000d+00 /
      data del( 71) /  1.0000000000d+01 /
      data npts / 72/
 

      if(xi .le. x(1)) then
        ii=1
      else if( xi .ge. x(npts)) then
        ii=npts-1
      else
        call hunt(x,npts,xi,ii)
      endif
      
 20   fi = fpp(ii)   * (x(ii+1)-xi)**3 / (6.0d0*del(ii)) + 
     x     fpp(ii+1) * (xi-x(ii))**3   / (6.0d0*del(ii)) +
     x   ((f(ii+1)/del(ii))-(fpp(ii+1)*del(ii)/6.0d0)) * (xi-x(ii)) + 
     x   ((f(ii)  /del(ii))-(fpp(ii)  *del(ii)/6.0d0)) * (x(ii+1)-xi)
      return
      end
