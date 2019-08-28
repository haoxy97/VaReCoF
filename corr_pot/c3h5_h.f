
      real*8 function c3h5_h(r, rpar, ipar)

c
c
c     analytic 3d fit for h+allyl
c
c     coordinates:            hi
c                            /
c                        ha /
c                        | /
c                        c1
c                       / \
c                   hb-c2  c3-hd
c                      |   |
c                      hc  he
c
c       r     = r(hi-c1)                    (au)     4.0 < r     < 12.0
c       theta = hi-c1-ha angle            (degrees)  0.0 < theta < 180
c       phi   = hi-c1-hac2 dihedral angle (degrees)       periodic
c
c     output: energy (kcal/mole) relative to h+allyl
c
c     note: the potential has four-fold symmetry in phi
c           i.e. energy(r,theta,phi)=energy(r,theta,phi+90) 
c
c     cartesian coordinates for allyl radical (au):
c
c                x       y                    z
c
c     c1         0.0     0.0                  0.0
c     c2         0.0     2.3398111211d+00    -1.2440996414d+00 
c     c3         0.0    -2.3398111211d+00    -1.2440996414d+00 
c     ha         0.0     0.0                  2.05
c     hb         0.0     4.0783097182d+00    -1.5776514970d-01 
c     hc         0.0     2.4113550893d+00    -3.2928508368d+00 /
c     hd         0.0    -4.0783097182d+00    -1.5776514970d-01 
c     he         0.0    -2.4113550893d+00    -3.2928508368d+00 /
c
c
c     energy - potential energy relative to h+c3h5 in au

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

c     ichab = ipotpar(iel,ieff,ich,1)
      phidvd = 90.0d0*pi/180.0d0
c     ibsc = ipotpar(iel,ieff,ich,2)
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
      rhh = sqrt((r(2,1,1)-r(1,3,1))**2 + (r(2,1,2)-r(1,3,2))**2 + 
     $ (r(2,1,3)-r(1,3,3))**2)
	rtest = sqrt((r(2,1,1)-r(1,4,1))**2 + (r(2,1,2)-r(1,4,2))**2 + 
     $ (r(2,1,3)-r(1,4,3))**2)

      cthe = (rr**2 + rch**2 - rhh**2)/(2.0d0*rr*rch)
      cthep = (rchp**2 + rch**2 - r122)/(2.0d0*rchp*rch)
!      cthep = (rcc**2 + rch**2 - r132)/(2.0d0*rcc*rch)
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
      call cross(r21,r23,tmp1)
      call cross(r23,r13,tmp2)
      tdot = 0.0d0
      do 1200 idim = 1 , ndim
         tdot = tdot + tmp1(idim)*tmp2(idim)
1200  continue
      rdenom = r232*sqrt(r212*r132)
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
      call cross(r23,tmp2,tmp3)
      do 2570 idim = 1 , ndim
         sum = sum - tmp1(idim)*tmp3(idim)
 2570 continue
      ctaub = sum/(r232*sqrt(r212*r232*r132)*sin(the)*sin(thep))
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
	theta = the*180.0d0/pi
      ichtmp = 1
      if (abs(phi).gt.90.0d0) ichtmp = 2
	if (tau.lt.0.0d0) ichtmp = ichtmp + 2
c ichtmp = 1 ==> l,+
c ichtmp = 2 ==> r,+
c ichtmp = 3 ==> l,-
c ichtmp = 4 ==> r,-
c ichab = 0 ==> all
c ichab = 1 ==> l,+ and l,-
c ichab = 2 ==> l,+ and r,+
c ichab = 3 ==> l,+ only
c      if ((ichab.eq.1).and.((ichtmp.eq.2).or.(ichtmp.eq.4))) then
c   vtot = 1.0e20
c         go to 5000
c	endif
c      if ((ichab.eq.2).and.((ichtmp.eq.3).or.(ichtmp.eq.4))) then
c	   vtot = 1.0e20
c         go to 5000
c	endif
c	if ((ichab.eq.3).and.(ichtmp.ne.1)) then
c	   vtot = 1.0e20
c         go to 5000
c	endif
      if (rr.lt.3.5d0) then
	   rr = 3.50001
c        call hpc3h5(rr,theta,phi,energy)
c        vtot = energy*1.59360124d-03
c   if (vtot.lt.0.0d0) write (6,*) 
c    $    'error in rr',rr*cautoang,vtot*cautoicm
	endif
      if (rr.gt.12.0d0) then
c   write (6,*) 'error in rr',rr*cautoang
	   rr = 11.9999
	endif
c	if (((rmepi.lt.rtest).and.(phi.gt.90.0d0)).or.
c     $ ((rmepi.gt.rtest).and.(phi.lt.90.0d0))) then
c	  write (6,*) 'pot test',rmepi,rtest,theta,phi,vtot*cautoicm
c	write (6,*) 'r1 geom',(r(1,1,idim)*cautoang,idim=1,3)
c	write (6,*) 'r2 geom',(r(1,2,idim)*cautoang,idim=1,3)
c	write (6,*) 'r3 geom',(r(1,3,idim)*cautoang,idim=1,3)
c	write (6,*) 'r4 geom',(r(1,4,idim)*cautoang,idim=1,3)
c	write (6,*) 'r5 geom',(r(2,1,idim)*cautoang,idim=1,3)
c	endif
      call hpc3h5(rr,theta,phi,energy)
      vtot = energy*1.59360124d-03
c	if (vtot.lt.vtotmin) then 
c	   vtotmin = vtot
c	   write (6,*) 'vtot test',vtot*cautoicm,rr,theta,phi,
c     $   rmepi*cautoang,rtest*cautoang
c	endif
c      call hpc3h5(rr,theta,-phi,energyp)
5000  continue  
c	write (6,*) 'pot test',rr,rtest,theta,phi,vtot*cautoicm
c      vtotp = energyp*1.59360124d-03
c	write (6,*) 'pot test',vtot*cautoicm,vtotp*cautoicm,ichab

      c3h5_h = vtot
      return
      end


      block data
      implicit real*8 (a-h,o-z)

      common/ri/x(12)
      common/es1/y1(12)
      common/es2/y2(12)
      common/info/nri

c
c     nri = number of points in x, y1 and y2
c
c     x values of r for 1d interpolations
c     y1(i) = mep energy for x(i)
c     y2(i) = scale point energy for x(i)
c
c     in this case the scale point was chosen to
c     be theta = 90.0 and phi = 90.0
c
      data x  /   3.5d0,  4.0d0,
     x            4.5d0,  5.0d0,  5.5d0,  6.0d0, 
     x            6.5d0,  7.0d0,  8.0d0,  9.0d0,
     x           10.0d0, 12.0d0  /

      data y1 / -59.29d0, -56.25d0,
     x          -42.62d0, -26.52d0,  -15.25d0, -8.219d0,
     x           -4.278d0, -2.212d0, -0.613d0, -0.181d0,
     x           -0.0571d0,-0.0098d0 /

      data y2 /  21.019d0, 10.880d0,
     x            5.435d0,  2.611d0,   1.190d0, 0.501d0,
     x            0.182d0,  0.0451d0, -0.0212d0, -0.0148d0,
     x           -0.0063d0,-0.0015d0 /

      data nri / 12 /
      end

      subroutine hpc3h5(rin,thetain,phiin,energy)
      implicit real*8 (a-h,o-z)
      common/ri/x(12)
      common/es1/y1(12)
      common/es2/y2(12)
      data degrad / 57.29577951308232d 00/
c
c     cartesian coordinates for the h atoms of c3h5
c
      data haz / 2.0500000000d+00 /
      data hby, hbz / 4.0783097182d+00, -1.5776514970d-01 /
      data hcy, hcz / 2.4113550893d+00, -3.2928508368d+00 /
c
c     cartesian coordinates for the c atoms of c3h5
c
      data cy, cz / 2.3398111211d+00,   -1.2440996414d+00 /


      xi = rin*sin(thetain/degrad) * sin(phiin/degrad)
      yi = rin*sin(thetain/degrad) * cos(phiin/degrad)
      zi = rin*cos(thetain/degrad)

c
c     ra,rb,rc,rd,re are the distances between the
c     active hydrogen atom and each of the hydrogens
c     of c3h5
c
      ra = dsqrt(xi**2 + yi**2 + (zi-haz)**2)
      rb = dsqrt(xi**2 + (yi-hby)**2 + (zi-hbz)**2)
      rc = dsqrt(xi**2 + (yi-hcy)**2 + (zi-hcz)**2)
      rd = dsqrt(xi**2 + (yi+hby)**2 + (zi-hbz)**2)
      re = dsqrt(xi**2 + (yi+hcy)**2 + (zi-hcz)**2)

      call hhspl09(ra,ea)
      call hhspl09(rb,eb)
      call hhspl09(rc,ec)
      call hhspl09(rd,ed)
      call hhspl09(re,ee)


c
c     rc1, rc2 rc3 are the distances between the
c     active hydrogen atom and each of the carbons
c     of c3h5
c
      rc1 = dsqrt(xi**2 + yi**2 + (zi)**2)
      rc2 = dsqrt(xi**2 + (yi-cy)**2 + (zi-cz)**2)
      rc3 = dsqrt(xi**2 + (yi+cy)**2 + (zi-cz)**2)

      call escale(rin,es1,es2)

      if(rin .ge. x(12)) then
         call r2dp12p0(thetain,phiin,y1(12),y2(12),fi)
      else if(rin .le. x(3) ) then
         call r2dp4p5(thetain,phiin,y1(3),y2(3),fi)
      else 
        if(rin .gt. x(3) .and. rin .lt. x(4) ) then
           call r2dp4p5(thetain,phiin,y1(3),y2(3),f1)
           call r2dp5p0(thetain,phiin,y1(4),y2(4),f2)
           r1= x(3)
           r2= x(4)
        else if(rin .ge. x(4) .and. rin .lt. x(5) ) then
           call r2dp5p0(thetain,phiin,y1(4),y2(4),f1)
           call r2dp5p5(thetain,phiin,y1(5),y2(5),f2)
           r1= x(4)
           r2= x(5)
        else if(rin .ge. x(5) .and. rin .lt. x(6) ) then
           call r2dp5p5(thetain,phiin,y1(5),y2(5),f1)
           call r2dp6p0(thetain,phiin,y1(6),y2(6),f2)
           r1= x(5)
           r2= x(6)
        else if(rin .ge. x(6) .and. rin .lt. x(7) ) then
           call r2dp6p0(thetain,phiin,y1(6),y2(6),f1)
           call r2dp6p5(thetain,phiin,y1(7),y2(7),f2)
           r1= x(6)
           r2= x(7)
        else if(rin .ge. x(7) .and. rin .lt. x(8) ) then
           call r2dp6p5(thetain,phiin,y1(7),y2(7),f1)
           call r2dp7p0(thetain,phiin,y1(8),y2(8),f2)
           r1= x(7)
           r2= x(8)
        else if(rin .ge. x(8) .and. rin .lt. x(9) ) then
           call r2dp7p0(thetain,phiin,y1(8),y2(8),f1)
           call r2dp8p0(thetain,phiin,y1(9),y2(9),f2)
           r1= x(8)
           r2= x(9)
        else if(rin .ge. x(9) .and. rin .lt. x(10) ) then
           call r2dp8p0(thetain,phiin,y1( 9),y2( 9),f1)
           call r2dp9p0(thetain,phiin,y1(10),y2(10),f2)
           r1= x(9)
           r2= x(10)
        else if(rin .ge. x(10) .and. rin .lt. x(11) ) then
           call r2dp9p0 (thetain,phiin,y1(10),y2(10),f1)
           call r2dp10p0(thetain,phiin,y1(11),y2(11),f2)
           r1=  x(10)
           r2=  x(11)
        else if(rin .ge. x(11) .and. rin .lt. x(12) ) then
           call r2dp10p0(thetain,phiin,y1(11),y2(11),f1)
           call r2dp12p0(thetain,phiin,y1(12),y2(12),f2)
           r1=  x(11)
           r2=  x(12)
        endif
        fi = (f1*rin-f2*rin+r1*f2-r2*f1)/(r1-r2)
      endif

      energy =  (fi)*(es2 -es1) + es1 + ea + eb + ec + ed + ee

c     add basis set corrections based on the c-h distances

      call ebasis09(rc1,ec1)
      call ebasis09(rc2,ec2)
      call ebasis09(rc3,ec3)
c     energy = energy    + ec2 + ec3
c	write (6,*) 'energy test',energy,ec2,ec3
      energy = energy + ec3

      if (ec2.lt.ec3) energy = energy + ec2 - ec3

      return
      end

      subroutine escale(xi,e1,e2)
      implicit real*8 (a-h,o-z)
      call escale1(xi,e1)
      call escale2(xi,e2)
      return
      end  

      subroutine escale1(xi,fi)
c
c     1d spline for h+c3h5 mep energy
c
      implicit real*8 (a-h,o-z)
      dimension x( 13), f( 13),fpp( 13),del( 12)
 
      data fpp(  1), fpp(  2) /  7.5048374002d+01 ,  4.1943251996d+01 /
      data fpp(  3), fpp(  4) /  1.1338618012d+01 , -2.8017724046d+01 /
      data fpp(  5), fpp(  6) / -1.5187721828d+01 , -1.2967388642d+01 /
      data fpp(  7), fpp(  8) / -7.1027236059d+00 , -3.6217169349d+00 /
      data fpp(  9), fpp( 10) / -7.8148739222d-01 , -2.5433349619d-01 /
      data fpp( 11), fpp( 12) / -4.9778623030d-02 , -2.4247382817d-02 /
      data fpp( 13) /  7.3243691408d-02 /
 
      data f(  1), f(  2) / -5.9290000000d+01 , -5.6250000000d+01 /
      data f(  3), f(  4) / -4.2620000000d+01 , -2.6520000000d+01 /
      data f(  5), f(  6) / -1.5250000000d+01 , -8.2190000000d+00 /
      data f(  7), f(  8) / -4.2780000000d+00 , -2.2120000000d+00 /
      data f(  9), f( 10) / -6.1300000000d-01 , -1.8100000000d-01 /
      data f( 11), f( 12) / -5.7100000000d-02 , -9.8000000000d-03 /
      data f( 13) /  0.0000000000d+00 /
 
      data x(  1), x(  2) /  3.5000000000d+00 ,  4.0000000000d+00 /
      data x(  3), x(  4) /  4.5000000000d+00 ,  5.0000000000d+00 /
      data x(  5), x(  6) /  5.5000000000d+00 ,  6.0000000000d+00 /
      data x(  7), x(  8) /  6.5000000000d+00 ,  7.0000000000d+00 /
      data x(  9), x( 10) /  8.0000000000d+00 ,  9.0000000000d+00 /
      data x( 11), x( 12) /  1.0000000000d+01 ,  1.2000000000d+01 /
      data x( 13) /  1.5000000000d+01 /
 
      data del(  1), del(  2) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  3), del(  4) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  5), del(  6) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  7), del(  8) /  5.0000000000d-01 ,  1.0000000000d+00 /
      data del(  9), del( 10) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 11), del( 12) /  2.0000000000d+00 ,  3.0000000000d+00 /
      data npts / 13/

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
      subroutine escale2(xi,fi)
c
c     1d spline for h+c3h5 scale (90,90) energy
c
      implicit real*8 (a-h,o-z)
      dimension x( 15), f( 15),fpp( 15),del( 14)
 
      data fpp(  1), fpp(  2) /  6.7946075351d+01 ,  5.0368249298d+01 /
      data fpp(  3), fpp(  4) /  3.0981727458d+01 ,  1.7925640872d+01 /
      data fpp(  5), fpp(  6) /  9.9381090562d+00 ,  5.2523229035d+00 /
      data fpp(  7), fpp(  8) /  2.7269993296d+00 ,  1.3788797780d+00 /
      data fpp(  9), fpp( 10) /  6.6868155849d-01 ,  2.9279398807d-01 /
      data fpp( 11), fpp( 12) /  3.7077256545d-02 , -4.9030142484d-03 /
      data fpp( 13), fpp( 14) / -4.8651995508d-03 , -1.2528942233d-03 /
      data fpp( 15) /  3.6197804450d-03 /
 
      data f(  1), f(  2) /  6.9830300000d+01 ,  3.9166100000d+01 /
      data f(  3), f(  4) /  2.1018600000d+01 ,  1.0880300000d+01 /
      data f(  5), f(  6) /  5.4346000000d+00 ,  2.6110000000d+00 /
      data f(  7), f(  8) /  1.1905000000d+00 ,  5.0080000000d-01 /
      data f(  9), f( 10) /  1.8240000000d-01 ,  4.5100000000d-02 /
      data f( 11), f( 12) / -2.1200000000d-02 , -1.4800000000d-02 /
      data f( 13), f( 14) / -6.3000000000d-03 , -1.5000000000d-03 /
      data f( 15) /  0.0000000000d+00 /
 
      data x(  1), x(  2) /  2.5000000000d+00 ,  3.0000000000d+00 /
      data x(  3), x(  4) /  3.5000000000d+00 ,  4.0000000000d+00 /
      data x(  5), x(  6) /  4.5000000000d+00 ,  5.0000000000d+00 /
      data x(  7), x(  8) /  5.5000000000d+00 ,  6.0000000000d+00 /
      data x(  9), x( 10) /  6.5000000000d+00 ,  7.0000000000d+00 /
      data x( 11), x( 12) /  8.0000000000d+00 ,  9.0000000000d+00 /
      data x( 13), x( 14) /  1.0000000000d+01 ,  1.2000000000d+01 /
      data x( 15) /  1.5000000000d+01 /
 
      data del(  1), del(  2) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  3), del(  4) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  5), del(  6) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  7), del(  8) /  5.0000000000d-01 ,  5.0000000000d-01 /
      data del(  9), del( 10) /  5.0000000000d-01 ,  1.0000000000d+00 /
      data del( 11), del( 12) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 13), del( 14) /  2.0000000000d+00 ,  3.0000000000d+00 /
      data npts / 15/
 

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
      subroutine ebasis09(xi,fi)
c
c     1d spline for h+ch3 
c     difference between aug-pvtz and pvdz calculations
c     correction scaled by mep(h+c3h5)/mep(h+ch3)
c
      implicit real*8 (a-h,o-z)
      dimension x( 48), f( 48),fpp( 48),del( 47)
 
      data fpp(  1), fpp(  2) /  4.5556716179d-03 , -9.1113432356d-03 /
      data fpp(  3), fpp(  4) /  3.1889701325d-02 , -1.1844746206d-01 /
      data fpp(  5), fpp(  6) /  4.2326330849d-01 ,  3.5475569117d-01 /
      data fpp(  7), fpp(  8) /  1.5981571401d+00 ,  1.1639730891d+00 /
      data fpp(  9), fpp( 10) /  1.6146296887d-01 ,  5.1670880292d-02 /
      data fpp( 11), fpp( 12) /  1.9695046287d-01 , -1.9127328578d-01 /
      data fpp( 13), fpp( 14) /  2.0249171072d-01 , -8.0990357440d-01 /
      data fpp( 15), fpp( 16) / -3.4125874274d-01 , -3.8461824133d-01 /
      data fpp( 17), fpp( 18) / -3.3079460772d-01 , -3.6976565467d-01 /
      data fpp( 19), fpp( 20) / -2.5108460184d-01 , -1.5498181054d-01 /
      data fpp( 21), fpp( 22) / -2.0932056597d-01 , -2.7089382033d-01 /
      data fpp( 23), fpp( 24) / -1.1984653498d-01 , -2.3305000784d-01 /
      data fpp( 25), fpp( 26) / -1.4803806140d-02 , -9.0006235747d-02 /
      data fpp( 27), fpp( 28) / -7.3924713476d-02 , -6.3048372952d-02 /
      data fpp( 29), fpp( 30) / -5.6153262859d-02 , -4.4748547912d-02 /
      data fpp( 31), fpp( 32) / -3.0780523333d-02 , -3.1575342137d-02 /
      data fpp( 33), fpp( 34) / -1.8011598534d-02 , -4.4903763943d-02 /
      data fpp( 35), fpp( 36) /  9.7903662618d-02 ,  2.3195843858d-03 /
      data fpp( 37), fpp( 38) / -7.4590084715d-03 , -5.7245477296d-03 /
      data fpp( 39), fpp( 40) / -2.8837978401d-03 , -5.5038004504d-03 /
      data fpp( 41), fpp( 42) / -2.8430388828d-03 , -1.1240440186d-03 /
      data fpp( 43), fpp( 44) / -1.6607850430d-03 ,  2.3671841905d-03 /
      data fpp( 45), fpp( 46) / -3.6079517189d-03 ,  6.6462268506d-04 /
      data fpp( 47), fpp( 48) / -1.8989219573d-04 ,  9.4946097865d-05 /
 
      data f(  1), f(  2) / -1.6343000000d+00 , -1.6343000000d+00 /
      data f(  3), f(  4) / -1.6343000000d+00 , -1.6343000000d+00 /
      data f(  5), f(  6) / -1.6343000000d+00 , -1.6220000000d+00 /
      data f(  7), f(  8) / -1.5890000000d+00 , -1.5084000000d+00 /
      data f(  9), f( 10) / -1.3892000000d+00 , -1.2588000000d+00 /
      data f( 11), f( 12) / -1.1250000000d+00 , -9.8730000000d-01 /
      data f( 13), f( 14) / -8.5180000000d-01 , -7.2440000000d-01 /
      data f( 15), f( 16) / -6.0980000000d-01 , -5.1060000000d-01 /
      data f( 17), f( 18) / -4.2470000000d-01 , -3.5130000000d-01 /
      data f( 19), f( 20) / -2.9030000000d-01 , -2.3850000000d-01 /
      data f( 21), f( 22) / -1.9320000000d-01 , -1.5550000000d-01 /
      data f( 23), f( 24) / -1.2630000000d-01 , -1.0410000000d-01 /
      data f( 25), f( 26) / -8.6900000000d-02 , -7.2000000000d-02 /
      data f( 27), f( 28) / -5.9800000000d-02 , -5.0300000000d-02 /
      data f( 29), f( 30) / -4.3100000000d-02 , -3.7900000000d-02 /
      data f( 31), f( 32) / -3.4300000000d-02 , -3.1900000000d-02 /
      data f( 33), f( 34) / -3.0600000000d-02 , -3.0100000000d-02 /
      data f( 35), f( 36) / -3.0200000000d-02 , -2.8200000000d-02 /
      data f( 37), f( 38) / -2.5600000000d-02 , -2.3200000000d-02 /
      data f( 39), f( 40) / -2.1000000000d-02 , -1.5400000000d-02 /
      data f( 41), f( 42) / -8.8000000000d-03 , -5.2000000000d-03 /
      data f( 43), f( 44) / -3.1000000000d-03 , -1.9000000000d-03 /
      data f( 45), f( 46) /  0.0000000000d+00 ,  0.0000000000d+00 /
      data f( 47), f( 48) /  0.0000000000d+00 ,  0.0000000000d+00 /
 
      data x(  1), x(  2) /  2.2000000000d+00 ,  2.4000000000d+00 /
      data x(  3), x(  4) /  2.6000000000d+00 ,  2.8000000000d+00 /
      data x(  5), x(  6) /  3.0200000000d+00 ,  3.2100000000d+00 /
      data x(  7), x(  8) /  3.4000000000d+00 ,  3.5900000000d+00 /
      data x(  9), x( 10) /  3.7800000000d+00 ,  3.9700000000d+00 /
      data x( 11), x( 12) /  4.1600000000d+00 ,  4.3500000000d+00 /
      data x( 13), x( 14) /  4.5400000000d+00 ,  4.7200000000d+00 /
      data x( 15), x( 16) /  4.9100000000d+00 ,  5.1000000000d+00 /
      data x( 17), x( 18) /  5.2900000000d+00 ,  5.4800000000d+00 /
      data x( 19), x( 20) /  5.6700000000d+00 ,  5.8600000000d+00 /
      data x( 21), x( 22) /  6.0500000000d+00 ,  6.2400000000d+00 /
      data x( 23), x( 24) /  6.4300000000d+00 ,  6.6100000000d+00 /
      data x( 25), x( 26) /  6.8000000000d+00 ,  6.9900000000d+00 /
      data x( 27), x( 28) /  7.1800000000d+00 ,  7.3700000000d+00 /
      data x( 29), x( 30) /  7.5600000000d+00 ,  7.7500000000d+00 /
      data x( 31), x( 32) /  7.9400000000d+00 ,  8.1300000000d+00 /
      data x( 33), x( 34) /  8.3100000000d+00 ,  8.5000000000d+00 /
      data x( 35), x( 36) /  8.6900000000d+00 ,  8.8800000000d+00 /
      data x( 37), x( 38) /  9.0700000000d+00 ,  9.2600000000d+00 /
      data x( 39), x( 40) /  9.4500000000d+00 ,  1.0000000000d+01 /
      data x( 41), x( 42) /  1.1000000000d+01 ,  1.2000000000d+01 /
      data x( 43), x( 44) /  1.3000000000d+01 ,  1.4000000000d+01 /
      data x( 45), x( 46) /  1.5000000000d+01 ,  1.6000000000d+01 /
      data x( 47), x( 48) /  1.8000000000d+01 ,  2.0000000000d+01 /
 
      data del(  1), del(  2) /  2.0000000000d-01 ,  2.0000000000d-01 /
      data del(  3), del(  4) /  2.0000000000d-01 ,  2.2000000000d-01 /
      data del(  5), del(  6) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del(  7), del(  8) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del(  9), del( 10) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 11), del( 12) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 13), del( 14) /  1.8000000000d-01 ,  1.9000000000d-01 /
      data del( 15), del( 16) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 17), del( 18) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 19), del( 20) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 21), del( 22) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 23), del( 24) /  1.8000000000d-01 ,  1.9000000000d-01 /
      data del( 25), del( 26) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 27), del( 28) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 29), del( 30) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 31), del( 32) /  1.9000000000d-01 ,  1.8000000000d-01 /
      data del( 33), del( 34) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 35), del( 36) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 37), del( 38) /  1.9000000000d-01 ,  1.9000000000d-01 /
      data del( 39), del( 40) /  5.5000000000d-01 ,  1.0000000000d+00 /
      data del( 41), del( 42) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 43), del( 44) /  1.0000000000d+00 ,  1.0000000000d+00 /
      data del( 45), del( 46) /  1.0000000000d+00 ,  2.0000000000d+00 /
      data del( 47) /  2.0000000000d+00 /
      data npts / 48/
 

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

      subroutine hhspl09(xi,fi)
c
c     h-h nonbonded potential from:
c     ch-h 3b1 cas+1+2/cc-pvdz
c     xi = hh bond distance in au
c     binding energy returned in kcal/mole
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
      fi = fi *627.51d 00
      return
      end  

      subroutine r2dp12p0(thetain,phiin,es1,es2,fi)
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(20)
      dimension cc(8)
c
      data fk(   1), fk(   2) /  1.2751977758d-03 , -1.0603110583d-03 /
      data fk(   3), fk(   4) /  7.9888379711d-05 ,  6.0396900059d-04 /
      data fk(   5), fk(   6) /  1.7055329072d-04 , -1.7307415079d-04 /
      data fk(   7), fk(   8) / -6.4771659396d-05 ,  1.4391074431d-05 /
      data fk(   9), fk(  10) /  8.7605898946d-06 , -2.3340626082d-06 /
      data fk(  11), fk(  12) /  1.2735688962d-05 , -2.5292965900d-06 /
      data fk(  13), fk(  14) / -4.4928160855d-05 , -4.1129501963d-06 /
      data fk(  15), fk(  16) /  3.3399076870d-07 ,  2.2314889204d-08 /
      data fk(  17), fk(  18) / -1.8826572348d-05 , -4.3650461340d-07 /
      data fk(  19), fk(  20) /  3.1793142944d-08 , -6.6563389157d-09 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )

      en = 0.0d 0
      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp10p0(thetain,phiin,es1,es2,fi)
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(25)
      dimension cc(8)
c
      data fk(   1), fk(   2) /  2.9308907233d-03 , -3.1184441622d-03 /
      data fk(   3), fk(   4) /  8.7274662093d-03 ,  4.0094006216d-03 /
      data fk(   5), fk(   6) / -7.5524083022d-03 , -1.0314301783d-03 /
      data fk(   7), fk(   8) / -5.0470860902d-04 ,  2.3576362964d-04 /
      data fk(   9), fk(  10) /  1.7058368651d-04 ,  3.3208588841d-03 /
      data fk(  11), fk(  12) / -5.9547189155d-05 , -4.1502034993d-05 /
      data fk(  13), fk(  14) / -2.7147847373d-03 ,  2.4539274006d-05 /
      data fk(  15), fk(  16) /  5.7429093441d-06 ,  6.8665717816d-07 /
      data fk(  17), fk(  18) /  9.8959753508d-04 ,  4.2402077637d-06 /
      data fk(  19), fk(  20) / -3.2665136880d-07 , -1.3420349329d-07 /
      data fk(  21), fk(  22) /  8.8467681496d-04 , -2.6432765949d-05 /
      data fk(  23), fk(  24) /  1.2070126398d-07 ,  1.5267209376d-08 /
      data fk(  25) /  8.7232425533d-10 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )

      en = 0.0d 00
      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 


      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp9p0(thetain,phiin,es1,es2,fi)      
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(25)
      dimension cc(8)
c
      data fk(   1), fk(   2) / -8.0861934304d-04 , -4.9047560669d-03 /
      data fk(   3), fk(   4) /  3.6916556535d-02 ,  1.2521616871d-02 /
      data fk(   5), fk(   6) / -3.0148894355d-02 , -3.9444205423d-03 /
      data fk(   7), fk(   8) / -2.3345397234d-03 ,  9.6503979015d-04 /
      data fk(   9), fk(  10) /  7.5259917972d-04 ,  1.4080531040d-02 /
      data fk(  11), fk(  12) / -4.1215814631d-04 , -1.8258017759d-04 /
      data fk(  13), fk(  14) / -1.3456308111d-02 ,  2.3302740195d-04 /
      data fk(  15), fk(  16) /  2.7042672047d-05 ,  3.2550522649d-06 /
      data fk(  17), fk(  18) /  3.9547146590d-03 ,  4.5882456424d-05 /
      data fk(  19), fk(  20) / -2.1072190470d-06 , -6.4446964237d-07 /
      data fk(  21), fk(  22) /  4.9785787188d-03 , -1.2906317109d-04 /
      data fk(  23), fk(  24) /  4.4567521483d-07 ,  8.6492682237d-08 /
      data fk(  25) /  5.2930051075d-09 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp8p0(thetain,phiin,es1,es2,fi)      
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(36)
      dimension cc(10)
c
      data fk(   1), fk(   2) / -5.0627793274d-03 , -2.6824657263d-02 /
      data fk(   3), fk(   4) /  1.3239347024d-01 ,  4.9763465536d-02 /
      data fk(   5), fk(   6) / -1.1184721690d-01 , -1.8912038411d-02 /
      data fk(   7), fk(   8) / -7.7177614191d-03 ,  4.8561561586d-03 /
      data fk(   9), fk(  10) /  3.3259304329d-03 ,  7.0953786701d-02 /
      data fk(  11), fk(  12) / -2.1397865619d-03 , -7.9507211607d-04 /
      data fk(  13), fk(  14) / -6.6640882752d-02 ,  1.2522594086d-03 /
      data fk(  15), fk(  16) /  1.1336178898d-04 ,  1.5666820106d-05 /
      data fk(  17), fk(  18) /  2.0447509918d-02 , -1.3999417395d-04 /
      data fk(  19), fk(  20) / -1.1002117565d-05 , -2.9950787459d-06 /
      data fk(  21), fk(  22) /  1.7081067057d-02 , -4.8219468476d-04 /
      data fk(  23), fk(  24) /  3.4252806466d-06 ,  3.9137200445d-07 /
      data fk(  25), fk(  26) /  2.2583010562d-08 , -1.8571151474d-02 /
      data fk(  27), fk(  28) /  3.5071028722d-04 , -1.1917996729d-06 /
      data fk(  29), fk(  30) / -3.4279522765d-08 , -3.3277177287d-09 /
      data fk(  31), fk(  32) /  3.7347747953d-03 , -6.8099754996d-05 /
      data fk(  33), fk(  34) /  9.0665152509d-08 ,  5.9281221531d-09 /
      data fk(  35), fk(  36) /  4.3473621513d-10 ,  1.2460809429d-11 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp7p0(thetain,phiin,es1,es2,fi)      
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(42)
      dimension cc(10)
c
      data fk(   1), fk(   2) /  6.3247682293d-02 , -2.0180625127d-01 /
      data fk(   3), fk(   4) /  4.7059966929d-01 ,  2.3755998130d-01 /
      data fk(   5), fk(   6) / -4.3413064804d-01 , -9.2854901562d-02 /
      data fk(   7), fk(   8) / -2.8509672400d-03 ,  2.4340115212d-02 /
      data fk(   9), fk(  10) /  1.4795182430d-02 ,  3.6550131834d-01 /
      data fk(  11), fk(  12) / -1.5314949799d-02 , -3.5456648294d-03 /
      data fk(  13), fk(  14) / -3.6786374256d-01 ,  9.7224879661d-03 /
      data fk(  15), fk(  16) /  5.3214586860d-04 ,  7.8451686359d-05 /
      data fk(  17), fk(  18) /  9.7669305375d-02 , -8.9476478867d-04 /
      data fk(  19), fk(  20) / -8.4831374991d-05 , -1.5052508455d-05 /
      data fk(  21), fk(  22) /  8.3831810625d-02 , -2.9502259470d-03 /
      data fk(  23), fk(  24) /  4.6970162881d-05 ,  1.9297307305d-06 /
      data fk(  25), fk(  26) /  1.4321891912d-07 , -8.0531212835d-02 /
      data fk(  27), fk(  28) /  1.8472065788d-03 , -1.6982539616d-05 /
      data fk(  29), fk(  30) / -2.5545288283d-07 , -2.2508175694d-08 /
      data fk(  31), fk(  32) /  1.4309633714d-02 , -2.4349407567d-04 /
      data fk(  33), fk(  34) / -7.3036609414d-07 ,  1.1295915738d-07 /
      data fk(  35), fk(  36) /  2.2960212382d-09 ,  1.4472498349d-10 /
      data fk(  37), fk(  38) /  8.7971932623d-03 , -1.4801555529d-04 /
      data fk(  39), fk(  40) /  2.4137967085d-06 , -4.1162983453d-08 /
      data fk(  41), fk(  42) / -2.5283845749d-10 , -1.8306557737d-11 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )

      en = 0.0d 00
      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 


      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp6p5(thetain,phiin,es1,es2,fi)      
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(69)
      dimension cc(14)
c
      data fk(   1), fk(   2) /  1.9726863362d-01 , -4.5601093056d-01 /
      data fk(   3), fk(   4) /  9.0818537259d-01 ,  4.9797146510d-01 /
      data fk(   5), fk(   6) / -9.0796020764d-01 , -1.9243850140d-01 /
      data fk(   7), fk(   8) /  5.2251615422d-02 ,  5.0805908016d-02 /
      data fk(   9), fk(  10) /  3.0450030564d-02 ,  7.9131673690d-01 /
      data fk(  11), fk(  12) / -3.7102632909d-02 , -7.2800214775d-03 /
      data fk(  13), fk(  14) / -8.3683238467d-01 ,  2.4318087612d-02 /
      data fk(  15), fk(  16) /  1.1029113919d-03 ,  1.6960009766d-04 /
      data fk(  17), fk(  18) /  2.3659727986d-01 , -2.7657261368d-03 /
      data fk(  19), fk(  20) / -2.1084998913d-04 , -3.2447476557d-05 /
      data fk(  21), fk(  22) /  1.6401097586d-01 , -6.2668040864d-03 /
      data fk(  23), fk(  24) /  1.1775600468d-04 ,  4.2677582842d-06 /
      data fk(  25), fk(  26) /  3.3239276174d-07 , -1.6013050646d-01 /
      data fk(  27), fk(  28) /  3.8741645001d-03 , -3.8963362627d-05 /
      data fk(  29), fk(  30) / -6.8593691531d-07 , -5.2372344258d-08 /
      data fk(  31), fk(  32) /  3.5365414763d-02 , -6.1062249382d-04 /
      data fk(  33), fk(  34) / -1.4883626222d-06 ,  3.0081257030d-07 /
      data fk(  35), fk(  36) /  5.7870667620d-09 ,  3.4926850131d-10 /
      data fk(  37), fk(  38) /  7.5036647234d-03 , -1.9981518811d-04 /
      data fk(  39), fk(  40) /  5.2334127999d-06 , -1.0380673942d-07 /
      data fk(  41), fk(  42) / -8.8386273128d-10 , -4.5857192674d-11 /
      data fk(  43), fk(  44) / -2.7860577608d-03 ,  8.1132723334d-05 /
      data fk(  45), fk(  46) / -1.6259258777d-06 ,  1.5271783715d-08 /
      data fk(  47), fk(  48) /  3.6671997236d-10 ,  4.3213493040d-12 /
      data fk(  49), fk(  50) /  2.4691934591d-13 ,  1.0633423374d-04 /
      data fk(  51), fk(  52) / -1.4628285674d-05 ,  2.3388300686d-07 /
      data fk(  53), fk(  54) /  1.7186045136d-09 , -1.2488322047d-10 /
      data fk(  55), fk(  56) / -6.2038615623d-13 , -3.0126023507d-14 /
      data fk(  57), fk(  58) / -8.4131704855d-04 ,  8.2705459868d-06 /
      data fk(  59), fk(  60) / -1.0635796232d-08 , -1.3315065121d-09 /
      data fk(  61), fk(  62) /  2.5343384501d-11 ,  2.4027425856d-13 /
      data fk(  63), fk(  64) /  2.8168215669d-15 ,  5.8018995800d-04 /
      data fk(  65), fk(  66) / -1.6909011520d-06 , -1.4853141957d-08 /
      data fk(  67), fk(  68) /  3.1133386409d-10 , -1.8273306884d-12 /
      data fk(  69) / -7.6552946723d-14 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )
      cc(12) = cos (12.0d 00 * theta )
      cc(14) = cos (14.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 
      en=en+fk(43)*plgndr(12,0,cxj) 
      en=en+fk(44)*plgndr(12,2,cxj)*cc(2) 
      en=en+fk(45)*plgndr(12,4,cxj)*cc(4) 
      en=en+fk(46)*plgndr(12,6,cxj)*cc(6) 
      en=en+fk(47)*plgndr(12,8,cxj)*cc(8) 
      en=en+fk(48)*plgndr(12,10,cxj)*cc(10) 
      en=en+fk(49)*plgndr(12,12,cxj)*cc(12) 
      en=en+fk(50)*plgndr(13,0,cxj) 
      en=en+fk(51)*plgndr(13,2,cxj)*cc(2) 
      en=en+fk(52)*plgndr(13,4,cxj)*cc(4) 
      en=en+fk(53)*plgndr(13,6,cxj)*cc(6) 
      en=en+fk(54)*plgndr(13,8,cxj)*cc(8) 
      en=en+fk(55)*plgndr(13,10,cxj)*cc(10) 
      en=en+fk(56)*plgndr(13,12,cxj)*cc(12) 
      en=en+fk(57)*plgndr(14,0,cxj) 
      en=en+fk(58)*plgndr(14,2,cxj)*cc(2) 
      en=en+fk(59)*plgndr(14,4,cxj)*cc(4) 
      en=en+fk(60)*plgndr(14,6,cxj)*cc(6) 
      en=en+fk(61)*plgndr(14,8,cxj)*cc(8) 
      en=en+fk(62)*plgndr(14,10,cxj)*cc(10) 
      en=en+fk(63)*plgndr(14,12,cxj)*cc(12) 
      en=en+fk(64)*plgndr(15,0,cxj) 
      en=en+fk(65)*plgndr(15,2,cxj)*cc(2) 
      en=en+fk(66)*plgndr(15,4,cxj)*cc(4) 
      en=en+fk(67)*plgndr(15,6,cxj)*cc(6) 
      en=en+fk(68)*plgndr(15,8,cxj)*cc(8) 
      en=en+fk(69)*plgndr(15,10,cxj)*cc(10) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp6p0(thetain,phiin,es1,es2,fi)      
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(69)
      dimension cc(14)
c
      data fk(   1), fk(   2) /  3.9846149962d-01 , -8.2678799208d-01 /
      data fk(   3), fk(   4) /  1.8071391922d+00 ,  9.1608947541d-01 /
      data fk(   5), fk(   6) / -1.9945669254d+00 , -3.4682360250d-01 /
      data fk(   7), fk(   8) /  2.2808349622d-01 ,  9.9865270518d-02 /
      data fk(   9), fk(  10) /  5.8545947775d-02 ,  1.6442902979d+00 /
      data fk(  11), fk(  12) / -8.4741887841d-02 , -1.3888148874d-02 /
      data fk(  13), fk(  14) / -1.7258943975d+00 ,  5.3026178389d-02 /
      data fk(  15), fk(  16) /  2.1815168154d-03 ,  3.3442057559d-04 /
      data fk(  17), fk(  18) /  4.6795130119d-01 , -5.1893352007d-03 /
      data fk(  19), fk(  20) / -4.9370272817d-04 , -6.3507812861d-05 /
      data fk(  21), fk(  22) /  3.0804585498d-01 , -1.2349651555d-02 /
      data fk(  23), fk(  24) /  2.5269436408d-04 ,  8.7722956173d-06 /
      data fk(  25), fk(  26) /  6.7559818569d-07 , -2.8510210024d-01 /
      data fk(  27), fk(  28) /  7.0908564644d-03 , -7.1944245774d-05 /
      data fk(  29), fk(  30) / -1.6228665657d-06 , -1.0712147655d-07 /
      data fk(  31), fk(  32) /  5.6446719326d-02 , -9.6452657894d-04 /
      data fk(  33), fk(  34) / -6.3040065609d-06 ,  6.3755351518d-07 /
      data fk(  35), fk(  36) /  1.3273430819d-08 ,  7.0336769435d-10 /
      data fk(  37), fk(  38) /  2.9172886634d-02 , -6.3580345755d-04 /
      data fk(  39), fk(  40) /  1.2044184862d-05 , -1.9585125981d-07 /
      data fk(  41), fk(  42) / -2.2246749871d-09 , -9.5500472395d-11 /
      data fk(  43), fk(  44) / -2.8750029444d-02 ,  4.1551144054d-04 /
      data fk(  45), fk(  46) / -4.0945695422d-06 ,  2.3049334902d-08 /
      data fk(  47), fk(  48) /  7.3476178518d-10 ,  1.1330482961d-11 /
      data fk(  49), fk(  50) /  4.9057144968d-13 ,  1.1666096833d-03 /
      data fk(  51), fk(  52) / -3.8220546810d-05 ,  7.9914778054d-09 /
      data fk(  53), fk(  54) /  1.0115605358d-08 , -2.4368477321d-10 /
      data fk(  55), fk(  56) / -1.7958862788d-12 , -6.4544264012d-14 /
      data fk(  57), fk(  58) /  6.0686066213d-03 , -1.0371806394d-04 /
      data fk(  59), fk(  60) /  6.7171677203d-07 , -6.2119452239d-09 /
      data fk(  61), fk(  62) /  5.2461327310d-11 ,  4.9204245687d-13 /
      data fk(  63), fk(  64) /  8.3089567712d-15 , -9.1627537080d-03 /
      data fk(  65), fk(  66) /  5.7584533332d-05 , -2.6915738121d-07 /
      data fk(  67), fk(  68) /  1.0365497801d-09 ,  6.9697334597d-13 /
      data fk(  69) / -1.5833746927d-13 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )
      cc(12) = cos (12.0d 00 * theta )
      cc(14) = cos (14.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 
      en=en+fk(43)*plgndr(12,0,cxj) 
      en=en+fk(44)*plgndr(12,2,cxj)*cc(2) 
      en=en+fk(45)*plgndr(12,4,cxj)*cc(4) 
      en=en+fk(46)*plgndr(12,6,cxj)*cc(6) 
      en=en+fk(47)*plgndr(12,8,cxj)*cc(8) 
      en=en+fk(48)*plgndr(12,10,cxj)*cc(10) 
      en=en+fk(49)*plgndr(12,12,cxj)*cc(12) 
      en=en+fk(50)*plgndr(13,0,cxj) 
      en=en+fk(51)*plgndr(13,2,cxj)*cc(2) 
      en=en+fk(52)*plgndr(13,4,cxj)*cc(4) 
      en=en+fk(53)*plgndr(13,6,cxj)*cc(6) 
      en=en+fk(54)*plgndr(13,8,cxj)*cc(8) 
      en=en+fk(55)*plgndr(13,10,cxj)*cc(10) 
      en=en+fk(56)*plgndr(13,12,cxj)*cc(12) 
      en=en+fk(57)*plgndr(14,0,cxj) 
      en=en+fk(58)*plgndr(14,2,cxj)*cc(2) 
      en=en+fk(59)*plgndr(14,4,cxj)*cc(4) 
      en=en+fk(60)*plgndr(14,6,cxj)*cc(6) 
      en=en+fk(61)*plgndr(14,8,cxj)*cc(8) 
      en=en+fk(62)*plgndr(14,10,cxj)*cc(10) 
      en=en+fk(63)*plgndr(14,12,cxj)*cc(12) 
      en=en+fk(64)*plgndr(15,0,cxj) 
      en=en+fk(65)*plgndr(15,2,cxj)*cc(2) 
      en=en+fk(66)*plgndr(15,4,cxj)*cc(4) 
      en=en+fk(67)*plgndr(15,6,cxj)*cc(6) 
      en=en+fk(68)*plgndr(15,8,cxj)*cc(8) 
      en=en+fk(69)*plgndr(15,10,cxj)*cc(10) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp5p5(thetain,phiin,es1,es2,fi)
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(69)
      dimension cc(14)
c
      data fk(   1), fk(   2) /  3.1384913315d-01 , -8.0229168222d-01 /
      data fk(   3), fk(   4) /  3.6553776308d+00 ,  1.2310913033d+00 /
      data fk(   5), fk(   6) / -4.5982739457d+00 , -4.2803058166d-01 /
      data fk(   7), fk(   8) /  8.3031068651d-01 ,  1.5473798595d-01 /
      data fk(   9), fk(  10) /  9.8417688269d-02 ,  3.0327035308d+00 /
      data fk(  11), fk(  12) / -1.6845676821d-01 , -2.2619420165d-02 /
      data fk(  13), fk(  14) / -2.9614362286d+00 ,  9.3494821835d-02 /
      data fk(  15), fk(  16) /  3.6840087541d-03 ,  5.5656642076d-04 /
      data fk(  17), fk(  18) /  7.1128861512d-01 , -6.0274384921d-03 /
      data fk(  19), fk(  20) / -9.0100708709d-04 , -1.0325731786d-04 /
      data fk(  21), fk(  22) /  3.5573648619d-01 , -1.6162599819d-02 /
      data fk(  23), fk(  24) /  3.4842280683d-04 ,  1.5319050747d-05 /
      data fk(  25), fk(  26) /  1.0933369867d-06 , -2.4820961483d-01 /
      data fk(  27), fk(  28) /  6.8489917120d-03 , -6.0725314762d-05 /
      data fk(  29), fk(  30) / -2.9417474259d-06 , -1.7528693086d-07 /
      data fk(  31), fk(  32) /  2.7112686332d-02 , -4.6612535455d-04 /
      data fk(  33), fk(  34) / -1.5906540284d-05 ,  7.5857354975d-07 /
      data fk(  35), fk(  36) /  2.5985772858d-08 ,  1.0585160880d-09 /
      data fk(  37), fk(  38) /  6.7388123975d-02 , -1.2140169108d-03 /
      data fk(  39), fk(  40) /  1.6418460046d-05 , -1.4887157019d-07 /
      data fk(  41), fk(  42) / -4.3767711613d-09 , -1.5395950640d-10 /
      data fk(  43), fk(  44) / -9.4972359240d-02 ,  1.1559174306d-03 /
      data fk(  45), fk(  46) / -6.9410335191d-06 , -5.8128090180d-10 /
      data fk(  47), fk(  48) /  7.7736673217d-10 ,  2.4614327674d-11 /
      data fk(  49), fk(  50) /  6.3907868332d-13 ,  1.8160253329d-02 /
      data fk(  51), fk(  52) / -9.9031151799d-05 , -1.1331683342d-06 /
      data fk(  53), fk(  54) /  3.2310680627d-08 , -2.3482966305d-10 /
      data fk(  55), fk(  56) / -3.8015461468d-12 , -1.0165131437d-13 /
      data fk(  57), fk(  58) /  3.8551180020d-02 , -4.6275398650d-04 /
      data fk(  59), fk(  60) /  2.6021560531d-06 , -1.8320041125d-08 /
      data fk(  61), fk(  62) /  5.5460675498d-11 ,  6.2511324449d-13 /
      data fk(  63), fk(  64) /  1.8689954397d-14 , -2.5183380976d-02 /
      data fk(  65), fk(  66) /  2.3799772987d-04 , -8.7309445681d-07 /
      data fk(  67), fk(  68) /  2.6942156650d-09 ,  2.0622165717d-11 /
      data fk(  69) / -2.2373008293d-13 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )
      cc(12) = cos (12.0d 00 * theta )
      cc(14) = cos (14.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 
      en=en+fk(43)*plgndr(12,0,cxj) 
      en=en+fk(44)*plgndr(12,2,cxj)*cc(2) 
      en=en+fk(45)*plgndr(12,4,cxj)*cc(4) 
      en=en+fk(46)*plgndr(12,6,cxj)*cc(6) 
      en=en+fk(47)*plgndr(12,8,cxj)*cc(8) 
      en=en+fk(48)*plgndr(12,10,cxj)*cc(10) 
      en=en+fk(49)*plgndr(12,12,cxj)*cc(12) 
      en=en+fk(50)*plgndr(13,0,cxj) 
      en=en+fk(51)*plgndr(13,2,cxj)*cc(2) 
      en=en+fk(52)*plgndr(13,4,cxj)*cc(4) 
      en=en+fk(53)*plgndr(13,6,cxj)*cc(6) 
      en=en+fk(54)*plgndr(13,8,cxj)*cc(8) 
      en=en+fk(55)*plgndr(13,10,cxj)*cc(10) 
      en=en+fk(56)*plgndr(13,12,cxj)*cc(12) 
      en=en+fk(57)*plgndr(14,0,cxj) 
      en=en+fk(58)*plgndr(14,2,cxj)*cc(2) 
      en=en+fk(59)*plgndr(14,4,cxj)*cc(4) 
      en=en+fk(60)*plgndr(14,6,cxj)*cc(6) 
      en=en+fk(61)*plgndr(14,8,cxj)*cc(8) 
      en=en+fk(62)*plgndr(14,10,cxj)*cc(10) 
      en=en+fk(63)*plgndr(14,12,cxj)*cc(12) 
      en=en+fk(64)*plgndr(15,0,cxj) 
      en=en+fk(65)*plgndr(15,2,cxj)*cc(2) 
      en=en+fk(66)*plgndr(15,4,cxj)*cc(4) 
      en=en+fk(67)*plgndr(15,6,cxj)*cc(6) 
      en=en+fk(68)*plgndr(15,8,cxj)*cc(8) 
      en=en+fk(69)*plgndr(15,10,cxj)*cc(10) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp5p0(thetain,phiin,es1,es2,fi)
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(72)
      dimension cc(14)
c
      data fk(   1), fk(   2) / -1.3416052255d+00 ,  1.5408680964d+00 /
      data fk(   3), fk(   4) /  7.0204776189d+00 ,  4.2689007871d-01 /
      data fk(   5), fk(   6) / -1.0307175298d+01 ,  1.1494357827d-02 /
      data fk(   7), fk(   8) /  2.7023439525d+00 ,  1.2281605040d-01 /
      data fk(   9), fk(  10) /  1.3198245497d-01 ,  3.8062719733d+00 /
      data fk(  11), fk(  12) / -2.3414309823d-01 , -2.8367867807d-02 /
      data fk(  13), fk(  14) / -2.2708092271d+00 ,  7.3333823408d-02 /
      data fk(  15), fk(  16) /  5.0367500938d-03 ,  6.8890851264d-04 /
      data fk(  17), fk(  18) / -1.4802424275d-02 ,  5.7597339090d-03 /
      data fk(  19), fk(  20) / -7.9657548794d-04 , -1.2659862666d-04 /
      data fk(  21), fk(  22) / -1.1940704509d+00 ,  2.7505118803d-02 /
      data fk(  23), fk(  24) / -4.2361697673d-04 ,  2.4977665817d-05 /
      data fk(  25), fk(  26) /  1.1247242289d-06 ,  1.7728155168d+00 /
      data fk(  27), fk(  28) / -3.2670394868d-02 ,  3.4261446115d-04 /
      data fk(  29), fk(  30) / -3.3088271593d-06 , -2.1213798875d-07 /
      data fk(  31), fk(  32) / -1.1677818344d+00 ,  6.8949318845d-03 /
      data fk(  33), fk(  34) / -1.2111384521d-07 , -1.8714769502d-06 /
      data fk(  35), fk(  36) /  5.6459460428d-08 ,  7.4913451761d-10 /
      data fk(  37), fk(  38) /  4.8525675670d-01 ,  6.6872204184d-04 /
      data fk(  39), fk(  40) / -5.1790480617d-05 ,  1.1201880921d-06 /
      data fk(  41), fk(  42) / -8.0092890590d-09 , -1.7730792086d-10 /
      data fk(  43), fk(  44) /  5.6830469144d-01 ,  2.4740873407d-03 /
      data fk(  45), fk(  46) /  5.7576482575d-06 , -1.7209562497d-07 /
      data fk(  47), fk(  48) / -2.7204654283d-09 ,  6.4391350833d-11 /
      data fk(  49), fk(  50) / -9.0464137481d-14 , -6.9349052067d-01 /
      data fk(  51), fk(  52) /  7.3082224429d-04 , -1.1725770735d-05 /
      data fk(  53), fk(  54) /  7.8535288841d-08 ,  1.0130223432d-09 /
      data fk(  55), fk(  56) / -9.0255035048d-12 , -1.0870160929d-13 /
      data fk(  57), fk(  58) /  1.2429782710d-01 , -4.5209632174d-03 /
      data fk(  59), fk(  60) /  2.0681676447d-05 , -8.2668668314d-08 /
      data fk(  61), fk(  62) / -4.4954342073d-11 , -1.3749247492d-12 /
      data fk(  63), fk(  64) /  7.9364567463d-14 , -2.1014445130d-02 /
      data fk(  65), fk(  66) /  2.7072954968d-03 , -8.4672222188d-06 /
      data fk(  67), fk(  68) /  1.3201541040d-08 ,  1.3471072425d-10 /
      data fk(  69), fk(  70) /  4.1361266372d-13 ,  2.1664588465d-12 /
      data fk(  71), fk(  72) / -8.2810181066d-13 ,  7.4557179870d-03 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )
      cc(12) = cos (12.0d 00 * theta )
      cc(14) = cos (14.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 
      en=en+fk(43)*plgndr(12,0,cxj) 
      en=en+fk(44)*plgndr(12,2,cxj)*cc(2) 
      en=en+fk(45)*plgndr(12,4,cxj)*cc(4) 
      en=en+fk(46)*plgndr(12,6,cxj)*cc(6) 
      en=en+fk(47)*plgndr(12,8,cxj)*cc(8) 
      en=en+fk(48)*plgndr(12,10,cxj)*cc(10) 
      en=en+fk(49)*plgndr(12,12,cxj)*cc(12) 
      en=en+fk(50)*plgndr(13,0,cxj) 
      en=en+fk(51)*plgndr(13,2,cxj)*cc(2) 
      en=en+fk(52)*plgndr(13,4,cxj)*cc(4) 
      en=en+fk(53)*plgndr(13,6,cxj)*cc(6) 
      en=en+fk(54)*plgndr(13,8,cxj)*cc(8) 
      en=en+fk(55)*plgndr(13,10,cxj)*cc(10) 
      en=en+fk(56)*plgndr(13,12,cxj)*cc(12) 
      en=en+fk(57)*plgndr(14,0,cxj) 
      en=en+fk(58)*plgndr(14,2,cxj)*cc(2) 
      en=en+fk(59)*plgndr(14,4,cxj)*cc(4) 
      en=en+fk(60)*plgndr(14,6,cxj)*cc(6) 
      en=en+fk(61)*plgndr(14,8,cxj)*cc(8) 
      en=en+fk(62)*plgndr(14,10,cxj)*cc(10) 
      en=en+fk(63)*plgndr(14,12,cxj)*cc(12) 
      en=en+fk(64)*plgndr(15,0,cxj) 
      en=en+fk(65)*plgndr(15,2,cxj)*cc(2) 
      en=en+fk(66)*plgndr(15,4,cxj)*cc(4) 
      en=en+fk(67)*plgndr(15,6,cxj)*cc(6) 
      en=en+fk(68)*plgndr(15,8,cxj)*cc(8) 
      en=en+fk(69)*plgndr(15,10,cxj)*cc(10) 
      en=en+fk(70)*plgndr(16,10,cxj) 
      en=en+fk(71)*plgndr(17,10,cxj) 
      en=en+fk(72)*plgndr(19,1,cxj) 

      fi = (en-es1)/(es2-es1)
      return
      end
      subroutine r2dp4p5(thetain,phiin,es1,es2,fi)
      implicit real*8 (a-h,o-z)
      parameter (degrad=57.29577951308232d 00)

      dimension fk(72)
      dimension cc(14)
c
      data fk(   1), fk(   2) / -4.4559935732d+00 ,  6.4627170003d+00 /
      data fk(   3), fk(   4) /  1.1983910502d+01 , -1.7174290014d+00 /
      data fk(   5), fk(   6) / -1.7855125064d+01 ,  9.2802378384d-01 /
      data fk(   7), fk(   8) /  3.9391076711d+00 ,  2.1204867268d-02 /
      data fk(   9), fk(  10) /  1.5663313557d-01 ,  2.9666741363d+00 /
      data fk(  11), fk(  12) / -1.2135828370d-01 , -3.3858462902d-02 /
      data fk(  13), fk(  14) /  3.3118171662d+00 , -1.9741344262d-01 /
      data fk(  15), fk(  16) /  8.2793382697d-03 ,  6.7356465333d-04 /
      data fk(  17), fk(  18) / -3.0158136785d+00 ,  9.3117829674d-02 /
      data fk(  19), fk(  20) /  1.7094643701d-04 , -1.3881828611d-04 /
      data fk(  21), fk(  22) / -4.2344560826d+00 ,  1.4434014353d-01 /
      data fk(  23), fk(  24) / -3.0440595849d-03 ,  4.7577241099d-05 /
      data fk(  25), fk(  26) /  6.5059414195d-07 ,  5.2605495047d+00 /
      data fk(  27), fk(  28) / -1.3031058379d-01 ,  1.5381851333d-03 /
      data fk(  29), fk(  30) / -3.7918742760d-06 , -2.1148192758d-07 /
      data fk(  31), fk(  32) / -1.5370237596d+00 ,  2.4958179840d-02 /
      data fk(  33), fk(  34) / -9.0605472579d-06 , -9.1883320919d-06 /
      data fk(  35), fk(  36) /  1.2173625146d-07 , -3.4050131964d-10 /
      data fk(  37), fk(  38) /  2.8970244242d-01 , -9.0672318587d-03 /
      data fk(  39), fk(  40) / -1.3667451455d-04 ,  4.3998407965d-06 /
      data fk(  41), fk(  42) / -1.6515007313d-08 , -9.6916573106d-11 /
      data fk(  43), fk(  44) / -2.1212464295d+00 ,  1.9443888298d-02 /
      data fk(  45), fk(  46) / -1.5077011567d-05 , -5.8064415270d-07 /
      data fk(  47), fk(  48) / -1.2195826841d-08 ,  1.3900214319d-10 /
      data fk(  49), fk(  50) / -1.4059235198d-12 ,  1.3476098818d+00 /
      data fk(  51), fk(  52) /  3.4235196299d-04 , -4.6922494367d-05 /
      data fk(  53), fk(  54) /  3.2745878511d-07 ,  3.9354915790d-09 /
      data fk(  55), fk(  56) / -2.1599286975d-11 ,  1.3415863663d-13 /
      data fk(  57), fk(  58) /  1.2312253181d+00 , -1.2476082309d-02 /
      data fk(  59), fk(  60) /  7.0470030698d-05 , -3.5978048266d-07 /
      data fk(  61), fk(  62) / -2.1177322374d-10 , -7.6407324680d-12 /
      data fk(  63), fk(  64) /  1.9241038999d-13 , -1.0383715491d+00 /
      data fk(  65), fk(  66) /  4.7175435191d-03 , -1.6287863399d-05 /
      data fk(  67), fk(  68) /  5.1383916392d-08 ,  9.2783770265d-10 /
      data fk(  69), fk(  70) /  3.3751023522d-12 , -3.3951540936d-13 /
      data fk(  71), fk(  72) /  1.1883527294d-12 , -6.6080696967d-03 /
c
      theta = thetain/degrad
      cxj = cos ( theta)
c
      theta = phiin/degrad
      cc( 2) = cos ( 2.0d 00 * theta )
      cc( 4) = cos ( 4.0d 00 * theta )
      cc( 6) = cos ( 6.0d 00 * theta )
      cc( 8) = cos ( 8.0d 00 * theta )
      cc(10) = cos (10.0d 00 * theta )
      cc(12) = cos (12.0d 00 * theta )
      cc(14) = cos (14.0d 00 * theta )

      en = 0.0d 00

      en=en+fk(1) 
      en=en+fk(2)*plgndr(1,0,cxj) 
      en=en+fk(3)*plgndr(2,0,cxj) 
      en=en+fk(4)*plgndr(2,2,cxj)*cc(2) 
      en=en+fk(5)*plgndr(3,0,cxj) 
      en=en+fk(6)*plgndr(3,2,cxj)*cc(2) 
      en=en+fk(7)*plgndr(4,0,cxj) 
      en=en+fk(8)*plgndr(4,2,cxj)*cc(2) 
      en=en+fk(9)*plgndr(4,4,cxj)*cc(4) 
      en=en+fk(10)*plgndr(5,0,cxj) 
      en=en+fk(11)*plgndr(5,2,cxj)*cc(2) 
      en=en+fk(12)*plgndr(5,4,cxj)*cc(4) 
      en=en+fk(13)*plgndr(6,0,cxj) 
      en=en+fk(14)*plgndr(6,2,cxj)*cc(2) 
      en=en+fk(15)*plgndr(6,4,cxj)*cc(4) 
      en=en+fk(16)*plgndr(6,6,cxj)*cc(6) 
      en=en+fk(17)*plgndr(7,0,cxj) 
      en=en+fk(18)*plgndr(7,2,cxj)*cc(2) 
      en=en+fk(19)*plgndr(7,4,cxj)*cc(4) 
      en=en+fk(20)*plgndr(7,6,cxj)*cc(6) 
      en=en+fk(21)*plgndr(8,0,cxj) 
      en=en+fk(22)*plgndr(8,2,cxj)*cc(2) 
      en=en+fk(23)*plgndr(8,4,cxj)*cc(4) 
      en=en+fk(24)*plgndr(8,6,cxj)*cc(6) 
      en=en+fk(25)*plgndr(8,8,cxj)*cc(8) 
      en=en+fk(26)*plgndr(9,0,cxj) 
      en=en+fk(27)*plgndr(9,2,cxj)*cc(2) 
      en=en+fk(28)*plgndr(9,4,cxj)*cc(4) 
      en=en+fk(29)*plgndr(9,6,cxj)*cc(6) 
      en=en+fk(30)*plgndr(9,8,cxj)*cc(8) 
      en=en+fk(31)*plgndr(10,0,cxj) 
      en=en+fk(32)*plgndr(10,2,cxj)*cc(2) 
      en=en+fk(33)*plgndr(10,4,cxj)*cc(4) 
      en=en+fk(34)*plgndr(10,6,cxj)*cc(6) 
      en=en+fk(35)*plgndr(10,8,cxj)*cc(8) 
      en=en+fk(36)*plgndr(10,10,cxj)*cc(10) 
      en=en+fk(37)*plgndr(11,0,cxj) 
      en=en+fk(38)*plgndr(11,2,cxj)*cc(2) 
      en=en+fk(39)*plgndr(11,4,cxj)*cc(4) 
      en=en+fk(40)*plgndr(11,6,cxj)*cc(6) 
      en=en+fk(41)*plgndr(11,8,cxj)*cc(8) 
      en=en+fk(42)*plgndr(11,10,cxj)*cc(10) 
      en=en+fk(43)*plgndr(12,0,cxj) 
      en=en+fk(44)*plgndr(12,2,cxj)*cc(2) 
      en=en+fk(45)*plgndr(12,4,cxj)*cc(4) 
      en=en+fk(46)*plgndr(12,6,cxj)*cc(6) 
      en=en+fk(47)*plgndr(12,8,cxj)*cc(8) 
      en=en+fk(48)*plgndr(12,10,cxj)*cc(10) 
      en=en+fk(49)*plgndr(12,12,cxj)*cc(12) 
      en=en+fk(50)*plgndr(13,0,cxj) 
      en=en+fk(51)*plgndr(13,2,cxj)*cc(2) 
      en=en+fk(52)*plgndr(13,4,cxj)*cc(4) 
      en=en+fk(53)*plgndr(13,6,cxj)*cc(6) 
      en=en+fk(54)*plgndr(13,8,cxj)*cc(8) 
      en=en+fk(55)*plgndr(13,10,cxj)*cc(10) 
      en=en+fk(56)*plgndr(13,12,cxj)*cc(12) 
      en=en+fk(57)*plgndr(14,0,cxj) 
      en=en+fk(58)*plgndr(14,2,cxj)*cc(2) 
      en=en+fk(59)*plgndr(14,4,cxj)*cc(4) 
      en=en+fk(60)*plgndr(14,6,cxj)*cc(6) 
      en=en+fk(61)*plgndr(14,8,cxj)*cc(8) 
      en=en+fk(62)*plgndr(14,10,cxj)*cc(10) 
      en=en+fk(63)*plgndr(14,12,cxj)*cc(12) 
      en=en+fk(64)*plgndr(15,0,cxj) 
      en=en+fk(65)*plgndr(15,2,cxj)*cc(2) 
      en=en+fk(66)*plgndr(15,4,cxj)*cc(4) 
      en=en+fk(67)*plgndr(15,6,cxj)*cc(6) 
      en=en+fk(68)*plgndr(15,8,cxj)*cc(8) 
      en=en+fk(69)*plgndr(15,10,cxj)*cc(10) 
      en=en+fk(70)*plgndr(16,10,cxj) 
      en=en+fk(71)*plgndr(17,10,cxj) 
      en=en+fk(72)*plgndr(19,1,cxj) 

      fi = (en-es1)/(es2-es1)
      return
      end
